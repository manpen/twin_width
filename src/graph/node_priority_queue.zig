const std = @import("std");
const bitset = @import("../util/two_level_bitset.zig");
const Graph = @import("graph.zig").Graph;

pub fn NodePriorityQueue(comptime T: type) type {
    return struct {
        const Self = @This();
        priority_queue: std.PriorityQueue(NodePriority, void, NodePriority.compareFn),
        postponed_queue: std.PriorityQueue(PostponedNode, void, PostponedNode.postponedQueueOrder),
        monotonic_tic_counter: T,
        max_postponed: T,
        postpone_set: bitset.FastBitSet,
        graph: *Graph(T),

        pub const PostponedNode = struct {
            node: T,
            priority: T,
            pub fn postponedQueueOrder(ctx: void, first: PostponedNode, second: PostponedNode) std.math.Order {
                _ = ctx;
                return std.math.order(first.priority, second.priority);
            }
        };

        pub inline fn init(graph: *Graph(T), max_postponed: T) !Self {
            var node_prio_queue = Self{ .priority_queue = std.PriorityQueue(NodePriority, void, NodePriority.compareFn).init(graph.allocator, {}), .graph = graph, .monotonic_tic_counter = 0, .postponed_queue = std.PriorityQueue(PostponedNode, void, PostponedNode.postponedQueueOrder).init(graph.allocator, {}), .max_postponed = max_postponed, .postpone_set = try bitset.FastBitSet.initEmpty(graph.number_of_nodes,graph.allocator) };
						try node_prio_queue.priority_queue.ensureTotalCapacity(graph.number_of_nodes);
						try node_prio_queue.postponed_queue.ensureTotalCapacity(max_postponed+1);
						return node_prio_queue;
        }

        pub inline fn deinit(self: *Self) void {
            self.postponed_queue.deinit();
            self.priority_queue.deinit();
            self.postpone_set.deinit(self.graph.allocator);
        }

				pub inline fn checkSmallest(self: *Self, node: T) bool {
					if(self.priority_queue.peek()) |item| {

            const prio_new = NodePriority{ .red_degree = self.graph.node_list[node].red_edges.cardinality(), .black_degree = self.graph.node_list[node].black_edges.cardinality(), .id = node };

						if(NodePriority.compareFn({},prio_new,item) == std.math.Order.lt) {
							return true;
						}
					}
					return false;
				}

        pub inline fn add(self: *Self, node: T) !void {
            if (self.graph.erased_nodes.get(node)) {
                return;
            }
            try self.priority_queue.add(NodePriority{ .red_degree = self.graph.node_list[node].red_edges.cardinality(), .black_degree = self.graph.node_list[node].black_edges.cardinality(), .id = node });
        }

        pub inline fn postponedIterator(self: *Self) std.PriorityQueue(PostponedNode, void, PostponedNode.postponedQueueOrder).Iterator {
            return self.postponed_queue.iterator();
        }

        pub inline fn shouldPostpone(self: *Self, current_tww: T, induced_tww: T) bool {
            _ = self;
						// Adaptive scheduling might improve this further
            return induced_tww > current_tww and (induced_tww * 5) > (current_tww * 6);
						//return induced_tww > current_tww;
        }

        pub inline fn postpone(self: *Self, node: T, by: T) !bool {
            if (self.postpone_set.get(node)) {
                return false;
            }
            // Deny postpone
            else if (self.max_postponed <= self.postponed_queue.len) {
                if (self.postponed_queue.peek()) |item| {
                    if (item.priority < (self.monotonic_tic_counter + by)) {
                        const removed_item = self.postponed_queue.remove();
                        _ = self.postpone_set.unset(removed_item.node);
                        try self.add(removed_item.node);
                    } else {
                        return false;
                    }
                } else {
                    return false;
                }
            }
            try self.postponed_queue.add(PostponedNode{ .priority = self.monotonic_tic_counter + by, .node = node });
            self.postpone_set.set(node);
            return true;
        }

        pub inline fn increaseTick(self: *Self) void {
            self.monotonic_tic_counter += 1;
        }

        pub inline fn addTick(self: *Self, ticks: T) void {
            self.monotonic_tic_counter += ticks;
        }

        pub inline fn reconcilePostponed(self: *Self, current_tww: T) !void {
            while (self.postponed_queue.peek()) |item| {
                if (item.priority < (self.monotonic_tic_counter + current_tww)) {
                    const removed_item = self.postponed_queue.remove();
                    _ = self.postpone_set.unset(removed_item.node);
                    try self.add(removed_item.node);
                } else {
                    break;
                }
            }
        }

				pub inline fn clear(self: *Self) void {
					self.postponed_queue.len = 0;
					self.priority_queue.len = 0;

					self.postpone_set.unsetAll();
					self.monotonic_tic_counter = 0;
				}

        pub inline fn removeNext(self: *Self, current_tww: T) !T {
            try self.reconcilePostponed(current_tww);
            while (true) {
                if (self.priority_queue.removeOrNull()) |prio_item| {
                    if (self.graph.erased_nodes.get(prio_item.id) or self.postpone_set.get(prio_item.id)) {
                        continue;
                    } else if (self.graph.node_list[prio_item.id].red_edges.cardinality() != prio_item.red_degree) {
                        try self.add(prio_item.id);
                        continue;
                    }
                    return prio_item.id;
                } else {
                    // Ok we ran out of items force reconcile!
                    const item = self.postponed_queue.remove().node;
                    _ = self.postpone_set.unset(item);
                    try self.add(item);
                }
            }
        }

        const NodePriority = struct {
            red_degree: T,
            black_degree: T,
            id: T,
            pub fn compareFn(ctx: void, lhs: NodePriority, rhs: NodePriority) std.math.Order {
                _ = ctx;
								//Sort leafes to the end
								if(lhs.red_degree+lhs.black_degree == 1 or rhs.red_degree+rhs.red_degree == 1) {
									return std.math.order(rhs.black_degree+rhs.red_degree,lhs.red_degree+lhs.black_degree);
								}
                if (lhs.red_degree == rhs.red_degree) {
                    return std.math.order(lhs.black_degree, rhs.black_degree);
                }
                return std.math.order(lhs.red_degree, rhs.red_degree);
            }
        };
    };
}
