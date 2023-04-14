const std = @import("std");
const bitset = @import("../util/two_level_bitset.zig");
const edge_list = @import("edge_list.zig");
const red_edge_stack = @import("red_edge_stack.zig");
const contraction = @import("../tww/contraction_sequence.zig");
const comptime_util = @import("../util/comptime_checks.zig");
const RetraceableContractionSequence = @import("../tww/retraceable_contraction_sequence.zig").RetraceableContractionSequence;
const connected_components = @import("connected_component.zig");
const Node = @import("node.zig").Node;
const bfs_mod = @import("bfs.zig");
const compressed_bitset = @import("../util/compressed_bitmap.zig");
const solver_resources = @import("solver.zig");

const pace_2023 = @import("../pace_2023/pace_fmt.zig");

pub const GraphError = error{ FileNotFound, NotPACEFormat, GraphTooLarge, MisformedEdgeList, InvalidContractionOneNodeErased, ContractionOverflow, NoContractionLeft, RetraceNoParent };


pub fn Graph(comptime T: type) type {
    comptime if (!comptime_util.checkIfIsCompatibleInteger(T)) {
        @compileError("T must either be u8,u16 or u32!");
    };

    const promote_thresh = comptime if (T == u8) 0 else if (T == u16) 200 else 300;
    const degrade_tresh = comptime if (T == u8) 0 else if (T == u16) 100 else 200;

    const NodeType = Node(T, promote_thresh, degrade_tresh);

    return struct {
        const Self = @This();
        number_of_nodes: u32,
        number_of_edges: u32,
        node_list: []NodeType,
        erased_nodes: bitset.FastBitSet,
        scratch_bitset: bitset.FastBitSet,
        connected_components: std.ArrayListUnmanaged(connected_components.ConnectedComponent(T)),
        connected_components_min_heap: std.PriorityQueue(connected_components.ConnectedComponentIndex(T), void, connected_components.ConnectedComponentIndex(T).compareComponentIndexDesc),
        allocator: std.mem.Allocator,
        contraction: RetraceableContractionSequence(T),

        pub const LargeListStorageType = compressed_bitset.FastCompressedBitmap(T, promote_thresh, degrade_tresh);

        pub inline fn addEdge(self: *Self, u: T, v: T) !void {
            self.connected_components.clearRetainingCapacity();
            self.connected_components_min_heap.shrinkAndFree(0);
            std.debug.assert(u < self.node_list.len);
            std.debug.assert(v < self.node_list.len);

            const u_node = &self.node_list[u];
            const v_node = &self.node_list[v];

            _ = try u_node.black_edges.add(self.allocator, v);
            _ = try v_node.black_edges.add(self.allocator, u);

            self.number_of_edges += 1;
        }

        pub inline fn memoryUsage(self: *Self) usize {
            return self.allocator.end_index;
        }

        pub inline fn getCurrentTwinWidth(self: *const Self) u32 {
            return self.contraction.getTwinWidth();
        }

        pub fn reverseLastContraction(self: *Self) !void {
            _ = self;
        }

        pub const InducedTwinWidth = struct { 
					tww: T,
					delta_red_edges: i32,
					pub inline fn default() InducedTwinWidth {
						return InducedTwinWidth {
							.tww = std.math.maxInt(T),
							.delta_red_edges = std.math.maxInt(i32),
						};
					}

					pub inline fn reset(self: *InducedTwinWidth) void {
						self.tww = std.math.maxInt(T);
						self.delta_red_edges = std.math.maxInt(i32);
					}

					pub inline fn lessThan(self: InducedTwinWidth, other: InducedTwinWidth) bool {
						if(self.tww < other.tww or (self.tww == other.tww and self.delta_red_edges < other.delta_red_edges)) {
							return true;
						}
						else {
							return false;
						}
					}
				};

        pub fn calculateInducedTww(self: *Self, erased: T, survivor: T, upper_bound: ?T) InducedTwinWidth {
            const erased_cardinality_red = self.node_list[erased].red_edges.cardinality();
            const survivor_cardinality_red = self.node_list[survivor].red_edges.cardinality();
            if (upper_bound) |ub| {
                const erased_cardinality_black = self.node_list[erased].black_edges.cardinality();
                const survivor_cardinality_black = self.node_list[survivor].black_edges.cardinality();

                const erased_cardinality = erased_cardinality_black + erased_cardinality_red;
                const survivor_cardinality = survivor_cardinality_black + survivor_cardinality_red;

                //heuristic_024.gr
                const min_dist_black = std.math.max(erased_cardinality_black, survivor_cardinality_black) - std.math.min(erased_cardinality_black, survivor_cardinality_black);

                if (min_dist_black > ub) {
                    return InducedTwinWidth{
                        .tww = min_dist_black,
                        .delta_red_edges = std.math.maxInt(i32),
                    };
                }

                const min_dist_total = std.math.max(erased_cardinality, survivor_cardinality) - std.math.min(erased_cardinality, survivor_cardinality);

                // Fast exit
                if (min_dist_total > ub) {
                    return InducedTwinWidth{ .tww = min_dist_total, .delta_red_edges = std.math.maxInt(i32) };
                }
            }

            const upper_bound_all = upper_bound orelse std.math.maxInt(T);

            var delta_red: T = 0;
            var correction_factor: T = 0;
            if (erased_cardinality_red == 0 or survivor_cardinality_red == 0) {
                delta_red = std.math.max(erased_cardinality_red, survivor_cardinality_red);
            } else {
                var red_iter = self.node_list[erased].red_edges.xorIterator(&self.node_list[survivor].red_edges);

                while (red_iter.next()) |item| {
                    if (item != survivor and item != erased) {
                        delta_red += 1;
                    } else {
                        correction_factor = 1;
                    }
                }
            }

            var delta_red_edges = @intCast(i32, (delta_red + correction_factor)) - @intCast(i32, erased_cardinality_red + survivor_cardinality_red);

            if (delta_red > upper_bound_all) {
                return InducedTwinWidth{ .tww = delta_red, .delta_red_edges = delta_red_edges };
            }

            var tww: T = 0;
            var black_iter = self.node_list[erased].black_edges.xorIterator(&self.node_list[survivor].black_edges);

            while (black_iter.next()) |item| {
                if (item == survivor or item == erased) {
                    continue;
                }
                // Came from erased
                if (black_iter.first) {
                    if (!self.node_list[survivor].red_edges.contains(item)) {
                        delta_red += 1;
                        delta_red_edges += 1;
                        tww = std.math.max(self.node_list[item].red_edges.cardinality() + 1, tww);
                    }
                }
                // Came from survivor
                else {
                    delta_red += 1;
                    delta_red_edges += 1;
                    tww = std.math.max(self.node_list[item].red_edges.cardinality() + 1, tww);
                }

                if (tww > upper_bound_all or delta_red > upper_bound_all) {
                    return InducedTwinWidth{ .tww = std.math.max(tww, delta_red), .delta_red_edges = delta_red_edges };
                }
            }

            tww = std.math.max(tww, delta_red);
            return InducedTwinWidth{ .tww = tww, .delta_red_edges = delta_red_edges };
        }

				pub fn addContractionNewRedEdges(self: *Self, erased: T, survivor: T, new_red_edges: *std.ArrayListUnmanaged(T), first_level_merge: *bool ) !T {
            if (self.erased_nodes.get(erased) or self.erased_nodes.get(survivor)) {
                return GraphError.InvalidContractionOneNodeErased;
            }
						first_level_merge.* = false;
            self.erased_nodes.set(erased);
            var red_iter = self.node_list[erased].red_edges.iterator();
            while (red_iter.next()) |item| {
                if (item != survivor) {
                    if (!try self.node_list[survivor].addRedEdgeExists(self.allocator, item)) {
												try new_red_edges.append(self.allocator,item);
                        try self.node_list[item].addRedEdge(self.allocator, survivor);
                    }
                }
								else {
									first_level_merge.* = true;
								}
                try self.node_list[item].removeRedEdge(self.allocator, erased);
            }

            var tww: T = 0;
            var black_iter = self.node_list[erased].black_edges.xorIterator(&self.node_list[survivor].black_edges);

            var remove_list = std.ArrayList(T).init(self.allocator);
            defer remove_list.deinit();

            while (black_iter.next()) |item| {
                if (item == survivor or item == erased) {
										first_level_merge.* = true;
                    continue;
                }
                // Came from erased
                if (black_iter.first) {
                    if (!try self.node_list[survivor].addRedEdgeExists(self.allocator, item)) {
												try new_red_edges.append(self.allocator,item);
                        try self.node_list[item].addRedEdge(self.allocator, survivor);
                    }
                }
                // Came from survivor
                else {
                    try self.node_list[survivor].addRedEdge(self.allocator, item);
                    try self.node_list[item].addRedEdge(self.allocator, survivor);

                    try self.node_list[item].removeBlackEdge(self.allocator, survivor);
                    try remove_list.append(item);
                }
                tww = std.math.max(tww, @intCast(T, self.node_list[item].red_edges.cardinality()));
            }

            for (remove_list.items) |item| {
                try self.node_list[survivor].removeBlackEdge(self.allocator, item);
            }

            var black_remove_iter = self.node_list[erased].black_edges.iterator();
            while (black_remove_iter.next()) |item| {
                try self.node_list[item].removeBlackEdge(self.allocator, erased);
            }

            tww = std.math.max(tww, @intCast(T, self.node_list[survivor].red_edges.cardinality()));
            try self.contraction.addContraction(erased, survivor, std.math.max(tww, self.contraction.getTwinWidth()));
            return tww;
        }

        pub fn addContraction(self: *Self, erased: T, survivor: T) !T {
            if (self.erased_nodes.get(erased) or self.erased_nodes.get(survivor)) {
                return GraphError.InvalidContractionOneNodeErased;
            }
            self.erased_nodes.set(erased);
            var red_iter = self.node_list[erased].red_edges.iterator();
            while (red_iter.next()) |item| {
                if (item != survivor) {
                    if (!try self.node_list[survivor].addRedEdgeExists(self.allocator, item)) {
                        try self.node_list[item].addRedEdge(self.allocator, survivor);
                    }
                }
                try self.node_list[item].removeRedEdge(self.allocator, erased);
            }

            var tww: T = 0;
            var black_iter = self.node_list[erased].black_edges.xorIterator(&self.node_list[survivor].black_edges);

            var remove_list = std.ArrayList(T).init(self.allocator);
            defer remove_list.deinit();

            while (black_iter.next()) |item| {
                if (item == survivor or item == erased) {
                    continue;
                }
                // Came from erased
                if (black_iter.first) {
                    if (!try self.node_list[survivor].addRedEdgeExists(self.allocator, item)) {
                        try self.node_list[item].addRedEdge(self.allocator, survivor);
                    }
                }
                // Came from survivor
                else {
                    try self.node_list[survivor].addRedEdge(self.allocator, item);
                    try self.node_list[item].addRedEdge(self.allocator, survivor);

                    try self.node_list[item].removeBlackEdge(self.allocator, survivor);
                    try remove_list.append(item);
                }
                tww = std.math.max(tww, @intCast(T, self.node_list[item].red_edges.cardinality()));
            }

            for (remove_list.items) |item| {
                try self.node_list[survivor].removeBlackEdge(self.allocator, item);
            }

            var black_remove_iter = self.node_list[erased].black_edges.iterator();
            while (black_remove_iter.next()) |item| {
                try self.node_list[item].removeBlackEdge(self.allocator, erased);
            }

            tww = std.math.max(tww, @intCast(T, self.node_list[survivor].red_edges.cardinality()));
            try self.contraction.addContraction(erased, survivor, std.math.max(tww, self.contraction.getTwinWidth()));
            return tww;
        }

        pub fn solveGreedy(self: *Self) !T {
						var solver = try solver_resources.SolverResources(T,10,100).init(self);
						defer solver.deinit(self.allocator);
						var cc_iter = self.connected_components_min_heap.iterator();
						var tww:T = 0;
						var last_node:?T = null;
						while(cc_iter.next()) |cc| {
							tww = std.math.max(try self.connected_components.items[cc.index].solveGreedyTopK(10,100,&solver),tww);
							if(self.connected_components.items[cc.index].subgraph.nodes.cardinality() == 1) {
								const survivor = self.connected_components.items[cc.index].subgraph.nodes.edges.items[0];
								if(last_node) |ln| {
									_ = try self.addContraction(survivor,ln);
								}
								else {
									last_node = survivor;
								}
							}
							else if(self.contraction.lastContraction()) |ctr| {
								if(last_node) |ln| {
									_ = try self.addContraction(ctr.survivor,ln);
								}
								else {
									last_node = ctr.survivor;
								}
							}
							solver.reset();
						}
						//std.debug.print("Contraction len {} graph len {}\n",.{self.contraction.seq.isComplete(),self.number_of_nodes});
						return tww;
        }

        pub fn new(number_of_nodes: T, allocator: std.mem.Allocator) !Self {
            var node_list = try allocator.alloc(NodeType, number_of_nodes);

            for (node_list) |*node| {
                node.black_edges = compressed_bitset.FastCompressedBitmap(T, promote_thresh, degrade_tresh).init(number_of_nodes);
                node.red_edges = compressed_bitset.FastCompressedBitmap(T, promote_thresh, degrade_tresh).init(number_of_nodes);
								node.high_degree_node = null;
            }

            //TODO: Add some errdefer's here

            var graph = Self{
                .number_of_nodes = number_of_nodes,
                .number_of_edges = 0,
                .node_list = node_list,
                .allocator = allocator,
                .scratch_bitset = try bitset.FastBitSet.initEmpty(number_of_nodes, allocator),
                .erased_nodes = try bitset.FastBitSet.initEmpty(number_of_nodes, allocator),
                .connected_components = std.ArrayListUnmanaged(connected_components.ConnectedComponent(T)){},
                .connected_components_min_heap = std.PriorityQueue(connected_components.ConnectedComponentIndex(T), void, connected_components.ConnectedComponentIndex(T).compareComponentIndexDesc).init(allocator, {}),
                .contraction = try RetraceableContractionSequence(T).init(allocator, number_of_nodes),
            };

            return graph;
        }

        pub inline fn density(self: *const Self) f64 {
            return (2 * @intToFloat(f64, self.number_of_edges)) / (@intToFloat(f64, self.number_of_nodes) * @intToFloat(f64, self.number_of_nodes - 1));
        }

        pub fn loadFromPace(allocator: std.mem.Allocator, pace: *pace_2023.Pace2023Fmt(T)) !Self {
            var node_list = try allocator.alloc(NodeType, pace.number_of_nodes);

            for (0..pace.number_of_nodes) |index| {
                node_list[index].black_edges = try compressed_bitset.FastCompressedBitmap(T, promote_thresh, degrade_tresh).fromUnsorted(allocator, &pace.nodes[index].edges, @intCast(T, pace.number_of_nodes));
                node_list[index].red_edges = compressed_bitset.FastCompressedBitmap(T, promote_thresh, degrade_tresh).init(@intCast(T, pace.number_of_nodes));
								node_list[index].high_degree_node = null;
            }

            // Remove allocations on failure
            errdefer {
                for (node_list) |*node| {
                    node.black_edges.deinit(allocator);
                    node.red_edges.deinit(allocator);
										node.high_degree_node = null;
                }
                allocator.free(node_list);
            }

            return Self {
                .number_of_nodes = pace.number_of_nodes,
                .number_of_edges = pace.number_of_edges, //ATM
                .node_list = node_list,
                .allocator = allocator,
                .erased_nodes = try bitset.FastBitSet.initEmpty(pace.number_of_nodes, allocator),
                .scratch_bitset = try bitset.FastBitSet.initEmpty(pace.number_of_nodes, allocator),
                .connected_components = std.ArrayListUnmanaged(connected_components.ConnectedComponent(T)){},
                .connected_components_min_heap = std.PriorityQueue(connected_components.ConnectedComponentIndex(T), void, connected_components.ConnectedComponentIndex(T).compareComponentIndexDesc).init(allocator, {}),
                .contraction = try RetraceableContractionSequence(T).init(allocator, pace.number_of_nodes),
            };
        }

        pub fn findAllConnectedComponents(self: *Self) !void {
            self.scratch_bitset.unsetAll();

            var bfs_stack = try bfs_mod.BfsQueue(T).init(self.allocator, self.number_of_nodes);

            var unsetIter = self.scratch_bitset.iterUnset();
            var components: u32 = 0;
            var total_sum: u32 = 0;
            var largest: u32 = 0;
            var node_list = edge_list.ParametrizedUnsortedArrayList(T).init();

            while (unsetIter.next()) |item| {
                if (self.scratch_bitset.get(item)) {
                    continue;
                }
                var iterator = bfs_mod.bfs(T, @intCast(T,item), self, &self.scratch_bitset, &bfs_stack, .{ .max_level = std.math.maxInt(T), .kind = .black });

                while (iterator.next()) |node| {
                    self.scratch_bitset.set(node);
                    try node_list.add(self.allocator, node);
                }

                components += 1;
                total_sum += node_list.cardinality();
                largest = std.math.max(largest, node_list.cardinality());
                try self.connected_components.append(self.allocator, try connected_components.ConnectedComponent(T).init(self.allocator, node_list, iterator.level,self));
                node_list = edge_list.ParametrizedUnsortedArrayList(T).init();
            }

            try self.connected_components_min_heap.ensureTotalCapacity(self.connected_components.items.len);
            for (0..self.connected_components.items.len) |index| {
                try self.connected_components_min_heap.add(connected_components.ConnectedComponentIndex(T){ .tww = self.connected_components.items[index].tww, .index = @intCast(T, index) });
            }

            var tww = self.connected_components_min_heap.remove();
            //std.debug.print("Found {} components largest {} and tww {} density {}\n", .{ components, largest, tww.tww, self.density() });
            try self.connected_components_min_heap.add(tww);
        }

        pub fn deinit(self: *Self) void {
            for (self.node_list) |*node| {
                node.black_edges.deinit(self.allocator);
                node.red_edges.deinit(self.allocator);
            }
            for (self.connected_components.items) |*connected_component| {
                connected_component.deinit(self.allocator);
            }
            self.connected_components.deinit(self.allocator);
            self.contraction.deinit();
            self.scratch_bitset.deinit(self.allocator);
            self.erased_nodes.deinit(self.allocator);
            self.allocator.free(self.node_list);
        }
    };
}

test "Check simple loading" {
    var allocator = std.heap.GeneralPurposeAllocator(.{}){};
    defer std.debug.assert(!allocator.deinit());

    var graph = try Graph(u8).loadFromPace(allocator.allocator(), "instances/tiny/tiny001.gr");
    // Free ressources
    defer graph.deinit();

    try std.testing.expectEqual(graph.number_of_nodes, 10);
    try std.testing.expectEqual(graph.number_of_edges, 9);
}

test "Test contraction Tiny 1" {
    var allocator = std.heap.GeneralPurposeAllocator(.{}){};
    defer std.debug.assert(!allocator.deinit());

    var graph = try Graph(u8).loadFromPace(allocator.allocator(), "instances/tiny/tiny001.gr");
    // Free ressources
    defer graph.deinit();

    try std.testing.expectEqual(graph.number_of_nodes, 10);
    try std.testing.expectEqual(graph.number_of_edges, 9);

    var i: u8 = 9;
    while (i > 0) {
        var tww = try graph.addContraction(i - 1, 9);
        if (i > 1) {
            try std.testing.expectEqual(@as(u32, 1), tww);
        } else {
            try std.testing.expectEqual(@as(u32, 0), tww);
        }
        i -= 1;
    }
}

test "Test contraction Tiny 2" {
    var allocator = std.heap.GeneralPurposeAllocator(.{}){};
    defer std.debug.assert(!allocator.deinit());

    var graph = try Graph(u8).loadFromPace(allocator.allocator(), "instances/tiny/tiny002.gr");
    // Free ressources
    defer graph.deinit();

    try std.testing.expectEqual(graph.number_of_nodes, 10);
    try std.testing.expectEqual(graph.number_of_edges, 10);

    var i: u8 = 9;
    while (i > 0) {
        var tww = try graph.addContraction(i - 1, 9);
        if (i > 2) {
            try std.testing.expectEqual(@as(u32, 2), tww);
        } else if (i > 1) {
            try std.testing.expectEqual(@as(u32, 1), tww);
        } else {
            try std.testing.expectEqual(@as(u32, 0), tww);
        }
        i -= 1;
    }
}

test "Test contraction retrace Tiny 2" {
    var allocator = std.heap.GeneralPurposeAllocator(.{}){};
    defer std.debug.assert(!allocator.deinit());

    var graph = try Graph(u8).loadFromPace(allocator.allocator(), "instances/tiny/tiny002.gr");
    // Free ressources
    defer graph.deinit();

    try std.testing.expectEqual(graph.number_of_nodes, 10);
    try std.testing.expectEqual(graph.number_of_edges, 10);

    var i: u8 = 9;
    while (i > 0) {
        var tww = try graph.addContraction(i - 1, 9);
        if (i > 2) {
            try std.testing.expectEqual(@as(u32, 2), tww);
        } else if (i > 1) {
            try std.testing.expectEqual(@as(u32, 1), tww);
        } else {
            try std.testing.expectEqual(@as(u32, 0), tww);
        }
        i -= 1;
    }
    i = 9;
    while (i > 0) {
        try graph.reverseLastContraction();

        const tww = graph.getCurrentTwinWidth();

        if (i == 9) {
            try std.testing.expectEqual(@as(u32, 2), tww);
        } else if (i == 1) {
            try std.testing.expectEqual(@as(u32, 0), tww);
        } else {
            try std.testing.expectEqual(@as(u32, 2), tww);
        }
        for (graph.node_list) |*node| {
            try std.testing.expect(node.red_edges.cardinality() <= tww);
        }
        i -= 1;
    }
}

test "Test contraction Tiny 3" {
    var allocator = std.heap.GeneralPurposeAllocator(.{}){};
    defer std.debug.assert(!allocator.deinit());

    var graph = try Graph(u8).loadFromPace(allocator.allocator(), "instances/tiny/tiny003.gr");
    // Free ressources
    defer graph.deinit();

    try std.testing.expectEqual(graph.number_of_nodes, 10);
    try std.testing.expectEqual(graph.number_of_edges, 45);

    var i: u8 = 9;
    while (i > 0) {
        var tww = try graph.addContraction(i - 1, 9);
        try std.testing.expectEqual(@as(u32, 0), tww);
        i -= 1;
    }
}

test "Test contraction retrace Tiny 3" {
    var allocator = std.heap.GeneralPurposeAllocator(.{}){};
    defer std.debug.assert(!allocator.deinit());

    var graph = try Graph(u8).loadFromPace(allocator.allocator(), "instances/tiny/tiny003.gr");
    // Free ressources
    defer graph.deinit();

    try std.testing.expectEqual(graph.number_of_nodes, 10);
    try std.testing.expectEqual(graph.number_of_edges, 45);

    var i: u8 = 9;
    while (i > 0) {
        var tww = try graph.addContraction(i - 1, 9);
        try std.testing.expectEqual(@as(u32, 0), tww);
        i -= 1;
    }
    i = 9;
    while (i > 0) {
        try graph.reverseLastContraction();
        var tww = graph.getCurrentTwinWidth();
        try std.testing.expectEqual(@as(u32, 0), tww);
        for (graph.node_list) |*node| {
            try std.testing.expect(node.red_edges.cardinality() <= tww);
        }
        i -= 1;
    }
}

test "Test contraction Tiny 4" {
    var allocator = std.heap.GeneralPurposeAllocator(.{}){};
    defer std.debug.assert(!allocator.deinit());

    var graph = try Graph(u8).loadFromPace(allocator.allocator(), "instances/tiny/tiny004.gr");
    // Free ressources
    defer graph.deinit();

    try std.testing.expectEqual(graph.number_of_nodes, 10);
    try std.testing.expectEqual(graph.number_of_edges, 9);

    var i: u8 = 9;
    while (i > 0) {
        var tww = try graph.addContraction(i - 1, 9);
        try std.testing.expectEqual(@as(u32, 0), tww);
        i -= 1;
    }
}

test "Test contraction retrace Tiny 4" {
    var allocator = std.heap.GeneralPurposeAllocator(.{}){};
    defer std.debug.assert(!allocator.deinit());

    var graph = try Graph(u8).loadFromPace(allocator.allocator(), "instances/tiny/tiny004.gr");
    // Free ressources
    defer graph.deinit();

    try std.testing.expectEqual(graph.number_of_nodes, 10);
    try std.testing.expectEqual(graph.number_of_edges, 9);

    var i: u8 = 9;
    while (i > 0) {
        var tww = try graph.addContraction(i - 1, 9);
        try std.testing.expectEqual(@as(u32, 0), tww);
        i -= 1;
    }
    i = 9;
    while (i > 0) {
        try graph.reverseLastContraction();
        var tww = graph.getCurrentTwinWidth();
        try std.testing.expectEqual(@as(u32, 0), tww);
        for (graph.node_list) |*node| {
            try std.testing.expect(node.red_edges.cardinality() <= tww);
        }
        i -= 1;
    }
}

test "Test contraction Tiny 6" {
    var allocator = std.heap.GeneralPurposeAllocator(.{}){};
    defer std.debug.assert(!allocator.deinit());

    var graph = try Graph(u8).loadFromPace(allocator.allocator(), "instances/tiny/tiny006.gr");
    // Free ressources
    defer graph.deinit();

    try std.testing.expectEqual(graph.number_of_nodes, 10);
    try std.testing.expectEqual(graph.number_of_edges, 5);

    var i: u8 = 9;
    var tww: u32 = 0;
    while (true) {
        tww = std.math.max(tww, try graph.addContraction(i - 1, i));
        if (i == 1) {
            break;
        }
        i -= 2;
    }
    try std.testing.expectEqual(@as(u32, 0), tww);
}

test "Check simple failed loading" {
    var allocator = std.heap.GeneralPurposeAllocator(.{}){};
    // Should never fail since the ressources should be managed
    defer std.debug.assert(!allocator.deinit());

    var graph = Graph(u8).loadFromPace(allocator.allocator(), "instances/tiny/tiny001.graph") catch |err| {
        // Loading should fail with file not found
        try std.testing.expectEqual(err, error.FileNotFound);
        return;
    };
    _ = graph;
    // Should never be reached
    std.debug.assert(false);
}

test "Check simple graph" {
    var allocator = std.heap.GeneralPurposeAllocator(.{}){};
    // Should never fail since the ressources should be managed
    //defer std.debug.assert(!allocator.deinit());

    var graph = try Graph(u8).new(5, allocator.allocator());

    try graph.addEdge(0, 1);
    try graph.addEdge(1, 2);
    try graph.addEdge(1, 3);
    try graph.addEdge(1, 4);

    _ = try graph.addContraction(1, 0);

    try std.testing.expectEqual(graph.node_list[0].red_edges.cardinality(), 3);
}
