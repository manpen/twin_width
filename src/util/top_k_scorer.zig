const std = @import("std");
const comptime_util = @import("comptime_checks.zig");
const set = @import("../util/two_level_bitset.zig");
const Graph = @import("../graph/graph.zig").Graph;

pub const TopKError = error {
	ArrayTooSmall,
};

pub fn TopKScorer(comptime T: type, comptime K: u32) type {
	comptime if (!comptime_util.checkIfIsCompatibleInteger(T)) {
		@compileError("Type must be either u8,u16 or u32");
	};

	return struct {
		const Self = @This();

		pub const TopKQueueEntry = struct {
			node: T,
			score: i32
		};

		pub fn entryCompareFn(context: void, lhs: TopKQueueEntry,rhs: TopKQueueEntry) std.math.Order {
			_ = context;
			return std.math.order(rhs.score,lhs.score);
		}

		pub fn lessThan(context: void, lhs: TopKQueueEntry,rhs: TopKQueueEntry) bool {
			_ = context;
			return lhs.score < rhs.score;
		}

		backing_storage: []i32,
		priority_queue: std.PriorityQueue(TopKQueueEntry,void, entryCompareFn),

		node_list: []T,
		write_ptr: u32,

		sorted_targets: std.BoundedArray(TopKQueueEntry,K),

		pub fn init(allocator: std.mem.Allocator, total_number_of_nodes: T) !Self {
			var memory = try allocator.alloc(i32, total_number_of_nodes);
			std.mem.set(i32, memory,0);

			var nodel = try allocator.alloc(T, total_number_of_nodes+1);
			std.mem.set(T, nodel,0);
			
			var priority_queue = std.PriorityQueue(TopKQueueEntry,void, entryCompareFn).init(allocator,{});

			try priority_queue.ensureTotalCapacity(K);

			return Self {
				.backing_storage = memory,
				.priority_queue = priority_queue,
				.node_list = nodel,
				.write_ptr = 0,
				.sorted_targets = try std.BoundedArray(TopKQueueEntry,K).init(0),
			};
		}

		pub inline fn deinit(self: *Self, allocator: std.mem.Allocator) void {
			self.priority_queue.deinit();
			allocator.free(self.backing_storage);
			allocator.free(self.node_list);
		}

		pub const TopKScorerIterator = struct {
			bounded_arr: *std.BoundedArray(TopKQueueEntry,K),
			index: u32,
			pub inline fn next(self: *TopKScorerIterator) ?T {
				if(self.index >= self.bounded_arr.len) return null;
				const item = self.bounded_arr.buffer[self.index].node;
				self.index+=1;
				return item;
			}

			pub fn reset(self: *TopKScorerIterator) void {
				self.index = 0;
			}
		};

		pub fn reset(self: *Self) void {
			self.priority_queue.len = 0;
			self.write_ptr = 0;
		}

		pub inline fn addVisit(self: *Self, node: T, exists: bool) void {
			if(!exists) {
				self.backing_storage[node] = 0;
			}
			self.backing_storage[node] += 2;
			self.node_list[self.write_ptr] = node;
			self.write_ptr += @boolToInt(!exists);
		}

		pub inline fn iterator(self: *Self, graph: *Graph(T)) !TopKScorerIterator {
			var index:u32 = 0;
			for(self.node_list[0..self.write_ptr]) |item| {
				// Compact node list
				if(graph.erased_nodes.get(item)) {
					@panic("There should never be erased nodes here!");
				}
				index+=1;
				if(graph.node_list[item].cardinality() > 10000 and !graph.node_list[item].isLargeNode()) {
					try graph.node_list[item].promoteToLargeDegreeNode(graph);
				}
				const item_score = @intCast(i32,graph.node_list[item].cardinality()) - self.backing_storage[item];


				if(self.priority_queue.count() < K) {
					try self.priority_queue.add(
						TopKQueueEntry {
							.node = @intCast(T,item),
							.score = item_score
						}
					);
				}
				else if(self.priority_queue.peek()) |prio_item| {
					if(prio_item.score > item_score) {
						_ = self.priority_queue.remove();
						try self.priority_queue.add(TopKQueueEntry {
							.node = @intCast(T,item),
							.score = item_score
						});
					}
				}
			}

			var iter_queue = self.priority_queue.iterator();
			self.sorted_targets.len = 0;
			while(iter_queue.next()) |item| {
				try self.sorted_targets.append(item);
			}

			std.sort.sort(TopKQueueEntry,self.sorted_targets.buffer[0..self.sorted_targets.len],{}, Self.lessThan);

			return TopKScorerIterator {
				.bounded_arr = &self.sorted_targets,
				.index = 0
			};
		}

		pub inline fn cachedIterator(self: *Self) TopKScorerIterator {
			return TopKScorerIterator {
				.bounded_arr = &self.sorted_targets,
				.index = 0
			};
		}

		pub inline fn unsetVisitedBitset(self: *Self, bs: *set.FastBitSet) void {
			if (bs.cardinality < 200) {
				for(self.node_list[0..self.write_ptr]) |item| {
					_ = bs.unset(item);
				}
			}
			else {
				bs.unsetAll();
			}
			if (bs.cardinality != 0) {
				@panic("Bitset cardinality is not zero!");
			}
		}


		pub inline fn copyResultsToArray(self: *Self, array: []i32, visited: *set.FastBitSet) TopKError!void {
			if(array.len < self.backing_storage.len) {
				return TopKError.ArrayTooSmall;
			}

			var vis = visited.iter();
			while(vis.next()) |item| {
				array[item] = self.backing_storage[item];
			}
		}
	};
}


test "TopKScorer: Check top k scorer basics" {
	var gpa = std.heap.GeneralPurposeAllocator(.{}){};
	var topk = try TopKScorer(u16,1).init(gpa.allocator(),100);
	var graph = try Graph(u16).new(10, gpa.allocator());
	try graph.addEdge(2, 3);
	try graph.addEdge(3, 4);
	topk.reset();

	topk.addVisit(2,false);
	topk.addVisit(3,false);
	topk.addVisit(3,true);

	var topkiter = try topk.iterator(&graph);
	const item = topkiter.next() orelse unreachable;
	try std.testing.expectEqual(item,3);
	try std.testing.expectEqual(topkiter.next(),null);
}
