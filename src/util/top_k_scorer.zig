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

		backing_storage: []T,
		
		node_list: []T,
		write_ptr: u32,

		priority_queue: std.PriorityQueue(TopKQueueEntry,void, entryCompareFn),

		pub fn init(allocator: std.mem.Allocator, total_number_of_nodes: T) !Self {
			var memory = try allocator.alloc(T, total_number_of_nodes);
			std.mem.set(T, memory,0);

			var nodel = try allocator.alloc(T, total_number_of_nodes);
			std.mem.set(T, nodel,0);
			
			var priority_queue = std.PriorityQueue(TopKQueueEntry,void, entryCompareFn).init(allocator,{});

			try priority_queue.ensureTotalCapacity(K);

			return Self {
				.backing_storage = memory,
				.priority_queue = priority_queue,
				.node_list = nodel,
				.write_ptr = 0
			};
		}

		pub inline fn deinit(self: *Self, allocator: std.mem.Allocator) void {
			self.priority_queue.deinit();
			allocator.free(self.backing_storage);
		}

		pub const TopKScorerIterator = struct {
			iterator: std.PriorityQueue(TopKQueueEntry,void, entryCompareFn).Iterator,
			pub inline fn next(self: *TopKScorerIterator) ?T {
				if(self.iterator.next()) |item| {
					return item.node;
				}
				return null;
			}
		};

		pub fn reset(self: *Self) void {
			self.write_ptr = 0;
			while(self.priority_queue.removeOrNull()) |_| {
			}
		}

		pub fn addNewNode(self: *Self, node: T) void {
			self.node_list[self.write_ptr] = node;
			self.write_ptr+=1;
		}

		pub inline fn iterator(self: *Self, graph: *Graph(T)) !TopKScorerIterator {
			for(self.node_list[0..self.write_ptr]) |item| {
				const score = self.backing_storage[item];

				const item_score = @intCast(i32,graph.node_list[item].cardinality())-@intCast(i32,score);
				
				self.backing_storage[item] = 0;

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

			return TopKScorerIterator {
				.iterator = self.priority_queue.iterator()
			};
		}


		pub inline fn addVisit(self: *Self, visited: T) void {
			self.backing_storage[visited]+=2;
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
	topk.setNextTargetDegree(5);

	topk.addVisit(2,2);
	topk.addVisit(3,1);

	var topkiter = topk.iterator();
	const item = topkiter.next() orelse unreachable;
	try std.testing.expectEqual(item,3);
	try std.testing.expectEqual(topkiter.next(),null);
}
