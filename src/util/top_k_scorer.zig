const std = @import("std");
const comptime_util = @import("comptime_checks.zig");

pub fn TopKScorer(comptime T: type, comptime K: u32) type {
	comptime if (!comptime_util.checkIfIsCompatibleInteger(T)) {
		@compileError("Type must be either u8,u16 or u32");
	};

	return struct {
		const Self = @This();
		pub const TopKScorerEntry  = packed struct {
			unique_timestamp_of_iteration: u32,
			score: T
		};

		pub const TopKScorerContext = struct {
			target_degree: T,
			memory: []TopKScorerEntry,
		};

		pub fn entryCompareFn(context: *TopKScorerContext, lhs: T,rhs: T) std.math.Order {
			return std.math.order(context.memory[rhs].score,context.memory[lhs].score);
		}

		backing_storage: []TopKScorerEntry,
		unique_id: u32,
		top_k_context: *TopKScorerContext,
		priority_queue: std.PriorityQueue(T,*TopKScorerContext, entryCompareFn),
		in_prio: std.bit_set.DynamicBitSetUnmanaged,
		add_visit_calls: u32,

		pub fn init(allocator: std.mem.Allocator, total_number_of_nodes: T) !Self {
			var memory = try allocator.alloc(TopKScorerEntry, total_number_of_nodes);
			for(0..total_number_of_nodes) |i| {
				memory[i].unique_timestamp_of_iteration = 0;
			}
			var top_k_entry = try allocator.create(TopKScorerContext);
			top_k_entry.memory = memory;
			top_k_entry.target_degree = 0;

			var priority_queue = std.PriorityQueue(T,*TopKScorerContext, entryCompareFn).init(allocator,top_k_entry);

			try priority_queue.ensureTotalCapacity(K);

			return Self {
				.backing_storage = memory,
				.unique_id = 1,
				.priority_queue = priority_queue,
				.top_k_context = top_k_entry,
				.in_prio = try std.bit_set.DynamicBitSetUnmanaged.initEmpty(allocator,total_number_of_nodes),
				.add_visit_calls = 0,
			};
		}

		pub inline fn deinit(self: *Self, allocator: std.mem.Allocator) void {
			self.priority_queue.deinit();
			allocator.destroy(self.top_k_context);
			allocator.free(self.backing_storage);
		}

		inline fn increaseUniqueId(self: *Self) void {
			self.unique_id+=1;
			self.add_visit_calls = 0;

			// Clear the queue
			while(self.priority_queue.removeOrNull()) |item| {
				self.in_prio.unset(item);
			}
		}

		pub inline fn setNextTargetDegree(self: *Self, target_degree: T) void {
			self.increaseUniqueId();
			self.top_k_context.target_degree = target_degree;
		}

		pub const TopKScorerIterator = struct {
			scorer: *Self,
			pub inline fn next(self: *TopKScorerIterator) ?T {
				if(self.scorer.priority_queue.removeOrNull()) |item| {
					self.scorer.in_prio.unset(item);
					return item;
				}
				return null;
			}
		};

		pub inline fn iterator(self: *Self) TopKScorerIterator {
			return TopKScorerIterator {
				.scorer = self
			};
		}

		pub inline fn addVisit(self: *Self, visited: T, total_degree: T) void {
			self.add_visit_calls+=1;
			if(self.backing_storage[visited].unique_timestamp_of_iteration != self.unique_id) {
				self.backing_storage[visited].unique_timestamp_of_iteration = self.unique_id;
				const new_score = (self.top_k_context.target_degree-1)+(total_degree-1);
				self.backing_storage[visited].score = new_score;
			}
			else {
				self.backing_storage[visited].score-=2;
			}

			
			if(self.in_prio.isSet(visited)) return;

			if(self.priority_queue.count() < K) {
				self.priority_queue.add(visited) catch unreachable;
				self.in_prio.set(visited);
			}
			else {
				const highest_prio = self.priority_queue.peek().?;
				if(self.backing_storage[highest_prio].score > self.backing_storage[visited].score) {
					self.priority_queue.update(highest_prio,visited) catch unreachable;
					self.in_prio.unset(highest_prio);
					self.in_prio.set(visited);
				}
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
