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
			score: i32
		};

		pub const TopKBfsIncreasor  = packed struct {
			unique_timestamp_of_iteration: u32,
			degree_added: u64,
		};

		pub const TopKScorerContext = struct {
			target_degree: T,
			memory: []TopKScorerEntry,
		};

		pub const TopKQueueEntry = struct {
			node: T,
			score: i32
		};

		pub fn entryCompareFn(context: *TopKScorerContext, lhs: TopKQueueEntry,rhs: TopKQueueEntry) std.math.Order {
			_ = context;
			return std.math.order(rhs.score,lhs.score);
		}

		backing_storage: []TopKScorerEntry,
		unique_id: u32,
		top_k_context: *TopKScorerContext,
		priority_queue: std.PriorityQueue(TopKQueueEntry,*TopKScorerContext, entryCompareFn),
		in_prio: std.bit_set.DynamicBitSetUnmanaged,
		add_visit_calls: u32,

		largest_bfs_increasor: []TopKBfsIncreasor,
		bfs_unique_id: u32,
		largest_bfs_node_id: T,

		pub fn init(allocator: std.mem.Allocator, total_number_of_nodes: T) !Self {
			var memory = try allocator.alloc(TopKScorerEntry, total_number_of_nodes);
			for(0..total_number_of_nodes) |i| {
				memory[i].unique_timestamp_of_iteration = 0;
			}
			var bfsinc = try allocator.alloc(TopKBfsIncreasor, total_number_of_nodes);
			for(0..total_number_of_nodes) |i| {
				bfsinc[i].unique_timestamp_of_iteration = 0;
			}
			bfsinc[0].degree_added = 0;
			var top_k_entry = try allocator.create(TopKScorerContext);
			top_k_entry.memory = memory;
			top_k_entry.target_degree = 0;

			var priority_queue = std.PriorityQueue(TopKQueueEntry,*TopKScorerContext, entryCompareFn).init(allocator,top_k_entry);

			try priority_queue.ensureTotalCapacity(K);

			return Self {
				.backing_storage = memory,
				.largest_bfs_increasor = bfsinc,
				.unique_id = 1,
				.bfs_unique_id = 1,
				.largest_bfs_node_id = 0,
				.priority_queue = priority_queue,
				.top_k_context = top_k_entry,
				.in_prio = try std.bit_set.DynamicBitSetUnmanaged.initEmpty(allocator,total_number_of_nodes),
				.add_visit_calls = 0,
			};
		}

		pub inline fn deinit(self: *Self, allocator: std.mem.Allocator) void {
			self.in_prio.deinit(allocator);
			self.priority_queue.deinit();
			allocator.destroy(self.top_k_context);
			allocator.free(self.largest_bfs_increasor);
			allocator.free(self.backing_storage);
		}

		inline fn increaseUniqueId(self: *Self) void {
			self.unique_id+=1;
			self.add_visit_calls = 0;

			// Clear the queue
			while(self.priority_queue.removeOrNull()) |item| {
				self.in_prio.unset(item.node);
			}
		}

		pub inline fn setNextTargetDegree(self: *Self, target_degree: T) void {
			self.increaseUniqueId();
			self.top_k_context.target_degree = target_degree;
		}

		pub const TopKScorerIterator = struct {
			iterator: std.PriorityQueue(TopKQueueEntry,*TopKScorerContext, entryCompareFn).Iterator,
			pub inline fn next(self: *TopKScorerIterator) ?T {
				if(self.iterator.next()) |item| {
					return item.node;
				}
				return null;
			}
		};

		pub inline fn iterator(self: *Self) TopKScorerIterator {
			return TopKScorerIterator {
				.iterator = self.priority_queue.iterator()
			};
		}

		pub inline fn consumeLargestBfsIncreasor(self: *Self) T {
			self.bfs_unique_id += 1;
			self.largest_bfs_increasor[self.largest_bfs_node_id].degree_added = 0;
			return self.largest_bfs_node_id;
		}

		pub inline fn addBfsIncreasor(self: *Self, visited: T, total_degree: T) void {
			if(self.largest_bfs_increasor[visited].unique_timestamp_of_iteration != self.bfs_unique_id) {
				self.largest_bfs_increasor[visited].unique_timestamp_of_iteration = self.bfs_unique_id;
				self.largest_bfs_increasor[visited].degree_added = total_degree;
			}
			else {
				self.largest_bfs_increasor[visited].degree_added += total_degree;
			}

			if(self.largest_bfs_increasor[visited].degree_added > self.largest_bfs_increasor[self.largest_bfs_node_id].degree_added) {
				self.largest_bfs_node_id = visited;
			}
		}

		pub inline fn addVisit(self: *Self, visited: T, total_degree: T) void {
			self.add_visit_calls+=1;
			if(self.backing_storage[visited].unique_timestamp_of_iteration != self.unique_id) {
				self.backing_storage[visited].unique_timestamp_of_iteration = self.unique_id;
				//const new_score = (self.top_k_context.target_degree-1)+(total_degree-1);
				self.backing_storage[visited].score = @intCast(i32,total_degree)-2;
			}
			else {
				self.backing_storage[visited].score-=2;
			}

			
			if(self.in_prio.isSet(visited)) return;

			if(self.priority_queue.count() < K) {
				self.priority_queue.add(TopKQueueEntry {
					.node = visited,
					.score = self.backing_storage[visited].score
				}) catch unreachable;
				self.in_prio.set(visited);
			}
			else {
				var score_now_highest:i32 = undefined;
				while(true) {
					const highest_prio = self.priority_queue.peek().?;
					score_now_highest = self.backing_storage[highest_prio.node].score;
					if(score_now_highest != highest_prio.score) {
						var hp = self.priority_queue.remove();
						hp.score = score_now_highest;
						self.priority_queue.add(hp) catch unreachable;
					}
					else {
						break;
					}
				}
				if(score_now_highest > self.backing_storage[visited].score) {
					const removed = self.priority_queue.remove();
					self.priority_queue.add(TopKQueueEntry {
						.node = visited,
						.score = self.backing_storage[visited].score
					}) catch unreachable;
					self.in_prio.unset(removed.node);
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