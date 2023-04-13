const bitset = @import("../util/two_level_bitset.zig");
const bfs_mod = @import("bfs.zig");
const topk = @import("../util/top_k_scorer.zig");
const comptime_util = @import("../util/comptime_checks.zig");
const node_priority_queue = @import("node_priority_queue.zig");
const Graph = @import("graph.zig").Graph;
const std = @import("std");

pub fn SolverResources(comptime T: type, comptime K: u32, comptime P: u32) type {
		if(!comptime_util.checkIfIsCompatibleInteger(T)) {
			@compileError("T must be an integer type u8,u16 or u32!\n");
		}
    return struct {
				const Self = @This();
        bfs_stack: bfs_mod.BfsQueue(T),
        scratch_bitset: bitset.FastBitSet,
				scorer: topk.TopKScorer(T,K),
				priority_queue: node_priority_queue.NodePriorityQueue(T),

				pub inline fn init(graph: *Graph(T)) !Self {
					var bfs = try bfs_mod.BfsQueue(T).init(graph.allocator,graph.number_of_nodes);
					
					var scratch_bitset = try bitset.FastBitSet.initEmpty(graph.number_of_nodes, graph.allocator);
					
					var scorer = try topk.TopKScorer(T, K).init(graph.allocator, @intCast(T, graph.number_of_nodes));

					var priority_queue = try node_priority_queue.NodePriorityQueue(T).init(graph,@intCast(T,P));
					

					return Self {
						.bfs_stack = bfs,
						.scratch_bitset = scratch_bitset,
						.scorer = scorer,
						.priority_queue = priority_queue,
					};
				}

				pub inline fn reset(self: *Self) void {
					self.scratch_bitset.unsetAll();
					self.bfs_stack.clear();
					self.priority_queue.clear();
					//Scorer does not need clearing
				}

				pub inline fn deinit(self: *Self, allocator: std.mem.Allocator) void {
					self.priority_queue.deinit();
					self.scorer.deinit(allocator);
					self.scratch_bitset.deinit(allocator);
					self.bfs_stack.deinit(allocator);
				}
    };
}
