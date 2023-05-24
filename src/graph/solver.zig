const bitset = @import("../util/two_level_bitset.zig");
const bfs_mod = @import("bfs.zig");
const topk = @import("../util/top_k_scorer.zig");
const comptime_util = @import("../util/comptime_checks.zig");
const node_priority_queue = @import("node_priority_queue.zig");
const Graph = @import("graph.zig").Graph;
const std = @import("std");
const dfs_mod = @import("dfs.zig");

pub fn NodePairs(comptime T: type) type {
    return struct {
        first: T,
        second: T,
        potential: Graph(T).InducedTwinWidth,

        pub fn sortPairs(ctx: T, lhs: NodePairs(T), rhs: NodePairs(T)) bool {
            return lhs.potential.isLess(rhs.potential, ctx);
        }

        pub fn compare(ctx: T, lhs: NodePairs(T), rhs: NodePairs(T)) std.math.Order {
            _ = ctx;
            if (lhs.potential.isLess(rhs.potential)) {
                return std.math.Order.gt;
            } else {
                return std.math.Order.lt;
            }
        }
    };
}

pub fn NodeTuple(comptime T: type) type {
    _ = T;
    return struct { first: u32, second: u32 };
}

pub fn ArticulationPointResources(comptime T: type) type {
    if (!comptime_util.checkIfIsCompatibleInteger(T)) {
        @compileError("T must be an integer type u8,u16 or u32!\n");
    }

    return struct {
        const Self = @This();
        parent_array: []T,
        depth_array: []T,
        low_array: []T,
        node_count: []T,

        articulation_points: bitset.FastBitSet,

        node_stack: std.ArrayListUnmanaged(dfs_mod.NodeStackEntry(T)),

        pub fn init(graph: *Graph(T)) !Self {
            return Self{
                .parent_array = try graph.allocator.alloc(T, graph.number_of_nodes),
                .depth_array = try graph.allocator.alloc(T, graph.number_of_nodes),
                .low_array = try graph.allocator.alloc(T, graph.number_of_nodes),
                .node_count = try graph.allocator.alloc(T, graph.number_of_nodes),
                .articulation_points = try bitset.FastBitSet.initEmpty(graph.number_of_nodes, graph.allocator),
                .node_stack = try std.ArrayListUnmanaged(dfs_mod.NodeStackEntry(T)).initCapacity(graph.allocator, graph.number_of_edges * 2),
            };
        }

        pub fn deinit(self: *Self, allocator: std.mem.Allocator) void {
            allocator.free(self.parent_array);
            allocator.free(self.depth_array);
            allocator.free(self.low_array);
            allocator.free(self.node_count);
            self.node_stack.deinit(allocator);
            self.articulation_points.deinit(allocator);
        }
    };
}

pub fn SolverResources(comptime T: type, comptime K: u32, comptime P: u32) type {
    if (!comptime_util.checkIfIsCompatibleInteger(T)) {
        @compileError("T must be an integer type u8,u16 or u32!\n");
    }
    return struct {
        const Self = @This();
        bfs_stack: bfs_mod.BfsQueue(T),
        scratch_bitset: bitset.FastBitSet,
        scorer: topk.TopKScorer(T, K),
        priority_queue: node_priority_queue.NodePriorityQueue(T),

        potential_ordered_nodes: std.PriorityQueue(NodePairs(T), T, NodePairs(T).compare),
        node_tuple: []NodeTuple(T),

        node_mask_bitset: bitset.FastBitSet,
        scratch_node_list: []T,

        articulation_point_resources: ArticulationPointResources(T),

        pub inline fn init(graph: *Graph(T)) !Self {
            var bfs = try bfs_mod.BfsQueue(T).init(graph.allocator, graph.number_of_nodes);

            var scratch_bitset = try bitset.FastBitSet.initEmpty(graph.number_of_nodes, graph.allocator);

            var node_mask_bitset = try bitset.FastBitSet.initEmpty(graph.number_of_nodes, graph.allocator);

            var scratch_node_list = try graph.allocator.alloc(T, graph.number_of_nodes);

            var scorer = try topk.TopKScorer(T, K).init(graph.allocator, @intCast(T, graph.number_of_nodes));

            var priority_queue = try node_priority_queue.NodePriorityQueue(T).init(graph, @intCast(T, P));

            var potential_ordered_node = std.PriorityQueue(NodePairs(T), T, NodePairs(T).compare).init(graph.allocator, 0);

            try potential_ordered_node.ensureTotalCapacity(K);

            var node_tuple = try graph.allocator.alloc(NodeTuple(T), graph.number_of_nodes);

            var art = try ArticulationPointResources(T).init(graph);

            return Self{
                .bfs_stack = bfs,
                .scratch_bitset = scratch_bitset,
                .scratch_node_list = scratch_node_list,
                .node_mask_bitset = node_mask_bitset,
                .scorer = scorer,
                .priority_queue = priority_queue,
                .articulation_point_resources = art,
                .potential_ordered_nodes = potential_ordered_node,
                .node_tuple = node_tuple,
            };
        }

        pub inline fn reset(self: *Self) void {
            self.scratch_bitset.unsetAll();
            self.bfs_stack.clear();
            self.priority_queue.clear();
            self.node_mask_bitset.unsetAll();
            self.potential_ordered_nodes.context = 0;
            self.potential_ordered_nodes.len = 0;
            //Scorer does not need clearing
        }

        pub inline fn deinit(self: *Self, allocator: std.mem.Allocator) void {
            self.priority_queue.deinit();
            self.scorer.deinit(allocator);
            self.scratch_bitset.deinit(allocator);
            self.bfs_stack.deinit(allocator);
            allocator.free(self.scratch_node_list);
            self.articulation_point_resources.deinit(allocator);
            self.node_mask_bitset.deinit(allocator);
            self.potential_ordered_nodes.deinit();
            allocator.free(self.node_tuple);
        }
    };
}
