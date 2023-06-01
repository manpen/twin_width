const std = @import("std");
const edge_list = @import("edge_list.zig");
const contraction = @import("../tww/contraction_sequence.zig");
const graph_mod = @import("graph.zig");
const comptime_util = @import("../util/comptime_checks.zig");
const bfs_mod = @import("bfs.zig");
const CompactField = graph_mod.CompactField;
const bench_timer = @import("../util/benchmark_helper.zig");
const topk = @import("../util/top_k_scorer.zig");
const InducedSubGraph = @import("subgraph.zig").InducedSubGraph;
const solver_resources = @import("solver.zig");
const retraceable_contraction = @import("../tww/retraceable_contraction_sequence.zig");
const Graph = @import("graph.zig").Graph;

const NodeScorerPotential = @import("../scorer/node_scorer.zig").NodeScorerPotential;
const NodeScorerSimple = @import("../scorer/node_scorer.zig").NodeScorerSimple;
const NodeScorerInduced = @import("../scorer/node_scorer.zig").NodeScorerInduced;
const NodeScorerNewReds = @import("../scorer/node_scorer.zig").NodeScorerNewReds;
const NodeScorerTotalRedDeg = @import("../scorer/node_scorer.zig").NodeScorerTotalRedDeg;
const GraphScorer = @import("../scorer/graph_scorer.zig").GraphScorer;
const GraphScorerWeightedMajor = @import("../scorer/graph_scorer.zig").GraphScorerWeightedMajor;
const GraphScorerWeightedMinor = @import("../scorer/graph_scorer.zig").GraphScorerWeightedMinor;
const GraphScorerMinSim = @import("../scorer/graph_scorer.zig").GraphScorerMinSim;
const GraphScorerMinor = @import("../scorer/graph_scorer.zig").GraphScorerMinor;
const GraphScorerMinBlacks = @import("../scorer/graph_scorer.zig").GraphScorerMinBlack;
const graph_scorers = @import("../scorer/graph_scorer.zig");
const node_scorers = @import("../scorer/node_scorer.zig");

pub fn ConnectedComponentIndex(comptime T: type) type {
    comptime if (!comptime_util.checkIfIsCompatibleInteger(T)) {
        @compileError("The type of ConnectedComponent must either be u8,u16 or u32!");
    };
    return struct {
        const Self = @This();
        tww: T,
        index: T,
        pub fn compareComponentIndexDesc(ctx: void, lhs: Self, rhs: Self) std.math.Order {
            _ = ctx;
            return std.math.order(rhs.tww, lhs.tww);
        }
    };
}

pub fn ConnectedComponent(comptime T: type) type {
    comptime if (!comptime_util.checkIfIsCompatibleInteger(T)) {
        @compileError("The type of ConnectedComponent must either be u8,u16 or u32!");
    };

    return struct {
        const Self = @This();
        subgraph: InducedSubGraph(T),
        best_contraction_sequence: contraction.ContractionSequence(T),
        current_contraction_seq: retraceable_contraction.RetraceableContractionSequence(T),
        tww: T,
        iteration: u32,
        contraction_slice: []u32,
        backup_contraction_slice: []u32,

        pub fn deinit(self: *Self, allocator: std.mem.Allocator) void {
            allocator.free(self.backup_contraction_slice);
            allocator.free(self.contraction_slice);
            self.best_contraction_sequence.deinit(allocator);
            self.current_contraction_seq.deinit(allocator);
        }

        pub fn solveSweepingTopK(self: *Self, comptime K: u32, comptime P: u32, solver: *solver_resources.SolverResources(T, K, P), probing: bool, seed: u64) !T {
            const result = try self.subgraph.solveSweepingSolverTopK(K, P, &self.current_contraction_seq, solver, probing, seed);
            if (result < self.tww) {
                try self.update_best();
            }
            return result;
        }

        pub fn solveSweepingTopKPrecision(self: *Self, comptime K: u32, comptime P: u32, solver: *solver_resources.SolverResources(T, K, P), seed: u64) !T {
            const result = try self.subgraph.solveSweepingSolverTopKPrecision(K, P, &self.current_contraction_seq, solver, self.iteration, seed);
            if (result < self.tww) {
                try self.update_best();
            }
            return result;
        }

        fn update_best(self: *Self) !void {
            self.tww = self.current_contraction_seq.getTwinWidth();
            self.current_contraction_seq.seq.write_to_slice(self.backup_contraction_slice);
            // swap ensure that the global slice is updated by a simple redirection of a pointer.
            // TECHNICALLY this is not atomic, as it consists of two assignments (length and pointer).
            // If this is a problem, we should use a mutex. should add negligible overhead
            var tmp = self.contraction_slice;
            self.contraction_slice = self.backup_contraction_slice;
            self.backup_contraction_slice = tmp;
            try self.best_contraction_sequence.copyInto(&self.current_contraction_seq.seq);
        }

        pub fn solveGreedyTopK(self: *Self, comptime K: u32, comptime P: u32, solver: *solver_resources.SolverResources(T, K, P), seed: u64) !T {
            var result: T = undefined;
            if (self.iteration == 0) {
                result = try self.subgraph.solveGreedyTopK(K, P, &self.current_contraction_seq, solver, true, seed, self.tww);

                if (result < self.tww) {
                    try self.update_best();
                }

                if (self.tww < 200 or self.subgraph.nodes.len < 10_000) {
                    try self.resetGraph();
                    solver.reset();
                    const sweeping = try self.solveSweepingTopK(K, P, solver, false, seed);
                    result = std.math.min(sweeping, result);
                    try self.resetGraph();
                } else if (self.tww < 500) {
                    try self.resetGraph();
                    solver.reset();
                    const sweeping = try self.solveSweepingTopK(K, P, solver, true, seed);
                    result = std.math.min(sweeping, result);
                }
            } else {
                if (self.tww < 200 or self.subgraph.nodes.len < 10_000) {
                    try self.resetGraph();
                    solver.reset();
                    const sweeping = try self.solveSweepingTopKPrecision(K, P, solver, seed);
                    result = std.math.min(sweeping, result);
                } else if (self.subgraph.nodes.len >= 500) {
                    try self.resetGraph();
                    solver.reset();
                    result = try self.subgraph.solveGreedyTopKMinHash(K, P, &self.current_contraction_seq, solver, self.iteration, seed, self.tww);
                    if (result < self.tww) {
                        try self.update_best();
                    }
                }
            }
						
						const time = try std.time.Instant.now();
						std.debug.print("Finished iteration with tww {} and time {}ms\n",.{self.tww, time.since(self.subgraph.graph.started_at)/(1000*1000)});

            self.iteration += 1;
            return result;
        }

        pub inline fn resetGraph(self: *Self) !void {
            while (self.current_contraction_seq.lastContraction()) |_| {
                try self.subgraph.graph.revertLastContraction(&self.current_contraction_seq);
            }
        }

        pub fn solveGreedy(self: *Self, comptime K: u32, comptime P: u32, solver: *solver_resources.SolverResources(T, K, P), seed: u64) !T {
            var ctx_node = NodeScorerInduced(T){};
            var ctx_score = GraphScorer(T){};
            self.iteration += 1;

            var final_tww: T = self.tww;
            for (0..8) |i| {
                const result = try self.subgraph.solveGreedyLookahead(@TypeOf(ctx_node), GraphScorer(T), &ctx_node, &ctx_score, K, @intCast(T, i), &self.current_contraction_seq, &solver.bfs_stack, &solver.scratch_bitset, final_tww, seed);
                if (result < final_tww) {
                    try self.update_best();
                }
                final_tww = std.math.min(result, final_tww);

                try self.resetGraph();
            }

            var ctx_min_sim = GraphScorerMinSim(T, .cumulative){
                .visited = &solver.scratch_bitset,
                .bfs = &solver.bfs_stack,
                .subgraph = &self.subgraph,
                .node_tuples = solver.node_tuple,
            };
            ctx_min_sim.initAll();
            for (0..3) |i| {
                const result = try self.subgraph.solveGreedyLookahead(@TypeOf(ctx_node), @TypeOf(ctx_min_sim), &ctx_node, &ctx_min_sim, K, @intCast(T, i), &self.current_contraction_seq, &solver.bfs_stack, &solver.scratch_bitset, final_tww, seed);
                if (result < final_tww) {
                    try self.update_best();
                }
                final_tww = std.math.min(result, final_tww);

                try self.resetGraph();
            }
            var ctx_new_reds = NodeScorerNewReds(T){};

            var node_edge_reducer = node_scorers.NodeScorerEdgeReducer(T){ .seq = &self.current_contraction_seq };
            for (0..5) |i| {
                const result = try self.subgraph.solveGreedyLookahead(@TypeOf(node_edge_reducer), @TypeOf(ctx_score), &node_edge_reducer, &ctx_score, K, @intCast(T, i), &self.current_contraction_seq, &solver.bfs_stack, &solver.scratch_bitset, final_tww, seed);
                if (result < final_tww) {
                    try self.update_best();
                }
                final_tww = std.math.min(result, final_tww);

                try self.resetGraph();
            }

            var ctx_min_sim_min = GraphScorerMinSim(T, .min){
                .visited = &solver.scratch_bitset,
                .bfs = &solver.bfs_stack,
                .subgraph = &self.subgraph,
                .node_tuples = solver.node_tuple,
            };
            ctx_min_sim_min.initAll();
            for (0..3) |i| {
                const result = try self.subgraph.solveGreedyLookahead(@TypeOf(ctx_node), @TypeOf(ctx_min_sim_min), &ctx_node, &ctx_min_sim_min, K, @intCast(T, i), &self.current_contraction_seq, &solver.bfs_stack, &solver.scratch_bitset, final_tww, seed);
                if (result < final_tww) {
                    try self.update_best();
                }
                final_tww = std.math.min(result, final_tww);

                try self.resetGraph();
            }

            var ctx_weighted_scorer_min = GraphScorerWeightedMinor(T){};

            var ctx_weighted_scorer = GraphScorerWeightedMajor(T){};
            for (0..4) |i| {
                const result = try self.subgraph.solveGreedyLookahead(@TypeOf(ctx_new_reds), @TypeOf(ctx_score), &ctx_new_reds, &ctx_score, K, @intCast(T, i), &self.current_contraction_seq, &solver.bfs_stack, &solver.scratch_bitset, final_tww, seed);
                if (result < final_tww) {
                    try self.update_best();
                }
                final_tww = std.math.min(result, final_tww);

                try self.resetGraph();
            }

            for (0..2) |i| {
                const result = try self.subgraph.solveGreedyLookahead(@TypeOf(ctx_node), @TypeOf(ctx_weighted_scorer_min), &ctx_node, &ctx_weighted_scorer_min, K, @intCast(T, i), &self.current_contraction_seq, &solver.bfs_stack, &solver.scratch_bitset, final_tww, seed);
                if (result < final_tww) {
                    try self.update_best();
                }
                final_tww = std.math.min(result, final_tww);

                try self.resetGraph();
            }

            for (0..2) |i| {
                const result = try self.subgraph.solveGreedyLookahead(@TypeOf(ctx_node), @TypeOf(ctx_weighted_scorer), &ctx_node, &ctx_weighted_scorer, K, @intCast(T, i), &self.current_contraction_seq, &solver.bfs_stack, &solver.scratch_bitset, final_tww, seed);
                if (result < final_tww) {
                    try self.update_best();
                }
                final_tww = std.math.min(result, final_tww);

                try self.resetGraph();
            }

            var selector_simple = NodeScorerSimple(T){};
            for (0..2) |i| {
                const result = try self.subgraph.solveGreedyLookahead(@TypeOf(selector_simple), @TypeOf(ctx_weighted_scorer), &selector_simple, &ctx_weighted_scorer, K, @intCast(T, i), &self.current_contraction_seq, &solver.bfs_stack, &solver.scratch_bitset, final_tww, seed);
                if (result < final_tww) {
                    try self.update_best();
                }
                final_tww = std.math.min(result, final_tww);

                try self.resetGraph();
            }

            return final_tww;
        }

        pub fn init(allocator: std.mem.Allocator, nodes: []T, graph: *graph_mod.Graph(T)) !Self {
            var contraction_seq = try contraction.ContractionSequence(T).init(allocator, @intCast(u32, nodes.len));
            var previous: ?T = null;
            var number_of_edges: u32 = 0;

            for (nodes) |item| {
                number_of_edges += graph.node_list[item].cardinality();
                if (previous) |into| {
                    try contraction_seq.addContraction(contraction.Contraction(T){ .erased = item, .survivor = into });
                } else {
                    previous = item;
                }
            }

            // initialize output solution slice
            var num_entries = nodes.len;
            if (num_entries > 2) {
                // 1 =>
                // 1
                // 1 2 =>
                // 1 2
                // 1 2 3 4 =>
                // 1 2, 2 3, 3 4
                // 1 2 3 4 5 6 7 =>
                // 1 2, 2 3, 3 4, 4 5, 5 6, 6 7
                num_entries = (num_entries - 1) * 2;
            }
            var slice: []u32 = try allocator.alloc(u32, num_entries);
            var backup_slice: []u32 = try allocator.alloc(u32, num_entries);
            if (nodes.len > 2) {
                var j: usize = 0;
                for (0..(nodes.len - 1)) |i| {
                    slice[j] = @as(u32, nodes[i + 1]);
                    slice[j + 1] = @as(u32, nodes[i]);
                    j += 2;
                }
            } else {
                slice[0] = nodes[0];
            }

            number_of_edges = number_of_edges >> 1;
            var retraceable = try retraceable_contraction.RetraceableContractionSequence(T).init(graph.allocator, @intCast(T, nodes.len), number_of_edges);
            return .{
                .iteration = 0,
                .subgraph = InducedSubGraph(T).fromSlice(graph, nodes),
                .tww = @intCast(T, nodes.len - 1),
                .best_contraction_sequence = contraction_seq,
                .current_contraction_seq = retraceable,
                .contraction_slice = slice,
                .backup_contraction_slice = backup_slice,
            };
        }
    };
}
