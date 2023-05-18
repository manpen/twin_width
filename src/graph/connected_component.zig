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
        bfs_levels: u32,

        pub fn deinit(self: *Self, allocator: std.mem.Allocator) void {
            self.best_contraction_sequence.deinit(allocator);
						self.current_contraction_seq.deinit(allocator);
        }

        pub fn solveSweepingTopK(self: *Self, comptime K: u32, comptime P: u32, solver: *solver_resources.SolverResources(T,K,P), probing: bool) !T {
					var now = try std.time.Instant.now();
					const time_passed = (now.since(self.subgraph.graph.started_at)/(1000*1000*1000));
					// Graph 194 takes around ~16 Seconds to finish up 2x safety factor
					if(time_passed >= 267) {
						return std.math.maxInt(T);
					}
					const result = try self.subgraph.solveSweepingSolverTopK(K,P,&self.current_contraction_seq,solver,(267-time_passed), probing);
					if(result < self.tww) {
						try self.best_contraction_sequence.copyInto(&self.current_contraction_seq.seq);
						self.tww = result;
					}
					return result;
        }

        pub fn solveGreedyTopK(self: *Self, comptime K: u32, comptime P: u32, solver: *solver_resources.SolverResources(T,K,P)) !T {
					var result = try self.subgraph.solveGreedyTopK(K,P,&self.current_contraction_seq,solver,true);
					try self.best_contraction_sequence.copyInto(&self.current_contraction_seq.seq);
					self.tww = result;
					if(self.tww > 0) return self.tww;

					if(self.tww < 200 or self.subgraph.nodes.len < 3700) {
						if(self.subgraph.nodes.len > 2_000_000 and self.tww >= 90) {}
						else {
							try self.resetGraph();
							solver.reset();
							const sweeping = try self.solveSweepingTopK(K,P,solver, false);
							result = std.math.min(sweeping,result);
						}
					}
					else if(self.tww < 500) {
							try self.resetGraph();
							solver.reset();
							const sweeping = try self.solveSweepingTopK(K,P,solver, true);
							result = std.math.min(sweeping,result);
					}

					return result;
        }

				pub inline fn resetGraph(self: *Self) !void {
					while(self.current_contraction_seq.lastContraction()) |_| {
						try self.subgraph.graph.revertLastContraction(&self.current_contraction_seq);
					}
				}

        pub fn solveGreedy(self: *Self, comptime K: u32, comptime P: u32, solver: *solver_resources.SolverResources(T,K,P)) !T {
					
					var ctx_node = NodeScorerInduced(T){};
					var ctx_score = GraphScorer(T){};

					var final_tww:?T = null;
					for(0..8) |i| {
						const result = try self.subgraph.solveGreedyLookahead(@TypeOf(ctx_node),GraphScorer(T), &ctx_node, &ctx_score, K, @intCast(T,i),&self.current_contraction_seq, &solver.bfs_stack, &solver.scratch_bitset, final_tww);
						self.tww = result;
						if(final_tww==null) {
							final_tww = result;
							try self.best_contraction_sequence.copyInto(&self.current_contraction_seq.seq);
						}
						else {
							if(result < final_tww.?) {
								try self.best_contraction_sequence.copyInto(&self.current_contraction_seq.seq);
							}
							final_tww = std.math.min(result,final_tww.?);
						}

						try self.resetGraph();
					}


					var ctx_min_sim = GraphScorerMinSim(T, .cumulative) {
						.visited = &solver.scratch_bitset,
						.bfs = &solver.bfs_stack,
						.subgraph = &self.subgraph,
						.node_tuples = solver.node_tuple,
					};
					ctx_min_sim.initAll();
					for(0..3) |i| {
						const result = try self.subgraph.solveGreedyLookahead(@TypeOf(ctx_node),@TypeOf(ctx_min_sim), &ctx_node, &ctx_min_sim, K, @intCast(T,i),&self.current_contraction_seq, &solver.bfs_stack, &solver.scratch_bitset, final_tww);
						self.tww = result;
						if(final_tww==null) {
							final_tww = result;
							try self.best_contraction_sequence.copyInto(&self.current_contraction_seq.seq);
						}
						else {
							if(result < final_tww.?) {
								try self.best_contraction_sequence.copyInto(&self.current_contraction_seq.seq);
							}
							final_tww = std.math.min(result,final_tww.?);
						}

						try self.resetGraph();
					}
					var ctx_new_reds = NodeScorerNewReds(T){};

					
					var node_edge_reducer = node_scorers.NodeScorerEdgeReducer(T){
						.seq = &self.current_contraction_seq
					};
					for(0..5) |i| {
						const result = try self.subgraph.solveGreedyLookahead(@TypeOf(node_edge_reducer),@TypeOf(ctx_score), &node_edge_reducer, &ctx_score, K, @intCast(T,i),&self.current_contraction_seq, &solver.bfs_stack, &solver.scratch_bitset, final_tww);
						self.tww = result;
						if(final_tww==null) {
							final_tww = result;
							try self.best_contraction_sequence.copyInto(&self.current_contraction_seq.seq);
						}
						else {
							if(result < final_tww.?) {
								try self.best_contraction_sequence.copyInto(&self.current_contraction_seq.seq);
							}
							final_tww = std.math.min(result,final_tww.?);
						}

						try self.resetGraph();
					}
					

					var ctx_min_sim_min = GraphScorerMinSim(T, .min) {
						.visited = &solver.scratch_bitset,
						.bfs = &solver.bfs_stack,
						.subgraph = &self.subgraph,
						.node_tuples = solver.node_tuple,
					};
					ctx_min_sim_min.initAll();
					for(0..3) |i| {
						const result = try self.subgraph.solveGreedyLookahead(@TypeOf(ctx_node),@TypeOf(ctx_min_sim_min), &ctx_node, &ctx_min_sim_min, K, @intCast(T,i),&self.current_contraction_seq, &solver.bfs_stack, &solver.scratch_bitset, final_tww);
						self.tww = result;
						if(final_tww==null) {
							final_tww = result;
							try self.best_contraction_sequence.copyInto(&self.current_contraction_seq.seq);
						}
						else {
							if(result < final_tww.?) {
								try self.best_contraction_sequence.copyInto(&self.current_contraction_seq.seq);
							}
							final_tww = std.math.min(result,final_tww.?);
						}

						try self.resetGraph();
					}


					var ctx_weighted_scorer_min = GraphScorerWeightedMinor(T){};
					
					var ctx_weighted_scorer = GraphScorerWeightedMajor(T){};
					for(0..4) |i| {
						const result = try self.subgraph.solveGreedyLookahead(@TypeOf(ctx_new_reds),@TypeOf(ctx_score), &ctx_new_reds, &ctx_score, K, @intCast(T,i),&self.current_contraction_seq, &solver.bfs_stack, &solver.scratch_bitset, final_tww);
						self.tww = result;
						if(final_tww==null) {
							final_tww = result;
							try self.best_contraction_sequence.copyInto(&self.current_contraction_seq.seq);
						}
						else {
							if(result < final_tww.?) {
								try self.best_contraction_sequence.copyInto(&self.current_contraction_seq.seq);
							}
							final_tww = std.math.min(result,final_tww.?);
						}

						try self.resetGraph();
					}

					for(0..2) |i| {
						const result = try self.subgraph.solveGreedyLookahead(@TypeOf(ctx_node),@TypeOf(ctx_weighted_scorer_min), &ctx_node, &ctx_weighted_scorer_min, K, @intCast(T,i),&self.current_contraction_seq, &solver.bfs_stack, &solver.scratch_bitset, final_tww);
						self.tww = result;
						if(final_tww==null) {
							final_tww = result;
							try self.best_contraction_sequence.copyInto(&self.current_contraction_seq.seq);
						}
						else {
							if(result < final_tww.?) {
								try self.best_contraction_sequence.copyInto(&self.current_contraction_seq.seq);
							}
							final_tww = std.math.min(result,final_tww.?);
						}

						try self.resetGraph();
					}
					
					for(0..2) |i| {
						const result = try self.subgraph.solveGreedyLookahead(@TypeOf(ctx_node),@TypeOf(ctx_weighted_scorer), &ctx_node, &ctx_weighted_scorer, K, @intCast(T,i),&self.current_contraction_seq, &solver.bfs_stack, &solver.scratch_bitset, final_tww);
						self.tww = result;
						if(final_tww==null) {
							try self.best_contraction_sequence.copyInto(&self.current_contraction_seq.seq);
							final_tww = result;
						}
						else {
							if(result < final_tww.?) {
								try self.best_contraction_sequence.copyInto(&self.current_contraction_seq.seq);
							}
							final_tww = std.math.min(result,final_tww.?);
						}

						try self.resetGraph();
					}

					var selector_simple = NodeScorerSimple(T){
					};
					for(0..2) |i| {
						const result = try self.subgraph.solveGreedyLookahead(@TypeOf(selector_simple),@TypeOf(ctx_weighted_scorer), &selector_simple, &ctx_weighted_scorer, K, @intCast(T,i),&self.current_contraction_seq, &solver.bfs_stack, &solver.scratch_bitset, final_tww);
						self.tww = result;
						if(final_tww==null) {
							try self.best_contraction_sequence.copyInto(&self.current_contraction_seq.seq);
							final_tww = result;
						}
						else {
							if(result < final_tww.?) {
								try self.best_contraction_sequence.copyInto(&self.current_contraction_seq.seq);
							}
							final_tww = std.math.min(result,final_tww.?);
						}

						try self.resetGraph();
					}


					return final_tww.?;
        }

        pub fn init(allocator: std.mem.Allocator, nodes: []T, bfs_level: u32, graph: *graph_mod.Graph(T)) !Self {
            var contraction_seq = try contraction.ContractionSequence(T).init(allocator, @intCast(u32,nodes.len));
            var previous: ?T = null;
						var number_of_edges:u32 = 0;

            for (nodes) |item| {
								number_of_edges += graph.node_list[item].cardinality();
                if (previous) |into| {
                    try contraction_seq.addContraction(contraction.Contraction(T){ .erased = item, .survivor = into });
                } else {
                    previous = item;
                }
            }
						number_of_edges = number_of_edges>>1;
            var retraceable = try retraceable_contraction.RetraceableContractionSequence(T).init(graph.allocator, @intCast(T,nodes.len), number_of_edges);
            return .{ 
						.bfs_levels = bfs_level, .subgraph = InducedSubGraph(T).fromSlice(graph,nodes), .tww = @intCast(T, nodes.len - 1), .best_contraction_sequence = contraction_seq,
						.current_contraction_seq = retraceable
						};
        }
    };
}
