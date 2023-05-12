const std = @import("std");
const Graph = @import("../graph/graph.zig").Graph;
const retraceable_contraction = @import("../tww/retraceable_contraction_sequence.zig");
const bitset = @import("../util/two_level_bitset.zig");
const bfs_mod = @import("../graph/bfs.zig");
const subgraph = @import("../graph/subgraph.zig");
const NodeTuple = @import("../graph/solver.zig").NodeTuple;

pub fn GraphScorer(comptime T: type) type {
    return struct {
        upper_bound: T = 0,
        pub const ScoreType = struct { T, u32 };

        pub fn setUpperBound(own: *@This(), bound: T) void {
            own.upper_bound = bound;
        }

        pub fn compare(ctx: *@This(), first: ScoreType, second: ScoreType) std.math.Order {
            _ = ctx;
            if (first.@"0" == second.@"0") {
                return std.math.order(second.@"1", first.@"1");
            }
            // Gt if tww_first < tww_second
            return std.math.order(second.@"0", first.@"0");
        }

        pub fn evaluate(ctx: @This(), graph: *Graph(T), seq: *retraceable_contraction.RetraceableContractionSequence(T)) ?ScoreType {
            if (seq.getTwinWidth() > ctx.upper_bound) return null;
            var total_reds: u32 = 0;
            for (graph.node_list) |node| {
                total_reds += node.red_edges.cardinality();
            }
            return .{ @intCast(T, seq.getTwinWidth()), total_reds };
        }
    };
}

pub fn GraphScorerMaxRed(comptime T: type) type {
    return struct {
        upper_bound: T = 0,
        pub const ScoreType = struct { T, u32 };

        pub fn setUpperBound(own: *@This(), bound: T) void {
            own.upper_bound = bound;
        }

        pub fn compare(ctx: *@This(), first: ScoreType, second: ScoreType) std.math.Order {
            if (first.@"0" <= ctx.upper_bound and second.@"0" <= ctx.upper_bound) {
                return std.math.order(first.@"1", second.@"1");
            }
            return std.math.order(second.@"0", first.@"0");
        }

        pub fn evaluate(ctx: @This(), graph: *Graph(T), seq: *retraceable_contraction.RetraceableContractionSequence(T)) ?ScoreType {
            if (seq.getTwinWidth() > ctx.upper_bound) return null;
            var total_reds: u32 = 0;
            for (graph.node_list) |node| {
                total_reds += node.red_edges.cardinality();
            }
            return .{ @intCast(T, seq.getTwinWidth()), total_reds };
        }
    };
}

pub fn GraphScorerMinBlack(comptime T: type) type {
    return struct {
        upper_bound: T = 0,
        pub const ScoreType = struct { T, u32 };

        pub fn setUpperBound(own: *@This(), bound: T) void {
            own.upper_bound = bound;
        }

        pub fn compare(ctx: *@This(), first: ScoreType, second: ScoreType) std.math.Order {
            _ = ctx;
            if (first.@"0" == second.@"0") {
                return std.math.order(second.@"1", first.@"1");
            }
            // Gt if tww_first < tww_second
            return std.math.order(second.@"0", first.@"0");
        }

        pub fn evaluate(ctx: @This(), graph: *Graph(T), seq: *retraceable_contraction.RetraceableContractionSequence(T)) ?ScoreType {
            if (seq.getTwinWidth() > ctx.upper_bound) return null;
            var total_blacks: u32 = 0;
            for (graph.node_list) |node| {
                total_blacks+= node.black_edges.cardinality();
            }
            return .{ @intCast(T, seq.getTwinWidth()), total_blacks};
        }
    };
}

pub fn GraphScorerMinor(comptime T: type) type {
    return struct {
        upper_bound: T = 0,
        pub const ScoreType = struct { T, u32 };

        pub fn setUpperBound(own: *@This(), bound: T) void {
            own.upper_bound = bound;
        }

        pub fn compare(ctx: *@This(), first: ScoreType, second: ScoreType) std.math.Order {
            if (first.@"0" <= ctx.upper_bound and second.@"0" <= ctx.upper_bound) {
                return std.math.order(second.@"1", first.@"1");
            }
            // Gt if tww_first < tww_second
            return std.math.order(second.@"0", first.@"0");
        }

        pub fn evaluate(ctx: @This(), graph: *Graph(T), seq: *retraceable_contraction.RetraceableContractionSequence(T)) ?ScoreType {
            if (seq.getTwinWidth() > ctx.upper_bound) return null;
            var total_reds: u32 = 0;
            for (graph.node_list) |node| {
                total_reds += node.red_edges.cardinality();
            }
            return .{ @intCast(T, seq.getTwinWidth()), total_reds };
        }
    };
}

pub fn GraphScorerWeightedMajor(comptime T: type) type {
    return struct {
        upper_bound: T = 0,
        pub const ScoreType = struct { T, u32 };

        pub fn setUpperBound(own: *@This(), bound: T) void {
            own.upper_bound = bound;
        }

        pub fn compare(ctx: *@This(), first: ScoreType, second: ScoreType) std.math.Order {
            if (first.@"0" <= ctx.upper_bound and second.@"0" <= ctx.upper_bound) {
            	return std.math.order(first.@"1", second.@"1");
						}
            return std.math.order(second.@"0", first.@"0");
        }

        pub fn evaluate(ctx: @This(), graph: *Graph(T), seq: *retraceable_contraction.RetraceableContractionSequence(T)) ?ScoreType {
            if (seq.getTwinWidth() > ctx.upper_bound) return null;
            var total_reds: u32 = 0;
            for (graph.node_list) |node| {
                total_reds += (ctx.upper_bound - node.red_edges.cardinality());
            }
            return .{ @intCast(T, seq.getTwinWidth()), total_reds };
        }
    };
}

pub fn GraphScorerWeightedMinor(comptime T: type) type {
    return struct {
        upper_bound: T = 0,
        pub const ScoreType = struct { T, u32 };

        pub fn setUpperBound(own: *@This(), bound: T) void {
            own.upper_bound = bound;
        }

        pub fn compare(ctx: *@This(), first: ScoreType, second: ScoreType) std.math.Order {
            _ = ctx;
            if (first.@"0" == second.@"0") {
            	return std.math.order(first.@"1", second.@"1");
						}
            return std.math.order(second.@"0", first.@"0");
        }

        pub fn evaluate(ctx: @This(), graph: *Graph(T), seq: *retraceable_contraction.RetraceableContractionSequence(T)) ?ScoreType {
            if (seq.getTwinWidth() > ctx.upper_bound) return null;
            var total_reds: u32 = 0;
            for (graph.node_list) |node| {
                total_reds += (ctx.upper_bound - node.red_edges.cardinality());
            }
            return .{ @intCast(T, seq.getTwinWidth()), total_reds };
        }
    };
}

pub const MinSimType = enum {
	cumulative,
	min
};

pub fn GraphScorerMinSim(comptime T: type, comptime options: MinSimType) type {
    return struct {
        upper_bound: T = 0,
        visited: *bitset.FastBitSet,
        bfs: *bfs_mod.BfsQueue(T),
        subgraph: *subgraph.InducedSubGraph(T),
        node_tuples: []NodeTuple(T),
        round: u32 = 3,

        pub const ScoreType = struct { T, u32 };

				pub fn initAll(self: *@This()) void {
					@memset(self.node_tuples, NodeTuple(T) {.first = 0,.second = 0});
				}
        pub fn setUpperBound(own: *@This(), bound: T) void {
            own.upper_bound = bound;
        }

        pub fn compare(ctx: *@This(), first: ScoreType, second: ScoreType) std.math.Order {
            if (first.@"0" <= ctx.upper_bound and second.@"0" <= ctx.upper_bound) {
            	return std.math.order(second.@"1", first.@"1");
						}
            return std.math.order(second.@"0", first.@"0");
        }

        pub fn evaluate(ctx: *@This(), graph: *Graph(T), seq: *retraceable_contraction.RetraceableContractionSequence(T)) ?ScoreType {
            if (seq.getTwinWidth() > ctx.upper_bound) return null;
						ctx.round+=1;

            for (ctx.subgraph.nodes) |node| {
                if (graph.erased_nodes.get(node)) continue;
                ctx.visited.unsetAll();
                var bf = bfs_mod.bfs(T, node, graph, ctx.visited, ctx.bfs, .{ .max_level = 2, .kind = .both });
                while (bf.next()) |other_node| {
                    if (other_node >= node) continue;
                    if (graph.erased_nodes.get(other_node)) {
                        @panic("A BFS should never return an erased node!");
                    }

                    const cal_induced_tww = graph.calculateTwwOfMergeSurvivor(other_node, node);
										if(ctx.round != ctx.node_tuples[node].second) {
											ctx.node_tuples[node].second = ctx.round;
											if(options == .cumulative) {
												ctx.node_tuples[node].first = 0;
											}
											else if(options == .min) {
												ctx.node_tuples[node].first = std.math.maxInt(u32);
											}
										}
										if(ctx.round != ctx.node_tuples[other_node].second) {
											ctx.node_tuples[other_node].second = ctx.round;
											if(options == .cumulative) {
												ctx.node_tuples[other_node].first = 0;
											}
											else if(options == .min) {
												ctx.node_tuples[other_node].first = std.math.maxInt(u32);
											}
										}
										if(options == .cumulative) {
                    	ctx.node_tuples[node].first += cal_induced_tww;
                    	ctx.node_tuples[other_node].first += cal_induced_tww;
										}
										else if(options == .min) {
                    	ctx.node_tuples[node].first = std.math.min(ctx.node_tuples[node].first, cal_induced_tww);
                    	ctx.node_tuples[other_node].first = std.math.min(ctx.node_tuples[other_node].first,cal_induced_tww);
										}

                }
            }

            var total_sim_exceeds: u32 = 0;

            for (0..ctx.node_tuples.len) |i| {
								if(ctx.node_tuples[i].second == ctx.round) {
									if(options == .min) {
										if(ctx.node_tuples[i].first > ctx.upper_bound) {
            					total_sim_exceeds += ctx.node_tuples[i].first;
										}
									}
									else {
            				total_sim_exceeds += ctx.node_tuples[i].first;
									}
								}
            }
            return .{ @intCast(T, seq.getTwinWidth()), total_sim_exceeds };
        }
    };
}


pub fn GraphScorerRemainingMoves(comptime T: type, comptime NodeCtx: type) type {
    return struct {
        upper_bound: T = 0,
        visited: *bitset.FastBitSet,
        bfs: *bfs_mod.BfsQueue(T),
        subgraph: *subgraph.InducedSubGraph(T),
				node_ctx: *NodeCtx,

        pub const ScoreType = struct { T, u32 };

        pub fn setUpperBound(own: *@This(), bound: T) void {
            own.upper_bound = bound;
        }

        pub fn compare(ctx: *@This(), first: ScoreType, second: ScoreType) std.math.Order {
            if (first.@"0" <= ctx.upper_bound and second.@"0" <= ctx.upper_bound) {
            	return std.math.order(first.@"1", second.@"1");
						}
            return std.math.order(second.@"0", first.@"0");
        }

        pub fn evaluate(ctx: *@This(), graph: *Graph(T), seq: *retraceable_contraction.RetraceableContractionSequence(T)) ?ScoreType {
            if (seq.getTwinWidth() > ctx.upper_bound) return null;

						var score:u32 = 0;
            for (ctx.subgraph.nodes) |node| {
                if (graph.erased_nodes.get(node)) continue;
                ctx.visited.unsetAll();
                var bf = bfs_mod.bfs(T, node, graph, ctx.visited, ctx.bfs, .{ .max_level = 2, .kind = .both });
                while (bf.next()) |other_node| {
                    if (other_node >= node) continue;
                    if (graph.erased_nodes.get(other_node)) {
                        @panic("A BFS should never return an erased node!");
                    }

										if(ctx.node_ctx.map(node, other_node, graph) != null) {
											score+=1;
										}
                }
            }

            return .{ @intCast(T, seq.getTwinWidth()), score};
        }
    };
}
