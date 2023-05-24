const std = @import("std");
const Graph = @import("../graph/graph.zig").Graph;
const Subgraph = @import("../graph/subgraph.zig").InducedSubGraph;
const retraceable_contraction = @import("../tww/retraceable_contraction_sequence.zig");
const bfs_mod = @import("../graph/bfs.zig");
const set = @import("../util/two_level_bitset.zig");

pub fn NodeScorerPotential(comptime T: type) type {
    return struct {
        const NodePriorityType = struct {
            survivor: T,
            erased: T,
            potential: Graph(T).InducedTwinWidthPotential,
        };

        pub const MapType = NodePriorityType;
        upper_bound: T = 0,
        sequence: *retraceable_contraction.RetraceableContractionSequence(T),

        pub fn setUpperBound(own: *@This(), bound: T) void {
            own.upper_bound = bound;
        }

        pub fn map(ctx: *@This(), first: T, second: T, graph: *Graph(T)) ?MapType {
            var def = Graph(T).InducedTwinWidthPotential.default();
            const result = graph.calculateInducedTwwPotential(first, second, &def, ctx.upper_bound);
            if (result.tww > ctx.upper_bound) return null;
            return NodePriorityType{ .survivor = second, .erased = first, .potential = result };
        }

        pub fn compare(ctx: *@This(), first: MapType, second: MapType) std.math.Order {
            return first.potential.order(second.potential, ctx.sequence.getTwinWidth());
        }
    };
}

pub fn NodeScorerEdgeReducer(comptime T: type) type {
    return struct {
        const NodePriorityType = struct { survivor: T, erased: T, reduced_edges: i32, tww: T };

        pub const MapType = NodePriorityType;
        upper_bound: T = 0,
        seq: *retraceable_contraction.RetraceableContractionSequence(T),

        pub fn setUpperBound(own: *@This(), bound: T) void {
            own.upper_bound = bound;
        }

        pub fn map(ctx: *@This(), first: T, second: T, graph: *Graph(T)) ?MapType {
            const result = graph.calculateTwwOfMergeSurvivor(first, second);
            if (result > ctx.upper_bound) return null;

            const total = graph.node_list[first].red_edges.cardinality() + graph.node_list[second].red_edges.cardinality();

            return NodePriorityType{ .survivor = second, .erased = first, .tww = result, .reduced_edges = @intCast(i32, total) - @intCast(i32, result) };
        }

        pub fn compare(ctx: *@This(), first: MapType, second: MapType) std.math.Order {
            if (first.tww > ctx.seq.getTwinWidth() or second.tww > ctx.seq.getTwinWidth()) {
                return std.math.order(second.tww, first.tww);
            }
            return std.math.order(first.reduced_edges, second.reduced_edges);
        }
    };
}

pub fn NodeScorerSimple(comptime T: type) type {
    return struct {
        const NodePriorityType = struct {
            survivor: T,
            erased: T,
            merge_tww: T,
        };

        pub const MapType = NodePriorityType;
        upper_bound: T = 0,

        pub fn setUpperBound(own: *@This(), bound: T) void {
            own.upper_bound = bound;
        }

        pub fn map(ctx: *@This(), first: T, second: T, graph: *Graph(T)) ?MapType {
            const result = graph.calculateTwwOfMergeSurvivor(first, second);
            if (result > ctx.upper_bound) return null;
            return NodePriorityType{ .survivor = second, .erased = first, .merge_tww = result };
        }

        pub fn compare(ctx: *@This(), first: MapType, second: MapType) std.math.Order {
            _ = ctx;
            return std.math.order(second.merge_tww, first.merge_tww);
        }
    };
}

pub fn NodeScorerInduced(comptime T: type) type {
    return struct {
        const NodePriorityType = struct {
            survivor: T,
            erased: T,
            merge_tww: T,
        };

        pub const MapType = NodePriorityType;
        upper_bound: T = 0,

        pub fn setUpperBound(own: *@This(), bound: T) void {
            own.upper_bound = bound;
        }

        pub fn map(ctx: *@This(), first: T, second: T, graph: *Graph(T)) ?MapType {
            const result = graph.calculateInducedTww(first, second, null);
            if (result.tww > ctx.upper_bound) return null;
            return NodePriorityType{ .survivor = second, .erased = first, .merge_tww = result.tww };
        }

        pub fn compare(ctx: *@This(), first: MapType, second: MapType) std.math.Order {
            _ = ctx;
            return std.math.order(second.merge_tww, first.merge_tww);
        }
    };
}

pub fn NodeScorerMaximizeMoves(comptime T: type) type {
    return struct {
        const NodePriorityType = struct {
            survivor: T,
            erased: T,
            moves: u64,
        };

        pub const MapType = NodePriorityType;
        upper_bound: T = 0,
        seq: *retraceable_contraction.RetraceableContractionSequence(T),
        subgraph: *Subgraph(T),
        visited: *set.FastBitSet,
        bfs: *bfs_mod.BfsQueue(T),

        pub fn setUpperBound(own: *@This(), bound: T) void {
            own.upper_bound = bound;
        }

        pub fn map(ctx: *@This(), first: T, second: T, graph: *Graph(T)) ?MapType {
            const result = graph.calculateInducedTww(first, second, null);
            if (result.tww > ctx.upper_bound) return null;

            _ = graph.addContraction(first, second, ctx.seq) catch unreachable;
            var moves: u64 = 0;

            for (ctx.subgraph.nodes) |node| {
                if (graph.erased_nodes.get(node)) continue;
                ctx.visited.unsetAll();
                var bf = bfs_mod.bfs(T, node, graph, ctx.visited, ctx.bfs, .{ .max_level = 2, .kind = .both });
                while (bf.next()) |other_node| {
                    if (other_node >= node) continue;
                    if (graph.erased_nodes.get(other_node)) {
                        @panic("A BFS should never return an erased node!");
                    }

                    const cal_induced_tww = graph.calculateInducedTww(other_node, node, null);
                    if (cal_induced_tww.tww > ctx.upper_bound) continue;
                    moves += 1;
                }
            }
            graph.revertLastContraction(ctx.seq) catch @panic("Should never happen");
            if (moves == 0) return null;

            return NodePriorityType{ .erased = first, .survivor = second, .moves = moves };
        }

        pub fn compare(ctx: *@This(), first: MapType, second: MapType) std.math.Order {
            _ = ctx;
            return std.math.order(first.moves, second.moves);
        }
    };
}

pub fn NodeScorerNewReds(comptime T: type) type {
    return struct {
        const NodePriorityType = struct {
            survivor: T,
            erased: T,
            merge_tww: T,
        };

        pub const MapType = NodePriorityType;
        upper_bound: T = 0,

        pub fn setUpperBound(own: *@This(), bound: T) void {
            own.upper_bound = bound;
        }

        pub fn map(ctx: *@This(), first: T, second: T, graph: *Graph(T)) ?MapType {
            const result = graph.calculateMaxTwwOfNewNeighbors(first, second).@"0";
            if (result > ctx.upper_bound) return null;
            return NodePriorityType{ .survivor = second, .erased = first, .merge_tww = result };
        }

        pub fn compare(ctx: *@This(), first: MapType, second: MapType) std.math.Order {
            _ = ctx;
            return std.math.order(second.merge_tww, first.merge_tww);
        }
    };
}

pub fn NodeScorerTotalRedDeg(comptime T: type) type {
    return struct {
        const NodePriorityType = struct {
            survivor: T,
            erased: T,
            total_red_deg: u32,
            tww: u32,
        };

        pub const MapType = NodePriorityType;
        upper_bound: T = 0,
        seq: *retraceable_contraction.RetraceableContractionSequence(T),
        subgraph: *Subgraph(T),

        pub fn setUpperBound(own: *@This(), bound: T) void {
            own.upper_bound = bound;
        }

        pub fn map(ctx: *@This(), first: T, second: T, graph: *Graph(T)) ?MapType {
            const result = graph.calculateInducedTww(first, second, null).tww;
            if (result > ctx.upper_bound) return null;
            _ = graph.addContraction(first, second, ctx.seq) catch unreachable;
            var total_reds: u32 = 0;
            for (ctx.subgraph.nodes) |node| {
                if (graph.erased_nodes.get(node)) continue;
                total_reds += graph.node_list[node].red_edges.cardinality();
            }
            graph.revertLastContraction(ctx.seq) catch unreachable;

            return NodePriorityType{ .survivor = second, .erased = first, .total_red_deg = total_reds, .tww = result };
        }

        pub fn compare(ctx: *@This(), first: MapType, second: MapType) std.math.Order {
            _ = ctx;
            if (second.tww == first.tww) {
                return std.math.order(second.total_red_deg, first.total_red_deg);
            }
            return std.math.order(second.tww, first.tww);
        }
    };
}
