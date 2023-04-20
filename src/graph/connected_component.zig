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

        pub fn solveGreedyTopK(self: *Self, comptime K: u32, comptime P: u32, solver: *solver_resources.SolverResources(T,K,P)) !T {
					const result = try self.subgraph.solveGreedyTopK(K,P,&self.current_contraction_seq,solver,true);
					self.tww = result;
					return result;
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
