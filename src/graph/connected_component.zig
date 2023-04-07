const std = @import("std");
const edge_list = @import("edge_list.zig");
const contraction = @import("../tww/contraction_sequence.zig");
const graph_mod = @import("graph.zig");
const comptime_util = @import("../util/comptime_checks.zig");
const Graph = graph_mod.Graph;
const CompactField = graph_mod.CompactField;

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
            return std.math.order(rhs.tww,lhs.tww);
        }
    };
}

pub fn ConnectedComponent(comptime T: type) type {
    comptime if (!comptime_util.checkIfIsCompatibleInteger(T)) {
        @compileError("The type of ConnectedComponent must either be u8,u16 or u32!");
    };

    return struct {
        const Self = @This();
        nodes: edge_list.ParametrizedUnsortedArrayList(T),
        best_contraction_sequence: contraction.ContractionSequence(T),
        tww: T,
        bfs_levels: u32,

				pub fn deinit(self: *Self, allocator: std.mem.Allocator) void {
					self.nodes.deinit(allocator);
					self.best_contraction_sequence.deinit(allocator);
				}

        pub fn init(allocator: std.mem.Allocator, nodes: edge_list.ParametrizedUnsortedArrayList(T), bfs_level: u32) !Self {
            var contraction_seq = try contraction.ContractionSequence(T).init(allocator, nodes.cardinality());
            var iter = nodes.iterator();
            var previous: ?T = null;
            while (iter.next()) |item| {
                if (previous) |into| {
                    try contraction_seq.addContraction(contraction.Contraction(T){ .erased = item, .survivor = into });
                } else {
                    previous = item;
                }
            }
            return .{ .bfs_levels = bfs_level, .nodes = nodes, .tww = @intCast(T,nodes.cardinality() - 1), .best_contraction_sequence = contraction_seq };
        }
    };
}
