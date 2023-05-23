const mg = @import("matrix_graph.zig");
const Graph = @import("../graph/graph.zig").Graph;
const sg = @import("../graph/subgraph.zig");
const std = @import("std");

const MatrixDispatchError = error{InputGraphTooLarge};

pub const MatrixGraphUnion = union(enum) {
    g8: mg.MatrixGraph(8),
    g16: mg.MatrixGraph(16),
    g32: mg.MatrixGraph(32),
    g64: mg.MatrixGraph(64),
    g128: mg.MatrixGraph(128),
    g256: mg.MatrixGraph(256),
    g512: mg.MatrixGraph(512),
    g1024: mg.MatrixGraph(1024),
    g2048: mg.MatrixGraph(2048),
    g4096: mg.MatrixGraph(4096),
    g8192: mg.MatrixGraph(8192),
    empty: void,
};

pub const MatrixGraphFromInducedSubGraph = struct {
    const Self = @This();
    allocator: std.mem.Allocator,
    graph: MatrixGraphUnion,
    mapping: std.ArrayList(u32),

    pub fn init(allocator: std.mem.Allocator, nodes: u32) !Self {
        var mapping = try std.ArrayList(u32).initCapacity(allocator, nodes);
        errdefer mapping.deinit();

        return Self{ .allocator = allocator, .graph = undefined, .mapping = mapping };
    }

    pub fn deinit(self: *Self) void {
        self.dispatch(.{}, deinit_graph);
        self.mapping.deinit();
    }

    fn deinit_graph(_: anytype, graph: anytype, _: anytype) void {
        graph.deinit();
    }

    pub fn dispatch(self: *Self, context: anytype, callback: anytype) void {
        switch (self.graph) {
            .g8 => return callback(context, self.graph.g8, self.mapping),
            .g16 => return callback(context, self.graph.g16, self.mapping),
            .g32 => return callback(context, self.graph.g32, self.mapping),
            .g64 => return callback(context, self.graph.g64, self.mapping),
            .g128 => return callback(context, self.graph.g128, self.mapping),
            .g256 => return callback(context, self.graph.g256, self.mapping),
            .g512 => return callback(context, self.graph.g512, self.mapping),
            .g1024 => return callback(context, self.graph.g1024, self.mapping),
            .g2048 => return callback(context, self.graph.g2048, self.mapping),
            .g4096 => return callback(context, self.graph.g4096, self.mapping),
            .g8192 => return callback(context, self.graph.g8192, self.mapping),
            else => {},
        }
        unreachable;
    }

    pub fn fromInducedSubgraph(self: *Self, comptime G: type, subgraph: anytype) !void {
        var graph = try G.new(self.allocator);

        for (subgraph.nodes) |org_node| {
            try self.mapping.append(@intCast(u32, org_node));
        }

        std.sort.sort(u32, self.mapping.items, {}, std.sort.asc(u32));

        // TODO: maybe faster if we use orderedIterator and avoid binary search
        for (subgraph.nodes) |org_node| {
            var org_idx = self.newIdOf(org_node);

            var neigh_iter = subgraph.graph.node_list[org_node].unorderedIterator();
            while (neigh_iter.next()) |org_neigh| {
                var nei_idx = self.newIdOf(org_neigh);
                var color = if (neigh_iter.red) mg.Color.Red else mg.Color.Black;
                _ = graph.addEdge(org_idx, nei_idx, color);
            }
        }

        inline for (std.meta.fields(MatrixGraphUnion)) |f| {
            if (G == f.type) {
                self.graph = @unionInit(MatrixGraphUnion, f.name, graph);

                return;
            }
        }

        return MatrixDispatchError.InputGraphTooLarge;
    }

    fn newIdOf(self: *const Self, org: u32) u32 {
        var idx_into_mapping = std.sort.binarySearch(u32, @intCast(u32, org), self.mapping.items, {}, order_u32).?;

        return @intCast(u32, idx_into_mapping);
    }
};

pub fn matrixGraphUnionFromInducedSubGraph(comptime T: type, subgraph: *const sg.InducedSubGraph(T), allocator: std.mem.Allocator) !MatrixGraphFromInducedSubGraph {
    var context = .{ .subgraph = subgraph, .result = try MatrixGraphFromInducedSubGraph.init(allocator, @intCast(u32, subgraph.nodes.len)) };

    try matrixGraphDispatch(@intCast(u32, subgraph.nodes.len), &context, matrixGraphUnionFromInducedSubGraphWorker);

    return context.result;
}

fn matrixGraphUnionFromInducedSubGraphWorker(context: anytype, comptime G: type) void {
    context.result.fromInducedSubgraph(G, context.subgraph) catch @panic("dispatch failed");
}

fn order_u32(context: void, lhs: u32, rhs: u32) std.math.Order {
    _ = context;
    return std.math.order(lhs, rhs);
}

pub fn matrixGraphDispatch(num_nodes: u32, context: anytype, callback: anytype) !void {
    inline for (std.meta.fields(MatrixGraphUnion)) |f| {
        if (f.type != void and num_nodes <= f.type.NumNodes) {
            return callback(context, mg.MatrixGraph(f.type.NumNodes));
        }
    }
    return MatrixDispatchError.InputGraphTooLarge;
}

fn callback_test(required_size: u32, comptime G: type) void {
    var graph = G.new();
    std.debug.assert(required_size <= graph.number_of_nodes());
    std.debug.assert(required_size < 4 or graph.number_of_nodes() <= 2 * required_size);
}

test "test" {
    var n: u32 = 1;
    while (n < 50) : (n += 1) {
        matrixGraphDispatch(n, n, callback_test);
    }
}
