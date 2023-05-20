const std = @import("std");
const bitset = @import("../util/two_level_bitset.zig");
const Graph = @import("graph.zig").Graph;
const comptime_util = @import("../util/comptime_checks.zig");
const solver = @import("solver.zig");

pub const DfsKind = enum { red, black, both };

pub const DfsOptions = struct {
    max_level: u32 = std.math.maxInt(u32),
    kind: DfsKind = DfsKind.both,
};

pub fn DfsIterator(comptime T: type, comptime options: DfsOptions) type {
    comptime if (options.max_level < 1) {
        @compileError("BFS max level must be at least 1");
    };
    comptime if (!comptime_util.checkIfIsCompatibleInteger(T)) {
        @compileError("Type must either be u8, u16 or u32");
    };

    return struct {
        const Self = @This();
        visited: *bitset.FastBitSet,
        stack: std.ArrayListUnmanaged(T),
        depth_stack: std.ArrayListUnmanaged(T),
        parent_stack: std.ArrayListUnmanaged(T),
        graph: *const Graph(T),
        level: T,
        parent: T,

        pub inline fn next(self: *Self) ?T {
            while (self.stack.popOrNull()) |current| {
                self.level = self.depth_stack.pop();
                self.parent = self.parent_stack.pop();

                if (!self.visited.setExists(current)) {
                    if (options.kind == .black or options.kind == .both) {
                        var black_iter = self.graph.node_list[current].black_edges.iterator();
                        while (black_iter.next()) |item| {
                            self.stack.append(self.graph.allocator, item);
                            self.depth_stack.append(self.graph.allocator, self.level + 1);
                            self.parent_stack.append(current);
                        }
                    }
                    if (options.kind == .red or options.kind == .both) {
                        var red_iter = self.graph.node_list[current].red_edges.iterator();
                        while (red_iter.next()) |item| {
                            self.stack.append(self.graph.allocator, item);
                            self.depth_stack.append(self.graph.allocator, self.level + 1);
                            self.parent_stack.append(current);
                        }
                    }
                    return current;
                }
            }
            return null;
        }

        pub inline fn deinit(self: *Self) void {
            self.stack.deinit(self.graph.allocator);
            self.depth_stack.deinit(self.graph.allocator);
            self.parent_stack.deinit(self.parent_stack);
        }
    };
}

pub inline fn dfs(comptime T: type, start_node: T, graph: *const Graph(T), visited: *bitset.FastBitSet, comptime options: DfsOptions) !DfsIterator(T, options) {
    var stack = try std.ArrayListUnmanaged(T).initCapacity(graph.allocator, graph.number_of_edges);
    var depth_stack = try std.ArrayListUnmanaged(T).initCapacity(graph.allocator, graph.number_of_edges);
    var parent_stack = try std.ArrayListUnmanaged(T).initCapacity(graph.allocator, graph.number_of_edges);

    try stack.append(graph.allocator, start_node);
    try depth_stack.append(graph.allocator, 0);
    try parent_stack.append(graph.allocator, start_node);

    return DfsIterator(T, options){
        .graph = graph,
        .stack = stack,
        .visited = visited,
        .depth_stack = depth_stack,
        .parent_stack = parent_stack,
        .level = 0,
    };
}

pub fn NodeStackEntry(comptime T: type) type {
    return struct { node: T, parent: T, recurse_update: bool };
}

pub fn findArticulationPoints(comptime T: type, root_node: T, graph: *const Graph(T), visited: *bitset.FastBitSet, resources: *solver.ArticulationPointResources(T)) !void {
    var failing_allocator = std.heap.FixedBufferAllocator.init(&[_]u8{});
    var allc = failing_allocator.allocator();

    resources.articulation_points.unsetAll();
    resources.parent_array[root_node] = root_node;
    resources.depth_array[root_node] = 0;
    resources.low_array[root_node] = 0;
    resources.node_stack.shrinkRetainingCapacity(0);

    try resources.node_stack.append(allc, NodeStackEntry(T){ .node = root_node, .parent = root_node, .recurse_update = false });

    var root_children: T = 0;

    while (resources.node_stack.popOrNull()) |item| {

        // Item does not exist yet
        if (!visited.setExists(item.node)) {
            if (item.parent == root_node and item.node != root_node) {
                root_children += 1;
            }
            resources.parent_array[item.node] = item.parent;
            resources.depth_array[item.node] = resources.depth_array[item.parent] + 1;
            resources.low_array[item.node] = resources.depth_array[item.parent] + 1;
            resources.node_count[item.node] = 1;

            var iter_nodes = graph.node_list[item.node].orderedIterator();
            while (iter_nodes.next()) |next_node| {
                if (!visited.get(next_node)) {
                    try resources.node_stack.append(allc, NodeStackEntry(T){ .node = item.node, .parent = next_node, .recurse_update = true });
                    try resources.node_stack.append(allc, NodeStackEntry(T){ .node = next_node, .parent = item.node, .recurse_update = false });
                } else if (next_node != resources.parent_array[item.node]) {
                    resources.low_array[item.node] = std.math.min(resources.low_array[item.node], resources.depth_array[next_node]);
                }
            }
        } else {
            if (item.recurse_update and resources.parent_array[item.parent] == item.node) {
                // parent are child are swapped and this is after the depth first recursion
                resources.low_array[item.node] = std.math.min(resources.low_array[item.node], resources.low_array[item.parent]);

                // Calculate the nodes in the dfs tree below us
                if (resources.parent_array[item.parent] == item.node) {
                    resources.node_count[item.node] += resources.node_count[item.parent];
                }
                // Propagate upwards
                if (resources.parent_array[item.node] != item.node and resources.low_array[item.parent] >= resources.depth_array[item.node]) {
                    resources.articulation_points.set(item.node);
                }
            } else if (item.node != resources.parent_array[item.parent]) {
                resources.low_array[item.parent] = std.math.min(resources.low_array[item.parent], resources.depth_array[item.node]);
            }
        }
    }
    if (root_children > 1) {
        resources.articulation_points.set(root_node);
    }
}

test "Check articulation points I" {
    var allocator = std.heap.GeneralPurposeAllocator(.{}){};
    var graph = try Graph(u8).new(5, allocator.allocator());
    try graph.addEdge(0, 1);
    try graph.addEdge(0, 2);
    try graph.addEdge(2, 1);
    try graph.addEdge(0, 3);
    try graph.addEdge(3, 4);

    var res = try solver.ArticulationPointResources(u8).init(&graph);

    try findArticulationPoints(u8, 0, &graph, &graph.scratch_bitset, &res);

    var iter = res.articulation_points.iter();
    try std.testing.expectEqual(iter.next().?, 0);
    try std.testing.expectEqual(iter.next().?, 3);
    try std.testing.expectEqual(iter.next(), null);
}

test "Check articulation points II" {
    var allocator = std.heap.GeneralPurposeAllocator(.{}){};
    var graph = try Graph(u8).new(4, allocator.allocator());
    try graph.addEdge(0, 1);
    try graph.addEdge(1, 2);
    try graph.addEdge(2, 3);

    var res = try solver.ArticulationPointResources(u8).init(&graph);
    try findArticulationPoints(u8, 0, &graph, &graph.scratch_bitset, &res);

    var iter = res.articulation_points.iter();
    try std.testing.expectEqual(iter.next().?, 1);
    try std.testing.expectEqual(iter.next().?, 2);
    try std.testing.expectEqual(iter.next(), null);
}

test "Check articulation points III" {
    var allocator = std.heap.GeneralPurposeAllocator(.{}){};
    var graph = try Graph(u8).new(7, allocator.allocator());
    try graph.addEdge(0, 1);
    try graph.addEdge(1, 2);
    try graph.addEdge(2, 0);

    try graph.addEdge(1, 3);
    try graph.addEdge(1, 4);
    try graph.addEdge(1, 6);

    try graph.addEdge(3, 5);
    try graph.addEdge(4, 5);
    var res = try solver.ArticulationPointResources(u8).init(&graph);

    try findArticulationPoints(u8, 0, &graph, &graph.scratch_bitset, &res);

    var iter = res.articulation_points.iter();
    try std.testing.expectEqual(iter.next().?, 1);
    try std.testing.expectEqual(iter.next(), null);
}
