const std = @import("std");
const bitset = @import("../util/two_level_bitset.zig");
const graph_mode = @import("graph.zig");
const comptime_util = @import("../util/comptime_checks.zig");
const topk = @import("../util/top_k_scorer.zig");
const Graph = graph_mode.Graph;
const CompactField = graph_mode.CompactField;

pub const BfsKind = enum { red, black, both };

pub const BfsOptions = struct {
    max_level: u32 = std.math.maxInt(u32),
    kind: BfsKind = BfsKind.both,
};

pub fn BfsIterator(comptime T: type, comptime options: BfsOptions) type {
    comptime if (options.max_level < 1) {
        @compileError("BFS max level must be at least 1");
    };
    comptime if (!comptime_util.checkIfIsCompatibleInteger(T)) {
        @compileError("Type must either be u8, u16 or u32");
    };
    return struct {
        visited: *bitset.FastBitSet,
        graph: *const Graph(T),
        stack: *BfsQueue(T),
        iter_current: BfsQueue(T).BfsQueueIterator,
        level: T,
        pub inline fn visited(self: @This()) *bitset.FastBitSet {
            return &self.visited;
        }
        pub inline fn next(self: *@This()) ?T {
            while (self.level < options.max_level) {
                while (self.iter_current.next()) |id| {
                    if (options.kind == .black or options.kind == .both) {
                        var black_iter = self.graph.node_list[id].black_edges.iterator();
                        while (black_iter.next()) |item| {
                            if (!self.visited.setExists(item)) {
                                self.stack.addNext(@intCast(T, item));
                            }
                        }
                    }
                    if (options.kind == .red or options.kind == .both) {
                        var red_iter = self.graph.node_list[id].red_edges.iterator();
                        while (red_iter.next()) |item| {
                            if (!self.visited.setExists(item)) {
                                self.stack.addNext(@intCast(T, item));
                            }
                        }
                    }
                    return id;
                }

                self.level += 1;
                if (self.stack.nextFrontierSize() == 0) {
                    return null;
                }
                self.stack.swapFrontiers();
                self.iter_current = self.stack.iterator();
            }
            while (self.iter_current.next()) |id| {
                self.visited.set(id);
                return id;
            }
            return null;
        }
    };
}

pub fn BfsQueue(comptime T: type) type {
    comptime if (T != u8 and T != u16 and T != u32) {
        @compileError("Type must either be u8, u16 or u32");
    };

    return struct {
        const Self = @This();
        current: []T,
        current_write_ptr: u32,
        next: []T,
        next_write_ptr: u32,
        pub inline fn init(allocator: std.mem.Allocator, graph_size: u32) !BfsQueue(T) {
            const current = try allocator.alloc(T, graph_size);
            const next = try allocator.alloc(T, graph_size);
            return .{
                .current = current,
                .next = next,
                .current_write_ptr = 0,
                .next_write_ptr = 0,
            };
        }

        pub inline fn addNext(self: *Self, next: T) void {
            self.next[self.next_write_ptr] = next;
            self.next_write_ptr += 1;
        }

        pub inline fn nextFrontierSize(self: *const Self) u32 {
            return self.next_write_ptr;
        }

        pub inline fn swapFrontiers(self: *Self) void {
            const next_tmp = self.next;
            self.next = self.current;

            self.current = next_tmp;
            self.current_write_ptr = self.next_write_ptr;

            self.next_write_ptr = 0;
        }

        pub inline fn clear(self: *Self) void {
            self.current_write_ptr = 0;
            self.next_write_ptr = 0;
        }

        pub const BfsQueueIterator = struct {
            stack: *const BfsQueue(T),
            index: u32,
            pub inline fn next(self: *BfsQueueIterator) ?T {
                if (self.index >= self.stack.current_write_ptr) {
                    return null;
                }
                const item = self.stack.current[self.index];
                self.index += 1;
                return item;
            }
        };

        pub inline fn iterator(self: *const Self) BfsQueueIterator {
            return BfsQueueIterator{
                .stack = self,
                .index = 0,
            };
        }

        pub inline fn deinit(self: *const Self, allocator: std.mem.Allocator) void {
            allocator.free(self.next);
            allocator.free(self.current);
        }
    };
}

pub inline fn bfs(comptime T: type, start_node: T, graph: *const Graph(T), visited: *bitset.FastBitSet, stack: *BfsQueue(T), comptime options: BfsOptions) BfsIterator(T, options) {
    stack.clear();
    stack.addNext(start_node);
    stack.swapFrontiers();

    visited.set(start_node);
    return BfsIterator(T, options){
        .graph = graph,
        .stack = stack,
        .visited = visited,
        .iter_current = stack.iterator(),
        .level = 0,
    };
}

pub inline fn bfs_topk_high_performance_inner(comptime T: type, comptime K: u32, start_node: T, graph: *const Graph(T), visited: *bitset.FastBitSet, scorer: *topk.TopKScorer(T, K), comptime large_root: bool) void {
    const own_card = graph.node_list[start_node].cardinality();
    {
        var black_iter = graph.node_list[start_node].black_edges.iterator();
        while (black_iter.next()) |item| {
            visited.set(item);
            scorer.addVisit(item, false);
        }
        var red_iter = graph.node_list[start_node].red_edges.iterator();
        while (red_iter.next()) |item| {
            visited.set(item);
            scorer.addVisit(item, false);
        }
    }

    const copy_frontier = scorer.write_ptr;

    // Prune large degree direct neighbors
    if (own_card > 100) {
        var generator = std.rand.DefaultPrng.init(start_node +% own_card);
        // Prune at most ~17% of direct neighbors for high degree nodes
        for (scorer.node_list[0..copy_frontier]) |id| {
            if (generator.next() > std.math.maxInt(u64) / 3) {
                if (!large_root and graph.node_list[id].isLargeNode()) {
                    var iter_large = graph.node_list[id].takeIterator(40);
                    while (iter_large.next()) |item| {
                        scorer.addVisit(item, visited.setExists(item));
                    }
                } else {
                    var black_iter = graph.node_list[id].black_edges.iterator();
                    while (black_iter.next()) |item| {
                        scorer.addVisit(item, visited.setExists(item));
                    }
                    var red_iter = graph.node_list[id].red_edges.iterator();
                    while (red_iter.next()) |item| {
                        scorer.addVisit(item, visited.setExists(item));
                    }
                }
            }
        }
    } else {
        for (scorer.node_list[0..copy_frontier]) |id| {
            if (!large_root and graph.node_list[id].isLargeNode()) {
                var iter_large = graph.node_list[id].takeIterator(40);
                while (iter_large.next()) |item| {
                    scorer.addVisit(item, visited.setExists(item));
                }
            } else {
                var black_iter = graph.node_list[id].black_edges.iterator();
                while (black_iter.next()) |item| {
                    scorer.addVisit(item, visited.setExists(item));
                }
                var red_iter = graph.node_list[id].red_edges.iterator();
                while (red_iter.next()) |item| {
                    scorer.addVisit(item, visited.setExists(item));
                }
            }
        }
    }
}

pub inline fn bfs_topk_high_performance(comptime T: type, comptime K: u32, start_node: T, graph: *const Graph(T), visited: *bitset.FastBitSet, scorer: *topk.TopKScorer(T, K)) void {
    visited.set(start_node);
    scorer.reset();

    if (graph.node_list[start_node].cardinality() > 10000 and !graph.node_list[start_node].isLargeNode()) {
        graph.node_list[start_node].promoteToLargeDegreeNode(graph) catch unreachable;
    }
    var root_large = graph.node_list[start_node].isLargeNode();
    if (root_large) {
        bfs_topk_high_performance_inner(T, K, start_node, graph, visited, scorer, true);
    } else {
        bfs_topk_high_performance_inner(T, K, start_node, graph, visited, scorer, false);
    }
    _ = visited.unset(start_node);
    scorer.backing_storage[start_node] = -1_000_000;
}

pub inline fn bfs_topk_high_performance_incremental(comptime T: type, comptime K: u32, start_node: T, erased_node: T, level_one_merge: bool, new_reds_erased: []T, graph: *const Graph(T), visited: *bitset.FastBitSet, scorer: *topk.TopKScorer(T, K)) bool {
    if (level_one_merge) {
        // Decrease score for all neighbors of erased
        var nb_iter = graph.node_list[erased_node].unorderedIterator();
        while (nb_iter.next()) |item| {
            scorer.backing_storage[item] -= 2;
        }
    } else {
        var nb_iter = graph.node_list[erased_node].unorderedIterator();
        var visited_already: u32 = 0;
        while (nb_iter.next()) |item| {
            if (visited.get(item)) {
                visited_already += 1;
            }
        }

        if (visited_already == 0) return false;
    }

    scorer.priority_queue.len = 0;

    visited.set(start_node);
    var large_root = graph.node_list[start_node].isLargeNode();

    for (new_reds_erased) |red| {
        // Visit ourselfes
        scorer.addVisit(red, visited.setExists(red));

        // Visit neighbors
        if (!large_root and graph.node_list[red].isLargeNode()) {
            var iter_large = graph.node_list[red].takeIterator(40);
            while (iter_large.next()) |item| {
                scorer.addVisit(item, visited.setExists(item));
            }
        } else {
            var black_iter = graph.node_list[red].black_edges.iterator();
            while (black_iter.next()) |item| {
                scorer.addVisit(item, visited.setExists(item));
            }
            var red_iter = graph.node_list[red].red_edges.iterator();
            while (red_iter.next()) |item| {
                scorer.addVisit(item, visited.setExists(item));
            }
        }
    }

    _ = visited.unset(start_node);
    scorer.backing_storage[start_node] = -1_000_000;
    return true;
}

pub inline fn bfs_topk(comptime T: type, comptime K: u32, start_node: T, graph: *const Graph(T), visited: *bitset.FastBitSet, stack: *BfsQueue(T), scorer: *topk.TopKScorer(T, K), comptime options: BfsOptions) void {
    stack.clear();
    stack.addNext(start_node);
    stack.swapFrontiers();

    visited.set(start_node);
    scorer.reset();

    var iter_current = stack.iterator();
    if (graph.node_list[start_node].cardinality() > 10000 and !graph.node_list[start_node].isLargeNode()) {
        graph.node_list[start_node].promoteToLargeDegreeNode(graph) catch unreachable;
    }
    var root_large = graph.node_list[start_node].isLargeNode();

    for (0..(options.max_level - 1)) |_| {
        while (iter_current.next()) |id| {
            if (graph.node_list[id].cardinality() > 10000 and !graph.node_list[id].isLargeNode()) {
                graph.node_list[id].promoteToLargeDegreeNode(graph) catch unreachable;
            }
            if (!root_large and graph.node_list[id].isLargeNode()) {
                // Was 30 before
                var iter_large = graph.node_list[id].takeIterator(40);
                while (iter_large.next()) |item| {
                    if (!visited.setExists(item)) {
                        stack.addNext(@intCast(T, item));
                        scorer.addNewNode(item);
                    }
                    scorer.addVisit(item);
                }
            } else {
                if (options.kind == .black or options.kind == .both) {
                    var black_iter = graph.node_list[id].black_edges.iterator();
                    while (black_iter.next()) |item| {
                        if (!visited.setExists(item)) {
                            stack.addNext(@intCast(T, item));
                            scorer.addNewNode(item);
                        }
                        scorer.addVisit(item);
                    }
                }
                if (options.kind == .red or options.kind == .both) {
                    var red_iter = graph.node_list[id].red_edges.iterator();
                    while (red_iter.next()) |item| {
                        if (!visited.setExists(item)) {
                            stack.addNext(@intCast(T, item));
                            scorer.addNewNode(item);
                        }
                        scorer.addVisit(item);
                    }
                }
            }
        }

        if (stack.nextFrontierSize() == 0) {
            return;
        }
        stack.swapFrontiers();
        iter_current = stack.iterator();
    }
    while (iter_current.next()) |id| {
        if (!root_large and graph.node_list[id].isLargeNode()) {
            var iter_large = graph.node_list[id].takeIterator(40);
            while (iter_large.next()) |item| {
                if (!visited.setExists(item)) {
                    scorer.addNewNode(item);
                }
                scorer.addVisit(item);
            }
        } else {
            if (options.kind == .black or options.kind == .both) {
                var black_iter = graph.node_list[id].black_edges.iterator();
                while (black_iter.next()) |item| {
                    if (!visited.setExists(item)) {
                        scorer.addNewNode(item);
                    }
                    scorer.addVisit(item);
                }
            }
            if (options.kind == .red or options.kind == .both) {
                var red_iter = graph.node_list[id].red_edges.iterator();
                while (red_iter.next()) |item| {
                    if (!visited.setExists(item)) {
                        scorer.addNewNode(item);
                    }
                    scorer.addVisit(item);
                }
            }
        }
    }
}

pub inline fn bfs_topk_incremental(comptime T: type, comptime K: u32, start_node: T, graph: *const Graph(T), scorer: *topk.TopKScorer(T, K), new_level_one_nodes: *std.ArrayListUnmanaged(T), level_one_merge: bool, comptime options: BfsOptions) void {
    _ = start_node;
    for (new_level_one_nodes.items) |level_one_node| {
        if (!level_one_merge) {
            scorer.addVisit(level_one_node, graph.node_list[level_one_node].cardinality());
        }
        if (options.kind == .black or options.kind == .both) {
            var black_iter = graph.node_list[level_one_node].black_edges.iterator();
            while (black_iter.next()) |item| {
                const total_deg = graph.node_list[item].cardinality();
                scorer.addVisit(item, total_deg);
            }
        }
        if (options.kind == .red or options.kind == .both) {
            var red_iter = graph.node_list[level_one_node].red_edges.iterator();
            while (red_iter.next()) |item| {
                const total_deg = graph.node_list[item].cardinality();
                scorer.addVisit(item, total_deg);
            }
        }
    }
}

pub inline fn bfs_topk_ignore_large(comptime T: type, comptime K: u32, start_node: T, graph: *const Graph(T), visited: *bitset.FastBitSet, stack: *BfsQueue(T), scorer: *topk.TopKScorer(T, K), large_threshold: T, comptime options: BfsOptions) void {
    stack.clear();
    stack.addNext(start_node);
    stack.swapFrontiers();

    visited.set(start_node);
    scorer.setNextTargetDegree(graph.node_list[start_node].cardinality());

    var iter_current = stack.iterator();
    var current_level: u32 = 0;

    while (current_level < options.max_level) {
        while (iter_current.next()) |id| {
            if (current_level == 1 and graph.node_list[id].cardinality() > large_threshold) {
                continue;
            }
            if (options.kind == .black or options.kind == .both) {
                var black_iter = graph.node_list[id].black_edges.iterator();
                while (black_iter.next()) |item| {
                    const total_deg = graph.node_list[item].cardinality();
                    if (!visited.setExists(item)) {
                        stack.addNext(@intCast(T, item));
                    }
                    scorer.addVisit(item, total_deg);
                }
            }
            if (options.kind == .red or options.kind == .both) {
                var red_iter = graph.node_list[id].red_edges.iterator();
                while (red_iter.next()) |item| {
                    const total_deg = graph.node_list[item].cardinality();
                    if (!visited.setExists(item)) {
                        stack.addNext(@intCast(T, item));
                    }
                    scorer.addVisit(item, total_deg);
                }
            }
        }

        current_level += 1;
        if (stack.nextFrontierSize() == 0) {
            return;
        }
        stack.swapFrontiers();
        iter_current = stack.iterator();
    }
    while (iter_current.next()) |id| {
        const total_deg = graph.node_list[id].black_edges.cardinality() + graph.node_list[id].red_edges.cardinality();
        scorer.addVisit(id, total_deg);
    }
}

test "BFS: Find specific level" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer std.debug.assert(!gpa.deinit());
    var graph = try Graph(u8).new(10, gpa.allocator());
    defer graph.deinit();
    try graph.addEdge(0, 1);
    try graph.addEdge(0, 2);
    try graph.addEdge(0, 5);
    try graph.addEdge(5, 6);
    try graph.addEdge(6, 7);

    var visited = try bitset.FastBitSet.initEmpty(10, gpa.allocator());
    defer visited.deinit(gpa.allocator());
    var stack = try BfsQueue(u8).init(gpa.allocator(), 10);
    defer stack.deinit(gpa.allocator());

    var bfs_search = bfs(u8, 0, &graph, &visited, &stack, BfsOptions{ .max_level = 1, .kind = .black });

    try std.testing.expectEqual(bfs_search.next(), 0);
    try std.testing.expectEqual(bfs_search.next(), 1);
    try std.testing.expectEqual(bfs_search.next(), 2);
    try std.testing.expectEqual(bfs_search.next(), 5);
    try std.testing.expectEqual(bfs_search.next(), null);

    try std.testing.expectEqual(visited.cardinality, 4);

    visited.unsetAll();

    var bfs_search_red = bfs(u8, 0, &graph, &visited, &stack, BfsOptions{ .max_level = 1, .kind = .red });

    try std.testing.expectEqual(bfs_search_red.next(), 0);
    try std.testing.expectEqual(bfs_search_red.next(), null);
    try std.testing.expectEqual(visited.cardinality, 1);

    visited.unsetAll();

    var bfs_search_two_levels = bfs(u8, 0, &graph, &visited, &stack, BfsOptions{ .max_level = 2, .kind = .black });
    try std.testing.expectEqual(bfs_search_two_levels.next(), 0);
    try std.testing.expectEqual(bfs_search_two_levels.next(), 1);
    try std.testing.expectEqual(bfs_search_two_levels.next(), 2);
    try std.testing.expectEqual(bfs_search_two_levels.next(), 5);
    try std.testing.expectEqual(bfs_search_two_levels.next(), 6);
    try std.testing.expectEqual(bfs_search_two_levels.next(), null);
    try std.testing.expectEqual(visited.cardinality, 5);
}

test "BfsStack: Basic functions" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer std.debug.assert(!gpa.deinit());

    var stack = try BfsQueue(u8).init(gpa.allocator(), 100);
    defer stack.deinit(gpa.allocator());

    {
        var stack_iter = stack.iterator();
        try std.testing.expectEqual(stack_iter.next(), null);
    }

    stack.addNext(10);
    {
        var stack_iter = stack.iterator();
        try std.testing.expectEqual(stack_iter.next(), null);
    }

    stack.swapFrontiers();
    // Expect next to be cleared
    try std.testing.expectEqual(stack.nextFrontierSize(), 0);
    {
        var stack_iter = stack.iterator();
        try std.testing.expectEqual(stack_iter.next(), 10);
        try std.testing.expectEqual(stack_iter.next(), null);
    }
    stack.swapFrontiers();
    try std.testing.expectEqual(stack.nextFrontierSize(), 0);
    {
        var stack_iter = stack.iterator();
        try std.testing.expectEqual(stack_iter.next(), null);
    }
}

test "BfsStack: Clearing functions" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer std.debug.assert(!gpa.deinit());

    var stack = try BfsQueue(u8).init(gpa.allocator(), 100);
    defer stack.deinit(gpa.allocator());

    {
        var stack_iter = stack.iterator();
        try std.testing.expectEqual(stack_iter.next(), null);
    }

    stack.addNext(10);
    stack.addNext(20);
    stack.addNext(30);
    {
        var stack_iter = stack.iterator();
        try std.testing.expectEqual(stack_iter.next(), null);
    }

    stack.clear();
    stack.swapFrontiers();
    {
        var stack_iter = stack.iterator();
        try std.testing.expectEqual(stack_iter.next(), null);
    }
}
