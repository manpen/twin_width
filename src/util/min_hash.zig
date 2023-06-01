const std = @import("std");
const comptime_util = @import("comptime_checks.zig");
const compressed_bitmap = @import("compressed_bitmap.zig");
const two_level_bitset = @import("../util/two_level_bitset.zig");
const node_mod = @import("../graph/node.zig");
const contraction = @import("../tww/contraction_sequence.zig");
const retraceable_contraction = @import("../tww/retraceable_contraction_sequence.zig");
const graph_mod = @import("../graph/graph.zig");
const subgraph = @import("../graph/subgraph.zig");
const dart_hash = @import("dart_hash.zig");
const updateable_pq = @import("../util/updateable_priority_queue.zig");
const bfs_mod = @import("../graph/bfs.zig");

pub const MinHashError = error{ OutOfMemory, KeyIsMissing };

pub inline fn fisher_yates_shuffle(comptime T: type, data: []T, generator: *std.rand.DefaultPrng) void {
    comptime if (!comptime_util.checkIfIsCompatibleInteger(T)) {
        @compileError("Can only use the fisher yates shuffle on integer types u8,u16 or u32");
    };

    for (0..(data.len - 1)) |i| {
        const random = generator.next();
        const target = i + random % (data.len - i);

        const tmp = data[target];
        data[target] = data[i];
        data[i] = tmp;
    }
}

pub inline fn fisher_yates_sample_first_n(comptime T: type, data: []T, n: u32, generator: *std.rand.DefaultPrng) void {
    comptime if (!comptime_util.checkIfIsCompatibleInteger(T)) {
        @compileError("Can only use the fisher yates shuffle on integer types u8,u16 or u32");
    };

    for (0..(std.math.min(n, data.len - 1))) |i| {
        const random = generator.next();
        const target = i + random % (data.len - i);

        const tmp = data[target];
        data[target] = data[i];
        data[i] = tmp;
    }
}

pub inline fn circular_permutation_shift(comptime T: type, data: []T) void {
    // See https://openreview.net/pdf?id=NrkAAcMpRoT

    const item = data[data.len - 1];

    var next = data[0];

    for (0..(data.len - 1)) |i| {
        const tmp = data[i + 1];
        data[i + 1] = next;
        next = tmp;
    }

    data[0] = item;
}

pub const MinHash = struct {
    const Self = @This();

    const PRIME: u32 = 4_294_967_291; // = 2^32 - 5
    //const PRIME: u32 = 2_147_483_647; // = 2^31 - 1

    shift: u32,
    mult: u32,

    pub fn sample_hash_function(generator: *std.rand.Random, n: u32) Self {
        _ = n; // may be needed for other hash functions

        var shift = generator.intRangeAtMost(u32, 0, PRIME - 1);
        var mult = generator.intRangeAtMost(u32, 1, PRIME - 1);
        return Self{ .shift = shift, .mult = mult };
    }

    pub fn invalid() Self {
        return Self{ .shift = 0, .mult = 0 };
    }

    pub inline fn getHash(self: *const Self, item: u32) u32 {
        std.debug.assert(self.mult != 0);
        const hashed = (@intCast(u64, self.mult) * item + self.shift) % PRIME;
        return @intCast(u32, hashed);
    }

    pub inline fn hash(self: *const Self, comptime InputIteratorType: type, iter: InputIteratorType) u32 {
        var min: u32 = std.math.maxInt(u32);
        while (iter.next()) |item| {
            min = std.math.min(min, self.getHash(item));
        }
        return min;
    }
};

const TestIterator = struct {
    const Self = @This();
    base: *std.BoundedArray(u32, 10),
    index: u32,
    offset: u32,
    pub fn next(self: *Self) ?u32 {
        if (self.base.len <= self.index) return null;
        const item = self.base.buffer[self.index];
        self.index += 1;
        return item + self.offset;
    }
};

pub fn hashFn(ctx: MinHashBandContext, value: []u32) u64 {
    _ = ctx;
    return std.hash.Wyhash.hash(0, std.mem.sliceAsBytes(value));
}

pub fn eqFn(ctx: MinHashBandContext, rhs: []u32, lhs: []u32) bool {
    _ = ctx;
    return std.mem.eql(u32, rhs, lhs);
}

pub const MinHashBandContext = struct {
    pub const hash = hashFn;
    pub const eql = eqFn;
};

pub fn MinHashBand(comptime B: u32) type {
    return struct {
        const Self = @This();
        hash_functions: []MinHash,
        hash_cache: []u32,
        collisions: std.HashMapUnmanaged([]u32, std.AutoArrayHashMapUnmanaged(u32, void), MinHashBandContext, 80),
        allocator: std.mem.Allocator,
        downgraded_level: u32,
        initialized_nodes: std.bit_set.DynamicBitSet,

        pub fn clear(self: *Self) void {
            var iter = self.collisions.iterator();
            while (iter.next()) |item| {
                item.value_ptr.deinit(self.allocator);
                self.allocator.free(item.key_ptr.*);
            }
            for (0..self.initialized_nodes.capacity()) |index| {
                self.initialized_nodes.unset(index);
            }
            self.collisions.clearRetainingCapacity();
        }

        pub fn randomize_functions(self: *Self, generator: *std.rand.Random) void {
            for (self.hash_functions) |*h| {
                h.* = MinHash.sample_hash_function(generator, 0);
            }
        }

        pub inline fn init(allocator: std.mem.Allocator, cache_size: u32) !Self {
            var hashes = try allocator.alloc(MinHash, B);
            for (hashes) |*h| {
                h.* = MinHash.invalid();
            }

            var collisions = std.HashMapUnmanaged([]u32, std.AutoArrayHashMapUnmanaged(u32, void), MinHashBandContext, 80){};

            var hash_cache = try allocator.alloc(u32, cache_size * B);
            var bitset = try std.bit_set.DynamicBitSet.initEmpty(allocator, cache_size);

            return Self{ .downgraded_level = 0, .hash_functions = hashes, .collisions = collisions, .allocator = allocator, .hash_cache = hash_cache, .initialized_nodes = bitset };
        }

        pub inline fn rehashItem(self: *Self, comptime IterType: type, comptime Context: type, key: u32, iterator: IterType, context: Context, callback_removed: fn (Context, u32, ?*std.AutoArrayHashMapUnmanaged(u32, void)) MinHashError!void, callback_added: fn (Context, u32, ?*std.AutoArrayHashMapUnmanaged(u32, void)) MinHashError!void) !void {
            const first = try self.removeItem(key);
            try callback_removed(context, key, first);

            const next = try self.addItem(IterType, key, iterator);
            try callback_added(context, key, next);
        }

        pub inline fn rehashItemNoRemove(self: *Self, comptime IterType: type, comptime Context: type, key: u32, iterator: IterType, context: Context, callback_added: fn (Context, u32, ?*std.AutoArrayHashMapUnmanaged(u32, void)) MinHashError!void) !void {
            const next = try self.addItem(IterType, key, iterator);
            try callback_added(context, key, next);
        }

        pub inline fn updateItem(self: *Self, comptime IterType: type, comptime Context: type, key: u32, iterator: IterType, context: Context, callback_removed: fn (Context, u32, ?*std.AutoArrayHashMapUnmanaged(u32, void)) MinHashError!void, callback_added: fn (Context, u32, ?*std.AutoArrayHashMapUnmanaged(u32, void)) MinHashError!void, removed_feature: ?u32, added_feature: ?u32) !void {
            if (!self.initialized_nodes.isSet(key)) {
                return try callback_added(context, key, try self.addItem(IterType, key, iterator));
            }
            var removed: bool = false;
            for (0..B) |i| {
                const min_before = self.hash_cache[key * B + i];

                if (removed_feature) |rm| {
                    if (self.hash_functions[i].getHash(rm) == min_before) {
                        if (!removed) {
                            return self.rehashItem(IterType, Context, key, iterator, context, callback_removed, callback_added);
                        } else {
                            return self.rehashItemNoRemove(IterType, Context, key, iterator, context, callback_added);
                        }
                    }
                }

                if (added_feature) |feature| {
                    const hash = self.hash_functions[i].getHash(feature);
                    if (!removed and min_before > hash) {
                        removed = true;
                        const first = try self.removeItem(key);
                        try callback_removed(context, key, first);
                    }
                    self.hash_cache[key * B + i] = std.math.min(self.hash_cache[key * B + i], hash);
                }
            }

            if (removed) {
                try callback_added(context, key, try self.addItemFromCache(key));
            }
        }

        pub inline fn downgradeBand(self: *Self, hit_map: *std.AutoHashMapUnmanaged(u64, u32)) !void {
            if (self.downgraded_level == B - 1) return;
            self.downgraded_level += 1;
            var second_collisions_map = std.HashMapUnmanaged([]u32, std.AutoArrayHashMapUnmanaged(u32, void), MinHashBandContext, 80){};
            try second_collisions_map.ensureTotalCapacity(self.allocator, self.collisions.count());
            var iterator = self.collisions.iterator();
            while (iterator.next()) |item| {
                {
                    var hits_iter_outer = item.value_ptr.iterator();
                    while (hits_iter_outer.next()) |inner_item| {
                        var hits_iter_inner = item.value_ptr.iterator();
                        while (hits_iter_inner.next()) |inner_item_second| {
                            if (inner_item.key_ptr.* == inner_item_second.key_ptr.*) continue;
                            const unique_key = MinHashSimilarity(u32, B).calculateUniqueKey(inner_item.key_ptr.*, inner_item_second.key_ptr.*);
                            if (hit_map.getPtr(unique_key)) |v| {
                                v.* += 1;
                            } else {
                                try hit_map.put(self.allocator, unique_key, 1);
                            }
                        }
                    }
                }
                if (second_collisions_map.getPtr(item.key_ptr.*[0..(B - self.downgraded_level)])) |hits| {
                    var hits_iter_outer = hits.iterator();
                    while (hits_iter_outer.next()) |inner_item| {
                        var hits_iter_inner = item.value_ptr.iterator();
                        while (hits_iter_inner.next()) |inner_item_second| {
                            if (inner_item.key_ptr.* == inner_item_second.key_ptr.*) continue;
                            const unique_key = MinHashSimilarity(u32, B).calculateUniqueKey(inner_item.key_ptr.*, inner_item_second.key_ptr.*);
                            if (hit_map.getPtr(unique_key)) |v| {
                                v.* += 1;
                            } else {
                                try hit_map.put(self.allocator, unique_key, 1);
                            }
                        }
                    }

                    var iter_inner = item.value_ptr.iterator();
                    while (iter_inner.next()) |inner_item| {
                        try hits.put(self.allocator, inner_item.key_ptr.*, {});
                    }

                    item.value_ptr.deinit(self.allocator);
                } else {
                    try second_collisions_map.put(self.allocator, item.key_ptr.*[0..(B - self.downgraded_level)], item.value_ptr.*);
                }
            }
            self.collisions.deinit(self.allocator);
            self.collisions = second_collisions_map;
        }

        inline fn addItemFromCache(self: *Self, key: u32) !?*std.AutoArrayHashMapUnmanaged(u32, void) {
            if (self.collisions.getPtr(self.getSignature(key))) |result| {
                try result.put(self.allocator, key, {});
                return result;
            } else {
                var copy = try self.allocator.alloc(u32, B);
                @memcpy(copy, self.hash_cache[key * B .. key * B + B]);

                var new_list = std.AutoArrayHashMapUnmanaged(u32, void){};
                try new_list.put(self.allocator, key, {});
                try self.collisions.put(self.allocator, copy[0..(B - self.downgraded_level)], new_list);
            }
            return null;
        }

        // Returns the collisions if there were any
        pub inline fn addItem(self: *Self, comptime IterType: type, key: u32, iterator: IterType) !?*std.AutoArrayHashMapUnmanaged(u32, void) {
            self.initialized_nodes.set(key);
            @memset(self.hash_cache[key * B .. key * B + B], std.math.maxInt(u32));
            while (iterator.next()) |item| {
                for (0..B) |i| {
                    self.hash_cache[key * B + i] = std.math.min(self.hash_cache[key * B + i], self.hash_functions[i].getHash(item));
                }
            }
            return self.addItemFromCache(key);
        }

        pub inline fn getSignature(self: *const Self, key: u32) []u32 {
            return self.hash_cache[key * B .. key * B + (B - self.downgraded_level)];
        }

        pub inline fn removeItem(self: *Self, key: u32) MinHashError!?*std.AutoArrayHashMapUnmanaged(u32, void) {
            if (!self.initialized_nodes.isSet(key)) return null;
            if (self.collisions.getPtr(self.getSignature(key))) |ptr| {
                if (ptr.count() > 1) {
                    _ = ptr.swapRemove(key);
                    return ptr;
                }
                // Release owned memory
                ptr.deinit(self.allocator);
                if (self.collisions.fetchRemove(self.getSignature(key))) |v| {
                    self.allocator.free(v.key);
                }
            } else {
                return MinHashError.KeyIsMissing;
            }
            return null;
        }

        pub fn deinit(self: *Self) void {
            self.initialized_nodes.deinit();
            var iter = self.collisions.iterator();
            while (iter.next()) |item| {
                item.value_ptr.deinit(self.allocator);
                self.allocator.free(item.key_ptr.*);
            }
            self.collisions.clearAndFree(self.allocator);
            self.allocator.free(self.hash_cache);
            self.allocator.free(self.hash_functions);
        }
    };
}

pub const SimilarityPriority = struct {
    cardinality: u32,
    score: u32,

    pub fn fromScoreAndNodes(comptime T: type, a: u32, b: u32, score: u32, graph: *graph_mod.Graph(T)) SimilarityPriority {
        var max: u32 = std.math.max(graph.node_list[a].red_edges.cardinality(), graph.node_list[b].red_edges.cardinality());
        return .{
            .cardinality = max,
            .score = score,
        };
    }

    pub fn compare(ctx: void, lhs: SimilarityPriority, rhs: SimilarityPriority) std.math.Order {
        _ = ctx;
        if (rhs.score == lhs.score) {
            return std.math.order(lhs.cardinality, rhs.cardinality);
        }
        return std.math.order(rhs.score, lhs.score);
    }
};

pub fn MinHashSimilarity(comptime T: type, comptime B: u32) type {
    return struct {
        const Self = @This();
        pub const PriorityItem = struct {
            key: u64,
            prio: SimilarityPriority,
        };

        bands: []MinHashBand(B),
        hit_map: std.AutoHashMapUnmanaged(u64, u32),
        sim_pq: updateable_pq.UpdateablePriorityQueue(u64, SimilarityPriority, void, SimilarityPriority.compare),

        number_of_nodes: u32,
        allocator: std.mem.Allocator,
        graph: *graph_mod.Graph(T),
        canidate_count: u32,
        canidate_list: std.ArrayListUnmanaged(PriorityItem),

        const NodeIteratorSplitRedBlack = struct {
            number_of_nodes: u32,
            iter: graph_mod.Graph(T).NodeType.UnorderedNodeEdgeIterator,
            final_node: ?u32 = null,
            pub fn next(self: *NodeIteratorSplitRedBlack) ?u32 {
                if (self.iter.next()) |item| {
                    if (self.iter.red) return item + self.number_of_nodes;
                    return item;
                }
                if (self.final_node) |item| {
                    const copy = item;
                    self.final_node = null;
                    return copy;
                }
                return null;
            }
        };

        pub inline fn keyIntoMove(key: u64) contraction.Contraction(T) {
            var first = @intCast(u32, key >> 32); // / number_of_nodes;
            var second = @truncate(u32, key);

            return contraction.Contraction(T){ .erased = @intCast(T, first), .survivor = @intCast(T, second) };
        }

        const NodeIterator = struct {
            number_of_nodes: u32,
            iter: graph_mod.Graph(T).NodeType.EdgeIterType,
            pub fn next(self: *NodeIterator) ?u32 {
                if (self.iter.next()) |item| {
                    return item;
                }
                return null;
            }
        };

        pub inline fn calculateUniqueKey(a: u32, b: u32) u64 {
            if (a < b) {
                return (@intCast(u64, b) << 32) | a;
            } else {
                return (@intCast(u64, a) << 32) | b;
            }
        }

        fn removedCallback(self: *Self, key: u32, list: ?*std.AutoArrayHashMapUnmanaged(u32, void)) MinHashError!void {
            if (list) |l| {
                var i: usize = 0;
                while (i < l.count()) {
                    const partner = l.keys()[i];
                    if (partner == key) {
                        i += 1;
                        continue;
                    }
                    const unique_key = Self.calculateUniqueKey(partner, key);
                    if (self.graph.erased_nodes.get(partner) or self.graph.erased_nodes.get(key)) {
                        _ = self.hit_map.remove(unique_key);
                        _ = self.sim_pq.removeKey(unique_key);
                        i += 1;
                        continue;
                    }

                    if (self.hit_map.getPtr(unique_key)) |pt| {
                        i += 1;
                        pt.* -= 1;
                        if (pt.* == 0) {
                            _ = self.hit_map.remove(unique_key);
                            _ = self.sim_pq.removeKey(unique_key);
                        } else {
                            _ = self.sim_pq.updateOrInsert(unique_key, SimilarityPriority.fromScoreAndNodes(T, key, partner, pt.*, self.graph)) catch return error.OutOfMemory;
                        }
                        // Untrack this node now!
                    } else {
                        return error.OutOfMemory;
                    }
                }
            }
        }

        fn addedCallback(self: *Self, key: u32, list: ?*std.AutoArrayHashMapUnmanaged(u32, void)) MinHashError!void {
            if (list) |l| {
                var i: usize = 0;
                while (i < l.count()) {
                    const partner = l.keys()[i];
                    i += 1;
                    if (partner == key) {
                        continue;
                    }
                    if (self.graph.erased_nodes.get(partner)) continue;
                    const unique_key = Self.calculateUniqueKey(partner, key);
                    if (self.hit_map.getPtr(unique_key)) |pt| {
                        pt.* += 1;
                        _ = self.sim_pq.updateOrInsert(unique_key, SimilarityPriority.fromScoreAndNodes(T, key, partner, pt.*, self.graph)) catch return error.OutOfMemory;
                    } else {
                        self.hit_map.put(self.allocator, unique_key, 1) catch return error.OutOfMemory;
                        _ = self.sim_pq.updateOrInsert(unique_key, SimilarityPriority.fromScoreAndNodes(T, key, partner, 1, self.graph)) catch return error.OutOfMemory;
                    }
                }
            }
        }

        pub fn removeNode(self: *Self, node: u32) !void {
            for (0..self.bands.len) |index| {
                const result = try self.bands[index].removeItem(node);
                try Self.removedCallback(self, node, result);
            }
        }

        pub fn rehashNode(self: *Self, node: u32, graph: *graph_mod.Graph(T)) !void {
            const split = self.bands.len / 3;
            for (0..split) |index| {
                var iter = NodeIteratorSplitRedBlack{ .number_of_nodes = graph.number_of_nodes, .iter = graph.node_list[node].unorderedIterator() };

                try self.bands[index].rehashItem(*NodeIteratorSplitRedBlack, @TypeOf(self), node, &iter, self, Self.removedCallback, Self.addedCallback);
            }
            for (split..split * 2) |index| {
                var iter = NodeIteratorSplitRedBlack{ .number_of_nodes = graph.number_of_nodes, .iter = graph.node_list[node].unorderedIterator(), .final_node = node };

                try self.bands[index].rehashItem(*NodeIteratorSplitRedBlack, @TypeOf(self), node, &iter, self, Self.removedCallback, Self.addedCallback);
            }
            for (split * 2..self.bands.len) |index| {
                var iter = NodeIteratorSplitRedBlack{ .number_of_nodes = graph.number_of_nodes, .iter = graph.node_list[node].unorderedIterator(), .final_node = node + self.number_of_nodes };

                try self.bands[index].rehashItem(*NodeIteratorSplitRedBlack, @TypeOf(self), node, &iter, self, Self.removedCallback, Self.addedCallback);
            }
        }

        pub const ChangedEdge = struct {
            removed: u32,
            red: bool = false,
            added: ?u32 = null,
            added_red: bool = false,
        };

        pub fn changedEdge(self: *Self, node: u32, graph: *graph_mod.Graph(T), changed: ChangedEdge) !void {
            const removed = if (changed.red) changed.removed + self.number_of_nodes else changed.removed;
            const added = if (changed.added) |a| if (changed.added_red) a + self.number_of_nodes else changed.added else changed.added;

            const split = self.bands.len / 3;
            for (0..split) |index| {
                var iter = NodeIteratorSplitRedBlack{ .number_of_nodes = graph.number_of_nodes, .iter = graph.node_list[node].unorderedIterator() };
                try self.bands[index].updateItem(*NodeIteratorSplitRedBlack, @TypeOf(self), node, &iter, self, Self.removedCallback, Self.addedCallback, removed, added);
            }
            for (split..split * 2) |index| {
                var iter = NodeIteratorSplitRedBlack{ .number_of_nodes = graph.number_of_nodes, .iter = graph.node_list[node].unorderedIterator(), .final_node = node };
                try self.bands[index].updateItem(*NodeIteratorSplitRedBlack, @TypeOf(self), node, &iter, self, Self.removedCallback, Self.addedCallback, removed, added);
            }
            for (split * 2..self.bands.len) |index| {
                var iter = NodeIteratorSplitRedBlack{ .number_of_nodes = graph.number_of_nodes, .iter = graph.node_list[node].unorderedIterator(), .final_node = node + self.number_of_nodes };
                try self.bands[index].updateItem(*NodeIteratorSplitRedBlack, @TypeOf(self), node, &iter, self, Self.removedCallback, Self.addedCallback, removed, added);
            }
        }

        pub fn bootstrapNodes(self: *Self, nodes: []T, graph: *graph_mod.Graph(T), seed: u64, visited: *two_level_bitset.FastBitSet, stack: *bfs_mod.BfsQueue(T)) !void {
            var rng = std.rand.DefaultPrng.init(seed);
            var random = rng.random();

            self.graph = graph;
            for (self.bands) |*b| {
                b.clear();
                b.randomize_functions(&random);
            }

            self.hit_map.clearRetainingCapacity();
            while (self.sim_pq.removeOrNull()) |_| {}

            const split = self.bands.len / 3;

            visited.unsetAll();
            stack.clear();
            if (nodes.len == 0) return;
            var bfs = bfs_mod.bfs(T, nodes[0], graph, visited, stack, .{ .kind = .both, .max_level = std.math.maxInt(T) });
            while (bfs.next()) |node| {
                if (graph.erased_nodes.get(node)) continue;
                if (self.hit_map.count() > 10_000) break;
                for (0..split) |index| {
                    var iter = NodeIteratorSplitRedBlack{ .number_of_nodes = graph.number_of_nodes, .iter = graph.node_list[node].unorderedIterator() };
                    try self.bands[index].rehashItemNoRemove(*NodeIteratorSplitRedBlack, @TypeOf(self), node, &iter, self, Self.addedCallback);
                }
                for (split..split * 2) |index| {
                    var iter = NodeIteratorSplitRedBlack{ .number_of_nodes = graph.number_of_nodes, .iter = graph.node_list[node].unorderedIterator(), .final_node = node };
                    try self.bands[index].rehashItemNoRemove(*NodeIteratorSplitRedBlack, @TypeOf(self), node, &iter, self, Self.addedCallback);
                }
                for (split * 2..self.bands.len) |index| {
                    var iter = NodeIteratorSplitRedBlack{ .number_of_nodes = graph.number_of_nodes, .iter = graph.node_list[node].unorderedIterator(), .final_node = node + self.number_of_nodes };
                    try self.bands[index].rehashItemNoRemove(*NodeIteratorSplitRedBlack, @TypeOf(self), node, &iter, self, Self.addedCallback);
                }
            }
            visited.unsetAll();
            stack.clear();
        }

        pub const TwinWidthPriorityContraction = struct {
            tww: T,
            merge: contraction.Contraction(T),
            pub fn compare(ctx: void, lhs: TwinWidthPriorityContraction, rhs: TwinWidthPriorityContraction) std.math.Order {
                _ = ctx;
                return std.math.order(lhs.tww, rhs.tww);
            }
            pub fn compareInverse(ctx: void, lhs: T, rhs: T) std.math.Order {
                _ = ctx;
                return std.math.order(rhs, lhs);
            }
        };

        pub fn getBestMove(self: *Self, graph: *graph_mod.Graph(T), current_tww: T) !?contraction.Contraction(T) {
            var min_cont: ?contraction.Contraction(T) = null;

            self.canidate_list.clearRetainingCapacity();

            if (self.sim_pq.len < self.canidate_count * 10) {
                self.hit_map.clearRetainingCapacity();
                // Clear queue
                while (self.sim_pq.removeOrNull()) |_| {}

                for (0..self.bands.len) |index| {
                    try self.bands[index].downgradeBand(&self.hit_map);
                }

                var hit_iter = self.hit_map.iterator();
                while (hit_iter.next()) |hit| {
                    const mv = Self.keyIntoMove(hit.key_ptr.*);
                    const prio = SimilarityPriority.fromScoreAndNodes(T, mv.erased, mv.survivor, hit.value_ptr.*, self.graph);
                    _ = try self.sim_pq.updateOrInsert(hit.key_ptr.*, prio);
                }
            }

            var best_tww = graph_mod.Graph(T).InducedTwinWidthPotential.default();
            while (self.sim_pq.removeOrNull()) |it| {
                const mv = Self.keyIntoMove(it.key);
                if (self.graph.erased_nodes.get(mv.erased) or self.graph.erased_nodes.get(mv.survivor)) continue;
                if (self.hit_map.get(it.key)) |k| {
                    const prio = SimilarityPriority.fromScoreAndNodes(T, mv.erased, mv.survivor, k, self.graph);
                    if (prio.cardinality != it.priority.cardinality) {
                        _ = try self.sim_pq.updateOrInsert(it.key, prio);
                        continue;
                    }
                    var tww = graph.calculateInducedTwwPotential(mv.erased, mv.survivor, &best_tww, current_tww);

                    try self.canidate_list.append(self.allocator, .{ .key = it.key, .prio = it.priority });

                    if (tww.isLess(best_tww, current_tww)) {
                        min_cont = .{ .erased = mv.erased, .survivor = mv.survivor };
                        best_tww = tww;
                    }
                    if (self.canidate_list.items.len == self.canidate_count) break;
                }
            }

            for (self.canidate_list.items) |it| {
                _ = try self.sim_pq.updateOrInsert(it.key, it.prio);
            }
            return min_cont;
        }

        pub fn init(allocator: std.mem.Allocator, comptime num_bands: u32, number_of_nodes: u32) !Self {
            var bands = try allocator.alloc(MinHashBand(B), num_bands);

            comptime if (num_bands % 3 != 0) @compileError("Number of bands must be divisible by 3!");
            var counter: u32 = 0;
            for (bands) |*band| {
                band.* = try MinHashBand(B).init(allocator, number_of_nodes);
                counter += number_of_nodes * 2;
            }

            var hit_map = std.AutoHashMapUnmanaged(u64, u32){};

            var sim_pq = updateable_pq.UpdateablePriorityQueue(u64, SimilarityPriority, void, SimilarityPriority.compare).init(allocator, {});

            return Self{
                .bands = bands, //
                .hit_map = hit_map,
                .allocator = allocator,
                .number_of_nodes = number_of_nodes,
                .graph = undefined,
                .sim_pq = sim_pq,
                .canidate_count = 20,
                .canidate_list = try std.ArrayListUnmanaged(PriorityItem).initCapacity(allocator, number_of_nodes),
            };
        }

        pub fn deinit(self: *Self) void {
            self.canidate_list.deinit(self.allocator);

            self.hit_map.deinit(self.allocator);
            for (0..self.bands.len) |i| {
                self.bands[self.bands.len - i - 1].deinit();
            }

            self.allocator.free(self.bands);
        }
    };
}

test "MinHash: Check is permutation" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer std.debug.assert(gpa.deinit() == .ok);
    var perm = try MinHash.generatePermutation(gpa.allocator(), 200, 19);
    defer gpa.allocator().free(perm);
    var min_hash = MinHash.init(perm, 0);
    var min_hash2 = MinHash.init(perm, 1);

    var hash_map = std.AutoHashMap(u32, void).init(gpa.allocator());
    defer hash_map.deinit();

    for (min_hash.permutation) |item| {
        if (hash_map.contains(item)) try std.testing.expectEqual(true, false);
        try hash_map.put(item, {});
    }

    try std.testing.expectEqual(min_hash.permutation.len, 200);

    try std.testing.expectEqual(min_hash.getHash(101), min_hash2.getHash(100));
}

test "MinHash: Check hashing of complex types" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer std.debug.assert(gpa.deinit() == .ok);

    var perm = try MinHash.generatePermutation(gpa.allocator(), 200, 19);
    defer gpa.allocator().free(perm);
    var min_hash = MinHash.init(perm, 0);

    var hash_map = std.AutoHashMap(u32, void).init(gpa.allocator());
    defer hash_map.deinit();

    var initial = try std.BoundedArray(u32, 10).init(0);
    try initial.append(10);
    try initial.append(20);

    var iter_black = TestIterator{ .base = &initial, .index = 0, .offset = 0 };
    var iter_red = TestIterator{ .base = &initial, .index = 0, .offset = 100 };

    var hash = min_hash.hash(@TypeOf(&iter_black), &iter_black);
    hash = std.math.min(min_hash.hash(@TypeOf(&iter_red), &iter_red), hash);

    var min: u32 = std.math.maxInt(u32);
    for (initial.buffer[0..initial.len]) |item| {
        min = std.math.min(min, min_hash.permutation[item]);
    }
    for (initial.buffer[0..initial.len]) |item| {
        min = std.math.min(min, min_hash.permutation[item + 100]);
    }

    try std.testing.expectEqual(min, hash);
}

test "MinHashBand: Basics" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer std.debug.assert(gpa.deinit() == .ok);

    var perm = try MinHash.generatePermutation(gpa.allocator(), 200, 19);
    defer gpa.allocator().free(perm);

    var band = try MinHashBand(4).init(gpa.allocator(), 100, perm, 0);
    defer band.deinit();

    var initial = try std.BoundedArray(u32, 10).init(0);
    try initial.append(10);
    try initial.append(20);

    var iter_black = TestIterator{ .base = &initial, .index = 0, .offset = 0 };
    try std.testing.expectEqual(try band.addItem(@TypeOf(&iter_black), 1, &iter_black), null);

    var iter_black2 = TestIterator{ .base = &initial, .index = 0, .offset = 0 };
    var result = try band.addItem(@TypeOf(&iter_black2), 2, &iter_black2);

    try std.testing.expect(result != null);
    try std.testing.expectEqual(result.?.items.len, 2);
    try std.testing.expectEqual(result.?.items[0], 1);
    try std.testing.expectEqual(result.?.items[1], 2);

    var result2 = try band.removeItem(2);

    try std.testing.expect(result2 != null);
    try std.testing.expectEqual(result2.?.items.len, 1);
    try std.testing.expectEqual(result2.?.items[0], 1);
}
