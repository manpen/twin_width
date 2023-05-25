const std = @import("std");
const comptime_util = @import("comptime_checks.zig");
const compressed_bitmap = @import("compressed_bitmap.zig");
const two_level_bitset = @import("../util/two_level_bitset.zig");
const node_mod = @import("../graph/node.zig");
const contraction = @import("../tww/contraction_sequence.zig");
const graph_mod = @import("../graph/graph.zig");

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
    permutation: []u32,

    pub inline fn permutate(slice: []u32, seed: u64) void {
        var gen = std.rand.DefaultPrng.init(seed);
        fisher_yates_shuffle(u32, slice, &gen);
    }

    pub inline fn generatePermutation(allocator: std.mem.Allocator, length: u32) ![]u32 {
        var memory = try allocator.alloc(u32, length);
        for (0..memory.len) |i| {
            memory[i] = @intCast(u32, i);
        }
        return memory;
    }

    pub inline fn init(permutation: []u32) Self {
        return Self{
            .permutation = permutation,
        };
    }

    pub inline fn getHash(self: *const Self, item: u32) u32 {
        return self.permutation[item];
    }

    pub inline fn hash(self: *const Self, comptime InputIteratorType: type, iter: InputIteratorType) u32 {
        var min: u32 = std.math.maxInt(u32);
        while (iter.next()) |item| {
            min = std.math.min(min, self.getHash(item));
        }
        return min;
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
        erased_set: std.AutoHashMapUnmanaged(u32, void),

        pub fn clear(self: *Self) void {
            var iter = self.collisions.iterator();
            while (iter.next()) |item| {
                item.value_ptr.deinit(self.allocator);
                self.allocator.free(item.key_ptr.*);
            }
            self.collisions.clearRetainingCapacity();
            self.erased_set.clearRetainingCapacity();
        }

        pub inline fn init(allocator: std.mem.Allocator, cache_size: u32, permutation: []u32) !Self {
            var hashes = try allocator.alloc(MinHash, B);

            for (hashes) |*hash| {
                hash.* = MinHash.init(permutation);
            }

            var collisions = std.HashMapUnmanaged([]u32, std.AutoArrayHashMapUnmanaged(u32, void), MinHashBandContext, 80){};

            var hash_cache = try allocator.alloc(u32, cache_size * B);

            var erased_set = std.AutoHashMapUnmanaged(u32, void){};

            return Self{ .hash_functions = hashes, .collisions = collisions, .allocator = allocator, .hash_cache = hash_cache, .erased_set = erased_set };
        }

        pub inline fn rehashItem(self: *Self, comptime IterType: type, comptime Context: type, key: u32, iterator: IterType, context: Context, callback_removed: fn (Context, u32, ?*std.AutoArrayHashMapUnmanaged(u32, void)) void, callback_added: fn (Context, u32, ?*std.AutoArrayHashMapUnmanaged(u32, void)) void) !void {
            const first = try self.removeItem(key);
            callback_removed(context, key, first);

            const next = try self.addItem(IterType, key, iterator);
            callback_added(context, key, next);
        }

        pub inline fn rehashItemNoRemove(self: *Self, comptime IterType: type, comptime Context: type, key: u32, iterator: IterType, context: Context, callback_added: fn (Context, u32, ?*std.AutoArrayHashMapUnmanaged(u32, void)) void) !void {
            const next = try self.addItem(IterType, key, iterator);
            callback_added(context, key, next);
        }

        pub inline fn updateItem(self: *Self, comptime IterType: type, comptime Context: type, key: u32, iterator: IterType, context: Context, callback_removed: fn (Context, u32, ?*std.AutoArrayHashMapUnmanaged(u32, void)) void, callback_added: fn (Context, u32, ?*std.AutoArrayHashMapUnmanaged(u32, void)) void, removed_feature: ?u32, added_feature: ?u32) !void {
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
                        callback_removed(context, key, first);
                    }
                    self.hash_cache[key * B + i] = std.math.min(self.hash_cache[key * B + i], hash);
                }
            }

            if (removed) {
                callback_added(context, key, try self.addItemFromCache(key));
            }
        }

        inline fn addItemFromCache(self: *Self, key: u32) !?*std.AutoArrayHashMapUnmanaged(u32, void) {
            if (self.collisions.getPtr(self.hash_cache[key * B .. key * B + B])) |result| {
                try result.put(self.allocator, key, {});
                return result;
            } else {
                var copy = try self.allocator.alloc(u32, B);
                @memcpy(copy, self.hash_cache[key * B .. key * B + B]);

                var new_list = std.AutoArrayHashMapUnmanaged(u32, void){};
                try new_list.put(self.allocator, key, {});
                try self.collisions.put(self.allocator, copy, new_list);
            }
            return null;
        }

        // Returns the collisions if there were any
        pub inline fn addItem(self: *Self, comptime IterType: type, key: u32, iterator: IterType) !?*std.AutoArrayHashMapUnmanaged(u32, void) {
            @memset(self.hash_cache[key * B .. key * B + B], std.math.maxInt(u32));
            while (iterator.next()) |item| {
                for (0..B) |i| {
                    self.hash_cache[key * B + i] = std.math.min(self.hash_cache[key * B + i], self.hash_functions[i].getHash(item));
                }
            }
            return self.addItemFromCache(key);
        }

        pub inline fn getSignature(self: *const Self, key: u32) []u32 {
            return self.hash_cache[key * B .. key * B + B];
        }

        pub inline fn removeItem(self: *Self, key: u32) !?*std.AutoArrayHashMapUnmanaged(u32, void) {
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
                std.debug.print("Key {}\n", .{key});
                @panic("Every key should be represented in the map!");
            }
            return null;
        }

        pub fn deinit(self: *Self) void {
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

pub const MinHashEntry = struct {
    const Self = @This();
    key: u64,
    score: u32,

    pub fn fromKey(key: u64) MinHashEntry {
        return MinHashEntry{
            .key = key,
            .score = 0,
        };
    }

    pub inline fn calculateScore(comptime T: type, key: u64, score: u32, graph: *graph_mod.Graph(T)) u32 {
        const mv = MinHashEntry.keyIntoMove(T, key, graph.number_of_nodes);
        const combined_card = @intCast(u32, graph.node_list[mv.erased].cardinality()) + @intCast(u32, graph.node_list[mv.survivor].cardinality());
        const tww_nb = graph.calculateMaxTwwScore(mv.erased, mv.survivor).tww_nb;
        const score_calc: u32 = (combined_card - (score * combined_card) / @intCast(u32, graph.min_hash.bands.len));
        _ = score_calc;
        return tww_nb;
    }

    pub fn from(comptime T: type, key: u64, score: u32, graph: *graph_mod.Graph(T)) MinHashEntry {
        return MinHashEntry{ .key = key, .score = Self.calculateScore(T, key, score, graph) };
    }

    pub inline fn intoMove(self: *const Self, comptime T: type, number_of_nodes: u32) contraction.Contraction(T) {
        var first = self.key / number_of_nodes;
        var second = self.key % number_of_nodes;
        return contraction.Contraction(T){ .erased = @intCast(T, first), .survivor = @intCast(T, second) };
    }

    pub inline fn keyIntoMove(comptime T: type, key: u64, number_of_nodes: u32) contraction.Contraction(T) {
        var first = key / number_of_nodes;
        var second = key % number_of_nodes;
        return contraction.Contraction(T){ .erased = @intCast(T, first), .survivor = @intCast(T, second) };
    }

    pub fn compare(ctx: void, lhs: Self, rhs: Self) std.math.Order {
        _ = ctx;
        if (lhs.key == rhs.key) return std.math.Order.eq;
        if (lhs.score < rhs.score) return std.math.Order.lt;
        return std.math.Order.gt;
    }
};

pub const SmallDegreeNode = struct {
    id: u32,
    cardinality: u32,

    pub fn compare(ctx: void, lhs: SmallDegreeNode, rhs: SmallDegreeNode) std.math.Order {
        _ = ctx;
        if (lhs.id == rhs.id) return std.math.Order.eq;
        if (lhs.cardinality < rhs.cardinality) return std.math.Order.lt;
        return std.math.Order.gt;
    }
};

pub fn MinHashSimilarity(comptime T: type, comptime B: u32) type {
    return struct {
        const Self = @This();

        bands: []MinHashBand(B),
        permutation: []u32,
        hit_map: std.AutoHashMapUnmanaged(u64, u32),

        tww_nb: std.AutoArrayHashMapUnmanaged(u64, u32),
        update_list: std.AutoHashMapUnmanaged(u32, std.AutoArrayHashMapUnmanaged(u32, void)),

        similarity_cardinality: []std.AutoArrayHashMapUnmanaged(u64, void),
        similarity_threshold: u32,

        number_of_nodes: u32,
        allocator: std.mem.Allocator,
        graph: *graph_mod.Graph(T),

        const NodeIteratorWithSelf = struct {
            number_of_nodes: u32,
            iter: graph_mod.Graph(T).NodeType.EdgeIterType,
            self: u32,
            done: bool = false,
            pub fn next(self: *NodeIteratorWithSelf) ?u32 {
                if (self.iter.next()) |item| {
                    return item;
                }
                if (!self.done) {
                    self.done = true;
                    return self.self;
                }

                return null;
            }
        };

        const NodeIteratorSplitRedBlack = struct {
            number_of_nodes: u32,
            iter: graph_mod.Graph(T).NodeType.UnorderedNodeEdgeIterator,
            pub fn next(self: *NodeIteratorSplitRedBlack) ?u32 {
                if (self.iter.next()) |item| {
                    if (self.iter.red) return item + self.number_of_nodes;
                    return item;
                }
                return null;
            }
        };

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

        pub inline fn calculateUniqueKey(a: u32, b: u32, number_of_nodes: u32) u64 {
            if (a < b) return @intCast(u64, b) * number_of_nodes + @intCast(u64, a);
            return @intCast(u64, a) * number_of_nodes + @intCast(u64, b);
        }

        fn removedCallback(self: *Self, key: u32, list: ?*std.AutoArrayHashMapUnmanaged(u32, void)) void {
            if (list) |l| {
                var iter = l.iterator();
                while (iter.next()) |partner_ptr| {
                    const partner = partner_ptr.key_ptr.*;
                    if (partner == key) continue;
                    const unique_key = Self.calculateUniqueKey(partner, key, self.number_of_nodes);
                    if (self.hit_map.getPtr(unique_key)) |pt| {
                        _ = self.similarity_cardinality[pt.* - 1].swapRemove(unique_key);
                        pt.* -= 1;
                        if (pt.* == 0) {
                            _ = self.hit_map.remove(unique_key);
                        } else {
                            _ = self.similarity_cardinality[pt.* - 1].put(self.allocator, unique_key, {}) catch @panic("Could not change similarity");
                        }
                        // Untrack this node now!
                        if ((pt.* + 1) == self.similarity_threshold) {
                            // Remove from similarity list
                            self.untrackNodePair(partner, key) catch @panic("Error untrack node pair!");
                        }
                    } else {
                        @panic("Should not happen!");
                    }
                }
            }
        }

        fn addedCallback(self: *Self, key: u32, list: ?*std.AutoArrayHashMapUnmanaged(u32, void)) void {
            if (list) |l| {
                var iter = l.iterator();
                while (iter.next()) |partner_ptr| {
                    const partner = partner_ptr.key_ptr.*;
                    if (partner == key) continue;
                    const unique_key = Self.calculateUniqueKey(partner, key, self.number_of_nodes);
                    if (self.hit_map.getPtr(unique_key)) |pt| {
                        _ = self.similarity_cardinality[pt.* - 1].swapRemove(unique_key);
                        pt.* += 1;
                        _ = self.similarity_cardinality[pt.* - 1].put(self.allocator, unique_key, {}) catch @panic("Could not change similarity");
                        // Track this move now
                        if (self.similarity_threshold == pt.* and self.tww_nb.count() < 100) {
                            self.trackNodePair(partner, key) catch @panic("Error while tracking node");
                        }
                    } else {
                        self.hit_map.put(self.allocator, unique_key, 1) catch @panic("Out of memory!");
                    }
                }
            }
        }

        fn trackNodePair(self: *Self, erased: u32, survivor: u32) !void {
            var unique_key = Self.calculateUniqueKey(erased, survivor, self.number_of_nodes);
            if (!self.tww_nb.contains(unique_key)) {
                var tww = self.graph.calculateMaxTwwScore(@intCast(T, erased), @intCast(T, survivor));
                try self.tww_nb.put(self.allocator, unique_key, tww.tww);
                if (self.update_list.getPtr(erased)) |ptr| {
                    try ptr.put(self.allocator, survivor, {});
                } else {
                    var set = std.AutoArrayHashMapUnmanaged(u32, void){};
                    try set.put(self.allocator, survivor, {});
                    try self.update_list.put(self.allocator, erased, set);
                }
                if (self.update_list.getPtr(survivor)) |ptr| {
                    try ptr.put(self.allocator, erased, {});
                } else {
                    var set = std.AutoArrayHashMapUnmanaged(u32, void){};
                    try set.put(self.allocator, erased, {});
                    try self.update_list.put(self.allocator, survivor, set);
                }
            }
        }

        fn untrackNodePair(self: *Self, erased: u32, survivor: u32) !void {
            if (self.update_list.getPtr(erased)) |pt| {
                _ = pt.swapRemove(survivor);
                if (self.update_list.getPtr(survivor)) |sr| {
                    _ = sr.swapRemove(erased);
                    const key = Self.calculateUniqueKey(erased, survivor, self.number_of_nodes);
                    _ = self.tww_nb.swapRemove(key);
                }
            }
        }

        fn untrackNode(self: *Self, node: u32) !void {
            if (self.update_list.getPtr(node)) |pt| {
                var iter = pt.iterator();
                while (iter.next()) |it| {
                    _ = self.tww_nb.swapRemove(Self.calculateUniqueKey(node, it.key_ptr.*, self.number_of_nodes));
                    var update_list = self.update_list.getPtr(it.key_ptr.*);
                    _ = update_list.?.swapRemove(node);
                }
                pt.deinit(self.allocator);
            }
            _ = self.update_list.remove(node);
        }

        pub fn removeNode(self: *Self, node: u32) !void {
            const split = @intCast(u32, self.bands.len >> 1);
            for (0..split) |index| {
                const result = try self.bands[index].removeItem(node);
                Self.removedCallback(self, node, result);
            }
            try self.untrackNode(node);
        }

        pub fn rehashNode(self: *Self, node: u32, graph: *graph_mod.Graph(T)) !void {
            for (0..self.bands.len) |index| {
                var iter = NodeIteratorSplitRedBlack{ .number_of_nodes = graph.number_of_nodes, .iter = graph.node_list[node].unorderedIterator() };

                try self.bands[index].rehashItem(*NodeIteratorSplitRedBlack, @TypeOf(self), node, &iter, self, Self.removedCallback, Self.addedCallback);
            }
        }

        pub fn changedEdge(self: *Self, node: u32, removed: u32, comptime removed_is_red: bool, added: ?u32, comptime added_is_red: bool, graph: *graph_mod.Graph(T)) !void {
            if (added != null and !added_is_red) {
                std.debug.panic("Added black edge!\n", .{});
            }

            for (0..self.bands.len) |index| {
                var iter = NodeIteratorSplitRedBlack{ .number_of_nodes = graph.number_of_nodes, .iter = graph.node_list[node].unorderedIterator() };
                try self.bands[index].updateItem(*NodeIteratorSplitRedBlack, @TypeOf(self), node, &iter, self, Self.removedCallback, Self.addedCallback, if (removed_is_red) removed + self.number_of_nodes else removed, if (added) |a| if (added_is_red) a + self.number_of_nodes else added else added);
            }
        }

        pub fn bootstrapNodes(self: *Self, nodes: []T, graph: *graph_mod.Graph(T), seed: u64) !void {
            MinHash.permutate(self.permutation, seed);
            self.graph = graph;
            for (0..self.bands.len) |j| {
                self.bands[j].clear();
            }
            self.hit_map.clearRetainingCapacity();
            self.similarity_threshold = 0;
            for (self.similarity_cardinality) |*n| {
                n.clearRetainingCapacity();
            }

            self.tww_nb.clearRetainingCapacity();
            var iter_v = self.update_list.valueIterator();
            while (iter_v.next()) |item| {
                item.deinit(self.allocator);
            }

            for (0..nodes.len) |index| {
                const i = nodes[index];

                if (graph.erased_nodes.get(i)) continue;
                for (0..self.bands.len) |j| {
                    var iter = NodeIteratorSplitRedBlack{ .number_of_nodes = graph.number_of_nodes, .iter = graph.node_list[i].unorderedIterator() };

                    if (try self.bands[j].addItem(*NodeIteratorSplitRedBlack, i, &iter)) |hits| {
                        var iterator = hits.iterator();
                        while (iterator.next()) |hit_pt| {
                            const hit = hit_pt.key_ptr.*;
                            if (hit != i) {
                                const key = Self.calculateUniqueKey(@intCast(u32, i), hit, graph.number_of_nodes);
                                if (self.hit_map.getPtr(key)) |it| {
                                    it.* += 1;
                                } else {
                                    try self.hit_map.put(graph.allocator, key, 1);
                                }
                            }
                        }
                    }
                }
            }

            var entries = self.hit_map.iterator();
            while (entries.next()) |value| {
                try self.similarity_cardinality[value.value_ptr.* - 1].put(self.allocator, value.key_ptr.*, {});
            }

            var threshold = self.bands.len - 1;
            var cum_cardinality: u32 = 0;
            while (threshold > 0) {
                cum_cardinality += @intCast(u32, self.similarity_cardinality[threshold].count());
                if (cum_cardinality > 100) break;
                threshold -= 1;
            }
            std.debug.print("Cardinality {}\n", .{cum_cardinality});

            self.similarity_threshold = @intCast(u32, threshold);

            for (threshold..self.bands.len) |index| {
                var iter = self.similarity_cardinality[index].iterator();
                while (iter.next()) |it| {
                    if (self.tww_nb.count() > 100) break;
                    const mv = MinHashEntry.keyIntoMove(T, it.key_ptr.*, graph.number_of_nodes);
                    try self.trackNodePair(mv.erased, mv.survivor);
                }
            }
        }

        pub fn updatedNode(self: *Self, node: T, graph: *graph_mod.Graph(T)) !void {
            if (self.update_list.getPtr(node)) |pt| {
                var iter = pt.iterator();
                while (iter.next()) |it| {
                    if (graph.erased_nodes.get(it.key_ptr.*)) continue;
                    try self.tww_nb.put(self.allocator, Self.calculateUniqueKey(node, it.key_ptr.*, graph.number_of_nodes), graph.calculateMaxTwwScore(node, @intCast(T, it.key_ptr.*)).tww_nb);
                }
            }
        }

        pub fn getDeltaMergeSimilarity(self: *Self, erased: u32, survivor: u32) !i32 {
            var max: i32 = 0;
            self.sim_nodes_after.clearRetainingCapacity();
            for (0..self.bands.len) |index| {
                if (try self.bands[index].getMergeSignatureHits(T, erased, survivor, self.graph)) |l| {
                    for (l.items) |it| {
                        if (self.sim_nodes_after.getPtr(it)) |x| {
                            x.* += 1;
                            max = std.math.max(max, @intCast(i32, x.*));
                        } else {
                            try self.sim_nodes_after.put(self.allocator, it, 1);
                        }
                    }
                }
            }
            return @intCast(i32, self.bands.len) - max;
        }

        pub fn getBestMove(self: *Self, graph: *graph_mod.Graph(T), current_tww: T) !?contraction.Contraction(T) {
            var min_cont: ?contraction.Contraction(T) = null;
            var best_tww = graph_mod.Graph(T).InducedTwinWidthPotential.default();

            // Fill up
            outer: while (self.tww_nb.count() < 100) {
                for (1..self.bands.len + 1) |id| {
                    const index = @intCast(u32, self.bands.len - id);
                    self.similarity_threshold = index;
                    var iter = self.similarity_cardinality[index].iterator();
                    while (iter.next()) |it| {
                        if (self.tww_nb.count() >= 100) break :outer;
                        const mv = MinHashEntry.keyIntoMove(T, it.key_ptr.*, graph.number_of_nodes);
                        try self.trackNodePair(mv.erased, mv.survivor);
                    }
                }
            }

            var iter = self.tww_nb.iterator();
            while (iter.next()) |it| {
                const mv = MinHashEntry.keyIntoMove(T, it.key_ptr.*, self.number_of_nodes);
								if(graph.erased_nodes.get(mv.erased) or graph.erased_nodes.get(mv.survivor)) {
									_ = self.tww_nb.swapRemove(it.key_ptr.*);
									continue;
								}
                var tww = graph.calculateInducedTwwPotential(mv.erased, mv.survivor, &best_tww, current_tww);

                if (tww.isLess(best_tww, current_tww)) {
                    min_cont = .{ .erased = mv.erased, .survivor = mv.survivor };
                    best_tww = tww;
                }
            }
            return min_cont;
        }

        pub fn init(allocator: std.mem.Allocator, num_bands: u32, number_of_nodes: u32) !Self {
            var bands = try allocator.alloc(MinHashBand(B), num_bands);
            var permutation = try MinHash.generatePermutation(allocator, number_of_nodes * num_bands * 2);
            var counter: u32 = 0;
            for (bands) |*band| {
                band.* = try MinHashBand(B).init(allocator, number_of_nodes, permutation[counter..(counter + 2 * number_of_nodes)]);
                counter += number_of_nodes * 2;
            }

            var hit_map = std.AutoHashMapUnmanaged(u64, u32){};

            var tww_nb = std.AutoArrayHashMapUnmanaged(u64, u32){};
            var update_list = std.AutoHashMapUnmanaged(u32, std.AutoArrayHashMapUnmanaged(u32, void)){};

            var similarity_cardinality = try allocator.alloc(std.AutoArrayHashMapUnmanaged(u64, void), num_bands);
            for (similarity_cardinality) |*card| {
                card.* = std.AutoArrayHashMapUnmanaged(u64, void){};
            }
            var similarity_threshold: u32 = 0;

            return Self{ .bands = bands, .permutation = permutation, .hit_map = hit_map, .allocator = allocator, .number_of_nodes = number_of_nodes, .graph = undefined, .tww_nb = tww_nb, .update_list = update_list, .similarity_cardinality = similarity_cardinality, .similarity_threshold = similarity_threshold };
        }
    };
}

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
