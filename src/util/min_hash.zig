const std = @import("std");
const comptime_util = @import("comptime_checks.zig");
const compressed_bitmap = @import("compressed_bitmap.zig");
const two_level_bitset = @import("../util/two_level_bitset.zig");
const node_mod = @import("../graph/node.zig");
const contraction = @import("../tww/contraction_sequence.zig");
const graph_mod = @import("../graph/graph.zig");
const dart_hash = @import("dart_hash.zig");

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

pub fn hashFnWeighted(ctx: MinHashBandContextWeighted, value: []u8) u64 {
    _ = ctx;
    return std.hash.Wyhash.hash(0, value);
}

pub fn eqFnWeighted(ctx: MinHashBandContextWeighted, rhs: []u8, lhs: []u8) bool {
    _ = ctx;
    return std.mem.eql(u8, rhs, lhs);
}
pub const MinHashBandContextWeighted = struct {
    pub const hash = hashFnWeighted;
    pub const eql = eqFnWeighted;
};
pub fn MinHashBandWeighted(comptime B: u32) type {
    return struct {
        const Self = @This();
        dart_hash: dart_hash.DartMinHash,
        hash_cache: []u8,
        collisions: std.HashMapUnmanaged([]u8, std.AutoArrayHashMapUnmanaged(u32, void), MinHashBandContextWeighted, 80),
        allocator: std.mem.Allocator,

        pub fn clear(self: *Self) void {
            var iter = self.collisions.iterator();
            while (iter.next()) |item| {
                item.value_ptr.deinit(self.allocator);
                self.allocator.free(item.key_ptr.*);
            }
            self.collisions.clearRetainingCapacity();
        }

        pub inline fn init(allocator: std.mem.Allocator, cache_size: u32, seed: u64) !Self {
            var rng = std.rand.DefaultPrng.init(seed);
            var hash = try dart_hash.DartMinHash.init(allocator, &rng, B);

            var collisions = std.HashMapUnmanaged([]u8, std.AutoArrayHashMapUnmanaged(u32, void), MinHashBandContextWeighted, 80){};

            var hash_cache = try allocator.alloc(u8, cache_size * B);

            return Self{ .dart_hash = hash, .collisions = collisions, .allocator = allocator, .hash_cache = hash_cache };
        }

        pub inline fn rehashItem(self: *Self, comptime IterType: type, comptime Context: type, key: u32, total_weight: f64, iterator: IterType, context: Context, callback_removed: fn (Context, u32, ?*std.AutoArrayHashMapUnmanaged(u32, void)) void, callback_added: fn (Context, u32, ?*std.AutoArrayHashMapUnmanaged(u32, void)) void) !void {
            const hash = try self.dart_hash.hash(IterType, iterator, total_weight);
            if (!std.mem.eql(u8, self.hash_cache[key * B .. key * B + B], hash)) {
                const first = try self.removeItem(key);
                callback_removed(context, key, first);

                const next = try self.addItem(IterType, key, iterator, total_weight);
                callback_added(context, key, next);
            }
        }

        inline fn addItemFromCache(self: *Self, key: u32) !?*std.AutoArrayHashMapUnmanaged(u32, void) {
            if (self.collisions.getPtr(self.hash_cache[key * B .. key * B + B])) |result| {
                try result.put(self.allocator, key, {});
                return result;
            } else {
                var copy = try self.allocator.alloc(u8, B);
                @memcpy(copy, self.hash_cache[key * B .. key * B + B]);

                var new_list = std.AutoArrayHashMapUnmanaged(u32, void){};
                try new_list.put(self.allocator, key, {});
                try self.collisions.put(self.allocator, copy, new_list);
            }
            return null;
        }

        // Returns the collisions if there were any
        pub inline fn addItem(self: *Self, comptime IterType: type, key: u32, iterator: IterType, total_weight: f64) !?*std.AutoArrayHashMapUnmanaged(u32, void) {
            const hash = try self.dart_hash.hash(IterType, iterator, total_weight);
            @memcpy(self.hash_cache[key * B .. key * B + B], hash);
            return self.addItemFromCache(key);
        }

        pub inline fn getSignature(self: *const Self, key: u32) []u8 {
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

pub fn MinHashSimilarityWeighted(comptime T: type, comptime B: u32) type {
    return struct {
        const Self = @This();

        bands: []MinHashBandWeighted(B),
        hit_map: std.AutoHashMapUnmanaged(u64, u32),

        tww_nb: std.AutoArrayHashMapUnmanaged(u64, u32),

        similarity_cardinality: []std.AutoArrayHashMapUnmanaged(u64, void),
        similarity_threshold: u32,

        number_of_nodes: u32,
        allocator: std.mem.Allocator,
        graph: *graph_mod.Graph(T),

        pub const CANIDATE_COUNT: u32 = 1000;

        pub const WeightedItem = struct {
            index: u32,
            weight: f64,
        };

        pub inline fn keyIntoMove(key: u64, number_of_nodes: u32) contraction.Contraction(T) {
            var first = key / number_of_nodes;
            var second = key % number_of_nodes;
            return contraction.Contraction(T){ .erased = @intCast(T, first), .survivor = @intCast(T, second) };
        }

        const NodeIteratorSplitRedBlack = struct {
            number_of_nodes: u32,
            graph: *graph_mod.Graph(T),
            iter: graph_mod.Graph(T).NodeType.UnorderedNodeEdgeIterator,
            final_node: ?u32 = null,
            final_node_sent: bool = false,
            pub inline fn reset(self: *NodeIteratorSplitRedBlack) void {
                self.iter.reset();
                self.final_node_sent = false;
            }

            pub fn next(self: *NodeIteratorSplitRedBlack) ?WeightedItem {
                if (self.iter.next()) |item| {
                    const weight = @intToFloat(f64, self.graph.node_list[item].red_edges.cardinality() + 1);
                    //if (self.iter.red) return .{.index=item + self.number_of_nodes, .weight = weight};
                    return .{ .index = item, .weight = weight };
                }
                if (self.final_node) |item| {
                    if (self.final_node_sent) return null;
                    self.final_node_sent = true;
                    if (item >= self.number_of_nodes) {
                        const weight = @intToFloat(f64, self.graph.node_list[item - self.number_of_nodes].red_edges.cardinality() + 1);
                        return .{ .index = item - self.number_of_nodes, .weight = weight };
                    } else {
                        const weight = @intToFloat(f64, self.graph.node_list[item].red_edges.cardinality() + 1);
                        return .{ .index = item, .weight = weight };
                    }
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
            if (a >= number_of_nodes or b >= number_of_nodes) @panic("Wrong key");
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
                        if (self.similarity_threshold == pt.* and self.tww_nb.count() < Self.CANIDATE_COUNT) {
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
                //var tww = self.graph.calculateMaxTwwScore(@intCast(T, erased), @intCast(T, survivor));
                try self.tww_nb.put(self.allocator, unique_key, 0);
            }
        }

        fn untrackNodePair(self: *Self, erased: u32, survivor: u32) !void {
            const key = Self.calculateUniqueKey(erased, survivor, self.number_of_nodes);
            _ = self.tww_nb.swapRemove(key);
        }

        pub fn removeNode(self: *Self, node: u32) !void {
            for (0..self.bands.len) |index| {
                const result = try self.bands[index].removeItem(node);
                Self.removedCallback(self, node, result);
            }
        }

        pub const ChangedEdge = struct {
            removed: u32,
            red: bool = false,
            added: ?u32 = null,
            added_red: bool = false,
        };

        pub fn changedEdge(self: *Self, node: u32, graph: *graph_mod.Graph(T), changed: ChangedEdge) !void {
            {
                var delta_weight: f32 = @intToFloat(f32, self.graph.node_list[changed.removed].red_edges.cardinality());
                graph.node_list[node].total_weight -= delta_weight;
                delta_weight = delta_weight / @floatCast(f32, graph.node_list[node].total_weight);

                graph.node_list[node].delta_potential_weighted_jaccard += delta_weight;
            }

            if (changed.added) |add| {
                var delta_weight: f32 = @intToFloat(f32, self.graph.node_list[add].red_edges.cardinality());
                graph.node_list[node].total_weight += delta_weight;
                delta_weight = delta_weight / @floatCast(f32, graph.node_list[node].total_weight);
                graph.node_list[node].delta_potential_weighted_jaccard += delta_weight;
            }
            if (graph.node_list[node].delta_potential_weighted_jaccard > 0.05) {
                return self.rehashNode(node, graph);
            }
        }

        pub fn rehashNode(self: *Self, node: u32, graph: *graph_mod.Graph(T)) !void {
            const split = self.bands.len / 3;
            var total_weight: f64 = 0.0;
            graph.node_list[node].delta_potential_weighted_jaccard = 0.0;

            var iterator = NodeIteratorSplitRedBlack{ .number_of_nodes = graph.number_of_nodes, .iter = graph.node_list[node].unorderedIterator(), .graph = graph };
            while (iterator.next()) |item| {
                total_weight += item.weight;
            }
            iterator.reset();
            graph.node_list[node].total_weight = @floatCast(f32, total_weight);
            for (0..split) |index| {
                var iter = NodeIteratorSplitRedBlack{ .number_of_nodes = graph.number_of_nodes, .iter = graph.node_list[node].unorderedIterator(), .graph = graph };

                try self.bands[index].rehashItem(*NodeIteratorSplitRedBlack, @TypeOf(self), node, total_weight, &iter, self, Self.removedCallback, Self.addedCallback);
            }
            total_weight += @intToFloat(f64, graph.node_list[node].red_edges.cardinality());
            for (split..split * 2) |index| {
                var iter = NodeIteratorSplitRedBlack{ .number_of_nodes = graph.number_of_nodes, .iter = graph.node_list[node].unorderedIterator(), .final_node = node, .graph = graph };

                try self.bands[index].rehashItem(*NodeIteratorSplitRedBlack, @TypeOf(self), node, total_weight, &iter, self, Self.removedCallback, Self.addedCallback);
            }
            for (split * 2..self.bands.len) |index| {
                var iter = NodeIteratorSplitRedBlack{ .number_of_nodes = graph.number_of_nodes, .iter = graph.node_list[node].unorderedIterator(), .final_node = node + self.number_of_nodes, .graph = graph };

                try self.bands[index].rehashItem(*NodeIteratorSplitRedBlack, @TypeOf(self), node, total_weight, &iter, self, Self.removedCallback, Self.addedCallback);
            }
        }

        pub fn bootstrapNodes(self: *Self, nodes: []T, graph: *graph_mod.Graph(T)) !void {
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

            const split = self.bands.len / 3;

            for (0..nodes.len) |index| {
                const i = nodes[index];

                if (graph.erased_nodes.get(i)) continue;
                for (0..split) |j| {
                    var iter = NodeIteratorSplitRedBlack{ .number_of_nodes = graph.number_of_nodes, .iter = graph.node_list[i].unorderedIterator(), .graph = graph };

                    if (try self.bands[j].addItem(*NodeIteratorSplitRedBlack, i, &iter, @intToFloat(f64, graph.node_list[i].cardinality()))) |hits| {
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
                for (split..split * 2) |j| {
                    var iter = NodeIteratorSplitRedBlack{ .number_of_nodes = graph.number_of_nodes, .iter = graph.node_list[i].unorderedIterator(), .final_node = i, .graph = graph };

                    if (try self.bands[j].addItem(*NodeIteratorSplitRedBlack, i, &iter, @intToFloat(f64, graph.node_list[i].cardinality()))) |hits| {
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
                for (split * 2..self.bands.len) |j| {
                    var iter = NodeIteratorSplitRedBlack{ .number_of_nodes = graph.number_of_nodes, .iter = graph.node_list[i].unorderedIterator(), .final_node = i + self.number_of_nodes, .graph = graph };

                    if (try self.bands[j].addItem(*NodeIteratorSplitRedBlack, i, &iter, @intToFloat(f64, graph.node_list[i].cardinality()))) |hits| {
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
                if (cum_cardinality > Self.CANIDATE_COUNT) break;
                threshold -= 1;
            }
            std.debug.print("Cardinality {}\n", .{cum_cardinality});

            self.similarity_threshold = @intCast(u32, threshold);

            for (threshold..self.bands.len) |index| {
                var iter = self.similarity_cardinality[index].iterator();
                while (iter.next()) |it| {
                    if (self.tww_nb.count() > Self.CANIDATE_COUNT) break;
                    const mv = Self.keyIntoMove(it.key_ptr.*, graph.number_of_nodes);
                    try self.trackNodePair(mv.erased, mv.survivor);
                }
            }
        }

        pub fn getBestMove(self: *Self, graph: *graph_mod.Graph(T), current_tww: T) !?contraction.Contraction(T) {
            var min_cont: ?contraction.Contraction(T) = null;

            var erased_list = try std.BoundedArray(u64, 1000).init(0);
            var iternb = self.tww_nb.iterator();
            while (iternb.next()) |it| {
                const mv = Self.keyIntoMove(it.key_ptr.*, graph.number_of_nodes);
                if (graph.erased_nodes.get(mv.erased) or graph.erased_nodes.get(mv.survivor)) {
                    try erased_list.append(it.key_ptr.*);
                }
            }
            for (erased_list.buffer[0..erased_list.len]) |it| {
                _ = self.tww_nb.swapRemove(it);
            }

            // Fill up
            outer: for (1..self.bands.len + 1) |id| {
                const index = @intCast(u32, self.bands.len - id);
                self.similarity_threshold = index;
                var iter = self.similarity_cardinality[index].iterator();
                while (iter.next()) |it| {
                    if (self.tww_nb.count() >= Self.CANIDATE_COUNT) break :outer;
                    const mv = Self.keyIntoMove(it.key_ptr.*, graph.number_of_nodes);
                    try self.trackNodePair(mv.erased, mv.survivor);
                }
            }

            if (self.similarity_threshold == 0 and self.tww_nb.count() == 0) {
                std.debug.print("No moves left\n", .{});
            }

            var best_tww = graph_mod.Graph(T).InducedTwinWidthPotential.default();
            var iter = self.tww_nb.iterator();
            while (iter.next()) |it| {
                const mv = Self.keyIntoMove(it.key_ptr.*, self.number_of_nodes);
                if (graph.erased_nodes.get(mv.erased) or graph.erased_nodes.get(mv.survivor)) {
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
            var bands = try allocator.alloc(MinHashBandWeighted(B), num_bands);
            if (num_bands % 3 != 0) @panic("Number of bands must be divisible by 3!");
            var seed: u64 = 129;
            for (bands) |*band| {
                band.* = try MinHashBandWeighted(B).init(allocator, number_of_nodes, seed);
                seed += 12362;
            }

            var hit_map = std.AutoHashMapUnmanaged(u64, u32){};

            var tww_nb = std.AutoArrayHashMapUnmanaged(u64, u32){};

            var similarity_cardinality = try allocator.alloc(std.AutoArrayHashMapUnmanaged(u64, void), num_bands);
            for (similarity_cardinality) |*card| {
                card.* = std.AutoArrayHashMapUnmanaged(u64, void){};
            }
            var similarity_threshold: u32 = 0;

            return Self{ .bands = bands, .hit_map = hit_map, .allocator = allocator, .number_of_nodes = number_of_nodes, .graph = undefined, .tww_nb = tww_nb, .similarity_cardinality = similarity_cardinality, .similarity_threshold = similarity_threshold };
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

        pub fn clear(self: *Self) void {
            var iter = self.collisions.iterator();
            while (iter.next()) |item| {
                item.value_ptr.deinit(self.allocator);
                self.allocator.free(item.key_ptr.*);
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

            return Self{ .hash_functions = hashes, .collisions = collisions, .allocator = allocator, .hash_cache = hash_cache };
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

pub fn MinHashSimilarity(comptime T: type, comptime B: u32) type {
    return struct {
        const Self = @This();

        bands: []MinHashBand(B),
        hit_map: std.AutoHashMapUnmanaged(u64, u32),

        tww_nb: std.AutoArrayHashMapUnmanaged(u64, u32),

        similarity_cardinality: []std.AutoArrayHashMapUnmanaged(u64, void),
        similarity_threshold: u32,

        number_of_nodes: u32,
        allocator: std.mem.Allocator,
        graph: *graph_mod.Graph(T),

        pub const CANIDATE_COUNT: u32 = 1000;

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

        fn removedCallback(self: *Self, key: u32, list: ?*std.AutoArrayHashMapUnmanaged(u32, void)) void {
            if (list) |l| {
                var iter = l.iterator();
                while (iter.next()) |partner_ptr| {
                    const partner = partner_ptr.key_ptr.*;
                    if (partner == key) continue;
                    const unique_key = Self.calculateUniqueKey(partner, key);
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
                    const unique_key = Self.calculateUniqueKey(partner, key);
                    if (self.hit_map.getPtr(unique_key)) |pt| {
                        _ = self.similarity_cardinality[pt.* - 1].swapRemove(unique_key);
                        pt.* += 1;
                        _ = self.similarity_cardinality[pt.* - 1].put(self.allocator, unique_key, {}) catch @panic("Could not change similarity");
                        // Track this move now
                        if (self.similarity_threshold == pt.* and self.tww_nb.count() < Self.CANIDATE_COUNT) {
                            self.trackNodePair(partner, key) catch @panic("Error while tracking node");
                        }
                    } else {
                        self.hit_map.put(self.allocator, unique_key, 1) catch @panic("Out of memory!");
                    }
                }
            }
        }

        fn trackNodePair(self: *Self, erased: u32, survivor: u32) !void {
            var unique_key = Self.calculateUniqueKey(erased, survivor);
            if (!self.tww_nb.contains(unique_key)) {
                //var tww = self.graph.calculateMaxTwwScore(@intCast(T, erased), @intCast(T, survivor));
                try self.tww_nb.put(self.allocator, unique_key, 0);
            }
        }

        fn untrackNodePair(self: *Self, erased: u32, survivor: u32) !void {
            const key = Self.calculateUniqueKey(erased, survivor);
            _ = self.tww_nb.swapRemove(key);
        }

        pub fn removeNode(self: *Self, node: u32) !void {
            for (0..self.bands.len) |index| {
                const result = try self.bands[index].removeItem(node);
                Self.removedCallback(self, node, result);
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
            if (changed.added != null and !changed.added_red) {
                std.debug.panic("Added black edge!\n", .{});
            }

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

        pub fn bootstrapNodes(self: *Self, nodes: []T, graph: *graph_mod.Graph(T), seed: u64) !void {
            var rng = std.rand.DefaultPrng.init(seed);
            var random = rng.random();

            self.graph = graph;
            for (self.bands) |*b| {
                b.clear();
                b.randomize_functions(&random);
            }

            self.hit_map.clearRetainingCapacity();
            self.similarity_threshold = 0;
            for (self.similarity_cardinality) |*n| {
                n.clearRetainingCapacity();
            }

            self.tww_nb.clearRetainingCapacity();

            const split = self.bands.len / 3;

            for (0..nodes.len) |index| {
                const i = nodes[index];

                if (graph.erased_nodes.get(i)) continue;
                for (0..split) |j| {
                    var iter = NodeIteratorSplitRedBlack{ .number_of_nodes = graph.number_of_nodes, .iter = graph.node_list[i].unorderedIterator() };

                    if (try self.bands[j].addItem(*NodeIteratorSplitRedBlack, i, &iter)) |hits| {
                        var iterator = hits.iterator();
                        while (iterator.next()) |hit_pt| {
                            const hit = hit_pt.key_ptr.*;
                            if (hit != i) {
                                const key = Self.calculateUniqueKey(@intCast(u32, i), hit);
                                if (self.hit_map.getPtr(key)) |it| {
                                    it.* += 1;
                                } else {
                                    try self.hit_map.put(graph.allocator, key, 1);
                                }
                            }
                        }
                    }
                }
                for (split..split * 2) |j| {
                    var iter = NodeIteratorSplitRedBlack{ .number_of_nodes = graph.number_of_nodes, .iter = graph.node_list[i].unorderedIterator(), .final_node = i };

                    if (try self.bands[j].addItem(*NodeIteratorSplitRedBlack, i, &iter)) |hits| {
                        var iterator = hits.iterator();
                        while (iterator.next()) |hit_pt| {
                            const hit = hit_pt.key_ptr.*;
                            if (hit != i) {
                                const key = Self.calculateUniqueKey(@intCast(u32, i), hit);
                                if (self.hit_map.getPtr(key)) |it| {
                                    it.* += 1;
                                } else {
                                    try self.hit_map.put(graph.allocator, key, 1);
                                }
                            }
                        }
                    }
                }
                for (split * 2..self.bands.len) |j| {
                    var iter = NodeIteratorSplitRedBlack{ .number_of_nodes = graph.number_of_nodes, .iter = graph.node_list[i].unorderedIterator(), .final_node = i + self.number_of_nodes };

                    if (try self.bands[j].addItem(*NodeIteratorSplitRedBlack, i, &iter)) |hits| {
                        var iterator = hits.iterator();
                        while (iterator.next()) |hit_pt| {
                            const hit = hit_pt.key_ptr.*;
                            if (hit != i) {
                                const key = Self.calculateUniqueKey(@intCast(u32, i), hit);
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
                if (cum_cardinality > Self.CANIDATE_COUNT) break;
                threshold -= 1;
            }

            self.similarity_threshold = @intCast(u32, threshold);

            for (threshold..self.bands.len) |index| {
                var iter = self.similarity_cardinality[index].iterator();
                while (iter.next()) |it| {
                    if (self.tww_nb.count() > Self.CANIDATE_COUNT) break;
                    const mv = Self.keyIntoMove(it.key_ptr.*);
                    try self.trackNodePair(mv.erased, mv.survivor);
                }
            }
        }

        pub fn getBestMove(self: *Self, graph: *graph_mod.Graph(T), current_tww: T) !?contraction.Contraction(T) {
            var min_cont: ?contraction.Contraction(T) = null;

            {
                var i: usize = 0;
                while (i < self.tww_nb.count()) {
                    const key = self.tww_nb.keys()[i];
                    const mv = Self.keyIntoMove(key);
                    if (graph.erased_nodes.get(mv.erased) or graph.erased_nodes.get(mv.survivor)) {
                        _ = self.tww_nb.swapRemove(key);
                    } else {
                        i += 1;
                    }
                }
            }

            // Fill up
            outer: for (1..self.bands.len + 1) |id| {
                const index = @intCast(u32, self.bands.len - id);
                self.similarity_threshold = index;
                var iter = self.similarity_cardinality[index].iterator();
                while (iter.next()) |it| {
                    if (self.tww_nb.count() >= Self.CANIDATE_COUNT) break :outer;
                    const mv = Self.keyIntoMove(it.key_ptr.*);
                    try self.trackNodePair(mv.erased, mv.survivor);
                }
            }

            if (self.similarity_threshold == 0 and self.tww_nb.count() == 0) {
                std.debug.print("No moves left\n", .{});
            }

            var best_tww = graph_mod.Graph(T).InducedTwinWidthPotential.default();
            var iter = self.tww_nb.iterator();
            while (iter.next()) |it| {
                const mv = Self.keyIntoMove(it.key_ptr.*);
                if (graph.erased_nodes.get(mv.erased) or graph.erased_nodes.get(mv.survivor)) {
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

            if (num_bands % 3 != 0) @panic("Number of bands must be divisible by 3!");
            var counter: u32 = 0;
            for (bands) |*band| {
                band.* = try MinHashBand(B).init(allocator, number_of_nodes);
                counter += number_of_nodes * 2;
            }

            var hit_map = std.AutoHashMapUnmanaged(u64, u32){};

            var tww_nb = std.AutoArrayHashMapUnmanaged(u64, u32){};

            var similarity_cardinality = try allocator.alloc(std.AutoArrayHashMapUnmanaged(u64, void), num_bands);
            for (similarity_cardinality) |*card| {
                card.* = std.AutoArrayHashMapUnmanaged(u64, void){};
            }
            var similarity_threshold: u32 = 0;

            return Self{ .bands = bands, .hit_map = hit_map, .allocator = allocator, .number_of_nodes = number_of_nodes, .graph = undefined, .tww_nb = tww_nb, .similarity_cardinality = similarity_cardinality, .similarity_threshold = similarity_threshold };
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

    var band = try MinHashBandWeighted(4).init(gpa.allocator(), 100, perm, 0);
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
