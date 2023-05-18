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
    shift: u32,
    pub inline fn generatePermutation(allocator: std.mem.Allocator, length: u32, seed: u64) ![]u32 {
        var memory = try allocator.alloc(u32, length);
        for (0..memory.len) |i| {
            memory[i] = i;
        }

        var gen = std.rand.DefaultPrng.init(seed);
        fisher_yates_shuffle(u32, memory, &gen);
        return memory;
    }

    pub inline fn init(permutation: []u32, shift: u32) Self {
        return Self{
            .permutation = permutation,
            .shift = shift,
        };
    }

    pub inline fn getHash(self: *const Self, item: u32) u32 {
        return self.permutation[(item + self.shift) % self.permutation.len];
    }

    pub inline fn hash(self: *const Self, comptime InputIteratorType: type, iter: InputIteratorType) u32 {
        var min: u32 = std.math.maxInt(u32);
        while (iter.next()) |item| {
            min = std.math.min(min, self.getHash(item));
        }
        return min;
    }
};

pub fn MinHashBand(comptime B: u32) type {
    return struct {
        const Self = @This();
        hash_functions: []MinHash,
				permutation: []u32,
				hash_cache: []u32,
				collision_set: std.AutoHashMapUnmanaged([B]u32,void),
        collisions: std.AutoHashMapUnmanaged([B]u32, std.ArrayListUnmanaged(u32)),
				allocator: std.mem.Allocator,

        pub inline fn init(allocator: std.mem.Allocator, permutation: []u32, shift: u32) !Self {
            var hashes = try allocator.alloc(MinHash, B);

            for (hashes) |hash| {
                hash.* = MinHash.init(permutation,shift);
            }

            var collisions_set = std.AutoArrayHashMapUnmanaged(u128,void){};
            var collisions = std.AutoArrayHashMapUnmanaged(u128,std.ArrayListUnmanaged(u32)){};


            return Self {
                .hashes = hashes,
								.collisions_set = collisions_set,
                .collisions = collisions,
								.allocator = allocator
            };
        }

				pub inline fn rehashItem(self: *Self, comptime IterType: type, key: u32, hash: u128, iterator: IterType) !struct{?*std.ArrayListUnmanaged(u32),?*std.ArrayListUnmanaged} {
				    _ = iterator;
					const first = try self.removeItem(key,hash);
					_ = first;

				}

				// Returns the collisions if there were any
				pub inline fn addItem(self: *Self, comptime IterType: type, key: u32, iterator: IterType) !?*std.ArrayListUnmanaged(u32) {
						var bounded_arr = try std.BoundedArray(u32,B).init(0);
						{
							var i:u32 = 0;
							inline while(i<B) : (i+=1) {
								try bounded_arr.append(std.math.maxInt(u32));
							}
						}

						while(iterator.next()) |item| {
							for(0..B) |i| {
								bounded_arr.buffer[i] = std.math.min(bounded_arr.buffer[i],self.hash_functions[i].getHash(item));
							}
						}

						{
							var i:u32 = 0;
							var hash:u128 = 0;
							inline while(i<B) : (i+=1) {
								hash |= (@intCast(u128,bounded_arr.buffer[i])<<(i*32))	;
							}

							if(try self.collision_set.fetchPut(self.allocator, hash,{}) != null) {
								var result = self.collisions.getPtr(hash).?;
								try result.append(self.allocator, key);
								return result;
							}
						}
						return null;
				}

				pub inline fn removeItem(self: *Self, item: u32, key: u128) !?*std.ArrayListUnmanaged(u32) {
					if(self.collisions.getPtr(key)) |ptr| {
						for(0..ptr.items.len) |index| {
							if(ptr.items[index] == item) {
								ptr.swapRemove(index);
								break;
							}
						}
						if(ptr.items.len != 1) return ptr;
						
						ptr.deinit(self.allocator);
						_ = self.collisions.remove(key);
					}
					else {
						_ = self.collision_set.remove(key);
					}
					return null;
				}

        pub fn deinit(self: *Self, allocator: std.mem.Allocator) void {
            self.sim_nodes.deinit(allocator);
            var iter = self.collisions.iterator();
            while (iter.next()) |item| {
                item.value_ptr.deinit(allocator);
            }
            self.collisions.deinit(allocator);
            allocator.free(self.node_hashes);
            for (self.hashes) |*hash| {
                hash.deinit(allocator);
            }
        }
    };
}

pub const MinHashEntry = struct {
    const Self = @This();
    key: u64,
    score: u32,
    pub fn from(key: u64, score: u32, max_score: u32, max_cardinality: u32) MinHashEntry {
        return MinHashEntry{ .key = key, .score = Self.calculateScore(score, max_score, max_cardinality) };
    }

    pub inline fn calculateScore(score: u32, max_score: u32, max_cardinality: u32) u32 {
        return (max_score - score) * max_cardinality;
    }

    pub inline fn intoMove(self: Self, comptime T: type, number_of_nodes: u32) contraction.Contraction(T) {
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

pub fn MinHashSimiliarity(comptime T: type) type {
    return struct {
        const Self = @This();
        const TransferedEdge = struct { T, bool, bool };
        const RemovedEdge = struct { T, bool };

        bands: []MinHashBand(4),
        score_table: std.AutoArrayHashMapUnmanaged(u64, u32),
        pq_moves: std.PriorityQueue(MinHashEntry, void, MinHashEntry.compare),
        score_counter: []u32,
        cutoff: u32,
        allocator: std.mem.Allocator,

        changed_color_survivor: std.ArrayListUnmanaged(T),

        // From erased to survivor
        transfered_edges: std.ArrayListUnmanaged(TransferedEdge),

        // When erased and survivor have both the edge
        removed_edges: std.ArrayListUnmanaged(RemovedEdge),

        fetched_moves: std.ArrayListUnmanaged(MinHashEntry),

        pub inline fn init(num_bands: u32, band_width: u32, graph: *graph_mod.Graph(T)) !Self {
            var bands = try graph.allocator.alloc(MinHashBand(4), num_bands);
            var score_table = std.AutoArrayHashMapUnmanaged(u64, u32){};
            var moves = std.PriorityQueue(MinHashEntry, void, MinHashEntry.compare).init(graph.allocator, {});

            var score_counter = try graph.allocator.alloc(u32, num_bands + 1);
            @memset(score_counter, 0);

            var transfered_edges = try std.ArrayListUnmanaged(TransferedEdge).initCapacity(graph.allocator, graph.number_of_nodes);
            var removed_edges = try std.ArrayListUnmanaged(RemovedEdge).initCapacity(graph.allocator, graph.number_of_nodes);
            var changed_color_survivor = try std.ArrayListUnmanaged(T).initCapacity(graph.allocator, graph.number_of_nodes);

            var seed: u64 = 19;

            var mhs = Self{ .bands = bands, .score_table = score_table, .pq_moves = moves, .cutoff = num_bands >> 2, .allocator = graph.allocator, .transfered_edges = transfered_edges, .removed_edges = removed_edges, .score_counter = score_counter, .fetched_moves = try std.ArrayListUnmanaged(MinHashEntry).initCapacity(graph.allocator, graph.number_of_nodes), .changed_color_survivor = changed_color_survivor };

            for (mhs.bands) |*band| {
                band.* = try MinHashBand(T).init(band_width, graph, seed);
                seed += 12395321;
                var iter = band.sim_nodes.iterator();
                while (iter.next()) |item| {
                    if (mhs.score_table.getPtr(item.key_ptr.*)) |score| {
                        score_counter[score.*] -= 1;
                        score.* += 1;
                        score_counter[score.*] += 1;
                    } else {
                        try mhs.score_table.put(graph.allocator, item.key_ptr.*, 1);
                        score_counter[1] += 1;
                    }
                }
            }

            try mhs.adjustCutoffValue(graph);
            return mhs;
        }

        pub fn adjustCutoffValue(self: *Self, graph: *graph_mod.Graph(T)) !void {
            var total_counter: u32 = 0;
            var cutoff: u32 = 0;
            for (0..self.score_counter.len) |i| {
                total_counter += self.score_counter[self.score_counter.len - (i + 1)];
                cutoff = @intCast(u32, self.score_counter.len - (i + 1));
                if (total_counter >= 10000) break;
            }

            // Clear pq
            self.cutoff = cutoff;
            self.pq_moves.len = 0;

            var iter = self.score_table.iterator();
            while (iter.next()) |item| {
                if (item.value_ptr.* >= self.cutoff) {
                    const move = MinHashEntry.keyIntoMove(T, item.key_ptr.*, graph.number_of_nodes);
                    try self.pq_moves.add(MinHashEntry.from(item.key_ptr.*, item.value_ptr.*, @intCast(u32, self.bands.len), std.math.max(graph.node_list[move.erased].cardinality(), graph.node_list[move.survivor].cardinality())));
                }
            }
        }

        pub inline fn decreaseScore(self: *Self, key: u64) !void {
            if (self.score_table.getPtr(key)) |k| {
                self.score_counter[k.*] -= 1;
                k.* -= 1;
                self.score_counter[k.*] += 1;
                if (k.* == 0) {
                    _ = self.score_table.swapRemove(key);
                }
            }
        }

        pub inline fn increaseScore(self: *Self, key: u64, graph: *graph_mod.Graph(T)) !void {
            if (self.score_table.getPtr(key)) |k| {
                self.score_counter[k.*] -= 1;
                k.* += 1;
                self.score_counter[k.*] += 1;
                if (k.* >= self.cutoff) {
                    const move = MinHashEntry.keyIntoMove(T, key, graph.number_of_nodes);
                    try self.pq_moves.add(MinHashEntry.from(key, k.*, @intCast(u32, self.bands.len), std.math.max(graph.node_list[move.erased].cardinality(), graph.node_list[move.survivor].cardinality())));
                }
            }
        }

        pub inline fn addTransferedEdge(self: *Self, node: T, from_color: bool, to_color: bool) void {
            self.transfered_edges.append(self.allocator, .{ node, from_color, to_color }) catch unreachable;
        }

        pub inline fn addChangeColorSurvivorEdge(self: *Self, node: T) void {
            self.changed_color_survivor.append(self.allocator, node) catch unreachable;
        }

        pub inline fn addRemovedEdge(self: *Self, node: T, color: bool) void {
            self.removed_edges.append(self.allocator, .{ node, color }) catch unreachable;
        }

        pub inline fn batchUpdateRehashNodes(self: *Self, erased: T, survivor: T, graph: *graph_mod.Graph(T)) !void {
            for (self.bands) |*band| {
                for (self.removed_edges.items) |node_index| {
                    try band.rehashNodeRemovedEdge(node_index.@"0", node_index.@"1", graph, self, erased);
                }
                for (self.transfered_edges.items) |node_index| {
                    try band.rehashNodeTransferedEdge(node_index.@"0", node_index.@"1", node_index.@"2", graph, self, erased, survivor);
                }
                for (self.changed_color_survivor.items) |node_index| {
                    try band.rehashNodeChangedColor(node_index, graph, self, survivor);
                }
                try band.removeNode(erased, self);
            }
            self.removed_edges.clearRetainingCapacity();
            self.transfered_edges.clearRetainingCapacity();
            self.changed_color_survivor.clearRetainingCapacity();
        }

        pub inline fn rehashNode(self: *Self, node_id: T, graph: *graph_mod.Graph(T)) !void {
            for (self.bands) |*band| {
                try band.rehashNode(node_id, graph, self);
            }
        }

        pub inline fn removeNode(self: *Self, node_id: T) !void {
            for (self.bands) |*band| {
                try band.removeNode(node_id, self);
            }
        }

        pub inline fn reinsertMove(self: *Self, move: MinHashEntry) !void {
            try self.pq_moves.add(move);
        }

        pub inline fn reinsertFetchedMoves(self: *Self) !void {
            try self.pq_moves.addSlice(self.fetched_moves.items);
            self.fetched_moves.clearRetainingCapacity();
        }

        pub inline fn getNextMove(self: *Self, graph: *graph_mod.Graph(T)) !?MinHashEntry {
            while (true) {
                while (self.pq_moves.removeOrNull()) |move| {
                    var mv = move.intoMove(T, graph.number_of_nodes);

                    if (graph.erased_nodes.get(mv.erased) or graph.erased_nodes.get(mv.survivor)) continue;

                    if (self.score_table.get(move.key)) |score| {
                        if (score < self.cutoff) continue;

                        const calculated_score = MinHashEntry.calculateScore(score, @intCast(u32, self.bands.len), std.math.max(graph.node_list[mv.erased].cardinality(), graph.node_list[mv.survivor].cardinality()));

                        if (calculated_score > move.score) {
                            try self.pq_moves.add(MinHashEntry.from(move.key, score, @intCast(u32, self.bands.len), std.math.max(graph.node_list[mv.erased].cardinality(), graph.node_list[mv.survivor].cardinality())));
                            continue;
                        }
                        return move;
                    }
                }

                // Adjust cutoff value
                try self.adjustCutoffValue(graph);
                if (self.pq_moves.len == 0) return null;
            }
        }

        pub inline fn fetchNextMoves(self: *Self, k: u32, graph: *graph_mod.Graph(T)) !void {
            std.debug.assert(self.fetched_moves.items.len == 0);
            while (try self.getNextMove(graph)) |item| {
                try self.fetched_moves.append(self.allocator, item);
                if (self.fetched_moves.items.len >= k) return;
            }
        }
    };
}

test "MinHash: Check is permutation" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer std.debug.assert(gpa.deinit() == .ok);
    var perm = MinHash(u16).generatePermutation(gpa.allocator(), 200, 19);
    var min_hash = try MinHash(u16).init(perm, 0);
    var min_hash2 = try MinHash(u16).init(perm, 1);

    var hash_map = std.AutoHashMap(u32, void).init(gpa.allocator());
    defer hash_map.deinit();

    for (min_hash.permutation) |item| {
        if (hash_map.contains(item)) try std.testing.expectEqual(true, false);
        try hash_map.put(item, {});
    }

    try std.testing.expectEqual(min_hash.permutation.len, 200);

    try std.testing.expectEqual(min_hash.getHash(100), min_hash2.getHash(101));
    try std.testing.expectEqual(min_hash.getValueOfNode(99, false), min_hash.permutation[199]);
}

test "MinHash: Check hashing of complex types" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer std.debug.assert(gpa.deinit() == .ok);

    var min_hash = try MinHash(u16).init(gpa.allocator(), 100, 19);
    defer min_hash.deinit(gpa.allocator());

    var hash_map = std.AutoHashMap(u32, void).init(gpa.allocator());
    defer hash_map.deinit();

    var initial = try std.BoundedArray(u16, 10).init(0);
    try initial.append(10);
    try initial.append(20);

    const TestIterator = struct {
        const Self = @This();
        base: *std.BoundedArray(u16, 10),
        index: u32,
        pub fn next(self: *Self) ?u16 {
            if (self.base.len <= self.index) return null;
            const item = self.base.buffer[self.index];
            self.index += 1;
            return item;
        }
    };

    var iter_black = TestIterator{ .base = &initial, .index = 0 };
    var iter_red = TestIterator{ .base = &initial, .index = 0 };

    var hash = min_hash.hash(@TypeOf(&iter_black), &iter_black, &iter_red);

    var min: u32 = std.math.maxInt(u16);
    for (initial.buffer[0..initial.len]) |item| {
        min = std.math.min(min, min_hash.permutation[item]);
    }
    for (initial.buffer[0..initial.len]) |item| {
        min = std.math.min(min, min_hash.permutation[item + 100]);
    }

    try std.testing.expectEqual(min, hash);
}
