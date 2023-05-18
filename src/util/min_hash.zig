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
            memory[i] = @intCast(u32,i);
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

pub fn MinHashSimilarity(comptime T: type, comptime B:u32) type {
    _ = B;
    _ = T;
	return struct {

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
    var min_hash = MinHash.init(perm,0);

    var hash_map = std.AutoHashMap(u32, void).init(gpa.allocator());
    defer hash_map.deinit();

    var initial = try std.BoundedArray(u32, 10).init(0);
    try initial.append(10);
    try initial.append(20);

    const TestIterator = struct {
        const Self = @This();
        base: *std.BoundedArray(u32, 10),
        index: u32,
				offset: u32,
        pub fn next(self: *Self) ?u32 {
            if (self.base.len <= self.index) return null;
            const item = self.base.buffer[self.index];
            self.index += 1;
            return item+self.offset;
        }
    };

    var iter_black = TestIterator{ .base = &initial, .index = 0, .offset = 0 };
    var iter_red = TestIterator{ .base = &initial, .index = 0, .offset = 100 };

    var hash = min_hash.hash(@TypeOf(&iter_black), &iter_black);
    hash = std.math.min(min_hash.hash(@TypeOf(&iter_red), &iter_red),hash);

    var min: u32 = std.math.maxInt(u32);
    for (initial.buffer[0..initial.len]) |item| {
        min = std.math.min(min, min_hash.permutation[item]);
    }
    for (initial.buffer[0..initial.len]) |item| {
        min = std.math.min(min, min_hash.permutation[item + 100]);
    }

    try std.testing.expectEqual(min, hash);
}
