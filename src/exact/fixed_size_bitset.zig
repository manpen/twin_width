const std = @import("std");

pub fn FixedSizeSet(comptime size: usize) type {
    const storeCardByDefault = true;
    _ = storeCardByDefault;

    if (size <= 8) {
        return FixedSizeSetImpl(u8, size, false);
    }

    if (size <= 16) {
        return FixedSizeSetImpl(u16, size, false);
    }

    if (size <= 32) {
        return FixedSizeSetImpl(u32, size, false);
    }

    if (size <= 64) {
        return FixedSizeSetImpl(u64, size, false);
    }

    return FixedSizeSetImpl(u64, size, false or (size > 128));
}

fn log2Type(comptime T: type) type {
    return switch (T) {
        u8 => u3,
        u16 => u4,
        u32 => u5,
        u64 => u6,
        else => unreachable,
    };
}

fn FixedSizeSetImpl(comptime T: type, comptime size: usize, comptime store_cardinality: bool) type {
    const N = comptime (size + @bitSizeOf(T) - 1) / @bitSizeOf(T);
    const Zero = @as(T, 0);
    const One = @as(T, 1);

    const maskLastWord = if (size % @bitSizeOf(T) == 0)
        ~Zero
    else
        (One << @intCast(log2Type(T), size % @bitSizeOf(T))) - One;

    const Container = comptime if (store_cardinality) struct { data: [N]T, card: u32 } else struct { data: [N]T };

    return struct {
        const Self = @This();

        container: Container,

        pub fn new() Self {
            return Self{ .container = comptime if (store_cardinality) Container{ .data = [_]T{0} ** N, .card = 0 } else Container{ .data = [_]T{0} ** N } };
        }

        pub fn newAllSet() Self {
            var me = Self{ .container = comptime if (store_cardinality) Container{ .data = [_]T{~Zero} ** N, .card = size } else Container{ .data = [_]T{~Zero} ** N } };
            me.container.data[N - 1] = maskLastWord;
            return me;
        }

        pub fn new_with_bits_set(bits: []const u32) Self {
            var result = Self.new();
            for (bits) |i| {
                _ = result.setBit(i);
            }
            return result;
        }

        pub fn is_equal(self: *const Self, other: *const Self) bool {
            if (store_cardinality) {
                if (self.container.card != other.container.card) {
                    return false;
                }
            }

            for (self.container.data, other.container.data) |a, b| {
                if (a != b) {
                    return false;
                }
            }

            return true;
        }

        pub fn is_subset_of(self: *const Self, other: *const Self) bool {
            if (store_cardinality) {
                if (self.container.card > other.container.card) {
                    return false;
                }
            }

            for (self.container.data, other.container.data) |a, b| {
                if (a & ~b != 0) {
                    return false;
                }
            }

            return true;
        }

        pub inline fn cardinality(self: *const Self) u32 {
            if (store_cardinality) {
                return self.container.card;
            }

            return self.compute_cardinality();
        }

        pub inline fn areAllUnset(self: *const Self) bool {
            if (store_cardinality) {
                return 0 == self.container.card;
            }

            for (self.container.data) |x| {
                if (x != 0) {
                    return false;
                }
            }

            return true;
        }

        fn compute_cardinality(self: *const Self) u32 {
            var sum: u32 = 0;

            for (self.container.data) |x| {
                sum += @popCount(x);
            }

            return sum;
        }

        fn update_cardinality_if_necessary(self: *Self) void {
            if (store_cardinality) {
                self.container.card = self.compute_cardinality();
            }
        }

        pub inline fn isSet(self: *const Self, idx: u32) bool {
            const word_idx = idx / @bitSizeOf(T);

            const L = log2Type(T);
            const bit_idx = @intCast(L, idx % @bitSizeOf(T));

            const mask: T = One << bit_idx;

            return 0 != (self.container.data[word_idx] & mask);
        }

        pub inline fn assignBit(self: *Self, idx: u32, value: bool) bool {
            const word_idx = idx / @bitSizeOf(T);

            const L = log2Type(T);
            const bit_idx = @intCast(L, idx % @bitSizeOf(T));

            const mask: T = One << bit_idx;

            var word = self.container.data[word_idx];

            const previous = (word & mask) != 0;
            word &= ~mask;
            word |= @as(T, @boolToInt(value)) << bit_idx;
            self.container.data[word_idx] = word;

            if (store_cardinality) {
                self.container.card -= @boolToInt(previous);
                self.container.card += @boolToInt(value);
            }

            return previous;
        }

        pub inline fn moveBitAndUnset(self: *Self, src: u32, dst: u32) void {
            const src_word: *T = &self.container.data[src / @bitSizeOf(T)];
            var dst_word: *T = &self.container.data[dst / @bitSizeOf(T)];

            const L = log2Type(T);
            const src_bit_idx = @intCast(L, src % @bitSizeOf(T));
            const dst_bit_idx = @intCast(L, dst % @bitSizeOf(T));

            var tmp: T = src_word.*;

            if (src_bit_idx < dst_bit_idx) {
                tmp <<= (dst_bit_idx - src_bit_idx);
            } else {
                tmp >>= (src_bit_idx - dst_bit_idx);
            }

            // unset original bit
            const src_bit_mask = One << src_bit_idx;
            src_word.* &= ~src_bit_mask;

            // shift bit
            const dst_bit_mask = One << dst_bit_idx;
            dst_word.* = (dst_word.* & ~dst_bit_mask) | (tmp & dst_bit_mask);
        }

        pub fn unsetAll(self: *Self) void {
            for (&self.container.data) |*s| {
                s.* = 0;
            }
            if (store_cardinality) {
                self.container.card = 0;
            }
        }

        pub inline fn setBit(self: *Self, idx: u32) bool {
            return self.assignBit(idx, true);
        }

        pub inline fn unsetBit(self: *Self, idx: u32) bool {
            return self.assignBit(idx, false);
        }

        pub fn assignOr(self: *Self, other: *const Self) void {
            for (&self.container.data, other.container.data) |*s, o| {
                s.* |= o;
            }
            self.update_cardinality_if_necessary();
        }

        pub fn assignXor(self: *Self, other: *const Self) void {
            for (&self.container.data, other.container.data) |*s, o| {
                s.* ^= o;
            }
            self.update_cardinality_if_necessary();
        }

        pub fn assignAnd(self: *Self, other: *const Self) void {
            for (&self.container.data, other.container.data) |*s, o| {
                s.* &= o;
            }
            self.update_cardinality_if_necessary();
        }

        pub fn assignSub(self: *Self, other: *const Self) void {
            for (&self.container.data, other.container.data) |*s, o| {
                s.* &= ~o;
            }
            self.update_cardinality_if_necessary();
        }

        pub fn copyWithOr(self: *const Self, other: *const Self) Self {
            var result = Self{ .container = undefined };
            for (&result.container.data, self.container.data, other.container.data) |*r, s, o| {
                r.* = s | o;
            }
            result.update_cardinality_if_necessary();
            return result;
        }

        pub fn copyWithXor(self: *const Self, other: *const Self) Self {
            var result = Self{ .container = undefined };
            for (&result.container.data, self.container.data, other.container.data) |*r, s, o| {
                r.* = s ^ o;
            }
            result.update_cardinality_if_necessary();
            return result;
        }

        pub fn copyWithAnd(self: *const Self, other: *const Self) Self {
            var result = Self{ .container = undefined };
            for (&result.container.data, self.container.data, other.container.data) |*r, s, o| {
                r.* = s & o;
            }
            result.update_cardinality_if_necessary();
            return result;
        }

        pub fn copyWithSub(self: *const Self, other: *const Self) Self {
            var result = Self{ .container = undefined };
            for (&result.container.data, self.container.data, other.container.data) |*r, s, o| {
                r.* = s & ~o;
            }
            result.update_cardinality_if_necessary();
            return result;
        }

        pub fn iter_set(self: *const Self) FixedSizeSetIterator(T, true) {
            return FixedSizeSetIterator(T, true).new(self.container.data[0..], size);
        }

        pub fn iter_unset(self: *const Self) FixedSizeSetIterator(T, false) {
            return FixedSizeSetIterator(T, false).new(self.container.data[0..], size);
        }

        pub fn debugPrint(self: *const Self) void {
            for (self.container.data) |x| {
                std.debug.print("{b:0>32} ", .{@bitReverse(x)});
            }
            std.debug.print("\n", .{});
        }
    };
}

pub fn FixedSizeSetIterator(comptime T: type, comptime IterOnes: bool) type {
    return struct {
        const Self = @This();

        slice: []const T,

        current_word: T,
        word_index: u32,
        number_of_bits: u32,

        pub fn new(slice: []const T, number_of_bits: u32) Self {
            return Self{
                .slice = slice,
                .current_word = @as(T, 0),
                .word_index = 0,
                .number_of_bits = number_of_bits,
            };
        }

        pub inline fn next(self: *Self) ?u32 {
            while (self.current_word == 0) {
                if (self.word_index >= self.slice.len) {
                    return null;
                }

                self.current_word = self.slice[self.word_index];
                if (!IterOnes) {
                    self.current_word = ~self.current_word;
                }

                self.word_index += 1;
            }

            const bit_idx = @intCast(log2Type(T), @ctz(self.current_word));
            self.current_word &= ~(@intCast(T, 1) << bit_idx);

            const result = (self.word_index - 1) * @bitSizeOf(T) + bit_idx;
            if (!IterOnes and result >= self.number_of_bits) {
                return null;
            }

            return result;
        }
    };
}

const assert = std.debug.assert;

test "Size of FixedSizeSet" {
    assert(@sizeOf(FixedSizeSet(1)) == 1);
    assert(@sizeOf(FixedSizeSet(8)) == 1);

    assert(@sizeOf(FixedSizeSet(9)) == 2);
    assert(@sizeOf(FixedSizeSet(16)) == 2);

    assert(@sizeOf(FixedSizeSet(17)) == 4);
    assert(@sizeOf(FixedSizeSet(32)) == 4);

    assert(@sizeOf(FixedSizeSet(33)) == 8);
    assert(@sizeOf(FixedSizeSet(64)) == 8);
}

test "Size of FixedSizeSetImpl" {
    assert(@sizeOf(FixedSizeSetImpl(u32, 1, false)) == 4);
    assert(@sizeOf(FixedSizeSetImpl(u64, 100, false)) == 16);
    assert(@sizeOf(FixedSizeSetImpl(u64, 130, false)) == 24);
    assert(@sizeOf(FixedSizeSetImpl(u64, 1, true)) > 8); // cardinality needs storage
}

test "NewAllSet" {
    assert(37 == FixedSizeSetImpl(u32, 37, false).newAllSet().cardinality());
}

test "Bit assignment" {
    var bs = FixedSizeSetImpl(u64, 100, false).new();
    assert(bs.cardinality() == 0);

    assert(bs.setBit(1) == false);
    assert(bs.cardinality() == 1);

    assert(bs.setBit(1) == true);
    assert(bs.cardinality() == 1);

    assert(bs.unsetBit(1) == true);
    assert(bs.cardinality() == 0);

    assert(bs.unsetBit(1) == false);
    assert(bs.cardinality() == 0);
}

test "Bit assignment with card" {
    var bs = FixedSizeSetImpl(u64, 100, true).new();
    assert(bs.cardinality() == 0);

    assert(bs.setBit(1) == false);
    assert(bs.cardinality() == 1);

    assert(bs.setBit(1) == true);
    assert(bs.cardinality() == 1);

    assert(bs.unsetBit(1) == true);
    assert(bs.cardinality() == 0);

    assert(bs.unsetBit(1) == false);
    assert(bs.cardinality() == 0);
}

test "bitwise operators" {
    const a = [_]u32{ 1, 2, 5, 7 };
    const b = [_]u32{ 0, 3, 5, 6 };

    const res_or = [_]u32{ 0, 1, 2, 3, 5, 6, 7 };
    const res_xor = [_]u32{ 0, 1, 2, 3, 6, 7 };
    const res_and = [_]u32{5};
    const res_sub = [_]u32{ 1, 2, 7 };

    const Set = FixedSizeSet(16);

    // Or
    {
        var set_a = Set.new_with_bits_set(&a);
        const set_b = Set.new_with_bits_set(b[0..]);
        const ref = Set.new_with_bits_set(res_or[0..]);

        const c = set_a.copyWithOr(&set_b);
        assert(c.is_equal(&ref));

        assert(!set_a.is_equal(&ref));
        set_a.assignOr(&set_b);
        assert(set_a.is_equal(&ref));
    }

    // Xor
    {
        var set_a = Set.new_with_bits_set(&a);
        const set_b = Set.new_with_bits_set(b[0..]);
        const ref = Set.new_with_bits_set(res_xor[0..]);

        const c = set_a.copyWithXor(&set_b);
        assert(c.is_equal(&ref));

        assert(!set_a.is_equal(&ref));
        set_a.assignXor(&set_b);
        assert(set_a.is_equal(&ref));
    }

    // And
    {
        var set_a = Set.new_with_bits_set(&a);
        const set_b = Set.new_with_bits_set(b[0..]);
        const ref = Set.new_with_bits_set(res_and[0..]);

        const c = set_a.copyWithAnd(&set_b);
        assert(c.is_equal(&ref));

        assert(!set_a.is_equal(&ref));
        set_a.assignAnd(&set_b);
        assert(set_a.is_equal(&ref));
    }

    // Sub
    {
        var set_a = Set.new_with_bits_set(&a);
        const set_b = Set.new_with_bits_set(b[0..]);
        const ref = Set.new_with_bits_set(res_sub[0..]);

        const c = set_a.copyWithSub(&set_b);
        assert(c.is_equal(&ref));

        assert(!set_a.is_equal(&ref));
        set_a.assignSub(&set_b);
        assert(set_a.is_equal(&ref));
    }
}

test "subset" {
    const bits = [_]u32{ 1, 2, 5, 7 };
    var sub = FixedSizeSet(32).new_with_bits_set(&bits);
    var sup = sub;

    assert(sub.is_equal(&sup));
    assert(sub.is_subset_of(&sup));

    _ = sup.setBit(3);
    assert(!sub.is_equal(&sup));
    assert(sub.is_subset_of(&sup));

    _ = sub.setBit(0);
    assert(!sub.is_subset_of(&sup));
}

test "MoveBitsLeft" {
    const bits = [_]u32{ 1, 2, 5, 7, 8 };
    var sub = FixedSizeSet(32).new_with_bits_set(&bits);
    var sup = sub;
    _ = sup;

    sub.moveBitAndUnset(7, 6); // 1 -> 0
    assert(!sub.isSet(7));
    assert(sub.isSet(6));

    std.debug.print("\n>", .{});
    sub.debugPrint();
    sub.moveBitAndUnset(7, 6); // 0 -> 1
    std.debug.print("<", .{});

    sub.debugPrint();
    assert(!sub.isSet(7));
    assert(!sub.isSet(6));

    sub.moveBitAndUnset(7, 6); // 0 -> 0
    assert(!sub.isSet(7));
    assert(!sub.isSet(6));

    sub.moveBitAndUnset(5, 2); // 1 -> 1
    assert(!sub.isSet(5));
    assert(sub.isSet(2));
}

test "MoveBitsRight" {
    const bits = [_]u32{ 1, 2, 5, 7 };
    var sub = FixedSizeSet(32).new_with_bits_set(&bits);
    var sup = sub;
    _ = sup;

    sub.moveBitAndUnset(2, 3); // 1 -> 0
    assert(!sub.isSet(2));
    assert(sub.isSet(3));

    sub.moveBitAndUnset(2, 3); // 0 -> 1
    assert(!sub.isSet(2));
    assert(!sub.isSet(3));

    sub.moveBitAndUnset(2, 3); // 0 -> 0
    assert(!sub.isSet(2));
    assert(!sub.isSet(3));

    sub.moveBitAndUnset(1, 5); // 1 -> 1
    assert(!sub.isSet(1));
    assert(sub.isSet(5));
}
