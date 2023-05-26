const std = @import("std");
const bitset = @import("../util/two_level_bitset.zig");
const comptime_util = @import("../util/comptime_checks.zig");

pub fn ParametrizedUnsortedArrayList(comptime T: type) type {
    if (!comptime_util.checkIfIsCompatibleInteger(T)) {
        @compileError("Type must either be u8,16 or u32!");
    }

    return struct {
        const Self = @This();
        edges: std.ArrayListUnmanaged(T),
        pub inline fn initCapacity(allocator: std.mem.Allocator, size: u32) !Self {
            std.debug.assert(size > 0);
            var list = try std.ArrayListUnmanaged(T).initCapacity(allocator, size);

            return Self{ .edges = list };
        }

        pub inline fn init() Self {
            var list = std.ArrayListUnmanaged(T){};

            return Self{ .edges = list };
        }

        pub const ParametrizedUnsortedArrayListIterator = struct {
            index: u32,
            list: *const Self,
            pub inline fn next(self: *ParametrizedUnsortedArrayListIterator) ?T {
                if (self.index >= self.list.edges.items.len) return null;
                const item = self.list.edges.items[self.index];
                self.index += 1;
                return item;
            }
        };

        pub inline fn intoSorted(self: *Self) ParametrizedSortedArrayList(T) {
            std.sort.sort(T, self.edges.items, {}, comptime std.sort.asc(T));
            return ParametrizedSortedArrayList(T){ .edges = self.edges };
        }

        pub inline fn iterator(self: *const Self) ParametrizedUnsortedArrayListIterator {
            return ParametrizedUnsortedArrayListIterator{ .index = 0, .list = self };
        }

        pub inline fn add(self: *Self, allocator: std.mem.Allocator, item: T) !void {
            try self.edges.append(allocator, item);
        }

        pub inline fn cardinality(self: *const Self) u32 {
            return @intCast(u32, self.edges.items.len);
        }

        pub fn deinit(self: *Self, allocator: std.mem.Allocator) void {
            self.edges.deinit(allocator);
        }

        pub inline fn nodePosition(self: *const Self, item: T) ?u32 {
            var i: u32 = 0;
            while (i < self.edges.items.len) : (i += 1) {
                if (self.edges.items[i] == item) {
                    return i;
                }
            }
            return null;
        }

        pub inline fn contains(self: *const Self, item: T) bool {
            return self.nodePosition(item) != null;
        }

        pub inline fn remove(self: *Self, item: T) bool {
            const result = self.nodePosition(item);
            if (result) |position| {
                _ = self.edges.swapRemove(position);
                return true;
            }
            return false;
        }
    };
}

pub fn ParametrizedSortedArrayList(comptime T: type) type {
    if (!comptime_util.checkIfIsCompatibleInteger(T)) {
        @compileError("Type must either be u8,16 or u32!");
    }

    return struct {
        const Self = @This();
        edges: std.ArrayListUnmanaged(T),
        pub inline fn initCapacity(allocator: std.mem.Allocator, size: u32) !Self {
            var list = try std.ArrayListUnmanaged(T).initCapacity(allocator, size);

            return Self{ .edges = list };
        }

        pub inline fn init() Self {
            var list = std.ArrayListUnmanaged(T){};
            return Self{ .edges = list };
        }

        pub const ParametrizedSortedArrayListIterator = struct {
            index: u32,
            list: *const Self,
						pub inline fn reset(self: *ParametrizedSortedArrayListIterator) void {
							self.index = 0;
						}

            pub inline fn next(self: *ParametrizedSortedArrayListIterator) ?T {
                if (self.index >= self.list.edges.items.len) return null;
                const item = self.list.edges.items[self.index];
                self.index += 1;
                return item;
            }
        };

        pub inline fn iterator(self: *const Self) ParametrizedSortedArrayListIterator {
            return ParametrizedSortedArrayListIterator{ .index = 0, .list = self };
        }

        pub inline fn xorIterator(self: *const Self, other: *const Self) ParametrizedSortedArrayListXorIterator {
            var iterator_first = self.iterator();
            var iterator_second = other.iterator();
            const item_first = iterator_first.next();
            const item_second = iterator_second.next();
            return ParametrizedSortedArrayListXorIterator{ .iterator_first = iterator_first, .iterator_second = iterator_second, .item_first = item_first, .item_second = item_second, .first = false };
        }

        pub const ParametrizedSortedArrayListXorIterator = struct {
            iterator_first: ParametrizedSortedArrayListIterator,
            iterator_second: ParametrizedSortedArrayListIterator,
            item_first: ?T,
            item_second: ?T,
            first: bool,
            pub inline fn next(self: *ParametrizedSortedArrayListXorIterator) ?T {
                while (self.item_first != null and self.item_second != null) {
                    if (self.item_first.? < self.item_second.?) {
                        const return_item_first = self.item_first;
                        self.item_first = self.iterator_first.next();
                        self.first = true;
                        return return_item_first;
                    } else if (self.item_first.? > self.item_second.?) {
                        const return_item_second = self.item_second;
                        self.item_second = self.iterator_second.next();
                        self.first = false;
                        return return_item_second;
                    }

                    self.item_first = self.iterator_first.next();
                    self.item_second = self.iterator_second.next();
                }

                if (self.item_second != null) {
                    const ret = self.item_second;
                    self.item_second = self.iterator_second.next();
                    self.first = false;
                    return ret;
                } else if (self.item_first != null) {
                    const ret = self.item_first;
                    self.item_first = self.iterator_first.next();
                    self.first = true;
                    return ret;
                }
                return null;
            }
        };

        pub inline fn add(self: *Self, allocator: std.mem.Allocator, item: T) !bool {
            const result = self.nodePosition(item);
            if (result.@"0" == false) {
                try self.edges.insert(allocator, result.@"1", item);
            }
            return result.@"0";
        }

        pub inline fn cardinality(self: *const Self) u32 {
            return @intCast(u32, self.edges.items.len);
        }

        pub fn deinit(self: *Self, allocator: std.mem.Allocator) void {
            self.edges.deinit(allocator);
        }

        pub inline fn contains(self: *const Self, item: T) bool {
            return self.nodePosition(item).@"0";
        }

        pub inline fn nodePosition(self: *const Self, item: T) struct { bool, u32 } {
            // Larger than 4 cache lines
            //const threshold = comptime (64/@sizeOf(T))*4;
            //TODO: Change this to a sensible number again!
            if (self.edges.items.len >= comptime (64 / @sizeOf(T)) * 8) {
                return self.binarySearch(item);
            } else {
                var i: u32 = 0;
                while (i < self.edges.items.len) : (i += 1) {
                    if (self.edges.items[i] == item) {
                        return .{ true, i };
                    } else if (self.edges.items[i] > item) {
                        return .{ false, i };
                    }
                }
                return .{ false, @intCast(u32, self.edges.items.len) };
            }
        }

        pub fn binarySearch(self: *const Self, target: T) struct { bool, u32 } {
            var left: usize = 0;
            var right = self.edges.items.len;

            while (left < right) {
                const mid = left + (right - left) / 2; // Avoid overflow.
                if (self.edges.items[mid] == target) {
                    return .{ true, @intCast(u32, mid) };
                } else if (self.edges.items[mid] < target) {
                    left = mid + 1;
                } else {
                    right = mid;
                }
            }
            return .{ false, @intCast(u32, left) };
        }

        pub inline fn removeMask(self: *Self, set: *bitset.FastBitSet) void {
            var i: u32 = 0;
            var write_ptr: u32 = 0;
            while (i < self.edges.items.len) : (i += 1) {
                const item = self.edges.items[i];
                self.edges.items[write_ptr] = item;
                if (!set.get(self.edges.items[i])) {
                    write_ptr += 1;
                }
            }
            self.edges.shrinkRetainingCapacity(write_ptr);
        }

        pub inline fn remove(self: *Self, item: T) bool {
            const result = self.nodePosition(item);
            if (result.@"0" == true) {
                _ = self.edges.orderedRemove(result.@"1");
            }
            return result.@"0";
        }
    };
}

test "ParametrizedSortedArrayList: Add" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer std.debug.assert(!gpa.deinit());

    var list = ParametrizedSortedArrayList(u32).init();
    defer list.deinit(gpa.allocator());

    _ = try list.add(gpa.allocator(), 1);

    try std.testing.expectEqual(list.contains(1), true);
    try std.testing.expectEqual(list.contains(2), false);

    // Should not add anything since the item is in the set already
    try std.testing.expectEqual(try list.add(gpa.allocator(), 1), true);

    var iter = list.iterator();
    var counter: u32 = 0;
    while (iter.next()) |item| {
        try std.testing.expectEqual(item, 1);
        counter += 1;
    }
    try std.testing.expectEqual(counter, 1);
    try std.testing.expectEqual(counter, list.cardinality());
}

test "ParametrizedSortedArrayList: Remove" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer std.debug.assert(!gpa.deinit());

    var list = ParametrizedSortedArrayList(u32).init();
    defer list.deinit(gpa.allocator());

    for (0..500) |item| {
        _ = try list.add(gpa.allocator(), @intCast(u32, item) * 2);
    }

    try std.testing.expectEqual(list.contains(1), false);
    try std.testing.expectEqual(list.contains(2), true);
    try std.testing.expectEqual(list.binarySearch(300).@"0", true);
    try std.testing.expectEqual(list.binarySearch(301).@"0", false);
    try std.testing.expectEqual(list.binarySearch(301).@"1", 151);

    var iter = list.iterator();
    var counter: u32 = 0;
    while (iter.next()) |_| {
        counter += 1;
    }
    try std.testing.expectEqual(counter, 500);
    try std.testing.expectEqual(counter, list.cardinality());

    var result = list.remove(100);
    try std.testing.expectEqual(result, true);
    result = list.remove(100);
    try std.testing.expectEqual(result, false);
}

test "ParametrizedSortedArrayList: Remove Mask" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer std.debug.assert(!gpa.deinit());

    var list = ParametrizedSortedArrayList(u32).init();
    defer list.deinit(gpa.allocator());

    var mask = try bitset.FastBitSet.initEmpty(5000, gpa.allocator());
    defer mask.deinit(gpa.allocator());

    for (0..500) |item| {
        _ = try list.add(gpa.allocator(), @intCast(u32, item) * 2);
        mask.set(@intCast(u32, item) * 4);
    }

    list.removeMask(&mask);

    var iter = list.iterator();
    var counter: u32 = 0;

    var last_item: ?u32 = null;
    while (iter.next()) |item| {
        counter += 1;
        if (last_item) |it| {
            // Expect ordered set
            try std.testing.expect(it < item);
        } else {
            last_item = item;
        }
    }
    try std.testing.expectEqual(counter, 250);
    try std.testing.expectEqual(counter, list.cardinality());
}
