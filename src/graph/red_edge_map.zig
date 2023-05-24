const std = @import("std");
const comptime_util = @import("../util/comptime_checks.zig");

pub fn RedEdgeMap(comptime T: type) type {
    if (!comptime_util.checkIfIsCompatibleInteger(T)) {
        @compileError("Only supports u8,u16 and u32");
    }

    if (T == u8 or T == u16) {
        return struct {
            const Self = @This();
            set: std.bit_set.DynamicBitSetUnmanaged,
            pub fn init(allocator: std.mem.Allocator) !Self {
                var bitset = try std.bit_set.DynamicBitSetUnmanaged.initEmpty(allocator, std.math.maxInt(T) * std.math.maxInt(T));
                return Self{ .set = bitset };
            }

            pub fn add(self: *Self, allocator: std.mem.Allocator, item: u64) !void {
                _ = allocator;
                self.set.set(item);
            }

            pub fn remove(self: *Self, item: u64) void {
                self.set.unset(item);
            }

            pub fn contains(self: *Self, item: u64) bool {
                return self.set.isSet(item);
            }

            pub fn deinit(self: *Self, allocator: std.mem.Allocator) void {
                self.set.deinit(allocator);
            }
        };
    } else if (T == u32) {
        return struct {
            const Self = @This();
            set: std.AutoHashMapUnmanaged(u64, u8),
            pub fn init(allocator: std.mem.Allocator) !Self {
                _ = allocator;
                var hashset = std.AutoHashMapUnmanaged(u64, u8){};
                return Self{ .set = hashset };
            }

            pub fn add(self: *Self, allocator: std.mem.Allocator, item: u64) !void {
                try self.set.put(allocator, item, 1);
            }

            pub fn remove(self: *Self, item: u64) void {
                _ = self.set.remove(item);
            }

            pub fn contains(self: *Self, item: u64) bool {
                return self.set.contains(item);
            }

            pub fn deinit(self: *Self, allocator: std.mem.Allocator) void {
                self.set.deinit(allocator);
            }
        };
    } else {
        @compileError("Unsupported type!");
    }
}
