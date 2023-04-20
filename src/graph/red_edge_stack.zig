const std = @import("std");
const comptime_util = @import("../util/comptime_checks.zig");

pub fn NewRedEdgeType(comptime T: type) type {
    return enum(T) {
				// Converted own black edge to red edge due to last merge
        black_to_red_own,
				// Converted erased black edge to red edge due to last merge
				black_to_red_other,
				// Deleted a black edge from the surviving node since the erased node had the black edge too
        black_to_deleted,
				// Deleted a red edge from the surviving node since the erased node had the red edge too
				red_to_deleted,
    };
}

pub fn NewRedEdge(comptime T: type) type {
    comptime if (!comptime_util.checkIfIsCompatibleInteger(T)) {
        @compileError("Type must either be u8,16 or u32");
    };

    // This structure only stores turned red edges meaning red edges which were created due to the merge with the erased node together with the erased node and this information we can piece together from who came which edge at which time
    // new red edges
    // what edge turned (black -> red)
    // override (black -> deleted)
    return struct {
				const Self = @This();
        target: T,
        edge_type: NewRedEdgeType(T),
				
				pub inline fn blackToRedOwn(target: T) Self {
					return Self {
						.target = target,
						.edge_type = .black_to_red_own,
					};
				}

				pub inline fn blackToRedOther(target: T) Self {
					return Self {
						.target = target,
						.edge_type = .black_to_red_other,
					};
				}

				pub inline fn blackToDeleted(target: T) Self {
					return Self {
						.target = target,
						.edge_type = .black_to_deleted,
					};
				}

				pub inline fn redToDeleted(target: T) Self {
					return Self {
						.target = target,
						.edge_type = .red_to_deleted,
					};
				}
    };
}

pub const RedEdgeStackError = error{ StackNotSealed, StackNoLevelLeft, StackOverflow };

pub fn RedEdgeStack(comptime T: type) type {
    return struct {
        const Self = @This();
        edge_stack: std.ArrayListUnmanaged(NewRedEdge(T)),
        size: std.ArrayListUnmanaged(u32),

        pub fn init(allocator: std.mem.Allocator) !Self {
            var size = std.ArrayListUnmanaged(u32){};
            try size.append(allocator, 0);

            return Self { .edge_stack = std.ArrayListUnmanaged(NewRedEdge(T)){}, .size = size};
        }

        pub fn initCapacity(allocator: std.mem.Allocator, capacity: u32) !Self {
            std.debug.assert(capacity > 0);
            var size = std.ArrayListUnmanaged(u32){};
            try size.append(allocator, 0);

            return Self { .edge_stack = try std.ArrayListUnmanaged(NewRedEdge(T)).initCapacity(allocator, capacity), .size = size };
        }

        pub fn deinit(self: *Self, allocator: std.mem.Allocator) void {
            self.edge_stack.deinit(allocator);
            self.size.deinit(allocator);
        }

        pub const RedEdgeStackIterator = struct {
            index: u32,
            stack: *const Self,
            pub inline fn next(self: *RedEdgeStackIterator) ?NewRedEdge(T) {
                if (self.index >= self.stack.edge_stack.items.len) {
                    return null;
                }

                const item = self.stack.edge_stack.items[self.index];
                self.index += 1;
                return item;
            }
        };

        pub fn iterateLastLevel(self: *const Self) RedEdgeStackError!RedEdgeStackIterator {
            if (self.size.items[self.size.items.len - 1] != self.edge_stack.items.len) {
                return RedEdgeStackError.StackNotSealed;
            } else if (self.size.items.len < 2) {
                return RedEdgeStackError.StackNoLevelLeft;
            }
            return RedEdgeStackIterator{ .index = self.size.items[self.size.items.len - 2], .stack = self };
        }

        pub fn revertLastContraction(self: *Self) RedEdgeStackError!void {
            if (self.size.items.len >= 2) {
                // Reset to start of this level
                _ = self.size.pop();
                const resize = self.size.items[self.size.items.len - 1];
                self.edge_stack.shrinkRetainingCapacity(resize);
            } else {
                return RedEdgeStackError.StackNoLevelLeft;
            }
        }

        //TODO: Add exception here!
        pub inline fn addEdge(self: *Self, allocator: std.mem.Allocator, edge: NewRedEdge(T)) !void {
            try self.edge_stack.append(allocator, edge);
        }

        pub inline fn sealLevel(self: *Self, allocator: std.mem.Allocator) !void {
            try self.size.append(allocator, @intCast(u32,self.edge_stack.items.len));
        }
    };
}

test "RedEdgeStack: New red edges" {
    {
        const edge = NewRedEdge(u16).blackToRedOwn(100);
        try std.testing.expectEqual(edge.target, 100);
        try std.testing.expectEqual(edge.edge_type, .black_to_red_own);
    }

    {
        const edge = NewRedEdge(u16).blackToDeleted(102);
        try std.testing.expectEqual(edge.target, 102);
        try std.testing.expectEqual(edge.edge_type, .black_to_deleted);
    }
}

test "RedEdgeStack: Red edge stack add and iterate" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer {
        std.debug.assert(!gpa.deinit());
    }
    var stack = try RedEdgeStack(u16).init(gpa.allocator());
    defer stack.deinit(gpa.allocator());

    stack.addEdge(gpa.allocator(),NewRedEdge(u16).blackToRedOwn(1)) catch unreachable;
    stack.addEdge(gpa.allocator(),NewRedEdge(u16).blackToDeleted(2)) catch unreachable;
    try stack.sealLevel(gpa.allocator());
    stack.addEdge(gpa.allocator(),NewRedEdge(u16).blackToRedOwn(3)) catch unreachable;
    stack.addEdge(gpa.allocator(),NewRedEdge(u16).blackToDeleted(4)) catch unreachable;
    try stack.sealLevel(gpa.allocator());

    var iter = stack.iterateLastLevel() catch unreachable;
    var next = iter.next().?;
    var next_2 = iter.next().?;

    try std.testing.expectEqual(next.target, 3);
    try std.testing.expectEqual(next.edge_type, .black_to_red_own);

    try std.testing.expectEqual(next_2.target, 4);
    try std.testing.expectEqual(next_2.edge_type, .black_to_deleted);

    try stack.revertLastContraction();

    var iter_level_2 = stack.iterateLastLevel() catch unreachable;
    var next_level2 = iter_level_2.next().?;
    var next_level2_2 = iter_level_2.next().?;

    try std.testing.expectEqual(next_level2.target, 1);
    try std.testing.expectEqual(next_level2.edge_type, .black_to_red_own);

    try std.testing.expectEqual(next_level2_2.target, 2);
    try std.testing.expectEqual(next_level2_2.edge_type, .black_to_deleted);

    try stack.revertLastContraction();

    // Nothing left to iterate
    try std.testing.expectError(error.StackNoLevelLeft, stack.iterateLastLevel());
}

test "RedEdgeStack: Red edge stack add and empty stack" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer {
        std.debug.assert(!gpa.deinit());
    }
    var stack = try RedEdgeStack(u16).init(gpa.allocator());
    defer stack.deinit(gpa.allocator());
    stack.addEdge(gpa.allocator(),NewRedEdge(u16).blackToDeleted(1)) catch unreachable;

    try stack.sealLevel(gpa.allocator());
    try stack.sealLevel(gpa.allocator());

    var iterator = try stack.iterateLastLevel();
    try std.testing.expectEqual(iterator.next(), null);

    try stack.revertLastContraction();
    var iterator_now = try stack.iterateLastLevel();
    const item = iterator_now.next().?;
    try std.testing.expectEqual(item.target, 1);
    try std.testing.expectEqual(item.edge_type, .black_to_deleted);
    try std.testing.expectEqual(iterator_now.next(), null);
}

test "RedEdgeStack: Red edge stack revert and add again" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer {
        std.debug.assert(!gpa.deinit());
    }
    var stack = try RedEdgeStack(u16).init(gpa.allocator());
    defer stack.deinit(gpa.allocator());
    stack.addEdge(gpa.allocator(),NewRedEdge(u16).blackToDeleted(1)) catch unreachable;

    try stack.sealLevel(gpa.allocator());
    stack.addEdge(gpa.allocator(),NewRedEdge(u16).blackToDeleted(2)) catch unreachable;
    try stack.sealLevel(gpa.allocator());

    {
        var iterator = try stack.iterateLastLevel();
        const item_now = iterator.next().?;
        try std.testing.expectEqual(item_now.target, 2);
        try std.testing.expectEqual(item_now.edge_type, .black_to_deleted);
        try std.testing.expectEqual(iterator.next(), null);
    }

    try stack.revertLastContraction();
    stack.addEdge(gpa.allocator(),NewRedEdge(u16).blackToRedOwn(3)) catch unreachable;
    try stack.sealLevel(gpa.allocator());

    {
        var iterator = try stack.iterateLastLevel();
        const item_now = iterator.next().?;
        try std.testing.expectEqual(item_now.target, 3);
        try std.testing.expectEqual(item_now.edge_type, .black_to_red_own);

        //Next item must be null
        try std.testing.expectEqual(iterator.next(), null);
    }
    try stack.revertLastContraction();
    {
        var iterator = try stack.iterateLastLevel();
        const item_now = iterator.next().?;
        try std.testing.expectEqual(item_now.target, 1);
        try std.testing.expectEqual(item_now.edge_type, .black_to_deleted);
        try std.testing.expectEqual(iterator.next(), null);
    }
    try stack.revertLastContraction();
    try std.testing.expectError(error.StackNoLevelLeft, stack.iterateLastLevel());
    try std.testing.expectError(error.StackNoLevelLeft, stack.revertLastContraction());
}

test "RedEdgeStack: Red edge stack check safety guards" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer {
        std.debug.assert(!gpa.deinit());
    }
    var stack = try RedEdgeStack(u16).init(gpa.allocator());
    defer stack.deinit(gpa.allocator());
    stack.addEdge(gpa.allocator(),NewRedEdge(u16).blackToRedOwn(1)) catch unreachable;
    try std.testing.expectError(error.StackNotSealed, stack.iterateLastLevel());
}

test "RedEdgeStack: Check capacity" {
    // Space for 3 items
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer {
        std.debug.assert(!gpa.deinit());
    }
    var stack = try RedEdgeStack(u16).initCapacity(gpa.allocator(), 100);
    defer stack.deinit(gpa.allocator());

    try std.testing.expectEqual(stack.edge_stack.capacity, 100);
}
