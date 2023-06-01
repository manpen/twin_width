const std = @import("std");
const Allocator = std.mem.Allocator;
const assert = std.debug.assert;
const Order = std.math.Order;
const testing = std.testing;
const expect = testing.expect;
const expectEqual = testing.expectEqual;
const expectError = testing.expectError;

pub fn UpdateablePriorityQueueNode(comptime T: type, comptime P: type) type {
    return struct {
        const Self = @This();
        key: T,
        priority: P,

        pub fn init(key: T, priority: P) Self {
            return Self{
                .key = key,
                .priority = priority,
            };
        }
    };
}

/// Updateable Priority queue with up-to 64 bit keys.
pub fn UpdateablePriorityQueue(comptime T: type, comptime P: type, comptime Context: type, comptime compareFn: fn (context: Context, a: P, b: P) Order) type {
    return struct {
        const Self = @This();
        const KEY_NOT_FOUND = @as(T, std.math.maxInt(T));
        id_to_heap_pos: std.AutoHashMapUnmanaged(u64, u64),
        heap: []UpdateablePriorityQueueNode(T, P),
        len: usize,
        allocator: Allocator,
        context: Context,

        /// Initialize and return a priority queue.
        pub fn init(allocator: Allocator, context: Context) Self {
            return Self{
                .id_to_heap_pos = std.AutoHashMapUnmanaged(u64, u64){},
                .heap = &[_]UpdateablePriorityQueueNode(T, P){},
                .len = 0,
                .allocator = allocator,
                .context = context,
            };
        }

        pub fn empty(self: Self) bool {
            return (self.len == 0);
        }

        pub fn size(self: Self) usize {
            return self.len;
        }

        pub fn peek(self: *Self) ?*UpdateablePriorityQueueNode(T, P) {
            if (self.empty()) {
                return null;
            } else {
                return &self.heap[0];
            }
        }

        pub fn pop(self: *Self) !UpdateablePriorityQueueNode(T, P) {
            return self.removeAtHeapIdx(0);
        }

        pub fn removeOrNull(self: *Self) ?UpdateablePriorityQueueNode(T, P) {
            if (self.empty()) {
                return null;
            } else {
                return self.removeAtHeapIdx(0);
            }
        }

        pub fn removeKey(self: *Self, key: T) ?UpdateablePriorityQueueNode(T, P) {
            if (self.id_to_heap_pos.get(key)) |pos| {
                return self.removeAtHeapIdx(pos);
            }
            return null;
        }

        fn removeAtHeapIdx(self: *Self, index: T) UpdateablePriorityQueueNode(T, P) {
            const last = self.heap[self.len - 1];
            const item = self.heap[index];
            _ = self.id_to_heap_pos.remove(item.key);

            if (self.id_to_heap_pos.getPtr(last.key)) |ptr| {
                ptr.* = 0;
            }

            self.heap[index] = last;
            self.len -= 1;
            if (index == 0) {
                siftDown(self, index);
            } else {
                const parent_index = ((index - 1) >> 1);
                const parent = self.heap[parent_index];
                if (compareFn(self.context, last.priority, parent.priority) == .gt) {
                    self.siftDown(index);
                } else {
                    self.siftUp(index);
                }
            }

            return item;
        }

        pub fn updateOrInsert(self: *Self, key: T, priority: P) !bool {
            if (!self.id_to_heap_pos.contains(key)) { // new key
                try self.ensureUnusedCapacity(1);
                self.addUnchecked(key, priority);
            } else { // old key, just update
                try self.update(key, priority);
                return false;
            }
            return true;
        }

        fn capacity(self: Self) usize {
            return self.heap.len;
        }

        fn addUnchecked(self: *Self, key: T, priority: P) void {
            const node = UpdateablePriorityQueueNode(T, P).init(key, priority);
            self.heap[self.len] = node;
            const idx = @intCast(u32, self.len);
            self.id_to_heap_pos.put(self.allocator, key, idx) catch @panic("Out of memory");
            siftUp(self, idx);
            self.len += 1;
        }

        fn update(self: *Self, key: T, priority: P) !void {
            const idx = self.id_to_heap_pos.get(key).?;
            const old_priority = self.heap[idx].priority;
            switch (compareFn(self.context, priority, old_priority)) {
                .lt => siftUp(self, idx),
                .gt => siftDown(self, idx),
                .eq => {}, // Nothing to do as the items have equal priority
            }
        }

        fn siftUp(self: *Self, start_index: u64) void {
            var child_index = start_index;
            while (child_index > 0) {
                var parent_index = ((child_index - 1) >> 1);
                const child = self.heap[child_index];
                const parent = self.heap[parent_index];

                const child_key = child.key;
                const parent_key = parent.key;

                if (compareFn(self.context, child.priority, parent.priority) != .lt) break;

                self.heap[parent_index] = child;
                self.heap[child_index] = parent;

                self.id_to_heap_pos.getPtr(child_key).?.* = parent_index;
                self.id_to_heap_pos.getPtr(parent_key).?.* = child_index;

                child_index = parent_index;
            }
        }

        fn siftDown(self: *Self, start_index: u64) void {
            var index = start_index;
            const half = self.len >> 1;
            while (true) {
                var left_index = (index << 1) + 1;
                var right_index = left_index + 1;
                var left = if (left_index < self.len) self.heap[left_index] else null;
                var right = if (right_index < self.len) self.heap[right_index] else null;

                var smallest_index = index;
                var smallest = self.heap[index];

                if (left) |e| {
                    if (compareFn(self.context, e.priority, smallest.priority) == .lt) {
                        smallest_index = left_index;
                        smallest = e;
                    }
                }

                if (right) |e| {
                    if (compareFn(self.context, e.priority, smallest.priority) == .lt) {
                        smallest_index = right_index;
                        smallest = e;
                    }
                }

                if (smallest_index == index) return;

                const smallest_key = self.heap[smallest_index].key;
                const key = self.heap[index].key;
                self.heap[smallest_index] = self.heap[index];
                self.heap[index] = smallest;

                self.id_to_heap_pos.getPtr(smallest_key).?.* = index;
                self.id_to_heap_pos.getPtr(key).?.* = smallest_index;

                index = smallest_index;

                if (index >= half) return;
            }
        }

        pub fn getPriority(self: *Self, key: T) ?P {
            if (self.id_to_heap_pos.get(key)) |idx| {
                return self.heap[idx].priority;
            } else {
                return null;
            }
        }

        /// Free memory used by the queue.
        pub fn deinit(self: Self) void {
            self.allocator.free(self.heap);
            self.id_to_heap_pos.deinit(self.allocator);
        }

        /// Ensure that the queue can fit at least `new_capacity` items.
        pub fn ensureTotalCapacity(self: *Self, new_capacity: usize) !void {
            var better_capacity = self.capacity();
            if (better_capacity >= new_capacity) return;
            while (true) {
                better_capacity += better_capacity / 2 + 8;
                if (better_capacity >= new_capacity) break;
            }
            self.heap = try self.allocator.realloc(self.heap, better_capacity);
        }

        /// Ensure that the queue can fit at least `additional_count` **more** item.
        pub fn ensureUnusedCapacity(self: *Self, additional_count: usize) !void {
            try self.ensureTotalCapacity(self.len + additional_count);
        }
    };
}

fn lessThan(context: void, a: u32, b: u32) Order {
    _ = context;
    return std.math.order(a, b);
}

fn greaterThan(context: void, a: u32, b: u32) Order {
    return lessThan(context, a, b).invert();
}

const PQlt = UpdateablePriorityQueue(u32, u32, void, lessThan);
const PQgt = UpdateablePriorityQueue(u32, u32, void, greaterThan);

test "UpdateablePriorityQueue: add and remove min heap" {
    var queue = PQlt.init(std.testing.allocator, {});
    defer queue.deinit();

    const priorities = [_]u32{ 54, 12, 7, 23, 25, 13 };
    const expected_prio = [_]u32{ 7, 12, 13, 23, 25, 54 };
    const expected_keys = [_]u32{ 2, 1, 5, 3, 4, 0 };

    var i: u32 = 0;
    while (i < 6) {
        _ = try queue.updateOrInsert(i, priorities[i]);
        i += 1;
    }
    i = 0;
    while (i < 6) {
        const expected_key = expected_keys[i];
        const expected_priorty = expected_prio[i];
        var res = try queue.pop();
        const key = res.key;
        const priority = res.priority;
        try expectEqual(@as(u32, expected_key), key);
        try expectEqual(@as(u32, expected_priorty), priority);
        i += 1;
    }
}
//
test "UpdateablePriorityQueue: add and then update same key" {
    var queue = PQlt.init(std.testing.allocator, {});
    defer queue.deinit();

    var updated = try queue.updateOrInsert(0, 54);
    try expectEqual(updated, true);
    updated = try queue.updateOrInsert(1, 12);
    try expectEqual(updated, true);
    updated = try queue.updateOrInsert(2, 7);
    try expectEqual(updated, true);
    updated = try queue.updateOrInsert(3, 23);
    try expectEqual(updated, true);
    updated = try queue.updateOrInsert(4, 25);
    try expectEqual(updated, true);
    updated = try queue.updateOrInsert(5, 13);
    try expectEqual(updated, true);

    updated = try queue.updateOrInsert(0, 55);
    try expectEqual(updated, false);
    updated = try queue.updateOrInsert(1, 13);
    try expectEqual(updated, false);
    updated = try queue.updateOrInsert(2, 14);
    try expectEqual(updated, false);
    updated = try queue.updateOrInsert(3, 24);
    try expectEqual(updated, false);
    updated = try queue.updateOrInsert(4, 26);
    try expectEqual(updated, false);
    updated = try queue.updateOrInsert(5, 14);
    try expectEqual(updated, false);
}
//
test "UpdateablePriorityQueue: removeOrNull on empty" {
    var queue = PQlt.init(testing.allocator, {});
    defer queue.deinit();

    try expect(queue.removeOrNull() == null);
}

test "UpdateablePriorityQueue: edge case 3 elements" {
    var queue = PQlt.init(testing.allocator, {});
    defer queue.deinit();

    _ = try queue.updateOrInsert(0, 3);
    _ = try queue.updateOrInsert(1, 2);
    _ = try queue.updateOrInsert(2, 9);
    const first = try queue.pop();
    try expectEqual(@as(u32, 2), first.priority);
    const second = try queue.pop();
    try expectEqual(@as(u32, 3), second.priority);
    const third = try queue.pop();
    try expectEqual(@as(u32, 9), third.priority);
}
//
//
test "UpdateablePriorityQueue: sift up with odd indices" {
    var queue = PQlt.init(testing.allocator, {});
    defer queue.deinit();
    const items = [_]u32{ 15, 7, 21, 14, 13, 22, 12, 6, 7, 25, 5, 24, 11, 16, 15, 24, 2, 1 };
    var i: u32 = 0;
    for (items) |e| {
        _ = try queue.updateOrInsert(i, e);
        i += 1;
    }

    const sorted_items = [_]u32{ 1, 2, 5, 6, 7, 7, 11, 12, 13, 14, 15, 15, 16, 21, 22, 24, 24, 25 };
    for (sorted_items) |e| {
        const node = try queue.pop();
        try expectEqual(e, node.priority);
    }
}
