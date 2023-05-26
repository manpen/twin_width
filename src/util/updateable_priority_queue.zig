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
            return Self {
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
        id_to_heap_pos: []T,
        heap: []UpdateablePriorityQueueNode(T, P),
        len: usize,
        allocator: Allocator,
        context: Context,

        /// Initialize and return a priority queue.
        pub fn init(allocator: Allocator, context: Context) Self {
            return Self {
                .id_to_heap_pos = &[_]T{},
                .heap = &[_]T{},
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

        pub fn peek(self: *Self) ?UpdateablePriorityQueueNode {
            if (self.empty()) {
                return null;
            } else {
                return self.heap[0];
            }
        }

        pub fn pop(self: *Self) !UpdateablePriorityQueueNode {
            return self.removeAtHeapIdx(0);
        }

        pub fn removeOrNull(self: *Self) ?UpdateablePriorityQueueNode {
            if (self.empty()) {
                return null;
            } else {
                return self.removeAtHeapIdx(0);
            }
        }

        pub fn removeKey(self: *Self, key: T) ?UpdateablePriorityQueueNode {
            if ((self.empty()) or (self.key >= self.id_to_heap_pos.len) or (self.id_to_heap_pos[key]==KEY_NOT_FOUND)) {
                return null;
            } else {
                return self.removeAtHeapIdx(self.id_to_heap_pos[key]);
            }
        }

        fn removeAtHeapIdx(self: *Self, index: T) T {
            const last = self.heap[self.len - 1];
            const item = self.heap[index];
            self.id_to_heap_pos[item.key] = KEY_NOT_FOUND;
            self.id_to_heap_pos[last.key] = 0;

            self.heap[index] = last;
            self.len -= 1;
            if (index == 0) {
                siftDown(self, index);
            } else {
                const parent_index = ((index - 1) >> 1);
                const parent = self.items[parent_index];
                if (compareFn(self.context, last.priority, parent.priority) == .gt) {
                    self.siftDown(index);
                } else {
                    self.siftUp(index);
                }
            }

            return item;
        }


        fn extend_ids(self: *Self, new_max_key: T) !bool {
            const new_capacity = new_max_key+1;
            var old_capacity = self.id_to_heap_pos.len;
            var better_capacity = self.id_to_heap_pos.len;
            if (better_capacity >= new_capacity) {
                return false;
            }
            while (true) {
                better_capacity += better_capacity / 2 + 8;
                if (better_capacity >= new_capacity) break;
            }
            self.id_to_heap_pos = try self.allocator.realloc(self.id_to_heap_pos, better_capacity);
            const N = self.id_to_heap_pos.len;
            for (old_capacity..N) |i| { // set all new indices to a null value
                self.id_to_heap_pos[i] = KEY_NOT_FOUND;
            }
            return true;
        }

        pub fn updateOrInsert(self: *Self, key: T, priority: P) !bool {
            const did_extend = try self.extend_ids(key);
            if (did_extend or (self.id_to_heap_pos[key] == KEY_NOT_FOUND)) { // new key
                try self.ensureUnusedCapacity(1);
                self.addUnchecked(key, priority);
            } else { // old key, just update
                self.update(key, priority);
            }
        }

        fn capacity(self: Self) usize {
            return self.heap.len;
        }

        fn addUnchecked(self: *Self, elem: T) void {
            self.heap[self.len] = elem;
            siftUp(self, self.len);
            self.len += 1;
        }

        fn update(self: *Self, key: T, priority: P) !void {
            const idx = self.id_to_heap_pos[key];
            const old_priority = self.heap[self.id_to_heap_pos[key]].priority;
            switch (compareFn(self.context, priority, old_priority)) {
                .lt => siftUp(self, idx),
                .gt => siftDown(self, idx),
                .eq => {}, // Nothing to do as the items have equal priority
            }
        }

        fn siftUp(self: *Self, start_index: usize) void {
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

                self.id_to_heap_pos[child_key] = parent_index;
                self.id_to_heap_pos[parent_key] = child_index;
                
                child_index = parent_index;
            }
        }
        
        fn siftDown(self: *Self, start_index: usize) void {
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

                const smallest_key = self.self.heap[smallest_index].key;
                const key = self.self.heap[index].key;
                self.heap[smallest_index] = self.heap[index];
                self.heap[index] = smallest;

                self.id_to_heap_pos[smallest_key] = index;
                self.id_to_heap_pos[key] = smallest_index;

                index = smallest_index;

                if (index >= half) return;
            }
        }


        pub fn getPriority(self: *Self, key: T) ?P {
            if ((self.id_to_heap_pos.len >= key) or (self.id_to_heap_pos[key] == KEY_NOT_FOUND)) {
                return null;
            } else {
                const idx = self.id_to_heap_pos[key];
                return self.heap[idx].priority;
            }
        }

        /// Free memory used by the queue.
        pub fn deinit(self: Self) void {
            self.allocator.free(self.heap);
            self.allocator.free(self.id_to_heap_pos);
        }

        pub fn ensureMaxKey(self: *Self, max_key: T) !void {
            _ = try self.extend_ids(max_key);
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

    const priorities = [_]u32{54, 12, 7, 23, 25, 13};
    const expected_prio = [_]u32{7, 12, 13, 23, 25, 54};
    const expected_keys = [_]u32{2, 1, 5, 3, 4, 0};
    
    var i: u32 =0;
     while (i < 6) {
        _ = try queue.updateOrInsert(i, priorities[i]);
    }
    i=0;
    while (i < 6) {
        const expected_key = expected_keys[i];
        const expected_priorty = expected_prio[i];
        var res =  queue.pop();
        try expectEqual(@as(u32, expected_key), res.key);
        try expectEqual(@as(u32, expected_priorty), res.priority);
        i+=1;
    }
}

test "UpdateablePriorityQueue: add and remove same min heap" {
    var queue = PQlt.init(std.testing.allocator, {});
    defer queue.deinit();

    const priorities = [_]u32{54, 12, 7, 23, 25, 13};
    const expected_prio = [_]u32{7, 12, 13, 23, 25, 54};
    const expected_keys = [_]u32{2, 1, 5, 3, 4, 0};

    var i: u32 = 0;
    while (i < 6) {
        _ = try queue.updateOrInsert(i, priorities[i]);
        i+=1;
    }
    i=0;
    while (i < 6) {
        _ = try queue.updateOrInsert(i, priorities[i]+1);
        i+=1;
    }
    try expectEqual(queue.size(), 6);
    i=0;
    while (i < 6) {
        const expected_key = expected_keys[i];
        const expected_priorty = expected_prio[i];
        var res =  queue.pop();
        try expectEqual(@as(u32, expected_key), res.key);
        try expectEqual(@as(u32, expected_priorty+1), res.priority);
        i+=1;
    }
    try expectEqual(queue.empty(), true);
}

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
    _ = try  queue.updateOrInsert(2, 9);
    try expectEqual(@as(u32, 2), queue.pop().priority);
    try expectEqual(@as(u32, 3), queue.pop().priority);
    try expectEqual(@as(u32, 9), queue.pop().priority);
}


test "UpdateablePriorityQueue: sift up with odd indices" {
    var queue = PQlt.init(testing.allocator, {});
    defer queue.deinit();
    const items = [_]u32{ 15, 7, 21, 14, 13, 22, 12, 6, 7, 25, 5, 24, 11, 16, 15, 24, 2, 1 };
    var i: u32 = 0;
    for (items) |e| {
        _ = try queue.updateOrInsert(i, e);
        i+=1;
    }

    const sorted_items = [_]u32{ 1, 2, 5, 6, 7, 7, 11, 12, 13, 14, 15, 15, 16, 21, 22, 24, 24, 25 };
    for (sorted_items) |e| {
        try expectEqual(e, queue.pop().priority);
    }
}

test "UpdateablePriorityQueue: add and remove max heap" {
    var queue = PQgt.init(testing.allocator, {});
    defer queue.deinit();

    _ = try queue.updateOrInsert(0, 54);
    _ = try queue.updateOrInsert(1, 12);
    _ = try queue.updateOrInsert(2, 7);
    _ = try queue.updateOrInsert(3, 23);
    _ = try queue.updateOrInsert(4, 25);
    _ = try queue.updateOrInsert(5, 13);
    try expectEqual(@as(u32, 54), queue.pop().priority);
    try expectEqual(@as(u32, 25), queue.pop().priority);
    try expectEqual(@as(u32, 23), queue.pop().priority);
    try expectEqual(@as(u32, 13), queue.pop().priority);
    try expectEqual(@as(u32, 12), queue.pop().priority);
    try expectEqual(@as(u32, 7), queue.pop().priority);

    try expectEqual(true, queue.empty());
}

test "UpdateablePriorityQueue: add and remove same max heap" {
    var queue = PQgt.init(testing.allocator, {});
    defer queue.deinit();

    _ = try queue.updateOrInsert(0, 1);
    _ = try queue.updateOrInsert(0, 1);
    _ = try queue.updateOrInsert(2, 2);
    _ = try queue.updateOrInsert(2, 2);
    _ = try queue.updateOrInsert(1, 1);
    _ = try queue.updateOrInsert(1, 1);
    try expectEqual(@as(u32, 1), queue.pop().priority);
    try expectEqual(@as(u32, 1), queue.pop().priority);
    try expectEqual(@as(u32, 3), queue.pop().priority);
    
    try expectEqual(true, queue.empty());
}


test "UpdateablePriorityQueue: update min heap" {
    var queue = PQlt.init(testing.allocator, {});
    defer queue.deinit();

    _ = try queue.updateOrInsert(55, 55);
    _ = try queue.updateOrInsert(44, 44);
    _ = try queue.updateOrInsert(11, 11);
    _ = try queue.updateOrInsert(55, 4);
    _ = try queue.updateOrInsert(44, 1);
    try expectEqual(@as(u32, 44), queue.pop().key);
    try expectEqual(@as(u32, 55), queue.pop().key);
    try expectEqual(@as(u32, 11), queue.pop().key);

    try expectEqual(true, queue.empty());
}


test "UpdateablePriorityQueue: update max heap" {
    var queue = PQgt.init(testing.allocator, {});
    defer queue.deinit();

    defer queue.deinit();

    _ = try queue.updateOrInsert(55, 55);
    _ = try queue.updateOrInsert(44, 44);
    _ = try queue.updateOrInsert(11, 11);
    _ = try queue.updateOrInsert(55, 4);
    _ = try queue.updateOrInsert(44, 1);
    try expectEqual(@as(u32, 11), queue.pop().key);
    try expectEqual(@as(u32, 55), queue.pop().key);
    try expectEqual(@as(u32, 44), queue.pop().key);

    try expectEqual(true, queue.empty());
}
