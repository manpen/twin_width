pub const CompactField = struct {
    const Self = @This();
    bitset: *bitset.FastBitSet,
    field: std.BoundedArray(u32, 128),
    dense: bool,
    pub const CompactFieldIterator = struct {
        iterator: ?bitset.FastBitSetIterator,
        index: u32,
        field: *const CompactField,
        pub inline fn next(self: *CompactFieldIterator) ?u32 {
            if (self.iterator) |*iter| {
                return iter.next();
            } else if (self.index >= self.field.field.len) {
                return null;
            }
            const item = self.field.field.buffer[self.index];
            self.index += 1;
            return item;
        }
    };
    pub inline fn iterator(self: *const Self) CompactFieldIterator {
        if (self.dense) {
            return CompactFieldIterator{ .iterator = self.bitset.iter(), .index = 0, .field = self };
        } else {
            return CompactFieldIterator{ .iterator = null, .index = 0, .field = self };
        }
    }

    pub inline fn cardinality(self: *const Self) u32 {
        if (self.dense) {
            return self.bitset.cardinality;
        } else {
            return @intCast(u32, self.field.len);
        }
    }

    pub inline fn toBitset(self: *Self) void {
        for (self.field.buffer[0..self.field.len]) |item| {
            self.bitset.set(item);
        }
        self.dense = true;
        self.field.resize(0) catch {};
    }

    pub inline fn clear(self: *Self) void {
        if (self.dense) {
            self.bitset.unsetAll();
        } else {
            self.field.resize(0) catch {};
        }
        self.dense = false;
    }

    pub inline fn setUnionOnParam(self: *const Self, other: *bitset.FastBitSet) void {
        if (self.dense) {
            other.setUnion(self.bitset);
        } else {
            for (self.field.buffer[0..self.field.len]) |item| {
                other.set(item);
            }
        }
    }

    pub inline fn get(self: *Self, item: u32) bool {
        if (self.dense) {
            return self.bitset.get(item);
        } else {
            var i: u32 = 0;
            while (i < self.field.len) : (i += 1) {
                if (self.field.buffer[i] == item) {
                    return true;
                } else if (self.field.buffer[i] > item) {
                    return false;
                }
            }
            return false;
        }
    }

    pub inline fn unset(self: *Self, item: u32) void {
        if (self.dense) {
            _ = self.bitset.unset(item);
        } else {
            var i: u32 = 0;
            while (i < self.field.len) : (i += 1) {
                if (self.field.buffer[i] == item) {
                    _ = self.field.orderedRemove(i);
                    return;
                } else if (self.field.buffer[i] > item) {
                    return;
                }
            }
        }
    }

    pub inline fn add(self: *Self, item: u32) void {
        if (self.dense) {
            self.bitset.set(item);
        } else {
            var i: u32 = 0;
            while (i < self.field.len) : (i += 1) {
                if (self.field.buffer[i] == item) {
                    return;
                } else if (self.field.buffer[i] > item) {
                    self.field.insert(i, item) catch {
                        self.toBitset();
                        self.bitset.set(item);
                    };
                    return;
                }
            }
            self.field.append(item) catch {
                self.toBitset();
                self.bitset.set(item);
            };
        }
    }
};
