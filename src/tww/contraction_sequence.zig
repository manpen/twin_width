const std = @import("std");
const comptime_util = @import("../util/comptime_checks.zig");

pub const ContractionError = error{ ContractionOverflow, NoContractionLeft, Incomplete };

pub fn Contraction(comptime T: type) type {
    comptime if (!comptime_util.checkIfIsCompatibleInteger(T)) {
        @compileError("The type must either be u8,u16 or u32!");
    };
    return struct { erased: T, survivor: T };
}

pub fn ContractionSequence(comptime T: type) type {
    comptime if (!comptime_util.checkIfIsCompatibleInteger(T)) {
        @compileError("The type must either be u8,u16 or u32!");
    };
    return struct {
        const Self = @This();
        list: []Contraction(T),
        write_ptr: u32,
        graph_size: u32,

        pub inline fn init(allocator: std.mem.Allocator, graph_size: u32) !Self {
            std.debug.assert(graph_size > 0);
            const list = try allocator.alloc(Contraction(T), graph_size - 1);
            return .{ .list = list, .write_ptr = 0, .graph_size = graph_size };
        }

        pub const ContractionSequenceIterator = struct {
            index: u32,
            sequence: *const ContractionSequence(T),
            pub inline fn next(self: *ContractionSequenceIterator) ?Contraction(T) {
                if (self.index >= self.sequence.write_ptr) {
                    return null;
                }
                const item = self.sequence.list[self.index];
                self.index += 1;
                return item;
            }
        };

        pub inline fn iterator(self: *Self) ContractionSequenceIterator {
            return ContractionSequenceIterator{ .index = 0, .sequence = self };
        }

        pub fn deinit(self: *const Self, allocator: std.mem.Allocator) void {
            allocator.free(self.list);
        }

        pub inline fn reset(self: *Self) void {
            self.write_ptr = 0;
        }

        pub inline fn append(self: *Self, other: *const Self) ContractionError!void {
            if (self.write_ptr + other.write_ptr > self.list.len) {
                return ContractionError.ContractionOverflow;
            }
            std.mem.copy(Contraction(T), self.list[self.write_ptr..], other.list[0..other.write_ptr]);
            self.write_ptr += other.write_ptr;
        }

        pub inline fn copyInto(self: *Self, other: *const Self) ContractionError!void {
            std.debug.assert(self.graph_size >= other.graph_size);
            self.reset();
            try self.append(other);
        }

        pub inline fn isComplete(self: *const Self) bool {
            return self.graph_size - 1 == self.write_ptr;
        }

        pub inline fn addContraction(self: *Self, contraction: Contraction(T)) ContractionError!void {
            if (self.write_ptr >= self.list.len) {
                return ContractionError.ContractionOverflow;
            }

            self.list[self.write_ptr] = contraction;
            self.write_ptr += 1;
        }

        pub inline fn getLastContraction(self: *Self) ?Contraction(T) {
            if (self.write_ptr == 0) {
                return null;
            }
            return self.list[self.write_ptr - 1];
        }

        pub inline fn removeLastContraction(self: *Self) ContractionError!Contraction(T) {
            if (self.write_ptr == 0) {
                return ContractionError.NoContractionLeft;
            }
            self.write_ptr -= 1;
            return self.list[self.write_ptr];
        }

        pub inline fn writeSolution(self: *Self, filename: []const u8) !void {
            var file = try std.fs.cwd().createFile(filename, .{});
            defer file.close();
            var writer = file.writer();

            var iterate_seq = self.iterator();
            while (iterate_seq.next()) |seq| {
                try std.fmt.format(writer, "{d} {d}\n", .{ seq.survivor + 1, seq.erased + 1 });
            }
        }

        pub inline fn write_to_slice(self: *Self, slice: []u32) void {
            var iterate_seq = self.iterator();
            var j: usize = 0;
            while (iterate_seq.next()) |seq| {
                slice[j] = @as(u32, seq.survivor);
                slice[j + 1] = @as(u32, seq.erased);
                j += 2;
            }
        }

        pub inline fn writeSolutionToStdout(self: *Self) !void {
            var file = std.io.getStdOut();
            var buffered = std.io.BufferedWriter(32768, @TypeOf(file.writer())){ .unbuffered_writer = file.writer() };
            var writer = buffered.writer();

            var iterate_seq = self.iterator();
            while (iterate_seq.next()) |seq| {
                try std.fmt.format(writer, "{d} {d}\n", .{ seq.survivor + 1, seq.erased + 1 });
            }
            try buffered.flush();
        }
    };
}

const expectEqual = std.testing.expectEqual;

test "ContractionSequence: Check init and deinit" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer {
        std.debug.assert(!gpa.deinit());
    }
    var sequence = try ContractionSequence(u8).init(gpa.allocator(), 100);
    defer sequence.deinit(gpa.allocator());

    try expectEqual(sequence.list.len, 99);
}

test "ContractionSequence: Append until overflow" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer {
        std.debug.assert(!gpa.deinit());
    }
    var sequence = try ContractionSequence(u8).init(gpa.allocator(), 2);
    defer sequence.deinit(gpa.allocator());

    try sequence.addContraction(.{ .erased = 0, .survivor = 1 });
    try std.testing.expectError(ContractionError.ContractionOverflow, sequence.addContraction(.{ .erased = 1, .survivor = 0 }));
}

test "ContractionSequence: Append until end" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer {
        std.debug.assert(!gpa.deinit());
    }
    var sequence = try ContractionSequence(u8).init(gpa.allocator(), 2);
    defer sequence.deinit(gpa.allocator());

    try sequence.addContraction(.{ .erased = 0, .survivor = 1 });
    try expectEqual(sequence.isComplete(), true);
}

test "ContractionSequence: Append multiple" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer {
        std.debug.assert(!gpa.deinit());
    }
    var sequence = try ContractionSequence(u8).init(gpa.allocator(), 4);
    defer sequence.deinit(gpa.allocator());

    var sequence_2 = try ContractionSequence(u8).init(gpa.allocator(), 3);
    defer sequence_2.deinit(gpa.allocator());

    try sequence.addContraction(.{ .erased = 1, .survivor = 2 });
    try sequence_2.addContraction(.{ .erased = 0, .survivor = 1 });
    try expectEqual(sequence.isComplete(), false);
    try expectEqual(sequence_2.isComplete(), false);

    try sequence.append(&sequence_2);
    try sequence.addContraction(.{ .erased = 2, .survivor = 1 });

    try std.testing.expectError(ContractionError.ContractionOverflow, sequence.addContraction(.{ .erased = 1, .survivor = 0 }));

    try expectEqual(sequence.isComplete(), true);
    try expectEqual(sequence.list[0], .{ .erased = 1, .survivor = 2 });
    try expectEqual(sequence.list[1], .{ .erased = 0, .survivor = 1 });
    try expectEqual(sequence.list[2], .{ .erased = 2, .survivor = 1 });
}

test "ContractionSequence: Copy into" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer {
        std.debug.assert(!gpa.deinit());
    }
    var sequence = try ContractionSequence(u8).init(gpa.allocator(), 4);
    defer sequence.deinit(gpa.allocator());

    var sequence_2 = try ContractionSequence(u8).init(gpa.allocator(), 4);
    defer sequence_2.deinit(gpa.allocator());

    try sequence.addContraction(.{ .erased = 1, .survivor = 2 });
    try sequence_2.addContraction(.{ .erased = 0, .survivor = 1 });

    try sequence_2.copyInto(&sequence);

    try expectEqual(sequence_2.write_ptr, sequence.write_ptr);
    try expectEqual(sequence_2.list[0], sequence.list[0]);
}

test "ContractionSequence: Iterator" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer {
        std.debug.assert(!gpa.deinit());
    }
    var sequence = try ContractionSequence(u8).init(gpa.allocator(), 4);
    defer sequence.deinit(gpa.allocator());

    try sequence.addContraction(.{ .erased = 1, .survivor = 2 });
    try sequence.addContraction(.{ .erased = 2, .survivor = 3 });

    {
        var iterator = sequence.iterator();
        try expectEqual(iterator.next(), .{ .erased = 1, .survivor = 2 });
        try expectEqual(iterator.next(), .{ .erased = 2, .survivor = 3 });
        try expectEqual(iterator.next(), null);
    }

    sequence.reset();
    {
        var iterator = sequence.iterator();
        try expectEqual(iterator.next(), null);
    }
}
