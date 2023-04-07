const std = @import("std");
const edge_list = @import("../graph/edge_list.zig");
const comptime_util = @import("../util/comptime_checks.zig");

const PaceError = error{ NotPACEFormat, TooManyEdges, TooFewEdges, NodeIdsTooLarge };

pub fn loadPaceGraphVerticesSizeFromString(string: []const u8) !u32 {
    var line_splitter = std.mem.split(u8, string, "\n");
    while (line_splitter.next()) |line| {
        if (line.len > 0 and line[0] == 'c') {
            continue;
        }

        var splits = std.mem.split(u8, line, " ");

        const p = splits.next() orelse return PaceError.NotPACEFormat;
        const tww = splits.next() orelse return PaceError.NotPACEFormat;

        if (!std.mem.eql(u8, p, "p")) {
            return PaceError.NotPACEFormat;
        }

        if (!std.mem.eql(u8, tww, "tww")) {
            return PaceError.NotPACEFormat;
        }

        const number_of_nodes: u32 = std.fmt.parseInt(u32, splits.next() orelse return PaceError.NotPACEFormat, 10) catch {
            return PaceError.NotPACEFormat;
        };
        return number_of_nodes;
    }
    return PaceError.NotPACEFormat;
}

pub fn loadPaceGraphVerticesSizeFromFile(filename: []const u8) !u32 {
    const file = try std.fs.cwd().openFile(filename, .{});
    defer file.close();

    var reader = file.reader();
    var buffer: [1024]u8 = undefined;
    while (try reader.readUntilDelimiterOrEof(buffer[0..], '\n')) |line| {
        if (line.len > 0 and line[0] == 'c') {
            continue;
        }

        var splits = std.mem.split(u8, line, " ");

        const p = splits.next() orelse return PaceError.NotPACEFormat;
        const tww = splits.next() orelse return PaceError.NotPACEFormat;

        if (!std.mem.eql(u8, p, "p")) {
            return PaceError.NotPACEFormat;
        }

        if (!std.mem.eql(u8, tww, "tww")) {
            return PaceError.NotPACEFormat;
        }

        const number_of_nodes: u32 = std.fmt.parseInt(u32, splits.next() orelse return PaceError.NotPACEFormat, 10) catch {
            return PaceError.NotPACEFormat;
        };
        return number_of_nodes;
    }
    return PaceError.NotPACEFormat;
}

pub fn PaceNode(comptime T: type) type {
    comptime if (!comptime_util.checkIfIsCompatibleInteger(T)) {
        @compileError("T must either be u8,u16 or u32!");
    };

    return struct { edges: edge_list.ParametrizedUnsortedArrayList(T) };
}

pub fn Pace2023Fmt(comptime T: type) type {
    comptime if (!comptime_util.checkIfIsCompatibleInteger(T)) {
        @compileError("T must either be u8,u16 or u32!");
    };

    return struct {
        const Self = @This();
        number_of_nodes: u32,
        number_of_edges: u32,
        nodes: []PaceNode(T),
        pub fn fromFile(allocator: std.mem.Allocator, filename: []const u8) !Self {
            const file = try std.fs.cwd().openFile(filename, .{});
            defer file.close();

            const file_size = @intCast(usize, (try file.stat()).size);

            const buffer = try allocator.alloc(u8, file_size + 1);
            // Free the buffer if we exit with an error
            defer allocator.free(buffer);

            const read_size = try file.read(buffer);

            std.debug.assert(read_size == file_size);
            return try Self.fromString(allocator, buffer);
        }

        pub fn fromString(allocator: std.mem.Allocator, string: []const u8) !Self {
            var line_splitter = std.mem.split(u8, string, "\n");

            var problem_line = line_splitter.next() orelse return PaceError.NotPACEFormat;

            var splits = std.mem.split(u8, problem_line, " ");

            const p = splits.next() orelse return PaceError.NotPACEFormat;
            const tww = splits.next() orelse return PaceError.NotPACEFormat;

            if (!std.mem.eql(u8, p, "p")) {
                return PaceError.NotPACEFormat;
            }

            if (!std.mem.eql(u8, tww, "tww")) {
                return PaceError.NotPACEFormat;
            }

            const number_of_nodes: u32 = std.fmt.parseInt(u32, splits.next() orelse return PaceError.NotPACEFormat, 10) catch {
                return PaceError.NotPACEFormat;
            };

            const number_of_edges: u32 = std.fmt.parseInt(u32, splits.next() orelse return PaceError.NotPACEFormat, 10) catch {
                return PaceError.NotPACEFormat;
            };

            var nodes: []PaceNode(T) = try allocator.alloc(PaceNode(T), number_of_nodes);
            errdefer allocator.free(nodes);

            for (nodes) |*node| {
                node.edges = edge_list.ParametrizedUnsortedArrayList(T).init();
            }

            while (line_splitter.next()) |line| {
                if (line.len == 0 or line[0] == 'c') {
                    // Skip comments
                    continue;
                }

                var find_first = if (std.mem.indexOf(u8, line, " ")) |split_index| split_index else break;
                const node_id = try std.fmt.parseInt(u32, line[0..find_first], 10) - 1;
                const second_node_id = try std.fmt.parseInt(u32, line[find_first + 1 ..], 10) - 1;

                if (node_id >= number_of_nodes or second_node_id >= number_of_nodes) {
                    return PaceError.NodeIdsTooLarge;
                } else if (node_id >= std.math.maxInt(T) or second_node_id >= std.math.maxInt(T)) {
                    return PaceError.NodeIdsTooLarge;
                }

                try nodes[node_id].edges.add(allocator, @intCast(T,second_node_id));
                try nodes[second_node_id].edges.add(allocator, @intCast(T,node_id));
            }

            return .{ .number_of_nodes = number_of_nodes, .number_of_edges = number_of_edges, .nodes = nodes };
        }

        pub fn deinit(self: *const Self, allocator: std.mem.Allocator) void {
            allocator.free(self.nodes);
        }

        pub fn deinitAllIncludingEdges(self: *const Self, allocator: std.mem.Allocator) void {
            for (self.nodes) |node| {
                node.edges.deinit(allocator);
            }
            allocator.free(self.nodes);
        }
    };
}

const expectEqual = std.testing.expectEqual;

test "Check pace format" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer std.debug.assert(!gpa.deinit());

    const pace = try Pace2023Fmt(u8).fromString(gpa.allocator(), "p tww 10 2\n2 1\n2 3\n");
    defer pace.deinit(gpa.allocator());
    try expectEqual(pace.number_of_nodes, 10);
    try expectEqual(pace.number_of_edges, 2);
}

test "Check pace edge count" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer std.debug.assert(!gpa.deinit());

    const pace = try Pace2023Fmt(u8).fromString(gpa.allocator(), "p tww 10 2\n2 1\n2 3\n");
    defer pace.deinit(gpa.allocator());
    try expectEqual(pace.number_of_nodes, 10);
    try expectEqual(pace.number_of_edges, 2);

    try expectEqual(pace.nodes[0].edges.edges.items.len, 1);
    try expectEqual(pace.nodes[1].edges.edges.items.len, 2);
    try expectEqual(pace.nodes[2].edges.edges.items.len, 1);
}

test "Check pace edges" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer std.debug.assert(!gpa.deinit());

    const pace = try Pace2023Fmt(u8).fromString(gpa.allocator(), "p tww 10 2\n2 1\n2 3\n");
    defer pace.deinit(gpa.allocator());
    try expectEqual(pace.number_of_nodes, 10);
    try expectEqual(pace.number_of_edges, 2);

    try expectEqual(pace.nodes[0].edges.edges.items[0], 1);
    try expectEqual(pace.nodes[0].edges.edges.items.len, 1);

    try expectEqual(pace.nodes[1].edges.edges.items[1], 0);
    try expectEqual(pace.nodes[1].edges.edges.items[2], 2);
    try expectEqual(pace.nodes[1].edges.edges.items.len, 2);

    try expectEqual(pace.nodes[2].edges.edges.items[0], 1);
    try expectEqual(pace.nodes[2].edges.edges.items.len, 1);
}

test "Check pace too few edges" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer std.debug.assert(!gpa.deinit());

    try std.testing.expectError(PaceError.TooFewEdges, Pace2023Fmt(u8).fromString(gpa.allocator(), "p tww 10 2\n2 1\n"));
}

test "Check pace too many edges" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer std.debug.assert(!gpa.deinit());

    try std.testing.expectError(PaceError.TooManyEdges, Pace2023Fmt(u8).fromString(gpa.allocator(), "p tww 10 2\n2 1\n1 2\n2 2\n"));
}

test "Check pace node id too large" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer std.debug.assert(!gpa.deinit());

    try std.testing.expectError(PaceError.NodeIdsTooLarge, Pace2023Fmt(u8).fromString(gpa.allocator(), "p tww 10 2\n20 1\n1 2\n"));
}

test "Check pace load tiny 005" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer std.debug.assert(!gpa.deinit());

    const pace = try Pace2023Fmt(u8).fromFile(gpa.allocator(), "test/pace_2023/tiny005.gr");
    defer pace.deinit(gpa.allocator());

    try expectEqual(pace.number_of_nodes, 25);
    try expectEqual(pace.number_of_edges, 40);
}
