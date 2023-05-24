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

pub fn loadPaceProblemHeaderFromFile(filename: []const u8) !PaceProblem2023 {
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

        const number_of_edges: u32 = std.fmt.parseInt(u32, splits.next() orelse return PaceError.NotPACEFormat, 10) catch {
            return PaceError.NotPACEFormat;
        };

        return PaceProblem2023{ .nodes = number_of_nodes, .edges = number_of_edges };
    }
    return PaceError.NotPACEFormat;
}

pub fn loadPaceProblemHeaderFromStdin() !PaceProblem2023 {
    var file = std.io.getStdIn();
    var line_splitter = file.reader();

    var local_buffer: [2048]u8 = undefined;
    var problem_line = try line_splitter.readUntilDelimiterOrEof(local_buffer[0..], '\n') orelse return PaceError.NotPACEFormat;

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

    return PaceProblem2023{ .nodes = number_of_nodes, .edges = number_of_edges };
}

pub const PaceProblem2023 = struct { nodes: u32, edges: u32 };

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
            buffer[read_size] = 0;
            return try Self.fromString(allocator, buffer);
        }

        pub fn fromStdin(allocator: std.mem.Allocator, problem: PaceProblem2023) !Self {
            var file = std.io.getStdIn();
            var unbuffered_line_splitter = file.reader();
            var buffered_reader = std.io.BufferedReader(32768, @TypeOf(unbuffered_line_splitter)){ .unbuffered_reader = unbuffered_line_splitter };
            var line_splitter = buffered_reader.reader();
            var local_buffer: [2048]u8 = undefined;

            var nodes: []PaceNode(T) = try allocator.alloc(PaceNode(T), problem.nodes);
            errdefer allocator.free(nodes);

            for (nodes) |*node| {
                node.edges = edge_list.ParametrizedUnsortedArrayList(T).init();
            }

            errdefer {
                for (nodes) |*node| {
                    node.edges.deinit(allocator);
                }
            }

            var edge_count: u32 = 0;
            while (try line_splitter.readUntilDelimiterOrEof(local_buffer[0..], '\n')) |line| {
                if (line.len == 0 or line[0] == 'c') {
                    // Skip comments
                    continue;
                }

                var find_first = if (std.mem.indexOf(u8, line, " ")) |split_index| split_index else break;
                const node_id = try std.fmt.parseInt(u32, line[0..find_first], 10) - 1;
                const second_node_id = try std.fmt.parseInt(u32, line[find_first + 1 ..], 10) - 1;

                if (node_id >= problem.nodes or second_node_id >= problem.nodes) {
                    return PaceError.NodeIdsTooLarge;
                } else if (node_id >= std.math.maxInt(T) or second_node_id >= std.math.maxInt(T)) {
                    return PaceError.NodeIdsTooLarge;
                }

                edge_count += 1;

                try nodes[node_id].edges.add(allocator, @intCast(T, second_node_id));
                try nodes[second_node_id].edges.add(allocator, @intCast(T, node_id));
                if (edge_count == problem.edges) {
                    break;
                }
            }

            return .{ .number_of_nodes = problem.nodes, .number_of_edges = problem.edges, .nodes = nodes };
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

            errdefer {
                for (nodes) |*node| {
                    node.edges.deinit(allocator);
                }
            }

            var edge_count: u32 = 0;
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

                if (edge_count >= number_of_edges) {
                    return PaceError.TooManyEdges;
                }
                edge_count += 1;

                try nodes[node_id].edges.add(allocator, @intCast(T, second_node_id));
                try nodes[second_node_id].edges.add(allocator, @intCast(T, node_id));
            }
            if (edge_count < number_of_edges) {
                return PaceError.TooFewEdges;
            }

            return .{ .number_of_nodes = number_of_nodes, .number_of_edges = number_of_edges, .nodes = nodes };
        }

        pub fn deinit(self: *const Self, allocator: std.mem.Allocator) void {
            allocator.free(self.nodes);
        }

        pub fn deinitAllIncludingEdges(self: *const Self, allocator: std.mem.Allocator) void {
            for (self.nodes) |*node| {
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
    defer pace.deinitAllIncludingEdges(gpa.allocator());
    try expectEqual(pace.number_of_nodes, 10);
    try expectEqual(pace.number_of_edges, 2);
}

test "Check pace edge count" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer std.debug.assert(!gpa.deinit());

    const pace = try Pace2023Fmt(u8).fromString(gpa.allocator(), "p tww 10 2\n2 1\n2 3\n");
    defer pace.deinitAllIncludingEdges(gpa.allocator());
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
    defer pace.deinitAllIncludingEdges(gpa.allocator());
    try expectEqual(pace.number_of_nodes, 10);
    try expectEqual(pace.number_of_edges, 2);

    try expectEqual(pace.nodes[0].edges.edges.items[0], 1);
    try expectEqual(pace.nodes[0].edges.edges.items.len, 1);

    try expectEqual(pace.nodes[1].edges.edges.items[0], 0);
    try expectEqual(pace.nodes[1].edges.edges.items[1], 2);
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
    defer pace.deinitAllIncludingEdges(gpa.allocator());

    try expectEqual(pace.number_of_nodes, 25);
    try expectEqual(pace.number_of_edges, 40);
}
