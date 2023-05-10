const std = @import("std");
const mg = @import("exact/matrix_graph.zig");
const bnb = @import("exact/branch_and_bound.zig");

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    var allocator = gpa.allocator();

    const files = [_][]const u8{
        "instances/tiny/tiny001.gr",
        "instances/tiny/tiny002.gr",
        "instances/tiny/tiny003.gr",
        "instances/tiny/tiny004.gr",
        "instances/tiny/tiny005.gr",
        "instances/tiny/tiny006.gr",
        "instances/tiny/tiny007.gr",
        "instances/tiny/tiny008.gr",
        "instances/tiny/tiny009.gr",
        "instances/tiny/tiny010.gr",
    };

    const twws = [_]u32{ 1, 2, 0, 0, 3, 0, 2, 4, 1, 2 };

    for (files, twws) |path, correct_tww| {
        const file = try std.fs.cwd().openFile(path, .{});
        defer file.close();

        var graph = mg.MatrixGraph(32).new();
        var reader = file.reader();
        var buffer: [1024]u8 = undefined;

        while (try reader.readUntilDelimiterOrEof(buffer[0..], '\n')) |line| {
            if (line[0] == 'p') {
                continue;
            }

            var find_first = if (std.mem.indexOf(u8, line, " ")) |split_index| split_index else break;
            const u = try std.fmt.parseInt(u32, line[0..find_first], 10) - 1;
            const v = try std.fmt.parseInt(u32, line[find_first + 1 ..], 10) - 1;
            _ = graph.addEdge(u, v, mg.Color.Black);
        }

        std.debug.print("Opened {s} with {d} edges\n", .{ path, graph.number_of_edges() });

        if (true) {
            var solver = try bnb.ExactBranchAndBound.new(allocator, graph.number_of_nodes());
            defer solver.deinit();

            var tww = try solver.solve(@TypeOf(graph), graph, 0);
            std.debug.print("Final solution size {d} afert {d} recursive calls\n", .{ tww, solver.number_of_recursive_calls() });
            std.debug.assert(tww == correct_tww);
        }
    }
}
