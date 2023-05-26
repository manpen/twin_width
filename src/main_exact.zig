const std = @import("std");
const graph = @import("graph/graph.zig");
const bitset = @import("util/two_level_bitset.zig");
const pace = @import("pace_2023/pace_fmt.zig");
comptime {
    _ = @import("graph/graph.zig");
}
comptime {
    _ = @import("util/two_level_bitset.zig");
}
comptime {
    _ = @import("util/compressed_bitmap.zig");
}
comptime {
    _ = @import("graph/red_edge_stack.zig");
}
comptime {
    _ = @import("tww/contraction_sequence.zig");
}
comptime {
    _ = @import("tww/retraceable_contraction_sequence.zig");
}
comptime {
    _ = @import("pace_2023/pace_fmt.zig");
}
comptime {
    _ = @import("graph/bfs.zig");
}
comptime {
    _ = @import("graph/dfs.zig");
}
comptime {
    _ = @import("util/top_k_scorer.zig");
}
const builtin = @import("builtin");

const BestKnownError = error{
    NotFound,
};

fn load_best_known(filename: []const u8) !?u32 {
    var file = try std.fs.cwd().openFile("instances/best_known_solutions.csv", .{});
    defer file.close();

    var buf_reader = std.io.bufferedReader(file.reader());
    var in_stream = buf_reader.reader();

    var buf: [1024]u8 = undefined;
    while (try in_stream.readUntilDelimiterOrEof(&buf, '\n')) |line| {
        if (line.len < filename.len + 2) {
            continue;
        }

        if (std.mem.eql(u8, line[0..filename.len], filename)) {
            var tww = try std.fmt.parseInt(u32, line[filename.len + 1 ..], 10);
            return tww;
        }
    }

    return null;
}

pub fn inner_initial_solver(comptime T: type, allocator: std.mem.Allocator, filename: []const u8, short_name: []const u8) !T {
    std.debug.print("Start solver for {s}\n", .{filename});

    var timer = try std.time.Instant.now();
    var pace_part = try pace.Pace2023Fmt(T).fromFile(allocator, filename);
    defer pace_part.deinit(allocator);

    var loaded_graph = graph.Graph(T).loadFromPace(allocator, &pace_part) catch |err| {
        //Print error message if the graph could not be loaded std.debug.print("Could not load graph: {}", .{err});
        return err;
    };
    defer loaded_graph.deinit();

    var ccs = try loaded_graph.findAllConnectedComponents();
    defer {
        allocator.free(ccs);
    }

    var best_known: ?u32 = null;
    var test_mode = false;
    {
        var u: usize = 0;
        while (u + 8 < filename.len) : (u += 1) {
            if (filename[u] == '_' and filename[u + 1] == 't' and filename[u + 2] == 'w' and filename[u + 3] == 'w' and filename[u + 7] == '_') {
                best_known = try std.fmt.parseInt(u32, filename[u + 4 .. u + 7], 10);
                test_mode = true;
                break;
            }
        }
    }

    if (test_mode) {
        loaded_graph.forceExactSolverToProduceSolution();
    }

    if (best_known == null) {
        best_known = load_best_known(short_name) catch null;
    }

    const tww = loaded_graph.solveExact() catch |err| {
        std.debug.print("FAILED: Error {} in file {s}\n", .{ err, filename });
        return err;
    };

    const formatted = try std.fmt.allocPrint(loaded_graph.allocator, "{s}.solution", .{filename});
    defer loaded_graph.allocator.free(formatted);

    var now = try std.time.Instant.now();
    var elapsed = now.since(timer) / (1000 * 1000);

    var cmp: u8 = '?';
    if (best_known) |best| {
        cmp = '!';
        if (tww < best) {
            cmp = '<';
        } else if (tww == best) {
            cmp = '=';
        }
    }

    std.debug.print("{s:<25} | {:>8} | {:>8} | {:>4} {c} {?:>4} | {:>6}ms ({:>3} min)\n", .{ short_name, loaded_graph.number_of_nodes, loaded_graph.number_of_edges, tww, cmp, best_known, elapsed, elapsed / 60_000 });
    try loaded_graph.contraction.writeSolution(formatted);
    std.debug.print("Wrote solution to {s}\n", .{formatted});

    if (best_known) |value| {
        if (tww > value) {
            std.debug.print("FAILED: {s}\n", .{filename});
            @panic("INVALID SOLUTION");
        }
        std.debug.assert(tww <= value);
    }
    return tww;
}

pub fn initial_solver(allocator: std.mem.Allocator, filename: []const u8, short_name: []const u8) !u32 {
    var problem_size = try pace.loadPaceProblemHeaderFromFile(filename);
    if (problem_size.nodes <= std.math.maxInt(u8)) {
        return try inner_initial_solver(u8, allocator, filename, short_name);
    } else if (problem_size.nodes <= std.math.maxInt(u16)) {
        return try inner_initial_solver(u16, allocator, filename, short_name);
    } else {
        return try inner_initial_solver(u32, allocator, filename, short_name);
    }
}

pub fn lessThanU8(context: void, lhs: []u8, rhs: []u8) bool {
    _ = context;
    for (0..lhs.len) |i| {
        if (lhs[i] < rhs[i]) {
            return true;
        } else if (lhs[i] > rhs[i]) {
            return false;
        }
    }
    return false;
}

pub fn main() !void {
    //PACE 2023 submission
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();

    const target_directory = "instances/tiny";

    var allocator = gpa.allocator();

    const args = try std.process.argsAlloc(allocator);
    defer std.process.argsFree(allocator, args);

    if (args.len > 1) {
        _ = try initial_solver(allocator, args[1], args[1]);
        return;
    }

    std.debug.print("{s:<25} | {s:>8} | {s:>8} | {s:>4} {s} {s:>4} | {s:>6}\n", .{ "filename", "nodes", "edges", "tww", "?", "best", "time (ms)" });

    var dirIter = try std.fs.cwd().openIterableDir(target_directory, .{});
    defer dirIter.close();
    var dirit = dirIter.iterate();

    var file_list = try std.ArrayListUnmanaged([]u8).initCapacity(allocator, 100);
    while (try dirit.next()) |item| {
        if (item.kind == .File) {
            if (item.name.len >= 3 and std.mem.eql(u8, item.name[item.name.len - 3 ..], ".gr")) { // and std.mem.eql(u8, item.name[item.name.len - 12 .. item.name.len - 6], "exact_")) {
                try file_list.append(allocator, try std.fmt.allocPrint(allocator, "{s}", .{item.name}));
            }
        }
    }

    std.sort.sort([]u8, file_list.items, {}, lessThanU8);

    var cumulative: u32 = 0;
    for (file_list.items) |name| {
        var complete_name = try std.fmt.allocPrint(allocator, "{s}/{s}", .{ target_directory, name });
        defer allocator.free(complete_name);

        cumulative += try initial_solver(allocator, complete_name, name);
    }

    for (file_list.items) |name| {
        allocator.free(name);
    }
    file_list.deinit(allocator);
}
