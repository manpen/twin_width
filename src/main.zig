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

pub fn inner_initial_solver(comptime T: type, allocator: std.mem.Allocator, filename: []const u8, short_name: []const u8) !T {
    // 138 OOM
    var timer = try std.time.Instant.now();
    var pace_part = pace.Pace2023Fmt(T).fromFile(allocator, filename) catch |err| {
        std.debug.print("Error load from file {}\n", .{err});
        return err;
    };
    var loaded_graph = graph.Graph(T).loadFromPace(allocator, &pace_part) catch |err| {
        //Print error message if the graph could not be loaded std.debug.print("Could not load graph: {}", .{err});
        std.debug.print("Error load graph {}\n", .{err});
        return err;
    };
    defer loaded_graph.deinit();

    _ = try loaded_graph.findAllConnectedComponents();
    const tww = loaded_graph.solveGreedy(.{ .single_pass = true }) catch |err| {
        std.debug.print("Error {}\n", .{err});
        return err;
    };

    const formatted = try std.fmt.allocPrint(loaded_graph.allocator, "{s}.solution", .{filename});
    defer loaded_graph.allocator.free(formatted);
    try loaded_graph.contraction.writeSolution(formatted);
    std.debug.print("Wrote solution to {s}\n", .{formatted});

    //TODO: Check these instances again
    //REALLY SLOW: heuristic_122.gr better but still ~550 sec heuristic_136.gr slow too.
    //SLOW: heuristic_052.gr
    //BAD: heuristic_116.gr

    var now = try std.time.Instant.now();
    var elapsed = now.since(timer) / (1000 * 1000);
    std.debug.print("{s:<25} | {:>8} | {:>8} | {:>8} | {:>6}ms\n", .{ short_name, loaded_graph.number_of_nodes, loaded_graph.number_of_edges, tww, elapsed });
    return tww;
}

pub fn initial_solver(allocator: std.mem.Allocator, filename: []const u8, short_name: []const u8) !u32 {
    var problem_size = pace.loadPaceProblemHeaderFromFile(filename) catch |err| {
        std.debug.print("Error {}\n", .{err});
        @panic("Caught error!");
    };
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

    const target_directory = "instances/heuristic-public";

    var allocator = gpa.allocator();
    std.debug.print("{s:<25} | {s:>8} | {s:>8} | {s:>8} | {s:>6}\n", .{ "filename", "nodes", "edges", "tww", "time (ms)" });

    var dirIter = try std.fs.cwd().openIterableDir(target_directory, .{});
    defer dirIter.close();
    var dirit = dirIter.iterate();

    var large_buffer = try allocator.alloc(u8, 1024 * 1024 * 7500);
    defer allocator.free(large_buffer);

    var hpa_allocator = std.heap.FixedBufferAllocator.init(large_buffer);
    var hpa = hpa_allocator.allocator();

    {
        const args = try std.process.argsAlloc(allocator);
        defer std.process.argsFree(allocator, args);

        if (args.len > 1) {
            _ = try initial_solver(hpa, args[1], args[1]);
            return;
        }
    }

    // HARD COLLECTION
    // SOLVABLE!
    //try initial_solver("instances/heuristic-public/heuristic_034.gr","heuristic_034.gr",&fixed_alloc);
    //try initial_solver("instances/heuristic-public/heuristic_052.gr","heuristic_052.gr",&fixed_alloc);
    //try initial_solver("instances/heuristic-public/heuristic_112.gr","heuristic_112.gr",&fixed_alloc);
    //try initial_solver("instances/heuristic-public/heuristic_138.gr","heuristic_138.gr",&fixed_alloc);
    //try initial_solver("instances/heuristic-public/heuristic_140.gr","heuristic_140.gr",&fixed_alloc);
    //try initial_solver("instances/heuristic-public/heuristic_144.gr","heuristic_144.gr",&fixed_alloc);
    //try initial_solver("instances/heuristic-public/heuristic_170.gr","heuristic_170.gr",&fixed_alloc);
    //fixed_alloc.reset();
    //try initial_solver("instances/heuristic-public/heuristic_184.gr","heuristic_184.gr",&fixed_alloc);
    //fixed_alloc.reset();
    //try initial_solver("instances/heuristic-public/heuristic_198.gr","heuristic_198.gr",&fixed_alloc);
    //fixed_alloc.reset();

    // NOT SOLVABLE!
    //try initial_solver("instances/heuristic-public/heuristic_122.gr","heuristic_122.gr",&fixed_alloc);
    //fixed_alloc.reset();
    //try initial_solver("instances/heuristic-public/heuristic_156.gr","heuristic_156.gr",&fixed_alloc);
    //fixed_alloc.reset();

    // TAKES AGES!
    // Update: Not anymore
    //try initial_solver(hpa,"instances/heuristic-public/heuristic_162.gr","heuristic_162.gr");
    //hpa_allocator.reset();

    //166 important
    //172 important

    //try initial_solver(hpa,"instances/heuristic-public/heuristic_172.gr","heuristic_172.gr");
    //hpa_allocator.reset();

    //try initial_solver(hpa,"instances/heuristic-public/heuristic_186.gr","heuristic_186.gr");
    //hpa_allocator.reset();

    //_ = try initial_solver(hpa,"instances/heuristic-public/heuristic_196.gr","heuristic_196.gr");
    //hpa_allocator.reset();

    //_ = try initial_solver(hpa,"instances/heuristic-public/heuristic_192.gr","heuristic_192.gr");
    //hpa_allocator.reset();

    //_ = try initial_solver(hpa,"instances/heuristic-public/heuristic_128.gr","heuristic_128.gr");
    //hpa_allocator.reset();

    //_ = try initial_solver(hpa,"instances/heuristic-public/heuristic_174.gr","heuristic_174.gr");
    //hpa_allocator.reset();

    //_ = try initial_solver(hpa,"instances/heuristic-public/heuristic_180.gr","heuristic_180.gr");
    //hpa_allocator.reset();

    // Not best instances/heuristic-public/heuristic_174.gr.solution

    //_ = try initial_solver(hpa,"instances/heuristic-public/heuristic_174.gr","heuristic_174.gr");
    //hpa_allocator.reset();

    //_ = try initial_solver(hpa,"instances/heuristic-public/heuristic_192.gr","heuristic_192.gr");
    //hpa_allocator.reset();

    //_ = try initial_solver(hpa,"instances/heuristic-public/heuristic_138.gr","heuristic_138.gr");
    //hpa_allocator.reset();

    //_ = try initial_solver(hpa,"instances/heuristic-public/heuristic_114.gr","heuristic_114.gr");
    //hpa_allocator.reset();
    //_ = try initial_solver(hpa,"instances/heuristic-public/heuristic_008.gr","heuristic_008.gr");
    //hpa_allocator.reset();

    var file_list = try std.ArrayListUnmanaged([]u8).initCapacity(allocator, 100);
    while (try dirit.next()) |item| {
        if (item.kind == .File) {
            if (item.name.len > 3 and std.mem.eql(u8, item.name[item.name.len - 3 ..], ".gr")) {
                try file_list.append(allocator, try std.fmt.allocPrint(allocator, "{s}", .{item.name}));
            }
        }
    }

    std.sort.sort([]u8, file_list.items, {}, lessThanU8);

    var cumulative: u32 = 3;
    var skip: u32 = 0;
    for (file_list.items) |name| {
        if (skip > 0) {
            skip -= 1;
            continue;
        }
        var complete_name = try std.fmt.allocPrint(allocator, "{s}/{s}", .{ target_directory, name });
        defer allocator.free(complete_name);

        cumulative += try initial_solver(hpa, complete_name, name);
        hpa_allocator.reset();
    }

    for (file_list.items) |name| {
        allocator.free(name);
    }
    file_list.deinit(allocator);
}
