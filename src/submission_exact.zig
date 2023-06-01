const std = @import("std");
const graph = @import("graph/graph.zig");
const bitset = @import("util/two_level_bitset.zig");
const pace = @import("pace_2023/pace_fmt.zig");
const builtin = @import("builtin");
const signal_handler = @import("util/signal_handler.zig");

pub fn inner_initial_solver_memory(comptime T: type, allocator: std.mem.Allocator, pace_header: pace.PaceProblem2023) !void {
    var pace_inst = try pace.Pace2023Fmt(T).fromStdin(allocator, pace_header);
    defer pace_inst.deinit(allocator);
    var loaded_graph = graph.Graph(T).loadFromPace(allocator, &pace_inst) catch |err| {
        //Print error message if the graph could not be loaded
        std.debug.print("Could not load graph: {}", .{err});
        return err;
    };
    defer loaded_graph.deinit();

    var ccs = try loaded_graph.findAllConnectedComponents();
    signal_handler.initialize_signal_handler_slices_only(ccs);

    defer allocator.free(ccs);
    _ = loaded_graph.solveExact() catch |err| {
        std.debug.print("Error {}\n", .{err});
        return err;
    };
    try loaded_graph.contraction.writeSolutionToStdout();
}

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    var allocator = gpa.allocator();

    var problem = try pace.loadPaceProblemHeaderFromStdin();
    if (problem.nodes <= std.math.maxInt(u8)) {
        try inner_initial_solver_memory(u8, allocator, problem);
    } else if (problem.nodes <= std.math.maxInt(u16)) {
        try inner_initial_solver_memory(u16, allocator, problem);
    } else {
        try inner_initial_solver_memory(u32, allocator, problem);
    }
}
