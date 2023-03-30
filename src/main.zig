const std = @import("std");
const graph = @import("graph/graph.zig");
comptime { _ = @import("graph/graph.zig"); }

pub fn main() !void {
		// Generate new allocator to hold the graph
		var allocator = std.heap.GeneralPurposeAllocator(.{}){};
		// Check if the allocator is does not have any leaks or double free's at the end
		defer std.debug.assert(!allocator.deinit());


		// Load the graph from the file
		const loaded_graph = graph.Graph.load_from_pace(allocator.allocator(),"instances/tiny/tiny001.gr") catch |err| {
			//Print error message if the graph could not be loaded
			std.log.info("Could not load graph: {}\n", .{err});
			return;
		};
		// Free the graph at the end
		defer loaded_graph.deinit(allocator.allocator());
}

test "simple test" {
    var list = std.ArrayList(i32).init(std.testing.allocator);
    defer list.deinit(); // try commenting this out and see if zig detects the memory leak!
    try list.append(42);
    try std.testing.expectEqual(@as(i32, 42), list.pop());
}
