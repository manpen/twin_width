const std = @import("std");
const graph = @import("graph/graph.zig");
const bitset = @import("util/two_level_bitset.zig");
comptime { _ = @import("graph/graph.zig"); }
comptime { _ = @import("util/two_level_bitset.zig"); }
comptime { _ = @import("util/compressed_bitmap.zig"); }
const builtin = @import("builtin");

pub fn main() !void {
		var gpa = std.heap.GeneralPurposeAllocator(.{}){};
		defer _ = gpa.deinit();

		var allocator = gpa.allocator();
		var allocation = try allocator.alloc(u8, 500*1024*1024);
		//var allocation = try std.os.mmap(null, 1000*1024*1024, std.os.PROT.READ | std.os.PROT.WRITE, std.os.MAP.PRIVATE | std.os.MAP.ANONYMOUS | 0x40, -1,0);
		//defer std.os.munmap(allocation);

		var make_fixed_out_of_it = std.heap.FixedBufferAllocator.init(allocation);
		std.debug.print("Allocated memory!\n",.{});

		//Load the graph from the file
		var loaded_graph = graph.Graph.loadFromPace(make_fixed_out_of_it.allocator(),"instances/heuristic-public/heuristic_200.gr") catch |err| {
			//Print error message if the graph could not be loaded
			std.log.info("Could not load graph: {}", .{err});
			return err;
		};
		std.debug.print("Loaded graph!\n",.{});

		var i: u32 = loaded_graph.number_of_nodes-1;
		var tww: u32 = 0;
		while(i > 0) {
			tww = std.math.max(try loaded_graph.addContraction(i-1,loaded_graph.number_of_nodes-1),tww);
			i-=1;
			std.debug.print("I {} twin width {}\n",.{i,tww});
		}
		std.debug.print("tww {}\n",.{tww});
}

test "simple test" {
    var list = std.ArrayList(i32).init(std.testing.allocator);
    defer list.deinit(); // try commenting this out and see if zig detects the memory leak!
    try list.append(42);
    try std.testing.expectEqual(@as(i32, 42), list.pop());
}
