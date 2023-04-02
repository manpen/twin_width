const std = @import("std");
const graph = @import("graph/graph.zig");
const bitset = @import("util/two_level_bitset.zig");
comptime { _ = @import("graph/graph.zig"); }
comptime { _ = @import("util/two_level_bitset.zig"); }
const builtin = @import("builtin");

pub fn main() !void {
		var gpa = std.heap.GeneralPurposeAllocator(.{}){};
		//var std_bitset = try std.bit_set.DynamicBitSet.initEmpty(gpa.allocator(),3_000_000);
		var std_bitset = try bitset.FastBitSet.initEmpty(3_000_000,gpa.allocator());


//		var i:u32 = 0;
		//while(i < 1000) : (i+=1) {
		//	std_bitset.set(@intCast(u32,rnd.next())%3_000_000);
		//	var iter = std_bitset.consuming_iter();
		//	while(iter.next()) |item| {
		//		std.debug.print("Item {}\n",.{item});
		//	}
		//}
		
		//var i:u32 = 0;
		//while(i < 100_000_000) : (i+=1) {
		//	std_bitset.set(1_500_000+i%1_400_000);
		//	std_bitset.unset(1_500_000+i%1_400_000);
		//}

		var i:u32 = 0;
		var calc:u64 = 0;
		while(i < 3_000_000) : (i+=1) {
			std_bitset.set(100);
			std_bitset.set(2_999_999);
			std_bitset.set(1_999_999);
			var iter = std_bitset.iter();
			//var iter = std_bitset.iterator(.{});
			while(iter.next()) |item| {
				calc+=@intCast(u64,item);
			}
			std_bitset.unset(100);
			std_bitset.unset(2_999_999);
			std_bitset.unset(1_999_999);
		}

		std.debug.print("Calc: {}\n",.{calc});

		//var allocation = try std.os.mmap(null, 512*1024*1024, std.os.PROT.READ | std.os.PROT.WRITE, std.os.MAP.PRIVATE | std.os.MAP.ANONYMOUS | std.os.MAP.POPULATE | 0x40, -1,0);
		//defer std.os.munmap(allocation);

		//var make_fixed_out_of_it = std.heap.FixedBufferAllocator.init(allocation);

		// Load the graph from the file
		//var loaded_graph = graph.Graph.loadFromPace(make_fixed_out_of_it.allocator(),"instances/heuristic-public/heuristic_002.gr") catch |err| {
			//Print error message if the graph could not be loaded
			//std.log.info("Could not load graph: {}", .{err});
			//return err;
		//};

		//var i: u32 = loaded_graph.number_of_nodes-1;
		//var tww: u32 = 0;
		//while(i > 0) {
		//	tww = std.math.max(try loaded_graph.addContraction(i-1,loaded_graph.number_of_nodes-1),tww);
		//	i-=1;
		//}
}

test "simple test" {
    var list = std.ArrayList(i32).init(std.testing.allocator);
    defer list.deinit(); // try commenting this out and see if zig detects the memory leak!
    try list.append(42);
    try std.testing.expectEqual(@as(i32, 42), list.pop());
}
