const std = @import("std");
const graph = @import("graph/graph.zig");
const bitset = @import("util/two_level_bitset.zig");
const pace = @import("pace_2023/pace_fmt.zig");
comptime { _ = @import("graph/graph.zig"); }
comptime { _ = @import("util/two_level_bitset.zig"); }
comptime { _ = @import("util/compressed_bitmap.zig"); }
comptime { _ = @import("graph/red_edge_stack.zig"); }
comptime { _ = @import("tww/contraction_sequence.zig"); }
comptime { _ = @import("tww/retraceable_contraction_sequence.zig"); }
comptime { _ = @import("pace_2023/pace_fmt.zig"); }
comptime { _ = @import("graph/bfs.zig"); }
const builtin = @import("builtin");


pub fn inner_initial_solver(comptime T: type, filename: []const u8, short_name: []const u8, fixed: *std.heap.FixedBufferAllocator) !void {
	var allocator = fixed.allocator();
	var timer = try std.time.Instant.now();
	var loaded_graph = graph.Graph(T).loadFromPace(allocator,filename) catch |err| {
			//Print error message if the graph could not be loaded
			std.debug.print("Could not load graph: {}", .{err});
			return err;
	};
	defer loaded_graph.deinit();

	try loaded_graph.findAllConnectedComponents();

	//_ = loaded_graph.splitIntoConnectedComponents() catch |err| {
	//	std.debug.print("Could not split graph: {}", .{err});
	//	return err;
	//};

	//const tww = try loaded_graph.connected_components.items[loaded_graph.connected_components.items.len-1].solveGreedy(&loaded_graph);

	//var tww: u32 = 0;
	//var i:u32 = loaded_graph.number_of_nodes-1;
	//while(i > 0) {
	//	tww = std.math.max(try loaded_graph.addContraction(i-1,loaded_graph.number_of_nodes-1),tww);
	//	i-=1;
	//}

	var now = try std.time.Instant.now();
	var elapsed = now.since(timer)/(1000*1000);
	std.debug.print("{s:<25} | {:>8}MB | {:>8} | {:>8} | {:>8} | {:>6}ms\n",.{short_name,fixed.end_index/(1024*1024),loaded_graph.number_of_nodes,loaded_graph.number_of_edges,100,elapsed});

}

pub fn initial_solver(filename: []const u8, short_name: []const u8, fixed: *std.heap.FixedBufferAllocator) !void {
	var problem_size = try pace.loadPaceGraphVerticesSizeFromFile(filename);
	if(problem_size <= std.math.maxInt(u8)) {
		try inner_initial_solver(u8,filename, short_name, fixed);
	}
	else if(problem_size <= std.math.maxInt(u16)) {
		try inner_initial_solver(u16,filename, short_name, fixed);
	}
	else {
		try inner_initial_solver(u32,filename, short_name, fixed);
	}
}

pub fn main() !void {
		var gpa = std.heap.GeneralPurposeAllocator(.{}){};
		defer _ = gpa.deinit();

		var allocator = gpa.allocator();
		std.debug.print("{s:<25} | {s:>10} | {s:>8} | {s:>8} | {s:>8} | {s:>6}\n",.{"filename","memory","nodes","edges","tww","time (ms)"});

		var dirIter = try std.fs.cwd().openIterableDir("instances/heuristic-public",.{});
		defer dirIter.close();
		var dirit = dirIter.iterate();

		var hpa_alloc_buffer = try allocator.alloc(u8,1400*1024*1024);
		defer allocator.free(hpa_alloc_buffer);
		var fixed_alloc = std.heap.FixedBufferAllocator.init(hpa_alloc_buffer);

		while(try dirit.next()) |item| {
			if(item.kind == .File) {
					const completeName = try std.fmt.allocPrint(allocator,"instances/heuristic-public/{s}",.{item.name});
					defer allocator.free(completeName);

					try initial_solver(completeName,item.name,&fixed_alloc);
					fixed_alloc.reset();
			}
		}
}
