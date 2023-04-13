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
comptime { _ = @import("graph/dfs.zig"); }
comptime { _ = @import("util/top_k_scorer.zig"); }
const builtin = @import("builtin");


//pub fn inner_initial_solver(comptime T: type, filename: []const u8, short_name: []const u8) !void {
//	var node_allocator = std.heap.c_allocator;

//	var timer = try std.time.Instant.now();
//	var loaded_graph = graph.Graph(T).loadFromPace(node_allocator, node_allocator,filename) catch |err| {
			//Print error message if the graph could not be loaded
//			std.debug.print("Could not load graph: {}", .{err});
//			return err;
//	};
//	defer loaded_graph.deinit();

//	try loaded_graph.findAllConnectedComponents();
//	const tww = loaded_graph.solveGreedy() catch |err| {
//		std.debug.print("Error {}\n",.{err});
//		return err;
//	};


//	const formatted = try std.fmt.allocPrint(node_allocator,"{s}.solution",.{filename});
//	defer node_allocator.free(formatted);
//	try loaded_graph.contraction.seq.writeSolution(formatted);
//	std.debug.print("Wrote solution to {s}\n",.{formatted});

	//TODO: Check these instances again
	//REALLY SLOW: heuristic_122.gr better but still ~550 sec heuristic_136.gr slow too.
	//SLOW: heuristic_052.gr
	//BAD: heuristic_116.gr

	//_ = loaded_graph.splitIntoConnectedComponents() catch |err| {
	//	std.debug.print("Could not split graph: {}", .{err});
	//	return err;
	//};

	//const tww = try loaded_graph.connected_components.items[loaded_graph.connected_components.items.len-1].solveGreedy(&loaded_graph);

	//var tww: u32 = 0;
	//var i:T = @intCast(T,loaded_graph.number_of_nodes)-1;
	//while(i > 0) {
	//	tww = std.math.max(try loaded_graph.addContraction(i-1,@intCast(T,loaded_graph.number_of_nodes)-1),tww);
	//	i-=1;
	//}

//	var now = try std.time.Instant.now();
//	var elapsed = now.since(timer)/(1000*1000);
//	std.debug.print("{s:<25} | {:>8} | {:>8} | {:>8} | {:>6}ms\n",.{short_name,loaded_graph.number_of_nodes,loaded_graph.number_of_edges,tww,elapsed});

//}

pub fn inner_initial_solver_memory(comptime T: type, pace_header: pace.PaceProblem2023) !void {
	var node_allocator = std.heap.c_allocator;

	var pace_inst = try pace.Pace2023Fmt(T).fromStdin(node_allocator,pace_header);
	defer pace_inst.deinit(node_allocator);
	var loaded_graph = graph.Graph(T).loadFromPace(node_allocator,node_allocator,&pace_inst) catch |err| {
			//Print error message if the graph could not be loaded
			std.debug.print("Could not load graph: {}", .{err});
			return err;
	};
	defer loaded_graph.deinit();

	try loaded_graph.findAllConnectedComponents();
	_ = loaded_graph.solveGreedy() catch |err| {
		std.debug.print("Error {}\n",.{err});
		return err;
	};
	try loaded_graph.contraction.seq.writeSolutionToStdout();

	//pace_2023_submission
}

//pub fn initial_solver(filename: []const u8, short_name: []const u8) !void {
//	var problem_size = try pace.loadPaceGraphVerticesSizeFromFile(filename);
//	if(problem_size <= std.math.maxInt(u8)) {
//		try inner_initial_solver(u8,filename, short_name);
//	}
//	else if(problem_size <= std.math.maxInt(u16)) {
//		try inner_initial_solver(u16,filename, short_name);
//	}
//	else {
//		try inner_initial_solver(u32,filename, short_name);
//	}
//}

//pub fn lessThanU8(context: void, lhs: []u8, rhs: []u8) bool {
//	_ = context;
//	for(0..lhs.len) |i| {
//		if(lhs[i] < rhs[i]) {
//			return true;
//		}
//		else if(lhs[i] > rhs[i]) {
//			return false;
//		}
//	}
//	return false;
//}

pub fn main() !void {
	var problem = try pace.loadPaceProblemHeaderFromStdin();
	if(problem.nodes <= std.math.maxInt(u8)) {
		try inner_initial_solver_memory(u8,problem);
	}
	else if(problem.nodes <= std.math.maxInt(u16)) {
		try inner_initial_solver_memory(u16,problem);
	}
	else {
		try inner_initial_solver_memory(u32,problem);
	}
}

//pub fn main() !void {
		//PACE 2023 submission
//		var gpa = std.heap.GeneralPurposeAllocator(.{}){};
//		defer _ = gpa.deinit();

//		var allocator = gpa.allocator();
//		std.debug.print("{s:<25} | {s:>8} | {s:>8} | {s:>8} | {s:>6}\n",.{"filename","nodes","edges","tww","time (ms)"});

//		var dirIter = try std.fs.cwd().openIterableDir("instances/heuristic-public",.{});
//		defer dirIter.close();
//		var dirit = dirIter.iterate();


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
		//try initial_solver("instances/heuristic-public/heuristic_186.gr","heuristic_186.gr",&fixed_alloc);
		//fixed_alloc.reset();

//		var file_list = try std.ArrayListUnmanaged([]u8).initCapacity(allocator,100);
//		while(try dirit.next()) |item| {
//			if(item.kind == .File) {
//					try file_list.append(allocator, try std.fmt.allocPrint(allocator,"{s}",.{item.name}));
//			}
//		}

//		std.sort.sort([]u8, file_list.items, {}, lessThanU8);

//		for(file_list.items) |name| {
//			var complete_name = try std.fmt.allocPrint(allocator,"instances/heuristic-public/{s}",.{name});
//			defer allocator.free(complete_name);

//			try initial_solver(complete_name,name);
//		}

//		for(file_list.items) |name| {
//			allocator.free(name);
//		}
//		file_list.deinit(allocator);
//}
