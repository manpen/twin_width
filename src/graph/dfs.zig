const std = @import("std");
const bitset = @import("../util/two_level_bitset.zig");
const Graph = @import("graph.zig").Graph;
const comptime_util = @import("../util/comptime_checks.zig");

pub const DfsKind = enum { red, black, both };

pub const DfsOptions = struct {
    max_level: u32 = std.math.maxInt(u32),
    kind: DfsKind = DfsKind.both,
};

pub fn DfsIterator(comptime T: type, comptime options: DfsOptions) type {
    comptime if (options.max_level < 1) {
        @compileError("BFS max level must be at least 1");
    };
    comptime if (!comptime_util.checkIfIsCompatibleInteger(T)) {
        @compileError("Type must either be u8, u16 or u32");
    };

		return struct {
			const Self = @This();
			visited: *bitset.FastBitSet,
			stack: std.ArrayListUnmanaged(T),
			depth_stack: std.ArrayListUnmanaged(T),
			parent_stack: std.ArrayListUnmanaged(T),
			graph: *const Graph(T),
			level: T,
			parent: T,

			pub inline fn next(self: *Self) ?T {
				while(self.stack.popOrNull()) |current| {
					self.level = self.depth_stack.pop();
					self.parent = self.parent_stack.pop();

					if(!self.visited.setExists(current)) {
						if (options.kind == .black or options.kind == .both) {
							var black_iter = self.graph.node_list[current].black_edges.iterator();
							while (black_iter.next()) |item| {
								self.stack.append(self.graph.allocator,item);
								self.depth_stack.append(self.graph.allocator,self.level+1);
								self.parent_stack.append(current);
							}
						}
						if (options.kind == .red or options.kind == .both) {
							var red_iter = self.graph.node_list[current].red_edges.iterator();
							while (red_iter.next()) |item| {
								self.stack.append(self.graph.allocator,item);
								self.depth_stack.append(self.graph.allocator,self.level+1);
								self.parent_stack.append(current);
							}
						}
						return current;
					}
				}
				return null;
			}

			pub inline fn deinit(self: *Self) void {
				self.stack.deinit(self.graph.allocator);
				self.depth_stack.deinit(self.graph.allocator);
				self.parent_stack.deinit(self.parent_stack);
			}
		};
}

pub inline fn dfs(comptime T: type, start_node: T, graph: *const Graph(T), visited: *bitset.FastBitSet, comptime options: DfsOptions) !DfsIterator(T, options) {
		var stack = try std.ArrayListUnmanaged(T).initCapacity(graph.allocator,graph.number_of_edges);
		var depth_stack = try std.ArrayListUnmanaged(T).initCapacity(graph.allocator,graph.number_of_edges);
		var parent_stack = try std.ArrayListUnmanaged(T).initCapacity(graph.allocator,graph.number_of_edges);

		try stack.append(graph.allocator,start_node);
		try depth_stack.append(graph.allocator,0);
		try parent_stack.append(graph.allocator,start_node);

    return DfsIterator(T, options){
        .graph = graph,
        .stack = stack,
        .visited = visited,
        .depth_stack = depth_stack,
				.parent_stack = parent_stack,
        .level = 0,
    };
}

pub const ArticulationPoints = struct {

};

pub fn findArticulationPoints(comptime T: type, root_node: T, graph: *const Graph(T), visited: *bitset.FastBitSet) !void {
		var parent_array = try graph.allocator.alloc(T,graph.number_of_nodes);
		defer graph.allocator.free(parent_array);

		var depth_array = try graph.allocator.alloc(T,graph.number_of_nodes);
		defer graph.allocator.free(depth_array);

		var low_array = try graph.allocator.alloc(T,graph.number_of_nodes);
		defer graph.allocator.free(low_array);

		var aritculation_points = try bitset.FastBitSet.initEmpty(graph.number_of_nodes,graph.allocator);
		defer aritculation_points.deinit(graph.allocator);

		var node_stack = try std.ArrayListUnmanaged(T).initCapacity(graph.allocator,2*graph.number_of_edges);
		var depth_stack = try std.ArrayListUnmanaged(T).initCapacity(graph.allocator,2*graph.number_of_edges);

		try node_stack.append(graph.allocator,root_node);
		try depth_stack.append(graph.allocator,0);

		var root_children:T = 0;

		var node_before:T = root_node;
		while(node_stack.popOrNull()) |item| {
			const current_depth = depth_stack.pop();

			// Item does not exist yet
			if(!visited.setExists(item)) {
				parent_array[item] = node_before;
				low_array[item] = current_depth;
				depth_array[item] = current_depth;

				if(node_before == root_node and item != root_node) {
					root_children+=1;
				}
				var black_iter = graph.node_list[item].black_edges.iterator();
				while (black_iter.next()) |child| {
					// Revisit ourselves after
					if(!visited.get(child)) {
						try node_stack.append(graph.allocator,item);
						try depth_stack.append(graph.allocator,current_depth);

						try node_stack.append(graph.allocator,child);
						try depth_stack.append(graph.allocator,current_depth+1);
					}
					else {
						low_array[item] = std.math.min(low_array[item],depth_array[node_before]);
					}
				}
			}
			else {
				// Ok we walked upwards
				if(parent_array[node_before] == item) {
					std.debug.print("\nCurrent {} is parent of {}\n",.{item,node_before});
					if(low_array[node_before] >= depth_array[item] and item != root_node) {
						aritculation_points.set(item);
					}
					low_array[item] = std.math.min(low_array[item],low_array[node_before]);
				}
			}
			node_before = item;
		}
		if(root_children > 1) {
			aritculation_points.set(root_node);
		}

		var articulationsp = aritculation_points.iter();
		while(articulationsp.next()) |item| {
			std.debug.print("Articulation point {}\n",.{item});
		}
}


test "Check articulation points" {
	var allocator = std.heap.GeneralPurposeAllocator(.{}){};
	var graph = try Graph(u8).new(5,allocator.allocator());
	try graph.addEdge(0,1);
	try graph.addEdge(0,2);
	try graph.addEdge(2,1);
	try graph.addEdge(0,3);
	try graph.addEdge(3,4);

	try findArticulationPoints(u8,0,&graph,&graph.scratch_bitset);
}
