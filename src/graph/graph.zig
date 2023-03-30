const std = @import("std");
pub const GraphError = error {
	FileNotFound,
	NotPACEFormat
};

pub const Graph = struct {
	number_of_nodes: u32,
	number_of_edges: u32,
	graph: []u8,

	pub fn load_from_pace(allocator: std.mem.Allocator, filename: []const u8) !Graph {
		const file = try std.fs.cwd().openFile(filename, .{});
		defer file.close();

		const file_size = (try file.stat()).size;

		const buffer = try allocator.alloc(u8, file_size+1);
		// Free the buffer if we exit with an error
		errdefer allocator.free(buffer);

		const read_size = try file.read(buffer);

		std.debug.assert(read_size == file_size);
		// Zero terminate the buffer
		buffer[read_size] = 0;
		
		var splits = std.mem.split(u8,buffer," ");

		const p = splits.next() orelse return GraphError.NotPACEFormat;
		const tww = splits.next() orelse return GraphError.NotPACEFormat;


		if(!std.mem.eql(u8,p,"p")) {
			return GraphError.NotPACEFormat;
		}

		if(!std.mem.eql(u8,tww,"tww")) {
			return GraphError.NotPACEFormat;
		}

		const number_of_nodes: u32 = std.fmt.parseInt(u32,splits.next() orelse return GraphError.NotPACEFormat,10) catch {
			return GraphError.NotPACEFormat;
		};
		var final_split = std.mem.split(u8,splits.next() orelse return GraphError.NotPACEFormat,"\n");
		const number_of_edges = if(final_split.next()) |edges| std.fmt.parseInt(u32,edges,10) catch {
				return GraphError.NotPACEFormat;
		} else return GraphError.NotPACEFormat;
		
		std.log.info("Loaded graph {s} with bytes len {}", .{filename,buffer.len});
		
		return Graph {
			.number_of_nodes = number_of_nodes,
			.number_of_edges = number_of_edges,
			.graph = buffer,
		};
	}

	pub fn deinit(self: Graph, allocator: std.mem.Allocator) void {
		allocator.free(self.graph);
	}
};

test "Check simple loading" {
	var allocator = std.heap.GeneralPurposeAllocator(.{}){};
	defer std.debug.assert(!allocator.deinit());

	var graph = try Graph.load_from_pace(allocator.allocator(),"instances/tiny/tiny001.gr");
	// Free ressources
	defer graph.deinit(allocator.allocator());

	try std.testing.expectEqual(graph.graph.len,49);
	try std.testing.expectEqual(graph.number_of_nodes,10);
	try std.testing.expectEqual(graph.number_of_edges,9);
}

test "Check simple failed loading" {
	var allocator = std.heap.GeneralPurposeAllocator(.{}){};
	// Should never fail since the ressources should be managed
	defer std.debug.assert(!allocator.deinit());

	var graph = Graph.load_from_pace(allocator.allocator(),"instances/tiny/tiny001.graph") catch |err| {
		// Loading should fail
		std.debug.assert(err == error.FileNotFound);
		return;
	};
	_ = graph;
	// Should never be reached
	std.debug.assert(false);
}
