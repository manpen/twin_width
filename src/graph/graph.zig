const std = @import("std");
const bitset = @import("../util/two_level_bitset.zig");

pub const GraphError = error {
	FileNotFound,
	NotPACEFormat,
	GraphTooLarge,
	MisformedEdgeList,
	InvalidContractionOneNodeErased,
	ContractionOverflow
};


// Id is omitted must be stored outside of the struct
pub const Node = struct {
	const Self = @This();

	// Index into the edge list (Smaller than pointer since pointer has 8 Bytes)
	black_edges: []u32,
	red_edges: std.ArrayList(u32),
	scratch: u32,

	pub inline fn delete_black_edge(self: *Self, edge_id: u32) void {
		for(self.black_edges) |*edge| {
			if(edge.* == edge_id) {
				edge.* |= 0x80000000;
				return;
			}
		}
	}
	

	pub inline fn removeRedEdge(self: *Self, node_id: u32) void {
		var i:u32 = 0;
		while(i < self.red_edges.items.len) : (i+=1) {
			if(self.red_edges.items[i] == node_id) {
				_ = self.red_edges.swapRemove(i);
				return;
			}
		}
	}

	pub inline fn addRedEdge(self: *Self, node_id: u32) !bool {
		for(self.red_edges.items) |edge| {
			if(edge==node_id) {
				return false;
			}
		}
		try self.red_edges.append(node_id);
		return true;
	}

	pub inline fn addRedEdgeNoCheck(self: *Self, node_id: u32) !void {
		try self.red_edges.append(node_id);
	}

	pub inline fn is_deleted(edge_id: u32) bool {
		return (edge_id & 0x80000000) != 0;
	}
};

pub const Contraction = struct {
	erased: u32,
	survivor: u32
};

pub const ContractionSequence = struct {
	list: []Contraction,
	write_ptr: u32,
	allocator: std.mem.Allocator,

	pub inline fn init(graphSize: u32, allocator: std.mem.Allocator) !ContractionSequence {
		const mem = try allocator.alloc(Contraction,graphSize-1);
		return ContractionSequence {
			.list = mem,
			.write_ptr = 0,
			.allocator = allocator
		};
	}

	pub inline fn finished(self: *ContractionSequence) bool {
		return self.write_ptr == self.list.len-1;
	}

	pub inline fn addContraction(self: *ContractionSequence, deleted: u32, survivor: u32) GraphError!void {
		if (self.write_ptr == self.list.len) {
			return GraphError.ContractionOverflow;
		}
		self.list[self.write_ptr] = Contraction{.erased = deleted,.survivor= survivor};
		self.write_ptr += 1;
	}

	pub inline fn lastContraction(self: *ContractionSequence) ?Contraction {
		if(self.write_ptr == 0) {
			return null;
		}
		return self.list[self.write_ptr-1];
	}

	pub inline fn release(self: ContractionSequence) void {
		self.allocator.free(self.list);
	}
};

const Pace2023Fmt = struct {
	number_of_nodes: u32,
	number_of_edges: u32,
	pub fn fromString(string: []const u8) !Pace2023Fmt {
		var splits = std.mem.split(u8,string," ");

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

		const number_of_edges: u32 = std.fmt.parseInt(u32,splits.next() orelse return GraphError.NotPACEFormat,10) catch {
			return GraphError.NotPACEFormat;
		};

		return Pace2023Fmt { 
			.number_of_nodes = number_of_nodes,
			.number_of_edges = number_of_edges
		};
	}
};

pub const Graph = struct {
	number_of_nodes: u32,
	number_of_edges: u32,
	graph: []u8,
	node_list: []Node,
	scratch_bitset: bitset.FastBitSet,
	scratch_bitset_2: bitset.FastBitSet,
	allocator: std.mem.Allocator,
	contraction: ContractionSequence,
	
	pub inline fn addEdge(self: *Graph, u: u32, v: u32) void {
		std.debug.assert(u < self.node_list.len);
		std.debug.assert(v < self.node_list.len);

		const u_node = &self.node_list[u];
		const v_node = &self.node_list[v];

		std.debug.assert(u_node.scratch < u_node.black_edges.len);
		u_node.black_edges[u_node.scratch] = v;
		u_node.scratch+=1;

		std.debug.assert(v_node.scratch < v_node.black_edges.len);
		v_node.black_edges[v_node.scratch] = u;
		v_node.scratch+=1;
	}

	pub inline fn reverseLastContraction(self: *Graph) !void {
		_ = self;
	}

	pub inline fn addContraction(self: *Graph, erased: u32, survivor: u32) !u32 {
		// Add all red edges from the erased node which we do not have yet!

		for(self.node_list[erased].red_edges.items) |red| {
			if(red != survivor and try self.node_list[survivor].addRedEdge(red)) {
				_ = try self.node_list[red].addRedEdge(survivor);
			}
			self.node_list[red].removeRedEdge(erased);
		}


		var erasedPtr:u32 = 0;
		var survivorPtr:u32 = 0;

		const black_edges_erased = self.node_list[erased].black_edges;
		const black_edges_survivor = self.node_list[survivor].black_edges;


		var tww:u32 = 0;

		while(erasedPtr < black_edges_erased.len and survivorPtr < black_edges_survivor.len) {
			if(Node.is_deleted(black_edges_erased[erasedPtr])) {
				erasedPtr+=1;
				continue;
			}
			else if(Node.is_deleted(black_edges_survivor[survivorPtr])) {
				survivorPtr+=1;
				continue;
			}
			else if(black_edges_erased[erasedPtr] == survivor) {
				erasedPtr+=1;
				continue;
			}
			else if(black_edges_survivor[survivorPtr] == erased) {
				//Remove edge
				black_edges_survivor[survivorPtr] |= 0x80000000;
				survivorPtr+=1;
				continue;
			}
			if(black_edges_erased[erasedPtr] < black_edges_survivor[survivorPtr]) {
				try self.node_list[survivor].red_edges.append(black_edges_erased[erasedPtr]);
				self.node_list[black_edges_erased[erasedPtr]].delete_black_edge(erased);
				try self.node_list[black_edges_erased[erasedPtr]].addRedEdgeNoCheck(survivor);
				tww = std.math.max(@intCast(u32,self.node_list[black_edges_erased[erasedPtr]].red_edges.items.len),tww);
				erasedPtr+=1;
				// Advance to the first not erased edge
				continue;
			}
			else if(black_edges_erased[erasedPtr] > black_edges_survivor[survivorPtr]) {
				try self.node_list[survivor].red_edges.append(black_edges_survivor[survivorPtr]);
				self.node_list[black_edges_survivor[survivorPtr]].delete_black_edge(survivor);
				try self.node_list[black_edges_survivor[survivorPtr]].addRedEdgeNoCheck(survivor);
				tww = std.math.max(@intCast(u32,self.node_list[black_edges_survivor[survivorPtr]].red_edges.items.len),tww);
				black_edges_survivor[survivorPtr] |= 0x80000000;
				survivorPtr+=1;

				continue;
			}

			survivorPtr+=1;
			erasedPtr+=1;
		}

		while(survivorPtr < black_edges_survivor.len) : (survivorPtr+=1) {
			if(Node.is_deleted(black_edges_survivor[survivorPtr]) or black_edges_survivor[survivorPtr] == erased) {
				black_edges_survivor[survivorPtr] |= 0x80000000;
				continue;
			}
			try self.node_list[survivor].red_edges.append(black_edges_survivor[survivorPtr]);
			self.node_list[black_edges_survivor[survivorPtr]].delete_black_edge(survivor);
			try self.node_list[black_edges_survivor[survivorPtr]].addRedEdgeNoCheck(survivor);
			tww = std.math.max(@intCast(u32,self.node_list[black_edges_survivor[survivorPtr]].red_edges.items.len),tww);
			black_edges_survivor[survivorPtr] |= 0x80000000;
		}

		while(erasedPtr < black_edges_erased.len) : (erasedPtr+=1) {
			if(Node.is_deleted(black_edges_erased[erasedPtr]) or black_edges_erased[erasedPtr] == survivor) {
				continue;
			}
			try self.node_list[survivor].red_edges.append(black_edges_erased[erasedPtr]);
			self.node_list[black_edges_erased[erasedPtr]].delete_black_edge(erased);
			try self.node_list[black_edges_erased[erasedPtr]].addRedEdgeNoCheck(survivor);
			tww = std.math.max(@intCast(u32,self.node_list[black_edges_erased[erasedPtr]].red_edges.items.len),tww);
		}

		return std.math.max(@intCast(u32,self.node_list[survivor].red_edges.items.len),tww);
	}

	pub fn loadFromPace(allocator: std.mem.Allocator, filename: []const u8) !Graph {
		const file = try std.fs.cwd().openFile(filename, .{});
		defer file.close();

		const file_size = @intCast(usize,(try file.stat()).size);

		const buffer = try allocator.alloc(u8, file_size+1);
		// Free the buffer if we exit with an error
		errdefer allocator.free(buffer);

		const read_size = try file.read(buffer);

		std.debug.assert(read_size == file_size);

		// For the moment expect the graph to not exceed 50MB in size
		if(read_size>1024*1024*500) {
			return GraphError.GraphTooLarge;
		}

		// Zero terminate the buffer
		buffer[read_size] = 0;

		var line_splitter = std.mem.split(u8,buffer,"\n");

		var problem_line = while(line_splitter.next()) |line| {
			if(line[0] == 'p') {
				break line;
			}
		} else return GraphError.NotPACEFormat;
		var problem = try Pace2023Fmt.fromString(problem_line);

		var node_list = try allocator.alloc(Node,problem.number_of_nodes);
		
		for(node_list) |*node| {
			node.red_edges = std.ArrayList(u32).init(allocator);
			node.scratch = 0;
		}
		// Remove allocations on failure
		errdefer {
			for(node_list) |*node| {
				node.red_edges.deinit();
			}
			allocator.free(node_list);
		}

		std.log.info("Loaded graph {s} with bytes len {} nodes {} and edges {}", .{filename,buffer.len,problem.number_of_nodes,problem.number_of_edges});

		var graph = Graph {
			.number_of_nodes = problem.number_of_nodes,
			.number_of_edges = problem.number_of_edges,
			.graph = buffer,
			.node_list = node_list,
			.allocator = allocator,
			.scratch_bitset = try bitset.FastBitSet.initEmpty(problem.number_of_nodes,allocator),
			.scratch_bitset_2 = try bitset.FastBitSet.initEmpty(problem.number_of_nodes,allocator),
			.contraction = try ContractionSequence.init(problem.number_of_nodes,allocator)
		};
		errdefer {
			graph.scratch_bitset.deinit();
			graph.scratch_bitset_2.deinit();
			graph.contraction.release();
		}
		var line_splitter_first_pass = std.mem.split(u8,buffer,"\n");
		// Skip problem line
		while(line_splitter_first_pass.next()) |line| {
			if(line[0] == 'p') {
				break;
			}
		}

		while(line_splitter_first_pass.next()) |line| {
			if(line[0] == 'c') {
				// Skip comments
				continue;
			}

			var find_first = if(std.mem.indexOf(u8,line," ")) |split_index| split_index else break;
			const node_id = try std.fmt.parseInt(u32,line[0..find_first],10)-1;
			const second_node_id = try std.fmt.parseInt(u32,line[find_first+1..],10)-1;

			std.debug.assert(node_id <= problem.number_of_nodes);
			std.debug.assert(second_node_id <= problem.number_of_nodes);
			graph.node_list[node_id].scratch+=1;
			graph.node_list[second_node_id].scratch+=1;
		}
		for(graph.node_list) |*node| {
			node.black_edges = try allocator.alloc(u32,node.scratch);
			node.scratch = 0;
		}

		errdefer {
			for(graph.node_list) |*node| {
				allocator.free(node.black_edges);
			}
		}


		while(line_splitter.next()) |line| {
			if(line[0] == 'c') {
				// Skip comments
				continue;
			}

			var find_first = if(std.mem.indexOf(u8,line," ")) |split_index| split_index else break;
			const node_id = try std.fmt.parseInt(u32,line[0..find_first],10)-1;
			const second_node_id = try std.fmt.parseInt(u32,line[find_first+1..],10)-1;

			std.debug.assert(node_id <= problem.number_of_nodes);
			std.debug.assert(second_node_id <= problem.number_of_nodes);

			graph.addEdge(node_id,second_node_id);
		}
		// Sort the black edges
		for(graph.node_list) |*node| {
			std.sort.sort(u32,node.black_edges,{},comptime std.sort.asc(u32));
		}
		
		return graph;
	}

	pub fn release(self: *Graph) void {
		self.allocator.free(self.graph);
		for(self.node_list) |*node| {
			self.allocator.free(node.black_edges);
			node.red_edges.deinit();
		}
		self.contraction.release();
		self.scratch_bitset.deinit();
		self.scratch_bitset_2.deinit();
		self.allocator.free(self.node_list);
	}
};

test "Check simple loading" {
	var allocator = std.heap.GeneralPurposeAllocator(.{}){};
	defer std.debug.assert(!allocator.deinit());

	var graph = try Graph.loadFromPace(allocator.allocator(),"instances/tiny/tiny001.gr");
	// Free ressources
	defer graph.release();

	try std.testing.expectEqual(graph.graph.len,49);
	try std.testing.expectEqual(graph.number_of_nodes,10);
	try std.testing.expectEqual(graph.number_of_edges,9);
}

test "Test contraction" {
	var allocator = std.heap.GeneralPurposeAllocator(.{}){};
	defer std.debug.assert(!allocator.deinit());

	var graph = try Graph.loadFromPace(allocator.allocator(),"instances/tiny/tiny001.gr");
	// Free ressources
	defer graph.release();

	try std.testing.expectEqual(graph.graph.len,49);
	try std.testing.expectEqual(graph.number_of_nodes,10);
	try std.testing.expectEqual(graph.number_of_edges,9);

	var i: u32 = 9;
	while(i > 0) {
		var tww = try graph.addContraction(i-1,9);
		if(i>1) {
			try std.testing.expectEqual(@as(u32,1),tww);
		}
		else {
			try std.testing.expectEqual(@as(u32,0),tww);
		}
		i-=1;
	}
}

test "Check simple failed loading" {
	var allocator = std.heap.GeneralPurposeAllocator(.{}){};
	// Should never fail since the ressources should be managed
	defer std.debug.assert(!allocator.deinit());

	var graph = Graph.loadFromPace(allocator.allocator(),"instances/tiny/tiny001.graph") catch |err| {
		// Loading should fail with file not found
		try std.testing.expectEqual(err,error.FileNotFound);
		return;
	};
	_ = graph;
	// Should never be reached
	std.debug.assert(false);
}
