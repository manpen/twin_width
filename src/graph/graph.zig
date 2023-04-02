const std = @import("std");
pub const GraphError = error {
	FileNotFound,
	NotPACEFormat,
	GraphTooLarge,
	MisformedEdgeList,
	InvalidContractionOneNodeErased,
	ContractionOverflow
};


const EdgeColor = enum(u32) {
	black = 0,
	_,
	pub inline fn isRed(self: EdgeColor) bool {
		return switch(self) {
			.black => false,
			_ => true
		};
	}
};

pub const Edge = struct {
	target: u32,
	color: EdgeColor
};

pub fn edgeCmp(context: usize, a: Edge, b: Edge) bool {
	_ = context;
	return a.target < b.target;
}

pub const EdgeList = struct {
	const Self = @This();

	deleted: bool,
	red_degree: u32,
	// Target node and color
	list: std.ArrayList(Edge),
	pub inline fn removeEdge(self: *Self, target: u32) u32 {
		var i: u32 = 0;
		while(i < self.list.items.len): (i+=1) {
			if(self.list.items[i].target == target) {
				if(self.list.items[i].color.isRed()) {
					self.red_degree-=1;
				}
				_ = self.list.orderedRemove(i);
				return self.red_degree;
			}
		}
		// Should never happen
		std.debug.assert(false);
		return self.red_degree;
	}
	
	pub inline fn addEdge(self: *Self, edge: Edge) !u32 {
		var i: u32 = 0;
		if(edge.color.isRed()) {
			self.red_degree+=1;
		}

		while(i < self.list.items.len): (i+=1) {
			if(self.list.items[i].target > edge.target) {
				try self.list.insert(i,edge);
				return self.red_degree;
			}
		}
		try self.list.append(edge);
		return self.red_degree;
	}

	pub inline fn replaceColor(self: *Self, edge: Edge) u32 {
		var i: u32 = 0;

		while(i < self.list.items.len): (i+=1) {
			if(self.list.items[i].target == edge.target) {
				if(!self.list.items[i].color.isRed()) {
					// Increase red degree if the edge is already red this is a no op
					self.red_degree+=1;
					self.list.items[i].color = edge.color;
				}
				return self.red_degree;
			}
		}
		return self.red_degree;
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
		//TODO: Throw an error if this is overflowing!
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
	adj_list: []EdgeList,
	alive_nodes: std.ArrayList(u32),
	allocator: std.mem.Allocator,
	contraction: ContractionSequence,
	
	pub inline fn addEdge(self: *Graph, from: u32, to: u32, color: EdgeColor) !void {
		try self.adj_list[from].list.append(Edge {
			.target = to,
			.color = color
		});
		try self.adj_list[to].list.append(Edge {
			.target = from,
			.color = color
		});
	}

	pub inline fn reverseLastContraction(self: *Graph) !void {
		_ = self;
	}

	pub inline fn addContraction(self: *Graph, erased: u32, survivor: u32) !u32 {
		try self.contraction.addContraction(erased,survivor);
		if (self.adj_list[erased].deleted or self.adj_list[survivor].deleted) {
			return GraphError.InvalidContractionOneNodeErased;
		}
		self.adj_list[erased].deleted = true;


		var ptrErased : u32 = 0;
		var ptrSurvivor: u32 = 0;
		var listErased = &self.adj_list[erased].list;
		var listSurvivor = &self.adj_list[survivor].list;
		
		var tww: u32 = 0;
		while((ptrErased < listErased.items.len) and (ptrSurvivor < listSurvivor.items.len)) {
			// Drop the edge
			if (listErased.items[ptrErased].target == survivor) {
				ptrErased+=1;
				continue;
			}
			else if(listSurvivor.items[ptrSurvivor].target == erased) {
				if(listSurvivor.orderedRemove(ptrSurvivor).color.isRed()) {
					self.adj_list[survivor].red_degree-=1;
				}
				continue;
			}
			else if(listErased.items[ptrErased].target < listSurvivor.items[ptrSurvivor].target) {
				// Add new red edge
				const new_edge = Edge {.target = listErased.items[ptrErased].target, .color = @intToEnum(EdgeColor,erased+1)};
				try listSurvivor.insert(ptrSurvivor,new_edge);

				const target = listErased.items[ptrErased].target;
				_ = self.adj_list[target].removeEdge(erased);

				tww = std.math.max(try self.adj_list[target].addEdge(Edge {
					.target = survivor,
					.color = @intToEnum(EdgeColor,erased+1)
				}),tww);

				// Increase red edge
				self.adj_list[survivor].red_degree+=1;
				ptrErased+=1;
				ptrSurvivor+=1;
			}
			else if(listErased.items[ptrErased].target > listSurvivor.items[ptrSurvivor].target) {
				if(!listSurvivor.items[ptrSurvivor].color.isRed()) {
					// Color old edge red
					listSurvivor.items[ptrSurvivor].color = @intToEnum(EdgeColor,erased+1);
					self.adj_list[survivor].red_degree+=1;
					// Colors the edge red if it is not yet.
					tww = std.math.max(self.adj_list[listSurvivor.items[ptrSurvivor].target].replaceColor(Edge {
								.target = survivor,
								.color = @intToEnum(EdgeColor,erased+1)
								}),tww);
				}
				ptrSurvivor+=1;
			}
			else {
				const target = listErased.items[ptrErased].target;
				// Remove one edge
				_ = self.adj_list[target].removeEdge(erased);

				// Recolor our own edge only if it was black before
				if(listErased.items[ptrErased].color.isRed() and !listSurvivor.items[ptrSurvivor].color.isRed()) {
					// Color the other side of the edge red
					listSurvivor.items[ptrSurvivor].color = @intToEnum(EdgeColor,erased+1);
					self.adj_list[survivor].red_degree+=1;
					tww = std.math.max(self.adj_list[target].replaceColor(listSurvivor.items[ptrSurvivor]),tww);
				}
				ptrSurvivor+=1;
				ptrErased+=1;
			}
		}

		if(ptrErased < listErased.items.len) {
			while(ptrErased < listErased.items.len) {
				if(listErased.items[ptrErased].target == survivor) {
					// Just skip it
					ptrErased+=1;
					continue;
				}

				const target = listErased.items[ptrErased].target;
				const new_edge = Edge {.target = target, .color = @intToEnum(EdgeColor,erased+1)};
				try listSurvivor.append(new_edge);

				_ = self.adj_list[target].removeEdge(erased);

				tww = std.math.max(try self.adj_list[target].addEdge(Edge {
					.target = survivor,
					.color = @intToEnum(EdgeColor,erased+1)
				}),tww);

				// Increase red edge
				self.adj_list[survivor].red_degree+=1;
				ptrErased+=1;
			}
		}
		else if (ptrSurvivor < listSurvivor.items.len) {
			while(ptrSurvivor < listSurvivor.items.len) {
				if(listSurvivor.items[ptrSurvivor].target == erased) {
					if(listSurvivor.orderedRemove(ptrSurvivor).color.isRed()) {
						self.adj_list[survivor].red_degree-=1;
					}
					continue;
				}
				// Color old edge red
				const target = listSurvivor.items[ptrSurvivor].target;
				listSurvivor.items[ptrSurvivor].color = @intToEnum(EdgeColor,erased+1);
				self.adj_list[survivor].red_degree+=1;
				// Colors the edge red if it is not yet.
				tww = std.math.max(self.adj_list[target].replaceColor(Edge {
					.target = survivor,
					.color = @intToEnum(EdgeColor,erased+1)
				}),tww);
				ptrSurvivor+=1;
			}
		}

		return std.math.max(tww,self.adj_list[survivor].red_degree);
	}

	pub fn loadFromPace(allocator: std.mem.Allocator, filename: []const u8) !Graph {
		const file = try std.fs.cwd().openFile(filename, .{});
		defer file.close();

		const file_size = (try file.stat()).size;

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

		var edge_list = try allocator.alloc(EdgeList,problem.number_of_nodes);
		
		for(edge_list) |*edge| {
			edge.list = std.ArrayList(Edge).init(allocator);
			edge.red_degree = 0;
			edge.deleted = false;
		}
		// Remove allocations on failure
		errdefer {
			for(edge_list) |edge| {
				edge.list.deinit();
			}
			allocator.free(edge_list);
		}

		std.log.info("Loaded graph {s} with bytes len {}", .{filename,buffer.len});
		
		var alive_nodes = try std.ArrayList(u32).initCapacity(allocator,problem.number_of_nodes);
		errdefer {
			alive_nodes.deinit();
		}
	
		var i: u32 = 0;
		while(i < problem.number_of_nodes): (i+=1) {
			try alive_nodes.append(i);
		}

		var graph = Graph {
			.number_of_nodes = problem.number_of_nodes,
			.number_of_edges = problem.number_of_edges,
			.graph = buffer,
			.adj_list = edge_list,
			.alive_nodes = alive_nodes,
			.allocator = allocator,
			.contraction = try ContractionSequence.init(problem.number_of_nodes,allocator)
		};

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

			try graph.addEdge(node_id,second_node_id,EdgeColor.black);
		}

		for(graph.adj_list) |*edge| {
			std.sort.sort(Edge,edge.list.items,@as(usize,0),edgeCmp);
		}
		
		return graph;
	}

	pub fn release(self: Graph) void {
		self.allocator.free(self.graph);
		for(self.adj_list) |edge| {
			edge.list.deinit();
		}
		self.contraction.release();
		self.alive_nodes.deinit();
		self.allocator.free(self.adj_list);
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
