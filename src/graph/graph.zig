const std = @import("std");
const bitset = @import("../util/two_level_bitset.zig");
const edge_list = @import("edge_list.zig");
const red_edge_stack = @import("red_edge_stack.zig");
const contraction = @import("../tww/contraction_sequence.zig");
const comptime_util = @import("../util/comptime_checks.zig");
const RetraceableContractionSequence = @import("../tww/retraceable_contraction_sequence.zig").RetraceableContractionSequence;
const connected_components = @import("connected_component.zig");
const Node = @import("node.zig").Node;
const bfs_mod = @import("bfs.zig");

const pace_2023 = @import("../pace_2023/pace_fmt.zig");

pub const GraphError = error {
	FileNotFound,
	NotPACEFormat,
	GraphTooLarge,
	MisformedEdgeList,
	InvalidContractionOneNodeErased,
	ContractionOverflow,
	NoContractionLeft,
	RetraceNoParent
};

pub fn Graph(comptime T: type) type {
	comptime if (!comptime_util.checkIfIsCompatibleInteger(T)) {
		@compileError("T must either be u8,u16 or u32!");
	};

	return struct {
		const Self = @This();
		number_of_nodes: u32,
		number_of_edges: u32,
		node_list: []Node(T),
		erased_nodes: bitset.FastBitSet,
		scratch_bitset: bitset.FastBitSet,
		connected_components: std.ArrayListUnmanaged(connected_components.ConnectedComponent(T)),
		connected_components_min_heap: std.PriorityQueue(connected_components.ConnectedComponentIndex(T),void,connected_components.ConnectedComponentIndex(T).compareComponentIndexDesc),
		allocator: std.mem.Allocator,
		contraction: RetraceableContractionSequence(T),
	
		pub inline fn addEdge(self: *Self, u: T, v: T) !void {
			self.connected_components.clearRetainingCapacity();
			self.connected_components_min_heap.shrinkAndFree(0);
			std.debug.assert(u < self.node_list.len);
			std.debug.assert(v < self.node_list.len);

			const u_node = &self.node_list[u];
			const v_node = &self.node_list[v];

			_ = try u_node.black_edges.add(self.allocator,v);
			_ = try v_node.black_edges.add(self.allocator,u);

			self.number_of_edges+=1;
		}

	pub inline fn memoryUsage(self: *Self) usize {
		return self.allocator.end_index;
	}

	pub inline fn getCurrentTwinWidth(self: *const Self) u32 {
		return self.contraction.getTwinWidth();
	}

	pub fn reverseLastContraction(self: *Self) !void {
		_ = self;
	}
	pub fn addContraction(self: *Self, erased: T, survivor: T) !T {
		_ = self;
		_ = survivor;
		_ = erased;
		return 0;
	}

	pub fn new(number_of_nodes: T, allocator: std.mem.Allocator) !Self {
		var node_list = try allocator.alloc(Node(T),number_of_nodes);
		
		for(node_list) |*node| {
			node.black_edges = try edge_list.ParametrizedSortedArrayList(T).initCapacity(allocator,2);
			node.red_edges = try edge_list.ParametrizedUnsortedArrayList(T).initCapacity(allocator,2);
		}

		//TODO: Add some errdefer's here

		var graph = Self {
			.number_of_nodes = number_of_nodes,
			.number_of_edges = 0,
			.node_list = node_list,
			.allocator = allocator,
			.scratch_bitset = try bitset.FastBitSet.initEmpty(number_of_nodes,allocator),
			.erased_nodes = try bitset.FastBitSet.initEmpty(number_of_nodes,allocator),
			.connected_components = std.ArrayListUnmanaged(connected_components.ConnectedComponent(T)){},
			.connected_components_min_heap = std.PriorityQueue(connected_components.ConnectedComponentIndex(T),void,connected_components.ConnectedComponentIndex(T).compareComponentIndexDesc).init(allocator,{}),
			.contraction = try RetraceableContractionSequence(T).init(allocator,number_of_nodes)
		};

		return graph;
	}

	pub fn loadFromPace(allocator: std.mem.Allocator, filename: []const u8) !Self {
		const pace = try pace_2023.Pace2023Fmt(T).fromFile(allocator, filename);
		defer pace.deinit(allocator);

		var node_list = try allocator.alloc(Node(T),pace.number_of_nodes);
		
		for(0..pace.number_of_nodes) |index| {
			node_list[index].black_edges = pace.nodes[index].edges.intoSorted();
			node_list[index].red_edges = edge_list.ParametrizedUnsortedArrayList(T).init();
		}

		// Remove allocations on failure
		errdefer {
			for(node_list) |*node| {
				node.black_edges.deinit(allocator);
				node.red_edges.deinit(allocator);
			}
			allocator.free(node_list);
		}


		return Self {
			.number_of_nodes = pace.number_of_nodes,
			.number_of_edges = pace.number_of_edges, //ATM
			.node_list = node_list,
			.allocator = allocator,
			.erased_nodes = try bitset.FastBitSet.initEmpty(pace.number_of_nodes,allocator),
			.scratch_bitset = try bitset.FastBitSet.initEmpty(pace.number_of_nodes,allocator),
			.connected_components = std.ArrayListUnmanaged(connected_components.ConnectedComponent(T)){},
			.connected_components_min_heap = std.PriorityQueue(connected_components.ConnectedComponentIndex(T),void,connected_components.ConnectedComponentIndex(T).compareComponentIndexDesc).init(allocator,{}),
			.contraction = try RetraceableContractionSequence(T).init(allocator,pace.number_of_nodes)
		};
	}

	pub fn findAllConnectedComponents(self: *Self) !void {
		self.scratch_bitset.unsetAll();

		var bfs_stack = try bfs_mod.BfsQueue(T).init(self.allocator, self.number_of_nodes);

		var unsetIter = self.scratch_bitset.iterUnset();
		var components:u32 = 0;
		var total_sum:u32 = 0;
		var largest:u32 = 0;
		var node_list = edge_list.ParametrizedUnsortedArrayList(T).init();

		while(unsetIter.next()) |item| {
			if(self.scratch_bitset.get(item)) {
				continue;
			}
			var iterator = bfs_mod.bfs(T,0,self,&self.scratch_bitset,&bfs_stack,.{.max_level = std.math.maxInt(T), .kind=.black});

			while(iterator.next()) |node| {
				self.scratch_bitset.set(node);
				try node_list.add(self.allocator,node);
			}

			components+=1;
			total_sum+=node_list.cardinality();
			largest = std.math.max(largest,node_list.cardinality());
			try self.connected_components.append(self.allocator,try connected_components.ConnectedComponent(T).init(self.allocator,node_list,iterator.level));
			node_list = edge_list.ParametrizedUnsortedArrayList(T).init();
		}

		try self.connected_components_min_heap.ensureTotalCapacity(self.connected_components.items.len);
		for(0..self.connected_components.items.len) |index| {
			try self.connected_components_min_heap.add(connected_components.ConnectedComponentIndex(T) {
				.tww = self.connected_components.items[index].tww,
				.index = @intCast(T,index)
			});
		}

		var tww = self.connected_components_min_heap.remove();
		std.debug.print("Found {} components largest {} and tww {}\n",.{components, largest, tww});
	}

	pub fn deinit(self: *Self) void {
		for(self.node_list) |*node| {
			node.black_edges.deinit(self.allocator);
			node.red_edges.deinit(self.allocator);
		}
		for(self.connected_components.items) |*connected_component| {
			connected_component.deinit(self.allocator);
		}
		self.connected_components.deinit(self.allocator);
		self.contraction.deinit();
		self.scratch_bitset.deinit();
		self.erased_nodes.deinit();
		self.allocator.free(self.node_list);
	}
	};
}


test "Check simple loading" {
	var allocator = std.heap.GeneralPurposeAllocator(.{}){};
	defer std.debug.assert(!allocator.deinit());

	var graph = try Graph(u8).loadFromPace(allocator.allocator(),"instances/tiny/tiny001.gr");
	// Free ressources
	defer graph.deinit();

	try std.testing.expectEqual(graph.number_of_nodes,10);
	try std.testing.expectEqual(graph.number_of_edges,9);
}

test "Test contraction Tiny 1" {
	var allocator = std.heap.GeneralPurposeAllocator(.{}){};
	defer std.debug.assert(!allocator.deinit());

	var graph = try Graph(u8).loadFromPace(allocator.allocator(),"instances/tiny/tiny001.gr");
	// Free ressources
	defer graph.deinit();

	try std.testing.expectEqual(graph.number_of_nodes,10);
	try std.testing.expectEqual(graph.number_of_edges,9);

	var i: u8 = 9;
	while(i > 0) {
		var tww = try graph.addContraction(i-1,9);
		if(i>1) {
			try std.testing.expectEqual(@as(u32,1),tww);
		}
		else {
			try std.testing.expectEqual(@as(u32,1),tww);
		}
		i-=1;
	}
}

test "Test contraction Tiny 2" {
	var allocator = std.heap.GeneralPurposeAllocator(.{}){};
	defer std.debug.assert(!allocator.deinit());

	var graph = try Graph(u8).loadFromPace(allocator.allocator(),"instances/tiny/tiny002.gr");
	// Free ressources
	defer graph.deinit();

	try std.testing.expectEqual(graph.number_of_nodes,10);
	try std.testing.expectEqual(graph.number_of_edges,10);

	var i: u8 = 9;
	while(i > 0) {
		var tww = try graph.addContraction(i-1,9);
		if(i>2) {
			try std.testing.expectEqual(@as(u32,2),tww);
		}
		else if(i>1) {
			try std.testing.expectEqual(@as(u32,2),tww);
		}
		else {
			try std.testing.expectEqual(@as(u32,2),tww);
		}
		i-=1;
	}
}

test "Test contraction retrace Tiny 2" {
	var allocator = std.heap.GeneralPurposeAllocator(.{}){};
	defer std.debug.assert(!allocator.deinit());

	var graph = try Graph(u8).loadFromPace(allocator.allocator(),"instances/tiny/tiny002.gr");
	// Free ressources
	defer graph.deinit();

	try std.testing.expectEqual(graph.number_of_nodes,10);
	try std.testing.expectEqual(graph.number_of_edges,10);

	var i: u8 = 9;
	while(i > 0) {
		var tww = try graph.addContraction(i-1,9);
		if(i>2) {
			try std.testing.expectEqual(@as(u32,2),tww);
		}
		else if(i>1) {
			try std.testing.expectEqual(@as(u32,2),tww);
		}
		else {
			try std.testing.expectEqual(@as(u32,2),tww);
		}
		i-=1;
	}
	i = 9;
	while(i > 0) {
		try graph.reverseLastContraction();

		const tww = graph.getCurrentTwinWidth();

		if(i==9) {
			try std.testing.expectEqual(@as(u32,2),tww);
		}
		else if(i==1) {
			try std.testing.expectEqual(@as(u32,0),tww);
		}
		else {
			try std.testing.expectEqual(@as(u32,2),tww);
		}
		for(graph.node_list) |*node| {
			try std.testing.expect(node.red_edges.cardinality() <= tww);
		}
		i-=1;
	}
}

test "Test contraction Tiny 3" {
	var allocator = std.heap.GeneralPurposeAllocator(.{}){};
	defer std.debug.assert(!allocator.deinit());

	var graph = try Graph(u8).loadFromPace(allocator.allocator(),"instances/tiny/tiny003.gr");
	// Free ressources
	defer graph.deinit();

	try std.testing.expectEqual(graph.number_of_nodes,10);
	try std.testing.expectEqual(graph.number_of_edges,45);

	var i: u8 = 9;
	while(i > 0) {
		var tww = try graph.addContraction(i-1,9);
		try std.testing.expectEqual(@as(u32,0),tww);
		i-=1;
	}
}

test "Test contraction retrace Tiny 3" {
	var allocator = std.heap.GeneralPurposeAllocator(.{}){};
	defer std.debug.assert(!allocator.deinit());

	var graph = try Graph(u8).loadFromPace(allocator.allocator(),"instances/tiny/tiny003.gr");
	// Free ressources
	defer graph.deinit();

	try std.testing.expectEqual(graph.number_of_nodes,10);
	try std.testing.expectEqual(graph.number_of_edges,45);

	var i: u8 = 9;
	while(i > 0) {
		var tww = try graph.addContraction(i-1,9);
		try std.testing.expectEqual(@as(u32,0),tww);
		i-=1;
	}
	i = 9;
	while(i > 0) {
		try graph.reverseLastContraction();
		var tww = graph.getCurrentTwinWidth();
		try std.testing.expectEqual(@as(u32,0),tww);
		for(graph.node_list) |*node| {
			try std.testing.expect(node.red_edges.cardinality() <= tww);
		}
		i-=1;
	}
}

test "Test contraction Tiny 4" {
	var allocator = std.heap.GeneralPurposeAllocator(.{}){};
	defer std.debug.assert(!allocator.deinit());

	var graph = try Graph(u8).loadFromPace(allocator.allocator(),"instances/tiny/tiny004.gr");
	// Free ressources
	defer graph.deinit();

	try std.testing.expectEqual(graph.number_of_nodes,10);
	try std.testing.expectEqual(graph.number_of_edges,9);

	var i: u8 = 9;
	while(i > 0) {
		var tww = try graph.addContraction(i-1,9);
		try std.testing.expectEqual(@as(u32,0),tww);
		i-=1;
	}
}

test "Test contraction retrace Tiny 4" {
	var allocator = std.heap.GeneralPurposeAllocator(.{}){};
	defer std.debug.assert(!allocator.deinit());

	var graph = try Graph(u8).loadFromPace(allocator.allocator(),"instances/tiny/tiny004.gr");
	// Free ressources
	defer graph.deinit();

	try std.testing.expectEqual(graph.number_of_nodes,10);
	try std.testing.expectEqual(graph.number_of_edges,9);

	var i: u8 = 9;
	while(i > 0) {
		var tww = try graph.addContraction(i-1,9);
		try std.testing.expectEqual(@as(u32,0),tww);
		i-=1;
	}
	i = 9;
	while(i > 0) {
		try graph.reverseLastContraction();
		var tww = graph.getCurrentTwinWidth();
		try std.testing.expectEqual(@as(u32,0),tww);
		for(graph.node_list) |*node| {
			try std.testing.expect(node.red_edges.cardinality() <= tww);
		}
		i-=1;
	}
}

test "Test contraction Tiny 6" {
	var allocator = std.heap.GeneralPurposeAllocator(.{}){};
	defer std.debug.assert(!allocator.deinit());

	var graph = try Graph(u8).loadFromPace(allocator.allocator(),"instances/tiny/tiny006.gr");
	// Free ressources
	defer graph.deinit();

	try std.testing.expectEqual(graph.number_of_nodes,10);
	try std.testing.expectEqual(graph.number_of_edges,5);

	var i: u8 = 9;
	var tww:u32 = 0;
	while(true) {
		tww = std.math.max(tww,try graph.addContraction(i-1,i));
		if(i==1) {
			break;
		}
		i-=2;
	}
	try std.testing.expectEqual(@as(u32,0),tww);
}

test "Check simple failed loading" {

	var allocator = std.heap.GeneralPurposeAllocator(.{}){};
	// Should never fail since the ressources should be managed
	defer std.debug.assert(!allocator.deinit());

	var graph = Graph(u8).loadFromPace(allocator.allocator(),"instances/tiny/tiny001.graph") catch |err| {
		// Loading should fail with file not found
		try std.testing.expectEqual(err,error.FileNotFound);
		return;
	};
	_ = graph;
	// Should never be reached
	std.debug.assert(false);
}

test "Check simple graph" {

	var allocator = std.heap.GeneralPurposeAllocator(.{}){};
	// Should never fail since the ressources should be managed
	//defer std.debug.assert(!allocator.deinit());
	
	var graph = try Graph(u8).new(5,allocator.allocator());

	try graph.addEdge(0,1);
	try graph.addEdge(1,2);
	try graph.addEdge(1,3);
	try graph.addEdge(1,4);

	_ = try graph.addContraction(1,0);

	try std.testing.expectEqual(graph.node_list[0].red_edges.cardinality(),3);
}
