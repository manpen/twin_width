const std = @import("std");
const comptime_util = @import("../util/comptime_checks.zig");
const compressed_bitmap = @import("../util/compressed_bitmap.zig");
const two_level_bitmap = @import("../util/two_level_bitset.zig");
const Graph = @import("graph.zig").Graph;

pub const NodeError = error {
	AlreadyPromoted
};


pub fn LargeNodeQueryProcessor(comptime T: type) type {
	comptime if (!comptime_util.checkIfIsCompatibleInteger(T)) {
		@compileError("T must be either u8,u16 or u32");
	};
	return struct {
		const Self = @This();

		pub const NodeDegreeEntry = struct {
			id: T,
			degree: T
		};

		lowest_degree_nodes: std.PriorityQueue(NodeDegreeEntry,void, Self.compareFunction),
		added_nodes: std.bit_set.DynamicBitSetUnmanaged,
		graph: *const Graph(T),
		pub fn compareFunction(ctx: void, lhs: NodeDegreeEntry, rhs: NodeDegreeEntry) std.math.Order {
			_ = ctx;
			return std.math.order(lhs.degree,rhs.degree);
		}

		pub fn init(graph: *const Graph(T)) !Self {
			var bt = std.PriorityQueue(NodeDegreeEntry,void,Self.compareFunction).init(graph.allocator, {});
			var added_nodes = try std.bit_set.DynamicBitSetUnmanaged.initEmpty(graph.allocator,graph.number_of_nodes);
			try bt.ensureTotalCapacity(graph.number_of_nodes);
			return Self {
				.lowest_degree_nodes = bt,
				.graph = graph,
				.added_nodes = added_nodes
			};
		}

		pub fn addNode(self: *Self, node: T) !void {
			if(!self.added_nodes.isSet(node)) {
				self.added_nodes.set(node);
				try self.lowest_degree_nodes.add(NodeDegreeEntry {
					.id = node,
					.degree = self.graph.node_list[node].cardinality()
				});
			}
		}

		pub const TakingLargeNodeIterator = struct {
			index: T,
			items_left: T,
			ptr: *Self,
			pub inline fn next(self: *TakingLargeNodeIterator) ?T {
				while(self.items_left > 0) {
					if(self.index >= self.ptr.lowest_degree_nodes.len) {
						return null;
					}
					const item = self.ptr.lowest_degree_nodes.items[self.index];
					if(self.ptr.graph.erased_nodes.get(item.id)) {
						_ = self.ptr.lowest_degree_nodes.removeIndex(self.index);
						continue;
					}
					else if(self.ptr.graph.node_list[item.id].cardinality() != item.degree) {
						var updated_item = self.ptr.lowest_degree_nodes.removeIndex(self.index);
						self.ptr.added_nodes.unset(updated_item.id);
						self.ptr.addNode(updated_item.id) catch unreachable;
						continue;
					}
					self.index+=1;
					self.items_left-=1;
					return item.id;
				}
				return null;
			}
		};

		pub fn takingIterator(self: *Self, num_items: T) TakingLargeNodeIterator {
			return TakingLargeNodeIterator {
				.ptr = self,
				.index = 0,
				.items_left = num_items
			};
		}
	};
}

// Id is omitted must be stored outside of the struct
pub fn Node(comptime T: type, comptime promote_threshold: u32, comptime degrade_threshold: u32) type {
	comptime if (!comptime_util.checkIfIsCompatibleInteger(T)) {
		@compileError("T must be either u8,u16 or u32");
	};
	
	return struct {
		const Self = @This();
		pub const EdgeIterType = compressed_bitmap.FastCompressedBitmap(T,promote_threshold,degrade_threshold).FastCompressedBitmapIterator;
		black_edges: compressed_bitmap.FastCompressedBitmap(T,promote_threshold,degrade_threshold),
		red_edges: compressed_bitmap.FastCompressedBitmap(T,promote_threshold,degrade_threshold),
		high_degree_node: ?*LargeNodeQueryProcessor(T),
		num_leafes: T,
		largest_red_one_nb: T,

		pub fn promoteToLargeDegreeNode(self: *Self, graph: *const Graph(T)) !void {
			if(self.high_degree_node == null) {
				var created = try graph.allocator.create(LargeNodeQueryProcessor(T));
				created.* = try LargeNodeQueryProcessor(T).init(graph);
				var iter = self.unorderedIterator();
				while(iter.next()) |item| {
					try created.addNode(item);
				}
				self.high_degree_node = created;
			}
			else {
				return error.AlreadyPromoted;
			}
		}

		pub inline fn takeIterator(self: *Self, num_items: T) LargeNodeQueryProcessor(T).TakingLargeNodeIterator {
			return self.high_degree_node.?.takingIterator(num_items);
		}

		pub inline fn isLargeNode(self: *const Self) bool {
			if(self.high_degree_node != null) {
				return true;
			}
			return false;
		}

		pub inline fn addRedEdgeExists(self: *Self, allocator: std.mem.Allocator, node: T) !bool {
			if(self.high_degree_node) |hd| {
				try hd.addNode(node);
			}
			return try self.red_edges.addExists(allocator,node);
		}


		pub inline fn addRedEdge(self: *Self, allocator: std.mem.Allocator, node: T) !void {
			if(self.high_degree_node) |hd| {
				try hd.addNode(node);
			}
			try self.red_edges.add(allocator,node);
		}

		pub inline fn addBlackEdge(self: *Self, allocator: std.mem.Allocator, node: T) !void {
			if(self.high_degree_node) |hd| {
				try hd.addNode(node);
			}
			try self.black_edges.add(allocator,node);
		}

		pub inline fn removeBlackEdge(self: *Self, allocator: std.mem.Allocator, node: T) !void {
			_ = try self.black_edges.remove(allocator, node);
		}

		pub inline fn removeRedEdge(self: *Self, allocator: std.mem.Allocator, node: T) !void {
			_ = try self.red_edges.remove(allocator, node);
		}

		pub inline fn cardinality(self: *const Self) T {
			return self.black_edges.cardinality() + self.red_edges.cardinality();
		}

		pub inline fn getFirstNeighboor(self: *const Self) T {
			if(self.black_edges.cardinality() > 0) {
				var iter = self.black_edges.iterator();
				return iter.next().?;
			}
			else {
				var iter = self.red_edges.iterator();
				return iter.next().?;
			}
		}

		pub inline fn isLeaf(self: *const Self) bool {
			return self.cardinality() == 1;
		}

		pub inline fn unorderedIterator(self: *const Self) UnorderedNodeEdgeIterator {
			var iterator_first = self.red_edges.iterator();
			var iterator_second = self.black_edges.iterator();
			return UnorderedNodeEdgeIterator {
				.red_iter = iterator_first,
				.black_iter = iterator_second,
				.red = false
			};
		}

		pub inline fn orderedIterator(self: *const Self) OrderedNodeEdgeIterator {
			var iterator_first = self.red_edges.iterator();
			var iterator_second = self.black_edges.iterator();
			const item_first = iterator_first.next();
			const item_second = iterator_second.next();
			return OrderedNodeEdgeIterator {
				.red_iter = iterator_first,
				.black_iter = iterator_second,
				.red_next = item_first,
				.black_next = item_second,
				.red = false
			};
		}

		pub const UnorderedNodeEdgeIterator = struct {
			red_iter: compressed_bitmap.FastCompressedBitmap(T,promote_threshold,degrade_threshold).FastCompressedBitmapIterator,
			black_iter: compressed_bitmap.FastCompressedBitmap(T,promote_threshold,degrade_threshold).FastCompressedBitmapIterator,
			red: bool,
			pub inline fn next(self: *UnorderedNodeEdgeIterator) ?T {
				if(self.black_iter.next()) |item| {
					return item;
				}
				self.red = true;
				if(self.red_iter.next()) |item| {
					return item;
				}
				return null;
			}
		};

		pub const OrderedNodeEdgeIterator = struct {
			red_iter: compressed_bitmap.FastCompressedBitmap(T,promote_threshold,degrade_threshold).FastCompressedBitmapIterator,
			black_iter: compressed_bitmap.FastCompressedBitmap(T,promote_threshold,degrade_threshold).FastCompressedBitmapIterator,
			red_next: ?T,
			black_next: ?T,
			red: bool,
			pub inline fn next(self: *OrderedNodeEdgeIterator) ?T {
				while(self.red_next!=null and self.black_next!=null) {
					if(self.red_next.? < self.black_next.?) {
						const return_item_first = self.red_next;
						self.red_next = self.red_iter.next();
						self.red = true;
						return return_item_first;
					}
					else if(self.red_next.? > self.black_next.?) {
						const return_item_second = self.black_next;
						self.black_next = self.black_iter.next();
						self.red = false;
						return return_item_second;
					}
					unreachable;
				}

				if(self.black_next != null) {
					const ret = self.black_next;
					self.black_next = self.black_iter.next();
					self.red = false;
					return ret;
				}
				else if(self.red_next != null) {
					const ret = self.red_next;
					self.red_next = self.red_iter.next();
					self.red = true;
					return ret;
				}
				return null;
			}
		};
	};
}
