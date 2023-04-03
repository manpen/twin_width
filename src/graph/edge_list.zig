const std = @import("std");
const bitset = @import("../util/two_level_bitset.zig");

pub const EdgeListType = enum {
	bitmap,
	list
};

pub const BitmapEdgeList = struct {
	edges: bitset.FastBitSet
};

pub const ListEdgeList = struct {
	const Self = @This();
	edges: std.ArrayList(u32),
	trailing_edges: bool,


	pub const ListEdgeListIterator = struct {
		index: u32,
		edge_list: *const ListEdgeList,
		pub fn next(self: *ListEdgeListIterator) ?u32 {
			while(true) {
				if(self.index>=self.edge_list.edges.items.len) {
					return null;
				}
				const item = self.edge_list.edges.items[self.index];
				self.index+=1;
				if((item&0x80000000) != 0) {
					continue;
				}
				return item;
			}
		}
	};

	pub inline fn compactList(self: *Self) void {
		if(self.trailing_edges) {
			var i:u32 = 0;
			var writePtr:u32 = 0;
			while(i < self.edges.items.len) : (i+=1) {
				self.edges.items[writePtr] = self.edges.items[i];
				if((self.edges.items[writePtr]&0x80000000)==0) {
					writePtr+=1;
				}
			}
			self.trailing_edges = false;
			self.edges.shrinkRetainingCapacity(writePtr);
		}
	}

	pub inline fn cardinality(self: *Self) u32 {
		self.compactList();
		return self.edges.items.len;
	}

	pub inline fn iterator(self: *const Self) ListEdgeListIterator {
		return ListEdgeListIterator {
			.index = 0,
			.edge_list = self
		};
	}

	pub inline fn markForRemove(self: *Self, node_id: u32) void {
		var i:u32 = 0;
		while(i < self.edges.items.len) : (i+=1) {
			if(self.edges.items[i] == node_id) {
				self.edges.items[i] |= 0x80000000;
				self.trailing_edges = true;
			}
			else if(self.edges.items[i] > node_id) {
				return;
			}
		}
	}

	pub inline fn add(self: *Self, node_id: u32) !bool {
		var i:u32 = 0;
		while(i < self.edges.items.len) {
			if(self.edges.items[i] == node_id) {
				return false;
			}
			else if(self.edges.items[i] > node_id) {
				try self.edges.insert(i,node_id);
				return true;
			}
			i+=1;
		}
	}
};

pub fn EdgeListXorIterator(comptime X: type, comptime Y: type) type {
	return struct {
		const currentType = @This();
iterator: X,
						iterator_2: Y,
						item: ?u32,
						item_2: ?u32,
						pub inline fn next(self: *currentType) ?u32 {
							while(self.item!=null and self.item_2!=null) {
								if(self.item.? < self.item_2.?) {
									const return_item = self.item;
									self.item = self.iterator.next();
									return return_item;
								}
								else if(self.item.? > self.item_2.?) {
									const return_item = self.item_2;
									self.item_2 = self.iterator_2.next();
									return return_item;					
								}

								self.item = self.iterator.next();
								self.item_2 = self.iterator_2.next();
							}
							if(self.item != null) {
								const ret = self.item;
								self.item = null;
								return ret;
							}
							if(self.item_2 != null) {
								const ret = self.item_2;
								self.item_2 = null;
								return ret;
							}

							while(self.iterator.next()) |item| {
								return item;
							}
							while(self.iterator_2.next()) |item| {
								return item;
							}
							return null;
						}
	};
}

pub const XorIteratorTag = enum {
	bitmap_bitmap,
	list_bitmap,
	list_list
};

pub const XorIterator = union(XorIteratorTag) {
	const Self = @This();
	bitmap_bitmap: EdgeListXorIterator(bitset.FastBitSetIterator, bitset.FastBitSetIterator),
	list_bitmap: EdgeListXorIterator(ListEdgeList.ListEdgeListIterator, bitset.FastBitSetIterator),
	list_list: EdgeListXorIterator(ListEdgeList.ListEdgeListIterator, ListEdgeList.ListEdgeListIterator),

	pub inline fn next(self: *Self) ?u32 {
		switch(self.*) {
			.bitmap_bitmap => return self.bitmap_bitmap.next(),
			.list_bitmap => return self.list_bitmap.next(),
			.list_list => return self.list_list.next(),
		}
	}
};

pub const EdgeListIterator = union(EdgeListType) {
	const Self = @This();
	bitmap: bitset.FastBitSetIterator,
	list: ListEdgeList.ListEdgeListIterator,

	pub inline fn next(self: *Self) ?u32 {
		switch(self.*) {
			.bitmap => return self.bitmap.next(),
			.list => return self.list.next()
		}
	}
};

pub const EdgeList = union(EdgeListType) {
	const Self = @This();
	bitmap: BitmapEdgeList,
	list: ListEdgeList,
	
	var thresholdPromote:u32 = 0;
	var thresholdDegrade:u32 = 0;
	var problemSize:u32 = std.math.maxInt(u32);

	pub fn initCapacity(capacity: u32, allocator: std.mem.Allocator) EdgeList {
		std.debug.assert(EdgeList.problemSize != std.math.maxInt(u32));

		if(capacity >= EdgeList.thresholdPromote) {
			return EdgeList {
				.bitmap = BitmapEdgeList {
					.edges = bitset.FastBitSet.initCapacity(EdgeList.problemSize, allocator)
				}
			};
		}
		else {
			return EdgeList {
				.list = ListEdgeList {
					.edges = std.ArrayList(u32).initCapacity(allocator,capacity)
				}
			};
		}
	}

	pub fn deinit(self: *Self) void {
		switch(self.*) {
			.bitmap => self.bitmap.edges.deinit(),
			.list => self.list.edges.deinit()
		}
	}

	pub inline fn cardinality(self: *Self) u32 {
		switch(self.*) {
			.bitmap => return self.bitmap.cardinality(),
			.list => return self.list.cardinality()
		}
	}

	pub inline fn iterator(self: *Self) EdgeListIterator {
		switch(self.*) {
			.bitmap => return EdgeListIterator {
				.bitmap = self.bitmap.edges.iterator()
			},
			.list => return EdgeListIterator {
				.list = self.list.iterator()
			}
		}
	}

	pub inline fn xorIterator(self: *Self, other: *Self) XorIterator {
		switch(self.*) {
			.bitmap => switch(other.*) {
				.bitmap => {
					var iter_1 = self.bitmap.edges.iterator();
					var iter_2 = other.bitmap.edges.iterator();
					const item_1 = iter_1.next();
					const item_2 = iter_2.next();
					return XorIterator {
						.bitmap_bitmap = EdgeListXorIterator(bitset.FastBitSetIterator,bitset.FastBitSetIterator) {
							.iterator = iter_1,
							.iterator_2 = iter_2,
							.item = item_1,
							.item_2 = item_2
						}
					};
				},
				.list => {
					var iter_1 = self.bitmap.edges.iterator();
					var iter_2 = other.list.iterator();
					const item_1 = iter_1.next();
					const item_2 = iter_2.next();
					return XorIterator {
						.list_bitmap = EdgeListXorIterator(ListEdgeList.ListEdgeListIterator,bitset.FastBitSetIterator) {
							.iterator = iter_1,
							.iterator_2 = iter_2,
							.item = item_1,
							.item_2 = item_2
						}
					};

				},
			},
			.list => switch(other.*) {
				.list => {
					var iter_1 = self.list.iterator();
					var iter_2 = other.list.iterator();
					const item_1 = iter_1.next();
					const item_2 = iter_2.next();
					return XorIterator {
						.list_bitmap = EdgeListXorIterator(ListEdgeList.ListEdgeListIterator,ListEdgeList.ListEdgeListIterator) {
							.iterator = iter_1,
							.iterator_2 = iter_2,
							.item = item_1,
							.item_2 = item_2
						}
					};
				},
				.bitmap => {
					var iter_1 = self.list.edges.iterator();
					var iter_2 = other.bitmap.iterator();
					const item_1 = iter_1.next();
					const item_2 = iter_2.next();
					return XorIterator {
						.list_bitmap = EdgeListXorIterator(ListEdgeList.ListEdgeListIterator,bitset.FastBitSetIterator) {
							.iterator = iter_1,
							.iterator_2 = iter_2,
							.item = item_1,
							.item_2 = item_2
						}
					};
				}
			}
		}
	}


	pub inline fn remove(self: *Self, node_id: u32) void {
		switch(self.*) {
			.bitmap => self.bitmap.edges.unset(node_id),
			.list => self.list.edges.markForRemove(node_id),
		}
	}
	
	pub inline fn add(self: *Self, node_id: u32) !void {
		switch(self.*) {
			.bitmap => self.bitmap.edges.set(node_id),
			.list => {
				self.list.edges.add(node_id);
				// Promote to bitmap if becoming quite dense
				if(self.list.edges.items.len >= EdgeList.thresholdPromote) {
					var bitmap = BitmapEdgeList {
						.edges = bitset.FastBitSet.initCapacity(EdgeList.problemSize, self.list.edges.allocator)
					};

					var iter = self.list.edges.iterator();
					while(iter.next()) |item| {
						bitmap.set(item);
					}
					self.list.edges.deinit();
					self.* = EdgeList {
						.bitmap = bitmap
					};
				}
			},
		}
	}
};
