const std = @import("std");
const comptime_util = @import("../util/comptime_checks.zig");
const edge_list = @import("../graph/edge_list.zig");
const two_level_bitset = @import("../util/two_level_bitset.zig");

pub const CompressedBitmapBucketTag = enum(u8) {
	bitset,
	array
};


pub fn FastCompressedBitmap(comptime T: type, comptime promote_threshold: u32, comptime degrade_threshold: u32) type {
	comptime if(!comptime_util.checkIfIsCompatibleInteger(T)) {
		@compileError("Type must either be u8,u16 or u32");
	};

	comptime if(promote_threshold < degrade_threshold) {
		@compileError("Promote threshold must be larger than degrade threshold!");
	};
	

	return struct {
		const Self = @This();
		storage: StorageRepresentation,
		max_size: T,

		pub const StorageRepresentation = union(CompressedBitmapBucketTag) {
			bitset: two_level_bitset.FastBitSet,
			array: edge_list.ParametrizedSortedArrayList(T),

			pub inline fn deinit(self: *StorageRepresentation, allocator: std.mem.Allocator) void {
				switch(self.*) {
					.bitset => self.bitset.deinit(allocator),
					.array => self.array.deinit(allocator),
				}
			}


		pub inline fn add(self: *StorageRepresentation, allocator: std.mem.Allocator, item: T) !void {
			switch(self.*) {
				.bitset => {
					self.bitset.set(item);
				},
				.array => {
					_ = try self.array.add(allocator,item);
				}
			}
		}

		pub inline fn addExists(self: *StorageRepresentation, allocator: std.mem.Allocator, item: T) !bool {
			switch(self.*) {
				.bitset => {
					return self.bitset.setExists(item);
				},
				.array => {
					return try self.array.add(allocator,item);
				}
			}
		}


		pub inline fn contains(self: *StorageRepresentation, item: T) bool {
			switch(self.*) {
				.bitset => return self.bitset.get(item),
				.array => return self.array.contains(item),
			}
		}
		pub inline fn remove(self: *StorageRepresentation, item: T) bool {
			switch(self.*) {
				.bitset => return self.bitset.unset(item),
				.array => return self.array.remove(item),
			}
		}

		pub inline fn cardinality(self: *StorageRepresentation) T {
			switch(self.*) {
				.bitset => return @intCast(T,self.bitset.cardinality),
					.array => return @intCast(T,self.array.cardinality()),
			}
		}
		};

		pub fn init(max_size: T) Self {
			return Self {
				.storage = StorageRepresentation {.array = edge_list.ParametrizedSortedArrayList(T).init()},
				.max_size = max_size
			};
		}

		pub fn fromUnsorted(allocator: std.mem.Allocator, list: *edge_list.ParametrizedUnsortedArrayList(T),max_size: T) !Self {
			if(list.cardinality() >= promote_threshold) {
				var bitset = try two_level_bitset.FastBitSet.initEmpty(max_size,allocator);

				var iter = list.iterator();
				while(iter.next()) |item| {
					bitset.set(item);
				}
				const container = Self {
						.storage = StorageRepresentation {.bitset = bitset},
						.max_size = max_size
				};

				return container;

			}
			else {
				return Self {
						.storage = StorageRepresentation {.array = list.intoSorted()},
						.max_size = max_size
				};
			}
		}

		pub fn deinit(self: *Self, allocator: std.mem.Allocator) void {
			self.storage.deinit(allocator);
		}

		pub const FastCompressedBitmapIterator = union(CompressedBitmapBucketTag) {
			bitset: two_level_bitset.FastBitSetIterator,
			array: edge_list.ParametrizedSortedArrayList(T).ParametrizedSortedArrayListIterator,
			pub inline fn next(self: *FastCompressedBitmapIterator) ?T {
				switch(self.*) {
					.bitset => {
						if(self.bitset.next()) |result| {
							return @intCast(T,result);
						}
						return null;
					},
					.array => return self.array.next(),
				}
			}
		};

		pub inline fn iterator(self: *const Self) FastCompressedBitmapIterator {
			switch(self.storage) {
				.bitset => return FastCompressedBitmapIterator {.bitset = self.storage.bitset.iter()},
				.array => return FastCompressedBitmapIterator {.array = self.storage.array.iterator()},
			}
		}

		pub inline fn xorIterator(self: *const Self, other: *const Self) FastCompressedBitmapXorIterator {
			var iterator_first = self.iterator();
			var iterator_second = other.iterator();
			const item_first = iterator_first.next();
			const item_second = iterator_second.next();
			return FastCompressedBitmapXorIterator {
				.iterator_first = iterator_first,
				.iterator_second = iterator_second,
				.item_first = item_first,
				.item_second = item_second,
				.first = false
			};
		}

		pub const FastCompressedBitmapXorIterator = struct {
			iterator_first: FastCompressedBitmapIterator,
			iterator_second: FastCompressedBitmapIterator,
			item_first: ?T,
			item_second: ?T,
			first: bool,
			pub inline fn next(self: *FastCompressedBitmapXorIterator) ?T {
				while(self.item_first!=null and self.item_second!=null) {
					if(self.item_first.? < self.item_second.?) {
						const return_item_first = self.item_first;
						self.item_first= self.iterator_first.next();
						self.first = true;
						return return_item_first;
					}
					else if(self.item_first.? > self.item_second.?) {
						const return_item_second = self.item_second;
						self.item_second = self.iterator_second.next();
						self.first = false;
						return return_item_second;
					}

					self.item_first= self.iterator_first.next();
					self.item_second = self.iterator_second.next();
				}

				if(self.item_second != null) {
					const ret = self.item_second;
					self.item_second = self.iterator_second.next();
					self.first = false;
					return ret;
				}
				else if(self.item_first != null) {
					const ret = self.item_first;
					self.item_first= self.iterator_first.next();
					self.first = true;
					return ret;
				}
				return null;
			}
		};

		pub inline fn forAll(self: *const Self, comptime ctx: type, comptime fnc: fn (ctx,T) callconv(.Inline) void, context: ctx) void {
			switch(self.storage) {
				.bitset => {
					var iter = self.storage.bitset.iter();
					while(iter.next()) |result| {
						fnc(context,result);
					}
				},
				.array => {
					var iter = self.storage.array.iterator();
					while(iter.next()) |result| {
						fnc(context,result);
					}
				},
			}
		}

		pub inline fn add(self: *Self, allocator: std.mem.Allocator, item: T) !void {
			try self.storage.add(allocator,item);
			try self.checkPromote(allocator);
		}

		inline fn checkPromote(self: *Self, allocator: std.mem.Allocator) !void {
			if(self.storage.cardinality() >= promote_threshold and self.storage == .array) {
				var iter = self.storage.array.iterator();
				var bitset = try two_level_bitset.FastBitSet.initEmpty(self.max_size,allocator);
				while(iter.next()) |value| {
					bitset.set(value);
				}
				self.storage.array.deinit(allocator);

				self.storage = StorageRepresentation {.bitset = bitset};
			}
		}

		pub inline fn addExists(self: *Self, allocator: std.mem.Allocator, item: T) !bool {
			const result = try self.storage.addExists(allocator,item);
			try self.checkPromote(allocator);
			return result;
		}

		pub inline fn contains(self: *Self, item: T) bool {
			return self.storage.contains(item);
		}

		pub inline fn cardinality(self: *Self) T {
			return self.storage.cardinality();
		}

		pub inline fn remove(self: *Self, allocator: std.mem.Allocator, item: T) !bool {
			if(self.storage.remove(item)) {
				// Degrade storage to array 
				if(degrade_threshold > self.cardinality() and self.storage == .bitset) {
					var bitset = self.storage.bitset.iter();
					var array = try edge_list.ParametrizedSortedArrayList(T).initCapacity(allocator,self.storage.bitset.cardinality);
					while(bitset.next()) |value| {
						array.edges.append(allocator,@intCast(T,value)) catch unreachable;
					}
					self.storage = StorageRepresentation {.array = array};
				}
				return true;
			}
			return false;
		}
	};
}

test "CompressedBitmap: Check get" {
	var gpa = std.heap.GeneralPurposeAllocator(.{}){};
	defer std.debug.assert(!gpa.deinit());



	var empty = FastCompressedBitmap(u16,4,6).init(2000);
	defer empty.deinit(gpa.allocator());	

	var i:u16 = 0;
	while(i<2000) : (i+=1) {
		std.debug.assert(empty.contains(i) == false);
	}
	try std.testing.expectEqual(empty.cardinality(),0);
}

test "CompressedBitmap: Check set" {
	var gpa = std.heap.GeneralPurposeAllocator(.{}){};
	defer std.debug.assert(!gpa.deinit());

	var empty = FastCompressedBitmap(u16,4,6).init(2000);
	defer empty.deinit(gpa.allocator());	

	var i:u16 = 0;
	while(i<2000) : (i+=1) {
		try empty.add(gpa.allocator(),i);
	}
	i = 0;
	while(i<2000) : (i+=1) {
		std.debug.assert(empty.contains(i) == true);
	}
	try std.testing.expectEqual(empty.cardinality(),2000);
}

test "CompressedBitmap: Check large" {
	var gpa = std.heap.GeneralPurposeAllocator(.{
		.enable_memory_limit = true
	}){};
	gpa.requested_memory_limit = 4000;
	defer std.debug.assert(!gpa.deinit());

	var empty = FastCompressedBitmap(u32,9,11).init(3_000_000);
	defer empty.deinit(gpa.allocator());	

	var i:u16 = 0;
	while(i<20) : (i+=1) {
		try empty.add(gpa.allocator(),i);
	}
	i = 0;
	while(i<20) : (i+=1) {
		std.debug.assert(empty.contains(i) == true);
	}
	_ = try empty.remove(gpa.allocator(),0);
	i = 0;
	while(i<20) : (i+=1) {
		if(i==0) {
			std.debug.assert(empty.contains(i) == false);
		}
		else {
			std.debug.assert(empty.contains(i) == true);
		}
	}
	try std.testing.expectEqual(empty.cardinality(),19);
	std.debug.print("Bytes {}\n",.{gpa.total_requested_bytes});
}

test "CompressedBitmap: Check large iterator" {
	var gpa = std.heap.GeneralPurposeAllocator(.{
		.enable_memory_limit = true
	}){};
	gpa.requested_memory_limit = 4000;
	defer std.debug.assert(!gpa.deinit());

	var empty = FastCompressedBitmap(u32,9,11).init(3_000_000);
	defer empty.deinit(gpa.allocator());	

	var i:u16 = 0;
	while(i<20) : (i+=2) {
		try empty.add(gpa.allocator(),i);
	}
	var iterator = empty.iterator();
	i=0;
	while(iterator.next()) |item| {
		std.debug.assert(item == i);
		i+=2;
	}

	try std.testing.expectEqual(empty.cardinality(),10);
	std.debug.print("Bytes {}\n",.{gpa.total_requested_bytes});
}

pub inline fn checkInline(ctx: *u64, result: u32) void {
	ctx.*+=result;
}

test "bench: CompressedBitmap: Check large iterator" {
	var gpa = std.heap.GeneralPurposeAllocator(.{
	}){};
	defer std.debug.assert(!gpa.deinit());

	var empty = FastCompressedBitmap(u32,10,12).init(10_000_000);
	defer empty.deinit(gpa.allocator());	

	var i:u32 = 0;
	while(i<10_000_000) : (i+=100) {
		try empty.add(gpa.allocator(),i);
	}
	var iterator = empty.iterator();

	var count:u64 = 0;
	while(iterator.next()) |item| {
		count+=item;
	}

	try std.testing.expectEqual(count,499995000000);
}

test "bench: CompressedBitmap: Check large iterator inline" {
	var gpa = std.heap.GeneralPurposeAllocator(.{
	}){};
	defer std.debug.assert(!gpa.deinit());

	var empty = FastCompressedBitmap(u32,10,12).init(10_000_000);
	defer empty.deinit(gpa.allocator());	

	var i:u32 = 0;
	while(i<10_000_000) : (i+=100) {
		try empty.add(gpa.allocator(),i);
	}

	var count:u64 = 0;
	empty.forAll(*u64,checkInline,&count);

	try std.testing.expectEqual(count,499995000000);
}
