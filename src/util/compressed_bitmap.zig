const std = @import("std");


const MAX_SPARSE_BUCKET_SIZE = 100;

pub const CompressedBitmapBucketTag = enum(u8) {
	bitmap,
	array,
	empty
};

pub const BucketBitmap = struct {
	const Self = @This();
	storage: []u64,
	cardinality: u32,
	pub const BucketBitmapIterator = struct {
		const Self = @This();
		coerced_ptr: [*]u64,
		length: usize,
		current: u64,
		index: u32,
		pub inline fn next(self: *BucketBitmapIterator) ?u32 {
			while(self.current == 0) {
				self.index+=1;
				if(self.index>=self.length) {
					return null;
				}
				self.current = self.coerced_ptr[self.index];
			}
			const next_higher_index = @ctz(self.current) + (self.index<<6);
			self.current&= self.current - 1;
			return next_higher_index;
		}
	};

	pub inline fn initEmpty(allocator: std.mem.Allocator) !BucketBitmap {
		var mem  = try allocator.alloc(u64,(65565>>6));
		std.mem.set(u64,mem,0);
		return BucketBitmap {
			.storage = mem,
			.cardinality = 0
		};
	}

	pub inline fn get(self: *const Self, index: u32) bool {
		const u64_index = index>>6;
		const bit_index:u6 = @intCast(u6,(index&0x3F));
		return (self.storage[u64_index] & (@as(u64,1)<<bit_index)) != 0;
	}

	pub inline fn iterator(self: *const Self) BucketBitmapIterator {
		return BucketBitmapIterator {
			.coerced_ptr = @ptrCast([*]u64,self.storage.ptr),
			.length = self.storage.len,
			.current = self.storage[0],
			.index = 0,
		};
	}

	pub inline fn set(self: *Self, index: u32) bool {
		const u64_index = index>>6;
		const bit_index:u6 = @intCast(u6,(index&0x3F));
		const before = self.storage[u64_index];
		self.storage[u64_index] |= (@as(u64,1)<<bit_index);
		if((before&(@as(u64,1)<<bit_index))>0) {
			return false;
		}
		self.cardinality+=1;
		return true;
	}
	
	pub inline fn unset(self: *Self, index: u32) bool {
		const u64_index = index>>6;
		const bit_index:u6 = @intCast(u6,(index&0x3F));
		const before = self.storage[u64_index];
		self.storage[u64_index] &= ~(@as(u64,1)<<bit_index);
		if((before&((@as(u64,1)<<bit_index))>0)) {
			self.cardinality-=1;
			return true;
		}
		return false;

	}
};

pub const BucketArray = struct {
	const Self = @This();
	storage: std.ArrayList(u16),
	pub const BucketArrayIterator = struct {
		current: u32,
		bucket: *const BucketArray,
		pub inline fn next(self: *BucketArrayIterator) ?u32 {
			if(self.current>=self.bucket.cardinality()) {
				return null;
			}
			const item = self.bucket.storage.items[self.current];
			self.current+=1;
			return @as(u32,item);
		}
	};
	pub inline fn get(self: *const Self, index: u32) bool {
		for(self.storage.items) |item| {
			if(item==index) {
				return true;
			}
			else if(item>index) {
				return false;
			}
		}
		return false;
	}

	pub inline fn set(self: *Self, index: u32) !bool {
		var i:u32 = 0;
		while(i<self.storage.items.len): (i+=1) {
			if(self.storage.items[i] > index)	{
				try self.storage.insert(i,@intCast(u16,index));
				return true;
			}
			else if(self.storage.items[i]==index) {
				return false;
			}
		}
		try self.storage.append(@intCast(u16,index));
		return true;
	}

	pub inline fn unset(self: *Self, index: u32) bool {
		var i:u32 = 0;
		while(i<self.storage.items.len): (i+=1) {
			if(self.storage.items[i] == index)	{
				_ = self.storage.orderedRemove(i);
				return true;
			}
			else if(self.storage.items[i] > index) {
				return false;
			}
		}
		return false;
	}

	pub inline fn iterator(self: *const Self) BucketArrayIterator {
		return BucketArrayIterator {
			.current = 0,
			.bucket = self
		};
	}

	pub inline fn cardinality(self: *const Self) u32 {
		return @intCast(u32,self.storage.items.len);
	}

};


pub const CompressedBucket = union(CompressedBitmapBucketTag) {
	bitmap: BucketBitmap,
	array: BucketArray,
	empty: void
};

pub const BucketEmptyIterator = struct {
	pub inline fn next(self: *const BucketEmptyIterator) ?u32 {
		_ = self;
		return null;
	}
};

pub const CompressedBucketIterator = union(CompressedBitmapBucketTag) {
	const Self = @This();
	bitmap: BucketBitmap.BucketBitmapIterator,
	array: BucketArray.BucketArrayIterator,
	empty: BucketEmptyIterator,
	pub inline fn next(self: *Self) ?u32 {
		return switch(self.*) {
			.bitmap => |*bitmap| bitmap.next(),
			.array => |*array| array.next(),
			.empty => |empty| empty.next(),
		};
	}
};

pub const CompressedBitmap = struct {
		const Self = @This();
		buckets: []CompressedBucket,
		allocator: std.mem.Allocator,
		cardinality: u32,

		pub const CompressedBitmapIterator = struct {
			bucket: u32,
			bitmap: *const CompressedBitmap,
			iterator: CompressedBucketIterator,
			pub inline fn next(self: *CompressedBitmapIterator) ?u32 {
				while(true) {
					if(self.iterator.next()) |item| {
						return item+(self.bucket<<16);
					}
					else {
						self.bucket+=1;
						if(self.bucket >= self.bitmap.buckets.len) {
							return null;
						}
						self.iterator = switch(self.bitmap.buckets[self.bucket]) {
							.bitmap => |bucket_b| CompressedBucketIterator {
								.bitmap = bucket_b.iterator()
							},
							.array => |bucket_a| CompressedBucketIterator {
								.array = bucket_a.iterator(),
							},
							.empty => {
								continue;
							}
						};
					}
				}
			}
		};
		pub const CompressedBitmapXORIterator = struct {
			iterator: CompressedBitmapIterator,
			iterator_2: CompressedBitmapIterator,
			item: ?u32,
			item_2: ?u32,
			pub inline fn next(self: *CompressedBitmapXORIterator) ?u32 {
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

		pub inline fn xor_iterator(self: *const Self, other: *const Self) CompressedBitmapXORIterator {
			var iter = self.iterator();
			const item = iter.next();
			var iter2 = other.iterator();
			var item_2 = iter2.next();
			return CompressedBitmapXORIterator {
				.item = item,
				.iterator = iter,
				.item_2 = item_2,
				.iterator_2 = iter2
			};
		}

		pub inline fn iterator(self: *const Self) CompressedBitmapIterator {
			switch(self.buckets[0]) {
				.bitmap => |bucket| {
					return CompressedBitmapIterator {
						.bucket = 0,
						.bitmap = self,
						.iterator = CompressedBucketIterator {
							.bitmap = bucket.iterator()
						}
					};
				},
				.array => |bucket| {
					return CompressedBitmapIterator {
						.bucket = 0,
						.bitmap = self,
						.iterator = CompressedBucketIterator {
							.array = bucket.iterator()
						}
					};
				},
				.empty => {
						return CompressedBitmapIterator {
						.bucket = 0,
						.bitmap = self,
						.iterator = CompressedBucketIterator {
							.empty = BucketEmptyIterator{}
						}
					};				
				}
			}
		}

		pub inline fn initEmpty(size: u32, allocator: std.mem.Allocator) !CompressedBitmap {
			const num_buckets = (size>>16)+1;
			var buckets = try allocator.alloc(CompressedBucket, num_buckets);
			var i:u32 = 0;
			if(num_buckets > 1) {
				while(i < num_buckets) : (i+=1) {
					buckets[i] = CompressedBucket {
						.empty = {}
					};
				}
			}
			else {
				buckets[0] = CompressedBucket {
					.bitmap = try BucketBitmap.initEmpty(allocator)
				};
			}

			return CompressedBitmap {
				.buckets = buckets,
				.allocator = allocator,
				.cardinality = 0
			};
		}

		pub inline fn get(self: *Self, index: u32) bool {
			const bucket_index = index&0xFFFF;
			const bucket = index>>16;
			switch(self.buckets[bucket]) {
				.bitmap => |bitmap_bucket| {
					return bitmap_bucket.get(bucket_index);
				},
				.array => |array_bucket| {
					return array_bucket.get(bucket_index);
				},
				.empty => {
					return false;
				}
			}
		}

		pub inline fn unset(self: *Self, index: u32) void {
			const bucket_index = index&0xFFFF;
			const bucket = index>>16;
			switch(self.buckets[bucket]) {
				.bitmap => |*bitmap_bucket| {
					if(bitmap_bucket.unset(bucket_index)) {
						self.cardinality-=1;
					}
				},
				.array => |*array_bucket| {
					if(array_bucket.unset(bucket_index)) {
						self.cardinality-=1;
					}
				},
				.empty => {}
			}
		}
		
		pub inline fn set(self: *Self, index: u32) !void {
			const bucket_index = index&0xFFFF;
			const bucket = index>>16;
			switch(self.buckets[bucket]) {
				.bitmap => |*bitmap_bucket| {
					if(bitmap_bucket.set(bucket_index)) {
						self.cardinality+=1;
					}
				},
				.array => |*array_bucket| {
					if(try array_bucket.set(bucket_index)) {
						self.cardinality+=1;
					}
					if(array_bucket.cardinality() > MAX_SPARSE_BUCKET_SIZE) {
						var bitmap_bucket = try BucketBitmap.initEmpty(self.allocator);

						var iterator_arr = array_bucket.iterator();
						while(iterator_arr.next()) |item| {
							_ = bitmap_bucket.set(item);
						}
						array_bucket.storage.deinit();
						self.buckets[bucket] = CompressedBucket {
							.bitmap = bitmap_bucket
						};
					}
				},
				.empty => {
					var bucket_arr = BucketArray {
						.storage = std.ArrayList(u16).init(self.allocator)
					};
					_ = try bucket_arr.set(bucket_index);
					self.cardinality+=1;
					self.buckets[bucket] = CompressedBucket {
						.array = bucket_arr
					};
				}
			}
		}

		pub fn deinit(self: *Self) void {
			for(self.buckets) |*bucket| {
				switch(bucket.*) {
					.bitmap => |bitmap_bucket| {
						self.allocator.free(bitmap_bucket.storage);
					},
					.array => |array_bucket| {
						array_bucket.storage.deinit();
					},
					.empty => {}
				}
			}
			self.allocator.free(self.buckets);
		}
};


test "CompressedBitmap: Check get" {
	var gpa = std.heap.GeneralPurposeAllocator(.{}){};
	defer std.debug.assert(!gpa.deinit());



	var empty = try CompressedBitmap.initEmpty(2000,gpa.allocator());
	defer empty.deinit();	

	var i:u32 = 0;
	while(i<2000) : (i+=1) {
		std.debug.assert(empty.get(i) == false);
	}
	try std.testing.expectEqual(empty.cardinality,0);
}

test "CompressedBitmap: Check set" {
	var gpa = std.heap.GeneralPurposeAllocator(.{}){};
	defer std.debug.assert(!gpa.deinit());



	var empty = try CompressedBitmap.initEmpty(2000,gpa.allocator());
	defer empty.deinit();	

	var i:u32 = 0;
	while(i<2000) : (i+=1) {
		try empty.set(i);
	}
	i = 0;
	while(i<2000) : (i+=1) {
		std.debug.assert(empty.get(i) == true);
	}
	try std.testing.expectEqual(empty.cardinality,2000);
}

test "CompressedBitmap: Check iterator" {
	var gpa = std.heap.GeneralPurposeAllocator(.{}){};
	defer std.debug.assert(!gpa.deinit());



	var empty = try CompressedBitmap.initEmpty(2000,gpa.allocator());
	defer empty.deinit();	

	var i:u32 = 0;
	var check:u32 = 0;
	while(i<2000) : (i+=100) {
		try empty.set(i);
		check+=i;
	}

	var iter = empty.iterator();
	var calc: u32 = 0;
	while(iter.next()) |item| {
		calc+=item;
	}

	try std.testing.expectEqual(calc,check);
	try std.testing.expectEqual(empty.cardinality,20);
}

const builtin = @import("builtin");
test "bench: CompressedBitmap: 3M iterator" {
	if(builtin.mode != std.builtin.OptimizeMode.ReleaseFast) {
		return error.SkipZigTest;
	}
	var gpa = std.heap.GeneralPurposeAllocator(.{
		.enable_memory_limit = true,
	}){};
	gpa.requested_memory_limit = 5000;
	defer _ = gpa.deinit();
	var std_bitset = try CompressedBitmap.initEmpty(3_000_000,gpa.allocator());

	var calc:u64 = 0;
	var i:u64 = 0;
	while(i<1000) : (i+=1) {
		try std_bitset.set(100);
		try std_bitset.set(900);
		try std_bitset.set(1_500_000);
		var iter = std_bitset.iterator();
		while(iter.next()) |item| {
			calc+=@intCast(u64,item);
		}
	}

	try std.testing.expectEqual(calc,(1000+1_500_000)*1000);
	try std.testing.expectEqual(std_bitset.cardinality,3);
}

test "bench: CompressedBitmap: 3M xor iterator" {
	if(builtin.mode != std.builtin.OptimizeMode.ReleaseFast) {
		return error.SkipZigTest;
	}
	var gpa = std.heap.GeneralPurposeAllocator(.{
		.enable_memory_limit = true,
	}){};
	gpa.requested_memory_limit = 10000;
	defer _ = gpa.deinit();
	var std_bitset = try CompressedBitmap.initEmpty(3_000_000,gpa.allocator());
	var second_bitset = try CompressedBitmap.initEmpty(3_000_000,gpa.allocator());

	var calc:u64 = 0;
	var i:u64 = 0;
	while(i<1000) : (i+=1) {
		try std_bitset.set(100);
		try second_bitset.set(101);
		try std_bitset.set(900);
		try second_bitset.set(901);
		try second_bitset.set(900);
		try std_bitset.set(1_500_000);
		try second_bitset.set(1_500_000);

		var iter = std_bitset.xor_iterator(&second_bitset);
		while(iter.next()) |item| {
			calc+=@intCast(u64,item);
		}
	}

	try std.testing.expectEqual(calc,(101+100+901)*1000);
	try std.testing.expectEqual(std_bitset.cardinality,3);
	try std.testing.expectEqual(second_bitset.cardinality,4);
}
