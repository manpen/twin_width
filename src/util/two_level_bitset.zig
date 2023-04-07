const std = @import("std");

pub const FastBitSetIterator = struct {
	const Self = @This();
	coerced_ptr: [*]u64,
	coerced_higher_ptr: [*]u64,
	current_higher: u64,
	length_higher: usize,
	index_higher: u32,
	current: u64,
	index: u32,
	pub inline fn next(self: *Self) ?u32 {
				 while(self.current == 0) {
					 while(self.current_higher == 0) {
						 self.index_higher+=1;
						 if(self.index_higher>=self.length_higher) {
							 return null;
						 }
						 self.current_higher = self.coerced_higher_ptr[self.index_higher];
					 }
					 const next_higher_index = @ctz(self.current_higher) + (self.index_higher<<6);
					 self.current_higher &= self.current_higher - 1;

					 self.current = self.coerced_ptr[next_higher_index];
					 self.index = next_higher_index;
				 }

				 const next_index = @ctz(self.current) + (self.index<<6);
				 self.current &= self.current - 1;
				 return next_index;
			 }
};

pub const FastBitSetZeroIterator = struct {
	const Self = @This();
	coerced_ptr: [*]u64,
	length: usize,
	current: u64,
	index: u32,
	max_index: u32,
	pub inline fn next(self: *Self) ?u32 {
		while(self.current==0) {
			self.index+=1;
			if(self.index>=self.length) {
				return null;
			}
			self.current = ~self.coerced_ptr[self.index];
		}

		const next_index = @ctz(self.current) + (self.index<<6);
		self.current &= self.current - 1;
		if(next_index>=self.max_index) {
			return null;
		}
		return next_index;
	}
};

pub const FastBitSet = struct {
		const Self = @This();
		storage: []u64,
		storage_higher: []u64,
		allocator: std.mem.Allocator,
		cardinality: u32,
		total_len: u32,

		pub inline fn initEmpty(size: u32, allocator: std.mem.Allocator) !Self {
			var ensure_aligned = (size>>6)+1;
			var ensure_higher = (ensure_aligned>>6)+1;
			var storage = try allocator.alignedAlloc(u64,32,ensure_aligned);
			errdefer allocator.free(storage);
			var storage_higher = try allocator.alignedAlloc(u64,32,ensure_higher);
			errdefer allocator.free(storage_higher);
			std.mem.set(u64,storage,0);
			std.mem.set(u64,storage_higher,0);

			return Self {
				.allocator = allocator,
				.storage = storage,
				.storage_higher = storage_higher,
				.cardinality = 0,
				.total_len = size
			};
		}

		pub inline fn copy(self: *const Self, allocator: std.mem.Allocator) !Self {
			var ensure_aligned = self.storage.len;
			var ensure_higher = self.storage_higher.len;
			var storage = try allocator.alignedAlloc(u64,32,ensure_aligned);
			errdefer allocator.free(storage);
			var storage_higher = try allocator.alignedAlloc(u64,32,ensure_higher);
			errdefer allocator.free(storage_higher);
			std.mem.copy(u64,storage,self.storage);
			std.mem.copy(u64,storage_higher,self.storage_higher);

			return Self {
				.allocator = allocator,
				.storage = storage,
				.storage_higher = storage_higher,
				.cardinality = self.cardinality,
				.total_len = self.total_len
			};
		}

		pub fn deinit(self: *Self) void {
			self.allocator.free(self.storage);
			self.allocator.free(self.storage_higher);
		}

		pub inline fn iter(self: *Self) FastBitSetIterator {
			return FastBitSetIterator {
				.coerced_ptr = self.storage.ptr,
				.coerced_higher_ptr = self.storage_higher.ptr,
				.length_higher = self.storage_higher.len,
				.current_higher = self.storage_higher[0],
				.index_higher = 0,
				.current = 0,
				.index = 0,
			};
		}

		pub inline fn iterUnset(self: *Self) FastBitSetZeroIterator {
			return FastBitSetZeroIterator {
				.coerced_ptr = self.storage.ptr,
				.length = self.storage.len,
				.max_index = self.total_len,
				.current = ~self.storage[0],
				.index = 0,
			};
		}


		pub inline fn set(self: *Self, index: u32) void {
			const overflow:u6 = @intCast(u6,index&0x3F);
			const new_index = index>>6;
			const overflow_higher:u6 = @intCast(u6,new_index&0x3F);
			const new_index_higher = index>>12;
			const before = self.storage[new_index];
			self.storage[new_index] |= @as(u64,1)<<overflow;
			self.storage_higher[new_index_higher] |= @as(u64,1)<<overflow_higher;
			if((before&@as(u64,1)<<overflow)==0) {
				self.cardinality+=1;
			}
		}


		pub inline fn setExists(self: *Self, index: u32) bool {
			const overflow:u6 = @intCast(u6,index&0x3F);
			const new_index = index>>6;
			const overflow_higher:u6 = @intCast(u6,new_index&0x3F);
			const new_index_higher = index>>12;
			const before = self.storage[new_index];
			self.storage[new_index] |= @as(u64,1)<<overflow;
			self.storage_higher[new_index_higher] |= @as(u64,1)<<overflow_higher;
			if((before&@as(u64,1)<<overflow)==0) {
				self.cardinality+=1;
				return false;
			}
			return true;
		}

		pub inline fn unset(self: *Self, index: u32) bool {
			const overflow:u6 = @intCast(u6,index&0x3F);
			const new_index = index>>6;
			const before = self.storage[new_index];
			self.storage[new_index] &= ~(@as(u64,1)<<overflow);
			if((before&(@as(u64,1)<<overflow))!=0) {
				self.cardinality-=1;
				return true;
			}
			return false;
		}

		pub inline fn get(self: *Self, index: u32) bool {
			const overflow:u6 = @intCast(u6,index&0x3F);
			const new_index = index>>6;
			return (self.storage[new_index] & @as(u64,1)<<overflow)!=0;
		}

		pub inline fn unsetAll(self: *Self) void {
			if(self.cardinality > 0) {
				std.mem.set(u64,self.storage,0);
				std.mem.set(u64,self.storage_higher,0);
				self.cardinality = 0;
			}
		}

		pub inline fn setUnion(self: *Self, other: *Self) void {
			if(other.cardinality < (self.storage.len>>8)) {
				var iterator = other.iter();
				while(iterator.next()) |item| {
					self.set(item);
				}
			}
			else {
				const length = std.math.min(self.storage.len,other.storage.len);
				var i:u32 = 0;
				while(i < length) : (i+=1) {
					self.storage[i] |= other.storage[i];
				}
				const length_higher = std.math.min(self.storage_higher.len,other.storage_higher.len);
				i=0;
				while(i < length_higher) : (i+=1) {
					self.storage_higher[i] |= other.storage_higher[i];
				}
			}
		}
};

const builtin = @import("builtin");

test "Check TwoLevelBitset iterator" {
	var gpa = std.heap.GeneralPurposeAllocator(.{}){};
	defer std.debug.assert(!gpa.deinit());
	var std_bitset = try FastBitSet.initEmpty(50_000,gpa.allocator());
	defer std_bitset.deinit();

	var random_gen = std.rand.DefaultPrng.init(@intCast(u64,std.time.timestamp()));

	var total_steps:u32 = 0;
	while(total_steps < 10) : (total_steps+=1) {
		var i:u32 = 0;
		var step_size:u64 = std.rand.limitRangeBiased(u64,random_gen.next(),49_000);

		var check: u64 = 0;
		while(i<50_000) : (i+=@intCast(u32,step_size)) {
			std_bitset.set(i);
			check+=i;
		}
		var iter = std_bitset.iter();
		while(iter.next()) |item| {
			check-=item;
		}
		try std.testing.expectEqual(check,0);

		std_bitset.unsetAll();
	}
}


test "bench: TwoLevelBitset 3M iterator" {
	if(builtin.mode != std.builtin.OptimizeMode.ReleaseFast) {
		return error.SkipZigTest;
	}
	var gpa = std.heap.GeneralPurposeAllocator(.{}){};
	defer _ = gpa.deinit();
	var std_bitset = try FastBitSet.initEmpty(3_000_000,gpa.allocator());

	var calc:u64 = 0;
	var i:u64 = 0;
	while(i<1000) : (i+=1) {
		std_bitset.set(100);
		std_bitset.set(900);
		std_bitset.set(1_500_000);
		var iter = std_bitset.iter();
		//var iter = std_bitset.iterator(.{});
		while(iter.next()) |item| {
			calc+=@intCast(u64,item);
		}
	}

	try std.testing.expectEqual(calc,(1000+1_500_000)*1000);
}

test "bench: std.bit_set.DynamicBitset 3M iterate" {
	if(builtin.mode != std.builtin.OptimizeMode.ReleaseFast) {
		return error.SkipZigTest;
	}
	var gpa = std.heap.GeneralPurposeAllocator(.{}){};
	defer _ = gpa.deinit();
	var std_bitset = try std.bit_set.DynamicBitSet.initEmpty(gpa.allocator(),3_000_000);

	var calc:u64 = 0;
	var i:u64 = 0;
	while(i<1000) : (i+=1) {
		std_bitset.set(100);
		std_bitset.set(900);
		std_bitset.set(1_500_000);
		var iter = std_bitset.iterator(.{});
		while(iter.next()) |item| {
			calc+=@intCast(u64,item);
		}
	}

	try std.testing.expectEqual(calc,(1000+1_500_000)*1000);
}
