const std = @import("std");

pub const FastBitSetIterator = struct {
	const Self = @This();
	coerced_ptr: [*]u64,
	coerced_higher_ptr: [*]u64,
	current_higher: u64,
	length_higher: usize,
	index_higher: u32,
	current: u64,
	length: usize,
	index: u32,
	pub inline fn next(self: *Self) ?u32 {
		again: while(true) {
				 if(self.index >= self.length) {
					 while(self.current_higher == 0) {
						 self.index_higher+=1;
						 if(self.index_higher>=self.length_higher) {
							 return null;
						 }
						 self.current_higher = self.coerced_higher_ptr[self.index_higher];
					 }
					 const next_higher_index = @ctz(self.current_higher) + (self.index_higher<<6);
					 self.current_higher &= self.current_higher - 1;

					 self.index = next_higher_index;
					 self.length = self.index+1;
					 self.current = self.coerced_ptr[self.index];
				 }

				 while(self.current==0) {
					 self.index+=1;
					 if(self.index>=self.length) {
							continue :again;
					 }
					 self.current = self.coerced_ptr[self.index];
				 }

				 const next_index = @ctz(self.current) + (self.index<<6);
				 self.current &= self.current - 1;
				 return next_index;
			 }
	}
};

pub const FastBitSet = struct {
		const Self = @This();
		storage: []u64,
		storage_higher: []u64,
		allocator: std.mem.Allocator,
		pub inline fn initEmpty(size: u32, allocator: std.mem.Allocator) !Self {
			var ensure_aligned = (size>>6)+1;
			var ensure_higher = (ensure_aligned>>6)+1;
			var storage = try allocator.alignedAlloc(u64,32,ensure_aligned);
			var storage_higher = try allocator.alignedAlloc(u64,32,ensure_higher);

			std.mem.set(u64,storage,0);
			std.mem.set(u64,storage_higher,0);

			return Self {
				.allocator = allocator,
				.storage = storage,
				.storage_higher = storage_higher
			};
		}

		pub inline fn iter(self: *Self) FastBitSetIterator {
			return FastBitSetIterator {
				.coerced_ptr = self.storage.ptr,
				.coerced_higher_ptr = self.storage_higher.ptr,
				.length_higher = self.storage_higher.len,
				.current_higher = self.storage_higher[0],
				.index_higher = 0,
				.length = 0,
				.current = undefined,
				.index = 0,
			};
		}


		pub inline fn set(self: *Self, index: u32) void {
			const overflow:u6 = @intCast(u6,index&0x3F);
			const new_index = index>>6;
			const overflow_higher:u6 = @intCast(u6,new_index&0x3F);
			const new_index_higher = index>>12;
			self.storage[new_index] |= @as(u64,1)<<overflow;
			self.storage_higher[new_index_higher] |= @as(u64,1)<<overflow_higher;
		}

		pub inline fn unset(self: *Self, index: u32) void {
			const overflow:u6 = @intCast(u6,index&0x3F);
			const new_index = index>>6;
			self.storage[new_index] &= ~(@as(u64,1)<<overflow);
		}

		pub inline fn get(self: *Self, index: u32) bool {
			const overflow:u6 = @intCast(u6,index&0x3F);
			const new_index = index>>6;
			return (self.storage[new_index] & @as(u64,1)<<overflow)!=0;
		}
};

const builtin = @import("builtin");

test "bench: MultiLevelBitset 3M iterator" {
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
