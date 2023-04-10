const std = @import("std");
const comptime_util = @import("comptime_checks.zig");
const compressed_bitmap = @import("compressed_bitmap.zig");

pub inline fn fisher_yates_shuffle(comptime T: type, data: []T) void {
	comptime if( !comptime_util.checkIfIsCompatibleInteger(T) ) {
		@compileError("Can only use the fisher yates shuffle on integer types u8,u16 or u32");
	};
	
	var default = std.rand.DefaultPrng{};
	default.init(@intCast(u64,std.time.timestamp()));

	for(0..(data.len-1)) |i| {
		const random = default.next();
		const target = i+random%(data.len-i);
		
		const tmp = data[target];
		data[target] = data[i];
		data[i] = tmp;
	}
}


pub inline fn circular_permutation_shift(comptime T: type, data: []T) void {
	// See https://openreview.net/pdf?id=NrkAAcMpRoT

	const item = data[data.len-1];

	var next = data[0];

	for(0..(data.len-1)) |i| {
		const tmp = data[i+1];
		data[i+1] = next;
		next = tmp;
	}

	data[0] = item;
}


pub fn MinHash(comptime T: type) type {
	return struct {
		const Self = @This();
		permutation: []T,
		pub inline fn init(number_of_nodes: T, allocator: std.mem.Allocator) !Self {
			var memory = try allocator.alloc(T,number_of_nodes);
			for(0..number_of_nodes) |i| {
				memory[i] = @intCast(T,i);
			}
			
			fisher_yates_shuffle(T,memory);

			return Self {
				.permutation = memory,
			};
		}

		pub inline fn deinit(self: *Self, allocator: std.mem.Allocator) void {
			allocator.free(self.permutation);
		}

		pub inline fn hash(self: *const Self, input: *const compressed_bitmap.FastCompressedBitmap) T {
			var min = std.math.maxInt(T);
			var iter = input.iterator();
			while(iter.next()) |item| {
				min = std.math.min(min, self.permutation[item]);
			}
			return min;
		}
	};
}
