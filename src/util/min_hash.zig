const std = @import("std");
const comptime_util = @import("comptime_checks.zig");
const compressed_bitmap = @import("compressed_bitmap.zig");
const two_level_bitset = @import("../util/two_level_bitset.zig");

pub inline fn fisher_yates_shuffle(comptime T: type, data: []T) void {
	comptime if( !comptime_util.checkIfIsCompatibleInteger(T) ) {
		@compileError("Can only use the fisher yates shuffle on integer types u8,u16 or u32");
	};
	
	var default = std.rand.DefaultPrng.init(@intCast(u64,std.time.timestamp()));

	for(0..(data.len-1)) |i| {
		const random = default.next();
		const target = i+random%(data.len-i);
		
		const tmp = data[target];
		data[target] = data[i];
		data[i] = tmp;
	}
}

pub inline fn fisher_yates_sample_first_n(comptime T: type, data: []T, n: u32) void {
	comptime if( !comptime_util.checkIfIsCompatibleInteger(T) ) {
		@compileError("Can only use the fisher yates shuffle on integer types u8,u16 or u32");
	};
	
	var default = std.rand.DefaultPrng.init(@intCast(u64,std.time.timestamp()));

	for(0..(std.math.min(n,data.len-1))) |i| {
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


pub fn MinHash(comptime T: type, comptime promote_thresh:u32, comptime degrade_thresh:u32) type {
	return struct {
		const Self = @This();
		permutation: []T,
		pub inline fn init(allocator: std.mem.Allocator, number_of_nodes: T) !Self {
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

		pub inline fn hash(self: *const Self, input: *const compressed_bitmap.FastCompressedBitmap(T,promote_thresh,degrade_thresh)) T {
			var min:T = std.math.maxInt(T);
			var iter = input.iterator();
			while(iter.next()) |item| {
				min = std.math.min(min, self.permutation[item]);
			}
			return min;
		}
	};
}

pub fn MinHashBand(comptime T: type, comptime N: u32, comptime promote_thresh:u32, comptime degrade_thresh:u32) type {
	comptime if(N > 7) {
		@compileError("At the moment N must be smaller than 8!");
	};

	return struct {
		const Self = @This();
		hashes: [N]MinHash(T,promote_thresh,degrade_thresh),
		tables: [N]std.AutoHashMapUnmanaged(T,std.ArrayListUnmanaged(T)),
		sim_nodes: two_level_bitset.FastBitSet, 
		node_num_sim: []u3,

		pub inline fn init(allocator: std.mem.Allocator, number_of_nodes: u32) !Self {
			var alloc = try allocator.alloc(u3,number_of_nodes);
			std.mem.set(u3,alloc,0);

			var hb = Self {
				.hashes = undefined,
				.tables = undefined,
				.sim_nodes = try two_level_bitset.FastBitSet.initEmpty(number_of_nodes,allocator),
				.node_num_sim = alloc,
			};

			for(&hb.hashes) |*hash| {
				hash.* = try MinHash(T,promote_thresh,degrade_thresh).init(allocator, @intCast(T,number_of_nodes));
			}
			for(&hb.tables) |*hash| {
				hash.* = std.AutoHashMapUnmanaged(T,std.ArrayListUnmanaged(T)){};
			}

			return hb;
		}

		pub fn deinit(self: *Self, allocator: std.mem.Allocator) void {
			allocator.free(self.node_num_sim);
			self.sim_nodes.deinit(allocator);
			for(self.tables) |*hash| {
				hash.deinit(allocator);
			}
		}

		pub inline fn removeNode(self: *Self, node: T, nbh: *compressed_bitmap.FastCompressedBitmap(T,promote_thresh,degrade_thresh)) !void {
			for(0..self.hashes.len) |index| {
				const hash = self.hashes[index].hash(nbh);
				if(self.tables[index].getPtr(hash)) |table| {
					var i:u32 = 0;
					while(i < table.items.len) : (i+=1) {
						if(table.items[i] == node) {
							_ = table.swapRemove(i);
							break;
						}
					}
					_ = self.sim_nodes.unset(node);

					if (table.items.len == 1) {
						self.node_num_sim[node] = 0;
						self.node_num_sim[table.items[0]] -= 1;

						if(self.node_num_sim[table.items[0]] == 0) {
							_ = self.sim_nodes.unset(table.items[0]);
						}
					}
				}
			}
		}

		pub inline fn hashNodeNeighborhood(self: *Self, allocator: std.mem.Allocator, node: T, nbh: *compressed_bitmap.FastCompressedBitmap(T,promote_thresh,degrade_thresh)) !void {
			for(0..self.hashes.len) |index| {
				const hash = self.hashes[index].hash(nbh);
				if(self.tables[index].getPtr(hash)) |table| {
					try table.append(allocator, node);
					if (table.items.len == 2) {
						self.sim_nodes.set(table.items[0]);
						self.node_num_sim[table.items[0]] += 1;
					}
					
					self.sim_nodes.set(node);
					self.node_num_sim[node] += 1;
				}
				else {
					var arr = try std.ArrayListUnmanaged(T).initCapacity(allocator,8);
					arr.append(allocator,node) catch unreachable;
					try self.tables[index].put(allocator, hash, arr);
				}
			}
		}
	};
}

pub const MinHashEntry = struct {
	counter: u32,	
	timestamp: u32
};

pub fn MinHashSimiliarity(comptime T: type, comptime N: u32, comptime B: u32, comptime promote_thresh:u32, comptime degrade_thresh:u32) type {
	return struct {
		const Self = @This();
		bands: [B]MinHashBand(T,N,promote_thresh,degrade_thresh),
		score_table: []MinHashEntry,
		unique_id: u32,

		pub inline fn init(allocator: std.mem.Allocator, number_of_nodes: T) !Self {
			var mhs = Self {
				.bands = undefined,
				.score_table = try allocator.alloc(MinHashEntry,number_of_nodes),
				.unique_id = 1,
			};

			for(&mhs.bands) |*band| {
				band.* = try MinHashBand(T,N,promote_thresh,degrade_thresh).init(allocator, number_of_nodes);
			}

			for(mhs.score_table) |*x| {
				x.timestamp = 0;
			}
			
			return mhs;
		}

		pub inline fn hashNodeNeighborhood(self: *Self, allocator: std.mem.Allocator, node: T, nbh: *compressed_bitmap.FastCompressedBitmap(T,promote_thresh,degrade_thresh)) !void {
			for(&self.bands) |*band| {
				try band.hashNodeNeighborhood(allocator,node,nbh);
			}
		}

		pub inline fn removeNode(self: *Self, allocator: std.mem.Allocator, node: T, nbh: *compressed_bitmap.FastCompressedBitmap(T,promote_thresh,degrade_thresh)) !void {
			_ = allocator;
			for(&self.bands) |*band| {
				try band.removeNode(node,nbh);
			}
		}

		pub inline fn getMostSimNode(self: *Self) ?T {
			self.unique_id += 1;
			var max_node:?T = null;
			var max_score:u32 = 0;

			for(0..B) |i| {
				var iter = self.bands[i].sim_nodes.iter();
				while(iter.next()) |item| {
					if(self.score_table[item].timestamp != self.unique_id) {
						self.score_table[item].counter = 1;
						self.score_table[item].timestamp = self.unique_id;
					}
					else {
						self.score_table[item].counter += 1;
					}

					if(self.score_table[item].counter > max_score) {
						max_score = self.score_table[item].counter;
						max_node = @intCast(T,item);
					}
				}
			}
			return max_node;
		}
	};
}
