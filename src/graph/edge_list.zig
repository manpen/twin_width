const std = @import("std");
const bitset = @import("../util/two_level_bitset.zig");
const comptime_util = @import("../util/comptime_checks.zig");

pub const EdgeListType = enum(u8) {
	bitmap,
	list
};

pub fn ParametrizedUnsortedArrayList(comptime T: type) type {
	if (!comptime_util.checkIfIsCompatibleInteger(T)) {
		@compileError("Type must either be u8,16 or u32!");
	}

	return struct {
		const Self = @This();
		edges: std.ArrayListUnmanaged(T),
		pub inline fn initCapacity(allocator: std.mem.Allocator, size: u32) !Self {
			std.debug.assert(size>0);
			var list = try std.ArrayListUnmanaged(T).initCapacity(allocator,size);

			return Self {
				.edges = list
			};
		}

		pub inline fn init() Self {
			var list = std.ArrayListUnmanaged(T){};

			return Self {
				.edges = list
			};
		}

		pub const ParametrizedUnsortedArrayListIterator = struct {
			index: u32,
			list: *const Self,
			pub inline fn next(self: *ParametrizedUnsortedArrayListIterator) ?T {
				if(self.index>=self.list.edges.items.len) return null;
				const item = self.list.edges.items[self.index];
				self.index+=1;
				return item;
			}
		};

		pub inline fn intoSorted(self: *Self) ParametrizedSortedArrayList(T) {
			std.sort.sort(T,self.edges.items,{}, comptime std.sort.asc(T));
			return ParametrizedSortedArrayList(T) {
				.edges = self.edges
			};
		}

		pub inline fn iterator(self: *const Self) ParametrizedUnsortedArrayListIterator {
			return ParametrizedUnsortedArrayListIterator {
				.index = 0,
				.list = self
			};
		}

		pub inline fn add(self: *Self, allocator: std.mem.Allocator, item: T) !void {
			try self.edges.append(allocator,item);
		}

		pub inline fn cardinality(self: *const Self) u32 {
			return @intCast(u32,self.edges.items.len);
		}

		pub fn deinit(self: *Self, allocator: std.mem.Allocator) void {
			self.edges.deinit(allocator);
		}

		pub inline fn nodePosition(self: *const Self, item: T) ?u32 {
			var i: u32 = 0;
			while(i < self.edges.items.len) : (i+=1) {
				if(self.edges.items[i] == item) {
					return i;
				}
			}
			return null;
		}

		pub inline fn contains(self: *const Self, item: T) bool {
			return self.nodePosition(item) != null;
		}

		pub inline fn remove(self: *Self, item: T) bool {
			const result = self.nodePosition(item);
			if(result) |position| {
				_ = self.edges.swapRemove(position);
				return true;
			}
			return false;
		}
	};
}

pub fn ParametrizedSortedArrayList(comptime T: type) type {
	if (!comptime_util.checkIfIsCompatibleInteger(T)) {
		@compileError("Type must either be u8,16 or u32!");
	}

	return struct {
		const Self = @This();
		edges: std.ArrayListUnmanaged(T),
		pub inline fn initCapacity(allocator: std.mem.Allocator, size: u32) !Self {
			std.debug.assert(size>0);
			var list = try std.ArrayListUnmanaged(T).initCapacity(allocator,size);

			return Self {
				.edges = list
			};
		}


		pub inline fn init() Self {
			var list = std.ArrayListUnmanaged(T){};
			return Self {
				.edges = list
			};
		}

		pub const ParametrizedSortedArrayListIterator = struct {
			index: u32,
			list: *const Self,
			pub inline fn next(self: *ParametrizedSortedArrayListIterator) ?T {
				if(self.index>=self.list.edges.items.len) return null;
				const item = self.list.edges.items[self.index];
				self.index+=1;
				return item;
			}
		};

		pub inline fn iterator(self: *const Self) ParametrizedSortedArrayListIterator {
			return ParametrizedSortedArrayListIterator {
				.index = 0,
				.list = self
			};
		}

		pub inline fn xorIterator(self: *const Self, other: *const Self) ParametrizedSortedArrayListXorIterator {
			var iterator_first = self.iterator();
			var iterator_second = other.iterator();
			const item_first = iterator_first.next();
			const item_second = iterator_second.next();
			return ParametrizedSortedArrayListXorIterator {
				.iterator_first = iterator_first,
				.iterator_second = iterator_second,
				.item_first = item_first,
				.item_second = item_second,
				.first = false
			};
		}

		pub const ParametrizedSortedArrayListXorIterator = struct {
			iterator_first: ParametrizedSortedArrayListIterator,
			iterator_second: ParametrizedSortedArrayListIterator,
			item_first: ?T,
			item_second: ?T,
			first: bool,
			pub inline fn next(self: *ParametrizedSortedArrayListXorIterator) ?T {
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

		pub inline fn add(self: *Self, allocator: std.mem.Allocator, item: T) !bool {
			const result = self.nodePosition(item);
			if(result.@"0" == false) {
				try self.edges.insert(allocator,result.@"1",item);
			}
			return result.@"0";
		}

		pub inline fn cardinality(self: *const Self) u32 {
			return @intCast(u32,self.edges.items.len);
		}

		pub fn deinit(self: *Self, allocator: std.mem.Allocator) void {
			self.edges.deinit(allocator);
		}

		pub inline fn contains(self: *const Self, item: T) bool {
			return self.nodePosition(item).@"0";
		}

		pub inline fn nodePosition(self: *const Self, item: T) struct{bool,u32} {
			// Larger than 4 cache lines
			//const threshold = comptime (64/@sizeOf(T))*4;
			//TODO: Change this to a sensible number again!
			if(self.edges.items.len >= 1_000_000) {
				return self.binarySearch(item);
			}
			else {
				var i: u32 = 0;
				while(i < self.edges.items.len) : (i+=1) {
					if(self.edges.items[i] == item) {
						return .{true,i};
					}
					else if(self.edges.items[i] > item) {
						return .{false,i};
					}
				}
				return .{false,@intCast(u32,self.edges.items.len)};
			}
		}

		pub fn binarySearch(self: *const Self, target: T) struct{bool,u32} {
			var left: usize = 0;
			var right = self.edges.items.len;

			while (left < right) {
				const mid = left + (right - left) / 2; // Avoid overflow.
				if (self.edges.items[mid] == target) {
					return .{false,@intCast(u32,mid)};
				} else if (self.edges.items[mid] < target) {
					left = mid + 1;
				} else {
					right = mid;
				}
			}
			return .{false,@intCast(u32,left)};
		}

		pub inline fn remove(self: *Self, item: T) bool {
			const result = self.nodePosition(item);
			if(result.@"0" == true) {
				_ = self.edges.orderedRemove(result.@"1");
			}
			return result.@"0";
		}
	};
}
