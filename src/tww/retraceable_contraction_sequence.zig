const std = @import("std");
const contraction = @import("contraction_sequence.zig");
const comptime_util = @import("../util/comptime_checks.zig");

pub const RetraceableContractionError = error {
	NodeNotErasedOrNotParent,
	NodeIdTooLarge
};

pub fn RetraceableContractionSequence(comptime T: type) type {
	comptime if(!comptime_util.checkIfIsCompatibleInteger(T)) {
		@compileError("The type must either be u8,u16 or u32!");
	};

	return struct {
	const Self = @This();
	seq: contraction.ContractionSequence(T),
	tww: []T,
	parent_array: []T,
	write_ptr: u32,
	allocator: std.mem.Allocator,

	pub inline fn init(allocator: std.mem.Allocator, graph_size: u32) !Self {
		std.debug.assert(graph_size <= std.math.maxInt(T));

		const seq = try contraction.ContractionSequence(T).init(allocator,graph_size);
		const tww_mem = try allocator.alloc(T,graph_size-1);
		const parent_array = try allocator.alloc(T,graph_size);
		for(parent_array,0..) |*p,index| {
			p.* = @intCast(T,index);
		}

		return Self {
			.seq = seq,
			.allocator = allocator,
			.write_ptr = 0,
			.parent_array = parent_array,
			.tww = tww_mem
		};
	}

	pub inline fn getTwinWidth(self: *const Self) T {
		if(self.write_ptr == 0) {
			return 0;
		}
		return self.tww[self.write_ptr-1];
	}

	pub inline fn removeLast(self: *Self) contraction.ContractionError!void {
		const cont = try self.seq.removeLastContraction();
		self.write_ptr -= 1;
		// Restore parent array
		self.parent_array[cont.erased] = cont.erased;
	}

	pub inline fn isComplete(self: *Self) bool {
		return self.seq.isComplete();
	}

	pub inline fn addContraction(self: *Self, deleted: T, survivor: T, tww: T) !void {
		if (self.write_ptr == self.tww.len) {
			return error.ContractionOverflow;
		}
		else if(deleted >= self.parent_array.len or survivor >= self.parent_array.len) {
			return error.NodeIdTooLarge;
		}
		self.parent_array[deleted] = survivor;
		self.tww[self.write_ptr] = tww;
		try self.seq.addContraction(.{.erased = deleted, .survivor = survivor});
		self.write_ptr += 1;
	}

	pub inline fn retraceParentOfRedEdge(self: *const Self, parent_canidate_1: T, parent_canidate_2: T, red_edge_induced_by_erasure_of: T) RetraceableContractionError!T {
		var start = red_edge_induced_by_erasure_of;
		if(start == parent_canidate_1) {
			return parent_canidate_1;
		}
		else if(start == parent_canidate_2) {
			return parent_canidate_2;
		}

		while(true) {
			if(start == parent_canidate_1) {
				return parent_canidate_1;
			}
			else if(start == parent_canidate_2) {
				return parent_canidate_2;
			}

			const new_parent = self.parent_array[start];

			if(new_parent == start) {
				return RetraceableContractionError.NodeNotErasedOrNotParent;
			}
			start = new_parent;
		}
	}

	pub inline fn lastContraction(self: *Self) ?contraction.Contraction(T) {
		return self.seq.getLastContraction();
	}

	pub inline fn deinit(self: *const Self) void {
		self.seq.deinit(self.allocator);
		self.allocator.free(self.tww);
		self.allocator.free(self.parent_array);
	}
};

}


const expectEqual = std.testing.expectEqual;

test "RetraceableContractionSequence: Check normal add" {
	var gpa = std.heap.GeneralPurposeAllocator(.{}){};
	defer std.debug.assert(!gpa.deinit());
	
	var traceable = try RetraceableContractionSequence(u32).init(gpa.allocator(),10);
	defer traceable.deinit();

	try std.testing.expectEqual(traceable.parent_array.len,traceable.tww.len+1);
	try std.testing.expectEqual(traceable.seq.list.len,traceable.tww.len);

	try traceable.addContraction(8,7,100);

	try expectEqual(traceable.getTwinWidth(),100);
	try expectEqual(traceable.lastContraction(),.{.erased=8,.survivor=7});
	try expectEqual(traceable.parent_array[8],7);
}

test "RetraceableContractionSequence: Remove contractions" {
	var gpa = std.heap.GeneralPurposeAllocator(.{}){};
	defer std.debug.assert(!gpa.deinit());
	
	var traceable = try RetraceableContractionSequence(u32).init(gpa.allocator(),10);
	defer traceable.deinit();

	try std.testing.expectEqual(traceable.parent_array.len,traceable.tww.len+1);
	try std.testing.expectEqual(traceable.seq.list.len,traceable.tww.len);

	try traceable.addContraction(8,7,100);

	try expectEqual(traceable.getTwinWidth(),100);
	try expectEqual(traceable.lastContraction(),.{.erased=8,.survivor=7});
	try expectEqual(traceable.parent_array[8],7);

	try traceable.removeLast();
	try expectEqual(traceable.getTwinWidth(),0);
	try expectEqual(traceable.seq.write_ptr,0);
	try expectEqual(traceable.lastContraction(),null);
	try expectEqual(traceable.parent_array[8],8);
}

test "RetraceableContractionSequence: Retrace parent" {
	var gpa = std.heap.GeneralPurposeAllocator(.{}){};
	defer std.debug.assert(!gpa.deinit());
	
	var traceable = try RetraceableContractionSequence(u32).init(gpa.allocator(),10);
	defer traceable.deinit();

	try std.testing.expectEqual(traceable.parent_array.len,traceable.tww.len+1);
	try std.testing.expectEqual(traceable.seq.list.len,traceable.tww.len);

	try traceable.addContraction(8,7,80);
	try traceable.addContraction(7,9,90);
	try traceable.addContraction(9,2,100);
	try traceable.addContraction(5,2,110);

	try expectEqual(traceable.getTwinWidth(),110);
	try expectEqual(traceable.lastContraction(),.{.erased=5,.survivor=2});

	try expectEqual(traceable.retraceParentOfRedEdge(2,5,8), 2);
	try expectEqual(traceable.retraceParentOfRedEdge(2,9,8), 9);
	try expectEqual(traceable.retraceParentOfRedEdge(2,9,9), 9);
	try expectEqual(traceable.retraceParentOfRedEdge(2,9,2), 2);

	try traceable.removeLast();
	try traceable.removeLast();
	try traceable.removeLast();
	try expectEqual(traceable.lastContraction(),.{.erased=8,.survivor=7});
	try expectEqual(traceable.getTwinWidth(),80);
	try std.testing.expectError(error.NodeNotErasedOrNotParent,traceable.retraceParentOfRedEdge(2,9,8));
}

test "RetraceableContractionSequence: Overflow" {
	var gpa = std.heap.GeneralPurposeAllocator(.{}){};
	defer std.debug.assert(!gpa.deinit());
	
	var traceable = try RetraceableContractionSequence(u32).init(gpa.allocator(),2);
	defer traceable.deinit();

	try std.testing.expectEqual(traceable.parent_array.len,traceable.tww.len+1);
	try std.testing.expectEqual(traceable.seq.list.len,traceable.tww.len);

	try std.testing.expectError(RetraceableContractionError.NodeIdTooLarge,traceable.addContraction(8,7,80));
	try traceable.addContraction(0,1,90);
	try std.testing.expectError(contraction.ContractionError.ContractionOverflow,traceable.addContraction(1,0,80));

	try expectEqual(traceable.isComplete(), true);
}
