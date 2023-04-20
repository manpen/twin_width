const std = @import("std");
const contraction = @import("contraction_sequence.zig");
const comptime_util = @import("../util/comptime_checks.zig");
const red_edge_stack = @import("../graph/red_edge_stack.zig");

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
	red_edge_stack: red_edge_stack.RedEdgeStack(T),
	tww: []T,
	write_ptr: u32,

	pub inline fn init(allocator: std.mem.Allocator, graph_size: u32, number_of_edges: u32) !Self {
		std.debug.assert(graph_size <= std.math.maxInt(T));

		const seq = try contraction.ContractionSequence(T).init(allocator,graph_size);
		const tww_mem = try allocator.alloc(T,graph_size-1);

		// Should suffice for any contraction
		const res = try red_edge_stack.RedEdgeStack(T).initCapacity(allocator,2*number_of_edges);

		return Self {
			.seq = seq,
			.red_edge_stack = res,
			.write_ptr = 0,
			.tww = tww_mem
		};
	}

	pub inline fn getTwinWidth(self: *const Self) T {
		if(self.write_ptr == 0) {
			return 0;
		}
		return self.tww[self.write_ptr-1];
	}

	pub inline fn removeLast(self: *Self) !void {
		_ = try self.seq.removeLastContraction();
		try self.red_edge_stack.revertLastContraction();
		self.write_ptr -= 1;
	}

	pub inline fn isComplete(self: *Self) bool {
		return self.seq.isComplete();
	}

	pub inline fn addContraction(self: *Self, allocator: std.mem.Allocator, deleted: T, survivor: T, tww: T) !void {
		if (self.write_ptr == self.tww.len) {
			return error.ContractionOverflow;
		}
		self.tww[self.write_ptr] = tww;
		try self.seq.addContraction(.{.erased = deleted, .survivor = survivor});
		self.write_ptr += 1;
		try self.red_edge_stack.sealLevel(allocator);
	}

	pub inline fn lastContraction(self: *Self) ?contraction.Contraction(T) {
		return self.seq.getLastContraction();
	}

	pub inline fn deinit(self: *Self, allocator: std.mem.Allocator) void {
		self.seq.deinit(allocator);
		allocator.free(self.tww);
		self.red_edge_stack.deinit(allocator);
	}
};

}


const expectEqual = std.testing.expectEqual;

test "RetraceableContractionSequence: Check normal add" {
	var gpa = std.heap.GeneralPurposeAllocator(.{}){};
	defer std.debug.assert(!gpa.deinit());
	
	var traceable = try RetraceableContractionSequence(u32).init(gpa.allocator(),10,5);
	defer traceable.deinit(gpa.allocator());

	try std.testing.expectEqual(traceable.seq.list.len,traceable.tww.len);

	try traceable.addContraction(gpa.allocator(),8,7,100);

	try expectEqual(traceable.getTwinWidth(),100);
	try expectEqual(traceable.lastContraction(),.{.erased=8,.survivor=7});
}

test "RetraceableContractionSequence: Remove contractions" {
	var gpa = std.heap.GeneralPurposeAllocator(.{}){};
	defer std.debug.assert(!gpa.deinit());
	
	var traceable = try RetraceableContractionSequence(u32).init(gpa.allocator(),10,5);
	defer traceable.deinit(gpa.allocator());

	try std.testing.expectEqual(traceable.seq.list.len,traceable.tww.len);

	try traceable.addContraction(gpa.allocator(),8,7,100);

	try expectEqual(traceable.getTwinWidth(),100);
	try expectEqual(traceable.lastContraction(),.{.erased=8,.survivor=7});

	try traceable.removeLast();
	try expectEqual(traceable.getTwinWidth(),0);
	try expectEqual(traceable.seq.write_ptr,0);
	try expectEqual(traceable.lastContraction(),null);
}


test "RetraceableContractionSequence: Overflow" {
	var gpa = std.heap.GeneralPurposeAllocator(.{}){};
	defer std.debug.assert(!gpa.deinit());
	
	var traceable = try RetraceableContractionSequence(u32).init(gpa.allocator(),2,5);
	defer traceable.deinit(gpa.allocator());

	try std.testing.expectEqual(traceable.seq.list.len,traceable.tww.len);

	try traceable.addContraction(gpa.allocator(),0,1,90);
	try std.testing.expectError(contraction.ContractionError.ContractionOverflow,traceable.addContraction(gpa.allocator(),1,0,80));

	try expectEqual(traceable.isComplete(), true);
}
