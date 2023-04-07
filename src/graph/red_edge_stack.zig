const std = @import("std");

pub const NewRedEdge = packed struct {
	target: u31,
	turned_or_overwrite: bool,
	pub inline fn new(target: u32) NewRedEdge {
		return NewRedEdge {
			.target = @intCast(u31,target),
			.turned_or_overwrite = false
		};
	}
	pub inline fn newOverwritten(target: u32) NewRedEdge {
		return NewRedEdge {
			.target = @intCast(u31,target),
			.turned_or_overwrite = true 
		};
	}
};


pub const RedEdgeStackError = error {
	StackNotSealed,
	StackNoLevelLeft,
	StackOverflow
};

pub const RedEdgeStack = struct {
	const Self = @This();
	edge_stack: std.ArrayListUnmanaged(NewRedEdge),
	size: std.ArrayListUnmanaged(u32),
	allocator: std.mem.Allocator,
	pub fn init(allocator: std.mem.Allocator) !RedEdgeStack {
		var size = std.ArrayListUnmanaged(u32){};
		try size.append(allocator,0);

		return RedEdgeStack {
			.edge_stack = std.ArrayListUnmanaged(NewRedEdge){},
			.size = size,
			.allocator = allocator
		};
	}

	pub fn initCapacity(allocator: std.mem.Allocator, capacity: u32) !RedEdgeStack {
		std.debug.assert(capacity>0);
		var size = std.ArrayListUnmanaged(u32){};
		try size.append(allocator,0);

		return RedEdgeStack {
			.edge_stack = try std.ArrayListUnmanaged(NewRedEdge).initCapacity(allocator,capacity),
			.size = size,
			.allocator = allocator
		};
	}

	pub fn deinit(self: *Self) void {
		self.edge_stack.deinit(self.allocator);
		self.size.deinit(self.allocator);
	}

	pub const RedEdgeStackIterator = struct {
		index: u32,
		length: u32,
		stack: *const RedEdgeStack,
		pub inline fn next(self: *RedEdgeStackIterator) ?NewRedEdge {
			if(self.index >= self.length) {
				return null;
			}

			const item = self.stack.edge_stack.items[self.index];
			self.index+=1;
			return item;
		}

	};

	pub fn iterateLastLevel(self: *const Self) RedEdgeStackError!RedEdgeStackIterator {
		if(self.size.items[self.size.items.len-1] != self.edge_stack.items.len) {
			return RedEdgeStackError.StackNotSealed;
		}
		else if(self.size.items.len < 2) {
			return RedEdgeStackError.StackNoLevelLeft;
		}
		return RedEdgeStackIterator {
			.index = self.size.items[self.size.items.len-2],
			.length = @intCast(u32,self.edge_stack.items.len),
			.stack = self
		};
	}

	pub fn revertLastContraction(self: *Self) RedEdgeStackError!void {
		if(self.size.items.len >= 2) {
			// Reset to start of this level
			_ = self.size.pop();
			const resize = self.size.items[self.size.items.len-1];
			self.edge_stack.shrinkRetainingCapacity(resize);
		}
		else {
			return RedEdgeStackError.StackNoLevelLeft;
		}
	}


	//TODO: Add exception here!
	pub fn addEdge(self: *Self, edge: NewRedEdge) !void {
		try self.edge_stack.append(self.allocator,edge);
	}

	pub fn sealLevel(self: *Self) !void {
		try self.size.append(self.allocator,@intCast(u32,self.edge_stack.items.len));
	}
};


test "RedEdgeStack: New red edges" {
	{
		const edge = NewRedEdge.new(100);
		try std.testing.expectEqual(edge.target,100);
		try std.testing.expectEqual(edge.turned_or_overwrite, false);
	}


	{
		const edge = NewRedEdge.newOverwritten(102);
		try std.testing.expectEqual(edge.target,102);
		try std.testing.expectEqual(edge.turned_or_overwrite, true);
	}
}

test "RedEdgeStack: Red edge stack add and iterate" {	
	var gpa = std.heap.GeneralPurposeAllocator(.{}){};
	defer {
		std.debug.assert(!gpa.deinit());
	}
	var stack = try RedEdgeStack.init(gpa.allocator());
	defer stack.deinit();

	stack.addEdge(NewRedEdge.new(1)) catch unreachable;
	stack.addEdge(NewRedEdge.newOverwritten(2)) catch unreachable;
	try stack.sealLevel();
	stack.addEdge(NewRedEdge.new(3)) catch unreachable;
	stack.addEdge(NewRedEdge.newOverwritten(4)) catch unreachable;
	try stack.sealLevel();

	var iter = stack.iterateLastLevel() catch unreachable;
	var next = iter.next().?;
	var next_2 = iter.next().?;

	try std.testing.expectEqual(next.target, 3);
	try std.testing.expectEqual(next.turned_or_overwrite, false);

	try std.testing.expectEqual(next_2.target, 4);
	try std.testing.expectEqual(next_2.turned_or_overwrite, true);

	try stack.revertLastContraction();

	var iter_level_2 = stack.iterateLastLevel() catch unreachable;
	var next_level2 = iter_level_2.next().?;
	var next_level2_2 = iter_level_2.next().?;

	try std.testing.expectEqual(next_level2.target, 1);
	try std.testing.expectEqual(next_level2.turned_or_overwrite, false);

	try std.testing.expectEqual(next_level2_2.target, 2);
	try std.testing.expectEqual(next_level2_2.turned_or_overwrite, true);


	try stack.revertLastContraction();

	// Nothing left to iterate
	try std.testing.expectError(error.StackNoLevelLeft,stack.iterateLastLevel());
}

test "RedEdgeStack: Red edge stack add and empty stack" {	
	var gpa = std.heap.GeneralPurposeAllocator(.{}){};
	defer {
		std.debug.assert(!gpa.deinit());
	}
	var stack = try RedEdgeStack.init(gpa.allocator());
	defer stack.deinit();
	stack.addEdge(NewRedEdge.new(1)) catch unreachable;

	try stack.sealLevel();
	try stack.sealLevel();

	var iterator = try stack.iterateLastLevel();
	try std.testing.expectEqual(iterator.next(),null);

	try stack.revertLastContraction();
	var iterator_now = try stack.iterateLastLevel();
	const item = iterator_now.next().?;
	try std.testing.expectEqual(item.target, 1);
	try std.testing.expectEqual(item.turned_or_overwrite, false);
	try std.testing.expectEqual(iterator_now.next(),null);
}


test "RedEdgeStack: Red edge stack revert and add again" {	
	var gpa = std.heap.GeneralPurposeAllocator(.{}){};
	defer {
		std.debug.assert(!gpa.deinit());
	}
	var stack = try RedEdgeStack.init(gpa.allocator());
	defer stack.deinit();
	stack.addEdge(NewRedEdge.new(1)) catch unreachable;

	try stack.sealLevel();
	stack.addEdge(NewRedEdge.new(2)) catch unreachable;
	try stack.sealLevel();

	{
		var iterator = try stack.iterateLastLevel();
		const item_now = iterator.next().?;
		try std.testing.expectEqual(item_now.target,2);
		try std.testing.expectEqual(item_now.turned_or_overwrite,false);
		try std.testing.expectEqual(iterator.next(),null);
	}

	try stack.revertLastContraction();
	stack.addEdge(NewRedEdge.newOverwritten(3)) catch unreachable;
	try stack.sealLevel();

	{
		var iterator = try stack.iterateLastLevel();
		const item_now = iterator.next().?;
		try std.testing.expectEqual(item_now.target,3);
		try std.testing.expectEqual(item_now.turned_or_overwrite,true);

		//Next item must be null
		try std.testing.expectEqual(iterator.next(),null);

	}
	try stack.revertLastContraction();
	{
		var iterator = try stack.iterateLastLevel();
		const item_now = iterator.next().?;
		try std.testing.expectEqual(item_now.target,1);
		try std.testing.expectEqual(item_now.turned_or_overwrite,false);
		try std.testing.expectEqual(iterator.next(),null);
	}
	try stack.revertLastContraction();
	try std.testing.expectError(error.StackNoLevelLeft,stack.iterateLastLevel());
	try std.testing.expectError(error.StackNoLevelLeft,stack.revertLastContraction());
}

test "RedEdgeStack: Red edge stack check safety guards" {	
	var gpa = std.heap.GeneralPurposeAllocator(.{}){};
	defer {
		std.debug.assert(!gpa.deinit());
	}
	var stack = try RedEdgeStack.init(gpa.allocator());
	defer stack.deinit();
	stack.addEdge(NewRedEdge.new(1)) catch unreachable;
	try std.testing.expectError(error.StackNotSealed,stack.iterateLastLevel());
}

test "RedEdgeStack: Check capacity" {	
	// Space for 3 items 
	var gpa = std.heap.GeneralPurposeAllocator(.{}){};
	defer {
		std.debug.assert(!gpa.deinit());
	}
	var stack = try RedEdgeStack.initCapacity(gpa.allocator(),100);
	defer stack.deinit();

	try std.testing.expectEqual(stack.edge_stack.capacity, 100);
}
