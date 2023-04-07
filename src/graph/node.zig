const std = @import("std");
const edge_list = @import("edge_list.zig");
const comptime_util = @import("../util/comptime_checks.zig");



// Id is omitted must be stored outside of the struct
pub fn Node(comptime T: type) type {
	comptime if (!comptime_util.checkIfIsCompatibleInteger(T)) {
		@compileError("T must be either u8,u16 or u32");
	};
	
	return struct {
		black_edges: edge_list.ParametrizedSortedArrayList(T),
		red_edges: edge_list.ParametrizedUnsortedArrayList(T)
	};
}
