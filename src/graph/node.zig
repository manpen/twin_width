const std = @import("std");
const comptime_util = @import("../util/comptime_checks.zig");
const compressed_bitmap = @import("../util/compressed_bitmap.zig");



// Id is omitted must be stored outside of the struct
pub fn Node(comptime T: type, comptime promote_threshold: u32, comptime degrade_threshold: u32) type {
	comptime if (!comptime_util.checkIfIsCompatibleInteger(T)) {
		@compileError("T must be either u8,u16 or u32");
	};
	
	return struct {
		black_edges: compressed_bitmap.FastCompressedBitmap(T,promote_threshold,degrade_threshold),
		red_edges: compressed_bitmap.FastCompressedBitmap(T,promote_threshold,degrade_threshold),
	};
}
