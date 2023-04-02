const std = @import("std");
const builtin = @import("builtin");
comptime { _ = @import("main.zig"); }


pub fn main() !void {
    for (builtin.test_functions) |test_fn| {
				if(std.mem.eql(u8,test_fn.name[0..11],"test.bench:")) {
					var calc:u64 = 0;
					var i: u32 = 0;
					while(i<10):(i+=1) {
						var timer = try std.time.Instant.now();
						try test_fn.func();
						var end = try std.time.Instant.now();
						calc+=end.since(timer);
					}
					calc = calc/10;
					if(calc>=1000_000) {
						std.debug.print("bench: {s} took {d}ms\n",.{test_fn.name,calc/1000_000});
					}
					else if(calc>=1000) {
						std.debug.print("bench: {s} took {d}us\n",.{test_fn.name,calc/1000});
					}
					else {
						std.debug.print("bench: {s} took {d}ns\n",.{test_fn.name,calc});
					}
				}
    }
}
