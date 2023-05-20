const std = @import("std");

pub const BenchmarkError = error{ TimerAlreadyStarted, TimerNotStarted };

pub const BenchmarkHelper = struct {
    const Self = @This();
    total_time: u64,
    start_time: ?std.time.Instant,
    pub fn init() Self {
        return BenchmarkHelper{ .total_time = 0, .start_time = null };
    }

    pub fn start(self: *Self) !void {
        if (self.start_time != null) {
            return BenchmarkError.TimerAlreadyStarted;
        }
        self.start_time = try std.time.Instant.now();
    }

    pub fn stop(self: *Self) !void {
        if (self.start_time) |start_tp| {
            const now = try std.time.Instant.now();
            self.total_time += now.since(start_tp);
            self.start_time = null;
            return;
        }
        return BenchmarkError.TimerNotStarted;
    }
};

test "BenchmarkHelper: Check start and end error" {
    var bh = BenchmarkHelper.init();

    var count = 0;
    bh.start();
    try std.testing.expectError(bh.start(), BenchmarkError.TimerAlreadyStarted);
    for (0..1000) |i| {
        count += i;
    }
    bh.stop();
    try std.testing.expectError(bh.stop(), BenchmarkError.TimerNotStarted);
}
