const std = @import("std");
const os = std.os.linux;
const Sigaction = os.Sigaction;
const sigaction = os.sigaction;
const empty_sigset = os.empty_sigset;
const SIGTERM: u6 = 15;
const BUF_SIZE = 32768;

var cc_slice: [][]u32 = undefined;

// Writes the current solution of cc_slice into stdout. all errors all discarded, due to callconv(.C)
// requirement they can not be propagated as a return type (!void).
// any better idea is appreciated.
fn handle_sigterm(_: c_int) callconv(.C) void {
    var file = std.io.getStdOut();
    var buffered = std.io.BufferedWriter(BUF_SIZE, @TypeOf(file.writer())){ .unbuffered_writer = file.writer() };
    var writer = buffered.writer();

    const n = cc_slice.len;
    for (0..n) |i| {
        const m = cc_slice[i].len;
        // isolated vertex/empty slice
        if (m < 2) {
            continue;
        }
        var j: usize = 0;
        while (j < m - 1) {
            std.fmt.format(writer, "{d} {d}\n", .{ cc_slice[i][j] + 1, cc_slice[i][j + 1] + 1 }) catch {};
            j += 2;
        }
    }
    // each cc must be contracted with each other cc.
    for (0..(n - 1)) |i| {
        var j = cc_slice[i].len;
        var k = cc_slice[i + 1].len;
        if (j == 0 or k == 0) {
            continue;
        }
        if (j == 1) {
            j += 1;
        }
        if (k == 1) {
            k += 1;
        }
        std.fmt.format(writer, "{d} {d}\n", .{ cc_slice[i + 1][k - 2] + 1, cc_slice[i][j - 2] + 1 }) catch {};
    }
    buffered.flush() catch {};
    os.exit(0);
}

const act = Sigaction{
    .handler = .{ .handler = handle_sigterm },
    .mask = empty_sigset,
    .flags = 0,
};

// initializes the signal handler. it expects a slice of slices of u32s.
// each slice S represents the solution of a connected component.
// each slice S for a connected component is expected to have a pair of contractions stored at S[i..i+1]
// for all i < S.len - 1 and i%2==0.
// the first entry in the slice S[i..i+1] is the vertex that remains after the contraction.
// vertices are expected to start at 0 and end at n-1.
// special case: a connected component consisting of a single isolated vertex.
// in this case, the respective slice contains only the single vertex.
pub fn initialize_signal_handler(init_cc_slice: [][]u32) void {
    cc_slice = init_cc_slice;
    _ = sigaction(SIGTERM, &act, null);
}
