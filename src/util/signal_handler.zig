const std = @import("std");
const os = std.os.linux;
const Sigaction = os.Sigaction;
const sigaction = os.sigaction;
const empty_sigset = os.empty_sigset;
const SIGTERM: u6 = 15;
const BUF_SIZE = 32768;

var cc_slice: [][]u32  = undefined;

// Writes the current solution of cc_slice into stdout. errors all all discarded, due to callconv(.C)
// requirement they can not be propagated as a return type (!void).
// any better idea is appreciated.
fn handle_sigterm(_: c_int) callconv(.C) void {
    var file = std.io.getStdOut();
    var buffered = std.io.BufferedWriter(BUF_SIZE, @TypeOf(file.writer())){
        .unbuffered_writer = file.writer()
    };
    var writer = buffered.writer();

    const n = cc_slice.len;
    var i: usize = 0;
    while (i<n) {
        const m = cc_slice[i].len;
        var j: usize = 0;
        while (j < m-1) {
            std.fmt.format(writer,"{d} {d}\n",.{cc_slice[i][j]+1,cc_slice[i][j+1]+1}) catch {};
            j+=2;
        }
        i+=1;
    }
    buffered.flush() catch {};
    os.exit(0);
}

const act = Sigaction {
    .handler = .{. handler = handle_sigterm},
    .mask = empty_sigset,
    .flags = 0,
};

// initializes the signal handler. it expects a slice of slices of u32s.
// each slice S represents the solution of a connected component.
// each slice S for a connected component is expected to have a pair of contractions stored at S[i..i+1]
// for all i < S.len - 1 and i%2==0.
// vertices are expected to start at 0 and end at n-1.
pub fn initialize_signal_handler(init_cc_slice: [][]u32) void {
    cc_slice = init_cc_slice;
    _ = sigaction(SIGTERM, &act, null);
}

