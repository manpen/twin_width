const std = @import("std");

// Although this function looks imperative, note that its job is to
// declaratively construct a build graph that will be executed by an external
// runner.
pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});

    const sub = b.addExecutable(.{
        .name = "solver_heuristic",
        .root_source_file = .{ .path = "src/submission.zig" },
        .target = target,
        .optimize = optimize,
    });
    sub.single_threaded = true;
    sub.addIncludePath("src/tww");
    b.installArtifact(sub);

    const run_sub = b.addRunArtifact(sub);
    run_sub.step.dependOn(b.getInstallStep());
    const solver_step = b.step("solver_heuristic", "Compile heuristic solver for submission");
    solver_step.dependOn(b.getInstallStep());

    /////////////////////////////

    const sub_exact = b.addExecutable(.{
        .name = "solver_exact",
        .root_source_file = .{ .path = "src/submission_exact.zig" },
        .target = target,
        .optimize = optimize,
    });
    sub_exact.single_threaded = true;
    sub_exact.addIncludePath("src/tww");
    b.installArtifact(sub_exact);

    const run_sub_exact = b.addRunArtifact(sub_exact);
    run_sub_exact.step.dependOn(b.getInstallStep());
    const solver_step_exact = b.step("solver_exact", "Compile exact solver for submission");
    solver_step_exact.dependOn(b.getInstallStep());
    
}
