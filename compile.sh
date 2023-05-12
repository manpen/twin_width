#!/usr/bin/env bash
zig build solver -Doptimize=ReleaseFast -Dcpu=haswell
zig build solver_exact -Doptimize=ReleaseFast -Dcpu=haswell
strip zig-out/bin/solver
strip zig-out/bin/solver_exact
