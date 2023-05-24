#!/usr/bin/env bash
zig build -Doptimize=ReleaseFast -Dcpu=haswell
strip zig-out/bin/solver
strip zig-out/bin/solver_exact
