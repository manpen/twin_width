# Twin Width Solver

This repository contains exact and heuristic twin width solvers that were developed for the [PACE 2023](https://pacechallenge.org/2023/) challenge.

The solvers are implemented in the Zig programming language. Since this language (and its standard library) are still rapidly evolving, we highly recommend to compile the code with Zig version `0.11.0-dev.3045`:

- [64-bit linux](https://ziglang.org/builds/zig-linux-x86_64-0.11.0-dev.3045+526065723.tar.xz)
- [64-bit windows](https://ziglang.org/builds/zig-windows-x86_64-0.11.0-dev.3045+526065723.zip)

The only explicit dependency is the Zig compiler.
After the installation of the compiler ---([essentially unpacking of the packages linked before](https://ziglang.org/learn/getting-started/#installing-zig))--- run the following command to build the solvers and expect a build time of roughly a minute.

```
zig build -Doptimize=ReleaseFast -Dcpu=native
```

If everything works fine (if not please check the installed Zig version!) the binaries are generated in the directory `zig-out/bin`.
They use the same semantics as on the [optil.io platform](https://www.optil.io/optilion/problem/3205), i.e.:

- The problem is piped in via STDIN
- The solution is emitted via STDOUT
- Additional information are printed via STDERR
- If the program panics (non-zero return code), the output shall be disregarded
- The exact solver terminates once optimality of the solution was established
- The heuristic solver keeps running until it receives a SIGTERM and then dumps the best solution found
