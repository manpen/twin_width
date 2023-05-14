const std = @import("std");
const mg = @import("matrix_graph.zig");
const bs = @import("bootstrapping.zig");
const subgraph = @import("../graph/subgraph.zig");
const cc = @import("../graph/connected_component.zig");
const contr = @import("../tww/contraction_sequence.zig");

const assert = std.debug.assert;

const Node = u32;
const Contraction = struct { rem: Node, sur: Node };

pub const SolverError = error{
    Infeasable,
    Timeout,
};

const ContraScore = struct { score: Node, distTwo: bool, rem: Node, sur: Node };
const CandidateListRef = struct {
    num_nodes: u32,
    candidates: []ContraScore,
    invalid: *void,
};

const FeatureReuseCandidates: bool = true;
const FeatureUseInfeasibleCache: bool = true;
const FeatureShrinkGraph: bool = true;
const FeatureLeafNodePruning: bool = true;
const FeatureTwinPruning: bool = true;

pub const ExactBranchAndBound = struct {
    const Self = @This();
    const InfeasibleCache = std.AutoHashMap(mg.Digest, void);

    allocator: std.mem.Allocator,
    lower: Node,
    upper_excl: Node,
    number_of_nodes: Node,

    working_seq: std.ArrayList(Contraction),
    best_seq: std.ArrayList(Contraction),

    num_calls: u64,
    recursion_depth: usize,

    smallest_depth_last_report: usize,
    time_last_report: std.time.Instant,

    infeasible_cache: InfeasibleCache,
    candidate_lists: std.ArrayList(CandidateListRef),

    pub fn new(allocator: std.mem.Allocator, number_of_nodes: Node) !Self {
        var working_seq = try std.ArrayList(Contraction).initCapacity(allocator, number_of_nodes);
        errdefer working_seq.deinit();

        var best_seq = try std.ArrayList(Contraction).initCapacity(allocator, number_of_nodes);
        errdefer best_seq.deinit();

        var candidate_lists = try std.ArrayList(CandidateListRef).initCapacity(allocator, @boolToInt(FeatureReuseCandidates) * number_of_nodes);
        errdefer candidate_lists.deinit();

        var cache = InfeasibleCache.init(allocator);
        errdefer cache.deinit();
        if (FeatureUseInfeasibleCache) {
            try cache.ensureTotalCapacity((2 << 30) / @sizeOf(mg.Digest));
        }

        var solver = Self{
            .allocator = allocator, //
            .lower = 0,
            .upper_excl = number_of_nodes,
            .number_of_nodes = number_of_nodes,
            .working_seq = working_seq,
            .best_seq = best_seq,
            .num_calls = 0,
            .recursion_depth = 0,
            .time_last_report = try std.time.Instant.now(),
            .smallest_depth_last_report = 0,
            .infeasible_cache = cache,
            .candidate_lists = candidate_lists,
        };

        return solver;
    }

    pub fn deinit(self: *Self) void {
        self.infeasible_cache.deinit();
        self.best_seq.deinit();
        self.working_seq.deinit();
    }

    pub fn setLowerBound(self: *Self, lb: Node) void {
        self.lower = lb;
    }

    /// Sets exclusive upper bound, i.e. a solution of size `ub - 1` is acceptable.
    pub fn setUpperBound(self: *Self, ub: Node) void {
        std.debug.assert(ub > 0);
        self.upper_excl = ub;
    }

    pub fn numberOfRecursiveCalls(self: *Self) u64 {
        return self.num_calls;
    }

    fn shrinkGraph(self: *Self, comptime Graph: type, graph: *const Graph, slack: Node) anyerror!Node {
        const Smaller = mg.MatrixGraph(Graph.NumNodes / 2);

        var newNodeIdOfOld = try self.allocator.alloc(u32, Graph.NumNodes);
        defer self.allocator.free(newNodeIdOfOld);

        var oldNodeIdOfNew = try self.allocator.alloc(u32, Smaller.NumNodes);
        defer self.allocator.free(oldNodeIdOfNew);

        {
            var u: u32 = 0;
            var i: u32 = 0;
            while (i < Graph.NumNodes) : (i += 1) {
                newNodeIdOfOld[i] = u;
                if (graph.has_neighbors.isSet(i)) {
                    oldNodeIdOfNew[u] = i;
                    u += 1;
                }
            }
        }

        var smaller_storage = try self.allocator.alloc(Smaller, 1);
        defer self.allocator.free(smaller_storage);

        var smaller = &smaller_storage[0];
        smaller.* = Smaller.new();

        // copy edges; TODO: this can be implement with parallel bit shifts
        {
            smaller.beginBatchUpdates();
            defer smaller.endBatchUpdates();

            var iter = graph.has_neighbors.iter_set();
            while (iter.next()) |u| {
                var mapped_u = newNodeIdOfOld[u];
                // neigh
                var niter = graph.constNeighbors(u).iter_set();
                while (niter.next()) |v| {
                    if (u < v) continue;
                    var color = if (graph.constRedNeighbors(u).isSet(v)) mg.Color.Red else mg.Color.Black;
                    _ = smaller.addEdge(mapped_u, newNodeIdOfOld[v], color);
                }
            }
        }

        const old_contractions = self.working_seq.items.len;
        const upper_before = self.upper_excl;
        var result = try self.solve(Smaller, smaller, slack);
        std.debug.assert(self.upper_excl < upper_before);

        if (false) {
            var i: u32 = 0;
            for (self.best_seq.items) |ctr| {
                std.debug.print("{d}: {d}-{d}\n", .{ i, ctr.rem, ctr.sur });
                i += 1;
            }

            std.debug.print("G: {d} S: {d} old: {d}\n", .{ Graph.NumNodes, Smaller.NumNodes, old_contractions });
        }

        for (self.best_seq.items[old_contractions..]) |*ctr| {
            ctr.rem = oldNodeIdOfNew[ctr.rem];
            ctr.sur = oldNodeIdOfNew[ctr.sur];
        }

        return result;
    }

    pub fn solve(self: *Self, comptime Graph: type, graph: *const Graph, slack: Node) anyerror!Node {
        if (FeatureShrinkGraph and Graph.NumNodes > 8 and 2 * graph.has_neighbors.cardinality() < Graph.NumNodes) {
            return self.shrinkGraph(Graph, graph, slack);
        }

        self.num_calls += 1;
        self.recursion_depth += 1;
        defer self.smallest_depth_last_report = @min(self.smallest_depth_last_report, self.recursion_depth);
        defer self.recursion_depth -= 1;

        if (self.num_calls == 1) {
            var local_graph = try self.allocator.alloc(Graph, 1);
            defer self.allocator.free(local_graph);
            local_graph[0] = graph.*;

            try self.initial_kernelization(Graph, &local_graph[0], slack);
            return self.solve(Graph, &local_graph[0], slack);
        }

        if (self.num_calls % 100000 == 0) {
            var now = try std.time.Instant.now();
            var elapsed = now.since(self.time_last_report);
            self.time_last_report = now;

            std.debug.print("Call: {d} Depth: {d} (smallest: {d}) Graph({d}) n={d} slack={d} upper={d} time={d}ms \n", .{ self.num_calls, self.recursion_depth, self.smallest_depth_last_report, Graph.NumNodes, graph.has_neighbors.cardinality(), slack, self.upper_excl, elapsed / 1000_000 });
            self.smallest_depth_last_report = self.recursion_depth;
        }

        var digest: mg.Digest = undefined;
        if (FeatureUseInfeasibleCache) {
            digest = graph.hash();
            if (self.infeasible_cache.contains(digest)) {
                return SolverError.Infeasable;
            }
        }

        var frame = &(try self.allocator.alloc(Frame(Graph), 1))[0];
        defer self.allocator.free(@ptrCast(*[1]Frame(Graph), frame));

        frame.* = try Frame(Graph).new(self, graph, slack);
        defer frame.deinit();
        var result = frame.run();

        if (FeatureUseInfeasibleCache and result == SolverError.Infeasable) {
            try self.infeasible_cache.put(digest, {});
        }

        return result;
    }

    fn initial_kernelization(self: *Self, comptime Graph: type, graph: *Graph, slack: Node) !void {
        _ = slack;
        while (FeatureLeafNodePruning and try self.eliminate_leafs(Graph, graph)) {}
    }

    fn eliminate_leafs(
        self: *Self,
        comptime Graph: type,
        graph: *Graph,
    ) !bool {
        var leafs = @TypeOf(graph.*).BitSet.new();
        var u: Node = 0;
        while (u < graph.numberOfNodes()) : (u += 1) {
            if (graph.deg(u) == 1) {
                _ = leafs.setBit(u);
            }
        }

        if (leafs.cardinality() < 2) {
            return false;
        }

        u = 0;
        var change = false;
        while (u < graph.numberOfNodes()) : (u += 1) {
            if (graph.deg(u) < 2) {
                continue;
            }
            var leafs_at_node = leafs.copyWithAnd(graph.constNeighbors(u));

            if (leafs_at_node.cardinality() < 2) {
                continue;
            }

            var leafs_iter = leafs_at_node.iter_set();
            const survivor = leafs_iter.next().?;
            while (leafs_iter.next()) |rem| {
                try self.working_seq.append(Contraction{ .rem = rem, .sur = survivor });
                graph.mergeNodes(rem, survivor, null); // TODO: this can be made faster by removing the Edge
                change = true;
            }
        }

        return change;
    }
};

fn Frame(comptime Graph: type) type {
    return struct {
        const Self = @This();

        context: *ExactBranchAndBound,
        candidates: std.ArrayList(ContraScore),

        input_graph: *const Graph,

        work_graph: *Graph,
        work_tww: Node,
        work_contracted: *Graph.BitSet,
        slack: Node,

        fn new(solver: *ExactBranchAndBound, graph: *const Graph, slack: Node) !Self {
            var work_graph = try solver.allocator.alloc(Graph, 1);
            errdefer solver.allocator.free(work_graph);

            var work_contracted = try solver.allocator.alloc(Graph.BitSet, 1);
            errdefer solver.allocator.free(work_contracted);

            var n = @max(4, @as(usize, graph.has_neighbors.cardinality()));

            var candidates = try std.ArrayList(ContraScore).initCapacity(solver.allocator, (n + 1) / 2 * n);
            errdefer candidates.deinit();

            return Self{
                .context = solver,

                .input_graph = graph,
                .candidates = candidates,

                .work_graph = &work_graph[0],
                .work_contracted = &work_contracted[0],
                .slack = slack,
                .work_tww = slack,
            };
        }

        fn deinit(self: *Self) void {
            self.candidates.deinit();
            self.context.allocator.free(@ptrCast(*[1]Graph.BitSet, self.work_contracted));
            self.context.allocator.free(@ptrCast(*[1]Graph, self.work_graph));
        }

        fn run(self: *Self) !Node {
            if (false) {
                for (self.context.working_seq.items) |x| {
                    _ = x;
                    std.debug.print(" ", .{});
                }
                std.debug.print("LB: {d} Slack: {d} UB: {d} n={d} d={d}\n", .{ self.context.lower, self.slack, self.context.upper_excl, self.input_graph.has_neighbors.cardinality(), self.context.candidate_lists.items.len });
            }

            // It's wasteful to reach this point with a slack that is too large. We should have pruned it earlier.
            std.debug.assert(self.slack < self.context.upper_excl);

            if (self.input_graph.numberOfEdges() == 0) {
                std.debug.print("Found a solution with tww={d}\n", .{self.slack});

                self.context.upper_excl = self.slack;
                self.context.best_seq.clearRetainingCapacity();
                try self.context.best_seq.appendSlice(self.context.working_seq.items);

                return self.slack;
            }

            const initial_seq_len = self.context.working_seq.items.len;

            self.computeCandidates();

            var result: SolverError!Node = SolverError.Infeasable;

            for (self.candidates.items) |c| {
                try self.context.working_seq.resize(initial_seq_len);

                self.work_graph.* = self.input_graph.*;
                self.work_tww = self.slack;
                self.work_contracted.unsetAll();

                try self.contractNodes(c.rem, c.sur, !c.distTwo);
                if (self.work_tww >= self.context.upper_excl) {
                    continue;
                }

                if (FeatureReuseCandidates) {
                    self.context.candidate_lists.appendAssumeCapacity(CandidateListRef{
                        .num_nodes = Graph.NumNodes, //
                        .candidates = self.candidates.items[0..],
                        .invalid = @ptrCast(*void, &self.candidates),
                    });
                }
                defer if (FeatureReuseCandidates) {
                    _ = self.context.candidate_lists.pop();
                };
                const tww_from_recursion = self.context.solve(Graph, self.work_graph, self.work_tww) catch continue;

                std.debug.assert(tww_from_recursion >= self.work_tww);

                result = tww_from_recursion;

                if (tww_from_recursion <= @max(self.slack, self.context.lower)) {
                    break;
                }
            }

            return result;
        }

        fn contractNodes(self: *Self, rem: Node, sur: Node, atleast_dist_three: bool) !void {
            if (atleast_dist_three) {
                self.work_graph.mergeDistantNodes(rem, sur);
            } else {
                self.work_graph.*.mergeNodes(rem, sur, null);
            }

            self.work_tww = @max(self.work_tww, self.work_graph.redDeg(sur));
            self.work_tww = @max(self.work_tww, self.work_graph.maxRedDegreeIn(self.work_graph.constRedNeighbors(sur)));
            self.work_contracted.assignOr(self.work_graph.twoNeighbors(sur));

            if (self.work_tww >= self.context.upper_excl) {
                return;
            }

            try self.context.working_seq.append(Contraction{ .rem = rem, .sur = sur });

            if (!atleast_dist_three) {
                _ = try self.eliminate_leafs_at(sur);
            }

            _ = try self.eliminate_twins_at(sur);
        }

        fn computeCandidates(self: *Self) void {
            const depth = self.context.candidate_lists.items.len;
            if (!FeatureReuseCandidates or depth == 0 or self.context.candidate_lists.items[depth - 1].num_nodes != Graph.NumNodes) {
                return self.computeCandidatesFromScratch();
            }

            var prev: *const CandidateListRef = &self.context.candidate_lists.items[depth - 1];
            var invalid = @ptrCast(*Graph.BitSet, @alignCast(@alignOf(Graph.BitSet), prev.invalid));

            for (prev.candidates) |c| {
                const rem = c.rem;
                const sur = c.sur;

                if (!self.input_graph.has_neighbors.isSet(rem) or !self.input_graph.has_neighbors.isSet(sur)) {
                    continue;
                }

                if (invalid.isSet(rem) or invalid.isSet(sur)) {
                    self.candidates.appendAssumeCapacity(c);
                    continue;
                }

                if (self.input_graph.constTwoNeighbors(rem).isSet(sur)) {
                    var reds = self.input_graph.redNeighborsAfterMerge(rem, sur);
                    var red_degree = reds.cardinality();

                    if (!FeatureReuseCandidates and red_degree >= self.context.upper_excl) {
                        continue;
                    }

                    self.candidates.appendAssumeCapacity(ContraScore{ .score = red_degree, .distTwo = true, .rem = rem, .sur = sur });
                } else {
                    var red_degree = self.input_graph.deg(rem) + self.input_graph.deg(sur);

                    if (!FeatureReuseCandidates and red_degree >= self.context.upper_excl) {
                        continue;
                    }

                    self.candidates.appendAssumeCapacity(ContraScore{ .score = red_degree, .distTwo = false, .rem = rem, .sur = sur });
                }
            }

            std.sort.sort(ContraScore, self.candidates.items, {}, cmpByValue);
        }

        fn computeCandidatesFromScratch(self: *Self) void {
            var nodes = self.input_graph.has_neighbors;
            var u: Node = 0;

            while (u < self.input_graph.numberOfNodes()) : (u += 1) {
                if (!nodes.unsetBit(u)) {
                    std.debug.assert(self.input_graph.deg(u) == 0);
                    continue;
                }
                std.debug.assert(self.input_graph.deg(u) > 0);

                var two_neighbors = self.input_graph.constTwoNeighbors(u);
                var deg_u = self.input_graph.deg(u);

                var iter = nodes.iter_set();
                while (iter.next()) |v| {
                    std.debug.assert(u < v);

                    if (two_neighbors.isSet(v)) {
                        var reds = self.input_graph.redNeighborsAfterMerge(u, v);
                        var red_degree = reds.cardinality();

                        // simplifies reusing of the candidate list
                        //if (red_degree >= self.context.upper_excl) {
                        //    continue;
                        //}

                        self.candidates.appendAssumeCapacity(ContraScore{ .score = red_degree, .distTwo = true, .rem = v, .sur = u });
                    } else {
                        var deg_v = self.input_graph.deg(v);

                        var red_degree = deg_u + deg_v;

                        // simplifies reusing of the candidate list
                        //if (red_degree >= self.context.upper_excl) {
                        //    continue;
                        //}

                        self.candidates.appendAssumeCapacity(ContraScore{ .score = red_degree, .distTwo = false, .rem = v, .sur = u });
                    }
                }
            }

            std.sort.sort(ContraScore, self.candidates.items, {}, cmpByValue);
        }

        noinline fn eliminate_leafs_at(self: *Self, u: Node) !bool {
            if (!FeatureLeafNodePruning) {
                return false;
            }

            var leafs = Graph.BitSet.new();
            var neighs = self.work_graph.constNeighbors(u).iter_set();
            while (neighs.next()) |v| {
                if (self.work_graph.deg(v) == 1) {
                    _ = leafs.setBit(v);
                }
            }

            if (leafs.cardinality() < 2) {
                return false;
            }

            var liter = leafs.iter_set();
            var sur = liter.next().?;

            var any_red: bool = false;
            while (liter.next()) |v| {
                any_red = any_red or self.work_graph.removeEdge(u, v).?.isRed();

                // no need to call contractNodes, since removing leafs does not warrant new checks
                try self.context.working_seq.append(Contraction{ .rem = v, .sur = sur });
            }

            if (any_red) {
                _ = self.work_graph.addEdge(u, sur, mg.Color.Red);
            }

            return true;
        }

        noinline fn eliminate_twins_at(self: *Self, u: Node) anyerror!bool {
            if (!FeatureTwinPruning) {
                return false;
            }

            var two_neighbors = self.work_graph.constTwoNeighbors(u);
            if (two_neighbors.cardinality() > 10) {
                return false;
            }

            var it1 = two_neighbors.iter_set();
            while (it1.next()) |v| {
                var it2 = it1;
                var deg_v = self.work_graph.deg(v);
                std.debug.assert(deg_v > 0);

                while (it2.next()) |w| {
                    std.debug.assert(v != w);
                    var a = w; // superset
                    var b = v; // subset
                    if (deg_v >= self.work_graph.deg(w)) {
                        a = v;
                        b = w;
                    }
                    std.debug.assert(self.work_graph.deg(w) > 0);

                    var prev = self.work_graph.addEdge(a, b, mg.Color.Black);

                    if (self.work_graph.constNeighbors(b).is_subset_of(self.work_graph.constNeighbors(a)) and
                        self.work_graph.constRedNeighbors(b).is_subset_of(self.work_graph.constRedNeighbors(a)))
                    {
                        var red = self.work_graph.redNeighborsAfterMerge(a, b);
                        std.debug.assert(red.is_equal(self.work_graph.constRedNeighbors(a)));

                        _ = try self.contractNodes(a, b, false);
                        return true;
                    }

                    if (prev) |p| {
                        if (p.isRed()) {
                            _ = self.work_graph.addEdge(a, b, mg.Color.Red);
                        }
                    } else {
                        _ = self.work_graph.removeEdge(a, b);
                    }
                }
            }

            return false;
        }

        fn cmpByValue(context: void, a: ContraScore, b: ContraScore) bool {
            _ = context;
            const a_score = (a.score << @boolToInt(!a.distTwo));
            const b_score = (b.score << @boolToInt(!b.distTwo));

            return (a_score < b_score) or (a_score == b_score and (a.rem < b.rem or (a.rem == b.rem and a.sur < b.sur)));
        }
    };
}

fn solveCCContext(comptime T: type) type {
    return struct {
        const IntType = T;

        cc: *cc.ConnectedComponent(T),
        allocator: std.mem.Allocator,
        lower: Node,
        upper: Node,
        result: anyerror!Node,
    };
}

pub fn solveCCExactly(comptime T: type, component: *cc.ConnectedComponent(T), allocator: std.mem.Allocator, graph: *bs.MatrixGraphFromInducedSubGraph, lower: Node, upper: Node) !Node {
    const ContextType = solveCCContext(T);
    var context = ContextType{ .cc = component, .allocator = allocator, .lower = lower, .upper = upper, .result = upper + 1 };

    graph.dispatch(&context, solveCCExactlyHandler);

    return context.result;
}

fn solveCCExactlyHandler(context: anytype, graph: anytype, mapping: anytype) void {
    const T = @TypeOf(context.*).IntType;

    std.debug.print("Invoke exact solver on n={d}, m={d}, lower={d}, upper={d}\n", .{ graph.has_neighbors.cardinality(), graph.numberOfEdges(), context.lower, context.upper });

    var solver = ExactBranchAndBound.new(context.allocator, graph.numberOfNodes()) catch |e| {
        context.result = e;
        return;
    };
    defer solver.deinit();

    solver.setLowerBound(context.lower);
    solver.setUpperBound(context.upper);

    var tww = solver.solve(@TypeOf(graph), &graph, 0) catch |e| {
        std.debug.print("Catch\n", .{});
        context.result = e;
        return;
    };

    context.result = tww;

    // update CC
    context.cc.tww = @intCast(T, tww);
    var cs = &context.cc.best_contraction_sequence;
    cs.reset();

    const OuterContraction = contr.Contraction(T);
    for (solver.best_seq.items) |ctr| {
        cs.addContraction(OuterContraction{ .erased = @intCast(T, mapping.items[ctr.rem]), .survivor = @intCast(T, mapping.items[ctr.sur]) }) catch unreachable;
    }
}

test "Init BB solver" {
    var allocator = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = allocator.deinit();
    var solver = try ExactBranchAndBound.new(allocator.allocator(), 100);
    defer solver.deinit();
}
