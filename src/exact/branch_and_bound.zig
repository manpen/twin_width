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
    UnknownError,
};

const ContraScore = struct {
    score: Node,
    redDeg: Node,
    rem: Node,
    sur: Node,
    processed: bool,
};

const CandidateListRef = struct {
    num_nodes: u32,
    candidates: []ContraScore,
    invalid: *void,
};

const INITIAL_CACHE_SIZE_MB: usize = 100;

const FeatureParanoia: bool = false;

const FeatureReuseCandidates: bool = true;
const FeatureParanoidCandidates: bool = FeatureReuseCandidates and false; // additional expensive assertions
const FeatureSkipProcessed: bool = FeatureReuseCandidates and true;

const FeatureShrinkGraph: bool = true;

const FeatureSkipInfeasibleCandidatesEarly: bool = true;
const FeatureSkipPathNodes: bool = true;
const FeatureCriticalCandidates: bool = true; // skips dist>3 pairs with critical black neighbor

const FeaturePessimiticCandidateAllocation: bool = false;

const FeatureUseInfeasibleCache: bool = true;
const FeatureUseSolutionCache: bool = false;
const FeatureUseWFHash: bool = false;

const FeatureReportProgress: u32 = 0 * 1000_000; // set to 0 to disable

const FeaturePruning: bool = true;

const FeatureInitialPruning: bool = FeaturePruning and true;
const FeatureLeafNodePruning: bool = FeaturePruning and true;
const FeatureTwinPruning: bool = FeaturePruning and true;
const FeatureTinyGraphBelowSlack: bool = FeaturePruning and true; // almost no effect, as it only helps, if we're about to find a solution
const FeaturePruneTinyRedBridges: bool = FeaturePruning and false;
const FeaturePruneGeneralizedTwins: bool = FeaturePruning and false;
const FeaturePruneBlackCCs: bool = FeaturePruning and false;

const FeatureComplements: bool = false; // BROKEN
const FeatureOuterPathPruning: bool = FeaturePruning and false;

const FeatureTryLBFirst: bool = true;
const BudgetSubgraphLB: u64 = 20_000; // ms

const SolSummary = struct {
    tww: Node,
    slack: Node,
};

pub const ExactBranchAndBound = struct {
    const Self = @This();
    const InfeasibleCache = std.AutoHashMap(mg.Digest, void);
    const SolutionCache = std.AutoHashMap(mg.Digest, SolSummary);

    allocator: std.mem.Allocator,
    lower: Node,
    upper_excl: Node,
    number_of_nodes: Node,

    working_seq: std.ArrayList(Contraction),
    best_seq: std.ArrayList(Contraction),

    num_calls: u64,
    num_cache_hits: u64,
    recursion_depth: usize,

    smallest_depth_last_report: usize,
    time_last_report: std.time.Instant,

    infeasible_cache: InfeasibleCache,
    solution_cache: SolutionCache,
    candidate_lists: std.ArrayList(CandidateListRef),

    start_time: std.time.Instant,
    timeout_ms: ?u64,

    lower_bound_mode: bool,
    rng: ?*std.rand.DefaultPrng,
    silent: bool,

    single_pass: bool,

    pub fn new(allocator: std.mem.Allocator, number_of_nodes: Node) !Self {
        assert(number_of_nodes > 0);

        var working_seq = try std.ArrayList(Contraction).initCapacity(allocator, number_of_nodes - 1);
        errdefer working_seq.deinit();

        var best_seq = try std.ArrayList(Contraction).initCapacity(allocator, number_of_nodes - 1);
        errdefer best_seq.deinit();

        var candidate_lists = try std.ArrayList(CandidateListRef).initCapacity(allocator, @boolToInt(FeatureReuseCandidates) * (20 + number_of_nodes));
        errdefer candidate_lists.deinit();

        var infeasible_cache = InfeasibleCache.init(allocator);
        errdefer infeasible_cache.deinit();
        if (FeatureUseInfeasibleCache) {
            try infeasible_cache.ensureTotalCapacity((INITIAL_CACHE_SIZE_MB << 10) / @sizeOf(mg.Digest));
        }

        var solution_cache = SolutionCache.init(allocator);
        errdefer solution_cache.deinit();
        if (FeatureUseSolutionCache) {
            try solution_cache.ensureTotalCapacity((INITIAL_CACHE_SIZE_MB << 10) / @sizeOf(mg.Digest));
        }

        var now = try std.time.Instant.now();

        var solver = Self{
            .allocator = allocator, //
            .lower = 0,
            .upper_excl = number_of_nodes,
            .number_of_nodes = number_of_nodes,
            .working_seq = working_seq,
            .best_seq = best_seq,
            .num_calls = 0,
            .num_cache_hits = 0,
            .recursion_depth = 0,
            .time_last_report = now,
            .smallest_depth_last_report = 0,
            .infeasible_cache = infeasible_cache,
            .solution_cache = solution_cache,
            .candidate_lists = candidate_lists,
            .start_time = now,
            .timeout_ms = null,
            .lower_bound_mode = false,
            .rng = null,
            .silent = false,
            .single_pass = false,
        };

        return solver;
    }

    pub fn deinit(self: *Self) void {
        self.solution_cache.deinit();
        self.infeasible_cache.deinit();
        self.candidate_lists.deinit();
        self.best_seq.deinit();
        self.working_seq.deinit();
    }

    pub fn setLowerBound(self: *Self, lb: Node) void {
        self.lower = lb;
    }

    /// Sets exclusive upper bound, i.e. a solution of size `ub - 1` is acceptable.
    pub fn setUpperBound(self: *Self, ub: Node) void {
        assert(ub > 0);
        self.upper_excl = ub;
    }

    pub fn setTimeout(self: *Self, duration_ms: u64) void {
        self.timeout_ms = duration_ms;
    }

    pub fn enableLowerBoundMode(self: *Self, rng: *std.rand.DefaultPrng) void {
        self.lower_bound_mode = true;
        self.rng = rng;
    }

    pub fn numberOfRecursiveCalls(self: *Self) u64 {
        return self.num_calls;
    }

    noinline fn shrinkGraph(self: *Self, comptime Graph: type, graph: *const Graph, slack: Node) anyerror!Node {
        const Smaller = mg.MatrixGraph(Graph.NumNodes / 2);
        var n = graph.has_neighbors.cardinality();

        var newNodeIdOfOld = try self.allocator.alloc(u32, Graph.NumNodes);
        defer self.allocator.free(newNodeIdOfOld);

        var oldNodeIdOfNew = try self.allocator.alloc(u32, Smaller.NumNodes);
        defer self.allocator.free(oldNodeIdOfNew);

        // compute mapping between old and new node ids
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

        var candidates = std.ArrayList(ContraScore).init(self.allocator);
        defer candidates.deinit();

        var contracted = Smaller.BitSet.new();

        if (FeatureReuseCandidates and self.candidate_lists.items.len > 0) {
            var last = &self.candidate_lists.items[self.candidate_lists.items.len - 1];
            try candidates.ensureTotalCapacity(last.candidates.len);
            for (last.candidates) |c| {
                if (!graph.has_neighbors.isSet(c.rem) or !graph.has_neighbors.isSet(c.sur)) {
                    continue;
                }

                assert(oldNodeIdOfNew[newNodeIdOfOld[c.rem]] == c.rem);
                assert(oldNodeIdOfNew[newNodeIdOfOld[c.sur]] == c.sur);
                candidates.appendAssumeCapacity(ContraScore{ .rem = newNodeIdOfOld[c.rem], .sur = newNodeIdOfOld[c.sur], .score = c.score, .redDeg = c.redDeg, .processed = c.processed });
            }

            var invalid = @ptrCast(*Graph.BitSet, @alignCast(@alignOf(Graph.BitSet), last.invalid));

            var u: Node = 0;
            while (u < n) : (u += 1) {
                if (invalid.isSet(oldNodeIdOfNew[u])) {
                    _ = contracted.setBit(u);
                }
            }

            self.candidate_lists.appendAssumeCapacity(CandidateListRef{
                .candidates = candidates.items,
                .num_nodes = Smaller.NumNodes,
                .invalid = @ptrCast(*void, &contracted),
            });
        }

        defer if (FeatureReuseCandidates and self.candidate_lists.items.len > 1) {
            _ = self.candidate_lists.pop();
        };

        var smaller_storage = try self.allocator.alloc(Smaller, 1);
        defer self.allocator.free(smaller_storage);

        var smaller = try Smaller.new(self.allocator);
        defer smaller.deinit();

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
        var result = try self.solve(Smaller, &smaller, slack);
        assert(self.upper_excl < upper_before);

        for (self.best_seq.items[old_contractions..]) |*ctr| {
            ctr.rem = oldNodeIdOfNew[ctr.rem];
            ctr.sur = oldNodeIdOfNew[ctr.sur];
        }

        return result;
    }

    pub fn solve(self: *Self, comptime Graph: type, graph: *const Graph, slack: Node) anyerror!Node {
        var n: Node = graph.has_neighbors.cardinality();

        if (n > 0 and n + self.working_seq.items.len != self.number_of_nodes) {
            std.debug.print("n={d} |cs|={d} initial={d}\n", .{ n, self.working_seq.items.len, self.number_of_nodes });
            @panic("incomplete cs");
        }

        if (self.timeout_ms) |dur_ms| {
            if ((try std.time.Instant.now()).since(self.start_time) > dur_ms * 1000_000) {
                std.debug.print(" Timeout of {d} ms reached\n", .{dur_ms});
                return SolverError.Timeout;
            }
        }

        if (FeatureShrinkGraph and Graph.NumNodes > 8 and 2 * n < Graph.NumNodes) {
            return self.shrinkGraph(Graph, graph, slack);
        }

        self.num_calls += 1;
        self.recursion_depth += 1;
        defer self.smallest_depth_last_report = @min(self.smallest_depth_last_report, self.recursion_depth);
        defer self.recursion_depth -= 1;

        if (FeatureInitialPruning and self.num_calls == 1) {
            var local_graph = try graph.copy();
            defer local_graph.deinit();

            try self.initial_kernelization(Graph, &local_graph, slack);
            return self.solve(Graph, &local_graph, slack);
        }

        if (FeatureReportProgress > 0 and self.num_calls % FeatureReportProgress == 0) {
            var now = try std.time.Instant.now();
            var elapsed = now.since(self.time_last_report);
            self.time_last_report = now;

            std.debug.print("Call: {d} Depth: {d} (smallest: {d}) Graph({d}) n={d} slack={d} upper={d} time={d}ms\n", //
                .{ //
                self.num_calls,
                self.recursion_depth,
                self.smallest_depth_last_report,
                Graph.NumNodes,
                graph.has_neighbors.cardinality(),
                slack,
                self.upper_excl,
                elapsed / 1000_000,
            });
            self.smallest_depth_last_report = self.recursion_depth;
        }

        var digest: mg.Digest = undefined;
        if ((FeatureUseInfeasibleCache or FeatureUseSolutionCache) and !self.single_pass) {
            if (FeatureUseWFHash) {
                digest = graph.hashWF();
            } else {
                digest = graph.hashMatrix();
            }

            if (self.infeasible_cache.contains(digest)) {
                self.num_cache_hits += 1;
                return SolverError.Infeasable;
            }

            if (self.solution_cache.get(digest)) |sol| {
                if (sol.tww < self.upper_excl) {
                    self.num_cache_hits += 1;
                    return sol.tww;
                }

                if (sol.slack < sol.tww and sol.tww >= self.upper_excl) {
                    return SolverError.Infeasable;
                }
            }
        }

        var frame = try self.allocator.alloc(Frame(Graph), 1);
        defer self.allocator.free(frame);

        frame[0] = try Frame(Graph).new(self, graph, slack);
        defer frame[0].deinit();
        var result = frame[0].run();

        if (FeatureUseInfeasibleCache and result == SolverError.Infeasable) {
            try self.infeasible_cache.put(digest, {});
        }

        if (FeatureUseSolutionCache and digest != null) {
            if (result) |tww| {
                try self.solution_cache.put(digest, SolSummary{ .slack = slack, .tww = tww });
            } else |_| {}
        }

        return result;
    }

    pub fn initial_kernelization(self: *Self, comptime Graph: type, graph: *Graph, slack: Node) !void {
        _ = slack;

        if (!FeatureLeafNodePruning) {
            return;
        }

        assert(self.working_seq.items.len == 0);

        if (!self.silent) {
            std.debug.print("Before Kernel: n={d} m={d}\n", .{ graph.has_neighbors.cardinality(), graph.numberOfEdges() });
        }

        while (true) {
            _ = try self.pruneLeafs(Graph, graph);
            if (try self.pruneTwins(Graph, graph)) {
                continue;
            }

            break;
        }

        for (self.working_seq.items) |ctr| {
            assert(graph.deg(ctr.rem) == 0);
        }

        if (!self.silent) {
            std.debug.print("After Kernel: n={d} m={d}\n", .{ graph.has_neighbors.cardinality(), graph.numberOfEdges() });
        }
    }

    noinline fn pruneTwins(
        self: *Self,
        comptime Graph: type,
        graph: *Graph,
    ) anyerror!bool {
        if (!FeatureTwinPruning) {
            return false;
        }

        var u_iter = graph.has_neighbors.iter_set();
        while (u_iter.next()) |u| {
            var v_iter = graph.constTwoNeighbors(u).iter_set();

            while (v_iter.next()) |v| {
                if (v >= u) {
                    break;
                }

                var red = graph.redNeighborsAfterMerge(u, v);
                if (red.is_equal(graph.constRedNeighbors(u)) or
                    red.is_equal(graph.constRedNeighbors(v)))
                {
                    //std.debug.print("Initial prune twin {d} -> {d}\n", .{ u, v });
                    graph.mergeNodes(v, u, &red);
                    self.working_seq.appendAssumeCapacity(Contraction{ .rem = v, .sur = u });

                    assert(graph.deg(v) == 0);
                    return true;
                }
            }
        }

        return false;
    }

    noinline fn pruneLeafs(
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
                //std.debug.print("Initial leaf {d} -> {d}\n", .{ rem, survivor });
                assert(graph.deg(survivor) == 1);
                assert(graph.deg(rem) == 1);

                self.working_seq.appendAssumeCapacity(Contraction{ .rem = rem, .sur = survivor });
                graph.mergeNodes(rem, survivor, null); // TODO: this can be made faster by removing the Edge

                assert(graph.deg(survivor) == 1);

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
        mergables: Graph.BitSet,
        has_critical_neighbor: Graph.BitSet,

        work_graph: Graph,
        work_tww: Node,
        work_contracted: Graph.BitSet,
        slack: Node,
        lower_bound_mode: bool,

        fn new(solver: *ExactBranchAndBound, graph: *const Graph, slack: Node) !Self {
            var n = @max(4, @as(usize, graph.has_neighbors.cardinality()));

            var cand_size: usize = (n + 1) * n / 2;
            if (!FeaturePessimiticCandidateAllocation) {
                var n_prev = solver.candidate_lists.items.len;
                if (solver.candidate_lists.items.len > 0) {
                    cand_size = @min(cand_size, 2 * solver.candidate_lists.items[n_prev - 1].candidates.len);
                }
            }

            var candidates = try std.ArrayList(ContraScore).initCapacity(solver.allocator, cand_size);
            errdefer candidates.deinit();

            var work_graph = try Graph.new(solver.allocator);
            errdefer work_graph.deinit();

            return Self{
                .context = solver,

                .input_graph = graph,
                .mergables = undefined,
                .has_critical_neighbor = undefined,

                .candidates = candidates,

                .work_graph = work_graph,
                .work_contracted = Graph.BitSet.new(),
                .slack = slack,
                .work_tww = slack,

                .lower_bound_mode = solver.lower_bound_mode,
            };
        }

        fn deinit(self: *Self) void {
            self.work_graph.deinit();
            self.candidates.deinit();
        }

        fn run(self: *Self) !Node {
            if (false and self.input_graph.has_neighbors.cardinality() > 43) {
                for (self.context.working_seq.items) |x| {
                    _ = x;
                    std.debug.print(" ", .{});
                }
                std.debug.print("LB: {d} Slack: {d} UB: {d} n={d} d={d}\n", .{ self.context.lower, self.slack, self.context.upper_excl, self.input_graph.has_neighbors.cardinality(), self.context.candidate_lists.items.len });
            }

            // It's wasteful to reach this point with a slack that is too large. We should have pruned it earlier.
            assert(self.slack < self.context.upper_excl);

            // It's not optimal to invoke the solver on an empty graph.
            // However, we only reach this point if we found a strictly better solution, which happens rarely.
            // Doing it this way simplifies the code and avoids corner cases.
            if (self.input_graph.numberOfEdges() == 0) {
                self.registerSolution();
                return self.slack;
            }

            self.computeCandidates();
            const initial_seq_len = self.context.working_seq.items.len;
            var result: SolverError!Node = SolverError.Infeasable;

            for (self.candidates.items) |*c| {
                if (FeatureSkipProcessed and c.processed) {
                    continue;
                } else {
                    c.processed = true;
                }

                if (c.redDeg >= self.context.upper_excl) {
                    continue;
                }

                assert(self.input_graph.has_neighbors.isSet(c.rem) and self.input_graph.has_neighbors.isSet(c.sur));

                // undo all changes possibly made by an earlier iteration
                self.context.working_seq.shrinkRetainingCapacity(initial_seq_len);
                self.work_graph.copyFrom(self.input_graph);
                self.work_tww = self.slack;
                self.work_contracted.unsetAll();

                try self.contractNodes(c.rem, c.sur);

                if (!((self.lower_bound_mode and self.work_tww == 0) or (c.redDeg <= self.work_tww))) {
                    std.debug.print("lbm: {d} red: {d} tww: {d}", .{ @boolToInt(self.lower_bound_mode), c.redDeg, self.work_tww });
                }

                assert((self.lower_bound_mode and self.work_tww == 0) or (c.redDeg <= self.work_tww));

                if (self.work_tww >= self.context.upper_excl) {
                    continue;
                }

                if (FeatureReuseCandidates and !self.lower_bound_mode) {
                    self.context.candidate_lists.appendAssumeCapacity(CandidateListRef{
                        .num_nodes = Graph.NumNodes, //
                        .candidates = self.candidates.items[0..],
                        .invalid = @ptrCast(*void, &self.work_contracted),
                    });
                }
                defer if (FeatureReuseCandidates and !self.lower_bound_mode) {
                    _ = self.context.candidate_lists.pop();
                };

                if (FeatureComplements and (Graph.NumNodes > 8 and self.work_graph.numEdgesInComplement() * 5 < self.work_graph.numberOfEdges() * 4)) {
                    // std.debug.print("Complement!\n", .{});
                    self.work_graph.complement();
                }

                const tww_from_rec = self.context.solve(Graph, &self.work_graph, self.work_tww) catch |e| {
                    if (e == SolverError.Infeasable) {
                        if (self.context.lower >= self.context.upper_excl) {
                            // we may have found a better LB in the mean time
                            return SolverError.Infeasable;
                        } else {
                            continue;
                        }
                    } else {
                        return e;
                    }
                };
                result = tww_from_rec;

                assert(tww_from_rec >= self.work_tww);

                if (tww_from_rec <= @max(self.slack, self.context.lower)) {
                    break;
                }
            }

            if (result == SolverError.Infeasable) {}

            return result;
        }

        fn contractNodes(self: *Self, rem: Node, sur: Node) !void {
            if (rem < sur) {
                return self.contractNodes(sur, rem);
            }

            self.contractNodesWithoutPruning(rem, sur);

            if (self.work_tww >= self.context.upper_excl) {
                return;
            }

            while (true) {
                _ = try self.pruneLeafsAt(sur);

                if (try self.pruneTwinsAt(sur)) {
                    continue;
                }

                if (try self.pruneGeneralizedTwinsAt(sur)) {
                    continue;
                }

                if (try self.pruneOuterPath()) {
                    continue;
                }

                if (try self.tinyRedBridgesAt(sur)) {
                    continue;
                }

                if (try self.pruneBlackCCs()) {
                    continue;
                }

                self.tinyGraphBelowSlack();

                break;
            }
        }

        fn contractNodesWithoutPruning(self: *Self, rem: Node, sur: Node) void {
            if (rem < sur) {
                return self.contractNodesWithoutPruning(sur, rem);
            }

            self.work_graph.mergeNodes(rem, sur, null);

            self.work_tww = @max(self.work_tww, self.work_graph.redDeg(sur));
            self.work_tww = @max(self.work_tww, self.work_graph.maxRedDegreeIn(self.work_graph.constRedNeighbors(sur)));
            {
                var three = self.work_graph.closedNeighborsOfSet(self.work_graph.constTwoNeighbors(sur));
                self.work_contracted.assignOr(&three);
            }
            _ = self.work_contracted.setBit(sur);

            self.context.working_seq.appendAssumeCapacity(Contraction{ .rem = rem, .sur = sur });
        }

        fn pruneBlackCCs(self: *Self) anyerror!bool {
            if (!FeaturePruneBlackCCs) {
                return false;
            }

            var iter = self.work_graph.cc_iter();
            while (iter.next()) |bcc| {
                assert(bcc.cardinality() > 1);

                if (bcc.cardinality() >= self.context.upper_excl) {
                    continue;
                }

                var reds = Graph.BitSet.new();
                var max_red: Node = 0;
                {
                    var riter = bcc.iter_set();
                    while (riter.next()) |u| {
                        reds.assignOr(self.work_graph.constRedNeighbors(u));
                        max_red = @max(max_red, self.work_graph.redDeg(u));
                    }
                }
                reds.assignSub(&bcc);

                if (reds.cardinality() > max_red) {
                    continue;
                }

                var card = bcc.cardinality();

                if (card == 2) {
                    var it = bcc.iter_set();
                    var u = it.next().?;
                    var v = it.next().?;

                    try self.contractNodes(u, v);
                    return true;
                } else if (card == 3 and max_red < self.slack) {
                    var it = bcc.iter_set();
                    var u = it.next().?;
                    var v = it.next().?;
                    var w = it.next().?;

                    self.contractNodesWithoutPruning(w, v);
                    try self.contractNodes(v, u);
                    return true;
                } else {
                    std.debug.print("pruneBlackCC with {d} nodes with m {d} and s {d}\n", .{ bcc.cardinality(), max_red, self.slack });
                }
            }

            return false;
        }

        fn tinyRedBridgesAt(self: *Self, host: Node) anyerror!bool {
            _ = host;
            if (!FeaturePruneTinyRedBridges) {
                return false;
            }

            var u: u32 = 0;
            while (u < Graph.NumNodes) : (u += 1) {
                var deg = self.work_graph.deg(u);

                if (deg > 2 or deg != self.work_graph.redDeg(u)) {
                    continue;
                }

                var neigh = self.work_graph.constNeighbors(u).iter_set();
                while (neigh.next()) |v| {
                    var rdeg_v = self.work_graph.redDeg(v);
                    if (self.work_graph.deg(v) != rdeg_v) {
                        continue;
                    }

                    try self.contractNodes(u, v);
                    return true;
                }
            }

            return false;
        }

        fn tinyGraphBelowSlack(self: *Self) void {
            if (!FeatureTinyGraphBelowSlack or self.work_graph.has_neighbors.cardinality() > self.work_tww + 1) {
                return;
            }

            if (self.work_graph.has_neighbors.areAllUnset()) {
                return;
            }

            var iter = self.work_graph.has_neighbors.iter_set();
            var sur = iter.next().?;
            if (!self.context.silent) {
                std.debug.print("Collaps tiny graph with n={d} tww={d}\n", .{ self.work_graph.has_neighbors.cardinality(), self.work_tww });
            }

            while (iter.next()) |rem| {
                self.work_graph.removeEdgesAtNode(rem);
                self.context.working_seq.appendAssumeCapacity(Contraction{ .rem = rem, .sur = sur });
            }

            assert(self.work_graph.numberOfEdges() == 0);
        }

        fn computeCriticalNeighbors(self: *Self) void {
            if (!FeatureCriticalCandidates) {
                return;
            }

            self.has_critical_neighbor = Graph.BitSet.new();

            if (self.slack + 1 < self.context.upper_excl) {
                return;
            }

            var u: u32 = 0;
            while (u < Graph.NumNodes) : (u += 1) {
                if (self.input_graph.redDeg(u) != self.slack) {
                    continue;
                }

                self.has_critical_neighbor.assignOr(&self.input_graph.blackNeighborsCopied(u));
            }
        }

        fn computeCandidates(self: *Self) void {
            var has_neighbors: *const Graph.BitSet = undefined;
            if (FeatureSkipPathNodes) {
                self.computeMergeables();
                has_neighbors = &self.mergables;
                assert(self.mergables.is_subset_of(&self.input_graph.has_neighbors));
            } else {
                has_neighbors = &self.input_graph.has_neighbors;
            }

            self.computeCriticalNeighbors();

            const depth = self.context.candidate_lists.items.len;
            if (!FeatureReuseCandidates or self.lower_bound_mode or depth == 0) {
                // compute from scratch
                var outer_iter = has_neighbors.iter_set();
                while (outer_iter.next()) |u| {
                    var inner_iter = outer_iter;
                    while (inner_iter.next()) |v| {
                        self.appendCandidatePair(u, v);
                    }
                }
            } else {
                // reuse candidates from previous step
                var prev: *const CandidateListRef = &self.context.candidate_lists.items[depth - 1];
                assert(prev.num_nodes == Graph.NumNodes);

                var invalid: *const Graph.BitSet = @ptrCast(*Graph.BitSet, @alignCast(@alignOf(Graph.BitSet), prev.invalid));

                for (prev.candidates) |c| {
                    const rem = c.rem;
                    const sur = c.sur;

                    if (!has_neighbors.isSet(rem) or !has_neighbors.isSet(sur)) {
                        continue;
                    }

                    if (invalid.isSet(rem) or invalid.isSet(sur)) {
                        continue;
                    }

                    if (FeatureParanoidCandidates) {
                        // debugging only
                        var redDeg = self.input_graph.redDegreeInNeighborhoodAfterMerge(c.rem, c.sur);
                        if (c.redDeg > redDeg) {
                            std.debug.print("FAILED: old-red={d} new-red={d}", .{ c.redDeg, redDeg });
                        }

                        assert(c.redDeg == redDeg);
                    }

                    if (c.redDeg < self.context.upper_excl) {
                        self.candidates.appendAssumeCapacity(c);
                    }
                }

                var outer_iter = invalid.iter_set();
                while (outer_iter.next()) |u| {
                    if (!self.input_graph.has_neighbors.isSet(u)) {
                        continue;
                    }

                    var inner_iter = has_neighbors.iter_set();
                    while (inner_iter.next()) |v| {
                        if (v >= u) {
                            break;
                        }

                        self.appendCandidatePair(u, v);
                    }
                }
            }

            std.sort.sort(ContraScore, self.candidates.items, {}, cmpByValue);
        }

        fn computeMergeables(self: *Self) void {
            self.mergables = Graph.BitSet.new();
            var iter = self.input_graph.has_neighbors.iter_set();
            while (iter.next()) |u| {
                if (self.input_graph.deg(u) != 2 or !self.input_graph.constRedNeighbors(u).areAllUnset()) {
                    self.mergables.assignOr(self.input_graph.constNeighbors(u));
                    _ = self.mergables.setBit(u);
                }
            }

            if (self.mergables.areAllUnset()) {
                // happens for a graph consisting of disjoint cycles
                self.mergables = self.input_graph.has_neighbors;
            }
        }

        fn appendCandidatePair(self: *Self, u: Node, v: Node) void {
            assert(u != v);

            var distAtmostTwo = self.input_graph.constTwoNeighbors(u).isSet(v);
            var red_degree = Graph.NumNodes;

            if (distAtmostTwo or FeatureParanoidCandidates) {
                //if (FeatureCriticalCandidates and (self.has_critical_neighbor.isSet(u) != self.has_critical_neighbor.isSet(v))) {
                //    return;
                //}
                red_degree = self.input_graph.redDegreeInNeighborhoodAfterMerge(u, v);
            } else {
                if (FeatureCriticalCandidates and (self.has_critical_neighbor.isSet(u) or self.has_critical_neighbor.isSet(v))) {
                    assert(!FeatureParanoia or self.input_graph.redDegreeInNeighborhoodAfterMerge(u, v) >= self.context.upper_excl);
                    return;
                }
                red_degree = self.input_graph.deg(u) + self.input_graph.deg(v);

                assert(!FeatureParanoia or red_degree <= self.input_graph.redDegreeInNeighborhoodAfterMerge(u, v));
            }

            if (FeatureSkipInfeasibleCandidatesEarly and red_degree >= self.context.upper_excl) {
                return;
            }

            var score = red_degree;
            var red_before = @min(red_degree, @max(self.input_graph.redDeg(u), self.input_graph.redDeg(v)));

            score -= red_before / 2;
            score += @boolToInt(distAtmostTwo);

            var candidate =
                ContraScore{
                .score = score, //
                .redDeg = red_degree,
                .rem = @max(u, v),
                .sur = @min(u, v),
                .processed = false,
            };

            if (FeaturePessimiticCandidateAllocation) {
                self.candidates.appendAssumeCapacity(candidate);
            } else {
                self.candidates.append(candidate) catch @panic("OOM");
            }
        }

        fn registerSolution(self: *Self) void {
            assert(self.context.upper_excl > self.work_tww);
            if (!self.context.silent) {
                std.debug.print("Found a solution with tww={d}\n", .{self.work_tww});
            }

            self.context.upper_excl = self.work_tww;
            self.context.best_seq.clearRetainingCapacity();
            self.context.best_seq.appendSliceAssumeCapacity(self.context.working_seq.items);
        }

        noinline fn pruneOuterPath(self: *Self) anyerror!bool {
            if (!FeatureOuterPathPruning or self.slack < 2) {
                return false;
            }

            var u: Node = 0;
            while (u < Graph.NumNodes) : (u += 1) {
                if (self.work_graph.deg(u) != 1) {
                    continue;
                }

                var v: Node = undefined;
                {
                    var iter = self.work_graph.constNeighbors(u).iter_set();
                    v = iter.next().?;
                }

                if (self.work_graph.deg(v) != 2) {
                    continue;
                }

                var w: Node = undefined;
                {
                    var iter = self.work_graph.constNeighbors(v).iter_set();
                    w = iter.next().?;
                    if (w == u) {
                        w = iter.next().?;
                    }
                }

                if (!self.work_graph.constRedNeighbors(v).isSet(w) and self.work_graph.deg(w) != 2) {
                    continue;
                }

                try self.contractNodes(u, v);
                return true;
            }

            return false;
        }

        noinline fn pruneLeafsAt(self: *Self, u: Node) anyerror!bool {
            if (!FeatureLeafNodePruning) {
                return false;
            }

            var neighs = self.work_graph.constNeighbors(u).iter_set();

            var l1 = while (neighs.next()) |v| {
                if (self.work_graph.deg(v) == 1) {
                    break v;
                }
            } else {
                return false;
            };

            while (neighs.next()) |v| {
                if (self.work_graph.deg(v) == 1) {
                    try self.contractNodes(l1, v);
                    return true;
                }
            }

            return false;
        }

        noinline fn pruneTwinsAt(self: *Self, u: Node) anyerror!bool {
            if (!FeatureTwinPruning or FeaturePruneGeneralizedTwins) {
                return false;
            }

            var two_neighbors = self.work_graph.constTwoNeighbors(u);
            if (two_neighbors.cardinality() > 50) {
                return false;
            }

            var it1 = two_neighbors.iter_set();
            while (it1.next()) |v| {
                var it2 = it1;
                while (it2.next()) |w| {
                    var red = self.work_graph.redNeighborsAfterMerge(v, w);
                    if (red.is_equal(self.work_graph.constRedNeighbors(v)) or
                        red.is_equal(self.work_graph.constRedNeighbors(w)))
                    {
                        try self.contractNodes(v, w);
                        return true;
                    }
                }
            }

            return false;
        }

        noinline fn pruneGeneralizedTwinsAt(self: *Self, u: Node) anyerror!bool {
            if (!FeaturePruneGeneralizedTwins) {
                return false;
            }

            var two_neighbors = self.work_graph.constTwoNeighbors(u);
            if (two_neighbors.cardinality() > 100) {
                return false;
            }

            var it1 = two_neighbors.iter_set();
            while (it1.next()) |v| {
                var it2 = it1;
                var red_deg_v = self.work_graph.redDeg(v);
                while (it2.next()) |w| {
                    var red_deg_w = self.work_graph.redDeg(w);
                    var red = self.work_graph.redNeighborsAfterMerge(v, w);

                    var red_card = red.cardinality();
                    if (red_card > red_deg_v or red_card > red_deg_w) {
                        continue;
                    }

                    red.assignSub(self.work_graph.constRedNeighbors(v));
                    red.assignSub(self.work_graph.constRedNeighbors(w));

                    if (!red.areAllUnset()) {
                        continue;
                    }

                    try self.contractNodes(v, w);
                    return true;
                }
            }

            return false;
        }

        fn cmpByValue(_: void, a: ContraScore, b: ContraScore) bool {
            return (a.score < b.score) or (a.score == b.score and (a.rem < b.rem or (a.rem == b.rem and a.sur < b.sur)));
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

pub fn findLowerBoundSubgraph(comptime N: u32, input_graph: *const mg.MatrixGraph(N), allocator: std.mem.Allocator, budget_ms: u64, lower: Node, upper: Node) !Node {
    const Graph = mg.MatrixGraph(N);
    var lower_bound = lower;

    var graph = try input_graph.copy();
    defer graph.deinit();
    {
        var solver = try ExactBranchAndBound.new(allocator, graph.numberOfNodes());
        defer solver.deinit();
        try solver.initial_kernelization(Graph, &graph, lower);
    }

    var n = graph.has_neighbors.cardinality();
    var max_size = @max(5 * upper, 18);
    if (n < max_size) {
        // if the graph is too small, we solve it directly
        // if its way too large, no chance we can solve it within the budget
        return lower;
    }
    var start_time: std.time.Instant = try std.time.Instant.now();
    var rng = std.rand.DefaultPrng.init(0x1273_5412_327a_8431);
    var random = rng.random();

    var nodes = Graph.BitSet.new();

    while (true) {
        if ((try std.time.Instant.now()).since(start_time) > budget_ms * 1000_000) {
            return lower_bound;
        }

        if (lower_bound >= upper) {
            return lower_bound;
        }

        nodes.unsetAll();
        while (true) {
            var node = random.intRangeLessThan(u32, 0, graph.numberOfNodes());
            if (graph.has_neighbors.isSet(node)) {
                _ = nodes.setBit(node);
                break;
            }
        }

        while (nodes.cardinality() < max_size) {
            var frontier = Graph.BitSet.new();
            var iter = nodes.iter_set();
            while (iter.next()) |u| {
                frontier.assignOr(graph.constNeighbors(u));
            }
            frontier.assignSub(&nodes);

            var nNodes = nodes.cardinality();
            while (frontier.cardinality() + nNodes > max_size) {
                var node = random.intRangeLessThan(u32, 0, graph.numberOfNodes());
                _ = frontier.unsetBit(node);
            }

            nodes.assignOr(&frontier);
        }

        lower_bound = subgraphLB(N, allocator, &graph, &nodes, budget_ms / 2, lower_bound, upper) catch continue;
    }

    return lower_bound;
}

fn subgraphLB(comptime N: u32, allocator: std.mem.Allocator, graph: *const mg.MatrixGraph(N), nodes: *const mg.MatrixGraph(N).BitSet, budget_ms: u64, lower: Node, upper: Node) !Node {
    const Graph = mg.MatrixGraph(N);

    var local_graph = try graph.copy();
    defer local_graph.deinit();
    {
        var it = nodes.iter_unset();
        while (it.next()) |i| {
            local_graph.removeEdgesAtNode(i);
        }
    }

    //    std.debug.print("Attempt LB with subgraph n={d} m={d} lower={d} upper={d}\n", .{ nodes.cardinality(), local_graph.numberOfEdges(), lower, upper });

    var solver = ExactBranchAndBound.new(allocator, local_graph.has_neighbors.cardinality()) catch {
        return lower;
    };
    defer solver.deinit();

    solver.setLowerBound(lower);
    solver.setUpperBound(upper);
    solver.setTimeout(budget_ms);
    solver.silent = true;

    var tww = solver.solve(Graph, &local_graph, 0) catch |e| blk: {
        if (e == SolverError.Infeasable) {
            std.debug.print(" Found optimal LB via sugraph!\n", .{});
            return upper;
        }
        break :blk lower;
    };

    if (tww > lower) {
        std.debug.print(" Found LB via subgraph with n={d}: {d}!\n", .{ nodes.cardinality(), tww });
    }

    return @max(lower, tww);
}

fn solveCCExactlyHandler(context: anytype, graph: anytype, mapping: anytype) void {
    const T = @TypeOf(context.*).IntType;
    var cs: *contr.ContractionSequence(T) = &context.cc.best_contraction_sequence;

    std.debug.print("Invoke exact solver on n={d}, m={d}, lower={d}, upper={d}, maxRed={d}\n", .{ graph.has_neighbors.cardinality(), graph.numberOfEdges(), context.lower, context.upper, graph.maxRedDegree() });

    var lower = context.lower;
    if (FeatureTryLBFirst) {
        var boost: u64 = if (graph.has_neighbors.cardinality() > 128) 6 else 1;

        if (BudgetSubgraphLB > 0) {
            lower = findLowerBoundSubgraph(@TypeOf(graph).NumNodes, &graph, context.allocator, BudgetSubgraphLB * boost, context.lower, context.upper) catch context.lower;
            assert(lower >= context.lower);
            if (lower >= context.upper) {
                context.result = SolverError.Infeasable;
                return;
            }
        }
    }

    var solver = ExactBranchAndBound.new(context.allocator, graph.has_neighbors.cardinality()) catch |e| {
        context.result = e;
        return;
    };
    defer solver.deinit();

    solver.setLowerBound(lower);
    solver.setUpperBound(context.upper);

    var tww = solver.solve(@TypeOf(graph), &graph, 0) catch |e| {
        std.debug.print("Iterations: {d} Cache-Hits: {d} ({d:.1} %)", .{
            solver.num_calls, //
            solver.num_cache_hits,
            100.0 * @intToFloat(f64, solver.num_cache_hits) / @intToFloat(f64, solver.num_calls),
        });
        std.debug.print("Catch {!}\n", .{e});
        context.result = e;
        return;
    };

    std.debug.print("Iterations: {d} Cache-Hits: {d} ({d:.1} %)\n", .{
        solver.num_calls, //
        solver.num_cache_hits,
        100.0 * @intToFloat(f64, solver.num_cache_hits) / @intToFloat(f64, solver.num_calls),
    });

    context.result = tww;

    // update CC
    context.cc.tww = @intCast(T, tww);
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
