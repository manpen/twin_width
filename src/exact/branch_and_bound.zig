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

const ContraScore = struct {
    score: Node,
    redDeg: Node,
    rem: Node,
    sur: Node,
    processed: bool,
    distTwo: bool,
};

const CandidateListRef = struct {
    num_nodes: u32,
    candidates: []ContraScore,
    invalid: *void,
};

const FeatureReuseCandidates: bool = true;
const FeatureUseInfeasibleCache: bool = true;
const FeatureShrinkGraph: bool = true;
const FeatureSkipProcessed: bool = true;
const FeatureSkipInfeasibleCandidatesEarly: bool = true;
const FeatureSkipPathNodes: bool = true;
const FeatureComplements: bool = false; // BROKEN

const FeatureReportProgress: u32 = 1 * 100000; // set to 0 to disable

const FeatureInitialPruning: bool = true;
const FeatureLeafNodePruning: bool = true;
const FeatureTwinPruning: bool = true;
const FeatureTinyGraphBelowSlack: bool = true; // almost no effect, as it only helps, if we're about to find a solution

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

    start_time: std.time.Instant,
    timeout_ms: ?u64,

    lower_bound_mode: bool,
    rng: ?*std.rand.DefaultPrng,

    pub fn new(allocator: std.mem.Allocator, number_of_nodes: Node) !Self {
        var working_seq = try std.ArrayList(Contraction).initCapacity(allocator, number_of_nodes);
        errdefer working_seq.deinit();

        var best_seq = try std.ArrayList(Contraction).initCapacity(allocator, number_of_nodes);
        errdefer best_seq.deinit();

        var candidate_lists = try std.ArrayList(CandidateListRef).initCapacity(allocator, @boolToInt(FeatureReuseCandidates) * (20 + number_of_nodes));
        errdefer candidate_lists.deinit();

        var cache = InfeasibleCache.init(allocator);
        errdefer cache.deinit();
        if (FeatureUseInfeasibleCache) {
            try cache.ensureTotalCapacity((2 << 30) / @sizeOf(mg.Digest));
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
            .recursion_depth = 0,
            .time_last_report = now,
            .smallest_depth_last_report = 0,
            .infeasible_cache = cache,
            .candidate_lists = candidate_lists,
            .start_time = now,
            .timeout_ms = null,
            .lower_bound_mode = false,
            .rng = null,
        };

        return solver;
    }

    pub fn deinit(self: *Self) void {
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
        std.debug.assert(ub > 0);
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

                std.debug.assert(oldNodeIdOfNew[newNodeIdOfOld[c.rem]] == c.rem);
                std.debug.assert(oldNodeIdOfNew[newNodeIdOfOld[c.sur]] == c.sur);
                candidates.appendAssumeCapacity(ContraScore{ .rem = newNodeIdOfOld[c.rem], .sur = newNodeIdOfOld[c.sur], .score = c.score, .redDeg = c.redDeg, .distTwo = c.distTwo, .processed = c.processed });
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
        std.debug.assert(self.upper_excl < upper_before);

        for (self.best_seq.items[old_contractions..]) |*ctr| {
            ctr.rem = oldNodeIdOfNew[ctr.rem];
            ctr.sur = oldNodeIdOfNew[ctr.sur];
        }

        return result;
    }

    pub fn solve(self: *Self, comptime Graph: type, graph: *const Graph, slack: Node) anyerror!Node {
        var n: Node = graph.has_neighbors.cardinality();

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

        if (self.lower_bound_mode and n < self.upper_excl * 4) {
            std.debug.assert(slack == 0);
            std.debug.print("Disable lower bound mode\n", .{});
            self.lower_bound_mode = false;
        }

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
        if (FeatureUseInfeasibleCache) {
            digest = graph.hash();
            if (self.infeasible_cache.contains(digest)) {
                return SolverError.Infeasable;
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

        return result;
    }

    fn initial_kernelization(self: *Self, comptime Graph: type, graph: *Graph, slack: Node) !void {
        _ = slack;
        if (!FeatureLeafNodePruning) {
            return;
        }

        while (true) {
            _ = try self.pruneLeafs(Graph, graph);
            if (try self.pruneTwins(Graph, graph)) {
                continue;
            }

            break;
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
                if (v <= u) {
                    continue;
                }

                var red = graph.redNeighborsAfterMerge(u, v);
                if (red.is_equal(graph.constRedNeighbors(u)) or
                    red.is_equal(graph.constRedNeighbors(v)))
                {
                    std.debug.print("Initial prune twin\n", .{});
                    graph.mergeNodes(u, v, &red);
                    self.working_seq.appendAssumeCapacity(Contraction{ .rem = u, .sur = v });
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
                std.debug.print("Initial leaf {d} -> {d}\n", .{ rem, survivor });
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
        mergables: Graph.BitSet,

        work_graph: Graph,
        work_tww: Node,
        work_contracted: Graph.BitSet,
        slack: Node,
        lower_bound_mode: bool,

        fn new(solver: *ExactBranchAndBound, graph: *const Graph, slack: Node) !Self {
            var n = @max(4, @as(usize, graph.has_neighbors.cardinality()));

            var candidates = try std.ArrayList(ContraScore).initCapacity(solver.allocator, (n + 1) * n / 2);
            errdefer candidates.deinit();

            var work_graph = try Graph.new(solver.allocator);
            errdefer work_graph.deinit();

            return Self{
                .context = solver,

                .input_graph = graph,
                .mergables = undefined,

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
            std.debug.assert(self.slack < self.context.upper_excl);

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

                std.debug.assert(self.input_graph.has_neighbors.isSet(c.rem) and self.input_graph.has_neighbors.isSet(c.sur));

                // undo all changes possibly made by an earlier iteration
                try self.context.working_seq.resize(initial_seq_len);
                self.work_graph.copyFrom(self.input_graph);
                self.work_tww = self.slack;
                self.work_contracted.unsetAll();

                try self.contractNodes(c.rem, c.sur, !c.distTwo);

                if (!((self.lower_bound_mode and self.work_tww == 0) or (c.redDeg <= self.work_tww))) {
                    std.debug.print("lbm: {d} red: {d} tww: {d} dist2: {any}", .{ @boolToInt(self.lower_bound_mode), c.redDeg, self.work_tww, c.distTwo });
                }

                std.debug.assert((self.lower_bound_mode and self.work_tww == 0) or (c.redDeg <= self.work_tww));

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
                    std.debug.print("Complement!", .{});
                    self.work_graph.complement();
                }

                const tww_from_rec = self.context.solve(Graph, &self.work_graph, self.work_tww) catch |e| {
                    if (e == SolverError.Infeasable) {
                        continue;
                    } else {
                        return e;
                    }
                };
                result = tww_from_rec;

                std.debug.assert(tww_from_rec >= self.work_tww);

                if (tww_from_rec <= @max(self.slack, self.context.lower)) {
                    break;
                }
            }

            return result;
        }

        fn contractNodes(self: *Self, rem: Node, sur: Node, atleast_dist_three: bool) !void {
            std.debug.assert(self.work_graph.deg(rem) > 0);
            std.debug.assert(self.work_graph.deg(sur) > 0);

            if (self.lower_bound_mode) {
                self.work_contracted.assignOr(self.work_graph.constTwoNeighbors(rem));

                self.work_graph.removeEdgesAtNode(rem);
            } else {
                self.work_graph.mergeNodes(rem, sur, null);
                self.work_tww = @max(self.work_tww, self.work_graph.redDeg(sur));
                self.work_tww = @max(self.work_tww, self.work_graph.maxRedDegreeIn(self.work_graph.constRedNeighbors(sur)));
                self.work_contracted.assignOr(self.work_graph.constTwoNeighbors(sur));
                _ = self.work_contracted.setBit(sur);
            }

            if (self.work_tww >= self.context.upper_excl) {
                return;
            }

            self.context.working_seq.appendAssumeCapacity(Contraction{ .rem = rem, .sur = sur });

            while (true) {
                if (!atleast_dist_three) {
                    _ = try self.pruneLeafsAt(sur);
                }

                if (try self.pruneTwinsAt(sur)) {
                    continue;
                }

                self.tinyGraphBelowSlack();

                break;
            }
        }

        fn tinyGraphBelowSlack(self: *Self) void {
            if (!FeatureTinyGraphBelowSlack or self.work_graph.has_neighbors.cardinality() > self.work_tww) {
                return;
            }

            if (self.work_graph.has_neighbors.areAllUnset()) {
                return;
            }

            var iter = self.work_graph.has_neighbors.iter_set();
            var sur = iter.next().?;
            std.debug.print("Collaps tiny graph with n={d} tww={d}\n", .{ self.work_graph.has_neighbors.cardinality(), self.work_tww });

            while (iter.next()) |rem| {
                self.work_graph.removeEdgesAtNode(rem);
                self.context.working_seq.appendAssumeCapacity(Contraction{ .rem = rem, .sur = sur });
            }
        }

        fn computeCandidates(self: *Self) void {
            var has_neighbors: *const Graph.BitSet = undefined;
            if (FeatureSkipPathNodes) {
                self.computeMergeables();
                has_neighbors = &self.mergables;
                std.debug.assert(self.mergables.is_subset_of(&self.input_graph.has_neighbors));
            } else {
                has_neighbors = &self.input_graph.has_neighbors;
            }

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
                std.debug.assert(prev.num_nodes == Graph.NumNodes);

                var invalid: *const Graph.BitSet = @ptrCast(*Graph.BitSet, @alignCast(@alignOf(Graph.BitSet), prev.invalid));

                for (prev.candidates) |c| {
                    const rem = c.rem;
                    const sur = c.sur;

                    if (!has_neighbors.isSet(rem) or !has_neighbors.isSet(sur)) {
                        continue;
                    }

                    if (!invalid.isSet(rem) and !invalid.isSet(sur)) {
                        if (c.redDeg < self.context.upper_excl) {
                            self.candidates.appendAssumeCapacity(c);
                        }
                        continue;
                    }
                }

                var outer_iter = invalid.iter_set();
                while (outer_iter.next()) |u| {
                    if (!self.input_graph.has_neighbors.isSet(u)) {
                        continue;
                    }

                    var inner_iter = has_neighbors.iter_set();
                    while (inner_iter.next()) |v| {
                        if (u < v) {
                            self.appendCandidatePair(u, v);
                        }
                    }
                }
            }

            std.sort.sort(ContraScore, self.candidates.items, {}, cmpByValue);

            if (self.lower_bound_mode) {
                var rand = self.context.rng.?.random();
                var i: usize = 1;

                rand.shuffle(ContraScore, self.candidates.items[0..@min(self.candidates.items.len, 10)]);

                while (i < @min(3, self.candidates.items.len)) : (i += 1) {
                    if (rand.int(u8) > 16) {
                        break;
                    }
                }

                if (i < self.candidates.items.len) {
                    self.candidates.shrinkRetainingCapacity(i);
                }
            }
        }

        fn computeMergeables(self: *Self) void {
            self.mergables = Graph.BitSet.new();
            var iter = self.input_graph.has_neighbors.iter_set();
            while (iter.next()) |u| {
                if (self.input_graph.deg(u) != 2 or !self.input_graph.constRedNeighbors(u).areAllUnset()) {
                    self.mergables.assignOr(self.input_graph.constTwoNeighbors(u));
                }
            }

            if (self.mergables.areAllUnset()) {
                // happens for a graph consisting of disjoint cycles
                self.mergables = self.input_graph.has_neighbors;
            }
        }

        fn appendCandidatePair(self: *Self, u: Node, v: Node) void {
            std.debug.assert(u != v);

            var distAtmostTwo = self.input_graph.constTwoNeighbors(u).isSet(v);
            var red_degree = Graph.NumNodes;

            if (distAtmostTwo) {
                red_degree = self.input_graph.redDegreeInNeighborhoodAfterMerge(u, v);
            } else {
                red_degree = self.input_graph.deg(u) + self.input_graph.deg(v);
            }

            if (FeatureSkipInfeasibleCandidatesEarly and red_degree >= self.context.upper_excl) {
                return;
            }

            self.candidates.appendAssumeCapacity(ContraScore{
                .score = red_degree, //
                .redDeg = red_degree,
                .rem = v,
                .sur = u,
                .distTwo = distAtmostTwo,
                .processed = false,
            });
        }

        fn registerSolution(self: *Self) void {
            std.debug.assert(self.context.upper_excl > self.work_tww);
            std.debug.print("Found a solution with tww={d}\n", .{self.work_tww});

            self.context.upper_excl = self.work_tww;
            self.context.best_seq.clearRetainingCapacity();
            self.context.best_seq.appendSliceAssumeCapacity(self.context.working_seq.items);
        }

        noinline fn pruneLeafsAt(self: *Self, u: Node) !bool {
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

        noinline fn pruneTwinsAt(self: *Self, u: Node) anyerror!bool {
            if (!FeatureTwinPruning) {
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
                        //std.debug.print("Twin be gone", .{});
                        try self.contractNodes(v, w, false);
                        return true;
                    }
                }
            }

            return false;
        }

        fn cmpByValue(_: void, a: ContraScore, b: ContraScore) bool {
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

pub fn findLowerBound(comptime N: u32, graph: *const mg.MatrixGraph(N), allocator: std.mem.Allocator, budget_ms: u64, lower: Node, upper: Node) !Node {
    var lower_bound = lower;

    var n = graph.has_neighbors.cardinality();
    if (n < @min(20, 4 * upper)) {
        // if the graph is too small, we solve it directly
        // if its way too large, no chance we can solve it within the budget
        return lower;
    }

    var start_time: std.time.Instant = try std.time.Instant.now();
    var rng = std.rand.DefaultPrng.init(0x1273_5412_327a_8431);

    for (0..100) |_| {
        if ((try std.time.Instant.now()).since(start_time) > budget_ms * 1000_000) {
            return lower_bound;
        }

        if (lower_bound >= upper) {
            return lower_bound;
        }

        var solver = ExactBranchAndBound.new(allocator, graph.numberOfNodes()) catch continue;
        defer solver.deinit();

        solver.setLowerBound(lower_bound);
        solver.setUpperBound(upper);
        solver.setTimeout(@max(budget_ms / 3, 1000));
        solver.enableLowerBoundMode(&rng);

        var this_lower = solver.solve(mg.MatrixGraph(N), graph, 0) catch |e| {
            if (e == SolverError.Infeasable) {
                std.debug.print(" Found optimal LB!\n", .{});
                return upper;
            }
            continue;
        };
        var now = try std.time.Instant.now();
        std.debug.print(" Established LB: {d:>8} (best: {d:>8}) gap: {d:>8} in {d:>6}ms\n", .{ this_lower, lower_bound, upper - this_lower, now.since(start_time) / 1000_000 });
        lower_bound = @max(this_lower, lower_bound);
    }

    return lower_bound;
}

pub fn findLowerBoundSubgraph(comptime N: u32, graph: *const mg.MatrixGraph(N), allocator: std.mem.Allocator, budget_ms: u64, lower: Node, upper: Node) !Node {
    const Graph = mg.MatrixGraph(N);
    var lower_bound = lower;

    var n = graph.has_neighbors.cardinality();
    if (n < @min(20, 4 * upper)) {
        // if the graph is too small, we solve it directly
        // if its way too large, no chance we can solve it within the budget
        return lower;
    }

    var start_time: std.time.Instant = try std.time.Instant.now();
    var rng = std.rand.DefaultPrng.init(0x1273_5412_327a_8431);
    var random = rng.random();

    var nodes = Graph.BitSet.new();
    var max_size = @max(4 * upper, n / 100);

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

        while (nodes.cardinality() * 3 < max_size * 2) {
            var prev = nodes;
            var iter = prev.iter_set();
            while (iter.next()) |u| {
                nodes.assignOr(graph.constNeighbors(u));
            }
        }

        while (nodes.cardinality() * 2 > max_size * 3) {
            var node = random.intRangeLessThan(u32, 0, graph.numberOfNodes());
            _ = nodes.unsetBit(node);
        }

        lower_bound = subgraphLB(N, allocator, graph, &nodes, budget_ms / 10, lower_bound, upper) catch continue;
    }

    return lower_bound;
}

pub fn findLowerBoundWithCS(comptime N: u32, graph: *const mg.MatrixGraph(N), cs: []Contraction, allocator: std.mem.Allocator, budget_ms: u64, lower: Node, upper: Node) !Node {
    const Graph = mg.MatrixGraph(N);

    var local_graph = try graph.copy();
    defer local_graph.deinit();
    var nodes_used = Graph.BitSet.new();

    {
        var tww: Node = 0;
        for (cs) |c| {
            local_graph.mergeNodes(c.rem, c.sur, null);
            _ = nodes_used.setBit(c.rem);
            _ = nodes_used.setBit(c.sur);
            nodes_used.assignOr(local_graph.constRedNeighbors(c.sur));

            tww = @max(tww, local_graph.maxRedDegree());
            if (tww == upper) {
                break;
            }
        } else {
            std.debug.print("Supplied contraction sequence is better than upper bound. Expected: {d} Found: {d}\n", .{ upper, tww });
            @panic("");
        }
    }

    return subgraphLB(N, allocator, graph, &nodes_used, budget_ms, lower, upper);
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

    std.debug.print("Attempt LB with subgraph n={d} m={d} lower={d} upper={d}\n", .{ nodes.cardinality(), local_graph.numberOfEdges(), lower, upper });

    var solver = ExactBranchAndBound.new(allocator, local_graph.numberOfNodes()) catch {
        return lower;
    };
    defer solver.deinit();

    solver.setLowerBound(lower);
    solver.setUpperBound(upper);
    solver.setTimeout(budget_ms);

    var tww = solver.solve(Graph, &local_graph, 0) catch |e| blk: {
        if (e == SolverError.Infeasable) {
            std.debug.print(" Found optimal LB via CS!\n", .{});
            return upper;
        }
        break :blk lower;
    };

    std.debug.print(" Found LB via subgraph with n={d}: {d}!\n", .{ nodes.cardinality(), tww });

    return @max(lower, tww);
}

fn solveCCExactlyHandler(context: anytype, graph: anytype, mapping: anytype) void {
    const T = @TypeOf(context.*).IntType;
    var cs: *contr.ContractionSequence(T) = &context.cc.best_contraction_sequence;

    std.debug.print("Invoke exact solver on n={d}, m={d}, lower={d}, upper={d}, maxRed={d}\n", .{ graph.has_neighbors.cardinality(), graph.numberOfEdges(), context.lower, context.upper, graph.maxRedDegree() });

    //var large_buffer = context.allocator.alloc(u8, 1024 * 1024 * 3000) catch @panic("OOM");
    //defer context.allocator.free(large_buffer);

    //var hpa_allocator = std.heap.FixedBufferAllocator.init(large_buffer);
    //var hpa = hpa_allocator.allocator();

    var lower = context.lower;
    if (true) {
        if (false) {
            var inner_cs = context.allocator.alloc(Contraction, cs.list.len) catch return;
            defer context.allocator.free(inner_cs);

            for (cs.list, inner_cs) |o, *i| {
                i.* = Contraction{ .rem = o.erased, .sur = o.survivor };
            }

            lower = findLowerBoundWithCS(@TypeOf(graph).NumNodes, &graph, inner_cs, context.allocator, 10_000, context.lower, context.upper) catch context.lower;
        }
        if (lower >= context.upper) {
            context.result = SolverError.Infeasable;
            return;
        }

        lower = findLowerBoundSubgraph(@TypeOf(graph).NumNodes, &graph, context.allocator, 10_000, context.lower, context.upper) catch context.lower;
        std.debug.assert(lower >= context.lower);
        if (lower >= context.upper) {
            context.result = SolverError.Infeasable;
            return;
        }

        lower = findLowerBound(@TypeOf(graph).NumNodes, &graph, context.allocator, 10_000, context.lower, context.upper) catch context.lower;
        std.debug.assert(lower >= context.lower);
        if (lower >= context.upper) {
            context.result = SolverError.Infeasable;
            return;
        }
    }

    var solver = ExactBranchAndBound.new(context.allocator, graph.numberOfNodes()) catch |e| {
        context.result = e;
        return;
    };
    defer solver.deinit();

    solver.setLowerBound(lower);
    solver.setUpperBound(context.upper);

    var tww = solver.solve(@TypeOf(graph), &graph, 0) catch |e| {
        std.debug.print("Catch {!}\n", .{e});
        context.result = e;
        return;
    };

    std.debug.print("Nodes: {d} CS: {d}", .{ graph.has_neighbors.cardinality(), solver.best_seq.items.len });

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
