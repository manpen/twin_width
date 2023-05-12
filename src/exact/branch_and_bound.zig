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

pub const ExactBranchAndBound = struct {
    const Self = @This();

    allocator: std.mem.Allocator,
    lower: Node,
    upper_excl: Node,
    number_of_nodes: Node,

    working_seq: std.ArrayList(Contraction),
    best_seq: std.ArrayList(Contraction),

    num_calls: u64,
    recursion_depth: usize,

    pub fn new(allocator: std.mem.Allocator, number_of_nodes: Node) !Self {
        var working_seq = try std.ArrayList(Contraction).initCapacity(allocator, number_of_nodes);
        errdefer working_seq.deinit();

        var best_seq = try std.ArrayList(Contraction).initCapacity(allocator, number_of_nodes);
        errdefer best_seq.deinit();

        var solver = Self{ .allocator = allocator, .lower = 0, .upper_excl = number_of_nodes, .number_of_nodes = number_of_nodes, .working_seq = working_seq, .best_seq = best_seq, .num_calls = 0, .recursion_depth = 0 };

        return solver;
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

    pub fn solve(self: *Self, comptime Graph: type, graph: *const Graph, slack: Node) anyerror!Node {
        self.num_calls += 1;
        self.recursion_depth += 1;
        std.debug.print("Depth: {d}\n", .{self.recursion_depth});
        var frame = Frame(Graph).new(self, graph.*, slack);
        var result = frame.run();
        self.recursion_depth -= 1;
        return result;
    }

    pub fn deinit(self: *Self) void {
        self.working_seq.deinit();
        self.best_seq.deinit();
    }
};

fn Frame(comptime Graph: type) type {
    return struct {
        const Self = @This();
        const ContraScore = struct { score: Node, rem: Node, sur: Node };

        context: *ExactBranchAndBound,
        graph: Graph,
        slack: Node,

        fn new(solver: *ExactBranchAndBound, graph: Graph, slack: Node) Self {
            return Self{
                .context = solver,
                .graph = graph,
                .slack = slack,
            };
        }

        fn run(self: *Self) !Node {
            //for (self.context.working_seq.items) |x| {_=x;std.debug.print(" ", .{});}
            //std.debug.print("LB: {d} Slack: {d} UB: {d}\n", .{self.context.lower, self.slack, self.context.upper_excl});

            if (self.slack >= self.context.upper_excl) {
                return SolverError.Infeasable;
            }

            if (self.graph.numberOfEdges() == 0) {
                std.debug.print("Found a solution with tww={d}\n", .{self.slack});

                self.context.upper_excl = self.slack;
                self.context.best_seq.clearRetainingCapacity();
                try self.context.best_seq.appendSlice(self.context.working_seq.items);

                return self.slack;
            }

            const initial_seq_len = self.context.working_seq.items.len;

            var candidates = try self.compute_candidates();
            defer candidates.deinit();

            var result: SolverError!Node = SolverError.Infeasable;

            for (candidates.items) |c| {
                try self.context.working_seq.resize(initial_seq_len);

                var local_graph = self.graph;
                local_graph.mergeNodes(c.rem, c.sur, null);
                var tww = @max(self.slack, local_graph.maxRedDegree());

                if (tww >= self.context.upper_excl) {
                    continue;
                }

                try self.context.working_seq.append(Contraction{ .rem = c.rem, .sur = c.sur });

                _ = try self.eliminate_leafs(&local_graph);

                const tww_from_recursion = self.context.solve(@TypeOf(local_graph), &local_graph, tww) catch continue;
                std.debug.assert(tww_from_recursion >= tww);

                result = tww_from_recursion;

                if (tww_from_recursion <= @max(self.slack, self.context.lower)) {
                    break;
                }
            }

            return result;
        }

        fn compute_candidates(self: *Self) !std.ArrayList(ContraScore) {
            var result = try std.ArrayList(ContraScore).initCapacity(self.context.allocator, self.graph.numberOfNodes());
            errdefer result.deinit();

            var nodes = self.graph.has_neighbors;
            var u: Node = 0;
            while (u < self.graph.numberOfNodes()) : (u += 1) {
                if (!nodes.unsetBit(u)) {
                    std.debug.assert(self.graph.deg(u) == 0);
                    continue;
                }
                std.debug.assert(self.graph.deg(u) > 0);

                var iter = nodes.iter_set();
                while (iter.next()) |v| {
                    assert(self.graph.deg(u) > 0);
                    assert(self.graph.deg(v) > 0);

                    var reds = self.graph.redNeighborsAfterMerge(u, v);
                    var red_degree = reds.cardinality();

                    if (red_degree >= self.context.upper_excl) {
                        continue;
                    }

                    try result.append(ContraScore{ .score = red_degree, .rem = v, .sur = u });
                }
            }

            std.sort.sort(ContraScore, result.items, {}, cmpByValue);

            return result;
        }

        fn cmpByValue(context: void, a: ContraScore, b: ContraScore) bool {
            _ = context;
            return a.score < b.score or (a.score == b.score and (a.rem < b.rem or (a.rem == b.rem and a.sur < b.sur)));
        }

        fn eliminate_leafs(self: *Self, graph: *Graph) !bool {
            var leafs = Graph.BitSet.new();
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
                if (self.graph.deg(u) < 2) {
                    continue;
                }
                var leafs_at_node = leafs.copyWithAnd(self.graph.constNeighbors(u));

                if (leafs_at_node.cardinality() < 2) {
                    continue;
                }

                var leafs_iter = leafs_at_node.iter_set();
                const survivor = leafs_iter.next().?;
                while (leafs_iter.next()) |rem| {
                    try self.context.working_seq.append(Contraction{ .rem = rem, .sur = survivor });
                    graph.mergeNodes(rem, survivor, null); // TODO: this can be made faster by removing the Edge
                    change = true;
                }
            }

            return change;
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

pub fn solveCCExactly(comptime T: type, component: *cc.ConnectedComponent(T), allocator: std.mem.Allocator, lower: Node, upper: Node) !Node {
    const ContextType = solveCCContext(T);
    var context = ContextType{ .cc = component, .allocator = allocator, .lower = lower, .upper = upper, .result = upper + 1 };
    try bs.matrixGraphFromInducedSubGraph(T, &component.subgraph, &context, solveCCExactlyHandler);
    return context.result;
}

fn solveCCExactlyHandler(context: anytype, graph: anytype, mapping: anytype) void {
    const T = @TypeOf(context.*).IntType;

    std.debug.print("Invoke exact solver on n={d}, m={d}, lower={d}, upper={d}", .{ graph.has_neighbors.cardinality(), graph.numberOfEdges(), context.lower, context.upper });

    var solver = ExactBranchAndBound.new(context.allocator, graph.numberOfNodes()) catch |e| {
        context.result = e;
        return;
    };
    defer solver.deinit();

    solver.setLowerBound(context.lower);
    solver.setUpperBound(context.upper);

    var tww = solver.solve(@TypeOf(graph), &graph, 0) catch |e| {
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
        cs.addContraction(OuterContraction{ .erased = @intCast(T, mapping[ctr.rem]), .survivor = @intCast(T, mapping[ctr.sur]) }) catch unreachable;
    }
}

test "Init BB solver" {
    var allocator = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = allocator.deinit();
    var solver = try ExactBranchAndBound.new(allocator.allocator(), 100);
    defer solver.deinit();
}
