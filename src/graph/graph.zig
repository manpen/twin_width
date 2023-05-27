const std = @import("std");
const bitset = @import("../util/two_level_bitset.zig");
const edge_list = @import("edge_list.zig");
const red_edge_stack = @import("red_edge_stack.zig");
const contraction = @import("../tww/contraction_sequence.zig");
const comptime_util = @import("../util/comptime_checks.zig");
const RetraceableContractionSequence = @import("../tww/retraceable_contraction_sequence.zig").RetraceableContractionSequence;
const connected_components = @import("connected_component.zig");
const Node = @import("node.zig").Node;
const bfs_mod = @import("bfs.zig");
const compressed_bitset = @import("../util/compressed_bitmap.zig");
const solver_resources = @import("solver.zig");
const min_hash_mod = @import("../util/min_hash.zig");

const solver_bnb = @import("../exact/branch_and_bound.zig");
const exact_bs = @import("../exact/bootstrapping.zig");

const pace_2023 = @import("../pace_2023/pace_fmt.zig");

pub const GraphError = error{ FileNotFound, //
NotPACEFormat, GraphTooLarge, MisformedEdgeList, InvalidContractionOneNodeErased, ContractionOverflow, NoContractionLeft, NegativeNumberOfLeafes, ExactSuboptimal };

pub fn Graph(comptime T: type) type {
    comptime if (!comptime_util.checkIfIsCompatibleInteger(T)) {
        @compileError("T must either be u8,u16 or u32!");
    };

    const promote_thresh = comptime if (T == u8) 0 else if (T == u16) 200 else 300;
    const degrade_tresh = comptime if (T == u8) 0 else if (T == u16) 100 else 200;

    return struct {
        const Self = @This();
        pub const IntType = T;
        pub const NodeType = Node(T, promote_thresh, degrade_tresh);
        pub const promote_thresh_inner = promote_thresh;
        number_of_nodes: u32,
        number_of_edges: u32,
        // Stores the all nodes the index is also the id of the node
        node_list: []NodeType,

        // Keep track of erased nodes since the node list will not shrink!
        erased_nodes: bitset.FastBitSet,

        // Scratch bitset for dfs/bfs visited etc.
        scratch_bitset: bitset.FastBitSet,

        // The current best contraction sequence
        contraction: contraction.ContractionSequence(T),

        // An allocator which always report OutOfMemory used to verify assumptions about
        // memory consumption being O(2*m) etc.
        failing_allocator: std.heap.FixedBufferAllocator,

        // Stores all connected components that are not trivial (< 2 vertices)
        connected_components: std.ArrayListUnmanaged(connected_components.ConnectedComponent(T)),

        // Stores a trivial merge sequence that merges all graphs of order < 2
        trivial_connected_component_contraction_sequence: std.ArrayListUnmanaged(u32),

        // Contains all nodes but permuted so that each slice can will store the nodes in
        // a connected component consecutive in memory
        connected_components_node_list_slice: []T,

        // Max heap keeping track of the component with the largest twin width so that the solver spends the majority of the time there
        connected_components_min_heap: std.PriorityQueue(connected_components.ConnectedComponentIndex(T), void, connected_components.ConnectedComponentIndex(T).compareComponentIndexDesc),

        // Main allocator is used to allocate everything that is needed at runtime.
        allocator: std.mem.Allocator,

        started_at: std.time.Instant,

        last_merge_first_level_merge: bool,
        last_merge_red_edges_erased: std.ArrayListUnmanaged(T),

        min_hash: min_hash_mod.MinHashSimilarity(T, 7),

        // Usually the exact tries to find a solution that is strictly better than
        // the heuristic one (oftentimes only proving a lower bound). If set to true,
        // we force the solver to produce a solution (by increasing the upper bound)
        // and report an error if it fails to do so. Useful for testing.
        force_exact_solver_to_solve: bool,

        pub const LargeListStorageType = compressed_bitset.FastCompressedBitmap(T, promote_thresh, degrade_tresh);

        pub inline fn addEdge(self: *Self, u: T, v: T) !void {
            @setEvalBranchQuota(3000);
            // At the moment clear the connected component since we have no efficient way to find out
            // to which connected component a node belongs
            self.connected_components.clearRetainingCapacity();
            // Clear out the min heap
            self.connected_components_min_heap.shrinkAndFree(0);
            std.debug.assert(u < self.node_list.len);
            std.debug.assert(v < self.node_list.len);

            const u_node = &self.node_list[u];
            const v_node = &self.node_list[v];

            // Normally addContraction keeps track of the number of leafes in each node the initialization
            // needs to keep track of that by itself
            if (u_node.cardinality() == 0) {
                v_node.num_leafes += 1;
            } else if (u_node.isLeaf()) {
                self.node_list[u_node.getFirstNeighboor()].num_leafes -= 1;
            }

            if (v_node.cardinality() == 0) {
                u_node.num_leafes += 1;
            } else if (v_node.isLeaf()) {
                self.node_list[v_node.getFirstNeighboor()].num_leafes -= 1;
            }

            _ = try u_node.addBlackEdge(self.allocator, v);
            _ = try v_node.addBlackEdge(self.allocator, u);

            self.number_of_edges += 1;
        }

        pub inline fn getCurrentTwinWidth(self: *Self) u32 {
            // Tww is the largest CC Tww
            if (self.connected_components_min_heap.peek()) |largest_cc| {
                return largest_cc.tww;
            }
            return 0;
        }

        pub inline fn checkUpdateNewLeaf(self: *Self, new_leaf: T, comptime increase: bool) void {
            if (self.node_list[new_leaf].isLeaf()) {
                const parent = self.node_list[new_leaf].getFirstNeighboor();
                if (increase) {
                    self.node_list[parent].num_leafes += 1;
                } else {
                    self.node_list[parent].num_leafes -= 1;
                }
            }
        }

        pub fn revertLastContraction(self: *Self, seq: *RetraceableContractionSequence(T)) !void {
            // WARNING: Is not usable yet since it will not restore num leafes at the moment this function will not restore exactly the same state as before!
            if (seq.lastContraction()) |last| {
                self.scratch_bitset.unsetAll();
                _ = self.erased_nodes.unset(last.erased);

                var black_iter = self.node_list[last.erased].black_edges.iterator();
                if (self.node_list[last.survivor].isLeaf()) {
                    const parent = self.node_list[last.survivor].getFirstNeighboor();
                    self.node_list[parent].num_leafes -= 1;
                }

                // Add black edges
                while (black_iter.next()) |black_edge| {
                    try self.node_list[black_edge].addBlackEdge(self.allocator, last.erased);
                }

                var red_iter = self.node_list[last.erased].red_edges.iterator();
                // Add black edges
                while (red_iter.next()) |red_edge| {
                    self.scratch_bitset.set(red_edge);
                    // Reduce survivor num leafes
                    try self.node_list[red_edge].addRedEdge(self.allocator, last.erased);
                    // noop

                    // noop
                    try self.node_list[red_edge].removeRedEdge(self.allocator, last.survivor);
                    // increase erased num leafes
                }
                self.checkUpdateNewLeaf(last.erased, true);

                // Remove all red edges which are set in the erased set
                self.node_list[last.survivor].red_edges.removeMask(&self.scratch_bitset);

                // Add back all red edges that were deleted by the merge and all black edges that were deleted by the merge
                var iter_last = try seq.red_edge_stack.iterateLastLevel();
                while (iter_last.next()) |red_edge| {
                    switch (red_edge.edge_type) {
                        .black_to_red_own => {
                            // Remove red edge
                            _ = try self.node_list[last.survivor].red_edges.remove(self.allocator, red_edge.target);
                            _ = try self.node_list[red_edge.target].red_edges.remove(self.allocator, last.survivor);

                            // Readd black edge
                            try self.node_list[last.survivor].black_edges.add(self.allocator, red_edge.target);
                            try self.node_list[red_edge.target].black_edges.add(self.allocator, last.survivor);
                        },
                        .black_to_red_other => {
                            // Remove red edge
                            _ = try self.node_list[last.survivor].red_edges.remove(self.allocator, red_edge.target);

                            _ = try self.node_list[red_edge.target].red_edges.remove(self.allocator, last.survivor);
                        },
                        .black_to_deleted => {
                            try self.node_list[last.survivor].black_edges.add(self.allocator, red_edge.target);

                            try self.node_list[red_edge.target].black_edges.add(self.allocator, last.survivor);
                        },
                        .red_to_deleted => {
                            try self.node_list[last.survivor].red_edges.add(self.allocator, red_edge.target);

                            try self.node_list[red_edge.target].red_edges.add(self.allocator, last.survivor);
                        },
                    }
                }
                self.updateLeafCount(last.survivor);

                // Inform the red edge stack about the revert
                try seq.removeLast();
            } else {
                @panic("No last contraction to revert!");
            }
        }

        pub const InducedTwinWidthPotential = struct {
            tww: T,
            cumulative_red_edges: i64,
            delta_red_edges: i32,
            pub inline fn default() InducedTwinWidthPotential {
                return InducedTwinWidthPotential{
                    .tww = std.math.maxInt(T),
                    .cumulative_red_edges = std.math.maxInt(i64),
                    .delta_red_edges = std.math.maxInt(i32),
                };
            }

            pub inline fn reset(self: *InducedTwinWidthPotential) void {
                self.tww = std.math.maxInt(T);
                self.cumulative_red_edges = std.math.maxInt(i64);
                self.delta_red_edges = std.math.maxInt(i32);
            }

            pub fn compare(ctx: void, self: InducedTwinWidthPotential, other: InducedTwinWidthPotential) std.math.Order {
                _ = ctx;
                if (self.tww == other.tww) {
                    return std.math.order(self.cumulative_red_edges, other.cumulative_red_edges);
                }
                return std.math.order(self.tww, other.tww);
            }

            pub inline fn order(self: InducedTwinWidthPotential, other: InducedTwinWidthPotential, current_twin_width: T) std.math.Order {
                if (self.tww >= current_twin_width or other.tww >= current_twin_width) {
                    if (self.tww < other.tww) {
                        return .lt;
                    } else if (self.tww == other.tww) {
                        return std.math.order(self.cumulative_red_edges, other.cumulative_red_edges);
                    } else {
                        return .gt;
                    }
                } else {
                    if (self.cumulative_red_edges < other.cumulative_red_edges) {
                        return .lt;
                    } else if (self.cumulative_red_edges == other.cumulative_red_edges) {
                        return std.math.order(other.delta_red_edges, self.delta_red_edges);
                    } else {
                        return .gt;
                    }
                }
            }

            pub inline fn isLessOrEqual(self: InducedTwinWidthPotential, other: InducedTwinWidthPotential, current_twin_width: T) bool {
                if (self.tww >= current_twin_width or other.tww >= current_twin_width) {
                    if (self.tww <= other.tww) {
                        return true;
                    } else {
                        return false;
                    }
                } else {
                    if (self.cumulative_red_edges <= other.cumulative_red_edges) {
                        return true;
                    }
                }
                return false;
            }

            pub inline fn isLess(self: InducedTwinWidthPotential, other: InducedTwinWidthPotential, current_twin_width: T) bool {
                if (self.tww >= current_twin_width or other.tww >= current_twin_width) {
                    if (self.tww < other.tww or (self.tww == other.tww and self.cumulative_red_edges < other.cumulative_red_edges)) {
                        return true;
                    } else {
                        return false;
                    }
                } else {
                    if (self.cumulative_red_edges < other.cumulative_red_edges) {
                        return true;
                    } else if (self.cumulative_red_edges == other.cumulative_red_edges) {
                        return self.delta_red_edges > other.delta_red_edges;
                    } else {
                        return false;
                    }
                }
            }

            pub inline fn isLessDeltaRedEdgesMajor(self: InducedTwinWidthPotential, other: InducedTwinWidthPotential, current_twin_width: T) bool {
                if (self.tww >= current_twin_width or other.tww >= current_twin_width) {
                    if (self.tww < other.tww or (self.tww == other.tww and self.delta_red_edges < other.delta_red_edges)) {
                        return true;
                    } else {
                        return false;
                    }
                } else {
                    if (self.delta_red_edges < other.delta_red_edges) {
                        return true;
                    } else if (self.delta_red_edges == other.delta_red_edges) {
                        return self.cumulative_red_edges < other.cumulative_red_edges;
                    } else {
                        return false;
                    }
                }
            }
        };

        pub const TwwScorer = struct {
            tww: T,
            tww_nb: T,
            delta_red_edges: i32,

            pub fn default() TwwScorer {
                return .{ .tww = std.math.maxInt(T), .tww_nb = std.math.maxInt(T), .delta_red_edges = std.math.maxInt(i32) };
            }

            pub fn betterTww(self: *TwwScorer, other: *TwwScorer, current_tww: T) bool {
                if (self.tww_nb <= current_tww) {
                    return self.tww < other.tww;
                } else {
                    return self.tww_nb <= other.tww_nb and self.tww < other.tww;
                }
            }

            pub fn better(self: *const TwwScorer, other: *const TwwScorer, current_tww: T) bool {
                if (self.tww_nb <= current_tww) {
                    if (self.delta_red_edges == other.delta_red_edges) return self.tww < other.tww;
                    return self.delta_red_edges < other.delta_red_edges;
                } else {
                    return self.tww_nb < other.tww_nb or (self.tww_nb == other.tww_nb and self.tww < other.tww);
                }
            }
        };

        pub fn calculateMaxTwwScore(self: *Self, erased: T, survivor: T) TwwScorer {
            var delta_red: T = 0;
            var tww: T = 0;
            var red_iter = self.node_list[erased].red_edges.iterator();
            var delta_red_edges: i32 = 0;
            var correction_factor: T = 0;

            while (red_iter.next()) |item| {
                if (item == survivor) {
                    delta_red_edges -= 1;
                    correction_factor = 1;
                    continue;
                }

                if (!self.node_list[survivor].red_edges.contains(item)) {
                    delta_red += 1;
                } else {
                    delta_red_edges -= 1;
                }
            }

            delta_red += self.node_list[survivor].red_edges.cardinality();
            delta_red -= correction_factor;

            var black_iter = self.node_list[erased].black_edges.xorIterator(&self.node_list[survivor].black_edges);

            while (black_iter.next()) |item| {
                if (item == survivor or item == erased) {
                    continue;
                }
                // Came from erased
                if (black_iter.first) {
                    if (!self.node_list[survivor].red_edges.contains(item)) {
                        delta_red += 1;
                        delta_red_edges += 1;
                        tww = std.math.max(self.node_list[item].red_edges.cardinality() + 1, tww);
                    }
                }
                // Came from survivor
                else {
                    if (!self.node_list[erased].red_edges.contains(item)) {
                        delta_red += 1;
                        delta_red_edges += 1;
                        tww = std.math.max(self.node_list[item].red_edges.cardinality() + 1, tww);
                    }
                }
            }
            return .{ .tww = delta_red, .tww_nb = std.math.max(delta_red, tww), .delta_red_edges = delta_red_edges };
        }

        pub fn calculateMaxTwwOfNewNeighbors(self: *Self, erased: T, survivor: T) struct { T, T } {
            var delta_red: T = 0;
            var tww: T = 0;
            var red_iter = self.node_list[erased].red_edges.iterator();

            while (red_iter.next()) |item| {
                if (item == survivor) continue;

                if (!self.node_list[survivor].red_edges.contains(item)) {
                    delta_red += 1;
                    tww = std.math.max(self.node_list[item].red_edges.cardinality() + 1, tww);
                }
            }

            delta_red += self.node_list[survivor].red_edges.cardinality();

            var black_iter = self.node_list[erased].black_edges.xorIterator(&self.node_list[survivor].black_edges);

            while (black_iter.next()) |item| {
                if (item == survivor or item == erased) {
                    continue;
                }
                // Came from erased
                if (black_iter.first) {
                    if (!self.node_list[survivor].red_edges.contains(item)) {
                        delta_red += 1;
                        tww = std.math.max(self.node_list[item].red_edges.cardinality() + 1, tww);
                    }
                }
                // Came from survivor
                else {
                    if (!self.node_list[erased].red_edges.contains(item)) {
                        delta_red += 1;
                        tww = std.math.max(self.node_list[item].red_edges.cardinality() + 1, tww);
                    }
                }
            }
            return .{ std.math.max(delta_red, tww), delta_red };
        }

        pub fn calculateTwwOfMergeSurvivor(self: *Self, erased: T, survivor: T) T {
            var delta_red: T = 0;
            var red_iter = self.node_list[erased].red_edges.iterator();

            while (red_iter.next()) |item| {
                if (item == survivor) continue;

                if (!self.node_list[survivor].red_edges.contains(item)) {
                    delta_red += 1;
                }
            }

            delta_red += self.node_list[survivor].red_edges.cardinality();

            var black_iter = self.node_list[erased].black_edges.xorIterator(&self.node_list[survivor].black_edges);

            while (black_iter.next()) |item| {
                if (item == survivor or item == erased) {
                    continue;
                }
                // Came from erased
                if (black_iter.first) {
                    if (!self.node_list[survivor].red_edges.contains(item)) {
                        delta_red += 1;
                    }
                }
                // Came from survivor
                else {
                    if (!self.node_list[erased].red_edges.contains(item)) {
                        delta_red += 1;
                    }
                }
            }
            return delta_red;
        }

        pub fn calculateInducedTwwPotential(self: *Self, erased: T, survivor: T, ub: *InducedTwinWidthPotential, current_tww: T) InducedTwinWidthPotential {
            // NOTICE: This function performs better than calculateInducedTww at the moment
            var red_potential: i64 = 0;

            var delta_red: T = 0;

            var new_red_edges: i32 = 0;
            var red_iter = self.node_list[erased].red_edges.iterator();
            var before_red_card = self.node_list[survivor].red_edges.cardinality();

            // Intuition is as following:
            // New red edges to nodes with a lot of red edges should decrease
            // the probability of taking  this move

            while (red_iter.next()) |item| {
                if (item == survivor) {
                    continue;
                }

                if (!self.node_list[survivor].red_edges.contains(item)) {
                    delta_red += 1;
                } else {
                    // We destroyed a red edge update potential
                    red_potential -= (@intCast(i64, self.node_list[item].red_edges.cardinality()) * @intCast(i64, self.node_list[item].red_edges.cardinality()));
                    new_red_edges -= 1;
                }
            }

            delta_red += self.node_list[survivor].red_edges.cardinality();

            var tww: T = 0;
            var black_iter = self.node_list[erased].black_edges.xorIterator(&self.node_list[survivor].black_edges);

            while (black_iter.next()) |item| {
                if (item == survivor or item == erased) {
                    continue;
                }
                // Came from erased
                if (black_iter.first) {
                    if (!self.node_list[survivor].red_edges.contains(item)) {
                        delta_red += 1;
                        new_red_edges += 1;
                        red_potential += (@intCast(i64, (self.node_list[item].red_edges.cardinality() + 1)) * @intCast(i64, (self.node_list[item].red_edges.cardinality() + 1)));
                        tww = std.math.max(self.node_list[item].red_edges.cardinality() + 1, tww);
                    }
                }
                // Came from survivor
                else {
                    if (!self.node_list[erased].red_edges.contains(item)) {
                        delta_red += 1;
                        new_red_edges += 1;
                        red_potential += (@intCast(i64, (self.node_list[item].red_edges.cardinality() + 1)) * @intCast(i64, (self.node_list[item].red_edges.cardinality() + 1)));
                        tww = std.math.max(self.node_list[item].red_edges.cardinality() + 1, tww);
                    }
                }

                const current = InducedTwinWidthPotential{ .tww = tww, .cumulative_red_edges = red_potential, .delta_red_edges = std.math.maxInt(i32) };
                if (!current.isLess(ub.*, current_tww)) {
                    return current;
                }
            }

            if (before_red_card < delta_red) {
                // replace by small gaussian
                for (before_red_card + 1..delta_red + 1) |c| {
                    red_potential += (@intCast(i64, c * c));
                }
            } else {
                for (delta_red + 1..before_red_card + 1) |c| {
                    red_potential -= (@intCast(i64, c * c));
                }
            }

            tww = std.math.max(tww, delta_red);
            return InducedTwinWidthPotential{ .tww = tww, .cumulative_red_edges = red_potential, .delta_red_edges = new_red_edges };
        }

        pub const InducedTwinWidth = struct {
            tww: T,
            delta_red_edges: i32,
            pub inline fn default() InducedTwinWidth {
                return InducedTwinWidth{
                    .tww = std.math.maxInt(T),
                    .delta_red_edges = std.math.maxInt(i32),
                };
            }

            pub inline fn reset(self: *InducedTwinWidth) void {
                self.tww = std.math.maxInt(T);
                self.delta_red_edges = std.math.maxInt(i32);
            }

            pub inline fn isLess(self: InducedTwinWidth, other: InducedTwinWidth) bool {
                return self.lessThan(other);
            }
            pub inline fn lessThan(self: InducedTwinWidth, other: InducedTwinWidth) bool {
                if (self.tww < other.tww or (self.tww == other.tww and self.delta_red_edges < other.delta_red_edges)) {
                    return true;
                } else {
                    return false;
                }
            }
        };

        pub fn calculateInducedTww(self: *Self, erased: T, survivor: T, upper_bound: ?T) InducedTwinWidth {
            const erased_cardinality_red = self.node_list[erased].red_edges.cardinality();
            const survivor_cardinality_red = self.node_list[survivor].red_edges.cardinality();
            if (upper_bound) |ub| {
                const erased_cardinality_black = self.node_list[erased].black_edges.cardinality();
                const survivor_cardinality_black = self.node_list[survivor].black_edges.cardinality();

                const erased_cardinality = erased_cardinality_black + erased_cardinality_red;
                const survivor_cardinality = survivor_cardinality_black + survivor_cardinality_red;

                //heuristic_024.gr
                const min_dist_black = std.math.max(erased_cardinality_black, survivor_cardinality_black) - std.math.min(erased_cardinality_black, survivor_cardinality_black);

                if (min_dist_black > ub) {
                    return InducedTwinWidth{
                        .tww = min_dist_black,
                        .delta_red_edges = std.math.maxInt(i32),
                    };
                }

                const min_dist_total = std.math.max(erased_cardinality, survivor_cardinality) - std.math.min(erased_cardinality, survivor_cardinality);

                // Fast exit
                if (min_dist_total > ub) {
                    return InducedTwinWidth{ .tww = min_dist_total, .delta_red_edges = std.math.maxInt(i32) };
                }
            }

            const upper_bound_all = upper_bound orelse std.math.maxInt(T);

            var tww: T = 0;

            var delta_red: T = 0;
            var correction_factor: T = 0;
            if (erased_cardinality_red == 0 or survivor_cardinality_red == 0) {
                delta_red = std.math.max(erased_cardinality_red, survivor_cardinality_red);
            } else {
                var red_iter = self.node_list[erased].red_edges.xorIterator(&self.node_list[survivor].red_edges);

                while (red_iter.next()) |item| {
                    if (item != survivor and item != erased) {
                        delta_red += 1;
                        tww = std.math.max(self.node_list[item].red_edges.cardinality(), tww);
                    } else {
                        correction_factor = 1;
                    }
                }
            }

            var delta_red_edges = @intCast(i32, (delta_red + correction_factor)) - @intCast(i32, erased_cardinality_red + survivor_cardinality_red);

            if (delta_red > upper_bound_all) {
                return InducedTwinWidth{ .tww = delta_red, .delta_red_edges = delta_red_edges };
            }

            var black_iter = self.node_list[erased].black_edges.xorIterator(&self.node_list[survivor].black_edges);

            while (black_iter.next()) |item| {
                if (item == survivor or item == erased) {
                    continue;
                }
                // Came from erased
                if (black_iter.first) {
                    if (!self.node_list[survivor].red_edges.contains(item)) {
                        delta_red += 1;
                        delta_red_edges += 1;
                        tww = std.math.max(self.node_list[item].red_edges.cardinality() + 1, tww);
                    }
                }
                // Came from survivor
                else {
                    if (!self.node_list[erased].red_edges.contains(item)) {
                        delta_red += 1;
                        delta_red_edges += 1;
                        tww = std.math.max(self.node_list[item].red_edges.cardinality() + 1, tww);
                    }
                }

                if (tww > upper_bound_all or delta_red > upper_bound_all) {
                    return InducedTwinWidth{ .tww = std.math.max(tww, delta_red), .delta_red_edges = delta_red_edges };
                }
            }

            tww = std.math.max(tww, delta_red);
            return InducedTwinWidth{ .tww = tww, .delta_red_edges = delta_red_edges };
        }

        pub inline fn updateLeafCount(self: *Self, node: T) void {
            if (self.node_list[node].isLeaf()) {
                const parent = self.node_list[node].getFirstNeighboor();
                self.node_list[parent].num_leafes += 1;
                self.node_list[node].num_leafes = if (self.node_list[parent].isLeaf()) 1 else 0;
            } else {
                var nb_iter = self.node_list[node].unorderedIterator();
                var count: T = 0;
                while (nb_iter.next()) |item| {
                    if (self.node_list[item].isLeaf()) {
                        count += 1;
                    }
                }
                self.node_list[node].num_leafes = count;
            }
        }

        pub fn addContractionNoMinHash(self: *Self, erased: T, survivor: T, seq: *RetraceableContractionSequence(T)) !T {
            if (self.erased_nodes.get(erased) or self.erased_nodes.get(survivor)) {
                std.debug.print("Result {} {}\n", .{ erased, survivor });
                return GraphError.InvalidContractionOneNodeErased;
            } else if (erased == survivor) {
                return GraphError.MisformedEdgeList;
            }

            self.last_merge_first_level_merge = false;
            self.last_merge_red_edges_erased.shrinkRetainingCapacity(0);

            self.erased_nodes.set(erased);
            if (self.node_list[survivor].isLeaf()) {
                const parent = self.node_list[survivor].getFirstNeighboor();
                self.node_list[parent].num_leafes -= 1;
            }
            if (self.node_list[erased].isLeaf()) {
                const parent = self.node_list[erased].getFirstNeighboor();
                self.node_list[parent].num_leafes -= 1;
            }

            var red_iter = self.node_list[erased].red_edges.iterator();
            while (red_iter.next()) |item| {
                try self.node_list[item].removeRedEdge(self.allocator, erased);
                if (item != survivor) {
                    if (!try self.node_list[survivor].addRedEdgeExists(self.allocator, item)) {
                        try self.node_list[item].addRedEdge(self.allocator, survivor);
                        self.last_merge_red_edges_erased.append(self.allocator, item) catch unreachable;
                    } else {
                        // Inform about the removal of the red edge
                        try seq.red_edge_stack.addEdge(self.failing_allocator.allocator(), red_edge_stack.NewRedEdge(T).redToDeleted(item));
                    }
                } else {
                    self.last_merge_first_level_merge = true;
                }
            }

            var tww: T = 0;
            var black_iter = self.node_list[erased].black_edges.xorIterator(&self.node_list[survivor].black_edges);

            var remove_list = std.ArrayList(T).init(self.allocator);
            defer remove_list.deinit();

            while (black_iter.next()) |item| {
                if (item == survivor or item == erased) {
                    self.last_merge_first_level_merge = true;
                    continue;
                }
                // Came from erased
                if (black_iter.first) {
                    if (!try self.node_list[survivor].addRedEdgeExists(self.allocator, item)) {
                        try seq.red_edge_stack.addEdge(self.failing_allocator.allocator(), red_edge_stack.NewRedEdge(T).blackToRedOther(item));

                        try self.node_list[item].addRedEdge(self.allocator, survivor);

                        self.last_merge_red_edges_erased.append(self.allocator, item) catch unreachable;
                    }
                }
                // Came from survivor
                else {
                    if (try self.node_list[survivor].addRedEdgeExists(self.allocator, item)) {
                        // If it existed before we inherited from erased

                        try seq.red_edge_stack.addEdge(self.failing_allocator.allocator(), red_edge_stack.NewRedEdge(T).blackToDeleted(item));
                    } else {
                        // Did not exist therefore turned
                        try seq.red_edge_stack.addEdge(self.failing_allocator.allocator(), red_edge_stack.NewRedEdge(T).blackToRedOwn(item));
                    }
                    try self.node_list[item].addRedEdge(self.allocator, survivor);

                    try self.node_list[item].removeBlackEdge(self.allocator, survivor);
                    try remove_list.append(item);
                }
                tww = std.math.max(tww, @intCast(T, self.node_list[item].red_edges.cardinality()));
            }

            for (remove_list.items) |item| {
                // Batch remove?
                try self.node_list[survivor].removeBlackEdge(self.allocator, item);
            }

            var black_remove_iter = self.node_list[erased].black_edges.iterator();
            while (black_remove_iter.next()) |item| {
                try self.node_list[item].removeBlackEdge(self.allocator, erased);
            }

            tww = std.math.max(tww, @intCast(T, self.node_list[survivor].red_edges.cardinality()));
            try seq.addContraction(self.allocator, erased, survivor, std.math.max(tww, seq.getTwinWidth()));

            self.updateLeafCount(survivor);
            return tww;
        }

        pub fn addContraction(self: *Self, erased: T, survivor: T, seq: *RetraceableContractionSequence(T)) !T {
            if (self.erased_nodes.get(erased) or self.erased_nodes.get(survivor)) {
                std.debug.print("Result {} {}\n", .{ erased, survivor });
                return GraphError.InvalidContractionOneNodeErased;
            } else if (erased == survivor) {
                return GraphError.MisformedEdgeList;
            }

            self.last_merge_first_level_merge = false;
            self.last_merge_red_edges_erased.shrinkRetainingCapacity(0);

            self.erased_nodes.set(erased);
            if (self.node_list[survivor].isLeaf()) {
                const parent = self.node_list[survivor].getFirstNeighboor();
                self.node_list[parent].num_leafes -= 1;
            }
            if (self.node_list[erased].isLeaf()) {
                const parent = self.node_list[erased].getFirstNeighboor();
                self.node_list[parent].num_leafes -= 1;
            }

            var red_iter = self.node_list[erased].red_edges.iterator();
            while (red_iter.next()) |item| {
                try self.node_list[item].removeRedEdge(self.allocator, erased);
                if (item != survivor) {
                    if (!try self.node_list[survivor].addRedEdgeExists(self.allocator, item)) {
                        try self.node_list[item].addRedEdge(self.allocator, survivor);
                        self.last_merge_red_edges_erased.append(self.allocator, item) catch unreachable;
                        try self.min_hash.changedEdge(item, self, .{ .removed = erased, .red = true, .added = survivor, .added_red = true });
                    } else {
                        // Inform about the removal of the red edge
                        try seq.red_edge_stack.addEdge(self.failing_allocator.allocator(), red_edge_stack.NewRedEdge(T).redToDeleted(item));
                        try self.min_hash.changedEdge(item, self, .{ .removed = erased, .red = true });
                    }
                } else {
                    self.last_merge_first_level_merge = true;
                }
            }

            var tww: T = 0;
            var black_iter = self.node_list[erased].black_edges.xorIterator(&self.node_list[survivor].black_edges);

            var remove_list = std.ArrayList(T).init(self.allocator);
            defer remove_list.deinit();

            while (black_iter.next()) |item| {
                if (item == survivor or item == erased) {
                    self.last_merge_first_level_merge = true;
                    continue;
                }
                // Came from erased
                if (black_iter.first) {
                    if (!try self.node_list[survivor].addRedEdgeExists(self.allocator, item)) {
                        try seq.red_edge_stack.addEdge(self.failing_allocator.allocator(), red_edge_stack.NewRedEdge(T).blackToRedOther(item));

                        try self.node_list[item].addRedEdge(self.allocator, survivor);

                        self.last_merge_red_edges_erased.append(self.allocator, item) catch unreachable;
                        //self.min_hash.addTransferedEdge(item,true,false);
                        try self.min_hash.changedEdge(item, self, .{ .removed = erased, .red = false, .added = survivor, .added_red = true });
                    } else {
                        try self.min_hash.changedEdge(item, self, .{ .removed = erased, .red = false });
                    }
                }
                // Came from survivor
                else {
                    if (try self.node_list[survivor].addRedEdgeExists(self.allocator, item)) {
                        // If it existed before we inherited from erased

                        try seq.red_edge_stack.addEdge(self.failing_allocator.allocator(), red_edge_stack.NewRedEdge(T).blackToDeleted(item));

                        try self.min_hash.changedEdge(item, self, .{ .removed = survivor, .red = false });
                    } else {
                        // Did not exist therefore turned
                        try seq.red_edge_stack.addEdge(self.failing_allocator.allocator(), red_edge_stack.NewRedEdge(T).blackToRedOwn(item));
                        try self.min_hash.changedEdge(item, self, .{ .removed = survivor, .red = false, .added = survivor, .added_red = true });
                    }
                    try self.node_list[item].addRedEdge(self.allocator, survivor);

                    try self.node_list[item].removeBlackEdge(self.allocator, survivor);
                    //self.min_hash.addChangeColorSurvivorEdge(item);
                    try remove_list.append(item);
                }
                tww = std.math.max(tww, @intCast(T, self.node_list[item].red_edges.cardinality()));
            }

            for (remove_list.items) |item| {
                // Batch remove?
                try self.node_list[survivor].removeBlackEdge(self.allocator, item);
            }

            var black_remove_iter = self.node_list[erased].black_edges.iterator();
            while (black_remove_iter.next()) |item| {
                try self.node_list[item].removeBlackEdge(self.allocator, erased);
                //self.min_hash.addRemovedEdge(item,true);
            }

            var black_iter_sur = self.node_list[survivor].black_edges.iterator();
            while (black_iter_sur.next()) |t| {
                try self.min_hash.changedEdge(t, self, .{ .removed = erased, .red = false });
            }

            try self.min_hash.rehashNode(survivor, self);
            try self.min_hash.removeNode(erased);

            tww = std.math.max(tww, @intCast(T, self.node_list[survivor].red_edges.cardinality()));
            try seq.addContraction(self.allocator, erased, survivor, std.math.max(tww, seq.getTwinWidth()));

            self.updateLeafCount(survivor);
            return tww;
        }

        pub fn forceExactSolverToProduceSolution(self: *Self) void {
            std.debug.print("Force exact solver to produce solution", .{});
            self.force_exact_solver_to_solve = true;
        }

        pub fn solveExact(self: *Self) !T {
            var org_graph = try std.ArrayList(exact_bs.MatrixGraphFromInducedSubGraph).initCapacity(self.allocator, self.connected_components.items.len);
            defer org_graph.deinit();
            defer for (org_graph.items) |*item| item.deinit();

            for (self.connected_components.items) |cc| {
                try org_graph.append(try exact_bs.matrixGraphUnionFromInducedSubGraph(T, &cc.subgraph, self.allocator));
            }

            var heu_tww = @intCast(u32, try self.solveGreedy(30));

            if (true) {
                std.debug.print("Heuristic TWW: {d}\n", .{heu_tww});
                for (self.connected_components.items) |*cc| {
                    if (cc.subgraph.nodes.len > 10) {
                        std.debug.print("  -> CC: n={d} tww={d}\n", .{ cc.subgraph.nodes.len, cc.tww });
                    }
                }
            }

            var lower: T = 0;
            while (true) {
                var upper_bound: u32 = 0;
                for (self.connected_components.items) |cc| {
                    if (cc.tww > heu_tww) {
                        return GraphError.ExactSuboptimal;
                    }
                    upper_bound = @max(upper_bound, cc.tww);
                }

                if (upper_bound <= lower) {
                    return lower;
                }

                var cc_iter = self.connected_components_min_heap.iterator();
                while (cc_iter.next()) |cc| {
                    var component = &self.connected_components.items[cc.index];
                    if (component.tww < upper_bound) {
                        if (component.subgraph.nodes.len > 10) {
                            std.debug.print("Skip CC for the moment n={d} tww={d} | lower={d} upper={d}\n", .{ component.subgraph.nodes.len, component.tww, lower, upper_bound });
                        }
                        continue;
                    }
                    var subgraph = &org_graph.items[cc.index];

                    // the exact solver treats the upper as exclusive; i.e. by setting it to the heuristic tww,
                    // we force the solver to produce a better solution or fail.
                    std.debug.print("Invoke exact solver with |CC|={d} lower={d} upper={d}\n", .{ component.subgraph.nodes.len, lower, upper_bound });
                    var improved_result = solver_bnb.solveCCExactly(T, component, self.allocator, subgraph, lower, upper_bound + @boolToInt(self.force_exact_solver_to_solve)) catch |e| {
                        if (e == solver_bnb.SolverError.Infeasable) {
                            std.debug.print(" ... infeasable\n", .{});

                            if (self.force_exact_solver_to_solve) {
                                return GraphError.ExactSuboptimal;
                            }

                            lower = @intCast(T, upper_bound);
                            break;
                        }
                        std.debug.print(" ... failed\n", .{});
                        return e;
                    };
                    std.debug.print(" ... solve with tww={d}\n", .{improved_result});

                    lower = @max(lower, @intCast(T, improved_result));
                }
            }

            try self.combineContractionSequences();

            return lower;
        }

        pub fn combineContractionSequences(self: *Self) !void {
            self.contraction.reset();

            var cc_iter = self.connected_components_min_heap.iterator();
            var tww: T = 0;
            _ = tww;
            var last_node: ?T = null;

            while (cc_iter.next()) |cc| {
                if (self.connected_components.items[cc.index].subgraph.nodes.len == 1) {
                    const survivor = self.connected_components.items[cc.index].subgraph.nodes[0];
                    if (last_node) |ln| {
                        try self.contraction.addContraction(contraction.Contraction(T){
                            .erased = survivor,
                            .survivor = ln,
                        });
                    } else {
                        last_node = survivor;
                    }
                } else if (self.connected_components.items[cc.index].best_contraction_sequence.getLastContraction()) |ctr| {
                    try self.contraction.append(&self.connected_components.items[cc.index].best_contraction_sequence);
                    if (last_node) |ln| {
                        try self.contraction.addContraction(contraction.Contraction(T){
                            .erased = ctr.survivor,
                            .survivor = ln,
                        });
                    } else {
                        last_node = ctr.survivor;
                    }
                }
            }
        }

        pub fn solveGreedy(self: *Self, timeout_seconds: ?u64) !T {
            const K = 20;
            const P = 100;

            // NOTICE: This function is single pass at the moment!
            var solver = try solver_resources.SolverResources(T, K, P).init(self);
            defer solver.deinit(self.allocator);

						var time = try std.time.Instant.now();

						var seed:u64 = 19;

            self.contraction.reset();
            while (self.connected_components_min_heap.removeOrNull()) |cc| {
                // Only for small exact graphs

                var cc_inst = &self.connected_components.items[cc.index];
								try cc_inst.resetGraph();

                if (T == u8 and cc_inst.subgraph.nodes.len < 128) {
                    var cc_tww = try cc_inst.solveGreedy(K, P, &solver,seed);
                    cc_inst.tww = cc_tww;
                } else {
                    var cc_tww = try cc_inst.solveGreedyTopK(K, P, &solver,seed);
                    cc_inst.tww = cc_tww;
                }
								seed += 1009;
								try self.connected_components_min_heap.add(connected_components.ConnectedComponentIndex(T) {
									.index = cc.index,
									.tww = cc_inst.tww,
								});

								if(timeout_seconds) |ts| {
									const now = try std.time.Instant.now();
									if(now.since(time) > 1000_000_000*ts) {
										break;
									}
								}
            }

						const index = self.connected_components_min_heap.peek().?;
						const tww = self.connected_components.items[index.index].tww;

            try self.combineContractionSequences();

            return tww;
        }

        pub fn new(number_of_nodes: T, allocator: std.mem.Allocator) !Self {
            var node_list = try allocator.alloc(NodeType, number_of_nodes);

            for (node_list) |*node| {
                node.black_edges = compressed_bitset.FastCompressedBitmap(T, promote_thresh, degrade_tresh).init(number_of_nodes);
                node.red_edges = compressed_bitset.FastCompressedBitmap(T, promote_thresh, degrade_tresh).init(number_of_nodes);
                node.high_degree_node = null;
                node.num_leafes = 0;
            }

            //TODO: Add some errdefer's here

            var graph = Self{
                .number_of_nodes = number_of_nodes,
                .number_of_edges = 0,
                .node_list = node_list,
                .allocator = allocator,
                .contraction = try contraction.ContractionSequence(T).init(allocator, number_of_nodes),
                .scratch_bitset = try bitset.FastBitSet.initEmpty(number_of_nodes, allocator),
                .erased_nodes = try bitset.FastBitSet.initEmpty(number_of_nodes, allocator),
                .connected_components = std.ArrayListUnmanaged(connected_components.ConnectedComponent(T)){},
                .trivial_connected_component_contraction_sequence = std.ArrayListUnmanaged(u32){},
                .connected_components_min_heap = std.PriorityQueue(connected_components.ConnectedComponentIndex(T), void, connected_components.ConnectedComponentIndex(T).compareComponentIndexDesc).init(allocator, {}),
                .failing_allocator = std.heap.FixedBufferAllocator.init(&[_]u8{}),
                .connected_components_node_list_slice = try allocator.alloc(T, number_of_nodes),
                .started_at = try std.time.Instant.now(),
                .min_hash = undefined,
                .last_merge_first_level_merge = false,
                .last_merge_red_edges_erased = try std.ArrayListUnmanaged(T).initCapacity(allocator, number_of_nodes),
            };
            graph.min_hash = try min_hash_mod.MinHashSimiliarity(T, 1).init(allocator, 36, graph.number_of_nodes);

            return graph;
        }

        pub inline fn density(self: *const Self) f64 {
            return (2 * @intToFloat(f64, self.number_of_edges)) / (@intToFloat(f64, self.number_of_nodes) * @intToFloat(f64, self.number_of_nodes - 1));
        }

        pub fn loadFromPace(allocator: std.mem.Allocator, pace: *pace_2023.Pace2023Fmt(T)) !Self {
            var node_list = try allocator.alloc(NodeType, pace.number_of_nodes);

            for (0..pace.number_of_nodes) |index| {
                node_list[index].black_edges = try compressed_bitset.FastCompressedBitmap(T, promote_thresh, degrade_tresh).fromUnsorted(allocator, &pace.nodes[index].edges, @intCast(T, pace.number_of_nodes));
                node_list[index].red_edges = compressed_bitset.FastCompressedBitmap(T, promote_thresh, degrade_tresh).init(@intCast(T, pace.number_of_nodes));
                node_list[index].high_degree_node = null;
                node_list[index].num_leafes = 0;
            }

            for (0..pace.number_of_nodes) |index| {
                if (node_list[index].cardinality() == 1) {
                    const parent = node_list[index].getFirstNeighboor();
                    node_list[parent].num_leafes += 1;
                }
            }

            // Remove allocations on failure
            errdefer {
                for (node_list) |*node| {
                    node.black_edges.deinit(allocator);
                    node.red_edges.deinit(allocator);
                    node.high_degree_node = null;
                    node.total_weight = @intToFloat(f32, node.cardinality());
                    node.delta_potential_weighted_jaccard = 0.0;
                }
                allocator.free(node_list);
            }

            var graph_instance = Self{
                .number_of_nodes = pace.number_of_nodes,
                .number_of_edges = pace.number_of_edges, //ATM
                .node_list = node_list,
                .allocator = allocator,
                .contraction = try contraction.ContractionSequence(T).init(allocator, pace.number_of_nodes),
                .erased_nodes = try bitset.FastBitSet.initEmpty(pace.number_of_nodes, allocator),
                .scratch_bitset = try bitset.FastBitSet.initEmpty(pace.number_of_nodes, allocator),
                .connected_components = std.ArrayListUnmanaged(connected_components.ConnectedComponent(T)){},
                .trivial_connected_component_contraction_sequence = std.ArrayListUnmanaged(u32){},
                .connected_components_node_list_slice = try allocator.alloc(T, pace.number_of_nodes),
                .connected_components_min_heap = std.PriorityQueue(connected_components.ConnectedComponentIndex(T), void, connected_components.ConnectedComponentIndex(T).compareComponentIndexDesc).init(allocator, {}),
                .failing_allocator = std.heap.FixedBufferAllocator.init(&[_]u8{}),
                .started_at = try std.time.Instant.now(),
                .min_hash = undefined,
                .last_merge_first_level_merge = false,
                .last_merge_red_edges_erased = try std.ArrayListUnmanaged(T).initCapacity(allocator, pace.number_of_nodes),
                .force_exact_solver_to_solve = false,
            };

            var hash = try min_hash_mod.MinHashSimilarity(T, 7).init(graph_instance.allocator, 120, graph_instance.number_of_nodes);
            graph_instance.min_hash = hash;
            return graph_instance;
        }

        pub fn findAllConnectedComponents(self: *Self) ![][]u32 {
            self.scratch_bitset.unsetAll();

            var bfs_stack = try bfs_mod.BfsQueue(T).init(self.allocator, self.number_of_nodes);
            defer bfs_stack.deinit(self.allocator);

            var unsetIter = self.scratch_bitset.iterUnset();
            var non_trivial_components: u32 = 0;

            var current_slice_start: u32 = 0;
            var current_slice_ptr: u32 = 0;

            while (unsetIter.next()) |item| {
                if (self.scratch_bitset.get(item)) {
                    continue;
                }
                var iterator = bfs_mod.bfs(T, @intCast(T, item), self, &self.scratch_bitset, &bfs_stack, .{ .max_level = std.math.maxInt(T), .kind = .black });

                while (iterator.next()) |node| {
                    self.scratch_bitset.set(node);
                    self.connected_components_node_list_slice[current_slice_ptr] = node;
                    current_slice_ptr += 1;
                }

                var tmp = self.connected_components_node_list_slice[current_slice_start..current_slice_ptr];
                if (tmp.len < 2) { // trivial component
                    for (tmp) |node| {
                        try self.trivial_connected_component_contraction_sequence.append(self.allocator, @as(u32, node));
                    }
                } else {
                    non_trivial_components += 1;
                    try self.connected_components.append(self.allocator, try connected_components.ConnectedComponent(T).init(self.allocator, tmp, iterator.level, self));
                }
                current_slice_start = current_slice_ptr;
            }

            // construct solution tracking
            var cc_solutions = try self.allocator.alloc([]u32, non_trivial_components + 1);
            for (0..non_trivial_components) |i| {
                cc_solutions[i] = self.connected_components.items[i].contraction_slice;
            }
            const items = self.trivial_connected_component_contraction_sequence.items;
            const last_index = non_trivial_components;
            if (items.len == 0) {
                cc_solutions[last_index] = &[_]u32{};
            } else if (items.len == 1) {
                cc_solutions[last_index] = try self.allocator.alloc(u32, 1);
                cc_solutions[last_index][0] = self.trivial_connected_component_contraction_sequence.items[0];
            } else {
                cc_solutions[last_index] = try self.allocator.alloc(u32, (items.len - 1) * 2);
                var j: usize = 0;
                for (0..(items.len - 1)) |i| {
                    cc_solutions[last_index][j] = items[i + 1];
                    cc_solutions[last_index][j + 1] = items[i];
                    j += 2;
                }
            }

            try self.connected_components_min_heap.ensureTotalCapacity(self.connected_components.items.len);
            for (0..self.connected_components.items.len) |index| {
                try self.connected_components_min_heap.add(connected_components.ConnectedComponentIndex(T){ .tww = self.connected_components.items[index].tww, .index = @intCast(T, index) });
            }

            var tww = self.connected_components_min_heap.remove();
            //std.debug.print("Found {} components largest {} and tww {} density {}\n", .{ components, largest, tww.tww, self.density() });
            try self.connected_components_min_heap.add(tww);
            return cc_solutions;
        }

        pub fn deinit(self: *Self) void {
            for (self.node_list) |*node| {
                node.black_edges.deinit(self.allocator);
                node.red_edges.deinit(self.allocator);
            }
            for (self.connected_components.items) |*connected_component| {
                connected_component.deinit(self.allocator);
            }
            self.connected_components.deinit(self.allocator);
            self.connected_components_min_heap.deinit();
            self.trivial_connected_component_contraction_sequence.deinit(self.allocator);

            self.contraction.deinit(self.allocator);
            self.allocator.free(self.connected_components_node_list_slice);
            self.scratch_bitset.deinit(self.allocator);
            self.erased_nodes.deinit(self.allocator);
            self.last_merge_red_edges_erased.deinit(self.allocator);
            self.allocator.free(self.node_list);
        }
    };
}

test "Check simple loading" {
    var allocator = std.heap.GeneralPurposeAllocator(.{}){};
    defer std.debug.assert(!allocator.deinit());

    var pace_fmt = try pace_2023.Pace2023Fmt(u8).fromFile(allocator.allocator(), "instances/tiny/tiny001.gr");
    defer pace_fmt.deinit(allocator.allocator());
    var graph = try Graph(u8).loadFromPace(allocator.allocator(), &pace_fmt);
    // Free ressources
    defer graph.deinit();

    try std.testing.expectEqual(graph.number_of_nodes, 10);
    try std.testing.expectEqual(graph.number_of_edges, 9);
}

test "Test contraction Tiny 1" {
    var allocator = std.heap.GeneralPurposeAllocator(.{}){};
    defer std.debug.assert(!allocator.deinit());

    var pace_fmt = try pace_2023.Pace2023Fmt(u8).fromFile(allocator.allocator(), "instances/tiny/tiny001.gr");
    defer pace_fmt.deinit(allocator.allocator());

    var ret = try RetraceableContractionSequence(u8).init(allocator.allocator(), 10, 9);
    defer ret.deinit(allocator.allocator());

    var graph = try Graph(u8).loadFromPace(allocator.allocator(), &pace_fmt);
    // Free ressources
    defer graph.deinit();

    try std.testing.expectEqual(graph.number_of_nodes, 10);
    try std.testing.expectEqual(graph.number_of_edges, 9);
    try std.testing.expectEqual(graph.node_list[8].num_leafes, 1);
    try std.testing.expectEqual(graph.node_list[1].num_leafes, 1);

    var i: u8 = 9;
    while (i > 0) {
        if (i > 2) {
            try std.testing.expectEqual(graph.node_list[i - 1].num_leafes, 1);
        } else if (i == 2) {
            try std.testing.expectEqual(graph.node_list[i - 1].num_leafes, 2);
        } else if (i == 1) {
            try std.testing.expectEqual(graph.node_list[9].num_leafes, 1);
            try std.testing.expectEqual(graph.node_list[0].num_leafes, 1);
        }
        var tww = try graph.addContraction(i - 1, 9, &ret);
        if (i > 1) {
            try std.testing.expectEqual(@as(u32, 1), tww);
        } else {
            try std.testing.expectEqual(@as(u32, 0), tww);
        }
        i -= 1;
    }
}

test "Test contraction Tiny 2" {
    var allocator = std.heap.GeneralPurposeAllocator(.{}){};
    defer std.debug.assert(!allocator.deinit());

    var pace_fmt = try pace_2023.Pace2023Fmt(u8).fromFile(allocator.allocator(), "instances/tiny/tiny002.gr");
    defer pace_fmt.deinit(allocator.allocator());

    var ret = try RetraceableContractionSequence(u8).init(allocator.allocator(), 10, 10);
    defer ret.deinit(allocator.allocator());

    var graph = try Graph(u8).loadFromPace(allocator.allocator(), &pace_fmt);
    // Free ressources
    defer graph.deinit();

    try std.testing.expectEqual(graph.number_of_nodes, 10);
    try std.testing.expectEqual(graph.number_of_edges, 10);

    var i: u8 = 9;
    while (i > 0) {
        var tww = try graph.addContraction(i - 1, 9, &ret);
        if (i > 2) {
            try std.testing.expectEqual(@as(u32, 2), tww);
        } else if (i > 1) {
            try std.testing.expectEqual(@as(u32, 1), tww);
        } else {
            try std.testing.expectEqual(@as(u32, 0), tww);
        }
        i -= 1;
    }
}

test "Test contraction custom" {
    var allocator = std.heap.GeneralPurposeAllocator(.{}){};
    defer std.debug.assert(!allocator.deinit());

    var ret = try RetraceableContractionSequence(u8).init(allocator.allocator(), 10, 10);
    defer ret.deinit(allocator.allocator());

    var graph = try Graph(u8).new(10, allocator.allocator());
    // Free ressources
    defer graph.deinit();

    try std.testing.expectEqual(graph.number_of_nodes, 10);

    try graph.addEdge(0, 1);
    try graph.addEdge(1, 2);
    try graph.addEdge(2, 3);
    try graph.addEdge(3, 4);
    try graph.addEdge(4, 5);
    try graph.addEdge(5, 6);
    _ = try graph.addContraction(0, 6, &ret);
    var i: u8 = 9;
    while (i > 0) {
        try std.testing.expectEqual(graph.node_list[i].num_leafes, 0);
        i -= 1;
    }
}

test "Test contraction retrace Tiny 2" {
    var allocator = std.heap.GeneralPurposeAllocator(.{}){};
    defer std.debug.assert(!allocator.deinit());

    var pace_fmt = try pace_2023.Pace2023Fmt(u8).fromFile(allocator.allocator(), "instances/tiny/tiny002.gr");
    defer pace_fmt.deinit(allocator.allocator());

    var ret = try RetraceableContractionSequence(u8).init(allocator.allocator(), 10, 10);
    defer ret.deinit(allocator.allocator());

    var graph = try Graph(u8).loadFromPace(allocator.allocator(), &pace_fmt);
    // Free ressources
    defer graph.deinit();

    try std.testing.expectEqual(graph.number_of_nodes, 10);
    try std.testing.expectEqual(graph.number_of_edges, 10);

    var i: u8 = 9;
    while (i > 0) {
        if (i == 1) {
            try std.testing.expectEqual(graph.node_list[9].num_leafes, 1);
            try std.testing.expectEqual(graph.node_list[0].num_leafes, 1);
        }
        var tww = try graph.addContraction(i - 1, 9, &ret);
        if (i > 2) {
            try std.testing.expectEqual(@as(u32, 2), tww);
        } else if (i > 1) {
            try std.testing.expectEqual(@as(u32, 1), tww);
        } else {
            try std.testing.expectEqual(@as(u32, 0), tww);
        }
        i -= 1;
    }
    i = 9;
    while (i > 0) {
        try graph.revertLastContraction(&ret);

        const tww = ret.getTwinWidth();

        if (i == 9) {
            try std.testing.expectEqual(@as(u32, 2), tww);
        } else if (i == 1) {
            try std.testing.expectEqual(@as(u32, 0), tww);
        } else {
            try std.testing.expectEqual(@as(u32, 2), tww);
        }
        for (graph.node_list) |*node| {
            try std.testing.expect(node.red_edges.cardinality() <= tww);
        }
        i -= 1;
    }
    var count: u32 = 0;
    for (graph.node_list) |node| {
        try std.testing.expectEqual(node.red_edges.cardinality(), 0);
        count += node.black_edges.cardinality();
    }
    try std.testing.expectEqual(count, 20);
}

test "Test contraction Tiny 3" {
    var allocator = std.heap.GeneralPurposeAllocator(.{}){};
    defer std.debug.assert(!allocator.deinit());

    var pace_fmt = try pace_2023.Pace2023Fmt(u8).fromFile(allocator.allocator(), "instances/tiny/tiny003.gr");
    defer pace_fmt.deinit(allocator.allocator());

    var ret = try RetraceableContractionSequence(u8).init(allocator.allocator(), 10, 10);
    defer ret.deinit(allocator.allocator());

    var graph = try Graph(u8).loadFromPace(allocator.allocator(), &pace_fmt);
    // Free ressources
    defer graph.deinit();

    try std.testing.expectEqual(graph.number_of_nodes, 10);
    try std.testing.expectEqual(graph.number_of_edges, 45);

    var i: u8 = 9;
    while (i > 0) {
        var tww = try graph.addContraction(i - 1, 9, &ret);
        try std.testing.expectEqual(@as(u32, 0), tww);
        i -= 1;
    }
}

test "Test contraction retrace Tiny 3" {
    var allocator = std.heap.GeneralPurposeAllocator(.{}){};
    defer std.debug.assert(!allocator.deinit());

    var pace_fmt = try pace_2023.Pace2023Fmt(u8).fromFile(allocator.allocator(), "instances/tiny/tiny003.gr");
    defer pace_fmt.deinit(allocator.allocator());

    var ret = try RetraceableContractionSequence(u8).init(allocator.allocator(), 10, 45);
    defer ret.deinit(allocator.allocator());
    var graph = try Graph(u8).loadFromPace(allocator.allocator(), &pace_fmt);
    // Free ressources
    defer graph.deinit();

    try std.testing.expectEqual(graph.number_of_nodes, 10);
    try std.testing.expectEqual(graph.number_of_edges, 45);

    var i: u8 = 9;
    while (i > 0) {
        var tww = try graph.addContraction(i - 1, 9, &ret);
        try std.testing.expectEqual(@as(u32, 0), tww);
        i -= 1;
    }
    i = 9;
    while (i > 0) {
        try graph.revertLastContraction(&ret);
        var tww = graph.getCurrentTwinWidth();
        try std.testing.expectEqual(@as(u32, 0), tww);
        for (graph.node_list) |*node| {
            try std.testing.expect(node.red_edges.cardinality() <= tww);
        }
        i -= 1;
    }

    i = 9;
    while (i > 0) {
        var tww = try graph.addContraction(i - 1, 9, &ret);
        try std.testing.expectEqual(@as(u32, 0), tww);
        i -= 1;
    }
}

test "Test contraction Tiny 4" {
    var allocator = std.heap.GeneralPurposeAllocator(.{}){};
    defer std.debug.assert(!allocator.deinit());

    var pace_fmt = try pace_2023.Pace2023Fmt(u8).fromFile(allocator.allocator(), "instances/tiny/tiny004.gr");
    defer pace_fmt.deinit(allocator.allocator());

    var ret = try RetraceableContractionSequence(u8).init(allocator.allocator(), 10, 9);
    defer ret.deinit(allocator.allocator());

    var graph = try Graph(u8).loadFromPace(allocator.allocator(), &pace_fmt);
    // Free ressources
    defer graph.deinit();

    try std.testing.expectEqual(graph.number_of_nodes, 10);
    try std.testing.expectEqual(graph.number_of_edges, 9);

    var i: u8 = 9;
    while (i > 0) {
        try std.testing.expectEqual(graph.node_list[0].num_leafes, i);
        var tww = try graph.addContraction(i - 1, 9, &ret);
        try std.testing.expectEqual(@as(u32, 0), tww);
        i -= 1;
    }
}

test "Test contraction retrace Tiny 4" {
    var allocator = std.heap.GeneralPurposeAllocator(.{}){};
    defer std.debug.assert(!allocator.deinit());

    var pace_fmt = try pace_2023.Pace2023Fmt(u8).fromFile(allocator.allocator(), "instances/tiny/tiny004.gr");
    defer pace_fmt.deinit(allocator.allocator());
    var ret = try RetraceableContractionSequence(u8).init(allocator.allocator(), 10, 9);
    defer ret.deinit(allocator.allocator());
    var graph = try Graph(u8).loadFromPace(allocator.allocator(), &pace_fmt);
    // Free ressources
    defer graph.deinit();

    try std.testing.expectEqual(graph.number_of_nodes, 10);
    try std.testing.expectEqual(graph.number_of_edges, 9);

    var i: u8 = 9;
    while (i > 0) {
        var tww = try graph.addContraction(i - 1, 9, &ret);
        try std.testing.expectEqual(@as(u32, 0), tww);
        i -= 1;
    }
    i = 9;
    while (i > 0) {
        try graph.revertLastContraction(&ret);
        var tww = graph.getCurrentTwinWidth();
        try std.testing.expectEqual(@as(u32, 0), tww);
        for (graph.node_list) |*node| {
            try std.testing.expect(node.red_edges.cardinality() <= tww);
        }
        i -= 1;
    }
}

test "Test contraction Tiny 6" {
    var allocator = std.heap.GeneralPurposeAllocator(.{}){};
    defer std.debug.assert(!allocator.deinit());

    var pace_fmt = try pace_2023.Pace2023Fmt(u8).fromFile(allocator.allocator(), "instances/tiny/tiny006.gr");
    defer pace_fmt.deinit(allocator.allocator());

    var ret = try RetraceableContractionSequence(u8).init(allocator.allocator(), 10, 5);
    defer ret.deinit(allocator.allocator());
    var graph = try Graph(u8).loadFromPace(allocator.allocator(), &pace_fmt);
    // Free ressources
    defer graph.deinit();

    try std.testing.expectEqual(graph.number_of_nodes, 10);
    try std.testing.expectEqual(graph.number_of_edges, 5);

    var i: u8 = 9;
    var tww: u32 = 0;
    while (true) {
        tww = std.math.max(tww, try graph.addContraction(i - 1, i, &ret));
        if (i == 1) {
            break;
        }
        i -= 2;
    }
    try std.testing.expectEqual(@as(u32, 0), tww);
}

test "Check simple failed loading" {
    var allocator = std.heap.GeneralPurposeAllocator(.{}){};
    // Should never fail since the ressources should be managed
    defer std.debug.assert(!allocator.deinit());

    var pace_fmt = pace_2023.Pace2023Fmt(u8).fromFile(allocator.allocator(), "instances/tiny/tiny001.grfailure") catch |err| {
        try std.testing.expectEqual(err, error.FileNotFound);
        return;
    };

    var graph = try Graph(u8).loadFromPace(allocator.allocator(), &pace_fmt);
    _ = graph;
    // Should never be reached
    std.debug.assert(false);
}

test "Check simple graph" {
    var allocator = std.heap.GeneralPurposeAllocator(.{}){};
    // Should never fail since the ressources should be managed
    //defer std.debug.assert(!allocator.deinit());

    var graph = try Graph(u8).new(5, allocator.allocator());

    var ret = try RetraceableContractionSequence(u8).init(allocator.allocator(), 10, 4);
    defer ret.deinit(allocator.allocator());

    try graph.addEdge(0, 1);
    try graph.addEdge(1, 2);
    try graph.addEdge(1, 3);
    try graph.addEdge(1, 4);

    _ = try graph.addContraction(1, 0, &ret);

    try std.testing.expectEqual(graph.node_list[0].red_edges.cardinality(), 3);
}
