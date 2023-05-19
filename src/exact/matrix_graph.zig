const std = @import("std");
const print = std.debug.print;
const FixedSizeSet = @import("fixed_size_bitset.zig").FixedSizeSet;
const builtin = @import("builtin");
const assert = std.debug.assert;

const ExpensiveConsistencyChecks = false;

pub const Color = enum {
    Black,
    Red,

    pub fn isRed(self: Color) bool {
        return self == Color.Red;
    }

    pub fn isBlack(self: Color) bool {
        return self == Color.Black;
    }
};

pub const RowType = enum(usize) {
    Neighbor = 0,
    Red = 1,
    Two = 2,
};

pub const DigestAlgo = std.crypto.hash.Sha1;
pub const Digest = [DigestAlgo.digest_length]u8;

pub fn MatrixGraph(comptime num_nodes: u32) type {
    return struct {
        const Self = @This();
        pub const NumNodes: u32 = num_nodes;
        pub const Node = u32;
        pub const BitSet = FixedSizeSet(num_nodes);

        allocator: std.mem.Allocator,
        num_edges: u32,
        matrix: []BitSet,
        has_neighbors: BitSet,
        invalid_two_neighbors: BitSet,
        batch_updates: bool,

        pub fn new(allocator: std.mem.Allocator) !Self {
            var matrix = try allocator.alloc(BitSet, 3 * num_nodes);
            errdefer allocator.free(matrix);

            const empty = FixedSizeSet(num_nodes).new();
            for (matrix) |*r| r.* = empty;

            return Self{
                .allocator = allocator,
                .num_edges = 0,
                .matrix = matrix,
                .has_neighbors = BitSet.new(),
                .invalid_two_neighbors = BitSet.newAllSet(),
                .batch_updates = false,
            };
        }

        pub fn deinit(self: *const Self) void {
            self.allocator.free(self.matrix);
        }

        pub fn copy(self: *const Self) !Self {
            var result = try Self.new(self.allocator);
            result.copyFrom(self);
            return result;
        }

        pub fn copyFrom(self: *Self, other: *const Self) void {
            self.num_edges = other.num_edges;
            self.has_neighbors = other.has_neighbors;
            self.invalid_two_neighbors = other.invalid_two_neighbors;
            self.batch_updates = other.batch_updates;
            std.mem.copy(BitSet, self.matrix, other.matrix);
        }

        pub inline fn numberOfNodes(self: *const Self) Node {
            _ = self;
            return num_nodes;
        }

        pub inline fn numberOfEdges(self: *const Self) u32 {
            return self.num_edges;
        }

        pub inline fn deg(self: *const Self, u: Node) Node {
            return @as(Node, self.constNeighbors(u).cardinality());
        }

        pub inline fn redDeg(self: *const Self, u: Node) Node {
            return @as(Node, self.constRedNeighbors(u).cardinality());
        }

        pub inline fn blackDeg(self: *const Self, u: Node) Node {
            return self.deg(u) - self.redDeg(u);
        }

        pub fn beginBatchUpdates(self: *Self) void {
            self.batch_updates = true;
        }

        pub fn endBatchUpdates(self: *Self) void {
            assert(self.batch_updates);
            self.batch_updates = false;
            var u: Node = 0;
            while (u < num_nodes) : (u += 1) {
                self.recomputeTwoNeighbors(u);
            }
            self.assertIsConsistent();
        }

        pub fn addEdge(self: *Self, u: Node, v: Node, color: Color) ?Color {
            var prev: ?Color = null;

            if (self.neighbors(u).setBit(v)) {
                prev = Color.Black;
            } else {
                _ = self.neighbors(v).setBit(u);
                self.num_edges += 1;

                inline for ([_]Node{ u, v }, [_]Node{ v, u }) |a, b| {
                    _ = self.has_neighbors.setBit(a);

                    if (!self.batch_updates) {
                        _ = self.twoNeighbors(a).setBit(b);
                        self.setColumn(RowType.Two, a, self.constNeighbors(b).iter_set());
                        self.twoNeighbors(a).assignOr(self.constNeighbors(b));
                    }
                }
            }

            if (color.isRed()) {
                if (self.redNeighbors(u).setBit(v)) {
                    prev = Color.Red;
                }
                _ = self.redNeighbors(v).setBit(u);
            } else {
                if (self.redNeighbors(u).unsetBit(v)) {
                    prev = Color.Red;
                    _ = self.redNeighbors(v).unsetBit(u);
                }
            }

            self.assertIsConsistent();

            return prev;
        }

        pub fn removeEdge(self: *Self, u: Node, v: Node) ?Color {
            if (!self.neighbors(u).unsetBit(v)) {
                return null;
            }

            var prev = Color.Black;
            _ = self.neighbors(v).unsetBit(u);

            if (self.redNeighbors(u).unsetBit(v)) {
                _ = self.redNeighbors(v).unsetBit(u);
                prev = Color.Red;
            }

            self.num_edges -= 1;

            if (self.neighbors(u).areAllUnset()) {
                _ = self.has_neighbors.unsetBit(u);
            }

            if (self.neighbors(v).areAllUnset()) {
                _ = self.has_neighbors.unsetBit(v);
            }

            if (!self.batch_updates and u != v) {
                var nodes_to_update = self.constNeighbors(u).copyWithOr(self.constNeighbors(v));
                _ = nodes_to_update.setBit(u);
                _ = nodes_to_update.setBit(v);

                var iter = nodes_to_update.iter_set();
                while (iter.next()) |i| {
                    self.recomputeTwoNeighbors(i);
                }
            }

            self.assertIsConsistent();

            return prev;
        }

        pub fn removeEdgesAtNode(self: *Self, u: Node) void {
            // neigbors
            {
                var neighs = self.neighbors(u).iter_set();
                while (neighs.next()) |v| {
                    _ = self.neighbors(v).unsetBit(u);

                    if (self.neighbors(v).areAllUnset()) {
                        _ = self.has_neighbors.unsetBit(v);
                    }
                }
            }

            // red neighbors
            {
                var neighs = self.redNeighbors(u).iter_set();
                while (neighs.next()) |v| {
                    _ = self.redNeighbors(v).unsetBit(u);
                }
            }

            // two neighbors
            if (!self.batch_updates) {
                var iter = self.twoNeighbors(u).iter_set();
                while (iter.next()) |i| {
                    _ = self.twoNeighbors(i).unsetBit(u);
                }
                self.twoNeighbors(u).unsetAll();
            }

            // neighbors
            if (!self.batch_updates) {
                var neighs = self.neighbors(u).iter_set();
                while (neighs.next()) |v| {
                    self.recomputeTwoNeighbors(v);
                }
            }

            self.num_edges -= self.neighbors(u).cardinality();
            self.neighbors(u).unsetAll();
            self.redNeighbors(u).unsetAll();
            _ = self.has_neighbors.unsetBit(u);

            self.assertIsConsistent();
        }

        pub fn redNeighborsAfterMerge(self: *const Self, rem: Node, sur: Node) BitSet {
            // compute currently black neighbors
            var result: BitSet = self.blackNeighborsCopied(sur);
            var black_rem: BitSet = self.blackNeighborsCopied(rem);

            // if a black neighbor appears exactly at one node, it becomes red
            result.assignXor(&black_rem);

            // also already red neighbors stay red
            result.assignOr(self.constRedNeighbors(sur));
            result.assignOr(self.constRedNeighbors(rem));

            _ = result.unsetBit(rem);
            _ = result.unsetBit(sur);

            return result;
        }

        pub fn redDegreeInNeighborhoodAfterMerge(self: *const Self, rem: Node, sur: Node) Node {
            var reds = self.redNeighborsAfterMerge(rem, sur);
            var red_degree = reds.cardinality();
            reds.assignSub(self.constRedNeighbors(rem));
            reds.assignSub(self.constRedNeighbors(sur));

            // remaining reds were previously black, so this merge increases their degree
            if (reds.cardinality() > 0) {
                red_degree = @max(red_degree, self.maxRedDegreeIn(&reds) + 1);
            }

            return red_degree;
        }

        pub fn mergeNodes(self: *Self, rem: Node, sur: Node, red_neighbors: ?*const BitSet) void {
            if (red_neighbors) |reds| {
                // TODO: Implement this by bit magic
                var new_reds = reds.copyWithSub(self.redNeighbors(sur));
                _ = new_reds.unsetBit(sur);

                var iter = new_reds.iter_set();
                while (iter.next()) |v| {
                    _ = self.addEdge(v, sur, Color.Red);
                }

                self.removeEdgesAtNode(rem);

                self.assertIsConsistent();
            } else {
                const locally_computed_red = self.redNeighborsAfterMerge(rem, sur);
                return self.mergeNodes(rem, sur, &locally_computed_red);
            }
        }

        pub fn mergeDistantNodes(self: *Self, rem: Node, sur: Node) void {
            self.mergeNodes(rem, sur, null);
        }

        pub fn mergeDistantNodesOpt(self: *Self, rem: Node, sur: Node) void {
            var newReds = self.neighbors(sur).*;
            newReds.assignOr(self.constNeighbors(rem));
            newReds.assignSub(self.constRedNeighbors(sur));

            var iter = newReds.iter_set();
            while (iter.next()) |v| {
                _ = self.addEdge(sur, v, Color.Red);
            }

            _ = self.removeEdge(sur, sur);
            _ = self.removeEdgesAtNode(rem);

            self.assertIsConsistent();
        }

        pub fn neighbors(self: *Self, u: Node) *BitSet {
            return &self.matrix[3 * u + @enumToInt(RowType.Neighbor)];
        }

        pub fn constNeighbors(self: *const Self, u: Node) *const BitSet {
            return &self.matrix[3 * u + @enumToInt(RowType.Neighbor)];
        }

        pub fn redNeighbors(self: *Self, u: Node) *BitSet {
            return &self.matrix[3 * u + @enumToInt(RowType.Red)];
        }

        pub fn constRedNeighbors(self: *const Self, u: Node) *const BitSet {
            return &self.matrix[3 * u + @enumToInt(RowType.Red)];
        }

        pub fn twoNeighbors(self: *Self, u: Node) *BitSet {
            assert(!self.batch_updates);
            return &self.matrix[3 * u + @enumToInt(RowType.Two)];
        }

        pub fn constTwoNeighbors(self: *const Self, u: Node) *const BitSet {
            assert(!self.batch_updates);
            return &self.matrix[3 * u + @enumToInt(RowType.Two)];
        }

        fn recomputeTwoNeighbors(self: *Self, u: Node) void {
            var two_neighbors: *BitSet = self.twoNeighbors(u);
            two_neighbors.* = self.constNeighbors(u).*;

            var iter = self.constNeighbors(u).iter_set();
            while (iter.next()) |v| {
                two_neighbors.assignOr(self.constNeighbors(v));
            }
        }

        pub fn blackNeighborsCopied(self: *const Self, u: Node) BitSet {
            return self.constNeighbors(u).copyWithSub(self.constRedNeighbors(u));
        }

        pub fn maxRedDegree(self: *const Self) Node {
            var max_deg: Node = 0;
            var u: Node = 0;
            while (u < num_nodes) : (u += 1) {
                max_deg = @max(max_deg, self.redDeg(u));
            }
            return max_deg;
        }

        pub fn maxRedDegreeIn(self: *const Self, nodes: *const BitSet) Node {
            var max_deg: Node = 0;
            var iter = nodes.iter_set();
            while (iter.next()) |u| {
                max_deg = @max(max_deg, self.redDeg(u));
            }
            return max_deg;
        }

        pub fn assertIsConsistent(self: *const Self) void {
            if (builtin.mode != std.builtin.Mode.Debug) {
                return;
            }

            var num_edges: u32 = 0;
            var u: Node = 0;

            while (u < self.numberOfNodes()) : (u += 1) {
                var v: Node = u;

                assert(self.constRedNeighbors(u).is_subset_of(self.constNeighbors(u)));

                if (ExpensiveConsistencyChecks) {
                    while (v < self.numberOfNodes()) : (v += 1) {
                        assert(self.constNeighbors(u).isSet(v) == self.constNeighbors(v).isSet(u));
                        assert(self.constRedNeighbors(u).isSet(v) == self.constRedNeighbors(v).isSet(u));
                    }
                }

                num_edges += self.constNeighbors(u).cardinality();

                assert(self.has_neighbors.isSet(u) == (self.deg(u) > 0));

                if (!self.batch_updates) {
                    var two = self.constNeighbors(u).*;
                    var iter = self.constNeighbors(u).iter_set();
                    while (iter.next()) |i| {
                        two.assignOr(self.constNeighbors(i));
                    }
                    if (self.deg(u) > 0) {
                        _ = two.setBit(u);
                    }

                    var is_equal = self.constTwoNeighbors(u).is_equal(&two);

                    if (!is_equal) {
                        print("u: {d}\nExp: ", .{u});
                        two.debugPrint();
                        print("Act: ", .{});
                        self.constTwoNeighbors(u).debugPrint();

                        assert(false);
                    }
                }
            }
        }

        fn setColumn(self: *Self, comptime rowType: RowType, column: u32, iter: anytype) void {
            var myiter = iter;
            while (myiter.next()) |i| {
                _ = self.matrix[3 * i + @enumToInt(rowType)].setBit(column);
            }
        }

        fn unsetColumn(self: *Self, comptime rowType: RowType, column: u32, iter: anytype) void {
            var myiter = iter;
            while (myiter.next()) |i| {
                _ = self.matrix[3 * i + @enumToInt(rowType)].unsetBit(column);
            }
        }

        pub fn hash(self: *const Self) Digest {
            var digest = DigestAlgo.init(.{});
            digest.update(std.mem.sliceAsBytes(self.matrix[0..]));
            var output: Digest = undefined;
            digest.final(&output);
            return output;
        }

        pub fn numEdgesInComplement(self: *const Self) u32 {
            var totalNodes = self.has_neighbors.cardinality();
            if (totalNodes < 2) {
                return 0;
            }

            var redDegreeSum: u32 = 0;
            var u: u32 = 0;
            while (u < NumNodes) : (u += 1) {
                redDegreeSum = self.redDeg(u);
            }

            var blackEdges = self.numberOfEdges() - redDegreeSum / 2;
            var totalEdges = totalNodes * (totalNodes - 1) / 2;
            return totalEdges - blackEdges;
        }

        pub fn complement(self: *Self) void {
            self.beginBatchUpdates();

            var u: u32 = 0;
            var degreeSum: u32 = 0;
            while (u < NumNodes) : (u += 1) {
                if (!self.has_neighbors.isSet(u)) {
                    continue;
                }

                self.neighbors(u).assignXor(&self.has_neighbors);
                self.neighbors(u).assignOr(self.constRedNeighbors(u));

                var d = self.deg(u);

                if (d == 0) {
                    _ = self.has_neighbors.unsetBit(u);
                }

                degreeSum += d;
            }

            self.num_edges = degreeSum / 2;

            self.endBatchUpdates();
            self.assertIsConsistent();
        }
    };
}

test "New Matrix Graph" {
    var graph = MatrixGraph(100).new();
    assert(graph.numberOfNodes() == 100);
    assert(graph.numberOfEdges() == 0);
}

test "Add Edge" {
    var graph = MatrixGraph(8).new();
    assert(graph.numberOfEdges() == 0);
    assert(graph.addEdge(0, 1, Color.Black) == null);
    assert(graph.addEdge(1, 0, Color.Black) == Color.Black);

    assert(graph.neighbors(0).cardinality() == 1);
    assert(graph.neighbors(1).cardinality() == 1);
    assert(graph.redNeighbors(0).areAllUnset());
    assert(graph.redNeighbors(1).areAllUnset());

    assert(graph.addEdge(0, 1, Color.Red) == Color.Black);

    assert(graph.neighbors(0).cardinality() == 1);
    assert(graph.neighbors(1).cardinality() == 1);
    assert(graph.redNeighbors(0).cardinality() == 1);
    assert(graph.redNeighbors(1).cardinality() == 1);

    assert(graph.addEdge(0, 1, Color.Red) == Color.Red);
    assert(graph.numberOfEdges() == 1);

    assert(graph.addEdge(1, 1, Color.Black) == null);
    assert(graph.numberOfEdges() == 2);
}

test "Remove Edge" {
    var graph = MatrixGraph(100).new();

    assert(graph.numberOfEdges() == 0);
    assert(graph.addEdge(0, 1, Color.Black) == null);
    assert(graph.removeEdge(0, 3) == null);
    assert(graph.numberOfEdges() == 1);
    assert(graph.removeEdge(1, 0) == Color.Black);
    assert(graph.numberOfEdges() == 0);

    assert(graph.removeEdge(0, 2) == null);
    assert(graph.numberOfEdges() == 0);
}

test "Merge" {
    var graph = MatrixGraph(10).new();

    _ = graph.addEdge(0, 1, Color.Black);
    _ = graph.addEdge(0, 2, Color.Black);
    _ = graph.addEdge(0, 3, Color.Black);
    _ = graph.addEdge(0, 4, Color.Red);

    _ = graph.addEdge(9, 1, Color.Black);
    _ = graph.addEdge(9, 2, Color.Red);
    _ = graph.addEdge(9, 4, Color.Red);
    _ = graph.addEdge(9, 5, Color.Red);
    _ = graph.addEdge(9, 6, Color.Black);

    graph.mergeNodes(0, 9, null);

    assert(graph.deg(0) == 0);
    assert(graph.deg(9) == 6);
    assert(graph.redDeg(9) == 5);
}
