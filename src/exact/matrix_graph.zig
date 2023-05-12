const std = @import("std");
const print = std.debug.print;
const FixedSizeSet = @import("fixed_size_bitset.zig").FixedSizeSet;
const builtin = @import("builtin");
const assert = std.debug.assert;

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

pub fn MatrixGraph(comptime num_nodes: u32) type {
    return struct {
        const Self = @This();
        pub const numNodes: u32 = num_nodes;
        pub const Node = u32;
        pub const BitSet = FixedSizeSet(num_nodes);

        num_edges: u32,
        matrix: [2 * num_nodes]BitSet,
        has_neighbors: BitSet,

        pub fn new() Self {
            var graph = Self{ .num_edges = 0, .matrix = undefined, .has_neighbors = BitSet.new() };

            const empty = FixedSizeSet(num_nodes).new();
            for (&graph.matrix) |*r| r.* = empty;

            return graph;
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

        pub fn addEdge(self: *Self, u: Node, v: Node, color: Color) ?Color {
            var prev: ?Color = null;

            if (self.neighbors(u).setBit(v)) {
                prev = Color.Black;
            } else {
                _ = self.neighbors(v).setBit(u);

                self.num_edges += 1;
                _ = self.has_neighbors.setBit(u);
                _ = self.has_neighbors.setBit(v);
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

        pub fn neighbors(self: *Self, u: Node) *BitSet {
            return &self.matrix[2 * u];
        }

        pub fn constNeighbors(self: *const Self, u: Node) *const BitSet {
            return &self.matrix[2 * u];
        }

        pub fn redNeighbors(self: *Self, u: Node) *BitSet {
            return &self.matrix[2 * u + 1];
        }

        pub fn constRedNeighbors(self: *const Self, u: Node) *const BitSet {
            return &self.matrix[2 * u + 1];
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

        pub fn assertIsConsistent(self: *Self) void {
            if (builtin.mode != std.builtin.Mode.Debug) {
                return;
            }

            var num_edges: u32 = 0;
            var u: Node = 0;

            while (u < self.numberOfNodes()) : (u += 1) {
                var v: Node = u;

                assert(self.constRedNeighbors(u).is_subset_of(self.constNeighbors(u)));
                while (v < self.numberOfNodes()) : (v += 1) {
                    assert(self.constNeighbors(u).isSet(v) == self.constNeighbors(v).isSet(u));
                    assert(self.constRedNeighbors(u).isSet(v) == self.constRedNeighbors(v).isSet(u));
                }

                num_edges += self.constNeighbors(u).cardinality();

                assert(self.has_neighbors.isSet(u) == (self.deg(u) > 0));
            }
        }
    };
}

test "New Matrix Graph" {
    var graph = MatrixGraph(100).new();
    assert(graph.numberOfNodes() == 100);
    assert(graph.numberOfEdges() == 0);
}

test "Add Edge" {
    var graph = MatrixGraph(100).new();
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
