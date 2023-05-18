const std = @import("std");
const Graph = @import("graph.zig").Graph;
const edge_list = @import("edge_list.zig");
const bitset = @import("../util/two_level_bitset.zig");
const bfs_mod = @import("bfs.zig");
const topk = @import("../util/top_k_scorer.zig");
const contraction = @import("../tww/contraction_sequence.zig");
const node_priority_queue = @import("node_priority_queue.zig");
const solver_resources = @import("solver.zig");
const retraceable_contraction = @import("../tww/retraceable_contraction_sequence.zig");
const dfs_mod = @import("dfs.zig");
const benchmark_helper = @import("../util/benchmark_helper.zig");
const min_hash = @import("../util/min_hash.zig");
const NodePairs = @import("solver.zig").NodePairs;
const NodeTuple = @import("solver.zig").NodeTuple;

pub const SubGraphError = error{
    SolverResourcesNotInitialized,
    SolverResourcesAlreadyInitialized,
};

pub fn InducedSubGraph(comptime T: type) type {
    return struct {
        const Self = @This();
        graph: *Graph(T),
        // Must be this type to support fast lookups!
        nodes: []T,

        pub fn fromSlice(graph: *Graph(T), nodes: []T) Self {
            return Self{
                .nodes = nodes,
                .graph = graph,
            };
        }

        pub const TargetMinimalInducedTww = struct {
            potential: Graph(T).InducedTwinWidthPotential,
            target: T,
        };

        pub fn compareNodeDegrees(ctx: *const Graph(T), lhs: T, rhs: T) std.math.Order {
            return std.math.order(ctx.node_list[rhs].cardinality(), ctx.node_list[lhs].cardinality());
        }

        pub inline fn selectBestMoveOfIter(self: *Self, comptime Iter: type, iter: *Iter, first_node: T, current_tww: T) TargetMinimalInducedTww {
            var induced_tww = Graph(T).InducedTwinWidthPotential.default();
            var min_target: T = 0;
            while (iter.next()) |item| {
                if (item == first_node) continue;
                if (self.graph.erased_nodes.get(item)) continue;

                const cal_induced_tww = self.graph.calculateInducedTwwPotential(item, first_node, &induced_tww, current_tww);

                if (cal_induced_tww.isLess(induced_tww, current_tww)) {
                    induced_tww = cal_induced_tww;
                    min_target = item;
                }
            }
            return TargetMinimalInducedTww{ .potential = induced_tww, .target = min_target };
        }

        pub inline fn selectBestMoveIncremental(self: *Self, comptime K: u32, comptime P: u32, first_node: T, solver: *solver_resources.SolverResources(T, K, P), current_tww: T, first_level_merge: bool, erased: T, red_edges_erased: []T) !TargetMinimalInducedTww {
            var induced_tww = Graph(T).InducedTwinWidthPotential.default();

            var min_target: T = 0;
            const result = bfs_mod.bfs_topk_high_performance_incremental(T, K, first_node, erased, first_level_merge, red_edges_erased, self.graph, &solver.scratch_bitset, &solver.scorer);
            if (!result) return TargetMinimalInducedTww{ .potential = induced_tww, .target = min_target };

            var iterator = try solver.scorer.iterator(self.graph);

            const result_best_move = self.selectBestMoveOfIter(@TypeOf(iterator), &iterator, first_node, current_tww);

            if (result_best_move.potential.cumulative_red_edges == std.math.maxInt(i64)) {
                @panic("Error cannot find contraction partner!");
            }

            return result_best_move;
        }

        pub inline fn selectBestMove(self: *Self, comptime K: u32, comptime P: u32, first_node: T, solver: *solver_resources.SolverResources(T, K, P), current_tww: T) !TargetMinimalInducedTww {
            bfs_mod.bfs_topk_high_performance(T, K, first_node, self.graph, &solver.scratch_bitset, &solver.scorer);

            var iterator = try solver.scorer.iterator(self.graph);

            const result = self.selectBestMoveOfIter(@TypeOf(iterator), &iterator, first_node, current_tww);

            if (result.potential.cumulative_red_edges == std.math.maxInt(i64)) {
                @panic("Error cannot find contraction partner!");
            }

            return result;
        }

        pub inline fn selectBestMoveExhaustive(self: *Self, comptime N: u32, remaining_nodes: *std.BoundedArray(T, N), current_tww: T) contraction.Contraction(T) {
            var induced_tww = Graph(T).InducedTwinWidthPotential.default();
            var min_contraction = contraction.Contraction(T){ .erased = 0, .survivor = 0 };
            var min_index_erased: T = 0;

            for (0..remaining_nodes.len) |first| {
                for (0..first) |second| {
                    const first_node = remaining_nodes.buffer[first];
                    const second_node = remaining_nodes.buffer[second];
                    var ind = self.graph.calculateInducedTwwPotential(first_node, second_node, &induced_tww, current_tww);
                    if (ind.isLess(induced_tww, current_tww)) {
                        induced_tww = ind;
                        min_contraction = .{ .erased = first_node, .survivor = second_node };
                        min_index_erased = @intCast(T, first);
                    }
                }
            }

            _ = remaining_nodes.swapRemove(min_index_erased);

            return min_contraction;
        }

        pub fn compareTValue(ctx: void, lhs: T, rhs: T) std.math.Order {
            _ = ctx;
            return std.math.order(rhs, lhs);
        }

        pub const LargestArticulationPoint = struct {
            point: T,
            second_largest_component_size: T,
        };

        pub fn findLargestArticulationPointCut(self: *Self, comptime K: u32, comptime P: u32, solver: *solver_resources.SolverResources(T, K, P)) !?T {
            solver.scratch_bitset.unsetAll();
            try dfs_mod.findArticulationPoints(T, self.nodes[0], self.graph, &solver.scratch_bitset, &solver.articulation_point_resources);
            solver.scratch_bitset.unsetAll();

            // Reset bitset we need it to use it as a scratch for the visited field
            if (solver.articulation_point_resources.articulation_points.cardinality == 0) {
                return null;
            }

            var priority_queue_components = std.PriorityQueue(T, void, Self.compareTValue).init(self.graph.allocator, {});
            defer priority_queue_components.deinit();

            try priority_queue_components.ensureTotalCapacity(self.graph.number_of_nodes);

            var largest = LargestArticulationPoint{ .point = 0, .second_largest_component_size = 0 };

            var iter_art = solver.articulation_point_resources.articulation_points.iter();
            while (iter_art.next()) |item| {
                var nb_iter = self.graph.node_list[item].unorderedIterator();
                var total_count: T = 0;

                while (nb_iter.next()) |nb| {
                    // Ok new connected component
                    // Check root node
                    if (solver.articulation_point_resources.parent_array[nb] == item) {
                        if (item == self.nodes[0]) {
                            // Root nodes will cover everything with that
                            try priority_queue_components.add(solver.articulation_point_resources.node_count[nb]);
                        } else if (solver.articulation_point_resources.low_array[nb] >= solver.articulation_point_resources.depth_array[item]) {
                            try priority_queue_components.add(solver.articulation_point_resources.node_count[nb]);
                            total_count += solver.articulation_point_resources.node_count[nb];
                        }
                    }
                }

                if (total_count != 0) {
                    // Add the left over connected component
                    try priority_queue_components.add(@intCast(T, self.nodes.len) - (total_count + 1));
                }

                _ = priority_queue_components.remove();
                const second_largest = priority_queue_components.remove();
                if (second_largest > largest.second_largest_component_size) {
                    largest.point = @intCast(T, item);
                    largest.second_largest_component_size = second_largest;
                }

                // Delete all items
                priority_queue_components.len = 0;
            }

            return largest.point;
        }

        pub inline fn addContractionAndLeafReduction(self: *Self, enable_follow_up_merge: *bool, seq: *retraceable_contraction.RetraceableContractionSequence(T), min_contraction: contraction.Contraction(T), contractions_left: *u32) !T {
            var total_tww: T = 0;
            enable_follow_up_merge.* = true;

            // Adjust the contractions left and maximal number of postpones allowed (Note this is only problematic in the case of contractions_left < P)
            contractions_left.* -= 1;

            var last_evoked_twin_width = try self.graph.addContraction(min_contraction.erased, min_contraction.survivor, seq);
            total_tww = std.math.max(last_evoked_twin_width, total_tww);

            // Reduce leafes if one was created
            if (self.graph.node_list[min_contraction.survivor].isLeaf()) {
                // Ok is leaf check for other leaf
                var parent = self.graph.node_list[min_contraction.survivor].getFirstNeighboor();
                if (self.graph.node_list[parent].num_leafes > 1) {
                    var parent_node_iter = self.graph.node_list[parent].unorderedIterator();
                    while (parent_node_iter.next()) |item| {
                        if (item == min_contraction.survivor) continue;
                        if (self.graph.node_list[item].isLeaf()) {
                            last_evoked_twin_width = try self.graph.addContraction(item, min_contraction.survivor, seq);
                            total_tww = std.math.max(last_evoked_twin_width, total_tww);
                            contractions_left.* -= 1;
                            enable_follow_up_merge.* = false;
                            // Only merge one other leaf since we do this for every leaf and there should not be any more leafes left
                            // WARNING: The iterator probably break if we try to go on since addContraction will update the edge list
                            break;
                        }
                    }
                }
            }
            // Reduce leaves if we are the parent of at least 2 leafes
            else if (self.graph.node_list[min_contraction.survivor].num_leafes > 1) {
                var parent_node_iter = self.graph.node_list[min_contraction.survivor].unorderedIterator();
                var first_leaf: ?T = null;

                while (parent_node_iter.next()) |item| {
                    if (self.graph.node_list[item].isLeaf()) {
                        if (first_leaf) |k| {
                            last_evoked_twin_width = try self.graph.addContraction(k, item, seq);
                            total_tww = std.math.max(last_evoked_twin_width, total_tww);
                            contractions_left.* -= 1;
                            enable_follow_up_merge.* = false;
                            break;
                        } else {
                            first_leaf = item;
                        }
                    }
                }
            }

            return total_tww;
        }

        pub fn LookaheadResult(comptime NodeCtx: type) type {
            return struct {
                const Own = @This();
                result: ?NodeCtx.ScoreType,
                failed_score: ?NodeCtx.ScoreType,
                depth: u32,

                pub inline fn default() Own {
                    return Own{
                        .result = null,
                        .failed_score = null,
                        .depth = 0,
                    };
                }

                pub inline fn assignIfBetter(self: *Own, other: Own, ctx: *NodeCtx, generator: *std.rand.DefaultPrng) bool {
                    if (self.result) |r_own| {
                        if (other.result) |r_other| {
                            switch (ctx.compare(r_other, r_own)) {
                                .gt => {
                                    self.* = other;
                                    return true;
                                },
                                .eq => {
                                    // On equal score chose at random if we overwrite or leave it
                                    if (generator.next() % 2 == 0) {
                                        self.* = other;
                                        return true;
                                    } else {
                                        return false;
                                    }
                                },
                                .lt => {
                                    return false;
                                },
                            }
                        } else {
                            return false;
                        }
                    } else {
                        // Asign if depth of other is greater
                        if (other.result != null or self.depth <= other.depth) {
                            if (self.depth == other.depth and other.result == null) {
                                if (self.failed_score != null and other.failed_score != null) {
                                    switch (ctx.compare(other.failed_score.?, self.failed_score.?)) {
                                        .gt => {
                                            self.* = other;
                                            return true;
                                        },
                                        .eq => {
                                            // On equal score chose at random if we overwrite or leave it
                                            if (generator.next() % 2 == 0) {
                                                self.* = other;
                                                return true;
                                            } else {
                                                return false;
                                            }
                                        },
                                        .lt => {
                                            return false;
                                        },
                                    }
                                }
                            }
                            self.* = other;
                            return true;
                        }
                        return false;
                    }
                }
            };
        }

				pub fn findOneMoveBelowUpperBound(self: *Self, bfs_queue: *bfs_mod.BfsQueue(T), visited: *bitset.FastBitSet, upper_bound: T, seed: u64, less: std.PriorityQueue()) !?contraction.Contraction(T) {
				    _ = less;
						var bounded_arr = try std.BoundedArray(contraction.Contraction(T),100).init(0);


				    
						for (self.nodes) |node| {
							if (self.graph.erased_nodes.get(node)) continue;
							visited.unsetAll();
							var bf = bfs_mod.bfs(T, node, self.graph, visited, bfs_queue, .{ .max_level = 2, .kind = .both });
							while (bf.next()) |other_node| {
								if (other_node >= node) continue;

								const result = self.graph.calculateInducedTww(other_node,node,null);

								if(result.tww >= upper_bound) continue;
								if(bounded_arr.len < 100) {
									try bounded_arr.append(contraction.Contraction(T) {
											.erased = other_node,
											.survivor = node,
									});
								}
								else {
									break;
								}
							}
						}

						if(bounded_arr.len == 0) return null;

						var gen = std.rand.DefaultPrng.init(seed);
						return bounded_arr.buffer[gen.next() % bounded_arr.len];
				}

				pub fn solveGreedyBacktrackToLastIncrease(self: *Self, best_contraction_seq: *contraction.ContractionSequence(T),bfs_queue: *bfs_mod.BfsQueue(T), visited: *bitset.FastBitSet, seq: *retraceable_contraction.RetraceableContractionSequence(T), upper_bound: T) !T {
					var check = best_contraction_seq.iterator();
					while(check.next()) |item| {
						_ = try self.graph.addContraction(item.erased,item.survivor,seq);
					}
					
					// Ok we have the best known current sequence
					var current_tww: T = seq.getTwinWidth();

					var current_level = self.nodes.len;
					while(seq.lastContraction()) |_| {
						try self.graph.revertLastContraction(seq);
						current_level-=1;
						if(seq.getTwinWidth() != current_tww) {
							break;
						}
					}

					var seed:u64 = 17;

					var max_backtracks:u32 = 10000;
					var generator = std.rand.DefaultPrng.init(seed);

					// Here is somewehere another error
					while(!seq.isComplete()) {
						seed+=17;
						var move = try self.findOneMoveBelowUpperBound(bfs_queue,visited,upper_bound, seed);
						if(move == null) {
							max_backtracks-=1;
							if(max_backtracks == 0) return std.math.maxInt(T);

							while(seq.lastContraction()) |_| {
								if(generator.next()%4 != 0) {
									try self.graph.revertLastContraction(seq);
								}
								else {
									break;
								}
							}
							continue;
						}

						_ = try self.graph.addContraction(move.?.erased, move.?.survivor, seq);
						current_level+=1;
					}

					return seq.getTwinWidth();
				}

        pub fn solveGreedyConstraint(self: *Self, bfs_queue: *bfs_mod.BfsQueue(T), visited: *bitset.FastBitSet, seq: *retraceable_contraction.RetraceableContractionSequence(T), upper_bound: T) !T {
            var constraints = std.AutoHashMap(T, T).init(self.graph.allocator);

            var contractions_left: u32 = @intCast(u32, self.nodes.len - 1);

            var total_max_reverts:u32 = 1_000;

            while (contractions_left > 0) {
                var total_nodes_checked: u32 = 0;
                var min_tww: ?T = null;
                var min_contraction = contraction.Contraction(T){
                    .erased = 0,
                    .survivor = 0,
                };

                for (self.nodes) |node| {
                    if (self.graph.erased_nodes.get(node)) continue;
                    visited.unsetAll();
                    var bf = bfs_mod.bfs(T, node, self.graph, visited, bfs_queue, .{ .max_level = 2, .kind = .both });
                    outer: while (bf.next()) |other_node| {
                        if (other_node >= node) continue;

                        const result = try self.graph.addContraction(other_node, node, seq);
                        for (self.nodes) |current_node| {
                            if (self.graph.erased_nodes.get(current_node)) continue;
                            if (constraints.get(current_node)) |less_than| {
                                if (self.graph.node_list[current_node].red_edges.cardinality() >= less_than) {
                                    try self.graph.revertLastContraction(seq);
                                    continue :outer;
                                }
                            }
                        }
                        const compareable_result = self.graph.node_list[node].red_edges.cardinality();
                        try self.graph.revertLastContraction(seq);
                        if (result >= upper_bound) continue;

                        if (min_tww == null) {
                            min_tww = compareable_result;
                            min_contraction = contraction.Contraction(T){
                                .erased = other_node,
                                .survivor = node,
                            };
                        } else if (min_tww.? >= compareable_result) {
                            min_tww = compareable_result;
                            min_contraction = contraction.Contraction(T){
                                .erased = other_node,
                                .survivor = node,
                            };
                        }

                        total_nodes_checked += 1;
                    }
                }
                if (total_nodes_checked == 0) {
                    total_max_reverts -= 1;
                    if (total_max_reverts == 0) return std.math.maxInt(T);
										constraints.clearRetainingCapacity();

                    for (self.nodes) |current_node| {
                        if (self.graph.erased_nodes.get(current_node)) continue;
												if(self.graph.node_list[current_node].red_edges.cardinality() == upper_bound-1) {
													try constraints.put(current_node,upper_bound-1);
													break;
												}
                    }

                    while (seq.lastContraction()) |_| {
                        try self.graph.revertLastContraction(seq);
                    }
            				contractions_left = @intCast(u32, self.nodes.len - 1);
										continue;
                }

                _ = try self.graph.addContraction(min_contraction.erased, min_contraction.survivor, seq);

                if (seq.getTwinWidth() >= upper_bound) {
                    return std.math.maxInt(T);
                }
                contractions_left -= 1;
            }

            return seq.getTwinWidth();
        }

        fn solveGreedyLookaheadGeneric(self: *Self, comptime NodeCtx: type, comptime ScoreCtx: type, ctx: *NodeCtx, ctx_score: *ScoreCtx, level: u32, bfs_queue: *bfs_mod.BfsQueue(T), visited: *bitset.FastBitSet, seq: *retraceable_contraction.RetraceableContractionSequence(T), seed: u64) !LookaheadResult(ScoreCtx) {
            const initial_level = level;
            var modifiable_level = level;

            const max_trials_per_level: u32 = 1;
            var priority_queue = std.PriorityQueue(NodeCtx.MapType, *NodeCtx, NodeCtx.compare).init(self.graph.allocator, ctx);
            try priority_queue.ensureTotalCapacity(max_trials_per_level);
            defer priority_queue.deinit();

            defer {
                for (0..(initial_level - modifiable_level)) |_| {
                    self.graph.revertLastContraction(seq) catch {
                        @panic("Error in defer revert last contraction!");
                    };
                }
            }

            var generator = std.rand.DefaultPrng.init(seq.lastContraction().?.erased + seed);

            while (modifiable_level > 0) {
                priority_queue.len = 0;

                var nodes_left: u32 = 0;
                for (self.nodes) |node| {
                    if (self.graph.erased_nodes.get(node)) continue;
                    visited.unsetAll();
                    nodes_left += 1;
                    var bf = bfs_mod.bfs(T, node, self.graph, visited, bfs_queue, .{ .max_level = 2, .kind = .both });
                    while (bf.next()) |other_node| {
                        if (other_node >= node) continue;
                        if (self.graph.erased_nodes.get(other_node)) {
                            @panic("A BFS should never return an erased node!");
                        }

                        const mapped = ctx.map(other_node, node, self.graph);
                        if (mapped == null) continue;
                        if (priority_queue.count() < max_trials_per_level) {
                            try priority_queue.add(mapped.?);
                        } else if (priority_queue.peek()) |top_prio| {
                            switch (ctx.compare(mapped.?, top_prio)) {
                                .gt => {
                                    _ = priority_queue.remove();
                                    try priority_queue.add(mapped.?);
                                },
                                .eq => {
                                    // On equal score chose at random if we overwrite or leave it
                                    if (generator.next() % 2 == 0) {
                                        _ = priority_queue.remove();
                                        try priority_queue.add(mapped.?);
                                    }
                                },
                                .lt => {},
                            }
                        }
                    }
                }
                if (nodes_left <= 1) {
                    return LookaheadResult(ScoreCtx){
                        .result = ctx_score.evaluate(self.graph, seq),
                        .failed_score = null,
                        .depth = level - modifiable_level,
                    };
                }
                if (priority_queue.count() == 0) {
                    return LookaheadResult(ScoreCtx){
                        .result = null,
                        .failed_score = ctx_score.evaluate(self.graph, seq),
                        .depth = level - modifiable_level,
                    };
                }

                var selected = generator.next() % priority_queue.count();
                const item = priority_queue.items[selected];
                _ = try self.graph.addContraction(item.erased, item.survivor, seq);
                modifiable_level -= 1;
            }

            return LookaheadResult(ScoreCtx){
                .result = ctx_score.evaluate(self.graph, seq),
                .failed_score = null,
                .depth = level,
            };
        }

        fn collectNodePairsIntoListGeneric(self: *Self, comptime Ctx: type, ctx: *Ctx, k: u32, bfs_queue: *bfs_mod.BfsQueue(T), visited: *bitset.FastBitSet, node_priority: *std.PriorityQueue(Ctx.MapType, *Ctx, Ctx.compare)) !void {
            for (self.nodes) |node| {
                if (self.graph.erased_nodes.get(node)) continue;
                visited.unsetAll();
                var bf = bfs_mod.bfs(T, node, self.graph, visited, bfs_queue, .{ .max_level = 2, .kind = .both });
                while (bf.next()) |other_node| {
                    if (other_node >= node) continue;

                    const result = ctx.map(other_node, node, self.graph);
                    if (result) |r| {
                        if (node_priority.count() < k) {
                            try node_priority.add(r);
                        } else if (node_priority.peek()) |top_prio| {
                            if (ctx.compare(r, top_prio) == .gt) {
                                _ = node_priority.remove();
                                try node_priority.add(r);
                            }
                        }
                    }
                }
            }
        }

        pub fn solveGreedyLookahead(self: *Self, comptime NodeCtx: type, comptime ScoreCtx: type, ctx: *NodeCtx, ctx_score: *ScoreCtx, k: u32, lookahead: u32, seq: *retraceable_contraction.RetraceableContractionSequence(T), bfs_queue: *bfs_mod.BfsQueue(T), visited: *bitset.FastBitSet, upper_bound: ?T) !T {
            var priority_queue = std.PriorityQueue(NodeCtx.MapType, *NodeCtx, NodeCtx.compare).init(self.graph.allocator, ctx);
            defer priority_queue.deinit();

            var contractions_left: u32 = @intCast(u32, self.nodes.len - 1);
            // Hard graphs exact 6,14,18,32,38

            var upper_bound_set = upper_bound orelse std.math.maxInt(T);
						if(upper_bound_set == 0) return std.math.maxInt(T);
            ctx.setUpperBound(upper_bound_set - 1);
            ctx_score.setUpperBound(upper_bound_set - 1);

            var seed: u64 = 0;

            var generator = std.rand.DefaultPrng.init(17);
            while (contractions_left > 0) {
                priority_queue.len = 0;
                try self.collectNodePairsIntoListGeneric(NodeCtx, ctx, k, bfs_queue, visited, &priority_queue);

                var max_result: LookaheadResult(ScoreCtx) = LookaheadResult(ScoreCtx).default();
                var max_contraction: contraction.Contraction(T) = undefined;
                var total_egligible_moves: u32 = 0;

                while (priority_queue.removeOrNull()) |i| {
                    total_egligible_moves += 1;
                    _ = try self.graph.addContraction(i.erased, i.survivor, seq);
                    seed += 19;
                    const result = try self.solveGreedyLookaheadGeneric(NodeCtx, ScoreCtx, ctx, ctx_score, lookahead, bfs_queue, visited, seq, seed);

                    try self.graph.revertLastContraction(seq);

                    if (max_result.assignIfBetter(result, ctx_score, &generator)) {
                        max_contraction = .{ .survivor = i.survivor, .erased = i.erased };
                    }
                }
                if (total_egligible_moves == 0) {
                    // Could not find solution
                    return std.math.maxInt(T);
                }
                _ = try self.graph.addContraction(max_contraction.erased, max_contraction.survivor, seq);
                contractions_left -= 1;
            }
            return seq.getTwinWidth();
        }

				pub const TwinWidthPriorityWithContraction = struct {
					tww: T,
					pub fn compare(ctx: void, lhs: TwinWidthPriorityWithContraction, rhs: TwinWidthPriorityWithContraction) std.math.Order {
					    _ = ctx;
						return std.math.order(rhs.tww,lhs.tww);
					}
				};

				pub fn solveFoldingSingleNode(self: *Self, comptime K: u32, comptime P: u32, seq: *retraceable_contraction.RetraceableContractionSequence(T), solver: *solver_resources.SolverResources(T, K, P)) !T {
				    
            var contractions_left: u32 = @intCast(u32, self.nodes.len - 1);
            if (contractions_left == 0) return 0;

						var generator = std.rand.DefaultPrng.init(23);
						min_hash.fisher_yates_shuffle(T,self.nodes[0..self.nodes.len], &generator);
						
						var min_node:T = 0;
						var min_tww:T = std.math.maxInt(T);


						for(0..std.math.min(self.nodes.len,1000)) |i| {
							const item = self.nodes[i];

              solver.scorer.unsetVisitedBitset(&solver.scratch_bitset);

							// Select greedy
              var selection = try self.selectBestMove(K, P, item, solver, 0);
							if(selection.potential.tww < min_tww) {
								min_node = item;
								min_tww = selection.potential.tww;
							}
						}

						while(contractions_left > 0) {
              solver.scorer.unsetVisitedBitset(&solver.scratch_bitset);

							// Select greedy
              var selection = try self.selectBestMove(K, P, min_node, solver, 0);
							_ = try self.graph.addContraction(selection.target, min_node, seq);
							contractions_left-=1;
							std.debug.print("Folding contractions left {} tww {}\n",.{contractions_left, seq.getTwinWidth()});
						}

						return seq.getTwinWidth();

				}


				pub fn solveMinHash(self: *Self, comptime K: u32, comptime P: u32, seq: *retraceable_contraction.RetraceableContractionSequence(T), solver: *solver_resources.SolverResources(T, K, P)) !T {
            solver.scratch_bitset.unsetAll();
            var contractions_left: u32 = @intCast(u32, self.nodes.len - 1);
						var diff_nodes = std.AutoHashMap(u32,u32).init(self.graph.allocator);

            if (contractions_left == 0) return 0;
            var min_contraction = contraction.Contraction(T){ .erased = 0, .survivor = 0 };
            while (contractions_left > 0) {
								//try self.graph.min_hash.fetchNextMoves(120,self.graph);

								var potential = Graph(T).InducedTwinWidthPotential.default();

								var min_mv:u32 = std.math.maxInt(u32);
								var max_mv:u32 = 0;
								diff_nodes.clearRetainingCapacity();
								var def = Graph(T).InducedTwinWidthPotential.default();
								for(0..self.graph.min_hash.fetched_moves.items.len) |mv_index| {
									const mv = self.graph.min_hash.fetched_moves.items[mv_index];

									const move = mv.intoMove(T,self.graph.number_of_nodes);
									var global_score = self.graph.calculateInducedTwwPotential(move.erased,move.survivor,&def, seq.getTwinWidth());
									if(diff_nodes.getPtr(move.erased)) |er| {
										er.* += 1;
										max_mv = std.math.max(max_mv,er.*);
									}
									else {
										try diff_nodes.put(move.erased,1);
									}
									if(diff_nodes.getPtr(move.survivor)) |sr| {
										sr.* += 1;
										max_mv = std.math.max(max_mv,sr.*);
									}
									else {
										try diff_nodes.put(move.survivor,1);
									}
									min_mv = std.math.min(min_mv,global_score.tww);
									if(global_score.tww < potential.tww) {
										min_contraction = move;
										potential = global_score;
									}
								}

								std.debug.print("Min tww {} and selected {} total number of different nodes {} max_mv {}\n",.{min_mv,potential.tww, diff_nodes.count(), max_mv});
								const length = self.graph.min_hash.fetched_moves.items.len;
								try self.graph.min_hash.reinsertFetchedMoves();
								_ = self.graph.addContraction(min_contraction.erased,min_contraction.survivor,seq) catch |err| {
									std.debug.print("Contractions left {} length {} tww {}\n",.{contractions_left, length, seq.getTwinWidth()});
									return err;
								};
								contractions_left-=1;
						}
						return seq.getTwinWidth();
				}

				pub fn solveSweepingSolverTopK(self: *Self, comptime K: u32, comptime P: u32, seq: *retraceable_contraction.RetraceableContractionSequence(T), solver: *solver_resources.SolverResources(T, K, P), budget_secs: u64, probing: bool) !T {
            solver.scratch_bitset.unsetAll();
            var contractions_left: u32 = @intCast(u32, self.nodes.len - 1);
            if (contractions_left == 0) return 0;
            var min_contraction = contraction.Contraction(T){ .erased = 0, .survivor = 0 };
            _ = min_contraction;
						// Should be zero which leads in the best move selection to always greedily choose the better twin width rather than the better potential
            var total_tww: T = seq.getTwinWidth();
						
						var start_time = try std.time.Instant.now();
						

						var priority_queue = std.PriorityQueue(TwinWidthPriorityWithContraction,void,TwinWidthPriorityWithContraction.compare).init(self.graph.allocator,{});
						defer priority_queue.deinit();

						// This is an overestimation can be cut down to N/(2*log(N))
						try priority_queue.ensureTotalCapacity(self.nodes.len);

						//try self.reduceLeafesAndPaths(K,P,seq,solver,false);
						var remaining_nodes:T = 0;
						for(self.nodes) |n| {
							if(self.graph.erased_nodes.get(n)) @panic("Should never happen!");
							solver.scratch_node_list[remaining_nodes] = n;
							remaining_nodes+=1;
						}
						
						var sweeping_thresh:u32 = 0;

						var generator = std.rand.DefaultPrng.init(19);
						//heuristic_128.gr

						
						var first_iteration: bool = true;
            outer: while (contractions_left > 0) {
							priority_queue.len = 0;
							var total_contractions:u32 = 0;
							solver.node_mask_bitset.unsetAll();


							const remaining_nodes_before = remaining_nodes;

							var sqrt = remaining_nodes_before/std.math.log2(remaining_nodes_before);

							var sample_amount = 2*sqrt;
							if(sample_amount > 200) {
								min_hash.fisher_yates_sample_first_n(T, solver.scratch_node_list[0..remaining_nodes], sample_amount, &generator);
							}
							else {
								sample_amount = remaining_nodes+1;
							}

							var min_tww: T = std.math.maxInt(T);
							var max_tww:T = 0;
							var cumulative_tww:u64 = 0;
							var visited:u32 = 0;

							// Need overflow
							var index:i32 = 0;
							while (index < remaining_nodes): (index+=1) {
								const item = solver.scratch_node_list[@intCast(u32,index)];
								if(solver.node_mask_bitset.get(item) and self.nodes.len > 10000) continue;
								if(self.graph.erased_nodes.get(item)) @panic("Should never happen!");
                solver.scorer.unsetVisitedBitset(&solver.scratch_bitset);
                var selection = try self.selectBestMove(K, P, item, solver, total_tww);
									
								if(selection.potential.tww <= sweeping_thresh) {
									_ = try self.graph.addContraction(item,selection.target,seq);
									contractions_left-=1;
									if(contractions_left==0) break :outer;

									remaining_nodes-=1;
									solver.scratch_node_list[@intCast(u32,index)] = solver.scratch_node_list[remaining_nodes];
									// Swap remove node
									total_contractions += 1;
								}
								solver.node_mask_bitset.set(item);
								solver.node_mask_bitset.set(selection.target);
								min_tww = std.math.min(min_tww,selection.potential.tww);
								max_tww = std.math.max(max_tww,selection.potential.tww);
								cumulative_tww += selection.potential.tww;
								visited+=1;
								if(priority_queue.count() < sqrt) {
									try priority_queue.add(TwinWidthPriorityWithContraction {
											.tww = selection.potential.tww,
											});
								}
								else if(priority_queue.peek().?.tww > selection.potential.tww) {
									try priority_queue.update(priority_queue.peek().?, TwinWidthPriorityWithContraction {
											.tww = selection.potential.tww,
											});
								}

								if(index == sample_amount) {
									sweeping_thresh = std.math.max(sweeping_thresh,min_tww);
									if(first_iteration) {
										if(probing and (((cumulative_tww/visited) > 50) or max_tww > 50)) return std.math.maxInt(T);
									}
								}

								if(index&0x800 > 0) {
									var time = try std.time.Instant.now();
									if(time.since(start_time)/(1000*1000*1000) > budget_secs) {
										return std.math.maxInt(T);	
									}
								}
							}
							first_iteration = false;

							
							if(total_contractions < sqrt) {
								var new_thresh:u32 = priority_queue.peek().?.tww;
								// Reset it since it may be bad due to the sampling
								sweeping_thresh = new_thresh;
							}
						}
						return seq.getTwinWidth();
				}

        pub fn solveGreedyTopK(self: *Self, comptime K: u32, comptime P: u32, seq: *retraceable_contraction.RetraceableContractionSequence(T), solver: *solver_resources.SolverResources(T, K, P), find_articulation_points: bool) !T {
            _ = find_articulation_points;
            // Use once at the beginning later we will only consider merges which involve the survivor of the last merge
            // Paths might be dangerous
            try self.reduceLeafesAndPaths(K, P, seq, solver, true);

            // Reset bitset we need it to use it as a scratch for the visited field
            solver.scratch_bitset.unsetAll();

            // Sadly only small components can be split of while removing only one articulation point
            // Removing all articulation points induces a large regression in the solver score
            //if(find_articulation_points) {
            //	_ = try self.findLargestArticulationPointCut(K,P,solver);
            //}

            // Total contractions left to do since we do not modify the nodes in the list
            var contractions_left: u32 = @intCast(u32, self.nodes.len - 1);
            if (contractions_left == 0) return 0;
            var min_contraction = contraction.Contraction(T){ .erased = 0, .survivor = 0 };
            var total_tww: T = seq.getTwinWidth();

            var node_counter: u32 = 0;
            for (self.nodes) |item| {
                // Remove all nodes that have been erased already before
                if (self.graph.erased_nodes.get(item)) {
                    contractions_left -= 1;
                    continue;
                }
                // Promote to large degree node
                if (self.graph.node_list[item].cardinality() > 10000 and !self.graph.node_list[item].isLargeNode()) {
                    try self.graph.node_list[item].promoteToLargeDegreeNode(self.graph);
                    //std.debug.print("Promoted degree {} to large node.\n",.{self.graph.node_list[item].cardinality()});
                    node_counter += 1;
                }
                try solver.priority_queue.add(item);
            }
            if (contractions_left == 0) return @intCast(T, seq.getTwinWidth());

            solver.bfs_stack.clear();

            // Keep track of the best move we had if we reach the postpone limit
            var best_contraction_potential_postponed: Graph(T).InducedTwinWidthPotential = Graph(T).InducedTwinWidthPotential.default();
            var best_contraction_postponed: contraction.Contraction(T) = .{ .erased = 0, .survivor = 0 };

            var max_consecutive_postpones: u32 = contractions_left;
            var current_postpones: u32 = 0;

            // Disable exhaustive solving for now
            const exhaustive_solving_thresh = 240;

            var enable_follow_up_merge: bool = true;

            while (contractions_left > exhaustive_solving_thresh) {
                var prio_item = try solver.priority_queue.removeNext(total_tww);
                var first_node = prio_item;

                solver.scorer.unsetVisitedBitset(&solver.scratch_bitset);
                var selection = try self.selectBestMove(K, P, first_node, solver, total_tww);

                // Store the best move up to now in the case of multiple postpones
                if (selection.potential.isLess(best_contraction_potential_postponed, total_tww)) {
                    best_contraction_potential_postponed = selection.potential;
                    best_contraction_postponed = contraction.Contraction(T){ .survivor = first_node, .erased = selection.target };
                }

                // If the priority queue proposes to postpone the node we will try to do so
                if (current_postpones < max_consecutive_postpones and solver.priority_queue.shouldPostpone(total_tww, selection.potential.tww) and try solver.priority_queue.postpone(first_node, selection.potential.tww)) {
                    current_postpones += 1;
                    continue;
                }
                // If we could not postpone due to exceeding the maximal postpone limit and the current move is known to be not optimal take
                // the optimal move tried before
                else if (current_postpones >= max_consecutive_postpones and !selection.potential.isLess(best_contraction_potential_postponed, total_tww)) {
                    min_contraction = best_contraction_postponed;
                    selection.potential = best_contraction_potential_postponed;
                } else {
                    // No case before therefore just use the selected move as the next move
                    min_contraction = contraction.Contraction(T){ .survivor = first_node, .erased = selection.target };
                }

								//try self.graph.min_hash.fetchNextMoves(40,self.graph);
								//for(0..self.graph.min_hash.fetched_moves.items.len) |mv_index| {
								//	const mv = self.graph.min_hash.fetched_moves.items[mv_index];

								//	const move = mv.intoMove(T,self.graph.number_of_nodes);
								//	var global_score = self.graph.calculateInducedTwwPotential(move.erased,move.survivor,&selection.potential, seq.getTwinWidth());
								//	std.debug.print("Move {any} tww {}\n",.{mv,global_score.tww});
								//	if(global_score.isLess(selection.potential,seq.getTwinWidth())) {
								//		min_contraction = move;
								//		selection.potential = global_score;
								//	}
								//}
								//try self.graph.min_hash.reinsertFetchedMoves();

                // Reset all variables which were set to select the best postponed move
                current_postpones = 0;
                best_contraction_potential_postponed.reset();

                const contractions_left_copy = contractions_left;

                total_tww = std.math.max(try self.addContractionAndLeafReduction(&enable_follow_up_merge, seq, min_contraction, &contractions_left), total_tww);

                solver.priority_queue.addTick(@intCast(T, contractions_left_copy - contractions_left));

                // Adjust the contractions left and maximal number of postpones allowed (Note this is only problematic in the case of contractions_left < P)
                max_consecutive_postpones = contractions_left + 1;

                // Readd the surviving node to the priority queue
                try solver.priority_queue.add(prio_item);

                // If we exceeded the twin width we can use the fast exit without calculating the moves
                if (total_tww >= contractions_left) {
                    // Contract all remaining vertices in any order
                    var first_left_node: ?T = null;
                    while (contractions_left > 0) {
                        if (first_left_node) |f_node| {
                            const sur = try solver.priority_queue.removeNext(total_tww);
                            if (self.graph.erased_nodes.get(sur)) continue;
                            _ = try self.graph.addContraction(sur, f_node, seq);
                            contractions_left -= 1;
                        } else {
                            first_left_node = try solver.priority_queue.removeNext(total_tww);
                            if (self.graph.erased_nodes.get(first_left_node.?)) {
                                first_left_node = null;
                            }
                        }
                    }

                    //std.debug.print("Time in bfs {} and in calculate tww {}\n",.{select_move_perf_bfs.total_time/(1000*1000),select_move_perf_calc.total_time/(1000*1000)});
                    return total_tww;
                }
            }

            var bounded_array = try std.BoundedArray(T, exhaustive_solving_thresh + 1).init(0);
            for (0..(contractions_left + 1)) |_| {
                try bounded_array.append(try solver.priority_queue.removeNext(total_tww));
            }

            // Do an exhaustive search for the best contraction sequence
            while (contractions_left > 0) {
                var move = self.selectBestMoveExhaustive(exhaustive_solving_thresh + 1, &bounded_array, total_tww);

                total_tww = std.math.max(try self.graph.addContraction(move.erased, move.survivor, seq), total_tww);
                contractions_left -= 1;
                if (total_tww >= contractions_left) {
                    // Contract all remaining vertices in any order
                    var first_left_node: ?T = null;
                    while (contractions_left > 0) {
                        if (first_left_node) |f_node| {
                            const sur = bounded_array.pop();
                            _ = try self.graph.addContraction(sur, f_node, seq);
                            contractions_left -= 1;
                        } else {
                            first_left_node = bounded_array.pop();
                        }
                    }

                    //std.debug.print("Time in bfs {} and in calculate tww {}\n",.{select_move_perf_bfs.total_time/(1000*1000),select_move_perf_calc.total_time/(1000*1000)});
                    return total_tww;
                }
            }

            return total_tww;
        }

        pub fn reduceLeafesAndPaths(self: *Self, comptime K: u32, comptime P: u32, seq: *retraceable_contraction.RetraceableContractionSequence(T), solver: *solver_resources.SolverResources(T, K, P), comptime reduce_paths: bool) !void {
            solver.bfs_stack.clear();
            for (self.nodes) |item| {
                solver.bfs_stack.addNext(item);
            }

            var leafes_reduced: u32 = 0;
            var reduced_paths: u32 = 0;
            while (solver.bfs_stack.nextFrontierSize() > 0) {
                solver.bfs_stack.swapFrontiers();
                var iterator = solver.bfs_stack.iterator();
                while (iterator.next()) |item| {
                    if (self.graph.erased_nodes.get(item)) continue;
                    if (self.graph.node_list[item].cardinality() == 1) {
                        var parent = self.graph.node_list[item].getFirstNeighboor();
                        var iterparent = self.graph.node_list[parent].unorderedIterator();
                        while (iterparent.next()) |childp| {
                            if (item != childp and self.graph.node_list[childp].cardinality() == 1) {
                                _ = try self.graph.addContraction(childp, item, seq);
                                leafes_reduced += 1;
                            }
                        }
                        if (reduce_paths) {
                            var new_item = item;
                            // Ok no leaf here check path and the corresponding length
                            while (self.graph.node_list[parent].cardinality() <= 2) {
                                if (self.graph.node_list[parent].cardinality() == 1) {
                                    _ = try self.graph.addContraction(new_item, parent, seq);
                                    reduced_paths += 1;
                                    break;
                                }
                                var iternext = self.graph.node_list[parent].unorderedIterator();
                                var next_parent_p1 = iternext.next().?;
                                var next_parent_p2 = iternext.next().?;

                                var next_parent = if (next_parent_p1 == new_item) next_parent_p2 else next_parent_p1;
                                if (self.graph.node_list[next_parent].cardinality() == 2) {
                                    _ = try self.graph.addContraction(new_item, parent, seq);
                                    solver.bfs_stack.addNext(parent);
                                    new_item = parent;
                                    parent = next_parent;
                                    reduced_paths += 1;
                                } else {
                                    break;
                                }
                            }
                        }
                    }
                }
            }

            solver.bfs_stack.clear();
        }
    };
}
