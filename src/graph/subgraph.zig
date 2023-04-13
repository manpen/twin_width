const std = @import("std");
const Graph = @import("graph.zig").Graph;
const edge_list = @import("edge_list.zig");
const bitset = @import("../util/two_level_bitset.zig");
const bfs_mod = @import("bfs.zig");
const topk = @import("../util/top_k_scorer.zig");
const contraction = @import("../tww/contraction_sequence.zig");
const node_priority_queue = @import("node_priority_queue.zig");
const solver_resources = @import("solver.zig");

pub const SubGraphError = error{
    SolverResourcesNotInitialized,
    SolverResourcesAlreadyInitialized,
};

pub fn InducedSubGraph(comptime T: type) type {
    return struct {
        const Self = @This();
        graph: *Graph(T),
        // Must be this type to support fast lookups!
        nodes: edge_list.ParametrizedUnsortedArrayList(T),

        pub fn fromUnsortedArrayList(graph: *Graph(T), nodes: edge_list.ParametrizedUnsortedArrayList(T)) !Self {
            return Self {
                .nodes = nodes,
                .graph = graph,
            };
        }

        pub fn deinit(self: *Self) void {
            self.nodes.deinit(self.graph.allocator);
        }

				pub const TargetMinimalInducedTww = struct {
					tww: T,
					target: T
				};

				pub fn compareNodeDegrees(ctx: *const Graph(T), lhs: T, rhs: T) std.math.Order {
					return std.math.order(ctx.node_list[rhs].cardinality(),ctx.node_list[lhs].cardinality());
				}

        pub fn solveGreedyTopK(self: *Self, comptime K: u32, comptime P:u32, solver: *solver_resources.SolverResources(T,K,P)) !T {
            try self.reduceLeafesAndPaths(K,P,solver);
							// Reset bitset we need it to use it as a scratch for the visited field
							solver.scratch_bitset.unsetAll();


							var contractions_left: u32 = self.nodes.cardinality() - 1;
							if(contractions_left == 0) return 0;
							var min_contraction = contraction.Contraction(T){ .erased = 0, .survivor = 0 };
							var induced_tww = Graph(T).InducedTwinWidth.default();
							var total_tww: T = 0;

							var node_iterator = self.nodes.iterator();
							var node_counter:u32 = 0;
							while (node_iterator.next()) |item| {
								if (self.graph.erased_nodes.get(item)) {
									contractions_left-=1;
									continue;
								}
								if(self.graph.node_list[item].cardinality() > 10000) {
									try self.graph.node_list[item].promoteToLargeDegreeNode(self.graph);
									//std.debug.print("Promoted degree {} to large node.\n",.{self.graph.node_list[item].cardinality()});
									node_counter+=1;
								}
								try solver.priority_queue.add(item);
							}
							if(contractions_left == 0) return @intCast(T,self.graph.getCurrentTwinWidth());

							solver.bfs_stack.clear();


							var max_consecutive_postpones:u32 = contractions_left;
							var current_postpones:u32 = 0;
							while (contractions_left > 0) {
								
								var prio_item = try solver.priority_queue.removeNext(total_tww);
								//if(prio_item == largest_degree_id) {
								//	prio_item = try priority_queue.removeNext();
								//	try priority_queue.add(@intCast(T,largest_degree_id));
								//}
								var first_node = prio_item;
								//solver.scratch_bitset.set(largest_degree_id);
								bfs_mod.bfs_topk(T, K, first_node, self.graph, &solver.scratch_bitset, &solver.bfs_stack, &solver.scorer, .{ .max_level = 2, .kind = .both });

								var iterator = solver.scorer.iterator();
								var items_in_loop: u32 = 0;
								while (iterator.next()) |item| {
									if (item == first_node) continue;
									items_in_loop += 1;
									const cal_induced_tww = self.graph.calculateInducedTww(item, first_node, induced_tww.tww);

									if (cal_induced_tww.lessThan(induced_tww)) {
										induced_tww = cal_induced_tww;
										min_contraction = contraction.Contraction(T){ .erased = item, .survivor = first_node };
									}
								}
								solver.scratch_bitset.unsetAll();
								// Extremely large tww 
								//if(total_tww < (contractions_left>>3) and induced_tww.tww > total_tww and induced_tww.tww > (contractions_left>>3)) {
									// Optimal time to not increase twin width beyond the current tww.
							//		if(try priority_queue.postpone(first_node,induced_tww.tww)) {
							//				induced_tww.reset();
							//				continue;
						//			}
							//	}

								if(current_postpones < max_consecutive_postpones and solver.priority_queue.shouldPostpone(total_tww,induced_tww.tww) and try solver.priority_queue.postpone(first_node,induced_tww.tww)) {
									induced_tww.reset();
									current_postpones+=1;
									continue;
								}
								current_postpones = 0;

								contractions_left -= 1;
								max_consecutive_postpones = contractions_left+1;
								induced_tww.reset();

								var last_evoked_twin_width = try self.graph.addContraction(min_contraction.erased, min_contraction.survivor);
								total_tww = std.math.max(last_evoked_twin_width,total_tww);
								solver.priority_queue.increaseTick();

								try solver.priority_queue.add(first_node);
								if (contractions_left % 100000 == 0) {
									var node_promoter_iterator = self.nodes.iterator();
									while(node_promoter_iterator.next()) |item| {
										if(self.graph.erased_nodes.get(item)) continue;
										if((!self.graph.node_list[item].isLargeNode()) and self.graph.node_list[item].cardinality() > 10000) {
											try self.graph.node_list[item].promoteToLargeDegreeNode(self.graph);
											//std.debug.print("Promoted degree {} to large node.\n",.{self.graph.node_list[item].cardinality()});
										}
									}
								}

								//if (contractions_left % 2000 == 0) {
								//	var bfs_max = scorer.consumeLargestBfsIncreasor();
							//		std.debug.print("Largest bfs increasor {} is postponed? {}\n",.{bfs_max,priority_queue.postpone_set.isSet(bfs_max)});
									//if(!priority_queue.postpone_set.isSet(bfs_max)) {
									//	var tww = try self.calculateMinimalInducedTww(K,&scorer,bfs_max);
									//	solver.scratch_bitset.unsetAll();
									//	if(total_tww >= tww.tww) {
									//		// Node seems to be egligable
									//		std.debug.print("Does not increase TWW do contraction!\n",.{});
									//		contractions_left -= 1;
									//		induced_tww.reset();

//											total_tww = std.math.max(try self.graph.addContraction(tww.target,bfs_max), total_tww);
	//										priority_queue.increaseTick();
		//								}
			//						}
							//	}

								if (total_tww >= contractions_left) {
									//std.debug.print("Early stopping due to high twin width!\n", .{});
									var first_left_node:?T = null;
									while(contractions_left > 0) {
										if(first_left_node) |f_node| {
											const sur = try solver.priority_queue.removeNext(std.math.maxInt(T));
											_ = try self.graph.addContraction(sur,f_node);
											contractions_left-=1;
										}
										else {
											first_left_node = try solver.priority_queue.removeNext(std.math.maxInt(T));
										}
									}
									return total_tww;
								}
							}
							return total_tww;
        }

        pub fn reduceLeafesAndPaths(self: *Self, comptime K: u32, comptime P:u32, solver: *solver_resources.SolverResources(T,K,P)) !void {
                solver.bfs_stack.clear();
                var node_iterator = self.nodes.iterator();
                while (node_iterator.next()) |item| {
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
                            var neighbors = self.graph.node_list[item].orderedIterator();
                            var parent = neighbors.next().?;
                            var iterparent = self.graph.node_list[parent].orderedIterator();
                            while (iterparent.next()) |childp| {
                                if (item != childp and self.graph.node_list[childp].cardinality() == 1) {
                                    _ = try self.graph.addContraction(childp, item);
                                    leafes_reduced += 1;
                                }
                            }
                            var new_item = item;
                            // Ok no leaf here check path and the corresponding length
                            while (self.graph.node_list[parent].cardinality() <= 2) {
                                if (self.graph.node_list[parent].cardinality() == 1) {
                                    _ = try self.graph.addContraction(new_item, parent);
                                    reduced_paths += 1;
                                    break;
                                }
                                var iternext = self.graph.node_list[parent].orderedIterator();
                                var next_parent_p1 = iternext.next().?;
                                var next_parent_p2 = iternext.next().?;

                                var next_parent = if (next_parent_p1 == new_item) next_parent_p2 else next_parent_p1;
                                if (self.graph.node_list[next_parent].cardinality() == 2) {
                                    _ = try self.graph.addContraction(new_item, parent);
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

                solver.bfs_stack.clear();
        }
    };
}
