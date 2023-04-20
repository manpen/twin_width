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
            return Self {
                .nodes = nodes,
                .graph = graph,
            };
        }

				pub const TargetMinimalInducedTww = struct {
					tww: T,
					target: T
				};

				pub fn compareNodeDegrees(ctx: *const Graph(T), lhs: T, rhs: T) std.math.Order {
					return std.math.order(ctx.node_list[rhs].cardinality(),ctx.node_list[lhs].cardinality());
				}


				pub inline fn selectBestMove(self: *Self, comptime K: u32, comptime P: u32, first_node: T, solver: *solver_resources.SolverResources(T,K,P), current_tww: T) !TargetMinimalInducedTww {
					var induced_tww = Graph(T).InducedTwinWidthPotential.default();
					var min_target: T = 0;
					bfs_mod.bfs_topk(T, K, first_node, self.graph, &solver.scratch_bitset, &solver.bfs_stack, &solver.scorer, .{ .max_level = 2, .kind = .both });

					var iterator = solver.scorer.iterator();
					var items_in_loop: u32 = 0;
					while (iterator.next()) |item| {
						if (item == first_node) continue;
						items_in_loop += 1;
						const cal_induced_tww = self.graph.calculateInducedTwwPotential(item, first_node);

						if (cal_induced_tww.isLess(induced_tww, current_tww)) {
							induced_tww = cal_induced_tww;
							min_target = item;
						}
					}

					if(items_in_loop == 0) {
						std.debug.print("Error cannot find contraction partner!\n",.{});
					}
					
					solver.scratch_bitset.unsetAll();
					return TargetMinimalInducedTww {
						.tww = induced_tww.tww,
						.target = min_target
					};
				}

				pub inline fn selectBestMoveExhaustive(self: *Self, comptime N:u32, remaining_nodes: *std.BoundedArray(T,N)) contraction.Contraction(T) {
					var induced_tww = Graph(T).InducedTwinWidth.default();
					var min_contraction = contraction.Contraction(T) {.erased = 0, .survivor = 0};
					var min_index_erased:T = 0;

					for(0..remaining_nodes.len) |first| {
						for(0..first) |second| {
							const first_node = remaining_nodes.buffer[first];
							const second_node = remaining_nodes.buffer[second];
							var ind = self.graph.calculateInducedTww(first_node,second_node, induced_tww.tww);
							if (ind.lessThan(induced_tww)) {
								induced_tww = ind;
								min_contraction = .{.erased = first_node, .survivor = second_node};
								min_index_erased = @intCast(T,first);
							}
						}
					}

					_ = remaining_nodes.swapRemove(min_index_erased);

					return min_contraction;
				}

				
				pub fn compareTValue(ctx: void, lhs: T, rhs: T) std.math.Order {
				    _ = ctx;
					return std.math.order(rhs,lhs);
				}

				pub const LargestArticulationPoint = struct {
					point: T,
					second_largest_component_size: T,
				};

				pub fn findLargestArticulationPointCut(self: *Self, comptime K:u32, comptime P:u32, solver: *solver_resources.SolverResources(T,K,P)) !?T {
            solver.scratch_bitset.unsetAll();
						try dfs_mod.findArticulationPoints(T, self.nodes[0], self.graph, &solver.scratch_bitset, &solver.articulation_point_resources);
						solver.scratch_bitset.unsetAll();

						// Reset bitset we need it to use it as a scratch for the visited field
						if(solver.articulation_point_resources.articulation_points.cardinality == 0) {
							return null;
						}

						var priority_queue_components = std.PriorityQueue(T,void,Self.compareTValue).init(self.graph.allocator,{});
						defer priority_queue_components.deinit();

						try priority_queue_components.ensureTotalCapacity(self.graph.number_of_nodes);


						var largest = LargestArticulationPoint {
							.point = 0,
							.second_largest_component_size = 0
						};

						var iter_art = solver.articulation_point_resources.articulation_points.iter();
						while(iter_art.next()) |item| {
							var nb_iter = self.graph.node_list[item].unorderedIterator();
							var total_count:T = 0;

							while(nb_iter.next()) |nb| {
								// Ok new connected component
								// Check root node
								if(solver.articulation_point_resources.parent_array[nb] == item) {
									if(item == self.nodes[0]) {
										// Root nodes will cover everything with that
										try priority_queue_components.add(solver.articulation_point_resources.node_count[nb]);
									}
									else if(solver.articulation_point_resources.low_array[nb] >= solver.articulation_point_resources.depth_array[item]) {
										try priority_queue_components.add(solver.articulation_point_resources.node_count[nb]);
										total_count += solver.articulation_point_resources.node_count[nb];
									}
								}
							}

							if(total_count != 0) {
								// Add the left over connected component
								try priority_queue_components.add(@intCast(T,self.nodes.len)-(total_count+1));
							}

							_ = priority_queue_components.remove();
							const second_largest = priority_queue_components.remove();
							if(second_largest > largest.second_largest_component_size) {
								largest.point = @intCast(T,item);
								largest.second_largest_component_size = second_largest;
							}
							
							// Delete all items
							priority_queue_components.len = 0;
						}
						
						return largest.point;
        }

        pub fn solveGreedyTopK(self: *Self, comptime K: u32, comptime P:u32, seq: *retraceable_contraction.RetraceableContractionSequence(T), solver: *solver_resources.SolverResources(T,K,P), find_articulation_points: bool) !T {
            	_ = find_articulation_points;
							// Use once at the beginning later we will only consider merges which involve the survivor of the last merge
            	try self.reduceLeafesAndPaths(K,P,seq,solver);

							// Reset bitset we need it to use it as a scratch for the visited field
							solver.scratch_bitset.unsetAll();

							// Sadly only small components can be split of while removing only one articulation point
							// Removing all articulation points induces a large regression in the solver score
							//if(find_articulation_points) {
							//	_ = try self.findLargestArticulationPointCut(K,P,solver);
							//}
							

							// Total contractions left to do since we do not modify the nodes in the list
							var contractions_left: u32 = @intCast(u32,self.nodes.len - 1);
							if(contractions_left == 0) return 0;
							var min_contraction = contraction.Contraction(T){ .erased = 0, .survivor = 0 };
							var total_tww: T = seq.getTwinWidth();

							var node_counter:u32 = 0;
							for(self.nodes) |item| {
								if (self.graph.erased_nodes.get(item)) {
									contractions_left-=1;
									continue;
								}
								// Promote to large degree node
								if(self.graph.node_list[item].cardinality() > 10000 and !self.graph.node_list[item].isLargeNode()) {
									try self.graph.node_list[item].promoteToLargeDegreeNode(self.graph);
									//std.debug.print("Promoted degree {} to large node.\n",.{self.graph.node_list[item].cardinality()});
									node_counter+=1;
								}
								try solver.priority_queue.add(item);
							}
							if(contractions_left == 0) return @intCast(T,seq.getTwinWidth());

							solver.bfs_stack.clear();

							// Keep track of the best move we had if we reach the postpone limit
							var best_contraction_tww_postponed:T = std.math.maxInt(T);
							var best_contraction_postponed: contraction.Contraction(T)  = .{.erased = 0, .survivor = 0};


							var max_consecutive_postpones:u32 = contractions_left;
							var current_postpones:u32 = 0;

							while (contractions_left > 0) {
								var prio_item = try solver.priority_queue.removeNext(total_tww);
								var first_node = prio_item;

								var selection = try self.selectBestMove(K,P,first_node,solver,total_tww);


								// Store the best move up to now in the case of multiple postpones
								if(selection.tww < best_contraction_tww_postponed) {
									best_contraction_tww_postponed = selection.tww;
									best_contraction_postponed = contraction.Contraction(T) {.survivor = first_node, .erased = selection.target};
								}

								// If the priority queue proposes to postpone the node we will try to do so
								if(current_postpones < max_consecutive_postpones and solver.priority_queue.shouldPostpone(total_tww,selection.tww) and try solver.priority_queue.postpone(first_node,selection.tww)) {
									current_postpones+=1;
									continue;
								}
								// If we could not postpone due to exceeding the maximal postpone limit and the current move is known to be not optimal take
								// the optimal move tried before
								else if(current_postpones >= max_consecutive_postpones and selection.tww > best_contraction_tww_postponed) {
									min_contraction = best_contraction_postponed;	
								}
								else {
									// No case before therefore just use the selected move as the next move
									min_contraction = contraction.Contraction(T) {.survivor = first_node, .erased = selection.target};
								}
								// Reset all variables which were set to select the best postponed move
								current_postpones = 0;
								best_contraction_tww_postponed = std.math.maxInt(T);

								// Adjust the contractions left and maximal number of postpones allowed (Note this is only problematic in the case of contractions_left < P)
								contractions_left -= 1;
								max_consecutive_postpones = contractions_left+1;

								var last_evoked_twin_width = try self.graph.addContraction(min_contraction.erased, min_contraction.survivor, seq);
								total_tww = std.math.max(last_evoked_twin_width,total_tww);

								// Increase the tick to inform the priority queue about the current time step so that it can
								// release postponed nodes
								solver.priority_queue.increaseTick();
								
								// Reduce leafes if one was created
								if(self.graph.node_list[min_contraction.survivor].isLeaf()) {
									// Ok is leaf check for other leaf
									var parent = self.graph.node_list[min_contraction.survivor].getFirstNeighboor();
									if(self.graph.node_list[parent].num_leafes > 1) {
										var parent_node_iter = self.graph.node_list[parent].unorderedIterator();
										while(parent_node_iter.next()) |item| {
											if(item == min_contraction.survivor) continue;
											if(self.graph.node_list[item].isLeaf()) {
												last_evoked_twin_width = try self.graph.addContraction(item, min_contraction.survivor, seq);
												total_tww = std.math.max(last_evoked_twin_width,total_tww);
												solver.priority_queue.increaseTick();
												contractions_left-=1;
												// Only merge one other leaf since we do this for every leaf and there should not be any more leafes left
												// WARNING: The iterator probably break if we try to go on since addContraction will update the edge list
												break;
											}
										}
									}
								}
								// Reduce leaves if we are the parent of at least 2 leafes
								else if(self.graph.node_list[min_contraction.survivor].num_leafes > 1) {
									var parent_node_iter = self.graph.node_list[min_contraction.survivor].unorderedIterator();
									var first_leaf: ?T = null;

									while(parent_node_iter.next()) |item| {
										if(self.graph.node_list[item].isLeaf()) {
											if(first_leaf) |k| {
												last_evoked_twin_width = try self.graph.addContraction(k, item, seq);
												total_tww = std.math.max(last_evoked_twin_width,total_tww);
												solver.priority_queue.increaseTick();
												contractions_left-=1;
												break;
											}
											else {
												first_leaf = item;
											}
										}
									}
								}

								// Readd the surviving node to the priority queue
								try solver.priority_queue.add(prio_item);

								// If we exceeded the twin width we can use the fast exit without calculating the moves
								if (total_tww >= contractions_left) {
									// Contract all remaining vertices in any order
									var first_left_node:?T = null;
									while(contractions_left > 0) {
										if(first_left_node) |f_node| {
											const sur = try solver.priority_queue.removeNext(total_tww);
											if(self.graph.erased_nodes.get(sur)) continue;
											_ = try self.graph.addContraction(sur,f_node, seq);
											contractions_left-=1;
										}
										else {
											first_left_node = try solver.priority_queue.removeNext(total_tww);
											if(self.graph.erased_nodes.get(first_left_node.?)) {
												first_left_node = null;	
											}
										}
									}
									return total_tww;
								}
							}

							return total_tww;
        }

        pub fn reduceLeafesAndPaths(self: *Self, comptime K: u32, comptime P:u32, seq: *retraceable_contraction.RetraceableContractionSequence(T), solver: *solver_resources.SolverResources(T,K,P)) !void {
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

                solver.bfs_stack.clear();
        }
    };
}
