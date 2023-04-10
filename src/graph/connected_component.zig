const std = @import("std");
const edge_list = @import("edge_list.zig");
const contraction = @import("../tww/contraction_sequence.zig");
const graph_mod = @import("graph.zig");
const comptime_util = @import("../util/comptime_checks.zig");
const bfs_mod = @import("bfs.zig");
const CompactField = graph_mod.CompactField;
const bench_timer = @import("../util/benchmark_helper.zig");
const topk = @import("../util/top_k_scorer.zig");

pub fn ConnectedComponentIndex(comptime T: type) type {
    comptime if (!comptime_util.checkIfIsCompatibleInteger(T)) {
        @compileError("The type of ConnectedComponent must either be u8,u16 or u32!");
    };
    return struct {
        const Self = @This();
        tww: T,
        index: T,
        pub fn compareComponentIndexDesc(ctx: void, lhs: Self, rhs: Self) std.math.Order {
            _ = ctx;
            return std.math.order(rhs.tww,lhs.tww);
        }
    };
}

pub fn ConnectedComponent(comptime T: type) type {
    comptime if (!comptime_util.checkIfIsCompatibleInteger(T)) {
        @compileError("The type of ConnectedComponent must either be u8,u16 or u32!");
    };

    return struct {
        const Self = @This();
        nodes: edge_list.ParametrizedUnsortedArrayList(T),
				erased_nodes: std.bit_set.DynamicBitSetUnmanaged,
        best_contraction_sequence: contraction.ContractionSequence(T),
        tww: T,
        bfs_levels: u32,

				pub fn deinit(self: *Self, allocator: std.mem.Allocator) void {
					self.nodes.deinit(allocator);
					self.best_contraction_sequence.deinit(allocator);
				}

				pub const NodePriority = struct {
					red_degree: T,
					black_degree: T,
					id: T,
					pub fn compareFn(ctx: void, lhs: NodePriority, rhs: NodePriority) std.math.Order {
					    _ = ctx;
						if(lhs.red_degree == rhs.red_degree) {
							return std.math.order(lhs.black_degree,rhs.red_degree);
						}
						return std.math.order(lhs.red_degree,rhs.red_degree);
					}
				};

				fn reduceLeafes(self: *Self, graph: *graph_mod.Graph(T), bfs_queue: *bfs_mod.BfsQueue(T)) !void {
					bfs_queue.clear();
					var node_iterator = self.nodes.iterator();
					while(node_iterator.next()) |item| {
						bfs_queue.addNext(item);
					}

					var leafes_reduced:u32 = 0;
					while(bfs_queue.nextFrontierSize() > 0) {
						bfs_queue.swapFrontiers();
						var iterator = bfs_queue.iterator();
outer: while(iterator.next()) |item| {
				 if(self.erased_nodes.isSet(item)) continue;
				 if(graph.node_list[item].black_edges.cardinality() == 1) {
					 var neighbors = graph.node_list[item].black_edges.iterator();
					 var parent = neighbors.next().?;
					 var iterparent = graph.node_list[parent].black_edges.iterator();
					 while(iterparent.next()) |childp| {
						 if(item != childp and graph.node_list[childp].black_edges.cardinality() == 1) {
							 _ = try graph.addContraction(item,childp);
							 bfs_queue.addNext(childp);
							 self.erased_nodes.set(item);
							 leafes_reduced+=1;
continue :outer;
						 }
					 }
				 }
			 }

					}

					bfs_queue.clear();	
					std.debug.print("Total leafes reduced {}\n",.{leafes_reduced});
				}

							 pub fn solveGreedyTopK(self: *Self, graph: *graph_mod.Graph(T)) !T {
								 var bfs_stack = try bfs_mod.BfsQueue(T).init(graph.allocator, graph.number_of_nodes);
							 		try self.reduceLeafes(graph,&bfs_stack);

								 graph.scratch_bitset.unsetAll();
								 var scorer = try topk.TopKScorer(T,10).init(graph.allocator,@intCast(T,graph.number_of_nodes));

								 var contractions:u32 = self.nodes.cardinality()-(1+@intCast(T,self.erased_nodes.count()));

								 var min_contraction = contraction.Contraction(T){ .erased = 0, .survivor = 0 };
								 var min_tww:?T = null;
								 var delta_red_edges:i32 = std.math.maxInt(i32);
								 var total_tww:T = 0;


								 var priority_queue = std.PriorityQueue(NodePriority,void,NodePriority.compareFn).init(graph.allocator,{});
								 defer priority_queue.deinit();

								 var node_iterator = self.nodes.iterator();
								 var largest_deg:T = 0;
								 var deg1_nodes:T = 0;
								 var leafes:T = 0;
								 while(node_iterator.next()) |item| {
								 	if(self.erased_nodes.isSet(item)) continue;
								 	try priority_queue.add(NodePriority {
										.red_degree = graph.node_list[item].red_edges.cardinality(),
										.black_degree = graph.node_list[item].black_edges.cardinality(),
										.id = item
									});
									largest_deg = std.math.max(graph.node_list[item].black_edges.cardinality(),largest_deg);
									if(graph.node_list[item].black_edges.cardinality()<=1) {
										var neighbors = graph.node_list[item].black_edges.iterator();
										var parent = neighbors.next().?;
										var iterparent = graph.node_list[parent].black_edges.iterator();
										while(iterparent.next()) |childp| {
											if(item != childp and graph.node_list[childp].black_edges.cardinality() == 1) {
												leafes+=1;
											}
										}
										deg1_nodes+=1;
									}
								 }
								 std.debug.print("largest deg: {d} deg one {} leafes {}",.{largest_deg,deg1_nodes,leafes});

								var timer_bfs = bench_timer.BenchmarkHelper.init();
								var timer_induced = bench_timer.BenchmarkHelper.init();
								var timer_add_cont = bench_timer.BenchmarkHelper.init();

								 while(contractions > 0) {
								 	 var prio_item = priority_queue.remove();
									 if(graph.node_list[prio_item.id].red_edges.cardinality() != prio_item.red_degree) {
									 	prio_item.red_degree = graph.node_list[prio_item.id].red_edges.cardinality();
										prio_item.black_degree = graph.node_list[prio_item.id].black_edges.cardinality();
										try priority_queue.add(prio_item);
										continue;
									 }
									 else if(graph.erased_nodes.get(prio_item.id)) {
									 	continue;
									 }
									 const first_node = prio_item.id;

										 scorer.setNextTargetDegree(graph.node_list[prio_item.id].black_edges.cardinality()+graph.node_list[prio_item.id].red_edges.cardinality());

										 try timer_bfs.start();
										 bfs_mod.bfs_topk(T,10,first_node,graph,&graph.scratch_bitset,&bfs_stack,&scorer,.{.max_level = 2, .kind=.both});
										 try timer_bfs.stop();

										 try timer_induced.start();
										 var iterator = scorer.iterator();
										 var items_in_loop:u32 = 0;
										 while(iterator.next()) |item| {
											 if(item == first_node) continue;
											 items_in_loop+=1;
											 const induced_tww = graph.calculateInducedTww(item,first_node,min_tww);

											 if(min_tww == null or induced_tww.tww < min_tww.? or (induced_tww.tww == min_tww.? and induced_tww.delta_red_edges < delta_red_edges)) {
												 min_tww = induced_tww.tww;
												 delta_red_edges = induced_tww.delta_red_edges;
												 min_contraction = contraction.Contraction(T){ .erased = item, .survivor = first_node };
												 if(induced_tww.tww == 0) {
													 break;
												 }
											 }
										 }
										 try timer_induced.stop();
									 	 graph.scratch_bitset.unsetAll();

									 contractions-=1;
									 min_tww = null;
									 delta_red_edges = std.math.maxInt(i32);

									 //if(contractions%100 == 0) {
									 //	std.debug.print("Size of black {}::{} and red {}::{}\n",.{graph.node_list[min_contraction.erased].black_edges.cardinality(),graph.node_list[min_contraction.survivor].black_edges.cardinality(), graph.node_list[min_contraction.erased].red_edges.cardinality(), graph.node_list[min_contraction.survivor].red_edges.cardinality()});
									 //}
									 try timer_add_cont.start();
									 total_tww = std.math.max(try graph.addContraction(min_contraction.erased,min_contraction.survivor),total_tww);
									 try timer_add_cont.stop();
									 prio_item.red_degree = graph.node_list[first_node].red_edges.cardinality();
									 prio_item.black_degree = graph.node_list[first_node].black_edges.cardinality();
									 try priority_queue.add(prio_item);
									 if(contractions%1000 == 0) {
										 std.debug.print("Remaining nodes {}\n",.{contractions});
										 //std.debug.print("Time BFS {}ms time induced {}ms and time add contractions {}ms total edges followed {}\n",.{timer_bfs.total_time/(1000*1000),timer_induced.total_time/(1000*1000),timer_add_cont.total_time/(1000*1000),scorer.add_visit_calls});
									 }

									 if(total_tww >= contractions) {
									 	std.debug.print("Early stopping due to high twin width!\n",.{});
										self.tww = total_tww;
										return total_tww;
									 }
								 }
								 self.tww = total_tww;
								 return total_tww;
							 }

							 pub fn solveGreedy(self: *Self, graph: *graph_mod.Graph(T)) !T {
								 graph.scratch_bitset.unsetAll();
								 var bfs_stack = try bfs_mod.BfsQueue(T).init(graph.allocator, graph.number_of_nodes);
								 var contractions:u32 = self.nodes.cardinality()-1;

								 var min_contraction = contraction.Contraction(T){ .erased = 0, .survivor = 0 };
								 var min_tww:?T = null;
								 var delta_red_edges:i32 = std.math.maxInt(i32);
								 var total_tww:T = 0;

								 var priority_queue = std.PriorityQueue(NodePriority,void,NodePriority.compareFn).init(graph.allocator,{});
								 defer priority_queue.deinit();

								 var node_iterator = self.nodes.iterator();
								 while(node_iterator.next()) |item| {
								 	try priority_queue.add(NodePriority {
										.red_degree = graph.node_list[item].red_edges.cardinality(),
										.black_degree = graph.node_list[item].black_edges.cardinality(),
										.id = item
									});
								 }

									var timer_bfs = bench_timer.BenchmarkHelper.init();
									var timer_induced = bench_timer.BenchmarkHelper.init();
									var timer_add_cont = bench_timer.BenchmarkHelper.init();
									 var total_bfs_size:u64 = 0;
								 while(contractions > 0) {
								 	 var prio_item = priority_queue.remove();
									 if(graph.node_list[prio_item.id].red_edges.cardinality() != prio_item.red_degree) {
									 	prio_item.red_degree = graph.node_list[prio_item.id].red_edges.cardinality();
										prio_item.black_degree = graph.node_list[prio_item.id].black_edges.cardinality();
										try priority_queue.add(prio_item);
										continue;
									 }
									 else if(graph.erased_nodes.get(prio_item.id)) {
									 	continue;
									 }
									 const first_node = prio_item.id;

									 try timer_bfs.start();
									 var iterator = bfs_mod.bfs(T,first_node,graph,&graph.scratch_bitset,&bfs_stack,.{.max_level = 2, .kind=.both});
									 while(iterator.next()) |item| {
										 if(item == first_node) continue;
										 total_bfs_size+=1;
										 try timer_bfs.stop();
										 try timer_induced.start();
										 const induced_tww = graph.calculateInducedTww(item,first_node,min_tww);
										 try timer_induced.stop();
										 try timer_bfs.start();
										 if(min_tww == null or induced_tww.tww < min_tww.? or (induced_tww.tww == min_tww.? and induced_tww.delta_red_edges < delta_red_edges)) {
											 min_tww = induced_tww.tww;
											 delta_red_edges = induced_tww.delta_red_edges;
											 min_contraction = contraction.Contraction(T){ .erased = item, .survivor = first_node };
											 if(induced_tww.tww == 0) {
												 break;
											 }
										 }
									 }
									 try timer_bfs.stop();

									 contractions-=1;
									 min_tww = null;
									 delta_red_edges = std.math.maxInt(i32);

									 //if(contractions%100 == 0) {
									 //	std.debug.print("Size of black {}::{} and red {}::{}\n",.{graph.node_list[min_contraction.erased].black_edges.cardinality(),graph.node_list[min_contraction.survivor].black_edges.cardinality(), graph.node_list[min_contraction.erased].red_edges.cardinality(), graph.node_list[min_contraction.survivor].red_edges.cardinality()});
									 //}
									 try timer_add_cont.start();
									 total_tww = std.math.max(try graph.addContraction(min_contraction.erased,min_contraction.survivor),total_tww);
									 try timer_add_cont.stop();
									 prio_item.red_degree = graph.node_list[first_node].red_edges.cardinality();
									 prio_item.black_degree = graph.node_list[first_node].black_edges.cardinality();
									 try priority_queue.add(prio_item);
									 if(contractions%1000 == 0) {
									 		std.debug.print("Total bfs size {}\n",.{total_bfs_size/1000});
									 		total_bfs_size = 0;
										 //std.debug.print("Remaining nodes {}\n",.{contractions});
										 //std.debug.print("Time BFS {}ms time induced {}ms and time add contractions {}ms\n",.{timer_bfs.total_time/(1000*1000),timer_induced.total_time/(1000*1000),timer_add_cont.total_time/(1000*1000)});
									 }
									 graph.scratch_bitset.unsetAll();
								 }
								 self.tww = total_tww;
								 return total_tww;
							 }

        pub fn init(allocator: std.mem.Allocator, nodes: edge_list.ParametrizedUnsortedArrayList(T), bfs_level: u32) !Self {
            var contraction_seq = try contraction.ContractionSequence(T).init(allocator, nodes.cardinality());
            var erased_nodes = try std.bit_set.DynamicBitSetUnmanaged.initEmpty(allocator, nodes.cardinality());
            var iter = nodes.iterator();
            var previous: ?T = null;
            while (iter.next()) |item| {
                if (previous) |into| {
                    try contraction_seq.addContraction(contraction.Contraction(T){ .erased = item, .survivor = into });
                } else {
                    previous = item;
                }
            }
            return .{ .bfs_levels = bfs_level, .nodes = nodes, .tww = @intCast(T,nodes.cardinality() - 1), .best_contraction_sequence = contraction_seq, .erased_nodes = erased_nodes };
        }
    };
}
