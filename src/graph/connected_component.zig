const std = @import("std");
const edge_list = @import("edge_list.zig");
const contraction = @import("../tww/contraction_sequence.zig");
const graph_mod = @import("graph.zig");
const comptime_util = @import("../util/comptime_checks.zig");
const bfs_mod = @import("bfs.zig");
const CompactField = graph_mod.CompactField;
const bench_timer = @import("../util/benchmark_helper.zig");

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

							 pub fn solveGreedy(self: *Self, graph: *graph_mod.Graph(T)) !T {
								 graph.scratch_bitset.unsetAll();
								 var bfs_stack = try bfs_mod.BfsQueue(T).init(graph.allocator, graph.number_of_nodes);
								 var contractions:u32 = self.nodes.cardinality()-1;

								 var min_contraction = contraction.Contraction(T){ .erased = 0, .survivor = 0 };
								 var min_tww:?T = null;
								 var reduced_red_edges:T = 0;
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
										 try timer_bfs.stop();
										 try timer_induced.start();
										 const induced_tww = graph.calculateInducedTww(item,first_node,min_tww);
										 try timer_induced.stop();
										 try timer_bfs.start();
										 if(min_tww == null or induced_tww.tww < min_tww.? or (induced_tww.tww == min_tww.? and induced_tww.erased_red_edges > reduced_red_edges)) {
											 min_tww = induced_tww.tww;
											 reduced_red_edges = induced_tww.erased_red_edges;
											 min_contraction = contraction.Contraction(T){ .erased = item, .survivor = first_node };
										 }
										 if(induced_tww.tww == 0) {
											 break;
										 }
									 }
									 try timer_bfs.stop();

									 contractions-=1;
									 min_tww = null;
									 reduced_red_edges = 0;

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
										 //std.debug.print("Time BFS {}ms time induced {}ms and time add contractions {}ms\n",.{timer_bfs.total_time/(1000*1000),timer_induced.total_time/(1000*1000),timer_add_cont.total_time/(1000*1000)});
									 }
									 graph.scratch_bitset.unsetAll();
								 }
								 self.tww = total_tww;
								 return total_tww;
							 }

        pub fn init(allocator: std.mem.Allocator, nodes: edge_list.ParametrizedUnsortedArrayList(T), bfs_level: u32) !Self {
            var contraction_seq = try contraction.ContractionSequence(T).init(allocator, nodes.cardinality());
            var iter = nodes.iterator();
            var previous: ?T = null;
            while (iter.next()) |item| {
                if (previous) |into| {
                    try contraction_seq.addContraction(contraction.Contraction(T){ .erased = item, .survivor = into });
                } else {
                    previous = item;
                }
            }
            return .{ .bfs_levels = bfs_level, .nodes = nodes, .tww = @intCast(T,nodes.cardinality() - 1), .best_contraction_sequence = contraction_seq };
        }
    };
}
