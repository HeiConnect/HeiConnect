#ifndef GREEDY_HPP
#define GREEDY_HPP

#include "graph.hpp"

namespace solver {
std::list<graph::Edge> greedy_weak(graph::GraphPair &graph);
std::list<graph::Edge> greedy_strong(graph::GraphPair &graph);
std::list<graph::Edge> greedy_heuristic_strong(graph::GraphPair &g);
std::list<graph::Edge> greedy_heuristic_sampling(graph::GraphPair &g,
                                                 int buckets);
std::list<graph::Edge> greedy_dynamic(graph::DynamicCactus &g);
std::list<graph::Edge> greedy_dynamic_sampling(graph::DynamicCactus &g);
std::list<graph::Edge> greedy_dynamic_bounds(graph::DynamicCactus &g);
std::list<graph::Edge> greedy_mst(graph::GraphPair &g);
std::list<graph::Edge> greedy_mst2(graph::GraphPair &g);
std::list<graph::Edge> greedy_mst_ilp(graph::GraphPair &g, int trees);
std::list<graph::Edge> greedy_mst_improved(graph::GraphPair &g,
                                           graph::DynamicCactus cactus);
std::pair<std::list<graph::Edge>, std::list<graph::Edge>>
greedy_mst_max_flow(graph::GraphPair &g);
std::list<graph::Edge>
greedy_mst_max_flow_heuristic(graph::GraphPair &g,
                              graph::DynamicCactus &dyn_cactus);
std::list<graph::Edge>
greedy_mst_max_flow_order_heuristic(graph::GraphPair &g,
                                    graph::DynamicCactus &dyn_cactus);
    std::list<graph::Edge> full_mst(graph::GraphPair &g);
std::list<graph::Edge> greedy_2mst_localsearch(graph::GraphPair &g,
                                               int depth_limit);
std::list<graph::Edge> greedy_2mst_localsearch_flow(graph::GraphPair &g,
                                                    int depth_limit, bool cache,
                                                    int trees = 2);
} // namespace solver

#endif // GREEDY_HPP