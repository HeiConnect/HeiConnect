#ifndef APPROX_HPP
#define APPROX_HPP

#include "graph.hpp"
#include <list>

namespace solver {

std::list<graph::Edge> approximate_2_lp(graph::GraphPair &graph);

std::list<graph::Edge> approximate_1_ln2_e(graph::GraphPair &graph, double epsilon);

std::list<graph::Edge> approximate_1_5_e(graph::GraphPair &graph, double epsilon);

} // namespace solver

#endif // APPROX_HPP
