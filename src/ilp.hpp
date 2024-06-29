#ifndef ILP_HPP
#define ILP_HPP

#include <list>

#include "graph.hpp"

namespace solver {
std::list<graph::Edge> ilp(graph::GraphPair &graph, bool use_initial = false,
                           int presolve = 0);
}

#endif // ILP_HPP
