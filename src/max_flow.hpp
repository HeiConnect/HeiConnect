#ifndef MAX_FLOW_HPP
#define MAX_FLOW_HPP

#include <algorithm>
#include <iostream>
#include <queue>
#include <vector>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/boykov_kolmogorov_max_flow.hpp>
#include <boost/graph/edmonds_karp_max_flow.hpp>
#include <boost/graph/push_relabel_max_flow.hpp>

#include "graph.hpp"
#include "util.hpp"

namespace max_flow {

typedef boost::adjacency_list<
    boost::vecS, boost::vecS, boost::directedS, boost::no_property,
    boost::property<
        boost::edge_capacity_t, int,
        boost::property<
            boost::edge_residual_capacity_t, int,
            boost::property<boost::edge_reverse_t,
                            boost::adjacency_list<>::edge_descriptor>>>>
    Graph;
typedef boost::property_map<Graph, boost::edge_capacity_t>::type
    EdgeCapacityMap;
typedef boost::property_map<Graph, boost::edge_residual_capacity_t>::type
    ResidualCapacityMap;
typedef boost::property_map<Graph, boost::edge_reverse_t>::type ReverseEdgeMap;

inline int max_flow(int n, std::list<graph::Edge> links,
                    std::list<graph::Edge> edges, graph::Edge removed_link) {

  Graph g(n);
  EdgeCapacityMap capacity_map = get(boost::edge_capacity, g);
  ReverseEdgeMap reverse_map = get(boost::edge_reverse, g);
  ResidualCapacityMap residual_map = get(boost::edge_residual_capacity, g);

  for (auto &edge : edges) {
    auto e = add_edge(edge.first, edge.second, g).first;
    auto rev_e = add_edge(edge.second, edge.first, g).first;
    capacity_map[e] = (int)edge.weight;
    capacity_map[rev_e] = (int)edge.weight; // Reverse edge
    reverse_map[e] = rev_e;
    reverse_map[rev_e] = e;
  }
  for (auto &edge : links) {
    if (edge.first == removed_link.first &&
        edge.second == removed_link.second) {
      continue;
    }
    auto e = add_edge(edge.first, edge.second, g).first;
    auto rev_e = add_edge(edge.second, edge.first, g).first;
    capacity_map[e] = 1;
    capacity_map[rev_e] = 1; // Reverse edge
    reverse_map[e] = rev_e;
    reverse_map[rev_e] = e;
  }

  int maxFlow =
      boost::edmonds_karp_max_flow(g, removed_link.first, removed_link.second);
  //   boost::boykov_kolmogorov_max_flow(g, removed_link.first,
  //   removed_link.second);
  //   boost::push_relabel_max_flow(g, removed_link.first, removed_link.second);

  return maxFlow;
}

class ReusableMaxFlow {
  Graph g;
  EdgeCapacityMap capacity_map;
  ReverseEdgeMap reverse_map;
  ResidualCapacityMap residual_map;

public:
  void init(int n, std::list<graph::Edge> links, std::list<graph::Edge> edges,
            std::list<graph::Edge> removed_links,
            std::list<graph::Edge> additional_links) {
    g = Graph(n);
    capacity_map = get(boost::edge_capacity, g);
    reverse_map = get(boost::edge_reverse, g);
    residual_map = get(boost::edge_residual_capacity, g);

    for (auto &edge : edges) {
      auto e = add_edge(edge.first, edge.second, g).first;
      auto rev_e = add_edge(edge.second, edge.first, g).first;
      capacity_map[e] = (int)edge.weight;
      capacity_map[rev_e] = (int)edge.weight; // Reverse edge
      reverse_map[e] = rev_e;
      reverse_map[rev_e] = e;
    }
    for (auto &edge : links) {
      if (std::any_of(removed_links.begin(), removed_links.end(),
                      [&edge](auto &removed_link) {
                        return (edge.first == removed_link.first &&
                                edge.second == removed_link.second) ||
                               (edge.second == removed_link.first &&
                                edge.first == removed_link.second);
                      })) {
        continue;
      }
      auto e = add_edge(edge.first, edge.second, g).first;
      auto rev_e = add_edge(edge.second, edge.first, g).first;
      capacity_map[e] = 1;
      capacity_map[rev_e] = 1; // Reverse edge
      reverse_map[e] = rev_e;
      reverse_map[rev_e] = e;
    }
    for (auto &edge : additional_links) {
      auto e = add_edge(edge.first, edge.second, g).first;
      auto rev_e = add_edge(edge.second, edge.first, g).first;
      capacity_map[e] = 1;
      capacity_map[rev_e] = 1; // Reverse edge
      reverse_map[e] = rev_e;
      reverse_map[rev_e] = e;
    }
  }

  int max_flow(int s, int t) {
    int flow = boost::edmonds_karp_max_flow(g, s, t);
    // "undo" flow
    // augment(t, s, flow);
    return flow;
  }

  void augment(int s, int t, int value) {
    boost::graph_traits<Graph>::edge_descriptor edge;
    bool exists;
    boost::tie(edge, exists) = boost::edge(s, t, g);
    if (exists) {
      capacity_map[edge] += value;
    } else {
      auto e = add_edge(s, t, g).first;
      auto rev_e = add_edge(t, s, g).first;
      capacity_map[e] = value;
      capacity_map[rev_e] = 0; // Reverse edge
      reverse_map[e] = rev_e;
      reverse_map[rev_e] = e;
    }
  }
};

inline int max_flow(int n, int s, int t, std::list<graph::Edge> links,
                    std::list<graph::Edge> edges,
                    std::list<graph::Edge> removed_links,
                    std::list<graph::Edge> additional_links) {
  ReusableMaxFlow g;
  g.init(n, links, edges, removed_links, additional_links);

  int maxFlow = g.max_flow(s, t);

  return maxFlow;
}

} // namespace max_flow

#endif // MAX_FLOW_HPP
