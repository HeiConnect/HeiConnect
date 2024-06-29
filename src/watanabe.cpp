#include "watanabe.hpp"
#include "graph.hpp"
#include "lemon/core.h"
#include "lemon/smart_graph.h"
#include "util.hpp"
#include <algorithm>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/maximum_weighted_matching.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <cassert>
#include <functional>
#include <limits>
#include <queue>
#include <random>
#include <type_traits>
#include <unordered_map>
#include <vector>

#include <lemon/concepts/graph.h>
#include <lemon/concepts/maps.h>
#include <lemon/lgf_reader.h>
#include <lemon/matching.h>
#include <lemon/math.h>
#include <lemon/smart_graph.h>

namespace watanabe {

typedef lemon::concepts::Graph LGraph;

using Weight = // boost::multiprecision::cpp_dec_float_50;
    boost::multiprecision::number<boost::multiprecision::cpp_dec_float<50>,
                                  boost::multiprecision::et_off>;

typedef boost::property<boost::edge_weight_t, Weight,
                        boost::property<boost::edge_index_t, int>>
    EdgeProperty;
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
                              boost::no_property, EdgeProperty>
    BoostGraph;

std::list<graph::Edge> smc_iteration(graph::DynamicCactus &gp) {
  WGraph g = remake(gp);
  g.compute_dist();

  std::list<graph::Edge> result;
  for (int v = 0; v < g.N; ++v) {
    if (g.adj_list[v].size() > 1) {
      // Not a leaf
      continue;
    }
    graph::Edge best = {0, 0, std::numeric_limits<double>::infinity()};
    for (int u = 0; u < g.N; ++u) {
      // assert(g.d[u][v] == g.d[v][u]);
      if (g.d[u][v] > 0. && g.d[u][v] < best.weight) {
        best = {g.b[u][v].first, g.b[u][v].second, g.d[u][v]};
      }
    }
    assert(best.first != 0 || best.second != 0);
    DEBUG((v - g.cycles + 1)
         << ": " << (best.first - g.cycles + 1) << "-"
         << (best.second - g.cycles + 1) << " with weight " << best.weight);
    result.push_back({g.w_to_cactus[best.first - g.cycles + 1],
                      g.w_to_cactus[best.second - g.cycles + 1], best.weight});
  }

  return result;
}

std::list<graph::Edge> dijkstra(WGraph &g, int src, int target) {
  std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>,
                      std::greater<std::pair<int, int>>>
      pq;
  std::vector<double> dist(g.N, std::numeric_limits<double>::infinity());
  std::vector<std::pair<int, graph::Edge>> parent(
      g.N, {-1, {-1, -1, std::numeric_limits<double>::infinity()}});

  pq.push({0, src});
  dist[src] = 0;
  parent[src] = {src, {0, 0, 0.}};

  while (!pq.empty()) {
    int u = pq.top().second;
    pq.pop();

    for (int v = 0; v < g.d[u].size(); ++v) {
      double weight = g.d[u][v];
      if (dist[u] != std::numeric_limits<double>::infinity() &&
          dist[u] + weight < dist[v]) {
        dist[v] = dist[u] + weight;
        parent[v] = {u, {u, v, weight}};
        pq.push({dist[v], v});
      }
    }
  }

  std::list<graph::Edge> path;
  std::pair<int, graph::Edge> current = parent[target];
  do {
    path.push_back(current.second);
    current = parent[current.first];
  } while (current.first != parent[current.first].first);
  return path;
}

std::vector<std::list<graph::Edge>> dijkstra(WGraph &g, int src) {
  std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>,
                      std::greater<std::pair<int, int>>>
      pq;
  std::vector<double> dist(g.N, std::numeric_limits<double>::infinity());
  std::vector<std::pair<int, graph::Edge>> parent(
      g.N, {-1, {-1, -1, std::numeric_limits<double>::infinity()}});

  pq.push({0, src});
  dist[src] = 0;
  parent[src] = {src, {0, 0, 0.}};

  while (!pq.empty()) {
    int u = pq.top().second;
    pq.pop();

    for (int v = 0; v < g.d[u].size(); ++v) {
      double weight = g.d[u][v];
      if (dist[u] != std::numeric_limits<double>::infinity() &&
          dist[u] + weight < dist[v]) {
        dist[v] = dist[u] + weight;
        parent[v] = {u, {u, v, weight}};
        pq.push({dist[v], v});
      }
    }
  }

  std::vector<std::list<graph::Edge>> paths;
  paths.resize(g.N, {});
  for (int target = 0; target < g.N; ++target) {
    if (target == src)
      continue;
    std::pair<int, graph::Edge> current = parent[target];
    do {
      paths[target].push_back(current.second);
      current = parent[current.first];
    } while (current.first != parent[current.first].first);
  }
  return paths;
}

std::list<graph::Edge> fsm_iteration(graph::DynamicCactus &gp) {
  WGraph g = remake(gp);
  g.compute_dist();

  std::vector<int> leafes;
  for (int v = 0; v < g.N; ++v) {
    if (g.adj_list[v].size() <= 1) {
      leafes.push_back(v);
    }
  }

  // compute d'' using shortest path distances
  std::vector<std::vector<double>> dpp(g.N);
  for (auto &v : dpp) {
    v.resize(g.N, std::numeric_limits<double>::infinity());
  }
  std::vector<std::vector<std::list<graph::Edge>>> Ep(g.N);
  for (auto &v : Ep) {
    v.resize(g.N, {});
  }

  double max_dpp = 0.;
  for (int ui = 0; ui < leafes.size(); ++ui) {
    int u = leafes[ui];
    std::vector<std::list<graph::Edge>> paths = dijkstra(g, u);
    for (int vi = ui + 1; vi < leafes.size(); ++vi) {
      int v = leafes[vi];
      double sum = 0.;
      auto &path = paths[v];
      for (auto &e : path) {
        sum += e.weight;
      }
      dpp[u][v] = sum;
      dpp[v][u] = sum;
      for (auto &e : path) {
        Ep[u][v].push_back(e);
        Ep[v][u].push_back(e);
      }

      if (g.d[u][v] != std::numeric_limits<double>::infinity() &&
          max_dpp < g.d[u][v]) {
        max_dpp = g.d[u][v];
      }
    }
  }
  max_dpp += 1.;

  std::list<graph::Edge> result;
  // Create edges with altered cost function
  std::list<graph::Edge> edges;
  for (int i = 0; i < leafes.size(); ++i) {
    int u = leafes[i];
    for (int j = i + 1; j < leafes.size(); ++j) {
      int v = leafes[j];
      if (g.d[u][v] > 0. &&
          g.d[u][v] != std::numeric_limits<double>::infinity()) {
        edges.push_back({u, v, max_dpp - dpp[u][v]});
      }
    }
  }

  // construct boost graph to find matching among edges

  // new vertex ids
  std::unordered_map<int, int> map_to;
  std::unordered_map<int, int> map_from;
  int num_boost_nodes = 0;
  for (auto &e : edges) {
    if (map_to.find(e.first) == map_to.end()) {
      map_to[e.first] = num_boost_nodes;
      map_from[num_boost_nodes] = e.first;
      ++num_boost_nodes;
    }
    if (map_to.find(e.second) == map_to.end()) {
      map_to[e.second] = num_boost_nodes;
      map_from[num_boost_nodes] = e.second;
      ++num_boost_nodes;
    }
  }

  // deduplicate edges
  std::vector<std::vector<double>> weights(
      num_boost_nodes, std::vector<double>(num_boost_nodes, 0.));
  for (auto &e : edges) {
    if (e.weight > weights[map_to[e.first]][map_to[e.second]]) {
      weights[map_to[e.first]][map_to[e.second]] = e.weight;
      weights[map_to[e.second]][map_to[e.first]] = e.weight;
    }
  }

  // LEMON Maximum Weighted Matching
  lemon::SmartGraph lemon_graph;
  lemon::SmartGraph::EdgeMap<int> lemon_weight(lemon_graph);
  std::vector<lemon::SmartGraph::Node> lemon_nodes;

  for (int i = 0; i < num_boost_nodes; ++i) {
    lemon_nodes.push_back(lemon_graph.addNode());
  }

  // add edges
  for (int u = 0; u < weights.size(); ++u) {
    for (int v = u + 1; v < weights.size(); ++v) {
      if (weights[u][v] != 0.) {
        assert(weights[u][v] != std::numeric_limits<double>::infinity());
        auto lemon_edge = lemon_graph.addEdge(lemon_nodes[u], lemon_nodes[v]);
        lemon_weight.set(lemon_edge, (int)(100000 * weights[u][v]));
      }
    }
  }

  lemon::MaxWeightedMatching<lemon::SmartGraph> matching(lemon_graph,
                                                         lemon_weight);
  matching.run();
  for (lemon::SmartGraph::Node &node : lemon_nodes) {
    auto mate = matching.mate(node);
    if (mate != lemon::INVALID) {
      int u = lemon_graph.id(node);
      int v = lemon_graph.id(mate);
      if (u < v) {
        int uw = map_from[u];
        int vw = map_from[v];
        for (auto &[u, v, weight] : Ep[uw][vw]) {
          result.push_back({g.w_to_cactus[g.b[u][v].first - g.cycles + 1],
                            g.w_to_cactus[g.b[u][v].second - g.cycles + 1],
                            g.d[u][v]});
        }
      }
    }
  }

  return result;
}

std::list<graph::Edge> hbd_iteration(graph::DynamicCactus &gp) {
  WGraph g = remake(gp);
  g.compute_dist();

  std::vector<int> leafes;
  for (int v = 0; v < g.N; ++v) {
    if (g.adj_list[v].size() <= 1) {
      leafes.push_back(v);
    }
  }

  // compute d' using shortest path distances
  std::vector<std::vector<double>> dp(g.N);
  for (auto &v : dp) {
    v.resize(g.N, std::numeric_limits<double>::infinity());
  }
  std::vector<std::vector<std::list<graph::Edge>>> Ep(g.N);
  for (auto &v : Ep) {
    v.resize(g.N, {});
  }

  double max_dp = 0.;
  for (int ui = 0; ui < leafes.size(); ++ui) {
    int u = leafes[ui];
    std::vector<std::list<graph::Edge>> paths = dijkstra(g, u);
    for (int vi = ui + 1; vi < leafes.size(); ++vi) {
      int v = leafes[vi];
      double sum = 0.;
      auto &path = paths[v];
      for (auto &e : path) {
        sum += e.weight;
      }
      dp[u][v] = sum;
      dp[v][u] = sum;
      for (auto &e : path) {
        Ep[u][v].push_back(e);
        Ep[v][u].push_back(e);
      }

      if (g.d[u][v] != std::numeric_limits<double>::infinity() &&
          max_dp < g.d[u][v]) {
        max_dp = g.d[u][v];
      }
    }
  }
  max_dp += 1.;

  std::list<graph::Edge> result;
  // Create edges with altered cost function
  std::list<graph::Edge> edges;
  std::vector<double> maximum_weights(g.N, 0.);
  for (int i = 0; i < leafes.size(); ++i) {
    int u = leafes[i];
    for (int j = i + 1; j < leafes.size(); ++j) {
      int v = leafes[j];
      if (g.d[u][v] > 0. &&
          g.d[u][v] != std::numeric_limits<double>::infinity()) {
        edges.push_back({u, v, max_dp - dp[u][v]});
        if (maximum_weights[u] < max_dp - dp[u][v]) {
          maximum_weights[u] = max_dp - dp[u][v];
        }
        if (maximum_weights[v] < max_dp - dp[u][v]) {
          maximum_weights[v] = max_dp - dp[u][v];
        }
      }
    }
  }
  edges.remove_if([&maximum_weights](graph::Edge &e) {
    return e.weight < maximum_weights[e.first] &&
           e.weight < maximum_weights[e.second];
  });

  // construct boost graph to find matching among edges

  // new vertex ids
  std::unordered_map<int, int> map_to;
  std::unordered_map<int, int> map_from;
  int num_boost_nodes = 0;
  for (auto &e : edges) {
    if (map_to.find(e.first) == map_to.end()) {
      map_to[e.first] = num_boost_nodes;
      map_from[num_boost_nodes] = e.first;
      ++num_boost_nodes;
    }
    if (map_to.find(e.second) == map_to.end()) {
      map_to[e.second] = num_boost_nodes;
      map_from[num_boost_nodes] = e.second;
      ++num_boost_nodes;
    }
  }
  // deduplicate edges
  std::vector<std::vector<double>> weights(
      num_boost_nodes, std::vector<double>(num_boost_nodes, 0.));
  for (auto &e : edges) {
    if (e.weight > weights[map_to[e.first]][map_to[e.second]]) {
      weights[map_to[e.first]][map_to[e.second]] = e.weight;
      weights[map_to[e.second]][map_to[e.first]] = e.weight;
    }
  }

  // LEMON Maximum Weighted Matching
  lemon::SmartGraph lemon_graph;
  lemon::SmartGraph::EdgeMap<int> lemon_weight(lemon_graph);
  std::vector<lemon::SmartGraph::Node> lemon_nodes;

  for (int i = 0; i < num_boost_nodes; ++i) {
    lemon_nodes.push_back(lemon_graph.addNode());
  }

  // add edges
  for (int u = 0; u < weights.size(); ++u) {
    for (int v = u + 1; v < weights.size(); ++v) {
      if (weights[u][v] != 0.) {
        assert(weights[u][v] != std::numeric_limits<double>::infinity());
        auto lemon_edge = lemon_graph.addEdge(lemon_nodes[u], lemon_nodes[v]);
        lemon_weight.set(lemon_edge, (int)(100000 * weights[u][v]));
      }
    }
  }

  lemon::MaxWeightedMatching<lemon::SmartGraph> matching(lemon_graph,
                                                         lemon_weight);
  matching.run();
  for (lemon::SmartGraph::Node &node : lemon_nodes) {
    auto mate = matching.mate(node);
    if (mate != lemon::INVALID) {
      int u = lemon_graph.id(node);
      int v = lemon_graph.id(mate);
      if (u < v) {
        int uw = map_from[u];
        int vw = map_from[v];
        for (auto &[u, v, weight] : Ep[uw][vw]) {
          result.push_back({g.w_to_cactus[g.b[u][v].first - g.cycles + 1],
                            g.w_to_cactus[g.b[u][v].second - g.cycles + 1],
                            g.d[u][v]});
        }
      }
    }
  }

  return result;
}

std::list<graph::Edge>
algorithm_frame(graph::DynamicCactus &gp,
                std::function<std::list<graph::Edge>(graph::DynamicCactus &)>
                    algorithm_iteration) {
  std::list<graph::Edge> result;
  std::set<std::pair<int, int>> duplicate_check;
  while (gp.size()) {
    std::list<graph::Edge> partial_augmentation = algorithm_iteration(gp);

    // remove duplicates
    partial_augmentation.remove_if([&duplicate_check](graph::Edge &e) {
      std::pair<int, int> pair = {e.first, e.second};
      if (pair.first > pair.second) {
        std::swap(pair.first, pair.second);
      }
      if (duplicate_check.find(pair) == duplicate_check.end()) {
        duplicate_check.insert(pair);
        return false;
      }
      return true;
    });

    for (auto &e : partial_augmentation) {
      int u = gp.pi(e.first);
      int v = gp.pi(e.second);
      if (u != v) {
        gp.contract(u, v);
        result.push_back(e);
      } else {
        DEBUG("Unnecessary link: " << e.first << "-" << e.second << " (" << u
                                   << "-" << v << ")");
      }
    }
  }
  return result;
}

std::list<graph::Edge> smc(graph::DynamicCactus &g) {
  return algorithm_frame(g, smc_iteration);
}

std::list<graph::Edge> fsm(graph::DynamicCactus &g) {
  return algorithm_frame(g, fsm_iteration);
}

std::list<graph::Edge> hbd(graph::DynamicCactus &g) {
  return algorithm_frame(g, hbd_iteration);
}

WGraph remake(graph::DynamicCactus &g) {
  WGraph wg(g.num_nodes(), g.size());
  int cnt = 1;
  for (auto &cycle : g.node_list) {
    for (auto &nd : cycle) {
      if (wg.cactus_to_w.find(nd) == wg.cactus_to_w.end()) {
        wg.cactus_to_w[nd] = cnt;
        wg.w_to_cactus[cnt] = nd;
        ++cnt;
      }
    }
  }
  for (int i = 0; i < g.node_list.size(); ++i) {
    std::list<CactusNodeID> &cycle = g.node_list[i];
    for (CactusNodeID nd : cycle) {
      wg.add_tree_edge(i, wg.cycles + wg.cactus_to_w[nd] - 1);
    }
    if (cycle.size() >= 2) {
      wg.add_edge_cost(cycle.front(), cycle.back(),
                       std::numeric_limits<double>::infinity());
      if (cycle.size() >= 3) {
        auto it = cycle.begin();
        int cur = *it;
        ++it;
        while (it != cycle.end()) {
          wg.add_edge_cost(cur, *it, std::numeric_limits<double>::infinity());
          ++it;
        }
      }
    }
    for (auto &link : g.links) {
      if (g.pi(link.first) != g.pi(link.second)) {
        // FIXME: output graph will be augmented, but does not have the correct
        // links (instead pi(link))
        wg.add_edge_cost(g.pi(link.first), g.pi(link.second), link.weight);
      }
    }
  }

  return wg;
}

} // namespace watanabe
