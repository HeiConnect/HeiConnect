#include <algorithm>
#include <boost/pending/disjoint_sets.hpp>
#include <list>
#include <vector>

#include "graph.hpp"

namespace kruskal {

using graph::Edge;

class KruskalGraph {
  int n, m;
  std::vector<std::pair<Edge, double>> edges;

public:
  KruskalGraph(int _n) : n(_n), m(0) {}

  KruskalGraph(graph::Cactus cactus) : n(cactus.num_nodes()), m(0) {
    for (int u = 1; u <= cactus.num_nodes(); ++u) {
      for (int v = u + 1; v <= cactus.num_nodes(); ++v) {
        auto &link = cactus.links[u][v];
        if (link.first) {
          add_edge({u, v, link.weight}, link.weight);
        }
      }
    }
  }

  KruskalGraph(graph::Cactus cactus, graph::DynamicCactus &dyn_cactus)
      : n(cactus.num_nodes()), m(0) {
    for (int u = 1; u <= cactus.num_nodes(); ++u) {
      for (int v = u + 1; v <= cactus.num_nodes(); ++v) {
        auto &link = cactus.links[u][v];
        if (link.first) {
          int num_covered_cuts =
              dyn_cactus.num_increased_cuts(link.first, link.second);
          if (num_covered_cuts > 0) {
            add_edge({u, v, link.weight}, link.weight / num_covered_cuts);
          }
        }
      }
    }
  }

  void add_edge(Edge edge, double heuristic) {
    ++m;
    edges.push_back({edge, heuristic});
  }

  std::list<Edge> kruskal_mst() {
    DEBUG("Computing MST");
    std::list<Edge> minimum_spanning_tree;
    std::sort(
        edges.begin(), edges.end(),
        [](const std::pair<Edge, double> &a, const std::pair<Edge, double> &b) {
          return a.second < b.second;
        });

    boost::disjoint_sets_with_storage<> ds(n);

    for (int edge_idx = 0; edge_idx < edges.size(); ++edge_idx) {
      int src_root = ds.find_set(edges[edge_idx].first.first - 1);
      int dest_root = ds.find_set(edges[edge_idx].first.second - 1);

      if (src_root != dest_root) {
        DEBUG("Link in MST: " << edges[edge_idx].first.first << "-"
                              << edges[edge_idx].first.second << ": "
                              << edges[edge_idx].second << " / "
                              << edges[edge_idx].first.weight);
        minimum_spanning_tree.push_back(edges[edge_idx].first);
        ds.union_set(src_root, dest_root);
      }
    }

    return minimum_spanning_tree;
  }
};

} // namespace kruskal
