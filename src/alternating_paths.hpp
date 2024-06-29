#ifndef ALTERNATING_PATHS_HPP
#define ALTERNATING_PATHS_HPP

#include "graph.hpp"
#include "util.hpp"
#include <boost/functional/hash.hpp>
#include <cassert>
#include <stack>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace alternating_paths {

struct AlternatingEdge {
  bool active;
  int source;
  int target;
  double weight;
  AlternatingEdge *other;

  double potential() { return active ? -weight : weight; }
};

class AlternatingGraph {
  int n;
  int depth_limit = 6;
  std::vector<std::list<AlternatingEdge>> adj_list;
  // cache invalid paths
  // stores the number of iterations it should not be checked
  std::unordered_set<size_t> path_cache;

public:
  AlternatingGraph(int _n, int _depth_limit = 5)
      : n(_n), depth_limit(_depth_limit) {
    adj_list.resize(n + 1, {});
  }
  void add_edge(graph::Edge edge, bool active) {
    adj_list[edge.second].push_back(
        {active, edge.second, edge.first, edge.weight, nullptr});
    adj_list[edge.first].push_back({active, edge.first, edge.second,
                                    edge.weight,
                                    &adj_list[edge.second].back()});
    adj_list[edge.second].back().other = &adj_list[edge.first].back();
  }

  void add_solution(std::list<graph::Edge> edges) {
    for (auto &e : edges) {
      add_edge(e, true);
    }
  }

  void add_links(std::list<graph::Edge> edges) {
    for (auto &e : edges) {
      add_edge(e, false);
    }
  }

  double potential(std::list<AlternatingEdge *> path) {
    bool last = !path.front()->active;
    double potential = 0.;
    for (auto edge : path) {
      assert(last != edge->active);
      last = !last;
      potential += edge->potential();
    }
    return potential;
  }

  std::list<std::pair<std::list<AlternatingEdge *>, double>>
  find_alternating_path() {
    std::list<std::pair<std::list<AlternatingEdge *>, double>> result;
    for (int i = 1; i <= n; ++i) {
      find_alternating_path(i, result);
    }
    return result;
  }

  std::list<std::pair<std::list<AlternatingEdge *>, double>>
  find_alternating_path(int start_node) {
    std::list<std::pair<std::list<AlternatingEdge *>, double>> result;
    find_alternating_path(start_node, result);
    return result;
  }

  void find_alternating_path(
      int start_node,
      std::list<std::pair<std::list<AlternatingEdge *>, double>> &result) {
    std::vector<bool> visited;
    visited.resize(n + 1, false);
    visited[start_node] = true;
    bool start_with_active = adj_list[start_node].size() > 1;
    std::list<AlternatingEdge *> path;
    for (auto &e : adj_list[start_node]) {
      if (!e.active || start_with_active) {
        // actual dfs
        rec_find_alternating_path(e, path, 1, visited, result);
      }
    }
  }

  void rec_find_alternating_path(
      AlternatingEdge &edge, std::list<AlternatingEdge *> &path, int depth,
      std::vector<bool> &visited,
      std::list<std::pair<std::list<AlternatingEdge *>, double>> &paths) {
    visited[edge.target] = true;
    path.push_back(&edge);

    if (path.size() > 1) {
      if (!path.back()->active || adj_list[path.back()->target].size() > 1) {
        // check if path decreases cost
        double pot = potential(path);
        if (pot < 0) {
          // check cache if path should be skipped
          size_t h = hash_path(path);
          if (path_cache.find(h) == path_cache.end() &&
              path.front()->source < path.back()->target) { // symmetry
            std::list<AlternatingEdge *> path_copy;
            for (const auto &e : path) {
              path_copy.push_back(e);
            }
            paths.push_back({path_copy, pot});
          }
        }
      }
    }

    if (depth < depth_limit) {
      for (auto &e : adj_list[edge.target]) {
        if (!visited[e.target] && edge.active != e.active) {
          rec_find_alternating_path(e, path, depth + 1, visited, paths);
        }
      }
    }

    visited[path.back()->target] = false;
    path.pop_back();
  }

  void apply_path(std::list<AlternatingEdge *> &path) {
    for (auto &e : path) {
      e->active = !e->active;
      e->other->active = !e->other->active;
    }
  }

  size_t hash_path(const std::list<AlternatingEdge *> &path) {
    std::vector<int> nodes(path.size() + 1);
    if (path.front()->source < path.back()->target) {
      nodes.push_back(path.front()->source);
      for (auto &e : path) {
        nodes.push_back(e->target);
      }
    } else {
      auto it = path.rbegin();
      while (it != path.rend()) {
        nodes.push_back((*it)->target);
        ++it;
      }
      nodes.push_back(path.front()->source);
    }
    // hash node sequence
    size_t seed = 0;
    boost::hash_range(seed, nodes.begin(), nodes.end());
    return seed;
  }

  void invalidate_path(const std::list<AlternatingEdge *> &path) {
    size_t h = hash_path(path);
    path_cache.insert(h);
  }
};

} // namespace alternating_paths

#endif // ALTERNATING_PATHS_HPP