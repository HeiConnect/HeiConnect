#ifndef WATANABE_HPP
#define WATANABE_HPP

#include "graph.hpp"
#include <limits>
#include <stack>
#include <unordered_map>

namespace watanabe {

std::list<graph::Edge> smc(graph::DynamicCactus &g);
std::list<graph::Edge> fsm(graph::DynamicCactus &g);
std::list<graph::Edge> hbd(graph::DynamicCactus &g);

/**
 * Data structure
 */

class WGraph {
public:
  std::vector<std::vector<double>> c;
  std::vector<std::vector<double>> d;
  std::vector<std::vector<std::pair<int, int>>> b;
  std::vector<std::list<int>> adj_list;
  std::vector<std::vector<bool>> in_subgraph;
  std::unordered_map<int, int> cactus_to_w;
  std::unordered_map<int, int> w_to_cactus;
  int n, cycles, N;

  WGraph(int _n, int _cycles) : n(_n), cycles(_cycles) {
    N = _n + _cycles;
    c.resize(n + cycles,
             std::vector(n + cycles, std::numeric_limits<double>::infinity()));
    d.resize(n + cycles,
             std::vector(n + cycles, std::numeric_limits<double>::infinity()));
    b.resize(n + cycles, std::vector(n + cycles, std::pair<int, int>{0, 0}));
    in_subgraph.resize(n + cycles, std::vector(n + cycles, false));
    adj_list.resize(n + cycles, {});
  }

  void add_tree_edge(int u, int v) {
    adj_list[u].push_back(v);
    adj_list[v].push_back(u);
  }

  void add_edge_cost(int u, int v, double cost) {
    u = cactus_to_w[u];
    v = cactus_to_w[v];
    if (cost < c[cycles + u - 1][cycles + v - 1]) {
      c[cycles + u - 1][cycles + v - 1] = cost;
      c[cycles + v - 1][cycles + u - 1] = cost;
    }
  }

  void compute_dist() {
    // DIST 1.
    std::vector<std::vector<std::pair<int, int>>> a;
    a.resize(N, std::vector<std::pair<int, int>>(N, {-1, 0}));
    for (int u = 0; u < N; ++u) {
      a[u][u] = {0, u};
      d[u][u] = std::numeric_limits<double>::infinity();
      b[u][u] = {u, u};
      for (int v = u + 1; v < N; ++v) {
        auto path = dist(u, v);
        int distance = path.size() - 1;
        a[u][v] = {distance, *(++path.begin())};
        a[v][u] = {distance, *(++path.rbegin())};
        d[u][v] = c[u][v];
        d[v][u] = c[v][u];
        b[u][v] = {u, v};
        b[v][u] = {v, u};
      }
    }
    // DIST 2.
    std::list<graph::Edge> edge_list;
    for (int u = 0; u < N; ++u) {
      for (int v = u + 1; v < N; ++v) {
        if (c[u][v] != 0.) {
          edge_list.push_back({u, v, c[u][v]});
          edge_list.push_back({v, u, c[u][v]});
        }
      }
    }
    edge_list.sort([&a](graph::Edge &e1, graph::Edge &e2) {
      return a[e1.first][e1.second] > a[e2.first][e2.second];
    });
    for (auto &[u, v, weight] : edge_list) {
      int us = a[u][v].second;
      int vs = a[v][u].second;
      if (d[u][v] < d[u][vs]) {
        d[u][vs] = d[u][v];
        b[u][vs] = b[u][v];
        // d[us][u] = d[u][v];
        // b[us][u] = b[u][v];
      }
      if (d[u][v] < d[us][v]) {
        d[us][v] = d[u][v];
        b[us][v] = b[u][v];
        // d[v][vs] = d[u][v];
        // b[u][vs] = b[u][v];
      }
    }
  }

private:
  std::list<int> dist(int u, int v) {
    // DFS to find distance
    std::stack<int> stack;
    stack.push(u);
    std::vector<bool> visited;
    visited.resize(N, false);
    std::vector<int> parent;
    parent.resize(N, -1);
    parent[u] = u;
    while (!stack.empty()) {
      int current = stack.top();
      stack.pop();

      if (current == v) {
        // reconstruct and return path
        std::list<int> path = {v};
        int c = v;
        while (parent[c] != c) {
          path.push_front(parent[c]);
          c = parent[c];
        }
        return path;
      }

      if (!visited[current]) {
        visited[current] = true;
        for (auto n : adj_list[current]) {
          if (!visited[n]) {
            parent[n] = current;
            stack.push(n);
          }
        }
      }
    }
    throw std::invalid_argument("DFS cannot find path, something is wrong");
  }
};

WGraph remake(graph::DynamicCactus &g);

} // namespace watanabe

#endif
