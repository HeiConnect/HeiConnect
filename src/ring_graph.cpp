#include "ring_graph.hpp"
#include <iterator>
#include <stdexcept>

namespace ring_graph {

/**
 * RingGraph
 */

RingGraph::RingGraph(graph::Cactus cactus) {
  std::list<int> path;
  cactus.eulerian_path(1, path);
  std::vector<int> first_appearance;
  first_appearance.resize(cactus.num_nodes() + 1, -1);

  // close path in the end instead of adding new node
  path.pop_back();

  for (int nd : path) {
    int new_id = mapping.size();
    mapping.push_back(nd); // inserts new implicit node
    adj_lists.push_back({});
    if (first_appearance[nd] >= 0) {
      // insert new link of weight zero
      add_link(first_appearance[nd], new_id, 0.);
    } else {
      first_appearance[nd] = new_id;
    }
  }
  // set new number of nodes
  n = mapping.size();

  // copy previous links
  for (int u = 1; u <= cactus.num_nodes(); ++u) {
    for (int v = u + 1; v <= cactus.num_nodes(); ++v) {
      if (cactus.links[u][v].first) {
        add_link(first_appearance[u], first_appearance[v],
                 cactus.links[u][v].weight, cactus.links[u][v].first,
                 cactus.links[u][v].second);
      }
    }
  }
}

void RingGraph::add_link(int u, int v, double weight, int original_u,
                         int original_v) {
  links.push_back({u, v, weight, original_u, original_v});
  adj_lists[u].push_back(&links.back());
  adj_lists[v].push_back(&links.back());
  ++num_links;
}

/**
 * LinkIntersectionGraph
 */

LinkIntersectionGraph::LinkIntersectionGraph(std::vector<Link> &K) {
  n = K.size();
  link_of_node.resize(n);
  adj_lists.resize(n, {});
  int num_nodes = 0;
  for (int i = 0; i < K.size(); ++i) {
    link_of_node[i] = &K[i];
    if (K[i].u > num_nodes) {
      num_nodes = K[i].u;
    }
    if (K[i].v > num_nodes) {
      num_nodes = K[i].v;
    }
    for (int j = i + 1; j < K.size(); ++j) {
      if (K[i].crossing(K[j])) {
        add_link(i, j);
      }
    }
  }
  adjacent_links.resize(num_nodes + 1, {});
  link_of_node.resize(num_nodes + 1, nullptr);
  for (int i = 0; i < K.size(); ++i) {
    adjacent_links[K[i].u].insert(i);
    adjacent_links[K[i].v].insert(i);
  }
}

LinkIntersectionGraph::LinkIntersectionGraph(std::vector<Link *> &K) {
  n = K.size();
  link_of_node.resize(n);
  adj_lists.resize(n, {});
  int num_nodes = 0;
  for (int i = 0; i < K.size(); ++i) {
    link_of_node[i] = K[i];
    if (K[i]->u > num_nodes) {
      num_nodes = K[i]->u;
    }
    if (K[i]->v > num_nodes) {
      num_nodes = K[i]->v;
    }
    for (int j = i + 1; j < K.size(); ++j) {
      if (K[i]->crossing(*K[j])) {
        add_link(i, j);
      }
    }
  }
  adjacent_links.resize(num_nodes + 1, {});
  link_of_node.resize(num_nodes + 1, nullptr);
  for (int i = 0; i < K.size(); ++i) {
    adjacent_links[K[i]->u].insert(i);
    adjacent_links[K[i]->v].insert(i);
  }
}

bool LinkIntersectionGraph::is_connected_to_v_good(int v, RingInterval &v_bad) {
  // not connected if link intersection graph is empty or there is no link
  // adjacent to v
  if (n == 0 || adjacent_links.size() <= v) {
    return false;
  }
  // bfs to find a vertex not in v_bad
  std::list<int> queue;
  // start points: adjacent links to v
  std::copy(adjacent_links[v].begin(), adjacent_links[v].end(),
            std::back_inserter(queue));
  std::vector<bool> visited;
  visited.resize(n);
  for (auto &x : queue) {
    visited[x] = true;
    // target nodes: adjacent links to v-good nodes
    if (!v_bad.in_range(link_of_node[x]->u) ||
        !v_bad.in_range(link_of_node[x]->v)) {
      // found a v-good node
      return true;
    }
  }
  while (!queue.empty()) {
    int current = queue.front();
    queue.pop_front();
    for (int next : adj_lists[current]) {
      if (!visited[next]) {
        visited[next] = true;
        queue.push_back(next);
        // target nodes: adjacent links to v-good nodes
        if (!v_bad.in_range(link_of_node[next]->u) ||
            !v_bad.in_range(link_of_node[next]->v)) {
          // found a v-good node
          return true;
        }
      }
    }
  }
  // no v-good node found during bfs
  return false;
}

std::pair<std::vector<int>, int> LinkIntersectionGraph::connected_components() {
  if (!n) {
    // WARN("Link intersection graph is empty");
    return {{}, 0};
  }
  // BFS to find connected components
  int size = 0;
  std::vector<int> partition(n);
  std::vector<bool> visited;
  visited.resize(n, false);
  for (int i = 0; i < n; ++i) {
    if (visited[i])
      continue;
    visited[i] = true;
    std::list<int> queue = {i};

    while (!queue.empty()) {
      int nd = queue.front();
      queue.pop_front();
      partition[nd] = size;
      for (int next : adj_lists[nd]) {
        if (!visited[next]) {
          visited[next] = true;
          queue.push_back(next);
        }
      }
    }
    ++size;
  }

  return {partition, size};
}

/**
 * Solution
 */

RingInterval Solution::get_v_bad(int v) {
  assert(shortening_adj.size() > 0);
  int min = v;
  int max = v;
  // dfs from v, remembering borders
  std::stack<int> queue;
  queue.push(v);
  while (!queue.empty()) {
    int next = queue.top();
    queue.pop();
    for (auto &nd : shortening_adj[next]) {
      queue.push(nd);
      if (nd < min) {
        min = nd;
      } else if (nd > max) {
        max = nd;
      }
    }
  }
  return {min, max};
}

void Solution::compute_shortening_adj(int n) {
  shortening_adj.resize(n, {});
  for (auto &l : shortening) {
    shortening_adj[l.u].push_back(l.v);
  }
}

int Solution::lca(int u, int v) {
  std::list<int> ru_path = dfs(0, u);
  std::list<int> rv_path = dfs(0, v);
  int current_lca = 0;
  while (!ru_path.empty() && !rv_path.empty() &&
         ru_path.front() == rv_path.front()) {
    current_lca = ru_path.front();
    ru_path.pop_front();
    rv_path.pop_front();
  }
  return current_lca;
}

std::list<int> Solution::dfs(int source, int target) {
  std::stack<int> stack;
  stack.push(source);
  std::vector<bool> visited;
  visited.resize(directed_solution.size() + 1, false);
  std::vector<int> parent;
  parent.resize(directed_solution.size() + 1, -1);
  parent[source] = source;
  while (!stack.empty()) {
    int current = stack.top();
    stack.pop();

    if (current == target) {
      // reconstruct and return path
      std::list<int> path = {target};
      int c = target;
      while (parent[c] != c) {
        path.push_front(parent[c]);
        c = parent[c];
      }
      return path;
    }

    if (!visited[current]) {
      visited[current] = true;
      for (auto n : shortening_adj[current]) {
        if (!visited[n]) {
          parent[n] = current;
          stack.push(n);
        }
      }
    }
  }
  throw std::invalid_argument("DFS cannot find path, something is wrong");
}

} // namespace ring_graph
