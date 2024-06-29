#include "graph.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <numeric>
#include <omp.h>
#include <random>
#include <stdexcept>
#include <string>
#include <unordered_set>

#include "pugixml.hpp"
#include "util.hpp"

namespace graph {

/* GraphPair */

void GraphPair::read_graph(const std::filesystem::path &graph_file,
                           const std::filesystem::path &cactus_file) {
  original_graph.read_from_file(graph_file);
  cactus.read_from_file(cactus_file, *this);
}

int GraphPair::cactus_id(int original_graph_id) {
  return map_to_cactus[original_graph_id];
}
std::list<int> GraphPair::original_graph_ids(int cactus_id) {
  return map_to_original_graph[cactus_id];
}

/* Graph */

unsigned int Graph::num_nodes() const { return n; }
unsigned int Graph::num_edges() const { return m; }

bool Graph::is_edge(int u, int v) {
  for (auto target : adj_list[u]) {
    if (target.first == v) {
      return true;
    }
  }
  return false;
}

void Graph::init(int _n, int _m) {
  n = _n;
  m = _m;
  if (m == n * (n - 1) / 2) {
    WARN("The graph is complete");
  }
  adj_list.resize(n + 1, {});
}

void Graph::finalize() {
  // sort adjacency lists
  for (unsigned int i = 1; i < adj_list.size(); ++i) {
    adj_list[i].sort([](std::pair<int, double> a, std::pair<int, double> b) {
      return a.first < b.first;
    });
  }
}

void Graph::add_edge(int u, int v, double weight) {
  assert(u > 0);
  assert(v > 0);
  adj_list[u].push_back({v, weight});
  adj_list[v].push_back({u, weight});
}

void Graph::add_directed_edge(int u, int v, double weight) {
  adj_list[u].push_back({v, weight});
}

const std::list<std::pair<int, double>> &Graph::get_adj_list(int node) {
  return adj_list[node];
}

/* OriginalGraph */

void OriginalGraph::read_from_file(const std::filesystem::path &path) {
  if (!std::filesystem::exists(path)) {
    throw std::runtime_error("File not found: " + path.string());
  }
  std::ifstream file{path};
  std::string line;
  do {
    std::getline(file, line);
  } while (line[0] == '%');
  std::istringstream iss(line);
  iss >> n >> m;
  DEBUG("Parsing original graph with " << n << " nodes and " << m << " edges");
  init(n, m);
  unsigned int u = 1, v;
  while (std::getline(file, line)) {
    if (line[0] == '%') {
      continue;
    }
    std::istringstream iss(line);
    while ((iss >> v)) {
      if (u < v)
        add_edge(u, v);
    }

    ++u;
  }
  finalize();
}

/* Cactus */

void Cactus::init(int _n, int _m) {
  Graph::init(_n, _m);
  links.resize(n + 1, {});
  for (int i = 0; i <= n; ++i) {
    links[i].resize(n + 1, {0, 0, -1.});
  }
}

void Cactus::read_from_file(const std::filesystem::path &path,
                            GraphPair &graph) {
  if (!std::filesystem::exists(path)) {
    throw std::runtime_error("File not found: " + path.string());
  }
  graph.map_to_cactus.resize(graph.original_graph.num_nodes() + 1, 0);

  pugi::xml_document doc;
  pugi::xml_parse_result result = doc.load_file(path.c_str());
  if (!result) {
    throw std::runtime_error("Could not open file " + path.string());
  }
  pugi::xml_node xml_graph = doc.child("graphml").child("graph");
  n = 0;
  m = 0;
  for (auto &node : xml_graph.children("node")) {
    ++n;
  }
  for (auto &edge : xml_graph.children("edge")) {
    ++m;
  }
  graph.map_to_original_graph.resize(n + 1, {});
  DEBUG("Parsing cactus with " << n << " nodes and " << m << " edges");

  init(n, m);

  // parse contained vertices
  contained_nodes.resize(n, {});
  for (auto &node : xml_graph.children("node")) {
    int cactus_node_id = node.attribute("id").as_int();
    std::string nodes_string =
        node.find_child_by_attribute("key", "containedVertices")
            .text()
            .as_string();
    std::vector<int> nodes = split(nodes_string, ',');
    for (int nd : nodes) {
      contained_nodes[cactus_node_id].insert(nd);
      graph.map_to_original_graph[cactus_node_id + 1].push_back(nd + 1);
      graph.map_to_cactus[nd + 1] = cactus_node_id + 1;
    }
  }

  // parse edges
  for (auto &edge : xml_graph.children("edge")) {
    double weight =
        edge.find_child_by_attribute("key", "weight").text().as_double();
    unsigned int source = edge.attribute("source").as_int() + 1;
    unsigned int target = edge.attribute("target").as_int() + 1;
    add_edge(source, target, weight);
  }

  finalize();
  compute_min_cut();
}

double round_places(double value, int places) {
  double factor = std::pow(10., places);
  int tmp = value * factor;
  return tmp / factor;
}

void GraphPair::generate_links(int seed, double fraction, int distribution) {
  DEBUG("Dist: " << distribution);
  std::mt19937 mt(seed);
  std::mt19937 mt_fraction(seed);
  // uniform
  std::uniform_real_distribution<double> dist{0., 1.};
  // normal
  std::normal_distribution<double> normal_dist(1., 0.5);
  // int (watanabe)
  std::uniform_int_distribution<int> w2(1, 2);
  std::uniform_int_distribution<int> w9(1, 9);
  std::uniform_int_distribution<int> w99(1, 99);
  std::uniform_int_distribution<int> w999(1, 999);
  // for link fraction
  std::uniform_real_distribution<> flip{0, 1};
  for (int i = 1; i <= cactus.num_nodes(); ++i) {
    for (int j = i + 1; j <= cactus.num_nodes(); ++j) {
      double weight;
      switch (distribution) {
      case 0:
        weight = 1.;
        break;
      case 1:
        weight = dist(mt);
        break;
      case 2:
        weight = std::abs(normal_dist(
            mt)); // FIXME: negative weights are invalid, regenerate instead?
        break;
      case 3:
        weight = w2(mt) / 2.;
        break;
      case 4:
        weight = w9(mt) / 9.;
        break;
      case 5:
        weight = w99(mt) / 99.;
        break;
      case 6:
        weight = w999(mt) / 999.;
        break;
      case 7:
        weight = dist(mt);
        weight = round_places(weight, 5);
        break;
      default:
        throw std::invalid_argument("Invalid distribution");
      }
      // double weight = random ? dist(mt) : 1.;
      for (auto &f : original_graph_ids(i)) {
        for (auto &s : original_graph_ids(j)) {
          if (!std::any_of(original_graph.get_adj_list(f).begin(),
                           original_graph.get_adj_list(f).end(),
                           [&](auto &el) { return el.first == s; })) {
            if (flip(mt_fraction) <= fraction) {
              add_link(i, j, f, s, weight);
            }
            goto out;
          }
        }
      }
    out:;
    }
  }
}

void GraphPair::load_links_from_file(const std::filesystem::path &link_file) {
  if (!std::filesystem::exists(link_file)) {
    throw std::runtime_error("Link file not found: " + link_file.string());
  }

  std::ifstream file{link_file};
  int u, v;
  double weight;
  while ((file >> u >> v >> weight)) {
    int cactus_u = cactus_id(u), cactus_v = cactus_id(v);
    if (cactus.links[cactus_u][cactus_v].first == 0 ||
        weight < cactus.links[cactus_u][cactus_v].weight) {
      add_link(cactus_u, cactus_v, u, v, weight);
    }
  }
}

void GraphPair::add_links(const std::filesystem::path &link_file,
                          double fraction, int distribution) {
  if (std::filesystem::exists(link_file)) {
    load_links_from_file(link_file);
  } else {
    WARN("Link file does not exists, generating uniform links with seed 0");
    add_links(0, fraction, distribution);
  }
}

void GraphPair::add_links(int seed, double fraction, int distribution) {
  generate_links(seed, fraction, distribution);
}

void GraphPair::add_link(int u, int v, int original_u, int original_v,
                         double weight) {
  cactus.links[u][v] = {original_u, original_v, weight};
  cactus.links[v][u] = {original_u, original_v, weight};
}

double Cactus::compute_min_cut() {
  if (min_cut >= 0) {
    return min_cut;
  }
  // Disconnected
  if (m == 0) {
    min_cut = 0.;
    return min_cut;
  }
  // Tree
  if (m == n - 1) {
    min_cut = adj_list[1].front().second;
    return min_cut;
  }
  // Cactus with at least one cycle
  double min = adj_list[1].front().second;
  for (int i = 1; i < adj_list.size(); ++i) {
    for (auto &pair : adj_list[i]) {
      if (pair.second < min) {
        min = pair.second;
      }
    }
  }
  min_cut = 2. * min;
  return min_cut;
}

std::list<CactusCut> GraphPair::get_min_cuts() {
  double min_cut = cactus.compute_min_cut();
  INFO("Min cut: " << min_cut);
  // Assume a cactus
  double threshold = min_cut - 0.5; // for double comparison
  if (cactus.num_edges() == 0) {
    throw std::runtime_error("Disconnected graph not supported");
  }
  std::list<CactusCut> cuts;
  std::set<std::pair<int, int>> processed;
  for (int u = 1; u <= cactus.num_nodes(); ++u) {
    for (const auto &edge : cactus.get_adj_list(u)) {
      int v = edge.first;
      double weight = edge.second;
      if (u > v) {
        {
          continue;
        }
      }
      if (processed.find({u, v}) == processed.end()) { // contains({u, v})) {
        processed.insert({u, v});
        Edge start_edge = {u, v, weight};
        if (start_edge.weight > threshold) { // tree edge
          cuts.push_back(CactusCut{*this, start_edge});
        } else { // cycle edge
          std::vector<std::pair<int, int>> cycle_edges;
          // BFS to find cycle edges
          std::vector<bool> visited;
          visited.resize(cactus.num_nodes() + 1, false);
          std::vector<int> parent;
          parent.resize(cactus.num_nodes() + 1);
          visited[start_edge.first] = true;
          parent[start_edge.first] = 0;
          std::list<int> queue;
          for (const auto &target : cactus.get_adj_list(start_edge.first)) {
            if (target.first != start_edge.second) {
              queue.push_back(target.first);
              parent[target.first] = start_edge.first;
            }
          }
          while (!queue.empty()) {
            int current = queue.front();
            queue.pop_front();
            visited[current] = true;
            for (auto &next : cactus.get_adj_list(current)) {
              if (!visited[next.first]) {
                parent[next.first] = current;
                if (next.first == start_edge.second) { // cycle found
                  queue.clear();
                  break;
                }
                queue.push_back(next.first);
              }
            }
          }
          // compute cycle edges
          cycle_edges.push_back({u, v});
          int current = start_edge.second;
          while (parent[current]) {
            int edge_u = current, edge_v = parent[current];
            if (edge_u > edge_v) {
              std::swap(edge_u, edge_v);
            }
            cycle_edges.push_back({edge_u, edge_v});
            processed.insert({edge_u, edge_v});
            current = parent[current];
          }
          // add all pairs of edges to cuts
          DEBUG("Adding cycle of size " << cycle_edges.size());
          // #pragma omp parallel
          {
            std::list<CactusCut> local_cuts;
            // #pragma omp for
            for (int i = 0; i < cycle_edges.size(); ++i) {
              for (int j = i + 1; j < cycle_edges.size(); ++j) {
                local_cuts.push_back(
                    CactusCut{*this, cycle_edges[i], cycle_edges[j]});
              }
            }
            // #pragma omp critical
            { cuts.splice(cuts.end(), local_cuts); }
          }
        }
      }
    }
  }

  INFO("Number of min cuts: " << cuts.size());
  return cuts;
}

void Cactus::expand_tree_edges() {
  double cut_threshold = compute_min_cut() - 0.5;
  for (int nd = 0; nd < num_nodes(); ++nd) {
    for (auto &e : adj_list[nd]) {
      if (e.second > cut_threshold) {
        e.second /= 2;
        m += 1;
        add_edge(nd, e.first, e.second);
      }
    }
  }
}

void Cactus::eulerian_path(int u, std::list<int> &circuit) {
  double cut_threshold = compute_min_cut() - 0.5;
  while (!adj_list[u].empty()) {
    auto [v, weight] = adj_list[u].back();
    adj_list[u].pop_back();
    // remove reverse edge if it is not a tree edge
    if (weight < cut_threshold)
      adj_list[v].remove_if(
          [u](std::pair<int, double> &val) { return val.first == u; });
    eulerian_path(v, circuit);
  }
  circuit.push_back(u);
}

/* CactusCut */

void CactusCut::compute_partition(GraphPair &g) {
  partition.resize(g.cactus.num_nodes() + 1, false);
  // do BFS to one side
  std::vector<bool> visited;
  visited.resize(g.cactus.num_nodes() + 1, false);
  visited[e1.first] = true;
  visited[e1.second] = true;
  std::list<int> queue;
  queue.push_back(e1.first);
  while (!queue.empty()) {
    int target = queue.front();
    queue.pop_front();
    partition[target] = true;
    visited[target] = true;
    for (auto &next : g.cactus.get_adj_list(target)) {
      if (!visited[next.first] && !is_same_edge(e1, target, next.first) &&
          (!e2.has_value() || !is_same_edge(*e2, target, next.first))) {
        queue.push_back(next.first);
      }
    }
  }
}

void OriginalGraph::write_to_file(const std::filesystem::path &path) {
  std::ofstream file;
  file.open(path);
  file << n << " " << m << std::endl;
  for (int i = 1; i <= n; ++i) {
    for (auto &target : adj_list[i]) {
      file << target.first << " ";
    }
    file << std::endl;
  }
  file.close();
}

void OriginalGraph::add_link(Edge &e) {
  if (std::any_of(adj_list[e.first].begin(), adj_list[e.first].end(),
                  [&](auto &pair) { return pair.first == e.second; })) {
    WARN("Tried to add double edge: " << e.first << "-" << e.second
                                      << " of weight " << e.weight);
  }
  ++m;
  add_edge(e.first, e.second);
  DEBUG("Adding " << e.first << "-" << e.second << " " << e.weight
                  << " to original graph");
}

std::vector<std::vector<int>> Graph::all_pair_distances() {
  std::vector<std::vector<int>> distances;
  distances.resize(n + 1, {});
  for (int i = 1; i <= n; ++i) {
    distances[i].resize(n + 1, std::numeric_limits<int>::max());
    distances[i][i] = 0;
  }
  for (int i = 1; i <= n; ++i) {
    // bfs from each node
    std::list<int> queue;
    std::list<int> next_queue;
    std::vector<bool> visited;
    visited.resize(n + 1, false);
    queue.push_back(i);
    visited[i] = true;
    int level = 1;
    do {
      while (!queue.empty()) {
        int current = queue.front();
        queue.pop_front();
        for (auto &target : adj_list[current]) {
          if (!visited[target.first]) {
            next_queue.push_back(target.first);
            visited[target.first] = true;
            distances[i][target.first] = level;
          }
        }
      }
      ++level;
      std::swap(queue, next_queue);
    } while (!queue.empty());
  }
  return distances;
}

// Dynamic cactus

void DynamicCactus::read_from_file(const std::filesystem::path &path) {
  if (!std::filesystem::exists(path)) {
    throw std::runtime_error("File not found: " + path.string());
  }

  pugi::xml_document doc;
  pugi::xml_parse_result result = doc.load_file(path.c_str());
  if (!result) {
    throw std::runtime_error("Could not open file " + path.string());
  }
  pugi::xml_node xml_graph = doc.child("graphml").child("graph");
  n = 0;
  m = 0;
  for (auto &node : xml_graph.children("node")) {
    ++n;
  }
  for (auto &edge : xml_graph.children("edge")) {
    ++m;
  }
  DEBUG("Parsing cactus with " << n << " nodes and " << m << " edges");

  init(n, m);
  node_to_cycle.resize(n + 1, {});

  // parse edges
  for (auto &edge : xml_graph.children("edge")) {
    double weight =
        edge.find_child_by_attribute("key", "weight").text().as_double();
    unsigned int source = edge.attribute("source").as_int() + 1;
    unsigned int target = edge.attribute("target").as_int() + 1;
    add_edge(source, target, weight);
  }

  // find cycles
  // each edge is on exactly one cycle -> go trough edges
  // compute edge list
  std::vector<std::pair<int, int>> cactus_edge_list;
  for (int i = 1; i <= n; ++i) {
    for (auto &e : adj_list[i]) {
      if (i < e.first) {
        cactus_edge_list.push_back({i, e.first});
      }
    }
  }
  // dfs to find cycle
  CycleID cycle_id;
  std::vector<bool> edge_visited;
  edge_visited.resize(cactus_edge_list.size(), false);
  for (int edge_idx = 0; edge_idx < cactus_edge_list.size(); ++edge_idx) {
    if (edge_visited[edge_idx]) {
      continue;
    }
    auto [source, target] = cactus_edge_list[edge_idx];
    std::list<int> queue;
    queue.push_back(source);
    std::vector<bool> visited;
    visited.resize(n + 1, false);
    visited[source] = true;
    std::vector<int> parent;
    parent.resize(n + 1, -1);
    parent[source] = source;
    while (!queue.empty()) {
      int next = queue.front();
      queue.pop_front();
      //   cycle_nodes.push_back(next);
      for (auto [t, _] : adj_list[next]) {
        // ignore initial edge
        if ((next == source && t == target) ||
            (next == target && t == source)) {
          continue;
        }

        if (!visited[t]) {
          queue.push_back(t);
          visited[t] = true;
          parent[t] = next;
        }
        if (t == target) {
          cycle_id = node_list.size();
          node_list.push_back({target});
          int current = target;
          while (current != source) {
            // mark edges of cycle as visited
            // FIXME: improve, binary search?
            for (int i = 0; i < cactus_edge_list.size(); ++i) {
              if ((current == cactus_edge_list[i].first &&
                   parent[current] == cactus_edge_list[i].second) ||
                  (current == cactus_edge_list[i].second &&
                   parent[current] == cactus_edge_list[i].first)) {
                edge_visited[i] = true;
              }
            }
            current = parent[current];
            node_list[cycle_id].push_back(current);
          }
          goto cycle_finished;
        }
      }
    }
    // edge not part of a cycle, edge is own cycle
    cycle_id = node_list.size();
    node_list.push_back({source, target});
  cycle_finished:
    for (CactusNodeID nd : node_list[cycle_id]) {
      node_to_cycle[nd].insert(cycle_id);
    }
    DEBUG("Add cycle of size " << node_list[cycle_id].size() << " for edge "
                               << source << "-" << target);
  }

  tree_n = node_list.size();
  edge_list.resize(tree_n);

  // add cycle tree edges
  for (int nd = 1; nd < node_to_cycle.size(); ++nd) {
    if (node_to_cycle[nd].size() > 1) {
      for (auto it_i = node_to_cycle[nd].begin();
           it_i != node_to_cycle[nd].end(); ++it_i) {
        auto it_j = it_i;
        ++it_j;
        for (; it_j != node_to_cycle[nd].end(); ++it_j) {
          // for (int i = 0; i < node_to_cycle[nd].size(); ++i) {
          // for (int j = i + 1; j < node_to_cycle[nd].size(); ++j) {
          edge_list[*it_i].push_back({nd, *it_j});
          edge_list[*it_j].push_back({nd, *it_i});
        }
      }
    }
  }

  DEBUG("Cycles: " << tree_n);

  finalize();
}

std::list<CactusNodeID> DynamicCactus::bfs(CactusNodeID u, CactusNodeID v) {
  std::unordered_set<CycleID> &cycles_u = node_to_cycle[u];
  std::unordered_set<CycleID> &cycles_v = node_to_cycle[v];

  // check if u and v are in the same cycle
  if (cycles_u.size() < cycles_v.size()) {
    if (std::any_of(cycles_u.begin(), cycles_u.end(), [&cycles_v](CycleID c) {
          return std::find(cycles_v.begin(), cycles_v.end(), c) !=
                 cycles_v.end();
        })) {
      return {v, u};
    }
  } else {
    if (std::any_of(cycles_v.begin(), cycles_v.end(), [&cycles_u](CycleID c) {
          return std::find(cycles_u.begin(), cycles_u.end(), c) !=
                 cycles_u.end();
        })) {
      return {v, u};
    }
  }

  // check if neither u nor v are articulation points
  if (cycles_u.size() == 1 && cycles_v.size() == 1) {
    CycleID cycle_u = *cycles_u.begin();
    CycleID cycle_v = *cycles_v.begin();
    // if they have the same articulation point, there is an edge => path found
    if (edge_list[cycle_u].front().first == edge_list[cycle_v].front().first) {
      return {v, edge_list[cycle_u].front().first, u};
    }
  }

  // bfs from multiple sources to multiple targets
  std::vector<bool> visited;
  visited.resize(tree_n, false);
  std::list<int> queue;
  std::vector<std::pair<CycleID, CactusNodeID>> parent;
  parent.resize(tree_n, {-1, -1});
  for (CycleID c : cycles_u) {
    visited[c] = true;
    queue.push_back(c);
    parent[c] = {c, u};
  }
  while (!queue.empty()) {
    CycleID next = queue.front();
    queue.pop_front();
    for (auto &[cactus_node, target] : edge_list[next]) {
      if (!visited[target]) {
        parent[target] = {next, cactus_node};
        if (std::find(cycles_v.begin(), cycles_v.end(), target) !=
            cycles_v.end()) {
          // reconstruct path
          std::pair<CycleID, CactusNodeID> current = {target, cactus_node};
          std::list<CactusNodeID> path = {v};
          while (current.first != parent[current.first].first) {
            path.push_back(parent[current.first].second);
            current = parent[current.first];
          }
          path.push_back(u);
          return path;
        }
        visited[target] = true;
        queue.push_back(target);
      }
    }
  }
  throw std::invalid_argument("no path");
}

void DynamicCactus::contract(CactusNodeID u, CactusNodeID v) {
  DEBUG("Contracting " << u << "-" << v);
  assert(exists_node(u) && exists_node(v));
  std::list<int> path = bfs(u, v);
  int last = path.front();
  // contract path
  std::for_each(std::next(path.begin()), path.end(), [&](CactusNodeID v) {
    DEBUG("Contracting in cycle: " << v << "-" << last);
    contract_in_cycle(v, last);
    last = v;
  });
  // Update mapping of link node ids to cactus nodes
  std::unordered_set<int> contracted_nodes;
  for (auto &x : path) {
    contracted_nodes.insert(x);
  }
  for (auto &target : original_to_cactus) {
    if (contracted_nodes.find(target) != contracted_nodes.end()) {
      target = last;
    }
  }
}

bool DynamicCactus::adjacent(CycleID cycle, CactusNodeID u, CactusNodeID v) {
  if ((node_list[cycle].front() == u && node_list[cycle].back() == v) ||
      (node_list[cycle].front() == v && node_list[cycle].back() == u)) {
    return true;
  }
  auto u_it = std::find(node_list[cycle].begin(), node_list[cycle].end(), u);
  auto u_next = u_it;
  u_next++;
  auto u_prev = u_it;
  u_prev--;
  return (u_it != node_list[cycle].begin() && *u_prev == v) ||
         (u_next != node_list[cycle].end() && *u_next == v);
}

// contract two nodes in the same cycle, node v will exist afterwards
void DynamicCactus::contract_in_cycle(CactusNodeID u, CactusNodeID v) {
  assert(u != v);
  merge_links(v, u);

  int cycle_id = common_cycle(u, v);
  // check if u and v are adjacent
  if (adjacent(cycle_id, u, v)) {
    // no new cycle, contraction may remove a cycle
    remove_cycle_edge(cycle_id, u, v);
  } else {
    // creates new cycle
    contract_articulation_points(cycle_id, u, v);
  }
  // update node_to_cycle mapping
  for (auto x : node_to_cycle[v]) {
    if (std::find(node_to_cycle[u].begin(), node_to_cycle[u].end(), x) ==
        node_to_cycle[u].end()) {
      node_to_cycle[u].insert(x);
    }
  }
  node_to_cycle[v].clear();
  --n;
}

void DynamicCactus::contract_articulation_points(CycleID cycle_id, int u,
                                                 int v) {
  // assumption: u, v not adjacent
  // make u the new articulation point
  // therefore, collect incident cycles
  std::list<CycleID> incident_cycles = {cycle_id};
  for (auto &[ap, cycle] : edge_list[cycle_id]) {
    if (ap == u || ap == v) {
      incident_cycles.push_back(cycle);
    }
  }
  // remove old incident edges
  for (auto &l : edge_list) {
    auto new_end =
        std::remove_if(l.begin(), l.end(), [u, v, cycle_id](auto &edge) {
          return (edge.first == u || edge.first == v);
        });
    l.resize(new_end - l.begin());
  }
  // Add new cycle
  node_list.push_back({});
  edge_list.push_back({});
  incident_cycles.push_back(tree_n);
  node_to_cycle[u].insert(tree_n);
  tree_n++;
  // Move nodes to new cycle
  auto it = node_list[cycle_id].begin();
  bool moving = true;
  while (it != node_list[cycle_id].end()) {
    CactusNodeID nd = *it;
    auto old_it = it;
    ++it;
    if (nd == u || nd == v) {
      moving = !moving;
      if (nd == u) {
        node_list[tree_n - 1].push_back(u);
      } else {
        // remove v
        node_list[cycle_id].erase(old_it);
      }
    } else {
      if (moving) {
        node_list[cycle_id].erase(old_it);
        node_list[tree_n - 1].push_back(nd);
        // fix pointers to new cycle
        // for (auto &l : edge_list) {
        for (CycleID i = 0; i < edge_list.size(); ++i) {
          for (auto &[ap, cycle] : edge_list[i]) {
            if (cycle == cycle_id && ap == nd) {
              cycle = tree_n - 1;
              // move edge
              edge_list[tree_n - 1].push_back({nd, i});
              auto it = std::find(edge_list[cycle_id].begin(),
                                  edge_list[cycle_id].end(),
                                  std::pair<int, int>{nd, i});
              if (it != edge_list[cycle_id].end()) {
                edge_list[cycle_id].erase(it);
              }
            }
          }
        }

        // Update node to cycle mapping
        node_to_cycle[nd].erase(cycle_id);
        node_to_cycle[nd].insert(tree_n - 1);
      }
    }
  }
  // change articulation point v to u
  for (auto &l : node_list) {
    for (auto &nd : l) {
      if (nd == v) {
        nd = u;
      }
    }
  }
  for (auto &l : edge_list) {
    for (auto &e : l) {
      if (e.first == v) {
        e.first = u;
      }
    }
  }
  // add new incident edges regarding articulation point u
  for (int first : incident_cycles) {
    for (int second : incident_cycles) {
      if (first < second) {
        edge_list[first].push_back({u, second});
        edge_list[second].push_back({u, first});
      }
    }
  }
}

void DynamicCactus::remove_cycle_edge(CycleID cycle_id, int u, int v) {
  if (node_list[cycle_id].size() == 2) {
    // delete cycle u -- v
    // remove cycle from node_to_cycle mapping
    node_to_cycle[u].erase(cycle_id);
    node_to_cycle[v].erase(cycle_id);
    // make u the new articulation point
    // therefore, collect incident edges
    std::set<CycleID> incident_cycles;
    for (auto &cycle_edge : edge_list[cycle_id]) {
      incident_cycles.insert(cycle_edge.second);
    }
    incident_cycles.erase(cycle_id);
    std::vector<CycleID> ic_vec;
    std::copy(incident_cycles.begin(), incident_cycles.end(),
              std::back_inserter(ic_vec));
    incident_cycles.clear();
    // remove old incident edges
    for (auto &l : edge_list) {
      auto new_end = std::remove_if(l.begin(), l.end(), [u, v](auto &edge) {
        return edge.first == u || edge.first == v;
      });
      l.resize(new_end - l.begin());
    }
    // add new incident edges
    for (int first = 0; first < ic_vec.size(); ++first) {
      for (int second = first + 1; second != ic_vec.size(); ++second) {
        edge_list[ic_vec[first]].push_back({u, ic_vec[second]});
        edge_list[ic_vec[second]].push_back({u, ic_vec[first]});
      }
    }

    // delete cycle
    std::swap(node_list[cycle_id], node_list.back());
    node_list.back().clear();
    node_list.pop_back();
    std::swap(edge_list[cycle_id], edge_list.back());
    edge_list.back().clear();
    edge_list.pop_back();
    --tree_n;
    // fix pointers to last cycle
    for (auto &l : node_list) {
      for (auto &val : l) {
        if (val == v)
          val = u;
      }
    }
    for (auto &l : edge_list) {
      for (auto &[_, val] : l) {
        if (val == tree_n)
          val = cycle_id;
      }
    }
    // fix node_to_cycle mapping
    for (auto &l : node_to_cycle) {
      if (l.find(tree_n) != l.end()) {
        l.erase(tree_n);
        l.insert(cycle_id);
      }
    }
  } else {
    // make cycle smaller
    // remove v from cycle
    node_list[cycle_id].remove(v);
    // fix pointers to v
    for (auto &el : edge_list) {
      for (auto &e : el) {
        if (e.first == v) {
          e.first = u;
        }
      }
    }
    // Add new links between cycles connected to u and v
    for (auto &x : node_to_cycle[u]) {
      if (x == cycle_id)
        continue;
      for (auto &y : node_to_cycle[v]) {
        if (y == cycle_id)
          continue;
        edge_list[x].push_back({u, y});
        edge_list[y].push_back({u, x});
      }
    }
    for (auto &nl : node_list) {
      for (CactusNodeID &articulation_point : nl) {
        if (articulation_point == v) {
          articulation_point = u;
        }
      }
    }
  }
}

int DynamicCactus::num_increased_cuts(CactusNodeID u, CactusNodeID v) {
  std::list<CactusNodeID> path = bfs(u, v);
  int sum = 0;
  std::reduce(std::next(path.begin()), path.end(), path.front(),
              [&sum, this](CactusNodeID a, CactusNodeID b) {
                sum += num_increased_cuts_cycle(a, b);
                return b;
              });
  return sum;
}

int DynamicCactus::num_increased_cuts_cycle(CactusNodeID u, CactusNodeID v) {
  CycleID cycle = common_cycle(u, v);
  int pos_u = 0, pos_v = 0;
  int cnt = 0;
  for (auto &nd : node_list[cycle]) {
    if (nd == u) {
      pos_u = cnt;
    } else if (nd == v) {
      pos_v = cnt;
    }
    ++cnt;
  }

  if (pos_u > pos_v) {
    std::swap(pos_u, pos_v);
  }

  return (pos_v - pos_u) * (node_list[cycle].size() - pos_v + pos_u);
}

CycleID DynamicCactus::common_cycle(CactusNodeID u, CactusNodeID v) {
  if (node_to_cycle[u].size() > node_to_cycle[v].size()) {
    std::swap(u, v);
  }
  if (node_to_cycle[u].size() == 1) {
    assert(std::find(node_to_cycle[v].begin(), node_to_cycle[v].end(),
                     *node_to_cycle[u].begin()) != node_to_cycle[v].end());
    return *node_to_cycle[u].begin();
  }
  // compute intersection
  std::vector<int> cycles_intersection;
  for (const int &element : node_to_cycle[u]) {
    if (node_to_cycle[v].find(element) != node_to_cycle[v].end()) {
      cycles_intersection.push_back(element);
    }
  }
  assert(cycles_intersection.size() == 1);
  return cycles_intersection.front();
}

void DynamicCactus::copy_links(const GraphPair &g) {
  // Initially: each node maps to itself
  // (links are already between cactus nodes)
  original_to_cactus.resize(num_nodes() + 1);
  std::iota(original_to_cactus.begin(), original_to_cactus.end(), 0);
  // copy all links into link array
  for (int u = 1; u <= g.cactus.num_nodes(); ++u) {
    for (int v = u + 1; v <= g.cactus.num_nodes(); ++v) {
      if (g.cactus.links[u][v].first) {
        links.push_back({u, v, g.cactus.links[u][v].weight});
      }
    }
  }
  if (links_per_node) {
    // sort all links into adjacency array
    adj_links.resize(num_nodes() + 1);
    for (Edge &l : links) {
      adj_links[l.first].push_back(l);
      adj_links[l.second].push_back(l);
    }
    for (auto &list : adj_links) {
      if (list.empty())
        continue;
      double best = std::numeric_limits<double>::max();
      for (auto &e : list) {
        if (e.weight < best)
          best = e.weight;
      }
      int old = list.size();
      list.remove_if([best](Edge &e) { return e.weight > 2 * best; });
    }
  }
}

void DynamicCactus::check_integrity() {
  for (CactusNodeID i = 1; i < node_to_cycle.size(); ++i) {
    if (node_to_cycle[i].size() > 1) {
      for (CycleID x : node_to_cycle[i]) {
        for (CycleID y : node_to_cycle[i]) {
          if (x != y) {
            assert(std::find(edge_list[x].begin(), edge_list[x].end(),
                             std::pair<int, int>{i, y}) != edge_list[x].end());
          }
        }
      }
    }
  }
  for (int i = 0; i < edge_list.size(); ++i) {
    for (auto &e : edge_list[i]) {
      assert(std::find(node_list[i].begin(), node_list[i].end(), e.first) !=
             node_list[i].end());
      assert(std::find(node_list[e.second].begin(), node_list[e.second].end(),
                       e.first) != node_list[e.second].end());
    }
  }
  DEBUG("Integrity success");
}

void DynamicCactus::reduce_links() {
  int old_size = links.size();
  std::vector<std::vector<double>> best;
  best.resize(num_nodes() + 1, {});
  for (int i = 0; i < best.size(); ++i) {
    best[i].resize(num_nodes() + 1, -1.);
  }
  for (auto &e : links) {
    int u = pi(e.first);
    int v = pi(e.second);
    if (u == v)
      continue;
    if (u > v)
      std::swap(u, v);
    if (best[u][v] < 0. || best[u][v] > e.weight) {
      best[u][v] = e.weight;
    }
  }
  links.remove_if([this, &best](Edge &l) {
    int u = pi(l.first);
    int v = pi(l.second);
    if (u > v)
      std::swap(u, v);
    return u == v || best[u][v] < l.weight;
  });
}

// DiGraph
void DiGraph::read_links_from_cactus(Cactus &cactus) {
  init(cactus.num_nodes(), cactus.num_edges() * 2);
  // for each undirected link, add two directed edges
  for (int u = 1; u < cactus.links.size(); ++u) {
    for (int v = u + 1; v < cactus.links[u].size(); ++v) {
      if (cactus.links[u][v].first) {
        add_edge(u, v, cactus.links[u][v].weight);
      }
    }
  }
  finalize();
}

void DynamicCactus::merge_links(int v, int u) {
  if (!links_per_node)
    return;
  adj_links[u].insert(adj_links[u].end(), adj_links[v].begin(),
                      adj_links[v].end());
  adj_links[u].remove_if(
      [this](Edge &e) { return pi(e.first) == pi(e.second); });
  adj_links[v].clear();
}

std::list<Edge> &DynamicCactus::links_of_node(int node) {
  return adj_links[node];
}

} // namespace graph
