#ifndef GRAPH_HPP
#define GRAPH_HPP

#include "util.hpp"

#include <filesystem>
#include <list>
#include <optional>
#include <set>
#include <unordered_set>
#include <vector>

#define PRINT_NODE_LIST                                                        \
  DEBUG("node_list: [");                                                       \
  for (auto &l : node_list) {                                                  \
    DEBUG("  [");                                                              \
    for (auto &val : l) {                                                      \
      DEBUG("    " << val << ", ");                                            \
    }                                                                          \
    DEBUG("  ],");                                                             \
  }                                                                            \
  DEBUG("]");
#define PRINT_EDGE_LIST                                                        \
  DEBUG("edge_list: [");                                                       \
  for (auto &l : edge_list) {                                                  \
    DEBUG("  [");                                                              \
    for (auto &val : l) {                                                      \
      DEBUG("    (" << val.first << " " << val.second << "), ");               \
    }                                                                          \
    DEBUG("  ],");                                                             \
  }                                                                            \
  DEBUG("]");

namespace graph {

// defined classes
class OriginalGraph;
class Cactus;
class GraphPair;
class CactusCut;

struct Edge {
  int first;
  int second;
  double weight;

  Edge(std::pair<int, int> edge)
      : first(edge.first), second(edge.second), weight(-1.) {}
  Edge(int u, int v, double _weight) : first(u), second(v), weight(_weight) {}
  bool operator==(const Edge &other) const {
    return first == other.first && second == other.second;
  }
};

class Graph {
public:
  unsigned int num_nodes() const;
  unsigned int num_edges() const;
  bool is_edge(int u, int v);
  std::vector<std::vector<int>> all_pair_distances();
  const std::list<std::pair<int, double>> &get_adj_list(int node);

protected:
  int n, m;
  std::vector<std::list<std::pair<int, double>>> adj_list;

  virtual void init(int n, int m);
  virtual void finalize();
  void add_edge(int u, int v, double weight = 1.);
  void add_directed_edge(int u, int v, double weight = 1.);
};

class OriginalGraph : public Graph {
public:
  void read_from_file(const std::filesystem::path &graph);
  void write_to_file(const std::filesystem::path &path);
  void add_link(Edge &e);
};

class Cactus : public Graph {
  friend class CactusCut;
  friend class GraphPair;

public:
  std::vector<std::vector<Edge>> links;
  std::vector<std::vector<double>> heuristic;

  void read_from_file(const std::filesystem::path &path, GraphPair &graph);
  double compute_min_cut();
  void expand_tree_edges();
  void eulerian_path(int u, std::list<int> &circuit); // destroys graph

private:
  std::vector<std::set<int>> contained_nodes;
  double min_cut = -1.;

  virtual void init(int n, int m);
};

class GraphPair {
  friend class Cactus;

public:
  OriginalGraph original_graph;
  Cactus cactus;

  void read_graph(const std::filesystem::path &graph,
                  const std::filesystem::path &cactus);

  int cactus_id(int original_graph_id);
  std::list<int> original_graph_ids(int cactus_id);

  void load_links_from_file(const std::filesystem::path &link_file);
  void generate_links(int seed = 0, double fraction = 1., int distribution = 1);
  void add_links(const std::filesystem::path &link_file, double fraction = 1.,
                 int distribution = 1);
  void add_links(int seed, double fraction, int distribution);
  void add_link(int u, int v, int original_u, int original_v, double weight);
  std::list<CactusCut> get_min_cuts();

private:
  std::vector<int> map_to_cactus;
  std::vector<std::list<int>> map_to_original_graph;
};

class CactusCut {
private:
  Edge e1;
  std::optional<Edge> e2;

public:
  std::vector<bool> partition;
  CactusCut(GraphPair &g, Edge e) : e1(e), e2(std::nullopt) {
    compute_partition(g);
  }
  CactusCut(GraphPair &g, Edge _e1, Edge _e2) : e1(_e1), e2({_e2}) {
    compute_partition(g);
  }

  bool is_edge_cut(int u, int v) { return partition[u] != partition[v]; }
  const bool operator[](int u) const { return partition[u]; }

  const Edge get_first() { return e1; }
  const std::optional<Edge> get_second() { return e2; }

private:
  /**
   * @brief Compute for all nodes of the graph on which side of the cut they are
   *
   * @param g the graph
   */
  void compute_partition(GraphPair &g);
};

inline bool is_same_edge(Edge &e, int u, int v) {
  return (e.first == u && e.second == v) || (e.first == v && e.second == u);
}

// Dynamic cactus graph data structure
#define CactusNodeID int
#define CycleID int

class DynamicCactus : public Graph {
public:
  // original graph related
  std::vector<int> original_to_cactus;
  // cactus graph related
  std::vector<std::vector<CycleID>>
      node_to_cycle_old; // Indexed by CactusNodeID
  std::vector<std::unordered_set<CycleID>> node_to_cycle;
  // cycle graph related
  int tree_n;
  std::vector<std::list<CactusNodeID>> node_list; // Indexed by CycleID
  std::vector<std::vector<std::pair<CactusNodeID, CycleID>>>
      edge_list; // Indexed by CycleID
  std::vector<std::list<Edge>> adj_links;

public:
  std::list<Edge> links;
  bool links_per_node = false;
  void read_from_file(const std::filesystem::path &graph);
  void copy_links(const GraphPair &g);
  void contract(CactusNodeID u, CactusNodeID v);
  void check_integrity();
  void reduce_links();
  int num_increased_cuts(CactusNodeID u, CactusNodeID v);
  std::list<CactusNodeID> bfs(CactusNodeID u, CactusNodeID v);
  int size() { return tree_n; }
  bool exists_node(CactusNodeID nd) { return !node_to_cycle[nd].empty(); }
  void debug() {
    PRINT_NODE_LIST;
    PRINT_EDGE_LIST;
    DEBUG("node_to_cycle: [");
    for (auto &l : node_to_cycle) {
      DEBUG("  [");
      for (auto &x : l) {
        DEBUG(x << ", ");
      }
      DEBUG("  ],");
    }
    DEBUG("]");
  }
  int pi(int node) { return original_to_cactus[node]; }
  std::list<Edge> &links_of_node(int node);
  int cycles_of_cactus_node(int node) { return node_to_cycle[pi(node)].size(); }

private:
  void contract_in_cycle(CactusNodeID u, CactusNodeID v);
  void remove_cycle_edge(CycleID cycle, int u, int v);
  void contract_articulation_points(CycleID cycle, int u, int v);
  int num_increased_cuts_cycle(CactusNodeID u, CactusNodeID v);
  CycleID common_cycle(CactusNodeID u, CactusNodeID v);
  bool adjacent(CycleID cycle, CactusNodeID u, CactusNodeID v);
  void merge_links(int u, int v);
};

// Directed graph
class DiGraph : public Graph {
public:
  void read_links_from_cactus(Cactus &cactus);
};

class AdjMatrixGraph {
public:
  std::vector<std::vector<double>> adj_matrix;
  int n;

  AdjMatrixGraph(int _n) : n(_n) {
    adj_matrix.resize(n + 1, {});
    for (auto &v : adj_matrix) {
      v.resize(n + 1, std::numeric_limits<double>::infinity());
    }
  }
};

} // namespace graph
#endif
