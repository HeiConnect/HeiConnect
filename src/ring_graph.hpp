#ifndef APX_GRAPH_HPP
#define APX_GRAPH_HPP

#include "graph.hpp"
#include "util.hpp"
#include <algorithm>
#include <boost/functional/hash.hpp>
#include <cassert>
#include <iterator>
#include <limits>
#include <pstl/glue_algorithm_defs.h>
#include <stack>
#include <unordered_map>
#include <vector>

namespace ring_graph {

struct Link {
  int u, v;
  double weight;
  int original_u, original_v;

  bool crossing(Link &link) {
    // check if they share an endpoint
    if (u == link.u || u == link.v || v == link.u || v == link.v) {
      return true;
    }
    // check if they are strictly crossing
    return in_range(link.u) != in_range(link.v);
  }

private:
  bool in_range(int nd) { return nd > u && nd < v; }
};

struct DirLink : Link {};

struct Shadow : DirLink {
  DirLink *original;
};

struct RingInterval {
  int first, last;

  bool is_cut(Link l) { return is_cut(l.u, l.v); }
  bool is_cut_directed(Link l) { return is_cut_directed(l.u, l.v); }
  bool is_cut(int u, int v) { return in_range(u) != in_range(v); }
  bool is_cut_directed(int u, int v) { return !in_range(u) && in_range(v); }

  bool in_range(int u) { return u >= first && u <= last; }
};

struct Solution {
  std::vector<DirLink> directed_solution;
  std::vector<Shadow> shortening;
  std::vector<std::list<int>> shortening_adj;

  RingInterval get_v_bad(int v);
  void compute_shortening_adj(int n);

  int lca(int u, int v);
  std::list<int> dfs(int source, int target);
};

// Represents a cut; stores the first and last node in the cut
struct Cut : RingInterval {
  static Cut combine_unchecked(Cut c1, Cut c2) {
    return {std::min(c1.first, c2.first), std::max(c1.last, c2.last)};
  }
};

// implicit ring graph, node ids range from 0 to n-1 and are connected in order
class RingGraph {

  // mapping to previous nodes
  std::vector<int> mapping;
  // links
  std::vector<std::list<Link *>> adj_lists;

public:
  std::vector<Link> links;
  int n, num_links = 0;

  RingGraph(graph::Cactus cactus);

private:
  void add_link(int u, int v, double weight = 0., int original_u = 0,
                int original_v = 0);
};

class LinkIntersectionGraph {
  int n;
  std::vector<std::list<int>> adj_lists;
  std::vector<Link *> link_of_node;
  // stores adjacent links for each node
  std::vector<std::set<int>> adjacent_links;

public:
  LinkIntersectionGraph(std::vector<Link> &K);
  LinkIntersectionGraph(std::vector<Link *> &K);
  bool is_connected_to_v_good(int v, RingInterval &v_bad);
  void add_link(int u, int v) {
    adj_lists[u].push_back(v);
    adj_lists[v].push_back(u);
  }
  std::pair<std::vector<int>, int> connected_components();
};

struct Pattern {
  Cut C;
  std::vector<Link *> B; // δ_S(C)
  int Tao_size;
  std::vector<int> Tao;  // partition of B
  std::vector<int> phi;  // φ[i] = lca(Tao[i])
  std::vector<bool> psi; // ψ[i] = φ[i] ∈ V(S_i)
  Pattern *parent1;
  Pattern *parent2;

  size_t hash_B() {
    std::sort(B.begin(), B.end());
    size_t seed = 0;
    boost::hash_range(seed, B.begin(), B.end());
    return seed;
  }

  static Pattern create(RingGraph &g, Solution &solution, std::vector<Link *> S,
                        int cut_node) {
    // create C-pattern with C = {cut_node}
    Pattern pattern = {{cut_node, cut_node}, S};

    // All links in B are in the same component,
    // because they are all connected to cut_node
    pattern.Tao_size = pattern.B.size() ? 1 : 0;
    pattern.Tao.resize(pattern.B.size(), 0);

    pattern.compute_phi(solution);
    pattern.compute_psi();
    return pattern;
  }

  static bool compatible(Pattern &p1, Pattern &p2, int alpha) {
    // cuts must be neighboring, so that combination is also a cut
    bool neighboring =
        (p1.C.first == p2.C.last + 1) || (p1.C.last + 1 == p2.C.first);
    if (!neighboring)
      return false;

    // links between C1 and C2 must be the same
    // FIXME: may be simpler
    std::set<Link *> d_B1_C2;
    for (Link *b1 : p1.B) {
      if (p2.C.in_range(b1->u) != p2.C.in_range(b1->v)) {
        d_B1_C2.insert(b1);
      }
    }
    for (Link *b2 : p2.B) {
      if (p1.C.in_range(b2->u) != p1.C.in_range(b2->v)) {
        // in d_{B_2}(C_1)
        auto found = d_B1_C2.find(b2);
        if (found == d_B1_C2.end()) {
          // but not in d_{B_1}(C_2)
          return false;
        }
        d_B1_C2.erase(found);
      }
    }
    if (!d_B1_C2.empty()) {
      // there are elements in d_{B_1}(C_2) that are not in d_{B_2}(C_1)
      return false;
    }

    // new pattern must be alpha thin
    std::vector<Link *> B12 = p1.B;
    std::copy(p2.B.begin(), p2.B.end(), std::back_inserter(B12));
    Cut C12 = Cut::combine_unchecked(p1.C, p2.C);
    // elements up to this iterator have exactly one endpoint in C12
    std::vector<Link *>::iterator d_B12_C12_it_end =
        std::remove_if(B12.begin(), B12.end(), [&C12](Link *link) {
          return C12.in_range(link->u) != C12.in_range(link->v);
        });
    return (d_B12_C12_it_end - B12.begin()) <= alpha;
  }

  static Pattern merge_compatible(Pattern &p1, Pattern &p2, Solution &solution,
                                  Pattern &tao_bar) {
    /* C */
    Cut C = Cut::combine_unchecked(p1.C, p2.C);

    /* B, Tao, Tao_size */
    std::vector<Link *> B;
    std::vector<int> Tao;
    B.reserve(tao_bar.B.size());
    Tao.reserve(tao_bar.B.size());
    int Tao_size = 0;
    std::vector<int> partition_relabel(tao_bar.Tao_size, -1);
    for (int i = 0; i < tao_bar.B.size(); ++i) {
      if (!(C.in_range(tao_bar.B[i]->u) && C.in_range(tao_bar.B[i]->v))) {
        if (partition_relabel[tao_bar.Tao[i]] < 0) {
          partition_relabel[tao_bar.Tao[i]] = Tao_size++;
        }
        B.push_back(tao_bar.B[i]);
        Tao.push_back(partition_relabel[tao_bar.Tao[i]]);
      }
    }

    /* phi, psi */
    Pattern new_pattern = {C, B, Tao_size, Tao};
    new_pattern.compute_phi(solution);
    new_pattern.compute_psi();

    new_pattern.parent1 = &p1;
    new_pattern.parent2 = &p2;

    return new_pattern;
  }

  void compute_phi(Solution &solution) {
    phi.resize(Tao_size);
    // lca(V(S)) = lca(leftmost(V(S)), rightmost(V(S)))
    for (int i = 0; i < Tao_size; ++i) {
      int min = std::numeric_limits<int>::max();
      int max = std::numeric_limits<int>::min();
      for (auto &idx : Tao) {
        Link &l = *B[idx];
        if (idx == i) {
          if (l.u < min) {
            min = l.u;
          }
          if (l.v < min) {
            min = l.v;
          }
          if (l.u > max) {
            max = l.u;
          }
          if (l.v > max) {
            max = l.v;
          }
        }
      }
      phi[i] = solution.lca(min, max);
    }
  }

  void compute_psi() {
    psi.resize(Tao_size, false);
    for (int idx = 0; idx < Tao.size(); ++idx) {
      Link &l = *B[idx];
      int partition = Tao[idx];
      if (l.u == phi[partition] || l.v == phi[partition]) {
        psi[partition] = true;
      }
    }
  }
};

struct DPEntry {
  Pattern pattern;
  std::vector<Link *> realizer;
  double pi;
};

// [cut size][cut start node][B]
using DPTable = std::vector<std::vector<std::unordered_map<long, DPEntry>>>;

struct Apx15Solution {
  std::vector<Link> F;
  std::vector<std::list<Shadow>> W;
  std::unordered_map<int, int> shadow_to_index;

  Apx15Solution(Solution &solution) {
    for (auto &link : solution.shortening) {
      int idx = hash(link);
      if (shadow_to_index.find(idx) != shadow_to_index.end()) {
        W[shadow_to_index[idx]].push_back(link);
      } else {
        shadow_to_index[idx] = F.size();
        F.push_back(*link.original);
        W.push_back({link});
      }
    }
  }

  int witness_set_size(Shadow &shadow) {
    assert(shadow_to_index.find(hash(shadow)) != shadow_to_index.end());
    return W[shadow_to_index[hash(shadow)]].size();
  }

  double Phi() {
    double sum = 0.;
    for (int i = 0; i < F.size(); ++i) {
      sum += F[i].weight * (W[i].size() == 2 ? 1.5 : 1.);
    }
    return sum;
  }

  void drop(Shadow &shadow) {}

private:
  int hash(Shadow &shadow) {
    if (shadow.original->u < shadow.original->v) {
      return shadow.original->u * 100000 + shadow.original->v;
    } else {
      return shadow.original->u + shadow.original->v * 100000;
    };
  }
};

} // namespace ring_graph

#endif
