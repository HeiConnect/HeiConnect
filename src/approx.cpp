#include "approx.hpp"
#include "graph.hpp"
#include "gurobi_util.hpp"
#include "ring_graph.hpp"
#include "util.hpp"

#include <algorithm>
#include <cmath>
#include <exception>
#include <iterator>
#include <limits>
#include <numeric>
#include <stdexcept>

#include <fstream>
#include <gurobi_c++.h>
#include <string>
#include <vector>

namespace solver {
using namespace ring_graph;

// internal functions
std::vector<Cut> get_rooted_cuts_of_ring(RingGraph &g);
Solution directed_WRAP(RingGraph &g, std::vector<Cut> &cuts);
template <bool convert, typename LinkList>
std::list<graph::Edge> directed_to_undirected(LinkList &dir_solution);
void shorten(RingGraph &g, Solution &solution, std::vector<Cut> &cuts);
template <typename LinkList>
void write_solution_to_file(std::string path, LinkList &link_list, int n);
std::vector<Link> optimal_thin_component(RingGraph &g, std::vector<Cut> cuts,
                                         Solution &solution,
                                         std::list<Shadow> &F, int alpha);
template <bool apx15 = false>
std::pair<double, std::vector<Link *>>
max_slack_dynamic_program(RingGraph &g, std::vector<Cut> cuts,
                          Solution &solution, std::list<Shadow> &F, double rho,
                          int alpha, Apx15Solution *apx15Solution);
bool is_in_drop(LinkIntersectionGraph &h_k, Solution &solution, DirLink &l);
template <typename L> double c(std::vector<L> S);

std::list<graph::Edge> approximate_2_lp(graph::GraphPair &graph) {
  // reduce to ring graph
  RingGraph g(graph.cactus);
  // create directed instance
  std::vector<Cut> cuts = get_rooted_cuts_of_ring(g);
  DEBUG("Nubmer of rooted cuts: " << cuts.size());
  // solve directed version using LP
  Solution dWRAP_solution = directed_WRAP(g, cuts);

  return directed_to_undirected<true>(dWRAP_solution.directed_solution);
}

std::list<graph::Edge> approximate_1_5_e(graph::GraphPair &graph,
                                         double epsilon) {
  // reduce to ring graph
  RingGraph g(graph.cactus);
  // create directed instance
  std::vector<Cut> cuts = get_rooted_cuts_of_ring(g);
  DEBUG("Number of rooted cuts: " << cuts.size());
  // Compute constants
  int alpha = 4 * std::ceil(4. / epsilon);
  double factor = 1 - (epsilon / (6 * g.n));
  // solve directed version using LP
  Solution solution = directed_WRAP(g, cuts);
  // shorten links of directed solution
  shorten(g, solution, cuts);
  // Compute F and witness sets from shortened links
  Apx15Solution sol(solution);
  double Phi = sol.Phi();
  double old_Phi;
  do {
    std::list<Shadow> F;
    std::copy(solution.shortening.begin(), solution.shortening.end(),
              std::back_inserter(F));
    auto [_pi, K] =
        max_slack_dynamic_program<true>(g, cuts, solution, F, 1., alpha, &sol);
    // compute drop
    LinkIntersectionGraph h_k(K); // construct H[K]
    F.remove_if([&h_k, &solution](Shadow &l) {
      // check if v is connected to a v-good vertex in H[K]
      return is_in_drop(h_k, solution, l);
    });
    solution.directed_solution.clear();
    // Add links not dropped
    for (auto &s : F) {
      solution.directed_solution.push_back(*s.original);
    }
    // Add K
    for (auto *l : K) {
      solution.directed_solution.push_back(
          {l->u, l->v, l->weight, l->original_u, l->original_v});
      solution.directed_solution.push_back(
          {l->v, l->u, l->weight, l->original_u, l->original_v});
    }
    shorten(g, solution, cuts);
    sol = Apx15Solution(solution);

    old_Phi = Phi;
    Phi = sol.Phi();
  } while (Phi / old_Phi <= factor);

  // convert solution to links in the original graph, removing duplicates
  std::list<graph::Edge> solution_original_graph;
  std::set<std::pair<int, int>> already_added;
  for (auto &l : sol.F) {
    if (l.weight > 0. && already_added.find({l.original_u, l.original_v}) ==
                             already_added.end()) {
      already_added.insert({l.original_u, l.original_v});
      solution_original_graph.push_back({l.original_u, l.original_v, l.weight});
    }
  }
  return solution_original_graph;
}

std::list<graph::Edge> approximate_1_ln2_e(graph::GraphPair &graph,
                                           double epsilon) {
  // reduce to ring graph
  RingGraph g(graph.cactus);
  // create directed instance
  std::vector<Cut> cuts = get_rooted_cuts_of_ring(g);
  DEBUG("Number of rooted cuts: " << cuts.size());
  // solve directed version using LP
  Solution solution = directed_WRAP(g, cuts);
  // shorten links of directed solution
  shorten(g, solution, cuts);

  std::list<Shadow> F; // non-shortenable solution
  std::copy(solution.shortening.begin(), solution.shortening.end(),
            back_inserter(F));
  std::vector<Link> S; // undirected replacement
  // add all links with weight 0
  for (auto &l : g.links) {
    if (l.weight == 0.) {
      S.push_back(l);
    }
  }

  for (auto &x : F) {
    std::cout << x.u << "-" << x.v << " ";
  }
  std::cout << std::endl;
  // eventually drop links
  LinkIntersectionGraph h_k(S); // construct H[K]
  // remove drop (all covered directed links)
  F.remove_if([&h_k, &solution](Shadow &l) {
    // check if v is connected to a v-good vertex in H[K]
    return is_in_drop(h_k, solution, l);
  });

  int alpha = 4 * std::ceil(2 / epsilon);
  while (!F.empty()) {
    for (auto &x : F) {
      std::cout << x.u << "-" << x.v << " ";
    }
    std::cout << std::endl;
    // compute a alpha-thin minimizer K
    std::vector<Link> K = optimal_thin_component(g, cuts, solution, F, alpha);
    // apply K
    for (auto &l : K) {
      DEBUG("Add link: " << l.u << "-" << l.v);
    }
    std::copy(K.begin(), K.end(), back_inserter(S));
    double cost = c(K);
    double gain = 0.;
    LinkIntersectionGraph h_k(K); // construct H[K]
    // remove drop (all covered directed links)
    F.remove_if([&h_k, &solution, &gain](Shadow &l) {
      // check if v is connected to a v-good vertex in H[K]
      bool drop = is_in_drop(h_k, solution, l);
      if (drop) {
        gain += l.weight;
      }
      return drop;
    });
    DEBUG("Cost: " << cost << ", gain: " << gain
                   << ", ratio: " << (cost / gain))
  }

  // convert solution to links in the original graph, removing duplicates
  std::list<graph::Edge> solution_original_graph;
  std::set<std::pair<int, int>> already_added;
  for (auto &l : S) {
    if (l.weight > 0. && already_added.find({l.original_u, l.original_v}) ==
                             already_added.end()) {
      already_added.insert({l.original_u, l.original_v});
      solution_original_graph.push_back({l.original_u, l.original_v, l.weight});
    }
  }
  return solution_original_graph;
}

std::vector<Cut> get_rooted_cuts_of_ring(RingGraph &g) {
  // each range between 1 and n-1 is a cut
  std::vector<Cut> cuts;

  for (int u = 1; u < g.n; ++u) {
    for (int v = u; v < g.n; ++v) {
      cuts.push_back({u, v});
    }
  }

  return cuts;
}

Solution directed_WRAP(RingGraph &g, std::vector<Cut> &cuts) {
  try {
    INFO("Creating directed WRAP LP");
    // Create Gurobi LP model
    GRBEnv env = GRBEnv(true);
    env.set("LogFile", "dWRAP.log");
    env.start();
    GRBModel model = GRBModel(env);

    // get list of directed links
    std::vector<DirLink> dir_links;
    dir_links.reserve(g.links.size() * 2);
    for (auto &l : g.links) {
      dir_links.push_back({l.u, l.v, l.weight, l.original_u, l.original_v});
      dir_links.push_back({l.v, l.u, l.weight, l.original_u, l.original_v});
    }
    // variables, one for each directed link
    std::vector<GRBVar> vars;
    vars.reserve(dir_links.size());
    for (int i = 0; i < dir_links.size(); ++i) {
      vars.push_back(
          model.addVar(0.0, 1.0, 1.0, GRB_CONTINUOUS, "v" + std::to_string(i)));
    }
    model.update();
    INFO("Number of variables = " << model.get(GRB_IntAttr_NumVars));

    // constraints, one for each cut
    std::vector<GRBConstr> constraints;
    constraints.reserve(cuts.size());
    int count = 0;
    for (auto &cut : cuts) {
      GRBLinExpr expr = 0;
      for (int i = 0; i < dir_links.size(); ++i) {
        if (cut.is_cut_directed(dir_links[i])) {
          expr += vars[i];
        }
      }
      GRBConstr con = model.addConstr(expr >= 1, "c" + std::to_string(count++));
      constraints.push_back(con);
    }
    model.update();
    INFO("Number of constraints = " << model.get(GRB_IntAttr_NumConstrs));

    // objective
    GRBLinExpr obj = 0;
    for (int i = 0; i < dir_links.size(); ++i) {
      obj += vars[i] * dir_links[i].weight;
    }
    model.setObjective(obj, GRB_MINIMIZE);
    model.update();

    // solve
    solve_lp(model);

    // get solution
    std::vector<DirLink> augmentation;
    for (int i = 0; i < vars.size(); ++i) {
      int val = vars[i].get(GRB_DoubleAttr_X);
      if (val > 0.) {
        augmentation.push_back(dir_links[i]);
      }
    }

    return {augmentation, {}};
  } catch (GRBException e) {
    ERROR("GRBException, error code = " << e.getErrorCode());
    WARN(e.getMessage());
    throw e;
  }
}

// get undirected solution from directed one, not adding the same link twice
template <bool convert, typename LinkList>
std::list<graph::Edge> directed_to_undirected(LinkList &directed_solution) {
  std::list<graph::Edge> solution;
  std::set<std::pair<int, int>> already_added;
  for (auto &directed_link : directed_solution) {
    if (directed_link.weight > 0.) {
      std::pair<int, int> e =
          directed_link.u < directed_link.v
              ? std::pair<int, int>{directed_link.u, directed_link.v}
              : std::pair<int, int>{directed_link.v, directed_link.u};
      if (!(already_added.find(e) != already_added.end())) {
        already_added.insert(e);
        if constexpr (convert) {
          solution.push_back({directed_link.original_u,
                              directed_link.original_v, directed_link.weight});
        } else {
          solution.push_back(
              {directed_link.u, directed_link.v, directed_link.weight});
        }
      }
    }
  }
  return solution;
}

void shorten(RingGraph &g, Solution &solution, std::vector<Cut> &cuts) {
  solution.shortening.clear();
  solution.shortening_adj = {};
  // count cut coverages
  std::vector<int> num_coverages;
  num_coverages.resize(cuts.size());
  for (auto &link : solution.directed_solution) {
    for (int i = 0; i < cuts.size(); ++i) {
      if (cuts[i].is_cut_directed(link)) {
        ++num_coverages[i];
      }
    }
  }

  // go through all links, only adding necessary and shortest possible links to
  // shortening
  for (auto &link : solution.directed_solution) {
    // check if link can be deleted
    bool link_required = false;
    for (int i = 0; i < cuts.size(); ++i) {
      if (num_coverages[i] <= 1 && cuts[i].is_cut_directed(link)) {
        link_required = true;
        break;
      }
    }
    if (!link_required) {
      DEBUG("Remove " << link.u << "-" << link.v);
      // do not add to shortened solution, but remove coverages
      for (int i = 0; i < cuts.size(); ++i) {
        if (cuts[i].is_cut_directed(link)) {
          --num_coverages[i];
        }
      }
      continue;
    }
    // add shortest shortening
    int target = link.v;
    int dir = target > link.u ? -1 : 1;
    for (int source = target + dir; source - dir != link.u; source += dir) {
      // check if shortening is valid
      bool valid = true;
      for (int i = 0; i < cuts.size(); ++i) {
        if (num_coverages[i] <= 1 && cuts[i].is_cut_directed(link) &&
            !cuts[i].is_cut_directed(source, target)) {
          valid = false;
          break;
        }
      }
      if (valid) {
        // add shortened link, updating coverages
        solution.shortening.push_back({source, target, link.weight,
                                       link.original_u, link.original_v,
                                       &link});
        for (int i = 0; i < cuts.size(); ++i) {
          if (cuts[i].is_cut_directed(link) &&
              !cuts[i].is_cut_directed(source, target)) {
            --num_coverages[i];
          }
        }
        break;
      }
    }
  }
  // set compute adjacency lists for shortening
  solution.compute_shortening_adj(g.n);
}

template <typename LinkList>
void write_solution_to_file(std::string path, LinkList &link_list, int n) {
  std::vector<std::vector<int>> adj_lists;
  std::set<std::pair<int, int>> already_added;
  adj_lists.resize(n, {});
  for (int i = 0; i < n; ++i) {
    adj_lists[i].push_back((i + 1) % n);
    adj_lists[(i + 1) % n].push_back(i);
  }
  for (auto &l : link_list) {
    int u = l.u, v = l.v;
    if (v < u)
      std::swap(u, v);
    std::pair<int, int> uv = {u, v};
    if (already_added.find(uv) == already_added.end()) {
      adj_lists[l.u].push_back(l.v);
      adj_lists[l.v].push_back(l.u);
      already_added.insert(uv);
    }
  }

  std::ofstream file;
  file.open(path);
  file << n << " " << link_list.size() << std::endl;
  for (int i = 0; i < n; ++i) {
    for (int nd : adj_lists[i]) {
      file << (nd + 1) << " ";
    }
    file << std::endl;
  }
  file.close();
}

// Lemma 2.10: a shadow (u, v) is in the drop if it is connected to a v-good
// vertex in H[S]
bool is_in_drop(LinkIntersectionGraph &h_k, Solution &solution, DirLink &l) {
  RingInterval v_bad = solution.get_v_bad(l.v);
  return h_k.is_connected_to_v_good(l.v, v_bad);
}

std::vector<Link> optimal_thin_component(RingGraph &g, std::vector<Cut> cuts,
                                         Solution &solution,
                                         std::list<Shadow> &F, int alpha) {
  // find rho* along with K using binary search w.r.t. rho
  // compute the maximum interval size in which we need to find rho*
  // scale up the minimum cost to 1
  // FIXME: may not be sufficient in every case because the proof requires
  // integer weights. However, it was in examples tested
  double cost_F = 0.;
  double min_cost = std::numeric_limits<double>::infinity();
  for (auto &l : solution.directed_solution) {
    if (l.weight > 0. && l.weight < min_cost) {
      min_cost = l.weight;
    }
  }
  for (auto &l : F) {
    cost_F += l.weight / min_cost;
  }
  double required_interval_size =
      std::min<double>(0.5, 1. / (cost_F * cost_F)); // 1/c(F')^2

  std::vector<Link *> K;
  // binary search
  double lower_bound = 0., upper_bound = 1.;
  while (upper_bound - lower_bound > required_interval_size) {
    double middle = (lower_bound + upper_bound) / 2.;
    auto [eta, _K] =
        max_slack_dynamic_program(g, cuts, solution, F, middle, alpha, nullptr);
    K = std::move(_K);
    if (eta > 0.) {
      upper_bound = middle;
    } else {
      lower_bound = middle;
    }
  }

  std::vector<Link> K_deref;
  for (Link *l : K)
    K_deref.push_back(*l);
  return K_deref;
}

template <typename L> double c(std::vector<L> S) {
  double sum = 0.;
  for (L &l : S) {
    sum += l.weight;
  }
  return sum;
}

template <typename L> double c(std::vector<L *> S) {
  double sum = 0.;
  for (L *l : S) {
    sum += l->weight;
  }
  return sum;
}

template <bool apx15>
double pi(Solution &solution, Pattern &pattern, std::vector<Link *> &S,
          std::list<Shadow> &F, double rho, Apx15Solution *apx15solution) {
  std::vector<Shadow> set;
  LinkIntersectionGraph h(S);
  // FIXME: ensure shadow \in F \cap F_0?
  for (Shadow &shadow : F) {
    // check if link is incoming to C, i.e. it is in intersection of (10)
    if (pattern.C.in_range(shadow.v)) {
      if (is_in_drop(h, solution, shadow)) {
        set.push_back(shadow);
      }
    }
  }
  double cost = -c(S);
  if constexpr (apx15) {
    cost *= 1.5;
    for (auto &l : set) {
      cost += l.weight / apx15solution->witness_set_size(l);
    }
  } else {
    cost += rho * c(set);
  }
  return cost;
}

template <bool apx15>
double combined_pi(DPEntry &e1, DPEntry &e2, Solution &solution,
                   std::list<Shadow> &F, double rho, Pattern &tao_bar,
                   Apx15Solution *apx15solution) {
  // previous pi
  double pi = e1.pi + e2.pi;

  // c(B_1 \cap B_2)
  std::vector<Link *> B_intersection;
  std::set_intersection(e1.pattern.B.begin(), e1.pattern.B.end(),
                        e2.pattern.B.begin(), e2.pattern.B.end(),
                        std::back_inserter(B_intersection));
  if constexpr (apx15) {
    pi += 1.5 * c(B_intersection);
  } else {
    pi += c(B_intersection);
  }

  // \tilde{c}(\Cup_(u \in U_{Q_1, Q_2})\delta^-_{F_0}(u))
  std::set<int> U;
  for (int T = 0; T < e1.pattern.Tao_size; ++T) {
    if (e1.pattern.psi[T] == 1 && e1.pattern.C.in_range(e1.pattern.phi[T])) {
      U.insert(e1.pattern.phi[T]);
    }
  }
  for (int T = 0; T < e2.pattern.Tao_size; ++T) {
    if (e2.pattern.psi[T] == 1 && e2.pattern.C.in_range(e2.pattern.phi[T])) {
      U.insert(e2.pattern.phi[T]);
    }
  }
  // remove combined lca from U
  {
    for (auto u = U.begin(), last = U.end(); u != last;) {
      if (std::find(tao_bar.phi.begin(), tao_bar.phi.end(), *u) !=
          tao_bar.phi.end())
        u = U.erase(u);
      else
        ++u;
    }
    // C++20 version
    // std::remove_if(U.begin(), U.end(), [&](const int u) {
    //  return std::find(tao_bar.phi.begin(), tao_bar.phi.end(), u) !=
    //         tao_bar.phi.end();
    //});
  }

  // FIXME: ensure shadow \in F \cap F_0?
  double cost = 0.;
  for (Shadow &shadow : F) {
    if (U.find(shadow.v) != U.end()) {
      if constexpr (apx15) {
        cost += shadow.weight / apx15solution->witness_set_size(shadow);
      } else {
        cost += rho * shadow.weight;
      }
    }
  }
  pi += cost;

  return pi;
}

template <bool apx15>
void add_entries_for_1cut(DPTable &table, int v, RingGraph &g, int alpha,
                          Solution &solution, std::list<Shadow> &F, double rho,
                          Apx15Solution *apx15solution) {
  // compute all feasible links that have one endpoint in C = {v}
  std::vector<Link *> superset_B;
  for (auto &link : g.links) {
    if ((link.u == v) != (link.v == v)) {
      if constexpr (apx15) {
        superset_B.push_back(&link);
      } else {
        if (link.weight != 0.) {
          superset_B.push_back(&link);
        }
      }
    }
  }
  if (superset_B.size() == 0) {
    return;
  }
  // enumerate all index subsets of size at most alpha
  std::list<std::vector<int>> enumerations =
      enumerate(superset_B.size(), alpha);
  // Add the entries to the dp table
  for (auto &enumeration : enumerations) {
    std::vector<Link *> B = get_indices(superset_B, enumeration);
    std::vector<Link *> S;
    for (Link *x : B) {
      S.push_back(x);
    }
    Pattern pattern = Pattern::create(g, solution, B, v);
    DPEntry entry = {pattern, S,
                     pi<apx15>(solution, pattern, S, F, rho, apx15solution)};
    table[1][v][pattern.hash_B()] = entry;
  }
}

template <bool apx15>
std::pair<double, std::vector<Link *>>
max_slack_dynamic_program(RingGraph &g, std::vector<Cut> cuts,
                          Solution &solution, std::list<Shadow> &F, double rho,
                          int alpha, Apx15Solution *apx15Solution) {
  DPTable dp_table;
  dp_table.resize(g.n + 1, {});
  for (auto &row : dp_table) {
    row.resize(g.n + 1, {});
  }
  // Initial row (|C| = 1)
  for (int i = 1; i < g.n; ++i) {
    add_entries_for_1cut<apx15>(dp_table, i, g, alpha, solution, F, rho,
                                apx15Solution);
  }
  for (int i = 2; i < g.n; ++i) {
    // go through all combinations resulting in |C| = i
    int half = i / 2;
    for (int first_cut_size = 1; first_cut_size <= half; ++first_cut_size) {
      int second_cut_size = i - first_cut_size;
      for (int first_cut_size_idx = 0;
           first_cut_size_idx < dp_table[first_cut_size].size();
           ++first_cut_size_idx) {
        for (int second_cut_size_idx =
                 first_cut_size == second_cut_size ? first_cut_size_idx + 1 : 0;
             second_cut_size_idx < dp_table[second_cut_size].size();
             ++second_cut_size_idx) {
          for (auto &[_, e1] : dp_table[first_cut_size][first_cut_size_idx]) {
            for (auto &[_, e2] :
                 dp_table[second_cut_size][second_cut_size_idx]) {
              if (Pattern::compatible(e1.pattern, e2.pattern, alpha)) {
                /* Tao bar (Definition 8.9) */
                // (a)
                std::vector<Link *> B = set_union(e1.pattern.B, e2.pattern.B);
                LinkIntersectionGraph h(B);
                std::unordered_map<Link *, int> link_to_h_nodeid;
                for (int i = 0; i < B.size(); ++i) {
                  link_to_h_nodeid[B[i]] = i;
                }
                // (b)
                for (int i = 0; i < e1.pattern.B.size(); ++i) {
                  if (link_to_h_nodeid.find(e1.pattern.B[i]) !=
                      link_to_h_nodeid.end()) {
                    for (int j = i + 1; j < e1.pattern.B.size(); ++j) {
                      if (link_to_h_nodeid.find(e1.pattern.B[i]) !=
                          link_to_h_nodeid.end()) {
                        // going through all pairs in B1 that are also in B
                        if (e1.pattern.Tao[i] == e1.pattern.Tao[j]) {
                          // FIXME: may add multiple edges, only affects
                          // performance
                          h.add_link(link_to_h_nodeid[e1.pattern.B[i]],
                                     link_to_h_nodeid[e1.pattern.B[j]]);
                        }
                      }
                    }
                  }
                }
                // (c)
                for (int i = 0; i < e2.pattern.B.size(); ++i) {
                  if (link_to_h_nodeid.find(e2.pattern.B[i]) !=
                      link_to_h_nodeid.end()) {
                    for (int j = i + 1; j < e2.pattern.B.size(); ++j) {
                      if (link_to_h_nodeid.find(e2.pattern.B[i]) !=
                          link_to_h_nodeid.end()) {
                        // going through all pairs in B1 that are also in B
                        if (e2.pattern.Tao[i] == e2.pattern.Tao[j]) {
                          // FIXME: may add multiple edges, only affects
                          // performance
                          h.add_link(link_to_h_nodeid[e2.pattern.B[i]],
                                     link_to_h_nodeid[e2.pattern.B[j]]);
                        }
                      }
                    }
                  }
                }
                auto [Tao, Tao_size] = h.connected_components();
                Pattern tmp_pattern = {
                    Cut::combine_unchecked(e1.pattern.C, e2.pattern.C), B,
                    Tao_size, Tao};
                tmp_pattern.compute_phi(solution);

                Pattern new_pattern = Pattern::merge_compatible(
                    e1.pattern, e2.pattern, solution, tmp_pattern);
                size_t hash = new_pattern.hash_B();
                // FIXME: recomputing is faster - why?
                double new_pi = combined_pi<apx15>(e1, e2, solution, F, rho,
                                                   tmp_pattern, apx15Solution);
                // double new_pi = pi(solution, new_pattern, new_S, F, rho);
                // FIXME: overwrites best solution becaus pi is 0, but is
                // getting ignored later. This results in non-optimal solutions
                // being applied
                if (!e1.realizer.empty() || !e2.realizer.empty()) {
                  auto &map = dp_table[i][new_pattern.C.first];
                  if (map.find(hash) == map.end() || map[hash].pi < new_pi) {
                    std::vector<Link *> new_S;
                    std::set_union(e1.realizer.begin(), e1.realizer.end(),
                                   e2.realizer.begin(), e2.realizer.end(),
                                   std::back_inserter(new_S));
                    map[hash] = {new_pattern, new_S, new_pi};
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  // find max
  DPEntry *max_entry;
  double max_pi = -std::numeric_limits<double>::max();
  for (auto &cut_size_row : dp_table) {
    for (auto &cut_row : cut_size_row) {
      for (auto &[hash, entry] : cut_row) {
        if (!entry.realizer.empty() && entry.pi > max_pi) {
          max_pi = entry.pi;
          max_entry = &entry;
        }
      }
    }
  }

  return {max_pi, max_entry->realizer};
}

} // namespace solver
