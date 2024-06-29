#include "greedy.hpp"
#include "alternating_paths.hpp"
#include "graph.hpp"
#include "ilp.hpp"
#include "kruskal_mst.hpp"
#include "max_flow.hpp"
#include "util.hpp"

#include <algorithm>
// #include <bits/chrono.h>
#include <cassert>
#include <chrono>
#include <cmath>
#include <exception>
#include <functional>
#include <limits>
#include <numeric>
#include <queue>
#include <random>
#include <stdexcept>
#include <vector>

namespace solver {

std::list<graph::Edge> greedy_weak(graph::GraphPair &g) {
  INFO("Running weak greedy algorithm");
  std::list<graph::Edge> augmentation;
  std::list<graph::CactusCut> cuts = g.get_min_cuts();
  while (cuts.size() > 0) {
    graph::CactusCut cut = cuts.front();
    cuts.pop_front();
    // if (std::any_of(augmentation.begin(), augmentation.end(),
    //                 [&](graph::Edge &val) {
    //                   return cut.is_edge_cut(val.first, val.second);
    //                 })) {
    //   continue;
    // }
    // greedily add cheapest edge to increase cut
    graph::Edge best_edge = {0, 0, std::numeric_limits<double>::max()};
    for (int u = 1; u <= g.cactus.num_nodes(); ++u) {
      for (int v = u + 1; v <= g.cactus.num_nodes(); ++v) {
        if (g.cactus.links[u][v].first && cut.is_edge_cut(u, v) &&
            g.cactus.links[u][v].weight < best_edge.weight) {
          best_edge.weight = g.cactus.links[u][v].weight;
          best_edge.first = u;
          best_edge.second = v;
        }
      }
    }

    if (!best_edge.first) {
      INFO("Cut: " << cut.get_first().first << ", " << cut.get_first().second);
      for (int i = 0; i < g.cactus.num_nodes(); ++i) {
        DEBUG(i << ": " << cut.partition[i]);
      }
      FATAL("Links cannot improve cut");
    }
    augmentation.push_back(best_edge);
    // remove all cuts that were increased
    cuts.remove_if([&](graph::CactusCut cut) {
      return cut.is_edge_cut(best_edge.first, best_edge.second);
    });
    DEBUG("Cuts remaining: " << cuts.size());
  }

  return augmentation;
}

std::list<graph::Edge> greedy_strong(graph::GraphPair &g) {
  std::list<graph::Edge> augmentation;
  std::list<graph::CactusCut> cuts = g.get_min_cuts();
  while (cuts.size() > 0) {
    // greedily add cheapest edge of all cuts to increase any cut
    graph::Edge best_edge = {0, 0, std::numeric_limits<double>::max()};
    for (auto &cut : cuts) {
      for (int u = 1; u <= g.cactus.num_nodes(); ++u) {
        for (int v = u + 1; v <= g.cactus.num_nodes(); ++v) {
          if (g.cactus.links[u][v].first && cut.is_edge_cut(u, v) &&
              g.cactus.links[u][v].weight < best_edge.weight) {
            best_edge.weight = g.cactus.links[u][v].weight;
            best_edge.first = u;
            best_edge.second = v;
          }
        }
      }
    }

    if (!best_edge.first) {
      FATAL("Links cannot improve cut");
    }
    augmentation.push_back(best_edge);
    // remove all cuts that were increased
    cuts.remove_if([&](graph::CactusCut cut) {
      return cut.is_edge_cut(best_edge.first, best_edge.second);
    });
    DEBUG("Cuts remaining: " << cuts.size());
  }
  return augmentation;
}

std::list<graph::Edge> greedy_heuristic_strong(graph::GraphPair &g) {
  std::list<graph::Edge> augmentation;
  std::list<graph::CactusCut> cuts = g.get_min_cuts();
  while (cuts.size() > 0) {
    // greedily add cheapest edge to increase cut
    graph::Edge best_edge = {0, 0, 0.};
    double best_heuristic = std::numeric_limits<double>::max();
    for (int u = 1; u <= g.cactus.num_nodes(); ++u) {
      for (int v = u + 1; v <= g.cactus.num_nodes(); ++v) {
        if (g.cactus.links[u][v].first) {
          int num_increased_cuts =
              std::accumulate(cuts.begin(), cuts.end(), 0,
                              [u, v](int result, graph::CactusCut &b) {
                                return result + b.is_edge_cut(u, v);
                              });
          if (num_increased_cuts) {
            double heuristic = g.cactus.links[u][v].weight / num_increased_cuts;
            if (heuristic < best_heuristic) {
              best_heuristic = heuristic;
              best_edge.weight = g.cactus.links[u][v].weight;
              best_edge.first = u;
              best_edge.second = v;
            }
          }
        }
      }
    }

    if (!best_edge.first) {
      FATAL("Links cannot improve cut");
    }
    augmentation.push_back(best_edge);
    // remove all cuts that were increased
    cuts.remove_if([&](graph::CactusCut cut) {
      return cut.is_edge_cut(best_edge.first, best_edge.second);
    });
    DEBUG("Cuts remaining: " << cuts.size());
  }

  return augmentation;
}

std::list<graph::Edge> greedy_heuristic_sampling(graph::GraphPair &g,
                                                 int buckets) {
  std::list<graph::Edge> augmentation;
  std::list<graph::CactusCut> cuts = g.get_min_cuts();
  int num_cuts = cuts.size();
  if (!buckets) {
    buckets = std::sqrt(cuts.size());
  }
  INFO("Heuristic sampling using " << buckets << " buckets");
  std::vector<std::list<graph::CactusCut>> all_cuts;
  all_cuts.resize(buckets, {});
  int mod = 0;
  for (auto &cut : cuts) {
    all_cuts[mod].push_back(cut);
    mod = ++mod % buckets;
  }
  mod = 0;
  while (num_cuts > 5 * buckets) {
    if (all_cuts[mod].size() == 0) {
      mod = ++mod % buckets;
      continue;
    }
    // greedily add cheapest edge to increase cut
    graph::Edge best_edge = {0, 0, 0.};
    double best_heuristic = std::numeric_limits<double>::max();
    for (int u = 1; u <= g.cactus.num_nodes(); ++u) {
      for (int v = u + 1; v <= g.cactus.num_nodes(); ++v) {
        if (g.cactus.links[u][v].first) {
          int num_increased_cuts =
              std::accumulate(all_cuts[mod].begin(), all_cuts[mod].end(), 0,
                              [u, v](int result, graph::CactusCut &b) {
                                return result + b.is_edge_cut(u, v);
                              });
          if (num_increased_cuts) {
            double heuristic = g.cactus.links[u][v].weight / num_increased_cuts;
            if (heuristic < best_heuristic) {
              best_heuristic = heuristic;
              best_edge.weight = g.cactus.links[u][v].weight;
              best_edge.first = u;
              best_edge.second = v;
            }
          }
        }
      }
    }

    if (!best_edge.first) {
      FATAL("Links cannot improve cut");
    }
    augmentation.push_back(best_edge);
    num_cuts = 0;
    for (auto &c : all_cuts) { // to avoid duplicate edges
      c.remove_if([&](graph::CactusCut cut) {
        return cut.is_edge_cut(best_edge.first, best_edge.second);
      });
      num_cuts += c.size();
    }
    mod = ++mod % buckets;
    DEBUG("Cuts remaining: " << num_cuts);
  }
  cuts.clear();
  for (auto &c : all_cuts) {
    cuts.insert(cuts.end(), c.begin(), c.end());
  }

  while (cuts.size() > 0) {
    // greedily add cheapest edge to increase cut
    graph::Edge best_edge = {0, 0, 0.};
    double best_heuristic = std::numeric_limits<double>::max();
    for (int u = 1; u <= g.cactus.num_nodes(); ++u) {
      for (int v = u + 1; v <= g.cactus.num_nodes(); ++v) {
        if (g.cactus.links[u][v].first) {
          int num_increased_cuts =
              std::accumulate(cuts.begin(), cuts.end(), 0,
                              [u, v](int result, graph::CactusCut &b) {
                                return result + b.is_edge_cut(u, v);
                              });
          if (num_increased_cuts) {
            double heuristic = g.cactus.links[u][v].weight / num_increased_cuts;
            if (heuristic < best_heuristic) {
              best_heuristic = heuristic;
              best_edge.weight = g.cactus.links[u][v].weight;
              best_edge.first = u;
              best_edge.second = v;
            }
          }
        }
      }
    }

    if (!best_edge.first) {
      FATAL("Links cannot improve cut");
    }
    augmentation.push_back(best_edge);
    // remove all cuts that were increased
    cuts.remove_if([&](graph::CactusCut cut) {
      return cut.is_edge_cut(best_edge.first, best_edge.second);
    });
    DEBUG("Cuts remaining: " << cuts.size());
  }

  return augmentation;
}

std::list<graph::Edge> greedy_dynamic(graph::DynamicCactus &g) {
  DEBUG("Cactus size: " << g.size());
  std::list<graph::Edge> augmentation;
  while (g.size() > 0) {
    // greedily add cheapest edge to increase cut
    graph::Edge best_edge = {0, 0, 0.};
    graph::Edge original_edge = {0, 0, 0.};
    double best_heuristic = std::numeric_limits<double>::max();

    auto it = g.links.begin();

    while (it != g.links.end()) {
      auto [orig_u, orig_v, weight] = *it;
      int u = g.pi(orig_u);
      int v = g.pi(orig_v);
      if (u == v) {
        auto old_it = it;
        ++it;
        g.links.erase(old_it);
        continue;
      }
      int num_increased_cuts = g.num_increased_cuts(u, v);
      double heuristic = weight / num_increased_cuts;
      if (heuristic < best_heuristic) {
        best_heuristic = heuristic;
        best_edge.weight = weight;
        best_edge.first = u;
        best_edge.second = v;
        original_edge.first = orig_u;
        original_edge.second = orig_v;
        original_edge.weight = weight;
      }
      ++it;
    }

    if (!best_edge.first) {
      FATAL("Links cannot improve cut");
      return augmentation;
    }
    augmentation.push_back(original_edge);
    g.contract(best_edge.first, best_edge.second);
    DEBUG("Cactus size: " << g.size() << ", links: " << g.links.size());

    // if (cnt % 5 == 0)
    g.reduce_links();
  }

  return augmentation;
}

std::list<graph::Edge> greedy_dynamic_sampling(graph::DynamicCactus &g) {
  DEBUG("Cactus size: " << g.size());
  std::list<graph::Edge> augmentation;
  int n = g.links.size();
  int num_buckets = std::sqrt(n);
  std::mt19937 mt;
  std::uniform_int_distribution<> gen{0, num_buckets - 1};
  std::vector<std::list<graph::Edge>> buckets{static_cast<size_t>(num_buckets),
                                              std::list<graph::Edge>{}};
  for (auto &l : g.links) {
    buckets[gen(mt)].push_back(l);
  }
  int current_bucket = 0;
  while (g.size() > 0 && g.num_nodes() > num_buckets) {
    // greedily add cheapest edge to increase cut
    graph::Edge best_edge = {0, 0, 0.};
    graph::Edge original_edge = {0, 0, 0.};
    double best_heuristic = std::numeric_limits<double>::max();

    auto it = buckets[current_bucket].begin();

    while (it != buckets[current_bucket].end()) {
      auto [orig_u, orig_v, weight] = *it;
      int u = g.pi(orig_u);
      int v = g.pi(orig_v);
      if (u == v) {
        auto old_it = it;
        ++it;
        buckets[current_bucket].erase(old_it);
        continue;
      }
      int num_increased_cuts = g.num_increased_cuts(u, v);
      double heuristic = weight / num_increased_cuts;
      if (heuristic < best_heuristic) {
        best_heuristic = heuristic;
        best_edge.weight = weight;
        best_edge.first = u;
        best_edge.second = v;
        original_edge.first = orig_u;
        original_edge.second = orig_v;
        original_edge.weight = weight;
      }
      ++it;
    }

    if (!best_edge.first) {
      current_bucket = (current_bucket + 1) % num_buckets;
      continue;
    }
    augmentation.push_back(original_edge);
    g.contract(best_edge.first, best_edge.second);
    DEBUG("Cactus size: " << g.size() << ", links: " << g.links.size());
  }
  std::list<graph::Edge> links;
  for (auto &l : buckets) {
    links.insert(links.end(), l.begin(), l.end());
  }

  return augmentation;
}

std::list<graph::Edge> greedy_dynamic_bounds(graph::DynamicCactus &g) {
  DEBUG("Cactus size: " << g.size());
  std::list<graph::Edge> augmentation;

  // lower bound buckets
  int num_lists = 32;
  std::vector<double> list_bounds;
  std::vector<std::list<graph::Edge>> lists;
  list_bounds.resize(num_lists);
  lists.resize(num_lists, {});
  for (int i = 0; i < num_lists; ++i) {
    list_bounds[i] = 1. / std::pow(2, i);
  }
  int cnt = 0;
  while (!g.links.empty()) {
    auto &x = g.links.back();
    int u = g.pi(x.first);
    int v = g.pi(x.second);
    int num_increased_cuts = g.num_increased_cuts(u, v);
    double heuristic = x.weight / num_increased_cuts;
    int to_list =
        std::min<double>(std::log2(1. / heuristic), num_lists - 2) + 1;
    lists[to_list].push_back(x);
    g.links.pop_back();
  }
  for (int i = 0; i < num_lists; ++i) {
    DEBUG("Bucket " << i << ": " << lists[i].size());
  }

  while (g.size() > 0) {
    // greedily add cheapest edge to increase cut
    graph::Edge best_edge = {0, 0, 0.};
    graph::Edge original_edge = {0, 0, 0.};
    double best_heuristic = std::numeric_limits<double>::max();

    for (int i = num_lists - 1; i >= 0; --i) {
      auto &bucket = lists[i];
      if (bucket.empty())
        continue;
      if (best_edge.first && list_bounds[i] > best_heuristic) {
        break;
      }
      auto it = bucket.begin();

      while (it != bucket.end()) {
        auto [orig_u, orig_v, weight] = *it;
        int u = g.pi(orig_u);
        int v = g.pi(orig_v);
        if (u == v) {
          auto old_it = it;
          ++it;
          bucket.erase(old_it);
          continue;
        }
        int num_increased_cuts = g.num_increased_cuts(u, v);
        double heuristic = weight / num_increased_cuts;
        if (heuristic < best_heuristic) {
          best_heuristic = heuristic;
          best_edge.weight = weight;
          best_edge.first = u;
          best_edge.second = v;
          original_edge.first = orig_u;
          original_edge.second = orig_v;
          original_edge.weight = weight;
        }
        if (i > 0 && heuristic > list_bounds[i - 1]) {
          lists[std::log2(1. / heuristic) + 1].push_back(*it);
          auto old_it = it;
          ++it;
          bucket.erase(old_it);
        } else {
          ++it;
        }
      }
    }

    if (!best_edge.first) {
      return augmentation;
      FATAL("Links cannot improve cut");
    }
    augmentation.push_back(original_edge);
    g.contract(best_edge.first, best_edge.second);
    int cnt = 0;
    for (auto &x : lists)
      cnt += x.size();
    DEBUG("Cactus size: " << g.size() << ", links: " << cnt);
  }

  return augmentation;
}

std::list<graph::Edge> greedy_mst(graph::GraphPair &g) {
  // Create graph from links in cactus graph
  kruskal::KruskalGraph kruskal_graph(g.cactus);
  std::list<graph::Edge> mst = kruskal_graph.kruskal_mst();

  // sort mst by descending weight (we want to remove expensive ones first)
  mst.sort([](auto &a, auto &b) { return a.weight > b.weight; });

  // compute cuts
  INFO("Computing cuts");
  std::vector<std::pair<int, graph::CactusCut>> cuts;
  for (auto &cut : g.get_min_cuts()) {
    cuts.push_back({0, cut});
  }

  // Count number of coverages per cut
  INFO("Counting cut coverages");
  for (graph::Edge &link : mst) {
    for (auto &cut : cuts) {
      if (cut.second.is_edge_cut(link.first, link.second))
        ++cut.first;
    }
    // remove if link is not necessary
  }

  // Remove unnecessary links
  INFO("Removing unnecessary links from MST");
  mst.remove_if([&cuts](graph::Edge &link) {
    bool can_remove =
        std::none_of(cuts.begin(), cuts.end(),
                     [&link](std::pair<int, graph::CactusCut> &cut) {
                       return cut.second.is_edge_cut(link.first, link.second) &&
                              cut.first <= 1;
                     });
    if (can_remove) {
      for (auto &cut : cuts) {
        if (cut.second.is_edge_cut(link.first, link.second)) {
          --cut.first;
        }
      }
    }
    return can_remove;
  });

  return mst;
}

std::list<graph::Edge> greedy_mst_ilp(graph::GraphPair &g, int trees) {
  // store links that should be used
  std::vector<std::vector<graph::Edge>> new_links;
  new_links.resize(g.cactus.links.size(), {});
  for (auto &vec : new_links) {
    vec.resize(g.cactus.links[0].size(), {0, 0, 0.});
  }

  std::list<graph::Edge> mst;

  // Add msts
  for (int i = 0; i < trees; ++i) {
    kruskal::KruskalGraph kruskal_graph(g.cactus);
    std::list<graph::Edge> new_mst = kruskal_graph.kruskal_mst();
    for (auto &e : new_mst) {
      new_links[e.first][e.second] = g.cactus.links[e.first][e.second];
      g.cactus.links[e.first][e.second] = {0, 0, 0.};
    }
    mst.splice(mst.end(), new_mst);
  }

  g.cactus.links = new_links;

  return ilp(g);
}

std::pair<std::list<graph::Edge>, std::list<graph::Edge>>
greedy_mst_max_flow(graph::GraphPair &g) {
  // Create graph from links in cactus graph
  std::list<graph::Edge> mst;
  {
    kruskal::KruskalGraph kruskal_graph(g.cactus);
    mst = kruskal_graph.kruskal_mst();
  }
  // sort mst by descending weight (we want to remove expensive ones first)
  mst.sort([](auto &a, auto &b) { return a.weight > b.weight; });

  std::list<graph::Edge> edges;
  double cut_threshold = g.cactus.compute_min_cut() - 0.5;
  for (int u = 1; u <= g.cactus.num_nodes(); u++) {
    for (auto &[v, weight] : g.cactus.get_adj_list(u)) {
      if (u < v) {
        double w = weight > cut_threshold ? 2. : 1.;
        edges.push_back({u, v, w});
      }
    }
  }

  // use max flow algorithm to find out if we can drop edges
  std::list<graph::Edge> solution;
  auto it = mst.begin();
  while (it != mst.end()) {
    auto link = *it;
    int flow = max_flow::max_flow(g.cactus.num_nodes(), mst, edges, link);
    auto old_it = it;
    ++it;
    if (flow < 3) {
      solution.push_back(link);
    } else {
      mst.erase(old_it);
    }
  }

  return {solution, mst};
}

std::list<graph::Edge> full_mst(graph::GraphPair &g) {
  // Create graph from links in cactus graph
  kruskal::KruskalGraph kruskal_graph(g.cactus);
  return kruskal_graph.kruskal_mst();
}

std::list<graph::Edge> greedy_2mst_localsearch(graph::GraphPair &g,
                                               int depth_limit) {
  // Compute greedy removal solution
  auto [mst_solution, mst_rest] = greedy_mst_max_flow(g);

  // compute cuts
  INFO("Computing cuts");
  std::vector<std::pair<int, graph::CactusCut>> cuts;
  for (auto &cut : g.get_min_cuts()) {
    cuts.push_back({0, cut});
  }

  // Count number of coverages per cut
  INFO("Counting cut coverages");
  for (graph::Edge &link : mst_solution) {
    for (auto &cut : cuts) {
      if (cut.second.is_edge_cut(link.first, link.second))
        ++cut.first;
    }
    // remove if link is not necessary
  }

  // Compute 2 MSTs (MST from max flow got destroyed)
  std::list<std::pair<std::pair<int, int>, graph::Edge>> readd;
  {
    kruskal::KruskalGraph kruskal_graph(g.cactus);
    auto mst = kruskal_graph.kruskal_mst();
    for (auto &e : mst) {
      readd.push_back({{e.first, e.second}, g.cactus.links[e.first][e.second]});
      g.cactus.links[e.first][e.second] = {0, 0, 0.};
    }
  }
  {
    kruskal::KruskalGraph kruskal_graph(g.cactus);
    std::list<graph::Edge> new_mst = kruskal_graph.kruskal_mst();
    mst_rest.splice(mst_rest.end(), new_mst);
  }
  for (auto &[pos, e] : readd) {
    g.cactus.links[pos.first][pos.second] = e;
  }

  alternating_paths::AlternatingGraph alternating_graph(g.cactus.num_nodes(),
                                                        depth_limit);
  alternating_graph.add_solution(mst_solution);
  alternating_graph.add_links(mst_rest);
  bool one_path_valid = true;
  while (one_path_valid) {
    auto paths = alternating_graph.find_alternating_path();
    DEBUG("Found " << paths.size() << " paths");
    one_path_valid = false;
    paths.sort([](auto &a, auto &b) { return a.second < b.second; });
    for (auto &[path, pot] : paths) {
      std::vector<int> diff;
      diff.resize(cuts.size(), 0);
      // check if augmenting path is still a connectivity augmentation
      for (auto &e : path) {
        int sign = e->active ? -1. : 1.;
        for (int i = 0; i < cuts.size(); ++i) {
          if (cuts[i].second.is_edge_cut(e->source, e->target)) {
            diff[i] += sign;
          }
        }
      }
      bool valid = true;
      for (int i = 0; i < cuts.size(); ++i) {
        if (diff[i] + cuts[i].first <= 0) {
          valid = false;
          break;
        }
      }
      // if it is valid, use it
      if (valid) {
        DEBUG("Applying best path with potential " << pot);
        std::cout << "Path: ";
        for (auto &e : path) {
          std::cout << e->source << "-" << e->target << " ";
        }
        std::cout << std::endl;

        for (int i = 0; i < cuts.size(); ++i) {
          cuts[i].first += diff[i];
        }

        for (auto &e : path) {
          if (e->active) {
            mst_solution.remove({e->source, e->target, e->weight});
            mst_solution.remove({e->target, e->source, e->weight});
          } else {
            mst_solution.push_back({e->source, e->target, e->weight});
          }
        }
        alternating_graph.apply_path(path);
        one_path_valid = true;
        break;
      }
    }
  }

  return mst_solution;
}

std::list<graph::Edge> greedy_2mst_localsearch_flow(graph::GraphPair &g,
                                                    int depth_limit, bool cache,
                                                    int trees) {
  if (!trees)
    trees = 2;
  DEBUG("Using " << trees << " MSTs");
  // Compute greedy removal solution
  auto [mst_solution, mst_rest] = greedy_mst_max_flow(g);

  auto begin = std::chrono::steady_clock::now();
  double solution_weight =
      std::accumulate(mst_solution.begin(), mst_solution.end(), 0.,
                      [](double v, graph::Edge &e) { return v + e.weight; });

  std::list<graph::Edge> edges;
  double cut_threshold = g.cactus.compute_min_cut() - 0.5;
  for (int u = 1; u <= g.cactus.num_nodes(); u++) {
    for (auto &[v, weight] : g.cactus.get_adj_list(u)) {
      if (u < v) {
        double w = weight > cut_threshold ? 2. : 1.;
        edges.push_back({u, v, w});
      }
    }
  }

  // Compute 2 MSTs (MST from max flow got destroyed)
  std::list<std::pair<std::pair<int, int>, graph::Edge>> readd;
  {
    kruskal::KruskalGraph kruskal_graph(g.cactus);
    auto mst = kruskal_graph.kruskal_mst();
    // temporarily remove for second MST
    for (auto &e : mst) {
      readd.push_back({{e.first, e.second}, g.cactus.links[e.first][e.second]});
      g.cactus.links[e.first][e.second] = {0, 0, 0.};
    }
  }
  for (int i = 1; i < trees; ++i) {
    kruskal::KruskalGraph kruskal_graph(g.cactus);
    std::list<graph::Edge> new_mst = kruskal_graph.kruskal_mst();
    mst_rest.splice(mst_rest.end(), new_mst);
  }
  for (auto &[pos, e] : readd) {
    g.cactus.links[pos.first][pos.second] = e;
  }

  alternating_paths::AlternatingGraph alternating_graph(g.cactus.num_nodes(),
                                                        depth_limit);
  alternating_graph.add_solution(mst_solution);
  alternating_graph.add_links(mst_rest);
  bool one_path_valid = true;
  while (one_path_valid) {
    auto paths = alternating_graph.find_alternating_path();
    DEBUG("Found " << paths.size() << " paths");
    one_path_valid = false;
    paths.sort([](auto &a, auto &b) { return a.second < b.second; });
    for (auto &[path, pot] : paths) {
      std::list<graph::Edge> removed_links;
      std::list<graph::Edge> added_links;
      for (auto &e : path) {
        if (e->active) {
          removed_links.push_back({e->source, e->target, e->weight});
        } else {
          added_links.push_back({e->source, e->target, e->weight});
        }
      }

      bool valid = true;
      max_flow::ReusableMaxFlow flow_graph;
      flow_graph.init(g.cactus.num_nodes(), mst_solution, edges, removed_links,
                      added_links);
      for (auto &e : removed_links) {
        int flow =
            //  flow_graph.max_flow(e.first, e.second);
            max_flow::max_flow(g.cactus.num_nodes(), e.first, e.second,
                               mst_solution, edges, removed_links, added_links);
        if (flow < 3) {
          valid = false;
          break;
        }
      }

      // if it is valid, use it
      if (valid) {
        DEBUG("Applying best path with potential " << pot);
        solution_weight += pot;
        auto now = std::chrono::steady_clock::now();
        std::cout << "Path: ";
        for (auto &e : path) {
          std::cout << e->source << "-" << e->target << " ";
        }
        std::cout << std::endl;

        for (auto &e : path) {
          if (e->active) {
            mst_solution.remove({e->source, e->target, e->weight});
            mst_solution.remove({e->target, e->source, e->weight});
          } else {
            mst_solution.push_back({e->source, e->target, e->weight});
          }
        }
        alternating_graph.apply_path(path);
        one_path_valid = true;
        break;
      } else if (cache) {
        alternating_graph.invalidate_path(path);
      }
    }
  }

  auto now = std::chrono::steady_clock::now();

  return mst_solution;
}

std::list<graph::Edge>
greedy_mst_max_flow_heuristic(graph::GraphPair &g,
                              graph::DynamicCactus &dyn_cactus) {
  // Create graph from links in cactus graph
  std::list<graph::Edge> mst;
  {
    kruskal::KruskalGraph kruskal_graph(g.cactus, dyn_cactus);
    mst = kruskal_graph.kruskal_mst();
  }
  // sort mst by descending weight (we want to remove expensive ones first)
  mst.sort([](auto &a, auto &b) { return a.weight > b.weight; });

  std::list<graph::Edge> edges;
  double cut_threshold = g.cactus.compute_min_cut() - 0.5;
  for (int u = 1; u <= g.cactus.num_nodes(); u++) {
    for (auto &[v, weight] : g.cactus.get_adj_list(u)) {
      if (u < v) {
        double w = weight > cut_threshold ? 2. : 1.;
        edges.push_back({u, v, w});
      }
    }
  }

  // use max flow algorithm to find out if we can drop edges
  std::list<graph::Edge> solution;
  auto it = mst.begin();
  while (it != mst.end()) {
    auto link = *it;
    int flow = max_flow::max_flow(g.cactus.num_nodes(), mst, edges, link);
    auto old_it = it;
    ++it;
    if (flow < 3) {
      solution.push_back(link);
    } else {
      mst.erase(old_it);
    }
  }

  return solution;
}

std::list<graph::Edge>
greedy_mst_max_flow_order_heuristic(graph::GraphPair &g,
                                    graph::DynamicCactus &dyn_cactus) {
  // Create graph from links in cactus graph
  std::list<graph::Edge> mst;
  {
    kruskal::KruskalGraph kruskal_graph(g.cactus);
    mst = kruskal_graph.kruskal_mst();
  }
  // sort mst by descending weight (we want to remove expensive ones first)
  mst.sort([&dyn_cactus](graph::Edge &a, graph::Edge &b) {
    int ai = dyn_cactus.num_increased_cuts(a.first, a.second);
    int bi = dyn_cactus.num_increased_cuts(b.first, b.second);
    double ha =
        ai == 0 ? std::numeric_limits<double>::infinity() : a.weight / ai;
    double hb =
        bi == 0 ? std::numeric_limits<double>::infinity() : b.weight / bi;
    return ha > hb;
  });

  std::list<graph::Edge> edges;
  double cut_threshold = g.cactus.compute_min_cut() - 0.5;
  for (int u = 1; u <= g.cactus.num_nodes(); u++) {
    for (auto &[v, weight] : g.cactus.get_adj_list(u)) {
      if (u < v) {
        double w = weight > cut_threshold ? 2. : 1.;
        edges.push_back({u, v, w});
      }
    }
  }

  // use max flow algorithm to find out if we can drop edges
  std::list<graph::Edge> solution;
  auto it = mst.begin();
  while (it != mst.end()) {
    auto link = *it;
    int flow = max_flow::max_flow(g.cactus.num_nodes(), mst, edges, link);
    auto old_it = it;
    ++it;
    if (flow < 3) {
      solution.push_back(link);
      dyn_cactus.contract(dyn_cactus.pi(link.first),
                          dyn_cactus.pi(link.second));
    } else {
      mst.erase(old_it);
    }
  }

  return solution;
}

} // namespace solver
