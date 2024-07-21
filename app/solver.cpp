#include "approx.hpp"
#include "config.hpp"
#include "graph.hpp"
#include "greedy.hpp"
#include "ilp.hpp"
#include "util.hpp"
#include "watanabe.hpp"

#include <numeric>
#include <optional>
#include <stdexcept>

using namespace graph;

int main(int argc, char **argv) {
  config::Params params{argc, argv};

  GraphPair g;
  TIMEIT("Read graph", g.read_graph(params.original_graph, params.cactus));
  if (params.seed >= 0) {
    TIMEIT("Generate links",
           g.add_links(params.seed, params.link_fraction, params.distribution));
  } else {
    TIMEIT("Generate links", g.add_links(params.links));
  }

  std::list<graph::Edge> solution;
  bool convert_links = true;

  switch (params.algorithm) {
  case config::Algorithm::GWC:
  case config::Algorithm::GWC_SAMPLING: {
    convert_links = true;
    DynamicCactus g_dynamic;
    g_dynamic.read_from_file(params.cactus);
    g_dynamic.copy_links(g);
    if (params.algorithm == config::Algorithm::GWC_SAMPLING) {
      solution = solver::greedy_dynamic_sampling(g_dynamic);
    } else {
      solution = solver::greedy_dynamic_bounds(g_dynamic);
    }
  } break;
  case config::Algorithm::EILP:
    convert_links = false;
    solution = solver::ilp(g, params.use_initial, params.count);
    break;
  case config::Algorithm::GREEDY:
    solution = solver::greedy_weak(g);
    break;
  case config::Algorithm::GREEDY_GLOBAL:
    solution = solver::greedy_strong(g);
    break;
  case config::Algorithm::HEURISTIC:
    solution = solver::greedy_heuristic_strong(g);
    break;
  case config::Algorithm::APX2_LP:
    convert_links = false;
    solution = solver::approximate_2_lp(g);
    break;
  case config::Algorithm::MST:
    solution = solver::greedy_mst(g);
    break;
  case config::Algorithm::MST_CONNECT:
    solution = solver::greedy_mst_max_flow(g).first;
    break;
  case config::Algorithm::MST_CONNECT_HEURISTIC: {
    DynamicCactus g_dynamic;
    g_dynamic.read_from_file(params.cactus);
    g_dynamic.copy_links(g);
    solution = solver::greedy_mst_max_flow_heuristic(g, g_dynamic);
    break;
  }
  case config::Algorithm::MST_CONNECT_ORDER_HEURISTIC: {
    DynamicCactus g_dynamic;
    g_dynamic.read_from_file(params.cactus);
    g_dynamic.copy_links(g);
    solution = solver::greedy_mst_max_flow_order_heuristic(g, g_dynamic);
    break;
  }
  case config::Algorithm::FULL_MST:
    solution = solver::full_mst(g);
    break;
  case config::Algorithm::MST_CONNECT_LS:
    solution = solver::greedy_2mst_localsearch(g, params.depth);
    break;
  case config::Algorithm::MST_CONNECT_LS_FLOW:
    solution = solver::greedy_2mst_localsearch_flow(g, params.depth,
                                                    params.cache, params.trees);
    break;
  case config::Algorithm::MST_ILP:
    convert_links = false;
    solution = solver::greedy_mst_ilp(g, params.sampling);
    break;
  case config::Algorithm::APX1_LN2_E:
    convert_links = false;
    solution = solver::approximate_1_ln2_e(g, params.epsilon);
    break;
  case config::Algorithm::APX1_5_E:
    convert_links = false;
    solution = solver::approximate_1_5_e(g, params.epsilon);
    break;
  case config::Algorithm::SMC: {
    convert_links = true;
    DynamicCactus g_dynamic;
    g_dynamic.read_from_file(params.cactus);
    g_dynamic.copy_links(g);
    solution = watanabe::smc(g_dynamic);
  } break;
  case config::Algorithm::FSM: {
    convert_links = true;
    DynamicCactus g_dynamic;
    g_dynamic.read_from_file(params.cactus);
    g_dynamic.copy_links(g);
    solution = watanabe::fsm(g_dynamic);
  } break;
  case config::Algorithm::HBD: {
    convert_links = true;
    DynamicCactus g_dynamic;
    g_dynamic.read_from_file(params.cactus);
    g_dynamic.copy_links(g);
    solution = watanabe::hbd(g_dynamic);
  } break;
  default:
    throw std::invalid_argument("Algorithm not covered");
  }

  // convert to original edge
  if (convert_links) {
    for (auto &x : solution) {
      auto &orig = g.cactus.links[x.first][x.second];
      if (orig.first == 0 || orig.second == 0) {
        WARN("Could not convert link "
             << x.first << "-" << x.second
             << "to the original graph, this may be caused by a bug");
      }
      x.first = orig.first;
      x.second = orig.second;
    }
  }

  if (params.output.has_value()) {
    for (auto x : solution) {
      g.original_graph.add_link(x);
    }
    g.original_graph.write_to_file(*params.output);
  }

  INFO("Found solution with " << solution.size() << " links");
  INFO("Augmentation weight: " << std::accumulate(
           solution.begin(), solution.end(), 0.,
           [](double v, graph::Edge &e) { return v + e.weight; }));

  return 0;
}
