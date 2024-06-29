#include "ilp.hpp"
#include "gurobi_util.hpp"
#include "kruskal_mst.hpp"
#include "util.hpp"

#include <gurobi_c++.h>

namespace solver {

#define LINK_THRESHOLD 2000.

std::list<graph::Edge> greedy_mst_i(graph::GraphPair &g) {
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

std::list<graph::Edge> ilp(graph::GraphPair &graph, bool use_initial,
                           int presolve) {
  try {
    std::set<std::pair<int, int>> mst_edges;
    if (use_initial) {
      // compute initial heuristic solution

      // Create graph from links in cactus graph
      {
        auto heuristic_solution = greedy_mst_i(graph);
        for (auto &link : heuristic_solution) {
          std::pair<int, int> e;
          if (link.first < link.second)
            e = {link.first, link.second};
          else
            e = {link.second, link.first};
          mst_edges.insert(e);
        }
      }
    }

    INFO("Creating ILP");
    std::list<graph::CactusCut> cuts = graph.get_min_cuts();

    // Create ILP model
    GRBEnv env = GRBEnv(true);
    env.set("LogFile", "connectivity_augmentation.log");
    env.start();
    GRBModel model = GRBModel(env);
    model.getEnv().set(GRB_IntParam_Presolve, presolve);

    // variables
    std::vector<std::vector<GRBVar>> vars;
    vars.resize(graph.cactus.num_nodes() + 1, {});
    for (int i = 1; i <= graph.cactus.num_nodes(); ++i) {
      vars[i].resize(graph.cactus.num_nodes() + 1, {});
    }
    for (int i = 1; i <= graph.cactus.num_nodes(); ++i) {
      for (int j = i + 1; j <= graph.cactus.num_nodes(); ++j) {
        if (graph.cactus.links[i][j].first != 0 &&
            graph.cactus.links[i][j].weight < LINK_THRESHOLD) {
          vars[i][j] =
              model.addVar(0.0, 1.0, 0.0, GRB_BINARY,
                           "e" + std::to_string(i) + "-" + std::to_string(j));
          if (use_initial) {
            // set mst as initial solution
            bool is_initial_solution =
                mst_edges.find({i, j}) != mst_edges.end();
            vars[i][j].set(GRB_DoubleAttr_Start, is_initial_solution);
          }
        }
      }
    }
    model.update();
    INFO("Number of variables = " << model.get(GRB_IntAttr_NumBinVars));

    // constraints
    std::vector<GRBConstr> con_cut;
    con_cut.reserve(cuts.size());
    unsigned int count = 1;
    for (auto &cut : cuts) {
      GRBLinExpr expr = 0;
      for (unsigned int i = 1; i <= graph.cactus.num_nodes(); ++i) {
        for (unsigned int j = i + 1; j <= graph.cactus.num_nodes(); ++j) {
          if (graph.cactus.links[i][j].first != 0 &&
              graph.cactus.links[i][j].weight < LINK_THRESHOLD) {
            if (cut.is_edge_cut(i, j)) {
              expr += vars[i][j];
            }
          }
        }
      }
      GRBConstr con = model.addConstr(expr >= 1, "c" + std::to_string(count++));
      // con.set(GRB_IntAttr_Lazy, 1);

      con_cut.push_back(con);
    }
    model.update();
    INFO("Number of constraints = " << model.get(GRB_IntAttr_NumConstrs));

    // objective
    GRBLinExpr obj = 0;
    for (int i = 1; i <= graph.cactus.num_nodes(); ++i) {
      for (int j = i + 1; j <= graph.cactus.num_nodes(); ++j) {
        if (graph.cactus.links[i][j].first != 0 &&
            graph.cactus.links[i][j].weight < LINK_THRESHOLD) {
          obj += vars[i][j] * graph.cactus.links[i][j].weight;
        }
      }
    }
    model.setObjective(obj, GRB_MINIMIZE);
    model.update();

    // solve
    solve_lp(model);

    // get result
    int nSolutions = model.get(GRB_IntAttr_SolCount);
    INFO("Number of solutions found: " << nSolutions);
    double objVal = model.get(GRB_DoubleAttr_ObjVal);
    INFO("Objective: " << objVal);

    std::list<graph::Edge> links;
    for (int i = 1; i <= graph.cactus.num_nodes(); ++i) {
      for (int j = i + 1; j <= graph.cactus.num_nodes(); ++j) {
        if (graph.cactus.links[i][j].first != 0 &&
            graph.cactus.links[i][j].weight < LINK_THRESHOLD) {
          if (vars[i][j].get(GRB_DoubleAttr_X) > 0) {
            links.push_back({graph.cactus.links[i][j].first,
                             graph.cactus.links[i][j].second,
                             graph.cactus.links[i][j].weight});
            INFO("Add edge: " << vars[i][j].get(GRB_StringAttr_VarName) << ", "
                              << graph.cactus.links[i][j].first << " -> "
                              << graph.cactus.links[i][j].second << " (weight "
                              << graph.cactus.links[i][j].weight << ")");
          }
        }
      }
    }
    INFO("Number of links added: '" << links.size() << "'");
    return links;
  } catch (GRBException e) {
    ERROR("GRBException, error code = " << e.getErrorCode());
    WARN(e.getMessage());
    throw e;
  }
}

} // namespace solver
