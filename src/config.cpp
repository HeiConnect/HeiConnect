#ifndef CONFIG_HPP
#define CONFIG_HPP

#include "config.hpp"
#include "argtable3.h"
#include "util.hpp"
#include <cstddef>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <vector>

namespace config {

struct AlgorithmParam {
  Algorithm algorithm;
  std::string name;
  bool public_api;
  std::string description;
  std::vector<std::string> params;
};

static const std::vector<AlgorithmParam> algorithm_params = {
    {Algorithm::GREEDY, "greedy", false, "", {}},
    {Algorithm::GREEDY_GLOBAL, "greedy-global", false, "", {}},
    {Algorithm::HEURISTIC, "heuristic", false, "", {}},
    {Algorithm::GWC, "gwc", true, "GreedyWeightCoverage", {"sampling"}},
    {Algorithm::GWC_SAMPLING, "gwc-sampling", false, "", {}},
    {Algorithm::EILP, "eilp", true, "eILP", {"use-initial"}},
    {Algorithm::APX2_LP, "apx2lp", false, "", {}},
    {Algorithm::MST, "mst", false, "", {}},
    {Algorithm::MST_ILP, "mst-connect-ilp", false, "", {}},
    {Algorithm::APX1_LN2_E, "apx1ln2e", false, "", {}},
    {Algorithm::APX1_5_E, "apx1.5e", false, "", {}},
    {Algorithm::MST_CONNECT, "mst-connect", true, "MSTConnect", {}},
    {Algorithm::FULL_MST, "full-mst", false, "", {}},
    {Algorithm::MST_CONNECT_LS, "mst-ls", false, "", {}},
    {Algorithm::MST_CONNECT_LS_FLOW,
     "mst-connect-ls",
     true,
     "MSTConnect with Local Search",
     {"depth", "cache", "trees"}},
    {Algorithm::SMC, "smc", false, "", {}},
    {Algorithm::FSM, "fsm", false, "", {}},
    {Algorithm::HBD, "hbd", false, "", {}},
    {Algorithm::MST_CONNECT_HEURISTIC, "mst-heuristic", false, "", {}},
    {Algorithm::MST_CONNECT_ORDER_HEURISTIC,
     "mst-order-heuristic",
     false,
     "",
     {}},
};

Algorithm parse_algorithm(const std::string &arg) {
  for (auto &algo_param : algorithm_params) {
    if (arg == algo_param.name) {
      return algo_param.algorithm;
    }
  }
  throw std::invalid_argument("Algorithm is invalid");
}

int parse_distribution(const std::string &arg) {
  if (arg == "unit")
    return 0;
  if (arg == "uniform")
    return 1;
  if (arg == "normal")
    return 2;
  if (arg == "w2")
    return 3;
  if (arg == "w9")
    return 4;
  if (arg == "w99")
    return 5;
  if (arg == "w999")
    return 6;
  if (arg == "u100k")
    return 7;
  else
    throw std::invalid_argument("Invalid distribution");
}

Params::Params(int argc, char **argv) {
  struct arg_lit *help, *arg_cache, *arg_use_initial, *arg_verbose;
  struct arg_str *arg_algorithm, *arg_distribution;
  struct arg_file *arg_original_graph, *arg_cactus, *arg_output_file,
      *arg_link_file;
  struct arg_int *arg_cut, *arg_sampling, *arg_count, *arg_depth, *arg_seed,
      *arg_trees;
  struct arg_dbl *arg_e, *arg_link_fraction;
  struct arg_end *end;

  void *argtable[] = {
      help = arg_litn("h", "help", 0, 1, "display this help and exit"),
      arg_original_graph = arg_file1("g", "graph", "<graph>",
                                     "original graph (in METIS format)"),
      arg_cactus = arg_file1("c", "cactus", "<cactus>",
                             "cactus graph (GraphML xml format)"),
      arg_link_file = arg_file0(
          "l", "links", "<link graph>",
          "path to link file (one link per line, 'source target weight')"),
      arg_output_file = arg_file0("o", "output", "<output graph>",
                                  "file to write augmented graph to"),
      arg_cut = arg_int0(NULL, "cut", "<int>", "the min cut of the graph"),
      arg_algorithm =
          arg_str1("a", "algorithm", "<algorithm>", "Algorithm to use"),
      arg_sampling =
          arg_int0(NULL, "sampling", "<int>", "Divide cuts into n sample sets"),
      arg_count = arg_int0(NULL, "count", "<int>",
                           "Set a count which can be used by some algorithms"),
      arg_cache = arg_lit0(NULL, "cache", "Enable caching if supported"),
      arg_use_initial =
          arg_lit0(NULL, "use-initial", "Use initial solution for ILP"),
      arg_e = arg_dbl0("e", "epsilon", "<double>",
                       "Set epsilon for approximations. Defaults to 0.25"),
      arg_seed = arg_int0("s", "seed", "<int>",
                          "Set the seed for random graph generation"),
      arg_link_fraction =
          arg_dbl0(NULL, "link-fraction", "<double>",
                   "Fraction of links that should be generated"),
      arg_distribution = arg_str0("d", "distribution", "<unit|uniform|normal>",
                                  "Distribution for link weight generation"),
      arg_trees = arg_int0(NULL, "trees", "<int>", "Number of MSTs to use"),
      arg_depth =
          arg_int0(NULL, "depth", "<int>", "Search depth of the algorithm"),
      arg_verbose = arg_lit0("v", "verbose",
#ifdef DBG
                             "Enable verbose output"
#else
                             "Verbose output not available because program was "
                             "not compiled with DBG enabled"
#endif
                             ),
      end = arg_end(20),
  };

  int nerrors;
  nerrors = arg_parse(argc, argv, argtable);

  /* special case: '--help' takes precedence over error reporting */
  if (help->count > 0) {
    std::cout << "Usage: solver";
    arg_print_syntax(stdout, argtable, "\n");
    arg_print_glossary(stdout, argtable, "  %-28s %s\n\n");
    std::cout << "Possible algorithms are:" << std::endl;
    for (auto &algo_param : algorithm_params) {
      if (algo_param.public_api) {
        std::cout << "\t" << algo_param.name << " (" << algo_param.description;
        if (algo_param.params.size()) {
          std::cout << ", available parameters: ";
          for (auto &p : algo_param.params) {
            std::cout << p << ", ";
          }
          std::cout << "\b\b";
        }
        std::cout << ")" << std::endl;
      }
    }
    arg_freetable(argtable, sizeof(argtable) / sizeof(argtable[0]));
    exit(0);
  }

  /* If the parser returned any errors then display them and exit */
  if (nerrors > 0) {
    /* Display the error details contained in the arg_end struct.*/
    arg_print_errors(stdout, end, argv[0]);
    std::cout << "Try " << argv[0] << " --help' for more information."
              << std::endl;
    arg_freetable(argtable, sizeof(argtable) / sizeof(argtable[0]));
    exit(1);
  }

#ifdef DBG
  set_verbose(arg_verbose->count);
#else
  if (arg_verbose->count)
    WARN("Verbose logging is not available because DBG was not enabled at "
         "compile time");
#endif

  original_graph = std::filesystem::path{arg_original_graph->filename[0]};
  if (arg_link_file->count) {
    links = arg_link_file->filename[0];
  } else {
    links = {original_graph};
    links.replace_extension("links");
  }
  cactus = std::filesystem::path{arg_cactus->filename[0]};
  if (arg_output_file->count) {
    output = std::optional<std::filesystem::path>{arg_output_file->filename[0]};
  }
  if (arg_cut->count > 0) {
    cut = arg_cut->ival[0];
  }

  if (arg_algorithm->count) {
    algorithm = parse_algorithm(arg_algorithm->sval[0]);
  }

  if (arg_sampling->count) {
    sampling = arg_sampling->ival[0];
  }
  if (arg_count->count) {
    count = arg_count->ival[0];
  }
  if (arg_cache->count) {
    cache = true;
  }
  if (arg_use_initial->count) {
    use_initial = true;
  }
  if (arg_e->count) {
    epsilon = arg_e->dval[0];
  }
  if (arg_seed->count) {
    seed = arg_seed->ival[0];
  }
  if (arg_link_fraction->count) {
    link_fraction = arg_link_fraction->dval[0];
  }
  if (arg_trees->count) {
    trees = arg_trees->ival[0];
  }

  if (arg_distribution->count) {
    distribution = parse_distribution(arg_distribution->sval[0]);
  }

  if (arg_depth->count) {
    depth = arg_depth->ival[0];
  }
} // namespace config

} // namespace config

#endif // CONFIG_HPP
