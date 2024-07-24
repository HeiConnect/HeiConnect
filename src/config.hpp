#include <filesystem>
#include <optional>

namespace config {

enum class Algorithm {
  GREEDY,
  GREEDY_GLOBAL,
  HEURISTIC,
  GWC,
  GWC_SAMPLING,
  EILP,
  APX2_LP,
  MST,
  MST_ILP,
  APX1_LN2_E,
  APX1_5_E,
  MST_CONNECT,
  FULL_MST,
  MST_CONNECT_LS,
  MST_CONNECT_LS_FLOW,
  SMC,
  FSM,
  HBD,
  MST_CONNECT_HEURISTIC,
  MST_CONNECT_ORDER_HEURISTIC,
};

enum class Distribution {
  UNIT,
  UNIFORM,
  NORMAL,
  W2,
  W9,
  W99,
  W999,
  U100K,
};

class Params {
public:
  std::filesystem::path original_graph;
  std::filesystem::path links;
  std::filesystem::path cactus;
  std::optional<std::filesystem::path> output = {};

  int cut = -1;
  Algorithm algorithm;
  int sampling = 0;
  int trees = 0;
  int count = 1;
  int depth = 0;
  bool cache = false;
  bool use_initial = false;
  double epsilon = 0.25;
  int seed = -1;
  double link_fraction = 1.;
  int distribution = 1;
  bool output_links = false;

  Params(int argc, char **argv);
};

} // namespace config
