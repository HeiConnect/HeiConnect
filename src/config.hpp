#include <filesystem>
#include <optional>

namespace config {

enum class Algorithm {
  GREEDY,
  GREEDY_STRONG,
  HEURISTIC,
  DYNAMIC,
  DYNAMIC_BOUNDED,
  HEURISTIC_SAMPLING,
  ILP,
  APX2_LP,
  MST,
  MST_ILP,
  APX1_LN2_E,
  APX1_5_E,
  MST_FLOW,
  FULL_MST,
  MST_LS,
  MST_LS_FLOW,
  SMC,
  FSM,
  HBD,
  MST_FLOW_HEURISTIC,
  MST_FLOW_ORDER_HEURISTIC,
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
  std::optional<std::filesystem::path> output = std::nullopt;

  int cut = -1;
  Algorithm algorithm;
  int sampling = 0;
  int trees = 1;
  int count = 1;
  int depth = 1;
  bool cache = false;
  bool use_initial = false;
  double epsilon = 0.25;
  int seed = -1;
  double link_fraction = 1.;
  int distribution = 1;

  Params(int argc, char **argv);
};

} // namespace config
