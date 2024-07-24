#ifndef VIECUT_RUNNER_HPP
#define VIECUT_RUNNER_HPP

#include "util.hpp"

#include <cstdlib>
#include <exception>
#include <filesystem>
#include <iostream>
#include <limits.h>
#include <sstream>
#include <unistd.h>

inline std::filesystem::path executable_path() {
  char result[PATH_MAX];
  ssize_t count = readlink("/proc/self/exe", result, PATH_MAX);
  return {std::string(result, (count > 0) ? count : 0)};
}

inline void compute_cactus(std::filesystem::path graph,
                           std::filesystem::path cactus) {
  if (std::filesystem::exists(cactus)) {
    INFO("Found cactus file " << cactus);
    return;
  }

  INFO("Generating missing cactus file " << cactus);
  std::filesystem::path repo_root =
      executable_path().parent_path().parent_path();
  std::filesystem::path viecut =
      repo_root / "extern" / "VieCut" / "build" / "mincut";

  if (!std::filesystem::exists(viecut)) {
    std::stringstream ss;
    ss << "VieCut binary " << viecut
       << " does not exist, make sure VieCut submodule is built successfully";
    throw std::runtime_error(ss.str());
  }

  std::stringstream command;
  command << viecut.string() << " " << graph.string() << " cactus -t "
          << cactus.string();
  system(command.str().c_str());

  INFO("Cactus generated");
}

#endif // VIECUT_RUNNER_HPP
