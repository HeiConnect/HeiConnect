#ifndef UTIL_HPP
#define UTIL_HPP

#include <algorithm>
#include <chrono>
#include <iostream>
#include <iterator>
#include <list>
#include <numeric>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#ifdef DBG
void set_verbose(bool);
bool is_verbose();
#define DEBUG(msg)                                                             \
  if (is_verbose())                                                            \
    std::cout << "DEBUG: " << msg << std::endl;
#else
#define DEBUG(msg) ;
#endif
#define INFO(msg) std::cout << " INFO: " << msg << std::endl;
#define WARN(msg) std::cout << " WARN: " << msg << std::endl;
#define ERROR(msg) std::cerr << "ERROR: " << msg << std::endl;
#define FATAL(msg)                                                             \
  std::cerr << "ERROR:" << msg << std::endl;                                   \
  exit(1);

#define COUNT_SIZE(c)                                                          \
  (std::accumulate(c.begin(), c.end(), 0, [](int cnt, auto &container) {       \
    return cnt + container.size();                                             \
  }))

inline std::vector<int> split(const std::string &s, char delim) {
  std::vector<int> result;
  std::stringstream ss(s);
  std::string item;

  while (getline(ss, item, delim)) {
    result.push_back(std::stoi(item));
  }

  return result;
}

template <typename Edge> void print_solution(std::list<Edge> solution) {
  INFO("Added " << solution.size() << " edges");
  for (auto &x : solution) {
    DEBUG("Link added: " << x.first << "-" << x.second << ": " << x.weight);
  }
  INFO("Augmentation weight: "
       << std::accumulate(solution.begin(), solution.end(), 0.,
                          [](double v, Edge &e) { return v + e.weight; }));
}

#define TIMEIT(msg, cmd) timeit([&]() { cmd; }, msg);

template <typename F>
inline void timeit(F &&functor, std::string message = "Function") {
  auto begin = std::chrono::steady_clock::now();
  functor();
  auto end = std::chrono::steady_clock::now();

  auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin)
                .count();
  INFO(message << " took " << ms << " ms");
}

// returns all non-empty subsets of size up to alpha of {0,...,n-1}
inline std::list<std::vector<int>> enumerate(int n, int alpha) {
  alpha = std::min(alpha, n);
  std::list<std::vector<int>> result = {{}};

  std::vector<int> subset(alpha, -1);
  subset.back() = 0;

  int min_idx = alpha - 1;
  while (true) {
    // Add current combination
    std::vector<int> tmp;
    std::copy(subset.begin() + min_idx, subset.end(), std::back_inserter(tmp));
    result.push_back(tmp);

    // Find, if existing, the rightmost index that can be incremented
    int idx = alpha - 1;
    while (idx >= 0 && subset[idx] == n - alpha + idx) {
      idx--;
    }
    if (idx < 0) {
      break;
    }
    if (idx < min_idx) {
      min_idx = idx;
    }

    // Increment the found index and update the rest
    subset[idx]++;
    for (int j = idx + 1; j < alpha; ++j) {
      subset[j] = subset[j - 1] + 1;
    }
  }
  return result;
}

template <typename T>
std::vector<T> get_indices(std::vector<T> &data, std::vector<int> &indices) {
  std::vector<T> list;
  for (int idx : indices) {
    list.push_back(data[idx]);
  }
  return list;
}

template <template <typename...> class Container, typename T>
Container<T> set_union(Container<T> a, Container<T> b) {
  std::set<T> set;
  for (T &e : a)
    set.insert(e);
  for (T &e : b)
    set.insert(e);
  Container<T> result;
  std::copy(set.begin(), set.end(), std::back_inserter(result));
  return result;
}

#endif // UTIL_HPP
