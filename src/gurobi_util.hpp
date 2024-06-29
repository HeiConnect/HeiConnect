#ifndef GUROBI_UTIL_HPP
#define GUROBI_UTIL_HPP

#include "util.hpp"

#include <gurobi_c++.h>

inline void solve_lp(GRBModel &model) {
  // solve
  TIMEIT("Optimizing model", model.optimize());
  // Check that the problem has an optimal solution.
  int status = model.get(GRB_IntAttr_Status);
  if (status == GRB_INF_OR_UNBD || status == GRB_INFEASIBLE ||
      status == GRB_UNBOUNDED) {
    FATAL("The model cannot be solved because it is infeasible or unbounded");
  }
  if (status != GRB_OPTIMAL) {
    FATAL("Optimization was stopped with status " << status);
  }
}

#endif // GUROBI_UTIL_HPP
