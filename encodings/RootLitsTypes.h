#ifndef ROOTLITSTYPES_H
#define ROOTLITSTYPES_H
#include "core/SolverTypes.h"
#include <utility>
#include <vector>

namespace rootLits{
  using myType = std::vector<std::pair<uint64_t, Glucose::Lit>>;
  using value_t =std::pair<uint64_t, Glucose::Lit>; 
}

#endif
