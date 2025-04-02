#pragma once

#include "multiopt.hpp"

namespace rs {
namespace run {

template <typename CF, typename DG>
class PMinimal : public MultiOpt<CF, DG> {
  using CE = CePtr<ConstrExp<CF, DG>>;
  using CS = ConstrSimple<CF, DG>;
  using ND = NonDom<DG>;
  using RObj = ReifiedObj<CF, DG>;
  using MO = MultiOpt<CF, DG>;

 protected:
  std::vector<RObj> rObjectives;

  // Performs P-minimization based of a solution. The solution and cost parameters are both an input of a starting point
  // and an output of the minimized solution. If UNSAT was encoutered, the method returns true.
  bool pMinimize(std::vector<Lit> &solution, std::vector<DG> &costs);
  // Generates assumptions that enforce the next solution to be weakly dominating the given costs
  bool enforceDominating(const std::vector<DG> &costs, IntSet &assumps);

 public:
  virtual void solve() override;

  PMinimal(const std::vector<CeArb> objs) : MO(objs) {
    rObjectives.reserve(MultiOpt<CF, DG>::reformObjs.size());
    for (const auto &obj : MO::reformObjs) rObjectives.emplace_back(obj);
  }
};

}  // namespace run
}  // namespace rs
