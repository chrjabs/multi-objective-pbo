#pragma once

#include "coreguided.hpp"
#include "pminimal.hpp"

namespace rs {
namespace run {

template <typename CF, typename DG>
class PareTopK : public PMinimal<CF, DG> {
  using CE = CePtr<ConstrExp<CF, DG>>;
  using CS = ConstrSimple<CF, DG>;
  using ND = NonDom<DG>;
  using RObj = ReifiedObj<CF, DG>;
  using MO = MultiOpt<CF, DG>;
  using PMin = PMinimal<CF, DG>;
  using LazyEnc = cg::LazyEnc<CF, DG>;
  using SumEnc = cg::SumEnc<CF, DG>;

 protected:
  // Joint objective to perform CG search on
  CS origJointObj{};
  CE reformJointObj;
  DG upper_bound;
  CF coeflim = 0;
  // Lazily reformulated cores
  std::vector<LazyEnc> lazies;

 public:
  virtual void solve() override;

  PareTopK(const std::vector<CeArb> objs) : PMin(objs) {
    reformJointObj = solver.cePools.take<CF, DG>();
    for (auto& obj : MO::reformObjs) {
      for (Var v : obj->vars) {
        reformJointObj->addLhs(obj->getCoef(v), v);
      }
    }
    upper_bound = reformJointObj->absCoeffSum();
    coeflim = options.cgStrat ? reformJointObj->getLargestCoef() : 0;
    for (auto& obj : MO::origObjs) {
      for (auto& t : obj.terms) {
        origJointObj.terms.emplace_back(t);
      }
    }
  }
};

}  // namespace run
}  // namespace rs
