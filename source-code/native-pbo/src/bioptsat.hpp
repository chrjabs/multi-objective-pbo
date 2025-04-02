#pragma once

#include "coreguided.hpp"
#include "multiopt.hpp"

namespace rs {
namespace run {

template <typename CF, typename DG>
class BiOptSat : public MultiOpt<CF, DG> {
  using CE = CePtr<ConstrExp<CF, DG>>;
  using CS = ConstrSimple<CF, DG>;
  using ND = NonDom<DG>;
  using RObj = ReifiedObj<CF, DG>;
  using MO = MultiOpt<CF, DG>;

 protected:
  std::vector<RObj> rObjectives;
  std::vector<cg::LazyEnc<CF, DG>> ollLazies;

  ID incrLb = ID_Undef;

  inline bool minIncr(std::vector<Lit> &solution, std::vector<DG> &costs) {
    if (options.bosVariant.is("sat-unsat")) return minIncrSolImpr(solution, costs);
    if (options.bosVariant.is("oll")) return minIncrOll(solution, costs);
    assert(false);
    return false;
  }

  bool minDecr(std::vector<Lit> &solution, std::vector<DG> &costs);

  bool minIncrSolImpr(std::vector<Lit> &solution, std::vector<DG> &costs);
  bool minIncrOll(std::vector<Lit> &solution, std::vector<DG> &costs);

  bool addDecConstr(const std::vector<DG> &costs, std::vector<Lit> &witness);

 public:
  virtual void solve() override;

  BiOptSat(const std::vector<CeArb> objs) : MO(objs), rObjectives{RObj(MO::reformObjs[0]), RObj(MO::reformObjs[1])} {
    assert(objs.size() == 2);
  }
};

}  // namespace run
}  // namespace rs
