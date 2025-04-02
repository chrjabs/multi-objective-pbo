#pragma once

#include <memory>
#include "auxiliary.hpp"
#include "run.hpp"

namespace rs {
namespace cg {
template <typename DG>
using UbCb = const std::function<DG()>&;
using SolCb = const std::function<bool(const std::vector<Lit>&)>&;
using BoundCb = const std::function<void()>&;
template <typename CF>
using WCores = std::vector<std::pair<ConstrSimple32, CF>>;

// Lazy core reformulation
template <typename CF, typename DG>
struct LazyEnc {
  using CS = ConstrSimple<CF, DG>;
  using CE = CePtr<ConstrExp<CF, DG>>;

  // note: only valid if core is non-card
  std::unique_ptr<aux::subsum::IncComp<CF, DG>> subsumComp{};
  CS atLeast{};  // X >= k + y1 + ... + yi
  CS atMost{};   // k + y1 + ... + yi-1 + (1+n-k-i)yi >= X
  DG coveredVal;
  DG maxBound;
  CF mult;
  Var watchedVar;
  ID atLeastID = ID_Undef;
  ID atMostID = ID_Undef;
  enum Flags { ATLEAST = 0b1, ATMOST = 0b10, SYM = 0b100, FINAL = 0b1000, PBREF = 0b10000 };
  char flags = 0;  // pbref|finalized|sym|atmost|atleast

  LazyEnc(const CS& core, DG maxBound, CF mult);
  ~LazyEnc() {
    if (atLeastID != ID_Unsat) run::solver.dropExternal(atLeastID, false, false);
    if (atMostID != ID_Unsat) run::solver.dropExternal(atMostID, false, false);
    if (subsumComp) subsumComp.reset();
  }

  // move constructor
  LazyEnc(LazyEnc&& other) noexcept {
    subsumComp = std::move(other.subsumComp);
    atLeast = std::move(other.atLeast);
    atMost = std::move(other.atMost);
    coveredVal = std::move(other.coveredVal);
    maxBound = std::move(other.maxBound);
    mult = std::move(other.mult);
    watchedVar = other.watchedVar;
    atLeastID = other.atLeastID;
    atMostID = other.atMostID;
    flags = other.flags;
    // Remove IDs from old object so that the deconstructor does not try to free them
    other.atLeastID = ID_Undef;
    other.atMostID = ID_Undef;
  };
  LazyEnc& operator=(LazyEnc&& other) noexcept {
    subsumComp = std::move(other.subsumComp);
    atLeast = std::move(other.atLeast);
    atMost = std::move(other.atMost);
    coveredVal = std::move(other.coveredVal);
    maxBound = std::move(other.maxBound);
    mult = std::move(other.mult);
    watchedVar = other.watchedVar;
    atLeastID = other.atLeastID;
    atMostID = other.atMostID;
    flags = other.flags;
    // Remove IDs from old object so that the deconstructor does not try to free them
    other.atLeastID = ID_Undef;
    other.atMostID = ID_Undef;
    return *this;
  };

  inline bool isInSolver() const {
    return ((!options.cgRefLb || (flags & ATLEAST)) && (flags & ATMOST) && (!options.cgRefSym || (flags & SYM))) ||
           (flags & FINAL);
  }
  CF addVar(Var v, bool reified);
  bool addConstraints(bool reified, bool addAtLeast, Var prevvar);

  // Don't use these three directly
  ID addAtLeastConstraint(bool reified);
  ID addAtMostConstraint(bool reified);
  // Output implications: o_i -> o_{i-1}
  ID addSymBreakingConstraint(Var prevvar);

  bool isFinal() const { return isInSolver() && maxBound <= coveredVal; }
  void decreaseMaxBound(DG cardMaxBound) {
    assert(maxBound >= cardMaxBound);
    maxBound = cardMaxBound;
  }
  bool exhaust(CE reformObj, bool reified, Limits lims, const IntSet baseAssumps, SolCb solCb, BoundCb boundCb);
};

template <typename CF, typename DG>
struct SumEnc {
  using CS = ConstrSimple<CF, DG>;

  CS reform{};
  int firstAuxIdx;

  SumEnc() {}
  SumEnc(const CS& core, std::pair<Var, Var> auxRange);
  SumEnc(const CS& core, std::vector<Term<CF>> auxs);

  SumEnc(SumEnc&&) = default;
  SumEnc& operator=(SumEnc&&) = default;

  bool addFull();
};

template <typename CF, typename DG>
using ExhaustCb = const std::function<bool(LazyEnc<CF, DG>&)>&;

template <typename CF, typename DG>
using Lazies = std::vector<LazyEnc<CF, DG>>;

bool trimCore(CeSuper core);
bool minimizeCore(CeSuper core);

template <typename CF, typename DG>
CF populateAssumptions(IntSet& assumps, CePtr<ConstrExp<CF, DG>> reformObj, CF coeflim, bool refine,
                       const std::vector<Lit>& bestSol);

template <typename CF, typename DG>
bool finalizeWcePhase(WCores<CF>& cores, CePtr<ConstrExp<CF, DG>> reformObj, Lazies<CF, DG>& lazies, UbCb<DG> ubCb,
                      ExhaustCb<CF, DG> exCb);

// Build a PB reformulation from a cardinality core
template <typename CF, typename DG>
ConstrSimple<CF, DG> buildPbReform(const ConstrSimple32& core, const CePtr<ConstrExp<CF, DG>> reformObj, CF& mult,
                                   DG& upperBound);

// Round found core to a cardinality core
template <typename CF, typename DG>
std::pair<Ce32, int> reduceToCardinality(const CeSuper& core, const CePtr<ConstrExp<CF, DG>> reformObj);

template <typename CF, typename DG>
bool handleInconsistency(std::vector<CeSuper>& origCores, WCores<CF>& coreBuf, CePtr<ConstrExp<CF, DG>> reformObj,
                         Lit equalitySwitch);

template <typename CF, typename DG>
bool coreGuidedReformulate(CePtr<ConstrExp<CF, DG>> reformObj, const ConstrSimple<CF, DG>& origObj,
                           std::vector<Lit>& optSol);

template <typename CF, typename DG>
bool coreGuidedReformulateWithLazies(CePtr<ConstrExp<CF, DG>> reformObj, const ConstrSimple<CF, DG>& origObj,
                                     std::vector<Lit>& optSol, Lazies<CF, DG>& lazies);

template <typename CF, typename DG>
bool harden(CePtr<ConstrExp<CF, DG>> reformObj, DG upper_bound) {
  if (!options.cgHardening) return false;
  DG diff = upper_bound + reformObj->getDegree();
  for (Var v : reformObj->vars) {
    if (aux::abs(reformObj->coefs[v]) > diff) {
      if (run::solver.addUnitConstraint(-reformObj->getLit(v), Origin::HARDENEDBOUND).second == ID_Unsat) {
        return true;
      }
    }
  }
  reformObj->removeUnitsAndZeroes(run::solver, false);
  return false;
}

template <typename CF, typename DG>
bool checkLazyVariables(CePtr<ConstrExp<CF, DG>> reformObj, std::vector<LazyEnc<CF, DG>>& lazies) {
  if (options.cgEncoding.is("sum")) return false;

  bool reified = options.cgEncoding.is("reified");
  for (int i = 0; i < (int)lazies.size(); ++i) {
    LazyEnc<CF, DG>& abstraction = lazies[i];
    if (!abstraction.isInSolver()) continue;

    if (reformObj->getLit(abstraction.watchedVar) == 0) {
      if (abstraction.isFinal() ||
          run::solver.isUnit(-abstraction.watchedVar)) {  // binary constraints make all new auxiliary variables unit
        aux::swapErase(lazies, i--);                      // fully expanded, no need to keep in memory
      } else {                                            // add auxiliary variable
        int newN = run::solver.getNbVars() + 1;
        run::solver.setNbVars(newN);
        Var oldvar = abstraction.watchedVar;
        CF w = abstraction.addVar(newN, reified);
        // reformulate the objective
        reformObj->addLhs(w, newN);
        // add necessary lazy constraints
        if (abstraction.addConstraints(reified, options.cgRefLb.get(), oldvar)) return true;
        if (abstraction.isFinal()) aux::swapErase(lazies, i--);
      }
    }
  }
  return false;
}

template <typename CF, typename DG>
DG computeCost(const ConstrSimple<CF, DG>& obj, const std::vector<Lit>& sol) {
  DG val = -obj.rhs;
  for (auto t : obj.terms) val += t.c * (int)(sol[toVar(t.l)] == t.l);
  return val;
}

template <typename CF, typename DG>
DG computeCost(const CePtr<ConstrExp<CF, DG>>& obj, const std::vector<Lit>& sol) {
  DG val = -obj->getRhs();
  for (Var v : obj->vars) val += obj->getCoef(v) * (int)(sol[v] == v);
  return val;
}

template <typename CF, typename DG>
DG maxObjVal(const CePtr<ConstrExp<CF, DG>> obj) {
  DG result = -obj->getRhs();
  for (Var v : obj->vars)
    if (obj->coefs[v] > 0) result += obj->coefs[v];
  return result;
}
}  // namespace cg
}  // namespace rs
