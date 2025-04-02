#pragma once

#include "ConstrSimple.hpp"
#include "IntSet.hpp"
#include "Solver.hpp"
#include "auxiliary.hpp"
#include "coreguided.hpp"
#include "propgraph.hpp"
#include "run.hpp"
#include "typedefs.hpp"

namespace rs {
namespace run {

using rs::operator<<;
using CardE = Ce32;
using CardS = ConstrSimple32;

template <typename SMALL, typename LARGE>
class Optimization {
 public:
  using CE = CePtr<ConstrExp<SMALL, LARGE>>;
  using CS = ConstrSimple<SMALL, LARGE>;
  using LazyEnc = cg::LazyEnc<SMALL, LARGE>;
  using SumEnc = cg::SumEnc<SMALL, LARGE>;

 private:
  CE origObj;
  CE reformObj;

  std::vector<Lit> bestSol{};

  // Propagation graph of objective literals
  PropGraph pgraph;

  LARGE hard_lower_bound;
  LARGE upper_bound;
  LARGE hard_upper_bound;
  ID lastUpperBound = ID_Undef;
  ID lastUpperBoundUnprocessed = ID_Undef;
  ID lastLowerBound = ID_Undef;
  ID lastLowerBoundUnprocessed = ID_Undef;

  // enforce this to switch upper bounds from `<= best` to `< best`
  // 0 means not yet defined
  Var equalitySwitch = 0;

  bool ubStrict = true;

  // Cores extracted during one WCE phase
  std::vector<std::pair<CardS, SMALL>> wceCoreBuf;
  // Lazily reformulated cores
  std::vector<LazyEnc> lazies;
  IntSet assumps;

  // Check if lazily abstractions need to be extended
  bool checkLazyVariables();
  // Core trimming
  bool trimCore(CeSuper core) {
    if (!options.cgTrim) return false;
    return cg::trimCore(core);
  }
  // Core minimization
  bool minimizeCore(CeSuper core) {
    if (!options.cgMinimize) return false;
    return cg::minimizeCore(core);
  }
  bool addLowerBound(bool force = false);
  bool addUpperBound(bool softened);
  // End WCE phase. To be called after getting SAT.
  bool finalizeWcePhase() {
    return cg::finalizeWcePhase<SMALL, LARGE>(
        wceCoreBuf, reformObj, lazies, [this]() -> LARGE { return normalizedUpperBound(); },
        [this](LazyEnc& enc) -> bool {
          if (!options.cgExhaust) return false;
          assert(!options.cgEncoding.is("sum"));
          bool reified = options.cgEncoding.is("reified");
          IntSet assumps{};
          if (ubStrict && equalitySwitch) assumps.add(equalitySwitch);
          int propLim = options.cgExPropLim.get();
          int confLim = options.cgExConfLim.get();
          return enc.exhaust(
              reformObj, reified, Limits(propLim < 0 ? Limits::UNLIM : propLim, confLim < 0 ? Limits::UNLIM : confLim),
              assumps,
              [this](const std::vector<Lit>& sol) {
                return handleNewSolution(sol, !options.cgUbConstr.is("no"), !options.cgUbConstr.is("<"));
              },
              [this] { return printObjBounds(); });
        });
  }
  // Handle newly extracted core(s)
  bool handleInconsistency(std::vector<CeSuper>& cores) {
    return cg::handleInconsistency(cores, wceCoreBuf, reformObj, equalitySwitch);
  };
  // Handle newly discovered solution
  bool handleNewSolution(const std::vector<Lit>& sol, bool addConstr, bool softened);
  // Choose assumptions for next CG call
  SMALL populateAssumptions(SMALL coeflim, bool refine) {
    assumps.clear();
    return cg::populateAssumptions(assumps, reformObj, coeflim, refine, bestSol);
  }
  // Update the propagation graph
  bool updatePropGraph();
  // Run at-most-one detection, returns true if unsat
  bool atMostOnes();
  // Build a PB reformulation from a cardinality core
  CS buildPbReform(const CardS& core, SMALL& mult, LARGE& upperBound) const {
    return cg::buildPbReform(core, reformObj, mult, upperBound);
  }

  inline LARGE lowerBound() const { return -reformObj->getDegree(); }

  bool harden() {
    if (!options.cgHardening) return false;
    LARGE diff = upper_bound - lowerBound();
    for (Var v : reformObj->vars) {
      if (aux::abs(reformObj->coefs[v]) > diff) {
        if (solver.addUnitConstraint(-reformObj->getLit(v), Origin::HARDENEDBOUND).second == ID_Unsat) {
          return true;
        }
      }
    }
    reformObj->removeUnitsAndZeroes(solver, false);
    return false;
  }

 public:
  Optimization(CeArb obj) {
    origObj = solver.cePools.take<SMALL, LARGE>();
    obj->copyTo(origObj);
    assert(origObj->vars.size() > 0);
    // NOTE: -origObj->getDegree() keeps track of the offset of the reformulated objective (or after removing unit lits)
    upper_bound = origObj->absCoeffSum() - origObj->getRhs() + 1;
    hard_upper_bound = upper_bound;

    reformObj = solver.cePools.take<SMALL, LARGE>();
    reformObj->stopLogging();
    if (!options.optMode.is("linear")) {
      origObj->copyTo(reformObj);
    }

    hard_lower_bound = lowerBound();
  };

  inline LARGE normalizedLowerBound() { return lowerBound() + origObj->getDegree(); }
  inline LARGE normalizedUpperBound() { return upper_bound + origObj->getDegree(); }

  void printObjBounds() const {
    if (options.verbosity.get() == 0) return;
    std::cout << "c bounds ";
    if (solver.foundSolution()) {
      std::cout << upper_bound;
    } else {
      std::cout << "-";
    }
    std::cout << " >= " << lowerBound() << " @ " << stats.getTime() << "\n";
  }

  void logProof();

  void optimize();
};

}  // namespace run

}  // namespace rs
