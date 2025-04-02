#include "optimize.hpp"
#include <boost/math/tools/traits.hpp>
#include "ConstrExp.hpp"
#include "auxiliary.hpp"

namespace rs {
namespace run {

template <typename SMALL, typename LARGE>
CePtr<ConstrExp<SMALL, LARGE>> expand(const ConstrSimple<SMALL, LARGE>& simp) {
  CePtr<ConstrExp<SMALL, LARGE>> exp = solver.cePools.take<SMALL, LARGE>();
  exp->addRhs(simp.rhs);
  for (const Term<SMALL>& t : simp.terms) exp->addLhs(t.c, t.l);
  exp->orig = simp.orig;
  if (exp->plogger) {
    exp->proofBuffer.str(std::string());
    exp->proofBuffer << simp.proofLine;
  }
  return exp;
}

template <typename SMALL, typename LARGE>
bool Optimization<SMALL, LARGE>::checkLazyVariables() {
  if (options.cgEncoding.is("sum")) return false;

  bool reified = options.cgEncoding.is("reified");
  for (int i = 0; i < (int)lazies.size(); ++i) {
    LazyEnc& abstraction = lazies[i];
    if (!abstraction.isInSolver()) continue;

    if (reformObj->getLit(abstraction.watchedVar) == 0) {
      LARGE coreUpper = abstraction.maxBound;
      if (options.cgCoreUpper) {
        coreUpper = static_cast<int>(std::min<LARGE>(coreUpper, normalizedUpperBound() / abstraction.mult));
        // NOTE: integer division rounds down
      }
      assert(coreUpper >= 0);
      abstraction.decreaseMaxBound(coreUpper);
      if (abstraction.isFinal() ||
          solver.isUnit(-abstraction.watchedVar)) {  // binary constraints make all new auxiliary variables unit
        aux::swapErase(lazies, i--);                 // fully expanded, no need to keep in memory
      } else {                                       // add auxiliary variable
        int newN = solver.getNbVars() + 1;
        solver.setNbVars(newN);
        Var oldvar = abstraction.watchedVar;
        SMALL w = abstraction.addVar(newN, reified);
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

template <typename SMALL, typename LARGE>
bool Optimization<SMALL, LARGE>::addLowerBound(bool force) {
  if (lowerBound() <= hard_lower_bound && !force) return false;
  CE aux = solver.cePools.take<SMALL, LARGE>();
  origObj->copyTo(aux);
  aux->addRhs(lowerBound());
  solver.dropExternal(lastLowerBound, true, true);
  std::pair<ID, ID> res = solver.addConstraint(aux, Origin::LOWERBOUND);
  lastLowerBoundUnprocessed = res.first;
  lastLowerBound = res.second;
  hard_lower_bound = lowerBound();
  return (lastLowerBound == ID_Unsat);
}

template <typename SMALL, typename LARGE>
bool Optimization<SMALL, LARGE>::addUpperBound(bool softened) {
  CE aux = solver.cePools.take<SMALL, LARGE>();
  origObj->copyTo(aux);
  if (softened) {
    if (!equalitySwitch) {
      equalitySwitch = solver.getNbVars() + 1;
      solver.setNbVars(equalitySwitch);
    }
    aux->addLhs(1, equalitySwitch);
  }
  aux->invert();
  aux->addRhs(-upper_bound + (softened ? 0 : 1));
  solver.dropExternal(lastUpperBound, true, true);
  std::pair<ID, ID> res = solver.addConstraint(aux, Origin::UPPERBOUND);
  lastUpperBoundUnprocessed = res.first;
  lastUpperBound = res.second;
  hard_upper_bound = upper_bound;
  return (lastUpperBound == ID_Unsat);
}

template <typename SMALL, typename LARGE>
bool Optimization<SMALL, LARGE>::handleNewSolution(const std::vector<Lit>& sol, bool addConstr, bool softened) {
  LARGE new_val = -origObj->getRhs();
  for (Var v : origObj->vars) new_val += origObj->coefs[v] * (int)(sol[v] > 0);
#ifndef NDEBUG
  bool missing_abstractions = false;
  for (auto& abs : lazies)
    if (!abs.isInSolver() || (!abs.isFinal() && !reformObj->getCoef(abs.watchedVar))) {
      missing_abstractions = true;
      break;
    }
  LARGE new_val_ref = -reformObj->getRhs();
  for (Var v : reformObj->vars) new_val_ref += reformObj->coefs[v] * (int)(sol[v] > 0);
  assert(!wceCoreBuf.empty() || missing_abstractions || new_val >= new_val_ref);
#endif
  if (new_val < upper_bound) {
    bestSol = sol;
    upper_bound = new_val;
    if (addConstr && addUpperBound(softened)) return true;
  }
  return false;
}

template <typename SMALL, typename LARGE>
void Optimization<SMALL, LARGE>::optimize() {
  long long most_lsu_confls = 0, most_cg_confls = 0;
  double hybrid_factor = options.hybridLimitFactor.get();
  Limits limits = Limits::none();
  SolveState reply = SolveState::SAT;
  SMALL coeflim = options.cgStrat ? reformObj->getLargestCoef() : 0;
  bool nextstrat = true;
  bool cgPhase = !options.optMode.is("linear");
  long long start_confls = 0;

  while (lowerBound() < upper_bound) {
    if (options.time_limit.get() != -1.0 && stats.getTime() > options.time_limit.get()) {
      std::cout << "time limit reached\n";
      return;
    }

    if (reply != SolveState::INPROCESSED && reply != SolveState::RESTARTED) {
      printObjBounds();

      start_confls = stats.NCONFL;

      if (options.optMode.is("coreboosted") && stats.getRunTime() >= options.cgBoosted.get()) {
        cgPhase = false;
      } else if (options.optMode.is("hybrid")) {
        // old time-based hybrid strategy
        // cgPhase = cg_time < options.cgHybrid.get() * (lsu_time + cg_time);

        // new limit based hybrid
        // set conflict limit to `hybrid_factor * most_confls` and switch phase
        // if not successful within this limit
        if (options.verbosity.get() > 1) std::cout << "c " << (cgPhase ? "cg" : "lsu") << " phase" << std::endl;
        long long most_confls = cgPhase ? most_cg_confls : most_lsu_confls;
        limits = (most_confls > 0)
                     ? Limits::conflicts(
                           std::max(static_cast<unsigned long long>(options.hybridMinLimit.get()),
                                    static_cast<unsigned long long>(static_cast<double>(most_confls) * hybrid_factor)))
                     : Limits::none();
      }
    }

    assumps.clear();
    if (cgPhase && assumps.isEmpty()) {
      // use core-guided step by setting assumptions
      coeflim = populateAssumptions(coeflim, nextstrat);
      if (options.cgAm1s) {
        if (updatePropGraph()) atMostOnes();
      }
      // update
      if (ubStrict && equalitySwitch) assumps.add(equalitySwitch);
      nextstrat = false;
    } else if (!cgPhase) {
      assumps.clear();
      if (hard_upper_bound > upper_bound) addUpperBound(!options.cgUbConstr.is("<"));
      if (equalitySwitch) assumps.add(equalitySwitch);
    }

    // std::cout << "c nAssumptions for solve: " << assumps.size() << " @ " << stats.getTime() << " s\n";

    SolveAnswer out = aux::timeCall<SolveAnswer>(
        [&] {
          if (reply != SolveState::INPROCESSED && reply != SolveState::RESTARTED) {
            solver.setLimits(limits);
          }
          solver.setAssumptions(assumps);
          return solver.solve();
        },
        cgPhase ? stats.SOLVETIMECG : stats.SOLVETIME);
    reply = out.state;

    if (cgPhase) {
      if (reply != SolveState::INPROCESSED && reply != SolveState::RESTARTED && reply != SolveState::CONFLIMIT) {
        most_cg_confls = std::max(most_cg_confls, stats.NCONFL - start_confls);
      }
    } else {
      if (reply != SolveState::INPROCESSED && reply != SolveState::RESTARTED && reply != SolveState::CONFLIMIT) {
        most_lsu_confls = std::max(most_lsu_confls, stats.NCONFL - start_confls);
      }
    }

    switch (reply) {
      case SolveState::UNSAT:
        break;
      case SolveState::SAT:
        assumps.clear();
        assert(solver.foundSolution());
        ++stats.NSOLS;
        if (handleNewSolution(out.solution, !cgPhase || !options.cgUbConstr.is("no"), !options.cgUbConstr.is("<")))
          break;
        if (wceCoreBuf.empty())
          nextstrat = true;
        else if (finalizeWcePhase())
          break;
        if (harden()) break;
        // hardening might propagate things, so the lazies need to be updated
        if (checkLazyVariables()) break;
        continue;
      case SolveState::INCONSISTENT:
        if (equalitySwitch && solver.isUnit(-equalitySwitch)) break;
        assert(!options.optMode.is("linear"));
        if (!cgPhase) {
          assert(options.optMode.is("hybrid"));
          assert(equalitySwitch);
          break;
        }
        ++stats.NCORES;
        if (handleInconsistency(out.cores)) break;
        if (options.cgLbConstr && addLowerBound()) break;
        if (lowerBound() <= upper_bound) {
          if (harden()) break;
          // hardening might propagate things, so the lazies need to be updated
          if (checkLazyVariables()) break;
        }
        continue;
      case SolveState::RESTARTED:
      case SolveState::INPROCESSED:
        continue;
      case SolveState::CONFLIMIT:
        if (cgPhase && wceCoreBuf.empty()) {
          cgPhase = false;
          if (ubStrict) {
            ubStrict = false;
            if (options.verbosity.get() > 0) std::cout << "c switching to CG(<=)" << std::endl;
          } else {
            cgPhase = false;
            ubStrict = true;
            if (options.verbosity.get() > 0) std::cout << "c switching to LSU" << std::endl;
          }
        } else if (cgPhase) {
          if (finalizeWcePhase()) break;
        } else {
          cgPhase = true;
          // finalize WCE before starting new WCE (add reformulations for found cores)
          if (finalizeWcePhase()) break;
          if (options.verbosity.get() > 0) std::cout << "c switching to CG(<)" << std::endl;
          hybrid_factor *= options.hybridBloom.get();
        }
        assumps.clear();
        continue;
      case SolveState::PROPLIMIT:
        assert(false);  // never setting a propagation limit
    }

    break;
  }
  if (lowerBound() >= upper_bound) {
    printObjBounds();
    logProof();
  }
  quit::exit_UNSAT(solver, upper_bound, bestSol);
  return;
}

template <typename SMALL, typename LARGE>
bool Optimization<SMALL, LARGE>::updatePropGraph() {
  // determine new literals to propagate
  std::vector<Lit> new_lits{};
  for (unsigned i = 0; i < assumps.size(); i++) {
    Lit l = assumps.getKeys()[i];
    if (!pgraph.has_lit(-l)) new_lits.push_back(-l);
  }
  pgraph.add_lits(new_lits);
  // propagate new literals and extend graph
  bool new_edges = false;
  for (Lit l1 : new_lits) {
    if (solver.isUnit(l1) || solver.isUnit(-l1)) {
      assumps.remove(-l1);
      continue;
    }
    IntSet a(1, {-l1});
    PropAnswer out = aux::timeCall<PropAnswer>(
        [&] {
          solver.setAssumptions(a);
          return solver.propagateAssumps();
        },
        stats.PROPTIMEAM1);
    if (out.conflict) {
      // learn unit
#ifndef NDEBUG
      auto rs =
#endif
          solver.addUnitConstraint(l1, Origin::COREGUIDED);
      assert(rs.second != ID_Unsat);
      assumps.remove(-l1);
      continue;
    }
    for (Lit l2 : out.propagated)
      if (pgraph.has_lit(l2)) {
        pgraph.add_edge(l1, l2);
        new_edges = true;
      }
  }
  reformObj->removeUnitsAndZeroes(solver, false);
  return new_edges;
}

template <typename SMALL, typename LARGE>
bool Optimization<SMALL, LARGE>::atMostOnes() {
  // Find maximal cliques in propagation graph and reformulate
  std::vector<Lit> assumpsCopy = assumps.getKeys();
#ifndef NDEBUG
  LARGE lower_bound = lowerBound();
#endif
  for (unsigned i = 0; i < assumpsCopy.size(); i++) {
    Lit l = -assumpsCopy[i];
    if (!assumps.has(-l)) continue;
    if (solver.isUnit(l) || solver.isUnit(-l)) {
      assumps.remove(-l);
      continue;
    }
    auto order = [this](Lit& l1, Lit& l2) { return pgraph.degree(l1) > pgraph.degree(l2); };
    auto isObj = [this](Lit& l1) { return reformObj->getCoef(l1) != 0; };
    std::vector<Lit> clique = pgraph.greedy_max_clique(l, order, isObj);
    if (clique.size() <= 1) continue;
    // Compute minimum clique weight
    SMALL minWeight = reformObj->getCoef(clique[0]);
    for (Lit l : clique) minWeight = std::min(minWeight, reformObj->getCoef(l));
    assert(minWeight > 0);
    // Reformulate found clique
    Ce32 reform = solver.cePools.take32();
    reform->addRhs(1);
#ifndef NDEBUG
    lower_bound += minWeight * (clique.size() - 1);
#endif
    for (Lit l : clique) reform->addLhs(1, -l);
    Var newV = solver.getNbVars() + 1;
    solver.setNbVars(newV);
    reform->addLhs(1, newV);
    if (solver.addConstraint(reform, Origin::COREGUIDED).second == ID_Unsat) return true;
    // reformulate objective
    reformObj->addUp(reform, minWeight);
    assert(!reformObj->hasNoZeroes());
    // update assumptions
    assumps.add(-newV);
    reformObj->removeZeroesWithCallback([&](Var v) {
      assumps.remove(v);
      assumps.remove(-v);
    });
    // Update lazy reformulations if needed
    checkLazyVariables();
    assert(lower_bound == lowerBound());
    stats.NAM1S++;
    stats.AM1LB += static_cast<long long>(minWeight) * (clique.size() - 1);
  }
  return false;
}

template <typename SMALL, typename LARGE>
void Optimization<SMALL, LARGE>::logProof() {
  if (!solver.logger) return;
  assert(lastUpperBound != ID_Undef);
  assert(lastUpperBound != ID_Unsat);
  assert(lastLowerBound != ID_Undef);
  assert(lastLowerBound != ID_Unsat);
  CePtr<ConstrExp<SMALL, LARGE>> coreAggregate = solver.cePools.take<SMALL, LARGE>();
  CePtr<ConstrExp<SMALL, LARGE>> aux = solver.cePools.take<SMALL, LARGE>();
  origObj->copyTo(aux);
  aux->invert();
  aux->addRhs(1 - upper_bound);
  aux->proof.reset(lastUpperBoundUnprocessed);
  coreAggregate->addUp(aux);
  aux->reset();
  origObj->copyTo(aux);
  aux->addRhs(lowerBound());
  aux->proof.reset(lastLowerBoundUnprocessed);
  coreAggregate->addUp(aux);
  assert(coreAggregate->hasNegativeSlack(solver));
  assert(solver.decisionLevel() == 0);
  coreAggregate->logInconsistency(solver, stats);
}

INSTANTIATE_CLASS(Optimization)

}  // namespace run
}  // namespace rs
