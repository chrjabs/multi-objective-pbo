#include "coreguided.hpp"
#include <unordered_set>
#include "Logger.hpp"
#include "auxiliary.hpp"

namespace rs {
namespace cg {

template <typename CF, typename DG>
CePtr<ConstrExp<CF, DG>> expand(const ConstrSimple<CF, DG>& simp) {
  CePtr<ConstrExp<CF, DG>> exp = run::solver.cePools.take<CF, DG>();
  exp->addRhs(simp.rhs);
  for (const Term<CF>& t : simp.terms) exp->addLhs(t.c, t.l);
  exp->orig = simp.orig;
  if (exp->plogger) {
    exp->proof.copyFrom(simp.proofLine);
  }
  return exp;
}

template <typename CF, typename DG>
SumEnc<CF, DG>::SumEnc(const CS& core, std::pair<Var, Var> auxRange) {
  assert(core.isCardinality());
  core.copyTo(reform);
  firstAuxIdx = reform.terms.size();
  for (Var v = auxRange.first; v <= auxRange.second; v++) reform.terms.emplace_back(-1, v);
}

template <typename CF, typename DG>
SumEnc<CF, DG>::SumEnc(const CS& core, std::vector<Term<CF>> auxs) {
  core.copyTo(reform);
  firstAuxIdx = reform.terms.size();
  for (Term<CF>& term : auxs) reform.terms.emplace_back(-term.c, term.l);
}

template <typename CF, typename DG>
LazyEnc<CF, DG>::LazyEnc(const CS& core, DG maxBound, CF mult)
    : subsumComp(nullptr),
      coveredVal(core.rhs),
      maxBound(maxBound),
      mult(mult),
      watchedVar(0),
      atLeastID(ID_Undef),
      atMostID(ID_Undef),
      flags(SYM) {
  core.copyTo(atLeast);
  atLeast.proofLine = "";
  atLeast.toNormalFormLit();
  assert(atLeast.rhs == core.rhs);
  atMost.rhs = -atLeast.rhs;
  atMost.terms.reserve(atLeast.terms.size());
  for (auto& t : atLeast.terms) atMost.terms.emplace_back(-t.c, t.l);
  if (!core.isCardinality()) {
    // initialize subsum computation
    flags |= PBREF;
    std::vector<CF> cfs;
    cfs.reserve(core.terms.size());
    for (Term<CF> t : core.terms) {
      assert(t.c > 0);
      cfs.push_back(t.c);
    }
    subsumComp = std::make_unique<aux::subsum::IncComp<CF, DG>>(cfs);
  }
  // note: after this constructor addVar needs to first be called
}

template <typename CF, typename DG>
bool SumEnc<CF, DG>::addFull() {
  if (run::solver.addConstraint(reform, Origin::COREGUIDED).second == ID_Unsat) {
    return true;
  }
  CS second;
  second.rhs = -reform.rhs;
  second.terms.reserve(reform.terms.size());
  for (auto& t : reform.terms) second.terms.emplace_back(-t.c, t.l);
  if (run::solver.addConstraint(second, Origin::COREGUIDED).second == ID_Unsat) {
    return true;
  }
  for (size_t i = firstAuxIdx; i < reform.terms.size() - 1; i++) {  // add symmetry breaking constraints
    if (run::solver
            .addConstraint(ConstrSimple32({{1, reform.terms[i].l}, {1, -reform.terms[i + 1].l}}, 1), Origin::COREGUIDED)
            .second == ID_Unsat) {
      return true;
    }
  }
  return false;
}

template <typename CF, typename DG>
CF LazyEnc<CF, DG>::addVar(Var v, bool reified) {
  assert(!(flags & FINAL));
  flags &= ~(ATLEAST | ATMOST | SYM);
  // no symmetry breaking for first variables, so mark it as added
  if (!watchedVar) flags |= SYM;

  CF w = 1;
  if (flags & PBREF) w = static_cast<CF>(subsumComp->nextLarger(coveredVal) - coveredVal);
  assert(w > 0);

  // construct next reformulation constraints
  if (reified) {
    if (watchedVar) {
      Term<CF>& last = atLeast.terms.back();
      last = {last.c - w, v};
      atMost.rhs = -coveredVal;
      atMost.terms.back() = {static_cast<CF>(maxBound - coveredVal), v};
    } else {
      atLeast.terms.emplace_back(-w, v);
      atMost.terms.emplace_back(static_cast<CF>(maxBound - coveredVal), v);
    }
    // justification for redundance in proof
    std::stringstream subst1{};
    subst1 << "x" << v << " -> 1 ";
    atMost.proofLine = subst1.str();
    // TODO: consider providing explicity proof for #1 subgoal
    std::stringstream subst2{};
    subst2 << "x" << v << " -> 0 ";
    // subst2 << "x" << v << " -> 0 ; begin" << std::endl
    //        << "proofgoal #1 ; begin " << std::endl
    //        << "p -1 " << coreID << std::endl
    //        << "qed -1" << std::endl
    //        << "qed";
    atLeast.proofLine = subst2.str();
  } else {
    atLeast.terms.emplace_back(-w, v);
    if (watchedVar) {
      assert(atLeast.terms.end()[-2].l == watchedVar);
      CF last_w = -atLeast.terms.end()[-2].c;
      atMost.terms.back().c = last_w;
    }
    atMost.terms.emplace_back(static_cast<CF>(maxBound - coveredVal), v);
    // justification for redundance in proof
    std::stringstream subst1{atMost.proofLine};
    subst1 << " x" << v << " -> 1 ";
    atMost.proofLine = subst1.str();
    // TODO: consider providing explicity proof for #1 subgoal
    std::stringstream subst2{atLeast.proofLine};
    subst2 << " x" << v << " -> 0 ";
    atLeast.proofLine = subst2.str();
  }

  watchedVar = v;
  coveredVal += w;
  return w * mult;
}

template <typename CF, typename DG>
bool LazyEnc<CF, DG>::addConstraints(bool reified, bool addAtLeast, Var prevvar) {
  // Add all constraints at the same time by first adding them to the proof log and then to the solver, to avoid issues
  // from simplifications happening when the constraints are added to the solver
  if (run::solver.logger) {
    CeSuper least = nullptr;
    CeSuper most = nullptr;
    CeSuper sym = nullptr;
#ifndef NDEBUG
    run::solver.logger->logComment("LazyEnc CG constraints");
#endif
    if (addAtLeast) {
      assert(!(flags & ATLEAST));
      assert(atLeast.terms.back().l == watchedVar);
      run::solver.dropExternal(atLeastID, !reified, false);
      least = atLeast.toExpanded(run::solver.cePools);
      least->getProof().rule = ProofRule::REDUNDANT;
      least->logCG();
      flags |= ATLEAST;
    }
    assert(!(flags & ATMOST));
    assert(atMost.terms.back().l == watchedVar);
    run::solver.dropExternal(atMostID, !reified, false);
    most = atMost.toExpanded(run::solver.cePools);
    most->getProof().rule = ProofRule::REDUNDANT;
    most->logCG();
    flags |= ATMOST;
    if (prevvar) {
      assert(!(flags & SYM));
      assert(prevvar < watchedVar);
      // y-- + ~y >= 1 (equivalent to y-- >= y)
      flags |= SYM;
      ConstrSimple32 constr{{{1, prevvar}, {1, -watchedVar}}, 1};
      sym = constr.toExpanded(run::solver.cePools);
      sym->getProof().addToWitness(watchedVar, 0);
      sym->getProof().rule = ProofRule::REDUNDANT;
      sym->logCG();
    }
    if (least) {
      atLeastID = run::solver.addConstraint(least, Origin::COREGUIDED).second;
      if (atLeastID == ID_Unsat) return true;
    };
    atMostID = run::solver.addConstraint(most, Origin::COREGUIDED).second;
    if (atMostID == ID_Unsat) return true;
    if (sym && run::solver.addConstraint(sym, Origin::COREGUIDED).second == ID_Unsat) return true;
  } else {
    if (addAtLeast && addAtLeastConstraint(reified) == ID_Unsat) return true;
    if (addAtMostConstraint(reified) == ID_Unsat) return true;
    if (prevvar && addSymBreakingConstraint(prevvar) == ID_Unsat) return true;
  }
  return false;
}

template <typename CF, typename DG>
ID LazyEnc<CF, DG>::addAtLeastConstraint(bool reified) {
  assert(!(flags & ATLEAST));
  assert(atLeast.terms.back().l == watchedVar);
  run::solver.dropExternal(atLeastID, !reified, false);
  atLeastID = run::solver.addConstraint(atLeast, Origin::COREGUIDED).second;
  flags |= ATLEAST;
  return atLeastID;
}

template <typename CF, typename DG>
ID LazyEnc<CF, DG>::addAtMostConstraint(bool reified) {
  assert(!(flags & ATMOST));
  assert(atMost.terms.back().l == watchedVar);
  run::solver.dropExternal(atMostID, !reified, false);
  atMostID = run::solver.addConstraint(atMost, Origin::COREGUIDED).second;
  flags |= ATMOST;
  return atMostID;
}

template <typename CF, typename DG>
ID LazyEnc<CF, DG>::addSymBreakingConstraint(Var prevvar) {
  assert(!(flags & SYM));
  assert(prevvar < watchedVar);
  // y-- + ~y >= 1 (equivalent to y-- >= y)
  flags |= SYM;
  return run::solver.addConstraint(ConstrSimple32({{1, prevvar}, {1, -watchedVar}}, 1), Origin::COREGUIDED).second;
}

template <typename CF, typename DG>
bool LazyEnc<CF, DG>::exhaust(CE reformObj, bool reified, Limits lims, const IntSet baseAssumps,
                              const std::function<bool(const std::vector<Lit>&)>& solCb,
                              const std::function<void()>& boundCb) {
  CF w = reformObj->getCoef(watchedVar);
  assert(isInSolver());

#ifndef NDEBUG
  DG lb = -reformObj->getDegree();
#endif

  IntSet assumps;
  bool unsat = false;
  bool exhausted = false;
  DG bound_increase = 0;
  while (!exhausted) {
    assumps = baseAssumps;
    assumps.add(-watchedVar);
    SolveAnswer out = aux::timeCall<SolveAnswer>(
        [&] {
          run::solver.setLimits(lims);
          run::solver.setAssumptions(assumps);
          return run::solver.solve();
        },
        stats.SOLVETIMEEXHAUST);
    switch (out.state) {
      case SolveState::SAT:
        unsat = solCb(out.solution);
        [[fallthrough]];
      case SolveState::PROPLIMIT:
      case SolveState::CONFLIMIT:
        exhausted = true;
        break;
      case SolveState::UNSAT:
        unsat = exhausted = true;
        break;
      case SolveState::RESTARTED:
      case SolveState::INPROCESSED:
        continue;
      case SolveState::INCONSISTENT:
        stats.NEXHAUSTED++;
        bound_increase += w;
#ifndef NDEBUG
        lb += w;
#endif

        for (Lit l : baseAssumps.getKeys()) out.cores[0]->harden(l);
        if (out.cores[0]->hasNegativeSlack(run::solver)) {
          if (run::solver.logger) out.cores[0]->logInconsistency(run::solver, stats);
          unsat = exhausted = true;
        }

        // Add extracted core as constraint
        // Typically this is a unit, we add the core for proof production
        if (!run::solver.isUnit(watchedVar) &&
            run::solver.addConstraint(out.cores[0], Origin::COREGUIDED).second == ID_Unsat)
          exhausted = unsat = true;

        if (isFinal()) {
          exhausted = true;
          break;
        }

        // Extend lazy encoding
        int newN = run::solver.getNbVars() + 1;
        run::solver.setNbVars(newN);
        Var oldvar = watchedVar;
        w = addVar(newN, reified);
        // reformulate the objective
        reformObj->addLhs(w, newN);
        // add necessary lazy constraints
        if (addConstraints(reified, options.cgRefLb.get(), oldvar)) {
          unsat = exhausted = true;
        }
        boundCb();
    }
  }
  run::solver.setLimits(Limits::none());
  reformObj->removeUnitsAndZeroes(run::solver, false);
  assert(-reformObj->getDegree() >= lb);  // might have derived more units
  if (options.verbosity.get() > 1)
    std::cout << "c core exhaustion raised lower bound by " << bound_increase << std::endl;
  return unsat;
}

INSTANTIATE_STRUCT(LazyEnc)

INSTANTIATE_STRUCT(SumEnc)

bool trimCore(CeSuper core) {
  bool decCore = (bool)options.cgDecisionCore;
  options.cgDecisionCore.set(false);
  size_t firstSize = core->vars.size();

  bool unsat = false;
  bool done = false;
  while (!done) {
    size_t oldSize = core->vars.size();
    // shuffle assumptions and solve, hoping to get a smaller core
    IntSet assumps = {};
    for (Var v : core->vars) {
      assumps.add(-core->getLit(v));
    }
    assumps.shuffle();
    run::solver.setLimits(Limits::none());
    run::solver.setAssumptions(assumps);
    SolveAnswer out = run::solveFully([&] { return run::solver.solve(); }, stats.SOLVETIME);
    switch (out.state) {
      case SolveState::INCONSISTENT:
        assert(out.cores.size() == 1);
        core = out.cores[0];
        out.cores.clear();
        assert(!core->vars.empty());
        if (core->vars.size() >= oldSize) done = true;
        continue;
      case SolveState::UNSAT:
        unsat = true;
        break;
      case SolveState::SAT:
      case SolveState::PROPLIMIT:
      case SolveState::CONFLIMIT:
      case SolveState::RESTARTED:
      case SolveState::INPROCESSED:
        assert(false);
    }
  }

  stats.NTRIMMED += firstSize - core->vars.size();
  options.cgDecisionCore.set(decCore);
  return unsat;
}

bool internalMin(CeSuper core) {
  IntSet assumps = {};
#ifndef NDEBUG
  size_t origSize = core->vars.size();
#endif
  for (size_t i = 0; i < core->vars.size(); i++) {
    assumps.clear();
    for (size_t j = 0; j < core->vars.size(); j++) {
      if (j == i) continue;
      assumps.add(-core->getLit(core->vars[j]));
    }
    while (true) {
      SolveAnswer out = aux::timeCall<SolveAnswer>(
          [&] {
            run::solver.setLimits(Limits::none());
            run::solver.setAssumptions(assumps);
            return run::solver.solve();
          },
          stats.SOLVETIMEMIN);
      switch (out.state) {
        case SolveState::INCONSISTENT:
          assert(out.cores.size() == 1);
          core = out.cores[0];
          out.cores.clear();
          assert(core->vars.size() < origSize);
          return internalMin(core);
        case SolveState::UNSAT:
          return true;
        case SolveState::SAT:
          break;
        case SolveState::PROPLIMIT:
        case SolveState::CONFLIMIT:
          assert(false);
        case SolveState::RESTARTED:
        case SolveState::INPROCESSED:
          continue;
      }
      break;
    }
  }
  return false;
}

bool minimizeCore(CeSuper core) {
  bool decCore = (bool)options.cgDecisionCore;
  options.cgDecisionCore.set(false);
  size_t oldSize = core->vars.size();
  bool unsat = internalMin(core);
  stats.NMINIMIZED += oldSize - core->vars.size();
  options.cgDecisionCore.set(decCore);
  return unsat;
}

template <typename CF, typename DG>
CF populateAssumptions(IntSet& assumps, CePtr<ConstrExp<CF, DG>> reformObj, CF coeflim, bool refine,
                       const std::vector<Lit>& bestSol) {
  reformObj->removeUnitsAndZeroes(run::solver, false);

  if (refine && coeflim > 0) {
    // find a new coeflim
    reformObj->sortInDecreasingCoefOrder();
    int numLargerCoefs = 0;
    int numLargerUniqueCoefs = 0;
    CF prevCoef = -1;
    bool changed = false;
    for (Var v : reformObj->vars) {
      CF coef = aux::abs(reformObj->coefs[v]);
      ++numLargerCoefs;
      numLargerUniqueCoefs += prevCoef != coef;
      prevCoef = coef;
      if (coeflim > coef && (float)numLargerCoefs / (float)numLargerUniqueCoefs > 1.25) {
        changed = true;
        coeflim = coef;
        break;
      }
    }
    if (!changed) {
      coeflim = 0;
    }
    std::cout << "c stratification level: " << coeflim << std::endl;
  }

  std::vector<Term<double>> litcoefs;  // using double will lead to faster sort than bigint
  litcoefs.reserve(reformObj->vars.size());
  for (Var v : reformObj->vars) {
    assert(reformObj->getLit(v) != 0);
    CF cf = aux::abs(reformObj->coefs[v]);
    if (cf >= coeflim) litcoefs.emplace_back(static_cast<double>(cf), v);
  }
  std::sort(litcoefs.begin(), litcoefs.end(), [](const Term<double>& t1, const Term<double>& t2) {
    return t1.c > t2.c || (t1.l < t2.l && t1.c == t2.c);
  });
  for (const Term<double>& t : litcoefs) assumps.add(-reformObj->getLit(t.l));

  if (options.cgSkipStratLevels && !bestSol.empty() && coeflim > 0) {
    // check if current best solution already satisfies assumptions
    for (Lit l : bestSol)
      if (assumps.has(-l)) return coeflim;
    // assumptions satisfied, refine
    return populateAssumptions(assumps, reformObj, coeflim, true, bestSol);
  }

  return coeflim;
}

template <typename CF, typename DG>
bool addReform(const ConstrSimple<CF, DG>& reform, CePtr<ConstrExp<CF, DG>> reformObj, Lazies<CF, DG>& lazies, CF mult,
               DG upperBound, UbCb<DG> ubCb, ExhaustCb<CF, DG> exCb) {
  assert(reform.isNormalFormLit());
  if (!reform.isCardinality()) stats.NPBREFORMS++;

#ifndef NDEBUG
  DG lower_bound = -reformObj->getDegree() + reform.rhs * mult;
#endif

  if (options.cgCoreUpper) {
    upperBound = std::min<DG>(upperBound, ubCb() / mult);
    // NOTE: integer division rounds down
  }
  assert(upperBound >= 0);

  CePtr<ConstrExp<CF, DG>> ceReform = expand(reform);
  ceReform->removeUnitsAndZeroes(run::solver, false);
  if (ceReform->isTautology()) {
    // derived units since extracting this core that subsume the core
    return false;
  }
  if (ceReform->hasNegativeSlack(run::solver)) {
    // derived contradiction since extracting this core
    assert(run::solver.decisionLevel() == 0);
    if (run::solver.logger) ceReform->logInconsistency(run::solver, stats);
    return true;
  }

  // reformulate objective: step 1, decrease core lit coefficients
  ceReform->invert();
  reformObj->addUp(ceReform, mult);
  assert(!reformObj->hasNoZeroes());
  assert(-reformObj->getDegree() == lower_bound);

  if (options.cgEncoding.is("sum")) {
    if (reform.isCardinality()) {
      int oldN = run::solver.getNbVars();
      int newN = oldN + static_cast<int>(upperBound - reform.rhs);
      run::solver.setNbVars(newN);
      for (Var v = oldN + 1; v <= newN; ++v) reformObj->addLhs(mult, v);
      SumEnc enc(reform, std::make_pair(oldN + 1, newN));
      enc.addFull();
    } else {
      // calculate all possible values (these are exponential, which is why the
      // sum encoding is not advisible with pb reformulations)
      std::vector<CF> cfs;
      cfs.resize(reform.terms.size());
      for (Term<CF> t : reform.terms) cfs.push_back(t.c);
      std::vector<DG> subsums = aux::subsum::all<CF, DG>(cfs);
      std::sort(subsums.begin(), subsums.end());
      subsums.erase(std::unique(subsums.begin(), subsums.end()), subsums.end());
      subsums.shrink_to_fit();
      // compute weights of aux vars
      DG lastVal = reform.rhs;
      Var lastVar = run::solver.getNbVars();
      std::vector<Term<CF>> auxs;
      auxs.resize(subsums.size());
      for (DG val : subsums) {
        if (val <= lastVal) {
          assert(val <= reform.rhs);
          continue;
        }
        CF w = static_cast<CF>(val - lastVal);
        auxs.emplace_back(w, ++lastVar);
        reformObj->addLhs(w, lastVar);
      }
      run::solver.setNbVars(lastVar);
      SumEnc enc(reform, auxs);
      enc.addFull();
    }
  } else {
    bool reified = options.cgEncoding.is("reified");
    int newV = run::solver.getNbVars() + 1;
    run::solver.setNbVars(newV);
    // reformulate objective: step 2, lazily add one new literal
    lazies.emplace_back(reform, upperBound, mult);
    LazyEnc<CF, DG>& enc = lazies.back();
    CF w = lazies.back().addVar(newV, reified);
    reformObj->addLhs(w, newV);
    // no need to add symmetry breaking the first time
    if (enc.addConstraints(reified, options.cgRefLb.get(), 0)) return true;
    if (exCb(lazies.back())) return true;
    if (enc.isFinal()) {
      lazies.pop_back();
    }
  }
  return false;
}

template <typename CF, typename DG>
bool addCardReform(const ConstrSimple32& reform, CePtr<ConstrExp<CF, DG>> reformObj, Lazies<CF, DG>& lazies, CF mult,
                   UbCb<DG> ubCb, ExhaustCb<CF, DG> exCb) {
  assert(reform.isCardinality());
  assert(reform.isNormalFormLit());
  ConstrSimple<CF, DG> pbReform{};
  reform.copyTo(pbReform);
  return addReform(pbReform, reformObj, lazies, mult, static_cast<DG>(reform.terms.size()), ubCb, exCb);
}

// reformulate cores in `wceCoreBuf` with card reformulations
template <typename CF, typename DG>
bool reformulateCard(WCores<CF>& cores, CePtr<ConstrExp<CF, DG>> reformObj, Lazies<CF, DG>& lazies, UbCb<DG> ubCb,
                     ExhaustCb<CF, DG> exCb) {
  // reformulate the cores in exactly the order they were extracted in
  for (const auto& wcore : cores) {
    if (addCardReform(wcore.first, reformObj, lazies, wcore.second, ubCb, exCb)) return true;
  }
  cores.clear();
  return false;
}

// reformulate cores in `wceCoreBuf` with pb reformulations
template <typename CF, typename DG>
bool reformulatePb(WCores<CF>& cores, CePtr<ConstrExp<CF, DG>> reformObj, Lazies<CF, DG>& lazies, UbCb<DG> ubCb,
                   ExhaustCb<CF, DG> exCb) {
  ConstrSimple<CF, DG> pbReform;
  CF mult;
  DG ub;
  for (const auto& wcore : cores) {
    if (wcore.first.terms.size() <= options.cgMaxPbRef.get()) {
      pbReform = buildPbReform(wcore.first, reformObj, mult, ub);
      if (addReform(pbReform, reformObj, lazies, mult, ub, ubCb, exCb)) return true;
    } else {
      if (addCardReform(wcore.first, reformObj, lazies, wcore.second, ubCb, exCb)) return true;
    }
  }
  cores.clear();
  assert(cores.empty());
  return false;
}

// automatically choose how to reformulate cores in `wceCoreBuf`
template <typename CF, typename DG>
bool reformulateAuto(WCores<CF>& cores, CePtr<ConstrExp<CF, DG>> reformObj, Lazies<CF, DG>& lazies, UbCb<DG> ubCb,
                     ExhaustCb<CF, DG> exCb) {
  // initialize each core as a group
  std::vector<std::pair<std::unordered_set<Var>, std::vector<size_t>>> groups;
  for (size_t idx = 0; idx < cores.size(); idx++) {
    auto& wcore = cores[idx];
    std::unordered_set<Var> set = {};
    std::vector<size_t> cores = {idx};
    for (Term<int> t : wcore.first.terms) {
      set.insert(toVar(t.l));
    }
    groups.emplace_back(std::make_pair(set, cores));
  }
  // merge groups until fully disjoint
  size_t last_size = groups.size() + 1;
  while (groups.size() > 1 && groups.size() < last_size) {
    last_size = groups.size();
    // pairwise check overlap and merge
    for (size_t i = 0; i < groups.size(); i++) {
      for (size_t j = i + 1; j < groups.size(); j++) {
        bool overlaps = false;
        for (Var v : groups[j].first) {
          if (groups[i].first.count(v)) {
            overlaps = true;
            break;
          }
        }
        if (!overlaps) continue;
        // both groups overlap, merge
        groups[i].first.insert(groups[j].first.begin(), groups[j].first.end());
        groups[i].second.insert(groups[i].second.end(), groups[j].second.begin(), groups[j].second.end());
        aux::swapErase(groups, j);
        j--;
      }
    }
  }

  if (options.verbosity.get() > 1) std::cout << "c found " << groups.size() << " disjoint core groups" << std::endl;

  enum ReformType { CARD, PB };
  while (!groups.empty()) {
    auto& group = groups.back();
    // decide about how to reformulate each group
    // implemented cases:
    // CARD. cardinality reformulate the entire group (in the original order for now)
    // PB.   pb reformulate only the core leading the the maximum lower bound increase
    // The strategy is the following for now: choose the core that would lead to the highest lower bound increase,
    // then check whether any literals with coefficient > .5 * avg coefficient appear in any other cores. If yes, we
    // assume that extending the pb reformulation until that literal is too costly and therefore give up -> choose
    // cardinality reformulation.
    long long pbCoreIdx = -1;
    DG bestLbInc = 0;
    // choose core with highest pb lb increase
    for (size_t i = 0; i < group.second.size(); i++) {
      ConstrSimple32& core = cores[group.second[i]].first;
      if (core.terms.size() > options.cgMaxPbRef.get()) continue;
      CF mult;
      DG ub;
      ConstrSimple<CF, DG> pbRef = buildPbReform(core, reformObj, mult, ub);
      if (pbRef.rhs > bestLbInc) {
        pbCoreIdx = i;
        bestLbInc = pbRef.rhs;
      }
    }
    ReformType type;
    if (pbCoreIdx < 0) {
      // all cores are larger than the maximum core to be treated as a pb reformulations
      type = CARD;
    } else {
      ConstrSimple32& pbCore = cores[group.second[pbCoreIdx]].first;
      // compute threshold for coefficients of literals that appear in other cores
      double thr = 0;
      for (Term<int> t : pbCore.terms) thr += static_cast<double>(aux::abs(reformObj->getCoef(t.l)));
      thr = thr / pbCore.terms.size() * .5;
      type = PB;
      for (size_t i = 0; i < group.second.size(); i++) {
        if (i == static_cast<size_t>(pbCoreIdx)) continue;
        ConstrSimple32& core = cores[group.second[i]].first;
        for (Term<int> t : core.terms) {
          if (static_cast<double>(aux::abs(reformObj->getCoef(t.l))) > thr) {
            type = CARD;
            break;
          }
        }
      }
    }
    if (options.verbosity.get() > 1)
      std::cout << "c reformulating core group of size " << group.second.size() << " as ";
    // build reformulation
    ConstrSimple<CF, DG> pbReform;
    CF mult;
    DG ub;
    switch (type) {
      case CARD:
        // make sure to reformulate cores in the order they were extracted
        sort(group.second.begin(), group.second.end());
        if (options.verbosity.get() > 1) std::cout << "cardinalities" << std::endl;
        for (size_t idx = 0; idx < group.second.size(); idx++) {
          auto& wcore = cores[group.second[idx]];
          if (addCardReform(wcore.first, reformObj, lazies, wcore.second, ubCb, exCb)) return true;
        }
        break;
      case PB:
        if (options.verbosity.get() > 1)
          std::cout << "pb reformulation of core with index " << group.second[pbCoreIdx] << std::endl;
        stats.NDISCARDEDCORES += group.second.size() - 1;
        pbReform = buildPbReform(cores[group.second[pbCoreIdx]].first, reformObj, mult, ub);
        if (addReform(pbReform, reformObj, lazies, mult, ub, ubCb, exCb)) return true;
        break;
    }
    groups.pop_back();
  }
  assert(groups.empty());
  cores.clear();

  return false;
}

template <typename CF, typename DG>
bool finalizeWcePhase(WCores<CF>& cores, CePtr<ConstrExp<CF, DG>> reformObj, Lazies<CF, DG>& lazies, UbCb<DG> ubCb,
                      ExhaustCb<CF, DG> exCb) {
  if (options.verbosity.get() > 0) std::cout << "c finalizing WCE round: reformulating cores" << std::endl;
  if (options.verbosity.get() > 1)
    std::cout << "c lower bound from core extraction " << -reformObj->getDegree() << std::endl;
  // reset the reformulated objective to before the WCE phase to be able to choose how to actually reformulate
  for (auto& wcore : cores) {
    Ce32 exp = expand(wcore.first);
    reformObj->addUp(exp, wcore.second);
  }

  bool ret = false;
  if (options.cgObjCoreCoefs.is("no"))
    ret = reformulateCard(cores, reformObj, lazies, ubCb, exCb);
  else if (options.cgObjCoreCoefs.is("yes"))
    ret = reformulatePb(cores, reformObj, lazies, ubCb, exCb);
  else if (options.cgObjCoreCoefs.is("auto"))
    ret = reformulateAuto(cores, reformObj, lazies, ubCb, exCb);
  else
    assert(false);

  if (options.verbosity.get() > 1)
    std::cout << "c lower bound from core reformulation " << -reformObj->getDegree() << std::endl;
  return ret;
}

// Build a PB reformulation from a cardinality core
template <typename CF, typename DG>
ConstrSimple<CF, DG> buildPbReform(const ConstrSimple32& core, const CePtr<ConstrExp<CF, DG>> reformObj, CF& mult,
                                   DG& upperBound) {
  assert(!core.terms.empty());
  assert(core.isNormalFormLit());
  upperBound = 0;
  ConstrSimple<CF, DG> reform{};
  long long coreDegree = core.rhs;
  for (Term<int> t : core.terms) {
    assert(t.c == 1);
    if (run::solver.isUnit(t.l) || run::solver.isUnit(-t.l)) {
      // if unit for literal has been derived, it was removed from the objective
      if (run::solver.isUnit(t.l)) coreDegree--;
      continue;
    }
    assert(reformObj->getCoef(t.l) > 0);
    reform.terms.emplace_back(reformObj->getCoef(t.l), t.l);
    upperBound += reformObj->getCoef(t.l);
  }
  std::sort(reform.terms.begin(), reform.terms.end(), [](Term<CF> a, Term<CF> b) { return a.c < b.c; });
  DG newDegree = 0;
  for (int idx = 0; idx < coreDegree; idx++) newDegree += reform.terms[idx].c;
  reform.rhs = newDegree;
  mult = reform.divideByGcd();
  upperBound /= mult;
  return reform;
}

template <typename CF, typename DG>
std::pair<Ce32, int> reduceToCardinality(const CeSuper& core, const CePtr<ConstrExp<CF, DG>> reformObj) {
  int bestNbBlocksRemoved = 0;
  CeSuper card = core->clone(run::solver.cePools);
  if (options.cgReduction.is("clause")) {
    card->sortInDecreasingCoefOrder(
        [&](Var v1, Var v2) { return aux::abs(reformObj->coefs[v1]) > aux::abs(reformObj->coefs[v2]); });
    card->simplifyToClause();
  } else if (options.cgReduction.is("minslack")) {
    card->sortInDecreasingCoefOrder(
        [&](Var v1, Var v2) { return aux::abs(reformObj->coefs[v1]) > aux::abs(reformObj->coefs[v2]); });
    card->simplifyToMinLengthCardinality();
  } else {
    assert(options.cgReduction.is("bestbound"));
    CeSuper cloneCoefOrder = card->clone(run::solver.cePools);
    cloneCoefOrder->sortInDecreasingCoefOrder();
    std::reverse(cloneCoefOrder->vars.begin(), cloneCoefOrder->vars.end());  // *IN*creasing coef order
    card->sort([&](Var v1, Var v2) { return aux::abs(reformObj->coefs[v1]) > aux::abs(reformObj->coefs[v2]); });
    CeSuper clone = card->clone(run::solver.cePools);
    assert(clone->vars.size() > 0);
    DG bestLowerBound = 0;
    int bestNbVars = clone->vars.size();

    // find the optimum number of variables to weaken to
    int nbBlocksRemoved = 0;
    while (!clone->isTautology()) {
      int carddegree = cloneCoefOrder->getCardinalityDegreeWithZeroes();
      if (bestLowerBound < aux::abs(reformObj->coefs[clone->vars.back()]) * carddegree) {
        bestLowerBound = aux::abs(reformObj->coefs[clone->vars.back()]) * carddegree;
        bestNbVars = clone->vars.size();
        bestNbBlocksRemoved = nbBlocksRemoved;
      }
      CF currentObjCoef = aux::abs(reformObj->coefs[clone->vars.back()]);
      // weaken all lowest objective coefficient literals
      while (clone->vars.size() > 0 && currentObjCoef == aux::abs(reformObj->coefs[clone->vars.back()])) {
        ++nbBlocksRemoved;
        Var v = clone->vars.back();
        clone->weakenLast();
        cloneCoefOrder->weaken(v);
      }
    }

    // weaken to the optimum number of variables and generate cardinality constraint
    while ((int)card->vars.size() > bestNbVars) {
      card->weakenLast();
    }
    card->saturate();
    // sort in decreasing coef order to minimize number of auxiliary variables, but break ties so that *large*
    // objective coefficient literals are removed first, as the smallest objective coefficients are the ones that
    // will be eliminated from the objective function. Thanks Stephan!
    // TODO: check if this actually makes sense
    card->sortInDecreasingCoefOrder(
        [&](Var v1, Var v2) { return aux::abs(reformObj->coefs[v1]) < aux::abs(reformObj->coefs[v2]); });
    card->simplifyToCardinality(false);
  }

  Ce32 result = run::solver.cePools.take<int, long long>();
  card->copyTo(result);
  return {result, bestNbBlocksRemoved};
}

template <typename CF, typename DG>
bool handleInconsistency(std::vector<CeSuper>& origCores, WCores<CF>& coreBuf, CePtr<ConstrExp<CF, DG>> reformObj,
                         Lit equalitySwitch) {
  std::vector<std::pair<Ce32, int>> cores;
  for (CeSuper& core : origCores) {
    assert(core);
    CeSuper cpy = core->clone(run::solver.cePools);
    core->removeUnitsAndZeroes(run::solver, false);
    if (core->isTautology()) continue;
    if (core->hasNegativeSlack(run::solver)) {
      assert(run::solver.decisionLevel() == 0);
      if (run::solver.logger) core->logInconsistency(run::solver, stats);
      return true;
    }
    if (options.cgTrim && trimCore(core)) return true;
    if (options.cgMinimize && minimizeCore(core)) return true;

    // Filter out equality switch
    if (equalitySwitch) core->harden(equalitySwitch);
    core->removeUnitsAndZeroes(run::solver, false);

    // double check after core minimization
    if (core->isTautology()) continue;
    if (core->hasNegativeSlack(run::solver)) {
      assert(run::solver.decisionLevel() == 0);
      if (run::solver.logger) core->logInconsistency(run::solver, stats);
      return true;
    }

    if (core->isClause())
      stats.NCLAUSALCORES++;
    else if (core->isCardinality())
      stats.NCARDCORES++;
    else
      stats.NPBCORES++;

    // figure out an appropriate core
    cores.push_back(cg::reduceToCardinality(core, reformObj));
  }

#ifndef NDEBUG
  DG prev_lower = -reformObj->getDegree();
#endif
  reformObj->removeUnitsAndZeroes(run::solver, false);
  // clean up cores one more time after _all_ minimization which might have derived more units
  unsigned j = 0;
  for (unsigned i = 0; i < cores.size(); i++) {
    auto& core = cores[i].first;
    core->removeUnitsAndZeroes(run::solver, false);
    if (core->isTautology()) continue;
    if (core->hasNegativeSlack(run::solver)) {
      assert(run::solver.decisionLevel() == 0);
      if (run::solver.logger) core->logInconsistency(run::solver, stats);
      return true;
    }
    if (i != j) cores[j] = std::move(cores[i]);
    j++;
  }
  assert(j <= cores.size());
  cores.resize(j);

  if (cores.size() == 0) {
    // only violated unit assumptions were derived
    ++stats.UNITCORES;
    assert(-reformObj->getDegree() > prev_lower);
    return false;
  }

#ifndef NDEBUG
  prev_lower = -reformObj->getDegree();
#endif

  stats.SINGLECORES += cores.size() == 1;

  DG bestLowerBound = -1;
  Ce32& coreToUse = cores[0].first;
  int bestBlocksRemoved = 0;
  for (int i = 0; i < (int)cores.size(); ++i) {
    Ce32& core = cores[i].first;
    assert(core->hasNoZeroes());
    assert(core->vars.size() > 0);
    CF lowestCoef = aux::abs(reformObj->coefs[core->vars[0]]);
    for (Var v : core->vars) {
      if (aux::abs(reformObj->coefs[v]) < lowestCoef) lowestCoef = aux::abs(reformObj->coefs[v]);
    }
    DG lowerBound = lowestCoef * core->getDegree();
    if (i == 1) {
      stats.NOCOREBEST += lowerBound == bestLowerBound;
      stats.FIRSTCOREBEST += lowerBound < bestLowerBound;
      stats.DECCOREBEST += lowerBound > bestLowerBound;
    }
    if (lowerBound > bestLowerBound) {
      bestLowerBound = lowerBound;
      coreToUse = core;
      bestBlocksRemoved = cores[i].second;
    }
  }

  stats.REMOVEDBLOCKS += bestBlocksRemoved;
  stats.NCORECARDINALITIES += !coreToUse->isClause();
  stats.COREDEGSUM += static_cast<long long>(coreToUse->getDegree());
  stats.CORESLACKSUM += coreToUse->vars.size() - static_cast<long long>(coreToUse->getDegree());

  // compute the minimum core weight
  CF mult = 0;
  for (Var v : coreToUse->vars) {
    assert(reformObj->getLit(v) != 0);
    CF w = aux::abs(reformObj->coefs[v]);
    if (mult == 0) mult = w;
    if (w < mult) {
      mult = w;
    }
  }
  assert(mult > 0);

  if (coreToUse->isCardinality() && (int)coreToUse->vars.size() <= coreToUse->getDegree()) {
    assert((int)coreToUse->vars.size() == coreToUse->getDegree());
    // Immediately add assignment core
    for (Var v : coreToUse->vars)
      if (run::solver.addUnitConstraint(coreToUse->getLit(v), Origin::COREGUIDED).second == ID_Unsat) return true;
    reformObj->removeUnitsAndZeroes(run::solver);
    return false;
  }

  ConstrSimple32 coreSimp{};
  coreToUse->toSimple()->copyTo(coreSimp);
  coreSimp.toNormalFormLit();
  if (options.cgObjCoreCoefs.is("yes") && coreSimp.terms.size() <= options.cgMaxPbRef.get()) {
    DG ub;
    CF mult;
    ConstrSimple<CF, DG> pbRef = buildPbReform(coreSimp, reformObj, mult, ub);
    CePtr<ConstrExp<CF, DG>> pbExp = expand(pbRef);
    pbExp->invert();
    reformObj->addUp(pbExp, mult);
  } else {
    // for wce, pretend to be doing cardinality reformulation
    // in auto mode, we will actually revert this later in `finalizeWcePhase` and decide how to actually reformulate
    coreToUse->invert();
    reformObj->addUp(coreToUse, mult);
    coreToUse->invert();
  }
  assert(!reformObj->hasNoZeroes());
  // store the core to later reformulate
  coreBuf.emplace_back(std::make_pair(coreSimp, mult));

  reformObj->removeUnitsAndZeroes(run::solver, false);

  return false;
}

template <typename CF, typename DG>
bool coreGuidedReformulateWithLazies(CePtr<ConstrExp<CF, DG>> reformObj, const ConstrSimple<CF, DG>& origObj,
                                     std::vector<Lit>& optSol, Lazies<CF, DG>& lazies) {
  WCores<CF> wceCoreBuf{};
  DG norm_offset = 0;
  for (auto& t : origObj.terms)
    if (t.c < 0) norm_offset -= t.c;
  DG upper_bound = maxObjVal(reformObj) + 1;
  CF coeflim = reformObj->getLargestCoef();
  bool refine = true;
  IntSet assumps;

  // Callback implementations
  SolCb handleSol = [&](const std::vector<Lit>& sol) -> bool {
    DG val = computeCost(origObj, sol);
    if (val < upper_bound) {
      upper_bound = val;
      optSol = sol;
    }
    return false;
  };
  BoundCb printBounds = [&] {
    if (options.verbosity.get() >= 2)
      std::cout << "c " << -reformObj->getDegree() << " <= " << upper_bound << std::endl;
  };
  UbCb<DG> ubCb = [&]() -> DG { return upper_bound + norm_offset; };
  ExhaustCb<CF, DG> exhaustCb = [=](LazyEnc<CF, DG>& enc) -> bool {
    if (!options.cgExhaust) return false;
    assert(!options.cgEncoding.is("sum"));
    bool reified = options.cgEncoding.is("reified");
    IntSet assumps{};
    int propLim = options.cgExPropLim.get();
    int confLim = options.cgExConfLim.get();
    return enc.exhaust(reformObj, reified,
                       Limits(propLim < 0 ? Limits::UNLIM : propLim, confLim < 0 ? Limits::UNLIM : confLim), assumps,
                       handleSol, printBounds);
  };

  bool hasSol = false;

  // Main loop
  while (-reformObj->getDegree() < upper_bound) {
    assumps.clear();
    coeflim = populateAssumptions(assumps, reformObj, coeflim, refine, optSol);
    refine = false;

    SolveAnswer out = aux::timeCall<SolveAnswer>(
        [&] {
          run::solver.setAssumptions(assumps);
          return run::solver.solve();
        },
        stats.SOLVETIMECG);

    switch (out.state) {
      case SolveState::UNSAT:
        return true;
      case SolveState::SAT:
        hasSol = true;
        handleSol(out.solution);
        if (wceCoreBuf.empty()) {
          refine = true;
          if (coeflim <= 1) {
            assert(-reformObj->getDegree() == upper_bound);
            break;
          }
        }
        if (finalizeWcePhase(wceCoreBuf, reformObj, lazies, ubCb, exhaustCb)) return true;
        if (harden(reformObj, upper_bound)) return true;
        // hardening might propagate things, so the lazies need to be updated
        if (checkLazyVariables(reformObj, lazies)) return true;
        continue;
      case SolveState::INCONSISTENT:
        ++stats.NCORES;
        if (handleInconsistency(out.cores, wceCoreBuf, reformObj, 0)) return true;
        if (-reformObj->getDegree() <= upper_bound) {
          if (harden(reformObj, upper_bound)) return true;
          // hardening might propagate things, so the lazies need to be updated
          if (checkLazyVariables(reformObj, lazies)) return true;
        }
        printBounds();
        continue;
      case SolveState::RESTARTED:
      case SolveState::INPROCESSED:
        continue;
      case SolveState::CONFLIMIT:
      case SolveState::PROPLIMIT:
        assert(false);
        break;
    }
  }

  if (!hasSol) {
    if (run::solver.logger) {
      // prove final unsat
      run::solver.clearAssumptions();
      SolveAnswer out = run::solveFully([&] { return run::solver.solve(); }, stats.SOLVETIME);
      assert(out.state == SolveState::UNSAT);
    }
    return true;
  }

  //  since the solution is already optimal, all newly added variables can be assigend to 0
  for (Var v = toVar(optSol.back()) + 1; v <= run::solver.getNbVars(); v++) optSol.push_back(-v);

  if (wceCoreBuf.empty()) return false;

  if (options.verbosity.get() > 2)
    std::cout << "c CG building final reformulation after terminating by bounds" << std::endl;

  // terminated by bounds, build final _cardinality_ reformulation and extend solution
  auto refOpt = options.cgObjCoreCoefs.get();
  options.cgObjCoreCoefs.set("no");
  if (finalizeWcePhase(wceCoreBuf, reformObj, lazies, ubCb, exhaustCb)) return true;
  options.cgObjCoreCoefs.set(refOpt);
  //  since the solution is already optimal, all newly added variables can be assigend to 0
  for (Var v = toVar(optSol.back()) + 1; v <= run::solver.getNbVars(); v++) optSol.push_back(-v);

  return false;
}

template <typename CF, typename DG>
bool coreGuidedReformulate(CePtr<ConstrExp<CF, DG>> reformObj, const ConstrSimple<CF, DG>& origObj,
                           std::vector<Lit>& optSol) {
  Lazies<CF, DG> lazies;

  if (coreGuidedReformulateWithLazies(reformObj, origObj, optSol, lazies)) return true;

  // finalize objective reformulation by reformulating remaining lazies
  if (checkLazyVariables(reformObj, lazies)) return true;
  bool reified = options.cgEncoding.is("reified");
  for (auto& lazy : lazies) {
    while (!lazy.isFinal() && !run::solver.isUnit(-lazy.watchedVar)) {
      int newN = run::solver.getNbVars() + 1;
      run::solver.setNbVars(newN);
      Var oldvar = lazy.watchedVar;
      CF w = lazy.addVar(newN, reified);
      // reformulate the objective
      reformObj->addLhs(w, newN);
      // add necessary lazy constraints
      if (lazy.addConstraints(reified, options.cgRefLb.get(), oldvar)) return true;
    }
  }

  return false;
}

template bool coreGuidedReformulate(CePtr<ConstrExp<int, long long>>, const ConstrSimple<int, long long>&,
                                    std::vector<Lit>&);
template bool coreGuidedReformulate(CePtr<ConstrExp<long long, int128>>, const ConstrSimple<long long, int128>&,
                                    std::vector<Lit>&);
template bool coreGuidedReformulate(CePtr<ConstrExp<int128, int128>>, const ConstrSimple<int128, int128>&,
                                    std::vector<Lit>&);
template bool coreGuidedReformulate(CePtr<ConstrExp<int128, int256>>, const ConstrSimple<int128, int256>&,
                                    std::vector<Lit>&);
template bool coreGuidedReformulate(CePtr<ConstrExp<bigint, bigint>>, const ConstrSimple<bigint, bigint>&,
                                    std::vector<Lit>&);
}  // namespace cg
}  // namespace rs
