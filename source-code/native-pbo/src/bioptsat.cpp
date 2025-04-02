#include "bioptsat.hpp"
#include "Logger.hpp"
#include "coreguided.hpp"

namespace rs {
namespace run {

template <typename CF, typename DG>
void BiOptSat<CF, DG>::solve() {
  if (options.moCoreBoosting && MO::coreBoosting()) return;

  std::vector<Lit> solution;
  std::vector<DG> costs;
  while (true) {
    solution.clear();
    costs.clear();
    if (minIncr(solution, costs)) break;
    bool unsat = minDecr(solution, costs);
    MO::addSolution(solution, costs);
    if (unsat) break;
    // improving constraint has already been added by minDecr.
    // when getting to implement enumeration version, need to add the hard constraint here.
  }
}

template <typename CF, typename DG>
bool BiOptSat<CF, DG>::minIncrSolImpr(std::vector<Lit> &solution, std::vector<DG> &costs) {
  IntSet assumps{};
  while (true) {
    assumps.clear();
    Lit ol = 0;
    if (!costs.empty()) {
      if (costs[0] <= -MO::reformObjs[0]->getDegree()) break;
      ol = rObjectives[0].getOutput(costs[0] - 1);
      if (!ol) {
        // Adding the reified objective constraint was detected as UNSAT
        // This should probably never happen, and if it does it should produce a correct UNSAT proof
        solver.clearAssumptions();
        return true;
      }
      assumps.add(-ol);
    }
    SolveAnswer out = solveFully(
        [&] {
          solver.setAssumptions(assumps);
          return solver.solve();
        },
        stats.SOLVETIME);
    switch (out.state) {
      case SolveState::UNSAT:
        solver.clearAssumptions();
        return true;
      case SolveState::INCONSISTENT:
        if (run::solver.logger) incrLb = out.cores[0]->logProofLineWithInfo("minIncrSolImpr core", stats);
        solver.clearAssumptions();
        if (solver.logger && !assumps.isEmpty() && ol > toVar(solution.back())) {
          // when proof logging, make sure that the assignment is complete for all constraints in the proof
          assert(ol == toVar(solution.back()) + 1);
          solution.push_back(ol);
        }
        return false;
      case SolveState::SAT:
        ++stats.NSOLS;
        break;
      default:
        assert(false);
        break;
    }
    solution = out.solution;
    costs = MO::computeCosts(solution);
    if (options.verbosity.get() >= 3) {
      std::cout << "c bioptsat increasing candidate (";
      for (size_t idx = 0; idx < costs.size(); idx++)
        std::cout << costs[idx] << ((idx < costs.size() - 1) ? ", " : ")");
      std::cout << std::endl;
    }
  }
  solver.clearAssumptions();
  return false;
}

template <typename CF, typename DG>
bool BiOptSat<CF, DG>::minIncrOll(std::vector<Lit> &solution, std::vector<DG> &costs) {
  bool cgCoreUpper = (bool)options.cgCoreUpper;
  options.cgCoreUpper.set(false);
  bool cgHardening = (bool)options.cgHardening;
  options.cgHardening.set(false);

  bool unsat = cg::coreGuidedReformulateWithLazies(MO::reformObjs[0], MO::origObjs[0], solution, ollLazies);
  costs = MO::computeCosts(solution);

  options.cgCoreUpper.set(cgCoreUpper);
  options.cgHardening.set(cgHardening);
  return unsat;
}

template <typename CF, typename DG>
bool BiOptSat<CF, DG>::minDecr(std::vector<Lit> &solution, std::vector<DG> &costs) {
  assert(!solution.empty());
  assert(costs.size() == 2);
  // Solution improving search on the second objective while pinning the first.
  // Constraints are added as hard constraints.
  IntSet assumps{};
  if (options.bosVariant.is("sat-unsat")) {
    Lit ol = rObjectives[0].getOutput(costs[0]);
    assert(ol);
    assumps.add(-ol);
    if (solver.logger && ol > toVar(solution.back())) {
      // when proof logging, make sure that the assignment is complete for all constraints in the proof
      assert(ol == toVar(solution.back()) + 1);
      solution.push_back(-ol);
    }
  } else if (options.bosVariant.is("oll")) {
    cg::populateAssumptions(assumps, MO::reformObjs[0], static_cast<CF>(1), false, solution);
  }
  while (true) {
    if (costs[1] <= -MO::reformObjs[1]->getDegree()) {
      if (solver.logger) {
        // Prove UNSAT
        // (i) PD-cut
        auto pd_cut = MO::pdCut(costs, rObjectives);
        assert(pd_cut.terms.size() <= 1);
        MO::logPdCut(pd_cut, solution);
        auto pd_cut_e = pd_cut.toExpanded(solver.cePools);
        pd_cut_e->getProof().rule = ProofRule::RUP;
#ifndef NDEBUG
        run::solver.logger->logComment("BiOptSat preliminary PD cut");
#endif
        auto cut_id = pd_cut_e->logRup();
        if (!pd_cut.terms.empty()) {
          // (ii) Combine with lower bound on increasing objective for contradiction
          solver.logger->last_proofID++;
          solver.logger->proof_out << "p " << cut_id << " " << incrLb << " + " << std::endl;
        }
      }
      solver.clearAssumptions();
      return true;
    }
    if (addDecConstr(costs, solution)) {
      // Adding the unit detected UNSAT
      // If this happens, a correct UNSAT proof is already produced
      solver.clearAssumptions();
      return true;
    };
    SolveAnswer out = solveFully(
        [&] {
          solver.setAssumptions(assumps);
          return solver.solve();
        },
        stats.SOLVETIME);
    switch (out.state) {
      case SolveState::UNSAT:
        solver.clearAssumptions();
        return true;
      case SolveState::INCONSISTENT:
        solver.clearAssumptions();
        return false;
      case SolveState::SAT:
        ++stats.NSOLS;
        break;
      default:
        assert(false);
        break;
    }
    solution = out.solution;
    costs = MO::computeCosts(solution);
    if (options.verbosity.get() >= 3) {
      std::cout << "c bioptsat decreasing candidate (";
      for (size_t idx = 0; idx < costs.size(); idx++)
        std::cout << costs[idx] << ((idx < costs.size() - 1) ? ", " : ")");
      std::cout << std::endl;
    }
  }
}

template <typename CF, typename DG>
bool BiOptSat<CF, DG>::addDecConstr(const std::vector<DG> &costs, std::vector<Lit> &witness) {
  if (!solver.logger) return rObjectives[1].addHard(costs[1] - 1).first;

  auto pd_cut = MO::pdCut(costs, rObjectives);
  MO::logPdCut(pd_cut, witness);
  auto pd_cut_e = pd_cut.toExpanded(solver.cePools);
  pd_cut_e->getProof().rule = ProofRule::RUP;
#ifndef NDEBUG
  std::stringstream comment;
  comment << "BiOptSat preliminary PD cut (" << costs[0] << "," << costs[1] << ")";
  solver.logger->logComment(comment.str());
#endif
  ID cut_id = pd_cut_e->logRup();

  DG factor = rObjectives[1].objMaxVal + 1 - costs[1] - rObjectives[1].objective->getRhs();
  std::stringstream proof{};
  proof << cut_id << " " << factor << " * " << rObjectives[1].getPId(costs[1] - 1) << " + ";

  if (pd_cut.terms.size() > 1) {
    // shorten PD cut with increasing lower bound
    proof << incrLb << " " << factor << " * + ";
  }

  return rObjectives[1].addHard(costs[1] - 1, proof.str()).first;
}

INSTANTIATE_CLASS(BiOptSat)

}  // namespace run
}  // namespace rs
