#include "pminimal.hpp"
#include <ranges>
#include "Logger.hpp"
#include "auxiliary.hpp"

namespace rs {
namespace run {

template <typename CF, typename DG>
void PMinimal<CF, DG>::solve() {
  bool unsat = false;
  if (options.moCoreBoosting && MO::coreBoosting()) unsat = true;
  std::vector<Lit> solution;
  std::vector<DG> costs;
  while (!unsat) {
    solver.clearAssumptions();
    SolveAnswer out = solveFully([&] { return solver.solve(); }, stats.SOLVETIME);
    switch (out.state) {
      case SolveState::UNSAT:
        break;
      case SolveState::SAT:
        ++stats.NSOLS;
        solution = out.solution;
        costs = MO::computeCosts(solution);
        unsat = pMinimize(solution, costs);
        MO::addSolution(solution, costs);
        if (unsat) break;
        continue;
      default:
        assert(false);
        break;
    }
    break;
  }
}

template <typename CF, typename DG>
bool PMinimal<CF, DG>::pMinimize(std::vector<Lit> &solution, std::vector<DG> &costs) {
  IntSet assumps{};
  while (true) {
    auto bl_clause = MO::pdCut(costs, rObjectives);
    MO::logPdCut(bl_clause, solution);
    if (solver.addConstraint(bl_clause, Origin::PMINIMAL).second == ID_Unsat) return true;
    assumps.clear();
    if (enforceDominating(costs, assumps)) return true;
    SolveAnswer out = solveFully(
        [&] {
          solver.setAssumptions(assumps);
          return solver.solve();
        },
        stats.SOLVETIME);
    switch (out.state) {
      case SolveState::UNSAT:
        return true;
      case SolveState::INCONSISTENT:
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
      std::cout << "c p-minimal candidate (";
      for (size_t idx = 0; idx < costs.size(); idx++)
        std::cout << costs[idx] << ((idx < costs.size() - 1) ? ", " : ")");
      std::cout << std::endl;
    }
  }
  return false;
}

template <typename CF, typename DG>
bool PMinimal<CF, DG>::enforceDominating(const std::vector<DG> &costs, IntSet &assumps) {
  for (size_t idx = 0; idx < costs.size(); idx++) {
    const DG &cost = costs[idx];
    if (cost >= rObjectives[idx].objMaxVal) continue;
    Lit out = rObjectives[idx].getOutput(cost);
    if (!out) return true;
    assumps.add(-out);
  }
  return false;
}

INSTANTIATE_CLASS(PMinimal)

}  // namespace run
}  // namespace rs
