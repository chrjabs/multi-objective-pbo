#include "paretopk.hpp"
#include "auxiliary.hpp"
#include "coreguided.hpp"

namespace rs {
namespace run {

template <typename CF, typename DG>
void PareTopK<CF, DG>::solve() {
  if (options.moCoreBoosting && MO::coreBoosting()) return;
  bool cgCoreUpper = (bool)options.cgCoreUpper;
  options.cgCoreUpper.set(false);
  bool cgHardening = (bool)options.cgHardening;
  options.cgHardening.set(false);
  std::vector<Lit> solution;
  std::vector<DG> costs;
  while (true) {
    solver.clearAssumptions();
    // CG search over joint objective
    if (cg::coreGuidedReformulateWithLazies(reformJointObj, origJointObj, solution, lazies)) break;
    if (options.verbosity.get() > 2) std::cout << "c Joint objective value: " << -reformJointObj->getRhs() << std::endl;
    // add solution to pareto front
    costs = MO::computeCosts(solution);
    MO::addSolution(solution, costs);
    // block solution with P-min cut
    auto bl_clause = MO::pdCut(costs, PMin::rObjectives);
    if (solver.addConstraint(bl_clause, Origin::PMINIMAL).second == ID_Unsat) break;
  }
  options.cgCoreUpper.set(cgCoreUpper);
  options.cgHardening.set(cgHardening);
}

INSTANTIATE_CLASS(PareTopK)

}  // namespace run
}  // namespace rs
