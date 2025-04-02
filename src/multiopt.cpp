#include "multiopt.hpp"
#include <ranges>
#include <sstream>
#include <unordered_set>
#include "Logger.hpp"
#include "auxiliary.hpp"
#include "coreguided.hpp"

namespace rs {
namespace run {

template <typename SMALL, typename LARGE>
CePtr<ConstrExp<SMALL, LARGE>> expand(const ConstrSimple<SMALL, LARGE> &simp) {
  CePtr<ConstrExp<SMALL, LARGE>> exp = solver.cePools.take<SMALL, LARGE>();
  exp->addRhs(simp.rhs);
  for (const Term<SMALL> &t : simp.terms) exp->addLhs(t.c, t.l);
  exp->orig = simp.orig;
  if (exp->plogger) {
    exp->proofBuffer.str(std::string());
    exp->proofBuffer << simp.proofLine;
  }
  return exp;
}

template <typename CF, typename DG>
bool MultiOpt<CF, DG>::coreBoosting() {
  bool cgCoreUpper = (bool)options.cgCoreUpper;
  options.cgCoreUpper.set(false);
  bool cgHardening = (bool)options.cgHardening;
  options.cgHardening.set(false);
  for (size_t idx = 0; idx < reformObjs.size(); idx++) {
    if (options.verbosity.get() >= 1) std::cout << "c core boosting objective " << idx + 1 << std::endl;
    std::vector<Lit> sol{};  // ignored for now
    if (cg::coreGuidedReformulate(reformObjs[idx], origObjs[idx], sol)) return true;
    if (options.verbosity.get() >= 2)
      std::cout << "c derived lower bound " << -reformObjs[idx]->getDegree() << " for objective " << idx + 1
                << std::endl;
  }
  if (options.verbosity.get() >= 1) {
    std::cout << "c cb derived ideal point (";
    for (size_t idx = 0; idx < reformObjs.size(); idx++)
      std::cout << -reformObjs[idx]->getDegree() << ((idx < reformObjs.size() - 1) ? ", " : ")");
    std::cout << std::endl;
  }
  options.cgCoreUpper.set(cgCoreUpper);
  options.cgHardening.set(cgHardening);
  return false;
}

template <typename CF, typename DG>
std::vector<DG> MultiOpt<CF, DG>::computeCosts(const std::vector<Lit> &solution) const {
  std::vector<DG> costs{};
  costs.reserve(origObjs.size());
  for (size_t idx = 0; idx < origObjs.size(); idx++) {
    costs.push_back(cg::computeCost(origObjs[idx], solution));
#ifndef NDEBUG
    DG reformVal = -reformObjs[idx]->getRhs();
    for (Var v : reformObjs[idx]->vars) reformVal += reformObjs[idx]->coefs[v] * (int)(solution[v] > 0);
    assert(reformVal >= costs.back());
#endif
  }
  return costs;
}

template <typename CF, typename DG>
void MultiOpt<CF, DG>::addSolution(const std::vector<Lit> &solution, [[maybe_unused]] const std::vector<DG> &costs) {
  auto true_costs = computeCosts(solution);
#ifndef NDEBUG
  if (!costs.empty()) {
    for (size_t idx = 0; idx < true_costs.size() && idx < costs.size(); idx++) assert(true_costs[idx] == costs[idx]);
  }
#endif
  paretoFront.emplace_back(true_costs, solution);
  if (options.verbosity.get() >= 1) {
    std::cout << "c found non-dominated point (";
    for (size_t i = 0; i < true_costs.size(); i++)
      std::cout << true_costs[i] << ((i < true_costs.size() - 1) ? ", " : ")");
    std::cout << std::endl;
  }
}

template <typename CF, typename DG>
void MultiOpt<CF, DG>::printParetoFront(bool printSolutions) const {
  std::cout << "c ===[ Pareto Front ]===========================================================" << std::endl;
  for (auto &ndom : paretoFront) {
    std::cout << "c ---[ Non-dominated point ]----------------------------------------------------" << std::endl;
    std::cout << "o";
    for (DG cost : ndom.objValues) std::cout << " " << cost;
    std::cout << std::endl;
    if (printSolutions) quit::printSol(ndom.solution);
  }
}

template <typename CF, typename DG>
void MultiOpt<CF, DG>::printStats() const {
  std::cout << "c ===[ Statistics ]=============================================================" << std::endl;
  std::cout << "c n-non-dominated: " << paretoFront.size() << std::endl;
}

template <typename CF, typename DG>
void MultiOpt<CF, DG>::logMoOrder() {
  if (!solver.logger) return;

  // prepare data
  std::unordered_set<Var> objvars{};
  std::stringstream definition{};
  for (auto &obj : origObjs) {
    definition << "    ";
    for (auto term : obj.terms) {
      objvars.insert(toVar(term.l));
      definition << -term.c << (term.l < 0 ? " ~u" : " u") << toVar(term.l) << " " << term.c
                 << (term.l < 0 ? " ~v" : " v") << toVar(term.l) << " ";
    }
    definition << " >= 0 ;" << std::endl;
  }

  // write order definition
  solver.logger->proof_out << "def_order pareto" << std::endl;
  solver.logger->proof_out << "  vars" << std::endl;
  solver.logger->proof_out << "    left";
  for (Var var : objvars) solver.logger->proof_out << " u" << var;
  solver.logger->proof_out << std::endl << "    right";
  for (Var var : objvars) solver.logger->proof_out << " v" << var;
  solver.logger->proof_out << std::endl << "    aux" << std::endl;
  solver.logger->proof_out << "  end" << std::endl;
  solver.logger->proof_out << "  def" << std::endl;
  solver.logger->proof_out << definition.str();
  solver.logger->proof_out << "  end" << std::endl;
  solver.logger->proof_out << "  transitivity" << std::endl;
  solver.logger->proof_out << "    vars" << std::endl;
  solver.logger->proof_out << "      fresh_right";
  for (Var var : objvars) solver.logger->proof_out << " w" << var;
  solver.logger->proof_out << std::endl << "    end" << std::endl;
  solver.logger->proof_out << "    proof" << std::endl;
  for (size_t idx = 0; idx < origObjs.size(); idx++) {
    solver.logger->proof_out << "      proofgoal #" << idx + 1 << std::endl;
    solver.logger->proof_out << "        pol " << idx + 1 << " " << origObjs.size() + idx + 1 << " + -1 +" << std::endl;
    solver.logger->proof_out << "      qed -1" << std::endl;
  }
  solver.logger->proof_out << "    qed" << std::endl;
  solver.logger->proof_out << "  end" << std::endl;
  solver.logger->proof_out << "end" << std::endl;
}

template <typename CF, typename DG>
void MultiOpt<CF, DG>::logPdCut(ConstrSimple32 &cut, std::vector<Lit> &witness) {
  if (!solver.logger) return;

  // Make sure that all variables in the pd cut are in the witness, in case they were introduced just now
  std::sort(cut.terms.begin(), cut.terms.end(), [](Term<int> &a, Term<int> &b) { return toVar(a.l) < toVar(b.l); });
  for (const Term<int> &t : cut.terms) {
    if (-t.l > std::abs(witness.back())) {
      assert(-t.l == toVar(witness.back()) + 1);
      witness.push_back(-t.l);
    }
  }

  // STEP 1: map weakly dominated solutions to witness
  CE map_dom_c = solver.cePools.take<CF, DG>();
  map_dom_c->addRhs(witness.size() - 1);
  map_dom_c->getProof().rule = ProofRule::REDUNDANT;
  map_dom_c->getProof().buffer.clear();
  map_dom_c->getProof().buffer.str(std::string());
  for (Lit lit : witness | std::views::drop(1)) {
    map_dom_c->addLhs(1, lit);
    map_dom_c->getProof().buffer << " x" << toVar(lit) << " -> " << (lit > 0);
  }
  for (auto term : cut.terms) {
    map_dom_c->addLhs(witness.size() + 1, term.l);
    if (toVar(witness.back()) < -term.l) map_dom_c->getProof().buffer << " x" << toVar(term.l) << " -> 1";
  }
#ifndef NDEBUG
  run::solver.logger->logComment("PD cut map_dom");
#endif
  ID map_dom = map_dom_c->logRedundant();

  // STEP 2: exclude witness
  ID excl = solver.logger->solx(witness);

  // STEP 3: add justification of PD cut
  std::stringstream hints{};
  hints << map_dom << " " << excl << " ";
  cut.proofLine = hints.str();
}

template <typename CF, typename DG>
ConstrSimple32 MultiOpt<CF, DG>::pdCut(const std::vector<DG> &costs, std::vector<RObj> &rObjectives) {
  ConstrSimple32 bl_clause({}, 1, Origin::PMINIMAL, "");
  for (size_t idx = 0; idx < costs.size(); idx++) {
    const DG &cost = costs[idx];
    if (cost <= -reformObjs[idx]->getDegree()) continue;
    Lit out = rObjectives[idx].getOutput(cost - 1);
    if (!out) {
      // unsat
      ConstrSimple32 unsat({}, 1, Origin::PMINIMAL, "");
      return unsat;
    }
    bl_clause.terms.emplace_back(1, -out);
  }
  return bl_clause;
}

INSTANTIATE_CLASS(MultiOpt)

}  // namespace run
}  // namespace rs
