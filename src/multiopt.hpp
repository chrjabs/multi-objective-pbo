#pragma once

#include "ConstrSimple.hpp"
#include "Solver.hpp"
#include "coreguided.hpp"
#include "run.hpp"
#include "typedefs.hpp"

namespace rs {
namespace run {

template <typename DG>
struct NonDom {
  std::vector<DG> objValues;
  std::vector<Lit> solution;

  NonDom(const std::vector<DG> &objValues, const std::vector<Lit> &solution)
      : objValues(objValues), solution(solution) {}
};

class MultiOptBase {
 public:
  virtual void solve() = 0;
  virtual void printParetoFront(bool printSolutions = false) const = 0;
  virtual void printStats() const = 0;

  virtual ~MultiOptBase() {}
};

template <typename CF, typename DG>
struct ReifiedObj;

template <typename CF, typename DG>
class MultiOpt : virtual public MultiOptBase {
  using CE = CePtr<ConstrExp<CF, DG>>;
  using CS = ConstrSimple<CF, DG>;
  using ND = NonDom<DG>;
  using RObj = ReifiedObj<CF, DG>;

 public:
  const std::vector<ND> &getParetoFront() { return paretoFront; }
  bool coreBoosting();
  void printParetoFront(bool printSolutions = false) const final override;
  virtual void printStats() const override;

  void logPdCut(ConstrSimple32 &cut, std::vector<Lit> &witness);
  ConstrSimple32 pdCut(const std::vector<DG> &costs, std::vector<RObj> &rObjectives);

  MultiOpt(const std::vector<CeArb> objs) {
    origObjs.reserve(objs.size());
    reformObjs.reserve(objs.size());
    for (auto obj : objs) {
      reformObjs.push_back(solver.cePools.take<CF, DG>());
      obj->copyTo(reformObjs.back());
      origObjs.push_back(reformObjs.back()->toOwnedSimple());
    }
    logMoOrder();
  }

  virtual ~MultiOpt() {}

 protected:
  std::vector<CS> origObjs{};
  std::vector<CE> reformObjs;
  std::vector<ND> paretoFront{};

  std::vector<DG> computeCosts(const std::vector<Lit> &solution) const;
  void addSolution(const std::vector<Lit> &solution, const std::vector<DG> &costs = {});
  void logMoOrder();
};

template <typename DG>
struct RObjData {
  DG bound;
  Lit output;
  ID pid;
};

template <typename CF, typename DG>
struct ReifiedObj {
  using Obj = CePtr<ConstrExp<CF, DG>>;

  const DG objMaxVal;
  const Obj objective{};
  // Sorted list of already existing outputs
  std::vector<RObjData<DG>> outMap{};

  ReifiedObj(const Obj objective) : objMaxVal(cg::maxObjVal(objective)), objective(objective) {};

  RObjData<DG> newReification(DG bound, std::vector<RObjData<DG>>::const_iterator iter) {
    Lit out = solver.getNbVars() + 1;
    solver.setNbVars(out);
    Obj reified = solver.cePools.take<CF, DG>();
    objective->copyTo(reified);
    reified->addRhs(bound);
    reified->invert();
    reified->addLhs(static_cast<CF>(objMaxVal - bound), out);
    if (solver.logger) {
      reified->getProof().rule = ProofRule::REDUNDANT;
      reified->getProof().buffer.clear();
      reified->getProof().buffer.str(std::string());
      reified->getProof().buffer << "x" << out << " -> 1";
    }
    auto ids = solver.addConstraint(reified, Origin::REIFIEDOBJ);
    RObjData data = RObjData{.bound = bound, .output = ids.second == ID_Unsat ? 0 : out, .pid = ids.first};
    outMap.insert(iter, data);
    return data;
  }

  // Gets an `output`, i.e., a reified literal `a` where `-a -> (Obj <= bound)`.
  // If no reified output for the given bound exists yet, adds the reification to the solver.
  // If the solver detects UNSAT while adding the reification, will return 0.
  Lit getOutput(DG bound) {
    assert(bound >= -objective->getDegree());
    assert(bound < objMaxVal);
    auto iter = outMap.begin();
    for (; iter != outMap.end(); iter++) {
      if (iter->bound == bound) return iter->output;
      if (iter->bound > bound) break;
    }
    // Create new reification
    RObjData<DG> data = newReification(bound, iter);
    return data.output;
  }

  // Gets the proof ID of a given `output`
  ID getPId(DG bound) {
    assert(bound >= -objective->getDegree());
    assert(bound < objMaxVal);
    auto iter = outMap.begin();
    for (; iter != outMap.end(); iter++) {
      if (iter->bound == bound) return iter->pid;
      if (iter->bound > bound) break;
    }
    // Create new reification
    RObjData<DG> data = newReification(bound, iter);
    return data.pid;
  }

  // Adds a hard constraint enforcing a bound. Returns true if the solver detected UNSAT.
  std::pair<bool, ID> addHard(DG bound, std::string proof = "") {
    assert(bound >= -objective->getDegree());
    assert(bound < objMaxVal);
    // Create constraint
    Obj constr = solver.cePools.take<CF, DG>();
    objective->copyTo(constr);
    constr->addRhs(bound);
    constr->invert();
    constr->getProof().buffer.str(proof);
    auto ids = solver.addConstraint(constr, Origin::UPPERBOUND);
    return std::make_pair(ids.second == ID_Unsat, ids.first);
  }
};

}  // namespace run
}  // namespace rs
