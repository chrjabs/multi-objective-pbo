/***********************************************************************
Copyright (c) 2014-2020, Jan Elffers
Copyright (c) 2019-2021, Jo Devriendt
Copyright (c) 2020-2021, Stephan Gocht
Copyright (c) 2014-2021, Jakob Nordström

Parts of the code were copied or adapted from MiniSat.

MiniSat -- Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson
           Copyright (c) 2007-2010  Niklas Sorensson

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
***********************************************************************/

#pragma once

#include <memory>
#include <sstream>
#include "ConstrSimple.hpp"
#include "ProofBuffer.hpp"
#include "SolverStructs.hpp"
#include "globals.hpp"
#include "typedefs.hpp"

namespace rs {

struct IntSet;
struct Logger;

enum class AssertionStatus { NONASSERTING, ASSERTING, FALSIFIED };

// shared_ptr-like wrapper around ConstrExp, ensuring it gets released back to the pool when no longer needed.
template <typename CE>
struct CePtr {
  CE* ce;

  // default constructor
  CePtr() : ce(nullptr) {}
  // regular constructor
  CePtr(CE* c) : ce(c) {
    if (ce) ce->increaseUsage();
  }
  // copy constructor
  CePtr(const CePtr<CE>& other) : ce{other.ce} {
    if (ce) ce->increaseUsage();
  }
  // copy constructor allowing for polymorphism
  template <typename T, typename = std::enable_if_t<std::is_convertible_v<T&, CE&>>>
  CePtr(const CePtr<T>& other) : ce{other.ce} {
    if (ce) ce->increaseUsage();
  }
  // move constructor
  CePtr(CePtr<CE>&& other) : ce{other.ce} { other.ce = nullptr; }
  // move constructor allowing for polymorphism
  template <typename T, typename = std::enable_if_t<std::is_convertible_v<T&, CE&>>>
  CePtr(CePtr<T>&& other) : ce{other.ce} {
    other.ce = nullptr;
  }
  // destructor
  ~CePtr() {
    if (ce) ce->decreaseUsage();
  }
  // assignment operator
  CePtr<CE>& operator=(const CePtr<CE>& other) {
    if (this == &other) return *this;
    if (ce) ce->decreaseUsage();
    ce = other.ce;
    if (ce) ce->increaseUsage();
    return *this;
  }
  // move assignment operator
  CePtr<CE>& operator=(CePtr<CE>&& other) {
    if (this == &other) return *this;
    if (ce) ce->decreaseUsage();
    ce = other.ce;
    other.ce = nullptr;
    return *this;
  }

  CE& operator*() const { return *ce; }
  CE* operator->() const { return ce; }
  explicit operator bool() const { return ce; }
};

struct ConstraintAllocator;
class ConstrExpPools;
class Solver;

struct ConstrExpSuper {
  std::vector<Var> vars;
  Origin orig = Origin::UNKNOWN;

  virtual ~ConstrExpSuper() {}

  virtual void increaseUsage() = 0;
  virtual void decreaseUsage() = 0;

  virtual void copyTo(Ce32 ce) const = 0;
  virtual void copyTo(Ce64 ce) const = 0;
  virtual void copyTo(Ce96 ce) const = 0;
  virtual void copyTo(Ce128 ce) const = 0;
  virtual void copyTo(CeArb ce) const = 0;

  virtual CeSuper clone(ConstrExpPools& ce) const = 0;
  virtual CRef toConstr(ConstraintAllocator& ca, bool locked, ID id) const = 0;
  virtual std::unique_ptr<ConstrSimpleSuper> toSimple() const = 0;

  virtual void resize(size_t s) = 0;
  virtual void resetProof(ID proofID) = 0;
  virtual ProofBuffer& getProof() = 0;
  virtual void initializeLogging(std::shared_ptr<Logger>& l) = 0;
  virtual void stopLogging() = 0;
  virtual bool isReset() const = 0;
  virtual void reset() = 0;

  virtual Lit getLit(Lit l) const = 0;
  virtual bool hasLit(Lit l) const = 0;

  virtual void weaken(Var v) = 0;
  virtual void weakenLast() = 0;

  virtual void harden(Lit l) = 0;  // assume a literal will always be true

  virtual bool hasNegativeSlack(const Solver& solver) const = 0;
  virtual bool hasNegativeSlack(const IntSet& assumptions) const = 0;
  virtual bool isTautology() const = 0;
  virtual bool isInconsistency() const = 0;
  virtual bool isSatisfied(const Solver& solver) const = 0;

  virtual void removeUnitsAndZeroes(const Solver& solver, bool doSaturation = true) = 0;
  virtual bool hasNoUnits(const Solver& solver) const = 0;
  virtual void removeZeroes() = 0;
  virtual void removeZeroesWithCallback(const std::function<void(Var v)>& callback) = 0;
  virtual bool hasNoZeroes() const = 0;

  virtual void saturate(const std::vector<Var>& vs, bool check = true) = 0;
  virtual void saturate(bool check = true) = 0;
  virtual bool isSaturated() const = 0;
  virtual void saturateAndFixOverflow(const Solver& solver, bool fullWeakening, int bitOverflow, int bitReduce,
                                      Lit asserting) = 0;
  virtual void saturateAndFixOverflowRational(const std::vector<double>& lpSolution) = 0;
  virtual bool fitsInDouble() const = 0;
  virtual bool largestCoefFitsIn(int bits) const = 0;

  virtual void weakenDivideRound(const Solver& solver, Lit l, bool slackdiv, bool fullWeakening) = 0;

  virtual bool divideByGCD() = 0;
  virtual void postProcess(const Solver& solver, bool sortFirst, Stats& sts) = 0;
  virtual AssertionStatus isAssertingBefore(const Solver& solver, int lvl) const = 0;
  virtual int getAssertionLevel(const Solver& solver) const = 0;
  virtual void heuristicWeakening(const Solver& solver, Stats& sts) = 0;

  virtual bool simplifyToCardinality(bool equivalencePreserving) = 0;
  virtual void simplifyToClause() = 0;
  virtual bool isCardinality() const = 0;
  virtual int getCardinalityDegree() const = 0;
  virtual int getCardinalityDegreeWithZeroes() = 0;
  virtual void simplifyToMinLengthCardinality() = 0;
  virtual bool isClause() const = 0;
  virtual void sortInDecreasingCoefOrder(const std::function<bool(Var, Var)>& tiebreaker = [](Var, Var) {
    return false;
  }) = 0;
  virtual bool isSortedInDecreasingCoefOrder() const = 0;
  virtual void sort(const std::function<bool(Var, Var)>& comp) = 0;

  virtual ID logCG() = 0;
  virtual ID logRup() = 0;
  virtual ID logImplied() = 0;
  virtual ID logAsAssumption() = 0;
  virtual void logAssertEquals(ID id) = 0;
  virtual ID logInput() = 0;
  virtual ID logProofLine() = 0;
  virtual ID logRedundant() = 0;
  virtual ID logProofLineWithInfo([[maybe_unused]] std::string&& info, [[maybe_unused]] const Stats& sts) = 0;
  virtual void logUnit(const Solver& solver, Var v_unit, const Stats& sts) = 0;
  virtual ID logInconsistency(const Solver& solver, const Stats& sts) = 0;

  virtual void toStreamAsOPB(std::ostream& o) const = 0;
  virtual void toStreamWithAssignment(std::ostream& o, const Solver& solver) const = 0;

  virtual void resolveWith(Ce32 c, Lit l, IntSet* actSet, const Solver& solver) = 0;
  virtual void resolveWith(Ce64 c, Lit l, IntSet* actSet, const Solver& solver) = 0;
  virtual void resolveWith(Ce96 c, Lit l, IntSet* actSet, const Solver& solver) = 0;
  virtual void resolveWith(Ce128 c, Lit l, IntSet* actSet, const Solver& solver) = 0;
  virtual void resolveWith(CeArb c, Lit l, IntSet* actSet, const Solver& solver) = 0;
  virtual void resolveWith(const Clause& c, Lit l, IntSet* actSet, const Solver& solver) = 0;
  virtual void resolveWith(const Cardinality& c, Lit l, IntSet* actSet, const Solver& solver) = 0;
};

template <typename CE>
class ConstrExpPool;

template <typename SMALL, typename LARGE>  // LARGE should be able to fit sums of SMALL
struct ConstrExp final : public ConstrExpSuper {
 private:
  ConstrExpPool<ConstrExp<SMALL, LARGE>>& pool;
  long long usageCount = 0;

 public:
  LARGE degree = 0;
  LARGE rhs = 0;
  std::vector<SMALL> coefs;
  std::vector<bool> used;
  ProofBuffer proof;
  std::shared_ptr<Logger> plogger;

 private:
  void remove(Var v);
  bool increasesSlack(const Solver& solver, Var v) const;
  LARGE calcDegree() const;
  LARGE calcRhs() const;
  bool falsified(const Solver& solver, Var v) const;
  void logIfUnit(Lit l, const SMALL& c, const Solver& solver);

 public:
  ConstrExp(ConstrExpPool<ConstrExp<SMALL, LARGE>>& cep);
  void increaseUsage();
  void decreaseUsage();
  bool testConstraint() const;

  void copyTo(Ce32 ce) const { copyTo_(ce); }
  void copyTo(Ce64 ce) const { copyTo_(ce); }
  void copyTo(Ce96 ce) const { copyTo_(ce); }
  void copyTo(Ce128 ce) const { copyTo_(ce); }
  void copyTo(CeArb ce) const { copyTo_(ce); }

  CeSuper clone(ConstrExpPools& ce) const;
  CRef toConstr(ConstraintAllocator& ca, bool locked, ID id) const;
  std::unique_ptr<ConstrSimpleSuper> toSimple() const;

  void resize(size_t s);
  void resetProof(ID proofID) { proof.reset(proofID); }
  ProofBuffer& getProof() { return proof; }
  void initializeLogging(std::shared_ptr<Logger>& l);
  void stopLogging();

  bool isReset() const;
  void reset();

  LARGE getRhs() const;
  LARGE getDegree() const;
  SMALL getCoef(Lit l) const;
  SMALL getLargestCoef() const;
  SMALL getSmallestCoef() const;
  LARGE getCutoffVal() const;
  Lit getLit(Lit l) const;
  bool hasLit(Lit l) const;

  void setDegree(LARGE d) {
    degree = d;
    rhs = calcRhs();
  }

  void addRhs(const LARGE& r);
  void addLhs(const SMALL& cf, Lit l);  // TODO: Term?
  void weaken(const SMALL& m, Var v);
  void weaken(Var v);
  void weakenLast();

  void harden(Lit l);  // leads to a zero in the constraint, remove it later

  LARGE getSlack(const Solver& solver) const;
  bool hasNegativeSlack(const Solver& solver) const;
  LARGE getSlack(const IntSet& assumptions) const;
  bool hasNegativeSlack(const IntSet& assumptions) const;
  bool isTautology() const;
  bool isInconsistency() const;
  bool isSatisfied(const Solver& solver) const;

  /**
   * @brief Remove literals in reason with coeff 0 and weaken away literals that are unit constraints in our database.
   *
   * @post: preserves order of vars
   */
  void removeUnitsAndZeroes(const Solver& solver, bool doSaturation = true);
  bool hasNoUnits(const Solver& solver) const;
  // @post: preserves order of vars
  void removeZeroes();
  // @post: preserves order of vars
  void removeZeroesWithCallback(const std::function<void(Var v)>& callback);
  bool hasNoZeroes() const;

  /**
   * @brief Saturate the constraint.
   *
   * @tparam SMALL Coefficient type.
   * @tparam LARGE Degree type.
   * @param vs Vector of variable in the constraint.
   * @param check If true, check if constraint is already saturated and return directly if this is the case.
   *
   * @post preserves order of vars
   */
  void saturate(const std::vector<Var>& vs, bool check = true);
  /**
   * @brief Saturate the constraint.
   *
   * @tparam SMALL Coefficient type.
   * @tparam LARGE Degree type.
   * @param check If true, check if constraint is already saturated and return if this is the case.
   */
  void saturate(bool check = true);
  bool isSaturated() const;
  /**
   * @brief Saturate after resolution step and prevent overflow.
   *
   * @post saturated
   * @post nothing else if bitOverflow == 0
   * @post the largest coefficient is less than 2^bitOverflow
   * @post the degree and rhs are less than 2^bitOverflow * INF
   * @post if overflow happened, all division until 2^bitReduce happened
   * @post the constraint remains conflicting or propagating on asserting
   */
  void saturateAndFixOverflow(const Solver& solver, bool fullWeakening, int bitOverflow, int bitReduce, Lit asserting);
  /*
   * Fixes overflow for rationals
   * @post: saturated
   * @post: none of the coefficients, degree, or rhs exceed INFLPINT
   */
  void saturateAndFixOverflowRational(const std::vector<double>& lpSolution);
  bool fitsInDouble() const;
  bool largestCoefFitsIn(int bits) const;

  template <typename S, typename L>
  void addUp(CePtr<ConstrExp<S, L>> c, const SMALL& cmult = 1, const SMALL& thismult = 1) {
    assert(cmult >= 1);
    assert(thismult >= 1);
    if (plogger) proof.addConstraint(thismult, c->proof, cmult);
    if (thismult > 1) {
      degree *= thismult;
      rhs *= thismult;
      for (Var v : vars) coefs[v] *= thismult;
    }
    rhs += static_cast<LARGE>(cmult) * static_cast<LARGE>(c->rhs);
    degree += static_cast<LARGE>(cmult) * static_cast<LARGE>(c->degree);
    for (Var v : c->vars) {
      assert(v < (Var)coefs.size());
      assert(v > 0);
      SMALL val = cmult * static_cast<SMALL>(c->coefs[v]);
      if (!used[v]) {
        assert(coefs[v] == 0);
        vars.push_back(v);
        coefs[v] = val;
        used[v] = true;
      } else {
        if ((coefs[v] < 0) != (val < 0)) degree -= std::min(aux::abs(coefs[v]), aux::abs(val));
        coefs[v] += val;
      }
    }
  }

  void invert();
  void multiply(const SMALL& m);
  void divide(const LARGE& d);
  /**
   * @brief Divides the current constraint by `d` and rounds up the coefficients and degree.
   *
   * @tparam SMALL Coefficient type.
   * @tparam LARGE Degree type.
   * @param d Divisor for the division.
   */
  void divideRoundUp(const LARGE& d);

  /**
   * @brief Initial RoundingSAT refinement function. Choose a divisor and make all coefficient of non-falsified literals
   * divisible by the divisor.
   *
   * @tparam SMALL Coefficient type.
   * @tparam LARGE Degree type.
   * @param solver A reference to the solver object
   * @param l Literal to resolve over.
   * @param slackdiv If true, the divisor is reason constraint slack + 1. If false, divide by the literal coefficient in
   * the reason constraint.
   * @param fullWeakening If true, weaken literals completely. If false, weaken literals to biggest dividable
   * coefficient.
   */
  void weakenDivideRound(const Solver& solver, Lit l, bool slackdiv, bool fullWeakening);

  /**
   * @brief Weaken the constraint such that each coefficient of the non falsified literal in the constraint are
   * divisible by `div`.
   *
   * @param solver A reference to the solver object
   * @param div Divisor to divide the constraint by in the resolution step.
   * @param fullWeakening If true, weaken literals completely. If false, weaken literals to biggest dividable
   * coefficient.
   * @param asserting The literal to resolve over in the resolution step.
   */
  void weakenNonDivisibleNonFalsified(const Solver& solver, const LARGE& div, bool fullWeakening, Lit asserting);
  void applyMIR(const LARGE& d, std::function<Lit(Var)> toLit);

  SMALL getGCD();
  void divideBy(SMALL div);
  bool divideByGCD();
  // NOTE: only equivalence preserving operations!
  void postProcess(const Solver& solver, bool sortFirst, Stats& sts);
  AssertionStatus isAssertingBefore(const Solver& solver, int lvl) const;
  // @return: earliest decision level that propagates a variable
  int getAssertionLevel(const Solver& solver) const;
  // @post: preserves order after removeZeroes()
  void weakenNonImplied(const Solver& solver, const LARGE& slack, Stats& sts);
  /**
   * @brief Weaken literals that are falsified and have a small enough coefficient at the current level.
   *
   * @post: preserves order after removeZeroes()
   *
   * @todo: return modified slack?
   */
  bool weakenNonImplying(const Solver& solver, const SMALL& propCoef, const LARGE& slack, Stats& sts);
  // @post: preserves order after removeZeroes()
  void heuristicWeakening(const Solver& solver, Stats& sts);

  // @post: preserves order
  template <typename T>
  void weakenSmalls(const T& limit) {
    for (Var v : vars)
      if (aux::abs(coefs[v]) < limit) weaken(v);
    saturate();
  }

  LARGE absCoeffSum() const;

  // @post: preserves order of vars
  bool simplifyToCardinality(bool equivalencePreserving);
  bool isCardinality() const;
  int getCardinalityDegree() const;
  int getCardinalityDegreeWithZeroes();
  void simplifyToMinLengthCardinality();
  void simplifyToClause();
  bool isClause() const;

  void sortInDecreasingCoefOrder(const std::function<bool(Var, Var)>& tiebreaker = [](Var, Var) { return false; });
  bool isSortedInDecreasingCoefOrder() const;
  void sort(const std::function<bool(Var, Var)>& comp);

  ID logAsAssumption();
  ID logInput();
  ID logCG();
  ID logRup();
  ID logImplied();
  void logAssertEquals(ID id);
  ID logProofLine();
  ID logRedundant();
  ID logProofLineWithInfo([[maybe_unused]] std::string&& info, [[maybe_unused]] const Stats& sts);
  // @pre: reducible to unit over v
  void logUnit(const Solver& solver, Var v_unit, const Stats& sts);
  ID logInconsistency(const Solver& solver, const Stats& sts);

  void toStreamAsOPB(std::ostream& o) const;
  void toStreamWithAssignment(std::ostream& o, const Solver& solver) const;

  void resolveWith(Ce32 c, Lit l, IntSet* actSet, const Solver& solver);
  void resolveWith(Ce64 c, Lit l, IntSet* actSet, const Solver& solver);
  void resolveWith(Ce96 c, Lit l, IntSet* actSet, const Solver& solver);
  void resolveWith(Ce128 c, Lit l, IntSet* actSet, const Solver& solver);
  void resolveWith(CeArb c, Lit l, IntSet* actSet, const Solver& solver);
  void resolveWith(const Clause& c, Lit l, IntSet* actSet, const Solver& solver);
  void resolveWith(const Cardinality& c, Lit l, IntSet* actSet, const Solver& solver);

  ConstrSimple<SMALL, LARGE> toOwnedSimple() const {
    ConstrSimple<SMALL, LARGE> result{};
    result.rhs = rhs;
    result.terms.reserve(vars.size());
    for (Var v : vars)
      if (coefs[v] != 0) result.terms.emplace_back(coefs[v], v);
    if (plogger) result.proofLine = proof.buffer.str();
    result.orig = orig;
    return result;
  }

 private:
  void addUsedLitsToActiveSet(IntSet* actSet, Lit l, const Solver& solver);

  /**
   * @brief Refine the reason constraint such that the slack of the resolvent with the conflict constraint over `l` is
   * negative.
   *
   * @tparam CF Coefficient type.
   * @tparam DG Degree type.
   * @param reason Reason constraint.
   * @param l Literal to resolve over.
   * @param solver A reference to the solver object
   */
  template <typename CF, typename DG>
  void refineConstrToNegativeSlackResolvent(CePtr<ConstrExp<CF, DG>> reason, Lit l, const Solver& solver) {
    reason->removeUnitsAndZeroes(solver);
    if (options.weakenNonImplying)
      reason->weakenNonImplying(solver, reason->getCoef(l), reason->getSlack(solver), stats);
    reason->saturateAndFixOverflow(solver, (bool)options.weakenFull, options.bitsOverflow.get(),
                                   options.bitsReduced.get(), l);

    reason->weakenDivideRound(solver, l, (bool)options.slackdiv, (bool)options.weakenFull);
  }

  template <typename CF, typename DG>
  void generalizedResolution(CePtr<ConstrExp<CF, DG>> reason, Lit l) {
    SMALL reason_coef_l = static_cast<SMALL>(reason->getCoef(l));  // NOTE: SMALL >= CF
    SMALL confl_coef_l = getCoef(-l);
    SMALL gcd_coef_l = aux::gcd(reason_coef_l, confl_coef_l);
    addUp(reason, confl_coef_l / gcd_coef_l, reason_coef_l / gcd_coef_l);
  }

  template <typename CF, typename DG>
  void genericResolve(CePtr<ConstrExp<CF, DG>> reason, Lit l, IntSet* actSet, const Solver& solver) {
    assert(getCoef(-l) > 0);
    stats.NADDEDLITERALS += reason->vars.size();

    refineConstrToNegativeSlackResolvent(reason, l, solver);

    // Add used variables to active set.
    if (actSet != nullptr) {
      for (Var v : reason->vars) {
        Lit l = reason->getLit(v);
        addUsedLitsToActiveSet(actSet, l, solver);
      }
    }

    generalizedResolution(reason, l);
    saturateAndFixOverflow(solver, (bool)options.weakenFull, options.bitsOverflow.get(), options.bitsReduced.get(), 0);
    assert(getCoef(-l) == 0);
    assert(hasNegativeSlack(solver));
  }

  template <typename S, typename L>
  void copyTo_(CePtr<ConstrExp<S, L>> out) const {
    // TODO: assert whether S/L can fit SMALL/LARGE? Not always possible.
    assert(out->isReset());
    out->degree = static_cast<L>(degree);
    out->rhs = static_cast<L>(rhs);
    out->orig = orig;
    out->vars = vars;
    assert(out->coefs.size() == coefs.size());
    for (Var v : vars) {
      out->coefs[v] = static_cast<S>(coefs[v]);
      assert(used[v] == true);
      assert(out->used[v] == false);
      out->used[v] = true;
    }
    if (plogger) {
      out->proof.copyFrom(proof.buffer.str());
    }
  }

  template <typename S, typename L>
  std::unique_ptr<ConstrSimple<S, L>> toSimple_() const {
    std::unique_ptr<ConstrSimple<S, L>> result = std::make_unique<ConstrSimple<S, L>>();
    result->rhs = static_cast<L>(rhs);
    result->terms.reserve(vars.size());
    for (Var v : vars)
      if (coefs[v] != 0) result->terms.emplace_back(static_cast<S>(coefs[v]), v);
    if (plogger) result->proofLine = proof.buffer.str();
    result->orig = orig;
    return result;
  }
};

template <typename S, typename L>
std::ostream& operator<<(std::ostream& o, const ConstrExp<S, L>& C) {
  std::vector<Var> vars = C.vars;
  std::sort(vars.begin(), vars.end(), [](Var v1, Var v2) { return v1 < v2; });
  for (Var v : vars) {
    Lit l = C.coefs[v] < 0 ? -v : v;
    o << C.getCoef(l) << "x" << l << " ";
  }
  o << ">= " << C.degree;
  return o;
}

}  // namespace rs
