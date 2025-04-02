#ifndef UNSATSATMSU3MO
#define UNSATSATMSU3MO
// #define PARTIAL
#ifndef PARTIAL

#ifdef SIMP
#include "simp/SimpSolver.h"
#else
#include "core/Solver.h"
#endif

#include "../Encoder.h"

#include "../MaxSAT.h"
#include "./Alg_ServerMO.h"
#include "./Alg_UnsatSatMO.h"
#include "utils/System.h"
#include <utility>
#include <map>
#include <set>

#define MAXDIM 10

namespace openwbo {
  class UnsatSatMSU3MO : public virtual UnsatSatMO {
    
  public:
    UnsatSatMSU3MO(int verb = _VERBOSITY_MINIMAL_, int weight = _WEIGHT_NONE_, int strategy = _WEIGHT_NONE_, 
		   int enc = _CARD_MTOTALIZER_, int pb = _PB_SWC_, 
		   int pbobjf = _PB_GTE_, int apmode = _ap_outvars_, float eps = 1, 
		   int searchStrat=20, float redFact=-1) : 
      BLS( verb,  weight,  strategy, enc,  pb, pbobjf,  apmode,  eps, searchStrat,  redFact),
      UnsatSatMO( verb,  weight,  strategy, enc,  pb, pbobjf,  apmode,  eps, searchStrat,  redFact) {}
    
    bool searchUnsatSatMO() override;
    
  protected:
    vec<Lit> explanation{}; 	// unsat explanation
    bool extendUL(YPoint& ul) override;
    virtual void pushBlockedVars();
    bool popBlockedVars();
    bool buildWorkFormula() override;
  };

  class UnsatSatMSU3IncMO: public UnsatSatMSU3MO, public UnsatSatIncMO{
  public:
    UnsatSatMSU3IncMO(int verb = _VERBOSITY_MINIMAL_, int weight = _WEIGHT_NONE_, int strategy = _WEIGHT_NONE_, 
		      int enc = _CARD_MTOTALIZER_, int pb = _PB_SWC_, 
		      int pbobjf = _PB_GTE_, int apmode = _ap_outvars_, float eps = 1, 
		      int searchStrat=3, float redFact=-1) : 
      BLS( verb,  weight,  strategy, enc,  pb, pbobjf,  apmode,  eps, searchStrat,  redFact),
      UnsatSatMO( verb,  weight,  strategy, enc,  pb, pbobjf,  apmode,  eps, searchStrat,  redFact),
      UnsatSatMSU3MO(verb, weight, strategy, enc, pb, pbobjf, apmode, eps, searchStrat, redFact), 
      UnsatSatIncMO(verb, weight, strategy, enc, pb, pbobjf, apmode, eps, searchStrat, redFact), 
      lowerBound{this}{}
    StatusCode searchAgain() override {  setInitialTime(cpuTime()); searchUnsatSatMO(); return answerType;}
    bool buildWorkFormula() override {return UnsatSatMSU3MO::buildWorkFormula();}
    void checkSols() override;
    void blockStep(const YPoint& yp) override;
    void increment() override {return;}
    bool searchUnsatSatMO() override;

  protected:
    Solution lowerBound;
  };

  class UnsatSatMSU3IncObjMO: public UnsatSatMSU3MO, public UnsatSatIncObjMO{
  public:
    UnsatSatMSU3IncObjMO(int verb = _VERBOSITY_MINIMAL_, int weight = _WEIGHT_NONE_, int strategy = _WEIGHT_NONE_, 
		      int enc = _CARD_MTOTALIZER_, int pb = _PB_SWC_, 
		      int pbobjf = _PB_GTE_, int apmode = _ap_outvars_, float eps = 1, 
		      int searchStrat=3, float redFact=-1) : 
      BLS( verb,  weight,  strategy, enc,  pb, pbobjf,  apmode,  eps, searchStrat,  redFact),
      UnsatSatMO( verb,  weight,  strategy, enc,  pb, pbobjf,  apmode,  eps, searchStrat,  redFact),
      UnsatSatMSU3MO(verb, weight, strategy, enc, pb, pbobjf, apmode, eps, searchStrat, redFact), 
      UnsatSatIncObjMO(verb, weight, strategy, enc, pb, pbobjf, apmode, eps, searchStrat, redFact), 
      lowerBound{this}{}
    StatusCode searchAgain() override {  setInitialTime(cpuTime()); searchUnsatSatMO(); return answerType;}
    bool buildWorkFormula() override {return UnsatSatMSU3MO::buildWorkFormula();}
    void checkSols() override;
    void blockStep(const YPoint& yp) override;
    void increment() override {return;}
    bool searchUnsatSatMO() override;

  protected:
    Solution lowerBound;
  };
}
#endif
#undef PARTIAL
#endif
