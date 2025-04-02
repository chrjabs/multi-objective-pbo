#ifndef UNSATSATSTRATMO
#define UNSATSATSTRATMO
// #define PARTIAL
#include "Alg_StratMO.h"
#ifndef PARTIAL

#ifdef SIMP
#include "simp/SimpSolver.h"
#else
#include "core/Solver.h"
#endif

#include "../Encoder.h"
#include "../Pareto.h"
#include "../MaxSAT.h"
#include "./Alg_ServerMO.h"
#include "./Alg_UnsatSatMO.h"
#include "./Alg_UnsatSatMSU3MO.h"
#include "utils/System.h"
#include <utility>
#include "../partition.h"
#include <algorithm>    // std::max
#include "Alg_MasterMO.h"
#include "Alg_StratMO.h"
#define MAXDIM 10

namespace openwbo {
  using namespace partition;
  
  class UnsatSatStratMO : public virtual UnsatSatMO, public virtual StratMO{
    
  public:
    UnsatSatStratMO(int verb = _VERBOSITY_MINIMAL_, int weight = _WEIGHT_NONE_, int strategy = _WEIGHT_NONE_, 
	       int enc = _CARD_MTOTALIZER_, int pb = _PB_SWC_, 
	       int pbobjf = _PB_GTE_, int apmode = _ap_outvars_, float eps = 1, 
		    int searchStrat=3, int partition_parameter = 15, float redFact=-1) : 
      BLS(verb, weight, strategy, enc, pb, pbobjf, apmode, eps, searchStrat, redFact), UnsatSatMO(verb, weight, strategy, enc, pb, pbobjf, apmode, eps, searchStrat, redFact){
      optim = new UnsatSatMSU3IncObjMO(verb, weight, strategy, enc, pb, pbobjf, apmode, eps, searchStrat, redFact); // 
      _partition_parameter = partition_parameter;
    }
    
    void genLowerBoundSet();
    bool recycleLowerBoundSet();
    bool buildWorkFormula() override{return UnsatSatMO::buildWorkFormula();}
    void search_MO(int strategy) override{StratMO::search_MO(strategy);}
    bool extendUL(uint64_t * upperObjv, uint64_t * upperObix);
    virtual bool extendUL(YPoint& ul) override;
    void consolidateSolution() override;
    bool incorporate_approx() override;
    StatusCode partialSearch();

  protected:
    int _partition_parameter;
    void initializeOptimizer(Solver*, MaxSATFormula*) override;
    vec<Lit> explanation{}; 	// unsat explanation
  };
  void propagation(BLS::rootLits_t&, map<uint64_t,Lit>& vars, Solver* solver);
  void combination(BLS::rootLits_t&, BLS::rootLits_t&, map<uint64_t,Lit>& vars, Solver* solver);
}

#endif
#undef PARTIAL
#endif
