// #define PARTIAL
#ifndef PARTIAL

#ifdef SIMP
#include "simp/SimpSolver.h"
#else
#include "core/Solver.h"
#endif

#include "../Encoder.h"

#include "../MaxSAT.h"
#include "./Alg_BLS.h"
#include "utils/System.h"
#include <utility>
#include <map>
#include <set>

#define MAXDIM 10

namespace openwbo {
    
  
  class PMinimalMO : public BLS {
    
  public:
    PMinimalMO(int verb = _VERBOSITY_MINIMAL_, int weight = _WEIGHT_NONE_, int strategy = _WEIGHT_NONE_, 
	       int enc = _CARD_MTOTALIZER_, int pb = _PB_SWC_, 
	       int pbobjf = _PB_GTE_, int apmode = _ap_outvars_, float eps = 1, 
	       int searchStrat=3, float redFact=-1) : 
      BLS(verb, weight, strategy, enc, pb, pbobjf, apmode, eps, searchStrat, redFact){
    }
    
    ~PMinimalMO(){}
    
    bool searchPMinimalMO();
    void search_MO(int strategy);
    void assumeDominatingRegion(uint64_t * objix, int nObj);
  protected:
    vec<Lit> explanation; 	// unsat explanation


  };
}

#endif
#undef PARTIAL
