#ifndef ALG_HITTINGSETSMO_H
#define ALG_HITTINGSETSMO_H

#include "Alg_MasterMO.h"
#include "../Encoder.h"

#include "../MaxSAT.h"
#include "./Alg_BLS.h"
#include "./Alg_ServerMO.h"
#include "./Alg_UnsatSatMO.h"
#include "./Alg_UnsatSatMSU3MO.h"
#include "./Alg_UnsatSatStratMSU3MO.h"
#include "utils/System.h"
#include <utility>
#include <map>
#include <set>

#define MAXDIM 10

namespace openwbo {
    
  
  class HittingSetsMO : public virtual BLSMasterMO {
    class diagnosis: public std::pair<int, std::vector<Lit>>{
    public:
      int& id(){return first;}
      std::vector<Lit>& clause(){return second;}
};
  public:
    HittingSetsMO(int verb = _VERBOSITY_MINIMAL_, int weight = _WEIGHT_NONE_, int strategy = _WEIGHT_NONE_, 
	       int enc = _CARD_MTOTALIZER_, int pb = _PB_SWC_, 
	       int pbobjf = _PB_GTE_, int apmode = _ap_outvars_, float eps = 1, 
	       int searchStrat=3, float redFact=-1, bool objectiveCores = false) :
      BLS(verb, weight, strategy, enc, pb, pbobjf, apmode, eps, searchStrat, redFact), objectiveCores(objectiveCores)
{
  optim = new UnsatSatIncMO(verb, weight, strategy, enc, pb, pbobjf, apmode, eps, searchStrat, redFact);
  diagnoses = std::vector<diagnosis>{};
    }
    
    ~HittingSetsMO(){
      if (optim != NULL)
        delete optim;
      optim = NULL;
    }
    void vectorVec(const std::vector<Lit>& vector, vec<Lit>& vec);
    void genLowerBoundSet();
    bool buildWorkFormula() override;
    void incrementFormula();
    bool absorb(Solution::OneSolution& osol, int bvar);
    bool diagnose(Solution::OneSolution& osol, vec<Lit>&);
    bool virtual recycleLowerBoundSet();
    void initializeOptimizer(Solver* solv, MaxSATFormula* mxf) override;
    void consolidateSolution() override;
    bool setup_approx()override;
    bool incorporate_approx()override;
    void build() override;
    

  protected:
    std::vector<diagnosis> diagnoses;

    // CJ: option for extracting cores only over objective variables
    bool objectiveCores{};

  };
}

#endif
