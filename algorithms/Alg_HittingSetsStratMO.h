#ifndef ALG_HITTINGSETSSTRATMO_H
#define ALG_HITTINGSETSSTRATMO_H
#include "Alg_MasterMO.h"
#include "Alg_ServerMO.h"
#include "Alg_StratMO.h"
#ifdef SIMP
#include "simp/SimpSolver.h"
#else
#include "core/Solver.h"
#endif

#include "../Encoder.h"

#include "../MaxSAT.h"
#include "./Alg_BLS.h"
#include "./Alg_HittingSetsMO.h"
#include "./Alg_UnsatSatMO.h"
#include "./Alg_UnsatSatMSU3MO.h"
#include "./Alg_UnsatSatStratMSU3MO.h"
#include "utils/System.h"
#include <utility>
#include <map>
#include <set>

#define MAXDIM 10

namespace openwbo {
    
  
  class HittingSetsServerMO: public virtual HittingSetsMO, public virtual BLSServerMO{
  public:
    HittingSetsServerMO(int verb = _VERBOSITY_MINIMAL_, int weight = _WEIGHT_NONE_, int strategy = _WEIGHT_NONE_, 
			int enc = _CARD_MTOTALIZER_, int pb = _PB_SWC_, 
			int pbobjf = _PB_GTE_, int apmode = _ap_outvars_, float eps = 1, 
			int searchStrat=3, float redFact=-1)
      :BLS(verb, weight, strategy, enc, pb, pbobjf, apmode, eps, searchStrat, redFact),
       HittingSetsMO(verb,weight,strategy, enc,pb, pbobjf,apmode,eps, searchStrat,  redFact),
       BLSServerMO(verb,weight,strategy, enc,pb, pbobjf,apmode,eps, searchStrat,  redFact)
    {}
    void build() override;
    StatusCode searchAgain() override;
    void consolidateSolution() override;
    bool not_done() override;
    void increment() override;
    void checkSols() override;
    virtual void bootstrap(const Solution& sol) override;
    bool recycleLowerBoundSet() override;
  };

  class HittingSetsStratMO : public virtual StratMO{
    
  public:
    HittingSetsStratMO(int verb = _VERBOSITY_MINIMAL_, int weight = _WEIGHT_NONE_, int strategy = _WEIGHT_NONE_, 
		       int enc = _CARD_MTOTALIZER_, int pb = _PB_SWC_, 
		       int pbobjf = _PB_GTE_, int apmode = _ap_outvars_, float eps = 1, 
		       int searchStrat=3, int partition_parameter = 15, float redFact=-1) : 
      BLS(verb, weight, strategy, enc, pb, pbobjf, apmode, eps, searchStrat, redFact),
      StratMO(verb, weight, strategy, enc, pb, pbobjf, apmode, eps, searchStrat, partition_parameter, redFact)
    {
      optim_hs = new HittingSetsServerMO(verb, weight, strategy, enc, pb, pbobjf, apmode, eps, searchStrat, redFact);
      optim = optim_hs;
      diagnoses = std::vector<std::vector<Lit>>{};
      _partition_parameter = partition_parameter;
    

    }
    
    ~HittingSetsStratMO(){
      if (optim != NULL)
        delete optim;
      optim = NULL;
    }
    void vectorVec(const std::vector<Lit>& vector, vec<Lit>& vec);
    void genLowerBoundSet();
    void incrementFormula();
    void incrementEncoding(partition::MyPartition::part_t &p, int i) override;
    bool recycleLowerBoundSet();
    void initializeOptimizer(Solver* solv, MaxSATFormula* mxf) override;
    int part_par(){return _partition_parameter;}
    bool setupPSearch(std::vector<partition::MyPartition>& vps) ;
    StatusCode partialSearch(){return StatusCode::_OPTIMUM_;};
    void build() override;
    bool buildWorkFormula() override;
    bool incorporate_approx() override;
    bool setup_approx() override;
  protected:
    std::vector<std::vector<Lit>> diagnoses;
    int _partition_parameter;
    // set of variables from where cores can be extracted.
    std::vector<bool> x_current{};
  private:
    HittingSetsServerMO* optim_hs;
  };

}

#endif
