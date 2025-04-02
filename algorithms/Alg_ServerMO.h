#ifndef INCMO
#define INCMO
// #define PARTIAL
#include "Alg_DynamicMO.h"
#include <cstdint>
#include <memory>
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
  
  class ServerMO : public DynamicMO{
  public: 
    virtual void setSolver(Solver* sol) =0;
    virtual Solver* getSolver() =0;
    virtual StatusCode searchAgain() = 0;
    virtual bool not_done(){return true;};
    virtual void bootstrap(const Solution& sol) = 0;
  };
  class BLSServerMO :public virtual BLS, public ServerMO{
  public:
    BLSServerMO(int verb = _VERBOSITY_MINIMAL_, int weight = _WEIGHT_NONE_, 
	     int strategy = _WEIGHT_NONE_, int enc = _CARD_MTOTALIZER_, int pb = _PB_SWC_,
	     int pbobjf = _PB_GTE_, int apmode = _ap_outvars_, float eps = 1, int searchStrat=3, float redFact=-1)
      :BLS(verb, weight, strategy, enc, pb, pbobjf, apmode, eps, searchStrat, redFact){}
    Solver* getSolver() override {return solver;}
    void setSolver(Solver* sol) override{solver=sol;}
    void setFormula(MOCOFormula& mff){
      mf = make_shared<MOCOFormula>(mff);
      ubCost = mf->getSumWeights();
    }
    //used to transfer a maxsatFormula into the MOCOFormula. Use it at
    //the infancy of the solver, only once.
    void loadFormula(MaxSATFormula* mxf) override{
      MOCO::loadFormula(mxf);
      // for(int i = 0; i < getFormula()->nObjFunctions(); i++)
      // 	objRootLits.push_back(std::make_unique<rootLits::RootLits>(rootLits::RootLits{}));
};
    void build() override{
      init();
      nbMCS = 0;
      answerType = _UNKNOWN_;
    };
    //BLs algorithms use order variables. This algorithm is also
    //incremental, so blocked regions may need to be unblocked again. 
    virtual int  blockSoft(const YPoint& yp){
      //controlling variable
      int var = solver->newVar();
      //computing the indexes given the objective value. This should be abstracted away...
      tmpBlockDominatedRegion(yp, mkLit(var,true));
      return var;
    }
    //sometimes the status of some blocked region gets updated. Either because
    //the region can be blocked permanently or because it should no longer be blocked at all.
    void hardenSoft(int var) {
      if(blocking_vars.count(var)){
	solver->addClause(mkLit(var,false));
	blocking_vars[var] = var_type::hard;
      }
    } 
    void disableSoft(int var) {
      if(blocking_vars.count(var)){
	solver->addClause(mkLit(var,true));
	blocking_vars.erase(var);
      }
    } 
    // outsorce building the rootLits
    void copyObjRootLits(invRootLits_t& inv, std::vector<rootLits_t>& dir, int n){
      invObjRootLits = inv;
      for(auto& entry:objRootLits)
	entry->clear();
      for(int i = 0; i < (int) dir.size(); i++)
	for(const auto& el: *dir[i])
	  objRootLits[i]->insert(el.first, el);
    }
    void checkLower() override{}
    //checking solutions after incrementing the partial objective functions
    void checkSols() override {
      for(auto it = solution().begin(), end = solution().end();it != end;){
	Solution::OneSolution osol = it->second.first;
	Model mod = osol.model();
	int bvar = it->second.second;
	auto osol_n = Solution::OneSolution{&solution(),mod};
	if(osol.yPoint() != osol_n.yPoint()){
	  it = solution().remove(it);
	  blocking_vars.erase(bvar);
	  // disabling permanently clause counterpart to bvar
	  solver->addClause(mkLit(bvar, true));
	  solution().push(osol_n.model());
	  blockSoft(osol_n.yPoint());
	}else
	  ++it;
      }
      marked_sols.clear();
    }
    void mark_solution(uint id){
      marked_sols.insert(id);
    }
  protected:
    // blocking literals for clauses that block regions dominated by
    // solutions that are still not known to be optimal
    vec<Lit> block_assmpts{};

    enum class var_type {
      soft,
      hard
    };
    //set of variables used to block regions softly
    std::map<int, var_type> blocking_vars{};
    std::set<uint64_t> marked_sols{};
    
  };
}

#endif
#undef PARTIAL
#endif
