#ifndef UNSATSATMO
#define UNSATSATMO
// #define PARTIAL
#include "core/SolverTypes.h"
#include <memory>
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
#include "utils/System.h"
#include <utility>
#include "../partition.h"
#include <algorithm>    // std::max
#include <random>

#define MAXDIM 10

namespace openwbo {
  using namespace partition;
  
  class UnsatSatMO : public virtual BLS {
    
  public:
    UnsatSatMO(int verb = _VERBOSITY_MINIMAL_, int weight = _WEIGHT_NONE_, int strategy = _WEIGHT_NONE_, 
	       int enc = _CARD_MTOTALIZER_, int pb = _PB_SWC_, 
	       int pbobjf = _PB_GTE_, int apmode = _ap_outvars_, float eps = 1, 
	       int searchStrat=3, float redFact=-1) : 
      BLS(verb, weight, strategy, enc, pb, pbobjf, apmode, eps, searchStrat, redFact){
    }
    
    std::vector<MyPartition> generate();
    MyPartition mix(std::vector<MyPartition>);

    bool buildWorkFormula() override;
    StatusCode searchAgain();
    virtual bool searchUnsatSatMO();
    void search_MO(int strategy) override;
    bool extendUL(uint64_t * upperObjv, uint64_t * upperObix);
    virtual bool extendUL(YPoint& ul);
    const std::set<Lit>& blocked_vars(){return blockedVars;};
  protected:
    vec<Lit> explanation{}; 	// unsat explanation
    virtual bool rootedSearch(const YPoint& yp);
    std::set<Lit> blockedVars{};    
  };
  class UnsatSatIncMO: public virtual UnsatSatMO, public virtual BLSServerMO{
  public:
    UnsatSatIncMO(int verb = _VERBOSITY_MINIMAL_, int weight = _WEIGHT_NONE_, int strategy = _WEIGHT_NONE_, 
	       int enc = _CARD_MTOTALIZER_, int pb = _PB_SWC_, 
	       int pbobjf = _PB_GTE_, int apmode = _ap_outvars_, float eps = 1, 
	       int searchStrat=3, float redFact=-1) : 
      BLS(verb, weight, strategy, enc, pb, pbobjf, apmode, eps, searchStrat, redFact), 
      UnsatSatMO(verb, weight, strategy, enc, pb, pbobjf, apmode, eps, searchStrat, redFact), 
      BLSServerMO(verb, weight, strategy, enc, pb, pbobjf, apmode, eps, searchStrat, redFact){}
    bool buildWorkFormula() override {return updateMOFormulationIfSAT();}
    Solver* getSolver() override {return solver;}
    StatusCode searchAgain() override {
      setInitialTime(cpuTime());
      for(auto& el: solution())
	if(el.second.second < 0){
	  auto bvar = blockSoft(el.second.first.yPoint());
	  el.second.second = bvar;
	}
      searchUnsatSatMO();
      return answerType;
    }
    void checkSols() override;
    void increment() override {return;}
    bool searchUnsatSatMO() override;
    bool rootedSearch(const YPoint& yp) override;
    void blockStep(const YPoint& yp) override;
    void build() override{
      BLSServerMO::build();
      solver = buildSolverMO();
    }
    virtual void bootstrap(const Solution& sol) override {return;}
  };

  class UnsatSatIncXsMO: public UnsatSatIncMO{
  public:
    UnsatSatIncXsMO(int verb = _VERBOSITY_MINIMAL_, int weight = _WEIGHT_NONE_, int strategy = _WEIGHT_NONE_, 
		     int enc = _CARD_MTOTALIZER_, int pb = _PB_SWC_, 
		     int pbobjf = _PB_GTE_, int apmode = _ap_outvars_, float eps = 1, 
		     int searchStrat=3, float redFact=-1) : 
      BLS(verb, weight, strategy, enc, pb, pbobjf, apmode, eps, searchStrat, redFact), 
      UnsatSatMO(verb, weight, strategy, enc, pb, pbobjf, apmode, eps, searchStrat, redFact){}

    void build() override {
      BLSServerMO::build();
      solver = buildSolverMO();
      for(int i = 0; i < getFormula()->nInitialVars(); i++)
	blockedVars.insert({Glucose::mkLit(i,false)});
    }

    //RootLits are ready to be sliced
    bool buildWorkFormula() override  {
      bool result = UnsatSatIncMO::buildWorkFormula();
      // prepare slicing
      for(int i = 0; i < (int) objRootLits.size(); i++){
	auto old = *dynamic_cast<rootLits::RootLits*>(objRootLits[i].get());
	objRootLits[i] = 
	  std::make_unique<rootLits::RootLitsSliced>
	  (rootLits::RootLitsSliced{std::move(old), *getFormula()->getObjFunction(i)});
	// the complete objective functions were stored into
	// objRootLits. Their place is taken by the partial functions
	// from now on
	getFormula()->replaceObjFunction(i, 
					 std::make_unique<PBObjFunction>(PBObjFunction{}));
      }
      return result;
    }
    bool not_done() override {return true;}
    bool searchUnsatSatMO() override;
    // toggles sliced state of objRootLits. Requires that objRootLits
    //is a objRootLitsSliced subtype
    void toggleSlice(){
      for(auto& el: objRootLits){
	auto& sliced = *dynamic_cast<rootLits::RootLitsSliced*>(el.get());
	sliced.toggleSlice();
      } 
    }
    void forceSlice(bool b){
      for(auto& el: objRootLits){
	auto& sliced = *dynamic_cast<rootLits::RootLitsSliced*>(el.get());
	sliced.setSlice(b);
      } 
    }
    // updates the slice
    void incrementSlice(const partition::MyPartition::part_t&);
    void assumeDominatingRegion(const YPoint& yp) override;
    bool extendUL(YPoint& ul) override;
    void checkSols() override;
    void blockStep(const YPoint& yp) override;
    const PBObjFunction& slicedObjective(int i);
    void increment() override;
  private:
    std::vector<std::unique_ptr<PBObjFunction>> off{}; 
  };

    class UnsatSatIncObjMO: public UnsatSatIncMO{
  public:
    UnsatSatIncObjMO(int verb = _VERBOSITY_MINIMAL_, int weight = _WEIGHT_NONE_, int strategy = _WEIGHT_NONE_, 
		     int enc = _CARD_MTOTALIZER_, int pb = _PB_SWC_, 
		     int pbobjf = _PB_GTE_, int apmode = _ap_outvars_, float eps = 1, 
		     int searchStrat=3, float redFact=-1) : 
      BLS(verb, weight, strategy, enc, pb, pbobjf, apmode, eps, searchStrat, redFact), 
      UnsatSatMO(verb, weight, strategy, enc, pb, pbobjf, apmode, eps, searchStrat, redFact){}
    bool buildWorkFormula() override {return updateMOFormulationIfSAT();}
    void checkSols() override;
    void increment() override {return;}
  private:
    std::vector<std::unique_ptr<PBObjFunction>> off{}; 


  };

}

#endif
#undef PARTIAL
#endif
