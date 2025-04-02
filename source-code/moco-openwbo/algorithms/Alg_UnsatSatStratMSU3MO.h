#ifndef UNSATSATSTRATMSU3MO
#define UNSATSATSTRATMSU3MO
// #define PARTIAL
#include "Alg_UnsatSatMO.h"
#ifndef PARTIAL

#ifdef SIMP
#include "simp/SimpSolver.h"
#else
#include "core/Solver.h"
#endif

#include "../Encoder.h"
#include "../Pareto.h"
#include "./Alg_UnsatSatMSU3MO.h"
#include "./Alg_ServerMO.h"
#include "utils/System.h"
#include <utility>
#include "../partition.h"
#include <algorithm>    // std::max

#define MAXDIM 10

namespace openwbo {
  using namespace partition;
  
  class UnsatSatStratMSU3MO : public UnsatSatMSU3MO{
    
  public:
    UnsatSatStratMSU3MO(int verb = _VERBOSITY_MINIMAL_, int weight = _WEIGHT_NONE_, int strategy = _WEIGHT_NONE_, 
			int enc = _CARD_MTOTALIZER_, int pb = _PB_SWC_, 
			int pbobjf = _PB_GTE_, int apmode = _ap_outvars_, float eps = 1, 
			int searchStrat=3, int partition_parameter = 15, float redFact=-1) : 
      BLS(verb, weight, strategy, enc, pb, pbobjf, apmode, eps, searchStrat, redFact), UnsatSatMSU3MO(verb, weight, strategy, enc, pb, pbobjf, apmode, eps, searchStrat, redFact){
      part_par = partition_parameter;
    }
    
    void blockStep(const YPoint& yp) override{
      for(auto& el: objRootLits){
	auto& sliced = *dynamic_cast<rootLits::RootLitsSliced*>(el.get());
	sliced.toggleSlice();
      }
      blockDominatedRegion(yp);
      for(auto& el: objRootLits){
	auto& sliced = *dynamic_cast<rootLits::RootLitsSliced*>(el.get());
	sliced.toggleSlice();
      }

    }


    std::vector<MyPartition> generate();
    bool searchUnsatSatMO() override;
    bool buildWorkFormula() override {
      bool result = UnsatSatMSU3MO::buildWorkFormula();
      // prepare slicing
      for(int i = 0; i < (int) objRootLits.size(); i++){
	auto old = *dynamic_cast<rootLits::RootLits*>(objRootLits[i].get());
	objRootLits[i] = 
	  std::make_unique<rootLits::RootLitsSliced>
	  (rootLits::RootLitsSliced{std::move(old), *getFormula()->getObjFunction(i)});
      }
      return result;
    }

  protected:
    int part_par;
    
  };
  class UnsatSatStratMSU3IncMO: public virtual UnsatSatStratMSU3MO, public virtual UnsatSatIncMO{
  public:
    UnsatSatStratMSU3IncMO(int verb = _VERBOSITY_MINIMAL_, int weight = _WEIGHT_NONE_, int strategy = _WEIGHT_NONE_, int enc = _CARD_MTOTALIZER_, 
			   int pb = _PB_SWC_, int pbobjf = _PB_GTE_, int apmode = _ap_outvars_, float eps = 1, 
			   int searchStrat=3, int partition_parameter = 15, float redFact=-1) : 
      BLS(verb, weight, strategy, enc, pb, pbobjf, apmode, eps, searchStrat, redFact), 
      UnsatSatStratMSU3MO(verb, weight, strategy, enc, pb, pbobjf, apmode, eps, searchStrat, partition_parameter, redFact), 
      UnsatSatIncMO(verb, weight, strategy, enc, pb, pbobjf, apmode, eps, searchStrat, redFact){}
      StatusCode searchAgain() override;
    // void checkSols() override;
    bool buildWorkFormula() override {return UnsatSatStratMSU3MO::buildWorkFormula();}
    void blockStep(const YPoint& yp) override{
      for(auto& el: objRootLits){
	auto& sliced = *dynamic_cast<rootLits::RootLitsSliced*>(el.get());
	sliced.toggleSlice();
      }
      UnsatSatIncMO::blockStep(yp);
      for(auto& el: objRootLits){
	auto& sliced = *dynamic_cast<rootLits::RootLitsSliced*>(el.get());
	sliced.toggleSlice();
      }
    }
    void build() override;
    void increment() override {return;}
    bool searchUnsatSatMO() override {return UnsatSatIncMO::searchUnsatSatMO();};
    bool rootedSearch(const YPoint& yp) override {return UnsatSatIncMO::rootedSearch(yp);}

    bool not_done() override{
      while(parts.size()){
	auto p = parts.pop();
	if(p.size()){
	  printf("new part, size %lu\n", p.size());
	  partition::logPart(p);
	  int i = 0;
	  for(auto& el: objRootLits){
	    PBObjFunction new_pb = dynamic_cast<rootLits::RootLitsSliced*>(el.get())->slice(p); 
	    getFormula()->replaceObjFunction(i++, make_unique<PBObjFunction>(new_pb));
	  }
	  return true;
	}
      }
      return false;
    }
  private:
    MyPartition parts;
  };
}

#endif
#undef PARTIAL
#endif
