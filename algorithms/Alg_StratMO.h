#ifndef ALG_STRATMO_H
#define ALG_STRATMO_H
#include "Alg_BLS.h"
#include "../partition.h"
#include "Alg_MasterMO.h"
#include <algorithm>

class StratMO: public virtual BLSMasterMO {
public:
  StratMO(int verb = _VERBOSITY_MINIMAL_, int weight = _WEIGHT_NONE_, int strategy = _WEIGHT_NONE_, 
	  int enc = _CARD_MTOTALIZER_, int pb = _PB_SWC_, 
	  int pbobjf = _PB_GTE_, int apmode = _ap_outvars_, float eps = 1, 
	  int searchStrat=3, int partition_parameter = 15, float redFact=-1): 
    BLS(verb, weight, strategy, enc, pb, pbobjf, apmode, eps, searchStrat, redFact),
    _part_par{partition_parameter}{}
  int part_par(){return _part_par;}
  // reascribes meaning to objRootLits
  virtual void incrementEncoding(partition::MyPartition::part_t& p);
  virtual void incrementEncoding(partition::MyPartition::part_t& p, int i);
  std::pair<rootLits::RootLits, BLS::invRootLits_t> partialEncoding(int iObj, PBObjFunction& pb);

  void build()override{
    BLS::build();
    partitions = generate();
  };
  vector<partition::MyPartition> generate(){
    int nObj = getFormula()->nObjFunctions();
    partitions.reserve(nObj);
    for(int i = 0; i < nObj; i++ ){
      partition::MyPartition p = {getFormula()->getObjFunction(i).get(), part_par()};
      p.reverse();
      partitions.push_back(p);
	
    }
    return partitions;
  }
  bool setup_approx() override{
    auto new_part = false;
    auto i = 0;
    for(auto& partition: partitions){
    pop:
      if(partition.size()){
	auto p = partition.pop();
	if(!p.size())
	  goto pop;
	printf("new part, size %lu, objective %d\n", p.size(), i);
	partition::logPart(p);
	incrementEncoding(p,i++);
	new_part = true;
      }
    }
    vector<bool> slices{};
    for(auto& el: objRootLits){
      auto slice = dynamic_cast<rootLits::RootLitsSliced*>(el.get());
      slices.push_back(slice->setSlice(true));
    }
    reverse(slices.begin(), slices.end());
    optim->copyObjRootLits(invObjRootLits, objRootLits, 0);
    for(auto& el: objRootLits){
      auto slice = dynamic_cast<rootLits::RootLitsSliced*>(el.get());
      slice->setSlice(*--slices.end());
      slices.pop_back();
    }

    optim->increment();
    return new_part;
  }
protected:
  vector<partition::MyPartition> partitions{};
  int _part_par = 0;
};
#endif
