
// #define PARTIAL
#ifndef PARTIAL

#include <iostream>
#include "Alg_UnsatSatStratMSU3MO.h"	
#include <algorithm>    // std::max
#include "../Pareto.h"
using namespace openwbo;
using namespace std;
using NSPACE::toLit;

vector<MyPartition> UnsatSatStratMSU3MO::generate(){
  int nObj = getFormula()->nObjFunctions();
  vector<MyPartition> partitions{};
  partitions.reserve(nObj);
  for(int i = 0; i < nObj; i++ )
    partitions.emplace_back(getFormula()->getObjFunction(i).get(), part_par);
  return partitions;
}

bool UnsatSatStratMSU3MO::searchUnsatSatMO() {
  int nObj = getFormula()->nObjFunctions();
  YPoint ul(nObj);
  for(int i = 0; i < nObj; i++){
    ul[i] = getTighterLB(i);
  }
  MyPartition partition = partition::mix(generate());
  while(partition.size())
    {
      auto p = partition.pop();
      if(p.size()){
	partition::logPart(p);
	printf("new part, size %lu\n", p.size());
	for(auto& el: objRootLits){
	  PBObjFunction new_pb = dynamic_cast<rootLits::RootLitsSliced*>(el.get())->slice(p); 
	}
	for(auto& el: p)
	  blockedVars.insert(el);
	auto dom = pareto::dominator(solution());
	if(pareto::dominates(ul,dom))
	  ul = dom;
	rootedSearch(ul);
      }
      else
	continue;

    }

 if (solution().size() == 0) {
    answerType=_UNSATISFIABLE_;
    return false;
  }else{
    answerType = _OPTIMUM_;
  }
  return true;
}

StatusCode UnsatSatStratMSU3IncMO::searchAgain(){
  setInitialTime(cpuTime()); 
  searchUnsatSatMO(); 
  return answerType;
}

void UnsatSatStratMSU3IncMO::build(){
  BLS::build();
  parts = partition::mix(generate());
}

#endif
#undef PARTIAL

