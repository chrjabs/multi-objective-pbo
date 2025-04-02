
// #define PARTIAL
#include <vector>
#ifndef PARTIAL

#include <iostream>
#include "Alg_UnsatSatStratMO.h"	
#include <algorithm>    // std::max
#include "../Pareto.h"
using namespace openwbo;
using namespace std;
using NSPACE::toLit;

void UnsatSatStratMO::initializeOptimizer(Solver* solv, MaxSATFormula* mxf){
  optim->loadFormula(mxf);
  optim->build();
  optim->setSolver(solv);
  auto nObj = optim->getFormula()->nObjFunctions();
  // objectives start out as empty
  for(int i = 0; i < nObj; i++ ){
    optim->getFormula()->replaceObjFunction(i, std::make_unique<PBObjFunction>(PBObjFunction{}));
  }
  // build remanining things, for instance, inintialize rootLits
  optim->buildWorkFormula();
}

bool UnsatSatStratMO::incorporate_approx(){
  auto ret = recycleLowerBoundSet();
  BLS::consolidateSolution();
  return ret;
}

void UnsatSatStratMO::genLowerBoundSet(){
  optim->searchAgain();
}

bool UnsatSatStratMO::recycleLowerBoundSet(){
  vec<Lit> assmpts{maxsat_formula->nVars()};
  for(auto& el: optim->solution()){
    auto& osol = el.second.first;
    Solution::notes_t bvar = el.second.second;
    auto m = Model{osol.model()};
    solution().pushSafe(m, bvar, true, true);
    auto yp = solution().yPoint();
  }
  return false;
}



bool UnsatSatStratMO::extendUL(uint64_t * upperObjv, uint64_t * upperObjix){
  bool extend = false;
  vec<Lit> conflict;
  Lit lit;
  int iObj;
  solver->conflict.copyTo(conflict);
  while(conflict.size() > 0){
    lit = conflict.last();
    conflict.pop();
    iObj = getIObjFromLit(lit);
    if(iObj > -1){
      if(upperObjix[iObj] + 1 < (*objRootLits[iObj]).size()){
	extend = true;
	upperObjix[iObj]++;
	upperObjv[iObj] = (*objRootLits[iObj])[upperObjix[iObj]].first;
      }
    }
  }

  return extend;
}
bool UnsatSatStratMO::extendUL(YPoint& yp){
  int nObj=yp.size();
  uint64_t upperObjv[nObj];
  uint64_t upperObjix[nObj];
  
  for(uint64_t i = 0; i < yp.size(); i++){
    upperObjv[i] = yp[i];
  }
  evalToIndex(upperObjv, upperObjix);
  bool res =  extendUL(upperObjv, upperObjix);

  for(uint64_t i = 0; i < yp.size(); i++){
    yp[i] = upperObjv[i] ;
  }
  return res;
}



void UnsatSatStratMO::consolidateSolution(){
  recycleLowerBoundSet();
  BLS::consolidateSolution();
}



#endif
#undef PARTIAL

