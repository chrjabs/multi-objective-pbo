
// #define PARTIAL
#ifndef PARTIAL

#include <iostream>
#include "Alg_UnsatSatMSU3MO.h"	

using namespace openwbo;
//using namespace NSPACE;
using NSPACE::toLit;
using namespace std;

void UnsatSatMSU3IncMO::blockStep(const YPoint& yp){
  int var = blockSoft(solution().yPoint());
  block_assmpts.push(mkLit(var,false));
  assumptions.push(block_assmpts.last());
  blocking_vars[var]=var_type::soft;
  solution().note(-1) = var;
}


bool UnsatSatMSU3IncMO::searchUnsatSatMO() {
  if(!lowerBound.size())
    rootedSearch(YPoint{getFormula()->nObjFunctions()});
  else{
    YPoint dom = pareto::dominator(lowerBound);
    rootedSearch(dom);
  }
  if (solution().size() == 0) {
    return false;
  }else{
  }
  
  return true;
}

void UnsatSatMSU3IncObjMO::blockStep(const YPoint& yp){
  int var = blockSoft(solution().yPoint());
  block_assmpts.push(mkLit(var,false));
  assumptions.push(block_assmpts.last());
  blocking_vars[var]=var_type::soft;
  solution().note(-1) = var;
}


bool UnsatSatMSU3IncObjMO::searchUnsatSatMO() {
  if(!lowerBound.size())
    rootedSearch(YPoint{getFormula()->nObjFunctions()});
  else{
    YPoint dom = pareto::dominator(lowerBound);
    rootedSearch(dom);
  }
  if (solution().size() == 0) {
    return false;
  }else{
  }
  
  return true;
}

bool UnsatSatMSU3MO::searchUnsatSatMO() {
  int nObj = maxsat_formula->nObjFunctions();
  YPoint ul(nObj);
  for(int i = 0; i < nObj; i++){
    ul[i] = getTighterLB(i);
  }
  pushBlockedVars();
  auto dom = pareto::dominator(solution());
  if(pareto::dominates(ul,dom))
    ul = dom;
  rootedSearch(ul);
  if (solution().size() == 0) {
    return false;
  }else{
  return true;
  }
}

bool UnsatSatMSU3MO::extendUL(YPoint& yp){
  bool ret = false;
  ret = popBlockedVars();
  ret = ret || UnsatSatMO::extendUL(yp);
  return ret ;
}

bool UnsatSatMSU3MO::popBlockedVars(){
  bool extend = false;
  vec<Lit> conflict;
  Lit lit;
  solver->conflict.copyTo(conflict);
  int n = conflict.size();
  while(conflict.size() > 0){
    lit = conflict.last();
    conflict.pop();
    if(var(lit) < maxsat_formula->nInitialVars()){
      blockedVars.erase(lit);
      extend = true;
    }
  }
  if(extend)
    cout<<"c popped "<< n << " x vars\n";
  return extend;
}

void UnsatSatMSU3MO::pushBlockedVars(){
  for(int i = 0; i < getFormula()->nInitialVars(); i++)
    blockedVars.insert(mkLit(i));
}


void UnsatSatMSU3IncMO::checkSols() {
  vec<Lit> assmpts{maxsat_formula->nVars()};
  for(auto it = solution().begin(), end = solution().end();it != end;){
    int bvar = it->second.second;
    modelClause(Model{it->second.first.model()},assmpts);
    lbool sat = solver->solveLimited(assmpts);
    if(sat == l_False){
      lowerBound.push(std::move(it->second));
      it = solution().remove(it);
      blocking_vars.erase(bvar);
    }else
      ++it;
  }
}
void UnsatSatMSU3IncObjMO::checkSols() {
  lowerBound.clear();
  for(auto it = solution().begin(), end = solution().end();it != end;){
    lowerBound.push(it->second.first, it->second.second);
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
      blockStep(osol_n.yPoint());
    }else
      ++it;
  }

}
bool UnsatSatMSU3MO::buildWorkFormula(){
  // Init Structures
  return  updateMOFormulationIfSAT(); 
}


#endif
#undef PARTIAL

