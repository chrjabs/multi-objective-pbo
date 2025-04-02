#include <cstdint>
#include <iostream>
#include "Alg_UnsatSatMO.h"	
#include <algorithm>    // std::max
#include <memory>
#include "../Pareto.h"
#include "core/SolverTypes.h"
using namespace openwbo;
using namespace std;

vector<MyPartition> UnsatSatMO::generate(){
  int nObj = maxsat_formula->nObjFunctions();
  vector<MyPartition> partitions{};
  partitions.reserve(nObj);
  for(int i = 0; i < nObj; i++ )
    partitions.emplace_back(maxsat_formula->getObjFunction(i),30);
  return partitions;
}


MyPartition UnsatSatMO::mix(vector<MyPartition> partitions){
  MyPartition mix{};
  random_device rd{};
  mt19937 rng{rd()};
  uniform_int_distribution<int> uni(0,partitions.size()-1);
  int rnd{};
  while(partitions.size()){
    uniform_int_distribution<int> uni(0,partitions.size()-1);
    rnd = uni(rng);    
    mix.push(partitions[rnd].pop());

    if(partitions[rnd].size()==0){
      partitions[rnd]=partitions[partitions.size()-1];
      partitions.pop_back();
    }
  }
  return mix;   
}


bool UnsatSatMO::buildWorkFormula(){
  // Init Structures
  bool ret = updateMOFormulationIfSAT(); 
  return ret;
}

void UnsatSatMO::search_MO(int strategy){
  //Build solver
  build();
  //encode obj functions
  bool resform = buildWorkFormula();
        
  if(resform){
  
    printf("c search\n");
    switch(strategy){
    case 20: searchUnsatSatMO(); break;
    default:
      printf("Wrong strategy!\n");
      exit(EXIT_FAILURE);
    }
  }else{
    printf("c No more solutions!\n");
  }
  printf("c Done searching\n");
  BLS::consolidateSolution();
  if(solution().size() > 0){
    answerType = _OPTIMUM_;
  }else{
    answerType = _UNSATISFIABLE_;
  }
 
  printAnswer(answerType);
}

bool UnsatSatMO::rootedSearch(const YPoint& yp){
  double runtime = cpuTime();
  assumptions.clear();
  YPoint ul = yp;
  
 newHarvest:
  cout<<"c new harvest. upperLimit: "<< ul <<endl;
  assumptions.clear();
  //reinserts the MSU3 blocked vars
  
  for(const auto& el: blockedVars)
    assumptions.push(~el);
  assumeDominatingRegion(ul);

  while( solve() == l_True) {
    Model m =Model{solver->model};
    vec<lbool> proj;
    for(int i = 0, n = maxsat_formula->nInitialVars(); i < n; i++)
      proj.push(solver->model[i]);
    Model m1 = Model{proj};
    // Only block dominated region if m1 gets into the Solution
    if(solution().pushSafe(m1)){
      blockStep(solution().yPoint());
      printf("c o ");
      std::cout<< solution().yPoint()<<std::endl;
      runtime = cpuTime();
      printf("c new optimal solution (time: %.3f)\n", runtime - initialTime);
    }

  }
  if(extendUL(ul))
    goto newHarvest;

  return true;
}


//searches, starting at YPoint yp
bool UnsatSatIncMO::rootedSearch(const YPoint& yp){
  cout << "c rooted search\n";
  double runtime = cpuTime();
  YPoint ul = yp;
  //initialize blocking_vars
  blocking_vars.clear();
  for(const auto& el: solution())
    blocking_vars[el.second.second] = var_type::soft;

 newHarvest:
  cout<<"c new harvest. upperLimit: "<< ul <<endl;
  //rebuild the assumptions
  assumptions.clear();
  //block blocked vars. For instance, when doing a MSU3 search
  for(const auto& el: blockedVars)
    assumptions.push(~el);
  //block region dominated so far
  for(auto& el: blocking_vars)
    assumptions.push(Glucose::mkLit(el.first, false));
  assumeDominatingRegion(ul);
  while( solve() == l_True) {

    Model m =Model{solver->model};
    vec<lbool> proj;
    for(int i = 0, n = getFormula()->nInitialVars(); i < n; i++)
      proj.push(solver->model[i]);
    Model m1 = Model{proj};
    // Only block dominated region if m1 gets into the Solution

    if(solution().pushSafe(m1)){
      blockStep(solution().yPoint());
      // refresh blocking_vars, after pushSafe
      blocking_vars.clear();
      for(const auto& el: solution())
	blocking_vars[el.second.second] = var_type::soft;
      auto yp = solution().yPoint();

      runtime = cpuTime();
      cout <<"c o " << yp << endl;
      printf("c new inner optimal solution (time: %.3f)\n", runtime - initialTime);
    }
    else {
      auto one = Solution::OneSolution{&solution(),m1};
      cout <<"c o " << one.yPoint() << endl;
      runtime = cpuTime();
      printf("c new non-optimal solution (time: %.3f)\n", runtime - initialTime);
    }
    assumptions.clear();
    //block blocked vars. For instance, when doing a MSU3 search
    for(const auto& el: blockedVars)
      assumptions.push(~el);
    //block region dominated so far
    for(auto& el: blocking_vars)
      assumptions.push(Glucose::mkLit(el.first, false));
    assumeDominatingRegion(ul);
  }
  if(extendUL(ul))
    goto newHarvest;

  return true;
}
bool UnsatSatIncMO::searchUnsatSatMO() {
  int nObj = getFormula()->nObjFunctions();
  YPoint ul(nObj);
  for(int i = 0; i < nObj; i++){
    uint64_t lb = 0;
    // calc.setObjective(getFormula()->getObjFunction(i).get());
    // calc.buildObj();
    // auto lbl = calc.lowerBound();
    // if(lbl > lb) lb = lbl;
    ul[i] = lb;
  }
  if(!lowerBound.size())
    rootedSearch(ul);
  else{
    
    YPoint dom = pareto::dominator(lowerBound);
    if(pareto::dominates(ul,dom))
      ul = dom;
    rootedSearch(ul);
  }
  if (solution().size() == 0) {
    answerType=_UNSATISFIABLE_;
    return false;
  }else{
    answerType = _OPTIMUM_;
  }
  assumptions.clear();
  return true;
}

bool UnsatSatIncXsMO::searchUnsatSatMO() {
  int nObj = getFormula()->nObjFunctions();
  YPoint ul(nObj);
  if(!lowerBound.size())
    rootedSearch(ul);
  else{
    YPoint dom = pareto::dominator(lowerBound);
    if(pareto::dominates(ul,dom))
      ul = dom;
    rootedSearch(ul);
  }
  if (solution().size() == 0) {
    answerType=_UNSATISFIABLE_;
    return false;
  }else{
    answerType = _OPTIMUM_;
  }
  
  return true;
}
bool UnsatSatMO::searchUnsatSatMO() {
  int nObj = maxsat_formula->nObjFunctions();
  YPoint ul(nObj);
  // for(int i = 0; i < nObj; i++){
  //   ul[i] = getTighterLB(i);
  // }

  auto dom = pareto::dominator(solution());
  if(pareto::dominates(ul,dom))
    ul = dom;
  rootedSearch(ul);
  if (solution().size() == 0) {
    answerType=_UNSATISFIABLE_;
    return false;
  }else{
    answerType = _OPTIMUM_;
  }
  return true;
}


bool UnsatSatMO::extendUL(uint64_t * upperObjv, uint64_t * upperObjix){
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
      if(upperObjix[iObj] + 1 <= (*objRootLits[iObj]).size()){
	extend = true;
	upperObjix[iObj]++;
	upperObjv[iObj] = (*objRootLits[iObj])[upperObjix[iObj]].first;
      }
    }
  }

  return extend;
}
bool UnsatSatMO::extendUL(YPoint& yp){
  int nObj=yp.size();
  uint64_t upperObjv[nObj];
  uint64_t upperObjix[nObj];
  
  for(uint64_t i = 0; i < yp.size(); i++){
    upperObjv[i] = yp[i]+1;
  }
  evalToIndex(upperObjv, upperObjix);
  bool res =  extendUL(upperObjv, upperObjix);

  for(uint64_t i = 0; i < yp.size(); i++){
    yp[i] = upperObjv[i]-1;
  }
  return res;
}

void UnsatSatIncMO::blockStep(const YPoint& yp){
  int var = blockSoft(solution().yPoint());
  solution().note(-1) = var;
}
void UnsatSatIncXsMO::blockStep(const YPoint& yp){
  forceSlice(true);
  int var = blockSoft(solution().yPoint());
  solution().note(-1) = var;
  toggleSlice();
}

void UnsatSatIncXsMO::checkSols() {
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
      blockStep(osol_n.yPoint());
    }else
      ++it;
  }
}

void UnsatSatIncObjMO::checkSols() {
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
      blockStep(osol_n.yPoint());
    }else
      ++it;
  }
}

void UnsatSatIncMO::checkSols() {
  lowerBound.clear();
  vec<Lit> assmpts{maxsat_formula->nInitialVars()};
  for(auto it = solution().begin(), end = solution().end();it != end;){
    int bvar = it->second.second;
    if(marked_sols.count(it->first)){
      lowerBound.push(std::move(it->second));
      it = solution().remove(it);
      continue;
    }
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

void UnsatSatIncXsMO::incrementSlice(const partition::MyPartition::part_t& p){
  int i = 0;
  for(const auto& x: p)
    if(blockedVars.count(x))
      blockedVars.erase(x);
  
  for(auto& el: objRootLits){
    auto sliced = dynamic_cast<rootLits::RootLitsSliced*>(el.get());
    sliced->slice(p);
    // objective function is updated whenever the slice is incremented
    getFormula()->replaceObjFunction(i++, 
				     std::make_unique<PBObjFunction>(PBObjFunction{sliced->cur}));
  }

}
void UnsatSatIncXsMO::assumeDominatingRegion(const YPoint& yp) {
  forceSlice(true);
  UnsatSatMO::assumeDominatingRegion(yp);
  toggleSlice();
}

bool UnsatSatIncXsMO::extendUL(YPoint& ul){
  forceSlice(true);
  bool ret =  UnsatSatMO::extendUL(ul);
  toggleSlice();
  return ret;
}
const PBObjFunction& UnsatSatIncXsMO::slicedObjective(int i){
  auto sliced = dynamic_cast<rootLits::RootLitsSliced*>(objRootLits[i].get());
  return sliced->cur;
}
void UnsatSatIncXsMO::increment(){
  for(int i = 0, n = getFormula()->nObjFunctions(); i < n; i++){
    auto& lits = getFormula()->getObjFunction(i)->_lits;
    for(int j = 0, n = lits.size(); j < n; j++)
      blockedVars.insert(lits[j]);
  }
}
