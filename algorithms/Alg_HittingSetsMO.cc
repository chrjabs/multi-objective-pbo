
// #define PARTIAL
#include <utility>
#ifndef PARTIAL

#include <iostream>
#include "Alg_HittingSetsMO.h"	
#include <algorithm>    // std::max

using namespace openwbo;
//using namespace NSPACE;
using NSPACE::toLit;

//this clones every variable in the solver. Make sure it is called
//accordingly
void HittingSetsMO::initializeOptimizer(Solver* solv, MaxSATFormula* mxf){
  auto& f = *mxf;
  for(int i = 0, n = getFormula()->nObjFunctions(); i < n; ++i){
    auto pb = getFormula()->getObjFunction(i).get();
    f.addObjFunction(pb);
  }
  for(int i = 0, n = getFormula()->nInitialVars(); i < n; ++i)
    f.newVar();
  {
    f.setInitialVars(f.nVars());
    optim->loadFormula(&f);
    optim->build();
    auto formula = optim->getFormula();
    int64_t min = 0, max = 0;
    for (int i = 0; i < formula->nObjFunctions(); i++){
      max = formula->getObjFunction(i)->getUB();
      formula->setUB(i, max);
      formula->setTighterUB(i, max);
      min = formula->getObjFunction(i)->getLB();
      formula->setLB(i, min);
      formula->setTighterLB(i, min);
    };
    formula->setFormat(_FORMAT_PB_);
  }
}

void HittingSetsMO::genLowerBoundSet(){
  optim->searchAgain();
}

// CJ: save core and extract additional disjoint cores (if they exist)
bool HittingSetsMO::diagnose(Solution::OneSolution& osol, vec<Lit>& assmpts){
  auto sat =l_False;
  do{
    diagnoses.emplace_back();
    auto& conflict = (--diagnoses.end())->clause();
    for(int i = 0, n = solver->conflict.size(); i<n;i++){
      //i.e, diagnoses clause satisfied iff conflict is hit
      Lit lit = solver->conflict[i];
      conflict.push_back(lit);
      assmpts.remove(~lit);
    }
    sat = solver->solveLimited(assmpts);
  }while(sat == l_False);
  return true;
}

bool HittingSetsMO::absorb(Solution::OneSolution& osol, int bvar){
  {
      auto m = Model{osol.model()};
      //removes elements of solution that are dominated by m.
      if(solution().pushSafe(m, bvar, true, true)){
	auto yp = solution().yPoint();
	cout <<"c o " << yp << endl;
	auto runtime = cpuTime();
	printf("c new optimal solution (time: %.3f)\n", runtime - initialTime);
      }
    }
 return true;
}

//saves sat lower bound set
// CJ: extract cores based on found non-dominated set
bool HittingSetsMO::recycleLowerBoundSet(){
  const int nVars = getFormula()->nInitialVars();
  
  vec<Lit> assmpts{getFormula()->nInitialVars()};
  vec<Lit> partial_assmpts{};
  bool andf = true;
  for(auto& el: optim->solution()){
    auto& osol = el.second.first;
    int id = el.first;
    Solution::notes_t bvar = el.second.second;
    // CJ: generate assumptions
    if (objectiveCores) {
      // CJ: assumptions only include objective variables
      objAssumptions(osol.model(), getFormula()->maxsat_formula(), assmpts);
    } else {
      //checks satisfiability of complete model
      modelClause(modelEmbed(osol.model(), nVars),assmpts);
    }
    lbool sat = solver->solveLimited(assmpts);

    if(sat == l_True)
      // TODO CJ: absorb model of newest SAT call rather than solution, as the
      // model from the solution might actually be UNSAT
      absorb(osol, bvar);
    else{
      andf = false;
      optim->mark_solution(id);
      // CJ: save core and extract additional disjoint cores (if they exist)
      diagnose(osol, assmpts);
    }
  }
  return andf;
}


void HittingSetsMO::incrementFormula(){
  cout << "diagnoses size: " << diagnoses.size() << endl;
  for(auto& diag: diagnoses){
    vec<Lit> vecDiag(diag.clause().size());
    vectorVec(diag.clause(), vecDiag);
    optim->getSolver()->addClause(vecDiag);
  }
  optim->checkSols();
}
bool HittingSetsMO::setup_approx(){
  if(!diagnoses.size())
    return false; // CJ: terminates IHS algorithm
  incrementFormula();
  diagnoses.clear();
  optim->increment(); // CJ: does nothing
  return optim->not_done(); // CJ: returns true -> new iteration
}

bool HittingSetsMO::incorporate_approx(){
  consolidateSolution();
  return true;
}

void HittingSetsMO::consolidateSolution(){
  recycleLowerBoundSet();  
  BLS::consolidateSolution();
}

void HittingSetsMO::vectorVec(const std::vector<Lit>& vector, vec<Lit>& vec){
  for(int i = 0, n = vec.size();  i < n; i++)
    vec[i] = vector[i];
}
void HittingSetsMO::build(){
  BLS::build();
  MaxSATFormula* f = new MaxSATFormula{};
  initializeOptimizer(NULL, f);
}
bool HittingSetsMO::buildWorkFormula(){
  return optim->buildWorkFormula();
}

#endif
#undef PARTIAL

