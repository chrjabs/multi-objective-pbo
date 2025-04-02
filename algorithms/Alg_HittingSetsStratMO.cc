#include "Alg_HittingSetsStratMO.h"
#include "Alg_MasterMO.h"
#include "Alg_StratMO.h"
//after, optimizer is ready to work 
void HittingSetsStratMO::initializeOptimizer(Solver* solv, MaxSATFormula* mxf) {
  optim->setFormula(*getFormula());
  for(int i = 0,nObj = optim->getFormula()->nObjFunctions(); i < nObj; i++ ){
    optim->getFormula()->replaceObjFunction(i, std::make_unique<PBObjFunction>(PBObjFunction{}));
  }
  optim->setSolver(solv);
  optim->build();
  // objectives start out empty
  }

bool HittingSetsStratMO::incorporate_approx(){
  for(auto& el: optim->solution()){
    auto& osol = el.second.first;
    Solution::notes_t bvar = el.second.second;
    auto m = Model{osol.model()};
    solution().pushSafe(m, bvar, true, true);
    auto yp = solution().yPoint();
  }
  return false;
}
StatusCode HittingSetsServerMO::searchAgain(){
  return searchMasterMO();
}

bool HittingSetsServerMO::not_done(){return true;}
//will propagate changes innoculated by Master.
void HittingSetsServerMO::increment(){
  checkSols();
  optim->solution().clear();
  optim->lowerBound.clear();
  optim->copyObjRootLits(invObjRootLits, objRootLits, 0);
  for(int i = 0, n = getFormula()->nObjFunctions(); i < n; i++){
    auto& new_pb = *getFormula()->getObjFunction(i);
    optim->getFormula()->replaceObjFunction(i, make_unique<PBObjFunction>(new_pb));
  }
  for(auto& el: solution())
    optim->solution().pushSafe(el.second.first.model(), el.second.second);
  nbMCS = 0;
  answerType = _UNKNOWN_;
  optim->increment();
}
void HittingSetsServerMO::build(){
  init();
  nbMCS = 0;
  answerType = _UNKNOWN_;
  MaxSATFormula* f = new MaxSATFormula{};
  initializeOptimizer(solver, f);
}
void HittingSetsStratMO::build() {
  StratMO::build();
  initializeOptimizer(solver, NULL);
}
void HittingSetsServerMO::consolidateSolution(){
  recycleLowerBoundSet();  
}

void HittingSetsStratMO::incrementEncoding(partition::MyPartition::part_t &p, int i) {
  auto& el = objRootLits[i];
  auto slice = dynamic_cast<rootLits::RootLitsSliced*>(el.get());
  PBObjFunction new_pb = slice->slice(p); 
  optim->getFormula()->replaceObjFunction(i, make_unique<PBObjFunction>(new_pb));
  }
bool HittingSetsStratMO::buildWorkFormula(){
  Solver* tmp = solver;
  solver = optim_hs->optim->getSolver();
  auto res = StratMO::buildWorkFormula();
  solver = tmp;
	for(int i = 0; i < (int) objRootLits.size(); i++){
	  auto old = *dynamic_cast<rootLits::RootLits*>(objRootLits[i].get());
	  objRootLits[i] = 
	    std::make_unique<rootLits::RootLitsSliced>
	    (rootLits::RootLitsSliced{std::move(old), *getFormula()->getObjFunction(i)});
	}
  return res;
}
bool HittingSetsStratMO::setup_approx(){
  auto res = StratMO::setup_approx();
  return res;
}

void HittingSetsServerMO::checkSols(){
  for(auto it = solution().begin(), end = solution().end();it != end;){
    Solution::OneSolution osol = it->second.first;
    Model mod = osol.model();
    int bvar = it->second.second;
    auto osol_n = Solution::OneSolution{&solution(),mod};
    if(osol.yPoint() != osol_n.yPoint()){
      it = solution().remove(it);
      blocking_vars.erase(bvar);
      solution().push(osol_n.model());
    }else
      ++it;
  }
}
void HittingSetsServerMO::bootstrap(const Solution &usol){
  return;

}

bool HittingSetsServerMO::recycleLowerBoundSet(){
  const int nVars = getFormula()->nInitialVars();
  
  vec<Lit> assmpts{getFormula()->nInitialVars()};
  std::set<Lit> vars;
  for(int i = 0, n = getFormula()->nObjFunctions(); i < n; i++){
    auto& lits = getFormula()->getObjFunction(i)->_lits;
    for(int j = 0, m = lits.size() ;j < m; j++)
      vars.insert(lits[j]);
  }
  assmpts.clear();
  for(auto& el: vars)
    assmpts.push(el);
  if(!assmpts.size())
    return false;
  bool andf = true;
  for(auto& el: optim->solution()){
    auto& osol = el.second.first;
    int id = el.first;
    Solution::notes_t bvar = el.second.second;
    //checks satisfiability of complete model
    modelClause(modelEmbed(osol.model(), nVars),assmpts);
    lbool sat = solver->solveLimited(assmpts);

    if(sat == l_True)
      absorb(osol, bvar);
    else{
      andf = false;
      optim->mark_solution(id);
      diagnose(osol, assmpts);
    }
  }
  return andf;
}
