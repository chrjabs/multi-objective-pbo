#include "Alg_BLS.h"
#include <algorithm>    // std::max
#include <memory>

using namespace openwbo;
//using namespace NSPACE;
using NSPACE::toLit;

// #define __DEBUG__


/************************************************************************************************
 //
 // Utils for model management
 //
 ************************************************************************************************/ 

/*_________________________________________________________________________________________________
  |
  |  saveModel : (currentModel : vec<lbool>&)  ->  [void]
  |  
  |  Description:
  |		 
  |    Saves the current model found by the SAT solver.
  |
  |  Pre-conditions:
  |    * Assumes that 'nbInitialVariables' has been initialized.
  |    * Assumes that 'currentModel' is not empty.
  |
  |  Post-conditions:
  |    * 'model' is updated to the current model.
  |    * 'nbSatisfiable' is increased by 1.
  |
  |________________________________________________________________________________________________@*/
void BLS::saveModel(vec<lbool> &currentModel){
  //assert (n_initial_vars != 0);
  assert (currentModel.size() != 0);
  
  model.clear();
  // Only store the value of the variables that belong to the original MaxSAT formula.
  for (int i = 0; i < getFormula()->nVars(); i++){
    model.push(currentModel[i]);
  }
  
  nbSatisfiable++;
}




void BLS::saveSmallestModel(vec<lbool> &currentModel){
  //assert (n_initial_vars != 0);
  assert (currentModel.size() != 0);
  
  _smallestModel.clear();
  // Only store the value of the variables that belong to the original MaxSAT formula.
  for (int i = 0; i < getFormula()->nVars(); i++){
    _smallestModel.push(currentModel[i]);
  }
}


/************************************************************************************************
 //
 // BLS search
 //
 ************************************************************************************************/ 

// Disjoint cores

void BLS::identifyDisjointCores() {
  lbool res = l_True;
  std::map<Lit,int> assumptionMapping;
  uint64_t w = 0;
  
  assumptions.clear();
  _cores.clear();
  _coreSatClauses.clear();
  _coreUnsatClauses.clear();
  _satClauses.clear();
  _prevAssumptions.clear();
  
  for (int i = 0; i < getFormula()->nSoft(); i++) {
    assumptions.push(~(getFormula()->getSoftClause(i).assumption_var));
    assumptionMapping[getFormula()->getSoftClause(i).assumption_var] = i;
  }
  
  do {
    res = solve();
    if (res == l_False) {
      printf("c Core #%d found.\n", _cores.size());
      // Core found!
      // If core is empty, then hard clause set is unsat!
      if (solver->conflict.size() == 0) {
	if (nbMCS == 0) {
	  printAnswer(_UNSATISFIABLE_);
	  exit(_UNSATISFIABLE_);
	}
	else {
	  // It is not the first MCS. Hence, all MCS were found.
	  printAnswer(_OPTIMUM_);
	  exit(_OPTIMUM_);
	}
      }
      
      _cores.push();
      _prevAssumptions.push();
      _coreSatClauses.push();
      _coreUnsatClauses.push();
      
      for (int i = 0; i < solver->conflict.size(); i++) {
	int indexSoft = coreMapping[solver->conflict[i]];
	int assumptionPos = assumptionMapping[solver->conflict[i]];
	
	//Add to core
	_cores[_cores.size()-1].push(indexSoft);
	
	//remove from assumptions
	assumptions[assumptionPos] = assumptions.last();
	assumptionMapping[~(assumptions.last())] = assumptionPos;
	assumptions.pop();
      }
      
      w += coreMinCost(_cores.size()-1);
    }
    else if (res == l_True) {
      for (int i = 0; i < assumptions.size(); i++) {
	_satClauses.push(coreMapping[~(assumptions[i])]);
      }
    }
  } while (res == l_False);
  printf("c # Disjoint Cores: %d\n", _cores.size());
  if (w > _lbWeight) _lbWeight = w;
}




uint64_t BLS::coreMinCost(int c) {
  unsigned int w = UINT_MAX;
  for (int i = 0; i < _cores[c].size(); i++)
    if (getFormula()->getSoftClause(_cores[c][i]).weight < w)
      w = getFormula()->getSoftClause(_cores[c][i]).weight;
  return w;
}


// Checks if a soft clause is satisfied by saved model

bool BLS::satisfiedSoft(int i) {
  for (int j = 0; j < getFormula()->getSoftClause(i).clause.size(); j++){
    assert (var(getFormula()->getSoftClause(i).clause[j]) < model.size());
    if ((sign(getFormula()->getSoftClause(i).clause[j]) && model[var(getFormula()->getSoftClause(i).clause[j])] == l_False) || 
	(!sign(getFormula()->getSoftClause(i).clause[j]) && model[var(getFormula()->getSoftClause(i).clause[j])] == l_True)) {
      return true;
    }
  }
  return false;
}





// Call to the SAT solver...

lbool BLS::solve() {
  nbSatCalls++;  
   lbool res;
#ifdef SIMP
  res = ((SimpSolver*)solver)->solveLimited(assumptions);
#else
  res =  solver->solveLimited(assumptions);
#endif
  if(res == l_True)
    nbSatisfiable++;
  return res;
}


void BLS::addMCSClause(vec<int>& unsatClauses) {
    
  if(getFormula()->nObjFunctions() == 1){
    vec<Lit> lits;
    for (int i = 0; i < unsatClauses.size(); i++) {
      lits.push(~getFormula()->getSoftClause(unsatClauses[i]).assumption_var);
    }

    //addHardClause(lits);
    getFormula()->addHardClause(lits);
  
    solver->addClause(getFormula()->getHardClause(getFormula()->nHard()-1).clause);
    //solver->addClause(lits);
  }
}




void BLS::initUndefClauses(vec<int>& undefClauses) {
  for (int i = 0; i < getFormula()->nSoft(); i++) 
    undefClauses.push(i);
}


void BLS::basicCoreBasedSearch(int maxMCS = 10, bool maxSAT = false) {
  // Init Structures
  init();
  
  //Build solver
  solver = buildSolver();
  
  identifyDisjointCores();
  
  while (_nMCS < maxMCS) {
    bool foundMCS = findNextCoreBasedMCS();
    
    if (foundMCS) _nMCS++;
    // if (!foundMCS || (maxSAT && _nMCS % 100 == 0))
    //   identifyDisjointCores();
    if (!foundMCS) identifyDisjointCores();
    if (maxSAT && _smallestMCS <= _lbWeight) {
      printAnswer(_OPTIMUM_);
      exit(_OPTIMUM_);      
    }
  }
  printf("c All requested MCSs found\n");
}



void BLS::addBackboneLiterals(int index) {
  assumptions.push(getFormula()->getSoftClause(index).assumption_var);
  for (int j = 0; j < getFormula()->getSoftClause(index).clause.size(); j++){
    assumptions.push(~(getFormula()->getSoftClause(index).clause[j]));
  }
}



// Find next MCS.

bool BLS::findNextCoreBasedMCS() {
  vec<int> undefClauses;
  uint64_t costModel = 0;
  lbool res = l_True;
  int initCore = 0;
  
  // Initialize CostModel. This could be done incrementally!!!
  while (initCore < _cores.size() && _coreUnsatClauses[initCore].size() > 0) {
    for (int i = 0; i < _coreUnsatClauses[initCore].size(); i++) 
      costModel += getFormula()->getSoftClause(_coreUnsatClauses[initCore][i]).weight;
    initCore++;
  }
  
  do {
    res = solve();
    if (res == l_False) {
      // All MCS with these disjoint cores were enumerated. Must get new set of disjoint cores.
      if (initCore == 0) return false;
      
      initCore--;
      while (assumptions.size() > _prevAssumptions[initCore])
	assumptions.pop();
      
      _coreSatClauses[initCore].clear();
      
      for (int i = 0; i < _coreUnsatClauses[initCore].size(); i++)
	costModel -= getFormula()->getSoftClause(_coreUnsatClauses[initCore][i]).weight;
      
      _coreUnsatClauses[initCore].clear();
    }
  } while (res == l_False);

  for (int c = initCore; c < _cores.size(); c++) {
    for (int i = 0; i < _cores[c].size(); i++) 
      costModel += getFormula()->getSoftClause(_cores[c][i]).weight;
  }
  
  // Iterate all cores to find next MCS
  for (int c = initCore; c < _cores.size(); c++) {
    undefClauses.clear();
    
    for (int i = 0; i < _cores[c].size(); i++) undefClauses.push(_cores[c][i]);
    
    _prevAssumptions[c] = assumptions.size();
    
    while (undefClauses.size()) {
      int i = undefClauses.last();
      undefClauses.pop();
      _coreSatClauses[c].push(i);
      assumptions.push(~(getFormula()->getSoftClause(i).assumption_var));
      
      res = solve();
      
      if (res == l_False) {
	_coreUnsatClauses[c].push(_coreSatClauses[c].last());
	_coreSatClauses[c].pop();
	assumptions.pop();
	addBackboneLiterals(i);
      }
      else if (res == l_True) {
	costModel -= getFormula()->getSoftClause(_coreSatClauses[c].last()).weight;
	saveModel(solver->model);
	
	i = 0;
	//Remove satisfied soft clauses from undefClauses
	while (i < undefClauses.size()) {
	  if (satisfiedSoft(undefClauses[i])) {
	    //Add soft clause as Hard!
	    _coreSatClauses[c].push(undefClauses[i]);
	    assumptions.push(~(getFormula()->getSoftClause(undefClauses[i]).assumption_var));
	    costModel -= getFormula()->getSoftClause(undefClauses[i]).weight;
	    undefClauses[i] = undefClauses.last();
	    undefClauses.pop();
	  }
	  else i++;
	}
	if (costModel < _smallestMCS) {
	  saveSmallestModel(solver->model);
	  printf("o %ld\n", costModel);
	  fflush(stdout);
	}
      }
      else {
	printf("c SAT Solver unable to solve formula?!?!?\n");
      }
    }
  }
  
  if (costModel < _smallestMCS) {
    _smallestMCS = costModel;
    //Last saved model is smallest MCS...
  }
  nbMCS++;
  printf("c MCS #%d Weight: %ld\n", nbMCS, costModel);
  
  vec<int> unsatClauses;
  for (int c = 0; c < _cores.size(); c++) {
    for (int i = 0; i < _coreUnsatClauses[c].size(); i++)
      unsatClauses.push(_coreUnsatClauses[c][i]);
  }
  addMCSClause(unsatClauses);
  
  return true;
}




void BLS::search_MO(int strategy){
  // Init Structures
  init();
    
  //     printf("\t\tc strategic_search\n");
  //     printf("c eps: %f\n", epsilon);
    
  //Build solver
  double epsthreshold = 1 + 1e-4;
    
  solver = buildSolverMO();    
    
  bool resform, terminate=false;
  nbMCS = 0;
    
  answerType = _UNKNOWN_;
    
  bool permanentBlock = false;
    
    
  while(!terminate){
    //encode obj functions
        
    resform = updateMOFormulationIfSAT();
        
    if(resform){
  
      printf("c search\n");
      switch(strategy){
      case 20: MCSSearchMO(); break;
      case 21: permanentBlock = true;
	assert(approxMode == 0);
	MCSSearchMO(permanentBlock); break;
                
      default:
	printf("Wrong strategy!\n");
	exit(EXIT_FAILURE);
      }
    }else{
      printf("c No more solutions!\n");
    }
    printf("c Done searching\n");
    printf("c epsilon: %f\n", epsilon);
    printf("c reductionFactor: %f\n", redFactor);
    if(epsilon <= 1 
       || redFactor <= 0 
       || (!resform &&
	   (permanentBlock 
	    || (approxMode == _ap_coeffs_ && apObjectivesAreExact())))){
      //         if(epsilon <= 1 || redFactor < 0){
      terminate = true;
      printf("c time to terminate\n");
    }else{
      updateExpectedApprox(epsilon, nreencodes);
      if(epsilon <= epsthreshold)
	epsilon = 1;
      else
	epsilon = 1 + (epsilon-1)/redFactor;
      printf("c REENCODE epsilon = %f\n", epsilon);
    }
  }
    
  if(nondom.size() > 0){
        
    if(epsilon <= 1 || (permanentBlock && redFactor > 0)){
      printf("c LBset = PF\n");
      clearLowerBoundSet();
      for(size_t i = 0; i < nondom.size(); i++)
	updateLowerBoundSet(nondom[i], false, nreencodes);
    }
    if(redFactor > 0)
      updateExpectedApprox(1, nreencodes);
    else
      updateExpectedApprox(epsilon, nreencodes);

    answerType = _OPTIMUM_;
  }else{
    if(nreencodes == 0){
      answerType = _UNSATISFIABLE_;
    }else if(nreencodes == 1){
      clearLowerBoundSet();
    }
  }
 
  printAnswer(answerType);
}



void BLS::tmpBlockDominatedRegion(uint64_t * objix, int nObj, Lit tmplit){
  vec<Lit> d_cl;
  d_cl.push(tmplit);
  for(int di = 0; di < nObj; di++){
    int j = objix[di];
    if(j > 0){
      d_cl.push((*objRootLits[di])[j].second);
      //                 printf("block z%d\n", var(objRootLits[di][j].second)+1);
#ifdef __DEBUG__
      printf("c tmp-block %d [%lu][var: %d]\n", di, objRootLits[di][j].first, var(objRootLits[di][j].second));
#endif
    }
  }
  solver->addClause(d_cl);
}
void BLS::tmpBlockDominatedRegion(const YPoint& yp, Lit tmplit){
  int nObj = yp.size();
  //computing the indexes given the objective value. This should be abstracted away...
  uint64_t objix[nObj];
  evalToIndex(yp, objix);
  for(int i = 0; i < nObj; i++)
    objix[i]--;
  tmpBlockDominatedRegion(objix, nObj, tmplit);
}



void BLS::blockDominatedRegion(uint64_t * objix, int nObj){
  vec<Lit> d_cl;
  for(int di = 0; di < nObj; di++){
    int j = objix[di];
    if(j > 0){
      d_cl.push((*objRootLits[di])[j].second);
      //                 printf("block z%d\n", var(objRootLits[di][j].second)+1);
#ifdef __DEBUG__
      printf("c block %d [%lu][var: %d]\n", di, (*objRootLits[di])[j].first, var((*objRootLits[di])[j].second));
#endif
    }
  }
  solver->addClause(d_cl);
}

void BLS::blockDominatedRegion(const YPoint& yp){
  int nObj = yp.size();
  //computing the indexes given the objective value. This should be abstracted away...
  uint64_t objix[nObj];
  for(int iObj = 0; iObj < nObj  ; iObj++){
    objix[iObj] = 0;
    if(objRootLits[iObj])
      for(auto const& el: (*objRootLits[iObj]))
	if(yp[iObj] >= el.first) 
	  objix[iObj]++;
	else break;
  }

  blockDominatedRegion(objix, yp.size());
}


void BLS::activateLit(Lit lit){
  vec<Lit> d_cl;
  d_cl.push(lit);
  solver->addClause(d_cl);
}


void BLS::getObjRootLits(uint64_t * objv, uint64_t * objix, uint64_t * exact_objix, int nObj){
  Lit l;
  for(int i = 0; i < nObj; i++){
    if(objix[i] > 0 && (*objRootLits[i])[objix[i]].first != objv[i]){
      l = kps[i].getValRootLit(solver, objv[i]-1);
      exact_objix[i] = objix[i]+1;
      objRootLits[i]->insert(exact_objix[i], std::pair<uint64_t, Lit>(objv[i], l));
    }else{
      exact_objix[i] = objix[i];
    }
  }
  /*    
	printf("\n-> objix:");
	for(int di = 0; di < nObj; di++){
        printf(" %lu", objix[di]);
	}
	printf("\n-> objix2:");
	for(int di = 0; di < nObj; di++){
        printf(" %lu", exact_objix[di]);
	}
	printf("\n");
  */
    
  //   printRootLit(getFormula()->nObjFunctions());
    
}


void BLS::getNewObjRootLitsAndBlock(uint64_t * objv, uint64_t * objix, int nObj){
  Lit l;
  std::vector<std::pair<uint64_t, Lit>>::iterator it;
  vec<Lit> d_cl;
    
  for(int i = 0; i < nObj; i++){
    if(objix[i] > 0 && (*objRootLits[i])[objix[i]].first != objv[i]){
      l = kps[i].getValRootLit(solver, objv[i]-1);
            
            
      d_cl.push(l);
#ifdef __DEBUG__
      printf("c block2-1 %d [%lu][var: %d]\n", i, objv[i]-1, var(l));
#endif
            
    }else{
                
      if(objix[i] > 0){
	d_cl.push((*objRootLits[i])[objix[i]].second);
#ifdef __DEBUG__
	printf("c block2-2 %d [%lu][var: %d]\n", i,  objRootLits[i][objix[i]].first, var(objRootLits[i][objix[i]].second));
#endif
      }
    }
  }
        
  solver->addClause(d_cl);    
}




bool BLS::apObjectivesAreExact(){
  assert(approxMode == _ap_coeffs_);
  printf("c [BLS] is the approx function exact?\n");
    
  for(int di = 0; di < getFormula()->nObjFunctions(); di++){
    if(enc_is_kp_based() 
       && kps[di].hasCreatedEncoding() 
       && !kps[di].apFunctionIsExact()){
                
      printf("c No (%d)!\n", di);
      return false;
    }
  }
  printf("c Yes!\n");
  return true;
}



// really with MCSes (12/2020)
bool BLS::MCSSearchMO(bool permanentBlock) {
#ifdef __DEBUG__
  printf("c MCSSearchMO\n");
#endif
  lbool res = l_True;

  assumptions.clear();
  addAssumptMOFormulation();
  
  
  double runtime = cpuTime();
  
  //   uint64_t costModel = _maxWeight;
  int nObj = getFormula()->nObjFunctions();
  
  vec<Lit> d_cl;
  uint64_t objix[nObj];
  uint64_t lastobjix[nObj];
  uint64_t ap_objv[nObj];
  uint64_t objv[nObj];
  uint64_t last_ndom[nObj];
  
  bool terminate = false;
  bool foundSolution = false;
  
  for(int i = 0; i < nObj; i++){
    lastobjix[i] = objRootLits[i]->size()-1;
  }
    
  //   assumptions.clear();
  Lit blockLit_focus = mkLit(solver->nVars(), false);
  newSATVariable(solver);
  Lit blockLit_dom = mkLit(solver->nVars(), false);
  newSATVariable(solver);
  
  assumptions.push(~blockLit_focus);
  assumptions.push(~blockLit_dom);
  
  if(permanentBlock){
    //so bloqueia a ultima solucao que e a que foi encontrada na ultima sat call e ainda nao foi
    //bloqueada
    evalModel(solver->model, objv, ap_objv);
    //para actualizar o lower bound set
    evalToIndex(ap_objv, objix);
    getPoint(objix, nObj, last_ndom);
    fixLBpoint(last_ndom, objv, nObj);
    updateLowerBoundSet(last_ndom, true, nreencodes);
    // e bloquear regiao dentro da aprox
    tmpBlockDominatedRegion(objix, nObj, blockLit_dom);
        
    foundSolution = true; // para que o a regiao dominada por esta solucao seja bloqueada mais abaixo
  
  }else if(effsols.size() > 0){
    //if the objective functions were reencoded
    for(size_t i = 0; i < effsols.size(); i++){
      //TODO: Tornar isto mais eficiente! 
      vec<lbool> model;
      lbool * modeli = effsols[i].first;
      for(int j = 0; j < getFormula()->nVars(); j++)
	model.push(modeli[j]);
      //TODO: Evitar re-avaliar o modelo... isto pode ser demorado...
      for(int di = 0; di < getFormula()->nObjFunctions(); di++)
	kps[di].evalModel(model, objv, ap_objv, di);
      evalToIndex(ap_objv, objix);

      //para actualizar o lower bound set
      getPoint(objix, nObj, last_ndom);
      fixLBpoint(last_ndom, objv, nObj);
      updateLowerBoundSet(last_ndom, true, nreencodes);
          
      tmpBlockDominatedRegion(objix, nObj, blockLit_dom);
#ifdef __DEBUG__
      printf("\n");
#endif
    }
  }

  //   solver->my_print();
    
  while(!terminate) {
    

      
    while(res == l_True) {
        
      res = solve();
#ifdef __DEBUG__
      printf("c solved\n");
#endif
      if(res == l_Undef) {
	printf("c Warn1: SAT Solver exit due to conflict budget.\n");
	return false;
      }
      if(res == l_True){
#ifdef __DEBUG__
	printf("c sat\n");
#endif
	saveModel(solver->model); //include o nbSatisfiable++

	answerType = _SATISFIABLE_;
        //   solver->my_print();
	foundSolution = true;
            
	//eval
	evalModel(solver->model, objv, ap_objv);
	evalToIndex(ap_objv, objix);
            
#ifdef __DEBUG__
	printf("objv:");
	for(int di = 0; di < nObj; di++){
	  printf(" %lu", objv[di]);
	}
	printf("\nobjix:");
	for(int di = 0; di < nObj; di++){
	  printf(" %lu", objix[di]);
	}
	printf("\n");
                
	printf("lastbjix:");
	for(int di = 0; di < nObj; di++){
	  printf(" %lu", lastobjix[di]);
	}
	printf("\n");
            
#endif    
    
	printf("c o");
	for(int di = 0; di < nObj; di++){
	  printf(" %lu", objv[di]);
	}
	printf("\t[");
	for(int di = 0; di < nObj; di++){
	  printf(" %ld", objv[di]+getFormula()->getObjFunction(di)->_const);
	}
	printf("]\n");
	fflush(stdout);
            
	//TODO: Fazer de forma diferente, para adicionar apenas pontos nao dominados,
	//e, no maximo, adicionar um ponto possivelmente nao optimo (a ultima solucao gerada
	//antes da interrupcao)
	uint64_t last_ndom[nObj];
	getPoint(objix, nObj, last_ndom);
        
	fixLBpoint(last_ndom, objv, nObj);
	updateLowerBoundSet(last_ndom, true, nreencodes);
	saveEfficientSol(solver->model, objv, true);
                      
	update_assumpt_n_mcs(blockLit_focus, objix, lastobjix, nObj);
                            
	for(int i = 0; i < nObj; i++)
	  lastobjix[i] = objix[i];
            
      }else{
#ifdef __DEBUG__
	printf("c unsat\n");
#endif
      }
        
    }
    
    if(foundSolution){
        
      printf("o");
      for(int di = 0; di < nObj; di++){
	printf(" %ld", objv[di]+getFormula()->getObjFunction(di)->_const);
      }
      printf("\n");
      fflush(stdout);
        
      runtime = cpuTime();
      printf("c new optimal solution (time: %.3f)\n", runtime - initialTime);
        

      activateLit(blockLit_focus);
        
      //add new tmp clause
      blockLit_focus = mkLit(solver->nVars(), false);
      newSATVariable(solver);
      assumptions.clear();
      addAssumptMOFormulation();
      assumptions.push(~blockLit_focus);
      assumptions.push(~blockLit_dom);
        
      if(epsilon > 1 && approxMode == _ap_outvars_){
	//block region dominated by the last solution found
	tmpBlockDominatedRegion(objix, nObj, blockLit_dom);
            
	if(permanentBlock){
	  printf("c mandatory block\n");
	  getNewObjRootLitsAndBlock(objv, objix, nObj);
	}
      }
      else{ 
	blockDominatedRegion(objix, nObj);
      }
        
      res = l_True;
      terminate = false;
      foundSolution = false;
                
      for(int i = 0; i < nObj; i++){
	lastobjix[i] = (*objRootLits[i]).size()-1;
      }

    }else{
      terminate = true;
    }

  }
  
  assumptions.clear();
  activateLit(blockLit_dom);
  activateLit(blockLit_focus);
  

  //   epsilon = 1;
  if (nondom.size() == 0) {
    printf("c UNSAT\n");
    answerType=_UNSATISFIABLE_;
    return false;
  }else{
    answerType = _OPTIMUM_;
  }
  return true;
  //     return false;

}




void BLS::update_assumpt_n_mcs(Lit mcsesBlockLit, uint64_t * objix, uint64_t * lastobjix, int nObj){
  //update block
  vec<Lit> d_cl;
  d_cl.push(mcsesBlockLit);
  for(int i = 0; i < nObj; i++){
        
    if(lastobjix[i] > 0){
      uint64_t j;
      for(j = lastobjix[i]; j > 0 && j > objix[i]; j--){
	assumptions.push(((*objRootLits[i])[j].second));            
      }
      //             printf("i,j: %d %d\n", i,j);
      for(; j > 0; j--){
	d_cl.push((*objRootLits[i])[j].second);
      }
    }
  }
  solver->addClause(d_cl);
}




// ------------------------------------------------------------------- //
// ------------------------------------------------------------------- //


// Always consider all clauses by adding assumption var to each soft clause.

void BLS::basicSearch(int maxMCS = 50, bool maxSAT = false) {
  // Init Structures
  init();
  
  //Build solver
  solver = buildSolver();
  
  while (maxMCS) {
    if (!findNextMCS()) break;
    
    maxMCS--;
    if (maxSAT && _smallestMCS <= _lbWeight) {
      printAnswer(_OPTIMUM_);
      exit(_OPTIMUM_);      
    }
  }
  if (maxMCS == 0) 
    printf("c All requested MCSs found\n");
}


// Find next MCS.
// Returns false if the SAT solver was not able to finish. Otherwise, returns true.

bool BLS::findNextMCS() {
  vec<int> undefClauses;
  vec<int> satClauses;
  vec<int> unsatClauses;
  uint64_t costModel = _maxWeight;
  lbool res = l_True;
  int conflict_limit = 100000;
  
  initUndefClauses(undefClauses);
  assumptions.clear();
  
  solver->setConfBudget(conflict_limit);
  
  // make first call.
  res = solve();
  
  // Check outcome of first call
  if (res == l_False) {
    // Hard clause set in unsat!
    if (nbMCS == 0) {
      printAnswer(_UNSATISFIABLE_);
      exit(_UNSATISFIABLE_);
    }
    else {
      // It is not the first MCS. Hence, all MCS were found.
      printAnswer(_OPTIMUM_);
      exit(_OPTIMUM_);
    }
  }
  else if (res == l_True) saveModel(solver->model);
  else {
    printf("c Warn: SAT Solver exit due to conflict budget.\n");
    return false;
  }
  
  // Iterate to find next MCS
  while (undefClauses.size()) {
    int i = 0;
    if (res == l_True) {
      //Remove satisfied soft clauses from undefClauses
      
      while (i < undefClauses.size()) {
      	if (satisfiedSoft(undefClauses[i])) {
      	  //Add soft clause as Hard!
      	  satClauses.push(undefClauses[i]);
	  assumptions.push(~(getFormula()->getSoftClause(undefClauses[i]).assumption_var));
      	  costModel -= getFormula()->getSoftClause(undefClauses[i]).weight;
      	  undefClauses[i] = undefClauses.last();
      	  undefClauses.pop();
      	}
      	else i++;
      }
      
      if (costModel < _smallestMCS) {
	saveSmallestModel(solver->model);
	printf("o %ld\n", costModel);
	fflush(stdout);
      }
    }
    
    if (undefClauses.size() == 0) {
      break;
    }
    
    i = undefClauses.last();
    undefClauses.pop();
    satClauses.push(i);
    assumptions.push(~(getFormula()->getSoftClause(i).assumption_var));

    res = solve();
    
    
    if (res == l_False) {
      unsatClauses.push(satClauses.last());
      satClauses.pop();
      assumptions.pop();
    }
    else if (res == l_True) {
      saveModel(solver->model);
      costModel -= getFormula()->getSoftClause(satClauses.last()).weight;
      if (undefClauses.size() == 0 && costModel < _smallestMCS) {
        printf("o %ld\n", costModel);
        fflush(stdout);
      }
    }
    else {
      printf("c Warn: SAT Solver exit due to conflict budget.\n");
      return false;
    }
  }
  
  if (costModel < _smallestMCS) {
    _smallestMCS = costModel;
    //Last saved model is smallest MCS...
  }
  nbMCS++;
  printf("c MCS #%d Weight: %ld\n", nbMCS, costModel);
  
  // collect some information to compute hamming distance
  for (int i = 0; i < _assigned_true.size(); i++){
    if (solver->model[i] == l_True){
      _assigned_true[i]++;
    }
  }
  
  addMCSClause(unsatClauses);
  return true;
}


//fubs - function upper bounds (may be to restrictive and may exclude some opt sols)
void BLS::updateMOEncoding(){
  nreencodes++;
  if(enc_is_kp_based() || encoder.getPBEncoding() == _PB_GTE_ || encoder.getPBEncoding() == _PB_IGTE_){
    printf("c updateMOEncoding\n");
    //         solver->my_print();
    //         assumptions.clear();
    for(int i = 0; i < getFormula()->nObjFunctions(); i++){
      rootLits::RootLits ith_orl = rootLits::RootLits{};

      //             gtes.push(GTE());
      //             printf("c %d needs reencoding? RHS: %lu, fubs: %lu\n", i, igtes[i].getRHS(), fubs[i]);
      if(enc_is_kp_based()){
                    
	if(objRootLits[i] != NULL && (*objRootLits[i]).size() > 0){
	  //                     objRootLits[i].clear();
	  printf("c clear encoding of obj. funct. %d\n", i);
	  if(enc_is_kp_based()){
	    //                         kps[i].clearedEncoding(solver); //TODO
	  }
	}
                
	printf("\nc encode (function %d upper bound: %lu)\n", i, fubs[i]);
	getFormula()->getObjFunction(i)->my_print(getFormula()->getIndexToName());
                
	wlit_mapt rootLits;
	if(getFormula()->getObjFunction(i)->_lits.size() >= 1){
	  if(enc_is_kp_based()){
	    kps[i].setApproxRatio(solver, epsilon, approxMode);
	    //                         kps[i].setNameToIndex(getFormula()->getNameToIndex());
	    kps[i].encode(solver, getFormula()->getObjFunction(i)->_lits, getFormula()->getObjFunction(i)->_coeffs, fubs[i]);
	  }else if(encoder.getPBEncoding() == _PB_GTE_){
	    gtes[i].encode(solver, getFormula()->getObjFunction(i)->_lits, getFormula()->getObjFunction(i)->_coeffs, fubs[i]);
	  }
                    
	  //             printf("get root lits\n");
	  //TODO: get root literals (sorted), save them, and create order encoding
	  //                 printf("Root Literals (o.f. %d)\n", i+1);
	  //                 printRootLit(i);
	  if(enc_is_kp_based())
	    rootLits = kps[i].getRootLits();
	  else if(encoder.getPBEncoding() == _PB_GTE_)
	    rootLits = gtes[i].getRootLits();
                   
	}/* else{ */
	/*   rootLits.insert(wlit_pairt(getFormula()->getObjFunction(i)->_coeffs[0], getFormula()->getObjFunction(i)->_lits[0])); */
	/* } */
	int o = 0;
	vec<Lit> clause;
	//order encoding
	if(enc_is_kp_based()){
	  ith_orl.clear();
	  ith_orl.reserve(rootLits.size());
	  wlit_mapt::iterator prev;
	  for (wlit_mapt::iterator rit = rootLits.begin(); rit != rootLits.end(); rit++){
	    ith_orl.push_back(std::pair<uint64_t, Lit>(rit->first+1, rit->second));
	    invObjRootLits.insert({var(rit->second), i});
	  }
                    
                    
	}else{
	  wlit_mapt::reverse_iterator prev;
	  //se x_d = 1, entao a sol tera de ser <= d
	  for (wlit_mapt::reverse_iterator rit = rootLits.rbegin(); rit != rootLits.rend(); rit++){
            //                 objRootLits[i].push_back(std::pair<uint64_t, Lit>(rit->first, rit->second));
#ifdef __DEBUG__
	    printf("\t[%lu] ",rit->first);
	    if (sign(rit->second))
	      printf("~");
	    printf("y%d\n", var(rit->second)+1);
#endif
	    if(o > 0){
	      clause.push(~prev->second);
	      clause.push(rit->second);
	      solver->addClause(clause);
	      clause.clear();
	    }
	    prev = rit;
	    o++;
	  }
                    
                    
	  //             printf("Inverse literals:\n");
	  size_t j=0;
	  Lit p;
                    
	  for (wlit_mapt::iterator rit = rootLits.begin(); rit != rootLits.end(); rit++, j++){
	    if(ith_orl.size() > j+1){
	      //                         printf("Lit %lu %lu\n", rit->first, objRootLits[i][j+1].first);
	      p = ith_orl[j+1].second;
	    }else{
	      //                         printf("Lit %lu\n", rit->first);
	      p = mkLit(solver->nVars(), false);
	    }
	    newSATVariable(solver);
	    //~p \/ ~rit->second
	    clause.push(~rit->second);
	    clause.push(~p);
	    solver->addClause(clause);
	    clause.clear();
	    //p \/ rit->second)
	    clause.push(rit->second);
	    clause.push(p);
	    solver->addClause(clause);
	    clause.clear();
#ifdef __DEBUG__
	    printf("\t[%lu] ",rit->first);
	    printf("y%d - z%d\n", var(rit->second)+1, var(p)+1);
#endif                
	    //                 objRootLits[i].push_back(std::pair<uint64_t, Lit>(rit->first, rit->second));
	    if(ith_orl.size() <= j+1){
	      ith_orl.push_back(std::pair<uint64_t, Lit>(rit->first, p));
	    }
	  }
	}
                
      }
      objRootLits[i] = std::make_unique<rootLits::RootLits>(std::move(ith_orl));

    }
  }

#ifdef __DEBUG__
  printRootLit(getFormula()->nObjFunctions());
#endif                
  // exit(1);
}

shared_ptr<MOCOFormula> BLS::workFormula(){return wf;}
bool BLS::loadWorkFormula(){
  for(int i = 0, n = workFormula()->nHard(); i < n; i++)
    solver->addClause(workFormula()->getHardClause(i).clause);
  return true;
}

//adiciona as assumptions, todas as variaveis de relaxacao associadas as funcoes objectivo
// (que, por exemplo, permitem desactivar o encoding)
void BLS::addAssumptMOFormulation(bool omitRHS){
  //     printf("c addAssumptMOFormulation\n");

  if(enc_is_kp_based()){
    for(int i = 0; i < getFormula()->nObjFunctions(); i++){ 
      std::vector<Lit> assumpts;
      assumpts = kps[i].getAssumptions();
      for(size_t z = 0; z < assumpts.size(); z++){
	assumptions.push(assumpts[z]);
	//                 printf("a %lu: %d\n", i, var(assumpts[z]));
      }
        
    }
  }
}



void BLS::getPoint(uint64_t * ix, int d, uint64_t * out){
  //     uint64_t * pt = (out == NULL) ? new uint64_t[d] : out;
  //     uint64_t * pt ;
  for(int di = 0; di < d; di++){
    // //         pt[di] = objRootLits[di][ix[di]-1].first;
    //         printf("[BLS::getPoint] ix[%d]: %lu\n", di, ix[di]);
    //         printf("objRootLits[%d].size(): %lu\n", di, objRootLits[di].size());
    if(ix[di] >= (*objRootLits[di]).size()){
      out[di] = UINT_MAX;
    }else{
      out[di] = (*objRootLits[di])[ix[di]].first;
    }
        
    //         printf("objval: %lu\n", pt[di]);
  }
  //     return pt;
}

void BLS::fixLBpoint(uint64_t * lb, uint64_t * ptdom, int d){
  for(int di = 0; di < d; di++){
    //importante por causa dos pontos no LBset (e para calcular o epsilon indicator)
    if(lb[di] == 0 && enc_is_kp_based()){
      lb[di] = kps[di].getMinimum();
    }
    if(lb[di] == 0 && ptdom[di] > 0){
      lb[di] = ptdom[di];
    }
  }
}


void BLS::printRootLit(int d){
  size_t maxvs = 5;
  for(int di = 0; di < d; di++){
    printf("c [BLS] printRootLit[%d] (%lu)\n", di, (*objRootLits[di]).size());
    //         printf("obj %d:", di);
    for(size_t i = 0; i < (*objRootLits[di]).size(); i++){
      if(i < maxvs || i > (*objRootLits[di]).size() - maxvs)
	printf(" [%lu] z%d", (*objRootLits[di])[i].first, var((*objRootLits[di])[i].second)+1);
      else if (i == maxvs) printf(" ...");
    }
    printf("\n");
  }
}


void BLS::printRootLitModel(int d){
  printf("c [BLS] printRootLitModel\n");
  for(int di = 0; di < d; di++){
    size_t i = 0;
    printf("obj %d:", di);
    while(i < (*objRootLits[di]).size()){
      if(solver->modelValue(var((*objRootLits[di])[i].second)) == l_True)
	printf(" [%lu] z%d", (*objRootLits[di])[i].first, var((*objRootLits[di])[i].second));
      else
	printf(" [%lu] ~z%d", (*objRootLits[di])[i].first, var((*objRootLits[di])[i].second));
      i++;
    }
    printf("\n");
  }
}


// Public search method
StatusCode BLS::search() {
  
  search_MO(sstrategy);
  printStats();
  printf("done\n");
  // Make sure the conflict budget is turned off.
  solver->budgetOff();  
  return answerType;
}


/************************************************************************************************
 //
 // Utils for printing
 //
 ************************************************************************************************/ 



void BLS::evalModel(vec<lbool> &currentModel, uint64_t * objv, uint64_t *ap_objv){
  Lit l;
  for(int di = 0; di < getFormula()->nObjFunctions(); di++){
        
        
    if(enc_is_kp_based() && kps[di].hasCreatedEncoding()){
      kps[di].evalModel(currentModel, objv, ap_objv, di);
    }else{
      objv[di] = 0;
      auto& objf = getFormula()->getObjFunction(di);
      for(int i = 0; i < objf->_lits.size(); i++){
	l = objf->_lits[i];
                
	if ((sign(l) && currentModel[var(l)] == l_False) || 
	    (!sign(l) && currentModel[var(l)] == l_True)) {
	  objv[di] += objf->_coeffs[i];
	  //                 printf(" %cx%d", (sign(l)) ? '-':' ', var(l)+1);
	}
      }
      ap_objv[di] = objv[di];
      //         printf("\n");
    }
  }
}

void BLS::evalToIndex(const YPoint& yp, uint64_t * objix){
    
  for(uint64_t di = 0; di < yp.size(); di++){
    objix[di] = 0;
    size_t i = 1;
    while(i < (*objRootLits[di]).size() && yp[di]+1 >= (*objRootLits[di])[i].first){
      objix[di] = i;
      i++;
    }
        
  }
}

void BLS::evalToIndex(uint64_t * objv, uint64_t * objix){
    
  for(int di = 0; di <  getFormula()->nObjFunctions(); di++){
    objix[di] = 0;
    size_t i = 1;
    if(objRootLits[di]!=NULL)
    while(i <= (*objRootLits[di]).size() && objv[di] >= (*objRootLits[di])[i].first){
      objix[di] = i;
      i++;
    }
        
  }
}





// Prints the best satisfying model. Assumes that 'model' is not empty.
void BLS::printModel(){
  
  assert (model.size() != 0);
  
  printf("v ");
  for (int i = 0; i < model.size(); i++){
    if (model[i] == l_True) printf("%d ",i+1);
    else printf("%d ",-(i+1));
  }
  printf("\n");
}


// Prints the best satisfying model. Assumes that '_smalledtModel' is not empty.
void BLS::printSmallestModel(){
  assert (_smallestModel.size() != 0);
  
  printf("v ");
  for (int i = 0; i < _smallestModel.size(); i++){
    if (_smallestModel[i] == l_True) printf("%d ",i+1);
    else printf("%d ",-(i+1));
  }
  printf("\n");
}


// Prints search statistics.
void BLS::printStats(){
  //double totalTime = cpuTime();
  
  printf("c\n");
  //printf("c  Best solution:          %12"PRIu64"\n", ubCost);
  //printf("c  Total time:             %12.2f s\n",totalTime - initialTime);
  printf("c  Nb SAT solver calls:    %12d\n", nbSatCalls);
  printf("c  Nb SAT calls:           %12d\n", nbSatisfiable);
  printf("c  Nb UNSAT calls:         %12d\n", nbSatCalls-nbSatisfiable);
  printf("c  Nb MCS:                 %12d\n", nbMCS);
  printf("c  Smallest MCS:           %12ld\n", _smallestMCS);
  printf("c  Lower Bound:            %12ld\n", _lbWeight);
  //printf("c  Average core size:      %12.2f\n", (float)sumSizeCores/nbCores);
  printf("c\n"); 
}

/*
// Prints the corresponding answer.
void BLS::printAnswer(int type){

if (verbosity > 0)
printStats();

if (type == _UNKNOWN_ && model.size() > 0)
type = _SATISFIABLE_;
  
switch(type){
case _SATISFIABLE_:
printf("s SATISFIABLE\n");
printSmallestModel();
break;
case _OPTIMUM_:
printf("c All requested MCSs found\n");
printf("s OPTIMUM\n");
printSmallestModel();
break;
case _UNSATISFIABLE_:
printf("s UNSATISFIABLE\n");
break;  
case _UNKNOWN_:
printf("s UNKNOWN\n");
break;
default:
printf("c Error: Unknown answer type.\n");
}
}*/

/************************************************************************************************
 //
 // Other protected methods
 //
 ************************************************************************************************/ 

Solver *BLS::buildSolver() {

  vec<bool> seen;
  seen.growTo(getFormula()->nVars(), false);

  Solver *S = newSATSolver();
  
  for (int i = 0; i < getFormula()->nVars(); i++)
    newSATVariable(S);
  
  for (int i = 0; i < getFormula()->nHard(); i++)
    S->addClause(getFormula()->getHardClause(i).clause);

  vec<Lit> clause;
  for (int i = 0; i < getFormula()->nSoft(); i++) {
    
    clause.clear();
    getFormula()->getSoftClause(i).clause.copyTo(clause);
    
    for (int j = 0; j < getFormula()->getSoftClause(i).relaxation_vars.size();
	 j++) {
      clause.push(getFormula()->getSoftClause(i).relaxation_vars[j]);
    }
    
    S->addClause(clause);
  }
  
  return S;
}



Solver *BLS::buildSolverMO() {

  printf("c [BLS] buildSolverMO -- nVars: %d nHard: %d nSoft: %d nPB: %d nObj: %d\n",
	 getFormula()->nVars(), getFormula()->nHard(), getFormula()->nSoft(),
	 getFormula()->nPB(), getFormula()->nObjFunctions());
    
  vec<bool> seen;
  seen.growTo(getFormula()->nVars(), false);

  Solver *S = newSATSolver();
  
  for (int i = 0; i < getFormula()->nVars(); i++)
    newSATVariable(S);
  
  for (int i = 0; i < getFormula()->nHard(); i++)
    S->addClause(getFormula()->getHardClause(i).clause);

  //   printf("c Encode PB constraints\n");

  for (int i = 0; i < getFormula()->nPB(); i++) {
    Encoder *enc = new Encoder(_INCREMENTAL_NONE_, encoding,
                               _AMO_LADDER_, pb_encoding); //AG

    assert(getFormula()->getPBConstraint(i)->_sign);

    enc->encodePB(S, getFormula()->getPBConstraint(i)->_lits,
                  getFormula()->getPBConstraint(i)->_coeffs,
                  getFormula()->getPBConstraint(i)->_rhs);
    delete enc;
  }

  //   printf("c Encode cardinality constraints\n");
  for (int i = 0; i < getFormula()->nCard(); i++) {
    //     Encoder *enc = new Encoder(_INCREMENTAL_NONE_, _CARD_MTOTALIZER_,
    //                                _AMO_LADDER_, _PB_GTE_); //original
    Encoder *enc = new Encoder(_INCREMENTAL_NONE_, encoding,
                               _AMO_LADDER_, pb_encoding); //encoding

    if (getFormula()->getCardinalityConstraint(i)->_rhs == 1) {
      enc->encodeAMO(S, getFormula()->getCardinalityConstraint(i)->_lits);
    } else {
#ifdef __DEBUG__
      printf("c [BLS] buildSolverMO] encodeCardinality() constraint\n\t");
      getFormula()->getCardinalityConstraint(i)->my_print(getFormula()->getIndexToName());
#endif
      
      enc->encodeCardinality(S,
                             getFormula()->getCardinalityConstraint(i)->_lits,
                             getFormula()->getCardinalityConstraint(i)->_rhs);
    }

    delete enc;
  }

  vec<Lit> clause;
  
  //   printf("c Encode soft constraints\n");
  for (int i = 0; i < getFormula()->nSoft(); i++) {
    //       DPRINTF("c [Alg_LinearSU: rebuildSolver] add soft clause to SAT solver\n");
      
    //     if (getFormula()->getSoftClause(i).weight < min_weight)
    //       continue;

    clause.clear();
    getFormula()->getSoftClause(i).clause.copyTo(clause);

    for (int j = 0; j < getFormula()->getSoftClause(i).relaxation_vars.size();
         j++) {
      clause.push(getFormula()->getSoftClause(i).relaxation_vars[j]);
    }

    S->addClause(clause);
  }
  
  //   printf("c Encode objective functions\n");
  int nObj = getFormula()->nObjFunctions();
  
  for(int di = 0; di < nObj; di++) fubs[di] = 0;
  getMaxSATFormula()->sync_first(S);
  encoder.kpa_fixed_vars(getMaxSATFormula()->fixed_vars());
  return S;
}




uint64_t BLS::getTighterUB(int di){
  return getFormula()->getTighterUB(di);
}

uint64_t BLS::getTighterLB(int di){
  return  getFormula()->getTighterLB(di);
}
bool BLS::firstSolution(){
  
  double before_1stsol = cpuTime();
  //   bool res = findNextMCS();
  printf("c first call to solver\n");
  lbool res = solve();
  double totaltime_1stsol = cpuTime() - before_1stsol;
  printf("c time for %s: %f\n", (res==l_True) ? "SAT" : "UNSAT", totaltime_1stsol);

  if(res != l_True)
    return false;
  Model m {solver->model};
  first = Solution::OneSolution{&solution(), m};
  std::cout << "c oi " << first.yPoint() << std::endl;
  return true;
}

//encodes the objective functions only if there is, at least, one feasible solution (and return true)
bool BLS::updateMOFormulationIfSAT(){
  printf("c [updateMOFormulationIfSAT]\n");
  //     solver->my_print();
  
  if(!firstSolution())
    return false;
  updateMOFormulation();
  return true;
}

bool BLS::updateMOFormulation(){
  printf("c [updateMOFormulation]\n");
    if(nreencodes == 0){
  int nObj = getFormula()->nObjFunctions();
    for(int di = 0; di < nObj; di++) 
      fubs[di] = getFormula()->getObjFunction(di)->getUB();
}
  double before_enc = cpuTime();
  updateMOEncoding();   //aqui e esquecido o encoding anterior, caso seja refeito
  double total_time = cpuTime();
  printf("c encode time: %f\n\n", total_time - before_enc);
  
  return true;
}



void BLS::init() {
  vec<int> vars;
  vars.growTo(getFormula()->nVars(), 0);
  
  if (!_useAllVars) {
    for (int i = 0; i < getFormula()->nSoft(); i++) {
      for (int j = 0; j < getFormula()->getSoftClause(i).clause.size(); j++) {
	int v = var(getFormula()->getSoftClause(i).clause[j]);
	assert (v < vars.size());
	if (vars[v] == 0){
	  _soft_variables.push(v);
	  vars[v] = 1;
	}
      }
    }
  }
  else {
    for (int i = 0; i < getFormula()->nVars(); i++) {
      _soft_variables.push(i);
    }
  }
  _assigned_true.growTo(getFormula()->nInitialVars(),0);
  _varScore.growTo(getFormula()->nInitialVars(),0);
  
  _maxWeight = 0;
  for (int i = 0; i < getFormula()->nSoft(); i++) {
    Lit l = getFormula()->newLiteral();
    getFormula()->getSoftClause(i).relaxation_vars.push(l);
    getFormula()->getSoftClause(i).assumption_var = getFormula()->getSoftClause(i).relaxation_vars[0]; // Assumption Var is relaxation var
    
    _maxWeight += getFormula()->getSoftClause(i).weight;
  }
  printf("c Max. Weight: %ld\n", _maxWeight);
  // build the objRootLits vector
  for(int i = 0; i < getFormula()->nObjFunctions(); i++)
    objRootLits.push_back(std::make_unique<rootLits::RootLits>(rootLits::RootLits{}));


}

void BLS::updateStats(){
  printf("c lobjstats obj nvars, nclauses, nrootvars\n");
  int nev = 0, nec = 0, nerv = 0;
  runstats[_nencvars_] = 0;
  runstats[_nencclauses_] = 0;
  runstats[_nencrootvars_] = 0;
  runstats[_nreencodes_] = nreencodes;
        
  for(int i = 0; i < getFormula()->nObjFunctions(); i++){
        
    if(enc_is_kp_based())
      kps[i].getEncodeSizes(&nev, &nec, &nerv);
    else if(encoder.getPBEncoding() == _PB_GTE_)
      gtes[i].getEncodeSizes(&nev, &nec, &nerv);
    runstats[_nencvars_] += nev;
    runstats[_nencclauses_] += nec;
    runstats[_nencrootvars_] += nerv;
    printf("c objstats %d %d %d %d\n", i+1, nev, nec, nerv);
  }
}
void BLS::consolidateSolution(){
  if(first.model().size())
    solution().pushSafe(first.model());
}
void BLS::build(){
  init();
  nbMCS = 0;
  answerType = _UNKNOWN_;
  solver = buildSolverMO(); 
}

/*assumes the region that dominates the point given by objix*/
void BLS::assumeDominatingRegion(uint64_t * objix, int nObj) {
  for(int di = 0; di < nObj; di++){
    int j = objix[di];
    if(j > 0){
      assumptions.push((*objRootLits[di])[j].second);
#ifdef __DEBUG__
      printf("assume all %d [%lu][var: %d]\n", di, (*objRootLits[di])[j].first, var((*objRootLits[di])[j].second));
#endif

    }
  }
}
/*assumes the region that dominates the point given by yp*/
void BLS::assumeDominatingRegion(const YPoint& yp) {
  auto nObj=yp.size();
  YPoint yp1 = yp;
  uint64_t objix[nObj];
  
  evalToIndex(yp1, objix);
  for(YPoint::size_type di = 0; di < nObj; di++){
    int j = objix[di];
    if(j > 0 && j < (int) (*objRootLits[di]).size()){
      assumptions.push((*objRootLits[di]).at(j).second);
#ifdef __DEBUG__
      printf("assume all %d [%lu][var: %d]\n", di, (*objRootLits[di])[j].first, var((*objRootLits[di])[j].second));
#endif

    }
  }
}

PBObjFunction  openwbo::filter(PBObjFunction* pb, set<Lit>& vars){
  map<Lit, uint64_t> terms{};  
  vec<Lit> lits{};
  vec<uint64_t> coeffs{};

  if(pb==nullptr)
    return PBObjFunction{};
  
  for(int j = 0; j < pb->_lits.size(); j++)
    terms[pb->_lits[j]]+=pb->_coeffs[j];
  
  for(auto& el: terms){
    if(vars.count(el.first)){
      lits.push(el.first);
      coeffs.push(el.second);
    }
  }
  return PBObjFunction{lits,coeffs,(int64_t) pb->_const};
}

void BLS::blockStep(const YPoint& yp){
  blockDominatedRegion(yp);
}
