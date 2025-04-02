#include "Alg_Naive.h"
#include <algorithm>    // std::max

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
void Naive::saveModel(vec<lbool> &currentModel){
  //assert (n_initial_vars != 0);
  assert (currentModel.size() != 0);
  
  model.clear();
  // Only store the value of the variables that belong to the original MaxSAT formula.
  for (int i = 0; i < maxsat_formula->nVars(); i++){
    model.push(currentModel[i]);
  }
  
  nbSatisfiable++;
}




void Naive::saveSmallestModel(vec<lbool> &currentModel){
  //assert (n_initial_vars != 0);
  assert (currentModel.size() != 0);
  
  _smallestModel.clear();
  // Only store the value of the variables that belong to the original MaxSAT formula.
  for (int i = 0; i < maxsat_formula->nVars(); i++){
    _smallestModel.push(currentModel[i]);
  }
}


/************************************************************************************************
 //
 // Naive search
 //
 ************************************************************************************************/ 


// Checks if a soft clause is satisfied by saved model

bool Naive::satisfiedSoft(int i) {
  for (int j = 0; j < maxsat_formula->getSoftClause(i).clause.size(); j++){
    assert (var(maxsat_formula->getSoftClause(i).clause[j]) < model.size());
    if ((sign(maxsat_formula->getSoftClause(i).clause[j]) && model[var(maxsat_formula->getSoftClause(i).clause[j])] == l_False) || 
	(!sign(maxsat_formula->getSoftClause(i).clause[j]) && model[var(maxsat_formula->getSoftClause(i).clause[j])] == l_True)) {
      return true;
    }
  }
  return false;
}





// Call to the SAT solver...

lbool Naive::solve() {
  nbSatCalls++;  
#ifdef SIMP
  return ((SimpSolver*)solver)->solveLimited(assumptions);
#else
  return solver->solveLimited(assumptions);
#endif
}



void Naive::initUndefClauses(vec<int>& undefClauses) {
  for (int i = 0; i < maxsat_formula->nSoft(); i++) 
    undefClauses.push(i);
}


void Naive::addBackboneLiterals(int index) {
  assumptions.push(maxsat_formula->getSoftClause(index).assumption_var);
  for (int j = 0; j < maxsat_formula->getSoftClause(index).clause.size(); j++){
    assumptions.push(~(maxsat_formula->getSoftClause(index).clause[j]));
  }
}

bool Naive::searchMO(){
  Solution dirty_sol{this};
  vec<Lit> clause{};
  while(solve() == l_True) {
    Model mdl =Model{solver->model};
    dirty_sol.push(mdl);
    std::cout <<  "c v " << mdl << std::endl;
    std::cout <<  "c o " << dirty_sol.yPoint() << std::endl;

    // block last found model
    {
      clause.clear();
      int i = 0;
      for(auto& at: mdl){
	bool sign = at == l_True? true : false;
	Lit lit = mkLit(i, sign);
	clause.push(lit);
	++i;
      }}
    solver->addClause(clause);
  }

  for(auto& osol: dirty_sol){
    Model mdl = osol.second.first.model();
      solution().pushSafe(mdl);
  }

  if (solution().size() == 0) {
    printf("c UNSAT\n");
    answerType=_UNSATISFIABLE_;
    return false;
  }else{
    answerType = _OPTIMUM_;
  }
  
  return true;
}

void Naive::search_MO(int strategy){
  // Init Structures
  init();
  solver = buildSolverMO();    
  printf("c search\n");
  searchMO();
  printf("c Done searching\n");
  printf("c epsilon: %f\n", epsilon);
  printf("c reductionFactor: %f\n", redFactor);
  if(nondom.size() > 0)
    answerType = _OPTIMUM_;
  else
    answerType = _UNSATISFIABLE_;
 
 printAnswer(answerType);
}




void Naive::activateLit(Lit lit){
  vec<Lit> d_cl;
  d_cl.push(lit);
  solver->addClause(d_cl);
}







void Naive::getPoint(uint64_t * ix, int d, uint64_t * out){
  //     uint64_t * pt = (out == NULL) ? new uint64_t[d] : out;
  //     uint64_t * pt ;
  for(int di = 0; di < d; di++){
    // //         pt[di] = objRootLits[di][ix[di]-1].first;
    //         printf("[Naive::getPoint] ix[%d]: %lu\n", di, ix[di]);
    //         printf("objRootLits[%d].size(): %lu\n", di, objRootLits[di].size());
    if(ix[di] >= objRootLits[di].size()){
      out[di] = UINT_MAX;
    }else{
      out[di] = objRootLits[di][ix[di]].first;
    }
        
    //         printf("objval: %lu\n", pt[di]);
  }
  //     return pt;
}

void Naive::fixLBpoint(uint64_t * lb, uint64_t * ptdom, int d){
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


void Naive::printRootLit(int d){
  size_t maxvs = 5;
  for(int di = 0; di < d; di++){
    printf("c [Naive] printRootLit[%d] (%lu)\n", di, objRootLits[di].size());
    //         printf("obj %d:", di);
    for(size_t i = 0; i < objRootLits[di].size(); i++){
      if(i < maxvs || i > objRootLits[di].size() - maxvs)
	printf(" [%lu] z%d", objRootLits[di][i].first, var(objRootLits[di][i].second)+1);
      else if (i == maxvs) printf(" ...");
    }
    printf("\n");
  }
}


void Naive::printRootLitModel(int d){
  printf("c [Naive] printRootLitModel\n");
  for(int di = 0; di < d; di++){
    size_t i = 0;
    printf("obj %d:", di);
    while(i < objRootLits[di].size()){
      if(solver->modelValue(var(objRootLits[di][i].second)) == l_True)
	printf(" [%lu] z%d", objRootLits[di][i].first, var(objRootLits[di][i].second));
      else
	printf(" [%lu] ~z%d", objRootLits[di][i].first, var(objRootLits[di][i].second));
      i++;
    }
    printf("\n");
  }
}


// Public search method
StatusCode Naive::search() {
  
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



void Naive::evalModel(vec<lbool> &currentModel, uint64_t * objv, uint64_t *ap_objv){
  Lit l;
  for(int di = 0; di < maxsat_formula->nObjFunctions(); di++){
        
        
    if(enc_is_kp_based() && kps[di].hasCreatedEncoding()){
      kps[di].evalModel(currentModel, objv, ap_objv, di);
    }else{
      objv[di] = 0;
      PBObjFunction * objf = maxsat_formula->getObjFunction(di);
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

void Naive::evalToIndex(YPoint& yp, uint64_t * objix){
    
  for(YPoint::size_type di = 0; di < yp.size(); di++){
    objix[di] = 0;
    size_t i = 1;
    while(i < objRootLits[di].size() && yp[di] >= objRootLits[di][i].first){
      objix[di] = i;
      i++;
    }
        
  }
}

void Naive::evalToIndex(uint64_t * objv, uint64_t * objix){
    
  for(int di = 0; di <  maxsat_formula->nObjFunctions(); di++){
    objix[di] = 0;
    size_t i = 1;
    while(i < objRootLits[di].size() && objv[di] >= objRootLits[di][i].first){
      objix[di] = i;
      i++;
    }
        
  }
}





// Prints the best satisfying model. Assumes that 'model' is not empty.
void Naive::printModel(){
  
  assert (model.size() != 0);
  
  printf("v ");
  for (int i = 0; i < model.size(); i++){
    if (model[i] == l_True) printf("%d ",i+1);
    else printf("%d ",-(i+1));
  }
  printf("\n");
}


// Prints the best satisfying model. Assumes that '_smalledtModel' is not empty.
void Naive::printSmallestModel(){
  assert (_smallestModel.size() != 0);
  
  printf("v ");
  for (int i = 0; i < _smallestModel.size(); i++){
    if (_smallestModel[i] == l_True) printf("%d ",i+1);
    else printf("%d ",-(i+1));
  }
  printf("\n");
}


// Prints search statistics.
void Naive::printStats(){
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
void Naive::printAnswer(int type){

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

Solver *Naive::buildSolver() {

  vec<bool> seen;
  seen.growTo(maxsat_formula->nVars(), false);

  Solver *S = newSATSolver();
  
  for (int i = 0; i < maxsat_formula->nVars(); i++)
    newSATVariable(S);
  
  for (int i = 0; i < maxsat_formula->nHard(); i++)
    S->addClause(maxsat_formula->getHardClause(i).clause);

  vec<Lit> clause;
  for (int i = 0; i < maxsat_formula->nSoft(); i++) {
    
    clause.clear();
    maxsat_formula->getSoftClause(i).clause.copyTo(clause);
    
    for (int j = 0; j < maxsat_formula->getSoftClause(i).relaxation_vars.size();
	 j++) {
      clause.push(maxsat_formula->getSoftClause(i).relaxation_vars[j]);
    }
    
    S->addClause(clause);
  }
  
  return S;
}



Solver *Naive::buildSolverMO() {

  printf("c [Naive] buildSolverMO -- nVars: %d nHard: %d nSoft: %d nPB: %d nObj: %d\n",
	 maxsat_formula->nVars(), maxsat_formula->nHard(), maxsat_formula->nSoft(),
	 maxsat_formula->nPB(), maxsat_formula->nObjFunctions());
    
  vec<bool> seen;
  seen.growTo(maxsat_formula->nVars(), false);

  Solver *S = newSATSolver();
  
  for (int i = 0; i < maxsat_formula->nVars(); i++)
    newSATVariable(S);
  
  for (int i = 0; i < maxsat_formula->nHard(); i++)
    S->addClause(maxsat_formula->getHardClause(i).clause);
  solver = S;
  //   printf("c Encode PB constraints\n");

  vec<Lit> clause;
  
  //   printf("c Encode soft constraints\n");
  for (int i = 0; i < maxsat_formula->nSoft(); i++) {
    //       DPRINTF("c [Alg_LinearSU: rebuildSolver] add soft clause to SAT solver\n");
      
    //     if (maxsat_formula->getSoftClause(i).weight < min_weight)
    //       continue;

    clause.clear();
    maxsat_formula->getSoftClause(i).clause.copyTo(clause);

    for (int j = 0; j < maxsat_formula->getSoftClause(i).relaxation_vars.size();
         j++) {
      clause.push(maxsat_formula->getSoftClause(i).relaxation_vars[j]);
    }

    S->addClause(clause);
  }
  
  //   printf("c Encode objective functions\n");
  int nObj = maxsat_formula->nObjFunctions();
  
  for(int di = 0; di < nObj; di++) fubs[di] = 0;
  return S;
}




void Naive::init() {
  vec<int> vars;
  vars.growTo(maxsat_formula->nVars(), 0);
  
    for (int i = 0; i < maxsat_formula->nVars(); i++) {
      _soft_variables.push(i);
    }

    _assigned_true.growTo(maxsat_formula->nInitialVars(),0);
    _varScore.growTo(maxsat_formula->nInitialVars(),0);
  
  _maxWeight = 0;
  for (int i = 0; i < maxsat_formula->nSoft(); i++) {
    Lit l = maxsat_formula->newLiteral();
    maxsat_formula->getSoftClause(i).relaxation_vars.push(l);
    maxsat_formula->getSoftClause(i).assumption_var = maxsat_formula->getSoftClause(i).relaxation_vars[0]; // Assumption Var is relaxation var
    
    _maxWeight += maxsat_formula->getSoftClause(i).weight;
  }
  printf("c Max. Weight: %ld\n", _maxWeight);
}

void Naive::updateStats(){
  printf("c lobjstats obj nvars, nclauses, nrootvars\n");
  int nev = 0, nec = 0, nerv = 0;
  runstats[_nencvars_] = 0;
  runstats[_nencclauses_] = 0;
  runstats[_nencrootvars_] = 0;
  runstats[_nreencodes_] = nreencodes;
        
  for(int i = 0; i < maxsat_formula->nObjFunctions(); i++){
        
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

void Naive::build(){
  init();
  nbMCS = 0;
  answerType = _UNKNOWN_;
  solver = buildSolverMO(); 
}
