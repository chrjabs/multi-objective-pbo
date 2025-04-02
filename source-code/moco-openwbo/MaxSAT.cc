/*!
 * \author Ruben Martins - ruben@sat.inesc-id.pt
 *
 * @section LICENSE
 *
 * Open-WBO, Copyright (c) 2013-2017, Ruben Martins, Vasco Manquinho, Ines Lynce
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 */

#include "MaxSAT.h"
#include "Pareto.h"
#include <sstream>
#include <float.h>

using namespace openwbo;


//does a weakly dominates b?
bool wdominates(const uint64_t * a, const uint64_t * b, int d){
    for(int i = 0; i < d; i++)
        if(a[i] > b[i])
            return false;
    return true;
}

/************************************************************************************************
 //
 // Public methods
 //
 ************************************************************************************************/

StatusCode MaxSAT::search() {
  if(print) printf("Error: Invalid MaxSAT algoritm.\n");
  throw MaxSATException(__FILE__, __LINE__, "Did not implement MaxSAT search");
  return _ERROR_;
}

void MaxSAT::setInitialTime(double initial) {
  initialTime = initial;
} // Sets the initial time.

/************************************************************************************************
 //
 // SAT solver interface
 //
 ************************************************************************************************/

// Creates an empty SAT Solver.
Solver *MaxSAT::newSATSolver() {

#ifdef SIMP
  NSPACE::SimpSolver *S = new NSPACE::SimpSolver();
#else
  Solver *S = new Solver();
#endif

  return (Solver *)S;
}

// Creates a new variable in the SAT solver.
void MaxSAT::newSATVariable(Solver *S) {

#ifdef SIMP
  ((NSPACE::SimpSolver *)S)->newVar();
#else
  S->newVar();
#endif
}

// Makes sure the underlying SAT solver has the given amount of variables
// reserved.
void MaxSAT::reserveSATVariables(Solver *S, unsigned maxVariable) {
#ifdef SAT_HAS_RESERVATION
#ifdef SIMP
  ((NSPACE::SimpSolver *)S)->reserveVars(maxVariable);
#else
  S->reserveVars(maxVariable);
#endif
#endif
}

// Solve the formula that is currently loaded in the SAT solver with a set of
// assumptions and with the option to use preprocessing for 'simp'.
lbool MaxSAT::searchSATSolver(Solver *S, vec<Lit> &assumptions, bool pre) {

// Currently preprocessing is disabled by default.
// Variable elimination cannot be done on relaxation variables nor on variables
// that belong to soft clauses. To preprocessing to be used those variables
// should be frozen.

#ifdef SIMP
  lbool res = ((NSPACE::SimpSolver *)S)->solveLimited(assumptions, pre);
#else
  lbool res = S->solveLimited(assumptions);
#endif
  nbSatCalls++;
  return res;
}

// Solve the formula without assumptions.
lbool MaxSAT::searchSATSolver(Solver *S, bool pre) {
  vec<Lit> dummy; // Empty set of assumptions.
  return searchSATSolver(S, dummy, pre);
}

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
  |
  |________________________________________________________________________________________________@*/
void MaxSAT::saveModel(vec<lbool> &currentModel) {
  assert(maxsat_formula->nInitialVars() != 0);
  assert(currentModel.size() != 0);

  model.clear();
  // Only store the value of the variables that belong to the
  // original MaxSAT formula.
  for (int i = 0; i < maxsat_formula->nInitialVars(); i++)
    model.push(currentModel[i]);
}

/*_________________________________________________________________________________________________
  |
  |  computeCostModel : (currentModel : vec<lbool>&) (weight : int) ->
  |                     [uint64_t]
  |
  |  Description:
  |
  |    Computes the cost of 'currentModel'. The cost of a model is the sum of
  |    the weights of the unsatisfied soft clauses.
  |    If a weight is specified, then it only considers the sum of the weights
  |    of the unsatisfied soft clauses with the specified weight.
  |
  |  Pre-conditions:
  |    * Assumes that 'currentModel' is not empty.
  |
  |________________________________________________________________________________________________@*/
uint64_t MaxSAT::computeCostModel(vec<lbool> &currentModel, uint64_t weight) {

  assert(currentModel.size() != 0);
  uint64_t currentCost = 0;

  for (int i = 0; i < maxsat_formula->nSoft(); i++) {
    bool unsatisfied = true;
    for (int j = 0; j < maxsat_formula->getSoftClause(i).clause.size(); j++) {

      if (weight != UINT64_MAX &&
          maxsat_formula->getSoftClause(i).weight != weight) {
        unsatisfied = false;
        continue;
      }

      assert(var(maxsat_formula->getSoftClause(i).clause[j]) <
             currentModel.size());
      if ((sign(maxsat_formula->getSoftClause(i).clause[j]) &&
           currentModel[var(maxsat_formula->getSoftClause(i).clause[j])] ==
               l_False) ||
          (!sign(maxsat_formula->getSoftClause(i).clause[j]) &&
           currentModel[var(maxsat_formula->getSoftClause(i).clause[j])] ==
               l_True)) {
        unsatisfied = false;
        break;
      }
    }

    if (unsatisfied) {
      currentCost += maxsat_formula->getSoftClause(i).weight;
    }
  }

  return currentCost;
}

/*_________________________________________________________________________________________________
  |
  |  isBMO : (cache : bool)  ->  [void]
  |
  |  Description:
  |
  |    Tests if the MaxSAT formula has lexicographical optimization criterion.
  |
  |  For further details see:
  |    * Joao Marques-Silva, Josep Argelich, Ana Graça, Inês Lynce: Boolean
  |      lexicographic optimization: algorithms & applications. Ann. Math.
  |      Artif. Intell. 62(3-4): 317-343 (2011)
  |
  |  Post-conditions:
  |    * 'orderWeights' is updated with the weights in lexicographical order if
  |      'cache' is true.
  |
  |________________________________________________________________________________________________@*/
bool MaxSAT::isBMO(bool cache) {
  assert(orderWeights.size() == 0);
  bool bmo = true;
  std::set<uint64_t> partitionWeights;
  std::map<uint64_t, uint64_t> nbPartitionWeights;

  for (int i = 0; i < maxsat_formula->nSoft(); i++) {
    partitionWeights.insert(maxsat_formula->getSoftClause(i).weight);
    nbPartitionWeights[maxsat_formula->getSoftClause(i).weight]++;
  }

  for (std::set<uint64_t>::iterator iter = partitionWeights.begin();
       iter != partitionWeights.end(); ++iter) {
    orderWeights.push_back(*iter);
  }

  std::sort(orderWeights.begin(), orderWeights.end(), greaterThan);

  uint64_t totalWeights = 0;
  for (int i = 0; i < (int)orderWeights.size(); i++)
    totalWeights += orderWeights[i] * nbPartitionWeights[orderWeights[i]];

  for (int i = 0; i < (int)orderWeights.size(); i++) {
    totalWeights -= orderWeights[i] * nbPartitionWeights[orderWeights[i]];
    if (orderWeights[i] < totalWeights) {
      bmo = false;
      break;
    }
  }

  if (!cache)
    orderWeights.clear();

  return bmo;
}

/************************************************************************************************
 //
 // Utils for printing
 //
 ************************************************************************************************/

// Prints information regarding the AMO encoding.
void MaxSAT::print_AMO_configuration(int encoding) {
  switch (encoding) {
  case _AMO_LADDER_:
    printf("c |  AMO Encoding:         %12s                      "
           "                                             |\n",
           "Ladder");
    break;

  default:
    printf("c Error: Invalid AMO encoding.\n");
    printf("s UNKNOWN\n");
    break;
  }
}

// Prints information regarding the PB encoding.
void MaxSAT::print_PB_configuration(int encoding) {
  switch (encoding) {
  case _PB_SWC_:
    printf("c |  PB Encoding:         %13s                        "
           "                                           |\n",
           "SWC");
    break;

  case _PB_GTE_:
    printf("c |  PB Encoding:         %13s                        "
           "                                           |\n",
           "GTE");
    break;

      case _PB_IGTE_:
    printf("c |  PB Encoding:         %13s                        "
           "                                           |\n",
           "IGTE");
    break;
  default:
    printf("c Error: Invalid PB encoding.\n");
    printf("s UNKNOWN\n");
    break;
  }
}

// Prints information regarding the cardinality encoding.
void MaxSAT::print_Card_configuration(int encoding) {
  switch (encoding) {
  case _CARD_CNETWORKS_:
    printf("c |  Cardinality Encoding: %12s                                "
           "                                   |\n",
           "CNetworks");
    break;

  case _CARD_TOTALIZER_:
    printf("c |  Cardinality Encoding: %12s                                "
           "                                   |\n",
           "Totalizer");
    break;

  case _CARD_MTOTALIZER_:
    printf("c |  Cardinality Encoding:    %19s                             "
           "                            |\n",
           "Modulo Totalizer");
    break;

  default:
    printf("c Error: Invalid cardinality encoding.\n");
    printf("s UNKNOWN\n");
    break;
  }
}

void MaxSAT::blockModel(Solver *solver) {
  assert(model.size() != 0);

  vec<Lit> blocking;

  printf("v ");
  for (int i = 0; i < model.size(); i++) {
    indexMap::const_iterator iter = maxsat_formula->getIndexToName().find(i);
    if (iter != maxsat_formula->getIndexToName().end()) {
      if (model[i] == l_False)
        printf("-");
      printf("%s ", iter->second.c_str());
    }
  }
  printf("\n");

  for (int i = 0; i < model.size(); i++) {
    blocking.push((model[i] == l_True) ? ~mkLit(i) : mkLit(i));
  }

  solver->addClause(blocking);
}

void MaxSAT::printBound(int64_t bound)
{
  if(!print) return;

  printf("o %" PRId64 "\n", bound);
}

// Prints the best satisfying model. Assumes that 'model' is not empty.
void MaxSAT::printModel() {
    
  assert(model.size() != 0);

  std::stringstream s;
  s << "v ";

  if (maxsat_formula->getFormat() == _FORMAT_PB_) {
    for (int i = 0; i < model.size(); i++) {
      indexMap::const_iterator iter = maxsat_formula->getIndexToName().find(i);
      if (iter != maxsat_formula->getIndexToName().end()) {
        if (model[i] == l_False)
          s << "-";
        s << iter->second.c_str() << " ";
      }
    }
  } else {
    for (int i = 0; i < model.size(); i++) {
      if (model[i] == l_True)
        s << i+1 << " ";
      else
        s << -(i+1) << " ";
    }
  }

  printf("%s\n", s.str().c_str());
}


void MaxSAT::printEffSolutions(bool printLBset) {
    
  assert(effsols.size() != 0);

  int d = maxsat_formula->nObjFunctions();

  
  FILE * f = stdout;
  if(print_model){
//     std::stringstream s;
        f = (print_my_output) ? fopen (effsols_file,"w") : stdout;
//         printf("c print to %s? %d\n", effsols_file, print_my_output);
        for(size_t j = 0; j < effsols.size(); j++){
            lbool *modelj = effsols[j].first;
            
//             s << "v ";
            if(!print_my_output)
                fprintf(f, "v");
            
            if (maxsat_formula->getFormat() == _FORMAT_PB_) {
                for (int i = 0; i < maxsat_formula->nVars(); i++) {
                    indexMap::const_iterator iter = maxsat_formula->getIndexToName().find(i);
                    if (iter != maxsat_formula->getIndexToName().end()) {
                        fprintf(f, " ");
                        if (modelj[i] == l_False)
//                         s << "-";
                            fprintf(f, "-");
                        fprintf(f, "%s", iter->second.c_str());
//                         s << iter->second.c_str() << " ";
                    }
                }
            } else {
              // Reconstruct model with preprocessing
              vec<lbool> rec_model{};
              for (int i = 0; i < maxsat_formula->nVars(); i++) {
                rec_model.push(modelj[i]);
              }
              maxsat_formula->reconstruct(rec_model);

              for (int i = 0; i < rec_model.size(); i++) {
                if (rec_model[i] == l_True)
//                     s << i+1 << " ";
                    fprintf(f, " %d", i+1);
                else
//                     s << -(i+1) << " ";
                    fprintf(f, " -%d", i+1);
                }
            }
//             printf("%s", s.str().c_str());
            fprintf(f, "\n");
        }
        
        if(print_my_output) fclose(f);
    }
    printf("c %lu (efficient) solutions\n", effsols.size());
    
    
//     if(!print_my_output){
        printf("c ------- \n");
        printf("c pts of transformed prob\n");
        for(size_t i = 0; i < nondom.size(); i++){
            printf("c pt");
            for(int di = 0; di < d; di++){
                printf(" %ld", nondom[i][di]); // + maxsat_formula->getObjFunction(di)->_const);
            }
            printf("\n");
        }
        printf("c ------- \n");
        printf("c %lu points T\n", nondom.size());
        printf("c ------- \n");
        
        printf("c lower bound set of transformed prob\n");
        for(size_t i = 0; i < LBset.size(); i++){
            printf("c lb");
            for(int di = 0; di < d; di++){
                printf(" %ld", LBset[i][di]); // + maxsat_formula->getObjFunction(di)->_const);
            }
            printf("\n");
        }
        printf("c ------- \n");
        printf("c %lu lbs T\n", LBset.size());
        printf("c ------- \n");
//     }
    
    //Lower bound set
//     if(!print_my_output){
    if(printLBset && LBset.size() > 0){
        f = (print_my_output) ? fopen (lbset_file,"w") : stdout;
        for(size_t i = 0; i < LBset.size(); i++){
            if(!print_my_output) fprintf(f, "c LBs");
            for(int di = 0; di < d; di++){
                fprintf(f, " %ld", LBset[i][di] + maxsat_formula->getObjFunction(di)->_const); // + maxsat_formula->getObjFunction(di)->_const);
            }
            fprintf(f, "\n");
        }
//         fprintf(f, "c ------- \n");
        printf("c %lu points in lower bound set\n", LBset.size());
        printf("c ------- \n");
        if(print_my_output) fclose(f);
    }
 
    
    // Solutions (objective space)
    f = (print_my_output) ? fopen (objv_file,"w") : stdout;
    
    for(size_t i = 0; i < nondom.size(); i++){
        if(!print_my_output)
            fprintf(f, "o");
        for(int di = 0; di < d; di++){
            fprintf(f, " %ld", nondom[i][di] + maxsat_formula->getObjFunction(di)->_const);
        }
        fprintf(f, "\n");
    }
    printf("c %lu nondominated points\n", nondom.size());
    if(print_my_output) fclose(f);

    printf("c _consts:");
    
    for(int di = 0; di < d; di++){
        printf(" %ld", maxsat_formula->getObjFunction(di)->_const);
    }
    printf("\n");
}



std::string MaxSAT::printSoftClause(int id){
  assert (maxsat_formula->getFormat() == _FORMAT_MAXSAT_);
  assert (id < maxsat_formula->nSoft());

  std::stringstream ss;
  ss << maxsat_formula->getSoftClause(id).weight << " ";

  for (int j = 0; j < maxsat_formula->getSoftClause(id).clause.size(); j++){
    if (sign(maxsat_formula->getSoftClause(id).clause[j]))
      ss << "-";
    ss << (var(maxsat_formula->getSoftClause(id).clause[j])+1) << " ";
  }
  ss << "0\n";
  return ss.str();
}

void MaxSAT::printUnsatisfiedSoftClauses() {
  assert (model.size() != 0);

  std::stringstream s;
  int soft_size = 0;
  
  for (int i = 0; i < maxsat_formula->nSoft(); i++) {
    bool unsatisfied = true;
    for (int j = 0; j < maxsat_formula->getSoftClause(i).clause.size(); j++) {

      assert(var(maxsat_formula->getSoftClause(i).clause[j]) <
             model.size());
      if ((sign(maxsat_formula->getSoftClause(i).clause[j]) &&
           model[var(maxsat_formula->getSoftClause(i).clause[j])] ==
               l_False) ||
          (!sign(maxsat_formula->getSoftClause(i).clause[j]) &&
           model[var(maxsat_formula->getSoftClause(i).clause[j])] ==
               l_True)) {
        unsatisfied = false;
        break;
      }
    }

    if (unsatisfied) {
      s << printSoftClause(i);
      soft_size++;
    }
  }
  FILE * file = fopen (getPrintSoftFilename(),"w");
  fprintf(file,"p cnf %d %d\n",maxsat_formula->nInitialVars(),soft_size);
  fprintf(file,"%s", s.str().c_str());
}

// Prints search statistics.
void MaxSAT::printStats() {
  double totalTime = cpuTime();
  float avgCoreSize = 0;
  if (nbCores != 0)
    avgCoreSize = (float)sumSizeCores / nbCores;

  printf("c\n");
  if (model.size() == 0)
    printf("c  Best solution:          %12s\n", "-");
  else
    printf("c  Best solution:          %12" PRIu64 "\n", ubCost);
  printf("c  Total time:             %12.2f s\n", totalTime - initialTime);
  printf("c  Nb SAT calls:           %12d\n", nbSatisfiable);
  printf("c  Nb UNSAT calls:         %12d\n", nbCores);
  printf("c  Average core size:      %12.2f\n", avgCoreSize);
  printf("c  Nb symmetry clauses:    %12d\n", nbSymmetryClauses);
  printf("c lfinal satcalls usatcalls best totaltime\n");
  printf("c final\t%d\t%d\t%lu\t%12.2f\n", nbSatisfiable, nbCores, ubCost, totalTime - initialTime);

  printf("c\n");
}

// Prints the corresponding answer.
void MaxSAT::printAnswer(int type) {
  if (verbosity > 0 && print)
    printStats();
  bool isMO = maxsat_formula->nObjFunctions() > 1;
  printf("c ---------- OUTPUT ---------------\n");

  if (isMO && type == _UNKNOWN_ && effsols.size() > 0)
      type = _SATISFIABLE_;
  else if (type == _UNKNOWN_ && model.size() > 0)
    type = _SATISFIABLE_;

  // store type in member variable
  searchStatus = (StatusCode)type;
  if(!print) return;

  switch (type) {
  case _SATISFIABLE_:
    printf("s SATISFIABLE\n");
    if(isMO){			// 
        bool printlbset = lbseti_expeps >= 1;
        printf("c satisfiable, print lb set (? %s\n", (printlbset) ? "yes" : "no");   
        if(printlbset){
            printf("c lbset of it %lu\n", lbseti_expeps);
            clearLowerBoundSet(lbseti_expeps); //remove the last (incomplete) LBset
            
            for(size_t i = 0; i < nondom.size(); i++)
                updateLowerBoundSet(nondom[i], true, lbseti_expeps);
            
        }
        printEffSolutions(printlbset);
        if(printlbset) 
            printApproxRatio();
        
    }else{
        if (print_model)
        printModel();
        if (print_soft)
        printUnsatisfiedSoftClauses();
    }
    printMyStats();
    fflush(stdout);
    break;
  case _OPTIMUM_:
    printf("s OPTIMUM\n");
    if(isMO){ printEffSolutions(true);
        printApproxRatio();
    }
    else{
        if (print_model)
        printModel();
        if (print_soft)
        printUnsatisfiedSoftClauses();
         printEffSolutions(false);
    }
    printMyStats();
    fflush(stdout);
    break;
  case _UNSATISFIABLE_:
    printf("s UNSATISFIABLE\n");
    break;
  case _UNKNOWN_:
    printf("s UNKNOWN\n");
    break;
    
  case _MEMOUT_:
    if(isMO){ printEffSolutions(false);}
    printMyStats();
    printf("s MEMOUT\n");
    fflush(stdout);
    break;
  default:
    printf("c Error: Invalid answer type.\n");
  }
  exit(0); //AG - adicionei depois usar o setrlimit para ignorar o tempo de encoding (em Alg_BLS.cc)
}

uint64_t MaxSAT::getUB() {
  // only works for partial MaxSAT currently
  Solver *solver = newSATSolver();

  vec<Lit> relaxation_vars;
  for (int i = 0; i < maxsat_formula->nSoft(); i++) {
    Lit p = mkLit(maxsat_formula->nVars() + i, false);
    relaxation_vars.push(p);
  }

  for (int i = 0; i < maxsat_formula->nVars() + maxsat_formula->nSoft(); i++)
    newSATVariable(solver);

  for (int i = 0; i < maxsat_formula->nHard(); i++)
    solver->addClause(maxsat_formula->getHardClause(i).clause);

  vec<Lit> clause;
  for (int i = 0; i < maxsat_formula->nSoft(); i++) {
    clause.clear();
    maxsat_formula->getSoftClause(i).clause.copyTo(clause);

    for (int j = 0; j < maxsat_formula->getSoftClause(i).relaxation_vars.size();
         j++)
      clause.push(maxsat_formula->getSoftClause(i).relaxation_vars[j]);

    clause.push(relaxation_vars[i]);

    solver->addClause(clause);
  }

  int limit = 1000;
  solver->setConfBudget(limit);

  vec<Lit> dummy;
  lbool res = searchSATSolver(solver, dummy);
  if (res == l_True) {
    uint64_t ub = computeCostModel(solver->model);
    return ub;
  } else if (res == l_False) {
    printAnswer(_UNSATISFIABLE_);
    exit(_UNSATISFIABLE_);
  }

  return maxsat_formula->nSoft();
}

std::pair<uint64_t, int> MaxSAT::getLB() {
  // only works for partial MaxSAT currently
  Solver *solver = newSATSolver();

  vec<Lit> relaxation_vars;
  for (int i = 0; i < maxsat_formula->nSoft(); i++) {
    Lit p = mkLit(maxsat_formula->nVars() + i, false);
    relaxation_vars.push(p);
  }

  for (int i = 0; i < maxsat_formula->nVars() + maxsat_formula->nSoft(); i++)
    newSATVariable(solver);

  for (int i = 0; i < maxsat_formula->nHard(); i++)
    solver->addClause(maxsat_formula->getHardClause(i).clause);

  vec<Lit> clause;
  for (int i = 0; i < maxsat_formula->nSoft(); i++) {
    clause.clear();
    maxsat_formula->getSoftClause(i).clause.copyTo(clause);

    clause.push(relaxation_vars[i]);

    solver->addClause(clause);
  }

  std::map<Lit, int> core; // Mapping between the assumption literal and
                           // the respective soft clause.

  for (int i = 0; i < maxsat_formula->nSoft(); i++)
    core[relaxation_vars[i]] = i;

  int limit = 1000;
  lbool res = l_False;
  uint64_t lb = 0;

  vec<bool> active;
  active.growTo(relaxation_vars.size(), false);
  vec<Lit> assumptions;
  for (int i = 0; i < relaxation_vars.size(); i++) {
    if (!active[i]) {
      assumptions.push(~relaxation_vars[i]);
    }
  }

  while (res == l_False) {
    solver->setConfBudget(limit);
    res = searchSATSolver(solver, assumptions);
    if (res == l_False) {

      for (int i = 0; i < solver->conflict.size(); i++) {
        Lit p = solver->conflict[i];
        if (core.find(p) != core.end()) {
          assert(!active[core[p]]);
          active[core[p]] = true;
        }
      }

      assumptions.clear();
      for (int i = 0; i < relaxation_vars.size(); i++) {
        if (!active[i]) {
          assumptions.push(~relaxation_vars[i]);
        }
      }
      lb++;
    }
  }

  int nb_relaxed = 0;
  for (int i = 0; i < relaxation_vars.size(); i++) {
    if (active[i])
      nb_relaxed++;
  }

  return std::make_pair(lb, nb_relaxed);
}


void MaxSAT::saveEfficientSol(const vec<lbool> &currentModel, const uint64_t * point, bool filter){
    
    double runtime = cpuTime();
    lbool * model;
    uint64_t * pt;
    int d = maxsat_formula->nObjFunctions();
    
#ifdef __DEBUG__
    printf("c time (of new solution): %.3f\n", runtime-initialTime);
#endif
    if(nondom.size() == 0){
        timestats[_time1stSol_] = runtime-initialTime;
        runstats[_nsatcalls1stSol_] = nbSatisfiable;
    }
    //TODO FIX: Na versao de enumerar todas as eficientes, a relacao entre nondom e effsols sera de 1-N
    bool isdom = false;
    if(filter){
        //estou a assumir que nenhum ponto em nondom domina o point, mas o point pode dominar pontos
        //em nondom
        for(size_t i = 0; i < nondom.size() && !isdom; i++){
            if(wdominates(nondom[i], point, d)){
                printf("c point is dominated by:\n");
                for(int j = 0; j < d; j++)
                    printf(" %lu", nondom[i][j]);
                printf("\n");
                
//                 exit(-1);
                isdom = true;
            }else if(wdominates(point, nondom[i], d)){
                delete [] nondom[i];
                delete [] effsols[i].first;
                nondom.erase(nondom.begin()+i);
                effsols.erase(effsols.begin()+i);
                i--;
            }
        }
    }
    if(!isdom){
        pt = new uint64_t[d];
        for(int i = 0; i < d; i++) 
            pt[i] = point[i];
        nondom.push_back(pt);
        
        // Only store the value of the variables that belong to the original MaxSAT formula.
        model = new lbool[maxsat_formula->nVars()];
        for (int i = 0; i < maxsat_formula->nVars(); i++){
            model[i] = currentModel[i];
        //     printf(" %d", model[i]);
        }
        //   printf("\n"); 
        
        effsols.push_back(std::pair<lbool *, uint64_t>(model, nondom.size()-1));
    }
//   printEffSolutions();
//     printf("archive size: %lu\n", nondom.size());
}

//na pratica, isto calcula o multiplicative epsilon indicator, em que o LBset e o conjunto de
//referencia e o nondom e o conjunto a ser avaliado
void MaxSAT::printApproxRatio(){
    int d = maxsat_formula->nObjFunctions();
    double realeps = 0;
    double pteps, ptepsi;
    for(size_t j = 0; j < LBset.size(); j++){
        pteps = DBL_MAX;
        for(size_t i = 0; i < nondom.size(); i++){
            ptepsi = 1;
            for(int di = 0; di < d; di++){
                //se nao for 0 (se for, entao o valor acima e' o minimo da funcao objectivo di)
                //e se o racio for maior do que dos outros di
                //e se LB*epsilon > LB+1 (para contemplar os casos em que epsilon_approx < 2)
                if(LBset[j][di] >  0 && float(nondom[i][di])/LBset[j][di] > ptepsi && LBset[j][di]*(float(nondom[i][di])/LBset[j][di]) >= LBset[j][di]+1){
                    ptepsi = float(nondom[i][di])/LBset[j][di];
                }
            }
//             printf("%f < %f ?\n", ptepsi, pteps);
            if(ptepsi < pteps) pteps = ptepsi;
        }
        if(pteps > realeps) realeps = pteps;
    }
    
    printf("c ------- \n");
    printf("c observed and expected approx ratio\n");
    printf("c tapprox <= %.4f\n", realeps); //real approximation ratio
    printf("c eapprox <= %.4f\n", expepsilon);
    printf("c ------- \n");
    repsilon = realeps;
}
//returns true if 'point' was successfully added
bool MaxSAT::updateLowerBoundSet(const uint64_t * point, bool filter, uint64_t ireencode){
    uint64_t * pt;
    
    int d = maxsat_formula->nObjFunctions();
        /*
    printf("add LB (it: %lu): ", ireencode);
    for(int i = 0; i < d; i++)
        printf(" %lu", point[i]);
    printf("\n");*/
    
    bool isdom = false;
    if(filter){
        for(size_t i = 0; i < LBset.size(); i++){
            if(LBset[i][d] == ireencode && wdominates(LBset[i], point, d)){
                isdom = true;
            }else
            if((LBset[i][d] + 2 <= ireencode) || (LBset[i][d] == ireencode && wdominates(point, LBset[i], d))){
                delete [] LBset[i];
                LBset.erase(LBset.begin()+i);
                i--;
            }
        }
    }
    if(!isdom){
        pt = new uint64_t[d+1]; //a ultima posicao guarda a iteracao de reencode em que foi adicionado
        for(int i = 0; i < d; i++) 
            pt[i] = point[i];
        pt[d] = ireencode; //augmented
        LBset.push_back(pt);
    }
    return !isdom; 
}




void MaxSAT::clearLowerBoundSet(uint64_t excepti){
    if(excepti <= 0){
        for(size_t i = 0; i < LBset.size(); i++)
            delete [] LBset[i];
        LBset.clear();
    }else{
        int d = maxsat_formula->nObjFunctions();
        for(size_t i = 0; i < LBset.size(); i++){
            if(LBset[i][d] != excepti){
                delete [] LBset[i];
                LBset.erase(LBset.begin()+i);
                i--;
            }
        }
        
    }
}




void MaxSAT::printMyStats(){
    double totalTime = cpuTime();
//     Solver *S = getSolver();
//     printf("-> %d , nondom.size(): %d\n", _nnondom_, nondom.size());
    runstats[_nsatcalls_] = nbSatisfiable;
    runstats[_ncalls_] = nbSatCalls;
    runstats[_neffsols_] = effsols.size();
    runstats[_nnondom_] = nondom.size();
//     
//     if(S != NULL){
//         printf("Solver is not null!\n");
//         printf("clauses solver, formula: %d %d\n", S->nClauses(), maxsat_formula->nHard());
//         runstats[_nencvars_] = S->nVars() - maxsat_formula->nVars();
//         runstats[_nencclauses_] = S->nClauses() - maxsat_formula->nHard();  
//     }
    runstats[_nprobvars_] = maxsat_formula->nVars();
    runstats[_nprobclauses_] = maxsat_formula->nHard();
    
    updateStats(); //update info about encoding size

    printf("clrunstats %18s %12s %12s %12s %12s %12s %18s %12s %18s %18s %8s %8s %8s\n", "nsatcalls_1stSol", "nsatcalls", "ncalls", "n_eff_sols", "n_nondom", "n_prob_vars", "n_prob_clauses", "n_enc_vars(sum)", "n_enc_clauses(sum)", "n_enc_rootvars(sum)", "n_reencodes", "rapprox", "nobj");
    printf("crunstats %12d %12d %12d %12d %12d %12d %18d %12d %18d %18d %18d %18.4f %10d\n", runstats[_nsatcalls1stSol_], runstats[_nsatcalls_], runstats[_ncalls_], runstats[_neffsols_], runstats[_nnondom_], runstats[_nprobvars_], runstats[_nprobclauses_], runstats[_nencvars_], runstats[_nencclauses_], runstats[_nencrootvars_], runstats[_nreencodes_], repsilon, maxsat_formula->nObjFunctions());
//     
//     
    timestats[_totaltime_] = totalTime - initialTime;
    printf("cltimestats %12s %12s\n", "time_1stSol", "totaltime");
    printf("ctimestats %12.2f %12.2f\n", timestats[_time1stSol_], timestats[_totaltime_]);
    
    if(print_my_output){
        FILE * file = fopen (stats_file,"w");
        fprintf(file, "# %18s %12s %12s %12s %12s %12s %18s %12s %18s %18s", "nsatcalls_1stSol", "nsatcalls", "ncalls", "n_eff_sols", "n_nondom", "n_prob_vars", "n_prob_clauses", "n_enc_vars", "n_enc_clauses", "n_enc_rootvars");
        fprintf(file, " %12s %12s %12s\n", "rapprox", "time_1stSol", "totaltime");
        
        fprintf(file, "  %18d %12d %12d %12d %12d %12d %18d %12d %18d %18d", runstats[_nsatcalls1stSol_], runstats[_nsatcalls_], runstats[_ncalls_], runstats[_neffsols_], runstats[_nnondom_], runstats[_nprobvars_], runstats[_nprobclauses_], runstats[_nencvars_], runstats[_nencclauses_], runstats[_nencrootvars_]);
        fprintf(file, " %12.4f %12.2f %12.2f\n", repsilon, timestats[_time1stSol_], timestats[_totaltime_]);
    
        fclose(file);
    }
}

const YPoint& Solution::OneSolution::yPoint() {
  if(ev){
    if(s == nullptr) throw std::runtime_error("trying to get yPoint from OneSolution, without acess to objectives: missing Solution");
    if(s->maxs == nullptr) throw std::runtime_error("trying to get yPoint from OneSolution, without acess to objectives: missing MaxSAT");
    yp = s->maxs->evalModel(m);
    ev = false;
  }
  return yp;
}

void MaxSAT::updateStats(){}

void openwbo::blockModel(Solver* sol, Model& mod){
  vec<Lit> clause;
  int i = 0;
  for(auto& at: mod){
    bool sign = at == l_True? true : false;
    Lit lit = mkLit(i, sign);
    clause.push(lit);
    ++i;
  }
  sol->addClause(clause);
}

void openwbo::modelClause(Model&& mod, vec<Lit>& clause){
  int i = 0;
  clause.clear();
  for(auto& at: mod){
    bool sign = at == l_True? false : true;
    Lit lit = mkLit(i, sign);
    clause.push(lit);
    ++i;
  }
}
void openwbo::modelClause(Model&& mod, vec<Lit>& clause, const std::vector<bool>& filter){
  clause.clear();
  for(int i = 0; i < (int) mod.size(); i++){
    auto& at = mod[i];
    bool sign = at == l_True? false : true;
    if(filter[i]){
      Lit lit = mkLit(i, sign);
      clause.push(lit);
    }
  }
}

void openwbo::modelClausePlus(Model&& mod, vec<Lit>& clause){
  int i = 0;
  clause.clear();
  for(auto& at: mod){
    bool sign = at == l_True? false : true;
    if(!sign)
      continue;
    Lit lit = mkLit(i, sign);
    clause.push(lit);
    ++i;
  }
}
// extends mod, by attributing l_False to tail.
Model openwbo::modelEmbed(const Model& mod, unsigned int nVars){
  Model modl{};
  if(nVars == mod.size()){
    modl = mod;
  }
  else{
    modl = Model{(Model::size_type)nVars};
    int i = 0;
    for(int n = mod.size(); i < n; i++ )
      modl[i] = mod[i];
    for(int n = nVars; i < n; i++)
      modl[i] = l_False;
  }
  return modl;
}

void Solution::push(const Model& m, notes_t notes){mods[id++]=std::make_pair(OneSolution{this, m}, notes);}
void Solution::push(Model&& m, notes_t notes){mods[id++]=std::make_pair(OneSolution{this, std::move(m)}, notes);}
void Solution::push(OneSolution& os, notes_t notes){mods[id++]=std::make_pair(os, notes);}
// by default makes sure the result is free of inter-dominances. check
// toggles testing the new element, filter toggles testing the set. If
// Returns true if the model gets into the solution.
bool Solution::pushSafe(const Model& m, notes_t notes, bool check_new, bool check_old ){
  OneSolution osol_a{this, m};
  if(check_old)
    for(auto it = begin(), last = end();it != last;){
      Solution::OneSolution& osol = it->second.first;
      if((osol_a.yPoint() != osol.yPoint()) && pareto::dominates(osol_a, osol)){
	it = remove(it);
	dropped++;
      }else
	++it;
    }

  if(!check_new || (check_new && !pareto::dominates((*this), osol_a))){
    push(m, notes);
    return true;
  }
  barred++;
  return false;
}

std::ostream& openwbo::operator<<(std::ostream& os, const Model& mdl){
  auto print = [&os](const lbool& lb, int i){
    if(lb == l_False) os << "-";
    os << "x" << i;
  };
  if(mdl.size()){
    {
      Model::size_type i = 1;
      print(mdl[i],i);
      for(; i < mdl.size(); i++){
	os << " ";
	print(mdl[i],i);
      }}}
  return os;

}

std::ostream& openwbo::operator<<(std::ostream& os, const YPoint& yp){
  auto print = [&os](uint32_t v){
    os << v;
  };
  if(yp.size()){
    {YPoint::size_type i = 0;
      print(yp[i++]);
      for(; i < yp.size(); i++){
	os << " ";
	print(yp[i]);
      }}}
  return os;
}
