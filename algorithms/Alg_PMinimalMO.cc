
// #define PARTIAL
#ifndef PARTIAL

#include "Alg_PMinimalMO.h"	
#include <algorithm>    // std::max

using namespace openwbo;
//using namespace NSPACE;
using NSPACE::toLit;

void PMinimalMO::search_MO(int strategy){
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
      case 20: searchPMinimalMO(); break;
                
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
    if((permanentBlock && !resform) || epsilon <= 1 || redFactor < 0){
      terminate = true;
      printf("c time to terminate\n");
    }else{
      if(epsilon <= epsthreshold)
	epsilon = 1;
      else
	epsilon = 1 + (epsilon-1)/redFactor;
      printf("c REENCODE epsilon = %f\n", epsilon);
    }
  }
    
  if(nondom.size() > 0){
        
    if(epsilon <= 1){
      printf("c LBset = PF\n");
      clearLowerBoundSet();
      for(size_t i = 0; i < nondom.size(); i++)
	updateLowerBoundSet(nondom[i], false);
    }
    answerType = _OPTIMUM_;
  }else{
    if(nreencodes == 1)
      clearLowerBoundSet();
  }
 
  printAnswer(answerType);
 }


bool PMinimalMO::searchPMinimalMO() {

  double runtime = cpuTime();
  int nObj = maxsat_formula->nObjFunctions();
  
  YPoint ul (nObj);
  for(int i = 0; i < nObj; i++){
    ul[i] = 0;
  }
  for(auto sat = solve(); sat == l_True; sat = solve()) {
    for(;sat == l_True; sat = solve()) {
      Model m =Model{solver->model};
      vec<lbool> proj;
      for(int i = 0, n = maxsat_formula->nInitialVars(); i < n; i++)
	proj.push(solver->model[i]);
      Model m1 = Model{proj};
      // Only block dominated region if m1 gets into the Solution
      if(solution().pushSafe(m1)){
	ul = solution().yPoint();
	blockDominatedRegion(ul);
	printf("c o ");
	std::cout<< ul <<std::endl;
	runtime = cpuTime();
	printf("c new feasible solution (time: %.3f)\n", runtime - initialTime);
      }
      BLS::assumeDominatingRegion(ul);
    }

    runtime = cpuTime();
    printf("c new optimal solution (time: %.3f)\n", runtime - initialTime);
    blockDominatedRegion(ul);
    assumptions.clear();
  }

  if (solution().size() == 0) {
    answerType=_UNSATISFIABLE_;
    return false;
  }else{
    answerType = _OPTIMUM_;
  }
  return true;
}


/*assumes the region that dominates the point given by bojix*/
void PMinimalMO::assumeDominatingRegion(uint64_t * objix, int nObj) {
  assumptions.clear();
  for(int di = 0; di < nObj; di++){
        int j = objix[di];
        if(j > 0){
            assumptions.push((*objRootLits[di])[j].second);
#ifdef __DEBUG__
                printf("assume all %d [%lu][var: %d]\n", di, objRootLits[di][j].first, var(objRootLits[di][j].second));
#endif

        }
    }
}



#endif
#undef PARTIAL

