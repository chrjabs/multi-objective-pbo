#include "MOCO.h"
#include "MOCOFormula.h"
#include "bounds/BoundsCalculatorMaker.h"
#include <memory>

namespace openwbo {
  
  void MOCO::loadFormula(MaxSATFormula *maxsat){
    MaxSAT::loadFormula(maxsat);
    mf = std::make_shared<openwbo::MOCOFormula>(maxsat);
    calc = bounds::maker(maxsat);
    mf->loadFormula();
  }

  YPoint MOCO::evalModel(const Model& model){
      
    YPoint yp{getFormula()->nObjFunctions()};
    if(model.size()){
      Lit l;

      for(int di = 0; di < getFormula()->nObjFunctions(); di++){
	auto& objf = getFormula()->getObjFunction(di);
	for(int i = 0; i < objf->_lits.size(); i++){
	  l = objf->_lits[i];
                
	  if ((sign(l) && model[var(l)] == l_False) || 
	      (!sign(l) && model[var(l)] == l_True)) {
	    yp[di] += objf->_coeffs[i];
	    //                 printf(" %cx%d", (sign(l)) ? '-':' ', var(l)+1);
	  }
	}
      }
    }
    return yp;
  }

  void Solution::pop(){
    mods.erase(--mods.end());
  }

  void MOCO::transferToEffSols(){
    int no = getFormula()->nObjFunctions();
    uint64_t point[getFormula()->nObjFunctions()];  
    vec<lbool> vm;
    auto & sol= *(--solution().end());
    sol.second.first.model().copyTo(vm);
    YPoint yp = sol.second.first.yPoint();
    for(int i = 0; i < no; i++ )
      point[i] = yp[i];
    saveEfficientSol(vm, point);
  }

  void MOCO::transferToSolution(){
    solution().push(Model{(*(--effsols.end())).first, (Model::size_type)getFormula()->nVars()});
  }
  std::shared_ptr<MOCOFormula> Solution::maxsat_formula(){
    return maxs->getFormula();
  }

  // Prints the corresponding answer.
  void MOCO::printAnswer(int type) {
    consolidateSolution();
    while(solution().size()){
      transferToEffSols();
      solution().pop();
    }
    if (verbosity > 0 && print)
      printStats();
    printf("c ---------- OUTPUT ---------------\n");

    if (type == _UNKNOWN_ && effsols.size() > 0)
      type = _SATISFIABLE_;
    else if (type == _UNKNOWN_ && model.size() > 0)
      type = _SATISFIABLE_;

    // store type in member variable
    searchStatus = (StatusCode)type;
    if(!print) return;

    switch (type) {
    case _SATISFIABLE_:
      printf("s SATISFIABLE\n");
      clearLowerBoundSet(lbseti_expeps); //remove the last (incomplete) LBset
      printEffSolutions(true);
      printMyStats();
      fflush(stdout);
      break;
    case _OPTIMUM_:
      printf("s OPTIMUM\n");
      printEffSolutions(true);
      printApproxRatio();
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
      printEffSolutions(false);
      printMyStats();
      printf("s MEMOUT\n");
      fflush(stdout);
      break;
    default:
      printf("c Error: Invalid answer type.\n");
    }
    exit(0); //AG - adicionei depois usar o setrlimit para ignorar o tempo de encoding (em Alg_BLS.cc)
  }

  // CJ: generates assumptions for the objective variables from a model
  void objAssumptions(const Model &mod, MaxSATFormula *formula,
                      vec<Lit> &assumpts) {
    assumpts.clear();
    for (Var v : formula->objVars()) {
      assert(v < mod.size());
      bool sign = mod[v] == l_True ? false : true;
      assumpts.push(mkLit(v, sign));
    }
  }
}
