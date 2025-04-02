#ifndef ALG_MASTERMO_H
#define ALG_MASTERMO_H
#include "Alg_ServerMO.h"
#include "../partition.h"
#include <vector>

namespace openwbo{
  class MasterMO{
  public:
    BLSServerMO* optim;
  public:
    virtual void initializeOptimizer(Solver*, MaxSATFormula*) = 0;
    //issues solve on the optim auxiliar solver
    virtual StatusCode compute_approx(){return optim->searchAgain();}
    // changes data according to the results of compute_approx, stored
    // in solution
    virtual bool incorporate_approx()=0;
    //get optim ready for the next solve call. If there is no need to
    //go on, return false.
    virtual bool setup_approx()=0;
  };

  class  BLSMasterMO: public virtual BLS , public virtual MasterMO{
  public:
    bool buildWorkFormula() override;
    virtual void build() override =0;
    void search_MO(int) override{
      build();
      if(firstSolution()){
	buildWorkFormula();
	searchMasterMO();
	consolidateSolution();
	if(solution().size() > 0){
	  answerType = _OPTIMUM_;
	}else{
	  answerType = _UNSATISFIABLE_;
	}
	printAnswer(answerType);
      }
    }

    StatusCode searchMasterMO(){
      do{
	compute_approx();
	incorporate_approx();
      }while(setup_approx());
      if (solution().size() == 0) {
	answerType=_UNSATISFIABLE_;
      }else{
	answerType = _OPTIMUM_;
      }
      return answerType;
    }
  };

}

#endif
