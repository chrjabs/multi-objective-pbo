#ifndef BOUNDSCALCULATORCPLEX_H
#define BOUNDSCALCULATORCPLEX_H

#include "ilconcert/ilomodel.h"
#include "ilconcert/ilosmodel.h"
#include "ilcplex/cplex.h"
#include <ilcplex/ilocplexi.h> 
#include <ilconcert/ilolinear.h> 
#include "ilconcert/iloexpression.h"
#include "ilconcert/ilolinear.h"
#include "ilconcert/ilosys.h"
#include "BoundsCalculator_I.h"
#include "BoundsCalculator.h"
namespace bounds{
class BoundsCalculatorCPLEX: public bounds::BoundsCalculator_I {
    
  public:
    BoundsCalculatorCPLEX(){};
    BoundsCalculatorCPLEX(openwbo::MaxSATFormula * msff):msf{msff}{
      model = IloModel(env);
      var = IloNumVarArray(env,msf->nVars(),0,1);
      con = IloRangeArray(env);
      criteria = IloNumExprArray(env);
      cplex = IloCplex(model);
      obj = IloObjective(env);
      loadFormula();
    };
    
    ~BoundsCalculatorCPLEX(){}
    int64_t upperBound();
    int64_t lowerBound();
    void buildObj();
    void setObjective(openwbo::PBObjFunction* pbb){pb = pbb;buildObj();}
    void loadFormula();
  protected:
    openwbo::MaxSATFormula * msf;
    openwbo::PBObjFunction* pb;
    IloEnv env;
    IloModel model;
    IloObjective obj;
    IloNumVarArray var;
    IloRangeArray con;
    IloNumExprArray criteria;
    IloCplex cplex;
    bool done = false;
  };


  IloConstraint ilo_cons(openwbo::Hard& h, IloNumVarArray var, IloEnv& env);
  IloConstraint ilo_cons(openwbo::PB* pb, IloNumVarArray var, IloEnv& env);
  IloConstraint ilo_cons(openwbo::Card* c, IloNumVarArray var, IloEnv& env);


}
#endif
