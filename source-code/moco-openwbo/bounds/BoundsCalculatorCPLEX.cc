#include "BoundsCalculatorCPLEX.h"
#include "ilconcert/iloexpression.h"
#include "ilconcert/ilolinear.h"
#include "ilconcert/ilomodel.h"
#include "ilconcert/ilosys.h"

void bounds::BoundsCalculatorCPLEX::buildObj(){
  const auto& lits = pb->_lits;
  const auto& coeffs = pb->_coeffs;
  const auto& _const = pb->_const;
  IloExpr exp{env};
  auto n = lits.size();
  for(IloInt i = 0; i < n; i++){
    if(!Glucose::sign(lits[i]))
      exp += IloInt(coeffs[i]) * (var[Glucose::var(lits[i])]);
    else
      exp += IloInt(coeffs[i]) * (1 - var[Glucose::var(lits[i])]);
  }
  obj= IloObjective(env,exp + IloInt(_const));
}
  
int64_t bounds::BoundsCalculatorCPLEX::lowerBound(){
  model.add(obj);
  obj.setSense(IloObjective::Minimize);
  if(!cplex.solve()){
    env.error()<<"Failed to optimize LP"<<endl;
    throw(-1);}
  auto res = cplex.getObjValue();
  model.remove(obj);
  return res;
}

int64_t bounds::BoundsCalculatorCPLEX::upperBound(){
  model.add(obj);
  obj.setSense(IloObjective::Maximize);
  if(!cplex.solve()){
    env.error()<<"Failed to optimize LP"<<endl;
    throw(-1);}
  auto res = cplex.getObjValue();
  model.remove(obj);
  return res;
}
void bounds::BoundsCalculatorCPLEX::loadFormula(){
  for(int i = 0; i < msf->nHard(); i++)
    model.add(ilo_cons(msf->getHardClause(i), var, env));
  for(int i = 0; i < msf->nPB(); i++)
    model.add(ilo_cons(msf->getPBConstraint(i), var, env));
  for(int i = 0; i < msf->nCard(); i++)
    model.add(ilo_cons(msf->getCardinalityConstraint(i), var, env));
}

IloConstraint bounds::ilo_cons(openwbo::Hard& h, IloNumVarArray var, IloEnv& env ){
  const auto& lits = h.clause;
  IloExpr exp{env};
  if(lits.size()){
    for(int j = 0; j < lits.size(); j++){
      if(!Glucose::sign(lits[j]))
	exp += var[Glucose::var(lits[j])];
      else{
	exp += 1 - var[Glucose::var(lits[j])];
      }
    }
  }
  return IloConstraint(exp >= 1);
}

IloConstraint bounds::ilo_cons(openwbo::PB* pb, IloNumVarArray var, IloEnv& env){
  const auto& lits = pb->_lits;
  const auto& coeffs = pb->_coeffs;
  const auto& rhs = pb->_rhs;
  IloExpr exp{env};
  if(lits.size()){
    for(int j = 0; j < lits.size(); j++){
      if(!Glucose::sign(lits[j]))
	exp += IloInt(coeffs[j]) * (var[Glucose::var(lits[j])]);
      else
	exp += IloInt(coeffs[j]) * (1 - var[Glucose::var(lits[j])]);
    }
  }
  if(pb->_sign)
    return IloConstraint(exp <= rhs);
  else
    return IloConstraint(exp >= rhs);	
}

IloConstraint bounds::ilo_cons(openwbo::Card* c, IloNumVarArray var, IloEnv& env){
  const auto& lits = c->_lits;
  const auto& rhs = c->_rhs;
  IloExpr exp{env};
  if(lits.size()){
    for(int j = 0; j < lits.size(); j++){
      if(!Glucose::sign(lits[j]))
	exp += var[Glucose::var(lits[j])];
      else
	exp += 1 - var[Glucose::var(lits[j])];
    }
  }
  return IloConstraint(exp <= rhs);
}
