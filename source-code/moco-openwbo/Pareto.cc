#include "Pareto.h"
// point that tightly dominates sol
openwbo::YPoint pareto::dominator(Solution& sol){
  int nObj = sol.maxsat_formula()->nObjFunctions();
  YPoint yp(nObj);
  if(sol.size() == 0)
    return yp;
  yp = (*sol.begin()).second.first.yPoint();
  for(auto& os: sol){
    auto& yp1 = os.second.first.yPoint();
    for(int i = 0; i < nObj; i++)
      if(yp[i] > yp1[i])
	yp[i] = yp[i];
  }
  return yp;
}

bool pareto::dominates(Solution::OneSolution &osol_a, 
		       Solution::OneSolution &osol_b){
  if(osol_a.comparable(osol_b))
    return dominates(osol_a.yPoint(), osol_b.yPoint());
  throw std::runtime_error("cannot compare incomparable OneSolution pair");
}


bool pareto::dominates(const YPoint &ypa, const YPoint &ypb){
  assert(ypa.size() > 0 && ypa.size() == ypb.size());
  for(int i = 0, n = ypa.size(); i < n; i++)
    if(ypa[i] > ypb[i])
      return false;
  return true;
}

bool pareto::dominates(Solution &sol, Solution::OneSolution& osol){
  //no const iterator, evaluating yPoint changes OneSolution
  if(sol.size() == 0)
    return false;
  try{
    for(auto& el: sol)
      if(dominates(el.second.first,osol))
	return true;
  } catch(std::runtime_error& e) {std::cout << e.what();}
  return false;
}
bool pareto::dominates(Solution &sol_a, Solution &sol_b){
  //no const iterator, evaluating yPoint changes OneSolution
  auto it_b = sol_b.begin();
  try{
    while(dominates(sol_a, (it_b++)->second.first));
  if(it_b == sol_b.cend())
    return true;
  } catch(std::runtime_error& e) {std::cout << e.what();}
  return false;
}

bool pareto::dominates(Solution &sol_a, const Model& m){
  Solution::OneSolution osol = sol_a.wrap(m);
  return dominates(sol_a, osol);
}


