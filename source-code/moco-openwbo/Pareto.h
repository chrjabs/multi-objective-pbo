#include <iostream>
#include "MaxSAT.h"
#include "MOCO.h"
namespace pareto{
  using namespace openwbo;
  bool dominates(Solution::OneSolution& osol_a, 
		 Solution::OneSolution& mb);
  bool dominates(const YPoint& ypa, const YPoint& ypb);
  bool dominates(Solution& sola, Solution& solb);
  bool dominates(Solution& sol, const Model& m);
  bool dominates(Solution& sol, Solution::OneSolution& osol);
  bool dominates(Solution& sol, YPoint& yp);
  openwbo::YPoint dominator(Solution& sol);

}
