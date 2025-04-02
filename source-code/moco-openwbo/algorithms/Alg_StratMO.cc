#include "Alg_StratMO.h"

void StratMO::incrementEncoding(partition::MyPartition::part_t& p, int i){
  auto& entry = objRootLits[i];
  auto& pb = getFormula()->getObjFunction(i);

    map<Lit, uint64_t> old_terms{};  
    map<Lit, uint64_t> terms{};
    vec<Lit> lits{};
    vec<uint64_t> coeffs{};
    //get terms from complete objective function
    for(int j = 0; j < pb->_lits.size(); j++)
      if(p.count(pb->_lits[j]))
	terms[pb->_lits[j]]+=pb->_coeffs[j];

    auto old_pb = optim->getFormula()->getObjFunction(i).get();
    //get terms from the partial objective function that is going to
    //be fully updated at the end of this function
    for(int j = 0; j < old_pb->_lits.size(); j++)
      if(p.count(old_pb->_lits[j]))
	old_terms[old_pb->_lits[j]]+=old_pb->_coeffs[j];

    for(auto& el: terms){
      auto value  = el.second - old_terms[el.first];
      if(value){
	lits.push(el.first);
	coeffs.push(el.second);
      }
    }
    auto ith_orl = rootLits::RootLits{};
    for(const auto& el: *entry)
      ith_orl.insert(el.first, el);
    //only increment if the increment is not zero. This will be true almost every time.
    if(lits.size())
      {
	auto inc_pb = PBObjFunction(lits,coeffs,pb->_const);
	auto tmp = partialEncoding(i, inc_pb);
	std::map<uint64_t,Lit> vars;
	rootLits::combination(ith_orl, tmp.first, vars, solver);
	rootLits::propagation(ith_orl, vars, solver);
	rootLits::propagation(tmp.first, vars, solver);
	rootLits::propagationClauses(ith_orl, vars, solver);
	rootLits::propagationClauses(tmp.first, vars, solver);
	rootLits::combinationClauses(ith_orl, tmp.first, vars, solver);
	ith_orl.clear();
	ith_orl.insert(vars.begin(), vars.end());
	rootLits::orderEncoding(ith_orl, solver);
	// correct invObjRootLits
	// remove old values,
	for(auto it = invObjRootLits.begin(), end = invObjRootLits.begin();
	    it!=end; )
	  if(it->second == i)
	    invObjRootLits.erase(it++);
	  else
	    ++it;
	//add new values
	for(auto const& w_var: ith_orl)
	  invObjRootLits[var(w_var.second)]=i;
	// only required if the same variable appears in different parts
	// of the partition:
	auto new_pb = openwbo::add(old_pb, &inc_pb);
	optim->getFormula()->replaceObjFunction(i, std::move(new_pb));
	entry = std::make_unique<rootLits::RootLits>(std::move(ith_orl));
      }
}

void StratMO::incrementEncoding(partition::MyPartition::part_t& p){
  int nObj = maxsat_formula->nObjFunctions();
  for(int i = 0; i < nObj; i++){
    incrementEncoding(p, i);
  }
  optim->copyObjRootLits(invObjRootLits, objRootLits, nObj);
}




  std::pair<rootLits::RootLits, BLS::invRootLits_t> StratMO::partialEncoding(int iObj, PBObjFunction& pb){
    // check BLS::updateMOEncoding, if following code is confusing
    rootLits::RootLits rootLits{};
    BLS::invRootLits_t invRootLits;

    KPA kpa{};
    kpa.fixed_vars(getMaxSATFormula()->fixed_vars());
    calc->setObjective(&pb);
    auto ub = getTighterUB(iObj);
    auto ubl = calc->upperBound();
    if(ubl < ub) ub = ubl;
    kpa.encode(solver, pb._lits, pb._coeffs, ub);
    auto rootLits_1 = kpa.getRootLits();
    {
      for (const auto& el: rootLits_1){
	rootLits.push_back({el.first+1, el.second});
	invObjRootLits.insert({var(el.second), iObj});
      }
    }
    return make_pair(std::move(rootLits), std::move(invRootLits));
  }
