#include "RootLits.h"
#include <algorithm>

namespace rootLits{
  void allSums(vector<uint64_t>::iterator beg,
	       vector<uint64_t>::iterator end, 
	       set<uint64_t>& result){
    if(beg == end)
      return;
    auto first = *beg;
    allSums(++beg,end,result);
    vector<uint64_t> l_res{result.begin(), result.end()};
    result.insert(first);
    for(auto& el: l_res)
      result.insert(first+el);
    return;
  }

  PBObjFunction RootLitsSliced::slice(const std::set<Lit>& vars){
    vector<uint64_t> terms{};
    vec<Lit> lits;
    vec<uint64_t> coeffs;
    PBObjFunction result{};
    // 
    for(int j = 0; j < pb._lits.size(); j++)
      if(vars.count(pb._lits[j])){
	terms.push_back(pb._coeffs[j]);
	lits.push(pb._lits[j]);
	coeffs.push(pb._coeffs[j]);
      }
    PBObjFunction inc_pb{lits,coeffs,pb._const};
    cur = std::move(*openwbo::add(&inc_pb, &cur));
    clear();
    set<uint64_t> sums{};
    sums.insert(0);
    terms.clear();
    for(int j = 0; j < cur._lits.size(); j++)
      terms.push_back(cur._coeffs[j]);
    allSums(terms.begin(), terms.end(), sums);
    for(auto& el: sums)
      me.insert(el+1, {el+1, muse.at_key(el+1)->second});
    return cur;
  }



  Iterator RootLits::at_key(uint64_t val){
    const auto pair = lower_bound(
			      myType::begin(), myType::end(),
			      val,
			      [] (myType::value_type p, uint64_t val)
			      {return p.first < val;}
			      );
    return {pair};
  }

  void RootLits::insert(uint64_t n, std::pair<uint64_t, Lit> p){
    myType::push_back(p);
}
  void combination(RootLitsInt& rootLitsOld, RootLitsInt& rootLits, map<uint64_t,Lit>& vars, Solver* solver){
    for(auto const& cur: rootLitsOld){
      for(auto const& cur_1: rootLits){
	auto val = cur.first+cur_1.first-1;
	Lit lit;
	if(vars.count(val))
	  lit = vars.at(val);
	else{
	  lit = mkLit(solver->newVar());
	  vars[val]=lit;
	}	
      }
    }
  }

  void combinationClauses(RootLitsInt& rootLitsOld, RootLitsInt& rootLits, map<uint64_t,Lit>& vars, Solver* solver){
    auto prev = *(rootLitsOld.cbegin());
    for(auto const& cur: rootLitsOld) {
      if(cur == prev)
	continue;
      auto start = vars.cbegin();
      auto prev1 = *(rootLits.cbegin());
      for(auto const& cur1: rootLitsOld){
	if(cur1 == prev1)
	  continue;
	auto val = cur.first+cur1.first-1;
	auto it = --lower_bound(start,vars.cend(),val,
				[] (const pair<uint64_t, Lit>& p, uint64_t val)
				{return p.first < val;});
	if(it != vars.end()){
	  start = it;
	  Lit lit = it->second;
	  solver->addClause(prev.second, prev.second, ~lit);
	}
	prev1 = cur1;
      }
      prev = cur;
    }
  }

void propagation(RootLitsInt& rootLits, map<uint64_t,Lit>& vars, Solver* solver){
  for(auto const& el: rootLits){
    Lit lit;
    if(vars.count(el.first))
      lit = vars.at(el.first);
    else{
      lit = mkLit(solver->newVar());
      vars[el.first]=lit;
    }
  }
}

void propagationClauses(RootLitsInt& rootLits, map<uint64_t,Lit>& vars, Solver* solver){
  for(auto const& el: vars){
    auto oldt = rootLits.at_key(el.first);
    auto newt = el.second;
    if(oldt != rootLits.end())
      solver->addClause(~newt, oldt->second);
    }
  }


  void orderEncoding(const RootLitsInt& objRootLits, Solver* solver){
    if(objRootLits.size() >= 2)
      // removing the sentinel
      for(auto prev=objRootLits.cbegin(), 
	    end =objRootLits.cend(),
	    rit=++objRootLits.cbegin();
	  rit!=end; rit++, prev++){
	solver->addClause(~prev->second, rit->second);
      }
  }







}
