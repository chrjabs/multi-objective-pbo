#ifndef MYPARTITION
#define MYPARTITION


#include "MaxSAT.h"
#include <functional>
#include <map>
#include <set>
#include <iostream>
#include <random>

using namespace std;
using namespace openwbo;

namespace partition{
  class MyPartition;
  MyPartition partition(PBObjFunction pb);
  class MyPartition{
    std::vector<Lit> _data;
    // upper and lower indices delimiting each part. Start with 0,
    // ends with length of data.
    std::vector<int> _breaks{0};
    int param = 15;
  public:
    friend std::ostream& operator<<(std::ostream&, MyPartition& mp);
    using part_t= set<Lit>;
    MyPartition(): _data{}{}
    MyPartition(PBObjFunction* pb, int paramm){
      param = paramm;
      multimap<uint64_t, Lit,greater<uint64_t>> lits{};
      part_t prt{};
      for(int i = 0; i < pb->_lits.size(); i++) 
	lits.insert({pb->_coeffs[i], pb->_lits[i]});
      auto it = lits.begin();
      auto end = lits.end();
      int w_n = 0, n = 0;
      uint64_t last{};
      while(it!=end){
	last = it->first;
	while(it!=end && last == it->first){
	  prt.insert(it->second);
	  it++;
	  n++;
	}
	w_n++;
	if(n / w_n >= param){
	  w_n=0; n = 0;
	  push(prt);
	  prt.clear();
	}
      }
      if(prt.size()!=0)
	push(prt);
    }
    part_t part(part_t::size_type i) const{
      if(i + 1 > _breaks.size())
	throw std::out_of_range("there is no such part in this partition");
      int ux = _breaks[_breaks.size()-1];
      int lx = _breaks[_breaks.size() - 2];
      part_t pt{};
      for(int i = ux-1; i >= lx; i--){
	pt.insert(_data[i]);
      }
      return pt;}
    void push(const part_t pt){
      _data.insert(_data.end(),pt.begin(), pt.end());
      _breaks.push_back(_data.size());
    }
    part_t pop(){
      part_t pt{};
      if(_breaks.size()<=1)
	return pt;
      int ux = _breaks[_breaks.size()-1];
      int lx = _breaks[_breaks.size() - 2];

      for(int i = ux-1; i >= lx; i--){
	pt.insert(_data[i]);
	_data.pop_back();
      }
      _breaks.pop_back();
      return pt;
    }
    int size() const {return _breaks.size()-1;}
    void reverse() {
      std::reverse(_data.begin(), _data.end());
      std::reverse(_breaks.begin(), _breaks.end());      
      uint64_t len = _breaks[0];
      for (auto& el: _breaks)
	el  = len - el;
}
    void collapse(){
      _breaks.clear();
      _breaks.push_back(0);
      _breaks.push_back(_data.size());
}
  };

  //  keep only terms that intersect part p
  void filter(typename MyPartition::part_t& p, map<Lit, uint64_t>& terms);
  void logPart(typename MyPartition::part_t& p);
  MyPartition mix(vector<MyPartition>&& partitions);




}

#endif
