#include "partition.h"

namespace partition{
  void logPart(typename MyPartition::part_t& p){
    
  }
  void filter(typename MyPartition::part_t& p, map<Lit, uint64_t>& terms){
    auto iter = terms.begin();
    while(iter != terms.end()){
      if(!p.count(iter->first))
	terms.erase(iter++);
      else
	      ++iter;
    }
  }

MyPartition mix(vector<MyPartition>&& partitions){
  MyPartition mix{};
  random_device rd{};
  mt19937 rng{rd()};
  uniform_int_distribution<int> uni(0,partitions.size()-1);
  int rnd{};
  while(partitions.size()){
    uniform_int_distribution<int> uni(0,partitions.size()-1);
    rnd = uni(rng);    
    mix.push(partitions[rnd].pop());

    if(partitions[rnd].size()==0){
      partitions[rnd]=partitions[partitions.size()-1];
      partitions.pop_back();
    }
  }
  return mix;
}




  std::ostream& operator<<(std::ostream& os, MyPartition& mp){
    return os;
  }
}
