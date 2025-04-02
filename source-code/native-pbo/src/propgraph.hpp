#pragma once

#include <unordered_map>
#include <unordered_set>
#include <vector>
#include "Solver.hpp"

namespace rs {

// finds sets of literals out of which at most one can be false
std::vector<Lit> find_am1s(Solver slv, std::vector<Lit> lits);

class PropGraph {
 private:
  std::unordered_map<Lit, std::unordered_set<Lit>> cons{};

 public:
  PropGraph() = default;
  PropGraph(const std::unordered_set<Lit>& lits) { add_lits(lits); }
  PropGraph(const std::vector<Lit>& lits) { add_lits(lits); }

  inline void add_lits(const std::unordered_set<Lit>& lits) {
    for (Lit l : lits)
      if (!has_lit(l)) cons.insert({l, {}});
  }

  inline void add_lits(const std::vector<Lit>& lits) {
    for (Lit l : lits)
      if (!has_lit(l)) cons.insert({l, {}});
  }

  inline bool has_lit(Lit l) const { return cons.find(l) != cons.end(); }

  inline size_t degree(Lit l) const { return cons.at(l).size(); }

  inline void add_edge(Lit l1, Lit l2) {
    assert(cons.find(l1) != cons.end());
    cons[l1].insert(l2);
    assert(cons.find(l2) != cons.end());
    cons[l2].insert(l1);
  }

  inline bool has_edge(Lit l1, Lit l2) const { return cons.at(l1).find(l2) != cons.at(l1).end(); }

  std::vector<Lit> greedy_max_clique(Lit start_lit) const {
    assert(cons.find(start_lit) != cons.end());
    std::vector<Lit> clique = {start_lit};
    for (Lit l : cons.at(start_lit)) {
      bool can_add = true;
      for (auto it = ++clique.begin(); it != clique.end(); ++it) {
        if (!has_edge(l, *it)) {
          can_add = false;
          break;
        }
      }
      if (can_add) clique.push_back(l);
    }
    return clique;
  }

  template <class Compare, class Filter>
  std::vector<Lit> greedy_max_clique(Lit start_lit, Compare comp, Filter filt) const {
    assert(cons.find(start_lit) != cons.end());
    std::vector<Lit> start_cons(cons.at(start_lit).begin(), cons.at(start_lit).end());
    sort(start_cons.begin(), start_cons.end(), comp);
    std::vector<Lit> clique = {start_lit};
    for (Lit l : start_cons) {
      bool can_add = true;
      for (auto it = ++clique.begin(); it != clique.end(); ++it) {
        if (!has_edge(l, *it)) {
          can_add = false;
          break;
        }
      }
      if (can_add && filt(l)) clique.push_back(l);
    }
    return clique;
  }
};

}  // namespace rs
