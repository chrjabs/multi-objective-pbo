#ifndef ROOTLITS_h
#define ROOTLITS_h
#include "../partition.h"
#include "Enc_KPA.h"
#include "core/SolverTypes.h"
#include <algorithm>
#include <cstdint>
#include <iterator>
#include "RootLitsTypes.h"
#include "RootLitsInt.h"
#include "RootLitsIterator.h"
#ifdef SIMP
#include "simp/SimpSolver.h"
#else
#include "core/Solver.h"
#include <memory>
#endif

#include <vector>

using NSPACE::Lit;
using namespace std;
using namespace openwbo;

namespace rootLits{

  class RootLits: public myType, public RootLitsInt{
  
  public: 
    Iterator at_key(uint64_t) override;
    value_t at(uint64_t i) override {return myType::at(i);};
    value_t operator[](uint64_t i) override{return vector::operator[](i);};
    Iterator begin() override{return Iterator{++myType::begin()};}
    Iterator end() override{return Iterator{myType::end()};}
    CIterator cbegin() const override{return CIterator{++myType::cbegin()};}
    CIterator cend() const override{return CIterator{myType::cend()};}
    void insert(uint64_t n, value_t) override;
    void reserve(uint64_t n){myType::reserve(n + 1);}
    virtual void clear() override{myType::clear(); push_back({0, Glucose::toLit(0)});}//sentinel
    template<typename iterator1> void insert(iterator1 first, iterator1 last){
      myType::insert(++myType::begin(), first, last);
    }
    uint64_t size() const override {return myType::size();};//sentinel

    RootLits(){
      myType{};
      push_back({0, Glucose::toLit(0)});//sentinel
    }
  };
  
  class RootLitsSliced: public RootLitsInt{
  public:
    RootLitsSliced(RootLits&& _muse, PBObjFunction _pb): muse{_muse}, me{}, ref{&muse},pb{std::move(_pb)},cur{}{
    }
    RootLitsSliced(RootLitsSliced&& other):
      muse{std::move(other.muse)}, me{std::move(other.me)},
      // ref must be updated after moving, not just copied
      ref{&muse},
      pb{std::move(other.pb)},cur{std::move(other.cur)}
    {
 other.ref = nullptr;
}
       PBObjFunction slice(const std::set<Lit>& vars);
    bool setSlice(bool b){
      auto tmp = sliced;
      if(!(sliced == b))
	toggleSlice();
      return tmp;
}
    void toggleSlice(){
      if(sliced)
	ref = &muse;
      else
	ref = &me;
      sliced = !sliced;
    }
    Iterator at_key(uint64_t val) override{return (*ref).at_key(val);}
    value_t at(uint64_t i) override {return (*ref).at(i);}
    value_t operator[](uint64_t i) override{return (*ref)[i];}
    Iterator begin() override{return (*ref).begin();}
    Iterator end() override{return (*ref).end();}
    CIterator cbegin() const override{return (*ref).cbegin();}
    CIterator cend() const override{return (*ref).cend();}
    void clear() override{me.clear();}
    void insert(uint64_t n, value_t val) override{me.insert(n,val);}
    uint64_t size() const override{return (*ref).size();}
  private:
    // stores the completely encoded rootLits
    RootLits muse;
    RootLits me;
    RootLits* ref;
    PBObjFunction pb;
    bool sliced = false;
  public:
    PBObjFunction cur;

  };

  void allSums(vector<uint64_t>::iterator beg, vector<uint64_t>::iterator end, set<uint64_t>& result);
  void combination(RootLitsInt& rootLitsOld, RootLitsInt& rootLits, map<uint64_t,Lit>& vars, Solver* solver);
  void combinationClauses(RootLitsInt& rootLitsOld, RootLitsInt& rootLits, map<uint64_t,Lit>& vars, Solver* solver);
  void propagation(RootLitsInt& rootLits, map<uint64_t,Lit>& vars, Solver* solver);
  void propagationClauses(RootLitsInt& rootLits, map<uint64_t,Lit>& vars, Solver* solver);
  void orderEncoding(const RootLitsInt& objRootLits, Solver* solver);
}


#endif
