#include "FormulaPB.h"

namespace openwbo{
std::unique_ptr<PBObjFunction>  add(PBObjFunction* pb, PBObjFunction* pbb){
  vec<Lit> lits{};
  vec<uint64_t> coeffs{};
  std::map<Lit, uint64_t> terms{};  
  
  if(pb!=nullptr)
    for(int j = 0; j < pb->_lits.size(); j++)
      terms[pb->_lits[j]]+=pb->_coeffs[j];
  if(pbb!=nullptr)
    for(int j = 0; j < pbb->_lits.size(); j++)
      terms[pbb->_lits[j]]+=pbb->_coeffs[j];
  for(auto& el: terms){
    lits.push(el.first);
    coeffs.push(el.second);
  }

  return std::make_unique<PBObjFunction>(PBObjFunction{lits,coeffs,pb->_const + pbb->_const});
}



}

