#include "Enc_UnaryCounter.h"

namespace encoder {
  pair<uint64_t, Glucose::Lit> UnaryCounter::get(int iObj, uint64_t val){
    return _rootLits[iObj].at_key(val);};

  std::string UnaryCounter::notes(Lit){
    return "this is the meaning of Lit";
  }
}
