#ifndef ENC_UNARYCOUNTER_H
#define ENC_UNARYCOUNTER_H

#include "RootLits.h"
#include "core/SolverTypes.h"
#include "encoder.h"
#include <string>

#define MAXDIM 10

namespace encoder {
  class UnaryCounter{
    //one RootLits instance per objective.
    rootLits::RootLits _rootLits[MAXDIM]; 
    //given variable index, returns corresponding value. Inverts _rootLits
    std::map<int, int> _invRootLits;

  public:
    UnaryCounter():_rootLits{}{}
    // first pair <key, Lit> with key larger than val, for objective iObj
    pair<uint64_t, Glucose::Lit> get(int iObj, uint64_t val);
    //human readable string, with interpretation of Lit.
    std::string notes(Lit);
  };
}

#endif
