#ifndef BOUNDSCALCULATOR_I_H
#define BOUNDSCALCULATOR_I_H
#include <cstdint>
#include "../FormulaPB.h"

namespace bounds {
  class BoundsCalculator_I{
  public:
    int64_t virtual upperBound() = 0;
    int64_t virtual lowerBound() = 0;
    void virtual setObjective(openwbo::PBObjFunction* pbb)= 0;
    void virtual loadFormula()= 0;
  };
}


#endif
