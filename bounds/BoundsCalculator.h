#ifndef BOUNDS_H
#define BOUNDS_H
#include "../MaxSATFormula.h"
#include "core/SolverTypes.h"
#include <cstdint>
#include "BoundsCalculator_I.h"


#include <memory>
#include <vector>
using namespace std;
namespace bounds{
  class BoundsCalculatorPlain: public BoundsCalculator_I {
    
  public:
    BoundsCalculatorPlain(){};
    BoundsCalculatorPlain(openwbo::MaxSATFormula * msff):msf{msff}{
      loadFormula();
    };
    
    ~BoundsCalculatorPlain(){}
    int64_t upperBound();
    int64_t lowerBound();
    void setObjective(openwbo::PBObjFunction* pbb){pb = pbb;}
    void loadFormula();
  protected:
    openwbo::MaxSATFormula * msf;
    openwbo::PBObjFunction* pb;
    bool done = false;
  };

}

#endif
