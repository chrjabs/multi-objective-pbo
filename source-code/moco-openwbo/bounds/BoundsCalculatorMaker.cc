#include "BoundsCalculatorMaker.h"
#include "BoundsCalculator.h"
#ifdef CPLEX_BOUNDS_DEPENDENCY
#include "BoundsCalculatorCPLEX.h"
#endif


std::unique_ptr<bounds::BoundsCalculator_I> bounds::maker(openwbo::MaxSATFormula* msff){
#ifdef CPLEX_BOUNDS_DEPENDENCY
  return std::make_unique<BoundsCalculatorCPLEX>(BoundsCalculatorCPLEX{msff});
#endif
  return std::make_unique<BoundsCalculatorPlain>(BoundsCalculatorPlain{msff});
}
