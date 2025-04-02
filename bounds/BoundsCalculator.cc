#include "BoundsCalculator.h"
#include "BoundsCalculator_I.h"
#include "core/SolverTypes.h"
#include <cstdint>
#include <memory>


int64_t bounds::BoundsCalculatorPlain::lowerBound(){
  return 0 + pb->_const;
}

int64_t bounds::BoundsCalculatorPlain::upperBound(){
  return pb->total + pb->_const;
}
void bounds::BoundsCalculatorPlain::loadFormula(){
}





