#ifndef BOUNDSCALCULATORMAKER_H
#define BOUNDSCALCULATORMAKER_H
#include "BoundsCalculator_I.h"
#include "../MaxSATFormula.h"
namespace bounds{
std::unique_ptr<BoundsCalculator_I> maker(openwbo::MaxSATFormula* msff);
}


#endif
