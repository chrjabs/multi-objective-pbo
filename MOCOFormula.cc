#include "MOCOFormula.h"
#include "FormulaPB.h"
#include <memory>

namespace openwbo{

  void MOCOFormula::loadFormula(){
    of.clear();
    for(int i = 0; i < msf->nObjFunctions(); i++){
      auto ofo = msf->getObjFunction(i);
      addObjFunction(ofo);
    }
  }

  void MOCOFormula::addObjFunction(PBObjFunction *ofo){
    //     objective_function = new PBObjFunction(of->_lits, of->_coeffs, of->_const);
    of.push_back(std::make_unique<PBObjFunction>(*ofo));
  }

  std::unique_ptr<PBObjFunction>& MOCOFormula::getObjFunction(int i) { return of.at(i); } //AG
  void MOCOFormula::replaceObjFunction(int i, std::unique_ptr<PBObjFunction>&& new_pb) {
    of.at(i).swap(new_pb);
 } //AG


}
