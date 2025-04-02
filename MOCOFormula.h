/*!
 * \author Ruben Martins - ruben@sat.inesc-id.pt
 *
 * @section LICENSE
 *
 * Open-WBO, Copyright (c) 2013-2015, Ruben Martins, Vasco Manquinho, Ines Lynce
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 */

#ifndef MOCOFormula_h
#define MOCOFormula_h

//max number of objective functions
#include <memory>
#define MAXDIM 10

#ifdef SIMP
#include "simp/SimpSolver.h"
#else
#include "core/Solver.h"
#endif

#include "maxConsts.h"
#include "FormulaPB.h"
#include "MaxSATFormula.h"
#include "MaxTypes.h"
#include <map>
#include <vector>
#include <string>

namespace openwbo {

  class MOCOFormula {
    /*! This class contains the MaxSAT formula and methods for adding soft and
     * hard clauses. */
  private:
    MaxSATFormula* msf;
  public:
    MOCOFormula():msf{}, of{}{}
    MOCOFormula(MaxSATFormula* msff): msf{msff}, of{}{}
    // objectives are copied, formula is shared.
    MOCOFormula(const MOCOFormula& mf){
      msf = mf.msf;
      for(const auto& el: mf.of)
	of.push_back(std::make_unique<PBObjFunction>(*el));
    }
    
    void loadFormula();
    //objective functions
    std::vector<std::unique_ptr<PBObjFunction>> of{}; 
    void addObjFunction(PBObjFunction *of);
    //   PBObjFunction *getObjFunction() { return objective_function; }
    MaxSATFormula* maxsat_formula(){return msf;}
    std::unique_ptr<PBObjFunction>& getObjFunction(int i=0); //AG
    void replaceObjFunction(int i, std::unique_ptr<PBObjFunction>&& new_pb); //AG
    int nObjFunctions() { return of.size(); } //AG
    int nVars() {
      return msf->nVars();
    } // Returns the number of variables in the working MaxSAT formula.

    int nSoft() {
      return maxsat_formula()->nSoft();
    }

    Soft &getSoftClause(int pos) {
      return maxsat_formula()->getSoftClause(pos);
    }

    // Adds a new hard clause to the hard clause database.
    void addHardClause(vec<Lit> &lits) {
      maxsat_formula()->addHardClause(lits);
    }

    Hard &getHardClause(int pos) {
      return maxsat_formula()->getHardClause(pos);
    }

    int nHard() {
      return maxsat_formula()->nHard();
    } // Returns the number of hard clauses in the working MaxSAT formula.
    indexMap &getIndexToName() { return maxsat_formula()->getIndexToName(); }
    nameMap &getNameToIndex() {  return maxsat_formula()->getNameToIndex(); } 
    int nPB() { return maxsat_formula()->nPB(); }
    
    /*! Return i-PB constraint. */
    PB *getPBConstraint(int pos) {
      return maxsat_formula()->getPBConstraint(pos); }
    int nCard() { return maxsat_formula()->nCard(); }

    /*! Return i-card constraint. */
    Card *getCardinalityConstraint(int pos) {
      return maxsat_formula()->getCardinalityConstraint(pos);
    }
    int nInitialVars() {
      return maxsat_formula()->nInitialVars();
    } // Returns the number of variables in the working MaxSAT formula.


    int64_t getTighterUB(int di){
      return maxsat_formula()->bounds[di][2];
    }

    int64_t getTighterLB(int di){
      return maxsat_formula()->bounds[di][1];
    }

    void setTighterUB(int di, int64_t tub){
      maxsat_formula()->bounds[di][2] = tub;
    }

    void setTighterLB(int di, int64_t tlb){
      maxsat_formula()->bounds[di][1] = tlb;
    }
    int64_t getUB(int di){
      return maxsat_formula()->bounds[di][3];
    }

    int64_t getLB(int di){
      return maxsat_formula()->bounds[di][0];
    }

    void setUB(int di, int64_t ub){
      maxsat_formula()->bounds[di][3] = ub;
    }

    void setLB(int di, int64_t lb){
      maxsat_formula()->bounds[di][0] = lb;
    }

    // Makes a new literal to be used in the working MaxSAT formula.
    Lit newLiteral(bool sign=false) {
      return maxsat_formula()->newLiteral(sign);
    }

    void newVar(int v=-1) {
      return maxsat_formula()->newVar(v);
    } // Increases the number of variables in the working MaxSAT formula.
    void setFormat(int form) { maxsat_formula()->setFormat(form); }
    int getSumWeights(){
      return msf->getSumWeights();
    }
  };
}

// namespace openwbo

#endif
