/*!
 * \author Ruben Martins - ruben@sat.inesc-id.pt
 *
 * @section LICENSE
 *
 * MiniSat,  Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson
 *           Copyright (c) 2007-2010, Niklas Sorensson
 * Open-WBO, Copyright (c) 2013-2017, Ruben Martins, Vasco Manquinho, Ines Lynce
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

#ifndef MaxSAT_h
#define MaxSAT_h

#ifdef SIMP
#include "simp/SimpSolver.h"
#else
#include "core/Solver.h"
#endif
#include "stdexcept"

#include "MaxSATFormula.h"
#include "MaxTypes.h"
#include "utils/System.h"
#include <algorithm>
#include <map>
#include <set>
#include <utility>
#include <vector>

using NSPACE::vec;
using NSPACE::Lit;
using NSPACE::lit_Undef;
using NSPACE::mkLit;
using NSPACE::lbool;
using NSPACE::Solver;
using NSPACE::cpuTime;

namespace openwbo {
  // generalized MaxSAT instance. In particular, it might actually be
  // a MOCO instance
  class MaxSAT {
    
  public:
    
    MaxSAT(MaxSATFormula *mx) {
      maxsat_formula = mx;
      searchStatus = _UNKNOWN_;
      
      // 'ubCost' will be set to the sum of the weights of soft clauses
      //  during the parsing of the MaxSAT formula.
      ubCost = 0;
      lbCost = 0;
      
      off_set = 0;
      // Statistics
      nbSymmetryClauses = 0;
      nbCores = 0;
      nbSatisfiable = 0;
      nbSatCalls = 0;
      sumSizeCores = 0;

      print_model = false;
      print_soft = false;
      print = false;
      unsat_soft_file = NULL;
    
      print_my_output = false;
      objv_file = effsols_file = stats_file = info_file = lbset_file = NULL;
    
      expepsilon = -1;
      lbseti_expeps = 0;
    }

    MaxSAT() {
      maxsat_formula = NULL;
      searchStatus = _UNKNOWN_;

      // 'ubCost' will be set to the sum of the weights of soft clauses
      //  during the parsing of the MaxSAT formula.
      ubCost = 0;
      lbCost = 0;

      off_set = 0;

      // Statistics
      nbSymmetryClauses = 0;
      nbCores = 0;
      nbSatisfiable = 0;
      nbSatCalls = 0;
      sumSizeCores = 0;

      print_model = false;
      print_soft = false;
      print = false;
      unsat_soft_file = NULL;
    
      print_my_output = false;
      objv_file = effsols_file = stats_file = info_file = lbset_file = NULL;
    
    }

    virtual ~MaxSAT() {

      for(size_t i = 0; i < effsols.size(); i++)
        delete [] effsols[i].first;
    
      for(size_t i = 0; i < nondom.size(); i++)
        delete [] nondom[i];
        
      for(size_t i = 0; i < LBset.size(); i++)
        delete [] LBset[i];
    
      effsols.clear();
      nondom.clear();
      LBset.clear();
    }
    void setInitialTime(double initial); // Set initial time.

    // Print configuration of the MaxSAT solver.
    // virtual void printConfiguration();
    void printConfiguration();

    // Encoding information.
    void print_AMO_configuration(int encoding);
    void print_PB_configuration(int encoding);
    void print_Card_configuration(int encoding);

    // Incremental information.
    void print_Incremental_configuration(int incremental);

    virtual StatusCode search();      // MaxSAT search.
    virtual void printAnswer(int type); // Print the answer.
    virtual void consolidateSolution(){};
    // Tests if a MaxSAT formula has a lexicographical optimization criterion.
    bool isBMO(bool cache = true);

    virtual void loadFormula(MaxSATFormula *maxsat) {
      maxsat_formula = maxsat;

      if (_FORMAT_PB_ == maxsat_formula->getFormat()) {
        printf("c [MaxSAT.h: loadFormula] maxsat_formula->convertPBtoMaxSAT()\n");
        off_set = maxsat_formula->getObjFunction()->_const;
        //e aqui que o a funcao obj PB e usada para definit as soft clauses
        maxsat_formula->convertPBtoMaxSAT();
      } else if (_FORMAT_MCNF_ == maxsat_formula->getFormat()) {
        printf("c [MaxSAT.h: loadFormula] converting from MCNF\n");
        maxsat_formula->clausalToPbObjs();
        
        // Debugging
        // maxsat_formula->my_print();
      }
      maxsat_formula->setInitialVars(maxsat_formula->nVars());

      ubCost = maxsat_formula->getSumWeights();
    
      printf("c [MaxSAT.h: loadFormula] ubCost: %ld\n", ubCost);
    }

    void blockModel(Solver *solver);

    // Get bounds methods
    uint64_t getUB();
    std::pair<uint64_t, int> getLB();

    Soft &getSoftClause(int i) { return maxsat_formula->getSoftClause(i); }
    Hard &getHardClause(int i) { return maxsat_formula->getHardClause(i); }
    Lit getAssumptionLit(int soft) {
      return maxsat_formula->getSoftClause(soft).assumption_var;
    }
    Lit getRelaxationLit(int soft, int i = 0) {
      return maxsat_formula->getSoftClause(soft).relaxation_vars[i];
    }

    int64_t getOffSet() { return off_set; }

    MaxSATFormula *getMaxSATFormula() { return maxsat_formula; }

    void setPrintModel(bool model) { print_model = model; }
    bool getPrintModel() { return print_model; }

    void setPrint(bool doPrint) { print = doPrint; }
    bool getPrint() { return print; }

    void setPrintSoft(const char* file) { 
      if (file != NULL){
	unsat_soft_file = (char*)malloc(sizeof(char) * (sizeof(file)));
	strcpy(unsat_soft_file,file) ; 
	print_soft = true;
      }
    }
    bool isPrintSoft() { return print_soft; }
    char * getPrintSoftFilename() { return unsat_soft_file; }
  
  
    char * getFilename(const char* file, const char *ext){
      char * fname = (char*)malloc(sizeof(char) * (strlen(file)+strlen(ext)+1)); //+sizeof(ext)+1));
      strcpy(fname,file);
      //       int sz = strlen(file);
      //       printf("size fname: %d\n", sz);
      //       strcpy(fname + sz, ext);
      strcat(fname, ext);
      return fname;
    }
  
    void setMyOutputFiles(const char* file) { 
      printf("c set my output files\n");
      if (file != NULL){
	//solutions' image in objective space
	objv_file = getFilename(file, ".objv");
	printf("c objv_file saved to: %s\n", objv_file);
	//solutions (in decision space)
	effsols_file = getFilename(file, ".effs");
	printf("c effsols_file saved to: %s\n", effsols_file);
	//statistics (ncalls, runtime, nclauses, ...
	stats_file = getFilename(file, ".stats");
	printf("c stats_file saved to: %s\n", stats_file);
	//lower bound set (if epsilon > 1 and certificate is OPT)
	lbset_file = getFilename(file, ".lbset");
	printf("c lbset_file saved to: %s\n", lbset_file);
	//Pareto front (if epsilon == 1 and certificate is OPT)
	//       pf_file = getFilename(file, ".pf");
	//       printf("c pf_file saved to: %s\n", pf_file);
	//       info_file = getFilename(file, ".objv");
	//       printf("info_file saved to: %s\n", info_file);
	print_my_output = true;
      }
      //     printf("c done setting up output files\n");
    }
    /** return status of current search
     *
     *  This method helps to extract the status in case the solver is used as a
     *  library object without printing solutions.
     */
    StatusCode getStatus() { return searchStatus; }

    /** return truth values for variables
     *
     *  This method returns the truth value for a variable in the internal
     *  format. However, the return value reflects the external format, e.g.
     *  getValue(0) will return 1 or -1, depending on the sign of the variable
     *  in the model.
     */
    int getValue(const NSPACE::Var v)
    {
      if (v > model.size()) return 0;
      if (model[v] == l_True) return v+1;
      return -(int)v - 1;
    }

    //MO stuff
    void printEffSolutions(bool printLBset);
    void saveEfficientSol(const vec<lbool> &currentModel, const uint64_t * point, bool filter=false);
    bool updateLowerBoundSet(const uint64_t * point, bool filter, uint64_t ireencode=1);
    void clearLowerBoundSet(uint64_t excepti=0);
    void updateExpectedApprox(double eps, uint64_t lbi){printf("c updateExpectedApprox\n"); expepsilon = eps; lbseti_expeps = lbi;};
    void printApproxRatio();
    void printMyStats();
  
  protected:
    // Interface with the SAT solver
    //
    Solver *newSATSolver(); // Creates a SAT solver.
    // Solves the formula that is currently loaded in the SAT solver.
    lbool searchSATSolver(Solver *S, vec<Lit> &assumptions, bool pre = false);
    lbool searchSATSolver(Solver *S, bool pre = false);

    void newSATVariable(Solver *S); // Creates a new variable in the SAT solver.

    void reserveSATVariables(Solver *S, unsigned maxVariable); // Reserve space for multiple variables in the SAT solver.

    // Properties of the MaxSAT formula
    //
    vec<lbool> model; // Stores the best satisfying model.
    StatusCode searchStatus; // Stores the current state of the formula

    // Statistics
    //
    int nbCores;           // Number of cores.
    int nbSymmetryClauses; // Number of symmetry clauses.
    uint64_t sumSizeCores; // Sum of the sizes of cores.
    int nbSatisfiable;     // Number of satisfiable calls.
    int nbSatCalls; // Number of calls to sat solver

    // Bound values
    //
    uint64_t ubCost; // Upper bound value.
    uint64_t lbCost; // Lower bound value.
    int64_t off_set; // Offset of the objective function for PB solving.

    MaxSATFormula *maxsat_formula;

    // Others
    // int currentWeight;  // Initialized to the maximum weight of soft clauses.
    double initialTime; // Initial time.
    int verbosity;      // Controls the verbosity of the solver.
    bool print_model;   // Controls if the model is printed at the end.
    bool print;         // Controls if data should be printed at all
    bool print_soft;    // Controls if the unsatified soft clauses are printed at the end.
    char * unsat_soft_file;  // Name of the file where the unsatisfied soft clauses will be printed.
  
    char * objv_file;
    char * effsols_file;
    char * stats_file;
    char * lbset_file;
    //   char * pf_file;
    char * info_file; //not used

    // Different weights that corresponds to each function in the BMO algorithm.
    std::vector<uint64_t> orderWeights;

    // Utils for model management
    //
    void saveModel(vec<lbool> &currentModel); // Saves a Model.
    // Compute the cost of a model.
    uint64_t computeCostModel(vec<lbool> &currentModel,
			      uint64_t weight = UINT64_MAX);

    // Utils for printing
    //
    void printBound(int64_t bound); // Print the current bound.
    void printModel(); // Print the best satisfying model.
    void printStats(); // Print search statistics.
    std::string printSoftClause(int id); // Prints a soft clause.
    void printUnsatisfiedSoftClauses(); // Prints unsatisfied soft clauses.

    // Greater than comparator.
    bool static greaterThan(uint64_t i, uint64_t j) { return (i > j); }
  
    //MO stuff
    //objvalues
    std::vector<uint64_t *> nondom;
    std::vector<uint64_t *> LBset; //para as versoes approximadas
    //model (sol space) and index of nondom point to the corresponding point (in obj space)
    std::vector<std::pair<lbool *, uint64_t>> effsols;
    //MO MO stats
    bool print_my_output;
  


    //   std::string statslabels[8] = {"nsatcalls_1stSol", "nsatcalls", "n_eff_sols", "n_nondom", "time_1stSol", "totaltime", "n_enc_vars", "n_enc_clauses"};
    enum {_nsatcalls1stSol_=0, _nsatcalls_, _ncalls_, _neffsols_, _nnondom_, _nprobvars_, _nprobclauses_, _nencvars_, _nencclauses_, _nencrootvars_, _nreencodes_};
    enum {_time1stSol_, _totaltime_};
    int runstats[11] = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1};
    double timestats[2] = {-1, -1};
    double repsilon = -1;
    double expepsilon = -1; //expected approx ratio
    uint64_t lbseti_expeps = 0; //those with lbset[][d] == lbseti_expeps are the ones guaranteeing expepsilon
  
    virtual void updateStats();

  };
} // namespace openwbo

#endif
