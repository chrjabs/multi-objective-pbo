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


#include "utils/Options.h"
#include "utils/ParseUtils.h"
#include "utils/System.h"
#include <cstddef>
#include <errno.h>
#include <signal.h>
#include <zlib.h>

#include <fstream>
#include <iostream>
#include <map>
#include <stdlib.h>
#include <libgen.h>
#include <string>
#include <vector>

#ifdef SIMP
#include "simp/SimpSolver.h"
#else
#include "core/Solver.h"
#endif

#include "MaxSAT.h"
#include "MaxTypes.h"
#include "ParserMaxSAT.h"
#include "ParserPB.h"
#include "Mcnf.h"

// Algorithms
#include "algorithms/Alg_LinearSU.h"
#include "algorithms/Alg_MSU3.h"
#include "algorithms/Alg_OLL.h"
#include "algorithms/Alg_PartMSU3.h"
#include "algorithms/Alg_WBO.h"
#include "algorithms/Alg_BLS.h"
#include "algorithms/Alg_Naive.h"
#include "algorithms/Alg_UnsatSatMO.h"
#include "algorithms/Alg_UnsatSatStratMO.h"
#include "algorithms/Alg_UnsatSatStratMSU3MO.h"
#include "algorithms/Alg_UnsatSatMSU3MO.h"
#include "algorithms/Alg_PMinimalMO.h"
#include "algorithms/Alg_HittingSetsMO.h"
#include "algorithms/Alg_HittingSetsStratMO.h"
#include "bounds/BoundsCalculatorMaker.h"
#include "context.h"
#include "context.h"
#include "maxConsts.h"
// path used when calling the articact

#define VER1_(x) #x
#define VER_(x) VER1_(x)
#define SATVER VER_(SOLVERNAME)
#define VER VER_(VERSION)

using NSPACE::cpuTime;
using NSPACE::OutOfMemoryException;
using NSPACE::IntOption;
using NSPACE::DoubleOption;
using NSPACE::BoolOption;
using NSPACE::StringOption;
using NSPACE::IntRange;
using NSPACE::DoubleRange;
using NSPACE::parseOptions;
using namespace openwbo;

#define INSTANCE argv[1]


//=================================================================================================

static MaxSAT *mxsolver;

static void SIGINT_exit(int signum) {
  cout<<"c current results\n";
  if(mxsolver != NULL)
    mxsolver->printAnswer(_UNKNOWN_);
  else
    cout<<"c solver has nothing to report;\n";
  exit(_UNKNOWN_);
}

//=================================================================================================
#if !defined(_MSC_VER) && !defined(__MINGW32__)
void limitMemory(uint64_t max_mem_mb)
{
// FIXME: OpenBSD does not support RLIMIT_AS. Not sure how well RLIMIT_DATA works instead.
#if defined(__OpenBSD__)
#define RLIMIT_AS RLIMIT_DATA
#endif

    // Set limit on virtual memory:
    if (max_mem_mb != 0){
        rlim_t new_mem_lim = (rlim_t)max_mem_mb * 1024*1024;
        rlimit rl;
        getrlimit(RLIMIT_AS, &rl);
        if (rl.rlim_max == RLIM_INFINITY || new_mem_lim < rl.rlim_max){
            rl.rlim_cur = new_mem_lim;
            if (setrlimit(RLIMIT_AS, &rl) == -1)
                printf("c WARNING! Could not set resource limit: Virtual memory.\n");
        }
    }

#if defined(__OpenBSD__)
#undef RLIMIT_AS
#endif
}
#else
void limitMemory(uint64_t /*max_mem_mb*/)
{
    printf("c WARNING! Memory limit not supported on this architecture.\n");
}
#endif


#if !defined(_MSC_VER) && !defined(__MINGW32__)
void limitTime(uint32_t max_cpu_time)
{
    if (max_cpu_time != 0){
        rlimit rl;
        getrlimit(RLIMIT_CPU, &rl);
        if (rl.rlim_max == RLIM_INFINITY || (rlim_t)max_cpu_time < rl.rlim_max){
            rl.rlim_cur = max_cpu_time;
            if (setrlimit(RLIMIT_CPU, &rl) == -1)
                printf("c WARNING! Could not set resource limit: CPU-time.\n");
        }
    }
}
#else
void limitTime(uint32_t /*max_cpu_time*/)
{
    printf("c WARNING! CPU-time limit not supported on this architecture.\n");
}
#endif

//=================================================================================================
// Main:

namespace context
{
  
  char * cwd;

  static void readLinkRecurse(){
    int result = 1;
    char expanded[max::name + 1];
    result = readlink(cwd, expanded, max::name);
    while(result > 0 ){
      result = readlink(cwd, expanded, max::name);
      snprintf(cwd, 2*max::name, "%s", expanded);
      cwd[result]= '\0';
    }
  
  }
  void setCwd(char **argv)
  {
    // setting the name of cwd.
    cwd = new char[max::name];
    snprintf(cwd, max::name, "%s", argv[0]);      
    readLinkRecurse();
    cwd = dirname(cwd);
    printf("c cwd was expanded to %s\n", cwd);


  }
  
}

namespace options
{
    BoolOption printmodel("Open-WBO", "print-model", "Print model.\n", true);

    StringOption printsoft("Open-WBO", "print-unsat-soft", "Print unsatisfied soft claues in the optimal assignment.\n", NULL);
    
    StringOption saveoutput("Open-WBO", "save-my-output", "Print efficient solutions, nondominated points, the lower bound set, and other info to separate files.\n", NULL);

    IntOption verbosity("Open-WBO", "verbosity",
                        "Verbosity level (0=minimal, 1=more).\n", 0,
                        IntRange(0, 1));

    IntOption cpu_lim("Open-WBO", "cpu-lim",
                      "Limit on CPU time allowed in seconds.\n", 0,
                      IntRange(0, INT32_MAX));

    IntOption mem_lim("Open-WBO", "mem-lim",
                      "Limit on memory usage in megabytes.\n", 0,
                      IntRange(0, INT32_MAX));

    IntOption algorithm("Open-WBO", "algorithm",
                        "Search algorithm "
                        "(0=wbo,1=linear-su,2=msu3,3=part-msu3,4=oll,5=best,6=mcs-mo,7=mo-us, 8=mo-pm, 9=mo-hs, 10=mo-usMSU3,11=mo-naive, 12=mo-us_strat, 13=mo-us_strat_msu3)."
                        "\n",
                        5, IntRange(0, 13));

    IntOption partition_strategy("PartMSU3", "partition-strategy",
                                 "Partition strategy (0=sequential, "
                                 "1=sequential-sorted, 2=binary)"
                                 "(only for unsat-based partition algorithms).",
                                 2, IntRange(0, 2));

    IntOption graph_type("PartMSU3", "graph-type",
                         "Graph type (0=vig, 1=cvig, 2=res) (only for unsat-"
                         "based partition algorithms).",
                         2, IntRange(0, 2));

    BoolOption bmo("Open-WBO", "bmo", "BMO search.\n", false);

    IntOption cardinality("Encodings", "cardinality",
                          "Cardinality encoding (0=cardinality networks, "
                          "1=totalizer, 2=modulo totalizer).\n",
                          1, IntRange(0, 2));

    IntOption amo("Encodings", "amo", "AMO encoding (0=Ladder).\n", 0,
                  IntRange(0, 0));

    IntOption pb("Encodings", "pb", "PB encoding (0=SWC,1=GTE,2=Adder,3=IGTE, 4=KP).\n", 2,
                 IntRange(0, 4));
    
    IntOption pbobjf("Encodings", "pbobjf", "PB encoding (0=SWC,1=GTE,2=Adder,3=IGTE, 4=KP, 5=KP_minisatp).\n", 4,
                 IntRange(0, 5));
    
    IntOption blockingMode("blockingMode", "blockingMode", "IGTE blocking (0=one by one mode,1=aggressive mode).\n", 0,
                 IntRange(0, 1));
    
    IntOption sstrategy("Alg_BLS", "sstrategy", "Search strategy of Alg_BLS.\n", 20,
                 IntRange(0, 21));
    
    BoolOption clbounds("tight bounds computation", "clbounds", "Compute tight lower bounds?(requires CPLEX)\n", true);
    BoolOption cubounds("tight bounds computation", "cubounds", "Compute tight upper bounds?(requires CPLEX)\n", true);

    IntOption apmode("Encoding: KPA", "apmode", "approximation mode (0=outvars, 1=coeffs).\n", 0,
                 IntRange(0, 2));
    
    DoubleOption eps("Encoding: IGTE/KPA", "eps", "epsilon value (>= 1).\n", 1,
                 DoubleRange(1, true, INT32_MAX, false));
    
    DoubleOption redFactor("Encoding: IGTE", "redFactor", "epsilon's reduction factor\n", -1,
                 DoubleRange(2, true, 100, false));
    
    BoolOption fulligte("Encoding", "fulligte", "reduce approx error until eps=1 (run complete version).\n", true);

    IntOption formula("Open-WBO", "formula",
                      "Type of formula (0=WCNF, 1=OPB, 2=MCNF).\n", 1, IntRange(0, 2));

    IntOption weight(
        "WBO", "weight-strategy",
        "Weight strategy (0=none, 1=weight-based, 2=diversity-based).\n", 2,
        IntRange(0, 2));

    BoolOption symmetry("WBO", "symmetry", "Symmetry breaking.\n", true);

    IntOption symmetry_lim(
        "WBO", "symmetry-limit",
        "Limit on the number of symmetry breaking clauses.\n", 500000,
        IntRange(0, INT32_MAX));

  IntOption partition_parameter("partition parameter", "part_par",". Larger values result in larger parts\n", 15);

    BoolOption prepro("MO Preprocessing", "prepro", "MaxPre preprocessing (only with MCNF input).\n", false);

    StringOption prepro_tech("MO Preprocessing", "prepro-tech", "Techniques to use in MaxPre.\n", "[[uvsrgc]VRTG]");

    BoolOption objective_cores("IHS", "obj-cores", "Extract cores only over objective variables\n", false);
}

void setLimits()
{
  using namespace options;
  // Try to set resource limits:
  if (cpu_lim != 0) limitTime(cpu_lim);
  if (mem_lim != 0) limitMemory(mem_lim);
  
}

MaxSAT * buildSolver()
{
  using namespace options;
  MaxSAT *S = nullptr;
  
    switch ((int)algorithm) {
    case _ALGORITHM_WBO_:
      S = new WBO(verbosity, weight, symmetry, symmetry_lim);
      break;

    case _ALGORITHM_LINEAR_SU_:
      S = new LinearSU(verbosity, bmo, cardinality, pb);
      break;

    case _ALGORITHM_PART_MSU3_:
      S = new PartMSU3(verbosity, partition_strategy, graph_type, cardinality);
      break;

    case _ALGORITHM_MSU3_:
      S = new MSU3(verbosity);
      break;

    case _ALGORITHM_OLL_:
      S = new OLL(verbosity, cardinality);
      break;

    case _ALGORITHM_BEST_:
      break;
      
    case _ALGORITHM_MCSMO_:
      S = new BLS(verbosity, weight, partition_strategy, cardinality, pb, pbobjf, apmode, eps, sstrategy, redFactor);
      break;

    case _ALGORITHM_MOUNSATSAT_:
      S = new UnsatSatMO(verbosity, weight, partition_strategy, cardinality, pb, pbobjf, apmode, eps, sstrategy);
      break;

    case _ALGORITHM_MOPMINIMAL_:
      S = new PMinimalMO(verbosity, weight, partition_strategy, cardinality, pb, pbobjf, apmode, eps, sstrategy);
      break;

    case _ALGORITHM_MOHITTINGSETS_:
      S = new HittingSetsMO(verbosity, weight, partition_strategy, cardinality, pb, pbobjf, apmode, eps, sstrategy, -1, objective_cores);
      break;

    case _ALGORITHM_MOUNSATSATMSU3_:
      S = new UnsatSatMSU3MO(verbosity, weight, partition_strategy, cardinality, pb, pbobjf, apmode, eps, sstrategy);
      break;

    case _ALGORITHM_MONAIVE_:
      S = new Naive(verbosity, weight, partition_strategy, cardinality, pb, pbobjf, apmode, eps, sstrategy);
      break;
    case _ALGORITHM_MOUNSATSATSTRAT_:
      S = new UnsatSatStratMO(verbosity, weight, partition_strategy, cardinality, pb, pbobjf, apmode, eps, sstrategy, partition_parameter);
      break;
    case _ALGORITHM_MOUNSATSATSTRATMSU3_:
      S = new UnsatSatStratMSU3MO(verbosity, weight, partition_strategy, cardinality, pb, pbobjf, apmode, eps, sstrategy, partition_parameter);
      break;

    default:
      printf("c Error: Invalid MaxSAT algorithm.\n");
      printf("s UNKNOWN\n");
      exit(_ERROR_);
    }  
    printf("c solver typeid: %s\n", typeid(*S).name());
    return S;
}

void printHeader(){
  printf(
      "c\nc Open-WBO:\t a Modular MaxSAT Solver -- based on %s (%s version)\n",
      SATVER, VER);
  printf("c Version:\t September 2018 -- Release: 2.1\n");
  printf("c Authors:\t Ruben Martins, Vasco Manquinho, Ines Lynce\n");
  printf("c Contributors:\t Miguel Neves, Saurabh Joshi, Norbert Manthey, Mikolas Janota\n");
  printf("c Contact:\t open-wbo@sat.inesc-id.pt -- "
         "http://sat.inesc-id.pt/open-wbo/\nc\n");
}


MaxSATFormula* buildFormula(int argc, char **argv){
  using namespace options;

  MaxSATFormula *maxsat_formula = new MaxSATFormula();
  if ((int)formula == _FORMAT_MAXSAT_) {
    gzFile in = (argc == 1) ? gzdopen(0, "rb") : gzopen(INSTANCE, "rb");
    if (in == NULL)
      printf("c ERROR! Could not open file: %s\n",
             argc == 1 ? "<stdin>" : INSTANCE),
          printf("s UNKNOWN\n"), exit(_ERROR_);
    parseMaxSATFormula(in, maxsat_formula);
    maxsat_formula->setFormat(_FORMAT_MAXSAT_);
    gzclose(in);
  } else if ((int)formula == _FORMAT_PB_) {
    printf("c [Main] ParserPB()\n");
    ParserPB *parser_pb = new ParserPB();
    parser_pb->parsePBFormula(INSTANCE, maxsat_formula);

    auto tmp_calc = bounds::maker(maxsat_formula);
    bounds::BoundsCalculator_I &calc = *tmp_calc;
    for (int i = 0, n = maxsat_formula->nObjFunctions(); i < n; i++) {
      auto const &pb = maxsat_formula->getObjFunction(i);
      auto total = pb->total;
      auto _const = pb->_const;
      maxsat_formula->bounds[i][2] = total + _const;
      maxsat_formula->bounds[i][3] = total + _const;
      maxsat_formula->bounds[i][0] = 0 + _const;
      maxsat_formula->bounds[i][1] = 0 + _const;
      // Improve the bounds with cplex
      if (clbounds || cubounds) {
        calc.setObjective(pb);
        if (clbounds) {
          maxsat_formula->bounds[i][1] = calc.lowerBound();
        }
        if (cubounds) {
          maxsat_formula->bounds[i][2] = calc.upperBound();
        }
        printf("c [Main] bounds, obj %d: [%ld %ld %ld %ld]\n", i,
               maxsat_formula->bounds[i][0], maxsat_formula->bounds[i][1],
               maxsat_formula->bounds[i][2], maxsat_formula->bounds[i][3]);
      }
    };
    maxsat_formula->setFormat(_FORMAT_PB_);
  } else {
    parseMcnf(INSTANCE, maxsat_formula);

    if ((bool)prepro) {
      maxsat_formula->preprocess((const char *)prepro_tech);
    }

    for (int i = 0, n = maxsat_formula->nObjFunctions(); i < n; i++) {
      auto obj = maxsat_formula->getClausalObj(i);
      uint64_t sum = 0;
      for (int j = 0; j < obj->size(); j++)
        sum += (*obj)[j].weight;
      maxsat_formula->bounds[i][2] = sum + maxsat_formula->clausalOffset(i);
      maxsat_formula->bounds[i][3] = sum + maxsat_formula->clausalOffset(i);
      maxsat_formula->bounds[i][0] = maxsat_formula->clausalOffset(i);
      maxsat_formula->bounds[i][1] = maxsat_formula->clausalOffset(i);
    }

    maxsat_formula->setFormat(_FORMAT_MCNF_);
  }

  // DEBUGGING:
  // if (_FORMAT_MCNF_ == maxsat_formula->getFormat())
  //   maxsat_formula->clausalToPbObjs();
  // maxsat_formula->my_print();
  // exit(-1);

  return maxsat_formula;
}

void printConfig(MaxSATFormula* maxsat_formula, double initial_time){
  using namespace options;

    printf("c |                                                                "
           "                                       |\n");
    printf("c ========================================[ Problem Statistics "
           "]===========================================\n");
    printf("c |                                                                "
           "                                       |\n");

    if (maxsat_formula->getFormat() == _FORMAT_MAXSAT_)
      printf(
          "c |  Problem Format:  %17s                                         "
          "                          |\n",
          "MaxSAT");
    else
//       printf(
//           "c |  Problem Format:  %17s                                         "
//           "                          |\n",
//           "PB");
      printf(
          "c |  Problem Format:  %17s                                         "
          "                          |\n"
          "c |  Number of objectives:  %11d                                         "
          "                          |\n",
          "PB",
          maxsat_formula->nObjFunctions());

    if (maxsat_formula->getProblemType() == _UNWEIGHTED_)
      printf("c |  Problem Type:  %19s                                         "
             "                          |\n",
             "Unweighted");
    else
      printf("c |  Problem Type:  %19s                                         "
             "                          |\n",
             "Weighted");

    printf("c |  Number of variables:  %12d                                    "
           "                               |\n",
           maxsat_formula->nVars());
    printf("c |  Number of hard clauses:    %7d                                "
           "                                   |\n",
           maxsat_formula->nHard());
    printf("c |  Number of soft clauses:    %7d                                "
           "                                   |\n",
           maxsat_formula->nSoft());
    printf("c |  Number of cardinality:     %7d                                "
           "                                   |\n",
           maxsat_formula->nCard());
    printf("c |  Number of PB :             %7d                                "
           "                                   |\n",
           maxsat_formula->nPB());
    double parsed_time = cpuTime();

    printf("c |  Parse time:           %12.2f s                                "
           "                                 |\n",
           parsed_time - initial_time);
    printf("c |                                                                "
           "                                       |\n");

}

void loadSolver(MaxSATFormula* maxsat_formula, MaxSAT* S, double initial_time){
  using namespace options;
  if (algorithm == _ALGORITHM_BEST_) {
    assert(S == NULL);

    if (maxsat_formula->getProblemType() == _UNWEIGHTED_) {
      // Unweighted
      S = new PartMSU3(verbosity, _PART_BINARY_, RES_GRAPH,
		       cardinality);
      S->loadFormula(maxsat_formula);

      if (((PartMSU3 *)S)->chooseAlgorithm() == _ALGORITHM_MSU3_) {
	// FIXME: possible memory leak
	S = new MSU3(verbosity);
      }

    } else {
      // Weighted
      S = new OLL(verbosity, cardinality);
    }
  }

  if (S->getMaxSATFormula() == NULL){
    printf("c [Main] mxsolver->loadFormula()\n");
    S->loadFormula(maxsat_formula);
  }

  S->setPrintModel(printmodel);
  S->setPrintSoft((const char *)printsoft);
  S->setMyOutputFiles((const char *)saveoutput);
  S->setInitialTime(initial_time);
  S->setPrint(true);
}

int main(int argc, char **argv) {
    
  printHeader();
  context::setCwd(argv);

#if defined(__linux__)
  fpu_control_t oldcw, newcw;
  _FPU_GETCW(oldcw);
  newcw = (oldcw & ~_FPU_EXTENDED) | _FPU_DOUBLE;
  _FPU_SETCW(newcw);
  printf(
	 "c WARNING: for repeatability, setting FPU to use double precision\n");
#endif


  try {
    using namespace options;

    signal(SIGXCPU, SIGINT_exit);
    signal(SIGTERM, SIGINT_exit);
    signal(SIGINT, SIGINT_exit);

    NSPACE::setUsageHelp("c USAGE: %s [options] <input-file>\n\n");
    parseOptions(argc, argv, true);
    setLimits();
    MaxSAT *S = buildSolver();
    
    double initial_time = cpuTime();


    if (argc == 1) {
      printf("c Warning: no filename.\n");
      return 1;
    }
    
    MaxSATFormula *maxsat_formula = buildFormula(argc, argv);

    printConfig(maxsat_formula,  initial_time);

    loadSolver(maxsat_formula, S, initial_time);
    mxsolver = S;

    printf("c [Main] mxsolver->search()\n");
    int ret = (int)mxsolver->search();
    delete context::cwd;
    delete S;
    delete maxsat_formula;
   return ret;
  } catch (OutOfMemoryException &) {
    sleep(1);
    printf("c Error: Out of memory.\n");
    mxsolver->printAnswer(_MEMOUT_);
    exit(_ERROR_);
  } catch(MaxSATException &e) {
    sleep(1);
    printf("c Error: MaxSAT Exception: %s\n", e.getMsg());
    printf("s MAXSATERROR\n");
    exit(_ERROR_);
  }
}
