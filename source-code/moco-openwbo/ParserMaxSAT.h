/*!
 * \author Ruben Martins - ruben@sat.inesc-id.pt
 *
 * @section LICENSE
 *
 * MiniSat,  Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson
 *           Copyright (c) 2007-2010, Niklas Sorensson
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

#ifndef ParserMaxSAT_h
#define ParserMaxSAT_h

#include <stdio.h>

#include "MaxSATFormula.h"
#include "core/SolverTypes.h"
#include "utils/ParseUtils.h"

#ifdef HAS_EXTRA_STREAMBUFFER
#include "utils/StreamBuffer.h"
#endif

#include <vector>

using NSPACE::mkLit;
using NSPACE::StreamBuffer;

namespace openwbo {
  using Glucose::Solver;
//=================================================================================================
// DIMACS Parser:

template <class B> static uint64_t parseWeight(B &in) {
  uint64_t val = 0;
  while ((*in >= 9 && *in <= 13) || *in == 32)
    ++in;
  if (*in < '0' || *in > '9')
    fprintf(stderr, "PARSE ERROR! Unexpected char: %c\n", *in), exit(3);
  while (*in >= '0' && *in <= '9')
    val = val * 10 + (*in - '0'), ++in;
  return val;
}

template <class B, class MaxSATFormula>
static uint64_t readClause(B &in, MaxSATFormula *maxsat_formula,
                           vec<Lit> &lits) {
  int parsed_lit, var;
  int64_t weight = 1;
  lits.clear();
  if (maxsat_formula->getProblemType() == _WEIGHTED_)
    weight = parseWeight(in);
  assert(weight > 0);

  for (;;) {
    parsed_lit = parseInt(in);
    if (parsed_lit == 0)
      break;
    var = abs(parsed_lit) - 1;
    while (var >= maxsat_formula->nVars())
      maxsat_formula->newVar();
    lits.push((parsed_lit > 0) ? mkLit(var) : ~mkLit(var));
  }
  return weight;
}
template <class B, class MaxSATFormula>
static void readClauseH(B &in, MaxSATFormula *maxsat_formula,
                           vec<Lit> &lits) {
  int parsed_lit, var;
  lits.clear();

  for (;;) {
    parsed_lit = parseInt(in);
    if (parsed_lit == 0)
      break;
    var = abs(parsed_lit) - 1;
    while (var >= maxsat_formula->nVars())
      maxsat_formula->newVar();
    lits.push((parsed_lit > 0) ? mkLit(var) : ~mkLit(var));
  }
  return;
}

template <class B, class MaxSATFormula>
static void parseMaxSAT(B &in, MaxSATFormula *maxsat_formula) {
  vec<Lit> lits;
  uint64_t hard_weight = UINT64_MAX;
  for (;;) {
    skipWhitespace(in);
    if (*in == EOF)
      break;
    else if (*in == 'p') {
      if (eagerMatch(in, "p cnf")) {
        parseInt(in); // Variables
        parseInt(in); // Clauses
      } else if (eagerMatch(in, "wcnf")) {
        maxsat_formula->setProblemType(_WEIGHTED_);
        parseInt(in); // Variables
        parseInt(in); // Clauses
        if (*in != '\r' && *in != '\n') {
          hard_weight = parseWeight(in);
          maxsat_formula->setHardWeight(hard_weight);
        }
      } else
        printf("c PARSE ERROR! Unexpected char: %c\n", *in),
            printf("s UNKNOWN\n"), exit(_ERROR_);
    } else if (*in == 'c' || *in == 'p')
      skipLine(in);
    else if (*in == 'h'){
      ++in;
      readClauseH(in, maxsat_formula, lits);
    }    else {
      uint64_t weight = readClause(in, maxsat_formula, lits);
      if (weight < hard_weight ||
          maxsat_formula->getProblemType() == _UNWEIGHTED_) {
        assert(weight > 0);
        // Updates the maximum weight of soft clauses.
        maxsat_formula->setMaximumWeight(weight);
        // Updates the sum of the weights of soft clauses.
        maxsat_formula->updateSumWeights(weight);
        maxsat_formula->addSoftClause(weight, lits);
      } else
        maxsat_formula->addHardClause(lits);
    }
  }
}

// Inserts problem into solver.
//
template <class MaxSATFormula>
static void parseMaxSATFormula(gzFile input_stream,
                               MaxSATFormula *maxsat_formula) {
  StreamBuffer in(input_stream);
  parseMaxSAT(in, maxsat_formula);
  if (maxsat_formula->getMaximumWeight() == 1)
    maxsat_formula->setProblemType(_UNWEIGHTED_);
  else
    maxsat_formula->setProblemType(_WEIGHTED_);

  // maxsat_formula->setInitialVars(maxsat_formula->nVars());
}

//=================================================================================================
// Adaptation for mocnf

typedef struct _tmp_kp_t{
    int nbase;
    int ncoeffs;
    uint64_t * base;                
    uint64_t * mxb_coeffs;         
    int ** ls;
    uint64_t * ls_szs;
    int ** zs;
    uint64_t * zs_szs;
    bool haszs;
    
    uint64_t nvars;
    uint64_t nclauses;
} tmp_kp_t;

template <class B>
static void readLineVals(B &in, std::vector<int> &tmp){
    tmp.clear();
    while(*in == 32) ++in;
    while(*in != '\r' && *in != '\n' && *in != '0'){
        tmp.push_back(parseInt(in));
//         printf("->%d\n", tmp.back());
        while(*in == 32) ++in;
    }
    if(*in == '0') ++in;
    
    while ((*in >= 9 && *in <= 13) || *in == 32)
        ++in;
    
//     printf("tmp size: %d\n", tmp.size());
}

template <class B>
static void parseOriginalVar(B &in, char *varName, int *varNameSize) {
    int i = 0;
    while (isgraph(varName[i] = *in)) {
      i++;
      ++in;
    }
//     --in;
    varName[i] = '\0';
    *varNameSize = i;
  }


  
  
  // Creates a new variable in the SAT solver
void newSATVariable(Solver *S);
  
int get_new_Lit_id(Solver * S);

  
//Nota: Assume que so tem uma funcao objectivo
void parseMOCNF(char * fname, Solver *S, tmp_kp_t * tmp, std::map<std::string, int> name2id, std::map<int, bool> id2sign, Lit relax_var);
/*
void free_lsids(int ** ls_ids, int nls){
    for(int i = 0; i < nls; i++) free(ls_ids[i]);
    free(ls_ids);
}*/



//=================================================================================================
} // namespace openwbo

#endif
