/*!
 * \author Saurabh Joshi - sbjoshi@iith.ac.in
 *
 * @section LICENSE
 *
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

#include "Enc_GTE.h"
#include <algorithm>
#include <numeric>

// #define __DEBUG__
#ifdef __DEBUG__
    #define    DPRINTF(f,...)    fprintf(stderr,f,##__VA_ARGS__)
#else
    #define    DPRINTF(f,...)
#endif




void my_printPB( vec<Lit> &lits, vec<uint64_t> &coeffs,
                 uint64_t rhs){
#ifdef __DEBUG__
    for (int i = 0; i < coeffs.size(); i++) {
        printf("%d ", (int)coeffs[i]);
        if (sign(lits[i]))
            printf("~");
        printf("y%d ", var(lits[i]) + 1);
    }
    printf("<= %lu\n", rhs);
#endif
}

void my_printLit(Lit l){
#ifdef __DEBUG__
    printf(" ");
    if (sign(l))
        printf("~");
    printf("y%d", var(l) + 1);
#endif
}


void my_printLitNL(Lit l){
#ifdef __DEBUG__
    my_printLit(l);
    printf("\n");
#endif
}

using namespace openwbo;




struct less_than_wlitt {
  inline bool operator()(const wlitt &wl1, const wlitt &wl2) {
    return (wl1.weight < wl2.weight);
  }
};

Lit GTE::getNewLit(Solver *S) {
  Lit p = mkLit(S->nVars(), false);
  newSATVariable(S);
  nb_variables++;
  return p;
}

Lit GTE::get_var(Solver *S, wlit_mapt &oliterals, uint64_t weight) {
  wlit_mapt::iterator it = oliterals.find(weight);
  if (it == oliterals.end()) {
    Lit v = getNewLit(S);
    oliterals[weight] = v;
  }
  return oliterals[weight];
}

bool GTE::encodeLeq(uint64_t k, Solver *S, const weightedlitst &iliterals,
                    wlit_mapt &oliterals) {

  if (iliterals.size() == 0 || k == 0){
    return false;
  }

  if (iliterals.size() == 1) {

    oliterals.insert(
        wlit_pairt(iliterals.front().weight, iliterals.front().lit));
    return true;
  }

  unsigned int size = iliterals.size();

  // formulat lformula,rformula;
  weightedlitst linputs, rinputs;
  wlit_mapt loutputs, routputs;

  unsigned int lsize = size >> 1;
  // unsigned int rsize=size-lsize;
  weightedlitst::const_iterator myit = iliterals.begin();
  weightedlitst::const_iterator myit1 = myit + lsize;
  weightedlitst::const_iterator myit2 = iliterals.end();

  linputs.insert(linputs.begin(), myit, myit1);
  rinputs.insert(rinputs.begin(), myit1, myit2);

  /*wlitt init_wlit;
  init_wlit.lit = lit_Undef;
  init_wlit.weight=0;*/
  wlit_sumt wlit_sum;
  uint64_t lk = std::accumulate(linputs.begin(), linputs.end(), uint64_t(0), wlit_sum);
  uint64_t rk = std::accumulate(rinputs.begin(), rinputs.end(), uint64_t(0), wlit_sum);

  lk = k >= lk ? lk : k;
  rk = k >= rk ? rk : k;

  bool result = encodeLeq(lk, S, linputs, loutputs);
  if (!result)
    return result;
  result = result && encodeLeq(rk, S, rinputs, routputs);
  if (!result)
    return result;

  {
    assert(!loutputs.empty());

    for (wlit_mapt::iterator mit = loutputs.begin(); mit != loutputs.end();
         mit++) {

      if (mit->first > k) {
          
        DPRINTF("\tc b1:\t"); my_printLit(~mit->second);
        my_printLitNL(get_var(S, oliterals, k));
        addBinaryClause(S, ~mit->second, get_var(S, oliterals, k));
        nb_clauses++;
      } else {
        DPRINTF("\tc b2:\t"); my_printLit(~mit->second);
        my_printLitNL(get_var(S, oliterals, mit->first));
        addBinaryClause(S, ~mit->second, get_var(S, oliterals, mit->first));
        nb_clauses++;
        // clause.push_back(get_var(auxvars,oliterals,l.first));
      }

      // formula.push_back(std::move(clause));
    }
  }

  {
    assert(!routputs.empty());
    for (wlit_mapt::iterator mit = routputs.begin(); mit != routputs.end();
         mit++) {

      if (mit->first > k) {
        DPRINTF("\tc b3:\t"); my_printLit(~mit->second);
        my_printLitNL(get_var(S, oliterals, k));
        addBinaryClause(S, ~mit->second, get_var(S, oliterals, k));
        nb_clauses++;
        // clause.push_back(get_var(auxvars,oliterals,k));
      } else {
        DPRINTF("\tc b4:\t"); my_printLit(~mit->second);
        my_printLitNL(get_var(S, oliterals, mit->first));
        addBinaryClause(S, ~mit->second, get_var(S, oliterals, mit->first));
        nb_clauses++;
        // clause.push_back(get_var(auxvars,oliterals,r.first));
      }

      // formula.push_back(std::move(clause));
    }
  }

  // if(!lformula.empty() && !rformula.empty())
  {
    for (wlit_mapt::iterator lit = loutputs.begin(); lit != loutputs.end();
         lit++) {
      for (wlit_mapt::iterator rit = routputs.begin(); rit != routputs.end();
           rit++) {
        /*clauset clause;
        clause.push_back(-l.second);
        clause.push_back(-r.second);*/
        uint64_t tw = lit->first + rit->first;
        if (tw > k) {
            
          DPRINTF("\tc t1:\t"); my_printLit(~lit->second);
          my_printLit(~rit->second);
          my_printLitNL(get_var(S, oliterals, k));
      
          addTernaryClause(S, ~lit->second, ~rit->second,
                           get_var(S, oliterals, k));
          nb_clauses++;
          // clause.push_back(get_var(auxvars,oliterals,k));
        } else {
          DPRINTF("\tc t2:\t"); my_printLit(~lit->second);
          my_printLit(~rit->second);
          my_printLitNL(get_var(S, oliterals, tw));
          
          addTernaryClause(S, ~lit->second, ~rit->second,
                           get_var(S, oliterals, tw));
          nb_clauses++;
          // clause.push_back(get_var(auxvars,oliterals,tw));
        }

        // formula.push_back(std::move(clause));
      }
    }
  }

  return true;
}



void GTE::encode(Solver *S, vec<Lit> &lits, vec<uint64_t> &coeffs,
                 uint64_t rhs) {
  // FIXME: do not change coeffs in this method. Make coeffs const.

  // If the rhs is larger than INT32_MAX is not feasible to encode this
  // pseudo-Boolean constraint to CNF.
  if (rhs >= UINT64_MAX) {
    DPRINTF("c Overflow in the Encoding\n");
    DPRINTF("s UNKNOWN\n");
    exit(_ERROR_);
  }
  
  DPRINTF("c [Enc_GTE: encode()]\n");
//   DPRINTF("c pb:");
//   my_printPB(lits, coeffs, rhs);

  hasEncoding = false;
  nb_variables = 0;
  nb_clauses = 0;

  vec<Lit> simp_lits;
  vec<uint64_t> simp_coeffs;
  lits.copyTo(simp_lits);
  coeffs.copyTo(simp_coeffs);

  lits.clear();
  coeffs.clear();

  // Fix literals that have a coeff larger than rhs.
  for (int i = 0; i < simp_lits.size(); i++) {
    if (simp_coeffs[i] == 0)
      continue;

    if (simp_coeffs[i] >= UINT64_MAX) {
      DPRINTF("c Overflow in the Encoding\n");
      DPRINTF("s UNKNOWN\n");
      exit(_ERROR_);
    }

    if (simp_coeffs[i] <= (unsigned)rhs) {
      lits.push(simp_lits[i]);
      coeffs.push(simp_coeffs[i]);
    } else{
        DPRINTF("\tc u1:\t"); my_printLitNL(~simp_lits[i]);
      addUnitClause(S, ~simp_lits[i]);
    }
  }

  if (lits.size() == 1) {
        DPRINTF("\tc e1: size 1\n");
    // addUnitClause(S, ~lits[0]);
    return;
  }

  if (lits.size() == 0){
    DPRINTF("\tc e0: size 0");
    return;
  }

  weightedlitst iliterals;
  for (int i = 0; i < lits.size(); i++) {
    wlitt wl;
    wl.lit = lits[i];
    wl.weight = coeffs[i];
    iliterals.push_back(wl);
  }
  less_than_wlitt lt_wlit;
  std::sort(iliterals.begin(), iliterals.end(), lt_wlit);
  encodeLeq(rhs + 1, S, iliterals, pb_oliterals);

  for (wlit_mapt::reverse_iterator rit = pb_oliterals.rbegin();
       rit != pb_oliterals.rend(); rit++) {
    if (rit->first > rhs) {
      DPRINTF("\tc u2:\t"); my_printLitNL(~rit->second);
      addUnitClause(S, ~rit->second);
    } else {
      break;
    }
  }
  // addUnitClause(S,~pb_oliterals.rbegin()->second);
  /*
  if (pb_oliterals.rbegin()->first != rhs+1){
        DPRINTF("%d - %d\n",pb_oliterals.rbegin()->first,rhs);
        for(wlit_mapt::reverse_iterator
  rit=pb_oliterals.rbegin();rit!=pb_oliterals.rend();rit++)
  {
        DPRINTF("rit->first %d\n",rit->first);
  }
  }
  */
  // assert (pb_oliterals.rbegin()->first == rhs+1);
  // printLit(~pb_oliterals.rbegin()->second);
  /* ... PUT CODE HERE FOR CREATING THE ENCODING ... */
  /* ... do not forget to sort the coefficients so that GTE is more efficient
   * ... */

  current_pb_rhs = rhs;
  hasEncoding = true;
  printf("c lgte nvars nclauses nvars_root\n"); 
  printf("c gte\t%d\t%d\t%lu\n", nb_variables, nb_clauses, pb_oliterals.size());

  
  DPRINTF("\tc ncs, nvs:\t%d\t%d\n", nb_clauses, nb_variables);
  fflush(stdout);
}

void GTE::update(Solver *S, uint64_t rhs) {

  assert(hasEncoding);
  for (wlit_mapt::reverse_iterator rit = pb_oliterals.rbegin();
       rit != pb_oliterals.rend(); rit++) {
    if (rit->first > current_pb_rhs)
      continue;
    if (rit->first > rhs) {
        DPRINTF("\tc turn off1 "); my_printLitNL(~rit->second); DPRINTF("[%d]\n", rit->first);
      addUnitClause(S, ~rit->second);
    } else {
      break;
    }
  }
  /* ... PUT CODE HERE TO UPDATE THE RHS OF AN ALREADY EXISTING ENCODING ... */

  current_pb_rhs = rhs;
}

// TODO: refactor the code to reduce duplication for the predict methods

// predict number of variables and clauses that this encode will generate
Lit GTE::get_var_predict(Solver *S, wlit_mapt &oliterals, uint64_t weight) {
  wlit_mapt::iterator it = oliterals.find(weight);
  if (it == oliterals.end()) {
    Lit v = mkLit(nb_current_variables, false);
    nb_current_variables++;
    oliterals[weight] = v;
  }
  return oliterals[weight];
}

bool GTE::predictEncodeLeq(uint64_t k, Solver *S, const weightedlitst &iliterals,
                    wlit_mapt &oliterals) {


  if (nb_clauses_expected >= _MAX_CLAUSES_)
    return false;

  if (iliterals.size() == 0 || k == 0)
    return false;

  if (iliterals.size() == 1) {
    
    oliterals.insert(
        wlit_pairt(iliterals.front().weight, iliterals.front().lit));
    return true;
  }

  unsigned int size = iliterals.size();

  weightedlitst linputs, rinputs;
  wlit_mapt loutputs, routputs;

  unsigned int lsize = size >> 1;
  weightedlitst::const_iterator myit = iliterals.begin();
  weightedlitst::const_iterator myit1 = myit + lsize;
  weightedlitst::const_iterator myit2 = iliterals.end();

  linputs.insert(linputs.begin(), myit, myit1);
  rinputs.insert(rinputs.begin(), myit1, myit2);

  wlit_sumt wlit_sum;
  uint64_t lk = std::accumulate(linputs.begin(), linputs.end(), uint64_t(0), wlit_sum);
  uint64_t rk = std::accumulate(rinputs.begin(), rinputs.end(), uint64_t(0), wlit_sum);

  lk = k >= lk ? lk : k;
  rk = k >= rk ? rk : k;

  bool result = predictEncodeLeq(lk, S, linputs, loutputs);
  if (!result)
    return result;
  result = result && predictEncodeLeq(rk, S, rinputs, routputs);
  if (!result)
    return result;

  {
    assert(!loutputs.empty());
    for (wlit_mapt::iterator mit = loutputs.begin(); mit != loutputs.end();
         mit++) {
      if (mit->first > k) {
        get_var_predict(S, oliterals, k);
        nb_clauses_expected++;
        if (nb_clauses_expected >= _MAX_CLAUSES_)
          return false;
      } else {
        get_var_predict(S, oliterals, mit->first);
        nb_clauses_expected++;
        if (nb_clauses_expected >= _MAX_CLAUSES_)
          return false;
      }

    }
  }

  {
    assert(!routputs.empty());
    for (wlit_mapt::iterator mit = routputs.begin(); mit != routputs.end();
         mit++) {
      if (mit->first > k) {
        get_var_predict(S, oliterals, k);
        nb_clauses_expected++;
        if (nb_clauses_expected >= _MAX_CLAUSES_)
          return false;
      } else {
        get_var_predict(S, oliterals, mit->first);
        nb_clauses_expected++;
        if (nb_clauses_expected >= _MAX_CLAUSES_)
          return false;
      }

    }
  }

  {
    for (wlit_mapt::iterator lit = loutputs.begin(); lit != loutputs.end();
         lit++) {
      for (wlit_mapt::iterator rit = routputs.begin(); rit != routputs.end();
           rit++) {
        uint64_t tw = lit->first + rit->first;
        if (tw > k) {
          get_var_predict(S, oliterals, k);
          nb_clauses_expected++;
          if (nb_clauses_expected >= _MAX_CLAUSES_)
            return false;
        } else {
          get_var_predict(S, oliterals, tw);
          nb_clauses_expected++;
          if (nb_clauses_expected >= _MAX_CLAUSES_)
            return false;
        }

      }
    }
  }

  return true;
}

int GTE::predict(Solver *S, vec<Lit> &lits, vec<uint64_t> &coeffs,
                 uint64_t rhs) {

  vec<Lit> simp_lits;
  vec<uint64_t> simp_coeffs;
  lits.copyTo(simp_lits);
  coeffs.copyTo(simp_coeffs);
  lits.clear();
  coeffs.clear();

  nb_current_variables = S->nVars();

  for (int i = 0; i < simp_lits.size(); i++) {
    if (simp_coeffs[i] == 0)
      continue;

    if (simp_coeffs[i] >= INT64_MAX) {
      return _MAX_CLAUSES_;
    }

    if (simp_coeffs[i] <= (unsigned)rhs) {
      lits.push(simp_lits[i]);
      coeffs.push(simp_coeffs[i]);
    } else
      nb_clauses_expected++;
  }

  if (lits.size() <= 1)
    return nb_clauses_expected;
  
  weightedlitst iliterals;
  for (int i = 0; i < lits.size(); i++) {
    wlitt wl;
    wl.lit = lits[i];
    wl.weight = coeffs[i];
    iliterals.push_back(wl);
  }
  less_than_wlitt lt_wlit;
  std::sort(iliterals.begin(), iliterals.end(), lt_wlit);
  wlit_mapt oliterals;
  predictEncodeLeq(rhs + 1, S, iliterals, oliterals);
  return nb_clauses_expected;

}
