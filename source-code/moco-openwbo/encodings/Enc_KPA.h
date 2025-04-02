/*!
 * \author Andreia P. Guerreiro
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

#ifndef Enc_KP_h
#define Enc_KP_h

#include <memory>
#ifdef SIMP
#include "simp/SimpSolver.h"
#else
#include "core/Solver.h"
#endif


#include "../FormulaPB.h"

#include "Encodings.h"
#include "core/SolverTypes.h"
#include <map>
#include <utility>
#include <vector>
#include <set>
#include "../maxConsts.h"
using NSPACE::lbool;

namespace openwbo {
    
  typedef std::map<uint64_t, Lit> wlit_mapt;
  enum kp_enc_mode {_kp_sat4j_, _kp_minisatp_};
  enum kp_enc_version {_full_, _inc_};
  enum kp_approx_mode {_ap_outvars_, _ap_coeffs_};
    
  class KPA : public Encodings {

  public:
    KPA() {
    
      nb_clauses = 0;
      nb_variables = 0;
    
      current_pb_rhs = -1;
    
      enc_version = _full_;
      enc_mode = _kp_sat4j_;
    
      nbase = 0;
      ncoeffs = 0;
    
      f_offset = 0;
    
      hasEncoding = false;
    
      pb = NULL;
      pb_approx_coeffs = NULL;
    
      epsilon = 1;
      approxMode = _ap_outvars_;
    
      ntmp_vars = 0;
    
      hasRelaxVar = false;
      relax_var = lit_Undef;
    
      base = NULL;
      mxb_coeffs = NULL;
      ls = NULL;
      ls_szs = NULL;
    
    }
    ~KPA() {
      if(hasEncoding){
        //TODO: free stuff
        delete(pb);
        free(base);
        free(mxb_coeffs);
            
        //ls
        for(int i = 0; i < ncoeffs; i++) free(ls[i]);
        free(ls);
        free(ls_szs);
        //zs
        
        if(epsilon > 1 && approxMode == _ap_coeffs_)
	  delete(pb_approx_coeffs);
      }
    }

    // Encode constraint.
    void encode(Solver *S, vec<Lit> &lits, vec<uint64_t> &coeffs, uint64_t rhs, int mode=_kp_sat4j_, int version=_full_);
    
    // Update constraint.
    void update(Solver *S, uint64_t rhs);
    
    wlit_mapt getRootLits() { return oliterals; } 
    
    Lit getValRootLit(Solver *S, uint64_t val);
    
    
    std::vector<Lit>& getAssumptions() { return assumptions; }
    std::vector<Lit> assumptions; //o encode() actualiza as assumptions
    
    void getEncodeSizes(int *nvar, int *nclauses, int *nrootvars){
      *nvar = nb_variables;
      *nclauses = nb_clauses;
      *nrootvars = oliterals.size();
    }
    
    bool hasCreatedEncoding() { return hasEncoding; }
    
    
    void evalModel(vec<lbool> &currentModel, uint64_t * objv, uint64_t *ap_objv, int di);

    void setApproxRatio(Solver *S, float eps, int mode=_ap_outvars_);
    void forgetEncoding(Solver *S);
    
    uint64_t getMinimum(){return f_offset;} //o menor valor possivel da funcao
    
    bool apFunctionIsExact();
  void fixed_vars(const std::set<Lit>& fvo){fv = std::make_shared<const std::set<Lit>>(fvo);};
protected:

    // encoding options
    int enc_version = _full_;       // (full, incremental) 
    int enc_mode;                   // method used to do the actual encoding
                                    //(based on sat4j, original code, direct usage of original code)
    // encoding: additional info
    PB * pb;                          //pseudo boolean function
    
    PB * pb_approx_coeffs;            //pseudo boolean function - with coeffs rounded down (_ap_coeffs_ approx mode)
    
    // encoding: core data         
    int nbase;                      
    int ncoeffs;
    uint64_t * base;                //ex. base = [2,3,3]
    uint64_t * mxb_coeffs;          // ex. mxb_coeffs = [1,2,6,18]  - mixed base coeffs
    int ** ls;                      // ex. ls = [[l1_1],[l2_1], [l6_1,l6_2],[l18_1]]
    uint64_t * ls_szs;                   // number of ls per coeff, ex.: [1,1,2,1]
    int ** zs;                      // ex. ls = [[z1_1],[z2_1], [z6_1,z6_2],[z18_1]]

    
    // encoded objective values
    wlit_mapt oliterals;            // objvalues and corresponding literals
    
    std::map<std::pair<uint64_t,int>,Lit> olits_aux; //lits that encodes 'leq' than pair[0], up to bit pair[1]
    
    // Number of variables and clauses for statistics.
    uint64_t nb_variables;
    uint64_t nb_clauses;
	
    bool hasEncoding;
    uint64_t current_pb_rhs;
    std::map<std::string, int> name2id;
    std::map<int, bool> id2sign;
    
    //file names (sat4j, ...)
    char ifname[max::name];
    char ofname[max::name];
    
    uint64_t f_offset; // valor que deve ser somado a funcao codificada (que corresponde a parte da
                       // funcao original), que corresponde a soma dos valores das vars com valor fixo
    
    // -------------------------------- methods ------------------------------- //
    // ------------------------------------------------------------------------ //
    
    // converters
    uint64_t convert_kp_sat4j(Solver *S, uint64_t ub);
    uint64_t convert_kp_minisatp(Solver * S, uint64_t ub);
    //------------
    void force_sort_within_ls(Solver *S);
    //------------
    Lit less_or_equal_than(Solver * S, uint64_t rhs);
    
    void encode_one_var_func(Solver *S);
    
    
    //IO auxiliar functions of convert_kp_sat4j
    uint64_t printPBfunction(Solver *S);
    void run_sat4j(uint64_t ub);
    void run_kp_minisatp(uint64_t ub);
    void read_mocnf(Solver *S);
    
    Lit getNewLit(Solver *S);
    
    void feasibleObjectiveValuesDP(Solver * S, std::set<uint64_t> &vals, uint64_t sz);
    void order_enc_rootvars(Solver *S);
    
    // info functions
    void printLit(Lit l) { printf("%s%d\n", sign(l) ? "-" : "", var(l) + 1); }
    
    uint64_t ntmp_vars; //number of vars in PB that are not fixed by the solver
    
    // approximation part
    float epsilon; //approximation ratio of (1+epsilon)
    int approxMode; // outvars represent a range of values (_ap_outvars_), round coefficients (_ap_coeffs_)
    
    //para o reencoding do _ap_coeffs_
    Lit relax_var;
    bool hasRelaxVar;
    
    std::set<uint64_t> vals;
    std::shared_ptr<const std::set<Lit>> fv{};

  };

} // namespace openwbo


#endif
