/*!
 * \author Andreia P. Guerreiro - andreia.guerreiro@tecnico.ulisboa.pt
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

#include "../maxConsts.h"
#include "Enc_KPA.h"
#include <algorithm>
#include <numeric>

#include <stdio.h>
#include <stdlib.h>
#include <set>

#include "../context.h"
#include "../depsLib/sat4jEncoder.h"
#include "../ParserMaxSAT.h"

// #define __DEBUG__
#ifdef __DEBUG__
    #define    DPRINTF(f,...)    fprintf(stderr,f,##__VA_ARGS__)
#else
    #define    DPRINTF(f,...)
#endif



using namespace openwbo;


Lit KPA::getNewLit(Solver *S) {
  Lit p = mkLit(S->nVars(), false);
  newSATVariable(S);
  nb_variables++;
  return p;
}


int64_t MAXDIG = 7;

void KPA::encode(Solver *S, vec<Lit> &lits, vec<uint64_t> &coeffs,
                 uint64_t rhs, int mode, int version) {
    
    printf("c [KP::encode] rhs: %lu, epsilon: %f\n", rhs, epsilon);
    // If the rhs is larger than INT32_MAX is not feasible to encode this
    // pseudo-Boolean constraint to CNF.
    if (rhs >= UINT64_MAX) {
        DPRINTF("c Overflow in the Encoding\n");
        DPRINTF("s UNKNOWN\n");
        exit(_ERROR_);
    }
  
    uint64_t nb_clauses_prev = nb_clauses;

//     hasEncoding = false;
    enc_mode = mode;
    enc_version = version;
    
    if(!hasEncoding){
        nb_variables = 0;
        nb_clauses = 0;
        pb = new PB(lits, coeffs, rhs);
    }
    
    
    if(!hasEncoding || (approxMode == _ap_coeffs_)){ //no caso do _ap_coeffs_ so vale a pena se pb!=pb_approx_coeffs
        switch(enc_mode){
            case _kp_sat4j_:
                ntmp_vars = convert_kp_sat4j(S, rhs);
                break;
                
            case _kp_minisatp_:
                ntmp_vars = convert_kp_minisatp(S, rhs);
                break;
            default:
                printf("ERROR: Wrong encoding type (Enc_KP)");
                exit(1);
        }
    }
    
  current_pb_rhs = rhs;
  hasEncoding = true;
  printf("c lkpa nvars nclauses nvars_root\n"); 
  printf("c kpa\t%lu\t%lu\t%lu\n", nb_variables, nb_clauses, oliterals.size());

  fflush(stdout);
  
  uint64_t lessclauses = nb_clauses;
  int nobjvals = 0;
  
  if(ntmp_vars > 0){
    
    if(enc_version == _full_){
        if(epsilon <= 1 || approxMode == _ap_coeffs_){
            
            if(epsilon > 1){
                rhs = f_offset;
                for(int i = 0; i < pb_approx_coeffs->_coeffs.size(); i++){
		  if(true){
                        rhs += pb_approx_coeffs->_coeffs[i];
                    }
                }
                if(current_pb_rhs < rhs)
                    rhs = current_pb_rhs;
                printf("c new rhs: %lu\n", rhs); 
            }
            

            vals.clear();
            
            feasibleObjectiveValuesDP(S, vals, rhs+1-f_offset);
            
            printf("c [KPA] encode full version -- all obj values (rhs: %lu)\n", rhs);

            printf("c vals size: %lu\n", vals.size());
            for(std::set<uint64_t>::iterator it = vals.begin(); it != vals.end(); it++){
#ifdef __DEBUG__
                      printf("c obj val: %lu\n", *it);
#endif
                less_or_equal_than(S, *it+f_offset);
                
            }
                        
        }else{
            printf("c range obj values (f_offset: %lu)\n", f_offset);
            
            if(f_offset == 0){
#ifdef __DEBUG__
                printf("c r_lb %d\n", 0);
#endif
                less_or_equal_than(S, 0);
            }
            
            int ndigit = 2*log10(ntmp_vars);
            printf("c ndigit: %d\n", ndigit);

            if(vals.size() == 0 && (epsilon <= 1 || ndigit <= MAXDIG)){
                //so e preciso calcular uma vez
                printf("\n\nc feasibleObjectiveValuesDP!\n");
                feasibleObjectiveValuesDP(S, vals, rhs+1-f_offset);
            }
            
            uint64_t lu = rhs;
            for(int i = 0; i < pb->_coeffs.size(); i++){
                lu = (pb->_coeffs[i] < lu) ? pb->_coeffs[i] : lu;
            }
            lu = (f_offset > lu) ? f_offset : lu;
            uint64_t i = lu, inxt, onxt;
            
            wlit_mapt::iterator it;
            while(i < rhs){
#ifdef __DEBUG__
                printf("c r_lb %lu\n", i);
#endif
                
                if((ndigit > MAXDIG) || vals.find(i-f_offset) != vals.end()){
                    less_or_equal_than(S, i);
                    inxt = epsilon * i;
                    
                    
                    inxt = (inxt > i+1) ? inxt : i+1;
                    
                    it = oliterals.lower_bound(inxt);
                    if(it != oliterals.begin())
                        it--;
                    onxt = it->first;
                    if(it != oliterals.end() && onxt > i && onxt < inxt)
                        i = onxt;
                    else
                        i = inxt;
                    
                }else{
                    i++;
                }
            }
            
        }
        
        order_enc_rootvars(S); //force order encoding in rootVars (this used to be done in Alg_BLS)
        
                
        lessclauses = nb_clauses - lessclauses;
        
        nobjvals = oliterals.size();
        
        printf("c n new clauses: %lu / %lu\n", nb_clauses-nb_clauses_prev, nb_clauses);
        printf("c n new less clauses: %lu / %lu\n", lessclauses, nb_clauses-nb_clauses_prev);
        printf("c different obj values: %d (max: %lu)\n", nobjvals, rhs);

    }else{ // incremental
        
        update(S, rhs);
    }
  }
}


Lit KPA::getValRootLit(Solver * S, uint64_t val){
//     printf("\nval: %ld\n", val);
    if(oliterals.find(val) != oliterals.end()) return oliterals.at(val);
    
    less_or_equal_than(S, val);
    
    Lit l = oliterals.at(val);
    wlit_mapt::iterator it = oliterals.lower_bound(val);
    if(it != oliterals.begin()){
        it--;
        Lit prev = it->second;
//         printf("prev: %ld\n", it->first);
        clause.push(~prev);
        clause.push(l);
        if(hasRelaxVar)
            clause.push(relax_var);
        S->addClause(clause);
        clause.clear();
    }
    
    
    if(oliterals.end()->first > val){
        Lit nxt = oliterals.upper_bound(val)->second;

//         printf("next: %ld\n", oliterals.upper_bound(val)->first);
        clause.push(~l);
        clause.push(nxt);
        if(hasRelaxVar)
            clause.push(relax_var);
        S->addClause(clause);
        clause.clear();
    }
    return l;
}

void KPA::order_enc_rootvars(Solver *S){
    wlit_mapt::iterator prev;
    size_t j=0;
    for (wlit_mapt::iterator rit = oliterals.begin(); rit != oliterals.end(); rit++, j++){
        #ifdef __DEBUG__
        printf("\t[%lu] ",rit->first);
        if (sign(rit->second))
            printf("~");
        printf("y%d\n", var(rit->second)+1);
#endif
        if(j >= 1){
            clause.push(~prev->second);
            clause.push(rit->second);
            if(hasRelaxVar)
                clause.push(relax_var);
            S->addClause(clause);
            clause.clear();
        }
        prev = rit;
            
    }
    
}


void KPA::feasibleObjectiveValuesDP(Solver * S, std::set<uint64_t> &vals, uint64_t sz){
    int n = pb_approx_coeffs->_lits.size();
    std::set<uint64_t>::reverse_iterator rit;
    
    vals.insert(0);
    
    
    uint64_t c;
    for(int i = 0; i < n; i++){
        std::set<uint64_t> new_vals;
        if(true){
            c = pb_approx_coeffs->_coeffs[i];
            for(rit = vals.rbegin(); rit != vals.rend(); rit++){
                if(*rit + c < sz){
                    new_vals.insert(*rit + c);
                }
            }
        }
        
        vals.insert(new_vals.begin(), new_vals.end());
    }
    
}

void KPA::update(Solver *S, uint64_t rhs) {

//   if(enc_version == _full_){
    printf("c [KP::update]\n");
    assert(hasEncoding);
    Lit o = less_or_equal_than(S, rhs);
    addUnitClause(S, o);
    nb_clauses++;
    //   S->my_print();
    current_pb_rhs = rhs;
//   }
    
    //TODO: No modo de aproximacao, ir buscar apenas a proxima classe
}



/* ---------------------------- converters ------------------------------------- */

//retorna o numero de variaveis na funcao simplificada (sem vars com valores fixos)
uint64_t KPA::convert_kp_sat4j(Solver * S, uint64_t ub){
//     printf("[KP::convert_kp_sat4j]\n");
    uint64_t ntmpvars = printPBfunction(S);
    if(ntmpvars == 1){
        encode_one_var_func(S);
    
    }else if(ntmpvars > 1){
        run_sat4j(ub);
        read_mocnf(S);
        
        force_sort_within_ls(S);
    }else{
        if(hasEncoding)
            printf("c approximated objective function remains the same!\n");
        else
            printf("c no need to convert, all vars have assigned values!\n");
    }
    
    return ntmpvars;
}

uint64_t KPA::convert_kp_minisatp(Solver * S, uint64_t ub){
//     printf("[KP::convert_kp_minisatp]\n");
    uint64_t ntmpvars = printPBfunction(S);
    
    if(ntmpvars > 0){
        run_kp_minisatp(ub);
        read_mocnf(S);
        
        force_sort_within_ls(S);
    }else{
        printf("c no need to convert, all vars have assigned values!\n");
    }
    
    return ntmpvars;
    
}
/* ------------------ converters - obj encoding ----------------------------------- */

void KPA::encode_one_var_func(Solver * S){
    printf("c just one free var\n");
        ncoeffs = 2;
        int i = pb->_coeffs.size() - 1;
        //for (i = 0; i < pb->_coeffs.size() && true; i++);
        uint64_t c = pb->_coeffs[i];
        int li;
        
        ls_szs = (uint64_t *) malloc(2 * sizeof(uint64_t));
        
        if(c > 1){
            li = 1;
            ls_szs[0] = 0;
            ls_szs[1] = 1;
        }else{
            li = 0;
            ls_szs[0] = 1;
            ls_szs[1] = 0;
        }
        ls = (int **) malloc(2 * sizeof(int *));
        ls[li] = (int *) malloc(1 * sizeof(int));
        ls[1-li] = NULL;
        
        base = (uint64_t *) malloc(2 * sizeof(uint64_t));
        mxb_coeffs = (uint64_t *) malloc(2 * sizeof(uint64_t));
        base[0] = mxb_coeffs[0] = 1;
        base[1] = mxb_coeffs[1] = c;
//         printf("ls: %d\n", ls[0][0]);
        
        if(sign(pb->_lits[i])){ //se sinal negativo, adicionar var para inverter sinal
            int t = get_new_Lit_id(S);
            ls[li][0] = t;
            vec<Lit> lits;
            
            lits.clear();
            //t => ~ls[0][0]
            if(relax_var != lit_Undef)
                lits.push(relax_var);
            lits.push(~mkLit(t)); lits.push(pb->_lits[i]); S->addClause(lits);
    
            //~ls[0][0] => t
    
            lits.clear();
            if(relax_var != lit_Undef)
                lits.push(relax_var);
            lits.push(mkLit(t)); lits.push(~(pb->_lits[i])); S->addClause(lits);
            
            ls[li][0] = t;
            nb_clauses = 2;
            nb_variables = 1;
//             printf("ls: %d\n", ls[0][0]);
        }else{
            ls[li][0] = var(pb->_lits[i]);
        }
}


void KPA::force_sort_within_ls(Solver * S){
    for(int i = 0; i < ncoeffs; i++){
        for(int j = ls_szs[i]-1; j >= 1; j--){
            addBinaryClause(S, ~mkLit(ls[i][j]), mkLit(ls[i][j-1]), relax_var); // l_i_j => l_i_{j-1}
        }
    }
}


Lit KPA::less_or_equal_than(Solver * S, uint64_t rhs){
    
    if(oliterals.find(rhs) != oliterals.end()) return oliterals.at(rhs);

    uint64_t v = rhs, r;
    
    Lit orhs = getNewLit(S);
    Lit o = orhs;
    oliterals[v] = o;
    
    uint64_t cumsum[ncoeffs];
    cumsum[0] = ls_szs[0];
    for(int i = 1; i < ncoeffs; i++){
        cumsum[i] = cumsum[i-1] + ls_szs[i]*mxb_coeffs[i];
    }

    v = rhs = (rhs > f_offset) ? rhs - f_offset : 0; // o valor minimo da funcao (consequencia das variaveis fixas) e' f_offset, portanto so se pode restringir usando as restantes variaveis, indicando que teem de ser <= que rhs - f_offset 
    
    Lit o2;
    bool stop;
    
    for(int i = ncoeffs-1; i >= 0; i--){
        if(ls_szs[i] > 0){
            r = v / mxb_coeffs[i];
//             printf("mxb_coeffs[%d]: %lu\n", i, mxb_coeffs[i]);
            
            
            if(r == 0 || i == 0 || v % mxb_coeffs[i] == cumsum[i-1] || cumsum[i] == 0){
                if(ls_szs[i] > r){
                    addBinaryClause(S, ~o, ~mkLit(ls[i][r]), relax_var);
//                     printf("(%d => ~%d)\n", var(o), ls[i][r]);
                    nb_clauses++;
                    
                }
                if(i > 0 && v % mxb_coeffs[i] == cumsum[i-1]){
//                     printf("%lu == %lu\n", v % mxb_coeffs[i], cumsum[i-1]);
                    break;
                }
            }else{
                v = v % mxb_coeffs[i];
                if(r-1 < ls_szs[i]){
                    stop = olits_aux.find(std::make_pair(v, i)) != olits_aux.end();
                    if(!stop){
                        o2 = getNewLit(S);
                        olits_aux[std::make_pair(v, i)] = o2;
                    }else{
                        o2 = olits_aux.at(std::make_pair(v, i));
                    }
                    
                    if(r < ls_szs[i]){
                        addBinaryClause(S, ~o, ~mkLit(ls[i][r]), relax_var);
//                         printf("(%d => ~%d)\n", var(o), ls[i][r]);
                        nb_clauses++;
                    }
                    
                
                    addTernaryClause(S, ~o, ~mkLit(ls[i][r-1]), o2, relax_var);
                    nb_clauses++;
//                     printf("(%d => ~%d \\/ %d)\n", var(o), ls[i][r-1], var(o2));
                        
                    if(stop)
                        break;
                    o = o2;
                }else{
                    break;
                }
            }
        }
    }
    
//     printf("[leq:ok] %lu (%d)\n\n", rhs, var(orhs));
    return orhs;
}


/* -------------------------------------------------------------------------------- */
/* ---------------------------- I/O functions ------------------------------------- */

//retorna o numero de variaveis na funcao simplificada (sem vars com valores fixos)
uint64_t KPA::printPBfunction(Solver *S){
    
    char hostname[HOST_NAME_MAX];
    gethostname(hostname, HOST_NAME_MAX);
    sprintf(ifname, "/tmp/in-%s-%d.opb", hostname, ::getpid());
    
#ifdef __DEBUG__
    printf("c [KP::printPBfunction] %s\n", ifname);
#endif
    FILE * f = fopen(ifname, "w");
    fprintf(f, "min:");
    
    f_offset = 0;
    
    name2id.clear();
    id2sign.clear();
    
    char tmp[20];
    
    int j = 0;
    
    uint64_t c, tmpc;
    
    bool f_changed = false;
    PB * old_pb_approx_coeffs = pb_approx_coeffs;
    
//     printf("c Approx mode: %d\n", approxMode);
    
    bool useOriginalPB = (epsilon == 1 || approxMode == _ap_outvars_);
    
    if(useOriginalPB)
        pb_approx_coeffs = pb;
    else
        pb_approx_coeffs = new PB();
    
//     printf("c useOriginalPB? %d\n", useOriginalPB);
        
    
    for (int i = 0; i < pb->_coeffs.size(); i++) {
//         printf("%lu ", pb->_coeffs[i]);
//             fprintf(f, " %c%lu %s",((pb->_coeffs[i] >= 0) ? '+':'-'), pb->_coeffs[i], tmp);
//             fprintf(f, " %c%lu %s",(sign(pb->_lits[i]) ? '+':'-'), pb->_coeffs[i], tmp);
        if(true){
            sprintf(tmp, "x%d", j+1);
            
            if(useOriginalPB){
                fprintf(f, " +%lu %s", pb->_coeffs[i], tmp);
            
            }else{                
                //approximation mode where PB coeffs are rounded down
                c = pb->_coeffs[i];
                tmpc = 1;
                while(tmpc <= pb->_coeffs[i]){
                    c = tmpc;
                    tmpc = (tmpc*epsilon > tmpc+1) ? round(tmpc*epsilon) : tmpc+1;
//                     printf("c c %lu\n", c);
                }
                fprintf(f, " +%lu %s", c, tmp);
                pb_approx_coeffs->addProduct(pb->_lits[i], c);
                
                if(hasEncoding && old_pb_approx_coeffs->_coeffs[i] != pb_approx_coeffs->_coeffs[i])
                    f_changed = true;
            }
            
            name2id[tmp] = var(pb->_lits[i]); 
            id2sign[var(pb->_lits[i])] = !sign(pb->_lits[i]);
            j++;
            
        }else{
            if(!useOriginalPB){
                pb_approx_coeffs->addProduct(pb->_lits[i], pb->_coeffs[i]);
                
                if(hasEncoding && old_pb_approx_coeffs->_coeffs[i] != pb_approx_coeffs->_coeffs[i])
                        f_changed = true;
            }
            
            if(S->value(pb->_lits[i]) == l_True && !sign(pb->_lits[i])){
                f_offset += pb->_coeffs[i];
            }else if(S->value(pb->_lits[i]) == l_True && sign(pb->_lits[i])){
                f_offset += pb->_coeffs[i];
            }
        }
//             printf("%cx%d\n", (!sign(pb->_lits[i])) ? '+':'-', i+1);
//             printf("sv: "); S->printLit(pb->_lits[i]); printf("\n");
//             printf("%s ", maxsat_formula->getIndexToName().at(var(pb->_lits[i])).c_str());

    }
    
    
    fprintf(f, ";\n");
    fclose(f);
    
#ifdef __DEBUG__
    printf("PBo: "); pb->print();
    printf("PBa: "); pb_approx_coeffs->print();
    
    printf("c vars_fixas: %d\n", pb->_coeffs.size() - j);
    printf("c vars_free: %d\n", j);
    /*
    if(j == 0){
        printf("c funcao objectivo vazia!\n");
        exit(-1);
    }*/
#endif
    
    printf("c offset: %lu\n", f_offset);
    
    if(approxMode == _ap_coeffs_ && hasEncoding){
        if (!f_changed && old_pb_approx_coeffs->_coeffs.size() == pb_approx_coeffs->_coeffs.size()){
            return 0;
        }else{
            forgetEncoding(S);
        }
    }
    
    return j;
}

void KPA::run_sat4j(uint64_t ub){
    printf("\nc [KP::run_sat4j]\n");
    
        
    char hostname[HOST_NAME_MAX];
    gethostname(hostname, HOST_NAME_MAX);
//     sprintf(ifname, "/tmp/in-%s-%d.opb", hostname, ::getpid());
    
    sprintf(ofname, "/tmp/out-%s-%d.mocnf", hostname, ::getpid());

    const char* cmdbase = "java -jar %s/%s -ib 1 -ul 0:%d %s -o %s";
//      cmdbase = "java -jar  ./../../joaoGitHub/moco/target/org.sat4j.moco.threeAlgorithms-0.0.1-SNAPSHOT-jar-with-dependencies.jar -ib 1 -ul 0:";

    
    char cmd[strlen(cmdbase) + 4 * max::name + context::max_cwd + 100];
    snprintf(cmd, strlen(cmdbase) + 4 * max::name + context::max_cwd, cmdbase, context::cwd, SAT4J_ENC, ub, ifname, ofname);
    
    printf("c [run_sat4j] cmd: %s\n", cmd);
    
    system(cmd);
    printf("c [run_sat4j] finished\n");
}


//TODO: Fix in and out file names
void KPA::run_kp_minisatp(uint64_t ub){
    printf("[KP::run_kp_minisatp]\n");
    
    char st[200];
    sprintf(st, "./../../KP-minisat/kp-minisatp-master/build/release/bin/minisatp -gs in.opb -goal=%lu -mocnf=out.mocnf", ub);
    printf("run: %s\n", st);
    
    system(st);
}


void KPA::read_mocnf(Solver *S){
    printf("c [KP::read_mocnf] %s\n", ofname);
    tmp_kp_t tmp;
    
//     char * fname = "out.mocnf";
    parseMOCNF(ofname, S, &tmp, name2id, id2sign, relax_var);
    
    nbase = tmp.nbase;
    ncoeffs = tmp.ncoeffs;
    base = tmp.base;
    mxb_coeffs = tmp.mxb_coeffs;
    
    ls = tmp.ls;
    ls_szs = tmp.ls_szs;
    
    if(tmp.haszs)
        zs = tmp.zs;
    else
        zs = NULL;
    
    nb_clauses = tmp.nclauses;
    nb_variables = tmp.nvars;
    
    //clean files
    char cmd[2*max::name + 5];
//     sprintf(cmd, "cp %s /tmp/in.opb; cp %s /tmp/out.mocnf; rm %s %s", ifname, ofname, ifname, ofname);
    sprintf(cmd, "rm %s %s", ifname, ofname);

    printf("c cmd: %s\n", cmd);
    system(cmd);
    
    //save ls/zs
    //TODO
    
}


//assume que o novo eps e inferior ao anterior
void KPA::setApproxRatio(Solver *S, float eps, int mode){
    assert(eps >= 1);
    epsilon = eps;
    
    approxMode = mode;
//     printf("c hasEncoding %d\n", hasEncoding);
    if(hasEncoding){
        //TODO
    }
    
    if(approxMode == _ap_coeffs_){
        if(!hasRelaxVar && epsilon > 1){
            assumptions.clear();
            relax_var = getNewLit(S);
            hasRelaxVar = true;
            assumptions.push_back(~relax_var);
        }
    }
    
    
    
}

void KPA::forgetEncoding(Solver *S){
    printf("c forget previous encoding\n");
    if(approxMode == _ap_coeffs_){
        //TODO
        if(hasRelaxVar){
            //clean encoding (permanently turn off relax var)
            addUnitClause(S, relax_var);
            oliterals.clear();
            olits_aux.clear();
    
        }
        
        assumptions.clear();
            
        //create new relax var
        if(epsilon > 1){
            relax_var = getNewLit(S);
            hasRelaxVar = true;
            assumptions.push_back(~relax_var);
        }else{
            relax_var = lit_Undef;
            hasRelaxVar = false;
        }
    }
}
    
    
    
/* currentModel - evaluate this model
 * objv - real evaluation (according to original PB)
 * ap_objv - approximated evaluation (according to relaxed PB for approximation mode _ap_coeffs_,
 * in any other mode objv = ap_objv)
 */
void KPA::evalModel(vec<lbool> &currentModel, uint64_t * objv, uint64_t *ap_objv, int di){
    Lit l;
    objv[di] = 0;
    ap_objv[di] = 0;
    for(int i = 0; i < pb->_lits.size(); i++){
        l = pb->_lits[i];
        if ((sign(l) && currentModel[var(l)] == l_False) || 
            (!sign(l) && currentModel[var(l)] == l_True)) {
            objv[di] += pb->_coeffs[i];
            ap_objv[di] += pb_approx_coeffs->_coeffs[i];
//                 printf(" %cx%d", (sign(l)) ? '-':' ', var(l)+1);
        }
    }
//         printf("\n");
}




bool KPA::apFunctionIsExact(){
    
    for(int i = 0; i < pb->_lits.size(); i++){
        if (pb->_coeffs[i] != pb_approx_coeffs->_coeffs[i]){
            return false;
        }
    }
    return true;
}
