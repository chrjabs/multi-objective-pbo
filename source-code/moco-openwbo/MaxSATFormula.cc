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

#include <iostream>
#include <vector>
#include <algorithm>

#include "MaxSATFormula.h"

using namespace openwbo;

MaxSATFormula *MaxSATFormula::copyMaxSATFormula() {
  assert(format == _FORMAT_MAXSAT_);

  MaxSATFormula *copymx = new MaxSATFormula();
  copymx->setInitialVars(nVars());

  for (int i = 0; i < nVars(); i++)
    copymx->newVar();

  for (int i = 0; i < nSoft(); i++)
    copymx->addSoftClause(getSoftClause(i).weight, getSoftClause(i).clause);

  for (int i = 0; i < nHard(); i++)
    copymx->addHardClause(getHardClause(i).clause);

  copymx->setProblemType(getProblemType());
  copymx->updateSumWeights(getSumWeights());
  copymx->setMaximumWeight(getMaximumWeight());
  copymx->setHardWeight(getHardWeight());

  return copymx;
}

// Adds a new hard clause to the hard clause database.
void MaxSATFormula::addHardClause(vec<Lit> &lits) {
  hard_clauses.push();
  vec<Lit> copy_lits;
  lits.copyTo(copy_lits);
  new (&hard_clauses[hard_clauses.size() - 1]) Hard(copy_lits);
  n_hard++;
}

// Adds a new soft clause to the hard clause database.
void MaxSATFormula::addSoftClause(uint64_t weight, vec<Lit> &lits) {
  soft_clauses.push();
  vec<Lit> vars;
  Lit assump = lit_Undef;
  vec<Lit> copy_lits;
  lits.copyTo(copy_lits);

  new (&soft_clauses[soft_clauses.size() - 1])
      Soft(copy_lits, weight, assump, vars);
  n_soft++;
}

// Adds a new soft clause to the hard clause database with predefined relaxation
// variables.
void MaxSATFormula::addSoftClause(uint64_t weight, vec<Lit> &lits,
                                  vec<Lit> &vars) {
  soft_clauses.push();
  Lit assump = lit_Undef;
  vec<Lit> copy_lits;
  lits.copyTo(copy_lits);

  new (&soft_clauses[soft_clauses.size() - 1])
      Soft(copy_lits, weight, assump, vars);
  n_soft++;
}

void MaxSATFormula::addSoftClause(uint64_t weight, vec<Lit> &lits, size_t obj_idx) {
  while (clausal_objectives.size() <= obj_idx) {
    clausal_objectives.push();
    offsets.push(0);
    ++n_objf;
  }
  clausal_objectives[obj_idx].push();
  vec<Lit> vars;
  Lit assump = lit_Undef;
  vec<Lit> copy_lits;
  lits.copyTo(copy_lits);

  new (&clausal_objectives[obj_idx][clausal_objectives[obj_idx].size() - 1])
      Soft(copy_lits, weight, assump, vars);
}

int MaxSATFormula::nInitialVars() {
  return n_initial_vars;
} // Returns the number of variables in the working MaxSAT formula.

void MaxSATFormula::setInitialVars(int vars) { n_initial_vars = vars; }

int MaxSATFormula::nVars() {
  return n_vars;
} // Returns the number of variables in the working MaxSAT formula.

int MaxSATFormula::nSoft() {
  return n_soft;
} // Returns the number of soft clauses in the working MaxSAT formula.

int MaxSATFormula::nSoft(uint obj) {
  return clausal_objectives[obj].size();
}

int MaxSATFormula::nHard() {
  return n_hard;
} // Returns the number of hard clauses in the working MaxSAT formula.

void MaxSATFormula::newVar(int v) {
  if(v == -1) n_vars++;
  else if(v > n_vars) n_vars = v;
} // Increases the number of variables in the working MaxSAT formula.

// Makes a new literal to be used in the working MaxSAT formula.
Lit MaxSATFormula::newLiteral(bool sign) {
  Lit p = mkLit(nVars(), sign);
  newVar();
  return p;
}

void MaxSATFormula::setProblemType(int type) {
  problem_type = type;
} // Sets the problem type.

int MaxSATFormula::getProblemType() {
  return problem_type; // Return the problem type.
}

// 'ubCost' is initialized to the sum of weights of the soft clauses.
void MaxSATFormula::updateSumWeights(uint64_t weight) {
  if (weight != hard_weight)
    sum_soft_weight += weight;
}

// The initial 'currentWeight' corresponds to the maximum weight of the soft
// clauses.
void MaxSATFormula::setMaximumWeight(uint64_t weight) {
  if (weight > max_soft_weight && weight != hard_weight)
    max_soft_weight = weight;
}

uint64_t MaxSATFormula::getMaximumWeight() { return max_soft_weight; }

void MaxSATFormula::setHardWeight(uint64_t weight) {
  hard_weight = weight;
} // Sets the weight of hard clauses.

Soft &MaxSATFormula::getSoftClause(int pos) {
  assert(pos < nSoft());
  return soft_clauses[pos];
}

Soft &MaxSATFormula::getSoftClause(uint obj, int pos) {
  assert(pos < nSoft(obj));
  return clausal_objectives[obj][pos];
}

Hard &MaxSATFormula::getHardClause(int pos) {
  assert(pos < nHard());
  return hard_clauses[pos];
}

void MaxSATFormula::addPBConstraint(PB *p) {

  // Add constraint to formula data structure.
  if (p->isClause()) {
    addHardClause(p->_lits);
  } else if (p->isCardinality()) {
    if (!p->_sign) {
      p->changeSign();
    }
    cardinality_constraints.push(new Card(p->_lits, p->_rhs));

  } else {
    if (!p->_sign) {
      p->changeSign();
    }

    pb_constraints.push(new PB(p->_lits, p->_coeffs, p->_rhs, p->_sign));
  }
}

int MaxSATFormula::newVarName(char *varName) {
  int id = varID(varName);
  if (id == var_Undef) {
    id = nVars();
    newVar();
    std::string s(varName);
    std::pair<std::string, int> nv(s, id);
    std::pair<int, std::string> ni(id, s);
    _nameToIndex.insert(nv);
    _indexToName.insert(ni);
  }
  return id;
}

int MaxSATFormula::varID(char *varName) {
  std::string s(varName);

  nameMap::const_iterator iter = _nameToIndex.find(s);
  if (iter != _nameToIndex.end()) {
    return iter->second;
  }
  return var_Undef;
}

void MaxSATFormula::convertPBtoMaxSAT() {
// void MaxSATFormula::convertPBtoMaxSAT(Solver * S = NULL) { //AG
    
//   printf("MaxSATFormula::convertPBtoMaxSAT n_objf: %d\n", n_objf);
  assert(objective_functions != NULL);
  vec<Lit> unit_soft(1);

  if (n_objf == 1) {
    // Convert objective function to soft clauses
    for (int i = 0; i < objective_functions[0]->_lits.size(); i++) {
        assert(objective_functions[0]->_coeffs[i] > 0);
        unit_soft[0] = ~objective_functions[0]->_lits[i];
        addSoftClause(objective_functions[0]->_coeffs[i], unit_soft);

        // Updates the maximum weight of soft clauses.
        setMaximumWeight(objective_functions[0]->_coeffs[i]);
        // Updates the sum of the weights of soft clauses.
        updateSumWeights(objective_functions[0]->_coeffs[i]);
    }
  }else{ //AG - n_objf > 1 (Multiobjective case)
      //TODO: converter fobj para GTE
      // objective_function->_lits[i] ou ~objective_function->_lits[i]?
      for(int j = 0; j < n_objf; j++){
          int totalw = 0;
            for (int i = 0; i < objective_functions[j]->_lits.size(); i++, totalw += objective_functions[j]->_coeffs[i]);
            
//TODO: Fazer versao standalone to GTE (recebe pbc e devolve conjunto de clausulas e literais
//             na raiz)
//         encode(this, objective_functions[j]->_lits, objective_functions[j]->_coeffs, totalw);
      }
      //      obter os literais do no raiz, R

      //      usar o order encoding nos literais de R (adicionar hard clauses)
      
      //      adicionar os literais em R como soft clauses
      
  }
  
  
  if (getMaximumWeight() == 1)
    setProblemType(_UNWEIGHTED_);
  else
    setProblemType(_WEIGHTED_);
  
  
  //AG
//   printf("MaxSATFormula::convertPBtoMaxSAT -- formula after conversion\n");
//   my_print();
  
}

void MaxSATFormula::clausalToPbObjs() {
  assert(!objective_functions.size());
  for (int obj_idx = 0; obj_idx < clausal_objectives.size(); obj_idx++) {
    objective_functions.push(new PBObjFunction());
    for (int cl_idx = 0; cl_idx < clausal_objectives[obj_idx].size(); cl_idx++) {
      Soft &soft = clausal_objectives[obj_idx][cl_idx];
      if (1 == soft.clause.size()) {
        objective_functions[obj_idx]->addProduct(~soft.clause[0], soft.weight);
      } else {
        // Add clause as hard with relaxation variable
        if (lit_Undef == soft.assumption_var) {
          soft.assumption_var = newLiteral();
        }
        soft.clause.push(soft.assumption_var);
        addHardClause(soft.clause);
        objective_functions[obj_idx]->addProduct(soft.assumption_var, soft.weight);
      }
    }
    clausal_objectives[obj_idx].clear();
    objective_functions[obj_idx]->_const = offsets[obj_idx];
  }
  clausal_objectives.clear();
  offsets.clear();
}


void Hard:: my_print(indexMap indexToName){   
    printf("c\t");
    for(int i = 0; i < clause.size(); i++)
//         printf(" %s%d", ((sign(clause[i])) ? "~" : ""), var(clause[i]));
        if(indexToName.find(var(clause[i])) != indexToName.end())
            printf(" %s%s", ((sign(clause[i])) ? "~" : ""), indexToName.at(var(clause[i])).c_str());
        else
            printf(" %sx%d", ((sign(clause[i])) ? "~" : ""), var(clause[i]));
    printf("\n");
}

void Soft:: my_print(indexMap indexToName){   
    printf("c\t");
    for(int i = 0; i < clause.size(); i++)
//         printf(" %s%d", ((sign(clause[i])) ? "~" : ""), var(clause[i]));
        
        if(indexToName.find(var(clause[i])) != indexToName.end())
            printf(" %lu %s%s", weight, ((sign(clause[i])) ? "~" : ""), indexToName.at(var(clause[i])).c_str());
        else
            printf(" %sx%d", ((sign(clause[i])) ? "~" : ""), var(clause[i]));
    printf(" [");
    for(int i = 0; i < relaxation_vars.size(); i++)
        printf(" %s%d", ((sign(relaxation_vars[i])) ? "~" : ""), var(relaxation_vars[i]));
    printf("]");
    printf(" A(x%d)\n", var(assumption_var) + 1);
}



void MaxSATFormula::my_print(){
    printf("c -----------------------------------\n");
    printf("c --------- index to var ------------\n");
//     int k;
//     string v;
    /*
    for (auto const& pair: _indexToName) {
        std::cout << "c\t{" << pair.first+1 << ": " << pair.second << "}\n";
    }*/
    
//     int maxsize = 20;
    int maxcs = 10; //100
    printf("c ------- Objective Func (%d) -----------------\n", n_objf);
    for(int di = 0; di < objective_functions.size(); di++)
        objective_functions[di]->my_print(_indexToName);
    for (int obj_idx = 0; obj_idx < clausal_objectives.size(); obj_idx++) {
      printf("c objective %d\n", obj_idx);
      for (int cl_idx = 0; cl_idx < std::min(clausal_objectives[obj_idx].size(), maxcs); cl_idx++)
        clausal_objectives[obj_idx][cl_idx].my_print(_indexToName);
      if (clausal_objectives[obj_idx].size() > maxcs)
        printf("c ...\n");
    }
    
    printf("c ------- PB Constraints (%d) ------------\n", pb_constraints.size());
    for(int i = 0; i < std::min(pb_constraints.size(), maxcs); i++)
        pb_constraints[i]->my_print(_indexToName);
    if(pb_constraints.size() > maxcs)
      printf("c ...\n");
    
    printf("c ------- Cardinality Constr (%d) --------\n", cardinality_constraints.size());
    for(int i = 0; i < std::min(cardinality_constraints.size(), maxcs); i++)
        cardinality_constraints[i]->my_print(_indexToName);
    if(cardinality_constraints.size() > maxcs)
      printf("c ...\n");

    printf("c ---------------------------------------\n");
    printf("c ------- Hard clauses (%d) --------------\n", hard_clauses.size());
    for(int i = 0; i < std::min(hard_clauses.size(), maxcs); i++)
        hard_clauses[i].my_print(_indexToName);
    if(hard_clauses.size() > maxcs)
      printf("c ...\n");
    printf("c ------- Soft clauses (%d) --------------\n", soft_clauses.size());
    for(int i = 0; i < std::min(soft_clauses.size(), maxcs); i++)
        soft_clauses[i].my_print(_indexToName);
    if(soft_clauses.size() > maxcs)
      printf("c ...\n");
    
    
    
    printf("c ---------------------------------------\n");
    
}

void MaxSATFormula::sync_first(Glucose::Solver *s){
  for(int i = 0; i < nInitialVars(); i++)
    if(s->value(i) != l_Undef)
      fv.insert(mkLit(i,true));
}

#define IPASIR(l) (sign(l) ? -(var(l) + 1) : (var(l) + 1))
#define MSLIT(l) (mkLit(abs(l) - 1, l < 0))

void MaxSATFormula::preprocess(const char *tech) {
  assert(nullptr == prepro);
  assert(_FORMAT_MCNF_ == format);
  assert(0 == soft_clauses.size());
  assert(0 == objective_functions.size());

  std::vector<std::vector<int>> clauses{};
  std::vector<std::vector<uint64_t>> weights{};

  for (int cl_idx = 0; cl_idx < nHard(); cl_idx++) {
    Hard &cl = getHardClause(cl_idx);
    std::vector<int> ipasir_cl{};
    for (int l_idx = 0; l_idx < cl.clause.size(); l_idx++) {
      ipasir_cl.push_back(IPASIR(cl.clause[l_idx]));
    }
    clauses.push_back(ipasir_cl);
    weights.push_back({});
  }

  uint64_t top = 0;
  for (int obj_idx = 0; obj_idx < nObjFunctions(); obj_idx++) {
    for (int cl_idx = 0; cl_idx < nSoft(obj_idx); cl_idx++) {
      Soft &cl = getSoftClause(obj_idx, cl_idx);
      std::vector<int> ipasir_cl{};
      std::vector<uint64_t> weight{};
      for (int l_idx = 0; l_idx < cl.clause.size(); l_idx++) {
        ipasir_cl.push_back(IPASIR(cl.clause[l_idx]));
      }
      clauses.push_back(ipasir_cl);
      top += cl.weight;
      for (int i = 0; i < obj_idx; i++) {
        weight.push_back(0);
      }
      weight.push_back(cl.weight);
      weights.push_back(weight);
    }
  }

  prepro = new maxPreprocessor::PreprocessorInterface(clauses, weights, top);
  prepro->preprocess(tech);
  clauses.clear();
  weights.clear();
  std::vector<int> labels{};
  prepro->getInstance(clauses, weights, labels);
  top = prepro->getTopWeight();

  hard_clauses.clear();
  n_hard = 0;
  clausal_objectives.clear();
  n_objf = 0;
  n_vars = 0;

  for (size_t cl_idx = 0; cl_idx < clauses.size(); cl_idx++) {
    std::vector<int> &cl = clauses[cl_idx];
    std::vector<uint64_t> &ws = weights[cl_idx];

    vec<Lit> ms_cl{};
    for (auto l : cl) {
      if (abs(l) > nVars())
        newVar(abs(l));
      ms_cl.push(MSLIT(l));
    }

    bool is_hard = true;
    for (size_t obj_idx = 0; obj_idx < ws.size(); obj_idx++) {
      if (ws[obj_idx] == top)
        continue;
      if (ws[obj_idx] > 0) {
        addSoftClause(ws[obj_idx], ms_cl, obj_idx);
      }
      is_hard = false;
    }
    if (is_hard)
      addHardClause(ms_cl);
  }
  
  std::vector<uint64_t> removed_weight = prepro->getRemovedWeight();
  for (size_t obj_idx = 0; obj_idx < removed_weight.size(); obj_idx++) {
    offsets[obj_idx] = removed_weight[obj_idx];
  }
}

void MaxSATFormula::reconstruct(vec<lbool> &model) {
  if (nullptr == prepro) return;
  
  std::vector<int> true_lits{};
  for (int var_idx = 0; var_idx < model.size(); var_idx++) {
    if (model[var_idx] == l_True) {
      true_lits.push_back(var_idx + 1);
    } else if (model[var_idx] == l_False) {
      true_lits.push_back(-(var_idx + 1));
    }
  }
  
  true_lits = prepro->reconstruct(true_lits);
  model.clear();
  model.growTo(prepro->getOriginalVariables(), l_Undef);
  
  for (auto l : true_lits) {
    model[abs(l) - 1] = l < 0 ? l_False : l_True;
  }
}

void MaxSATFormula::extractObjVars() {
  obj_vars.clear();
  for (uint obj_idx = 0; obj_idx < objective_functions.size(); obj_idx++) {
    PBObjFunction *pb_fnc = objective_functions[obj_idx];
    for (uint lit_idx = 0; lit_idx < pb_fnc->_lits.size(); lit_idx++) {
      obj_vars.push_back(var(pb_fnc->_lits[lit_idx]));
    }
  }
  std::sort(obj_vars.begin(), obj_vars.end());
  auto end = std::unique(obj_vars.begin(), obj_vars.end());
  obj_vars.resize(std::distance(obj_vars.begin(), end));
}
