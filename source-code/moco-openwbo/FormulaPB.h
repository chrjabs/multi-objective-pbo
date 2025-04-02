/*!
 * \author Vasco Manquinho - vmm@sat.inesc-id.pt
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

#ifndef FormulaPB_h
#define FormulaPB_h

#ifdef SIMP
#include "simp/SimpSolver.h"
#else
#include "core/Solver.h"
#endif
#include <map>
#include <memory>
using NSPACE::vec;
using NSPACE::Lit;

namespace openwbo {

    
typedef std::map<int, std::string> indexMap;
    
// Cardinality constraint of the form atMostK
class Card {

public:
  Card(vec<Lit> &lits, int64_t rhs, bool sign = false) {
    lits.copyTo(_lits);
    _rhs = rhs;
    if (sign) {
      int s = 0;
      for (int i = 0; i < _lits.size(); i++) {
        s += 1;
        _lits[i] = ~_lits[i];
      }
      _rhs = s - _rhs;
    }
//     print();
  }

  Card() { _rhs = 0; }
  ~Card() {}

  void print() {
    printf("Card: ");

    for (int i = 0; i < _lits.size(); i++) {
      if (sign(_lits[i]))
        printf("~");
      printf("%d ", var(_lits[i]) + 1);
    }
    printf(" <= %d\n", (int)_rhs);
  }
  
  
    void my_print(indexMap indexToName, bool original_vars = true) {
    // Assume _sign == false...
    printf("Card: ");

    for (int i = 0; i < _lits.size(); i++) {
      if (sign(_lits[i]))
        printf("~");
      if(!original_vars)
        printf("X%d ", var(_lits[i]) + 1);
      else
          if(indexToName.find(var(_lits[i])) != indexToName.end())
            printf("%s ", indexToName.at(var(_lits[i])).c_str());
          else
            printf("X%d ", var(_lits[i]) + 1);
    }
    printf(" <= %d\n", (int)_rhs);
  }

  vec<Lit> _lits;
  int64_t _rhs;
};

// PB constraint. The constraint sign is encoded in the structure.
class PB {

public:
  PB(vec<Lit> &lits, vec<uint64_t> &coeffs, int64_t rhs, bool s = false) {
    lits.copyTo(_lits);
    coeffs.copyTo(_coeffs);
    _rhs = rhs;
    _sign = s;
    for(int i = 0; i< coeffs.size(); i++) ub += coeffs[i];
    
//     for(int i = 0; i< coeffs.size(); i++) printf("%lu\n", coeffs[i]);
//     printf("%d\n" % ((s) ? 1 : 0));
//     printf("rhs: %lu\n", rhs);
  }

  PB() {
    _rhs = 0;
    _sign = false;
    ub = 0;
  }
  ~PB() {}

  void addProduct(Lit l, int64_t c) {
    _coeffs.push();
    _lits.push();
    if (c >= 0) {
      _coeffs[_coeffs.size() - 1] = c;
      _lits[_lits.size() - 1] = l;
      ub += c;
    } else {
      _coeffs[_coeffs.size() - 1] = -c;
      _lits[_lits.size() - 1] = ~l;
      _rhs += -c;
      ub += -c;
    }
  }

  void addRHS(int64_t rhs) { _rhs += rhs; }

  void changeSign() {
    int s = 0;
    for (int i = 0; i < _coeffs.size(); i++) {
      s += _coeffs[i];
      _lits[i] = ~(_lits[i]);
    }
    _rhs = s - _rhs;
    _sign = !(_sign);
  }

  bool isClause() {
    // Assume _sign == false...
    bool sign = _sign;
    if (_sign)
      changeSign();
    if (_rhs != 1) {
      if (_sign != sign)
        changeSign();
      return false;
    }
    for (int i = 0; i < _coeffs.size(); i++) {
      if (_coeffs[i] != 1) {
        if (_sign != sign)
          changeSign();
        return false;
      }
    }
    return true;
  }

  bool isCardinality() {
    // Assume _sign == false...
    bool sign = _sign;
    if (_sign)
      changeSign();
    for (int i = 0; i < _coeffs.size(); i++) {
      if (_coeffs[i] != 1) {
        if (_sign != sign)
          changeSign();
        return false;
      }
    }
    return true;
  }

  void print() {
    // Assume _sign == false...
    if (isClause())
      printf("Clause: ");
    else if (isCardinality())
      printf("Card: ");
    else
      printf("PB: ");

    for (int i = 0; i < _coeffs.size(); i++) {
      printf("%d ", (int)_coeffs[i]);
      if (sign(_lits[i]))
        printf("~");
      printf("%d ", var(_lits[i]) + 1);
    }
    if(!_sign)
        printf(" >= %d\n", (int)_rhs);
    else
        printf(" <= %d\n", (int)_rhs);
  }

  void my_print(indexMap indexToName, bool original_vars = true) {
    // Assume _sign == false...
    if (isClause())
      printf("c\tClause:\t");
    else if (isCardinality())
      printf("c\tCard:\t");
    else
      printf("c\tPB:\t");

    for (int i = 0; i < _coeffs.size(); i++) {
      printf("%d ", (int)_coeffs[i]);
      if (sign(_lits[i]))
        printf("~");
      if(!original_vars)
        printf("y%d ", var(_lits[i]) + 1);
      else
        printf("%s ", indexToName.at(var(_lits[i])).c_str());
    }
    
    
    if(!_sign)
        printf(" >= %d\n", (int)_rhs);
    else
        printf(" <= %d\n", (int)_rhs);
  }
  
  
  uint64_t getUB(){ return ub;}
  
  vec<uint64_t> _coeffs;
  vec<Lit> _lits;
  int64_t _rhs;
  bool _sign; // atLeast: false; atMost: true
  uint64_t ub;
};

class PBObjFunction {

    public:
    PBObjFunction(vec<Lit> &lits, vec<uint64_t> &coeffs, int64_t c = 0) {
        lits.copyTo(_lits);
        coeffs.copyTo(_coeffs);
        _const = c;
        ub = 0;
	lb = 0;
	total = 0;
        for(int i = 0; i < coeffs.size(); i++){
	  if(coeffs[i] > 0) ub += coeffs[i];
	  else lb +=coeffs[i];
	  }
	total = ub + lb;
    }
    PBObjFunction(PBObjFunction&& other){
      _lits = std::move(other._lits);
      _coeffs = std::move(other._coeffs);
      ub = other.ub;
      lb = other.lb;
      total = other.total;
      _const = other._const;
      other.ub = 0;
      other.lb = 0;
      other._const = 0;
      other.total = 0;
    }
  
  PBObjFunction& operator=(PBObjFunction&& other){
    _lits.clear();
    _coeffs.clear();
    _lits = std::move(other._lits);
    _coeffs = std::move(other._coeffs);
    ub = other.ub;
    lb= other.lb;
    _const = other._const;
    other.ub = 0;
    other.lb = 0;
    other._const = 0;
    other.total = 0;
    return *this;
  }
    PBObjFunction(const PBObjFunction& other){
        other._lits.copyTo(_lits);
        other._coeffs.copyTo(_coeffs);
	ub = other.ub;
	lb = other.lb;
	total = other.total;
	_const = other._const;
    }


  PBObjFunction() { _const = 0; ub = 0;lb = 0; total = 0;}
    ~PBObjFunction() {}

    void addProduct(Lit l, int64_t c) {
        _coeffs.push();
        _lits.push();
        if (c >= 0) {
        _coeffs[_coeffs.size() - 1] = c;
	ub += c;
	total += c;
        _lits[_lits.size() - 1] = l;
        } else {
        _coeffs[_coeffs.size() - 1] = -c;
        _lits[_coeffs.size() - 1] = ~l;
	lb += -c;
	total += -c;
        _const += c;
        }
//         printf("update _const: %ld\n", _const);
    }

    vec<uint64_t> _coeffs;
    vec<Lit> _lits;
    int64_t _const;
    uint64_t ub;
    uint64_t lb;
  uint64_t total;
    uint64_t getUB(){ return ub; }
    uint64_t getLB(){ return lb; }

    void my_print(indexMap indexToName, bool original_vars = true, int maxsize=20) {
        // Assume _sign == false...
        printf("c\tmin (%d)\t", _coeffs.size());
        if(maxsize < 0)
            maxsize = _coeffs.size();
        

        for (int i = 0; i < maxsize && i < _coeffs.size(); i++) {
            printf("%lu ", _coeffs[i]);
            if (sign(_lits[i]))
                printf("~");
            if(!original_vars)
                printf("x%d ", var(_lits[i]));
            else
                if(indexToName.find(var(_lits[i])) != indexToName.end())
                    printf("%s ", indexToName.at(var(_lits[i])).c_str());
                else
                    printf("x%d ", var(_lits[i]));
        }
        if(_coeffs.size() > maxsize)
            printf("... (%d lits)", _coeffs.size());
        printf("\n");
        
        double tweight = 0;
        for (int i = 0; i < _coeffs.size(); i++)
            tweight += _coeffs[i];
        printf("c max obj value: %.0f\n", tweight);
        printf("c _const: %ld\n", _const);
    }
  bool empty(){
    return _coeffs.size() == _lits.size() && (_lits.size() == 1);
}
    
    };
  std::unique_ptr<PBObjFunction>  add(PBObjFunction* pb, PBObjFunction* pbb);
} // namespace openwbo

#endif
