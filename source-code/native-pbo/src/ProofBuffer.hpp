/***********************************************************************
Copyright (c) 2014-2020, Jan Elffers
Copyright (c) 2019-2021, Jo Devriendt
Copyright (c) 2020-2021, Stephan Gocht
Copyright (c) 2014-2021, Jakob Nordström

Parts of the code were copied or adapted from MiniSat.

MiniSat -- Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson
           Copyright (c) 2007-2010  Niklas Sorensson

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
***********************************************************************/

#pragma once

#include <sstream>
#include "typedefs.hpp"

namespace rs {

struct Logger;

// Syntactic sugar
template <typename T>
struct proofMult {
  const T& m;
  proofMult(const T& mult) : m(mult) {}
};
template <typename T>
std::ostream& operator<<(std::ostream& os, const proofMult<T>& mult) {
  if (mult.m != 1) os << mult.m << " * ";
  return os;
}

// Syntactic sugar
template <typename T>
struct proofDiv {
  const T& d;
  proofDiv(const T& div) : d(div) {}
};
template <typename T>
std::ostream& operator<<(std::ostream& os, const proofDiv<T>& div) {
  return os << div.d << " d ";
}

// Syntactic sugar
struct proofLit {
  Lit l;
  proofLit(Lit l) : l(l) {}
};
inline std::ostream& operator<<(std::ostream& o, proofLit l) { return o << (l.l < 0 ? " ~x" : " x") << toVar(l.l); }

enum class ProofRule { POLISH, REDUNDANT, RUP };

// Tree of operations used to derive a constraint,
// represented in reverse Polish notation
struct ProofBuffer {
  std::stringstream buffer;
  enum ProofRule rule = ProofRule::POLISH;
  size_t pad = 0;
  operator std::ostream&() { return buffer; }

  bool isTrivial() const;

  void reset(ID proofId);
  void copyFrom(const std::string& other);

  template <typename SMALL>
  void addVariable(Var v, const SMALL& m) {
    if (m > 0)
      buffer << "x" << v << " " << proofMult(m) << "+ ";
    else
      buffer << "~x" << v << " " << proofMult(-m) << "+ ";
  }
  template <typename SMALL>
  void addLiteral(Lit l, const SMALL& c) {
    buffer << proofLit(l) << " " << proofMult(c) << "+ ";
  }
  template <typename SMALL>
  void addClause(ID id, const SMALL& c) {
    buffer << id << " " << proofMult(c) << "+ ";
  }
  template <typename SMALL>
  void addConstraint(const SMALL& thismult, ID id, const SMALL& c) {
    buffer << proofMult(thismult) << " " << id << " " << proofMult(c) << "+ ";
  }
  template <typename SMALL>
  void addConstraint(const SMALL& thismult, const ProofBuffer& other, const SMALL& othermult) {
    buffer << proofMult(thismult) << other.buffer.view() << proofMult(othermult) << "+ ";
  }
  void saturate() { buffer << "s "; }
  void addToWitness(Var v, int sign);
  void addSubProof(ID id);

  ID logPolish(Logger& logger);
  ID logRedundant(Logger& logger, CeSuper ce);
  ID logRup(Logger& logger, CeSuper ce);
  ID logImplied(Logger& logger, CeSuper ce);
  ID logAssumption(Logger& logger, CeSuper ce);
  ID logExternal(const std::string& rule, Logger& logger, CeSuper ce);
  ID log(Logger& logger, CeSuper ce);
};

}  // namespace rs
