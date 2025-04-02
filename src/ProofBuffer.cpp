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

#include "ProofBuffer.hpp"
#include "ConstrExp.hpp"
#include "Logger.hpp"

namespace rs {

bool ProofBuffer::isTrivial() const {
  auto str = buffer.view();
  assert(str.back() == ' ');
  long long spacecount = 0;
  for (char const& c : str) {
    spacecount += (c == ' ');
    if (spacecount > 1) return false;
  }
  return true;
}

void ProofBuffer::reset(ID proofID) {
  assert(proofID != ID_Undef);
  assert(proofID != ID_Unsat);
  buffer.clear();
  buffer.str(std::string());
  buffer << proofID << " ";
  rule = ProofRule::POLISH;
  pad = 0;
}

void ProofBuffer::copyFrom(const std::string& other) {
  buffer.clear();
  buffer.str(std::string());
  buffer << other;
  rule = ProofRule::POLISH;
  pad = 0;
}

// TODO: support more than one literal
void ProofBuffer::addToWitness(Var v, int sign) {
  rule = ProofRule::REDUNDANT;
  buffer.clear();
  buffer.str(std::string());
  buffer << " x" << v << " -> " << sign << " ;";
}

// TODO: support more than one objective
void ProofBuffer::addSubProof(ID id) {
  buffer << " begin\nproofgoal #1\npol -1 " << id << " +\nend -1\nend";
  pad = 3;
}

ID ProofBuffer::logPolish(Logger& logger) {
  ID id;
  if (isTrivial()) {  // line is just one id, don't print it
    buffer >> id;
    buffer.seekg(0);
    assert(id != ID_Invalid);
    assert(id != logger.ID_Trivial);
  } else {  // non-trivial line
    id = ++logger.last_proofID;
    logger.proof_out << "p " << buffer.view() << std::endl;
    reset(id);
  }
  return id;
}

ID ProofBuffer::logRedundant(Logger& logger, CeSuper ce) {
  logger.proof_out << "red ";
  ce->toStreamAsOPB(logger.proof_out);
  logger.proof_out << " ; " << buffer.view() << std::endl;
  logger.last_proofID += pad;
  ID id = ++logger.last_proofID;
  reset(id);  // ensure consistent proofBuffer
  return id;
}

ID ProofBuffer::logExternal(const std::string& rule, Logger& logger, CeSuper ce) {
  logger.proof_out << rule;
  ce->toStreamAsOPB(logger.proof_out);
  logger.proof_out << std::endl;
  ID id = ++logger.last_proofID;
  reset(id);  // ensure consistent proofBuffer
  return id;
}

ID ProofBuffer::logRup(Logger& logger, CeSuper ce) {
  logger.proof_out << "rup ";
  ce->toStreamAsOPB(logger.proof_out);
  logger.proof_out << " ; " << buffer.view() << std::endl;
  logger.last_proofID += pad;
  ID id = ++logger.last_proofID;
  reset(id);  // ensure consistent proofBuffer
  return id;
}

ID ProofBuffer::logImplied(Logger& logger, CeSuper ce) { return logExternal("ia ", logger, ce); }

ID ProofBuffer::logAssumption(Logger& logger, CeSuper ce) { return logExternal("a ", logger, ce); }

ID ProofBuffer::log(Logger& logger, CeSuper ce) {
  if (rule == ProofRule::POLISH)
    return logPolish(logger);
  else if (rule == ProofRule::RUP)
    return logRup(logger, ce);
  else
    return logRedundant(logger, ce);
}

}  // namespace rs
