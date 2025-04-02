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

#include "Logger.hpp"
#include <ranges>
#if defined IOSTREAMS_WITH_ZLIB
#include <boost/iostreams/filter/gzip.hpp>
#endif

#include "ProofBuffer.hpp"

namespace rs {

Logger::Logger(std::ostream& proof_out) : proof_out(proof_out) { proof_out << "pseudo-Boolean proof version 2.0\n"; }

void Logger::init(uint num_constraints) {
  proof_out << "f " << num_constraints << "\n";
  last_proofID = num_constraints;
  // Dummy constraint for easier proof logging.
  proof_out << "rup >= 0 ;\n";
  ID_Trivial = ++last_proofID;
  default_proof_line = std::to_string(ID_Trivial) + " ";
}

ID Logger::solution(const std::vector<Lit>& sol) {
  proof_out << "soli";
  for (Lit x : sol | std::views::drop(1)) proof_out << proofLit(x);
  proof_out << "\n";
  return ++last_proofID;
}

ID Logger::solx(const std::vector<Lit>& sol) {
  proof_out << "solx";
  for (Lit x : sol | std::views::drop(1)) proof_out << proofLit(x);
  proof_out << "\n";
  return ++last_proofID;
}

void Logger::remove(ID id) {
  return;  // TODO: does not play well with the LP solver; disabling
  proof_out << "del id " << id << "\n";
}

void Logger::remove(const std::vector<ID>& ids) {
  return;  // TODO: does not play well with the LP solver; disabling
  if (ids.empty()) return;
  proof_out << "del id";
  for (ID id : ids) proof_out << " " << id;
  proof_out << "\n";
}

void Logger::sat() {
  // TODO: log bounds
  if (optimization) return bailout();
  proof_out << "output NONE\n";
  proof_out << "conclusion SAT\n";
  proof_out << "end pseudo-Boolean proof\n";
  flush();
}

void Logger::unsat() {
  proof_out << "output NONE\n";
  if (optimization) {
    proof_out << "conclusion BOUNDS INF INF\n";
  } else {
    proof_out << "conclusion UNSAT : " << last_proofID << "\n";
  }
  proof_out << "end pseudo-Boolean proof\n";
  flush();
}

template <typename LARGE>
void Logger::optimum(LARGE opt) {
  assert(optimization);
  proof_out << "output NONE\n";
  proof_out << "conclusion BOUNDS " << opt << " " << opt << "\n";
  proof_out << "end pseudo-Boolean proof\n";
  flush();
}

template void Logger::optimum<long long>(long long);
template void Logger::optimum<int128>(int128);
template void Logger::optimum<int256>(int256);
template void Logger::optimum<bigint>(bigint);

void Logger::bailout() {
  proof_out << "output NONE\n";
  proof_out << "conclusion NONE\n";
  proof_out << "end pseudo-Boolean proof\n";
  flush();
}

};  // namespace rs
