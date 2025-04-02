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

#include "run.hpp"
#include "ConstrExp.hpp"
#include "Logger.hpp"
#include "Solver.hpp"
#include "auxiliary.hpp"
#include "bioptsat.hpp"
#include "optimize.hpp"
#include "paretopk.hpp"
#include "pminimal.hpp"

namespace rs {

std::default_random_engine aux::rng;

Solver run::solver;
run::MultiOptBase* run::mosolver = nullptr;

void run::decide() {
  while (true) {
    SolveState reply = aux::timeCall<SolveState>([&] { return solver.solve().state; }, stats.SOLVETIME);
    assert(reply != SolveState::INCONSISTENT);
    if (reply == SolveState::SAT) {
      quit::exit_SAT(solver);
      return;
    } else if (reply == SolveState::UNSAT) {
      quit::exit_UNSAT(solver);
      return;
    }
    if (options.time_limit.get() != -1.0 && stats.getTime() > options.time_limit.get()) {
      quit::exit_INDETERMINATE(solver);
      return;
    }
  }
}

template <typename CF, typename DG>
run::MultiOptBase* initMultiOptSolver(std::vector<CeArb> objectives) {
  if (options.moAlg.is("p-minimal")) return new run::PMinimal<CF, DG>(objectives);
  if (options.moAlg.is("paretop-k")) return new run::PareTopK<CF, DG>(objectives);
  if (options.moAlg.is("bioptsat")) {
    if (objectives.size() != 2) {
      std::cerr << "c BiOptSat Algorithm requires exactly 2 objectives" << std::endl;
      exit(1);
    }
    return new run::BiOptSat<CF, DG>(objectives);
  }
  return nullptr;
}

void run::run(std::vector<CeArb> objectives) {
  stats.RUNSTARTTIME = aux::cpuTime();
  if (options.verbosity.get() > 0)
    std::cout << "c #variables " << solver.getNbOrigVars() << " #constraints " << solver.getNbConstraints()
              << std::endl;
  aux::rng.seed(options.seed.get());
  try {
    if (objectives.empty()) {
      decide();
    } else if (objectives.size() == 1) {
      CeArb objective = objectives[0];
      objective->stopLogging();
      objective->removeUnitsAndZeroes(solver, false);

      BigVal maxVal = objective->getCutoffVal();
      if (maxVal <= limit32) {
        if (options.verbosity.get() > 0) std::cout << "c SMALL=int, LARGE=long long" << std::endl;
        Optimization<int, long long> optim(objective);
        optim.optimize();
      } else if (maxVal <= limit64) {
        if (options.verbosity.get() > 0) std::cout << "c SMALL=long long, LARGE=int128" << std::endl;
        Optimization<long long, int128> optim(objective);
        optim.optimize();
      } else if (maxVal <= BigVal(limit96)) {
        if (options.verbosity.get() > 0) std::cout << "c SMALL=int128, LARGE=int128" << std::endl;
        Optimization<int128, int128> optim(objective);
        optim.optimize();
      } else if (maxVal <= BigVal(limit128)) {
        if (options.verbosity.get() > 0) std::cout << "c SMALL=int128, LARGE=int256" << std::endl;
        Optimization<int128, int256> optim(objective);
        optim.optimize();
      } else {
        if (options.verbosity.get() > 0) std::cout << "c SMALL=bigint, LARGE=bigint" << std::endl;
        Optimization<bigint, bigint> optim(objective);
        optim.optimize();
      }
    } else {
      std::vector<BigVal> maxVals{};
      maxVals.reserve(objectives.size());
      for (auto obj : objectives) {
        obj->stopLogging();
        obj->removeUnitsAndZeroes(solver, false);

        maxVals.push_back(obj->getCutoffVal());
      }
      BigVal maxVal = aux::max(maxVals);
      if (maxVal <= limit32) {
        if (options.verbosity.get() > 0) std::cout << "c SMALL=int, LARGE=long long" << std::endl;
        mosolver = initMultiOptSolver<int, long long>(objectives);
      } else if (maxVal <= limit64) {
        if (options.verbosity.get() > 0) std::cout << "c SMALL=long long, LARGE=int128" << std::endl;
        mosolver = initMultiOptSolver<long long, int128>(objectives);
      } else if (maxVal <= BigVal(limit96)) {
        if (options.verbosity.get() > 0) std::cout << "c SMALL=int128, LARGE=int128" << std::endl;
        mosolver = initMultiOptSolver<int128, int128>(objectives);
      } else if (maxVal <= BigVal(limit128)) {
        if (options.verbosity.get() > 0) std::cout << "c SMALL=int128, LARGE=int256" << std::endl;
        mosolver = initMultiOptSolver<int128, int256>(objectives);
      } else {
        if (options.verbosity.get() > 0) std::cout << "c SMALL=bigint, LARGE=bigint" << std::endl;
        mosolver = initMultiOptSolver<bigint, bigint>(objectives);
      }
      mosolver->solve();
      mosolver->printParetoFront((bool)options.printSol);
      mosolver->printStats();
      if (solver.logger) solver.logger->unsat();
      delete mosolver;
    }
  } catch (const AsynchronousInterrupt& ai) {
    std::cout << "c " << ai.what() << std::endl;
    if (mosolver) {
      mosolver->printParetoFront((bool)options.printSol);
      mosolver->printStats();
      if (solver.logger) solver.logger->bailout();
      delete mosolver;
      std::cout << "s UNKNOWN" << std::endl;
    } else
      quit::exit_INDETERMINATE(solver);
  }
}

}  // namespace rs
