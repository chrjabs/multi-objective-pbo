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

#include <boost/algorithm/string/predicate.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <csignal>
#include <fstream>
#include <memory>
#include "Logger.hpp"
#include "auxiliary.hpp"
#include "globals.hpp"
#include "io.hpp"
#include "parsing.hpp"
#include "run.hpp"

namespace rs {

bool asynch_interrupt;
Options options;
Stats stats;

}  // namespace rs

static void SIGINT_interrupt([[maybe_unused]] int signum) { rs::asynch_interrupt = true; }

static void SIGINT_exit([[maybe_unused]] int signum) {
  printf("\n*** INTERRUPTED ***\n");
  exit(1);
}

int main(int argc, char** argv) {
  rs::stats.STARTTIME = rs::aux::cpuTime();
  rs::asynch_interrupt = false;

  signal(SIGINT, SIGINT_exit);
  signal(SIGTERM, SIGINT_exit);
  signal(SIGXCPU, SIGINT_exit);

  rs::options.parseCommandLine(argc, argv);

  if (rs::options.verbosity.get() > 0) {
    std::cout << "c RoundingSat 2\n";
    std::cout << "c branch " << EXPANDED(GIT_BRANCH) << "\n";
    std::cout << "c commit " << EXPANDED(GIT_COMMIT_HASH) << std::endl;
  }

  std::shared_ptr<rs::Logger> logger;
  const std::string& proof_log_name = rs::options.proofLog.get();
  rs::io::ostream proof_out;

  if (!proof_log_name.empty()) {
    std::ostream& out = proof_out.open(proof_log_name);
    logger = std::make_shared<rs::Logger>(out);
  }

  rs::run::solver.init(logger);
  std::vector<rs::CeArb> objectives{};
  bool infeasible_or_error = false;

  if (!rs::options.formulaName.empty()) {
    std::ifstream fin(rs::options.formulaName, std::ifstream::in);
    if (!fin.is_open()) {
      std::cerr << "c Failed to open input file" << rs::options.formulaName << std::endl;
      exit(1);
    }
    boost::iostreams::filtering_istream in;

    if (boost::algorithm::ends_with(rs::options.formulaName, ".gz")) in.push(boost::iostreams::gzip_decompressor());
    if (boost::algorithm::ends_with(rs::options.formulaName, ".bz2")) in.push(boost::iostreams::bzip2_decompressor());
    in.push(fin);
    infeasible_or_error = rs::parsing::file_read(in, rs::run::solver, objectives);
  }
  //  else {
  //    if (rs::options.verbosity.get() > 0) std::cout << "c No filename given, reading from standard input" <<
  //    std::endl; infeasible_or_error = rs::parsing::file_read(std::cin, rs::run::solver, objective);
  //  }

  if (infeasible_or_error) return 0;

  signal(SIGINT, SIGINT_interrupt);
  signal(SIGTERM, SIGINT_interrupt);
  signal(SIGXCPU, SIGINT_interrupt);

  rs::run::solver.initLP(objectives);

  rs::run::run(objectives);
}
