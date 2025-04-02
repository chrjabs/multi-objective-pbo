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

#include "parsing.hpp"
#include <boost/iostreams/filtering_stream.hpp>
#include <sstream>
#include "Logger.hpp"
#include "Solver.hpp"
#include "quit.hpp"

namespace rs {

bigint parsing::read_number(const std::string& s) {
  bigint answer = 0;
  bool negate = false;
  for (char c : s) {
    if ('0' <= c && c <= '9') {
      answer *= 10;
      answer += c - '0';
    }
    negate = (negate || (c == '-'));
  }
  return negate ? -answer : answer;
}

bool parsing::opb_read(std::istream& in, Solver& solver, std::vector<CeArb>& objectives, std::string line) {
  assert(objectives.empty());
  CeArb input = solver.cePools.takeArb();
  [[maybe_unused]] bool first_constraint = true;
  for (; line.empty() ? getline(in, line) : in;) {
    if (line.empty() || line[0] == '*') {
      line.clear();
      continue;
    }
    bool opt_line = line.substr(0, 4) == "min:";
    std::string line0;
    if (opt_line)
      line = line.substr(4), assert(first_constraint);
    else {
      std::string symbol;
      if (line.find(">=") != std::string::npos)
        symbol = ">=";
      else
        symbol = "=";
      assert(line.find(symbol) != std::string::npos);
      line0 = line;
      line = line.substr(0, line.find(symbol));
    }
    std::istringstream is(line);
    input->reset();
    std::vector<std::string> tokens;
    std::string tmp;
    while (is >> tmp) tokens.push_back(tmp);
    if (tokens.back() == ";")
      tokens.pop_back();
    else if (tokens.back().back() == ';')
      tokens.back().pop_back();
    if (tokens.size() % 2 != 0) quit::exit_UNSUPPORTED("non-linear constraints");
    for (size_t i = 0; i < tokens.size(); i += 2)
      if (find(tokens[i].begin(), tokens[i].end(), 'x') != tokens[i].end())
        quit::exit_UNSUPPORTED("non-linear constraints");
    for (size_t i = 0; i < tokens.size(); i += 2) {
      std::string scoef = tokens[i];
      std::string var = tokens[i + 1];
      BigCoef coef = read_number(scoef);
      bool negated = false;
      std::string origvar = var;
      if (!var.empty() && var[0] == '~') {
        negated = true;
        var = var.substr(1);
      }
      if (var.empty() || var[0] != 'x') quit::exit_ERROR({"Invalid literal token: ", origvar});
      var = var.substr(1);
      Lit l = stoi(var);
      if (l < 1) quit::exit_ERROR({"Variable token less than 1: ", origvar});
      if (negated) l = -l;
      solver.setNbVars(std::abs(l), true);
      input->addLhs(coef, l);
    }
    if (opt_line) {
      CeArb obj = solver.cePools.takeArb();
      input->copyTo(obj);
      objectives.push_back(obj);
    } else {
      first_constraint = false;
      input->addRhs(read_number(line0.substr(line0.find("=") + 1)));
      if (solver.addConstraint(input, Origin::FORMULA).second == ID_Unsat) {
        quit::exit_UNSAT(solver);
        return true;
      }
      if (line0.find(" = ") != std::string::npos) {  // Handle equality case with second constraint
        input->invert();
        if (solver.addConstraint(input, Origin::FORMULA).second == ID_Unsat) {
          quit::exit_UNSAT(solver);
          return true;
        }
      }
    }
    line.clear();
  }
  return false;
}

bool parsing::wcnf_read(std::istream& in, BigCoef top, Solver& solver, CeArb objective) {
  assert(objective->isReset());
  CeArb input = solver.cePools.takeArb();
  for (std::string line; getline(in, line);) {
    if (line.empty() || line[0] == 'c')
      continue;
    else {
      std::istringstream is(line);
      std::string sweight;
      is >> sweight;
      BigCoef weight = read_number(sweight);
      if (weight == 0) continue;
      input->reset();
      input->addRhs(1);
      Lit l;
      while (is >> l, l) input->addLhs(1, l);
      if (weight < top) {  // soft clause
        std::stringstream s;
        s << "Negative clause weight: " << weight;
        if (weight < 0) quit::exit_ERROR({s.str()});
        if (input->vars.size() == 1) {
          // avoid adding constraint for unit soft clauses
          Lit l = input->getLit(input->vars[0]);
          objective->addLhs(weight, -l);
          continue;
        }
        solver.setNbVars(solver.getNbVars() + 1, true);  // increases n to n+1
        objective->addLhs(weight, solver.getNbVars());
        input->addLhs(1, solver.getNbVars());
      }  // else hard clause
      if (solver.addConstraint(input, Origin::FORMULA).second == ID_Unsat) {
        quit::exit_UNSAT(solver);
        return true;
      }
    }
  }
  return false;
}

bool parsing::cnf_read(std::istream& in, Solver& solver) {
  Ce32 input = solver.cePools.take32();
  for (std::string line; getline(in, line);) {
    if (line.empty() || line[0] == 'c')
      continue;
    else {
      std::istringstream is(line);
      input->reset();
      input->addRhs(1);
      Lit l;
      while (is >> l, l) {
        solver.setNbVars(std::abs(l), true);
        input->addLhs(1, l);
      }
      if (solver.addConstraint(input, Origin::FORMULA).second == ID_Unsat) {
        quit::exit_UNSAT(solver);
        return true;
      }
    }
  }
  return false;
}

bool parsing::file_read(boost::iostreams::filtering_istream& in, Solver& solver, std::vector<CeArb>& objectives) {
  for (std::string line; getline(in, line);) {
    if (line.empty() || line[0] == 'c') continue;
    if (line[0] == 'p') {
      std::istringstream is(line);
      is >> line;  // skip 'p'
      std::string type;
      long long num_variables, num_constraints;
      is >> type >> num_variables >> num_constraints;
      solver.setNbVars(num_variables);
      solver.reserveConstraints(num_constraints);
      if (solver.logger) {
        solver.logger->init(num_constraints);
      }
      if (type == "cnf") {
        return cnf_read(in, solver);
      } else if (type == "wcnf") {
        std::string stop;
        is >> stop;  // top
        BigCoef top = read_number(stop);
        objectives.push_back(solver.cePools.takeArb());
        return wcnf_read(in, top, solver, objectives.back());
      }
    } else if (line[0] == '*' || line[0] == 'm' || line[0] == '+' || line[0] == '-') {
      // parse hint line
      if (line.substr(0, 11) == "* #variable") {
        std::istringstream is(line);
        std::vector<std::string> tokens;
        std::string tmp;
        while (is >> tmp) tokens.push_back(tmp);
        if (tokens.size() > 2) {
          size_t offset = 2;
          // reserve variables
          if (tokens[2] == "=") {
            offset++;
          }
          solver.setNbVars(stoi(tokens[offset]), true);
          offset++;
          long long num_constraints = 0;
          for (; offset < tokens.size(); offset += 2) {
            if (tokens[offset] == "#constraint=" && tokens.size() > offset + 1) {
              num_constraints += stoul(tokens[offset + 1]);
            } else if (tokens[offset] == "#constraint" && tokens.size() > offset + 2 && tokens[offset + 1] == "=") {
              num_constraints += stoul(tokens[offset + 2]);
            } else if (tokens[offset] == "#equal=" && tokens.size() > offset + 1) {
              num_constraints += stoul(tokens[offset + 1]);
            } else if (tokens[offset] == "#equal=" && tokens.size() > offset + 2 && tokens[offset + 1] == "=") {
              num_constraints += stoul(tokens[offset + 2]);
            } else if (options.bitsInput.get() != 0 && tokens[offset] == "intsize=" && tokens.size() > offset + 1 &&
                       stoi(tokens[offset + 1]) > options.bitsInput.get()) {
              quit::exit_UNSUPPORTED("coefficient size");
            } else if (options.bitsInput.get() != 0 && tokens[offset] == "intsize" && tokens.size() > offset + 2 &&
                       tokens[offset + 2] == "=" && stoi(tokens[offset + 1]) > options.bitsInput.get()) {
              quit::exit_UNSUPPORTED("coefficient size");
            }
          }
          solver.reserveConstraints(num_constraints);
          if (num_constraints > 0 && solver.logger) {
            solver.logger->init(num_constraints);
          }
        }
      }
      if (solver.logger && solver.logger->last_proofID == 0) {
        quit::exit_ERROR({"Need a hint for the number of constraints when proof logging"});
      }
      return opb_read(in, solver, objectives, line);
    } else {
      quit::exit_ERROR({"No supported format [(m)opb, cnf, wcnf] detected."});
      return true;
    }
  }
  return true;
}

}  // namespace rs
