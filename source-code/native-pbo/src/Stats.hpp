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

#include "auxiliary.hpp"
#include "typedefs.hpp"

namespace rs {

struct Stats {
  long long NTRAILPOPS = 0, NWATCHLOOKUPS = 0, NWATCHLOOKUPSBJ = 0, NWATCHCHECKS = 0, NPROPCHECKS = 0,
            NADDEDLITERALS = 0;
  long long NCONFL = 0, NDECIDE = 0, NPROP = 0, NPROPCLAUSE = 0, NPROPCARD = 0, NPROPWATCH = 0, NPROPCOUNTING = 0,
            NRESOLVESTEPS = 0;
  long long NWATCHED = 0, NCOUNTING = 0;
  int128 EXTERNLENGTHSUM = 0, LEARNEDLENGTHSUM = 0;
  bigint EXTERNDEGREESUM = 0, LEARNEDDEGREESUM = 0;
  long long NCLAUSESEXTERN = 0, NCARDINALITIESEXTERN = 0, NGENERALSEXTERN = 0;
  long long NCLAUSESLEARNED = 0, NCARDINALITIESLEARNED = 0, NGENERALSLEARNED = 0;
  long long NGCD = 0, NCARDDETECT = 0, NCORECARDINALITIES = 0, NCORES = 0, NSOLS = 0;
  long long NWEAKENEDNONIMPLYING = 0, NWEAKENEDNONIMPLIED = 0;
  long long NRESTARTS = 0, NCLEANUP = 0;
  double STARTTIME = 0;
  long long NORIGVARS = 0, NAUXVARS = 0;
  long long NCONSFORMULA = 0, NCONSLEARNED = 0, NCONSBOUND = 0, NCONSCOREGUIDED = 0;
  long long NENCFORMULA = 0, NENCLEARNED = 0, NENCBOUND = 0,
            NENCCOREGUIDED;  // Number of times a reason constraint of this type was encountered.

  long long NLPADDEDROWS = 0, NLPDELETEDROWS = 0;
  long long NLPPIVOTSINTERNAL = 0, NLPPIVOTSROOT = 0, NLPNOPIVOT = 0, NLPRESETBASIS = 0;
  double LPSOLVETIME = 0, LPTOTALTIME = 0;
  long long NLPCALLS = 0, NLPOPTIMAL = 0, NLPINFEAS = 0, NLPFARKAS = 0;
  long long NLPCYCLING = 0, NLPNOPRIMAL = 0, NLPNOFARKAS = 0, NLPSINGULAR = 0, NLPOTHER = 0;
  long long NLPGOMORYCUTS = 0, NLPLEARNEDCUTS = 0, NLPLEARNEDFARKAS = 0, NLPDELETEDCUTS = 0;
  long long NLPENCGOMORY = 0, NLPENCFARKAS = 0,
            NLPENCLEARNEDFARKAS = 0;  // Number of times a reason constraint of this type was encountered.

  long long UNITCORES = 0, SINGLECORES = 0, REMOVEDBLOCKS = 0, FIRSTCOREBEST = 0, DECCOREBEST = 0, NOCOREBEST = 0,
            COREDEGSUM = 0, CORESLACKSUM = 0;
  long long NEXHAUSTED = 0, NTRIMMED = 0, NMINIMIZED = 0, NCLAUSALCORES = 0, NCARDCORES = 0, NPBCORES = 0;
  long long NPBREFORMS = 0, NAM1S = 0, AM1LB = 0, NDISCARDEDCORES = 0;

  double SOLVETIME = 0, SOLVETIMECG = 0, SOLVETIMEEXHAUST = 0, SOLVETIMETRIM = 0, SOLVETIMEMIN = 0, CATIME = 0,
         PROPTIME = 0;
  double RUNSTARTTIME = 0, PROPTIMEAM1 = 0;

  inline double getTime() const { return aux::cpuTime() - STARTTIME; }
  inline double getRunTime() const { return aux::cpuTime() - RUNSTARTTIME; }
  inline double getSolveTime() const { return SOLVETIME + SOLVETIMECG; }

  inline long long getDetTime() const {
    return 1 + NADDEDLITERALS + NWATCHLOOKUPS + NWATCHLOOKUPSBJ + NWATCHCHECKS + NPROPCHECKS + NPROP + NTRAILPOPS +
           NDECIDE + NLPPIVOTSROOT + NLPPIVOTSINTERNAL;
  }

  void print() const {
    printf("c cpu time %g s\n", getTime());
    printf("c deterministic time %lld %.2e\n", getDetTime(), (double)getDetTime());
    printf("c optimization time %g s\n", getRunTime() - getSolveTime());
    printf("c total solve time %g s\n", getSolveTime());
    printf("c core-guided solve time %g s\n", SOLVETIMECG);
    printf("c propagation time %g s\n", PROPTIME);
    printf("c conflict analysis time %g s\n", CATIME);
    printf("c propagations %lld\n", NPROP);
    printf("c resolve steps %lld\n", NRESOLVESTEPS);
    printf("c decisions %lld\n", NDECIDE);
    printf("c conflicts %lld\n", NCONFL);
    printf("c restarts %lld\n", NRESTARTS);
    printf("c inprocessing phases %lld\n", NCLEANUP);
    printf("c input clauses %lld\n", NCLAUSESEXTERN);
    printf("c input cardinalities %lld\n", NCARDINALITIESEXTERN);
    printf("c input general constraints %lld\n", NGENERALSEXTERN);
    long long nonLearneds = NCLAUSESEXTERN + NCARDINALITIESEXTERN + NGENERALSEXTERN;
    printf("c input average constraint length %.2f\n", nonLearneds == 0 ? 0 : (double)EXTERNLENGTHSUM / nonLearneds);
    printf("c input average constraint degree %.2f\n", nonLearneds == 0 ? 0 : (double)EXTERNDEGREESUM / nonLearneds);
    printf("c learned clauses %lld\n", NCLAUSESLEARNED);
    printf("c learned cardinalities %lld\n", NCARDINALITIESLEARNED);
    printf("c learned general constraints %lld\n", NGENERALSLEARNED);
    long long learneds = NCLAUSESLEARNED + NCARDINALITIESLEARNED + NGENERALSLEARNED;
    printf("c learned average constraint length %.2f\n", learneds == 0 ? 0 : (double)LEARNEDLENGTHSUM / learneds);
    printf("c learned average constraint degree %.2f\n", learneds == 0 ? 0 : (double)LEARNEDDEGREESUM / learneds);
    printf("c watched constraints %lld\n", NWATCHED);
    printf("c counting constraints %lld\n", NCOUNTING);
    printf("c gcd simplifications %lld\n", NGCD);
    printf("c detected cardinalities %lld\n", NCARDDETECT);
    printf("c weakened non-implied lits %lld\n", NWEAKENEDNONIMPLIED);
    printf("c weakened non-implying lits %lld\n", NWEAKENEDNONIMPLYING);
    printf("c original variables %lld\n", NORIGVARS);
    printf("c clausal propagations %lld\n", NPROPCLAUSE);
    printf("c cardinality propagations %lld\n", NPROPCARD);
    printf("c watched propagations %lld\n", NPROPWATCH);
    printf("c counting propagations %lld\n", NPROPCOUNTING);
    printf("c watch lookups %lld\n", NWATCHLOOKUPS);
    printf("c watch backjump lookups %lld\n", NWATCHLOOKUPSBJ);
    printf("c watch checks %lld\n", NWATCHCHECKS);
    printf("c propagation checks %lld\n", NPROPCHECKS);
    printf("c constraint additions %lld\n", NADDEDLITERALS);
    printf("c trail pops %lld\n", NTRAILPOPS);
    printf("c formula constraints %lld\n", NCONSFORMULA);
    printf("c learned constraints %lld\n", NCONSLEARNED);
    printf("c bound constraints %lld\n", NCONSBOUND);
    printf("c core-guided constraints %lld\n", NCONSCOREGUIDED);
    printf("c encountered formula constraints %lld\n", NENCFORMULA);
    printf("c encountered learned constraints %lld\n", NENCLEARNED);
    printf("c encountered bound constraints %lld\n", NENCBOUND);
    printf("c encountered core-guided constraints %lld\n", NENCCOREGUIDED);
    printf("c LP total time %g s\n", LPTOTALTIME);
    printf("c LP solve time %g s\n", LPSOLVETIME);
    printf("c LP constraints added %lld\n", NLPADDEDROWS);
    printf("c LP constraints removed %lld\n", NLPDELETEDROWS);
    printf("c LP pivots internal %lld\n", NLPPIVOTSINTERNAL);
    printf("c LP pivots root %lld\n", NLPPIVOTSROOT);
    printf("c LP calls %lld\n", NLPCALLS);
    printf("c LP optimalities %lld\n", NLPOPTIMAL);
    printf("c LP no pivot count %lld\n", NLPNOPIVOT);
    printf("c LP infeasibilities %lld\n", NLPINFEAS);
    printf("c LP valid Farkas constraints %lld\n", NLPFARKAS);
    printf("c LP learned Farkas constraints %lld\n", NLPLEARNEDFARKAS);
    printf("c LP basis resets %lld\n", NLPRESETBASIS);
    printf("c LP cycling count %lld\n", NLPCYCLING);
    printf("c LP singular count %lld\n", NLPSINGULAR);
    printf("c LP no primal count %lld\n", NLPNOPRIMAL);
    printf("c LP no farkas count %lld\n", NLPNOFARKAS);
    printf("c LP other issue count %lld\n", NLPOTHER);
    printf("c LP Gomory cuts %lld\n", NLPGOMORYCUTS);
    printf("c LP learned cuts %lld\n", NLPLEARNEDCUTS);
    printf("c LP deleted cuts %lld\n", NLPDELETEDCUTS);
    printf("c LP encountered Gomory constraints %lld\n", NLPENCGOMORY);
    printf("c LP encountered Farkas constraints %lld\n", NLPENCFARKAS);
    printf("c LP encountered learned Farkas constraints %lld\n", NLPENCLEARNEDFARKAS);
    printf("c CG auxiliary variables introduced %lld\n", NAUXVARS);
    printf("c CG solutions found %lld\n", NSOLS);
    printf("c CG cores constructed %lld\n", NCORES);
    printf("c CG core cardinality constraints returned %lld\n", NCORECARDINALITIES);
    printf("c CG unit cores %lld\n", UNITCORES);
    printf("c CG single cores %lld\n", SINGLECORES);
    printf("c CG blocks removed during cardinality reduction %lld\n", REMOVEDBLOCKS);
    printf("c CG first core best %lld\n", FIRSTCOREBEST);
    printf("c CG decision core best %lld\n", DECCOREBEST);
    printf("c CG core reduction tie %lld\n", NOCOREBEST);
    printf("c CG core degree average %.2f\n",
           (NCORES - UNITCORES) == 0 ? 0 : COREDEGSUM / (double)(NCORES - UNITCORES));
    printf("c CG core slack average %.2f\n",
           (NCORES - UNITCORES) == 0 ? 0 : CORESLACKSUM / (double)(NCORES - UNITCORES));
    printf("c CG exhausted %lld\n", NEXHAUSTED);
    printf("c CG exhaust time %g s\n", SOLVETIMEEXHAUST);
    printf("c CG trimmed %lld\n", NTRIMMED);
    printf("c CG trim time %g s\n", SOLVETIMETRIM);
    printf("c CG minimized %lld\n", NMINIMIZED);
    printf("c CG minimize time %g s\n", SOLVETIMEMIN);
    printf("c CG am1 prop time %g s\n", PROPTIMEAM1);
    printf("c CG # am1s reformulated %lld\n", NAM1S);
    printf("c CG lower bound derived from am1s %lld\n", AM1LB);
    printf("c CG # clausal cores %lld\n", NCLAUSALCORES);
    printf("c CG # cardinality cores %lld\n", NCARDCORES);
    printf("c CG # pb cores %lld\n", NPBCORES);
    printf("c CG # pb reforms %lld\n", NPBREFORMS);
    printf("c CG # cores discarded because of pb reform %lld\n", NDISCARDEDCORES);
  }
};

}  // namespace rs
