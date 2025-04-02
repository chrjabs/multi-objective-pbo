João O'Neill Cortes


Table of Contents
─────────────────

1. Version 0.1
.. 1. Compile
..... 1. Without CPLEX
..... 2. With CPLEX
.. 2. Usage
..... 1. Input format
..... 2. Output format
.. 3. Datasets
..... 1. TACAS2023


Authors: João Cortes, Vasco Manquinho, Ines Lynce

To contact the authors please send an email to:
joaooneillcortes@outlook.pt


1 Version 0.1
═════════════

  This tool solves the multi-objective optimization problem (MOCO).

  Given some instance (check the example below) you will get a set of
  /Pareto-optimal/ solutions, if you wait long enough. If the solver is
  interrupted, it will report back an approximation to the complete set
  of Pareto-optimal solutions.


1.1 Compile
───────────

  You can choose between installing CPLEX suite. You will get smaller
  encodings if you do, because CPLEX is used to compute a priori tighter
  bounds on the objective functions's ranges.


1.1.1 Without CPLEX
╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌

  Hopefully, you will be able to build everything with a single make
  call,
  ┌────
  │ make cplex_bounds_dependency=False
  └────


1.1.2 With CPLEX
╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌

  You must install IBM's CPLEX suite first. Check
  <https://www.ibm.com/products/ilog-cplex-optimization-studio/cplex-optimizer>
  for details. There is a free licence for academic purposes (as of Mon
  5 Dec 12:00:19 WET 2022).

  You will need to edit the value of the variables CPLEXDIR and
  CONCERTDIR, by editing cplex.make.

  Hopefully, you will be able to build everything with a single make
  call,

  ┌────
  │ make cplex_bounds_dependency=True
  └────


1.2 Usage
─────────

  You can run any of the following algorithms:
  1. P-Minimal (based on `sat4j` framework), label `p-minimal-sat4j`;
     <https://link.springer.com/chapter/10.1007/978-3-319-66158-2_38>
  2. P-Minimal (based on `openwbo` solver), label `p-minimal-openwbo`;
  3. ParetoMCS (sat4j), label `pareto-mcs`;
     <https://www.researchgate.net/publication/326204042_Multi-Objective_Optimization_Through_Pareto_Minimal_Correction_Subsets>
  4. Hitting-Sets (openwbo), label `hitting-sets`;
     <https://arxiv.org/abs/2204.10856>
  5. Core-Guided (openwbo), label `core-guided`;
     <https://arxiv.org/abs/2204.10856>

  For example:
  ┌────
  │ ./RUNME.sh pareto-mcs inputfile.pbmo
  └────


1.2.1 Input format
╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌

  The file should start with a list of Pseudo-Boolean functions to
  minimize, followed by a set of pseudo-Boolean restrictions. For
  instance,

        min: +1 ~x1 -5 x2 min: +2 x1 +9 x3 x1 + x2 < 1 x1 + x3 >=
        1

  The tilde operator '~' negates a boolean variable. The weights can be
  either positive or negative integers.


1.2.2 Output format
╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌

  Typical problems will take a long time to compute the full
  Pareto-front. In any case, the solver will report back a partial
  solution if interrupted by SIGINT.

  The solver follows the standard output of MaxSAT solvers:

  ⁃ Comments ("c " lines);
  ⁃ Solution Status ("s " line):
    • `s OPTIMUM` : the full Pareto-front was found;
    • `s UNSATISFIABLE` : the hard clauses are unsatisfiable;
    • `s SATISFIABLE` : some solutions were found but optimality was not
      proven;
  ⁃ Solution Cost Line ("o " lines):
    • This represents the cost of some solution found by the solver;
  ⁃ Solution Values (Truth Assignment) ("v " lines):
    • This represents the truth assignment (true/false) assigned to each
    variable;
  ⁃ Number of Solutions ("n " line): the cardinality of the current
    partial solution;


1.3 Datasets
────────────

1.3.1 TACAS2023
╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌

  You can find the full TACAS2023 dataset at [TACAS2023].


[TACAS2023] <http://sat.inesc-id.pt/~jcortes/artifacts/TACAS2023>
