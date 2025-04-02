using MoaMopb
MoaMopb.main()
import HiGHS
import MultiObjectiveAlgorithms as MOA
MoaMopb.run("tiny-inst.mopb", HiGHS.Optimizer, MOA.Lexicographic(), nothing, true, true, false, 1, true)
MoaMopb.run("tiny-inst.mopb", HiGHS.Optimizer, MOA.EpsilonConstraint(), nothing, true, true, true, 1, true)
MoaMopb.run("tiny-inst.mopb", HiGHS.Optimizer, MOA.Hierarchical(), nothing, true, true, true, 1, true)
MoaMopb.run("tiny-inst.mopb", HiGHS.Optimizer, MOA.TambyVanderpooten(), nothing, true, true, false, 1, true)
MoaMopb.run("tiny-inst.mopb", HiGHS.Optimizer, MOA.KirlikSayin(), nothing, true, true, false, 1, true)
MoaMopb.run("tiny-inst.mopb", HiGHS.Optimizer, MOA.DominguezRios(), nothing, true, true, false, 1, true)
import SCIP
MoaMopb.run("tiny-inst.mopb", SCIP.Optimizer, MOA.Lexicographic(), nothing, true, true, true, 1, true)
MoaMopb.run("tiny-inst.mopb", SCIP.Optimizer, MOA.EpsilonConstraint(), nothing, true, true, true, 1, true)
MoaMopb.run("tiny-inst.mopb", SCIP.Optimizer, MOA.Hierarchical(), nothing, true, true, true, 1, true)
MoaMopb.run("tiny-inst.mopb", SCIP.Optimizer, MOA.TambyVanderpooten(), nothing, true, true, false, 1, true)
MoaMopb.run("tiny-inst.mopb", SCIP.Optimizer, MOA.KirlikSayin(), nothing, true, true, false, 1, true)
MoaMopb.run("tiny-inst.mopb", SCIP.Optimizer, MOA.DominguezRios(), nothing, true, true, false, 1, true)
import Gurobi
MoaMopb.run("tiny-inst.mopb", Gurobi.Optimizer, MOA.Lexicographic(), nothing, true, true, true, 1, true)
MoaMopb.run("tiny-inst.mopb", Gurobi.Optimizer, MOA.EpsilonConstraint(), nothing, true, true, true, 1, true)
MoaMopb.run("tiny-inst.mopb", Gurobi.Optimizer, MOA.Hierarchical(), nothing, true, true, true, 1, true)
MoaMopb.run("tiny-inst.mopb", Gurobi.Optimizer, MOA.TambyVanderpooten(), nothing, true, true, false, 1, true)
MoaMopb.run("tiny-inst.mopb", Gurobi.Optimizer, MOA.KirlikSayin(), nothing, true, true, false, 1, true)
MoaMopb.run("tiny-inst.mopb", Gurobi.Optimizer, MOA.DominguezRios(), nothing, true, true, false, 1, true)
import CPLEX
MoaMopb.run("tiny-inst.mopb", CPLEX.Optimizer, MOA.Lexicographic(), nothing, true, true, true, 1, true)
MoaMopb.run("tiny-inst.mopb", CPLEX.Optimizer, MOA.EpsilonConstraint(), nothing, true, true, true, 1, true)
MoaMopb.run("tiny-inst.mopb", CPLEX.Optimizer, MOA.Hierarchical(), nothing, true, true, true, 1, true)
MoaMopb.run("tiny-inst.mopb", CPLEX.Optimizer, MOA.TambyVanderpooten(), nothing, true, true, false, 1, true)
MoaMopb.run("tiny-inst.mopb", CPLEX.Optimizer, MOA.KirlikSayin(), nothing, true, true, false, 1, true)
MoaMopb.run("tiny-inst.mopb", CPLEX.Optimizer, MOA.DominguezRios(), nothing, true, true, false, 1, true)
