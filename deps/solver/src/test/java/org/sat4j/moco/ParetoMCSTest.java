package org.sat4j.moco;

import org.sat4j.moco.algorithm.ParetoMCS;
import org.sat4j.moco.algorithmTest;


import static org.junit.Assert.assertTrue;

import org.junit.Test;
import org.sat4j.core.Vec;
import org.sat4j.core.VecInt;

import org.sat4j.moco.analysis.Result;
import org.sat4j.moco.pb.PBExpr;
import org.sat4j.moco.problem.DivObj;
import org.sat4j.moco.problem.Objective;
import org.sat4j.moco.util.Real;

public class ParetoMCSTest extends algorithmTest {

    public ParetoMCSTest(){};

    @Override
     protected ParetoMCS instateAlgorithm(){
	return new ParetoMCS(this.moco);

    }

    @Test
    public void testDivReduction() {
        PBExpr num = new PBExpr(new VecInt(new int[] { -2, 3 }),
                                new Vec<Real>(new Real[] { new Real(2), new Real(2) }));
        PBExpr den = new PBExpr(new VecInt(new int[] { 2, 3 }),
                                new Vec<Real>(new Real[] { Real.ONE, Real.ONE }));
        DivObj other_obj = new DivObj(new Vec<PBExpr>(new PBExpr[] { num }),
                                      new Vec<PBExpr>(new PBExpr[] { den }));
        this.moco.addObj(other_obj);
        this.solver = this.instateAlgorithm();
        this.solver.solve();
        Result result = this.solver.getResult();
        assertTrue(result.isParetoFront());
        assertTrue(result.nSolutions() == 2);
        boolean[][] front_sols = new boolean[][] { new boolean[] { false, true, true },
                                                   new boolean[] { true, true, false } };
        double[][] front_costs = new double[][] { new double[] { 1, 1 }, new double[] { 3, 0 } };
        Objective[] objs = new Objective[] { this.main_obj, other_obj };
        validateResult(result, objs, front_sols, front_costs);
    }
    
    @Test
    public void testSumOfDivReduction() {
        PBExpr num1 = new PBExpr(new VecInt(new int[] { -2, 3 }),
                                 new Vec<Real>(new Real[] { new Real(2), new Real(2) }));
        PBExpr den1 = new PBExpr(new VecInt(new int[] { 2, 3 }),
                                 new Vec<Real>(new Real[] { Real.ONE, Real.ONE }));
        PBExpr num2 = new PBExpr(new VecInt(new int[] { -1, 3 }),
                                 new Vec<Real>(new Real[] { new Real(2), new Real(2) }));
        PBExpr den2 = new PBExpr(new VecInt(new int[] { 1, 2 }),
                                 new Vec<Real>(new Real[] { Real.ONE, Real.ONE }));
        DivObj other_obj = new DivObj(new Vec<PBExpr>(new PBExpr[] { num1, num2 }),
                                      new Vec<PBExpr>(new PBExpr[] { den1, den2 }));
        this.moco.addObj(other_obj);
        this.solver = this.instateAlgorithm();
        this.solver.solve();
        Result result = this.solver.getResult();
        assertTrue(result.isParetoFront());
        assertTrue(result.nSolutions() == 2);
        boolean[][] front_sols = new boolean[][] { new boolean[] { false, true, true },
                                                   new boolean[] { true, true, false } };
        double[][] front_costs = new double[][] { new double[] { 1, 5 }, new double[] { 3, 0 } };
        Objective[] objs = new Objective[] { this.main_obj, other_obj };
        validateResult(result, objs, front_sols, front_costs);
    }
    

}
