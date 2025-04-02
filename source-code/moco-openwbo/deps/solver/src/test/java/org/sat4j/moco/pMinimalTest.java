package org.sat4j.moco;
import org.sat4j.moco.Launcher;
import org.sat4j.moco.util.IOUtils;
import java.io.IOException;
import org.sat4j.moco.parsing.OPBReader;

import org.sat4j.moco.algorithm.pMinimal;
import org.junit.AfterClass;
import org.junit.Test;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.sat4j.core.Vec;
import org.sat4j.core.VecInt;
import org.sat4j.moco.algorithm.algorithm;

import org.sat4j.moco.analysis.Result;
import org.sat4j.moco.pb.PBFactory;
import org.sat4j.moco.problem.LinearObj;
import org.sat4j.moco.problem.Instance;
import org.sat4j.moco.problem.Objective;
import org.sat4j.moco.util.Real;


public class pMinimalTest extends algorithmTest {
    public pMinimalTest(){};

    @Override
    public pMinimal instateAlgorithm(){
	return new pMinimal(this.moco);

    }
    @Test

    public void testBiOjective() {
	LinearObj other_obj = new LinearObj(new VecInt(new int[] { -2, 3 }),
                                            new Vec<Real>(new Real[] { new Real(2), new Real(2) }));
        this.moco.addObj(other_obj);
        assertTrue(this.moco.nObjs() == 2);
	assertTrue(this.moco.nConstrs() == 1);
        this.solver = this.instateAlgorithm();
        this.solver.solve();
        Result result = this.solver.getResult();
        assertTrue(result.isParetoFront());
        assertTrue(result.nSolutions() == 2);
        boolean[][] front_sols = new boolean[][] { new boolean[] { false, true, true },
                                                   new boolean[] { true, true, false } };
        double[][] front_costs = new double[][] { new double[] { 1, 2 }, new double[] { 3, 0 } };
        assertTrue(this.moco.nObjs() == 2);
        Objective[] objs = new Objective[] { this.main_obj, other_obj };
        validateResult(result, objs, front_sols, front_costs);
    }

@Test
    public void testBiOjectiveFromFile() {
	this.moco = new Instance();
	try{
	    OPBReader reader = IOUtils.mkFileReader("/home/joaooneillcortes/Mnemosyne/Aion/moco/mocoSource/solver/bugs/bugPMinimalUnit/minimal.opb");
	    this.moco  = reader.readMOCO();
	    
    }	
        catch (IOException e) {
            System.out.println("PARSER ERROR!");
            e.printStackTrace();
	    return;
        }

        this.solver = this.instateAlgorithm();
        this.solver.solve();
        Result result = this.solver.getResult();
        assertTrue(result.isParetoFront());
        assertTrue(result.nSolutions() == 2);
        boolean[][] front_sols = new boolean[][] { new boolean[] { false, true, true },
                                                   new boolean[] { true, true, false } };
        double[][] front_costs = new double[][] {  new double[] { 0, 1 } ,new double[] { -2, 3 } };
	assertTrue(this.moco.nObjs() == 2);
        Objective[] objs = new Objective[] { this.moco.getObj(0),  this.moco.getObj(1)};
        validateResult(result, objs, front_sols, front_costs);
    }
}
