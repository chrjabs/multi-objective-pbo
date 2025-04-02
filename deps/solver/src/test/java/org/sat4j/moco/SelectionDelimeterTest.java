/*******************************************************************************
 * SAT4J: a SATisfiability library for Java Copyright (C) 2004, 2012 Artois University and CNRS
 *
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the Eclipse Public License v1.0
 * which accompanies this distribution, and is available at
 *  http://www.eclipse.org/legal/epl-v10.html
 *
 * Alternatively, the contents of this file may be used under the terms of
 * either the GNU Lesser General Public License Version 2.1 or later (the
 * "LGPL"), in which case the provisions of the LGPL are applicable instead
 * of those above. If you wish to allow use of your version of this file only
 * under the terms of the LGPL, and not to allow others to use your version of
 * this file under the terms of the EPL, indicate your decision by deleting
 * the provisions above and replace them with the notice and other provisions
 * required by the LGPL. If you do not delete the provisions above, a recipient
 * may use your version of this file under the terms of the EPL or the LGPL.
 *
 * Contributors:
 *   CRIL - initial API and implementation
 *   Miguel Terra-Neves, Ines Lynce and Vasco Manquinho - MOCO solver
 *******************************************************************************/
package org.sat4j.moco;

import static org.junit.Assert.assertTrue;

import java.util.Iterator;

import org.junit.jupiter.api.Test;
import org.sat4j.core.Vec;
import org.sat4j.core.VecInt;
import org.sat4j.moco.pb.PBFactory;
import org.sat4j.moco.pb.PBSolver;
import org.sat4j.moco.problem.DigitalEnv;
import org.sat4j.moco.problem.DigitalEnv.DigitalNumber;
import org.sat4j.moco.problem.Instance;
import org.sat4j.moco.problem.LinearObj;
import org.sat4j.moco.problem.Objective;
import org.sat4j.moco.goal_delimeter.SelectionDelimeter;
import org.sat4j.moco.goal_delimeter.Circuit;
import org.sat4j.moco.goal_delimeter.Circuit.ControlledComponent;
import org.sat4j.moco.goal_delimeter.SelectionDelimeter.ObjManager;
import org.sat4j.moco.util.Log;
import org.sat4j.moco.util.MyModelIterator;
import org.sat4j.moco.util.Real;
import org.sat4j.specs.ContradictionException;
import org.sat4j.specs.IVecInt;

public class SelectionDelimeterTest {
    protected SelectionDelimeter sd = null;
    protected PBSolver pbSolver;
    protected Instance moco;
    protected IVecInt range;
    protected int[] upperBound;
    static{Log.setVerbosity(6);}

    public SelectionDelimeterTest(){};

    public void partialSetup(){
	try {
	    this.pbSolver = buildSolver();
	    
	}
	catch (ContradictionException e) {
	    Log.comment("Could not build the solver");
            return;
        }
	this.sd  = new SelectionDelimeter(moco, this.pbSolver,true){

		/**
		 *Adds the disjunction of setOfLiterals, and logs
		 *@param setOfliterals
		 */

		public boolean AddClause(IVecInt setOfLiterals){
		    this.prettyPrintVecInt(setOfLiterals, true);
		    try{
			this.solver.AddClause(setOfLiterals);
		    } catch (ContradictionException e) {
			Log.comment(2, "contradiction when adding clause: ");
			for(int j = 0; j < setOfLiterals.size(); ++j)
			    Log.comment(2, " " + setOfLiterals.get(j) + " " );
			return false;
		    }
		    return true;
		}

	    };
	this.range = new VecInt(this.pbSolver.nVars());
	for(int i = 0, n = this.pbSolver.nVars(); i < n; i++)
	    this.range.push(i + 1);

}

    /**
     * another incomplete instance, one goal function with 8 literals,
     * with weights equal to 1. No constraints
     */
    public void mocoSetup1() {
	this.moco = new Instance();
	// min: +1 x1 +1 x2 +1 x3 +1 x4 +1 x5 +1 x6 +1 x7 +1 x8 ;
	Objective main_obj = new LinearObj(new VecInt(new int[] { 1, 2, 3, 4, 5, 6, 7, 8 }),
					   new Vec<Real>(new Real[] {Real.ONE, Real.ONE, Real.ONE, Real.ONE, Real.ONE, Real.ONE, Real.ONE, Real.ONE }));
	this.moco.addObj(main_obj);
	this.partialSetup();
	this.upperBound = new int[]{3};
    }

    /**
     *Sets up the main moco instance. Only add the constraint if
     *{@code constraint}
     */
    public void mocoSetUp(boolean constraint) {
	this.moco = new Instance();

	// min: +2 x1 +9 x2 +5 x3 +7 x4
	// min: -8 x1 +3 x2 -1 x3 -1 x4
	// 4 x1 +1 x2 +3 x3 +2 x4 <= 5
	Objective main_obj = new LinearObj(new VecInt(new int[] { 1, 2, 3, 4 }),
				      new Vec<Real>(new Real[] { new Real(2), new Real(9), new Real(5), new Real(7) }));
	this.moco.addObj(main_obj);

	Objective second_obj =  new LinearObj(new VecInt(new int[] { 1,2,3,4}),
					      new Vec<Real>(new Real[] {new Real(-8), new Real(3), new Real(-1), new Real(-1) }));
	this.moco.addObj(second_obj);

	if(constraint)
	    this.moco.addConstr(PBFactory.instance().mkLE(new VecInt(new int[] { 1, 2, 3, 4 }),
							      new Vec<Real>(new Real[] { new Real(4), new Real(1), new Real(3), new Real(2) }),
							      new Real(5)));
	this.partialSetup();
    }

    @Test
    public void testDigits(){
	Integer[] ratios = new Integer[]{2};
	DigitalEnv digitalEnv = new DigitalEnv(ratios);
	DigitalEnv.DigitalNumber digitalNumber = digitalEnv.toDigital(9);
	int digit0 = digitalNumber.getDigit(1);
	assertTrue("digit0 fails", digit0 == 1);
	int digit1 = digitalNumber.getDigit(2);
	assertTrue("digit1 fails", digit1 == 0);
	int digit2 = digitalNumber.getDigit(4);
	assertTrue("digit2 fails", digit2 == 0);
	int digit3 = digitalNumber.getDigit(8);
	assertTrue("digit3 fails", digit3 == 1);
    }

    @Test
    public void testDigitsMulti(){
	Integer[] ratios = new Integer[]{2, 3};
	DigitalEnv digitalEnv = new DigitalEnv(ratios);
	DigitalEnv.DigitalNumber digitalNumber = digitalEnv.toDigital(8);
	int digit0 = digitalNumber.getDigitI(0);
	assertTrue("digit0 fails", digit0 == 0);
	int digit1 = digitalNumber.getDigitI(1);
	assertTrue("digit1 fails: " + digit1, digit1 == 1);
	int digit2 = digitalNumber.getDigitI(2);
	assertTrue("digit2 fails", digit2 == 1);
	int digit3 = digitalNumber.getDigitI(4);
	assertTrue("digit3 fails", digit3 == 0);
    }

    @Test
    public void testWeightExpansion(){
	Integer[] ratios = new Integer[]{2, 3};
	DigitalEnv digitalEnv = new DigitalEnv(ratios);
	DigitalEnv.DigitalNumber digitalNumber = digitalEnv.toDigital(8);
	int digit0 = digitalNumber.getDigitI(0);
	assertTrue("digit0 fails", digit0 == 0);
	int digit1 = digitalNumber.getDigitI(1);
	assertTrue("digit1 fails: " + digit1, digit1 == 1);
	int digit2 = digitalNumber.getDigitI(2);
	assertTrue("digit2 fails", digit2 == 1);
	int digit3 = digitalNumber.getDigitI(4);
	assertTrue("digit3 fails", digit3 == 0);
    }



    // @Test
    // public void digitalCounterTest(){
    // 	this.mocoSetUp(false);
    // 	int[] upperBound = new int[]{2, 2};
    // 	assertTrue(this.moco.nObjs() + " objectives and " + upperBound.length + "upper bounds", this.moco.nObjs() == upperBound.length);
    // 	IVecInt assumptions = this.sd.generateUpperBoundAssumptions(upperBound);
    // 	boolean[] inputValues = new boolean[this.moco.nVars()];
    // 	for(int i = 1, n = this.moco.nVars(); i <= n;i++)
    // 	    if(inputValues[i - 1])
    // 		assumptions.push(i);
    // 	    else
    // 		assumptions.push(-i);

    // 	Iterator<boolean[]> iterator = new MyModelIterator(this.pbSolver, assumptions);
    // 	boolean[] model;
    // 	boolean[]  saturatedComplete = new boolean[this.moco.nObjs()];
    // 	boolean[] saturated = new boolean[this.moco.nObjs()];
    // 	for(int i = 0, n = saturatedComplete.length; i < n; i++)
    // 	    saturatedComplete[i] = false;
    // 	int modelN = 0;
    // 	while(iterator.hasNext()){
    // 	    model = iterator.next();
    // 	    saturated = this.testDigitalValues(upperBound);
    // 	    for(int i = 0, n = saturated.length; i < n; i++)
    // 		saturatedComplete[i] = saturated[i] || saturatedComplete[i];
    // 	    modelN++;
    // 	    // this.sd.prettyFormatVecIntWithValues(this.range);
	
    // 	}
    // 	Log.comment(modelN + " models");
    // 	int obtained = 0;
    // 	// for(int i = 0, n = this.moco.nObjs(); i < n; i++){
    // 	//     assertTrue(i + "'th not saturated", saturated[i]);
    // 	// }
    // }
    private boolean[] testDigitalValues(int[] upperBound){
	this.mocoSetUp(false);
	boolean[] saturated = new boolean[upperBound.length];
	for(int i = 0, n = this.moco.nObjs(); i < n; i++){
	    ObjManager objManager  = this.sd.getIthObjManager(i);
	    Objective objective = this.moco.getObj(i);
	    int obtained = objective.evaluateDiff(this.pbSolver);
	    DigitalNumber digitalNumber = objManager.getDigitalEnv().toDigital(obtained);
	    DigitalNumber.IteratorJumps iterator = digitalNumber.iterator();
	    Circuit circuit = objManager.getCircuit();
	    int digit0 = 0;
	    int base = 0;
	    int lit = 0;
	    ControlledComponent contComp;
	    while(iterator.hasNext()){
		digit0 = iterator.next();
		base = iterator.currentBase();
		contComp = circuit.getControlledComponentBase(base);
		Log.comment("digit "+ base + " inputs:");
		System.out.print(this.sd.prettyFormatArrayWithValues(contComp.getInputs()));
		Log.comment("digit outputs:");
		System.out.println(this.sd.prettyFormatArrayWithValues(contComp.getOutputs()));
		assertTrue(obtained + " is > than " + upperBound[i],obtained <= upperBound[i]);
		if(obtained == upperBound[i])
		    saturated[i] = true;
	    }

	}
	return saturated;
    }
    // /**
    //  *Tests if all valid models respect the delimitation imposed by
    //  *{@code upperBound}.
    //  */
    // @Test
    // public void delimitationTest(){
    // 	this.mocoSetup1();
	
    // 	//last valid values
    // 	assertTrue(this.moco.nObjs() + " objectives and " + this.upperBound.length + "upper bounds", this.moco.nObjs() == this.upperBound.length);

    // 	IVecInt assumptions = this.sd.generateUpperBoundAssumptions(this.upperBound);
    // 	Iterator<boolean[]> iterator = new MyModelIterator(this.pbSolver, assumptions);
    // 	boolean[] model;
    // 	int modelN = 0;
    // 	Integer[] inputs = new Integer[this.moco.nVars()];
    // 	for(int i = 0; i < inputs.length; i++)
    // 	    inputs[i] = i + 1;
    // 	while(iterator.hasNext()){
    // 	    modelN++;
    // 	    // this.sd.prettyFormatVecIntWithValues(this.range);
    // 	    model = iterator.next();
    // 	    assertTrue("model \n" + this.sd.prettyFormatArrayWithValues(inputs) + "failed the test", this.testUpperBound(model, this.upperBound));}
    // 	Log.comment(modelN + " models");

    // }

    // private boolean testUpperBound(boolean[] model, int[] upperBound){
    // 	for(int i = 0, n = this.moco.nObjs(); i < n; i++)
    // 	    if(this.moco.getObj(i).evaluateDiff(model) <= upperBound[i])
    // 		continue;
    // 	    else
    // 		return false;
    // 	return true;	

    // }
    // @Test
    // public void InputAndDelimitationTest(){
    // 	this.mocoSetUp(true);
    // 	int[] upperBound = new int[]{1, 1};
    // 	IVecInt assumptions = this.sd.generateUpperBoundAssumptions(upperBound);
    // 	assertTrue(this.moco.nObjs() + " objectives and " + upperBound.length + "upper bounds", this.moco.nObjs() == upperBound.length);
    // 	assumptions.push(1);
    // 	assumptions.push(-2);
    // 	assumptions.push(-3);
    // 	assumptions.push(-4);
    // 	Iterator<boolean[]> iterator = new MyModelIterator(this.pbSolver, assumptions);
    // 	boolean[] model;
    // 	int modelN = 0;
    // 	while(iterator.hasNext()){
    // 	    modelN++;
    // 	    // this.sd.prettyFormatVecIntWithValues(this.range);
    // 	    model = iterator.next();
    // 	    assertTrue("model " + model + "failed the test", this.testUpperBound(model, upperBound));}
    // 	Log.comment(modelN + " models");

    // }


    /**
     * Creates a PB oracle initialized with the MOCO's constraints.
     * @return The oracle.
     * @throws ContradictionException if the oracle detects that the
     * MOCO's constraint set is unsatisfiable.
     */
    private PBSolver buildSolver() throws ContradictionException {
        PBSolver solver = new PBSolver();
        solver.newVars(this.moco.nVars());
        for (int i = 0; i < this.moco.nConstrs(); ++i) {
            solver.addConstr(this.moco.getConstr(i));
        }
        // Log.comment(5, "out UnsatSat.buildSolver");
	solver.setConstantID();
        return solver;
    }
}

