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
import org.sat4j.moco.util.MyModelIterator;
import java.util.List;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.Random;

import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Disabled;
import org.junit.jupiter.api.Test;
import org.sat4j.core.VecInt;
import org.sat4j.moco.algorithm.UnsatSat;
import org.sat4j.moco.pb.PBSolver;
import org.sat4j.moco.goal_delimeter.SelectionDelimeter;
import org.sat4j.moco.goal_delimeter.Circuit;
import org.sat4j.moco.goal_delimeter.Circuit.ControlledComponent;
import org.sat4j.moco.util.General;
import org.sat4j.moco.util.Log;
import org.sat4j.specs.ContradictionException;
import org.sat4j.specs.IVecInt;

public class CircuitTest {
    protected SelectionDelimeter sd = null;
    protected UnsatSat solver;
    protected PBSolver pbSolver = new PBSolver();
    static {Log.setVerbosity(6);}

    @BeforeEach
    public void partialSetUp() {
	pbSolver = new PBSolver();
	this.pbSolver.setConstantID();


    }
    
    /**
     *Given a sorted unary input, check if the ModComponent is
     *semi-correct. 
     */
    @Disabled
    @Test
    public void ModComponentTest(){
	int modN = 3;
	Integer[] inputs = new Integer[8];
	this.fillInputWithVars(inputs);

	Circuit circuit = new Circuit(pbSolver){
		public void buildCircuit(){
		    ModComponent modComponent = new ModComponent(inputs, modN);
		    modComponent.constitutiveClause();
		    new ControlledComponent(0, modComponent);
		}
		public int getFreshVar1(){pbSolver.newVar();return pbSolver.nVars();}

		public boolean AddClause1(IVecInt setOfLiterals){
		    return AddClause(setOfLiterals, true);
		}
	    };
	circuit.buildCircuit();

	
	ControlledComponent comp = circuit.getControlledComponentBase(0);
	Iterator<boolean[]> iterator = new MyModelIterator(this.pbSolver);
	boolean[] model;
	while(iterator.hasNext()){
	    model = iterator.next();
	    this.testModComponentModel(comp, modN);
	}
    }


    /**
     *Assertion helper of {@code  ModComponentTest()}.
     */

    private void testModComponentModel(ControlledComponent comp, int modN){
	Integer[] inputs = comp.getInputs();
	Integer[] outputs = comp.getOutputs();
	int value = 0;
	if(!this.valuesAreSorted(inputs))
	    return;
	for(int lit: inputs)
	    if(this.pbSolver.modelValue(lit))
		value++;
	    else
		break;
	value = value % modN;
	for(int i = 0; i < value; i++){
	    if(!this.pbSolver.modelValue( outputs[i])){
		Log.comment("inputs of ModComponent:");
		General.FormatArrayWithValues(inputs, pbSolver, true);
		Log.comment("outputs of ModComponent:");
		General.FormatArrayWithValues(outputs, pbSolver, true);
		assertTrue("Failing at " + i +"'th comparison", false);
	    }
	}
}



    /**
     *Checks if the i'th output of SelectionComponent is active if the
     *i+1'th value of the sorted inputs is active.
     */
    @Disabled
    @Test
    public void SelectionComponentTest(){
	int inputsLength = 8;
	int sortedPortionN = inputsLength;
	Random rand = new Random();
	Integer[] inputs = new Integer[inputsLength];
	for(int i = 0; i < inputs.length; i++){
	    this.pbSolver.newVar();
	    inputs[i] =  this.pbSolver.nVars();
	}

	Circuit circuit = new Circuit(this.pbSolver){
		public void buildCircuit(){
		    SelectionComponent comp = new SelectionComponent(inputs, sortedPortionN);
		    comp.constitutiveClause();
		    new ControlledComponent(0, comp);

		}
		public int getFreshVar1(){pbSolver.newVar();return pbSolver.nVars();}

		public boolean AddClause1(IVecInt setOfLiterals){
		    return AddClause(setOfLiterals, true);
		}
	    };

	circuit.buildCircuit();

	ControlledComponent controlledComp =  circuit.getControlledComponentBase(0);
	Integer[] outputs = controlledComp.getOutputs();
	// assumptions.push(-controlledComp.getIthOutput(2));

	MyModelIterator iterator = new MyModelIterator(pbSolver);
	while(iterator.hasNext()){
	    iterator.next();
	    this.SelectionComponentAssertion(inputs, outputs, sortedPortionN);
	}
    }

    /**
     *Assertion helper of {@code  SelectionComponentTest()}.
     */

    public void SelectionComponentAssertion(Integer[] inputs, Integer[] outputs, int sortedPortionN){
	int value = 0;

	for(int i = 0, n = inputs.length ; i < n; i++)
	    if(pbSolver.modelValue(inputs[i]))
		value++;

	if(value > sortedPortionN) value = sortedPortionN;
	for(int i = 0, n = value ; i < n; i++){
	    if(!pbSolver.modelValue(outputs[i])){
		Log.comment(3, "inputs:");
		General.FormatArrayWithValues(inputs, pbSolver,true);
		Log.comment(3, "outputs:");
		General.FormatArrayWithValues(outputs, pbSolver,true);
		assertTrue("Failing at " + i + "'th output", false);
	    }
	}
    }

    /**
     *Checks if negating the SelecComp.outputs[i] literal delimits the
     *number of active inputs to i+1;
     */
    @Disabled
    @Test
    public void SelectionComponentDelimetedTest(){

	int inputsLength = 8;
	int sortedPortionN = inputsLength;
	int value = 3;
	Random rand = new Random();
	Integer[] inputs = new Integer[inputsLength];
	Integer[] inputValues = new Integer[inputsLength];
	Arrays.fill(inputValues, 0);
	for(int i = 0; i < inputs.length; i++){
	    this.pbSolver.newVar();
	    inputs[i] =  this.pbSolver.nVars();
	}

	IVecInt assumptions = new VecInt();
	Circuit circuit = new Circuit(this.pbSolver){
		public void buildCircuit(){
		    SelectionComponent comp = new SelectionComponent(inputs, sortedPortionN);
		    comp.constitutiveClause();
		    new ControlledComponent(0, comp);

		}
		public int getFreshVar1(){pbSolver.newVar();return pbSolver.nVars();}

		public boolean AddClause1(IVecInt setOfLiterals){
		    return AddClause(setOfLiterals, true);
		}
	    };

	circuit.buildCircuit();
	ControlledComponent controlledComp =  circuit.getControlledComponentBase(0);
	assumptions.push(-controlledComp.getIthOutput(2));

	MyModelIterator iterator = new MyModelIterator(pbSolver, assumptions);
	while(iterator.hasNext()){
	    iterator.next();
	    this.SelectionComponentDelimetedAssertion(controlledComp, value);
	}
    }

    /**
     *Assertion helper of {@code SelectionComponentDelimetedTest()}.
     */

    private void SelectionComponentDelimetedAssertion(ControlledComponent comp, int value){
	Integer[] inputs = comp.getInputs();
	int result = 0;
	for(int i = 0, n = inputs.length ; i < n; i++){
	    if(pbSolver.modelValue(inputs[i]))
		result++;
	}
	if(!(result<=value)){
	    General.FormatArrayWithValues(inputs, pbSolver,true);
	    General.FormatArrayWithValues(comp.getOutputs(), pbSolver,true);
	    assertTrue("not delimeted", result < value);
	}
    }

    /**
     *Checks if the output of MergComponent is sorted. Notice each of
     *the 4 (sub)lits in inputsArray must be sorted.
     */
    @Disabled
    @Test void MergeComponentTest(){
	Integer[] inputs = new Integer[16];
	int sortedPortionN = inputs.length;
	this.fillInputWithVars(inputs);
	Integer[][] inputsArray = new Integer[4][];
	int n = inputs.length;
	{
	    int i = 0;
	    int m = 0;
	    int i0 = 0;
	    for(; m < 4 - 1; m++){
		inputsArray[m] = new Integer[n / 4];
		i0 = i;
		for(; i < (m + 1) * n / 4; i++)
		    inputsArray[m][i - i0 ] = inputs[i];
	    }
	    inputsArray[m] = new Integer[n - 3 * n / 4];
	    i0 = i;
	    for(; i < n; i++)
		inputsArray[m][i - i0] = inputs[i];
	}

	Circuit circuit = new Circuit(this.pbSolver){
		public void buildCircuit(){
		    MergeComponent comp = new MergeComponent(sortedPortionN);
		    comp.constitutiveClause(inputsArray);
		    new ControlledComponent(0, comp);
		}


		public int getFreshVar1(){pbSolver.newVar();return pbSolver.nVars();}

		public boolean AddClause1(IVecInt setOfLiterals){
		    return AddClause(setOfLiterals, true);
		}
	    };

	circuit.buildCircuit();
	ControlledComponent comp1 = circuit.getControlledComponentBase(0);

	Iterator<boolean[]> iterator = new MyModelIterator(this.pbSolver);
	boolean[] model;
	while(iterator.hasNext()){
	    model = iterator.next();
	    this.MergeComponentAssertion(comp1, inputsArray);
	}
	    return;

    }
    private void MergeComponentAssertion(ControlledComponent comp, Integer[][] inputsArray){
	int value = 0;
	for(Integer[] inputs: inputsArray)
	    value += this.activePrefixSize(inputs);
	if(!(this.activePrefixSize(comp.getOutputs()) >= value)){
	    Log.comment("inputs:");
	    General.FormatArrayWithValues(comp.getInputs(),this.pbSolver ,true);
	    Log.comment("outputs:");
	    General.FormatArrayWithValues(comp.getOutputs(),this.pbSolver ,true);
	assertTrue("output is not semi-sorted", false);
	}
}
    

    /**
     *Checks if the optimum component is correct. Does not enumerate all models
     */
    @Test
    void OptimumComponentTest(){
	Integer[] inputs = new Integer[10];
	boolean polarity = true;
	this.fillInputWithVars(inputs);
	Circuit circuit = new Circuit(pbSolver){
		public void buildCircuit(){
		    optimumComponent comp = new optimumComponent(inputs, polarity);
		    comp.constitutiveClause();
		    new ControlledComponent(0, comp);
		}


		public int getFreshVar1(){pbSolver.newVar();return pbSolver.nVars();}

		public boolean AddClause1(IVecInt setOfLiterals){
		    return AddClause(setOfLiterals, true);
		}
	    };
	circuit.buildCircuit();
	ControlledComponent comp = circuit.getControlledComponentBase(0);
	Iterator<boolean[]> iterator = new MyModelIterator(this.pbSolver);
	boolean[] model;
	assertTrue("length is not 1", comp.getOutputs().length == 1);
	while(iterator.hasNext()){
	    model = iterator.next();
	    this.optimumComponentAssertion(comp, polarity);
	}

    }

    public boolean optimumComponentAssertion(ControlledComponent comp, boolean polarity){
	boolean expected = !polarity;
	Integer[] inputs = comp.getInputs();
	int output = comp.getIthOutput(0);
	for(int i = 0, n = inputs.length; i < n; i++)
	    if((!this.pbSolver.modelValue(inputs[i]) || polarity) && (this.pbSolver.modelValue(inputs[i]) || !polarity)){
		expected = polarity;
		break;
	    }

	assertTrue("value is not correct", (!this.pbSolver.modelValue(output) || expected) && (this.pbSolver.modelValue(output) || !expected));
	return true;
	
}

    /**
     *Checks if the output of CombineComponent is sorted
     */
    @Test
    public void CombineComponentTest(){
	Random rand = new Random();
	Integer[] inputsSize = new Integer[]{4,8};
	int sortedPortionN = inputsSize[0] + inputsSize[1];
	Integer[][] inputsArray = new Integer[2][];
	for(int k = 0; k < 2; k++){
	    inputsArray[k] = new Integer[inputsSize[k]];
	    this.fillInputWithVars(inputsArray[k]);
}

	Circuit circuit = new Circuit(pbSolver){
		public void buildCircuit(){
		    Circuit.CombineComponent comp = new Circuit.CombineComponent(sortedPortionN);
		    comp.constitutiveClause(inputsArray[0], inputsArray[1]);
		    new ControlledComponent(0, comp);
		}


		public int getFreshVar1(){pbSolver.newVar();return pbSolver.nVars();}

		public boolean AddClause1(IVecInt setOfLiterals){
		    return AddClause(setOfLiterals, true);
		}
	    };

	circuit.buildCircuit();
	ControlledComponent comp1 = circuit.getControlledComponentBase(0);
	comp1.getOutputs();
	Iterator<boolean[]> iteratorModels = new MyModelIterator(this.pbSolver);
	boolean[] model;
	int modelN = 0;

	while(iteratorModels.hasNext()){
	    model = iteratorModels.next();
	    if(this.CombineComponentAssertion(comp1, inputsArray))
		modelN++;
	    // for(int i = 0, n = semiSorted.length; i < n; i++  ){
	    // 	int lit = comp1.getIthOutput(i);
	    // 	if(semiSorted[i] == 1)
	    // 	    assertTrue("failing " + i +"'th comparison", this.pbSolver.modelValue(lit));
	    // 	else
	    // 	    break;
	    // }
	}
	    Log.comment(modelN + " assertable models");

    }
    private boolean CombineComponentAssertion(ControlledComponent comp, Integer[][] inputsArray){
	for(Integer[] inputs: inputsArray)
	    if(!this.valuesAreSorted(inputs))
		return true;
	if(!(this.activePrefixSize(inputsArray[0]) > this.activePrefixSize(inputsArray[1])))
	    return true;
	int expected = this.activePrefixSize(inputsArray[0]) + this.activePrefixSize(inputsArray[1]);
	if(!(this.valuesAreSorted(comp.getOutputs()) && (this.activePrefixSize(comp.getOutputs()) >= expected))){
	    Log.comment("inputs:");
	    General.FormatArrayWithValues(comp.getInputs(),this.pbSolver ,true);
	    Log.comment("outputs:");
	    General.FormatArrayWithValues(comp.getOutputs(),this.pbSolver ,true);
	    assertTrue("output is not sorted", false);
	    return false;
	}
	return true;
    }

    /**
     * Checks if {@code upperValue} delimits the inputs accordingly,
     * that is, if n is the count of active inputs, then the unary
     * output cannot exceed mod2 upperValue.
     */
    @Disabled
    @Test void DigitComponentTest(){
	int modN = 2;
	Integer[] inputs = new Integer[6];
	
	for(int i = 0; i < inputs.length; i++){
	    this.pbSolver.newVar();
	    inputs[i] =  this.pbSolver.nVars();
	}
	Circuit circuit = new Circuit(pbSolver){
		public void buildCircuit(){
		    DigitComponent digitComp = new DigitComponent(inputs, modN);
		    digitComp.constitutiveClause();
		    new ControlledComponent(0, digitComp);
		}
		public int getFreshVar1(){pbSolver.newVar();return pbSolver.nVars();}

		public boolean AddClause1(IVecInt setOfLiterals){
		    // this.prettyPrintVecInt(setOfLiterals,true);
		    try{
			pbSolver.AddClause(setOfLiterals);
		    } catch (ContradictionException e) {return false;} return true;
		}
	    };
	circuit.buildCircuit();
	
	ControlledComponent comp = circuit.getControlledComponentBase(0);
	Iterator<boolean[]> iterator = new MyModelIterator(this.pbSolver);
	boolean[] model;
	while(iterator.hasNext()){
	    model = iterator.next();
	    this.testDigitComponentModel(comp, modN);
	}
    }

    /**
     *Assertion helper of {@code DigitComponentTest()}.
     */
    private void testDigitComponentModel(ControlledComponent comp , int modN){

	Integer[] inputs = comp.getInputs();
	Integer[] outputs = comp.getOutputs();
	int value = 0;

	for(int lit: inputs)
	    if(this.pbSolver.modelValue(lit))
		value++;

	value = value % modN;
	for(int i = 0; i < value; i++)
	    if(!this.pbSolver.modelValue(outputs[i])){
		Log.comment("inputs of DigitComponent:");
		General.FormatArrayWithValues(comp.getInputs(), pbSolver, true);
		Log.comment("outputs of DigitComponent:");
		General.FormatArrayWithValues(comp.getOutputs(), pbSolver, true);
		assertTrue("Failing at " + i +"'th comparison", false);
	    }
    }


    public IVecInt buildAssumption(Integer[] inputValues, Integer[] inputs){
	IVecInt assumptions = new VecInt();
	for(int i = 0, n = inputValues.length; i < n; i++){
	    if(inputValues[i] == 1)
	        assumptions.push(inputs[i]);
	    else
	        assumptions.push(-inputs[i]);
	}
	return assumptions;
    }

    public boolean AddClause(IVecInt setOfLiterals, boolean print){
	int lit;
	if(print)
	    {
		Log.comment("clause:");
		for(int i = 0, n = setOfLiterals.size(); i < n; i++){
		    lit = setOfLiterals.get(i);
		    Log.comment("lit:"+ lit);
		}
		
		Log.comment("");
	    }
	try{
	    pbSolver.AddClause(setOfLiterals);
	} catch (ContradictionException e) {return false;} return true;
    }
    public void fillInputWithVars(Integer[] inputs){
	for(int i = 0; i < inputs.length; i++){
	    this.pbSolver.newVar();
	    inputs[i] =  this.pbSolver.nVars();
	}

    }
    public boolean valuesAreSorted(Integer[] inputs){
	int i = 0, n =inputs.length;
	for(; i < n; i++)
	    if(!this.pbSolver.modelValue(inputs[i]))
		break;
	for(; i < n; i++)
	    if(this.pbSolver.modelValue(inputs[i]))
		return false;
	return true;
    }

    public int activePrefixSize(Integer[] literals){
	int result = 0;
	for(int lit: literals)
	    if(this.pbSolver.modelValue(lit))
		result++;
	return result;
    }

}


