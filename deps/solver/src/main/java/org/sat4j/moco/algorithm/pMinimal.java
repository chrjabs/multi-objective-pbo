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
package org.sat4j.moco.algorithm;

import java.util.Vector;

import org.sat4j.core.VecInt;
import org.sat4j.moco.analysis.Result;
import org.sat4j.moco.goal_delimeter.GoalDelimeterI;
import org.sat4j.moco.goal_delimeter.SeqEncoder;
import org.sat4j.moco.problem.Instance;
import org.sat4j.moco.problem.Objective;
import org.sat4j.moco.util.Log;
import org.sat4j.specs.ContradictionException;
import org.sat4j.specs.IVecInt;

/**
 * Class that implements the p-minimal inspired algorithm
 * @author João Cortes
 */

public class pMinimal extends algorithm implements IWithGoalDelimeter  {

    /**
     * IDs of the variables used int the sequential encoder. The first
     * index is the goal, the second is the first index of s from " On
     * using Incremental Encodings...".Remember that s(i,j) is an
     * indicator of the propositions of the form x_i>=j.
     */

    private GoalDelimeterI goalDelimeter = null;

    /**
     * Last explored differential k, for each objective function.
     */
    private int[] UpperKD = null;
    /**
     *  Last id of the real, non auxiliary,  variables 
     */  
    private int realVariablesN = 0;


    /**
     * Creates an instance of a MOCO solver, for a given instance,
     * that applies the Pareto-MCS algorithm.
     * @param m The MOCO instance.
     */
    
    public pMinimal(Instance m) {
	Log.comment(5, "In pMinimal.pMinimal");
	this.problem = m;
	this.result = new Result(m, true);
	try {
            this.solver = buildSolver();
        }
        catch (ContradictionException e) {
            Log.comment(3, "Contradiction in ParetoMCS.buildSolver");
	    Log.comment(5, "}");
            return;
        }
	this.realVariablesN = this.solver.nVars();
	this.UpperKD =  new int[(this.problem.nObjs())];
	Log.comment(5, "}");
    }

    
    
    /**
     * Applies the p-minimal algorithm to the MOCO instance provided
     */

    public void solve() {
	// IVecInt currentYModel = new VecInt(new int[] {});
	// IVecInt currentXModel = new VecInt(new int[] {});
	for(int iObj = 0, nObj = this.problem.nObjs(); iObj < nObj; ++iObj)
	    this.setUpperKD(iObj);


	boolean[] currentXModelValues = new boolean[this.problem.nVars()];
	IVecInt assumptions = new VecInt(new int[] {});
	boolean sat = false;
        Log.comment(3, "{ pMinimal.solve");
	this.solver.check(assumptions);
	sat = this.solver.isSat();
	while(sat){
	    while(sat){		
		// currentYModel = this.getYModel();
		currentXModelValues = this.getXModelValues();
		this.setAssumptions(assumptions, currentXModelValues);
		this.solver.check(assumptions);
		sat = this.solver.isSat();
	    }
	    this.result.saveThisModelUnsafe(currentXModelValues);
	    sat = this.blockDominatedRegion(currentXModelValues);
	    if(sat){
		assumptions = new VecInt(new int[] {});
		this.solver.check();
		sat = this.solver.isSat();
	    }
	}
	this.result.setParetoFrontFound();
    }


    private void setAssumptions(IVecInt assumptions, boolean[] XModelValues){
	int[] upperLimits = this.findUpperLimits(XModelValues);
	boolean allZero = true;
	for(int iObj = 0, n = this.problem.nObjs(); iObj < n; ++iObj){
	    int lit = -this.goalDelimeter.getY(iObj, upperLimits[iObj]);
	    if(upperLimits[iObj] != 0)
		allZero = false;
	    if(lit != 0)
		assumptions.push(-this.goalDelimeter.getY(iObj, upperLimits[iObj]));
}
	if(allZero)
	    assumptions.push(this.solver.constantLiteral(false));
    }



    /**
     *Sets the current upper limit of the explored value of the
     *differential k of the ithOjective to newKD
     *@param iObj
     */
    private void setUpperKD(int iObj){
	int newKD = this.problem.getObj(iObj).getWeightDiff();
	if(this.goalDelimeter.getCurrentKD(iObj) < newKD)
	    this.goalDelimeter.UpdateCurrentK(iObj, newKD);
	this.UpperKD[iObj] = newKD;
    }

    
    /**
     *Checks if literal is an STop variable
     *@param literal
     */

    public boolean isY(int literal){
	if(this.goalDelimeter.isY(literal))
	    return true;
	return false;
    }
    /**
     *Checks if literal is an STop variable
     *@param literal
     */

    public boolean isX(int literal){
	int id = (literal>0)? literal: -literal;
	return id <= this.realVariablesN && id >= 1;

    }

    /**
     *returns the model in DIMACS format, including only the real
     *variables and the STop variables of the sequential encoder
     *@return a filtered model
     */

    public boolean[] getXModelValues(){
	boolean[] modelValues = new boolean[this.problem.nVars()];
	for(int id = 1; id <= this.problem.nVars();++id){
	    modelValues[id - 1] = this.solver.modelValue(id);
	}
	return modelValues;
    }
    /**
     *returns the model in DIMACS format, including only the real
     *variables and the STop variables of the sequential encoder
     *@return a filtered model
     */

    public IVecInt getYModel(){
	IVecInt model = new VecInt(new int[] {});
	for(int id = 1; id <= this.solver.nVars();++id){
	    int literal = (this.solver.modelValue(id))? id: -id;
	    if(this.isY(literal))
		model.push(literal);
	}
	return model;
    }
    /**
     *returns the model in DIMACS format, including only the real
     *variables and the STop variables of the sequential encoder
     *@return a filtered model
     */

    public IVecInt getXModel(){
	IVecInt model = new VecInt(new int[] {});
	for(int id = 1; id <= this.solver.nVars();++id){
	    int literal = (this.solver.modelValue(id))? id: -id;
	    if(this.isX(literal))
		model.push(literal);
	}
	return model;
    }
    /**
     *returns the model in DIMACS format, including only the real
     *variables and the STop variables of the sequential encoder
     *@return a filtered model
     */

    public IVecInt getFullModel(){
	IVecInt model = new VecInt(new int[] {});
	for(int id = 1; id <= this.solver.nVars();++id){
	    int literal = (this.solver.modelValue(id))? id: -id;
	    model.push(literal);
	}
	return model;
    }

    /** 
     * Print the models 
     * @param models, the obtained models
     */
    public void printModels(Vector<IVecInt> models) {
	for(int i = 0; i <models.size(); ++i){
	    System.out.println("Model " + i);
	    this.printModel(models.get(i));
	    System.out.println();
	}
	return;
    }

    /** 
     * Print a model 
     * @param models, the obtained models
     */
    public void printModel(IVecInt model) {
	for(int j = 0; j <model.size(); ++j)
	    this.goalDelimeter.prettyPrintVariable(model.get(j));
	System.out.println();


	return;
    }


    /**
     * The attained value of objective  in the interpretation of model 
     @param model
    */
    private int attainedValue(Objective objective, boolean[] XModelValues){
	return 	objective.evaluate(XModelValues).asIntExact();
    }
    /**
     * Block the region dominated by the known models.
     */

    public int[] findUpperLimits(boolean[] XModelValues){
	Log.comment(5, "In pMinimal.findUpperLimits");
	int[] upperLimits = new int[this.problem.nObjs()];
	for(int i = 0; i < this.problem.nObjs(); ++i){
	    upperLimits[i] = this.attainedValue(this.problem.getObj(i), XModelValues);
	    upperLimits[i]-=this.problem.getObj(i).getMinValue();
	}
	Log.comment(5, "}");
	return upperLimits;
    }

    public boolean blockDominatedRegion(boolean[] XModelValues){
	Log.comment(5, "{ pMinimal.blockDominatedregion");
	int[] upperLimits = this.findUpperLimits(XModelValues);
	IVecInt newHardClause = new VecInt();
	for (int iObj = 0; iObj < this.problem.nObjs(); ++iObj){
	    int lit = -this.goalDelimeter.getY(iObj, upperLimits[iObj]);
	    if(lit != 0)
		newHardClause.push(-this.goalDelimeter.getY(iObj, upperLimits[iObj]));
	}
	Log.comment(5, "}");
	return this.AddClause(newHardClause);
    }


    public boolean blockModelX(IVecInt modelX){
	IVecInt notPreviousModel = new VecInt(new int[] {});
	for(int iX = 0; iX < modelX.size(); ++iX)
	    notPreviousModel.push(-modelX.get(iX));
	return this.AddClause(notPreviousModel);
    }

    public void prettyPrintVecInt(IVecInt vecInt, boolean clausing){
	if(clausing)
	    Log.clausing(this.goalDelimeter.prettyFormatVecInt(vecInt));
	else
	    Log.comment(6, this.goalDelimeter.prettyFormatVecInt(vecInt));
	return;
    }
    
    public String prettyFormatVecInt(IVecInt literals){return this.goalDelimeter.prettyFormatVecInt(literals);}
    public void printFlightRecordParticular(){

    }



	@Override
	public void setGoalDelimeter(GoalDelimeterI gd) {
	    this.goalDelimeter = gd;
		
	}



	@Override
	public GoalDelimeterI GetGoalDelimeter() {
		return goalDelimeter;
	}
}

