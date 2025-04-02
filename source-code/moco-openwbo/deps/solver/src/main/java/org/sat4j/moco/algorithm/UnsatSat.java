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
import java.nio.IntBuffer;
import java.util.Arrays;
import java.util.Comparator;
import org.sat4j.core.VecInt;
import org.sat4j.core.ReadOnlyVec;
import org.sat4j.core.ReadOnlyVecInt;
import org.sat4j.moco.analysis.Result;
import org.sat4j.moco.analysis.SubResult;
import org.sat4j.moco.goal_delimeter.GoalDelimeterI;
import org.sat4j.moco.util.Real;
import org.sat4j.moco.problem.Instance;
import org.sat4j.moco.problem.Objective;
import org.sat4j.moco.util.Log;
import org.sat4j.specs.ContradictionException;
import org.sat4j.specs.IVecInt;

/**
 * Class that implements UnsatSat, MSU3 flavoured
 * @author João Cortes
 */

public class UnsatSat extends algorithm implements IWithGoalDelimeter {

    /**
     * IDs of the variables used int the sequential encoder. The first
     * index is the goal, the second is the first index of s from " On
     * using Incremental Encodings...".Remember that s(i,j) is an
     * indicator of the propositions of the form x_i>=j.
     */

    private GoalDelimeterI goalDelimeter = null;

    public UnsatSat(){}
    /**
     * Creates an instance of a MOCO solver, for a given instance,
     * that applies the Pareto-MCS algorithm.
     * @param m The MOCO instance.
     */

    public UnsatSat(Instance m, boolean MSU3) {
        // Log.comment(3, "in UnsatSat constructor");
	this.problem = m;
	this.result = new Result(m, true);
	try {
            this.solver = buildSolver();
        }
        catch (ContradictionException e) {
            // Log.comment(3, "Contradiction in ParetoMCS.buildSolver");
            return;
        }
	
	this.subResult = new SubResult(this.problem);

    }

    public UnsatSat(Instance m) {
	this(m, false);
    }

    /**
     * Applies the UnsatSat algorithm to the MOCO instance provided
     */

    public void solve(){
	IVecInt currentExplanation = new VecInt(new int[] {});
	IVecInt currentAssumptions = new VecInt(new int[] {});

	boolean goOn = true;
	boolean goOn1 = true;
	Log.comment(2, "encoding setup completed.");
	currentAssumptions = this.generateUpperBoundAssumptions(currentExplanation, false);
	this.logUpperLimit();
	while(goOn){
	    Log.comment("");
	    Log.comment("new harvest cycle\n");
	    this.logUpperLimit();
	    Log.comment(6, "assumptions:");
	    Log.comment(6, this.prettyFormatVecInt(currentAssumptions));
	    solver.check(currentAssumptions);
	    this.solver.printStats();
	    if(goOn1 && solver.isSat()){
		this.saveModel();
		int[] diffAttainedValue = this.diffAttainedValue();
 		if(! this.blockDominatedRegion(diffAttainedValue)){
		    goOn1 = false;
		}

	    }else{
		finalizeHarvest();
		goOn = goOn1;
		if(goOn){
		    currentExplanation = solver.unsatExplanation();
		    Log.comment(2, "explanation length: " + currentAssumptions.size());
		    if(currentExplanation.size() == 0){
			goOn = false;
			Log.comment(6, "empty explanation");
		    }else{
			currentAssumptions = this.generateUpperBoundAssumptions(currentExplanation, true);
			// if currentAssumptions are null, then the
			// attainable domain did was not expanded and
			// there is no need to keep goind
			if(currentAssumptions == null){
			    Log.comment(6, "There was no expansion");
			    goOn = false;
			}else{
			    Log.comment(2, "explanation length: " + currentExplanation.size());

}
		    }
		}
	    }
	}
	this.result.setParetoFrontFound();

	return;
    }



    /**
     *Log the value of the upperLimit
     */

    private void logUpperLimit()    {
	String logUpperLimit = "diff upper limit: ["+this.getUpperKD(0);
	for(int iObj = 1; iObj < this.problem.nObjs(); ++iObj)
	    logUpperLimit +=", "+this.getUpperKD(iObj) ;//+ this.problem.getObj(iObj).getMinValue())
	//..log
	
	logUpperLimit +="]";
	Log.comment(2, logUpperLimit );
    }

    
    
    
    /**
     * Generate the upper limit assumptions
     */
    public IVecInt generateUpperBoundAssumptions(IVecInt explanation, boolean checkRange){
	IVecInt assumptions = new VecInt(new int[]{});
	assumptions = this.goalDelimeter.generateUpperBoundAssumptions(explanation, checkRange);
	return assumptions;
    }


    /**
     *gets the current upper limit of the explored value of the
     *differential k of the ithOjective
     *@param iObj
     */

    private int getUpperKD(int iObj){
	return this.goalDelimeter.getUpperKD(iObj);
    }


    /**
     *Checks if literal is an Y variable
     *@param literal
     */

    public boolean isY(int literal){
	if(this.goalDelimeter.isY(literal))
	    return true;
	return false;
    }
    /**
     *Checks if literal is an Y variable
     *@param literal
     */

    public boolean isX(int literal){
	return this.goalDelimeter.isX(literal);

    }

    /**
     *returns the model in DIMACS format, including only the real
     *variables and the Y variables of the sequential encoder
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
     *variables and the Y variables of the sequential encoder
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
     *variables and the Y variables of the sequential encoder
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
     *variables and the Y variables of the sequential encoder
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
	    Log.comment(5, "Model " + i);
	    this.printModel(models.get(i));
	}
	return;
    }

    // /**
    //  * Print an modelY
    //  * @param model,
    //  */

    public void printModelY(IVecInt modelY) {
    	int[][] convertedModel = new int[(modelY.size())][];
    	for(int i=0, n = modelY.size();i<n;i++){
    	    int yId = this.solver.idFromLiteral( modelY.get(i));
    	    int iObj = this.goalDelimeter.getIObjFromY(yId);
    	    int kD = this.goalDelimeter.getKDFromY(yId);
    	    convertedModel[i] = new int[]{  iObj, kD,modelY.get(i),};
    	}


    	Arrays.sort(convertedModel, Comparator.comparing(IntBuffer::wrap));

	String logYModel = "";
	int  currentIObj = convertedModel[0][0];
    	for(int i=0, n = convertedModel.length;i<n;i++){
	    if(convertedModel[i][0] == currentIObj)	 {
		// if(convertedModel[i][2] > 0)
		logYModel += this.goalDelimeter.prettyFormatVariable(convertedModel[i][2]) + " ";
	    }
	    else{
		Log.comment(3, logYModel);
		logYModel = "";
		currentIObj = convertedModel[i][0];
		Log.comment(5, "");
		i--;
	    }
	}
		Log.comment(3, logYModel);
    }

    /**
     * Print a model
     * @param models, the obtained models
     */

    public void printModel(IVecInt model) {
	this.goalDelimeter.prettyPrintVecInt(model, 3);
	}


    /**
     * The attained value of objective in the interpretation of the
     * last found model
     @param objective
    */
    private int attainedValue(Objective objective){
	int result = 0;
	int objectiveNLit = objective.getTotalLits();
	ReadOnlyVecInt objectiveLits = objective.getSubObjLits(0);
	ReadOnlyVec<Real> objectiveCoeffs = objective.getSubObjCoeffs(0);
	for(int iLit = 0; iLit < objectiveNLit; ++iLit  ){
	    int coeff = objectiveCoeffs.get(iLit).asInt();
	    int literal = objectiveLits.get(iLit);
	    if(this.solver.modelValue(literal))
		result += coeff;
	}
	return result;
    }

    public int[] diffAttainedValue(){
	int[] diffAttainedValue = new int[this.problem.nObjs()];
	for(int i = 0; i < this.problem.nObjs(); ++i){
	    diffAttainedValue[i] = this.attainedValue(this.problem.getObj(i));
	    diffAttainedValue[i]-=this.problem.getObj(i).getMinValue();
	}
	return diffAttainedValue;
    }

    /**
     * Block the region dominated by the last found model.
     */

    public boolean blockDominatedRegion(int[] diffAttainedValue ){
    
	String logDiffAttainedValue = "diff attained value: ["+ diffAttainedValue[0];
	for(int iObj = 1; iObj < this.problem.nObjs(); ++iObj)
	    logDiffAttainedValue +=", "+ diffAttainedValue[iObj];
	//..log
	
	logDiffAttainedValue +="]";
	Log.comment(2, logDiffAttainedValue );
	IVecInt newHardClause = new VecInt();
	for (int iObj = 0; iObj < this.problem.nObjs(); ++iObj){
	    if(diffAttainedValue[iObj] != 0){
		int possibleClause =- this.goalDelimeter.getY(iObj, diffAttainedValue[iObj]);
		//this better always be true.
		if(possibleClause != 0)
		    newHardClause.push(possibleClause);
	    }	}	    
	// Log.comment(6, "Blocking clause:");
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

    public void saveModel(){
	Log.comment(6, "model:");
	Log.comment(6, this.prettyFormatVecInt(this.getXModel()));
	this.subResult.saveModel(this.solver);
}
    public void finalizeHarvest(){
	this.transferSubResult();
}
	@Override
	public void setGoalDelimeter(GoalDelimeterI gd) {
	    this.goalDelimeter = gd;
		
	}
	@Override
	public GoalDelimeterI GetGoalDelimeter() {
	    return this.goalDelimeter;
	}

}
