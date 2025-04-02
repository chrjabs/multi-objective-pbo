package org.sat4j.moco.goal_delimeter;

import java.util.HashMap;

import org.sat4j.core.ReadOnlyVec;
import org.sat4j.core.ReadOnlyVecInt;
import org.sat4j.core.VecInt;
import org.sat4j.moco.pb.PBSolver;
import org.sat4j.moco.problem.Instance;
import org.sat4j.moco.problem.Objective;
import org.sat4j.moco.util.Log;
import org.sat4j.moco.util.Real;
import org.sat4j.specs.ContradictionException;
import org.sat4j.specs.IVecInt;


public abstract class GoalDelimeter<PIndex extends Index> implements GoalDelimeterI{

    /**
     * the instance to be solved
     */
    protected Instance instance = null;
    /**
     * the solver being used
     */
    protected PBSolver solver = null;

    protected Librarian<PIndex> librarian = null;

    /**
     * Last explored differential k, for each objective function.
     */
    protected int[] upperKD = null;

    /**
     *First solver variable that pertains to the goal delimeter
     *encoding
     */

    protected int firstVariable = 0;

    public GoalDelimeter(){
	this.librarian = new Librarian<PIndex>();
    };


    public GoalDelimeter(Instance instance, PBSolver solver) {
	this();
	this.instance = instance;	
	this.solver = solver;
	this.setFirstVariable();
	this.upperKD = new int[this.instance.nObjs()];
}

    public void setSolver(PBSolver solver){this.solver = solver;}
    public PBSolver getSolver(){return this.solver;}

    public void setInstance(Instance instance){this.instance = instance;}
    public Instance getInstance(){return this.instance;}

    public int getFirstVariable(){return this.firstVariable;}
    public void setFirstVariable(){this.firstVariable = this.getSolver().nVars() + 1;}
    
    /**
     * Generate the upper limit assumptions
     */

    public IVecInt generateUpperBoundAssumptions(IVecInt explanation, boolean checkChange){
	if(!this.preAssumptionsExtend(explanation) & checkChange)
	    return null;
	return this.generateUpperBoundAssumptions(explanation);
    }


    public IVecInt generateUpperBoundAssumptions(IVecInt explanation){
	
	IVecInt assumptions = new VecInt(new int[]{});
	
	for(int iObj = 0; iObj < this.instance.nObjs(); ++iObj){
	    int IthUpperBound = this.nextKDValue(iObj, getUpperKD(iObj));
	    Objective ithObjective = this.instance.getObj(iObj);
	    if(this.getUpperKD(iObj)  != IthUpperBound){
		int newY = -this.getY(iObj, IthUpperBound);
		if(newY!=0)
		    assumptions.push(newY);
	    }

	    ReadOnlyVecInt objectiveLits = ithObjective.getSubObjLits(0);
	    ReadOnlyVec<Real> objectiveCoeffs = ithObjective.getSubObjCoeffs(0);
	    int sign;
	    int ithAbsoluteWeight;

	    for(int iX = 0, nX = ithObjective.getTotalLits(); iX <nX; iX ++){
		int ithX = objectiveLits.get(iX);
		ithAbsoluteWeight = objectiveCoeffs.get(iX).asInt();
		sign = (ithAbsoluteWeight > 0? 1 : -1);
		ithAbsoluteWeight *= sign;
		if(ithAbsoluteWeight > getUpperKD(iObj))
		    assumptions.push(-sign * ithX);
	    }

	}

	return assumptions;
    }

    /**
     *Checks if a variable is an X(original) variable.
     */
    public boolean isX(int literal){
	int id = this.solver.idFromLiteral(literal);
	//TODO: the 1 is the constant literal, that the solver creates.
	if(id < this.firstVariable - 1)
	    return true;
	return false;
    }


    /**
     *Adds the disjunction of setOfLiterals
     *@param setOfliterals
     */

    public boolean AddClause(IVecInt setOfLiterals){
	// Log.comment(6, "clause to add:");
	// Log.comment(6, this.prettyFormatVecInt(setOfLiterals));

	try{
	    this.solver.AddClause(setOfLiterals);
	} catch (ContradictionException e) {
	    Log.comment(6, "contradiction when adding clause: ");
	    for(int j = 0; j < setOfLiterals.size(); ++j)
		Log.comment(3, " " + setOfLiterals.get(j) + " " );
	    return false;
	}
	return true;
    }
	public int generateNext(int iObj, int kD, int inclusiveMax) {
	    if(kD < inclusiveMax)
		return kD +1;
	    if(kD == inclusiveMax)
		return kD;
		return 0;
	}

	public int nextKDValue(int iObj, int kD) {
	    int max = this.instance.getObj(iObj).getWeightDiff();
	    if(kD < max)
		return  kD + 1;
	    if(kD == max)
		return kD;
	    else return -1;
	}

    public int getUpperKD(int iObj){
	return this.upperKD[iObj];
    }

    public void setUpperKD(int iObj, int newKD){
	if(this.getUpperKD(iObj)< newKD)
	    this.upperKD[iObj] = newKD;
    }

    /**
     *Increase the upperKD, taking care of the semantics
     *
     */

    public boolean preAssumptionsExtend(IVecInt currentExplanation){
	// Log.comment(0, "covered x variables: " + this.coveredLiterals.size());
	boolean change = false;
	HashMap<Integer,Boolean> objectivesToChange = new HashMap<Integer, Boolean>(this.instance.nObjs());
	for(int lit: currentExplanation.toArray()){
	    int id = this.solver.idFromLiteral(lit);
	    if(this.isX(id)){
		for(int iObj = 0; iObj < this.instance.nObjs(); ++iObj){
		    if(this.instance.getObj(iObj).getSubObj(0).weightFromLit(id) != null)
			objectivesToChange.put(iObj, null);
		}
	    }
	    else
		objectivesToChange.put(this.getIObjFromY(id), null);
	    for(int iObj :objectivesToChange.keySet()){
		// Log.comment(3, "changing upperlimit " + iObj);
		int upperKDBefore = this.getUpperKD(iObj);
		this.generateNext(iObj, upperKDBefore, this.instance.getObj(iObj).getWeightDiff());
		this.setUpperKD(iObj, this.nextKDValue(iObj, this.getUpperKD(iObj)));
		if(this.getUpperKD(iObj)!= upperKDBefore)
		    change = true;
	    }
	}
	return change;
    }

    /**
     *Pretty print the variable in literal. 
     */
    public String prettyFormatVariable(int literal){
	int sign =(literal>0)? 1: -1;
	int id =  literal * sign;

	if(isX(id)){
	    return (sign>0? "+":"-")+"X["+id+"] ";
	}
	if(this.isY(id)){
	    int iObj = this.getIObjFromY(id);
	    int kD = this.getKDFromY(id);
	    int k = kD; // + this.instance.getObj(iObj).getMinValue();
	    return "Y[" + iObj + ", " + k +"]"+ "::" + literal + " ";
	}
	return literal + " ";
    }

}

