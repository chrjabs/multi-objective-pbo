package org.sat4j.moco.goal_delimeter;

import org.sat4j.moco.util.Log;
import org.sat4j.specs.IVecInt;


import org.sat4j.moco.pb.PBSolver;
import org.sat4j.moco.problem.Instance;

public interface GoalDelimeterI{

    public void setSolver(PBSolver solver);
    public PBSolver getSolver();

    public void setInstance(Instance instance);
    public Instance getInstance();

    public int getFirstVariable();
    public void setFirstVariable();
    
    default public boolean UpdateCurrentK(int iObj, int upperKD){return true;};
    public boolean isY(int id);

    /**
     * Generate the upper limit assumptions
     */
    public IVecInt generateUpperBoundAssumptions(IVecInt explanation);
    public IVecInt generateUpperBoundAssumptions(IVecInt explanation, boolean whoKnows);

;

    /**
     *Checks if a variable is an X(original) variable.
     */
    default public boolean isX(int literal){
	PBSolver solver = this.getSolver();
	int id = solver.idFromLiteral(literal);
	//TODO: the 1 is the constant literal, that the solver creates.
	if(id < this.getFirstVariable() - 1)
	    return true;
	return false;
    }

    abstract public int getCurrentKD(int iObj);

    abstract public int getIObjFromY(int id);
    abstract public int getKDFromY(int id);

    abstract public int getY(int iObj, int iKD);

    abstract public String prettyFormatVariable(int literal);    

    default public String prettyFormatVecInt(IVecInt vecInt){
	 String result = "";
	 for(int j = 0; j < vecInt.size(); ++j)
	     result += this.prettyFormatVariable(vecInt.get(j));
	 return result;
     }
    default public String prettyFormatVecIntWithValues(IVecInt vecInt){
	 String result = "\n";
	 for(int j = 0; j < vecInt.size(); j++)
	     result += this.prettyFormatVariable(vecInt.get(j)) + " " +this.getSolver().modelValue(vecInt.get(j)) + "|\n";
	 return result;
     }
     
    default public String prettyFormatArrayWithValues(Integer[] literals){
	String result = "";
	for(int j = 0, n = literals.length; j < n ; ++j)
	    result += prettyFormatVariable(literals[j]) + " " +this.getSolver().modelValue(literals[j]) + "|\n";
	return result;
    }


    default public void prettyPrintVecInt(IVecInt vecInt, boolean clausing){
	if(clausing)
	    Log.clausing(this.prettyFormatVecInt(vecInt));
	else
	    Log.comment(6, this.prettyFormatVecInt(vecInt));
	return;
    }
    default public void prettyPrintVecInt(IVecInt vecInt){
	Log.comment(2, prettyFormatVecInt(vecInt));
	return;
     }
    default public void prettyPrintVecInt(IVecInt vecInt, int level){
	Log.comment(level, prettyFormatVecInt(vecInt));
	return;
     }
    default public void prettyPrintVariable(int literal){
	Log.comment(6,prettyFormatVariable(literal));
    }
    default public void prettyPrintVariable(int literal, int level){
	Log.comment(level,prettyFormatVariable(literal));
    }


    /**
     * Increase the upperKD, taking care of the semantics
     */
    public boolean preAssumptionsExtend(IVecInt currentExplanation);
    public int getUpperKD(int iObj);
    public void setUpperKD(int iObj, int upperKD);
}
