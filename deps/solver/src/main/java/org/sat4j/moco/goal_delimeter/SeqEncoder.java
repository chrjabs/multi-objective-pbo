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
 *   João O'Neill Cortes, INESC
 *******************************************************************************/
package org.sat4j.moco.goal_delimeter;

import org.sat4j.moco.util.Log;
import java.lang.Math;
import java.util.Arrays;
import java.util.ArrayList;
import java.util.List;
import java.util.Hashtable;
import org.sat4j.core.ReadOnlyVec;
import org.sat4j.core.ReadOnlyVecInt;
import org.sat4j.core.VecInt;
import org.sat4j.moco.util.Real;
import org.sat4j.moco.pb.PBSolver;
import org.sat4j.moco.problem.Instance;
import org.sat4j.moco.problem.Objective;
import org.sat4j.specs.IVecInt;
import org.sat4j.specs.ContradictionException;

/**
 * Class with the implementation of the sequetial encoder.
 * Notice that the differential is both a value and an index starting at 0
 * @author Joao O'Neill Cortes
 */

 public class SeqEncoder extends GoalDelimeter<Index> {

    /** 
     * IDs of the S(equential) variables used to enforce the semantics of the sequential encoder.
     * idsS[iObj][iX][iKD] 
     */
    private List<ArrayList<ArrayList<Integer>>>  idsS = null;

     /** 
      * IDs of the B(locking) variables used to allow incrementality
      */
     private int[][] idsB = null;
     
    // private ConstrID topConstraint = null;

    /** 
     * Current top initialized differential k's for all objective functions
     */

    private int[] currentKDs = null;
    /** 
     * Current top differential k's for all objective functions
     */
    private int[] initializedKDs = null;

    /**
     *The inverse index map for the S variables
     */
    private Hashtable<Integer,int[]> sTopVariablesInverseIndex  = new Hashtable<Integer, int[]>();

    /**
     *The inverse index map for the S variables
     */
    private Hashtable<Integer,int[]> auxVariablesInverseIndex  = new Hashtable<Integer, int[]>();
     /**
      *The inverse index map for the blocking variables
      */
     private Hashtable<Integer, Integer> bVariablesInverseIndex  = new Hashtable<Integer, Integer>();


    /**
     * Creates an Instance of the sequential encoder
     * @param instance, the pseudo boolean instance
     * @param solver, the solver to be updated
     */

    public SeqEncoder(Instance instance, PBSolver solver) {
	Log.comment(4, "{ SeqEncoder");
	this.instance = instance;	
	this.solver = solver;
	this.firstVariable = solver.nVars() + 1;
	this.initializeIdsS();
	this.initializeIdsB();
	this.initializedKDs = new int[this.instance.nObjs()];
	this.currentKDs = new int[this.instance.nObjs()];
	for(int iObj = 0;iObj< instance.nObjs(); ++iObj){
	    this.setInitializedKD(iObj,-1);
	    this.UpdateCurrentK(iObj,0);
	}
	Log.comment(4, "}");
	}



    /**
     *My little method. It should add the hard constraints needed for
     *the sequential encoding, using afterKD as the superior limit for
     *the differential k. I need to initialize the variables before I
     *can use them to build clauses.
     *@param iObj The objective index
     *@param afterKD The desired max value for the objective 
     */
    
     public boolean UpdateCurrentK(int iObj , int afterKD ){
	 // Log.comment(6, "Internal:");
	 // Log.comment(6, "["+this.getCurrentKD(0));
	 // for(int iObj1 = 1; iObj1 < this.instance.nObjs(); ++iObj1)
	 //     Log.comment(6, ", "+this.getCurrentKD(iObj1));
	 // Log.comment(6, "]");
 	Log.comment(4, "{ SeqEncoder.UpdatecurrentK");
	if(this.getInitializedKD(iObj) == this.instance.getObj(iObj).getWeightDiff())
	    return false;
	if(this.getInitializedKD(iObj)< afterKD ){
	     // Y variables are also extended at 
	     this.extendInitializedIdsSInK(iObj, afterKD); 
	     this.extendInitializedIdsBInK(iObj, afterKD); 
	     this.setInitializedKD(iObj, afterKD);
	 }
	 if(afterKD == 0);
	 this.largerThan0(iObj, afterKD);
	 if(this.getCurrentKD(iObj) < afterKD){
	     this.blockingVariableB(iObj, afterKD);
	     if(iObj == 0)
		 Log.comment(5, "Clauses -1 4 8 9");
	     this.ifNotLessNotMore(iObj,afterKD);
	     Log.comment(5, "");
	     this.IfXAtLeastW(iObj, afterKD);
	     Log.comment(5, "");
	     this.IfLessAlsoMore(iObj, afterKD);
	     Log.comment(5, "");
	     this.IfLessAndIthXAtLeastIthW(iObj, afterKD);
	     this.setCurrentKD(iObj, afterKD);
	 }
	 Log.comment(5, "}");
	 return true;
     }



    /**
     * Initialize the container of the Blocking variables
     */
     
     private void initializeIdsB(){
	this.idsB = new int[this.instance.nObjs()][];
	for(int iObj = 0;iObj< instance.nObjs(); ++iObj){
	    Objective ithObj = instance.getObj(iObj);
	    this.idsB[iObj] = new int[ithObj.getWeightDiff() + 1];

       }
     }

    /**
     * Initialize the container of the Sequential variables
     */

    private void initializeIdsS(){
	this.idsS = new ArrayList<ArrayList<ArrayList<Integer>>>();
	for(int iObj = 0;iObj< instance.nObjs(); ++iObj){
 	    Objective ithObj = instance.getObj(iObj);
	    ArrayList<ArrayList<Integer>> arrayX_dK =new ArrayList<ArrayList<Integer>>(ithObj.getTotalLits());
	    int nLits = ithObj.getTotalLits();
	    for(int x = 1; x <= nLits ;++x){
		// + 1 necessary: remember kd is simultaneously a
		// value and an index
		ArrayList<Integer> arrayDK = new ArrayList<Integer>(ithObj.getWeightDiff() + 1);
		arrayX_dK.add(arrayDK);
	    }
		this.idsS.add(arrayX_dK);
	}
    }

    /**
     * Get the ithObj, x literal, kD S variable ID.
     */

    public int getS(int iObj, int x, int kD){
	int iX = x - 1;
	int iKD = kD ;
	List<ArrayList<Integer>> arrayX_dK = this.idsS.get(iObj);
	List<Integer> arrayDK = arrayX_dK.get(iX);
	int result = arrayDK.get(iKD);
	return result;
    }
 

    /**
     * Set the ithObj, ith Literal ith kD S variable to
     *@param iObj, the index of the objective
     *@param iX, the index of the literal from the objective
     *@param iKD, the index of the current differential k
     *@param id, the id of the variable S created
     * */

    public void setS(int iObj, int x, int kD, int id){
	int iX = x - 1;
	int iKD = kD;

	List<ArrayList<Integer>> arrayX_dK = this.idsS.get(iObj);
	List<Integer> arrayDK = arrayX_dK.get(iX);
	if( iKD > arrayDK.size() - 1){
	    Integer[] graftArray =  new Integer[iKD - arrayDK.size() + 1];
	    arrayDK.addAll(Arrays.asList(graftArray));
	}
	arrayDK.set(iKD, id);
	this.auxVariablesInverseIndex.put(id, new int[] {iObj,x,kD});
    }


    /**
     * Get the ithObj, ithX obj, kD Y variable ID.
     *@param iObj, the index of the objective
     *@param iX, the index of the literal from the objective
     *@param iKD, the index of the current differential k
     */

    public int getY(int iObj, int kD){
	int nLits = this.instance.getObj(iObj).getTotalLits();
	return 	this.getS(iObj, nLits, kD);
    }

    /**
     * Set the ithObj,  ith kD Y variable to
     *@param iObj, the index of the objective
     *@param iKD, the index of the current differential k
     *@param id, the id of the variable Y
     * */

    public void setY(int iObj, int iKD, int id){
	this.sTopVariablesInverseIndex.put(id, new int[] {iObj,iKD});
    }

    /**
     * Pops the ithObj,  ith kD Y variable 
     *@param iObj, the index of the objective
     *@param iKD, the index of the current differential k
     *@param id, the id of the variable Y
     * */

    public void popY(int iObj, int id){
	return ;
	// this.sTopVariablesInverseIndex.remove(id);
    }


    /**
     * Get the ithObj current upper limit on the initialized differential k
     *@param iObj, the index of the objective
     *@return  ithObj 
     */

    public int getInitializedKD(int iObj){
	return	 this.initializedKDs[iObj];
    }

    /**
     * Set the ithObj,  ith kD upper differential k
     *@param iObj, the index of the objective
     *@param iKD, the index of the current differential k
     * */

    public void setInitializedKD(int iObj, int kD){
	this.initializedKDs[iObj] = kD;
    }

    /**
     * Get the ithObj current upper limit on the initialized differential k
     *@param iObj, the index of the objective
     *@return  ithObj 
     */

    public int getCurrentKD(int iObj){
	return	 this.currentKDs[iObj];
    }

    /**
     * Set the ithObj,  ith kD upper differential k
     *@param iObj, the index of the objective
     *@param iKD, the index of the current differential k
     * */

    public void setCurrentKD(int iObj, int kD){
	this.currentKDs[iObj] = kD;
    }

    /**
     * Get the ithObj, ith differential k, Blocking  variable
     *@param iObj, the index of the objective
     *@param iKD, the index of the current differential k
     * */

     public int getB(int iObj, int dK){
	return this.idsB[iObj][dK];
    }

    /**
     * Set the ithObj, ith kD id of a blocking variable to id
     *@param iObj, the index of the objective
     *@param iKD, the index of the current differential k
     *@param the new blocking variable id
     **/
     public void setB(int iObj, int kD, int id){
	 this.idsB[iObj][kD] = id;
	}
    

    // private ConstrID AddRemovableClause(IVecInt setOfLiterals){

    // 	ConstrID constrainId = null;

    // 	for(int i = 0; i < setOfLiterals.size(); ++i)
    // 	    this.prettyPrintVariable(setOfLiterals.get(i));
    // 	try{
    // 	    constrainId = this.solver.addRemovableConstr(PBFactory.instance().mkClause(setOfLiterals));
    // 	} catch (ContradictionException e) {
    // 	    Log.comment(6, "contradiction when adding clause: ");
    // 	    for(int j = 0; j < setOfLiterals.size(); ++j)
    // 		Log.comment(6, " " + setOfLiterals.get(j) + " " );
    // 	    return constrainId;
    // 	}
    // 	return null;
    // }

     
    // /**
    //  * Add Clauses of type 4 in "On Using Incremental Encodings in..'"
    //  * deprecated, as of Wed 30 Oct 17:33:23 WET 2019
    //  */
    // private void ClausesIndependentOfK(){

    // 	for(int iObj = 0;iObj< this.instance.nObjs(); ++iObj){
    // 	    Objective ithObj = this.instance.getObj(iObj);
    // 	    this.extendInitializedIdsSInK(iObj, ithObj.getMaxAbsCoeff());
    // 	    int ithObjNLit = ithObj.getTotalLits();
    // 	    ReadOnlyVecInt ithObjLits = ithObj.getSubObjLits(0);
    // 	    ReadOnlyVec<Real> ithObjCoeffs = ithObj.getSubObjCoeffs(0);
    // 	    //	     assert ithObjNLit ==ithObjLits.size();
    // 	    for (int iX = 0 ; iX < ithObjNLit-1; ++iX){
    // 		int ithXW = Math.round(ithObjCoeffs.get(iX).asInt());
    // 		int sign = (ithXW > 0)? 1: -1;
    // 		ithXW = sign * ithXW;
    // 		int x =  sign * ithObjLits.get(iX);
    // 		for (int kD  = this.getCurrentKD(iObj) ;  kD <= ithXW; ++kD){
    // 		    int s = this.getS(iObj, iX, kD);
    // 		    IVecInt clauseSet = new VecInt(2);
    // 		    clauseSet.push(-x);
    // 		    clauseSet.push(s);
    // 		    this.AddClause(clauseSet);
    // 		}
    // 	    }
    // 	}
	 
    // }
     /**
      *The sum is always larger or equal than 0
      */
     private void largerThan0(int iObj, int afterKD){
	 Objective ithObj = this.instance.getObj(iObj);
	 int ithObjNLit = ithObj.getTotalLits();
	 //	     assert ithObjNLit ==ithObjLits.size();
	 for (int iX = 1 ; iX <= ithObjNLit; ++iX){
	     int s = this.getS(iObj, iX, 0);
	     this.AddClause(new VecInt(new int[]{s}));
	 }
 }



     /**
     * Clause 4. If a literal is contributing positively, then the sum
     * must be at least either the literal's weight or the upper
     * differential k, whichever is less.
     * @param iObj, the index of the objective function
     * @param afterKD, the new upper value of differential k for which
     * the semantics of the sequential encoding is complete
     */
    private void IfXAtLeastW(int iObj, int afterKD){

	    Objective ithObj = this.instance.getObj(iObj);
	    int ithObjNLit = ithObj.getTotalLits();
	    ReadOnlyVecInt ithObjLits = ithObj.getSubObjLits(0);
	    ReadOnlyVec<Real> ithObjCoeffs = ithObj.getSubObjCoeffs(0);
	    //	     assert ithObjNLit ==ithObjLits.size();
	    for (int iX = 1 ; iX <= ithObjNLit; ++iX){
		int ithXW = Math.round(ithObjCoeffs.get(iX-1).asInt());
		int sign = (ithXW > 0)? 1: -1;
		ithXW = sign * ithXW;
		int literal =  sign * ithObjLits.get(iX-1);
		int upperLimit = ithXW;
		upperLimit = (upperLimit < afterKD)? upperLimit: afterKD;
		for (int kD  = this.getCurrentKD(iObj) +1;
		     kD <= upperLimit; ++kD){
		    int s = this.getS(iObj, iX, kD);
		    IVecInt clauseSet = new VecInt(2);
		    clauseSet.push(-literal);
		    clauseSet.push(s);
		    this.AddClause(clauseSet);
		}
	    }
	}
	 

    /**
     * Clause -1. TODO NOT CORRECT Encoding that the value of the sum
     * is inductive on the index of the last literal of the sum
     * @param iObj, the index of the objective function
     * @param afterKD, the new upper value of differential k for which
     * the semantics of the sequential encoding is complete
     */
    private void ifNotLessNotMore(int iObj,int afterKD){

	int nLit = this.instance.getObj(iObj).getTotalLits();
	int x = nLit;
	    for (int kD  = this.currentKDs[iObj]+1;  kD <= afterKD; ++kD){
		int s1 = this.getS(iObj, x, kD-1);
		int s2 = this.getS(iObj, x, kD);
		IVecInt clauseSet = new VecInt(new int[] {s1,-s2});
		this.AddClause(clauseSet);
	    }
    }



    /**
     * Clause 8. Encoding that the value of the sum is inductive on
     * the index of the last literal of the sum
     * @param iObj, the index of the objective function
     * @param afterKD, the new upper value of differential k for which
     * the semantics of the sequential encoding is complete
     */
    private void IfLessAlsoMore(int iObj,int afterKD){

	int nLit = this.instance.getObj(iObj).getTotalLits();
	for (int x = 2 ; x <= nLit; ++x){
	    for (int kD  = this.currentKDs[iObj]+1;  kD <= afterKD; ++kD){
		IVecInt clauseSet = new VecInt(2);
		int s1 = this.getS(iObj, x-1, kD);
		int s2 = this.getS(iObj, x, kD);
		clauseSet.push(-s1);
		clauseSet.push(s2);
		this.AddClause(clauseSet);
	    }
	}
    }

    /**
     * Clause 9. If a given literal is present and the sum until that
     * same literal is at least k, then the sum including the literal
     * is at least the known estimative plus the literal's absolute
     * weight
     * @param iObj, the index of the objective function
     * @param afterKD, the new upper value of differential k for which
     * the semantics of the sequential encoding is complete
     */
    private void IfLessAndIthXAtLeastIthW(int iObj,int afterKD){

	int ithObjNLit = this.instance.getObj(iObj).getTotalLits();
	ReadOnlyVecInt ithObjLits = this.instance.getObj(iObj).getSubObjLits(0);
	ReadOnlyVec<Real> ithObjCoeffs = this.instance.getObj(iObj).getSubObjCoeffs(0);
	assert ithObjNLit ==ithObjLits.size();

	for (int iX = 2 ; iX <= ithObjNLit; ++iX){
	    int ithXW = ithObjCoeffs.get(iX-1).asInt();
	    int sign = (ithXW > 0)? 1: -1;
	    ithXW = sign * ithXW;
	    int literal = sign * ithObjLits.get(iX-1);
	    int lowerLimit = this.currentKDs[iObj]- ithXW + 1;
	    lowerLimit = (lowerLimit > 1)? lowerLimit: 1;
	    int upperLimit = afterKD - ithXW ;
	    // upperLimit = (upperLimit >= 1)? upperLimit: 1;
	    for(int kD  = lowerLimit;  kD <= upperLimit ; ++kD){
		int s1 = this.getS(iObj, iX-1, kD);
		int s2 = this.getS(iObj, iX, kD + ithXW );
		IVecInt clauseSet = new VecInt(new int[] {-s1,-literal,s2});
		this.AddClause(clauseSet);
	    }
	}
    }

     /**
      *blocks the older non incremental clauses
      */
     private void blockingVariableB(int iObj,int afterK){
	 int beforeK = this.getCurrentKD(iObj);
	 if(beforeK >= 1){
	     IVecInt clauseSet = new VecInt(1);
	     clauseSet.push(this.getB(iObj, beforeK));
	     this.AddClause(clauseSet);
	 }
     }


     /**
      *Just a wrapper to PBSolver newVar, that returns the new id
      */

    private int newVar(){
	this.solver.newVars(1);
	return solver.nVars();
    }




    /** 
     *  extend the B variables in k
     */

    private void extendInitializedIdsBInK(int iObj, int afterKD){

	for(int kD = this.getInitializedKD(iObj)+1; kD <= afterKD; ++kD){
	    this.setB(iObj,kD , this.newVar());	     
	}
	 
    }

    /** 
     * extend the S variables in the differential k index
     */

    public void extendInitializedIdsSInK(int iObj, int afterKD){

	int nLit = this.instance.getObj(iObj).getTotalLits();
	for(int kd = this.initializedKDs[iObj]+1; kd <= afterKD ; ++kd){
	    for (int x = 1 ; x <= nLit; ++x){
		/* Log.comment(6, iObj + " " + iX + " " + kd + " " + this.newVar()); */
		this.setS(iObj, x, kd, this.newVar());
	    }
	}
	this.extendInitializedIdsYinK(iObj, afterKD);
    }


    /** 
     * extend the Y variables in the differential k index until afterKD. Assumes
     * the S variables are already extended accordingly
     * @param iObj
     * @param afterKD
     */

    private void extendInitializedIdsYinK(int iObj, int afterKD){
	int nLit = this.instance.getObj(iObj).getTotalLits();
	int beforeKD = this.getCurrentKD(iObj);
	for(int kd = this.initializedKDs[iObj]+1; kd <= afterKD  ; ++kd){
	    if( beforeKD != 0)
		this.popY(iObj, this.getS(iObj, nLit, beforeKD));
	    this.setY(iObj, kd, this.getS(iObj, nLit, afterKD ));
	    
	}
    }



    /**
     * Get the objective from an Y variable
     * @param literal
     */

    public int getIObjFromY(int literal){
	assert this.isY(literal);
	literal = (literal>0)? literal: -literal;
	return this.sTopVariablesInverseIndex.get(literal)[0] ;
    }

    /**
     * return the value of the differential k from an S variable
     * @param literal
     */

    public int getKDFromY(int literal){
	assert this.isY(literal);
	literal = (literal>0)? literal: -literal;
	return this.sTopVariablesInverseIndex.get(literal)[1] ;
    }

    /**
     * Get the objective from an Y variable
     * @param literal
     */

    public int getObjFromSVariable(int literal){
	assert this.isS(literal);
	literal = (literal>0)? literal: -literal;
	return this.auxVariablesInverseIndex.get(literal)[0] ;
    }

    /**
     * return the value of the differential k from an S variable
     * @param literal
     */

    public int getKDFromSVariable(int literal){
	assert this.isS(literal);
	literal = (literal>0)? literal: -literal;
	return this.auxVariablesInverseIndex.get(literal)[2] ;

    }
    /**
     * Get the objective from an B variable
     * @param literal
     */

    public int getObjFromBVariable(int literal){
	assert this.isB(literal);
	literal = (literal>0)? literal: -literal;
	return this.bVariablesInverseIndex.get(literal) ;
    }

    /**
     * return the value of the differential k from an B variable
     * @param literal
     */

    /**
     * return the value of the literal index
     * @param literal
     */

    public int getXFromSVariable(int literal){
	assert this.isS(literal);
	literal = (literal>0)? literal: -literal;
	return this.auxVariablesInverseIndex.get(literal)[1] ;
    }


    /**
     *Checks if literal is an Y variable
     *@param literal
     */

    public boolean isY(int literal){
	literal = (literal>0)? literal: -literal;
	if(this.sTopVariablesInverseIndex.containsKey(literal))
	    return true;
	return false;
    }

    /**
     *Checks if literal is an S variable
     *@param literal
     */

    public boolean isS(int literal){
	literal = (literal>0)? literal: -literal;
	if(this.auxVariablesInverseIndex.containsKey(literal))
	    return true;
	return false;
    }
    /**
     *Checks if literal is an B variable
     *@param literal
     */

    public boolean isB(int literal){
	literal = (literal>0)? literal: -literal;
	if(this.bVariablesInverseIndex.containsKey(literal))
	    return true;
	return false;
    }

     public String prettyFormatVecInt(IVecInt vecInt){
	 String result = "";
	 for(int j = 0; j < vecInt.size(); ++j)
	     result += this.prettyFormatVariable(vecInt.get(j));
	 return result;
     }

     public void prettyPrintVecInt(IVecInt vecInt){
	Log.comment(6,prettyFormatVecInt(vecInt));
	 return;
     }

    public void prettyPrintVariable(int literal){
	Log.comment(6,prettyFormatVariable(literal));
    }

    public String prettyFormatVariable(int literal){
	int sign =(literal>0)? 1: -1;
	int id =  literal * sign;

	if(this.isY(id)){
	    int iObj = this.getIObjFromY(id);
	    int kd = this.getKDFromY(id);
	    return "Y[" + iObj + ", " + kd +"]"+ "::" + literal + " ";
	}
	 
	if(this.isS(id)){
	    int iObj = this.getObjFromSVariable(id);
	    int iX = this.getXFromSVariable(id);
	    int kd = this.getKDFromSVariable(id);
	    return  "S[" + iObj + ", " + iX + ", " + kd +"]"+"::" + literal+ " ";

	}
	if(this.isB(id)){
	    int iObj = this.getObjFromBVariable(id);
	    return "B[" + ", " + iObj +"]"+ "::" + id + " ";
	}
	if(id < this.firstVariable){
	    return (sign>0? "+":"-")+"X["+id+"] ";
	}
	return "";
    }


    /**
     * Generate the upper limit assumptions
     */
    public IVecInt generateUpperBoundAssumptions(int[] UpperKD){

	IVecInt assumptions = new VecInt(new int[]{});
	for(int iObj = 0; iObj < this.instance.nObjs(); ++iObj){
	    Objective ithObjective = this.instance.getObj(iObj);
	    if(UpperKD[iObj]  < ithObjective.getWeightDiff())
		assumptions.push(-this.getY(iObj, UpperKD[iObj] + 1));
	    
	    ReadOnlyVecInt objectiveLits = ithObjective.getSubObjLits(0);
	    ReadOnlyVec<Real> objectiveCoeffs = ithObjective.getSubObjCoeffs(0);
	    int sign = 1;
	    int ithAbsoluteWeight;
	    for(int iX = 0, nX = ithObjective.getTotalLits(); iX <nX; iX ++){
		ithAbsoluteWeight = objectiveCoeffs.get(iX).asInt();
		sign = (ithAbsoluteWeight > 0? 1 : -1);
		ithAbsoluteWeight *= sign;
		if( ithAbsoluteWeight > UpperKD[iObj])
		    assumptions.push(- sign * objectiveLits.get(iX));
	    }

	}

	return assumptions;
    }

}

