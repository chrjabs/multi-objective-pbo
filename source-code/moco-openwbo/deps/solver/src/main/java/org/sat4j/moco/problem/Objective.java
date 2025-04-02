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
package org.sat4j.moco.problem;

import org.sat4j.core.ReadOnlyVec;
import org.sat4j.core.ReadOnlyVecInt;
import org.sat4j.core.Vec;
import org.sat4j.moco.pb.PBExpr;
import org.sat4j.moco.pb.PBSolver;
import org.sat4j.moco.util.Real;
import org.sat4j.specs.IVec;

/**
 * Superclass for representations of PB objective function.
 * @author Miguel Terra-Neves
 */
public abstract class Objective {
    
    /**
     * Stores the sub-objectives associated with the objective (e.g. reduced division objectives).
     */
    private IVec<PBExpr> sub_objs = new Vec<PBExpr>();
    
    /**
     * Creates an instance of a PB objective function.
     * @param e The objective's PB expression.
     */
    protected Objective(PBExpr e) { this.sub_objs.push(e); }
    
    /**
     * Creates an instance of a reduced PB objective function.
     * @param es The objective's sub-expressions.
     */
    protected Objective(IVec<PBExpr> es) { es.copyTo(this.sub_objs); }
    
    /**
     * Retrieves the number of sub-objectives in the PB objective function.
     * @return The objective's number of sub-objectives.
     */
    public int nSubObj() { return this.sub_objs.size(); }
    
    /**
     * Retrieves an expression that represents one of the PB objective's sub-objectives.
     * @param i The sub-objective index.
     * @return The {@code i}-th sub-objective's expression.
     */
    public PBExpr getSubObj(int i) {
        assert(i < nSubObj());
        return this.sub_objs.get(i);
    }
    
    /**
     * Retrieves the literals of one of the PB objective's sub-objectives.
     * @param i The sub-objective index.
     * @return The literals in the {@code i}-th sub-objective.
     */
    public ReadOnlyVecInt getSubObjLits(int i) { return getSubObj(i).getLits(); }
    
    /**
     *  returns the total number of literals for each objective
     */

    public int getTotalLits(){
	int total = 0;
	for (int i = 0; i < this.nSubObj() ; ++i)
	    total += this.getSubObjLits(i).size();
	return total;
    }

    /** Added by Joao Cortes at Wed 23 Oct 16:22:19 WEST 2019
     * get the Max Weight of an objective.
     * This only works well for integer weights
     */
    public int getMaxValue(){

	int total = 0;
	for (int i = 0; i < this.nSubObj() ; ++i){
	    try{	    
	    total+=this.getSubObj(i).getMaxSum().asIntExact();
	    }catch(ArithmeticException e){
		if(total > 0)
		    return Integer.MAX_VALUE;
		else
		    return -Integer.MAX_VALUE;
}
	}
	return total;
    }


    /**
     * get the Min Weight of an objective
     */
    public int getMinValue(){

	int total = 0;
	for (int i = 0; i < this.nSubObj() ; ++i){
	    try{
		total += this.getSubObj(i).getMinSum().asIntExact();
	    }catch(ArithmeticException e){
		if(total > 0)
		    return Integer.MAX_VALUE;
		else
		    return -Integer.MAX_VALUE;
		    
	    }
	}
	return total;
    }

    /**
     * get the range of the objective
     */
    public int getWeightDiff(){
	int result = this.getMaxValue() - this.getMinValue();
	if(result > 0)
	    return result;
	return -result;
    }

    /**
     * get the greatest of all coeffs
     */
    public int getMaxAbsCoeff(){

	int max = 0;
	for (int i = 0; i < this.nSubObj() ; ++i){
	    ReadOnlyVec<Real> ithCoeffs =  this.getSubObjCoeffs(i);
	    int ithCoeffsN = ithCoeffs.size();
	    for(int  k = 0; k < ithCoeffsN;++k){
		int ithKthCoeff = ithCoeffs.get(k).asInt();
		ithKthCoeff = (ithKthCoeff > 0)? ithKthCoeff: -ithKthCoeff;
		if( ithKthCoeff > max)
		    max = ithKthCoeff;
	    }
	}
	return max;
    }




    public ReadOnlyVec<Real> getSubObjCoeffs(int i) { return getSubObj(i).getCoeffs(); }
    
    /**
     * Computes the value of the PB objective under a given assignment.
     * @param a The assignment. The {@code i}-th position is the Boolean value assigned to variable
     * {@code i+1}.
     * @return The objective's value under assignment {@code a}.
     */
    public abstract Real evaluate(boolean[] a);
    
    /**
     * Computes the value of the PB objective under a model in a given PB solver.
     * @param s The solver.
     * @return The objective's value under the model stored in {@code s}.
     */
    public abstract Real evaluate(PBSolver s);
    
    public int evaluateDiff(PBSolver solver){
	return this.evaluate(solver).asIntExact() - this.getMinValue();
}

    public int evaluateDiff(boolean[] model){
	return this.evaluate(model).asIntExact() - this.getMinValue();
}
}
