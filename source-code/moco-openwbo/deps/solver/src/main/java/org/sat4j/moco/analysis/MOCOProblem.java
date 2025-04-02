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
package org.sat4j.moco.analysis;

import org.moeaframework.core.Solution;
import org.moeaframework.core.variable.EncodingUtils;
import org.moeaframework.problem.AbstractProblem;
import org.sat4j.core.Vec;
import org.sat4j.moco.pb.PBConstr;
import org.sat4j.moco.problem.Instance;
import org.sat4j.moco.problem.Objective;
import org.sat4j.moco.util.Real;
import org.sat4j.specs.IVecInt;

/**
 * Implementation of the MOCO problem as an {@link AbstractProblem} in the MOEA framework.
 * To be used for performance analysis.
 * @author Miguel Terra-Neves
 */
class MOCOProblem extends AbstractProblem {

    /**
     * Stores the actual MOCO instance.
     */
    private Instance instance = null;
    
    /**
     * Creates an instance of a MOEA framework MOCO problem with a given MOCO instance.
     * @param instance The instance.
     */
    public MOCOProblem(Instance instance) {
        super(instance.nVars(), instance.nObjs(), 1);
        this.instance = instance;
    }
    
    /**
     * Retrieves an assignment from a given MOEA framework {@link Solution} object.
     * @param sol The solution.
     * @return The assignment.
     */
    public boolean[] getAssignment(Solution sol) {
        assert(sol.getNumberOfVariables() == getNumberOfVariables());
        boolean[] a = new boolean[getNumberOfVariables()];
        for (int j = 0; j < sol.getNumberOfVariables(); ++j) {
            a[j] = EncodingUtils.getBoolean(sol.getVariable(j));
        }
        return a;
    }

    /**
     *returns the maximal point
     */

    public Vec<Integer> maxPoint(){
	Vec<Integer> result = new Vec<Integer>();
	for(int i = 0, n = this.instance.nObjs();i<n;i++)
	    result.push(this.instance.getObj(i).getMaxValue());
	return result;			
    }
    /**
     *returns the minimal point
     */

    public Vec<Integer> minPoint(){
	Vec<Integer> result = new Vec<Integer>();
	for(int i = 0, n = this.instance.nObjs();i<n;i++)
	    result.push(this.instance.getObj(i).getMinValue());
	return result;			
    }

    /**
     * Computes the constraint violation and objective function cost values of a given MOEA framework
     * {@link Solution} object.
     * @param sol The solution.
     */
    public void evaluate(Solution sol) {
        boolean[] a = getAssignment(sol);
        for (int i = 0; i < this.instance.nObjs(); ++i) {
            Objective obj = this.instance.getObj(i);
            sol.setObjective(i, obj.evaluate(a).asDouble());
        }
        double viol = 0.0;
        for (int i = 0; i < this.instance.nConstrs(); ++i) {
            PBConstr constr = this.instance.getConstr(i);
            Real lhs = constr.getLHS().evaluate(a);
            viol += constr.violatedBy(lhs) ? constr.getRHS().subtract(lhs).abs().asDouble() : 0.0;
        }
        sol.setConstraint(0, viol);
    }

    /**
     * Creates an empty MOEA framework {@link Solution} object.
     * @return The solution object.
     */
    public Solution newSolution() {
        Solution sol = new Solution(getNumberOfVariables(), getNumberOfObjectives(), getNumberOfConstraints());
        for (int i = 0; i < getNumberOfVariables(); ++i) {
            sol.setVariable(i, EncodingUtils.newBoolean());
        }
        return sol;
    }
    
    /**
       remove objective
*/
    public void removeObj(int i){this.instance.removeObj(i);}

    /**
     * Creates an empty MOEA framework {@link Solution} object purely for storing the objective function
     * cost values of some assignment in a memory efficient manner.
     * @return The solution object.
     */
    Solution newCostVec() { return new Solution(0, getNumberOfObjectives(), 0); }
    
}
