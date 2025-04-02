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
package org.sat4j.moco.pb;

import java.util.Map;
import java.util.HashMap;

import org.sat4j.core.ReadOnlyVec;
import org.sat4j.core.ReadOnlyVecInt;
import org.sat4j.moco.util.Real;
import org.sat4j.specs.IVec;
import org.sat4j.specs.IVecInt;

/**
 * Representation of a PB expression.
 * @author Miguel Terra-Neves
 */
public class PBExpr {
    
    /**
     * Stores the PB expression's literals.
     */
    private ReadOnlyVecInt lits = null;
    
    /**
     * Map ids to weights.
     */
    private Map<Integer, Real> litsToWeight = null;
    

    /**
     * Stores the PB expression's coefficients.
     */
    private ReadOnlyVec<Real> coeffs = null;
    
    /**
     * Creates an instance of a PB expression.
     * @param lits The expression's literals.
     * @param coeffs The expression's coefficients.
     */
    public PBExpr(IVecInt lits, IVec<Real> coeffs) {
        this.lits = new ReadOnlyVecInt(lits);
        this.coeffs = new ReadOnlyVec<Real>(coeffs);
	this.litsToWeight = new HashMap<Integer, Real>();
        assert(this.lits.size() == this.coeffs.size());
	for(int iLit = 0, n = lits.size();  iLit < n ;iLit++){
	    this.litsToWeight.put(this.lits.get(iLit), this.coeffs.get(iLit));
	}
    }
    
    /**
     * Retrieves the PB expression's literals.
     * @return The expression's literals.
     */
    public ReadOnlyVecInt getLits() { return this.lits; }
    
    /**
     * Retrieves the PB expression's coefficients.
     * @return The expression's coefficients.
     */
    public ReadOnlyVec<Real> getCoeffs() { return this.coeffs; }
    
    /**
     * Retrieves the number of terms in the PB expression.
     * @return The expression's number of terms.
     */
    public int nTerms() { return getLits().size(); }
    
    /**
     * Checks if a given literal is true under a given assignment.
     * @param a The assignment. The {@code i}-th position is the Boolean value assigned to variable
     * {@code i+1}.
     * @param lit The literal.
     * @return True if {@code lit} is true under assignment {@code a}, false otherwise.
     */
    private boolean isTrue(boolean[] a, int lit) {
        int var = Math.abs(lit)-1;
        return lit > 0 ? a[var] : !(a[var]);
    }
    
    /**
     * Computes the value of the PB expression under a given assignment.
     * @param a The assignment. The {@code i}-th position is the Boolean value assigned to variable
     * {@code i+1}.
     * @return The expression's value under assignment {@code a}.
     */
    public Real evaluate(boolean[] a) {
        Real val = Real.ZERO;
        for (int i = 0; i < nTerms(); ++i) {
            if (isTrue(a, getLits().get(i))) {
                val = val.add(getCoeffs().get(i));
            }
        }
        return val;
    }
    
    /**
     * Computes the value of the PB expression under a model in a given PB solver.
     * @param s The solver.
     * @return The expression's value under the model stored in {@code s}.
     */
    public Real evaluate(PBSolver s) {
        assert(s.isSolved() && s.isSat());
        Real val = Real.ZERO;
        for (int i = 0; i < nTerms(); ++i) {
            if (s.modelValue(getLits().get(i))) {
                val = val.add(getCoeffs().get(i));
            }
        }
        return val;
    }


    /**
     * Get the weight given the literal id
     */
    public Real weightFromLit(int id){
	Real weight = this.litsToWeight.get(id);
	return weight;
}

    
    /**
     * Auxiliary interface for strategy objects used by {@link PBExpr#checkAllTerms(ITermChecker)} to check
     * properties of the PB expression's terms (e.g. if all coefficients are positive).
     * @author Miguel Terra-Neves
     */
    private interface ITermChecker {
        
        /**
         * Checks some property of a given term.
         * @param l The term's literal.
         * @param c The term's coefficient.
         * @return True if the property is verified for the given term, false otherwise.
         */
        public boolean check(int l, Real c);
        
    }
    
    /**
     * Checks if a given property is verified by all terms.
     * @param c An {@link ITermChecker} object that implements the property check.
     * @return True if the check implemented by {@code c} passes for all terms, false otherwise.
     */
    private boolean checkAllTerms(ITermChecker c) {
        for (int i = 0; i < nTerms(); ++i) {
            if (!c.check(getLits().get(i), getCoeffs().get(i))) { return false; }
        }
        return true;
    }
    
    /**
     * Checks if all coefficients in the PB expression are 1.
     * @return True if all of the PB expression's coefficients are 1, false otherwise.
     */

    boolean allUnitCoeffs() {
        return checkAllTerms(new ITermChecker() {
            public boolean check(int l, Real c) { return c.equals(Real.ONE); }
        });
    }
    
    /**
     * Checks if all coefficients in the PB expression are greater than or equal to a given value.
     * @param d The value.
     * @return True if all coefficients are greater than or equal to {@code d}, false otherwise.
     */
    boolean allCoeffsGE(final Real d) {
        return checkAllTerms(new ITermChecker() {
            public boolean check(int l, Real c) { return c.greaterOrEqual(d); }
        });
    }
    
    /**
     * Checks if all coefficients in the PB expression are less than or equal to a given value.
     * @param d The value.
     * @return True if all coefficients are less than or equal to {@code d}, false otherwise.
     */
    boolean allCoeffsLE(final Real d) {
        return checkAllTerms(new ITermChecker() {
            public boolean check(int l, Real c) { return c.lessOrEqual(d); }
        });
    }
    
    /**
     * {@link ITermChecker} implementation that checks for terms with positive coefficients.
     */
    private final static ITermChecker POS_CHECKER = new ITermChecker() {
        public boolean check(int l, Real c) { return c.isPositive(); }
    };
    
    /**
     * {@link ITermChecker} implementation that checks for terms with negative coefficients.
     */
    private final static ITermChecker NEG_CHECKER = new ITermChecker() {
        public boolean check(int l, Real c) { return c.isNegative(); }
    };
    
    /**
     * Checks if all coefficients in the PB expression are positive.
     * @return True if all coefficients are positive, false otherwise.
     */
    boolean allPosCoeffs() { return checkAllTerms(POS_CHECKER); }
    
    /**
     * Checks if all coefficients in the PB expression are negative.
     * @return True if all coefficients are negative, false otherwise.
     */
    boolean allNegCoeffs() { return checkAllTerms(NEG_CHECKER); }
    
    /**
     * Computes the sum of the coefficients of the terms that satisfy a condition implemented by some
     * {@link ITermChecker} object.
     * @param checker The term checker.
     * @return The coefficient sum.
     */
    private Real getCoeffSum(ITermChecker checker) {
        Real sum = Real.ZERO;
        for (int i = 0; i < nTerms(); ++i) {
            Real coeff = getCoeffs().get(i);
            if (checker.check(getLits().get(i), coeff)) {
                sum = sum.add(coeff);
            }
        }
        return sum;
    }
    
    /**
     * Computes the maximum value possible for this expression.
     * @return The expression's maximum value.
     */
    public Real getMaxSum() { return getCoeffSum(POS_CHECKER); }
    
    /**
     * Computes the minimum value possible for this expression.
     * @return The expression's minimum value.
     */
    public Real getMinSum() { return getCoeffSum(NEG_CHECKER); }
    
}
