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

import java.math.BigInteger;

import org.sat4j.core.ReadOnlyVec;
import org.sat4j.core.ReadOnlyVecInt;
import org.sat4j.core.Vec;
import org.sat4j.core.VecInt;
import org.sat4j.moco.util.Real;
import org.sat4j.pb.IPBSolver;
import org.sat4j.specs.ContradictionException;
import org.sat4j.specs.IConstr;
import org.sat4j.specs.IVec;
import org.sat4j.specs.IVecInt;

/**
 * Superclass for representations of PB constraints.
 * @author Miguel Terra-Neves
 */
public abstract class PBConstr {
    
    /**
     * The left-hand side of the PB constraint.
     */
    private PBExpr lhs = null;
    
    /**
     * The right-hand side of the PB constraint.
     */
    private Real rhs = null;
    
    /**
     * Creates an instance of a PB constraint.
     * @param lits The literals in the constraint's left-hand side.
     * @param coeffs The coefficients in the constraint's left-hand side.
     * @param rhs The constraint's right-hand side.
     */
    public PBConstr(IVecInt lits, IVec<Real> coeffs, Real rhs) {
        this.lhs = new PBExpr(lits, coeffs);
        this.rhs = rhs;
    }
    
    /**
     * Retrieves the PB constraint's left-hand side.
     * @return The constraint's left-hand side.
     */
    public PBExpr getLHS() { return this.lhs; }
    
    /**
     * Retrieves the PB constraint's right-hand side.
     * @return The constraint's right-hand side.
     */
    public Real getRHS() { return this.rhs; }
    
    /**
     * Retrieves the literals in the PB constraint's left-hand side.
     * @return The literals in the constraint's left-hand side.
     */
    public ReadOnlyVecInt getLits() { return this.lhs.getLits(); }
    
    /**
     * Retrieves the coefficients in the PB constraint's right-hand side.
     * @return The coefficients in the constraint's right-hand side.
     */
    public ReadOnlyVec<Real> getCoeffs() { return this.lhs.getCoeffs(); }
    
    /**
     * Sets a given literal as the activator for this PB constraint.
     * @param act The activator literal.
     */
    public void setActivator(int act) {
        IVecInt lits = new VecInt(getLits().size() + 1);
        IVec<Real> coeffs = new Vec<Real>(getCoeffs().size() + 1);
        getLits().copyTo(lits);
        getCoeffs().copyTo(coeffs);
        lits.push(-act);
        coeffs.push(getActivatorCoeff());
        setLHS(lits, coeffs);
    }
    
    /**
     * {@link #setActivator(int)} is a template method that uses this method to retrieve the coefficient
     * that should be associated with the activator literal.
     * Must be implemented by sub-classes.
     * @return The coefficient for the activator literal.
     */
    protected abstract Real getActivatorCoeff();
    
    /**
     * Sets the left-hand side of this PB constraint to a given expression.
     * @param expr The new left-hand side.
     */
    protected void setLHS(PBExpr expr) { this.lhs = expr; }
    
    /**
     * Sets the left-hand side of this PB constraint to an expression with the given literals and
     * coefficients.
     * @param lits The literals for the new left-hand side.
     * @param coeffs The coefficients for the new left-hand side.
     */
    protected void setLHS(IVecInt lits, IVec<Real> coeffs) { setLHS(new PBExpr(lits, coeffs)); }

    /**
     * Template method that scales the PB constraint's coefficients and right-hand side to integer values
     * before adding the constraint to the PB solver object {@code solver}.
     * {@link #addScaledToSolver(IPBSolver, IVecInt, IVec, BigInteger)}, which must be implemented by
     * sub-classes, is called with the scaled result.
     * @param solver The PB solver.
     * @return A constraint ID object that can be used later for removal
     * ({@link PBSolver#removeConstr(ConstrID)} and {@link PBSolver#removeConstrs(IVec)}).
     * @throws ContradictionException if {@code solver} detects that the formula becomes unsatisfiable after
     * adding this constraint.
     */
    protected IConstr scaleToIntThenAddToSolver(IPBSolver solver) throws ContradictionException {
        int factor = getRHS().nDecimals();
        for (int i = 0; i < getCoeffs().size(); ++i) {
            factor = Math.max(factor, getCoeffs().get(i).nDecimals());
        }
        IVec<BigInteger> coeffs = new Vec<BigInteger>();
        for (int i = 0; i < getCoeffs().size(); ++i) {
            coeffs.push(getCoeffs().get(i).scaleByPowerOfTen(factor).asBigIntegerExact());
        }
        BigInteger rhs = getRHS().scaleByPowerOfTen(factor).asBigIntegerExact();
        return addScaledToSolver(solver, getLits(), coeffs, rhs);
    }
    
    /**
     * Checks if the PB constraint becomes violated if the left-hand side equals a given value.
     * @param val The value to be tested.
     * @return True if a left-hand side equal to {@code val} violates the constraint, false otherwise.
     */
    public abstract boolean violatedBy(Real val);
    
    /**
     * Adds the PB constraint to a given PB solver.
     * @param solver The PB solver.
     * @return A constraint ID object that can be used later for removal
     * ({@link PBSolver#removeConstr(ConstrID)} and {@link PBSolver#removeConstrs(IVec)}).
     * @throws ContradictionException if {@code solver} detects that the formula becomes unsatisfiable after
     * adding this constraint.
     */
    abstract IConstr addToSolver(IPBSolver solver) throws ContradictionException;
    
    /**
     * Method called by {@link #scaleToIntThenAddToSolver(IPBSolver)} after scaling the PB constraint's
     * coefficients and right-hand side to integer values.
     * Adds the scaled constraint to a given PB solver object.
     * Must be implemented by sub-classes.
     * @param solver The PB solver.
     * @param lits The literal in the constraint's left-hand side.
     * @param coeffs The coefficients in the constraint's left-hand side scaled to integer values.
     * @param rhs The constraint's right-hand side scaled to an integer value.
     * @return A constraint ID object that can be used later for removal
     * ({@link PBSolver#removeConstr(ConstrID)} and {@link PBSolver#removeConstrs(IVec)}).
     * @throws ContradictionException if {@code solver} detects that the formula becomes unsatisfiable after
     * adding this constraint.
     */
    protected abstract IConstr addScaledToSolver(
            IPBSolver solver, IVecInt lits, IVec<BigInteger> coeffs, BigInteger rhs) throws ContradictionException;
    
    /**
     * Retrieves a string representation of the PB constraint's operator (e.g. ">=" for greater-or-equal
     * constraints).
     * @return The string representation of the constraint's operator.
     */
    public abstract String getOpStrRep();
    
}
