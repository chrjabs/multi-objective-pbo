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
import org.sat4j.specs.IVecInt;

/**
 * Representation of a sum of divisions PB objective.
 * @author Miguel Terra-Neves
 */
public class DivObj extends Objective {

    /**
     * The PB objective's numerators.
     */
    private IVec<PBExpr> nums = null;
    
    /**
     * The PB objective's denominators.
     */
    private IVec<PBExpr> dens = null;

    /**
     * Creates the sub-objectives that correspond to a sum of divisions PB objective's reduction.
     * @param nums The objective's numerators.
     * @param dens The objective's denominators.
     * @return A vector with the PB expressions that correspond to the reduced objective's sub-objectives.
     */
    private static IVec<PBExpr> mkSubObjs(IVec<PBExpr> nums, IVec<PBExpr> dens) {
        assert(nums.size() == dens.size());
        IVec<PBExpr> objs = new Vec<PBExpr>(nums.size() + dens.size());
        for (int i = 0; i < nums.size(); ++i) {
            objs.unsafePush(nums.get(i));
            PBExpr den = dens.get(i);
            IVec<Real> neg_den_coeffs = new Vec<Real>(den.nTerms());
            for (int j = 0; j < den.nTerms(); ++j) {
                neg_den_coeffs.unsafePush(den.getCoeffs().get(j).negate());
            }
            objs.unsafePush(new PBExpr(den.getLits(), neg_den_coeffs));
        }
        return objs;
    }
    
    /**
     * Creates a set of PB expressions from given sets of literals and coefficients.
     * @param lits The expressions' literal vectors.
     * @param coeffs The expressions' coefficient vectors.
     * @return A vector of expressions, where the {@code i}-th expression contains the {@code i}-th literals
     * of {@code lits} and the {@code i}-th coefficients of {@code coeffs}.
     */
    private static IVec<PBExpr> mkExprs(IVec<IVecInt> lits, IVec<IVec<Real>> coeffs) {
        assert(lits.size() == coeffs.size());
        IVec<PBExpr> es = new Vec<PBExpr>(lits.size());
        for (int i = 0; i < lits.size(); ++i) {
            es.unsafePush(new PBExpr(lits.get(i), coeffs.get(i)));
        }
        return es;
    }
    
    /**
     * Creates an instance of a sum of divisions PB objective.
     * @param nums The objective's numerators.
     * @param dens The objective's denominators.
     */
    public DivObj(IVec<PBExpr> nums, IVec<PBExpr> dens) {
        super(mkSubObjs(nums, dens));
        this.nums = nums;
        this.dens = dens;
    }
    
    /**
     * Creates an instance of a sum of divisions PB objective.
     * @param num_lits The numerators' literal vectors.
     * @param num_coeffs The numerators' coefficient vectors.
     * @param den_lits The denominators' literal vectors.
     * @param den_coeffs The denominators' coefficient vectors.
     */
    public DivObj(IVec<IVecInt> num_lits, IVec<IVec<Real>> num_coeffs,
                  IVec<IVecInt> den_lits, IVec<IVec<Real>> den_coeffs) {
        this(mkExprs(num_lits, num_coeffs), mkExprs(den_lits, den_coeffs));
    }
    
    /**
     * Retrieves the number of division terms in the sum of divisions PB objective.
     * @return The objective's number of division terms.
     */
    public int nDivs() { return this.nums.size(); }
    
    /**
     * Retrieves the PB expression that corresponds to a given numerator of the sum of divisions PB objective.
     * @param i The division term index.
     * @return The numerator of the {@code i}-th division term.
     */
    public PBExpr getNum(int i) { return this.nums.get(i); }
    
    /**
     * Retrieves the PB expression that corresponds to a given denominator of the sum of divisions PB
     * objective.
     * @param i The division term index.
     * @return The denominator of the {@code i}-th division term.
     */
    public PBExpr getDen(int i) { return this.dens.get(i); }
    
    /**
     * Retrieves the literals in a given numerator of the sum of divisions PB objective.
     * @param i The division term index.
     * @return The literals in the numerator of the {@code i}-th division term.
     */
    public ReadOnlyVecInt getNumLits(int i) { return getNum(i).getLits(); }
    
    /**
     * Retrieves the literals in a given denominator of the sum of divisions PB objective.
     * @param i The division term index.
     * @return The literals in the denominator of the {@code i}-th division term.
     */
    public ReadOnlyVecInt getDenLits(int i) { return getDen(i).getLits(); }
    
    /**
     * Retrieves the coefficients in a given numerator of the sum of divisions PB objective.
     * @param i The division term index.
     * @return The coefficients in the numerator of the {@code i}-th division term.
     */
    public ReadOnlyVec<Real> getNumCoeffs(int i) { return getNum(i).getCoeffs(); }
    
    /**
     * Retrieves the coefficients in a given denominator of the sum of divisions PB objective.
     * @param i The division term index.
     * @return The coefficients in the denominator of the {@code i}-th division term.
     */
    public ReadOnlyVec<Real> getDenCoeffs(int i) { return getDen(i).getCoeffs(); }

    @Override
    public Real evaluate(boolean[] a) {
        Real val = Real.ZERO;
        for (int i = 0; i < nDivs(); ++i) {
            Real num_val = getNum(i).evaluate(a);
            val = num_val.equals(Real.ZERO) ? val : val.add(num_val.divide(getDen(i).evaluate(a)));
        }
        return val;
    }

    // FIXME: too similar to evaluate(boolean[] a); refactor
    @Override
    public Real evaluate(PBSolver s) {
        Real val = Real.ZERO;
        for (int i = 0; i < nDivs(); ++i) {
            Real num_val = getNum(i).evaluate(s);
            val = num_val.equals(Real.ZERO) ? val : val.add(num_val.divide(getDen(i).evaluate(s)));
        }
        return val;
    }
    
}
