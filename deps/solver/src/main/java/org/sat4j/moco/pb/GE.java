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

import org.sat4j.moco.util.Real;
import org.sat4j.pb.IPBSolver;
import org.sat4j.specs.ContradictionException;
import org.sat4j.specs.IConstr;
import org.sat4j.specs.IVec;
import org.sat4j.specs.IVecInt;

/**
 * Representation of a Greater-or-Equal (GE) constraint.
 * @author Miguel Terra-Neves
 */
public class GE extends PBConstr {
    
    /**
     * Creates an instance of a GE constraint.
     * @param lits The literals in the constraint's left-hand side.
     * @param coeffs The coefficients in the constraint's left-hand side.
     * @param rhs The constraint's right-hand side.
     */
    GE(IVecInt lits, IVec<Real> coeffs, Real rhs) { super(lits, coeffs, rhs); }
    
    @Override
    protected Real getActivatorCoeff() {
        Real min_sum = getLHS().getMinSum();
        assert(!min_sum.isPositive());
        return min_sum.negate().add(getRHS());
    }
    
    @Override
    public boolean violatedBy(Real val) { return val.lessThan(getRHS()); }

    @Override
    IConstr addToSolver(IPBSolver solver) throws ContradictionException {
        if (getRHS().isPositive() && getLHS().allCoeffsGE(getRHS())) {
            return solver.addClause(getLits());
        }
        else if (getLHS().allUnitCoeffs()) {
            try {
                return solver.addAtLeast(getLits(), getRHS().asIntExact());
            }
            catch (ArithmeticException ae) {        // may occur if RHS is not integer
                return solver.addAtLeast(getLits(), getRHS().asInt()+1);
            }
        }
        return scaleToIntThenAddToSolver(solver);
    }
    
    @Override
    protected IConstr addScaledToSolver(IPBSolver solver, IVecInt lits, IVec<BigInteger> coeffs, BigInteger rhs)
            throws ContradictionException {
        return solver.addPseudoBoolean(lits, coeffs, true, rhs);
    }
    
    @Override
    public String getOpStrRep() { return ">="; }

}
