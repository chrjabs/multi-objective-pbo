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
import org.sat4j.moco.pb.PBExpr;
import org.sat4j.moco.pb.PBSolver;
import org.sat4j.moco.util.Real;
import org.sat4j.specs.IVec;
import org.sat4j.specs.IVecInt;

/**
 * Representation of a linear PB objective.
 * @author Miguel Terra-Neves
 */
public class LinearObj extends Objective {

    /**
     * Creates an instance of a linear PB objective.
     * @param e The objective's PB expression.
     */
    public LinearObj(PBExpr e) { super(e); }
    
    /**
     * Creates an instance of a linear PB objective.
     * @param lits The objective's literals.
     * @param coeffs The objective's coefficients.
     */
    public LinearObj(IVecInt lits, IVec<Real> coeffs) { this(new PBExpr(lits, coeffs)); }
    
    /**
     * Retrieves the PB expression that represents the linear PB objective.
     * @return The objective's expression.
     */
    public PBExpr getExpr() { return getSubObj(0); }
    
    /**
     * Retrieves the PB objective's literals.
     * @return The objective's literals.
     */
    public ReadOnlyVecInt getLits() { return getExpr().getLits(); }
    
    /**
     * Retrieves the PB objective's coefficients.
     * @return The objective's coefficients.
     */
    public ReadOnlyVec<Real> getCoeffs() { return getExpr().getCoeffs(); }

    @Override
    public Real evaluate(boolean[] a) { return getExpr().evaluate(a); }

    @Override
    public Real evaluate(PBSolver s) { return getExpr().evaluate(s); }
    
}
