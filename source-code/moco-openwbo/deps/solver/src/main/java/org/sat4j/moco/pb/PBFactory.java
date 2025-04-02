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

import org.sat4j.core.Vec;
import org.sat4j.moco.util.Log;
import org.sat4j.moco.util.Real;
import org.sat4j.specs.IVec;
import org.sat4j.specs.IVecInt;

/**
 * Singleton class implementing a PB constraint factory.
 * The basic PB constraint types are supported by default.
 * Support for other constraint types can be added by providing new implementations of the {@link IPBBuilder}
 * interface through the {@link #addBuilder(IPBBuilder)} method.
 * Such constraints may then be instantiated using the generic
 * {@link #mkConstr(String, IVecInt, IVec, Real)} method.
 * @author Miguel Terra-Neves
 */
public class PBFactory {

    /**
     * The single instance of the PB constraint factory.
     */
    private static final PBFactory instance = new PBFactory();
    
    /**
     * Retrieves the singleton PB factory instance.
     * @return The PB factory.
     */
    public static PBFactory instance() { return instance; }
    
    /**
     * Stores the builders of the PB constraint types supported by the PB factory.
     */
    private IVec<IPBBuilder> builders = new Vec<IPBBuilder>();
    
    static {
        instance().addBuilder(new IPBBuilder() {
            public PBConstr mkConstr(String op, IVecInt lits, IVec<Real> coeffs, Real rhs) {
                if (op.equals("=")) { return instance().mkEQ(lits, coeffs, rhs); }
                else if (op.equals(">=") || op.equals("=>")) { return instance().mkGE(lits, coeffs, rhs); }
                else if (op.equals("<=") || op.equals("=<")) { return instance().mkLE(lits, coeffs, rhs); }
                return null;
            }
        });
        Log.comment(2, "standard constraint builders ready");
    }
    
    /**
     * Creates an instance of a PB factory.
     */
    private PBFactory() {}
    
    /**
     * Creates an instance of a PB expression.
     * @param lits The expression's literals.
     * @param coeffs The expression's coefficients.
     * @return The PB expression.
     */
    public PBExpr mkExpr(IVecInt lits, IVec<Real> coeffs) { return new PBExpr(lits, coeffs); }
    
    /**
     * Adds a custom PB constraint builder.
     * Use this method to add support for custom constraint types.
     * @param b The constraint builder.
     */
    public void addBuilder(IPBBuilder b) { this.builders.push(b); }
    
    /**
     * Creates an instance of a PB constraint.
     * @param op A string representation of the constraint's operator.
     * @param lits The literals in the constraint's left-hand side.
     * @param coeffs The coefficients in the constraint's left-hand side.
     * @param rhs The constraint's right-hand side.
     * @return The PB constraint object, or {@code null} if the factory does not support the constraint type
     * given by {@code op} ({@link #addBuilder(IPBBuilder)}).
     */
    public PBConstr mkConstr(String op, IVecInt lits, IVec<Real> coeffs, Real rhs) {
        PBConstr c = null;
        for (int i = this.builders.size()-1; i >= 0 && c == null; --i) {
            c = this.builders.get(i).mkConstr(op, lits, coeffs, rhs);
        }
        return c;
    }
    
    /**
     * Creates a vector of 1s with a given size.
     * @param size The vector's size.
     * @return A vector of 1s of size {@code size}.
     */
    private IVec<Real> mkUnitCoeffs(int size) {
        IVec<Real> c = new Vec<Real>(size);
        for (int i = 0; i < size; ++i) {
            c.unsafePush(Real.ONE);
        }
        return c;
    }
    
    /**
     * Creates an instance of an LE constraint.
     * @param lits The literals in the constraint's left-hand side.
     * @param coeffs The coefficients in the constraint's left-hand side.
     * @param rhs The constraint's right-hand side.
     * @return The LE constraint object.
     */
    public LE mkLE(IVecInt lits, IVec<Real> coeffs, Real rhs) { return new LE(lits, coeffs, rhs); }
    
    /**
     * Creates an instance of an LE constraint with all unit coefficients, also known as at-most-k
     * constraint.
     * @param lits The literals in the constraint's left-hand side.
     * @param rhs The constraint's right-hand side.
     * @return The LE constraint object.
     */
    public LE mkLE(IVecInt lits, int rhs) { return mkLE(lits, mkUnitCoeffs(lits.size()), new Real(rhs)); }
    
    /**
     * Creates an instance of a GE constraint.
     * @param lits The literals in the constraint's left-hand side.
     * @param coeffs The coefficients in the constraint's left-hand side.
     * @param rhs The constraint's right-hand side.
     * @return The GE constraint object.
     */
    public GE mkGE(IVecInt lits, IVec<Real> coeffs, Real rhs) { return new GE(lits, coeffs, rhs); }
    
    /**
     * Creates an instance of a GE constraint with all unit coefficients, also known as at-least-k
     * constraint.
     * @param lits The literals in the constraint's left-hand side.
     * @param rhs The constraint's right-hand side.
     * @return The GE constraint object.
     */
    public GE mkGE(IVecInt lits, int rhs) { return mkGE(lits, mkUnitCoeffs(lits.size()), new Real(rhs)); }
    
    /**
     * Creates an instance of an EQ constraint.
     * @param lits The literals in the constraint's left-hand side.
     * @param coeffs The coefficients in the constraint's left-hand side.
     * @param rhs The constraint's right-hand side.
     * @return The EQ constraint object.
     */
    public EQ mkEQ(IVecInt lits, IVec<Real> coeffs, Real rhs) { return new EQ(lits, coeffs, rhs); }
    
    /**
     * Creates an instance of a EQ constraint with all unit coefficients.
     * @param lits The literals in the constraint's left-hand side.
     * @param rhs The constraint's right-hand side.
     * @return The EQ constraint object.
     */
    public EQ mkEQ(IVecInt lits, int rhs) { return mkEQ(lits, mkUnitCoeffs(lits.size()), new Real(rhs)); }
    
    /**
     * Creates an instance of a clause constraint (disjunction of literals).
     * @param lits The clause's literals.
     * @return A {@link PBConstr} object representing the clause constraint.
     */
    public PBConstr mkClause(IVecInt lits) { return mkGE(lits, 1); }
    
}
