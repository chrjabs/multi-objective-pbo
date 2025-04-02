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
package org.sat4j.moco;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import org.junit.Before;
import org.junit.Test;
import org.sat4j.core.Vec;
import org.sat4j.core.VecInt;
import org.sat4j.moco.pb.ConstrID;
import org.sat4j.moco.pb.GE;
import org.sat4j.moco.pb.LE;
import org.sat4j.moco.pb.PBConstr;
import org.sat4j.moco.pb.PBFactory;
import org.sat4j.moco.pb.PBSolver;
import org.sat4j.moco.util.Real;
import org.sat4j.specs.ContradictionException;
import org.sat4j.specs.IVecInt;

public class PBSolverTest {

    private PBSolver solver;
    
    @Before
    public void setUp() {
        this.solver = new PBSolver();
        this.solver.newVars(3);
        this.solver.unsafeAddConstr(PBFactory.instance().mkLE(new VecInt(new int[] { 1, 2, 3 }), 2));
    }
    
    public void testPropagatedUnsatBug(PBConstr c) {
        this.solver.check();
        assertTrue(this.solver.isSolved() && this.solver.isSat() && !this.solver.isUnsat());
        ConstrID id = this.solver.unsafeAddRemovableConstr(c);
        this.solver.check();
        assertTrue(this.solver.isSolved() && !this.solver.isSat() && this.solver.isUnsat());
        this.solver.removeConstr(id);
        this.solver.check();
        assertTrue(this.solver.isSolved() && this.solver.isSat() && !this.solver.isUnsat());
    }
    
    @Test
    public void testPropagatedUnsatGE() {
        GE c = PBFactory.instance().mkGE(new VecInt(new int[] { 1, 2, 3 }),
                                         new Vec<Real>(new Real[] { Real.ONE, Real.ONE, new Real(2)}),
                                         new Real(4));
        testPropagatedUnsatBug(c);
    }
    
    @Test
    public void testPropagatedUnsatLE() {
        LE c = PBFactory.instance().mkLE(new VecInt(new int[] { 1, 2, 3 }),
                                         new Vec<Real>(new Real[] { Real.ONE.negate(),
                                                                    Real.ONE.negate(),
                                                                    new Real(2).negate()}),
                                         new Real(4).negate());
        testPropagatedUnsatBug(c);
    }
    
    @Test
    public void testUnsatExplanation() {
        IVecInt asms = new VecInt(new int[] { 1, 2, 3 });
        this.solver.check(asms);
        assertTrue(this.solver.isSolved() && !this.solver.isSat() && this.solver.isUnsat());
        IVecInt core = this.solver.unsatExplanation();
        assertTrue(core.size() == 3);
        assertTrue(core.contains(1) && core.contains(2) && core.contains(3));
    }
    
    @Test
    public void testEmptyUnsatExplanation() {
        GE c = PBFactory.instance().mkGE(new VecInt(new int[] { 1, 2, 3 }),
                                         new Vec<Real>(new Real[] { Real.ONE, Real.ONE, new Real(2)}),
                                         new Real(4));
        this.solver.unsafeAddConstr(c);
        this.solver.check();
        assertTrue(this.solver.isSolved() && !this.solver.isSat() && this.solver.isUnsat());
        IVecInt core = this.solver.unsatExplanation();
        assertNotNull(core);
        assertTrue(core.isEmpty());
    }
    
    @Test
    public void testRepeatedLits() {
        GE c = PBFactory.instance().mkGE(new VecInt(new int[] { -1, 2, 2, 3 }), 4);
        this.solver.unsafeAddConstr(c);
        this.solver.check();
        assertTrue(this.solver.isSolved() && this.solver.isSat() && !this.solver.isUnsat());
        assertFalse(this.solver.modelValue(1));
        assertTrue(this.solver.modelValue(2));
        assertTrue(this.solver.modelValue(3));
    }
    
    @Test
    public void testRemoval() {
        GE c1 = PBFactory.instance().mkGE(new VecInt(new int[] { 1, 2 }), 2);
        PBConstr c2 = PBFactory.instance().mkClause(new VecInt(new int[] { -2, -3 }));
        ConstrID id1 = this.solver.unsafeAddRemovableConstr(c1);
        this.solver.unsafeAddConstr(c2);
        this.solver.check();
        assertTrue(this.solver.isSolved() && this.solver.isSat() && !this.solver.isUnsat());
        assertTrue(this.solver.modelValue(1));
        assertTrue(this.solver.modelValue(2));
        assertFalse(this.solver.modelValue(3));
        this.solver.removeConstr(id1);
        GE c3 = PBFactory.instance().mkGE(new VecInt(new int[] { 2, 3 }), 2);
        try {
            this.solver.addConstr(c3);
        }
        catch (ContradictionException e) {
            return;     // PASSED TEST!!
        }
        this.solver.check();
        assertTrue(this.solver.isSolved() && !this.solver.isSat() && this.solver.isUnsat());
    }
    
}
