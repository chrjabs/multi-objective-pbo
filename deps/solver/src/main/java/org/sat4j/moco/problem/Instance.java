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

import java.io.IOException;
import java.io.Writer;

import org.sat4j.core.ReadOnlyVec;
import org.sat4j.core.ReadOnlyVecInt;
import org.sat4j.core.Vec;
import org.sat4j.moco.pb.PBConstr;
import org.sat4j.moco.pb.PBExpr;
import org.sat4j.moco.util.IOUtils;
import org.sat4j.moco.util.Real;
import org.sat4j.specs.IVec;

/**
 * Representation of a MOCO instance.
 * @author Miguel Terra-Neves
 */
public class Instance {
    
    /**
     * Stores the MOCO instance's constraints.
     */
    private IVec<PBConstr> constrs = null;
    
    /**
     * Stores the MOCO instance's objectives.
     */
    private IVec<Objective> objs = null;
    
    /**
     * The number of variables in the MOCO instance.
     */
    private int nvars = -1;

    /**
     * Creates an empty MOCO instance.
     */
    public Instance() {
        this.constrs = new Vec<PBConstr>();
        this.objs = new Vec<Objective>();
    }
    
    /**
     * Creates a MOCO instance initialized with a given set of constraints and objectives.
     * @param constrs The constraints.
     * @param objs The objectives.
     */
    public Instance(IVec<PBConstr> constrs, IVec<Objective> objs) {
        this();
        constrs.copyTo(this.constrs);
        objs.copyTo(this.objs);
    }
    
    /**
     * Retrieves the number of constraints in the MOCO instance.
     * @return The number of constraints in the MOCO instance.
     */
    public int nConstrs() { return this.constrs.size(); }
    
    /**
     * Retrieves the number of objectives in the MOCO instance.
     * @return The number of objectives in the MOCO instance.
     */
    public int nObjs() { return this.objs.size(); }
    
    /**
     * Counts the number of variables in the MOCO instance.
     * @return The number of variables in the MOCO instance.
     */
    private int countVars() {
        int n = 0;
        for (int i = 0; i < nConstrs(); ++i) {
            ReadOnlyVecInt lits = getConstr(i).getLHS().getLits();
            for (int j = 0; j < lits.size(); ++j) {
                int x = Math.abs(lits.get(j));
                n = x > n ? x : n;
            }
        }
	for (int i = 0; i < nObjs(); ++i) {
	    Objective ithObjective = getObj(i);
	    for (int k = 0; k < ithObjective.nSubObj(); ++k) {
		ReadOnlyVecInt lits = ithObjective.getSubObj(k).getLits();
		for (int j = 0; j < lits.size(); ++j) {
		    int x = Math.abs(lits.get(j));
		    n = x > n ? x : n;
		}
	    }
	}
return n;
    }
    
    /**
     * Retrieves the number of variables in the MOCO instance.
     * @return The number of variables in the MOCO instance.
     */
    public int nVars() { return nvars = (nvars >= 0 ? nvars : countVars()); }
    
    /**
     * Retrieves a constraint in the MOCO instance.
     * @param i The constraint's index.
     * @return The {@code i}-th constraint.
     */
    public PBConstr getConstr(int i) { return this.constrs.get(i); }
    
    /**
     * Retrieves an objective in the MOCO instance.
     * @param i The objective's index.
     * @return The {@code i}-th objective.
     */
    public Objective getObj(int i) { return this.objs.get(i); }
    
    /**
     * Adds a constraint to the MOCO instance.
     * @param c The constraint.
     */
    public void addConstr(PBConstr c) {
        this.nvars = -1;
        this.constrs.push(c);
    }
    
    /**
     * Adds an objective to the MOCO instance.
     * @param o The objective.
     */
    public void addObj(Objective o) {
        this.nvars = -1;
        this.objs.push(o);
    }
    
    /**
       Removes an objective
     */
       public void removeObj(int i){
	   this.objs.delete(i);
}


    /**
     * Writes the MOCO instance to a given {@link Writer} object.
     * @param w The writer object.
     * @throws IOException if an error occurs writing the MOCO instance.
     */
    public void dump(Writer w) throws IOException {
        for (int i = 0; i < nObjs(); ++i) {
            w.write("min: ");
            if (getObj(i) instanceof LinearObj) {
                writeExpr(w, ((LinearObj)getObj(i)).getExpr());
            }
            else {
                assert(getObj(i) instanceof DivObj);
                DivObj obj = (DivObj)getObj(i);
                for (int j = 0; j < obj.nDivs(); ++j) {
                    w.write("(");
                    writeExpr(w, obj.getNum(j));
                    w.write(") / (");
                    writeExpr(w, obj.getDen(j));
                    w.write(j == obj.nDivs()-1 ? ")" : ") + ");
                }
            }
            w.write(";" + IOUtils.NEWLINE);
        }
        for (int i = 0; i < nConstrs(); ++i) {
            PBConstr c = getConstr(i);
            writeExpr(w, c.getLHS());
            w.write(" " + c.getOpStrRep() + " " + c.getRHS().toPlainString() + ";" + IOUtils.NEWLINE);
        }
    }
    
    /**
     * Writes a given expression to a given {@link #Writer} object.
     * @param w The writer object.
     * @param e The expression.
     * @throws IOException if an error occurs writing the expression.
     */
    private void writeExpr(Writer w, PBExpr e) throws IOException {
        ReadOnlyVecInt lits = e.getLits();
        ReadOnlyVec<Real> coeffs = e.getCoeffs();
        for (int i = 0; i < lits.size(); ++i) {
            int l = lits.get(i);
            w.write(coeffs.get(i).toPlainString() + (l > 0 ? " x" : " ~x") + Math.abs(l) +
                    (i < lits.size()-1 ? " " : ""));
        }
    }
    
}
