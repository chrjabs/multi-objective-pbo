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
package org.sat4j.moco.mcs;

import org.sat4j.core.ReadOnlyVecInt;
import org.sat4j.moco.pb.PBSolver;
import org.sat4j.moco.util.Log;
import org.sat4j.specs.IVecInt;

/**
 * Superclass for MCS extraction algorithms.
 * @author Miguel Terra-Neves
 */
public abstract class MCSAlgorithm {

    /**
     * Stores the MCS found by the algorithm if an MCS was found on the last successful call to
     * {@link #extract(PBSolver, IVecInt)}.
     */
    private ReadOnlyVecInt mcs = null;
    
    /**
     * Stores the MSS found by the algorithm if an MCS was found on the last successful call to
     * {@link #extract(PBSolver, IVecInt)}.
     */
    private ReadOnlyVecInt mss = null;
    
    /**
     * Stores the last model obtained by the MCS algorithm.
     * When the algorithm terminates, this should be a witness of the computed MCS.
     */
    private boolean[] model = null;
    
    /**
     * Stores if current model exploitation is enabled.
     * @see #enableExploitModel()
     */
    private boolean exploit_model = false;
    
    /**
     * Boolean used to store if the PB instance was solved successfully on the last call to
     * {@link #extract(PBSolver, IVecInt)}, i.e., an MCS was found or the hard formula was proven
     * unsatisfiable.
     */
    private boolean is_solved = false;
    
    /**
     * Boolean used to store if an MCS was found on the last call to {@link #extract(PBSolver, IVecInt)}.
     */
    private boolean found_mcs = false;
    
    /**
     * Stores a listener object that listens for models found during the MCS extraction process.
     */
    private IModelListener listener = new IModelListener() {
        public void onModel(PBSolver s) { /* default model listener does nothing */ }
    };
    
    /**
     * If an MCS was found on the last successful call to {@link #extract(PBSolver, IVecInt)}, retrieves that
     * MCS.
     * @return The MCS.
     */
    ReadOnlyVecInt getMCS() { return this.mcs; }
    
    /**
     * If an MCS was found on the last successful call to {@link #extract(PBSolver, IVecInt)}, retrieves the
     * corresponding MSS.
     * @return The MSS.
     */
    ReadOnlyVecInt getMSS() { return this.mss; }
    
    /**
     * Sets the listener that will listen for models found during the MCS extraction process.
     * @param l The model listener.
     */
    void setModelListener(IModelListener l) { this.listener = l; }
    
    /**
     * Enables current model exploitation.
     * This is a packaged protected method to be used in MCS extraction to hint the MCS algorithm that the
     * current model is also a model of the next hard formula and can be safely exploited (e.g. for
     * initialization purposes).
     */
    void enableExploitModel() { this.exploit_model = true; }
    
    /**
     * Disables current model exploitation.
     * @see #enableExploitModel()
     */
    void disableExploitModel() { this.exploit_model = false; }
    
    /**
     * Checks if current model exploitation is enabled.
     * @return True if exploitation is enabled, false otherwise.
     * @see #enableExploitModel()
     */
    protected boolean exploitModelEnabled() { return this.exploit_model; }

    /**
     * Extracts a model from a given PB solver, stores it and invokes the model listener.
     * To be used by {@link #extract(PBSolver, IVecInt)} whenever a model is found by the algorithm.
     * @param s The solver containing the model.
     */
    protected void newModel(PBSolver s) {
        this.listener.onModel(s);
        if (this.model == null || this.model.length < s.nVars()) {
            this.model = new boolean[s.nVars()];
        }
        for (int x = 1; x <= s.nVars(); ++x) {
            this.model[x-1] = s.modelValue(x);
        }
    }
    
    /**
     * If a model was produced on the last call to {@link #extract(PBSolver, IVecInt)}, retrieves the value
     * of a given literal in that model.
     * If an MCS was found, the model is a witness of that MCS.
     * @param lit The literal.
     * @return True if {@code lit} has value 1 in the model, false otherwise.
     */
    public boolean modelValue(int lit) {
        assert(lit != 0);
        return lit > 0 ? this.model[lit-1] : !this.model[-lit-1];
    }
    
    /**
     * Stores an MCS and the corresponding MSS, updating state to indicate that an MCS was found.
     * To be used by {@link #extract(PBSolver, IVecInt)} when an MCS is found.
     * @param mcs The MCS.
     * @param mss The corresponding MSS.
     */
    protected void saveMCS(IVecInt mcs, IVecInt mss) {
        Log.comment(3, "in MCSAlgorithm.saveMCS");
        this.mcs = new ReadOnlyVecInt(mcs);
        this.mss = new ReadOnlyVecInt(mss);
        this.is_solved = true;
        this.found_mcs = true;
        Log.comment(3, "out MCSAlgorithm.saveMCS");
    }
    
    /**
     * Updates state to indicate that the hard formula is unsatisfiable.
     * To be used by {@link #extract(PBSolver, IVecInt)} when the formula is proven unsatisfiable.
     */
    protected void setUnsat() {
        this.is_solved = true;
        this.found_mcs = false;
    }
    
    /**
     * Updates state to indicate that the formula was not solved successfully.
     * To be used by {@link #extract(PBSolver, IVecInt)} when something prevents the formula from being
     * solved (e.g. timeout).
     */
    protected void setUnsolved() { this.is_solved = false; }
    
    /**
     * Checks if the MCS algorithm was able to find an MCS or prove unsatisfiability on the last call to
     * {@link #extract(PBSolver, IVecInt)}.
     * @return True if the algorithm was successful, false otherwise.
     */
    boolean isSolved() { return this.is_solved; }
    
    /**
     * Checks if the MCS algorithm found an MCS on the last successful call to
     * {@link #extract(PBSolver, IVecInt)}.
     * @return True if an MCS was found, false otherwise.
     */
    boolean foundMCS() { return this.found_mcs; }
    
    /**
     * Checks if the MCS algorithm proved unsatisfiability of the hard formula on the last successful call to
     * {@link #extract(PBSolver, IVecInt)}.
     * @return True if the hard formula is unsatisfiable, false otherwise.
     */
    boolean isUnsat() { return !foundMCS(); }
    
    /**
     * Extracts an MCS from a given set of soft literals.
     * @param s A PB solver containing the hard formula and to be used by the MCS algorithm as an oracle.
     * @param lits The soft literals.
     */
    void extract(PBSolver s, IVecInt lits) { run(s, lits); }
    
    /**
     * Runs the actual MCS algorithm.
     * Must be implemented by subclasses.
     * @param s A PB solver containing the hard formula and to be used by the MCS algorithm as an oracle.
     * @param lits The soft literals.
     */
    protected abstract void run(PBSolver s, IVecInt lits);
    
    /**
     * Checks if the last MCS extraction was trivial for this algorithm.
     * This method is used by the MCS extractor during stratified extraction to decide if the remaining
     * partitions are too easy and should be merged into a single one.
     * By default, no extractions are trivial.
     * Subclasses may implement this method as a hint for the MCS extractor.
     * @return True if the last extraction was trivial, false otherwise.
     */
    boolean trivialExtraction() { return false; }
    
}
