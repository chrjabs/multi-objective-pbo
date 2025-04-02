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
import org.sat4j.core.Vec;
import org.sat4j.core.VecInt;
import org.sat4j.moco.Params;
import org.sat4j.moco.pb.ConstrID;
import org.sat4j.moco.pb.PBFactory;
import org.sat4j.moco.pb.PBSolver;
import org.sat4j.moco.util.Clock;
import org.sat4j.moco.util.Log;
import org.sat4j.specs.IVec;
import org.sat4j.specs.IVecInt;

/**
 * Class that adds general MCS techniques (e.g. stratification) on top of the {@link MCSAlgorithm} class.
 * @author Miguel Terra-Neves
 */
// TODO: multiple MCS algorithms
// TODO: disjoint cores
public class MCSExtractor {

    /**
     * The MCS algorithm to be used for extraction.
     */
    private MCSAlgorithm alg = null;
    
    /**
     * The PB solver to be used as an oracle of the MCS algorithm.
     */
    private PBSolver solver = null;
    
    /**
     * Stores the MCS found on the last successful call to {@link #extract(IVecInt)} or
     * {@link #extract(IVec)}.
     */
    private ReadOnlyVecInt mcs = null;
    
    /**
     * Stores the MSS found on the last successful call to {@link #extract(IVecInt)} or
     * {@link #extract(IVec)}.
     */
    private ReadOnlyVecInt mss = null;
    
    /**
     * Stratified MCS extraction parameter that controls the maximum conflicts allowed before merging a
     * partition with the next one.
     */
    private int part_max_confl = 200000;
    
    /**
     * Stratified MCS exrtaction parameter that controls the number of trivially solved partitions that can
     * occur in a row before merging the remaining ones.
     */
    private int trivial_thres = 20;
    
    /**
     * Creates an instance of an MCS extractor.
     * @param s A PB solver, containing a hard formula, to be used as an oracle for the MCS algorithm.
     */
    public MCSExtractor(PBSolver s) {
        this.alg = new CLD();
        this.solver = s;
    }
    
    /**
     * If an MCS was found on the last successful call to {@link #extract(IVecInt)} or {@link #extract(IVec)},
     * retrieves that MCS.
     * @return The MCS.
     */
    public ReadOnlyVecInt getMCS() { return this.mcs; }
    
    /**
     * If an MCS was found on the last successful call to {@link #extract(IVecInt)} or {@link #extract(IVec)},
     * retrieves the corresponding MSS.
     * @return The MSS.
     */
    public ReadOnlyVecInt getMSS() { return this.mss; }
    
    /**
     * Sets the listener that will listen for models found during the execution of the MCS algorithm.
     * @param l The model listener.
     */
    public void setModelListener(IModelListener l) { this.alg.setModelListener(l); }
    
    /**
     * Extracts an MCS from a given set of soft literals.
     * @param undef The soft literals.
     */
    public void extract(IVecInt undef) {
        Log.comment(3, "in MCSExtractor.extract");
        this.alg.extract(this.solver, undef);
        if (isSolved() && foundMCS()) {
            saveMCS(this.alg.getMCS(), this.alg.getMSS());
        }
        Log.comment(3, "out MCSExtractor.extract");
    }
    
    /**
     * Performs stratified MCS extraction given a sequence of partitions of soft literals.
     * @param undef The soft literal partition sequence.
     */
    // TODO: does not re-use MCS algorithm's partial progress in case of conflict timeout; may be worth it to re-use?
    public void extract(IVec<IVecInt> undef) {
        if (undef.size() == 1) {
            extract(undef.get(0));
            return;
        }
        Log.comment(3, "in MCSExtractor.extract");
        IVecInt mcs = new VecInt();
        IVecInt mss = new VecInt();
        IVec<ConstrID> ids = new Vec<ConstrID>();
        IVec<IVecInt> undef_cpy = new Vec<IVecInt>();
        undef.copyTo(undef_cpy);    // ensure that changes to undef aren't visible outside
        int ntrivial = 0;
        getSolver().setMaxConflicts(this.part_max_confl);
        for (int i = 0; i < undef_cpy.size(); ++i) {
            IVecInt part = undef_cpy.get(i);
            if (ntrivial >= this.trivial_thres) {
                Log.comment(1, "trivial extraction threshold reached, merging last " + (undef_cpy.size() - i) +
                               " partitions");
                part = new VecInt();
                for (; i < undef_cpy.size(); ++i) {
                    undef_cpy.get(i).copyTo(part);
                }
                --i;
            }
            Log.comment(1, ":partition " + i + " :size " + part.size());
            if (i == undef_cpy.size()-1) { getSolver().resetMaxConflicts(); }
            if (i > 0 && (mcs.size() > 0 || mss.size() > 0)) { this.alg.enableExploitModel(); }
            this.alg.extract(getSolver(), part);
            if (isSolved() && foundMCS()) {
                ReadOnlyVecInt part_mcs = this.alg.getMCS();
                ReadOnlyVecInt part_mss = this.alg.getMSS();
                assert(allSatisfied(part_mss));
                part_mcs.copyTo(mcs);
                part_mss.copyTo(mss);
                ids.push(getSolver().unsafeAddRemovableConstr(PBFactory.instance().mkLE(part_mcs, 0)));
                ids.push(getSolver().unsafeAddRemovableConstr(PBFactory.instance().mkGE(part_mss, part_mss.size())));
                ntrivial = this.alg.trivialExtraction() ? ntrivial+1 : 0;
            }
            else if (!isSolved() && !Clock.instance().timedOut()) {
                assert(i < undef_cpy.size()-1);
                Log.comment(2, "conflict timeout at partition " + i);
                IVecInt new_part = new VecInt();
                part.copyTo(new_part);
                undef_cpy.get(i+1).copyTo(new_part);
                undef_cpy.set(i+1, new_part);
            }
            else {
                assert((!isSolved() && Clock.instance().timedOut()) || (isSolved() && isUnsat() && i == 0));
                break;
            }
        }
        if (isSolved() && foundMCS()) {
            saveMCS(mcs, mss);
        }
        Log.comment(3, ":to-remove " + ids.size());
        getSolver().removeConstrs(ids);
        getSolver().resetMaxConflicts();
        this.alg.disableExploitModel();
        Log.comment(3, "out MCSExtractor.extract");
    }
    
    /**
     * Checks if the all the literals in a given vector are satisfied by the MCS algorithm's current model.
     * @param lits The literals.
     * @return True if all the literals in {@code lits} are satisfied, false otherwise.
     */
    private boolean allSatisfied(IVecInt lits) {
        assert(isSolved() && foundMCS());
        for (int i = 0; i < lits.size(); ++i) {
            if (!this.alg.modelValue(lits.get(i))) {
                return false;
            }
        }
        return true;
    }
    
    /**
     * Stores an MCS and the corresponding MSS.
     * @param mcs The MCS.
     * @param mss The corresponding MSS.
     */
    private void saveMCS(IVecInt mcs, IVecInt mss) {
        Log.comment(1, ":mss-size " + mss.size() + " :mcs-size " + mcs.size());
        this.mcs = new ReadOnlyVecInt(mcs);
        this.mss = new ReadOnlyVecInt(mss);
    }
    
    /**
     * Checks if the MCS extractor was able to find an MCS or prove unsatisfiability on the last call to
     * {@link #extract(IVecInt)} or {@link #extract(IVec)}.
     * @return True if the extractor was successful, false otherwise.
     */
    public boolean isSolved() { return this.alg.isSolved(); }
    
    /**
     * Checks if the MCS extractor proved unsatisfiability of the hard formula on the last successful call to
     * {@link #extract(IVecInt)} or {@link #extract(IVec)}.
     * @return True if the hard formula is unsatisfiable, false otherwise.
     */
    public boolean isUnsat() { return this.alg.isUnsat(); }
    
    /**
     * Checks if the MCS extractor found an MCS on the last successful call to {@link #extract(IVecInt)} or
     * {@link #extract(IVec)}.
     * @return True if an MCS was found, false otherwise.
     */
    public boolean foundMCS() { return this.alg.foundMCS(); }
    
    /**
     * Retrieves the PB solver being used as an oracle for the MCS algorithm.
     * @return The oracle.
     */
    public PBSolver getSolver() { return this.solver; }
    
    /**
     * Sets the MCS extractor configuration to the one stored in a given set of parameters.
     * @param p The parameters object.
     */
    public void updtParams(Params p) {
        this.part_max_confl = p.getPartMaxConfl();
        this.trivial_thres = p.getTrivialThres();
        Log.comment(":part-max-confl " + this.part_max_confl);
        Log.comment(":trivial-thres " + this.trivial_thres);
    }
    
}
