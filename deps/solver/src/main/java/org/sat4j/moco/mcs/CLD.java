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

import org.sat4j.core.Vec;
import org.sat4j.core.VecInt;
import org.sat4j.moco.pb.ConstrID;
import org.sat4j.moco.pb.PBFactory;
import org.sat4j.moco.pb.PBSolver;
import org.sat4j.moco.util.Log;
import org.sat4j.specs.ContradictionException;
import org.sat4j.specs.IVec;
import org.sat4j.specs.IVecInt;

/**
 * Class that implements the CLD algorithm for MCS extraction, proposed in:<br>
 *      Marques-Silva, J., Heras, F., Janota, M., Previti, A., &amp; Belov, A. (2013, August).
 *      On Computing Minimal Correction Subsets. In IJCAI (pp. 615-622).
 * @author Miguel Terra-Neves
 */
public class CLD extends MCSAlgorithm {
    
    /**
     * Stores the number of models found in the last run of this algorithm.
     * Used to determine if the last extraction was trivial.
     */
    private int models_found = 0;

    @Override
    protected void newModel(PBSolver s) {
        this.models_found++;
        super.newModel(s);
    }
    
    @Override
    protected void run(PBSolver solver, IVecInt lits) {
        Log.comment(3, "in CLD.extract");
        this.models_found = 0;
        IVecInt mcs = new VecInt();
        IVecInt mss = new VecInt();
        lits.copyTo(mcs);
        IVec<ConstrID> ids = new Vec<ConstrID>();
        boolean first = true, is_sat = true;
        while (is_sat) {
            if (!first) { newModel(solver); }
            if (!first || exploitModelEnabled()) { updateSat(mss, mcs); }
            first = false;
            Log.comment(2, ":sat " + mss.size() + " :undef " + mcs.size());
            if (mcs.isEmpty()) { break; }
            try {
                ids.push(solver.addRemovableConstr(PBFactory.instance().mkClause(mcs)));
                ids.push(solver.addRemovableConstr(PBFactory.instance().mkGE(mss, mss.size())));
                solver.check();
                Log.comment(3, "CLD SAT check :is-solved " + solver.isSolved() + " :is-sat " + solver.isSat());
                is_sat = solver.isSolved() && solver.isSat();
            }
            catch (ContradictionException e) {
                Log.comment(3, "contradiction adding clause D");
                is_sat = false;
            }
        }
        Log.comment(3, ":to-remove " + ids.size());
        solver.removeConstrs(ids);
        if (!exploitModelEnabled() && solver.isSolved() && mss.isEmpty()) {     // do MCS existence check
            solver.check();
            Log.comment(3, "MCS exists check :is-solved " + solver.isSolved() + " :is-sat " + solver.isSat());
            if (solver.isSolved() && solver.isSat()) { newModel(solver); }
        }
        if (!exploitModelEnabled() && solver.isSolved() && mss.isEmpty() && solver.isUnsat()) {
            setUnsat();
        }
        else if (solver.isSolved()) {
            assert(!is_sat || mcs.isEmpty());
            assert(!mss.isEmpty() || solver.isSat() || exploitModelEnabled());
            saveMCS(mcs, mss);
        }
        else {
            Log.comment(1, "CLD timeout");
            setUnsolved();
        }
        Log.comment(3, "out CLD.extract");
    }
    
    /**
     * Checks the literals in a vector of undefined literals that are satisfied by the current model.
     * Satisfied literals are moved to a vector of satisfied literals.
     * @param sat The vector of satisfied literals.
     * @param undef The vector of undefined literals.
     */
    private void updateSat(IVecInt sat, IVecInt undef) {
        Log.comment(3, "in CLD.updateSat");
        for (int i = 0; i < undef.size();) {
            if (modelValue(undef.get(i))) {
                sat.push(undef.get(i));
                undef.set(i, undef.last());
                undef.pop();
            }
            else {
                ++i;
            }
        }
        Log.comment(3, "out CLD.updateSat");
    }

    @Override
    boolean trivialExtraction() { return this.models_found == 0; }
    
}
