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
package org.sat4j.moco.analysis;

import org.moeaframework.core.NondominatedPopulation;
import org.moeaframework.core.Population;
import org.moeaframework.core.Solution;
import org.moeaframework.core.Variable;
import org.moeaframework.core.variable.EncodingUtils;
import org.sat4j.moco.pb.PBSolver;
import org.sat4j.moco.problem.Instance;
import org.sat4j.moco.util.Clock;
import org.sat4j.moco.util.Log;

public class SubResult extends Result {
    /**
     * Extracts and stores the solution that corresponds to a model in
     * a given PB solver, without logging the costs
     * @param solver The solver
     */
     public SubResult(Instance problem){
	 super(problem);
    }

    @Override public void saveModel(PBSolver solver) {
        Solution sol = this.problem.newSolution();
        for (int lit = 1; lit <= sol.getNumberOfVariables(); ++lit) {
            Variable var = sol.getVariable(lit-1);
            EncodingUtils.setBoolean(var, solver.modelValue(lit));
        }
        this.problem.evaluate(sol);
        if (!sol.violatesConstraints() && !isWeaklyDominated(sol, this.solutions)) {
            this.solutions.add(sol);
            // Log.costs(sol.getObjectives());
            Log.comment(1, ":elapsed " + Clock.instance().getElapsed() + " :front-size " + nSolutions());
        }
    }

    /**
     * Extracts and stores the solution that corresponds to a model in
     * a given PB solver, without logging the costs
     * @param solver The solver.
     */
    @Override public void saveThisModel(boolean[] xModelValues ) {
        Solution sol = this.problem.newSolution();
        for (int lit = 1; lit <= sol.getNumberOfVariables(); ++lit) {
            Variable var = sol.getVariable(lit-1);
            EncodingUtils.setBoolean(var, xModelValues[lit-1]);
        }
        this.problem.evaluate(sol);
        if (!sol.violatesConstraints() && !isWeaklyDominated(sol, this.solutions)) {
            this.solutions.add(sol);
            // Log.costs(sol.getObjectives());
            Log.comment(1, ":elapsed " + Clock.instance().getElapsed() + " :front-size " + nSolutions());
        }
    }

}
