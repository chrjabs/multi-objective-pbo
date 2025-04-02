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

/**
 * Container class for the result of an execution of the MOCO solver.
 * Stores the nondominated solutions found by the solver and if that set of solutions is the Pareto front
 * of the MOCO instance.
 * @author Miguel Terra-Neves
 */
public class Result {

    /**
     * A MOEA framework representation of the MOCO instance.
     */
    protected MOCOProblem problem = null;

    /**
     * Stores the nondominated solutions found by the solver.
     */
    protected NondominatedPopulation solutions = null;

    /**
     * Boolean value indicating if the {@link #solutions} set is the Pareto front of the MOCO instance.
     */
    protected boolean is_opt = false;

    public Result(){}
    /**
     * Creates an instance of a container for the nondominated solutions found for a given MOCO instance.
     * @param m The instance.
     */
    public Result(Instance m) {
        this.problem = new MOCOProblem(m);
        this.solutions = new NondominatedPopulation();
    }

    /**
     * Creates an instance of a container for the nondominated solutions found for a given MOCO instance.
     * @param m The instance.
     */
    public Result(Instance m, boolean unsafe) {
	this(m);
	if(unsafe)
	    this.solutions = new NondominatedPopulationUnsafe();
}

    /**
     * Extracts and stores the solution that corresponds to a model in a given PB solver.
     * @param solver The solver.
     */
    public void saveModel(PBSolver solver) {
        Solution sol = this.problem.newSolution();
        for (int lit = 1; lit <= sol.getNumberOfVariables(); ++lit) {
            Variable var = sol.getVariable(lit-1);
            EncodingUtils.setBoolean(var, solver.modelValue(lit));
        }
        this.problem.evaluate(sol);
        if (!sol.violatesConstraints() && !isWeaklyDominated(sol, this.solutions)) {
            this.solutions.add(sol);
            Log.costs(sol.getObjectives());
            Log.comment(1, ":elapsed " + Clock.instance().getElapsed() + " :front-size " + nSolutions());
        }
    }

    /**
     * Extracts and stores the solution that corresponds to xModelValues
     * @param xModelValues, the model in (eval.x_1, ...,eval.x_n) format
     */
    public void saveThisModel(boolean[] xModelValues ) {
        Solution sol = this.problem.newSolution();
        for (int lit = 1; lit <= sol.getNumberOfVariables(); ++lit) {
            Variable var = sol.getVariable(lit-1);
            EncodingUtils.setBoolean(var, xModelValues[lit-1]);
        }
        this.problem.evaluate(sol);
        if (!sol.violatesConstraints() && !isWeaklyDominated(sol, this.solutions)) {
            this.solutions.add(sol);
            Log.costs(sol.getObjectives());
            Log.comment(1, ":elapsed " + Clock.instance().getElapsed() + " :front-size " + nSolutions());
        }
    }

    /**
     * Extracts and stores the solution that corresponds to a model in a given PB solver, without checking
     * @param solver The solver.
     */
    public void saveModelUnsafe(PBSolver solver) {
        Solution sol = this.problem.newSolution();
        for (int lit = 1; lit <= sol.getNumberOfVariables(); ++lit) {
            Variable var = sol.getVariable(lit-1);
            EncodingUtils.setBoolean(var, solver.modelValue(lit));
        }
        this.problem.evaluate(sol);
	this.solutions.add(sol);
	Log.costs(sol.getObjectives());
	Log.comment(1, ":elapsed " + Clock.instance().getElapsed() + " :front-size " + nSolutions());
    }

    /**
     * Extracts and stores the solution that corresponds to xModelValues, without checking.
     * @param  xModelValues
     */
    public void saveThisModelUnsafe(boolean[] xModelValues ) {
        Solution sol = this.problem.newSolution();
        for (int lit = 1; lit <= sol.getNumberOfVariables(); ++lit) {
            Variable var = sol.getVariable(lit-1);
            EncodingUtils.setBoolean(var, xModelValues[lit-1]);
        }
        this.problem.evaluate(sol);
            this.solutions.add(sol);
            Log.costs(sol.getObjectives());
            Log.comment(1, ":elapsed " + Clock.instance().getElapsed() + " :front-size " + nSolutions());
    }
    
    /**
     * Checks if a given solution is weakly dominated by some other solution in a given population.
     * @param sol The solution.
     * @param pop The population.
     * @return True if {@code sol} is weakly dominated by some solution in {@code pop}, false otherwise.
     */
    protected boolean isWeaklyDominated(Solution sol, Population pop) {
        for (int i = 0; i < pop.size(); ++i) {
            if (isWeaklyDominated(sol, pop.get(i))) {
                return true;
            }
        }
        return false;
    }

    /**
     * Checks if a given solution is weakly dominated by some other solution.
     * @param sol The solution to check.
     * @param other_sol The solution used for dominance check.
     * @return True if {@code sol} is weakly dominated by {@code other_sol}, false otherwise.
     */
    private boolean isWeaklyDominated(Solution sol, Solution other_sol) {
        assert(!sol.violatesConstraints() && !other_sol.violatesConstraints());
        for (int i = 0; i < other_sol.getNumberOfObjectives(); ++i) {
            if (sol.getObjective(i) < other_sol.getObjective(i)) {
                return false;
            }
        }
        return true;
    }

    /**
     * Method used to notify the container that its set of nondominated solutions is the Pareto front of the
     * MOCO instance.
     */
    public void setParetoFrontFound() { this.is_opt = true; }

    /**
     * Checks if the set of nondominated solutions stored in the container is the Pareto front of the MOCO
     * instance.
     * @return True if the solution set is the Pareto front, false otherwise.
     */
    public boolean isParetoFront() { return this.is_opt; }

    /**
     * Retrieves the number of nondominated solutions in the container.
     * @return The number of nondominated solutions.
     */
    public int nSolutions() { return this.solutions.size(); }

    /**
     * Retrieves the nondominated solutions in the container.
     * @return The nondominated solutions.
     */
    NondominatedPopulation getSolutions() { return this.solutions; }

    /**
     * Retrieves a given nondominated solution in the container.
     * @param i The solution index.
     * @return The {@code i}-th nondominated solution.
     */
    public Solution getSolution(int i) { return getSolutions().get(i); }

    /**
     * Adds a solution to the container and updates the set of nondominated solutions.
     * If the solution is dominated, it is discarded.
     * If not, it is added to the set and now dominated solutions are discarded.
     * @param s The solution.
     */
     void addSolution(Solution s) { this.solutions.add(s); }

    /**
     * Adds a solution to the container and updates the set of nondominated solutions.
     * If the solution is dominated, it is discarded.
     * If not, it is added to the set and now dominated solutions are discarded.
     * @param s The solution.
     */
    public void addSolutionUnsafe(Solution sol) { 
	this.problem.evaluate(sol);
	this.solutions.add(sol);
	Log.costs(sol.getObjectives());
	Log.comment(1, ":elapsed " + Clock.instance().getElapsed() + " :front-size " + nSolutions());
	; }

    /**
     * Retrieves the assignment of a given nondominated solution in the container.
     * @param i The solution index.
     * @return The {@code i}-th nominated solution's assignment.
     */
    public boolean[] getAssignment(int i) {
        assert(i < nSolutions());
        return this.problem.getAssignment(this.solutions.get(i));
    }

    /**
     * Retrieves the cost vector of a given nondominated solution in the container.
     * @param i The solution index.
     * @return The {@code i}-th nondominated solution's cost vector.
     */
    public double[] getCosts(int i) {
        Solution sol = this.solutions.get(i);
        assert(sol.getNumberOfObjectives() == this.problem.getNumberOfObjectives());
        double[] c = new double[this.problem.getNumberOfObjectives()];
        for (int j = 0; j < sol.getNumberOfObjectives(); ++j) {
            c[j] = sol.getObjective(j);
        }
        return c;
    }

    /**
     *Deletes every element in population
     */

   public void clearPopulation(){
	this.solutions.clear();
    }

}
