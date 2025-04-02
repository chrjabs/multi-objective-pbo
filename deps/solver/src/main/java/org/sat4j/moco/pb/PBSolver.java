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

import java.io.PrintWriter;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Timer;
import java.util.TimerTask;
import java.util.Map.Entry;

import org.sat4j.core.VecInt;
import org.sat4j.minisat.learning.ClauseOnlyLearning;
import org.sat4j.minisat.restarts.LubyRestarts;
import org.sat4j.minisat.restarts.NoRestarts;
import org.sat4j.moco.util.Clock;
import org.sat4j.pb.IPBSolver;
import org.sat4j.pb.SolverFactory;
import org.sat4j.pb.core.PBDataStructureFactory;
import org.sat4j.pb.core.PBSolverResolution;
import org.sat4j.specs.ContradictionException;
import org.sat4j.specs.IVec;
import org.sat4j.specs.IVecInt;
import org.sat4j.specs.TimeoutException;

/**
 * Solver wrapper around SAT4J for solving Pseudo-Boolean problems.
 * The wrapper handles some bugs and undocumented behavior, making SAT4J easier to use.
 * Works incrementally.
 * @author Miguel Terra-Neves
 */
public class PBSolver {
    
    /**
     * An instance of the underlying SAT4J PB solver.
     */
    private IPBSolver solver = null;

    /**
     * Boolean used to store the satisfiability of the PB instance of the last call to {@link #check()}
     * or {@link #check(IVecInt)}.
     * True if the instance is satisfiable, false otherwise.
     */
    private boolean is_sat = false;
    
    /**
     * Boolean used to store if the PB instance was solved successfully on the last call to {@link #check()}
     * or {@link #check(IVecInt)}.
     * True if so, false otherwise.
     */
    private boolean is_solved = false;
    
    /**
     * Stores mapping of removable constraint IDs to corresponding activator literals.
     */
    private Map<ConstrID, Integer> act_map = new HashMap<ConstrID, Integer>();
    
    /**
     * Maximum conflicts allowed on a single call to {@link #check()} or {@link #check(IVecInt)}.
     * If {@code max_conflicts} is smaller than 0, then no conflict limit is imposed.
     */
    private int max_conflicts = -1;
    
    /**
     *Number of clauses
     */

    private int clausesN = 0;
    
    /**
     *Constant variable;
     */
    private int constantID ;

    public int getClausesN(){
	return this.clausesN;
    }

    /**
     * Creates an instance of a PB solver.
     */
    public PBSolver() {
	PBSolverResolution solver = SolverFactory.newResolutionGlucose21();
	solver.setLearningStrategy(new ClauseOnlyLearning<PBDataStructureFactory>());
	solver.setRestartStrategy(new NoRestarts());
	solver.setSimplifier(solver.SIMPLE_SIMPLIFICATION);
        solver.setLearnedConstraintsDeletionStrategy(solver.activity_based_low_memory);
	this.solver = solver;
	this.newVar();

    }
    
    /**
     *Return the ID of a freshly created variable
     */

    public int getFreshVar(){
	this.newVar();
	return this.nVars();
    }



    /**
     * Creates a new Boolean variable in the PB solver.
     */
    public void newVar() { newVars(1); }
    
    /**
     * Creates multiple new Boolean variables in the PB solver.
     * @param nvars The number of variables to be created.
     */
    public void newVars(int nvars) { this.solver.newVar(this.solver.nVars() + nvars); }
    
    /**
     * Retrieves the number of variables created in the PB solver.
     * PB constraints cannot contain variables with identifiers larger than the number returned by this
     * method.
     * @return The number of variables in the solver.
     */
    public int nVars() { return this.solver.nVars(); }
    
    /**
     * Sets the maximum number of conflicts allowed for future calls to {@link #check()} or
     * {@link #check(IVecInt)}.
     * By default, no limit is imposed until this method is called.
     * @param conflicts The number of conflicts.
     */
    public void setMaxConflicts(int conflicts) { this.max_conflicts = conflicts; }
    
    /**
     * Disables the maximum conflicts limit set with {@link #setMaxConflicts(int)}.
     */
    public void resetMaxConflicts() { this.max_conflicts = -1; }
    
    /**
     * Checks if the PB solver has a maximum conflict limit set.
     * @return True if a conflict limit is set, false otherwise.
     * @see #setMaxConflicts(int)
     */
    private boolean hasMaxConflicts() { return this.max_conflicts >= 0; }
    
    /**
     * Retrieves the maximum number of conflicts allowed for a single call to {@link #check()} or
     * {@link #check(IVecInt)}.
     * @return The conflict limit, or a value smaller than 0 if no limit is set.
     */
    private int getMaxConflicts() { return this.max_conflicts; }
    
    /**
     * Adds an activator literal to a given constraint.
     * @param c The constraint.
     * @return The activator literal.
     */
    private int addActivator(PBConstr c) {
        newVar();
        int act = nVars();
        c.setActivator(act);
        return act;
    }
    
    /**
     * Adds a PB constraint to the PB solver.
     * @param c The constraint.
     * @throws ContradictionException if the solver detects that the addition of {@code c} would cause the
     * formula to become unsatisfiable.
     */

    public void addConstr(PBConstr c) throws ContradictionException { c.addToSolver(this.solver); }
    
    /**
     * Adds a removable PB constraint to the PB solver.
     * @param c The constraint.
     * @return A constraint ID object that can be used to remove the constraint in the future through
     * {@link #removeConstr(ConstrID)} or {@link #removeConstrs(IVec)}.
     * @throws ContradictionException if the solver detects that the addition of {@code c} would cause the
     * formula to become unsatisfiable.
     */
    public ConstrID addRemovableConstr(PBConstr c) throws ContradictionException {
        int act = addActivator(c);
        ConstrID id = ConstrID.mkFresh();
        this.act_map.put(id, act);
        c.addToSolver(this.solver);
        return id;
    }
    
    /**
     * Adds a PB constraint to the PB solver.
     * Only use this method instead of {@link #addConstr(PBConstr)} if known that the constraint will not
     * cause the formula to become unsatisfiable.
     * @param c The constraint.
     */
    public void unsafeAddConstr(PBConstr c) {
        try {
            addConstr(c);
        }
        catch (ContradictionException e) {
            throw new RuntimeException("ContradictionException thrown during unsafe constraint addition", e);
        }
    }
    
    /**
     * Adds a removable PB constraint to the PB solver.
     * Only use this method instead of {@link #addRemovableConstr(PBConstr)} if known that the constraint
     * will not cause the formula to become unsatisfiable.
     * @param c The constraint.
     * @return A constraint ID object that can be used to remove the constraint in the future through
     * {@link #removeConstr(ConstrID)} or {@link #removeConstrs(IVec)}.
     */
    public ConstrID unsafeAddRemovableConstr(PBConstr c) {
        try {
            return addRemovableConstr(c);
        }
        catch (ContradictionException e) {
            throw new RuntimeException("ContradictionException thrown during unsafe constraint addition", e);
        }
    }
    
    /**
     * Removes a PB constraint from the solver.
     * @param id The constraint's ID.
     */
    public void removeConstr(ConstrID id) {
        int act = this.act_map.get(id);
        this.act_map.remove(id);
        try {
            this.solver.addClause(new VecInt(new int[] { -act }));
        }
        catch (ContradictionException e) { /* only occurs if activator was added to empty clause */ }
    }
    
    /**
     * Removes a set of PB constraints from the solver.
     * @param ids A vector with the IDs of the constraints.
     */
    public void removeConstrs(IVec<ConstrID> ids) {
        for (int i = 0; i < ids.size(); ++i) {
            removeConstr(ids.get(i));
        }
    }
    
    /**
     * Checks the satisfiability of the formula in the PB solver under a given set of assumptions.
     * The assumptions are a set of literals that also must be satisfied in addition to the formula.
     * The result of the satisfiability check can be retrieved through the method {@link #isSolved()},
     * {@link #isSat()} and {@link #isUnsat()}.
     * If the formula is satisfiable, the satisfying assignment can be retrieved through
     * {@link #modelValue(int)}.
     * If the formula is unsatisfiable, a subset of the assumptions responsible for unsatisfiability can be
     * retrieved through {@link #unsatExplanation()}.
     * @param asms A vector of assumptions.
     */
    public void check(IVecInt asms) {
        this.is_solved = false;
        IVecInt act_asms = new VecInt(asms.size() + this.act_map.size());
        asms.copyTo(act_asms);
        for (Iterator<Integer> it = this.act_map.values().iterator(); it.hasNext();) {
            act_asms.unsafePush(it.next());
        }
        if (Clock.instance().timedOut()) { return; }
        this.solver.expireTimeout();
        this.solver.setTimeout(Clock.instance().getRemaining());
        Timer timer = new Timer();                  // FIXME: timer required because SAT4J does not support conflict
        if (hasMaxConflicts()) {                    // and time timeouts simultaneously
            timer.schedule(new TimerTask() {
                @Override
                public void run() {
                    solver.expireTimeout();
                }
            }, 1000*(long)Clock.instance().getRemaining());
            this.solver.setTimeoutOnConflicts(getMaxConflicts());
        }
        try {
            this.is_sat = this.solver.isSatisfiable(act_asms);
            this.is_solved = true;
        }
        catch (TimeoutException e) { /* intentionally left empty */ }
        catch (NullPointerException e) {    // FIXME: DIRTY HACK!!! SAT4J is buggy and might throw a
            /* intentionally left empty */  // NullPointerException when expiring timeout from outside;
        }                                   // might lead to unexpected behavior if solver is re-used
        timer.cancel();
        timer.purge();
    }
    
    /**
     * Checks the satisfiability of the formula in the PB solver.
     * The result of the satisfiability check can be retrieved through the method {@link #isSolved()},
     * {@link #isSat()} and {@link #isUnsat()}.
     * If the formula is satisfiable, the satisfying assignment can be retrieved through
     * {@link #modelValue(int)}.
     */
    public void check() { check(new VecInt()); }
    
    /**
     * Checks if the solver was able to check the satisfiability on the last call to {@link #check()} or
     * {@link #check(IVecInt)}.
     * @return True if the check was successful, false otherwise.
     */
    public boolean isSolved() { return this.is_solved; }
    
    /**
     * Checks if the formula checked on the last successful call to {@link #check()} or
     * {@link #check(IVecInt)} is satisfiable.
     * @return True if the formula is satisfiable, false otherwise.
     */
    public boolean isSat() { return this.is_sat; }
    
    /**
     * Checks if the formula checked on the last successful call to {@link #check()} or
     * {@link #check(IVecInt)} is unsatisfiable.
     * @return True if the formula is unsatisfiable, false otherwise.
     */
    public boolean isUnsat() { return !isSat(); }
    
    /**
     * If the formula checked on the last successful call to {@link #check()} or {@link #check(IVecInt)} is
     * satisfiable, retrieves the value of a given literal in the satisfying assignment found by the solver.
     * @param lit The literal.
     * @return True if {@code lit} has value 1 in the satisfying assignment, false otherwise.
     */
    public boolean modelValue(int lit) {
        assert(lit != 0);
        return lit > 0 ? this.solver.model(lit) : !this.solver.model(-lit);
    }
    
    /**
     * If the formula checked on the last successful call to {@link #check()} or {@link #check(IVecInt)} is
     * unsatisfiable, retrieves a subset of the assumptions that are responsible for unsatisfiability.
     * @return A vector with the assumption literals responsible for unsatisfiability.
     */
    public IVecInt unsatExplanation() {
        IVecInt explanation = this.solver.unsatExplanation();
        if (explanation != null) {
            for (int i = 0; i < explanation.size();) {
                if (this.act_map.containsValue(explanation.get(i))) {
                    explanation.set(i, explanation.last());
                    explanation.pop();
                }
                else {
                    ++i;
                }
            }
            return explanation;
        }
        return new VecInt();
    }
    /**
     * Is the literal positive?
     * @param literal
     */
    public boolean isLiteralPositive(int literal){
	return literal > 0;

}

    /**
     *Adds the disjunction of setOfLiterals
     *@param setOfliterals
     */

    public void AddClause(IVecInt setOfLiterals) throws ContradictionException{ 
	this.clausesN++;
	this.solver.addClause(setOfLiterals);
    }


    /**
     * Adds a removable clause to the PB solver.
     * @param setOfLiterals
     * @return A constraint ID object that can be used to remove the
     * clause in the future through {@link #removeConstr(ConstrID)} or
     * {@link #removeConstrs(IVec)}.
     * @throws ContradictionException if the solver detects that the
     * addition of {@code setOfLiterals} would cause the formula to
     * become unsatisfiable.
     */

    public ConstrID addRemovableClause(IVecInt setOfLiterals) throws ContradictionException {
	this.newVar();
        int act = this.nVars();
        ConstrID id = ConstrID.mkFresh();
        this.act_map.put(id, act);
	this.solver.addClause(setOfLiterals);
        return id;
    }

    

    /**
     * returns the literal that represents the logic {@code value}
     * @param value
     * @return A literal
     */
    public int constantLiteral(boolean value){
	if(value)
	    return  this.constantID;
	else 
	    return -this.constantID;}
    /**
     * returns the the id of the literal
     * @param literal
     * @return id
     */


    public int idFromLiteral(int literal){
	int id = this.isLiteralPositive(literal) ? literal: -literal;
	return id;
}

    /**
     *Prints solver stats
     */
    public void printStats(){
	PrintWriter writer = new PrintWriter(System.out);
	writer.print("c constraints: " + this.solver.nConstraints());
	
	writer.print(", variables: " + this.solver.nVars());
	int decisions = this.solver.getStat().get("decisions").intValue();
	writer.print(", decisions: " + decisions);
	int propagations = this.solver.getStat().get("propagations").intValue();
	writer.print(", propagations: " + propagations + "\n");
	writer.flush();

    }
    public String[] getStat(){
	Map<String, Number> stats = new HashMap<String, Number>();
	stats = this.solver.getStat();
	String[] statLog = new String[stats.size()];
	int i = 0;
	for(Entry<String, Number> entry: stats.entrySet()){
	    statLog[i] = entry.getKey();
	    statLog[i] += ": ";
	    statLog[i] += "" + entry.getValue().intValue();
	    i++;
}
	return statLog;
    }
    
    public void setConstantID(){
	this.constantID = this.nVars();
	try{
	    this.AddClause(new VecInt(new int[]{this.constantID}));
	}catch(ContradictionException exception){
	    throw new RuntimeException("Contradition exception while adding constant value", exception);
	};

}
}
