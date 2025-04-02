package org.sat4j.moco.algorithm;

import org.sat4j.moco.analysis.Result;
import org.sat4j.moco.pb.PBFactory;
import org.sat4j.moco.pb.PBSolver;
import org.sat4j.moco.problem.Instance;
import org.sat4j.moco.util.Log;
import org.sat4j.specs.ContradictionException;
import org.sat4j.specs.IVecInt;


abstract public class  algorithm{

    /**
     * An instance of a MOCO problem to be solved.
     */
    protected Instance problem = null;

    /**
     * Stores the result (e.g. nondominated solutions) of the execution of the Pareto-MCS algorithm.
     */
    protected Result result = null;

    /**
     * Stores the partial result (e.g. nondominated solutions) of the execution of the Pareto-MCS algorithm.
     */
    protected Result subResult = null;

    /**
     * Stores the PB solver to be used by the Pareto-MCS algorithm.
     */
    protected PBSolver solver = null;

    /**
     * Retrieves the result of the last call to {@link #solve()}.
     * @return The result.
     */
    public Result getResult() { return this.result; }

    /**
     *Stores any solutions that are in subResult.
     */

    public void transferSubResult() {
	if(this.subResult == null)
	    return;
	if(this.subResult != this.result){
	    for(int i = 0; i < this.subResult.nSolutions(); ++i)
		this.result.addSolutionUnsafe(this.subResult.getSolution(i));
	    this.subResult.clearPopulation();
	}
    }

    abstract public void prettyPrintVecInt(IVecInt vecInt, boolean clausing);
    abstract public String prettyFormatVecInt(IVecInt vecInt);
    abstract public void solve();



    /**
     *Adds the disjunction of setOfLiterals, and logs
     *@param setOfliterals
     */

    public boolean AddClause(IVecInt setOfLiterals){
	// Log.comment(6, "{ algorithm.AddClause");
	Log.comment(6, this.prettyFormatVecInt(setOfLiterals));
	try{
	    this.solver.AddClause(setOfLiterals);
	} catch (ContradictionException e) {
	    Log.comment(2, "contradiction when adding clause: ");
	    for(int j = 0; j < setOfLiterals.size(); ++j)
		Log.comment(2, " " + setOfLiterals.get(j) + " " );
	    return false;
	}
	// Log.comment(6, "}");
	return true;
    }


    public void printFlightRecord(){
	this.solver.printStats();
	this.printFlightRecordParticular();
    }
    public void printFlightRecordParticular(){};

    /**
     * Creates a PB oracle initialized with the MOCO's constraints.
     * @return The oracle.
     * @throws ContradictionException if the oracle detects that the
     * MOCO's constraint set is unsatisfiable.
     */
    protected PBSolver buildSolver() throws ContradictionException {
        // Log.comment(5, "in Algorithm.buildSolver");
        PBSolver solver = new PBSolver();
        solver.newVars(this.problem.nVars());
        for (int i = 0; i < this.problem.nConstrs(); ++i) {
            solver.addConstr(this.problem.getConstr(i));
        }
	solver.setConstantID();
        // Log.comment(5, "out UnsatSat.buildSolver");
        return solver;
    }

	public PBSolver getSolver() {
		return solver;
	}

	public void setSolver(PBSolver solver) {
		this.solver = solver;
	}

	public Instance getProblem() {
		return problem;
	}

	public void setProblem(Instance problem) {
		this.problem = problem;
	}

}

