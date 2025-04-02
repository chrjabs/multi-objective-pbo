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
package org.sat4j.moco.algorithm;

import java.util.Arrays;

import org.moeaframework.core.PRNG;
import org.sat4j.core.ReadOnlyVec;
import org.sat4j.core.ReadOnlyVecInt;
import org.sat4j.core.Vec;
import org.sat4j.core.VecInt;
import org.sat4j.moco.Params;
import org.sat4j.moco.analysis.Result;
import org.sat4j.moco.mcs.IModelListener;
import org.sat4j.moco.mcs.MCSExtractor;
import org.sat4j.moco.pb.PBFactory;
import org.sat4j.moco.pb.PBSolver;
import org.sat4j.moco.problem.Instance;
import org.sat4j.moco.problem.Objective;
import org.sat4j.moco.util.Log;
import org.sat4j.moco.util.Real;
import org.sat4j.specs.ContradictionException;
import org.sat4j.specs.IVec;
import org.sat4j.specs.IVecInt;

/**
 * Class that implements the Pareto-MCS based algorithm for MOCO, proposed in:<br>
 *      Terra-Neves, M., Lynce, I., &amp; Manquinho, V. (2017, August).
 *      Introducing Pareto minimal correction subsets. In International Conference on Theory and Applications
 *      of Satisfiability Testing (pp. 195-211). Springer, Cham.<br>
 * Includes MOCO stratification, proposed in:<br>
 *      Terra-Neves, M., Lynce, I., &amp; Manquinho, V. M. (2018).
 *      Stratification for Constraint-Based Multi-Objective Combinatorial Optimization. In IJCAI
 *      (pp. 1376-1382).
 * @author Miguel Terra-Neves
 */
public class ParetoMCS extends algorithm {
    
    
    /**
     * Stores the MCS extractor to be used by the Pareto-MCS algorithm.
     */
    private MCSExtractor extractor = null;
    
    /**
     * Creates an instance of a MOCO solver, for a given instance, that applies the Pareto-MCS algorithm.
     * @param m The MOCO instance.
     */
    public ParetoMCS(Instance m) {
        this.problem = m;
        this.result = new Result(m);
        try {
            this.solver = buildSolver();
        }
        catch (ContradictionException e) {
            Log.comment(3, "Contradiction in ParetoMCS.buildSolver");
            this.result.setParetoFrontFound();
            return;
        }
        this.extractor = new MCSExtractor(this.solver);
        this.extractor.setModelListener(new IModelListener() {
            public void onModel(PBSolver s) {
                result.saveModel(s);
            }
        });
    }
    
    /**
     * Retrieves the result of the last call to {@link #solve()}.
     * @return The result.
     */
    public Result getResult() { return this.result; }
    
    /**
     * Applies the Pareto-MCS algorithm to the MOCO instance provided in {@link #ParetoMCS(Instance)}.
     * If the instance has already been solved, nothing happens.
     */
    public void solve() {
        if (this.result.isParetoFront()) {
            Log.comment(1, "ParetoMCS.solve called on already solved instance");
            return;
        }
        Log.comment(3, "{ ParetoMCS.solve");
        int nmcs = 0;
        initUndefFmls();
        IVec<IVecInt> undef_fmls = buildUndefFmls();
        extractor.extract(undef_fmls);
        while (extractor.isSolved() && extractor.foundMCS()) {
            ++nmcs;
            if (extractor.getMSS().isEmpty()) { break; }    // if MSS is empty, then only 1 MCS exists
            try {
                solver.addConstr(PBFactory.instance().mkClause(extractor.getMCS()));
                if (this.stratify) { undef_fmls = buildUndefFmls(); }
                extractor.extract(undef_fmls);
            }
            catch (ContradictionException e) {
                Log.comment(3, "contradiction blocking MCS");
                break;
            }
        }
        Log.comment(1, ":mcs-found " + nmcs);
        if (extractor.isSolved()) {
            this.result.setParetoFrontFound();
        }
        else {
            Log.comment(1, "MCS extraction timeout");
        }
        Log.comment(3, "out ParetoMCS.solve");
    }
    
    
    /**
     * Builds a partition sequence of the literals in the objective functions to be used for stratified
     * MCS extraction.
     * If stratification is disabled, a single partition is returned.
     * @return The objective literals partition sequence.
     */
    private IVec<IVecInt> buildUndefFmls() {
        Log.comment(3, "{ ParetoMCS.buildUndefFmls");
        IVec<IVecInt> fmls = new Vec<IVecInt>();
        IVec<IVec<IVecInt>> p_stacks = new Vec<IVec<IVecInt>>(this.undef_parts.size());
        for (int i = 0; i < this.undef_parts.size(); ++i) {
            IVec<IVecInt> parts = this.undef_parts.get(i);
            IVec<IVecInt> p_stack = new Vec<IVecInt>(parts.size());
            for (int j = parts.size()-1; j >= 0; --j) {
                p_stack.unsafePush(parts.get(j));
            }
            p_stacks.unsafePush(p_stack);
        }
        while (!p_stacks.isEmpty()) {
            int rand_i = PRNG.nextInt(p_stacks.size());
            IVec<IVecInt> rand_stack = p_stacks.get(rand_i);
            fmls.push(new ReadOnlyVecInt(rand_stack.last()));
            rand_stack.pop();
            if (rand_stack.isEmpty()) {
                p_stacks.set(rand_i, p_stacks.last());
                p_stacks.pop();
            }
        }
        Log.comment(3, "out ParetoMCS.buildUndefFmls");
        return fmls;
    }
    
    /**
     * Representation of weighted literals.
     * Used to store and sort literals based on their coefficients in some objective function.
     * @author Miguel Terra-Neves
     */
    private class WeightedLit implements Comparable<WeightedLit> {

        /**
         * Stores the literal.
         */
        private int lit = 0;
        
        /**
         * Stores the weight.
         */
        private Real weight = Real.ZERO;
        
        /**
         * Creates an instance of a weighted literal with a given weight.
         * @param l The literal.
         * @param w The weight.
         */
        WeightedLit(int l, Real w) {
            this.lit = l;
            this.weight = w;
        }
        
        /**
         * Retrieves the literal part of the weighted literal.
         * @return The literal.
         */
        int getLit() { return this.lit; }
        
        /**
         * Retrieves the weight part of the weighted literal.
         * @return The weight.
         */
        Real getWeight() { return this.weight; }
        
        /**
         * Compares the weighted literal to another weighted literal.
         * The weighted literal order is entailed by their weights.
         * @param other The other weighted literal.
         * @return An integer smaller than 0 if this literal's weight is smaller than {@code other}'s, 0 if
         * the weight are equal, an integer greater than 0 if this literal's weight is larger than
         * {@code other}'s.
         */
        public int compareTo(WeightedLit other) {
            return getWeight().compareTo(other.getWeight());
        }
        
    }
    
    /**
     * Stores each objective's individual literal partition sequence.
     */
    private IVec<IVec<IVecInt>> undef_parts = null;
    
    /**
     * Boolean indicating if stratification is to be used.
     */
    private boolean stratify = false;
    
    /**
     * Stratification parameter that controls the literal-weight ratio used in the objective literals
     * partitioning process.
     */
    private double lwr = 15.0;
    
    /**
     * Initializes the objective literal partition sequences.
     * If stratification is disabled, a single partition is created with all objective literals for all
     * objective functions.
     * If stratification is enabled, an individual partition sequence is built for each objective function,
     * to later be mixed in during the search process.
     * @see #buildUndefFmls()
     */
    private void initUndefFmls() {
        Log.comment(3, "{ ParetoMCS.initUndefFmls");
        this.undef_parts = new Vec<IVec<IVecInt>>();
        for (int i = 0; i < this.problem.nObjs(); ++i) {
            Objective o = this.problem.getObj(i);
            if (this.stratify) {
                this.undef_parts.push(partition(o));
            }
            else if (this.undef_parts.isEmpty()) {
                this.undef_parts.push(new Vec<IVecInt>());
                this.undef_parts.get(0).push(singlePartition(o));
            }
            else {
                singlePartition(o).copyTo(this.undef_parts.get(0).get(0));
            }
        }
        logPartitions();
        Log.comment(3, "out ParetoMCS.initUndefFmls");
    }
    
    /**
     * Logs the number of partitions for each objective function and partition sizes.
     */
    private void logPartitions() {
        for (int i = 0; i < this.undef_parts.size(); ++i) {
            IVec<IVecInt> obj_parts = this.undef_parts.get(i);
            Log.comment(1, ":obj-idx " + i + " :partitions " + obj_parts.size());
            for (int j = 0; j < obj_parts.size(); ++j) {
                Log.comment(1, ":part-idx " + j + " :part-size " + obj_parts.get(j).size());
            }
        }
    }
    
    /**
     * Builds a literal partition sequence for a given objective.
     * @param o The objective.
     * @return A partition sequence for objective {@code o}.
     */
    private IVec<IVecInt> partition(Objective o) {
        IVec<IVecInt> parts = new Vec<IVecInt>();
        IVec<WeightedLit> w_lits = getWeightedLits(o);
        WeightedLit[] w_lits_array = new WeightedLit[w_lits.size()];
        w_lits.copyTo(w_lits_array);
        Arrays.sort(w_lits_array);
        IVecInt part = new VecInt();
        int w_count = 0;
        for (int i = w_lits_array.length-1; i >= 0; --i) {
            if (    i < w_lits_array.length-1 &&
                    !w_lits_array[i].getWeight().equals(w_lits_array[i+1].getWeight()) &&
                    (double)part.size() / w_count > this.lwr) {
                parts.push(part);
                part = new VecInt();
                w_count = 0;
            }
            part.push(-w_lits_array[i].getLit());
            if (w_count == 0 || !w_lits_array[i].getWeight().equals(w_lits_array[i+1].getWeight())) {
                w_count++;
            }
        }
        assert(!part.isEmpty());
        parts.push(part);
        return parts;
    }
    
    /**
     * Builds a single partition for a given objective containing all of the objective's literals.
     * @param o The objective.
     * @return A partition with all of objective {@code o}'s literals.
     */
    private IVecInt singlePartition(Objective o) {
        IVec<WeightedLit> w_lits = getWeightedLits(o);
        IVecInt part = new VecInt(w_lits.size());
        for (int i = 0; i < w_lits.size(); ++i) {
            part.unsafePush(-w_lits.get(i).getLit());
        }
        return part;
    }
    
    /**
     * Retrieves the literals and respective coefficients in an objective function as a vector of weighted
     * literals.
     * @param o The objective.
     * @return The objective's literals and coefficients as weighted literals.
     */
    private IVec<WeightedLit> getWeightedLits(Objective o) {
        IVec<WeightedLit> w_lits = new Vec<WeightedLit>();
        for (int i = 0; i < o.nSubObj(); ++i) {
            ReadOnlyVecInt lits = o.getSubObjLits(i);
            ReadOnlyVec<Real> coeffs = o.getSubObjCoeffs(i);
            for (int j = 0; j < lits.size(); ++j) {
                int lit = lits.get(j);
                Real coeff = coeffs.get(j);
                if (coeff.isPositive()) {
                    w_lits.push(new WeightedLit(lit, coeff));
                }
                else if (coeff.isNegative()) {
                    w_lits.push(new WeightedLit(-lit, coeff.negate()));
                }
                else {
                    Log.comment(2, "0 coefficient ignored");
                }
            }
        }
        return w_lits;
    }
    
    /**
     * Sets the algorithm configuration to the one stored in a given set of parameters.
     * @param p The parameters object.
     */
    public void updtParams(Params p) {
        this.stratify = p.getStratify();
        Log.comment(":stratify " + this.stratify);
        if (this.stratify) {
            this.lwr = p.getLWR();
            Log.comment(":lwr " + this.lwr);
        }
        this.extractor.updtParams(p);
    }
    /**
     *TODO
     */
    public void prettyPrintVecInt(IVecInt vecInt, boolean clausing){};
    public String prettyFormatVecInt(IVecInt literals){return "";}
    public void printFlightRecordParticular(){};


}
