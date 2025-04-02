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

import java.util.Collection;
import java.util.Iterator;

import org.moeaframework.core.Population;
import org.moeaframework.core.Settings;
import org.moeaframework.core.Solution;
import org.moeaframework.util.ReferenceSetMerger;
import org.sat4j.core.Vec;
import org.sat4j.moco.problem.Instance;
import org.sat4j.moco.util.Log;

/**
 * Class used for building reference sets to be used in performance indicator (e.g. inverted generational
 * distance) computation by the MOEA framework.
 * @author Miguel Terra-Neves
 */
class ReferenceSet {

    /**
     * A MOEA framework representation of the MOCO instance.
     */
    private MOCOProblem problem = null;
    
    /**
     * A MOEA framework object that merges several nondominated populations into a single reference set.
     */
    private ReferenceSetMerger merger = new ReferenceSetMerger();
    
    /**
     * The {@link #merger} object requires a different key for each nondominated population.
     * Used to store the last key generated for some population.
     */
    private int key_count = 0;
    
    /**
     * Creates an instance of an empty reference set for a given MOCO instance.
     * @param m The instance.
     */
    public ReferenceSet(Instance m) { this.problem = new MOCOProblem(m); }
    
    /**
     * Makes a fresh unused key for the {@link #merger}.
     * @return The key.
     */
    private String mkFreshKey() { return Integer.toString(this.key_count++); }
    
    /**
     * Updates the reference set with the nondominated solutions in a given result object.
     * @param r The result.
     */
    public void add(Result r) { this.merger.add(mkFreshKey(), r.getSolutions()); }

    /**
     *Ref point
     */
    private Solution refPoint = null;

    /**
     *Ideal point
     */
    private Solution idealPoint = null; 

    /**
     * Updates the reference set with the nondominated solutions in a collection of result objects.
     * @param c The result object collection.
     */
    public void addAll(Collection<Result> c) {
        for (Iterator<Result> it = c.iterator(); it.hasNext();) {
            add(it.next());
        }
    }
    
    /**
     * Retrieves the merged reference set as a MOEA framework {@link Population}.
     * @return The reference set.
     */
    private Population getMergedSet() { return this.merger.getCombinedPopulation(); }
    
    /**
     * Checks if the reference set is empty.
     * @return True if the reference set is empty, false otherwise.
     */
    public boolean isEmpty() { return getMergedSet().isEmpty(); }
    
    /**
     * Computes the ideal or nadir point of the MOCO instance based on the current reference set.
     * @param nadir True (false) if the nadir (ideal) point is to be computed.
     * @return A cost vector that corresponds to the requested point.
     */
    private Solution getPoint(boolean nadir) {
        Population p = getMergedSet();
        if (p.size() == 0) { return null; }
        Solution s = this.problem.newCostVec();
        for (int i = 0; i < this.problem.getNumberOfObjectives(); ++i) {
            double val = p.get(0).getObjective(i);
            for (int j = 1; j < p.size(); ++j) {
                double j_val = p.get(j).getObjective(i);
                val = nadir ? Math.max(val, j_val) : Math.min(val, j_val);
            }
            s.setObjective(i, val);
        }
        return s;
    }
    
    /**
     * Computes the nadir point based of the MOCO instance based on the current reference set.
     * @return The nadir point.
     */
    // TODO: multiple calls can become inefficient; use caching?
    public Solution getNadirPoint() { return getPoint(true); }
    
    /**
     * Computes the ideal point of the MOCO instance based on the current reference set.
     * @return The ideal point.
     */
    // TODO: multiple calls can become inefficient; use caching?
    
    public Solution getIdealPoint() { 
	if(this.idealPoint!=null)
	    return this.idealPoint;

	double[] idealObjInt = this.getPoint(false).getObjectives();
	double[] idealObjDouble = new double[idealObjInt.length];
	for(int iObj=0, nObj = idealObjInt.length; iObj < nObj; iObj++){
	    idealObjDouble[iObj] = idealObjInt[iObj];
	}
	this.idealPoint = new Solution(idealObjDouble);
	return this.idealPoint;
    }
    
    /**
     * Computes a reference point of the MOCO instance, for hypervolume computation, based on the current
     * reference set.
     * @return The reference point.
     */
    public Solution getRefPoint() {
	if(this.refPoint!=null)
	    return this.refPoint;

        this.refPoint = getNadirPoint();
        for (int i = 0, n = this.refPoint.getNumberOfObjectives(); i < n;++i) {
	    this.refPoint.setObjective(i, problem.maxPoint().get(i)); // using the max point as the reference point.
        }
        return this.refPoint;
    }
    
    /**
     * Retrieves the solutions in the reference set as a MOEA framework {@link Population}.
     * @return The solutions in the reference set.
     */
    // TODO: multiple calls can become inefficient; use caching?
    public Population getSolutions() {
        Population p = getMergedSet();
        if (p.size() == 0) { return p; }
        // MOEA Framework's analysis fails if the reference set has a single individual or if an objective's
        // range is empty; in order to avoid such scenarios, new solutions are injected based on the reference
        // and ideal points
        boolean single = p.size() == 1;
	Log.comment( "single  " + p.size());
        Solution nadir = getNadirPoint(), ideal = getIdealPoint(), ref = getRefPoint();
	// make sure all objectives attain a finite range inside the population. Add artificial points, if that is not the case.
        for (int i = 0; i < this.problem.getNumberOfObjectives(); ++i) {
            if (single || nadir.getObjective(i) - ideal.getObjective(i) < Settings.EPS) {
                Solution s = ref.copy();
                s.setObjective(i, s.getObjective(i) - Settings.EPS);    
                p.add(s);                                                   
		Log.comment(1, "Artificial point, objective " + i);
            }
        }
        Log.comment(1, ":ref-set-size " + p.size());
        return p;
    }
    
}
