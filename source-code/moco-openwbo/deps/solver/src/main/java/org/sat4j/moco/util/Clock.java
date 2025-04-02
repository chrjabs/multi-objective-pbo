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
package org.sat4j.moco.util;

import org.sat4j.moco.Params;

/**
 * Singleton class for a clock used to count elapsed time.
 * @author Miguel Terra-Neves
 */
public class Clock {

    /**
     * The single instance of the clock.
     */
    private static Clock instance = null;
    
    /**
     * Retrieves the clock instance.
     * @return The clock instance.
     */
    public static Clock instance() {
        if (instance == null) { instance = new Clock(); }
        return instance;
    }
    
    /**
     * The time instant to be used as the origin when counting elapsed time.
     */
    private long start;

    /**
     * Maximum time, in seconds, allowed for the MOCO solver to run.
     * If {@code timeout} is smaller than 0, then no time limit is imposed.
     */
    private int timeout = -1;
    
    /**
     * Creates an instance of a clock.
     */
    private Clock() { }
    
    /**
     * Sets the maximum time, in seconds, allowed for the MOCO solver to run.
     * By default, no limit is imposed until this method is called.
     * @param timeout The time limit.
     */
    public void setTimeout(int timeout) { this.timeout = timeout; }
    
    /**
     * Retrieves the maximum time, in seconds, allowed for the MOCO solver to run.
     * @return The time limit.
     */
    public int getTimeout() { return this.timeout; }
    
    /**
     * Checks if the MOCO solver has a maximum time limit set.
     * @return True if a time limit is set, false otherwise.
     * @see #setTimeout(int)
     */
    public boolean hasTimeout() { return this.timeout >= 0; }
    
    /**
     * Retrieves the remaining time, in seconds, allowed for the MOCO solver to run.
     * @return The remaining time.
     */
    public int getRemaining() {
        if (!hasTimeout()) return Integer.MAX_VALUE;
        assert(timeout >= 0);
        return Math.max(this.timeout - (int)getElapsed(), 0);
    }
    
    /**
     * Checks if the MOCO solver has expired its time limit.
     * @return True if the time limit has expired, false otherwise.
     */
    public boolean timedOut() { return getRemaining() <= 0; }
    
    /**
     * Resets the clock. The current time instant becomes 0.
     */
    public void reset() { this.start = System.nanoTime(); }
    
    /**
     * Retrieves the elapsed time since the last call to {@link #reset()}. Assumes that {@link #reset()} was
     * called at least once.
     * @return The elapsed time.
     */
    public double getElapsed() { return (double)(System.nanoTime() - this.start) / 1000000000.0; }

    /**
     * Sets the time limit to the one stored in a given set of parameters.
     * @param p The parameters object.
     */
    public void updtParams(Params p) {
        if (p.hasTimeout()) {
            setTimeout(p.getTimeout());
            Log.comment(":timeout " + timeout);
        }
    }
    
}
