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

import java.io.PrintStream;

import org.sat4j.moco.Params;

/**
 * Class with several static methods for logging.
 * @author Miguel Terra-Neves
 */
public class Log {

    /**
     * comment indent 
     */
    private int commentIndent = 0;
    /**
     * The output stream used by the logger ({@link System#out} by default).
     */
    private static PrintStream out = System.out;
    
    /**
     * The verbosity level of the logger.
     * Any value lower than 0 disables the logger.
     */
    private static int verb_lvl = -1; // all output disabled by default
    
    /**
     * Boolean indicating if assignment logging should be suppressed.
     * If true, then assignments are not logged.
     */
    private static boolean suppress_assign = false;
    public Log(int verbosity){
	setVerbosity(verbosity);
}
    
    /**
     * Sets the verbosity level of the logger.
     * Any value lower than 0 disables the logger.
     * @param lvl The verbosity level.
     */
    public static void setVerbosity(int lvl) { Log.verb_lvl = lvl; }
    
    /**
     * Sets if assignment logging should be suppressed.
     * @param sa True if assignment logging is to be suppressed, false otherwise.
     */
    public static void setSuppressAssignments(boolean sa) { Log.suppress_assign = sa; }
    
    /**
     * Sets the output stream of the logger.
     * @param ps The output stream.
     */
    public static void setStream(PrintStream ps) { Log.out = ps; }
    
    /**
     * Writes an entry in the log if the verbosity level is greater than or equal to a given level.
     * @param lvl The minimum verbosity level for the write to occur.
     * @param s The log entry.
     */
    private static void write(int lvl, String s) { if (lvl <= Log.verb_lvl) Log.out.println(s); }
    
    /**
     * Writes a comment in the log if the verbosity level is greater than or equal to a given level.
     * @param lvl The minimum verbosity level for the write to occur.
     * @param s The comment's content.
     */
    public static void comment(int lvl, String s  ) { write(lvl, "c " + s); }
    
    public static void clausing(String clause) { write(4, "p " + clause); }

    /**
     * Writes a comment in the log if logging is enabled.
     * @param s The comment's content.
     */
    public static void comment(String s) { comment(0, s); }
    
    /**
     * Writes an entry of cost values if logging is enabled.
     * @param costs An array of costs.
     */
    public static void costs(double[] costs) {
        String s = "o";
        for (int i = 0; i < costs.length; ++i) {
            s += " " + costs[i];
        }
        write(0, s);
    }
    
    /**
     * Logs the total time, and the number of solutions found
     */
    public static void numberOfSolutions(int n){write(0, "n " + n);}
    /**
     * Writes an 'optimum found' entry if logging is enabled.
     */
    public static void optimum() { write(0, "s OPTIMUM"); }
    
    /**
     * Writes an 'unsatisfiable instance' entry if logging is enabled.
     */
    public static void unsat() { write(0, "s UNSATISFIABLE"); }
    
    /**
     * Writes a 'satisfiable instance' entry if logging is enabled.
     */
    public static void sat() { write(0, "s SATISFIABLE"); }
    
    /**
     * Writes an 'unknown' entry if logging is enabled.
     */
    public static void unknown() { write(0, "s UNKNOWN"); }
    
    /**
     * Writes an assignment entry if logging is enabled and assignment logging is not suppressed.
     * @param s An assignment.
     * @see #setSuppressAssignments(boolean)
     */
    public static void assignment(boolean[] s) {
        if (Log.suppress_assign) return;
        String s_str = "v";
        for (int i = 0; i < s.length; ++i) {
            s_str += (s[i] ? " x" : " -x") + (i+1);
        }
        Log.write(0, s_str);
    }
    
    /**
     * Sets the logger's verbosity level and assignment suppression configuration to the ones stored in a
     * given set of parameters.
     * @param p The parameters object.
     */
    public static void updtParams(Params p) {
        setVerbosity(p.getVerbosity());
        comment(":verbosity " + Log.verb_lvl);
        setSuppressAssignments(p.getSuppressAssignments());
    }
    
}
