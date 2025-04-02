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

import java.io.BufferedReader;
import java.io.IOException;
import java.io.Reader;

import org.moeaframework.core.Solution;
import org.sat4j.moco.problem.Instance;

/**
 * Class for parsing the result of a MOCO solver execution from a log file.
 * @author Miguel Terra-Neves
 */
class OutputReader {

    /**
     * The reader object from which to the parse the result.
     */
    private BufferedReader reader = null;
    
    /**
     * The MOCO instance for which the log was generated.
     */
    private Instance moco = null;
    
    /**
     * A MOEA framework representation of the MOCO instance.
     */
    private MOCOProblem problem = null;
    
    /**
     * Creates a log output reader object that parses the result of a solver's execution for a given MOCO
     * instance.
     * @param reader The reader object from which to parse the result.
     * @param m The instance for which the log was generated.
     */
    public OutputReader(Reader reader, Instance m) {
        this.reader = new BufferedReader(reader);
        this.moco = m;
        this.problem = new MOCOProblem(m);
    }
    
    /**
     * Reads a MOCO solver execution result for the MOCO instance and from the {@link Reader} object provided
     * in {@link #OutputReader(Reader, Instance)}.
     * @return The solver execution's result.
     * @throws IOException if an error occurs reading the log.
     */
    public Result readResult() throws IOException {
        Result r = new Result(this.moco);
        boolean add_mode = false;
        for (String line = this.reader.readLine(); line != null; line = this.reader.readLine()) {
            line = line.trim();
            if (line.isEmpty() || (!add_mode && !line.startsWith("s "))) continue;
            if (line.startsWith("s ")) { add_mode = true; 
		if(parseState(line))
		    r.setParetoFrontFound();
		continue;
	    }
            if (add_mode && line.startsWith("o ")) {
                r.addSolution(parseCostVec(line));
            }
        }
        return r;
    }
    
    /**
     * Parses a cost vector from a given entry of cost values in the log.
     * @param line The log entry.
     * @return The cost vector as a MOEA framework {@link Solution} object.
     */
    private Solution parseCostVec(String line) {
        assert(line.startsWith("o "));
        String[] tokens = line.trim().split("\\s+");
        Solution s = this.problem.newCostVec();
        for (int i = 1; i < tokens.length; ++i) {
            String tok = tokens[i];
            s.setObjective(i-1, Double.parseDouble(tok));
        }
        return s;
    }
    /**
     * Parses the final state
     */
    private boolean parseState(String line) {
        assert(line.startsWith("s "));
        String[] tokens = line.trim().split("\\s+");
	if(tokens[1].equals("OPTIMUM"))
	    return true;
	else
	    return false;
    }
    
    /**
     * Frees the I/O resources associated with the log output reader.
     * @throws IOException if an error occurs closing the reader.
     */
    public void close() throws IOException { this.reader.close(); }
    
}
