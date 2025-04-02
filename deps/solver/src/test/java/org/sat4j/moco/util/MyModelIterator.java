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

import java.util.Iterator;

import org.sat4j.core.VecInt;
import org.sat4j.moco.pb.PBSolver;
import org.sat4j.specs.ContradictionException;
import org.sat4j.specs.IVecInt;


    /**
     * Enumerates all models of a DIMACS pseudo boolean formula. It
     * does not preserve the state of the solver
     *
     */


public class MyModelIterator implements Iterator<boolean[]>{
	private PBSolver pbSolver;
	boolean contradiction = false;
	final IVecInt assumptions;

	public MyModelIterator(PBSolver solver, IVecInt assumptions){
	    this.pbSolver =  solver;
	    this.assumptions = assumptions;
	}
	public MyModelIterator(PBSolver solver){
	    this(solver, new VecInt(new int[0]));
	}


	public boolean hasNext(){
	    this.pbSolver.check(assumptions);
	    if(this.contradiction)
		return false;
	    return this.pbSolver.isSat();
	}
	public boolean[] next(){
	    boolean[] currentAssignment = new boolean[this.pbSolver.nVars()];
	    IVecInt notCurrent = new VecInt();
	    int litId;
	    for(int i = 0, n = currentAssignment.length; i < n ; i++){
		litId = i + 1;
		currentAssignment[i] = this.pbSolver.modelValue(litId);
		if(currentAssignment[i])
		    notCurrent.push(-litId);
		else
		    notCurrent.push(litId);

	    }
	    if(notCurrent.size() > 0)
		try {
		    this.pbSolver.AddClause(notCurrent);
		}
		catch (ContradictionException e) {
		    Log.comment(3, "Contradiction detected!");
		    this.contradiction = true;
		}
	    

	    return currentAssignment;
	}

    }
