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
import org.sat4j.moco.goal_delimeter.GoalDelimeter;
import org.sat4j.specs.ContradictionException;
import org.sat4j.specs.IVecInt;


    /**
     * Enumerates all models of a DIMACS pseudo boolean formula. It
     * does not preserve the state of the solver
     *
     */

  public class General{
      static public String FormatArrayWithValues(Integer[] literals, PBSolver solver, boolean print){
	String result = "";
	for(int j = 0, n = literals.length; j < n ; ++j)
	    result += literals[j] + " " +solver.modelValue(literals[j]) + "|\n";
	if(print)
	    System.out.println(result);
	return result;
    }
      static public String FormatArrayWithValues(int[] literals, PBSolver solver, boolean print){
	String result = "";
	for(int j = 0, n = literals.length; j < n ; ++j)
	    result += literals[j] + " " +solver.modelValue(literals[j]) + "|\n";
	if(print)
	    System.out.println(result);
	return result;
    }

      static public String FormatModel(boolean[] model, boolean print){
	  String result = "";
	  for(int i = 0, n = model.length; i < n; i++){
	      result += i + " " + (model[i]? 1 : 0) + "\n";
}
	  if(print)
	      System.out.println(result);
	  return result;
      }
    }
