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
 *   João O'Neill Cortes, INESC
 *******************************************************************************/
package org.sat4j.moco.goal_delimeter;

import java.util.Map;

import org.sat4j.moco.Params;
import org.sat4j.moco.pb.PBSolver;
import org.sat4j.moco.problem.Instance;
import org.sat4j.moco.util.Log;


public class GoalDelimeterCreator {
    Params params;
    
    public GoalDelimeterCreator(Params params){
	this.params = params;
}
    public GoalDelimeter<?> create(String encoding, Instance instance, PBSolver solver, boolean MSU3){
	GoalDelimeter<?> gd = null;
	Map<Integer, Integer[]> allRatios = this.params.getAllRatios();
	switch(encoding) {
	case "SD":
	    if (MSU3) {
		SelectionDelimeterMSU3 sd = new SelectionDelimeterMSU3(instance, solver, false, params.getUpperLimits());
		sd.initializeObjectManagers();
		for(int i = 0, n = sd.getInstance().nObjs();i<n;i++){
		    if(allRatios.get(i) != null)
			sd.getIthObjManager(i).getDigitalEnv().setRatios(allRatios.get(i));
		}
		sd.buildCircuits();
		gd = sd;
	    }
	    else{
		SelectionDelimeter gd1 = new SelectionDelimeter(instance, solver, true);
		gd1.generateY();
		gd = gd1;

	    }	    
	    break;
	case "GTE":
	    if (MSU3) 
		gd = new GenTotalEncoderMSU3(instance, solver);
	    else 
		gd = new GenTotalEncoder(instance, solver);
	    break;
	    
	default:
	    Log.comment("Don't know what encoding to use");
	    break;
	}
	return gd;
    }




}
