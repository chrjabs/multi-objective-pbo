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

import org.sat4j.moco.analysis.SubResult;
import org.sat4j.moco.problem.Instance;
import org.sat4j.moco.util.Log;

/**
 * Class that implements UnsatSat, MSU3 flavoured
 * @author João Cortes
 */

public class UnsatSatMSU3 extends  UnsatSat {

    
    public UnsatSatMSU3(Instance m) {
        // Log.comment(3, "in UnsatSat constructor");
	super(m);
	this.subResult = new SubResult(this.problem);
    }




    public void saveModel(){
	Log.comment(6, "model:");
	Log.comment(6, this.prettyFormatVecInt(this.getXModel()));
	this.result.saveModel(this.solver);
    }

    public void finalizeHarvest(){}

}
