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

import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.ArrayList;

import java.util.HashMap;


import org.sat4j.core.ReadOnlyVec;
import org.sat4j.moco.util.Real;
import org.sat4j.moco.pb.PBSolver;
import org.sat4j.moco.problem.Instance;
import org.sat4j.moco.problem.Objective;
import org.sat4j.specs.IVecInt;

/**
 * Class with the implementation of the Selection network based encoder.
 * @author Joao O'Neill Cortes
 */

abstract public class SelectionDelimeterT<PObjManager extends IObjManager> extends GoalDelimeter<SelectionDelimeterT.SDIndex> {

    private PObjManager[] objManagers;
    private ArrayList<TreeMap<Integer, Integer>> yTable = null;
    
    public SelectionDelimeterT(Instance instance, PBSolver solver, boolean buildCircuit) {
	super(instance, solver);
	this.objManagers =  objManagersCreator();
	this.initializeObjectManagers();
	this.initializeYTable();
	if(buildCircuit){
	    this.buildCircuits();
	}
	// Log.comment(5, "}");
    }

    


    static class SDIndex extends Index{

	SDIndex(int iObj, int kD){
	    super(iObj, kD);
	}
    }

    abstract protected PObjManager[] objManagersCreator();
    abstract protected PObjManager objManagerCreator(int iObj);

public void initializeObjectManagers(){
	for(int iObj = 0, nObj = instance.nObjs() ;iObj< nObj; ++iObj){
	    this.objManagers[iObj] = this.objManagerCreator(iObj);
	}
    }

public void buildCircuits(){
	for(int iObj = 0, nObj = instance.nObjs() ;iObj< nObj; ++iObj){
	    objManagers[iObj].buildMyself();

	}
    }

    /**
     * Initialize the container of the Y variables
     */
    protected void initializeYTable(){
	this.yTable = new ArrayList<TreeMap<Integer, Integer>>();

	for(int iObj = 0;iObj< instance.nObjs(); ++iObj){
	    TreeMap<Integer, Integer> ithSortedMap = new TreeMap<Integer, Integer>();
	    this.yTable.add(ithSortedMap);
	}
    }

    public PObjManager getIthObjManager(int i){return this.objManagers[i];}

;




    public Integer[] concatenate(Integer[][] seq){
	List<Integer> result = new ArrayList<Integer>();
	for(Integer[] array: seq)
	    for(Integer value: array)
		result.add(value);
	return result.toArray(new Integer[0]);
    }    



    protected HashMap<Integer, Integer> getWeights(int iObj){
	HashMap<Integer, Integer> result = new HashMap<Integer, Integer>();
	Objective ithObjective = this.instance.getObj(iObj);
	ReadOnlyVec<Real> objectiveCoeffs = ithObjective.getSubObjCoeffs(0);
	IVecInt literals = ithObjective.getSubObjLits(0);
	for(int i = 0, n = objectiveCoeffs.size(); i < n; i++)
	    {
		int weight = objectiveCoeffs.get(i).asIntExact();
		int lit = literals.get(i);
		result.put(lit, weight);
	    }
	return result;

    }


    public boolean isY(int id){
	SDIndex index = librarian.getIndex(id);
	if(index == null) return false;
	return true;
    }
    public int getIObjFromY(int id){
	SDIndex index = this.librarian.getIndex(id);
	if(index!=null)
	    return index.getIObj();
	return 0;};

    public int getKDFromY(int id){
	SDIndex index = this.librarian.getIndex(id);
	if(index!=null)
	    return index.getKD();
	return 0;};

    public void  setY (int iObj, int kD, int id){
	this.yTable.get(iObj).put(kD, id);
}
    public int getY(int iObj, int kD){
	
	Integer y = this.yTable.get(iObj).get(kD);
	if(y == null)
	    return 0;
	return y;
    }


    public int unaryToIndex(int kD){
	return kD  - 1;

    }

	@Override
	public int getCurrentKD(int iObj) {
	    return this.upperKD[iObj];
	}


	@Override
	public int nextKDValue(int iObj, int kD) {
	    Integer value = this.yTable.get(iObj).ceilingKey(kD + 1);
	    if(value == null)
		return kD;
	    return value;
	}

    public TreeMap<Integer, Integer> getIthYTable(int iObj) {
	return yTable.get(iObj);
    }
}


