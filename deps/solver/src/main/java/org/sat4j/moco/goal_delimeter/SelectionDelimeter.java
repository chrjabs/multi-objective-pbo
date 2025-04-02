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
import java.util.SortedMap;
import java.util.TreeMap;

import java.util.ArrayList;

import java.util.Arrays;
import java.util.Map.Entry;
import java.util.HashMap;


import org.sat4j.core.ReadOnlyVec;
import org.sat4j.core.VecInt;
import org.sat4j.moco.util.Log;
import org.sat4j.moco.util.Real;
import org.sat4j.moco.goal_delimeter.Circuit.ControlledComponent;
import org.sat4j.moco.pb.PBSolver;
import org.sat4j.moco.problem.Instance;
import org.sat4j.moco.problem.Objective;
import org.sat4j.moco.problem.DigitalEnv;
import org.sat4j.moco.problem.DigitalEnv.DigitalNumber;
import org.sat4j.moco.problem.DigitalEnv.DigitalNumber.IteratorContiguous;
import org.sat4j.specs.IVecInt;

/**
 * Class with the implementation of the Selection network based encoder.
 * @author Joao O'Neill Cortes
 */

public class SelectionDelimeter extends SelectionDelimeterT<SelectionDelimeter.ObjManager> {
    
    public SelectionDelimeter(Instance instance, PBSolver solver, boolean buildCircuit) {
		super(instance, solver, buildCircuit);
		this.initializeYTable();
	}


	public class ObjManager implements IObjManager{
	int iObj;
	Circuit circuit;
	DigitalEnv digitalEnv;
	ObjManager(int iObj){
	    this.iObj = iObj;
	    this.digitalEnv = new DigitalEnv();
	}


	/**
	 *Generates the inputs created by the weights of the objective
	 *function iObj
	 */
	    protected SortedMap<Integer,ArrayList<Integer>> getInputsFromWeights(int iObj){
	    DigitalEnv digitalEnv = this.getDigitalEnv();
	    SortedMap<Integer,ArrayList<Integer>> baseInputs= new TreeMap<Integer, ArrayList<Integer>>();
	    HashMap<Integer, Integer> weights = getWeights(iObj);
	    List<DigitalNumber> digitsList = new ArrayList<DigitalNumber>();
	    // IVecInt digits = new VecInt(new int[]{});
	    for(Entry<Integer, Integer> entry: weights.entrySet()){
		int weight = entry.getValue();
		boolean weightSign = weight > 0;
		int lit = entry.getKey();
		int absWeight = weight > 0 ? weight: - weight;
		DigitalNumber digits = digitalEnv.toDigital(absWeight);
		digitsList.add(digits);
		// if(maxNDigits < nDigits) maxNDigits = nDigits;
		DigitalNumber.IteratorContiguous iterator = digits.iterator2();

		int ithDigit = 0;
		int ithBase = 1;
		while(iterator.hasNext())
		    {
			ithBase = iterator.currentBase();
			ithDigit = iterator.next();
			while( ithDigit > 0){
			    if(baseInputs.containsKey(ithBase))
				baseInputs.get(ithBase).add(weightSign? lit: -lit);
			    else
				baseInputs.put(ithBase, new ArrayList<Integer>(Arrays.asList(new Integer[]{weightSign? lit: -lit})));
			    ithDigit--;
			}

		    }
	    }
	    return baseInputs;

	}


	public DigitalEnv getDigitalEnv(){return this.digitalEnv;}
	public Circuit getCircuit(){return this.circuit;}
	public int digitalLiteral(int base, int value){
	    ControlledComponent component = circuit.getControlledComponentBase(base);
	    if( value <= 0 || component.getOutputsSize() == 0)
		return 0;
	    int index = unaryToIndex(value);
	    if(index > circuit.getControlledComponentBase(base).getOutputsSize())  
		index = circuit.getControlledComponentBase(base).getOutputsSize() - 1;
	    return circuit.getControlledComponentBase(base).getIthOutput(index);
	}
	/**
	 *Adds the upper bound clauses  that enforce the inclusive
	 *upper limit upperLimit. Returns the blocking variables.
	 *@return blocking variables
	 *@param upperLimit inclusive upper limit
	 *@param iObj the objective index
	 */


	public int LexicographicOrder(int upperLimit){
	    DigitalNumber digits = digitalEnv.toDigital(upperLimit);
	    IteratorContiguous iterator = digits.iterator3();
	    int activator = getSolver().getFreshVar();
	    IVecInt clause = new VecInt(new int[]{activator});
	    SDIndex sDIndex = new SDIndex(iObj, upperLimit);
	    librarian.putIndex(activator, sDIndex);
	    setY(iObj, upperLimit, activator);
	    Log.comment(6, "Lexicographic order");
	    this.LexicographicOrderRecurse(iterator, clause);
	    AddClause(clause);
	    return activator;
	}

	private void LexicographicOrderRecurse(IteratorContiguous iterator, IVecInt clause){
	    int base = iterator.currentBase();
	    int ratio = digitalEnv.getRatio(iterator.getIBase());
	    int digit = iterator.next();
	    int lit = 0;
	    if(digit + 1 < ratio){
		lit = digitalLiteral(base, digit + 1);
		if(lit != 0){
		    clause.push(-lit);
		    AddClause(clause);
		    clause.pop();
		}
	    }
	    if(digit > 0){
		lit = digitalLiteral(base, digit);
		if(lit != 0)
		    clause.push(-lit);
	    }
	    if(iterator.hasNext())
		this.LexicographicOrderRecurse(iterator, clause);
	}



	@Override
	public int getIObj() {
	    return this.iObj;
	}


	    /**
	     *Optimize the ratios. Notice that this.digitalEnv is
	     *reset by this.
	     */
	    public boolean optimizeRatios(int maxW, int maxSum){
		DigitalEnv digitalEnv = this.getDigitalEnv();
		int t;
		DigitalNumber digits =  digitalEnv.toDigital(maxW);
		t = this.getDigitalEnv().getBaseI(digits.getMSBase()) - 1 ;
		//TODO: return false if ratios stay the same
		int x = maxSum;
		for(int i = 0, n = t ; i < n; i++ ){
		    x -= digitalEnv.getBase(i)* (digitalEnv.getRatio(i) - 1);
		}

		//only change ratios if the t'th digit pushes at least
		//one carry
		if(x >=1){
		    Integer[] ratios = new Integer[t + 1];
		    int lastRatio = 0;
		    int lastBase = digitalEnv.getBase(t);
		    while(x > 0){
			for(int j = 0; j < lastBase; j++)
			    x--;
			lastRatio++;
		    }
		    lastRatio++;
		    for(int i = 0; i < t ; i++)
			ratios[i] = digitalEnv.getRatio(i);
		    ratios[t] = lastRatio;
		    this.digitalEnv = new DigitalEnv();
		    this.digitalEnv.setRatios(ratios);
		    getDigitalEnv().toDigital(maxSum);
		}
		return true;		
	    }

	@Override
	public void buildMyself() {
    	    this.circuit = new Circuit(getSolver()){
		    public void buildCircuit(){
			SortedMap<Integer, ArrayList<Integer>> baseInputs = getInputsFromWeights(iObj);
			ArrayList<Integer> inputs = new ArrayList<Integer>();

			int maxValue = getInstance().getObj(getIObj()).getWeightDiff();
			ReadOnlyVec<Real> coeffs = getInstance().getObj(getIObj()).getSubObj(0).getCoeffs();
			int maxCoeff = 0;
			for(int i = 0, n = coeffs.size(); i < n; i++ ){
			    int currentCoeff = coeffs.get(i).asIntExact();
			    if(currentCoeff < 0) currentCoeff = -currentCoeff;
			    if(currentCoeff > maxCoeff)
				maxCoeff = currentCoeff;
			}
			optimizeRatios(maxCoeff, maxValue);
			// to recover digitalEnv setup,
			baseInputs = getInputsFromWeights(iObj);
			//to make sure digitalEnv.basesN is correct,
			getDigitalEnv().toDigital(maxValue);
			// last base needed to expand the weights
			int ratioI = 0;
			int base = 1;
			int ratio = 1;
			int maxBase = baseInputs.lastKey();
			List<Integer> carryBits = null;
			int basesN = 1;
			do{
			    ratio = digitalEnv.getRatio(ratioI++);
			    inputs.clear();
			    ArrayList<Integer> inputsWeights = baseInputs.get(base);
			    if(carryBits != null)
				inputs.addAll(carryBits);		    
			    if(inputsWeights!=null)
				inputs.addAll(inputsWeights);
			    if(base <= maxBase || inputs.size() != 0){
				if(base > maxBase)
				    digitalEnv.setBasesN(basesN);
				carryBits =
				    buildControlledComponent(inputs.toArray(new Integer[0]), base, ratio);
			    } else{break;}
			    base *=ratio;
			    basesN++;
			}while(true);
		    

		    }

		    /**
		     *range is the exclusive upper value the unary output may represent.
		     */
		    public List<Integer> buildControlledComponent(Integer[] inputs, int base, int modN){
			DigitComponent digitComp = new DigitComponent(inputs, modN);
			digitComp.constitutiveClause();
			new ControlledComponent(base, digitComp);
			return digitComp.getCarryBits(modN);
		    }
		    public int getFreshVar1(){return getSolver().getFreshVar();}
		    public boolean AddClause1(IVecInt setOfLiterals){return AddClause(setOfLiterals);}
		};

	    this.circuit.buildCircuit();
	}



    }

	@Override
	protected ObjManager[] objManagersCreator() {
	    return new ObjManager[this.getInstance().nObjs()];
	}

	@Override
	protected ObjManager objManagerCreator(int iObj) {
	    return new ObjManager(iObj);
	}


    public void generateY(){
	for(int iObj = 0, nObj = instance.nObjs() ;iObj< nObj; ++iObj){
	    Objective ithObjective = this.getInstance().getObj(iObj);
	    int oldActivator;
	    int activator = 0;
	    int max = ithObjective.getWeightDiff();

	    int oldKD = this.nextKDValue(iObj, 0);
	    int kD = this.nextKDValue(iObj, oldKD);
	    oldActivator = this.getIthObjManager(iObj).LexicographicOrder(oldKD);
	    while(kD > oldKD){
		activator = this.getIthObjManager(iObj).LexicographicOrder(kD);
		if(kD > 1){
		    Log.comment(6, "sequential clause");
		    this.AddClause(new VecInt(new int[]{-activator, oldActivator}));
		}
		oldKD = kD;
		kD = this.nextKDValue(iObj, oldKD);		
		oldActivator = activator;
	    }
	}
    }

	@Override
	protected void initializeYTable() {
	    super.initializeYTable();
	    for(int iObj = 0;iObj< instance.nObjs(); ++iObj){
		Objective ithObj = instance.getObj(iObj);
		ReadOnlyVec<Real> objectiveCoeffs = ithObj.getSubObjCoeffs(0);	    
		TreeMap<Integer, Integer> ithSortedMap = this.getIthYTable(iObj);
		SortedMap<Integer, Integer> ithSortedMapClone = new TreeMap<Integer, Integer>(ithSortedMap);
		for(int j = 0, n = objectiveCoeffs.size(); j < n ; j++){
		    int jthCoeff = objectiveCoeffs.get(j).asIntExact();
		    if(jthCoeff < 0) jthCoeff = -jthCoeff;
		    for(int entry: ithSortedMapClone.keySet())
			ithSortedMap.put(jthCoeff + entry, null);
		    ithSortedMap.put(jthCoeff, null);
		    ithSortedMapClone.putAll(ithSortedMap);
		}
	    }

	}


}




