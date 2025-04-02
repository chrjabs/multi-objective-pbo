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
import java.util.NoSuchElementException;
import java.util.SortedMap;
import java.util.TreeMap;
import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;

import java.util.Arrays;
import java.util.Map.Entry;
import java.util.HashMap;


import org.sat4j.core.ReadOnlyVec;
import org.sat4j.core.ReadOnlyVecInt;
import org.sat4j.core.VecInt;
import org.sat4j.moco.util.Log;
import org.sat4j.moco.util.Real;
import org.sat4j.moco.goal_delimeter.Circuit.ControlledCompIterator;
import org.sat4j.moco.goal_delimeter.Circuit.ControlledComponent;
import org.sat4j.moco.pb.PBSolver;
import org.sat4j.moco.problem.Instance;
import org.sat4j.moco.problem.Objective;
import org.sat4j.moco.problem.DigitalEnv;
import org.sat4j.moco.problem.DigitalEnv.DigitalNumber;
import org.sat4j.specs.IVecInt;

/**
 * Class with the implementation of the Selection network based encoder.
 * @author Joao O'Neill Cortes
 */

public class SelectionDelimeterMSU3 extends SelectionDelimeterT<SelectionDelimeterMSU3.ObjManager>{

    /**
     *upper limits. Only encode the circuit using weights that are
     *lighter than this.
     */
    protected Map<Integer, Integer> upperLimits;

    private int[] uncoveredMaxKD = null;
    private HashMap<Integer, Boolean> coveredLiterals = null;

    /**
     * Upper bound, exclusive
     */
    private int[] UpperBound = null;


    public SelectionDelimeterMSU3(Instance instance, PBSolver solver, boolean buildCircuit) {
	super(instance, solver, false);
	this.uncoveredMaxKD = new int[this.instance.nObjs()];
	this.UpperBound =  new int[(this.instance.nObjs())];
	this.coveredLiterals = new HashMap<Integer, Boolean>(this.solver.nVars());
	this.upperLimits = new HashMap<Integer, Integer>();
	this.initializeCoveredLiterals();
	if(buildCircuit)
	    this.buildCircuits();
    }

    public SelectionDelimeterMSU3(Instance instance, PBSolver solver, boolean buildCircuit, Map<Integer, Integer> upperLimits) {
	this(instance, solver, false);
	this.upperLimits = upperLimits;
	if(buildCircuit)
	    this.buildCircuits();

}


    private void initializeCoveredLiterals(){
	for(int iObj = 0, nObj = this.instance.nObjs();iObj < nObj; iObj++){
	    Objective ithObjective = this.instance.getObj(iObj);
	    ReadOnlyVecInt objectiveLits = ithObjective.getSubObjLits(0);
	    ReadOnlyVec<Real> objectiveCoeffs = ithObjective.getSubObjCoeffs(0);
	    int sign = 1;
	    int ithAbsoluteWeight;
	    for(int iX = 0, nX = ithObjective.getTotalLits(); iX <nX; iX ++){
		int ithX = objectiveLits.get(iX);
		ithAbsoluteWeight = objectiveCoeffs.get(iX).asIntExact();
		sign = (ithAbsoluteWeight > 0? 1 : -1);
		ithAbsoluteWeight *= sign;
		if(this.upperLimits.get(iObj) != null &&  ithAbsoluteWeight <= this.upperLimits.get(iObj))
		    this.coveredLiterals.putIfAbsent(-sign * ithX, true);
	    }
	}

    }
    // static class SDIndex extends Index{

    // 	SDIndex(int iObj, int kD){
    // 	    super(iObj, kD);
    // 	}
    // }
    ;

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
		Integer upperLimit = upperLimits.get(iObj);
		if(upperLimit != null && upperLimit <= absWeight){
		    AddClause(new VecInt(new int[]{weight > 0? -lit: lit}));
		    continue;
		}

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
	    Circuit.ControlledComponent component = circuit.getControlledComponentBase(base);
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
	    DigitalNumber.IteratorContiguous iterator = digits.iterator3();
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

	private void LexicographicOrderRecurse(DigitalNumber.IteratorContiguous iterator, IVecInt clause){
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


			// last base needed to expand the weights
			int ratioI = 0;
			int base = 1;
			int ratio = 1;
			int maxBase = 0;
			try{
			    maxBase = baseInputs.lastKey();
			}catch(NoSuchElementException e){
}

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


	@Override
	public int getIObj() {
	    return this.iObj;
	}
	

    }

	public int getUncoveredMaxKD(int iObj) {
	    return this.uncoveredMaxKD[iObj];
	}


	public void setMaxUncoveredKD(int iObj, int a) {
	    this.uncoveredMaxKD[iObj] = a;
	}



	public int getMaxUncoveredKD(int iObj) {
	    return this.uncoveredMaxKD[iObj];
	}


	@Override
	protected ObjManager[] objManagersCreator() {
	    return new ObjManager[this.getInstance().nObjs()];
	}

	@Override
	protected ObjManager objManagerCreator(int iObj) {
	    return new ObjManager(iObj);
	}


	public HashMap<Integer, Boolean> getCoveredLiterals() {
		return coveredLiterals;
	}


	public void setCoveredLiterals(HashMap<Integer, Boolean> coveredLiterals) {
		this.coveredLiterals = coveredLiterals;
	}

    /**
     *If necessary for the construction of the current assumptions,
     *initialize more of the domain of the goal delimeter
     */
    @Override
    public boolean preAssumptionsExtend(IVecInt currentExplanation){
	Log.comment(6, "explanation: ");
	Log.comment(6, this.prettyFormatVecInt(currentExplanation));
	boolean change = false;
	// Log.comment(0, "covered x variables: " + this.coveredLiterals.size());
	IVecInt currentExplanationX = new VecInt(new int[] {});
	HashMap<Integer,Boolean> objectivesToChange = new HashMap<Integer, Boolean>(this.instance.nObjs());
	for(int lit: currentExplanation.toArray()){
	    int id = this.solver.idFromLiteral(lit);
	    if(this.isX(id)){
		currentExplanationX.push(lit);
		for(int iObj = 0; iObj < this.instance.nObjs(); ++iObj){
		    if(this.instance.getObj(iObj).getSubObj(0).weightFromLit(id) != null)
			objectivesToChange.put(iObj, null);
		}
	    }
	    else
		objectivesToChange.put(this.getIObjFromY(id), null);

	}
	change = this.uncoverXs(currentExplanationX);
	for(int iObj :objectivesToChange.keySet()){
	    // Log.comment(3, "changing upperlimit " + iObj);
	    int upperKDBefore = this.getUpperKD(iObj);
 	    if(this.getUpperKD(iObj) == this.getUpperBound(iObj))
		this.generateNext(iObj,this.getUpperKD(iObj), this.getMaxUncoveredKD(iObj));
	    this.setUpperKD(iObj, this.nextKDValue(iObj, this.getUpperKD(iObj)));
	    if(this.getUpperKD(iObj) >= this.getUpperBound(iObj))
		this.generateNext(iObj, this.getUpperKD(iObj), this.getMaxUncoveredKD(iObj));
	    this.setUpperBound(iObj, this.nextKDValue(iObj, this.getUpperKD(iObj)));
	    if(this.getUpperKD(iObj)!= upperKDBefore)
		change = true;
	}
	return change;

    }
	@Override
	public IVecInt generateUpperBoundAssumptions(IVecInt explanation, boolean checkChange) {
	    IVecInt assumptions =  super.generateUpperBoundAssumptions(explanation, checkChange);
	    for(Integer x: this.coveredLiterals.keySet())
		assumptions.push(x);
	    return assumptions;
	}

    /**
     *Uncover leafs
     */
    private boolean uncoverXs(IVecInt explanationX) {
	// Log.comment(5, "{ UnsatSatMSU3.uncoverXs");
	
	for(int iObj = 0;iObj< instance.nObjs(); ++iObj){
	    Objective ithObj = instance.getObj(iObj);
	    TreeMap<Integer, Integer> ithSortedMap = this.getIthYTable(iObj);
	    SortedMap<Integer, Integer> ithSortedMapClone = new TreeMap<Integer, Integer>(this.getIthYTable(iObj));
	    for(int j = 0, n = explanationX.size(); j < n ; j++){
		int id = explanationX.get(j);
		if(id < 0) id = - id;
		Real jthCoeffReal = ithObj.getSubObj(0).weightFromLit(id);
		if(jthCoeffReal == null)
		    continue;
		int jthCoeff = jthCoeffReal.asIntExact();
		if(jthCoeff < 0) jthCoeff = -jthCoeff;
		for(int entry: ithSortedMapClone.keySet())
		    ithSortedMap.putIfAbsent(jthCoeff + entry, null);
		ithSortedMap.putIfAbsent(jthCoeff, null);
		ithSortedMapClone.putAll(ithSortedMap);
		}
	    }

	
	int lit = 0;
	for(int iLit = 0, n = explanationX.size(); iLit < n; iLit++){
	    lit = explanationX.get(iLit);
	    this.coveredLiterals.remove(lit);
	}
	this.updateAllUncoveredMaxKD();
	this.logUncoveredMaxKD();
	// Log.comment(5, "}");
	return true;
    }

private void updateAllUncoveredMaxKD(){
    for(int i = 0, n = this.instance.nObjs(); i < n; i++)
	this.updateUncoveredMaxKD(i);
}
    public void logUncoveredMaxKD(){
	String logUpperLimit = "uncovered max: ["+this.getUncoveredMaxKD(0);
	for(int iObj = 1; iObj < this.instance.nObjs(); ++iObj)
	    logUpperLimit +=", "+this.getUncoveredMaxKD(iObj) ;//+ this.instance.getObj(iObj).getMinValue())
	//..log
	
	logUpperLimit +="]";
	Log.comment(2, logUpperLimit );
    }

	public void updateUncoveredMaxKD(int iObj){
	    int a = 0;
	    Objective ithObjective = instance.getObj(iObj); // 
	    ReadOnlyVecInt objectiveLits = ithObjective.getSubObjLits(0);
	    ReadOnlyVec<Real> objectiveCoeffs = ithObjective.getSubObjCoeffs(0);
	    int sign = 1;
	    int ithAbsoluteWeight;
	    for(int iX = 0, nX = ithObjective.getTotalLits(); iX < nX; iX ++){
		int ithX = objectiveLits.get(iX);
		ithAbsoluteWeight = objectiveCoeffs.get(iX).asInt();
		sign = (ithAbsoluteWeight > 0? 1 : -1);
		ithAbsoluteWeight *= sign;
		if(coveredLiterals.get(-sign * ithX) == null)
		    a += ithAbsoluteWeight;
	    }
	    this.setMaxUncoveredKD(iObj, a);
	}
    public void setMaxUncoveredKDs(int iObj, int a){this.uncoveredMaxKD[iObj] = a;}

    private int getUpperBound(int iObj){
	return this.UpperBound[iObj];
    }

    /**
     *Sets the current upper bound of iObj to nowKD
     *@param newKD
     *@param iObj
     */
    private void setUpperBound(int iObj, int newKD){
	    this.UpperBound[iObj] = newKD;
    }

    
    @Override
    public int generateNext(int iObj, int kD, int inclusiveMax) {
	Integer next = super.nextKDValue(iObj, kD);
	if( next > inclusiveMax)
	    return kD;
	this.generateOneY(kD, iObj);
	next = super.nextKDValue(iObj, kD);
	if( next > inclusiveMax)
	    return kD;
	return next;
    }

    public int nextKDValue(int iObj, int kD) {
	if(kD < this.getUncoveredMaxKD(iObj))
	    return super.nextKDValue(iObj, kD);
	if(kD == this.getUncoveredMaxKD(iObj))
	    return kD;
	else return 0;
    }

    public void generateOneY(int kD, int iObj){
	int oldY;
	int oldKD = kD;
	int newKD = this.nextKDValue(iObj, oldKD);
	int newY = this.getY(iObj, newKD);
	if(newY != 0)
	    return;
	newY = this.getIthObjManager(iObj).LexicographicOrder(newKD);
	oldY = this.getY(iObj, oldKD);
	
	if(newKD >  this.nextKDValue(iObj, 0) && newKD > oldKD){
	    if(oldY == 0)
		oldY = this.getIthObjManager(iObj).LexicographicOrder(oldKD);
	    this.AddClause(new VecInt(new int[]{-newY, oldY}));
	    
	    int oldNextKD = this.nextKDValue(iObj, newKD);
	    if(oldNextKD > newKD){
		int oldNextY = this.getY(iObj, oldNextKD);
		if(oldNextY == 0)
		    return;
		else
		    this.AddClause(new VecInt(new int[]{-oldNextY, newY}));
	    }	
	}

    }
    public void printBasis(BufferedWriter out) throws IOException{
	for(int i = 0, n = this.instance.nObjs(); i < n; i++){
	    out.write("b " + i + " ");
	    DigitalEnv digitalEnv = this.getIthObjManager(i).getDigitalEnv();
	    for(int j = 0, m = digitalEnv.getBasesN() ; j < m; j++)
		out.write(digitalEnv.getRatio(j) + " ");
	    out.write("\n");
	}
    }

    public void printOutVariables(BufferedWriter out) throws IOException{
	for(int i = 0, n = this.instance.nObjs(); i < n; i++){
	    Circuit circuit = this.getIthObjManager(i).getCircuit();
	    ControlledCompIterator iterator =  circuit.controlledCompIterator();
	    while(iterator.hasNext()){
		out.write("l " + i + " ");
		ControlledComponent contComp = iterator.next();
		out.write(contComp.getBase() + " ");
		for(int lit: contComp.getOutputs())
		    out.write(lit + " ");
		out.write("\n");
	    }

	}
    }
    public int numberOutVars(){
	int sum = 0;
	for(int i = 0, n = this.instance.nObjs(); i < n; i++){
	    Circuit circuit = this.getIthObjManager(i).getCircuit();
	    ControlledCompIterator iterator =  circuit.controlledCompIterator();
	    while(iterator.hasNext()){
		ControlledComponent contComp = iterator.next();
		sum += contComp.getOutputsSize();
	    }
	}
	return sum;

    }

}


