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

import org.sat4j.moco.util.Log;

import java.util.Iterator;
import java.util.PriorityQueue;
import java.util.ArrayList;
import java.util.TreeMap;
import java.util.Collection;
import java.util.Hashtable;

import org.sat4j.core.ReadOnlyVec;
import org.sat4j.core.ReadOnlyVecInt;
import org.sat4j.core.VecInt;
import org.sat4j.moco.util.Real;
import org.sat4j.moco.pb.PBSolver;
import org.sat4j.moco.problem.Instance;
import org.sat4j.moco.problem.Objective;
import org.sat4j.specs.IVecInt;

/**
 * Class with the implementation of the generalized totalizor encoder.
 * Notice that the differential is both a value and an index starting at 0
 * @author Joao O'Neill Cortes
 */


public class GenTotalEncoder extends GoalDelimeter<GenTotalEncoder.GTEIndex> {

    /**
     *The inverse index map for the partial sum variables. For each ID, a
     *value that is an array vector with the value of the goal and the
     *value of the sum
     */

    class SumTree {

	private int iObj = -1;
	/**
	 *must old the desired upperLimit, at any time. Desired is
	 *purposefully
	 */
	private int upperLimit = 0;
	/**
	 *Must old the last but effective upperLimit.
	 */
	private int olderUpperLimit = 0;

	/**
	 *Maximal possible value
	 */
	private int maxUpperLimit = 0;

	/**
	 *The root of the SumTree.
	 */
	private Node parent = null;
	/**
	 *List of nodes. Simple array.
	 */
	private ArrayList<Node> nodes = new ArrayList<Node>();

	/**
	 *List of unlinked nodes. Discardable.
	 */
	private PriorityQueue<Node> unlinkedNodes = new PriorityQueue<Node>((a,b) -> a.nodeSum - b.nodeSum);

	/**
	 *Node of a SumTree.
	 */
	class Node {

	    /**
	     *Container of the variables associated with the SumTree
	     */
	    class NodeVars{

		/**
		 * Ordered map.key is the nodeSum, value is the variable.
		 */
		private TreeMap<Integer, NodeVar> containerAll = null;
		 
		/**
		 *variable representation. 
		 */
		class NodeVar {
		    private int kD;
		    private Integer id= null;

		    private boolean iAmFresh = false;
		    public NodeVar(int kD){
			this.setKD(kD);
		    }

		    public Integer getId(){return this.id;}
		    public int getKD(){return this.kD;}
		    public void setKD(int newKD){
			this.kD = newKD;
		    }
		    private void cutValue(){
			if(this.kD > upperLimit )
			    this.setKD(upperLimit);
			
		    }
		    /**
		     *Return the ID of a freshly created auxiliar variable
		     */
		    
		    protected void setFreshId(){
			assert this.id == null;
			solver.newVar();
			this.id = solver.nVars();
		    }
                   protected boolean newValidVariable(){
		       if(olderUpperLimit < this.kD)
			   if(this.kD <= maxUpperLimit)
			       return true;
		       return false;
		    }

		    public boolean iAmFresh(){
			return iAmFresh;
		    }
		}
		public NodeVars(){
		    this.containerAll = new TreeMap<Integer, NodeVar>();
		}

		public NodeVar add(int kD, int id, boolean cutKD, boolean clausing){
		    cutKD = false;
		    NodeVar newNodeVar =  new NodeVar(kD);
		    if(cutKD){ newNodeVar.cutValue(); kD = newNodeVar.getKD();}
		    if(!clausing || newNodeVar.newValidVariable()){
		    this.containerAll.put(kD, newNodeVar);
		    if(id == 0) 
			newNodeVar.setFreshId();
		    else    newNodeVar.id = id;
		    librarian.putIndex(newNodeVar.getId(), new GTEIndex(iObj, kD, nodeName));
		    // Log.comment(6, "var " + prettyFormatVariable(newNodeVar.getId()));
		    if(kD == 0)
			AddClause( new VecInt(new int[] {newNodeVar.getId()}));
		    return newNodeVar;
		    }
		    return null;
		}


		public NodeVar get(int value){
		    return this.containerAll.get(value);
		}
		/**
		 * Only adds a new nodeVar if there isn't one
		 * already. The kD value is normalized by the current
		 * upperLimit. A value exceeding the upperLimit is
		 * mapped to the upperLimit.
		 * @param kD
		 */

		public NodeVar addOrRetrieve(int kD){
		    NodeVar nodeVar = this.containerAll.get(kD);
		    if(nodeVar == null){
			nodeVar = this.add(kD, 0, false, true);
			nodeVar.iAmFresh = true;
		    }
		    return nodeVar;
		}

		/**
		 * Given iKD, returns the Id of the ceiling nodevar,
		 * that is, the id of the entry with a key that is
		 * larger or equal to iKD.
		 */
		private int getCeilingId(int iKD){
		    Integer key = this.containerAll.ceilingKey(iKD);
		    if(key == null)return 0;
		    return this.containerAll.get(key).id;
		}

		/**
		 * Given iKD, returns the Id of the ceiling nodevar,
		 * that is, the id of the entry with a key that is
		 * larger or equal to iKD.
		 */

		public int getCeilingId (){
		    return this.containerAll.lastEntry().getValue().id;
		}

		/**
		 * Given iKD, returns the value of the ceiling
		 * nodevar, that is, the value iKD of the entry
		 * imediately above iKD.
		 * @param iKD
		 * @return value
		 */

		public int getCeilingKD (int iKD){
		    if(this.containerAll == null)return -1;
		    Integer key =this.containerAll.ceilingKey(iKD);
		    if(key == null)return -1;
		    return this.containerAll.get(key).getKD();
		}

		public Collection<NodeVars.NodeVar> currentUpper(){
		    return this.containerAll.tailMap(upperLimit + 1).values();
		}
		public Collection<NodeVars.NodeVar> currentTail(){
		    return this.containerAll.tailMap(olderUpperLimit).values();
		}

		public Collection<NodeVars.NodeVar> currentHead(){
		    return this.containerAll.headMap(upperLimit+1).values();
		}
		public Collection<NodeVars.NodeVar> olderHead(){
		    return this.containerAll.headMap(olderUpperLimit+1).values();
		}
		public Collection<NodeVars.NodeVar> currentNovel(){
		    return this.containerAll.subMap(olderUpperLimit+1, upperLimit+1).values();
		}
	    }	     

	    NodeVars nodeVars = null;
	    private int nodeSum = 0;
	    private Node left = null;
	    private Node right = null;
	    private int leafID = 0;
	    private int nodeName = 0;
	    

	    public Node(int weight, int X){
		nodes.add(this);
		int sign = 1;
		if(weight < 0)
		     sign = -1;
		this.nodeSum = sign * weight;
		this.left = null; 
		this.right = null;
		this.nodeVars = new NodeVars();
		this.leafID = sign * X;
		//TODO this is the only significant difference with GenTotalEncoder.
		this.nodeVars.add(this.nodeSum, leafID, false, false);
		


	    }
	     
	    public Node(Node left, Node right){
		nodes.add(this);
		this.left = left;
		this.right = right;
		this.nodeVars =  new NodeVars();
		this.nodeSum = left.nodeSum + right.nodeSum;
	    }



	}

	public void setOlderUpperLimit(){
	    this.olderUpperLimit = this.upperLimit;
	}
	public void setUpperLimit(int newUpperLimit){
	    this.upperLimit = newUpperLimit;
	}
	/**
	 *Links the SumTree, in such a fashion that at any time all
	 *unlinked nodes are lighter than any linked node.
	 */
	public void linkTreeNameNodes(){
	    int size = unlinkedNodes.size();
	    int name = 0;
	    while(size >=2){
	 	Node leftNode = unlinkedNodes.poll();
		leftNode.nodeName = name;
		name++;
		Node rightNode = unlinkedNodes.poll();
		rightNode.nodeName = name;
		name++;
		Node parentNode = new Node(leftNode, rightNode);
		unlinkedNodes.add(parentNode);
		size--;
	    }
	    this.parent = this.unlinkedNodes.poll();
	    this.parent.nodeName = name;

	}

	public SumTree(int iObj, int[] leafWeights){
	    this.iObj = iObj;
	    this.maxUpperLimit = instance.getObj(iObj).getWeightDiff();
	    int iX = 0;
	    for(int weight : leafWeights){
		int ithX = instance.getObj(iObj).getSubObjLits(0).get(iX);
		Node node =  new Node(weight, ithX);
		this.unlinkedNodes.add(node);
		iX ++;
	    }
	    linkTreeNameNodes();

	}
	 
    }

    static class GTEIndex extends Index{
	int nodeName = 0;

	GTEIndex(int iObj, int kD, int nodeName){
	    super(iObj, kD);
	    this.nodeName = nodeName;
	}
	int getNodeName(){
	    return this.nodeName;
	}
    }

    /**
     *Trees used to encode the goal limits
     */

    private SumTree[] sumTrees = null;

    /**
     * Creates an instance of the generalized totalizor encoder
     * @param instance
     * the pseudo boolean instance
     * @param solver
     * the solver to be updated
     */

    public GenTotalEncoder(Instance instance, PBSolver solver) {
	// Log.comment(5, "{ GenTotalEncoder");
	super(instance, solver);
	this.sumTrees = new SumTree[this.instance.nObjs()];
	for(int iObj = 0, nObj = instance.nObjs() ;iObj< nObj; ++iObj){
	    Objective ithObj = this.instance.getObj(iObj);
	    ReadOnlyVec<Real> ithObjCoeffsReal = ithObj.getSubObjCoeffs(0);
	    int[] ithObjCoeffsInt = new int[ithObjCoeffsReal.size()];

	    for(int iX = 0, nX = ithObjCoeffsReal.size(); iX < nX; ++iX){
		ithObjCoeffsInt[iX] = Math.round(ithObjCoeffsReal.get(iX).asInt());
	    }
	    this.sumTrees[iObj] = new SumTree(iObj ,ithObjCoeffsInt);
	}

	for(int iObj = 0, nObj = instance.nObjs() ;iObj< nObj; ++iObj){
	    for(SumTree.Node node: this.sumTrees[iObj].nodes)
		node.nodeVars.add(0, 0, false, false);
	}
	// Log.comment(5, "}");
    }

    public int getCurrentKD(int iObj){
	return this.sumTrees[iObj].upperLimit;
    }


    /**
     * get iObj from an Y variable
     */

    public int getIObjFromY(int id){
	assert this.isY(id);
	return getIObjFromS(id);
    }

    /**
     * get kD from an Y variable
     */
    public int getKDFromY(int id){
	assert this.isY(id);
	return this.getKDFromS(id);

    }

    /**
     * get iObj from an S variable
     */

    public int getIObjFromS(int id){
	assert this.isS(id);
	return this.librarian.getIndex(id).getIObj();
    }


    /**
     * get kD from an S variable
     */

    public int getKDFromS(int id){
	assert this.isS(id);
	return this.librarian.getIndex(id).getKD();
    }

    /**
     * get name of the node, from an S variable
     */

    public int getNameFromS(int id){
	assert this.isS(id);
	return this.librarian.getIndex(id).getNodeName();
    }


    /**
     *Tricky. This is not a real getter. Given kD, it returns the id
     *of the variable with a smaller kD, yet larger or equal to kD.
     */
    public int getY(int iObj, int kD){
	return this.sumTrees[iObj].parent.nodeVars.getCeilingId(kD);
    }


    /**
     *Checks if the variable in the literal is an Y variable.
     */
    public boolean isY(int literal){
	int id = this.solver.idFromLiteral(literal);
	if(!isS(literal))
	    return false;
	for(SumTree sumTree: this.sumTrees)
	    for(SumTree.Node.NodeVars.NodeVar nodeVar: sumTree.parent.nodeVars.containerAll.values())
		if(nodeVar.getId() == id)
		return true;
	return false;
    }


    /**
     *Checks if a variable is an S variable, given a literal.
     */

    public boolean isS(int literal){
	if(isX(literal))
	    return false;
	return true;
    }

    /**
     *Pretty print the variable in literal. 
     */
    public String prettyFormatVariable(int literal){
	int sign =(literal>0)? 1: -1;
	int id =  literal * sign;

	if(this.isY(id)){
	    int iObj = this.getIObjFromS(id);
	    int kD = this.getKDFromS(id);
	    int k = kD; // + this.instance.getObj(iObj).getMinValue();
	    return "Y[" + iObj + ", " + k +"]"+ "::" + literal + " ";
	}
	 
	if(isX(id)){
	    return (sign>0? "+":"-")+"X["+id+"] ";
	}
	/**
	 *Then, it is S!
	 */
	int iObj = this.getIObjFromS(id);
	int kD = this.getKDFromS(id);
	int k = kD ;// + this.instance.getObj(iObj).getMinValue();
	int name = this.getNameFromS(id);
	return "S[" + iObj +  ", " + name  +", "+ k +"]"+ "::" + literal + " ";
    }



    /**
     * Add the sequential clauses. This are the clauses of the form v1
     *=> v2,where v2 belongs to the same node and is associated to a
     *smaller value kD.
     */
    private boolean addClauseSequential(SumTree.Node root){
        // Log.comment(5, "{ GenTotalEncoder.addClauseSequential");
	boolean change = false;

	SumTree.Node.NodeVars.NodeVar past;
	SumTree.Node.NodeVars.NodeVar current;
	Collection<SumTree.Node.NodeVars.NodeVar> tail =
	    root.nodeVars.currentTail();
	Iterator<SumTree.Node.NodeVars.NodeVar> it = tail.iterator();
	if(it.hasNext()){
	    past = it.next();
	    if(it.hasNext()){
	    do{
		    current = it.next();
		    if(current.iAmFresh() || past.iAmFresh()){
			IVecInt clause = new VecInt(new int[] {-current.id, past.id});
			AddClause(clause);
			past.iAmFresh = false;
			change = true;
		    }
		    past = current;
	    }while(it.hasNext());
	    current.iAmFresh = false;
	    }
	}
        // Log.comment(5, "}");
	return change;
    }
    /**
     * This adds the clause that makes this an GTE. That is, v1 v2 =>
     * v3, where kD of v3 is the (corrected) sum kD of v1 and v2
     */
    
    private boolean addSumClauses(SumTree.Node parent, SumTree.Node first, SumTree.Node second){
	// Log.comment(5, "{ GenTotalEncoder.addSumClauses");
	boolean change = false;
	Collection<SumTree.Node.NodeVars.NodeVar> firstAll =
	    first.nodeVars.currentHead();
	
	Collection<SumTree.Node.NodeVars.NodeVar> secondPartial =
	    second.nodeVars.currentNovel();
	
	for(SumTree.Node.NodeVars.NodeVar firstVar : firstAll){
	    for(SumTree.Node.NodeVars.NodeVar secondVar : secondPartial ){
		SumTree.Node.NodeVars.NodeVar parentVar =
		    parent.nodeVars.addOrRetrieve(firstVar.kD + secondVar.kD);
		if(parentVar != null && parentVar.newValidVariable()) 
		    if(parentVar.getKD()!=0 ){
			IVecInt clause = new VecInt(new int[]{parentVar.getId()});
			if(firstVar.getKD()>0)
			    clause.push(-firstVar.getId());
			if(secondVar.getKD()>0)
			    clause.push(-secondVar.getId());
			AddClause(clause);
			change = true;
		    }
	    }
	}
	firstAll =
	    first.nodeVars.currentNovel();
	
	secondPartial =
	    second.nodeVars.olderHead();
	
	for(SumTree.Node.NodeVars.NodeVar firstVar : firstAll){
	    for(SumTree.Node.NodeVars.NodeVar secondVar : secondPartial ){
		SumTree.Node.NodeVars.NodeVar parentVar =
		    parent.nodeVars.addOrRetrieve(firstVar.kD + secondVar.kD);
		if(parentVar != null && parentVar.newValidVariable()) 
		    if(parentVar.getKD()!=0 ){
			IVecInt clause = new VecInt(new int[]{parentVar.getId()});
			if(firstVar.getKD()>0)
			    clause.push(-firstVar.getId());
			if(secondVar.getKD()>0)
			    clause.push(-secondVar.getId());
			AddClause(clause);
			change = true;
		    }
	    }
	}
	
	// Log.comment(5, "}");
	return change;
    }
    /**
     * Simple propagation of variables above the current limit
     * 
     */
    
    private boolean simplePropagation(SumTree.Node parent){
	ArrayList<SumTree.Node> children = new ArrayList<SumTree.Node>(2);
	children.add(parent.left);
	children.add(parent.right);
	// Log.comment(5, "{ GenTotalEncoder.simplePropagation");
	boolean change = false;
	
	for(SumTree.Node child: children){
	    Collection<SumTree.Node.NodeVars.NodeVar> childTail =
		child.nodeVars.currentUpper();
	    for(SumTree.Node.NodeVars.NodeVar childVar : childTail){
		if(childVar.iAmFresh){
		    childVar.iAmFresh = false;
		    SumTree.Node.NodeVars.NodeVar parentVar = parent.nodeVars.addOrRetrieve(childVar.getKD());
		    IVecInt clause = new VecInt(new int[]{-childVar.getId(), parentVar.getId()});
		    
		    AddClause(clause);
		    change = true;
		}
	    }
	}
	// Log.comment(5, "}");
	return change;
    }


    /**
     *Recursive helper of addClausesSumTree
     */

    public boolean addClausesSubSumTree(SumTree sumTree, SumTree.Node currentNode, boolean secondPhase){
	boolean change = false;
	SumTree.Node left = currentNode.left;
	SumTree.Node right = currentNode.right;
	if(currentNode.leafID != 0){
	    return false;
	}
	else{
	    change = addClausesSubSumTree(sumTree, left, secondPhase) || change;
	    change = addClausesSubSumTree(sumTree, right, secondPhase) || change;
	    change = addSumClauses(currentNode, left, right) || change;    
	    change = simplePropagation(currentNode) || change;

	    // else
	    // 	change = addBindingInternal(sumTree, currentNode, left, right);
	}
	return change;
    }

    /**
     *Adds all clauses, respecting the current upperLimit, that complete the semantics of the GTE 
     */
    public boolean addClausesSumTree(int iObj){
	boolean change = false;
	SumTree ithObjSumTree = this.sumTrees[iObj];
	change = addClausesSubSumTree(ithObjSumTree, ithObjSumTree.parent, false) || change;


	// if(change)
	//     addClausesSubSumTree(ithObjSumTree, ithObjSumTree.parent, true);
	return change;
    }




    /**
     *Updates the semantics, in such a way that everything is valid
     *for any value kD less or equal to upperKD. Notice that it may
     *well be extend more, depending on the possible sums
     */

    public boolean UpdateCurrentK(int iObj, int upperKD){
	boolean change = false;
	SumTree ithObjSumTree = this.sumTrees[iObj];
	if(upperKD > this.getCurrentKD(iObj)){
	ithObjSumTree.setOlderUpperLimit();
	
	while(!change && upperKD <= this.instance.getObj(iObj).getWeightDiff()){
	    // Log.comment(5, "{ GenTotalEncoder.UpdateCurrentK of "+ iObj + " to " + upperKD);
	    this.sumTrees[iObj].setUpperLimit(upperKD);
	    change = addClausesSumTree(iObj);
	    upperKD++;
	}
	if(change)
	    addClauseSequential(ithObjSumTree.parent );
}
    // Log.comment(5, "}");
    return change;
    }

    /**
     * Generate the upper limit assumptions
     */
    public IVecInt generateUpperBoundAssumptions(int[] UpperKD){

	IVecInt assumptions = new VecInt(new int[]{});
	for(int iObj = 0; iObj < this.instance.nObjs(); ++iObj){
	    Objective ithObjective = this.instance.getObj(iObj);
	    if(UpperKD[iObj]  < ithObjective.getWeightDiff())
		assumptions.push(-this.getY(iObj, UpperKD[iObj] + 1));
	    
	    ReadOnlyVecInt objectiveLits = ithObjective.getSubObjLits(0);
	    ReadOnlyVec<Real> objectiveCoeffs = ithObjective.getSubObjCoeffs(0);
	    int sign = 1;
	    int ithAbsoluteWeight;
	    for(int iX = 0, nX = ithObjective.getTotalLits(); iX <nX; iX ++){
		ithAbsoluteWeight = objectiveCoeffs.get(iX).asInt();
		sign = (ithAbsoluteWeight > 0? 1 : -1);
		ithAbsoluteWeight *= sign;
		if( ithAbsoluteWeight > UpperKD[iObj])
		    assumptions.push(- sign * objectiveLits.get(iX));
	    }

	}
	return assumptions;
    }


    /**
     *Finds the next valid kD value, starting in kD and extending
     *until newKD, inclusive. This will not repeat clauses only if the
     *intervale (kD, newKD] is empty of already computed kD values
     *
     */
    @Override
    public int generateNext(int iObj, int kD, int inclusiveMax){
	// Log.comment(5, "{ GenTotalEncoder.generateNext");

	boolean change = false;
	SumTree ithObjSumTree = this.sumTrees[iObj];
	// store  values of upperLimit and olderUpperLimit
	int upperLimit = ithObjSumTree.upperLimit;
	int olderUpperLimit = ithObjSumTree.olderUpperLimit;
	ithObjSumTree.olderUpperLimit = kD;

	ithObjSumTree.setUpperLimit(kD);
	int upperKD = kD;
	while(!change && upperKD < inclusiveMax){
	    upperKD++;
	    this.sumTrees[iObj].setUpperLimit(upperKD);
	    // Log.comment(5, "{ GenTotalEncoder.generateNext of "+ iObj + "from " + kD + " to " + upperKD);
	    change = addClausesSumTree(iObj);
	    // Log.comment(5, "}");
	}
	if(change)
	    addClauseSequential(ithObjSumTree.parent);

	// Log.comment(5, "}");
	ithObjSumTree.upperLimit = upperLimit;
	ithObjSumTree.olderUpperLimit = olderUpperLimit;
	if(change)
	    return ithObjSumTree.upperLimit;
	else
	    return kD;
    }

    /**
     *Finds the next valid kD value, starting from lastK and extending
     *until newKD, inclusive. This will not repeat clauses only if the
     *intervale (kD, newKD] is empty of already computed kD values
     *
     */

    public int nextKDValue(int iObj, int kD){
	SumTree ithObjSumTree= this.sumTrees[iObj];
	if(kD == ithObjSumTree.maxUpperLimit)
	    return kD;
	if(ithObjSumTree.parent == null)
	    return -1;
	int aproxNextKD = ithObjSumTree.parent.nodeVars.getCeilingKD(kD + 1);
	if(aproxNextKD == -1)
	    return kD;
	return aproxNextKD;
    }
}
