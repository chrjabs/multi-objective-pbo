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

import java.lang.Math;
import org.sat4j.moco.util.Log;
import java.util.Iterator;

import java.util.PriorityQueue;
import java.util.ArrayList;
import java.util.TreeMap;
import java.util.Collection;
import java.util.HashMap;

import org.sat4j.core.ReadOnlyVec;

import org.sat4j.core.ReadOnlyVecInt;
import org.sat4j.core.VecInt;
import org.sat4j.moco.util.Real;
import org.sat4j.moco.goal_delimeter.GenTotalEncoder.GTEIndex;
import org.sat4j.moco.pb.PBSolver;
import org.sat4j.moco.problem.Instance;
import org.sat4j.moco.problem.Objective;
import org.sat4j.specs.IVecInt;


/**
 * Class with the implementation of the generalized totalizor encoder.
 * Notice that the differential is both a value and an index starting at 0
 * @author Joao O'Neill Cortes
 */


public class GenTotalEncoderMSU3 extends GoalDelimeter<GTEIndex> {
    
    /**
     *The inverse index map for the partial sum variables. For each ID, a
     *value that is an array vector with the value of the goal and the
     *value of the sum, namely {kD, iObj, nodeName}
     */
    private HashMap<Integer, Boolean> coveredLiterals = null;

    /**
     * Upper bound, exclusive
     */
    private int[] UpperBound = null;


    /**
     *gets the current upper bound
     *@param iObj
     */

    /**
     * Last explored differential k, for each objective function.
     */
    private int[] upperKD = null;


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
	private int olderUpperLimit = -1;

	/**
	 *Maximal possible value
	 */
	private int maxUpperLimit = 0;

	/**
	 *Max uncovered KD value
	 */
	private int maxUncoveredKD = 0;

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
	 *The root of the SumTree added.
	 */

	private Node freshParent = null;

	private int depth = 0;
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
			librarian.putIndex(newNodeVar.getId(), new GTEIndex(iObj,kD,  nodeName));
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
		    
		    Integer key =  this.containerAll.ceilingKey(iKD);
		    if(key == null)
			return 0;
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
		    if(this.containerAll == null)
			return -1;
		    Integer key =  this.containerAll.ceilingKey(iKD);
		    if(key == null)
			return -1;
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
	    private int leafLit = 0;
	    private int nodeName = 0;
	    
	    public Node(){
		this.nodeName = nodes.size();
		nodes.add(this);
		this.nodeVars = new NodeVars();
		this.nodeVars.add(0, 0, false, false);
	    }
	    

	    public Node(int weight, int X){
		this();
		int sign = 1;
		if(weight < 0)
		    sign = -1;
		this.nodeSum = sign * weight;
		this.left = null; 
		this.right = null;
		this.leafLit = sign * X;
		// this.nodeVars.add(this.nodeSum, leafLit, false, false);
		


	    }
	     
	    public Node(Node left, Node right){
		this();
		this.left = left;
		this.right = right;
		this.nodeSum = left.nodeSum + right.nodeSum;
	    }
	    public boolean isLeaf(){
		return this.leafLit!=0;
	    }

	}

	public int getMaxUpperLimit(){return this.maxUpperLimit;}
	public int getMaxUncoveredKD(){return this.maxUncoveredKD;}
	public void setMaxUncoveredKD(int a){this.maxUncoveredKD = a;}
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
	public SumTree.Node linkTreeNameNodes(){
	    // Log.comment(5, " { GenTotalEncoderMSU3.linkTreeNameNodes");
	    int size = unlinkedNodes.size();
	    // int name = nodes.size()-size;
	    int depth = 0;
	    Node newParent = null;
	    if(size>=1){
		while(size >=2){
		    Node leftNode = unlinkedNodes.poll();
		    Node rightNode = unlinkedNodes.poll();
		    Node parentNode = new Node(leftNode, rightNode);
		    unlinkedNodes.add(parentNode);
		    size--;
		    depth++;
		}
		newParent = this.unlinkedNodes.poll();
		
		// newParent.nodeName = name;
		if(this.parent == null){
		    this.parent = newParent;
		    this.depth = depth;
		}
		else{
		    this.parent = new Node(this.parent, newParent);
		    if(this.depth < depth)
			this.depth = depth;
		    this.depth++;
		    // this.parent.nodeName = name + 1;
		}
	    }
	    // Log.comment(5, "}");
	    return newParent;
	}
	public SumTree(int iObj){
	    this.iObj = iObj;
	    this.maxUpperLimit = instance.getObj(iObj).getWeightDiff();

	}
	public boolean isLeafAlreadyHere(int lit){
	    boolean alreadyHere = false;
	    for(Node node: this.nodes){
		int leaflitTrue = node.leafLit;
		if(leaflitTrue == lit){
		    alreadyHere = true;
		    break;
		}
	    }	    
	    return alreadyHere;
	}

	
	/**
	 * push a set of literals into the sum tree.
	 */
	public boolean pushNewLeafs(IVecInt litsToAdd, boolean invertLits){
	    boolean alreadyHere = false;
	    for(int i = 0, n = litsToAdd.size(); i < n; i++){
		int lit = litsToAdd.get(i), lit1 = 0;
		Real weight = instance.getObj(iObj).getSubObj(0).weightFromLit(lit);
		if(invertLits){
		    lit1 = lit;
		    if(weight==null){
			weight = instance.getObj(iObj).getSubObj(0).weightFromLit(-lit);
			lit1 = -lit;
		    }		    
		}
		if(weight!=null){
		    if(invertLits){
			int sign = 1;
			if(weight.asInt() < 0) sign = -1;
			if(lit != - sign * lit1 )
			    continue;
			else
			    lit = lit1;
	    }
		    alreadyHere= this.isLeafAlreadyHere(lit);
		    if(!alreadyHere){
			Node node =  new Node(weight.asIntExact(), lit);
			// Log.comment(5, "new leaf: "+ node.leafLit);
			this.unlinkedNodes.add(node);
		    }
		}
	    }
	    


	    return this.unlinkedNodes.size()!=0;
	}
	/**
	 *Adds a new sub tree, with nodes associated to leafs in leafsXId
	 */
	public SumTree.Node linkNewNodes(){
	    Node newParent = linkTreeNameNodes();
	    return newParent;
	}

	/**
	 *report unbalance
	 */
	public double reportUnbalance(){
	    int leafsN = 0;
	    double expectedDepth = 0;
	    for(Node node: this.nodes){
		int leaflitTrue = solver.idFromLiteral(node.leafLit);
		if(leaflitTrue != 0){
		    leafsN++;
		}
		expectedDepth = Math.log(leafsN)/Math.log(2);
		// Log.comment("expected depth: " + expectedDepth );
		// Log.comment("true depth: " + this.depth);

	    }
	    return expectedDepth/this.depth;
	}


	/**
	 *Updates the current maxValues for each objective
	 */

	public void updateUncoveredMaxKD(){
	    int a = 0;
	    Objective ithObjective = instance.getObj(this.iObj); // 
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
	    this.setMaxUncoveredKD(a);
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

    public GenTotalEncoderMSU3(Instance instance, PBSolver solver) {
	// Log.comment(5, "{ GenTotalEncoder");
	this.instance = instance;	
	this.solver = solver;
	this.firstVariable = this.solver.nVars() + 1;
	this.sumTrees = new SumTree[this.instance.nObjs()];
	this.coveredLiterals = new HashMap<Integer, Boolean>(this.solver.nVars());
	this.upperKD =  new int[(this.instance.nObjs())];
	this.UpperBound =  new int[(this.instance.nObjs())];
	    
	for(int iObj = 0, nObj = instance.nObjs() ;iObj< nObj; ++iObj){
	    this.sumTrees[iObj] = new SumTree(iObj);
	}
	this.initializeCoveredLiterals();
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
		ithAbsoluteWeight = objectiveCoeffs.get(iX).asInt();
		sign = (ithAbsoluteWeight > 0? 1 : -1);
		ithAbsoluteWeight *= sign;
		this.coveredLiterals.putIfAbsent(-sign * ithX, true);
	    }
	}

    }

    public int getCurrentKD(int iObj){
	return this.sumTrees[iObj].upperLimit;
    }

    public void reportUnbalances(){
	String report = "unbalances: ";
	for(SumTree sumTree: this.sumTrees)
	    report += " " + sumTree.reportUnbalance();
	Log.comment(2, report);
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
     *of the variable with a smaller kD, yet larger or equal than kD.
     */
    public int getY(int iObj, int kD){
	if(this.sumTrees[iObj].parent!=null)
	    return this.sumTrees[iObj].parent.nodeVars.getCeilingId(kD);
	return 0;
    }


    /**
     *Checks if the variable in the literal is an Y variable.
     */
    public boolean isY(int literal){
	int id = this.solver.idFromLiteral(literal);
	// if(!isS(literal))
	//     return false;
	for(SumTree sumTree: this.sumTrees)
	    if(sumTree.parent!=null)
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

	if(isX(id)){
	    return (sign>0? "+":"-")+"X["+id+"] ";
	}
	if(this.isY(id)){
	    int iObj = this.getIObjFromS(id);
	    int kD = this.getKDFromS(id);
	    int k = kD; // + this.instance.getObj(iObj).getMinValue();
	    return "Y[" + iObj + ", " + k +"]"+ "::" + literal + " ";
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
			// past.iAmFresh = false;
			change = true;
		    }
		    past = current;
		}while(it.hasNext());
		// current.iAmFresh = false;
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
		if(parentVar != null) 
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
		if(parentVar != null) 
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
    
    private boolean simplePropagation(SumTree.Node parent, boolean ignoreFreshness){
	// Log.comment(5, "{ GenTotalEncoder.simplePropagation");
	ArrayList<SumTree.Node> children = new ArrayList<SumTree.Node>(2);
	children.add(parent.left);
	children.add(parent.right);
	boolean change = false;
	
	for(SumTree.Node child: children){
	    Collection<SumTree.Node.NodeVars.NodeVar> childTail =
		child.nodeVars.currentUpper();
	    for(SumTree.Node.NodeVars.NodeVar childVar : childTail){
		if(ignoreFreshness || childVar.iAmFresh){
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
    public boolean addClausesCurrentNode(SumTree sumTree, SumTree.Node currentNode, boolean iAmtheNewRoot){
	// Log.comment(5, "{ GenTotalEncoder.addClausesCurrentNode");
	boolean change = false;
	SumTree.Node left = currentNode.left;
	SumTree.Node right = currentNode.right;
	change = addSumClauses(currentNode, left, right) || change;    
	change = simplePropagation(currentNode, iAmtheNewRoot) || change;
	// Log.comment(5, "}");
	return change;
    }

    /**
     *Recursive helper of addClausesSumTree
     */

    public boolean addClausesSubSumTree(SumTree sumTree, SumTree.Node currentNode){
	// Log.comment(5, "{ GenTotalEncoder.addClausesSubSumTree");
	
	boolean change = false;
	SumTree.Node left = currentNode.left;
	SumTree.Node right = currentNode.right;
	SumTree.Node.NodeVars.NodeVar nodeVar = null;
	if(currentNode.isLeaf()){
	    if(sumTree.olderUpperLimit < currentNode.nodeSum &&  currentNode.nodeSum <= sumTree.upperLimit){
		nodeVar = currentNode.nodeVars.addOrRetrieve(currentNode.nodeSum);
		if((nodeVar != null) && nodeVar.iAmFresh){
		    AddClause(new VecInt(new int[]{-currentNode.leafLit, nodeVar.getId()}));
		    // Log.comment(5, "}");
		    return true;
		}
		return false;
	    }
	}
	else{
	    change = addClausesSubSumTree(sumTree, left) || change;
	    change = addClausesSubSumTree(sumTree, right) || change;
	    change = addClausesCurrentNode(sumTree, currentNode, false);
	    // else
	    // 	change = addBindingInternal(sumTree, currentNode, left, right);
	}
	// Log.comment(5, "}");

	return change;

    }

    /**
     *Adds all clauses, respecting the current upperLimit, that complete the semantics of the GTE 
     */
    public boolean addClausesSumTree(int iObj){
	boolean change = false;
	SumTree ithObjSumTree = this.sumTrees[iObj];
	if(ithObjSumTree.parent!=null)
	    change = addClausesSubSumTree(ithObjSumTree, ithObjSumTree.parent) || change;


	// if(change)
	//     addClausesSubSumTree(ithObjSumTree, ithObjSumTree.parent);
	return change;
    }


    /**
     *Adds leafs to tree
     */
    public boolean addLeafs(int iObj, IVecInt explanationX ){
	// Log.comment(5, "{ GenTotalEncoder.addLeafs: " + iObj);
	SumTree ithObjSumTree = this.sumTrees[iObj];

	boolean result = ithObjSumTree.pushNewLeafs(explanationX, true);
	ithObjSumTree.freshParent = ithObjSumTree.linkNewNodes();
	ithObjSumTree.updateUncoveredMaxKD();
	// Log.comment(5, "}");

	return result;
    }
    

    public boolean isLeafAlreadyHere(int iObj, int lit){
	return this.sumTrees[iObj].isLeafAlreadyHere(lit);
    }

    /**
     *clause fresh subTree, whose nodes were already linked by addLeafsTo
     */

    public boolean bindFreshSubTree(int iObj, int upperLimit){
	// Log.comment(5, "{ GenTotalEncoder.bindFreshSubTree: " + iObj + ", " + upperLimit );
	boolean change = false;
	SumTree ithObjSumTree = this.sumTrees[iObj];
	SumTree.Node newParent = ithObjSumTree.freshParent;
	if(newParent!=null){
	    int oldOldUpperLimit = ithObjSumTree.olderUpperLimit;
	    int oldUpperLimit = ithObjSumTree.upperLimit;

	    ithObjSumTree.olderUpperLimit = -1;
	    ithObjSumTree.upperLimit = upperLimit;

	    change = this.addClausesSubSumTree(ithObjSumTree, newParent);
	    if(newParent!=ithObjSumTree.parent){
		change = this.addClausesCurrentNode(ithObjSumTree, ithObjSumTree.parent, true) || change;
	    }
	    this.addClauseSequential(ithObjSumTree.parent);
	    ithObjSumTree.olderUpperLimit = oldOldUpperLimit;
	    ithObjSumTree.upperLimit = oldUpperLimit;
	}
	// Log.comment(5, "}");
	return change;
    }


    /**
     *Finds the next valid kD value, starting in kD and extending
     *until newKD, inclusive. This will not repeat clauses only if the
     *intervale (kD, newKD] is empty of already computed kD values
     *
     */

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


	public int getMaxUpperLimit(int iObj){return this.sumTrees[iObj].getMaxUpperLimit();}

	public int getMaxUncoveredKD(int iObj){return this.sumTrees[iObj].getMaxUncoveredKD();}



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

    /**
     *Uncover leafs
     */
    private boolean uncoverXs(IVecInt explanationX) {
	// Log.comment(5, "{ UnsatSatMSU3.uncoverXs");

	boolean change = false;
	for(int iObj = 0; iObj < this.instance.nObjs(); ++iObj){
	    change = this.addLeafs(iObj, explanationX) || change;
	    change = this.bindFreshSubTree(iObj, this.getUpperBound(iObj)) || change;
	}

	int[] explanationXarray = explanationX.toArray();
	for(int x : explanationXarray)
	    this.coveredLiterals.remove(x);
	this.updateAllUncoveredMaxKD();
	this.logUncoveredMaxKD();
	// Log.comment(5, "}");
	return change;
    }
private void updateAllUncoveredMaxKD(){
    for(int i = 0, n = this.instance.nObjs(); i < n; i++)
	this.sumTrees[i].updateUncoveredMaxKD();
}


    /**
     *Sets the current upper bound of iObj to nowKD
     *@param newKD
     *@param iObj
     */
    private void setUpperBound(int iObj, int newKD){
	    this.UpperBound[iObj] = newKD;
    }

    /**
     *gets the current upper limit of the explored value of the
     *differential k of the ithOjective
     *@param iObj
     */

    public int getUpperKD(int iObj){
	return this.upperKD[iObj];
    }

    /**
     *Sets the current upper limit of the explored value of the
     *differential k of the ithOjective to newKD
     *@param newKD
     *@param iObj
     */
    public void setUpperKD(int iObj, int newKD){
	if(this.getUpperKD(iObj)< newKD)
	    this.upperKD[iObj] = newKD;
    }

    /**
     *Log the current uncovered max values
     */

    public void logUncoveredMaxKD(){
	String logUpperLimit = "uncovered max: ["+this.getUncoveredMaxKD(0);
	for(int iObj = 1; iObj < this.instance.nObjs(); ++iObj)
	    logUpperLimit +=", "+this.getUncoveredMaxKD(iObj) ;//+ this.instance.getObj(iObj).getMinValue())
	//..log
	
	logUpperLimit +="]";
	Log.comment(2, logUpperLimit );
    }
    private int getUpperBound(int iObj){
	return this.UpperBound[iObj];
    }

    public void setUpperBound() {
	for(int i = 0, n = this.UpperBound.length; i < n; i++)
	    this.UpperBound[i] = this.nextKDValue(i, this.getUpperKD(i));
	this.setFirstVariable();
    }
    
    public int getUncoveredMaxKD(int iObj){return this.sumTrees[iObj].getMaxUncoveredKD();}
    public void setMaxUncoveredKDs(int iObj, int a){this.sumTrees[iObj].setMaxUncoveredKD(a);}
    public void updateUncoveredMaxKD(int iObj){this.sumTrees[iObj].updateUncoveredMaxKD();}



	@Override
	public IVecInt generateUpperBoundAssumptions(IVecInt explanation, boolean checkChange) {
	    IVecInt assumptions =  super.generateUpperBoundAssumptions(explanation, checkChange);
	    for(Integer x: this.coveredLiterals.keySet())
		assumptions.push(x);
	    return assumptions;
	}
}
