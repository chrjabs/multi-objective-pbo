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
package org.sat4j.moco.translator;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Map;

import org.sat4j.core.VecInt;
import org.sat4j.moco.goal_delimeter.SelectionDelimeterMSU3;
import org.sat4j.moco.pb.PBSolver;
import org.sat4j.moco.problem.Instance;
import org.sat4j.specs.ContradictionException;
import org.sat4j.specs.IVecInt;

/**
 * Class that implements UnsatSat, MSU3 flavoured
 * @author João Cortes
 */

public class Translator{
    protected PBSolver solver = null;
    private int nClauses = 0;
    private SelectionDelimeterMSU3 selectionDelimeter = null;
    private Instance problem = null;
    private int nVars = 0;
    BufferedWriter out = null;
    BufferedWriter tempOut = null;

    public Translator(Instance m, BufferedWriter out, BufferedWriter tempOut) {
	this.out = out;
	this.tempOut = tempOut;
	this.problem = m;
	try {
            this.solver = buildSolver();
        }
        catch (ContradictionException e) {
            // Log.comment(3, "Contradiction in ParetoMCS.buildSolver");
            return;
        }
	
	this.nVars = this.problem.nVars();
	    
    }

    public void translate(File tempFile, Map<Integer, Integer> upperLimits, Map<Integer, Integer[]> allRatios)  throws IOException{
	
	this.selectionDelimeter = new SelectionDelimeterMSU3(this.problem, solver, false, upperLimits){
		public boolean AddClause(IVecInt setOfLiterals){
		    AddClause1(setOfLiterals);
		    return true;
		    }
		};

	this.selectionDelimeter.initializeObjectManagers();
	for(int i = 0, n = this.problem.nObjs();i<n;i++){
	    if(allRatios.get(i) != null)
		this.selectionDelimeter.getIthObjManager(i).getDigitalEnv().setRatios(allRatios.get(i));
	}
	this.selectionDelimeter.buildCircuits();
	this.tempOut.close();
	this.printHeaderLine();
	this.selectionDelimeter.printBasis(out);
	this.selectionDelimeter.printOutVariables(this.out);
	this.printOriginalVariables();

	try{
	this.mergeFiles(tempFile);
	tempFile.delete();
	    
	}
	catch (IOException e){
	    System.out.println("Couldn't merge files");
	}
	this.out.close();
	
    }

    public void printHeaderLine(){
	int nOutVars = this.selectionDelimeter.numberOutVars();
	int nInVars = this.nVars;
	int nObj = this.problem.nObjs();
	int nTotalVars = this.solver.nVars();
	try{
	this.out.write("p mocnf " + nTotalVars + " " + this.nClauses + " " + nObj + " " + nInVars + " " + nOutVars + "\n");
	}
	catch(IOException e){}
    }

    public void printOriginalVariables() throws IOException{
	
	for(int i = 1, n =this.nVars; i <= n; i++){
	    //TODO: variables may not be contiguous...
	    this.out.write("x x"+ i + " " + i + "\n");
	}


    }

    /**
     * Creates a PB oracle initialized with the MOCO's constraints.
     * @return The oracle.
     * @throws ContradictionException if the oracle detects that the
     * MOCO's constraint set is unsatisfiable.
     */
    protected PBSolver buildSolver() throws ContradictionException {
        // Log.comment(5, "in Algorithm.buildSolver");
        PBSolver solver = new PBSolver();
        solver.newVars(this.problem.nVars());
        for (int i = 0; i < this.problem.nConstrs(); ++i) {
            solver.addConstr(this.problem.getConstr(i));
        }
	solver.setConstantID();
	this.AddClause1(new VecInt(new int[]{solver.constantLiteral(true)}));
        // Log.comment(5, "out UnsatSat.buildSolver");
        return solver;
    }
    void mergeFiles(File tempFile) throws IOException
    { 
	BufferedReader br = new BufferedReader(new FileReader(tempFile.getAbsolutePath())); 
          
	String line = br.readLine(); 
          
	// loop to copy each line of  
	// file1.txt to  file3.txt 
	while (line != null) 
	    { 
		this.out.write(line); 
		this.out.write("\n");
		line = br.readLine(); 
	    } 
          
	br.close(); 
    }
    public boolean AddClause1(IVecInt setOfLiterals){
    try{
	{int i = 0, n = setOfLiterals.size();
	    for(; i < n - 1; i++)
		tempOut.write(setOfLiterals.get(i) + " ");
	    if(n > 0)
		tempOut.write(setOfLiterals.get(i) + "\n");
	}
    }
	catch (IOException e){
	}
    this.nClauses++;
    return true;
    }
}

 
