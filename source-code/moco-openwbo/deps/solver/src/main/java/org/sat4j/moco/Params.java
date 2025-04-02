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
package org.sat4j.moco;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Scanner;
import java.util.regex.MatchResult;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.sat4j.moco.util.Real;

/**
 * Class used to store the solver's configuration.
 * @author Miguel Terra-Neves
 */
// TODO: implement setters?
public class Params {

    /**
     * default output for translation
     */
    private static final String DEFAULT_OUTPUT = "out.mocnf";

    /**
     *Default blankness
     */
    private static final String DEFAULT_ISBLANK = "0";

    /**
     * Default verbosity level.
     */
    private static final String DEFAULT_VERB = "0";
    
    /**
     * Default decimal scale.
     */
    private static final String DEFAULT_SCALE = "5";
    
    /**
     * Default literal-weight ratio for stratification.
     */
    private static final String DEFAULT_LWR = "15.0";
    
    /**
     * Default maximum conflicts allowed before partition merging in stratified algorithms.
     */
    private static final String DEFAULT_PMC = "200000";

    /**
     * Default algorithm
     */
    private static final String DEFAULT_ALGI = "0";

    /**
     * default encoding for the goal delimeter
     */
    private static final String DEFAULT_ENC = "GTE";
    
    /**
     * Default trivial threshold (number of trivially solved partitions in a row before merging the
     * remaining ones) for stratified algorithms.
     */
    private static final String DEFAULT_TT = "20";
    
    /**
     * Builds an {@link Options} object with the solver's configuration parameters to be used for parsing
     * command line options.
     * @return An {@link Options} object to be used by the command line interface.
     */
    public static Options buildOpts() {
        Options o = new Options();
	o.addOption("rb", "ratios", true, "ratios for each objective function. Only meaningful with SD encodings.");
	o.addOption("ul", "upperLimits", true, "upper limits for each objective function. Only meaningful with delimeted algorithms");
	o.addOption("ib", "isBlank", true, "Should I translate or solve the instance?");
	o.addOption("o", "output", true, "output file for translation.");
        o.addOption("v", "verbosity", true,
                    "Set the verbosity level (from 0 to 3). Default is " + DEFAULT_VERB + ".");
        o.addOption("t", "timeout", true, "Set the time limit in seconds. No limit by default.");
        o.addOption("sa", "suppress-assign", false, "Suppress assignment output.");
        o.addOption("ds", "decimal-scale", true,
                    "Set the maximum scale (number of digits to the right of the decimal point) for real numbers." +
                    "Default is " + DEFAULT_SCALE + ".");
        o.addOption("s", "stratify", false, "Enable stratification in MCS based algorithm.");
        o.addOption("lwr", "lit-weight-ratio", true,
                    "Set the literal-weight ratio for stratification. Default is " + DEFAULT_LWR + ".");
        o.addOption("pmc", "part-max-confl", true,
                    "Set the maximum conflicts allowed before merging with the next partition in stratified " +
                    "algorithms. Default is " + DEFAULT_PMC + ".");
        o.addOption("tt", "trivial-thres", true,
                    "Set the trivial threshold for stratified algorithms (number of trivially solved partitions " +
                    "in a row before merging the remaining ones). Default is " + DEFAULT_TT + ".");
	o.addOption("alg", "algorithm-index", true, "Choose the algorithm to use. options are 0-paretomcs, 1-unsatsat, 2-pminimal");
	o.addOption("enc", "GD encoding", true, "Choose the goal delimeter encoding to use. options are SWE, GTE");
        return o;
    }
    
    /**
     * map with upper limits for each objective
     */
    Map<Integer, Integer> upperLimits;

    /**
     * map with ratios for each objective
     */
    Map<Integer, Integer[]> allRatios;


    /** 
     * output file name for translation purposes
     */
    private String output = "";

    /**
     * Stores the verbosity level of the solver.
     */
    private int verb = 0;
    
    /**
     * Stores the maximum time, in seconds, allowed for the solver to run.
     * If less than 0, then no time limit is imposed.
     */
    private int timeout = -1;
    
    /**
     * Stores if assignment logging should be suppressed.
     * If true, then assignments should not be logged.
     */
    private boolean suppress_assign = false;
    
    /**
     * Stores the maximum scale (number of digits to the right of the decimal point) to be considered for
     * {@link Real} operations.
     */
    private int scale = 5;
    
    /**
     * Stores if stratification is to be enabled for Pareto-MCS based algorithms.
     */
    private boolean stratify = false;
    
    /**
     * Stores the literal-weight ratio to be used when partitioning objectives for stratification in
     * Pareto-MCS based algorithms.
     */
    private double lwr = 15.0;
    
    /**
     * Stores the maximum conflicts allowed in stratified algorithms before merging a partition with the next
     * one.
     */
    private int pmc = 200000;
    
    /**
     *Stores the algorithm to use.
     */
    private int algorithmI = 0;

    /**
     *Stores the GD encoding to use.
     */
    private String encodingGD = "";

    /**
     * Stores the trivial threshold (number of trivially solved partitions in a row before merging the
     * remaining ones) for stratified algorithms.
     */
    private int tt = 20;
    
    /**
     *Stores the isBlank value
     */
    private int isBlank;
    /**
     * Creates a parameters object with default configuration options.
     */
    public Params() {
	this.output = DEFAULT_OUTPUT;
	this.isBlank = Integer.parseInt(DEFAULT_ISBLANK);
        this.verb = Integer.parseInt(DEFAULT_VERB);
        this.scale = Integer.parseInt(DEFAULT_SCALE);
        this.lwr = Double.parseDouble(DEFAULT_LWR);
        this.pmc = Integer.parseInt(DEFAULT_PMC);
        this.tt = Integer.parseInt(DEFAULT_TT);
	this.algorithmI = Integer.parseInt(DEFAULT_ALGI);
	this.encodingGD = "";
    }
    
    /**
     * Creates a parameters object with the configuration options provided in the command line.
     * @param cl The command line object.
     */
    public Params(CommandLine cl) {
	this.upperLimits = parseUpperLimit(cl.getOptionValue("ul"));
	this.allRatios = parseAllRatios(cl.getOptionValue("rb"));
	this.output = cl.getOptionValue("o", DEFAULT_OUTPUT);
	this.isBlank = Integer.parseInt(cl.getOptionValue("ib", DEFAULT_ISBLANK));
        this.verb = Integer.parseInt(cl.getOptionValue("v", DEFAULT_VERB));
        this.suppress_assign = cl.hasOption("sa");
        this.scale = Integer.parseInt(cl.getOptionValue("ds", DEFAULT_SCALE));
        this.stratify = cl.hasOption("s");
        this.lwr = Double.parseDouble(cl.getOptionValue("lwr", DEFAULT_LWR));
        this.pmc = Integer.parseInt(cl.getOptionValue("pmc", DEFAULT_PMC));
        this.tt = Integer.parseInt(cl.getOptionValue("tt", DEFAULT_TT));
	this.algorithmI = Integer.parseInt(cl.getOptionValue("alg", DEFAULT_ALGI));
	this.encodingGD = cl.getOptionValue("enc", DEFAULT_ENC);
        if (cl.hasOption("t")) {
            this.timeout = Integer.parseInt(cl.getOptionValue("t"));
        }
    }

    private Map<Integer, Integer> parseUpperLimit(String string){
	Map<Integer, Integer> upperLimits = new HashMap<Integer, Integer>();
	if(string ==null)
	    return upperLimits;
	Scanner scanner = new Scanner(string);
	scanner.useDelimiter(",");
	String pattern = "(\\d+)\\:(\\d+)";
	while(scanner.hasNext(pattern)){
	    scanner.next(pattern);
	    MatchResult match = scanner.match();
	    int iObj = Integer.parseInt(match.group(1));
	    int limit = Integer.parseInt(match.group(2));
	    upperLimits.put(iObj, limit);
	}
	scanner.close();
	return upperLimits;
    }

    private Map<Integer, Integer[]> parseAllRatios(String string){
	Map<Integer, Integer[]> allRatios  = new HashMap<Integer, Integer[]>();
	if(string ==null)
	    return allRatios;
	Scanner scanner = new Scanner(string);
	String pattern = "(\\d+)\\:\\[([\\d+,]+)\\]";
	while(scanner.hasNext(pattern)){
	    scanner.next(pattern);
	    MatchResult match = scanner.match();
	    int iObj = Integer.parseInt(match.group(1));
	    Integer[] ratios = this.parseRatios(match.group(2));
	    allRatios.put(iObj, ratios);
	}
	scanner.close();
	return allRatios;
    }

    private Integer[] parseRatios(String string){
	Integer[] ratios;
	List<Integer> ratiosList  = new ArrayList<Integer>(0);
	if(string ==null)
	    return null;
	Scanner scanner = new Scanner(string);
	scanner.useDelimiter(",");
	while(scanner.hasNext())
	    ratiosList.add(scanner.nextInt());
	scanner.close();
	ratios = ratiosList.toArray(new Integer[0]);
	scanner.close();
	return ratios;
    }
    
    /**
     * Retrieves the desired verbosity level.
     * @return The verbosity level.
     */
    public int getVerbosity() { return this.verb; }
    
    /**
     * Checks if a time limit was provided by the user.
     * @return True if a time limit was provided, false otherwise.
     */
    public boolean hasTimeout() { return this.timeout >= 0; }
    
    /**
     * If a time limit was provided by the user, retrieves that time limit in seconds.
     * @return The time limit.
     */
    public int getTimeout() { return this.timeout; }
    
    /**
     * Checks if assignment logging should be suppressed.
     * @return True if assignments are to be suppressed, false otherwise.
     */
    public boolean getSuppressAssignments() { return this.suppress_assign; }
    
    /**
     * Retrieves the maximum scale (number of digits to the right of the decimal point) to be considered for
     * {@link Real} operations.
     * @return The scale.
     */
    public int getScale() { return this.scale; }
    
    /**
     * Checks if stratification is to be enabled for Pareto-MCS based algorithms.
     * @return True if stratification is enabled, false otherwise.
     */
    public boolean getStratify() { return this.stratify; }
    
    /**
     * Retrieves the literal-weight ratio to be used by the stratified Pareto-MCS algorithm in the
     * partitioning process.
     * @return The literal-weight ratio.
     */
    public double getLWR() { return this.lwr; }
    
    /**
     * Retrieves the maximum number of conflicts to be allowed in stratified algorithms before merging a
     * partition with the next one.
     * @return The maximum conflicts per partition.
     */
    public int getPartMaxConfl() { return this.pmc; }
    
    /**
     * Retrieves the trivial threshold (number of trivially solved partitions in a row before merging the
     * remaining ones) for stratified algorithms.
     * @return The trivial threshold.
     */
    public int getTrivialThres() { return this.tt; }
    
    /**
     *Returns the algorithm to be used
     */
    public int getAlgorithmI(){return this.algorithmI;}

    /**
     *returns the goal delimeter (GD) encoding to be used
     @return the name of the encoding of{@code GoalDelimeter} 
    */
    public String getEncodingGD(){return this.encodingGD;}

    public int getIsBlank() {return isBlank;}

    public String getOutput() {
	return output;
    }

    public Map<Integer, Integer> getUpperLimits() {
	return upperLimits;
	}

    public void setUpperLimits(Map<Integer, Integer> upperLimits) {
	this.upperLimits = upperLimits;
	}

	public Map<Integer, Integer[]> getAllRatios() {
		return allRatios;
	}

	public void setAllRatios(Map<Integer, Integer[]> allRatios) {
		this.allRatios = allRatios;
	}
}
