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
package org.sat4j.moco.analysis;
import org.moeaframework.core.Settings;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Arrays;
import java.util.Map.Entry;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.moeaframework.Analyzer.AlgorithmResult;
import org.moeaframework.Analyzer.AnalyzerResults;
import org.moeaframework.Analyzer.IndicatorResult;
import org.moeaframework.core.NondominatedPopulation;
import org.moeaframework.core.Population;
import org.moeaframework.core.PopulationIO;
import org.sat4j.core.VecInt;
import org.sat4j.moco.parsing.OPBReader;
import org.sat4j.moco.problem.Instance;
import org.sat4j.moco.util.IOUtils;
import org.sat4j.moco.util.Log;
import org.sat4j.specs.IVecInt;

import org.moeaframework.core.Solution;
import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Multimap;

/**
 * Class for performing performance analysis based on the results produced by the MOCO solver.
 * Contains a {@link #main(String[])} method, which can be used to perform analysis based on log outputs.
 * @author Miguel Terra-Neves
 */
public class Analyzer {

    /**
     * The MOCO instance being analyzed.
     */
    private Instance moco = null;
    
    /**
     * A MOEA framework representation of the MOCO instance.
     */
    private MOCOProblem problem = null;
    
    /**
     * Stores the container objects with the results of different executions of the solver.
     * A multimap is used because there may be multiple results for the same algorithm (multiple executions
     * of the same algorithm with different seeds).
     */
    private Multimap<String, Result> dataset = ArrayListMultimap.create();
    
    /**
     * Creates an instance of the analyzer for a given MOCO instance.
     * @param m The instance.
     */
    public Analyzer(Instance m) {
        this.moco = m;
        this.problem = new MOCOProblem(m);
	// Settings.PROPERTIES.setString("org.moeaframework.configuration.EPS", "0.0E-1");
	// System.out.println("EPS is:");
	// System.out.println(Settings.PROPERTIES.getDouble("org.moeaframework.configuration.EPS", 0));

    }
    
    /**
     * Adds a new result for a given algorithm key.
     * @param key The key.
     * @param r The result container.
     */
    public void addResult(String key, Result r) { this.dataset.put(key, r); }
    
    /**
     * Imports the result of a solver execution from a given log file.
     * @param key The algorithm key.
     * @param path The log file path.
     * @throws IOException if an error occurs reading from the log file.
     */
    public void importResult(String key, String path) throws IOException {
        Log.comment(2, "importing result from " + path);
        OutputReader reader = new OutputReader(new FileReader(path), this.moco);
        Result r = reader.readResult();
        this.dataset.put(key, r);
        reader.close();
        Log.comment(2, "import successful");
    }
    
    /**
     * Builds the reference set of the MOCO instance to be used in computation of performance indicators
     * that require a reference set (e.g. inverted generational distance).
     * @return The reference set.
     */
    private ReferenceSet mkRefSet() {
        Log.comment(3, "{ Analyzer.mkRefSet");
        ReferenceSet r = new ReferenceSet(this.moco);
        r.addAll(this.dataset.values());
        Log.comment(3, "out Analyzer.mkRefSet");
        return r;
    }
    /**
     * Builds the reference set of the MOCO instance to be used in computation of performance indicators
     * that require a reference set (e.g. inverted generational distance).
     * @return The reference set.
     */
    private void projectDataSet(IVecInt consObjectives) {
	for( String key: this.dataset.keySet().toArray(new String[0])) {
	    Result[] rs = this.dataset.get(key).toArray(Result[]::new);
	    Result[] newRs = new Result[rs.length];
	    int ir = 0;
		for(Result result: rs){
		    newRs[ir] = new Result(this.moco);
		    Iterator<Solution> it = result.getSolutions().iterator();
		    // project the solutions, and build a collection of
		    // results for the same key
		    while(it.hasNext()){
			Solution sol = it.next();
			Solution newSol = this.problem.newSolution();
			// build the new solution, and add it to the new result
			double[] objectives = new double[this.moco.nObjs()];
			int i = 0;
			for(int j  = 0, m = sol.getNumberOfObjectives(); j < m; j++){
			    if(!consObjectives.contains(j)){
				objectives[i] = sol.getObjective(j);
				i++;
			    }
			}

			newSol.setObjectives(objectives);
			newRs[ir].addSolution(newSol);
		    }
		}
	    // replace the dataset with the projected values
	    this.dataset.replaceValues(key, Arrays.asList(newRs));
	}
}

    
    /**
     * Writes the cost vectors of the solutions in a given population of solutions to a given file.
     * @param p The solutions.
     * @param f The file.
     */
    private void dumpObjVals(Population p, File f) {
        Log.comment(3, "{ Analyzer.dumpObjVals");
        try {
            PopulationIO.writeObjectives(f, p);
        }
        catch (IOException e) {
            throw new RuntimeException("Failed to write reference set to temporary file", e);
        }
        Log.comment(3, "out Analyzer.dumpObjVals");
    }
    
    /**
     * Logs the number of solutions in each of the results produced by a given algorithm.
     * @param alg The algorithm key.
     * @param rs The result collection.
     */
    private void logPopSizes(String alg, Collection<Result> rs) {
        String line = ":algorithm " + alg + " :pop-sizes [";
        for (Iterator<Result> it = rs.iterator(); it.hasNext();) {
            line += it.next().nSolutions() + (it.hasNext() ? ", " : "]");
        }
        Log.comment(1, line);
    }
    
    /**
     * Adds the contents of the dataset to a MOEA framework analyzer object.
     * @param analyzer The MOEA framework analyzer.
     */
    private void addDataset(org.moeaframework.Analyzer analyzer) {
        Log.comment(3, "{ Analyzer.addDataset");
        for (Iterator<String> it = this.dataset.keySet().iterator(); it.hasNext();) {
            String key = it.next();
            Collection<Result> rs = this.dataset.get(key);
            logPopSizes(key, rs);
            List<NondominatedPopulation> pops = new ArrayList<NondominatedPopulation>(rs.size());
            for (Iterator<Result> r_it = rs.iterator(); r_it.hasNext();) {
                pops.add(r_it.next().getSolutions());
            }
            analyzer.addAll(key, pops);
        }
        Log.comment(3, "out Analyzer.addDataset");
    }
    
    /**
     * Builds the MOEA framework analyzer object with the given reference set.
     * @param rs The reference set.
     * @return The MOEA framework analyzer.
     */
    private org.moeaframework.Analyzer buildMOEAAnalyzer(ReferenceSet rs) {
        Log.comment(3, "{ Analyzer.buildMOEAAnalyzer");
        // MOEA Framework's Analyzer requires the reference set to be specified through a file
        File r_file = IOUtils.mkTempFile("ref_set", ".pop", true);
        dumpObjVals(rs.getSolutions(), r_file);
        org.moeaframework.Analyzer analyzer = new org.moeaframework.Analyzer();
	
        analyzer = analyzer.withProblem(this.problem)
	    .withReferencePoint(rs.getRefPoint().getObjectives())
	    .withIdealPoint(rs.getIdealPoint().getObjectives())
	    // .includeInvertedGenerationalDistance()
	    .includeHypervolume()
	    .includeContribution()
	    // .includeInvertedGenerationalDistance()
	    // .showIndividualValues()
	    // .showStatisticalSignificance()
	    .withReferenceSet(r_file);
        addDataset(analyzer);
        Log.comment(3, "out Analyzer.buildMOEAAnalyzer");
        return analyzer;
    }
    /**
       Check for objective functions that do not change over the reference set
    */
    private IVecInt checkConstantObjectives(ReferenceSet rs) {
	Population p = rs.getSolutions();
	IVecInt consObjectives = new VecInt(new int[]{});
	// Objcectives need to be ordered in a descendent manner, to
	// make sure I can delete them one by one.
	for(int i = this.moco.nObjs()-1; i >= 0; i--)
	    consObjectives.push(i);
	for(int i = 1, n = p.size(); i < n; i++) {
	    for(int j  = 0; j < consObjectives.size(); j++)
		if(p.get(i).getObjective(consObjectives.get(j)) != p.get(0).getObjective(consObjectives.get(j))){
		    consObjectives.remove(consObjectives.get(j));
		    j--;
		}
	}
	return consObjectives;
    }

    
    /**
     project the reference set
     */
    private void projectReferenceSet(ReferenceSet rs, IVecInt consObjectives){
	int nObj = this.moco.nObjs() - consObjectives.size();
	Population p = rs.getSolutions();
	Population p1 = new NondominatedPopulation();
	for(int i = 0, n = p.size(); i < n; i++) {
	    double[] objectives = new double[nObj];
	    for(int j  = 0, m = nObj; j < m; j++){
		objectives[j] = p.get(i).getObjective(j);
		
	    }
}

};

    /**
     projects the problem
    */
    private void projectInstance(Instance instance, IVecInt consObjectives){
	int nObj = consObjectives.size();
	for(int j  = 0; j < nObj; j++)
	    instance.removeObj(consObjectives.get(j));
	

};

    /**
     projects tho data
    */
    private void projectInstance(Multimap<String, Result> dataset, IVecInt consObjectives){
	int nObj = consObjectives.size();
	for(Entry<String, Result> entryA: this.dataset.entries()){
	    NondominatedPopulation B = entryA.getValue().solutions;	
	    
	}
};


    /**
     *compare populations
     */
    private int compareAlgorithms(){

    	for(Entry<String, Result> entryA: this.dataset.entries()) {
	    NondominatedPopulation A = entryA.getValue().solutions;
	    for(Entry<String, Result> entryB: this.dataset.entries()) {
		NondominatedPopulation B = entryA.getValue().solutions;
		int iA = 0;
		int iB = 0;
		int [][] interDominance = new int [A.size()][ B.size()];
		Iterator<Solution> itA = A.iterator();
		Iterator<Solution> itB = B.iterator();
		while(itA.hasNext()){
		    iA++;
		    while(itB.hasNext()){
			iB++;
			Solution solA = itA.next();
			Solution solB = itB.next();
			switch(A.getComparator().compare(solA, solB)){
			case -1:
			    interDominance[iA][iB] = -1;
			case 1:
			    interDominance[iA][iB] = 1;
			    break;
			case 0:
			    interDominance[iA][iB] = 0;
			    break;
			}
		    }}
		System.out.print(interDominance);
	    }}
	return 0;
    }

    /**
     *compare populations(simple)
     */
    private int simpleCompareAlgorithms(){
	HashMap<String, Integer> interDominance = new HashMap<String, Integer>();
	for(Entry<String, Result> entryA: this.dataset.entries()) {
	    NondominatedPopulation A = entryA.getValue().solutions;
	    Iterator<Solution> itA = A.iterator();
	    int countA = 0;
	    for(Entry<String, Result> entryB: this.dataset.entries()){
		NondominatedPopulation B = entryA.getValue().solutions;
		Iterator<Solution> itB = B.iterator();
		while(itA.hasNext()){
		    Solution solA = itA.next();
		    while(itB.hasNext()){
			Solution solB = itB.next();
			if(A.getComparator().compare(solA, solB) == 1){
			    countA++;
			}
		    }
		}
		interDominance.put(entryA.getKey() + " < " + entryB.getKey(), countA);
	    }
	}
	String interDominanceLog ="interdominance: ";
	for(Entry<String, Integer> interDominanceEntry:interDominance.entrySet()){
	    interDominanceLog += interDominanceEntry.getKey();
	    interDominanceLog += ", ";
	    interDominanceLog += interDominanceEntry.getValue();
	    interDominanceLog += "| ";

	}
	Log.comment(interDominanceLog);
	return 0;
    }

    /**
     * Analyzes the execution results in the {@link #dataset} and prints the results of the analysis
     * (hypervolumes, inverted generational distances and respective statistical significances).
     */
    public void printAnalysis() {
        Log.comment(3, "{ Analyzer.printAnalysis");
        ReferenceSet ref = mkRefSet();
	IVecInt objectives =  checkConstantObjectives(ref);
	this.problem = new MOCOProblem(this.moco);
	int nObjs = moco.nObjs();
	if(nObjs != objectives.size()){
	    this.projectInstance(this.moco, objectives);
	    this.projectDataSet(objectives);
	}
	// project the ref set, after projecting the data set
	ref = this.mkRefSet();

        if (ref.isEmpty()) {
            Log.comment("0 known solutions, impossible to analyze");
        }
        else {
            org.moeaframework.Analyzer analyzer = buildMOEAAnalyzer(ref);
	    
	    AnalyzerResults results = analyzer.getAnalysis();
            List<String> algs = results.getAlgorithms();
            for (Iterator<String> a_it = algs.iterator(); a_it.hasNext();) {
                String alg = a_it.next();
                AlgorithmResult ar = results.get(alg);
                List<String> indicators = ar.getIndicators();
                for (Iterator<String> i_it = indicators.iterator(); i_it.hasNext();) {
                    String ind = i_it.next();
                    IndicatorResult ir = ar.get(ind);
                    String line = ":algorithm " + alg + " :indicator " + ind + " :values [";
                    double[] vals = ir.getValues();
                    for (int i = 0; i < vals.length; ++i) {
                        line += vals[i] + (i == vals.length-1 ? "] " : ", ");
                    }
                    line += ":min " + ir.getMin() + " :max " + ir.getMax() + " :median " + ir.getMedian();
                    Log.comment(line);
                }
            }
        }
        Log.comment(3, "out Analyzer.printAnalysis");
    }
    
    /**
     * Prints a help message to standard output.
     * @param options An object that represents the application's options.
     */
    private static void printHelpMessage(Options options) {
        HelpFormatter formatter = new HelpFormatter();
        formatter.printHelp("Analyzer instance output-description...", options);
    }
    
    /**
     * Builds an {@link Options} object with the analyzer's configuration parameters to be used for parsing
     * command line options.
     * @return An {@link Options} object to be used by the command line interface.
     */
    private static Options buildOpts() {
        Options options = new Options();
        options.addOption("v", "verbosity", true, "Set the verbosity level (from 0 to 3). Default is 0.");
        return options;
    }
    
    /**
     * Sets the analyzer's verbosity level to the one given in the command line.
     * @param cl The command line object.
     */
    private static void updtOpts(CommandLine cl) {
        int verb_lvl = Integer.parseInt(cl.getOptionValue("v", "0"));
        Log.setVerbosity(verb_lvl);
        Log.comment(":verbosity " + verb_lvl);
    }
    
    /**
     * Reads a MOCO instance from a file given as a command line argument.
     * @param cl The command line object that contains the input file.
     * @return The instance.
     * @throws IOException if an error occurs parsing the instance.
     */
    // TODO: same as CbMOCO.readMOCO; DIMACS will be supported in the future, so such code should be unified
    private static Instance readMOCO(CommandLine cl) throws IOException {
        OPBReader reader = IOUtils.mkFileReader(cl.getArgs()[0]);
        Instance moco = reader.readMOCO();
        reader.close();
        return moco;
    }
    
    /**
     * Parses a set of dataset description entries.
     * A dataset description entry has the following format:
     * '&lt;algorithm key&gt;:&lt;log file path 1&gt;:&lt;log file path 2&gt;:...:&lt;log file path n&gt;'.
     * Multiple log files for the same algorithm key correspond to multiple runs of the same algorithm with
     * different seeds.
     * @param descs An array of dataset description entries.
     * @return A multimap of algorithm keys to log file paths.
     */
    private static Multimap<String, String> parseDescs(String[] descs) {
        Multimap<String, String> paths = ArrayListMultimap.create();
        for (int i = 0; i < descs.length; ++i) {
            String[] tokens = descs[i].split(":");
            String key = tokens[0];
            if (tokens.length > 1) {
                for (int j = 1; j < tokens.length; ++j) {
                    paths.put(key, tokens[j]);
                }
            }
            else {
                System.out.println("WARNING: no output files provided for key " + key);
            }
        }
        return paths;
    }
    
    /**
     * Builds the analyzer for a given MOCO instance and a given dataset.
     * The dataset is described as an array of dataset description entries with the following format:
     * '&lt;algorithm key&gt;:&lt;log file path 1&gt;:&lt;log file path 2&gt;:...:&lt;log file path n&gt;'.
     * Multiple log files for the same algorithm key correspond to multiple runs of the same algorithm with
     * different seeds.
     * @param m The instance.
     * @param descs An array of dataset description entries.
     * @return The analyzer.
     * @throws IOException if an error occurs importing a log file's results.
     */
    private static Analyzer buildAnalyzer(Instance m, String[] descs) throws IOException {
        Log.comment(3, "{ Analyzer.buildAnalyzer");
        Multimap<String, String> dataset_paths = parseDescs(descs);
        Analyzer analyzer = new Analyzer(m);
        for (Iterator<String> it = dataset_paths.keySet().iterator(); it.hasNext();) {
            String key = it.next();
            Collection<String> paths = dataset_paths.get(key);
            for (Iterator<String> p_it = paths.iterator(); p_it.hasNext();) {
                analyzer.importResult(key, p_it.next());
            }
        }
        Log.comment(3, "out Analyzer.buildAnalyzer");
        return analyzer;
    }
    
    /**
     * Retrieves the dataset description entries given as command line arguments.
     * A dataset description entry has the following format:
     * '&lt;algorithm key&gt;:&lt;log file path 1&gt;:&lt;log file path 2&gt;:...:&lt;log file path n&gt;'.
     * Multiple log files for the same algorithm key correspond to multiple runs of the same algorithm with
     * different seeds.
     * @param cl The command line object.
     * @return An array with the dataset description entries.
     */
    private static String[] retrieveDescs(CommandLine cl) {
        String[] descs = new String[cl.getArgs().length-1];
        for (int i = 1; i < cl.getArgs().length; ++i) {
            descs[i-1] = cl.getArgs()[i];
        }
        return descs;
    }
    
    /**
     * The analyzer's entry point.
     * @param args The command line arguments.
     */
    public static void main(String[] args) {
        Options options = buildOpts();
        CommandLineParser cl_parser = new DefaultParser();
        CommandLine cl = null;
        try {
            cl = cl_parser.parse(options, args);
            if (cl.getArgs().length < 2) {
                System.out.println("Illegal number of arguments");
                printHelpMessage(options);
                return;
            }
            updtOpts(cl);
            Instance moco = readMOCO(cl);
            Analyzer analyzer = buildAnalyzer(moco, retrieveDescs(cl));
            analyzer.printAnalysis();
	    analyzer.simpleCompareAlgorithms();
        }
        catch (ParseException e) {
            printHelpMessage(options);
        }
        catch (IOException e) {
            System.out.println("PARSER ERROR!");
            e.printStackTrace();
        }
    }

}
