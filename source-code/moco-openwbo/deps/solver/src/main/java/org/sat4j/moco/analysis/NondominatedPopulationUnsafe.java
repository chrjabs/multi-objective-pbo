/*   João Cortes, Ines Lynce and Vasco Manquinho - MOCO solver
 *******************************************************************************/

package org.sat4j.moco.analysis;

import org.moeaframework.core.NondominatedPopulation;
import org.moeaframework.core.Population;
import org.moeaframework.core.Solution;
import org.moeaframework.core.Variable;
import org.moeaframework.core.variable.EncodingUtils;
import org.sat4j.moco.pb.PBSolver;
import org.sat4j.moco.problem.Instance;
import org.sat4j.moco.util.Clock;
import org.sat4j.moco.util.Log;

/**
 * 
 * 
 * 
 * @author João Cortes
 */

public class NondominatedPopulationUnsafe extends NondominatedPopulation {

    public void addUnsafe(Solution solution){
	this.forceAddWithoutCheck(solution);
    }
}
