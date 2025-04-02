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
package org.sat4j.moco.pb;

import java.math.BigInteger;

/**
 * An identifier for constraints added to a solver that supports constraint removal.
 * These solvers return a {@code ConstrID} object when adding a removable constraint, that can then be used
 * to remove it.
 * Contraint ids can also be used as timestamps, since these can be compared and are created in increasing
 * order by {@link #mkFresh()}.
 * @author Miguel Terra-Neves
 */
public class ConstrID implements Comparable<ConstrID> {

    /**
     * Counter used to generate unique constraint ids.
     */
    private static BigInteger id_gen = BigInteger.ZERO;

    /**
     * The {@code BigInteger} representation of the constraint id.
     */
    private BigInteger id = null;
    
    /**
     * Creates an instance of a constraint id.
     * @param id {@code BigInteger} representation of the id.
     */
    private ConstrID(BigInteger id) { this.id = id; }
    
    /**
     * Equals comparator method for constraint ids.
     * Checks if the constraint id and a given object are equal.
     * @param other The object to be compared with the constraint id.
     * @return True if {@code other} is a constraint id and references the same constraint as this
     * constraint id, false otherwise.
     */
    @Override
    public boolean equals(Object other) {
        if (other instanceof ConstrID) {
            ConstrID other_id = (ConstrID)other;
            return other_id.id.equals(id);
        }
        return false;
    }

    /**
     * Implementation of the {@code compareTo} method of the {@code Comparable} interface for
     * constraint ids.
     * Checks if the constraint id is smaller or larger than another given constraint id.
     * @param other The constraint id object to be compared with this constraint id.
     * @return An integer smaller, larger or equal to 0 if this constraint id is smaller, larger or
     * equal to {@code other} respectively.
     */
    public int compareTo(ConstrID other) { return id.compareTo(other.id); }
    
    /**
     * Produces and returns a string representation of the constraint id.
     * @return The string representation of the constraint id.
     */
    @Override
    public String toString() { return id.toString(); }
    
    /**
     * Generates a fresh unique constraint id, not equal to any constraint id generated previously.
     * @return A fresh unique constraint id.
     */
    public static ConstrID mkFresh() {
        ConstrID new_id = new ConstrID(id_gen);
        id_gen = id_gen.add(BigInteger.ONE);
        return new_id;
    }
    
}
