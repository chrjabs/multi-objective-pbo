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
package org.sat4j.moco.util;

import java.math.BigDecimal;
import java.math.BigInteger;
import java.math.RoundingMode;

import org.sat4j.moco.Params;

/**
 * {@link BigDecimal} wrapper that enforces a fixed maximum scale (number of digits to the right of the
 * decimal point).
 * @author Miguel Terra-Neves
 */
public class Real implements Comparable<Real> {

    /**
     * The rounding mode to use when the scale is not enough.
     */
    private static final RoundingMode ROUNDING_MODE = RoundingMode.HALF_EVEN;
    
    /**
     * A real with value 0.
     */
    public static final Real ZERO = new Real(BigDecimal.ZERO);
    
    /**
     * A real with value 1.
     */
    public static final Real ONE = new Real(BigDecimal.ONE);
    
    /**
     * A real with value 10.
     */
    public static final Real TEN = new Real(BigDecimal.TEN);
    
    /**
     * Maximum scale considered for real numbers.
     */
    private static int scale = 5;

    /**
     * Sets the maximum scale considered for real numbers.
     * @param s The scale.
     */
    public static void setScale(int s) { Real.scale = s; }
    
    /**
     * Sets the maximum scale for real numbers to the one stored in a given set of parameters.
     * @param p The parameters object.
     */
    public static void updtParams(Params p) {
        setScale(p.getScale());
        Log.comment(":decimal-scale " + Real.scale);
    }
    
    /**
     * An instance of the underlying {@link BigDecimal}.
     */
    private BigDecimal val = null;
    
    /**
     * Creates a real from a given {@link BigDecimal}.
     * @param val The {@link BigDecimal} value.
     */
    public Real(BigDecimal val) { this.val = val.setScale(Real.scale, ROUNDING_MODE); }
    
    /**
     * Creates a real from a given {@link BigInteger}.
     * @param val The {@link BigInteger} value.
     */
    public Real(BigInteger val) { this(new BigDecimal(val)); }
    
    /**
     * Creates a real from a given double.
     * @param val The double value.
     */
    public Real(double val) { this(new BigDecimal(val)); }
    
    /**
     * Creates a real from a given integer value.
     * @param val The integer value.
     */
    public Real(int val) { this(new BigDecimal(val)); }
    
    /**
     * Creates a real from a given long value.
     * @param val The long value.
     */
    public Real(long val) { this(new BigDecimal(val)); }
    
    /**
     * Creates a real from a given string representation.
     * @param val The string representation.
     * @see BigDecimal#BigDecimal(String)
     */
    public Real(String val) { this(new BigDecimal(val)); }
    
    /**
     * Converts the real to a {@link BigDecimal}.
     * @return The real as a {@link BigDecimal}.
     */
    public BigDecimal asBigDecimal() { return this.val; }
    
    /**
     * Converts the real to a {@link BigInteger}.
     * @return The real as a {@link BigInteger}.
     * @see {@link BigDecimal#toBigInteger()}.
     */
    public BigInteger asBigInteger() { return this.val.toBigInteger(); }

    /**
     * Converts the real to a {@link BigInteger}.
     * An exception is thrown if the real has a non-zero fractional part.
     * @return The real as a {@link BigInteger}.
     * @see {@link BigDecimal#toBigIntegerExact()}
     */
    public BigInteger asBigIntegerExact() { return this.val.toBigIntegerExact(); }
    
    /**
     * Converts the real to a double.
     * @return The real as a double.
     * @see {@link BigDecimal#doubleValue()}
     */
    public double asDouble() { return this.val.doubleValue(); }
    
    /**
     * Converts the real to an integer.
     * @return The real as an integer.
     * @see {@link BigDecimal#intValue()}
     */
    public int asInt() { return this.val.intValue(); }
    
    /**
     * Converts the real to an integer.
     * An exception is thrown if the real has a non-zero fractional part.
     * @return The real as an integer.
     * @see {@link BigDecimal#intValueExact()}
     */
    public int asIntExact() { return this.val.intValueExact(); }

    /**
     * Converts the real to a long.
     * @return The real as a long.
     * @see {@link BigDecimal#longValue()}
     */
    public long asLong() { return this.val.longValue(); }
    
    /**
     * Converts the real to a long.
     * An exception is thrown if the real has a non-zero fractional part.
     * @return The real as a long.
     * @see {@link BigDecimal#longValueExact()}
     */
    public long asLongExact() { return this.val.longValueExact(); }
    
    /**
     * Computes the real's absolute value.
     * @return The real's absolute value.
     */
    public Real abs() { return new Real(this.val.abs()); }
    
    /**
     * Computes the result of adding this real to another real.
     * @param augend The addition's augend.
     * @return {@code this + augend}.
     */
    public Real add(Real augend) { return new Real(this.val.add(augend.asBigDecimal())); }
    
    /**
     * Computes the result of subtracting another real from this real.
     * @param subtrahend The subtraction's subtrahend.
     * @return {@code this - subtrahend}.
     */
    public Real subtract(Real subtrahend) { return new Real(this.val.subtract(subtrahend.asBigDecimal())); }
    
    /**
     * Computes the result of multiplying this real by another.
     * @param multiplicand The multiplication's multiplicand.
     * @return {@code this * multiplicand}.
     */
    public Real multiply(Real multiplicand) { return new Real(this.val.multiply(multiplicand.asBigDecimal())); }
    
    /**
     * Computes the result of dividing this real by another.
     * @param divisor The division's divisor.
     * @return {@code this / divisor}.
     */
    public Real divide(Real divisor) {
        return new Real(this.val.divide(divisor.asBigDecimal(), Real.scale, ROUNDING_MODE));
    }
    
    /**
     * Checks if this real is equal to a given object.
     * @param other The other object.
     * @return True if {@code other} is a real equal to this one, false otherwise.
     */
    public boolean equals(Object other) {
        // compareTo is used because BigDecimal's equals method returns false if the BigDecimals have
        // different scales, even if representing the same value
        return other instanceof Real && this.compareTo((Real)other) == 0;
    }
    
    /**
     * Compares this real with another real.
     * @param other The other real.
     * @return 0 if {@code this} is equal to {@code other}, a value greater than 0 if {@code this} is greater
     * than {@code other}, a value smaller than 0 if {@code this} is smaller than {@code other}.
     */
    public int compareTo(Real other) { return this.val.compareTo(other.asBigDecimal()); }
    
    /**
     * Checks if the real is greater than another real.
     * @param other The other real.
     * @return True if {@code this} is greater than {@code other}, false otherwise.
     */
    public boolean greaterThan(Real other) { return this.val.compareTo(other.asBigDecimal()) > 0; }
    
    /**
     * Checks if the real is greater or equal to another real.
     * @param other The other real.
     * @return True if {@code this} is greater or equal to {@code other}, false otherwise.
     */
    public boolean greaterOrEqual(Real other) { return this.val.compareTo(other.asBigDecimal()) >= 0; }
    
    /**
     * Checks if the real is less than another real.
     * @param other The other real.
     * @return True if {@code this} is less than {@code other}, false otherwise.
     */
    public boolean lessThan(Real other) { return !greaterOrEqual(other); }
    
    /**
     * Checks if the real is less or equal to another real.
     * @param other The other real.
     * @return True if {@code this} is less or equal to {@code other}, false otherwise.
     */
    public boolean lessOrEqual(Real other) { return !greaterThan(other); }
    
    /**
     * Computes the signum function of the real.
     * @return 1, 0 or -1 if {@code this} is positive, 0 or negative respectively.
     */
    public int signum() { return this.val.signum(); }
    
    /**
     * Checks if the real is positive.
     * @return True if {@code this} is positive, false otherwise.
     */
    public boolean isPositive() { return signum() > 0; }
    
    /**
     * Checks if the real is negative.
     * @return True if {@code this} is negative, false otherwise.
     */
    public boolean isNegative() { return signum() < 0; }
    
    /**
     * Computes the number of digits in the real to the right of the decimal point.
     * For example, if {@code this} represents the value {@code 25.023}, then this method returns 3.
     * @return The number of digits in {@code this} to the right of the decimal point.
     */
    public int nDecimals() { return Math.max(this.val.stripTrailingZeros().scale(), 0); }
    
    /**
     * Multiples the real by a given power of ten.
     * @param n The exponent.
     * @return {@code this * 10^n}.
     */
    public Real scaleByPowerOfTen(int n) { return new Real(this.val.scaleByPowerOfTen(n)); }
    
    /**
     * Computes the negation of the real.
     * @return {@code -this}.
     */
    public Real negate() { return new Real(this.val.negate()); }
    
    /**
     * Returns the string representation of the real, using scientific notation if needed.
     * @return The real's string representation.
     * @see {@link BigDecimal#toString()}
     */
    public String toString() { return this.val.toString(); }
    
    /**
     * Returns the string representation of the real without an exponent field.
     * @return The real's string representation without an exponent field.
     * @see {@link BigDecimal#toPlainString()}
     */
    public String toPlainString() { return this.val.toPlainString(); }
    
}
