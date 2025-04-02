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
package org.sat4j.moco.parsing;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Reader;

import org.sat4j.core.Vec;
import org.sat4j.core.VecInt;
import org.sat4j.moco.pb.PBConstr;
import org.sat4j.moco.pb.PBFactory;
import org.sat4j.moco.problem.DivObj;
import org.sat4j.moco.problem.LinearObj;
import org.sat4j.moco.problem.Instance;
import org.sat4j.moco.util.Log;
import org.sat4j.moco.util.Real;
import org.sat4j.specs.IVec;
import org.sat4j.specs.IVecInt;

/**
 * Class for parsing MOCO instances in extended OPB format.
 * The extended OPB format supports multiple objectives and objectives expressed as sums of divisions.
 * @author Miguel Terra-Neves
 */
public class OPBReader {
    
    /**
     * The reader object from which to parse the MOCO instance.
     */
    private BufferedReader reader = null;
    
    /**
     * Stores the MOCO instance.
     */
    private Instance moco = null;
    
    /**
     * Creates an OPB format reader object that parses a MOCO instance from a given {@link Reader} object.
     * @param reader The reader object from which to parse the instance.
     */
    public OPBReader(Reader reader) { this.reader = new BufferedReader(reader); }
    
    /**
     * Creates an OPB format reader object that parse a MOCO instance from a given {@link InputStream} object.
     * @param stream The input stream from which to parse the instance.
     */
    public OPBReader(InputStream stream) { this(new InputStreamReader(stream)); }
    
    /**
     * Reads a MOCO instance from the {@link Reader} object provided in {@link #OPBReader(Reader)} or the
     * {@link InputStream} object provided in {@link #OPBReader(InputStream)}.
     * @return The MOCO instance.
     * @throws IOException if an error occurs reading the instance.
     */
    public Instance readMOCO() throws IOException {
        Log.comment(3, "in OPBReader.readMOCO");
        this.moco = new Instance();
        int lineno = 1;
        for (String line = this.reader.readLine(); line != null; line = this.reader.readLine(), lineno++) {
            line = line.trim();
            if (line.isEmpty() || line.startsWith("*")) continue;
            parseLine(line, lineno, line.startsWith("m"));
        }
        Log.comment(1, ":nvars " + this.moco.nVars() + " " +
                       ":constraints " + this.moco.nConstrs() + " " +
                       ":objectives " + this.moco.nObjs());
        Log.comment(3, "out OPBReader.readMOCO");
        return this.moco;
    }
    
    /**
     * Parses a line of the MOCO instance.
     * Updates the MOCO instance stored in {@link #moco}.
     * @param line The line.
     * @param lineno The line's number.
     * @param is_obj True if the line is an objective function, false otherwise.
     */
    private void parseLine(String line, int lineno, boolean is_obj) {
        String[] tokens = line.trim().split("\\s+");
        if (is_obj && !(tokens[0].equals("min:") || tokens[0].equals("max:"))) {
            throw new ParserException(lineno, "unknown objective type " + tokens[0]);
        }
        IVecInt lits = new VecInt();
        IVec<Real> coeffs = new Vec<Real>();
        IVec<IVecInt> num_lits = new Vec<IVecInt>();
        IVec<IVec<Real>> num_coeffs = new Vec<IVec<Real>>();
        IVec<IVecInt> den_lits = new Vec<IVecInt>();
        IVec<IVec<Real>> den_coeffs = new Vec<IVec<Real>>();
        Real rhs = null;
        String op = null;
        boolean neg = is_obj && tokens[0].equals("max:"), is_div_term = false;
        for (int i = is_obj ? 1 : 0; i < tokens.length; ++i) {
            String tok = tokens[i];
            if (tok.endsWith(";")) {
                tok = tok.substring(0, tok.length()-1);
            }
            if (!tok.startsWith(")") && tok.endsWith(")")) {
                tok = tok.substring(0, tok.length()-1);
                tokens[i--] = ")";
            }
            if (tok.startsWith("(")) {
                if (!is_obj) Log.comment("WARNING: ignoring '(' within constraint line");
                if (neg) throw new ParserException(lineno, "maximization of division objectives is not supported");
                if (is_div_term) throw new ParserException(lineno, "'(' within a division term");
                if (lits.size() > 0) throw new ParserException(lineno, "linear objective terms among division objective");
                is_div_term = true;
                tok = tok.substring(1);
            }
            if (tok.isEmpty()) continue;
            if (tok.equals(")") && !is_obj) {
                Log.comment("WARNING: ignoring ')' within constraint line");
                continue;
            }
            if (tok.equals(")") && !is_div_term) throw new ParserException(lineno, "')' without matching '('");
            if (tok.startsWith("/")) throw new ParserException(lineno, "missplaced '/'");
            if (tok.equals(")")) {
                for (++i; i < tokens.length && tokens[i].isEmpty(); ++i);
                if (i == tokens.length || tokens[i].equals(";") || tokens[i].equals("+")) {
                    den_lits.push(lits);
                    den_coeffs.push(coeffs);
                    if (num_lits.size() < den_lits.size()) {
                        throw new ParserException(lineno, "division missing numerator");
                    }
                }
                else if (tokens[i].equals("/")) {
                    num_lits.push(lits);
                    num_coeffs.push(coeffs);
                    if (den_lits.size() < num_lits.size()-1) {
                        throw new ParserException(lineno, "division missing denominator");
                    }
                }
                else {
                    throw new ParserException(lineno, "')' followed by unknown operator " + tokens[i]);
                }
                lits = new VecInt();
                coeffs = new Vec<Real>();
                is_div_term = false;
            }
            else if (tok.startsWith("x") || tok.startsWith("~")) {
                lits.push(parseLit(tok, lineno));
                if (lits.size() > coeffs.size()) {
                    throw new ParserException(lineno, "variable " + tok + " without coefficient");
                }
            }
            else if (isOp(tok)) {
                op = tok;       // TODO: support glued op and rhs
            }
            else {
                if (op == null && coeffs.size() > lits.size()) {
                    throw new ParserException(lineno, "coefficient without variable");
                }
                Real c = parseCoeff(tok, lineno);
                if (op == null) {
                    coeffs.push(neg ? c.negate() : c);
                }
                else {
                    assert(!neg);
                    rhs = c;
                }
            }
        }
        assert(lits.size() == coeffs.size());
        if (lits.size() > 0 && num_lits.size() > 0) {
            throw new ParserException(lineno, "linear objective terms among division objective");
        }
        else if (lits.size() == 0 && (!is_obj || num_lits.size() == 0)) {
            throw new ParserException(lineno, "empty " + (is_obj ? "objective" : "constraint"));
        }
        else if (is_obj && op != null) {
            throw new ParserException(lineno, "operator " + op + " in objective");
        }
        else if (!is_obj && op == null) {
            throw new ParserException(lineno, "missing operator in constraint");
        }
        else if (!is_obj && rhs == null) {
            throw new ParserException(lineno, "missing right-hand side in constraint");
        }
        if (is_obj && num_lits.size() > 0) {
            this.moco.addObj(new DivObj(num_lits, num_coeffs, den_lits, den_coeffs));
        }
        else if (is_obj) {
            this.moco.addObj(new LinearObj(lits, coeffs));
        }
        else {
            this.moco.addConstr(mkConstr(lits, coeffs, rhs, op, lineno));
        }
    }
    
    /**
     * Checks if a given token is the relational operator of a constraint (e.g. ">=").
     * @param tok The token.
     * @return True if token is a relational operator, false otherwise.
     */
    private boolean isOp(String tok) { return tok.equals(">=") || tok.equals("=") || tok.equals("<="); }
    
    /**
     * Parses a literal from a given token.
     * @param tok The token.
     * @param lineno The number of the line the token was extracted from.
     * @return The literal parsed from {@code tok}.
     */
    private int parseLit(String tok, int lineno) {
        assert(tok.startsWith("x") || tok.startsWith("~"));
        boolean neg = false;
        if (tok.startsWith("~")) {
            neg = true;
            tok = tok.substring(2);
        }
        else {
            tok = tok.substring(1);
        }
        try {
            return neg ? -Integer.parseInt(tok) : Integer.parseInt(tok);
        }
        catch (NumberFormatException e) {
            throw new ParserException(lineno, "invalid variable id " + tok);
        }
    }
    
    /**
     * Parse a coefficient from a given token.
     * @param tok The token.
     * @param lineno The number of the line the token was extracted from.
     * @return The coefficient parsed from {@code tok}.
     */
    private Real parseCoeff(String tok, int lineno) {
        try {
            return new Real(tok);
        }
        catch (NumberFormatException e) {
            throw new ParserException(lineno, "invalid coefficient " + tok);
        }
    }
    
    /**
     * Creates a PB constraint object.
     * @param lits The literals in the left-hand side of the constraint.
     * @param coeffs The coefficients in the right-hand side of the constraint.
     * @param rhs The right-hand side of the constraint.
     * @param op The string representation of the constraint's relational operator.
     * @param lineno The number of the line the constraint was extracted from.
     * @return The constraint object.
     */
    private PBConstr mkConstr(IVecInt lits, IVec<Real> coeffs, Real rhs, String op, int lineno) {
        PBConstr c = PBFactory.instance().mkConstr(op, lits, coeffs, rhs);
        if (c == null) throw new ParserException(lineno, "unknown operator " + op);
        return c;
    }
    
    /**
     * Frees the I/O resources associated with the OPB reader.
     * @throws IOException if an error occurs closing the reader.
     */
    public void close() throws IOException { this.reader.close(); }
    
}
