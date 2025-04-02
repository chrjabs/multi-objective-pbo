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

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.zip.GZIPInputStream;

import org.sat4j.moco.parsing.OPBReader;

/**
 * Class with several static utility methods for input/output operations.
 * @author Miguel Terra-Neves
 */
public class IOUtils {

    /**
     * This system's line separator.
     */
    public static final String NEWLINE = System.getProperty("line.separator");
    
    /**
     * Makes a fresh temporary file.
     * @param prefix The prefix for the temporary file. 
     * @param suffix The suffix for the temporary file.
     * @param unique True if the file should have a unique name, false otherwise.
     * @param tmp_dir The directory in which the temporary file should be created.
     * @return The temporary file.
     */
    public static File mkTempFile(String prefix, String suffix, boolean unique, File tmp_dir) {
        if (unique) {
            prefix += "_" + Long.toString(System.currentTimeMillis());
        }
        File tmp_file = null;
        try {
            if (tmp_dir != null) {
                tmp_file = File.createTempFile(prefix, suffix, tmp_dir);
            }
            else {
                tmp_file = File.createTempFile(prefix, suffix);
            }
        }
        catch (IOException e) {
            throw new RuntimeException("Failed to create temporary file " + prefix + suffix, e);
        }
        tmp_file.deleteOnExit();
        return tmp_file;
    }
    
    /**
     * Makes a fresh temporary file in the default temporary directory.
     * @param prefix The prefix for the temporary file. 
     * @param suffix The suffix for the temporary file.
     * @param unique True if the file should have a unique name, false otherwise.
     * @return The temporary file.
     */
    public static File mkTempFile(String prefix, String suffix, boolean unique) {
        return mkTempFile(prefix, suffix, unique, null);
    }
    
    /**
     * Builds a reader for parsing a MOCO instance from a given file.
     * @param fname The file path.
     * @return The MOCO reader.
     * @throws IOException if an error occurs building the reader.
     */
    public static OPBReader mkFileReader(String fname) throws IOException {
        InputStream in = new FileInputStream(fname);
        if (fname.endsWith(".gz")) {
            in = new GZIPInputStream(in);
        }
        return new OPBReader(in);
    }
    
}
