/*
	This file is part of OSPREY.

	OSPREY Protein Redesign Software Version 2.1 beta
	Copyright (C) 2001-2012 Bruce Donald Lab, Duke University

	OSPREY is free software: you can redistribute it and/or modify
	it under the terms of the GNU Lesser General Public License as
	published by the Free Software Foundation, either version 3 of
	the License, or (at your option) any later version.

	OSPREY is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
	GNU Lesser General Public License for more details.

	You should have received a copy of the GNU Lesser General Public
	License along with this library; if not, see:
	      <http://www.gnu.org/licenses/>.

	There are additional restrictions imposed on the use and distribution
	of this open-source code, including: (A) this header must be included
	in any modification or extension of the code; (B) you are required to
	cite our papers in any publications that use this code. The citation
	for the various different modules of our software, together with a
	complete list of requirements and restrictions are found in the
	document license.pdf enclosed with this distribution.

	Contact Info:
			Bruce Donald
			Duke University
			Department of Computer Science
			Levine Science Research Center (LSRC)
			Durham
			NC 27708-0129
			USA
			e-mail:   www.cs.duke.edu/brd/

	<signature of Bruce Donald>, Mar 1, 2012
	Bruce Donald, Professor of Computer Science
*/
////////////////////////////////////////////////////////////////////////////////////////////
// GenericResidueTemplates.java
//
//  Version:           2.1 beta
//
//
// authors:
//    initials     name                   organization                email
//   ---------    --------------     ------------------------    ------------------------------
//     RHL         Ryan Lilien         Dartmouth College           ryan.lilien@dartmouth.edu
//	   ISG		   Ivelin Georgiev	   Duke University			   ivelin.georgiev@duke.edu
//
////////////////////////////////////////////////////////////////////////////////////////////

/**
 * Written by Ryan Lilien (2001-2004);
 * 	minor modifications by Ivelin Georgiev (2004-2009)
 *
 */

import java.io.InputStream;
import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.io.FileInputStream;
import java.util.*;
import java.io.Serializable;

/**
 * This class reads from a generic residue file that includes element type, AMBER atom type, limited connectivity,
 * and partial charge. This file is analogous to AminoAcidTemplates.java; instead of amino acid parameters, parameters
 * for general compounds and nucleic acids (referred to as 'generic residues') are read.
 * By matching these generic residue templates to actual generic residues in a molecule, the corresponding
 * template atom types and partial charges can be assigned to the matched residues.
 * The format of the input parameter file is similar to the PARM AMBER datafiles, identical to the all_amino94.in.
 */
public class GenericResidueTemplates implements Serializable {

    public static final int MAX_NUM_RES = 50;

    String grFilename = "all_nuc94_and_gr.in";

    Residue grResidues[];   // arrays of residues
    int numGRs;

    GenericResidueTemplates() throws Exception {

        FileInputStream is = new FileInputStream( EnvironmentVars.getDataDir().concat(grFilename) );
        BufferedReader bufread = new BufferedReader(new InputStreamReader(is));
        String curLine = null, tmpName = null;
        int dumPresent = 0;
        numGRs = 0;
        grResidues = new Residue[MAX_NUM_RES];
        // Skip over first 2 lines of header info
        curLine = bufread.readLine();
        curLine = bufread.readLine();

        while( curLine != null ) {
            // Skip over first line which is the long residue name
            curLine = bufread.readLine();
            if (curLine.length() >= 4)
                if (curLine.substring(0,4).equalsIgnoreCase("stop")) {
                    curLine = bufread.readLine();
                    continue;
                }
            // Skip blank line
            curLine = bufread.readLine();
            // The next line contains the 3 letter amino acid name
            curLine = bufread.readLine();
            tmpName = getToken(curLine,1);
            Residue newRes = new Residue();
            newRes.name = tmpName;
            // Skip next 2 lines
            curLine = bufread.readLine();
            curLine = bufread.readLine();
            // Now we're into the section with atoms
            curLine = bufread.readLine();
            // Skip the dummy atoms
            dumPresent = 0;
            while (getToken(curLine,2).equalsIgnoreCase("DUMM")) {
                dumPresent++;
                curLine = bufread.readLine();
            }
            dumPresent++; // to adjust for 0-based
            while (!getToken(curLine,2).equals("")) {
                Atom at = new Atom();
                tmpName = getToken(curLine,2);
                at.changeType(tmpName);
                at.forceFieldType = getToken(curLine,3);
                at.charge = (float) (new Float(getToken(curLine,11)).floatValue());
                at.isBBatom = at.setIsBBatom();
                at.addBond(((new Integer(getToken(curLine,5))).intValue())-dumPresent);
                newRes.addAtom(at);  // add atom
                curLine = bufread.readLine();  // read next line
            }
            // Eventually we might want to be able to handle the improper
            //  torsions listed here

            // Add the residue to the array
            grResidues[numGRs++] = newRes;
            // Read until the end of the residue
            boolean atDone = false;
            if (curLine.length() >= 4)
                atDone = curLine.substring(0,4).equalsIgnoreCase("done");
            else
                atDone = false;
            while (!atDone) {
                curLine = bufread.readLine();
                if (curLine.length() >= 4)
                    atDone = curLine.substring(0,4).equalsIgnoreCase("done");
            }
        }
        bufread.close();
    }

    /******************************/
    // This function returns the xth token in string s
    private String getToken(String s, int x) {

        int curNum = 1;
        StringTokenizer st = new StringTokenizer(s," ,;\t\n\r\f");

        while (curNum < x) {
            curNum++;
            if (st.hasMoreTokens())
                st.nextToken();
            else {
                return(new String(""));
            }
        }

        if (st.hasMoreTokens())
            return(st.nextToken());
        return(new String(""));

    } // end getToken

}
