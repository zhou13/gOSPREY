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
// AminoAcidTemplates.java
//
//  Version:           2.1 beta
//
//
// authors:
//    initials    name            organization                email
//   ---------   --------------  ------------------------    ------------------------------
//     RHL        Ryan Lilien     Dartmouth College           ryan.lilien@dartmouth.edu
//     KER        Kyle E. Roberts       Duke University         ker17@duke.edu
//     PGC        Pablo Gainza C.       Duke University         pablo.gainza@duke.edu
////////////////////////////////////////////////////////////////////////////////////////////

/**
 * Writtten by Ryan Lilien (2001-2004)
 *
 */

import java.io.InputStream;
import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.io.FileInputStream;
import java.util.*;
import java.io.Serializable;

/**
 *
 * This class reads from three data files specifying amino acid residues, N-terminal amino acid residues,
 * and C-terminal amino acid residues. Information read includes element type, atom type, limited connectivity,
 * and partial charge. By matching these amino acid templates to actual residues in a molecule, the corresponding
 * template atom types and partial charges can be assigned to the matched residues.
 *
 */
public class AminoAcidTemplates implements Serializable {

    String aaFilename = "";
    String aaNTFilename = "";
    String aaCTFilename = "";
    String aaCoordsFilename = "all_amino_coords.in";

    Residue aaResidues[];   // arrays of residues
    Residue aaNTResidues[];
    Residue aaCTResidues[];
    int numAAs, numAANTs, numAACTs;

    // Constructor which reads amino acid information from the template files
    AminoAcidTemplates() throws Exception {

        switch(EnvironmentVars.forcefld) {
        case AMBER:
            //KER: This is for the standard amber parameters
            aaFilename =  "all_amino94.in";
            aaNTFilename =  "all_aminont94.in";
            aaCTFilename =  "all_aminoct94.in";
            break;
        case CHARMM22:
            //KER: This if for using the charmm22 parameters:
            aaFilename = "all_amino_charmm22.txt";
            aaNTFilename = "all_amino_charmm22_nt.txt";
            aaCTFilename = "all_amino_charmm22_ct.txt";
            break;
        case CHARMM19NEUTRAL:
            //KER: This is for CHARMM19 parameters:
            aaFilename =  "all_amino_charmm19_neutral.in";
            aaNTFilename =  "all_amino_charmm19_neutral_nt.in";
            aaCTFilename =  "all_amino_charmm19_neutral_ct.in";
            break;
        case CHARMM19:
            aaFilename =  "all_amino_charmm19.in";
            aaNTFilename =  "all_amino_charmm19_nt.in";
            aaCTFilename =  "all_amino_charmm19_ct.in";
            break;
        default:
            System.out.println("FORCEFIELD not recognized...Exiting");
            System.exit(0);
        }


        // HANDLE THE NORMAL AAs
        FileInputStream is = new FileInputStream( EnvironmentVars.getDataDir().concat(aaFilename) );
        BufferedReader bufread = new BufferedReader(new InputStreamReader(is));
        String curLine = null, tmpName = null;
        int dumPresent = 0;
        numAAs = numAANTs = numAACTs = 0;
        aaResidues = new Residue[30];
        aaNTResidues = new Residue[30];
        aaCTResidues = new Residue[30];

        // Skip over first 2 lines of header info
        curLine = bufread.readLine();
        curLine = bufread.readLine();

        while( curLine != null ) {
            // Skip over first line which is the long amino acid name
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
                //KER: The first atom is bonded to a dummy atom so we can't include that
                //KER: in the bond list, so check atom is >= 0
                int atomBondedTo = ((new Integer(getToken(curLine,5))).intValue())-dumPresent;
                if(atomBondedTo >=0)
                    at.addBond(atomBondedTo);
                newRes.addAtom(at);  // add atom
                curLine = bufread.readLine();  // read next line
            }

            //KER: Read LOOP data if any
            curLine = bufread.readLine();
            if (curLine.length() >= 4) {
                if(getToken(curLine, 1).equalsIgnoreCase("LOOP")) {
                    curLine = bufread.readLine();
                    while(!getToken(curLine,2).equals("")) {
                        //find atom1
                        for(Atom a : newRes.atom) {
                            if(a.name.equalsIgnoreCase(getToken(curLine,1))) {
                                //find atom2
                                for(Atom b : newRes.atom) {
                                    if(b.name.equalsIgnoreCase(getToken(curLine,2))) {
                                        a.addBond(b.residueAtomNumber);
                                    }
                                }
                            }
                        }
                        curLine = bufread.readLine();
                    }
                }
            }

            // Eventually we might want to be able to handle the improper
            //  torsions listed here

            // Add the residue to the array
            aaResidues[numAAs++] = newRes;
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

        //Read in the atom coords for the previous residues
        is = new FileInputStream( EnvironmentVars.getDataDir().concat(aaCoordsFilename) );
        bufread = new BufferedReader(new InputStreamReader(is));
        curLine = null;
        tmpName = null;
        int tmpCtr = 0;

        curLine = bufread.readLine();
        while ( curLine.startsWith("#") ) {
            curLine = bufread.readLine();
        }
        boolean foundRes = false;
        boolean foundAtom = false;
        while (curLine != null ) {
            String resName = getToken(curLine,1);
            int numAtoms = new Integer(getToken(curLine,2));
            foundRes = false;
            for(int j = 0; j<numAAs; j++) {
                Residue r = aaResidues[j];
                //the coord template matches the forcefield template
                if(r.name.equalsIgnoreCase(resName)) {
                    foundRes = true;
                    for(int i=0; i<numAtoms; i++) {
                        curLine = bufread.readLine();
                        //Find the current atom in the residue
                        foundAtom = false;
                        for(Atom a:r.atom) {
                            if(a.name.equalsIgnoreCase(getToken(curLine,1))) {
                                foundAtom = true;
                                a.setCoords(new Float(getToken(curLine,2)), new Float(getToken(curLine,3)), new Float(getToken(curLine,4)));
                                break;
                            }
                        }
                        if(!foundAtom) {
                            System.out.println("Residue coord and forcefield templates did not match up.");
                            System.exit(0);
                        }
                    }
                    break;
                }
            }
            //If we didn't find a match we need to get rid of those
            //lines from the file
            if(!foundRes) {
                for(int i=0; i<numAtoms; i++) {
                    curLine=bufread.readLine();
                }
            }
            //Read to catch the ENDRES line and then
            //get the start of the next AA
            curLine = bufread.readLine();
            curLine = bufread.readLine();
        }
        bufread.close();

        // HANDLE THE N-terminal AAs
        is = new FileInputStream( EnvironmentVars.getDataDir().concat(aaNTFilename) );
        bufread = new BufferedReader(new InputStreamReader(is));
        curLine = null;
        tmpName = null;

        // Skip over first 2 lines of header info
        curLine = bufread.readLine();
        curLine = bufread.readLine();

        while( curLine != null ) {
            // Skip over first line which is the long amino acid name
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
            aaNTResidues[numAANTs++] = newRes;
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

        // HANDLE THE C-terminal AAs
        is = new FileInputStream( EnvironmentVars.getDataDir().concat(aaCTFilename) );
        bufread = new BufferedReader(new InputStreamReader(is));
        curLine = null;
        tmpName = null;

        // Skip over first 2 lines of header info
        curLine = bufread.readLine();
        curLine = bufread.readLine();

        while( curLine != null ) {
            // Skip over first line which is the long amino acid name
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
            aaCTResidues[numAACTs++] = newRes;
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
