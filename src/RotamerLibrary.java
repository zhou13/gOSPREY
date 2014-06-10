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
// RotamerLibrary.java
//
//  Version:           2.1 beta
//
//
// authors:
//    initials    name                 organization                email
//   ---------   -----------------    ------------------------    ----------------------------
//     RHL        Ryan Lilien          Dartmouth College           ryan.lilien@dartmouth.edu
//	   ISG		  Ivelin Georgiev	   Duke University			   ivelin.georgiev@duke.edu
//	  KER        Kyle E. Roberts       Duke University         ker17@duke.edu
//    PGC        Pablo Gainza C.       Duke University         pablo.gainza@duke.edu
///////////////////////////////////////////////////////////////////////////////////////////////

/**
 * Written by Ryan Lilien (2001-2004) and Ivelin Georgiev (2004-2009)
 *
 */

import java.io.*;
import java.util.*;

/**
 * This class implements a rotamer library reader. It reads from an input file that contains rotamer information
 * for amino acid types or other generic residues.
 */
public class RotamerLibrary implements Serializable {

    // If the debug flag is set to true then additional debug statements are
    //  printed to standard out.
    public static final boolean debug = false;

    String aaNames[];       // Array of amino-acid names
    private int numAAallowed;              // Number of AA types read
    private int numDihedrals[];     // Number of dihedrals per AA
    private int numRotamers[];      // Number of rotamers per AA

    private String dihedralAtomNames[][][];  // Names of atoms involved in the dihedrals for each amino acid
    private int rotamerValues[][][];  // Actual angle values for each rotamer for each amino acid
    private float rotamerVolumes[][];	// Volumes of each rotamer for each amino acid

    private int totalNumRotamers; //AAs with 0 rotamers are counted as 1 rotamer
    private int rotamerIndexOffset[] = null; //the rotamer index offset for each amino acid (AAs with 0 rotamers are counted as 1 rotamer)

    private String rotFile;
    private String volFilename;

    private boolean canMutate = false;

    private boolean addedRotamers = false;

    public static HashMap<String,String> three2one = null;

    public boolean isAddedRotamers() {
        return addedRotamers;
    }

    public void setAddedRotamers(boolean addedRotamers) {
        this.addedRotamers = addedRotamers;
    }

    // Generic constructor
    RotamerLibrary(String rotFilename, boolean isProtein) {

        canMutate = isProtein;

        try {
            readRotLibrary(rotFilename);
        } catch (Exception e) {
            System.out.println("ERROR reading rotamer library file: "+e);
            System.exit(1);
        }

        if(three2one == null)
            initThree2One();
    }

    //Read in all of the rotamers for all amino acids from the rotFilename file
    private void readRotLibrary(String rotFilename) throws Exception {

        rotFile = rotFilename;
        volFilename = rotFile + ".vol";

        // HANDLE THE NORMAL AAs
        FileInputStream is = new FileInputStream( rotFilename );
        BufferedReader bufread = new BufferedReader(new InputStreamReader(is));
        String curLine = null;
        int curAA = 0;
        numAAallowed = 0;
        totalNumRotamers = 0;
        // Skip over comments (lines starting with !)
        curLine = bufread.readLine();
        while( curLine.charAt(0) == '!' )
            curLine = bufread.readLine();

        numAAallowed = (new Integer(getToken(curLine,1))).intValue(); //the first non-comment line is the number of AA types
        curLine = bufread.readLine();

        aaNames = new String[numAAallowed];
        numDihedrals = new int[numAAallowed];
        numRotamers = new int[numAAallowed];
        rotamerIndexOffset = new int[numAAallowed];
        dihedralAtomNames = new String[numAAallowed][][];
        rotamerValues = new int[numAAallowed][][];

        while( curLine != null ) {
            if(curLine.charAt(0) == '!') {
                curLine = bufread.readLine();
                continue;
            }

            aaNames[curAA] = getToken(curLine,1);
            numDihedrals[curAA] = (new Integer(getToken(curLine,2))).intValue();
            numRotamers[curAA] = (new Integer(getToken(curLine,3))).intValue();
            if (numRotamers[curAA]<0)
                numRotamers[curAA] = 0;

            dihedralAtomNames[curAA] = new String[numDihedrals[curAA]][4];
            rotamerValues[curAA] = new int[numRotamers[curAA]][numDihedrals[curAA]];

            // Read in the actual dihedrals
            for(int q=0; q<numDihedrals[curAA]; q++) {
                curLine = bufread.readLine();
                dihedralAtomNames[curAA][q][0] = getToken(curLine,1);
                dihedralAtomNames[curAA][q][1] = getToken(curLine,2);
                dihedralAtomNames[curAA][q][2] = getToken(curLine,3);
                dihedralAtomNames[curAA][q][3] = getToken(curLine,4);
            }
            // Read in the actual rotamers
            for(int q=0; q<numRotamers[curAA]; q++) {
                curLine = bufread.readLine();
                for(int w=0; w<numDihedrals[curAA]; w++) {
                    rotamerValues[curAA][q][w] = (new Integer(getToken(curLine,(w+1)))).intValue();
                }
            }

            rotamerIndexOffset[curAA] = totalNumRotamers;

            totalNumRotamers += numRotamers[curAA];
            if (numRotamers[curAA]<=0) //ALA or GLY
                totalNumRotamers += 1;

            curAA++;
            curLine = bufread.readLine();
        }

        bufread.close();

        if (curAA!=numAAallowed) {
            System.out.println("ERROR: not all amino acid types read from rotamer library");
            System.exit(1);
        }
    }

    //reads in the rotamer volume data from volFilename
    public void loadVolFile () {

        try {
            readRotVol(volFilename);
        } catch (Exception e) {
            System.out.println("Rotamer volumes file not found. Computing it..");
            computeAAVolumes(volFilename);
            try {
                readRotVol(volFilename);
            } catch (Exception e1) {
                System.out.println("ERROR: "+e1);
                System.exit(1);
            }
        }
    }

    //Read in all of the rotamer volumes for all amino acids from the volFilename file
    private void readRotVol(String volFilename) throws Exception {

        FileInputStream is = new FileInputStream( volFilename );
        BufferedReader bufread = new BufferedReader(new InputStreamReader(is));
        String curLine = null;

        rotamerVolumes = new float[numAAallowed][];

        is = new FileInputStream( volFilename );
        bufread = new BufferedReader(new InputStreamReader(is));

        // Skip over comments (lines starting with !)
        curLine = bufread.readLine();
        while( curLine.charAt(0) == '!' )
            curLine = bufread.readLine();

        int curResult = 0;
        while( curLine != null ) {

            int aaIndex = getAARotamerIndex(getToken(curLine,1));

            int numRotVol = numRotamers[aaIndex];
            if (numRotamers[aaIndex]==0)
                numRotVol++;

            rotamerVolumes[aaIndex] = new float[numRotVol];
            for (int j=0; j<numRotVol; j++)
                rotamerVolumes[aaIndex][j] = new Float(getToken(curLine,j+2)).floatValue();

            curLine = bufread.readLine();
            curResult++;
        }

        bufread.close();

        if (curResult!=numAAallowed) {
            System.out.println("ERROR: not all amino acid types read from rotamer volumes file");
            System.exit(1);
        }
    }

    // Uses the VolModule class to calculate the volume of each rotamer of each amino acid
    public void computeAAVolumes(String volFileName) {

        Molecule m = new Molecule();

        Amber96PolyPeptideResidue ppr = new Amber96PolyPeptideResidue();
        StrandRotamers LR = null;

        Residue res = ppr.getResidue("ala");
        //res.fullName = "ALA  ";
        m.addResidue(0,res);

        VolModule sm = new VolModule(m);
        sm.setResidueTreatment(0,1);

        LR = new StrandRotamers(rotFile,m.strand[0]);

        PrintStream printStream = setupOutputFile(volFileName);

        String aanames[] = getAAtypesAllowed();
        int numAAs = getNumAAallowed();

        for(int i=0; i<numAAs; i++) {
            if(canMutate)
                LR.changeResidueType(m,0,aanames[i],true);
            printStream.print(aanames[i] + " ");
            System.out.println(aanames[i] + " ");
            if(getNumRotamers(aanames[i])==0) {		// ALA or GLY
                float vol = sm.getMoleculeVolume(0.25f,0.0f);
                printStream.print(vol + " ");
                System.out.println(vol + " ");
            }
            for(int j=0; j<getNumRotamers(aanames[i]); j++) {
                LR.applyRotamer(m,0,j);
                float vol = sm.getMoleculeVolume(0.25f,0.0f);
                printStream.print(vol + " ");
                System.out.println(vol + " ");
            }
            printStream.println();
        }
        printStream.close();
    }

    // This function returns the rotamer index for rotamers of
    //  amino acid aaName; returns -1 if name not found
    public int getAARotamerIndex(String aaName) {

        // KER: Allow HID, HIE, and HIP to all be different now
        /*if (aaName.equalsIgnoreCase("HID") || aaName.equalsIgnoreCase("HIE") || aaName.equalsIgnoreCase("HIP")) {
        	aaName = "HIS";
        	if(debug)
        		System.out.println("ASSUMING HID/E/P is " + aaName + " for rotamer purposes.");
        }*/
        if (aaName.equalsIgnoreCase("CYX")) {
            aaName = "CYS";
            if(debug)
                System.out.println("ASSUMING CYX is " + aaName + " for rotamer purposes.");
        }

        for(int q=0; q<numAAallowed; q++) {
            if (aaNames[q].equalsIgnoreCase(aaName))
                return q;
        }
        return -1;
    }

    // Returns the name of the amino acid at array index aaIndex
    public String getAAName(int aaIndex) {
        return(aaNames[aaIndex]);
    }

    // Returns the number of rotamers for the specified amino acid type
    public int getNumRotForAAtype(int aaTypeInd) {
        return(numRotamers[aaTypeInd]);
    }

    // This function returns the number of rotamers for a given
    //  amino acid type (by name)
    public int getNumRotamers(String aaName) {
        int aaNum = getAARotamerIndex(aaName);
        return(numRotamers[aaNum]);
    }

    // Returns the number of dihedrals for the amino acid
    //  number aaTypeInd. NOTE: aaTypeInd is not an index
    //  into a molecule, it's 0..19
    public int getNumDihedrals(int aaTypeInd) {
        if (aaTypeInd != -1)
            return(numDihedrals[aaTypeInd]);
        else
            return(0);
    }


    // Returns the residue local atom numbers for dihedral dihedNum of
    //  residue resNum of strand strNum
    public int[] getDihedralInfo(Molecule m, int strNum, int resNum, int dihedNum) {

        Residue localResidue = m.strand[strNum].residue[resNum];
        String tmpName = null;
        if(localResidue.name.equalsIgnoreCase("CYX"))
            tmpName = "CYS";
        else
            tmpName = localResidue.name;

        int aaNum = getAARotamerIndex(tmpName);
        int atNum[] = new int[4];

        // Find atoms involved in the dihedral
        for(int q=0; q<localResidue.numberOfAtoms; q++) {
            for(int w=0; w<4; w++) {
                if (localResidue.atom[q].name.equalsIgnoreCase(dihedralAtomNames[aaNum][dihedNum][w])) {
                    atNum[w] = q;
                }
            }
        }
        return(atNum);
    }

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

    public int getNumAAallowed() {
        return numAAallowed;
    }

    public float [][] getRotVol() {
        return rotamerVolumes;
    }

    public String getDihedralAtomNames(int i, int j, int k) {
        return dihedralAtomNames[i][j][k];
    }

    public int getRotamerValues(int i, int j, int k) {
        return rotamerValues[i][j][k];
    }

    public int [] getRotamerIndexOffset() {
        return rotamerIndexOffset;
    }

    public int getTotalNumRotamers() {
        return totalNumRotamers;
    }

    public String [] getAAtypesAllowed() {
        return aaNames;
    }

    //Setup the file with name filename for output
    private PrintStream setupOutputFile(String fileName) {
        PrintStream logPS = null; //the output file for conf info
        try {
            FileOutputStream fileOutputStream = new FileOutputStream(fileName);
            BufferedOutputStream bufferedOutputStream = new BufferedOutputStream( fileOutputStream );
            logPS = new PrintStream( bufferedOutputStream );
        } catch (Exception ex) {
            System.out.println("ERROR: An exception occured while opening log file");
        }
        return logPS;
    }

    //KER: Adding this function so that I can dope the rotamers
    //	with the original rotamer from the input structure
    public void addRotamer(String AAname, int[] dihedVals) {
        //find AAnumber
        int aaNum = getAARotamerIndex(AAname);

        //Extend the rotamerValuesArray
        int[][] temp = new int[rotamerValues[aaNum].length+1][];
        System.arraycopy(rotamerValues[aaNum], 0, temp, 0, rotamerValues[aaNum].length);
        rotamerValues[aaNum] = temp;
        rotamerValues[aaNum][rotamerValues[aaNum].length-1] = new int[rotamerValues[aaNum][0].length];
        for(int w=0; w<dihedVals.length; w++)
            rotamerValues[aaNum][rotamerValues[aaNum].length-1][w] = dihedVals[w];

        //System.out.println(" Added for AAindex: "+aaNum+" rotamer: "+(rotamerValues[aaNum].length-1));
        totalNumRotamers++;
        numRotamers[aaNum]++;
        for(int i=aaNum+1; i<rotamerIndexOffset.length; i++)
            rotamerIndexOffset[i]++;

    }

    public static void addOrigRots(int[][] strandMut,RotamerLibrary rl, Molecule m) {
        //TODO: this assumes first strandRot is a protein
        if(! rl.isAddedRotamers()) {
            for(int str=0; str<strandMut.length; str++) {
                for(int res=0; res<strandMut[str].length; res++) {
                    Residue curRes = m.strand[str].residue[strandMut[str][res]];
                    //Get Num Dihedrals
                    int numDiheds = rl.getNumDihedrals(rl.getAARotamerIndex(curRes.name));
                    if(numDiheds<=0)
                        continue;
                    int[] diheds = new int[numDiheds];
                    int atoms[] = new int[4];
                    //get all dihedrals
                    for(int i=0; i<numDiheds; i++) {
                        atoms = rl.getDihedralInfo(m, curRes.strandNumber, curRes.strandResidueNumber, i);
                        diheds[i] = (int) Math.round(curRes.atom[atoms[3]].torsion(curRes.atom[atoms[0]], curRes.atom[atoms[1]], curRes.atom[atoms[2]]));
                    }
                    //System.out.print("Str: "+str+" Res: "+res+" ");
                    rl.addRotamer(curRes.name, diheds);
                    //System.out.println("");
                }

            }
            rl.setAddedRotamers(true);
        } else {
            System.out.println("DEBUG: ALREADY ADDED ROTAMERS");
        }
    }

    public static void initThree2One() {
        three2one = new HashMap<String,String>();
        three2one.put("ALA","A");
        three2one.put("CYS","C");
        three2one.put("ASP","D");
        three2one.put("GLU","E");
        three2one.put("PHE","F");
        three2one.put("GLY","G");
        three2one.put("HIS","H");
        three2one.put("ILE","I");
        three2one.put("LYS","K");
        three2one.put("LEU","L");
        three2one.put("MET","M");
        three2one.put("ASN","N");
        three2one.put("PRO","P");
        three2one.put("GLN","Q");
        three2one.put("ARG","R");
        three2one.put("SER","S");
        three2one.put("THR","T");
        three2one.put("VAL","V");
        three2one.put("TRP","W");
        three2one.put("TYR","Y");
    }

    public static String getOneLet(String aa3Name) {
        String res = three2one.get(aa3Name);
        if (res == null)
            res = "X";
        return res;
    }

}
