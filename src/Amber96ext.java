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
// Amber96ext.java
//
//  Version:           2.1 beta
//
//
// authors:
//    initials    name                  organization                email
//   ---------   --------------      ------------------------     ------------------------------
//     RHL        Ryan Lilien          Dartmouth College           ryan.lilien@dartmouth.edu
//     ISG		  Ivelin Georgiev	   Duke University			   ivelin.georgiev@duke.edu
//     KER        Kyle E. Roberts       Duke University         ker17@duke.edu
//     PGC        Pablo Gainza C.       Duke University         pablo.gainza@duke.edu
////////////////////////////////////////////////////////////////////////////////////////////

/*
 *
 * Major changes were made by Ryan Lilien (2000-2004) and Ivelin Georgiev (2004-2009)
 *  -Added forcefield constant reader (from parameter file). This eliminated
 *    the huge functions required to get AMBER constants
 *  -Updated the forcefield to AMBER94 with 1999 parameters (parm99.dat)
 *
 * Rewritten by Ryan Lilien based on code by Neill White
 * Many functions have been added, others removed, most have had 
 *  at least some parts rewritten. Code rewrites have often been
 *  major to fix bugs or add functionality.
 * 
 * Based on software copyrighted, 1999, by Neill White. 
 *  The author hereby grants permission to use, copy, modify, and re-distribute
 *  this software and its documentation for any purpose, provided
 *  that existing copyright notices are retained in all copies and that this
 *  notice is included verbatim in any distributions. No written agreement,
 *  license, or royalty fee is required for any of the authorized uses.
 *  Modifications to this software may be distributed provided that
 *  the nature of the modifications are clearly indicated.
 *
 */

import java.io.*;
import java.util.*;

/**
 * This class handles the energy computation for a given molecule. The Amber force field parameters
 * are read in and saved upon initialization. The EEF1 solvation parameters are also read in using
 * the EEF1 class functions. The energy (including electrostatic, vdW, and solvation) and gradient 
 * computation for the full molecule or a selected subset of the molecule atoms is performed by this class.
 */
public class Amber96ext implements ForceField, Serializable {

	// If the debug flag is set to true then additional debug statements are
	//  printed to standard out. I'm hoping that by making it a public static
	//  final variable that the compiler will be smart and not compile/include
	//  code that is unreachable.
	public static final boolean debug = false;
	
	private final int noMatchInt = 9999;
	
	boolean doSolvationE = false; //should solvation energies be computed
	
	double dielectric = 1.0;	
	boolean distDepDielect = true;
	final double constCoulomb = 332.0;

	Molecule m;
	int numberOfDihedralTerms = 0;
	int numberNonBonded = 0;
	int numHalfNonBondedTerms = 0;
	int number12Terms = 0;
	int number13Terms = 0;
	int numSolvationTerms = 0;
	int halfNBeval[], NBeval[];
	double bondStretchTerms[], angleBendTerms[], dihedralAngleTerms[];
	double nonBondedTerms[], halfNonBondedTerms[];
	double solvationTerms[];
	boolean solvExcludePairs[][];
	double D2R = 0.01745329251994329576;
	double R2D = 57.29577951308232090712;
	float vdwMultiplier = 1.0f;
	int ligandNum = -1;
	
	// The following was added by Ryan Lilien
	AminoAcidTemplates aat = null;
	GenericResidueTemplates grt = null;
	
	final int atomTypeX = -2; //the atom type number for the X wildcard atom type
	String[] atomTypeNames = null;
	float[] atomAtomicMasses = null;
	int[]   bondAtomType1 = null;
	int[]   bondAtomType2 = null;
	float[] bondHFC = null; // Harmonic Force Constant
	float[] bondEBL = null; // Equilibrium Bond Length
	int[]		angleAtomType1 = null;
	int[]		angleAtomType2 = null;
	int[]		angleAtomType3 = null;
	float[] angleHFC = null; // Harmonic Force Constant
	float[] angleEBA = null; // Equilibrium Bond Angle
	int numGeneralDihedParams = 0; //the number of generic dihedral parameters (that have wildcard X atoms)
	int[]		dihedAtomType1 = null;
	int[]		dihedAtomType2 = null;
	int[]		dihedAtomType3 = null;
	int[]		dihedAtomType4 = null;
	// Dihedral: (PK/IDIVF) * (1 + cos(PN*phi - PHASE))
	float[] dihedTerm1 = null; // (PK/IDIVF)
	int[] dihedPN = null; // Periodicity
	float[] dihedPhase = null; // Phase
	int[]		impDihedAtomType1 = null;
	int[]		impDihedAtomType2 = null;
	int[]		impDihedAtomType3 = null;
	int[]		impDihedAtomType4 = null;
	// Imprper Dihedral: PK * (1 + cos(PN*phi - PHASE))
	float[] impDihedTerm1 = null; // PK
	int[] impDihedPN = null; // Periodicity
	float[] impDihedPhase = null; // Phase
	int[]		vdwAtomType1 = null;
	float[] vdwR = null;
	float[] vdwE = null;
	int[][] equivAtoms = null;
	// For two atoms i & j
	//  B(ij) = 2 * (Ri + Rj)^6 * ei * ej
	//  A(ij) = (Ri + Rj)^12 * ei * ej
	// Those that are 1-4 separated are scaled 

	// Used to keep track of partial subsets of EV nonbonded terms
	int[] numPartHalfNonBonded = new int[0];
	int[] numPartNonBonded = new int[0];
	double[][] partHalfNonBonded = new double[0][0];
	double[][] partNonBonded = new double[0][0];
	int[][] partHalfNBeval = new int[0][0];
	int[][] partNBeval = new int[0][0];
	
	//Used to keep track of partial subsets of Dihed terms
	int numPartDihed[] = null;
	double partDihed[][] = null;
	
	//Used to keep track of partial subsets of solvation terms
	int numPartSolv[] = null;
	double partSolv[][] = null;
	
	//Solvation interactions for atoms more than 9.0A apart are already counted in dG(ref);
	//		Only count solvation interactions between atoms within 9.0A distance
	final double solvCutoff = 9.0;
	
	double solvScale = 1.0; //the scale factor for the solvation energies
	
	//The solvation parameters object
	EEF1 eef1parms = null;
	
	String amberDatInFile = "";
	

	Amber96ext(Molecule m, boolean ddDielect, double dielectConst, boolean doSolv, double solvScFactor, double vdwMult){
		
		this.m = m;
		setForcefieldInputs();
		// Read in AA and Generic templates
		try {
			aat = new AminoAcidTemplates();
			grt = new GenericResidueTemplates();
		}
		catch (FileNotFoundException e) {
			System.out.println("ERROR: Template File Not Found: "+e);
			System.exit(0);
		}
		catch ( Exception e ){
			System.out.println("ERROR: An error occurred while reading a template file: "+e);
			System.exit(0);
		}
		
		// Read in AMBER forcefield parameters
		// parm96a.dat
		try {
			readParm96();
		}
		catch (FileNotFoundException e) {
			System.out.println("ERROR: File Not Found: "+e);
			System.exit(0);
		}
		catch (IOException e) {
			System.out.println("ERROR: An error occurred while reading file: "+e);
			System.exit(0);
		}
		catch ( Exception e ){
			System.out.println("ERROR: An error occurred while reading file: "+e);
			System.exit(0);
		}
		
		// Read in the EEF1 solvation parameters
		try {
			eef1parms = new EEF1(m);
			eef1parms.readEEF1parm();
		}
		catch ( Exception e ){
			System.out.println("ERROR: An error occurred while reading file: "+e);
			System.exit(0);
		}
		
		distDepDielect = ddDielect;
		dielectric = dielectConst;
		
		vdwMultiplier = (float)vdwMult;
		
		doSolvationE = doSolv;
		solvScale = solvScFactor;
	}

	//************************************
	// This function reads the AMBER forcefield parameter file
	//  parm96a.dat
	// Although this function is meant to be relatively generic
	//  it is optimized for reading parm96a.dat. Slight changes
	//  will most likely be required to read other parameter
	//  files. Reading of other files should be done in other
	//  functions
	private void readParm96() throws Exception {
	
		FileInputStream is = new FileInputStream( EnvironmentVars.getDataDir().concat(amberDatInFile) );
		BufferedReader bufread = new BufferedReader(new InputStreamReader(is));
		String curLine = null, tmpStr = null;
		int tmpInt = 0;
		
		final int initSize = 10; //the initial size of the arrays to store the data that is read

		// Skip over the first line of header info
		curLine = bufread.readLine();
		
		// 1. Read atom names and atomic masses
		atomTypeNames = new String[initSize];
		atomAtomicMasses = new float[initSize];
		curLine = bufread.readLine();		
		tmpInt = 0; // temporary integer
		// Until we're at a blank line (or until we've read numAtomTypes)
		while (!(getToken(curLine,1).equals(""))) {
			
			if (tmpInt>=atomTypeNames.length){ //double the array sizes
				atomTypeNames = doubleArraySize(atomTypeNames);
				atomAtomicMasses = doubleArraySize(atomAtomicMasses);
			}
			
			atomTypeNames[tmpInt] = getToken(curLine,1);  // snag atom name
			atomAtomicMasses[tmpInt] = (new Float(getToken(curLine,2))).floatValue();
			tmpInt++;
			curLine = bufread.readLine();
		}
		atomTypeNames = reduceArraySize(atomTypeNames,tmpInt);
		atomAtomicMasses = reduceArraySize(atomAtomicMasses,tmpInt);
		

		// Skip unknown line
		curLine = bufread.readLine();

		// 2. Read Bonds
		bondAtomType1 = new int[initSize];
		bondAtomType2 = new int[initSize];
		bondHFC = new float[initSize];
		bondEBL = new float[initSize];
		curLine = bufread.readLine();		
		tmpInt = 0;
		while (!(getToken(curLine,1).equals(""))) {
			
			if (tmpInt>=bondAtomType1.length){
				bondAtomType1 = doubleArraySize(bondAtomType1);
				bondAtomType2 = doubleArraySize(bondAtomType2);
				bondHFC = doubleArraySize(bondHFC);
				bondEBL = doubleArraySize(bondEBL);
			}
			
			//tmpStr = curLine.substring(0,5);
			bondAtomType1[tmpInt] = atomTypeToInt(getDashedToken(curLine, 1));
			bondAtomType2[tmpInt] = atomTypeToInt(getDashedToken(curLine, 2));
			bondHFC[tmpInt] = (new Float(getDashedToken(curLine,3))).floatValue();
			bondEBL[tmpInt] = (new Float(getDashedToken(curLine,4))).floatValue();
			tmpInt++;
			curLine = bufread.readLine();
		}
		bondAtomType1 = reduceArraySize(bondAtomType1,tmpInt);
		bondAtomType2 = reduceArraySize(bondAtomType2,tmpInt);
		bondHFC = reduceArraySize(bondHFC,tmpInt);
		bondEBL = reduceArraySize(bondEBL,tmpInt);
		

		// 3. Read Angles
		angleAtomType1 = new int[initSize];
		angleAtomType2 = new int[initSize];
		angleAtomType3 = new int[initSize];
		angleHFC = new float[initSize];
		angleEBA = new float[initSize];
		curLine = bufread.readLine();		
		tmpInt = 0;
		while (!(getToken(curLine,1).equals(""))) {
			
			if (tmpInt>=angleAtomType1.length){
				angleAtomType1 = doubleArraySize(angleAtomType1);
				angleAtomType2 = doubleArraySize(angleAtomType2);
				angleAtomType3 = doubleArraySize(angleAtomType3);
				angleHFC = doubleArraySize(angleHFC);
				angleEBA = doubleArraySize(angleEBA);
			}
			
			//tmpStr = curLine.substring(0,8);
			angleAtomType1[tmpInt] = atomTypeToInt(getDashedToken(curLine,1));
			angleAtomType2[tmpInt] = atomTypeToInt(getDashedToken(curLine,2));
			angleAtomType3[tmpInt] = atomTypeToInt(getDashedToken(curLine,3));
			angleHFC[tmpInt] = (new Float(getDashedToken(curLine,4))).floatValue();
			angleEBA[tmpInt] = (new Float(getDashedToken(curLine,5))).floatValue();
			tmpInt++;
			curLine = bufread.readLine();
		}
		angleAtomType1 = reduceArraySize(angleAtomType1,tmpInt);
		angleAtomType2 = reduceArraySize(angleAtomType2,tmpInt);
		angleAtomType3 = reduceArraySize(angleAtomType3,tmpInt);
		angleHFC = reduceArraySize(angleHFC,tmpInt);
		angleEBA = reduceArraySize(angleEBA,tmpInt);
		
		
		// 4. Read Dihedrals
		numGeneralDihedParams = 0;
		dihedAtomType1 = new int[initSize];
		dihedAtomType2 = new int[initSize];
		dihedAtomType3 = new int[initSize];
		dihedAtomType4 = new int[initSize];
		dihedTerm1 = new float[initSize];
		dihedPhase = new float[initSize];
		dihedPN = new int[initSize];
		curLine = bufread.readLine();		
		tmpInt = 0;
		float tmpFlt = 0.0f;
		while (!(getToken(curLine,1).equals(""))) {
			
			if (tmpInt>=dihedAtomType1.length){
				dihedAtomType1 = doubleArraySize(dihedAtomType1);
				dihedAtomType2 = doubleArraySize(dihedAtomType2);
				dihedAtomType3 = doubleArraySize(dihedAtomType3);
				dihedAtomType4 = doubleArraySize(dihedAtomType4);
				dihedTerm1 = doubleArraySize(dihedTerm1);
				dihedPhase = doubleArraySize(dihedPhase);
				dihedPN = doubleArraySize(dihedPN);
			}
			
			//tmpStr = curLine.substring(0,11);
			dihedAtomType1[tmpInt] = atomTypeToInt(getDashedToken(curLine,1));
			dihedAtomType2[tmpInt] = atomTypeToInt(getDashedToken(curLine,2));
			dihedAtomType3[tmpInt] = atomTypeToInt(getDashedToken(curLine,3));
			dihedAtomType4[tmpInt] = atomTypeToInt(getDashedToken(curLine,4));
			
			if ( dihedAtomType1[tmpInt]==atomTypeX || dihedAtomType2[tmpInt]==atomTypeX || dihedAtomType3[tmpInt]==atomTypeX || dihedAtomType4[tmpInt]==atomTypeX ) //at least one of the atoms is a wildcard
				numGeneralDihedParams++;
			
			tmpFlt = (new Float(getDashedToken(curLine,5))).floatValue();
			dihedTerm1[tmpInt] = (new Float(getDashedToken(curLine,6))).floatValue() / tmpFlt;
			dihedPhase[tmpInt] = (new Float(getDashedToken(curLine,7))).floatValue();
			dihedPN[tmpInt] = (new Float(getDashedToken(curLine,8))).intValue();
			// If dihedPN is negative then there are one or more additional terms
			//  nothing fancy needs to be done because they will all be read in anyway
			//  but we do need to correct the sign
			if (dihedPN[tmpInt] < 0)
				dihedPN[tmpInt] = -dihedPN[tmpInt];
			tmpInt++;
			curLine = bufread.readLine();
		}
		dihedAtomType1 = reduceArraySize(dihedAtomType1,tmpInt);
		dihedAtomType2 = reduceArraySize(dihedAtomType2,tmpInt);
		dihedAtomType3 = reduceArraySize(dihedAtomType3,tmpInt);
		dihedAtomType4 = reduceArraySize(dihedAtomType4,tmpInt);
		dihedTerm1 = reduceArraySize(dihedTerm1,tmpInt);
		dihedPhase = reduceArraySize(dihedPhase,tmpInt);
		dihedPN = reduceArraySize(dihedPN,tmpInt);
		

		// 5. Read Improper Dihedrals
		impDihedAtomType1 = new int[initSize];
		impDihedAtomType2 = new int[initSize];
		impDihedAtomType3 = new int[initSize];
		impDihedAtomType4 = new int[initSize];
		impDihedTerm1 = new float[initSize];
		impDihedPhase = new float[initSize];
		impDihedPN = new int[initSize];
		curLine = bufread.readLine();		
		tmpInt = 0;
		while (!(getToken(curLine,1).equals(""))) {
			
			if (tmpInt>=impDihedAtomType1.length){
				impDihedAtomType1 = doubleArraySize(impDihedAtomType1);
				impDihedAtomType2 = doubleArraySize(impDihedAtomType2);
				impDihedAtomType3 = doubleArraySize(impDihedAtomType3);
				impDihedAtomType4 = doubleArraySize(impDihedAtomType4);
				impDihedTerm1 = doubleArraySize(impDihedTerm1);
				impDihedPhase = doubleArraySize(impDihedPhase);
				impDihedPN = doubleArraySize(impDihedPN);
			}
			
			//tmpStr = curLine.substring(0,11);
			impDihedAtomType1[tmpInt] = atomTypeToInt(getDashedToken(curLine,1));
			impDihedAtomType2[tmpInt] = atomTypeToInt(getDashedToken(curLine,2));
			impDihedAtomType3[tmpInt] = atomTypeToInt(getDashedToken(curLine,3));
			impDihedAtomType4[tmpInt] = atomTypeToInt(getDashedToken(curLine,4));
			impDihedTerm1[tmpInt] = (new Float(getDashedToken(curLine,5))).floatValue();
			impDihedPhase[tmpInt] = (new Float(getDashedToken(curLine,6))).floatValue();
			impDihedPN[tmpInt] = (new Float(getDashedToken(curLine,7))).intValue();
			tmpInt++;
			curLine = bufread.readLine();
		}
		impDihedAtomType1 = reduceArraySize(impDihedAtomType1,tmpInt);
		impDihedAtomType2 = reduceArraySize(impDihedAtomType2,tmpInt);
		impDihedAtomType3 = reduceArraySize(impDihedAtomType3,tmpInt);
		impDihedAtomType4 = reduceArraySize(impDihedAtomType4,tmpInt);
		impDihedTerm1 = reduceArraySize(impDihedTerm1,tmpInt);
		impDihedPhase = reduceArraySize(impDihedPhase,tmpInt);
		impDihedPN = reduceArraySize(impDihedPN,tmpInt);
		

		// Skip 2 lines (we might also be able to go until the keyword MOD4
		curLine = bufread.readLine();
		curLine = bufread.readLine();

		// Read the equivalence lines
		// The first atomnum in equivAtoms is the main atom and numbers with index 1..n are
		//  equivalent to the atom in index 0.
		equivAtoms = new int[initSize][];
		curLine = bufread.readLine();
		tmpInt = 0;
		while (!(getToken(curLine,1).equals(""))) {
			
			if (tmpInt>=equivAtoms.length){
				equivAtoms = doubleArraySize(equivAtoms);
			}
			
			int numEquivAtoms = (new StringTokenizer(curLine," ,;\t\n\r\f")).countTokens();
			equivAtoms[tmpInt] = new int[numEquivAtoms];
			
			for(int q=0;q<equivAtoms[tmpInt].length;q++)
				equivAtoms[tmpInt][q] = -noMatchInt;
			int tmpInt2=1;
			while (!getToken(curLine,tmpInt2).equalsIgnoreCase("")) { 
				equivAtoms[tmpInt][tmpInt2-1] = atomTypeToInt( getToken(curLine,tmpInt2) );
				tmpInt2++;
			}
			tmpInt++;
			curLine = bufread.readLine();
		}
		equivAtoms = reduceArraySize(equivAtoms,tmpInt);
		
				
		// Skip a line (we might also be able to go until the keyword MOD4
		curLine = bufread.readLine();		

		// 6. Read vdw
		vdwAtomType1 = new int[initSize];
		vdwR = new float[initSize];
		vdwE = new float[initSize];
		curLine = bufread.readLine();
		tmpInt = 0;
		while (!(getToken(curLine,1).equals(""))) {
			
			if (tmpInt>=vdwAtomType1.length){
				vdwAtomType1 = doubleArraySize(vdwAtomType1);
				vdwR = doubleArraySize(vdwR);
				vdwE = doubleArraySize(vdwE);
			}
			
			vdwAtomType1[tmpInt] = atomTypeToInt(getToken(curLine,1));
			vdwR[tmpInt] = (new Float(getToken(curLine,2))).floatValue();
			vdwE[tmpInt] = (new Float(getToken(curLine,3))).floatValue();
			tmpInt++;
			curLine = bufread.readLine();
		}
		vdwAtomType1 = reduceArraySize(vdwAtomType1,tmpInt);
		vdwR = reduceArraySize(vdwR,tmpInt);
		vdwE = reduceArraySize(vdwE,tmpInt);		
		
		bufread.close();

	// DEBUG START * Good to keep for when the parameter file changes to
	//  make sure you're reading it correctly
/*	System.out.println("ATOM TYPES");
	for(int q=0;q<5;q++) {
		System.out.println(q + " " + atomTypeNames[q] + " mass:" + atomAtomicMasses[q]);
	}
	System.out.println("BONDS");
	for(int q=0;q<5;q++) {
		System.out.println(q + " " + bondAtomType1[q] + " " + bondAtomType2[q] + " HFC:" + 
			bondHFC[q] + " EBL:" + bondEBL[q]);
	}
	System.out.println("ANGLES");
	for(int q=0;q<5;q++) {
		System.out.println(q + " " + angleAtomType1[q] + " " + angleAtomType2[q] + " " + angleAtomType3[q]
		 + " HFC:" + angleHFC[q] + " EBA:" + angleEBA[q]);
	}
	System.out.println("DIHED");
	for(int q=0;q<5;q++) {
		System.out.println(q + " " + dihedAtomType1[q] + " " + dihedAtomType2[q] + " " + dihedAtomType3[q]
		 + " " + dihedAtomType4[q] + " term1:" + dihedTerm1[q] + " PN:" + dihedPN[q] + " Phase:" + dihedPhase[q]);
	}
	System.out.println("IMPROP DIHED");
	for(int q=0;q<5;q++) {
		System.out.println(q + " " + impDihedAtomType1[q] + " " + impDihedAtomType2[q] + " " + impDihedAtomType3[q]
		 + " " + impDihedAtomType4[q] + " term1:" + impDihedTerm1[q] + " PN:" + impDihedPN[q] + " Phase:" + impDihedPhase[q]);
	}
	System.out.println("VDW");
	for(int q=0;q<5;q++) {
		System.out.println(q + " " + vdwAtomType1[q] + " " + " R:" + vdwR[q] + " E:" + vdwE[q]);
	}
	// DEBUG END
*/
	
	}


	private String getDashedToken(String s, int x) {
	
		int curNum = 1;	
		StringTokenizer st = new StringTokenizer(s," ,;\t\n\r\f-");
		
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
	

	// This function returns the numeric atom type based on the string atom type
	// If atom type is 'x' then return atomTypeX which means it's a wildcard
	private int atomTypeToInt(String s) {
		s = s.trim();
		if (s.equalsIgnoreCase("x"))
			return atomTypeX;
		for(int q=0;q<atomTypeNames.length;q++) {
			if (atomTypeNames[q].equalsIgnoreCase(s))
				return q;
		}
		return -1;
	}
	
	
	// This function searches the bond constants for the atoms specified
	//  and returns the approprate constants 
	public boolean getStretchParameters(int atomType1, int atomType2,
		double forceConstant[],	double equilibriumDistance[]){

		int numGeneric = 5, tmpInt = 0;
		double tmpFC = 0.0f, tmpED = 0.0f;
		boolean matched = false;
		for(int q=0;q<bondAtomType1.length;q++) {
			if (match2Atoms(atomType1, atomType2, bondAtomType1[q], bondAtomType2[q])) {
				tmpInt = 0;
				matched = true;
				if (bondAtomType1[q] == atomTypeX)
					tmpInt++;
				if (bondAtomType2[q] == atomTypeX)
					tmpInt++;
				if (tmpInt < numGeneric) {
					numGeneric = tmpInt;
					tmpFC = bondHFC[q];
					tmpED = bondEBL[q];
				}
			}
		}
		
		if (matched) {
			forceConstant[0] = tmpFC;
			equilibriumDistance[0] = tmpED;
			return(true);
		}
		else {
			forceConstant[0] = 317;
      equilibriumDistance[0] = 1.522;
			System.out.println("Ambstretch DEFAULTING TO C-CT");
			return(false);
		}
	}
	
	// This function searches the angle constants for the atoms specified
	//  and returns the approprate constants
	public boolean getBendParameters(int atomType1, int atomType2, int atomType3, 
		double forceConstant[], double equilibriumAngle[]){

		if (atomType3 < atomType1){
			int temp = atomType3;
			atomType3 = atomType1;
			atomType1 = temp;
		}

		int numGeneric = 5, tmpInt = 0;
		double tmpFC = 0.0f, tmpEA = 0.0f;
		boolean matched = false;
		for(int q=0;q<angleAtomType1.length;q++) {
			if (match3Atoms(atomType1, atomType2, atomType3, angleAtomType1[q], 
				angleAtomType2[q], angleAtomType3[q])) {
				tmpInt = 0;
				matched = true;
				if (angleAtomType1[q] == atomTypeX)
					tmpInt++;
				if (angleAtomType2[q] == atomTypeX)
					tmpInt++;
				if (angleAtomType3[q] == atomTypeX)
					tmpInt++;
				if (tmpInt < numGeneric) {
					numGeneric = tmpInt;
					tmpFC = angleHFC[q];
					tmpEA = angleEBA[q];
				}
			}
		}
		
		if (matched) {
			forceConstant[0] = tmpFC;
			equilibriumAngle[0] = tmpEA;
			return(true);
		}
		else {
			forceConstant[0] = 63.000000;
			equilibriumAngle[0] = 113.100000;
			System.out.println( "AmberBend: Could not find correct angle, defaulting to CC-CT-CT");
			return( false );
		}
	}

	// This function searches the dihedral constants for the atoms specified
	//  and returns the approprate constants
	public boolean getTorsionParameters(int atomType1, int atomType2, int atomType3, 
		int atomType4, double forceConstant[], double equilibriumAngle[], int terms[], 
		int multiplicity[]){

		if (atomType3 < atomType2){
			int temp = atomType3;
			atomType3 = atomType2;
			atomType2 = temp;
			temp = atomType4;
			atomType4 = atomType1;
			atomType1 = temp;
		}

		// First search through generic torsions (those terms containing wildcards)
		boolean matched = false;
		for(int q=0;q<numGeneralDihedParams;q++) {
			if (match2Atoms(atomType2, atomType3, dihedAtomType2[q], dihedAtomType3[q])) {
				matched = true;
				forceConstant[0] = dihedTerm1[q];
				equilibriumAngle[0] = dihedPhase[q];
				terms[0] = dihedPN[q];
				multiplicity[0] = 0;
			}
		}
		
		// According to the original paper "any specific parameter, such as OS-CH-CH-OS,
		//  overrides any general parameter."
		int forceCounter = 0, eqCounter = 0, multCounter = 0, termCounter = 0;
		for(int q=numGeneralDihedParams;q<dihedAtomType1.length;q++) {
			if (match4Atoms(atomType1, atomType2, atomType3, atomType4, dihedAtomType1[q],
						dihedAtomType2[q], dihedAtomType3[q], dihedAtomType4[q])) {
				matched = true;
				forceConstant[forceCounter++] = dihedTerm1[q];
				equilibriumAngle[eqCounter++] = dihedPhase[q];
				terms[termCounter++] = dihedPN[q];
				multiplicity[0] = multCounter++;
			}
		}
		
		if (matched) {
			return(true);
		}
		else {
			System.out.println("AmberDihed: Could not find correct torsion");
			return( false );
		}
	}
	
	// This function attempts to match two atoms for a bond. An atom type
	//  of atomTypeX is a generic term and can match anything
	private boolean match2Atoms(int atType1, int atType2, int known1, int known2) {
	
		if ((atType1 == known1) && (atType2 == known2))
			return(true);
		if ((atType1 == known2) && (atType2 == known1))
			return(true);
		if ((atType1 == known1) && (known2 == atomTypeX))
			return(true);
		if ((atType1 == known2) && (known1 == atomTypeX))
			return(true);
		if ((atType2 == known1) && (known2 == atomTypeX))
			return(true);
		if ((atType2 == known2) && (known1 == atomTypeX))
			return(true);
		return(false);
	}
	
	// This function attempts to match three atoms for an angle. An atom
	//  type of atomTypeX is a generic term and can match anything
	private boolean match3Atoms(int atType1, int atType2, int atType3,
		int known1, int known2, int known3) {
		
		if ((atType1 == known1) && (atType2 == known2) && (atType3 == known3))
			return(true);
		if ((atType1 == known3) && (atType2 == known2) && (atType3 == known1))
			return(true);
		if ((known1 == atomTypeX) && (atType2 == known2) && (atType3 == known3))
			return(true);
		if ((known1 == atomTypeX) && (atType2 == known2) && (atType1 == known3))
			return(true);
		if ((atType1 == known1) && (atType2 == known2) && (known3 == atomTypeX))
			return(true);
		if ((atType3 == known1) && (atType2 == known2) && (known3 == atomTypeX))
			return(true);
		if ((atType1 == known1) && (known2 == atomTypeX) && (atType3 == known3))
			return(true);
		if ((atType3 == known1) && (known2 == atomTypeX) && (atType1 == known3))
			return(true);
		return(false);
	}

	// This function attempts to match four atoms for a dihedral (no generic atoms)
	private boolean match4Atoms(int atType1, int atType2, int atType3, int atType4,
		int known1, int known2, int known3, int known4) {
	
		if ((atType1 == known1) && (atType2 == known2) &&
			(atType3 == known3) && (atType4 == known4))
			return(true);
		else if ((atType1 == known4) && (atType2 == known3) &&
			(atType3 == known2) && (atType4 == known1))
					return(true);

		return(false);
	}

	
	// This function returns the equivalent class for the given atomtype
	private int getEquivalentType(int atomType) {

		for(int i=0;i<equivAtoms.length;i++) {
			for(int j=1;j<equivAtoms[i].length;j++) {
				if (atomType == equivAtoms[i][j])
					return(equivAtoms[i][0]);
			}
		}
		return(-1);
	}


	// This function returns the r and epsilon paramters for a given atom type
	public boolean getNonBondedParameters(int atomType, double r[], double epsilon[]){

		for(int q=0;q<vdwAtomType1.length;q++) {
			if (vdwAtomType1[q]==atomType) {
				r[0]=vdwR[q];
				epsilon[0]=vdwE[q];
				return (true);
			}
		}

		// Check for equivalent atoms
		int equivType = getEquivalentType(atomType);
		for(int q=0;q<vdwAtomType1.length;q++) {
			if (vdwAtomType1[q]==equivType) {
				r[0]=vdwR[q];
				epsilon[0]=vdwE[q];
				return (true);
			}
		}
		
		return(false);
	}

	
	// Checks residues res against template templateRes for atom type assignment
	// If residues don't match, noMatchInt is returned
	// If residues do match, the difference between the number of atoms in the template
	//   and res is returned
	private int checkAAType(Molecule m, Residue templateRes, Residue res, int[] atArray) { 
		boolean[] atUsedArray = new boolean[35];
		for(int q=0; q<35; q++) {
			atArray[q] = -1;
			atUsedArray[q] = false;
		}			

		if(!m.connectivity12Valid)
			m.establishConnectivity(true);

		// Find backbone N
		for(int q=0; q<res.numberOfAtoms; q++) {
			Atom at = res.atom[q];
			if (at.name.equalsIgnoreCase("N")) {
				atUsedArray[q]=true;
				atArray[0]=q;
			}
		} // end for
					
		// Now chain off the N (starting at atom 1)
		// Require all hydrogens to match (last parameter of recurseAA() is 'true')
		if (recurseAA(m,templateRes,res,atArray,atUsedArray,1,true)) {
			int d = (templateRes.numberOfAtoms - res.numberOfAtoms);
			if (d < 0)
				return (-d);
			return (d);
		}

		return noMatchInt;
	}
	
	
	// Assigns AMBER types as well as default AMBER charge for atoms
	//  in molecule m, residue nr (the nr index is molecule based)
	public boolean assignAA(Molecule m, int nr) {
	
		Residue res = m.residue[nr];
		Residue templateRes = null, templateResNT = null, templateResCT = null;
		int[] atArray = new int[35];
		int[] atArrayNT = new int[35];
		int[] atArrayCT = new int[35];
		int d = noMatchInt, dNT = noMatchInt, dCT = noMatchInt;

		// Attempt to find the resiude name in the residue assigner
		for (int i=0;i<aat.numAAs;i++) {
			if (aat.aaResidues[i].name.equalsIgnoreCase(res.name)) {
				if (d!=noMatchInt)
					System.out.println("** Residue already matched to nonterm template, this should not happen **");
				templateRes = aat.aaResidues[i];
				d = checkAAType(m,templateRes,res,atArray);
			}
		}
		
		for (int i=0;i<aat.numAANTs;i++) {
			if (aat.aaNTResidues[i].name.equalsIgnoreCase(res.name)) {
				if (dNT!=noMatchInt)
					System.out.println("** Residue already matched to NT template, this should not happen **");
				templateResNT = aat.aaNTResidues[i];
				dNT = checkAAType(m,templateResNT,res,atArrayNT);
			}
		}
		
		for (int i=0;i<aat.numAACTs;i++) {
			if (aat.aaCTResidues[i].name.equalsIgnoreCase(res.name)) {
				if (dCT!=noMatchInt)
					System.out.println("** Residue already matched to CT template, this should not happen **");
				templateResCT = aat.aaCTResidues[i];
				dCT = checkAAType(m,templateResCT,res,atArrayCT);
			}
		}

		if ((d == noMatchInt) && (dNT == noMatchInt) && (dCT == noMatchInt))
			return false;
		
		// assignment worked so determine which template was the best match then
		//  copy over charges and atom types
		Residue tR = null;
		int[] atA = null;
		if ((d < dNT) && (d < dCT)) {
			tR = templateRes;
			tR.cterm = false;
			tR.nterm = false;
			atA = atArray;
			if (debug)
				System.out.println("Residue " + nr + " matched to non-terminal: " + res.name);
		}
		else if (dNT < dCT) {
			tR = templateResNT;
			atA = atArrayNT;
			tR.cterm = false;
			tR.nterm = true;
			if (debug)
				System.out.println("Residue " + nr + " matched to N-terminal: " + res.name);
		}
		else {
			tR = templateResCT;
			tR.cterm = true;
			tR.nterm = false;
			atA = atArrayCT;
			if (debug)
				System.out.println("Residue " + nr + " matched to C-terminal: " + res.name);
		}
		for(int q=0;q<tR.numberOfAtoms;q++) {
			if (atA[q] != -1) {
				res.atom[atA[q]].forceFieldType = tR.atom[q].forceFieldType;
//			System.out.println(" FFT: " + atA[q] + " " + tR.atom[q].forceFieldType + " charge: " + tR.atom[q].charge);
				res.atom[atA[q]].charge = tR.atom[q].charge;
				assignNumericalType(res.atom[atA[q]], res.atom[atA[q]].forceFieldType);
			}
		}
		res.cterm = tR.cterm;
		res.nterm = tR.nterm;
		res.ffAssigned = true;		
		return true;
	}

	/******************/
	private boolean recurseAA(Molecule m, Residue templateRes, Residue res,
		int[] atArray, boolean[] atUsedArray, int curAtResNum, boolean requireHydrogens) {
		
		if (curAtResNum >= templateRes.numberOfAtoms)
			return true;

		// the trick with the next line is that in the template residue
		//  each atom has one bond which points to a higher atom in the
		//  residue it's bonded to. this is used to get a handle to the
		//  this atom's parent node in the tree of atoms that is the
		//  template residue
		int rootNum = templateRes.atom[curAtResNum].bond[0];
		String templEleType = templateRes.atom[curAtResNum].elementType;
		for(int w=0; w<res.numberOfAtoms; w++) {
			if (atUsedArray[w]==false) {
				if (res.atom[w].elementType.equalsIgnoreCase(templEleType)) {
					int curAtMolNum = res.atom[w].moleculeAtomNumber;
					for(int k=0; k<m.connected[curAtMolNum][0]; k++) {
						Atom tmpAt = m.atom[m.connected[curAtMolNum][k+1]];
						if (tmpAt.residueAtomNumber == atArray[rootNum]) {
							atUsedArray[w]=true;
							atArray[curAtResNum]=w;
							if (recurseAA(m,templateRes,res,atArray,atUsedArray,curAtResNum+1,requireHydrogens)) {
								if (debug)
									System.out.println("FOUND MATCH, AT:"+(w+1)+" to template "+curAtResNum);
								return true;
							}
							atUsedArray[w]=false;
							atArray[curAtResNum]=-1;
						}
					}
				}
			} // end if atUsedArray[w]
		}

	// if we got here then no assignment worked, but we allow hydrogens to slip
	//  so if this is a hydrogen atom then we let it slide
	if (templEleType.equalsIgnoreCase("H") && !requireHydrogens)
		if (recurseAA(m,templateRes,res,atArray,atUsedArray,curAtResNum+1,requireHydrogens)) {
			return true;
		}
	
	return false;
	}

	// Calculates AMBER atom types using molecule templates
	//  for one amino acid
	public boolean calculateOneAAWithTemplates(Residue res){
		// Blank the atom types
		for(int z=0; z<res.numberOfAtoms; z++) {
			res.atom[z].forceFieldType = "";
			res.atom[z].type = -1;
		}
		return(assignAA(m,res.moleculeResidueNumber));
	}


	// Checks residues res against template templateRes for atom type assignment
	// If residues don't match, noMatchInt is returned
	// If residues do match, the difference between the number of atoms in the template
	//   and res is returned
	private int checkGRType(Molecule m, Residue templateRes, Residue res, int[] atArray) { 
		boolean[] atUsedArray = new boolean[atArray.length];
		for(int q=0; q<atUsedArray.length; q++) {
			atArray[q] = -1;
			atUsedArray[q] = false;
		}			

		if(!m.connectivity12Valid)
			m.establishConnectivity(true);
		
		
		// Find the anchor atom for each generic residue in the input PDB file: test all atoms;
		// The anchor atom should be the first atom of the template residue;
		// Require all hydrogens to match (last parameter of recurseAA() is 'true')
		for(int q=0; q<res.numberOfAtoms; q++) {
			atUsedArray[q]=true;
			atArray[0]=q;
			
			// Now chain off the current anchor atom (starting at atom 1)
			// Note that we can use recurseAA because that function is NOT AA specific
			if (recurseAA(m,templateRes,res,atArray,atUsedArray,1,true)) {
				int d = (templateRes.numberOfAtoms - res.numberOfAtoms);
				if (d < 0)
					return (-d);
				return (d);
			}
			
			atUsedArray[q]=false;
			atArray[0]=-1;
		}

		return noMatchInt;
	}
	
	// Assigns AMBER types as well as default AMBER charge for atoms
	//  in molecule m, residue nr (the nr index is molecule based)
	public boolean assignGR(Molecule m, int nr) {
	
		Residue res = m.residue[nr];
		Residue templateRes = null;
		int[] atArray = new int[m.residue[nr].numberOfAtoms*2];
		int d = noMatchInt;

		// Attempt to find the resiude name in the generic residue assigner
		for (int i=0;i<grt.numGRs;i++) {
			if (grt.grResidues[i].name.equalsIgnoreCase(res.name)) {
				if (d!=noMatchInt)
					System.out.println("** Residue already matched to template, this should not happen **");
				templateRes = grt.grResidues[i];
				d = checkGRType(m,templateRes,res,atArray);
			}
		}
		
		if (d == noMatchInt)
			return false;
		
		// assignment worked, copy over charges and atom types
		for(int q=0;q<templateRes.numberOfAtoms;q++) {
			if (atArray[q] != -1) {
				res.atom[atArray[q]].forceFieldType = templateRes.atom[q].forceFieldType;
//			System.out.println(" FFT: " + atArray[q] + " " + templateRes.atom[q].forceFieldType + " charge: " + templateRes.atom[q].charge);
				res.atom[atArray[q]].charge = templateRes.atom[q].charge;
				assignNumericalType(res.atom[atArray[q]], res.atom[atArray[q]].forceFieldType);
			}
		}
		res.ffAssigned = true;		
		return true;
	}

	// Calculates AMBER atom types using molecule templates
	//  for one residue using the generic residues class
	public boolean calculateOneGRWithTemplates(Residue res){
		// Blank the atom types
		for(int z=0; z<res.numberOfAtoms; z++) {
			res.atom[z].forceFieldType = "";
			res.atom[z].type = -1;
		}
		return(assignGR(m,res.moleculeResidueNumber));
	}


	// Calculates AMBER atom types using molecule templates
	// Currently the templates are amino acid based
	// To use the amino acid templates, each strand should
	//  have its isProtein boolean set to true
	public void calculateTypesWithTemplates() {
		//KER: Ensure molecule is connected
		if(!m.connectivityValid)
			m.establishConnectivity(false);
	
		// Modified to use AA templates if possible
		for(int q=0; q<m.numberOfStrands; q++) {
			if (m.strand[q].isProtein) {
				for(int w=0; w<m.strand[q].numberOfResidues; w++) {
					//if ( m.strand[q].residue[w].getEnergyEvalSC() || m.strand[q].residue[w].getEnergyEvalBB() ){
						Residue res = m.strand[q].residue[w];
						if(res.cofactor){ //KER: I allow cofactor residues in a protein strand cause they will never move/rotate
							if(!calculateOneGRWithTemplates(res)){
								System.out.println("WARNING: UNABLE TO FIND GENERIC TEMPLATE FOR: " + m.strand[q].residue[w].fullName);
							}
							else { // we were able to assign types with template
								// System.out.println("Assigned strand: " + q + " res: " + w + " to an AAT");
								// Make sure all atoms were assigned, if not assign them by standard methods
								for(int z=0; z<res.numberOfAtoms; z++) {
									if (res.atom[z].forceFieldType.equals("")) {
										System.out.println(" patching res: " + w + " atom: " + z + " " + res.atom[z].elementType + " " + res.atom[z].name);
										System.out.println("WARNING: Unable to patch residue\n");
									}
								}
							}
						}
						else{
							if(!calculateOneAAWithTemplates(res)){
								System.out.println("WARNING: UNABLE TO FIND AA TEMPLATE FOR: " + m.strand[q].residue[w].fullName);
							}
							else { // we were able to assign types with template
								// System.out.println("Assigned strand: " + q + " res: " + w + " to an AAT");
								// Make sure all atoms were assigned, if not assign them by standard methods
								for(int z=0; z<res.numberOfAtoms; z++) {
									if (res.atom[z].forceFieldType.equals("")) {
										System.out.println(" patching res: " + w + " atom: " + z + " " + res.atom[z].elementType + " " + res.atom[z].name);
										System.out.println("WARNING: Unable to patch residue\n");
									}
								}
							}
						}
					//}
				}
			}
			else { // strand is not a protein
				for(int w=0; w<m.strand[q].numberOfResidues; w++) {
					//if ( m.strand[q].residue[w].getEnergyEvalSC() || m.strand[q].residue[w].getEnergyEvalBB() ){
						Residue res = m.strand[q].residue[w];
						if(!calculateOneGRWithTemplates(res)){
							System.out.println("WARNING: UNABLE TO FIND GENERIC TEMPLATE FOR: " + m.strand[q].residue[w].fullName);
						}
						else { // we were able to assign types with template
							// System.out.println("Assigned strand: " + q + " res: " + w + " to an AAT");
							// Make sure all atoms were assigned, if not assign them by standard methods
							for(int z=0; z<res.numberOfAtoms; z++) {
								if (res.atom[z].forceFieldType.equals("")) {
									System.out.println(" patching res: " + w + " atom: " + z + " " + res.atom[z].elementType + " " + res.atom[z].name);
									System.out.println("WARNING: Unable to patch residue\n");
								}
							}
						}
					//}
				}
			} 
		}

	}

	// Calculates AMBER type for atoms in the molecule already
	//   assigned to this object
	public void calculateTypes(){

		calculateTypesWithTemplates();
	}


	// This function goes from forcefieldtype (string) to atomType (int)
	// If the atom type can not be found -1 is assigned
	public void assignNumericalType(Atom theAtom, String atomName){
  
		for(int q=0;q<atomTypeNames.length;q++) {
			if (atomTypeNames[q].equalsIgnoreCase(atomName)) {
				theAtom.type=q;
				return;
			}
		}
		theAtom.type = -1;
	}


	// This function goes through the nonbonded interactions and sets
	//  the halfNBeval and NBeval terms involving hydrogen to:
	//   0 - compute neither elect nor vdw terms
	//   1 - compute both elect and vdw terms
	//   2 - compute only elect term
	//   3 - compute only vdw term
	//  based on the input booleans
	public void setNBEval(boolean electEval, boolean vdwEval) {
	
		int evalNum = 0;
		
		if (electEval && vdwEval)
			evalNum = 1;
		else if (electEval)
			evalNum = 2;
		else if (vdwEval)
			evalNum = 3;
	
		halfNBeval = new int[numHalfNonBondedTerms];
		NBeval = new int[numberNonBonded];
		
		for(int i=0; i<numHalfNonBondedTerms; i++) {
			if (m.atom[(int)halfNonBondedTerms[i*4]].elementType.equalsIgnoreCase("H"))
				halfNBeval[i] = evalNum;
			else {
				if (m.atom[(int)halfNonBondedTerms[i*4+1]].elementType.equalsIgnoreCase("H"))
					halfNBeval[i] = evalNum;
				else
					halfNBeval[i] = 1;
			}
		}	

		for(int i = 0; i<numberNonBonded; i++) {
			if (m.atom[(int)nonBondedTerms[i*4]].elementType.equalsIgnoreCase("H"))
				NBeval[i] = evalNum;
			else {
				if (m.atom[(int)nonBondedTerms[i*4+1]].elementType.equalsIgnoreCase("H"))
					NBeval[i] = evalNum;
				else
					NBeval[i] = 1;
			}
		}	
	}
	
	// Sets the ligand molecule residue number
	public void setLigandNum(int lignum){
		ligandNum = lignum;
		
	}

//////////////////////////////////////////////////////////////////////////////////////////////////
//	This section initializes the energy calculation
//////////////////////////////////////////////////////////////////////////////////////////////////
	// This function sets up the arrays for energy evaluation
	//  it does lookup calls to getStretchParameters and such
	// It prepares terms for bond, angle, and dihedral, vdw, and electrostatic
	// Terms involving residues with energyEval == false
	//  are not included
	public void initializeCalculation(){
		
		if (debug)
			System.out.println("Starting initializeCalculation");
		
		// reset the gradient
		int natomsx3 = m.numberOfAtoms * 3;
		m.gradient = new double[natomsx3];
		for(int i=0; i<natomsx3; i++){
			m.gradient[i] = 0;
		}
		if (!m.connectivityValid)
			m.establishConnectivity(false);		
		
		initializeEVCalculation(); //initialize the calculation of the electrostatic and vdW terms
		
		if (doSolvationE) //initialize solvation energy calculation
			initializeSolvationCalculation();		
	}

	// This function sets up the arrays for energy evaluation
	//  for electrostatics and vdW only (EV)
	// Terms involving residues with energyEval == false
	//  are not included
	private void initializeEVCalculation(){

		int atom1, atom2, atom4, ix2, ix4, ix4b;
		int numberNonBondedx4;
		int atomType1, atomType2, atomType4;
		double equilibriumDistance[] = new double[1];
		double epsilon[] = new double[1];
		double smallerArray[];
		boolean evalAtom[];
		
		numberNonBonded = 0;
		numHalfNonBondedTerms = 0;

		if (debug)
			System.out.println("Starting initializeEVCalculation");

		// Build array of atom based energy evaluation booleans
		// If energy terms involving atom i are to be computed then
		//  evalAtom[i] == true
		evalAtom = new boolean[m.numberOfAtoms];
		for(int i=0;i<m.numberOfAtoms;i++){
			evalAtom[i] = false;
		}
		for(int i=0;i<m.numberOfResidues;i++){
			evalAtom = getEvalForRes(m.residue[i],evalAtom);
		}
		

		halfNonBondedTerms = new double[m.numberOf14Connections * 4];

		if (debug)
			System.out.println("Initial number of 1-4 pairs: " + m.numberOf14Connections);

		ix4 = -4;
		ix4b = -4;
		numHalfNonBondedTerms = 0;
		for(int i=0; i<m.numberOf14Connections; i++){
		  ix4 += 4;
			atom1 = m.atom[m.connected14[ix4]].moleculeAtomNumber;
			atom4 = m.atom[m.connected14[ix4 + 3]].moleculeAtomNumber;

			if (evalAtom[atom1] && evalAtom[atom4]) {
				atomType1 = m.atom[atom1].type;
				atomType4 = m.atom[atom4].type;

				ix4b += 4;
				double epsilonProduct = 0, ri = 0, rj = 0;
				if (!(getNonBondedParameters(atomType1, equilibriumDistance, epsilon)))
					System.out.println("WARNING: Could not find nb parameters for " + atom1 + " type: " + m.atom[atom1].forceFieldType);
				else {
					if((EnvironmentVars.forcefld == EnvironmentVars.FORCEFIELD.CHARMM19 
							|| EnvironmentVars.forcefld == EnvironmentVars.FORCEFIELD.CHARMM19NEUTRAL )
							&& m.atom[m.connected14[ix4]].elementType.equalsIgnoreCase("C")){
						//KER: if charmm19 then reduce C radii for 1-4 interactions
						epsilonProduct = 0.1;
						ri = 1.9;
					}
					else{
						epsilonProduct = epsilon[0];
						ri = equilibriumDistance[0];
					}
					if (!(getNonBondedParameters(atomType4, equilibriumDistance, epsilon)))
						System.out.println("WARNING: Could not find nb parameters for " + atom4 + " type: " + m.atom[atom4].forceFieldType);
					else {
						if((EnvironmentVars.forcefld == EnvironmentVars.FORCEFIELD.CHARMM19 
								|| EnvironmentVars.forcefld == EnvironmentVars.FORCEFIELD.CHARMM19NEUTRAL )
								&& m.atom[m.connected14[ix4]].elementType.equalsIgnoreCase("C")){
							//KER: if charmm19 then reduce C radii for 1-4 interactions
							epsilonProduct *= 0.1;
							rj = 1.9;
						}
						else{
							epsilonProduct *= epsilon[0];
							rj = equilibriumDistance[0];
						}
						epsilonProduct = Math.sqrt(epsilonProduct);
						// This part is 1-4 interactions which are scaled by 1/2
						double Bij = ( ri + rj ) * ( ri + rj );
						Bij = Bij * Bij * Bij;
						double Aij = Bij * Bij;
						switch(EnvironmentVars.forcefld){
							case AMBER:
								Aij *= epsilonProduct * 0.5;
								Bij *= epsilonProduct;
								break;
							case CHARMM19: 
							case CHARMM19NEUTRAL:
								Aij *= epsilonProduct;
								Bij *= epsilonProduct * 2.0;
								// Aij = (ri+rj)^12 * sqrt(ei*ej)
								// Bij = (ri+rj)^6 * sqrt(ei*ej) * 2.0
								break;
						}
						// Aij = (ri+rj)^12 * sqrt(ei*ej) * 0.5
						// Bij = (ri+rj)^6 * sqrt(ei*ej)
						halfNonBondedTerms[ix4b] = atom1;
						halfNonBondedTerms[ix4b + 1] = atom4;
						halfNonBondedTerms[ix4b + 2] = Aij;
						halfNonBondedTerms[ix4b + 3] = Bij;
						numHalfNonBondedTerms++;
					}
				}
			}
		}
		
		// Reduce the size of the halfNonBondedTerms to the size we actually used
		smallerArray = new double[numHalfNonBondedTerms * 4];
		System.arraycopy(halfNonBondedTerms, 0, smallerArray, 0, numHalfNonBondedTerms * 4);
		halfNonBondedTerms = smallerArray;
		if (debug)
			System.out.println("Final number of halfNonBondedTerms: " + numHalfNonBondedTerms);

		// make an array of 4 terms 1=atom1, 2=atom2, 3=Aij, 4=Bij
		nonBondedTerms = new double[m.numberNonBonded * 4];

		if (debug)
			System.out.println("Initial number of full nonbonded pairs: " + m.numberNonBonded);

		ix2 = -2;
		numberNonBondedx4 = -4;
		numberNonBonded = 0;
		for(int i=0; i<m.numberNonBonded; i++) {
			ix2 += 2;
			atom1 = m.atom[m.nonBonded[ix2]].moleculeAtomNumber;
			atom2 = m.atom[m.nonBonded[ix2 + 1]].moleculeAtomNumber;
			if (evalAtom[atom1] && evalAtom[atom2]) {
				atomType1 = m.atom[atom1].type;
				atomType2 = m.atom[atom2].type;
				if (!(getNonBondedParameters(atomType1, equilibriumDistance, epsilon)))
					System.out.println("WARNING: Could not find nb parameters for (at1) " + atom1 + " type: " + m.atom[atom1].forceFieldType);
				else {
					double epsilonProduct = epsilon[0];
					double ri = equilibriumDistance[0];
					if (!(getNonBondedParameters(atomType2, equilibriumDistance, epsilon)))
						System.out.println("WARNING: Could not find nb parameters for (at2) " + atom2 + " type: " + m.atom[atom2].forceFieldType);
					else {						
						epsilonProduct *= epsilon[0];
						double rj = equilibriumDistance[0];
						epsilonProduct = Math.sqrt(epsilonProduct);
						double Bij = ( ri + rj ) * ( ri + rj );
						Bij = Bij * Bij * Bij;
						double Aij = Bij * Bij * epsilonProduct;
						Bij *= epsilonProduct * 2.0;
						// Aij = (ri+rj)^12 * sqrt(ei*ej)
						// Bij = (ri+rj)^6 * sqrt(ei*ej) * 2
						numberNonBondedx4 += 4;
						nonBondedTerms[numberNonBondedx4] = atom1;
						nonBondedTerms[numberNonBondedx4 + 1] = atom2;
						nonBondedTerms[numberNonBondedx4 + 2] = Aij;
						nonBondedTerms[numberNonBondedx4 + 3] = Bij;
						numberNonBonded++;
					}
				}
			}
		}
		
		// Reduce the size of the nonBondedTerms to the size we actually used
		smallerArray = new double[numberNonBonded * 4];
		System.arraycopy(nonBondedTerms, 0, smallerArray, 0, numberNonBonded * 4);
		nonBondedTerms = smallerArray;
		
		if (debug)
			System.out.println("Final number of full nonbonded pairs: " + numberNonBonded);
	}
	
	// This function sets up the arrays for energy evaluation
	//  for solvation energies only
	// Terms involving residues with energyEval == false
	//  are not included
	// Since EEF1 handles only natural amino acids, we compute solvation 
	//	energies for proteins and ligands (natural amino acids) only, and
	//	not for cofactors. To determine if an atom belongs to a protein or a ligand,
	//	we use the isProtein flag of the Strand class. In KSParser, this flag is
	//	set to true for the protein and the ligand, but not for the cofactor
	private void initializeSolvationCalculation(){
		
		int atom1, ix6, numTerms;
		double smallerArray[];
		boolean evalAtom[];

		if (debug)
			System.out.println("Starting initializeSolvationCalculation");

		// Build array of atom based energy evaluation booleans
		// If energy terms involving atom i are to be computed then
		//  evalAtom[i] == true
		evalAtom = new boolean[m.numberOfAtoms];
		for(int i=0;i<m.numberOfAtoms;i++){
			evalAtom[i] = false;
		}
		for(int i=0;i<m.numberOfResidues;i++){			
			if (m.strand[m.residue[i].strandNumber].isProtein && !m.residue[i].cofactor) //only compute solvation energies for the protein and an AA ligand
				evalAtom = getEvalForRes(m.residue[i],evalAtom);
		}
		
		
		//Count the number of atoms to be evaluated
		numSolvationTerms = 0;
		for (int i=0; i<evalAtom.length; i++){
			if (evalAtom[i])
				numSolvationTerms++;
		}

		// Setup an array of 6 terms: 1=atom1(moleculeAtomNumber), 2=dG(ref), 3=dG(free), 4=volume, 5=lambda,
		//  6=vdW radius
		solvationTerms = new double[numSolvationTerms * 6];

		if (debug)
			System.out.println("Initial number of solvation terms: " + numSolvationTerms);

		ix6 = -6;
		numTerms = 0;
		for(int i=0; i<m.numberOfAtoms; i++){
			atom1 = m.atom[i].moleculeAtomNumber;

			if (evalAtom[atom1]) {
				
				if (!m.atom[atom1].elementType.equalsIgnoreCase("H")){//solvation terms do not include H
				
					ix6 += 6;
					
					double dGref[] = new double[1];
					double dGfree[] = new double[1];
					double atVolume[] = new double[1];
					double lambda[] = new double[1];
					double vdWradiusExt[] = new double[1]; //extended vdWradius (uses the EEF1 parameters)
					
					if (!(eef1parms.getSolvationParameters(atom1,dGref,dGfree,atVolume,lambda,vdWradiusExt))){
						System.out.println("WARNING: Could not find solvation parameters for atom: " + atom1+" ("+m.atom[atom1].name+") res: "+m.atom[atom1].moleculeResidueNumber+" ("+m.residue[m.atom[atom1].moleculeResidueNumber].name+")");
						System.exit(1);
					}
					else {
						
						solvationTerms[ix6] = atom1;
						solvationTerms[ix6 + 1] = dGref[0];
						solvationTerms[ix6 + 2] = dGfree[0];
						solvationTerms[ix6 + 3] = atVolume[0];
						solvationTerms[ix6 + 4] = lambda[0];
						solvationTerms[ix6 + 5] = vdWradiusExt[0];
						numTerms++;
					}
				}
			}
		}
		
		// Shrink the dihedralAngleTerms array down
		smallerArray = new double[numTerms * 6];
		System.arraycopy(solvationTerms,0,smallerArray,0,numTerms*6);
		solvationTerms = smallerArray;
		numSolvationTerms = numTerms;

		if (debug)
			System.out.println("Final number of solvation terms: " + numSolvationTerms);
		
		//Determine which pairs of atoms can be excluded from the solvation energy computation
		solvExcludePairs = new boolean[numSolvationTerms][numSolvationTerms];
		for (int i=0; i<solvExcludePairs.length; i++){
			
			int atomi = (int)solvationTerms[i*6];
			
			for (int j=0; j<solvExcludePairs.length; j++){
				
				if (i!=j){
					
					int atomj = (int)solvationTerms[j*6];
					if ( (!m.are12connected(atomi,atomj)) && (!m.are13connected(atomi,atomj)) ) //(not 1-2) and (not 1-3) connected
						solvExcludePairs[i][j] = false;
					else
						solvExcludePairs[i][j] = true;
				}
			}
		}
	}
	
	//Determines which atoms of residue res (molecule-relative numbering) should be included in the energy computation
	private boolean [] getEvalForRes(Residue res, boolean evalAtom[]){
		
		
		if (res.getEnergyEvalSC()){
			if (res.ffAssigned){
				for(int j=0;j<res.numberOfAtoms;j++){
					if (!res.atom[j].getIsBBatom())
						evalAtom[res.atom[j].moleculeAtomNumber] = true;
				}
			}
			else {
				System.out.println("ERROR: the force field for residue "+res.fullName+" is not assigned");
				System.exit(1);
			}
		}
		if (res.getEnergyEvalBB()){
			if (res.ffAssigned){
				for(int j=0;j<res.numberOfAtoms;j++){
					if (res.atom[j].getIsBBatom())
						evalAtom[res.atom[j].moleculeAtomNumber] = true;
				}
			}
			else {
				System.out.println("ERROR: the force field for residue "+res.fullName+" is not assigned");
				System.exit(1);
			}
		}
		
		return evalAtom;
	}
//////////////////////////////////////////////////////////////////////////////////////////////////
	

//////////////////////////////////////////////////////////////////////////////////////////////////
//	Sets up partial atom arrays: used by SimpleMaximizer to quickly compute energies on given subsets of
// 		the full atoms arrays (nonBonded, halfNonBonded, dihedralAngleTerms, etc.)	
//////////////////////////////////////////////////////////////////////////////////////////////////	
	//Sets up the partial arrays
	public void setupPartialArrays(int numRows, int maxNumColumns, int atomList[][],int numColumns[]){
		
		setupPartialNonBondedArrays(numRows, maxNumColumns, atomList, numColumns); //setup nonbonded arrays
		
		if (doSolvationE) //setup dihedral arrays
			setupPartialSolvationArrays(numRows, maxNumColumns, atomList, numColumns);
	}
	
	// Sets up local datastructures to hold lists of the halfNonBonded and nonBonded
	//  terms that include atoms in atomList. Each row of atomList is a different
	//  subset. atomList has numRows rows and numCol columns.
	// For example this is used by SimpleMinimizer to quickly compute nonbonded
	//  energies in conjunction with calculateEVEnergyPartWithArrays
	private void setupPartialNonBondedArrays(int numRows, int maxNumColumns, int atomList[][],
		int numColumns[]){
	
		int ix4 = -4;
		int tempCount = 0, tempIndx = 0;
		int atomi = 0, atomj = 0;
		
		numPartHalfNonBonded = new int[numRows];
		numPartNonBonded = new int[numRows];
		
		partHalfNonBonded = new double[numRows][];
		partNonBonded = new double[numRows][];
			// In the worst case each atom in a column of atomList is involved with
			//  every other atom in the molecule
		partHalfNBeval = new int[numRows][];
		partNBeval = new int[numRows][];
		
		for(int i=0; i<numRows;i++){
			partHalfNonBonded[i] = new double[numColumns[i]*m.numberOfAtoms*4];
			partNonBonded[i] = new double[numColumns[i]*m.numberOfAtoms*4];
				// In the worst case each atom in a column of atomList is involved with
				//  every other atom in the molecule
			partHalfNBeval[i] = new int[numColumns[i]*m.numberOfAtoms];
			partNBeval[i] = new int[numColumns[i]*m.numberOfAtoms];
		}
	
		for(int q=0;q<numRows;q++){
			int[] tempAtomList = new int[m.numberOfAtoms];
			for(int i=0;i<numColumns[q];i++)
				tempAtomList[atomList[q][i]] = 1;
			tempCount = 0;
			ix4 = -4;
			for(int i=0; i<numHalfNonBondedTerms; i++) {
				ix4 += 4;
				atomi = (int)halfNonBondedTerms[ix4];
				atomj = (int)halfNonBondedTerms[ix4 + 1];
				if ((tempAtomList[atomi] + tempAtomList[atomj]) > 0){
					tempIndx = tempCount * 4;
					partHalfNonBonded[q][tempIndx] = halfNonBondedTerms[ix4];
					partHalfNonBonded[q][tempIndx+1] = halfNonBondedTerms[ix4 + 1];
					partHalfNonBonded[q][tempIndx+2] = halfNonBondedTerms[ix4 + 2];
					partHalfNonBonded[q][tempIndx+3] = halfNonBondedTerms[ix4 + 3];
					partHalfNBeval[q][tempCount] = halfNBeval[i];
					tempCount++;
				}
			}
			numPartHalfNonBonded[q] = tempCount;
			
			tempCount = 0;
			ix4 = -4;
			for(int i=0; i<numberNonBonded; i++) {
				ix4 += 4;
				atomi = (int)nonBondedTerms[ix4];
				atomj = (int)nonBondedTerms[ix4 + 1];
				if ((tempAtomList[atomi] + tempAtomList[atomj]) > 0){
					tempIndx = tempCount * 4;
					partNonBonded[q][tempIndx] = nonBondedTerms[ix4];
					partNonBonded[q][tempIndx+1] = nonBondedTerms[ix4 + 1];
					partNonBonded[q][tempIndx+2] = nonBondedTerms[ix4 + 2];
					partNonBonded[q][tempIndx+3] = nonBondedTerms[ix4 + 3];
					partNBeval[q][tempCount] = NBeval[i];
					tempCount++;
				}
			}
			numPartNonBonded[q] = tempCount;
		}
	
	}
	
	// Sets up local datastructures to hold lists of the solvation
	//  terms that include atoms in atomList. Each row of atomList is a different
	//  subset. atomList has numRows rows and numCol columns.
	// See setupPartialNonBondedArrays() for additional comments
	private void setupPartialSolvationArrays(int numRows, int maxNumColumns, int atomList[][],
		int numColumns[]){
		
		int ix6 = -6;
		int tempCount = 0, tempIndx = 0;
		int atomi = 0;
		
		numPartSolv = new int[numRows];
		partSolv = new double[numRows][maxNumColumns*m.numberOfAtoms*7];		
	
		for(int q=0;q<numRows;q++){
			int[] tempAtomList = new int[m.numberOfAtoms];
			for(int i=0;i<numColumns[q];i++)
				tempAtomList[atomList[q][i]] = 1;
			tempCount = 0;
			ix6 = -6;
			for(int i=0; i<numSolvationTerms; i++) {
				ix6 += 6;
				atomi = (int)solvationTerms[ix6];
				if ((tempAtomList[atomi]) > 0){
					tempIndx = tempCount * 7;
					partSolv[q][tempIndx] = solvationTerms[ix6];
					partSolv[q][tempIndx+1] = solvationTerms[ix6 + 1];
					partSolv[q][tempIndx+2] = solvationTerms[ix6 + 2];
					partSolv[q][tempIndx+3] = solvationTerms[ix6 + 3];
					partSolv[q][tempIndx+4] = solvationTerms[ix6 + 4];
					partSolv[q][tempIndx+5] = solvationTerms[ix6 + 5];
					partSolv[q][tempIndx+6] = i;
					tempCount++;
				}
			}
			numPartSolv[q] = tempCount;
		}		
	}
///////////////////////////////////////////////////////////////////////////////////////////////
	

//////////////////////////////////////////////////////////////////////////////////////////////////
//	This section calculates the energy of the given system
//////////////////////////////////////////////////////////////////////////////////////////////////	
	//Calculates the total energy of the system specified by coordinates[];
	//Depending on the flags, different types of energies are included/excluded
	//If (curIndex!=-1), then precomputed partial arrays  are used to quickly compute the specified energy terms
	public double [] calculateTotalEnergy(float coordinates[], int curIndex){
		
		double energyTerms[] = new double[4]; //total, electrostatics, vdW, and solvation
		for (int i=0; i<energyTerms.length; i++)
			energyTerms[i] = 0.0;
		
		
		calculateEVEnergy(coordinates,curIndex,energyTerms); //compute electrostatic and vdW energies
		
		if (doSolvationE) //compute solvation energies
			 calculateSolvationEnergy(coordinates,curIndex,energyTerms);
		
		//compute total energy (electrostatics + vdW + solvation)
		energyTerms[0] = energyTerms[1] + energyTerms[2] + energyTerms[3];
		
		return energyTerms;
	}

	// This function calculates the electrostatic and vdw (EV) energy of a system
	// Energy values are 'R'eturned in EenergyR and VenergyR
	// vdwMultiplier is the soft-potential multiplier, if equal to 1.0 it has no effect,
	//  values <1.0 allow for slight overpacking
	// ligand number is the residue number that corresponds to the ligand
	//  this is used to break the energy terms down into P-P and P-L
	// ligand number is 1 based numbering
	// EenergyR[] and VenergyR[] have the following values
	//   XenergyR[0] = total energy
	//   XenergyR[1] = P-P energy
	//   XenergyR[2] = P-L energy
	//   XenergyR[3] = L-L energy
	// If ligandNum == -1 then only the total energy is returned
	// Makes use of the halfNBeval and NBeval arrays
	//   0 - compute neither elect nor vdw terms
	//   1 - compute both elect and vdw terms
	//   2 - compute only elect term
	//   3 - compute only vdw term
	private void calculateEVEnergy(float coordinates[], int curIndex, double energyTerms[]){

		int atomix3, atomjx3, atomi, atomj;
		int ix4;
		double rij, rij2, rij6, rij12, coulombTerm, vdwTerm;
		double rijx, rijy, rijz;
		double chargei, chargej, Aij, Bij;
		double coulombFactor;
		float Eenergy[], Venergy[];
		float Amult, Bmult;
		
		int numHalfNBterms = 0; int numNBterms = 0;
		double halfNBterms[] = null; double nbTerms[] = null;
		int halfNBev[] = null; int nbEv[] = null;
		
		if (curIndex==-1){ //full energy is computed
			numHalfNBterms = numHalfNonBondedTerms;
			halfNBterms = halfNonBondedTerms;
			halfNBev = halfNBeval;
			numNBterms = numberNonBonded;
			nbTerms = nonBondedTerms;
			nbEv = NBeval;
		}
		else { //partial energy is computed, based on flexible residue curIndex
			numHalfNBterms = numPartHalfNonBonded[curIndex];
			halfNBterms = partHalfNonBonded[curIndex];
			halfNBev = partHalfNBeval[curIndex];
			numNBterms = numPartNonBonded[curIndex];
			nbTerms = partNonBonded[curIndex];
			nbEv = partNBeval[curIndex];
		}
		
		Eenergy = new float[4];
		Venergy = new float[4];
		for(int i=0;i<4;i++) {
			Eenergy[i] = 0.0f;
			Venergy[i] = 0.0f;
		}

		// Note: Bmult = vdwMultiplier^6 and Amult = vdwMultiplier^12
		Bmult = vdwMultiplier * vdwMultiplier;
		Bmult = Bmult*Bmult*Bmult;
		Amult = Bmult*Bmult;

		// half non-bonded terms
		ix4 = -4;
		// 1-4 electrostatic terms are scaled by 1/1.2
		switch(EnvironmentVars.forcefld){
			case AMBER:
				coulombFactor = (constCoulomb/1.2) / (dielectric);
				break;
			case CHARMM19:
			case CHARMM19NEUTRAL:
				coulombFactor = (constCoulomb * 0.4) / (dielectric);
				break;
			default:
				coulombFactor = 0;
				System.out.println("FORCEFIELD NOT RECOGNIZED!!!");
				System.exit(0);
				break;
		}
		
		double tmpCoulFact;
		for(int i=0; i<numHalfNBterms; i++) {
			ix4 += 4;
			atomi = (int)halfNBterms[ix4];
			atomj = (int)halfNBterms[ix4 + 1];
			Aij = halfNBterms[ix4 + 2] * Amult;
			Bij = halfNBterms[ix4 + 3] * Bmult;
			chargei = m.atom[atomi].charge;
			chargej = m.atom[atomj].charge;
			atomix3 = atomi * 3;
			atomjx3 = atomj * 3;
			rijx = coordinates[atomix3] - coordinates[atomjx3];
			rijy = coordinates[atomix3 + 1] - coordinates[atomjx3 + 1];
			rijz = coordinates[atomix3 + 2] - coordinates[atomjx3 + 2];

			rij2 = rijx * rijx + rijy * rijy + rijz * rijz;
			rij = Math.sqrt(rij2);
			rij6 = rij2 * rij2 * rij2;
			rij12 = rij6 * rij6;
			
			//coulombFactor = (constCoulomb/1.2) / (dielectric);
			tmpCoulFact = coulombFactor;
			if (distDepDielect) //distance-dependent dielectric
				tmpCoulFact /= rij;
	
			coulombTerm = (chargei * chargej * tmpCoulFact) / rij;
			vdwTerm = Aij / rij12 - Bij / rij6;

			// This is not the fastest way to do this, but based on the
			//  halfNBeval array either the elect or vdw energies might
			//  not be counted
			if (halfNBev[i] == 2)
				vdwTerm = 0.0;
			else if (halfNBev[i] == 3)
				coulombTerm = 0.0;
			else if (halfNBev[i] == 0) {
				vdwTerm = 0.0;
				coulombTerm = 0.0;
			}
			Eenergy[0] += coulombTerm;
			Venergy[0] += vdwTerm;
			if (m.atom[atomi].moleculeResidueNumber == ligandNum) {
				if (m.atom[atomj].moleculeResidueNumber == ligandNum) {
					// L-L
					Eenergy[3] += coulombTerm;
					Venergy[3] += vdwTerm;
				}
				else {
					// P-L
					Eenergy[2] += coulombTerm;
					Venergy[2] += vdwTerm;
				}
			}
			else if (m.atom[atomj].moleculeResidueNumber == ligandNum) {
				// P-L
				Eenergy[2] += coulombTerm;
				Venergy[2] += vdwTerm;
			}
			else {
			  // P-P
				Eenergy[1] += coulombTerm;
				Venergy[1] += vdwTerm;
			}
		}

		ix4 = -4;
		// The full nonbonded electrostatic terms are NOT scaled down by 1/1.2
		coulombFactor = constCoulomb / (dielectric);
		for(int i=0; i<numNBterms; i++) {
			ix4 += 4;
			atomi = (int)nbTerms[ ix4 ];
			atomj = (int)nbTerms[ ix4 + 1 ];
				
			Aij = nbTerms[ ix4 + 2 ] * Amult;
			Bij = nbTerms[ ix4 + 3 ] * Bmult;
			chargei = m.atom[ atomi ].charge;
			chargej = m.atom[ atomj ].charge;
			atomix3 = atomi * 3;
			atomjx3 = atomj * 3;
			rijx = coordinates[ atomix3 ] - coordinates[ atomjx3 ];
			rijy = coordinates[ atomix3 + 1 ] - coordinates[ atomjx3 + 1 ];
			rijz = coordinates[ atomix3 + 2 ] - coordinates[ atomjx3 + 2 ];
			rij2 = rijx * rijx + rijy * rijy + rijz * rijz;
			rij = Math.sqrt( rij2 );
			rij6 = rij2 * rij2 * rij2;
			rij12 = rij6 * rij6;
			
			//coulombFactor = constCoulomb / (dielectric);
			tmpCoulFact = coulombFactor;
			if (distDepDielect) //distance-dependent dielectric
				tmpCoulFact /= rij;

			coulombTerm = (chargei * chargej * tmpCoulFact) / rij;
			vdwTerm = Aij / rij12 - Bij / rij6;

			// This is not the fastest way to do this, but based on the
			//  NBeval array either the elect or vdw energies might
			//  not be counted
			if (nbEv[i] == 2)
				vdwTerm = 0.0;
			else if (nbEv[i] == 3)
				coulombTerm = 0.0;
			else if (nbEv[i] == 0) {
				vdwTerm = 0.0;
				coulombTerm = 0.0;
			}
			Eenergy[0] += coulombTerm;
			Venergy[0] += vdwTerm;
			if (m.atom[atomi].moleculeResidueNumber == ligandNum) {
				if (m.atom[atomj].moleculeResidueNumber == ligandNum) {
					// L-L
					Eenergy[3] += coulombTerm;
					Venergy[3] += vdwTerm;
				}
				else {
					// P-L
					Eenergy[2] += coulombTerm;
					Venergy[2] += vdwTerm;
				}
			}
			else if (m.atom[atomj].moleculeResidueNumber == ligandNum) {
				// P-L
				Eenergy[2] += coulombTerm;
				Venergy[2] += vdwTerm;
			}
			else {
			  // P-P
				Eenergy[1] += coulombTerm;
				Venergy[1] += vdwTerm;
			}
		}
		
		//store computed energies
		energyTerms[1] = Eenergy[0]; //electrostatics
		energyTerms[2] = Venergy[0]; //vdW
	}
	
	//Calculates the solvation energies for the system with given coordinates[]
	private void calculateSolvationEnergy(float coordinates[], int curIndex, double energyTerms[]){
		
		if (curIndex==-1) //all residues included
			calculateSolvationEnergyFull(coordinates,energyTerms);
		else //only residue curIndex included (partial matrices used)
			calculateSolvationEnergyPart(coordinates,curIndex,energyTerms);
	}
	
	//Calculates the solvation energies for the system with given coordinates[]
	private void calculateSolvationEnergyFull(float coordinates[], double energyTerms[]){
		
		double energy = 0.0;
		int atomix3, atomjx3, atomi, atomj;
		double rij, rij2;
		double rijx, rijy, rijz;
		int indMult = 6;
		
		int numSolvTerms = 0;
		double solvTerms[] = null;
		
		numSolvTerms = numSolvationTerms;
		solvTerms = solvationTerms;
		
		for ( int i = 0; i < numSolvTerms; i++ ){

			atomi = (int)solvTerms[ i*indMult ];
			atomix3 = atomi * 3;
			
			energy += solvTerms[i*indMult+1]; //dGi(ref)
			
			double dGi_free = solvTerms[i*indMult+2]; //dGi(free)
			double V_i = solvTerms[i*indMult+3]; //Vi
			double lambda_i = solvTerms[i*indMult+4]; //lambdai
			double vdWr_i = solvTerms[i*indMult+5]; //vdWri
			
			for (int j=i+1; j<numSolvationTerms; j++){ //the pairwise solvation energies
				
				atomj = (int)solvationTerms[j*indMult];
				atomjx3 = atomj*3;
				
				//atoms 1 or 2 bonds apart are excluded from each other's calculation of solvation free energy
				if (!solvExcludePairs[i][j]){
					
					rijx = coordinates[ atomix3 ] - coordinates[ atomjx3 ];
					rijy = coordinates[ atomix3 + 1 ] - coordinates[ atomjx3 + 1 ];
					rijz = coordinates[ atomix3 + 2 ] - coordinates[ atomjx3 + 2 ];
					rij2 = rijx * rijx + rijy * rijy + rijz * rijz;
					rij = Math.sqrt( rij2 ); //distance between the two atoms
					
					if (rij < solvCutoff){
						
						double dGj_free = solvationTerms[j*indMult+2]; //dGj(free)
						double V_j = solvationTerms[j*indMult+3]; //Vj
						double lambda_j = solvationTerms[j*indMult+4]; //lambdaj
						double vdWr_j = solvationTerms[j*indMult+5]; //vdWrj
					
						double coeff = 1/(4*Math.PI*Math.sqrt(Math.PI));
						
						double Xij = (rij-vdWr_i)/lambda_i;
						double Xji = (rij-vdWr_j)/lambda_j;
						
						energy -= ( (2 * coeff * dGi_free * Math.exp(-Xij*Xij) * V_j) / (lambda_i * rij2)
									+ (2 * coeff * dGj_free * Math.exp(-Xji*Xji) * V_i) / (lambda_j * rij2) );
					}
				}
			}
		}
		
		//store computed energy
		energyTerms[3] = solvScale*energy; //solvation
	}
	
	//Calculates the solvation energies for the system with given coordinates[]
	private void calculateSolvationEnergyPart(float coordinates[], int curIndex, double energyTerms[]){
		
		double energy = 0.0;
		int atomix3, atomjx3, atomi, atomj;
		double rij, rij2;
		double rijx, rijy, rijz;
		int indMult = 7;
		int startInd;
		
		int numSolvTerms = 0;
		double solvTerms[] = null;
		
		numSolvTerms = numPartSolv[curIndex];
		solvTerms = partSolv[curIndex];
		
		for ( int i = 0; i < numSolvTerms; i++ ){

			atomi = (int)solvTerms[ i*indMult ];
			atomix3 = atomi * 3;
			
			energy += solvTerms[i*indMult+1]; //dGi(ref)
			
			double dGi_free = solvTerms[i*indMult+2]; //dGi(free)
			double V_i = solvTerms[i*indMult+3]; //Vi
			double lambda_i = solvTerms[i*indMult+4]; //lambdai
			double vdWr_i = solvTerms[i*indMult+5]; //vdWri
			
			startInd = (int)solvTerms[i*indMult+6];
			
			for (int j=0; j<numSolvationTerms; j++){ //the pairwise solvation energies
				
				atomj = (int)solvationTerms[j*6];
				atomjx3 = atomj*3;
				
				boolean comp = true;
				if (m.atom[atomi].moleculeResidueNumber==m.atom[atomj].moleculeResidueNumber){
					if (j<=startInd)
						comp = false;
				}
				
				if (comp){
					//atoms 1 or 2 bonds apart are excluded from each other's calculation of solvation free energy
					if (!solvExcludePairs[startInd][j]){
						
						rijx = coordinates[ atomix3 ] - coordinates[ atomjx3 ];
						rijy = coordinates[ atomix3 + 1 ] - coordinates[ atomjx3 + 1 ];
						rijz = coordinates[ atomix3 + 2 ] - coordinates[ atomjx3 + 2 ];
						rij2 = rijx * rijx + rijy * rijy + rijz * rijz;
						rij = Math.sqrt( rij2 ); //distance between the two atoms
						
						if (rij < solvCutoff){
							
							double dGj_free = solvationTerms[j*6+2]; //dGj(free)
							double V_j = solvationTerms[j*6+3]; //Vj
							double lambda_j = solvationTerms[j*6+4]; //lambdaj
							double vdWr_j = solvationTerms[j*6+5]; //vdWrj
						
							double coeff = 1/(4*Math.PI*Math.sqrt(Math.PI));
							
							double Xij = (rij-vdWr_i)/lambda_i;
							double Xji = (rij-vdWr_j)/lambda_j;
							
							energy -= ( (2 * coeff * dGi_free * Math.exp(-Xij*Xij) * V_j) / (lambda_i * rij2)
										+ (2 * coeff * dGj_free * Math.exp(-Xji*Xji) * V_i) / (lambda_j * rij2) );
						}
					}
				}
			}
		}
		
		//store computed energy
		energyTerms[3] = solvScale*energy; //solvation
	}
//////////////////////////////////////////////////////////////////////////////////////////////////
	

//////////////////////////////////////////////////////////////////////////////////////////////////
//	This section computes the energy gradient
//////////////////////////////////////////////////////////////////////////////////////////////////
	public void calculateGradient(){
		System.out.println("ERROR: the function calculateGradient(int index) must be used.");
		System.exit(1);
	}
	
	//Computes the gradient of the different energy terms;
	//The computed gradient is in the molecule's gradient member variable
	//The parameter curIndex specifies the row in the partial arrays
	//		(the corresponding flexible residue);
	//If (curIndex==-1), then the full gradient is computed
	public void calculateGradient(int curIndex){
		
		// clear the gradient		
		m.gradient = new double[m.numberOfAtoms * 3];
		for(int i=0; i<m.numberOfAtoms*3;i ++){
			m.gradient[i]=0;
		}
		
		calculateEVGradient(curIndex); //compute electrostatic and vdW energies
		
		if (doSolvationE) //compute solvation energies
			calculateSolvationGradient(curIndex);
	}
	
	// This code computes the gradient of the electrostatic and vdw energy terms
	// The computed gradient is in the molecule's gradient member variable
	private void calculateEVGradient(int curIndex){
		
		int ix4;
		double coulombFactor;
		int atomi, atomj, atomix3, atomjx3;
		double Aij, Bij, rij6, rij7, rij8, rij14;
		double chargei, chargej, coulombTerm;
		double rijx, rijy, rijz, rij2, rij, rij3;
		double term1, term2, term3;
		double forceix, forceiy, forceiz, forcejx, forcejy, forcejz;
		
		int numHalfNBterms = 0; int numNBterms = 0;
		double halfNBterms[] = null; double nbTerms[] = null;
		
		if (curIndex==-1){ //full gradient is computed
			numHalfNBterms = numHalfNonBondedTerms;
			halfNBterms = halfNonBondedTerms;
			numNBterms = numberNonBonded;
			nbTerms = nonBondedTerms;
		}
		else { //partial gradient is computed, based on flexible residue curIndex
			numHalfNBterms = numPartHalfNonBonded[curIndex];
			halfNBterms = partHalfNonBonded[curIndex];
			numNBterms = numPartNonBonded[curIndex];
			nbTerms = partNonBonded[curIndex];
		}

		// Note: Bmult = vdwMultiplier^6 and Amult = vdwMultiplier^12
		float Bmult; float Amult;
		Bmult = vdwMultiplier * vdwMultiplier;
		Bmult = Bmult*Bmult*Bmult;
		Amult = Bmult*Bmult;
		
		// compute gradient for 1/2 non-bonded terms
		ix4 = -4;
		//coulombFactor = (constCoulomb/2.0) / (dielectric);
		//KER: Made change to 1.2
		switch(EnvironmentVars.forcefld){
			case AMBER: 
				coulombFactor = (constCoulomb/1.2) / dielectric;
				break;
			case CHARMM19:
			case CHARMM19NEUTRAL:
				coulombFactor = (constCoulomb*0.4) / dielectric;
				break;
			default:
				coulombFactor = 0;
				System.out.println("FORCEFIELD NOT RECOGNIZED!!!");
				System.exit(0);
				break;
		}
		double tmpCoulFact;
		for (int i=0; i<numHalfNBterms; i++){
			ix4 += 4;
			atomi = (int)halfNBterms[ix4];
			atomj = (int)halfNBterms[ix4 + 1];
			Aij = halfNBterms[ix4 + 2]*Amult;
			Bij = halfNBterms[ix4 + 3]*Bmult;
			chargei = m.atom[atomi].charge;
			chargej = m.atom[atomj].charge;
			atomix3 = atomi * 3;
			atomjx3 = atomj * 3;
			rijx = m.actualCoordinates[atomix3] - 
				m.actualCoordinates[atomjx3];
			rijy = m.actualCoordinates[atomix3 + 1] - 
				m.actualCoordinates[atomjx3 + 1];
			rijz = m.actualCoordinates[atomix3 + 2] - 
				m.actualCoordinates[atomjx3 + 2];
			rij2 = rijx * rijx + rijy * rijy + rijz * rijz;
			if (rij2 < 1.0e-2)
				rij2 = 1.0e-2;	
			rij = Math.sqrt(rij2);
			rij3 = rij2 * rij;
			rij6 = rij3 * rij3;
			rij7 = rij6 * rij;
			rij8 = rij7 * rij;
			rij14 = rij7 * rij7;
			
			tmpCoulFact = coulombFactor;
			if (distDepDielect) //distance-dependent dielectric
				tmpCoulFact = (tmpCoulFact * 2) / rij;
			
			coulombTerm = (chargei * chargej * tmpCoulFact) / rij3;
			term1 = 12 * Aij / rij14;
			term2 = 6 * Bij / rij8;
			term3 = term1 - term2 + coulombTerm;
			forceix = term3 * rijx;
			forceiy = term3 * rijy;
			forceiz = term3 * rijz;
			forcejx = -forceix;
			forcejy = -forceiy;
			forcejz = -forceiz;
			m.gradient[atomix3] -= forceix;
			m.gradient[atomix3 + 1] -= forceiy;
			m.gradient[atomix3 + 2] -= forceiz;
			m.gradient[atomjx3] -= forcejx;
			m.gradient[atomjx3 + 1] -= forcejy;
			m.gradient[atomjx3 + 2] -= forcejz;
		}

		// Full non bonded terms
		ix4 = -4;
		coulombFactor = constCoulomb / (dielectric);
		for(int i=0; i<numNBterms; i++){
			ix4 += 4;
			atomi = (int)nbTerms[ix4];
			atomj = (int)nbTerms[ix4 + 1];
			Aij = nbTerms[ix4 + 2]*Amult;
			Bij = nbTerms[ix4 + 3]*Bmult;
			chargei = m.atom[atomi].charge;
			chargej = m.atom[atomj].charge;
			atomix3 = atomi * 3;
			atomjx3 = atomj * 3;
			rijx = m.actualCoordinates[atomix3] - 
				m.actualCoordinates[atomjx3];
			rijy = m.actualCoordinates[atomix3 + 1] - 
				m.actualCoordinates[atomjx3 + 1];
			rijz = m.actualCoordinates[atomix3 + 2] - 
				m.actualCoordinates[atomjx3 + 2];
			rij2 = rijx * rijx + rijy * rijy + rijz * rijz;
			if (rij2 < 1.0e-2)
				rij2 = 1.0e-2;	
			rij = Math.sqrt(rij2);
			rij3 = rij2 * rij;
			rij6 = rij3 * rij3;
			rij7 = rij6 * rij;
			rij8 = rij7 * rij;
			rij14 = rij7 * rij7;
			
			tmpCoulFact = coulombFactor;
			if (distDepDielect) //distance-dependent dielectric
				tmpCoulFact = (tmpCoulFact * 2) / rij;
			
			coulombTerm = (chargei * chargej * coulombFactor) / rij3;
			term1 = 12 * Aij / rij14;
			term2 = 6 * Bij / rij8;
			term3 = term1 - term2 + coulombTerm;
			forceix = term3 * rijx;
			forceiy = term3 * rijy;
			forceiz = term3 * rijz;
			forcejx = -forceix;
			forcejy = -forceiy;
			forcejz = -forceiz;
			m.gradient[atomix3] -= forceix;
			m.gradient[atomix3 + 1] -= forceiy;
			m.gradient[atomix3 + 2] -= forceiz;
			m.gradient[atomjx3] -= forcejx;
			m.gradient[atomjx3 + 1] -= forcejy;
			m.gradient[atomjx3 + 2] -= forcejz;
		}  
	}
	
	//Computes the gradient for the solvation energy term;
	//The computed gradient is in the molecules gradient member variable
	private void calculateSolvationGradient(int curIndex) {
		
		double forceix, forceiy, forceiz;
		int atomix3, atomjx3, atomi, atomj;
		double rij, rij2, rij3;
		double rijx, rijy, rijz;
		double tempTerm_i;
		int indMult = 0;
		
		int numSolvTerms = 0;
		double solvTerms[] = null;
		
		if (curIndex==-1){ //full energy is computed
			numSolvTerms = numSolvationTerms;
			solvTerms = solvationTerms;
			indMult = 6;
		}
		else { //partial energy is computed, based on flexible residue curIndex
			numSolvTerms = numPartSolv[curIndex];
			solvTerms = partSolv[curIndex];
			indMult = 7;
		}
		
		for ( int i = 0; i < numSolvTerms; i++ ){

			atomi = (int)solvTerms[ i*indMult ];
			atomix3 = atomi * 3;
			
			double dGi_free = solvTerms[i*indMult+2]; //dGi(free)
			double V_i = solvTerms[i*indMult+3]; //Vi
			double lambda_i = solvTerms[i*indMult+4]; //lambdai
			double vdWr_i = solvTerms[i*indMult+5]; //vdWri
			
			int startInd = i;
			if (curIndex!=-1)
				startInd = (int)solvTerms[i*indMult+6];
			
			forceix = 0.0;forceiy = 0.0; forceiz = 0.0;
			for (int j=0; j<numSolvationTerms; j++){ //the pairwise solvation energies
				
				if (j!=startInd){
				
					atomj = (int)solvTerms[j*6];
					atomjx3 = atomj*3;
					
					//atoms 1 or 2 bonds apart are excluded from each other's calculation of solvation free energy
					if (!solvExcludePairs[startInd][j]) {
						
						rijx = m.actualCoordinates[ atomix3 ] - m.actualCoordinates[ atomjx3 ];
						rijy = m.actualCoordinates[ atomix3 + 1 ] - m.actualCoordinates[ atomjx3 + 1 ];
						rijz = m.actualCoordinates[ atomix3 + 2 ] - m.actualCoordinates[ atomjx3 + 2 ];
						rij2 = rijx * rijx + rijy * rijy + rijz * rijz;
						rij = Math.sqrt( rij2 ); //distance between the two atoms
						rij3 = rij2 * rij;
						
						if (rij < solvCutoff){
							
							double dGj_free = solvationTerms[j*6+2]; //dGj(free)
							double V_j = solvationTerms[j*6+3]; //Vj
							double lambda_j = solvationTerms[j*6+4]; //lambdaj
							double vdWr_j = solvationTerms[j*6+5]; //vdWrj
						
							double coeff = 1/(Math.PI*Math.sqrt(Math.PI));
							
							double Xij = (rij-vdWr_i)/lambda_i;
							double Xji = (rij-vdWr_j)/lambda_j;
							
							double Vj_coeff = Xij/lambda_i + 1/rij;
							double Vi_coeff = Xji/lambda_j + 1/rij;
							
							tempTerm_i = ( (coeff * dGi_free * Math.exp(-Xij*Xij) * Vj_coeff * V_j) / (lambda_i * rij3)
										+ (coeff * dGj_free * Math.exp(-Xji*Xji) * Vi_coeff * V_i) / (lambda_j * rij3) ) ;
							
							forceix += tempTerm_i * rijx;
							forceiy += tempTerm_i * rijy;
							forceiz += tempTerm_i * rijz;
						}
					}
				}
			}		
			
			m.gradient[ atomix3 ] += solvScale*forceix;
			m.gradient[ atomix3 + 1 ] += solvScale*forceiy;
			m.gradient[ atomix3 + 2 ] += solvScale*forceiz;
		}
	}
//////////////////////////////////////////////////////////////////////////////////////////////////
	
	//Doubles the size of the a[] String array
	private String [] doubleArraySize(String a[]){
		String tmp[] = new String[a.length*2];
		System.arraycopy(a, 0, tmp, 0, a.length);
		return tmp;
	}
	
	//Doubles the size of the a[] float array
	private float [] doubleArraySize(float a[]){
		float tmp[] = new float[a.length*2];
		System.arraycopy(a, 0, tmp, 0, a.length);
		return tmp;
	}
	
	//Doubles the size of the a[] int array
	private int [] doubleArraySize(int a[]){
		int tmp[] = new int[a.length*2];
		System.arraycopy(a, 0, tmp, 0, a.length);
		return tmp;
	}
	
	//Doubles the size (first dimension only) of the a[] int array
	private int [][] doubleArraySize(int a[][]){
		int tmp[][] = new int[a.length*2][];
		for (int i=0; i<a.length; i++){
			tmp[i] = new int[a[i].length];
			System.arraycopy(a[i], 0, tmp[i], 0, a[i].length);
		}
		return tmp;
	}
	
	//Reduce the a[] String array to keep only the first newSize elements
	private String [] reduceArraySize(String a[], int newSize){
		String tmp[] = new String[newSize];
		System.arraycopy(a, 0, tmp, 0, tmp.length);
		return tmp;
	}
	
	//Reduce the a[] float array to keep only the first newSize elements
	private float [] reduceArraySize(float a[], int newSize){
		float tmp[] = new float[newSize];
		System.arraycopy(a, 0, tmp, 0, tmp.length);
		return tmp;
	}
	
	//Reduce the a[] int array to keep only the first newSize elements
	private int [] reduceArraySize(int a[], int newSize){
		int tmp[] = new int[newSize];
		System.arraycopy(a, 0, tmp, 0, tmp.length);
		return tmp;
	}
	
	private int [][] reduceArraySize(int a[][], int newSize){
		int tmp[][] = new int[newSize][];
		for (int i=0; i<newSize; i++){
			tmp[i] = new int[a[i].length];
			System.arraycopy(a[i], 0, tmp[i], 0, a[i].length);
		}
		return tmp;
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
	
	private void setForcefieldInputs(){
		// These values are specific to parm96a.dat
		//   parm96a.dat that I made that has Cl info
		switch(EnvironmentVars.forcefld){
			case AMBER:
				amberDatInFile = "parm96a.dat";
				break;
			case CHARMM22: 
				//KER: These numbers are specific to the charmm2Amber.dat file
				amberDatInFile = "parmcharmm22.dat";
				break;
			case CHARMM19:
			case CHARMM19NEUTRAL:
				//KER: These numbers are specific for charmm19
				amberDatInFile = "parmcharmm19.dat";
				break;
			default:
				System.out.println("DON'T RECOGNIZE FORCEFIELD: "+EnvironmentVars.forcefld.name());
				System.exit(0);
				break;
		}

	}
}
