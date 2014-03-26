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

///////////////////////////////////////////////////////////////////////////////////////////////
//SimpleMinimizer.java
//
//Version:           2.1 beta
//
//
//authors:
//    initials    name                 organization                email
//   ---------   -----------------    ------------------------    ----------------------------
//     RHL        Ryan Lilien          Dartmouth College           ryan.lilien@dartmouth.edu
//	   ISG		  Ivelin Georgiev	   Duke University			   ivelin.georgiev@duke.edu
//	  KER        Kyle E. Roberts       Duke University         ker17@duke.edu
//    PGC        Pablo Gainza C.       Duke University         pablo.gainza@duke.edu
///////////////////////////////////////////////////////////////////////////////////////////////

/* 
* Written by: Ryan Lilien (2001-2004);
* 		modified by Ivelin Georgiev (2004-2009)
*
* This class implements a simple energy minimization routine for the side-chains only
*  using the Amber96 forcefield, electrostatic, vdw, and dihedral energy
*  terms, and using the assumption that only certain torsions
*  can bend. Additionally there is a special residue, the
*  ligand, that can translate and globally rotate.
* 
* The system consists of one molecule containing: a protein strand, 
 *  a ligand strand (optional), and a cofactor strand (optional).
 * The protein strand does not have to contain sequential residues
 *  but it must be made of standard amino acids
 * The ligand strand can only be one 'thing'
 *  -if this 'thing' is an AA then the Penultimate Rotamer library is used
 *  -if this 'thing' is not an AA then a generic rotamer library is used
* 
*/



import java.io.Serializable;

/**
 * This class implements a simple energy minimization routine for the side-chains only. Handles the computation
 * of the Amber dihedral energy penalties (for minimizing away from the initial rotamer dihedrals).
 * Additionally there is a special residue, the ligand, that can translate and globally rotate.
 */
public class SimpleMinimizer implements Serializable {

	// If the debug flag is set to true then additional debug statements are
	//  printed to standard out.
	boolean debug = false;
	public int GLOBALNUM = 0; //Used to count how many times the minimizer
	  //has been called.
	
	//private int MAX_NUM_ATOMS_DISTAL = 30; // the maximum number of atoms distal to the dihedral; this is increased as needed for big ligands


	Molecule m = null;
	Amber96ext a96ff = null;
		// the molecule and forcefield (typical stuff)
	//int sysStrNum = -1;
	//int ligStrNum = -1;
		// the system and ligand strand numbers
		// if ligStrNum == -1 then there is no ligand
	StrandRotamers sysRH = null;
	StrandRotamers ligRH = null;
		// the system and ligand rotamer handlers, these
		//  are required so that we can get a handle on
		//  the dihedrals
	//int numSysDihedrals = 0;
	//int numLigDihedrals = 0;
		// the number of system and ligand dihedrals
	int strDihedralAtNums[][][] = null;
		// array of dihedrals, this array is n by 4
		//  and has the moleculeAtomNumbers for each
		//  of the 4 atoms of the dihedral
	int strDihedralDistal[][][] = null;
		// array of dihedrals, the first index is over all
		//  dihedrals, each row consists of an atomList
		//  that can be passed to the setTorsion
	int strNumAtomsDistal[][] = null;
		// the number of atoms distal for each dihedral,
		//  ie. the first x entries in a row of 
		//  sysDihedralDistal are valid atoms
	//int ligDihedralAtNums[][] = null;
		// array of dihedrals, this array is n by 4
		//  and has the moleculeAtomNumbers for each
		//  of the 4 atoms of the dihedral
	//int ligDihedralDistal[][] = null;
		// similar to sysDihedralArray but for the
		//  dihedrals of the ligand residue
	//int ligNumAtomsDistal[] = null;
		// the number of atoms distal for each dihedral,
		//  ie. the first x entries in a row of 
		//  ligDihedralDistal are valid atoms

	float initialAngleStepSize = 1.0f;
		// initial angular step size (in degrees)
	double strDihedDiff[][] = null;
	//double sysDihedDiff[] = null;
		// current finite difference step for the given dihedral
	//double ligDihedDiff[] = null;
		// current finite difference step for the given dihedral
	float tempCoords[] = null;
		// temporary holding location for coordinates as we
		//  muck with a dihedral, note that the size of this
		//  array is defined by MAX_NUM_ATOMS_DISTAL * 3
	double strCumulativeDihedStep[][] = null;
	//double sysCumulativeDihedStep[] = null;
		// the cumulative change in angle taken by the given
		//  dihedral. this is kept track of because we limit
		//  the total movement
	//double ligCumulativeDihedStep[] = null;
		// the cumulative change in angle taken by the given
		//  dihedral. this is kept track of because we limit
		//  the total movement
	double strStartCOM[][] = null;
		// the position of the ligand's initial COM
	double strCurCOM[][] = null;
		// the ligand's current COM

	double maxMovement = 9.0;
		// maximum degrees by which a torsion can
		//  cumulatively change
	float RotStep = 0.25f;
		// step size for the rigid rotation
		//  of the ligand
	float TransStep = 0.10f;
		// step size in � for the rigid ligand
		//  translation
	double MaxTrans = 1.2;
		// the maximum ligand translation allowed
	
	int numFlexRes[] = null;
		// the number of flexible residues
	int flexResAtomList[][] = new int[0][0];
		// each row is a residue and contains
		//  the moleucleatomnumbers located more
		//  distal than the 3rd atom of any dihedral
	int flexResListSize[] = new int[0];
		// the number of valid elements in each
		//  row of flexResAtomList
	int strDihedToResNum[][] = null;
		// mapping from dihedral number to residue
		//  number
	//int ligResNum = -1;
		// the residue index of the ligand in the
		//  flexResAtomList, flexResListSize
	int totalFlexRes = 0;
	int totalTransRotStrands = 0;
	StrandRotamers[] strandRot = null;
	int[] numStrDihedrals = null;
	
	// Section for dealing with dihedral energies
	// I know this shouldn't be in the minimizer but
	//  there's no other good way to do this.
	
	int strDihedPN[][] = new int[0][0];
	double strDihedWeight[][] = new double[0][0];
	int strDihedLocalNum[][] = new int[0][0];
	int strDihedNumTerms[] = new int[0];
	
	/*int sDihedPN[] = new int[0];
	double sDihedWeight[] = new double[0];
	int sDihedLocalNum[] = new int[0];
	int sDihedNumTerms = 0;*/
		// each of these arrays contains the PN and Weight for
		//  a dihedral term, each physical dihedral may contain
		//  multiple 'paths' (ie. set of four atoms that goes
		//  through the same central two). And each path may
		//  contain multiple terms (due to AMBER).
		// the local number is the number of the dihedral in
		//  the sysCumulativeDihedStep array
		// the dihednumterms is the total number of terms
	/*int lDihedPN[] = new int[0];
	double lDihedWeight[] = new double[0];
	int lDihedLocalNum[] = new int[0];
	int lDihedNumTerms = 0;*/
		// the same as the previous 4 variables, but for the
		//  ligand.
	boolean doDihedEnergy = false;
	// If true then dihedral energies are computed and
	//  used in minimization
	int numberOfStrands = 0;
	
	int MAX_NUM_ATOMS_DISTAL = 30;

	SimpleMinimizer(){
		// not much to do here
	
	}
	
		
	// Initializes the Minimizer for use with a system without a ligand
	public void initialize(Molecule theM, int numStrands, Amber96ext theA96ff,
		StrandRotamers strandRotamers[], int curAANum[], boolean doDihedral){
	
		
		
		
		// snag the local variables
		m = theM;
		a96ff = theA96ff;
		doDihedEnergy = doDihedral;
		//sysRH = theSysRH;
		//sysStrNum = theSysStrNum;
	
		strandRot = strandRotamers;
		
		//ligStrNum = -1;
		//numLigDihedrals = 0;
		//ligResNum = -1;
		/*for(int i=0;i<3;i++){
			ligCurCOM[i] = 0.0;
			ligStartCOM[i] = 0.0;
		}*/
		numberOfStrands = numStrands;
		numStrDihedrals = new int[numberOfStrands];
		
		for(int i=0;i<numberOfStrands;i++){
			for(int j=0;j<m.strand[i].numberOfResidues;j++){
				if(m.strand[i].residue[j].flexible)
					MAX_NUM_ATOMS_DISTAL = Math.max(MAX_NUM_ATOMS_DISTAL, m.strand[i].residue[j].numberOfAtoms);
			}
		}
		
		for(int i=0;i<numberOfStrands;i++){
			for(int j=0;j<m.strand[i].numberOfResidues;j++){
				if(m.strand[i].residue[j].flexible)
					numStrDihedrals[i] += strandRot[i].rl.getNumDihedrals(curAANum[m.strand[i].residue[j].moleculeResidueNumber]);
			}
		}		
			
		strDihedralAtNums = new int[numberOfStrands][][];
		strDihedralDistal = new int[numberOfStrands][][];
		strNumAtomsDistal = new int[numberOfStrands][];
		strDihedToResNum  = new int[numberOfStrands][];
		
		for(int i=0;i<numberOfStrands;i++){
			strDihedralAtNums[i] = new int[numStrDihedrals[i]][4];
			strDihedralDistal[i] = new int[numStrDihedrals[i]][MAX_NUM_ATOMS_DISTAL];
			strNumAtomsDistal[i] = new int[numStrDihedrals[i]];
			strDihedToResNum[i]  = new int[numStrDihedrals[i]];
		}		
		
		/*sysDihedralAtNums = new int[numSysDihedrals][4];
		sysDihedralDistal = new int[numSysDihedrals][MAX_NUM_ATOMS_DISTAL];
		sysNumAtomsDistal = new int[numSysDihedrals];*/
		tempCoords = new float[MAX_NUM_ATOMS_DISTAL * 3];

		// Count number of flexible residues
		numFlexRes = new int[numberOfStrands];
		for(int str=0; str<numberOfStrands;str++){
			for(int i=0;i<m.strand[str].numberOfResidues;i++){
				if(m.strand[str].residue[i].flexible)
					numFlexRes[str]++;
		}
		}
		
		// 2 is added to numFlexRes so that there is room for the ligand at the
		//  end if there is a ligand present, if there is no ligand then it
		//  doesn't really matter as we won't look at that row
		// The first ligand term includes nonbonded terms for computing dihedrals
		// The second ligand term includes nonbonded terms for computing the gradient
		//  and thus includes terms for all atoms
		totalFlexRes = 0;
		for(int i=0; i<numberOfStrands;i++)
			totalFlexRes += numFlexRes[i];
		
		totalTransRotStrands = 0;
		for(int i=0; i<numberOfStrands;i++)
			if(m.strand[i].rotTrans)
				totalTransRotStrands++;
				
		flexResAtomList = new int[totalFlexRes+totalTransRotStrands][MAX_NUM_ATOMS_DISTAL];
		flexResListSize = new int[totalFlexRes+totalTransRotStrands];
		
		// Now determine which atoms are involved with each system dihedral
		int curTransRotInd = -1;
		int strResNumPP = -1;
		int curNumFlex = 0;
		
		for(int str=0; str<numberOfStrands;str++){
		int atoms[] = new int[4];
		int tmpAtomList[] = new int[MAX_NUM_ATOMS_DISTAL];
		int at3num = -1, at2num = -1;
		int numDihed = 0;
			
			if(m.strand[str].rotTrans){
				curTransRotInd++;
				strResNumPP = totalFlexRes+curTransRotInd;
				flexResListSize[strResNumPP] = theM.strand[str].numberOfAtoms;
				flexResAtomList[strResNumPP] = new int[flexResListSize[strResNumPP]];
			}
			Residue localRes = null;
			StrandRotamers localRH = null;
			int prevNumAtoms = 0;
			for(int i=0;i<m.strand[str].numberOfResidues;i++){
				localRes = m.strand[str].residue[i];
				if(m.strand[str].rotTrans){
					for(int k=0;k<localRes.numberOfAtoms;k++){
						flexResAtomList[strResNumPP][k+prevNumAtoms] = localRes.atom[k].moleculeAtomNumber;
					}
				}
				prevNumAtoms += localRes.numberOfAtoms;
				if(localRes.flexible){
					
					for(int j=0;j<strandRot[str].rl.getNumDihedrals(curAANum[theM.strand[str].residue[i].moleculeResidueNumber]);j++){
						strDihedToResNum[str][numDihed] = curNumFlex;
						atoms = strandRot[str].rl.getDihedralInfo(m,str,i,j);
					// note: atoms are residue relative numbering
						strDihedralAtNums[str][numDihed][0] = localRes.atom[atoms[0]].moleculeAtomNumber;
						strDihedralAtNums[str][numDihed][1] = localRes.atom[atoms[1]].moleculeAtomNumber;
						strDihedralAtNums[str][numDihed][2] = localRes.atom[atoms[2]].moleculeAtomNumber;
						strDihedralAtNums[str][numDihed][3] = localRes.atom[atoms[3]].moleculeAtomNumber;
					
					at2num = localRes.atom[atoms[2]].moleculeAtomNumber;
					at3num = localRes.atom[atoms[3]].moleculeAtomNumber;					
					for(int k=0;k<MAX_NUM_ATOMS_DISTAL;k++){
						tmpAtomList[k]=1;
					}
					tmpAtomList[atoms[1]]=0;
					tmpAtomList[atoms[2]]=0;
						strandRot[str].getAtomsMoreDistal(m,localRes.moleculeResidueNumber,m.atom[at2num],tmpAtomList);
						strNumAtomsDistal[str][numDihed]=0;
					for(int k=0;k<MAX_NUM_ATOMS_DISTAL;k++){
						if ((tmpAtomList[k]==2) && (k != atoms[3])){
								strDihedralDistal[str][numDihed][strNumAtomsDistal[str][numDihed]]=localRes.atom[k].moleculeAtomNumber;
								strNumAtomsDistal[str][numDihed] += 1;
						}
					}
					if (j==0){
						// If this is the first dihedral for this residue snag information
						flexResListSize[curNumFlex] = strNumAtomsDistal[str][numDihed]+1;
						flexResAtomList[curNumFlex][0] = at3num;
							for(int k=0;k<strNumAtomsDistal[str][numDihed];k++){
								flexResAtomList[curNumFlex][k+1] = strDihedralDistal[str][numDihed][k];
						}
					}
					numDihed++;
				}
				curNumFlex++;
				}
			}
		}
	}	
	
	// Initializes the Minimizer for use with a system and a ligand
	/*public void initialize(Molecule theM, int theSysStrNum, int theLigStrNum,
		Amber96ext theA96ff, StrandRotamers theSysRH, StrandRotamers theLigRH, int
		curAANum[], int curLigAANum, boolean doDihedral, RotamerLibrary rl, RotamerLibrary grl){
		
		MAX_NUM_ATOMS_DISTAL = Math.max(MAX_NUM_ATOMS_DISTAL, theM.strand[theLigStrNum].residue[0].numberOfAtoms);
	
		// First call the other initialize interface to grab the
		//  system dihedrals
		initialize(theM,theSysStrNum,theA96ff,theSysRH,curAANum,doDihedral,rl);
		doDihedEnergy = doDihedral;
		
		ligRH = theLigRH;
		ligStrNum = theLigStrNum;

		// Now compute the ligand dihedrals
		// First determine how many there are
		numLigDihedrals = grl.getNumDihedrals(curLigAANum);
		
		ligDihedralAtNums = new int[numLigDihedrals][4];
		ligDihedralDistal = new int[numLigDihedrals][MAX_NUM_ATOMS_DISTAL];
		ligNumAtomsDistal = new int[numLigDihedrals];

		// numFlexRes, flexResAtomList, flexResListSize are setup in the system
		//  initialize code. The ligand information at numFlexRes are the dihedral
		//  relevant nonbonded terms, the ligand information at numFlexRes+1 are
		//  used in computing the gradient and include all atoms of the ligand not
		//  just those distal to the first dihedral
	
		// Now determine which atoms are involved with each dihedral
		int atoms[] = new int[4];
		int tmpAtomList[] = new int[MAX_NUM_ATOMS_DISTAL];
		int at3num = -1, at2num = -1;
		int numDihed = 0;
		Residue localRes = m.strand[ligStrNum].residue[0];
		// Snag all ligand atom information for flexResAtomList, flexResListSize
		ligResNum = numFlexRes;  // numFlexRes was set in the other initialize call
		int ligResNumPP = ligResNum + 1;
		flexResListSize[ligResNumPP] = localRes.numberOfAtoms;
		for(int k=0;k<flexResListSize[ligResNumPP];k++){
			flexResAtomList[ligResNumPP][k] = localRes.atom[k].moleculeAtomNumber;
		}

		if(m.strand[ligStrNum].residue[0].flexible){
			for(int j=0;j<grl.getNumDihedrals(curLigAANum);j++){
				// ligDihedToResNum[numDihed] = 0; // all ligand dihedrals are from residue 0
				atoms = grl.getDihedralInfo(m,ligStrNum,0,j);
				// note: atoms are residue relative numbering
				ligDihedralAtNums[numDihed][0] = localRes.atom[atoms[0]].moleculeAtomNumber;
				ligDihedralAtNums[numDihed][1] = localRes.atom[atoms[1]].moleculeAtomNumber;
				ligDihedralAtNums[numDihed][2] = localRes.atom[atoms[2]].moleculeAtomNumber;
				ligDihedralAtNums[numDihed][3] = localRes.atom[atoms[3]].moleculeAtomNumber;
				
				at2num = localRes.atom[atoms[2]].moleculeAtomNumber;
				at3num = localRes.atom[atoms[3]].moleculeAtomNumber;
					
				for(int k=0;k<MAX_NUM_ATOMS_DISTAL;k++){
					tmpAtomList[k]=1;
				}
				tmpAtomList[atoms[1]]=0;
				tmpAtomList[atoms[2]]=0;
				ligRH.getAtomsMoreDistal(m,localRes.moleculeResidueNumber,m.atom[at2num],tmpAtomList);
				ligNumAtomsDistal[numDihed]=0;
				for(int k=0;k<MAX_NUM_ATOMS_DISTAL;k++){
					if ((tmpAtomList[k]==2) && (k != atoms[3])){
						ligDihedralDistal[numDihed][ligNumAtomsDistal[numDihed]]=localRes.atom[k].moleculeAtomNumber;
						ligNumAtomsDistal[numDihed] += 1;
					}
				}
				if (j==0){
					// If this is the first dihedral for this residue snag information
					flexResListSize[ligResNum] = ligNumAtomsDistal[numDihed]+1;
					flexResAtomList[ligResNum][0] = at3num;
					for(int k=0;k<ligNumAtomsDistal[numDihed];k++){
						flexResAtomList[ligResNum][k+1] = ligDihedralDistal[numDihed][k];
					}
				}
				numDihed++;
			}
		}
	}*/
	
	public int getNumTotDihed(){
		int numDihedrals = 0;
		for(int str=0;str<numberOfStrands;str++)
			numDihedrals+=numStrDihedrals[str];
		
		return numDihedrals;
	}
	
	// Sets the initial angle step size
	public void setInitialAngleStepSize(float num){	
		initialAngleStepSize = num;
	}
	
	// Setup the dihedral angle terms for both the system and ligand
	public boolean setupDihedralTerms(){
		
		// At this point we assume sysDihedralAtNums, ligDihedralAtNums,
		//  numSysDihedrals, and numLigDihedrals are all correct.
		
		int tmpStrSize = 1000;
		//int tmpLSize = 1000;
		strDihedPN = new int[numberOfStrands][tmpStrSize];
		//sDihedPN = new int[tmpSSize];
		//lDihedPN = new int[tmpLSize];
		strDihedWeight = new double[numberOfStrands][tmpStrSize];
		//sDihedWeight = new double[tmpSSize];
		//lDihedWeight = new double[tmpLSize];
		strDihedLocalNum = new int[numberOfStrands][tmpStrSize];
		//sDihedLocalNum = new int[tmpSSize];
		//lDihedLocalNum = new int[tmpLSize];
		strDihedNumTerms = new int[numberOfStrands];
		for(int i=0; i<numberOfStrands;i++)
			strDihedNumTerms[i] = 0;
		//sDihedNumTerms = 0;
		//lDihedNumTerms = 0;
		int atom2 = 0, atom3 = 0;
		double fConst[] = new double[5];
		double eqAngle[] = new double[5];
		int terms[] = new int[5];
		int multiplic[] = new int[1];
		
		int newDihedPN[] = null;
		double newDihedWeight[] = null;
		int newDihedLocalNum[] = null;
		for(int str=0;str<numberOfStrands;str++){
		// SYSTEM DIHEDRALS
			for(int i=0;i<numStrDihedrals[str];i++) {
			// Get center two atoms
				atom2 = strDihedralAtNums[str][i][1];
				atom3 = strDihedralAtNums[str][i][2];
			for(int q=0;q<m.connected[atom2][0];q++) {
				if (m.connected[atom2][q+1] != atom3) {
					for(int w=0;w<m.connected[atom3][0];w++) {
						if (m.connected[atom3][w+1] != atom2) {
							// At this point 'q',atom2,atom3,'w' is a dihedral
							if (!(a96ff.getTorsionParameters(m.atom[m.connected[atom2][q+1]].type,m.atom[atom2].type,
															 m.atom[atom3].type,m.atom[m.connected[atom3][w+1]].type,
															 fConst,eqAngle,terms,multiplic))) {
								System.out.println("WARNING: Could not find torsion parameters for " + m.connected[atom2][q+1] + "," + atom2 + "," + atom3 +"," + m.connected[atom3][w+1]);
								return false;
							}
							for(int r=0;r<=multiplic[0];r++) {
									strDihedNumTerms[str]++;
									if (strDihedNumTerms[str] > tmpStrSize) {
									// increase the size of the arrays
										newDihedPN = new int[tmpStrSize + 500];
										newDihedWeight = new double[tmpStrSize + 500];
										newDihedLocalNum = new int[tmpStrSize + 500];
										System.arraycopy(strDihedPN[str],0,newDihedPN,0,tmpStrSize);
										System.arraycopy(strDihedWeight[str],0,newDihedWeight,0,tmpStrSize);
										System.arraycopy(strDihedLocalNum[str],0,newDihedLocalNum,0,tmpStrSize);
										strDihedPN[str] = newDihedPN;
										strDihedWeight[str] = newDihedWeight;
										strDihedLocalNum[str] = newDihedLocalNum;
										tmpStrSize += 500;
								}
									strDihedPN[str][strDihedNumTerms[str]] = terms[r];
									strDihedWeight[str][strDihedNumTerms[str]] = fConst[r];
									strDihedLocalNum[str][strDihedNumTerms[str]] = i;
							}
						}
					}
				}
			}
		}
		
		// Shrink the size of the array
			newDihedPN = new int[strDihedNumTerms[str]];
			newDihedWeight = new double[strDihedNumTerms[str]];
			newDihedLocalNum = new int[strDihedNumTerms[str]];
			System.arraycopy(strDihedPN[str],0,newDihedPN,0,strDihedNumTerms[str]);
			System.arraycopy(strDihedWeight[str],0,newDihedWeight,0,strDihedNumTerms[str]);
			System.arraycopy(strDihedLocalNum[str],0,newDihedLocalNum,0,strDihedNumTerms[str]);
			strDihedPN[str] = newDihedPN;
			strDihedWeight[str] = newDihedWeight;
			strDihedLocalNum[str] = newDihedLocalNum;
		}
		/*// LIGAND DIHEDRALS
		for(int i=0;i<numLigDihedrals;i++) {
			// Get center two atoms
			atom2 = ligDihedralAtNums[i][1];
			atom3 = ligDihedralAtNums[i][2];
			for(int q=0;q<m.connected[atom2][0];q++) {
				if (m.connected[atom2][q+1] != atom3) {
					for(int w=0;w<m.connected[atom3][0];w++) {
						if (m.connected[atom3][w+1] != atom2) {
							// At this point 'q',atom2,atom3,'w' is a dihedral
							if (!(a96ff.getTorsionParameters(m.atom[m.connected[atom2][q+1]].type,m.atom[atom2].type,
															 m.atom[atom3].type,m.atom[m.connected[atom3][w+1]].type,
															 fConst,eqAngle,terms,multiplic))) {
								System.out.println("WARNING: Could not find torsion parameters for " + m.connected[atom2][q+1] + "," + atom2 + "," + atom3 +"," + m.connected[atom3][w+1]);
								return false;
							}
							for(int r=0;r<=multiplic[0];r++) {
								lDihedNumTerms++;
								if (lDihedNumTerms > tmpLSize) {
									// increase the size of the arrays
									newDihedPN = new int[tmpLSize + 500];
									newDihedWeight = new double[tmpLSize + 500];
									newDihedLocalNum = new int[tmpLSize + 500];
									System.arraycopy(lDihedPN,0,newDihedPN,0,tmpLSize);
									System.arraycopy(lDihedWeight,0,newDihedWeight,0,tmpLSize);
									System.arraycopy(lDihedLocalNum,0,newDihedLocalNum,0,tmpLSize);
									lDihedPN = newDihedPN;
									lDihedWeight = newDihedWeight;
									lDihedLocalNum = newDihedLocalNum;
									tmpLSize += 500;
								}
								lDihedPN[lDihedNumTerms] = terms[r];
								lDihedWeight[lDihedNumTerms] = fConst[r];
								lDihedLocalNum[lDihedNumTerms] = i;
							}
						}
					}
				}
			}
		}
		
		// Shrink the size of the array
		newDihedPN = new int[lDihedNumTerms];
		newDihedWeight = new double[lDihedNumTerms];
		newDihedLocalNum = new int[lDihedNumTerms];
		System.arraycopy(lDihedPN,0,newDihedPN,0,lDihedNumTerms);
		System.arraycopy(lDihedWeight,0,newDihedWeight,0,lDihedNumTerms);
		System.arraycopy(lDihedLocalNum,0,newDihedLocalNum,0,lDihedNumTerms);
		lDihedPN = newDihedPN;
		lDihedWeight = newDihedWeight;
		lDihedLocalNum = newDihedLocalNum;
		*/

		return true;
	}
	
	
	// Computes and returns the dihedral energy for the dihedrals
	//  that are allowed to change during minimization (both for
	//  the system and the ligand). Each dihedral is assumed to
	//  start in a energy well of the appropriate dihedral energy
	//  terms from the AMBER forcefield.
	// Should be called only after calling minimize()
	public double computeDihedEnergy() {
	
		double dihedE = 0.0;
		double d2r = 3.14159265 / 180.0;

		// System dihedral terms first
		for(int str=0; str<numberOfStrands;str++){
			for(int i=0;i<strDihedNumTerms[str];i++){
				dihedE += strDihedWeight[str][i] * (1 - Math.cos(strDihedPN[str][i]*strCumulativeDihedStep[str][strDihedLocalNum[str][i]]*d2r));
			}
		}
		
		/*for(int i=0;i<sDihedNumTerms;i++){
			dihedE += sDihedWeight[i] * (1 - Math.cos(sDihedPN[i]*sysCumulativeDihedStep[sDihedLocalNum[i]]*d2r));
		}
		// Ligand dihedral terms
		for(int i=0;i<lDihedNumTerms;i++){
			dihedE += lDihedWeight[i] * (1 - Math.cos(lDihedPN[i]*ligCumulativeDihedStep[lDihedLocalNum[i]]*d2r));
		}*/
		
		//System.out.println("computeDihedEnergy "+dihedE);
		
		return(dihedE);
	}
	
	// This function sets up the dihedral energy computation for the
	//	given dihedral; the computation is done in computeOneDihedEnergyDiffHelper()
	public double computeOneDihedEnergyDiff(int strNumber, int dihedNumForCur, double dihedChange){
		int[] DNum = new int[numberOfStrands];
		double[] DChange = new double[numberOfStrands];
		
		for(int str=0;str<numberOfStrands;str++){
			if(str == strNumber){
				DNum[str] = dihedNumForCur;
				DChange[str] = dihedChange;
			}
			else{
				DNum[str] = -1;
				DChange[str] = 0.0;
			}
		}
		
		return computeOneDihedEnergyDiffHelper(DNum,DChange);
	}
	
	/*public double computeOneDihedEnergyDiff(boolean isLigDihed, int dihedNumForCur, double dihedChange){
		
		int sDNum, lDNum;
		double sDChange, lDChange;
		if (!isLigDihed) { //system dihedral
			sDNum = dihedNumForCur;
			sDChange = dihedChange;
			lDNum = -1;
			lDChange = 0.0;
		}
		else { //ligand dihedral
			sDNum = -1;
			sDChange = 0.0;
			lDNum = dihedNumForCur;
			lDChange = dihedChange;
		}
		
		return computeOneDihedEnergyDiffHelper(sDNum, sDChange, lDNum, lDChange);
	}	*/

	// Computes the difference in energy for specified dihedral
	//  initial - change. 
	// System dihedral sDNum is changed by sDChange and ligand
	//  dihedral lDNum is changed by lDChange. These changes are
	//  only temporary and used to compute an energy, that is the
	//  sys(lig)CumulativeDihedStep is not changed. If sDNum or
	//  lDNum are -1 then they are not changed.
	public double computeOneDihedEnergyDiffHelper(int[] DNum, double[] DChange){
		double dihedE = 0.0;
		double d2r = 3.14159265 / 180.0;
		
		for(int str=0; str<DNum.length;str++){
			if (DNum[str] >= 0) {
				for(int i=0;i<strDihedNumTerms[str];i++) {
					if (strDihedLocalNum[str][i] == DNum[str])
						dihedE -= strDihedWeight[str][i] * (1 - Math.cos(strDihedPN[str][i]*(DChange[str])*d2r));
				}
			}
		}
		
		return dihedE;
	}
	
	/*public double computeOneDihedEnergyDiffHelper(int sDNum, double sDChange,
										int lDNum, double lDChange) {
		
		double dihedE = 0.0;
		double d2r = 3.14159265 / 180.0;
		
		// System dihedral terms first
		if (sDNum >= 0) {
			for(int i=0;i<sDihedNumTerms;i++) {
				if (sDihedLocalNum[i] == sDNum)
					dihedE -= sDihedWeight[i] * (1 - Math.cos(sDihedPN[i]*(sDChange)*d2r));
			}
		}
		// Ligand dihedral terms
		if (lDNum >= 0) {
			for(int i=0;i<lDihedNumTerms;i++) {
				if (lDihedLocalNum[i] == lDNum)
					dihedE -= lDihedWeight[i] * (1 - Math.cos(lDihedPN[i]*(lDChange)*d2r));
			}
		}

		return(dihedE);
	}*/

	// Returns the cross product: a x b
	public double[] crossProduct(double a[], double b[]){

		double f[] = new double[3];
		f[0] = a[1]*b[2] - a[2]*b[1];
		f[1] = a[2]*b[0] - a[0]*b[2];
		f[2] = a[0]*b[1] - a[1]*b[0];
		return(f);
	}
	
	// Returns the dot product of a and b
	public double dotProduct(double a[], double b[]){
		return(a[0]*b[0] + a[1]*b[1] + a[2]*b[2]);
	}

	public double getMagnitude(double a[]){
		double sum = 0.0;
		for(int i=0;i<a.length;i++)
			sum += (a[i]*a[i]);
		return(Math.sqrt(sum));
	}
	
	// Computes the torque around the specified dihedral
	public double[] computeTorque(int diArray[]){
	
		// remember that the indices of the dihedralArray are
		// [0] residueNum
		// [1] numberOfAtoms (n)
		// [2] atom1
		// [3] atom2
		// [4] atom3
		// [5] atom4 (the distal atom of the dihedral)
		// ...
		// [n+1] atomn
		// the result is a float array where the first
		//  three values are the axis and the fourth
		//  is the torque
		// interestingly since our force is kcal/mol/A
		//  the torque is kcal*A/mol/A or kcal/mol
		
		double result[] = new double[4];
		double torque[] = new double[3];
		double temp[] = new double[3];
		double R[] = new double[3];
		double F[] = new double[3];
		int atBase = -1;
		Residue localRes = m.residue[diArray[0]];
		int at3Base = localRes.atom[diArray[4]].moleculeAtomNumber * 3;
		int at2Base = localRes.atom[diArray[3]].moleculeAtomNumber * 3;
		
		// compute the rotation (dihedral) axis
		result[0] = m.actualCoordinates[at3Base] - m.actualCoordinates[at2Base];
		result[1] = m.actualCoordinates[at3Base+1] - m.actualCoordinates[at2Base+1];
		result[2] = m.actualCoordinates[at3Base+2] - m.actualCoordinates[at2Base+2];
		// normalize the axis to unit length
		double axisMag = getMagnitude(result);
		result[0] /= axisMag;
		result[1] /= axisMag;
		result[2] /= axisMag;
		
		torque[0] = 0.0;
		torque[1] = 0.0;
		torque[2] = 0.0;

		// loop over all involved atoms
		for(int i=5;i<diArray[1]+2;i++){
			atBase = localRes.atom[diArray[i]].moleculeAtomNumber * 3;
			// R is the vector from atom3 to the atom of interest
			R[0] = m.actualCoordinates[atBase] - m.actualCoordinates[at3Base];
			R[1] = m.actualCoordinates[atBase+1] - m.actualCoordinates[at3Base+1];
			R[2] = m.actualCoordinates[atBase+2] - m.actualCoordinates[at3Base+2];
			// F is the current gradient on the atom of interest
			F[0] = m.gradient[atBase];
			F[1] = m.gradient[atBase+1];
			F[2] = m.gradient[atBase+2];
			temp = crossProduct(R,F);
			torque[0] += temp[0];
			torque[1] += temp[1];
			torque[2] += temp[2];
		}
		// Project the sum of the torques onto the axis of rotation
		//  which is stored in the first 3 elements of result
		result[3] = dotProduct(torque,result);

		return(result);
	}
	
	//Sets up and performs ligand translation/rotation
	/*private void doLigTransRot(double ligTorque[], double ligTrans[], double lligRotStep, double lligTransStep, double lligMaxTrans){
		
		float bckpLigCoords[] = backupLigCoord(); //backup actual ligand coordinates
		double initE[] = a96ff.calculateTotalEnergy(m.actualCoordinates, -1); //compute energy before translation/rotation
	
		a96ff.calculateGradient(-1); //calculate the gradient (perhaps eventually just calculate part of the gradient)
		//a96ff.calculateEVGradientPartWithArrays(ligResNumPP);
		computeLigTorqueTrans(ligTorque, ligTrans);
		applyLigTorqueTrans(ligTorque, lligRotStep, ligTrans, lligTransStep, lligMaxTrans);
		
		double minE[] = a96ff.calculateTotalEnergy(m.actualCoordinates, -1); //compute energy after translation/rotation
		if (initE[0]<minE[0]) { //restore initial ligand actual coordinates if trans/rot increses energy
			restoreLigCoord(bckpLigCoords);
		}
	}*/

	protected void doStrTransRot(int strNumber, double strTorque[], double strTrans[], double RotStep, double TransStep, double MaxTrans){
		
		float bckpStrCoords[] = backupStrCoord(strNumber); //backup actual ligand coordinates
		double initE[] = a96ff.calculateTotalEnergy(m.actualCoordinates, -1); //compute energy before translation/rotation
	
		a96ff.calculateGradient(-1); //calculate the gradient (perhaps eventually just calculate part of the gradient)
		//a96ff.calculateEVGradientPartWithArrays(ligResNumPP);
		computeStrTorqueTrans(strNumber, strTorque, strTrans);
		applyStrTorqueTrans(strNumber, strTorque, RotStep, strTrans, TransStep, MaxTrans);
		
		double minE[] = a96ff.calculateTotalEnergy(m.actualCoordinates, -1); //compute energy after translation/rotation
		//m.saveMolecule("minSingle.pdb", 0.0f);
		if (initE[0]<=minE[0]) { //restore initial ligand actual coordinates if trans/rot increses energy
			restoreStrCoord(strNumber, bckpStrCoords);	
		}
	}

	// Computes the torque around the ligand
	public void computeStrTorqueTrans(int strNumber, double strTorque[], double strTrans[]){
	
		Strand thisStrand = m.strand[strNumber];
		double centOfMass[] = m.getStrandCOM(strNumber);
		int atBase = -1;
		double R[] = new double[3];
		double F[] = new double[3];
		double temp[];
		strTorque[0] = 0.0;
		strTorque[1] = 0.0;
		strTorque[2] = 0.0;
		strTrans[0] = 0.0;
		strTrans[1] = 0.0;
		strTrans[2] = 0.0;

		// loop over all involved atoms
		for(int i=0;i<thisStrand.numberOfResidues;i++){
			for(int j=0;j<thisStrand.residue[i].numberOfAtoms;j++){
				atBase = thisStrand.residue[i].atom[j].moleculeAtomNumber * 3;
				// R is the vector from the COM to the current atom
				R[0] = m.actualCoordinates[atBase] - centOfMass[0];
				R[1] = m.actualCoordinates[atBase+1] - centOfMass[1];
				R[2] = m.actualCoordinates[atBase+2] - centOfMass[2];
				// F is the current gradient for the atom of interest
				F[0] = m.gradient[atBase];
				F[1] = m.gradient[atBase+1];
				F[2] = m.gradient[atBase+2];
				temp = crossProduct(R,F);
				strTorque[0] += temp[0];
				strTorque[1] += temp[1];
				strTorque[2] += temp[2];
				strTrans[0] += F[0];
				strTrans[1] += F[1];
				strTrans[2] += F[2];
			}
		}

	}
	
	/*// Computes the torque around the ligand
	public void computeLigTorqueTrans(double ligTorque[], double ligTrans[]){
	
		Strand thisStrand = m.strand[ligStrNum];
		double centOfMass[] = m.getStrandCOM(ligStrNum);
		int atBase = -1;
		double R[] = new double[3];
		double F[] = new double[3];
		double temp[];
		ligTorque[0] = 0.0;
		ligTorque[1] = 0.0;
		ligTorque[2] = 0.0;
		ligTrans[0] = 0.0;
		ligTrans[1] = 0.0;
		ligTrans[2] = 0.0;

		// loop over all involved atoms
		for(int i=0;i<thisStrand.numberOfResidues;i++){
			for(int j=0;j<thisStrand.residue[i].numberOfAtoms;j++){
				atBase = thisStrand.residue[i].atom[j].moleculeAtomNumber * 3;
				// R is the vector from the COM to the current atom
				R[0] = m.actualCoordinates[atBase] - centOfMass[0];
				R[1] = m.actualCoordinates[atBase+1] - centOfMass[1];
				R[2] = m.actualCoordinates[atBase+2] - centOfMass[2];
				// F is the current gradient for the atom of interest
				F[0] = m.gradient[atBase];
				F[1] = m.gradient[atBase+1];
				F[2] = m.gradient[atBase+2];
				temp = crossProduct(R,F);
				ligTorque[0] += temp[0];
				ligTorque[1] += temp[1];
				ligTorque[2] += temp[2];
				ligTrans[0] += F[0];
				ligTrans[1] += F[1];
				ligTrans[2] += F[2];
			}
		}

	}*/

	// Applies the specified torque and rigid body translation to the
	//  atoms in the ligand
	// The total translation can not be more than lligMaxTrans
	// The torque is applied with a magnitude of lligRotStep
	private void applyStrTorqueTrans(int strNumber, double strTorque[], double RotStep,
		double strTrans[], double TransStep, double MaxTrans){
	
		// normalize the ligand translation
		// ligTrans will now have magnitude lligTransStep
		// we change sign because the gradient points uphill
		double scale = getMagnitude(strTrans)/TransStep;
		if (scale != 0.0) {
			strTrans[0] /= -scale;
			strTrans[1] /= -scale;
			strTrans[2] /= -scale;
	}

		// determine the new COM if we took this step
		double tempCOM[] = new double[3];
		tempCOM[0] = strCurCOM[strNumber][0] + strTrans[0];
		tempCOM[1] = strCurCOM[strNumber][1] + strTrans[1];
		tempCOM[2] = strCurCOM[strNumber][2] + strTrans[2];

		// compute how large of a step this would be from the
		//  start point
		double totalMovement [] = new double[3];
		totalMovement[0] = tempCOM[0] - strStartCOM[strNumber][0];
		totalMovement[1] = tempCOM[1] - strStartCOM[strNumber][1];
		totalMovement[2] = tempCOM[2] - strStartCOM[strNumber][2];
		scale = getMagnitude(totalMovement);
		
		// if the step would take us too far away then
		//  scale it back
		if (scale > MaxTrans){
			scale = scale / MaxTrans;
			totalMovement[0] /= scale;
			totalMovement[1] /= scale;
			totalMovement[2] /= scale;
		}
		
		// compute the translation to get us to the new COM
		float theTranslation[] = new float[3];
		theTranslation[0] = (float)(strStartCOM[strNumber][0] + totalMovement[0] - strCurCOM[strNumber][0]);
		theTranslation[1] = (float)(strStartCOM[strNumber][1] + totalMovement[1] - strCurCOM[strNumber][1]);
		theTranslation[2] = (float)(strStartCOM[strNumber][2] + totalMovement[2] - strCurCOM[strNumber][2]);

		// update the current COM
		strCurCOM[strNumber][0] = strStartCOM[strNumber][0]+totalMovement[0];
		strCurCOM[strNumber][1] = strStartCOM[strNumber][1]+totalMovement[1];
		strCurCOM[strNumber][2] = strStartCOM[strNumber][2]+totalMovement[2];
		
		// Apply torque
		// Since our step size is fixed we just need to know the
		//  vector to rotate around
		// we change sign because the gradient points uphill
		double mag = getMagnitude(strTorque);
		if (mag != 0.0) {
			strTorque[0] /= -mag;
			strTorque[1] /= -mag;
			strTorque[2] /= -mag;
			// The last parameter says not to update the atom coordintes,
			//  only update the actualCoords
			if (Math.abs(RotStep)>0.0001)   // this should be in
				m.rotateStrandAroundCOM(strNumber,(float)strTorque[0],(float)strTorque[1],(float)strTorque[2],(float)RotStep,false);
		}
		// Apply translation, don't update atom coordinates
		// The last parameter says not to update the atom coordintes,
		//  only update the actualCoords
		m.translateStrand(strNumber,theTranslation[0],theTranslation[1],theTranslation[2],false);
		
	}
	
	// Applies the specified torque and rigid body translation to the
	//  atoms in the ligand
	// The total translation can not be more than lligMaxTrans
	// The torque is applied with a magnitude of lligRotStep
	/*private void applyLigTorqueTrans(double ligTorque[], double lligRotStep,
		double ligTrans[], double lligTransStep, double lligMaxTrans){
	
		// normalize the ligand translation
		// ligTrans will now have magnitude lligTransStep
		// we change sign because the gradient points uphill
		double scale = getMagnitude(ligTrans)/lligTransStep;
		if (scale != 0.0) {
			ligTrans[0] /= -scale;
			ligTrans[1] /= -scale;
			ligTrans[2] /= -scale;
		}
		
		// determine the new COM if we took this step
		double tempCOM[] = new double[3];
		tempCOM[0] = ligCurCOM[0] + ligTrans[0];
		tempCOM[1] = ligCurCOM[1] + ligTrans[1];
		tempCOM[2] = ligCurCOM[2] + ligTrans[2];
		
		// compute how large of a step this would be from the
		//  start point
		double totalMovement [] = new double[3];
		totalMovement[0] = tempCOM[0] - ligStartCOM[0];
		totalMovement[1] = tempCOM[1] - ligStartCOM[1];
		totalMovement[2] = tempCOM[2] - ligStartCOM[2];
		scale = getMagnitude(totalMovement);
		
		// if the step would take us too far away then
		//  scale it back
		if (scale > lligMaxTrans){
			scale = scale / lligMaxTrans;
			totalMovement[0] /= scale;
			totalMovement[1] /= scale;
			totalMovement[2] /= scale;
		}
		
		// compute the translation to get us to the new COM
		float theTranslation[] = new float[3];
		theTranslation[0] = (float)(ligStartCOM[0] + totalMovement[0] - ligCurCOM[0]);
		theTranslation[1] = (float)(ligStartCOM[1] + totalMovement[1] - ligCurCOM[1]);
		theTranslation[2] = (float)(ligStartCOM[2] + totalMovement[2] - ligCurCOM[2]);

		// update the current COM
		ligCurCOM[0] = ligStartCOM[0]+totalMovement[0];
		ligCurCOM[1] = ligStartCOM[1]+totalMovement[1];
		ligCurCOM[2] = ligStartCOM[2]+totalMovement[2];
		
		// Apply torque
		// Since our step size is fixed we just need to know the
		//  vector to rotate around
		// we change sign because the gradient points uphill
		double mag = getMagnitude(ligTorque);
		if (mag != 0.0) {
			ligTorque[0] /= -mag;
			ligTorque[1] /= -mag;
			ligTorque[2] /= -mag;
			// The last parameter says not to update the atom coordintes,
			//  only update the actualCoords
			if (Math.abs(lligRotStep)>0.0001)   // this should be in
				m.rotateStrandAroundCOM(ligStrNum,(float)ligTorque[0],(float)ligTorque[1],(float)ligTorque[2],(float)lligRotStep,false);
		}
		// Apply translation, don't update atom coordinates
		// The last parameter says not to update the atom coordintes,
		//  only update the actualCoords
		m.translateStrand(ligStrNum,theTranslation[0],theTranslation[1],theTranslation[2],false);
		
	}*/

	private void storeCoord(int index, int mAtNum){
	
		tempCoords[index*3] = m.actualCoordinates[mAtNum*3];
		tempCoords[index*3 + 1] = m.actualCoordinates[mAtNum*3 + 1];
		tempCoords[index*3 + 2] = m.actualCoordinates[mAtNum*3 + 2];
	}

	private void restoreCoord(int index, int mAtNum){
		
		m.actualCoordinates[mAtNum*3] = tempCoords[index*3];
		m.actualCoordinates[mAtNum*3 + 1] = tempCoords[index*3 + 1];
		m.actualCoordinates[mAtNum*3 + 2] = tempCoords[index*3 + 2];
	}
	
	protected void applyDihedStep(int diAtArray[], int atomList[], int alSize,
		double dihedDiff){
		
		// Perform the rotation
		// at1, at2, and at3 don't move, at4 and the atoms
		//  in the atomList rotate
		// The last parameter says not to update the atom coordintes,
		//  only update the actualCoords
		if (Math.abs(dihedDiff)>0.0001)
			m.changeTorsion(diAtArray[0],diAtArray[1],diAtArray[2],
				diAtArray[3],dihedDiff,atomList,alSize,false);
	}	
	
////////////////////////////////////////////////////////////////////////////////
//  Minimization Section
////////////////////////////////////////////////////////////////////////////////

	// For the given dihedral, compute initial energy, save coordinates,
	//  modify dihedral, compute new energy, compute difference in
	//  energy, restore coords, return step size and direction of lower energy
	protected float computeDihedDiff(int diAtArray[], int atomList[],
		int alSize, int AANum, float stepSize, int strNumber, int dihedNumForCur){
		
		double initialEnergy[], secondEnergy[];
		
		double cumulStep = 0.0;
		for(int str=0; str<numberOfStrands;str++){
			if (str==strNumber){
				cumulStep = strCumulativeDihedStep[str][dihedNumForCur];
				break;
			}
		}
		
		// Store coordinates
		storeCoord(0,diAtArray[3]);
		for(int i=0;i<alSize;i++)
			storeCoord(i+1,atomList[i]);
		
		// Compute first partial energy
		initialEnergy = a96ff.calculateTotalEnergy(m.actualCoordinates,AANum);
		if (doDihedEnergy)
			initialEnergy[0] += computeOneDihedEnergyDiff(strNumber,dihedNumForCur,0.0f);
		
		// Apply a rotation
		// at1, at2, and at3 don't move, at4 and the atoms
		//  in the atomList rotate
		m.changeTorsion(diAtArray[0],diAtArray[1],diAtArray[2],
			diAtArray[3],stepSize,atomList,alSize,false);
			
		// Compute second energy
		secondEnergy = a96ff.calculateTotalEnergy(m.actualCoordinates,AANum);
		if (doDihedEnergy){
			secondEnergy[0] += computeOneDihedEnergyDiff(strNumber,dihedNumForCur,cumulStep+stepSize);
		}
		
		// Restore coordinates
		restoreCoord(0,diAtArray[3]);
		for(int i=0;i<alSize;i++)
			restoreCoord(i+1,atomList[i]);
		
		
		double thirdEnergy[];
		//Store coordinates
		storeCoord(0,diAtArray[3]);
		for(int i=0;i<alSize;i++)
			storeCoord(i+1,atomList[i]);
		
		// Apply a rotation
		// at1, at2, and at3 don't move, at4 and the atoms
		//  in the atomList rotate
		m.changeTorsion(diAtArray[0],diAtArray[1],diAtArray[2],
			diAtArray[3],-stepSize,atomList,alSize,false);
			
		// Compute second energy
		thirdEnergy = a96ff.calculateTotalEnergy(m.actualCoordinates,AANum);
		if (doDihedEnergy)
			thirdEnergy[0] += computeOneDihedEnergyDiff(strNumber,dihedNumForCur,cumulStep-stepSize);
		
		//Restore coordinates
		restoreCoord(0,diAtArray[3]);
		for(int i=0;i<alSize;i++)
			restoreCoord(i+1,atomList[i]);
		
		if ((initialEnergy[0] > secondEnergy[0])&&(initialEnergy[0] > thirdEnergy[0])){
			if ((initialEnergy[0] - secondEnergy[0])>(initialEnergy[0] - thirdEnergy[0]))
				return stepSize;
			else
				return -stepSize;
		}
		else if(initialEnergy[0] > secondEnergy[0])
			return stepSize;
		else if(initialEnergy[0] > thirdEnergy[0])
			return -stepSize;
		else
			return 0.0f;
	}

	// This function updates the total amount (degrees) in which the given dihedral
	// has moved. It limits this movement to be within +/- lmaxMovement
	protected void updateCumulative(double cumulativeDihedStep[], double newDihedDiff[], int index,
		double lmaxMovement){
		
		if ((cumulativeDihedStep[index] + newDihedDiff[index]) > lmaxMovement)
			newDihedDiff[index] = lmaxMovement - cumulativeDihedStep[index];
		if ((cumulativeDihedStep[index] + newDihedDiff[index]) < -lmaxMovement)
			newDihedDiff[index] = -lmaxMovement - cumulativeDihedStep[index];
		cumulativeDihedStep[index] += newDihedDiff[index];		
	}
	
	//Clears the molecule gradient
	protected void clearMolGradient() {
		int natomsx3 = m.numberOfAtoms * 3;
		m.gradient = new double[natomsx3];
		for(int i=0; i<natomsx3; i++){
			m.gradient[i] = 0;
		}
	}
	
	// Performs a simple steepest descent minimization
	// Assumptions:
	//  -a96ff is current
	//     all atoms of appropriate residues have been assigned
	//     initializeEVCalculation has been called
	//  -m.actualCoords contains current atomic coordinates
	// The user should have set numMinimizationSteps and
	//  initialAngleStepSize (else defaults will be used)
	// Uses specific precomputed nonbonded arrays for each residue (this makes things
	//  run faster)
	// If ligandOnly is true, then only the ligand is allowed to minimize, while
	//		the system residues are fixed
	public void minimize(int numSteps){
	
		/*if ((ligandOnly)&&(ligStrNum == -1)){
			return;
		}*/
	
		float step = initialAngleStepSize;
		double lmaxMovement = maxMovement;
			// maximum degrees by which a torsion can
			//  cumulatively change
		float strRotStep = RotStep;
			// step size for the rigid rotation
			//  of the ligand
		float strTransStep = TransStep;
			// step size in � for the rigid ligand
			//  translation
		double strMaxTrans = MaxTrans;
			// the maximum ligand translation allowed
		
		int i=0;
		boolean done = false;
		double strTorque[] = new double[3];
		double strTrans[] = new double[3];

		strDihedDiff = new double[numberOfStrands][];
		strCumulativeDihedStep = new double[numberOfStrands][];
		for(int str=0; str<numberOfStrands; str++){
			strDihedDiff[str] = new double[numStrDihedrals[str]];
			strCumulativeDihedStep[str] = new double[numStrDihedrals[str]];
		}
		//sysDihedDiff = new double[numSysDihedrals];
		//ligDihedDiff = new double[numLigDihedrals];
		//sysCumulativeDihedStep = new double[numSysDihedrals];
		//ligCumulativeDihedStep = new double[numLigDihedrals];

		strStartCOM = new double[numberOfStrands][3];
		strCurCOM = new double[numberOfStrands][3];
		for(int str=0;str<numberOfStrands;str++){
			strStartCOM[str] = m.getStrandCOM(str);
			strCurCOM[str][0] = strStartCOM[str][0];
			strCurCOM[str][1] = strStartCOM[str][1];
			strCurCOM[str][2] = strStartCOM[str][2];
		}
		/*if(ligStrNum != -1){
			// get the staring COM
			ligStartCOM = m.getStrandCOM(ligStrNum);
			ligCurCOM[0] = ligStartCOM[0];
			ligCurCOM[1] = ligStartCOM[1];
			ligCurCOM[2] = ligStartCOM[2];
		}*/
				
		// Initialize the dihedral movement arrays
		for(int str=0; str<numberOfStrands;str++){
			for(int j=0; j<numStrDihedrals[str];j++){
				strDihedDiff[str][j] = 0.0;
				strCumulativeDihedStep[str][j] = 0.0;
			}
		}
		
		/*for(int j=0;j<numSysDihedrals;j++){
			sysDihedDiff[j] = 0.0;
			sysCumulativeDihedStep[j] = 0.0;
		}
		for(int j=0;j<numLigDihedrals;j++){
			ligDihedDiff[j] = 0.0;
			ligCumulativeDihedStep[j] = 0.0;
		}*/
		
		// If computing dihedral energies initialize them
		if (doDihedEnergy){
			if (!setupDihedralTerms()) //could not initialize dihed energies
				System.exit(1);
		}
		
		float deltaStep = step / numSteps;
		float deltaRotStep = strRotStep / numSteps;
		float deltaTransStep = strTransStep / numSteps;		
		
			// numFlexRes, flexResAtomList, and flexResListSize include the ligand if one exists
			/*if(ligStrNum != -1)
				a96ff.setupPartialArrays(numFlexRes+2,MAX_NUM_ATOMS_DISTAL,flexResAtomList,
					flexResListSize);
			else*/
		a96ff.setupPartialArrays(totalFlexRes+totalTransRotStrands,MAX_NUM_ATOMS_DISTAL,flexResAtomList,
					flexResListSize);

		while(!done){
			
			for(int str=0; str<numberOfStrands;str++){
				for(int j=0;j<numStrDihedrals[str];j++) {
					strDihedDiff[str][j] = computeDihedDiff(strDihedralAtNums[str][j],strDihedralDistal[str][j],
						strNumAtomsDistal[str][j],strDihedToResNum[str][j], step, str, j);
					updateCumulative(strCumulativeDihedStep[str],strDihedDiff[str],j,lmaxMovement);
					applyDihedStep(strDihedralAtNums[str][j],strDihedralDistal[str][j],strNumAtomsDistal[str][j],strDihedDiff[str][j]);
				}
			}
				
			/*for(int j=0;j<numLigDihedrals;j++) {
				ligDihedDiff[j] = computeDihedDiff(ligDihedralAtNums[j],ligDihedralDistal[j],
					ligNumAtomsDistal[j],ligResNumber,step,true,j);
				updateCumulative(ligCumulativeDihedStep,ligDihedDiff,j,lmaxMovement);
				applyDihedStep(ligDihedralAtNums[j],ligDihedralDistal[j],ligNumAtomsDistal[j],ligDihedDiff[j]);
			}*/

			//Translate and rotate the ligand
			/*if(ligStrNum != -1)				
				doLigTransRot(ligTorque, ligTrans, lligRotStep, lligTransStep, lligMaxTrans);*/

			for(int str=0;str<numberOfStrands;str++){
				if(m.strand[str].rotTrans){
					doStrTransRot(str, strTorque, strTrans, strRotStep, strTransStep, strMaxTrans);
				}
			}
			
			/*if(debug){
				String filename = String.format("badVdw/run%1$d/structure%2$02d.pdb",GLOBALNUM,i);
				//Write out pdbs of minimization movement
				m.saveMolecule(filename, 0.0f);
			}*/
			
			i++;
			if(i>=numSteps){
				done = true;
			}

			step -= deltaStep;
			strRotStep -= deltaRotStep;
			strTransStep -= deltaTransStep;
		}
		
		clearMolGradient(); //after minimization is done, clear the molecule gradient
		
		if(debug){
			// Display movement
			for(int str=0; str<numberOfStrands;str++){
				for(int j=0;j<numStrDihedrals[str];j++) {
					System.out.print(strCumulativeDihedStep[str][j] + " ");
				}
			}
			
			GLOBALNUM++;
			
			System.out.println();
		}
	}
////////////////////////////////////////////////////////////////////////////////
//	 End of Minmization Section
////////////////////////////////////////////////////////////////////////////////
	
	//Backup the actual ligand coordinates (assumes there is only 1 residue in the ligand strand)
	private float [] backupStrCoord(int strNumber){
		int index=0;
		float bckpStrCoords[] = new float[m.strand[strNumber].numberOfAtoms*3];
		for (int j=0; j<m.strand[strNumber].numberOfResidues;j++){
			for (int i=0; i<m.strand[strNumber].residue[j].numberOfAtoms; i++){
				int curStrAtom = m.strand[strNumber].residue[j].atom[i].moleculeAtomNumber;
				bckpStrCoords[index*3] = m.actualCoordinates[curStrAtom*3];
				bckpStrCoords[index*3+1] = m.actualCoordinates[curStrAtom*3+1];
				bckpStrCoords[index*3+2] = m.actualCoordinates[curStrAtom*3+2];
				index+=1;
			}
		}
		return bckpStrCoords;
	}
	
	//Restore the backup actual ligand coordinates (assumes there is only 1 residue in the ligand strand)
	private void restoreStrCoord(int strNumber, float bckpStrCoords[]){			
		int index=0;
		for(int j=0; j<m.strand[strNumber].numberOfResidues;j++){
			for (int i=0; i<m.strand[strNumber].residue[j].numberOfAtoms; i++){
				int curStrAtom = m.strand[strNumber].residue[j].atom[i].moleculeAtomNumber;
				m.actualCoordinates[curStrAtom*3] = bckpStrCoords[index*3];
				m.actualCoordinates[curStrAtom*3+1] = bckpStrCoords[index*3+1];
				m.actualCoordinates[curStrAtom*3+2] = bckpStrCoords[index*3+2];
				index+=1;
			}
		}
	}
	
	/*private float [] backupLigCoord(){		
		float bckpLigCoords[] = new float[m.strand[ligStrNum].numberOfAtoms*3];
		for (int i=0; i<m.strand[ligStrNum].numberOfAtoms; i++){
			int curLigAtom = m.strand[ligStrNum].residue[0].atom[i].moleculeAtomNumber;
			bckpLigCoords[i*3] = m.actualCoordinates[curLigAtom*3];
			bckpLigCoords[i*3+1] = m.actualCoordinates[curLigAtom*3+1];
			bckpLigCoords[i*3+2] = m.actualCoordinates[curLigAtom*3+2];
		}
		return bckpLigCoords;
	}
	
	//Restore the backup actual ligand coordinates (assumes there is only 1 residue in the ligand strand)
	private void restoreLigCoord(float bckpLigCoords[]){			
		for (int i=0; i<m.strand[ligStrNum].numberOfAtoms; i++){
			int curLigAtom = m.strand[ligStrNum].residue[0].atom[i].moleculeAtomNumber;
			m.actualCoordinates[curLigAtom*3] = bckpLigCoords[i*3];
			m.actualCoordinates[curLigAtom*3+1] = bckpLigCoords[i*3+1];
			m.actualCoordinates[curLigAtom*3+2] = bckpLigCoords[i*3+2];
		}
	}*/
}
