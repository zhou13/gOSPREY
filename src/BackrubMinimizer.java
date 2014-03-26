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
//	BackrubMinimizer.java
//
//	Version:           2.1 beta
//
//
//	  authors:
// 	  initials    name                 organization                email
//	 ---------   -----------------    ------------------------    ----------------------------
//	  ISG		 Ivelin Georgiev	  Duke University			  ivelin.georgiev@duke.edu
//     KER        Kyle E. Roberts       Duke University         ker17@duke.edu
//     PGC        Pablo Gainza C.       Duke University         pablo.gainza@duke.edu
///////////////////////////////////////////////////////////////////////////////////////////////

/**
* Written by Ivelin Georgiev (2004-2009)
* 
*/

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.ObjectInputStream;
import java.io.PrintStream;
import java.io.Serializable;
import java.util.StringTokenizer;

/** 
 * Handles two types of energy minimization:
 * 		(1) the minimization required for computing the pairwise energy matrices;
 * 		(2) the minimization of a full conformation: the side-chain dihedrals are kept rigid,
 * 				while the backbone is allowed to move using backrub motions
 * Currently, Backrub minimization can be applied only to the system strand of the molecule; 
 * 		the ligand (if present) is allowed to rotate and translate
 * 
 */
public class BackrubMinimizer implements Serializable {
	
	private int MAX_NUM_ATOMS_RES = 30; //the max number of atoms for a given residue; this is increased as needed for big ligands
	
	int mutRes2Strand[] = null;
	int mutRes2StrandMutIndex[] = null;
	
	Molecule m = null;
	Amber96ext a96ff = null;
	
	int numFlexRes[] = null; //the number of residues with flexible side-chains (only in the system strand, and not the ligand)
	int flexResAtomList[][] = new int[0][0]; 
		// each row is a residue and contains
		//  the moleucleatomnumbers located more
		//  distal than the 3rd atom of any dihedral
	int flexResListSize[] = new int[0]; // the number of valid elements in each row of flexResAtomList
	
	int totalFlexRes = 0;
	int totalTransRotStrands = 0;
	int totalNonTransRotStrands = 0;
	
	//int residueMap[] = null; //the mapping between AS and sysStrand-relative residue numbering
	int strandMut[][] = null;
	
	//int ligResNum = -1; // the residue index of the ligand in the flexResAtomList, flexResListSize
	//int ligStrNum = -1; // the ligand strand number (if ligStrNum == -1 then there is no ligand)
	
	//int sysStrNum = -1; //the system strand number
	
	Backrubs br = null; //handles the application of the backrub motions
	
	String backrubFile = null; //the file in which the allowed backrub angles are stored
	float brSamples[][] = null; //the backrub samples (for the big rotation) for each flexible residue (in the system strand only)
	boolean useBRsample[][] = null; //determines if the given backrub for the given residue should be considered
	float theta2[][] = null; //the rotation angles for the first small rotation
	float theta3[][] = null; //the rotation angles for the second small rotation
	
	final float ligRotSize = 0.5f; //initial rotation angle size (in degrees)
	final float ligTransSize = 0.5f; //initial translation size (in angstrom)
	final int numTransRotSteps = 10; //number translation/rotation steps for the ligand
	float ligMaxTrans = 1.2f; //the maximum ligand translation allowed (computed as allowed CA displacement)
	//double initLigCA[] = null;
	
	private double bigE = Math.pow(10, 38);
	
	private double minE = bigE; //the minimum computed energy
	private float minEbrSample[][] = null; //the backrub angles for each of the flexible residues (in the system strand) that give the minimum energy
	private double maxE = -bigE; //the max computed energy (among sterically-allowed conformations)
	private float maxEbrSample[][] = null; //the backrub angles for each of the flexible residues (in the system strand) that give the max energy
	
	//final float idealTauGly = 113.1f; //the ideal tau value for Gly
	//final float idealTauPro = 112.1f; //the ideal tau value for Pro
	final float idealTauOther = 111.0f; //the ideal tau value for all other amino acids (not Gly and not Pro)
	final float tauCutoff = 5.5f; //the cutoff value for tau deviation from the ideal values (if the initial tau values are within the ideal range)
	
	boolean hSteric = true; //determines if hydrogens are used in steric checks
	float overlapThresh = -10000.0f; //the allowed overlap threshold
	
	int numberOfStrands;
	int numMutRes; //KER: total number of mutable residues (i.e. what used to be residueMap.length)
	
	//constructor
	BackrubMinimizer(){
	}
	
	//Initialize for the system strand only
	public void initialize(Molecule mol, Amber96ext theA96ff, int strMut[][], String brFile, boolean hS, float stThresh, int numStrands, boolean readFile) {
		
		m = mol;
		a96ff = theA96ff;
		numberOfStrands = numStrands;		
		
		strandMut = strMut; //the numbering in residueMap[] is system-strand-relative
		initMutRes2Str(strMut);

		//NumMutRes (Don't count non-protein strands)
		numMutRes = 0;
		for(int i=0;i<strandMut.length;i++){
			if(m.strand[i].isProtein)
				numMutRes += strandMut[i].length;
		}
		
		backrubFile = brFile;
		
		hSteric = hS;
		overlapThresh = stThresh;

		// Count number of flexible residues
		numFlexRes = new int[numberOfStrands];;
		for(int str=0; str<numberOfStrands;str++){
			for(int i=0;i<m.strand[str].numberOfResidues;i++){
				if(m.strand[str].residue[i].flexible)
					numFlexRes[str]++;
			}
		}
		
		// 2 is added to numFlexRes so that there is room for the ligand at the
		//  end if there is a ligand present, if there is no ligand then it
		//  doesn't really matter as we won't look at that row
		// The first ligand term includes nonbonded terms for computing energies
		// The second ligand term includes nonbonded terms for computing the gradient
		//  and thus includes terms for all atoms
		totalFlexRes = 0;
		for(int i=0; i<numberOfStrands;i++)
			totalFlexRes += numFlexRes[i];
		
		totalTransRotStrands = 0;
		for(int i=0; i<numberOfStrands;i++)
			if(m.strand[i].rotTrans)
				totalTransRotStrands++;
		
		totalNonTransRotStrands = numberOfStrands-totalTransRotStrands;
		
		flexResAtomList = new int[totalFlexRes+totalTransRotStrands][];
		flexResListSize = new int[totalFlexRes+totalTransRotStrands];
		
		int curNumFlex = 0;
		Residue localRes = null;
		int curTransRotInd = -1;
		int strResNumPP = -1;
		for(int str=0; str<numberOfStrands;str++){
			if(m.strand[str].rotTrans){
				curTransRotInd++;
				strResNumPP = totalFlexRes+curTransRotInd;
				flexResListSize[strResNumPP] = mol.strand[str].numberOfAtoms;
				flexResAtomList[strResNumPP] = new int[flexResListSize[strResNumPP]];
			}
			
			int prevNumAtoms = 0;
			for(int i=0;i<m.strand[str].numberOfResidues;i++){
				localRes = m.strand[str].residue[i];
				
				for(int k=0;k<localRes.numberOfAtoms;k++){
					if(m.strand[str].rotTrans)
						flexResAtomList[strResNumPP][k+prevNumAtoms] = localRes.atom[k].moleculeAtomNumber;
				}
				prevNumAtoms += localRes.numberOfAtoms;
				
				if(localRes.flexible){
					flexResListSize[curNumFlex] = localRes.numberOfAtoms;
					flexResAtomList[curNumFlex] = new int[flexResListSize[curNumFlex]];
					for(int k=0;k<flexResListSize[curNumFlex];k++){
						flexResAtomList[curNumFlex][k] = localRes.atom[k].moleculeAtomNumber;
					}
					curNumFlex++;
				}
			}
		}
		
		if(readFile)
			readBackrubs(); //read in the allowed backrubs
	}
	
	//Initialize for a system and a ligand
	/*public void initialize(Molecule mol, Amber96ext theA96ff, int residueMap[], int sysStrand,
			int ligStrand, String brFile, boolean hS, float stThresh){
		
		MAX_NUM_ATOMS_RES = Math.max(MAX_NUM_ATOMS_RES, mol.strand[ligStrand].residue[0].numberOfAtoms);
	
		// First call the other initialize interface to grab the system information
		initialize(mol,theA96ff,residueMap,sysStrand,brFile,hS,stThresh);
		
		ligStrNum = ligStrand;
		ligResNum = numFlexRes;  // numFlexRes was set in the other initialize call
		
		Residue localRes = m.strand[ligStrNum].residue[0];
		
		flexResListSize[ligResNum + 1] = localRes.numberOfAtoms;
		for(int k=0;k<flexResListSize[ligResNum + 1];k++){
			flexResAtomList[ligResNum + 1][k] = localRes.atom[k].moleculeAtomNumber;
		}

		if(m.strand[ligStrNum].residue[0].flexible){
			flexResListSize[ligResNum] = localRes.numberOfAtoms;
			for(int k=0;k<flexResListSize[ligResNum];k++){
				flexResAtomList[ligResNum][k] = localRes.atom[k].moleculeAtomNumber;
			}
		}
	}*/
	
	//Calls the version below with (shellRun==false) and (templateOnly==false);
	//This version is only called for the pairwise matrix energy precomputation runs involving the shell (SHL-AS, LIG-SHL, and TEMPL)
	public void minimizeFull(){
		minimizeFull(false,false);
	}
	
	//Performs backrubs-based energy minimization for the given system;
	//If (shellRun==true), then this is a pairwise matrix energy computation involving the shell;
	//		if (templateOnly==true), then this is a pairwise matrix energy computation involving *only* the shell
	public void minimizeFull(boolean shellRun, boolean templateOnly){
		
		if (brSamples==null){
			System.out.println("ERROR: Backrub angles data not loaded. Use precomputeBackrubs() to compute the allowed backrub angle.");
			System.exit(1);
		}
		
		br = new Backrubs();
		
		
		
		minE = Math.pow(10, 38);
		minEbrSample = new float[numMutRes][];
		float curEbrSample[][] = new float[numMutRes][];
		maxE = -Math.pow(10, 38);
		maxEbrSample = new float[numMutRes][];
		
		//if (ligStrNum!=-1) //find the initial ligand CA position
		//	initLigCA = getCAcoord(m.strand[ligStrNum].residue[0].moleculeResidueNumber);
		
		setupPartialAmber();
		
		//useBRsample[][] determines if a backrub can be pruned for all rotamer;
		//here, we already have rotamer assignments, so we can determine which backrubs can be pruned for the given rotamers;
		//		this should results in additional pruning, since now we can check some flexible side-chain atoms for steric clashed
		boolean useBRforCurRot[][] = getUseBRforCurRot();
		if (allBRforCurRotPruned(useBRforCurRot)) //all backrubs are pruned for at least one of the residue positions
			return;
		
		
		//generate the excluded atom lists that will be used for partially-assigned backrub conformations during the search
		int excludeList[][][] = new int[numMutRes][][];
		for (int i=0; i<numMutRes; i++)
			excludeList[i] = generateExcludeListSteric(m,strandMut,i,shellRun,templateOnly);
		
		minimizeFullHelper(0,curEbrSample,shellRun,templateOnly,useBRforCurRot,excludeList); //find the best backrubs combination
		
		//Apply the best backrubs combination and translate/rotate the ligand (if present) accordingly
		for (int i=0; i<numMutRes; i++){
			// Check with allowed AAs
			int str = mutRes2Strand[i];
			int strResNum = strandMut[str][mutRes2StrandMutIndex[i]];
			
			Residue curRes = m.strand[str].residue[strResNum];
			if ( (minEbrSample[i]!=null)&&(minEbrSample[i][0]!=0.0f) )
				br.applyBackrub(m, str, curRes.strandResidueNumber, minEbrSample[i][0], false, minEbrSample[i][1], minEbrSample[i][2]);
		}
		for(int str=0;str<numberOfStrands;str++){
			if(m.strand[str].rotTrans){
				doStrTransRot(str);
			}
		}
	}
	
	//Called by minimizeFull()
	//Test all combinations of backrubs for the residues in residueMap[] (in the system strand)
	//		and determine the lowest energy among these backrub combinations
	private void minimizeFullHelper(int curDepth, float curEbrSample[][], boolean shellRun, boolean templateOnly, 
			boolean useBRforCurRot[][], int excludeList[][][]){
		
		if (curDepth==numMutRes){ //a fully-assigned conformation for the system strand
			
			float storedCoord[][][] = new float[numberOfStrands][][];
			for(int str=0;str<numberOfStrands;str++){
				if(m.strand[str].rotTrans){ //there is a ligand, so translate/rotate the ligand
					storedCoord[str] = new float[m.strand[str].numberOfResidues][];
					for(int i=0; i<m.strand[str].numberOfResidues;i++)
						storedCoord[str][i] = storeCoord(m.strand[str].residue[i].moleculeResidueNumber);
				//doLigTransRot(); //this should be un-commented if we want to translate/rotate the ligand for each backrub combination
			}}
			
			double curE[] = null;
			if ( shellRun && !templateOnly ) //shell run, so only compute the energy between the single flexible residue and the shell
				curE = a96ff.calculateTotalEnergy(m.actualCoordinates, 0);
			else 
				curE = a96ff.calculateTotalEnergy(m.actualCoordinates, -1);
			
			if (curE[0]<minE){ //new lowest energy, so update the energy and the corresponding backrub values
				minE = curE[0];
				for (int i=0; i<minEbrSample.length; i++)
					minEbrSample[i] = curEbrSample[i];
			}
			if (curE[0]>maxE){ //new max energy, so update the energy and the corresponding backrub values
				maxE = curE[0];
				for (int i=0; i<maxEbrSample.length; i++)
					maxEbrSample[i] = curEbrSample[i];
			}
			
			for(int str=0;str<numberOfStrands;str++){
				if(m.strand[str].rotTrans){
					for(int i=0; i<m.strand[str].numberOfResidues;i++)
						restoreCoord(m.strand[str].residue[i].moleculeResidueNumber,storedCoord[str][i]); //restore the ligand coordinates to those before minimization
				}
			}
		}
		else {
			// Check with allowed AAs
			int str = mutRes2Strand[curDepth];
			int strResNum = strandMut[str][mutRes2StrandMutIndex[curDepth]];
			
			Residue curRes = m.strand[str].residue[strResNum];
			
			//Backrubs are only applied for the residues currently marked as flexible;
			//If shellRun or templateOnly are true, then this is a pairwise matrix energy precomputation involving the shell, and so
			//		backrubs are applied for all residues in residueMap[], not just the ones marked as flexible
			if (curRes.flexible || shellRun || templateOnly){
				for (int i=0; i<brSamples[curDepth].length; i++){ //apply all backrub samples
					if (useBRsample[curDepth][i] && useBRforCurRot[curDepth][i]) { //allowed backrub
						
						float storedCoord[][] = new float[3][]; //backup the actualCoordinates[] for the three affected residues before the current backrub
						if (brSamples[curDepth][i]!=0.0f){
							storedCoord[0] = storeCoord(m.strand[str].residue[strResNum-1].moleculeResidueNumber); 
							storedCoord[1] = storeCoord(curRes.moleculeResidueNumber);
							storedCoord[2] = storeCoord(m.strand[str].residue[strResNum+1].moleculeResidueNumber);
						
							br.applyBackrub(m, str, curRes.strandResidueNumber, brSamples[curDepth][i], false, theta2[curDepth][i], theta3[curDepth][i]);
						}
						
						if ( ( curRes.flexible && checkSterics(m,str,strResNum,excludeList[curDepth],true) )
								|| ( (!curRes.flexible) && checkSterics(m,str,strResNum,excludeList[curDepth],false) ) ){ //allowed steric
							
							curEbrSample[curDepth] = new float[3];
							curEbrSample[curDepth][0] = brSamples[curDepth][i];
							curEbrSample[curDepth][1] = theta2[curDepth][i];
							curEbrSample[curDepth][2] = theta3[curDepth][i];
							
							minimizeFullHelper(curDepth+1,curEbrSample,shellRun,templateOnly,useBRforCurRot,excludeList); //move to the next residue
						}
						
						if (brSamples[curDepth][i]!=0.0f){
							//restore the actualCoordinates for the three residues to the ones before the current backrub
							restoreCoord(m.strand[str].residue[strResNum-1].moleculeResidueNumber,storedCoord[0]); 
							restoreCoord(curRes.moleculeResidueNumber,storedCoord[1]);
							restoreCoord(m.strand[str].residue[strResNum+1].moleculeResidueNumber,storedCoord[2]);
						}
					}
				}				
			}
			else { //residue currently not flexible, but in brSamples[], so do not apply any backrubs and call at next depth
				minimizeFullHelper(curDepth+1,curEbrSample,shellRun,templateOnly,useBRforCurRot,excludeList); //move to the next residue
			}
		}
	}
	
	//Setup the Amber partial arrays
	private void setupPartialAmber(){
		// numFlexRes, flexResAtomList, and flexResListSize include the ligand if one exists
		/*if(ligStrNum != -1)
			a96ff.setupPartialArrays(numFlexRes+2,MAX_NUM_ATOMS_RES,flexResAtomList,
				flexResListSize);
		else*/
			a96ff.setupPartialArrays(totalFlexRes+totalTransRotStrands,MAX_NUM_ATOMS_RES,flexResAtomList,
				flexResListSize);
	}
	
	//Applies the backrub conformation that corresponds to the max found energy
	//NOTE: minimizeFull() must be run before calling this function
	public void applyMaxBackrub(){
		for (int i=0; i<numMutRes; i++){
			// Check with allowed AAs
			int str = mutRes2Strand[i];
			int strResNum = strandMut[str][mutRes2StrandMutIndex[i]];
			
			Residue curRes = m.strand[str].residue[strResNum];
			if ( (maxEbrSample[i]!=null)&&(maxEbrSample[i][0]!=0.0f) )
				br.applyBackrub(m, str, curRes.strandResidueNumber, maxEbrSample[i][0], false, maxEbrSample[i][1], maxEbrSample[i][2]);
		}
		
		for(int i=0; i<numberOfStrands;i++)
			if(m.strand[i].rotTrans)
				doStrTransRot(i);
	}
	
	//Performs the ligand translation/rotation
	private void doStrTransRot(int strNumber){
		
		float rotStep = ligRotSize;
		float transStep = ligTransSize;
		
		float deltaRotStep = rotStep/numTransRotSteps;
		float deltaTransStep = transStep/numTransRotSteps;
		
		double[] initCOM = m.getStrandCOM(strNumber);
		
		int resNums[] = new int[m.strand[strNumber].numberOfResidues];
		for(int i=0; i<resNums.length;i++)
			resNums[i] = m.strand[strNumber].residue[i].moleculeResidueNumber;
		
		double centOfMass[];
		for (int j=0; j<numTransRotSteps; j++){
			centOfMass = m.getStrandCOM(strNumber);
			for (int curCoord=0; curCoord<3; curCoord++){
				float dTrans = compTrans(resNums,resNums.length,curCoord,transStep);
				if (Math.abs(dTrans)!=0.0)
					updateCumulativeTrans(strNumber, curCoord,dTrans,false, initCOM);
			}
			centOfMass = m.getStrandCOM(strNumber);
			for (int curCoord=0; curCoord<3; curCoord++){
				float axisToRot[] = getRotVector(curCoord); //determine the axis of rotation
				float dRot = compRot(resNums,rotStep,centOfMass,axisToRot);
				if (Math.abs(dRot)!=0.0)
					updateCumulativeRot(resNums,dRot,centOfMass,false,axisToRot);
			}
			
			rotStep -= deltaRotStep;
			transStep -= deltaTransStep;
		}
	}
	
	//Determines the direction for the translation of size transStep for
	//		the numRes number of residues in resNums[] (molecule-relative numbering) in the direction of coord
	private float compTrans(int resNums[], int numRes, int coord, float transStep){
		
		//determine the translation size
		float d[] = new float[3];
		for (int i=0; i<d.length; i++)
			d[i] = 0.0f;
		d[coord] = transStep;	
		
		double initialEnergy[], secondEnergy[], thirdEnergy[];
		float storedCoord[][] = new float[numRes][];
		
		//Store the actualCoordinates for resNum before any changes
		for (int i=0; i<numRes; i++)
			storedCoord[i] = storeCoord(resNums[i]);
		
		
		//Check at initial position
		//TODO: Making these "-1" instead of AAnum makes the calculation a lot slower
		//      which should be fixed by hashing rotatable strand arrays
		initialEnergy = a96ff.calculateTotalEnergy(m.actualCoordinates,-1);
		
		//Check at +transStep rotation
		for (int i=0; i<numRes; i++)
			m.translateResidue(resNums[i], d[0], d[1], d[2], false);	
		secondEnergy = a96ff.calculateTotalEnergy(m.actualCoordinates,-1);		
		for (int i=0; i<numRes; i++)
			restoreCoord(resNums[i],storedCoord[i]);
		
		//Check at -transStep rotation
		d[coord] = -transStep;
		for (int i=0; i<numRes; i++)
			m.translateResidue(resNums[i], d[0], d[1], d[2], false);		
		thirdEnergy = a96ff.calculateTotalEnergy(m.actualCoordinates,-1);		
		for (int i=0; i<numRes; i++)
			restoreCoord(resNums[i],storedCoord[i]);
		
		return getDir(initialEnergy[0],secondEnergy[0],thirdEnergy[0],transStep);
	}
	
	//Checks if a translation of transStep in the direction of coord will move the CA of the ligand further than the limit
	//		and applies the optimal translation if (onlyCheck==false)
	private float updateCumulativeTrans(int strNum, int coord, float transStep, boolean onlyCheck, double[] initCOM){
		
		//int resNum = m.strand[ligStrNum].residue[0].moleculeResidueNumber;
		
		float curCOM[] = m.strand[strNum].getCenterOfMass();
		
		//Determine the new CA if we took this step (the translation is only in the direction of the coord coordinate)
		double tmpCOM[] = new double[3];
		for (int i=0; i<tmpCOM.length; i++)
			tmpCOM[i] = curCOM[i];
		tmpCOM[coord] += transStep;
		
		//Compute how large of a step this would be from the start point
		double dist = getDist(tmpCOM,initCOM);
			
		if (dist > ligMaxTrans){ // if the step would take us too far away then scale it back
			double curDisp[] = new double[3];
			curDisp[0] = curCOM[0] - initCOM[0];
			curDisp[1] = curCOM[1] - initCOM[1];
			curDisp[2] = curCOM[2] - initCOM[2];
			
			double a[] = new double[2];
			int ind = 0;
			for (int i=0; i<3; i++){
				if (i!=coord)
					a[ind++] = curDisp[i] * curDisp[i];
			}
			
			double mt2 = ligMaxTrans * ligMaxTrans;
			double s = 0.0;
			if ((mt2-a[0]-a[1])>0.0) //due to rounding, this number can actually be negative
				s = Math.sqrt(mt2-a[0]-a[1]);
			
			if (curDisp[coord]>=0)
				transStep = (float)(-curDisp[coord] + s);
			else
				transStep = (float)(-curDisp[coord] - s);
				
		}
		
		// compute the translation to get us to the new CA
		float theTranslation[] = new float[3];
		for (int i=0; i<3; i++)
			theTranslation[i] = 0.0f;
		theTranslation[coord] = transStep;

		if (!onlyCheck) //actually apply the new translation
			for(int i=0; i<m.strand[strNum].numberOfResidues;i++)
				m.translateResidue(m.strand[strNum].residue[i].moleculeResidueNumber, theTranslation[0], theTranslation[1], theTranslation[2], false);
		
		return transStep;
	}
	
	//Determines the direction for the rotation of size rotStep centered at center[], for
	//		the numRes number of residues in resNums[]  (molecule-relative numbering) around the axisToRot axis
	private float compRot(int resNums[], float rotStep, double center[], float axisToRot[]){		
		
		double initialEnergy[], secondEnergy[], thirdEnergy[];
		float storedCoord[][] = new float[resNums.length][];
		
		//Store the actualCoordinates for resNum before any changes
		for (int i=0; i<resNums.length; i++)
			storedCoord[i] = storeCoord(resNums[i]);
		
		
		//Check at initial position
		//TODO: Making this "-1" instead of "AAnum" will make it slower. Fix by adding
		//      whole strands to the partial matrix computation
		initialEnergy = a96ff.calculateTotalEnergy(m.actualCoordinates,-1);
		
		//Check at +rotStep rotation
		for (int i=0; i<resNums.length; i++)
			m.rotateResidue(resNums[i], axisToRot[0], axisToRot[1], axisToRot[2], center[0], center[1], center[2], rotStep, false);		
		secondEnergy = a96ff.calculateTotalEnergy(m.actualCoordinates,-1);		
		for (int i=0; i<resNums.length; i++)
			restoreCoord(resNums[i],storedCoord[i]);
		
		//Check at -rotStep rotation
		for (int i=0; i<resNums.length; i++)
			m.rotateResidue(resNums[i], axisToRot[0], axisToRot[1], axisToRot[2], center[0], center[1], center[2], -rotStep, false);		
		thirdEnergy = a96ff.calculateTotalEnergy(m.actualCoordinates,-1);		
		for (int i=0; i<resNums.length; i++)
			restoreCoord(resNums[i],storedCoord[i]);
		
		return getDir(initialEnergy[0],secondEnergy[0],thirdEnergy[0],rotStep);
	}
	
	//Applies a rotation of rotStep degrees for residue resNum (molecule-relative numbering)
	//		around axis axisToRot[] and center center[];
	//If the rotation moves the CA for the given residue further than the allowed limit, 
	//		the rotation is decresed so that CA is moved to the limit;
	//If centerIsCA is true, then the center of rotation is the CA atom of residue resNum, so the rotation will not change the CA position;
	//Only if (checkOnly==false), then the rotation is actually performed
	private float updateCumulativeRot(int[] resNums, float rotStep, double center[], boolean checkOnly, float axisToRot[]){
		
		
		float storedCoord[][] = new float[resNums.length][];
		
		//Store the actualCoordinates for resNum before any changes
		for (int i=0; i<resNums.length; i++)
			storedCoord[i] = storeCoord(resNums[i]);
				
		for (int i=0; i<resNums.length; i++)
			m.rotateResidue(resNums[i], axisToRot[0], axisToRot[1], axisToRot[2], center[0], center[1], center[2], rotStep, false);
		
		//KER: Don't worry about rotation moving center too far since we rotate around the center
		/*double CAcoord[] = getCAcoord(resNum);
		double dist = getDist(CAcoord,initLigCA);
		
		if ( (!centerIsCA) && (dist>ligMaxTrans) ) { //further than allowed, so restore, find the optimal rotStep value, and rotate with that value, if necessary
			
			restoreCoord(resNum,storedCoord);
			
			if ( getDist(initLigCA,getCAcoord(resNum)) > (ligMaxTrans-0.01) ) //the current CA is already approximately at the limit
				return 0.0f;
			
			else {
				CAcoord = getCAcoord(resNum);
				rotStep = binSearchRotStep(CAcoord,axisToRot,center,rotStep,resNum); //find the optimal rotStep value that will rotate the CA to the max limit
				if ((rotStep!=0.0f)&&(!checkOnly)) //apply the rotation
					m.rotateResidue(resNum, axisToRot[0], axisToRot[1], axisToRot[2], center[0], center[1], center[2], rotStep, false);
			}
		}*/
		
		if (checkOnly) //only check, so restore the previous coordinates
			for (int i=0; i<resNums.length; i++)
				restoreCoord(resNums[i],storedCoord[i]);
		
		return rotStep;
	}
	
	//Performs a binary search to find a value for thetaDeg rotation that will rotate the CA for resNum to the maximum allowed limit;
	//The initial coordinates for CA are given in coord[], the rotation is done around axis axisToRot[] and around center center[]
	/*private float binSearchRotStep(double coord[], float axisToRot[], double center[], float origTheta, int resNum){
		
		double newDist;
		
		double newCoord[] = new double[3];
		for (int i=0; i<3; i++)
			newCoord[i] = coord[i];
		
		float thetaDeg = origTheta;
		float thetaStep = origTheta/2.0f;
		
		int curStep = 1;
		while (Math.abs(thetaStep)>0.01){ //thetaStep can be negative
		
			newCoord = rotatePoint(coord,thetaDeg,axisToRot,center);
			newDist = getDist(newCoord,initLigCA);
			if ( (newDist < ligMaxTrans) && (newDist >= (ligMaxTrans-0.01) ) ){
				return thetaDeg;
			}
			else {
				thetaStep = origTheta/(float)Math.pow(2.0,curStep);
				if (newDist>ligMaxTrans)
					thetaDeg -= thetaStep;
				else
					thetaDeg += thetaStep;
			}
			curStep++;
		}
		return 0.0f;
	}*/
	
	//Determines in which direction the best energy is (called by compRot() and compTrans() )
	private float getDir(double e1, double e2, double e3, float step){
		
		if ((e1 > e2)&&(e1 > e3)){
			if ((e1 - e2)>(e1 - e3))
				return step;
			else
				return -step;
		}
		else if (e1 > e2)
			return step;
		else if (e1 > e3)
			return -step;
		else
			return 0.0f;
	}
	
	public double getMagnitude(double a[]){
		double sum = 0.0;
		for(int i=0;i<a.length;i++)
			sum += (a[i]*a[i]);
		return(Math.sqrt(sum));
	}
	
	//Returns the CA coordinates for residue resNum (molecule-relative numbering)
	public double [] getCAcoord(int resNum){
		
		double coord[] = new double[3];
		
		Residue r = m.residue[resNum];
		if (m.strand[r.strandNumber].isProtein) { //use the CA atom
			for (int i=0; i<r.numberOfAtoms; i++){
				if (r.atom[i].name.equalsIgnoreCase("CA")){
					int curAtom = r.atom[i].moleculeAtomNumber;
					coord[0] = m.actualCoordinates[curAtom*3];
					coord[1] = m.actualCoordinates[curAtom*3+1];
					coord[2] = m.actualCoordinates[curAtom*3+2];
					break;
				}
			}
		}
		else { //compute and return the geometric center of the residue
			coord = m.getResidueGC(resNum);
		}
		return coord;
	}
	
	//Determines the vector around which the rotation is to be performed;
	//Currently, the vector is the coordinate axis specified by axisNum
	public float [] getRotVector(int axisNum){
		
		float axisToRot[] = new float[3];
		for (int i=0; i<axisToRot.length; i++)
			axisToRot[i] = 0.0f;
		axisToRot[axisNum] = 1.0f;
		
		return axisToRot;
	}
	
	//Returns the distance between the two points given by coordinates coord1[] and coord2[] (double version)
	private double getDist(double coord1[], double coord2[]){
		
		double rijx, rijy, rijz, rij, rij2;
		
		rijx = coord1[0] - coord2[0];
		rijy = coord1[1] - coord2[1];
		rijz = coord1[2] - coord2[2];
		rij2 = rijx * rijx + rijy * rijy + rijz * rijz;
		rij = Math.sqrt(rij2);
		
		return rij;
	}
	
	//Rotates the point with coordinates coord[] with thetaDeg degrees around the axis of rotation d[] and around center c[];
	//Returns the new coordinates in newCoord[] and does not modify the original coordinates
	private double [] rotatePoint(double coord[], float thetaDeg, float d[], double c[]){
		
		double newCoord[] = new double[3];
		
		double tx,ty,tz;

		float[][] rot_mtx = new float[3][3];
		RotMatrix rM = new RotMatrix();
		rM.getRotMatrix(d[0],d[1],d[2],thetaDeg,rot_mtx);
			
		tx=coord[0] - c[0];
		ty=coord[1] - c[1];
		tz=coord[2] - c[2];

		newCoord[0] = (float)(tx * rot_mtx[0][0] + ty * rot_mtx[0][1] + tz * rot_mtx[0][2] + c[0]);
		newCoord[1] = (float)(tx * rot_mtx[1][0] + ty * rot_mtx[1][1] + tz * rot_mtx[1][2] + c[1]);
		newCoord[2] = (float)(tx * rot_mtx[2][0] + ty * rot_mtx[2][1] + tz * rot_mtx[2][2] + c[2]);
		
		return newCoord;
	}
	
	//Returns the current actualCoordinates for residue resNum (molecule-relative numbering)
	private float [] storeCoord(int resNum){
		
		float storedCoord[] = new float[m.residue[resNum].numberOfAtoms * 3];
		for (int i=0; i<m.residue[resNum].numberOfAtoms; i++){
			
			int curAtom = m.residue[resNum].atom[i].moleculeAtomNumber;
			storedCoord[i*3] = m.actualCoordinates[curAtom*3];
			storedCoord[i*3+1] = m.actualCoordinates[curAtom*3+1];
			storedCoord[i*3+2] = m.actualCoordinates[curAtom*3+2];
		}
		return storedCoord;
	}
	
	//Sets the actualCoordinates for residue resNum (molecule-relative numbering) to the ones in storedCoord[]
	private void restoreCoord(int resNum, float storedCoord[]){
		
		for (int i=0; i<m.residue[resNum].numberOfAtoms; i++){
			
			int curAtom = m.residue[resNum].atom[i].moleculeAtomNumber;
			m.actualCoordinates[curAtom*3] = storedCoord[i*3];
			m.actualCoordinates[curAtom*3+1] = storedCoord[i*3+1];
			m.actualCoordinates[curAtom*3+2] = storedCoord[i*3+2];
		}
	}
	
	//Determines which of the possible backrubs can be pruned for the current rotamer assignment
	private boolean [][] getUseBRforCurRot(){
		
		//Generate the list of atoms to exclude from steric overlap checks (this only considers residues in the system strand)
		int excludeList[][] = generateExcludeListSteric(m,strandMut,-1,false,false);
		
		boolean useBRforCurRot[][] = new boolean[brSamples.length][];
		
		for (int curDepth=0; curDepth<useBRforCurRot.length; curDepth++){
			// Check with allowed AAs
			int str = mutRes2Strand[curDepth];
			int strResNum = strandMut[str][mutRes2StrandMutIndex[curDepth]];
			
			Residue curRes = m.strand[str].residue[strResNum];
			
			if (curRes.flexible) { //perform the check only for flexible residues, since rotamers are assigned only for the flexible residues
			
				useBRforCurRot[curDepth] = new boolean[brSamples[curDepth].length];
				
				for (int i=0; i<useBRforCurRot[curDepth].length; i++){ //apply all backrub samples
					if (useBRsample[curDepth][i]) {
						
						try{
						float storedCoord[][] = new float[3][]; //backup the actualCoordinates[] for the three affected residues before the current backrub
						if (brSamples[curDepth][i]!=0.0f){
							storedCoord[0] = storeCoord(m.strand[str].residue[strResNum-1].moleculeResidueNumber); 
							storedCoord[1] = storeCoord(curRes.moleculeResidueNumber);
							storedCoord[2] = storeCoord(m.strand[str].residue[strResNum+1].moleculeResidueNumber);
							
							br.applyBackrub(m, str, curRes.strandResidueNumber, brSamples[curDepth][i], false, theta2[curDepth][i], theta3[curDepth][i]);
						}
						
						
						if (checkSterics(m, str, strResNum, excludeList, true)) //allowed steric
							useBRforCurRot[curDepth][i] = true;
						else
							useBRforCurRot[curDepth][i] = false;
						
						//restore the actualCoordinates for the three residues to the ones before the current backrub
						if (brSamples[curDepth][i]!=0.0f){
							restoreCoord(m.strand[str].residue[strResNum-1].moleculeResidueNumber,storedCoord[0]); 
							restoreCoord(curRes.moleculeResidueNumber,storedCoord[1]);
							restoreCoord(m.strand[str].residue[strResNum+1].moleculeResidueNumber,storedCoord[2]);
						}
						}
						catch(Exception E){
							System.out.println("Delete me please!!");
							E.printStackTrace();
						}
						
					}
					else
						useBRforCurRot[curDepth][i] = false;
				}
			}
			else //not a flexible residue
				useBRforCurRot[curDepth] = useBRsample[curDepth];
		}
		
		return useBRforCurRot;
	}
	
	//Returns true if all backrubs for at least one flexible residue position are pruned
	private boolean allBRforCurRotPruned(boolean useBRforCurRot[][]){
		
		for (int i=0; i<useBRforCurRot.length; i++){
			// Check with allowed AAs
			int str = mutRes2Strand[i];
			int strResNum = strandMut[str][mutRes2StrandMutIndex[i]];
			
			if (m.strand[str].residue[strResNum].flexible){ //only check for the flexible residues
				boolean allPruned = true;
				for (int j=0; j<useBRforCurRot[i].length; j++){
					if (useBRforCurRot[i][j])
						allPruned = false;
				}
				if (allPruned)
					return true;
			}
		}
		return false;
	}
	
	//Performs a precomputation of the allowed backrub angles (brSamples[][], theta2[][], and theta3[][]);
	//		Outputs the computed angles to the backrubFile, which is then read upon each initialize
	public void precomputeBackrubs(int numBackrubSamples, float backrubStepSize){
		
		a96ff.calculateTypesWithTemplates(); //assign force field atom types and parameters
		
		//Compute the big rotation angles
		brSamples = new float[numMutRes][2*numBackrubSamples+1];
		useBRsample = new boolean[numMutRes][brSamples[0].length];
		for (int i=0; i<brSamples.length; i++){
			for (int j=0; j<brSamples[i].length; j++){
				brSamples[i][j] = (j-numBackrubSamples)*backrubStepSize;
				useBRsample[i][j] = true; // initially, all backrubs are allowed for all residues
			}
		}
		
		//Compute the two small rotation angles and prune the infeasible big/small angle choices
		br = new Backrubs();
		theta2 = new float[brSamples.length][brSamples[0].length];
		theta3 = new float[brSamples.length][brSamples[0].length];
		
		//Generate the list of atoms to exclude from steric overlap checks (this only considers residues in the system strand)
		int excludeList[][] = generateExcludeListSteric(m,strandMut,-1,false,false);
		
		for (int i=0; i<brSamples.length; i++){
			for (int j=0; j<brSamples[0].length; j++){
				if (brSamples[i][j]==0.0f){
					theta2[i][j] = 0.0f;
					theta3[i][j] = 0.0f;
				}
				else if (useBRsample[i][j]){ //for each allowed big rotation, determine the rotation angles for the two small rotations
					// Check with allowed AAs
					int str = mutRes2Strand[i];
					int strResNum = strandMut[str][mutRes2StrandMutIndex[i]];
					
					
					Residue curRes = m.strand[str].residue[strResNum];
					
					float initTaus[] = new float [3]; //get the initial tau values
					initTaus[0] = getTau(m,str,strResNum-1);
					initTaus[1] = getTau(m,str,strResNum);
					initTaus[2] = getTau(m,str,strResNum+1);
					
					float theta[] = br.applyBackrub(m, str, curRes.strandResidueNumber, brSamples[i][j], true, 0.0f, 0.0f);
					
					theta2[i][j] = theta[0];
					theta3[i][j] = theta[1];
					
					if (!isBRfeasible(m,str,strResNum,initTaus,excludeList))
						useBRsample[i][j] = false;
					
					m.updateCoordinates(); //restore the actualCoordinates[] to the values before the backrub
				}
			}
		}
		
		outputBackrubs(); //output the computed allowed backrub angles to backrubFile
	}
	
	//Checks if the conformation for residue resNum of strand strNum in molecule m is feasible;
	//Currently, checks for significant tau deviations for the tri-peptide between (resnum-1), resNum, (resNum+1)
	//		and for steric clashes with the fixed part of the molecule
	private boolean isBRfeasible(Molecule m, int strNum, int resNum, float initTaus[], int excludeList[][]){
		
		if ( (!checkSterics(m, strNum, resNum-1, excludeList, false))
				|| (!checkSterics(m, strNum, resNum, excludeList, false))
				|| (!checkSterics(m, strNum, resNum+1, excludeList, false)) ) //unallowed steric for one of the three residues involved in the current backrub
					return false;
		
		if ( ! checkAllowedTaus(m, strNum, resNum, initTaus) )
			return false;
		
		return true;
	}
	
	//Checks for steric overlap between:
	//		1) the backbone atoms of residue resNum (strand-relative numbering) in strand strNum in molecule m, and
	//		2) all atoms in the molecule not belonging to the residues in excludeRes[] (strand-relative numbering) in strand strNum
	//If hSteric is true, then hydrogens are used in the steric checks;
	//If a ligand is present, do not check for overlap with the ligand atoms, since the ligand translates/rotates in the last step of the minimization procedure
	private boolean checkSterics(Molecule m, int strNum, int resNum, int excludeList[][], boolean useResSC){
		
		ProbeStericCheck psc = new ProbeStericCheck();
		
		Residue res = m.strand[strNum].residue[resNum];
		
		for(int i=0;i<res.numberOfAtoms;i++) {
			Atom a1 = res.atom[i];
			if ( (hSteric) || (!a1.elementType.equalsIgnoreCase("H")) ) {
				if (a1.getIsBBatom() || a1.name.equalsIgnoreCase("CB") || useResSC){ //include the CB of he current residue for the steric check
					for(int q=0;q<m.numberOfStrands;q++) {
						if (!m.strand[q].rotTrans) {//not the ligand strand (i.e. strand does not rotate/translate)
							int resToCheck = m.strand[q].numberOfResidues;
							for(int w=0;w<resToCheck;w++) {
								if( !( (q==strNum) && (w==resNum) ) ) {
									for(int t=0;t<m.strand[q].residue[w].numberOfAtoms;t++) {
										Atom a2 = m.strand[q].residue[w].atom[t];
										if ( !isInList(a2.moleculeAtomNumber,excludeList) ){ //atom not excluded
											if ( (hSteric) || (!a2.elementType.equalsIgnoreCase("H")) )
												if (!psc.isAllowedSteric(m, a1, a2, overlapThresh))
													return false;
										}
									}
								}
							}
						}
					}
				}
			}		
		}
	
		// If you got here then everything passed
		return true;
		
	}
	
	//Checks if the tau angles for residues (resNum-1), resNum, and (resNum+1) (strand-relative numbering)
	//		of strand strNum in molecule m, are within the allowed range
	private boolean checkAllowedTaus(Molecule m, int strNum, int resNum, float initTaus[]){
		
		float curTaus[] = new float[3];
		curTaus[0] = getTau(m,strNum,resNum-1);
		curTaus[1] = getTau(m,strNum,resNum);
		curTaus[2] = getTau(m,strNum,resNum+1);
		
		for (int i=0; i<curTaus.length; i++) {
			
			float curIdealTau = 10000.0f;
			/*if (m.strand[strNum].residue[resNum-1+i].name.equalsIgnoreCase("GLY"))
				curIdealTau = idealTauGly;
			else if (m.strand[strNum].residue[resNum-1+i].name.equalsIgnoreCase("PRO"))
				curIdealTau = idealTauPro;
			else */
				curIdealTau = idealTauOther;
			
			if (Math.abs(curIdealTau-curTaus[i])>tauCutoff){ //new tau value not within ideal range
				
				if ( (Math.abs(curIdealTau-initTaus[i])<tauCutoff) ) //initial tau value within ideal range, so do not allow this backrub
					return false;
				else if ( (Math.abs(curIdealTau-initTaus[i])>tauCutoff) && ( Math.abs(curIdealTau-initTaus[i])<Math.abs(curIdealTau-curTaus[i]) ) )
					return false; //(initial not within ideal range) and (initial closer to ideal than current) and (change from initial to current greater than cutoff)
			}
		}
		return true;
	}
	
	//Computes the tau value (from the actualCoordinates[]) for residue resNum (strand-relative numbering)	of strand strNum in molecule m
	private float getTau(Molecule m, int strNum, int resNum){
		
		Residue res = m.strand[strNum].residue[resNum];
		Atom a[] = new Atom[3];
		for (int i=0; i<res.numberOfAtoms; i++){ //find the N, CA, and C atoms
			if (res.atom[i].name.equalsIgnoreCase("N"))
				a[0] = res.atom[i];
			else if (res.atom[i].name.equalsIgnoreCase("CA"))
				a[1] = res.atom[i];
			else if (res.atom[i].name.equalsIgnoreCase("C"))
				a[2] = res.atom[i];
		}
		
		float paC[][] = new float[3][3]; //get the actualCoordinates[] for the atoms
		for (int i=0; i<3; i++)
			paC[i] = m.getActualCoord(a[i].moleculeAtomNumber);
		
		Atom pa[] = new Atom[3];
		for (int i=0; i<3; i++)
			pa[i] = new Atom("",paC[i][0],paC[i][1],paC[i][2]);
		
		return ((float)pa[2].angle(pa[0], pa[1]));
	}
	
	//Generate an exclusion list of atoms not to be checked in steric checks, since their position is unknown due to possible backrubs;
	//The residue numbers in resMap[] are strand-relative
	private int [][] generateExcludeListSteric(Molecule m, int strandMut[][], int aboveInd, boolean shellRun, boolean templateOnly){
		
		int excludeList[][] = new int[numMutRes][];
		int ctr = -1;
		for (int str=0; str<strandMut.length;str++){
			if(m.strand[str].isProtein){
				for (int i=0; i<strandMut[str].length; i++){
					ctr++;
					
					//KER: if this is the most C-terminal or N-terminal then don't allow backrub to occur
					if(strandMut[str][i]-1 < 0 || strandMut[str][i]+1 >= m.strand[str].residue.length){
						//KER: Don't allow backrubs so turn all backrubs false except the one with 0 change
						for(int q=0; q<useBRsample[ctr].length;q++)
							useBRsample[ctr][q] = false;
						            
						int midPoint = (useBRsample[ctr].length-1)/2;
						useBRsample[ctr][midPoint] = true;
						excludeList[ctr] = new int[0];
						continue;
					}
					
					Residue prevRes = m.strand[str].residue[strandMut[str][i]-1];
					Residue curRes = m.strand[str].residue[strandMut[str][i]];
					Residue nextRes = m.strand[str].residue[strandMut[str][i]+1];
					
					
					
					if ( (ctr>aboveInd) || ( (ctr<aboveInd) && (!curRes.flexible) && (!(shellRun || templateOnly)) ) ) { //exclude all
						
						excludeList[ctr] = new int[2 + curRes.numberOfAtoms + 2];
						
						excludeList[ctr][0] = prevRes.getAtomNameToMolnum("C");
						excludeList[ctr][1] = prevRes.getAtomNameToMolnum("O");
						
						for (int j=0; j<curRes.numberOfAtoms; j++)
							excludeList[ctr][2+j] = curRes.atom[j].moleculeAtomNumber;
						
						excludeList[ctr][curRes.numberOfAtoms+2] = nextRes.getAtomNameToMolnum("N");
						excludeList[ctr][curRes.numberOfAtoms+3] = nextRes.getAtomNameToMolnum("H");
					}
					else if ( (ctr<aboveInd) && (!curRes.flexible) ) { //residue in resMap[], but not flexible, so exclude only the side-chain atoms of this residue
						
						excludeList[ctr] = new int[curRes.numberOfAtoms];
						for (int j=0; j<curRes.numberOfAtoms; j++){
							if ( ! (curRes.atom[j].getIsBBatom() || curRes.atom[j].name.equalsIgnoreCase("CB")) )
								excludeList[ctr][j] = curRes.atom[j].moleculeAtomNumber;
							else
								excludeList[ctr][j] = -1;
						}				
					}
					else //do not exclude anything, since no atoms will change their position
						excludeList[ctr] = new int[0];
				}
			}
		}
		
		return excludeList;
	}
	
	//Checks if a is in list[]
	private boolean isInList(int a, int list[][]){
		for (int i=0; i<list.length; i++){
			for (int j=0; j<list[i].length; j++){
				if (list[i][j]==a)
					return true;
			}
		}
		return false;
	}
	
	//Output the computed backrub angles to backrubFile
	private void outputBackrubs(){
		PrintStream logPS = setupOutputFile(backrubFile);
		logPS.println(brSamples.length+" "+brSamples[0].length);
		for (int i=0; i<brSamples.length; i++){
			for (int j=0; j<brSamples[i].length; j++){
				if (useBRsample[i][j]){
					logPS.println(i+" "+j+" "+brSamples[i][j]+" "+theta2[i][j]+" "+theta3[i][j]);
				}
			}
		}
		logPS.flush();
		logPS.close();
	}
	
	//Read in the allowed backrubs from backrubFile
	private void readBackrubs(){
		
		BufferedReader bufread = null;
		try {
			File file = new File(backrubFile);
			FileReader fr = new FileReader(file);
			bufread = new BufferedReader(fr);
		}
		catch (FileNotFoundException e) {
			System.out.println("Can't find backrubfile: "+backrubFile);
			System.out.println("Please fix backrubFile parameter or run precomputeBackrubs to get backrub angles");
			System.exit(0);
			return;
		}

		boolean done = false;
		String str = null;		
		
		int curLine = 0;
		while (!done) {
			try {
				str = bufread.readLine();
			}
			catch ( Exception e ){
			 System.out.println("ERROR: An error occurred while reading input");
			 System.exit(0);
			}

			if (str == null) // stop if we've reached EOF
				done = true;
			else if (getToken(str,1).charAt(0)=='%') //skip comment lines
				continue;
			else {	
				if (curLine==0) {//first line is the dimensions of brSamples[][], theta2[][], and theta3[][]
					int s1 = new Integer(getToken(str,1)).intValue();
					int s2 = new Integer(getToken(str,2)).intValue();
					brSamples = new float[s1][s2];
					theta2 = new float[s1][s2];
					theta3 = new float[s1][s2];
					useBRsample = new boolean[s1][s2];
					for (int i=0; i<useBRsample.length; i++){
						for (int j=0; j<useBRsample[i].length; j++){
							useBRsample[i][j] = false; //we will only set to true the entries that are read in below
						}
					}
				}
				else { //each line is in the format: ind1 ind2 brSamples[ind1][ind2] theta2[ind1][ind2] theta3[ind1][ind2]
					int ind1 = new Integer(getToken(str,1)).intValue();
					int ind2 = new Integer(getToken(str,2)).intValue();
					brSamples[ind1][ind2] = new Float(getToken(str,3)).floatValue();
					theta2[ind1][ind2] = new Float(getToken(str,4)).floatValue();
					theta3[ind1][ind2] = new Float(getToken(str,5)).floatValue();
					useBRsample[ind1][ind2] = true;
				}
				curLine++;
			}
		}		
		try {bufread.close();} catch(Exception e){} //done reading them in
	}
	
	//Setup the file with name filename for output
	private PrintStream setupOutputFile(String fileName){
		PrintStream logPS = null; //the output file for conf info
		try {			
			FileOutputStream fileOutputStream = new FileOutputStream(fileName);
			BufferedOutputStream bufferedOutputStream = new BufferedOutputStream( fileOutputStream );
			logPS = new PrintStream( bufferedOutputStream );
		}
		catch (Exception ex) {
			System.out.println("ERROR: An exception occured while opening log file");
		}
		return logPS;
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
				// System.out.println("ERROR: Unable to access parameter " + x);
				return(new String(""));
			}
		}
		
		if (st.hasMoreTokens())		
			return(st.nextToken());
		return(new String(""));

	} // end getToken
	
	//Computes the min and max intra-energy for resMapPos (index into residueMap[]);
	//Since the backbone atoms of a given residue are included in the INTRA energy precomputation, a backrub
	//		at that residue can change the intra-energy, so we need to compute the lower/upper bounds on the intra energies
	public float [] getMinMaxIntraEnergyBR(int resMapPos){
		
		br = new Backrubs();
		
		double minIntra = bigE;
		double maxIntra = -bigE;
		
		int str = mutRes2Strand[resMapPos];
		int strResNum = strandMut[str][mutRes2StrandMutIndex[resMapPos]];
		
		Residue curRes = m.strand[str].residue[strResNum];
		
		for (int i=0; i<brSamples[resMapPos].length; i++){ //apply all backrub samples
			if (useBRsample[resMapPos][i]) { //allowed backrub
		
				float storedCoord[][] = new float[3][]; //backup the actualCoordinates[] for the three affected residues before the current backrub
				if (brSamples[resMapPos][i]!=0.0f){
					storedCoord[0] = storeCoord(m.strand[str].residue[strResNum-1].moleculeResidueNumber); 
					storedCoord[1] = storeCoord(curRes.moleculeResidueNumber);
					storedCoord[2] = storeCoord(m.strand[str].residue[strResNum+1].moleculeResidueNumber);
				
					br.applyBackrub(m, str, curRes.strandResidueNumber, brSamples[resMapPos][i], false, theta2[resMapPos][i], theta3[resMapPos][i]);
				}
				
				double curE[] = a96ff.calculateTotalEnergy(m.actualCoordinates, -1);
				
				minIntra = Math.min(curE[0], minIntra);
				maxIntra = Math.max(curE[0], maxIntra);
				
				//restore the actualCoordinates for the three residues to the ones before the current backrub
				if (brSamples[resMapPos][i]!=0.0f){
					restoreCoord(m.strand[str].residue[strResNum-1].moleculeResidueNumber,storedCoord[0]); 
					restoreCoord(curRes.moleculeResidueNumber,storedCoord[1]);
					restoreCoord(m.strand[str].residue[strResNum+1].moleculeResidueNumber,storedCoord[2]);
				}
			}
		}
		
		float e[] = new float[2];
		e[0] = (float)minIntra;
		e[1] = (float)maxIntra;
		
		return e;
	}
	
	
	public void initMutRes2Str(int strandMut[][]){
		int totalLength=0;
		for(int i=0;i<strandMut.length;i++)
			for(int j=0;j<strandMut[i].length;j++)
				totalLength++;
		
		mutRes2Strand = new int[totalLength];
		mutRes2StrandMutIndex = new int[totalLength];
		int ctr=0;
		for(int i=0;i<strandMut.length;i++)
			for(int j=0;j<strandMut[i].length;j++){
				mutRes2Strand[ctr] = i;
				mutRes2StrandMutIndex[ctr] = j;
				ctr++;
			}
	}
	
}
