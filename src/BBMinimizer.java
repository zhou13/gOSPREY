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
//	BBMinimizer.java
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

import java.io.Serializable;

/** 
 * Handles two types of backbone energy minimization for BD:
 * 		(1) the minimization required for computing the pairwise energy matrices;
 * 		(2) the minimization of a full conformation: the side-chain dihedrals are kept rigid,
 * 				while the backbone dihedrals are allowed to move within given limits
 * 
 */
public class BBMinimizer implements Serializable {
	
	private int MAX_NUM_ATOMS_RES = 30; //the max number of atoms for a given residue; this is increased as needed for big ligands
	
	Molecule m = null;
	Amber96ext a96ff = null;
	
	int numFlexRes[] = null; //the number of residues with flexible side-chains (only in the system strand, and not the ligand)
	int totalFlexRes = 0;
	int totalTransRotStrands = 0;
	int totalNonTransRotStrands = 0;
	int flexResAtomList[][] = new int[0][0]; 
		// each row is a residue and contains
		//  the moleucleatomnumbers located more
		//  distal than the 3rd atom of any dihedral
	int flexResListSize[] = new int[0]; // the number of valid elements in each row of flexResAtomList
	
	//int residueMap[] = null; //the mapping between AS and molecule-relative residue numbering
	int strandMut[][] = null;
	
	//int ligResNum = -1; // the residue index of the ligand in the flexResAtomList, flexResListSize
	//int ligStrNum = -1; // the ligand strand number (if ligStrNum == -1 then there is no ligand)
	
	//int sysStrNum = -1; //the system strand number
	int numberOfStrands;
	final int numStepsFiPsi = 15; //the number of phi/psi minimization steps to be performed	
	final float sysMaxDihedRot = 3.0f; //the maximum rotation from the initial phi/sp values
	private float sysFiPsiStep = 0.8f; //initial step size for phi/psi changes; used for full conformation energy minimization
	
	float molCurFiPsiDisp[][] = null; //the current displacement from the initial phi/psi angle for each residue
	
	float ligRotSize = 0.5f; //initial rotation angle size (in degrees)
	float ligTransSize = 0.5f; //initial translation size (in angstrom)
	
	float sysMaxTrans = 1.5f;//2.0f; //the maximum displacement (in angstrom) from the initial CA
	private int molAtNumCA[] = null; //the molecule atom numbers for the CA atoms for each residue
	double molStartCA[][] = null; // the position of the initial CA for each residue
	private double molCurCA[][] = null; // the current CA coordinates for each residue
	
	
	//constructor
	BBMinimizer(){
	}
	
	//Initialize for the system strand only
	public void initialize(Molecule mol, Amber96ext theA96ff, int strMut[][], int numStrands) {
		
		m = mol;
		a96ff = theA96ff;
		numberOfStrands = numStrands;
		strandMut = strMut;
		//sysStrNum = sysStrand;
		//residueMap = resMap;

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
		int curTransRotInd = -1;
		int strResNumPP = -1;
		for(int str=0;str<numberOfStrands;str++){
			if(m.strand[str].rotTrans){
				curTransRotInd++;
				strResNumPP = totalFlexRes+curTransRotInd;
				flexResListSize[strResNumPP] = mol.strand[str].numberOfAtoms;
				flexResAtomList[strResNumPP] = new int[flexResListSize[strResNumPP]];
			}
			
			int prevNumAtoms = 0;
			for(int i=0; i<m.strand[str].numberOfResidues;i++){
				Residue localRes = m.strand[str].residue[i];
				
				for(int k=0;k<localRes.numberOfAtoms;k++){
					if(m.strand[str].rotTrans)
						flexResAtomList[strResNumPP][k+prevNumAtoms] = localRes.atom[k].moleculeAtomNumber;
				}
				prevNumAtoms += localRes.numberOfAtoms;
				if(m.strand[str].residue[i].flexible){
					flexResListSize[curNumFlex] = localRes.numberOfAtoms; 
					flexResAtomList[curNumFlex] = new int[flexResListSize[curNumFlex]];
					for(int k=0;k<flexResListSize[curNumFlex];k++){
						flexResAtomList[curNumFlex][k] = localRes.atom[k].moleculeAtomNumber;
					}
					curNumFlex++;
				}
			}
		}
		
		
		/*int curNumFlex = 0;
		Residue localRes = null;
		for(int i=0;i<m.strand[sysStrNum].numberOfResidues;i++){
			localRes = m.strand[sysStrNum].residue[i];
			if(localRes.flexible){
				flexResListSize[curNumFlex] = localRes.numberOfAtoms;
				for(int k=0;k<flexResListSize[curNumFlex];k++){
					flexResAtomList[curNumFlex][k] = localRes.atom[k].moleculeAtomNumber;
				}
				curNumFlex++;
			}
		}*/
		
		molAtNumCA = getMolAtNumCA(); //this must be set before setting molStartCA and molCurCA
		
		molStartCA = new double[m.numberOfResidues][3];
		molCurCA = new double[m.numberOfResidues][3];
		for (int i=0; i<m.numberOfResidues; i++){
			molStartCA[i] = getCAcoord(i);
			molCurCA[i][0] = molStartCA[i][0];
			molCurCA[i][1] = molStartCA[i][1];
			molCurCA[i][2] = molStartCA[i][2];
		}
		
		molCurFiPsiDisp = new float[m.numberOfResidues][2];
		for (int i=0; i<molCurFiPsiDisp.length; i++){
			molCurFiPsiDisp[i][0] = 0.0f;
			molCurFiPsiDisp[i][1] = 0.0f;
		}
	}
	
	//Initialize for a system and a ligand
	/*public void initialize(Molecule mol, Amber96ext theA96ff, int residueMap[], int sysStrand,
			int ligStrand){
		
		MAX_NUM_ATOMS_RES = Math.max(MAX_NUM_ATOMS_RES, mol.strand[ligStrand].residue[0].numberOfAtoms);
	
		// First call the other initialize interface to grab the system information
		initialize(mol,theA96ff,residueMap,sysStrand);
		
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
	
	//Minimizes the energy of a fully-assigned AA sequence with fully-assigned rotamer identities,
	//		starting with a given template
	public void minimizeFull(boolean pemComp){
		
		setupPartialAmber();
		updateCurCA();
		
		float fiPsiStep = sysFiPsiStep;		
		float rotStep = ligRotSize;
		float transStep = ligTransSize;
		
		float deltaFiPsiStep = fiPsiStep/numStepsFiPsi;
		float deltaRotStep = rotStep/numStepsFiPsi;
		float deltaTransStep = transStep/numStepsFiPsi;
		
		Backbone bb = new Backbone();
		
		for (int i=0; i<numStepsFiPsi; i++) {
			for (int str=0; str<numberOfStrands;str++){
				for (int j=0; j<strandMut[str].length; j++){ //apply the current phi/psi step for each of the flexible residues
					int asResNum = strandMut[str][j]; //strandMut is strandResidue specific
				for (int a=0; a<2; a++){ // (a==0) applies a phi change; (a==1) applies a psi change
						
						if( (!m.strand[str].residue[asResNum].nterm && !(asResNum==0) && a==0) ||
								(!m.strand[str].residue[asResNum].cterm && !(asResNum==(m.strand[str].numberOfResidues-1)) && a==1) ){
							float change = compFiPsi(asResNum,str,a,fiPsiStep,bb,pemComp);
					if (Math.abs(change)>0.001){ //smaller changes should not be applies
								bb.applyFiPsi(m,str,asResNum,change,a);
								molCurFiPsiDisp[m.strand[str].residue[strandMut[str][j]].moleculeResidueNumber][a] += change;
					}
				}
			}
				}
			}
			int curTransRotInd = totalFlexRes-1;
			for (int str=0; str<numberOfStrands;str++){
				if(m.strand[str].rotTrans){
					curTransRotInd++;
					int resNums[] = new int[m.strand[str].numberOfResidues];
					for(int q=0;q<resNums.length;q++)
						resNums[q] = m.strand[str].residue[q].moleculeResidueNumber;
				
				for (int curCoord=0; curCoord<3; curCoord++){
						float dTrans = compTrans(resNums,resNums.length,curCoord,transStep,curTransRotInd);
						if (Math.abs(dTrans)!=0.0){
							for(int k=0;k<resNums.length;k++)
								updateCumulativeTrans(resNums[k],curCoord,dTrans,false);
				}
					}
				for (int curCoord=0; curCoord<3; curCoord++){
						float axisToRot[] = getRotVector(curCoord); //determine the axis of rotation
						float[] myCenter = m.strand[str].getCenterOfMass();
						double[] center = new double[3];
						for(int j=0;j<3;j++)
							center[j] = myCenter[j];
						float dRot = compRot(resNums,resNums.length,rotStep,center,axisToRot,curTransRotInd);
					if (Math.abs(dRot)!=0.0)
							for(int k=0; k<resNums.length;k++)
								updateCumulativeRot(resNums[k],dRot,center,false,axisToRot,false);
				}
			}
			}
			
			//reduce the step size for the next minimization step
			fiPsiStep -= deltaFiPsiStep;
			rotStep -= deltaRotStep;
			transStep -= deltaTransStep;
		}
	}
	
	//Checks the energy if a phi (angleType==0) or psi (angleType==1) change of +fiPsiStep
	//		is applied to residue resNum (strand-relative numbering, only residues in strand strandNum are considered);
	//Checks the energy at -fiPsiStep and at the initial position and compares these three energies
	//		to determine the optimal direction of change;
	//If at least one residue is moved further than the maximum allowed limit for its CA, then the
	//		step in the corresponding direction is not considered
	private float compFiPsi(int resNum, int strandNum, int angleType, float fiPsiStep, Backbone bb, boolean pemComp){
		
		double initialEnergy[], secondEnergy[], thirdEnergy[];
		float storedCoord[][] = new float[m.numberOfResidues][];
		
		//Store the actualCoordinates[] before any changes
		for (int i=0; i<m.numberOfResidues; i++)
			storedCoord[i] = storeCoord(m.residue[i].moleculeResidueNumber);
		
		
		//Check at initial position
		initialEnergy = a96ff.calculateTotalEnergy(m.actualCoordinates,-1);
		
		//Check at +fiPsiStep rotation
		bb.applyFiPsi(m,strandNum,resNum,fiPsiStep,angleType);	
		secondEnergy = a96ff.calculateTotalEnergy(m.actualCoordinates,-1);	
		if ( (!checkCAallRes(strandNum)) ) //at least one residue moved too far away or unallowed steric for PEMcomp, so make this step impossible
			secondEnergy[0] = Math.pow(10, 38);
		for (int i=0; i<m.strand[strandNum].numberOfResidues; i++)
			restoreCoord(m.strand[strandNum].residue[i].moleculeResidueNumber,storedCoord[m.strand[strandNum].residue[i].moleculeResidueNumber]);
		
		//Check at -fiPsiStep rotation
		bb.applyFiPsi(m,strandNum,resNum,-fiPsiStep,angleType);
		thirdEnergy = a96ff.calculateTotalEnergy(m.actualCoordinates,-1);
		if ( (!checkCAallRes(strandNum)) ) //at least one residue moved too far away or unallowed steric for PEMcomp, so make this step impossible
			thirdEnergy[0] = Math.pow(10, 38);
		for (int i=0; i<m.strand[strandNum].numberOfResidues; i++)
			restoreCoord(m.strand[strandNum].residue[i].moleculeResidueNumber,storedCoord[m.strand[strandNum].residue[i].moleculeResidueNumber]);
		
		float step = getDir(initialEnergy[0],secondEnergy[0],thirdEnergy[0],fiPsiStep);
		
		int molResNum = m.strand[strandNum].residue[resNum].moleculeResidueNumber;
		if (step + molCurFiPsiDisp[molResNum][angleType] > sysMaxDihedRot)
			step = sysMaxDihedRot - molCurFiPsiDisp[molResNum][angleType];
		if (step + molCurFiPsiDisp[molResNum][angleType] < -sysMaxDihedRot)
			step = -sysMaxDihedRot - molCurFiPsiDisp[molResNum][angleType];
		
		return step;
	}
	
	//Updates the CA coordinates for all residues; should be called from within functions that use/modify molCurCA[]
	private void updateCurCA(){
		for (int i=0; i<m.numberOfResidues; i++)
			molCurCA[i] = getCAcoord(i);
	}
	
	//Setup the Amber partial arrays
	private void setupPartialAmber(){
		// numFlexRes, flexResAtomList, and flexResListSize include the ligand if one exists
		//if(ligStrNum != -1)
			a96ff.setupPartialArrays(totalFlexRes+totalTransRotStrands,MAX_NUM_ATOMS_RES,flexResAtomList,
				flexResListSize);
		/*else
			a96ff.setupPartialArrays(totalFlexRes,MAX_NUM_ATOMS_RES,flexResAtomList,
				flexResListSize);*/
	}
	
	//Determines the direction for the translation of size transStep for
	//		the numRes number of residues in resNums[] (molecule-relative numbering) in the direction of coord
	private float compTrans(int resNums[], int numRes, int coord, float transStep, int AAnum){
		
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
		initialEnergy = a96ff.calculateTotalEnergy(m.actualCoordinates,AAnum);
		
		//Check at +transStep rotation
		for (int i=0; i<numRes; i++)
			m.translateResidue(resNums[i], d[0], d[1], d[2], false);	
		secondEnergy = a96ff.calculateTotalEnergy(m.actualCoordinates,AAnum);		
		for (int i=0; i<numRes; i++)
			restoreCoord(resNums[i],storedCoord[i]);
		
		//Check at -transStep rotation
		d[coord] = -transStep;
		for (int i=0; i<numRes; i++)
			m.translateResidue(resNums[i], d[0], d[1], d[2], false);		
		thirdEnergy = a96ff.calculateTotalEnergy(m.actualCoordinates,AAnum);		
		for (int i=0; i<numRes; i++)
			restoreCoord(resNums[i],storedCoord[i]);
		
		return getDir(initialEnergy[0],secondEnergy[0],thirdEnergy[0],transStep);
	}
	
	//Checks if a translation of transStep in the direction of coord will move the CA of resNum further than the limit
	//		and applies the optimal translation if (onlyCheck==false)
	private float updateCumulativeTrans(int resNum, int coord, float transStep, boolean onlyCheck){
		
		updateCurCA(); //get the current CA for each residue
		
		//Determine the new CA if we took this step (the translation is only in the direction of the coord coordinate)
		double tmpCA[] = new double[3];
		for (int i=0; i<3; i++)
			tmpCA[i] = molCurCA[resNum][i];
		tmpCA[coord] = molCurCA[resNum][coord] + transStep;
		
		//Compute how large of a step this would be from the start point
		double dist = getDist(tmpCA,molStartCA[resNum]);
			
		if (dist > sysMaxTrans){ // if the step would take us too far away then scale it back
			double curDisp[] = new double[3];
			curDisp[0] = molCurCA[resNum][0] - molStartCA[resNum][0];
			curDisp[1] = molCurCA[resNum][1] - molStartCA[resNum][1];
			curDisp[2] = molCurCA[resNum][2] - molStartCA[resNum][2];
			
			double a[] = new double[2];
			int ind = 0;
			for (int i=0; i<3; i++){
				if (i!=coord)
					a[ind++] = curDisp[i] * curDisp[i];
			}
			
			double mt2 = sysMaxTrans * sysMaxTrans;
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

		if (!onlyCheck) {//actually apply the new translation
			
			m.translateResidue(resNum, theTranslation[0], theTranslation[1], theTranslation[2], false);
			
			// update the current CA
			molCurCA[resNum] = getCAcoord(resNum);
		}
		
		return transStep;
	}
	
	//Determines the direction for the rotation of size rotStep centered at center[], for
	//		the numRes number of residues in resNums[]  (molecule-relative numbering) around the axisToRot axis
	private float compRot(int resNums[], int numRes, float rotStep, double center[], float axisToRot[], int AAnum){		
		
		double initialEnergy[], secondEnergy[], thirdEnergy[];
		float storedCoord[][] = new float[numRes][];
		
		//Store the actualCoordinates for resNum before any changes
		for (int i=0; i<numRes; i++)
			storedCoord[i] = storeCoord(resNums[i]);
		
		
		//Check at initial position
		initialEnergy = a96ff.calculateTotalEnergy(m.actualCoordinates,AAnum);
		
		//Check at +rotStep rotation
		for (int i=0; i<numRes; i++)
			m.rotateResidue(resNums[i], axisToRot[0], axisToRot[1], axisToRot[2], center[0], center[1], center[2], rotStep, false);		
		secondEnergy = a96ff.calculateTotalEnergy(m.actualCoordinates,AAnum);		
		for (int i=0; i<numRes; i++)
			restoreCoord(resNums[i],storedCoord[i]);
		
		//Check at -rotStep rotation
		for (int i=0; i<numRes; i++)
			m.rotateResidue(resNums[i], axisToRot[0], axisToRot[1], axisToRot[2], center[0], center[1], center[2], -rotStep, false);		
		thirdEnergy = a96ff.calculateTotalEnergy(m.actualCoordinates,AAnum);		
		for (int i=0; i<numRes; i++)
			restoreCoord(resNums[i],storedCoord[i]);
		
		return getDir(initialEnergy[0],secondEnergy[0],thirdEnergy[0],rotStep);
	}
	
	//Applies a rotation of rotStep degrees for residue resNum (molecule-relative numbering)
	//		around axis axisToRot[] and center center[];
	//If the rotation moves the CA for the given residue further than the allowed limit, 
	//		the rotation is decresed so that CA is moved to the limit;
	//If centerIsCA is true, then the center of rotation is the CA atom of residue resNum, so the rotation will not change the CA position;
	//Only if (checkOnly==false), then the rotation is actually performed
	private float updateCumulativeRot(int resNum, float rotStep, double center[], boolean checkOnly, float axisToRot[], boolean centerIsCA){
		
		updateCurCA(); //get the current CA for each residue
		
		float storedCoord[] = storeCoord(resNum);
		
		m.rotateResidue(resNum, axisToRot[0], axisToRot[1], axisToRot[2], center[0], center[1], center[2], rotStep, false);
		
		double CAcoord[] = getCAcoord(resNum);
		double dist = getDist(CAcoord,molStartCA[resNum]);
		
		if ( (!centerIsCA) && (dist>sysMaxTrans) ) { //further than allowed, so restore, find the optimal rotStep value, and rotate with that value, if necessary
			
			restoreCoord(resNum,storedCoord);
			
			if ( getDist(molStartCA[resNum],molCurCA[resNum]) > (sysMaxTrans-0.01) ) //the current CA is already approximately at the limit
				return 0.0f;
			
			else {
				CAcoord = getCAcoord(resNum);
				rotStep = binSearchRotStep(CAcoord,axisToRot,center,rotStep,resNum); //find the optimal rotStep value that will rotate the CA to the max limit
				if ((rotStep!=0.0f)&&(!checkOnly)) //apply the rotation
					m.rotateResidue(resNum, axisToRot[0], axisToRot[1], axisToRot[2], center[0], center[1], center[2], rotStep, false);
			}
		}
		
		if (!checkOnly) //rotation is actually performed, so update the CA coordinates
			molCurCA[resNum] = getCAcoord(resNum);
		else
			restoreCoord(resNum,storedCoord);
		
		return rotStep;
	}
	
	//Performs a binary search to find a value for thetaDeg rotation that will rotate the CA for resNum to the maximum allowed limit;
	//The initial coordinates for CA are given in coord[], the rotation is done around axis axisToRot[] and around center center[]
	private float binSearchRotStep(double coord[], float axisToRot[], double center[], float origTheta, int resNum){
		
		double newDist;
		
		double newCoord[] = new double[3];
		for (int i=0; i<3; i++)
			newCoord[i] = coord[i];
		
		float thetaDeg = origTheta;
		float thetaStep = origTheta/2.0f;
		
		int curStep = 1;
		while (Math.abs(thetaStep)>0.01){ //thetaStep can be negative
		
			newCoord = rotatePoint(coord,thetaDeg,axisToRot,center);
			newDist = getDist(newCoord,molStartCA[resNum]);
			if ( (newDist < sysMaxTrans) && (newDist >= (sysMaxTrans-0.01) ) ){
				return thetaDeg;
			}
			else {
				thetaStep = origTheta/(float)Math.pow(2.0,curStep);
				if (newDist>sysMaxTrans)
					thetaDeg -= thetaStep;
				else
					thetaDeg += thetaStep;
			}
			curStep++;
		}
		return 0.0f;
	}
	
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
	
	//Returns the molecule atom numbers of the CA atoms for each residue in the molecule;
	//If a given residue does not have a CA atom, the corresponding entry is set to -1
	private int [] getMolAtNumCA(){
		
		int molAtomNumCA[] = new int[m.numberOfResidues];
		for (int resNum=0; resNum<m.numberOfResidues; resNum++){
			Residue r = m.residue[resNum];
			boolean found = false;
			for (int i=0; i<r.numberOfAtoms; i++){
				if (r.atom[i].name.equalsIgnoreCase("CA")){
					molAtomNumCA[resNum] = r.atom[i].moleculeAtomNumber;
					found = true;
					break;
				}
			}
			if (!found)
				molAtomNumCA[resNum] = -1;
		}
		return molAtomNumCA;
	}
	
	public double getMagnitude(double a[]){
		double sum = 0.0;
		for(int i=0;i<a.length;i++)
			sum += (a[i]*a[i]);
		return(Math.sqrt(sum));
	}
	
	public float getMaxDistCA(){
		return sysMaxTrans;
	}
	
	public float getMaxDihedRot(){
		return sysMaxDihedRot;
	}
	
	//Returns the CA coordinates for residue resNum (molecule-relative numbering);
	//If the given residue does not have a CA atom, then the coordinates of the geometric center of the residue are returned
	public double [] getCAcoord(int resNum){
		
		double coord[] = new double[3];
		
		if (molAtNumCA[resNum]>=0) {//use the CA atom
			int curAtom = molAtNumCA[resNum];
			coord[0] = m.actualCoordinates[curAtom*3];
			coord[1] = m.actualCoordinates[curAtom*3+1];
			coord[2] = m.actualCoordinates[curAtom*3+2];
		}
		else { //no CA atom, so compute and return the geometric center of the residue
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
	
	//Returns the distance between the two points given by coordinates coord1[] and coord2[]
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
	
	//Checks the distance between the current CA position and the initial CA position for all residues;
	//		If at least one residue has moved too far, return false; otherwise, return true
	public boolean checkCAallRes(int strandNum){
		
		for (int i=0; i<m.strand[strandNum].numberOfResidues; i++){
			int resNum = m.strand[strandNum].residue[i].moleculeResidueNumber;
			double tmpCA[] = getCAcoord(resNum);
			double dist = getDist(tmpCA,molStartCA[resNum]);
			if (dist>sysMaxTrans)
				return false;
		}
		return true;
	}
}
