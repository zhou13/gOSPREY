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
//	StericCheck.java
//
//	Version:           2.1 beta
//
//
//	  authors:
// 	  initials    name                 organization                email
//	 ---------   -----------------    ------------------------    ----------------------------
//	  ISG		 Ivelin Georgiev	  Duke University			  ivelin.georgiev@duke.edu
//	  KER        Kyle E. Roberts       Duke University         ker17@duke.edu
//    PGC        Pablo Gainza C.       Duke University         pablo.gainza@duke.edu
///////////////////////////////////////////////////////////////////////////////////////////////

/**
 * Written by Ivelin Georgiev (2004-2009)
 * 
 */

import java.io.BufferedOutputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Hashtable;
import java.math.*;

/**
 * Implements functions for checking the steric overlap for specified parts of a given structure at different
 * stages of the A* expansion of the conformation tree.
 *
 */
public class StericCheck {	
	
	//numAS residues
	//int numInAS = -1;
	int numMutable = -1;
	//the cur AA at each residue
	int curAANum[] = null;
	
	//the cur AA for the ligand
	//int ligAANum = -1;
	
	//the mapping from residue number to AS num and from AS num to res num
	//int curResToASMap[] = null;
	//int residueMap[] = null;
	int strandMut[][] = null;
	int numberOfStrands = 0;
	//the current molecule
	Molecule m = null;
	
	//the system and ligand strand numbers
	//int sysStrNum = -1;
	//int ligStrNum = -1;
	
	//the rotamer libraries
	//RotamerLibrary rl = null; //for the protein
	//RotamerLibrary grl = null; //for the ligand
	
	//the sys and lig rotamer handlers
	//StrandRotamers sysLR = null;
	//StrandRotamers ligROT = null;
	StrandRotamers strandRot[] = null;
	//the overlap threshold for the steric check
	double overlapThresh = 1.5;
	
	//determines whether hydrogens are used in the steric check
	boolean hSteric = false;
	
	//the MinDEE pruning matrix
	boolean eliminatedRot[] = null;
	
	//num rotamers for each residue
	int numRotForRes[] = null;

	//the number of conformations not examined yet
	BigInteger numConfsLeft = new BigInteger("0");
	
	//the number of conformations above each level (from level i+1 to the top level,
	//	which is the ligand)
	BigInteger numConfsAboveLevel[] = null;
	
	//the number of conformations pruned by the steric filter
	BigInteger numConfsPrunedByS = new BigInteger("0");
	
	int mutRes2Strand[] = null;
	int mutRes2StrandMutIndex[] = null;

        //Indicates use of DEEPer
        boolean doPerturbations = false;

	//determines whther a ligand is present
	//boolean ligPresent = false;
	
	//PrintStream logPS = null;
	
	/*StericCheck (int curAANumP[],int curResToASMapP[],int strandMutP[],boolean eliminatedRotAtPosRedP[],
			int numRotForResP[],Molecule mP, double overlapThreshP,	boolean hS, BigInteger numConfsLeftP,
			BigInteger numConfsAboveLevelP[], int numMutableP, StrandRotamers sysLRP, int ligStrNumP, 
			StrandRotamers ligROTP, int curLigAANum, RotamerLibrary rlP, RotamerLibrary grlP){
		
		curAANum = curAANumP;
		//curResToASMap = curResToASMapP;
		//residueMap = residueMapP;
		//numInAS = residueMap.length;
		strandMut = strandMutP;
		m = mP;
		overlapThresh = overlapThreshP;
		hSteric = hS;
		eliminatedRot = eliminatedRotAtPosRedP;
		numRotForRes = numRotForResP;
		numConfsLeft = numConfsLeftP;
		numConfsAboveLevel = numConfsAboveLevelP;
		numConfsPrunedByS = new BigInteger("0");
		sysStrNum = sysStrNumP;
		sysLR = sysLRP;
		ligStrNum = ligStrNumP;
		ligROT = ligROTP;
		ligAANum = curLigAANum;
		ligPresent = true;
		rl = rlP;
		grl = grlP;
		
		//logPS = logPSP;
	}*/
	
	StericCheck (int curAANumP[],/*int curResToASMapP[],*/int strandMutP[][],boolean eliminatedRotAtPosRedP[],
			int numRotForResP[],Molecule mP, double overlapThreshP, boolean hS, BigInteger numConfsLeftP, 
			BigInteger numConfsAboveLevelP[], int numMutableP, int numberOfStrandsP, StrandRotamers strandRotP[],
			int mutRes2StrandP[], int mutRes2StrandMutIndexP[], boolean doPerts){
		
		curAANum = curAANumP;
		//curResToASMap = curResToASMapP;
		//residueMap = residueMapP;
		//numInAS = residueMap.length;
		strandMut = strandMutP;
		numMutable = numMutableP;
		
		m = mP;
		overlapThresh = overlapThreshP;
		hSteric = hS;
		eliminatedRot = eliminatedRotAtPosRedP;
		numRotForRes = numRotForResP;
		numConfsLeft = numConfsLeftP;
		numConfsAboveLevel = numConfsAboveLevelP;
		numConfsPrunedByS = new BigInteger("0");
		//sysStrNum = sysStrNumP;
		//sysLR = sysLRP;
		//ligPresent = false;
		numberOfStrands = numberOfStrandsP;
		strandRot = strandRotP;
		mutRes2Strand = mutRes2StrandP;
		mutRes2StrandMutIndex = mutRes2StrandMutIndexP;
                doPerturbations = doPerts;

		
		//logPS = logPSP;
	}
	
	//Return the number of conformations not examined yet
	public BigInteger getNumConfsLeft() {
		return numConfsLeft;
	}
	
	//Return the number of conformations pruned by the steric filter
	public BigInteger getNumConfsPrunedByS(){
		return numConfsPrunedByS;
	}
	
	//Sets the number of conformations remaining to be examined
	public void setNumConfsLeft(BigInteger numCLeft){
		numConfsLeft = numCLeft;
	}
	
	//Sets up the steric check for the given partial (full) conformation: curConf[] has been
	//	assigned for levels 0...(curTopLevel-1) and curNode should be applied at curTopLevel;
	//Returns true if the partial conformation is sterically allowed
	//Called by A*
	public boolean checkAllowedSteric (int curTopLevel, int curConf[], int curNode){		
		
		//As the rotamers given to A* are only the non-pruned ones, there is a difference between the
		//	rotamer numbers returned by A* and the actual rotamer numbers for each residue (that is,
		//	A* may return rot 4 for res 3, but rot 3 for res 3 may be pruned, and so the actual number
		//	of the rot to be applied for res 3 is 5)
		int curPruningInd = 0;
		int curRotInd, compInd;
		int conf[] = new int[curTopLevel+1];
		for (int curRes=0; curRes<=curTopLevel; curRes++){//compute the actual rot numbers for levels 0...curTopLevel
			curRotInd = 0;
			for (int curRot=0; curRot<numRotForRes[curRes]; curRot++){
				if (!eliminatedRot[curPruningInd]){
					if (curRes==curTopLevel) //the rotamer (curNode) for curTopLevel is not a part of curConf[]
						compInd = curNode;
					else
						compInd = curConf[curRes];
					
					if (curRotInd==compInd)
						conf[curRes] = curRot;
					curRotInd++;
				}
				curPruningInd++;
			}
		}

		//Backup the atom coordinates, so that they can be restored after the steric check, as 
		//		applyRotamer() changes both the actualCoordinates[] and the atom coordinates,
		//		so we cannot restore the original position just using m.updateCoordinates()
		m.backupAtomCoord();
		
		//Apply the rotamers of the current partial conformation (up to level (curTopLevel-1))
		int curAS = 0;
		boolean applyLig = false;
		
		//If there is a ligand and curTopLevel is the ligand level (full conformation), apply the lig rotamer;
		//	otherwise, curTopLevel is an AS residue
		/*if ((ligPresent)&&(curTopLevel==numInAS)){ //apply the ligand rotamer
			if (grl.getNumRotForAAtype(ligAANum)!=0){//not GLY or ALA
				ligROT.applyRotamer(m, 0, conf[curTopLevel]);//the ligand level
			}
			applyLig = true;
		}*/
		for(int i=0; i<=curTopLevel;i++){
				//int curL = oldLevelToNew[i];
				int str = mutRes2Strand[i];
				int strResNum = strandMut[str][mutRes2StrandMutIndex[i]];
		//for (int str=0; str<numberOfStrands;str++){
		//	for(int i=0;i<strandMut[str].length;i++){
		//		if(curAS>curTopLevel)
		//			break;
				int molResNum = m.strand[str].residue[strResNum].moleculeResidueNumber;

                                if(doPerturbations)
                                        ((StrandRCs)strandRot[str]).applyRC(m, strResNum, conf[curAS]);
                                else if (strandRot[str].rl.getNumRotForAAtype(curAANum[molResNum])!=0){//not GLY or ALA
					strandRot[str].applyRotamer(m, strResNum, conf[curAS]);
                                }
		//		curAS++; //prepare the next AS residue
		//	}
		}
		
		/*for (int curRes=0; curRes<m.strand[sysStrNum].numberOfResidues; curRes++){
			if (curAS<curTopLevel){ //apply for AS res 0...(curTopLevel-1)
				if (curResToASMap[curRes]!=-1){//make a change only to the AS residues: use the native type for the other residues
										
					if (rl.getNumRotForAAtype(curAANum[residueMap[curAS]])!=0){//not GLY or ALA
						sysLR.applyRotamer(m, curRes, conf[curAS]);
					}
					curAS++; //prepare the next AS residue
				}
			}
			else if (!applyLig) { //we need to apply an AS rot at curTopLevel, as there is no ligand
				if (curResToASMap[curRes]!=-1){
					if (rl.getNumRotForAAtype(curAANum[residueMap[curAS]])!=0)//not GLY or ALA
						sysLR.applyRotamer(m, curRes, conf[curTopLevel]);
					break;
				}
			}
			else //we have already applied all of the rotamers for the given partial conformation
				break;
		}*/
		
		/*logPS.println("curTopLevel "+curTopLevel+" curNode "+curNode+" curConf ");
		for (int i=0;i<=curTopLevel;i++)logPS.print(conf[i]+" ");logPS.println();
		if (!ligPresent){
			for (int i=0;i<=curTopLevel;i++)logPS.print(sysLR.getCurRotNum(residueMap[i])+" ");logPS.println();logPS.flush();
		}
		else{
			for (int i=0;i<curTopLevel;i++)logPS.print(sysLR.getCurRotNum(residueMap[i])+" ");
			logPS.print(ligROT.getCurRotNum(0));logPS.println();logPS.flush();
		}*/
		
		boolean allowedSteric = true;
		//int curTL = oldLevelToNew[curTopLevel];
		int str = mutRes2Strand[curTopLevel];
		int strResNum = strandMut[str][mutRes2StrandMutIndex[curTopLevel]];
		//Do the steric checks
		/*if ((ligPresent)&&(curTopLevel==numInAS)) //check the ligand (which is at the top level) against all other residues
			allowedSteric = RS_CheckAllSterics(ligStrNum,0);
		else*/ 
			allowedSteric = RS_CheckAllSterics(str,strResNum);
		
		m.restoreAtomCoord(); //restore the atom coordinates
		m.updateCoordinates(); //restore the actualCoordinates
		
		if (!allowedSteric){ //decrease the number of remaining conformations
			numConfsLeft = numConfsLeft.subtract(numConfsAboveLevel[curTopLevel]);
			numConfsPrunedByS = numConfsPrunedByS.add(numConfsAboveLevel[curTopLevel]);
		}
		//logPS.println("allowedSteric "+allowedSteric+" confsLeft "+numConfsLeft+" confsPrunedByS "+numConfsPrunedByS);
		//logPS.println();logPS.flush();
		
		return allowedSteric;
	}
	
	// This function is similar to RS_CheckSterics
	// This version checks all residues against the target residue rather
	//  than just checking residues up to the specified one in the strand
	// The AS residues with res numbers greater than the resNum for the
	//	target residue are not checked, as they have not been assigned yet;
	// Also, when resNum is a residue in the AS, the ligand is not checked,
	//	as it has not been assigned yet
	private boolean RS_CheckAllSterics(int strandNum, int resNum) {
		
		ProbeStericCheck psc = new ProbeStericCheck();
		
		//The mapping from AS res to the system strand residue numbering
		boolean curASToResMap[][] = new boolean[numberOfStrands][];
		
		for (int str=0;str<numberOfStrands;str++){
			curASToResMap[str] = new boolean[m.strand[str].numberOfResidues];
		}
	
		for(int str=0;str<numberOfStrands;str++){
			for (int i=0; i<curASToResMap[str].length; i++)
				curASToResMap[str][i] = false;
		}
		for(int str=0;str<numberOfStrands;str++){
			for (int i=0; i<strandMut[str].length; i++)
				curASToResMap[str][strandMut[str][i]] = true;
	
		}		
		
		Residue res = m.strand[strandNum].residue[resNum];
	
		// It's important to divide below by 100.0 rather than
		//  just 100 as java assumes precision or does some
		//  strange rounding that causes an incorrect results
		//  when 100 is used
		
		for(int i=0;i<res.numberOfAtoms;i++) {
			Atom a1 = res.atom[i];
			if ( hSteric || (!a1.elementType.equalsIgnoreCase("H")) ) {
				for(int q=0;q<m.numberOfStrands;q++) {
					int resToCheck = m.strand[q].numberOfResidues;
					for(int w=0;w<resToCheck;w++) {
						if(!((q==strandNum) && (w==resNum))) {//not the same residue
							if( !( ((q>strandNum)&&curASToResMap[q][w]) || ((q==strandNum)&&(w>resNum)&&(curASToResMap[q][w])) ) ){
							//if ((!((strandNum==sysStrNum)&&(q==strandNum)&&(w>resNum)&&(curASToResMap[w])))&&
									//(!((q==ligStrNum)&&(strandNum==sysStrNum)))){//not an AS residue with a bigger res number AND not the ligand
								for(int t=0;t<m.strand[q].residue[w].numberOfAtoms;t++) {
									Atom a2 = m.strand[q].residue[w].atom[t];
									if ( hSteric || (!a2.elementType.equalsIgnoreCase("H")) ) {
										if (!psc.isAllowedSteric(m, a1, a2, (float)overlapThresh))
											return false;
									}
								}
							}
							else if ( ((q>strandNum)&&curASToResMap[q][w]) ||((q==strandNum)&&(w>resNum)&&(curASToResMap[q][w])) ) {
                                                            if(!doPerturbations){//Don't do this check in DEEPer because the backbone atoms may move still
							//else if ((strandNum==sysStrNum)&&(q==strandNum)&&(w>resNum)&&(curASToResMap[w])){ //an unassigned AS residue, so check only the backbone atoms of that residue
								for(int t=0;t<m.strand[q].residue[w].numberOfAtoms;t++) {
									Atom a2 = m.strand[q].residue[w].atom[t];
									if (a2.getIsBBatom()){
										if ( hSteric || (!a2.elementType.equalsIgnoreCase("H")) ) {
											if (!psc.isAllowedSteric(m, a1, a2, (float)overlapThresh))
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
		}
	
		// If you got here then everything passed
		return true;
	}
}
