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
//	DEEGoldsteinPairs.java
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
//     MAH        Mark A. Hallen	Duke University         mah43@duke.edu
///////////////////////////////////////////////////////////////////////////////////////////////

/**
* Written by Ivelin Georgiev (2004-2009)
* 
*/

/**
 * Performs DEE Goldstein pairs rotamer pruning
 * 
*/
public class DEEGoldsteinPairs extends DEE {
    
	
	//the single and pair interval terms in the DE pairs MinDEE criterion
	double indIntMinDEE2Pos[][] = null;
	double pairIntMinDEE2Pos[][] = null;
	
		
		
	//the max scaling factor for the interval terms
	float maxScale = 1.0f;
	

	//constructor
	DEEGoldsteinPairs(PairwiseEnergyMatrix arpMatrix, PairwiseEnergyMatrix arpMatrixMax, int numResMutable,
			int strMut[][], float initEw, 
			StrandRotamers strandLRot[], PrunedRotamers<Boolean> prunedRotAtRes, boolean residueMut[],
			boolean doMin, boolean spFlags[][][][][][], boolean useSF, boolean mb, boolean dDEE, boolean minBB, 
                        boolean scaleInt, float maxSc, int mutRes2StrandP[], int mutRes2MutIndexP[], boolean typeDep, boolean aIMinDEE,
                        float aIval, boolean tripFlags[][][][][][][][][], boolean doPerts) {


                init(arpMatrix, arpMatrixMax, numResMutable,
			strMut, initEw, strandLRot, prunedRotAtRes, doMin, null, null,
                        spFlags, useSF, minBB, mutRes2StrandP, mutRes2MutIndexP, typeDep, aIMinDEE, aIval,
                        mb, dDEE, residueMut, tripFlags, doPerts);

		
		if (!scaleInt) //no scaling of the interval terms performed
			maxScale = 1.0f;
		else
			maxScale = maxSc;
	}
	

	//Compute the conformations that can be eliminated
	//Return a boolean matrix in which an element is true if
	//the corresponding rotamer pair can be eliminated, and false otherwise
	public void ComputeEliminatedRotConf(){
		
		// 2010: No interval terms for useMinDEEPruningEw
		if (doMinimize && !doIMinDEE){ //compute the DE pairs MinDEE interval terms
			
			int numRes = numMutable;		
			/*if (numLigRot!=0) //ligand is present
				numRes++;*/
			
			indIntMinDEE2Pos = new double[numRes][numRes];
			pairIntMinDEE2Pos = new double[numRes][numRes];
			
			for (int posNum1=0; posNum1<numRes; posNum1++){
				for (int posNum2=posNum1+1; posNum2<numRes; posNum2++){
		
					indIntMinDEE2Pos[posNum1][posNum2] = SumMaxIndInt(posNum1,posNum2);				//formula term 3
					pairIntMinDEE2Pos[posNum1][posNum2] = SumSumMaxPairInt(posNum1,posNum2);		//formula term 4
					indIntMinDEE2Pos[posNum2][posNum1] = indIntMinDEE2Pos[posNum1][posNum2];		//formula term 3
					pairIntMinDEE2Pos[posNum2][posNum1] = pairIntMinDEE2Pos[posNum1][posNum2];		//formula term 4
				}
			}
		}
			
		//Check for pairs pruning
		int numRotForCurAAatPos1;
		
		int prunedCurRun = 0;
		boolean done = false;
		numRuns = 1;
		
		while (!done){
			
			prunedCurRun = 0;
			
			System.out.println("Current run: "+numRuns);
		
			//Compute for the AS residues first
			for (int curPos1=0; curPos1<numMutable; curPos1++){
				System.out.println("");
				System.out.print("Pos1 "+curPos1+": ");
				
				int str1=mutRes2Strand[curPos1];
				int strResNum1=strandMut[str1][mutRes2MutIndex[curPos1]];
				if ((magicBullet)||(!distrDEE)||(resInPair[curPos1])){ //mb-pairs or not distrDEE or cur res is in distr pair
				
					//System.out.print("Starting AS residue "+curPos1);
					//System.out.print("..");
					
					for (int curPos2=curPos1+1; curPos2<numMutable; curPos2++){
						System.out.print(curPos2+" ");
						int str2=mutRes2Strand[curPos2];
						int strResNum2=strandMut[str2][mutRes2MutIndex[curPos2]];
						if ((magicBullet)||(!distrDEE)||(resInPair[curPos2])){ //mb-pairs or not distrDEE or cur res is in distr pair
							
							for (int AA1=0; AA1<numAAtypes[curPos1]; AA1++){
								
								int curAA1 = strandRot[str1].getIndexOfNthAllowable(strResNum1,AA1);
								
								//find how many rotamers are allowed for the current AA type at the given residue;
								//note that ala and gly have 0 possible rotamers
								numRotForCurAAatPos1 = getNumRot(str1, strResNum1, curAA1);
								
								for(int curRot1=0; curRot1<numRotForCurAAatPos1; curRot1++){
									
									//int index_1 = curPos1*numTotalRot + rotIndOffset[curAA1] + curRot1;
								
									if (!eliminatedRotAtPos.get(curPos1,curAA1,curRot1)){//not already pruned
										
										for (int AA2=0; AA2<numAAtypes[curPos2]; AA2++){
											int curAA2 = strandRot[str2].getIndexOfNthAllowable(strResNum2,AA2);
											int numRotForCurAAatPos2 = getNumRot( str2, strResNum2, curAA2);
											
											for(int curRot2=0; curRot2<numRotForCurAAatPos2; curRot2++){
											
												//int index_2 = curPos2*numTotalRot + rotIndOffset[curAA2] + curRot2;
												
												if (!eliminatedRotAtPos.get(curPos2,curAA2,curRot2)){//not already pruned
												
													if (!splitFlags[curPos1][curAA1][curRot1][curPos2][curAA2][curRot2]){ //rotamer pair not already pruned
											
														if (CanEliminate(curPos1, curAA1, curRot1, curPos2, curAA2, curRot2)){
															splitFlags[curPos1][curAA1][curRot1][curPos2][curAA2][curRot2] = true;
															splitFlags[curPos2][curAA2][curRot2][curPos1][curAA1][curRot1] = true;
															
															prunedCurRun++;
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
					//System.out.println("done");
				}
			}
			
			//If there is a ligand, compute DEE for the lig rotamers as well
			/*if (numLigRot!=0){
				
				if ((magicBullet)||(!distrDEE)||(resInPair[numSiteResidues])){ //mb-pairs or not distrDEE or ligand is in distr pair
					
					System.out.print("Starting ligand run");
					System.out.print("..");
					
					for (int curPos2=0; curPos2<numSiteResidues; curPos2++){ //curPos2 is AS residue, always != ligand
						
						if ((magicBullet)||(!distrDEE)||(resInPair[curPos2])){ //mb-pairs or not distrDEE or cur res is in distr pair
							
							for (int curRotLig=0; curRotLig<numLigRot; curRotLig++){
								int indexLig = numSiteResidues*numTotalRot+curRotLig;
								
								if (!eliminatedRotAtPos[indexLig]){//not already pruned
									
									for (int AA2=0; AA2<numAAtypes[curPos2]; AA2++){
										int curAA2 = sysLR.getIndexOfNthAllowable(residueMap[curPos2],AA2);
										int numRotForCurAAatPos2 = rl.getNumRotForAAtype(curAA2);
										if (numRotForCurAAatPos2==0)	//ala or gly
											numRotForCurAAatPos2 = 1;
										
										for(int curRot2=0; curRot2<numRotForCurAAatPos2; curRot2++){
										
											int index_2 = curPos2*numTotalRot + rotIndOffset[curAA2] + curRot2;
											
											if (!eliminatedRotAtPos[index_2]){//not already pruned
											
												if (!splitFlags[indexLig][index_2]){ //rotamer pair not already pruned
								
													if (CanEliminateLig(curRotLig,curPos2,curAA2,curRot2)){
														splitFlags[indexLig][index_2] = true;
														splitFlags[index_2][indexLig] = true;
								
														prunedCurRun++;
													}
												}
											}
										}
									}
								}
							}
						}
					}
					System.out.println("done");
				}
			}*/
			
			System.out.println("Number of pairs pruned this run: "+prunedCurRun);
			System.out.println("DEE: The minimum difference is "+minDiff);
			System.out.println();
			
			if (prunedCurRun==0) //no rotamers pruned this run, so done
				done = true;
			else 
				numRuns++;
			
			if ((!magicBullet)&&(!distrDEE)) //full non-distributed pairs, so do not repeat (computationally-expensive)
				done = true;
		}
	}
	
	//Called only by ComputeEliminatedRotConf(.)
	private boolean CanEliminate (int posNum1, int AANumAtPos1, int rotNumAtPos1, int posNum2, 
			int AANumAtPos2, int rotNumAtPos2){
		
		double minIndVoxelE, maxIndVoxelE;
		double minPairE, maxPairE;
		double indVoxelInterval, pairVoxelInterval;
		double minDiffPairVoxelE;
				
		double checkSum;
		
		int str1=mutRes2Strand[posNum1];
		int strResNum1=strandMut[str1][mutRes2MutIndex[posNum1]];
		int str2=mutRes2Strand[posNum2];
		int strResNum2=strandMut[str2][mutRes2MutIndex[posNum2]];
		
		//In the energy matrix, column 0 gives the individual energies for each r at i;
		//skip row 0, as the individual energies start from row 1 (and are all in column 0)

		minIndVoxelE = pairwiseMinEnergyMatrix.getIntraAndShellE( posNum1, AANumAtPos1, rotNumAtPos1 )
                        + pairwiseMinEnergyMatrix.getIntraAndShellE( posNum2, AANumAtPos2, rotNumAtPos2 );//formula term 1.  Intra+shell energy
		minPairE = pairwiseMinEnergyMatrix.getPairwiseE( posNum1, AANumAtPos1, rotNumAtPos1, posNum2, AANumAtPos2, rotNumAtPos2 );
		
		
		if ((!eliminatedRotAtPos.get(posNum1,AANumAtPos1,rotNumAtPos1))&&(!eliminatedRotAtPos.get(posNum2,AANumAtPos2,rotNumAtPos2))
				&&(!splitFlags[posNum1][AANumAtPos1][rotNumAtPos1][posNum2][AANumAtPos2][rotNumAtPos2])){ //not pruned
		
			// 2010: Don't compute the intervals if useMinDEEPruningEw is true since they are not neces
			if (doMinimize && !doIMinDEE){ //MinDEE, so compute the interval terms
				indVoxelInterval = indIntMinDEE2Pos[posNum1][posNum2];							//formula term 3
				pairVoxelInterval = pairIntMinDEE2Pos[posNum1][posNum2];						//formula term 4
			}
			else { //traditional-DEE, so no interval terms // 2010: if useMinDEEPruningEw is set to true then it is the same as DEE
				indVoxelInterval = 0.0;
				pairVoxelInterval = 0.0;
			}
			
			//For the particular position, compare the energy performance (one by one)
			//of the remaining rotamer possibilities to that of the given rotamer:
			//given r at i, compare it to all t at i for pruning
			int numRotForAAatPos1;
			
			for (int AA1=0; AA1<numAAtypes[posNum1]; AA1++){
				
				int altAA1 = strandRot[str1].getIndexOfNthAllowable(strResNum1,AA1);
				
				numRotForAAatPos1 = getNumRot( str1, strResNum1, altAA1 );
				
				for (int altRot1=0; altRot1<numRotForAAatPos1; altRot1++){
							
					//if t and r are not actually the same rotamer of the same AA
					if (!((altAA1==AANumAtPos1)&&(altRot1==rotNumAtPos1))){
						
						//index_t1 = posNum1*numTotalRot + rotIndOffset[altAA1] + altRot1;
						
						if ((!eliminatedRotAtPos.get(posNum1,altAA1,altRot1))){ //not pruned
							
							int numRotForAAatPos2;
							
							for (int AA2=0; AA2<numAAtypes[posNum2]; AA2++){
								
								int altAA2 = strandRot[str2].getIndexOfNthAllowable(strResNum2,AA2);
								if(!typeDependent || (altAA1 == AANumAtPos1 && altAA2 == AANumAtPos2)){

                                                                    numRotForAAatPos2 = getNumRot( str2, strResNum2, altAA2 );
								
                                                                    for (int altRot2=0; altRot2<numRotForAAatPos2; altRot2++){
											
									//if t and r are not actually the same rotamer of the same AA
									if (!((altAA2==AANumAtPos2)&&(altRot2==rotNumAtPos2))){
										
											//index_t2 = posNum2*numTotalRot + rotIndOffset[altAA2] + altRot2;
										
											if ((!eliminatedRotAtPos.get(posNum2,altAA2,altRot2))){ //not pruned 


                                                                                            maxIndVoxelE = pairwiseMaxEnergyMatrix.getIntraAndShellE( posNum1, altAA1, altRot1 )
                                                                                                    + pairwiseMaxEnergyMatrix.getIntraAndShellE( posNum2, altAA2, altRot2 );//formula term 1.  Intra+shell energy
                                                                                            maxPairE = pairwiseMaxEnergyMatrix.getPairwiseE( posNum1, altAA1, altRot1, posNum2, altAA2, altRot2 );


                                                                                            minDiffPairVoxelE = SumMinDiffPVE(posNum1, AANumAtPos1, rotNumAtPos1, altAA1, altRot1,
                                                                                                            posNum2, AANumAtPos2, rotNumAtPos2, altAA2, altRot2);	//formula term 5

                                                                                            checkSum = -templateInt + (minIndVoxelE + minPairE) - (maxIndVoxelE + maxPairE)
                                                                                                                    - indVoxelInterval - pairVoxelInterval + minDiffPairVoxelE;

                                                                                            if (checkSum > curEw){
                                                                                                    //System.out.println(index_r+" "+index_t+" "+checkSum+" "+minIndVoxelE+" "+minShellResE+" "+maxIndVoxelE+" "+maxShellResE+" "+indVoxelInterval+" "+pairVoxelInterval+" "+minDiffPairVoxelE);

                                                                                                    return true;}//this rotamer can be pruned/eliminated
                                                                                            else {
                                                                                                    minDiff = Math.max(minDiff,checkSum);

                                                                                                    if (magicBullet) //magic bullet pairs, so no further checks
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
		}
		
		//We have tried all of the other rotamers at the current position and none
		//of them is able to prune the given rotamer, so we return false
		return false;
	}
	
	////////////////////////////////////////////////////////////////////////

	//////////////////////////////////////////////////////////////////////////////////
	
	//Called only by CanEliminate(.)
	private double SumMinDiffPVE (int atPos1, int withAA1, int withRot1, int altAA1, int altRot1,
			int atPos2, int withAA2, int withRot2, int altAA2, int altRot2){
		
		double sum = 0;
		
		//get the contribution from the active site residue rotamers
		for (int curPos=0; curPos<numMutable; curPos++){			
				
			if ((curPos != atPos1)&&(curPos!=atPos2))
			
				sum += IndMinDiffPVE(atPos1, withAA1, withRot1, altAA1, altRot1, 
						atPos2, withAA2, withRot2, altAA2, altRot2, curPos);
		}
		
		/*if (numLigRot!=0){ //there is a ligand
			//get the contribution from the ligand rotamers: there is only one ligand residue,
			//so there is only one position j here for which to add
			sum += LigandIndMinDiffPVE (atPos1, withAA1, withRot1, altAA1, altRot1, atPos2, withAA2, withRot2, altAA2, altRot2);
		}*/

		return sum;
	}
	
	//Called by SumMaxMaxPVE(.)
	private double IndMinDiffPVE (int firstPos, int firstAA1, int firstRot1, int firstAltAA1, int firstAltRot1, 
			int secondPos, int secondAA1, int secondRot1, int secondAltAA1, int secondAltRot1, int thirdPos){
		
		double minE = bigE;
		double curEmin, curEmax;
		
		Index3 index_r1, index_r2, index_t1, index_t2, index2;
		int numRotForAAatPos;
		
		//r at i
		index_r1 = new Index3(firstPos,firstAA1,firstRot1);//firstPos*numTotalRot + rotIndOffset[firstAA1] + firstRot1;
		index_r2 = new Index3(secondPos,secondAA1,secondRot1);//secondPos*numTotalRot + rotIndOffset[secondAA1] + secondRot1;
		
		//t at i
		index_t1 = new Index3(firstPos,firstAltAA1,firstAltRot1);//firstPos*numTotalRot + rotIndOffset[firstAltAA1] + firstAltRot1;
		index_t2 = new Index3(secondPos,secondAltAA1,secondAltRot1);//secondPos*numTotalRot + rotIndOffset[secondAltAA1] + secondAltRot1;
		
		int str3=mutRes2Strand[thirdPos];
		int strResNum3=strandMut[str3][mutRes2MutIndex[thirdPos]];
		
		
		boolean found = false;
		
		if (((!eliminatedRotAtPos.get(index_r1))&&(!eliminatedRotAtPos.get(index_r2)))&&
				((!eliminatedRotAtPos.get(index_t1))&&(!eliminatedRotAtPos.get(index_t2)))){ //not pruned
		
			for (int AA=0; AA<numAAtypes[thirdPos]; AA++){
				
				int curAA = strandRot[str3].getIndexOfNthAllowable(strResNum3,AA);
				
				numRotForAAatPos = getNumRot( str3, strResNum3, curAA );
				
				for (int curRot=0; curRot<numRotForAAatPos; curRot++){
						
					//There is a displacement: column 0 and row 0 have special entries, 
					//so pairwise energies start from row 1, column 1
					
					//s at j
					index2 = new Index3(thirdPos,curAA,curRot);//thirdPos*numTotalRot + rotIndOffset[curAA] + curRot;
					
					if ((!eliminatedRotAtPos.get(index2))){ //not pruned 
						
						if ( (!useFlags) || (
                                                        (!splitFlags[firstPos][firstAA1][firstRot1][thirdPos][curAA][curRot])
                                                        &&(!splitFlags[secondPos][secondAA1][secondRot1][thirdPos][curAA][curRot])
                                                        && (!isPrunedTriple(firstPos,firstAA1,firstRot1,secondPos,secondAA1,secondRot1,thirdPos,curAA,curRot))
                                                        )){ //not using split flags or not flagged
		
							curEmin = pairwiseMinEnergyMatrix.getPairwiseE( firstPos, firstAA1, firstRot1, thirdPos, curAA, curRot ) + pairwiseMinEnergyMatrix.getPairwiseE( secondPos, secondAA1, secondRot1, thirdPos, curAA, curRot );
							curEmax = pairwiseMaxEnergyMatrix.getPairwiseE( firstPos, firstAltAA1, firstAltRot1, thirdPos, curAA, curRot ) + pairwiseMaxEnergyMatrix.getPairwiseE( secondPos, secondAltAA1, secondAltRot1, thirdPos, curAA, curRot );
							
							if ((curEmin-curEmax) < minE)
								minE = curEmin-curEmax;
							
							found = true;
						}
					}					
				}
			}
			
			//if(minE==bigE)//make it contribute nothing to the sum
			//	minE = 0;
		}
		
		if (!found) //no possible pairs found
			minE = 0.0; //contributes nothing to the sum
		
		return minE;
	}
	
	//Called by SumMaxMaxPVE(.)
	/*private double LigandIndMinDiffPVE (int firstPos, int firstAA1, int firstRot1, int firstAltAA1, int firstAltRot1, 
			int secondPos, int secondAA1, int secondRot1, int secondAltAA1, int secondAltRot1){
		
		double minE = bigE;
		double curEmin, curEmax;
		
		int index_r1, index_r2, index_t1, index_t2, index2;
		
		//r at i
		index_r1 = firstPos*numTotalRot + rotIndOffset[firstAA1] + firstRot1;
		index_r2 = secondPos*numTotalRot + rotIndOffset[secondAA1] + secondRot1;
		
		//t at i
		index_t1 = firstPos*numTotalRot + rotIndOffset[firstAltAA1] + firstAltRot1;
		index_t2 = secondPos*numTotalRot + rotIndOffset[secondAltAA1] + secondAltRot1;
		
		boolean found = false;
		
		if (((!eliminatedRotAtPos[index_r1])&&(!eliminatedRotAtPos[index_r2]))&&
				((!eliminatedRotAtPos[index_t1])&&(!eliminatedRotAtPos[index_t2]))){ //not pruned
			
			for (int curLigPos=0; curLigPos<numLigRot; curLigPos++){
				
				//s at j (the ligand residue)
				index2 = numSiteResidues*numTotalRot + curLigPos;
				
				if ((!eliminatedRotAtPos[index2])){ //not pruned 
					
					if ((!useFlags)||((!splitFlags[index_r1][index2])&&(!splitFlags[index_r2][index2]))){ //not using split flags or not flagged
				
						curEmin = pairwiseMinEnergyMatrix[firstPos][firstAA1][firstRot1][numSiteResidues][ligAANum][curLigPos] + pairwiseMinEnergyMatrix[secondPos][secondAA1][secondRot1][numSiteResidues][ligAANum][curLigPos];
						curEmax = pairwiseMaxEnergyMatrix[firstPos][firstAltAA1][firstAltRot1][numSiteResidues][ligAANum][curLigPos] + pairwiseMaxEnergyMatrix[secondPos][secondAltAA1][secondAltRot1][numSiteResidues][ligAANum][curLigPos];
					
						if ((curEmin-curEmax) < minE)
							minE = curEmin-curEmax;
						
						found = true;
					}
				}
			}
			
			//if(minE==bigE)//make it contribute nothing to the sum
			//	minE = 0;
		}
		
		if (!found)
			minE = 0.0; //contributes nothing to the sum
		
		return minE;
	}*/
	//////////////////////////////////////////////////////////////////////////////////
	
	//////////////////////////////////////////////////////////////////////////////////		
	//Called only by CanEliminate() adn CanEliminateLig()
	private double SumMaxIndInt (int withoutPos1, int withoutPos2){
		
		double sum = 0;
		
		//get the contribution from the active site residue rotamers
		for (int curPos=0; curPos<numMutable; curPos++){
			
			if ((curPos != withoutPos1)&&(curPos != withoutPos2)){
				
				sum += maxScale*MaxIndInt(curPos);
			}
		}
		
		/*if (numLigRot!=0){ //there is a ligand
			//get the contribution from the ligand rotamers: there is only one ligand residue,
			//so there is only one position j here for which to add
			if ((withoutPos1!=numSiteResidues)&&(withoutPos2!=numSiteResidues)) //if we are not currently checking ligand rotamers for pruning
				sum += maxScale*LigandMaxIndInt();
		}*/
	
		return sum;
	}
	
	//Called by SumMaxIndInt(.)
	private double MaxIndInt (int atPos){
		
		double maxEInt = -999999.0;//this is OK, since E intervals are always positive
		double curEInt;
		
		int numRotForAAatPos;
		
		int str=mutRes2Strand[atPos];
		int strResNum=strandMut[str][mutRes2MutIndex[atPos]];
		
		for (int AA=0; AA<numAAtypes[atPos]; AA++){
			
			int curAA = strandRot[str].getIndexOfNthAllowable(strResNum,AA);
			
			numRotForAAatPos = getNumRot( str, strResNum, curAA );
			
			for (int curRot=0; curRot<numRotForAAatPos; curRot++){
			
				curEInt = IndInt(atPos, curAA, curRot);
				if (curEInt > maxEInt){
					maxEInt = curEInt;
				}
			}
		}
		
		return maxEInt;
	}

	//Called by MaxIndInt(.)
	private double IndInt (int atPos, int atAA, int atRot){
		
		//s at j
		Index3 index1 = new Index3(atPos,atAA,atRot);//atPos*numTotalRot + rotIndOffset[atAA] + atRot;
		
		if ((!eliminatedRotAtPos.get(index1))){ //not pruned 
		
			double maxE = pairwiseMaxEnergyMatrix.getIntraAndShellE( atPos, atAA, atRot);
                        double minE = pairwiseMinEnergyMatrix.getIntraAndShellE( atPos, atAA, atRot);
			
			return (maxE - minE);
		}
		else
			return 0.0;//contributes nothing
	}
	

	//////////////////////////////////////////////////////////////////////////////////
	
	//////////////////////////////////////////////////////////////////////////////////
	//Called only by CanEliminate() and CanEliminateLig()
	private double SumSumMaxPairInt(int withoutPos1, int withoutPos2){
		
		double sum = 0;
		
		//get the contribution from the active site residue rotamers
		for (int curPos1=0; curPos1<numMutable; curPos1++){
			if ((curPos1 != withoutPos1)&&(curPos1 != withoutPos2)){
				for (int curPos2=0; curPos2<curPos1; curPos2++){
					if ((curPos2 != withoutPos1)&&(curPos2 != withoutPos2)){
					
						sum += maxScale*MaxPairInt(curPos1,curPos2);
					}
				}
			}
		}
		
		/*if (numLigRot!=0){ //there is a ligand
			//get the contribution from the ligand rotamers: there is only one ligand residue,
			//so there is only one position k here for which to add;
			//the range of j is the number of active site residues
			if ((withoutPos1!=numSiteResidues)&&(withoutPos2!=numSiteResidues)){ //if we are not currently checking ligand rotamers for pruning
				for (int curPos=0; curPos<numSiteResidues; curPos++){
					if ((curPos != withoutPos1)&&(curPos != withoutPos2)){
						
						sum += maxScale*LigandMaxPairInt(curPos);
					}
				}
			}
		}*/
		
		return sum;
	}
	
	//Called by SumSumMaxPairInt(.)
	private double MaxPairInt (int atPos1, int atPos2){
		
		double maxEInt = -999999.0;//this is OK, since E intervals are always positive
		double curEInt;
		
		int numRotForAAatPos1;
		
		int str1=mutRes2Strand[atPos1];
		int strResNum1=strandMut[str1][mutRes2MutIndex[atPos1]];
		int str2=mutRes2Strand[atPos2];
		int strResNum2=strandMut[str2][mutRes2MutIndex[atPos2]];
		
		for (int AA1=0; AA1<numAAtypes[atPos1]; AA1++){
			
			int curAA1 = strandRot[str1].getIndexOfNthAllowable(strResNum1,AA1);
			
			numRotForAAatPos1 = getNumRot( str1, strResNum1, curAA1 );
                        
			for (int curRot1=0; curRot1<numRotForAAatPos1; curRot1++){
				
				int numRotForAAatPos2;
				
				for (int AA2=0; AA2<numAAtypes[atPos2]; AA2++){
					
					int curAA2 = strandRot[str2].getIndexOfNthAllowable(strResNum2,AA2);
					
					numRotForAAatPos2 = getNumRot( str2, strResNum2, curAA2 );
					
					for (int curRot2=0; curRot2<numRotForAAatPos2; curRot2++){
					
						curEInt = PairInt(atPos1, curAA1, curRot1, atPos2, curAA2, curRot2);
						if (curEInt > maxEInt){
							maxEInt = curEInt;
						}
					}
				}
			}
		}
		
		return maxEInt;
	}
	
	//Called by MaxPairInt(.)
	private double PairInt (int atPos1, int atAA1, int atRot1, int atPos2, int atAA2, int atRot2){
		
		//There is a displacement: colum 0 and row 0 have special entries, 
		//so pairwise energies start from row 1, column 1
		Index3 index1 = new Index3(atPos1,atAA1,atRot1);//atPos1*numTotalRot + rotIndOffset[atAA1] + atRot1;//u at k
		Index3 index2 = new Index3(atPos2,atAA2,atRot2);//atPos2*numTotalRot + rotIndOffset[atAA2] + atRot2;//s at j
		
		if (((!eliminatedRotAtPos.get(index1)))&&
				((!eliminatedRotAtPos.get(index2)))){ //not pruned 
		
			double maxE = pairwiseMaxEnergyMatrix.getPairwiseE( atPos1, atAA1, atRot1, atPos2, atAA2, atRot2 );
			double minE = pairwiseMinEnergyMatrix.getPairwiseE( atPos1, atAA1, atRot1, atPos2, atAA2, atRot2 );
		
			//if ((maxE<=stericEThreshPair)&&(minE<=stericEThreshPair))
				return (maxE - minE);
		}
		else
			return 0.0;//contributes nothing
	}
	
}
