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
//	DEESplit1f.java
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

import java.util.Random;

/**
 * Performs full split-DEE (conformational splitting) with 1 plit position
 *
*/
public class DEESplit1f extends DEE {



    //constructor
    DEESplit1f(PairwiseEnergyMatrix arpMatrix, PairwiseEnergyMatrix arpMatrixMax, int numResMutable, int strMut[][], float initEw,
               StrandRotamers strandLRot[], PrunedRotamers<Boolean> prunedRotAtRes, boolean residueMut[],
               boolean doMin, double indInt[], double pairInt[], boolean spFlags[][][][][][], boolean useSF, boolean dDEE, boolean minBB,
               int mutRes2StrandP[], int mutRes2MutIndexP[], boolean typeDep, boolean aIMinDEE, float aIval, boolean doPerts) {


        init(arpMatrix, arpMatrixMax, numResMutable,
             strMut, initEw, strandLRot, prunedRotAtRes, doMin, indInt, pairInt,
             spFlags, useSF, minBB, mutRes2StrandP, mutRes2MutIndexP, typeDep, aIMinDEE, aIval,
             false, dDEE, residueMut, null, doPerts);


        compNumRotForRes(); //compute the number of rotamers for each residue position
    }



    //return the split flags for all rotamer pairs
    public boolean[][][][][][] getSplitFlags() {
        return splitFlags;
    }

    //Compute the conformations that can be eliminated
    //Return a boolean matrix in which an element is true if
    //the corresponding r at i can be eliminated, and false otherwise
    public PrunedRotamers<Boolean> ComputeEliminatedRotConf() {

        int numRotForCurAAatPos;

        int prunedCurRun = 0;
        boolean done = false;
        numRuns = 1;

        while (!done) {

            prunedCurRun = 0;

            System.out.println("Current run: "+numRuns);

            //Compute for the AS residues first
            for (int curPos=0; curPos<numMutable; curPos++) {

                int str=mutRes2Strand[curPos];
                int strResNum=strandMut[str][mutRes2MutIndex[curPos]];

                if ((!distrDEE)||(resInPair[curPos])) {//not distributed DEE or the residue to be checked for distrDEE

                    //System.out.print("Starting AS residue "+curPos);

                    for (int AA=0; AA<numAAtypes[curPos]; AA++) {

                        //System.out.print(".");

                        int curAA = strandRot[str].getIndexOfNthAllowable(strResNum,AA);

                        //find how many rotamers are allowed for the current AA type at the given residue;
                        //note that ala and gly have 0 possible rotamers
                        numRotForCurAAatPos = getNumRot( str, strResNum, curAA );

                        for(int curRot=0; curRot<numRotForCurAAatPos; curRot++) {

                            if (!eliminatedRotAtPos.get(curPos,curAA,curRot)) { //not already pruned

                                if (CanEliminate(curPos, curAA, curRot)) {
                                    eliminatedRotAtPos.set(curPos,curAA,curRot,true);
                                    //System.out.println(curEc);

                                    prunedCurRun++;
                                } else
                                    eliminatedRotAtPos.set(curPos,curAA,curRot,false);
                            }
                        }
                    }
                    //System.out.println("done");
                }
            }//System.out.println("Ec: "+Ec);

            //If there is a ligand, compute MinDEE for the lig rotamers as well
            /*if (numLigRot!=0){

            	if ((!distrDEE)||(resInPair[numSiteResidues])) {//not distributed DEE or ligand to be checked for distrDEE

            		System.out.print("Starting ligand run");
            		System.out.print("..");
            		for (int curRot=0; curRot<numLigRot; curRot++){
            			if (!eliminatedRotAtPos[numSiteResidues*numTotalRot+curRot]){//not already pruned

            				if (CanEliminateLig(curRot)){
            					eliminatedRotAtPos[numSiteResidues*numTotalRot+curRot] = true;

            					prunedCurRun++;
            				}
            				else
            					eliminatedRotAtPos[numSiteResidues*numTotalRot+curRot] = false;
            			}
            		}
            		System.out.println("done");
            	}
            }*/

            System.out.println("Number of rotamers pruned this run: "+prunedCurRun);
            System.out.println("DEE: The minimum difference is "+minDiff);
            System.out.println();

            if (prunedCurRun==0) //no rotamers pruned this run, so done
                done = true;
            else
                numRuns++;
        }

        return eliminatedRotAtPos;
    }

    //Called only by ComputeEliminatedRotConf(.)
    /*
     * The logic is as follows:
     * 	for every AS residue
     * 		for every AA type at this residue
     * 			for every possible rotamer at the given AA type
     * 				check against every other possible rotamer for all AA types at the same residue;
     * 				eliminate if provable
     *
     * That is, for each residue: out of the 152 rotamer possibilities, choose 1
     * and compare it to the other 151 until elimination can be proven or there
     * are no more comparisons left. Repeat this for all 152 rotamers and all AS residues
     *
     */
    private boolean CanEliminate (int posNum, int AANumAtPos, int rotNumAtPos) {

        double minIndVoxelE, maxIndVoxelE;
        double indVoxelInterval, pairVoxelInterval;
        double minDiffPairVoxelE;


        double checkSum;

        int str=mutRes2Strand[posNum];
        int strResNum=strandMut[str][mutRes2MutIndex[posNum]];

        //In the energy matrix, column 0 gives the individual energies for each r at i;
        //skip row 0, as the individual energies start from row 1 (and are all in column 0)
        //index_r = posNum*numTotalRot + rotIndOffset[AANumAtPos] + rotNumAtPos;

        if ((!eliminatedRotAtPos.get(posNum,AANumAtPos,rotNumAtPos))) { //not already pruned

            minIndVoxelE = pairwiseMinEnergyMatrix.getIntraAndShellE( posNum, AANumAtPos, rotNumAtPos ); 	//formula term 1.  Intra+shell

            if ((minIndVoxelE)>=stericE) //rotamer incompatible with template, so prune
                return true;

            // 2010: Don't compute the intervals if useMinDEEPruningEw is true since they are not necessary
            if (doMinimize && !doIMinDEE) { //MinDEE, so compute the interval terms
                indVoxelInterval = indIntMinDEE[posNum];							//formula term 3
                pairVoxelInterval = pairIntMinDEE[posNum];							//formula term 4
            } else { //traditional-DEE, so no interval terms // 2010: if useMinDEEPruningEw is set to true then it is the same as DEE
                indVoxelInterval = 0.0;
                pairVoxelInterval = 0.0;
            }


            //int splitPos = chooseSplitPos(posNum,majorSplitPos,index_r); //choose the split position
            for (int splitPos=0; splitPos<numMutable; splitPos++) { //for each splitting position
                boolean found = true;

                int str2=mutRes2Strand[splitPos];
                int strResNum2=strandMut[str2][mutRes2MutIndex[splitPos]];


                /*for (int i=0; i<majorSplitPos.length; i++){ //split residue different from majorSplitPos
                	if (splitPos==majorSplitPos[i])
                		found = false;
                }*/
                if ((splitPos!=posNum)&&(found)) { //the split residue must be different from cur res

                    boolean partitionPruned[] = new boolean[numRotForRes[splitPos]]; //prune r at i for each partition
                    for (int partP=0; partP<partitionPruned.length; partP++)
                        partitionPruned[partP] = false;

                    //System.out.println("Splitting residue: "+splitPos);

                    //For the particular position, compare the energy performance (one by one)
                    //of the remaining rotamer possibilities to that of the given rotamer:
                    //given r at i, compare it to all t at i for pruning
                    int numRotForAAatPos;

                    for (int AA=0; AA<numAAtypes[posNum]; AA++) {

                        int altAA = strandRot[str].getIndexOfNthAllowable(strResNum,AA);

                        if(!typeDependent || (altAA == AANumAtPos )) {

                            numRotForAAatPos = getNumRot( str, strResNum, altAA );

                            for (int altRot=0; altRot<numRotForAAatPos; altRot++) {

                                //if t and r are not actually the same rotamer of the same AA
                                if (!((altAA==AANumAtPos)&&(altRot==rotNumAtPos))) {

                                    //at this point, we know what r at i and t at i are

                                    //index_t = posNum*numTotalRot + rotIndOffset[altAA] + altRot;

                                    maxIndVoxelE = pairwiseMaxEnergyMatrix.getIntraAndShellE( posNum, altAA, altRot );		//formula term 2

                                    //if ((maxIndVoxelE<=stericEThreshIntra)&&(maxShellResE<=stericEThreshPair)){//check only if not an unallowed steric
                                    if ((!eliminatedRotAtPos.get(posNum,altAA,altRot))) { //not pruned

                                        minDiffPairVoxelE = SumMinDiffPVE(posNum, AANumAtPos, rotNumAtPos, altAA, altRot, splitPos);	//formula term 5


                                        int partCount = 0; //the visited partitions count
                                        for (int spAA=0; spAA<numAAtypes[splitPos]; spAA++) { //for each AA at the splitting residue

                                            int partAA = strandRot[str2].getIndexOfNthAllowable(strResNum2,spAA);
                                            int numRotForPartAA = getNumRot( str2, strResNum2, partAA );

                                            for (int partRot=0; partRot<numRotForPartAA; partRot++) { //for each rot for the given partitioning AA

                                                if (!partitionPruned[partCount]) { //only if r at i not pruned for this partition yet

                                                    //int index_hv = splitPos*numTotalRot + rotIndOffset[partAA] + partRot; //the index of the partitioning rotamer

                                                    if ((!eliminatedRotAtPos.get(splitPos,partAA,partRot))&&((!useFlags)||(!splitFlags[posNum][AANumAtPos][rotNumAtPos][splitPos][partAA][partRot]))) { //not pruned

                                                        double splitPosDiffE = pairwiseMinEnergyMatrix.getPairwiseE( posNum, AANumAtPos, rotNumAtPos, splitPos, partAA, partRot )
                                                                               - pairwiseMaxEnergyMatrix.getPairwiseE( posNum, altAA, altRot, splitPos, partAA, partRot ); //formula term 6

                                                        checkSum = -templateInt + minIndVoxelE - maxIndVoxelE
                                                                   - indVoxelInterval - pairVoxelInterval + minDiffPairVoxelE
                                                                   + splitPosDiffE;

                                                        if (checkSum > curEw) {
                                                            //System.out.println(index_r+" "+index_t+" "+checkSum+" "+minIndVoxelE+" "+minShellResE+" "+maxIndVoxelE+" "+maxShellResE+" "+indVoxelInterval+" "+pairVoxelInterval+" "+minDiffPairVoxelE);
                                                            partitionPruned[partCount] = true; //r at i can be pruned for this partition

                                                            if (useFlags) {
                                                                splitFlags[posNum][AANumAtPos][rotNumAtPos][splitPos][partAA][partRot] = true; //flag as a dead-ending pair
                                                                splitFlags[splitPos][partAA][partRot][posNum][AANumAtPos][rotNumAtPos] = true;
                                                            }
                                                        }

                                                        else {
                                                            minDiff = Math.max(minDiff,checkSum);
                                                        }
                                                    } else //v at h is pruned, so we set the pruning of r at i for this partition to true
                                                        partitionPruned[partCount] = true;
                                                }
                                                partCount++; //the next partition
                                            }
                                        }

                                        //after checking all partitions, return if r at i is pruned for all of the partitions;
                                        //	otherwise, go to the next competitor t at i
                                        boolean canPrune = true;
                                        for (int curPartCheck=0; curPartCheck<partitionPruned.length; curPartCheck++) {
                                            if (!partitionPruned[curPartCheck]) {
                                                canPrune = false;
                                                break;
                                            }
                                        }
                                        if (canPrune)
                                            return true;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        } else //aready pruned
            return true;

        //We have tried all of the other rotamers at the current position and none
        //of them is able to prune the given rotamer, so we return false
        return false;
    }



    ////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////////

    //Called only by CanEliminate(.)
    private double SumMinDiffPVE (int atPos, int withAA1, int withRot1, int withAA2, int withRot2, int splitPos) {

        double sum = 0;

        //get the contribution from the active site residue rotamers
        for (int curPos=0; curPos<numMutable; curPos++) {

            if ((curPos != atPos)&&(curPos!=splitPos)) // j!=i and j!=k

                sum += IndMinDiffPVE(atPos, withAA1, withRot1, withAA2, withRot2, curPos);
        }

        /*if (numLigRot!=0){ //there is a ligand
        	//get the contribution from the ligand rotamers: there is only one ligand residue,
        	//so there is only one position j here for which to add
        	sum += LigandIndMinDiffPVE (atPos, withAA1, withRot1, withAA2, withRot2);
        }*/

        return sum;
    }

    //Called by SumMaxMaxPVE(.)
    private double IndMinDiffPVE (int firstPos, int firstAA1, int firstRot1, int firstAA2, int firstRot2, int secondPos) {

        double minE = bigE;
        double curEmin, curEmax;

        int numRotForAAatPos;

        //r at i
        //index1 = firstPos*numTotalRot + rotIndOffset[firstAA1] + firstRot1;

        int str2=mutRes2Strand[secondPos];
        int strResNum2=strandMut[str2][mutRes2MutIndex[secondPos]];


        //t at i
        //index3 = firstPos*numTotalRot + rotIndOffset[firstAA2] + firstRot2;

        boolean found = false;

        if (((!eliminatedRotAtPos.get(firstPos,firstAA1,firstRot1)))&&
                ((!eliminatedRotAtPos.get(firstPos,firstAA2,firstRot2)))) { //not pruned

            for (int AA=0; AA<numAAtypes[secondPos]; AA++) {

                int curAA = strandRot[str2].getIndexOfNthAllowable(strResNum2,AA);

                numRotForAAatPos = getNumRot( str2, strResNum2, curAA );

                for (int curRot=0; curRot<numRotForAAatPos; curRot++) {

                    //s at j
                    //index2 = secondPos*numTotalRot + rotIndOffset[curAA] + curRot;

                    if ((!eliminatedRotAtPos.get(secondPos,curAA,curRot))) { //not pruned

                        if ((!useFlags)||(!splitFlags[firstPos][firstAA1][firstRot1][secondPos][curAA][curRot])) { //not using split flags or not flagged

                            curEmin = pairwiseMinEnergyMatrix.getPairwiseE( firstPos, firstAA1, firstRot1, secondPos, curAA, curRot );
                            curEmax = pairwiseMaxEnergyMatrix.getPairwiseE( firstPos, firstAA2, firstRot2, secondPos, curAA, curRot );
                            //if (/*(curEmin<=stericEThreshPair)&&*/(curEmax<=stericEThreshPair)){//check only if not an unallowed steric
                            if ((curEmin-curEmax) < minE)
                                minE = curEmin-curEmax;
                            //}
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
    /*private double LigandIndMinDiffPVE (int firstPos, int firstAA1, int firstRot1, int firstAA2, int firstRot2){

    	double minE = bigE;
    	double curEmin, curEmax;

    	int index1, index2, index3;

    	//r at i
    	index1 = firstPos*numTotalRot + rotIndOffset[firstAA1] + firstRot1;

    	//t at i
    	index3 = firstPos*numTotalRot + rotIndOffset[firstAA2] + firstRot2;

    	boolean found = false;

    	if (((!eliminatedRotAtPos[index1]))&&
    			((!eliminatedRotAtPos[index3]))){ //not pruned

    		for (int curLigPos=0; curLigPos<numLigRot; curLigPos++){

    			//s at j (the ligand residue)
    			index2 = numSiteResidues*numTotalRot + curLigPos;

    			if ((!eliminatedRotAtPos[index2])){ //not pruned

    				if ((!useFlags)||(!splitFlags[index1][index2])){ //not using split flags or not flagged

    					curEmin = pairwiseMinEnergyMatrix[firstPos][firstAA1][firstRot1][numSiteResidues][ligAANum][curLigPos];
    					curEmax = pairwiseMaxEnergyMatrix[firstPos][firstAA2][firstRot2][numSiteResidues][ligAANum][curLigPos];
    					//if (/*(curEmin<=stericEThreshPair)&&*//*(curEmax<=stericEThreshPair)){//check only if not an unallowed steric
						if ((curEmin-curEmax) < minE)
							minE = curEmin-curEmax;
					//}
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

    ///////////////////////////////////////////////////////////////////////////////////////
    //Same as CanEliminate(), just checks the ligand rotamers for pruning
    //Called by ComputeEliminatedRotConf()
    /*private boolean CanEliminateLig (int curLigRot){

    	double minIndVoxelE, maxIndVoxelE;
    	double minShellResE, maxShellResE;
    	double indVoxelInterval, pairVoxelInterval;
    	double minDiffPairVoxelE;

    	int index_r, index_t;

    	double checkSum;

    	//In the energy matrix, column 0 gives the individual energies for each r at i;
    	//skip row 0, as the individual energies start from row 1 (and are all in column 0)
    	index_r = numSiteResidues*numTotalRot + curLigRot;

    	if ((!eliminatedRotAtPos[index_r])){ //not already pruned

    		minIndVoxelE = pairwiseMinEnergyMatrix[numSiteResidues][ligAANum][curLigRot][numSiteResidues][0][0]; //formula term 1
    		minShellResE = pairwiseMinEnergyMatrix[numSiteResidues][ligAANum][curLigRot][numSiteResidues][0][1];

    		if ((minIndVoxelE + minShellResE)>=stericE) //rotamer incompatible with template, so prune
    			return true;

    		if (doMinimize){ //MinDEE, so compute the interval terms
    			indVoxelInterval = indIntMinDEE[numSiteResidues];							//formula term 3
    			pairVoxelInterval = pairIntMinDEE[numSiteResidues];							//formula term 4
    		}
    		else { //traditional-DEE, so no interval terms
    			indVoxelInterval = 0.0;
    			pairVoxelInterval = 0.0;
    		}

    		//For the particular position, compare the energy performance (one by one)
    		//of the remaining rotamer possibilities to that of the given rotamer:
    		//given r at i, compare it to all t at i for pruning
    		for (int altRot=0; altRot<numLigRot; altRot++){

    			//if t and r are not actually the same lig rotamer
    			if (curLigRot!=altRot){

    				//at this point, we know what r at i and t at i are

    				index_t = numSiteResidues*numTotalRot + altRot;

    				maxIndVoxelE = pairwiseMaxEnergyMatrix[numSiteResidues][ligAANum][altRot][numSiteResidues][0][0]; //formula term 2
    				maxShellResE = pairwiseMaxEnergyMatrix[numSiteResidues][ligAANum][altRot][numSiteResidues][0][1];

    				//if ((maxIndVoxelE<=stericEThreshIntra)&&(maxShellResE<=stericEThreshPair)){//check only if not an unallowed steric
    				if ((!eliminatedRotAtPos[index_t])){ //not pruned

    					minDiffPairVoxelE = SumMinDiffPVELig(curLigRot, altRot);	//formula term 5

    					checkSum = -templateInt + (minIndVoxelE + minShellResE) - (maxIndVoxelE + maxShellResE)
    								- indVoxelInterval - pairVoxelInterval + minDiffPairVoxelE;

    					if (checkSum > curEw){
    						//System.out.println(checkSum+" "+minIndVoxelE+" "+minShellResE+" "+maxIndVoxelE+" "+maxShellResE+" "+indVoxelInterval+" "+pairVoxelInterval+" "+minDiffPairVoxelE);

    						return true;}//this rotamer can be pruned/eliminated
    					else {
    						minDiff = Math.max(minDiff,checkSum);
    					}
    				}
    			}
    		}
    	}
    	else //aready pruned
    		return true;

    	//We have tried all of the other rotamers at the current position and none
    	//of them is able to prune the given rotamer, so we return false
    	return false;
    }*/

    //Same as SumMinDiffPVE(), just checks the ligand rotamers for pruning;
    //Called by CanEliminateLig()
    /*private double SumMinDiffPVELig (int withRot1, int withRot2){

    	double sum = 0;

    	//get the contribution from the active site residue rotamers
    	for (int curPos=0; curPos<numSiteResidues; curPos++){

    			sum += IndMinDiffPVELig(withRot1, withRot2, curPos);
    	}

    	return sum;
    }*/

    //Same as IndMinDiffPVE(), just checks the ligand rotamers for pruning
    //Called by SumMinDiffPVELig()
    /*private double IndMinDiffPVELig (int firstRot1, int firstRot2, int secondPos){

    	double minE = bigE;
    	double curEmin, curEmax;

    	int index1, index2, index3;
    	int numRotForAAatPos;

    	//r at i
    	index1 = numSiteResidues*numTotalRot + firstRot1;

    	//t at i
    	index3 = numSiteResidues*numTotalRot + firstRot2;

    	boolean found = false;

    	if (((!eliminatedRotAtPos[index1]))&&
    			((!eliminatedRotAtPos[index3]))){ //not pruned

    		for (int AA=0; AA<numAAtypes[secondPos]; AA++){

    			int curAA = sysLR.getIndexOfNthAllowable(residueMap[secondPos],AA);

    			numRotForAAatPos = rl.getNumRotForAAtype(curAA);
    			if (numRotForAAatPos==0)	//ala or gly
    				numRotForAAatPos = 1;

    			for (int curRot=0; curRot<numRotForAAatPos; curRot++){

    				//There is a displacement: column 0 and row 0 have special entries,
    				//so pairwise energies start from row 1, column 1

    				//s at j
    				index2 = secondPos*numTotalRot + rotIndOffset[curAA] + curRot;

    				if ((!eliminatedRotAtPos[index2])){ //not pruned

    					if ((!useFlags)||(!splitFlags[index1][index2])){ //not using split flags or not flagged

    						curEmin = pairwiseMinEnergyMatrix[numSiteResidues][ligAANum][firstRot1][secondPos][curAA][curRot];
    						curEmax = pairwiseMaxEnergyMatrix[numSiteResidues][ligAANum][firstRot2][secondPos][curAA][curRot];
    						//if (/*(curEmin<=stericEThreshPair)&&*//*(curEmax<=stericEThreshPair)){//check only if not an unallowed steric
							if ((curEmin-curEmax) < minE)
								minE = curEmin-curEmax;
						//}
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
}*/
}
