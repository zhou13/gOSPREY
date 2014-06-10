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
//	DEEGoldstein.java
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
 * Performs simple Goldstein DEE rotamer pruning
 *
*/
public class DEEGoldstein extends DEE {

    //constructor
    DEEGoldstein(PairwiseEnergyMatrix arpMatrix, PairwiseEnergyMatrix arpMatrixMax, int numResMutable,
                 int strMut[][], float initEw,
                 StrandRotamers strandLRot[], PrunedRotamers<Boolean> prunedRotAtRes, boolean doMin, double indInt[],
                 double pairInt[], boolean spFlags[][][][][][], boolean useSF, boolean minBB,
                 int mutRes2StrandP[], int mutRes2MutIndexP[], boolean typeDep, boolean iMinDEE, float Ival, boolean doPerts) {


        init(arpMatrix, arpMatrixMax, numResMutable,
             strMut, initEw, strandLRot, prunedRotAtRes, doMin, indInt,
             pairInt, spFlags, useSF, minBB, mutRes2StrandP, mutRes2MutIndexP, typeDep, iMinDEE, Ival,
             false, false, null, null, doPerts);
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

            //System.out.println("Current run: "+numRuns);

            //Compute for the AS residues first
            for (int curPos=0; curPos<numMutable; curPos++) {

                int str=mutRes2Strand[curPos];
                int strResNum=strandMut[str][mutRes2MutIndex[curPos]];
                //System.out.print("Starting AS residue "+curPos);

                for (int AA=0; AA<numAAtypes[curPos]; AA++) {

                    //System.out.print(".");

                    int curAA = strandRot[str].getIndexOfNthAllowable(strResNum,AA);

                    //find how many rotamers are allowed for the current AA type at the given residue;
                    //note that ala and gly have 0 possible rotamers
                    numRotForCurAAatPos = getNumRot(str, strResNum, curAA);

                    for(int curRot=0; curRot<numRotForCurAAatPos; curRot++) {

                        if (!eliminatedRotAtPos.get(curPos,curAA,curRot)) { //not already pruned

                            if (CanEliminate(curPos, curAA, curRot)) {
                                eliminatedRotAtPos.set(curPos,curAA,curRot,true);
                                prunedCurRun++;
                            } else
                                eliminatedRotAtPos.set(curPos,curAA,curRot,false);
                        }
                    }
                }
                //System.out.println("done");
            }//System.out.println("Ec: "+Ec);

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

        //In the energy matrix, column 0 gives the individual energies for each r at i;
        //skip row 0, as the individual energies start from row 1 (and are all in column 0)
        //index_r = posNum*numTotalRot + rotIndOffset[AANumAtPos] + rotNumAtPos;

        if ((!eliminatedRotAtPos.get(posNum,AANumAtPos,rotNumAtPos))) { //not pruned

            minIndVoxelE = pairwiseMinEnergyMatrix.getIntraAndShellE( posNum, AANumAtPos, rotNumAtPos );//formula term 1
            //intra+shell energy


            if ( minIndVoxelE >=stericE) //rotamer incompatible with template, so prune
                return true;

            // 2010: Don't compute the intervals if doIMinDEE is true since they are not necessary
            if (doMinimize && !doIMinDEE) { //MinDEE, so compute the interval terms
                indVoxelInterval = indIntMinDEE[posNum];							//formula term 3
                pairVoxelInterval = pairIntMinDEE[posNum];							//formula term 4
            } else { //traditional-DEE, so no interval terms
                indVoxelInterval = 0.0;
                pairVoxelInterval = 0.0;
            }

            //For the particular position, compare the energy performance (one by one)
            //of the remaining rotamer possibilities to that of the given rotamer:
            //given r at i, compare it to all t at i for pruning
            int numRotForAAatPos;

            int str=mutRes2Strand[posNum];
            int strResNum=strandMut[str][mutRes2MutIndex[posNum]];

            for (int AA=0; AA<numAAtypes[posNum]; AA++) {

                int altAA = strandRot[str].getIndexOfNthAllowable(strResNum,AA);
                if(!typeDependent || altAA == AANumAtPos) {

                    numRotForAAatPos = getNumRot(str,strResNum,altAA);

                    for (int altRot=0; altRot<numRotForAAatPos; altRot++) {

                        //if t and r are not actually the same rotamer of the same AA
                        if (!((altAA==AANumAtPos)&&(altRot==rotNumAtPos))) {

                            //at this point, we know what r at i and t at i are

                            maxIndVoxelE = pairwiseMaxEnergyMatrix.getIntraAndShellE( posNum, altAA, altRot );//formula term 2

                            //if ((maxIndVoxelE<=stericEThreshIntra)&&(maxShellResE<=stericEThreshPair)){//check only if not an unallowed steric
                            if ((!eliminatedRotAtPos.get(posNum,altAA,altRot))) { //not pruned

                                minDiffPairVoxelE = SumMinDiffPVE(posNum, AANumAtPos, rotNumAtPos, altAA, altRot);	//formula term 5

                                checkSum = -templateInt + minIndVoxelE - maxIndVoxelE
                                           - indVoxelInterval - pairVoxelInterval + minDiffPairVoxelE;

                                if (checkSum > curEw) {
                                    return true;
                                }//this rotamer can be pruned/eliminated
                                else {
                                    minDiff = Math.max(minDiff,checkSum);
                                }
                            }
                        }
                    }
                }
            }
        } else //already pruned
            return true;

        //We have tried all of the other rotamers at the current position and none
        //of them is able to prune the given rotamer, so we return false
        return false;
    }

    ////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////////

    //Called only by CanEliminate(.)
    private double SumMinDiffPVE (int atPos, int withAA1, int withRot1, int withAA2, int withRot2) {

        double sum = 0;

        //get the contribution from the active site residue rotamers
        for (int curPos=0; curPos<numMutable; curPos++) {

            if (curPos != atPos)

                sum += IndMinDiffPVE(atPos, withAA1, withRot1, withAA2, withRot2, curPos);
        }

        return sum;
    }

    //Called by SumMaxMaxPVE(.)
    //min_s (E_\ominus(i_r,j_s)-E_\oplus(i_t,j_s))
    //secondPos = j
    private double IndMinDiffPVE (int firstPos, int firstAA1, int firstRot1, int firstAA2, int firstRot2, int secondPos) {

        double minE = bigE;
        double curEmin, curEmax;

        int numRotForAAatPos;

        int str2=mutRes2Strand[secondPos];
        int strResNum2=strandMut[str2][mutRes2MutIndex[secondPos]];


        boolean found = false;

        if (((!eliminatedRotAtPos.get(firstPos,firstAA1,firstRot1)))&&
                ((!eliminatedRotAtPos.get(firstPos,firstAA2,firstRot2)))) { //not pruned

            for (int AA=0; AA<numAAtypes[secondPos]; AA++) {

                int curAA = strandRot[str2].getIndexOfNthAllowable(strResNum2,AA);

                numRotForAAatPos = getNumRot(str2, strResNum2, curAA);

                for (int curRot=0; curRot<numRotForAAatPos; curRot++) {

                    //There is a displacement: column 0 and row 0 have special entries,
                    //so pairwise energies start from row 1, column 1

                    //s at j
                    //index2 = strandRot[str2].getIndex(strPos2, curAA, curRot, strandMut);//secondPos*numTotalRot + rotIndOffset[curAA] + curRot;

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

}