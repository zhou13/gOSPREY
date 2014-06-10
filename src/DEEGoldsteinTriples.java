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
//	DEEGoldsteinTriples.java
//
//	Version:           2.1 beta
//
//
//	  authors:
// 	  initials    name                 organization                email
//	 ---------   -----------------    ------------------------    ----------------------------
//	  MAH           Mark A. Hallen	  Duke University               mah43@duke.edu
///////////////////////////////////////////////////////////////////////////////////////////////

/**
 * Pruned triples of rotamers (for "renormalized residues" as referred to by Goldstein 1994)
 * Uses rigid DEE or iMinDEE
*/



public class DEEGoldsteinTriples extends DEE {

    int prunedCurRun;
    //Number of triples pruned this run

    //For distributed DEE
    boolean resInTriple[];

    int magicBulletNum = 5;//Number of magic bullets per triple of positions

    int numPrunedRegular = 0;
    int numPrunedMin = 0;

    //constructor
    DEEGoldsteinTriples (PairwiseEnergyMatrix arpMatrix, int numResMutable, int strMut[][], float initEw,
                         StrandRotamers strandLRot[], PrunedRotamers<Boolean> prunedRotAtRes, boolean residueMut[],
                         boolean doMin, boolean spFlags[][][][][][], boolean mb, int mbNum, boolean dDEE,
                         boolean minBB, int mutRes2StrandP[], int mutRes2MutIndexP[], boolean typeDep, boolean aIMinDEE,
                         float aIval, boolean tripFlags[][][][][][][][][], boolean doPerts) {


        init(arpMatrix, null, numResMutable,
             strMut, initEw, strandLRot, prunedRotAtRes, doMin, null, null,
             spFlags, true, minBB, mutRes2StrandP, mutRes2MutIndexP, typeDep, aIMinDEE, aIval,
             mb, dDEE, residueMut, tripFlags, doPerts);

        if(doMinimize && (!doIMinDEE) ) {
            System.err.println("ERROR: iMinDEE required for pruning triples");
            System.exit(1);
        }

        useTriples = true;//We are using triples even if tripleFlags was initialized to null (initializeTripleFlags can be called later)
        //This still needs to be set to true because it is used by DEE.isPrunedTriple

        resInTriple= residueMut;

        magicBulletNum = mbNum;
    }


    //return the triple flags for all rotamer pairs
    public boolean[][][][][][][][][] getTripleFlags() {
        return tripleFlags;
    }

    //Compute the conformations that can be eliminated
    //Saves results to tripleFlags
    public void ComputeEliminatedRotConf() {

        boolean done = false;
        numRuns = 1;

        int magicBulletAA[][] = new int[magicBulletNum][3];
        int magicBulletRot[][] = new int[magicBulletNum][3];

        while (!done) {

            prunedCurRun = 0;

            System.out.println("Current run: "+numRuns);


            int numRotForCurAAatPos1;

            for (int curPos1=0; curPos1<numMutable; curPos1++) {
                System.out.println("");
                System.out.print("Pos1 "+curPos1+": ");

                if ((magicBullet)||(!distrDEE)||(resInTriple[curPos1])) { //mb-pairs or not distrDEE or cur res is in distr pair

                    int str1=mutRes2Strand[curPos1];
                    int strResNum1=strandMut[str1][mutRes2MutIndex[curPos1]];

                    for (int curPos2=curPos1+1; curPos2<numMutable; curPos2++) {

                        if ((magicBullet)||(!distrDEE)||(resInTriple[curPos2])) {

                            int str2=mutRes2Strand[curPos2];
                            int strResNum2=strandMut[str2][mutRes2MutIndex[curPos2]];

                            for(int curPos3=curPos2+1; curPos3<numMutable; curPos3++) {

                                if((magicBullet)||(!distrDEE)||(resInTriple[curPos3])) {

                                    int str3=mutRes2Strand[curPos3];
                                    int strResNum3=strandMut[str3][mutRes2MutIndex[curPos3]];

                                    if( magicBullet )
                                        findMagicBullets(curPos1, curPos2, curPos3, magicBulletAA, magicBulletRot);

                                    for (int AA1=0; AA1<numAAtypes[curPos1]; AA1++) {

                                        int curAA1 = strandRot[str1].getIndexOfNthAllowable(strResNum1,AA1);
                                        numRotForCurAAatPos1 = getNumRot(str1, strResNum1, curAA1);

                                        for(int curRot1=0; curRot1<numRotForCurAAatPos1; curRot1++) {

                                            if (!eliminatedRotAtPos.get(curPos1, curAA1, curRot1)) { //not already pruned

                                                for (int AA2=0; AA2<numAAtypes[curPos2]; AA2++) {
                                                    int curAA2 = strandRot[str2].getIndexOfNthAllowable(strResNum2,AA2);
                                                    int numRotForCurAAatPos2 = getNumRot(str2, strResNum2, curAA2);

                                                    for(int curRot2=0; curRot2<numRotForCurAAatPos2; curRot2++) {

                                                        if ( (!eliminatedRotAtPos.get(curPos2, curAA2, curRot2)) && (!splitFlags[curPos1][curAA1][curRot1][curPos2][curAA2][curRot2]) ) {

                                                            for (int AA3=0; AA3<numAAtypes[curPos3]; AA3++) {
                                                                int curAA3 = strandRot[str3].getIndexOfNthAllowable(strResNum3,AA3);
                                                                int numRotForCurAAatPos3 = getNumRot(str3, strResNum3, curAA3);

                                                                for(int curRot3=0; curRot3<numRotForCurAAatPos3; curRot3++) {

                                                                    if( (!eliminatedRotAtPos.get(curPos3,curAA3,curRot3))
                                                                            && (!splitFlags[curPos2][curAA2][curRot2][curPos3][curAA3][curRot3])
                                                                            && (!splitFlags[curPos1][curAA1][curRot1][curPos3][curAA3][curRot3])
                                                                            && (!tripleFlags[curPos3][curAA3][curRot3][curPos2][curAA2][curRot2][curPos1][curAA1][curRot1]) ) {


                                                                        boolean pruned = false;

                                                                        if( magicBullet ) {
                                                                            for( int mb=0; ( mb<magicBulletNum ) && ( magicBulletRot[mb][0] != -1 ); mb++ ) {

                                                                                if( canEliminateUsing( curPos1, curAA1, curRot1,
                                                                                                       magicBulletAA[mb][0], magicBulletRot[mb][0],
                                                                                                       curPos2, curAA2, curRot2,
                                                                                                       magicBulletAA[mb][1], magicBulletRot[mb][1],
                                                                                                       curPos3, curAA3, curRot3,
                                                                                                       magicBulletAA[mb][2], magicBulletRot[mb][2] ) ) {

                                                                                    pruned = true;
                                                                                    break;
                                                                                }
                                                                            }
                                                                        } else {
                                                                            pruned = CanEliminate(curPos1, curAA1, curRot1, curPos2, curAA2, curRot2,
                                                                                                  curPos3, curAA3, curRot3);
                                                                        }

                                                                        if ( pruned ) {
                                                                            tripleFlags[curPos3][curAA3][curRot3][curPos2][curAA2][curRot2][curPos1][curAA1][curRot1] = true;
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
                                }
                            }
                        }
                    }
                }
                System.out.println("done");
            }

            System.out.println("Number of triples pruned this run: "+prunedCurRun);
            System.out.println("DEE: The minimum difference is "+minDiff);
            System.out.println();

            if ( prunedCurRun==0 ) //no rotamers or pairs pruned this run, so done
                done = true;
            else
                numRuns++;
        }
    }



    //Check if a triple can be eliminated
    private boolean CanEliminate (int posNum1, int AANumAtPos1, int rotNumAtPos1,
                                  int posNum2, int AANumAtPos2, int rotNumAtPos2,
                                  int posNum3, int AANumAtPos3, int rotNumAtPos3) {

        //This should not be called for already-pruned triples (checked in ComputeEliminatedRotConf)
        //So we will not check again but just go right into pruning

        //Get strand numbers and strand residue numbers for the positions being pruned at
        int str1=mutRes2Strand[posNum1];
        int strResNum1=strandMut[str1][mutRes2MutIndex[posNum1]];
        int str2=mutRes2Strand[posNum2];
        int strResNum2=strandMut[str2][mutRes2MutIndex[posNum2]];
        int str3=mutRes2Strand[posNum3];
        int strResNum3=strandMut[str3][mutRes2MutIndex[posNum3]];


        for (int AA1=0; AA1<numAAtypes[posNum1]; AA1++) {

            int altAA1 = strandRot[str1].getIndexOfNthAllowable(strResNum1,AA1);
            int numRotForAAatPos1 = getNumRot(str1, strResNum1, altAA1);

            for (int altRot1=0; altRot1<numRotForAAatPos1; altRot1++) {

                //We allow t and r to be the same rotamer of the same AA

                if (!eliminatedRotAtPos.get(posNum1, altAA1, altRot1)) { //not pruned

                    int numRotForAAatPos2;

                    for (int AA2=0; AA2<numAAtypes[posNum2]; AA2++) {

                        int altAA2 = strandRot[str2].getIndexOfNthAllowable(strResNum2,AA2);
                        numRotForAAatPos2 = getNumRot(str2, strResNum2, altAA2);

                        for (int altRot2=0; altRot2<numRotForAAatPos2; altRot2++) {

                            if ( (!eliminatedRotAtPos.get(posNum2, altAA2, altRot2))
                                    && ( ! splitFlags[posNum1][altAA1][altRot1][posNum2][altAA2][altRot2] ) ) { //not pruned

                                int numRotForAAatPos3;

                                for(int AA3=0; AA3<numAAtypes[posNum3]; AA3++) {

                                    int altAA3 = strandRot[str3].getIndexOfNthAllowable(strResNum3,AA3);
                                    numRotForAAatPos3 = getNumRot(str3, strResNum3, altAA3);

                                    for(int altRot3=0; altRot3<numRotForAAatPos3; altRot3++) {

                                        if ( (!eliminatedRotAtPos.get(posNum3, altAA3, altRot3))
                                                && (!splitFlags[posNum2][altAA2][altRot2][posNum3][altAA3][altRot3])
                                                && (!splitFlags[posNum1][altAA1][altRot1][posNum3][altAA3][altRot3]) ) { //not pruned (except we allow pruned triples to be used because some of them have proved useful)


                                            if( canEliminateUsing(posNum1, AANumAtPos1, rotNumAtPos1, altAA1, altRot1,
                                                                  posNum2, AANumAtPos2, rotNumAtPos2, altAA2, altRot2,
                                                                  posNum3, AANumAtPos3, rotNumAtPos3, altAA3, altRot3) )
                                                return true;
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



    //Check if a triple (indicated by r) can be eliminated using another triple (indicated by t)
    //(residues indicated by 1, 2, 3)
    private boolean canEliminateUsing (int pos1, int r1AA, int r1Rot, int t1AA, int t1Rot,
                                       int pos2, int r2AA, int r2Rot, int t2AA, int t2Rot,
                                       int pos3, int r3AA, int r3Rot, int t3AA, int t3Rot) {



        if( typeDependent ) { //Can only prune a rotamer using another rotamer of the same AA type
            if( ( r1AA != t1AA ) || ( r2AA != t2AA ) || ( r3AA != t3AA ) )
                return false;
        }


        double checkSum = 0;

        //Energies of the individual rotamers (candidate minus competitor)
        checkSum += pairwiseMinEnergyMatrix.getIntraAndShellE(pos1,r1AA,r1Rot)	//1st residue
                    - pairwiseMinEnergyMatrix.getIntraAndShellE(pos1,t1AA,t1Rot)
                    + pairwiseMinEnergyMatrix.getIntraAndShellE(pos2,r2AA,r2Rot)	//2nd residue
                    - pairwiseMinEnergyMatrix.getIntraAndShellE(pos2,t2AA,t2Rot)
                    + pairwiseMinEnergyMatrix.getIntraAndShellE(pos3,r3AA,r3Rot)	//3rd residue
                    - pairwiseMinEnergyMatrix.getIntraAndShellE(pos3,t3AA,t3Rot);


        //Pairwise energies within the triple
        checkSum += pairwiseMinEnergyMatrix.getPairwiseE( pos1, r1AA, r1Rot, pos2, r2AA, r2Rot )
                    - pairwiseMinEnergyMatrix.getPairwiseE( pos1, t1AA, t1Rot, pos2, t2AA, t2Rot )
                    + pairwiseMinEnergyMatrix.getPairwiseE( pos1, r1AA, r1Rot, pos3, r3AA, r3Rot )
                    - pairwiseMinEnergyMatrix.getPairwiseE( pos1, t1AA, t1Rot, pos3, t3AA, t3Rot )
                    + pairwiseMinEnergyMatrix.getPairwiseE( pos2, r2AA, r2Rot, pos3, r3AA, r3Rot )
                    - pairwiseMinEnergyMatrix.getPairwiseE( pos2, t2AA, t2Rot, pos3, t3AA, t3Rot );


        //Pairwise energies between rotamers in the triple and rotamers at other residues
        for(int posj=0; posj<numMutable; posj++) {
            //We are assuming pos1 < pos2 < pos3 as coded above

            if( (posj == pos1) || (posj == pos2) || (posj == pos3) )
                continue;

            int strj=mutRes2Strand[posj];
            int strResNumj=strandMut[strj][mutRes2MutIndex[posj]];


            float minTerm = Float.POSITIVE_INFINITY;//This will be the minimum over the witness rotamers at posj of E(candidate, witness) - E(competitor, witness)

            for(int AAs=0; AAs<numAAtypes[posj]; AAs++) {
                int sAA = strandRot[strj].getIndexOfNthAllowable(strResNumj,AAs);
                int sNumRot = getNumRot(strj, strResNumj, sAA);

                for (int sRot=0; sRot<sNumRot; sRot++) {

                    if( eliminatedRotAtPos.get(posj, sAA, sRot) )//Don't consider pruned rotamers
                        continue;

                    if( splitFlags[pos1][r1AA][r1Rot][posj][sAA][sRot]
                            || splitFlags[pos2][r2AA][r2Rot][posj][sAA][sRot]
                            || splitFlags[pos3][r3AA][r3Rot][posj][sAA][sRot]
                            || isPrunedTriple( pos1, r1AA, r1Rot, pos2, r2AA, r2Rot, posj, sAA, sRot )
                            || isPrunedTriple( pos1, r1AA, r1Rot, pos3, r3AA, r3Rot, posj, sAA, sRot )
                            || isPrunedTriple( pos2, r2AA, r2Rot, pos3, r3AA, r3Rot, posj, sAA, sRot ) )//Don't consider rotamers incompatible with the candidate triple
                        continue;


                    float checkMin = pairwiseMinEnergyMatrix.getPairwiseE( pos1, r1AA, r1Rot, posj, sAA, sRot )
                                     - pairwiseMinEnergyMatrix.getPairwiseE( pos1, t1AA, t1Rot, posj, sAA, sRot )
                                     + pairwiseMinEnergyMatrix.getPairwiseE( pos2, r2AA, r2Rot, posj, sAA, sRot )
                                     - pairwiseMinEnergyMatrix.getPairwiseE( pos2, t2AA, t2Rot, posj, sAA, sRot )
                                     + pairwiseMinEnergyMatrix.getPairwiseE( pos3, r3AA, r3Rot, posj, sAA, sRot )
                                     - pairwiseMinEnergyMatrix.getPairwiseE( pos3, t3AA, t3Rot, posj, sAA, sRot );

                    if ( checkMin < minTerm )
                        minTerm = checkMin;
                }
            }


            if(minTerm == Float.POSITIVE_INFINITY)
                return true;//The candidate triple can be pruned because it is incompatible with all rotamers at posj

            checkSum += minTerm;
        }


        if ( checkSum > curEw)
            return true;//this rotamer can be pruned/eliminated
        else
            minDiff = Math.max(minDiff,checkSum);

        return false;

    }





    private void findMagicBullets(int curPos1, int curPos2, int curPos3, int[][] magicBulletAA, int[][] magicBulletRot) {
        //Find the magicBulletNum "best" triples at the specified position with respect to internal energy plus maximum sum of pairwise energy with other positions
        //Similar to what Gordon and Mayo 1999 do for pairs but with magicBulletNum instead of just 1 magic bullets

        float bestEnergies[] = new float[magicBulletNum];

        for(int a=0; a<magicBulletNum; a++) {
            bestEnergies[a] = Float.POSITIVE_INFINITY;
            for(int b=0; b<3; b++) {
                magicBulletAA[a][b] = -1;
                magicBulletRot[a][b] = -1;
            }
        }


        int str1=mutRes2Strand[curPos1];
        int strResNum1=strandMut[str1][mutRes2MutIndex[curPos1]];
        int str2=mutRes2Strand[curPos2];
        int strResNum2=strandMut[str2][mutRes2MutIndex[curPos2]];
        int str3=mutRes2Strand[curPos3];
        int strResNum3=strandMut[str3][mutRes2MutIndex[curPos3]];


        //The structure of this loop is copied from computeEliminatedRotConf
        for (int AA1=0; AA1<numAAtypes[curPos1]; AA1++) {

            int curAA1 = strandRot[str1].getIndexOfNthAllowable(strResNum1,AA1);
            int numRotForCurAAatPos1 = getNumRot(str1, strResNum1, curAA1);

            for(int curRot1=0; curRot1<numRotForCurAAatPos1; curRot1++) {

                if (!eliminatedRotAtPos.get(curPos1, curAA1, curRot1)) { //not already pruned

                    for (int AA2=0; AA2<numAAtypes[curPos2]; AA2++) {
                        int curAA2 = strandRot[str2].getIndexOfNthAllowable(strResNum2,AA2);
                        int numRotForCurAAatPos2 = getNumRot(str2, strResNum2, curAA2);

                        for(int curRot2=0; curRot2<numRotForCurAAatPos2; curRot2++) {

                            if ( (!eliminatedRotAtPos.get(curPos2, curAA2, curRot2))
                                    && (!splitFlags[curPos1][curAA1][curRot1][curPos2][curAA2][curRot2]) ) {

                                for (int AA3=0; AA3<numAAtypes[curPos3]; AA3++) {
                                    int curAA3 = strandRot[str3].getIndexOfNthAllowable(strResNum3,AA3);
                                    int numRotForCurAAatPos3 = getNumRot(str3, strResNum3, curAA3);


                                    for(int curRot3=0; curRot3<numRotForCurAAatPos3; curRot3++) {

                                        if( (!eliminatedRotAtPos.get(curPos3,curAA3,curRot3))
                                                && (!splitFlags[curPos2][curAA2][curRot2][curPos3][curAA3][curRot3])
                                                && (!splitFlags[curPos1][curAA1][curRot1][curPos3][curAA3][curRot3])
                                                && (!tripleFlags[curPos3][curAA3][curRot3][curPos2][curAA2][curRot2][curPos1][curAA1][curRot1]) ) {
                                            //We don't need to call isPrunedTriples because we know curPos3 > curPos2 > curPos1


                                            //This quantity is similar to the pruning condition but is just being evaluated for the single triple rather than comparing two
                                            float checkSum = 0;

                                            //Energies of the individual rotamers
                                            checkSum += pairwiseMinEnergyMatrix.getIntraAndShellE( curPos1, curAA1, curRot1 ) 	//intra and shell energies: 1st residue
                                                        + pairwiseMinEnergyMatrix.getIntraAndShellE( curPos2, curAA2, curRot2 ) //2nd residue
                                                        + pairwiseMinEnergyMatrix.getIntraAndShellE( curPos3, curAA3, curRot3 ); //3rd residue


                                            //Pairwise energies within the triple
                                            checkSum += pairwiseMinEnergyMatrix.getPairwiseE( curPos1, curAA1, curRot1, curPos2, curAA2, curRot2 )
                                                        + pairwiseMinEnergyMatrix.getPairwiseE( curPos1, curAA1, curRot1, curPos3, curAA3, curRot3 )
                                                        + pairwiseMinEnergyMatrix.getPairwiseE( curPos2, curAA2, curRot2, curPos3, curAA3, curRot3 );


                                            //Pairwise energies between rotamers in the triple and rotamers at other residues
                                            for(int posj=0; posj<numMutable; posj++) {
                                                //We are assuming pos1 < pos2 < pos3 as coded above

                                                if( (posj == curPos1) || (posj == curPos2) || (posj == curPos3) )
                                                    continue;

                                                int strj=mutRes2Strand[posj];
                                                int strResNumj=strandMut[strj][mutRes2MutIndex[posj]];


                                                float maxTerm = Float.NEGATIVE_INFINITY;//This will be the maximum over the witness rotamers at posj of E(triple, witness)

                                                for(int AAs=0; AAs<numAAtypes[posj]; AAs++) {
                                                    int sAA = strandRot[strj].getIndexOfNthAllowable(strResNumj,AAs);
                                                    int sNumRot = getNumRot(strj, strResNumj, sAA);

                                                    for (int sRot=0; sRot<sNumRot; sRot++) {

                                                        if( eliminatedRotAtPos.get(posj, sAA, sRot) )//Don't consider pruned rotamers
                                                            continue;

                                                        if( splitFlags[curPos1][curAA1][curRot1][posj][sAA][sRot]
                                                                || splitFlags[curPos2][curAA2][curRot2][posj][sAA][sRot]
                                                                || splitFlags[curPos3][curAA3][curRot3][posj][sAA][sRot]
                                                                || isPrunedTriple( curPos1, curAA1, curRot1, curPos2, curAA2, curRot2, posj, sAA, sRot )
                                                                || isPrunedTriple( curPos1, curAA1, curRot1, curPos3, curAA3, curRot3, posj, sAA, sRot )
                                                                || isPrunedTriple( curPos2, curAA2, curRot2, curPos3, curAA3, curRot3, posj, sAA, sRot ) )//Don't consider rotamers incompatible with the candidate triple
                                                            continue;

                                                        float checkMax = pairwiseMinEnergyMatrix.getPairwiseE( curPos1, curAA1, curRot1, posj, sAA, sRot )
                                                                         + pairwiseMinEnergyMatrix.getPairwiseE( curPos2, curAA2, curRot2, posj, sAA, sRot )
                                                                         + pairwiseMinEnergyMatrix.getPairwiseE( curPos3, curAA3, curRot3, posj, sAA, sRot );

                                                        if ( ( checkMax > maxTerm ) && ( checkMax != Float.POSITIVE_INFINITY ) )//Parametrically incompatible conformations not counted
                                                            maxTerm = checkMax;
                                                    }
                                                }


                                                if(maxTerm == Float.NEGATIVE_INFINITY) {
                                                    checkSum = Float.POSITIVE_INFINITY;//This triple should not be used as a magic bullet
                                                    break;
                                                }

                                                checkSum += maxTerm;
                                            }


                                            int worst = argmax( bestEnergies );//Index of the highest energy in bestEnergies so far

                                            if( bestEnergies[worst] > checkSum ) {
                                                //If this is the case then checkSum is one of the magicBulletNum best energies so far
                                                //and this triple should go in bestEnergies

                                                bestEnergies[worst] = checkSum;
                                                magicBulletAA[worst][0] = curAA1;
                                                magicBulletAA[worst][1] = curAA2;
                                                magicBulletAA[worst][2] = curAA3;
                                                magicBulletRot[worst][0] = curRot1;
                                                magicBulletRot[worst][1] = curRot2;
                                                magicBulletRot[worst][2] = curRot3;
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




    public void initializeTripleFlags() { //Initializes the tripleFlags array
        //Only initializes for each real, unpruned conformation once
        //(permutations are handled by using descending order of position #)

        System.out.println("Initializing triple flags...");
        tripleFlags = new boolean[numMutable][][][][][][][][];


        for (int curPos1=0; curPos1<numMutable; curPos1++) { //curPos1 and curPos2 are used differently here than in ComputeEliminatedRotConf:
            //curPos1 > curPos2 > curPos3 here to compactly store the triple flags

            int str1=mutRes2Strand[curPos1];
            int strResNum1=strandMut[str1][mutRes2MutIndex[curPos1]];

            tripleFlags[curPos1] = new boolean[strandRot[str1].rl.getNumAAallowed()][][][][][][][];

            for (int AA1=0; AA1<numAAtypes[curPos1]; AA1++) {
                int curAA1 = strandRot[str1].getIndexOfNthAllowable(strResNum1,AA1);
                int numRotForCurAAatPos1 = getNumRot(str1, strResNum1, curAA1);

                tripleFlags[curPos1][curAA1] = new boolean[numRotForCurAAatPos1][][][][][][];

                for(int curRot1=0; curRot1<numRotForCurAAatPos1; curRot1++) {

                    if(!eliminatedRotAtPos.get(curPos1, curAA1, curRot1)) {

                        tripleFlags[curPos1][curAA1][curRot1] = new boolean[curPos1][][][][][];

                        for(int curPos2=0; curPos2<curPos1; curPos2++) {

                            int str2=mutRes2Strand[curPos2];
                            int strResNum2=strandMut[str2][mutRes2MutIndex[curPos2]];

                            tripleFlags[curPos1][curAA1][curRot1][curPos2] = new boolean[strandRot[str2].rl.getNumAAallowed()][][][][];

                            for (int AA2=0; AA2<numAAtypes[curPos2]; AA2++) {
                                int curAA2 = strandRot[str2].getIndexOfNthAllowable(strResNum2,AA2);
                                int numRotForCurAAatPos2 = getNumRot(str2, strResNum2, curAA2);

                                tripleFlags[curPos1][curAA1][curRot1][curPos2][curAA2] = new boolean[numRotForCurAAatPos2][][][];

                                for(int curRot2=0; curRot2<numRotForCurAAatPos2; curRot2++) {

                                    if(!eliminatedRotAtPos.get(curPos2, curAA2, curRot2)) {

                                        tripleFlags[curPos1][curAA1][curRot1][curPos2][curAA2][curRot2] = new boolean[curPos2][][];

                                        for(int curPos3=0; curPos3<curPos2; curPos3++) {

                                            int str3=mutRes2Strand[curPos3];
                                            int strResNum3=strandMut[str3][mutRes2MutIndex[curPos3]];

                                            tripleFlags[curPos1][curAA1][curRot1][curPos2][curAA2][curRot2][curPos3] =
                                                new boolean[strandRot[str3].rl.getNumAAallowed()][];

                                            for (int AA3=0; AA3<numAAtypes[curPos3]; AA3++) {
                                                int curAA3 = strandRot[str3].getIndexOfNthAllowable(strResNum3,AA3);
                                                int numRotForCurAAatPos3 = getNumRot(str3, strResNum3, curAA3);

                                                tripleFlags[curPos1][curAA1][curRot1][curPos2][curAA2][curRot2][curPos3][curAA3] = new boolean[numRotForCurAAatPos3];
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



    private int argmax(float[] x) {
        //Return the index of the biggest element in x
        //(If it's a tie return the first of the tied indices)

        float max = Float.NEGATIVE_INFINITY;
        int argMax = 0;

        for(int a=0; a<x.length; a++) {
            if(x[a] > max) {
                argMax = a;
                max = x[a];
            }
        }

        return argMax;
    }

}