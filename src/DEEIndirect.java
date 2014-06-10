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
//	DEEIndirect.java
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
 * Performs indirect DEE rotamer and pairs pruning
 * considering the entire protein as the system being pruned.
 * Uses rigid DEE or iMinDEE
*/


public class DEEIndirect extends DEE {

    //useFlags and useMinDEEPruningEw are not used
    //because we always use flags and iMinDEE is assumed if doing minimization

    private float pairK[][][][][];//Precomputed K(i_r,i_t) matrix
    boolean inZ[];//Indicates which of the flexible residue positions are in the pruning zone

    int prunedSinglesCurRun;
    int prunedPairsCurRun;

    //constructor
    DEEIndirect (PairwiseEnergyMatrix arpMatrix, int numResMutable, int strMut[][], float initEw,
                 StrandRotamers strandLRot[], PrunedRotamers<Boolean> prunedRotAtRes, boolean residueMut[],
                 boolean doMin, boolean spFlags[][][][][][], boolean mb, boolean dDEE, boolean minBB,
                 int mutRes2StrandP[], int mutRes2MutIndexP[], boolean typeDep, boolean aIMinDEE,
                 float aIval, boolean tripFlags[][][][][][][][][], boolean doPerts, boolean inZ[]) {

        init(arpMatrix, null, numResMutable,
             strMut, initEw, strandLRot, prunedRotAtRes, doMin, null, null,
             spFlags, true, minBB, mutRes2StrandP, mutRes2MutIndexP, typeDep, aIMinDEE, aIval,
             mb, dDEE, residueMut, tripFlags, doPerts);

        this.inZ = inZ;


        if(doMinimize && (!doIMinDEE) ) {
            System.err.println("ERROR: iMinDEE required for indirect pruning and thus for algOption >= 4");
            System.exit(1);
        }

        //Now allocate space for pairK
        pairK = new float[numMutable][][][][];

        for(int pos=0; pos<numMutable; pos++) {

            if(inZ[pos]) {

                int str=mutRes2Strand[pos];
                int strResNum=strandMut[str][mutRes2MutIndex[pos]];

                pairK[pos] = new float[strandRot[str].rl.getNumAAallowed()][][][];

                for(int AA=0; AA<numAAtypes[pos]; AA++) {
                    int curAA = strandRot[str].getIndexOfNthAllowable(strResNum,AA);
                    int numRot = getNumRot(str, strResNum, curAA);
                    pairK[pos][curAA] = new float[numRot][][];

                    for(int rot=0; rot<numRot; rot++) {
                        pairK[pos][curAA][rot] = new float[strandRot[str].rl.getNumAAallowed()][];

                        for(int AA2=0; AA2<numAAtypes[pos]; AA2++) {
                            int curAA2 = strandRot[str].getIndexOfNthAllowable(strResNum,AA2);
                            int numRot2 = getNumRot(str, strResNum, curAA2);
                            pairK[pos][curAA][rot][curAA2] = new float[numRot2];
                        }
                    }
                }
            }
        }
    }


    //Compute the conformations that can be eliminated
    //Return a boolean array in which an element is true if
    //the corresponding r at i can be eliminated, and false otherwise
    //Also saves "true" to splitFlags for pairs that can be eliminated
    public PrunedRotamers<Boolean> ComputeEliminatedRotConf() {

        boolean done = false;
        numRuns = 1;

        while (!done) {

            prunedSinglesCurRun = 0;
            prunedPairsCurRun = 0;

            System.out.println("Current run: "+numRuns);

            //First precompute K(i_r,i_t) for each rotamer pair in the pruning zone
            for (int curPos=0; curPos<numMutable; curPos++) {

                if( inZ[curPos] ) {

                    System.out.print("Starting precomputation: residue " + curPos);

                    int str=mutRes2Strand[curPos];
                    int strResNum=strandMut[str][mutRes2MutIndex[curPos]];

                    //There are two types of indices for amino acids:
                    //those among the allowed AAs at some residue
                    //(variable names here for this are "AA" or "AA[extra letter or number]")
                    //and those among all the 20 possible AAs
                    //(used in energy matrix, pairK, etc; names are like curAA or AANumAtPos)

                    for (int AA=0; AA<numAAtypes[curPos]; AA++) {

                        System.out.print(".");

                        int curAA = strandRot[str].getIndexOfNthAllowable(strResNum,AA);
                        int numRotForCurAAatPos = getNumRot(str, strResNum, curAA);

                        for(int curRot=0; curRot<numRotForCurAAatPos; curRot++) { //This is i_r

                            if ( !eliminatedRotAtPos.get(curPos, curAA, curRot) ) { //not already pruned

                                //Find i_t now
                                for(int AAt=0; AAt<numAAtypes[curPos]; AAt++ ) {

                                    int tAA = strandRot[str].getIndexOfNthAllowable(strResNum,AAt);
                                    int tNumRot = getNumRot(str, strResNum, tAA);

                                    for(int tRot=0; tRot<tNumRot; tRot++) {
                                        if(!eliminatedRotAtPos.get(curPos, tAA, tRot))
                                            pairK[curPos][curAA][curRot][tAA][tRot] = precomputePairK(curPos, curAA, curRot, tAA, tRot);
                                    }
                                }
                            }
                        }
                    }
                    System.out.println("done");
                }
            }


            //Now try to prune singles.
            for (int curPos=0; curPos<numMutable; curPos++) {

                if( inZ[curPos] ) {

                    System.out.print("Starting residue "+curPos);

                    int str=mutRes2Strand[curPos];
                    int strResNum=strandMut[str][mutRes2MutIndex[curPos]];

                    for (int AA=0; AA<numAAtypes[curPos]; AA++) {

                        System.out.print(".");

                        int curAA = strandRot[str].getIndexOfNthAllowable(strResNum,AA);
                        int numRotForCurAAatPos = getNumRot(str, strResNum, curAA);

                        for(int curRot=0; curRot<numRotForCurAAatPos; curRot++) {

                            if (!eliminatedRotAtPos.get(curPos, curAA, curRot)) { //not already pruned

                                if (CanEliminate(curPos, curAA, curRot)) {
                                    eliminatedRotAtPos.set(curPos, curAA, curRot, true);
                                    prunedSinglesCurRun++;
                                } else
                                    eliminatedRotAtPos.set(curPos, curAA, curRot, false);
                            }
                        }
                    }
                    System.out.println("done");

                }
            }


            //Now try to prune pairs

            for (int curPos1=0; curPos1<numMutable; curPos1++) {

                if ( inZ[curPos1] && ( (magicBullet)||(!distrDEE)||(resInPair[curPos1]) ) ) { //mb-pairs or not distrDEE or cur res is in distr pair; also, this res must be in the pruning zone


                    System.out.print("Starting residue "+curPos1);
                    System.out.print("..");

                    int str1=mutRes2Strand[curPos1];
                    int strResNum1=strandMut[str1][mutRes2MutIndex[curPos1]];

                    for (int curPos2=curPos1+1; curPos2<numMutable; curPos2++) {

                        if ( inZ[curPos2] && ( (magicBullet)||(!distrDEE)||(resInPair[curPos2]) ) ) { //mb-pairs or not distrDEE or cur res is in distr pair

                            int str2=mutRes2Strand[curPos2];
                            int strResNum2=strandMut[str2][mutRes2MutIndex[curPos2]];

                            for (int AA1=0; AA1<numAAtypes[curPos1]; AA1++) {

                                int curAA1 = strandRot[str1].getIndexOfNthAllowable(strResNum1,AA1);
                                int numRotForCurAAatPos1 = getNumRot(str1, strResNum1, curAA1);

                                for(int curRot1=0; curRot1<numRotForCurAAatPos1; curRot1++) {

                                    if (!eliminatedRotAtPos.get(curPos1, curAA1, curRot1)) { //not already pruned

                                        for (int AA2=0; AA2<numAAtypes[curPos2]; AA2++) {
                                            int curAA2 = strandRot[str2].getIndexOfNthAllowable(strResNum2,AA2);
                                            int numRotForCurAAatPos2 = getNumRot(str2, strResNum2, curAA2);

                                            for(int curRot2=0; curRot2<numRotForCurAAatPos2; curRot2++) {

                                                if (!eliminatedRotAtPos.get(curPos2, curAA2, curRot2)) { //not already pruned

                                                    if (!splitFlags[curPos1][curAA1][curRot1][curPos2][curAA2][curRot2]) { //rotamer pair not already pruned

                                                        if (CanEliminate(curPos1, curAA1, curRot1, curPos2, curAA2, curRot2)) {
                                                            splitFlags[curPos1][curAA1][curRot1][curPos2][curAA2][curRot2] = true;
                                                            splitFlags[curPos2][curAA2][curRot2][curPos1][curAA1][curRot1] = true;
                                                            prunedPairsCurRun++;
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
            }

            //Higher-order indirect pruning would go here

            System.out.println("Number of rotamers pruned this run: "+prunedSinglesCurRun);
            System.out.println("Number of pairs pruned this run: "+prunedPairsCurRun);
            System.out.println("DEE: The minimum difference is "+minDiff);
            System.out.println();

            if ( prunedSinglesCurRun==0 && prunedPairsCurRun==0 ) //no rotamers or pairs pruned this run, so done
                done = true;
            else
                numRuns++;
        }

        return eliminatedRotAtPos;
    }

    //Called only by ComputeEliminatedRotConf(.)
    private boolean CanEliminate (int posNum, int AANumAtPos, int rotNumAtPos) {

        double checkSum;

        int str=mutRes2Strand[posNum];
        int strResNum=strandMut[str][mutRes2MutIndex[posNum]];


        //In the energy matrix, column 0 gives the individual energies for each r at i;
        //skip row 0, as the individual energies start from row 1 (and are all in column 0)

        if ((!eliminatedRotAtPos.get(posNum, AANumAtPos, rotNumAtPos))) { //not pruned yet

            //For the particular position, compare the energy performance (one by one)
            //of the remaining rotamer possibilities to that of the given rotamer:
            //given r at i, compare it to all t at i for pruning

            for (int AA=0; AA<numAAtypes[posNum]; AA++) {

                int altAA = strandRot[str].getIndexOfNthAllowable(strResNum,AA);

                if( (!typeDependent) || (altAA==AANumAtPos) ) {

                    int numRotForAAatPos = getNumRot(str, strResNum, altAA);

                    for (int altRot=0; altRot<numRotForAAatPos; altRot++) {

                        //if t and r are not actually the same rotamer of the same AA
                        if (!((altAA==AANumAtPos)&&(altRot==rotNumAtPos))) {

                            //at this point, we know what r at i and t at i are

                            //if ((maxIndVoxelE<=stericEThreshIntra)&&(maxShellResE<=stericEThreshPair)){//check only if not an unallowed steric
                            if (!eliminatedRotAtPos.get(posNum, altAA, altRot)) { //not pruned

                                checkSum = 0;

                                int s_AAtypes[] = new int[numMutable];//AA types and rotamers used to prune
                                int s_rots[] = new int[numMutable];

                                for(int pos2=0; pos2<numMutable; pos2++) {

                                    if(pos2 == posNum) {
                                        checkSum += pairK[posNum][AANumAtPos][rotNumAtPos][altAA][altRot];
                                        continue;
                                    } else if( inZ[pos2] ) {

                                        int str2=mutRes2Strand[pos2];
                                        int strResNum2=strandMut[str2][mutRes2MutIndex[pos2]];

                                        float Kmaxmin = Float.NEGATIVE_INFINITY;

                                        for(int AAs=0; AAs<numAAtypes[pos2]; AAs++) {
                                            int sAA = strandRot[str2].getIndexOfNthAllowable(strResNum2, AAs);
                                            int sNumRot = getNumRot(str2, strResNum2, sAA);
                                            for(int sRot=0; sRot<sNumRot; sRot++) {

                                                if( eliminatedRotAtPos.get(pos2, sAA, sRot)
                                                        || splitFlags[pos2][sAA][sRot][posNum][altAA][altRot] )
                                                    continue;//Don't use this j_s

                                                float sKmin = Float.POSITIVE_INFINITY;

                                                for(int AAu=0; AAu<numAAtypes[pos2]; AAu++) {

                                                    int uAA = strandRot[str2].getIndexOfNthAllowable(strResNum2, AAu);
                                                    int uNumRot = getNumRot(str2, strResNum2, uAA);

                                                    for(int uRot=0; uRot<uNumRot; uRot++) {

                                                        if( eliminatedRotAtPos.get(pos2, uAA, uRot)
                                                                || splitFlags[pos2][uAA][uRot][posNum][AANumAtPos][rotNumAtPos] )
                                                            continue;

                                                        if( pairK[pos2][uAA][uRot][sAA][sRot] < sKmin )
                                                            sKmin = pairK[pos2][uAA][uRot][sAA][sRot];
                                                    }
                                                }

                                                if(sKmin == Float.POSITIVE_INFINITY)
                                                    return true;//Rotamer r can be pruned because it is incompatible with all unpruned rotamers at pos2

                                                if(sKmin > Kmaxmin) {
                                                    boolean isCompatible = true;//is index_s compatible with the rest of the s-indices used?
                                                    for(int a=0; a<pos2; a++) {
                                                        if( inZ[a] && (a != posNum) && splitFlags[a][s_AAtypes[a]][s_rots[a]][pos2][sAA][sRot] )
                                                            isCompatible = false;
                                                    }
                                                    if(isCompatible) {
                                                        Kmaxmin = sKmin;
                                                        s_AAtypes[pos2] = sAA;
                                                        s_rots[pos2] = sRot;
                                                    }
                                                }

                                            }
                                        }

                                        if(Kmaxmin == Float.NEGATIVE_INFINITY) {
                                            checkSum = Float.NEGATIVE_INFINITY;
                                            break;//r will not be pruned using t
                                        }

                                        checkSum += Kmaxmin;//add max min K(j_u,j_s) to checkSum where j_u is compatible with i_r and j_s with i_t
                                    }
                                }

                                if ( checkSum > curEw)
                                    return true;//this rotamer can be pruned/eliminated
                                else
                                    minDiff = Math.max(minDiff,checkSum);


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



    //Check if a pair can be eliminated
    private boolean CanEliminate (int posNum1, int AANumAtPos1, int rotNumAtPos1, int posNum2,
                                  int AANumAtPos2, int rotNumAtPos2) {


        int str1=mutRes2Strand[posNum1];
        int strResNum1=strandMut[str1][mutRes2MutIndex[posNum1]];
        int str2=mutRes2Strand[posNum2];
        int strResNum2=strandMut[str2][mutRes2MutIndex[posNum2]];


        if( eliminatedRotAtPos.get(posNum1, AANumAtPos1, rotNumAtPos1)
                || eliminatedRotAtPos.get(posNum2, AANumAtPos2, rotNumAtPos2)
                || splitFlags[posNum1][AANumAtPos1][rotNumAtPos1][posNum2][AANumAtPos2][rotNumAtPos2] )
            return true;//Already pruned


        for (int AA1=0; AA1<numAAtypes[posNum1]; AA1++) {

            int altAA1 = strandRot[str1].getIndexOfNthAllowable(strResNum1,AA1);
            int numRotForAAatPos1 = getNumRot(str1, strResNum1, altAA1);

            for (int altRot1=0; altRot1<numRotForAAatPos1; altRot1++) {

                //if t and r are not actually the same rotamer of the same AA
                //if (!((altAA1==AANumAtPos1)&&(altRot1==rotNumAtPos1))){

                if (!eliminatedRotAtPos.get(posNum1, altAA1, altRot1)) { //not pruned

                    for (int AA2=0; AA2<numAAtypes[posNum2]; AA2++) {

                        int altAA2 = strandRot[str2].getIndexOfNthAllowable(strResNum2,AA2);
                        int numRotForAAatPos2 = getNumRot(str2, strResNum2, altAA2);

                        if( (!typeDependent) || ( (altAA1==AANumAtPos1) && (altAA2==AANumAtPos2) ) ) {

                            for (int altRot2=0; altRot2<numRotForAAatPos2; altRot2++) {

                                //if t and r are not actually the same rotamer of the same AA
                                //if (!((altAA2==AANumAtPos2)&&(altRot2==rotNumAtPos2))){

                                if ( (!eliminatedRotAtPos.get(posNum2 ,altAA2, altRot2))
                                        && (!splitFlags[posNum1][altAA1][altRot1][posNum2][altAA2][altRot2]) ) { //not pruned

                                    if( canEliminateUsing( posNum1, AANumAtPos1, rotNumAtPos1, altAA1, altRot1,
                                                           posNum2, AANumAtPos2, rotNumAtPos2, altAA2, altRot2 ) )
                                        return true;
                                    else if (magicBullet) //magic bullet pairs, so no further checks
                                        return false;
                                }
                                //}
                            }
                        }
                    }
                }
                //}
            }
        }

        //We have tried all of the other rotamers at the current position and none
        //of them is able to prune the given rotamer, so we return false
        return false;
    }



    //Check if a pair (indicated by r) can be eliminated using another pair (indicated by t)
    //(residues indicated by 1 and 2)
    private boolean canEliminateUsing (int pos1, int r1AA, int r1Rot, int t1AA, int t1Rot,
                                       int pos2, int r2AA, int r2Rot, int t2AA, int t2Rot) {

        double checkSum = 0;

        int s_AAtypes[] = new int[numMutable];//AA types and rotamers used to prune (in the Kmaxmin terms)
        int s_rots[] = new int[numMutable];

        //First add up K's at other positions (these are stored in pairK)
        for(int posj=0; posj<numMutable; posj++) {

            if( inZ[posj] && (posj != pos1) && (posj != pos2) ) {

                int strj=mutRes2Strand[posj];
                int strResNumj=strandMut[strj][mutRes2MutIndex[posj]];


                float Kmaxmin = Float.NEGATIVE_INFINITY;

                for(int AAs=0; AAs<numAAtypes[posj]; AAs++) {
                    int sAA = strandRot[strj].getIndexOfNthAllowable(strResNumj, AAs);
                    int sNumRot = getNumRot(strj, strResNumj, sAA);

                    for(int sRot=0; sRot<sNumRot; sRot++) {

                        if( eliminatedRotAtPos.get(posj, sAA, sRot)
                                || splitFlags[posj][sAA][sRot][pos1][t1AA][t1Rot]
                                || splitFlags[posj][sAA][sRot][pos2][t2AA][t2Rot] )
                            continue;//Don't use this j_s

                        if( isPrunedTriple(posj,sAA,sRot,pos1,t1AA,t1Rot,pos2,t2AA,t2Rot) )
                            continue;
                        //This is checked after lower-order pruning to make sure we don't call this function for pruned rotamers


                        float sKmin = Float.POSITIVE_INFINITY;

                        for(int AAu=0; AAu<numAAtypes[posj]; AAu++) {
                            int uAA = strandRot[strj].getIndexOfNthAllowable(strResNumj, AAu);
                            int uNumRot = getNumRot(strj, strResNumj, uAA);

                            for(int uRot=0; uRot<uNumRot; uRot++) {

                                if( eliminatedRotAtPos.get(posj, uAA, uRot)
                                        || splitFlags[posj][uAA][uRot][pos1][r1AA][r1Rot]
                                        || splitFlags[posj][uAA][uRot][pos2][r2AA][r2Rot] )
                                    continue;

                                if( isPrunedTriple(posj,uAA,uRot,pos1,r1AA,r1Rot,pos2,r2AA,r2Rot) )
                                    continue;

                                if( pairK[posj][uAA][uRot][sAA][sRot] < sKmin )
                                    sKmin = pairK[posj][uAA][uRot][sAA][sRot];
                            }
                        }

                        if ( sKmin == Float.POSITIVE_INFINITY )
                            return true;
                        //In this case, there is no rotamer u at posj that is compatible with both r1 at pos1 and r2 at pos1
                        //So the rotamer pair (r1 at pos1, r2 at pos2) is impossible


                        if(sKmin > Kmaxmin) {
                            boolean isCompatible = true;//is index_s compatible with the rest of the s-indices used?
                            for(int a=0; a<posj; a++) {
                                if( inZ[a] && (a != pos1) && (a != pos2) && splitFlags[a][s_AAtypes[a]][s_rots[a]][posj][sAA][sRot] )
                                    isCompatible = false;
                            }
                            if(isCompatible) {
                                Kmaxmin = sKmin;
                                s_AAtypes[posj] = sAA;
                                s_rots[posj] = sRot;
                            }
                        }
                    }
                }

                if(Kmaxmin == Float.NEGATIVE_INFINITY)
                    return false;
                //This messes up the pruning of (r1 at pos1, r2 at pos2)
                //So go try to find some other competitor pair

                checkSum += Kmaxmin;//add max min K(j_u,j_s) to checkSum where j_u is compatible with i_r and j_s with i_t
            }
        }



        //Now add in K(i_r,i_t) (bold i,r,t)
        checkSum += pairwiseMinEnergyMatrix.getIntraAndShellE( pos1, r1AA, r1Rot ) 	//intra and shell energies: 1st residue
                    - pairwiseMinEnergyMatrix.getIntraAndShellE( pos1, t1AA, t1Rot )
                    + pairwiseMinEnergyMatrix.getIntraAndShellE( pos2, r2AA, r2Rot ) //2nd residue
                    - pairwiseMinEnergyMatrix.getIntraAndShellE( pos2, t2AA, t2Rot );


        //Pairwise energy within the pair
        checkSum +=  pairwiseMinEnergyMatrix.getPairwiseE( pos1, r1AA, r1Rot, pos2, r2AA, r2Rot )
                     - pairwiseMinEnergyMatrix.getPairwiseE( pos1, t1AA, t1Rot, pos2, t2AA, t2Rot );


        //Now sum the other pairwise energy differences within the pruning zone
        for(int posj=0; posj<pos2; posj++) {
            //We are assuming pos2 > pos1

            if( inZ[posj] && (posj != pos1) ) {

                int strj=mutRes2Strand[posj];
                int strResNumj=strandMut[strj][mutRes2MutIndex[posj]];

                float minTerm = Float.POSITIVE_INFINITY;//This will be the minimum E(i_r,j_s) over j_s in R_(j,i_r)
                float maxTerm = Float.NEGATIVE_INFINITY;//This will be the maximum E(i_t,j_s) over j_s in R(j,i_t)

                for(int AAs=0; AAs<numAAtypes[posj]; AAs++) {
                    int sAA = strandRot[strj].getIndexOfNthAllowable(strResNumj, AAs);
                    int sNumRot = getNumRot(strj, strResNumj, sAA);

                    for (int sRot=0; sRot<sNumRot; sRot++) {

                        if( eliminatedRotAtPos.get(posj, sAA, sRot) )//Don't consider pruned rotamers
                            continue;

                        if( ! ( splitFlags[pos1][r1AA][r1Rot][posj][sAA][sRot]
                                || splitFlags[pos2][r2AA][r2Rot][posj][sAA][sRot]
                                || isPrunedTriple(pos1,r1AA,r1Rot,pos2,r2AA,r2Rot,posj,sAA,sRot) ) ) {
                            float checkMin = pairwiseMinEnergyMatrix.getPairwiseE( pos2, r2AA, r2Rot, posj, sAA, sRot );
                            if(posj<pos1)
                                checkMin += pairwiseMinEnergyMatrix.getPairwiseE( pos1, r1AA, r1Rot, posj, sAA, sRot );
                            if ( checkMin < minTerm )
                                minTerm = checkMin;
                        }
                        if( ! ( splitFlags[pos1][t1AA][t1Rot][posj][sAA][sRot] || splitFlags[pos2][t2AA][t2Rot][posj][sAA][sRot]
                                || isPrunedTriple(pos1,t1AA,t1Rot,pos2,t2AA,t2Rot,posj,sAA,sRot) ) ) {
                            float checkMax = pairwiseMinEnergyMatrix.getPairwiseE( pos2, t2AA, t2Rot, posj, sAA, sRot );
                            if(posj<pos1)
                                checkMax += pairwiseMinEnergyMatrix.getPairwiseE( pos1, t1AA, t1Rot, posj, sAA, sRot );
                            if ( checkMax > maxTerm )
                                maxTerm = checkMax;
                        }
                    }
                }


                if(minTerm == Float.POSITIVE_INFINITY)
                    return true;//(r1, r2) can be pruned because it is incompatible with all rotamers at posj

                if(maxTerm == Float.NEGATIVE_INFINITY) {
                    splitFlags[pos1][t1AA][t1Rot][pos2][t2AA][t2Rot] = true;
                    splitFlags[pos2][t2AA][t2Rot][pos1][t1AA][t1Rot] = true;
                    prunedPairsCurRun++;
                    return false;
                }
                //(t1, t2) can be pruned because it is incompatible with all rotamers at posj
                //So, as above, we cannot meaningfully continue trying to prune (r1,r2) with it



                checkSum += minTerm - maxTerm;
            }
        }


        //Finally add terms for residues outside the pruning zone
        for(int posj=0; posj<numMutable; posj++) {

            if ( !inZ[posj] ) {//posj is not in the pruning zone

                int strj=mutRes2Strand[posj];
                int strResNumj=strandMut[strj][mutRes2MutIndex[posj]];

                float minTerm = Float.POSITIVE_INFINITY;//This will be the minimum E(i_r,j_s) - E(i_t,j_s)  over j_s in R_j

                for(int AAs=0; AAs<numAAtypes[posj]; AAs++) {

                    int sAA = strandRot[strj].getIndexOfNthAllowable(strResNumj, AAs);
                    int sNumRot = getNumRot(strj, strResNumj, sAA);

                    for (int sRot=0; sRot<sNumRot; sRot++) {

                        if( ! ( eliminatedRotAtPos.get(posj, sAA, sRot)
                                || splitFlags[pos1][r1AA][r1Rot][posj][sAA][sRot]
                                || splitFlags[pos2][r2AA][r2Rot][posj][sAA][sRot] ) ) {

                            if( ! isPrunedTriple(pos1,r1AA,r1Rot,pos2,r2AA,r2Rot,posj,sAA,sRot) ) {

                                float checkMin = + pairwiseMinEnergyMatrix.getPairwiseE( pos1, r1AA, r1Rot, posj, sAA, sRot )
                                                 + pairwiseMinEnergyMatrix.getPairwiseE( pos2, r2AA, r2Rot, posj, sAA, sRot )
                                                 - pairwiseMinEnergyMatrix.getPairwiseE( pos1, t1AA, t1Rot, posj, sAA, sRot )
                                                 - pairwiseMinEnergyMatrix.getPairwiseE( pos2, t2AA, t2Rot, posj, sAA, sRot );

                                if ( checkMin < minTerm )
                                    minTerm = checkMin;
                            }
                        }
                    }
                }


                if(minTerm == Float.POSITIVE_INFINITY)
                    return true;//(r1, r2) can be pruned because it is incompatible with all rotamers at posj

                checkSum += minTerm;
            }
        }


        if ( checkSum > curEw)
            return true;//this rotamer can be pruned/eliminated
        else
            minDiff = Math.max(minDiff,checkSum);

        return false;

    }


    private float precomputePairK(int curPos, int curAA, int curRot, int tAA, int tRot) {
        //Computes K(i_r,i_t); i = curPos; curAA & curRot give r; AAt and tRot give t
        //This should only be called if i_r and i_t have not been pruned

        //First set K = E(i_r) - E(i_t)
        float K = pairwiseMinEnergyMatrix.getIntraAndShellE( curPos, curAA, curRot ) 	//intra energy for i_r
                  - pairwiseMinEnergyMatrix.getIntraAndShellE( curPos, tAA, tRot );   // i_t now

        //Now sum pairwise energy differences within the pruning zone
        for(int pos2=0; pos2<curPos; pos2++) {

            if( inZ[pos2] ) {

                int str2=mutRes2Strand[pos2];
                int strResNum2=strandMut[str2][mutRes2MutIndex[pos2]];

                float minTerm = Float.POSITIVE_INFINITY;//This will be the minimum E(i_r,j_s) over j_s in R_(j,i_r)
                float maxTerm = Float.NEGATIVE_INFINITY;//This will be the maximum E(i_t,j_s) over j_s in R(j,i_t)

                for(int AAs=0; AAs<numAAtypes[pos2]; AAs++) {
                    int sAA = strandRot[str2].getIndexOfNthAllowable(strResNum2, AAs);
                    int sNumRot = getNumRot(str2, strResNum2, sAA);

                    for (int sRot=0; sRot<sNumRot; sRot++) {

                        if( eliminatedRotAtPos.get(pos2, sAA, sRot) )//Don't consider pruned rotamers
                            continue;

                        if( ! splitFlags[curPos][curAA][curRot][pos2][sAA][sRot] ) {
                            if ( pairwiseMinEnergyMatrix.getPairwiseE( curPos, curAA, curRot, pos2, sAA, sRot ) < minTerm )
                                minTerm = pairwiseMinEnergyMatrix.getPairwiseE( curPos, curAA, curRot, pos2, sAA, sRot );
                        }
                        if( ! splitFlags[curPos][tAA][tRot][pos2][sAA][sRot] ) {
                            if ( pairwiseMinEnergyMatrix.getPairwiseE( curPos, tAA, tRot, pos2, sAA, sRot ) > maxTerm )
                                maxTerm = pairwiseMinEnergyMatrix.getPairwiseE( curPos, tAA, tRot, pos2, sAA, sRot );
                        }
                    }
                }

                if( minTerm == Float.POSITIVE_INFINITY || maxTerm == Float.NEGATIVE_INFINITY ) {

                    if(minTerm == Float.POSITIVE_INFINITY) {
                        eliminatedRotAtPos.set(curPos,curAA,curRot,true);//ir can be pruned because it is incompatible with all rotamers at position j (pos2)
                        prunedSinglesCurRun++;
                    }

                    if(maxTerm == Float.NEGATIVE_INFINITY) {
                        if(!eliminatedRotAtPos.get(curPos,tAA,tRot)) {
                            eliminatedRotAtPos.set(curPos,tAA,tRot,true);//it can be pruned (same argument)
                            prunedSinglesCurRun++;
                        }
                    }

                    return Float.NaN;//Either way, the pruning makes K(ir,it) useless
                }



                K += minTerm - maxTerm;
            }
        }

        //and finally terms from outside the pruning zone
        for(int pos2=0; pos2<curPos; pos2++) {

            if( !inZ[pos2] ) { //pos2 is not in the pruning zone

                int str2=mutRes2Strand[pos2];
                int strResNum2=strandMut[str2][mutRes2MutIndex[pos2]];

                float minTerm = Float.POSITIVE_INFINITY;//This will be the minimum E(i_r,j_s) - E(i_t,j_s)  over j_s in R_j

                for(int AAs=0; AAs<numAAtypes[pos2]; AAs++) {

                    int sAA = strandRot[str2].getIndexOfNthAllowable(strResNum2, AAs);
                    int sNumRot = getNumRot(str2, strResNum2, sAA);

                    for (int sRot=0; sRot<sNumRot; sRot++) {

                        if( ! ( eliminatedRotAtPos.get(pos2, sAA, sRot)
                                || splitFlags[curPos][curAA][curRot][pos2][sAA][sRot] ) ) {

                            float checkMin = pairwiseMinEnergyMatrix.getPairwiseE( curPos, curAA, curRot, pos2, sAA, sRot )
                                             - pairwiseMinEnergyMatrix.getPairwiseE( curPos, tAA, tRot, pos2, sAA, sRot );

                            if ( checkMin < minTerm )
                                minTerm = checkMin;
                        }
                    }
                }


                if(minTerm == Float.POSITIVE_INFINITY) {
                    eliminatedRotAtPos.set(curPos, curAA, curRot, true);//ir can be pruned because it is incompatible with all rotamers at position j (pos2)
                    prunedSinglesCurRun++;
                    return Float.NaN;
                }

                K += minTerm;
            }
        }

        return K;
    }

}