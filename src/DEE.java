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
//	DEE.java
//
//	Version:           2.1 beta
//
//
//	  authors:
// 	  initials    name                 organization                email
//	 ---------   -----------------    ------------------------    ----------------------------
//	  MAH           Mark A. Hallen	  Duke University               mah43@duke.edu
///////////////////////////////////////////////////////////////////////////////////////////////

//This is a base class for the various DEE classes.  It contains shared functionality for different DEE methods
public abstract class DEE {

    //two pairwise energy matrices: one for the min energies and one for the max
	protected PairwiseEnergyMatrix pairwiseMinEnergyMatrix = null;
	protected PairwiseEnergyMatrix pairwiseMaxEnergyMatrix = null;

	//eliminated rotamers at position i, for all positions
	protected PrunedRotamers<Boolean> eliminatedRotAtPos = null;


        //number of positions (active site residues or ligand) to consider
        //protected int numPos;


	//the number of AA types allowed for each AS residue
	int numAAtypes[] = null;

	//number of rotamers for the current AA type at the given residue
	//int numRotForAAtypeAtRes[];

	//this value depends on the particular value specified in the pairwise energy matrices;
	//		in KSParser, this value is 10^38;
	//entries with this particular value will not be examined, as they are not allowed;
	//note that when computing E intervals, if a steric is not allowed, (maxE-minE)=0,
	//		so no comparison with stericE is necessary there
	protected float bigE = (float)Math.pow(10,38);

	//steric energy that determines incompatibility of a rotamer with the template
	float stericE = bigE;

	protected float curEw = 0.0f;	//the max allowable difference from the GMEC (checkSum<=curEw should not be pruned)

	//the minimum difference in the checkSum when a rotamer cannot be pruned
	protected double minDiff = -(float)Math.pow(10,30);

	//boolean for turning on type dependent DEE i.e. only rotamers of the same type can prune each other
	boolean typeDependent = false;

	//the number of runs
	int numRuns = 1;

	//determines if energy minimization is performed: either traditional-DEE or MinDEE is used
	boolean doMinimize = false;

	//the single and pair interval terms in the MinDEE criterion
	double indIntMinDEE[] = null;
	double pairIntMinDEE[] = null;

	//split flags for all rotamer pairs;
	//only the dead-ending pairs that contain i_r (the rotamer to be pruned) can be
	// 		discarded from the summation; the pairs with i_t and the pairs in the interval terms (if MinDEE/BD/BRDEE)
	// 		must still be a part of the summation.
	boolean splitFlags[][][][][][] = null;

	//determines if split flags are used
	boolean useFlags = false;

	//determines if backbone minimization is performed
	boolean minimizeBB = false;

	
	//the template interval energy (0.0 if fixed backbone)
	float templateInt = 0.0f;

	// 2010: iMinDEE
	boolean doIMinDEE = false;
	float Ival = 0.0f;


        protected int numMutable;
	StrandRotamers strandRot[] = null;
	int strandMut[][] = null;
	int mutRes2Strand[] = null;
	int mutRes2MutIndex[] = null;


        //Used for conformational-splitting DEE
        int numRotForRes[] = null;



        //Pair-related parameters
        //determines if magic bullet or full pairs is used
	boolean magicBullet = false;
        //determines which two residues are in the current pair (only for the distributed DEE)
	boolean resInPair[] = null;
	//determines if distributed DEE is performed
	boolean distrDEE = false;

        //True for triples that are pruned (use tripleFlags[a][b][c] where a>b>c)
        boolean tripleFlags[][][][][][][][][] = null;
        boolean useTriples = false;

        //Use DEEPer
        boolean doPerturbations = false;


        //The constructors for the different classes are very similar, so calling this init function avoids redundancy
        public void init(PairwiseEnergyMatrix arpMatrix, PairwiseEnergyMatrix arpMatrixMax, int numResMutable,
			int strMut[][], float initEw,
			StrandRotamers strandLRot[], PrunedRotamers<Boolean> prunedRotAtRes, boolean doMin, double indInt[],
			double pairInt[], boolean spFlags[][][][][][], boolean useSF, boolean minBB,
                        int mutRes2StrandP[], int mutRes2MutIndexP[], boolean typeDep, boolean iMinDEE, float Ival,
			boolean mb, boolean dDEE, boolean residueMut[], boolean tripFlags[][][][][][][][][], boolean doPerts) {


                doMinimize = doMin;
		typeDependent = typeDep;
		doIMinDEE = iMinDEE;

		pairwiseMinEnergyMatrix = arpMatrix;
		// 2010: No max matrix if doIMinDEE set
		if (doMinimize && !doIMinDEE) //max matrix is different
			pairwiseMaxEnergyMatrix = arpMatrixMax;
		else //no minimization, so the same matrix // 2010: if doIMinDEE is set to true then it is the same as DEE
			pairwiseMaxEnergyMatrix = pairwiseMinEnergyMatrix;

		splitFlags = spFlags;
		eliminatedRotAtPos = prunedRotAtRes;
		//rotIndOffset = rotamerIndexOffset;
		//residueMap = resMap;
		indIntMinDEE = indInt;
		pairIntMinDEE = pairInt;
		//sysLR = systemLRot;
		//rl = rlP;
		useFlags = useSF;
		minimizeBB = minBB;

		strandMut = strMut;
		strandRot = strandLRot;
		mutRes2Strand = mutRes2StrandP;
		mutRes2MutIndex = mutRes2MutIndexP;
		numMutable = numResMutable;

		//numSiteResidues = numResInActiveSite;		// tested with 9
		//numTotalRot = numTotalRotamers;				// ?152?
		//numLigRot = numLigRotamers;					// 0 if no ligand

		numAAtypes = new int[numMutable];

		int ctr=0;
		for(int str=0;str<strandMut.length;str++){ //the number of AAs allowed for each AS residue
			for(int i=0;i<strandMut[str].length;i++){
				numAAtypes[ctr] = strandRot[str].getNumAllowable(strandMut[str][i]);
				ctr++;
			}
		}

		curEw = initEw + Ival;

		numRuns = 1;

		templateInt = 0.0f;
		if (minimizeBB) //backbone minimization, so we need the template interval energy (otherwise, templateInt will be 0.0)
			templateInt = pairwiseMaxEnergyMatrix.getShellShellE()-pairwiseMinEnergyMatrix.getShellShellE();


                //Some stuff for pairs
                resInPair = residueMut;
		distrDEE = dDEE;
		magicBullet = mb;

                if(tripFlags != null){//Use triples...if they are not going to be used, tripFlags will be null
                    tripleFlags = tripFlags;
                    useTriples = true;
                }

                doPerturbations = doPerts;
        }





	//return the split flags for all rotamer pairs
	public boolean[][][][][][] getSplitFlags(){
		return splitFlags;
	}


        //find how many rotamers are allowed for the current AA type at the given residue
        //(strand-based numbering)
        //str=strand number
        //gly and ala are counted as having 1 rotamer each
        protected int getNumRot(int str, int strResNum, int curAA){

            if(doPerturbations)
                return ((StrandRCs)strandRot[str]).getNumRCs(strResNum, curAA);
            else{
                int ans = strandRot[str].rl.getNumRotForAAtype(curAA);
                if (ans==0)	//ala or gly
                    return 1;
                else
                    return ans;
            }
        }


        //Compute the number of rotamers for each residue position (assign to numRotForRes[])
	protected void compNumRotForRes(){

		//boolean ligPresent = (numLigRot==0); //ligand present
		int treeLevels = numMutable;
		/*if (ligPresent)
			treeLevels++;*/

		numRotForRes = new int[treeLevels];

		int curNumRot = 0;
		for (int curLevel=0; curLevel<treeLevels; curLevel++){
			/*if ((ligPresent)&&(curLevel==(treeLevels-1))){ //the ligand level
				curNumRot = numLigRot;
			}
			else {*/ //AS residue
				curNumRot = 0;
				int str=mutRes2Strand[curLevel];
				int strResNum=strandMut[str][mutRes2MutIndex[curLevel]];

				for (int i=0; i<strandRot[str].getNumAllowable(strResNum); i++){ //add the rot for all allowable AA at this residue
					int newRot = getNumRot( str, strResNum, strandRot[str].getIndexOfNthAllowable(strResNum,i) );

					curNumRot += newRot;
				}
			//}
			numRotForRes[curLevel] = curNumRot;
		}
	}


        boolean isPrunedTriple(int curPos1, int curAA1, int curRot1,
                int curPos2, int curAA2, int curRot2,
                int curPos3, int curAA3, int curRot3 ){
            //Checks if a given triple is pruned in tripleFlags
            //curPos1, curPos2, and curPos3 (positions among mutable residues of
            //residues 1, 2, and 3) are assumed to be all different (otherwise it's not a prunable triple)
            //If useTriples==true, calling this with pruned rotamers might yield a null-pointer exception
            //(because of the space-saving structure of tripleFlags)

            if(!useTriples)
                return false;

            if( curPos1 > curPos2 ){
                if( curPos2 > curPos3 )
                    return tripleFlags[curPos1][curAA1][curRot1][curPos2][curAA2][curRot2][curPos3][curAA3][curRot3];
                else if ( curPos1 > curPos3 )
                    return tripleFlags[curPos1][curAA1][curRot1][curPos3][curAA3][curRot3][curPos2][curAA2][curRot2];
                else//index3 > index1 > index2
                    return tripleFlags[curPos3][curAA3][curRot3][curPos1][curAA1][curRot1][curPos2][curAA2][curRot2];
            }
            else{
                if( curPos1 > curPos3 )
                    return tripleFlags[curPos2][curAA2][curRot2][curPos1][curAA1][curRot1][curPos3][curAA3][curRot3];
                else if ( curPos2 > curPos3 )
                    return tripleFlags[curPos2][curAA2][curRot2][curPos3][curAA3][curRot3][curPos1][curAA1][curRot1];
                else
                    return tripleFlags[curPos3][curAA3][curRot3][curPos2][curAA2][curRot2][curPos1][curAA1][curRot1];
            }
        }


}
