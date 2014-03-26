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
//	PairwiseEnergyMatrix.java
//
//	Version:           2.1 beta
//
//
//	  authors:
// 	  initials    name                 organization                email
//	 ---------   -----------------    ------------------------    ----------------------------
//	  MAH           Mark A. Hallen	  Duke University               mah43@duke.edu
///////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////
//	Adapted from: PEMHandler.java
//
//	Version:           2.0
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
 * Manages operations on the pairwise energy matrices.
 * 
*/
public class PairwiseEnergyMatrix {

        float eMatrix[][][][][][] = null;

	/* 
	 * Initialize the pairwise energy matrices: each matrix has 6 dimensions;
	 * The first three dimensions correspond to the residue position, amino acid type, and rotamer identity
	 * 		for the first rotamer in the rotamer pair; the last three dimensions define the second rotamer in
	 * 		the rotamer pair;
	 * The first and fourth dimensions are of size (numInAS+1) if the input system does not have a ligand, or (numInAS+2) otherwise;
	 * 		The (numInAS+1)^th row in the first/fourth dimensions corresponds to the ligand entries (if present);
	 * 		The last row in the first dimension corresponds to the template energy (the last row in the fourth dimension is never used);
	 * The intra-rotamer and rot-shell energies for a given rotamer identity 'r' for amino acid type 'a' at residue position 'p'
	 * 		are stored in the two entries in [p][a][r][p][0][] (since there are no rotamer pairwise energies for rotamers at the 
	 * 		same residue position);
	 * NOTE: the difference between ligPresent and useLig is the following:
	 * 		ligPresent determines if a ligand is present in the input structure;
	 * 		useLig determines if the ligand will be used in the current computation (e.g., useLig will be false for SHL-AS)
	 */
        public PairwiseEnergyMatrix(){
            
        }


	public PairwiseEnergyMatrix(int numMutable, int resMut[], int[][] strandMut,
			RotamerSearch rs, boolean shellRun, boolean intraRun, boolean initAll) {
		
		int numPos = numMutable+1;
		/*if (ligPresent)
			numPos++;*/
		eMatrix = new float[numPos][][][][][];
		int p1 = 0;
		for (int str1=0; str1<strandMut.length; str1++){
                    for (int i=0; i<strandMut[str1].length; i++){
			if (resMut[p1] == 1) {
                            eMatrix[p1] = new float[rs.strandRot[str1].rl.getNumAAallowed()][][][][];
                            for (int a1=0; a1<rs.strandRot[str1].getNumAllowable(strandMut[str1][i]); a1++){
                                int curAAind1 = rs.strandRot[str1].getIndexOfNthAllowable(strandMut[str1][i],a1);
                                int numRot1 = rs.getNumRot(str1, strandMut[str1][i], curAAind1);
                                eMatrix[p1][curAAind1] = new float[numRot1][][][];
                                for (int r1=0; r1<numRot1; r1++){
                                    eMatrix[p1][curAAind1][r1] = new float[numPos][][];
                                    int p2=0;
                                    for (int str2=0;str2<strandMut.length;str2++){
                                            for (int j=0; j<strandMut[str2].length; j++){
                                                    if ((initAll)||((!intraRun)&&(p2!=p1))){
                                                            eMatrix[p1][curAAind1][r1][p2] = new float[rs.strandRot[str2].rl.getNumAAallowed()][];
                                                            for (int a2=0; a2<rs.strandRot[str2].getNumAllowable(strandMut[str2][j]); a2++){
                                                                int curAAind2 = rs.strandRot[str2].getIndexOfNthAllowable(strandMut[str2][j],a2);
                                                                int numRot2 = rs.getNumRot( str2, strandMut[str2][j], curAAind2 );
                                                                eMatrix[p1][curAAind1][r1][p2][curAAind2] = new float[numRot2];
                                                                for (int r2=0; r2<numRot2; r2++){
                                                                        eMatrix[p1][curAAind1][r1][p2][curAAind2][r2] = 0.0f;
                                                                }
                                                            }
                                                    }
                                                    if (p2==p1) {
                                                            //Be able to store intra rotamer energy and energy with each strand template
                                                            eMatrix[p1][curAAind1][r1][p2] = new float[1][2];
                                                    }
                                                    p2++;
                                            }
                                    }
                                }
                            }
			}
			p1++;		
                    }
		}

		if (shellRun){
			eMatrix[numPos-1] = new float[1][1][1][1][1];
		}
		
	}
	
	
	//Returns a new independent six-dimensional matrix that is a copy of fromMatrix[][][][][][]
	public PairwiseEnergyMatrix copy(){
		
                PairwiseEnergyMatrix newM = new PairwiseEnergyMatrix();

		if (eMatrix==null)
			return null;
		
		newM.eMatrix = new float[eMatrix.length][][][][][];
		for (int p1=0; p1<newM.eMatrix.length; p1++){
			if (eMatrix[p1]!=null){
				newM.eMatrix[p1] = new float[eMatrix[p1].length][][][][];
				for (int a1=0; a1<newM.eMatrix[p1].length; a1++){
					if (eMatrix[p1][a1]!=null){
						newM.eMatrix[p1][a1] = new float[eMatrix[p1][a1].length][][][];
						for (int r1=0; r1<newM.eMatrix[p1][a1].length; r1++){
							if (eMatrix[p1][a1][r1]!=null){
								newM.eMatrix[p1][a1][r1] = new float[eMatrix[p1][a1][r1].length][][];
								for (int p2=0; p2<newM.eMatrix[p1][a1][r1].length; p2++){
									if (eMatrix[p1][a1][r1][p2]!=null){
										newM.eMatrix[p1][a1][r1][p2] = new float[eMatrix[p1][a1][r1][p2].length][];
										for (int a2=0; a2<newM.eMatrix[p1][a1][r1][p2].length; a2++){
											if (eMatrix[p1][a1][r1][p2][a2]!=null){
												newM.eMatrix[p1][a1][r1][p2][a2] = new float[eMatrix[p1][a1][r1][p2][a2].length];
												System.arraycopy(eMatrix[p1][a1][r1][p2][a2], 0, newM.eMatrix[p1][a1][r1][p2][a2], 0, eMatrix[p1][a1][r1][p2][a2].length);
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
		
		return newM;
	}




        //Accessing and changing elements of energy matrices






        //The following methods access elements based on residue number among the flexible residues,
        //AA index, and rotamer (or RC) number


        public void setPairwiseE(int res1, int AA1, int rot1, int res2, int AA2, int rot2, float val){
            eMatrix[res1][AA1][rot1][res2][AA2][rot2] = val;
            eMatrix[res2][AA2][rot2][res1][AA1][rot1] = val;
        }

        public float getPairwiseE(int res1, int AA1, int rot1, int res2, int AA2, int rot2){
            return eMatrix[res1][AA1][rot1][res2][AA2][rot2];
        }

        //Add the given value to the specified pairwise energy
        public void addToPairwiseE(int res1, int AA1, int rot1, int res2, int AA2, int rot2, float val){
            eMatrix[res1][AA1][rot1][res2][AA2][rot2] += val;
            eMatrix[res2][AA2][rot2][res1][AA1][rot1] += val;
        }


        //Shell energy
        public float getShellShellE(){
		float shlshlE = 0.0f;
		for(int i=0; i<eMatrix[eMatrix.length-1][0][0][0].length;i++){
			for(int j=i; j<eMatrix[eMatrix.length-1][0][0][0][0].length;j++){
				shlshlE += eMatrix[eMatrix.length-1][0][0][0][i][j];
			}
		}
		return shlshlE;
	}


        public void setShellShellE(float val){
                eMatrix[eMatrix.length-1][0][0][0][0][0] = val;
	}


        //Shell-rotamer energy
	public float getShellRotE(int pos, int AANum, int rot){
		float shlRotE = 0.0f;
		//Skip the first energy which is the intra-rot energy still
		for(int i=1; i<eMatrix[pos][AANum][rot][pos][0].length;i++){
			shlRotE += eMatrix[pos][AANum][rot][pos][0][i];
		}
		return shlRotE;
	}

	public void setShellRotE(int pos, int AANum, int rot, float val){
		eMatrix[pos][AANum][rot][pos][0][1] = val;
	}

        public void addToShellRotE(int pos, int AANum, int rot, float val){
		eMatrix[pos][AANum][rot][pos][0][1] += val;
	}



        //Intra + shell-interaction energy for a rotamer or RC
        public float getIntraAndShellE(int pos, int AANum, int rot){
            float E = 0.0f;
            for(int i=0; i<eMatrix[pos][AANum][rot][pos][0].length;i++){
                    E += eMatrix[pos][AANum][rot][pos][0][i];
            }
            return E;
        }



        //Intra energies
        public void setIntraE(int res1, int AA1, int rot1, float val){
            eMatrix[res1][AA1][rot1][res1][0][0] = val;
        }

        public float getIntraE(int res1, int AA1, int rot1){
            return eMatrix[res1][AA1][rot1][res1][0][0];
        }

        //Add the given value to the specified pairwise energy
        public void addToIntraE(int res1, int AA1, int rot1, float val){
            eMatrix[res1][AA1][rot1][res1][0][0] += val;
        }




        //Ways to reorganize the energies


        
        //Generate reduced energy matrix for A*
        //eliminatedRotAtRes is reduced to contain only the entries for the current sequence
        //(if singleSeq, for use in K*) or for sequences being considered if !singleSeq
        //the energy matrix is reduced to a 2D matrix containing only unpruned rotamers/RCs
        //Called from RotamerSearch
        //eliminatedRotAtPosRed is filled in
        //if rs.useFlagsAStar==true then the split and triple flags are reduced too
        public ReducedEnergyMatrix reduceMatrix(boolean eliminatedRotAtPosRed[],
			int numRotForRes[], int numRotForResNonPruned[], int treeLevels,
			int numTotalRotRedNonPruned, int numMutable, int strandMut[][],
                        int numTotalRotRed, RotamerSearch rs, boolean singleSeq,
                        boolean splitFlagsRed[][], boolean tripleFlagsRed[][][]){


                int indicesEMatrixPos[] = new int[numTotalRotRedNonPruned]; //original (in the non-reduced matrices) indices of non-pruned rot to be included
		int indicesEMatrixAA[] = new int[numTotalRotRedNonPruned];
		int indicesEMatrixRot[] = new int[numTotalRotRedNonPruned];
		float arpMatrixRed[][] = new float[numTotalRotRedNonPruned+eMatrix[eMatrix.length-1][0][0][0][0].length][numTotalRotRedNonPruned+1];//include the intra-energies in the last column
		

		int curIndexRed = 0;//index into the reduced matrices
		int pruningIndex = 0;//index into the reduced MinDEE matrix

		for (int curLevel=0; curLevel<treeLevels; curLevel++){
			int str = rs.mutRes2Strand[curLevel];
			int strResNum = strandMut[str][rs.mutRes2StrandMutIndex[curLevel]];
                        int molResNum = rs.m.strand[str].residue[strResNum].moleculeResidueNumber;

                        int numAA=1;
                        if(!singleSeq)//not a K* run
                            numAA = rs.strandRot[str].getNumAllowable(strResNum);

                        for (int curAA=0; curAA<numAA; curAA++){ //for all allowed AA's


                                int index, newRot;
                                
                                if(singleSeq){
                                    index = rs.curAANum[molResNum];
                                    newRot = numRotForRes[curLevel];
                                }
                                else{
                                    index = rs.strandRot[str].getIndexOfNthAllowable(strResNum,curAA);
                                    newRot = rs.getNumRot(str, strResNum, index);
                                }


                                for (int curRot=0; curRot<newRot; curRot++){ //for all rotamers for the given AA

                                        if (!rs.eliminatedRotAtRes.get(curLevel,index,curRot)){ //not pruned, so add its index
                                                indicesEMatrixPos[curIndexRed] = curLevel;
                                                indicesEMatrixAA[curIndexRed] = index;
                                                indicesEMatrixRot[curIndexRed] = curRot;
                                                curIndexRed++;
                                        }
                                        eliminatedRotAtPosRed[pruningIndex] = rs.eliminatedRotAtRes.get(curLevel,index,curRot);
                                        //logPS.println(pruningIndex+" "+curIndex+" "+eliminatedRotAtRes[curIndex]);logPS.flush();
                                        pruningIndex++;
                                }
                        }
		}


                
                if(!singleSeq){
                    
                    System.out.println("pruneIndex "+pruningIndex);
                    for (int i=0;i<curIndexRed;i++)System.out.print("("+indicesEMatrixPos[i]+" "+indicesEMatrixAA[i]+" "+indicesEMatrixRot[i]+") ");System.out.println("curIndexRed "+curIndexRed);
                }



		//Reduce the min energy matrix
		for (int curRot1=0; curRot1<numTotalRotRedNonPruned; curRot1++){
			int p1 = indicesEMatrixPos[curRot1];
			int a1 = indicesEMatrixAA[curRot1];
			int r1 = indicesEMatrixRot[curRot1];
			for (int curRot2=0; curRot2<numTotalRotRedNonPruned; curRot2++){
				if (curRot1!=curRot2){
					int p2 = indicesEMatrixPos[curRot2];
					if (p1!=p2) {//not the same residue position
						int a2 = indicesEMatrixAA[curRot2];
						int r2 = indicesEMatrixRot[curRot2];
						arpMatrixRed[curRot1][curRot2] = getPairwiseE(p1,a1,r1,p2,a2,r2);//pairwise
					}
				}
			}
			arpMatrixRed[curRot1][numTotalRotRedNonPruned] = getIntraE(p1,a1,r1);//store intra-energies in the last column
			for(int i=0; i<eMatrix[p1][a1][r1][p1][0].length-1;i++)
				arpMatrixRed[numTotalRotRedNonPruned+i][curRot1] = eMatrix[p1][a1][r1][p1][0][1+i];//store shell-rotamer E in the last rows
		}



                if(rs.useFlagsAStar){//Reduce the flags for pruned tuples
                    
                    for(int a=0;a<numTotalRotRedNonPruned;a++){
                        for(int b=0;b<numTotalRotRedNonPruned;b++){
                            if(indicesEMatrixPos[a]!=indicesEMatrixPos[b])//Pairs of rotamers at the same residue aren't meaningful
                                splitFlagsRed[a][b] = rs.splitFlags[indicesEMatrixPos[a]][indicesEMatrixAA[a]][indicesEMatrixRot[a]]
                                        [indicesEMatrixPos[b]][indicesEMatrixAA[b]][indicesEMatrixRot[b]];
                        }

                        if(rs.useTriples){

                            tripleFlagsRed[a] = new boolean[a][];//There are at most a-1 values of b such that a>b>=0

                            for(int b=0; indicesEMatrixPos[b] < indicesEMatrixPos[a]; b++){

                                tripleFlagsRed[a][b] = new boolean[b];

                                //tripleFlags doesn't store triples repeating a residue position; this criterion accounts for this
                                for(int c=0; indicesEMatrixPos[c] < indicesEMatrixPos[b]; c++)
                                    tripleFlagsRed[a][b][c] = rs.tripleFlags[indicesEMatrixPos[a]][indicesEMatrixAA[a]][indicesEMatrixRot[a]]
                                            [indicesEMatrixPos[b]][indicesEMatrixAA[b]][indicesEMatrixRot[b]]
                                            [indicesEMatrixPos[c]][indicesEMatrixAA[c]][indicesEMatrixRot[c]];
                            }
                        }
                    }
                }



                return new ReducedEnergyMatrix(arpMatrixRed);

        }
        



        
        

	
	//Called by slave nodes to generate cObj.compEE[] entries to return to the main node
        //Called as a method of the minimum-eneryg matrix with the maximum-eneryg matrix as an argument
	//The two matrices should have the same structure (i.e., a computed entry in one matrix should also be computed in the other)
	public SamplingEEntries [] generateCompEE(PairwiseEnergyMatrix maxM){
		
		if ( maxM == null ) {
			System.out.println("ERROR: cannot generate compEE[] entries from a null PEM matrix.");
			System.exit(1);
		}
		
		else {
			SamplingEEntries compEE[] = new SamplingEEntries[100];
			int curEntry = 0;
			for (int p1=0; p1<eMatrix.length; p1++){
				if (eMatrix[p1]!=null){
					for (int a1=0; a1<eMatrix[p1].length; a1++){
						if (eMatrix[p1][a1]!=null){
							for (int r1=0; r1<eMatrix[p1][a1].length; r1++){
								if (eMatrix[p1][a1][r1]!=null){
									for (int p2=0; p2<eMatrix[p1][a1][r1].length; p2++){
										if (eMatrix[p1][a1][r1][p2]!=null){
											for (int a2=0; a2<eMatrix[p1][a1][r1][p2].length; a2++){
												if (eMatrix[p1][a1][r1][p2][a2]!=null){
													for (int r2=0; r2<eMatrix[p1][a1][r1][p2][a2].length; r2++){
														if ( (eMatrix[p1][a1][r1][p2][a2][r2]!=0.0f) || (maxM.eMatrix[p1][a1][r1][p2][a2][r2]!=0.0f) ) {
														
															compEE[curEntry] = new SamplingEEntries();
															compEE[curEntry].i1 = p1;
															compEE[curEntry].i2 = a1;
															compEE[curEntry].i3 = r1;
															compEE[curEntry].i4 = p2;
															compEE[curEntry].i5 = a2;
															compEE[curEntry].i6 = r2;
															compEE[curEntry].minE = eMatrix[p1][a1][r1][p2][a2][r2];
															compEE[curEntry].maxE = maxM.eMatrix[p1][a1][r1][p2][a2][r2];
															
															curEntry++;
															
															if (curEntry>=compEE.length){
																SamplingEEntries tmp[] = new SamplingEEntries[compEE.length*2];
																System.arraycopy(compEE, 0, tmp, 0, compEE.length);
																compEE = tmp;
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
			
			SamplingEEntries tmp[] = new SamplingEEntries[curEntry];
			System.arraycopy(compEE, 0, tmp, 0, curEntry);
			compEE = tmp;
			
			return compEE;
		}
		
		return null;
	}





        //Initialing other matrices to contain information on rotamer pairs
        //Using the same sets of rotamer pairs at this pairwise energy matrix

        //Pairwise matrix of booleans (used for split flags)
        public boolean[][][][][][] initializePairwiseBooleanMatrix(){

                boolean toMatrix[][][][][][] = new boolean[eMatrix.length][][][][][];
		//KER: we only want the pairwise interactions and not the shell-shell in the last row
		for (int p1=0; p1<eMatrix.length-1; p1++){
			if (eMatrix[p1]!=null){
				toMatrix[p1] = new boolean[eMatrix[p1].length][][][][];
				for (int a1=0; a1<eMatrix[p1].length; a1++){
					if (eMatrix[p1][a1]!=null){
						toMatrix[p1][a1] = new boolean[eMatrix[p1][a1].length][][][];
						for (int r1=0; r1<eMatrix[p1][a1].length; r1++){
							if (eMatrix[p1][a1][r1]!=null){
								toMatrix[p1][a1][r1] = new boolean[eMatrix[p1][a1][r1].length][][];
								for (int p2=0; p2<eMatrix[p1][a1][r1].length; p2++){
									if (eMatrix[p1][a1][r1][p2]!=null){
										toMatrix[p1][a1][r1][p2] = new boolean[eMatrix[p1][a1][r1][p2].length][];
										for (int a2=0; a2<eMatrix[p1][a1][r1][p2].length; a2++){
											if (eMatrix[p1][a1][r1][p2][a2]!=null){
												toMatrix[p1][a1][r1][p2][a2] = new boolean[eMatrix[p1][a1][r1][p2][a2].length];
												for (int r2=0; r2<eMatrix[p1][a1][r1][p2][a2].length; r2++){
													toMatrix[p1][a1][r1][p2][a2][r2] = false;
												}
												//System.arraycopy(eMatrix[p1][a1][r1][p2][a2], 0, toMatrix[p1][a1][r1][p2][a2], 0, eMatrix[p1][a1][r1][p2][a2].length);
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

                return toMatrix;
        }




        //Pairwise matrix of doubles (used for Ec in BoundFlags)
        //Initialize with the given value
        public double[][][][][][] initializePairwiseDoubleMatrix(double val){

                double toMatrix[][][][][][] = new double[eMatrix.length][][][][][];
		//KER: we only want the pairwise interactions and not the shell-shell in the last row
		for (int p1=0; p1<toMatrix.length-1; p1++){
			if (eMatrix[p1]!=null){
				toMatrix[p1] = new double[eMatrix[p1].length][][][][];
				for (int a1=0; a1<toMatrix[p1].length; a1++){
					if (eMatrix[p1][a1]!=null){
						toMatrix[p1][a1] = new double[eMatrix[p1][a1].length][][][];
						for (int r1=0; r1<toMatrix[p1][a1].length; r1++){
							if (eMatrix[p1][a1][r1]!=null){
								toMatrix[p1][a1][r1] = new double[eMatrix[p1][a1][r1].length][][];
								for (int p2=0; p2<toMatrix[p1][a1][r1].length; p2++){
									if (eMatrix[p1][a1][r1][p2]!=null){
										toMatrix[p1][a1][r1][p2] = new double[eMatrix[p1][a1][r1][p2].length][];
										for (int a2=0; a2<toMatrix[p1][a1][r1][p2].length; a2++){
											if (eMatrix[p1][a1][r1][p2][a2]!=null){
												toMatrix[p1][a1][r1][p2][a2] = new double[eMatrix[p1][a1][r1][p2][a2].length];
												for (int r2=0; r2<toMatrix[p1][a1][r1][p2][a2].length; r2++){
													toMatrix[p1][a1][r1][p2][a2][r2] = val;
												}
												//System.arraycopy(eMatrix[p1][a1][r1][p2][a2], 0, toMatrix[p1][a1][r1][p2][a2], 0, eMatrix[p1][a1][r1][p2][a2].length);
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

		return toMatrix;


        }


}
