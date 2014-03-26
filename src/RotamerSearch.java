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
// RotamerSearch.java
//
//  Version:           2.1 beta
//
//
// authors:
//    initials    name                 organization                email
//   ---------   -----------------    ------------------------    ----------------------------
//     RHL        Ryan Lilien          Dartmouth College           ryan.lilien@dartmouth.edu
//	   ISG		  Ivelin Georgiev	   Duke University			   ivelin.georgiev@duke.edu
//	  KER        Kyle E. Roberts       Duke University         ker17@duke.edu
//    PGC        Pablo Gainza C.       Duke University         pablo.gainza@duke.edu
//     MAH        Mark A. Hallen	Duke University         mah43@duke.edu
///////////////////////////////////////////////////////////////////////////////////////////////

/** 
 * Written by Ryan Lilien (2001-2004) and Ivelin Georgiev (2004-2009)
 * 
 * This class provides a variety of tools and search algorithms for
 *  doing rotamer-based searching over molecular conformations
 *
 * The system consists of one molecule containing: a protein strand, 
 *  a ligand strand (optional), and a cofactor strand (optional).
 * The protein strand does not have to contain sequential residues
 *  but it must be made of standard amino acids
 * The ligand strand can only be one 'thing'
 *  -if this 'thing' is an AA then the Penultimate Rotamer library is used
 *  -if this 'thing' is not an AA then a generic rotamer library is used
 *
 */


import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintStream;
import java.io.BufferedOutputStream;
import java.io.FileOutputStream;
import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.io.FileInputStream;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.PrintWriter;
import java.io.Serializable;
import java.util.*;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.math.*;

/**
 * 
 * This class provides a variety of tools and search algorithms for
 *  doing rotamer-based searching over molecular conformations.
 * Contains functions for computing the pairwise energy matrices and
 *  for performing DEE/A* and K* (with different types of minimization)
 *  mutation searches.
 * The functions in this class are typically called after the necessary
 *  setup in KSParser.java is performed.
 *
 */
public class RotamerSearch implements Serializable
{

	// If the debug flag is set to true then additional debug statements are
	//  printed to standard out.
	public static final boolean debug = true;
	final double constRT = 1.9891/1000.0 * 298.15;   // in kCal / kelvin-mole (the T value here should be consistent with T in EEF1)
	
	Molecule m;		// the molecule
	Amber96ext a96ff;	// the forcefield and energy function to use for energy evaluation etc...
	Amber96ext templateA96ff; //KER: forcefield that keeps the energy computation for the template 
	Amber96ext res1A96ff; //KER: forcefields to store terms for single amino acid energy 
    Amber96ext res2A96ff;
    Amber96ext templateOffA96ff; //KER: templateOff is to compute the energy since we don't want to include backbone atoms in template
    							//for the matrix, only for the minimization
	
	SimpleMinimizer simpMin;	// the simple energy minimizer (side-chains)
	BBMinimizer bbMin = null;	//the backbone minimizer
	BackrubMinimizer brMin = null;	//the backrub minimizer
	//KER: eliminatedRotAtRes is now a 3D matrix for easier indexing (I don't think this will save space like the 6D array)
	PrunedRotamers<Boolean> eliminatedRotAtRes = null;	// rotamers pruned by MinDEE
	//KER: splitFlags will now be a 6D array to save space
	boolean splitFlags[][][][][][] = null;

        boolean tripleFlags[][][][][][][][][] = null;//Pruned triples
        boolean useTriples = false;             //Indicates whether triples should be used
        boolean useFlagsAStar = false;          //Should split and potentially triple flags be used in A*?

	PrunedRotamers<Boolean> prunedIsSteric = null;	// pruned rotamers due to unallowed sterics
	double indIntMinDEE[] = null;	//the single-residue interval term in the MinDEE criterion
	double pairIntMinDEE[] = null;	//the pairwise interval term in the MinDEE criterion
	boolean repeatSearch = false;	// determines if the search must be repeated to achieve the desired accuracy
	int curConf[] = null;		// the current conformation returned by A*
	boolean allPruned = false;	// determines if all the rotamers for a given residue have been pruned by MinDEE;
									// sends this information to the master node in the mutation search
	
	final float stericE = (float)Math.pow(10,38);	// the energy stored for an unallowed steric
	private double Ec_const = stericE;	// the minimum lower energy bound for a pruned conformation
	private double boundForPartition = stericE; //a lower bound on the conformational energy for a given partition (used by DACS)
	
	final int samplesBB = 1000;//number samples for intra-rotamer energy matrix computation with backbone flexibility
	
	boolean distDepDielect = true; //distance-dependent dielectric
	double dielectConst = 1.0; //the dielectric constant
	
	boolean doDihedE = false; 		// if true dihedral energies are computed and used
									//  during energy minimization. Note that if
									//  energy minimization is NOT used then dihedral
									//  energies are not explicitly computed. In
									//  reality the total energy values are the same
									//  because although we're using AMBER dihedral
									//  energy terms we assume that each dihedral of
									//  each rotamer is at the bottom of an energy
									//  well. Thus without minimization the total
									//  dihedral energy is zero.
	boolean doSolvationE = false; //determines if solvation energies should be computed
	
	double solvScale = 1.0; //the solvation energies scaling factor
	
	PriorityQueue<ConfPair> topConfs = null; //Stores the top conformations of the run
	//Right now this is used to print out the conformations at the end
	
	//RotamerLibrary rl = null; //the standard rotamer library for the protein (the AA library)
	//RotamerLibrary[] grl = null; //the rotamer library for the ligand (could be the AA or non-AA library)
	//StrandRotamers sysLR = null;// the rotamers object for the system strand
	//StrandRotamers ligROT = null;	// the rotamers object for the ligand strand
	//int sysStrNum = -1;	// the strand number of the system
	//int ligStrNum = -1;	// the strand number of the ligand
	float overlapThresh = -10000.0f;	// hard overlap threshold used for checking sterics (this should be used when the atom positions will not be allowed to change after the steric check)
	float softOverlapThresh = -10000.0f; //soft overlap threshold used for checking sterics (this should be used when the atom positions may be allowed to change after the steric check)
	int curAANum[] = null;	// for each residue in the system strand, the
			// index of the current amino acid type; if the residue is not
			// rotamerizable (ie. it's not flexible) then the curAANum entry
			// should be -1
	//int curLigAANum = -1;	// the index of the current ligand type
	boolean computeEVEnergy = false;
			// do we compute EV energies during a conformation search
	boolean doMinimization = false;
			// do we some EV minimization steps during a conformation search
	boolean hElect = true;
		// should hydrogens be used in electrostatic energy calculations
	boolean hVDW = true;
		// should hydrogens be used in vdw energy calculations
	boolean hSteric = false; //should hydrogens be used in steric checks
	double vdwMultiplier = 1.0f;
		// vdw multiplier used in energy evaluation
	boolean addHydrogens = true;
		// during a mutation, should hydrogens be included when
		//  changing residue type
	boolean connectResidues = true;
		// during a mutation, should a new residue be bonded to
		//  the prior and subsequent resiudes if the numbering
		//  is sequential
	int curConfNum = 0;
	int numMinSteps = 35; //140
		// number of minimization steps to perform by simpmin

        boolean doPerturbations;//Indicates DEEPer
        boolean minimizePerturbations;//When doing DEEPer, indicates perturbation-minimizing mode (i.e. DEE/A* step)
        String pertFile;//File with the perturbation information

        float switchedTemplateE[] = null;//The changes in template energies caused by full structure switch perturbations


	private PairwiseEnergyMatrix arpMatrix = null;
		// all rotamer pairs lower min energy bound matrix, created with
		//  simplePairwiseMutationAllRotamerSearch and loaded with
		//  loadPairwiseEnergyMatrices()
	private PairwiseEnergyMatrix arpMatrixMax = null;
		// all rotamer pairs lower max energy bound matrix, created with
		//  simplePairwiseMutationAllRotamerSearchMax and loaded with
		//  loadPairwiseEnergyMatricesMax()

	int ASAANums[] = null;
		// integer array containing the index for each AS residues
		//  that can be used with rotamerIndexOffset and for the arpMatrix
	int curStrRotNum[] = null;
		// integer array containing the currently assumed rotamer for
		//  each amino acid in the active site
		// is allocated during a rotamer search.
		// note that it is _not_ the same size as curAANum
	//int curLigRotNum = 0;
		// the current rotamer number of the ligand
	float bestEMin = 9999999.0f; //this should ONLY be accessed/modified using the synchronized methods below
	float bestEUnMin = 9999999.0f;
		// the best minimized and unminimized energy found thus far
	BigInteger numConfsTotal = new BigInteger("0");
		// the number of total conformations for the current configuration
		// this is created and computed in computeTotalNumConfs()
	BigInteger numConfsLeft = new BigInteger("0");
		// the number of remaining conformations for the current configuration
		//  updated as the search progresses
	BigInteger numConfsBelowLevel[] = null;
		// the number of conformations below the specified level, level 0
		//  refers to the ligand level, if there's no ligand then level 0
		//  and level 1 have the same value
		// this is created and computed in computeTotalNumConfs()
	BigInteger numConfsAboveLevel[] = null; //at level i, num confs from level i+1 to the last level
	BigInteger numConfsPrunedByE = new BigInteger("0");
		// the number of conformations not minimized because their 'best'
		//  energy (as computed from the arpMatrix) was too unfavorable
		//  based on the accuracy threshold below
	BigInteger numConfsPrunedByS = new BigInteger("0");
		// number of conformations pruned due to a steric clash
	BigInteger numConfsPrunedByMinDEE = new BigInteger("0");	//the number of confs pruned by MinDEE
	BigInteger numConfsEvaluated = new BigInteger("0");
		// number of conformations that got all the way down to the energy
		//  evaluation
		// Note that numConfsPrunedByE + numConfsPrunedByS + numConfsEvaluated
		//  should equal the total number of conformations
	float KSepsilon = 0.03f;
		// the accuracy for computing energies for K*
		// a value of 0.03 means the energies computed
		//  will allow for a calculation of K*_approx
		//  that's within 3% of the true K*
	BigDecimal partial_q = new BigDecimal(0.0);
		// the partially computed partition function (updated as we go)
	BigDecimal partial_p = new BigDecimal(0.0);
		// the bound on the partition function of the pruned conformations
	BigDecimal initial_q = new BigDecimal(0.0);
		// used in mutation search as an initial partial_q if we're
		//  bootstrapping the search
	StrandRotamers[] strandRot = null;
	int numberOfStrands = -1;
	int mutRes2Strand[] = null;
	int mutRes2StrandMutIndex[] = null;
	int numberMutable = 0;

	boolean isTemplateOn = false;
	
	// Note that the mutation search functions in this class are relatively
	//  messy as a result of changing the algorithms multiple times. They
	//  could be rewritten to be much tighter and more elegant.

	// the constructor if you also have a ligand
	RotamerSearch(Molecule theMolec, int numMut,int strandsPresent, boolean hE, 
			boolean hV, boolean hS, boolean addH, boolean conRes, float eps, 
			float stericThresh, float softStericThresh,	boolean ddDielect, 
			double dielectC, boolean doDihedral, boolean doSolv,double solvScFactor, 
			double vdwMult, RotamerLibrary[] rotamerLibraries,
                        boolean doPerts, String pFile, boolean minPerts, boolean useTrips, boolean flagsAStar) {
		
		hElect = hE;
		hVDW = hV;
		hSteric = hS;
		addHydrogens = addH;
		connectResidues = conRes;
		KSepsilon = eps;
		overlapThresh = stericThresh;
		softOverlapThresh = softStericThresh;
		vdwMultiplier = vdwMult;
		distDepDielect = ddDielect;
		dielectConst = dielectC;
		doDihedE = doDihedral;
		doSolvationE = doSolv;
		solvScale = solvScFactor;
		numberMutable = numMut;
		setBestE(stericE);

                doPerturbations=doPerts;
                pertFile = pFile;
                minimizePerturbations = minPerts;
                useTriples = useTrips;
                useFlagsAStar = flagsAStar;

		m=theMolec;
		a96ff = new Amber96ext(m, distDepDielect, dielectConst, doSolvationE, solvScale, vdwMultiplier);
		a96ff.calculateTypesWithTemplates(); //KER: Needed so that the N and C terminus flags are set correctly

                numberOfStrands = strandsPresent;
		strandRot = new StrandRotamers[numberOfStrands];
		
		curAANum = new int[m.numberOfResidues]; //

                if(doPerturbations){
                    simpMin = new PMinimizer(minimizePerturbations);
                    for(int i=0; i<numberOfStrands;i++){
			strandRot[i] = new StrandRCs(rotamerLibraries[i],m.strand[i]);
                    }
                }
                else{
                    simpMin = new SimpleMinimizer();
                    bbMin = new BBMinimizer();
                    brMin = new BackrubMinimizer();
                    for(int i=0; i<numberOfStrands;i++){
			strandRot[i] = new StrandRotamers(rotamerLibraries[i],m.strand[i]);
                    }
                }
		
		
		//sysLR = new StrandRotamers(rl,m.strand[sysStrNum]);
		/*if (ligStrNum>=0) { //there is a ligand
			m.strand[ligStrNum].residue[0].flexible = true;
			ligROT = new StrandRotamers(grl,m.strand[ligStrNum]);		
			ligROT.setAllowable(0,m.strand[ligStrNum].residue[0].name);
		}*/
	}


        //Set up perturbations and residue conformations (RCs) for DEEPer
        //This should be called after the allowable AAs are set
        public void setupRCs(boolean addWTRot){

                PertFileHandler.readPertFile(pertFile, m, strandRot);

                for (int str=0;str<numberOfStrands;str++){
                    
                    StrandRCs sr = (StrandRCs)strandRot[str];

                    //Could be a problem here!  Ligand needs storage regardless because of rotation/translation!!
                    //Should handle this by reversing translation/rotation
                    if(addWTRot)
                        sr.storeWTRotamers(m);
                    
                    sr.addUnperturbedRCs(addWTRot);
                    sr.countRCs();
                }

        }

	
	// This function adds the AA type named name to the list
	//  of allowable types for residue number resNum in
	//  the system strand (resNum is strand based numbering)
	public void setAllowable(int resNum, String name, int strNum) {
		strandRot[strNum].setAllowable(resNum,name);
		m.strand[strNum].residue[resNum].flexible = true;
	}

        
	public void setSplitFlags(){

		if(arpMatrix == null){
			System.err.println("setSplitFlags was using the arpMatrix but it is null");
			System.exit(0);
		}

                splitFlags = arpMatrix.initializePairwiseBooleanMatrix();
	}
	
	public void setSplitFlags(boolean fromMatrix[][][][][][]){
		
		boolean toMatrix[][][][][][] = new boolean[fromMatrix.length][][][][][];
		if(fromMatrix == null){
			System.err.println("setSplitFlags was given a null array");
			System.exit(0);
		}
		for (int p1=0; p1<fromMatrix.length; p1++){
			if (fromMatrix[p1]!=null){
				toMatrix[p1] = new boolean[fromMatrix[p1].length][][][][];
				for (int a1=0; a1<fromMatrix[p1].length; a1++){
					if (fromMatrix[p1][a1]!=null){
						toMatrix[p1][a1] = new boolean[fromMatrix[p1][a1].length][][][];
						for (int r1=0; r1<fromMatrix[p1][a1].length; r1++){
							if (fromMatrix[p1][a1][r1]!=null){
								toMatrix[p1][a1][r1] = new boolean[fromMatrix[p1][a1][r1].length][][];
								for (int p2=0; p2<fromMatrix[p1][a1][r1].length; p2++){
									if (fromMatrix[p1][a1][r1][p2]!=null){
										toMatrix[p1][a1][r1][p2] = new boolean[fromMatrix[p1][a1][r1][p2].length][];
										for (int a2=0; a2<fromMatrix[p1][a1][r1][p2].length; a2++){
											if (fromMatrix[p1][a1][r1][p2][a2]!=null){
												toMatrix[p1][a1][r1][p2][a2] = new boolean[fromMatrix[p1][a1][r1][p2][a2].length];
												for (int r2=0; r2<fromMatrix[p1][a1][r1][p2][a2].length; r2++){
													toMatrix[p1][a1][r1][p2][a2][r2] = fromMatrix[p1][a1][r1][p2][a2][r2];
												}
												//System.arraycopy(fromMatrix[p1][a1][r1][p2][a2], 0, toMatrix[p1][a1][r1][p2][a2], 0, fromMatrix[p1][a1][r1][p2][a2].length);
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
		
		splitFlags = toMatrix;
		/*for (int i=0; i<spFlags.length; i++){
			for (int j=0; j<spFlags.length; j++){
				splitFlags[i][j] = spFlags[i][j];
			}
		}*/
	}
	
	//Return the split flags
	//	This version makes sure that the returned matrix is not associated with the RotamerSearch instance variable
	public boolean[][][][][][] getSplitFlags(boolean toMatrix[][][][][][]){
		
		toMatrix = new boolean[splitFlags.length][][][][][];
		 
		for (int p1=0; p1<splitFlags.length; p1++){
			if (splitFlags[p1]!=null){
				toMatrix[p1] = new boolean[splitFlags[p1].length][][][][];
				for (int a1=0; a1<splitFlags[p1].length; a1++){
					if (splitFlags[p1][a1]!=null){
						toMatrix[p1][a1] = new boolean[splitFlags[p1][a1].length][][][];
						for (int r1=0; r1<splitFlags[p1][a1].length; r1++){
							if (splitFlags[p1][a1][r1]!=null){
								toMatrix[p1][a1][r1] = new boolean[splitFlags[p1][a1][r1].length][][];
								for (int p2=0; p2<splitFlags[p1][a1][r1].length; p2++){
									if (splitFlags[p1][a1][r1][p2]!=null){
										toMatrix[p1][a1][r1][p2] = new boolean[splitFlags[p1][a1][r1][p2].length][];
										for (int a2=0; a2<splitFlags[p1][a1][r1][p2].length; a2++){
											if (splitFlags[p1][a1][r1][p2][a2]!=null){
												toMatrix[p1][a1][r1][p2][a2] = new boolean[splitFlags[p1][a1][r1][p2][a2].length];
												for (int r2=0; r2<splitFlags[p1][a1][r1][p2][a2].length; r2++){
													toMatrix[p1][a1][r1][p2][a2][r2] = splitFlags[p1][a1][r1][p2][a2][r2];
												}
												//System.arraycopy(fromMatrix[p1][a1][r1][p2][a2], 0, toMatrix[p1][a1][r1][p2][a2], 0, fromMatrix[p1][a1][r1][p2][a2].length);
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
		/*savedSpF = new boolean[splitFlags.length][splitFlags.length];
		for (int i=0; i<splitFlags.length; i++){
			for (int j=0; j<splitFlags.length; j++){
				savedSpF[i][j] = splitFlags[i][j];
			}
		}
		return savedSpF;*/
	}
	
	public boolean[][][][][][] getSplitFlags(){
		return splitFlags;
	}

	public boolean[][][][][][][][][] getTripleFlags(){
		return tripleFlags;
	}

	public synchronized float getBestE(){
		return bestEMin;
	}
	
	public synchronized void setBestE(float newBestE){
		bestEMin = newBestE;
	}
	
	//updates bestEMin only if (ue<bestEMin)
	public synchronized void updateBestE(float ue){
		bestEMin = Math.min(bestEMin, ue);
	}
	
	// This function clears the list of allowable AA types
	//  for residue number resNum in the system strand
	public void clearAllowable(int resNum, int strandNum) {
		strandRot[strandNum].clearAllowable(resNum);
	}

	// Refreshes the system strand
	public void refreshStrand(int str){
		strandRot[str] = strandRot[str].reInit(m.strand[str]);
	}
	
	public double getEc_const(){
		return Ec_const;
	}
	
	public double getBoundForPartition(){
		return boundForPartition;
	}

	// This function computes one energy 
	private float calcTotalSnapshotEnergy(){
		
		double energyTerms[] = a96ff.calculateTotalEnergy(m.actualCoordinates,-1); //compute the energy
	
		return (float)energyTerms[0]; //the total energy is in energyTerms[0]
	}

//// BEGIN CHECK_STERICS CODE SECTION
	//This version checks all residues against residue resNum (strand-relative numbering) of strand strandNum;
	//This function is only called *before* minimization for the pairwise matrix energy computation;
	//The residue numbers (strand-relative numbering) in excludeRes[] (from the system strand) are not included in the steric check
		private boolean RS_CheckAllSterics(int strandNum, int resNum, int excludeList[][]) {
		
		ProbeStericCheck psc = new ProbeStericCheck();
	
		Residue res = m.strand[strandNum].residue[resNum];
		
		for(int i=0;i<res.numberOfAtoms;i++) {
			Atom a1 = res.atom[i];
			if ( hSteric || (!a1.elementType.equalsIgnoreCase("H")) ) {
				for(int q=0;q<m.numberOfStrands;q++) {
					int resToCheck = m.strand[q].numberOfResidues;
					for(int w=0;w<resToCheck;w++) {
						if(!((q==strandNum) && (w==resNum))) {
							boolean resInList = false;
							//if (q==sysStrNum) //check for excluded residues in the system strand
							resInList = isInList(w,excludeList[q]);
							for(int t=0;t<m.strand[q].residue[w].numberOfAtoms;t++) {
								Atom a2 = m.strand[q].residue[w].atom[t];
								if ( !resInList || ( a2.getIsBBatom() && (!doPerturbations) ) ) {
                                                                    //For methods other than DEEPer, check only the backbone atoms for the sysStrand residues that are in excludeList[]
                                                                    //for DEEPer, these backbone atoms may move a lot, so exclude residues in excludeList entirely
									if ( hSteric || (!a2.elementType.equalsIgnoreCase("H"))) {
										if (!psc.isAllowedSteric(m, a1, a2, softOverlapThresh))
											return false;
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
	//Checks if a is in list[]
	private boolean isInList(int a, int list[]){
		for (int i=0; i<list.length; i++){
			if (list[i]==a)
				return true;
		}
		return false;
	}

//// END CHECK_STERICS CODE SECTION

//// BEGIN HELPER FUNCTION SECTION
	
	// Loads the min (minMatrix==true) or max (minMatrix==false) pairwise energy matrix
	public void loadPairwiseEnergyMatrices(String allRotamerPairsEnergyName, boolean minMatrix) {


		try{
			ObjectInputStream in = new ObjectInputStream(new FileInputStream(allRotamerPairsEnergyName));
			if (minMatrix){
                            arpMatrix = new PairwiseEnergyMatrix();
                            arpMatrix.eMatrix = (float [][][][][][])in.readObject();
                        }
			else{
                            arpMatrixMax = new PairwiseEnergyMatrix();
                            arpMatrixMax.eMatrix = (float [][][][][][])in.readObject();
                        }
			in.close();
		}
		catch (Exception e){}		
	}
	
	public PairwiseEnergyMatrix getMinMatrix(){
		return arpMatrix;
	}
	
	public PairwiseEnergyMatrix getMaxMatrix(){
		return arpMatrixMax;
	}
	
	//Adds the reference energies to the intra-energies in arpMatrix;
	//If doMinimize is true, then arpMatrixMax is also updated appropriately
	//Adds the reference energies to the intra-energies in arpMatrix;
	//If doMinimize is true, then arpMatrixMax is also updated appropriately
	public void addEref(float eRef[][], boolean doMinimize, int strandMut[][]){
		
		int ind = 1; //skip the entry [0][0], since this is the fixed template energy
		
		int ctr=0;
		for (int str=0;str<numberOfStrands;str++){
			for(int i=0;i<strandMut[str].length;i++){
				for(int j=0; j<strandRot[str].getNumAllowable(strandMut[str][i]);j++){
					int aaInd = strandRot[str].getIndexOfNthAllowable(strandMut[str][i],j);
					int numRot = getNumRot(str, strandMut[str][i], aaInd);

					for (int k=0; k<numRot; k++){
							arpMatrix.addToIntraE(ctr,aaInd,k, -eRef[ctr][aaInd]);
						if (doMinimize)
                                                        arpMatrixMax.addToIntraE(ctr,aaInd,k, -eRef[ctr][aaInd]);
						ind++;
					}
				}			
				ctr++;
			}		
		}
	}

	// Computes the best energy (lower bound) using the arpMatrix
	// This energy is rotamer based, that is it computes the best energy
	//  for the current rotamer assignment of each amino-acid	
	public float computeBestRotEnergyBound(/*int numTotalRotamers, int rotamerIndexOffset[]*/) {

		float bestE = arpMatrix.getShellShellE(); // Add shell-shell energy
	
		for(int i=0;i<ASAANums.length;i++) {
			bestE += arpMatrix.getIntraAndShellE(i, ASAANums[i], curStrRotNum[i]); //Add the intra-rotamer and shell energies for each rotamer
			for(int j=i+1;j<ASAANums.length;j++){ // Add the pairwise energies
				bestE += arpMatrix.getPairwiseE( i, ASAANums[i], curStrRotNum[i], j, ASAANums[j], curStrRotNum[j] );
			}
		}

		return bestE;
	}
//// END HELPER FUNCTION SECTION

	
//// BEGIN MASTER MUTATION SEARCH SECTION

	// This function is the model for a mutation search
	// This is the function used by the master node to generate a list of
	//  mutations that it wishes to consider.
	// Utilizes a number of helper functions
	public int simpleMasterMutationSearch(int strandMut[][], int numMutable,
		int theCurConfNum,
		Set<OneMutation> mutArray, float minVol, float maxVol) {

		curConfNum = theCurConfNum;
		
		String curAAnames[] = new String[numMutable];
		float rotamerVolumes[][][] = new float[strandMut.length][][];
		//TODO: rl.getRotVol will fail for non-amino acid rotamers
		for(int str=0;str<strandMut.length;str++)
			rotamerVolumes[str] = strandRot[str].rl.getRotVol();
		
		initMutRes2Str(strandMut);
		
		masterMutationSearchHelper(0, numMutable, strandMut, mutArray, minVol, 
				maxVol, curAAnames, rotamerVolumes);
		
		return curConfNum;
	}
	// This function is similar to mutationSearchHelper
	//  the only difference is that we only compute volumes and an amino acid
	//  level energy approximation. I could have modified that function, but
	//  decided not to so as to keep that function fast (ie. this way the
	//  execution of a bunch of conditionals is saved in the normal search)
	public void masterMutationSearchHelper(int depth, int maxDepth,
		int strandMut[][], Set<OneMutation> mutSet, float minVol, float maxVol, String curAAnames[], 
		float rotamerVolumes[][][]) {
	
		if (depth >= maxDepth) {
			// If we've arrived here then we're ready to
			//  compute a volume and approximate energy
			if(debug){
				System.out.print(".");
			}
			float curVolume = 0.0f;
			for (int i=0; i<maxDepth; i++){
				int str = mutRes2Strand[i];
				curVolume += rotamerVolumes[str][strandRot[str].rl.getAARotamerIndex(curAAnames[i])][0]; //use the rotamer with index 0
			}
			if ((curVolume > minVol) && (curVolume < maxVol)) {
				// Add mutation to mutation array
				OneMutation tMut = new OneMutation();
				assignAANums(strandMut);
				tMut.score = new BigDecimal("0.0");  // Added when aap removed 6/23/03
				tMut.resTypes = new String[maxDepth];
				for(int q=0;q<maxDepth;q++) {
					tMut.resTypes[q] = curAAnames[q];
				}
				tMut.vol = curVolume;
				/*if (curConfNum >= mutArray.length) {
					// If there's no space left, make space in mutArray
					OneMutation newArray[] = new OneMutation[mutArray.length + 5000];
					System.arraycopy(mutArray, 0, newArray, 0, mutArray.length);
					mutArray = newArray;
				}*/
				getNumConfForMut(tMut);
				mutSet.add(tMut);
				curConfNum++;
			}
			return;
		}

		// Check with allowed AAs
		int str = mutRes2Strand[depth];
		int strResNum = strandMut[str][mutRes2StrandMutIndex[depth]];
		
		
		for(int q=0;q<strandRot[str].getNumAllowable(strResNum);q++) {
			curAAnames[depth] = strandRot[str].rl.getAAName(strandRot[str].getIndexOfNthAllowable(strResNum,q));
			masterMutationSearchHelper(depth+1,maxDepth,strandMut,mutSet,minVol,maxVol,curAAnames,rotamerVolumes);
		}
	}
	// Assigns elements of the ASAANums[] array
	private void assignAANums(int strandMut[][]) {
		
		ASAANums = new int[numberMutable];
		int ctr=0;
		for(int str=0;str<numberOfStrands;str++){
			for(int i=0;i<strandMut[str].length;i++){
				ASAANums[ctr] = strandRot[str].rl.getAARotamerIndex(m.strand[str].residue[strandMut[str][i]].name);
				ctr++;
			}
		}
	}
	//Compute the number of conformations for the given mutation sequence
	private void getNumConfForMut(OneMutation tMut){
		tMut.numConfUB = BigInteger.ONE;
		tMut.numConfB = BigInteger.ONE;
		for (int i=0; i<tMut.resTypes.length; i++){
			int str = mutRes2Strand[i];
			int numRot = strandRot[str].rl.getNumRotamers(tMut.resTypes[i]);
			if (numRot==0)
				numRot = 1;
			tMut.numConfUB = tMut.numConfUB.multiply(BigInteger.valueOf(numRot));
			tMut.numConfB = tMut.numConfUB;
		}
		/*int numLigRot = grl.getNumRotamers(ligROT.getCurRotType(0));
		if (numLigRot==0)
			numLigRot = 1;
		tMut.numConfB = tMut.numConfB.multiply(BigInteger.valueOf(numLigRot));*/
	}
//// END MASTER MUTATION SEARCH SECTION
///////////////////////////////////////////////////////////////////////////////////	

///////////////////////////////////////////////////////////////////////////////////
//// BEGIN PAIRWISE MUTATION _ALL_ ROTAMER SEARCH SECTION - ENERGY PRECOMPUTATION
	
	// Turns all residues off except the side-chain for the residue specified by
	//  molResNum. Computes the energy of the current system.
	//  It does not do minimization.
	// *Make sure that the forcefield types for the atoms in the
	//  residue of interest have been computed prior to calling
	//  this function as they are not computed in this function
	private float computeEnergyOfOnlyRes(int molResNum, int pos, int[][] strandMut, Amber96ext resA96ff) {
		//Get aa (int index) of the cur molResNum
		int AAind = strandRot[m.residue[molResNum].strandNumber].rl.getAARotamerIndex(m.residue[molResNum].name);
		
		if(resA96ff == null){ //create the template amber FF if it's null
			if(pos==0){
				res1A96ff = new Amber96ext(m, distDepDielect, dielectConst, doSolvationE, solvScale, vdwMultiplier);
                resA96ff = res1A96ff;
            }
            else if(pos==1){
            	res2A96ff = new Amber96ext(m, distDepDielect, dielectConst, doSolvationE, solvScale, vdwMultiplier);
                resA96ff = res2A96ff;
            }
			
			
			boolean savedEnergyEvalSC[] = new boolean[m.numberOfResidues];
			boolean savedEnergyEvalBB[] = new boolean[m.numberOfResidues];
				
			// Save the energy eval flag, clear them at the same time
			for(int i=0;i<m.numberOfResidues;i++){
				savedEnergyEvalSC[i] = m.residue[i].getEnergyEvalSC();
				savedEnergyEvalBB[i] = m.residue[i].getEnergyEvalBB();
				m.residue[i].setEnergyEval(false, false);
			}
							
            m.residue[molResNum].setEnergyEval(true, true);
			
			resA96ff.initializeCalculation();
			resA96ff.setNBEval(hElect,hVDW);
			
			// Restore the energy eval and flexibility flags
			for(int i=0;i<m.numberOfResidues;i++){
				m.residue[i].setEnergyEval(savedEnergyEvalSC[i], savedEnergyEvalBB[i]);
			}
			
		}
				
		double energyTerms[] = resA96ff.calculateTotalEnergy(m.actualCoordinates,-1); //compute the energy
		float curEnergy = (float) energyTerms[0];
			
		return curEnergy;
	}
	//Turns all template residues on and the side-chains for all flexible residues off;
	//Computes the template energy for the current template; does not do minimization
	// *Make sure that the forcefield types for the atoms in the
	//  residue of interest have been computed prior to calling
	//  this function as they are not computed in this function
	private float computeEnergyOfOnlyTemplate(int strandMut[][]) {
		if(templateA96ff == null){ //create the template amber FF if it's null
			templateA96ff = new Amber96ext(m, distDepDielect, dielectConst, doSolvationE, solvScale, vdwMultiplier);
			
			boolean savedEnergyEvalSC[] = new boolean[m.numberOfResidues];
			boolean savedEnergyEvalBB[] = new boolean[m.numberOfResidues];
			
				
			// Save the energy eval flag, clear them at the same time
			for(int i=0;i<m.numberOfResidues;i++){
				savedEnergyEvalSC[i] = m.residue[i].getEnergyEvalSC();
				savedEnergyEvalBB[i] = m.residue[i].getEnergyEvalBB();
				m.residue[i].setEnergyEval(true, true);
			}
							
			//Clear the energy evaluation flags for the AS residues and the ligand (if present), so that
			//		the energies only between template residues are computed
			for (int str=0; str<strandMut.length;str++){
				for (int i=0; i<strandMut[str].length; i++){
					if(isTemplateOn)
						m.strand[str].residue[strandMut[str][i]].setEnergyEval(false, false);
					else
						m.strand[str].residue[strandMut[str][i]].setEnergyEval(false, false);
				}
			}
			
			templateA96ff.initializeCalculation();
			templateA96ff.setNBEval(hElect,hVDW);
			
			// Restore the energy eval and flexibility flags
			for(int i=0;i<m.numberOfResidues;i++){
				m.residue[i].setEnergyEval(savedEnergyEvalSC[i], savedEnergyEvalBB[i]);
			}
			
		}
		
		double energyTerms[] = templateA96ff.calculateTotalEnergy(m.actualCoordinates,-1); //compute the energy
		float curEnergy = (float) energyTerms[0];
			
		return curEnergy;
	}
	
	public void setTemplateOff(int[][] strandMut, int res1, int res2, boolean shellRun, boolean templateOnly){
		templateOffA96ff = new Amber96ext(m, distDepDielect, dielectConst, doSolvationE, solvScale, vdwMultiplier);
		
		boolean savedEnergyEvalSC[] = new boolean[m.numberOfResidues];
		boolean savedEnergyEvalBB[] = new boolean[m.numberOfResidues];
		
		
			
		// Save the energy eval flag, clear them at the same time
		for(int i=0;i<m.numberOfResidues;i++){
			savedEnergyEvalSC[i] = m.residue[i].getEnergyEvalSC();
			savedEnergyEvalBB[i] = m.residue[i].getEnergyEvalBB();
			m.residue[i].setEnergyEval(false, false);
		}
						
		//Clear the energy evaluation flags for the AS residues and the ligand (if present), so that
		//		the energies only between template residues are computed
		
		/*if (ligStrNum>=0)
			m.strand[ligStrNum].residue[0].setEnergyEval(false, false);*/
		
		
		if(!shellRun){
			int str2 = mutRes2Strand[res2];
			int strResNum2 = strandMut[str2][mutRes2StrandMutIndex[res2]];
			m.strand[str2].residue[strResNum2].setEnergyEval(true,true);
		}
		else { //(shellRun){
			//SHL-AS run needs SHL on
			for(int i=0;i<m.numberOfResidues;i++){
				m.residue[i].setEnergyEval(true,true);
			}
			for (int str=0; str<strandMut.length;str++)
			for (int i=0; i<strandMut[str].length; i++){
				if(isTemplateOn)
					m.strand[str].residue[strandMut[str][mutRes2StrandMutIndex[i]]].setEnergyEval(false, false);
				else
					m.strand[str].residue[strandMut[str][mutRes2StrandMutIndex[i]]].setEnergyEval(false, false);
			}
		}
		
		if(!templateOnly){
			int str1 = mutRes2Strand[res1];
			int strResNum1 = strandMut[str1][mutRes2StrandMutIndex[res1]];
			//Do this after so the above function doesn't make it false.
			//if((ligPresent && ! shellRun) || !ligPresent)
		
			m.strand[str1].residue[strResNum1].setEnergyEval(true,true);
			//if(ligPresent)
			//	m.strand[ligStrNum].residue[0].setEnergyEval(true, true);
		}
		
		templateOffA96ff.initializeCalculation();
		templateOffA96ff.setNBEval(hElect,hVDW);
		
		// Restore the energy eval and flexibility flags
		for(int i=0;i<m.numberOfResidues;i++){
			m.residue[i].setEnergyEval(savedEnergyEvalSC[i], savedEnergyEvalBB[i]);
		}
	}

	
	// This function helps to compute the min/max pairwise energy matrix, two
	//  residues are allowed to mutate (as specified in residueMutatable)
	//  a steric check is done for each rotamer pair if the steric check
	//  passes then the min/max energy is computed for the pair and is
	//  saved. ALL rotamer pairs that pass the steric threshold are
	//  saved. If a pair doesn't pass the threshold then an energy of
	//  10^38 is assigned.
	// To compute the pairwise interactions of each rotameric position with
	//  the ligand only one residue should be "allowed" in residueMutatable[]
	//  and the ligPresent should be set to true.
	// A shellRun computes the energy of the "allowed" residue with all other
	//  residues that are NOT in strandMut
	// Utilizes a number of helper functions
	public void simplePairwiseMutationAllRotamerSearch(int strandMut[][], int mutableSpots, 
			boolean searchDoMinimize, boolean shellRun, boolean intraRun, 
			int residueMutatable[], PairwiseEnergyMatrix retEMatrixMin, PairwiseEnergyMatrix retEMatrixMax,
			boolean minimizeBB, boolean doBackrubs, boolean templateOnly, String backrubFile, boolean templateAlwaysOn) {


		doMinimization = searchDoMinimize;
		computeEVEnergy = true;

                if(doPerturbations)
                    ((PMinimizer)simpMin).minimizePerturbations=minimizePerturbations;//This field needs to be set properly when this function is called

		// Prepare Amber
		if(computeEVEnergy){
			// Amber should already be loaded
			// First turn off energy evaluation for all residues
			//   since we're only computing pairwise energies
			//   we only want specific residues on (will be turned on later)
			// If we're doing a shell run then turn the shell on
			if(!templateAlwaysOn){
				if (shellRun) {
					for(int i=0;i<m.numberOfResidues;i++){
						m.residue[i].setEnergyEval(true, true);
						m.residue[i].flexible = false;
					}
					for(int i=0;i<strandMut.length;i++){
						for(int j=0;j<strandMut[i].length;j++){
							m.strand[i].residue[strandMut[i][j]].setEnergyEval(false, false);
							m.strand[i].residue[strandMut[i][j]].flexible = false;
						}
					}
				}
				else {
					for(int i=0;i<m.numberOfResidues;i++){
						m.residue[i].setEnergyEval(false, false);
						m.residue[i].flexible = false;
					}
				}
				isTemplateOn = false;
			}
			else{
				System.out.println("Template on during minimization of bounds (this might take longer)");
				// 2010: If templateAlwaysOn then always have the shell on for pairwise calculations.
				for(int i=0;i<m.numberOfResidues;i++){
					m.residue[i].setEnergyEval(true,true);
					m.residue[i].flexible = false;
				}
				for(int i=0;i<strandMut.length;i++){
					for(int j=0;j<strandMut[i].length;j++){
                                            if(doPerturbations)
                                                m.strand[i].residue[strandMut[i][j]].setSCEnergyEval(false);
                                            else
						m.strand[i].residue[strandMut[i][j]].setEnergyEval(false,false);
                                            m.strand[i].residue[strandMut[i][j]].flexible = false;
					}
				}
				isTemplateOn = true;
			}
						
			if (doMinimization){
				if (!minimizeBB){ //side-chain minimization
					if (simpMin == null) {
						System.out.println("Error: simpMin not allocated and you are attempting minimization");
						System.exit(1);
					}
					bbMin = null;
					brMin = null;
				}
				else { //backbone minimization
					if (!doBackrubs){ // phi/psi minimization
						if (bbMin == null) {
							System.out.println("Error: bbMin not allocated and you are attempting backbone minimization");
							System.exit(1);
						}
						simpMin = null;
						brMin = null;
					}
					else { //minimization with backrubs
						if (brMin == null) {
							System.out.println("Error: brMin not allocated and you are attempting backrub minimization");
							System.exit(1);
						}
						simpMin = null;
						bbMin = null;
					}
				}
			}
                }


                        if (shellRun) { //compute the template energies
				if ((!minimizeBB)&&m.perts.length>0) {//side-chain minimization, so the template is fixed


                                    //The full structure switch perturbation can change the template
                                    //If we have one, it'll be in m.perts[0], and we need to calculate a template energy for each of its states
                                   if( m.perts[0].type.equalsIgnoreCase("FULL STRUCTURE SWITCH") ){
                                       
                                        int numTemplates = ((FullStructureSwitch)m.perts[0]).numStructs;
                                        switchedTemplateE = new float[numTemplates];

                                        for(int tempNum=0; tempNum<numTemplates; tempNum++){

                                            if(tempNum>0)
                                                m.perts[0].applyPerturbation(tempNum);

                                            for(int str=0; str<numberOfStrands; str++){
                                                for (int i=0; i<strandMut[str].length; i++)
                                                    m.strand[str].residue[strandMut[str][i]].setEnergyEval(false,false);
                                            }
                                            a96ff.calculateTypesWithTemplates();
                                            a96ff.initializeCalculation();
                                            a96ff.setNBEval(hElect,hVDW);
                                            float minE = calcTotalSnapshotEnergy();

                                            if(tempNum==0){
                                                retEMatrixMin.setShellShellE(minE);
                                                retEMatrixMax.setShellShellE(minE);
                                            }
                                            else
                                                switchedTemplateE[tempNum] = minE - retEMatrixMin.getShellShellE();
                                        }
                                    }
				}
			}
			
			if (templateOnly) { //the template energies are computed only once	
				
				if (!minimizeBB && !doMinimization) {//No minimization, so the template is fixed
					for(int i=0;i<strandMut.length;i++){
						for(int j=0;j<strandMut[i].length;j++){
											
							m.strand[i].residue[strandMut[i][j]].flexible = false;
							if(isTemplateOn)
								m.strand[i].residue[strandMut[i][j]].setEnergyEval(false,false);
							else
								m.strand[i].residue[strandMut[i][j]].setEnergyEval(false,false);
						}
					}
					
					a96ff.calculateTypesWithTemplates();
					a96ff.initializeCalculation();
					a96ff.setNBEval(hElect,hVDW);
					
					float minE = calcTotalSnapshotEnergy();
					retEMatrixMin.setShellShellE( minE );
					retEMatrixMax.setShellShellE( minE );
				}
				if (!minimizeBB) {//side-chain minimization, so the template is fixed	
					for(int i=0; i<strandMut.length; i++){
						for(int j=0;j<strandMut[i].length;j++){
							String a = strandRot[i].getCurRotType(strandMut[i][j]);
							if ((!a.equalsIgnoreCase("GLY"))&&(!a.equalsIgnoreCase("PRO"))){
								String tmpAA = "gly";
								//if (minimizeBB)
								//	tmpAA = "ala";
								if(m.strand[i].isProtein)
									strandRot[i].changeResidueType(m,strandMut[i][j],tmpAA,addHydrogens,connectResidues);
							}
							
							m.strand[i].residue[strandMut[i][j]].flexible = false;
							if(isTemplateOn)
								m.strand[i].residue[strandMut[i][j]].setEnergyEval(false,false);
							else
								m.strand[i].residue[strandMut[i][j]].setEnergyEval(false,false);	
						}
					}
					a96ff.calculateTypesWithTemplates();
					a96ff.initializeCalculation();
					a96ff.setNBEval(hElect,hVDW);
					
					simpMin.initialize(m,numberOfStrands,a96ff,strandRot,curAANum,doDihedE);

					float beginE = computeEnergyOfOnlyTemplate(strandMut);//calcTotalSnapshotEnergy();
					simpMin.minimize(35);
					float curE = computeEnergyOfOnlyTemplate(strandMut); //calcTotalSnapshotEnergy();
					
					float minE = Math.min(beginE, curE);
					
					retEMatrixMin.setShellShellE(minE);
					retEMatrixMax.setShellShellE(minE);
				}
				else if (minimizeBB){ //the template energies for backbone minimization are computed only once
					for(int i=0; i<strandMut.length; i++){
						for(int j=0;j<strandMut[i].length;j++){
							// Make this residue a "gly", if it is not already GLY or PRO;
							//		If minimizeBB, then change to ALA, if not already GLY or PRO
							String a = strandRot[i].getCurRotType(strandMut[i][j]);
							if ((!a.equalsIgnoreCase("GLY"))&&(!a.equalsIgnoreCase("PRO"))){
								String tmpAA = "gly";
								//if (minimizeBB)
								//	tmpAA = "ala";
							
								if(m.strand[i].isProtein)
									strandRot[i].changeResidueType(m,strandMut[i][j],tmpAA,addHydrogens,connectResidues);
							}
							
							m.strand[i].residue[strandMut[i][j]].flexible = false;
							if(isTemplateOn)
								m.strand[i].residue[strandMut[i][j]].setEnergyEval(false,false);
							else
								m.strand[i].residue[strandMut[i][j]].setEnergyEval(false,false);
						}
					}
					a96ff.calculateTypesWithTemplates();
					a96ff.initializeCalculation();
					a96ff.setNBEval(hElect,hVDW);
					if (!doBackrubs){
						setTemplateOff(strandMut, -1, -1, shellRun, templateOnly);
						bbMin.initialize(m, a96ff, strandMut, numberOfStrands);
						pairwiseRotamerEnergyBackboneHelper(-1,-1,-1,-1,-1,-1,strandMut,shellRun, retEMatrixMin, retEMatrixMax, true,-1,-1,-1,-1);					
						bbMin = new BBMinimizer(); //reset the backbone minimizer
					}
					else {
						setTemplateOff(strandMut, -1, -1, shellRun, templateOnly);
						brMin.initialize(m, a96ff, strandMut, backrubFile, hSteric, overlapThresh, numberOfStrands,true); //no ligand for template energy computation					
						pairwiseRotamerEnergyBackrubsHelper(-1,-1,-1,-1,-1,-1,strandMut, shellRun, retEMatrixMin, retEMatrixMax, true, -1,-1,-1,-1);					
						brMin = new BackrubMinimizer(); //reset the backbone minimizer
					}
				}
			}

		else if (intraRun) { //INTRA run
			computeIntraRotEnergies(strandMut,retEMatrixMin,retEMatrixMax,residueMutatable,minimizeBB,doBackrubs, backrubFile);			
			return;
		}
		
		else {	//pairwise rotamer or rot-shell run	
			// Initialize curAANum array, we have to do this because we only recurse through the 9 core residues of the active site
			//  and not all 40 residues, so we use a residueMap and we have to do some preinitialization.
			for(int i=0;i<m.numberOfResidues;i++)
				curAANum[i] = -1;
			
			// Note: In this search we only search over the key active site  residues rather than all residues in the molecule, 
			//		thus maxDepth should be numInAS
			initMutRes2Str(strandMut);
						
			pairwiseEnergyComp (mutableSpots, strandMut, residueMutatable, retEMatrixMin, retEMatrixMax, shellRun, minimizeBB, doBackrubs, backrubFile);

                        if( ( switchedTemplateE != null ) && ( residueMutatable[0] != 0 ) )//This transfers energies from template changes in full structure switch perturbations to the first residue shell energies
                            correctSwitchedTemplates(strandMut, retEMatrixMin, retEMatrixMax);

			return;
		}
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
			
	private void computeIntraRotEnergies(int strandMut[][],	PairwiseEnergyMatrix retEMatrixMin,
			PairwiseEnergyMatrix retEMatrixMax, int residueMutatable[],
			boolean minimizeBB, boolean doBackrubs, String backrubFile) {
	
		//Save the energy eval and flexibilty flags, clear them at the same time
		/*boolean savedEnergyEvalSC[] = new boolean[m.numberOfResidues];
		boolean savedEnergyEvalBB[] = new boolean[m.numberOfResidues];
		boolean savedFlexible[] = new boolean[m.numberOfResidues];
		for(int i=0;i<m.numberOfResidues;i++){
			savedEnergyEvalSC[i] = m.residue[i].getEnergyEvalSC();
			savedEnergyEvalBB[i] = m.residue[i].getEnergyEvalBB();
			savedFlexible[i] = m.residue[i].flexible;
			m.residue[i].setEnergyEval(false, false);
			m.residue[i].flexible = false;
		}*/
		
		// Go through each active site residue, each AA type they could be and all their rotamers, 
		// 	saving the computed energies to the appropriate place.
		int curPos;
		initMutRes2Str(strandMut);
		for(int i=0; i<residueMutatable.length;i++){
			//for(int str=0; str<numberOfStrands;str++){
			//for(int i=0;i<strandMut[str].length;`i++) { //i is strand residue specific
			if (residueMutatable[i] == 1){ //only compute for those residues that are mutable
				curPos = i;
				//curPos++;
				System.out.println();
				
				int str = mutRes2Strand[i];
				int strResNum = strandMut[str][mutRes2StrandMutIndex[i]];
				int molResNum = m.strand[str].residue[strResNum].moleculeResidueNumber;
				//int strResNum = strandMut[str][i];
					
				for(int j=0;j<strandRot[str].getNumAllowable(strResNum);j++){
					
					System.out.print(".");
				
					// Apply mutation
					curAANum[molResNum] = strandRot[str].getIndexOfNthAllowable(strResNum,j);
					if(m.strand[str].isProtein)
						strandRot[str].changeResidueType(m,strResNum,strandRot[str].rl.getAAName(curAANum[molResNum]),addHydrogens,connectResidues);
					m.strand[str].residue[strResNum].flexible = true;
					
					if(isTemplateOn&&(!doPerturbations))
						m.strand[str].residue[strResNum].setSCEnergyEval(true);
					else
						m.strand[str].residue[strResNum].setEnergyEval(true,true);
						
					// Setup Amber, Setup Minimizer
					if (computeEVEnergy){
						a96ff.calculateTypesWithTemplates();
						a96ff.initializeCalculation();
						a96ff.setNBEval(hElect,hVDW);
						if (doMinimization){
							if (!minimizeBB) //side-chain minimization
									simpMin.initialize(m,numberOfStrands,a96ff,strandRot,curAANum,doDihedE);
							else { //backbone minimization
								if (!doBackrubs)
										bbMin.initialize(m, a96ff, strandMut, numberOfStrands);
									else{
										brMin.initialize(m, a96ff, strandMut, backrubFile, hSteric, overlapThresh, numberOfStrands,true);
									}
								}
						}
					}
	
					
					boolean done = false;
					int curRot = 0;
					int totRotForCur;
                                        if(doPerturbations)
                                            totRotForCur = ((StrandRCs)strandRot[str]).getNumRCs(strResNum, curAANum[molResNum]);
                                        else
                                            totRotForCur = strandRot[str].rl.getNumRotForAAtype(curAANum[molResNum]);
					
					while ( (!done) && ((curRot<totRotForCur) || (totRotForCur==0)) ){
					
						if (totRotForCur!=0){ // loop through rotamers
							if(doPerturbations)
                                                            ((StrandRCs)strandRot[str]).applyRC(m, strResNum, curRot);
                                                        else
                                                            strandRot[str].applyRotamer(m, strResNum, curRot);
                                                }
						else // no rotamers, so only the current AA state has to be computed and we exit the while loop
							done = true;
						
							computeIntraRotEnergiesHelper(curPos, strResNum, curAANum[molResNum], curRot, str, totRotForCur, 
								retEMatrixMin, retEMatrixMax, strandMut, minimizeBB, doBackrubs);
						
						curRot++;
					}				
				}
					m.strand[str].residue[strResNum].flexible = false;
					if(isTemplateOn&&(!doPerturbations))
						m.strand[str].residue[strResNum].setSCEnergyEval(false);
					else
						m.strand[str].residue[strResNum].setEnergyEval(false,false);
			}
		}
		
		/*if (ligPresent) { //go through all of the ligand rotamers and save the computed intra energies
			curLigAANum = ligROT.getIndexOfNthAllowable(0,0);
			if (m.strand[ligStrNum].isProtein) // apply mutation
				ligROT.changeResidueType(m,0,grl.getAAName(curLigAANum),addHydrogens);
			m.strand[ligStrNum].residue[0].flexible = true;
			m.strand[ligStrNum].residue[0].setEnergyEval(true,true);

			// Setup Amber, Setup Minimizer
			if (computeEVEnergy){
				a96ff.calculateTypesWithTemplates();
				a96ff.initializeCalculation();
				a96ff.setNBEval(hElect,hVDW);
				if (doMinimization){
					if (!minimizeBB) //side-chain minimization
						simpMin.initialize(m,sysStrNum,ligStrNum,a96ff,sysLR,ligROT,curAANum,curLigAANum,doDihedE,rl,grl);
					else { //backbone minimization
						if (!doBackrubs)
							bbMin.initialize(m, a96ff, residueMap, sysStrNum, ligStrNum);
						else
							brMin.initialize(m, a96ff, residueMap, sysStrNum, ligStrNum, backrubFile, hSteric, overlapThresh);
					}
				}
			}

			boolean done = false;
			int curRot = 0;
			int totRotForCur = grl.getNumRotForAAtype(curLigAANum);
			
			while ( (!done) && ((curRot<totRotForCur) || (totRotForCur==0)) ){
				
				if (totRotForCur!=0) // loop through rotamers
					ligROT.applyRotamer(m, 0, curRot);
				else // no rotamers, so only the current AA state has to be computed and we exit the while loop
					done = true;
						
				computeIntraRotEnergiesHelper(-1, -1, curRot, totRotForCur,	retEMatrixMin, 
						retEMatrixMax, rotamerIndexOffset, residueMap, numTotalRotamers, numInAS, true, minimizeBB, doBackrubs);
				
				curRot++;
			}
		}*/
		
		// Restore the energy eval and flexibility flags
		/*for(int i=0;i<m.numberOfResidues;i++){
			m.residue[i].setEnergyEval(savedEnergyEvalSC[i], savedEnergyEvalBB[i]);
			m.residue[i].flexible = savedFlexible[i];
		}*/
	}
	//Computes the intra-residue min and max energies for a given rotamer
	private void computeIntraRotEnergiesHelper(int curPos/*mutation residue specific*/, int strResNum, int curAA, int curRot, int curStr, int totRotForCur, 
			PairwiseEnergyMatrix retEMatrixMin, PairwiseEnergyMatrix retEMatrixMax,
			int strandMut[][], boolean minimizeBB, boolean doBackrubs){
		
		float curEnergy = 0.0f;	
		float beginE = 0.0f;
		float minEnergy = (float)Math.pow(10,30);
		float maxEnergy = -(float)Math.pow(10,30);
		
		//Compute min and max energies
		if ( doMinimization ){								
			
			//Iniitialize Amber for the current residue
			a96ff.initializeCalculation();
			a96ff.setNBEval(hElect,hVDW);
			
			if ( (!minimizeBB) && (!doBackrubs) ) { //side-chain minimization
				
				//first, compute the initial energy at the new position
				beginE = calcTotalSnapshotEnergy();
								
				//minimize, starting at the initial position
				simpMin.minimize(numMinSteps);
				curEnergy = calcTotalSnapshotEnergy();
				
				
				//Compare to the min and max energies found so far and update, if necessary
				if (beginE<curEnergy)
					curEnergy = beginE;
				float lE = Math.min(beginE, curEnergy);
				float hE = Math.max(beginE, curEnergy);
				minEnergy = Math.min(minEnergy,lE);									
				maxEnergy = Math.max(maxEnergy,hE);

				m.updateCoordinates();//restore the actualCoordinates array to the initial values
                                if(doPerturbations)
                                    m.revertPertParamsToCurState();
			}
			else if (!doBackrubs) { //phi/psi backbone minimization, so only rotate the backbone O and HN		
				
				int at[] = new int[5];
				int numAtoms = 0;
				Residue r1 = null;
				r1 = m.strand[curStr].residue[strResNum];
				
				numAtoms = r1.numberOfAtoms;
				
				//get the atoms
				for (int i=0; i<numAtoms; i++){
					if (r1.atom[i].name.equalsIgnoreCase("CA"))
						at[0] = r1.atom[i].moleculeAtomNumber;
					else if (r1.atom[i].name.equalsIgnoreCase("C"))
						at[1] = r1.atom[i].moleculeAtomNumber;
					else if (r1.atom[i].name.equalsIgnoreCase("O"))
						at[2] = r1.atom[i].moleculeAtomNumber;
					else if (r1.atom[i].name.equalsIgnoreCase("N"))
						at[3] = r1.atom[i].moleculeAtomNumber;
					else if (r1.atom[i].name.equalsIgnoreCase("H"))
						at[4] = r1.atom[i].moleculeAtomNumber;
				}
				
				//get the CA-C bond
				float dx = m.actualCoordinates[at[1]*3] - m.actualCoordinates[at[0]*3];
				float dy = m.actualCoordinates[at[1]*3+1] - m.actualCoordinates[at[0]*3+1];
				float dz = m.actualCoordinates[at[1]*3+2] - m.actualCoordinates[at[0]*3+2];
				
				//get the N-CA bond
				float dxH = m.actualCoordinates[at[0]*3] - m.actualCoordinates[at[3]*3];
				float dyH = m.actualCoordinates[at[0]*3+1] - m.actualCoordinates[at[3]*3+1];
				float dzH = m.actualCoordinates[at[0]*3+2] - m.actualCoordinates[at[3]*3+2];
				
				//get the center of rotation for O (the actualCoordinates[] of C)
				double center[] = new double[3];
				center[0] = m.actualCoordinates[at[1]*3];
				center[1] = m.actualCoordinates[at[1]*3+1];
				center[2] = m.actualCoordinates[at[1]*3+2];
				
				//get the center of rotation for H (the actualCoordinates[] of N)
				double centerH[] = new double[3];
				centerH[0] = m.actualCoordinates[at[3]*3];
				centerH[1] = m.actualCoordinates[at[3]*3+1];
				centerH[2] = m.actualCoordinates[at[3]*3+2];
				
				float rotForInitPos = bbMin.getMaxDihedRot(); //get the max phi/psi rotation
				
				//Do the sampling and minimization
				for (int curSample=0; curSample<samplesBB; curSample++){
				
					//randomly generate the rotation angle
					float rotChange[] = new float[2];
					Random randNum = new Random();
					
					if (curSample!=0) {
						for (int i=0; i<2; i++)
							rotChange[i] = (randNum.nextFloat()-0.5f)*rotForInitPos*2.0f;
					}
					else {
						for (int i=0; i<2; i++)
							rotChange[i] = 0.0f;
					}
					
					//Compute the energy corresponding to the new positions
					if (curSample!=0){ //do not apply a change for the initial position
						m.rotateAtom(at[2], dx, dy, dz, center[0], center[1], center[2], rotChange[0], false);
						m.rotateAtom(at[4], dxH, dyH, dzH, centerH[0], centerH[1], centerH[2], rotChange[1], false);
					}
					
					//compute the initial energy at the new position
					curEnergy = calcTotalSnapshotEnergy();
					
					//Compare to the min and max energies found so far and update, if necessary;
					//For intra-energies with BB flexibility, the initial point is taken as the max energy,
					//		since the O and H positions are only sampled without minimization
					minEnergy = Math.min(minEnergy,curEnergy);
					if (curSample==0)
						maxEnergy = Math.max(maxEnergy,curEnergy);
	
					m.updateCoordinates();//restore the actualCoordinates array to the initial values
				}
			}
			else { //backrub minimization
				if (curPos>=0){ //not the ligand
					beginE = calcTotalSnapshotEnergy();
					float e[] = brMin.getMinMaxIntraEnergyBR(curPos);
					float lE = Math.min(beginE, e[0]);
					float hE = Math.min(beginE, e[1]);
					minEnergy = Math.min(minEnergy,lE);									
					maxEnergy = Math.max(maxEnergy,hE);
				}
				else { //the ligand (backrubs are not applied and have no effect on the intra-energy)
					minEnergy = calcTotalSnapshotEnergy();
					maxEnergy = minEnergy;
				}
                                    m.updateCoordinates();//restore the actualCoordinates array to the initial values
			}
		}
		else if (computeEVEnergy){
			a96ff.initializeCalculation();
			a96ff.setNBEval(hElect,hVDW);
			curEnergy = calcTotalSnapshotEnergy();
			m.updateCoordinates();
			minEnergy = curEnergy;
			maxEnergy = curEnergy;
		}

		// Store result
		int pos = -1;
		int aa = -1;
		//if (!isLig){ //AS residue
			pos = curPos;
			aa = curAA;
		//}
		/*else { //ligand
			pos = numInAS;
			aa = ligROT.getIndexOfNthAllowable(0,0);
		}*/
		
		retEMatrixMin.setIntraE(pos, aa, curRot, minEnergy);
                retEMatrixMax.setIntraE(pos, aa, curRot, maxEnergy);
	}
	// Sets up the computation for residue-to-template (SHL-AS, LIG-SHL) and rotamer-rotamer (AS-AS, LIG-AS) energies;
	// Searches among amino acid types for all mutatable residues (as determined by residueMutatable);
	//		Sets all non-mutatable residues in residueMap to either Gly or Ala (and leaves all Pro)
	public void pairwiseEnergyComp (final int maxDepth, int strandMut[][], int residueMutatable[],
		PairwiseEnergyMatrix retEMatrixMin, PairwiseEnergyMatrix retEMatrixMax, boolean shellRun, boolean minimizeBB, boolean doBackrubs, String backrubFile) {
		
		
		/*if (ligPresent) { //there is a ligand and its energies should be computed
			if(ligROT.getNumAllowable(0)==0) {  // this shouldn't happen
				System.out.println("ERROR: Ligand has no allowables but you are using a ligand?");
				System.exit(1);
			}
			else {
				// Because this is a ligand there can be only one allowable type, ie. it can not mutate;
				// Change the residue type to this one type (if protein)
				curLigAANum = ligROT.getIndexOfNthAllowable(0,0);
				if (m.strand[ligStrNum].isProtein)		
					ligROT.changeResidueType(m,0,grl.getAAName(curLigAANum),addHydrogens);
				m.strand[ligStrNum].residue[0].flexible = true;
				m.strand[ligStrNum].residue[0].setEnergyEval(true, true);
			}
		}*/
		
		//Check with all allowed amino acids for all of the mutatable positions in residueMutatable[];
		//	all non-mutatable positions in residueMutatable (these are necessarily in residueMap[]) are set ot Gly or Ala (or remain Pro)
		
		int numMut = 0; //should be either 1 (shellRun) or 2 (rotamer-rotamer energies)
		int mutDepth[] = new int[2];
		for (int i=0; i<mutDepth.length; i++)
			mutDepth[i] = -1;
		
		for (int depth=0; depth<maxDepth; depth++){	
			if (residueMutatable[depth]!=0){
				mutDepth[numMut] = depth;
				numMut++;
			}
		}
		
		if ( ( ((shellRun) && (numMut!=1)) || ((!shellRun) && (numMut!=2)) ) ) {
			System.out.println("ERROR: incorrect number of mutatable positions in residueMutatable[]");
			System.exit(1);
		}
		
		//KER: this needs to be called to initialize the cterm and nterm flags before
		//KER: anything is mutated.
		a96ff.calculateTypesWithTemplates();
			
		for (int depth=0; depth<maxDepth; depth++){				
			if (residueMutatable[depth]==0){ //not a mutatable residue, but in residueMap[]
				int str = mutRes2Strand[depth];
				int strMutIndex = mutRes2StrandMutIndex[depth];
				// Make this residue a "gly", if it is not already GLY or PRO;
				//		If doing backbone minimization, then change to ALA, if not already GLY or PRO
				String a = strandRot[str].getCurRotType(strandMut[str][strMutIndex]);
				if ( (!a.equalsIgnoreCase("GLY")) &&
                                        ( (!a.equalsIgnoreCase("PRO")) || doPerturbations ) ){
                                        //DEEPer can mutate prolines so it treats them like other residues
					String tmpAA = "gly";
					if (minimizeBB || doPerturbations)
						tmpAA = "ala";
				
					//TODO: KER: strandRot should have the changeResidueType should be able to mutate things other than proteins
					if(m.strand[str].isProtein)
						strandRot[str].changeResidueType(m,strandMut[str][strMutIndex],tmpAA,addHydrogens,connectResidues);
				}
				int molResNumber = m.strand[str].residue[strandMut[str][strMutIndex]].moleculeResidueNumber;
				curAANum[molResNumber] = -1;
				m.strand[str].residue[strandMut[str][strMutIndex]].flexible = false;
				if(isTemplateOn&&(!doPerturbations))
					m.strand[str].residue[strandMut[str][strMutIndex]].setEnergyEval(false, true);
				else
					m.strand[str].residue[strandMut[str][strMutIndex]].setEnergyEval(false, false);
			}
		}

		pairwiseEnergyCompAllMutatedResHelper(maxDepth,	strandMut,	-1, -1, residueMutatable, retEMatrixMin, retEMatrixMax, shellRun, 
				minimizeBB, doBackrubs, numMut, mutDepth, 0, backrubFile);
	}
	// Helper to pairwiseEnergyComp();
	// Searches among amino acid types for all mutatable residues (as determined by residueMutatable[]);
	//		the non-mutatable residues in residueMutatable[] have already been set to either Gly or Ala (all Pro remain)
	private void pairwiseEnergyCompAllMutatedResHelper(int maxDepth, int strandMut[][], int res1Num, int res2Num, int residueMutatable[],
			PairwiseEnergyMatrix retEMatrixMin, PairwiseEnergyMatrix retEMatrixMax, boolean shellRun, boolean minimizeBB, boolean doBackrubs, final int numMut,
			int mutDepth[], int curMut, String backrubFile){
		
		if (curMut>=numMut){// If we've arrived here, then we have assigned a full amino acid sequence and we are ready to generate all rotamer combinations
			
			if(debug){
				System.out.print(".");
			}
			if (computeEVEnergy){
				a96ff.calculateTypesWithTemplates();
				a96ff.initializeCalculation();
				a96ff.setNBEval(hElect,hVDW);

                //Residue got mutated so molecule indices changed and we need to reset the saved AmberFF
                //KER: if we had a molecule that didn't have indices changed we could store AmberFF terms
                templateA96ff = null;
                res1A96ff = null;
                res2A96ff = null;
                
                setTemplateOff(strandMut,mutDepth[0],mutDepth[1],shellRun,false);
                
				if (doMinimization){
					if (!minimizeBB){ //side-chain minimization
						simpMin.initialize(m,numberOfStrands,a96ff,strandRot,curAANum,doDihedE);
					}
					else { //backbone minimization
						if (!doBackrubs){
							bbMin.initialize(m, a96ff, strandMut, numberOfStrands);
						}
						else {
							brMin.initialize(m, a96ff, strandMut, backrubFile, hSteric, overlapThresh, numberOfStrands,true);
						}
					}
				}
			}
			int res1 = -1, res2 = -1;
			for(int i=0;i<maxDepth;i++){
				if (residueMutatable[i]==1) {
					if (res1 == -1)
						res1 = i;
					else
						res2 = i;
				}
			}
			
			/*if (ligPresent){
				// The first residue is the ligand (so -2 is passed), the second is numerically res1
				pairwiseMutationAllRotamerSearch(maxDepth, residueMap, -2, res1, -2,
					res1Num,-1,-1,rotamerIndexOffset,numTotalRotamers,
					ligPresent, retEMatrixMin, retEMatrixMax, minimizeBB, doBackrubs, numMut, mutDepth, 0);
			}
			else{*/
				pairwiseMutationAllRotamerSearch(maxDepth, strandMut, res1, res2, res1Num,
					res2Num,-1,-1,retEMatrixMin, retEMatrixMax, minimizeBB, doBackrubs, numMut, mutDepth, 0, true);
			//}
			
			return;
		}
		else {
			boolean isFirstAA = (res1Num == -1);
			int depth = mutDepth[curMut];
			int str   = mutRes2Strand[depth];
			int strResNum = strandMut[str][mutRes2StrandMutIndex[depth]];
			for(int q=0;q<strandRot[str].getNumAllowable(strResNum);q++) {
				
				int AAindex = strandRot[str].getIndexOfNthAllowable(strResNum,q);
				// Change the residue type

				if (isFirstAA){
					res1Num = AAindex;
				}
				else{
					res2Num = AAindex;
				}
				curAANum[m.strand[str].residue[strResNum].moleculeResidueNumber] = AAindex;


				// Apply mutation
				if(m.strand[str].isProtein)
					strandRot[str].changeResidueType(m,strResNum,strandRot[str].rl.getAAName(curAANum[m.strand[str].residue[strResNum].moleculeResidueNumber]),addHydrogens,connectResidues);


				m.strand[str].residue[strResNum].flexible = true;
				
				m.strand[str].residue[strResNum].setEnergyEval(true,true);
				pairwiseEnergyCompAllMutatedResHelper(maxDepth,	strandMut,res1Num,res2Num, residueMutatable, 
					retEMatrixMin, retEMatrixMax, shellRun, minimizeBB, doBackrubs, 
					numMut, mutDepth, curMut+1, backrubFile);
			}
		}
	}
	// Called by pairwiseEnergyCompAllMutatedResHelper();
	// Computes the energies among all rotamer combinations for residues res1 and res2 (if res2==1, then computes the res1-to-template energies)
	private void pairwiseMutationAllRotamerSearch(int maxDepth, int strandMut[][],
		int res1, int res2, int res1AANum, int res2AANum, int res1RotNum, int res2RotNum,
		PairwiseEnergyMatrix retEMatrixMin, PairwiseEnergyMatrix retEMatrixMax,
		boolean minimizeBB, boolean doBackrubs, final int numMut, int mutDepth[], int curMut, boolean validConf) {
//validConf indicates that all RC applications were successful (always true if not DEEPer)

		// If we're at the ligand depth
		/*if (isLigRun) {
			isLigRun = false;
			int numLigRot = grl.getNumRotForAAtype(curLigAANum);
			if (numLigRot==0) { //ligand is ALA or GLY
				pairwiseMutationAllRotamerSearch(false,maxDepth,residueMap,res1,res2,
						res1AANum, res2AANum, 0, res2RotNum, rotamerIndexOffset,
						numTotalRotamers, ligPresent, retEMatrixMin, retEMatrixMax, 
						minimizeBB, doBackrubs, numMut, mutDepth, curMut);
			}
			else {
				for(int w=0;w<grl.getNumRotForAAtype(curLigAANum);w++){
					ligROT.applyRotamer(m, 0, w);
					pairwiseMutationAllRotamerSearch(false,maxDepth,residueMap,res1,res2,
						res1AANum, res2AANum, w, res2RotNum, rotamerIndexOffset,
						numTotalRotamers, ligPresent, retEMatrixMin, retEMatrixMax, 
						minimizeBB, doBackrubs, numMut, mutDepth, curMut);
				}
			}
		}*/
		//else {
			if (curMut>=numMut) {
				
				int res1Strand = -1;
				int res1Resnum = -1;
				int res2Strand = -1;
				int res2Resnum = -1;
				int res1MutIndex = -1;
				int res2MutIndex = -1;
				
				boolean shellRun = false;
				if (res2 == -1)
					shellRun = true;
				
				//int pos1 = -1;
				//int aa1 = -1;
				
				res1Strand = mutRes2Strand[res1];
				res1Resnum = strandMut[res1Strand][mutRes2StrandMutIndex[res1]];
				res1MutIndex = mutRes2StrandMutIndex[res1];
				if(res2!=-1){
					res2Strand = mutRes2Strand[res2];
					res2Resnum = strandMut[res2Strand][mutRes2StrandMutIndex[res2]];
					res2MutIndex = mutRes2StrandMutIndex[res2];
				}


                                //First, if using DEEPer and calculating a pairwise energy,
                                //check if the pair is parametrically incompatible
                                //and set the pairwise energy to infinity if it is
                                //Do the same for residues with non-pre-Pro compatible BB dihedrals
                                //right before a proline
                                //and for RCs whose application failed
                                if(doPerturbations){

                                    Residue firstRes = m.strand[res1Strand].residue[res1Resnum], secondRes;

                                    //For a shell run, if the next residue is a shell proline, make sure the RC's backbone dihedrals are OK for pre-pro
                                    //Also, if this residue is a proline and the previous residue is in the shell, make sure that shell residue's backbone dihedrals are OK for pre-pro
                                    if( shellRun ){

                                        if(!firstRes.validConf)//Pro ring not closed
                                            validConf = false;

                                        if( m.checkCBonded(firstRes.moleculeResidueNumber) ){

                                            Residue nextRes = m.strand[res1Strand].residue[res1Resnum+1];
                                            if( nextRes.getEnergyEvalSC() && nextRes.name.equalsIgnoreCase("PRO") ){//We are only calculating nextRes' energy if it's in the shell
                                                if( ! RamachandranChecker.getInstance().checkPrePro(m, firstRes.moleculeResidueNumber) ){
                                                    validConf = false;
                                                }
                                            }
                                        }


                                        if( m.checkNBonded(firstRes.moleculeResidueNumber) && firstRes.name.equalsIgnoreCase("PRO") ){

                                            Residue prevRes = m.strand[res1Strand].residue[res1Resnum-1];
                                            if( prevRes.getEnergyEvalSC() ){//We are only calculating prevRes' energy if it's in the shell
                                                if( ! RamachandranChecker.getInstance().checkPrePro(m, prevRes.moleculeResidueNumber) ){
                                                    validConf = false;
                                                }
                                            }
                                        }

                                        if( !validConf ){
                                            retEMatrixMin.setShellRotE( res1, res1AANum, res1RotNum, Float.POSITIVE_INFINITY );
                                            retEMatrixMax.setShellRotE( res1, res1AANum, res1RotNum, Float.POSITIVE_INFINITY );
                                            return;
                                        }
                                    }
                                    else if (!shellRun){

                                        boolean incompatible = !validConf;

                                        secondRes = m.strand[res2Strand].residue[res2Resnum];

                                        if( res1Strand == res2Strand ){
                                            if( (res2Resnum == res1Resnum + 1) && ( secondRes.name.equalsIgnoreCase("PRO") ) ){//firstRes is right before secondRes, which is proline
                                                if( ! RamachandranChecker.getInstance().checkPrePro(m, firstRes.moleculeResidueNumber) )
                                                    incompatible = true;//the pair of RCs is incompatible, because secondRes has AA type proline and firstRes' backbone dihedrals are not OK for pre-pro
                                            }
                                            else if( (res2Resnum == res1Resnum - 1) && ( firstRes.name.equalsIgnoreCase("PRO") ) ){//secondRes if right before firstRes, which is proline
                                                if( ! RamachandranChecker.getInstance().checkPrePro(m, secondRes.moleculeResidueNumber) )
                                                    incompatible = true;
                                            }
                                        }

                                        if( ! (firstRes.validConf&&secondRes.validConf) )
                                            incompatible = true;

                                        if( isParametricallyIncompatible(firstRes,res1AANum,res1RotNum,secondRes,res2AANum,res2RotNum) )
                                            incompatible = true;


                                        if(incompatible){
                                            retEMatrixMin.setPairwiseE(res1, res1AANum, res1RotNum, res2, res2AANum, res2RotNum, Float.POSITIVE_INFINITY);
                                            retEMatrixMax.setPairwiseE(res1, res1AANum, res1RotNum, res2, res2AANum, res2RotNum, Float.POSITIVE_INFINITY);
                                            return;
                                        }
                                    }



                                }


				//If we've gotten here then check the sterics to make sure we're sterically allowable;
				//Exclude the other residues from residueMap[], since they have been mutated to Gly/Ala
				int excludeRes[][] = new int[strandMut.length][];
				
				
				for(int str=0;str<numberOfStrands;str++){
					excludeRes[str] = new int[strandMut[str].length];
					System.arraycopy(strandMut[str], 0, excludeRes[str], 0, strandMut[str].length);
				}
				excludeRes[res1Strand][res1MutIndex] = -1; //do not exclude res1, since it is fixed here
				if (res2!=-1)
					excludeRes[res2Strand][res2MutIndex] = -1; //do not exclude res2, since it is fixed here
				boolean stericallyGood = false;
				if (RS_CheckAllSterics(res1Strand,res1Resnum,excludeRes)) {
					if (res2 != -1) {
						if (RS_CheckAllSterics(res2Strand,res2Resnum,excludeRes))
							stericallyGood = true;
					}
					else 
						stericallyGood = true;
				}
				
				if ((stericallyGood)) { //good steric found or doing backbone minimization
					
					// After minimization do a m.updateCoordinates() to resync the actualCoordinates which were changed
					//  in the minimization procedure									
					
					if (doMinimization){
						if (!minimizeBB){ //side-chain minimization
							pairwiseRotamerEnergySidechainHelper (res1,res2,res1Strand,res1Resnum,res2Strand,res2Resnum,strandMut,
									shellRun, retEMatrixMin, retEMatrixMax, res1AANum, res1RotNum, res2AANum, res2RotNum);
						}
						else { //backbone minimization
							if (!doBackrubs) // phi/psi minimization
								pairwiseRotamerEnergyBackboneHelper (res1,res2,res1Strand,res1Resnum,res2Strand,res2Resnum,strandMut,
										shellRun, retEMatrixMin, retEMatrixMax, false, res1AANum, res1RotNum, res2AANum, res2RotNum);
							else { //backrubs
								pairwiseRotamerEnergyBackrubsHelper (res1,res2,res1Strand,res1Resnum,res2Strand,res2Resnum,strandMut,
										shellRun, retEMatrixMin, retEMatrixMax, false, res1AANum, res1RotNum, res2AANum, res2RotNum);
							}
						}
					}
					else if (computeEVEnergy){
						//KER: Only need to initialize when mutating or changing eval flags
						/*a96ff.initializeCalculation();
						a96ff.setNBEval(hElect,hVDW);*/
						double energyTerms[] = templateOffA96ff.calculateTotalEnergy(m.actualCoordinates,-1); //compute the energy
						float curEnergy = (float)energyTerms[0];//calcTotalSnapshotEnergy();
						
						if(!shellRun){
							//Remove the intra-residue energies here
							float tmpE = computeEnergyOfOnlyRes(m.strand[res1Strand].residue[res1Resnum].moleculeResidueNumber,0,strandMut,res1A96ff);
							curEnergy -= tmpE;
							
							tmpE = computeEnergyOfOnlyRes(m.strand[res2Strand].residue[res2Resnum].moleculeResidueNumber,1,strandMut,res2A96ff);
							curEnergy -= tmpE;
						
						}
						else {
							// If this is a shell run then subtract off the shell-to-shell
							//  interaction energies so that they are not over counted
							curEnergy -= computeEnergyOfOnlyTemplate(strandMut);
						}
						//Store the computed energy
						if (shellRun) {
							retEMatrixMin.setShellRotE(res1, res1AANum, res1RotNum, curEnergy);
                                                        retEMatrixMax.setShellRotE(res1, res1AANum, res1RotNum, curEnergy);
						}
						else {
							retEMatrixMin.setPairwiseE( res1, res1AANum, res1RotNum, res2, res2AANum, res2RotNum, curEnergy );
                                                        retEMatrixMax.setPairwiseE( res1, res1AANum, res1RotNum, res2, res2AANum, res2RotNum, curEnergy );
						}
					}
					else {
						System.out.println("This should not happen. No energy evaluation specified");
						System.exit(1);
					}					
															
					m.updateCoordinates();
                                        if(doPerturbations)
                                                m.revertPertParamsToCurState();
				}
				else {
					if (shellRun) {
						retEMatrixMin.setShellRotE( res1, res1AANum, res1RotNum, stericE );
                                                retEMatrixMax.setShellRotE( res1, res1AANum, res1RotNum, stericE );
					}
					else {
                                                retEMatrixMin.setPairwiseE( res1, res1AANum, res1RotNum, res2, res2AANum, res2RotNum, stericE );
                                                retEMatrixMax.setPairwiseE( res1, res1AANum, res1RotNum, res2, res2AANum, res2RotNum, stericE );
					}
				}

				return;
			}

			else {
				// Next check with different rotamers
				int depth = mutDepth[curMut];
				int str = mutRes2Strand[depth];
                                int strResNum = strandMut[str][mutRes2StrandMutIndex[depth]];
				int molResNum = m.strand[str].residue[strResNum].moleculeResidueNumber;
				
                                int numRot;
                                if(doPerturbations)
                                    numRot = ((StrandRCs)strandRot[str]).getNumRCs(strResNum, curAANum[molResNum]);
                                else
                                    numRot = strandRot[str].rl.getNumRotForAAtype(curAANum[molResNum]);
				
				if (numRot == 0) {//If there are no rotamers for this AA then allow the default conformation; this will only happen with Ala and Gly
					if (depth==res1)
						pairwiseMutationAllRotamerSearch(maxDepth,strandMut,res1,res2,
							res1AANum, res2AANum, 0, res2RotNum, retEMatrixMin, retEMatrixMax, minimizeBB, doBackrubs,
							numMut, mutDepth, curMut+1, validConf);
					else
						pairwiseMutationAllRotamerSearch(maxDepth,strandMut,res1,res2,
							res1AANum, res2AANum, res1RotNum, 0, retEMatrixMin, retEMatrixMax, minimizeBB, doBackrubs,
							numMut, mutDepth, curMut+1, validConf);
				}
				else {		
					for(int w=0;w<numRot;w++) {


                                                if(doPerturbations)
                                                    validConf = validConf && ((StrandRCs)strandRot[str]).applyRC(m, strandMut[str][mutRes2StrandMutIndex[depth]], w);
                                                else
                                                    strandRot[str].applyRotamer(m, strandMut[str][mutRes2StrandMutIndex[depth]], w);

						if (depth==res1)
							pairwiseMutationAllRotamerSearch(maxDepth,strandMut,res1,res2,
								res1AANum, res2AANum, w, res2RotNum, retEMatrixMin, retEMatrixMax, minimizeBB, doBackrubs,
								numMut, mutDepth, curMut+1, validConf);
						else
							pairwiseMutationAllRotamerSearch(maxDepth,strandMut,res1,res2,
								res1AANum, res2AANum, res1RotNum, w, retEMatrixMin, retEMatrixMax, minimizeBB, doBackrubs,
								numMut, mutDepth, curMut+1, validConf);
					}
				}
			}
		}	
	//This method computes the minimized energy for a given pair of rotamers (or rot-to-template)
	//Called by pairwiseMutationAllRotamerSearch(.)
	private void pairwiseRotamerEnergySidechainHelper (int res1,int res2,int res1Strand,int res1Resnum,int res2Strand,int res2Resnum,int strandMut[][],
			boolean shellRun, PairwiseEnergyMatrix retEMatrixMin, PairwiseEnergyMatrix retEMatrixMax,
			int res1AANum, int res1RotNum, int res2AANum, int res2RotNum){
		
		//Initialize Amber for the current pair
		//KER: only need to initialize when mutating or changing eval flags
		//a96ff.initializeCalculation();
		//a96ff.setNBEval(hElect,hVDW);
		
		float minEnergy = (float)Math.pow(10,30);
		float maxEnergy = -(float)Math.pow(10,30);
		
		float curEnergy = 0.0f;
		float beginE = 0.0f;
		
		
		/////////////////////////////
		//formally making sure that when minimization is performed,
		//		the minimized value is never greater than the initial value
		double energyTerms[] = templateOffA96ff.calculateTotalEnergy(m.actualCoordinates,-1); //compute the energy
		beginE = (float)energyTerms[0];//calcTotalSnapshotEnergy();
		
		//KER: Include intra res with shell energy so we don't subtract off the first residue energy for the
		//KER: shell run
		if (!shellRun) {
			//Remove the intra-residue energies here
			float tmpE1 = computeEnergyOfOnlyRes(m.strand[res1Strand].residue[res1Resnum].moleculeResidueNumber,0,strandMut,res1A96ff);
			beginE -= tmpE1;
			
			tmpE1 = computeEnergyOfOnlyRes(m.strand[res2Strand].residue[res2Resnum].moleculeResidueNumber,1,strandMut,res2A96ff);
			beginE -= tmpE1;
		}
		else {
			// If this is a shell run then subtract off the shell-to-shell
			//  interaction energies so that they are not over counted
			beginE -= computeEnergyOfOnlyTemplate(strandMut);
		}
		/////////////////////////////

		//Minimize
		simpMin.minimize(numMinSteps); 
		
		energyTerms = templateOffA96ff.calculateTotalEnergy(m.actualCoordinates,-1); //compute the energy
		curEnergy = (float)energyTerms[0];//calcTotalSnapshotEnergy();
		
		if (!shellRun) {
			//Remove the intra-residue energies here
			float tmpE2 = computeEnergyOfOnlyRes(m.strand[res1Strand].residue[res1Resnum].moleculeResidueNumber,0,strandMut,res1A96ff);
			curEnergy -= tmpE2;
			
			tmpE2 = computeEnergyOfOnlyRes(m.strand[res2Strand].residue[res2Resnum].moleculeResidueNumber,1,strandMut,res2A96ff);
			curEnergy -= tmpE2;
		}
		else {
			// If this is a shell run then subtract off the shell-to-shell
			//  interaction energies so that they are not over counted
			curEnergy -= computeEnergyOfOnlyTemplate(strandMut);
			
			if (doDihedE) //add dihedral energies
				curEnergy += simpMin.computeDihedEnergy();
		}
		/////////////////////////////

		//Compare to the min and max energies found so far and update, if necessary
		if(doDihedE){
			if (beginE<curEnergy)
				curEnergy = beginE;
			float lE = Math.min(beginE, curEnergy);
			float hE = Math.max(beginE, curEnergy);
		}
		minEnergy = Math.min(minEnergy,curEnergy);
		maxEnergy = Math.max(maxEnergy,beginE);		


                m.updateCoordinates();//restore the actualCoordinates array to the initial values
                if(doPerturbations)
                    m.revertPertParamsToCurState();

		if (shellRun) {
			retEMatrixMin.setShellRotE(res1, res1AANum, res1RotNum, minEnergy);
                        retEMatrixMax.setShellRotE(res1, res1AANum, res1RotNum, maxEnergy);
		}
		else {
			retEMatrixMin.setPairwiseE( res1, res1AANum, res1RotNum, res2, res2AANum, res2RotNum, minEnergy );
                        retEMatrixMax.setPairwiseE( res1, res1AANum, res1RotNum, res2, res2AANum, res2RotNum, maxEnergy );
		}
		
		return;
	}
	//This method helps compute energy bounds for a given pair of rotamers for backbone phi/psi minimization
	//Called by pairwiseMutationRotamerSearch(.)
	private void pairwiseRotamerEnergyBackboneHelper (int res1,int res2,int res1Strand,int res1Resnum,int res2Strand,int res2Resnum,int strandMut[][],
			boolean shellRun, PairwiseEnergyMatrix retEMatrixMin, PairwiseEnergyMatrix retEMatrixMax, boolean templateOnly,
			int res1AANum, int res1RotNum, int res2AANum, int res2RotNum){		

		//Initialize Amber for the current pair
		a96ff.initializeCalculation();
		a96ff.setNBEval(hElect,hVDW);
	
		float minEnergy = (float)Math.pow(10,30);
		float maxEnergy = -(float)Math.pow(10,30);
			
		float beginE = 0.0f;
		
		double energyTerms[] = templateOffA96ff.calculateTotalEnergy(m.actualCoordinates,-1); //compute the energy
		beginE = (float)energyTerms[0];//calcTotalSnapshotEnergy();
		beginE -= removeExtraE(templateOnly, res1Strand, res1Resnum, res2Strand, res2Resnum, strandMut);
		
		bbMin.minimizeFull(true); //minimize
		
		energyTerms = templateOffA96ff.calculateTotalEnergy(m.actualCoordinates,-1); //compute the energy
		float curEnergy = (float)energyTerms[0];//calcTotalSnapshotEnergy();
		curEnergy -= removeExtraE(templateOnly, res1Strand, res1Resnum, res2Strand, res2Resnum, strandMut);
		
		
		if (beginE<curEnergy)
			curEnergy = beginE;
		float lE = Math.min(beginE, curEnergy);
		float hE = Math.max(beginE, curEnergy);
		minEnergy = Math.min(minEnergy,lE);
		maxEnergy = Math.max(maxEnergy,hE);	
								
		m.updateCoordinates();//restore the actualCoordinates array to the initial values
		/////////////////////////////////////////////
		
		if (templateOnly){
			retEMatrixMin.setShellShellE(minEnergy);
                        retEMatrixMax.setShellShellE(maxEnergy);
		}
		else if (shellRun) {
			retEMatrixMin.setShellRotE(res1, res1AANum, res1RotNum, minEnergy);
                        retEMatrixMax.setShellRotE(res1, res1AANum, res1RotNum, maxEnergy);
		}
		else {
			retEMatrixMin.setPairwiseE( res1, res1AANum, res1RotNum, res2, res2AANum, res2RotNum, minEnergy );
                        retEMatrixMax.setPairwiseE( res1, res1AANum, res1RotNum, res2, res2AANum, res2RotNum, maxEnergy );
		}
		
		return;
	}
	//This method helps perform sampling for a given pair of rotamers for backrubs minimization
	//Called by pairwiseMutationRotamerSearch(.)
	private void pairwiseRotamerEnergyBackrubsHelper (int res1,int res2,int res1Strand,int res1Resnum,int res2Strand,int res2Resnum,
			int strandMut[][], boolean shellRun, PairwiseEnergyMatrix retEMatrixMin, PairwiseEnergyMatrix retEMatrixMax, boolean templateOnly,
			int res1AANum, int res1RotNum, int res2AANum, int res2RotNum){
		
		//Initialize Amber for the current pair
		//KER: only need to initialize when mutating or changing eval flags
		//a96ff.initializeCalculation();
		//a96ff.setNBEval(hElect,hVDW);
	
		float minEnergy = (float)Math.pow(10,30);
		float maxEnergy = -(float)Math.pow(10,30);
			
		double energyTerms[] = templateOffA96ff.calculateTotalEnergy(m.actualCoordinates,-1); //compute the energy
		float beginE = (float)energyTerms[0];//calcTotalSnapshotEnergy();
		
		beginE -= removeExtraE(templateOnly, res1Strand, res1Resnum, res2Strand, res2Resnum, strandMut);
		
		float initActualCoords[] = m.getActualCoords(); //this is necessary for performing brMin.applyMaxBackrub() after brMin.minimizeFull()
		brMin.minimizeFull(shellRun,templateOnly); //minimize
		
		energyTerms = templateOffA96ff.calculateTotalEnergy(m.actualCoordinates,-1); //compute the energy
		float curEnergy = (float)energyTerms[0];//calcTotalSnapshotEnergy();
		curEnergy -= removeExtraE(templateOnly, res1Strand, res1Resnum, res2Strand, res2Resnum, strandMut);
		m.setActualCoords(initActualCoords);
		
		brMin.applyMaxBackrub();
		energyTerms = templateOffA96ff.calculateTotalEnergy(m.actualCoordinates,-1); //compute the energy
		float maxE = (float)energyTerms[0];//calcTotalSnapshotEnergy();
		maxE -= removeExtraE(templateOnly, res1Strand, res1Resnum, res2Strand, res2Resnum, strandMut);
		beginE = (float)Math.min(beginE,maxE); //max energy will be the min of (initial energy) and (max sterically-allowed energy from brMin)
		
		if (beginE<curEnergy)
			curEnergy = beginE;
		float lE = Math.min(beginE, curEnergy);
		float hE = Math.max(beginE, curEnergy);
		minEnergy = Math.min(minEnergy,lE);
		maxEnergy = Math.max(maxEnergy,hE);	
								
		m.updateCoordinates();//restore the actualCoordinates array to the initial values from the atoms coord[]
		/////////////////////////////////////////////
		
		if (templateOnly){
			retEMatrixMin.setShellShellE(minEnergy);
                        retEMatrixMax.setShellShellE(maxEnergy);
		}
		else if (shellRun) {
			retEMatrixMin.setShellRotE(res1, res1AANum, res1RotNum, minEnergy);
                        retEMatrixMax.setShellRotE(res1, res1AANum, res1RotNum, maxEnergy);
		}
		else {
			retEMatrixMin.setPairwiseE( res1, res1AANum, res1RotNum, res2, res2AANum, res2RotNum, minEnergy );
                        retEMatrixMax.setPairwiseE( res1, res1AANum, res1RotNum, res2, res2AANum, res2RotNum, maxEnergy );
		}
		
		return;
	}
	
	//Computes the intra- and shell-based energies that have been double-counted in the pairwise
	//		energy computation and must be removed from the computed energy;
	//This is currently called for the backrubs and backbone pairwise minimization, but not for the side-chain dihedral minimization
	private float removeExtraE(boolean templateOnly, int res1Strand, int res1Resnum, int res2Strand, int res2Resnum, int strandMut[][]){
		
		float removeE = 0.0f;
		float tmpE1;
		if (!templateOnly){
			if (res2Resnum != -1) { //not a shell run
				tmpE1 = computeEnergyOfOnlyRes(m.strand[res1Strand].residue[res1Resnum].moleculeResidueNumber,0,strandMut,res1A96ff);
				removeE += tmpE1;
				
				tmpE1 = computeEnergyOfOnlyRes(m.strand[res2Strand].residue[res2Resnum].moleculeResidueNumber,1,strandMut,res2A96ff);
				removeE += tmpE1;
			}
			else { //shell run
				tmpE1 = computeEnergyOfOnlyTemplate(strandMut);
				removeE += tmpE1;
			}
		}
		return removeE;
	}



        private void correctSwitchedTemplates(int strandMut[][], PairwiseEnergyMatrix retEMatrixMin, PairwiseEnergyMatrix retEMatrixMax){
            //This function takes the changes in template energy due to a full structure switch's template changes,
            //which are stored in switchedTemplateE, and transfers them to the intra energies of RCs at the first flexible residue
            //that are in the corresponding perturbation states.  There, they are used for pruning and can be efficiently
            //incorporated into A*.  The two formulations of the energy matrix are biophysically equivalent.


            if( ! m.perts[0].type.equalsIgnoreCase("FULL STRUCTURE SWITCH") )
                return;

            //Figure out which residue's energies will be adjusted
            int str = mutRes2Strand[0];
            int strResNum = strandMut[str][mutRes2StrandMutIndex[0]];

            for(int AA=0; AA<strandRot[str].getNumAllowable(strResNum); AA++){

                int curAA = strandRot[str].getIndexOfNthAllowable(strResNum, AA);
                int numRCs = getNumRot(str, strResNum, curAA);

                for(int curRC=0; curRC<numRCs; curRC++ ){

                    int resPertState = ((StrandRCs)strandRot[str]).RCPertStates[strResNum][curAA][curRC];//Residue perturbation state in this RC
                    int pertState = m.strand[str].residue[strResNum].pertStates[resPertState][0];//Corresponding state of the full structure switch
                    float switchE = switchedTemplateE[pertState];

                    retEMatrixMin.addToShellRotE(0, curAA, curRC, switchE);
                    retEMatrixMax.addToShellRotE(0, curAA, curRC, switchE);
                }
            }

        }


//// END PAIRWISE MUTATION _ALL_ ROTAMER SEARCH SECTION - ENERGY PRECOMPUTATION
///////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////
//// BEGIN ROTAMER SEARCH SECTION
	
	//For the given mutation sequence,
	//Count the total number of conformations, the number of confs pruned by MinDEE, the
	//	number of remaining conformations; the number of conformations above each level
	//	(from level i+1 to the top level, which is the ligand); the number of rotamers for each flexible residue;
	//	the number of rotamers (total and non-pruned by MinDEE) for the flexible residues only;
	//The MinDEE matrix eliminatedRotAtRes[] should already be computed;
	//For each pruned rotamer at curLevel, count the number of *new* conformations that 
	//	are pruned as a result, then decrease the number of rotamers for curLevel to
	//	make sure we do not overcount the pruned conformations (the next pruned rotamer
	//	should not count the conformations that include the first rotamer, and therefore
	//	have already been pruned)
	private int countConfs(int numMutable, int numRes, /*boolean ligPresent,*/ 
			int strandMut[][], int numRotForRes[], int numRotForResNonPruned[]){
		
		int numTotalRotRed = 0;
		numConfsTotal = BigInteger.valueOf(1);
		numConfsAboveLevel = new BigInteger[numRes];

		int curNumRot;
		
		//Store the number of rotamers for each AA in the current mutation
		//The ligand (if present) is in the last level
		for (int curLevel=0; curLevel<numRes; curLevel++){
			int str = mutRes2Strand[curLevel];
			int strResNum = strandMut[str][mutRes2StrandMutIndex[curLevel]];
			int molResNum = m.strand[str].residue[strResNum].moleculeResidueNumber;
			/*if ((ligPresent)&&(curLevel==(numRes-1))) //the ligand level
				curNumRot = grl.getNumRotForAAtype(curLigAANum);
			else*/
			curNumRot = getNumRot( str, strResNum, curAANum[molResNum] );
			numRotForRes[curLevel] = curNumRot;
			numRotForResNonPruned[curLevel] = numRotForRes[curLevel];
			numTotalRotRed += numRotForRes[curLevel];
			numConfsTotal = numConfsTotal.multiply(BigInteger.valueOf(numRotForRes[curLevel]));
		}
		
		BigInteger numPruned = BigInteger.ZERO;
		int numPrunedThisLevel;
		
		//Count the number of rotamers pruned by MinDEE
		for (int curLevel=0; curLevel<numRes; curLevel++){
			numPrunedThisLevel = 0;
			int str = mutRes2Strand[curLevel];
			int strResNum = strandMut[str][mutRes2StrandMutIndex[curLevel]];
			int molResNum = m.strand[str].residue[strResNum].moleculeResidueNumber;
			for (int curRot=0; curRot<numRotForRes[curLevel]; curRot++){
				//int curIndex;
				/*if ((ligPresent)&&(curLevel==(numRes-1))) //the ligand level
					curIndex = numInAS*numTotalRotamers + curRot;
				else*/
				//TODO: fix rotamerIndexOffset so that it works for multiple strands/multiple rotamer libraries
				//curIndex = curLevel*numTotalRotamers + rotamerIndexOffset[curAANum[molResNum]] + curRot;
				
				if (eliminatedRotAtRes.get(curLevel,curAANum[molResNum],curRot)){
					numPruned = numPruned.add(compPrunedConfsByRot(numRotForResNonPruned,numRes,curLevel));
					numPrunedThisLevel++;
				}
			}
			numRotForResNonPruned[curLevel] -= numPrunedThisLevel;
		}
		
		//Count the number of conformations below each level (flexible residue)
		if (numRes>0)
		numConfsAboveLevel[numRes-1] = new BigInteger("1"); //the last level
		if (numRes>1){
			for (int curLevel=numRes-2; curLevel>=0; curLevel--){
				numConfsAboveLevel[curLevel] = numConfsAboveLevel[curLevel+1].multiply(BigInteger.valueOf(numRotForResNonPruned[curLevel+1]));
			}
		}
		
		numConfsPrunedByMinDEE = numPruned; //set the number of confs pruned by MinDEE
		numConfsLeft = numConfsTotal.subtract(numConfsPrunedByMinDEE);
		numConfsPrunedByS = BigInteger.valueOf(0);
		
		return numTotalRotRed;
	}
	//Computes the number of conformations pruned by a rotamer (eliminated by MinDEE) at curLevel:
	//	multiply the number of rotamers for each level different from curLevel
	private BigInteger compPrunedConfsByRot(int numRotForRes[], int numRes, int curLevel){
		
		BigInteger numPruned = BigInteger.ONE;
		
		for (int i=0; i<numRes; i++){
			if (i!=curLevel)
				numPruned = numPruned.multiply(BigInteger.valueOf(numRotForRes[i]));
		}
		
		return numPruned;
	}
	
	//Computes the number of conformations that are pruned by MinDEE due to steric clashes:
	//	for each level, count the number of rotamers that are not pruned by MinDEE due to steric clash, then 
	//	multiply for all levels; the result is the number of conformations that do not include
	//	pruned rotamers; return (totalNumConfs - this number) as the number of conformations pruned by MinDEE
	private BigInteger countPrunedByMinDEESteric(int numMutable, int numRes, /*boolean ligPresent,*/ 
			int strandMut[][], int numRotForRes[], PrunedRotamers<Boolean> prunedIsSteric){
		
		BigInteger numConfNotPruned = BigInteger.ONE;
		int numNonPruned[] = new int[numRes];
		for (int i=0; i<numRotForRes.length; i++)
			numNonPruned[i] = numRotForRes[i];
		
		for (int curLevel=0; curLevel<numRes; curLevel++){
			int str = mutRes2Strand[curLevel];
			int strResNum = strandMut[str][mutRes2StrandMutIndex[curLevel]];
			int molResNum = m.strand[str].residue[strResNum].moleculeResidueNumber;
			for (int curRot=0; curRot<numRotForRes[curLevel]; curRot++){
				//int curIndex;
				/*if ((ligPresent)&&(curLevel==(numRes-1))) //the ligand level
					curIndex = numInAS*numTotalRotamers + curRot;
				else*/
					//curIndex = curLevel*numTotalRotamers + rotamerIndexOffset[curAANum[molResNum]] + curRot;
					if ((eliminatedRotAtRes.get(curLevel,curAANum[molResNum],curRot))&&(prunedIsSteric.get(curLevel,curAANum[molResNum],curRot))){
					numNonPruned[curLevel]--;
				}
			}
		}
		
		for (int curLevel=0; curLevel<numRes; curLevel++)
			numConfNotPruned = numConfNotPruned.multiply(BigInteger.valueOf(numNonPruned[curLevel]));
		
		return (numConfsTotal.subtract(numConfNotPruned));
	}
	//Sets-up the repeat mutation search if the desired accuracy has not been achived, although
	//	all non-pruned conformations have been generated by A*;
	//Updates eliminatedRotAtRes[] to not prune a subset of the pruned rotamers, such that the
	//	number of conformations not pruned by this is at least numPrunedConfsToAllow
	private void setupRepeatRun(BigInteger numPrunedConfsToAllow, int numRotForResNonPruned[], 
			int numLevels, int numMutable){
		
		BigInteger totalNumConfsUnpruned = new BigInteger("0");
		Iterator<RotInfo<Boolean>> elimRotIter = eliminatedRotAtRes.iterator();
		Iterator<RotInfo<Boolean>> prunedStericIter = prunedIsSteric.iterator();
		while(elimRotIter.hasNext()){
		//for (int i=0; i<eliminatedRotAtRes.length; i++){
			RotInfo<Boolean> elimRotInfo = elimRotIter.next();
			RotInfo<Boolean> prunedStericInfo = prunedStericIter.next();
			//check to make sure we are always checking same index
			assert (elimRotInfo.curPos == prunedStericInfo.curPos && elimRotInfo.curAA == prunedStericInfo.curAA &&
					elimRotInfo.curRot == prunedStericInfo.curRot) : "ElimRot and PrunedSteric indexes don't match";
			
			if ((elimRotInfo.state)&&(!prunedStericInfo.state)){ //pruned non-steric rotamer
			//if ((eliminatedRotAtRes[i])&&(!prunedIsSteric[i])){ //pruned non-steric rotamer
				
				int curLevel = elimRotInfo.curPos;//(int)Math.floor(i/numTotalRotamers); //the residue number for the cur rot index
				//TODO: fix that all positions are assumed to have "numTotalRotamers" number of rotamers
				/*if (curLevel>=numInAS) //ligand level (the ligand may have more than numTotalRotamers rotamers)
					curLevel = numInAS; //the residue number for the cur rot index*/
				numRotForResNonPruned[curLevel]++;
				eliminatedRotAtRes.set(elimRotInfo,false);
				
				BigInteger numConfsAtLevelUnpruned = new BigInteger("1"); //count the unpruned confs by unpruning the cur rot index
				for (int j=0; j<numLevels; j++){
					if (j!=curLevel){
						numConfsAtLevelUnpruned = numConfsAtLevelUnpruned.multiply(BigInteger.valueOf(numRotForResNonPruned[j]));
					}
				}
				
				totalNumConfsUnpruned = totalNumConfsUnpruned.add(numConfsAtLevelUnpruned);
				if (totalNumConfsUnpruned.compareTo(numPrunedConfsToAllow)>=0) // num pruned confs reduced by the necessary number
					break;
			}
		}
	}
	// This function performs a rotamer search to compute
	//  a partition function for the slave node.
	// Utilizes a number of helper functions
	public void slaveDoRotamerSearch(int runNum, boolean searchComputeEVEnergy, boolean searchDoMinimization, int numInAS, 
		int strandMut[][], boolean usingInitialBest, BigDecimal initialBest,
		CommucObj cObj, boolean minimizeBB, boolean saveConfs, String fName, boolean doBackrubs, String backrubFile,
		boolean saveTopConfs, boolean printTopConfs, int numTopConfs, int curMut,boolean useMaxKSconfs, BigInteger maxKSconfs) {

		// A rotamer search is performed. For each residue,
		//  every allowable rotamer is tried in combination
		//  with every other allowable rotamer

		computeEVEnergy = searchComputeEVEnergy;
		doMinimization = searchDoMinimization;
		ASAANums = new int[numberMutable];
		curStrRotNum = new int[numberMutable];
		//curASRotNum = new int[numInAS];
		//int curResToASMap[] = new int[m.strand[sysStrNum].numberOfResidues];
			// This map maps the system strand residues back to the AS numbering
			// So something like 8 -> 0, 10 -> 1, 11 -> 2, ...
		//curLigRotNum = 0;
		numConfsPrunedByE = BigInteger.ZERO;
		numConfsPrunedByS = BigInteger.ZERO;
		numConfsEvaluated = BigInteger.ZERO;
		numConfsPrunedByMinDEE = BigInteger.ZERO;
		allPruned = false;
		//confEnergies = new Vector(2048,1024);
		setBestE(9999999.0f);
		bestEUnMin = 9999999.0f;
		if (usingInitialBest)
			initial_q = initialBest.multiply(new BigDecimal((double)(1-KSepsilon)));
		else
			initial_q = new BigDecimal(0.0);
		partial_q = new BigDecimal(0.0);
		partial_p = new BigDecimal(0.0);
		
		//boolean ligPresent = (ligStrNum>=0); //determine if there is a ligand
		
		for(int str=0;str<numberOfStrands;str++){
			for (int i=0; i<strandMut[str].length; i++) //the AS residues are flexible - this is used by simpMin to set up the minimizer
				m.strand[str].residue[strandMut[str][i]].flexible = true;
		}
		

		if (searchDoMinimization && !searchComputeEVEnergy){
			System.out.println("Warning: In order to do minimization computeEVEnergy must be true");
			return;
		}

		// Prepare Amber
		if(searchComputeEVEnergy){
			// Amber should already be loaded
			/*if(ligPresent) {
				a96ff.setLigandNum(m.strand[ligStrNum].residue[0].moleculeResidueNumber);
			}*/
			if (doMinimization){
				if (!minimizeBB){ //side-chain minimization or DEEPer
					if (simpMin == null) {
						System.out.println("Warning: Attempting minimization run but simpMin not allocated, RotamerSearch aborting");
						return;
					}
					bbMin = null;
					//brMin = null;
				}
				else { //backbone minimization
					if (!doBackrubs){
						if (bbMin == null) {
							System.out.println("Warning: Attempting minimization run but bbMin not allocated, RotamerSearch aborting");
							return;
						}
						simpMin = null;
						//brMin = null;
					}
					/*else {
						if (brMin == null) {
							System.out.println("Warning: Attempting minimization run but brMin not allocated, RotamerSearch aborting");
							return;
						}
						simpMin = null;
						bbMin = null;
					}*/
					}
				}
			}

		// Make sure the allRotamerPairsEnergyName matrices exist
		if (arpMatrix == null){
			System.out.println("Warning: allRotamerPairsEnergy matrix not loaded");
			return;
		}
		
		if (eliminatedRotAtRes == null){
			System.out.println("Warning: MinDEE matrix not computed");
			return;
		}

		// Setup the residue number to AS number map
		/*for(int i=0;i<m.strand[sysStrNum].numberOfResidues;i++){
			curResToASMap[i] = -1;
		}
		for(int i=0;i<residueMap.length;i++){
			curResToASMap[residueMap[i]] = i;
		}*/

		/*if (ligPresent)
			slaveMutationRotamerSearch(-1, numInAS, strandMut, AAList, numAAAllowed,
				numTotalRotamers, rotamerIndexOffset,	curResToASMap, ligPresent, minimizeBB, saveConfs, fName, doBackrubs, backrubFile);
		else*/
			slaveMutationRotamerSearch(runNum, 0, numberMutable, strandMut, minimizeBB, 
				saveConfs, fName, doBackrubs, backrubFile,saveTopConfs, printTopConfs, numTopConfs, 
				curMut,useMaxKSconfs, maxKSconfs);
	
		// Store results to the communication object
		if (cObj!=null){
			/*
			if (ligPresent) {
				cObj.EL_searchNumConfsTotal = numConfsTotal.intValue();
				if (allPruned) { //all of the conformations were pruned by MinDEE, as there is a residue with no remaining rotamers
					cObj.EL_searchNumConfsPrunedByS = 0;
					cObj.EL_searchNumPrunedMinDEE = numConfsTotal.intValue();
					cObj.EL_searchNumConfsEvaluated = 0;
					cObj.EL_searchNumConfsLeft = 0;
				}
				else {
					cObj.EL_searchNumConfsPrunedByS = numConfsPrunedByS.intValue();
					cObj.EL_searchNumPrunedMinDEE = numConfsPrunedByMinDEE.intValue();
					cObj.EL_searchNumConfsEvaluated = numConfsEvaluated.intValue();
					cObj.EL_searchNumConfsLeft = numConfsLeft.intValue();
				}
			} else {*/
				cObj.searchNumConfsTotal[runNum] = numConfsTotal.intValue();
				if (allPruned) { //all of the conformations were pruned by MinDEE, as there is a residue with no remaining rotamers
					cObj.searchNumConfsPrunedByS[runNum] = 0;
					cObj.searchNumPrunedMinDEE[runNum] = numConfsTotal.intValue();
					cObj.searchNumConfsEvaluated[runNum] = 0;
					cObj.searchNumConfsLeft[runNum] = 0;
				}
				else {
					cObj.searchNumConfsPrunedByS[runNum] = numConfsPrunedByS.intValue();
					cObj.searchNumPrunedMinDEE[runNum] = numConfsPrunedByMinDEE.intValue();
					cObj.searchNumConfsEvaluated[runNum] = numConfsEvaluated.intValue();
					cObj.searchNumConfsLeft[runNum] = numConfsLeft.intValue();
				}
			//}
	
			/*if (ligPresent)
				cObj.EL_allPruned = allPruned;
			else*/
				cObj.allPruned[runNum] = allPruned;
			
			// Compute q_X
			/*if (ligPresent) {
				cObj.q_EL = partial_q;
				cObj.bestBoundE = (double)bestEUnMin;
				cObj.bestBoundEMin = (double)getBestE();
			} else {*/
				cObj.q[runNum] = partial_q;
				cObj.bestE[runNum] = (double)bestEUnMin;
				cObj.bestEMin[runNum] = (double)getBestE();
				//}
		}
		else {
			/*if (ligPresent)
				System.out.println("Statistics (bound):");
			else*/
				System.out.println("Statistics (unbound):");
			System.out.println("Best Energy:  " + (double)getBestE());
			System.out.println("partial_q: " + partial_q);
			System.out.println("partial_p: " + partial_p);
			System.out.println("NumConfsTotal:      	" + numConfsTotal);
			System.out.println("NumConfsPrunedByMinDEE: " + numConfsPrunedByMinDEE);
			System.out.println("NumConfsPrunedByS:  	" + numConfsPrunedByS);
			System.out.println("NumConfsEvaluated:  	" + numConfsEvaluated);
			System.out.println("NumConfsLeft:       	" + numConfsLeft);
		}
		
	}

	// This function does a mutation search then for each allowable mutation
	//  a simple rotamer search is performed
	private void slaveMutationRotamerSearch(int runNum, int depth, int maxDepth, int 
		strandMut[][], boolean minimizeBB,
		boolean saveConfs, String fName, boolean doBackrubs, String backrubFile, 
		boolean saveTopConfs, boolean printTopConfs, int numTopConfs, int curMut,
		boolean useMaxKSconfs, BigInteger maxKSconfs) {
	
		// If we're at the ligand depth
		/*if (depth == -1) {
			if(ligROT.getNumAllowable(0)==0) {
				curLigAANum = -1;  // this shouldn't happen
				System.out.println("ERROR: Ligand has no allowables but you are using a ligand?");
				System.exit(1);
			}
			else {
				// Because this is a ligand there can be only one
				//  allowable type, ie. it can not mutate
				// Change the residue type to this one type (if protein)
				curLigAANum = ligROT.getIndexOfNthAllowable(0,0);
				if (m.strand[ligStrNum].isProtein)
					ligROT.changeResidueType(m,0,grl.getAAName(curLigAANum),addHydrogens);
				
				slaveMutationRotamerSearch(depth+1,maxDepth,strandMut,AAList,numAAAllowed,
					numTotalRotamers,rotamerIndexOffset,/*curResToASMap,ligPresent,minimizeBB,saveConfs,fName,doBackrubs,backrubFile);
			}
		}
		else {*/
			if (depth >= maxDepth) {
				// If we've arrived here then we're ready to
				//  do a rotamerSearch
				if(debug)
					System.out.println("One Mutation Conf Found and is being tested");
				if (computeEVEnergy){
					a96ff.calculateTypesWithTemplates();
					a96ff.initializeCalculation();
					a96ff.setNBEval(hElect,hVDW);
					if (doMinimization){
						if (!minimizeBB){ //side-chain minimization
							//if(ligPresent)
							simpMin.initialize(m,numberOfStrands,a96ff,strandRot,curAANum,doDihedE);
							//else
							//	simpMin.initialize(m,sysStrNum,a96ff,sysLR,curAANum,doDihedE,rl);
						}
						else { //backbone minimization
							if (!doBackrubs){
								//if (ligPresent)
								//	bbMin.initialize(m, a96ff, residueMap, sysStrNum, ligStrNum);
								//else
									bbMin.initialize(m, a96ff, strandMut, numberOfStrands);
							}
							else {
								brMin.initialize(m, a96ff, strandMut, backrubFile, hSteric, overlapThresh, numberOfStrands,true);
							}
						}
					}
				}

				// Determine the following
				// -AANums of all AS residues (so that we can index properly into rotamerIndexOffset, others)
				//   residue i (i=0..8) has 3 letter code AAList(AANums[i])
				// -total number of conformations for this mutation -> assign to numConfsTotal, numConfsLeft
				// -number of conformations below each level (ie. if we prune one of the 3rd AAs rotamers
				//   because of steric overlap, how many confs are we pruning)
				assignAANums(strandMut);
				//int numConfs = computeTotalNumConfs(residueMap,curResToASMap,ligPresent);
				
				  // also computes num of confs at each level
				//numConfsLeft = (numConfsTotal-numConfsPrunedByMinDEE); //MinDEE should have already been applied
				
				//Perform A* search: compute the partial partition function q*
				slaveRotamerSearchAStar(maxDepth, strandMut, minimizeBB, saveConfs, fName, doBackrubs,
						saveTopConfs, printTopConfs, numTopConfs, useMaxKSconfs, maxKSconfs);
				
				//Print the top structures
				if((saveTopConfs || printTopConfs) && runNum >= 0){
					ConfPair confPairs[] = new ConfPair[topConfs.size()];
					confPairs = topConfs.toArray(confPairs);
					Arrays.sort(confPairs);
					int ctr = 0;
					PrintStream printStream = null;
					String outName = ""+curMut+"_";
					if(printTopConfs){
						for (int str = 0; str<strandMut.length; str++)
							for(int i=0; i<strandMut[str].length; i++){
								outName += RotamerLibrary.getOneLet(m.strand[str].residue[strandMut[str][i]].name);
							}
						
						outName += "_"+runNum;
						try{
							FileOutputStream fileOutputStream = new FileOutputStream(EnvironmentVars.ksConfDir+"/"+outName);
							BufferedOutputStream bufferedOutputStream = new BufferedOutputStream( 
									fileOutputStream);
							printStream = new PrintStream(bufferedOutputStream);
						}
						catch(Exception e){
							e.printStackTrace();
							System.out.println("Couldn't write confOut file.");
						}
					}
					for(int i = confPairs.length-1; i>=0; i--){
						//String filename = String.format("%1$s_%2$03d_%3$d.pdb",fName,ctr++,runNum);
						String filename = String.format("%1$s_%2$03d.pdb",outName,ctr++);
						if(saveTopConfs){
							saveConf(confPairs[i].conf, confPairs[i].energy[0], EnvironmentVars.ksConfDir+"/"+filename,strandMut,minimizeBB, doBackrubs);
							System.out.println("Saved Conf " + i);
						}
						if(printStream != null){ //printTopConfs
							printStream.print(""+ctr+" ");
							for (int str = 0; str<strandMut.length; str++)
								for(int j=0; j<strandMut[str].length; j++){
									printStream.print(m.strand[str].residue[strandMut[str][j]].name+" ");
								}
							for(int j: confPairs[i].conf)
								printStream.print(j+" ");
							
							printStream.print("unMin: "+confPairs[i].energy[1]+" minE: "+confPairs[i].energy[0]+" ");
							
							printStream.println("");
						}
					}
					if(printTopConfs || printStream != null){
						printStream.close();
					}	
				}
				
				return;
			}

			// If there are no allowables then test with 'native' form of the enzyme
			int str = mutRes2Strand[depth];
			int strResNum = strandMut[str][mutRes2StrandMutIndex[depth]];
			int molResNum = m.strand[str].residue[strResNum].moleculeResidueNumber;
			if (strandRot[str].getNumAllowable(strResNum) == 0) {
				curAANum[molResNum] = -1;
				slaveMutationRotamerSearch(runNum, depth+1,maxDepth,strandMut,
					minimizeBB,saveConfs,fName,doBackrubs,backrubFile,
					saveTopConfs, printTopConfs, numTopConfs, curMut,
					useMaxKSconfs, maxKSconfs);
			}
			else {
				// Otherwise check with different AAs
				for(int q=0;q<strandRot[str].getNumAllowable(strResNum);q++) {
					curAANum[molResNum] = strandRot[str].getIndexOfNthAllowable(strResNum,q);
					// Change the residue type
					if(m.strand[str].isProtein)
						strandRot[str].changeResidueType(m,strResNum,strandRot[str].rl.getAAName(curAANum[molResNum]),addHydrogens,connectResidues);
					slaveMutationRotamerSearch(runNum, depth+1,maxDepth,strandMut,
						minimizeBB,saveConfs,fName,doBackrubs,backrubFile,
						saveTopConfs, printTopConfs, numTopConfs, curMut,
						useMaxKSconfs, maxKSconfs);
				}
			}
		}
	// Calls AStar repeatedly while the returned conformations still have energy below the threshold;
	//		computes the partial partition function q*;
	// Called by slaveMutationRotamerSearch(.)
	private void slaveRotamerSearchAStar(int numMutable, int strandMut[][], boolean minimizeBB,
			boolean saveConfs, String fName, boolean doBackrubs,boolean saveTopConfs, 
			boolean printTopConfs, int numTopConfs, boolean useMaxKSconfs, BigInteger maxKSconfs){

		/*///////////////////////////////////////////////////////////////////////
		PrintStream logPS = null;
		try {			
			FileOutputStream fileOutputStream = new FileOutputStream("/net/grad/shaqbuzz/i+"+numRotSamples+".txt");
			BufferedOutputStream bufferedOutputStream = new BufferedOutputStream( fileOutputStream );
			logPS = new PrintStream( bufferedOutputStream );
			numRotSamples++;
		}
		catch (Exception ex) {
			System.out.println("ERROR: An exception occured while opening log file");
		}
		//for (int i=0;i<9;i++){logPS.println(curAANum[residueMap[i]]);}logPS.flush();
		///////////////////////////////////////////////////////////////////////*/

                if(doPerturbations)//Make sure minimizer is set properly
                    ((PMinimizer)simpMin).minimizePerturbations = minimizePerturbations;



		if(saveTopConfs || printTopConfs)
			topConfs = new PriorityQueue<ConfPair>(numTopConfs);
		
		ExpFunction ef = new ExpFunction();
		
		
		int treeLevels; // total num levels in the conformation tree
		/*if (ligPresent) //determine num tree levels: if ligPresent, then numInAS+1
			treeLevels = numInAS+1;
		else*/
		treeLevels = numMutable;
		
		int numRotForRes[] = new int[treeLevels]; //the number of rotamers for each flexible residue (AS+lig) during a mutation search
		int numRotForResNonPruned[] = new int[treeLevels]; //the number of non-pruned (by MinDEE) rotamers for each flexible residue
		int numTotalRotRed = 0;		//the total number of rotamers for the flexible residues only (during a mutation search)
		int numTotalRotRedNonPruned = 0; //the total num of non-pruned rotamers for the flexible residues
		
		
		//Count the total number of conformations, the number of conformations pruned by MinDEE,
		//	the remaining conformations; the total num rotamers for the flexible residues, the num
		//	rotamers (total and non-pruned by MinDEE) for each fllexible residue;
		numTotalRotRed = countConfs(numMutable, treeLevels, strandMut, numRotForRes, numRotForResNonPruned);
		
		BigInteger k_const = numConfsPrunedByMinDEE;
	
		BigInteger numConfsPrunedMinDEESteric = countPrunedByMinDEESteric(numMutable, treeLevels, /*ligPresent,*/ 
				strandMut, numRotForRes, prunedIsSteric);
		
		k_const = k_const.subtract(numConfsPrunedMinDEESteric); //only the non-steric prunings are used in the computation of k_const
		
		//Bound the contribution of the conformations pruned by MinDEE
		BigDecimal pStar = ef.exp(-Ec_const/constRT).multiply(new BigDecimal(k_const));
		
		final double ro = (double)KSepsilon /(double)(1-KSepsilon);
		
		System.out.println("k_const: "+k_const+" pStar: "+printBigNum(pStar,5)+" numConfsPrunedMinDEESteric: "+numConfsPrunedMinDEESteric);
		
		//Count the number of non-pruned rotamers
		for (int curLevel=0; curLevel<treeLevels; curLevel++){
			if (numRotForResNonPruned[curLevel]==0){ //no non-pruned rotamers for curLevel, so no possible conformations
				allPruned = true;
				BigDecimal e = ef.exp(-Ec_const/constRT);
				if ( (!k_const.equals(BigInteger.ZERO)) && (e.compareTo(BigDecimal.ZERO)!=0) ) { //some non-sterics pruned but accuracy not achieved, so search must be repeated
					BigDecimal psi = initial_q.max(partial_q.multiply(new BigDecimal(ro)));
					
					BigDecimal f = psi.divide(e,4); //rounding is ROUND_HALF_UP
					BigInteger l_const = k_const.subtract( BigInteger.valueOf( (long)Math.ceil(f.doubleValue()) ) );
					setupRepeatRun(l_const, numRotForResNonPruned, treeLevels, numMutable); //accuracy not achieved, so repeat the search with reduced num pruned conf MinDEE
					repeatSearch = true;
				}
				
				return;
			}
			else
				numTotalRotRedNonPruned += numRotForResNonPruned[curLevel];
		}
		
		//logPS.println(Ec_const+" "+numConfsPrunedByMinDEE+" "+numConfsPrunedMinDEESteric+" "+k_const+" "+pStar+" "+numTotalRotRedNonPruned);
		

                boolean[] eliminatedRotAtPosRed = new boolean[numTotalRotRed];//reduced MinDEE matrix


                boolean splitFlagsRed[][] = null, tripleFlagsRed[][][] = null;
                if( useFlagsAStar ){
                    splitFlagsRed = new boolean[numTotalRotRedNonPruned][numTotalRotRedNonPruned];
                    if( useTriples )
                        tripleFlagsRed = new boolean[numTotalRotRedNonPruned][][];//As with the normal tripleFlags, tripleFlagsRed[a][b][c] is only defined for a>b>c
                }



                //Get the reduced min energy matrix
                ReducedEnergyMatrix arpMatrixRed = arpMatrix.reduceMatrix( eliminatedRotAtPosRed,
			numRotForRes, numRotForResNonPruned, treeLevels,
			numTotalRotRedNonPruned, numMutable, strandMut,
                        numTotalRotRed, this, true, splitFlagsRed, tripleFlagsRed );



		//Set-up the A* search
		StericCheck stericF = null;
		if (!minimizeBB) {//do not do a steric check if backbone minimization
			/*if (ligPresent)
				stericF = new StericCheck(curAANum,curResToASMap,residueMap,eliminatedRotAtPosRed,numRotForRes,
						m,softOverlapThresh,hSteric,numConfsLeft,numConfsAboveLevel,sysStrNum,sysLR,ligStrNum,ligROT,curLigAANum,rl,grl);
			else*/
				stericF = new StericCheck(curAANum,/*curResToASMap,*/strandMut,eliminatedRotAtPosRed,numRotForRes,
						m,softOverlapThresh,hSteric,numConfsLeft,numConfsAboveLevel,numMutable,numberOfStrands,
                                                strandRot,mutRes2Strand,mutRes2StrandMutIndex,doPerturbations);
		}
		
		MSAStar AStarSearch = new MSAStar(treeLevels,numRotForResNonPruned,arpMatrixRed,stericF,
                        splitFlagsRed,tripleFlagsRed,doPerturbations);
		
		curConf = new int[treeLevels]; //the rotamer sequence				
		boolean run1 = true;

		while (numConfsLeft.compareTo(BigInteger.ZERO)==1){

			curConf = AStarSearch.doAStar(run1,numMutable,null,eliminatedRotAtPosRed,strandRot,null,numRotForRes,strandMut,true,mutRes2Strand,mutRes2StrandMutIndex); //the current rotamer sequence); //the current rotamer sequence
			run1 = false;

			//Update the number of remaining conformations
			if (stericF!=null){
				numConfsLeft = stericF.getNumConfsLeft();
				numConfsPrunedByS = stericF.getNumConfsPrunedByS();
			}
			
			for (int curRotCheck=0; curRotCheck<treeLevels; curRotCheck++){//check if the conformation is valid
				if (curConf[curRotCheck]==-1){ // no valid conformations remaining
					if (partial_q.multiply(new BigDecimal(ro)).compareTo(pStar)<0){ //approximation accuracy not achieved
						BigDecimal e = ef.exp(-Ec_const/constRT);
						if ( (!k_const.equals(BigInteger.ZERO)) && (e.compareTo(BigDecimal.ZERO)!=0) ){ //some non-sterics pruned but accuracy not achieved, so search must be repeated
							BigDecimal psi = initial_q.max(partial_q.multiply(new BigDecimal(ro)));
							
							BigDecimal f = psi.divide(e,4); //rounding is ROUND_HALF_UP
							BigInteger l_const = k_const.subtract(BigInteger.valueOf((long)Math.ceil(f.doubleValue())));
							setupRepeatRun(l_const, numRotForResNonPruned, treeLevels, numMutable); //accuracy not achieved, so repeat the search with reduced num pruned conf MinDEE
							repeatSearch = true;
						}
					}
					
					return;
				}
			}
			
			//As the rotamers given to A* are only the non-pruned ones, there is a difference between the
			//	rotamer numbers returned by A* and the actual rotamer numbers for each residue (that is,
			//	A* may return rot 4 for res 3, but rot 3 for res 3 may be pruned, and so the actual number
			//	of the rot to be applied for res 3 is 5)
			int conf[] = new int[curConf.length]; //the conformation with the actual rotamer numbers
			conf = getActualConf(curConf,eliminatedRotAtPosRed,treeLevels,numRotForRes,conf);
			
			m.backupAtomCoord();
			applyRotamers(strandMut, conf);
			
			/***
			 * for (int curLevel=0; curLevel<m.strand[sysStrNum].numberOfResidues; curLevel++){
				if (curResToASMap[curLevel]!=-1){//make a change only to the AS residues: use the native type for the other residues
										
					if (rl.getNumRotForAAtype(curAANum[residueMap[curAS]])!=0){//not GLY or ALA
						sysLR.applyRotamer(m, curLevel, conf[curAS]);
						curASRotNum[curResToASMap[curLevel]] = conf[curAS];
					}
					else { //GLY or ALA
						curASRotNum[curResToASMap[curLevel]] = 0;
					}
					curAS++; //prepare the next AS residue
				}
			}		
			if (ligPresent){ //apply the ligand rotamer
				if (grl.getNumRotForAAtype(curLigAANum)!=0){//not GLY or ALA
					ligROT.applyRotamer(m, 0, conf[treeLevels-1]);//the ligand level
					curLigRotNum = conf[treeLevels-1];
				}
				else { //GLY or ALA
					curLigRotNum = 0;
				}
			}***/
			
			for(int i=0;i<treeLevels;i++)System.out.print(conf[i]+" ");System.out.println();
			
			//Check the energy of the conformation and compute the score if necessary
			double minELowerBound = (double)(computeBestRotEnergyBound(/*numTotalRotamers,rotamerIndexOffset*/));
			//double psi = Math.max(initial_q,partial_q);
			BigDecimal psi = initial_q.max(partial_q.multiply(new BigDecimal(ro)));
			
			//double curThreshold = -constRT * (Math.log(psi)+Math.log(ro/numConfsLeft));
			double curThreshold = stericE;
			BigDecimal diff_qp = psi.subtract(pStar);
			if (diff_qp.compareTo(new BigDecimal(0.0))<0) { //the contribution of the pruned confs is bigger than ro*partial_q, so the search cannot be halted
				
				BigDecimal qBound = partial_q.add(ef.exp(-minELowerBound/constRT).multiply(new BigDecimal(numConfsLeft))); //an upper bound on what partial_q can be
				
				if (pStar.compareTo(initial_q.max(qBound.multiply(new BigDecimal(ro))))>0){ //approximation accuracy cannot be achieved
					
					BigDecimal e = ef.exp(-Ec_const/constRT);
					if ( (!k_const.equals(BigInteger.ZERO)) && (e.compareTo(BigDecimal.ZERO)!=0) ){ //some non-sterics pruned but accuracy not achieved, so search must be repeated						
						
						BigDecimal f = psi.divide(e,4); //rounding is ROUND_HALF_UP
						BigInteger l_const = k_const.subtract(BigInteger.valueOf((long)Math.ceil(f.doubleValue())));
						setupRepeatRun(l_const, numRotForResNonPruned, treeLevels, numMutable); //accuracy not achieved, so repeat the search with reduced num pruned conf MinDEE
						repeatSearch = true;
						return;
					}
					else { //approximation accuracy not achieved and nothing to un-prune
						System.out.println("ERROR: Approximation accuracy not achieved, but no rotamers to unprune..");
						System.exit(1);
					}
				}
				else //it may be possible to achieve approximation accuracy
					curThreshold = stericE;
			}
			else
				curThreshold = -constRT * (ef.log(diff_qp).doubleValue()-Math.log(numConfsLeft.doubleValue()));

			System.out.println("conf: "+numConfsEvaluated.add(BigInteger.ONE)+" minELowerBound: "+minELowerBound+" curThreshold: "+curThreshold);
			System.out.println("pStar: "+printBigNum(pStar,3)+" qStar: "+printBigNum(partial_q,3)+" rho*qStar: "+printBigNum(partial_q.multiply(new BigDecimal(ro)),3));
			
			if (minELowerBound > curThreshold) 
				return;
			else if(useMaxKSconfs && numConfsEvaluated.compareTo(maxKSconfs) > 0)
				return;
			
			else { //the energy of the new conformation is still below the threshold
				
				// After minimization do a m.updateCoordinates() to
				//  resync the actualCoordinates which were changed
				//  in the minimization procedure
				numConfsEvaluated = numConfsEvaluated.add(BigInteger.ONE);
				numConfsLeft = numConfsLeft.subtract(BigInteger.ONE);
				if (stericF!=null)
					stericF.setNumConfsLeft(numConfsLeft); //update the number of conformations still to examine
				float energy = 0.0f;
				float unMinE = 0.0f;
				if (doMinimization){
					if (!minimizeBB){ //side-chain minimization
						energy = calcTotalSnapshotEnergy();
						unMinE = energy;
						if (energy < bestEUnMin) {
							bestEUnMin = energy;
						}						
						simpMin.minimize(numMinSteps);
						energy = calcTotalSnapshotEnergy();
						if (doDihedE) //add dihedral energies
							energy += simpMin.computeDihedEnergy();
					}
					else { //backbone minimization
						energy = calcTotalSnapshotEnergy();
						unMinE = energy;
						if (energy < bestEUnMin) {
							bestEUnMin = energy;
						}
						if (!doBackrubs)
							bbMin.minimizeFull(false);
						else
							brMin.minimizeFull();
						energy = calcTotalSnapshotEnergy();
					}
					
					if (saveConfs)
						saveMolecule(m,(fName+numConfsEvaluated+".pdb"),energy);
				}
				else if (computeEVEnergy){
					energy = (float)minELowerBound;
					unMinE = energy;
					if (energy < bestEUnMin) {
						bestEUnMin = energy;
					}
				}
				
                                m.updateCoordinates();
                                if(doPerturbations)
                                    m.revertPertParamsToCurState();
				updateBestE(energy);
				
				partial_q = partial_q.add(ef.exp(-((double)(energy)) / constRT));
				
				if ((saveTopConfs || printTopConfs)){
					float[] myEnergy = new float[2];
					myEnergy[0] = energy;
					myEnergy[1] = unMinE;
					ConfPair cp = new ConfPair(conf, myEnergy);
					ConfPair head = topConfs.peek();
					if(topConfs.size() >= numTopConfs){
						if(cp.energy[0] < head.energy[0]){
							topConfs.add(cp);
							topConfs.remove(); //Will remove head of queue, so want to traverse
								//Backwards when generating conformations
						}
					}
					else
						topConfs.add(cp);
				}
				
				System.out.println("energy: "+energy);
			}
		}
		if ((numConfsLeft.equals(BigInteger.ZERO))&&(!k_const.equals(BigInteger.ZERO))){ //no conformations remaining, non-sterics pruned			
			if (partial_q.multiply(new BigDecimal(ro)).compareTo(pStar)<0){ //approximation accuracy not achieved, so repeat search
				BigDecimal psi = initial_q.max(partial_q.multiply(new BigDecimal(ro)));
				BigDecimal e = ef.exp(-Ec_const/constRT);
				if (e.compareTo(BigDecimal.ZERO)!=0){
					BigDecimal f = psi.divide(e,4); //rounding is ROUND_HALF_UP				
					BigInteger l_const = k_const.subtract(BigInteger.valueOf((long)Math.ceil(f.doubleValue())));
					setupRepeatRun(l_const, numRotForResNonPruned, treeLevels, numMutable); //accuracy not achieved, so repeat the search with reduced num pruned conf MinDEE
					
					repeatSearch = true;
				}
			}
		}
	}
	
	public void applyRotamers(int[][] strandMut, int[] conf) {
		//Apply the rotamers of the current conformation
		int curAS = 0;	
		
		for (int str=0; str<numberOfStrands;str++){
			for(int i=0; i<strandMut[str].length;i++){

                                if(doPerturbations){
                                    boolean validRC = ((StrandRCs)strandRot[str]).applyRC(m, strandMut[str][i], conf[curAS]);
                                    if(!validRC)
                                        System.err.println("Error: invalid RC " + conf[curAS] + " at residue " + 
                                                strandMut[str][i] + " of strand " + str );

                                    curStrRotNum[curAS] = conf[curAS];
                                }
                                else{
                                    
                                    int molResNum = m.strand[str].residue[strandMut[str][i]].moleculeResidueNumber;

                                    if (strandRot[str].rl.getNumRotForAAtype(curAANum[molResNum])!=0){//not GLY or ALA
                                            strandRot[str].applyRotamer(m, strandMut[str][i], conf[curAS]);
                                            curStrRotNum[curAS] = conf[curAS];
                                    }
                                    else { //GLY or ALA
                                            curStrRotNum[curAS] = 0;
                                    }
                                }
                                
				curAS++; //prepare the next AS residue
			}
		}

	}
	
	
	//Get the actual numbers for the rotamers of a conformation that is returned by A* by including the information
	//		for the pruned rotamers
	private int [] getActualConf(int curConf[], boolean eliminatedRotAtPosRed[], 
			int treeLevels, int numRotForRes[], int conf[]){
		
		int curPruningInd = 0;
		for (int curLevel=0; curLevel<treeLevels; curLevel++){
			int curRotInd = 0;
			for (int curRot=0; curRot<numRotForRes[curLevel]; curRot++){
				if (!eliminatedRotAtPosRed[curPruningInd]){
					if (curRotInd==curConf[curLevel])
						conf[curLevel] = curRot;
					curRotInd++;
				}
				curPruningInd++;
			}
		}
		return conf;
	}
	
//// END ROTAMER SEARCH SECTION
/////////////////////////////////////////////////////////////////////////////

	
/*
 * 
 * DEE section
 * 
 */		
	//Compute the two interval terms in the summation of the MinDEE criteria
	public void doCompMinDEEIntervals(int numMutable, int strandMut[][], PrunedRotamers<Boolean> prunedRotAtRes,
			boolean scaleInt, float maxIntScale){
		
		System.out.print("Computing MinDEE interval terms..");	
		
		float dist[][] = null;
		if (scaleInt)
			dist = doCompDistSC(numMutable,strandMut);
		
		MinDEEIntervals compInt = new MinDEEIntervals(arpMatrix, arpMatrixMax, numMutable, 	strandMut, 
				strandRot, prunedRotAtRes, scaleInt, dist, maxIntScale, mutRes2Strand,mutRes2StrandMutIndex, doPerturbations);
		
		compInt.compMinDEEIntervals();
		indIntMinDEE = compInt.getIndInt();
		pairIntMinDEE = compInt.getPairInt();
		
		System.out.println("done.");
		
		if (debug){
			System.out.println();
			System.out.print(" ind: ");
			for (int curPos=0; curPos<indIntMinDEE.length; curPos++){
				System.out.print(indIntMinDEE[curPos]+" ");
			}
			System.out.println();
			System.out.print(" pair: ");
			for (int curPos=0; curPos<pairIntMinDEE.length; curPos++){
				System.out.print(pairIntMinDEE[curPos]+" ");
			}
			System.out.println();
		}
	}
	//Compute the distance between the side-chains for each active site residue pair (the ligand is not considered here);
	//	Returns the minimum distance between a pair of atoms in the two side-chains, for each side-chain pair
	private float[][] doCompDistSC(int numMutable, int strandMut[][]){
		
		float dist[][] = new float[numMutable][numMutable];
		
		for (int str1=0;str1<numberOfStrands;str1++){
			for(int i=0;i<strandMut[str1].length;i++){
				Residue r1 = m.strand[str1].residue[strandMut[str1][i]];
				for(int str2=str1;str2<numberOfStrands;str2++){
					for (int j=i+1; j<strandMut[str2].length; j++){
						Residue r2 = m.strand[str2].residue[strandMut[str2][j]];
						dist[i][j] = r1.getDist(r2,false);
						dist[j][i] = dist[i][j];
					}
				}
			}
		}
		
		/*for (int i=0; i<numResInActiveSite; i++){
			Residue r1 = m.strand[sysStrNum].residue[residueMap[i]];
			for (int j=i+1; j<numResInActiveSite; j++){
				Residue r2 = m.strand[sysStrNum].residue[residueMap[j]];
				dist[i][j] = r1.getDist(r2,false);
				dist[j][i] = dist[i][j];
			}
		}*/
		return dist;
	}
	//Prune all rotamers that are incompatible with the template (intra E + res-to-template E >= stericE)
	public PrunedRotamers<Boolean> DoPruneStericTemplate(int numMutable, int strandMut[][],
			PrunedRotamers<Boolean> prunedRotAtRes, double stericE){
		
		eliminatedRotAtRes = prunedRotAtRes;
		
		if (prunedIsSteric==null){ //no pruning runs done yet
			prunedIsSteric = new PrunedRotamers<Boolean>(numMutable,strandMut, this, false);//boolean[prunedRotAtRes.length];
			//for (int i=0; i<prunedIsSteric.length; i++)
			//	prunedIsSteric[i] = false;
		}
		int numAAtypes[] = new int[numMutable];
		
		int ctr=0;
		for(int str=0;str<strandMut.length;str++){ //the number of AAs allowed for each AS residue
			for(int i=0;i<strandMut[str].length;i++){
				numAAtypes[ctr] = strandRot[str].getNumAllowable(strandMut[str][i]);
				ctr++;
			}
		}
		//int numAAtypes[] = new int[numMutable];
		//for (int i=0; i<numAAtypes.length; i++) //the number of AAs allowed for each AS residue
		//	numAAtypes[i] = sysLR.getNumAllowable(residueMap[i]);
		
		int numPruned = 0;
		
		//Compute for the AS residues first
		for (int curPos=0; curPos<numMutable; curPos++){
			int str=mutRes2Strand[curPos];
			int strResNum=strandMut[str][mutRes2StrandMutIndex[curPos]];
			for (int AA=0; AA<numAAtypes[curPos]; AA++){
				int curAA = strandRot[str].getIndexOfNthAllowable(strResNum,AA);
				
				//find how many rotamers are allowed for the current AA type at the given residue;
				//note that ala and gly have 0 possible rotamers
				int numRotForCurAAatPos = getNumRot( str, strResNum, curAA );
				
				for(int curRot=0; curRot<numRotForCurAAatPos; curRot++){
					try{
					if ( arpMatrix.getIntraAndShellE(curPos, curAA, curRot) >= stericE ){
						
						//int index_r = curPos*totalNumRotamers + rotamerIndexOffset[curAA] + curRot;
						
						eliminatedRotAtRes.set(curPos,curAA,curRot, true);
						
						prunedIsSteric.set(curPos, curAA, curRot, true);
						numPruned++;
					}
				}
					catch(Exception e){
						System.out.println(curPos + " "+ curAA+" "+curRot);
			}
		}
			}
		}
		
		//If there is a ligand, compute for the lig rotamers as well
		/*if (numLigRotamers!=0){
			for (int curRot=0; curRot<numLigRotamers; curRot++){
				
				int ligAAind = ligROT.getIndexOfNthAllowable(0,0);
				
				if ((arpMatrix[numInAS][ligAAind][curRot][numInAS][0][0]+arpMatrix[numInAS][ligAAind][curRot][numInAS][0][1])>=stericE){
					
					int index_r = numInAS*totalNumRotamers+curRot;
					
					eliminatedRotAtRes[index_r] = true;
					
					prunedIsSteric[index_r] = true;
					numPruned++;
				}
			}
		}*/
		
		System.out.println("Number of rotamers pruned due to incompatibility with the template: "+numPruned);		
		
		return eliminatedRotAtRes;
	}



        //Prune rotamer or residue-conformation pairs that are definitely impossible due to a steric clash or to parametric incompatibility
        //as measured by their interaction energy being greater than the cutoff
        //Particularly useful for DEEPer
        public void pruneRidiculousPairs(int numMutable, int strandMut[][],
			PrunedRotamers<Boolean> prunedRotAtRes, float cutoff){

		eliminatedRotAtRes = prunedRotAtRes;

		int numAAtypes[] = new int[numMutable];

		int ctr=0;
		for(int str=0;str<strandMut.length;str++){ //the number of AAs allowed for each AS residue
			for(int i=0;i<strandMut[str].length;i++){
				numAAtypes[ctr] = strandRot[str].getNumAllowable(strandMut[str][i]);
				ctr++;
			}
		}
                
		int numPruned = 0;

                for (int curPos=0; curPos<numMutable; curPos++){
			int str=mutRes2Strand[curPos];
			int strResNum=strandMut[str][mutRes2StrandMutIndex[curPos]];
			for (int AA=0; AA<numAAtypes[curPos]; AA++){
				int curAA = strandRot[str].getIndexOfNthAllowable(strResNum,AA);

				//find how many rotamers are allowed for the current AA type at the given residue;
				//note that ala and gly have 0 possible rotamers
				int numRotForCurAAatPos = getNumRot( str, strResNum, curAA );


				for(int curRot=0; curRot<numRotForCurAAatPos; curRot++){

                                    if( !eliminatedRotAtRes.get(curPos, curAA, curRot)){

                                        for (int altPos=0; altPos<curPos; altPos++){
                                            int str2=mutRes2Strand[altPos];
                                            int strResNum2=strandMut[str2][mutRes2StrandMutIndex[altPos]];
                                            for (int AA2=0; AA2<numAAtypes[altPos]; AA2++){
                                                    int altAA = strandRot[str2].getIndexOfNthAllowable(strResNum2,AA2);

                                                    //find how many rotamers are allowed for the current AA type at the given residue;
                                                    //note that ala and gly have 0 possible rotamers
                                                    int numRotForCurAAatPos2 = getNumRot( str2, strResNum2, altAA );


                                                    for(int altRot=0; altRot<numRotForCurAAatPos2; altRot++){


                                                        if( !eliminatedRotAtRes.get(altPos, altAA, altRot)){

                                                            if( arpMatrix.getPairwiseE(curPos, curAA, curRot, altPos, altAA, altRot) >= cutoff ){

                                                                splitFlags[curPos][curAA][curRot][altPos][altAA][altRot] = true;
                                                                splitFlags[altPos][altAA][altRot][curPos][curAA][curRot] = true;

                                                                numPruned++;
                                                            }
                                                        }
                                                    }
                                            }
                                        }
                                    }
                                }
                        }
                }


                if(doPerturbations)
                    System.out.println("Number of pairs pruned due to steric clashes or parametric incompatibility: "+numPruned);
                else
                    System.out.println("Number of pairs pruned due to steric clashes: "+numPruned);
	}


        //Remove RCs that are parametrically incompatible with all RCs at some other position (i.e. are impossible)
        //This can be done before energy precomputation to avoid precomputing useless energies
        //The RCs are actually removed from sysLR and ligROT
        public void removeImpossibleRCs(int numMutable, int strandMut[][]){

            boolean done = false;

            while( !done ){//We iterate until no more RCs can be removed

                    done = true;
                    
                    int numAAtypes[] = new int[numMutable];
                    
                    boolean prunedRCs[][][][] = new boolean[numberOfStrands][][][];//Which RCs should be removed
                    for(int str=0; str<numberOfStrands; str++)
                        prunedRCs[str] = new boolean[m.strand[str].numberOfResidues][strandRot[str].rl.getNumAAallowed()][];

                    int ctr=0;
                    for(int str=0;str<strandMut.length;str++){ //the number of AAs allowed for each AS residue
                            for(int i=0;i<strandMut[str].length;i++){
                                    numAAtypes[ctr] = strandRot[str].getNumAllowable(strandMut[str][i]);
                                    ctr++;
                            }
                    }

                    for (int curPos=0; curPos<numMutable; curPos++){
                            int str=mutRes2Strand[curPos];
                            int strResNum=strandMut[str][mutRes2StrandMutIndex[curPos]];
                            Residue res = m.strand[str].residue[strResNum];

                            for (int AA=0; AA<numAAtypes[curPos]; AA++){
                                    int curAA = strandRot[str].getIndexOfNthAllowable(strResNum,AA);

                                    //find how many rotamers are allowed for the current AA type at the given residue;
                                    //note that ala and gly have 0 possible rotamers
                                    int numRCForCurAAatPos = getNumRot( str, strResNum, curAA );

                                    prunedRCs[str][strResNum][curAA] = new boolean[numRCForCurAAatPos];

                                    for(int curRC=0; curRC<numRCForCurAAatPos; curRC++){

                                        prunedRCs[str][strResNum][curAA][curRC] = false;

                                            for (int altPos=0; altPos<numMutable; altPos++){

                                                if( altPos != curPos ){

                                                    int str2=mutRes2Strand[altPos];
                                                    int strResNum2=strandMut[str2][mutRes2StrandMutIndex[altPos]];

                                                    boolean prune = true;
                                                    Residue res2 = m.strand[str2].residue[strResNum2];
                                                    
                                                    for (int AA2=0; AA2<numAAtypes[altPos]; AA2++){
                                                            int altAA = strandRot[str2].getIndexOfNthAllowable(strResNum2,AA2);

                                                            //find how many rotamers are allowed for the current AA type at the given residue;
                                                            //note that ala and gly have 0 possible rotamers
                                                            int numRCForCurAAatPos2 = getNumRot( str2, strResNum2, altAA );


                                                            for(int altRC=0; altRC<numRCForCurAAatPos2; altRC++){

                                                                if( ! isParametricallyIncompatible(res, curAA, curRC, res2, altAA, altRC))
                                                                    prune = false;
                                                            }
                                                    }

                                                    prunedRCs[str][strResNum][curAA][curRC] = prunedRCs[str][strResNum][curAA][curRC] || prune;
                                                }
                                            }

                                        if(prunedRCs[str][strResNum][curAA][curRC])
                                            done = false;
                                    }
                            }
                    }


                    for(int str=0; str<numberOfStrands; str++)
                        ((StrandRCs)strandRot[str]).removeRCs(prunedRCs[str]);
            }
        }



	//Marks all rotamer pairs for which the min energy matrix entry is greater than cutoff as having a steric clash;
	//		If the max energy matrix exists, the corresponding entries are marked in the same way
	public void preprocessPairs(float cutoff, int numMutable, int strandMut[][]){
		
		int numLevels = numMutable;
		/*if (ligPresent)
			numLevels++;*/
		
		if (arpMatrix==null)
			return;
		else {			
			for (int curLevel1=0; curLevel1<numMutable; curLevel1++){ //the ligand level is taken into account in curLevel2
				int str1=mutRes2Strand[curLevel1];
				int strResNum1=strandMut[str1][mutRes2StrandMutIndex[curLevel1]];
				for (int curAA1=0; curAA1<strandRot[str1].getNumAllowable(strResNum1); curAA1++){ //for all allowed AA's
					int index1 = strandRot[str1].getIndexOfNthAllowable(strResNum1,curAA1);
					int newRot1 = getNumRot(str1, strResNum1, index1);
                                        
					for (int curRot1=0; curRot1<newRot1; curRot1++){ //for all rotamers for the given AA
						
						for (int curLevel2=(curLevel1+1); curLevel2<numLevels; curLevel2++){
							int str2=mutRes2Strand[curLevel2];
							int strResNum2=strandMut[str2][mutRes2StrandMutIndex[curLevel2]];	
							if (curLevel1!=curLevel2) {
								
								/*if ((ligPresent)&&(curLevel2==numInAS)){ //the ligand level
									int ligAAind = ligROT.getIndexOfNthAllowable(0,0);
									int newRot2 = grl.getNumRotForAAtype(ligAAind);
									if (newRot2==0)
										newRot2 = 1;
									for (int curRot2=0; curRot2<newRot2; curRot2++){ //for all rotamers for the given AA
										if (arpMatrix[curLevel1][index1][curRot1][curLevel2][ligAAind][curRot2]>cutoff){
											arpMatrix[curLevel1][index1][curRot1][curLevel2][ligAAind][curRot2] = stericE;
											arpMatrix[curLevel2][ligAAind][curRot2][curLevel1][index1][curRot1] = stericE;
											if (arpMatrixMax!=null){
												arpMatrixMax[curLevel1][index1][curRot1][curLevel2][ligAAind][curRot2] = stericE;
												arpMatrixMax[curLevel2][ligAAind][curRot2][curLevel1][index1][curRot1] = stericE;
											}
										}
									}
								}
								else {*/							
									for (int curAA2=0; curAA2<strandRot[str2].getNumAllowable(strResNum2); curAA2++){ //for all allowed AA's
										int index2 = strandRot[str2].getIndexOfNthAllowable(strResNum2,curAA2);
										int newRot2 = getNumRot( str2, strResNum2, index2 );
										for (int curRot2=0; curRot2<newRot2; curRot2++){ //for all rotamers for the given AA
											if ( arpMatrix.getPairwiseE( curLevel1, index1, curRot1, curLevel2, index2, curRot2 ) > cutoff ){
												arpMatrix.setPairwiseE( curLevel1, index1, curRot1, curLevel2, index2, curRot2, stericE );
												if (arpMatrixMax!=null)
                                                                                                        arpMatrix.setPairwiseE( curLevel1, index1, curRot1, curLevel2, index2, curRot2, stericE );
											}
										}
									}
								//}
							}
						}
					}
				}
			}
		}
	}
	//Do Bounds Pruning
	public PrunedRotamers<Boolean> DoMinBounds(int numInAS, int strandMut[][],
			double pruningE, PrunedRotamers<Boolean> prunedRotAtRes, float initEw, boolean useSF, boolean boundKS){
		
		return DoMinBounds(numInAS, strandMut, 
					pruningE, prunedRotAtRes, initEw, useSF, boundKS, false); //(onlyBounds==false)
	}
	
	public PrunedRotamers<Boolean> DoMinBounds(int numMutable, int strandMut[][],
			double pruningE, PrunedRotamers<Boolean> prunedRotAtRes, float initEw, boolean useSF, 
			boolean boundKS, boolean onlyBounds){
		
		if (boundKS && onlyBounds){ //cannot be both set at the same time
			System.out.println("ERROR (RotamerSearch, DoMinBounds): boundKS and onlyBounds cannot be simultaneously 'true'.");
			System.exit(1);
		}
		
		MSMinBounds minBoundsRun = new MSMinBounds(arpMatrix,numMutable,strandMut,strandRot,
				pruningE,prunedRotAtRes,splitFlags,useSF,initEw,boundKS,onlyBounds, mutRes2Strand, mutRes2StrandMutIndex, doPerturbations);
		
		if ( (!boundKS) && (!onlyBounds) ) { //use Bounds to prune new rotamers
			eliminatedRotAtRes = minBoundsRun.ComputeEliminatedRotConf();
		}
		else if (boundKS) { //compute Ec
			minBoundsRun.ComputeEliminatedRotConf();
			prunedIsSteric = minBoundsRun.getPrunedSteric();
			Ec_const = minBoundsRun.getEc();
		}
		else { //(onlyBounds==true)
			minBoundsRun.ComputeEliminatedRotConf();
			boundForPartition = minBoundsRun.getBoundForPartition();
		}
		
		minBoundsRun = null;
		
		return eliminatedRotAtRes;
	}
	//Do Bounding Flags
	public void DoBoundFlags(int numMutable, int strandMut[][],
			double pruningE, PrunedRotamers<Boolean> prunedRotAtRes, float initEw, boolean useSF){	
		
		BoundFlags bFlagsRun = new BoundFlags(arpMatrix,numMutable,strandMut,strandRot,
				pruningE,prunedRotAtRes,splitFlags,useSF,initEw, mutRes2Strand, mutRes2StrandMutIndex, doPerturbations);
		

		splitFlags = bFlagsRun.ComputeEliminatedRotConf();
		
		bFlagsRun = null;
	}
	//Do simple Goldstein DEE
	public PrunedRotamers<Boolean> DoDEEGoldstein(int numMutable,
			int strandMut[][], float initEw, PrunedRotamers<Boolean> prunedRotAtRes,
			boolean doMinimize, boolean useSF, boolean minimizeBB, boolean typeDep, 
			boolean useMinDEEPruningEw, float Ival){
		
		//arpmatrix has a row/column for the backbone energies, so we just need
		//the # remaining entries, which is length-1 (only the AS and ligand rotamers)
		DEEGoldstein DEERun = new DEEGoldstein(arpMatrix, arpMatrixMax, numMutable, strandMut, initEw, 
				strandRot, prunedRotAtRes, doMinimize, 
				indIntMinDEE, pairIntMinDEE, splitFlags, useSF, minimizeBB, mutRes2Strand, mutRes2StrandMutIndex,
				typeDep, useMinDEEPruningEw, Ival, doPerturbations);

		eliminatedRotAtRes = DEERun.ComputeEliminatedRotConf();
		
		DEERun = null;
		
		return eliminatedRotAtRes;
	}

        
        public PrunedRotamers<Boolean> DoDEEIndirect(int numMutable, int strandMut[][],
			float initEw, PrunedRotamers<Boolean> prunedRotAtRes, boolean resInPair[], boolean doMinimize,
			boolean magicBullet, boolean distrDEE, boolean minimizeBB, boolean typeDep,
                        boolean doIMinDEE, float Ival){


            DEEIndirect DEERun;

            int numPos = numMutable;

            ArrayList<boolean[]> mpz = new ArrayList<boolean[]>();
            if(doPerturbations)
                mpz = getMinimalPruningZones(numPos, strandMut);

            eliminatedRotAtRes = prunedRotAtRes;

            for(int zoneNum=-1; zoneNum<mpz.size(); zoneNum++){

                boolean inZ[];//Indicates if a residue is in Z

                if(zoneNum == -1){//All active-site residues and the ligand are in the pruning zone
                    inZ = new boolean[numPos];
                    for(int pos=0; pos<numPos; pos++)
                        inZ[pos] = true;
                }
                else
                    inZ = mpz.get(zoneNum);

                if(useTriples)
                    DEERun = new DEEIndirect(arpMatrix, numMutable, strandMut, initEw,
				strandRot, prunedRotAtRes, resInPair, doMinimize,
				splitFlags, magicBullet, distrDEE, minimizeBB, mutRes2Strand, mutRes2StrandMutIndex,
				typeDep, doIMinDEE, Ival, tripleFlags, doPerturbations, inZ);

                else
                    DEERun = new DEEIndirect(arpMatrix, numMutable, strandMut, initEw,
				strandRot, prunedRotAtRes, resInPair, doMinimize,
				splitFlags, magicBullet, distrDEE, minimizeBB, mutRes2Strand, mutRes2StrandMutIndex,
				typeDep, doIMinDEE, Ival, null, doPerturbations, inZ);


                eliminatedRotAtRes = DEERun.ComputeEliminatedRotConf();
                splitFlags = DEERun.getSplitFlags();

                DEERun = null;

            }

            return eliminatedRotAtRes;
	}


        public void DoDEETriples(int numMutable, int strandMut[][],
			float initEw, PrunedRotamers<Boolean> prunedRotAtRes, boolean resInTriple[], boolean doMinimize,
			boolean magicBullet, int magicBulletNum, boolean distrDEE, boolean minimizeBB, boolean typeDep,
                        boolean doIMinDEE, float Ival){

                DEEGoldsteinTriples DEERun = new DEEGoldsteinTriples(arpMatrix, numMutable, strandMut, initEw, strandRot, prunedRotAtRes, 
                        resInTriple, doMinimize, splitFlags, magicBullet, magicBulletNum, distrDEE, minimizeBB,
                        mutRes2Strand, mutRes2StrandMutIndex, typeDep, doIMinDEE, Ival, tripleFlags, doPerturbations );

                if(tripleFlags == null)//At this point tripleFlags == DEERun.tripleFlags
                    DEERun.initializeTripleFlags();

		DEERun.ComputeEliminatedRotConf();
                tripleFlags = DEERun.getTripleFlags();

		DEERun = null;
	}


	//SplitDEE (conformational splitting) with 1 or 2 split positions
	public PrunedRotamers<Boolean> DoDEEConfSplitting(int numMutable, 
			int strandMut[][], float initEw, PrunedRotamers<Boolean> prunedRotAtRes, boolean resInMut[],
			boolean doMinimize, boolean useSF, int numSplits, 
			boolean distrDEE, boolean minimizeBB, boolean typeDep, boolean doIMinDEE, float Ival){
		
		//arpmatrix has a row/column for the backbone energies, so we just need
		//the # remaining entries, which is length-1 (only the AS and ligand rotamers)
		if (numSplits==1){ //1 split position
			DEESplit1f DEERunConfSplitting = new DEESplit1f(arpMatrix, arpMatrixMax, 
					numMutable, strandMut, initEw, 
					strandRot, prunedRotAtRes, resInMut, doMinimize, indIntMinDEE, pairIntMinDEE, 
					splitFlags, useSF, distrDEE, minimizeBB, mutRes2Strand, mutRes2StrandMutIndex, typeDep,
					doIMinDEE, Ival, doPerturbations);
			
			eliminatedRotAtRes = DEERunConfSplitting.ComputeEliminatedRotConf();
			splitFlags = DEERunConfSplitting.getSplitFlags();
			
			DEERunConfSplitting = null;
		}
		else{ // 2 split positions
			DEESplit2f DEERunConfSplitting = new DEESplit2f(arpMatrix, arpMatrixMax, 
					numMutable, strandMut, initEw, 
					strandRot, prunedRotAtRes, resInMut, doMinimize, indIntMinDEE, pairIntMinDEE, 
					splitFlags, useSF, distrDEE, minimizeBB, mutRes2Strand, mutRes2StrandMutIndex, typeDep,
                                        doIMinDEE, Ival, doPerturbations);
			
			eliminatedRotAtRes = DEERunConfSplitting.ComputeEliminatedRotConf();
			splitFlags = DEERunConfSplitting.getSplitFlags();
			
			DEERunConfSplitting = null;
		}
		
		return eliminatedRotAtRes;
	}
	//Do Goldstein DEE pairs
	public void DoDEEPairs(int numMutable, int strandMut[][], 
			float initEw, PrunedRotamers<Boolean> prunedRotAtRes, boolean resInPair[], boolean doMinimize, 
			boolean useSF, boolean magicBullet, boolean distrDEE, boolean minimizeBB, boolean scaleInt, 
			float maxScale, boolean typeDep, boolean doIMinDEE, float Ival){
		
		//arpmatrix has a row/column for the backbone energies, so we just need
		//the # remaining entries, which is length-1 (only the AS and ligand rotamers)
		DEEGoldsteinPairs DEERunPairs = new DEEGoldsteinPairs(arpMatrix, arpMatrixMax, numMutable, strandMut, initEw, strandRot, prunedRotAtRes, resInPair, doMinimize, 
				splitFlags, useSF, magicBullet, distrDEE, minimizeBB, scaleInt, maxScale, 
				mutRes2Strand, mutRes2StrandMutIndex, typeDep, doIMinDEE, Ival, tripleFlags, doPerturbations );
		
		DEERunPairs.ComputeEliminatedRotConf();		
		splitFlags = DEERunPairs.getSplitFlags();
		
		DEERunPairs = null;
	}

        
        public ArrayList<boolean[]> getMinimalPruningZones(int numMutable, int strandMut[][]){//Get the minimal pruning zones for indirect pruning
            //This is for DEEPer
            //answer.get(a)[b] indicates if flexible residue b is in pruning zone a


            ArrayList<boolean[]> mpz = new ArrayList<boolean[]>();

            if(m.perts.length==0)
                return mpz;//Return an empty ArrayList because there are no minimal pruning zones

            int start = 0;
            if(m.perts[0].resAffected.length == numMutable)//If there is an initial perturbation affecting everything (e.g. a full structure switch) don't count it
                start = 1;

            for(int pertNum=start;pertNum<m.perts.length;pertNum++){

                Perturbation pert = m.perts[pertNum];

                int flexResAffected[] = new int[pert.resAffected.length];//Residues affected by the perturbation (numbering with respect to all flexible residues)
                for(int a=0;a<pert.resAffected.length;a++){
                    int molResNum = pert.resAffected[a];//The molecule residue number
                    int strResNum = m.residue[molResNum].strandResidueNumber;
                    int str = m.residue[molResNum].strandNumber;
                    for(int b=0;b<numMutable;b++){
                        int strb = mutRes2Strand[b];

                        //If the b'th flexible residue has the indicated molResNum
                        if( (strb==str) && (strandMut[str][mutRes2StrandMutIndex[b]]==strResNum) )
                            flexResAffected[a] = b;
                    }
                }

                ArrayList<Integer> zones = new ArrayList<Integer>();//Existing pruning zones that residues affected by this perturbation are found in
                for( int zoneNum=0; zoneNum<mpz.size(); zoneNum++){
                    for(int a=0;a<flexResAffected.length;a++){
                        if( mpz.get(zoneNum)[flexResAffected[a]] && (!zones.contains(zoneNum)) )
                            zones.add(zoneNum);
                    }
                }

                if( zones.isEmpty() ){//Need to create a new pruning zone
                    boolean[] new_mpz = new boolean[numMutable];
                    for(int a=0;a<flexResAffected.length;a++)
                        new_mpz[flexResAffected[a]] = true;
                    mpz.add(new_mpz);
                }
                else if ( zones.size() == 1 ){//Put this perturbation's affected residues in an existing pruning zone
                    for(int a=0;a<flexResAffected.length;a++)
                        mpz.get(zones.get(0))[flexResAffected[a]] = true;
                }
                else{//zones.size() > 1: Merge the pruning zones and put this perturbation's affected residues in it
                    boolean[] new_mpz = new boolean[numMutable];
                    for(int a=0;a<zones.size();a++){
                        for(int b=0;b<numMutable;b++){
                            if( mpz.get(zones.get(a))[b] )
                                new_mpz[b] = true;
                        }
                    }
                    for(int a=0;a<flexResAffected.length;a++)
                        new_mpz[flexResAffected[a]] = true;

                    mpz.add(new_mpz);
                    for(int a=zones.size()-1;a>=0;a--)//zones is in ascending order...delete in descending order to avoid messing up the numbering
                        mpz.remove(zones.get(a).intValue());
                }
            }

            return mpz;
        }


	/*
	 * 
	 * End of DEE section
	 * 
	 */	
	
//////////////////////////////////////////////////////////////////////////
//	Compute Min GMEC section
//////////////////////////////////////////////////////////////////////////
	public float doAStarGMEC(String fileName, boolean searchComputeEVEnergy, 
			boolean searchDoMinimization,int numMutable, 
			int strandMut[][], String strandDefault[][],int numMut, float Ew, double bestScore, 
			CommucObj cObj, boolean approxMinGMEC, float lambda, boolean minimizeBB, boolean useEref, 
			float eRef[][], boolean doBackrubs, String backrubFile, boolean useMinDEEPruningEw, float Ival) {
	
		// A rotamer search is performed. For each residue,
		//  every allowable rotamer is tried in combination
		//  with every other allowable rotamer
		// If we're doing a mutation search then residues
		//  are allowed to mutate


                if(doPerturbations)//Make sure minimizer is set properly
                    ((PMinimizer)simpMin).minimizePerturbations = minimizePerturbations;//which should generally be true


		numConfsEvaluated = BigInteger.ZERO;
		computeEVEnergy = searchComputeEVEnergy;
		doMinimization = searchDoMinimization;
		ASAANums = new int[numberMutable];
		curStrRotNum = new int[numberMutable];
		//curASRotNum = new int[numInAS];
		//int curResToASMap[] = new int[m.strand[sysStrNum].numberOfResidues];
			// This map maps the system strand residues back to the AS numbering
			// So something like 8 -> 0, 10 -> 1, 11 -> 2, ...
		//curLigRotNum = 0;
		
		//boolean ligPresent = (ligStrNum>=0); //determine if a ligand is present
		
		for (int str=0;str<numberOfStrands;str++){
			for(int i=0; i<strandMut[str].length;i++){
				m.strand[str].residue[strandMut[str][i]].flexible = true;
			}
		}
		
		/*for (int i=0; i<numInAS; i++) //the AS residues are flexible - this is used by simpMin to set up the minimizer
			m.residue[residueMap[i]].flexible = true;
		if (ligPresent) //the ligand is also flexible
			m.strand[ligStrNum].residue[0].flexible = true;*/
		
	
		if (searchDoMinimization && !searchComputeEVEnergy){
			System.out.println("Warning: In order to do minimization computeEVEnergy must be true");
			return -1;
		}
	
		// Prepare Amber
		if(searchComputeEVEnergy){
			// Amber should already be loaded
			/*if(ligPresent) {
				a96ff.setLigandNum(m.strand[ligStrNum].residue[0].moleculeResidueNumber);
			}*/
			if (doMinimization){
				if (!minimizeBB){ //side-chain minimization
					if (simpMin == null) {
						System.out.println("Warning: Attempting minimization run but simpMin not allocated, RotamerSearch aborting");
						return -1;
					}
					bbMin = null;
					brMin = null;
				}
				else { //backbone minimization
					if (!doBackrubs){
						if (bbMin == null) {
							System.out.println("Warning: Attempting minimization run but bbMin not allocated, RotamerSearch aborting");
							return -1;
						}
						simpMin = null;
						brMin = null;
					}
					else {
						if (brMin == null) {
							System.out.println("Warning: Attempting minimization run but brMin not allocated, RotamerSearch aborting");
							return -1;
						}
						simpMin = null;
						bbMin = null;
					}
				}
			}
		}
	
		// Make sure the allRotamerPairsEnergyName matrices exist
		if (arpMatrix == null) {
			System.out.println("Warning: allRotamerPairsEnergy matrix not loaded");
			return -1;
		}
	
		// Setup the residue number to AS number map
		/*for(int i=0;i<m.strand[sysStrNum].numberOfResidues;i++){
			curResToASMap[i] = -1;
		}
		for(int i=0;i<residueMap.length;i++){
			curResToASMap[residueMap[i]] = i;
		}*/
		
		return doAStarGMECHelper(numMutable, strandMut, fileName, strandDefault, numMut, Ew, bestScore, cObj, 
				approxMinGMEC, lambda, minimizeBB, useEref, eRef, doBackrubs, backrubFile, useMinDEEPruningEw, Ival);
	}
	// Calls AStar repeatedly while the returned conformations still have energy below the threshold;
	//		computes the partial partition function q*;
	// Called by mutationRotamerSearch(.)
	private float doAStarGMECHelper(int numMutable, int strandMut[][], String fileName, 
			String strandDefault[][], int numMut, float Ew, double bestScore, CommucObj cObj, 
			boolean approxMinGMEC, float lambda, boolean minimizeBB, boolean useEref, float eRef[][], 
			boolean doBackrubs, String backrubFile, boolean useMinDEEPruningEw, float Ival){
		
		/*/////////////////////////////////////////////////////////////////////
		PrintStream debugPS = null;
		try {			
			FileOutputStream fileOutputStream = new FileOutputStream("debug.txt");
			BufferedOutputStream bufferedOutputStream = new BufferedOutputStream( fileOutputStream );
			debugPS = new PrintStream( bufferedOutputStream );
		}
		catch (Exception ex) {
			System.out.println("ERROR: An exception occured while opening log file");
		}
		///////////////////////////////////////////////////////////////////////*/
		
		boolean outputFile = (fileName!=null); //output to file
	
		int treeLevels; // total num levels in the conformation tree
		/*if (ligPresent) //determine num tree levels: if ligPresent, then numInAS+1
			treeLevels = numInAS+1;
		else*/
			treeLevels = numMutable;
		
		int numRotForRes[] = new int[treeLevels]; //the number of rotamers for each flexible residue (AS+lig) during a mutation search
		int numRotForResNonPruned[] = new int[treeLevels]; //the number of non-pruned (by MinDEE) rotamers for each flexible residue
		int numTotalRotRed = 0;		//the total number of rotamers for the flexible residues only (during a mutation search)
		int numTotalRotRedNonPruned = 0; //the total num of non-pruned rotamers for the flexible residues
		
		int numTotalConf = 1;
		int curNumRot = 0;
		for (int curLevel=0; curLevel<treeLevels; curLevel++){
			int str = mutRes2Strand[curLevel];
			int strResNum = strandMut[str][mutRes2StrandMutIndex[curLevel]];
			int molResNum = m.strand[str].residue[strResNum].moleculeResidueNumber;
			/*if ((ligPresent)&&(curLevel==(treeLevels-1))){ //the ligand level
				curNumRot = grl.getNumRotForAAtype(ligROT.getIndexOfNthAllowable(0,0));
				if (curNumRot==0) //GLY or ALA
					curNumRot = 1;
			}
			else {*/ //AS residue				
				curNumRot = 0;
				for (int i=0; i<strandRot[str].getNumAllowable(strResNum); i++){ //add the rot for all allowable AA at this residue
					int newRot = getNumRot( str, strResNum, strandRot[str].getIndexOfNthAllowable(strResNum,i) );
					curNumRot += newRot; 
				}
			//}
			numRotForRes[curLevel] = curNumRot;
			numRotForResNonPruned[curLevel] = numRotForRes[curLevel];
			numTotalRotRed += numRotForRes[curLevel];
			
			numTotalConf *= numRotForRes[curLevel];
		}
		
		int numPrunedThisLevel;
		
		for (int curLevel=0; curLevel<treeLevels; curLevel++){ //for each residue
			int str = mutRes2Strand[curLevel];
			int strResNum = strandMut[str][mutRes2StrandMutIndex[curLevel]];
			int molResNum = m.strand[str].residue[strResNum].moleculeResidueNumber;
	
			int curIndex;		System.out.println("curLevel "+curLevel);
			numPrunedThisLevel = 0;
			int numWTrots = 0;
			/*if ((ligPresent)&&(curLevel==(treeLevels-1))){ //the ligand level
				for (int curRot=0; curRot<numRotForRes[curLevel]; curRot++){ //for all rotamers for the given AA
					curIndex = numInAS*numTotalRotamers + curRot;
					if (eliminatedRotAtRes[curIndex]){
						numPrunedThisLevel++;
					}
				}
			}
			else {*/ //AS residue
				for (int curAA=0; curAA<strandRot[str].getNumAllowable(strResNum); curAA++){ //for all allowed AA's
					int index = strandRot[str].getIndexOfNthAllowable(strResNum,curAA);
					int newRot = getNumRot( str, strResNum, index );
					String AAname = strandRot[str].rl.getAAName(index);
					int numPrunedThisAA = 0;
					System.out.print(AAname + "(");
                                        System.out.print(newRot+" Rot): ");
                                        
					for (int curRot=0; curRot<newRot; curRot++){ //for all rotamers for the given AA
						//curIndex = curLevel*numTotalRotamers + rotamerIndexOffset[index]+curRot;
						System.out.print(curLevel+" "+index+" "+curRot+" "+eliminatedRotAtRes.get(curLevel,index,curRot)+" ");
						if (eliminatedRotAtRes.get(curLevel,index,curRot)){
							numPrunedThisLevel++;
							numPrunedThisAA++;
						}
					}System.out.println();
					if(AAname.equals(strandDefault[str][mutRes2StrandMutIndex[curLevel]])){
						numWTrots = newRot - numPrunedThisAA;
				}
			}
			numRotForResNonPruned[curLevel] -= numPrunedThisLevel;
			System.out.println("NumWTRots: "+numWTrots);
		}
		
		int numConfNonPrunedMinDEE = 1;
		for (int curLevel=0; curLevel<treeLevels; curLevel++){
			if (numRotForResNonPruned[curLevel]==0){ //no non-pruned rotamers for curLevel, so no possible conformations
				System.out.println("MinDEE has pruned all possible conformations: position "+curLevel); //at least the minGMEC should be remaining
				if (cObj!=null) {//output num conf info to cObj
					//TODO: confirm program never gets here
					System.out.println("PROGRAM SHOULD NEVER REACH HERE");
					System.exit(0);
					/*cObj.EL_searchNumConfsTotal = numTotalConf;
					cObj.EL_searchNumPrunedMinDEE = numTotalConf;
					cObj.EL_searchNumConfsEvaluated = 0; //no confs evaluated
					cObj.bestBoundEMin = stericE;*/
				}
				
				return -1;
			}
			else {
				numTotalRotRedNonPruned += numRotForResNonPruned[curLevel];
				
				numConfNonPrunedMinDEE *= numRotForResNonPruned[curLevel];
			}
		}
		
		if (cObj!=null) {//output num conf info to cObj
			System.out.println("PROGRAM SHOULD NEVER REACH HERE");
			System.exit(0);
			/*cObj.EL_searchNumConfsTotal = numTotalConf;
			cObj.EL_searchNumPrunedMinDEE = (numTotalConf - numConfNonPrunedMinDEE);*/
		}
		
		System.out.println("ASTAR PRUNING INFO: Total Rots before pruning for each residue: ");
		for(int i=0;i<treeLevels;i++)System.out.print(numRotForRes[i]+" ");System.out.println();
		System.out.println("Total number of rotamers before pruning: "+numTotalRotRed+" ");
		System.out.println("ASTAR PRUNING INFO: Total Rots non-pruned for each residue ");
		for(int i=0;i<treeLevels;i++)System.out.print(numRotForResNonPruned[i]+" ");System.out.println();
		System.out.println("Total number of rotamers after pruning: "+numTotalRotRedNonPruned+" ");

                int AAdefault[] = new int[numMutable];
                int ctr=0;
                for(int str=0;str<numberOfStrands;str++){
                        for (int i=0; i<strandDefault[str].length; i++){
                                AAdefault[ctr] = strandRot[str].rl.getAARotamerIndex(strandDefault[str][i]);
                                ctr++;
                        }
                }


                for (int i=0; i<AAdefault.length; i++)System.out.print(AAdefault[i]+" ");System.out.println();

                //Allocate this so reduceMatrix can fill it in
                boolean[] eliminatedRotAtPosRed = new boolean[numTotalRotRed];//reduced MinDEE matrix

                boolean splitFlagsRed[][] = null, tripleFlagsRed[][][] = null;
                if( useFlagsAStar ){
                    splitFlagsRed = new boolean[numTotalRotRedNonPruned][numTotalRotRedNonPruned];
                    if( useTriples )
                        tripleFlagsRed = new boolean[numTotalRotRedNonPruned][][];//As with the normal tripleFlags, tripleFlagsRed[a][b][c] is only defined for a>b>c
                }



                //Get the reduced min energy matrix
                ReducedEnergyMatrix arpMatrixRed = arpMatrix.reduceMatrix( eliminatedRotAtPosRed,
			numRotForRes, numRotForResNonPruned, treeLevels,
			numTotalRotRedNonPruned, numMutable, strandMut,
                        numTotalRotRed, this, false, splitFlagsRed, tripleFlagsRed);


		
		//Declaring the logPS output here prevents opening an empty file and returning
		//	(for example, if all conformations are pruned by MinDEE above)
		PrintStream logPS = null; //the output file for conf info
		if (outputFile){
			try {			
				FileOutputStream fileOutputStream = new FileOutputStream(fileName,true); //append file if more than 1 partition
				BufferedOutputStream bufferedOutputStream = new BufferedOutputStream( fileOutputStream );
				logPS = new PrintStream( bufferedOutputStream );
			}
			catch (Exception ex) {
				System.out.println("ERROR: An exception occured while opening log file");
			}
		}
		
		
		//Set-up the A* search					
		MSAStar MSAStarSearch = new MSAStar(treeLevels,numRotForResNonPruned,arpMatrixRed,null,
                        splitFlagsRed,tripleFlagsRed,doPerturbations);
		
		curConf = new int[treeLevels]; //the rotamer sequence				
		boolean run1 = true;
		int numConfsOutput = 0;//num confs output to file
		float lowestBound = stericE;
		
		/*if (!doMinimization)
			approxMinGMEC = false; //reset approxMinGMEC, since it is valid only for MinDEE, and not for traditional DEE*/

		updateBestE((float)bestScore);
	
		while (true){
			
			//debugPS.println("curMinE: "+curMinE);debugPS.flush();
			
			if (cObj!=null){				
				/*if (numConfsEvaluated>=cObj.confSeq.length){
					CommucObj.ConfInfo tempConf[] = new CommucObj.ConfInfo[2*cObj.confSeq.length];
					System.arraycopy(cObj.confSeq,0,tempConf,0,cObj.confSeq.length);
					cObj.confSeq = tempConf;
				}
				cObj.confSeq[numConfsEvaluated] = cObj.new ConfInfo(treeLevels);*/
				System.out.println("PROGRAM SHOULD NEVER REACH HERE");
				System.exit(0);
				/*cObj.EL_searchNumConfsEvaluated = numConfsEvaluated.intValue();*/
			}
			
			//clear the values from the previous run
			for(int i=0; i<curAANum.length;i++)
				curAANum[i] = -1;
				//curASRotNum[i] = -1;
			//}
			//curLigAANum = -1;
			//curLigRotNum = -1;
			
			
			curConf = MSAStarSearch.doAStar(run1,numMut,AAdefault,eliminatedRotAtPosRed,strandRot,strandDefault,numRotForRes,strandMut,false,mutRes2Strand, mutRes2StrandMutIndex); //the current rotamer sequence
			
			System.out.println("confNum: "+(numConfsEvaluated.add(BigInteger.ONE)));
			System.out.print("curConf: ");for(int i=0;i<treeLevels;i++)System.out.print(curConf[i]+" ");System.out.println();
			
			//debugPS.println("confNum: "+(numConfsEvaluated+1));
			//debugPS.print("curConf: ");for(int i=0;i<treeLevels;i++)debugPS.print(curConf[i]+" ");debugPS.println();
			//debugPS.flush();
			
			
			
			for (int curRotCheck=0; curRotCheck<treeLevels; curRotCheck++){//check if the conformation is valid
				if (curConf[curRotCheck]==-1){ // no valid conformations remaining
					/*CommucObj.ConfInfo tempConf[] = new CommucObj.ConfInfo[numConfsEvaluated];
					System.arraycopy(cObj.confSeq,0,tempConf,0,numConfsEvaluated);
					cObj.confSeq = tempConf;*/
					
					if (cObj!=null)
						cObj.bestScore = new BigDecimal(getBestE()); //update the best score so far to supply to the next partition
					
					MSAStarSearch = null;
					
					if (outputFile){
						logPS.flush(); //there may still be results to output
					}
					return (getBestE() - lowestBound);
				}
			}
			
			//As the rotamers given to A* are only the non-pruned ones, there is a difference between the
			//	rotamer numbers returned by A* and the actual rotamer numbers for each residue (that is,
			//	A* may return rot 4 for res 3, but rot 3 for res 3 may be pruned, and so the actual number
			//	of the rot to be applied for res 3 is 5)
			int conf[] = new int[curConf.length]; //the conformation with the actual rotamer numbers
			conf = getActualConf(curConf,eliminatedRotAtPosRed,treeLevels,numRotForRes,conf);
			
			System.out.print("actualConf: ");for(int i=0;i<treeLevels;i++)System.out.print(conf[i]+" ");System.out.println();
			
			//debugPS.print("actualConf: ");for(int i=0;i<treeLevels;i++)debugPS.print(conf[i]+" ");debugPS.println();debugPS.flush();
			
			//KER: Need to calculate Types with template before we mutate anything to set the
			//KER: nterm and cterm flags for the residues
			a96ff.calculateTypesWithTemplates();
			//Extract the AA numbers for the current conformation and appply the corresponding AA types
			for (int i=0; i<treeLevels; i++){
				int str = mutRes2Strand[i];
				int strResNum = strandMut[str][mutRes2StrandMutIndex[i]];
				int molResNum = m.strand[str].residue[strResNum].moleculeResidueNumber;
				
				/*if ((ligPresent)&&(i==(treeLevels-1))){ //the ligand level: only one type is allowed
					curLigAANum = ligROT.getIndexOfNthAllowable(0,0);
					if (m.strand[ligStrNum].isProtein)
						ligROT.changeResidueType(m,0,grl.getAAName(curLigAANum),addHydrogens);
				}
				else {*/ // AS residue
					curAANum[molResNum] = getAAIndex(conf[i],i,strandDefault,strandMut);
					//TODO: allow all rotamer libraries the ability to mutate
					if(m.strand[str].isProtein)
						strandRot[str].changeResidueType(m,strResNum,strandRot[str].rl.getAAName(curAANum[molResNum]),addHydrogens,connectResidues);
					
					
				//}
			}
			
			
			
			
			ctr=0;
			for(int str=0;str<numberOfStrands;str++ ){
				for(int i=0; i<strandMut[str].length;i++){
					ASAANums[ctr] = curAANum[m.strand[str].residue[strandMut[str][i]].moleculeResidueNumber];
					ctr++;
				}
			}
			/*for (int i=0; i<ASAANums.length; i++)
				ASAANums[i] = curAANum[residueMap[i]];*/
			
			System.out.print("curAANum: ");for(int i=0;i<numMutable;i++)System.out.print(ASAANums[i]+" ");System.out.println("");
						
			
			//debugPS.print("curAANum: ");for(int i=0;i<numInAS;i++)debugPS.print(ASAANums[i]+" ");debugPS.println(curLigAANum);debugPS.flush();
			
			
			//Extract and apply the rotamers of the current conformation
			int curAS = 0;	
			int curRot = 0;
			for(int str=0;str<numberOfStrands;str++){
				for (int i=0; i<strandMut[str].length; i++){
					//if (curResToASMap[curLevel]!=-1){//make a change only to the AS residues: use the native type for the other residues
					int strResNum = strandMut[str][i];
					int molResNum = m.strand[str].residue[strResNum].moleculeResidueNumber;

                                        if(doPerturbations){
                                                curRot = conf[curAS] - getRotSum(conf[curAS],curAS,strandDefault,strandMut);
                                                //For DEEPer curRot is the the current RC
                                                //Since A* counts all RCs together for all AA types, subtract out the previous AA types' RCs
                                                boolean validRC = ((StrandRCs)strandRot[str]).applyRC(m, strResNum, curRot);
                                                if(!validRC)
                                                    System.err.println("Error: invalid RC " + curRot + " at residue " + strResNum +
                                                            " of strand " + str );

                                                curStrRotNum[curAS] = curRot;
                                        }
					else if (strandRot[str].rl.getNumRotForAAtype(curAANum[molResNum])!=0){//not GLY or ALA
						curRot = conf[curAS] - getRotSum(conf[curAS],curAS,strandDefault,strandMut);
						strandRot[str].applyRotamer(m, strResNum, curRot);
						curStrRotNum[curAS] = curRot;
					}
					else { //GLY or ALA
						curStrRotNum[curAS] = 0;
					}

					curAS++; //prepare the next AS residue
				}
			}		
			/*if (ligPresent){ //apply the ligand rotamer
				if (grl.getNumRotForAAtype(curLigAANum)!=0){//not GLY or ALA
					ligROT.applyRotamer(m, 0, conf[treeLevels-1]);//the ligand level
					curLigRotNum = conf[treeLevels-1];
				}
				else { //GLY or ALA
					curLigRotNum = 0;
				}
			}*/

                        if(doPerturbations){
                            //Output the RC numbers (stored in curStrRotNum) first
                            System.out.print("curRC: ");
                            for(int i=0;i<numMutable;i++)
                                System.out.print(curStrRotNum[i]+" ");//System.out.println(curLigRotNum);


                            //Now output the corresponding sidechain rotamers
                            System.out.println("curSCRot: ");
                            for(int i=0;i<numMutable;i++){
                                int str = mutRes2Strand[i];
                                int strResNum = strandMut[str][mutRes2StrandMutIndex[i]];
                                System.out.print( ((StrandRCs)strandRot[str]).RCRots[strResNum][ASAANums[i]][curStrRotNum[i]] + " " );
                            }
                        }
                        else{
                            System.out.print("curRot: ");
                            for(int i=0;i<numMutable;i++)
                                System.out.print(curStrRotNum[i]+" ");//System.out.println(curLigRotNum);
                        }

			//debugPS.print("curRot: ");for(int i=0;i<numInAS;i++)debugPS.print(curASRotNum[i]+" ");debugPS.println(curLigRotNum);debugPS.flush();
			
			
			
			// After minimization do a m.updateCoordinates() to
			//  resync the actualCoordinates which were changed
			//  in the minimization procedure
			
			//////////////////////////////////////////////////////////////////////////////////////////
			// Energy computation
			float unMinE = 0.0f;
			float minE = 0.0f;
				
			float minELowerBound = computeBestRotEnergyBound(/*numTotalRotamers,rotamerIndexOffset*/);
			
			//debugPS.println("minELowerBound: "+minELowerBound);debugPS.flush();
			
			
			if (run1) //this is the first extracted conformation, so it has the lowest energy lower bound, so store it
				lowestBound = minELowerBound;
			
			
			boolean done = false;
			
			if ((minELowerBound>(getBestE()+Ew)) && (!run1)){ //we already have all confs within Ew of the minGMEC
				done = true;
			}
			else if (approxMinGMEC){ //running the heuristic halting condition
				
				if ((minELowerBound>(lowestBound+lambda)) && (!run1)) //compare the current bound to the lowest bound
					done = true;
			}
			// 2010: This is a debatable issue: Our "new" useMinDEEPruningEw method doesn't prune anything
			//  withing Ew of the lowestBound, while the old MinDEE doesn't prune anything within Ew of
			//   the minGMEC.  In any case, if we pruned too much, returned the new energy window so 
			//   the process can be repeated.  
			if(useMinDEEPruningEw){  
				if(lowestBound+Ival < minELowerBound || getBestE() < minELowerBound){ // We are not done and we pruned too much, repeat search
					return (getBestE()-lowestBound);
				}			  
			}
			
			
			if (done){ //we already have all required conformations
				
				if (cObj!=null)
					cObj.bestScore = new BigDecimal(getBestE()); //update the best score so far to supply to the next partition
				
				MSAStarSearch = null;
				
				/*CommucObj.ConfInfo tempConf[] = new CommucObj.ConfInfo[numConfsEvaluated];				
				System.arraycopy(cObj.confSeq,0,tempConf,0,numConfsEvaluated);				
				cObj.confSeq = tempConf;*/
				
				if (outputFile){
					logPS.flush();
				}
				if(useMinDEEPruningEw && !approxMinGMEC)
					return (getBestE()-lowestBound);
				else
					return 0; //stop the search
			}
			
			else{
				if (computeEVEnergy){
					a96ff.calculateTypesWithTemplates();
					a96ff.initializeCalculation();
					a96ff.setNBEval(hElect,hVDW);
					if (doMinimization){
						if (!minimizeBB){
								/*if(ligPresent){
								simpMin.initialize(m,sysStrNum,ligStrNum,a96ff,sysLR,ligROT,curAANum,curLigAANum,doDihedE,rl,grl);
							}
								else*/
								simpMin.initialize(m,numberOfStrands,a96ff,strandRot,curAANum,doDihedE);
						}
						else { //backbone minimization
							if (!doBackrubs){
									/*if (ligPresent)
									bbMin.initialize(m, a96ff, residueMap, sysStrNum, ligStrNum);
									else*/
									bbMin.initialize(m, a96ff, strandMut, numberOfStrands);
							}
							else {
									brMin.initialize(m, a96ff, strandMut, backrubFile, hSteric, overlapThresh, numberOfStrands,true);
							}
						}
					}
				}
				if (doMinimization){
					if (!minimizeBB) {//side-chain minimization
						unMinE = calcTotalSnapshotEnergy();
						simpMin.minimize(numMinSteps);
						minE = calcTotalSnapshotEnergy();
						if (doDihedE) //add dihedral energies
							minE += simpMin.computeDihedEnergy();
					}
					else {//backbone minimization
						unMinE = calcTotalSnapshotEnergy();
						if (!doBackrubs)
							bbMin.minimizeFull(false);
						else{
							brMin.minimizeFull();
						}
						minE = calcTotalSnapshotEnergy();
					}
                                            m.updateCoordinates();
                                            if(doPerturbations)
                                                m.revertPertParamsToCurState();
				}
				else if (computeEVEnergy){ //no minimization, so traditional DEE
					minE = calcTotalSnapshotEnergy(); //the sum of the precomputed energy terms
					unMinE = minE;
					m.updateCoordinates();
				}

					float totEref = 0.0f;
					float totEntropy = 0.0f;
					if (useEref)
						totEref = getTotSeqEref(eRef,numMutable,strandMut);
					if (EnvironmentVars.useEntropy)
						totEntropy = getTotSeqEntropy(strandMut);
					unMinE -= totEref;
					minE -= totEref;
				
				updateBestE(minE); //for the halting condition
				////////////////////////////////////////////////////////////////////////////////////////////
				
				
				System.out.println(minELowerBound+" "+minE+" "+getBestE());				
				
				//Since we only need to save the information for conformations whose energy is within Ew of
				//	the minGMEC, we can compare the minimized energy for the current conformation to the
				//	lowest minimized energy in the search so far and output only the conformations that
				//	are within Ew. This optimization is important, since writing to the output file is
				//	relatively very exomputationally expensive; moreover, the output file can become very
				//	big, so the optimization reduces the storage requirement. This approach is most beneficial
				//	when low minimized energies are returned early in thse serach, so there will be only a small
				//	number of extra output conformations.
				if (outputFile){//Output to file
					
					if ((approxMinGMEC)||(minE<=(getBestE()+Ew))){ //heuristic stopping condition or minE within Ew of the current lowest energy
						
						numConfsOutput++;
						logPS.print(numConfsOutput+" ");
						for (int i=0; i<treeLevels; i++){
								int str = mutRes2Strand[i];
								int strResNum = strandMut[str][mutRes2StrandMutIndex[i]];
								int molResNum = m.strand[str].residue[strResNum].moleculeResidueNumber;
								/*if ((ligPresent)&&(i==(treeLevels-1)))
								logPS.print(grl.getAAName(curLigAANum)+" ");
								else*/
								logPS.print(strandRot[str].rl.getAAName(curAANum[molResNum])+" ");
						}
						for (int i=0; i<treeLevels; i++){
								/*if ((ligPresent)&&(i==(treeLevels-1)))
								logPS.print(curLigRotNum+" ");
								else*/
									logPS.print(curStrRotNum[i]+" ");
						}


                                                if(doPerturbations){
                                                    //Output the sidechain rotamers in DEEPer
                                                    //(we just outputted the RCs if running DEEPer)
                                                    logPS.print("SC rotamers: ");

                                                    for (int i=0; i<treeLevels; i++){
                                                        int str = mutRes2Strand[i];
                                                        int strResNum = strandMut[str][mutRes2StrandMutIndex[i]];
                                                        int molResNum = m.strand[str].residue[strResNum].moleculeResidueNumber;

                                                        logPS.print( ((StrandRCs)strandRot[str]).RCRots[strResNum][curAANum[molResNum]][curStrRotNum[i]] + " " );
                                                    }
                                                }



						//logPS.println();
						logPS.print("unMinE: "+unMinE+" ");
						logPS.print("minE: "+minE+" ");
						if (doMinimization)
							logPS.print("minBound: "+minELowerBound+" ");
						logPS.print("bestE: "+getBestE());
						logPS.println();
						//if (numConfsOutput%100==0) //flush only every 100 confs, since output is relatively *very* computationally expensive
							logPS.flush();
					}
				}
			}
			
			numConfsEvaluated = numConfsEvaluated.add(BigInteger.ONE);
			
			run1 = false;
		}
	}
	
	
	//Returns the reference energy for the current amino acid sequence assignment (for the mutatable positions only)
	private float getTotSeqEntropy(int [][] strandMut){
		
		double[] entropy = EnvironmentVars.getEntropyTerms();
		
		float totEref = 0.0f;
		
		
		for(int str=0; str<numberOfStrands;str++){
			if(m.strand[str].isProtein){
				for (int i=0; i<strandMut[str].length; i++){
					int strResNum = strandMut[str][i];
					int molResNum = m.strand[str].residue[strResNum].moleculeResidueNumber;
					totEref += entropy[curAANum[molResNum]];
				}
			}
		}
			
		
		return totEref;
	}

	//Returns the reference energy for the current amino acid sequence assignment (for the mutatable positions only)
	public float getTotSeqEref(float eRef[][], int numMutable, int strandMut[][]){
		
		float totEref = 0.0f;
				
		int ctr = 0;
		for(int str=0; str<numberOfStrands;str++){
			for (int i=0; i<strandMut[str].length; i++){
				int strResNum = strandMut[str][i];
				int molResNum = m.strand[str].residue[strResNum].moleculeResidueNumber;					
				totEref += eRef[ctr][strandRot[str].rl.getAARotamerIndex(m.residue[molResNum].name)];//curAANum[molResNum]];
				ctr++;
			}
		}
		/*if (ligPresent)
			totEref += eRef[curLigAANum];*/
		
		return totEref;
	}
	
	private int getAAIndex(int rotIndex, int curRes, String strandDefault[][], int strandMut[][]){
		
		int rotSum = 0;
		
		int str = mutRes2Strand[curRes];
		int strResNum = strandMut[str][mutRes2StrandMutIndex[curRes]];
		int molResNum = m.strand[str].residue[strResNum].moleculeResidueNumber;
		
		for (int i=0; i<strandRot[str].getNumAllowable(strResNum); i++){
                        int curAA = strandRot[str].getIndexOfNthAllowable(strResNum,i);
			int curRot = getNumRot( str, strResNum, curAA );
                        
			rotSum += curRot;
			if (rotSum>rotIndex)
				return curAA;
		}
		return -1;
	}
	
	private int getRotSum(int rotIndex, int curRes, String strandDefault[][], int strandMut[][]){
		
		int rotSum = 0;
		int str = mutRes2Strand[curRes];
		int strResNum = strandMut[str][mutRes2StrandMutIndex[curRes]];
		for (int i=0; i<strandRot[str].getNumAllowable(strResNum); i++){
			int curRot = getNumRot( str, strResNum, strandRot[str].getIndexOfNthAllowable(strResNum,i) );

			if ((rotSum+curRot)>rotIndex)
				return rotSum;
			else
				rotSum += curRot;
		}
		return -1;
	}


        //find how many rotamers or RCs are allowed for the given AA type index at the given
        //residue of the given strand
        // giving 1 for gly or ala
        public int getNumRot(int str, int strResNum, int AAind){

            if(doPerturbations){
                return ((StrandRCs)strandRot[str]).getNumRCs(strResNum, AAind);
            }
            else{
                int ans = strandRot[str].rl.getNumRotForAAtype(AAind);
                if (ans==0)	//ala or gly
                    return 1;
                else
                    return ans;
            }
        }

//////////////////////////////////////////////////////////////////////////
//	End of Compute Min GMEC section
//////////////////////////////////////////////////////////////////////////
	
	
//////////////////////////////////////////////////////////////////////////
//	Generate Backbones section
//////////////////////////////////////////////////////////////////////////
	public void doGenBackbones(String runName, int numMutable, int strandMut[][], double theta, double alpha,
			int numSamples, boolean systematicSampling){
	
		//setup the log file to store information about the generated backbones
		PrintStream logPS = null;
		try {			
			FileOutputStream fileOutputStream = new FileOutputStream((runName+".bb.all"));
			BufferedOutputStream bufferedOutputStream = new BufferedOutputStream( fileOutputStream );
			logPS = new PrintStream( bufferedOutputStream );
		}
		catch (Exception ex) {
			System.out.println("ERROR: An exception occured while opening log file");
			System.exit(1);
		}
		
		Backbone bb = new Backbone();
		
		//the initial phi and psi angles for the residues with flexible backbones
		double initFi[] = new double[numMutable];
		double initPsi[] = new double[numMutable];
		
		for (int i=0; i<numMutable; i++){
			int str = mutRes2Strand[i];
			int strResNum = strandMut[str][mutRes2StrandMutIndex[i]];
			int molResNum = m.strand[str].residue[strResNum].moleculeResidueNumber;
			
			initFi[i] = bb.getFiPsi(m, str, strResNum, 0);
			initPsi[i] = bb.getFiPsi(m, str, strResNum, 1);
			
			System.out.println("AS residue: "+i+" ("+m.strand[str].residue[strResNum].getResNumber()+") phi: "+initFi[i]+" psi: "+initPsi[i]);
			
			if ((initFi[i]==0.0)||(initPsi[i]==0.0)){
				System.out.println("ERROR: Strand "+str+" Residue "+strResNum+" does not have a valid (phi,psi) pair.");
				System.exit(1);
			}
		}
		
		if (systematicSampling){ //systematic sampling (generates a large number of backbones)
			//the number of steps (from -theta to +theta)
			int numSteps = 1;
			if (alpha!=0.0)
				numSteps = 2*(int)(theta/alpha) + 1; //count the initial phi/psi values as well
			
			doGenBackbonesHelper(runName, numMutable, strandMut, theta, alpha, numSteps, bb, initFi, initPsi, logPS, 0);
		}
		else { //random sampling (generates a small number of backbones)
			outputBB(logPS, numMutable, runName, strandMut, bb); //output the initial backbone first
			doGenBackbonesRandHelper(runName, numMutable, strandMut, theta, numSamples, bb, logPS);
		}
	}

	//Generates up to numSamples backbones by applying random phi/psi changes within theta degrees 
	//		This version is used to generate a small number of backbones (the systematic version below
	//		will generate a very large number of backbones, even for a very small number of steps)
	private void doGenBackbonesRandHelper(String runName, int numMutable, int strandMut[][], double theta,
			int numSamples, Backbone bb, PrintStream logPS){
		
		Random randNum = new Random();
		
		for (int curSample=0; curSample<numSamples; curSample++){ //generate up to numSamples backbones
			
			float curFiChange[] = new float[numMutable]; //the phi changes for each residue
			float curPsiChange[] = new float[numMutable]; //the psi changes for each residue
			
			for(int curRes=0; curRes<numMutable; curRes++){ //apply the random (phi,psi) changes for each residue
				int str = mutRes2Strand[curRes];
				int strResNum = strandMut[str][mutRes2StrandMutIndex[curRes]];
				int molResNum = m.strand[str].residue[strResNum].moleculeResidueNumber;
				//Get the random (phi,psi) change
				curFiChange[curRes] = 2*(randNum.nextFloat()-0.5f)*(float)theta; //phi change within (-theta,theta)
				curPsiChange[curRes] = 2*(randNum.nextFloat()-0.5f)*(float)theta; //psi change within (-theta,theta)
				
				//Apply the (phi,psi) change
				bb.applyFiPsi(m,str,strResNum,curFiChange[curRes],0);
				bb.applyFiPsi(m,str,strResNum,curPsiChange[curRes],1);
			}
			
			//we have a full backbone conformation, so output
			outputBB(logPS, numMutable, runName, strandMut, bb);
			
			for(int curRes=0; curRes<numMutable; curRes++){ //restore the original (phi,psi) for each residue
				int str = mutRes2Strand[curRes];
				int strResNum = strandMut[str][mutRes2StrandMutIndex[curRes]];
				bb.applyFiPsi(m,str,strResNum,-curFiChange[curRes],0);
				bb.applyFiPsi(m,str,strResNum,-curPsiChange[curRes],1);
			}
		}
	}

	//Recursively generates backbones by applying phi/psi changes within theta degrees 
	//		at a step of alpha degrees to each residue with a flexible backbone (systematic sampling)
	private void doGenBackbonesHelper(String runName, int numMutable, int strandMut[][], double theta, double alpha, 
			int numSteps, Backbone bb, double initFi[], double initPsi[], PrintStream logPS, int curRes){
		
		if (curRes==numMutable){//we have a full backbone conformation, so output
			outputBB(logPS, numMutable, runName, strandMut, bb);
		}
		else { //we only have a partial backbone, so apply changes to the current residue
			
			if (curRes==0) 
				System.out.println("Starting..");
			
			//First, move the phi/psi angle to -theta, then apply the changes to +theta, at alpha steps
			
			int str = mutRes2Strand[curRes];
			int strResNum = strandMut[str][mutRes2StrandMutIndex[curRes]];
			
			//move phi to -(theta+alpha), so that the first move in the *for* statement below goes to -theta
			bb.applyFiPsi(m,str,strResNum,-(theta+alpha),0);
			
			//apply phi changes up to +theta
			for (int curStepFi=0; curStepFi<numSteps; curStepFi++){ //apply each alpha step up to a displacement of +theta from the initial phi
				
				if (curRes==0)
					System.out.print("*");
				
				bb.applyFiPsi(m,str,strResNum,alpha,0);
				
				//move psi to -(theta+alpha), so that the first move in the *for* statement below goes to -theta
				bb.applyFiPsi(m,str,strResNum,-(theta+alpha),1);
				
				for (int curStepPsi=0; curStepPsi<numSteps; curStepPsi++){ //apply each alpha step up to a displacement of +theta from the initial psi
				
					if (curRes==0)
						System.out.print(".");
					
					bb.applyFiPsi(m,str,strResNum,alpha,1);
					if (checkStericsBBonly(str,strResNum)) {//passed steric test, so move to the next residue
						doGenBackbonesHelper(runName,numMutable,strandMut,theta,alpha,numSteps,bb,initFi,initPsi,logPS,curRes+1);
					}
				}
				
				//restore initial psi
				bb.applyFiPsi(m,str,strResNum,-theta,1);
				
				if (curRes==0)
					System.out.println();
			}
			
			//restore initial phi
			bb.applyFiPsi(m,str,strResNum,-theta,0);
			
			if (curRes==0)
				System.out.println("done");
		}
	}
	
	//Checks if the current backbone is sterically allowed and outputs the pdb and log information
	private void outputBB(PrintStream logPS, int numMutable, String runName, int strandMut[][], Backbone bb){
		
		//Check all residues for sterics against the whole strand since all residues
		//		are already assigned (up to this point, we have checked for sterics only
		//		against the residues up to a given AS residue);
		//The backbone movement may have actually caused unallowed sterics between rigid
		//		parts of the molecule, so we check for this possibility
		boolean stericAllowed = true;
		for(int str=0; str<numberOfStrands;str++){
			for (int i=0; i<m.strand[str].numberOfResidues; i++){
				if (!checkAllStericsBBonly(str,i)){
				stericAllowed = false;
				break;
				}
			}
		}
		
		if (stericAllowed){ //all sterics are allowed
		
			String fileName = (runName+System.currentTimeMillis()+".pdb");
			
			saveMolecule(m,fileName,0.0f);//save the molecule
			
			//output the molecule file name and all (phi,psi) pairs for the residues with flexible backbones
			logPS.print(fileName+" ");
			for(int str=0;str<numberOfStrands;str++){
				for (int i=0; i<strandMut[str].length; i++){
					logPS.print("( "+bb.getFiPsi(m, str, strandMut[str][i], 0)+" , "+bb.getFiPsi(m, str, strandMut[str][i], 1)+" ) ");
			}
			}
			logPS.println();
			logPS.flush();
		}
	}
	// This function checks the sterics of the given conformation;
	//  it checks backbone atoms of the given residue resNum against all
	//  backbone atoms that are in residues 0..resNum-1 for the given strand only
	//  If any two atoms overlap by more than overlapThresh then
	//  false is returned
	// strandNum is the number of the strand containing resNum
	// Hydrogens are NOT used in checking sterics in this function
	private boolean checkStericsBBonly(int strandNum, int resNum) {
	
		Residue res = m.strand[strandNum].residue[resNum];
	
		Atom tmpAtm = null;
		int resToCheck = 0;
		
		for(int i=0;i<res.numberOfAtoms;i++) {
			
			if (isBBatom(res.atom[i])){ //backbone atom
				
				resToCheck = resNum;
				for(int w=0;w<resToCheck;w++) {
					for(int t=0;t<m.strand[strandNum].residue[w].numberOfAtoms;t++) {
						tmpAtm = m.strand[strandNum].residue[w].atom[t];
						
						if (isBBatom(tmpAtm)){
							if (!(tmpAtm.elementType.equalsIgnoreCase("H"))) {
								if ((res.atom[i].distance(tmpAtm) < ((tmpAtm.radius + res.atom[i].radius)/100.0) - softOverlapThresh)){
									if (!(res.atom[i].bondedTo(tmpAtm.moleculeAtomNumber))) {
										return false;
									}
								}
							}
						}
					}
				}	
			}
		}
	
		return true; //if we are here, then everything passed the steric test
	}
	
	//This function is similar to checkStericsBBonly(), but it checks all residues in strandNum against resNum, 
	//		instead of just checking the residues up to resNum
	private boolean checkAllStericsBBonly(int strandNum, int resNum) {
		
		Residue res = m.strand[strandNum].residue[resNum];
	
		Atom tmpAtm = null;
		int resToCheck = 0;
		
		for(int i=0;i<res.numberOfAtoms;i++) {
			
			if (isBBatom(res.atom[i])){ //backbone atom
				
				resToCheck = m.strand[strandNum].numberOfResidues;
				for(int w=0;w<resToCheck;w++) {
					if (w!=resNum){
						for(int t=0;t<m.strand[strandNum].residue[w].numberOfAtoms;t++) {
							tmpAtm = m.strand[strandNum].residue[w].atom[t];
							
							if (isBBatom(tmpAtm)){
								if (!(tmpAtm.elementType.equalsIgnoreCase("H"))) {
									if ((res.atom[i].distance(tmpAtm) < ((tmpAtm.radius + res.atom[i].radius)/100.0) - softOverlapThresh)){
										if (!(res.atom[i].bondedTo(tmpAtm.moleculeAtomNumber))) {
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
	
		return true; //if we are here, then everything passed the steric test
	}
	
	//Determines if the given atom is a backbone atom
	private boolean isBBatom(Atom at){
		
		return ((at.name.equalsIgnoreCase("N"))||(at.name.equalsIgnoreCase("CA"))
				||(at.name.equalsIgnoreCase("C"))||(at.name.equalsIgnoreCase("O")));
	}
//////////////////////////////////////////////////////////////////////////
//	End of Generate Backbones section
//////////////////////////////////////////////////////////////////////////
	
	
//////////////////////////////////////////////////////////////////////////
	
	/*//Determines how many residue positions in the system strand (pos is strand-relative numbering)
	//		are within dist from residue position pos
	public boolean [] getProxAS(int pos, float dist, boolean as[]){
		
		Residue r1 = m.strand[sysStrNum].residue[pos];
		for (int i=0; i<m.strand[sysStrNum].numberOfResidues; i++){
			if (i!=pos){
				Residue r2 = m.strand[sysStrNum].residue[i];
				if (r1.getDistSC(r2)<=dist)
					as[i] = true;
				else
					as[i] = false;
			}
		}
		
		return as;
	}*/
	
	
//////////////////////////////////////////////////////////////////////////////////
	public void saveMolecule(Molecule m, String fname, float energy){
		m.backupAtomCoord();
		boolean printSegID = false;
		try{
			FileOutputStream fileOutputStream = new FileOutputStream(fname);
			BufferedOutputStream bufferedOutputStream = new BufferedOutputStream( 
					fileOutputStream);
			PrintStream printStream = new PrintStream(bufferedOutputStream);
			Hashtable params = new Hashtable(7);
			params.put("printSegID",new Boolean(printSegID));
			params.put("comment","");
			params.put("energy", energy);
			params.put("showConnect",new Boolean(false));					
			new SaveMolecule(m, printStream, params); 
			printStream.close();
		}
		catch (IOException e) {
			System.out.println("ERROR: An io exception occurred while writing file");
			System.exit(0);
		}
		catch ( Exception e ){
			System.out.println(e.toString());
			System.out.println("ERROR: An exception occurred while writing file");
			System.exit(0);
		}
		m.restoreAtomCoord();
	}

	//Computes and stores the pairwise energies for the current molecule; this can be used to compare
	//		the actual energies for each pair of residues (computed here) against the respective lower bounds computed
	//		during the pairwise energy minimization
	/*private void pairwiseEnergySidechainMolecule (int residueMap[]){
		
		boolean savedEnergyEvalSC[] = new boolean[m.numberOfResidues];
		boolean savedEnergyEvalBB[] = new boolean[m.numberOfResidues];
		boolean savedFlexible[] = new boolean[m.numberOfResidues];
			
		// Save the energy eval flag, clear them at the same time
		for(int i=0;i<m.numberOfResidues;i++){
			savedEnergyEvalSC[i] = m.residue[i].getEnergyEvalSC();
			savedEnergyEvalBB[i] = m.residue[i].getEnergyEvalBB();
			savedFlexible[i] = m.residue[i].flexible;
		}
		
		PrintStream logPS = setupOutputFile("energies.out");
		
		for(int i=0;i<m.numberOfResidues;i++){
			m.residue[i].setEnergyEval(false, false);
			m.residue[i].flexible = false;
		}
		
		for (int i=0; i<residueMap.length; i++){
			
			m.strand[sysStrNum].residue[residueMap[i]].setEnergyEval(true, true);
			m.strand[sysStrNum].residue[residueMap[i]].flexible = true;
			
			float intraEi = computeEnergyOfOnlyRes(m.strand[sysStrNum].residue[residueMap[i]].moleculeResidueNumber);
			
			for (int j=i+1; j<residueMap.length; j++){
				
				m.strand[sysStrNum].residue[residueMap[j]].setEnergyEval(true, true);
				m.strand[sysStrNum].residue[residueMap[j]].flexible = true;
				
				float intraEj = computeEnergyOfOnlyRes(m.strand[sysStrNum].residue[residueMap[j]].moleculeResidueNumber);
				
				float pairE = calcTotalSnapshotEnergy() - intraEi - intraEj;
				
				logPS.println(i+" "+j+" "+pairE+" "+intraEj);
				
				m.strand[sysStrNum].residue[residueMap[j]].setEnergyEval(false, false);
				m.strand[sysStrNum].residue[residueMap[j]].flexible = false;
			}
			
			for(int j=0;j<m.numberOfResidues;j++){
				m.residue[j].setEnergyEval(true, true);
				m.residue[j].flexible = false;
			}
			if (ligStrNum != -1) {
				for(int j=0;j<m.strand[ligStrNum].numberOfResidues;j++)
					m.strand[ligStrNum].residue[j].setEnergyEval(false, false);
			}
			for(int j=0;j<residueMap.length;j++)
				m.strand[sysStrNum].residue[residueMap[j]].setEnergyEval(false, false);
			m.strand[sysStrNum].residue[residueMap[i]].setEnergyEval(true, true);
			m.strand[sysStrNum].residue[residueMap[i]].flexible = true;
			
			float shellE = computeEnergyOfOnlyTemplate(residueMap);
			
			float resShellE = calcTotalSnapshotEnergy() - intraEi - shellE;
			
			for(int j=0;j<m.numberOfResidues;j++){
				m.residue[j].setEnergyEval(false, false);
				m.residue[j].flexible = false;
			}
			
			logPS.println(i+" "+intraEi+" "+resShellE+" "+shellE);
			logPS.flush();
		}
		logPS.close();
		
		// Restore the energy eval and flexibility flags
		for(int i=0;i<m.numberOfResidues;i++){
			m.residue[i].setEnergyEval(savedEnergyEvalSC[i], savedEnergyEvalBB[i]);
			m.residue[i].flexible = savedFlexible[i];
		}
	}*/
	
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
	
	
	public String printBigNum(BigDecimal bd, int sigDigits){
        int numDigits = (bd.toPlainString()).indexOf('.');
        if(numDigits == -1)
                numDigits = (bd.toPlainString().length());

        int newScale = -1*(numDigits - sigDigits);

        return bd.setScale(newScale,BigDecimal.ROUND_HALF_UP).toString();
	}

	public void resetMatrices() {
		arpMatrix = null;
		arpMatrixMax = null;
	}

	public static float[][] loadErefMatrix(String fileName){
		float[][] eref = null;
		try{
			ObjectInputStream in = new ObjectInputStream(new FileInputStream(fileName));
			eref = (float [][])in.readObject();
			in.close();
		}
		catch (Exception e){}
		return eref;
	}

	public void addEntropyTerm(boolean doMinimize, int[][] strandMut){
		
		double[] entropy = EnvironmentVars.getEntropyTerms();

		//NumMutRes (Don't count non-protein strands)
		int numMutRes = 0;
		for(int i=0;i<strandMut.length;i++){
			numMutRes += strandMut[i].length;
		}
		
		int ind = 1; //skip the entry [0][0], since this is the fixed template energy
		for (int i=0; i<numMutRes; i++){
			int str = mutRes2Strand[i];
			int strResNum = strandMut[str][mutRes2StrandMutIndex[i]];
			if(!m.strand[str].isProtein){
				for (int j=0; j<strandRot[str].getNumAllowable(strResNum); j++){
					int aaInd = strandRot[str].getIndexOfNthAllowable(strResNum,j);
					int numRot = getNumRot(str, strResNum, aaInd);
					for (int k=0; k<numRot; k++){
						arpMatrix.addToIntraE( i, aaInd, k, (float) entropy[aaInd] );
						if (doMinimize)
                                                        arpMatrixMax.addToIntraE( i, aaInd, k, (float) entropy[aaInd] );
						ind++;
					}
				}
			}
			else{
				System.out.println("Entropy term not added for non-amino acid residue: "+m.strand[str].residue[strResNum].fullName);
			}
		}
			
	}
	
	private void saveConf(int[] conf, float energy, String filename, int[][] strandMut,
			boolean minimizeBB, boolean doBackrubs) {
		// TODO Auto-generated method stub
		applyRotamers(strandMut, conf);
                m.backupAtomCoord();//The backup coordinates will now be in the unminimized state so we can return to it

		if(doMinimization){
			if(!minimizeBB)
				simpMin.minimize(numMinSteps);
			else{
				if(!doBackrubs)
					bbMin.minimizeFull(false);
				else{
					System.out.println("BACKRUBS NOT IMPLEMENTED YET");
					//brMin.minimizeFull();
				}
			}
		}
		m.saveMolecule(filename, energy);
                if(doPerturbations)
                    m.revertPertParamsToCurState();
                m.restoreAtomCoord();
                m.updateCoordinates();
                
	}


        //Check parametric incompatibility of a pair of RCs, with their residues and AA types also given
        public boolean isParametricallyIncompatible(Residue firstRes, int curAA1, int RC1, Residue secondRes, int curAA2, int RC2){

            //Get the residues' residue perturbations states
            int pertState1 = ((StrandRCs)strandRot[firstRes.strandNumber]).RCPertStates[firstRes.strandResidueNumber][curAA1][RC1];
            int pertState2 = ((StrandRCs)strandRot[secondRes.strandNumber]).RCPertStates[secondRes.strandResidueNumber][curAA2][RC2];
            //Now see if the residue perturbation states are compatible
                    
            for(int p1=0;p1<firstRes.perts.length;p1++){//Check all the perturbations of firstRes in this RC for compatibility
                for(int p2=0;p2<secondRes.perts.length;p2++){
                    if(firstRes.perts[p1] == secondRes.perts[p2]){//The residues share a perturbation
                        if(firstRes.pertStates[pertState1][p1] != secondRes.pertStates[pertState2][p2])//The states of this perturbation don't match for the given RCs
                            return true;//So the RCs are parametrically incompatible
                    }
                }
            }

            return false;//No incompatibility found, so we're good
        }
	
}//end of RotamerSearch class
