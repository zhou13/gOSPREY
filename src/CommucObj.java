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

////////////////////////////////////////////////////////////////////////////////////////////////
// CommucObj.java
//
//  Version:           2.1 beta
//
//
// authors:
//    initials    name                 organization                email
//   ---------   -----------------    ------------------------    ----------------------------
//     RHL        Ryan Lilien          Dartmouth College           ryan.lilien@dartmouth.edu
//	   ISG		  Ivelin Georgiev	   Duke University			   ivelin.georgiev@duke.edu
//     KER        Kyle E. Roberts       Duke University         ker17@duke.edu
//     PGC        Pablo Gainza C.       Duke University         pablo.gainza@duke.edu
//     MAH        Mark A. Hallen	Duke University         mah43@duke.edu
///////////////////////////////////////////////////////////////////////////////////////////////

/**
 *
 * Written by Ryan Lilien (2002-2004) and Ivelin Georgiev (2004-2009)
 *
 */

import java.io.Serializable;
import java.util.Hashtable;
import java.util.Vector;
import java.math.BigDecimal;
import java.math.BigInteger;

/**
 * The CommucObj class is a data structure used in communication between the master
 *  and slave nodes. It is basically just a data container.
 * It allows the master node to specify what type of search the slave should perform
 *  and it allows the slave to return the result of the computation to the master.
 */
public class CommucObj implements Serializable
{
	public class ConfInfo implements Serializable {
		String AA[] = null;
		int AAnums[] = null;
		int rot[] = null;
		float minBound = 0.0f;
		float unMinE = 0.0f;
		float minE = 0.0f;
		
		ConfInfo (int treeLevels){
			rot = new int[treeLevels];
			AA = new String[treeLevels];
			AAnums = new int[treeLevels];
		}
	}
	
	ConfInfo confSeq[] = null;
	/*int E_searchNumConfsTotal = 0;
	int E_searchNumConfsPrunedByE = 0;
	int E_searchNumConfsPrunedByS = 0;
	int E_searchNumConfsEvaluated = 0;
	int E_searchNumConfsLeft = 0;
	int E_searchNumPrunedMinDEE = 0;
	int EL_searchNumConfsTotal = 0;
	int EL_searchNumConfsPrunedByE = 0;
	int EL_searchNumConfsPrunedByS = 0;
	int EL_searchNumConfsEvaluated = 0;
	int EL_searchNumConfsLeft = 0;
	int EL_searchNumPrunedMinDEE = 0;
	float EL_searchBestEnergyFound = 99999.0f;*/
	int searchNumConfsTotal[] = null;
	int searchNumConfsPrunedByE[] = null;
	int searchNumConfsPrunedByS[] = null;
	int searchNumConfsEvaluated[] = null;
	int searchNumConfsLeft[] = null;
	int searchNumPrunedMinDEE[] = null;
	float searchBestEnergyFound[] = null;
	BigDecimal q[] = null;
	/*BigDecimal q_E = new BigDecimal(0.0);
	BigDecimal q_EL = new BigDecimal(0.0);
	BigDecimal q_L = new BigDecimal(0.0);*/
	BigDecimal bestScore = new BigDecimal(0.0);
	double bestEMin[] = null;
	double bestE[] = null;
	/*double bestBoundEMin = 0.0;		// The best minimized bound energy
	double bestUnBoundEMin = 0.0;	// The best minimized unbound energy
	double bestBoundE = 0.0;		// The best bound energy (no minimization)
	double bestUnBoundE = 0.0;*/		// The best unbound energy (no minimization)
		// passed from the master to the slave
	
boolean typeDep = false;	

	String currentMutation[] = null;
	//int resMap[] = null;
	/*String resDefault[] = null;*/
	int[][] strandMut = null;
	String[][] strandDefault = null;
	boolean[] strandPresent = null;
	String[][] strandLimits = null;
	int strandsPresent = 0;
	int mutableSpots = 0;
	//String ligType = null;
	//boolean ligPresent = false;
	int numMutations = 0;
	String arpFilenameMin = null;
	String arpFilenameMax = null;
	boolean minDEEtypeMS = false;
	int algOption = 0;
	int numSplits = 0;
	String AAallowed[][] = null;
	boolean resMutatable[] = null;
	String minDEEfile = null;
	float initEw = 0.0f;
	float pruningE = (float)Math.pow(10,38);
	boolean repeatEW[] = null;
	boolean allPruned[] = null;
	/*boolean E_repeatEw = false;
	boolean EL_repeatEw = false;
	boolean EL_allPruned = false;
	boolean E_allPruned = false;*/
	int mutationNumber = -1;
	ParamSet params = null;
	boolean PEMcomp = false;
	boolean entropyComp = false; //this *must* be false for the pairwise matrix energy computation
	boolean compASdist = false;
	boolean asDist[] = null;
	float dist = 0.0f;
	int mutRes2Strand[] = null;
	int mutRes2StrandMutIndex[] = null;
	
	double gamma = 0.01;
	float epsilon = 0.03f;
	int numResidues = 0;
	float stericThresh = -10000.0f;
	float softStericThresh = -10000.0f;
	//int numInAS = 0;
	
	boolean computeEVEnergy = true;
	boolean doMinimization = true;
	boolean minimizeBB = false;
	boolean doBackrubs = false;
	String backrubFile = null;
	boolean repeatSearch = true;
	boolean calculateVolumes = true;
	boolean approxMinGMEC = false;
	float lambda = (float)Math.pow(10, 38);
	boolean distDepDielect = true;
	double dielectConst = 1.0;
	boolean doDihedE = false;
	boolean doSolvationE = false;
	double solvScale = 1.0;
	double vdwMult = 1.0;
	boolean scaleInt = false;
	float maxIntScale = 1.0f;
	double stericE = Math.pow(10, 38);
	boolean useEref = false;
	float eRef[][] = null;
	int numberOfStrands = 0;
	// Timing info (in seconds)
	int q_Time[] = null;
	/*int q_E_Time = 0;
	int q_EL_Time = 0;*/
	int numComplexes = 0;

	// Identification information
	int slaveNum = -1;
		// the number of the slave with which this packet communicates 1..n
	int portNum = -1;
		// the port number at which the slave can be reached perhaps 10000...
	String machineName = null;
		// the machine name at which the slave can be reached

	// Stopping information
	boolean workToDo = true;
		// when workToDo is false the slave exits
	
	
	//Variables specific to PEM computation
	int resMut[] = null;
	String flagMutType = null;	
	SamplingEEntries compEE[] = null;//initialized by the slave node, not by the master
	//int numLigRotamers = 0;	
	int elapsedTime = 0; // timing info (in seconds)
	
	//Variables specific to distributed DACS and distributed DEE computations
	PrunedRotamers<Boolean> prunedRot = null;
	boolean useSF = false;
	boolean distrDACS = false;
	boolean distrDEE = false;
	boolean splitFlags[][] = null;
	String rotFileIn = null;
	String sfFileIn = null;
	String sfFileOut = null;
	int numSpPos = -1;
	int msp[] = null;
	int typeDEE = -1;
	int initDepth = -1;
	int subDepth = -1;
	int diffFact = -1;
	double minRatioDiff = 0.0;
	BigInteger numInitUnprunedConf = null;
	String outputPruneInfo = null;
	String outputConfInfo = null;
	Index3 partIndex[] = null;
	boolean KSGMEC = false;
	boolean KSCONFTHRESH = false;
	public int curMut = -1;
	public boolean templateAlwaysOn = false;
	public boolean addOrigRots = false;
	RotamerLibrary rl;
	
	public boolean saveTopConfs = false;
	public boolean printTopConfs = false;
	public int numTopConfs;
	public String pdbName = "out";
	public boolean useMaxKSconfs;
	public BigInteger maxKSconfs;
	public int curStrForMatrix;

        boolean useTriples;
        boolean useFlagsAStar;
        boolean magicBulletTriples;
        int magicBulletNumTriples;
        //Variables specific to DEEPer
        boolean doPerturbations;
        boolean minimizePerts;
        //Also, in this context doMinimization && !minimizePerts implies just minimizing dihedrals
        String pertFile;
        boolean addWTRot;
        boolean idealizeSC;


	CommucObj() {
	}
	
	public void initConfSeq(int numConfs,int treeLevels){
		confSeq = new ConfInfo[numConfs];
		for(int j=0; j<numConfs;j++){
			confSeq[j] = new ConfInfo(treeLevels);
		}
	}
	
	/*public void setPartitionProperties(int runNum, PartitionMessage p){
		this.q[runNum] = p.q;
		this.bestEMin[runNum] = p.bestEMin;
		this.bestE[runNum]=p.bestE;
		this.q_Time[runNum]=p.q_Time;
		this.repeatEW[runNum]=p.repeatEW;
		this.allPruned[runNum]=p.allPruned;
		this.searchNumConfsTotal[runNum]=p.searchNumConfsTotal;
		this.searchNumConfsPrunedByE[runNum]=p.searchNumConfsPrunedByE;
		this.searchNumConfsPrunedByS[runNum]=p.searchNumConfsPrunedByS;
		this.searchNumConfsEvaluated[runNum]=p.searchNumConfsEvaluated;
		this.searchNumConfsLeft[runNum]=p.searchNumConfsLeft;
		this.searchNumPrunedMinDEE[runNum]=p.searchNumPrunedMinDEE;
		this.searchBestEnergyFound[runNum]=p.searchBestEnergyFound;
	}*/

}
