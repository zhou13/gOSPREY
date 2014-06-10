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

////////////////////////////////////////////////////////////////////////////////////////////
// MutationManager.java
//
//  Version:           2.1 beta
//
//
// authors:
//    initials    name                 organization                email
//   ---------   --------------      ------------------------     ------------------------------
//     RHL        Ryan Lilien          Dartmouth College           ryan.lilien@dartmouth.edu
//	   ISG		  Ivelin Georgiev	   Duke University			   ivelin.georgiev@duke.edu
//	  KER        Kyle E. Roberts       Duke University         ker17@duke.edu
//    PGC        Pablo Gainza C.       Duke University         pablo.gainza@duke.edu
//     MAH        Mark A. Hallen	Duke University         mah43@duke.edu
////////////////////////////////////////////////////////////////////////////////////////////

/**
 * Written by Ryan Lilien (2002-2004) and Ivelin Georgiev (2004-2009)
 *
 */

import java.io.*;
import java.nio.channels.FileChannel;
import java.nio.channels.FileLock;
import java.util.Hashtable;
import java.util.Iterator;
import java.math.*;

import mpi.*;

/**
 * The MutationManager class maintains a list of mutations to be tested, maintains
 *  their scores, prints a log file, and generally manages the mutations to test.
 */
public class MutationManager {

    //the algorithm options that define what pruning criteria will be applied
    //	NOTE!!! These must be the same as in KSParser.java
    final int optSimple = 1;
    final int optBounds = 2;
    final int optSplit = 3;
    final int optPairs = 4;

    final double constRT = 1.9891/1000.0 * 298.15;   // in kCal / kelvin-mole (the T value here should be consistent with T in EEF1)

    // Information needed by all mutations
    CommucObj cObjArray[] = null;
    //int residueMap[] = null;
    //String resDefault[] = null;
    int[][] strandMut = null;
    String[][] strandDefault = null;
    boolean[] strandPresent = null;
    String[][] strandLimits = null;
    int strandsPresent = 0;
    boolean typeDep = false;
    //String ligType = null;
    boolean ligPresent = false;
    int numMutations = 0;
    String arpFilenameMin = null;
    String arpFilenameMax = null;
    boolean resMutatable[][] = null;
    int algOption = 0;
    int numSplits = 0;
    String AAallowed[][] = null;
    String minDEEfile = null;
    float initEw = 0.0f;
    float pruningE = (float)Math.pow(10,38);
    double gamma = 0.01;  // The gamma used in inter-mutation pruning
    float epsilon = 0.03f;  // The epsilon used in intra-mutation pruning
    float stericThresh = -10000.0f;
    float softStericThresh = -10000.0f;
    //int numInAS = 0;
    int mutableSpots = 0;
    boolean computeEVEnergy = true;
    boolean doMinimization = true;
    boolean minimizeBB = false;
    boolean doBackrubs = false;
    String backrubFile = null;
    boolean repeatSearch = true;
    boolean calculateVolumes = true;
    BigDecimal bestScore = new BigDecimal("0.0");
    BigDecimal q_L = new BigDecimal("0.0");
    ParamSet sParams = null;
    boolean approxMinGMEC = false;
    float lambda = (float)Math.pow(10,38);
    double stericE = Math.pow(10,38);
    boolean distDepDielect = true;
    double dielectConst = 1.0;
    boolean doDihedE = false;
    boolean doSolvationE = false;
    double solvScale = 1.0;
    double vdwMult = 1.0;
    boolean scaleInt = false;
    float maxIntScale = 1.0f;
    boolean useEref = false;
    float[][] eRef = null;
    boolean entropyComp = false; //this *must* be false for the pairwise matrix energy computation
    float asasE[][][][] = null;
    boolean compASdist = false;
    boolean asDist[][] = null;
    float dist = 0.0f;
    int numberOfStrands = 0;
    PrintStream logPS = null;
    OneMutation mutArray[] = null;
    int mutRes2Strand[] = null;
    int mutRes2StrandMutIndex[] = null;

    private boolean saveTopConfs;
    private boolean printTopConfs;
    private int numTopConfs;

    //Variables specific to PEM computation
    PairwiseEnergyMatrix pairEMatrixMin = null;
    PairwiseEnergyMatrix pairEMatrixMax = null;
    //float[][] eRefMatrix = null;
    float curMaxE = -(float)Math.pow(10,30);
    //int numLigRotamers = 0;

    float pairEMatrixMinEntropy[][] = null;
    //KER: This is only used for the SCMF entropy run.
    RotamerLibrary rotLib = null;

    //Variables specific to distributed DACS and distributed DEE computations
    PrunedRotamers<Boolean> prunedRot = null;
    String rotFile = null;
    boolean useSF = false;
    boolean splitFlags[][][][][][] = null;
    String sfFile = null;
    boolean distrDACS = false;
    boolean distrDEE = false;
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

    String pdbName;

    String eRefMatrix;
    boolean PEMcomp = false; //true if PEM computation is performed; false if mut search is performed
    private boolean templateAlwaysOn;
    private boolean addOrigRots = false;
    private boolean useMaxKSconfs;
    private BigInteger maxKSconfs;
    private int curStrForMatrix;

    boolean useTriples;
    boolean useFlagsAStar;
    boolean magicBulletTriples;
    int magicBulletNumTriples;
    //DEEPer
    boolean doPerturbations;
    boolean minimizePerts;
    String pertFile;
    boolean addWTRot;
    boolean idealizeSC;

    // Generic constructor
    MutationManager(String logName, OneMutation mArray[], boolean PEMcomputation) {
        pdbName = logName;
        if (logName!=null) { //open log file for writing
            try {
                FileOutputStream fileOutputStream = new FileOutputStream(logName);
                BufferedOutputStream bufferedOutputStream = new BufferedOutputStream( fileOutputStream );
                logPS = new PrintStream( bufferedOutputStream );
            } catch (Exception ex) {
                System.out.println("ERROR: An exception occured while opening log file");
            }
        }

        mutArray = mArray;
        cObjArray = new CommucObj[mutArray.length];

        PEMcomp = PEMcomputation;

        distrDACS = false; //if this is a distributed DACS computation, the flag will be set with setDistrDACS()
        distrDEE = false; //if this is a distributed DEE computation, the flag will be set with setDistrDEE()
    }

    // Returns the next mutation packaged in a communication object
    public synchronized CommucObj getNextComObj(int curMutIndex) {

        CommucObj cObj = new CommucObj();
        cObj.numberOfStrands = numberOfStrands;
        cObj.strandMut = strandMut;
        cObj.strandDefault = strandDefault;
        cObj.typeDep = typeDep;
        cObj.addOrigRots = addOrigRots;
        //cObj.ligPresent = ligPresent;
        //cObj.ligType = ligType;
        cObj.distrDEE = distrDEE;
        cObj.distrDACS = distrDACS;
        cObj.arpFilenameMin = arpFilenameMin;
        cObj.arpFilenameMax = arpFilenameMax;
        cObj.params = sParams;
        cObj.stericThresh = stericThresh;
        cObj.softStericThresh = softStericThresh;
        //cObj.numInAS = numInAS;
        cObj.saveTopConfs = saveTopConfs;
        cObj.printTopConfs = printTopConfs;
        cObj.numTopConfs = numTopConfs;
        cObj.rl = rotLib;
        cObj.computeEVEnergy = computeEVEnergy;
        cObj.doMinimization = doMinimization;
        cObj.minimizeBB = minimizeBB;
        cObj.doBackrubs = doBackrubs;
        cObj.backrubFile = backrubFile;
        cObj.calculateVolumes = calculateVolumes;
        cObj.distDepDielect = distDepDielect;
        cObj.dielectConst = dielectConst;
        cObj.doDihedE = doDihedE;
        cObj.doSolvationE = doSolvationE;
        cObj.solvScale = solvScale;
        cObj.vdwMult = vdwMult;
        cObj.PEMcomp = PEMcomp;
        cObj.mutableSpots = mutableSpots;
        cObj.mutRes2Strand = mutRes2Strand;
        cObj.mutRes2StrandMutIndex = mutRes2StrandMutIndex;
        cObj.strandPresent = strandPresent;
        cObj.curMut = curMutIndex;
        cObj.strandLimits = strandLimits;
        cObj.strandsPresent = strandsPresent;
        cObj.templateAlwaysOn = templateAlwaysOn;

        cObj.useFlagsAStar = useFlagsAStar;
        cObj.useTriples = useTriples;
        cObj.doPerturbations = doPerturbations;
        cObj.magicBulletTriples = magicBulletTriples;
        cObj.magicBulletNumTriples = magicBulletNumTriples;
        if(doPerturbations) {
            cObj.minimizePerts = minimizePerts;
            cObj.pertFile = pertFile;
            cObj.addWTRot = addWTRot;
            cObj.idealizeSC = idealizeSC;
        }

        if (PEMcomp) {//PEM computation

            if (!entropyComp) {
                //cObj.numLigRotamers = numLigRotamers;
                cObj.curStrForMatrix = curStrForMatrix;
                cObj.flagMutType = mutArray[curMutIndex].flagMutType;
                cObj.curMut = mutArray[curMutIndex].mutNum;
                cObj.resMut = new int[mutableSpots];
                cObj.currentMutation = new String[mutableSpots];
                for(int i=0; i<mutableSpots; i++) {
                    cObj.resMut[i] = mutArray[curMutIndex].resMut[i];
                    //cObj.currentMutation[i] = resDefault[i];
                }
                int ctr = 0;
                for(int i=0; i<strandDefault.length; i++)
                    for(int j=0; j<strandDefault[i].length; j++) {
                        cObj.currentMutation[ctr] = strandDefault[i][j];
                        ctr++;
                    }

            } else { //entropy energy computation run
                cObj.entropyComp = entropyComp;
                if (compASdist) { //AS-AS distance computation
                    cObj.compASdist = compASdist;
                    cObj.asDist = new boolean[mutableSpots];
                    cObj.dist = dist;
                } else { //AS-AS or SHL-AS energy matrix computation
                    //TODO: fix this code so it makes sense for multiple strands
                    assert false==true;
                    cObj.flagMutType = mutArray[curMutIndex].flagMutType;
                    if (cObj.flagMutType.equalsIgnoreCase("AS-AS")) { //AS-AS run
                        cObj.mutableSpots = 2;
                        cObj.strandMut = new int[1][2];
                        cObj.strandMut[0][0] = mutArray[curMutIndex].resMut[0];
                        cObj.strandMut[0][1] = mutArray[curMutIndex].resMut[1];
                    } else if (cObj.flagMutType.equalsIgnoreCase("INTRA")) { //INTRA run
                        cObj.mutableSpots = 1;
                        cObj.strandMut = new int[1][1];
                        cObj.strandMut[0][0] = curMutIndex;
                    } else {
                        System.out.println("ERROR: only AS-AS and INTRA runs allowed for the pairwise entropy matrix precomputation.");
                        System.exit(1);
                    }
                }
            }
        } else { //mutation search

            if ((!distrDEE)&&(!distrDACS)) { //mutation search run, not (distributed DACS or distributed DEE)

                //cObj.q_L = q_L;
                cObj.pdbName = pdbName;
                cObj.numMutations = numMutations;
                cObj.repeatSearch = repeatSearch;
                cObj.algOption = algOption;
                cObj.numSplits = numSplits;
                cObj.initEw = initEw;
                cObj.scaleInt = scaleInt;
                cObj.maxIntScale = maxIntScale;
                cObj.pruningE = pruningE;
                cObj.stericE = stericE;
                cObj.gamma = gamma;
                cObj.epsilon = epsilon;
                cObj.useMaxKSconfs = useMaxKSconfs;
                cObj.maxKSconfs = maxKSconfs;

                cObj.currentMutation = new String[mutableSpots];
                for(int i=0; i<mutableSpots; i++) {
                    cObj.currentMutation[i] = mutArray[curMutIndex].resTypes[i];
                }
                cObj.bestScore = bestScore;
            } else { //distributed DACS or distributed DEE

                cObj.initEw = initEw;
                cObj.scaleInt = scaleInt;
                cObj.maxIntScale = maxIntScale;
                cObj.prunedRot = prunedRot;
                cObj.useSF = useSF;
                cObj.sfFileIn = sfFile;
                //cObj.numLigRotamers = numLigRotamers;
                cObj.useEref = useEref;
                cObj.eRef = eRef;
                cObj.curStrForMatrix = curStrForMatrix;

                if (distrDACS) { //distributed DACS
                    cObj.numMutations = numMutations;
                    cObj.rotFileIn = rotFile;
                    cObj.pruningE = pruningE;
                    cObj.approxMinGMEC = approxMinGMEC;
                    cObj.lambda = lambda;
                    cObj.algOption = algOption;
                    cObj.initDepth = initDepth;
                    cObj.subDepth = subDepth;
                    cObj.diffFact = diffFact;
                    cObj.minRatioDiff = minRatioDiff;
                    cObj.msp = msp;
                    cObj.numInitUnprunedConf = numInitUnprunedConf;
                    cObj.currentMutation = new String[mutableSpots];
                    cObj.outputPruneInfo = outputPruneInfo;
                    cObj.outputConfInfo = outputConfInfo;
                    cObj.partIndex = new Index3[initDepth];
                    for (int i=0; i<initDepth; i++)
                        cObj.partIndex[i] = mutArray[curMutIndex].index[i];
                    int ctr = 0;
                    for(int i=0; i<strandDefault.length; i++)
                        for(int j=0; j<strandDefault[i].length; j++) {
                            cObj.currentMutation[ctr] = strandDefault[i][j];
                            ctr++;
                        }
                    cObj.bestScore = bestScore;
                } else { //distributed DEE
                    cObj.resMut = new int[mutArray[curMutIndex].resMut.length];
                    for (int i=0; i<cObj.resMut.length; i++) {
                        cObj.resMut[i] = mutArray[curMutIndex].resMut[i];
                    }
                    cObj.AAallowed = new String[numberOfStrands][];
                    for(int str=0; str<numberOfStrands; str++) {
                        cObj.AAallowed[str] = new String[strandMut[str].length];
                    }
                    cObj.currentMutation = new String[mutableSpots];
                    for (int str=0; str<numberOfStrands; str++) {
                        for(int i=0; i<strandMut[str].length; i++) {
                            cObj.AAallowed[str][i] = AAallowed[str][i];
                        }
                    }

                    int ctr = 0;
                    for(int i=0; i<strandDefault.length; i++)
                        for(int j=0; j<strandDefault[i].length; j++) {
                            cObj.currentMutation[ctr] = strandDefault[i][j];
                            ctr++;
                        }
                    cObj.numSpPos = numSpPos;
                    cObj.typeDEE = typeDEE;
                }
            }
        }

        cObj.mutationNumber = curMutIndex;
        curMutIndex++;

        return(cObj);
    }


    // Output a finished mutation to the results file
    public synchronized void processFinishedMutation(CommucObj cObj) {

        if (PEMcomp) { //energy matrix computation
            int countNewEntries = cObj.compEE.length;
            //Update the E matrices computed so far with the new computations supplied
            //	by the current cObj
            if (!entropyComp) { //PEM computation
                if(cObj.flagMutType.equals("INTRA")) {
                    //KER: initialize eref to big E
                    for(int i=0; i<eRef.length; i++)
                        for(int j=0; j<eRef[i].length; j++)
                            eRef[i][j]=Float.MAX_VALUE;
                    //KER: I now output a separate Eref matrix
                    for (int i=0; i<countNewEntries; i++) {
                        eRef[cObj.compEE[i].i1][cObj.compEE[i].i2] = Math.min(eRef[cObj.compEE[i].i1][cObj.compEE[i].i2],cObj.compEE[i].minE);
                    }
                    //KER: Remove all of the bigE
                    for(int i=0; i<eRef.length; i++)
                        for(int j=0; j<eRef[i].length; j++)
                            if(eRef[i][j]==Float.MAX_VALUE)
                                eRef[i][j] = 0.0f;
                    outputObject(eRef,eRefMatrix+".dat");
                } else {
                    for (int i=0; i<countNewEntries; i++) {
                        //We do not use the PairwiseEnergyMatrix.setPairwiseE method because we are setting all kinds of energies (pairwise, template, etc.)
                        //by index
                        pairEMatrixMin.eMatrix[cObj.compEE[i].i1][cObj.compEE[i].i2][cObj.compEE[i].i3][cObj.compEE[i].i4][cObj.compEE[i].i5][cObj.compEE[i].i6] = cObj.compEE[i].minE;
                        if (doMinimization)
                            pairEMatrixMax.eMatrix[cObj.compEE[i].i1][cObj.compEE[i].i2][cObj.compEE[i].i3][cObj.compEE[i].i4][cObj.compEE[i].i5][cObj.compEE[i].i6] = cObj.compEE[i].maxE;
                    }
                }
            } else { //entropy E matrix computation
                if (compASdist) { //AS-AS distance computation
                    asDist[cObj.mutationNumber] = cObj.asDist;
                } else {
                    if (cObj.flagMutType.equalsIgnoreCase("INTRA")) {
                        for (int i=0; i<countNewEntries; i++) {
                            int index1 = 1 + rotLib.getRotamerIndexOffset()[cObj.compEE[i].i2] + cObj.compEE[i].i3;
                            pairEMatrixMinEntropy[cObj.mutationNumber*rotLib.getTotalNumRotamers()+index1][0] = cObj.compEE[i].minE;
                        }
                    } else { //AS-AS run
                        for (int i=0; i<countNewEntries; i++) {
                            int ind1 = -1;
                            int ind2 = -1;
                            int index1 = cObj.compEE[i].i1*rotLib.getTotalNumRotamers() + rotLib.getRotamerIndexOffset()[cObj.compEE[i].i2] + cObj.compEE[i].i3;
                            int index2 = cObj.compEE[i].i4*rotLib.getTotalNumRotamers() + rotLib.getRotamerIndexOffset()[cObj.compEE[i].i5] + cObj.compEE[i].i6;
                            if (index1<rotLib.getTotalNumRotamers()) {
                                ind1 = index1;
                                ind2 = index2-rotLib.getTotalNumRotamers();
                            } else {
                                ind1 = index2;
                                ind2 = index1-rotLib.getTotalNumRotamers();
                            }
                            asasE[cObj.strandMut[0][0]][cObj.strandMut[0][1]][ind1][ind2] = cObj.compEE[i].minE;
                        }
                    }
                }

            }
            //Output mutation information to results file (for resume)
            /*System.out.println("MutNUM: "+cObj.curMut+" produced "+countNewEntries+" new entries.");
            logPS.print("Completed mutation "+cObj.curMut);
            	logPS.print(" SlaveNum "+cObj.slaveNum);
            	logPS.print(" Time "+(cObj.elapsedTime/60.0));
            	if (!entropyComp){
            		for(int i=0;i<cObj.mutableSpots;i++)
            			logPS.print(" "+cObj.resMut[i]);
            		logPS.print(" "+cObj.flagMutType);
            	}
            	logPS.println();
            	logPS.flush();*/
        } else { //mutation search
            if ((!distrDEE)&&(!distrDACS)) { //Hybrid MinDEE-K*, not (distributed DACS or distributed DEE)
                System.out.println("MutNUM: "+cObj.mutationNumber);
                logPS.print("Completed mutation "+cObj.mutationNumber);
                BigDecimal score = new BigDecimal("0.0");
                BigDecimal denom = new BigDecimal("1.0");

                ExpFunction e = new ExpFunction();

                for(int i=0; i<cObj.numComplexes-1; i++) {
                    denom = denom.multiply(cObj.q[i]);
                }
                if (denom.compareTo(new BigDecimal("0.0")) != 0)
                    score = cObj.q[cObj.numComplexes-1].divide(denom,4);
                logPS.print(" Score "+score);

                logPS.print(" Volume "+mutArray[cObj.mutationNumber].vol);
                logPS.print(" SlaveNum "+cObj.slaveNum);
                logPS.print(" Time ");
                for(int i=0; i<cObj.numComplexes; i++)
                    logPS.print((cObj.q_Time[i]/60.0)+" ");
                logPS.print(" InitBest "+cObj.bestScore);
                BigDecimal bs = cObj.bestScore;
                if (score.compareTo(cObj.bestScore) >0)
                    bs = score;
                logPS.print(" FinalBest "+bs);
                for(int i=0; i<cObj.mutableSpots; i++)
                    logPS.print(" "+cObj.currentMutation[i]);

                for(int i=0; i<cObj.numComplexes; i++) {
                    logPS.print(" "+i+"ConfInfo "+cObj.searchNumConfsEvaluated[i]+" "+cObj.searchNumPrunedMinDEE[i]+" "
                                +cObj.searchNumConfsPrunedByS[i]+" "+cObj.searchNumConfsLeft[i]);
                }
                logPS.print(" MinEMinimized ");
                for(int i=0; i<cObj.numComplexes; i++)
                    logPS.print(cObj.bestEMin[i]+" ");
                logPS.print(" MinEUnMinimized ");
                for(int i=0; i<cObj.numComplexes; i++)
                    logPS.print(cObj.bestE[i]+" ");
                logPS.print(" EffectiveEpsilon: ");
                /*for(int i=0;i<cObj.numComplexes;i++)
                	logPS.print(cObj.effEpsilon[i]+" ");*/
                logPS.print(" Partial_q_E ");
                for(int i=0; i<cObj.numComplexes; i++)
                    logPS.print(cObj.q[i]+" ");
                logPS.print(" E_total ");
                for(int i=0; i<cObj.numComplexes; i++)
                    logPS.print(cObj.searchNumConfsTotal[i]+" ");
                logPS.print(" SecondEw ");
                for(int i=0; i<cObj.numComplexes; i++)
                    logPS.print(cObj.repeatEW[i]+" ");
                logPS.print(" E_allPruned ");
                for(int i=0; i<cObj.numComplexes; i++)
                    logPS.print(cObj.allPruned[i]+" ");
                logPS.println();
                if (score.compareTo(bestScore) >0) {
                    logPS.println("BestScoreChange "+bestScore+" to "+score);
                    bestScore = score;
                }
                logPS.flush();
            } else if (distrDACS) { //distributed DACS
                bestScore = bestScore.min(cObj.bestScore);
                pruningE = bestScore.floatValue();
                logPS.print("Completed mutation "+cObj.mutationNumber);
                logPS.print(" Score "+cObj.bestScore);
                logPS.print(" BestScore "+bestScore);
                logPS.print(" PartitionIndices");
                for (int i=0; i<initDepth; i++)
                    logPS.print(" "+cObj.partIndex[i]);
                logPS.print(" Time "+(cObj.elapsedTime/60.0));
                logPS.println();
                logPS.flush();
                System.out.println("Partition "+cObj.mutationNumber+" done; best energy: "+cObj.bestScore);
            } else { //distributed DEE

                if (typeDEE==optPairs) { //pairs DEE
                    boolean tmpSF[][][][][][] = null;
                    tmpSF = (boolean [][][][][][])readFromFile(tmpSF,cObj.sfFileOut);

                    for (int p1=0; p1<tmpSF.length; p1++) {
                        if (tmpSF[p1]!=null) {
                            for (int a1=0; a1<tmpSF[p1].length; a1++) {
                                if (tmpSF[p1][a1]!=null) {
                                    for (int r1=0; r1<tmpSF[p1][a1].length; r1++) {
                                        if (tmpSF[p1][a1][r1]!=null) {
                                            for (int p2=0; p2<tmpSF[p1][a1][r1].length; p2++) {
                                                if (tmpSF[p1][a1][r1][p2]!=null) {
                                                    for (int a2=0; a2<tmpSF[p1][a1][r1][p2].length; a2++) {
                                                        if (tmpSF[p1][a1][r1][p2][a2]!=null) {
                                                            for (int r2=0; r2<tmpSF[p1][a1][r1][p2][a2].length; r2++) {
                                                                if(tmpSF[p1][a1][r1][p2][a2][r2])
                                                                    splitFlags[p1][a1][r1][p2][a2][r2] = true;
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

                    /*for (int i1=0; i1<tmpSF.length; i1++){
                    	for (int i2=0; i2<tmpSF[0].length; i2++){
                    		if (tmpSF[i1][i2]) //get the DE pairs from this computation
                    			splitFlags[i1][i2] = true;
                    	}
                    }*/
                    deleteFile(cObj.sfFileOut);
                    outputObject(splitFlags,sfFile); //must be output after each result read
                } else { //singles DEE
                    Iterator<RotInfo<Boolean>> iter = prunedRot.iterator();
                    while(iter.hasNext()) {
                        //for (int i=0; i<prunedRot.length; i++){
                        RotInfo<Boolean> ri = iter.next();
                        if (cObj.prunedRot.get(ri)) //get the pruned rotamers from this computation
                            prunedRot.set(ri, true);
                    }
                    outputObject(prunedRot,rotFile); //must be output after each result read
                }
            }
        }
    }

    private synchronized Object readFromFile(Object inObj, String inFile) {
        try {
            ObjectInputStream in = new ObjectInputStream(new FileInputStream(inFile));
            inObj = in.readObject();
            in.close();
        } catch (Exception e) {
            System.out.println(e.toString());
            System.out.println("ERROR: An exception occurred while reading from object file");
            System.exit(0);
        }
        return inObj;
    }

    private void outputObject(Object outObj, String outFile) {
        FileOutputStream fout = null;
        FileChannel ch = null;
        FileLock lock = null;
        try {
            fout = new FileOutputStream(outFile);
            ch = fout.getChannel();
            lock = ch.lock();
            ObjectOutputStream out = new ObjectOutputStream(fout);
            out.writeObject(outObj);
            //out.close();
        } catch (Exception e) {
            System.out.println(e.toString());
            System.out.println("ERROR: An exception occurred while writing object file");
            System.exit(0);
        } finally {
            try {
                lock.release();
            } catch (Exception e) {
                System.out.println(e.toString());
                System.out.println("ERROR: unable to release lock on file "+outFile);
                System.exit(1);
            }
        }
    }

    private void deleteFile(String file) {

        Runtime rt = Runtime.getRuntime();
        try {
            File tmpSFMat = new File(file);
            if(tmpSFMat.exists())
                tmpSFMat.delete();
        } catch(Exception ex) {
            System.out.println("Exception: runtime");
            System.out.println(ex.getMessage());
        }
    }

    public void closeLog() {
        if (logPS!=null) {
            logPS.flush();
            logPS.close();
        }
    }

    public void setNumOfStrands(int s) {
        numberOfStrands = s;
    }
    public void setStrandMut(int sm[][]) {
        strandMut = sm;
    }
    public void setStrandDefault(String sd[][]) {
        strandDefault = sd;
    }
    public void setStrandPresent(boolean sp[]) {
        strandPresent = sp;
    }
    public void setStrandLimits(String sl[][]) {
        strandLimits = sl;
    }
    public void setStrandsPresent(int sp) {
        strandsPresent = sp;
    }

    /*public void setLigType(String lt) {
    	ligType = lt;
    }*/
    public void setLigPresent(boolean lp) {
        ligPresent = lp;
    }
    public void setNumMutations(int nm) {
        numMutations = nm;
    }
    public void setResMutatable (int resMut[][]) {
        resMutatable = new boolean[resMut.length][resMut[0].length];
        for (int i=0; i<resMutatable.length; i++) {
            for (int j=0; j<resMutatable[0].length; j++) {
                if (resMut[i][j]==1)
                    resMutatable[i][j] = true;
                else
                    resMutatable[i][j] = false;
            }
        }
    }
    public void setAAallowed(String aal[][]) {
        AAallowed = aal;
    }
    public void setarpFilenameMin(String afnm) {
        arpFilenameMin = afnm;
    }
    public void setarpFilenameMax(String afnm) {
        arpFilenameMax = afnm;
    }
    public void setAlgOption(int ao) {
        algOption = ao;
    }
    public void setNumSplits(int ns) {
        numSplits = ns;
    }
    public void setMinDEEFileName(String mdf) {
        minDEEfile = mdf;
    }
    public void setInitEw(float iew) {
        initEw = iew;
    }
    public void setPruningE(float pe) {
        pruningE = pe;
    }
    public float getPruningE() {
        return pruningE;
    }
    public void setGamma(double g) {
        gamma = g;
    }
    public void setEpsilon(float g) {
        epsilon = g;
    }
    public void setEpsilon(double g) {
        epsilon = (new Double(g)).floatValue();
    }
    public void setParams(ParamSet theParams) {
        sParams = theParams;
    }
    public void setStericThresh(float st) {
        stericThresh = st;
    }
    public void setSoftStericThresh(float st) {
        softStericThresh = st;
    }
    /*public void setNumInAS(int nas) {
    	numInAS = nas;
    }*/

    public void setComputeEVEnergy(boolean ceve) {
        computeEVEnergy = ceve;
    }
    public void setDoMinimization(boolean dm) {
        doMinimization = dm;
    }
    public void setMinimizeBB(boolean mbb) {
        minimizeBB = mbb;
    }
    public void setDoBackrubs(boolean br) {
        doBackrubs = br;
    }
    public void setBackrubFile(String brf) {
        backrubFile = brf;
    }
    public void setRepeatSearch(boolean rs) {
        repeatSearch = rs;
    }
    public void setCalculateVolumes(boolean cv) {
        calculateVolumes = cv;
    }
    /*public void setnumLigRotamers(int nlr) {
    	numLigRotamers = nlr;
    }*/
    public float[][] getErefMatrix() {
        return eRef;
    }

    public void setErefMatrix(float eRef[][]) {
        this.eRef = eRef;
    }
    public void setPairEMatrixMin(PairwiseEnergyMatrix pemMin) {
        pairEMatrixMin = pemMin;
    }
    public void setPairEMatrixMax(PairwiseEnergyMatrix pemMax) {
        pairEMatrixMax = pemMax;
    }
    public void setPrunedRot(PrunedRotamers<Boolean> pr) {
        prunedRot = pr;
    }
    public void setRotFile(String rf) {
        rotFile = rf;
    }
    public void setUseSF(boolean usf) {
        useSF = usf;
    }
    public void setSpFlags(boolean spFlags[][][][][][]) {
        splitFlags = spFlags;
    }
    public void setSfFile(String sff) {
        sfFile = sff;
    }
    public void setDistrDACS(boolean dDACS) {
        distrDACS = dDACS;
    }
    public boolean getDistrDACS() {
        return distrDACS;
    }
    public void setDistrDEE(boolean dDEE) {
        distrDEE = dDEE;
    }
    public void setBestScore(BigDecimal bs) {
        bestScore = bs;
    }
    public void setNumSpPos(int spp) {
        numSpPos = spp;
    }
    public void setMSP(int m[]) {
        msp = m;
    }
    public void setTypeDEE(int t) {
        typeDEE = t;
    }
    public void setInitDepth(int id) {
        initDepth = id;
    }
    public void setSubDepth(int sd) {
        subDepth = sd;
    }
    public void setDiffFact(int df) {
        diffFact = df;
    }
    public void setMinRatioDiff(double mrd) {
        minRatioDiff = mrd;
    }
    public void setNumInitUnprunedConf(BigInteger niuc) {
        numInitUnprunedConf = niuc;
    }
    public void setOutputPruneInfo(String opi) {
        outputPruneInfo = opi;
    }
    public void setOutputConfInfo(String oci) {
        outputConfInfo = oci;
    }
    public void setApproxMinGMEC(boolean amg) {
        approxMinGMEC = amg;
    }
    public void setLambda(float l) {
        lambda = l;
    }
    public void setStericE(double se) {
        stericE = se;
    }
    public void setDistDepDielect(boolean ddd) {
        distDepDielect = ddd;
    }
    public void setDielectConst(double dc) {
        dielectConst = dc;
    }
    public void setDoDihedE(boolean dde) {
        doDihedE = dde;
    }
    public void setDoSolvationE(boolean dse) {
        doSolvationE = dse;
    }
    public void setSolvScale(double ss) {
        solvScale = ss;
    }
    public void setVdwMult(double vm) {
        vdwMult = vm;
    }
    public void setScaleInt(boolean si) {
        scaleInt = si;
    }
    public void setMaxIntScale(float is) {
        maxIntScale = is;
    }
    public void setUseEref(boolean uer) {
        useEref = uer;
    }
    public void setEref(float er[][]) {
        eRef = er;
    }
    public void setLigPartFn(BigDecimal ql) {
        q_L = ql;
    }
    public void setEntropyComp(boolean ec) {
        entropyComp = ec;
    }
    public PairwiseEnergyMatrix getMinEmatrix() {
        return pairEMatrixMin;
    }
    public PairwiseEnergyMatrix getMaxEmatrix() {
        return pairEMatrixMax;
    }
    public void setPairEntropyMatrix(float aae[][][][]) {
        asasE = aae;
    }
    public float [][][][] getPairEntropyEmatrix() {
        return asasE;
    }
    public void setIntraEntropyMatrixMin(float pemMin[][]) {
        pairEMatrixMinEntropy = pemMin;
    }
    public void setASdistMatrix(boolean ad[][]) {
        asDist = ad;
    }
    public void setASdist(float d) {
        dist = d;
    }
    public void setCompASdist(boolean ad) {
        compASdist = ad;
    }
    public boolean [][] getASdistMatrix() {
        return asDist;
    }
    public float [][] getMinEmatrixEntropy() {
        return pairEMatrixMinEntropy;
    }
    public void setMutableSpots(int ms) {
        mutableSpots = ms;
    }
    public void setMut2Strand(int m2s[]) {
        mutRes2Strand = m2s;
    }
    public void setMut2StrandMutIndex(int m2sMi[]) {
        mutRes2StrandMutIndex = m2sMi;
    }
    public void setTypeDep(boolean tD) {
        // TODO Auto-generated method stub
        typeDep = tD;
    }
    public void setRotamerLibrary(RotamerLibrary rl) {
        rotLib = rl;
    }
    public void setErefMatrixName(String erm) {
        eRefMatrix  = erm;
    }

    public void setTemplateAlwaysOn(boolean templAlwaysOn) {
        templateAlwaysOn = templAlwaysOn;
    }

    public void setAddOrigRots(boolean aor) {
        addOrigRots   = aor;
    }

    public void setSaveTopConfs(boolean stc) {
        // TODO Auto-generated method stub
        saveTopConfs = stc;
    }

    public void setPrintTopConfs(boolean ptc) {
        // TODO Auto-generated method stub
        printTopConfs = ptc;
    }

    public void setNumTopConfs(int ntc) {
        // TODO Auto-generated method stub
        numTopConfs = ntc;
    }

    public void setUseMaxKSconfs(boolean useKSconfs) {
        useMaxKSconfs = useKSconfs;
    }

    public void setMaxKSconfs(BigInteger mkc) {
        maxKSconfs = mkc;
    }

    public void setcurStrForMatrix(int cSFM) {
        curStrForMatrix = cSFM;
    }

    public void setUseFlagsAStar(boolean useFlagsAStar) {
        this.useFlagsAStar = useFlagsAStar;
    }
    public void setUseTriples(boolean useTriples) {
        this.useTriples = useTriples;
    }
    public void setDoPerturbations(boolean dp) {
        doPerturbations=dp;
    }
    public void setMinimizePerts(boolean mp) {
        minimizePerts=mp;
    }
    public void setPertFile(String pf) {
        pertFile = pf;
    }
    public void setAddWTRot(boolean awr) {
        addWTRot = awr;
    }
    public void setIdealizeSC(boolean idealizeSC) {
        this.idealizeSC = idealizeSC;
    }
    public void setMagicBulletNumTriples(int magicBulletNumTriples) {
        this.magicBulletNumTriples = magicBulletNumTriples;
    }
    public void setMagicBulletTriples(boolean magicBulletTriples) {
        this.magicBulletTriples = magicBulletTriples;
    }
}
