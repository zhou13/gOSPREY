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
// KSParser.java
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
 */

import java.io.*;
import java.nio.channels.*;
import java.lang.Runtime;
import java.util.*;
import java.lang.Integer;
import java.math.*;
// import com.neva.*;   // Not compatible with linux

import mpi.MPI;
import mpi.MPIException;
import mpi.Status;

/**
 *
 * The main class that sets up and handles the basic OSPREY computation and related functions.
 *
 * The OSPREY functions include:
 * 		doDEE - perform DEE/A* redesign (this includes MinDEE, BD, and BRDEE);
 * 		genStructDEE - generate structures for a selected set of the top doDEE conformations;
 * 		precomputeBackrubs - precompute a list of allowed backrubs for each flexible residue position (used by BRDEE);
 * 		KSMaster - perform K* redesign;
 * 		doSinglePartFn - generate (bound or unbound) structures for the K* ensemble of a given protein-ligand complex;
 * 		doResEntropy - use SCMF to compute the residue entropy for each (non-Pro) residue in a protein.
 *
 */
public class KSParser {

    boolean printSegID = false;

    boolean hElect = true; // should hydrogens be used in electrostatic energy calculations
    boolean hVDW = true; // should hydrogens be used in vdw energy calculations
    boolean hSteric = false; // should hydrogens be used in steric checks

    final double constRT = 1.9891/1000.0 * 298.15;   // in kCal / kelvin-mole (the T value here should be consistent with T in EEF1)

    // Config file name
    String cfgName = "KStar.cfg";

    ParamSet rParams = null; //the main KStar parameters

    // For soft vdw potential
    double softvdwMultiplier = 1.0;
    // For electrostatics
    boolean distDepDielect = true;
    double dielectConst = 1.0;

    //Determine if dihedral and solvation energies should be computed
    boolean doDihedE = false;
    boolean doSolvationE = false;
    double solvScale = 1.0;
    float stericThresh = -10000.0f; // allowed overlap between the vdW radii of two atoms, steric clash if larger overlap
    float softStericThresh = -10000.0f; // soft steric overlap threshold

    //the rotamer library files
    //String rl; //KER: Moved rl to a static rotamer library in EnvironmentVars
    RotamerLibrary[] grl;
    //RotamerLibrary rl = null; //the standard rotamer library for the protein (the Aa library)
    //RotamerLibrary[] grl = null; //the rotamer library for the ligand (could be the AA or non-AA library)

    int numThreads = 1;

    //int numAAallowed = -1; //number of allowed AA types
    //String resAllowed[] = null; //the type of allowed AA
    //int rotamerIndexOffset[] = null; //the rotamer index offset for each allowed AA
    //int totalNumRotamers = -1; //the total number of rotamers for the Lovell rotamer library


    final static int regTag = 1; //regular tag for MPI messages
    final static int updateTag = 2; //used in DACS for updating the best energy found for the different partitions
    static int numProc = 1; //number of processors for MPI

    final static int COMPLEX = -1;

    static boolean mpiRun = false; //determines if this is an MPI run

    //the algorithm options that define what pruning criteria will be applied
    //	NOTE!!! These must be the same as in MutationManager.java
    final int optSimple = 1;
    final int optBounds = 2;
    final int optSplit = 3;
    final int optPairs = 4;

    //the assigned protein, ligand, and cofactor strand numbers
    //final int sysStrNum = 0; //the protein strand number is always 0
    //int ligStrNum = -1;
    //int cofStrNum = -1;

    Molecule[] mols = new Molecule[3];
    RotamerLibrary[][] rotLibs = new RotamerLibrary[3][];

    /**
     * Checks if this is an MPI run and calls the respective functions
     *
     */
    public void checkMPI(String[] args) {

        if ((args.length>0)&&(args[0].equalsIgnoreCase("mpi"))) { //MPI run

            mpiRun = true;
            MPItoThread.initialize(mpiRun, 0);

            String tmp[] = new String[args.length-1]; //remove the mpi argument
            System.arraycopy(args, 1, tmp, 0, tmp.length);
            args = tmp;

            args = parseArgs(args);



            try {
                handleDoMPI(args);
            } catch (Exception e) {};
        } else {
            mpiRun = false;

            args = parseArgs(args);

            MPItoThread.initialize(mpiRun, numThreads);
            //Store all the threads that are available for Kyle's "thread mpi"
            MPItoThread.threadEle.put(Thread.currentThread(), new ThreadElement(0));

            outputProgInfo(); //output program information
            setConfigPars(); //set the parameters from the configuration file

            parse(args); //parse the arguments
        }
    }

    private String[] parseArgs(String[] args) {
        while(args.length>0 && args[0].startsWith("-")) {

            if (args[0].equalsIgnoreCase("-c")) {
                cfgName = args[1];
                String temp []= new String[args.length-2];
                System.arraycopy(args,2,temp,0,args.length-2);
                args = temp;
            } else if(args[0].equalsIgnoreCase("-t")) {
                numThreads = new Integer(args[1]);
                String temp []= new String[args.length-2];
                System.arraycopy(args,2,temp,0,args.length-2);
                args = temp;
            }
        }

        return args;

    }

    /**
    * The main function which handles the OSPREY commands
    */
    public void parse(String[] args) {

        boolean commandLineScript = false;
        boolean firstCommandLine = false;
        byte bytebuff[];
        String s = new String("");  // line being parsed

        if (args.length > 0) {
            commandLineScript = true;
            firstCommandLine = true;
        }

        bytebuff = new byte[150];
        if (!commandLineScript) {
            System.out.print("> ");
            try {
                System.in.read(bytebuff);
            } catch ( Exception e ) {
                System.out.println("ERROR: An error occurred while reading input");
                System.exit(0);
            }
            s = new String(bytebuff);  // create a string from bytebuff
        } else if (commandLineScript && !firstCommandLine) {
            // If you were running a command line script and the file is over then quit
            s = new String("quit");
        } else if (firstCommandLine) {
            s = new String("");
            for(int i=0; i<args.length; i++)
                s = s.concat(args[i] + " ");
            firstCommandLine = false;
        }

        s = s.trim();  // remove whitespace from beginning and end of line

        StringTokenizer st = new StringTokenizer(s," ;\t\n\r\f");
        String firstToken = new String("");
        if (st.hasMoreTokens())
            firstToken = st.nextToken();  // snag a copy of the first token

        if (firstToken.equalsIgnoreCase("doSinglePairE"))
            doSinglePairE(s,null);
        else if (firstToken.equalsIgnoreCase("doResEntropy"))
            handleDoResEntropy(s,null);
        else if (firstToken.equalsIgnoreCase("selectResidues"))
            selectResidues(s);
        else if (firstToken.equalsIgnoreCase("compStericOverlap"))
            handleCompStericOverlap(s);
        else if (firstToken.equalsIgnoreCase("precomputeBackrubs"))
            handlePrecomputeBackrubs(s);

        else if (firstToken.equalsIgnoreCase("doDEE"))
            handleDoDEE(s);
        else if (firstToken.equalsIgnoreCase("genStructDEE"))
            handleMinDEEApplyRot(s);
        else if (firstToken.equalsIgnoreCase("generateRandConfs"))
            generateRandConfs(s);
        else if (firstToken.equalsIgnoreCase("fitEparams"))
            fitEparams(s);

        //KER: It is better to do a threaded debug than continually
        //update this function IMO
        /*else if (firstToken.equalsIgnoreCase("doSinglePartFn"))
        	handleKSTest(s);*/
        else if (firstToken.equalsIgnoreCase("computeEnergyMol"))
            handleComputeEnergyMol(s);
        else if (firstToken.equalsIgnoreCase("KSMaster"))
            handleKSMaster(s);
        else if (firstToken.equalsIgnoreCase("computeEmats"))
            handleComputeAllPairwiseRotamerEnergies(s);

        else if (firstToken.equalsIgnoreCase("genBackbones"))
            generateBackbones(s);
        else if (firstToken.equalsIgnoreCase("identifyRots"))
            identifyRotamers(s);
        else if (firstToken.equalsIgnoreCase("fixStruct"))
            fixStruct(s);

        //exit from all slave nodes
        cleanUpNodes();


    } // End parse function

    public static void cleanUpNodes() {
        if (mpiRun || MPItoThread.exe != null) { //exit from all slave nodes
            CommucObj cObj[] = new CommucObj[1];
            cObj[0] = null;
            for (int curProc=1; curProc<numProc; curProc++) {
                try {
                    MPItoThread.Send(cObj, 0, 1, ThreadMessage.OBJECT, curProc, regTag);
                } catch (Exception e) {}
            }
            if(MPItoThread.exe!=null)
                MPItoThread.exe.shutdown();
        }
    }


    /**
    * Displays the program version and citations
    */
    public void outputProgInfo() {

        System.out.println();
        System.out.println("OSPREY Protein Redesign Software Version 2.0");
        System.out.println("Copyright (C) 2001-2012 Bruce Donald Lab, Duke University");
        System.out.println("");
        System.out.println("This program is free software: you can redistribute it and/or modify");
        System.out.println("it under the terms of the GNU Lesser General Public License as");
        System.out.println("published by the Free Software Foundation, either version 3 of the");
        System.out.println("License, or (at your option) any later version.");
        System.out.println("");
        System.out.println("This program is distributed in the hope that it will be useful,");
        System.out.println("but WITHOUT ANY WARRANTY; without even the implied warranty of");
        System.out.println("MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the");
        System.out.println("GNU Lesser General Public License for more details.");
        System.out.println("");
        System.out.println("There are additional restrictions imposed on the use and distribution");
        System.out.println("of this open-source code, including: (A) this header must be included");
        System.out.println("in any modification or extension of the code; (B) you are required to");
        System.out.println("cite our papers in any publications that use this code.  The citation");
        System.out.println("for the various different modules of our software, together with a");
        System.out.println("complete list of requirements and restrictions are found in the");
        System.out.println("document license.pdf enclosed with this distribution.");
        System.out.println("");

        if(mpiRun)
            System.out.println("OSPREY running on "+numProc+" processor(s)");
        else
            System.out.println("OSPREY running on "+numThreads+" thread(s)");
        System.out.println();

    }

    //Sets the parameters from the configuration file
    public void setConfigPars() {

        rParams = new ParamSet();
        rParams.addParamsFromFile(cfgName);

        // PGC Added default values for all variables.
        hElect = (new Boolean((String)rParams.getValue("HELECT", "true"))).booleanValue();
        hVDW = (new Boolean((String)rParams.getValue("HVDW", "true"))).booleanValue();
        hSteric = (new Boolean((String)rParams.getValue("HSTERIC","false"))).booleanValue();
        distDepDielect = (new Boolean((String)rParams.getValue("DISTDEPDIELECT","true"))).booleanValue();
        dielectConst = (new Double((String)rParams.getValue("DIELECTCONST","6.0"))).doubleValue();
        doDihedE = (new Boolean((String)rParams.getValue("DODIHEDE","false"))).booleanValue();
        doSolvationE = (new Boolean((String)rParams.getValue("DOSOLVATIONE","true"))).booleanValue();
        solvScale = (new Double((String)rParams.getValue("SOLVSCALE","0.5"))).doubleValue();
        softvdwMultiplier = (new Double((String)rParams.getValue("VDWMULT","0.95"))).doubleValue();
        stericThresh = (new Float((String)rParams.getValue("STERICTHRESH","0.4"))).floatValue();
        softStericThresh = (new Float((String)rParams.getValue("SOFTSTERICTHRESH","1.5"))).floatValue();
        EnvironmentVars.setDataDir(rParams.getValue("DATADIR","./"));
        EnvironmentVars.setForcefld(rParams.getValue("FORCEFIELD","AMBER"));
        double entropyScale = (new Double((String)rParams.getValue("ENTROPYSCALE","0.0"))).doubleValue();
        EnvironmentVars.setEntropyScale(entropyScale);

        EnvironmentVars.setAArotLib(EnvironmentVars.getDataDir().concat(rParams.getValue("ROTFILE","LovellRotamer.dat")));
        /*numAAallowed = rl.getNumAAallowed();
        resAllowed = rl.getAAtypesAllowed();
        rotamerIndexOffset = rl.getRotamerIndexOffset();
        totalNumRotamers = rl.getTotalNumRotamers();*/

        EnvironmentVars.autoFix = new Boolean((String)rParams.getValue("AUTOFIX","true")).booleanValue();

        String ramaGlyFile = (String)rParams.getValue("RAMAGLYFILE","rama500-gly-sym.data");

        if( ! ramaGlyFile.equalsIgnoreCase("none") ) {
            String ramaFiles[] = {
                EnvironmentVars.dataDir + ramaGlyFile,
                EnvironmentVars.dataDir + (String)rParams.getValue("RAMAPROFILE","rama500-pro.data"),
                EnvironmentVars.dataDir + (String)rParams.getValue("RAMAGENFILE","rama500-general.data"),
                EnvironmentVars.dataDir + (String)rParams.getValue("RAMAPREPROFILE","rama500-prepro.data")
            };
            RamachandranChecker.getInstance().readInputFiles( ramaFiles );
        }

        AStarConfig.enableJava    = Boolean.parseBoolean(rParams.getValue("ENABLEASTARJAVA", "true"));
        AStarConfig.enableNativeC = Boolean.parseBoolean(rParams.getValue("ENABLEAStarNATIVEC", "false"));
        AStarConfig.enableCUDA    = Boolean.parseBoolean(rParams.getValue("ENABLEAStarCUDA", "false"));
        AStarConfig.maxCPUMemory  = Long.parseLong(rParams.getValue("MAXNATIVECPUMEMORY", "1073741824"));
        AStarConfig.maxGPUMemory  = Long.parseLong(rParams.getValue("MAXNATIVEGPUMEMORY", "1073741824"));
        AStarConfig.numWorkGroup  = Integer.parseInt(rParams.getValue("NUMGPUWORKGROUP", "8"));
        AStarConfig.numWorkItem   = Integer.parseInt(rParams.getValue("NUMGPUWORKITEM", "192"));
        AStarConfig.numWorkItem2  = Integer.parseInt(rParams.getValue("NUMGPUWORKITEM2", "192"));
        AStarConfig.shrinkRatio   = Double.parseDouble(rParams.getValue("SHRINKRATIO", "1"));
    }

    /******************************/
    // This function returns the number of tokens in string s
    private int numTokens(String s) {

        int curNum = 0;
        StringTokenizer st = new StringTokenizer(s," ,;\t\n\r\f");

        while (st.hasMoreTokens()) {
            curNum++;
            st.nextToken();
        }
        return(curNum);
    }


    /******************************/
    // This function returns the xth token in string s
    private String getToken(String s, int x) {

        int curNum = 1;
        StringTokenizer st = new StringTokenizer(s," ,;\t\n\r\f");

        while (curNum < x) {
            curNum++;
            if (st.hasMoreTokens())
                st.nextToken();
            else {
                // System.out.println("ERROR: Unable to access parameter " + x);
                return(new String(""));
            }
        }

        if (st.hasMoreTokens())
            return(st.nextToken());
        return(new String(""));

    } // end getToken

    // Computes the factorial of input n
    public static BigInteger factorial(int n) {

        if (n==0)
            return BigInteger.valueOf(1);
        if (n<0) {
            System.out.println("ERROR: trying to find factorial of: "+n);
            System.exit(0);
        }

        return (factorial(n-1).multiply(BigInteger.valueOf(n)));
    }

    public void saveMolecule(Molecule m, String fname, float energy) {

        try {
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
        } catch (IOException e) {
            System.out.println("ERROR: An io exception occurred while writing file: "+fname);
            System.exit(0);
        } catch ( Exception e ) {
            System.out.println(e.toString());
            System.out.println("ERROR: An exception occurred while writing file: "+fname);
            System.exit(0);
        }
    }

    // This function generates all possible combinations of n choose m
    public void generateCombinations(int residueMutatable[][], int n, int m) {

        int curIndex[] = new int[1];
        int curComb[] = new int[n];
        curIndex[0] = 0;
        generateCombHelper(0,n,curIndex,residueMutatable,curComb,0,m);
    }
    private void generateCombHelper(int depth, int maxDepth, int curIndex[], int
                                    residueMutatable[][], int curComb[], int numUsed, int maxToUse) {

        if (depth >= maxDepth) {
            if (numUsed == maxToUse) {
                for (int i=0; i<maxDepth; i++) {
                    residueMutatable[curIndex[0]][i] = curComb[i];
                }
                curIndex[0]++;
            }
            return;
        }

        curComb[depth] = 0;
        generateCombHelper(depth+1,maxDepth,curIndex,residueMutatable,curComb,numUsed,maxToUse);

        if (numUsed < maxToUse) {
            curComb[depth] = 1;
            generateCombHelper(depth+1,maxDepth,curIndex,residueMutatable,curComb,numUsed+1,maxToUse);
        }
    }
    // end combination code

    //Sets the allowables for AS residue position curPos
    private void setAllowablesHelper(RotamerSearch rs, ParamSet sParams, boolean addWT, int strNum, int molStrand, int curPos, int[][] strandMut, String strandDefault[][]) {
        String tempResAllow = (String)sParams.getValue("RESALLOWED"+strNum+"_"+curPos, "");
        if(numTokens(tempResAllow) <= 0 && !addWT) {
            System.out.println("Error: resAllowed not set for strand");
            System.exit(1);
        }
        for(int q=0; q<numTokens(tempResAllow); q++)
            rs.setAllowable(strandMut[molStrand][curPos],getToken(tempResAllow,q+1),molStrand);
        if (addWT)
            rs.setAllowable(strandMut[molStrand][curPos],strandDefault[molStrand][curPos],molStrand); //the default type is set last
    }

    //Sets up the molecule system and returns the number of ligand rotamers
    private Molecule setupMolSystem(Molecule m, ParamSet sParams, boolean[] strandPresent, String[][] strandLimits) {
        return setupMolSystem(m,sParams,strandPresent,strandLimits,true);
    }

    //Sets up the molecule system and returns the number of ligand rotamers
    private Molecule setupMolSystem(Molecule m, ParamSet sParams, boolean[] strandPresent, String[][] strandLimits, boolean keepCofactor) {

        int ligStrNum = -1;
        int cofStrNum = -1;

        try {
            FileInputStream is = new FileInputStream((String)sParams.getValue("PDBNAME"));
            new PDBChemModel(m, is);
        } catch (Exception e) {
            System.out.println("WARNING: An error occurred while reading file");
            System.out.println(e);
            System.exit(1);
        }

        /*m.strand[sysStrNum].isProtein = strandPresent[sysStrNum];		// main protein
        m.strand[sysStrNum].rotTrans = (new Boolean((String)sParams.getValue("STRANDROTTRANS"+sysStrNum))).booleanValue();
        strandLength[sysStrNum] = m.mapPDBresNumToMolResNum(strandLimits[0][1])-m.mapPDBresNumToMolResNum(strandLimits[0][0])+1;*/
        int strNum = 0; //the current strand number; 0 is reserved for the protein strand, the ligand strand is 1 (if present)

        //get the number of strands that are present
        int numPresent = 0;
        for(int str=0; str<strandPresent.length; str++)
            if(strandPresent[str])
                numPresent++;

        Molecule newMol = new Molecule();

        grl = new RotamerLibrary[numPresent];
        //int curPresStr = 0;
        for(int i=0; i<strandLimits.length; i++) {
            String pdbStart = strandLimits[i][0];
            String pdbEnd   = strandLimits[i][1];
            //if (pdbEnd>=0){ //with ligand in PDB
            int molStartNum = m.mapPDBresNumToMolResNum(pdbStart);
            int molEndNum   = m.mapPDBresNumToMolResNum(pdbEnd);
            if(molStartNum <0 || molEndNum < 0) {
                System.out.println("Please make sure strand "+i+"'s begin and end are set properly.");
                System.exit(0);
            }
            int numInStrand = molEndNum-molStartNum+1;
            //strandLength[i] = numInStrand;
            Residue myStrand[] = new Residue[numInStrand];
            for (int j=(numInStrand-1); j>=0; j--) {
                myStrand[j] = m.residue[molStartNum+j]; // pull out the ligand
                m.deleteResidue(molStartNum+j);
                myStrand[j].renumberResidue();
            }
            if (strandPresent[i]) { //ligand will be used in design

                newMol.addStrand(""+i);
                strNum = newMol.strand.length-1;
                newMol.strand[strNum].rotTrans = (new Boolean((String)sParams.getValue("STRANDROTTRANS"+i))).booleanValue();

                for(int j=0; j<numInStrand; j++)
                    newMol.addResidue(strNum,myStrand[j],false);

                newMol.strand[strNum].isProtein = (new Boolean((String)sParams.getValue("STRANDAA"+i))).booleanValue();
                if (newMol.strand[strNum].isProtein) //use the AA rotamer library for the ligand
                    grl[strNum] = EnvironmentVars.aaRotLib;
                else //use the non-AA rotamer library for the ligand
                    grl[strNum] = new RotamerLibrary(sParams.getValue("GROTFILE"+i,"GenericRotamers.dat"),false);

                //change ligand to the specified residue type
                //if ( m.strand[ligStrNum].isProtein && !m.strand[ligStrNum].residue[0].name.equalsIgnoreCase(ligType) ) //not the same ligand type
                //	(new StrandRotamers(grl,m.strand[ligStrNum])).changeResidueType(m,0,ligType,true);

                strNum++;
                //curPresStr++;

            }
            //}
            /*else if (strandPresent[i]){
            System.out.println("ERROR: Attempting to use a ligand, but ligand not found in system config file");
            System.exit(1);
            }*/
        }


        //Get the cofactors (if present)
        if(keepCofactor) {
            int str=-1;
            for(int i=0; i<strandPresent.length; i++) {
                if(strandPresent[i]) {
                    str++;
                    String cofMapString = sParams.getValue("COFMAP"+i, "-1");
                    //int numCofactorRes = (new Integer((String)sParams.getValue("NUMCOFRES"))).intValue();

                    Residue cof;
                    for(int res=0; res<numTokens(cofMapString); res++) {
                        int cofactorRes = m.mapPDBresNumToMolResNum(getToken(cofMapString, res));
                        if(cofactorRes >= 0) {
                            cof = m.residue[cofactorRes];
                            cof.cofactor = true;
                            m.deleteResidue(cofactorRes);
                            newMol.addResidue(str,cof,false);
                        }
                    }
                }
            }
        }

        newMol.determineBonds();
        newMol.establishConnectivity(false);

        return newMol;

        //Determine the number of rotamers for the ligand (if used)
        /*int numLigRotamers = 0;
        if (useLig) {
        	numLigRotamers = grl.getNumRotamers(m.strand[ligStrNum].residue[0].name);
        	if (numLigRotamers == 0)
        		numLigRotamers = 1;
        }
        return numLigRotamers;*/
    }

    // Computes the bound or unbound partition function and
    //  can compute the energyminimized structure for a specified
    //  rotameric conformation
    /*public void handleKSTest(String s) {

    	// Takes the following parameters
    	// 1: System parameter filename (string)
    	// 2: Mutataion search parameter filename (string)

    	ParamSet sParams = new ParamSet();
    	sParams.addParamsFromFile(getToken(s,2)); //read system parameters
    	sParams.addParamsFromFile(getToken(s,3)); //read mutation search parameters

    	// Pull search parameters
    	//int numInAS = (new Integer((String)sParams.getValue("NUMINAS"))).intValue();
    	String eMatrixNameMin = (String)sParams.getValue("MINENERGYMATRIXNAME");
    	String eMatrixNameMax = (String)sParams.getValue("MAXENERGYMATRIXNAME");
    	boolean doMinimize = (new Boolean((String)sParams.getValue("DOMINIMIZE"))).booleanValue();
    	boolean minimizeBB = (new Boolean((String)sParams.getValue("MINIMIZEBB"))).booleanValue();
    	boolean doBackrubs = (new Boolean((String)sParams.getValue("DOBACKRUBS"))).booleanValue();
    	boolean repeatSearch = (new Boolean((String)sParams.getValue("REPEATSEARCH"))).booleanValue();
    	String backrubFile = (String)sParams.getValue("BACKRUBFILE");
    	boolean scaleInt = (new Boolean((String)sParams.getValue("SCALEINT"))).booleanValue();
    	float maxIntScale = (new Float((String)sParams.getValue("MAXINTSCALE"))).floatValue();
    	float initEw = (new Float((String)sParams.getValue("INITEW"))).floatValue();
    	double pruningE = (new Double((String)sParams.getValue("PRUNINGE"))).doubleValue();
    	double stericE = (new Double((String)sParams.getValue("STERICE"))).doubleValue();
    	//boolean ligPresent = (new Boolean((String)sParams.getValue("LIGPRESENT"))).booleanValue();
    	//String ligType = (String)sParams.getValue("LIGTYPE");
    	boolean saveConfs = (new Boolean((String)sParams.getValue("OUTPUTPDBS"))).booleanValue();
    	boolean KSGMEC = (new Boolean((String)sParams.getValue("KSGMEC"))).booleanValue();
    	boolean KSCONFTHRESH = (new Boolean((String)sParams.getValue("KSCONFTHRESH"))).booleanValue();
    	String tmpNumKSconfs = (String)sParams.getValue("numKSconfs");
    	BigInteger numKSconfs = new BigInteger(tmpNumKSconfs);

    	String fName = (String)sParams.getValue("PDBPREFIX");
    	String resMut = (String)sParams.getValue("RESMUT");
    	//String mutNum = (String)sParams.getValue("mutNum");
    	//float epsilon = (new Float((String)sParams.getValue("EPSILON"))).floatValue();

    	if (!doMinimize)
    		minimizeBB = false;
    	if (!minimizeBB)
    		doBackrubs = false;

    	if ( (!ligPresent) && ((new Boolean((String)sParams.getValue("USEUNBOUNDSTRUCT"))).booleanValue()) ) { //ligPresent, or a different input structure is used for the unbound partition function computation
    		sParams.setValue("PDBNAME",sParams.getValue("UNBOUNDPDBNAME"));
    		sParams.setValue("PDBLIGNUM","-1");
    		eMatrixNameMin = sParams.getValue("MINENERGYMATRIXNAMEUNBOUND");
    		eMatrixNameMax = sParams.getValue("MAXENERGYMATRIXNAMEUNBOUND");
    	}

    	//Setup the molecule system
    	Molecule m = new Molecule();
    	int numLigRotamers = setupMolSystem(m,sParams,ligPresent,ligType);

    	int residueMap[] = new int[numInAS];
    	String resDefault[] = new String[numInAS];
    	String resMapString = (String)sParams.getValue("RESIDUEMAP");
    	System.out.print("ResidueMap:");
    	for(int i=0;i<numInAS;i++){
    		int pdbResNum = (new Integer(getToken(resMapString,i+1))).intValue();
    		residueMap[i] = m.strand[sysStrNum].mapPDBresNumToStrandResNum(pdbResNum);
    		resDefault[i] = m.strand[sysStrNum].residue[residueMap[i]].name;
    		System.out.print(" "+residueMap[i]+"("+m.strand[sysStrNum].residue[residueMap[i]].fullName+")");
    	}
    	System.out.println();

    	RotamerSearch rs = new RotamerSearch(m,sysStrNum,ligStrNum,hElect,hVDW,hSteric,true,true,epsilon,stericThresh,softStericThresh,distDepDielect,dielectConst,doDihedE,doSolvationE,solvScale,softvdwMultiplier,rl,grl);

    	// Define the mutation amino-acid sequence
    	System.out.print("Mutation Sequence:");
    	String curSeq[] = new String[numInAS];
    	for(int i=0;i<numInAS;i++){
    		curSeq[i] = getToken(resMut,i+1);
    		System.out.print(" "+curSeq[i]);
    	}
    	System.out.println();

    	System.out.println("Beginning setAllowables");
    	for(int i=0;i<numInAS;i++){
    		rs.setAllowable(residueMap[i],curSeq[i]);
    	}

    	System.out.print("Loading precomputed min energy matrix...");
    	loadPairwiseEnergyMatrices(sParams,rs,eMatrixNameMin+".dat",doMinimize,eMatrixNameMax+".dat");
    	System.out.println("done");

    	BigDecimal q_L = BigDecimal.ZERO;
    	if (ligPresent)
    		q_L = getLigPartFn(m,numInAS,ligType,eMatrixNameMin+".dat"); //compute the ligand partition function

    	System.out.println("Before start");

    	boolean prunedRotAtRes[] = new boolean[numInAS*totalNumRotamers+numLigRotamers];
    	for (int i=0; i<prunedRotAtRes.length; i++)
    		prunedRotAtRes[i] = false;

    	//Prune all rotamers that are incompatible with the template (intra E + re-to-template E >= stericE)
    	prunedRotAtRes = rs.DoPruneStericTemplate(numInAS, totalNumRotamers, numLigRotamers,
    			residueMap, rotamerIndexOffset, prunedRotAtRes, stericE);

    	if (doMinimize) //precompute the interval terms in the MinDEE criterion
    		rs.doCompMinDEEIntervals(numInAS, totalNumRotamers, numLigRotamers, residueMap,
    				rotamerIndexOffset, prunedRotAtRes, scaleInt, maxIntScale);

    	prunedRotAtRes = rs.DoDEEGoldstein(numInAS, totalNumRotamers, numLigRotamers, residueMap,
    			rotamerIndexOffset, initEw, prunedRotAtRes, doMinimize, false, minimizeBB);

    	//Prune with MinBounds
    	prunedRotAtRes = rs.DoMinBounds(numInAS,totalNumRotamers,numLigRotamers,
    			residueMap,rotamerIndexOffset,pruningE,prunedRotAtRes,initEw, false, false);

    	//Compute the Ec value and prunedIsSteric[]
    	rs.DoMinBounds(numInAS,totalNumRotamers,numLigRotamers,
    			residueMap,rotamerIndexOffset,pruningE,prunedRotAtRes,initEw, false, true);

    	BigDecimal initialBest = BigDecimal.ZERO;
    	if (ligPresent)
    		initialBest =  q_E.multiply(bestScore.multiply(q_L)).multiply(new BigDecimal(gamma * epsilon));

    	//Do the rotamer search
    	rs.slaveDoRotamerSearch(true,doMinimize,numInAS,numAAallowed,totalNumRotamers,rotamerIndexOffset,resAllowed,
    			residueMap,ligPresent,initialBest,null,minimizeBB,saveConfs,fName,doBackrubs,backrubFile);

    	if ((repeatSearch)&&(rs.repeatSearch)){ //the desired accuracy was not achieved, so repeat the search: the setup is already done

    		System.out.println();
    		System.out.println("Repeating search..");
    		rs.repeatSearch = false; //reset the flag
    		rs.slaveDoRotamerSearch(true,doMinimize,numInAS,numAAallowed,totalNumRotamers,rotamerIndexOffset,resAllowed,
    				residueMap,ligPresent,initialBest,null,minimizeBB,saveConfs,fName,doBackrubs,backrubFile);
    	}
    }*/


    // Finds the energy for a given input system (a molecule with specified flexible residues)
    public void handleComputeEnergyMol(String s) {

        // Takes the following parameters
        // 1: System parameter filename (string)
        // 2: Ligand (boolean), is true if present
        // 3: Amino acid type for ligand (if ligand is absent, write none or anything)

        ParamSet sParams = new ParamSet();
        sParams.addParamsFromFile(getToken(s,2)); //read system parameters


        MolParameters mp = new MolParameters();
        loadStrandParams(sParams, mp, COMPLEX);

        //Check to see if the pdb name is a directory
        Object pdbName = (String)sParams.getValue("PDBNAME");
        File f = new File((String) pdbName);
        Object[] pdbFiles = null;
        if(f.isDirectory()) {
            pdbFiles = getPdbFiles(f);
        } else {
            pdbFiles = new String[1];
            pdbFiles[0] = pdbName;
        }

        System.out.println("Starting energy computation");
        for(int q = 0; q<pdbFiles.length ; q++) {
            //Change the pdbfile to be looked at;
            sParams.setValue("PDBNAME", (String) pdbFiles[q]);
            //Setup the molecule system
            Molecule m = new Molecule();
            m = setupMolSystem(m,sParams,mp.strandPresent,mp.strandLimits);

            Amber96ext a96ff = new Amber96ext(m, distDepDielect, dielectConst, doSolvationE, solvScale, softvdwMultiplier);
            a96ff.calculateTypesWithTemplates();
            a96ff.initializeCalculation();
            a96ff.setNBEval(hElect,hVDW);
            //TODO: Fix this so that ligand numbers can get set in the energy function
            /*if (ligPresent)
            	a96ff.setLigandNum((new Integer((String)sParams.getValue("PDBLIGNUM"))).intValue());
            */

            /*boolean specificInt = true;
            if(specificInt){

            	loadMutationParams(sParams, mp);

            	String flag = "SHL-AS";int str1 = 1; int res1 = 0;int str2 = -1; int res2 = -1;
            	RotamerSearch rs = new RotamerSearch(m,-1, mp.strandsPresent, hElect, hVDW, hSteric, true,
            			true, 0.0f, stericThresh, softStericThresh, distDepDielect, dielectConst, doDihedE,
            			doSolvationE, solvScale, softvdwMultiplier, rl, grl, forcefield);
            	//rs.strandMut = strandMut;
            	float energy = rs.getPairE(flag, str1, res1, str2, res2,mp.strandMut);
            	System.out.println(""+pdbFiles[q]+": Energy: "+energy);
            }
            else{*/
            double energy[] = a96ff.calculateTotalEnergy(m.actualCoordinates,-1);
            System.out.println("System energy: " + energy[0]+" (elect: "+energy[1]+" vdW: "+energy[2]+" solvation: "+energy[3]+")");
            //}

        }
    }

    //input should be a directory
    public Object[] getPdbFiles(File f) {
        Vector<String> pdbFiles = new Vector<String>();
        File allMyFolderObjects[]  = f.listFiles();
        for(int i =0; i<allMyFolderObjects.length; i++) {
            String filename = allMyFolderObjects[i].getName();
            String ext = (filename.lastIndexOf(".")==-1)?"":filename.substring(filename.lastIndexOf(".")+1,filename.length());
            if(ext.equals("pdb"))
                pdbFiles.add(""+f.getPath()+"\\"+filename);
        }

        return pdbFiles.toArray();

    }

    /**
     * Performs K* redesign; sets up the K* computation from the input model and configuration files and distributes the
     * candidate mutants for evaluation by the set of available processors.
    */
    public void handleKSMaster(String s) {

        // Takes the following parameters
        // 1: System parameter filename (string)
        // 2: Mutation search parameter filename (string)

        ParamSet sParams = new ParamSet();
        sParams.addParamsFromFile(getToken(s,2)); //read system parameters
        sParams.addParamsFromFile(getToken(s,3)); //read mutation search parameters


        //int numInAS = (new Integer((String)sParams.getValue("NUMINAS"))).intValue();
        int numMutations = (new Integer((String)sParams.getValue("NUMMUTATIONS"))).intValue();
        String runName = (String)sParams.getValue("RUNNAME");
        String mutFileName = (String)sParams.getValue("MUTFILENAME",runName+".mut");
        String eMatrixNameMin = (String)sParams.getValue("MINENERGYMATRIXNAME",runName+"minM");
        String eMatrixNameMax = (String)sParams.getValue("MAXENERGYMATRIXNAME",runName+"maxM");
        boolean doMinimize = (new Boolean((String)sParams.getValue("DOMINIMIZE","false"))).booleanValue();
        boolean minimizeBB = (new Boolean((String)sParams.getValue("MINIMIZEBB","false"))).booleanValue();
        boolean doBackrubs = (new Boolean((String)sParams.getValue("DOBACKRUBS","false"))).booleanValue();
        String backrubFile = "";
        if(doBackrubs) {
            backrubFile = (String)sParams.getValue("BACKRUBFILE");
        }
        boolean repeatSearch = (new Boolean((String)sParams.getValue("REPEATSEARCH","true"))).booleanValue();
        boolean scaleInt = (new Boolean((String)sParams.getValue("SCALEINT","false"))).booleanValue();
        float maxIntScale =0.0f;
        if(scaleInt) {
            maxIntScale = (new Float((String)sParams.getValue("MAXINTSCALE"))).floatValue();
        }
        float initEw = (new Float((String)sParams.getValue("INITEW","6.0"))).floatValue();
        float pruningE = (new Float((String)sParams.getValue("PRUNINGE","100.0"))).floatValue();
        double stericE = (new Double((String)sParams.getValue("STERICE","100.0"))).doubleValue();
        float targetVol = (new Float((String)sParams.getValue("TARGETVOLUME","0.0"))).floatValue();
        float volWindow = (new Float((String)sParams.getValue("VOLUMEWINDOW","50000000"))).floatValue();
        boolean resumeSearch = (new Boolean((String)sParams.getValue("RESUMESEARCH","false"))).booleanValue();
        String resumeFilename = "";
        if(resumeSearch) {
            resumeFilename = (String)sParams.getValue("RESUMEFILENAME");
        }
        double gamma = (new Double((String)sParams.getValue("GAMMA", "0"))).doubleValue();
        float epsilon = (new Float((String)sParams.getValue("EPSILON","0.3"))).floatValue();

        boolean saveTopConfs = (new Boolean((String)sParams.getValue("SAVETOPCONFSASPDB","false"))).booleanValue();
        boolean printTopConfs = (new Boolean((String)sParams.getValue("SAVETOPCONFSROTS","false"))).booleanValue();
        int numTopConfs = (new Integer((String)sParams.getValue("NUMTOPCONFSTOSAVE","0"))).intValue();

        if(printTopConfs || saveTopConfs) {
            //KER: make directory for the confs to be printed to.
            File ksConfDir = new File(EnvironmentVars.ksConfDir);
            if(!ksConfDir.exists())
                ksConfDir.mkdir();
        }

        /*KER: Set environment variables*/
        boolean useMaxKSconfs = new Boolean((String)sParams.getValue("useMaxKSconfs","false")).booleanValue();
        BigInteger maxKSconfs = BigInteger.ZERO;
        if(useMaxKSconfs)
            maxKSconfs = new BigInteger(sParams.getValue("maxKSconfs"));


        // DEEPer parameters
        boolean doPerturbations = (new Boolean((String)sParams.getValue("DOPERTURBATIONS","false"))).booleanValue();//Triggers DEEPer
        String pertFile = (String)sParams.getValue("PERTURBATIONFILE","defaultPerturbationFileName.pert");//Input file giving perturbation information
        boolean minimizePerts = (new Boolean((String)sParams.getValue("MINIMIZEPERTURBATIONS","false"))).booleanValue();//Allow continuous minimization with respect to perturbation parameters
        boolean selectPerturbations = (new Boolean((String)sParams.getValue("SELECTPERTURBATIONS","false"))).booleanValue();//Should perturbations be automatically selected?

        Perturbation.idealizeSC = (new Boolean((String)sParams.getValue("IDEALIZESIDECHAINS","true"))).booleanValue();


        boolean addWTRotsSomehow = (new Boolean((String)sParams.getValue("ADDWTROTS","false"))).booleanValue();
        boolean addOrigRots = false, addWTRot = false;
        if( addWTRotsSomehow ) {
            if(doPerturbations)//DEEPer allows adding WT rotamers only to the positions they came from, so do that if possible
                addWTRot = true;
            else//Otherwise just add them to the rotamer library
                addOrigRots = true;
        }


        /*if (!mpiRun){
        	System.out.println("ERROR: Distributed computation requires MPI");
        	System.exit(1);
        }*/

        System.out.println("Run Name: "+runName);
        System.out.println("Precomputed Min Energy Matrix: "+eMatrixNameMin);
        System.out.println("Precomputed Max Energy Matrix: "+eMatrixNameMax);
        //System.out.println("Ligand Type: "+ligType);
        System.out.println("Volume Center: "+targetVol);
        System.out.println("Volume Window Size: "+volWindow);
        System.out.println("Num Residues Allowed to Mutate: "+numMutations);

        if(resumeSearch) {
            System.out.println("** Resuming Search **");
            System.out.println("     resuming from file: "+resumeFilename);
        }

        MolParameters mp = loadMolecule(sParams, COMPLEX);
        //KER: This is a placeholder so I don't have to change all the variables in the code
        Molecule m = mp.m;
        int numberMutable = mp.numberMutable;
        int strandsPresent = mp.strandsPresent;
        String[][] strandLimits = mp.strandLimits;
        boolean[] strandPresent = mp.strandPresent;
        int[][] strandMut = mp.strandMut;
        String[][] strandDefault = mp.strandDefault;

        if(addOrigRots)
            RotamerLibrary.addOrigRots(strandMut, EnvironmentVars.aaRotLib, m);

        // Create the mutation list with estimated energies
        Set<OneMutation> mutSet = new TreeSet<OneMutation>();


        // Generate all combinations (include (n choose m), (n choose m-1), ... , (n choose 1), and (n choose 0) )
        int numCombAll = 0;
        if(numMutations > mp.numberMutable)
            numMutations = mp.numberMutable;
        for (int i=numMutations; i>=0; i--)
            numCombAll += factorial(mp.numberMutable).divide(factorial(mp.numberMutable-i).multiply(factorial(i))).intValue();
        int residueMutatableAll[][] = new int[numCombAll][mp.numberMutable];
        int curInd = 0;
        for (int i=numMutations; i>=0; i--) {
            int numCombCur = factorial(mp.numberMutable).divide(factorial(mp.numberMutable-i).multiply(factorial(i))).intValue();
            int residueMutatableCur[][] = new int[numCombCur][mp.numberMutable];
            generateCombinations(residueMutatableCur,mp.numberMutable,i);
            for (int j=0; j<numCombCur; j++) {
                residueMutatableAll[curInd] = residueMutatableCur[j];
                curInd++;
            }
        }

        // At this point each row of residueMutatble is a 0/1 array, 1 indicates
        //  that that residues can mutate

        if(selectPerturbations)//Need to run the automatic perturbation selection
            //This only needs to be done once though: after that the perturbations can be read from pertFile
            selectPerturbations(mp, doPerturbations, pertFile, minimizePerts, addWTRot, sParams);


        System.out.print("Checking if precomputed energy matrix is already computed...");
        RotamerSearch rs = new RotamerSearch(m,numberMutable, strandsPresent,hElect,hVDW,hSteric,true,true,
                                             epsilon,stericThresh,softStericThresh,distDepDielect,dielectConst,doDihedE,doSolvationE,solvScale,
                                             softvdwMultiplier,grl, doPerturbations, pertFile, minimizePerts, false, false);


        if(doPerturbations)
            rs.setupRCs(addWTRot);

        for(int i=0; i<m.numberOfStrands; i++) {
            if ((new Boolean((String)sParams.getValue("USEUNBOUNDSTRUCT"+i, "false"))).booleanValue()) { //a different input structure is used for the unbound partition function computation
                ParamSet ubParams = new ParamSet(); //create a new parameter set, just for the unbound-case matrix computation; sParams must not be changed here
                ubParams.setParamsValues(sParams.getParams(), sParams.getValues(), sParams.getCurNum());
                ubParams.setValue("PDBNAME",sParams.getValue("UNBOUNDPDBNAME"+i));
                ubParams.setValue("NUMOFSTRANDS", "1");
                ubParams.setValue("STRAND0", sParams.getValue("STRAND"+i));
                ubParams.setValue("STRANDAA0", sParams.getValue("STRANDAA"+i));
                ubParams.setValue("STRANDROTRANS0", "FALSE");
                ubParams.setValue("STRANDMUTNUMS", getToken(sParams.getValue("STRANDMUTNUMS"),i+1));
                ubParams.setValue("STRANDMUT0", sParams.getValue("STRANDMUT"+i));
                ubParams.setValue("COFMAP0", sParams.getValue("COFMAP"+i,"-1"));

                //KER: Fix the ResAllowed Flags
                for(int q=0; q<strandMut[i].length; q++) {
                    ubParams.setValue("RESALLOWED0_"+q, sParams.getValue("RESALLOWED"+i+"_"+q));
                }

                System.out.print("Checking if precomputed energy matrix (unbound) is already computed...");
                rs = new RotamerSearch(m,numberMutable, strandsPresent,hElect,hVDW,hSteric,true,true,
                                       epsilon,stericThresh,softStericThresh,distDepDielect,dielectConst,doDihedE,doSolvationE,solvScale,
                                       softvdwMultiplier,grl, doPerturbations, pertFile, minimizePerts, false, false);


                loadUnboundPairwiseEnergyMatrices(ubParams,rs,eMatrixNameMin,true,eMatrixNameMax,i);
                rs = null;
                ubParams = null;
                System.out.println("done");
                //BAD CODE //setupMolSystem(m,sParams,ligPresent,ligType); //re-initialize, since some molecule-relative variables have changed (e.g., ligStrNum)
                //VERY BAD CODE
                m = setupMolSystem(m,sParams,strandPresent,strandLimits);
            } else {
                rs.resetMatrices();
                loadPairwiseEnergyMatrices(sParams,rs,eMatrixNameMin,doMinimize,eMatrixNameMax,i);
            }
        }

        rs = null;
        System.out.println("done");

        m = new Molecule();

        m = setupMolSystem(m,sParams,strandPresent,strandLimits);
        rs = new RotamerSearch(m,numberMutable, strandsPresent,hElect,hVDW,hSteric,true,true,
                               epsilon,stericThresh,softStericThresh,distDepDielect,dielectConst,doDihedE,doSolvationE,solvScale,
                               softvdwMultiplier,grl, doPerturbations, pertFile, minimizePerts, false, false);
        //Compute the matrices for all strands (-1 is for the complex)
        //KER: Put this last so the correct Eref matrix is kept
        rs.resetMatrices();
        loadPairwiseEnergyMatrices(sParams,rs,eMatrixNameMin,doMinimize,eMatrixNameMax,COMPLEX);

        //Load mutation list for distribution
        mutSet = handleHybridKSLoadMutList(mutFileName, numberMutable, m, numCombAll, residueMutatableAll,
                                           sParams, strandMut, strandDefault, rs.strandRot, targetVol, volWindow,strandsPresent,strandPresent);

        BigDecimal bestScore = new BigDecimal("0.0"); //for the resume results
        // If doing a resume, read the initial results into a bunch of OneMutations
        OneMutation resumeResults[] = new OneMutation[mutSet.size()];
        if (resumeSearch) {

            bestScore = new BigDecimal("0.0"); //high scores are better

            //OneMutation resumeResults[] = new OneMutation[mutArray.length];
            for(int q=0; q<mutSet.size(); q++)
                resumeResults[q] = new OneMutation();
            resumeResults = readResumeFile(resumeResults,resumeFilename,numberMutable,false,false,-1);
            if(resumeResults != null) {
                System.out.println("Read "+resumeResults.length+" completed mutations");

                // Now filter removed mutations (already computed results
                //  are NOT written to file since you already have them)
                // We do need to maintain the best score
                int newIndex = 0;
                Vector<OneMutation> newVector2 = new Vector<OneMutation>();
                Iterator<OneMutation> i = mutSet.iterator();
                while(i.hasNext()) {
                    //for(int q=0;q<mutArray.length;q++) {
                    OneMutation curMut = i.next();
                    int w = findMutationIndex(resumeResults,curMut.resTypes);
                    if (w>=0)
                        bestScore = bestScore.max(resumeResults[w].score); //higher scores are better for Hybrid MinDEE-K*
                    else {
                        newIndex++;
                        newVector2.add(curMut);
                    }
                }
                mutSet = new TreeSet<OneMutation>();
                Iterator<OneMutation> q = newVector2.iterator();
                while(q.hasNext()) {
                    //for(int q=0; q<newArray2.length; q++){
                    mutSet.add(q.next());
                }
                //System.arraycopy(newArray2,0,mutArray,0,newIndex);
                System.out.println("Length of mutArray after removing already computed mutations: "+mutSet.size());
                if(mutSet.size() == 0) {
                    System.out.println("Length of mutArray is 0, so there are no mutations to compute.");
                    System.out.println("Quitting.....");
                    System.exit(0);
                }
            }
        }
        OneMutation[] mutArray = mutSet.toArray(new OneMutation[1]);
        //BigDecimal q_L = getLigPartFn(m,numInAS,ligType,eMatrixNameMin+".dat");

        MutationManager mutMan = new MutationManager(runName,mutArray,false);

        mutMan.setSaveTopConfs(saveTopConfs);
        mutMan.setPrintTopConfs(printTopConfs);
        mutMan.setNumTopConfs(numTopConfs);
        mutMan.setStrandMut(strandMut);
        mutMan.setStrandDefault(strandDefault);
        mutMan.setStrandPresent(strandPresent);
        mutMan.setStrandsPresent(strandsPresent);
        mutMan.setStrandLimits(strandLimits);
        mutMan.setMutableSpots(numberMutable);
        //mutMan.setLigPresent(ligPresent);
        mutMan.setAddOrigRots(addOrigRots);
        mutMan.setNumMutations(numMutations);
        mutMan.setarpFilenameMin(eMatrixNameMin+".dat");
        mutMan.setarpFilenameMax(eMatrixNameMax+".dat");
        mutMan.setDoMinimization(doMinimize);
        mutMan.setMinimizeBB(minimizeBB);
        mutMan.setDoBackrubs(doBackrubs);
        mutMan.setBackrubFile(backrubFile);
        mutMan.setRepeatSearch(repeatSearch);
        mutMan.setInitEw(initEw);
        mutMan.setGamma(gamma);
        mutMan.setEpsilon(epsilon);
        mutMan.setStericE(stericE);
        mutMan.setPruningE(pruningE);
        mutMan.setParams(sParams);
        mutMan.setStericThresh(stericThresh);
        mutMan.setSoftStericThresh(softStericThresh);
        mutMan.setComputeEVEnergy(true);
        mutMan.setCalculateVolumes(false);
        mutMan.setDistDepDielect(distDepDielect);
        mutMan.setDielectConst(dielectConst);
        mutMan.setDoDihedE(doDihedE);
        mutMan.setDoSolvationE(doSolvationE);
        mutMan.setSolvScale(solvScale);
        mutMan.setVdwMult(softvdwMultiplier);
        mutMan.setScaleInt(scaleInt);
        mutMan.setMaxIntScale(maxIntScale);
        mutMan.setRotamerLibrary(EnvironmentVars.aaRotLib);
        mutMan.setUseMaxKSconfs(useMaxKSconfs);
        mutMan.setMaxKSconfs(maxKSconfs);
        //DEEPer
        mutMan.setDoPerturbations(doPerturbations);
        mutMan.setPertFile(pertFile);
        mutMan.setMinimizePerts(minimizePerts);
        mutMan.setAddWTRot(addWTRot);
        mutMan.setIdealizeSC(Perturbation.idealizeSC);
        mutMan.setUseFlagsAStar(false);

        if (resumeSearch)
            mutMan.setBestScore(bestScore);	// Set the current best score from the partial results
        else
            mutMan.setBestScore(new BigDecimal("0.0")); //the initial best score is 0.0

        try {
            handleDoMPIMaster(mutMan,mutArray.length);
        } catch (Exception e) {
            System.out.println("ERROR: "+e);
            e.printStackTrace();
            System.exit(1);
        }

        System.out.println("DONE: K* computation");
    }

    /**
     * Computes the partition function for the ligand using the rotamers from the (ligand) rotamer library
     */
    /*private BigDecimal getLigPartFn(Molecule m, int numInAS, String ligType, String eMatrixNameMin){

    	float minMatrix[][][][][][] = (float [][][][][][])readObject(eMatrixNameMin);

    	int numRot = grl.getNumRotamers(m.strand[ligStrNum].residue[0].name);
    	if (numRot==0) //ALA or GLY
    		numRot = 1;

    	ExpFunction ef = new ExpFunction();

    	BigDecimal q_L = new BigDecimal("0.0");

    	for (int i=0; i<numRot; i++)
    		q_L = q_L.add(ef.exp(-minMatrix[numInAS][grl.getAARotamerIndex(ligType)][i][numInAS][0][0]/constRT));

    	System.out.println("Ligand partition function (double): "+q_L.doubleValue());

    	return q_L;
    }*/

    //Loads the mutation sequence list for Hybrid MinDEE-K*; computes a list if one cannot be loaded
    private Set<OneMutation> handleHybridKSLoadMutList (String mutFileName, int numMutable,
            Molecule m, int numComb, int residueMutatable[][], ParamSet sParams,int strandMut[][], String strandDefault[][],
            StrandRotamers[] strandRot, float targetVol,float volWindow,int strandsPresent, boolean strandPresent[]) {

        // Look for previous mutation file
        System.out.println();
        System.out.print("Looking for mutation list file ");
        Set<OneMutation> mutSet = loadMutationList(mutFileName,numMutable,false);

        if (mutSet == null) {

            // Create the mutation list with estimated energies
            mutSet = new TreeSet<OneMutation>();
            RotamerSearch rs = new RotamerSearch(m, numMutable, strandsPresent, hElect, hVDW, hSteric, true,
                                                 true, 0.0f, stericThresh, softStericThresh, distDepDielect, dielectConst,doDihedE,
                                                 doSolvationE,solvScale,softvdwMultiplier,grl,
                                                 false,null,false,false,false);//Not going to consider perturbations, etc. at this level of approximation

            //KER: load vol file for each rotamer library
            for (int strNum=0; strNum<rs.strandRot.length; strNum++) {
                //KER: Only load volfile if residue is allowed to mutate
                //if(rs.strandRot[strNum].rl.getNumAAallowed() > 1)
                rs.strandRot[strNum].rl.loadVolFile(); //load the rotamer volume file

            }

            int curNumSeq = 0;
            boolean valid;
            for(int i=0; i<numComb; i++) {
                valid = true;
                // Reset each amino acid type
                System.out.print("Starting mutation combination " + i + " ... ");
                for(int str=0; str<strandMut.length; str++)
                    rs.refreshStrand(str); // clears allowables and does some other stuff

                boolean addWT = (new Boolean((String)sParams.getValue("ADDWT", "true"))).booleanValue();
                if(!addWT)
                    checkWT(strandDefault, strandPresent, sParams);
                int molStrand = 0;
                int mutIndex = 0;
                for (int strNum=0; strNum<strandPresent.length; strNum++) {
                    if(strandPresent[strNum]) {
                        for (int k=0; k<strandMut[molStrand].length; k++) {
                            if (residueMutatable[i][mutIndex] == 1)
                                setAllowablesHelper(rs, sParams, addWT, strNum, molStrand, k, strandMut, strandDefault);
                            else {
                                valid = false;
                                String tempResAllow = (String)sParams.getValue("RESALLOWED"+strNum+"_"+k,"");
                                if(numTokens(tempResAllow) <= 0 && !addWT) {
                                    System.out.println("Error: resAllowed not set for strand");
                                    System.exit(1);
                                }
                                for(int q=0; q<numTokens(tempResAllow); q++) {
                                    if(getToken(tempResAllow,q+1).equalsIgnoreCase(strandDefault[molStrand][k])) {
                                        rs.setAllowable(strandMut[molStrand][k],strandDefault[molStrand][k],molStrand); //the default type is set last
                                        valid = true;
                                    }
                                }
                                if(addWT) {
                                    rs.setAllowable(strandMut[molStrand][k],strandDefault[molStrand][k],molStrand); //the default type is set last
                                    valid = true;
                                }
                            }
                            mutIndex++;
                        }
                        molStrand++;
                    }
                }

                // Perform simple mutation search for this set of mutatable residues
                if(valid) {
                    curNumSeq = rs.simpleMasterMutationSearch(strandMut,numMutable,
                                curNumSeq,mutSet,targetVol-volWindow,
                                targetVol+volWindow);
                }
                System.out.println("finished");
            }

            System.out.println("Sequences remaining after volume filter "+curNumSeq);

            // We now have all the mutations in mutArray, collapse the mutArray
            //  to the actual number of mutations we have.

            //KER: mutArray is now a set so we don't need to reallocate a smaller array
            //OneMutation newArray[] = new OneMutation[curNumSeq];
            //System.out.println("Allocated newArray");
            //System.out.println("Initial Length of mutArray: "+mutArray.length);
            //System.arraycopy(mutArray,0,newArray,0,curNumSeq);
            //mutArray = newArray;
            System.out.println("Trimmed Length of mutArray: "+mutSet.size());

            //KER: mutArray is now a set so there should be no duplicates to begin with
            /*System.out.print("Removing duplicates...");
            mutArray = removeDuplicates(mutArray);
            System.out.println("done");*/

            System.out.println(mutSet.size()+" unique mutation sequences found in volume range "+(targetVol-volWindow)+" to "+(targetVol+volWindow));
            BigInteger numConfs = BigInteger.ZERO;
            Iterator<OneMutation> i = mutSet.iterator();
            while(i.hasNext()) {
                OneMutation tmp = i.next();
                numConfs = numConfs.add(tmp.numConfUB.add(tmp.numConfB));
            }
            System.out.println("Total number of conformations (bound and unbound) for all sequences: "+numConfs);// Save mutation list
            saveMutationList(mutSet,mutFileName,false);
        }


        // Sort the mutation list
        // System.out.print("Sorting mutation list ... ");
        //KER: mutSet is a TreeSet which is sorted so no need to sort
        //RyanQuickSort rqs = new RyanQuickSort();
        //rqs.Sort(mutArray);
        //rqs = null;
        // System.out.println("done");

        return mutSet;
    }


    /**
     * Reads the results of a partially completed run into an array of CommucObj. The MutationManager then queries
     * this array before sending out a task.
    */
    public OneMutation[] readResumeFile(OneMutation resumeResults[], String resumeFilename, int numMutable, boolean distrDACS, boolean PEMcomp, int initDepth) {

        BufferedReader bufread = null;
        try {
            File file = new File(resumeFilename);
            FileReader fr = new FileReader(file);
            bufread = new BufferedReader(fr);
        } catch (FileNotFoundException e) {
            System.out.println("ERROR: Resume File Not Found");
            return(null);
        }

        boolean done = false;
        String str = null;
        int resultNum = 0;

        while (!done) {
            try {
                str = bufread.readLine();
            } catch ( Exception e ) {
                System.out.println("ERROR: An error occurred while reading input");
                System.exit(0);
            }

            if (str == null) // stop if we've reached EOF
                done = true;
            else if (getToken(str,1).equalsIgnoreCase("Completed")) {
                if (PEMcomp) {//PEM computation
                    resumeResults[resultNum].resMut = new int[numMutable];
                    for(int q=0; q<numMutable; q++)
                        resumeResults[resultNum].resMut[q] = new Integer(getToken(str,8+q)).intValue();
                    resumeResults[resultNum].flagMutType = getToken(str,8+numMutable);
                } else { //mutation search
                    if (!distrDACS) { //Hybrid-K* or MinDEE/A* resume
                        resumeResults[resultNum].score = new BigDecimal(getToken(str,5));
                        resumeResults[resultNum].resTypes = new String[numMutable];
                        for(int q=0; q<numMutable; q++) {
                            resumeResults[resultNum].resTypes[q] = getToken(str,21+q);
                        }
                    } else { //distributed DACS resume
                        resumeResults[resultNum].mutNum = new Integer(getToken(str,3)).intValue();
                        resumeResults[resultNum].score = new BigDecimal(getToken(str,7));
                        resumeResults[resultNum].resMut = new int[initDepth];
                        for(int q=0; q<initDepth; q++)
                            resumeResults[resultNum].resMut[q] = new Integer(getToken(str,9+q));
                    }
                }
                resultNum++;
            }
        }

        // We're done reading them in
        try {
            bufread.close();
        } catch(Exception e) {
        }

        // Resize completed mutation array
        OneMutation temp[] = new OneMutation[resultNum];
        System.arraycopy(resumeResults,0,temp,0,resultNum);
        resumeResults = temp;
        return(resumeResults);
    }


    // Finds the index of the mutation in resumeResults with the same
    //  mutation sequence as the targetMutation. If none are found, -1
    //  is returned.
    public int findMutationIndex(OneMutation resumeResults[],
                                 String targetMutation[]) {

        for(int q=0; q<resumeResults.length; q++) {
            if (resumeResults[q].isSame(targetMutation))
                return(q);
        }
        return(-1);
    }


    // Attempts to read a list of mutations from file
    public Set<OneMutation> loadMutationList(String fName, int numMutable, boolean PEMcomp) {

        BufferedReader bufread = null;
        try {
            File file = new File(fName);
            FileReader fr = new FileReader(file);
            bufread = new BufferedReader(fr);
        } catch (FileNotFoundException e) {
            System.out.println(" ... no mutation list file found. Computing one.");
            return(null);
        }

        boolean done = false;
        String str = null;
        int resultNum = 0;
        Set<OneMutation> mutList = new TreeSet<OneMutation>();

        while (!done) {
            try {
                str = bufread.readLine();
            } catch ( Exception e ) {
                System.out.println("ERROR: An error occurred while reading input");
                System.exit(0);
            }

            if (str == null) // stop if we've reached EOF
                done = true;
            else {
                if (PEMcomp) {//PEM computation
                    OneMutation tmp = new OneMutation();
                    //mutList[resultNum] = new OneMutation();
                    //mutList[resultNum].resMut = new int[numMutable];
                    tmp.resMut = new int[numMutable];
                    for(int q=0; q<numMutable; q++)
                        tmp.resMut[q] = new Integer(getToken(str,1+q)).intValue();

                    tmp.flagMutType = getToken(str,1+numMutable);
                    mutList.add(tmp);
                } else { //mutation search
                    OneMutation tmp = new OneMutation();
                    //mutList[resultNum] = new OneMutation();
                    tmp.score = new BigDecimal(getToken(str,1));
                    tmp.vol = new Float(getToken(str,2)).floatValue();
                    tmp.resTypes = new String[numMutable];
                    for(int q=0; q<numMutable; q++) {
                        tmp.resTypes[q] = getToken(str,3+q);
                    }
                    mutList.add(tmp);
                }

                resultNum++;
                /*if (resultNum >= mutList.length){
                	OneMutation newArray[] = new OneMutation[mutList.length+1000];
                	System.arraycopy(mutList,0,newArray,0,resultNum);
                	mutList = newArray;
                }*/
            }
        }

        // We're done reading them in
        try {
            bufread.close();
        } catch(Exception e) {
        }

        // Resize completed mutation array
        /*OneMutation temp[] = new OneMutation[resultNum];
        System.arraycopy(mutList,0,temp,0,resultNum);
        mutList = temp;*/
        System.out.println(" ... read "+mutList.size()+" mutations from mutation list "+fName);
        return(mutList);
    }

    // Saves the list of mutations so that a PEM computation/mutation search
    //  doesn't need to recompute these during a resume. Thus
    //  the resume can go more quickly.
    public void saveMutationList(Set<OneMutation> mutList, String fName, boolean PEMcomp) {

        if (mutList.size() == 0)
            return;

        int numInAS = 0;
        Iterator<OneMutation> i = mutList.iterator();
        if (PEMcomp)
            numInAS = i.next().resMut.length;
        else
            numInAS = i.next().resTypes.length;

        PrintStream printStream = setupOutputFile(fName);
        i = mutList.iterator();
        while(i.hasNext()) {
            //for(int q=0;q<mutList.length;q++) {
            OneMutation tmp = i.next();
            if (PEMcomp) {//PEM computation
                for(int w=0; w<numInAS; w++) {
                    printStream.print(" "+tmp.resMut[w]);
                }
                printStream.print(" "+tmp.flagMutType);
                printStream.println();
            } else { //mutation search
                printStream.print(tmp.score + " " + tmp.vol);
                for(int w=0; w<numInAS; w++) {
                    printStream.print(" "+tmp.resTypes[w]);
                }
                printStream.println();
            }
        }
        printStream.close();
    }

    //Removes duplicate mutations (for which the mutation sequence is the same) from a given list
    public OneMutation [] removeDuplicates(OneMutation mutArray[]) {

        //First, sort the list alphabetically, according to the mutation sequence
        RyanQuickSort rqs = new RyanQuickSort();
        rqs.Sort(mutArray);
        rqs = null;

        //Copy mutArray into nArray, excluding duplicate entries
        OneMutation nArray[] = new OneMutation[mutArray.length];

        //Copy the first element
        nArray[0] = mutArray[0];
        int nAIndex = 1;

        //Compare each mutation with the previous one in the list
        for (int i=1; i<mutArray.length; i++) { //for each mutation
            if (!(mutArray[i].isSame(mutArray[i-1].resTypes))) { //different sequence
                nArray[nAIndex] = mutArray[i];
                nAIndex++;
            }
        }

        mutArray = new OneMutation[nAIndex];
        System.arraycopy(nArray,0,mutArray,0,mutArray.length);

        return mutArray;//return the reduced list
    }

    // Mutation search Slave function
    public CommucObj handleKSSlave(CommucObj cObj) {
        EnvironmentVars.aaRotLib = cObj.rl;
        if (cObj.PEMcomp) { //PEM computation
            if (!cObj.entropyComp) //PEM computation
                cObj = handleComputeAllPairwiseRotamerEnergiesSlave(cObj);
            else //entropy E matrix computation
                cObj = handleDoResEntropySlave(cObj);
        } else { //distributed mutation search
            if (cObj.distrDACS) { //running distributed DACS
                cObj = doDistrDACSSlave(cObj);
            } else if (cObj.distrDEE) { //running distributed DEE
                cObj = doDistrDEESlave(cObj);
            } else { //running Hybrid MinDEE-K*
                cObj = hybridKScompute(cObj);
            }
        }
        return cObj;
    }

    /**
     * Handles the computation of the K* score for a single mutation sequence with the target ligand.
     * The 'cObj' parameter contains the mutation search input distributed by the main processor.
     * Returns the results of the computation to the main processor.
     */
    private CommucObj hybridKScompute(CommucObj cObj) {

        //System.out.println("Start of hybridKScompute");

        Molecule m;


        //KER: We only need to find the partition function for each unbound
        //KER: strand and then the whole complex
        //int[][] allCombos = generateCombos(cObj.strandMut.length);
        boolean[][] allCombos = new boolean[cObj.strandMut.length+1][];
        for(int com=0; com<allCombos.length; com++) {
            allCombos[com] = new boolean[cObj.strandMut.length];
            for(int i=0; i<allCombos[com].length; i++)
                if(com==allCombos.length-1)
                    allCombos[com][i] = true;
                else if(i==com)
                    allCombos[com][i] = true;
                else
                    allCombos[com][i] = false;
        }

        cObj.numComplexes = allCombos.length;
        cObj.repeatEW = new boolean[allCombos.length];
        cObj.allPruned = new boolean[allCombos.length];
        cObj.searchNumConfsTotal = new int[allCombos.length];
        cObj.searchNumConfsPrunedByE = new int[allCombos.length];
        cObj.searchNumConfsPrunedByS = new int[allCombos.length];
        cObj.searchNumConfsEvaluated = new int[allCombos.length];
        cObj.searchNumConfsLeft = new int[allCombos.length];
        cObj.searchNumPrunedMinDEE = new int[allCombos.length];
        cObj.searchBestEnergyFound = new float[allCombos.length];
        cObj.q = new BigDecimal[allCombos.length];
        cObj.bestEMin = new double[allCombos.length];
        cObj.bestE = new double[allCombos.length];
        cObj.q_Time = new int[allCombos.length];
        //cObj.effEpsilon = new double[allCombos.length];

        System.out.print("## CurMut: "+cObj.curMut+" Starting Sequence: ");
        for(int i=0; i<cObj.mutableSpots; i++)
            System.out.print(" "+cObj.currentMutation[i]);
        System.out.println(" &&");

        //KER: Run through all of the partition function calculations
        for(int runNum = 0; runNum<allCombos.length; runNum++) {
            long startTime = System.currentTimeMillis();


            boolean strandPresent[] = allCombos[runNum];

            String unboundStr = Integer.toString(runNum);
            int strandsPresent = 1;
            boolean notFullComplex = true;
            if(runNum == allCombos.length-1) {
                notFullComplex = false;
                strandsPresent = allCombos[runNum].length;
                System.out.println("Calculating partition function for entire complex");
            } else {
                System.out.println("Calculating partition function for strand: "+runNum);
            }

            ParamSet params = null;
            String minEmatrixFile = null;
            String maxEmatrixFile = null;
            //TODO:fix this hack....
            if ( notFullComplex && ((new Boolean((String)cObj.params.getValue("USEUNBOUNDSTRUCT"+unboundStr,"false"))).booleanValue()) ) { //use a different input PDB structure for the unbound case
                params = new ParamSet();
                params.setParamsValues(cObj.params.getParams(), cObj.params.getValues(), cObj.params.getCurNum());
                params.setValue("PDBNAME",params.getValue("UNBOUNDPDBNAME"+unboundStr));
                //params.setValue("PDBLIGNUM","-1");
                minEmatrixFile = cObj.arpFilenameMin;
                maxEmatrixFile = cObj.arpFilenameMax;
                minEmatrixFile = minEmatrixFile.replace(".dat", "_"+unboundStr+".dat");
                maxEmatrixFile = maxEmatrixFile.replace(".dat", "_"+unboundStr+".dat");
                //minEmatrixFile = params.getValue("MINENERGYMATRIXNAMEUNBOUND"+unboundStr)+".dat";
                //maxEmatrixFile = params.getValue("MAXENERGYMATRIXNAMEUNBOUND"+unboundStr)+".dat";
            } else { //a single input PDB structure is used for the bound and unbound computations
                params = cObj.params;
                minEmatrixFile = cObj.arpFilenameMin;
                maxEmatrixFile = cObj.arpFilenameMax;
                if(notFullComplex) {
                    minEmatrixFile = minEmatrixFile.replace(".dat", "_"+unboundStr+".dat");
                    maxEmatrixFile = maxEmatrixFile.replace(".dat", "_"+unboundStr+".dat");
                } else {
                    minEmatrixFile = minEmatrixFile.replace(".dat", "_COM.dat");
                    maxEmatrixFile = maxEmatrixFile.replace(".dat", "_COM.dat");
                }
            }

            if(cObj.doPerturbations)
                Perturbation.idealizeSC = cObj.idealizeSC;

            //Setup the molecule system
            MolParameters mp = new MolParameters();
            mp.numOfStrands = cObj.numberOfStrands;
            mp.strandLimits = cObj.strandLimits;
            mp.strandPresent = strandPresent;
            mp.strandsPresent = strandsPresent;

            //KER: to match up with the Ematrix code
            //we index into mols the same way we do there
            int curStr = runNum+1;
            if(curStr == allCombos.length)
                curStr = 0;
            if(mols[curStr] == null) {
                mp.m = new Molecule();
                mp.m = setupMolSystem(mp.m,cObj.params,mp.strandPresent,mp.strandLimits);
                mols[curStr] = mp.m;
                rotLibs[curStr] = new RotamerLibrary[strandsPresent];
                for(int i=0; i<strandsPresent; i++)
                    rotLibs[curStr][i] = grl[i];
            } else {
                mp.m = mols[curStr];
                grl = rotLibs[curStr];
            }

            loadMutationParams(cObj.params, mp);
            mp.numberMutable = getNumberMutable(mp.strandMut);
            System.out.println("NumMutable "+mp.numberMutable);

            //KER: Create local variables so I don't have to add mp to everything...
            m = mp.m;
            int numberMutable = mp.numberMutable;
            int[][] strandMut = mp.strandMut;

            double minEBound = 0.0;
            BigInteger numConfsPrunedMinDEESteric = null;

            RotamerSearch rs = new RotamerSearch(m,numberMutable,strandsPresent, hElect, hVDW, hSteric, true,
                                                 true, cObj.epsilon, cObj.stericThresh, cObj.softStericThresh, cObj.distDepDielect,
                                                 cObj.dielectConst,cObj.doDihedE,cObj.doSolvationE,cObj.solvScale,cObj.vdwMult,grl,
                                                 cObj.doPerturbations,cObj.pertFile, cObj.minimizePerts, false, false);

            System.out.print("Loading precomputed min energy matrix...");
            rs.loadPairwiseEnergyMatrices(minEmatrixFile,true);
            System.out.println("done");

            if (cObj.doMinimization) { //MinDEE, so load the max matrix
                System.out.print("MinDEE: Loading precomputed max energy matrix...");
                rs.loadPairwiseEnergyMatrices(maxEmatrixFile,false);
                System.out.println("done");
            }

            rs.initMutRes2Str(strandMut);
            //KER: currentMutOffset is needed because if certain strands aren't present we need
            //to index differently into the current mutable array. For example, in a str
            //lig system if we're only looking at the strand all the mutables in the str
            //need to be skipped over
            int currentMutOffset[] = null;

            /********Shrink the Matrices Such that they are the right size for the current mutation****/
            currentMutOffset = getCurrentMutOffset(cObj.strandMut,strandMut,strandPresent,cObj.strandPresent);
            String[] currentMutation = shrinkCurrentMutation(cObj.currentMutation, cObj.strandMut, strandMut, cObj.strandPresent, strandPresent);
            //Set the allowable AAs for each AS residue
            int ctr=0;
            for(int str=0; str<strandMut.length; str++) {
                for(int i=0; i<strandMut[str].length; i++) {
                    rs.setAllowable(strandMut[str][i], currentMutation[ctr], str);
                    ctr++;
                }
            }

            if(cObj.doPerturbations)
                rs.setupRCs(cObj.addWTRot);

            //Initially, no rotamers have been pruned
            PrunedRotamers<Boolean> prunedRotAtRes = new PrunedRotamers<Boolean>(mp.numberMutable, strandMut, rs, false);
            //boolean prunedRotAtRes[] = new boolean[numberMutable*cObj.numTotalRotamers];
            //for (int i=0; i<prunedRotAtRes.length; i++)
            //	prunedRotAtRes[i] = false;


            //Prune all rotamers that are incompatible with the template (intra E + re-to-template E >= stericE)
            prunedRotAtRes = rs.DoPruneStericTemplate(numberMutable, strandMut, prunedRotAtRes, cObj.stericE);

            //Perform the DEE pruning
            if (cObj.doMinimization) //compute the MinDEE interval terms
                rs.doCompMinDEEIntervals(numberMutable, strandMut,
                                         prunedRotAtRes, cObj.scaleInt, cObj.maxIntScale);

            //Perform Goldstein Pruning
            prunedRotAtRes = rs.DoDEEGoldstein(numberMutable, strandMut, cObj.initEw, prunedRotAtRes,
                                               cObj.doMinimization, false,	cObj.minimizeBB, cObj.typeDep, false, 0.0f);

            //Prune with MinBounds (last parameter is false)
            prunedRotAtRes = rs.DoMinBounds(numberMutable,strandMut,
                                            cObj.pruningE,prunedRotAtRes,cObj.initEw, false, false);

            //Compute the Ec value and prunedIsSteric[] (last parameter is true)
            rs.DoMinBounds(numberMutable,strandMut,cObj.pruningE,prunedRotAtRes,cObj.initEw, false, true);


            //TODO: Fix usingInitialBest and setting initialBest
            boolean usingInitialBest = (!notFullComplex);
            BigDecimal initialBest = (new BigDecimal("0.0"));

            //Inter-mutation pruning
            if (usingInitialBest) {
                initialBest = cObj.bestScore.multiply(new BigDecimal(cObj.gamma * cObj.epsilon));
                for(int i=0; i<allCombos.length-1; i++)
                    initialBest =  initialBest.multiply(cObj.q[i]);
            }

            rs.slaveDoRotamerSearch(runNum, cObj.computeEVEnergy,cObj.doMinimization,numberMutable,
                                    strandMut,usingInitialBest,initialBest,cObj,cObj.minimizeBB,false,cObj.pdbName,cObj.doBackrubs,cObj.backrubFile,
                                    cObj.saveTopConfs, cObj.printTopConfs, cObj.numTopConfs, cObj.curMut,
                                    cObj.useMaxKSconfs, cObj.maxKSconfs);

            if ((cObj.repeatSearch)&&(rs.repeatSearch)) { //the desired accuracy was not achieved, so repeat the search: the setup is already done

                rs.repeatSearch = false; //reset the flag
                cObj.repeatEW[runNum] = true;

                rs.slaveDoRotamerSearch(runNum, cObj.computeEVEnergy,cObj.doMinimization,numberMutable,
                                        strandMut,usingInitialBest,initialBest,cObj,cObj.minimizeBB,false,null,cObj.doBackrubs,
                                        cObj.backrubFile, cObj.saveTopConfs, cObj.printTopConfs, cObj.numTopConfs, cObj.curMut,
                                        cObj.useMaxKSconfs, cObj.maxKSconfs);
            }

            long stopTime = System.currentTimeMillis();
            cObj.q_Time[runNum] = Math.round((stopTime - startTime) / 1000.0f);
        } // end for(runNum)

        System.out.print("## CurMut: "+cObj.curMut+" Finished Sequence: ");
        for(int i=0; i<cObj.mutableSpots; i++)
            System.out.print(" "+cObj.currentMutation[i]);
        System.out.println(" &&");

        return cObj;
    }

    //Load the pairwise energy matrices; if not computed, compute, and the load
    //KER: strandPresent represents the strand that is present for this matrix (complex is -1)
    private void loadPairwiseEnergyMatrices(ParamSet sParams, RotamerSearch rs, String minMatrixFile, boolean doMinimize, String maxMatrixFile,
                                            int strandPresent) {
        //StrandNumbers are appended to PEM matrix, Complex is denoted with "_COM"
        // PGC
        String runName = sParams.getValue("RUNNAME");
        String suffix = "_" + strandPresent+".dat";
        if(strandPresent == COMPLEX )
            suffix = "_COM.dat";
        rs.loadPairwiseEnergyMatrices(minMatrixFile + suffix,true);
        if (doMinimize)
            rs.loadPairwiseEnergyMatrices(maxMatrixFile + suffix,false);

        if ( (rs.getMinMatrix()==null) || ( doMinimize && rs.getMaxMatrix()==null ) ) { //at least one of the matrices not computed, so compute

            System.out.println("Precomputed energy matrices not available..");

            long startTime = System.currentTimeMillis();
            ParamSet newParams = new ParamSet(); //create a new parameter set, just for the unbound-case matrix computation; sParams must not be changed here
            newParams.setParamsValues(sParams.getParams(), sParams.getValues(), sParams.getCurNum());
            newParams.setValue("MINENERGYMATRIXNAME", sParams.getValue("MINENERGYMATRIXNAME",runName+"minM")+suffix);
            newParams.setValue("MAXENERGYMATRIXNAME", sParams.getValue("MAXENERGYMATRIXNAME",runName+"maxM")+suffix);
            newParams.setValue("ONLYSINGLESTRAND", new Integer(strandPresent).toString());
            handleComputeAllPairwiseRotamerEnergiesMaster(newParams);

            long stopTime = System.currentTimeMillis();
            System.out.println("PEM execution time: "+((stopTime-startTime)/(60.0*1000.0)));
            System.out.println("DONE: Pairwise energy matrix precomputation");

            rs.loadPairwiseEnergyMatrices(minMatrixFile+suffix,true);
            if (doMinimize)
                rs.loadPairwiseEnergyMatrices(maxMatrixFile+suffix,false);
        }
    }

    //Load the pairwise energy matrices; if not computed, compute, and the load
    //KER: strandPresent represents the strand that is present for this matrix (complex is -1)
    private void loadUnboundPairwiseEnergyMatrices(ParamSet sParams, RotamerSearch rs, String minMatrixFile, boolean doMinimize, String maxMatrixFile,
            int unboundStrand) {
        //StrandNumbers are appended to PEM matrix, Complex is denoted with "_COM"
        // PGC
        String runName = sParams.getValue("RUNNAME");

        String suffix = "_" + unboundStrand+".dat";

        rs.loadPairwiseEnergyMatrices(minMatrixFile + suffix,true);
        if (doMinimize)
            rs.loadPairwiseEnergyMatrices(maxMatrixFile + suffix,false);

        if ( (rs.getMinMatrix()==null) || ( doMinimize && rs.getMaxMatrix()==null ) ) { //at least one of the matrices not computed, so compute

            System.out.println("Precomputed energy matrices not available..");

            long startTime = System.currentTimeMillis();
            ParamSet newParams = new ParamSet(); //create a new parameter set, just for the unbound-case matrix computation; sParams must not be changed here
            newParams.setParamsValues(sParams.getParams(), sParams.getValues(), sParams.getCurNum());
            newParams.setValue("MINENERGYMATRIXNAME", sParams.getValue("MINENERGYMATRIXNAME",runName+"minM")+suffix);
            newParams.setValue("MAXENERGYMATRIXNAME", sParams.getValue("MAXENERGYMATRIXNAME",runName+"maxM")+suffix);
            newParams.setValue("ONLYSINGLESTRAND", "-1");
            handleComputeAllPairwiseRotamerEnergiesMaster(newParams);

            long stopTime = System.currentTimeMillis();
            System.out.println("PEM execution time: "+((stopTime-startTime)/(60.0*1000.0)));
            System.out.println("DONE: Pairwise energy matrix precomputation");

            rs.loadPairwiseEnergyMatrices(minMatrixFile+suffix,true);
            if (doMinimize)
                rs.loadPairwiseEnergyMatrices(maxMatrixFile+suffix,false);
        }
    }


/////////////////////////////////////////////////////////////////////////
// MIN and MAX pairwise energy matrices computation
/////////////////////////////////////////////////////////////////////////

    //KER: wrapper for computing the pairwise energy matrices
    public void handleComputeAllPairwiseRotamerEnergies(String s) {

        ParamSet sParams = new ParamSet();
        sParams.addParamsFromFile(getToken(s,2)); //read system parameters
        sParams.addParamsFromFile(getToken(s,3)); //read mutation search parameters

        // PGC
        String runName = sParams.getValue("RUNNAME");
        int curStrForMatrix = new Integer(sParams.getValue("ONLYSINGLESTRAND","-1"));

        String suffix = "_" + curStrForMatrix + ".dat";
        if(curStrForMatrix == COMPLEX )
            suffix = "_COM.dat";

        sParams.setValue("MINENERGYMATRIXNAME", sParams.getValue("MINENERGYMATRIXNAME",runName+"minM")+suffix);
        sParams.setValue("MAXENERGYMATRIXNAME", sParams.getValue("MAXENERGYMATRIXNAME",runName+"maxM")+suffix);

        handleComputeAllPairwiseRotamerEnergiesMaster(sParams);
    }


    /**
     * This function sets up the computation for all min and max pairwise rotamer interaction energies and stores
     * these energies into user-specified precomputed energy matrices. The computation is distributed to the available processors.
     * For each rotamer, for each allowed amino acid type, for each flexible residue position, the following energies
     * are computed: intra-rotamer, rotamer-to-template, rotamer-rotamer, and template-only energies. If a ligand is
     * present, all energies involving the ligand are also computed.
     */
    public void handleComputeAllPairwiseRotamerEnergiesMaster(ParamSet sParams) {

        int numMutations = 2; //pairwise energies are computed
        boolean doMinimize = (new Boolean((String)sParams.getValue("DOMINIMIZE", "false"))).booleanValue();
        boolean minimizeBB = (new Boolean((String)sParams.getValue("MINIMIZEBB", "false"))).booleanValue();
        boolean doBackrubs = (new Boolean((String)sParams.getValue("DOBACKRUBS", "false"))).booleanValue();
        String backrubFile ="";
        if(doBackrubs) {
            backrubFile = (String)sParams.getValue("BACKRUBFILE");
        }
        String runName = (String)sParams.getValue("RUNNAME");
        String minEMatrixName = (String)sParams.getValue("MINENERGYMATRIXNAME",runName+"minM");
        String maxEMatrixName = (String)sParams.getValue("MAXENERGYMATRIXNAME",runName+"maxM");
        String eRefMatrix = (String)(sParams.getValue("EREFMATRIXNAME","Eref"));
        boolean templateAlwaysOn = (new Boolean((String)sParams.getValue("TEMPLATEALWAYSON","false"))).booleanValue();

        int curStrForMatrix = (new Integer((String)sParams.getValue("ONLYSINGLESTRAND","-1"))).intValue();
        boolean typeDep = (new Boolean((String)sParams.getValue("TYPEDEP","false"))).booleanValue();


        boolean doPerturbations = (new Boolean((String)sParams.getValue("DOPERTURBATIONS","false"))).booleanValue();//Triggers DEEPer
        String pertFile = (String)sParams.getValue("PERTURBATIONFILE","defaultPerturbationFileName.pert");//Input file giving perturbation information
        boolean minimizePerts = (new Boolean((String)sParams.getValue("MINIMIZEPERTURBATIONS","false"))).booleanValue();//Allow continuous minimization with respect to perturbation parameters
        Perturbation.idealizeSC = (new Boolean((String)sParams.getValue("IDEALIZESIDECHAINS","true"))).booleanValue();


        boolean addWTRotsSomehow = (new Boolean((String)sParams.getValue("ADDWTROTS","false"))).booleanValue();
        boolean addOrigRots = false, addWTRot = false;
        if( addWTRotsSomehow ) {
            if(doPerturbations)//DEEPer allows adding WT rotamers only to the positions they came from, so do that if possible
                addWTRot = true;
            else//Otherwise just add them to the rotamer library
                addOrigRots = true;
        }


        if (!doMinimize) //no minimization
            minimizeBB = false;
        if (!minimizeBB) //not backbone minimization
            doBackrubs = false;

        //Setup the molecule system
        MolParameters mp = loadMolecule(sParams, curStrForMatrix);
        Molecule m = mp.m;
        int numberMutable = mp.numberMutable;
        int strandsPresent = mp.strandsPresent;
        String[][] strandLimits = mp.strandLimits;
        boolean[] strandPresent = mp.strandPresent;
        int[][] strandMut = mp.strandMut;
        String[][] strandDefault = mp.strandDefault;

        if(addOrigRots)
            RotamerLibrary.addOrigRots(strandMut, EnvironmentVars.aaRotLib, m);

        System.out.println("Run Name: "+runName);
        System.out.println("Precomputed Minimum Energy Matrix: "+minEMatrixName);
        System.out.println("Precomputed Maximum Energy Matrix: "+maxEMatrixName);
        System.out.println("Num Residues Allowed to Mutate: "+numMutations);


        System.out.println("Computing _All_ Rotamer-Rotamer Energies");

        System.out.println("Starting minimum and maximum bound energy computation");

        RotamerSearch rs = new RotamerSearch(m, numberMutable, strandsPresent, hElect, hVDW, hSteric, true,
                                             true, 0.0f, stericThresh, softStericThresh, distDepDielect, dielectConst, doDihedE,
                                             doSolvationE, solvScale, softvdwMultiplier, grl, doPerturbations, pertFile, minimizePerts, false, false);

        rs.initMutRes2Str(strandMut);

        int resMut[] = new int[numberMutable];
        for (int i=0; i<resMut.length; i++)
            resMut[i] = 1;

        //System.out.println("Beginning setAllowables");
        //Set the allowable AAs for each AS residue
        boolean addWT = (new Boolean((String)sParams.getValue("ADDWT", "true"))).booleanValue();
        if(!addWT)
            checkWT(strandDefault, strandPresent, sParams);
        int molStrand = 0;
        for (int strNum=0; strNum<strandPresent.length; strNum++) {
            if(strandPresent[strNum]) {
                for (int k=0; k<strandMut[molStrand].length; k++) {
                    setAllowablesHelper(rs, sParams, addWT, strNum, molStrand, k, strandMut, strandDefault);
                }
                molStrand++;
            }
        }

        if(doPerturbations)
            rs.setupRCs(addWTRot);

        //initialize the pairwise energy matrices (full initialization - for all residues in residueMap[], the ligand, and the template)
        PairwiseEnergyMatrix mutationEnergiesMin = new PairwiseEnergyMatrix(numberMutable,resMut,strandMut,rs,true,true,true);
        PairwiseEnergyMatrix mutationEnergiesMax = null;
        if(doMinimize)
            mutationEnergiesMax = mutationEnergiesMin.copy();
        float eRef[][] = new float[numberMutable][];
        int ctr = 0;
        for (int str=0; str<strandMut.length; str++) {
            for (int j=0; j<strandMut[str].length; j++) {
                eRef[ctr] = new float[rs.strandRot[str].rl.getNumAAallowed()];
                ctr++;
            }
        }


        OneMutation mutArray[] = getMutArrayPairEcomp(numberMutable,minimizeBB);

        //Sort the mutation list
        System.out.print("Sorting mutation list ... ");
        //RyanQuickSort rqs = new RyanQuickSort();
        Arrays.sort(mutArray);
        for(int i=0; i<mutArray.length; i++)
            mutArray[i].mutNum = i;
        //rqs = null;
        System.out.println("done");

        MutationManager mutMan = new MutationManager(null,mutArray,true);
        mutMan.setTemplateAlwaysOn(templateAlwaysOn);
        mutMan.setStrandMut(strandMut);
        mutMan.setStrandDefault(strandDefault);
        mutMan.setStrandPresent(strandPresent);
        mutMan.setStrandLimits(strandLimits);
        mutMan.setStrandsPresent(strandsPresent);
        mutMan.setAddOrigRots(addOrigRots);
        mutMan.setTypeDep(typeDep);
        mutMan.setMutableSpots(numberMutable);
        mutMan.setarpFilenameMin(minEMatrixName);
        mutMan.setPairEMatrixMin(mutationEnergiesMin);
        mutMan.setEref(eRef);
        mutMan.setErefMatrixName(eRefMatrix);
        if (doMinimize) {
            mutMan.setarpFilenameMax(maxEMatrixName);
            mutMan.setPairEMatrixMax(mutationEnergiesMax);
        }
        mutMan.setParams(sParams);
        mutMan.setStericThresh(stericThresh);
        mutMan.setSoftStericThresh(softStericThresh);
        mutMan.setComputeEVEnergy(true);
        mutMan.setDoMinimization(doMinimize);
        mutMan.setMinimizeBB(minimizeBB);
        mutMan.setDoBackrubs(doBackrubs);
        mutMan.setBackrubFile(backrubFile);
        mutMan.setCalculateVolumes(false);
        mutMan.setDistDepDielect(distDepDielect);
        mutMan.setDielectConst(dielectConst);
        mutMan.setDoDihedE(doDihedE);
        mutMan.setDoSolvationE(doSolvationE);
        mutMan.setSolvScale(solvScale);
        mutMan.setVdwMult(softvdwMultiplier);
        mutMan.setRotamerLibrary(EnvironmentVars.aaRotLib);
        mutMan.setcurStrForMatrix(curStrForMatrix);

        mutMan.setDoPerturbations(doPerturbations);
        mutMan.setMinimizePerts(minimizePerts);
        mutMan.setPertFile(pertFile);
        mutMan.setIdealizeSC(Perturbation.idealizeSC);
        mutMan.setAddWTRot(addWTRot);


        try {
            handleDoMPIMaster(mutMan,mutArray.length);
        } catch (Exception e) {
            System.out.println("ERROR: "+e);
            e.printStackTrace();
            System.exit(1);
        }

        outputObject(mutMan.getMinEmatrix().eMatrix,minEMatrixName);
        if (doMinimize)
            outputObject(mutMan.getMaxEmatrix().eMatrix,maxEMatrixName);

        System.out.println("DONE: Pairwise energy matrix precomputation..");
    }

    /**
     * Computes a specific part of the pairwise energy matrices, as specified by the parameters in the 'cObj' parameter,
     * distributed by the main processor. Returns the results of the computation to the main processor.
     */
    public CommucObj handleComputeAllPairwiseRotamerEnergiesSlave(CommucObj cObj) {

        long startTime = System.currentTimeMillis();

        //Setup the molecule system
        Molecule m = new Molecule();
        //CurStrForMatrix starts at -1 so push everything up one;
        int curComplex = cObj.curStrForMatrix+1;
        if(mols[curComplex] == null || cObj.doPerturbations ) { //The molecule must be reloaded for DEEPer because the perturbation states are unknown and we might need WT rotamers
            m = setupMolSystem(m,cObj.params,cObj.strandPresent,cObj.strandLimits);
            mols[curComplex] = m;
            rotLibs[curComplex] = grl;
        } else {
            m = mols[curComplex];
            grl = rotLibs[curComplex];
        }

        RotamerSearch rs = new RotamerSearch(m,cObj.mutableSpots, cObj.strandsPresent, hElect, hVDW, hSteric, true,
                                             true, cObj.epsilon, cObj.stericThresh, cObj.softStericThresh, cObj.distDepDielect, cObj.dielectConst,
                                             cObj.doDihedE,cObj.doSolvationE,cObj.solvScale,cObj.vdwMult,grl,
                                             cObj.doPerturbations,cObj.pertFile,cObj.minimizePerts,false,false);

        rs.initMutRes2Str(cObj.strandMut);

        //System.out.println("Beginning setAllowables");
        //Set the allowable AAs for each AS residue
        boolean addWT = (new Boolean((String)cObj.params.getValue("ADDWT", "true"))).booleanValue();
        if(!addWT)
            checkWT(cObj.strandDefault, cObj.strandPresent, cObj.params);
        int molStrand = 0;
        for (int strNum=0; strNum<cObj.strandPresent.length; strNum++) {
            if(cObj.strandPresent[strNum]) {
                for (int k=0; k<cObj.strandMut[molStrand].length; k++) {
                    setAllowablesHelper(rs, cObj.params, addWT, strNum, molStrand, k, cObj.strandMut, cObj.strandDefault);
                }
                molStrand++;
            }
        }


        if(cObj.doPerturbations) {
            rs.setupRCs(cObj.addWTRot);
            Perturbation.idealizeSC = cObj.idealizeSC;
        }


        boolean shellRun = false;
        boolean intraRun = false;
        boolean templateOnly = false;

        if (cObj.flagMutType.compareTo("TEMPL")==0) {

            // **** Normal Active Site residue runs ****
            // Computes active site residue to active site residue pair energies
            shellRun = true;/*ligPresent = true*/;
            intraRun = false;
            templateOnly = true;
        } else if (cObj.flagMutType.compareTo("AS-AS")==0) {

            // **** Normal Active Site residue runs ****
            // Computes active site residue to active site residue pair energies
            shellRun = false;/*ligPresent = true*/;
            intraRun = false;
            templateOnly = false;
        } else if (cObj.flagMutType.compareTo("SHL-AS")==0) {

            // Then shell runs for the active site residues
            // Computes the active site residue rotamers to shell energies
            shellRun = true;/*ligPresent = true*/;
            intraRun = false;
            templateOnly = false;
        } else if (cObj.flagMutType.compareTo("INTRA")==0) {

            // Compute all intra-residue energies
            shellRun = false;
            intraRun = true;
            templateOnly = false;
        } else if (cObj.flagMutType.compareTo("LIG-AS")==0) {

            // **** Ligand present runs ****
            // This section computes the inter-residue energies between
            //  active site residues and the ligand
            shellRun = false;
            intraRun = false;
            templateOnly = false;
        } else { //(cObj.flagMutType.compareTo("LIG-SHL")==0)

            // Computes ligand rotamer to shell energies
            shellRun = true;
            intraRun = false;
            templateOnly = false;
        }

        // The goal is that the total energy of a system can be bounded by the sum of
        //  all pairwise active site residue entries plus the entry for each active site
        //  residue's shell run plus each active site residue's self intra-residue energy.
        //  If a ligand is present then one should add the ligand to shell energy, the
        //  ligand to each active site residue pairwise energy, and the ligand self intra-
        //  residue energy.

        //initialize the pairwise energy matrices (partial initialization - only for the residues involved in this computation, e.g., AS-AS)
        PairwiseEnergyMatrix minEmatrix = new PairwiseEnergyMatrix(cObj.mutableSpots,cObj.resMut,cObj.strandMut,
                rs,shellRun,intraRun,false);
        PairwiseEnergyMatrix maxEmatrix = null;
        if(cObj.doMinimization)
            maxEmatrix = minEmatrix.copy();
        else  //PGC
            maxEmatrix = minEmatrix;


        //Compute the corresponding matrix entries
        rs.simplePairwiseMutationAllRotamerSearch(cObj.strandMut,cObj.mutableSpots,cObj.doMinimization,shellRun,intraRun,
                cObj.resMut,minEmatrix,maxEmatrix,cObj.minimizeBB,cObj.doBackrubs,
                templateOnly,cObj.backrubFile, cObj.templateAlwaysOn);


        long stopTime = System.currentTimeMillis();
        cObj.elapsedTime = Math.round((stopTime - startTime) / 1000.0f);


        //Store the information in less space to allow the master node to buffer several cObj at once
        cObj.compEE = minEmatrix.generateCompEE(maxEmatrix);

        //KER: reset evalAtoms
        m.setAllEval();

        return cObj;

    }

    //Generates and saves to file the mutation list for the pairwise energy matrix computation
    private OneMutation[] getMutArrayPairEcomp(int numberMutable, boolean minimizeBB) {

        int numMutations = 2; //pairwise energy computation

        //KER: sometimes there will only be one thing flexible (for example if strand only has one residue flexible)
        if(numberMutable < numMutations)
            numMutations = numberMutable;

        // Generate all combinations
        int numComb = factorial(numberMutable).divide(factorial(numberMutable-numMutations).multiply(factorial(numMutations))).intValue();
        int residueMutatable[][] = new int[numComb][numberMutable];
        generateCombinations(residueMutatable,numberMutable,numMutations);
        // At this point each row of residueMutatble is a 0/1 array which specifies a mutation
        //  pair combination, 1 indicates that that residue can mutate in the specified combination

        System.out.println("Number of possible mutation combinations: "+numComb);

        // Create the mutation list with estimated energies
        OneMutation mutArray[] = new OneMutation[numComb];

        //Set the AS-AS mutations
        int curMutNum = 0;
        //KER: if we only have one mutable position there can't be pairwise interactions
        if(numMutations >= 2) {
            for(int i=0; i<numComb; i++) {

                mutArray[i] = new OneMutation();
                mutArray[i].flagMutType = "AS-AS";
                mutArray[i].resMut = new int[numberMutable];
                for(int j=0; j<numberMutable; j++) {
                    mutArray[i].resMut[j] = residueMutatable[i][j];
                }

                // Perform simple mutation search for this set of mutatable residues
                curMutNum++;
            }
        }

        //Add the runs for template only, AS-shell, AS-ligand, intra-residue energies, and ligand-shell
        int numOtherMut;
        /*if (ligPresent)
        	numOtherMut = 2+t+2*numInAS;
        else
        	numOtherMut = 1+t+numInAS;*/
        numOtherMut = 2+numberMutable;
        OneMutation otherMutArray[] = new OneMutation[numOtherMut];

        for (int i=0; i<numOtherMut; i++) {
            otherMutArray[i] = new OneMutation();
            otherMutArray[i].resMut = new int[numberMutable];
        }

        //Set the AS-shell mutations
        for (int i=0; i<numberMutable; i++) {
            otherMutArray[i].flagMutType = "SHL-AS";
            for (int j=0; j<numberMutable; j++) {
                if (i==j)
                    otherMutArray[i].resMut[j] = 1;
                else
                    otherMutArray[i].resMut[j] = 0;
            }
        }

        //Set the intra-residue energies run
        otherMutArray[numberMutable].flagMutType = "INTRA";
        for (int j=0; j<numberMutable; j++)
            otherMutArray[numberMutable].resMut[j] = 1;

        //Set the template energy run
        otherMutArray[numberMutable+1].flagMutType = "TEMPL";
        for (int j=0; j<numberMutable; j++)
            otherMutArray[numberMutable+1].resMut[j] = 0;


        /*if (ligPresent){//if the ligand is present, set the corresponding runs

        	//Set the AS-ligand mutations
        	for (int i=1+t+numInAS; i<=2*numInAS+t; i++){
        		otherMutArray[i].flagMutType = "LIG-AS";
        		for (int j=0; j<numInAS; j++){
        			if ((i-1-t-numInAS)==j)
        				otherMutArray[i].resMut[j] = 1;
        			else
        				otherMutArray[i].resMut[j] = 0;
        		}
        	}

        	//Set the ligand-shell run
        	otherMutArray[1+t+2*numInAS].flagMutType = "LIG-SHL";
        	for (int j=0; j<numInAS; j++)
        		otherMutArray[1+t+2*numInAS].resMut[j] = 0;
        }*/

        // We now have all the mutations in mutArray, collapse the mutArray
        //  to the actual number of mutations we have.
        OneMutation newArray[] = new OneMutation[curMutNum+numOtherMut];
        System.arraycopy(otherMutArray,0,newArray,0,numOtherMut);//add the other mutations first
        System.arraycopy(mutArray,0,newArray,numOtherMut,curMutNum);//then add the AS-AS mutations

        mutArray = newArray;
        System.out.println("Length of mutArray: "+mutArray.length);

        return mutArray;
    }

    // Finds the index of the mutation in resumeResults with the same
    //  mutation sequence as resMut. If none are found, -1 is returned.
    public int sampFindMutationIndex(OneMutation resumeResults[], String flMutType, int mutResidues[]) {

        for(int q=0; q<resumeResults.length; q++)
            if ((resumeResults[q].flagMutType.compareTo(flMutType)==0) && (sameSeq(resumeResults[q].resMut,mutResidues)))
                return(q);

        return(-1);
    }

    //Determines if the residues that can mutate are the same for two mutation sequences
    private boolean sameSeq (int computedResMut[], int allResMut[]) {

        boolean found = true;
        for (int i=0; i<computedResMut.length; i++) {
            if (computedResMut[i]!=allResMut[i])
                found = false;
        }
        return found;
    }
///////////////////////////////////////////////////////////////////////////
//	End of MIN and MAX Pairwise Energy Precomputation
///////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////
//	 Compute minimized-GMEC section
////////////////////////////////////////////////////////////////
    /**
     * Computes the (energy-minimized) structure for the rotameric conformations specified by the input file.
     * This function is used to generate structures for a set of output conformations from a DEE/A* search.
     */
    public void handleMinDEEApplyRot(String s) {

        // Takes the following parameters
        // 1: System parameter filename (string)
        // 2: Mutation search parameter filename (string)

        ParamSet sParams = new ParamSet();
        sParams.addParamsFromFile(getToken(s,2)); //read system parameters
        sParams.addParamsFromFile(getToken(s,3)); //read mutation search parameters

        // Pull search parameters
        String runName = (String)sParams.getValue("RUNNAME");
        String confResFile = (String)sParams.getValue("CONFRESFILE", "c_"+runName);
        int numResults = (new Integer((String)sParams.getValue("NUMRESULTS"))).intValue();
        boolean doMinimize = (new Boolean((String)sParams.getValue("DOMINIMIZE", "false"))).booleanValue();
        boolean minimizeBB = (new Boolean((String)sParams.getValue("MINIMIZEBB","false"))).booleanValue();
        boolean doBackrubs = (new Boolean((String)sParams.getValue("DOBACKRUBS", "false"))).booleanValue();
        String backrubFile ="";
        if(doBackrubs) {
            backrubFile = (String)sParams.getValue("BACKRUBFILE");
        }
        boolean outputPDB = (new Boolean((String)sParams.getValue("OUTPUTPDBS", "true"))).booleanValue();

        boolean doPerturbations = (new Boolean((String)sParams.getValue("DOPERTURBATIONS","false"))).booleanValue();//Triggers DEEPer
        String pertFile = (String)sParams.getValue("PERTURBATIONFILE","defaultPerturbationFileName.pert");//Input file giving perturbation information
        boolean minimizePerts = (new Boolean((String)sParams.getValue("MINIMIZEPERTURBATIONS","false"))).booleanValue();//Allow continuous minimization with respect to perturbation parameters
        Perturbation.idealizeSC = (new Boolean((String)sParams.getValue("IDEALIZESIDECHAINS","true"))).booleanValue();


        boolean addWTRotsSomehow = (new Boolean((String)sParams.getValue("ADDWTROTS","false"))).booleanValue();
        boolean addOrigRots = false, addWTRot = false;
        if( addWTRotsSomehow ) {
            if(doPerturbations)//DEEPer allows adding WT rotamers only to the positions they came from, so do that if possible
                addWTRot = true;
            else//Otherwise just add them to the rotamer library
                addOrigRots = true;
        }



        MolParameters mp = loadMolecule(sParams, COMPLEX);
        Molecule m = mp.m;
        int numberMutable = mp.numberMutable;
        int strandsPresent = mp.strandsPresent;
        String[][] strandLimits = mp.strandLimits;
        boolean[] strandPresent = mp.strandPresent;
        int[][] strandMut = mp.strandMut;
        String[][] strandDefault = mp.strandDefault;
        int numOfStrands = strandMut.length;
        int mutRes2Strand[] = mp.mutRes2Strand;
        int mutRes2StrandMutIndex[] = mp.mutRes2StrandMutIndex;

        if(addOrigRots)
            RotamerLibrary.addOrigRots(strandMut, EnvironmentVars.aaRotLib, m);

        int numResidues = numberMutable;

        if(outputPDB) {
            //KER: make pdbs directory if it doesn't exist
            File pdbDir = new File("pdbs");
            if(!pdbDir.exists())
                pdbDir.mkdir();
        }

        //Read the results file into the AA and rot matrices
        String AAtypes[] = new String[numResults*numResidues];
        int rotNums[] = new int[numResults*numResidues];
        readRotResFile(confResFile,AAtypes,rotNums,numResults,numResidues);
        String curSeq[] = new String[numberMutable];
        int curSeqInd[] = new int[numberMutable];
        int curRot[] = new int[numberMutable];

        int numSaved = 0;
        for (int curResult=0; curResult<numResults; curResult++) {


            if( doPerturbations && curResult>0 ) { //Reload the molecule
                mp = loadMolecule(sParams, COMPLEX);
                m = mp.m;
            }


            System.out.print("Starting minimization of result "+(curResult+1)+"..");

            Amber96ext a96ff = new Amber96ext(m, distDepDielect, dielectConst, doSolvationE, solvScale, softvdwMultiplier);

            StrandRotamers[] strandRot = new StrandRotamers[numOfStrands];
            for(int i=0; i<numOfStrands; i++) {
                if(doPerturbations)
                    strandRot[i] = new StrandRCs(grl[i],m.strand[i]);
                else
                    strandRot[i] = new StrandRotamers(grl[i],m.strand[i]);
            }

            //the starting index for the current result within AAtypes[] and rotNums[]
            int startInd = numResidues*curResult;
            a96ff.calculateTypesWithTemplates();

            int ctr = 0;
            for(int str=0; str<numOfStrands; str++) {
                for(int j=0; j<strandMut[str].length; j++) {
                    curSeq[ctr] = AAtypes[startInd+ctr];
                    curSeqInd[ctr] = strandRot[str].rl.getAARotamerIndex(curSeq[ctr]);
                    strandRot[str].setAllowable(strandMut[str][j],curSeq[ctr]);

                    if( doPerturbations && addWTRot )//Store this before mutating anything
                        ((StrandRCs)strandRot[str]).storeWTRotamer(m, strandMut[str][j]);

                    if(m.strand[str].isProtein)
                        strandRot[str].changeResidueType(m,strandMut[str][j],curSeq[ctr],true,true);
                    m.strand[str].residue[strandMut[str][j]].flexible = true;
                    ctr++;
                }
            }


            if(doPerturbations) { //This has to be done after setting the allowables
                PertFileHandler.readPertFile(pertFile, m, strandRot);
                for(int str=0; str<numOfStrands; str++) {
                    ((StrandRCs)strandRot[str]).addUnperturbedRCs(addWTRot);
                    ((StrandRCs)strandRot[str]).countRCs();
                }
            }

            a96ff.calculateTypesWithTemplates();
            a96ff.initializeCalculation();
            a96ff.setNBEval(hElect,hVDW);
            //TODO: Fix the energy function so we can set each strand
            //a96ff.setLigandNum((new Integer((String)sParams.getValue("PDBLIGNUM"))).intValue());

            int curAA[] = new int[m.numberOfResidues];
            for(int str=0; str<numOfStrands; str++) {
                for(int j=0; j<m.strand[str].numberOfResidues; j++) {
                    int molResNum = m.strand[str].residue[j].moleculeResidueNumber;
                    curAA[molResNum] = strandRot[str].getIndexOfNthAllowable(j,0);
                }
            }

            SimpleMinimizer simpMin = null;
            BBMinimizer bbMin = null;
            BackrubMinimizer brMin = null;
            if ( doMinimize && (!minimizeBB) ) {
                if(doPerturbations)
                    simpMin = new PMinimizer(minimizePerts);
                else
                    simpMin = new SimpleMinimizer();
                /*if(ligPresent)
                	simpMin.initialize(m,0,1,a96ff,sysLR,ligLR2,curAA,ligLR2.getIndexOfNthAllowable(0,0),doDihedE,rl,grl);
                else*/
                simpMin.initialize(m,numOfStrands,a96ff,strandRot,curAA,doDihedE);
            } else if (minimizeBB) {
                if (!doBackrubs) {
                    bbMin = new BBMinimizer();
                    /*if (ligPresent)
                    	bbMin.initialize(m, a96ff, residueMap, sysStrNum, ligStrNum);
                    else*/
                    bbMin.initialize(m, a96ff, strandMut, numOfStrands);
                } else {
                    brMin = new BackrubMinimizer();
                    /*if (ligPresent)
                    	brMin.initialize(m, a96ff, residueMap, sysStrNum, ligStrNum, backrubFile, hSteric, stericThresh);
                    else*/
                    brMin.initialize(m, a96ff, strandMut, backrubFile, hSteric, stericThresh,numOfStrands, true);
                }
            }

            //In DEEPer we will need to reload the molecule for each result instead of backing up
            //(the backup restoration would mess up the perturbations)
            if(!doPerturbations)
                m.backupAtomCoord();

            //Apply the corresponding rotamers
            ctr = 0;
            for(int str=0; str<numOfStrands; str++)
                for (int j=0; j<strandMut[str].length; j++) {
                    int molResNum = m.strand[str].residue[strandMut[str][j]].moleculeResidueNumber;

                    if(doPerturbations) {
                        curRot[ctr] = rotNums[startInd+ctr];
                        boolean outcome = ((StrandRCs)strandRot[str]).applyRC(m, strandMut[str][j], curRot[ctr]);
                        if( ! outcome )
                            System.err.println("Invalid RC:" + curRot[ctr] + " at position " + strandMut[str][j] + " of strand " + str );
                    } else if(strandRot[str].rl.getNumRotForAAtype(curAA[molResNum]) != 0) { //not GLY or ALA
                        curRot[ctr] = rotNums[startInd+ctr];
                        strandRot[str].applyRotamer(m, strandMut[str][j], curRot[ctr]);
                    }
                    ctr++;
                }

            double unMinE[] = a96ff.calculateTotalEnergy(m.actualCoordinates,-1); //the energy before minimization
            if (outputPDB) { //save molecule
                String fName = runName;
                String filename = String.format("pdbs/%1$s_%2$03d_unmin.pdb",fName,numSaved);
                saveMolecule(m,filename,(float)unMinE[0]);
            }
            if ( doMinimize && (!minimizeBB) )
                simpMin.minimize(35);
            else if (minimizeBB) {
                if (!doBackrubs)
                    bbMin.minimizeFull(false);
                else {
                    brMin.minimizeFull();
                }
            }

            double minE[] = a96ff.calculateTotalEnergy(m.actualCoordinates,-1); //the energy after minimization
            if ((doMinimize)&&(!minimizeBB)&&(doDihedE)) //add dihedral energies
                minE[0] += simpMin.computeDihedEnergy();

            if (outputPDB) { //save molecule
                String fName = runName;
                String filename = String.format("pdbs/%1$s_%2$03d_min.pdb",fName,numSaved);
                saveMolecule(m,filename,(float)minE[0]);
                numSaved++;
            }

            if(!doPerturbations) {
                m.restoreAtomCoord();
                m.updateCoordinates();
            }

            if (minE[0]>unMinE[0])
                minE = unMinE;

            System.out.println("done");
        }
        System.out.println("done");
    }

    private void readRotResFile (String resFile, String AAtypes[], int rotNums[], int numResults, int numResidues) {

        BufferedReader bufread = null;
        try {
            File file = new File(resFile);
            FileReader fr = new FileReader(file);
            bufread = new BufferedReader(fr);
        } catch (FileNotFoundException e) {
            System.out.println("ERROR: Results File Not Found");
        }

        boolean done = false;
        String str = null;
        int resultNum = 0;

        while (!done) {
            try {
                str = bufread.readLine();
            } catch ( Exception e ) {
                System.out.println("ERROR: An error occurred while reading input");
                System.exit(0);
            }

            if (str == null) // stop if we've reached EOF
                done = true;
            else if(resultNum >= numResults)
                done = true;
            else {
                for (int i=0; i<numResidues; i++) {
                    AAtypes[resultNum*numResidues+i] = getToken(str,2+i);
                    rotNums[resultNum*numResidues+i] = new Integer((String)getToken(str,2+numResidues+i)).intValue();
                }
                resultNum++;
            }
        }

        if (numResults!=resultNum) {
            System.out.println("Error: Not all results available for reading");
            System.exit(0);
        }

        // We're done reading them in
        try {
            bufread.close();
        } catch(Exception e) {
        }
    }



    //Identify the current rotamers in the flexible residues
    //This function is sort of the inverse of handleMinDEEApplyRot
    public void identifyRotamers(String s) {

        ParamSet sParams = new ParamSet();
        sParams.addParamsFromFile(getToken(s,2)); //read system parameters
        sParams.addParamsFromFile(getToken(s,3)); //read mutation search parameters


        //Setup the molecule system
        MolParameters mp = loadMolecule(sParams, COMPLEX);
        int numMutable = mp.numberMutable;

        //Data to output
        String AAtypes[] = new String[numMutable];
        int curRot[] = new int[numMutable];
        float maxDev[] = new float[numMutable];//Maximum deviation of a dihedral from the rigid-rotamer value


        //Find the rotamer at each position with the minimax deviation from the real structure
        for(int j=0; j<numMutable; j++) {

            maxDev[j] = Float.POSITIVE_INFINITY;

            //Get information on the residue
            int str = mp.mutRes2Strand[j];
            int strResNum = mp.strandMut[str][mp.mutRes2StrandMutIndex[j]];
            Residue res = mp.m.strand[str].residue[strResNum];
            AAtypes[j] = res.name;

            //Choose the correct rotamer library for this residue
            RotamerLibrary curRL = grl[str];


            int AATypeInd = curRL.getAARotamerIndex( AAtypes[j] );
            int numDihedrals = curRL.getNumDihedrals( AATypeInd );
            int symmCount = ResSymmetry.getSymmetry( AAtypes[j] );//Get the number of symmetry states of the residue

            double curDih[][] = new double[numDihedrals][symmCount];//Current values for the dihedrals
            //in each of the symmetry states

            for(int a=0; a<numDihedrals; a++) {

                for(int symm=0; symm<symmCount; symm++) {

                    int molAtNum[] = new int[4];//Molecule atoms numbers for atoms involved in the dihedral

                    for( int b=0; b<4; b++ ) {
                        String atName = curRL.getDihedralAtomNames(AATypeInd, a, b);
                        atName = ResSymmetry.getPermutedAtomName(atName, AAtypes[j], symm);
                        molAtNum[b] = res.getAtomNameToMolnum(atName);
                    }

                    curDih[a][symm] = mp.m.getTorsion( molAtNum[0], molAtNum[1], molAtNum[2], molAtNum[3] );
                }
            }



            int numRot = curRL.getNumRotamers( AAtypes[j] );

            for(int rot=0; rot<numRot; rot++) { //Check all the possible rotamers

                for(int symm=0; symm<symmCount; symm++) { //Check each symmetry state against them

                    float rotMaxDev = Float.NEGATIVE_INFINITY;
                    //Maximum dihedral deviation between the structure (in the current symmetry)
                    //and the current rotamer

                    for(int k=0; k<numDihedrals; k++) {

                        float idealDih = curRL.getRotamerValues(AATypeInd, rot, k);

                        float dev = (float)Math.abs( idealDih - curDih[k][symm] );//the angles' absolute values can't exceed 180 so this is at most 360
                        dev = Math.min(dev, 360-dev);//dev is the absolute difference between the angles, which are mod 360

                        if( dev > rotMaxDev )
                            rotMaxDev = dev;
                    }

                    if( rotMaxDev < maxDev[j] ) { //This rotamer is the best so far
                        maxDev[j] = rotMaxDev;
                        curRot[j] = rot;
                    }
                }
            }

            if(numRot == 0) {
                curRot[j] = 0;
                maxDev[j] = 0;
            }

        }


        //Now write out the answer
        //First all the amino-acid types
        System.out.print( "AA types: " );
        for(int j=0; j<numMutable; j++)
            System.out.print( AAtypes[j] + " " );
        System.out.println();

        System.out.print( "Rotamers: " );
        for(int j=0; j<numMutable; j++)
            System.out.print( curRot[j] + " " );
        System.out.println();

        System.out.print( "Max dihedral deviations: " );
        for(int j=0; j<numMutable; j++)
            System.out.print( maxDev[j] + " " );
        System.out.println();
    }

////////////////////////////////////////////////////////////////
//	 End of Compute minimized-GMEC section
////////////////////////////////////////////////////////////////




///////////////////////////////////////////////////////////////////////////
//	DEE section
///////////////////////////////////////////////////////////////////////////
    /**
     * Performs a DEE (Traditional, MinDEE, BD, BRDEE, or DEEPer) pruning with A* search, with or without DACS;
     * the only parameter 's' (String) includes the command-line arguments specifying the filenames of the two input configuration files.
     * If distributed DACS is performed, the computation is distributed to the available processors for evaluation.
    */
    public void handleDoDEE(String s) {

        // Takes the following parameters
        // 1: System parameter filename (string)
        // 2: DEE config filename (string)

        long startTimeAll = System.currentTimeMillis();

        System.out.println("Performing DEE");

        ParamSet sParams = new ParamSet();
        sParams.addParamsFromFile(getToken(s,2)); //read system parameters
        sParams.addParamsFromFile(getToken(s,3)); //read mutation search parameters

        // Pull search parameters

        String runName = ((String)sParams.getValue("RUNNAME"));
        int numMaxMut = (new Integer((String)sParams.getValue("NUMMAXMUT", "1000"))).intValue();
        int algOption = (new Integer((String)sParams.getValue("ALGOPTION", "3"))).intValue();
        boolean doDACS = (new Boolean((String)sParams.getValue("DODACS", "false"))).booleanValue();
        boolean useFlags = (new Boolean((String)sParams.getValue("SPLITFLAGS", "false"))).booleanValue();
        boolean distrDACS = (new Boolean((String)sParams.getValue("DISTRDACS", "false"))).booleanValue();
        //boolean distrDEE = (new Boolean((String)sParams.getValue("DISTRDEE"))).booleanValue();
        boolean distrDEE = false; //the distributed DEE section is outdated and must be carefully checked before called
        boolean doMinimize = (new Boolean((String)sParams.getValue("DOMINIMIZE", "false"))).booleanValue();
        boolean minimizeBB = (new Boolean((String)sParams.getValue("MINIMIZEBB", "false"))).booleanValue();
        boolean doBackrubs = (new Boolean((String)sParams.getValue("DOBACKRUBS", "false"))).booleanValue();
        String backrubFile = "";
        if(doBackrubs) {
            backrubFile = ((String)sParams.getValue("BACKRUBFILE"));
        }
        boolean approxMinGMEC = (new Boolean((String)sParams.getValue("APPROXMINGMEC", "false"))).booleanValue();
        boolean preprocPairs = (new Boolean((String)sParams.getValue("PREPROCPAIRS", "true"))).booleanValue();
        boolean scaleInt = (new Boolean((String)sParams.getValue("SCALEINT", "false"))).booleanValue();
        String runNameEMatrixMin = (String)(sParams.getValue("MINENERGYMATRIXNAME",runName+"minM" ));
        String runNameEMatrixMax = (String)(sParams.getValue("MAXENERGYMATRIXNAME",runName+"maxM" ));
        float initEw = (new Float((String)sParams.getValue("INITEW","0"))).floatValue();
        float lambda = (new Float((String)sParams.getValue("LAMBDA", "0"))).floatValue();
        float pruningE = (new Float((String)sParams.getValue("PRUNINGE", "100.0"))).floatValue();
        double stericE = (new Double((String)sParams.getValue("STERICE","30.0"))).doubleValue();
        float pairSt = (new Float((String)sParams.getValue("PAIRST", "100.0"))).floatValue();
        float maxIntScale =0;
        if(scaleInt) {
            maxIntScale = (new Float((String)sParams.getValue("MAXINTSCALE","0"))).floatValue();
        }
        double minRatioDiff=0;
        int initDepth = 0;
        int subDepth = 0;
        int diffFact =0;
        if(doDACS) {
            minRatioDiff = (new Double((String)sParams.getValue("MINRATIODIFF"))).doubleValue();
            initDepth = (new Integer((String)sParams.getValue("INITDEPTH"))).intValue();
            subDepth = (new Integer((String)sParams.getValue("SUBDEPTH"))).intValue();
            diffFact = (new Integer((String)sParams.getValue("DIFFFACT"))).intValue();
        }
        boolean genInteractionGraph = (new Boolean((String)sParams.getValue("GENINTERACTIONGRAPH","false"))).booleanValue();
        float distCutoff=0;
        float eInteractionCutoff=0;
        if(genInteractionGraph) {
            distCutoff = (new Float((String)sParams.getValue("DISTCUTOFF"))).floatValue();
            eInteractionCutoff = (new Float((String)sParams.getValue("EINTERACTIONCUTOFF"))).floatValue();
        }
        String outputConfInfo = (String)(sParams.getValue("OUTPUTCONFINFO","c_"+runName));
        String outputPruneInfo = (String)(sParams.getValue("OUTPUTPRUNEINFO","p_"+runName));

        int curStrForMatrix = (new Integer((String)sParams.getValue("ONLYSINGLESTRAND","-1"))).intValue();
        boolean typeDep = (new Boolean((String)sParams.getValue("TYPEDEP","false"))).booleanValue();

        boolean useEref = (new Boolean((String)sParams.getValue("USEEREF","true"))).booleanValue();
        boolean resumeSearch = (new Boolean((String)sParams.getValue("RESUMESEARCH","false"))).booleanValue();
        String resumeFilename ="";
        if(resumeSearch) {
            resumeFilename = ((String)sParams.getValue("RESUMEFILENAME"));
        }

        // 2010: Use energy window MinDEE method.  If this is set to true,
        //   MinDEE will use traditional DEE with an energy window (initEw)
        //   for pruning.  Max terms will be ignored and only the min terms for pruning and
        boolean useMinDEEPruningEw = (new Boolean((String)sParams.getValue("imindee", "false"))).booleanValue();
        float Ival = 0.0f;
        float interval = 0;
        if(useMinDEEPruningEw)
            Ival = (new Float((String)sParams.getValue("IVAL"))).floatValue();
        float difference = Ival;


        boolean useTriples = (new Boolean((String)sParams.getValue("USETRIPLES","false"))).booleanValue();
        boolean magicBulletTriples = (new Boolean((String)sParams.getValue("MAGICBULLETTRIPLES","true"))).booleanValue();//Use only "magic bullet" competitor triples
        int magicBulletNumTriples = (new Integer((String)sParams.getValue("MAGICBULLETNUMTRIPLES","5"))).intValue();//Number of magic bullet triples to use
        boolean useFlagsAStar = (new Boolean((String)sParams.getValue("USEFLAGSASTAR","false"))).booleanValue();


        // DEEPer parameters
        boolean doPerturbations = (new Boolean((String)sParams.getValue("DOPERTURBATIONS","false"))).booleanValue();//Triggers DEEPer
        boolean pertScreen = (new Boolean((String)sParams.getValue("PERTURBATIONSCREEN","false"))).booleanValue();//Triggers perturbation screen: pruning-only run with rigid perturbations
        String pertFile = (String)sParams.getValue("PERTURBATIONFILE","defaultPerturbationFileName.pert");//Input file giving perturbation information
        boolean minimizePerts = (new Boolean((String)sParams.getValue("MINIMIZEPERTURBATIONS","false"))).booleanValue();//Allow continuous minimization with respect to perturbation parameters
        String screenOutFile = ((String)sParams.getValue("SCREENOUTFILE","screenOutFileDefaultName.pert"));//Name of file for outputting results of screen (same format as PERTURBATIONFILE)
        boolean selectPerturbations = (new Boolean((String)sParams.getValue("SELECTPERTURBATIONS","false"))).booleanValue();//Should perturbations be automatically selected?
        Perturbation.idealizeSC = (new Boolean((String)sParams.getValue("IDEALIZESIDECHAINS","true"))).booleanValue();


        boolean addWTRotsSomehow = (new Boolean((String)sParams.getValue("ADDWTROTS","false"))).booleanValue();
        boolean addOrigRots = false, addWTRot = false;
        if( addWTRotsSomehow ) {
            if(doPerturbations)//DEEPer allows adding WT rotamers only to the positions they came from, so do that if possible
                addWTRot = true;
            else//Otherwise just add them to the rotamer library
                addOrigRots = true;
        }

        if( !useFlags && (useFlagsAStar || useTriples || algOption >=4) ) { //These all rely heavily on the split flags
            System.err.println("ERROR: Options requiring split flags (flags in A*, triples pruning, or algOption>=4) are set but split flags are turned off");
            System.exit(1);
        }
        if( doMinimize && (!useMinDEEPruningEw) && ( useTriples || algOption >=4 ) ) {
            System.err.println("ERROR: Options requiring iMinDEE (triples pruning and/or algOption >=4) are set but iMinDEE is turned off");
            System.exit(1);
        }

        if ((!mpiRun)&&((distrDACS)||distrDEE)) {
            System.out.println("ERROR: Distributed computation requires MPI");
            System.exit(1);
        }

        if (!doMinimize) //no minimization
            minimizeBB = false;
        if (!minimizeBB) //not backbone minimization
            doBackrubs = false;

        if (genInteractionGraph) //DACS is not performed when generating the interaction graph
            doDACS = false;

        //Setup the molecule system
        MolParameters mp = loadMolecule(sParams,curStrForMatrix);

        if(addOrigRots)
            RotamerLibrary.addOrigRots(mp.strandMut, EnvironmentVars.aaRotLib, mp.m);


        boolean reload=false;


        if(selectPerturbations) { //Need to run the automatic perturbation selection
            //This only needs to be done once though: after that the perturbations can be read from pertFile
            selectPerturbations(mp, doPerturbations, pertFile, minimizePerts, addWTRot, sParams);
            reload = true;
        }


        // 2010: If useMinDEEPruningEw is set to false, this cycle never repeats itself.
        //  If it is set to true, it can repeat at most once: if none of the rotamer vectors
        //  between the conformation of lowest energy (i.e. lowestBound)
        //  and lowestBound+InitEw can minimize to a lower energy than lowestBound+InitEw,
        //   then let minimumEnergy be the minimum nergy found among the enumerated conformations,
        //   we set a new Ew equal to minimumEnergy - lowestBount and repeat this cycle.  We
        //   only have to do it at most twice.
        do {

            if( reload )
                mp = loadMolecule(sParams,curStrForMatrix);
            else if(doPerturbations)
                reload = true;


            RotamerSearch rs = new RotamerSearch(mp.m,mp.numberMutable, mp.strandsPresent, hElect, hVDW, hSteric, true,
                                                 true, 0.0f, stericThresh, softStericThresh, distDepDielect, dielectConst, doDihedE, doSolvationE, solvScale, softvdwMultiplier, grl,
                                                 doPerturbations, pertFile, minimizePerts, useTriples, useFlagsAStar);

            rs.initMutRes2Str(mp.strandMut);

            System.out.print("Loading precomputed energy matrix...");
            loadPairwiseEnergyMatrices(sParams,rs,runNameEMatrixMin,doMinimize,runNameEMatrixMax,curStrForMatrix);
            System.out.println("done");


            /////////////////////////////////////////////////////////////
            // DEE section

            long startTime = System.currentTimeMillis();

            //Set the allowable AAs for each AS residue
            boolean addWT = (new Boolean((String)sParams.getValue("ADDWT", "true"))).booleanValue();
            if(!addWT)
                checkWT(mp.strandDefault, mp.strandPresent, sParams);
            int molStrand = 0;
            for (int strNum=0; strNum<mp.numOfStrands; strNum++) {
                if(mp.strandPresent[strNum]) {
                    for (int k=0; k<mp.strandMut[molStrand].length; k++) {
                        setAllowablesHelper(rs, sParams, addWT, strNum, molStrand, k, mp.strandMut, mp.strandDefault);
                    }
                    molStrand++;
                }
            }

            if(doPerturbations)
                rs.setupRCs(addWTRot);

            float eRef[][] = null;
            if (useEref) { //add the reference energies to the min (and max) intra-energies
                String eRefName = sParams.getValue("EREFMATRIXNAME", "Eref");
                eRef = RotamerSearch.loadErefMatrix(eRefName+".dat");
                //eRef = getResEntropyEmatricesEref(useEref,rs.getMinMatrix(),rs.strandRot,mp.strandMut,null,mp.numberMutable,mp.mutRes2Strand,mp.mutRes2StrandMutIndex);
                rs.addEref(eRef, doMinimize, mp.strandMut);
            }
            if(EnvironmentVars.useEntropy)
                rs.addEntropyTerm(doMinimize, mp.strandMut);


            PrunedRotamers<Boolean> prunedRotAtRes = new PrunedRotamers<Boolean>(mp.numberMutable,mp.strandMut,rs,false);

            final String rotFile = ("rot_out"+System.currentTimeMillis()); //output the pruned rotamers (for distributed DEE)
            final String sfFile = ("sf_matrix"+System.currentTimeMillis()); //output the split flags (for distributed DEE)

            int[] numRotForRes = compNumRotForRes(mp.numberMutable, rs, mp.strandMut, mp.mutRes2Strand,mp.mutRes2StrandMutIndex);
            for(int q:numRotForRes) {
                System.out.print(q+" ");
            };
            System.out.println("");

            //first prune all rotamers that are incompatible with the template (intra E + re-to-template E >= stericE)
            System.out.println();
            System.out.println("Pruning all rotamers incompatible with the template..");
            prunedRotAtRes = rs.DoPruneStericTemplate(mp.numberMutable, mp.strandMut,
                             prunedRotAtRes, stericE);
            System.out.println();

            int[] unprunedRotForRes = rotamersRemaining(numRotForRes, prunedRotAtRes);
            for(int q:unprunedRotForRes) {
                System.out.print(q+" ");
            }
            System.out.println();

            //preprocess pairs of rotamers (mark pairwise energies greater than the cutoff as steric clashes)
            //This is only useful if these pairs will not be pruned by rs.prunedRidiculousPairs (i.e. if !useFlags)
            if (preprocPairs && (!useFlags) ) {
                System.out.println("Preprocessing pairs of rotamers, cutoff of "+pairSt);
                rs.preprocessPairs(pairSt, mp.numberMutable,mp.strandMut);
                System.out.println();
            }


            //Setup and do the DEE pruning
            if ((useFlags)||(algOption>=3))
                rs.setSplitFlags();//initialize the split flags

            if(useFlags) {
                float cutoff = (float)stericE;
                if ( preprocPairs && (pairSt < stericE) )
                    cutoff = pairSt;
                rs.pruneRidiculousPairs(mp.numberMutable, mp.strandMut,
                                        prunedRotAtRes, cutoff);
            }

            int numPrunedRot = 0;
            int numPrunedPairs = 0;
            int numPrunedRotThisRun = 0;
            int numPrunedPairsThisRun = 0;
            boolean done = false;
            int numRuns = 1;

            boolean localUseMinDEEPruningEw = useMinDEEPruningEw;
            if(doDACS) { //If we are doing dacs we don't want to prune stuff too early so turn off
                localUseMinDEEPruningEw = false;   //iMinDEE until the last depth is reached

                //Without iMinDEE we can't do triples or indirect pruning though
                if(useTriples) {
                    System.out.println("Warning: Can't prune triples with DACS.  Turning off triples pruning.");
                    useTriples = false;
                    rs.useTriples = false;
                }
                if(algOption >= 4) {
                    System.out.println("Warning: Can't do indirect pruning (algOption 4) with DACS.  Reverting to algOption 3");
                    algOption = 3;
                }
            }


            while (!done) { //repeat the pruning cycle until no more rotamers are pruned

                numPrunedRotThisRun = 0;
                numPrunedPairsThisRun = 0;

                System.out.println("Starting DEE cycle run: "+numRuns);

                if (doMinimize && !localUseMinDEEPruningEw) //precompute the interval terms in the MinDEE criterion
                    rs.doCompMinDEEIntervals(mp.numberMutable, mp.strandMut, prunedRotAtRes, scaleInt, maxIntScale);

                //Depending on the chosen algorithm option, apply the corresponding pruning criteria;
                System.out.println("Starting pruning with DEE (simple Goldstein)");
                prunedRotAtRes = rs.DoDEEGoldstein(mp.numberMutable, mp.strandMut,
                                                   initEw, prunedRotAtRes, doMinimize, useFlags, minimizeBB,typeDep,localUseMinDEEPruningEw,Ival);
                System.out.println();

                if ((algOption>=3)) { //simple Goldstein pairs
                    System.out.println("Starting pruning with DEE (mb pairs)");
                    rs.DoDEEPairs(mp.numberMutable, mp.strandMut,
                                  initEw, prunedRotAtRes, null, doMinimize, useFlags, true,
                                  false, minimizeBB, scaleInt, maxIntScale,typeDep,localUseMinDEEPruningEw,Ival);
                    System.out.println();
                }

                if ((useFlags)||(algOption>=3)) {
                    System.out.println("Starting pruning with Bounding Flags");
                    rs.DoBoundFlags(mp.numberMutable, mp.strandMut,
                                    pruningE, prunedRotAtRes, initEw, useFlags);
                    System.out.println();
                }

                System.out.println("Starting pruning with DEE (1-sp split-DEE)");
                prunedRotAtRes = rs.DoDEEConfSplitting(mp.numberMutable, mp.strandMut,
                                                       initEw, prunedRotAtRes, null, doMinimize, useFlags, 1, false, minimizeBB,
                                                       typeDep,localUseMinDEEPruningEw,Ival);
                System.out.println();

                System.out.println("Starting pruning with Bounds");
                prunedRotAtRes = rs.DoMinBounds(mp.numberMutable, mp.strandMut,
                                                pruningE, prunedRotAtRes, initEw, useFlags, false);
                System.out.println();

                //check how many rotamers/pairs are pruned this run
                int numTotalPrunedRot = countPrunedRot(prunedRotAtRes);
                int numTotalPrunedPairs = 0;
                if ((useFlags)||(algOption>=3)) //pairs pruning is performed
                    numTotalPrunedPairs = countPrunedPairs(rs.getSplitFlags());

                if ((numTotalPrunedRot!=numPrunedRot)||(numTotalPrunedPairs!=numPrunedPairs)) { //new rotamers/pairs pruned this run

                    numPrunedRotThisRun += numTotalPrunedRot - numPrunedRot;
                    numPrunedPairsThisRun += numTotalPrunedPairs - numPrunedPairs;

                    numPrunedRot = numTotalPrunedRot; //update the number of pruned rotamers
                    numPrunedPairs = numTotalPrunedPairs;
                    numRuns++;
                } else { //no more rotamers pruned, so perform the computationally-expensive 2-sp split-DEE and pairs

                    if (!doDACS) { //DACS will not be performed

                        if ((algOption>=3)) { //simple Goldstein pairs
                            System.out.println("Starting pruning with DEE (full pairs)");

                            if (distrDEE) { //distributed full Goldstein pairs DEE
                                doDistrDEEMaster(mp.numberMutable, mp.strandMut,
                                                 sParams, mp.strandDefault,
                                                 runNameEMatrixMin, runNameEMatrixMax, initEw, prunedRotAtRes, doMinimize,
                                                 rs, sfFile, rotFile, useFlags, -1, optPairs, minimizeBB, scaleInt,
                                                 maxIntScale, useEref, eRef, addWTRot);

                                rs.setSplitFlags((boolean [][][][][][])readObject(sfFile)); //get the DE pairs from the distributed run
                            } else { //perform on a single processor
                                rs.DoDEEPairs(mp.numberMutable, mp.strandMut,
                                              initEw, prunedRotAtRes, null, doMinimize, useFlags,
                                              false, false, minimizeBB, scaleInt, maxIntScale, typeDep,localUseMinDEEPruningEw,Ival);
                            }
                            System.out.println();
                        }

                        if ((algOption>=2)) { //2-sp conf splitting
                            System.out.println("Starting pruning with DEE (2-sp split-DEE)");

                            if (distrDEE) { //distributed 2-sp split-DEE
                                doDistrDEEMaster(mp.numberMutable, mp.strandMut,
                                                 sParams, mp.strandDefault,
                                                 runNameEMatrixMin, runNameEMatrixMax, initEw, prunedRotAtRes, doMinimize,
                                                 rs, sfFile, rotFile, useFlags, 2, optSplit, minimizeBB, scaleInt,
                                                 maxIntScale, useEref, eRef, addWTRot);

                                prunedRotAtRes = (PrunedRotamers<Boolean>)readObject(rotFile);
                            } else { //perform on a single processor
                                prunedRotAtRes = rs.DoDEEConfSplitting(mp.numberMutable, mp.strandMut,
                                                                       initEw, prunedRotAtRes, null, doMinimize, useFlags, 2,
                                                                       false, minimizeBB, typeDep,localUseMinDEEPruningEw,Ival);
                            }
                            System.out.println();
                        }

                        if( useTriples ) {
                            //Prune Goldstein triples

                            System.out.println("Starting triples pruning");
                            rs.DoDEETriples(mp.numberMutable, mp.strandMut, initEw, prunedRotAtRes, null,
                                            doMinimize, magicBulletTriples, magicBulletNumTriples, false,
                                            minimizeBB, typeDep, localUseMinDEEPruningEw, Ival);

                            System.out.println();
                        }

                        if(algOption >= 4) {
                            //Indirect pruning goes at the end because it benefits strongly from prior pruned pairs and triples
                            System.out.println("Starting indirect pruning");

                            prunedRotAtRes = rs.DoDEEIndirect(mp.numberMutable, mp.strandMut, initEw,
                                                              prunedRotAtRes, null, doMinimize, false, false,
                                                              minimizeBB, typeDep, localUseMinDEEPruningEw, Ival);

                            System.out.println();
                            //This prunes full rather than magic-bullet pairs
                        }


                        //check if 2-sp split-DEE and pairs pruned new rotamers
                        numTotalPrunedRot = countPrunedRot(prunedRotAtRes);
                        numTotalPrunedPairs = 0;
                        if ((useFlags)||(algOption>=3)) //pairs pruning is performed
                            numTotalPrunedPairs = countPrunedPairs(rs.getSplitFlags());

                        if ((numTotalPrunedRot==numPrunedRot)&&(numTotalPrunedPairs==numPrunedPairs)) //no more rotamers/pairs pruned
                            done = true;
                        else { //new rotamers pruned this run
                            numPrunedRotThisRun += numTotalPrunedRot - numPrunedRot;
                            numPrunedPairsThisRun += numTotalPrunedPairs - numPrunedPairs;

                            numPrunedRot = numTotalPrunedRot; //update the number of pruned rotamers
                            numPrunedPairs = numTotalPrunedPairs;
                            numRuns++;
                        }
                    } else //DACS will be performed
                        done = true;
                }
                System.out.println("Num pruned rot this run: "+numPrunedRotThisRun);
                System.out.println("Num pruned pairs this run: "+numPrunedPairsThisRun);
                System.out.println();
            }

            long pruneTime = System.currentTimeMillis();

            if(pertScreen) { //No A*, so we're done
                PertFileHandler.writePertFile(screenOutFile, rs.m, prunedRotAtRes, rs.strandRot, mp.strandMut, true);
                System.out.println("Screening time: "+((pruneTime-startTime)/(60.0*1000.0)));
                System.out.println("Screening done");
                return;
            }

            if (!doDACS) { //DACS will not be performed

                if (genInteractionGraph) //generate interaction graph
                    genInteractionGraph(mp.numberMutable, rs, prunedRotAtRes, runName, mp.strandMut, eInteractionCutoff, distCutoff, mp.m, preprocPairs, pairSt,mp.mutRes2Strand,mp.mutRes2StrandMutIndex);

                else { //perform A* search to enumerate conformations
                    double bestScore = Math.pow(10,38); //for DEE/A*, the initial best score is the highest possible

                    // 2010: A* now returns a new value for Ew.  Note that right now useMinDEEPruningEw
                    //   is incompatible with the traditional usage of initEw.  This can be easily
                    //   fixed by adding a different Ew for this method.  If A* returns an initEw of
                    //    -1 that means that an error occured somewhere.  If it returns 0 it means
                    //    that the GMEC or minGMEC was found.  If it returns 'Ew'>0 then it means that
                    //    useMinDEEPruningEw = true and that the energy window must be "enlarged"
                    //    to 'Ew'
                    //  useTopKHeuristic: Only the top Kvalue conformations are enumerated.  If
                    //     not enough, a new initEw is returned. Note that we may have to do this
                    //     several times.

                    interval = rs.doAStarGMEC(outputConfInfo,true,doMinimize,
                                              mp.numberMutable,mp.strandMut,mp.strandDefault,numMaxMut,initEw,
                                              bestScore,null,approxMinGMEC,lambda,minimizeBB,useEref,eRef,doBackrubs,backrubFile,
                                              localUseMinDEEPruningEw, Ival);

                    difference = interval - Ival;
                    Ival = interval;

                }
            } else { //DACS
                numRotForRes = compNumRotForRes(mp.numberMutable, rs, mp.strandMut, mp.mutRes2Strand,mp.mutRes2StrandMutIndex);
                BigInteger numInitUnprunedConfs = compNumUnprunedConfs(mp.numberMutable, prunedRotAtRes, numRotForRes);

                int msp[] = new int[mp.numberMutable]; //split positions (for DACS)
                for (int i=0; i<msp.length; i++)
                    msp[i] = -1;

                if (distrDACS) { //distributed DACS (only for level 0)

                    if (initDepth<=0) {
                        System.out.println("ERROR: distributed DACS called with 'initDepth="+initDepth+"' partitioning positions; use a positive value");
                        System.exit(1);
                    }
                    if (subDepth<0)
                        subDepth = 0;

                    distrDEE = false; //do not perform both distributed DACS and distributed DEE

                    //choose the major splitting positions
                    for (int i=0; i<initDepth; i++)
                        msp[i] = chooseSplitPos(mp.numberMutable,prunedRotAtRes,numRotForRes, msp, i, minRatioDiff);

                    int maxNumPartitions = 1;
                    for (int i=0; i<initDepth; i++)
                        maxNumPartitions *= numRotForRes[msp[i]];

                    OneMutation resumeResults[] = null;
                    if (resumeSearch) { //read resume results
                        System.out.println("Reading resume results..");
                        resumeResults = new OneMutation[maxNumPartitions];
                        for(int q=0; q<resumeResults.length; q++)
                            resumeResults[q] = new OneMutation();
                        resumeResults = readResumeFile(resumeResults,resumeFilename,mp.numberMutable,true,false,initDepth);
                        System.out.println("Read "+resumeResults.length+" completed partitions.");
                        System.out.println();
                    }

                    doDistrDACSMaster(runName, mp.numberMutable, rs, mp.strandMut, mp.strandDefault,
                                      rotFile, prunedRotAtRes, algOption, sfFile, useFlags, initEw, pruningE, initDepth, msp,
                                      numInitUnprunedConfs, diffFact, outputPruneInfo, outputConfInfo, minRatioDiff, doMinimize,
                                      runNameEMatrixMin, runNameEMatrixMax, sParams, approxMinGMEC,
                                      lambda, numRotForRes, resumeResults, resumeFilename, minimizeBB, numMaxMut, scaleInt, maxIntScale,
                                      useEref, eRef, doBackrubs, backrubFile, subDepth,mp.mutRes2Strand,mp.mutRes2StrandMutIndex,mp.strandPresent,
                                      mp.strandLimits,mp.strandsPresent,addWTRot);
                } else { //single-processor DACS

                    initDepth = 0; //only used for distributed DACS
                    if (subDepth<=0) {
                        System.out.println("ERROR: single-processor DACS called with 'subDepth="+subDepth+"' partitioning positions; use a positive value");
                        System.exit(1);
                    }

                    PrintStream logPS = setupOutputFile(outputPruneInfo);

                    doDACS(mp.numberMutable, rs, mp.strandMut, mp.strandDefault,
                           rotFile, prunedRotAtRes, algOption, sfFile, useFlags, initEw, pruningE, initDepth, 0, logPS, msp,
                           numInitUnprunedConfs, diffFact, outputConfInfo, minRatioDiff, doMinimize, runNameEMatrixMin,
                           runNameEMatrixMax, distrDEE, sParams, approxMinGMEC, lambda, null,
                           null, minimizeBB, numMaxMut, scaleInt, maxIntScale, useEref, eRef, doBackrubs, backrubFile, subDepth,
                           mp.mutRes2Strand,mp.mutRes2StrandMutIndex,typeDep, addWTRot);
                }
            }

            long stopTime = System.currentTimeMillis();

            System.out.println("Pruning time: "+((pruneTime-startTime)/(60.0*1000.0)));
            if (genInteractionGraph)
                System.out.println("Graph generation time: "+((stopTime-pruneTime)/(60.0*1000.0)));
            else
                System.out.println("Enumeration/DACS time: "+((stopTime-pruneTime)/(60.0*1000.0)));
            System.out.println("DEE execution time: "+((stopTime-startTime)/(60.0*1000.0)));
            System.out.println("DEE done");
            System.out.println("Total execution time: "+((stopTime-startTimeAll)/(60.0*1000.0)));
            //end of DEE section
            /////////////////////////////////////////////////////////////
        } while(difference > 0.001 && useMinDEEPruningEw && doMinimize); // 2010: if I1-I0 >0 we must repeat the cycle with the new energy
        // window.  This can only happen if useMinDEEPruningEw is true
        //   and not topK


    }

    //KER: If addWT isn't true and there is only one AA type for a mutable
    //position we treat that position as WT
    private void checkWT(String[][] strandDefault,boolean[] strandPresent, ParamSet sParams) {
        int ctr = 0;
        for(int str=0; str<strandPresent.length; str++) {
            if(strandPresent[str]) {
                for(int i=0; i<strandDefault[ctr].length; i++) {
                    String tempResAllow = (String)sParams.getValue("RESALLOWED"+str+"_"+i);
                    if(numTokens(tempResAllow)==1)
                        strandDefault[ctr][i] = getToken(tempResAllow, 1);
                }
                ctr++;
            }
        }
    }


    private class MolParameters {
        String[][] strandLimits;
        int numOfStrands;
        boolean[] strandPresent;
        int strandsPresent;
        int[][] strandMut;
        String[][] strandDefault;
        Molecule m;
        int numberMutable;
        int mutRes2Strand[];
        int mutRes2StrandMutIndex[];
    }

    public MolParameters loadMolecule(ParamSet sParams, int curStrForMatrix) {

        MolParameters mp = new MolParameters();

        loadStrandParams(sParams, mp, curStrForMatrix);

        //Setup the molecule system
        mp.m = new Molecule();
        mp.m = setupMolSystem(mp.m,sParams,mp.strandPresent,mp.strandLimits);

        loadMutationParams(sParams, mp);

        int totalLength=0;
        for(int i=0; i<mp.strandMut.length; i++)
            for(int j=0; j<mp.strandMut[i].length; j++)
                totalLength++;

        mp.mutRes2Strand = new int[totalLength];
        mp.mutRes2StrandMutIndex = new int[totalLength];
        int ctr=0;
        for(int i=0; i<mp.strandMut.length; i++)
            for(int j=0; j<mp.strandMut[i].length; j++) {
                mp.mutRes2Strand[ctr] = i;
                mp.mutRes2StrandMutIndex[ctr] = j;
                ctr++;
            }

        mp.numberMutable = getNumberMutable(mp.strandMut);

        return mp;

    }

    private int getNumberMutable(int[][] strandMut) {
        int numberMutable = 0;
        for (int i=0; i<strandMut.length; i++)
            numberMutable += strandMut[i].length;

        return numberMutable;
    }

    private void loadStrandParams(ParamSet sParams, MolParameters mp, int curStrForMatrix) {
        mp.numOfStrands = (new Integer((String)sParams.getValue("NUMOFSTRANDS"))).intValue();
        mp.strandLimits = new String[mp.numOfStrands][2];
        mp.strandPresent = new boolean[mp.numOfStrands];
        mp.strandsPresent = 0;
        for (int i=0; i<mp.numOfStrands; i++) {
            String strandLimit = (String)sParams.getValue("STRAND"+i);
            String limit1 = getToken(strandLimit,1);
            String limit2 = getToken(strandLimit,2);
            mp.strandLimits[i][0] = limit1;
            mp.strandLimits[i][1] = limit2;
            if(curStrForMatrix == COMPLEX)
                mp.strandPresent[i] = true;
            else if(curStrForMatrix == i)
                mp.strandPresent[i] = true;
            else
                mp.strandPresent[i] = false;
            //strandPresent[i] = (new Boolean((String)sParams.getValue("STRANDPRESENT" + i))).booleanValue();
            if(mp.strandPresent[i]) {
                mp.strandsPresent++;
            }
        }
    }



    //KER: mp must already have the following terms set:
    //mp.numOfStrands
    //mp.strandLimits
    //mp.strandPresent
    //mp.strandsPresent
    private void loadMutationParams(ParamSet sParams, MolParameters mp) {
        //Setup the molecule system
        //mp.m = new Molecule();
        //mp.m = setupMolSystem(mp.m,sParams,mp.strandPresent,mp.strandLimits);


        boolean addWT = (new Boolean((String)sParams.getValue("ADDWT","true"))).booleanValue();
        /**********Get the regions of each strand that are mutable****/
        mp.strandMut = new int[mp.strandsPresent][];  //taking the place of resMap and ligMap
        mp.strandDefault = new String[mp.strandsPresent][];
        String strandMutNums = (String)sParams.getValue("STRANDMUTNUMS");
        int ctr = 0;
        try {
            for(int i=0; i<mp.strandPresent.length; i++) {
                if(mp.strandPresent[i]) {
                    int numberOfMutables = (new Integer(getToken(strandMutNums,i+1))).intValue();
                    mp.strandMut[ctr] = new int[numberOfMutables];
                    String strandMutResNum = (String)sParams.getValue("STRANDMUT"+i);
                    for(int j=0; j<numberOfMutables; j++) {
                        String strandMutRes = getToken(strandMutResNum,j+1);
                        mp.strandMut[ctr][j] = mp.m.residue[mp.m.mapPDBresNumToMolResNum(strandMutRes)].strandResidueNumber;
                    }
                    mp.strandDefault[ctr] = new String[mp.strandMut[ctr].length];
                    for(int j=0; j<mp.strandMut[ctr].length; j++) {
                        mp.strandDefault[ctr][j] = mp.m.strand[ctr].residue[mp.strandMut[ctr][j]].name;
                    }
                    ctr++;
                }
            }
            if(!addWT)
                checkWT(mp.strandDefault, mp.strandPresent, sParams);
        } catch(Exception E) {
            System.out.println("PROBLEM Loading the strand mut numbers (Check System.cfg) ");
            E.printStackTrace();
            return;
        }
    }

    //Implements threads for DACS: allows the current best score among all partitions to be distributed to every partition;
    //This thread performs the DACS computation, while the main thread monitors for updates of the best energy from the other partitions;
    //The communication is performed via the common RotamerSearch object (synchronized access to the bestEMin variable)


    /**
     * Performs the DACS partition-specific computation. The parameters 'rs' and 'rs.sysLR' must be valid.
     */
    private void doDACS(int numMutable, RotamerSearch rs, int strandMut[][],
                        String strandDefault[][], String rotFile,
                        PrunedRotamers<Boolean> prunedRotAtRes, int algOption, String sfFile, boolean useFlags, float initEw, float pruningE,
                        int initDepth, int curDepth, PrintStream logPS, int majorSplitPos[], BigInteger numInitUnprunedConfs,
                        int diffFact, String outputConfInfo, double minRatioDiff, boolean doMinimize, String minPEM,
                        String maxPEM, boolean distrDEE, ParamSet sParams,
                        boolean approxMinGMEC, float lambda, Index3 partIndex[], CommucObj cObj,
                        boolean minimizeBB, int numMaxMut, boolean scaleInt, float maxIntScale, boolean useEref, float eRef[][],
                        boolean doBackrubs, String backrubFile, int subDepth,int mutRes2Strand[],int mutRes2StrandMutIndex[],boolean typeDep, boolean addWTRot) {

        if (curDepth>=(initDepth+subDepth))
            return;

        System.out.println("Starting pruning with DEE (DACS).");

        //prunedRotAtRes[] should not be modified here, in order to be able to distinguish
        //	newly pruned rotamers and rotamers pruned by Bounds or DEE

        //the num rotamers for each AS residue and the ligand (if present);
        //	sysLR in rs must be valid (with all the possible AA's for each residue position)
        int numRotForRes[] = compNumRotForRes(numMutable, rs, strandMut,mutRes2Strand,mutRes2StrandMutIndex);

        if (curDepth>=initDepth) {//sub-partition, so majorSplitPos[curDepth] is unknown; compute it
            majorSplitPos[curDepth] = chooseSplitPos(numMutable,prunedRotAtRes,numRotForRes,
                                      majorSplitPos, curDepth, minRatioDiff);//the splitting position
        }

        int numPartitions = numRotForRes[majorSplitPos[curDepth]]; //the number of partitions
        int numPrunedPartitions = 0; //the number of partitions with no unpruned confs

        System.out.println("Current depth: "+curDepth);
        System.out.println("Major splitting residue number: "+majorSplitPos[curDepth]);
        System.out.println();

        //map rotamer index to rot num for the splitting residue
        Index3[] indexMap = getIndexMap(numPartitions,rs,majorSplitPos[curDepth],strandMut,mutRes2Strand,mutRes2StrandMutIndex);

        //count the total num confs and the num unpruned confs
        BigInteger numConfsTotalForPartition[] = new BigInteger[numPartitions];
        BigInteger numConfsUnprunedForPartition[] = new BigInteger[numPartitions];
        //BigInteger numInitConfsUnprunedForPartition[] = new BigInteger[numPartitions];
        BigInteger numTotalConfs = new BigInteger("0");
        BigInteger numUnprunedConfs = new BigInteger("0");
        BigInteger numEvaluatedConfs = new BigInteger("0");

        //update the best energy found
        pruningE = (float)Math.min(pruningE, rs.getBestE());
        double bestScore = pruningE;

        boolean savedSpFlags[][][][][][] = null;
        if (useFlags||(algOption>=3))
            savedSpFlags = rs.getSplitFlags(savedSpFlags);

        //determine the prunings for each of the sub-solutions (the partitions)
        PrunedRotamers<Boolean> prunedForPartition[] = new PrunedRotamers[numPartitions];//[prunedRotAtRes.length];
        for(int i=0; i<numPartitions; i++)
            prunedForPartition[i] = new PrunedRotamers<Boolean>(prunedRotAtRes,false);
        for (int i=0; i<prunedForPartition.length; i++) {

            if ((curDepth>=initDepth)||(indexMap[i]==partIndex[curDepth])) { //sub-partitions or current partition is the partition distributed for computation

                //copy the prunings from before conf splitting (from Bounds and simple traditional DEE)
                //System.arraycopy(prunedRotAtRes,0,prunedForPartition[i],0,prunedRotAtRes.length);
                Iterator<RotInfo<Boolean>> iter1 = prunedRotAtRes.iterator();
                while(iter1.hasNext()) {
                    RotInfo<Boolean> ri = iter1.next();
                    prunedForPartition[i].set(ri, new Boolean(ri.state));
                }

                Index3 curPartIndex = indexMap[i]; //the rotamer index of the partitioning rotamer

                //artificially set all rotamers at the splitting position, other than the rotamer
                //	for the current partition, to pruned, so that there will be only one rotamer at
                //	that residue position when the conf splitting criterion is applied;
                //	when done, subtract the artifcially pruned rotamers from the total number of pruned
                boolean indToUnprune[] = new boolean[numPartitions];
                for (int j=0; j<indToUnprune.length; j++)
                    indToUnprune[j] = false;

                //check the partition only if the current partitioning rotamer is not already pruned
                if (!prunedRotAtRes.get(curPartIndex)) {

                    for (int j=0; j<numPartitions; j++) {
                        if (j!=i) {//not the rotamer for the current partition
                            Index3 curInd = indexMap[j]; //the index of the current rotamer
                            if (!prunedForPartition[i].get(curInd)) { //not pruned by the other DEE methods
                                prunedForPartition[i].set(curInd, true);
                                indToUnprune[j] = true;
                            }
                        }
                    }

                    if ((useFlags)||(algOption>=3))
                        rs.setSplitFlags(savedSpFlags); //for each partition, reset the flags to the globally valid ones

                    int numPrunedRot = countPrunedRot(prunedForPartition[i]);
                    int numPrunedPairs = 0;
                    if ((useFlags)||(algOption>=3))
                        numPrunedPairs = countPrunedPairs(rs.getSplitFlags());
                    int numPrunedRotThisRun = 0;
                    int numPrunedPairsThisRun = 0;
                    boolean done = false;
                    int numRuns = 1;

                    while (!done) { //repeat the pruning cycle until no more rotamers are pruned

                        numPrunedRotThisRun = 0;
                        numPrunedPairsThisRun = 0;

                        System.out.println("Starting DEE cycle run: "+numRuns);

                        if (doMinimize) //precompute the interval terms in the MinDEE criterion
                            rs.doCompMinDEEIntervals(numMutable, strandMut,
                                                     prunedForPartition[i], scaleInt, maxIntScale);

                        System.out.println("Starting pruning with DEE (simple Goldstein)");
                        prunedForPartition[i] = rs.DoDEEGoldstein(numMutable, strandMut,
                                                initEw, prunedForPartition[i], doMinimize, useFlags,
                                                minimizeBB, typeDep, false, 0.0f);
                        System.out.println();

                        if ((algOption>=3)) { //simple Goldstein pairs
                            System.out.println("Starting pruning with DEE (mb pairs)");
                            rs.DoDEEPairs(numMutable, strandMut,
                                          initEw, prunedForPartition[i], null, doMinimize,
                                          useFlags, true, false, minimizeBB, scaleInt, maxIntScale,typeDep, false, 0.0f);
                            System.out.println();
                        }

                        if ((useFlags)||(algOption>=3)) {
                            System.out.println("Starting pruning with Bounding Flags");
                            rs.DoBoundFlags(numMutable, strandMut,
                                            pruningE, prunedForPartition[i], initEw, useFlags);
                            System.out.println();
                        }

                        System.out.println("Starting pruning with DEE (1-sp split-DEE)");
                        prunedForPartition[i] = rs.DoDEEConfSplitting(numMutable, strandMut,
                                                initEw, prunedForPartition[i], null, doMinimize,
                                                useFlags, 1, false, minimizeBB,typeDep, false, 0.0f);
                        System.out.println();

                        System.out.println("Starting pruning with Bounds");
                        prunedForPartition[i]= rs.DoMinBounds(numMutable, strandMut,
                                                              pruningE, prunedForPartition[i], initEw, useFlags, false);
                        System.out.println();

                        //check how many rotamers/pairs are pruned this run
                        int numTotalPrunedRot = countPrunedRot(prunedForPartition[i]);
                        int numTotalPrunedPairs = 0;
                        if ((useFlags)||(algOption>=3))
                            numTotalPrunedPairs = countPrunedPairs(rs.getSplitFlags());

                        if ((numTotalPrunedRot!=numPrunedRot)||(numTotalPrunedPairs!=numPrunedPairs)) { //new rotamers/pairs pruned this run

                            numPrunedRotThisRun += numTotalPrunedRot - numPrunedRot;
                            numPrunedPairsThisRun += numTotalPrunedPairs - numPrunedPairs;

                            numPrunedRot = numTotalPrunedRot; //update the number of pruned rotamers
                            numPrunedPairs = numTotalPrunedPairs;
                            numRuns++;
                        } else { //no more rotamers pruned, so perform the computationally-expensive 2-sp split-DEE and pairs

                            if ((algOption>=3)&&(curDepth>=(initDepth+subDepth-1))) { //simple Goldstein pairs
                                System.out.println("Starting pruning with DEE (full pairs)");

                                if (distrDEE) { //distributed full Goldstein pairs DEE
                                    doDistrDEEMaster(numMutable, strandMut,
                                                     sParams, strandDefault,
                                                     minPEM, maxPEM, initEw, prunedForPartition[i], doMinimize, rs, sfFile, rotFile,
                                                     useFlags, -1, optPairs, minimizeBB, scaleInt, maxIntScale, useEref, eRef, addWTRot);

                                    rs.setSplitFlags((boolean [][][][][][])readObject(sfFile)); //get the DE pairs from the distributed run
                                } else { //perform on a single processor
                                    rs.DoDEEPairs(numMutable, strandMut,
                                                  initEw, prunedForPartition[i], null, doMinimize,
                                                  useFlags, false, false, minimizeBB, scaleInt, maxIntScale,typeDep, false, 0.0f);
                                }
                                System.out.println();
                            }

                            if ((algOption>=2)) { //2-sp conf splitting
                                System.out.println("Starting pruning with DEE (2-sp split-DEE)");

                                if (distrDEE) { //distributed 2-sp split-DEE
                                    doDistrDEEMaster(numMutable, strandMut,
                                                     sParams, strandDefault,
                                                     minPEM, maxPEM, initEw, prunedForPartition[i], doMinimize, rs, sfFile, rotFile,
                                                     useFlags, 2, optSplit, minimizeBB, scaleInt, maxIntScale, useEref, eRef, addWTRot);

                                    prunedForPartition[i] = ((PrunedRotamers<Boolean>)readObject(rotFile));
                                } else {
                                    prunedForPartition[i] = rs.DoDEEConfSplitting(numMutable, strandMut,
                                                            initEw, prunedForPartition[i], null, doMinimize,
                                                            useFlags, 2, false, minimizeBB,typeDep, false, 0.0f);
                                }
                                System.out.println();
                            }

                            //check if 2-sp split-DEE and pairs pruned new rotamers/pairs
                            numTotalPrunedRot = countPrunedRot(prunedForPartition[i]);
                            if ((useFlags)||(algOption>=3))
                                numTotalPrunedPairs = countPrunedPairs(rs.getSplitFlags());

                            if ((numTotalPrunedRot==numPrunedRot)&&(numTotalPrunedPairs==numPrunedPairs)) //no more rotamers/pairs pruned
                                done = true;
                            else { //new rotamers pruned this run
                                numPrunedRotThisRun += numTotalPrunedRot - numPrunedRot;
                                numPrunedPairsThisRun += numTotalPrunedPairs - numPrunedPairs;

                                numPrunedRot = numTotalPrunedRot; //update the number of pruned rotamers
                                numPrunedPairs = numTotalPrunedPairs;
                                numRuns++;
                            }
                        }
                        System.out.println("Num pruned rot this run: "+numPrunedRotThisRun);
                        System.out.println("Num pruned pairs this run: "+numPrunedPairsThisRun);
                        System.out.println();
                    }
                }

                //count the number of pruned rotamers for each residue (except for the splitting residue,
                //	which always has only 1 available rotamer for the current partition)
                int numPrunedRotForRes[] = new int[numMutable]; //after pruning for this partition
                //int numInitPrunedRotForRes[] = new int[numInAS]; //initial prunings, before pruning for this partition
                for (int j=0; j<numMutable; j++) {
                    numPrunedRotForRes[j] = 0;
                    if (j!=majorSplitPos[curDepth]) {
                        Iterator<RotInfo<Boolean>> iter = prunedForPartition[i].iterator(j);
                        RotInfo<Boolean> ri = iter.next();
                        while(ri.curPos == j && iter.hasNext()) {
                            //for (int k=0; k<totalNumRotamers; k++){ //pruned rot are true and must be in the current set of allowed AA
                            if (prunedForPartition[i].get(ri))
                                numPrunedRotForRes[j]++;
                            //if (prunedRotAtRes[j*totalNumRotamers + k])
                            //numInitPrunedRotForRes[j]++;
                            ri = iter.next();
                        }
                    } else { // j==majorSplitPos
                        numPrunedRotForRes[j] = 0;
                        //numInitPrunedRotForRes[j] = 0;
                    }
                }
                /*int numPrunedLigRot = 0;
                //int numInitPrunedLigRot = 0;
                for (int k=0; k<numLigRotamers; k++){
                	if (prunedForPartition[i][numInAS*totalNumRotamers + k])
                		numPrunedLigRot++;
                	//if (prunedRotAtRes[numInAS*totalNumRotamers + k])
                	//	numInitPrunedLigRot++;
                }*/

                //count the total num confs and the num unpruned confs
                numConfsTotalForPartition[i] = new BigInteger("1");
                numConfsUnprunedForPartition[i] = new BigInteger("1");
                //numInitConfsUnprunedForPartition[i] = new BigInteger("1");
                if (prunedRotAtRes.get(curPartIndex)) { //current partitioning rotamer already pruned, so no unpruned confs for this partition
                    numConfsUnprunedForPartition[i] = new BigInteger("0");
                    numPrunedPartitions++;
                }

                for (int j=0; j<numMutable; j++) {
                    if (!(isSplitRes(j,majorSplitPos,curDepth))) { //the split residues contribute only 1 rotamer
                        numConfsTotalForPartition[i] = numConfsTotalForPartition[i].multiply(BigInteger.valueOf(numRotForRes[j]));
                        numConfsUnprunedForPartition[i] = numConfsUnprunedForPartition[i].multiply(BigInteger.valueOf(numRotForRes[j]-numPrunedRotForRes[j]));
                        //numInitConfsUnprunedForPartition[i] = numInitConfsUnprunedForPartition[i].multiply(BigInteger.valueOf(numRotForRes[j]-numInitPrunedRotForRes[j]));
                    }
                }
                /*if(ligPresent){
                	numConfsTotalForPartition[i] = numConfsTotalForPartition[i].multiply(BigInteger.valueOf(numLigRotamers));
                	numConfsUnprunedForPartition[i] = numConfsUnprunedForPartition[i].multiply(BigInteger.valueOf(numLigRotamers-numPrunedLigRot));
                	//numInitConfsUnprunedForPartition[i] = numInitConfsUnprunedForPartition[i].multiply(BigInteger.valueOf(numLigRotamers-numInitPrunedLigRot));
                }*/

                numTotalConfs = numTotalConfs.add(numConfsTotalForPartition[i]);
                numUnprunedConfs = numUnprunedConfs.add(numConfsUnprunedForPartition[i]);

                BigInteger pruneDiffFact = BigInteger.valueOf(10).pow(diffFact);

                System.out.println("Num unpruned confs: "+numConfsUnprunedForPartition[i]+" diffFact: "+pruneDiffFact);
                System.out.println();
                System.out.println();

                //output pruning info to file
                logPS.print("curDepth: "+curDepth+" curPartition: "+i+" majorSplitPos: "+majorSplitPos[curDepth]+" ");
                logPS.print(numConfsTotalForPartition[i]+" "+numInitUnprunedConfs+" "+numConfsUnprunedForPartition[i]);
                logPS.println();
                logPS.println();
                logPS.flush();

                //if ((curDepth+1<maxDepth)&&(numConfsUnprunedForPartition[i].compareTo(numInitUnprunedConfs.divide(pruneDiffFact))==1)){ //not enough pruned, so partition at new depth
                if ((curDepth+1<(initDepth+subDepth))&&(numConfsUnprunedForPartition[i].compareTo(pruneDiffFact)==1)) { //not enough pruned, so partition at new depth
                    doDACS(numMutable, rs, strandMut,strandDefault,
                           rotFile, prunedForPartition[i], algOption, sfFile, useFlags, initEw, pruningE, initDepth, curDepth+1,
                           logPS, majorSplitPos, numConfsUnprunedForPartition[i], diffFact, outputConfInfo, minRatioDiff,
                           doMinimize, minPEM, maxPEM, distrDEE, sParams, approxMinGMEC,
                           lambda, partIndex, cObj, minimizeBB, numMaxMut, scaleInt, maxIntScale, useEref, eRef, doBackrubs,
                           backrubFile, subDepth,mutRes2Strand,mutRes2StrandMutIndex,typeDep,addWTRot);
                } else if (!prunedRotAtRes.get(curPartIndex)) { //if enough pruned or maxDepth partitioning reached, do the rotamer search

                    bestScore = Math.min(bestScore,rs.getBestE());//best E for the partitions so far

                    //Do the rotamer search
                    rs.doAStarGMEC(outputConfInfo,true,doMinimize,numMutable,strandMut,
                                   strandDefault,numMaxMut,initEw,bestScore,null,approxMinGMEC,lambda,minimizeBB,useEref,eRef,doBackrubs,
                                   backrubFile, false, 0.0f);

                    numEvaluatedConfs = numEvaluatedConfs.add(rs.numConfsEvaluated); //add the evaluated confs for this partition
                    pruningE = (float)Math.min(pruningE,rs.getBestE());//update cutoff energy for MinBounds
                }

                //unprune the artificially pruned indices
                for (int j=0; j<indToUnprune.length; j++) {
                    if (indToUnprune[j]) {
                        prunedForPartition[i].set(indexMap[j], false);
                    }
                }
            }
        }

        if (cObj==null) { //not distributed DACS

            System.out.println("numTotalConfs: "+numTotalConfs+"; numUnprunedConfs: "+numUnprunedConfs+"; numEvaluatedConfs: "+numEvaluatedConfs);
            for (int i=0; i<numPartitions; i++)System.out.print(numConfsTotalForPartition[i]+" ");
            System.out.println();
            for (int i=0; i<numPartitions; i++)System.out.print(numConfsUnprunedForPartition[i]+" ");
            System.out.println();
            System.out.println("Major splitting residue number: "+majorSplitPos[curDepth]);
            System.out.println("Number of partitions: "+numPartitions);
            System.out.println("Number of non-zero partitions: "+(numPartitions-numPrunedPartitions));
            System.out.println("Additional pruned rotamers: ");

            //count the number of partitions for which a rotamer is pruned (counting only the rotamers
            //	not pruned by Bounds or simple traditional DEE)
            PrunedRotamers<Integer> countNumPartitionsPrunedRot = new PrunedRotamers<Integer>(prunedRotAtRes,0);
            Iterator<RotInfo<Boolean>> iter = prunedRotAtRes.iterator();
            while(iter.hasNext()) {
                RotInfo<Boolean> ri = iter.next();
                //for (int i=0; i<prunedRotAtRes.length; i++){//for each rotamer
                //countNumPartitionsPrunedRot[i] = 0;
                if (!ri.state) { //only if not pruned by the other two methods
                    for (int j=0; j<numPartitions; j++) { //check for each partition
                        if (prunedForPartition[j].get(ri))
                            countNumPartitionsPrunedRot.set(ri, countNumPartitionsPrunedRot.get(ri)+1);
                    }
                }

                //output information
                if (countNumPartitionsPrunedRot.get(ri)>0)
                    System.out.println("index: "+ri.printCoord()+"; num partitions in which pruned: "+countNumPartitionsPrunedRot.get(ri));
            }
        } else { //distributed DACS
            cObj.bestScore = cObj.bestScore.min(BigDecimal.valueOf(rs.getBestE()));
        }
    }

    /**
     * Implements threads for DACS: allows the current best score among all partitions to be distributed to every partition;
     * This thread performs the DACS computation, while the main thread monitors for updates of the best energy from the other partitions;
     * The communication is performed via the common RotamerSearch object (synchronized access to the bestEMin variable)
    */
    private class DACSthread implements Runnable {
        private RotamerSearch rs = null;
        private CommucObj cObj = null;
        private PrintStream logPS = null;
        String outputConfInfo = null;
        DACSthread(RotamerSearch rsP, CommucObj cObjP, PrintStream lP, String ociP) {
            rs = rsP;
            cObj = cObjP;
            logPS = lP;
            outputConfInfo = ociP;
        }
        public void run() {
            //Perform DACS
            doDACS(cObj.mutableSpots, rs, cObj.strandMut,
                   cObj.strandDefault, null,
                   cObj.prunedRot, cObj.algOption, cObj.sfFileIn, cObj.useSF, cObj.initEw, cObj.pruningE,
                   cObj.initDepth, 0, logPS, cObj.msp, cObj.numInitUnprunedConf,
                   cObj.diffFact, outputConfInfo, cObj.minRatioDiff, cObj.doMinimization, null,
                   null, false, cObj.params,
                   cObj.approxMinGMEC, cObj.lambda, cObj.partIndex, cObj, cObj.minimizeBB, cObj.numMutations,
                   cObj.scaleInt, cObj.maxIntScale, cObj.useEref, cObj.eRef, cObj.doBackrubs, cObj.backrubFile,
                   cObj.subDepth,cObj.mutRes2Strand,cObj.mutRes2StrandMutIndex,cObj.typeDep,cObj.addWTRot);
        }
    }

    //Compute the number of unpruned conformations
    //Compute the number of unpruned conformations
    private BigInteger compNumUnprunedConfs(int numMutable, PrunedRotamers<Boolean> prunedRotAtRes, int numRotForRes[]) {

        int numPrunedRotForRes[] = new int[numMutable]; //after pruning for this partition
        for(int i=0; i<numPrunedRotForRes.length; i++)
            numPrunedRotForRes[i] = 0;

        Iterator<RotInfo<Boolean>> iter = prunedRotAtRes.iterator();
        while(iter.hasNext()) {
            //for (int j=0; j<numMutable; j++){
            RotInfo<Boolean> ri = iter.next();
            if(ri.state) {
                numPrunedRotForRes[ri.curPos]++;
            }
            //for (int k=0; k<totalNumRotamers; k++){ //pruned rot are true and must be in the current set of allowed AA
            //if (prunedRotAtRes[j*totalNumRotamers + k])
            //	numPrunedRotForRes[j]++;
            //}
        }
        /*int numPrunedLigRot = 0;
        for (int k=0; k<numLigRotamers; k++){
        	if (prunedRotAtRes[numInAS*totalNumRotamers + k])
        		numPrunedLigRot++;
        }*/

        //count the total num confs and the num unpruned confs
        BigInteger numConfsTotal = new BigInteger("1");
        BigInteger numConfsUnpruned = new BigInteger("1");

        for (int j=0; j<numMutable; j++) {
            numConfsTotal = numConfsTotal.multiply(BigInteger.valueOf(numRotForRes[j]));
            numConfsUnpruned = numConfsUnpruned.multiply(BigInteger.valueOf(numRotForRes[j]-numPrunedRotForRes[j]));
        }
        /*if(ligPresent){
        	numConfsTotal = numConfsTotal.multiply(BigInteger.valueOf(numLigRotamers));
        	numConfsUnpruned = numConfsUnpruned.multiply(BigInteger.valueOf(numLigRotamers-numPrunedLigRot));
        }*/

        return numConfsUnpruned;
    }

    //Choose the splitting position (use only AS positions; the ligand is chosen not to be a splitting position)
    private int chooseSplitPosRandom(int numMutable) {
        Random randNum = new Random();
        return randNum.nextInt(numMutable);
    }

    //Choose the splitting position (use only AS positions; the ligand is chosen not to be a splitting position);
    //	Choose the AS residue position with the smallest fraction of pruned rotamers (from MinBounds and simple MinDEE)
    private int chooseSplitPos(int numMutable, PrunedRotamers<Boolean> prunedRotAtRes, int numRotForRes[],
                               int majorSplitPos[], int curDepth, double minRatioDiff) {

        final int minPartitions = 5; //the min number of rotamers that a splitting residue can have

        double pruneRatio[] = new double[numMutable];
        int minPos = -1;
        double minRatio = (double)Math.pow(10,38);
        for (int curRes=0; curRes<numMutable; curRes++) {

            if (!(isSplitRes(curRes,majorSplitPos,curDepth))) {
                if (numRotForRes[curRes]>=minPartitions) { //do not split at residues with very small number of rotamers
                    int curPruned = 0;
                    Iterator<RotInfo<Boolean>> iter = prunedRotAtRes.iterator(curRes);
                    RotInfo<Boolean> ri = iter.next();
                    while(ri.curPos == curRes && iter.hasNext()) {
                        //for (int curRot=0; curRot<numTotalRotamers; curRot++){
                        if (ri.state) { //prunedRotAtRes[curRes*numTotalRotamers+curRot]){ //cur rot is pruned (pruned rotamers are necessarily in the cur set of allowed AAs)
                            curPruned++;
                        }
                        ri = iter.next();
                    }
                    pruneRatio[curRes] = (double)curPruned/numRotForRes[curRes];
                    if (minRatio>=pruneRatio[curRes]) {
                        if ((minPos==-1)||(curRes<minPos)||(minRatio>=pruneRatio[curRes]+minRatioDiff)) {//preference to split at lower-numbered residues
                            minRatio = pruneRatio[curRes];
                            minPos = curRes;
                        }
                    }
                }
            }
        }

        if (minPos!=-1) {
            //System.out.println("minPos: "+minPos);
            //for (int i=0;i<numInAS;i++)System.out.print(pruneRatio[i]+" ");System.out.println();
            return minPos;
        } else //if split position not chosen, choose randomly
            return chooseSplitPosRandom(numMutable);
    }

    //Check if the residue curRes is one of the splitRes
    private boolean isSplitRes(int curRes, int majorSplitPos[], int curDepth) {
        for (int i=0; i<=curDepth; i++) {
            if (curRes==majorSplitPos[i])
                return true;
        }
        return false;
    }

    //Compute the number of rotamers for each residue position (assign to numRotForRes[])
    private int [] compNumRotForRes(int numberMutable, RotamerSearch rs, int strandMut[][],int mutRes2Strand[],int mutRes2StrandMutIndex[]) {

        int numRotForRes[] = new int[numberMutable];
        //boolean ligPresent = (numLigRot==0); //ligand present
        int treeLevels = numberMutable;
        /*if (ligPresent)
        	treeLevels++;*/

        numRotForRes = new int[treeLevels];

        int curNumRot = 0;
        for (int curLevel=0; curLevel<treeLevels; curLevel++) {
            int str = mutRes2Strand[curLevel];
            int strResNum = strandMut[str][mutRes2StrandMutIndex[curLevel]];
            /*if ((ligPresent)&&(curLevel==(treeLevels-1))){ //the ligand level
            	curNumRot = numLigRot;
            }
            else {*/ //AS residue
            curNumRot = 0;
            for (int i=0; i<rs.strandRot[str].getNumAllowable(strResNum); i++) { //add the rot for all allowable AA at this residue
                int newRot = rs.getNumRot(str, strResNum, rs.strandRot[str].getIndexOfNthAllowable(strResNum,i));
                curNumRot += newRot;
            }
            //}
            numRotForRes[curLevel] = curNumRot;
        }
        return numRotForRes;
    }


    //Get the mapping between rotamer indices (into the pruning matrix) and the number of the
    //	current rotamer for the giveen residue; assumes sysLR in rs is valid (all allowables for the AS residues)
    private Index3 [] getIndexMap(int numPartitions, RotamerSearch rs, int curRes, int strandMut[][],
                                  int mutRes2Strand[],int mutRes2StrandMutIndex[]) {

        int str = mutRes2Strand[curRes];
        int strResNum = strandMut[str][mutRes2StrandMutIndex[curRes]];

        Index3 indexMap[] = new Index3[numPartitions];
        int indNum = 0;
        for (int AA=0; AA<rs.strandRot[str].getNumAllowable(strResNum); AA++) { //for each AA for the given AS residue
            int curAA = rs.strandRot[str].getIndexOfNthAllowable(strResNum,AA);
            int numRotForAA = rs.getNumRot(str, strResNum, curAA);

            for (int curRot=0; curRot<numRotForAA; curRot++) { //for each rot for the given AA
                indexMap[indNum] = new Index3(curRes,curAA,curRot);//curRes*numTotalRot + rotamerIndexOffset[curAA] + curRot;
                indNum++;
            }
        }
        return indexMap;
    }

    //Setup the file with name filename for output
    private PrintStream setupOutputFile(String fileName) {
        PrintStream logPS = null; //the output file for conf info
        try {
            FileOutputStream fileOutputStream = new FileOutputStream(fileName);
            BufferedOutputStream bufferedOutputStream = new BufferedOutputStream( fileOutputStream );
            logPS = new PrintStream( bufferedOutputStream );
        } catch (Exception ex) {
            System.out.println("ERROR: An exception occured while opening log file");
        }
        return logPS;
    }

    ///////////////////////////////////////////////////////////////////////////
//	End of DEE section
///////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////
//	Distributed DACS section
///////////////////////////////////////////////////////////////////////////
    /**
     * Handles the distribution of the DACS computation to the set of available processors
    */
    private void doDistrDACSMaster(String runName, int numMutable, RotamerSearch rs, int strandMut[][],
                                   String strandDefault[][], String rotFile,
                                   PrunedRotamers<Boolean> prunedRotAtRes, int algOption, String sfFile, boolean useFlags, float initEw, float pruningE,
                                   int initDepth, int majorSplitPos[], BigInteger numInitUnprunedConfs,
                                   int diffFact, String outputPruneInfo, String outputConfInfo, double minRatioDiff, boolean doMinimize, String minPEM,
                                   String maxPEM, ParamSet sParams,
                                   boolean approxMinGMEC, float lambda, int numRotForRes[], OneMutation resumeResults[], String resumeFileName, boolean minimizeBB, int numMaxMut,
                                   boolean scaleInt, float maxIntScale, boolean useEref, float eRef[][], boolean doBackrubs, String backrubFile, int subDepth,
                                   int mutRes2Strand[],int mutRes2StrandMutIndex[],
                                   boolean [] strandPresent, String[][] strandLimits, int strandsPresent, boolean addWTRot) {

        System.out.println("Starting DACS (distributed)");
        System.out.println("Forming DACS partitions..");

        Index3 indexMap[][] = new Index3[initDepth][];
        int numPartitions[] = new int[indexMap.length];
        int maxNumPartitions = 1;

        //map rotamer index to rot num for the splitting residues
        for (int i=0; i<indexMap.length; i++) {
            numPartitions[i] = numRotForRes[majorSplitPos[i]]; //the number of partitions
            indexMap[i] = getIndexMap(numPartitions[i],rs,majorSplitPos[i],strandMut,mutRes2Strand,mutRes2StrandMutIndex);
            maxNumPartitions *= numPartitions[i];
        }

        //get the partitions
        OneMutation mutArray[] = formDACSpartitions(maxNumPartitions, initDepth, indexMap, numPartitions, prunedRotAtRes, majorSplitPos);

        if (resumeResults!=null) { //remove completed partitions
            int curMut = 0;
            OneMutation tmpArray2[] = new OneMutation[mutArray.length];
            for (int i=0; i<mutArray.length; i++) {
                boolean partFound = true;
                for (int j=0; j<resumeResults.length; j++) {
                    partFound = true;
                    for (int k=0; k<initDepth; k++) {
                        if (mutArray[i].resMut[k]!=resumeResults[j].resMut[k]) {
                            partFound = false;
                            break;
                        }
                    }
                    if (partFound) //partition already computed
                        break;
                }
                if (!partFound) {
                    tmpArray2[curMut] = new OneMutation();
                    tmpArray2[curMut].mutNum = mutArray[i].mutNum;
                    tmpArray2[curMut].resMut = mutArray[i].resMut;
                    curMut++;
                }
            }
            OneMutation tmpArray[] = new OneMutation[curMut]; //trim the size of the partition array
            System.arraycopy(tmpArray2, 0, tmpArray, 0, curMut);
            mutArray = tmpArray;

            System.out.println("Number non-zero partitions after removing completed results: "+mutArray.length);

            //update the best energy so far
            for (int j=0; j<resumeResults.length; j++) {
                pruningE = (float)Math.min(pruningE, resumeResults[j].score.doubleValue());
            }
        }

        //output the rs object
        outputObject(rs,rotFile);

        //sort the partitions
        sortDACSpartitions(mutArray,initDepth,majorSplitPos,numPartitions,prunedRotAtRes,indexMap,numMutable,
                           strandMut,pruningE,initEw,rs);

        System.out.println();

        MutationManager mutMan = new MutationManager(runName,mutArray,false);
        mutMan.setStrandMut(strandMut);
        mutMan.setStrandDefault(strandDefault);
        mutMan.setStrandPresent(strandPresent);
        mutMan.setStrandLimits(strandLimits);
        mutMan.setStrandsPresent(strandsPresent);
        mutMan.setMutableSpots(numMutable);
        mutMan.setMut2Strand(rs.mutRes2Strand);
        mutMan.setMut2StrandMutIndex(rs.mutRes2StrandMutIndex);
        //mutMan.setLigType(ligType);
        mutMan.setarpFilenameMin(minPEM);
        mutMan.setarpFilenameMax(maxPEM);
        mutMan.setParams(sParams);
        mutMan.setStericThresh(stericThresh);
        mutMan.setSoftStericThresh(softStericThresh);
        mutMan.setMutableSpots(numMutable);
        //mutMan.setnumLigRotamers(numLigRotamers);
        mutMan.setComputeEVEnergy(true);
        mutMan.setDoMinimization(doMinimize);
        mutMan.setMinimizeBB(minimizeBB);
        mutMan.setDoBackrubs(doBackrubs);
        mutMan.setBackrubFile(backrubFile);
        mutMan.setCalculateVolumes(false);
        mutMan.setNumMutations(numMaxMut);
        mutMan.setInitEw(initEw);
        mutMan.setPruningE(pruningE);
        //mutMan.setLigPresent(ligPresent);
        mutMan.setUseSF(useFlags);
        mutMan.setPrunedRot(prunedRotAtRes);
        mutMan.setSfFile(sfFile);
        mutMan.setRotFile(rotFile);
        mutMan.setDistrDACS(true);
        mutMan.setDistrDEE(false);
        mutMan.setBestScore(new BigDecimal(pruningE)); //the best E initially is the pruningE read from the parameter file
        mutMan.setAlgOption(algOption);
        mutMan.setInitDepth(initDepth);
        mutMan.setSubDepth(subDepth);
        mutMan.setDiffFact(diffFact);
        mutMan.setMinRatioDiff(minRatioDiff);
        mutMan.setNumInitUnprunedConf(numInitUnprunedConfs);
        mutMan.setOutputPruneInfo(outputPruneInfo);
        mutMan.setOutputConfInfo(outputConfInfo);
        mutMan.setMSP(majorSplitPos);
        mutMan.setApproxMinGMEC(approxMinGMEC);
        mutMan.setLambda(lambda);
        mutMan.setDistDepDielect(distDepDielect);
        mutMan.setDielectConst(dielectConst);
        mutMan.setDoDihedE(doDihedE);
        mutMan.setDoSolvationE(doSolvationE);
        mutMan.setSolvScale(solvScale);
        mutMan.setVdwMult(softvdwMultiplier);
        mutMan.setScaleInt(scaleInt);
        mutMan.setMaxIntScale(maxIntScale);
        mutMan.setUseEref(useEref);
        mutMan.setEref(eRef);
        mutMan.setRotamerLibrary(EnvironmentVars.aaRotLib);

        mutMan.setIdealizeSC(Perturbation.idealizeSC);
        mutMan.setAddWTRot(addWTRot);

        try {
            handleDoMPIMaster(mutMan,mutArray.length);
        } catch (Exception e) {
            System.out.println("ERROR: "+e);
            System.exit(1);
        }

        //Delete the temporary rs file
        if (!(new File(rotFile)).delete()) {
            System.out.println("ERROR: cannot delete file "+rotFile);
            System.exit(1);
        }
    }

    // Distributed DACS Slave function
    private CommucObj doDistrDACSSlave(CommucObj cObj) {

        long startTime = System.currentTimeMillis();

        //boolean ligPresent = cObj.ligPresent;

        //Setup the molecule system
        Molecule m = new Molecule();
        m = setupMolSystem(m,cObj.params,cObj.strandPresent,cObj.strandLimits);

        RotamerSearch rs = (RotamerSearch)readObject(cObj.rotFileIn); //load the saved rs from master

        Perturbation.idealizeSC = cObj.idealizeSC;

        String fn = "";
        for (int i=0; i<cObj.partIndex.length; i++)
            fn += ("_"+cObj.partIndex[i]);
        String outputConfInfo = ("./conf_info/"+cObj.outputConfInfo+fn);
        PrintStream logPS = setupOutputFile("./conf_info/"+cObj.outputPruneInfo+fn);

        if(logPS == null) {
            System.out.println("ERROR: Please make folder conf_info!!!");
            return null;
        }

        float otherBestE = cObj.pruningE; //the best energy from other partitions

        Thread t = new Thread(new DACSthread(rs,cObj,logPS,outputConfInfo)); //create a new thread for performing the DACS search
        t.start();
        long waitTime = 300000; //five minutes
        while (t.isAlive()) { //DACS search thread is still running

            try {
                t.join(waitTime);   //wait for waitTime before the next interruption of the DACS search
            } catch (Exception e) {}

            float rsBestE = rs.getBestE(); //the best energy from the current partition

            //System.out.println("partition "+cObj.partIndex+": curBestE "+bestE+" rsBestE "+rsBestE);

            if (rsBestE<otherBestE) { //new best energy for the current partition; update

                CommucObj c[] = new CommucObj[1];
                c[0] = new CommucObj();
                c[0].pruningE = rsBestE;

                //System.out.println("partition "+cObj.partIndex+": sending update to main node..");

                try {
                    MPItoThread.Send(c, 0, 1, ThreadMessage.OBJECT, 0, updateTag);
                } catch (Exception e) {}; //send back updated best energy

                otherBestE = rsBestE;
            }

            //check if there are updates for the best energy from the other partitions
            try {
                //System.out.println("partition "+cObj.partIndex+": checking for update from main node..");

                float c[] = new float[1];
                while (MPItoThread.Iprobe(0, updateTag)!=null) { //new update message received

                    //System.out.println("partition "+cObj.partIndex+": update from main node received..");

                    MPItoThread.Recv(c, 0, 1, ThreadMessage.FLOAT, 0, updateTag);

                    //System.out.println("partition "+cObj.partIndex+": updateE "+c[0]+" curBestE: "+bestE);

                    rs.updateBestE(c[0]);
                    otherBestE = Math.min(otherBestE,c[0]);
                }
            } catch (Exception e) {};
        }

        logPS.flush();
        logPS.close();

        rs = null;
        cObj.prunedRot = null;//smaller object, for sending back

        long stopTime = System.currentTimeMillis();
        cObj.elapsedTime = Math.round((stopTime - startTime) / 1000.0f);

        //System.out.println("Partition "+cObj.partIndex+" done, time: "+cObj.elapsedTime/60.0f+" minutes; sending results to main node..");

        return cObj;
    }

    //Forms the DACS partitions based on the splitting residue positions
    private OneMutation [] formDACSpartitions(int maxNumPartitions, int initDepth, Index3 indexMap[][], int numPartitions[],
            PrunedRotamers<Boolean> prunedRotAtRes, int majorSplitPos[]) {

        OneMutation mutArray[] = new OneMutation[maxNumPartitions];
        Index3 curInd[] = new Index3[initDepth];
        int curMut[] = new int[1];
        curMut[0] = 0;

        formDACSpartitionsHelper(initDepth, indexMap, numPartitions, prunedRotAtRes, mutArray, curMut, 0, curInd);

        OneMutation tmpArray[] = new OneMutation[curMut[0]]; //trim the size of the partition array
        System.arraycopy(mutArray, 0, tmpArray, 0, curMut[0]);
        mutArray = tmpArray;

        System.out.print("Partitioning residues: ");
        for (int i=0; i<initDepth; i++)
            System.out.print(majorSplitPos[i]+" ");
        System.out.println();
        System.out.println("Number of non-zero partitions: "+curMut[0]);

        return mutArray;
    }

    //Determines all non-pruned partitions for DACS deistribution
    //Called by formDACSpartitions(.)
    private void formDACSpartitionsHelper(int initDepth, Index3 indexMap[][], int numPartitions[],
                                          PrunedRotamers<Boolean> prunedRotAtRes, OneMutation mutArray[], int curMut[], int curDepth, Index3 curInd[]) {

        if (curDepth>=initDepth) { //new partition
            mutArray[curMut[0]] = new OneMutation();
            mutArray[curMut[0]].mutNum = curMut[0];
            mutArray[curMut[0]].index = new Index3[initDepth];
            for (int i=0; i<initDepth; i++)
                mutArray[curMut[0]].index[i] = curInd[i];
            curMut[0]++;
        } else {
            for (int i=0; i<numPartitions[curDepth]; i++) {
                if (!prunedRotAtRes.get(indexMap[curDepth][i])) { //only consider non-pruned partitions for distribution
                    curInd[curDepth] = indexMap[curDepth][i];
                    formDACSpartitionsHelper(initDepth, indexMap, numPartitions, prunedRotAtRes, mutArray, curMut, curDepth+1, curInd);
                }
            }
        }
    }

    //Sorts the DACS partitions by lower energy bounds;
    //The sorted array is returned in mutArray[]
    private void sortDACSpartitions(OneMutation mutArray[], int initDepth, int majorSplitPos[], int numPartitions[],
                                    PrunedRotamers<Boolean> prunedRotAtRes, Index3 indexMap[][], int numMutable, int strandMut[][],
                                    float pruningE, float initEw, RotamerSearch rsP) {

        RotamerSearch rs = rsP; //no changes should be made to the original RotamerSearch object
        rsP = null;

        System.out.print("Computing a lower bound on the conformational energy for each partition..");
        for (int m=0; m<mutArray.length; m++) { //for each partition

            if(m%100 == 0)
                System.out.println("Starting Partition.. "+m+" out of "+mutArray.length+" ");

            PrunedRotamers<Boolean> prunedForPartition = new PrunedRotamers<Boolean>(prunedRotAtRes,false); //no changes should be made to prunedRotAtRes[]
            Iterator<RotInfo<Boolean>> iter1 = prunedRotAtRes.iterator();
            while(iter1.hasNext()) {
                RotInfo<Boolean> ri = iter1.next();
                prunedForPartition.set(ri, new Boolean(ri.state));
            }
            //System.arraycopy(prunedRotAtRes, 0, prunedForPartition, 0, prunedRotAtRes.length);

            mutArray[m].setSortScores(true); //sort by lower bounds

            //artificially set to pruned all rotamers at the splitting position, other than the current partitioning rotamer for the given partitioning position
            for (int i=0; i<initDepth; i++) { //first, set all other rotamers for the partitioning positions to pruned
                Index3 curPart = mutArray[m].index[i];
                for (int j=0; j<numPartitions[i]; j++) {
                    Index3 curInd = indexMap[i][j];
                    if (curInd!=curPart) { //not the rotamer for the current partition
                        if (!prunedForPartition.get(curInd)) //rotamer not already pruned
                            prunedForPartition.set(curInd, true);
                    }
                }
            }

            //compute a lower bound on the conformational energies for this partition
            rs.DoMinBounds(numMutable,strandMut,pruningE,prunedForPartition,initEw, false, false, true);
            mutArray[m].score = new BigDecimal(rs.getBoundForPartition());
        }
        System.out.println("done");

        //sort the partitions
        System.out.print("Sorting partitions by their lower energy bounds..");
        //RyanQuickSort rqs = new RyanQuickSort();
        //rqs.Sort(mutArray);
        //rqs = null;
        Arrays.sort(mutArray);
        System.out.println("done");
        System.out.println("MinBoundOfMinPartition: "+mutArray[0].score);
    }

///////////////////////////////////////////////////////////////////////////
//	End of Distributed DACS section
///////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////////////////////
//	Distributed DEE section
// WARNING: the Distributed DEE section is outdated and must be carefully checked before called
///////////////////////////////////////////////////////////////////////////////////////////////
    //Handles the distribution of the DEE computation to the slave nodes
    private void doDistrDEEMaster(int numMutable, int strandMut[][],
                                  ParamSet sParams, String strandDefault[][],
                                  String minPEM, String maxPEM, float initEw, PrunedRotamers<Boolean> prunedRotAtRes,
                                  boolean doMinimize, RotamerSearch rs, String sfFile, String rotFile, boolean useSF,
                                  int numSpPos, int typeDEE, boolean minimizeBB, boolean scaleInt, float maxIntScale,
                                  boolean useEref, float eRef[][], boolean addWTRot) {

        //the total number of residues (active site + ligand, if present)
        int totalNumRes = numMutable;
        /*if (ligPresent) //ligand is present
        	totalNumRes++;*/

        // Generate all combinations
        int numMutRes = 1; //singles
        if (typeDEE==optPairs)
            numMutRes = 2; //pairs criterion
        int numComb = factorial(totalNumRes).divide(factorial(totalNumRes-numMutRes).multiply(factorial(numMutRes))).intValue();
        int residueMutatable[][] = new int[numComb][totalNumRes];
        generateCombinations(residueMutatable,totalNumRes,numMutRes);

        OneMutation mutArray[] = new OneMutation[numComb];
        for (int curMut=0; curMut<mutArray.length; curMut++) {
            mutArray[curMut] = new OneMutation();
            mutArray[curMut].mutNum = curMut;
            mutArray[curMut].resMut = new int[totalNumRes];
            for (int curRes=0; curRes<totalNumRes; curRes++)
                mutArray[curMut].resMut[curRes] = residueMutatable[curMut][curRes];
        }

        boolean splitFlags[][][][][][] = null;
        splitFlags = rs.getSplitFlags(splitFlags);//get the current DE pairs
        outputObject(splitFlags,sfFile);
        int numberOfStrands = strandMut.length;
        String AAallowed[][] = new String[numberOfStrands][]; //the AA's to which each residue can mutate
        for(int str=0; str<numberOfStrands; str++) {
            AAallowed[str] = new String[strandMut[str].length];
        }
        for(int str=0; str<numberOfStrands; str++) {
            for (int i=0; i<strandMut[str].length; i++)
                AAallowed[str][i] = (String)sParams.getValue("RESALLOWED"+str+"_"+i);
        }
        MutationManager mutMan = new MutationManager(null,mutArray,false);
        mutMan.setNumOfStrands(strandMut.length);
        mutMan.setStrandMut(strandMut);
        mutMan.setStrandDefault(strandDefault);
        //mutMan.setLigType(ligType);
        mutMan.setarpFilenameMin(minPEM);
        mutMan.setarpFilenameMax(maxPEM);
        mutMan.setParams(sParams);
        mutMan.setStericThresh(stericThresh);
        mutMan.setSoftStericThresh(softStericThresh);
        mutMan.setMutableSpots(numMutable);
        //mutMan.setnumLigRotamers(numLigRotamers);
        mutMan.setComputeEVEnergy(true);
        mutMan.setDoMinimization(doMinimize);
        mutMan.setMinimizeBB(minimizeBB);
        mutMan.setCalculateVolumes(false);
        mutMan.setInitEw(initEw);
        //mutMan.setLigPresent(ligPresent);
        mutMan.setUseSF(useSF);
        mutMan.setSpFlags(splitFlags);
        mutMan.setPrunedRot(prunedRotAtRes);
        mutMan.setSfFile(sfFile);
        mutMan.setRotFile(rotFile);
        mutMan.setDistrDACS(false);
        mutMan.setDistrDEE(true);
        mutMan.setAAallowed(AAallowed);
        mutMan.setNumSpPos(numSpPos);
        mutMan.setTypeDEE(typeDEE);
        mutMan.setDistDepDielect(distDepDielect);
        mutMan.setDielectConst(dielectConst);
        mutMan.setDoDihedE(doDihedE);
        mutMan.setDoSolvationE(doSolvationE);
        mutMan.setSolvScale(solvScale);
        mutMan.setVdwMult(softvdwMultiplier);
        mutMan.setScaleInt(scaleInt);
        mutMan.setMaxIntScale(maxIntScale);
        mutMan.setUseEref(useEref);
        mutMan.setEref(eRef);
        mutMan.setRotamerLibrary(EnvironmentVars.aaRotLib);

        mutMan.setDoPerturbations(rs.doPerturbations);
        mutMan.setPertFile(rs.pertFile);
        mutMan.setMinimizePerts(rs.minimizePerturbations);
        mutMan.setAddWTRot(addWTRot);
        mutMan.setIdealizeSC(Perturbation.idealizeSC);

        try {
            handleDoMPIMaster(mutMan,mutArray.length);
        } catch (Exception e) {
            System.out.println("ERROR: "+e);
            System.exit(1);
        }
    }

    // Distributed DEE Slave function
    private CommucObj doDistrDEESlave(CommucObj cObj) {

        long startTime = System.currentTimeMillis();

        //boolean ligPresent = cObj.ligPresent;

        //Setup the molecule system
        Molecule m = new Molecule();
        m = setupMolSystem(m,cObj.params,cObj.strandPresent,cObj.strandLimits);

        RotamerSearch rs = new RotamerSearch(m,cObj.mutableSpots,cObj.strandsPresent, hElect, hVDW, hSteric, true,
                                             true, cObj.epsilon, cObj.stericThresh, cObj.softStericThresh, cObj.distDepDielect,
                                             cObj.dielectConst, cObj.doDihedE, cObj.doSolvationE, cObj.solvScale, cObj.vdwMult, grl,
                                             cObj.doPerturbations, cObj.pertFile, cObj.minimizePerts, false, false);

        System.out.print("Loading precomputed energy matrix...");

        rs.loadPairwiseEnergyMatrices(new String(cObj.arpFilenameMin+".dat"),true);
        if (cObj.doMinimization)
            rs.loadPairwiseEnergyMatrices(new String(cObj.arpFilenameMax+".dat"),false);
        System.out.println("done");

        if (cObj.useEref) { //add the reference energies to the min (and max) intra-energies
            rs.addEref(cObj.eRef, cObj.doMinimization, cObj.strandMut);
        }


        int totalNumRes = cObj.mutableSpots;
        /*if (ligPresent)
        	totalNumRes++;*/

        //determine the two residues in the pair (for pairs DE) or the one residue (split-DEE)
        boolean resInMut[] = new boolean[totalNumRes];
        for (int i=0; i<totalNumRes; i++) {
            if (cObj.resMut[i]==1)
                resInMut[i] = true;
            else
                resInMut[i] = false;
        }


        //Set the allowable AA for each residue;
        // 		the ligand allowable set in the RotamerSearch() constructor
        boolean addWT = (new Boolean((String)cObj.params.getValue("ADDWT", "true"))).booleanValue();
        if(!addWT)
            checkWT(cObj.strandDefault, cObj.strandPresent, cObj.params);
        int molStrand = 0;
        for (int strNum=0; strNum<cObj.strandPresent.length; strNum++) {
            if(cObj.strandPresent[strNum]) {
                for (int k=0; k<cObj.strandMut[molStrand].length; k++) {
                    setAllowablesHelper(rs, cObj.params, addWT, strNum, molStrand, k, cObj.strandMut, cObj.strandDefault);
                }
                molStrand++;
            }
        }

        if(cObj.doPerturbations)
            rs.setupRCs(cObj.addWTRot);


        //Perform DEE pairs
        boolean splitFlags[][][][][][] = (boolean [][][][][][])readObject(cObj.sfFileIn);
        rs.setSplitFlags();
        rs.setSplitFlags(splitFlags);//initialize the split flags

        if (cObj.typeDEE==optPairs) { //simple Goldstein pairs
            rs.DoDEEPairs(cObj.mutableSpots, cObj.strandMut,
                          cObj.initEw, cObj.prunedRot, resInMut, cObj.doMinimization,
                          cObj.useSF, false, true, cObj.minimizeBB, cObj.scaleInt, cObj.maxIntScale,cObj.typeDep, false, 0.0f);
        } else if (cObj.typeDEE==optSplit) { //1- or 2-sp split-DEE
            //Precompute the MinDEE intervals
            rs.doCompMinDEEIntervals(cObj.mutableSpots, cObj.strandMut,
                                     cObj.prunedRot, cObj.scaleInt, cObj.maxIntScale);

            cObj.prunedRot = rs.DoDEEConfSplitting(cObj.mutableSpots,
                                                   cObj.strandMut, cObj.initEw, cObj.prunedRot, resInMut,
                                                   cObj.doMinimization, true, cObj.numSpPos, true, cObj.minimizeBB,cObj.typeDep, false, 0.0f);
        }

        long stopTime = System.currentTimeMillis();
        cObj.elapsedTime = Math.round((stopTime - startTime) / 1000.0f);


        //We cannot send back the whole splitFlags matrix, since it is too big; even just sending the
        //		newly-identified DE pairs can be infeasible for large systems, so we output the
        //		new splitFlags matrix to a temp file, which is then read by the master node
        if (cObj.typeDEE==optPairs) {
            cObj.sfFileOut = ("tmp_"+cObj.mutationNumber+"_"+stopTime);
            splitFlags = rs.getSplitFlags(splitFlags);
            outputObject(splitFlags,cObj.sfFileOut);
        }

        return cObj;
    }

///////////////////////////////////////////////////////////////////////////
//	End of Distributed DEE section
///////////////////////////////////////////////////////////////////////////

    private Object readObject(String inFile) {
        return readObject(inFile,true);
    }

    private Object readObject(String inFile, boolean repeat) {
        Object inObj = null;
        boolean done = false;
        while (!done) {
            try {
                ObjectInputStream in = new ObjectInputStream(new FileInputStream(inFile));
                inObj = in.readObject();
                in.close();
                done = true;
            } catch (Exception e) {
                //System.out.println(e.toString());
                //System.out.println("ERROR: An exception occurred while reading from object file");
                if (repeat)
                    done = false;
                else
                    done = true;
            }
        }
        return inObj;
    }

    static void outputObject(Object outObj, String outFile) {
        try {
            FileOutputStream fout = new FileOutputStream(outFile);
            ObjectOutputStream out = new ObjectOutputStream(fout);
            out.writeObject(outObj);
            out.close();
        } catch (Exception e) {
            System.out.println(e.toString());
            System.out.println("ERROR: An exception occurred while writing object file");
            System.exit(0);
        }
    }

    private int countPrunedRot(PrunedRotamers<Boolean> prunedRot) {
        int countPruned = 0;
        Iterator<RotInfo<Boolean>> iter = prunedRot.iterator();
        while(iter.hasNext()) {
            //for (int i=0; i<prunedRot.length; i++){
            RotInfo<Boolean> ri = iter.next();
            if (ri.state)
                countPruned++;
        }
        return countPruned;
    }

    private int countPrunedPairs(boolean prunedPairs[][][][][][]) {
        int countPruned = 0;

        boolean[][][][][][] fromMatrix = prunedPairs;
        for (int p1=0; p1<fromMatrix.length; p1++) {
            if (fromMatrix[p1]!=null) {
                for (int a1=0; a1<fromMatrix[p1].length; a1++) {
                    if (fromMatrix[p1][a1]!=null) {
                        for (int r1=0; r1<fromMatrix[p1][a1].length; r1++) {
                            if (fromMatrix[p1][a1][r1]!=null) {
                                for (int p2=0; p2<fromMatrix[p1][a1][r1].length; p2++) {
                                    if (fromMatrix[p1][a1][r1][p2]!=null) {
                                        for (int a2=0; a2<fromMatrix[p1][a1][r1][p2].length; a2++) {
                                            if (fromMatrix[p1][a1][r1][p2][a2]!=null) {
                                                for (int r2=0; r2<fromMatrix[p1][a1][r1][p2][a2].length; r2++) {
                                                    if(fromMatrix[p1][a1][r1][p2][a2][r2])
                                                        countPruned++;
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
        /*for (int i=0; i<prunedPairs.length; i++){
        	for (int j=i+1; j<prunedPairs.length; j++){
        		if (prunedPairs[i][j])
        			countPruned++;
        	}
        }*/
        return countPruned;
    }

    //Function taken from: http://www.javaworld.com/javaworld/javatips/jw-javatip76.html?page=2
    //Java Tip 76: An alternative to the deep copy technique
    //Author: Dave Miller
    static public Object deepCopy(Object oldObj) throws Exception {
        ObjectOutputStream oos = null;
        ObjectInputStream ois = null;
        try {
            ByteArrayOutputStream bos =
                new ByteArrayOutputStream(); // A
            oos = new ObjectOutputStream(bos); // B
            // serialize and pass the object
            oos.writeObject(oldObj);   // C
            oos.flush();               // D
            ByteArrayInputStream bin =
                new ByteArrayInputStream(bos.toByteArray()); // E
            ois = new ObjectInputStream(bin);                  // F
            // return the new object
            return ois.readObject(); // G
        } catch(Exception e) {
            System.out.println("Exception in ObjectCloner = " + e);
            throw(e);
        } finally {
            oos.close();
            ois.close();
        }
    }


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// MPI section
//		The tutorial "One-step Tutorial: MPI: It's easy to get started"
//			(http://www.lam-mpi.org/tutorials/one-step/ezstart.php ; accessed Oct 23, 2006) was used as MPI code reference
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //Setup and do MPI
    //	The parameter is only used by the master node
    public void handleDoMPI(String args[]) throws MPIException, InterruptedException {

        MPI.Init(args);

        int procRank = MPItoThread.Rank();
        numProc = MPItoThread.Size();

        //System.out.println("Node rank: "+procRank+" of "+numProc);
        MPI.COMM_WORLD.Barrier();

        setConfigPars();
        MPI.COMM_WORLD.Barrier();

        if (procRank==0) { //master node
            outputProgInfo();
            parse(args);

        } else { //slave node
            for(int i=0; i<mols.length; i++)
                mols[i] = null;

            handleDoMPISlave();
        }

        MPI.Finalize();
    }

    //Do MPI for the master node
    public void handleDoMPIMaster(MutationManager mutMan, int size) throws MPIException, InterruptedException {

        //KER: If it isn't an mpiRun start extra threads so that we can simulate
        //an mpiRun
        if(!mpiRun)
            MPItoThread.startThreads(this,Thread.currentThread());

        CommucObj cObjArray[] = new CommucObj[size];
        int numFinished = 0;

        int curMut = 0;
        for (int curProc=1; curProc<numProc; curProc++) { //distribute a single mutation per processor, for all processors

            if (curMut<cObjArray.length) { //more mutations to distribute

                System.out.println("Retrieving "+curMut+" of "+(cObjArray.length));
                cObjArray[curMut] = mutMan.getNextComObj(curMut);

                MPItoThread.Send(cObjArray, curMut, 1, ThreadMessage.OBJECT, curProc, regTag);
                curMut++;

                System.out.println("Sent to proc "+curProc);
                System.out.println();
            } else
                break;
        }

        boolean distrDACS = mutMan.getDistrDACS(); //distributed DACS computation

        while (numFinished<cObjArray.length) { //distribute and receive all remaining mutations

            CommucObj cObj[] = new CommucObj[1];
            cObj[0] = new CommucObj();

            //System.out.println("Receiving message on main node..");

            Object s = MPItoThread.Recv(cObj, 0, 1, ThreadMessage.OBJECT, ThreadMessage.ANY_SOURCE, ThreadMessage.ANY_TAG);

            //System.out.println("Received message on main node: tag "+s.tag+" source "+s.source);

            if (distrDACS) { //DACS computation, so check if the new energy is better than the best energy so far

                float curBestE = mutMan.getPruningE();

                if (cObj[0].pruningE<curBestE) { //the new energy is better

                    //System.out.println("Updating best energy on main node from "+curBestE+" to "+cObj[0].pruningE+", source (partition): "+cObj[0].partIndex);

                    mutMan.setBestScore(new BigDecimal(cObj[0].pruningE));
                    mutMan.setPruningE(cObj[0].pruningE);

                    float c[] = new float[1];
                    c[0] = cObj[0].pruningE;
                    for (int curProc=1; curProc<numProc; curProc++) { //update the best energy in each partition
                        MPItoThread.Isend(c, 0, 1, ThreadMessage.FLOAT, curProc, updateTag);
                    }
                }
            }

            if (MPItoThread.getStatusTag(s)==regTag) { //completed job
                mutMan.processFinishedMutation(cObj[0]);
                numFinished++;

                System.out.println("Finished: "+cObj[0].mutationNumber+", Time: "+(cObj[0].elapsedTime/60.0));

                if (curMut<cObjArray.length) {

                    System.out.print("Retrieving "+curMut+" of "+(cObjArray.length));
                    cObjArray[curMut] = mutMan.getNextComObj(curMut);

                    MPItoThread.Send(cObjArray, curMut, 1, ThreadMessage.OBJECT, MPItoThread.getStatusSource(s), regTag);
                    curMut++;

                    System.out.println(", Sent to proc "+MPItoThread.getStatusSource(s));
                    System.out.println();
                }
            }
        }
    }

    //Do MPI for a slave node
    public void handleDoMPISlave() throws MPIException, InterruptedException {

        int rank = MPItoThread.Rank();

        while (true) {

            Object s = MPItoThread.Probe(0, ThreadMessage.ANY_TAG);
            if (s!=null) {
                if (MPItoThread.getStatusTag(s)==regTag) { //new computation
                    CommucObj cObj[] = new CommucObj[1];
                    cObj[0] = new CommucObj();

                    //System.out.println("node "+rank+" receiving message from main node..");

                    MPItoThread.Recv(cObj, 0, 1, ThreadMessage.OBJECT, 0, regTag);

                    //System.out.println("node "+rank+" received message from main node..");

                    if (cObj[0]==null) //computation is done
                        return;

                    cObj[0] = handleKSSlave(cObj[0]); //perform computation
                    MPItoThread.Send(cObj, 0, 1, ThreadMessage.OBJECT, 0, regTag); //send back result
                } else { //(s.tag==updateTag), so discard
                    float c[] = new float[1];
                    MPItoThread.Recv(c, 0, 1, ThreadMessage.FLOAT, 0, updateTag);
                }
            }
        }
    }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	 End of MPI section
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
//	Backbone Flexibility Section
///////////////////////////////////////////////////////////////////////////
    public void generateBackbones(String s) {
        // Takes the following parameters
        // 1: System parameter filename (string)
        // 2: Backbone config filename (string)

        System.out.println("Performing Backbone Generation");

        ParamSet sParams = new ParamSet();
        sParams.addParamsFromFile(getToken(s,2)); //read system parameters
        sParams.addParamsFromFile(getToken(s,3)); //read mutation search parameters

        // Pull search parameters
        int numInAS = (new Integer((String)sParams.getValue("NUMINAS"))).intValue();
        String runName = ((String)sParams.getValue("RUNNAME"));
        boolean sysSampling = (new Boolean((String)sParams.getValue("SYSSAMPLING"))).booleanValue();
        double theta = (new Float((String)sParams.getValue("THETA"))).floatValue();
        double alpha = (new Float((String)sParams.getValue("ALPHA"))).floatValue();
        int numSamples = (new Integer((String)sParams.getValue("NUMSAMPLES"))).intValue();
        boolean ligPresent = (new Boolean((String)sParams.getValue("LIGPRESENT"))).booleanValue();
        String ligType = null;
        if (ligPresent)
            ligType = (String)(sParams.getValue("LIGTYPE"));

        if (theta%alpha!=0) {
            System.out.println("ERROR: Choose theta = k*alpha, for k - an integer.");
            System.exit(1);
        }


        MolParameters mp = loadMolecule(sParams, COMPLEX);
        //Setup the molecule system
        mp.m = new Molecule();
        mp.m = setupMolSystem(mp.m,sParams,mp.strandPresent,mp.strandLimits);

        int numberMutable = mp.numberMutable;
        int strandsPresent = mp.strandsPresent;
        int[][] strandMut = mp.strandMut;

        RotamerSearch rs = new RotamerSearch(mp.m,numberMutable,strandsPresent, hElect, hVDW, hSteric, true,
                                             true, 0.0f, stericThresh, softStericThresh, distDepDielect, dielectConst, doDihedE,
                                             doSolvationE, solvScale, softvdwMultiplier, grl, false, null, false, false, false);

        rs.doGenBackbones(runName, numberMutable, strandMut, theta, alpha, numSamples, sysSampling);
    }
///////////////////////////////////////////////////////////////////////////
//	End of Backbone Flexibility Section
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
//	Generate Random Conformations Section
///////////////////////////////////////////////////////////////////////////

    //Generates a random set of mutations/conformations for a given system
    public void generateRandConfs(String s) {
        System.out.println("RANDOM CONFS IS NOT IMPLEMENTED RIGHT NOW");
    }


    //Generates a random set of mutations/conformations for a given system
    /*public void generateRandConfs1(String s) {

    	// Takes the following parameters
    	// 1: System parameter filename (string)
    	// 2: Run name (for output files)
    	// 3: Ligand is present (boolean)
    	// 4: Ligand type (if present)
    	// 5: Number of conformations to be generated

    	ParamSet sParams = new ParamSet();
    	sParams.addParamsFromFile(getToken(s,2)); //read system parameters

    	// Pull search parameters
    	int numInAS = (new Integer((String)sParams.getValue("NUMINAS"))).intValue();
    	String runName = getToken(s,3);
    	boolean ligPresent = (new Boolean(getToken(s,4))).booleanValue();
    	String ligType = null;
    	if (ligPresent)
    		ligType = getToken(s,5);
    	int num = (new Integer(getToken(s,6))).intValue();

    	//Setup the molecule system
    	Molecule m = new Molecule();
    	setupMolSystem(m,sParams,ligPresent,ligType);

    	int residueMap[] = new int[numInAS];
    	String resDefault[] = new String[numInAS];
    	String resMapString = (String)sParams.getValue("RESIDUEMAP");
    	System.out.print("ResidueMap:");
    	for(int i=0;i<numInAS;i++){
    		int pdbResNum = (new Integer(getToken(resMapString,i+1))).intValue();
    		residueMap[i] = m.strand[sysStrNum].mapPDBresNumToStrandResNum(pdbResNum);
    		resDefault[i] = m.strand[sysStrNum].residue[residueMap[i]].name;
    		System.out.print(" "+residueMap[i]+"("+m.strand[sysStrNum].residue[residueMap[i]].fullName+")");
    	}
    	System.out.println();

    	PrintStream logPS = setupOutputFile(runName);

    	Random r = new Random();
    	for (int i=0; i<num; i++){
    		String AAs[] = new String[numInAS];
    		int rot[] = new int[numInAS+1];
    		for (int a=0; a<numInAS; a++){
    			AAs[a] = rl.getAAName(r.nextInt(numAAallowed));
    			int n = rl.getNumRotamers(AAs[a]);
    			if (n<=0)
    				n = 1;
    			rot[a] = r.nextInt(n);
    		}
    		if (ligPresent){
    			int n = grl.getNumRotamers(ligType);
    			if (n<=0)
    				n = 1;
    			rot[numInAS] = r.nextInt(n);
    		}

    		logPS.print(i+" ");
    		for (int a=0; a<numInAS; a++)
    			logPS.print(AAs[a]+" ");
    		if (ligPresent)
    			logPS.print(ligType+" ");

    		for (int a=0; a<numInAS; a++)
    			logPS.print(rot[a]+" ");
    		if (ligPresent)
    			logPS.print(rot[numInAS]+" ");

    		logPS.println();
    		logPS.flush();
    	}
    	logPS.close();
    }*/
///////////////////////////////////////////////////////////////////////////
//	End of Generate Random Conformations Section
///////////////////////////////////////////////////////////////////////////

    //KER: It is better to debug a threaded version of doDEE
    //than to continue to keep this function up to date.
    public void doSinglePairE(String s, ParamSet sParams) {

        // Takes the following parameters
        // 1: System parameter filename (string)
        // 2: Mutation search parameter filename (string)

        //Only read system and mutation files if sParams is null

        if (sParams==null) { //parameter files not read yet, so read them in
            sParams = new ParamSet();
            sParams.addParamsFromFile(getToken(s,2)); //read system parameters
            sParams.addParamsFromFile(getToken(s,3)); //read mutation search parameters
        }

        //int numInAS = (new Integer((String)sParams.getValue("NUMINAS"))).intValue();
        int numMutations = 2; //pairwise energies are computed
        boolean minimizeBB = (new Boolean((String)sParams.getValue("MINIMIZEBB"))).booleanValue();
        boolean doBackrubs = (new Boolean((String)sParams.getValue("DOBACKRUBS"))).booleanValue();
        String backrubFile = (String)sParams.getValue("BACKRUBFILE");
        String runName = (String)sParams.getValue("RUNNAME");
        String minEMatrixName = (String)sParams.getValue("MINENERGYMATRIXNAME",runName+"minM");
        String maxEMatrixName = (String)sParams.getValue("MAXENERGYMATRIXNAME",runName+"maxM");
        boolean templateAlwaysOn = (new Boolean((String)sParams.getValue("TEMPLATEALWAYSON"))).booleanValue();

        boolean ligPresent = (new Boolean((String)sParams.getValue("LIGPRESENT"))).booleanValue();
        boolean inputSysWithLig = ligPresent;
        String ligType = null;
        if (ligPresent)
            ligType = (String)sParams.getValue("LIGTYPE");

        boolean resumeSearch = (new Boolean((String)sParams.getValue("RESUMESEARCH"))).booleanValue();
        String resumeFilename = (String)sParams.getValue("RESUMEFILENAME");


        if( (new Boolean((String)sParams.getValue("DOPERTURBATIONS","false"))).booleanValue() ) {
            System.err.println("ERROR: DEEPer not supported in doSinglePairE");
            System.exit(1);
        }

        System.out.println("Run Name: "+runName);
        System.out.println("Precomputed Minimum Energy Matrix: "+minEMatrixName);
        System.out.println("Precomputed Maximum Energy Matrix: "+maxEMatrixName);
        System.out.println("Ligand Type: "+ligType);
        System.out.println("Num Residues Allowed to Mutate: "+numMutations);


        System.out.println("Computing _All_ Rotamer-Rotamer Energies");

        System.out.println("Starting minimum and maximum bound energy computation");

        if(resumeSearch) {
            System.out.println("** Resuming Search **");
            System.out.println("     resuming from file: "+resumeFilename);
        }

        MolParameters mp = loadMolecule(sParams, COMPLEX);
        Molecule m = mp.m;
        int numberMutable = mp.numberMutable;
        int strandsPresent = mp.strandsPresent;
        String[][] strandLimits = mp.strandLimits;
        boolean[] strandPresent = mp.strandPresent;
        int[][] strandMut = mp.strandMut;
        String[][] strandDefault = mp.strandDefault;

        RotamerSearch rs = new RotamerSearch(m,numberMutable,strandsPresent, hElect, hVDW, hSteric, true,
                                             true, 0.0f, stericThresh, softStericThresh, distDepDielect, dielectConst,
                                             doDihedE, doSolvationE, solvScale, softvdwMultiplier, grl,
                                             false, null, false, false, false);



        int resMut[] = new int[numberMutable];
        for (int i=0; i<resMut.length; i++)
            resMut[i] = 0;

        String flagMutType = "AS-AS";
        resMut[0] = 1;
        resMut[1] = 1;



        /*boolean useLig = ligPresent;

        if (((flagMutType.compareTo("AS-AS")==0)||(flagMutType.compareTo("SHL-AS")==0)||(flagMutType.compareTo("TEMPL")==0))&&(ligPresent)){
        	useLig = false;
        	m.deleteStrand(1);//we do not need the ligand for these runs
        	rs = new RotamerSearch(m, sysStrNum, -1, hElect, hVDW, hSteric, true,
        			true, 0.0f, stericThresh, softStericThresh, distDepDielect, dielectConst, doDihedE, doSolvationE, solvScale, softvdwMultiplier, rl, grl);
        }*/

        //System.out.println("Beginning setAllowables");
        // Ligand allowable set in the RotamerSearch() constructor
        boolean addWT = (new Boolean((String)sParams.getValue("ADDWT","true"))).booleanValue();
        if(!addWT)
            checkWT(strandDefault, strandPresent, sParams);
        int molStrand = 0;
        for (int strNum=0; strNum<strandMut.length; strNum++) {
            if(strandPresent[strNum]) {
                for (int k=0; k<strandMut[molStrand].length; k++) {
                    setAllowablesHelper(rs, sParams, addWT, strNum, molStrand, k, strandMut, strandDefault);
                }
                molStrand++;
            }
        }


        boolean shellRun = false;
        boolean intraRun = false;
        boolean templateOnly = false;

        if (flagMutType.compareTo("TEMPL")==0) {

            // **** Normal Active Site residue runs ****
            // Computes active site residue to active site residue pair energies
            shellRun = true;
            ligPresent = true;
            intraRun = false;
            templateOnly = true;
        } else if (flagMutType.compareTo("AS-AS")==0) {

            // **** Normal Active Site residue runs ****
            // Computes active site residue to active site residue pair energies
            shellRun = false;
            ligPresent = true;
            intraRun = false;
            templateOnly = false;
        } else if (flagMutType.compareTo("SHL-AS")==0) {

            // Then shell runs for the active site residues
            // Computes the active site residue rotamers to shell energies
            shellRun = true;
            ligPresent = true;
            intraRun = false;
            templateOnly = false;
        } else if (flagMutType.compareTo("INTRA")==0) {

            // Compute all intra-residue energies
            shellRun = false;
            intraRun = true;
            templateOnly = false;
        } else if (flagMutType.compareTo("LIG-AS")==0) {

            // **** Ligand present runs ****
            // This section computes the inter-residue energies between
            //  active site residues and the ligand
            shellRun = false;
            intraRun = false;
            templateOnly = false;
        } else { //(cObj.flagMutType.compareTo("LIG-SHL")==0)

            // Computes ligand rotamer to shell energies
            shellRun = true;
            intraRun = false;
            templateOnly = false;
        }

        // The goal is that the total energy of a system can be bounded by the sum of
        //  all pairwise active site residue entries plus the entry for each active site
        //  residue's shell run plus each active site residue's self intra-residue energy.
        //  If a ligand is present then one should add the ligand to shell energy, the
        //  ligand to each active site residue pairwise energy, and the ligand self intra-
        //  residue energy.

        //initialize the pairwise energy matrices (partial initialization - only for the residues involved in this computation, e.g., AS-AS)
        PairwiseEnergyMatrix minEmatrix = new PairwiseEnergyMatrix(numberMutable,resMut,strandMut,
                rs,shellRun,intraRun,false);
        PairwiseEnergyMatrix maxEmatrix = minEmatrix.copy();

        //Compute the corresponding matrix entries
        rs.simplePairwiseMutationAllRotamerSearch(strandMut,numberMutable,true,shellRun,intraRun,
                resMut,minEmatrix,maxEmatrix,minimizeBB,doBackrubs,templateOnly,backrubFile, templateAlwaysOn);

        return;
    }

    //Computes conformation energies for different combinations of the energy function parameters
    private void fitEparams(String s) {

        String firstParam = getToken(s,1);

        String sysFile = getToken(s,2);

        // Pull search parameters
        String confResFile = getToken(s,3);
        String runName = getToken(s,4);
        boolean ligPresent = (new Boolean(getToken(s,5))).booleanValue();
        String ligType = null;
        if (ligPresent)
            ligType = getToken(s,6);
        int numResults = (new Integer(getToken(s,7))).intValue();
        boolean minimizeBB = (new Boolean(getToken(s,8))).booleanValue();

        int numSteps = 10;
        double maxVdwMult = 1.05;
        double maxSolvScale = 0.3;

        solvScale = 0.0;
        double initVdwMult = 0.63;

        double vdwDelta = (maxVdwMult-initVdwMult)/numSteps;
        double solvDelta = (maxSolvScale-solvScale)/numSteps;

        for (int i1=0; i1<numSteps; i1++) {
            solvScale += solvDelta;
            for (int i2=0; i2<2; i2++) {
                if (i2==0)
                    distDepDielect = true;
                else
                    distDepDielect = false;

                for (int i3=0; i3<=numSteps; i3++) {
                    if (i3==0)
                        dielectConst = 1.0;
                    else
                        dielectConst = 4*i3;

                    softvdwMultiplier = initVdwMult-vdwDelta;
                    for (int i4=0; i4<=numSteps; i4++) {
                        softvdwMultiplier += vdwDelta;

                        String runNameParams = (runName+"_"+solvScale+"_"+distDepDielect+"_"+dielectConst+"_"+softvdwMultiplier);

                        String s1 = (firstParam+" "+sysFile+" "+confResFile+" "+runNameParams+" false none "+numResults+" "+minimizeBB);
                        handleMinDEEApplyRot(s1);

                        if (ligPresent) {
                            runNameParams = (runNameParams+"_lig");
                            s1 = (firstParam+" "+sysFile+" "+(confResFile+"_lig")+" "+runNameParams+" true "+ligType+" "+numResults+" "+minimizeBB);
                            handleMinDEEApplyRot(s1);
                        }
                    }
                }
            }
        }
    }


////////////////////////////////////////////////////////////////
// Compute Residue Entropy Section
////////////////////////////////////////////////////////////////
    /**
     * Handles the SCMF residue entropy computation.
     */
    public void handleDoResEntropy(String s, ParamSet sParams) {

        // Takes the following parameters
        // 1: System parameter filename (string)
        // 2: Mutation search parameter filename (string)

        //Only read system and mutation files if sParams is null

        if (sParams==null) { //parameter files not read yet, so read them in
            sParams = new ParamSet();
            sParams.addParamsFromFile(getToken(s,2)); //read system parameters
            sParams.addParamsFromFile(getToken(s,3)); //read mutation search parameters
        }

        String runName = (String)sParams.getValue("RUNNAME");

        String matrixName = (String)sParams.getValue("MATRIXNAME");
        boolean useEref = (new Boolean((String)sParams.getValue("USEEREF", "true"))).booleanValue();
        float dist = (new Float((String)sParams.getValue("DIST", "8.0"))).floatValue();
        String rotProbFile = (String)sParams.getValue("ROTPROBFILE", runName+"rotProb");
        float stericE = (new Float((String)sParams.getValue("STERICE", "30"))).floatValue();
        float maxPairE = (new Float((String)sParams.getValue("MAXPAIRE","1000.0"))).floatValue();

        MolParameters mp = new MolParameters();
        loadStrandParams(sParams, mp, COMPLEX);

        //Set nonprotein strands as not present
        for(int i=0; i<mp.numOfStrands; i++) {
            boolean isProtein = (new Boolean((String)sParams.getValue("STRANDAA"+i))).booleanValue();
            if(!isProtein) {
                mp.strandPresent[i] = false;
            }
        }

        //Setup the molecule system
        mp.m = new Molecule();
        mp.m = setupMolSystem(mp.m,sParams,mp.strandPresent,mp.strandLimits,false);
        Molecule m = mp.m;



        /*MolParameters mp = loadMolecule(sParams, COMPLEX);
        Molecule m = mp.m;
        int numberMutable = mp.numberMutable;
        int strandsPresent = mp.strandsPresent;
        String[][] strandLimits = mp.strandLimits;
        boolean[] strandPresent = mp.strandPresent;
        int[][] strandMut = mp.strandMut;
        String[][] strandDefault = mp.strandDefault;
        int[] mutRes2Strand = mp.mutRes2Strand;
        int [] mutRes2StrandMutIndex = mp.mutRes2StrandMutIndex;*/
        String[][] strandDefault = new String [mp.m.strand.length][];
        int[][] strandPDBnum = new int[mp.m.strand.length][];
        for(int str=0; str<strandDefault.length; str++) {
            strandDefault[str] = new String[m.strand[str].numberOfResidues];
            strandPDBnum[str] = new int[m.strand[str].numberOfResidues];
            for(int i=0; i<m.strand[str].numberOfResidues; i++) {
                strandDefault[str][i] = m.strand[str].residue[i].name;
                strandPDBnum[str][i] = m.strand[str].residue[i].getResNumber();
            }
        }

        int numRes = m.numberOfResidues;

        //KER: this function assumes we're only looking at one strand
        //KER: for now I'm assuming that's a protein strand
        RotamerLibrary rotLib = EnvironmentVars.aaRotLib;


        rotProbFile = (rotProbFile+".dat");
        float rotProb[][] = (float [][])readObject(rotProbFile,false);
        if (rotProb==null) { //Perform SCMF to compute the rotamer probabilities

            //read in or compute all of the energy matrices;
            //NOTE: backbone-to-backbone and rotamer-to-backbone energies are included in the pairwise rotamer energies
            float asasE[][][][] = getResEntropyEmatricesPair(matrixName, sParams, numRes, strandDefault, mp, runName, dist);
            float intraEnergies[][] = getResEntropyEmatricesIntra(matrixName, sParams, numRes, mp, runName);
            float eRef[] = getResEntropyEmatricesEref(useEref, null, null, null, intraEnergies, numRes,null,null,rotLib);

            rotProb = compRotProbSCMF(numRes, intraEnergies, asasE, eRef, rotProbFile, strandDefault, stericE, maxPairE,m,rotLib);
        }

        int numProx[] = new int[numRes]; //get the number of proximate residues for each residue position
        for (int i=0; i<numProx.length; i++)
            numProx[i] = 0;
        String asDistFile = matrixName+"_dist.dat";
        boolean as[][] = (boolean [][])readObject(asDistFile,false);
        for (int i=0; i<numRes; i++) {
            for (int j=i+1; j<numRes; j++) {
                if (as[i][j]) {
                    numProx[i]++;
                    numProx[j]++;
                }
            }
        }

        //m = null;

        PrintStream logPS = setupOutputFile(runName);

        logPS.print("resNum pdbResNum resDefault entropy"+" ");
        for (int j=0; j<rotLib.getAAtypesAllowed().length; j++)
            logPS.print(rotLib.getAAtypesAllowed()[j]+" ");
        logPS.println("numProx");
        logPS.flush();


        //Compute the AA probabilities for each residue position (as a function of the rotamer probabilities);
        //Compute the entropy at each position as a function of the amino acid probabilities for that position
        final double kB = 1.0;
        for (int i=0; i<numRes; i++) {
            int str=m.residue[i].strandNumber;
            int strResNum=m.residue[i].strandResidueNumber;
            if ( !strandDefault[m.residue[i].strandNumber][m.residue[i].strandResidueNumber].equalsIgnoreCase("PRO") ) {

                float aaProbBBE[] = new float[rotLib.getAAtypesAllowed().length];
                int curInd = 0;
                for (int j=0; j<aaProbBBE.length; j++) {

                    int numCurRot = rotLib.getNumRotamers(rotLib.getAAtypesAllowed()[j]);
                    if (numCurRot==0)
                        numCurRot = 1;

                    aaProbBBE[j] = 0.0f;
                    for (int r=0; r<numCurRot; r++) {
                        aaProbBBE[j] += rotProb[i][curInd];
                        curInd++;
                    }
                }

                //Compute the unnormalized AA probabilities as a weighted function of energies and PDB stats (if included)
                double aaProbUnNorm[] = new double[rotLib.getAAtypesAllowed().length];
                double aaNorm = 0.0;
                for (int j=0; j<rotLib.getAAtypesAllowed().length; j++) {

                    aaProbUnNorm[j] = 1.0;
                    aaProbUnNorm[j] *= aaProbBBE[j];

                    aaNorm += aaProbUnNorm[j];
                }

                //Normalize the probabilities
                double aaProb[] = new double[rotLib.getAAtypesAllowed().length];
                for (int j=0; j<rotLib.getAAtypesAllowed().length; j++) {
                    if (aaNorm!=0.0)
                        aaProb[j] = aaProbUnNorm[j]/aaNorm;
                    else
                        aaProb[j] = 0.0;
                }

                //Compute the entropy for the current residue position
                double sumAA = 0.0;
                for (int j=0; j<aaProb.length; j++) {
                    if (aaProb[j]>0.0)
                        sumAA += aaProb[j]*Math.log(aaProb[j]);
                }

                double entropy = -kB * sumAA;


                logPS.print(i+" "/*+defResNum[str][strResNum]*/+" "+strandDefault[str][strResNum]+" "+entropy+" ");
                for (int j=0; j<aaProb.length; j++)
                    logPS.print(aaProb[j]+" ");

                logPS.println(numProx[i]);
                logPS.flush();
            } else {
                logPS.println(i+" "/*+defResNum[str][strResNum]*/+" "+strandDefault[str][strResNum]+" "+0.0); //only for residue positions with wildtype Pro
            }
        }
        logPS.close();
    }

    //Computes the rotamer probabilities for all rotamers at all residue positions using SCMF
    private float[][] compRotProbSCMF(int numRes, float intraEnergies[][],
                                      float asasE[][][][], float eRef[], String rotProbFile, String strandDefault[][], float stericE, float maxPairE,Molecule m,
                                      RotamerLibrary rotLib) {

        final float constR = (float)(1.9891/1000.0);//the gas constant
        float T = 50000; //initial temperature
        final float endT = 298.15f; //the minimum temperature for annealing
        float tStepSize = 100.0f; //the temperature step size for annealing
        final float eps = 0.0001f; //the convergence threshold
        final float lambda = 0.5f; //scaling factor for updating the rotamer probabilities

        int totalNumRotamers = rotLib.getTotalNumRotamers();

        for (int i=0; i<asasE.length; i++) { //Set the max energy for any element in asasE[][][][] to maxPairE
            if (asasE[i]!=null) {
                for (int j=0; j<asasE[i].length; j++) {
                    if (asasE[i][j]!=null) {
                        for (int k=0; k<asasE[i][j].length; k++) {
                            if (asasE[i][j][k]!=null) {
                                for (int l=0; l<asasE[i][j][k].length; l++) {
                                    if (asasE[i][j][k][l]>maxPairE)
                                        asasE[i][j][k][l] = maxPairE;
                                }
                            }
                        }
                    }
                }
            }
        }

        int numPrunedRot = 0;
        boolean prunedRot[][] = new boolean[numRes][totalNumRotamers];
        for (int i=0; i<numRes; i++) {
            for (int j=0; j<totalNumRotamers; j++) {
                if ( (intraEnergies[1+i*totalNumRotamers+j][0]) > stericE) {
                    prunedRot[i][j] = true;
                    numPrunedRot++;
                } else
                    prunedRot[i][j] = false;
            }
        }
        System.out.println("Num rotamers pruned due to incompatibility with the template: "+numPrunedRot);


        //For each residue, compute the probability of each rotamer for that residue
        float Emf[][] = new float[numRes][totalNumRotamers];
        float rotProb[][] = new float[numRes][totalNumRotamers];
        float oldProb[][] = new float[numRes][totalNumRotamers];
        for (int i=0; i<numRes; i++) {
            for (int j=0; j<totalNumRotamers; j++) {
                if (!prunedRot[i][j])
                    rotProb[i][j] = 1.0f/totalNumRotamers;
                else
                    rotProb[i][j] = 0.0f;

                oldProb[i][j] = rotProb[i][j];
            }
        }

        while (T>=endT) { //perform annealing

            System.out.println("Starting run at T = "+T);

            boolean done = false;
            while (!done) {

                //Compute the new mean-field energy for each rotamer
                for (int i=0; i<numRes; i++) {

                    if ( !strandDefault[m.residue[i].strandNumber][m.residue[i].strandResidueNumber].equalsIgnoreCase("PRO") ) {

                        for (int j=0; j<totalNumRotamers; j++) {

                            if (!prunedRot[i][j]) {

                                Emf[i][j] = intraEnergies[1+i*totalNumRotamers+j][0] - eRef[getAAindFromRotNum(j,rotLib)];

                                if (asasE[i]!=null) {

                                    for (int k=0; k<numRes; k++) {
                                        if ( (k!=i) && (!strandDefault[m.residue[k].strandNumber][m.residue[k].strandResidueNumber].equalsIgnoreCase("PRO")) ) { //for all residues with which i_j has contact
                                            if ( (i<k) && (asasE[i][k]!=null) ) {
                                                for (int l=0; l<totalNumRotamers; l++) {
                                                    if (!prunedRot[k][l])
                                                        Emf[i][j] += asasE[i][k][j][l]*rotProb[k][l];
                                                }
                                            } else if ( (i>k) && (asasE[k][i]!=null) ) {
                                                for (int l=0; l<totalNumRotamers; l++) {
                                                    if (!prunedRot[k][l])
                                                        Emf[i][j] += asasE[k][i][l][j]*rotProb[k][l];
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

                //Update the rotamer probabilities
                for (int i=0; i<numRes; i++) {

                    if ( !strandDefault[m.residue[i].strandNumber][m.residue[i].strandResidueNumber].equalsIgnoreCase("PRO") ) {

                        float normFactor = 0.0f;
                        for (int j=0; j<totalNumRotamers; j++) {

                            if (!prunedRot[i][j])
                                normFactor += (float)Math.exp( -Emf[i][j] / (constR*T));
                        }

                        for (int j=0; j<totalNumRotamers; j++) {

                            if (!prunedRot[i][j]) {

                                oldProb[i][j] = rotProb[i][j]; //the probability before the update

                                if (normFactor!=0.0f)
                                    rotProb[i][j] = lambda*((float)Math.exp( -Emf[i][j] / (constR*T)) / normFactor) + (1-lambda)*oldProb[i][j];
                                else
                                    rotProb[i][j] = 0.0f;
                            }
                        }
                    }
                }

                float rms = checkRotProbConvEntropy(rotProb,oldProb,strandDefault,prunedRot,m);

                if (rms>eps)
                    done = false;
                else
                    done = true;
            }

            T -= tStepSize;
        }

        outputObject(rotProb,rotProbFile);

        return rotProb;
    }

    //Checks if the rotamer probabilities for the entropy computation have converged
    private float checkRotProbConvEntropy(float rotProb[][], float oldProb[][], String strandDefault[][], boolean prunedRot[][],Molecule m) {

        float sum = 0.0f;
        for (int i=0; i<rotProb.length; i++) {
            if ( !strandDefault[m.residue[i].strandNumber][m.residue[i].strandResidueNumber].equalsIgnoreCase("PRO") ) {
                for (int j=0; j<rotProb[i].length; j++) {
                    if (!prunedRot[i][j])
                        sum += (float)Math.pow( (rotProb[i][j]-oldProb[i][j]) , 2.0);
                }
            }
        }
        float rms = (float)Math.sqrt(sum);

        System.out.println("RMS: "+rms);

        return rms;
    }

    //Reads in (if computed) or computes the energy matrices for rot-to-rot pairwise energies;
    //This is not setup to handle energy minimization properly (mainly due to the approach used for computing only a small
    //		subset of the pairwise rot-to-rot energies), so doMinimize and minimizeBB should be false
    private float [][][][] getResEntropyEmatricesPair(String matrixName, ParamSet sParams, int numRes, String strandDefault[][], MolParameters mp,
            String runName, float dist) {

        String origPDB = (String)sParams.getValue("PDBNAME");

        //Check for the pairwise rot-to-rot file;
        String pairName = matrixName+"_pair.dat";

        float asasE[][][][] = readPairMatrixEntropy(pairName,numRes);
        if (asasE==null) { //compute the rot-to-template energy matrix

            //For each residue position i, get all residue positions that are within dist
            String asDistFile = matrixName+"_dist.dat";
            boolean as[][] = (boolean [][])readObject(asDistFile,false);
            if (as==null) {
                as = new boolean[numRes][numRes];
                computeEntropyEmatrixMaster(numRes,runName+".log",asDistFile,sParams,false,false,false,null,
                                            null,null,true,dist,as,false, mp); //the distances are returned in as[][]
            }

            int numPairs = 0;
            asasE = new float[numRes][numRes][][];

            for (int i=0; i<numRes; i++) {

                if (!strandDefault[mp.m.residue[i].strandNumber][mp.m.residue[i].strandResidueNumber].equalsIgnoreCase("PRO")) {
                    for (int j=i+1; j<numRes; j++) {
                        if (!strandDefault[mp.m.residue[j].strandNumber][mp.m.residue[j].strandResidueNumber].equalsIgnoreCase("PRO")) {
                            if (as[i][j]) {
                                asasE[i][j] = new float[1][];
                                numPairs++;
                            }
                        }
                    }
                }
            }

            computeEntropyEmatrixMaster(numRes,runName+".log",pairName,sParams,false,false,false,null,
                                        null,asasE,false,0.0f,null,false,mp);

            sParams.setValue("PDBNAME", origPDB);



            asasE = readPairMatrixEntropy(pairName,numRes);
        }

        return asasE;
    }

    //Reads in (if computed) or computes the energy matrices for intra-rot energies;
    //This is not setup to handle energy minimization properly (mainly due to the approach used for computing only a small
    //		subset of the pairwise rot-to-rot energies), so doMinimize and minimizeBB should be false
    private float [][] getResEntropyEmatricesIntra(String matrixName, ParamSet sParams, int numRes, MolParameters mp, String runName) {

        String origPDB = (String)sParams.getValue("PDBNAME");


        //Check for the intra energies file
        String intraName = matrixName+"_intra.dat";

        float intraEnergies[][] = (float [][])readObject(intraName,false); //check if already computed

        if (intraEnergies==null) //compute the intra-rotamer energy matrix
            computeEntropyEmatrixMaster(numRes,runName+".log",intraName,sParams,false,false,false,null,
                                        intraEnergies,null,false,0.0f,null,true,mp);

        sParams.setValue("PDBNAME", origPDB);


        /*if (ligStrNum>=0) //the ligand is not used here
        	m.deleteStrand(ligStrNum);*/

        intraEnergies = (float [][])readObject(intraName,false);

        return intraEnergies;
    }

    //Reads in the amino acid reference energies (if used);
    private float [] getResEntropyEmatricesEref(boolean useEref, float intraEnergies[][][][][][], StrandRotamers strandRot[], int strandMut[][], float intraEnergiesEntropy[][], int numRes,
            int mutRes2Strand[],int mutRes2StrandMutIndex[], RotamerLibrary rotLib) {

        if ( (intraEnergies==null && intraEnergiesEntropy==null) || (intraEnergies!=null && intraEnergiesEntropy!=null) ) { //exactly one of the two matrices should be non-null
            System.out.println("ERROR: exactly one matrix can be used for the reference energy computation.");
            System.exit(1);
        }

        float eRef[] = new float[rotLib.getAAtypesAllowed().length];

        if (useEref) { //use AA reference energies
            eRef = compEref(intraEnergies,strandRot,strandMut,intraEnergiesEntropy,numRes,mutRes2Strand, mutRes2StrandMutIndex,rotLib);
        } else {
            for(int j=0; j<eRef.length; j++)
                //for (int i=0; i<eRef[j].length; i++)
                eRef[j] = 0.0f;
        }

        return eRef;
    }

    //Reads in the pairwise energy matrix for the entropy computation
    private float [][][][] readPairMatrixEntropy(String fName, int numRes) {

        BufferedReader bufread = null;
        try {
            File file = new File(fName);
            FileReader fr = new FileReader(file);
            bufread = new BufferedReader(fr);
        } catch (FileNotFoundException e) {
            return null;
        }

        float asasE[][][][] = new float[numRes][][][];

        boolean done = false;
        String str = null;

        while (!done) {
            try {
                str = bufread.readLine();
            } catch ( Exception e ) {
                System.out.println("ERROR: An error occurred while reading input");
                System.exit(0);
            }

            if (str == null) // stop if we've reached EOF
                done = true;
            else if (getToken(str,1).charAt(0)=='%') //skip comment lines
                continue;
            else {
                int i = new Integer(getToken(str,1)).intValue();
                String name = getToken(str,2);
                asasE[i] = (float [][][])readObject(name,false);
                if (asasE[i]==null) {
                    System.out.println("ERROR: Could not read data from file "+name);
                    System.exit(1);
                }
            }
        }

        // We're done reading them in
        try {
            bufread.close();
        } catch(Exception e) {
        }

        return asasE;
    }

    //Determines how many residue positions in the molecule numbered from (pos+1...) (pos is strand-relative numbering)
    //		are within dist from residue position pos
    public boolean [] getProxAS(Molecule m, int pos, float dist, boolean as[]) {

        int strandMap[] = new int[2];
        strandMap[0] = pos;

        for (int i=pos+1; i<m.numberOfResidues; i++) {

            if (i!=pos) {

                boolean done = false;

                strandMap[1] = i;

                Molecule m1 = getASASMolEntropy(m, strandMap);

                StrandRotamers strandRot = new StrandRotamers(EnvironmentVars.aaRotLib,m1.strand[0]);

                for (int j=0; j<m1.numberOfResidues; j++) {
                    for (int r=0; r<strandRot.rl.getAAtypesAllowed().length; r++) {
                        strandRot.setAllowable(j,strandRot.rl.getAAtypesAllowed()[r]);
                        m1.residue[j].flexible = true;
                    }
                }

                for(int q1=0; q1<strandRot.getNumAllowable(0); q1++) {

                    int AAindex1 = strandRot.getIndexOfNthAllowable(0,q1);
                    strandRot.changeResidueType(m1,0,strandRot.rl.getAAName(AAindex1),true,true);

                    for(int q2=0; q2<strandRot.getNumAllowable(1); q2++) {

                        int AAindex2 = strandRot.getIndexOfNthAllowable(1,q2);
                        strandRot.changeResidueType(m1,1,strandRot.rl.getAAName(AAindex2),true,true);

                        int numRot1 = strandRot.rl.getNumRotForAAtype(AAindex1);

                        int w1 = 0;
                        if (numRot1<=0)
                            w1 = -1;

                        while ((w1<numRot1)&&(!done)) {

                            if (w1!=-1)
                                strandRot.applyRotamer(m1, 0, w1);

                            int numRot2 = strandRot.rl.getNumRotForAAtype(AAindex2);

                            int w2 = 0;
                            if (numRot2<=0)
                                w2= -1;

                            while ((w2<numRot2)&&(!done)) {

                                if (w2!=-1)
                                    strandRot.applyRotamer(m1, 1, w2);

                                Residue r1 = m1.residue[0];
                                Residue r2 = m1.residue[1];

                                if (r1.getDist(r2,true)<=dist) {
                                    as[i] = true;
                                    done = true;
                                }

                                w2++;
                            }

                            w1++;
                        }
                        if (done)
                            break;
                    }
                    if (done)
                        break;
                }
                if (!done)
                    as[i] = false;
            }
        }

        return as;
    }

    //Returns the AA index into rotamerIndexOffset to which rotNum belongs
    private int getAAindFromRotNum(int rotNum, RotamerLibrary rotLib) {
        int[] rotamerIndexOffset = rotLib.getRotamerIndexOffset();
        for (int i=0; i<rotamerIndexOffset.length-1; i++) {
            if ( (rotNum>=rotamerIndexOffset[i]) && (rotNum<rotamerIndexOffset[i+1]) )
                return i;
        }
        if (!(rotNum>=rotLib.getTotalNumRotamers()))
            return (rotamerIndexOffset.length-1);
        else
            return -1;
    }

    //Distributes the different types of energy computation for the entropy calculation
    private void computeEntropyEmatrixMaster(int numRes, String runName, String matrixName, ParamSet sParams, boolean doMinimize, boolean minimizeBB, boolean doBackrubs, String backrubFile,
            float bbEnergies[][], float asasE[][][][], boolean compASASdist, float dist, boolean asDist[][], boolean intraRun, MolParameters mp) {

        //MolParameters mp = loadMolecule(sParams, COMPLEX);
        //KER: This function assumes using the protein rotamer library
        RotamerLibrary rotLib = EnvironmentVars.aaRotLib;

        int mutEnerMatrixSize = 0;
        int residueMap[] = null;
        OneMutation mutArray[] = null;
        int numMutable = numRes;

        int numMut = 0;

        if (compASASdist) { //compute the min distance between any pair of rotamers for each pair of residue positions
            mutArray = new OneMutation[numMutable];
            for (int i=0; i<asDist.length; i++)
                mutArray[i] = new OneMutation();

            numMut = numRes;
        } else {
            if (intraRun) { //computing intra energies

                System.out.println("Starting intra-rot energy computation..");

                mutEnerMatrixSize = 1 + rotLib.getTotalNumRotamers()*numMutable;

                bbEnergies = new float[mutEnerMatrixSize][1];
                for(int i=0; i<mutEnerMatrixSize; i++) {
                    for(int j=0; j<1; j++) {
                        bbEnergies[i][j] = 0.0f;
                    }
                }

                numMutable = 1;
                residueMap = new int[numMutable];

                mutArray = new OneMutation[numRes];
                for (int i=0; i<mutArray.length; i++) {
                    mutArray[i] = new OneMutation();
                    mutArray[i].flagMutType = "INTRA";
                }

                numMut = numRes;
            }

            else { //AS-AS energies

                System.out.println("Starting rot-to-rot energy computation..");

                numMutable = 2;

                int numPairs = 0;
                for (int i=0; i<asasE.length; i++) {
                    for (int j=i+1; j<asasE[0].length; j++) {
                        if (asasE[i][j]!=null)
                            numPairs++;
                    }
                }
                mutArray = new OneMutation[numPairs];

                int curPair = 0;
                for (int i=0; i<asasE.length; i++) {
                    for (int j=i+1; j<asasE[0].length; j++) {
                        if (asasE[i][j]!=null) {
                            mutArray[curPair] = new OneMutation();
                            mutArray[curPair].flagMutType = "AS-AS";
                            mutArray[curPair].resMut = new int[2];
                            mutArray[curPair].resMut[0] = i; //molecule-relative numbering
                            mutArray[curPair].resMut[1] = j;
                            curPair++;

                            asasE[i][j] = new float[rotLib.getTotalNumRotamers()][rotLib.getTotalNumRotamers()];
                        }
                    }
                }

                numMut = numPairs;
            }
        }


        MutationManager mutMan = new MutationManager(runName,mutArray,true);
        //mutMan.setStrandMut(strandMut);
        mutMan.setStrandLimits(mp.strandLimits);
        mutMan.setStrandPresent(mp.strandPresent);
        mutMan.setStrandsPresent(mp.strandsPresent);
        //mutMan.setLigType(null);
        mutMan.setarpFilenameMin(matrixName);
        mutMan.setIntraEntropyMatrixMin(bbEnergies);
        mutMan.setParams(sParams);
        mutMan.setStericThresh(stericThresh);
        mutMan.setSoftStericThresh(softStericThresh);
        mutMan.setRotamerLibrary(rotLib);
        mutMan.setMutableSpots(numMutable);
        mutMan.setComputeEVEnergy(true);
        mutMan.setDoMinimization(doMinimize);
        mutMan.setMinimizeBB(minimizeBB);
        mutMan.setDoBackrubs(doBackrubs);
        mutMan.setBackrubFile(backrubFile);
        mutMan.setCalculateVolumes(false);
        mutMan.setLigPresent(false);
        mutMan.setDistDepDielect(distDepDielect);
        mutMan.setDielectConst(dielectConst);
        mutMan.setDoDihedE(doDihedE);
        mutMan.setDoSolvationE(doSolvationE);
        mutMan.setSolvScale(solvScale);
        mutMan.setVdwMult(softvdwMultiplier);
        mutMan.setEntropyComp(true);
        mutMan.setPairEntropyMatrix(asasE);
        mutMan.setASdistMatrix(asDist);
        mutMan.setASdist(dist);
        mutMan.setCompASdist(compASASdist);
        mutMan.setRotamerLibrary(EnvironmentVars.aaRotLib);


        try {
            handleDoMPIMaster(mutMan,numMut);
        } catch (Exception e) {
            System.out.println("ERROR: "+e);
            System.exit(1);
        }

        if (compASASdist) {
            asDist = mutMan.getASdistMatrix();
            outputObject(asDist,matrixName);
        } else {
            if (intraRun) {
                bbEnergies = mutMan.getMinEmatrixEntropy();

                int numCompEntries = 0;
                for(int i1=0; i1<mutEnerMatrixSize; i1++) {
                    if ((bbEnergies[i1][0]!=0.0f))
                        numCompEntries++;
                }
                System.out.println("Num computed entries: "+numCompEntries);

                outputObject(bbEnergies,matrixName);
            } else {
                asasE = mutMan.getPairEntropyEmatrix();

                PrintStream logPS = setupOutputFile(matrixName);
                for (int i=0; i<asasE.length; i++) {
                    if (asasE[i]!=null) {
                        String fn = ("peme/pem_entr_"+i);
                        logPS.println(i+" "+fn);
                        outputObject(asasE[i],fn);
                        logPS.flush();
                    }
                }
                logPS.close();
            }
        }
    }

    private CommucObj handleDoResEntropySlave(CommucObj cObj) {
        //TODO: understand what this function does and fix it for multiple strands
        long startTime = System.currentTimeMillis();

        //Set nonprotein strands as not present
        for(int i=0; i<cObj.strandPresent.length; i++) {
            boolean isProtein = (new Boolean((String)cObj.params.getValue("STRANDAA"+i))).booleanValue();
            if(!isProtein) {
                cObj.strandPresent[i] = false;
            }
        }
        //Setup the molecule system
        Molecule m = new Molecule();
        //BAD CODE setupMolSystem(m,cObj.params,false,null); //the ligand is not used here
        m = setupMolSystem(m,cObj.params,cObj.strandPresent,cObj.strandLimits, false);

        if (cObj.compASdist) { //AS-AS distance computation
            cObj.asDist = getProxAS(m,cObj.mutationNumber,cObj.dist,cObj.asDist);
            cObj.compEE = new SamplingEEntries[0];
        } else { //AS-AS or INTRA energy computation

            int numMutable = cObj.mutableSpots;

            int resMut[] = new int[numMutable];
            int strandMut[][] = new int[cObj.strandMut.length][];

            for(int str=0; str<cObj.strandMut.length; str++)
                strandMut[str] = new int[cObj.strandMut[str].length];

            boolean shellRun = false;
            boolean ligPresent = false;
            boolean intraRun = false;
            boolean templateOnly = false;

            if (cObj.flagMutType.compareTo("AS-AS")==0) { //AS-AS run

                for (int i=0; i<cObj.mutableSpots; i++)
                    resMut[i] = 1;

                m = getASASMolEntropy(m,cObj.strandMut[0]);

                strandMut[0][0] = 0;
                strandMut[0][1] = 1;

            }

            else if (cObj.flagMutType.compareTo("INTRA")==0) {

                intraRun = true;

                m = getASASMolEntropy(m,cObj.strandMut[0]);

                strandMut = new int[1][1];
                strandMut[0][0] = 0;

                resMut = new int[1];
                resMut[0] = 1;
            }

            else {
                System.out.println("ERROR: only AS-AS and INTRA runs allowed for the pairwise entropy matrix precomputation.");
                System.exit(1);
            }

            //KER: there will only be one strand with the remade molecule
            int strandsPresent = 1;

            RotamerSearch rs = new RotamerSearch(m, cObj.mutableSpots, strandsPresent, hElect, hVDW, hSteric, true,
                                                 true, cObj.epsilon, cObj.stericThresh, cObj.softStericThresh, cObj.distDepDielect,
                                                 cObj.dielectConst, cObj.doDihedE,cObj.doSolvationE,cObj.solvScale,cObj.vdwMult, grl,
                                                 false, null, false, false, false);

            for(int j=0; j<cObj.mutableSpots; j++) {
                for(int q=0; q<rs.strandRot[0].rl.getAAtypesAllowed().length; q++)
                    rs.setAllowable(strandMut[0][j],rs.strandRot[0].rl.getAAtypesAllowed()[q],0);
            }

            //initialize the pairwise energy matrices (partial initialization - only for the residues involved in this computation, e.g., AS-AS)
            PairwiseEnergyMatrix minEmatrix = new PairwiseEnergyMatrix(numMutable,resMut,strandMut,rs,shellRun,intraRun,false);
            PairwiseEnergyMatrix maxEmatrix = minEmatrix.copy();

            rs.simplePairwiseMutationAllRotamerSearch(strandMut,numMutable,cObj.doMinimization,shellRun,intraRun,
                    resMut,minEmatrix,maxEmatrix,cObj.minimizeBB,cObj.doBackrubs,templateOnly,cObj.backrubFile, false);

            long stopTime = System.currentTimeMillis();
            cObj.elapsedTime = Math.round((stopTime - startTime) / 1000.0f);

            //Store the information in less space to allow the master node to buffer several cObj at once
            cObj.compEE = minEmatrix.generateCompEE(maxEmatrix);
        }

        return cObj;
    }

    //Returns a molecule that contains only the residues in the system strand (sysStrNum) of molecule m that are specified by residueMap[];
    //	This function is used for the pairwise energy matrix computation in the residue entropy calculations
    private Molecule getASASMolEntropy (Molecule m, int strandMap[]) {

        Molecule m1 = new Molecule();

        for (int i=0; i<strandMap.length; i++) {

            Residue oldResidue = m.residue[strandMap[i]];

            Residue newResidue = new Residue();
            newResidue.name = oldResidue.name;
            newResidue.fullName = oldResidue.fullName;

            for (int j=0; j<oldResidue.numberOfAtoms; j++) {

                Atom oldAtom = oldResidue.atom[j];

                Atom newAtom = new Atom(oldAtom.name,oldAtom.coord[0],oldAtom.coord[1],oldAtom.coord[2]);
                newAtom.modelAtomNumber = oldAtom.modelAtomNumber;
                newAtom.strandNumber = oldAtom.strandNumber;
                newAtom.elementType = oldAtom.elementType;
                newResidue.addAtom(newAtom);
            }

            m1.addResidue(0,newResidue);
        }
        //Determine the bonds between the atoms in the molecule
        m1.determineBonds();

        // Assign the molecule relative atom numbers
        m1.updateMoleculeAtomNumbers();

        m1.strand[0].isProtein = true;

        return m1;
    }

    //Computes the amino acid reference energies using the intra-rotamer energies from intraEnergies or intraEnergiesEntropy;
    //For each amino acid type, takes the min energy among all rotamers for that amino acid type, for all numRes residues
    private float [] compEref(float intraEnergies[][][][][][], StrandRotamers strandRot[], int strandMut[][], float intraEnergiesEntropy[][], int numRes,
                              int mutRes2Strand[],int mutRes2StrandMutIndex[], RotamerLibrary rotLib) {

        int numAAallowed = rotLib.getNumAAallowed();

        float bigE = (float)Math.pow(10,38);
        float eRef[] = new float[numAAallowed];
        for(int j=0; j<numAAallowed; j++)
            //for (int i=0; i<eRef[j].length; i++)
            eRef[j] = bigE;

        int ind = 1; //skip the entry [0][0], since this is the fixed template energy



        for (int i=0; i<numRes; i++) {
            int str = -1;
            int strResNum = -1;
            if(intraEnergies!=null) {
                str = mutRes2Strand[i];
                strResNum = strandMut[str][mutRes2StrandMutIndex[i]];
            }
            int numAA = numAAallowed;
            if (intraEnergies!=null) { //the six-dimensional, so the energies for only a subset of the amino acid types are available
                numAA = strandRot[str].getNumAllowable(strResNum);
            }
            for (int j=0; j<numAA; j++) {
                int aaInd = j;
                if (intraEnergies!=null)
                    aaInd = strandRot[str].getIndexOfNthAllowable(strResNum,j);
                int numRot = rotLib.getNumRotForAAtype(aaInd);
                if (numRot==0) //ALA or GLY
                    numRot = 1;
                float curMin = bigE;
                for (int k=0; k<numRot; k++) {
                    if (intraEnergies!=null)
                        curMin = Math.min(curMin,intraEnergies[i][aaInd][k][i][0][0]);
                    else
                        curMin = Math.min(curMin,intraEnergiesEntropy[ind][0]);
                    ind++;
                }
                eRef[aaInd] = Math.min(eRef[aaInd],curMin);
            }
        }

        //for(int j=0; j<numRes;j++)
        for (int i=0; i<eRef.length; i++)
            if (eRef[i]==bigE)
                eRef[i] = 0.0f;

        return eRef;
    }
//////////////////////////////////////////////////////
// End Compute Residue Entropy Section
//////////////////////////////////////////////////////

    private void selectResidues(String s) {

        // Takes the following parameters
        // 1: System parameter filename (string)
        // 2: Residue search filename (string)

        ParamSet sParams = new ParamSet();
        sParams.addParamsFromFile(getToken(s,2)); //read system parameters
        sParams.addParamsFromFile(getToken(s,3)); //read mutation search parameters

        String runName = (String)sParams.getValue("RUNNAME");
        int numRes = (new Integer((String)sParams.getValue("NUMRES"))).intValue();
        int pdbRes[] = new int[numRes];
        float dist[] = new float[numRes];
        boolean ligPresent = (new Boolean((String)sParams.getValue("LIGPRESENT"))).booleanValue();
        String ligType = null;
        if (ligPresent)
            ligType = (String)(sParams.getValue("LIGTYPE"));

        String resString = (String)sParams.getValue("RESIDUES");
        String distString = ((String)sParams.getValue("DIST"));
        for (int i=0; i<numRes; i++) {
            pdbRes[i] = new Integer((String)getToken(resString,i+1)).intValue();
            dist[i] = new Float((String)getToken(distString,i+1)).floatValue();
        }

        MolParameters mp = new MolParameters();
        loadStrandParams(sParams, mp, COMPLEX);

        Molecule m = new Molecule();
        m = setupMolSystem(m,sParams,mp.strandPresent,mp.strandLimits);

        //Map from the pdb residue numbers to the residue index in m.residue[]
        int residues[] = new int[numRes];
        int curRes = 0;
        for (int i=0; i<m.numberOfResidues; i++) {
            for (int j=0; j<numRes; j++) {
                if (m.residue[i].getResNumber()==pdbRes[j]) {
                    residues[curRes] = i;
                    curRes++;
                    break;
                }
            }
        }


        boolean asProx[] = new boolean[m.numberOfResidues];
        for (int i=0; i<asProx.length; i++)
            asProx[i] = false;

        for (int res=0; res<numRes; res++) {
            Residue r1 = m.residue[residues[res]];

            for (int i=0; i<asProx.length; i++) {

                if (i!=residues[res]) {

                    Residue r2 = m.residue[i];

                    if (r1.getDist(r2,true)<=dist[res])
                        asProx[i] = true;
                } else
                    asProx[i] = true;
            }
        }

        Molecule m1 = new Molecule(); //add all proximate residues; the connectivity/bonds will not be valid
        for (int i=0; i<asProx.length; i++) {
            if (asProx[i]) {
                m1.addResidue(0, m.residue[i]);
            }
        }
        saveMolecule(m1,runName+".pdb",0.0f);
    }

    //Computes the information necessary to generate a residue interaction graph for the given system;
    //		computes the minimum distance and minimum energy (absolute value) for each residue pair in residueMap[],
    //		considering all possible unpruned rotamers for the given residues;
    //		ligand interactions are also computed if a ligand is present
    private void genInteractionGraph(int numMutable, RotamerSearch rs, PrunedRotamers<Boolean> prunedRotAtRes, String runName, int strandMut[][], float eInteractionCutoff, float distCutoff, Molecule m,
                                     boolean usePairSt, float pairSt, int mutRes2Strand[],int mutRes2StrandMutIndex[]) {

        if (eInteractionCutoff<0.0f) //the cutoff should be non-negative, since we are comparing absolute values of energies against it
            eInteractionCutoff = 0.0f;

        float dist[][] = new float[numMutable][numMutable];
        float eInteraction[][] = new float[numMutable][numMutable];
        for (int i=0; i<numMutable; i++) {
            for (int j=0; j<numMutable; j++) {
                dist[i][j] = (float)Math.pow(10, 38);
                eInteraction[i][j] = 0.0f;
            }
        }

        /*float ligDist[] = null;
        float ligE[] = null;
        if (ligPresent){
        	ligDist = new float[numInAS];
        	ligE = new float[numInAS];
        	for (int i=0; i<numInAS; i++){
        		ligDist[i] = (float)Math.pow(10, 38);
        		ligE[i] = 0.0f;
        	}
        }*/

        for (int i=0; i<numMutable; i++) {
            int stri = mutRes2Strand[i];
            int strResNumi = strandMut[stri][mutRes2StrandMutIndex[i]];
            for(int q1=0; q1<rs.strandRot[stri].getNumAllowable(strResNumi); q1++) {

                int AAindex1 = rs.strandRot[stri].getIndexOfNthAllowable(strResNumi,q1);
                rs.strandRot[stri].changeResidueType(m,0,rs.strandRot[stri].rl.getAAName(AAindex1),true,true);

                int numRot1 = rs.getNumRot( stri, strResNumi, AAindex1);

                for (int r1=0; r1<numRot1; r1++) {

                    if (!prunedRotAtRes.get(i,AAindex1,r1)) { //rotamer not pruned

                        if(rs.doPerturbations)
                            ((StrandRCs)rs.strandRot[stri]).applyRC(m, strResNumi, r1);
                        else
                            rs.strandRot[stri].applyRotamer(m, strResNumi, r1);

                        for (int j=i+1; j<numMutable; j++) {
                            int strj = mutRes2Strand[j];
                            int strResNumj = strandMut[strj][mutRes2StrandMutIndex[j]];
                            for(int q2=0; q2<rs.strandRot[strj].getNumAllowable(strResNumj); q2++) {

                                int AAindex2 = rs.strandRot[strj].getIndexOfNthAllowable(strResNumj,q2);
                                rs.strandRot[strj].changeResidueType(m,strResNumj,rs.strandRot[strj].rl.getAAName(AAindex2),true,true);

                                int numRot2 = rs.getNumRot( strj, strResNumj, AAindex2 );

                                for (int r2=0; r2<numRot2; r2++) {

                                    if (!prunedRotAtRes.get(j,AAindex2,r2)) {

                                        float pairE = rs.getMinMatrix().getPairwiseE( i, AAindex1, r1, j, AAindex2, r2 );

                                        if(rs.doPerturbations)
                                            ((StrandRCs)rs.strandRot[strj]).applyRC( m, strResNumj, r2 );
                                        else
                                            rs.strandRot[strj].applyRotamer(m, strResNumj, r2);

                                        float d = m.strand[stri].residue[strResNumi].getDist(m.strand[strj].residue[strResNumj],true);
                                        dist[i][j] = Math.min(dist[i][j],d);
                                        dist[j][i] = Math.min(dist[j][i],d);

                                        if ( (!usePairSt) || (pairE<=pairSt) ) {
                                            eInteraction[i][j] = Math.max(eInteraction[i][j],Math.abs(pairE));
                                            eInteraction[j][i] = Math.max(eInteraction[j][i],Math.abs(pairE));
                                        }
                                    }
                                }
                            }
                        }

                        /*if (ligPresent) {

                        	int AAindex2 = grl.getAARotamerIndex(ligType);

                        	int numRot2 = grl.getNumRotForAAtype(AAindex2);
                        	if (numRot2==0)
                        		numRot2 = 1;

                        	for (int r2=0; r2<numRot2; r2++){

                        		if (!prunedRotAtRes[numInAS*totalNumRotamers + r2]){

                        			rs.ligROT.applyRotamer(m, 0, r2);

                        			float d = m.strand[sysStrNum].residue[residueMap[i]].getDist(m.strand[ligStrNum].residue[0],true);
                        			ligDist[i] = Math.min(ligDist[i],d);

                        			float pairE = rs.getMinMatrix()[numInAS][AAindex2][r2][i][AAindex1][r1];
                        			if ( (!usePairSt) || (pairE<=pairSt) )
                        				ligE[i] = Math.max(ligE[i],Math.abs(pairE));
                        		}
                        	}
                        }*/
                    }
                }
            }
        }

        PrintStream logPS = setupOutputFile(runName+".log");
        PrintStream logPS2 = setupOutputFile(runName);

        logPS2.println("PIG:0 "+runName); //output in Pigale-compatible ASCII format

        //Output data
        for (int i=0; i<numMutable; i++) {
            int stri = mutRes2Strand[i];
            int strResNumi = strandMut[stri][mutRes2StrandMutIndex[i]];

            int pdbResNum1 = m.strand[stri].residue[strResNumi].getResNumber();
            for (int j=i+1; j<numMutable; j++) {
                int strj = mutRes2Strand[j];
                int strResNumj = strandMut[strj][mutRes2StrandMutIndex[j]];
                int pdbResNum2 = m.strand[strj].residue[strResNumj].getResNumber();

                logPS.println(pdbResNum1+" "+pdbResNum2+" "+dist[i][j]+" "+eInteraction[i][j]);
                if ( (dist[i][j]<=distCutoff) && (eInteraction[i][j]>eInteractionCutoff) ) //these two residues interact
                    logPS2.println(pdbResNum1+" "+pdbResNum2);
            }
            /*if (ligPresent){
            	int pdbResNum2 = m.strand[ligStrNum].residue[0].getResNumber();

            	logPS.println(pdbResNum1+" "+pdbResNum2+" "+ligDist[i]+" "+ligE[i]);
            	if ( (ligDist[i]<=distCutoff) && (ligE[i]>eInteractionCutoff) )
            		logPS2.println(pdbResNum1+" "+pdbResNum2);
            }*/
        }
        logPS2.println("0 0");

        logPS.flush();
        logPS.close();
        logPS2.flush();
        logPS2.close();

        outputObject(prunedRotAtRes,runName+"_pruneInfo.obj");
    }

    //Returns a molecule m1 that contains only the residues in molecule m that are specified by residueMap[] (molecul-relative residue indexing);
    private Molecule getMolRes (Molecule m, int strandMap[][]) {

        Molecule m1 = new Molecule();

        for (int i=0; i<m.numberOfStrands; i++) { //create the same number of strands
            m1.addStrand(m.strand[i].name);
            m1.strand[i].isProtein = m.strand[i].isProtein;
        }

        for(int str=0; str<strandMap.length; str++) {
            for (int i=0; i<strandMap[str].length; i++) {

                Residue oldResidue = m.residue[strandMap[str][i]];

                Residue newResidue = new Residue();
                newResidue.name = oldResidue.name;
                newResidue.fullName = oldResidue.fullName;

                for (int j=0; j<oldResidue.numberOfAtoms; j++) {

                    Atom oldAtom = oldResidue.atom[j];

                    Atom newAtom = new Atom(oldAtom.name,oldAtom.coord[0],oldAtom.coord[1],oldAtom.coord[2]);
                    newAtom.modelAtomNumber = oldAtom.modelAtomNumber;
                    newAtom.strandNumber = oldAtom.strandNumber;
                    newAtom.elementType = oldAtom.elementType;
                    newResidue.addAtom(newAtom);
                }

                m1.addResidue(oldResidue.strandNumber,newResidue);
            }
        }

        //Determine the bonds between the atoms in the molecule
        m1.determineBonds();

        // Assign the molecule relative atom numbers
        m1.updateMoleculeAtomNumbers();

        return m1;
    }

//////////////////////////////////////////////////////
// Begin Steric Overlap Check Section
//////////////////////////////////////////////////////
    //Compute the amount of overlap between a set of structures and a reference structure
    public void handleCompStericOverlap (String s) {

        // Takes the following parameters
        // 1: System parameter filename for the reference structure (string)
        // 2: System parameter filename for the set of structures to compare (string)
        // 2: Mutation search parameter filename (string)

        // Read System parameters for the reference structure
        ParamSet sParams = new ParamSet();
        sParams.addParamsFromFile(getToken(s,2)); //read system parameters
        Molecule mRef = new Molecule();

        MolParameters mp = new MolParameters();
        loadStrandParams(sParams, mp, COMPLEX);

        mRef = setupMolSystem(mRef,sParams,mp.strandPresent,mp.strandLimits);
        int numInASref = (new Integer((String)sParams.getValue("NUMINAS"))).intValue();
        int posMapRef[] = new int[numInASref]; //molecule-relative numbering
        String posMapString = (String)sParams.getValue("RESIDUEMAP");
        for(int i=0; i<numInASref; i++)
            posMapRef[i] = (new Integer(getToken(posMapString,i+1))).intValue();



        // Read System parameters for the set of structures to compare
        sParams = new ParamSet();
        sParams.addParamsFromFile(getToken(s,3)); //read system parameters
        sParams.addParamsFromFile(getToken(s,4)); //read system parameters

        String runName = (String)sParams.getValue("RUNNAME");

        String protPDBname = (String)sParams.getValue("PROTPDBNAME");
        int numPDBfiles = (new Integer((String)sParams.getValue("NUMPDBFILES"))).intValue();
        int numInAS2 = (new Integer((String)sParams.getValue("NUMINAS"))).intValue();

        int posMap2[] = new int[numInAS2]; //molecule-relative numbering
        posMapString = (String)sParams.getValue("RESIDUEMAP");
        for(int i=0; i<numInAS2; i++)
            posMap2[i] = (new Integer(getToken(posMapString,i+1))).intValue();

        String pdbFiles[] = getPDBfiles(protPDBname,numPDBfiles); //get the PDB filenames
        numPDBfiles = pdbFiles.length;

        double minMaxOverlap[] = new double[numInAS2];
        int minMaxOverlapStruct[] = new int[numInAS2];
        for (int i=0; i<minMaxOverlap.length; i++) {
            minMaxOverlap[i] = (float)Math.pow(10, 10);
            minMaxOverlapStruct[i] = -1;
        }

        // The idea is: given a structure from the set, for each of the residues in posMap2[], determine the largest
        //		steric overlap with an atom in the posMapRef[] residues from the reference structure;
        //		then, for each residue in posMap2[], find the minimum such largest overlap among all structures in the set
        for (int i=0; i<numPDBfiles; i++) {

            System.out.println("Starting structure "+pdbFiles[i]+" ("+i+")");

            sParams.setValue("PDBNAME",pdbFiles[i]);

            //Setup the molecule system
            Molecule m2 = new Molecule();
            m2 = setupMolSystem(m2,sParams,mp.strandPresent,mp.strandLimits);

            for (int res2=0; res2<numInAS2; res2++) { //for each included residue in the given structure
                double maxOverlap = 0.0;
                for (int at2=0; at2<m2.residue[posMap2[res2]].numberOfAtoms; at2++) { //for each atom in that residue
                    Atom a2 = m2.residue[posMap2[res2]].atom[at2];
                    if ( hSteric || (!a2.elementType.equalsIgnoreCase("H"))) {
                        for (int resRef=0; resRef<numInASref; resRef++) { //for each included residue in the reference structure
                            for (int atRef=0; atRef<mRef.residue[posMapRef[resRef]].numberOfAtoms; atRef++) { //for each atom
                                Atom aRef = mRef.residue[posMapRef[resRef]].atom[atRef];
                                if ( hSteric || (!aRef.elementType.equalsIgnoreCase("H"))) {
                                    double overlap = ((a2.radius + aRef.radius)/100.0) - a2.distance(aRef);
                                    if (overlap<0.0)
                                        overlap = 0.0;
                                    maxOverlap = Math.max(maxOverlap, overlap);
                                }
                            }
                        }
                    }
                }
                if (minMaxOverlap[res2]>maxOverlap) {
                    minMaxOverlap[res2] = maxOverlap;
                    minMaxOverlapStruct[res2] = i;
                }
            }
        }


        //Output the computed distances
        PrintStream logPS = setupOutputFile(runName);
        for (int i=0; i<numInAS2; i++) {
            logPS.println(posMap2[i]+" "+minMaxOverlap[i]+" "+pdbFiles[minMaxOverlapStruct[i]]);
        }
        logPS.flush();
        logPS.close();
    }

    //Reads the pdb filenames
    private String [] getPDBfiles(String fName, int numFiles) {

        BufferedReader bufread = null;
        try {
            File file = new File(fName);
            FileReader fr = new FileReader(file);
            bufread = new BufferedReader(fr);
        } catch (FileNotFoundException e) {
            System.out.println(" ... pdbs config file not found");
            System.exit(1);
        }

        String pdbFiles[] = new String[numFiles];

        boolean done = false;
        String str = null;
        int curFile = 0;

        while (!done) {
            try {
                str = bufread.readLine();
            } catch ( Exception e ) {
                System.out.println("ERROR: An error occurred while reading input");
                System.exit(0);
            }

            if (str == null) // stop if we've reached EOF
                done = true;
            else if (getToken(str,1).charAt(0)=='%') //skip comment lines
                continue;
            else {
                if (curFile>=numFiles)
                    break;
                else {
                    pdbFiles[curFile] = getToken(str,1);
                    curFile++;
                }
            }
        }

        if (curFile<numFiles) {
            String tmp[] = new String[curFile];
            System.arraycopy(pdbFiles, 0, tmp, 0, tmp.length);
            pdbFiles = tmp;
        }

        // We're done reading them in
        try {
            bufread.close();
        } catch(Exception e) {
        }

        return pdbFiles;
    }

//////////////////////////////////////////////////////
// End Steric Overlap Check Section
//////////////////////////////////////////////////////

//////////////////////////////////////////////////////
// Begin Backrub Precomputation Section
//////////////////////////////////////////////////////
    public void handlePrecomputeBackrubs (String s) {

        // Takes the following parameters
        // 1: System parameter filename
        // 2: Number of backrub samples in each direction
        // 3: Backrub step size
        // 4: Output file name

        // Read System parameters for the reference structure
        ParamSet sParams = new ParamSet();
        sParams.addParamsFromFile(getToken(s,2)); //read system parameters
        int numBackrubSamples = new Integer(getToken(s,3)).intValue();
        float backrubStepSize = new Float(getToken(s,4)).floatValue();
        String backrubFile = getToken(s,5);


        MolParameters mp = new MolParameters();
        loadStrandParams(sParams, mp, COMPLEX);

        //Setup the molecule system
        mp.m = new Molecule();
        mp.m = setupMolSystem(mp.m,sParams,mp.strandPresent,mp.strandLimits);


        /**********Get the regions of each strand that are mutable****/
        mp.strandMut = new int[mp.strandsPresent][];  //taking the place of resMap and ligMap
        mp.strandDefault = new String[mp.strandsPresent][];
        String strandMutNums = (String)sParams.getValue("STRANDMUTNUMS");
        int ctr = 0;
        for(int i=0; i<mp.strandPresent.length; i++) {
            if(mp.strandPresent[i]) {
                int numberOfMutables = (new Integer(getToken(strandMutNums,i+1))).intValue();
                mp.strandMut[ctr] = new int[numberOfMutables];
                String strandMutResNum = (String)sParams.getValue("STRANDMUT"+i);
                for(int j=0; j<numberOfMutables; j++) {
                    String strandMutRes = getToken(strandMutResNum,j+1);
                    mp.strandMut[ctr][j] = mp.m.residue[mp.m.mapPDBresNumToMolResNum(strandMutRes)].strandResidueNumber;
                }
                mp.strandDefault[ctr] = new String[mp.strandMut[ctr].length];
                for(int j=0; j<mp.strandMut[ctr].length; j++) {
                    mp.strandDefault[ctr][j] = mp.m.strand[ctr].residue[mp.strandMut[ctr][j]].name;
                }
                ctr++;
            }
        }

        Amber96ext a96ff = new Amber96ext(mp.m, distDepDielect, dielectConst, doSolvationE, solvScale, softvdwMultiplier);
        a96ff.calculateTypesWithTemplates();

        BackrubMinimizer brMin = new BackrubMinimizer();
        brMin.initialize(mp.m, a96ff, mp.strandMut, backrubFile, hSteric, stericThresh, mp.numOfStrands, false);
        brMin.precomputeBackrubs(numBackrubSamples, backrubStepSize);
        System.out.println("DONE: Backrub angle precomputation..");
    }

    //Compute the amount of overlap between a set of structures and a reference structure
    /*public void handlePrecomputeBackrubs1 (String s) {

    	// Takes the following parameters
    	// 1: System parameter filename
    	// 2: Number of backrub samples in each direction
    	// 3: Backrub step size
    	// 4: Output file name

    	// Read System parameters for the reference structure
    	ParamSet sParams = new ParamSet();
    	sParams.addParamsFromFile(getToken(s,2)); //read system parameters
    	int numBackrubSamples = new Integer(getToken(s,3)).intValue();
    	float backrubStepSize = new Float(getToken(s,4)).floatValue();
    	String backrubFile = getToken(s,5);

    	Molecule m = new Molecule();
    	setupMolSystem(m,sParams,false,null);
    	int numInAS = (new Integer((String)sParams.getValue("NUMINAS"))).intValue();
    	int residueMap[] = new int[numInAS];
    	String resDefault[] = new String[numInAS];
    	String resMapString = (String)sParams.getValue("RESIDUEMAP");
    	System.out.print("ResidueMap:");
    	for(int i=0;i<numInAS;i++){
    		int pdbResNum = (new Integer(getToken(resMapString,i+1))).intValue();
    		residueMap[i] = m.strand[sysStrNum].mapPDBresNumToStrandResNum(pdbResNum);
    		resDefault[i] = m.strand[sysStrNum].residue[residueMap[i]].name;
    		System.out.print(" "+residueMap[i]+"("+m.strand[sysStrNum].residue[residueMap[i]].fullName+")");
    	}
    	System.out.println();

    	Amber96ext a96ff = new Amber96ext(m, distDepDielect, dielectConst, doSolvationE, solvScale, softvdwMultiplier);
    	a96ff.calculateTypesWithTemplates();

    	BackrubMinimizer brMin = new BackrubMinimizer();
    	brMin.initialize(m, a96ff, residueMap, sysStrNum, backrubFile, hSteric, stericThresh);
    	brMin.precomputeBackrubs(numBackrubSamples, backrubStepSize);

    	System.out.println("DONE: Backrub angle precomputation..");
    }*/

//////////////////////////////////////////////////////
// End Backrub Precomputation Section
//////////////////////////////////////////////////////

    int[] getCurrentMutOffset(int[][] origStrandMut, int[][] currStrandMut, boolean[] strandPresent, boolean[] oldStrandPresent) {
        int offsetMat[] = new int[origStrandMut.length];
        int strCtr=0;
        int offset=0;
        for(int str=0; str<strandPresent.length; str++) {
            if(strandPresent[str] && oldStrandPresent[str]) {
                offsetMat[strCtr] = offset;
                strCtr++;
            } else if(oldStrandPresent[str]) {
                offsetMat[strCtr] = offset;
                strCtr++;
                offset += origStrandMut[str].length;
            }
        }
        return offsetMat;
    }

    String[] shrinkCurrentMutation(String currentMut[], int origStrandMut[][] ,int newStrandMut[][],boolean origStrandPresent[], boolean newStrandPresent[]) {
        //Find the length of the new mutable array
        int newSize=0;
        for(int i=0; i<newStrandMut.length; i++)
            newSize+=newStrandMut[i].length;
        String [] newCurMut = new String[newSize];

        int NewCtr=0;
        int OldCtr=0;
        for(int str=0; str<origStrandMut.length; str++) {
            for(int i=0; i<origStrandMut[str].length; i++) {
                if(newStrandPresent[str]&&origStrandPresent[str]) {
                    newCurMut[NewCtr] = currentMut[OldCtr];
                    NewCtr++;
                    OldCtr++;
                } else if(origStrandPresent[str]) {
                    OldCtr++;
                }
            }
        }
        return newCurMut;
    }

    private int[] rotamersRemaining(int numRotForRes[], PrunedRotamers<Boolean> prunedRotAtRes) {


        //initialize numPrunedRotForRes[]
        int numNotPrunedForRes[] = new int[numRotForRes.length]; //after pruning for this partition
        for(int i=0; i<numNotPrunedForRes.length; i++) {
            numNotPrunedForRes[i] = numRotForRes[i];
        }



        Iterator<RotInfo<Boolean>> i = prunedRotAtRes.iterator();

        while(i.hasNext()) {
            RotInfo<Boolean> ri = i.next();
            if(ri.state) {
                numNotPrunedForRes[ri.curPos]--;
                if(numNotPrunedForRes[ri.curPos] == 0) {
                    System.err.println("ALL OF THE ROTAMERS HAVE BEEN PRUNED at site: "+ri.curPos);
                    Exception e = new Exception();
                    e.printStackTrace();
                    System.exit(0);
                }
            }
        }


        /*for (int j=0; j<numRotForRes.length; j++){
        	for (int k=0; k<totalNumRotamers; k++){ //pruned rot are true and must be in the current set of allowed AA
        		if (prunedRotAtRes[j*totalNumRotamers + k]){
        			numNotPrunedForRes[j]--;
        			if(numNotPrunedForRes[j] == 0){
        				System.err.println("ALL OF THE ROTAMERS HAVE BEEN PRUNED at site: "+j);
        				Exception e = new Exception();
        				e.printStackTrace();
        				System.exit(0);
        			}
        		}
        	}
        }*/

        return numNotPrunedForRes;

    }


    private void selectPerturbations(MolParameters mp, boolean doPerturbations, String pertFile, boolean minimizePerts,
                                     boolean addWTRot, ParamSet sParams ) {


        //Parameters for perturbation selection
        float min_rmsd = (new Float((String)sParams.getValue("MINRMSD","0"))).floatValue();//Minimum backbone heavy-atom RMSD for perturbation states: used to select perturbations
        Shear.setParams( ((String)sParams.getValue("SHEARPARAMS","none")) );//Default shear parameters (for perturbation selection): State 1 min, state 1 max, state 2 min, state 2 max,...
        Backrub.setParams( ((String)sParams.getValue("BACKRUBPARAMS","none")) );//Default backrub parameters (same format)
        RamachandranChecker.denCutoff = (new Float((String)sParams.getValue("RAMACHANDRANCUTOFF","0.02"))).floatValue();//Density cutoff for Ramachandran filter


        RotamerSearch rs = new RotamerSearch(mp.m,mp.numberMutable, mp.strandsPresent, hElect, hVDW, hSteric, true,
                                             true, 0.0f, stericThresh, softStericThresh, distDepDielect, dielectConst, doDihedE, doSolvationE,
                                             solvScale, softvdwMultiplier, grl, doPerturbations, pertFile, minimizePerts, false, false);

        rs.initMutRes2Str(mp.strandMut);


        boolean addWT = (new Boolean((String)sParams.getValue("ADDWT", "true"))).booleanValue();
        if(!addWT)
            checkWT(mp.strandDefault, mp.strandPresent, sParams);
        int molStrand = 0;
        for (int strNum=0; strNum<mp.numOfStrands; strNum++) {
            if(mp.strandPresent[strNum]) {
                for (int k=0; k<mp.strandMut[molStrand].length; k++) {
                    setAllowablesHelper(rs, sParams, addWT, strNum, molStrand, k, mp.strandMut, mp.strandDefault);
                }
                molStrand++;
            }
        }

        //Read the starting perturbation file, if any
        //The perturbations in this file will be included along with any the selector generates
        //and their perturbation states will be selected as appropriate
        String startingPertFile = (String)sParams.getValue("STARTINGPERTURBATIONFILE","none");//Input file giving perturbation information
        boolean onlyStarting = (new Boolean((String)sParams.getValue("ONLYSTARTINGPERTURBATIONS", "false"))).booleanValue();//Use only the perturbations specified in the startingPerturbationFile
        //(select RCs for them without selecting any more perturbations)

        if( startingPertFile.equalsIgnoreCase("none") && onlyStarting ) {
            System.err.println("ERROR: Perturbation selector can't use only starting perturbations if startingPerturbationFile is set to 'none'");
        }

        PerturbationSelector ps = new PerturbationSelector(mp.numberMutable, mp.mutRes2Strand,
                mp.mutRes2StrandMutIndex, mp.strandMut, addWTRot, mp.m,
                rs.strandRot, min_rmsd, startingPertFile,onlyStarting);

        ps.selectPerturbations();
        rs.removeImpossibleRCs(mp.numberMutable, mp.strandMut);
        PertFileHandler.writePertFile(pertFile, mp.m, null, rs.strandRot, null, false);
    }



    public void fixStruct(String s) {
        //This simple command just reads in a structure with the autofix feature on
        //and then outputs the fixed structure
        //the first argument is the PDB file to read in and the second is the filename to write to
        //(note that this is different from the usual use of configuration files)
        EnvironmentVars.autoFix = true;
        String fileIn = getToken(s,2);
        String fileOut = getToken(s,3);

        //Stop if configuration files are provided to avoid overwriting them or something
        if( fileIn.endsWith(".cfg") || fileOut.endsWith(".cfg") ) {
            System.err.println("ERROR: Just provide PDB file names for fixStruct");
            System.exit(1);
        }

        Molecule m = new Molecule();

        try {
            FileInputStream is = new FileInputStream(fileIn);
            new PDBChemModel(m, is);
        } catch (Exception e) {
            System.out.println("WARNING: An error occurred while reading file");
            System.out.println(e);
            System.exit(1);
        }

        saveMolecule(m, fileOut, Float.NaN);
    }





} // end of KSParser class
