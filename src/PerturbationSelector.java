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
//	PerturbationSelector.java
//
//	Version:           2.1 beta
//
//
//	  authors:
// 	  initials    name                 organization                email
//	 ---------   -----------------    ------------------------    ----------------------------
//	  MAH           Mark A. Hallen	  Duke University               mah43@duke.edu
///////////////////////////////////////////////////////////////////////////////////////////////
import java.util.ArrayList;
import java.io.*;

public class PerturbationSelector {
    //This class selects perturbations to apply to a particular molecule, and RCs to allow with those perturbations

    boolean addWTRot;//Should the wild-type "rotamers" be considered?

    int flexMolResNum[];//Molecule residue numbers of flexible residues
    int numMutable;//Total number of flexible residues


    StrandRCs sRC[];
    Molecule m;
    float min_rmsd;

    String startingPertFile;
    boolean onlyStarting = false;


    //Now some parameters

    public PerturbationSelector(int numMut, int mutRes2Strand[], int mutRes2StrandMutIndex[], 
            int strandMut[][], boolean addWTRotamers, Molecule molec,
            StrandRotamers[] strandRot, float minrmsd, String startPF, boolean onlyStart){

        addWTRot = addWTRotamers;

        numMutable = numMut;
        flexMolResNum = new int[numMutable];

        for(int pos=0; pos<numMutable; pos++){
            int str = mutRes2Strand[pos];
            int strResNum = strandMut[str][mutRes2StrandMutIndex[pos]];
            flexMolResNum[pos] = molec.strand[str].residue[strResNum].moleculeResidueNumber;
        }


        m=molec;
        sRC=new StrandRCs[m.numberOfStrands];
        for(int str=0; str<m.numberOfStrands; str++)
            sRC[str] = (StrandRCs)strandRot[str];

        min_rmsd = minrmsd;
        startingPertFile = startPF;
        onlyStarting = onlyStart;
    }


    //This function generates the set of possible perturbations for the molecule
    //and their states
    public void generatePerturbations(Molecule m){

        final int numPertTypes = 8;

        Perturbation[][] pertsByType = new Perturbation[numPertTypes][];

        pertsByType[0] = new Perturbation[0];
        //Not currently auto-generating partial or full structure switches
        //except in starting perturbation file (if present)

        //The starting perturbations will probably be partial structure switches
        //so they will be put before the other perturbations as the 2nd perturbation type
        if( startingPertFile.equalsIgnoreCase("none") )
            pertsByType[1] = new Perturbation[0];
        else{

            try{

                BufferedReader br = new BufferedReader( new FileReader(startingPertFile) );
                br.readLine();//Skip title line

                pertsByType[1] = PertFileHandler.readPerts(br, m);
                br.close();
            }
            catch(Exception e){
                System.err.println(e.getMessage());
                System.exit(1);
            }

        }

        int pertCount = 0;

        if(onlyStarting){
            m.perts = pertsByType[1];
            pertCount = m.perts.length;
        }
        else {//generate other perturbations unless only using those in the starting perturbation file
            pertsByType[2] = Perturbation.generateAll("SSNE", m, flexMolResNum, 4);
            pertsByType[3] = Perturbation.generateAll("SSCE", m, flexMolResNum, 4);
            pertsByType[4] = Perturbation.generateAll("LOOP CLOSURE ADJUSTMENT", m, flexMolResNum, 3);
            pertsByType[5] = Perturbation.generateAll("SHEAR", m, flexMolResNum, 4);
            pertsByType[6] = Perturbation.generateAll("BACKRUB", m, flexMolResNum, 3);
            pertsByType[7] = ProlineFlip.generateAll(m, sRC);

            //Copy them all into m.perts
            for(int a=0;a<numPertTypes;a++)
                pertCount += pertsByType[a].length;

            m.perts = new Perturbation[pertCount];
            int start = 0;
            for(int a=0;a<numPertTypes;a++){
                System.arraycopy(pertsByType[a], 0, m.perts, start, pertsByType[a].length);
                start += pertsByType[a].length;
            }
        }

        System.out.println( pertCount + " perturbations to be used");

    }


    
    //Select the perturbations to be used from the set of possible perturbations
    //and load them into the molecule m, its residues, and the RC handlers
    public void selectPerturbations(){

        //First define all the feasible perturbations and their states and successors
        generatePerturbations(m);
        m.findPerturbationSuccessors();

      
        for(int str=0; str<m.numberOfStrands; str++){
            if(addWTRot)
                sRC[str].storeWTRotamers(m);
            sRC[str].addUnperturbedRCs(addWTRot);
        }


        RamachandranChecker rcheck = RamachandranChecker.getInstance();

        boolean proFlip[] = new boolean[numMutable];//Indicates with residues have a proline flip

        for(int curPos=0; curPos<numMutable; curPos++){

            //First find the perturbations affecting each residue
            Residue res = getResidue(curPos);

            res.pertStates = new int[1][];
            
            ArrayList<Integer> pertList = new ArrayList<Integer>();
            //Perturbations affecting res

            for( int pertNum=0; pertNum<m.perts.length; pertNum++){
                Perturbation pert = m.perts[pertNum];
                boolean affecting=false;

                for( int b=0 ; b<pert.resAffected.length; b++){
                    if( pert.resAffected[b] == res.moleculeResidueNumber )
                        affecting = true;
                }

                for( int p2 : pert.successors ){
                    Perturbation pert2 = m.perts[p2];
                    for( int b=0 ; b<pert2.resAffected.length; b++){
                        if( pert2.resAffected[b] == res.moleculeResidueNumber )
                            affecting = true;
                    }
                }

                if(affecting)
                    pertList.add(pertNum);
            }

            int resNumPerts = pertList.size();

            if(resNumPerts > 0){
            
                res.perts = new int[resNumPerts];
                for(int b=0; b<resNumPerts; b++)
                    res.perts[b] = pertList.get(b).intValue();

                //Also assign the affectedPerts
                m.assignAffectedPerts(res);

                proFlip[curPos] = m.perts[res.perts[res.perts.length-1]].type.equalsIgnoreCase("PROLINE FLIP");


                //Now find their valid and Ramachandran-allowed states
                //We let the unperturbed states be allowed by default
                res.pertStates = new int[1][];

                int state[] = new int[resNumPerts];//Start with the unperturbed state (denoted by all 0's)
                
                boolean isUnperturbed = true;//The first state is unperturbed (all 0's for perturbation states)
                ArrayList<int[]> pertStateList = new ArrayList<int[]>();//Residue perturbation states


                boolean done = false;

                int badState[] = {-1};
                //If a residue perturbation state whose last perturbations applied are in the unperturbed state
                //is invalid, then applying those last perturbations will not make it valid
                //So if we encounter such a state, we record it up to the last perturbation applied with nonzero parameter,
                //call this partial state badState, and then don't consider any perturbation states that start with badState

                while( ! done ){

                    boolean valid = false;

                    for(int a=0; a<badState.length; a++){
                        if( state[a] != badState[a] ){
                            valid = true;
                            break;
                        }
                    }


                    if(isUnperturbed){
                        pertStateList.add(state.clone());

                        //If there's a proline flip the second state will be unperturbed except for that flip; we make sure to keep it too
                        if( ! ( proFlip[curPos] && pertStateList.size() == 1 ) )
                            isUnperturbed = false;//Otherwise there are no more unperturbed state

                    }
                    else if(valid){

                        //Try to apply the perturbation state
                        res.pertStates[0] = state;

                        boolean ramaCheck = true;

                        valid = sRC[res.strandNumber].applyPertState(m, res.strandResidueNumber, 0);

                        if(valid){
                            boolean rama[] = rcheck.checkByAAType(m, res.moleculeResidueNumber);

                            if( ! ( rama[0] || rama[1] || rama[2] ) )//Bad state regardless of AA type
                                ramaCheck = false;
                            //rama is not stored in the interest of saving memory...we have to reapply the perturbation states later again anyway
                            //if there is an RMSD filter
                        }


                        if(valid && ramaCheck)
                            pertStateList.add(state.clone());
                        
                        else if (!valid){

                            int lastPert = resNumPerts-1;
                            while( state[lastPert]==0 )
                                lastPert--;

                            badState = new int[lastPert+1];

                            for(int a=0; a<=lastPert; a++)
                                badState[a] = state[a];
                        }

                    }

                    //Now find the next state
                    //We go through the states in order of the last perturbation, then the second to last, etc.
                    for(int pertNum=resNumPerts-1; pertNum>=0; pertNum--){

                        state[pertNum]++;

                        if( state[pertNum] < m.perts[res.perts[pertNum]].minParams.length )
                            break;
                        else if( pertNum > 0 )
                            state[pertNum] = 0;

                    }

                    if( state[0] >= m.perts[res.perts[0]].minParams.length )
                        done = true;
                }

                //Store the perturbation states
                res.pertStates = new int[pertStateList.size()][resNumPerts];

                for( int stateNum=0; stateNum<pertStateList.size(); stateNum++)
                    res.pertStates[stateNum] = pertStateList.get(stateNum);

                System.out.println( pertStateList.size() + " perturbation states considered for flexible residue " + curPos );
            }
            else//No perturbations at this residue
                res.pertStates = new int[1][0];//Unperturbed state

        }

        //Remove states that are impossible because they're incompatible with everything at some other residue
        boolean[][] impossiblePertStates = getImpossiblePertStates();

        for(int curPos=0; curPos<numMutable; curPos++)
            removePertStates(getResidue(curPos), countTrue(impossiblePertStates[curPos]), impossiblePertStates[curPos], null);

        impossiblePertStates = null;//Delete this because it might be big


        boolean ramaAllowed[][][] = new boolean[numMutable][][];

        //Now apply the RMSD filter if desired (get rid of perturbation states whose
        //conformations are all within backbone heavy-atom RMSD min_rmsd of another's)
        //Also calculate Ramachandran acceptability for residue types
        for(int curPos=0; curPos<numMutable; curPos++){

            Residue res = getResidue(curPos);

            int resNumStates = res.pertStates.length;

            boolean prunedPertStates[] = new boolean[resNumStates];
            boolean ramaAllowedLong[][] = new boolean[resNumStates][3];//Indicates whether each perturbation state is Ramachandran-allowed for gly, pro, or other residues
            int numPruned = 0;

            int startState = 0;
            if( proFlip[curPos] )
                startState = 1;


            for(int state1=startState; state1<resNumStates; state1++){


                if( !prunedPertStates[state1] ){//Apply perturbation state and do Ramachandran check

                    sRC[res.strandNumber].applyPertState(m, res.strandResidueNumber, state1);
                    ramaAllowedLong[state1] = rcheck.checkByAAType(m, res.moleculeResidueNumber);

                    if( min_rmsd > 0 ){//Perturbation state valid and not removed by Ramachandran filter, so apply RMSD filter

                        m.backupAtomCoord();

                        for(int state2=state1+1; state2<resNumStates; state2++){

                            if( !prunedPertStates[state2] ){

                                sRC[res.strandNumber].applyPertState(m, res.strandResidueNumber, state2);

                                float rmsd = 0;
                                String heavyAtoms[] = {"N","CA","C","O"};

                                Perturbation lastPert = m.perts[res.perts[res.perts.length - 1]];
                                for(int q=0; q<lastPert.resAffected.length; q++){
                                    //We take the RMSD over all heavy atoms in residues affected by the last perturbation applied,
                                    //because removing the perturbation state will remove the relevant flexibility from all these residues
                                    Residue res2 = m.residue[lastPert.resAffected[q]];
                                    for(int a=0;a<4;a++){
                                        int molAtNum = res2.getAtomNameToMolnum(heavyAtoms[a]);
                                        for(int b=0;b<3;b++){
                                            float dev = m.actualCoordinates[3*molAtNum + b] - m.backupCoordinates[3*molAtNum + b];
                                            rmsd += dev*dev;
                                        }
                                    }
                                }

                                rmsd = (float)Math.sqrt(rmsd/(4*lastPert.resAffected.length));

                                if(rmsd < min_rmsd){
                                    prunedPertStates[state2] = true;
                                    numPruned++;
                                }
                            }
                        }
                    }
                }
            }

            if( min_rmsd > 0 )
                System.out.println( "Flexible residue " + curPos + ": " + numPruned + " perturbation states pruned by RMSD filter" );

            ramaAllowed[curPos] = removePertStates(res, numPruned, prunedPertStates, ramaAllowedLong);

        }

        //Remove states that have newly become impossible
        impossiblePertStates = getImpossiblePertStates();

        for(int curPos=0; curPos<numMutable; curPos++)
            ramaAllowed[curPos] = removePertStates( getResidue(curPos), countTrue(impossiblePertStates[curPos]), impossiblePertStates[curPos], ramaAllowed[curPos]);

        impossiblePertStates = null;


        //Now generate all the perturbed RCs: all combinations of unpruned perturbation states
        //and rotamers (i.e. the unperturbed RCs already in place)
        //Add them to the unperturbed RCs in the RC handlers
        for(int curPos=0; curPos<numMutable; curPos++){

            Residue res = getResidue(curPos);
            int resNumStates = res.pertStates.length;


            if(resNumStates > 1){//There is at least one perturbed state
                
                int[][] oldRCRots = sRC[res.strandNumber].RCRots[res.strandResidueNumber];
                //These are the unperturbed RCs

                int newRCRots[][] = new int[oldRCRots.length][];
                int newRCPertStates[][] = new int[oldRCRots.length][];

                for(int AA=0;AA<oldRCRots.length;AA++){

                    if( oldRCRots[AA] != null ){//If this AA type is allowed at this position

                        boolean goodState[] = new boolean[resNumStates];//Indicates which perturbation states pass the Ramachandran filter for this residue type
                        int resGoodStates = 0;

                        for( int stateNum=0; stateNum<resNumStates; stateNum++ ){

                            boolean flippedPro = false;//This perturbation state is a flipped-proline state
                            if( m.perts[res.perts[res.perts.length-1]].type.equalsIgnoreCase("PROLINE FLIP") ){
                                if( res.pertStates[stateNum][res.perts.length-1] == 1 )//A proline flip would be the last perturbation
                                    flippedPro = true;
                            }

                            goodState[stateNum] = true;

                            String AAType = sRC[res.strandNumber].rl.getAAName(AA);


                            if( /*( curPos < numSiteResidues ) &&*/ ( stateNum > 0 ) &&
                                    ( (!proFlip[curPos]) || (!AAType.equalsIgnoreCase("pro")) ||  (stateNum>1)  ) ){
                                //Active-site-residue perturbed RCs are Ramachandran filtered on a residue-type basis
                                //We don't filter the unperturbed or any proline-flipped unperturbed states

                                if( AAType.equalsIgnoreCase("gly") ){
                                    if( ( ! ramaAllowed[curPos][stateNum][0] ) || flippedPro )
                                        goodState[stateNum] = false;
                                }
                                else if( AAType.equalsIgnoreCase("pro") ){
                                    if( ! ramaAllowed[curPos][stateNum][1] )
                                        goodState[stateNum] = false;
                                }
                                else if ( ( ! ramaAllowed[curPos][stateNum][2] ) || flippedPro )
                                    goodState[stateNum] = false;
                            }

                            if(goodState[stateNum])
                                resGoodStates++;
                        }

                        newRCRots[AA] = new int[ oldRCRots[AA].length * resGoodStates ];
                        newRCPertStates[AA] = new int[ oldRCRots[AA].length * resGoodStates ];
                        for( int rotNum=0; rotNum<oldRCRots[AA].length; rotNum++ ){
                            int b=0;
                            for( int stateNum=0; stateNum<resNumStates; stateNum++){
                                if( goodState[stateNum] ){
                                    newRCRots[AA][b*oldRCRots[AA].length + rotNum] = oldRCRots[AA][rotNum];//This ordering of RCs puts unperturbed ones first
                                    newRCPertStates[AA][b*oldRCRots[AA].length + rotNum] = stateNum;
                                    b++;
                                }
                            }
                        }
                    }
                }

                sRC[res.strandNumber].RCRots[res.strandResidueNumber] = newRCRots;
                sRC[res.strandNumber].RCPertStates[res.strandResidueNumber] = newRCPertStates;
            }
        }


        for(int str=0; str<m.numberOfStrands; str++)
            sRC[str].countRCs();
    }


    


    private boolean[][] getImpossiblePertStates(){

        boolean done = false;
        boolean prunedStates[][] = new boolean[numMutable][];//Which RCs are pruned


        while( !done ){//We iterate until no more perturbation states can be removed

            done = true;

            for (int curPos=0; curPos<numMutable; curPos++){

                    Residue res = getResidue(curPos);

                    int resNumStates = res.pertStates.length;

                    if(prunedStates[curPos] == null)//This will happen in the first iteration
                        prunedStates[curPos] = new boolean[resNumStates];

                    for(int curState=0; curState<resNumStates; curState++){

                        if( ! prunedStates[curPos][curState] ){

                            for (int altPos=0; altPos<numMutable; altPos++){

                                if( altPos != curPos ){

                                    boolean prune = true;//Indicates no perturbation state compatible with curRC has been found yet at altPos
                                    Residue res2 = getResidue(altPos);

                                    for(int altState=0; altState<res2.pertStates.length; altState++){

                                        if( ! arePertStatesIncompatible(res, curState, res2, altState) )
                                            prune = false;
                                    }

                                    prunedStates[curPos][curState] = prunedStates[curPos][curState] || prune;
                                }
                            }

                            if( prunedStates[curPos][curState] )//Iterate again if anything was pruned
                                done = false;

                        }
                    }

            }

        }

        return prunedStates;

    }




    private boolean[][] removePertStates(Residue res, int numPruned, boolean prunedPertStates[], boolean ramaAllowedLong[][]){
        //Removes the perturbation states indicated by true values in prunedPertStates from res
        //There are numPruned pruned states
        //We also remove the corresponding arrays from ramaAllowedLong and return the answer if ramaAllowedLong != null

        boolean ramaAllowed[][] = null;

        int resNumStates = res.pertStates.length;
        int resNumPerts = res.perts.length;

        //Now get rid of the unnecessary perturbation states
        resNumStates -= numPruned;
        int newPertStates[][] = new int[resNumStates][resNumPerts];

        if(ramaAllowedLong != null)
            ramaAllowed = new boolean[resNumStates][3];
        
        int newState = 0;
        for(int oldState=0; oldState<res.pertStates.length; oldState++){

            if( ! prunedPertStates[oldState] ){
                
                System.arraycopy( res.pertStates[oldState], 0, newPertStates[newState], 0, resNumPerts );
                if(ramaAllowedLong != null)
                    System.arraycopy( ramaAllowedLong[oldState], 0, ramaAllowed[newState], 0, 3);
                newState++;
            }
        }

        res.pertStates = newPertStates;

        return ramaAllowed;
    }



    private Residue getResidue(int curPos){//Get the residue at a given flexible position
        return m.residue[flexMolResNum[curPos]];
    }



    //Check parametric incompatibility of a pair of perturbation states
    public boolean arePertStatesIncompatible(Residue firstRes, int pertState1, Residue secondRes, int pertState2 ){

        //This loop is the same as in RotamerSearch.isParametricallyIncompatible
        for(int p1=0;p1<firstRes.perts.length;p1++){//Check all the perturbations of firstRes in this perturbation state for compatibility
            for(int p2=0;p2<secondRes.perts.length;p2++){
                if(firstRes.perts[p1] == secondRes.perts[p2]){//The residues share a perturbation
                    if(firstRes.pertStates[pertState1][p1] != secondRes.pertStates[pertState2][p2])//The states of this perturbation don't match for the given RCs
                        return true;//So the perturbation states are parametrically incompatible
                }
            }
        }

        return false;//No incompatibility found, so we're good
    }


    private int countTrue(boolean[] x){//Count how many values in the boolean array are true

        int count=0;

        for(int a=0; a<x.length;a++){
            if(x[a])
                count++;
        }

        return count;
    }
        
   
        /* Might steric-check perturbation states (BB only) to reduce number
         * Consider looping through positions, AA types, rotamers, perturbations states and checking for feasible RCs:
         * boolean goodRCs[][][][];//Indices: position (in active site), amino-acid type, perturbation state, rotamer
         * However, if this is just a steric check that can be done during energy precomputations
         *
         * If needed, a more stringent sorting system can
         * 1. Use trainable parameters to score amino acid-perturbation combinations according to known mutant-WT structure pairs
         * 2. Look at energies (quick contact score or actual rigid energy), and/or
         * 3. Incorporate these into a score that is calculated for each RC (e.g. in float RCscores[][][][])
         * RCs with the same perturbation states can be partially averaged (perturbation states can be scored
         * too) and then the top n RCs can be taken where n is specified by the user
         * The issue with energetics is it'll get rid of conformations already well pruned by DEE
         * RMSD checks are complementary to DEE
         * For example score RCs based on energies
         * Then weighted-sum with perturbation scores
         * Generate a perturbation file
         */

}
