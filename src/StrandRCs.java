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
//	StrandRCs.java
//
//	Version:           2.1 beta
//
//
//	  authors:
// 	  initials    name                 organization                email
//	 ---------   -----------------    ------------------------    ----------------------------
//	  MAH           Mark A. Hallen	  Duke University               mah43@duke.edu
///////////////////////////////////////////////////////////////////////////////////////////////

import java.util.HashMap;
import java.util.ArrayDeque;

/**
 * This class handles residue conformation (RC) assignment and maintenance for a given strand,
 * extending the capabilities of StrandRotamers to handle rotamers
 * 		Performs RC swaps and amino acid mutations for this strand.
 */
public class StrandRCs extends StrandRotamers {


    private int numRCs[][];//The number of residue conformations available to each residue in each AA type
    //(Same size, indexing as allowableAAs)

    int totNumRCs;//The total number of RCs for the strand

    private int RCOffsets[][];//Cumulative sum of numRCs, first over the residue number and then the AA type dimension,
    //indicating the first index in RotamerSearch.eliminatedRotAtRes, etc. of an RC with the given position and AA type


    //Probably need an RC count/offset per residue too

    int RCRots[][][];//Rotameric state for each RC.  -2 denotes the wild-type rotamer (only allowed for the WT residue type)
    int RCPertStates[][][];//Perturbation state for each RC
    //Indexing: RCRots[strandResidueNumber][AA type][RCIndex]

    //Perturbation states for each residue are denoted by their indices in Residue.pertStates

    private int curRC[];//Which RC each residue is in
    private int curPertState[];//Which perturbation state each residue is in

    private Residue[] WTRes;//Wild-type residues for restoration of WT rotamers
    private double[] WTGenChi1;//Generalized chi1 angles (as in Perturbation.java) for restoration of WT rotamers


    boolean isProtein;

    //Need to find a way to make numRCs and RCOffsetss
    // Generic constructor
    StrandRCs(RotamerLibrary rlP, Strand s) {
        super(rlP,s);

        numRCs = new int[numberOfResidues][rl.getNumAAallowed()];
        RCOffsets = new int[numberOfResidues][rl.getNumAAallowed()];

        curRC = new int[numberOfResidues];
        curPertState = new int[numberOfResidues];
        for(int i=0; i<numberOfResidues; i++) {
            curRC[i] = -1;
            curPertState[i] = -1;
        }

        WTRes = new Residue[numberOfResidues];
        WTGenChi1 = new double[numberOfResidues];

        isProtein = s.isProtein;
    }


    public void countRCs() { //Set up numRCs and RCOffsets
        for(int resNum=0; resNum<numberOfResidues; resNum++) {
            for(int AA=0; AA<getNumAllowable(resNum); AA++) {
                int curAA = getIndexOfNthAllowable(resNum,AA);
                numRCs[resNum][curAA] = RCRots[resNum][curAA].length;
            }

            for(int curAA=0; curAA<rl.getNumAAallowed(); curAA++) { //Get RCOffsets for all amino acids (allowed at this residue or not) to facilitate counting
                //numRCs will just be 0 for AAs not allowed at this residue
                if(curAA == 0) {
                    if( resNum > 0 )
                        RCOffsets[resNum][curAA] = RCOffsets[resNum-1][rl.getNumAAallowed()-1] + numRCs[resNum-1][rl.getNumAAallowed()-1];
                } else
                    RCOffsets[resNum][curAA] = RCOffsets[resNum][curAA-1] + numRCs[resNum][curAA-1];
            }
        }

        totNumRCs = RCOffsets[numberOfResidues-1][rl.getNumAAallowed()-1] + numRCs[numberOfResidues-1][rl.getNumAAallowed()-1];
    }



    public void addUnperturbedRCs(boolean addWTRot) {
        //Adds in RCs that are unperturbed rotamers
        //If addWTRot is true then creates an unperturbed rotamer for the WT rotameric state

        if(RCRots == null) { //RCRots and RCPertStates are expected to be both null or both not null; same with all sub-arrays within them
            //since they are formed together
            RCRots = new int[numberOfResidues][rl.getNumAAallowed()][];
            RCPertStates = new int[numberOfResidues][rl.getNumAAallowed()][];
        }
        for(int resNum=0; resNum<numberOfResidues; resNum++) { //Residue position

            if(RCRots[resNum]==null && getNumAllowable(resNum) > 0) { //We need to fill these arrays for flexible residues
                RCRots[resNum] = new int[rl.getNumAAallowed()][];
                RCPertStates[resNum] = new int[rl.getNumAAallowed()][];
            }

            for(int AA=0; AA<getNumAllowable(resNum); AA++) { //AA type
                int curAA = getIndexOfNthAllowable(resNum,AA);
                int numRot = rl.getNumRotForAAtype(curAA);
                if(numRot == 0)//Gly or Ala: Just create the single rotameric state
                    numRot = 1;
                if( addWTRot && isProtein  ) {
                    if( ( curAA == rl.getAARotamerIndex(WTRes[resNum].name) ) && ( ! WTRes[resNum].name.equalsIgnoreCase("gly") )
                            && ( ! WTRes[resNum].name.equalsIgnoreCase("pro") )  )
                        numRot++;
                    //Because glycine has no dihedrals there's no point having a wild-type rotamer for it.  Same for proline, since we don't adjust its dihedrals
                    //
                    //Alanine has only dihedrals involving the beta hydrogens but the capability is maintained for it just in case
                }

                if(RCRots[resNum][curAA] == null) {
                    RCRots[resNum][curAA] = new int[numRot];
                    RCPertStates[resNum][curAA] = new int[numRot];
                } else {
                    int oldRots[] = RCRots[resNum][curAA];
                    int oldPertStates[] = RCPertStates[resNum][curAA];
                    RCRots[resNum][curAA] = new int[oldRots.length + numRot];
                    RCPertStates[resNum][curAA] = new int[oldRots.length + numRot];
                    System.arraycopy(oldRots, 0, RCRots[resNum][curAA], numRot, oldRots.length);
                    System.arraycopy(oldPertStates, 0, RCPertStates[resNum][curAA], numRot, oldPertStates.length);
                }

                int startRot = 0;
                int shift = 0;
                if( addWTRot && isProtein ) {
                    if( ( curAA == rl.getAARotamerIndex(WTRes[resNum].name) ) && ( ! WTRes[resNum].name.equalsIgnoreCase("gly") )
                            && ( ! WTRes[resNum].name.equalsIgnoreCase("pro") )  ) {
                        RCRots[resNum][curAA][0] = -2;
                        startRot = 1;
                        shift = 1;
                    }
                }

                for(int curRot=startRot; curRot<numRot; curRot++) {
                    RCRots[resNum][curAA][curRot] = curRot-shift;
                    //The corresponding perturbation states are all 0
                    //Backbone-dependent rotamers could be easily accommodated here by checking the backbone
                    //and only creating the unperturbed rotameric state if appropriate
                }
            }
        }
    }


    public void storeWTRotamers(Molecule m) {
        //This can be called before mutations or dihedral adjustments have been performed
        //to store the WT rotameric states.  Also should be called before addUnperturbedRCs if
        //the unperturbed WT rotamer is to be added.
        //We just need to store the coordinates and residue and atom names

        m.resolveCoordinates();//Get the actualCoordinates into the Atom.coord arrays to be read

        for(int resNum=0; resNum<numberOfResidues; resNum++) { //Residue position
            if( getNumAllowable(resNum) > 0 ) {
                if( isProtein )
                    storeWTRotamer(m, resNum);
                else
                    storeWTCoord(m, resNum);
            }
        }
    }


    public void storeWTRotamer(Molecule m, int resNum) {
        //Store the WT rotamer for a single residue
        //Stored coordinates use the CA as the origin


        Residue localRes=m.strand[strandNumber].residue[resNum];
        WTRes[resNum] = new Residue(localRes.name);
        WTRes[resNum].numberOfAtoms = localRes.atom.length;
        WTRes[resNum].atom = new Atom[localRes.atom.length];



        HashMap<String, ArrayDeque<String> > HMap = new HashMap<String, ArrayDeque<String> >();
        //Maps heavy atom names to a queue of H names to use for them
        //There seem to be multiple naming conventions for sidechain hydrogens
        //So we are going to standardize them using the idealRes

        Residue idealRes = new Amber96PolyPeptideResidue().getResidue(localRes.name);
        //An ideal version of the residue, for standardizing names of sidechain hydrogens
        for(int a=0; a<idealRes.atom.length; a++) {
            Atom idealAt = idealRes.atom[a];
            if( idealAt.elementType.equalsIgnoreCase("H") ) {
                String heavyAtomName = idealRes.atom[idealAt.bond[0]].name;
                if( !HMap.containsKey(heavyAtomName) )
                    HMap.put(heavyAtomName, new ArrayDeque<String>());
                HMap.get(heavyAtomName).addLast(idealAt.name);
            }
        }

        float CACoord[] = localRes.getAtomByName("CA").coord;

        for(int a=0; a<localRes.atom.length; a++) {
            Atom localAt = localRes.atom[a];
            String WTAtName;
            if( localAt.elementType.equalsIgnoreCase("H") ) {

                int ha=0;
                Atom heavyAtom = m.atom[localAt.bond[0]];
                while( heavyAtom.moleculeResidueNumber != localRes.moleculeResidueNumber ) {
                    ha++;//A steric clash with something on another residue might have been interpreted as a bond
                    heavyAtom = m.atom[localAt.bond[ha]];//So we make sure we have the actual heavy atom
                }

                String heavyAtomName = heavyAtom.name;//Each H is bonded to one heavy atom

                if( heavyAtomName.equalsIgnoreCase("N") )//To avoid trouble with N-terminal residues
                    WTAtName = localAt.name;
                else {
                    if( HMap.get(heavyAtomName) == null ) {
                        System.err.println( "Error parsing atom names for hydrogens attached to " + heavyAtomName + " in " + localRes.fullName + " while adding WT rotamer" );
                        System.exit(1);
                    }
                    WTAtName = HMap.get(heavyAtomName).pollFirst();//Keep the ordering of hydrogen names at a given heavy atom
                    localAt.name = WTAtName;//If we restore the atom before performing a mutation we need the right name in place
                }

                if(WTAtName == null) {
                    System.err.println( "Error parsing atom names for " + localRes.fullName + " while adding WT rotamer" );
                    System.exit(1);
                }

            } else
                WTAtName = localAt.name;

            WTRes[resNum].atom[a] = new Atom( WTAtName, localAt.coord[0] - CACoord[0],
                                              localAt.coord[1] - CACoord[1], localAt.coord[2] - CACoord[2] );
        }

        WTGenChi1[resNum] = Perturbation.getGenChi1(m, localRes.moleculeResidueNumber);
    }



    public void storeWTCoord(Molecule m, int resNum) {
        //Stores the original coordinates of a residue in WTRes
        //This is used to store the untranslated and unrotated ligand coordinates
        //If the rotations and translations are implemented as a perturbation this becomes unnecessary

        Residue localRes=m.strand[strandNumber].residue[resNum];
        WTRes[resNum] = new Residue(localRes.name);
        WTRes[resNum].numberOfAtoms = localRes.atom.length;
        WTRes[resNum].atom = new Atom[localRes.atom.length];


        for(int a=0; a<localRes.atom.length; a++) {

            Atom localAt = localRes.atom[a];
            int offset = 3*localAt.moleculeAtomNumber;

            WTRes[resNum].atom[a] = new Atom( localAt.name, m.actualCoordinates[offset],
                                              m.actualCoordinates[offset+1], m.actualCoordinates[offset+2] );
        }

    }



    public void restoreWTCoord(Molecule m, int resNum) {
        //Restore the residue's coordinates from WTRes

        Residue localRes=m.strand[strandNumber].residue[resNum];

        for(int a=0; a<localRes.atom.length; a++) {

            int offset = 3*localRes.atom[a].moleculeAtomNumber;

            float WTCoord[] = WTRes[resNum].atom[a].coord;

            System.arraycopy(WTCoord, 0, m.actualCoordinates, offset, 3);
        }

    }


    // This function converts residue resNum to the residue conformation specified by RCNum
    //It will start with the actualCoordinates but will change both the molecule actualCoordinates and the atom coord arrays
    //(like applyRotamer)
    public boolean applyRC(Molecule m, int resNum, int RCNum) {


        int AANum = rl.getAARotamerIndex(m.strand[strandNumber].residue[resNum].name);
        int rotNum = RCRots[resNum][AANum][RCNum];
        int pertState = RCPertStates[resNum][AANum][RCNum];

        boolean outcome = true;

        if(rotNum == -2) { //WT rotamer

            Residue localRes=m.strand[strandNumber].residue[resNum];

            if( WTRes[resNum] == null ) {
                System.err.println("Error: WT rotamer being called but AddWTRotamers option is off");
                System.exit(1);
            }

            if( !isProtein ) {
                System.err.println("Error: WT rotamer not supported for non-protein: " + localRes.fullName);
                System.exit(1);
            }

            //MH 2012: It's an error to be assigning a WT rotamer to a mutant residue
            //But sometimes CYX gets called
            if( !localRes.name.equalsIgnoreCase( WTRes[resNum].name ) ) {
                System.err.println("Error: WT rotamer assigned to mutant residue "+localRes.fullName+", wild type " + WTRes[resNum].name );
                System.exit(1);
            }

            RotMatrix r = new RotMatrix();
            float CACoord[] = m.getActualCoord( localRes.getAtomNameToMolnum("CA") );

            //First put the atoms in in the original PDB orientation
            for(int a=0; a<localRes.numberOfAtoms; a++) {
                Atom localAt = localRes.atom[a];
                if( !localAt.isBBatom ) { //BB atoms shouldn't be restored
                    Atom WTAt = WTRes[resNum].atom[a];
                    if( ! WTAt.name.equalsIgnoreCase( localAt.name ) )//Different ordering of atoms
                        WTAt  = WTRes[resNum].getAtomByName(localAt.name);
                    localAt.coord = r.add(WTAt.coord, CACoord);
                }
            }


            for(int a=0; a<localRes.numberOfAtoms; a++) {
                if( !localRes.atom[a].isBBatom )
                    m.updateCoordinates(localRes.atom[a].moleculeAtomNumber);
            }
            //Update the actualCoordinates


            //BB conformational changes might make this orientation inappropriate though
            //So idealize the sidechain to fix it
            //This won't return false because it's not a proline
            outcome = outcome && m.idealizeResSidechain(localRes);

            //Finally make sure the (generalized) chi1 is correct
            Perturbation.setGenChi1(m, localRes.moleculeResidueNumber, WTGenChi1[resNum] );

            curRotNum[resNum] = -2;
        } else
            applyRotamer(m, resNum, rotNum);
        //applyRotamer returns false for conditions like the AA having no rotamers,
        //e.g. applying "rotamer 0" of glycine, which do not actually mean RC application failed
        //So we do not store its return value

        outcome = outcome && applyPertState(m, resNum, pertState);

        m.resolveCoordinates();//Get the perturbed coordinates into the Atom.coord arrays


        curRC[resNum] = RCNum;

        return outcome;
    }


    //Apply a residue perturbation state to the m.actualCoordinates
    //The false return value here is used to detect bad perturbation states
    public boolean applyPertState(Molecule m, int resNum, int pertState) {

        Residue localRes=m.strand[strandNumber].residue[resNum];

        boolean outcome = true;//Indicates success or failure

        if( localRes.perts.length>0 ) { //Apply perturbation state

            int pertInd = localRes.perts.length - 1;//Index in localRes.perts
            int affectedInd = localRes.affectedPerts.length - 1;//Index in localRes.affectedPerts

            //First see what perturbations we have to change
            int firstPert = m.perts.length;

            for(int a=0; a<localRes.perts.length; a++) {
                Perturbation pert = m.perts[localRes.perts[a]];
                int thisState = localRes.pertStates[pertState][a];
                if( pert.curParam != ( pert.minParams[thisState] + pert.maxParams[thisState] ) / 2 ) { //The perturbation is not at the default parameter value for the desired state
                    firstPert = localRes.perts[a];
                    break;
                }
            }

            for(int pertNum=m.perts.length-1; pertNum>=firstPert; pertNum--) { //Undoing perturbations, in reverse order
                if(pertNum == localRes.perts[pertInd]) {
                    m.perts[pertNum].undo();
                    pertInd--;
                } else if(affectedInd > -1) { //Still have affected perturbations to look for
                    if( pertNum == localRes.affectedPerts[affectedInd] ) {
                        m.perts[pertNum].undo();
                        affectedInd--;
                    }
                }
            }

            if(pertInd < localRes.perts.length - 1)
                pertInd++;//It's now the first perturbation to reapply
            if(affectedInd < localRes.affectedPerts.length - 1 )
                affectedInd++;//The first affected perturbation to reapply


            for(int pertNum=firstPert; pertNum<m.perts.length; pertNum++) { //Redoing them in the correct state
                if(pertNum==localRes.perts[pertInd]) {
                    Perturbation pert = m.perts[pertNum];

                    if(outcome) {
                        int newState = localRes.pertStates[pertState][pertInd];
                        outcome = pert.applyPerturbation( ( pert.minParams[newState] + pert.maxParams[newState] ) / 2 );
                        if(outcome)
                            pert.curState = newState;//The actual state of the perturbation (not the residue perturbation state)
                        else
                            pert.curState = 0;
                    } else {
                        pert.applyPerturbation(0);//Don't waste time applying this invalid perturbation state
                        pert.curState = 0;//Put it into an unperturbed state
                    }

                    if(pertInd < localRes.perts.length - 1)
                        pertInd++;
                } else if(affectedInd > -1) {
                    if(pertNum==localRes.affectedPerts[affectedInd]) {

                        if(outcome)
                            outcome = m.perts[pertNum].applyPerturbation(m.perts[pertNum].curParam);
                        else {
                            m.perts[pertNum].applyPerturbation(0);
                            m.perts[pertNum].curState = 0;
                        }


                        if( affectedInd < localRes.affectedPerts.length - 1 )
                            affectedInd++;
                    }
                }
            }

        }

        curPertState[resNum] = pertState;

        return outcome;
    }



    @Override
    public void changeResidueType(Molecule m, int resNum, String newResType, boolean addHydrogens, boolean connectResidue, boolean useOldBBatoms) {

        Residue localRes = m.strand[strandNumber].residue[resNum];


        //We will always mutate the unperturbed state because mutations change the backbone conformation
        //and this might not commute with perturbations
        //So we undo perturbations affecting or affected by this residue
        //This part is like applyRC
        if( localRes.perts.length>0 ) {

            int pertInd = localRes.perts.length - 1;//Index in localRes.perts
            int affectedInd = localRes.affectedPerts.length - 1;//Index in localRes.affectedPerts

            for(int pertNum=m.perts.length-1; pertNum>=localRes.perts[0]; pertNum--) { //Undoing perturbations, in reverse order
                if(pertNum == localRes.perts[pertInd]) {
                    m.perts[pertNum].undo();
                    pertInd--;
                } else if(affectedInd > -1) { //Still have affected perturbations to look for
                    if( pertNum == localRes.affectedPerts[affectedInd] ) {
                        m.perts[pertNum].undo();
                        affectedInd--;
                    }
                }
            }
        }


        //Get the actualCoordinates into the Atom.coord arrays to avoid losing perturbation information
        m.resolveCoordinates();

        //This uses Atom.coord values but alters both Atom.coord and m.actualCoordinates
        super.changeResidueType(m, resNum, newResType, addHydrogens, connectResidue, useOldBBatoms);

        curRC[resNum] = -1;
        //The Perturbation.curState is left the same for each perturbation affecting the residue
        //ans so the residue's perturbation state is the same as before

        //Reapplying perturbations
        if( localRes.perts.length>0 ) {

            int pertInd = 0;
            int affectedInd = 0;

            for(int pertNum=0; pertNum<m.perts.length; pertNum++) {
                if(pertNum==localRes.perts[pertInd]) {
                    m.perts[pertNum].applyPerturbation(m.perts[pertNum].curParam);
                    if(pertInd < localRes.perts.length - 1)
                        pertInd++;
                } else if(localRes.affectedPerts.length > 0) {
                    if(pertNum==localRes.affectedPerts[affectedInd]) {
                        m.perts[pertNum].applyPerturbation(m.perts[pertNum].curParam);
                        if(affectedInd < localRes.affectedPerts.length - 1 )
                            affectedInd++;
                    }
                }
            }
        }

    }



    public void removeRCs(boolean[][][] toRemove) { //Remove the RCs with true values in toRemove (indices: position in strand, AA type, RC number)


        for(int pos=0; pos<numberOfResidues; pos++) {

            for(int AA=0; AA<rl.getNumAAallowed(); AA++) {

                int RCCount=0;

                for(int RC=0; RC<getNumRCs(pos,AA); RC++) {

                    if( ! toRemove[pos][AA][RC] )
                        RCCount++;
                }

                int newRots[] = new int[RCCount];
                int newPertStates[] = new int[RCCount];

                int place=0;

                for(int RC=0; RC<getNumRCs(pos,AA); RC++) {

                    if( ! toRemove[pos][AA][RC] ) {
                        newRots[place] = RCRots[pos][AA][RC];
                        newPertStates[place] = RCPertStates[pos][AA][RC];

                        place++;
                    }

                }

                RCRots[pos][AA] = newRots;
                RCPertStates[pos][AA] = newPertStates;
            }
        }

        countRCs();//The numbering has changed so recount the RCs
    }


    public int getNumRCs(int resNum, int AANum) {
        return numRCs[resNum][AANum];
    }

    public int getTotNumRCs(int resNum) { //Total number of RCs at a given residue position
        int count = 0;

        for(int curAA=0; curAA<rl.getNumAAallowed(); curAA++)
            count+=numRCs[resNum][curAA];

        return count;
    }


    public int getTotNumPerturbedRCs(int resNum) { //Total number of perturbed RCs at a given residue position
        int count = 0;

        for(int curAA=0; curAA<rl.getNumAAallowed(); curAA++) {
            for(int RC=0; RC<numRCs[resNum][curAA]; RC++) {
                if( RCPertStates[resNum][curAA][RC] != 0 )
                    count++;
            }
        }

        return count;
    }


    public int getRCOffset(int resNum, int AANum) {
        return RCOffsets[resNum][AANum];
    }

    public int getCurPertState(int resNum) {
        return curPertState[resNum];
    }

    public void setRCPertStates(int[][][] RCPertStates) {
        this.RCPertStates = RCPertStates;
    }

    public void setRCRots(int[][][] RCRots) {
        this.RCRots = RCRots;
    }

    public int getNumAATypes() {
        return rl.getNumAAallowed();
    }

}
