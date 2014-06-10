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
//	LoopClosureAdjustment.java
//
//	Version:           2.1 beta
//
//
//	  authors:
// 	  initials    name                 organization                email
//	 ---------   -----------------    ------------------------    ----------------------------
//	  MAH           Mark A. Hallen	  Duke University               mah43@duke.edu
///////////////////////////////////////////////////////////////////////////////////////////////

/* This perturbation involves finding alternate ways to close loops
 * A plain "loop closure adjustment" does this by finding alternate sets of dihedrals
 * for a three-residue segment without altering the chain elsewhere
 * This is a discrete type of perturbation

 * There are also some specialized versions:
 * 1. A helix or sheet reduction is a loop closure adjustment containing one or more non-loop residues
 * (which will be removed from their helix or sheet)
 * 2. A helix or sheet expansion is formed by adding a loop residue to an adjacent helix or sheet
 * and then trying to close the loop
 * 3. A partial structure switch uses experimental data (maybe its own class??)
 *
 */
import java.util.HashMap;
import java.util.ArrayList;


public class LoopClosureAdjustment extends Perturbation {

    // HashMap<float[],float[][][][]> confMap;
    //Key: Perturbation parameters for the predecessors.  Value: Rotation matrices to use;


    float predParams[][] = new float[0][];//Perturbation parameters for predecessors.
    //Indices: 1. Which set of predecessor parameters this is (state index) 2. Which perturbation

    float[][][][][] rotMatrices = new float[0][][][][];//Rotation matrices
    //indices: 1. State index (like predParams) 2. # of solution 3. number of matrix
    //(0 for the middle sidechain & CA & everything before it; 1 for everything after) 4-5. Matrix indices

    //Should the first and/or last residue copy the dihedrals (i.e. the secondary structure) of the one
    //before or after it respectively?
    boolean addSSStart = false;
    boolean addSSEnd = false;
    //Possibly also add the option of some out-of-plane motion before and/or after?
    //(Added the same way)

//Partial structure switch: here or separate?  Not sure

    public LoopClosureAdjustment(String t, Molecule molec, int resList[]) {

        type=t;
        m=molec;
        resAffected=resList;

        if(type.equalsIgnoreCase("SSNE"))//Secondary structure N-terminal expansion
            addSSEnd = true;
        else if(type.equalsIgnoreCase("SSCE"))//Secondary structure C-terminal expansion
            addSSStart = true;

    }

    public void setDefaultParams() { //There can be up to 16 solutions so we create that many states
        //They're all discrete though
        minParams = new float[16];
        maxParams = new float[16];
        for(int a=0; a<16; a++)
            minParams[a] = maxParams[a] = a;
    }



    public boolean doPerturbationMotion(float param) {

        int predState = -1;//If the predecessors have been in this state before this will indicate which set of known parameters this is

        for( int a=0; ( a<predParams.length ) && ( predState == -1 ); a++ ) { //Look for a match in predParams to the current predecessor parameters

            predState = a;

            for( int b=0; b<predecessors.length; b++ ) {
                if( predParams[predState][b] != m.perts[predecessors[b]].curParam ) { //Predecessor parameters
                    predState = -1;
                    break;
                }
            }

        }

        if(predState == -1) { //We need to calculate the solutions
            calcSolns();
            predState = predParams.length - 1;
        }


        if( (param < 0) || param >= ( rotMatrices[predState].length) )
            return false;//The parameter value is invalid for this state of the predecessors
        //So we fail to apply the perturbation

        applySSChanges();

        float[][][] rm = rotMatrices[predState][(int)param];

        if(addSSStart) {
            int tripepRes[] = {resAffected[1], resAffected[2], resAffected[3]};
            applyTripeptideRot(rm, tripepRes);
        } else
            applyTripeptideRot(rm, resAffected);

        return true;
    }





    public void applyTripeptideRot(float rm[][][], int res[]) {
        //Perform the loop closure adjustment rotations on three consecutive residues,
        //given the rotation matrices and the molecule residue numbers


        //All the rotation are about anchor-point CAs, so order is not important here like for the shear
        //The second peptide plane is rotated about the ending anchor CA; the other stuff is rotated about the starting one
        int CA0AtNum = m.residue[res[0]].getAtomNameToMolnum("CA");
        int CA2AtNum = m.residue[res[2]].getAtomNameToMolnum("CA");


        //We rotate the middle residue except the carbonyl, plus the carbonyl group of the first, with the 1st matrix
        int pep1sc[] = concatenateArrays( m.residue[res[0]].getAtomList(false, false, false, true),
                                          m.residue[res[1]].getAtomList(true, true, true, false) );
        m.rotateAtomList( pep1sc, rm[0], m.actualCoordinates[3*CA0AtNum],
                          m.actualCoordinates[3*CA0AtNum+1], m.actualCoordinates[3*CA0AtNum+2], false);


        //Rotate the amide group of the last residue and carbonyl of the middle residue
        //with the second matrix
        int pep2[] = concatenateArrays( m.residue[res[1]].getAtomList(false, false, false, true),
                                        m.residue[res[2]].getAtomList(true,false,false,false) );
        m.rotateAtomList( pep2, rm[1], m.actualCoordinates[3*CA2AtNum],
                          m.actualCoordinates[3*CA2AtNum+1], m.actualCoordinates[3*CA2AtNum+2], false);

    }






    public int calcSolns() { //Solve the loop closure equations given the current state of predecessor perturbations
        //Store conformational change info and return the number of solutions

        //Get necessary lengths, angles, dihedrals

        RotMatrix rm = new RotMatrix();
        float[][][][] matrices = null;


        double[] genChi1 = null;
        int tripepRes[] = null;//Residues in the tripeptide to be closed

        //Backup coordinates for any residues to be moved (for secondary-structure adjustment purposes)
        //Would need to do the other same for other motions outside the closed tripeptide
        //We use oldCoords, like in undo; since we have not applied the perturbation yet this does not interfere with undoing
        if(addSSStart || addSSEnd) {
            oldCoords = new ArrayList<HashMap<String, float[]>>(resAffected.length);
            genChi1 = getGenChi1();//Note the generalized chi1 for each residue so we can restore it after the changes

            for(int a=0; a<resAffected.length; a++)
                oldCoords.add(a, storeResBB(m.residue[resAffected[a]]) );

            tripepRes = new int[3];

            if(addSSStart)
                System.arraycopy(resAffected, 1, tripepRes, 0, 3);
            else
                System.arraycopy(resAffected, 0, tripepRes, 0, 3);
        } else
            tripepRes = resAffected;


        TripeptideClosure tc = new TripeptideClosure(m, tripepRes);
        //This grabs the bond lengths, angles, and omega dihedrals before we start changing them

        applySSChanges();


        float r_soln_n[][][] = new float[16][3][3];
        float r_soln_a[][][] = new float[16][3][3];
        float r_soln_c[][][] = new float[16][3][3];

        float firstN[] = m.getActualCoord( m.residue[ tripepRes[0] ].getAtomNameToMolnum("N") );//Coordinates of the first residue's N
        float firstCA[] = m.getActualCoord( m.residue[ tripepRes[0] ].getAtomNameToMolnum("CA") );//Its CA
        float firstC[] = m.getActualCoord( m.residue[ tripepRes[0] ].getAtomNameToMolnum("C") );

        float midN[] = m.getActualCoord( m.residue[ tripepRes[1] ].getAtomNameToMolnum("N") );//Starting coordinates of the middle N
        float midCA[] = m.getActualCoord( m.residue[ tripepRes[1] ].getAtomNameToMolnum("CA") );
        float midC[] = m.getActualCoord( m.residue[ tripepRes[1] ].getAtomNameToMolnum("C") );

        float lastN[] = m.getActualCoord( m.residue[ tripepRes[2] ].getAtomNameToMolnum("N") );
        float lastCA[] = m.getActualCoord( m.residue[ tripepRes[2] ].getAtomNameToMolnum("CA") );
        float lastC[] = m.getActualCoord( m.residue[ tripepRes[2] ].getAtomNameToMolnum("C") );


        int numSoln = tc.solve_3pep_poly( firstN, firstCA, lastCA, lastC, r_soln_n, r_soln_a, r_soln_c);


        int shift = 0;
        if( addSSStart || addSSEnd )
            shift = 1;
        //This is to make room for the unperturbed state at matrices[0]

        matrices = new float[numSoln+shift][2][][];

        int unperturbed = -1;//The unperturbed state: will be found by least-squares comparison to the original state
        float lowestSum = Float.POSITIVE_INFINITY;


        for(int s=0; s<numSoln; s++) {

            //First rotation: Based on midCA and midN
            matrices[s+shift][0] = rm.getSuperposingRotMatrix( rm.subtract( midCA, firstCA ), rm.subtract( r_soln_a[s][1], firstCA),
                                   rm.subtract( midN, firstCA ), rm.subtract( r_soln_n[s][1], firstCA) );

            //This rotation is about the last CA instead of the first one
            matrices[s+shift][1] = rm.getSuperposingRotMatrix( rm.subtract( midC, lastCA ), rm.subtract( r_soln_c[s][1], lastCA),
                                   rm.subtract( lastN, lastCA ), rm.subtract( r_soln_n[s][2], lastCA) );


            if( ! (addSSStart || addSSEnd) ) { //See if this might be the unperturbed state
                float checkSum = rm.normsq( rm.subtract(firstC, r_soln_c[s][0]) )
                                 + rm.normsq( rm.subtract(midN, r_soln_n[s][1]) )
                                 + rm.normsq( rm.subtract(midCA, r_soln_a[s][1]) )
                                 + rm.normsq( rm.subtract(midC, r_soln_c[s][1]) )
                                 + rm.normsq( rm.subtract(lastN, r_soln_n[s][2]) );

                if(checkSum < lowestSum) {
                    lowestSum = checkSum;
                    unperturbed = s;
                }
            }

        }

        if(addSSStart || addSSEnd) { //This is like undo.  In this case we do not expect the unperturbed state to have been generated
            for(int a=0; a<resAffected.length; a++)
                restoreResBB( m.residue[resAffected[a]], oldCoords.get(a) );

            setGenChi1(genChi1);
        } else {

            if(numSoln == 0) { //At least the unperturbed state should have been generated
                System.err.println("Error applying tripeptide closure to residues " + resAffected[0] + " through " + resAffected[2]);
                System.exit(1);
            }

            //The unperturbed state needs to have parameter 0, so switch it to that position
            //We don't need matrices for the unperturbed state so those can be overwritten
            matrices[unperturbed] = matrices[0];
        }


        //Store the matrices and current perturbation parameters at the end of rotMatrices and predParams respectively
        float newPredParams[][] = new float[predParams.length+1][];
        float newRotMatrices[][][][][] = new float[rotMatrices.length+1][][][][];
        System.arraycopy(predParams,0,newPredParams,0,predParams.length);
        System.arraycopy(rotMatrices,0,newRotMatrices,0,rotMatrices.length);
        predParams = newPredParams;
        rotMatrices = newRotMatrices;

        predParams[predParams.length-1] = new float[predecessors.length];
        for(int b=0; b<predecessors.length; b++)
            predParams[predParams.length-1][b] = m.perts[predecessors[b]].curParam;

        rotMatrices[rotMatrices.length-1] = matrices;

        return numSoln;

    }



    public void applySSChanges() {
        //Add loop residues to adjacent secondary structures by copying the dihedrals of the nearest residues in those structures

        if( ! (addSSStart || addSSEnd) )
            return;//Nothing to do

        RamachandranChecker rcheck = RamachandranChecker.getInstance();
        RotMatrix rm = new RotMatrix();

        int closeStart=0;//Start of tripeptide to close

        float anchor1Old[] = null, anchor2Old[] = null;//Positions of anchor CAs before SS changes

        if(addSSStart) { //Add the first residue

            Residue res0 = m.residue[resAffected[0]];
            if(!m.checkNBonded(resAffected[0])) {
                System.err.println("ERROR: Trying to add " + res0.fullName +
                                   " to a preceding secondary structure, but it is not bonded to any preceding residue");
                System.exit(1);
            }

            Residue res1 = m.residue[resAffected[1]];

            anchor1Old = m.getActualCoord( res1.getAtomNameToMolnum("CA") );

            float phiPsi[] = rcheck.getPhiPsi(m, resAffected[0]-1);//phi, psi for residue before perturbation starts

            //First set phi
            int atomList1[] = res0.getAtomList(false, true, true, true);
            int atomList2[] = res1.getAtomList(true, true, true, false);
            int[] atomList = concatenateArrays(atomList1, atomList2);

            m.setTorsion( m.residue[resAffected[0]-1].getAtomNameToMolnum("C") ,
                          res0.getAtomNameToMolnum("N"), res0.getAtomNameToMolnum("CA"),
                          res0.getAtomNameToMolnum("C"), phiPsi[0], atomList, atomList.length, false);


            //Now psi
            atomList1 = res0.getAtomList(false, false, false, true);
            atomList2 = res1.getAtomList(true, true, true, false);
            atomList = concatenateArrays(atomList1, atomList2);

            m.setTorsion( res0.getAtomNameToMolnum("N"), res0.getAtomNameToMolnum("CA"),
                          res0.getAtomNameToMolnum("C"), res1.getAtomNameToMolnum("N"),
                          phiPsi[1], atomList, atomList.length, false);

            closeStart=1;
        }

        if(addSSEnd) {

            Residue resEnd = m.residue[resAffected[resAffected.length-1]];
            if(!m.checkCBonded(resAffected[resAffected.length-1])) {
                System.err.println("ERROR: Trying to add " + resEnd.fullName +
                                   " to an ensuing secondary structure, but it is not bonded to any ensuing residue");
                System.exit(1);
            }

            Residue resPrior = m.residue[resAffected[resAffected.length-2]];

            anchor2Old = m.getActualCoord( resPrior.getAtomNameToMolnum("CA") );

            float phiPsi[] = rcheck.getPhiPsi(m, resAffected[resAffected.length-1]+1 );
            //phi, psi for residue after perturbation ends

            //First set psi
            int atomList1[] = resEnd.getAtomList(true, true, true, false);
            int atomList2[] = resPrior.getAtomList(false, true, true, true);
            int[] atomList = concatenateArrays(atomList1, atomList2);

            m.setTorsion( m.residue[resAffected[resAffected.length-1]+1].getAtomNameToMolnum("N") ,
                          resEnd.getAtomNameToMolnum("C"), resEnd.getAtomNameToMolnum("CA"),
                          resEnd.getAtomNameToMolnum("N"), phiPsi[1], atomList, atomList.length, false);


            //Now phi
            atomList1 = resEnd.getAtomList(true, false, false, false);
            atomList2 = resPrior.getAtomList(false, true, true, true);
            atomList = concatenateArrays(atomList1, atomList2);

            m.setTorsion( resEnd.getAtomNameToMolnum("C"), resEnd.getAtomNameToMolnum("CA"),
                          resEnd.getAtomNameToMolnum("N"), resPrior.getAtomNameToMolnum("C"),
                          phiPsi[0], atomList, atomList.length, false);
        } else
            anchor2Old = m.getActualCoord( m.residue[resAffected.length-1].getAtomNameToMolnum("CA") );


        Residue anchorRes1 = m.residue[resAffected[closeStart]];
        Residue middleRes = m.residue[resAffected[closeStart+1]];
        Residue anchorRes2 = m.residue[resAffected[closeStart+2]];

        float anchor1New[] = m.getActualCoord( anchorRes1.getAtomNameToMolnum("CA") );
        float anchor2New[] = m.getActualCoord( anchorRes2.getAtomNameToMolnum("CA") );



        //Finally translate everything in between--the stuff to be closed--into place

        //The middle CA and everything before it is translated to avoid a discontinuity at the first CA of the tripeptide
        //Everything after is translated to avoid a discontinuity at the last CA of the tripeptide
        //But each translation only needs to happen if that CA has been messed with during secondary structure application


        if(addSSStart) {
            int atomList[] = concatenateArrays( anchorRes1.getAtomList(false, false, false, true),
                                                middleRes.getAtomList(true, true, true, false) );

            float tr[] = rm.subtract(anchor1New, anchor1Old);
            m.translateAtomList(atomList, tr, false, false);
        }

        if(addSSEnd) {
            int atomList[] = concatenateArrays( middleRes.getAtomList(false, false, false, true),
                                                anchorRes2.getAtomList(true, false, false, false) );

            float tr[] = rm.subtract(anchor2New, anchor2Old);
            m.translateAtomList(atomList, tr, false, false);
        }


        //Now the two rotations used to apply the loop closure adjustment will put the backbone in right
        //and sidechain idealization will fix the sidechains
    }


    private int[] concatenateArrays(int[]... toMerge) {
        //Concatenates integer arrays (e.g. atom lists)

        int fullSize = 0;

        for ( int[] arr : toMerge )
            fullSize += arr.length;

        int[] ans = new int[fullSize];

        int cursor = 0;

        for( int[] arr : toMerge ) {
            System.arraycopy(arr, 0, ans, cursor, arr.length);
            cursor += arr.length;
        }

        return ans;
    }






}
