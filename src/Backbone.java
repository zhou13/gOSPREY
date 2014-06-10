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
//	Backbone.java
//
//	Version:           2.1 beta
//
//
//	  authors:
// 	  initials    name                 organization                email
//	 ---------   -----------------    ------------------------    ----------------------------
//	  ISG		 Ivelin Georgiev	  Duke University			  ivelin.georgiev@duke.edu
//     KER        Kyle E. Roberts       Duke University         ker17@duke.edu
//     PGC        Pablo Gainza C.       Duke University         pablo.gainza@duke.edu
///////////////////////////////////////////////////////////////////////////////////////////////

/**
* Written by Ivelin Georgiev (2004-2009)
*
*/

import java.io.Serializable;

/**
 * Handles the backbone representation for the protein;
 * 		Applies (phi,psi) changes to the molecule;
 * 		Assumes that the order of the atoms for the phi angle is C(i-1), N(i), CA(i), C(i),
 * 			and for the psi angle: N(i), CA(i), C(i), N(i+1)
 *
 */
public class Backbone implements Serializable {

    //constructor
    Backbone() {

    }

    //Get the current phi (angleType==0) or psi (angleType=1) value for the given residue resNum of strand strandNum in molecule m
    public double getFiPsi(Molecule m, int strandNum, int resNum, int angleType) {

        if ((angleType==0)
                &&((resNum==0)||(m.strand[strandNum].residue[resNum-1].getResNumber()!=m.strand[strandNum].residue[resNum].getResNumber()-1))) { //phi angle needed and (N-terminus or no prev residue), so no phi angle
            return 0.0;
        } else if ((angleType==1)
                   &&((resNum==(m.strand[strandNum].numberOfResidues-1))||(m.strand[strandNum].residue[resNum+1].getResNumber()!=m.strand[strandNum].residue[resNum].getResNumber()+1))) { //psi angle needed and (C-terminus or no next residue), so no psi angle
            return 0.0;
        } else {

            Residue curRes = m.strand[strandNum].residue[resNum];
            Residue prevRes = null;
            Residue nextRes = null;

            Atom at[] = new Atom[4];

            if (angleType==0) { //find the atoms involved in the phi angle
                prevRes = m.strand[strandNum].residue[resNum-1];
                getFiAtoms(curRes,prevRes,at);
            } else { //find the atoms involved in the psi angle
                nextRes = m.strand[strandNum].residue[resNum+1];
                getPsiAtoms(curRes,nextRes,at);
            }

            return (at[3].torsion(at[0], at[1], at[2]));
        }
    }

    // Find the atoms involved in the phi angle and return in at[] and atNum[]
    private void getFiAtoms(Residue curRes, Residue prevRes, Atom at[]) {
        getAtforRes("C",prevRes,at,3); //first find the C atom of the previous residue
        getAtforRes("N",curRes,at,2); //then find the N, CA, and C atoms of the current residue
        getAtforRes("CA",curRes,at,1);
        getAtforRes("C",curRes,at,0);
    }

    // Find the atoms involved in the psi angle and return in at[] and atNum[]
    private void getPsiAtoms(Residue curRes, Residue nextRes, Atom at[]) {
        getAtforRes("N",curRes,at,0); //first find the N, CA, and C atoms of the current residue
        getAtforRes("CA",curRes,at,1);
        getAtforRes("C",curRes,at,2);
        getAtforRes("N",nextRes,at,3); //then find the N atom of the next residue
    }

    //Find the atName atom of the given residue and store in at[ind] and atNum[ind]
    private void getAtforRes (String atName, Residue res, Atom at[], int ind) {
        for(int q=0; q<res.numberOfAtoms; q++) {
            if (res.atom[q].name.equalsIgnoreCase(atName)) {
                at[ind] = res.atom[q];
            }
        }
    }

    public void applyFiPsi(Molecule m, int strandNum, int resNum, double alpha, int angleType) {

        if ((angleType==0)
                &&((resNum==0)||(m.strand[strandNum].residue[resNum-1].getResNumber()!=m.strand[strandNum].residue[resNum].getResNumber()-1))) { //phi angle needed and (N-terminus or no prev residue), so no phi angle
            System.out.println("Warning: invalid phi angle.");
            return;
        } else if ((angleType==1)
                   &&((resNum==(m.strand[strandNum].numberOfResidues-1))||(m.strand[strandNum].residue[resNum+1].getResNumber()!=m.strand[strandNum].residue[resNum].getResNumber()+1))) { //psi angle needed and (C-terminus or no next residue), so no psi angle
            System.out.println("Warning: invalid psi angle.");
            return;
        } else {

            Residue curRes = m.strand[strandNum].residue[resNum];
            Residue prevRes = null;
            Residue nextRes = null;
            if(angleType == 0)
                prevRes = m.strand[strandNum].residue[resNum-1];
            if(angleType == 1)
                nextRes = m.strand[strandNum].residue[resNum+1];

            int atomList[] = new int[m.numberOfAtoms]; //list of atoms to be rotated

            Atom at[] = new Atom[4];
            at[0] = at[1] = at[2] = at[3] = null;

            if (angleType==0) //find the atoms involved in the phi angle
                getFiAtoms(curRes,prevRes,at);
            else //find the atoms involved in the psi angle
                getPsiAtoms(curRes,nextRes,at);


            for(int q=0; q<atomList.length; q++) {
                atomList[q]=1;
            }
            atomList[at[1].moleculeAtomNumber]=0;
            atomList[at[2].moleculeAtomNumber]=0;

            // Now find all atoms that need to be rotated
            getAtomsToRotate(m,strandNum,at[2],atomList,angleType);

            // Copy atoms over in atomList, ie. if current atomList[q]==2 then
            //  the new atomList[] should include q (ie. atom q counts)
            // Skip the 4th atom of the torsion as it's automatically rotated
            //  by m.setTorsion
            int alSize=0;
            for(int q=0; q<atomList.length; q++) {
                if ((atomList[q]==2) && (q != at[3].moleculeAtomNumber))
                    atomList[alSize++] = q;
            }

            double curTorsion = at[3].torsion(at[0], at[1], at[2]); //the current torsion

            // Perform the rotation
            // at[0], at[1], and at[2] don't move, at[3] and the atoms
            //  in the atomList rotate
            m.setTorsion(at[0].moleculeAtomNumber,at[1].moleculeAtomNumber,
                         at[2].moleculeAtomNumber,at[3].moleculeAtomNumber,
                         (curTorsion+alpha),atomList,alSize,false);
        }
    }

    //Determines which atoms need to be rotated when the phi or psi change is applied, for a given strand strandNum in molecule m;
    //Only certain atoms for the residue to which atAtom belongs will be rotated;
    //All atoms belonging to the residues before this one will be rotated
    private void getAtomsToRotate(Molecule m, int strandNum, Atom atAtom, int atomList[], int angleType) {

        int curRes = atAtom.moleculeResidueNumber;
        getAtomsToRotateHelper(m,curRes,atAtom,atomList);

        for (int i=0; i<m.strand[strandNum].numberOfResidues; i++) {

            int resNum = m.strand[strandNum].residue[i].moleculeResidueNumber;

            boolean toRot = false;
            if ((angleType==0)&&(resNum<curRes)) //phi rotation
                toRot = true;
            else if ((angleType==1)&&(resNum>curRes)) //psi rotation
                toRot = true;

            if (toRot) {
                for (int j=0; j<m.residue[resNum].numberOfAtoms; j++) {
                    int atNum = m.residue[resNum].atom[j].moleculeAtomNumber;
                    if (atomList[atNum] == 1)
                        atomList[atNum] = 2;
                }
            }
        }
    }

    //Determines which atoms need to be rotated when the phi or psi change is applied for a given residue
    private void getAtomsToRotateHelper(Molecule m, int curRes, Atom atAtom, int atomList[]) {

        if (atAtom.bond != null) {
            for(int q=0; q<atAtom.bond.length; q++) {
                Atom secondAtom = m.atom[atAtom.bond[q]];
                if (secondAtom.moleculeResidueNumber==curRes) {
                    if (atomList[secondAtom.moleculeAtomNumber] == 1) {
                        atomList[secondAtom.moleculeAtomNumber] = 2;
                        getAtomsToRotateHelper(m,curRes,secondAtom,atomList);
                    }
                }
            }
        }
    }
}
