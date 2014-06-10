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
// Strand.java
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
///////////////////////////////////////////////////////////////////////////////////////////////

/*
 *
 * Major changes were made by Ryan Lilien (2001-2004);
 * 		minor changes by Ivelin Georgiev (2004-2009)
 *
 * Rewritten by Ryan Lilien based on code by Neill White
 * Many functions have been added, others removed, most have had
 *  at least some parts rewritten. Code rewrites have often been
 *  major to fix bugs or add functionality.
 *
 * Based on software copyrighted, 1999, by Neill White.
 *  The author hereby grants permission to use, copy, modify, and re-distribute
 *  this software and its documentation for any purpose, provided
 *  that existing copyright notices are retained in all copies and that this
 *  notice is included verbatim in any distributions. No written agreement,
 *  license, or royalty fee is required for any of the authorized uses.
 *  Modifications to this software may be distributed provided that
 *  the nature of the modifications are clearly indicated.
 *
 */

import java.io.Serializable;

/**
 * Handles functions and data associated with strands.
 */
public class Strand implements Serializable {

    String name = "";		// The strand name
    int	numberOfResidues=0;		// Number of residues
    int	numberOfAtoms=0;	// Number of atoms
    int	number=-1;			// The number of the current strand
    Residue	residue[];	// Array of residues in the strand
    boolean isProtein = false;	// Is this strand a protein?
    boolean rotTrans = false;

    Strand() {
        residue = new Residue[1];
    }

    Strand(Residue firstResidue) {
        residue = new Residue[1];
        residue[0] = firstResidue;
        numberOfResidues = 1;
    }

    Strand(String strandName) {
        numberOfResidues = 0;
        numberOfAtoms = 0;
        number = 0;
        residue = new Residue[1];
        isProtein = false;
        name = strandName;
    }

    // Displays some strand info to System.out for debugging
    public void printStrandInfo() {
        System.out.println("name = *"+name+"*");
        System.out.println("numberOfResides = *"+numberOfResidues+"*");
        System.out.println("numberOfAtoms = *"+numberOfAtoms+"*");
        System.out.println("number = *"+number+"*");
    }

    // Adds a residue to the strand
    // Best called by molecule.addResidue(.) but can be called
    //  directly. If called directly the following fields are not
    //  properly updated (as they can't be without molecule information)
    //  Residue.moleculeResidueNumber
    //  Atom.moleculeResidueNumber
    //  Atom.moleculeAtomNumber
    public int addResidue(Residue newResidue) {

        numberOfAtoms += newResidue.numberOfAtoms;

        // handle residue stuff
        Residue largerResidueArray[] = new Residue[numberOfResidues + 1];
        System.arraycopy(residue, 0, largerResidueArray, 0, residue.length);
        residue = largerResidueArray;

        newResidue.strandResidueNumber = numberOfResidues;
        newResidue.strandNumber = number;
        for(int i=0; i<newResidue.numberOfAtoms; i++) {
            newResidue.atom[i].strandNumber = number;
            newResidue.atom[i].strandResidueNumber = numberOfResidues;
        }
        residue[numberOfResidues] = newResidue;
        return(numberOfResidues++);
    }

    // Deletes a residue from the strand
    // All possible appropriate bookeeping is done including
    //  updates of bonds in this strand
    public int deleteResidue(int residueNumber) {

        int oldResNumAtoms = residue[residueNumber].numberOfAtoms;
        int lowAtomIndex = residue[residueNumber].atom[0].moleculeAtomNumber;
        int highAtomIndex = lowAtomIndex + oldResNumAtoms - 1;

        numberOfAtoms -= oldResNumAtoms;

        for(int i=residueNumber; i<numberOfResidues; i++) {
            residue[i].strandResidueNumber -= 1;
            for(int j=0; j<residue[i].numberOfAtoms; j++) {
                residue[i].atom[j].strandResidueNumber -= 1;
            }
        }
        Residue smallerResidueArray[] = new Residue[numberOfResidues-1];
        System.arraycopy(residue, 0, smallerResidueArray, 0, residueNumber);
        if (residueNumber < numberOfResidues-1)
            System.arraycopy(residue, residueNumber+1, smallerResidueArray,	residueNumber, residue.length-residueNumber-1);
        residue = smallerResidueArray;
        // update atom bond indices
        for(int i=0; i<numberOfResidues-1; i++) {
            for(int j=0; j<residue[i].numberOfAtoms; j++) {
                for(int m=0; m<residue[i].atom[j].numberOfBonds; m++) {
                    if (residue[i].atom[j].bond[m] >= lowAtomIndex) {
                        if (residue[i].atom[j].bond[m] <= highAtomIndex) {
                            residue[i].atom[j].deleteBond(m);
                            m--;
                        } else
                            residue[i].atom[j].bond[m] -= oldResNumAtoms;
                    }
                }
            }
        }

        return(--numberOfResidues);
    }

    // Adds an atom to the strand and the specified residue
    // All appropriate strand, residue, and atom fields are updated
    //  except for Atom.moleculeAtomNumber which we can't update
    //  It should be updated in the molecule class
    public int addAtom(int residueNumber, Atom newAtom) {
        // allow the residue to do its bookkeeping
        residue[residueNumber].addAtom(newAtom);
        return(numberOfAtoms++);
    }

    // Deletes an atom from the specified residue
    // Appropriate bookkeeping is done
    public int deleteAtom(int residueNumber, int atomNumber) {

        int moleculeAtomNumber = residue[residueNumber].atom[atomNumber].moleculeAtomNumber;

        if(residue[residueNumber].numberOfAtoms == 1)
            return(deleteResidue(residueNumber));

        residue[residueNumber].deleteAtom(atomNumber);

        // update atom bond indices
        for(int i=0; i<numberOfResidues; i++) {
            if (!(residue[i].strandResidueNumber==residueNumber)) {
                for(int j=0; j<residue[i].numberOfAtoms; j++) {
                    for(int m=0; m<residue[i].atom[i].numberOfBonds; m++) {
                        if (residue[i].atom[i].bond[m] == moleculeAtomNumber) {
                            residue[i].atom[i].deleteBond(m);
                            m--;
                        } else if (residue[i].atom[i].bond[m] > moleculeAtomNumber) {
                            residue[i].atom[i].bond[m] -= 1;
                        }
                    }
                }
            }
        }

        return(numberOfAtoms--);
    }

    // Returns the center of mass
    public float[] getCenterOfMass() {
        float centOfMass[] = new float[3];
        float xCenter = 0.0f, yCenter = 0.0f, zCenter = 0.0f;
        double totalMass = 0.0;
        double tmpMass = 0.0;

        for(int j=0; j<numberOfResidues; j++) {
            for(int i=0; i<residue[j].numberOfAtoms; i++) {
                tmpMass = residue[j].atom[i].mass;
                if (tmpMass==0)
                    System.out.println("** Zero mass atom detected???");
                totalMass += tmpMass;
                xCenter += (residue[j].atom[i].coord[0] * tmpMass);
                yCenter += (residue[j].atom[i].coord[1] * tmpMass);
                zCenter += (residue[j].atom[i].coord[2] * tmpMass);
            }
        }
        xCenter /= totalMass;
        yCenter /= totalMass;
        zCenter /= totalMass;
        centOfMass[0] = xCenter;
        centOfMass[1] = yCenter;
        centOfMass[2] = zCenter;
        return(centOfMass);
    }


    // This function rotates the entire strand by thetaDeg
    //  degrees around axis dx, dy, dz (around the center of mass)
    // If you don�t have a molecule pass null
    public void rotateAroundCOM(double dx, double dy, double dz,
                                double thetaDeg, Molecule m) {

        // Get center of mass
        float[] centOfMass = getCenterOfMass();

        rotateStrand(dx,dy,dz,thetaDeg,m,centOfMass[0],centOfMass[1],centOfMass[2]);
    }

    // This function rotates the entire strand by thetaDeg
    //  degrees around axis dx, dy, dz (around the point cx,cy,cz)
    // If you don�t have a molecule pass null
    public void rotateStrand(double dx, double dy, double dz,
                             double thetaDeg, Molecule m, float cx, float cy, float cz) {

        float fx,fy,fz, tx,ty,tz;
        fx = (new Double(dx)).floatValue();
        fy = (new Double(dy)).floatValue();
        fz = (new Double(dz)).floatValue();

        float[][] rot_mtx = new float[3][3];
        RotMatrix rM = new RotMatrix();
        rM.getRotMatrix(fx,fy,fz,(float) thetaDeg,rot_mtx);

        for(int w=0; w<numberOfResidues; w++) {
            for(int q=0; q<residue[w].numberOfAtoms; q++) {
                tx=residue[w].atom[q].coord[0] - cx;
                ty=residue[w].atom[q].coord[1] - cy;
                tz=residue[w].atom[q].coord[2] - cz;

                residue[w].atom[q].coord[0] = tx * rot_mtx[0][0] + ty * rot_mtx[0][1] + tz * rot_mtx[0][2] + cx;
                residue[w].atom[q].coord[1] = tx * rot_mtx[1][0] + ty * rot_mtx[1][1] + tz * rot_mtx[1][2] + cy;
                residue[w].atom[q].coord[2] = tx * rot_mtx[2][0] + ty * rot_mtx[2][1] + tz * rot_mtx[2][2] + cz;

                if (m!=null)
                    m.updateCoordinates(residue[w].atom[q]);
            }
        }
    }

    //Returns the strand-relative residue number of the residue with PDB residue number pdbResNum
    public int mapPDBresNumToStrandResNum(int pdbResNum) {
        for (int i=0; i<numberOfResidues; i++) {
            if (residue[i].getResNumber()==pdbResNum) {
                return residue[i].strandResidueNumber;
            }
        }
        return -1;
    }
}
