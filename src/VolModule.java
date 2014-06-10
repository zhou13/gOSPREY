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
// VolModule.java
//
//  Version:           2.1 beta
//
//
// authors:
//    initials    name            organization                email
//   ---------   --------------  ------------------------    ------------------------------
//     RHL        Ryan Lilien     Dartmouth College           ryan.lilien@dartmouth.edu
//
////////////////////////////////////////////////////////////////////////////////////////////

/*
 *
 * Written by Ryan Lilien (2001-2004)
 *
 * The Vol Module class can crudely compute a molecular volume for a
 *  specified molecule.
 * This class snags molecule information upon construction, if the
 *  molecule changes outside the class a new SASA object should be
 *  created.
 * To use this class one needs to set the molecule and each residue to
 *  either be on, off, or partially on (see setResidueTreatment and
 *  comment by residueTreat[])
 *
 * As all helper classes should, we deal with the molecule's actualCoordinates
 *  and we leave the actual atom.coord alone
 */

import java.lang.Math;
import java.util.Hashtable;
import java.util.Vector;
import java.io.PrintStream;
import java.io.BufferedInputStream;

/**
 * This class computes a crude molecular volume for a specified molecule.
 */
public class VolModule {

    Molecule m;					// Molecule for which this surface exists
    int numRes;					// Number of residues in molecule m
    int residueTreat[];			// Array of residues, 1 = surface should be computed
    //  for this residue, 2 = atoms should be treated as
    //  existing atoms in terms of computing surfaces for
    //  residues with a 1 but no molecular surface should
    //  be computed for them. 0 = no molecular surface
    //  should be generated and these atoms are invisible
    //  to other atoms.
    Hashtable atomBuckets;	// Will become a hashtable of Vectors of atoms
    //  in the search area
    Vector spacePoints;			// Vector of surface space (double[3]) the x,y,z coordinate
    double spaceRes;				// Resolution of the space field (in Angstroms)
    //  This is the spacing between points
    double bbScale;					// The bounding box range. If this was 80% of full
    //  size then if the involved residues were from -5 to
    //  5 this would return -4 to 4. This is done so that
    //  we don't get extra space points just outside the
    //  molecule
    double solventRadius;		// Radius of solvent sphere (in Angstroms)
    double bucketResolution;		// Resolution of buckets datastructure
    //  one shouldn't need to change this
    // 9 is good because the largest atoms have a
    //  radii of 2.4A and a typical solvent radius is
    //  1.5A for a total of 4A radius or 8A spacing


    // ****************************************
    // Constructor specifying the active molecule
    VolModule(Molecule mol) {

        solventRadius = 1.4;
        bucketResolution = 9.0;
        spaceRes = 0.5;
        bbScale = 0.80;
        m=mol;
        numRes=m.numberOfResidues;
        residueTreat = new int[numRes];
        for(int q=0; q<numRes; q++) {
            residueTreat[q]=0;
        }
    }


    // *************************************************************
    // Functions to set parameters for surface and space generation
    public void setResidueTreatment (int resNum, int treatment) {
        if ((resNum < numRes) && (treatment >= 0) && (treatment <=2))
            residueTreat[resNum]=treatment;
    }
    public void setSpaceResolution (double resol) {
        if (resol > 0)
            spaceRes = resol;
    }
    public void setBoundingBoxScale (double scale) {
        if (scale > 0)
            bbScale = scale;
    }
    public void setSolventRadius (double radius) {
        if (radius > 0.0)
            solventRadius = radius;
    }
    public void setBucketResolution (double spacing) {
        if (spacing > 0.0)
            bucketResolution = spacing;
    }

    // *************************************************************
    // Generates an atom space field, which is a dotted space corresponding
    //  to the non-SAS space. One can make the solventRadius = 0 to examine
    //  only the vdw space.
    public void generateSpace() {

        double bb[];
        int xpos,ypos,zpos;
        Vector vec;

        // Generate atomBuckets for rapid atom steric indexing
        // I could check if a buckets hashTable already exists but
        //  for simplicity I'll just make a new one for now.
        createAtomBuckets();

        // Compute the bounding box range. If this was 80% of full
        //  size then if the molecule extent was from -5 to 5 this
        //  would return -4 to 4. This is done so that we don't
        //  get extra space points just outside the molecule
        bb = getBoundingBox(bbScale);

        // Search the bounding region and check points against appropriate
        //  neighboring atoms for steric overlap
        double newspacept[];
        int numAtomsToCheck;
        spacePoints = new Vector();
        int neighborAt[] = new int[m.numberOfAtoms];
        for(double q=bb[0]; q<bb[1]; q+=spaceRes) {
            for(double w=bb[2]; w<bb[3]; w+=spaceRes) {
                for(double j=bb[4]; j<bb[5]; j+=spaceRes) {
                    xpos = (int)Math.round(q / bucketResolution);
                    ypos = (int)Math.round(w / bucketResolution);
                    zpos = (int)Math.round(j / bucketResolution);
                    numAtomsToCheck = 0;
                    for(int xi=xpos-1; xi<=xpos+1; xi++) {
                        for(int yi=ypos-1; yi<=ypos+1; yi++) {
                            for(int zi=zpos-1; zi<=zpos+1; zi++) {
                                String key = new String(xi+","+yi+","+zi);
                                vec = (Vector)atomBuckets.get(key);
                                if (!(vec==null))
                                    for(int r=0; r<vec.size(); r++)
                                        neighborAt[numAtomsToCheck++] = ((Integer)vec.elementAt(r)).intValue();
                            }
                        }
                    }
                    newspacept = new double[4];
                    newspacept[0] = q;
                    newspacept[1] = w;
                    newspacept[2] = j;
                    if (!tooClose(newspacept,neighborAt,numAtomsToCheck,-1,0.95f)) {
                        newspacept[3] = 0.0;
                        spacePoints.addElement(newspacept);
                    }
                }
            }
        }

    } // end generateSpace


    // Returns the internal molecular volume of all residues
    //  selected as either type 1 or 2
    // This is done in a manner very similar to generateSpace
    // A defulat high resolution is used
    public float getMoleculeVolume() {
        return(getMoleculeVolume(0.25f,0.0f));
    }

    // Returns the internal molecular volume of all residues
    //  selected as either type 1 or 2
    // This is done in a manner very similar to generateSpace
    // The grid spacing is volRes
    // theSolvRadius is the solvent radius to use for this run
    public float getMoleculeVolume(float volRes, float theSolvRadius) {

        double bb[];
        int xpos,ypos,zpos;
        Vector vec;
        int numPoints = 0;  // the number of points inside the molecule
        double heldSolventRadius = solventRadius;

        // Save and update the solventRadius to be 0 for this calculation
        solventRadius = theSolvRadius;

        // Generate atomBuckets for rapid atom steric indexing
        // I could check if a buckets hashTable already exists but
        //  for simplicity I'll just make a new one for now.
        createAtomBuckets();

        // Compute the bounding box range. We use an extend of
        //  105% so that we're just bigger than the molecule
        bb = getBoundingBox(1.05);
        // Extend bond a bit so that we get all parts of the
        //  vdw sphere, since no atom has a vdw radii of over
        //  5A this should be fine.
        bb[0] -= 5.0;
        bb[1] += 5.0;
        bb[2] -= 5.0;
        bb[3] += 5.0;
        bb[4] -= 5.0;
        bb[5] += 5.0;

        // Search the bounding region and check points against appropriate
        //  neighboring atoms for steric overlap
        int numAtomsToCheck;
        double tempPt[] = new double[3];
        int neighborAt[] = new int[m.numberOfAtoms];
        for(double q=bb[0]; q<bb[1]; q+=volRes) {
            for(double w=bb[2]; w<bb[3]; w+=volRes) {
                for(double j=bb[4]; j<bb[5]; j+=volRes) {
                    xpos = (int)Math.round(q / bucketResolution);
                    ypos = (int)Math.round(w / bucketResolution);
                    zpos = (int)Math.round(j / bucketResolution);
                    numAtomsToCheck = 0;
                    for(int xi=xpos-1; xi<=xpos+1; xi++) {
                        for(int yi=ypos-1; yi<=ypos+1; yi++) {
                            for(int zi=zpos-1; zi<=zpos+1; zi++) {
                                String key = new String(xi+","+yi+","+zi);
                                vec = (Vector)atomBuckets.get(key);
                                if (!(vec==null))
                                    for(int r=0; r<vec.size(); r++)
                                        neighborAt[numAtomsToCheck++] = ((Integer)vec.elementAt(r)).intValue();
                            }
                        }
                    }
                    tempPt[0] = q;
                    tempPt[1] = w;
                    tempPt[2] = j;
                    if (tooClose(tempPt,neighborAt,numAtomsToCheck,-1,1.0f)) {
                        numPoints++;
                    }
                }
            }
        }

        // Restore the solventRadius
        solventRadius = heldSolventRadius;

        // determine the total number of potential points in
        //  the bounding region
        double possibleVolume = (bb[1]-bb[0]) * (bb[3]-bb[2]) * (bb[5]-bb[4]);
        double possiblePoints = possibleVolume / (volRes * volRes * volRes);
        return((float)(numPoints / possiblePoints * possibleVolume));
    }


    // *************************************************************
    // Gets x,y,z extremes and then scales them by the scale factor
    private double[] getBoundingBox(double scale) {

        double bbox[] = new double[6];
        double axd,ayd,azd;
        float coords[] = new float[3];
        int tmpNum = 0;

        bbox[0]=999999.0;		// xmin
        bbox[1]=-999999.0;	// xmax
        bbox[2]=999999.0;		// ymin
        bbox[3]=-999999.0;	// ymax
        bbox[4]=999999.0;		// zmin
        bbox[5]=-999999.0;	// zmax

        for(int q=0; q<m.numberOfAtoms; q++) {
            tmpNum = m.atom[q].moleculeAtomNumber * 3;
            coords[0] = m.actualCoordinates[tmpNum];
            coords[1] = m.actualCoordinates[tmpNum+1];
            coords[2] = m.actualCoordinates[tmpNum+2];
            if (coords[0] < bbox[0])
                bbox[0] = coords[0];
            else if (coords[0] > bbox[1])
                bbox[1] = coords[0];
            if (coords[1] < bbox[2])
                bbox[2] = coords[1];
            else if (coords[1] > bbox[3])
                bbox[3] = coords[1];
            if (coords[2] < bbox[4])
                bbox[4] = coords[2];
            else if (coords[2] > bbox[5])
                bbox[5] = coords[2];
        }

        axd = (bbox[1]+bbox[0])/2;
        ayd = (bbox[3]+bbox[2])/2;
        azd = (bbox[5]+bbox[4])/2;

        bbox[0]=axd-(axd-bbox[0])*scale;
        bbox[1]=axd+(bbox[1]-axd)*scale;
        bbox[2]=ayd-(ayd-bbox[2])*scale;
        bbox[3]=ayd+(bbox[3]-ayd)*scale;
        bbox[4]=azd-(azd-bbox[4])*scale;
        bbox[5]=azd+(bbox[5]-azd)*scale;

        return(bbox);
    }


    // *************************************************************
    // Creates the atom buckets hashtable for use with generateSpace
    private void createAtomBuckets() {

        atomBuckets = new Hashtable();
        int xpos,ypos,zpos,tmpNum;
        Vector vec;

        // Fill buckets with atoms for rapid indexing in later stages
        for(int q=0; q<numRes; q++) {
            if ((residueTreat[q] == 1) || (residueTreat[q] == 2)) {
                for(int w=0; w<m.residue[q].numberOfAtoms; w++) {
                    tmpNum = m.residue[q].atom[w].moleculeAtomNumber * 3;
                    xpos = (int)Math.round(m.actualCoordinates[tmpNum] / bucketResolution);
                    ypos = (int)Math.round(m.actualCoordinates[tmpNum+1] / bucketResolution);
                    zpos = (int)Math.round(m.actualCoordinates[tmpNum+2] / bucketResolution);
                    String key = new String(xpos+","+ypos+","+zpos);
                    vec = (Vector)atomBuckets.get(key);
                    if (vec == null)
                        vec = new Vector();
                    vec.addElement(new Integer(m.residue[q].atom[w].moleculeAtomNumber));
                    atomBuckets.put(key,vec);
                }
            }
        }
    }


    // *************************************************************
    // Returns true if there is an atom too close to a specified new surface point
    private boolean tooClose(double newsurfpt[], int neighborAt[], int numAtomsToCheck,
                             int moleculeAtomNumber, float threshMultiplier) {

        int atNum,tmpNum;
        float atCoords[] = new float[3];
        double thresh;
        for(int i=0; i<numAtomsToCheck; i++) {
            atNum = neighborAt[i];
            if (atNum != moleculeAtomNumber) {
                tmpNum = m.atom[atNum].moleculeAtomNumber * 3;
                atCoords[0] = m.actualCoordinates[tmpNum];
                atCoords[1] = m.actualCoordinates[tmpNum+1];
                atCoords[2] = m.actualCoordinates[tmpNum+2];
                // Threshold is threshMultiplier * the distance of the atomic radius
                //  plus the solvent radius
                thresh = threshMultiplier * (m.atom[atNum].radius/100.0 + solventRadius);
                thresh = thresh * thresh;
                if ( ((newsurfpt[0]-atCoords[0])*(newsurfpt[0]-atCoords[0]) +
                        (newsurfpt[1]-atCoords[1])*(newsurfpt[1]-atCoords[1]) +
                        (newsurfpt[2]-atCoords[2])*(newsurfpt[2]-atCoords[2])) < thresh )
                    return (true);
            }
        }
        return(false);
    }


} // end VolModule class
