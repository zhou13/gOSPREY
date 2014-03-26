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
//	Backrubs.java
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
 * Handles the application of the Richardsons' Backrub motions for a given residue in a protein
 * 
 */
public class Backrubs implements Serializable {
	
	final float thetaSmallScale = 0.7f; //the scaling factor for the rotation angles for the small rotations
	
	//constructor
	Backrubs(){
		
	}
	
	//Applies a Backrub motion to residue resNum (strand-relative numbering) of strand strandNum in molecule m;
	//This motion includes three rotations:
	//	1) rotate resNum and its adjacent peptides by theta1 degrees around the CA(resNum-1)-CA(resNum+1) axis;
	//	2) rotate the peptide with resNum-1 by theta2 degrees around the CA(resNum-1)-CA(resNum) axis;
	//  3) rotate the peptide with resNum+1 by theta3 degrees around the CA(resNum)-CA(resNum+1) axis;
	//	The only atoms that move from their initial positions belong to residue resNum and the adjacent peptides
	//	If computeSmallRot is true, then the values of theta2 and theta3 are computed after the first rotation, 
	//		so that the H-bonding positions of CO and NH are roughly preserved
	public float[] applyBackrub(Molecule m, int strandNum, int resNum, float theta1, boolean computeSmallRot, float theta2, float theta3){
		
		if ( ((resNum==0)||(m.strand[strandNum].residue[resNum-1].getResNumber()!=m.strand[strandNum].residue[resNum].getResNumber()-1))
				|| ((resNum==(m.strand[strandNum].numberOfResidues-1))||(m.strand[strandNum].residue[resNum+1].getResNumber()!=m.strand[strandNum].residue[resNum].getResNumber()+1)) ) {
			System.out.println("ERROR: Backrubs cannot be applied for residue "+m.strand[strandNum].residue[resNum].getResNumber());
			System.exit(1);
		}
		
		Residue curRes = m.strand[strandNum].residue[resNum];
		Residue prevRes = m.strand[strandNum].residue[resNum-1];
		Residue nextRes = m.strand[strandNum].residue[resNum+1];
		
		//Get the three lists of atoms to rotate (the atom lists are in molecule-relative numbering): 
		// 	1) C, O (prevRes); all atoms (curRes); N, H (nextRes);
		// 	2) C, O (prevRes); N, H (curRes);
		// 	3) C, O (curRes); N, H (nextRes)
		int atomList1[] = new int[curRes.numberOfAtoms+4]; //list of atoms to be rotated
		int atomList2[] = new int[4];
		int atomList3[] = new int[4];
		
		Atom Cprev = getAtForRes("C",prevRes);
		Atom Oprev = getAtForRes("O",prevRes);
		Atom Nnext = getAtForRes("N",nextRes);
		Atom Hnext = getAtForRes("H",nextRes);
		if(nextRes.name.equalsIgnoreCase("PRO"))
			Hnext = getAtForRes("CD",nextRes);
		Atom Ccur = getAtForRes("C",curRes);
		Atom Ocur = getAtForRes("O",curRes);
		Atom Ncur = getAtForRes("N",curRes);
		Atom Hcur = getAtForRes("H",curRes);
		Atom CAprev = getAtForRes("CA",prevRes);
		Atom CAcur = getAtForRes("CA",curRes);
		Atom CAnext = getAtForRes("CA",nextRes);
		
		for (int i=0; i<curRes.numberOfAtoms; i++)
			atomList1[i] = curRes.atom[i].moleculeAtomNumber;
		atomList1[curRes.numberOfAtoms] = Cprev.moleculeAtomNumber;
		atomList1[curRes.numberOfAtoms + 1] = Oprev.moleculeAtomNumber;
		atomList1[curRes.numberOfAtoms + 2] = Nnext.moleculeAtomNumber;
		atomList1[curRes.numberOfAtoms + 3] = Hnext.moleculeAtomNumber;
		
		atomList2[0] = Cprev.moleculeAtomNumber;
		atomList2[1] = Oprev.moleculeAtomNumber;
		atomList2[2] = Ncur.moleculeAtomNumber;
		atomList2[3] = Hcur.moleculeAtomNumber;
		
		atomList3[0] = Ccur.moleculeAtomNumber;
		atomList3[1] = Ocur.moleculeAtomNumber;
		atomList3[2] = Nnext.moleculeAtomNumber;
		atomList3[3] = Hnext.moleculeAtomNumber;
		
		float OprevOldCoord[] = getActualCoord(m, Oprev);
		Atom pseudoOprevOld = new Atom(Oprev.name,OprevOldCoord[0],OprevOldCoord[1],OprevOldCoord[2]);
		float OcurOldCoord[] = getActualCoord(m, Ocur);
		Atom pseudoOcurOld = new Atom(Ocur.name,OcurOldCoord[0],OcurOldCoord[1],OcurOldCoord[2]);

		//Perform the big rotation
		rotateList(m,atomList1,CAprev,CAnext,CAprev,theta1);
		
		//Determine the rotation angles for the two smaller rotations, in order to roughly preserve the H-bonding positions of CO and NH
		if (computeSmallRot){
			theta2 = getSmallRotAngle(m,Oprev,CAcur,CAprev,pseudoOprevOld);
			theta2 *= thetaSmallScale;
			if (Math.signum(theta2)==Math.signum(theta1))
				theta2 = -theta2;
			
			theta3 = getSmallRotAngle(m,Ocur,CAcur,CAnext,pseudoOcurOld);
			theta3 *= thetaSmallScale;
			if (Math.signum(theta3)==Math.signum(theta1))
				theta3 = -theta3;
		}
		
		//Performs the two smaller rotations
		rotateList(m,atomList2,CAprev,CAcur,CAprev,theta2);
		rotateList(m,atomList3,CAcur,CAnext,CAcur,theta3);
		
		if (computeSmallRot){
			float theta23[] = new float[2];
			theta23[0] = theta2;
			theta23[1] = theta3;
			return theta23;
		}
		return null;
	}
	
	//Rotates all atoms in atomList[] around the axis between atoms a1-a2, centered at atom c, by theta degrees
	private void rotateList(Molecule m, int atomList[], Atom a1, Atom a2, Atom c, float theta){
		
		int a1num = a1.moleculeAtomNumber;
		int a2num = a2.moleculeAtomNumber;
		int cnum = c.moleculeAtomNumber;
		
		float dx = m.actualCoordinates[a2num*3] - m.actualCoordinates[a1num*3];
		float dy = m.actualCoordinates[a2num*3+1] - m.actualCoordinates[a1num*3+1];
		float dz = m.actualCoordinates[a2num*3+2] - m.actualCoordinates[a1num*3+2];
		
		double center[] = new double[3];
		center[0] = m.actualCoordinates[cnum*3];
		center[1] = m.actualCoordinates[cnum*3+1];
		center[2] = m.actualCoordinates[cnum*3+2];
		
		m.rotateAtomList(atomList, dx, dy, dz, center[0], center[1], center[2], theta, false);
	}
	
	//Find the atName atom of the given residue and store in at[ind] and atNum[ind]
	private Atom getAtForRes (String atName, Residue res){
		for(int q=0;q<res.numberOfAtoms;q++) {
			if (res.atom[q].name.equalsIgnoreCase(atName)) {
				return res.atom[q];
			}	
		}
		return null;
	}
	
	//Returns the coordinates in actualCoordinates[] for atom a from molecule m
	private float [] getActualCoord(Molecule m, Atom a){
		
		int anum = a.moleculeAtomNumber;
		
		float coord[] = new float[3];
		for (int i=0; i<3; i++)
			coord[i] = m.actualCoordinates[anum*3+i];
		
		return coord;
	}
	
	//Get the small rotation angle that will rotate atom a1 around the axis defined by atoms (a2,a3), so that a1 will be as close as possible to atom a4
	private float getSmallRotAngle(Molecule m, Atom a1, Atom a2, Atom a3, Atom a4){
		
		float a1Coord[] = getActualCoord(m, a1);
		Atom pp1 = new Atom(a1.name,a1Coord[0],a1Coord[1],a1Coord[2]);
		float a2Coord[] = getActualCoord(m,a2);
		Atom pp2 = new Atom(a2.name,a2Coord[0],a2Coord[1],a2Coord[2]);
		Atom pp3 = projectPointLine(a3, pp2, pp1);
		Atom pp4 = projectPointPlane(a3, pp2, pp3, a4);
		Atom closestPoint = getClosestPoint(pp3,pp1,pp4);
		return (float)closestPoint.angle(pp1, pp3);
	}


        //The next three functions have been changed from private to public so Backrub can use them

	
	//Returns the projection (a pseudo-atom) of atom p3 onto the line between atoms l1 and l2
	public Atom projectPointLine(Atom l1, Atom l2, Atom p3){
		
		float c[] = new float[3];
		double d12sq = Math.pow(l2.distance(l1),2);
		double u = ( (p3.coord[0]-l1.coord[0])*(l2.coord[0]-l1.coord[0]) + (p3.coord[1]-l1.coord[1])*(l2.coord[1]-l1.coord[1]) + (p3.coord[2]-l1.coord[2])*(l2.coord[2]-l1.coord[2]) );
		u = u / d12sq;
		
		c[0] = (float)(l1.coord[0] + u * (l2.coord[0]-l1.coord[0]));
		c[1] = (float)(l1.coord[1] + u * (l2.coord[1]-l1.coord[1]));
		c[2] = (float)(l1.coord[2] + u * (l2.coord[2]-l1.coord[2]));
		
		return (new Atom("",c[0],c[1],c[2]));
	}
	
	//Returns the projection (a pseudo-atom) of atom p4 onto the plane with normal defined by atoms (l1,l2) and a point on that plane c3
	public Atom projectPointPlane(Atom l1, Atom l2, Atom c3, Atom p4){
		
		float l[] = new float[3];
		for (int i=0; i<3; i++)
			l[i] = l2.coord[i] - l1.coord[i];
		
		float d = 0.0f;
		for (int i=0; i<3; i++)
			d += l[i]*c3.coord[i];
		
		float t1 = 0.0f;
		for (int i=0; i<3; i++)
			t1 += l[i]*p4.coord[i];
		t1 -= d;
		
		float t2 = 0.0f;
		for (int i=0; i<3; i++)
			t2 += l[i]*l[i];
		
		float t = t1/t2;
		
		float r[] = new float[3];
		for (int i=0; i<3; i++)
			r[i] = p4.coord[i] - t*l[i];
		
		return (new Atom("",r[0],r[1],r[2]));
	}
	
	//Compute the closest point on the circle defined by the Atom c (center) and radius (c,p1) to Atom q1 (this assumes c, p1, and q1 are coplanar)
	public Atom getClosestPoint (Atom c, Atom p1, Atom q1){
		
		float r = (float)c.distance(p1);
		
		float t[] = new float[3];
		for (int i=0; i<3; i++)
			t[i] = q1.coord[i] - c.coord[i];
		
		float d = 0.0f;
		for (int i=0; i<3; i++)
			d += t[i]*t[i];
		d = (float)Math.sqrt(d);
		
		float a[] = new float[3];
		for (int i=0; i<3; i++)
			a[i] = c.coord[i] + r*(t[i]/d);
		
		return (new Atom("",a[0],a[1],a[2]));
	}
}
