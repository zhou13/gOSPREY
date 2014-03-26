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
// StrandRotamers.java
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

/** 
 * Written by Ryan Lilien (2001-2004) and Ivelin Georgiev (2004-2009)
 * 
 */

import java.io.*;
import java.util.*;

/**
 * This class handles rotamer assignment and maintenance for a given strand;
 * 		Performs rotamer swaps and amino acid mutations fot this strand.
 */
public class StrandRotamers implements Serializable {
	
	// If the debug flag is set to true then additional debug statements are
	//  printed to standard out.
	public static final boolean debug = false;
	
	public RotamerLibrary rl = null; //the rotamer library object
	
	protected int strandNumber = -1; 		//the strand number
	private int totalNumDihedrals = 0;			//number of dihedrals in this strand
	protected int numberOfResidues=0;			// Number of residues in this strand
	private int allowableAAs[][];		// Which AA types are allowable at each position
	private int numAllowable[];			// Number of allowables
	private String curAAType[] = null;    //the three-letter code currently assumed by each residue in the strand
	protected int curRotNum[] = null;     //the rotamer number currently assumed by each residue in the strand
	
	// Generic constructor
	StrandRotamers(String rotFilename, Strand s) {
		
		rl = new RotamerLibrary(rotFilename,s.isProtein);
		setupStrand(s);
		
		
	}
	
	// Generic constructor
	StrandRotamers(RotamerLibrary rotLib, Strand s) {
		
		rl = rotLib;
		setupStrand(s);
	}
	
	public StrandRotamers reInit(Strand s){
		StrandRotamers sr = new StrandRotamers(rl,s);
		return sr;
	}

	// This function sets up the allowableAAs[][] array
	//  as well as the numAllowable[], curRotNum[], and
	//  curAAType[] arrays.
	private void setupStrand(Strand s) {
		
		strandNumber = s.number;
		
		numberOfResidues = s.numberOfResidues;
		allowableAAs = new int[numberOfResidues][rl.getNumAAallowed()];
		numAllowable = new int[numberOfResidues];
		curRotNum = new int[numberOfResidues];
		curAAType = new String[numberOfResidues];

		// Clear allowable amino acid arrays
		for(int i=0;i<numberOfResidues;i++) {
			for(int q=0;q<rl.getNumAAallowed();q++) {
				allowableAAs[i][q] = 0;
			}
			numAllowable[i] = 0;
			curAAType[i] = s.residue[i].name;
			curRotNum[i] = -1;
		}
	}

	// This function adds the AA type named name to the list
	//  of allowable types for residue number resNum
	public void setAllowable(int resNum, String name) {

		int aaNum = rl.getAARotamerIndex(name);
		if (aaNum>=0) {
			if (allowableAAs[resNum][aaNum] == 0) {
				allowableAAs[resNum][aaNum] = 1;
				numAllowable[resNum]++;
			}
		}
		
		if(debug){
			if (aaNum>=0)
				System.out.println("setAllow: res " + resNum + " type " + name);
		}
	}

        // This function checks if the AA type named name is allowed for residue number resNum
	public boolean checkAllowable(int resNum, String name) {

		int aaNum = rl.getAARotamerIndex(name);
		if (aaNum>=0) {
			if (allowableAAs[resNum][aaNum] == 1)
                            return true;
                        else
                            return false;
		}
                else{
                    System.err.println("Warning: AA type " + name + " not recognized");
                    return false;
                }

	}

	// This function clears the list of allowable AA types
	//  for residue number resNum
	public void clearAllowable(int resNum) {

		for(int i=0;i<rl.getNumAAallowed();i++)
			allowableAAs[resNum][i] = 0;
		numAllowable[resNum] = 0;
	}

	// Returns the number of allowable residue types
	//  for the specified residue
	public int getNumAllowable(int resIndex){
		return(numAllowable[resIndex]);
	}

	// Returns the numberOfResidues
	public int getNumResidues(){
		return(numberOfResidues);
	}
	
	// Returns the current three-letter code for the specified residue index strResPos (strand-relative numbering)
	public String getCurRotType(int strResPos){
		return(curAAType[strResPos]);
	}

	// Returns the current rotamer number for the specified residue index strResPos (strand-relative numbering)
	public int getCurRotNum(int strResPos){
		return(curRotNum[strResPos]);
	}
	
	// Returns the sum volume of the currently assumed
	//  rotamers for all residues marked with at least
	//  one allowable amino acid type
	// If no rotamer has been assigned (ie. if a mutation
	//  was just performed but no 'rotamer' was assigned
	//  and curRotNum[] is -1) then the volume of the first
	//  rotamer is used.
	public float getCurAllowableVolume(){
	
		float vol = 0.0f;
		int localAANum = -1;
		String tmpstg = null;
		
		float rotamerVolumes[][] = rl.getRotVol();
		
		for(int i=0;i<numberOfResidues;i++){
			if(numAllowable[i] > 0){
				localAANum = rl.getAARotamerIndex(curAAType[i]);
				int tRot = curRotNum[i];
				if (tRot == -1)
					tRot = 0;
				vol += rotamerVolumes[localAANum][tRot];
			}
		}	
		return(vol);
	}
	

	// This function returns the index of num-th allowable
	//  AA for the specified residue
	// allowableAAs is a 0/1 array
	public int getIndexOfNthAllowable(int resNum, int num) {
		
		int tmp = -1;
		for(int q=0;q<rl.getNumAAallowed();q++) {
			tmp += allowableAAs[resNum][q];
			if (num == tmp)
				return q;
		}
		return -1;
	}

	// This function returns a list of atoms from residue resNum that are
	//  further distal than at2, atomList for a2 == 0 so we ignore it
	// atomList elements are 0=ignore, 1=unchecked, 2=include
	public void getAtomsMoreDistal(Molecule m, int resNum, Atom at2, int atomList[]) {
		
		if (at2.bond != null) {
			for(int q=0;q<at2.bond.length;q++) {
					Atom at = m.atom[at2.bond[q]];
					if (at.moleculeResidueNumber==resNum){
					if (atomList[at.residueAtomNumber] == 1) {
						atomList[at.residueAtomNumber] = 2;
						getAtomsMoreDistal(m,resNum,at,atomList);
					}
				}
			}
		}
	}


	// This function is called when NNew is aligned with NOld. We want
	//  the the axis and angle that will rotate CANew into CAOld
	// In the function call at1 = CAOld, at2 = NOld = NNew,
	//  and at3 = CANew
	// The axis is returned in axis and the angle is in thetaAngle
	// This is one of the more obfuscated functions
	private void getRotationInfoA(Atom at1, Atom at2, Atom at3,
		double thetaAngle[], double axis[]){
	
		// Theta is the angle between at1-at2-at3
		double	R2D = 57.29577951308232090712;
		double x12,y12,z12,x32,y32,z32,l12,l32,dp,dx,dy,dz;
		double fx,fy,fz,fp,thetaDeg2;
		double thetaDeg;

		x12 = at1.coord[0] - at2.coord[0];
		y12 = at1.coord[1] - at2.coord[1];
		z12 = at1.coord[2] - at2.coord[2];
		x32 = at3.coord[0] - at2.coord[0];
		y32 = at3.coord[1] - at2.coord[1];
		z32 = at3.coord[2] - at2.coord[2];		
		l12 = Math.sqrt(x12*x12 + y12*y12 + z12*z12);
		l32 = Math.sqrt(x32*x32 + y32*y32 + z32*z32);
		if ((l12 == 0.0) || (l32 == 0.0)) {
			thetaDeg = 0;
		}
		else {
			dp = (x12*x32 + y12*y32 + z12*z32) / (l12*l32);
			if (dp < -1.0)
				dp = -1.0;
			else if (dp > 1.0)
				dp = 1.0;
			thetaDeg = R2D * Math.acos(dp);
		}
		// To exactly pin down the angle we need to take the
		//  cross product of 12 and 32
		fx = y12*z32 - z12*y32;
		fy = z12*x32 - x12*z32;
		fz = x12*y32 - y12*x32;
		fp = (Math.sqrt(fx*fx + fy*fy + fz*fz))/(l12*l32);
		if (fp < -1.0)
		  fp = -1.0;
		else if (fp > 1.0)
		  fp = 1.0;
		thetaDeg2 = R2D * Math.asin(fp);
		// If the sign of the angle from the asin and acos do not agree
		//  then make the sign of the cos the same as that of the sin
		if (((thetaDeg2 > 0) && (thetaDeg < 0)) || 
				((thetaDeg2 < 0) && (thetaDeg > 0)))
			thetaDeg = -thetaDeg;

		// The rotation axis is the cross product of 12 cross 32
		// this is fx,fy,fz
			
		thetaAngle[0]=thetaDeg;
		axis[0]=fx;
		axis[1]=fy;
		axis[2]=fz;
	}

	// This function is called when CANew is already aligned with CAOld and NNew is
	//  aligned with NOld. We want the axis and angle that will rotate
	//  CBOld into CBNew
	// In the function call at1 = CBOld, at2 = CAOld = CANew,
	//  at3 = NOld = NNew, and at4 = CBNew
	// The axis is returned in axis and the angle is in thetaAngle
	// This is one of the more obfuscated functions
	private void getRotationInfoB(Atom at1, Atom at2, Atom at3, Atom at4,
		double thetaAngle[], double axis[]){

		double R2D = 57.29577951308232090712;
		double x12,y12,z12,x32,y32,z32,l12,l32,dp,dx,dy,dz;
		double x12b,y12b,z12b,x32b,y32b,z32b,ex,ey,ez;
		double fx,fy,fz,fp,thetaDeg2;
		double thetaDeg;

		x12 = at1.coord[0] - at2.coord[0];
		y12 = at1.coord[1] - at2.coord[1];
		z12 = at1.coord[2] - at2.coord[2];
		x32 = at3.coord[0] - at2.coord[0];
		y32 = at3.coord[1] - at2.coord[1];
		z32 = at3.coord[2] - at2.coord[2];		
		// d is cross product of vectors 12 and 32
		dx = y12*z32 - z12*y32;
		dy = z12*x32 - x12*z32;
		dz = x12*y32 - y12*x32;
		x12b = at4.coord[0] - at2.coord[0];
		y12b = at4.coord[1] - at2.coord[1];
		z12b = at4.coord[2] - at2.coord[2];
		x32b = at3.coord[0] - at2.coord[0];
		y32b = at3.coord[1] - at2.coord[1];
		z32b = at3.coord[2] - at2.coord[2];		
		// e is cross product of vectors 12b and 32b
		ex = y12b*z32b - z12b*y32b;
		ey = z12b*x32b - x12b*z32b;
		ez = x12b*y32b - y12b*x32b;

		// Now determine dot product between vectors d and e
		//  to get cos of angle
		l12 = Math.sqrt(dx*dx + dy*dy + dz*dz);
		l32 = Math.sqrt(ex*ex + ey*ey + ez*ez);
		if ((l12 == 0.0) || (l32 == 0.0)) {
			thetaDeg = 0;
		}
		else {
			dp = (dx*ex + dy*ey + dz*ez) / (l12*l32);
			if (dp < -1.0)
				dp = -1.0;
			else if (dp > 1.0)
				dp = 1.0;
			thetaDeg = R2D * Math.acos(dp);
		}
		// To exactly pin down the angle we need to take the
		//  cross product of d and e
		fx = dy*ez - dz*ey;
		fy = dz*ex - dx*ez;
		fz = dx*ey - dy*ex;
		fp = (Math.sqrt(fx*fx + fy*fy + fz*fz))/(l12*l32);

		if (fp < -1.0)
		  fp = -1.0;
		else if (fp > 1.0)
		  fp = 1.0;
		thetaDeg2 = R2D * Math.asin(fp);

		// If the sign of the angle from the asin and acos do not agree
		//  then make the sign of the cos the same as that of the sin
		if (((thetaDeg2 > 0) && (thetaDeg < 0)) ||
				((thetaDeg2 < 0) && (thetaDeg > 0)))
			thetaDeg = -thetaDeg;
		thetaAngle[0]=thetaDeg;
		
		if(debug)
			if ((Math.abs(fp) < 0.00001) && (Math.abs(thetaDeg-180) > 0.1))
				System.out.println("Warning: LovellRotamers: fp is small but rotation is not 180!");
				
		// if we're doing a 180 degree rotation then we can't
		//  get the axis from the cross product above as the
		//  cross product is zero, let the axis be the 32
		//  vector
		if (Math.abs(fp) < 0.00001) {
			System.out.println("Warning: GetRotationInfoB, crossproduct = 0");
			axis[0]=x32;
			axis[1]=y32;
			axis[2]=z32;
		}
		else {		 
			axis[0]=fx;
			axis[1]=fy;
			axis[2]=fz;
		}
	}

	// This function converts residue (resNum) to the conformation specified by rotamer rotNum
	public boolean applyRotamer(Molecule m, int resNum, int rotNum) {

		int aaNum = rl.getAARotamerIndex(m.strand[strandNumber].residue[resNum].name);
		// If this AA has no rotamers return false
		if (aaNum == -1)
			return false;
		// If the rotamer number is out of range return false
		if (rotNum >= rl.getNumRotForAAtype(aaNum))
			return false;

		int aaDihedrals = rl.getNumDihedrals(aaNum);
		Residue localResidue = m.strand[strandNumber].residue[resNum];
		int atomList[] = new int[localResidue.numberOfAtoms];
		int alSize = 0;

		for(int i=0;i<aaDihedrals;i++) {
	
			Atom at[] = new Atom[4];
			at[0] = at[1] = at[2] = at[3] = null;
			int atNum[] = new int[4];

			// Find atoms involved in the dihedral
			for(int q=0;q<localResidue.numberOfAtoms;q++) {
				for(int w=0;w<4;w++) {
					if (localResidue.atom[q].name.equalsIgnoreCase(rl.getDihedralAtomNames(aaNum,i,w))) {
						at[w] = localResidue.atom[q];
						atNum[w] = q;
					}
				}		
			}

			// If we didn't find all the atoms return false
			if ((at[0] == null) || (at[1] == null) || (at[2] == null) ||
				(at[3] == null))
				return false;

			for(int q=0;q<localResidue.numberOfAtoms;q++) {
				atomList[q]=1;
			}
			atomList[atNum[1]]=0;
			atomList[atNum[2]]=0;
				
			// Now find all atoms in the rotamer that are more distal than at[2]
			getAtomsMoreDistal(m,localResidue.moleculeResidueNumber,at[2],atomList);

			// Copy atoms over in atomList, ie. if current atomList[q]==2 then
			//  the new atomList[] should include q (ie. atom q counts)
			// Skip the 4th atom of the torsion as it's automatically rotated
			//  by m.setTorsion
			alSize=0;
			for(int q=0;q<localResidue.numberOfAtoms;q++){
				if ((atomList[q]==2) && (q != atNum[3]))
					atomList[alSize++]=localResidue.atom[q].moleculeAtomNumber;
			}

			// Perform the rotation
			// at[0], at[1], and at[2] don't move, at[3] and the atoms
			//  in the atomList rotate
			m.setTorsion(at[0].moleculeAtomNumber,at[1].moleculeAtomNumber,
				at[2].moleculeAtomNumber,at[3].moleculeAtomNumber,
				rl.getRotamerValues(aaNum,rotNum,i),atomList,alSize);
		}
	
		curRotNum[resNum] = rotNum;
		return true;
	}

	// This function changes the residue type (it performs a mutation)
	// The residue resNum of strand standNumber of molecule m is changed
	//  to a residue with name newResType
	// If addHydrogens is true then hydrogens are added to the new
	//  residue
	public void changeResidueType(Molecule m, int resNum, String newResType, boolean addHydrogens) {

		// call with connectResidue = true and useOldBBatoms = false
		changeResidueType(m,resNum,newResType,addHydrogens,true,false);
	}
	
	public void changeResidueType(Molecule m, int resNum, String newResType, boolean addHydrogens, boolean connectResidue) {

		// call with useOldBBatoms = false
		changeResidueType(m,resNum,newResType,addHydrogens,connectResidue,false);
	}

	// This function changes the residue type (it performs a mutation)
	// The residue resNum of strand strandNumber of molecule m is changed
	//  to a residue with name newResType
	// If addHydrogens is true then hydrogens are added to the new
	//  residue, otherwise hydrogens are stripped
	// If connectResidue is true then AA's are connected if their residue
	//  numbers are sequential.
	// If useOldBBatoms is true, then the *exact* coordinates for the old backbone
	//	CA, C, and CB atoms are copied for the new residue
	//
	// This function makes use of the CB carbon, unfortunately if we are changing
	//  to or from Gly we don't have a CB, in that case we match with HA3
	// So the overall scheme is:
	//  1) Translate N's to overlap
	//  2) Rotate CA's to overlap
	//  3) Rotate around N-CA axis to get CB's to overlap
	//
	// In step 3), we could align the backbone C's, but due to differences between 
	//		the PPR amino acid template and actual residue conformations, the CB may 
	//		differ significantly in such a case. To preserve the CA-CB orientation, 
	//		the CB's are aligned instead.
	//
	// This function is rather complicated due to the bookkeeping
	//  required to maintain our messy molecule datastructure.
	public void changeResidueType(Molecule m, int resNum, String newResType, boolean addHydrogens, boolean connectResidue, boolean useOldBBatoms) {

		Residue localResidue = m.strand[strandNumber].residue[resNum];

                RotMatrix rm = new RotMatrix();

		boolean newResGly = false;
		boolean oldResGly = false;
		boolean glyMutation = false;

                boolean newResPro = false;
                boolean oldResPro = false;
                boolean proMutation = false;

		/* KER: Histidine can now take on the forms HIP, HIE, and HID so it should
		 * be allowed to change between the different types
		if (isResHis(localResidue.name)&&isResHis(newResType)) // checks if the new and old types are both forms of histidine, if so a mutation is not done
			return;														// This prevents say HID from reverting into HIS
		*/

		// If the old or new residues are glycine then we must do special things as we treat the H as CB
                //Proline is also special
		if(newResType.equalsIgnoreCase("gly"))
                    newResGly = true;
                else if(newResType.equalsIgnoreCase("pro"))
                    newResPro = true;

		if(localResidue.name.equalsIgnoreCase("gly"))
                    oldResGly = true;
                else if(localResidue.name.equalsIgnoreCase("pro"))
                    oldResPro = true;

		if (oldResGly && newResGly)// Nothing to do here, a null mutation
			return;

		glyMutation = newResGly || oldResGly;
                proMutation = newResPro || oldResPro;

		// We assume a standard backbone ordering of: N,CA,C,O

		int savedMoleculeResidueNumber = localResidue.atom[0].moleculeResidueNumber;
		int savedStrandResidueNumber = localResidue.atom[0].strandResidueNumber;
		int savedStrandNumber = localResidue.atom[0].strandNumber;
		String savedSegID = localResidue.atom[0].segID;

		// Get the new residue from the templates
 		Amber96PolyPeptideResidue ppr = new Amber96PolyPeptideResidue();
		Residue r = ppr.getResidue(newResType);

		// Residue r = ppr.getResidue("Lala");
		Molecule m2 = new Molecule();
		m2.addResidue(0,r);
		m2.establishConnectivity(false);


		// First get a handle on the backbone N, CA, C, O, H, and CB atoms
		Atom at[] = null;
		at = getBBatoms(r); //for the new residue
		Atom NNew = at[0]; Atom CANew = at[1]; Atom CNew = at[2]; Atom ONew = at[3]; Atom HNew = at[4]; Atom CBNew = at[5];

		at = getBBatoms(localResidue); //for the old residue
		Atom NOld = at[0]; Atom CAOld = at[1]; Atom COld = at[2]; Atom OOld = at[3]; Atom HOld = at[4]; Atom CBOld = at[5];

		if (oldResGly){ // we didn't find CBOld as Gly doesn't have it; find HA3 and point CBOld to it
			for(int q=0;q<localResidue.numberOfAtoms;q++) {
				if ( (localResidue.atom[q].name.equalsIgnoreCase("HA3")) || (localResidue.atom[q].name.equalsIgnoreCase("3HA")) )
					CBOld = localResidue.atom[q];
			}
			if(CBOld == null){
				System.out.println("HA3 on "+localResidue.fullName+" couldn't be found. Please check labelling.");
				System.out.println("Going to assume we want HA2");
				for(int q=0;q<localResidue.numberOfAtoms;q++) {
					if ( (localResidue.atom[q].name.equalsIgnoreCase("HA2")) || (localResidue.atom[q].name.equalsIgnoreCase("2HA")) )
						CBOld = localResidue.atom[q];
                                }
			}
		}

                float newNHLength = rm.norm( rm.subtract( HNew.coord, NNew.coord ) );//New amide NH or N-CD bond length

			// START ALIGNMENT
		// Translate N's to overlap
		try{
                    float Ntrans[] = new float[3];
                    Ntrans[0] = NNew.coord[0] - NOld.coord[0];
                    Ntrans[1] = NNew.coord[1] - NOld.coord[1];
                    Ntrans[2] = NNew.coord[2] - NOld.coord[2];
                    for(int q=0;q<r.numberOfAtoms;q++) {
                            r.atom[q].coord[0] -= Ntrans[0];
                            r.atom[q].coord[1] -= Ntrans[1];
                            r.atom[q].coord[2] -= Ntrans[2];
                    }
                }
		catch(Exception E){
			System.out.println("1");
			E.printStackTrace();
		}

		int numAtoms = -1;
		int atomList[] = null;
		double thetaDeg[] = new double[1];
		double rotAxis[] = new double[3];
		numAtoms = r.numberOfAtoms;
		atomList = new int[numAtoms];

		// Rotate CAs to overlap
		if (CAOld.distance(CANew) > 0.0001) {
			getRotationInfoA(CAOld, NOld, CANew, thetaDeg, rotAxis);
			for(int q=0;q<r.numberOfAtoms;q++)
				atomList[q] = q;
			r.rotateResidue(NNew,rotAxis[0],rotAxis[1],rotAxis[2],-thetaDeg[0],atomList,numAtoms);
		}

		// Rotate CBs to overlap - now the residue backbones should be aligned
		if ( (!glyMutation) && (CBOld.distance(CBNew) > 0.0001) ) { //not a to- or from- Gly mutation
			getRotationInfoB(CBOld, CAOld, NOld, CBNew, thetaDeg, rotAxis);
			for(int q=0;q<r.numberOfAtoms;q++)
				atomList[q] = q;
			r.rotateResidue(CANew,rotAxis[0],rotAxis[1],rotAxis[2],-thetaDeg[0],atomList,numAtoms);
		}

		if (oldResGly) //mutation from Gly
			alignCBOldGly(CBOld,CBNew,CAOld,CANew,NOld,r);
		if (newResGly) //mutation to Gly
			alignCBNewGly(CBOld,CBNew,CAOld,CANew,NOld,r,localResidue);


		// Remove hydrogens if we don't want them
		if (!addHydrogens) {
			for(int q=0;q<r.numberOfAtoms;q++) {
				if (r.atom[q].elementType.equalsIgnoreCase("H"))
					m2.deleteAtom(r.atom[q].moleculeAtomNumber);
			}
		}
		else {
			// ELSE KEEP ALL HYDROGENS
		}

		//Make the positions of the new backbone atoms coincide *exactly* with the old ones;
		// 	NNew already has the same coordinates as NOld;
		// 	the CAs, CBs, and CBs are already aligned, but may differ due to differences between the template PPR and the residue backbone in the PDB
		if (useOldBBatoms){
			setAtomCoord(CANew,CAOld);
			if ( (CBNew!=null) && (CBOld!=null) ) //not a from/to Gly mutation (Gly does not have a CB)
				setAtomCoord(CBNew,CBOld);
		}
		setAtomCoord(CNew,COld);
		setAtomCoord(ONew,OOld);

                if( newResPro || ( oldResPro && addHydrogens ) ){

                    float NVec[] = rm.subtract(HOld.coord, NOld.coord);//N- to H or CD bond vector
                    NVec = rm.scale(NVec, newNHLength/rm.norm(NVec) );
                    HNew.coord = rm.add(NVec, NNew.coord);
                }
                else if(addHydrogens && (HOld != null))
			setAtomCoord(HNew,HOld);


		//////////////////////////////////////////////////////////////////////////////////
		if (localResidue.nterm){ //we are changing the nterm residue
			//we must add the H1, H2, H3 atoms (and delete the HN atom) to the new residue,
			//	since the PPR templates only handle polypeptide residues

			Residue r1 = getNtermRes(r,localResidue);

			m2 = new Molecule();
			m2.addResidue(0, r1);
			m2.establishConnectivity(false);
		}
		if (localResidue.cterm){ //we are changing the cterm residue
			//we must add the OXT or H1, H2, H3 atoms (and delete the HN atom) to the new residue,
			//	since the PPR templates only handle polypeptide residues

			Residue r1 = getCtermRes(r,localResidue);

			m2 = new Molecule();
			m2.addResidue(0, r1);
			m2.establishConnectivity(false);
		}
		////////////////////////////////////////////////////////////////////////////////


		// Copy the new residue information into the old residue
		int changeInAtoms = r.numberOfAtoms - localResidue.numberOfAtoms;
		int baseChangedAtom = localResidue.atom[0].moleculeAtomNumber;
		// first atomnum in next residue
		int nextResidueBase = -1;
		// the first atomnum in this residue
		int thisResidueBase = -1;
		// atomnumber of the C in the last residue
		int lastResidueC = -1;
		// Determine if the residue before and the one after are have sequential
		//  residue numbers. If they do and if we're interested in connecting
		//  sequential residues then gather information so we can make the
		//  appropriate bonds.
		boolean connectedResidue = false;
		boolean connectNextResidue = false;
		boolean connectLastResidue = false;
		thisResidueBase = m.strand[strandNumber].residue[resNum].atom[0].moleculeAtomNumber;
		nextResidueBase = thisResidueBase + localResidue.numberOfAtoms;

		if ((resNum+1) < m.strand[strandNumber].numberOfResidues) {
			if (connectResidue)
				connectNextResidue = (localResidue.getResNumber() + 1 == m.strand[strandNumber].residue[resNum+1].getResNumber());
		}
		if (resNum > 0) {
			if (connectResidue)
				connectLastResidue = (localResidue.getResNumber() - 1 == m.strand[strandNumber].residue[resNum-1].getResNumber());
			if (connectLastResidue) {
				for(int q=0;q<m.strand[strandNumber].residue[resNum-1].numberOfAtoms;q++) {
					if (m.strand[strandNumber].residue[resNum-1].atom[q].name.equalsIgnoreCase("C"))
						lastResidueC = m.strand[strandNumber].residue[resNum-1].atom[q].moleculeAtomNumber;
				}
			}
			if (lastResidueC == -1)
				connectLastResidue = false;
		}


		localResidue.name = r.name;

		if (r.name.length() > 3)
			localResidue.fullName = new String (r.name.substring(0,3) + " " + localResidue.fullName.substring(4));
		else if (localResidue.fullName.length()  > 4)
			localResidue.fullName = r.name.substring(0,3) + " " + localResidue.fullName.substring(4);
		else
			localResidue.fullName = r.name.substring(0,3);
		localResidue.numberOfAtoms = r.numberOfAtoms;
		localResidue.atom = r.atom;

		for(int j=0;j<localResidue.numberOfAtoms;j++){
			localResidue.atom[j].moleculeResidueNumber = savedMoleculeResidueNumber;
			localResidue.atom[j].strandResidueNumber = savedStrandResidueNumber;
			localResidue.atom[j].strandNumber = savedStrandNumber;
			localResidue.atom[j].segID = new String(savedSegID);
		}

		// Update atoms in residues of this strand after this residue
		//  as well as the bookkeeping in the molecule itself
		int curAtom = 0;
		int linkfrom=-1, linkto=-1;
		m.numberOfAtoms += changeInAtoms;
		m.numberOfAtomsx3 = m.numberOfAtoms * 3;
		m.atom = new Atom[m.numberOfAtoms];
		m.actualCoordinates = new float[m.numberOfAtomsx3];
		for(int j=0;j<m.numberOfStrands;j++) {
			for(int q=0;q<m.strand[j].numberOfResidues;q++) {
				for(int w=0;w<m.strand[j].residue[q].numberOfAtoms;w++) {
					m.strand[j].residue[q].atom[w].moleculeAtomNumber = curAtom;
					m.atom[curAtom++] = m.strand[j].residue[q].atom[w];
					int tmpIntAry[];
					tmpIntAry = m.strand[j].residue[q].atom[w].bond;
					if ((q==resNum) && (j==strandNumber)) {
						if (tmpIntAry != null) {
							for(int i=0;i<m.strand[j].residue[q].atom[w].bond.length;i++) {
								tmpIntAry[i] += baseChangedAtom;
							}
						}
						if (m.strand[j].residue[q].atom[w].name.equalsIgnoreCase("C") && (nextResidueBase != -1) && connectNextResidue) {
							// mark the atoms to bond so we can join when we're done
							//  we can't join them now as the second atom doesn't
							//  yet exist in the atom[] array
							connectedResidue = true;
							linkfrom = nextResidueBase+changeInAtoms;
							linkto = m.strand[j].residue[q].atom[w].moleculeAtomNumber;
						}
						if (m.strand[j].residue[q].atom[w].name.equalsIgnoreCase("N") && (lastResidueC != -1) && connectLastResidue) {
							m.addBondBetween(lastResidueC,m.strand[j].residue[q].atom[w].moleculeAtomNumber);
						}
					}
					else {
						if ((tmpIntAry != null) && (nextResidueBase != -1)) {
							for(int i=0;i<m.strand[j].residue[q].atom[w].bond.length;i++) {
								// Because the position of this residue's backbone C atom can move
								//  we have to be careful about simply updating the bond numbers.
								// The backbone N of numRes+1 if it bonds back to this residue needs
								//  special attention
								// Basically if there is a bond to the old residue we'll delete that
								//  bond. The only types of bonds into the new residue should be
								//  the previous and next peptide bond which we handle ourselves.
								if (tmpIntAry[i] >= nextResidueBase)
									tmpIntAry[i] += changeInAtoms;
								else if (tmpIntAry[i] >= thisResidueBase){
									m.strand[j].residue[q].atom[w].deleteBond(i);
									tmpIntAry = m.strand[j].residue[q].atom[w].bond;
									i--;
								}
							}
						}
					}

				}
			}
		}
		if (connectedResidue)
			m.addBondBetween(linkfrom,linkto);

		// Establish all connectivity including non-bonded interactions
		//KER: Remember that if you try to get rid of this establishConnectivity, you have to check for connectivity in Amber96ext.calculateTypesWithTemplates
		m.connectivityValid = false;
		//m.establishConnectivity(false);
		m.updateNumAtoms();

		curAAType[resNum] = newResType;
		curRotNum[resNum] = -1;
		localResidue.ffAssigned = false;

		// Copy atom coordinates back into actualCoordinates
		for(int q=0;q<m.numberOfAtoms;q++)
			m.updateCoordinates(q);


                if(proMutation){//Idealize the sidechain since proline has an unusual geometry
                    //This is especially important if we are mutating to proline (some bonds are probably way off length then:
                    //the idealization reconstructs the ideal ring given the backbone, CB, and CD coordinates)

                    m.idealizeResSidechain(localResidue);
                    int firstAtom = localResidue.atom[0].moleculeAtomNumber;
                    for(int a=0; a<localResidue.numberOfAtoms; a++)
                        m.resolveCoordinates(firstAtom+a);//Copy the idealized coordinates back into the Atom.coord arrays

                    if(!newResPro)//No longer proline, so ring closure is no longer an issue
                        localResidue.validConf = true;
                }

	}
	//Checks if residue res is a form of histidine
	/*private boolean isResHis(String res){
		
		if (res.equalsIgnoreCase("HID") || res.equalsIgnoreCase("HIE") || res.equalsIgnoreCase("HIP") || res.equalsIgnoreCase("HIS"))
			return true;
		else
			return false;
	}*/
	
	//Returns a handle on the backbone N, CA, C, O, H, and CB atoms for residue res
	private Atom [] getBBatoms(Residue res){
		
		Atom at[] = new Atom[6];
		
		for(int q=0;q<res.numberOfAtoms;q++) {
			if (res.atom[q].name.equalsIgnoreCase("N"))
				at[0] = res.atom[q];
			else if (res.atom[q].name.equalsIgnoreCase("CA"))
				at[1] = res.atom[q];
			else if (res.atom[q].name.equalsIgnoreCase("C"))
				at[2] = res.atom[q];
			else if (res.atom[q].name.equalsIgnoreCase("O")||res.atom[q].name.equalsIgnoreCase("O1"))
				at[3] = res.atom[q];
			else if ( res.atom[q].name.equalsIgnoreCase("H") || ( res.atom[q].name.equalsIgnoreCase("CD") && res.name.equalsIgnoreCase("PRO") ) )
				at[4] = res.atom[q];
			else if (res.atom[q].name.equalsIgnoreCase("CB"))
				at[5] = res.atom[q];
		}
		
		return at;
	}
	
	//Sets the coordinates of atom toAt to be the same as the coordinates of atom fromAt
	private void setAtomCoord(Atom toAt, Atom fromAt){
		toAt.coord[0] = fromAt.coord[0];
		toAt.coord[1] = fromAt.coord[1];
		toAt.coord[2] = fromAt.coord[2];
	}
	
	//This function returns a modified version of residue r by replacing the single backbone H atom
	//		with the H1, H2, and H3 atoms from localResidue, and updaing the bond information
	private Residue getNtermRes(final Residue r, final Residue localResidue){
		
		Residue r1 = r;
		
		//Delete the HN atom from the new residue
		for(int q=0;q<r1.numberOfAtoms;q++) {
			if (r1.atom[q].name.equalsIgnoreCase("H")){
				r1.deleteAtom(q);
			}
		}
		
		//Find the old H1, H2, H3 atoms and add them to the new residue
		for(int q=0;q<localResidue.numberOfAtoms;q++) {
			if ((localResidue.atom[q].name.equalsIgnoreCase("H1"))||
					(localResidue.atom[q].name.equalsIgnoreCase("H2"))||
					(localResidue.atom[q].name.equalsIgnoreCase("H3"))){
				
				Atom lAtom = localResidue.atom[q];					
				Atom at1 = new Atom(lAtom.name,lAtom.coord[0],lAtom.coord[1],lAtom.coord[2],lAtom.charge,lAtom.forceFieldType);					
				r1.addAtom(at1);
			}
		}
		
		//Add the new bonds between the H1, H2, and H3 atoms and the N atom
		for(int q=0;q<r1.numberOfAtoms;q++) {
			if ((r1.atom[q].name.equalsIgnoreCase("H1"))||(r.atom[q].name.equalsIgnoreCase("H2"))||
					(r1.atom[q].name.equalsIgnoreCase("H3"))){
				
				r1.atom[q].addBond(0); //add the bond for this H atom to the N atom
				r1.atom[0].addBond(q); //add the bond for the N atom to this H atom
			}
		}		
		
		return r1;
	}
	
	private Residue getCtermRes(final Residue r, final Residue localResidue){
		
		Residue r1 = r;
		
		
		//Find the old H1, H2, H3 atoms and add them to the new residue
		for(int q=0;q<localResidue.numberOfAtoms;q++) {
			if (localResidue.atom[q].name.equalsIgnoreCase("OXT")||
					localResidue.atom[q].name.equalsIgnoreCase("O2")){
				
				Atom lAtom = localResidue.atom[q];					
				Atom at1 = new Atom(lAtom.name,lAtom.coord[0],lAtom.coord[1],lAtom.coord[2],lAtom.charge,lAtom.forceFieldType);					
				r1.addAtom(at1);
			}
		}
		
		int cAtom = -1;
		//find which atom is the C to connect the OXT to
		for(int q=0; q<r1.numberOfAtoms;q++){
			if(r1.atom[q].name.equalsIgnoreCase("C")){
				cAtom = q;
				break;
			}
		}
		
		//update the bond information for the new atom
		for(int q=0; q<r1.numberOfAtoms;q++){
			if (r1.atom[q].name.equalsIgnoreCase("OXT")||r1.atom[q].name.equalsIgnoreCase("O2")){
				r1.atom[q].addBond(cAtom); //add the bond to the new C atom
				r1.atom[cAtom].addBond(q);
			}
		}		
		
		return r1;
	}
	
	//Performs the CB alignment (changes the coordinates of the new residue r) for the new residue when the old residue is Gly
	private void alignCBOldGly(Atom CBOld, Atom CBNew, Atom CAOld, Atom CANew, Atom NOld, Residue r){
		
	  try{
		int numAtoms = -1;
		int atomList[] = null;
		double thetaDeg[] = new double[1];
		double rotAxis[] = new double[3];
		numAtoms = r.numberOfAtoms;
		atomList = new int[numAtoms];
		
		// Compute the angle between CBOld-CA and CBNew-CA; here, CBOld is actually HA3
		double magOld = CBOld.distance(CAOld);
		double magNew = CBNew.distance(CAOld);
		double dotProd = ((CBNew.coord[0] - CAOld.coord[0]) * (CBOld.coord[0] - CAOld.coord[0])) +
			((CBNew.coord[1] - CAOld.coord[1]) * (CBOld.coord[1] - CAOld.coord[1])) +
			((CBNew.coord[2] - CAOld.coord[2]) * (CBOld.coord[2] - CAOld.coord[2]));
		double angle = 180.0 / 3.1415 * Math.acos(dotProd / magOld / magNew);
		if (Math.abs(angle) > 0.1){
			getRotationInfoB(CBOld, CAOld, NOld, CBNew, thetaDeg, rotAxis);
			for(int q=0;q<r.numberOfAtoms;q++)
				atomList[q] = q;
			r.rotateResidue(CANew,rotAxis[0],rotAxis[1],rotAxis[2],-thetaDeg[0],atomList,numAtoms);
		}
	  }
	  catch(Exception e){
		  System.out.println("Why are you failing?");
	  }
	}
	
	//Performs the CB alignment (changes the coordinates of the new residue r) when the new residue is Gly
	private void alignCBNewGly(Atom CBOld, Atom CBNew, Atom CAOld, Atom CANew, Atom NOld, Residue r, Residue localResidue){
		
		//If the new residue is a Gly, we snap the HA2 and HA3 atoms
		
		// Snap the HA2
		int localH = -1;
		for(int q=0;q<localResidue.numberOfAtoms;q++){
			if(localResidue.atom[q].name.equalsIgnoreCase("HA"))
				localH = q;
		}
		for(int q=0;q<r.numberOfAtoms;q++){
			if( r.atom[q].name.equalsIgnoreCase("HA2") || r.atom[q].name.equalsIgnoreCase("2HA") ){
				r.atom[q].coord[0] = localResidue.atom[localH].coord[0];
				r.atom[q].coord[1] = localResidue.atom[localH].coord[1];
				r.atom[q].coord[2] = localResidue.atom[localH].coord[2];
			}
		}
		
		// Snap the HA3 (but note that we can't just copy coordinates
		//  because the C-H bond length is different than the C-C
		//  bond length, so we have to scale
		localH = -1;
		for(int q=0;q<r.numberOfAtoms;q++){
			if( r.atom[q].name.equalsIgnoreCase("HA3") || r.atom[q].name.equalsIgnoreCase("3HA") )
				localH = q;
		}
		Atom tempH = r.atom[localH];
		double magNew = Math.sqrt(((tempH.coord[0] - CANew.coord[0]) * (tempH.coord[0] - CANew.coord[0])) +
			((tempH.coord[1] - CANew.coord[1]) * (tempH.coord[1] - CANew.coord[1])) +
			((tempH.coord[2] - CANew.coord[2]) * (tempH.coord[2] - CANew.coord[2])));
		double magOld = Math.sqrt(((CBOld.coord[0] - CANew.coord[0]) * (CBOld.coord[0] - CANew.coord[0])) +
			((CBOld.coord[1] - CANew.coord[1]) * (CBOld.coord[1] - CANew.coord[1])) +
			((CBOld.coord[2] - CANew.coord[2]) * (CBOld.coord[2] - CANew.coord[2])));
		r.atom[localH].coord[0] = (float)(CANew.coord[0]+(CBOld.coord[0]-CANew.coord[0])*magNew/magOld);
		r.atom[localH].coord[1] = (float)(CANew.coord[1]+(CBOld.coord[1]-CANew.coord[1])*magNew/magOld);
		r.atom[localH].coord[2] = (float)(CANew.coord[2]+(CBOld.coord[2]-CANew.coord[2])*magNew/magOld);
	}
	
	public void addOrigRots(int[][] strandMut,RotamerLibrary rl, Molecule m){
		
		if(! rl.isAddedRotamers()){
			for(int str=0; str<strandMut.length;str++){
				for(int res=0; res<strandMut[str].length;res++){
					Residue curRes = m.strand[str].residue[strandMut[str][res]];
					//Get Num Dihedrals
					int numDiheds = rl.getNumDihedrals(rl.getAARotamerIndex(curRes.name));
					if(numDiheds<=0)
						continue;
					int[] diheds = new int[numDiheds]; 
					int atoms[] = new int[4];
					//get all dihedrals
					for(int i=0; i<numDiheds; i++){
						atoms = rl.getDihedralInfo(m, curRes.strandNumber, curRes.strandResidueNumber, i);
						diheds[i] = (int) Math.round(curRes.atom[atoms[3]].torsion(curRes.atom[atoms[0]], curRes.atom[atoms[1]], curRes.atom[atoms[2]]));
					}
					//System.out.print("Str: "+str+" Res: "+res+" ");
					rl.addRotamer(curRes.name, diheds);
					//System.out.println("");
				}
			
			}
			rl.setAddedRotamers(true);
		}
		else{
			System.out.println("DEBUG: ALREADY ADDED ROTAMERS");
		}
	}
	
}
