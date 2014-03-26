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
// Atom.java
//
//  Version:           2.1 beta
//
//
// authors:
//    initials    name                 organization                email
//   ---------   -----------------    ------------------------    ----------------------------
//     RHL        Ryan Lilien          Dartmouth College           ryan.lilien@dartmouth.edu
//	   ISG		  Ivelin Georgiev	   Duke University			   ivelin.georgiev@duke.edu
//     KER        Kyle E. Roberts       Duke University         ker17@duke.edu
//     PGC        Pablo Gainza C.       Duke University         pablo.gainza@duke.edu
///////////////////////////////////////////////////////////////////////////////////////////////

/*
 *
 * Major changes were made by Ryan Lilien (2001-2004)
 * Minor changes were made by Ivelin Georgiev (2004-2009)
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
 * Handles functions and data associated with atoms. Example functions include adding a bond between atoms,
 * computing torsion, computing atom distance, etc. Some of the data members include the atom name, radius,
 * mass, coordinates, and bond information.
 */
public class Atom implements Serializable {

	int	moleculeAtomNumber=-1;		// the atom number in this molecule
	int	residueAtomNumber=-1;			// the atom number within this residue
	int	modelAtomNumber=-1;				// a temp variable used on writing to file
	int	moleculeResidueNumber=-1;		// the number of the residue containing this atom
	int	strandResidueNumber=-1;		// the strand relative number of the residue continaing atom
	int	strandNumber=-1;					// the number of the strand containing atom
	int	bond[];		// arrays of indices for bonded atoms
	int	elementNumber=0;					// element number
	int	numberOfBonds=0;					// the number of bonds
	String elementType="";		// the element symbol
	String forceFieldType="";	// the force field type symbol
	int type=-1;							// the force field type number
	boolean	selected = true;
	String	name;									// the atom name ie. CA, CB, CG, ...
	String  segID = "";			// the atom segment id, specifically included for the molecular replacement proj
	float	charge=0.0f;	// atomic charge
	int	radius = 170;		// vdw radii in pm (100pm=1A)
	double	mass;				// atomic mass
	float	coord[] = new float[3];		// atomic coordinates
	boolean isBBatom = false;
        float BFactor;


	Atom(){
	}

	Atom(String atomName, float xpos, float ypos, float zpos){
		name = atomName;
		setProperties();
		coord[0] = xpos;
		coord[1] = ypos;
		coord[2] = zpos;
		
		isBBatom = setIsBBatom();
	}

	Atom(String atomName, float xpos, float ypos, float zpos, String ffType){
		name = atomName;
		setProperties();
		coord[0] = xpos;
		coord[1] = ypos;
		coord[2] = zpos;
		forceFieldType = ffType;
		
		isBBatom = setIsBBatom();
	}

	Atom(String atomName, float xpos, float ypos, float zpos, float thecharge){
		name = atomName;
		setProperties();
		coord[0] = xpos;
		coord[1] = ypos;
		coord[2] = zpos;
		charge = thecharge;
		
		isBBatom = setIsBBatom();
	}
	
	Atom(String atomName, float xpos, float ypos, float zpos, float thecharge, String ffType){	
		name = atomName;
		setProperties();
		coord[0] = xpos;
		coord[1] = ypos;
		coord[2] = zpos;
		forceFieldType = ffType;
		charge = thecharge;
		
		isBBatom = setIsBBatom();
	}

        Atom(float coor[]){
            coord = coor;
        }

	// Changes the atom name and sets the atom properties element number,
	//  element type, radii, and mass
	public void changeType(String newName){
		name = newName;	
		setProperties();
	}

	// Sets default properties such as radii, element number, mass, element type
	// RHL I don't know where the radii values come from but they are good estimates
	//  usually atoms will be assigned better information from a forcefield
	private void setProperties(){
		forceFieldType = "?";
		if ( name.indexOf("Ac") != -1 ){
			elementNumber = 89;
			radius = 295;
			mass = 227;
			elementType = "Ac";
		}
		else if( name.indexOf("Ag") != -1 ){
			radius = 398;
			elementNumber = 47;
			mass = 107.9;
			elementType = "Ag";
		}
		else if( name.indexOf("Al") != -1 ){
			radius = 338;
			elementNumber = 13;
			mass = 26.98;
			elementType = "Al";
		}
		else if( name.indexOf("Am") != -1 ){
			elementNumber = 95;
			radius = 230;
			mass = 243;
			elementType = "Am";
		}
		else if( name.indexOf("Ar") != -1 ){
			elementNumber = 18;
			radius = 392;
			mass = 39.95;
			elementType = "Ar";
		}
		else if( name.indexOf("As") != -1 ){
			elementNumber = 33;
			radius = 302;
			mass = 74.92;
			elementType = "As";
		}
		else if( name.indexOf("At") != -1 ){
			elementNumber = 85;
			radius = 302;
			mass = 210;
			elementType = "At";
		}
		else if( name.indexOf("Au") != -1 ){
			radius = 375;
			elementNumber = 79;
			mass = 197;
			elementType = "Au";
		}
		else if( name.indexOf("Ba") != -1 ){
			radius = 335;
			elementNumber = 56;
			mass = 137.3;
			elementType = "Ba";
		}
		else if( name.indexOf("Be") != -1 ){
			elementNumber = 4;
			radius = 88;
			mass = 9.012;
			elementType = "Be";
		}
		else if( name.indexOf("Bi") != -1 ){
			elementNumber = 83;
			radius = 385;
			mass = 209;
			elementType = "Bi";
		}
		else if( name.indexOf("Bk") != -1 ){
			elementNumber = 97;
			radius = 225;
			mass = 247;
			elementType = "Bk";
		}
		else if( name.indexOf("Br") != -1 ){
			radius = 302;
			elementNumber = 97;
			mass = 247;
			elementType = "Br";
		}
		else if( name.indexOf("Ca") != -1 ){
			radius = 248;
			elementNumber = 20;
			mass = 40.08;
			elementType = "Ca";
		}
		else if( name.indexOf("Cd") != -1 ){
			elementNumber = 48;
			radius = 422;
			mass = 112.4;
			elementType = "Cd";
		}
		else if( name.indexOf("Ce") != -1 ){
			elementNumber = 58;
			radius = 458;
			mass = 140.1;
			elementType = "Ce";
		}
		else if( name.indexOf("Cf") != -1 ){
			elementNumber = 98;
			radius = 222;
			mass = 249;
			elementType = "Cf";
		}
		else if( name.indexOf("Cl") != -1 ){
			elementNumber = 17;
			radius = 250;
			mass = 35.45;
			elementType = "Cl";
		}
		else if( name.indexOf("Cm") != -1 ){
			elementNumber = 96;
			radius = 228;
			mass = 247;
			elementType = "Cm";
		}
		else if( name.indexOf("Co") != -1 ){
			radius = 332;
			elementNumber = 27;
			mass = 58.93;
			elementType = "Co";
		}
		else if( name.indexOf("Cr") != -1 ){
			radius = 338;
			elementNumber = 24;
			mass = 52;
			elementType = "Cr";
		}
		else if( name.indexOf("Cs") != -1 ){
			elementNumber = 55;
			radius = 418;
			mass = 132.9;
			elementType = "Cs";
		}
		else if( name.indexOf("Cu") != -1 ){
			radius = 380;
			elementNumber = 29;
			mass = 63.55;
			elementType = "Cu";
		}
		else if( name.indexOf("Dy") != -1 ){
			elementNumber = 66;
			radius = 438;
			mass = 162.5;
			elementType = "Dy";
		}
		else if( name.indexOf("Er") != -1 ){
			elementNumber = 68;
			radius = 432;
			mass = 167.3;
			elementType = "Er";
		}
		else if( name.indexOf("Es") != -1 ){
			elementNumber = 99;
			radius = 220;
			mass = 254;
			elementType = "Es";
		}
		else if( name.indexOf("Eu") != -1 ){
			elementNumber = 63;
			radius = 498;
			mass = 152;
			elementType = "Eu";
		}
		else if( name.indexOf("Fe") != -1 ){
			radius = 335;
			elementNumber = 26;
			mass = 55.85;
			elementType = "Fe";
		}
		else if( name.indexOf("Fm") != -1 ){
			elementNumber = 100;
			radius = 218;
			mass = 250;
			elementType = "Fm";
		}
		else if( name.indexOf("Fr") != -1 ){
			elementNumber = 87;
			radius = 450;
			mass = 223;
			elementType = "Fr";
		}
		else if( name.indexOf("Ga") != -1 ){
			radius = 305;
			elementNumber = 31;
			mass = 69.72;
			elementType = "Ga";
		}
		else if( name.indexOf("Gd") != -1 ){
			elementNumber = 64;
			radius = 448;
			mass = 157.3;
			elementType = "Gd";
		}
		else if( name.indexOf("Ge") != -1 ){
			elementNumber = 32;
			radius = 292;
			mass = 72.59;
			elementType = "Ge";
		}
		else if( name.indexOf("He") != -1 ){
			radius = 400;
			elementNumber = 2;
			mass = 4.003;
			elementType = "He";
		}
		else if( name.indexOf("Hf") != -1 ){
			elementNumber = 72;
			radius = 392;
			mass = 178.5;
			elementType = "Hf";
		}
		else if( name.indexOf("Hg") != -1 ){
			elementNumber = 80;
			radius = 425;
			mass = 200.6;
			elementType = "Hg";
		}
		else if( name.indexOf("Ho") != -1 ){
			elementNumber = 67;
			radius = 435;
			mass = 164.9;
			elementType = "Ho";
		}
		else if( name.indexOf("In") != -1 ){
			elementNumber = 49;
			radius = 408;
			mass = 114.8;
			elementType = "In";
		}
		else if( name.indexOf("Ir") != -1 ){
			elementNumber = 77;
			radius = 330;
			mass = 192.2;
			elementType = "Ir";
		}
		else if( name.indexOf("Kr") != -1 ){
			elementNumber = 36;
			radius = 400;
			mass = 83.8;
			elementType = "Kr";
		}
		else if( name.indexOf("La") != -1 ){
			elementNumber = 57;
			radius = 468;
			mass = 138.9;
			elementType = "La";
		}
		else if( name.indexOf("Li") != -1 ){
			radius = 170;
			elementNumber = 3;
			mass = 6.941;
			elementType = "Li";
		}
		else if( name.indexOf("Lr") != -1 ){
			elementNumber = 103;
			radius = 210;
			mass = 257;
			elementType = "Lr";
		}
		else if( name.indexOf("Lu") != -1 ){
			elementNumber = 71;
			radius = 430;
			mass = 175;
			elementType = "Lu";
		}
		else if( name.indexOf("Md") != -1 ){
			elementNumber = 101;
			radius = 215;
			mass = 256;
			elementType = "Md";
		}
		else if( name.indexOf("Mg") != -1 ){
			radius = 275;
			elementNumber = 12;
			mass = 24.31;
			elementType = "Mg";
		}
		else if( name.indexOf("Mn") != -1 ){
			radius = 338;
			elementNumber = 25;
			mass = 54.94;
			elementType = "Mn";
		}
		else if( name.indexOf("Mo") != -1 ){
			elementNumber = 42;
			radius = 368;
			mass = 95.94;
			elementType = "Mo";
		}
		else if( name.indexOf("Na") != -1 ){
			elementNumber = 11;
			radius = 243;
			mass = 22.99;
			elementType = "Na";
		}
		else if( name.indexOf("Nb") != -1 ){
			elementNumber = 41;
			radius = 370;
			mass = 92.91;
			elementType = "Nb";
		}
		else if( name.indexOf("Nd") != -1 ){
			elementNumber = 60;
			radius = 452;
			mass = 144.2;
			elementType = "Nd";
		}
		else if( name.indexOf("Ne") != -1 ){
			elementNumber = 10;
			radius = 280;
			mass = 20.18;
			elementType = "Ne";
		}
		else if( name.indexOf("Ni") != -1 ){
			radius = 405;
			elementNumber = 28;
			mass = 58.69;
			elementType = "Ni";
		}
		else if( name.indexOf("No") != -1 ){
			elementNumber = 102;
			radius = 212;
			mass = 253;
			elementType = "No";
		}
		else if( name.indexOf("Np") != -1 ){
			elementNumber = 93;
			radius = 238;
			mass = 237;
			elementType = "Np";
		}
		else if( name.indexOf("Os") != -1 ){
			elementNumber = 76;
			radius = 342;
			mass = 190.2;
			elementType = "Os";
		}
		else if( name.indexOf("Pa") != -1 ){
			elementNumber = 91;
			radius = 222;
			mass = 231;
			elementType = "Pa";
		}
		else if( name.indexOf("Pb") != -1 ){
			elementNumber = 82;
			radius = 385;
			mass = 207.2;
			elementType = "Pb";
		}
		else if( name.indexOf("Pd") != -1 ){
			elementNumber = 46;
			radius = 375;
			mass = 106.4;
			elementType = "Pd";
		}
		else if( name.indexOf("Pm") != -1 ){
			elementNumber = 61;
			radius = 450;
			mass = 147;
			elementType = "Pm";
		}
		else if( name.indexOf("Po") != -1 ){
			elementNumber = 84;
			radius = 420;
			mass = 210;
			elementType = "Po";
		}
		else if( name.indexOf("Pr") != -1 ){
			elementNumber = 59;
			radius = 455;
			mass = 140.9;
			elementType = "Pr";
		}
		else if( name.indexOf("Pt") != -1 ){
			elementNumber = 78;
			radius = 375;
			mass = 195.1;
			elementType = "Pt";
		}
		else if( name.indexOf("Pu") != -1 ){
			elementNumber = 94;
			radius = 232;
			mass = 242;
			elementType = "Pu";
		}
		else if( name.indexOf("Ra") != -1 ){
			elementNumber = 88;
			radius = 358;
			mass = 226;
			elementType = "Ra";
		}
		else if( name.indexOf("Rb") != -1 ){
			elementNumber = 37;
			radius = 368;
			mass = 85.47;
			elementType = "Rb";
		}
		else if( name.indexOf("Re") != -1 ){
			elementNumber = 75;
			radius = 338;
			mass = 186.2;
			elementType = "Re";
		}
		else if( name.indexOf("Rh") != -1 ){
			elementNumber = 45;
			radius = 362;
			mass = 102.9;
			elementType = "Rh";
		}
		else if( name.indexOf("Rn") != -1 ){
			elementNumber = 88;
			radius = 475;
			mass = 222;
			elementType = "Rn";
		}
		else if( name.indexOf("Ru") != -1 ){
			elementNumber = 44;
			radius = 350;
			mass = 101.1;
			elementType = "Ru";
		}
		else if( name.indexOf("Sb") != -1 ){
			elementNumber = 51;
			radius = 365;
			mass = 121.8;
			elementType = "Sb";
		}
		else if( name.indexOf("Sc") != -1 ){
			radius = 360;
			elementNumber = 21;
			mass = 44.96;
			elementType = "Sc";
		}
		else if( name.indexOf("Se") != -1 ){
			elementNumber = 34;
			radius = 305;
			mass = 78.96;
			elementType = "Se";
		}
		else if( name.indexOf("Si") != -1 ){
			radius = 300;
			elementNumber = 14;
			mass = 28.09;
			elementType = "Si";
		}
		else if( name.indexOf("Sm") != -1 ){
			elementNumber = 62;
			radius = 450;
			mass = 150.4;
			elementType = "Sm";
		}
		else if( name.indexOf("Sn") != -1 ){
			elementNumber = 50;
			radius = 365;
			mass = 118.7;
			elementType = "Sn";
		}
		else if( name.indexOf("Sr") != -1 ){
			elementNumber = 38;
			radius = 280;
			mass = 87.62;
			elementType = "Sr";
		}
		else if( name.indexOf("Ta") != -1 ){
			elementNumber = 73;
			radius = 358;
			mass = 180.9;
			elementType = "Ta";
		}
		else if( name.indexOf("Tb") != -1 ){
			elementNumber = 65;
			radius = 440;
			mass = 158.9;
			elementType = "Tb";
		}
		else if( name.indexOf("Tc") != -1 ){
			elementNumber = 43;
			radius = 338;
			mass = 99;
			elementType = "Tc";
		}
		else if( name.indexOf("Te") != -1 ){
			elementNumber = 52;
			radius = 368;
			mass = 127.6;
			elementType = "Te";
		}
		else if( name.indexOf("Th") != -1 ){
			elementNumber = 90;
			radius = 255;
			mass = 232;
			elementType = "Th";
		}
		else if( name.indexOf("Ti") != -1 ){
			radius = 368;
			elementNumber = 22;
			mass = 47.88;
			elementType = "Ti";
		}
		else if( name.indexOf("Tl") != -1 ){
			elementNumber = 81;
			radius = 388;
			mass = 204.4;
			elementType = "Tl";
		}
		else if( name.indexOf("Tm") != -1 ){
			elementNumber = 69;
			radius = 430;
			mass = 168.9;
			elementType = "Tm";
		}
		else if( name.indexOf("Une") != -1 ){
			elementNumber = 109;
			mass = 266;
			elementType = "Une";
		}
		else if( name.indexOf("Unh") != -1 ){
			elementNumber = 106;
			mass = 263;
			elementType = "Unh";
		}
		else if( name.indexOf("Uno") != -1 ){
			elementNumber = 108;
			mass = 265;
			elementType = "Uno";
		}
		else if( name.indexOf("Unp") != -1 ){
			elementNumber = 105;
			mass = 260;
			elementType = "Unp";
		}
		else if( name.indexOf("Unq") != -1 ){
			elementNumber = 104;
			mass = 257;
			elementType = "Unq";
		}
		else if( name.indexOf("Uns") != -1 ){
			elementNumber = 107;
			mass = 262;
			elementType = "Uns";
		}
		else if( name.indexOf("Xe") != -1 ){
			elementNumber = 54;
			radius = 425;
			mass = 131.3;
			elementType = "Xe";
		}
		else if( name.indexOf("Yb") != -1 ){
			elementNumber = 70;
			radius = 485;
			mass = 173;
			elementType = "Yb";
		}
		else if( name.indexOf("Zn") != -1 ){
			radius = 362;
			elementNumber = 30;
			mass = 65.39;
			elementType = "Zn";
		}
		else if( name.indexOf("Zr") != -1 ){
			elementNumber = 40;
			radius = 390;
			mass = 91.22;
			elementType = "Zr";
		}
		//The radii for C, H, N, O, P, S were changed to those used by the Richardsons' PROBE;
		//		The old radii are commented out to the side
		else if( name.toUpperCase().indexOf("C") == 0 ){
			radius = 165; //radius = 180;
			elementNumber = 6;
			mass = 12.01;
			elementType = "C";
		}
		else if( (name.toUpperCase().indexOf("H") == 0) || 
				( (name.toUpperCase().indexOf("H") == 1) && ((name.charAt(0)>='0') && (name.charAt(0)<='9')) ) ){
			radius = 100; //radius = 80;
			elementNumber = 1;
			mass = 1;
			elementType = "H";
		}
		else if( name.toUpperCase().indexOf("N") == 0 ){
			radius = 155; //radius = 170;
			elementNumber = 7;
			mass = 14.01;
			elementType = "N";
		}
		else if( name.toUpperCase().indexOf("O") == 0 ){
			radius = 140; //radius = 170;
			elementNumber = 8;
			mass = 16;
			elementType = "O";
		}
		else if( name.toUpperCase().indexOf("B") == 0 ){
			radius = 208;
			elementNumber = 5;
			mass = 10.81;
			elementType = "B";
		}
		else if( name.toUpperCase().indexOf("I") == 0 ){
			radius = 350;
			elementNumber = 53;
			mass = 126.9;
			elementType = "I";
		}
		else if( name.toUpperCase().indexOf("F") == 0 ){
			radius = 160;
			elementNumber = 9;
			mass = 19.0;
			elementType = "F";
		}
		else if( name.toUpperCase().indexOf("P") == 0 ){
			radius = 180; //radius = 259;
			elementNumber = 15;
			mass = 30.97;
			elementType = "P";
		}
		else if( name.toUpperCase().indexOf("K") == 0 ){
			radius = 332;
			elementNumber = 19;
			mass = 39.1;
			elementType = "K";
		}
		else if( name.toUpperCase().indexOf("S") == 0 ){
			radius = 180; //radius = 255;
			elementNumber = 16;
			mass = 32.07;
			elementType = "S";
		}
		else if( name.toUpperCase().indexOf("U") == 0 ){
			radius = 242;
			elementNumber = 92;
			mass = 238;
			elementType = "U";
		}
		else if( name.toUpperCase().indexOf("V") == 0 ){
			radius = 332;
			elementNumber = 23;
			mass = 50.94;
			elementType = "V";
		}
		else if( name.toUpperCase().indexOf("W") == 0 ){
			radius = 342;
			elementNumber = 74;
			mass = 183.9;
			elementType = "W";
		}
		else if( name.toUpperCase().indexOf("Y") == 0 ){
			radius = 445;
			elementNumber = 39;
			mass = 88.91;
			elementType = "Y";
		}
		else {  //unrecognized atom type
			radius = 10000;
			elementType = "DU"; //dummy atom
		}
	}

	// Adds a bond from this atom to atom number j
	// If we're in a molecule then j should be molecule based
	public void addBond(int j){
		if (bond==null){
			bond = new int[++numberOfBonds];
			bond[0] = j;
		}
		else{
			// If we already have a bond to this atom we're done
			if (bondedTo(j))
				return;
			int newBond[] = new int[++numberOfBonds];
			System.arraycopy(bond,0,newBond,0,bond.length);
			bond = newBond;
			bond[numberOfBonds-1] = j;
		}
	}

	// Deletes the jth bond of this atom
	// *NOTE j is not an atom number, it is an index into the
	//  bond array
	public void deleteBond(int j){
		if (numberOfBonds == 0)
			return;
		int smallerBondArray[] = new int[numberOfBonds-1];
		if (numberOfBonds>1){
			System.arraycopy(bond, 0, smallerBondArray, 0, j);
			if (j<numberOfBonds-1)
				System.arraycopy(bond,j+1,smallerBondArray,j,bond.length-j-1);
		}
		bond = smallerBondArray;
		numberOfBonds--;
	}

	// Deletes the bond to atom number atomj
	// If we're in a molecule atomj should be molecule based
	public void deleteBondTo(int atomj){
		int j;
		for(j=0; j<numberOfBonds; j++)
			if (bond[j] == atomj){
				deleteBond(j);
				j--;
			}
	}

	// Adds two bonds, one each to atoms j and k
	// Nothing fancy, this just makes it easier to add multiple bonds
	public void addBond(int j, int k){
		addBond(j);
		addBond(k);		
 	}

	// Adds two bonds, one each to atoms j, k, and l
	// Nothing fancy, this just makes it easier to add multiple bonds
	public void addBond(int j, int k, int l){
		addBond(j);
		addBond(k);		
		addBond(l);
	}

	// Adds four bonds, one each to atoms j, k, l, and m
	// Nothing fancy, this just makes it easier to add multiple bonds
	public void addBond(int j, int k, int l, int m){
		addBond(j);
		addBond(k);		
		addBond(l);
		addBond(m);
	}

	// Returns true/false, is this atom bonded to atom j
	public boolean bondedTo(int j){
		for(int i=0; i<numberOfBonds; i++){
			if (bond[i] == j)
				return true;
		}
		return false;
	}

	// Returns the distance from this atom to the specified atom
	public double distance(Atom atom){
		return (Math.sqrt((coord[0] - atom.coord[0]) * (coord[0] - atom.coord[0]) + 
			(coord[1] - atom.coord[1]) * (coord[1] - atom.coord[1]) +
			(coord[2] - atom.coord[2]) * (coord[2] - atom.coord[2])));
	}

	// This returns the right handed rotation angle around the +x
	//  axis where the -z axis is the 0. (0..360)
	// This was not written by me, but I have checked it
	public double angleAboutXAxis(){
		double	R2D = 57.29577951308232090712;
		double distance = Math.sqrt( coord[ 0 ] * coord[ 0 ] +
			coord[ 1 ] * coord[ 1 ] + coord[ 2 ] * coord[ 2 ] );
		double yComponent = coord[ 1 ] / distance;
		double zComponent = coord[ 2 ] / distance;
		double theta = 0.0;
		if ( ( zComponent == 0 ) && ( yComponent > 0 ) )
			theta = 90.0;
		else if ( ( zComponent == 0 ) && ( yComponent < 0 ) )
			theta = 270.0;
		else if ( ( yComponent == 0 ) && ( zComponent > 0 ) )
			theta = 180.0;
		else if ( ( yComponent == 0 ) && ( zComponent <= 0 ) )
			theta = 0.0;
		else{ 
			theta = Math.atan( Math.abs( yComponent / zComponent ) ) * R2D;
			if ( ( yComponent > 0 ) && ( zComponent > 0 ) )
				theta = 90.0 + Math.atan( Math.abs( zComponent / 
					yComponent ) ) * R2D;
			else if ( ( yComponent < 0 ) && ( zComponent > 0 ) )
				theta += 180.0;
			else if ( ( yComponent < 0 ) && ( zComponent < 0 ) )
				theta = 270.0 + Math.atan( Math.abs( zComponent / 
					yComponent ) ) * R2D;
		}
		return( theta );
	}

	// This returns the right handed rotation angle around the +y
	//  axis where the -x axis is the 0
	// This was not written by me, but I have checked it
	public double angleAboutYAxis(){
		double	R2D = 57.29577951308232090712;
		double distance = Math.sqrt( coord[ 0 ] * coord[ 0 ] +
			coord[ 1 ] * coord[ 1 ] + coord[ 2 ] * coord[ 2 ] );
		double xComponent = coord[ 0 ] / distance;
		double zComponent = coord[ 2 ] / distance;
		double theta = 0.0;
		if ( ( zComponent == 0 ) && ( xComponent > 0 ) )
			theta = 180.0;
		else if ( ( zComponent == 0 ) && ( xComponent <= 0 ) )
			theta = 0.0;
		else if ( ( xComponent == 0 ) && ( zComponent > 0 ) )
			theta = 90.0;
		else if ( ( xComponent == 0 ) && ( zComponent < 0 ) )
			theta = 270.0;
		else{ 
			theta = Math.atan( Math.abs( zComponent / xComponent ) ) * R2D;
			if ( ( zComponent > 0 ) && ( xComponent > 0 ) )
				theta = 90.0 + Math.atan( Math.abs( xComponent / 
					zComponent ) ) * R2D;
			else if ( ( zComponent < 0 ) && ( xComponent > 0 ) )
				theta += 180.0;
			else if ( ( zComponent < 0 ) && ( xComponent < 0 ) )
				theta = 270.0 + Math.atan( Math.abs( xComponent / 
					zComponent ) ) * R2D;
		}
		return( theta );
	}

	// This returns the right handed rotation angle around the +z
	//  axis where the +x axis is the 0
	// This was not written by me, but I have checked it
	public double angleAboutZAxis(){
		double	R2D = 57.29577951308232090712;
		double distance = Math.sqrt( coord[ 0 ] * coord[ 0 ] +
			coord[ 1 ] * coord[ 1 ] + coord[ 2 ] * coord[ 2 ] );
		double xComponent = coord[ 0 ] / distance;
		double yComponent = coord[ 1 ] / distance;
		double theta = 0.0;
		if ( ( xComponent == 0 ) && ( yComponent > 0 ) )
			theta = 90.0;
		else if ( ( xComponent == 0 ) && ( yComponent < 0 ) )
			theta = 270.0;
		else if ( ( yComponent == 0 ) && ( xComponent >= 0 ) )
			theta = 0.0;
		else if ( ( yComponent == 0 ) && ( xComponent < 0 ) )
			theta = 180.0;
		else{ 
			theta = Math.atan( Math.abs( yComponent / xComponent ) ) * R2D;
			if ( ( yComponent > 0 ) && ( xComponent < 0 ) )
				theta = 90.0 + Math.atan( Math.abs( xComponent / 
					yComponent ) ) * R2D;
			else if ( ( yComponent < 0 ) && ( xComponent < 0 ) )
				theta += 180.0;
			else if ( ( yComponent < 0 ) && ( xComponent > 0 ) )
				theta = 270.0 + Math.atan( Math.abs( xComponent / 
					yComponent ) ) * R2D;
		}
		return( theta );
	}

	// Returns the angle (in degrees) made between atom1-atom2-thisatom
	// This was not written by me, but I have checked it
	public double angle(Atom atom1, Atom atom2){	
		double R2D = 57.29577951308232090712;
		return(R2D * angleInRadians(atom1,atom2));
	}

	// Returns the angle (in radians) made between atom1-atom2-thisatom
	// This was not written by me, but I have checked it
	public double angleInRadians(Atom atom1, Atom atom2){	
		double x12, x32, y12, y32, z12, z32, l12, l32, dp;
		double	R2D = 57.29577951308232090712;
	
    x12 = atom1.coord[0] - atom2.coord[0];
    y12 = atom1.coord[1] - atom2.coord[1];
    z12 = atom1.coord[2] - atom2.coord[2];
    x32 = coord[0] - atom2.coord[0];
    y32 = coord[1] - atom2.coord[1];
    z32 = coord[2] - atom2.coord[2];
    l12 = Math.sqrt(x12 * x12 + y12 * y12 + z12 * z12);
    l32 = Math.sqrt(x32 * x32 + y32 * y32 + z32 * z32);
    if(l12 == 0.0){
			return(0.0);
    }
    if(l32 == 0.0){
			return(0.0);
    }
    dp = (x12 * x32 + y12 * y32 + z12 * z32) / (l12 * l32);
		if (dp < -1.0)
			dp = -1.0;
		else if (dp > 1.0)
			dp = 1.0;
   	return(Math.acos(dp));
  }

	// Returns the torsion angle (in degrees) made between
	//  atom1-atom2-atom3-thisatom
	// This was not written by me, but I have checked it
	// If all 4 atoms lie in a plane and the first and fourth
	//  atoms are trans then 180 is returned, if they are cis
	//  then 0 is returned.
	// Between these two extremes, the angle of the right
	//  handed rotation where the axis is the vector from
	//  atom2 to atom3 (ie thumb points to atom3) is returned.
	// The returned angle is between -180-epsilon .. +180
	public double torsion(Atom atom1, Atom atom2, Atom atom3){
		double xij, yij, zij;
		double xkj, ykj, zkj;
		double xkl, ykl, zkl;
		double dx, dy, dz;
		double gx, gy, gz;
		double bi, bk;
		double ct, d, ap, app, bibk;

    xij = atom1.coord[0] - atom2.coord[0];
    yij = atom1.coord[1] - atom2.coord[1];
    zij = atom1.coord[2] - atom2.coord[2];
    xkj = atom3.coord[0] - atom2.coord[0];
    ykj = atom3.coord[1] - atom2.coord[1];
    zkj = atom3.coord[2] - atom2.coord[2];
    xkl = atom3.coord[0] - coord[0];
    ykl = atom3.coord[1] - coord[1];
    zkl = atom3.coord[2] - coord[2];

		// d = ij cross kj
		// g = kl cross kj
    dx = yij * zkj - zij * ykj;
    dy = zij * xkj - xij * zkj;
    dz = xij * ykj - yij * xkj;
    gx = zkj * ykl - ykj * zkl;
    gy = xkj * zkl - zkj * xkl;
    gz = ykj * xkl - xkj * ykl;

    bi = dx * dx + dy * dy + dz * dz;  // magnitude of d
    bk = gx * gx + gy * gy + gz * gz;  // magnitude of g
    ct = dx * gx + dy * gy + dz * gz;  // d dot g
		bibk = bi * bk;
		if (bibk < 1.0e-6)	
			return 0;
    ct = ct / Math.sqrt(bibk);
    if(ct < -1.0)
      ct = -1.0;
    else if(ct > 1.0)
      ct = 1.0;

    ap = Math.acos(ct);
    d  = xkj*(dz*gy-dy*gz) + ykj*(dx*gz-dz*gx) + zkj*(dy*gx-dx*gy);
    if(d < 0.0)
      ap = -ap;
    ap = Math.PI - ap;
    app = 180.0 * ap / Math.PI;
    if(app > 180.0)
      app = app - 360.0;
    return(app);
	}

	// Returns the torsion angle (in radians) made between
	//  atom1-atom2-atom3-thisatom
	// This was not written by me, but I have checked it
	public double torsionInRadians(Atom atom1, Atom atom2, Atom atom3){
		double D2R = 0.01745329251994329576;
    return(D2R*torsion(atom1,atom2,atom3));
	}

	//Checks if this is a backbone atom (the atom name is one of N, CA, C, O, OXT, H, H1, H2, H3)
	public boolean setIsBBatom(){
		if ( name.equalsIgnoreCase("N") || name.equalsIgnoreCase("CA") || name.equalsIgnoreCase("C") || name.equalsIgnoreCase("O") 
				|| name.equalsIgnoreCase("OXT") || name.equalsIgnoreCase("H")  
				 || name.equalsIgnoreCase("H1") || name.equalsIgnoreCase("H2") || name.equalsIgnoreCase("H3") )
			return true;
		else
			return false;
	}
	
	/*//Checks if this is a backbone atom (the atom name is one of N, CA, C, O, CB, OXT, H, H1, H2, H3)
	public boolean setIsBBatom(){
		if ( name.equalsIgnoreCase("N") || name.equalsIgnoreCase("CA") || name.equalsIgnoreCase("C") || name.equalsIgnoreCase("O")
				|| name.equalsIgnoreCase("CB") || name.equalsIgnoreCase("OXT") || name.equalsIgnoreCase("H")  
				 || name.equalsIgnoreCase("H1") || name.equalsIgnoreCase("H2") || name.equalsIgnoreCase("H3") )
			return true;
		else
			return false;
	}*/
	
	//Checks if this is a backbone atom
	public boolean getIsBBatom(){
		return isBBatom;
	}

	//Hopefully these two functions will provide a deep 
	//copy of the current atom instead of just pointing 
	//to it
	public Atom copy(){
		return new Atom(this);
	}
	
	protected Atom (Atom a){
		moleculeAtomNumber = a.moleculeAtomNumber;
		residueAtomNumber = a.residueAtomNumber;
		modelAtomNumber = a.modelAtomNumber;
		moleculeResidueNumber = a.moleculeResidueNumber;
		strandResidueNumber = a.strandResidueNumber;
		strandNumber = a.strandNumber;
		elementNumber = a.elementNumber;
		numberOfBonds = a.numberOfBonds;
		elementType = a.elementType;
		forceFieldType = a.forceFieldType;
		type = a.type;
		selected = a.selected;
		name = a.name;
		segID = a.segID;
		charge = a.charge;
		radius = a.radius;
		mass = a.mass;
		isBBatom = a.isBBatom;
		if(a.bond != null){
			bond = new int[a.bond.length];
			for(int i=0; i<bond.length;i++)
				bond[i] = a.bond[i];
		}
		
		coord = new float[3];
		for(int i=0; i<3;i++)
			coord[i] = a.coord[i];
		
	}

	public void setCoords(float xpos, float ypos, float zpos){
		coord[0] = xpos;
		coord[1] = ypos;
		coord[2] = zpos;
	}
}
