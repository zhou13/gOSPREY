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
//	ProbeStericCheck.java
//
//	Version:           2.1 beta
//
//
//	  authors:
// 	  initials    name                 organization                email
//	 ---------   -----------------    ------------------------    ----------------------------
//	  ISG		 Ivelin Georgiev	  Duke University			  ivelin.georgiev@duke.edu
//
///////////////////////////////////////////////////////////////////////////////////////////////

/**
* Written by Ivelin Georgiev (2004-2009)
* 
*/

/**
 * Implements a steric check for atom pairs based on the Richardsons' Probe approach.
 */
public class ProbeStericCheck {
	
	final float hBondCutoff = 0.8f; //the allowed vdW overlap for H-bonding atoms (used for steric checks)

	//constructor
	ProbeStericCheck (){		
	}
	
	//Checks whether there is steric overlap (above the overlapThresh cutoff) between atoms a1 and a2 in molecule m;
	//A special cutoff (hBondCutoff) is used if the two atoms could be forming a hydrogen bond;
	//NOTE: the steric checks is performed using the actualCoordinates[] in the molecule, and not the atoms coord[]
	public boolean isAllowedSteric(Molecule m, Atom a1, Atom a2, float overlapThresh){
		
		float a1coord[] = m.getActualCoord(a1.moleculeAtomNumber);
		float a2Coord[] = m.getActualCoord(a2.moleculeAtomNumber);
		float dist = getDist(a1coord,a2Coord);
		if ((dist < ((a2.radius + a1.radius)/100.0) - overlapThresh)){ //the two atoms overlap more than allowed
			if (areNonBondedSteric(m,a1,a2)) { //check for overlap betwen nonbonded pairs only			
				if (checkHbondTypes(a1,a2)){ //possible H bonding atoms, so allow bigger overlap
					if (dist < ((a2.radius + a1.radius)/100.0) - hBondCutoff)
						return false;
					else
						return true;
				}
				else
					return false;
			}
			else 
				return true;
		}
		else
			return true;
	}
	
	//Checks if atoms a1 and a2 are non-bonded (used in steric checks only; should not be used for the Amber energy computation);
	//In the Richardsons' Probe, for steric check purposes, two atoms are considered non-bonded if:
	//		1) They are more than 3 bonds apart if neither of the atoms is a hydrogen, or
	//		2) They are more than 4 bonds apart if at least one of the atoms is a hydrogen
	private boolean areNonBondedSteric(Molecule m, Atom a1, Atom a2){
		
		if (!m.connectivityValid) //valid molecule connectivity is required here
			m.establishConnectivity(false);
		
		int a1num = a1.moleculeAtomNumber;
		int a2num = a2.moleculeAtomNumber;
		
		if ( (!a1.elementType.equalsIgnoreCase("H")) && (!a2.elementType.equalsIgnoreCase("H")) ) //not hydrogens
			return (m.areNonBonded(a1num, a2num));
		
		else { //one of the atoms is a hydrogen
			
			if (!m.areNonBonded(a1num, a2num)) //the two atoms are at most 3 bonds apart
				return false;
			
			else { //the two atoms are either exactly 4 bonds apart or more than 4 bonds apart
				
				int i = 0;
				while ( i < m.connected14.length ){
					
					if ( (m.connected14[i]==a1num) && (m.atom[m.connected14[i+3]].bondedTo(a2num)) ) //exactly 4 bonds apart
						return false;
					
					else if ( (m.connected14[i+3]==a1num) && (m.atom[m.connected14[i]].bondedTo(a2num)) ) //exactly 4 bonds apart
						return false;
					
					else if ( (m.connected14[i]==a2num) && (m.atom[m.connected14[i+3]].bondedTo(a1num)) ) //exactly 4 bonds apart
						return false;
					
					else if ( (m.connected14[i+3]==a2num) && (m.atom[m.connected14[i]].bondedTo(a1num)) ) //exactly 4 bonds apart
						return false;
					
					i = i+4;
				}
				
				//if reached here, then the two atoms are not exactly 4 bonds apart, so they must be more than 4 bonds apart
				return true;
			}
		}
	}
	
	//Checks if atoms a1 and a2 are a possible H bonding pair (based only on the atom types)
	private boolean checkHbondTypes(Atom a1, Atom a2){
		
		String a1type = a1.forceFieldType;
		String a2type = a2.forceFieldType;
		
		if ( a1type.equalsIgnoreCase("H") || a1type.equalsIgnoreCase("HO") || a1type.equalsIgnoreCase("HS") || a1type.equalsIgnoreCase("HW")) {
			if ( a2type.equalsIgnoreCase("N") || a2type.equalsIgnoreCase("O") || a2type.equalsIgnoreCase("S") )
				return true;
		}
		else if ( a2type.equalsIgnoreCase("H") || a2type.equalsIgnoreCase("HO") || a2type.equalsIgnoreCase("HS") || a2type.equalsIgnoreCase("HW")) {
			if ( a1type.equalsIgnoreCase("N") || a1type.equalsIgnoreCase("O") || a1type.equalsIgnoreCase("S") )
				return true;
		}
		
		return false;
	}
	
	//Returns the distance between the two points given by coordinates coord1[] and coord2[] (float version)
	private float getDist(float coord1[], float coord2[]){
		float rijx, rijy, rijz, rij, rij2;
		
		rijx = coord1[0] - coord2[0];
		rijy = coord1[1] - coord2[1];
		rijz = coord1[2] - coord2[2];
		rij2 = rijx * rijx + rijy * rijy + rijz * rijz;
		rij = (float)Math.sqrt(rij2);
		
		return rij;
	}
}
