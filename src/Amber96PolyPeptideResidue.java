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
// Amber96PolyPeptideResidue.java
//
//  Version:           2.1 beta
//
//
// authors:
//    initials    name            organization                email
//   ---------   --------------  ------------------------    ------------------------------
//     RHL        Ryan Lilien     Dartmouth College           ryan.lilien@dartmouth.edu
//     KER        Kyle E. Roberts       Duke University         ker17@duke.edu
//     PGC        Pablo Gainza C.       Duke University         pablo.gainza@duke.edu
////////////////////////////////////////////////////////////////////////////////////////////

/*
 * Major changes were made by Ryan Lilien (2001-2004)
 * This is an extension to the polyPeptideResidue class
 *  It is hardcoded to the Amber96 atom types for each residue as specified
 *   in all_amino94.in **CURRENTLY FOR THE L FORMS ONLY**
 *  This is helpful in using Rotamers where you often switch residue type
 *   and don't want to have to recompute atom types
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

import java.awt.*;

/**
 * This class contains hard-coded templates for the different amino acid types. These templates are
 * used when performing residue mutations.
 */
public class Amber96PolyPeptideResidue {

    AminoAcidTemplates aat = null;
    Residue residue;

    public Amber96PolyPeptideResidue() {

        loadAminoAcidTemplates();


    }

    private void loadAminoAcidTemplates() {
        if(aat == null) {
            try {
                aat = new AminoAcidTemplates();
            } catch (Exception E) {
                System.out.println("Problem Loading Amino Acid Templates");
            }
        }

    }

    Residue getResidue( String residueName ) {

        int numAtoms = 0;
        int aaResIndex = -1;
        for(int i=0; i<aat.numAAs; i++) {
            if(aat.aaResidues[i].name.equalsIgnoreCase(residueName)) {
                numAtoms = aat.aaResidues[i].numberOfAtoms;
                aaResIndex = i;
                break;
            }
        }

        createResidue(residueName, numAtoms);
        //KER: Need to make sure that when we copy the residue we are
        //doing a deep copy and not a shallow copy
        for(int i=0; i<numAtoms; i++) {
            residue.atom[i] = aat.aaResidues[aaResIndex].atom[i].copy();
        }


        // RHL: If only a three character residue name is given, assume
        //  we want the L form
        /*if (residueName.length() == 3) {
        	String tmp = new String("L");
        	residueName = tmp.concat(residueName);
        }

        if ( residueName.equalsIgnoreCase( "Lace" ) ){
        	createResidue( "ACE", 6 );
        residue.atom[ 0 ] = new Atom( "H1", 0.000f, 0.000f, 0.000f );
        residue.atom[ 1 ] = new Atom( "CH3", 0.000f, 1.090f, 0.000f );
        residue.atom[ 2 ] = new Atom( "H2", 1.028f, 1.453f, 0.000f );
        residue.atom[ 3 ] = new Atom( "H3", -0.514f, 1.453f, -0.890f );
        residue.atom[ 4 ] = new Atom( "C", -0.721f, 1.600f, 1.249f, 0.616f );
        residue.atom[ 5 ] = new Atom( "O", -0.839f, 2.806f, 1.453f, -0.504f );
        	residue.atom[ 0 ].addBond( 1 );
        	residue.atom[ 1 ].addBond( 2, 3 );
        	residue.atom[ 1 ].addBond( 4 );
        	residue.atom[ 4 ].addBond( 5 );
        }
        else if ( residueName.equalsIgnoreCase( "Lala" ) ){
        	createResidue( "ALA", 10 );
        residue.atom[ 0 ] = new Atom( "N", 2.044f, 1.352f, 0.000f, -0.4157f, "N" );
        residue.atom[ 1 ] = new Atom( "H", 1.721f, 1.837f, 0.824f, 0.2719f, "H" );
        residue.atom[ 2 ] = new Atom( "CA", 1.522f, 0.000f, 0.000f, 0.03370f, "CT" );
        residue.atom[ 3 ] = new Atom( "HA", 1.896f, -0.481f, -0.904f, 0.0823f, "H1" );
        residue.atom[ 4 ] = new Atom( "C", 0.000f, 0.000f, 0.000f, 0.5973f, "C" );
        residue.atom[ 5 ] = new Atom( "O", -0.624f, 1.060f, 0.000f, -0.5679f, "O" );
        residue.atom[ 6 ] = new Atom( "CB", 1.988f, -0.769f, 1.233f, -0.1825f, "CT");
        residue.atom[ 7 ] = new Atom( "HB1", 1.580f, -1.780f, 1.205f, 0.0603f, "HC" );
        residue.atom[ 8 ] = new Atom( "HB2", 3.077f, -0.816f, 1.241f, 0.0603f, "HC" );
        	residue.atom[ 9 ] = new Atom( "HB3", 1.638f, -0.260f, 2.131f, 0.0603f, "HC" );
        	addBackboneBonds();
        	residue.atom[ 4 ].addBond( 5 );
        	residue.atom[ 2 ].addBond( 6 );
        	residue.atom[ 6 ].addBond( 7, 8, 9 );
        }
        else if ( residueName.equalsIgnoreCase( "Larg" ) ){
        	createResidue( "ARG", 24 );
        residue.atom[ 0 ] = new Atom( "N", 2.044f, 1.352f, 0.000f, -0.3479f, "N" );
        residue.atom[ 1 ] = new Atom( "H", 1.721f, 1.837f, 0.824f, 0.2747f, "H" );
        residue.atom[ 2 ] = new Atom( "CA", 1.522f, 0.000f, 0.000f, -0.2637f, "CT" );
        residue.atom[ 3 ] = new Atom( "HA", 1.896f, -0.481f, -0.904f, 0.1560f, "H1" );
        residue.atom[ 4 ] = new Atom( "C", 0.000f, 0.000f, 0.000f, 0.7341f, "C" );
        	residue.atom[ 5 ] = new Atom( "O", -0.624f, 1.060f, 0.000f, -0.5894f, "O" );
        	residue.atom[ 6 ] = new Atom( "CB", 1.988f, -0.769f, 1.233f, -0.0007f, "CT" );
        	residue.atom[ 7 ] = new Atom( "HB1", 3.077f, -0.816f, 1.241f, 0.0327f, "HC" );
        	residue.atom[ 8 ] = new Atom( "HB2", 1.638f, -0.260f, 2.131f, 0.0327f, "HC" );
        	residue.atom[ 9 ] = new Atom( "CG", 1.418f, -2.182f, 1.194f, 0.0390f, "CT" );
        	residue.atom[ 10 ] = new Atom( "HG1", 0.329f, -2.136f, 1.186f, 0.0285f, "HC" );
        	residue.atom[ 11 ] = new Atom( "HG2", 1.767f, -2.692f, 0.296f, 0.0285f, "HC" );
        	residue.atom[ 12 ] = new Atom( "CD", 1.884f, -2.951f, 2.427f, 0.0486f, "CT" );
        	residue.atom[ 13 ] = new Atom( "HD1", 2.973f, -2.984f, 2.445f, 0.0687f, "H1" );
        	residue.atom[ 14 ] = new Atom( "HD2", 1.511f, -2.438f, 3.314f, 0.0687f, "H1" );
        	residue.atom[ 15 ] = new Atom( "NE", 1.349f, -4.332f, 2.424f, -0.5295f, "N2" );
        	residue.atom[ 16 ] = new Atom( "HE", 0.761f, -4.619f, 1.655f, 0.3456f, "H" );
        	residue.atom[ 17 ] = new Atom( "CZ", 1.606f, -5.209f, 3.389f, 0.8076f, "CA" );
        	residue.atom[ 18 ] = new Atom( "NH1", 2.372f, -4.905f, 4.434f, -0.8627f, "N2" );
        	residue.atom[ 19 ] = new Atom( "HH11", 2.774f, -3.982f, 4.509f, 0.4478f, "H" );
        	residue.atom[ 20 ] = new Atom( "HH12", 2.546f, -5.597f, 5.148f, 0.4478f, "H" );
        	residue.atom[ 21 ] = new Atom( "NH2", 1.074f, -6.424f, 3.287f, -0.8627f, "N2" );
        	residue.atom[ 22 ] = new Atom( "HH21", 0.494f, -6.653f, 2.493f, 0.4478f, "H" );
        	residue.atom[ 23 ] = new Atom( "HH22", 1.252f, -7.113f, 4.004f, 0.4478f, "H" );
        	addBackboneBonds();
        	residue.atom[ 4 ].addBond( 5 );
        	residue.atom[ 2 ].addBond( 6 );
        	residue.atom[ 6 ].addBond( 7, 8, 9 );
        	residue.atom[ 9 ].addBond( 10, 11, 12 );
        	residue.atom[ 12 ].addBond( 13, 14, 15 );
        	residue.atom[ 15 ].addBond( 16, 17 );
        	residue.atom[ 17 ].addBond( 18, 21 );
        	residue.atom[ 18 ].addBond( 19, 20 );
        	residue.atom[ 21 ].addBond( 22, 23 );
        }
        else if ( residueName.equalsIgnoreCase( "Lasn" ) ){
        	createResidue( "ASN", 14 );
        residue.atom[ 0 ] = new Atom( "N", 2.044f, 1.352f, 0.000f, -0.4157f, "N" );
        residue.atom[ 1 ] = new Atom( "H", 1.721f, 1.837f, 0.824f, 0.2719f, "H" );
        residue.atom[ 2 ] = new Atom( "CA", 1.522f, 0.000f, 0.000f, 0.0143f, "CA" );
        residue.atom[ 3 ] = new Atom( "HA", 1.896f, -0.481f, -0.904f, 0.1048f, "H1" );
        residue.atom[ 4 ] = new Atom( "C", 0.000f, 0.000f, 0.000f, 0.5973f, "C" );
        	residue.atom[ 5 ] = new Atom( "O", -0.624f, 1.060f, 0.000f, -0.5679f, "O" );
        residue.atom[ 6 ] = new Atom( "CB", 1.988f, -0.769f, 1.233f, -0.2041f, "CT" );
        residue.atom[ 7 ] = new Atom( "HB1", 3.077f, -0.800f, 1.251f, 0.0797f, "HC" );
        residue.atom[ 8 ] = new Atom( "HB2", 1.613f, -0.256f, 2.118f, 0.0797f, "HC" );
        residue.atom[ 9 ] = new Atom( "CG", 1.439f, -2.188f, 1.232f, 0.713f, "C" );
        residue.atom[ 10 ] = new Atom( "OD1", 0.720f, -2.579f, 0.315f, -0.5931f, "O" );
        residue.atom[ 11 ] = new Atom( "ND2", 1.780f, -2.961f, 2.266f, -0.9191f, "N" );
        residue.atom[ 12 ] = new Atom( "HD21", 2.374f, -2.591f, 2.995f, 0.4196f, "H" );
        residue.atom[ 13 ] = new Atom( "HD22", 1.443f, -3.912f, 2.315f, 0.4196f, "H" );
        	addBackboneBonds();
        	residue.atom[ 4 ].addBond( 5 );
        	residue.atom[ 2 ].addBond( 6 );
        	residue.atom[ 6 ].addBond( 7, 8, 9 );
        	residue.atom[ 9 ].addBond( 10, 11 );
        	residue.atom[ 11 ].addBond( 12, 13 );
        }
        else if ( residueName.equalsIgnoreCase( "Lasp" ) ){
        	createResidue( "ASP", 12 );
        residue.atom[ 0 ] = new Atom( "N", 2.044f, 1.352f, 0.000f, -0.5163f, "N" );
        residue.atom[ 1 ] = new Atom( "H", 1.721f, 1.837f, 0.824f, 0.2936f, "H" );
        residue.atom[ 2 ] = new Atom( "CA", 1.522f, 0.000f, 0.000f, 0.0381f, "CT" );
        residue.atom[ 3 ] = new Atom( "HA", 1.896f, -0.481f, -0.904f, 0.0880f, "H1" );
        residue.atom[ 4 ] = new Atom( "C", 0.000f, 0.000f, 0.000f, 0.5366f, "C" );
        	residue.atom[ 5 ] = new Atom( "O", -0.624f, 1.060f, 0.000f, -0.5819f, "O" );
        residue.atom[ 6 ] = new Atom( "CB", 1.988f, -0.769f, 1.233f, -0.0303f, "CT" );
        residue.atom[ 7 ] = new Atom( "HB1", 3.077f, -0.816f, 1.241f, -0.0122f, "HC" );
        residue.atom[ 8 ] = new Atom( "HB2", 1.638f, -0.260f, 2.131f, -0.0122f, "HC" );
        residue.atom[ 9 ] = new Atom( "CG", 1.417f, -2.185f, 1.194f, 0.7994f, "C" );
        residue.atom[ 10 ] = new Atom( "OD1", 2.107f, -3.069f, 0.620f, -0.8014f, "O2" );
        residue.atom[ 11 ] = new Atom( "OD2", 0.297f, -2.369f, 1.741f, -0.8014f, "O2" );
        	addBackboneBonds();
        	residue.atom[ 4 ].addBond( 5 );
        	residue.atom[ 2 ].addBond( 6 );
        	residue.atom[ 6 ].addBond( 7, 8, 9 );
        	residue.atom[ 9 ].addBond( 10, 11 );
        }
        else if ( residueName.equalsIgnoreCase( "Lcys" ) ){
        	createResidue( "CYS", 11 );
        residue.atom[ 0 ] = new Atom( "N", 2.044f, 1.352f, 0.000f, -0.4157f, "N" );
        residue.atom[ 1 ] = new Atom( "H", 1.721f, 1.837f, 0.824f, 0.2719f, "H" );
        residue.atom[ 2 ] = new Atom( "CA", 1.522f, 0.000f, 0.000f, 0.0213f, "CT" );
        residue.atom[ 3 ] = new Atom( "HA", 1.896f, -0.481f, -0.904f, 0.1124f, "H1" );
        residue.atom[ 4 ] = new Atom( "C", 0.000f, 0.000f, 0.000f, 0.5973f, "C" );
        	residue.atom[ 5 ] = new Atom( "O", -0.624f, 1.060f, 0.000f, -0.5679f, "O" );
        residue.atom[ 6 ] = new Atom( "CB", 1.988f, -0.769f, 1.233f, -0.1231f, "CT" );
        residue.atom[ 7 ] = new Atom( "HB1", 3.076f, -0.748f, 1.283f, 0.1112f, "H1" );
        residue.atom[ 8 ] = new Atom( "HB2", 1.528f, -0.252f, 2.075f, 0.1112f, "H1" );
        residue.atom[ 9 ] = new Atom( "SG", 1.409f, -2.479f, 1.367f, -0.3119f, "SH" );
        residue.atom[ 10 ] = new Atom( "HG", 1.890f, -3.023f, 2.481f, 0.1933f, "HS" );
        	addBackboneBonds();
        	residue.atom[ 4 ].addBond( 5 );
        	residue.atom[ 2 ].addBond( 6 );
        	residue.atom[ 6 ].addBond( 7, 8, 9 );
        	residue.atom[ 9 ].addBond( 10 );
        }
        else if ( residueName.equalsIgnoreCase( "Lcyx" ) ){
        	createResidue( "CYX", 10 );
        residue.atom[ 0 ] = new Atom( "N", 2.044f, 1.352f, 0.000f, -0.4157f, "N" );
        residue.atom[ 1 ] = new Atom( "H", 1.721f, 1.837f, 0.824f, 0.2719f, "H" );
        residue.atom[ 2 ] = new Atom( "CA", 1.522f, 0.000f, 0.000f, 0.0429f, "CT" );
        residue.atom[ 3 ] = new Atom( "HA", 1.896f, -0.481f, -0.904f, 0.0766f, "H1" );
        residue.atom[ 4 ] = new Atom( "C", 0.000f, 0.000f, 0.000f, 0.5973f, "C" );
        	residue.atom[ 5 ] = new Atom( "O", -0.624f, 1.060f, 0.000f, -0.5679f, "O" );
        residue.atom[ 6 ] = new Atom( "CB", 1.988f, -0.769f, 1.233f, -0.0790f, "CT" );
        residue.atom[ 7 ] = new Atom( "HB1", 3.076f, -0.748f, 1.283f, 0.0910f, "H1" );
        residue.atom[ 8 ] = new Atom( "HB2", 1.528f, -0.252f, 2.075f, 0.0910f, "H1" );
        residue.atom[ 9 ] = new Atom( "SG", 1.409f, -2.479f, 1.367f, -0.1081f, "S" );
        	addBackboneBonds();
        	residue.atom[ 4 ].addBond( 5 );
        	residue.atom[ 2 ].addBond( 6 );
        	residue.atom[ 6 ].addBond( 7, 8, 9 );
        }
        else if ( residueName.equalsIgnoreCase( "Lgln" ) ){
        	createResidue( "GLN", 17 );
        	residue.atom[ 0 ] = new Atom( "N", 2.044f, 1.352f, 0.000f, -0.4157f, "N" );
        	residue.atom[ 1 ] = new Atom( "H", 1.721f, 1.837f, 0.824f, 0.2719f, "H" );
        	residue.atom[ 2 ] = new Atom( "CA", 1.522f, 0.000f, 0.000f, -0.0031f, "CT" );
        	residue.atom[ 3 ] = new Atom( "HA", 1.896f, -0.481f, -0.904f, 0.0850f, "H1" );
        	residue.atom[ 4 ] = new Atom( "C", 0.000f, 0.000f, 0.000f, 0.5973f, "C" );
        	residue.atom[ 5 ] = new Atom( "O", -0.624f, 1.060f, 0.000f, -0.5679f, "O" );
        residue.atom[ 6 ] = new Atom( "CB", 1.988f, -0.769f, 1.233f, -0.0036f, "CT" );
        residue.atom[ 7 ] = new Atom( "HB1", 3.077f, -0.816f, 1.241f, 0.0171f, "HC" );
        residue.atom[ 8 ] = new Atom( "HB2", 1.638f, -0.260f, 2.131f, 0.0171f, "HC" );
        residue.atom[ 9 ] = new Atom( "CG", 1.418f, -2.182f, 1.194f, -0.0645f, "CT" );
        residue.atom[ 10 ] = new Atom( "HG1", 0.330f, -2.135f, 1.168f, 0.0352f, "HC" );
        residue.atom[ 11 ] = new Atom( "HG2", 1.792f, -2.681f, 0.300f, 0.0352f, "HC" );
        residue.atom[ 12 ] = new Atom( "CD", 1.861f, -2.984f, 2.410f, 0.6951f, "C" );
        residue.atom[ 13 ] = new Atom( "OE1", 2.585f, -2.476f, 3.263f, -0.6086f, "O" );
        residue.atom[ 14 ] = new Atom( "NE2", 1.422f, -4.243f, 2.489f, -0.9407f, "N" );
        residue.atom[ 15 ] = new Atom( "HE21", 0.828f, -4.614f, 1.761f, 0.4251f, "H" );
        residue.atom[ 16 ] = new Atom( "HE22", 1.687f, -4.819f, 3.275f, 0.4251f, "H" );
        	addBackboneBonds();
        	residue.atom[ 4 ].addBond( 5 );
        	residue.atom[ 2 ].addBond( 6 );
        	residue.atom[ 6 ].addBond( 7, 8, 9 );
        	residue.atom[ 9 ].addBond( 10, 11, 12 );
        	residue.atom[ 12 ].addBond( 13, 14 );
        	residue.atom[ 14 ].addBond( 15, 16 );
        }
        else if ( residueName.equalsIgnoreCase( "Lglu" ) ){
        	createResidue( "GLU", 15 );
        residue.atom[ 0 ] = new Atom( "N", 2.044f, 1.352f, 0.000f, -0.5163f, "N" );
        residue.atom[ 1 ] = new Atom( "H", 1.721f, 1.837f, 0.824f, 0.2936f, "H" );
        residue.atom[ 2 ] = new Atom( "CA", 1.522f, 0.000f, 0.000f, 0.0397f, "CT" );
        residue.atom[ 3 ] = new Atom( "HA", 1.896f, -0.481f, -0.904f, 0.1105f, "H1" );
        residue.atom[ 4 ] = new Atom( "C", 0.000f, 0.000f, 0.000f, 0.5366f, "C" );
        	residue.atom[ 5 ] = new Atom( "O", -0.624f, 1.060f, 0.000f, -0.5819f, "O" );
        residue.atom[ 6 ] = new Atom( "CB", 1.988f, -0.769f, 1.233f, 0.0560f, "CT" );
        residue.atom[ 7 ] = new Atom( "HB1", 3.077f, -0.816f, 1.241f, -0.0173f, "HC" );
        residue.atom[ 8 ] = new Atom( "HB2", 1.638f, -0.260f, 2.131f, -0.0173f, "HC" );
        residue.atom[ 9 ] = new Atom( "CG", 1.423f, -2.168f, 1.195f, 0.0136f, "CT" );
        residue.atom[ 10 ] = new Atom( "HG1", 0.334f, -2.122f, 1.187f, -0.0425f, "HC" );
        residue.atom[ 11 ] = new Atom( "HG2", 1.772f, -2.678f, 0.296f, -0.0425f, "HC" );
        residue.atom[ 12 ] = new Atom( "CD", 1.890f, -2.939f, 2.429f, 0.8054f, "C" );
        residue.atom[ 13 ] = new Atom( "OE1", 1.161f, -2.878f, 3.455f, -0.8188f, "O2" );
        residue.atom[ 14 ] = new Atom( "OE2", 2.971f, -3.578f, 2.334f, -0.8188f, "O2" );
        	addBackboneBonds();
        	residue.atom[ 4 ].addBond( 5 );
        	residue.atom[ 2 ].addBond( 6 );
        	residue.atom[ 6 ].addBond( 7, 8, 9 );
        	residue.atom[ 9 ].addBond( 10, 11, 12 );
        	residue.atom[ 12 ].addBond( 13, 14 );
        }
        else if ( residueName.equalsIgnoreCase( "Lgly" ) ){
        	createResidue( "GLY", 7 );
        residue.atom[ 0 ] = new Atom( "N", 2.027f, 1.358f, 0.000f, -0.4157f, "N" );
        residue.atom[ 1 ] = new Atom( "H", 1.697f, 1.839f, 0.824f, 0.2719f, "H" );
        residue.atom[ 2 ] = new Atom( "CA", 1.522f, 0.000f, 0.000f, -0.0252f, "CT" );
        residue.atom[ 3 ] = new Atom( "HA2", 1.886f, -0.523f, -0.885f, 0.0698f, "H1" );
        residue.atom[ 4 ] = new Atom( "HA3", 1.874f, -0.506f, 0.899f, 0.0698f, "H1" );
        residue.atom[ 5 ] = new Atom( "C", 0.000f, 0.000f, 0.000f, 0.5973f, "C" );
        residue.atom[ 6 ] = new Atom( "O", -0.624f, 1.058f, 0.000f, -0.5679f, "O" );
        	addBackboneBonds();
        	residue.atom[ 2 ].addBond( 5 );
        	residue.atom[ 5 ].addBond( 6 );
        }
        else if ( residueName.equalsIgnoreCase( "Lhid" ) ){
        	createResidue( "HID", 17 );
        residue.atom[ 0 ] = new Atom( "N", 2.044f, 1.352f, 0.000f, -0.4157f, "N" );
        residue.atom[ 1 ] = new Atom( "H", 1.721f, 1.837f, 0.824f, 0.2719f, "H" );
        residue.atom[ 2 ] = new Atom( "CA", 1.522f, 0.000f, 0.000f, 0.0188f, "CT" );
        residue.atom[ 3 ] = new Atom( "HA", 1.896f, -0.481f, -0.904f, 0.0881f, "H1" );
        residue.atom[ 4 ] = new Atom( "C", 0.000f, 0.000f, 0.000f, 0.5973f, "C" );
        	residue.atom[ 5 ] = new Atom( "O", -0.624f, 1.060f, 0.000f, -0.5679f, "O" );
        	residue.atom[ 6 ] = new Atom( "CB", 1.988f, -0.769f, 1.233f, -0.0462f, "CT" );
        	residue.atom[ 7 ] = new Atom( "HB1", 3.076f, -0.759f, 1.276f, 0.0402f, "HC" );
        	residue.atom[ 8 ] = new Atom( "HB2", 1.547f, -0.251f, 2.084f, 0.0402f, "HC" );
        	residue.atom[ 9 ] = new Atom( "CG", 1.492f, -2.192f, 1.321f, -0.0266f, "CC" );
        	residue.atom[ 10 ] = new Atom( "ND1", 1.829f, -3.024f, 2.383f, -0.3811f, "NA" );
        	residue.atom[ 11 ] = new Atom( "HD1", 2.411f, -2.777f, 3.169f, 0.3649f, "H" );
        	residue.atom[ 12 ] = new Atom( "CD2", 0.698f, -2.921f, 0.493f, 0.1292f, "CV" );
        	residue.atom[ 13 ] = new Atom( "HD2", 0.306f, -2.473f, -0.421f, 0.1147f, "H4" );
        	residue.atom[ 14 ] = new Atom( "CE1", 1.252f, -4.194f, 2.183f, 0.2057f, "CR" );
        	residue.atom[ 15 ] = new Atom( "HE1", 1.387f, -5.000f, 2.905f, 0.1392f, "H5" );
        	residue.atom[ 16 ] = new Atom( "NE2", 0.576f, -4.150f, 1.061f, -0.5727f, "NB" );
        	addBackboneBonds();
        	residue.atom[ 4 ].addBond( 5 );
        	residue.atom[ 2 ].addBond( 6 );
        	residue.atom[ 6 ].addBond( 7, 8, 9 );
        	residue.atom[ 9 ].addBond( 10, 12 );
        	residue.atom[ 10 ].addBond( 11, 14 );
        	residue.atom[ 12 ].addBond( 13, 16 );
        	residue.atom[ 14 ].addBond( 15, 16 );
        }
        else if ( residueName.equalsIgnoreCase( "Lhie" ) ){
        	createResidue( "HIE", 17 );
        residue.atom[ 0 ] = new Atom( "N", 2.044f, 1.352f, 0.000f, -0.4157f, "N" );
        residue.atom[ 1 ] = new Atom( "H", 1.721f, 1.837f, 0.824f, 0.2719f, "H" );
        residue.atom[ 2 ] = new Atom( "CA", 1.522f, 0.000f, 0.000f, -0.0581f, "CT" );
        residue.atom[ 3 ] = new Atom( "HA", 1.896f, -0.481f, -0.904f, 0.1360f, "H1" );
        residue.atom[ 4 ] = new Atom( "C", 0.000f, 0.000f, 0.000f, 0.5973f, "C" );
        	residue.atom[ 5 ] = new Atom( "O", -0.624f, 1.060f, 0.000f, -0.5679f, "O" );
        residue.atom[ 6 ] = new Atom( "CB", 1.988f, -0.769f, 1.233f, -0.0074f, "CT" );
        residue.atom[ 7 ] = new Atom( "HB1", 3.076f, -0.759f, 1.276f, 0.0367f, "HC" );
        residue.atom[ 8 ] = new Atom( "HB2", 1.547f, -0.251f, 2.084f, 0.0367f, "HC" );
        residue.atom[ 9 ] = new Atom( "CG", 1.492f, -2.192f, 1.321f, 0.1868f, "CC" );
        residue.atom[ 10 ] = new Atom( "ND1", 1.829f, -3.024f, 2.383f, -0.5432f, "NB" );
        residue.atom[ 11 ] = new Atom( "CD2", 0.698f, -2.921f, 0.493f, -0.2207f, "CW" );
        residue.atom[ 12 ] = new Atom( "HD2", 0.306f, -2.473f, -0.421f, 0.1862f, "H4" );
        residue.atom[ 13 ] = new Atom( "CE1", 1.252f, -4.194f, 2.183f, 0.1635f, "CR" );
        residue.atom[ 14 ] = new Atom( "HE1", 1.387f, -5.000f, 2.905f, 0.1435f, "H5" );
        residue.atom[ 15 ] = new Atom( "NE2", 0.576f, -4.150f, 1.061f, -0.2795f, "NA" );
        residue.atom[ 16 ] = new Atom( "HE2", 0.041f, -4.916f, 0.677f, 0.3339f, "H" );
        	addBackboneBonds();
        	residue.atom[ 4 ].addBond( 5 );
        	residue.atom[ 2 ].addBond( 6 );
        	residue.atom[ 6 ].addBond( 7, 8, 9 );
        	residue.atom[ 9 ].addBond( 10, 11 );
        	residue.atom[ 10 ].addBond( 13 );
        	residue.atom[ 11 ].addBond( 12, 15 );
        	residue.atom[ 13 ].addBond( 14, 15 );
        	residue.atom[ 15 ].addBond( 16 );
        }
        else if ( residueName.equalsIgnoreCase( "Lhip" ) ){
        	createResidue( "HIP", 18 );
        residue.atom[ 0 ] = new Atom( "N", 2.044f, 1.352f, 0.000f, -0.3479f, "N" );
        residue.atom[ 1 ] = new Atom( "H", 1.721f, 1.837f, 0.824f, 0.2747f, "H" );
        residue.atom[ 2 ] = new Atom( "CA", 1.522f, 0.000f, 0.000f, -0.1354f, "CT" );
        residue.atom[ 3 ] = new Atom( "HA", 1.896f, -0.481f, -0.904f, 0.1212f, "H1" );
        residue.atom[ 4 ] = new Atom( "C", 0.000f, 0.000f, 0.000f, 0.7341f, "C" );
        	residue.atom[ 5 ] = new Atom( "O", -0.624f, 1.060f, 0.000f, -0.5894f, "O" );
        residue.atom[ 6 ] = new Atom( "CB", 1.988f, -0.769f, 1.233f, -0.0414f, "CT" );
        residue.atom[ 7 ] = new Atom( "HB1", 3.076f, -0.759f, 1.276f, 0.0810f, "HC" );
        residue.atom[ 8 ] = new Atom( "HB2", 1.547f, -0.251f, 2.084f, 0.0810f, "HC" );
        residue.atom[ 9 ] = new Atom( "CG", 1.492f, -2.192f, 1.321f, -0.0012f, "CC" );
        residue.atom[ 10 ] = new Atom( "ND1", 1.829f, -3.024f, 2.383f, -0.1513f, "NA" );
        residue.atom[ 11 ] = new Atom( "HD1", 2.411f, -2.777f, 3.169f, 0.3866f, "H" );
        residue.atom[ 12 ] = new Atom( "CD2", 0.698f, -2.921f, 0.493f, -0.1141f, "CW" );
        residue.atom[ 13 ] = new Atom( "HD2", 0.306f, -2.473f, -0.421f, 0.2317f, "H4" );
        residue.atom[ 14 ] = new Atom( "CE1", 1.252f, -4.194f, 2.183f, -0.0170f, "CR" );
        residue.atom[ 15 ] = new Atom( "HE1", 1.387f, -5.000f, 2.905f, 0.2681f, "H5" );
        residue.atom[ 16 ] = new Atom( "NE2", 0.576f, -4.150f, 1.061f, -0.1718f, "NA" );
        residue.atom[ 17 ] = new Atom( "HE2", 0.041f, -4.916f, 0.677f, 0.3911f, "H" );
        	addBackboneBonds();
        	residue.atom[ 4 ].addBond( 5 );
        	residue.atom[ 2 ].addBond( 6 );
        	residue.atom[ 6 ].addBond( 7, 8, 9 );
        	residue.atom[ 9 ].addBond( 10, 12 );
        	residue.atom[ 10 ].addBond( 11, 14 );
        	residue.atom[ 12 ].addBond( 13, 16 );
        	residue.atom[ 14 ].addBond( 15, 16 );
        	residue.atom[ 16 ].addBond( 17 );
        }
        else if ( residueName.equalsIgnoreCase( "Lhis" ) ){
        	createResidue( "HIS", 18 );
        	// RHL charge and type assigned as if it here HIP
        residue.atom[ 0 ] = new Atom( "N", 2.044f, 1.352f, 0.000f, -0.3479f, "N" );
        residue.atom[ 1 ] = new Atom( "H", 1.721f, 1.837f, 0.824f, 0.2747f, "H" );
        residue.atom[ 2 ] = new Atom( "CA", 1.522f, 0.000f, 0.000f, -0.1354f, "CT" );
        residue.atom[ 3 ] = new Atom( "HA", 1.896f, -0.481f, -0.904f, 0.1212f, "H1" );
        residue.atom[ 4 ] = new Atom( "C", 0.000f, 0.000f, 0.000f, 0.7341f, "C" );
        	residue.atom[ 5 ] = new Atom( "O", -0.624f, 1.060f, 0.000f, -0.5894f, "O" );
        residue.atom[ 6 ] = new Atom( "CB", 1.988f, -0.769f, 1.233f, -0.0414f, "CT" );
        residue.atom[ 7 ] = new Atom( "HB1", 3.076f, -0.759f, 1.276f, 0.0810f, "HC" );
        residue.atom[ 8 ] = new Atom( "HB2", 1.547f, -0.251f, 2.084f, 0.0810f, "HC" );
        residue.atom[ 9 ] = new Atom( "CG", 1.492f, -2.192f, 1.321f, -0.0012f, "CC" );
        residue.atom[ 10 ] = new Atom( "ND1", 1.829f, -3.024f, 2.383f, -0.1513f, "NA" );
        residue.atom[ 11 ] = new Atom( "HD1", 2.411f, -2.777f, 3.169f, 0.3866f, "H" );
        residue.atom[ 12 ] = new Atom( "CD2", 0.698f, -2.921f, 0.493f, -0.1141f, "CW" );
        residue.atom[ 13 ] = new Atom( "HD2", 0.248f, -2.587f, -0.442f, 0.2317f, "H4" );
        residue.atom[ 14 ] = new Atom( "CE1", 1.252f, -4.194f, 2.183f, -0.0170f, "CR" );
        residue.atom[ 15 ] = new Atom( "HE1", 1.327f, -5.058f, 2.843f, 0.2681f, "H5" );
        residue.atom[ 16 ] = new Atom( "NE2", 0.576f, -4.150f, 1.061f, -0.1718f, "NA" );
        residue.atom[ 17 ] = new Atom( "HE2", 0.041f, -4.916f, 0.677f, 0.3911f, "H" );
        	addBackboneBonds();
        	residue.atom[ 4 ].addBond( 5 );
        	residue.atom[ 2 ].addBond( 6 );
        	residue.atom[ 6 ].addBond( 7, 8, 9 );
        	residue.atom[ 9 ].addBond( 10, 12 );
        	residue.atom[ 10 ].addBond( 11, 14 );
        	residue.atom[ 12 ].addBond( 13, 16 );
        	residue.atom[ 14 ].addBond( 15, 16 );
        	residue.atom[ 16 ].addBond( 17 );
        }
        else if ( residueName.equalsIgnoreCase( "Lile" ) ){
        	createResidue( "ILE", 19 );
        residue.atom[ 0 ] = new Atom( "N", 2.044f, 1.352f, 0.000f, -0.4157f, "N" );
        residue.atom[ 1 ] = new Atom( "H", 1.721f, 1.837f, 0.824f, 0.2719f, "H" );
        residue.atom[ 2 ] = new Atom( "CA", 1.522f, 0.000f, 0.000f, -0.0597f, "CT" );
        residue.atom[ 3 ] = new Atom( "HA", 1.875f, -0.500f, -0.902f, 0.0869f, "H1" );
        residue.atom[ 4 ] = new Atom( "C", 0.000f, 0.000f, 0.000f, 0.5973f, "C" );
        	residue.atom[ 5 ] = new Atom( "O", -0.624f, 1.060f, 0.000f, -0.5679f, "O" );
        residue.atom[ 6 ] = new Atom( "CB", 2.009f, -0.733f, 1.245f, 0.1303f, "CT" );
        residue.atom[ 7 ] = new Atom( "HB", 3.098f, -0.765f, 1.245f, 0.0187f, "HC" );
        residue.atom[ 8 ] = new Atom( "CG1", 1.459f, -2.156f, 1.245f, -0.0430f, "CT" );
        residue.atom[ 9 ] = new Atom( "HG11", 0.370f, -2.124f, 1.245f, 0.0236f, "HC" );
        residue.atom[ 10 ] = new Atom( "HG12", 1.807f, -2.680f, 0.355f, 0.0236f, "HC" );
        residue.atom[ 11 ] = new Atom( "CG2", 1.522f, 0.000f, 2.491f, -0.3204f, "CT" );
        residue.atom[ 12 ] = new Atom( "HG21", 1.870f, -0.524f, 3.381f, 0.0882f, "HC" );
        residue.atom[ 13 ] = new Atom( "HG22", 1.914f, 1.017f, 2.490f, 0.0882f, "HC" );
        residue.atom[ 14 ] = new Atom( "HG23", 0.432f, 0.032f, 2.491f, 0.0882f, "HC" );
        residue.atom[ 15 ] = new Atom( "CD1", 1.947f, -2.889f, 2.491f, -0.0660f, "CT" );
        residue.atom[ 16 ] = new Atom( "HD11", 1.554f, -3.906f, 2.490f, 0.0186f, "HC" );
        residue.atom[ 17 ] = new Atom( "HD12", 3.036f, -2.921f, 2.491f, 0.0186f, "HC" );
        residue.atom[ 18 ] = new Atom( "HD13", 1.599f, -2.365f, 3.381f, 0.0186f, "HC" );
        	addBackboneBonds();
        	residue.atom[ 4 ].addBond( 5 );
        	residue.atom[ 2 ].addBond( 6 );
        	residue.atom[ 6 ].addBond( 7, 8, 11 );
        	residue.atom[ 8 ].addBond( 9, 10, 15 );
        	residue.atom[ 11 ].addBond( 12, 13, 14 );
        	residue.atom[ 15 ].addBond( 16, 17, 18 );
        }

        else if ( residueName.equalsIgnoreCase( "Lleu" ) ){
        	createResidue( "LEU", 19 );
        residue.atom[ 0 ] = new Atom( "N", 2.044f, 1.352f, 0.000f, -0.4157f, "N" );
        residue.atom[ 1 ] = new Atom( "H", 1.721f, 1.837f, 0.824f, 0.2719f, "H" );
        residue.atom[ 2 ] = new Atom( "CA", 1.522f, 0.000f, 0.000f, -0.0518f, "CT" );
        residue.atom[ 3 ] = new Atom( "HA", 1.896f, -0.481f, -0.904f, 0.0922f, "H1" );
        residue.atom[ 4 ] = new Atom( "C", 0.000f, 0.000f, 0.000f, 0.5973f, "C" );
        	residue.atom[ 5 ] = new Atom( "O", -0.624f, 1.060f, 0.000f, -0.5679f, "O" );
        	residue.atom[ 6 ] = new Atom( "CB", 1.988f, -0.769f, 1.233f, -0.1102f, "CT" );
        	residue.atom[ 7 ] = new Atom( "HB1", 3.077f, -0.816f, 1.241f, 0.0457f, "HC" );
        	residue.atom[ 8 ] = new Atom( "HB2", 1.638f, -0.260f, 2.131f, 0.0457f, "HC" );
        	residue.atom[ 9 ] = new Atom( "CG", 1.418f, -2.182f, 1.194f, 0.3531f, "CT" );
        	residue.atom[ 10 ] = new Atom( "HG", 0.329f, -2.136f, 1.186f, -0.0361f, "HC" );
        	residue.atom[ 11 ] = new Atom( "CD1", 1.906f, -2.894f, -0.063f, -0.4121f, "CT" );
        	residue.atom[ 12 ] = new Atom( "HD11", 1.499f, -3.905f, -0.090f, 0.1000f, "HC" );
        	residue.atom[ 13 ] = new Atom( "HD12", 1.573f, -2.345f, -0.943f, 0.1000f, "HC" );
        	residue.atom[ 14 ] = new Atom( "HD13", 2.995f, -2.941f, -0.055f, 0.1000f, "HC" );
        	residue.atom[ 15 ] = new Atom( "CD2", 1.884f, -2.951f, 2.427f, -0.4121f, "CT" );
        	residue.atom[ 16 ] = new Atom( "HD21", 1.476f, -3.962f, 2.400f, 0.1000f, "HC" );
        	residue.atom[ 17 ] = new Atom( "HD22", 2.973f, -2.998f, 2.436f, 0.1000f, "HC" );
        	residue.atom[ 18 ] = new Atom( "HD23", 1.534f, -2.443f, 3.325f, 0.1000f, "HC" );
        	addBackboneBonds();
        	residue.atom[ 4 ].addBond( 5 );
        	residue.atom[ 2 ].addBond( 6 );
        	residue.atom[ 6 ].addBond( 7, 8, 9 );
        	residue.atom[ 9 ].addBond( 10, 11, 15 );
        	residue.atom[ 11 ].addBond( 12, 13, 14 );
        	residue.atom[ 15 ].addBond( 16, 17, 18 );
        }
        else if ( residueName.equalsIgnoreCase( "Llys" ) ){
        	createResidue( "LYS", 22 );
        	residue.atom[ 0 ] = new Atom( "N", 2.044f, 1.352f, 0.000f, -0.3479f, "N" );
        	residue.atom[ 1 ] = new Atom( "H", 1.721f, 1.837f, 0.824f, 0.2747f, "H" );
        	residue.atom[ 2 ] = new Atom( "CA", 1.522f, 0.000f, 0.000f, -0.2400f, "CT" );
        	residue.atom[ 3 ] = new Atom( "HA", 1.896f, -0.481f, -0.904f, 0.1426f, "H1" );
        	residue.atom[ 4 ] = new Atom( "C", 0.000f, 0.000f, 0.000f, 0.7341f, "C" );
        	residue.atom[ 5 ] = new Atom( "O", -0.624f, 1.060f, 0.000f, -0.5894f, "O" );
        	residue.atom[ 6 ] = new Atom( "CB", 1.988f, -0.769f, 1.233f, -0.0094f, "CT" );
        	residue.atom[ 7 ] = new Atom( "HB1", 3.077f, -0.816f, 1.241f, 0.0362f, "HC" );
        	residue.atom[ 8 ] = new Atom( "HB2", 1.638f, -0.260f, 2.131f, 0.0362f, "HC" );
        	residue.atom[ 9 ] = new Atom( "CG", 1.418f, -2.182f, 1.194f, 0.0187f, "CT" );
        	residue.atom[ 10 ] = new Atom( "HG1", 0.329f, -2.136f, 1.186f, 0.0103f, "HC" );
        	residue.atom[ 11 ] = new Atom( "HG2", 1.767f, -2.692f, 0.296f, 0.0103f, "HC" );
        	residue.atom[ 12 ] = new Atom( "CD", 1.884f, -2.951f, 2.427f, -0.0479f, "CT" );
        	residue.atom[ 13 ] = new Atom( "HD1", 2.973f, -2.998f, 2.436f, 0.0621f, "HC" );
        	residue.atom[ 14 ] = new Atom( "HD2", 1.534f, -2.443f, 3.325f, 0.0621f, "HC" );
        	residue.atom[ 15 ] = new Atom( "CE", 1.314f, -4.366f, 2.389f, -0.0143f, "CT" );
        	residue.atom[ 16 ] = new Atom( "HE1", 0.225f, -4.318f, 2.381f, 0.1135f, "HP" );
        	residue.atom[ 17 ] = new Atom( "HE2", 1.663f, -4.874f, 1.491f, 0.1135f, "HP" );
        	residue.atom[ 18 ] = new Atom( "NZ", 1.763f, -5.107f, 3.577f, -0.3854f, "N3" );
        	residue.atom[ 19 ] = new Atom( "HZ1", 1.385f, -6.042f, 3.552f, 0.3400f, "H" );
        	residue.atom[ 20 ] = new Atom( "HZ2", 2.772f, -5.150f, 3.585f, 0.3400f, "H" );
        	residue.atom[ 21 ] = new Atom( "HZ3", 1.440f, -4.635f, 4.409f, 0.3400f, "H" );
        	addBackboneBonds();
        	residue.atom[ 4 ].addBond( 5 );
        	residue.atom[ 2 ].addBond( 6 );
        	residue.atom[ 6 ].addBond( 7, 8, 9 );
        	residue.atom[ 9 ].addBond( 10, 11, 12 );
        	residue.atom[ 12 ].addBond( 13, 14, 15 );
        	residue.atom[ 15 ].addBond( 16, 17, 18 );
        	residue.atom[ 18 ].addBond( 19, 20, 21 );
        }
        else if ( residueName.equalsIgnoreCase( "Lmet" ) ){
        	createResidue( "MET", 17 );
        	residue.atom[ 0 ] = new Atom( "N", 2.044f, 1.352f, 0.000f, -0.4157f, "N" );
        	residue.atom[ 1 ] = new Atom( "H", 1.721f, 1.837f, 0.824f, 0.2719f, "H" );
        	residue.atom[ 2 ] = new Atom( "CA", 1.522f, 0.000f, 0.000f, -0.0237f, "CT" );
        	residue.atom[ 3 ] = new Atom( "HA", 1.896f, -0.481f, -0.904f, 0.0880f, "H1" );
        	residue.atom[ 4 ] = new Atom( "C", 0.000f, 0.000f, 0.000f, 0.5973f, "C" );
        	residue.atom[ 5 ] = new Atom( "O", -0.624f, 1.060f, 0.000f, -0.5679f, "O" );
        	residue.atom[ 6 ] = new Atom( "CB", 1.988f, -0.769f, 1.233f, 0.0342f, "CT" );
        	residue.atom[ 7 ] = new Atom( "HB1", 3.077f, -0.816f, 1.241f, 0.0241f, "HC" );
        	residue.atom[ 8 ] = new Atom( "HB2", 1.638f, -0.260f, 2.131f, 0.0241f, "HC" );
        	residue.atom[ 9 ] = new Atom( "CG", 1.418f, -2.182f, 1.194f, 0.0018f, "CG" );
        	residue.atom[ 10 ] = new Atom( "HG1", 0.329f, -2.135f, 1.179f, 0.0440f, "H1" );
        	residue.atom[ 11 ] = new Atom( "HG2", 1.775f, -2.688f, 0.297f, 0.0440f, "H1" );
        	residue.atom[ 12 ] = new Atom( "SD", 1.962f, -3.109f, 2.652f, -0.2737f, "S" );
        	residue.atom[ 13 ] = new Atom( "CE", 1.167f, -4.670f, 2.341f, -0.0536f, "CT" );
        	residue.atom[ 14 ] = new Atom( "HE1", 1.399f, -5.364f, 3.149f, 0.0684f, "H1" );
        	residue.atom[ 15 ] = new Atom( "HE2", 0.088f, -4.523f, 2.287f, 0.0684f, "H1" );
        	residue.atom[ 16 ] = new Atom( "HE3", 1.525f, -5.079f, 1.396f, 0.0684f, "H1" );
        	addBackboneBonds();
        	residue.atom[ 4 ].addBond( 5 );
        	residue.atom[ 2 ].addBond( 6 );
        	residue.atom[ 6 ].addBond( 7, 8, 9 );
        	residue.atom[ 9 ].addBond( 10, 11, 12 );
        	residue.atom[ 12 ].addBond( 13 );
        	residue.atom[ 13 ].addBond( 14, 15, 16 );
        }
        else if ( residueName.equalsIgnoreCase( "Lnme" ) ){
        	createResidue( "NME", 6 );
        	residue.atom[ 0 ] = new Atom( "N", -1.227f, 0.728f, 2.125f, -0.463f );
        	residue.atom[ 1 ] = new Atom( "H", -1.124f, -0.261f, 1.947f, 0.252f );
        	residue.atom[ 2 ] = new Atom( "CA", -1.918f, 1.159f, 3.323f, 0.03502f );
        	residue.atom[ 3 ] = new Atom( "HA1", -1.939f, 2.249f, 3.358f, 0.0586f );
        	residue.atom[ 4 ] = new Atom( "HA2", -2.939f, 0.777f, 3.311f, 0.0586f );
        	residue.atom[ 5 ] = new Atom( "HA3", -1.398f, 0.777f, 4.201f, 0.0586f );
        	addBackboneBonds();
        	residue.atom[ 2 ].addBond( 5 );
        }
        else if ( residueName.equalsIgnoreCase( "Lphe" ) ){
        	createResidue( "PHE", 20 );
        	residue.atom[ 0 ] = new Atom( "N", 2.044f, 1.352f, 0.000f, -0.4157f, "N" );
        	residue.atom[ 1 ] = new Atom( "H", 1.721f, 1.837f, 0.824f, 0.2719f, "H" );
        	residue.atom[ 2 ] = new Atom( "CA", 1.522f, 0.000f, 0.000f, -0.0024f, "CT" );
        	residue.atom[ 3 ] = new Atom( "HA", 1.896f, -0.481f, -0.904f, 0.0978f, "H1" );
        	residue.atom[ 4 ] = new Atom( "C", 0.000f, 0.000f, 0.000f, 0.5973f, "C" );
        	residue.atom[ 5 ] = new Atom( "O", -0.624f, 1.060f, 0.000f, -0.5679f, "O" );
        	residue.atom[ 6 ] = new Atom( "CB", 1.988f, -0.769f, 1.233f, -0.0343f, "CT" );
        	residue.atom[ 7 ] = new Atom( "HB1", 3.076f, -0.759f, 1.276f, 0.0295f, "HC" );
        	residue.atom[ 8 ] = new Atom( "HB2", 1.547f, -0.251f, 2.084f, 0.0295f, "HC" );
        	residue.atom[ 9 ] = new Atom( "CG", 1.492f, -2.192f, 1.321f, 0.0118f, "CA" );
        	residue.atom[ 10 ] = new Atom( "CD1", 1.856f, -2.993f, 2.410f, -0.1256f, "CA" );
        	residue.atom[ 11 ] = new Atom( "HD1", 2.497f, -2.589f, 3.194f, 0.1330f, "HA" );
        	residue.atom[ 12 ] = new Atom( "CD2", 0.669f, -2.712f, 0.315f, -0.1256f, "CA" );
        	residue.atom[ 13 ] = new Atom( "HD2", 0.385f, -2.087f, -0.533f, 0.1330f, "HA" );
        	residue.atom[ 14 ] = new Atom( "CE1", 1.398f, -4.313f, 2.492f, -0.1704f, "CA" );
        	residue.atom[ 15 ] = new Atom( "HE1", 1.681f, -4.937f, 3.340f, 0.1430f, "HA" );
        	residue.atom[ 16 ] = new Atom( "CE2", 0.210f, -4.031f, 0.397f, -0.1704f, "CA" );
        	residue.atom[ 17 ] = new Atom( "HE2", -0.431f, -4.435f, -0.386f, 0.1430f, "HA" );
        	residue.atom[ 18 ] = new Atom( "CZ", 0.575f, -4.833f, 1.486f, -0.1072f, "CA" );
        	residue.atom[ 19 ] = new Atom( "HZ", 0.217f, -5.860f, 1.550f, 0.1297f, "HA" );
        	addBackboneBonds();
        	residue.atom[ 4 ].addBond( 5 );
        	residue.atom[ 2 ].addBond( 6 );
        	residue.atom[ 6 ].addBond( 7, 8, 9 );
        	residue.atom[ 9 ].addBond( 10, 12 );
        	residue.atom[ 10 ].addBond( 11, 14 );
        	residue.atom[ 12 ].addBond( 13, 16 );
        	residue.atom[ 14 ].addBond( 15, 18 );
        	residue.atom[ 16 ].addBond( 17, 18 );
        	residue.atom[ 18 ].addBond( 19 );
        }
        else if ( residueName.equalsIgnoreCase( "Lpro" ) ){
        	createResidue( "PRO", 14 );
        	residue.atom[ 0 ] = new Atom( "N", 2.068f, 1.413f, 0.000f, -0.2548f, "N" );
        	residue.atom[ 1 ] = new Atom( "CA", 1.523f, 0.000f, 0.000f, -0.0266f, "CT" );
        	residue.atom[ 2 ] = new Atom( "HA", 2.373f, -0.566f, -0.381f, 0.0641f, "H1" );
        	residue.atom[ 3 ] = new Atom( "C", 0.000f, 0.000f, 0.000f, 0.5896f, "C" );
        	residue.atom[ 4 ] = new Atom( "O", -0.623f, 1.059f, -0.001f, -0.5748f, "O" );
        	residue.atom[ 5 ] = new Atom( "CB", 0.361f, 0.011f, 0.950f, -0.0070f, "CT" );
        	residue.atom[ 6 ] = new Atom( "HB1", -0.547f, -0.274f, 0.417f, 0.0253f, "HC" );
        	residue.atom[ 7 ] = new Atom( "HB2", 0.547f, -0.695f, 1.758f, 0.0253f, "HC" );
        	residue.atom[ 8 ] = new Atom( "CG", 0.192f, 1.398f, 1.523f, 0.0189f, "CT" );
        	residue.atom[ 9 ] = new Atom( "HG1", -0.781f, 1.794f, 1.235f, 0.0213f, "HC" );
        	residue.atom[ 10 ] = new Atom( "HG2", 0.257f, 1.396f, 2.611f, 0.0213f, "HC" );
        	residue.atom[ 11 ] = new Atom( "CD", 1.274f, 2.235f, 0.905f, 0.0192f, "CT" );
        	residue.atom[ 12 ] = new Atom( "HD1", 1.916f, 2.636f, 1.689f, 0.0391f, "H1" );
        	residue.atom[ 13 ] = new Atom( "HD2", 0.824f, 3.057f, 0.349f, 0.0391f, "H1" );
        	residue.atom[ 0 ].addBond( 1, 11 );
        	residue.atom[ 1 ].addBond( 2, 3, 5 );
        	residue.atom[ 3 ].addBond( 4 );
        	residue.atom[ 5 ].addBond( 6, 7, 8 );
        	residue.atom[ 8 ].addBond( 9, 10, 11 );
        	residue.atom[ 11 ].addBond( 12, 13 );
        }
        else if ( residueName.equalsIgnoreCase( "Lser" ) ){
        	createResidue( "SER", 11 );
        	residue.atom[ 0 ] = new Atom( "N", 2.044f, 1.352f, 0.000f, -0.4157f, "N" );
        	residue.atom[ 1 ] = new Atom( "H", 1.721f, 1.837f, 0.824f, 0.2719f, "H" );
        	residue.atom[ 2 ] = new Atom( "CA", 1.522f, 0.000f, 0.000f, -0.0249f, "CT" );
        	residue.atom[ 3 ] = new Atom( "HA", 1.896f, -0.481f, -0.904f, 0.0843f, "H1" );
        	residue.atom[ 4 ] = new Atom( "C", 0.000f, 0.000f, 0.000f, 0.5973f, "C" );
        	residue.atom[ 5 ] = new Atom( "O", -0.624f, 1.060f, 0.000f, -0.5679f, "O" );
        	residue.atom[ 6 ] = new Atom( "CB", 1.988f, -0.769f, 1.233f, 0.2117f, "CT" );
        	residue.atom[ 7 ] = new Atom( "HB1", 3.077f, -0.816f, 1.241f, 0.0352f, "H1" );
        	residue.atom[ 8 ] = new Atom( "HB2", 1.638f, -0.260f, 2.131f, 0.0352f, "H1" );
        	residue.atom[ 9 ] = new Atom( "OG", 1.453f, -2.094f, 1.197f, -0.6546f, "OH" );
        	residue.atom[ 10 ] = new Atom( "HG", 1.746f, -2.579f, 1.973f, 0.4275f, "HO" );
        	addBackboneBonds();
        	residue.atom[ 4 ].addBond( 5 );
        	residue.atom[ 2 ].addBond( 6 );
        	residue.atom[ 6 ].addBond( 7, 8, 9 );
        	residue.atom[ 9 ].addBond( 10 );
        }
        else if ( residueName.equalsIgnoreCase( "Lthr" ) ){
        	createResidue( "THR", 14 );
        	residue.atom[ 0 ] = new Atom( "N", 2.044f, 1.352f, 0.000f, -0.4157f, "N" );
        	residue.atom[ 1 ] = new Atom( "H", 1.721f, 1.837f, 0.824f, 0.2719f, "H" );
        	residue.atom[ 2 ] = new Atom( "CA", 1.522f, 0.000f, 0.000f, -0.0389f, "CT" );
        	residue.atom[ 3 ] = new Atom( "HA", 1.896f, -0.481f, -0.904f, 0.1007f, "H1" );
        	residue.atom[ 4 ] = new Atom( "C", 0.000f, 0.000f, 0.000f, 0.5973f, "C" );
        	residue.atom[ 5 ] = new Atom( "O", -0.624f, 1.060f, 0.000f, -0.5679f, "O" );
        	residue.atom[ 6 ] = new Atom( "CB", 1.988f, -0.769f, 1.233f, 0.3654f, "CT" );
        	residue.atom[ 7 ] = new Atom( "HB", 1.580f, -1.780f, 1.205f, 0.0043f, "H1" );
        	residue.atom[ 8 ] = new Atom( "OG1", 1.530f, -0.101f, 2.411f, -0.6761f, "OH" );
        	residue.atom[ 9 ] = new Atom( "HG1", 1.823f, -0.585f, 3.187f, 0.4102f, "HO" );
        	residue.atom[ 10 ] = new Atom( "CG2", 3.510f, -0.835f, 1.244f, -0.2438f, "CT" );
        	residue.atom[ 11 ] = new Atom( "HG21", 3.844f, -1.384f, 2.125f, 0.0642f, "HC" );
        	residue.atom[ 12 ] = new Atom( "HG22", 3.860f, -1.343f, 0.346f, 0.0642f, "HC" );
        	residue.atom[ 13 ] = new Atom( "HG23", 3.918f, 0.177f, 1.271f, 0.0642f, "HC" );
        	addBackboneBonds();
        	residue.atom[ 4 ].addBond( 5 );
        	residue.atom[ 2 ].addBond( 6 );
        	residue.atom[ 6 ].addBond( 7, 8, 10 );
        	residue.atom[ 8 ].addBond( 9 );
        	residue.atom[ 10 ].addBond( 11, 12, 13 );
        }
        else if ( residueName.equalsIgnoreCase( "Ltrp" ) ){
        	createResidue( "TRP", 24 );
        residue.atom[ 0 ] = new Atom( "N", 2.044f, 1.352f, 0.000f, -0.4157f, "N" );
        residue.atom[ 1 ] = new Atom( "H", 1.721f, 1.837f, 0.824f, 0.2719f, "H" );
        residue.atom[ 2 ] = new Atom( "CA", 1.522f, 0.000f, 0.000f, -0.0275f, "CT" );
        residue.atom[ 3 ] = new Atom( "HA", 1.896f, -0.481f, -0.904f, 0.1123f, "H1" );
        residue.atom[ 4 ] = new Atom( "C", 0.000f, 0.000f, 0.000f, 0.5973f, "C" );
        	residue.atom[ 5 ] = new Atom( "O", -0.624f, 1.060f, 0.000f, -0.5679f, "O" );
        residue.atom[ 6 ] = new Atom( "CB", 1.988f, -0.769f, 1.233f, -0.0050f, "CT" );
        residue.atom[ 7 ] = new Atom( "HB1", 3.076f, -0.759f, 1.276f, 0.0339f, "HC" );
        residue.atom[ 8 ] = new Atom( "HB2", 1.547f, -0.251f, 2.084f, 0.0339f, "HC" );
        residue.atom[ 9 ] = new Atom( "CG", 1.492f, -2.192f, 1.321f, -0.1415f, "C*" );
        residue.atom[ 10 ] = new Atom( "CD1", 1.752f, -3.078f, 2.293f, -0.1638f, "CW" );
        residue.atom[ 11 ] = new Atom( "HD1", 2.365f, -2.906f, 3.178f, 0.2062f, "H4" );
        residue.atom[ 12 ] = new Atom( "CD2", 0.669f, -2.868f, 0.412f, 0.1243f, "CB" );
        residue.atom[ 13 ] = new Atom( "NE1", 1.072f, -4.288f, 1.950f, -0.3418f, "NA" );
        residue.atom[ 14 ] = new Atom( "HE1", 1.079f, -5.139f, 2.493f, 0.3412f, "H" );
        residue.atom[ 15 ] = new Atom( "CE2", 0.438f, -4.113f, 0.817f, 0.1380f, "CN" );
        residue.atom[ 16 ] = new Atom( "CE3", 0.103f, -2.412f, -0.785f, -0.2387f, "CA" );
        residue.atom[ 17 ] = new Atom( "HE3", 0.273f, -1.397f, -1.145f, 0.1700f, "HA" );
        residue.atom[ 18 ] = new Atom( "CZ2", -0.350f, -5.037f, 0.120f, -0.2601f, "CA" );
        residue.atom[ 19 ] = new Atom( "HZ2", -0.515f, -6.050f, 0.487f, 0.1572f, "HA" );
        residue.atom[ 20 ] = new Atom( "CZ3", -0.694f, -3.325f, -1.505f, -0.1972f, "CA" );
        residue.atom[ 21 ] = new Atom( "HZ3", -1.150f, -3.005f, -2.442f, 0.1447f, "HA" );
        	residue.atom[ 22 ] = new Atom( "CH2", -0.912f, -4.584f, -1.069f, -0.1134f, "CA" );
        	residue.atom[ 23 ] = new Atom( "HH2", -1.535f, -5.257f, -1.658f, 0.1417f, "HA" );
        	addBackboneBonds();
        residue.atom[ 4 ].addBond( 5 );
        residue.atom[ 2 ].addBond( 6 );
        residue.atom[ 6 ].addBond( 7, 8, 9 );
        residue.atom[ 9 ].addBond( 10, 12 );
        residue.atom[ 10 ].addBond( 11, 13 );
        residue.atom[ 13 ].addBond( 14, 15 );
        residue.atom[ 15 ].addBond( 18 );
        residue.atom[ 18 ].addBond( 19, 22 );
        residue.atom[ 22 ].addBond( 23 );
        residue.atom[ 12 ].addBond( 15, 16 );
        residue.atom[ 16 ].addBond( 17, 20 );
        residue.atom[ 20 ].addBond( 21, 22 );
        }
        else if ( residueName.equalsIgnoreCase( "Ltyr" ) ){
        	createResidue( "TYR", 21 );
        	residue.atom[ 0 ] = new Atom( "N", 2.044f, 1.352f, 0.000f, -0.4157f, "N" );
        	residue.atom[ 1 ] = new Atom( "H", 1.721f, 1.837f, 0.824f, 0.2719f, "H" );
        	residue.atom[ 2 ] = new Atom( "CA", 1.522f, 0.000f, 0.000f, -0.0014f, "CT" );
        	residue.atom[ 3 ] = new Atom( "HA", 1.896f, -0.481f, -0.904f, 0.0876f, "H1" );
        	residue.atom[ 4 ] = new Atom( "C", 0.000f, 0.000f, 0.000f, 0.5973f, "C" );
        	residue.atom[ 5 ] = new Atom( "O", -0.624f, 1.060f, 0.000f, -0.5679f, "O" );
        	residue.atom[ 6 ] = new Atom( "CB", 1.988f, -0.769f, 1.233f, -0.0152f, "CT" );
        	residue.atom[ 7 ] = new Atom( "HB1", 3.077f, -0.816f, 1.241f, 0.0295f, "HC" );
        	residue.atom[ 8 ] = new Atom( "HB2", 1.638f, -0.260f, 2.131f, 0.0295f, "HC" );
        	residue.atom[ 9 ] = new Atom( "CG", 1.423f, -2.168f, 1.195f, -0.0011f, "CA" );
        	residue.atom[ 10 ] = new Atom( "CD1", 1.715f, -3.068f, 2.227f, -0.1906f, "CA" );
        	residue.atom[ 11 ] = new Atom( "HD1", 2.348f, -2.758f, 3.058f, 0.1699f, "HA" );
        	residue.atom[ 12 ] = new Atom( "CD2", 0.609f, -2.567f, 0.128f, -0.1906f, "CA" );
        	residue.atom[ 13 ] = new Atom( "HD2", 0.382f, -1.867f, -0.676f, 0.1699f, "HA" );
        	residue.atom[ 14 ] = new Atom( "CE1", 1.192f, -4.367f, 2.193f, -0.2341f, "CA" );
        	residue.atom[ 15 ] = new Atom( "HE1", 1.418f, -5.067f, 2.996f, 0.1656f, "HA" );
        	residue.atom[ 16 ] = new Atom( "CE2", 0.086f, -3.866f, 0.093f, -0.2341f, "CA" );
        	residue.atom[ 17 ] = new Atom( "HE2", -0.548f, -4.176f, -0.737f, 0.1656f, "HA" );
        	residue.atom[ 18 ] = new Atom( "CZ", 0.378f, -4.766f, 1.126f, 0.3226f, "C" );
        	residue.atom[ 19 ] = new Atom( "OH", -0.131f, -6.026f, 1.092f, -0.5579f, "OH" );
        	residue.atom[ 20 ] = new Atom( "HH", 0.132f, -6.557f, 1.849f, 0.3992f, "HO" );
        	addBackboneBonds();
        	residue.atom[ 4 ].addBond( 5 );
        	residue.atom[ 2 ].addBond( 6 );
        	residue.atom[ 6 ].addBond( 7, 8, 9 );
        	residue.atom[ 9 ].addBond( 10, 12 );
        	residue.atom[ 10 ].addBond( 11, 14 );
        	residue.atom[ 14 ].addBond( 15, 18 );
        	residue.atom[ 18 ].addBond( 19 );
        	residue.atom[ 19 ].addBond( 20 );
        	residue.atom[ 12 ].addBond( 13, 16 );
        	residue.atom[ 16 ].addBond( 17, 18 );
        }
        else if ( residueName.equalsIgnoreCase( "Lval" ) ){
        	createResidue( "VAL", 16 );
        residue.atom[ 0 ] = new Atom( "N", 2.044f, 1.352f, 0.000f, -0.4157f, "N" );
        residue.atom[ 1 ] = new Atom( "H", 1.721f, 1.837f, 0.824f, 0.2719f, "H" );
        residue.atom[ 2 ] = new Atom( "CA", 1.522f, 0.000f, 0.000f, -0.0875f, "CT" );
        residue.atom[ 3 ] = new Atom( "HA", 1.896f, -0.481f, -0.904f, 0.0969f, "H1" );
        residue.atom[ 4 ] = new Atom( "C", 0.000f, 0.000f, 0.000f, 0.5973f, "C" );
        	residue.atom[ 5 ] = new Atom( "O", -0.624f, 1.060f, 0.000f, -0.5679f, "O" );
        residue.atom[ 6 ] = new Atom( "CB", 1.988f, -0.769f, 1.233f, 0.2985f, "CT" );
        residue.atom[ 7 ] = new Atom( "HB", 3.077f, -0.816f, 1.241f, -0.0297f, "HC" );
        residue.atom[ 8 ] = new Atom( "CG1", 1.499f, -0.057f, 2.490f, -0.3192f, "CT" );
        residue.atom[ 9 ] = new Atom( "HG11", 1.832f, -0.606f, 3.370f, 0.0791f, "HC" );
        residue.atom[ 10 ] = new Atom( "HG12", 1.906f, 0.953f, 2.516f, 0.0791f, "HC" );
        residue.atom[ 11 ] = new Atom( "HG13", 0.410f, -0.010f, 2.481f, 0.0791f, "HC" );
        residue.atom[ 12 ] = new Atom( "CG2", 1.418f, -2.182f, 1.194f, -0.3192f, "CT" );
        residue.atom[ 13 ] = new Atom( "HG21", 1.751f, -2.732f, 2.075f, 0.0791f, "HC" );
        residue.atom[ 14 ] = new Atom( "HG22", 0.329f, -2.136f, 1.186f, 0.0791f, "HC" );
        residue.atom[ 15 ] = new Atom( "HG23", 1.767f, -2.692f, 0.296f, 0.0791f, "HC" );
        	addBackboneBonds();
        residue.atom[ 4 ].addBond( 5 );
        residue.atom[ 2 ].addBond( 6 );
        residue.atom[ 6 ].addBond( 7, 8, 12 );
        residue.atom[ 8 ].addBond( 9, 10, 11 );
        residue.atom[ 12 ].addBond( 13, 14, 15 );
        }
        // If they're looking for a D amino acid, print an error message and return
        // But keep the templates in the file below so that we might add D amino
        //  acids one day
        else if ((residueName.charAt(0) == 'D') || (residueName.charAt(0) == 'd'))
        	{
        		System.out.println("*** Warning: Amber96PolyPeptideResidue does not yet have D amino acids implemented");
        		return new Residue();
        	}
        else if ( residueName.equalsIgnoreCase( "Dace" ) ){
        	createResidue( "ACE", 6 );
        	residue.atom[ 0 ] = new Atom( "H1", 0.000f, 0.000f, 0.000f );
        	residue.atom[ 1 ] = new Atom( "CH3", 0.000f, 1.090f, 0.000f );
        	residue.atom[ 2 ] = new Atom( "H2", 1.028f, 1.453f, 0.000f );
        	residue.atom[ 3 ] = new Atom( "H3", -0.514f, 1.453f, -0.890f );
        	residue.atom[ 4 ] = new Atom( "C", -0.721f, 1.600f, 1.249f, 0.616f );
        	residue.atom[ 5 ] = new Atom( "O", -0.839f, 2.806f, 1.453f, -0.504f );
        	residue.atom[ 0 ].addBond( 1 );
        	residue.atom[ 1 ].addBond( 2, 3 );
        	residue.atom[ 1 ].addBond( 4 );
        	residue.atom[ 4 ].addBond( 5 );
        }
        else if ( residueName.equalsIgnoreCase( "Dala" ) ){
        	createResidue( "ALA", 10 );
        	residue.atom[ 0 ] = new Atom( "N", 2.044f, 1.352f, 0.000f, -0.463f );
        	residue.atom[ 1 ] = new Atom( "H", 1.721f, 1.837f, 0.824f, 0.252f );
        	residue.atom[ 2 ] = new Atom( "CA", 1.522f, 0.000f, 0.000f, 0.035f );
        	residue.atom[ 3 ] = new Atom( "HA", 1.896f, -0.481f, 0.904f, 0.048f );
        	residue.atom[ 4 ] = new Atom( "C", 0.000f, 0.000f, 0.000f, 0.616f );
        	residue.atom[ 5 ] = new Atom( "O", -0.624f, 1.060f, 0.000f, -0.504f );
        	residue.atom[ 6 ] = new Atom( "CB", 1.988f, -0.769f, -1.232f, -0.098f );
        	residue.atom[ 7 ] = new Atom( "HB1", 1.580f, -1.779f, -1.205f, 0.036f );
        	residue.atom[ 8 ] = new Atom( "HB2", 1.638f, -0.260f, -2.131f, 0.036f );
        	residue.atom[ 9 ] = new Atom( "HB3", 3.077f, -0.816f, -1.241f, 0.036f );
        	addBackboneBonds();
        	residue.atom[ 4 ].addBond( 5 );
        	residue.atom[ 2 ].addBond( 6 );
        	residue.atom[ 6 ].addBond( 7, 8, 9 );
        }
        else if ( residueName.equalsIgnoreCase( "Darg" ) ){
        	createResidue( "ARG", 24 );
        	residue.atom[ 0 ] = new Atom( "N", 2.044f, 1.352f, 0.000f, -0.463f );
        	residue.atom[ 1 ] = new Atom( "H", 1.721f, 1.837f, 0.824f, 0.252f );
        	residue.atom[ 2 ] = new Atom( "CA", 1.522f, 0.000f, 0.000f, 0.035f );
        	residue.atom[ 3 ] = new Atom( "HA", 1.896f, -0.481f, 0.904f, 0.048f );
        	residue.atom[ 4 ] = new Atom( "C", 0.000f, 0.000f, 0.000f, 0.616f );
        	residue.atom[ 5 ] = new Atom( "O", -0.624f, 1.060f, 0.000f, -0.504f );
        	residue.atom[ 6 ] = new Atom( "CB", 1.988f, -0.769f, -1.232f, -0.080f );
        	residue.atom[ 7 ] = new Atom( "HB1", 1.638f, -0.260f, -2.131f, 0.056f );
        	residue.atom[ 8 ] = new Atom( "HB2", 3.077f, -0.816f, -1.241f, 0.056f );
        	residue.atom[ 9 ] = new Atom( "CG", 1.418f, -2.183f, -1.195f, -0.103f );
        	residue.atom[ 10 ] = new Atom( "HG1", 1.767f, -2.692f, -0.296f, 0.074f );
        	residue.atom[ 11 ] = new Atom( "HG2", 0.328f, -2.136f, -1.186f, 0.074f );
        	residue.atom[ 12 ] = new Atom( "CD", 1.884f, -2.951f, -2.427f, -0.228f );
        	residue.atom[ 13 ] = new Atom( "HD1", 1.549f, -2.433f, -3.326f, 0.133f );
        	residue.atom[ 14 ] = new Atom( "HD2", 2.972f, -3.004f, -2.410f, 0.133f );
        	residue.atom[ 15 ] = new Atom( "NE", 1.348f, -4.332f, -2.424f, -0.324f );
        	residue.atom[ 16 ] = new Atom( "HE", 0.761f, -4.619f, -1.655f, 0.269f );
        	residue.atom[ 17 ] = new Atom( "CZ", 1.606f, -5.210f, -3.390f, 0.760f );
        	residue.atom[ 18 ] = new Atom( "NH1", 2.371f, -4.905f, -4.434f, -0.624f );
        	residue.atom[ 19 ] = new Atom( "HH11", 2.774f, -3.982f, -4.509f, 0.361f );
        	residue.atom[ 20 ] = new Atom( "HH12", 2.545f, -5.597f, -5.148f, 0.361f );
        	residue.atom[ 21 ] = new Atom( "NH2", 1.074f, -6.424f, -3.287f, -0.624f );
        	residue.atom[ 22 ] = new Atom( "HH21", 0.494f, -6.652f, -2.492f, 0.361f );
        	residue.atom[ 23 ] = new Atom( "HH22", 1.252f, -7.112f, -4.004f, 0.361f );
        	addBackboneBonds();
        	residue.atom[ 4 ].addBond( 5 );
        	residue.atom[ 2 ].addBond( 6 );
        	residue.atom[ 6 ].addBond( 7, 8, 9 );
        	residue.atom[ 9 ].addBond( 10, 11, 12 );
        	residue.atom[ 12 ].addBond( 13, 14, 15 );
        	residue.atom[ 15 ].addBond( 16, 17 );
        	residue.atom[ 17 ].addBond( 18, 21 );
        	residue.atom[ 18 ].addBond( 19, 20 );
        	residue.atom[ 21 ].addBond( 22, 23 );
        }
        else if ( residueName.equalsIgnoreCase( "Dasn" ) ){
        	createResidue( "ASN", 14 );
        	residue.atom[ 0 ] = new Atom( "N", 2.044f, 1.352f, 0.000f, -0.463f );
        	residue.atom[ 1 ] = new Atom( "H", 1.721f, 1.837f, 0.824f, 0.252f );
        	residue.atom[ 2 ] = new Atom( "CA", 1.522f, 0.000f, 0.000f, 0.035f );
        	residue.atom[ 3 ] = new Atom( "HA", 1.896f, -0.481f, 0.904f, 0.048f );
        	residue.atom[ 4 ] = new Atom( "C", 0.000f, 0.000f, 0.000f, 0.616f );
        	residue.atom[ 5 ] = new Atom( "O", -0.624f, 1.060f, 0.000f, -0.504f );
        	residue.atom[ 6 ] = new Atom( "CB", 1.988f, -0.769f, -1.232f, -0.086f );
        	residue.atom[ 7 ] = new Atom( "HB1", 1.653f, -0.251f, -2.131f, 0.038f );
        	residue.atom[ 8 ] = new Atom( "HB2", 3.076f, -0.821f, -1.214f, 0.038f );
        	residue.atom[ 9 ] = new Atom( "CG", 1.438f, -2.188f, -1.232f, 0.675f );
        	residue.atom[ 10 ] = new Atom( "OD1", 0.720f, -2.579f, -0.316f, -0.470f );
        	residue.atom[ 11 ] = new Atom( "ND2", 1.780f, -2.961f, -2.265f, -0.867f );
        	residue.atom[ 12 ] = new Atom( "HD21", 2.374f, -2.591f, -2.994f, 0.344f );
        	residue.atom[ 13 ] = new Atom( "HD22", 1.443f, -3.912f, -2.315f, 0.344f );
        	addBackboneBonds();
        	residue.atom[ 4 ].addBond( 5 );
        	residue.atom[ 2 ].addBond( 6 );
        	residue.atom[ 6 ].addBond( 7, 8, 9 );
        	residue.atom[ 9 ].addBond( 10, 11 );
        	residue.atom[ 11 ].addBond( 12, 13 );
        }
        else if ( residueName.equalsIgnoreCase( "Dasp" ) ){
        	createResidue( "ASP", 12 );
        	residue.atom[ 0 ] = new Atom( "N", 2.044f, 1.352f, 0.000f, -0.463f );
        	residue.atom[ 1 ] = new Atom( "H", 1.721f, 1.837f, 0.824f, 0.252f );
        	residue.atom[ 2 ] = new Atom( "CA", 1.522f, 0.000f, 0.000f, 0.035f );
        	residue.atom[ 3 ] = new Atom( "HA", 1.896f, -0.481f, 0.904f, 0.048f );
        	residue.atom[ 4 ] = new Atom( "C", 0.000f, 0.000f, 0.000f, 0.616f );
        	residue.atom[ 5 ] = new Atom( "O", -0.624f, 1.060f, 0.000f, -0.504f );
        	residue.atom[ 6 ] = new Atom( "CB", 1.988f, -0.769f, -1.232f, -0.398f );
        	residue.atom[ 7 ] = new Atom( "HB1", 1.638f, -0.260f, -2.131f, 0.071f );
        	residue.atom[ 8 ] = new Atom( "HB2", 3.077f, -0.816f, -1.241f, 0.071f );
        	residue.atom[ 9 ] = new Atom( "CG", 1.417f, -2.184f, -1.195f, 0.714f );
        	residue.atom[ 10 ] = new Atom( "OD1", 0.297f, -2.369f, -1.741f, -0.721f );
        	residue.atom[ 11 ] = new Atom( "OD2", 2.107f, -3.069f, -0.620f, -0.721f );
        	addBackboneBonds();
        	residue.atom[ 4 ].addBond( 5 );
        	residue.atom[ 2 ].addBond( 6 );
        	residue.atom[ 6 ].addBond( 7, 8, 9 );
        	residue.atom[ 9 ].addBond( 10, 11 );
        }
        else if ( residueName.equalsIgnoreCase( "Dcys" ) ){
        	createResidue( "CYS", 11 );
        	residue.atom[ 0 ] = new Atom( "N", 2.044f, 1.352f, 0.000f, -0.463f );
        	residue.atom[ 1 ] = new Atom( "H", 1.721f, 1.837f, 0.824f, 0.252f );
        	residue.atom[ 2 ] = new Atom( "CA", 1.522f, 0.000f, 0.000f, 0.035f );
        	residue.atom[ 3 ] = new Atom( "HA", 1.896f, -0.481f, 0.904f, 0.048f );
        	residue.atom[ 4 ] = new Atom( "C", 0.000f, 0.000f, 0.000f, 0.616f );
        	residue.atom[ 5 ] = new Atom( "O", -0.624f, 1.060f, 0.000f, -0.504f );
        	residue.atom[ 6 ] = new Atom( "CB", 1.988f, -0.769f, -1.232f, -0.060f );
        	residue.atom[ 7 ] = new Atom( "HB1", 1.707f, -0.219f, -2.130f, 0.038f );
        	residue.atom[ 8 ] = new Atom( "HB2", 3.069f, -0.846f, -1.122f, 0.038f );
        	residue.atom[ 9 ] = new Atom( "SG", 1.409f, -2.479f, -1.365f, -0.135f );
        	residue.atom[ 10 ] = new Atom( "HG", 1.889f, -3.023f, -2.481f, 0.135f );
        	addBackboneBonds();
        	residue.atom[ 4 ].addBond( 5 );
        	residue.atom[ 2 ].addBond( 6 );
        	residue.atom[ 6 ].addBond( 7, 8, 9 );
        	residue.atom[ 9 ].addBond( 10 );
        }
        else if ( residueName.equalsIgnoreCase( "Dcyx" ) ){
        	createResidue( "CYX", 10 );
        	residue.atom[ 0 ] = new Atom( "N", 2.044f, 1.352f, 0.000f, -0.463f );
        	residue.atom[ 1 ] = new Atom( "H", 1.721f, 1.837f, 0.824f, 0.252f );
        	residue.atom[ 2 ] = new Atom( "CA", 1.522f, 0.000f, 0.000f, 0.035f );
        	residue.atom[ 3 ] = new Atom( "HA", 1.896f, -0.481f, 0.904f, 0.048f );
        	residue.atom[ 4 ] = new Atom( "C", 0.000f, 0.000f, 0.000f, 0.616f );
        	residue.atom[ 5 ] = new Atom( "O", -0.624f, 1.060f, 0.000f, -0.504f );
        	residue.atom[ 6 ] = new Atom( "CB", 1.988f, -0.769f, -1.232f, -0.098f );
        	residue.atom[ 7 ] = new Atom( "HB1", 1.707f, -0.219f, -2.130f, 0.0495f );
        	residue.atom[ 8 ] = new Atom( "HB2", 3.069f, -0.846f, -1.122f, 0.0495f );
        	residue.atom[ 9 ] = new Atom( "SG", 1.409f, -2.479f, -1.365f, 0.015f );
        	addBackboneBonds();
        	residue.atom[ 4 ].addBond( 5 );
        	residue.atom[ 2 ].addBond( 6 );
        	residue.atom[ 6 ].addBond( 7, 8, 9 );
        }
        else if ( residueName.equalsIgnoreCase( "Dgln" ) ){
        	createResidue( "GLN", 17 );
        	residue.atom[ 0 ] = new Atom( "N", 2.044f, 1.352f, 0.000f, -0.463f );
        	residue.atom[ 1 ] = new Atom( "H", 1.721f, 1.837f, 0.824f, 0.252f );
        	residue.atom[ 2 ] = new Atom( "CA", 1.522f, 0.000f, 0.000f, 0.035f );
        	residue.atom[ 3 ] = new Atom( "HA", 1.896f, -0.481f, 0.904f, 0.048f );
        	residue.atom[ 4 ] = new Atom( "C", 0.000f, 0.000f, 0.000f, 0.616f );
        	residue.atom[ 5 ] = new Atom( "O", -0.624f, 1.060f, 0.000f, -0.504f );
        	residue.atom[ 6 ] = new Atom( "CB", 1.988f, -0.769f, -1.232f, -0.098f );
        	residue.atom[ 7 ] = new Atom( "HB1", 1.638f, -0.260f, -2.131f, 0.038f );
        	residue.atom[ 8 ] = new Atom( "HB2", 3.077f, -0.816f, -1.241f, 0.038f );
        	residue.atom[ 9 ] = new Atom( "CG", 1.418f, -2.183f, -1.195f, -0.102f );
        	residue.atom[ 10 ] = new Atom( "HG1", 1.752f, -2.685f, -0.287f, 0.057f );
        	residue.atom[ 11 ] = new Atom( "HG2", 0.330f, -2.115f, -1.205f, 0.057f );
        	residue.atom[ 12 ] = new Atom( "CD", 1.861f, -2.984f, -2.410f, 0.675f );
        	residue.atom[ 13 ] = new Atom( "OE1", 2.585f, -2.476f, -3.263f, -0.470f );
        	residue.atom[ 14 ] = new Atom( "NE2", 1.422f, -4.243f, -2.489f, -0.867f );
        	residue.atom[ 15 ] = new Atom( "HE21", 0.828f, -4.613f, -1.760f, 0.344f );
        	residue.atom[ 16 ] = new Atom( "HE22", 1.687f, -4.819f, -3.275f, 0.344f );
        	addBackboneBonds();
        	residue.atom[ 4 ].addBond( 5 );
        	residue.atom[ 2 ].addBond( 6 );
        	residue.atom[ 6 ].addBond( 7, 8, 9 );
        	residue.atom[ 9 ].addBond( 10, 11, 12 );
        	residue.atom[ 12 ].addBond( 13, 14 );
        	residue.atom[ 14 ].addBond( 15, 16 );
        }
        else if ( residueName.equalsIgnoreCase( "Dglu" ) ){
        	createResidue( "GLU", 15 );
        	residue.atom[ 0 ] = new Atom( "N", 2.044f, 1.352f, 0.000f, -0.463f );
        	residue.atom[ 1 ] = new Atom( "H", 1.721f, 1.837f, 0.824f, 0.252f );
        	residue.atom[ 2 ] = new Atom( "CA", 1.522f, 0.000f, 0.000f, 0.035f );
        	residue.atom[ 3 ] = new Atom( "HA", 1.896f, -0.481f, 0.904f, 0.048f );
        	residue.atom[ 4 ] = new Atom( "C", 0.000f, 0.000f, 0.000f, 0.616f );
        	residue.atom[ 5 ] = new Atom( "O", -0.624f, 1.060f, 0.000f, -0.504f );
        	residue.atom[ 6 ] = new Atom( "CB", 1.988f, -0.769f, -1.232f, -0.184f );
        	residue.atom[ 7 ] = new Atom( "HB1", 1.638f, -0.260f, -2.131f, 0.092f );
        	residue.atom[ 8 ] = new Atom( "HB2", 3.077f, -0.816f, -1.241f, 0.092f );
        	residue.atom[ 9 ] = new Atom( "CG", 1.423f, -2.168f, -1.195f, -0.398f );
        	residue.atom[ 10 ] = new Atom( "HG1", 1.772f, -2.677f, -0.296f, 0.071f );
        	residue.atom[ 11 ] = new Atom( "HG2", 0.334f, -2.122f, -1.186f, 0.071f );
        	residue.atom[ 12 ] = new Atom( "CD", 1.890f, -2.938f, -2.429f, 0.714f );
        	residue.atom[ 13 ] = new Atom( "OE1", 2.971f, -3.579f, -2.334f, -0.721f );
        	residue.atom[ 14 ] = new Atom( "OE2", 1.160f, -2.879f, -3.455f, -0.721f );
        	addBackboneBonds();
        	residue.atom[ 4 ].addBond( 5 );
        	residue.atom[ 2 ].addBond( 6 );
        	residue.atom[ 6 ].addBond( 7, 8, 9 );
        	residue.atom[ 9 ].addBond( 10, 11, 12 );
        	residue.atom[ 12 ].addBond( 13, 14 );
        }
        else if ( residueName.equalsIgnoreCase( "Dgly" ) ){
        	createResidue( "GLY", 7 );
        	residue.atom[ 0 ] = new Atom( "N", 2.027f, 1.358f, 0.000f, -0.463f );
        	residue.atom[ 1 ] = new Atom( "H", 1.697f, 1.839f, 0.824f, 0.252f );
        	residue.atom[ 2 ] = new Atom( "CA", 1.522f, 0.000f, 0.000f, 0.035f );
        	residue.atom[ 3 ] = new Atom( "HA2", 1.886f, -0.523f, -0.885f, 0.032f );
        	residue.atom[ 4 ] = new Atom( "HA3", 1.874f, -0.506f, 0.899f, 0.032f );
        	residue.atom[ 5 ] = new Atom( "C", 0.000f, 0.000f, 0.000f, 0.616f );
        	residue.atom[ 6 ] = new Atom( "O", -0.624f, 1.058f, 0.000f, -0.504f );
        	addBackboneBonds();
        	residue.atom[ 2 ].addBond( 5 );
        	residue.atom[ 5 ].addBond( 6 );
        }
        else if ( residueName.equalsIgnoreCase( "Dhid" ) ){
        	createResidue( "HID", 17 );
        	residue.atom[ 0 ] = new Atom( "N", 2.044f, 1.352f, 0.000f, -0.463f );
        	residue.atom[ 1 ] = new Atom( "H", 1.721f, 1.837f, 0.824f, 0.252f );
        	residue.atom[ 2 ] = new Atom( "CA", 1.522f, 0.000f, 0.000f, 0.035f );
        	residue.atom[ 3 ] = new Atom( "HA", 1.896f, -0.481f, 0.904f, 0.048f );
        	residue.atom[ 4 ] = new Atom( "C", 0.000f, 0.000f, 0.000f, 0.616f );
        	residue.atom[ 5 ] = new Atom( "O", -0.624f, 1.060f, 0.000f, -0.504f );
        	residue.atom[ 6 ] = new Atom( "CB", 1.988f, -0.769f, -1.232f, -0.098f );
        	residue.atom[ 7 ] = new Atom( "HB1", 1.695f, -0.225f, -2.131f, 0.038f );
        	residue.atom[ 8 ] = new Atom( "HB2", 3.071f, -0.840f, -1.141f, 0.038f );
        	residue.atom[ 9 ] = new Atom( "CG", 1.492f, -2.192f, -1.321f, -0.032f );
        	residue.atom[ 10 ] = new Atom( "ND1", 1.828f, -3.024f, -2.383f, -0.146f );
        	residue.atom[ 11 ] = new Atom( "HD1", 2.411f, -2.776f, -3.170f, 0.228f );
        	residue.atom[ 12 ] = new Atom( "CD2", 0.698f, -2.921f, -0.492f, 0.195f );
        	residue.atom[ 13 ] = new Atom( "HD2", 0.306f, -2.472f, 0.421f, 0.018f );
        	residue.atom[ 14 ] = new Atom( "CE1", 1.251f, -4.195f, -2.182f, 0.241f );
        	residue.atom[ 15 ] = new Atom( "HE1", 1.387f, -4.999f, -2.905f, 0.036f );
        	residue.atom[ 16 ] = new Atom( "NE2", 0.575f, -4.150f, -1.061f, -0.502f );
        	addBackboneBonds();
        	residue.atom[ 4 ].addBond( 5 );
        	residue.atom[ 2 ].addBond( 6 );
        	residue.atom[ 6 ].addBond( 7, 8, 9 );
        	residue.atom[ 9 ].addBond( 10, 12 );
        	residue.atom[ 10 ].addBond( 11, 14 );
        	residue.atom[ 12 ].addBond( 13, 16 );
        	residue.atom[ 14 ].addBond( 15, 16 );
        }
        else if ( residueName.equalsIgnoreCase( "Dhie" ) ){
        	createResidue( "HIE", 17 );
        	residue.atom[ 0 ] = new Atom( "N", 2.044f, 1.352f, 0.000f, -0.463f );
        	residue.atom[ 1 ] = new Atom( "H", 1.721f, 1.837f, 0.824f, 0.252f );
        	residue.atom[ 2 ] = new Atom( "CA", 1.522f, 0.000f, 0.000f, 0.035f );
        	residue.atom[ 3 ] = new Atom( "HA", 1.896f, -0.481f, 0.904f, 0.048f );
        	residue.atom[ 4 ] = new Atom( "C", 0.000f, 0.000f, 0.000f, 0.616f );
        	residue.atom[ 5 ] = new Atom( "O", -0.624f, 1.060f, 0.000f, -0.504f );
        	residue.atom[ 6 ] = new Atom( "CB", 1.988f, -0.769f, -1.232f, -0.098f );
        	residue.atom[ 7 ] = new Atom( "HB1", 1.695f, -0.225f, -2.131f, 0.038f );
        	residue.atom[ 8 ] = new Atom( "HB2", 3.071f, -0.840f, -1.141f, 0.038f );
        	residue.atom[ 9 ] = new Atom( "CG", 1.492f, -2.192f, -1.321f, 0.251f );
        	residue.atom[ 10 ] = new Atom( "ND1", 1.828f, -3.024f, -2.383f,  -0.502f );
        	residue.atom[ 11 ] = new Atom( "CD2", 0.698f, -2.921f, -0.492f, -0.184f );
        	residue.atom[ 12 ] = new Atom( "HD2", 0.306f, -2.472f, 0.421f, 0.114f );
        	residue.atom[ 13 ] = new Atom( "CE1", 1.251f, -4.195f, -2.182f, 0.241f );
        	residue.atom[ 14 ] = new Atom( "HE1", 1.387f, -4.999f, -2.905f, 0.036f );
        	residue.atom[ 15 ] = new Atom( "NE2", 0.575f, -4.150f, -1.061f, -0.146f );
        	residue.atom[ 16 ] = new Atom( "HE2", 0.041f, -4.917f, -0.677f, 0.228f );
        	addBackboneBonds();
        	residue.atom[ 4 ].addBond( 5 );
        	residue.atom[ 2 ].addBond( 6 );
        	residue.atom[ 6 ].addBond( 7, 8, 9 );
        	residue.atom[ 9 ].addBond( 10, 11 );
        	residue.atom[ 10 ].addBond( 13 );
        	residue.atom[ 11 ].addBond( 12, 15 );
        	residue.atom[ 13 ].addBond( 14, 15 );
        	residue.atom[ 15 ].addBond( 16 );
        }
        else if ( residueName.equalsIgnoreCase( "Dhip" ) ){
        	createResidue( "HIP", 18 );
        	residue.atom[ 0 ] = new Atom( "N", 2.044f, 1.352f, 0.000f, -0.463f );
        	residue.atom[ 1 ] = new Atom( "H", 1.721f, 1.837f, 0.824f, 0.252f );
        	residue.atom[ 2 ] = new Atom( "CA", 1.522f, 0.000f, 0.000f, 0.035f );
        	residue.atom[ 3 ] = new Atom( "HA", 1.896f, -0.481f, 0.904f, 0.048f );
        	residue.atom[ 4 ] = new Atom( "C", 0.000f, 0.000f, 0.000f, 0.616f );
        	residue.atom[ 5 ] = new Atom( "O", -0.624f, 1.060f, 0.000f, -0.504f );
        	residue.atom[ 6 ] = new Atom( "CB", 1.988f, -0.769f, -1.232f, -0.098f );
        	residue.atom[ 7 ] = new Atom( "HB1", 1.695f, -0.225f, -2.131f, 0.086f );
        	residue.atom[ 8 ] = new Atom( "HB2", 3.071f, -0.840f, -1.141f, 0.086f );
        	residue.atom[ 9 ] = new Atom( "CG", 1.492f, -2.192f, -1.321f, 0.058f );
        	residue.atom[ 10 ] = new Atom( "ND1", 1.828f, -3.024f, -2.383f, -0.058f );
        	residue.atom[ 11 ] = new Atom( "HD1", 2.411f, -2.776f, -3.170f, 0.306f );
        	residue.atom[ 12 ] = new Atom( "CD2", 0.698f, -2.921f, -0.492f, -0.037f );
        	residue.atom[ 13 ] = new Atom( "HD2", 0.306f, -2.472f, 0.421f, 0.153f );
        	residue.atom[ 14 ] = new Atom( "CE1", 1.251f, -4.195f, -2.182f, 0.114f );
        	residue.atom[ 15 ] = new Atom( "HE1", 1.387f, -4.999f, -2.905f, 0.158f );
        	residue.atom[ 16 ] = new Atom( "NE2", 0.575f, -4.150f, -1.061f, -0.058f );
        	residue.atom[ 17 ] = new Atom( "HE2", 0.041f, -4.917f, -0.677f, 0.306f );
        	addBackboneBonds();
        	residue.atom[ 4 ].addBond( 5 );
        	residue.atom[ 2 ].addBond( 6 );
        	residue.atom[ 6 ].addBond( 7, 8, 9 );
        	residue.atom[ 9 ].addBond( 10, 12 );
        	residue.atom[ 10 ].addBond( 11, 14 );
        	residue.atom[ 12 ].addBond( 13, 16 );
        	residue.atom[ 14 ].addBond( 15, 16 );
        	residue.atom[ 16 ].addBond( 17 );
        }
        else if ( residueName.equalsIgnoreCase( "Dhis" ) ){
        	createResidue( "HIS", 18 );
        	residue.atom[ 0 ] = new Atom( "N", 2.044f, 1.352f, 0.000f, -0.463f );
        	residue.atom[ 1 ] = new Atom( "H", 1.721f, 1.837f, 0.824f, 0.252f );
        	residue.atom[ 2 ] = new Atom( "CA", 1.522f, 0.000f, 0.000f, 0.035f );
        	residue.atom[ 3 ] = new Atom( "HA", 1.896f, -0.481f, 0.904f, 0.048f );
        	residue.atom[ 4 ] = new Atom( "C", 0.000f, 0.000f, 0.000f, 0.616f );
        	residue.atom[ 5 ] = new Atom( "O", -0.624f, 1.060f, 0.000f, -0.504f );
        	residue.atom[ 6 ] = new Atom( "CB", 1.988f, -0.769f, -1.232f, -0.098f );
        	residue.atom[ 7 ] = new Atom( "HB1", 1.695f, -0.225f, -2.131f, 0.086f );
        	residue.atom[ 8 ] = new Atom( "HB2", 3.071f, -0.840f, -1.141f, 0.086f );
        	residue.atom[ 9 ] = new Atom( "CG", 1.492f, -2.192f, -1.321f, 0.058f );
        	residue.atom[ 10 ] = new Atom( "ND1", 1.828f, -3.024f, -2.383f, -0.058f );
        	residue.atom[ 11 ] = new Atom( "HD1", 2.411f, -2.776f, -3.170f, 0.306f );
        	residue.atom[ 12 ] = new Atom( "CD2", 0.698f, -2.921f, -0.492f, -0.037f );
        	residue.atom[ 13 ] = new Atom( "HD2", 0.248f, -2.586f, 0.442f, 0.153f );
        	residue.atom[ 14 ] = new Atom( "CE1", 1.251f, -4.195f, -2.182f, 0.114f );
        	residue.atom[ 15 ] = new Atom( "HE1", 1.326f, -5.058f, -2.843f, 0.158f );
        	residue.atom[ 16 ] = new Atom( "NE2", 0.575f, -4.150f, -1.061f, -0.058f );
        	residue.atom[ 17 ] = new Atom( "HE2", 0.041f, -4.917f, -0.677f, 0.306f );
        	addBackboneBonds();
        	residue.atom[ 4 ].addBond( 5 );
        	residue.atom[ 2 ].addBond( 6 );
        	residue.atom[ 6 ].addBond( 7, 8, 9 );
        	residue.atom[ 9 ].addBond( 10, 12 );
        	residue.atom[ 10 ].addBond( 11, 14 );
        	residue.atom[ 12 ].addBond( 13, 16 );
        	residue.atom[ 14 ].addBond( 15, 16 );
        	residue.atom[ 16 ].addBond( 17 );
        }
        else if ( residueName.equalsIgnoreCase( "Dile" ) ){
        	createResidue( "ILE", 19 );
        	residue.atom[ 0 ] = new Atom( "N", 2.044f, 1.352f, 0.000f, -0.463f );
        	residue.atom[ 1 ] = new Atom( "H", 1.721f, 1.837f, 0.824f, 0.252f );
        	residue.atom[ 2 ] = new Atom( "CA", 1.522f, 0.000f, 0.000f, 0.035f );
        	residue.atom[ 3 ] = new Atom( "HA", 1.875f, -0.500f, 0.902f, 0.048f );
        	residue.atom[ 4 ] = new Atom( "C", 0.000f, 0.000f, 0.000f, 0.616f );
        	residue.atom[ 5 ] = new Atom( "O", -0.624f, 1.060f, 0.000f, -0.504f );
        	residue.atom[ 6 ] = new Atom( "CB", 2.009f, -0.733f, -1.245f, -0.012f );
        	residue.atom[ 7 ] = new Atom( "HB", 1.661f, -0.209f, -2.135f, 0.022f );
        	residue.atom[ 8 ] = new Atom( "CG1", 1.459f, -2.156f, -1.246f, -0.049f );
        	residue.atom[ 9 ] = new Atom( "HG11", 1.807f, -2.680f, -0.355f, 0.027f );
        	residue.atom[ 10 ] = new Atom( "HG12", 0.370f, -2.124f, -1.246f, 0.027f );
        	residue.atom[ 11 ] = new Atom( "CG2", 3.533f, -0.777f, -1.245f, -0.085f );
        	residue.atom[ 12 ] = new Atom( "HG21", 3.882f, -1.301f, -2.135f, 0.029f );
        	residue.atom[ 13 ] = new Atom( "HG22", 3.927f, 0.239f, -1.245f, 0.029f );
        	residue.atom[ 14 ] = new Atom( "HG23", 3.882f, -1.301f, -0.355f, 0.029f );
        	residue.atom[ 15 ] = new Atom( "CD1", 1.946f, -2.889f, -2.490f, -0.085f );
        	residue.atom[ 16 ] = new Atom( "HD11", 1.554f, -3.905f, -2.490f, 0.028f );
        	residue.atom[ 17 ] = new Atom( "HD12", 1.598f, -2.365f, -3.380f, 0.028f );
        	residue.atom[ 18 ] = new Atom( "HD13", 3.036f, -2.920f, -2.490f, 0.028f );
        	addBackboneBonds();
        	residue.atom[ 4 ].addBond( 5 );
        	residue.atom[ 2 ].addBond( 6 );
        	residue.atom[ 6 ].addBond( 7, 8, 11 );
        	residue.atom[ 8 ].addBond( 9, 10, 15 );
        	residue.atom[ 11 ].addBond( 12, 13, 14 );
        	residue.atom[ 15 ].addBond( 16, 17, 18 );
        }
        else if ( residueName.equalsIgnoreCase( "Dleu" ) ){
        	createResidue( "LEU", 19 );
        	residue.atom[ 0 ] = new Atom( "N", 2.044f, 1.352f, 0.000f, -0.463f );
        	residue.atom[ 1 ] = new Atom( "H", 1.721f, 1.837f, 0.824f, 0.252f );
        	residue.atom[ 2 ] = new Atom( "CA", 1.522f, 0.000f, 0.000f, 0.035f );
        	residue.atom[ 3 ] = new Atom( "HA", 1.896f, -0.481f, 0.904f, 0.048f );
        	residue.atom[ 4 ] = new Atom( "C", 0.000f, 0.000f, 0.000f, 0.616f );
        	residue.atom[ 5 ] = new Atom( "O", -0.624f, 1.060f, 0.000f, -0.504f );
        	residue.atom[ 6 ] = new Atom( "CB", 1.988f, -0.769f, -1.232f, -0.061f );
        	residue.atom[ 7 ] = new Atom( "HB1", 1.638f, -0.260f, -2.131f, 0.033f );
        	residue.atom[ 8 ] = new Atom( "HB2", 3.077f, -0.816f, -1.241f, 0.033f );
        	residue.atom[ 9 ] = new Atom( "CG", 1.418f, -2.183f, -1.195f, -0.010f );
        	residue.atom[ 10 ] = new Atom( "HG", 1.766f, -2.692f, -0.296f, 0.031f );
        	residue.atom[ 11 ] = new Atom( "CD1", -0.106f, -2.117f, -1.182f, -0.107f );
        	residue.atom[ 12 ] = new Atom( "HD11", -0.513f, -3.128f, -1.155f, 0.034f );
        	residue.atom[ 13 ] = new Atom( "HD12", -0.438f, -1.567f, -0.302f, 0.034f );
        	residue.atom[ 14 ] = new Atom( "HD13", -0.455f, -1.608f, -2.081f, 0.034f );
        	residue.atom[ 15 ] = new Atom( "CD2", 1.884f, -2.951f, -2.427f, -0.107f );
        	residue.atom[ 16 ] = new Atom( "HD21", 1.476f, -3.962f, -2.400f, 0.034f );
        	residue.atom[ 17 ] = new Atom( "HD22", 1.534f, -2.443f, -3.326f, 0.034f );
        	residue.atom[ 18 ] = new Atom( "HD23", 2.973f, -2.999f, -2.436f, 0.034f );
        	addBackboneBonds();
        	residue.atom[ 4 ].addBond( 5 );
        	residue.atom[ 2 ].addBond( 6 );
        	residue.atom[ 6 ].addBond( 7, 8, 9 );
        	residue.atom[ 9 ].addBond( 10, 11, 15 );
        	residue.atom[ 11 ].addBond( 12, 13, 14 );
        	residue.atom[ 15 ].addBond( 16, 17, 18 );
        }
        else if ( residueName.equalsIgnoreCase( "Dlys" ) ){
        	createResidue( "LYS", 22 );
        	residue.atom[ 0 ] = new Atom( "N", 2.044f, 1.352f, 0.000f, -0.463f );
        	residue.atom[ 1 ] = new Atom( "H", 1.721f, 1.837f, 0.824f, 0.252f );
        	residue.atom[ 2 ] = new Atom( "CA", 1.522f, 0.000f, 0.000f, 0.035f );
        	residue.atom[ 3 ] = new Atom( "HA", 1.896f, -0.481f, 0.904f, 0.048f );
        	residue.atom[ 4 ] = new Atom( "C", 0.000f, 0.000f, 0.000f, 0.616f );
        	residue.atom[ 5 ] = new Atom( "O", -0.624f, 1.060f, 0.000f, -0.504f );
        	residue.atom[ 6 ] = new Atom( "CB", 1.988f, -0.769f, -1.232f, -0.098f );
        	residue.atom[ 7 ] = new Atom( "HB1", 1.638f, -0.260f, -2.131f, 0.038f );
        	residue.atom[ 8 ] = new Atom( "HB2", 3.077f, -0.816f, -1.241f, 0.038f );
        	residue.atom[ 9 ] = new Atom( "CG", 1.418f, -2.183f, -1.195f, -0.160f );
        	residue.atom[ 10 ] = new Atom( "HG1", 1.767f, -2.692f, -0.296f, 0.116f );
        	residue.atom[ 11 ] = new Atom( "HG2", 0.328f, -2.136f, -1.186f, 0.116f );
        	residue.atom[ 12 ] = new Atom( "CD", 1.884f, -2.951f, -2.427f, -0.180f );
        	residue.atom[ 13 ] = new Atom( "HD1", 1.534f, -2.443f, -3.326f, 0.122f );
        	residue.atom[ 14 ] = new Atom( "HD2", 2.973f, -2.999f, -2.436f, 0.122f );
        	residue.atom[ 15 ] = new Atom( "CE", 1.313f, -4.366f, -2.389f, -0.038f );
        	residue.atom[ 16 ] = new Atom( "HE1", 1.663f, -4.874f, -1.491f, 0.098f );
        	residue.atom[ 17 ] = new Atom( "HE2", 0.224f, -4.318f, -2.380f, 0.098f );
        	residue.atom[ 18 ] = new Atom( "NZ", 1.762f, -5.106f, -3.577f, -0.138f );
        	residue.atom[ 19 ] = new Atom( "HZ1", 1.385f, -6.042f, -3.552f, 0.294f );
        	residue.atom[ 20 ] = new Atom( "HZ2", 1.440f, -4.636f, -4.410f, 0.294f );
        	residue.atom[ 21 ] = new Atom( "HZ3", 2.771f, -5.149f, -3.585f, 0.294f );
        	addBackboneBonds();
        	residue.atom[ 4 ].addBond( 5 );
        	residue.atom[ 2 ].addBond( 6 );
        	residue.atom[ 6 ].addBond( 7, 8, 9 );
        	residue.atom[ 9 ].addBond( 10, 11, 12 );
        	residue.atom[ 12 ].addBond( 13, 14, 15 );
        	residue.atom[ 15 ].addBond( 16, 17, 18 );
        	residue.atom[ 18 ].addBond( 19, 20, 21 );
        }
        else if ( residueName.equalsIgnoreCase( "Dmet" ) ){
        	createResidue( "MET", 17 );
        	residue.atom[ 0 ] = new Atom( "N", 2.044f, 1.352f, 0.000f, -0.463f );
        	residue.atom[ 1 ] = new Atom( "H", 1.721f, 1.837f, 0.824f, 0.252f );
        	residue.atom[ 2 ] = new Atom( "CA", 1.522f, 0.000f, 0.000f, 0.035f );
        	residue.atom[ 3 ] = new Atom( "HA", 1.896f, -0.481f, 0.904f, 0.048f );
        	residue.atom[ 4 ] = new Atom( "C", 0.000f, 0.000f, 0.000f, 0.616f );
        	residue.atom[ 5 ] = new Atom( "O", -0.624f, 1.060f, 0.000f, -0.504f );
        	residue.atom[ 6 ] = new Atom( "CB", 1.988f, -0.769f, -1.232f, -0.151f );
        	residue.atom[ 7 ] = new Atom( "HB1", 1.638f, -0.260f, -2.131f, 0.027f );
        	residue.atom[ 8 ] = new Atom( "HB2", 3.077f, -0.816f, -1.241f, 0.027f );
        	residue.atom[ 9 ] = new Atom( "CG", 1.418f, -2.183f, -1.195f, -0.054f );
        	residue.atom[ 10 ] = new Atom( "HG1", 1.762f, -2.690f, -0.292f, 0.0652f );
        	residue.atom[ 11 ] = new Atom( "HG2", 0.329f, -2.129f, -1.191f, 0.0652f );
        	residue.atom[ 12 ] = new Atom( "SD", 1.962f, -3.109f, -2.652f, -0.025f );
        	residue.atom[ 13 ] = new Atom( "CE", 1.167f, -4.670f, -2.341f, -0.134f );
        	residue.atom[ 14 ] = new Atom( "HE1", 1.399f, -5.363f, -3.150f, 0.0652f );
        	residue.atom[ 15 ] = new Atom( "HE2", 1.525f, -5.079f, -1.397f, 0.0652f );
        	residue.atom[ 16 ] = new Atom( "HE3", 0.087f, -4.523f, -2.287f, 0.0652f );
        	addBackboneBonds();
        	residue.atom[ 4 ].addBond( 5 );
        	residue.atom[ 2 ].addBond( 6 );
        	residue.atom[ 6 ].addBond( 7, 8, 9 );
        	residue.atom[ 9 ].addBond( 10, 11, 12 );
        	residue.atom[ 12 ].addBond( 13 );
        	residue.atom[ 13 ].addBond( 14, 15, 16 );
        }
        else if ( residueName.equalsIgnoreCase( "Dnme" ) ){
        	createResidue( "NME", 6 );
        	residue.atom[ 0 ] = new Atom( "N", -1.227f, 0.728f, 2.125f, -0.463f );
        	residue.atom[ 1 ] = new Atom( "H", -1.124f, -0.261f, 1.947f, 0.252f );
        	residue.atom[ 2 ] = new Atom( "CA", -1.918f, 1.159f, 3.323f, 0.03502f );
        	residue.atom[ 3 ] = new Atom( "HA1", -1.939f, 2.249f, 3.358f, 0.0586f );
        	residue.atom[ 4 ] = new Atom( "HA2", -2.939f, 0.777f, 3.311f, 0.0586f );
        	residue.atom[ 5 ] = new Atom( "HA3", -1.398f, 0.777f, 4.201f, 0.0586f );
        	addBackboneBonds();
        	residue.atom[ 2 ].addBond( 5 );
        }
        else if ( residueName.equalsIgnoreCase( "Dphe" ) ){
        	createResidue( "PHE", 20 );
        	residue.atom[ 0 ] = new Atom( "N", 2.044f, 1.352f, 0.000f, -0.463f );
        	residue.atom[ 1 ] = new Atom( "H", 1.721f, 1.837f, 0.824f, 0.252f );
        	residue.atom[ 2 ] = new Atom( "CA", 1.522f, 0.000f, 0.000f, 0.035f );
        	residue.atom[ 3 ] = new Atom( "HA", 1.896f, -0.481f, 0.904f, 0.048f );
        	residue.atom[ 4 ] = new Atom( "C", 0.000f, 0.000f, 0.000f, 0.616f );
        	residue.atom[ 5 ] = new Atom( "O", -0.624f, 1.060f, 0.000f, -0.504f );
        	residue.atom[ 6 ] = new Atom( "CB", 1.988f, -0.769f, -1.232f, -0.100f );
        	residue.atom[ 7 ] = new Atom( "HB1", 1.695f, -0.225f, -2.131f, 0.108f );
        	residue.atom[ 8 ] = new Atom( "HB2", 3.071f, -0.840f, -1.141f, 0.108f );
        	residue.atom[ 9 ] = new Atom( "CG", 1.492f, -2.192f, -1.321f, -0.100f );
        	residue.atom[ 10 ] = new Atom( "CD1", 1.856f, -2.993f, -2.410f, -0.150f );
        	residue.atom[ 11 ] = new Atom( "HD1", 2.497f, -2.589f, -3.194f, 0.150f );
        	residue.atom[ 12 ] = new Atom( "CD2", 0.669f, -2.711f, -0.315f, -0.150f );
        	residue.atom[ 13 ] = new Atom( "HD2", 0.385f, -2.087f, 0.533f, 0.150f );
        	residue.atom[ 14 ] = new Atom( "CE1", 1.398f, -4.313f, -2.492f, -0.150f );
        	residue.atom[ 15 ] = new Atom( "HE1", 1.681f, -4.937f, -3.340f, 0.150f );
        	residue.atom[ 16 ] = new Atom( "CE2", 0.210f, -4.031f, -0.397f, -0.150f );
        	residue.atom[ 17 ] = new Atom( "HE2", -0.431f, -4.436f, 0.387f, 0.150f );
        	residue.atom[ 18 ] = new Atom( "CZ", 0.575f, -4.833f, -1.486f, -0.150f );
        	residue.atom[ 19 ] = new Atom( "HZ", 0.217f, -5.861f, -1.550f, 0.150f );
        	addBackboneBonds();
        	residue.atom[ 4 ].addBond( 5 );
        	residue.atom[ 2 ].addBond( 6 );
        	residue.atom[ 6 ].addBond( 7, 8, 9 );
        	residue.atom[ 9 ].addBond( 10, 12 );
        	residue.atom[ 10 ].addBond( 11, 14 );
        	residue.atom[ 12 ].addBond( 13, 16 );
        	residue.atom[ 14 ].addBond( 15, 18 );
        	residue.atom[ 16 ].addBond( 17, 18 );
        	residue.atom[ 18 ].addBond( 19 );
        }
        else if ( residueName.equalsIgnoreCase( "Dpro" ) ){
        	createResidue( "PRO", 14 );
        	residue.atom[ 0 ] = new Atom( "N", 2.067f, 1.413f, 0.000f, -0.229f );
        	residue.atom[ 1 ] = new Atom( "CA", 1.522f, 0.000f, 0.000f, 0.035f );
        	residue.atom[ 2 ] = new Atom( "HA", 1.308f, -0.765f, 0.745f, 0.048f );
        	residue.atom[ 3 ] = new Atom( "C", 0.000f, 0.000f, 0.000f, 0.526f );
        	residue.atom[ 4 ] = new Atom( "O", -0.624f, 1.059f, 0.000f, -0.500f );
        	residue.atom[ 5 ] = new Atom( "CB", 2.632f, -0.865f, -0.521f, -0.115f );
        	residue.atom[ 6 ] = new Atom( "HB1", 2.901f, -1.604f, 0.234f, 0.061f );
        	residue.atom[ 7 ] = new Atom( "HB2", 2.302f, -1.372f, -1.426f, 0.061f );
        	residue.atom[ 8 ] = new Atom( "CG", 3.834f, -0.007f, -0.835f, -0.121f );
        	residue.atom[ 9 ] = new Atom( "HG1", 4.671f, -0.309f, -0.206f, 0.063f );
        	residue.atom[ 10 ] = new Atom( "HG2", 4.135f, -0.100f, -1.878f, 0.063f );
        	residue.atom[ 11 ] = new Atom( "CD", 3.438f, 1.400f, -0.496f, -0.0012f );
        	residue.atom[ 12 ] = new Atom( "HD1", 3.507f, 2.022f, -1.389f, 0.060f );
        	residue.atom[ 13 ] = new Atom( "HD2", 4.105f, 1.791f, 0.271f, 0.060f );
        	residue.atom[ 0 ].addBond( 1, 11 );
        	residue.atom[ 1 ].addBond( 2, 3, 5 );
        	residue.atom[ 3 ].addBond( 4 );
        	residue.atom[ 5 ].addBond( 6, 7, 8 );
        	residue.atom[ 8 ].addBond( 9, 10, 11 );
        	residue.atom[ 11 ].addBond( 12, 13 );
        }
        else if ( residueName.equalsIgnoreCase( "Dser" ) ){
        	createResidue( "SER", 11 );
        	residue.atom[ 0 ] = new Atom( "N", 2.044f, 1.352f, 0.000f, -0.463f );
        	residue.atom[ 1 ] = new Atom( "H", 1.721f, 1.837f, 0.824f, 0.252f );
        	residue.atom[ 2 ] = new Atom( "CA", 1.522f, 0.000f, 0.000f, 0.035f );
        	residue.atom[ 3 ] = new Atom( "HA", 1.896f, -0.481f, 0.904f, 0.048f );
        	residue.atom[ 4 ] = new Atom( "C", 0.000f, 0.000f, 0.000f, 0.616f );
        	residue.atom[ 5 ] = new Atom( "O", -0.624f, 1.060f, 0.000f, -0.504f );
        	residue.atom[ 6 ] = new Atom( "CB", 1.988f, -0.769f, -1.232f, 0.018f );
        	residue.atom[ 7 ] = new Atom( "HB1", 1.638f, -0.260f, -2.131f, 0.119f );
        	residue.atom[ 8 ] = new Atom( "HB2", 3.077f, -0.816f, -1.241f, 0.119f );
        	residue.atom[ 9 ] = new Atom( "OG", 1.453f, -2.095f, -1.196f, -0.550f );
        	residue.atom[ 10 ] = new Atom( "HG", 1.746f, -2.578f, -1.973f, 0.310f );
        	addBackboneBonds();
        	residue.atom[ 4 ].addBond( 5 );
        	residue.atom[ 2 ].addBond( 6 );
        	residue.atom[ 6 ].addBond( 7, 8, 9 );
        	residue.atom[ 9 ].addBond( 10 );
        }
        else if ( residueName.equalsIgnoreCase( "Dthr" ) ){
        	createResidue( "THR", 14 );
        	residue.atom[ 0 ] = new Atom( "N", 2.044f, 1.352f, 0.000f, -0.463f );
        	residue.atom[ 1 ] = new Atom( "H", 1.721f, 1.837f, 0.824f, 0.252f );
        	residue.atom[ 2 ] = new Atom( "CA", 1.522f, 0.000f, 0.000f, 0.035f );
        	residue.atom[ 3 ] = new Atom( "HA", 1.896f, -0.481f, 0.904f, 0.048f );
        	residue.atom[ 4 ] = new Atom( "C", 0.000f, 0.000f, 0.000f, 0.616f );
        	residue.atom[ 5 ] = new Atom( "O", -0.624f, 1.060f, 0.000f, -0.504f );
        	residue.atom[ 6 ] = new Atom( "CB", 1.988f, -0.769f, -1.232f, 0.170f );
        	residue.atom[ 7 ] = new Atom( "HB", 1.580f, -1.779f, -1.205f, 0.082f );
        	residue.atom[ 8 ] = new Atom( "OG1", 3.416f, -0.831f, -1.244f, -0.550f );
        	residue.atom[ 9 ] = new Atom( "HG1", 3.710f, -1.314f, -2.020f, 0.310f );
        	residue.atom[ 10 ] = new Atom( "CG2", 1.499f, -0.057f, -2.490f, -0.191f );
        	residue.atom[ 11 ] = new Atom( "HG21", 1.832f, -0.606f, -3.370f, 0.065f );
        	residue.atom[ 12 ] = new Atom( "HG22", 0.410f, -0.010f, -2.480f, 0.065f );
        	residue.atom[ 13 ] = new Atom( "HG23", 1.906f, 0.953f, -2.516f, 0.065f );
        	addBackboneBonds();
        	residue.atom[ 4 ].addBond( 5 );
        	residue.atom[ 2 ].addBond( 6 );
        	residue.atom[ 6 ].addBond( 7, 8, 10 );
        	residue.atom[ 8 ].addBond( 9 );
        	residue.atom[ 10 ].addBond( 11, 12, 13 );
        }
        else if ( residueName.equalsIgnoreCase( "Dtrp" ) ){
        	createResidue( "TRP", 24 );
        	residue.atom[ 0 ] = new Atom( "N", 2.044f, 1.352f, 0.000f, -0.463f );
        	residue.atom[ 1 ] = new Atom( "H", 1.721f, 1.837f, 0.824f, 0.252f );
        	residue.atom[ 2 ] = new Atom( "CA", 1.522f, 0.000f, 0.000f, 0.035f );
        	residue.atom[ 3 ] = new Atom( "HA", 1.896f, -0.481f, 0.904f, 0.048f );
        	residue.atom[ 4 ] = new Atom( "C", 0.000f, 0.000f, 0.000f, 0.616f );
        	residue.atom[ 5 ] = new Atom( "O", -0.624f, 1.060f, 0.000f, -0.504f );
        	residue.atom[ 6 ] = new Atom( "CB", 1.988f, -0.769f, -1.232f, -0.098f );
        	residue.atom[ 7 ] = new Atom( "HB1", 1.695f, -0.225f, -2.131f, 0.038f );
        	residue.atom[ 8 ] = new Atom( "HB2", 3.071f, -0.840f, -1.141f, 0.038f );
        	residue.atom[ 9 ] = new Atom( "CG", 1.492f, -2.192f, -1.321f, -0.135f );
        	residue.atom[ 10 ] = new Atom( "CD1", 1.752f, -3.078f, -2.293f, 0.044f );
        	residue.atom[ 11 ] = new Atom( "HD1", 2.365f, -2.907f, -3.178f, 0.093f );
        	residue.atom[ 12 ] = new Atom( "CD2", 0.669f, -2.868f, -0.412f, 0.146f );
        	residue.atom[ 13 ] = new Atom( "NE1", 1.072f, -4.288f, -1.950f, -0.352f );
        	residue.atom[ 14 ] = new Atom( "HE1", 1.079f, -5.139f, -2.493f, 0.271f );
        	residue.atom[ 15 ] = new Atom( "CE2", 0.437f, -4.114f, -0.817f, 0.154f );
        	residue.atom[ 16 ] = new Atom( "CE3", 0.103f, -2.412f, 0.785f, -0.173f );
        	residue.atom[ 17 ] = new Atom( "HE3", 0.273f, -1.398f, 1.145f, 0.086f );
        	residue.atom[ 18 ] = new Atom( "CZ2", -0.350f, -5.037f, -0.120f, -0.168f );
        	residue.atom[ 19 ] = new Atom( "HZ2", -0.515f, -6.050f, -0.487f, 0.084f );
        	residue.atom[ 20 ] = new Atom( "CZ3", -0.694f, -3.326f, 1.506f, -0.066f );
        	residue.atom[ 21 ] = new Atom( "HZ3", -1.150f, -3.005f, 2.442f, 0.057f );
        	residue.atom[ 22 ] = new Atom( "CH2", -0.912f, -4.585f, 1.069f, -0.168f );
        	residue.atom[ 23 ] = new Atom( "HH2", -1.535f, -5.257f, 1.658f, 0.084f );
        	addBackboneBonds();
        	residue.atom[ 4 ].addBond( 5 );
        	residue.atom[ 2 ].addBond( 6 );
        	residue.atom[ 6 ].addBond( 7, 8, 9 );
        	residue.atom[ 9 ].addBond( 10, 12 );
        	residue.atom[ 10 ].addBond( 11, 13 );
        	residue.atom[ 13 ].addBond( 14, 15 );
        	residue.atom[ 15 ].addBond( 18 );
        	residue.atom[ 18 ].addBond( 19, 22 );
        	residue.atom[ 22 ].addBond( 23 );
        	residue.atom[ 12 ].addBond( 15, 16 );
        	residue.atom[ 16 ].addBond( 17, 20 );
        	residue.atom[ 20 ].addBond( 21, 22 );
        	residue.atom[ 22 ].addBond( 23 );
        }
        else if ( residueName.equalsIgnoreCase( "Dtyr" ) ){
        	createResidue( "TYR", 21 );
        	residue.atom[ 0 ] = new Atom( "N", 2.044f, 1.352f, 0.000f, -0.463f );
        	residue.atom[ 1 ] = new Atom( "H", 1.721f, 1.837f, 0.824f, 0.252f );
        	residue.atom[ 2 ] = new Atom( "CA", 1.522f, 0.000f, 0.000f, 0.035f );
        	residue.atom[ 3 ] = new Atom( "HA", 1.896f, -0.481f, 0.904f, 0.048f );
        	residue.atom[ 4 ] = new Atom( "C", 0.000f, 0.000f, 0.000f, 0.616f );
        	residue.atom[ 5 ] = new Atom( "O", -0.624f, 1.060f, 0.000f, -0.504f );
        	residue.atom[ 6 ] = new Atom( "CB", 1.988f, -0.769f, -1.232f, -0.098f );
        	residue.atom[ 7 ] = new Atom( "HB1", 1.638f, -0.260f, -2.131f, 0.038f );
        	residue.atom[ 8 ] = new Atom( "HB2", 3.077f, -0.816f, -1.241f, 0.038f );
        	residue.atom[ 9 ] = new Atom( "CG", 1.423f, -2.168f, -1.195f, -0.030f );
        	residue.atom[ 10 ] = new Atom( "CD1", 1.715f, -3.069f, -2.228f, -0.002f );
        	residue.atom[ 11 ] = new Atom( "HD1", 2.348f, -2.758f, -3.058f, 0.064f );
        	residue.atom[ 12 ] = new Atom( "CD2", 0.609f, -2.567f, -0.128f, -0.002f );
        	residue.atom[ 13 ] = new Atom( "HD2", 0.382f, -1.866f, 0.676f, 0.064f );
        	residue.atom[ 14 ] = new Atom( "CE1", 1.191f, -4.367f, -2.193f, -0.264f );
        	residue.atom[ 15 ] = new Atom( "HE1", 1.418f, -5.067f, -2.996f, 0.102f );
        	residue.atom[ 16 ] = new Atom( "CE2", 0.086f, -3.865f, -0.094f, -0.264f );
        	residue.atom[ 17 ] = new Atom( "HE2", -0.548f, -4.176f, 0.737f, 0.102f );
        	residue.atom[ 18 ] = new Atom( "CZ", 0.378f, -4.765f, -1.126f, 0.462f );
        	residue.atom[ 19 ] = new Atom( "OH", -0.131f, -6.027f, -1.093f, -0.528f );
        	residue.atom[ 20 ] = new Atom( "HH", 0.132f, -6.557f, -1.848f, 0.334f );
        	addBackboneBonds();
        	residue.atom[ 4 ].addBond( 5 );
        	residue.atom[ 2 ].addBond( 6 );
        	residue.atom[ 6 ].addBond( 7, 8, 9 );
        	residue.atom[ 9 ].addBond( 10, 12 );
        	residue.atom[ 10 ].addBond( 11, 14 );
        	residue.atom[ 14 ].addBond( 15, 18 );
        	residue.atom[ 18 ].addBond( 19 );
        	residue.atom[ 19 ].addBond( 20 );
        	residue.atom[ 12 ].addBond( 13, 16 );
        	residue.atom[ 16 ].addBond( 17, 18 );
        }
        else if ( residueName.equalsIgnoreCase( "Dval" ) ){
        	createResidue( "VAL", 16 );
        	residue.atom[ 0 ] = new Atom( "N", 2.044f, 1.352f, 0.000f, -0.463f );
        	residue.atom[ 1 ] = new Atom( "H", 1.721f, 1.837f, 0.824f, 0.252f );
        	residue.atom[ 2 ] = new Atom( "CA", 1.522f, 0.000f, 0.000f, 0.035f );
        	residue.atom[ 3 ] = new Atom( "HA", 1.896f, -0.481f, 0.904f, 0.048f );
        	residue.atom[ 4 ] = new Atom( "C", 0.000f, 0.000f, 0.000f, 0.616f );
        	residue.atom[ 5 ] = new Atom( "O", -0.624f, 1.060f, 0.000f, -0.504f );
        	residue.atom[ 6 ] = new Atom( "CB", 1.988f, -0.769f, -1.232f, -0.012f );
        	residue.atom[ 7 ] = new Atom( "HB", 1.638f, -0.260f, -2.131f, 0.024f );
        	residue.atom[ 8 ] = new Atom( "CG1", 3.510f, -0.834f, -1.245f, -0.091f );
        	residue.atom[ 9 ] = new Atom( "HG11", 3.844f, -1.383f, -2.125f, 0.031f );
        	residue.atom[ 10 ] = new Atom( "HG12", 3.918f, 0.176f, -1.271f, 0.031f );
        	residue.atom[ 11 ] = new Atom( "HG13", 3.860f, -1.343f, -0.346f, 0.031f );
        	residue.atom[ 12 ] = new Atom( "CG2", 1.418f, -2.183f, -1.195f, -0.091f );
        	residue.atom[ 13 ] = new Atom( "HG21", 1.751f, -2.732f, -2.075f, 0.031f );
        	residue.atom[ 14 ] = new Atom( "HG22", 1.767f, -2.692f, -0.296f, 0.031f );
        	residue.atom[ 15 ] = new Atom( "HG23", 0.328f, -2.136f, -1.186f, 0.031f );
        	addBackboneBonds();
        	residue.atom[ 4 ].addBond( 5 );
        	residue.atom[ 2 ].addBond( 6 );
        	residue.atom[ 6 ].addBond( 7, 8, 12 );
        	residue.atom[ 8 ].addBond( 9, 10, 11 );
        	residue.atom[ 12 ].addBond( 13, 14, 15 );
        }*/

        // Make the bonds bidirectional
        completeBonds();
        for(int i=0; i<residue.numberOfAtoms; i++) {
            residue.atom[i].residueAtomNumber = i;
        }
        return(residue);
    }

    public void createResidue( String residueName, int numberOfResidueAtoms ) {
        residue = new Residue();
        residue.name = residueName;
        residue.numberOfAtoms = numberOfResidueAtoms;
        residue.atom = new Atom[ numberOfResidueAtoms ];
    }

    public void addBackboneBonds() {
        residue.atom[ 0 ].addBond( 1, 2 );
        residue.atom[ 2 ].addBond( 3, 4 );
    }

    // Complete bonds created by amino acid template. Bonds
    //  that are only unidirectional are made bidirectional
    // Important for using Amber96PolyPeptideResiude
    public void completeBonds() {
        for(int i=0; i<residue.numberOfAtoms; i++) {
            for(int j=0; j<residue.atom[i].numberOfBonds; j++) {
                if (!(residue.atom[residue.atom[i].bond[j]].bondedTo(i)))
                    residue.atom[residue.atom[i].bond[j]].addBond(i);
            }
        }
    }



}
