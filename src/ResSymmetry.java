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
//	ResSymmetry.java
//
//	Version:           2.1 beta
//
//
//	  authors:
// 	  initials    name                 organization                email
//	 ---------   -----------------    ------------------------    ----------------------------
//	  MAH           Mark A. Hallen	  Duke University               mah43@duke.edu
///////////////////////////////////////////////////////////////////////////////////////////////

//MH 1/3/11: This file was adapted from my specbas energy-function package
//It is used here for the identifyRots functionality
//because rotamers for symmetric residues can have two different sets of dihedrals
//depending on how the equivalent atoms are named (e.g. PHE CD1 and CD2)


//This class is for dealing with symmetries in residues
//Particularly those causing ambiguities in atom naming
public class ResSymmetry {


    public static int getSymmetry(String resName){
        //For the given residue type, specify the number of symmetric forms that could cause
        //heavy-atom name ambiguity

        if ( resName.equalsIgnoreCase("ASP") ||
                resName.equalsIgnoreCase("GLU") ||
                resName.equalsIgnoreCase("PHE") ||
                resName.equalsIgnoreCase("ARG") ||
                resName.equalsIgnoreCase("VAL") ||
                resName.equalsIgnoreCase("TYR") )
            
            return 2;
        
        else
            return 1;
    }



    public static String getPermutedAtomName( String atomName, String resName, int s ){
        //In residue resName,
        //give the name of the atom to which atomName is permuted in symmetry state # s

        if(s==0)//0'th symmetry state = identity permutation
            return atomName;


        if( s == 1 ){

            String stem = getStem(atomName);

            if( resName.equalsIgnoreCase("ASP") ){
                //Switch OD1, OD2

                if( atomName.equalsIgnoreCase("OD1") )
                    return "OD2";
                else if(atomName.equalsIgnoreCase("OD2"));
                    return "OD1";
            }


            if( resName.equalsIgnoreCase("GLU") ){
                //Switch OE1, OE2

                if( atomName.equalsIgnoreCase("OE1") )
                    return "OE2";
                else if(atomName.equalsIgnoreCase("OE2"));
                    return "OE1";
            }




            if( resName.equalsIgnoreCase("PHE") || resName.equalsIgnoreCase("TYR") ){
                //Ring-flipping symmetry
                //Involves gamma and delta hydrogens and carbons

                if( stem.endsWith("D") || stem.endsWith("E") ){

                    if( getFirstDigit( atomName ) == '1' )
                        return atomName.replaceFirst("1", "2");
                    else if( getFirstDigit( atomName ) == '2' )
                        return atomName.replaceFirst("2", "1");
                    //This should be robust to 1HG versus HG1 format differences
                }
            }




            if( resName.equalsIgnoreCase("ARG") ){
                //Switch eta nitrogens and their hydrogens

                if( stem.equalsIgnoreCase("NH") || stem.equalsIgnoreCase("HH") ){

                    if( getFirstDigit( atomName ) == '1' )
                        return atomName.replaceFirst("1", "2");
                    else if ( getFirstDigit( atomName ) == '2' )
                        return atomName.replaceFirst("2", "1");
                }
            }



            if( resName.equalsIgnoreCase("VAL") ){
                //Switch gamma carbons and their hydrogens

                if( stem.endsWith("G") ){

                    if( getFirstDigit( atomName ) == '1' )
                        return atomName.replaceFirst("1", "2");
                    else if( getFirstDigit( atomName ) == '2' )
                        return atomName.replaceFirst("2", "1");
                }
            }


            return atomName;//Other atoms' names stay the same

        }

        //No symmetry numbers over 1 are recognized as of now
        System.err.println("ERROR: Symmetry number " + s + " not recognized for residue " + resName );
        System.exit(1);
        return "ERROR";
    }





    public static String getStem(String s){
        //Get the "stem" of an atom name (the part consisting of letters rather than digits)
        String stem = "";
        for(int a=0; a<s.length(); a++){
            char c = s.charAt(a);
            if( Character.isLetter(c) )
                stem = stem + c;
        }
        return stem;
    }


    public static char getFirstDigit(String s){

        for(int a=0; a<s.length(); a++){
            char c = s.charAt(a);
            if( Character.isDigit(c ) )
                return c;
        }

        //No digit found
        System.err.println("ERROR: No digit found in string " + s );
        System.exit(1);
        return '!';
    }



    //Omitting the permutation math functions that are in specbas.ResSymmetry: not needed in OSPREY


}
