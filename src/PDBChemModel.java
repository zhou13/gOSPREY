/*
	This file is part of OSPREY.

	OSPREY Protein Redesign Software Version 2.1 beta
	Copyright (C) 2001-2009 Bruce Donald Lab, Duke University

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

	<signature of Bruce Donald>, 12 Apr, 2009
	Bruce Donald, Professor of Computer Science
*/

////////////////////////////////////////////////////////////////////////////////////////////////
// PDBChemModel.java
//
//  Version:           2.1 beta
//
//
// authors:
//    initials    name                 organization                email
//   ---------   -----------------    ------------------------    ----------------------------
//     RHL        Ryan Lilien          Dartmouth College           ryan.lilien@dartmouth.edu
//	   ISG		  Ivelin Georgiev	   Duke University			   ivelin.georgiev@duke.edu
//     MAH        Mark A. Hallen	Duke University         mah43@duke.edu
//
///////////////////////////////////////////////////////////////////////////////////////////////

/*
 * Written by Ryan Lilien (2000-2004); some modifications by Ivelin Georgiev (2004-2009)
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

import java.io.InputStream;
import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

/**
 * Reads a pdb file and creates a molecule object.
 */
class PDBChemModel {

    ArrayList<Integer> sortedStarts, sortedEnds;//for use by sortSS: keeps track of secondary structures to sort
    ArrayList<Character> sortedStrands;

    float overlapThresh = 0.4f;//Steric overlap threshold for auto-fixing residues

    // Reads a pdb file and creates a molecule object
    PDBChemModel (Molecule m, InputStream is) throws Exception {

        int modelAtomNumber = 0;
        int residueNumber = 0;
        int strandNumber = 0;
        String atomName = "ZZZ", residueName = "ZZZ", strandName = "0";
        String lastResidueName = "ZZZ", fullResidueName = "ZZZZZZZZZ";
        String elementType = "  ";
        String curLine = null;
        String tmpStg = null;
        int tmpInt;
        char[] tmpChr = new char[15];
        float	x = 0f, y = 0f, z = 0f, B = 0f;
        Atom	newAtom;
        boolean newStrandPending = true;  // Will add the first strand with the first atom
        Residue newResidue = null;

        BufferedReader bufread = new BufferedReader(new InputStreamReader(is));

        ArrayList<Integer> helixStarts = new ArrayList<Integer>();//Residues where helices start
        ArrayList<Integer> helixEnds = new ArrayList<Integer>();//Residues where they end
        ArrayList<Character> helixStrands = new ArrayList<Character>();
        ArrayList<Integer> sheetStarts = new ArrayList<Integer>();
        ArrayList<Integer> sheetEnds = new ArrayList<Integer>();
        ArrayList<Character> sheetStrands = new ArrayList<Character>();


        curLine = bufread.readLine();
        if (curLine == null) { // stop if we've reached EOF
            bufread.close();
        }

        boolean emptyStrand = true;//Indicates no atoms have been added to the current strand yet

        while(curLine != null) {
            // First pad line to 80 characters
            tmpInt = curLine.length();
            for (int i=0; i < (80-tmpInt); i++)
                curLine += " ";

            if ((curLine.regionMatches(true,0,"ATOM  ",0,6)) || (curLine.regionMatches(true,0,"HETATM",0,6))) {

                if( EnvironmentVars.autoFix ) { //Ignore alternates other than A; treat A alternates as the real structure
                    char alt = curLine.charAt(16);//This specifies which alternate the atom is (space if not an alternate)
                    if( ( alt != ' ' ) && ( alt != 'A' ) ) {
                        curLine = bufread.readLine();
                        continue;//skip the line and  go to the next one
                    }
                }


                emptyStrand = false;

                // Is an ATOM line
                tmpStg = curLine.substring(6,11);  // Snag atom serial number
                tmpStg = tmpStg.trim();
                modelAtomNumber = (new Integer(tmpStg)).intValue();
                atomName = curLine.substring(12,16);  // Snag atom name
                atomName = atomName.trim();
                residueName = curLine.substring(17,20);  // Snag short residue name
                residueName = residueName.trim();
                fullResidueName = curLine.substring(17,26);  // Snag full residue atom name
                fullResidueName = fullResidueName.trim();

                if (!(fullResidueName.equals(lastResidueName))) {
                    if (newResidue != null) {
                        standardizeHNamesBB(newResidue);
                        if (newStrandPending)
                            m.addResidue(strandNumber-1, newResidue);
                        else
                            m.addResidue(strandNumber, newResidue);
                    }
                    newResidue = new Residue();
                    newResidue.name = residueName;
                    newResidue.fullName = fullResidueName;
                    lastResidueName = fullResidueName;
                    residueNumber++;
                }

                if (newStrandPending) {
                    strandName = curLine.substring(21,22);  // Snag strand name
                    m.addStrand(strandName);
                    newStrandPending = false;
                }
                tmpStg = curLine.substring(30,38);  // Snag x coord
                x = (float) new Double(tmpStg).doubleValue();
                tmpStg = curLine.substring(38,46);  // Snag y coord
                y = (float) new Double(tmpStg).doubleValue();
                tmpStg = curLine.substring(46,54);  // Snag z coord
                z = (float) new Double(tmpStg).doubleValue();

                tmpStg = curLine.substring(60,66).trim();  // Snag B-factor
                if(!tmpStg.isEmpty())
                    B = (float) new Double(tmpStg).doubleValue();

                elementType = curLine.substring(76,78);  // Snag atom elementType
                elementType = elementType.trim();
                // If we can't get element type from substring(76,78) snag
                //  the first character of the atom name
                if (elementType.equalsIgnoreCase(""))
                    elementType = getEleType(curLine.substring(12,15));
                newAtom = new Atom(atomName,x,y,z);
                newAtom.modelAtomNumber = modelAtomNumber;
                newAtom.strandNumber = strandNumber;
                newAtom.elementType = elementType;
                newAtom.BFactor = B;
                newResidue.addAtom(newAtom);
            } // end ATOM line
            else if (curLine.regionMatches(true,0,"TER   ",0,6)) {
                // Is the end of a strand
                if(!emptyStrand) {
                    lastResidueName = "ZZZ";
                    residueNumber = 0;
                    strandNumber++;
                    newStrandPending = true;
                }
                emptyStrand = true;
            } // end TER line
            else if(curLine.regionMatches(true,0,"HELIX  ",0,7)) { //Read helix records
                helixStarts.add(new Integer(curLine.substring(21,25).trim()));
                helixEnds.add(new Integer(curLine.substring(33,37).trim()));
                helixStrands.add(curLine.charAt(19));
            } else if(curLine.regionMatches(true,0,"SHEET  ",0,7)) {
                sheetStarts.add(new Integer(curLine.substring(22,26).trim()));
                sheetEnds.add(new Integer(curLine.substring(33,37).trim()));
                sheetStrands.add(curLine.charAt(21));
            } else // is a line we skip
                ;
            curLine = bufread.readLine();  // attempt to read next line
        }  // end while (curLine != null)

        standardizeHNamesBB(newResidue);

        if (newStrandPending)
            m.addResidue(strandNumber-1,newResidue);
        else
            m.addResidue(strandNumber,newResidue);
        bufread.close();  // close the buffer

        assignSecStruct(m, helixStarts, helixEnds, helixStrands, sheetStarts, sheetEnds, sheetStrands);

        //Determine the bonds between the atoms in the molecule
        m.determineBonds();

        // Assign the molecule relative atom numbers
        m.updateMoleculeAtomNumbers();


        for( Residue res : m.residue ) { //Identify the pucker associated with all the proline rings
            if( res.name.equalsIgnoreCase("PRO") )
                assignPucker(m, res);
        }


        if(EnvironmentVars.autoFix)
            fixResidues(m);

    }



    // This function pulls the element type from
    //  the atom name
    private String getEleType(String str) {

        int start=0, end=-1;
        int i=0;
        while( (str.charAt(i)==' ') || ((str.charAt(i)>='0') && (str.charAt(i)<='9')) ) {
            i++;
        }
        start = i;
        end = i++;
        if (i<str.length())
            if((str.charAt(i)>='a') && (str.charAt(i)<='z'))
                end = i;
        return(str.substring(start,end+1));
    }


    public static void standardizeHNamesBB(Residue newResidue) {
        //Change some nonstandard backbone hydrogen names so the rest of the program can understand them

        //The amide hydrogen is sometimes known as HN or HN1
        Atom HN = newResidue.getAtomByName("HN");
        if( HN != null ) {
            HN.name = "H";
            HN.isBBatom = true;//This would not have been set with the old name
        } else {
            Atom HN1 = newResidue.getAtomByName("HN1");
            if( HN1 != null ) {
                HN1.name = "H";
                HN1.isBBatom = true;//This would not have been set with the old name
            }
        }

        Atom HA1 = newResidue.getAtomByName("HA1");
        if( HA1 != null ) { //Convert HA1/HA2 numbering system for glycine to HA2/HA3
            newResidue.getAtomByName("HA2").name = "HA3";
            HA1.name = "HA2";
        }
    }

    public void assignSecStruct(Molecule m, ArrayList<Integer> helixStarts, ArrayList<Integer> helixEnds,
                                ArrayList<Character> helixStrands, ArrayList<Integer> sheetStarts, ArrayList<Integer>sheetEnds,
                                ArrayList<Character> sheetStrands) throws Exception {


        char[] strandList = new char[m.strand.length];//Chain IDs in order
        for(int str=0; str<m.strand.length; str++)
            strandList[str] = m.strand[str].name.charAt(0);

        //First need to sort lists of secondary structure elements, using StrandList for strand order
        sortSS(helixStarts, helixEnds, helixStrands, strandList);
        helixStarts = sortedStarts;
        helixEnds = sortedEnds;
        helixStrands = sortedStrands;

        sortSS(sheetStarts, sheetEnds, sheetStrands, strandList);
        sheetStarts = sortedStarts;
        sheetEnds = sortedEnds;
        sheetStrands = sortedStrands;


        int curHelix = 0;//position in the helix ArrayLists
        int curSheet = 0;//Same for sheets
        int curSecondaryStruct = Residue.LOOP;//Current secondary structure.  Loop unless otherwise specified.

        for( int resNum=0; resNum<m.residue.length; resNum++ ) {

            Residue res = m.residue[resNum];

            while (true) {//Get curSheet, curHelix, curSecondaryStruct up to date

                if( curSecondaryStruct == Residue.LOOP ) { //Check for starting of a helix or sheet
                    if( curHelix < helixStarts.size() ) { //If there are more helices to look for
                        if( (res.getResNumber() >= helixStarts.get(curHelix) )
                                && ( strandList[res.strandNumber] == helixStrands.get(curHelix)) ) { //If we reached the start of the next one
                            curSecondaryStruct = Residue.HELIX;
                            continue;
                        }
                    }
                    if( curSheet < sheetStarts.size() ) {
                        if( (res.getResNumber() >= sheetStarts.get(curSheet) )
                                && ( strandList[res.strandNumber] == sheetStrands.get(curSheet)) ) {
                            curSecondaryStruct = Residue.SHEET;
                            continue;
                        }
                    }
                }

                else if(curSecondaryStruct == Residue.HELIX ) { //Check for the end of a helix
                    if( ( res.getResNumber() > helixEnds.get(curHelix) ) ||
                            ( strandList[res.strandNumber] != helixStrands.get(curHelix)) ) {
                        curSecondaryStruct = Residue.LOOP;
                        curHelix++;
                        continue;//
                    }
                }

                else if(curSecondaryStruct == Residue.SHEET ) { //Check for the end of a sheet
                    if( ( res.getResNumber() > sheetEnds.get(curSheet) ) ||
                            ( strandList[res.strandNumber] != sheetStrands.get(curSheet)) ) {
                        curSecondaryStruct = Residue.LOOP;
                        curSheet++;
                        continue;
                    }
                }

                break;

            }

            res.secondaryStruct = curSecondaryStruct;
        }
    }


    private void sortSS(ArrayList<Integer> starts, ArrayList<Integer> ends, ArrayList<Character> strands,
                        char[] strandList) {
        //Used to sort lists of starts, ends, and strand IDs of secondary structure elements
        //Called by the constructor above
        //Sorting is first by strand, order given by strandList
        sortedStarts = new ArrayList<Integer>(starts.size());
        sortedEnds = new ArrayList<Integer>(ends.size());
        sortedStrands = new ArrayList<Character>(strands.size());

        int strandOffset = 0;

        for(int a=0; a<strandList.length; a++) { //There will tend to be only a few strands so this loop is inexpensive
            char curStrand = strandList[a];
            int ssCount = 0;

            for(int b=0; b<strands.size(); b++) {
                if(strands.get(b) == curStrand) {
                    ssCount++;
                    sortedStarts.add( starts.get(b) );
                    sortedEnds.add( ends.get(b) );
                    sortedStrands.add(curStrand);
                }
            }

            Collections.sort(sortedStarts.subList(strandOffset, strandOffset + ssCount));
            Collections.sort(sortedEnds.subList(strandOffset, strandOffset + ssCount));

            strandOffset += ssCount;

        }

    }



    static public void assignPucker(Molecule m, Residue res) {
        //Given a molecule and a proline residue in it, assign the ring pucker

        int CANum = res.getAtomNameToMolnum("CA");
        int CBNum = res.getAtomNameToMolnum("CB");
        int CGNum = res.getAtomNameToMolnum("CG");
        int CDNum = res.getAtomNameToMolnum("CD");

        double chi2 = m.getTorsion(CANum, CBNum, CGNum, CDNum);

        if( chi2 > 0 )
            res.pucker = ProlineFlip.UP;
        else
            res.pucker = ProlineFlip.DOWN;

    }





    public void fixResidues(Molecule m) {
        //This function deals with residues provided in the PDB file
        //in a form unsuitable for use by OSPREY with the current amino-acid and generic-residue
        //templates
        //This function is triggered by EnvironmentVars.autoFix == true
        //Residues with unrecognized names are deleted
        //His or Cys residues are renamed to Hid, Hie, Hip, or Cyx as needed
        //Incomplete residues are mutated if enough backbone information is available
        //and deleted otherwise


        AminoAcidTemplates aat = null;
        GenericResidueTemplates grt = null;

        //We use the amino-acid and generic-residue templates
        try {
            aat = new AminoAcidTemplates();
            grt = new GenericResidueTemplates();
        } catch (FileNotFoundException e) {
            System.out.println("ERROR: Template File Not Found: "+e);
            System.exit(0);
        } catch ( Exception e ) {
            System.out.println("ERROR: An error occurred while reading a template file: "+e);
            System.exit(0);
        }


        ArrayList<Integer> resToDelete = new ArrayList<Integer>(); //A list of molecule residue numbers
        //for which the residues should be deleted
        //resToDelete needs to be in order

        for( int molResNum=0; molResNum<m.numberOfResidues; molResNum++ ) {

            Residue res = m.residue[molResNum];


            //First deal with the residue name
            if( res.name.equalsIgnoreCase("HIS") ) {

                String protName = null;

                int HDCount = 0;//Count delta and epsilon hydrogens
                int HECount = 0;

                for( Atom at : res.atom ) {
                    if( at.name.contains("HD") )
                        HDCount++;
                    else if(at.name.contains("HE"))
                        HECount++;
                }

                //Determine protonation state
                if( ( HDCount == 1 ) && ( HECount == 2 ) )
                    protName = "HIE";
                else if( ( HDCount == 2 ) && ( HECount == 1 ) )
                    protName = "HID";
                else if( ( HDCount == 2 ) && ( HECount == 2 ) )
                    protName = "HIP";
                else {
                    System.err.println("ERROR: Invalid protonation state for " + res.fullName );
                    System.exit(1);
                }

                //Rename residue accordingly
                res.fullName = res.fullName.replaceFirst("HIS", protName);
                res.name = protName;
            } else if( res.name.equalsIgnoreCase("CYS") ) {

                //See if the residue has a gamma hydrogen, if not then rename it CYX
                if( res.getAtomByName("HG") == null ) {
                    res.fullName = res.fullName.replaceFirst("CYS", "CYX");
                    res.name = "CYX";
                }

            }



            //Now look for a template with the same name as the current residue
            //It can be either an AA or generic-residue template

            //Use N- or C-terminal templates if the H1 and OXT atoms are present respectively
            Residue[] AATemplRes = aat.aaResidues;
            Residue[] GRTemplRes = grt.grResidues;

            if( res.getAtomByName("H1") != null )
                AATemplRes = aat.aaNTResidues;
            else if( res.getAtomByName("OXT") != null )
                AATemplRes = aat.aaCTResidues;


            Residue templ = null;//Template residue
            boolean isTemplAA = false;

            for( Residue cand : AATemplRes ) {
                if( cand != null ) { //Some of the templates may be null
                    if( cand.name.equalsIgnoreCase(res.name) ) {
                        templ=cand;
                        isTemplAA = true;
                        break;
                    }
                }
            }

            if( ! isTemplAA ) {

                for( Residue cand : GRTemplRes ) {
                    if( cand != null ) { //Some of the templates may be null
                        if( cand.name.equalsIgnoreCase(res.name) ) {
                            templ=cand;
                            break;
                        }
                    }
                }
            }

            //Flag residue for deletion if unrecognized
            if( templ == null ) {
                resToDelete.add( molResNum );
                System.out.println("Warning: Deleting unrecognized residue " + res.fullName );
            }
            //If residue has less atoms than its template, restore residue from template if possible; otherwise delete
            //This is needed when the end of a Lys or Arg sidechain is disordered, etc.
            else if( res.numberOfAtoms < templ.numberOfAtoms ) {
                //Note we leave residues with more atoms than their templates, or the wrong atoms, alone
                //Some cases like this (e.g. ligand residues that are both N- and C-terminal) are OK
                //And we don't really have a good way of dealing with them...they will lead to errors on the assignment of AMBER templates later


                if( ! isTemplAA ) { //Can't mutate a non-AA residue
                    resToDelete.add( molResNum );
                    System.out.println("Warning: Deleting incomplete residue " + res.fullName );
                } else {

                    //See if all atoms needed for mutating are needed
                    boolean mutPossible = true;

                    if( ( res.getAtomByName("CB") == null ) &&
                            ( ! res.name.equalsIgnoreCase("GLY") ) )
                        mutPossible = false;

                    if( ( res.getAtomByName("CD") == null ) &&
                            ( res.name.equalsIgnoreCase("PRO") ) )
                        mutPossible = false;

                    //We require an intact backbone for mutating
                    for( Atom at : templ.atom ) {
                        if( at.isBBatom ) {
                            if( res.getAtomByName( at.name ) == null ) {
                                mutPossible = false;
                                break;
                            }
                        }
                    }


                    if( mutPossible ) {
                        //Fix the residue by changing it to the template if possible
                        //We will fix it using the AA rotamer library (from EnvironmentVars)

                        StrandRotamers sr = new StrandRotamers( EnvironmentVars.aaRotLib, m.strand[res.strandNumber] );
                        ProbeStericCheck psc = new ProbeStericCheck();
                        RotMatrix r = new RotMatrix();

                        HashMap<String, float[]> oldCoord = new HashMap<String, float[]>();//Record known coordinates
                        for( Atom at : res.atom )
                            oldCoord.put(at.name, at.coord);

                        sr.changeResidueType(m, res.strandResidueNumber, res.name, true);


                        int numRot = EnvironmentVars.aaRotLib.getNumRotamers( res.name );

                        float minMSDGood = Float.POSITIVE_INFINITY;//Lowest MSD of a sterically good rotamer to the template
                        float minMSDBad = Float.POSITIVE_INFINITY;//Lowest MSD of a sterically bad rotamer to the template
                        //(These are based on atoms of the same name in the residue and the template)

                        int argminMSDGood = -1;//The corresponding rotamer indices
                        int argminMSDBad = -1;

                        for(int rot=0; rot<numRot; rot++) {

                            sr.applyRotamer(m, res.strandResidueNumber, rot);

                            int firstMolAtNum = res.atom[0].moleculeAtomNumber;
                            int lastMolAtNum = res.atom[res.atom.length-1].moleculeAtomNumber;

                            float msd = 0;
                            boolean stericallyGood = true;

                            for( Atom at : res.atom ) {

                                int molAtNum1 = at.moleculeAtomNumber;

                                if( oldCoord.containsKey(at.name) ) {
                                    float dev = r.norm( r.subtract( m.getActualCoord(molAtNum1),
                                                                    oldCoord.get(at.name) ) );
                                    msd += dev*dev;
                                }


                                for( Atom at2 : m.atom ) {

                                    if( ( at2.moleculeAtomNumber < firstMolAtNum ) ||
                                            ( at2.moleculeAtomNumber > lastMolAtNum ) ) { //If at2 isn't in res
                                        //(rotamers shouldn't have interior steric clashes, so we check for clashes with other residues)
                                        //Note we do check for clashes against residues that will be deleted later (e.g. waters w/ no hydrogens)
                                        //since even if we aren't going to model them they can still inform the placement of the residues that are modeled
                                        if( ! psc.isAllowedSteric(m, at, at2, overlapThresh ) ) {
                                            stericallyGood = false;
                                        }
                                    }
                                }
                            }


                            if( stericallyGood && ( msd < minMSDGood ) ) {
                                minMSDGood = msd;
                                argminMSDGood = rot;
                            } else if( ( ! stericallyGood ) && ( msd < minMSDBad ) ) {
                                minMSDBad = msd;
                                argminMSDBad = rot;
                            }
                        }


                        if( argminMSDGood > -1 )//If there are any sterically good rotamers
                            //Apply the one closest (in MSD) to the original geometry
                            sr.applyRotamer(m, res.strandResidueNumber, argminMSDGood);
                        else if( argminMSDBad > -1 ) { //All rotamers are sterically bad so apply the one closest (in MSD) to the original geometry
                            //If argminMSDGood == argminMSDBad == -1 there are no rotamers so we don't need to assign one
                            sr.applyRotamer(m, res.strandResidueNumber, argminMSDBad);
                            System.out.println( "Warning: " + res.fullName + " incomplete but no sterically good rotamers found; applying clashing rotamer" );
                        }
                    } else { //Otherwise get rid of the residue
                        resToDelete.add( molResNum );
                        System.out.println("Warning: Deleting incomplete residue " + res.fullName );
                    }
                }
            }
        }

        //Now delete the residues that need to be deleted
        //Since the residues are in sorted order, doing this backwards makes sure the numbers are right
        for( int a=resToDelete.size()-1; a>=0; a-- ) {
            m.deleteResidue( resToDelete.get(a) );
        }



    }//End of fixResidues

}
