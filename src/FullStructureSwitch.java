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
//	FullStructureSwitch.java
//
//	Version:           2.1 beta
//
//
//	  authors:
// 	  initials    name                 organization                email
//	 ---------   -----------------    ------------------------    ----------------------------
//	  MAH           Mark A. Hallen	  Duke University               mah43@duke.edu
///////////////////////////////////////////////////////////////////////////////////////////////

import java.io.*;
import java.util.*;


//This perturbation changes the entire structure, including the template,
//to one specified by a different PDB file
//Because this perturbation changes the whole structure, it should be applied before any others
//(it will overwrite all backbone changes made before its application)
public class FullStructureSwitch extends Perturbation {


    float[][][] templateCoord;
    //templateCoord[a][b] holds the template sections of actualCoordinates in each alternate sructure
    //Each of the structures must have the same chains, residues, and atoms in the same order, except for active-site mutations
    //(having differences would make the energy comparison invalid, but removing differing parts
    //from the experimental structure before computation may be fine)

    int numProteinFlexRes;//Number of flexible protein residues (the last affected residue
    //may be a non-protein ligand but the others will count)

    int numStructs;//Number of structures

    String PDBList;//List of alternate structures

    ArrayList< ArrayList < HashMap < String, float[]>>> flexBB = new ArrayList<ArrayList<HashMap<String,float[]>>>();//Alternate flexible-residue backbones



    public FullStructureSwitch(Molecule m, int resList[], BufferedReader br) {

        resAffected = resList;//All flexible residues
        type = "FULL STRUCTURE SWITCH";
        this.m = m;


        if( m.strand[ m.residue[ resAffected [resAffected.length - 1] ].strandNumber ].isProtein )
            numProteinFlexRes = resAffected.length;
        else
            numProteinFlexRes = resAffected.length - 1;


        //Count sections of template
        int numSections = 0;
        if( resAffected[0] != 0 )
            numSections++;

        for( int a=1; a<numProteinFlexRes; a++) {
            if( resAffected[a] > resAffected[a-1] + 1 )
                numSections++;
        }

        if( resAffected[numProteinFlexRes-1] != m.residue.length - 1 )
            numSections++;


        try {
            PDBList = br.readLine();//This line will contain the list of alternate PDB files
            StringTokenizer st = new StringTokenizer(PDBList," ");
            numStructs = st.countTokens() + 1;

            //Store the different structures (everything in the template and non-protein residues; backbones for flexible residues)
            templateCoord = new float[numStructs][numSections][];
            Molecule m2;

            for( int structIndex=0; structIndex<numStructs; structIndex++) {

                ArrayList < HashMap < String, float[]>> structBB = new ArrayList< HashMap < String, float[]>>();

                if(structIndex == 0)
                    m2 = m;
                else {
                    FileInputStream is = new FileInputStream( st.nextToken() );
                    m2 = new Molecule();
                    new PDBChemModel(m2,is);
                }

                int section=0;

                if( resAffected[0] != 0 ) {
                    int secLength = 3*m2.residue[resAffected[0]].atom[0].moleculeAtomNumber;
                    templateCoord[structIndex][0] = new float[secLength];
                    System.arraycopy( m2.actualCoordinates, 0, templateCoord[structIndex][0], 0, secLength );
                    section++;
                }

                structBB.add( storeResBB( m2.residue[resAffected[0]], m2 ) );
                int curIndex = -1;
                if( resAffected[0] < m2.residue.length - 1)
                    curIndex = 3*m2.residue[resAffected[0]+1].atom[0].moleculeAtomNumber;//current index in m.actualCoordinates


                for( int a=1; a<numProteinFlexRes; a++) {

                    int molResNum = resAffected[a];

                    if( molResNum > resAffected[a-1] + 1 ) { //There is a template section before this residue
                        int secLength = 3*m2.residue[resAffected[a]].atom[0].moleculeAtomNumber - curIndex;
                        templateCoord[structIndex][section] = new float[secLength];
                        System.arraycopy( m2.actualCoordinates, curIndex, templateCoord[structIndex][section], 0, secLength );
                        section++;
                    }

                    //Now store the residue backbone itself
                    structBB.add( storeResBB( m2.residue[molResNum] , m2 ) );

                    if(molResNum < m2.residue.length - 1)
                        curIndex = 3*m2.residue[molResNum + 1].atom[0].moleculeAtomNumber;
                    else
                        curIndex = -1;//We're at the end of the molecule
                }

                //Change the final template section
                if( curIndex != -1) {
                    int secLength = m2.actualCoordinates.length - curIndex;
                    templateCoord[structIndex][section] = new float[secLength];
                    System.arraycopy( m2.actualCoordinates, curIndex, templateCoord[structIndex][section], 0, secLength );
                }

                flexBB.add(structBB);
            }
        } catch(Exception e) {
            System.err.println();
            e.printStackTrace();
        }

    }



    //Override applyPerturbation and undo in order to handle the template changes

    @Override
    public boolean applyPerturbation(float param) {
        //The parameter need to be an integer less than templateCoord.length

        if( param == 0 ) {
            curParam = 0;
            return true;
        }

        boolean outcome = doPerturbationMotion(param);

        if(outcome)
            curParam = param;
        else
            curParam=0;

        return outcome;

    }


    //This contains the sidechain idealization and chi1 adjustments, unlike the default version
    public boolean doPerturbationMotion(float param) {

        int structIndex = (int)param;
        if( structIndex < 0 || structIndex >= templateCoord.length ) {
            System.err.println("Bad full structure switch parameter: " + structIndex + ".  Perturbation state ignored." );
            return false;
        }

        double genChi1[] = getGenChi1();//Note the generalized chi1 for each residue so we can restore it after the perturbation

        int section = 0;

        if( resAffected[0] != 0 ) {
            System.arraycopy( templateCoord[structIndex][0], 0, m.actualCoordinates, 0, templateCoord[structIndex][0].length );
            section++;
        }

        restoreResBB( m.residue[resAffected[0]], flexBB.get(structIndex).get(0) );
        int curIndex = -1;
        if( resAffected[0] < m.residue.length - 1)
            curIndex = 3*m.residue[resAffected[0] + 1].atom[0].moleculeAtomNumber;//current index in m.actualCoordinates


        for( int a=1; a<numProteinFlexRes; a++) {

            int molResNum = resAffected[a];

            if( molResNum > resAffected[a-1] + 1 ) { //There is a template section before this residue
                System.arraycopy( templateCoord[structIndex][section], 0, m.actualCoordinates, curIndex, templateCoord[structIndex][section].length  );
                section++;
            }

            //Now change the residue backbone itself
            restoreResBB( m.residue[molResNum], flexBB.get(structIndex).get(a) );

            if(molResNum < m.residue.length - 1)
                curIndex = 3*m.residue[molResNum + 1].atom[0].moleculeAtomNumber;
            else
                curIndex = -1;//We're at the end of the molecule
        }

        //Change the final template section
        if( curIndex != -1)
            System.arraycopy( templateCoord[structIndex][section], 0, m.actualCoordinates, curIndex,  templateCoord[structIndex][section].length );


        if(idealizeSC)
            idealizeSidechains();

        setGenChi1(genChi1);

        return true;

    }


    @Override
    public void undo() {

        if(curParam==0)//Nothing to do
            return;

        doPerturbationMotion(0);
    }


    public void setDefaultParams() {

        minParams = new float[numStructs];
        maxParams = new float[numStructs];

        for(int a=0; a<numStructs; a++)
            minParams[a] = maxParams[a] = a;
    }


    public void writeExtraInfo( BufferedWriter bw ) {
        try {
            bw.append(PDBList);
            bw.newLine();
        } catch(IOException e) {
            System.err.println(e.getMessage());
        }
    }


}
