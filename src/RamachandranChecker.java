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
//	RamachandranChecker.java
//
//	Version:           2.1 beta
//
//
//	  authors:
// 	  initials    name                 organization                email
//	 ---------   -----------------    ------------------------    ----------------------------
//	  MAH           Mark A. Hallen	  Duke University               mah43@duke.edu
///////////////////////////////////////////////////////////////////////////////////////////////
/**
 * This class checks backbone dihedrals against experimental Ramachandran-plot contours.
 *
 */
import java.io.*;
import java.util.StringTokenizer;

public class RamachandranChecker {


    float[][][] tables;//tables[a][b][c] is density at phi=b+/-1, psi=c+/-1
    //with a=0 is for gly, a=1 for pro, a=2 for general, and a=3 for pre-pro

//This class is designed to read the Richardsons' top500 Ramachandran plot density data:
//(from top500-angles/pct/rama)

    //2-degree bins for each angle, centered at -179, -177,..
    //Density (0-to-1 scale) in 3rd column; sorted by phi and then by psi;
    //columns separated by spaces; starts with some comment lines marked with a #



    static float denCutoff = 0.02f;//Cutoff density for being allowed

    private RamachandranChecker() {

    }

    public static RamachandranChecker getInstance() {
        return RamachandranCheckerHolder.INSTANCE;
    }

    private static class RamachandranCheckerHolder {
        private static final RamachandranChecker INSTANCE = new RamachandranChecker();
    }


    public void readInputFiles(String[] fileNames) { //Filenames in order: gly, pro, general, pre-pro

        tables = new float[4][180][180];

        for(int a=0; a<4; a++) {

            try {
                BufferedReader br = new BufferedReader( new FileReader( fileNames[a] ) );
                String line = br.readLine();

                while( line.charAt(0) == '#' )
                    line = br.readLine();

                for(int phiBin=0; phiBin<180; phiBin++) {
                    for(int psiBin=0; psiBin<180; psiBin++) {

                        StringTokenizer st = new StringTokenizer(line," ");

                        st.nextToken();
                        st.nextToken();

                        tables[a][phiBin][psiBin] = Float.valueOf(st.nextToken());

                        line = br.readLine();
                    }
                }

                br.close();
            } catch(IOException e) {
                System.err.println( "Error reading Ramachandran plot file " + fileNames[a] );
                System.err.println(e.getMessage());
                System.exit(1);
            }

        }
    }


    public boolean[] checkByAAType(Molecule m, int resNum) {
        //Return the acceptability at {gly, pro, other AA types} of a given residue's BB dihedrals
        //resNum = molecule residue number in m
        //Looks at actualCoordinates

        boolean ans[] = new boolean[3];

        Residue res = m.residue[resNum];

        if( ! ( m.checkNBonded(resNum) && m.checkCBonded(resNum)
                && m.strand[res.strandNumber].isProtein ) ) { //Phi or psi will not be defined so we cannot rule against them
            for(int a=0; a<3; a++)
                ans[a] = true;
            return ans;
        }

        float phiPsi[] = getPhiPsi(m, resNum);

        for(int a=0; a<3; a++)
            ans[a] = checkAngles(phiPsi[0], phiPsi[1], a);

        return ans;
    }


    //Same but for prePro
    public boolean checkPrePro(Molecule m, int resNum) {

        if(!m.checkNBonded(resNum))
            return true;
        //If we're checking for pre-pro the residue should not be at a C-terminus

        float phiPsi[] = getPhiPsi(m, resNum);

        return checkAngles(phiPsi[0], phiPsi[1], 3);
    }


    //Returns {phi,psi}.  resNum is a molecule residue number.
    public float[] getPhiPsi(Molecule m, int resNum) {

        Residue res = m.residue[resNum];
        float ans[] = new float[2];

        //Get coordinates of relevant atoms
        Atom CLast = new Atom( m.getActualCoord( m.residue[ resNum-1 ].getAtomNameToMolnum("C") ) );
        Atom NCur = new Atom( m.getActualCoord( res.getAtomNameToMolnum("N") ) );
        Atom CACur = new Atom( m.getActualCoord( res.getAtomNameToMolnum("CA") ) );
        Atom CCur = new Atom( m.getActualCoord( res.getAtomNameToMolnum("C") ) );
        Atom NNext = new Atom( m.getActualCoord( m.residue[ resNum+1 ].getAtomNameToMolnum("N") ) );
        ans[0] = (float)CCur.torsion(CLast, NCur, CACur);//phi
        ans[1] = (float)NNext.torsion(NCur, CACur, CCur);//psi

        return ans;
    }



    public boolean checkAngles(float phi, float psi, int plotNum) {

        int phiBin = (int)((phi+180)/2);
        int psiBin = (int)((psi+180)/2);
        float den = tables[plotNum][phiBin][psiBin];
        if(den > denCutoff)
            return true;
        else
            return false;
    }



}
