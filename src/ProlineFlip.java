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
//	ProlineFlip.java
//
//	Version:           2.1 beta
//
//
//	  authors:
// 	  initials    name                 organization                email
//	 ---------   -----------------    ------------------------    ----------------------------
//	  MAH           Mark A. Hallen	  Duke University               mah43@duke.edu
///////////////////////////////////////////////////////////////////////////////////////////////

import java.util.ArrayList;

//This discrete perturbation flips the pucker of proline
//It should be applied to all flexible prolines but not elsewhere

//Note: perturbations should treat the amide nitrogen and the CD and its hydrogens together
//as a rigid body, as they would the amide hydorgen and nitrogen
//The rest of the sidechain is treated as usual, except during idealization
//and we not attempt to maintain constant chi1 when perturbing
//(chi1 follows from phi and bond angles according to Ho et al's formula)

public class ProlineFlip extends Perturbation {

    static boolean UP = true;
    static boolean DOWN = false;//Pucker options


    public ProlineFlip(Molecule molec,int resList[]) {

        type="PROLINE FLIP";
        m=molec;
        resAffected=resList;
    }


    public boolean doPerturbationMotion(float param) {

        Residue res = m.residue[resAffected[0]];

        if( ( param == 1 ) && (res.name.equalsIgnoreCase("PRO")) ) { //Do the flip

            res.pucker = !(res.pucker);//Change the pucker
            return m.idealizeProRing(res);//Idealize the ring with the new pucker
        }

        return true;
    }


    public void setDefaultParams() {
        float mp[] = {0f,1f};
        minParams = mp.clone();
        maxParams = mp.clone();
    }



    public static ProlineFlip[] generateAll( Molecule m, StrandRotamers strandRot[] ) {
        //Generate all the proline flips possible for a given molecule at the positions where proline is allowable
        //strandRot is an array of StrandRotamers objects for the strands in m
        ArrayList<Integer> resList = new ArrayList<Integer>();

        for( int resNum=0; resNum<m.residue.length; resNum++ ) {
            Residue res = m.residue[resNum];

            if( strandRot[res.strandNumber].checkAllowable( res.strandResidueNumber , "PRO" ) )
                resList.add(resNum);
        }

        ProlineFlip[] ans = new ProlineFlip[resList.size()];

        for( int a=0; a<resList.size(); a++ ) {
            int aff[] = { resList.get(a) };
            ProlineFlip pf = new ProlineFlip( m, aff );
            pf.setDefaultParams();
            ans[a] = pf;
        }

        return ans;
    }


}
