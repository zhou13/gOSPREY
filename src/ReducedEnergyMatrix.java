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
//	ReducedEnergyMatrix.java
//
//	Version:           2.1 beta
//
//
//	  authors:
// 	  initials    name                 organization                email
//	 ---------   -----------------    ------------------------    ----------------------------
//	  MAH           Mark A. Hallen	  Duke University               mah43@duke.edu
///////////////////////////////////////////////////////////////////////////////////////////////
//This class represents the energy bounds for unpruned rotamers in a compact format
//for use in A*
public class ReducedEnergyMatrix {

    private float RedEmat[][];
    // the min energy matrix: the last column contains the intra-energy for each rotamer; the last row
    // contains the shell-residue energy for each rotamer


    int numTotalNodes;//the total number of possible rotamers for the given mutation


    public float[][] getRedEmat() {
        return RedEmat;
    }

    public ReducedEnergyMatrix(float[][] a) {
        RedEmat = a;
    }

    public float getPairwiseE(int index1, int index2) {
        return RedEmat[index1][index2];
    }

    public float getIntraE(int index) { //the intra-energy is in the last column
        return RedEmat[index][numTotalNodes];
    }

    public float getShellRotE(int index) {
        float shlRotE = 0.0f;
        //Skip the first energy which is the intra-rot energy still
        for(int i=numTotalNodes; i<RedEmat.length; i++) {
            shlRotE += RedEmat[i][index];
        }
        return shlRotE;
    }

    public float getIntraAndShellE(int index) {
        return getIntraE(index) + getShellRotE(index);
    }


}
