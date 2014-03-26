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
//	Backrub.java
//
//	Version:           2.1 beta
//
//
//	  authors:
// 	  initials    name                 organization                email
//	 ---------   -----------------    ------------------------    ----------------------------
//	  MAH           Mark A. Hallen	  Duke University               mah43@duke.edu
///////////////////////////////////////////////////////////////////////////////////////////////
import java.util.StringTokenizer;


public class Backrub extends Perturbation {

    //This is a backrub, formulated as a perturbation.  It affects three residues.
    //Like the shear, it moves a section of chain between two alpha carbons,
    //and the implementation is based on the shear's.

    //There are three matrices used: m=0 for the first peptide link, m=1 for the middle CA and sidechain, m=2 for the second peptide link


    static float defaultMaxParams[] = {2.5f};
    static float defaultMinParams[] = {-2.5f};

    
    public Backrub(Molecule molec,int resList[]){

        type="Backrub";
        m=molec;
        resAffected=resList;

    }


    public void setDefaultParams(){
        minParams = defaultMinParams.clone();
        maxParams = defaultMaxParams.clone();
    }



    public boolean doPerturbationMotion(float param){//Apply the perturbation
        //Use an arbitrary param (primary backrub angle in degrees)
        //Don't store rotation matrices or translations

        float rm[][][]=new float[3][3][3];//Rotation matrices


        calcTransRot(param,rm);

        applyTransRot(rm);

        return true;
    }



    public void calcTransRot(float param, float rm[][][]){//calculate rotation matrices for a given perturbation parameter, put them in rm
//rm should be 3X3X3

        Atom[] Calphas=new Atom[3];
        for(int a=0;a<3;a++){
            Calphas[a]=m.residue[resAffected[a]].getAtomByName("CA");
        }

        Atom O1 = m.residue[resAffected[0]].getAtomByName("O");
        Atom O2 = m.residue[resAffected[1]].getAtomByName("O");


        float x[][]=new float[3][3];//Calpha coordinates
        for(int a=0;a<3;a++)
            x[a] = m.getActualCoord(Calphas[a].moleculeAtomNumber);

        float O1Coord[] = m.getActualCoord(O1.moleculeAtomNumber);
        float O2Coord[] = m.getActualCoord(O2.moleculeAtomNumber);

        RotMatrix r=new RotMatrix();//This is actually an object for performing rot. matrix-related calculations
        Backrubs b = new Backrubs();//And this is an object, created for BRDEE, that performs backrub-related calculations.
        //Used for the secondary peptide rotations)
        float rotax[]=r.subtract(x[2],x[0]);//vector from 1st to second Calpha: primary rotation axis

        //Create the primary rotation matrix.  This can be used to rotate about either anchor CA.  
        r.getRotMatrix(rotax[0],rotax[1],rotax[2],param,rm[1]);

        //Get first corrective peptide rotation matrix
        Atom anchor1 = new Atom("CA",x[0][0],x[0][1],x[0][2]);
        Atom anchor2 = new Atom("CA",x[2][0],x[2][1],x[2][2]);
        float midCACoord[] = r.add( r.applyRotMatrix( rm[1], r.subtract( x[1], x[0] ) ), x[0] );
        //Rotated middle CA coordinates
        Atom midCA = new Atom("CA",midCACoord[0],midCACoord[1],midCACoord[2]);
        Atom oldO1 = new Atom("O",O1Coord[0],O1Coord[1],O1Coord[2]);
        Atom oldO2 = new Atom("O",O2Coord[0],O2Coord[1],O2Coord[2]);
        O1Coord = r.add( r.applyRotMatrix( rm[1], r.subtract( O1Coord, x[0] ) ), x[0] );
        O2Coord = r.add( r.applyRotMatrix( rm[1], r.subtract( O2Coord, x[0] ) ), x[0] );
        //Rotated oxygen coordinates
        Atom newO1 = new Atom("O",O1Coord[0],O1Coord[1],O1Coord[2]);
        Atom newO2 = new Atom("O",O2Coord[0],O2Coord[1],O2Coord[2]);

        //We need these atoms to put into getSmallRotAngle, which uses the Atom.coord array
        
        float theta;
        float M[][] = new float[3][3];

        theta = getSmallRotAngle(newO1,midCA,anchor1,oldO1,b);
        theta *= b.thetaSmallScale;
        if (Math.signum(theta)==Math.signum(param))
            theta = -theta;
        rotax = r.subtract( x[1], x[0] );
        r.getRotMatrix( rotax[0], rotax[1], rotax[2], theta, M );
        rm[0] = r.multiplyMatrices( M, rm[1] );


        theta = getSmallRotAngle(newO2,midCA,anchor2,oldO2,b);
        theta *= b.thetaSmallScale;
        if (Math.signum(theta)==Math.signum(param))
            theta = -theta;
        rotax = r.subtract( x[2], x[1] );
        r.getRotMatrix( rotax[0], rotax[1], rotax[2], theta, M );
        rm[2] = r.multiplyMatrices( M, rm[1] );
    }


    public void applyTransRot(float rm[][][]){


        //All the rotation are about anchor-point CAs, so order is not important here like for the shear
        //The second peptide plane is rotated about the ending anchor CA; the other stuff is rotated about the starting one
        int CA0AtNum = m.residue[resAffected[0]].getAtomNameToMolnum("CA");
        int CA2AtNum = m.residue[resAffected[2]].getAtomNameToMolnum("CA");



        //Only rotate the amide group
        int amide[] = m.residue[resAffected[2]].getAtomList(true,false,false,false);
        m.rotateAtomList( amide, rm[2], m.actualCoordinates[3*CA2AtNum],
                m.actualCoordinates[3*CA2AtNum+1], m.actualCoordinates[3*CA2AtNum+2], false);


        //Rotate the carbonyl
        int carbonyl[] = m.residue[resAffected[1]].getAtomList(false, false, false, true);
        m.rotateAtomList( carbonyl, rm[2], m.actualCoordinates[3*CA2AtNum],
                m.actualCoordinates[3*CA2AtNum+1], m.actualCoordinates[3*CA2AtNum+2], false);

        //Apply the primary rotation only to the alpha carbon and sidechain
        int nonAmideCarbonyl[] = m.residue[resAffected[1]].getAtomList(false, true, true, false);
        m.rotateAtomList( nonAmideCarbonyl, rm[1], m.actualCoordinates[3*CA0AtNum],
                m.actualCoordinates[3*CA0AtNum+1], m.actualCoordinates[3*CA0AtNum+2], false);

        //Rotate the amide group
        amide = m.residue[resAffected[1]].getAtomList(true,false,false,false);
        m.rotateAtomList( amide, rm[0], m.actualCoordinates[3*CA0AtNum],
                m.actualCoordinates[3*CA0AtNum+1], m.actualCoordinates[3*CA0AtNum+2], false);


        //Just rotate the carbonyl group
        carbonyl = m.residue[resAffected[0]].getAtomList(false, false, false, true);
        m.rotateAtomList( carbonyl, rm[0], m.actualCoordinates[3*CA0AtNum],
                m.actualCoordinates[3*CA0AtNum+1], m.actualCoordinates[3*CA0AtNum+2], false);

    }


    public float getStepSizeForMinimizer(){//Return a step size for the minimizer to use when optimizing this perturbation
        //This could potentially be improved
        return 0.5f;
    }


    //Get the small rotation angle that will rotate atom pp1 around the axis defined by atoms (pp2,a3), so that pp1 will be as close as possible to atom a4
    //Adapted from Backrubs to use only Atom.coord arrays (needed for the setup of the calculations here)
    private float getSmallRotAngle(Atom pp1, Atom pp2, Atom a3, Atom a4, Backrubs b){

            Atom pp3 = b.projectPointLine(a3, pp2, pp1);
            Atom pp4 = b.projectPointPlane(a3, pp2, pp3, a4);
            Atom closestPoint = b.getClosestPoint(pp3,pp1,pp4);
            return (float)closestPoint.angle(pp1, pp3);
    }


    //This is almost the same as Shear.setParams, but it uses backrub static fields so it's separate
    public static void setParams(String s){

        if(!s.equalsIgnoreCase("none")){
            StringTokenizer st = new StringTokenizer(s);
            int numStates = st.countTokens()/2;

            if( ( numStates*2 != st.countTokens() ) || numStates == 0 ){//odd number of tokens or no parameter values provided
                System.err.println("Badly formulated backrub parameters: using single default state, -2.5 to 2.5.");
                return;
            }


            float min1 = Float.valueOf(st.nextToken());
            float max1 = Float.valueOf(st.nextToken());

            if( min1 + max1 != 0 ){//The first (unperturbed) state must be centered at 0
                System.err.println("First backrub state (unperturbed) must be centered at 0: using -2.5 to 2.5");
                min1 = -2.5f;
                max1 = 2.5f;
            }

            defaultMinParams = new float[numStates];
            defaultMaxParams = new float[numStates];

            defaultMinParams[0] = min1;
            defaultMaxParams[0] = max1;

            for(int state=1; state<numStates; state++){
                defaultMinParams[state] = Float.valueOf(st.nextToken());
                defaultMaxParams[state] = Float.valueOf(st.nextToken());
            }

        }
    }

}
