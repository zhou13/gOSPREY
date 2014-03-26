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
//	Shear.java
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

public class Shear extends Perturbation {

    //This motion comes from KiNG's chiropraxis.mc.CaShear
    //It it mostly found in helices and affects four residues

    //There are three matrices used: m=0 is for the first peptide link and sidechain, m=1 for the middle peptide, m=2 for the last sidechain and peptide
    //The sidechains on the anchor carbons will get their own matrices later
    //There are two translation vectors (first peptide is not translated): m=0 for 1st sidechain & middle peptide, m=1 for last sidechain & peptide

    public Atom Calphas[];//The four Calphas
    public Atom midC, midO;//the carbonyl carbon and oxygen of the middle peptide link



    static float defaultMaxParams[] = {2.5f};
    static float defaultMinParams[] = {-2.5f};

    public Shear(Molecule molec,int resList[]){

        type="Shear";
        m=molec;
        resAffected=resList;
    }


    public void setDefaultParams(){
        minParams = defaultMinParams.clone();
        maxParams = defaultMaxParams.clone();
    }




    public boolean doPerturbationMotion(float param){//Apply the perturbation
        //Use an arbitrary param (primary shear angle in degrees)

        float rm[][][]=new float[3][3][3];//single-use rotation matrices
        float tr[][]=new float[2][3];//single-use translation vectors


        Calphas=new Atom[4];
        for(int a=0;a<4;a++){
            Calphas[a]=m.residue[resAffected[a]].getAtomByName("CA");
        }

        Residue theRes=m.residue[resAffected[1]];
        midC = theRes.getAtomByName("C");
        midO = theRes.getAtomByName("O");


        calcTransRot(param,rm,tr);

        applyTransRot(rm,tr);

        return true;
    }



    public void calcTransRot(float param, float rm[][][], float tr[][]){//calculate rotation matrices and translation vectors for a given perturbation parameter, put them in rm and tr
//rm and tr should be 3X3X3 and 2X3 arrays respectively

        float x[][]=new float[4][3];//Calpha coordinates
        for(int a=0;a<4;a++){
            int CAAtNum = Calphas[a].moleculeAtomNumber;
            for(int b=0;b<3;b++)
                x[a][b]=m.actualCoordinates[3*CAAtNum + b];
        }

        float w[] = new float[3];
        for(int a=0;a<3;a++)
            w[a] = m.actualCoordinates[3*midO.moleculeAtomNumber + a] - m.actualCoordinates[3*midC.moleculeAtomNumber + a];


        RotMatrix r=new RotMatrix();//This is actually an object for performing rot. matrix-related calculations
        float x01[]=r.subtract(x[1],x[0]);//vector from 1st to second Calpha
        float x12[]=r.subtract(x[2],x[1]);
        float x23[]=r.subtract(x[3],x[2]);

        float rotax1[]=r.cross(x01,x12);//First peptide and sidechain
        r.getRotMatrix(rotax1[0],rotax1[1],rotax1[2],param,rm[0]);//Note getRotMatrix takes angles in degrees

        tr[0]=r.subtract(r.applyRotMatrix(rm[0],x01),x01);

        float y1[]=r.add(x[1],tr[0]);//Last peptide and sidechain 2.  y1=new 2nd CA coordinates
        float y13[]=r.subtract(x[3],y1);
        float a=r.norm(x23);
        float b=r.norm(y13);
        float d=r.norm(x12);
        float beta=(float)(Math.acos(r.dot(x23,y13)/(a*b)) - Math.acos((b*b+a*a-d*d)/(2*a*b)));//Second rotation angle
        float rotax2[]=r.cross(x23,y13);
        r.getRotMatrixRad(rotax2[0],rotax2[1],rotax2[2],beta,rm[2]);//beta is in radians
        tr[1]=r.subtract(x23,r.applyRotMatrix(rm[2],x23));


        float y2[]=r.add(x[2],tr[1]);//Calculating middle peptide rotation, rm[1].  y2=new 3rd CA coordinates
        float y12[]=r.subtract(y2,y1);
        float theta=r.getAngle(x12,y12);
        float srotax[]=r.cross(x12,y12);//Rotation axis to superimpose

        if( r.norm(srotax) == 0 )//This will happen if x12, y12 are the same.  No feasible shear will have x12 and y12 in opposite directions
            rm[1] = r.identity();
        else
            r.getRotMatrixRad(srotax[0],srotax[1],srotax[2],theta,rm[1]);//Rotate middle peptide to superimpose Calphas
        float u[]=r.applyRotMatrix(rm[1],w);
        
        float rmalpha[][]=new float[3][3];//Additional rotation to correct carbonyl orientation
        float vhat[]=r.scale(y12,1/r.norm(y12));//unit vector along line from 2nd to 3rd Calpha
        float alpha=(float)Math.atan2( r.dot(r.cross(u,w),vhat), r.dot(w,u)-r.dot(w,vhat)*r.dot(u,vhat));
        r.getRotMatrixRad(vhat[0],vhat[1],vhat[2],alpha,rmalpha);

        rm[1]=r.multiplyMatrices(rmalpha,rm[1]);

    }



    public void applyTransRot(float rm[][][], float tr[][]){

        //Movements are applied in reverse of the atom order in the peptide chain because the motion
        //of some atoms is based on the original position of CAs earlier in the chain
        //(For the same reason, if some of the residues do not move because they're
        //not flexible, this does not affect the motion of other residues)


            int CA2AtNum = Calphas[2].moleculeAtomNumber;//Just rotate and translate the amide group
            int amide[] = m.residue[resAffected[3]].getAtomList(true,false,false,false);
            m.rotateAtomList( amide, rm[2], m.actualCoordinates[3*CA2AtNum],
                    m.actualCoordinates[3*CA2AtNum+1], m.actualCoordinates[3*CA2AtNum+2], false);
            m.translateAtomList( amide, tr[1], false, false);


            int nonAmide[] = m.residue[resAffected[2]].getAtomList(false, true, true, true);
            m.rotateAtomList( nonAmide, rm[2], m.actualCoordinates[3*CA2AtNum],
                    m.actualCoordinates[3*CA2AtNum+1], m.actualCoordinates[3*CA2AtNum+2], false);
            m.translateAtomList( nonAmide, tr[1], false, false);

            int CA1AtNum = Calphas[1].moleculeAtomNumber;//Rotate and translate the amide group
            amide =  m.residue[resAffected[2]].getAtomList(true,false,false,false);
            m.rotateAtomList( amide, rm[1], m.actualCoordinates[3*CA1AtNum],
                    m.actualCoordinates[3*CA1AtNum+1], m.actualCoordinates[3*CA1AtNum+2], false);
            m.translateAtomList( amide, tr[0], false, false);


            int carbonyl[] = m.residue[resAffected[1]].getAtomList(false,false,false,true);
            m.rotateAtomList( carbonyl, rm[1], m.actualCoordinates[3*CA1AtNum],
                    m.actualCoordinates[3*CA1AtNum+1], m.actualCoordinates[3*CA1AtNum+2], false);
            m.translateAtomList( carbonyl, tr[0], false, false);

            int CA0AtNum = Calphas[0].moleculeAtomNumber;//Rotate all the residue except the carbonyl group about the first anchor CA
            int nonCarbonyl[] = m.residue[resAffected[1]].getAtomList(true,true,true,false);
            m.rotateAtomList( nonCarbonyl, rm[0], m.actualCoordinates[3*CA0AtNum],
                    m.actualCoordinates[3*CA0AtNum+1], m.actualCoordinates[3*CA0AtNum+2], false);

        

            m.rotateAtomList( m.residue[resAffected[0]].getAtomList(false, false, false, true), rm[0],
                    m.actualCoordinates[3*CA0AtNum], m.actualCoordinates[3*CA0AtNum+1],
                    m.actualCoordinates[3*CA0AtNum+2], false);
    }



    public float getStepSizeForMinimizer(){//Return a step size for the minimizer to use when optimizing this perturbation
        //This could potentially be improved
        return 0.5f;
    }


    //This is almost the same as Backrub.setParams, but it uses shear static fields so it's separate
    public static void setParams(String s){

        if(!s.equalsIgnoreCase("none")){
            StringTokenizer st = new StringTokenizer(s);
            int numStates = st.countTokens()/2;

            if( ( numStates*2 != st.countTokens() ) || numStates == 0 ){//odd number of tokens or no parameter values provided
                System.err.println("Badly formulated shear parameters: using single default state, -2.5 to 2.5.");
                return;
            }


            float min1 = Float.valueOf(st.nextToken());
            float max1 = Float.valueOf(st.nextToken());

            if( min1 + max1 != 0 ){//The first (unperturbed) state must be centered at 0
                System.err.println("First shear state (unperturbed) must be centered at 0: using -2.5 to 2.5");
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

   /*
    To help make default states:
    public void fillDefaultParams(){//Fill in default parameter values and where they apply
        //0 and +/-7 degrees are allowable defaults for all mutations and rotamers
        //Also +/- 3 degrees if there are no overlapping perturbations and no AA type in the shear has >10 rotamers
        
    }*/

}
