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

////////////////////////////////////////////////////////////////////////////////////////////
// RotMatrix.java
//
//  Version:           2.1 beta
//
//
// authors:
//    initials    name            organization                email
//   ---------   --------------  ------------------------    ------------------------------
//     RHL        Ryan Lilien     Dartmouth College           ryan.lilien@dartmouth.edu
//     MAH        Mark A. Hallen	Duke University         mah43@duke.edu
//
////////////////////////////////////////////////////////////////////////////////////////////

/**
 * This class implements rotation matricies.
 * Written by: Ryan Lilien  (2001-2004)
 */

/**
 * This class implements rotation matricies.
 */
public class RotMatrix
{


	RotMatrix()
	{

	}

	// Rotates an array of points around the +x axis
	//  using a right handed rotation of thetaDeg degrees
	// theCoords are the coordinates where the ith coodinates
	//  are in position 3*i .. (3*(i+1))-1
	public void xAxisRotate(float thetaDeg, float theCoords[],int numCoords) {

		axisRotate(1.0f,0.0f,0.0f,thetaDeg,theCoords,numCoords);
	}

	// Rotates an array of points around the +x axis
	//  using a right handed rotation of thetaDeg degrees
	// theCoords are the coordinates where the ith coodinates
	//  are in position 3*i .. (3*(i+1))-1
	public void xAxisRotate(double thetaDeg, float theCoords[],int numCoords) {

		axisRotate(1.0f,0.0f,0.0f,new Double(thetaDeg).floatValue(),theCoords,numCoords);
	}

	// Rotates an array of points around the +y axis
	//  using a right handed rotation of thetaDeg degrees
	// theCoords are the coordinates where the ith coodinates
	//  are in position 3*i .. (3*(i+1))-1
	public void yAxisRotate(float thetaDeg, float theCoords[],int numCoords) {

		axisRotate(0.0f,1.0f,0.0f,thetaDeg,theCoords,numCoords);
	}

	// Rotates an array of points around the +y axis
	//  using a right handed rotation of thetaDeg degrees
	// theCoords are the coordinates where the ith coodinates
	//  are in position 3*i .. (3*(i+1))-1
	public void yAxisRotate(double thetaDeg, float theCoords[],int numCoords) {

		axisRotate(0.0f,1.0f,0.0f,new Double(thetaDeg).floatValue(),theCoords,numCoords);
	}

	// Rotates an array of points around the +z axis
	//  using a right handed rotation of thetaDeg degrees
	// theCoords are the coordinates where the ith coodinates
	//  are in position 3*i .. (3*(i+1))-1
	public void zAxisRotate(float thetaDeg, float theCoords[],int numCoords) {

		axisRotate(0.0f,0.0f,1.0f,thetaDeg,theCoords,numCoords);
	}

	// Rotates an array of points around the +z axis
	//  using a right handed rotation of thetaDeg degrees
	// theCoords are the coordinates where the ith coodinates
	//  are in position 3*i .. (3*(i+1))-1
	public void zAxisRotate(double thetaDeg, float theCoords[],int numCoords) {

		axisRotate(0.0f,0.0f,1.0f,new Double(thetaDeg).floatValue(),theCoords,numCoords);
	}


	// Rotates an array of points around the axis ax, ay, az
	//  using a right handed rotation of thetaDeg degrees
	// theCoords are the coordinates where the ith coodinates
	//  are in position 3*i .. (3*(i+1))-1
	public void axisRotate(float ax, float ay, float az, float thetaDeg,
	 float theCoords[],int numCoords) {

			float tx,ty,tz;

			float[][] rot_mtx = new float[3][3];
			getRotMatrix(ax,ay,az,(float) thetaDeg,rot_mtx);

			int ix3 = 0;
			for(int i=0;i<numCoords;i++){
				tx=theCoords[ix3];
				ty=theCoords[ix3+1];
				tz=theCoords[ix3+2];

				theCoords[ix3] = tx * rot_mtx[0][0] + ty * rot_mtx[0][1] + tz * rot_mtx[0][2];
				theCoords[ix3+1] = tx * rot_mtx[1][0] + ty * rot_mtx[1][1] + tz * rot_mtx[1][2];
				theCoords[ix3+2] = tx * rot_mtx[2][0] + ty * rot_mtx[2][1] + tz * rot_mtx[2][2];

				ix3 += 3;
			}
	}

	// Translates an array of points by the specified amount
	public void translate(float tx, float ty, float tz, float theCoords[],
	 int numCoords) {

		int ix3 = 0;
		for(int i=0;i<numCoords;i++){
			theCoords[ix3] += tx;
			theCoords[ix3+1] += ty;
			theCoords[ix3+2] += tz;
			ix3 += 3;
		}
	}

	// This function constructs a rotation matrix from a rotation in
	//  axis-angle notation
	public void getRotMatrix(float fx, float fy, float fz, float angle,
		float[][] rot_mtx) {

		// First convert the axisangle to a quaternion
		float sin_a = (float) Math.sin(angle*3.14159265/180/2);
		float cos_a = (float) Math.cos(angle*3.14159265/180/2);
		float tmp = (float) Math.sqrt(fx*fx + fy*fy + fz*fz);
		float qx = fx / tmp * sin_a;
		float qy = fy / tmp * sin_a;
		float qz = fz / tmp * sin_a;
		float qw = cos_a;
		tmp = (float) Math.sqrt(qx*qx + qy*qy + qz*qz + qw*qw);
		qx /= tmp;
		qy /= tmp;
		qz /= tmp;
		qw /= tmp;
		float xx = qx * qx;
		float xy = qx * qy;
		float xz = qx * qz;
		float xw = qx * qw;

		float yy = qy * qy;
		float yz = qy * qz;
		float yw = qy * qw;

		float zz = qz * qz;
		float zw = qz * qw;

		rot_mtx[0][0] = 1 - 2 * (yy + zz);
		rot_mtx[0][1] = 2 * (xy - zw);
		rot_mtx[0][2] = 2 * (xz + yw);

		rot_mtx[1][0] = 2 * (xy + zw);
		rot_mtx[1][1] = 1 - 2 * (xx + zz);
		rot_mtx[1][2] = 2 * (yz - xw);

		rot_mtx[2][0] = 2 * (xz - yw);
		rot_mtx[2][1] = 2 * (yz + xw);
		rot_mtx[2][2] = 1 - 2 * (xx + yy);
	}


        public void getRotMatrixRad(float fx, float fy, float fz, float angle,
		float[][] rot_mtx) {//Same as getRotMatrix but takes the angle in radians

            getRotMatrix(fx,fy,fz,angle*180/(float)Math.PI,rot_mtx);
        }


        float[] applyRotMatrix(float[][] rm, float[] vec){//Apply rotation matrix rm to vector vec, i.e. compute rm*vec
            float ans[]=new float[3];
            float val;
            for(int a=0;a<3;a++){
                val=0;
                for(int b=0;b<3;b++){
                    val+=rm[a][b]*vec[b];
                }
                ans[a]=val;
            }
            return ans;
        }


        //Double-precision versions
        // This function constructs a rotation matrix from a rotation in
	//  axis-angle notation
	public void getRotMatrix(double fx, double fy, double fz, double angle,
		double[][] rot_mtx) {

		// First convert the axisangle to a quaternion
		double sin_a = Math.sin(angle*3.14159265/180/2);
		double cos_a = Math.cos(angle*3.14159265/180/2);
		double tmp = Math.sqrt(fx*fx + fy*fy + fz*fz);
		double qx = fx / tmp * sin_a;
		double qy = fy / tmp * sin_a;
		double qz = fz / tmp * sin_a;
		double qw = cos_a;
		tmp = Math.sqrt(qx*qx + qy*qy + qz*qz + qw*qw);
		qx /= tmp;
		qy /= tmp;
		qz /= tmp;
		qw /= tmp;
		double xx = qx * qx;
		double xy = qx * qy;
		double xz = qx * qz;
		double xw = qx * qw;

		double yy = qy * qy;
		double yz = qy * qz;
		double yw = qy * qw;

		double zz = qz * qz;
		double zw = qz * qw;

		rot_mtx[0][0] = 1 - 2 * (yy + zz);
		rot_mtx[0][1] = 2 * (xy - zw);
		rot_mtx[0][2] = 2 * (xz + yw);

		rot_mtx[1][0] = 2 * (xy + zw);
		rot_mtx[1][1] = 1 - 2 * (xx + zz);
		rot_mtx[1][2] = 2 * (yz - xw);

		rot_mtx[2][0] = 2 * (xz - yw);
		rot_mtx[2][1] = 2 * (yz + xw);
		rot_mtx[2][2] = 1 - 2 * (xx + yy);
	}


        public void getRotMatrixRad(double fx, double fy, double fz, double angle,
		double[][] rot_mtx) {//Same as getRotMatrix but takes the angle in radians

            getRotMatrix(fx,fy,fz,angle*180/Math.PI,rot_mtx);
        }


        double[] applyRotMatrix(double[][] rm, double[] vec){//Apply rotation matrix rm to vector vec, i.e. compute rm*vec
            double ans[]=new double[3];
            double val;
            for(int a=0;a<3;a++){
                val=0;
                for(int b=0;b<3;b++){
                    val+=rm[a][b]*vec[b];
                }
                ans[a]=val;
            }
            return ans;
        }



        float[] unapplyRotMatrix(float[][] rm, float[] vec){//Reverse rotation specified by matrix rm on vector vec, i.e. compute inv(rm)*vec=transpose(rm)*vec
            float ans[]=new float[3];
            float val;
            for(int a=0;a<3;a++){
                val=0;
                for(int b=0;b<3;b++){
                    val+=rm[b][a]*vec[b];
                }
                ans[a]=val;
            }
            return ans;
        }



        float[][] multiplyMatrices(float[][] M1, float[][] M2){
            float[][] ans=new float[3][3];
            for(int a=0;a<3;a++){
                for(int b=0;b<3;b++){
                    ans[a][b]=0;
                    for(int c=0;c<3;c++){
                        ans[a][b]+=M1[a][c]*M2[c][b];
                    }
                }
            }
            return ans;
        }


        float[][] getSuperposingRotMatrix(float uold[], float unew[], float vold[], float vnew[]){
        //Returns a rotation matrix that rotates vector uold to point in the direction of unew, and vold to point in the direction of vnew
        //Relies on uold-vold and unew-vnew angles being basically the same (no rotation exactly satisfies the requirements if they are different)
        //If the angles are close to equal the rotation matrix will exactly superimpose uold, unew and be close to right for vold, vnew

            float th;

            //First create a matrix rotating uold to unew
            //Axis for this will be uold X unew (sign for this will give an angle < 180 degrees)
            float mtx1[][]=new float[3][3];
            float axis1[] = cross(uold, unew);
            if( norm(axis1) == 0 ){//uold and unew are collinear
                if( dot(uold,unew) > 0 )//uold and unew are already superimposed
                    mtx1 = identity();
                else{//Need a 180-degree rotation.
                    float normal[] = getPerpendicular(uold);
                    getRotMatrix( normal[0], normal[1], normal[2], 180, mtx1 );
                }
            }
            else{
                th = getAngle(uold,unew);//angle of this first rotation
                getRotMatrixRad(axis1[0], axis1[1], axis1[2], th, mtx1);
            }

            //Now create a matrix to rotate about unew, rotating mtx1*vold to vnew
            //The angle assumption comes in here, because we are using unew as the axis
            //So we need vnew-unew angle = applyMatrix(mtx1,vold)-unew angle (which = vold-uold angle)
            float mtx2[][]=new float[3][3];
            float w1[] = perpendicularComponent( applyRotMatrix(mtx1,vold), unew );
            float w2[] = perpendicularComponent( vnew, unew );
            th = getAngle(w1,w2);//angle of second rotation
            if( dot( cross(w1,w2), unew) > 0)
                getRotMatrixRad( unew[0], unew[1], unew[2], th, mtx2 );
            else
                getRotMatrixRad( unew[0], unew[1], unew[2], -th, mtx2 );

            return multiplyMatrices(mtx2, mtx1);//Return the product of these two rotations

        }


        public float[][] identity(){//Return the 3X3 identity matrix

            float M[][] = new float[3][3];

            for(int a=0;a<3;a++){
                for(int b=0;b<3;b++){
                    if(a==b)
                        M[a][b] = 1;
                    else
                        M[a][b] = 0;
                }
            }

            return M;
        }




        //A couple of quick vector operations, all in 3-D like the above

        float dot(float[] vec1, float[] vec2){
            return vec1[0]*vec2[0]+vec1[1]*vec2[1]+vec1[2]*vec2[2];
        }

        float norm(float[] vec){
            return (float)Math.sqrt(dot(vec,vec));
        }

        float normsq(float[] vec){
            return dot(vec, vec);
        }

        float[] add(float[] vec1,float[] vec2){
            float ans[]=new float[3];
            for(int a=0;a<3;a++){
                ans[a]=vec1[a]+vec2[a];
            }
            return ans;
        }

        float[] subtract(float[] vec1,float[] vec2){
            float ans[]=new float[3];
            for(int a=0;a<3;a++){
                ans[a]=vec1[a]-vec2[a];
            }
            return ans;
        }

        float[] scale(float[] vec, float scalar){
            float ans[]=new float[3];
            for(int a=0;a<3;a++){
                ans[a]=vec[a]*scalar;
            }
            return ans;
        }

        float[] cross(float[] vec1,float[] vec2){
            float ans[]=new float[3];
            ans[0]=vec1[1]*vec2[2]-vec2[1]*vec1[2];
            ans[1]=vec1[2]*vec2[0]-vec2[2]*vec1[0];
            ans[2]=vec1[0]*vec2[1]-vec2[0]*vec1[1];
            return ans;
        }


        float[] parallelComponent(float[] vec1, float[] vec2){//Component of vec1 parallel to vec2
            return scale( vec2, dot(vec1,vec2) / (norm(vec2)*norm(vec2)) );
        }

        float[] perpendicularComponent(float[] vec1, float[] vec2){//Component of vec1 perpendicular to vec2
            return subtract(vec1, parallelComponent(vec1,vec2) );
        }


        float getAngle(float vec1[], float vec2[]){//Get the angle, in radians, between two vectors

            float costh = dot(vec1,vec2) / ( norm(vec1) * norm(vec2) );
            if( costh > 1)//It might be slightly over due to numerical error...this means the angle is basically 0
                costh = 1;
            else if(costh < -1 )
                costh = -1;

            return (float) Math.acos( costh );
        }

        float getAngle(float A[], float B[], float C[]){//Get the angle ABC
            float BA[] = subtract(A,B);
            float BC[] = subtract(C,B);
            return getAngle(BA,BC);
        }

        float[] getPerpendicular(float vec[]){//Generate a vector perpendicular to the argument
            if( vec[1]==0 && vec[2] == 0 ){
                float yhat[] = {0,1,0};
                return cross(vec, yhat);
            }
            else{
                float xhat[] = {1,0,0};
                return cross(vec, xhat);
            }
        }

         public float[] average(float v1[], float v2[]){//Average two vectors
             return scale( add(v1,v2), 0.5f );
         }

         //Some double-precision versions
         double dot(double[] vec1, double[] vec2){
            return vec1[0]*vec2[0]+vec1[1]*vec2[1]+vec1[2]*vec2[2];
        }

        double norm(double[] vec){
            return Math.sqrt(dot(vec,vec));
        }

        double normsq(double[] vec){
            return dot(vec, vec);
        }

        double[] add(double[] vec1,double[] vec2){
            double ans[]=new double[3];
            for(int a=0;a<3;a++){
                ans[a]=vec1[a]+vec2[a];
            }
            return ans;
        }

        double[] subtract(double[] vec1,double[] vec2){
            double ans[]=new double[3];
            for(int a=0;a<3;a++){
                ans[a]=vec1[a]-vec2[a];
            }
            return ans;
        }

        double[] scale(double[] vec, double scalar){
            double ans[]=new double[3];
            for(int a=0;a<3;a++){
                ans[a]=vec[a]*scalar;
            }
            return ans;
        }

        double[] cross(double[] vec1,double[] vec2){
            double ans[]=new double[3];
            ans[0]=vec1[1]*vec2[2]-vec2[1]*vec1[2];
            ans[1]=vec1[2]*vec2[0]-vec2[2]*vec1[0];
            ans[2]=vec1[0]*vec2[1]-vec2[0]*vec1[1];
            return ans;
        }


        double getAngle(double vec1[], double vec2[]){//Get the angle, in radians, between two vectors

            double costh = dot(vec1,vec2) / ( norm(vec1) * norm(vec2) );
            if( costh > 1)//It might be slightly over due to numerical error...this means the angle is basically 0
                costh = 1;
            else if(costh < -1 )
                costh = -1;

            return Math.acos( costh );
        }




        /**
         * From KiNG's driftwood.r3.Builder:
    * Given three points A, B, and C,
    * construct a line segment from C to D
    * of length len
    * at angle ang to BC (in degrees, 0 to 180)
    * and with a dihedral angle dihe to ABC (in degrees)
    * return D
    * Used in sidechain idealization
    */
    public float[] get4thPoint(float[] a, float[] b, float[] c, float len, float ang, float dihe)
    {
        float d[] = subtract(b,c);
        d = scale(d, len/norm(d) );

        // Not robust to a/b/c colinear
        // Doesn't matter since that makes dihe undef.
        float x1[] = subtract(a,b);
        float x2[] = subtract(c,b);
        x1 = cross(x1,x2);

        float[][] rot1  = new float[3][3];
        getRotMatrix(x1[0], x1[1], x1[2], ang, rot1);
        d = applyRotMatrix(rot1,d);

        getRotMatrix(x2[0], x2[1], x2[2], dihe, rot1);
        d = applyRotMatrix(rot1,d);

        return add(d,c);
    }

}
