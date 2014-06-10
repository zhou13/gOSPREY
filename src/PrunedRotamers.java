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
//	PrunedRotamers.java
//
//	Version:           2.1 beta
//
//
//	  authors:
// 	  initials    name                 organization                email
//	 ---------   -----------------    ------------------------    ----------------------------
//	  KER        Kyle E. Roberts       Duke University         ker17@duke.edu
///////////////////////////////////////////////////////////////////////////////////////////////

import java.io.Serializable;
import java.util.Iterator;
import java.util.StringTokenizer;


public class PrunedRotamers<T> implements Iterable<RotInfo<T>>, Serializable {

    private T[][][] prunedRot;

    public PrunedRotamers(int numMutable, int[][] strandMut, RotamerSearch rs, T initVal) {
        prunedRot = initializePrunedRot(numMutable, strandMut, rs, initVal);
    }

    public PrunedRotamers(PrunedRotamers<?> pr, Object initVal) {

        prunedRot = initFromCopy(pr.prunedRot,(T) initVal);


    }

    public T [][][] initializePrunedRot(int numMutable, int[][] strandMut, RotamerSearch rs, T initVal) {

        T prunedRot[][][] = null;
        int numPos = numMutable;
        /*if (ligPresent)
        	numPos++;*/
        prunedRot = (T[][][]) new Object[numPos][][];
        int p1 = 0;
        for (int str1=0; str1<strandMut.length; str1++) {
            for (int i=0; i<strandMut[str1].length; i++) {
                prunedRot[p1] = (T[][]) new Object[rs.strandRot[str1].rl.getNumAAallowed()][];
                for (int a1=0; a1<rs.strandRot[str1].getNumAllowable(strandMut[str1][i]); a1++) {
                    int curAAind1 = rs.strandRot[str1].getIndexOfNthAllowable(strandMut[str1][i],a1);

                    int numRot1;
                    if( rs.doPerturbations )
                        numRot1 = ((StrandRCs)rs.strandRot[str1]).getNumRCs( strandMut[str1][i], curAAind1 );
                    else {
                        numRot1 = rs.strandRot[str1].rl.getNumRotForAAtype(curAAind1);
                        if (numRot1==0) //ALA or GLY
                            numRot1 = 1;
                    }

                    prunedRot[p1][curAAind1] = (T[]) new Object[numRot1];
                    for (int r1=0; r1<numRot1; r1++) {
                        prunedRot[p1][curAAind1][r1] = initVal;

                    }
                }
                p1++;
            }
        }


        return prunedRot;
    }

    //Returns a new independent six-dimensional matrix that is a copy of fromMatrix[][][][][][]
    public T [][][] initFromCopy(Object fromMatrix[][][], T initVal) {

        if (fromMatrix==null)
            return null;

        T toMatrix[][][] = (T[][][]) new Object[fromMatrix.length][][];
        for (int p1=0; p1<toMatrix.length; p1++) {
            if (fromMatrix[p1]!=null) {
                toMatrix[p1] = (T[][]) new Object[fromMatrix[p1].length][];
                for (int a1=0; a1<toMatrix[p1].length; a1++) {
                    if (fromMatrix[p1][a1]!=null) {
                        toMatrix[p1][a1] = (T[]) new Object[fromMatrix[p1][a1].length];
                        for (int r1=0; r1<toMatrix[p1][a1].length; r1++) {
                            toMatrix[p1][a1][r1] = initVal;
                        }
                    }
                }
            }
        }

        return toMatrix;
    }

    public Iterator<RotInfo<T>> iterator() {
        return new PrunedRotIterator<T>(prunedRot);
    }

    public Iterator<RotInfo<T>> iterator(int startRes) {
        return new PrunedRotIterator<T>(prunedRot, startRes);
    }

    public T get(int curPos, int curAA, int curRot) {
        return prunedRot[curPos][curAA][curRot];
    }

    public void set(int curPos, int curAA, int curRot, T val) {
        prunedRot[curPos][curAA][curRot] = val;
    }

    public T get(RotInfo<?> ri) {
        return prunedRot[ri.curPos][ri.curAA][ri.curRot];
    }

    public void set(RotInfo<?> ri, T val) {
        prunedRot[ri.curPos][ri.curAA][ri.curRot] = val;
    }

    public T get(Index3 i) {
        return prunedRot[i.pos][i.aa][i.rot];
    }

    public void set(Index3 i, T val) {
        prunedRot[i.pos][i.aa][i.rot] = val;
    }

}

class PrunedRotIterator<T> implements Iterator<RotInfo<T>> {

    private final T[][][] prunedRot;
    private boolean hasNextItem = false;
    private RotInfo<T> nextItem = null;

    public PrunedRotIterator(T[][][] pr) {
        prunedRot = pr;

        if(prunedRot != null) {
            hasNextItem = true;
            int aaCtr=0;
            while(prunedRot[0][aaCtr] == null) {
                aaCtr++;
            }
            nextItem = new RotInfo<T>(0,aaCtr,0,prunedRot[0][aaCtr][0]);
        } else {
            hasNextItem = false;
        }

    }

    public PrunedRotIterator(T[][][] pr, int startPos) {
        prunedRot = pr;

        if(prunedRot != null) {
            hasNextItem = true;
            int aaCtr=0;
            while(prunedRot[startPos][aaCtr] == null) {
                aaCtr++;
            }
            nextItem = new RotInfo<T>(startPos,aaCtr,0,prunedRot[startPos][aaCtr][0]);
        } else {
            hasNextItem = false;
        }

    }

    public boolean hasNext() {
        return hasNextItem;
    }

    public RotInfo<T> next() {
        RotInfo<T> ret = nextItem;
        calcNext(ret);
        return ret;
    }

    public void calcNext(RotInfo<T> curItem) {
        hasNextItem = false;

        int[] ctr = {curItem.curPos,curItem.curAA,curItem.curRot};
        int[] max = {prunedRot.length,prunedRot[curItem.curPos].length,prunedRot[curItem.curPos][curItem.curAA].length};
        incrementCtr(ctr,max);
        if(ctr[0] == prunedRot.length)
            return;
        while(prunedRot[ctr[0]][ctr[1]] == null) {
            max[2] = 1;
            int tmpAA = ctr[1];
            incrementCtr(ctr,max);
            //Gone past the end of the matrix
            if(ctr[0] == prunedRot.length)
                return;

            max[1] = prunedRot[ctr[0]].length;
        }

        hasNextItem = true;
        nextItem = new RotInfo<T>(ctr[0],ctr[1],ctr[2],prunedRot[ctr[0]][ctr[1]][ctr[2]]);

        //First try to add to the rotamer, then the aa, then the position
        /*if(curItem.curRot < prunedRot[curItem.curPos][curItem.curAA].length-1){
        	nextItem = new RotInfo<T>(curItem.curPos,curItem.curAA,curItem.curRot+1,prunedRot[curItem.curPos][curItem.curAA][curItem.curRot+1]);
        	hasNextItem = true;
        }
        else if(curItem.curAA < prunedRot[curItem.curPos].length-1){
        	nextItem = new RotInfo<T>(curItem.curPos,curItem.curAA+1,0,prunedRot[curItem.curPos][curItem.curAA+1][0]);
        	hasNextItem = true;
        }
        else if(curItem.curAA < prunedRot[curItem.curPos].length-1){
        	nextItem = new RotInfo<T>(curItem.curPos+1,0,0,prunedRot[curItem.curPos+1][0][0]);
        	hasNextItem = true;
        }*/

        //there are no more
    }


    //This isn't quite correct because the max should be updated if the curAA changes,
    //but this works due to the matrix structure that if the curAA is defined rotamer 0 will be defined.
    private void incrementCtr(int[] ctr, int[] max) {
        ctr[ctr.length-1]++;
        for(int i=ctr.length-1; i>0; i--) {
            if(ctr[i] == max[i]) {
                ctr[i] = 0;
                ctr[i-1]++;
            }
        }
    }




    @Override
    public void remove() {
        // TODO Auto-generated method stub
    }
}

class RotInfo<T> {
    int curPos;
    int curAA;
    int curRot;
    T state;

    public RotInfo(int curPos, int curAA, int curRot, T s) {
        this.curPos = curPos;
        this.curAA = curAA;
        this.curRot = curRot;
        this.state = s;
    }

    public String printCoord() {
        return "("+curPos+","+curAA+","+curRot+")";
    }

}

class Index3 implements Comparable<Index3> {
    int pos;
    int aa;
    int rot;

    public Index3(int curPos, int curAA, int curRot) {
        this.pos = curPos;
        this.aa = curAA;
        this.rot = curRot;
    }

    public Index3(String s) {
        //Trim parens
        String sSub = s.substring(1, s.length()-1);
        StringTokenizer st = new StringTokenizer(sSub, ",");
        pos = new Integer(st.nextToken());
        aa = new Integer(st.nextToken());
        rot = new Integer(st.nextToken());
    }

    public String toString() {
        return "("+pos+","+aa+","+rot+")";
    }

    @Override
    public int compareTo(Index3 o) {
        if(pos < o.pos) {
            return -1;
        } else if(pos == o.pos) { //pos either = or great
            if(aa < o.aa) {
                return -1;
            } else if(aa == o.aa) { //pos either = or great
                if(rot < o.rot) {
                    return -1;
                } else if(rot == o.rot) { //pos either = or great
                    return 0;
                } else
                    return 1;
            } else
                return 1;
        } else
            return 1;


    }


}
