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

////////////////////////////////////////////////////////////////////////////////////////////
// RyanQuickSort.java
//
//  Version:           2.1 beta
//
//
// authors:
//    initials    name            organization                email
//   ---------   --------------  ------------------------    ------------------------------
//     RHL        Ryan Lilien     Dartmouth College           ryan.lilien@dartmouth.edu
//
////////////////////////////////////////////////////////////////////////////////////////////

/*
 *
 * This class was written by Ryan Lilien (2001-2004)
 *
 * For some reason Microsoft doesn't implement any sorting functions
 *  so I have to do it myself. This doesn't make any sense, especially
 *  because the standard Java specifications include sorting.
 *
 */

/**
 * This class implements quick sort; sorted objects must implement RyanComparable.
 */
public class RyanQuickSort {

    Object dArray[] = null;

    RyanQuickSort() {}

    // Sorts the array a. Objects in array a must implement
    //  RyanComparable
    public void Sort(Object[] a) {

        RyanQuickSort r = new RyanQuickSort();
        r.dArray = a;
        if (r.isAlreadySorted()) return;
        r.quickSort(0, r.dArray.length-1);
        if (r.isAlreadySorted())
            return;
        else {
            System.out.println("SORT FAILED!");
        }
        return;
    }

    private void quickSort(int q, int w) {
        if (q < w) {
            int p = partition(q,w);
            if (p == w)
                p--;
            quickSort(q,p);
            quickSort(p+1,w);
        }
    }

    private int partition(int b, int t) {
        Object pivot = dArray[b];
        while(true) {
            while ( (((RyanComparable)dArray[t]).compareTo(pivot) >= 0 ) && (b < t) )
                t--;
            while ( (((RyanComparable)dArray[b]).compareTo(pivot) < 0 ) && (b < t) )
                b++;
            if (b < t) {
                // exchange
                Object tmp = dArray[b];
                dArray[b] = dArray[t];
                dArray[t] = tmp;
            } else
                return t;
        }
    }

    // Returns true if array is sorted
    private boolean isAlreadySorted() {
        for(int i=1; i<dArray.length; i++) {
            if (((RyanComparable)dArray[i-1]).compareTo(dArray[i]) > 0)
                return false;
        }
        return true;
    }

}
