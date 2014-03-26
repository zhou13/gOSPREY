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
// OneMutation.java
//
//  Version:           2.1 beta
//
//
// authors:
//    initials    name            organization                email
//   ---------   --------------  ------------------------    ------------------------------
//     RHL        Ryan Lilien     Dartmouth College           ryan.lilien@dartmouth.edu
//	  KER        Kyle E. Roberts       Duke University         ker17@duke.edu
//    PGC        Pablo Gainza C.       Duke University         pablo.gainza@duke.edu
////////////////////////////////////////////////////////////////////////////////////////////

/* 
 * Written by Ryan Lilien (2002-2004) and Ivelin Georgiev (2004-2009)
 * 
 */

import java.math.BigDecimal;
import java.math.BigInteger;

/**
 * Handles the data for a single mutation sequence. Contains the amino acid identities for the given sequence
 * and can contain the computed score. Implements a method for comparing two sequences that is used for sorting
 * all sequences with respect to different criteria.
 */
public class OneMutation implements RyanComparable, Comparable
{

	int mutNum = -1;
	BigDecimal score = new BigDecimal("0.0");
	float vol = 0.0f;
	String resTypes[] = null;
	int resMut[] = null;
	Index3 index[] = null;
	String flagMutType = null;
	BigInteger numConfUB = null; //num conformations for the unbound sequence
	BigInteger numConfB = null; //num conformations for the bound sequence
	
	private boolean sortScores = false; //if true, then sort by score (default is false; should only be set using the method below)
	
	OneMutation() {
	}
	
	/*public int compareTo(Object otherObject) {
		OneMutation mut = (OneMutation)otherObject;
		if (score > mut.score) return -1;
		if (score < mut.score) return 1;
		return 0;
	}*/
	
	public int compareTo(Object otherObject){
		OneMutation mut = (OneMutation)otherObject;
		if (!sortScores) {
			String seq1 = "";
			String seq2 = "";
			if (resTypes!=null){
				for (int i=0; i<resTypes.length; i++){
					seq1 += resTypes[i];
					seq2 += mut.resTypes[i];
				}
			}
			else {
				for (int i=0; i<resMut.length; i++){
					seq1 += resMut[i];
					seq2 += mut.resMut[i];
				}
			}
			return (seq1.compareTo(seq2));
		}
		else {// (sortScores==true)
			return score.compareTo(mut.score);
		}
		
	}
	
	//Determines if entries should be sorted by their scores
	public void setSortScores(boolean ss){
		sortScores = ss;
	}
	
	// Returns true if the passed mutation sequence and this
	//  mutation sequence are the same. Otherwise returns 0.
	public boolean isSame(String anotherMutation[]) {
		for(int i=0;i<anotherMutation.length;i++) {
			if (!anotherMutation[i].equalsIgnoreCase(resTypes[i]))
				return(false);
		}
		return(true);
	}	

}
