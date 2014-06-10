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
// KStar.java
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

/**
 * Written by Ryan Lilien (2001-2004) and Ivelin Georgiev (2004-2009)
 *
 */

import java.io.*;
import java.io.PrintStream;
import java.util.*;

/**
 * This is the main class for the KStar program; essentially just a wrapper for the KSParser class.
 */
public class KStar {

    public static void main (String[] args) {
        KSParser parser = new KSParser();
        parser.checkMPI(args);
        System.out.println("KStar Finished.");

    }  // End main()

} // end class
