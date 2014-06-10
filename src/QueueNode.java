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
//	QueueNode.java
//
//	Version:           2.1 beta
//
//
//	  authors:
// 	  initials    name                 organization                email
//	---------   -----------------    ------------------------    ----------------------------
//	  ISG		 Ivelin Georgiev	  Duke University			  ivelin.georgiev@duke.edu
//	  KER        Kyle E. Roberts       Duke University         ker17@duke.edu
//    PGC        Pablo Gainza C.       Duke University         pablo.gainza@duke.edu
///////////////////////////////////////////////////////////////////////////////////////////////

/**
* Written by Ivelin Georgiev (2004-2009)
*
*/

import java.io.BufferedOutputStream;
import java.io.FileOutputStream;
import java.io.PrintStream;

/**
 * Handles the data for a single node in the A* queue.
 *
 */
public class QueueNode implements Comparable {

    //the f(n) score associated with the current node
    public double fScore;

    //corresponding level
    public int level;

    //the numbers of the nodes in the considered conformation up to the current level
    public int confSoFar[];

    //the number of the corresponding node at that level
    public int nodeNum;

    //a pointer to the previous and next nodes in the expansion list
    public QueueNode prevNode;
    public QueueNode nextNode;

    //constructor
    QueueNode (int curNode, int curLevel, int curConf[], double fn) {

        nodeNum = curNode;
        level = curLevel;

        confSoFar = new int[level+1];

        for (int i=0; i<=level; i++) {
            confSoFar[i] = curConf[i];
        }
        fScore = fn;

        prevNode = null;
        nextNode = null;
    }

    public int compareTo(Object otherObject) throws ClassCastException {
        if(!(otherObject instanceof QueueNode)) {
            throw new ClassCastException("A QueueNode object expected.");
        }
        QueueNode other = (QueueNode)otherObject;
        if(this.fScore > other.fScore) {
            return 1;
        } else if(this.fScore < other.fScore) {
            return -1;
        } else { // Nodes have the same score
            if(this.level > other.level || this.nodeNum != other.nodeNum) {
                return -1; // but different levels, this is larger by default
            } else if(checkConf(this, other)) {
                return 0;
            } else { // Two distinct nodes have the same fScore, say this one is larger.
                return 1;
            }
        }
    }

    //Checks if the two given nodes have the same partially assigned conformation
    private boolean checkConf(QueueNode node1, QueueNode node2) {

        if (node1.level!=node2.level) //different level
            return false;

        for (int l=0; l<node1.confSoFar.length; l++) {
            if (node1.confSoFar[l]!=node2.confSoFar[l])
                return false;
        }

        //The partially assigned conformations are the same
        return true;
    }

}
