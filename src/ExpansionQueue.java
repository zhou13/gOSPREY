import java.util.concurrent.PriorityBlockingQueue;

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
//	ExpansionQueue.java
//
//	Version:           2.1 beta
//
//
//	  authors:
// 	  initials    name                 organization                email
//	 ---------   -----------------    ------------------------    ----------------------------
//	  ISG		 Ivelin Georgiev	  Duke University			  ivelin.georgiev@duke.edu
//	  KER        Kyle E. Roberts       Duke University         ker17@duke.edu
//    PGC        Pablo Gainza C.       Duke University         pablo.gainza@duke.edu
///////////////////////////////////////////////////////////////////////////////////////////////

/**
* Written by Ivelin Georgiev (2004-2009)
* 
*/

/**
 * This queue is ordered in terms of increasing f(n) values of the nodes in the A* expansion tree;
 * 		only the visible nodes are contained in the queue.
 */
public class ExpansionQueue {
	
	 private PriorityBlockingQueue<QueueNode> thequeue;
	
	//a pointer to the first node in the expansion list
	public QueueNode curFront;
	
	//number of nodes in the list
	public int numNodes;
	
	//the unique id of each node in the queue
	public int idNUM;
	
	//constructor
	ExpansionQueue () {
		curFront = null;
		numNodes = 0;
		
		 thequeue = new PriorityBlockingQueue<QueueNode>(50000);
		
	}
	/*
	//inserts a new node into the expansion queue
	public void insert (QueueNode newNode){
		
		boolean done;
		int i;
		QueueNode curNode = null;
		
		if (numNodes==0){//first node to be inserted
			curFront = newNode;
			numNodes++;
		}
		else {//expansion list is not empty
			curNode = curFront;
			done = false;
			i=0;
			while ((!done)&&(i<numNodes)){
				if (curNode.fScore > newNode.fScore){//insert the new node right before curNode
					
					newNode.nextNode = curNode;
					newNode.prevNode = curNode.prevNode;
					
					if (i!=0)//curNode!=curMinNode
						curNode.prevNode.nextNode = newNode;
					else //we have a new minimum node
						curFront = newNode;
					curNode.prevNode = newNode;					
					
					numNodes++;
					done = true;
				}
				else {
					if (i==numNodes-1){
						newNode.prevNode = curNode;
						curNode.nextNode = newNode;
						numNodes++;
						done = true;
					}
					else {
						curNode = curNode.nextNode;
						i++;
					}
				}
			}
		}
	}

	//Deletes the specified node from the queue: determines which node to delete
	//	based on the fScore, the node level, the node number, and confSoFar[] (only
	//	the last one is sufficient, but it is faster not to check the full length
	//	of confSoFar[] for each node in the queue)
	public void delete (QueueNode delNode){
		
		QueueNode curNode = curFront;
		int i = 0;
		boolean done = false;
		
		while (!done){
			if(delNode.fScore==curNode.fScore){
				if((delNode.level==curNode.level)&&(delNode.nodeNum==curNode.nodeNum)){
					if (checkConf(delNode,curNode)){ //found the node to delete
						if(i!=numNodes-1)
							curNode.nextNode.prevNode = curNode.prevNode;
						if(i!=0)
							curNode.prevNode.nextNode = curNode.nextNode;
						if(i==0)//move the front pointer if we are deleting the min node
							curFront = curNode.nextNode;
						
						numNodes--;
						done = true;
						break;
					}
				}
			}
			if (!done){ //go to the next node in the queue
				curNode = curNode.nextNode;
				i++;
			}
		}
	}*/
	
	 public void insert(QueueNode newNode){
         thequeue.add(newNode);
         curFront = thequeue.peek();
	 }
	 public void delete(QueueNode delNode){
         thequeue.remove(delNode);
         curFront = thequeue.peek();
	 }
	 public QueueNode getMin(){
         return thequeue.poll();
	}
	
	//Checks if the two given nodes have the same partially assigned conformation
	private boolean checkConf(QueueNode node1, QueueNode node2){
	
		if (node1.level!=node2.level) //different level
			return false;
		
		for (int l=0; l<node1.confSoFar.length; l++){
			if (node1.confSoFar[l]!=node2.confSoFar[l])
				return false;
		}
		
		//The partially assigned conformations are the same
		return true;		
	}
}
