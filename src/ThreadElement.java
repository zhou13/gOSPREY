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
//	ThreadElement.java
//
//	Version:           2.1 beta
//
//
//	  authors:
// 	  initials    name                 organization                email
//	 ---------   -----------------    ------------------------    ----------------------------
//	  KER        Kyle E. Roberts       Duke University         ker17@duke.edu
///////////////////////////////////////////////////////////////////////////////////////////////

import java.util.concurrent.BlockingDeque;
import java.util.concurrent.LinkedBlockingDeque;
import java.lang.Runtime;
import mpi.Datatype;
import mpi.Status;

public class ThreadElement {
	/***** Thread Variables *****/
	BlockingDeque<ThreadMessage> newMessages = null;
	BlockingDeque<ThreadMessage> oldMessages = null;
	//KSParser[] kstarArray = null;
	int rank = -1;
	
	public ThreadElement(int rank) {
		// TODO Auto-generated constructor stub
		this.rank = rank;
		//kstarArray = ksa;
		newMessages = new LinkedBlockingDeque<ThreadMessage>();
		oldMessages = new LinkedBlockingDeque<ThreadMessage>();
	}
	
	public int getRank(){
		return rank;
	}
	
	public int getSize(){
		return Runtime.getRuntime().availableProcessors() + 1;
	}
	
	public ThreadStatus Probe(int source, int tag) throws InterruptedException{
		ThreadStatus s = null;
		while(s == null){
			ThreadMessage t = newMessages.takeFirst();
			if((source == -1 || t.source == source) && (tag==-1 || t.tag == tag)){
				s = new ThreadStatus();
				s.source = t.source;
				s.tag    = t.tag;
				oldMessages.addFirst(t);
			}
			else{
				oldMessages.addFirst(t);
			}
		}
		
		resetMessages();
		
		return s;
	}
	
	public ThreadStatus Iprobe(int source, int tag) throws InterruptedException{
		ThreadStatus s = null;
		for(ThreadMessage t : newMessages){
			if((source == -1 || t.source == source) && (tag==-1 || t.tag == tag)){
				s = new ThreadStatus();
				s.source = t.source;
				s.tag    = t.tag;
			}
		}
		
		return s;
	}
	
	
	public ThreadStatus Recv(Object buf, int offset, int count, int type,
			int source, int tag) throws InterruptedException{
		ThreadStatus s = null;
		while(s == null){
			ThreadMessage t = newMessages.takeFirst();
			if((source == -1 || t.source == source) && (tag==-1 || t.tag == tag)){
				s = new ThreadStatus();
				s.source = t.source;
				s.tag    = t.tag;
				Object[] objs = (Object[])t.getObj(offset, count);
				//Need to special case int since it isn't an object
				if(type == ThreadMessage.INT){
					for(int i=0; i<count;i++)
						((int[])buf)[i] = (Integer)objs[i];
				}
				else if(type == ThreadMessage.BOOLEAN){
					for(int i=0; i<count;i++)
						((boolean[])buf)[i] = (Boolean)objs[i];
				}
				else if(type == ThreadMessage.FLOAT){
					for(int i=0; i<count;i++)
						((float[])buf)[i] = (Float)objs[i];
				}
				else{
					for(int i=0; i<count;i++)
						((Object[])buf)[i] = objs[i];
				}
				//((CommucObj[])buf)[0].arpFilenameMin = ((CommucObj[])(t.getObj(offset,count)))[0].arpFilenameMin;
			}
			else{
				oldMessages.addFirst(t);
			}
		}
		
		resetMessages();
		
		return s;
	}
	
	public void Send(Object buf, int offset, int count, int type,
			int dest, int tag) {
		
		ThreadMessage t = new ThreadMessage(buf, type, rank, tag);
		t.obj = t.determineObj(buf, type, offset, count);
		try {
			t.obj = MPItoThread.deepCopy(t.obj);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		(MPItoThread.getThreadElement(dest)).setMessage(t);
		
	}
	
	public void resetMessages(){
		while(!oldMessages.isEmpty()){
			ThreadMessage t = null;
			try {
				t = oldMessages.takeFirst();
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			newMessages.addFirst(t);
		}
	}
	
	public void setMessage(ThreadMessage t){
		newMessages.addLast(t);
	}
	
	public class ThreadStatus {
		public int source;
		public int tag;
		
		public ThreadStatus(){
			
		}

	}

}
