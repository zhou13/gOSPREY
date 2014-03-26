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
//	MPItoThread.java
//
//	Version:           2.1 beta
//
//
//	  authors:
// 	  initials    name                 organization                email
//	 ---------   -----------------    ------------------------    ----------------------------
//	  KER        Kyle E. Roberts       Duke University         ker17@duke.edu
///////////////////////////////////////////////////////////////////////////////////////////////


import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.Hashtable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import mpi.Datatype;
import mpi.MPI;
import mpi.MPIException;
import mpi.Status;


public class MPItoThread {
	static Hashtable<Thread,ThreadElement> threadEle = null;
	static ThreadElement[] threadEleArray = null;
	static ExecutorService exe = null;
	//static Hashtable<Thread,KSParser> parsers = new Hashtable<Thread,KSParser>();
	static boolean mpiRun = false;
	static int numThreads;
	static int numProc;
	

	
	
	public MPItoThread(boolean mpiRun, int numThreads) {
		super();
		MPItoThread.mpiRun = mpiRun;
		MPItoThread.numThreads = numThreads;
		MPItoThread.numProc = numThreads+1;
	}
	
	public static void initialize(boolean mpiRun, int numThreads) {
		MPItoThread.mpiRun = mpiRun;
		MPItoThread.numThreads = numThreads;
		MPItoThread.numProc = numThreads+1;
		
		if(!mpiRun){
			threadEleArray = new ThreadElement[numProc];
			threadEle = new Hashtable<Thread,ThreadElement>();
		}
	}
	
	public static ThreadElement getThreadElement(int i){
		return threadEleArray[i];
	}

	public static void startThreads(KSParser mainKSP,Thread mainThread){
		if(exe == null){
			//numThreads =  1; //Runtime.getRuntime().availableProcessors();
			//numProc = numThreads+1;
			mainKSP.numProc = numProc;
			exe = Executors.newFixedThreadPool(numThreads);
			
			KSthread[] kst = new KSthread[numThreads];
			//KSParser[] ksp = new KSParser[numThreads+1];
			//ksp[0] = mainKSP;
			MPItoThread.threadEleArray[0] = threadEle.get(mainThread);
			for(int i=0; i<numThreads; i++){
				threadEleArray[i+1] = new ThreadElement(i+1);
				kst[i] = new KSthread(new KSParser(),i+1,threadEleArray[i+1]);
				kst[i].ksp.cfgName = mainKSP.cfgName;
				//ksp[i+1] = kst[i].ksp;
			}
			
			//threadEle.kstarArray = ksp;
			for(int i=0; i<numThreads;i++){
				//kst[i].ksp.threadEle.kstarArray = ksp;
				exe.execute(kst[i]);
			}
		}
	}
	
	/******** Threaded Functions *********/
	public static int Rank() throws MPIException{
		if(mpiRun)
			return MPI.COMM_WORLD.Rank();
		else
			return threadEle.get(Thread.currentThread()).getRank();
	}
	
	public static int Size() throws MPIException{
		if(mpiRun)
			return MPI.COMM_WORLD.Size();
		else
			return threadEle.get(Thread.currentThread()).getSize();
	}
	
	public static Object Probe(int source, int tag) throws MPIException, InterruptedException{
		if(mpiRun){
			if(source == -1)
				source = MPI.ANY_SOURCE;
			if(tag == -1)
				tag = MPI.ANY_TAG;
			return MPI.COMM_WORLD.Probe(source, tag);
		}
		else
			return threadEle.get(Thread.currentThread()).Probe(source, tag);
		
	}
	
	public static Object Iprobe(int source, int tag) throws MPIException, InterruptedException{
		if(mpiRun){
			if(source == -1)
				source = MPI.ANY_SOURCE;
			if(tag == -1)
				tag = MPI.ANY_TAG;
			return MPI.COMM_WORLD.Iprobe(source, tag);
		}
		else
			return threadEle.get(Thread.currentThread()).Iprobe(source, tag);
		
	}
	
	//Either going to return a status or ThreadStatus
	public static Object Recv(Object buf, int offset, int count, int type,
			int source, int tag) throws MPIException, InterruptedException{
		if(mpiRun){
			if(source == -1)
				source = MPI.ANY_SOURCE;
			if(tag == -1)
				tag = MPI.ANY_TAG;
			Datatype dtype = getDatatype(type);
			return MPI.COMM_WORLD.Recv(buf, offset, count, dtype, source, tag);
		}
		else
			return threadEle.get(Thread.currentThread()).Recv(buf, offset, count, type, source, tag);
		
	}
	
	public static void Send(Object buf, int offset, int count, int type,
			int dest, int tag) throws MPIException, InterruptedException{
		if(mpiRun){
			Datatype dtype = getDatatype(type);
			MPI.COMM_WORLD.Send(buf, offset, count, dtype, dest, tag);
		}
		else
			threadEle.get(Thread.currentThread()).Send(buf, offset, count, type, dest, tag);
		
	}
	
	//ISEND not supported yet for threaded functions
	//will just call the normal send command
	public static void Isend(Object buf, int offset, int count, int type,
			int dest, int tag) throws MPIException, InterruptedException{
		if(mpiRun){
			Datatype dtype = getDatatype(type);
			MPI.COMM_WORLD.Isend(buf, offset, count, dtype, dest, tag);
		}
		else
			threadEle.get(Thread.currentThread()).Send(buf, offset, count, type, dest, tag);
		
	}
	
	public static Datatype getDatatype(int type){
		switch(type){
			case ThreadMessage.BOOLEAN :
				return MPI.BOOLEAN;
			case ThreadMessage.INT :
				return MPI.INT;
			case ThreadMessage.OBJECT :
				return MPI.OBJECT;
			case ThreadMessage.DOUBLE :
				return MPI.DOUBLE;
			case ThreadMessage.FLOAT :
				return MPI.FLOAT;
			default:
				System.out.println("MPI OBJECT NOT CODED YET");
				System.exit(0);
		}
		return null;
		
		
	}
	
	public static int getStatusTag(Object obj){
		if(obj instanceof ThreadElement.ThreadStatus){
			return ((ThreadElement.ThreadStatus)obj).tag;
		}
		else if(obj instanceof mpi.Status){
			return ((Status)obj).tag;
		}
		else{
			System.out.println("Status not defined");
			System.exit(0);
		}
		return -1;
	}
	
	public static int getStatusSource(Object obj){
		if(obj instanceof ThreadElement.ThreadStatus){
			return ((ThreadElement.ThreadStatus)obj).source;
		}
		else if(obj instanceof mpi.Status){
			return ((Status)obj).source;
		}
		else{
			System.out.println("Status not defined");
			System.exit(0);
		}
		return -1;
	}
	/********* End Threaded Functions *****************/
	
	//Function taken from: http://www.javaworld.com/javaworld/javatips/jw-javatip76.html?page=2
	//Java Tip 76: An alternative to the deep copy technique
	//Author: Dave Miller
	static public Object deepCopy(Object oldObj) throws Exception {
		ObjectOutputStream oos = null;
	      ObjectInputStream ois = null;
	      try
	      {
	         ByteArrayOutputStream bos = 
	               new ByteArrayOutputStream(); // A
	         oos = new ObjectOutputStream(bos); // B
	         // serialize and pass the object
	         oos.writeObject(oldObj);   // C
	         oos.flush();               // D
	         ByteArrayInputStream bin = 
	               new ByteArrayInputStream(bos.toByteArray()); // E
	         ois = new ObjectInputStream(bin);                  // F
	         // return the new object
	         return ois.readObject(); // G
	      }
	      catch(Exception e)
	      {
	         System.out.println("Exception in ObjectCloner = " + e);
	         throw(e);
	      }
	      finally
	      {
	         oos.close();
	         ois.close();
	      }
	   }
	
}
