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
//	ParamSet.java
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

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.util.StringTokenizer;
import java.io.Serializable;

/**
 * Handles reading in and managing parameter/value pairs from the input configuration files
 */
public class ParamSet implements Serializable {
	
	private String params[] = null; //the parameter names
	private String values[] = null; //the parameter values
	private int curNum = 0; //the number of used positions in the arrays
	
	//constructor
	ParamSet(){
		int initSize = 50; //the initial size of the arrays; increased if necessary
		params = new String[initSize];
		values = new String[initSize];
	}
	
	//Reads in all parameter pairs from the file fName and updates the params[] and values[] arrays
	public void addParamsFromFile(String fName){
		
		BufferedReader bufread = null;
		String curLine = null;
		boolean done = false;
		
		// First attempt to open and read the config file
		try{
			FileInputStream is = new FileInputStream(fName);
			bufread = new BufferedReader(new InputStreamReader(is));

			curLine = bufread.readLine();

			while (curLine != null){
				done = false;
				while (!done) {
					if (curLine.charAt(0) == '%'){
						curLine = bufread.readLine();
					}
					else {
						done = true;
					}
					if (curLine == null)
						done = true;
				}
				if (curLine != null){
					String paramName = getToken(curLine,1);
					int ind = getParamInd(paramName);
					if (ind>=0){ //parameter already read, but only one definition should exist
						System.out.println("ERROR: parameter "+paramName+" already read");
						System.exit(1);
					}
					else { // new parameter
						if (curNum>=params.length){ //double the size of the arrays
							
							String tmp[] = new String[2*params.length];
							System.arraycopy(params, 0, tmp, 0, params.length);
							params = tmp;
							
							tmp = new String[2*values.length];
							System.arraycopy(values,0,tmp,0,values.length);
							values = tmp;
						}
						
						params[curNum] = getToken(curLine,1);
						values[curNum] = curLine.substring(params[curNum].length()+1);
						curNum++;
						curLine = bufread.readLine();
					}
				}
			}
			bufread.close();
		}
		catch(Exception ex)
		{
			System.out.println("ERROR: An error occurred reading configuration file "+fName);
			System.exit(1);
		}
	}
	
	//Returns the index into params[] of the parameter paramName; return -1 if not found
	private int getParamInd(String paramName){
		for (int i=0; i<curNum; i++){
			if (params[i].equalsIgnoreCase(paramName))
				return i;
		}
		return -1;
	}
	
	//Returns the value from values[] that corresponds to the parameter paramName from params[];
	public String getValue(String paramName){
		int ind = getParamInd(paramName);
		if (ind>=0)
			return values[ind];
		else {		
			System.out.println("ERROR: Parameter "+paramName+" not found");
			System.exit(1);
			return null;
		}
	}
	
	//Returns the value from values[] that corresponds to the parameter paramName from params[];
	public String getValue(String paramName,String defaultVal){
		int ind = getParamInd(paramName);
		int processRank = 0;
		try{
			processRank = MPItoThread.Rank();
		}
		catch(Exception e){
			e.printStackTrace();
			System.exit(1);			
		}
		
		if (ind>=0){
			if(processRank == 0){
				System.out.println("Parameter "+paramName+" set to "+values[ind]);
			}
			return values[ind];
		}
		else {
			if(processRank == 0){
				System.out.println("Parameter "+paramName+" not set. Using default value "+defaultVal);
			}
			return defaultVal;
		}
	}
	
	//Sets the value of parameter paramName to newValue
	public void setValue(String paramName, String newValue){
		int ind = getParamInd(paramName);
		if (ind>=0)
			values[ind] = newValue;
		else {
			if (curNum>=params.length){ //double the size of the arrays
				
				String tmp[] = new String[2*params.length];
				System.arraycopy(params, 0, tmp, 0, params.length);
				params = tmp;
				
				tmp = new String[2*values.length];
				System.arraycopy(values,0,tmp,0,values.length);
				values = tmp;
			}
			
			params[curNum] = paramName;
			values[curNum] = newValue;
			curNum++;
			
			//System.out.println("ERROR: Parameter "+paramName+" not found");
			//System.exit(1);
		}
	}
	
	public String [] getParams(){
		return params;
	}
	
	public String [] getValues(){
		return values;
	}
	
	public int getCurNum(){
		return curNum;
	}
	
	public void setParamsValues(String inP[], String inV[], int nP){
		curNum = nP;
		params = new String[curNum];
		values = new String[curNum];
		for (int i=0; i<curNum; i++){
			params[i] = inP[i];
			values[i] = inV[i];
		}
	}
	
	// This function returns the xth token in string s
	private String getToken(String s, int x) {
	
		int curNum = 1;	
		StringTokenizer st = new StringTokenizer(s," ,;\t\n\r\f");
		
		while (curNum < x) {
			curNum++;
			if (st.hasMoreTokens())
			  st.nextToken();
			else {
				// System.out.println("ERROR: Unable to access parameter " + x);
				return(new String(""));
			}
		}
		
		if (st.hasMoreTokens())		
			return(st.nextToken());
		return(new String(""));

	} // end getToken
}
