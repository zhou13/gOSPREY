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
//	EnvironmentVars.java
//
//	Version:           2.1 beta
//
//
//	  authors:
// 	  initials    name                 organization                email
//	 ---------   -----------------    ------------------------    ----------------------------
//     KER        Kyle E. Roberts       Duke University         ker17@duke.edu
///////////////////////////////////////////////////////////////////////////////////////////////


import java.math.BigInteger;


public class EnvironmentVars {
	
	public enum FORCEFIELD {
		AMBER, CHARMM22, CHARMM19NEUTRAL, CHARMM19
	}
	
	static String dataDir = "./";
	static FORCEFIELD forcefld; 

	//static BigInteger maxKSconfs;
	//static boolean useMaxKSconfs = false;
	
	public static boolean useEntropy = false;
	private static double[] entropy = null;
	public static double entropyScaling = 0;
	
	public static String ksConfDir = "ksConfs";
	
	public static RotamerLibrary aaRotLib;

        public static boolean autoFix;//Should structures being read in be autofixed?
	
	public static String getDataDir() {
		return dataDir;
	}

	public static void setDataDir(String dd) {
		if(!dd.endsWith("/") || !dd.endsWith("\\")){
			dd = dd.concat("/");
		}
		EnvironmentVars.dataDir = dd;
	}
	
	public static void setForcefld(String frcefld) {
		forcefld = FORCEFIELD.valueOf(frcefld.toUpperCase());
	}
	
	/*public static void setMaxKSconfs(BigInteger confs){
		maxKSconfs = confs;
	}
	
	public static void setUseMaxKSConfs(boolean val){
		useMaxKSconfs = val;
	}*/
	
	public static double[] getEntropyTerms() {
		if(entropy == null){
			double[] entr = {0,0.5,0.75,0.75,0.58,0.99,0.97,1.14,1.53,1.19,1.12,2.21,2.13,0.99,0.99,0.99,0.61,1.65,0.81,2.02,0};
			entropy = entr;	
			for(int i=0;i<entropy.length;i++)
				entropy[i] *= entropyScaling;
		}
		
		return entropy;
	}
	
	public static void setEntropyScale(double es){
		if(es>0)
			useEntropy = true;
		entropyScaling = es;
	}
	
	public static void setAArotLib(String rl){
		aaRotLib = new RotamerLibrary(rl, true);
	}


}
