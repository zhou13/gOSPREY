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
//	PertFileHandler.java
//
//	Version:           2.1 beta
//
//
//	  authors:
// 	  initials    name                 organization                email
//	 ---------   -----------------    ------------------------    ----------------------------
//	  MAH           Mark A. Hallen	  Duke University               mah43@duke.edu
///////////////////////////////////////////////////////////////////////////////////////////////
import java.io.*;
import java.util.*;


public class PertFileHandler {

    public static void readPertFile(String fileName, Molecule m, StrandRotamers strandRot[]){
        //Read the perturbation file and sets up the perturbations and residue conformations indicated in it

       /* int sysStrNum = sysLR.strandNumber;
        int ligStrNum = -1;
        int ligAANum = -1;

        if( ligROT != null ){
            ligStrNum = ligROT.strandNumber;
            ligAANum = ligROT.getIndexOfNthAllowable(0, 0);
        }*/

        int numStrands = m.numberOfStrands;//This should also be strandRot.length


        //Perturbations must be listed in the order they'll be applied.
        //They should be in as few non-overlapping layers as possible
        //with larger ones coming first, in order to minimize extra RCs from indirect effects
        try{

            BufferedReader br=new BufferedReader(new FileReader(fileName));
            StringTokenizer st;
            br.readLine();//Title
            m.perts = readPerts( br, m );

            /*
            int numResToRead = 0;//Number of residues whose information needs to be read (residues affected by perturbations)
            int numResRead = 0;//Number of residues whose information has been read

            
             This code built Residue.perts only from perturbations directly affecting each residue:
             *
             *
             for(int a=0;a<m.residue.length;a++)//Creating the Residue.perts arrays
                m.residue[a].perts = new int[resNumPerts[a]];
            int pertInd[] = new int[m.residue.length];
            for(int a=0;a<numPerts;a++){
                Perturbation pert = m.perts[a];
                for(int b=0;b<pert.resAffected.length;b++){
                    int rn = pert.resAffected[b];
                    m.residue[rn].perts[pertInd[rn]] = a;
                    pertInd[rn]++;
                    if(pertInd[rn]==1)
                        numResToRead++;
                }
             *
             *
            }*/



            //int numAATypes = sysLR.getNumAATypes();//How many AA types there are
            int RCRots[][][][] = new int[numStrands][][][];
            int RCPertStates[][][][] = new int[numStrands][][][];
            //int RCRotsLig[][][] = null, RCPertStatesLig[][][] = null, resRCCount[] = null, RCInfo[][] = null;

            for(int str=0; str<numStrands; str++){
                RCRots[str] = new int[m.strand[str].numberOfResidues][((StrandRCs)strandRot[str]).getNumAATypes()][];
                RCPertStates[str] = new int[m.strand[str].numberOfResidues][((StrandRCs)strandRot[str]).getNumAATypes()][];
            }

            while(br.readLine() != null){//Read residue perturbation states and RCs.  Skipping the "RES" line

                String inputNumber = br.readLine().trim();
                int molResNum = m.mapPDBresNumToMolResNum(inputNumber);

                if(molResNum != -1){//Don't read perturbation information for residues that aren't in the system
                    //(for example if we're only considering one member of a complex and the perturbation file is for the whole complex)

                    Residue res = m.residue[molResNum];
                    int strResNum = res.strandResidueNumber;
                    int strandNumber = res.strandNumber;
                    StrandRCs curStrandRCs = (StrandRCs)strandRot[strandNumber];
                    int numAATypes = curStrandRCs.getNumAATypes();

                    st = new StringTokenizer(br.readLine()," ");
                    int numStates = Integer.valueOf(st.nextToken());
                    st.nextToken();
                    int numRCs = Integer.valueOf(st.nextToken());

                    //Read the perturbations affecting this residue
                    //This section is included in the file because it specifies if some indirect effects should be neglected
                    st = new StringTokenizer(br.readLine()," ");
                    res.perts = new int[st.countTokens() - 1];
                    st.nextToken();//"PERTURBATIONS"
                    for(int a=0;a<res.perts.length;a++)
                        res.perts[a] = Integer.valueOf(st.nextToken());

                    res.pertStates = new int[numStates][res.perts.length];

                    for(int a=0;a<numStates;a++){//Read perturbation states
                        st = new StringTokenizer(br.readLine()," ");
                        for(int b=0;b<res.perts.length;b++)
                            res.pertStates[a][b] = Integer.valueOf(st.nextToken());
                    }

                    if( !br.readLine().equalsIgnoreCase("RCs") )
                        throw new Exception("Bad formatting of residue "+inputNumber);

                    int[] resRCCount = new int[numAATypes];
                    int[][] RCInfo = new int[numRCs][3];//Each row has the AA type, rotamer, and perturbation state for an RC


                    int numRCsRead = 0;

                    for(int a=0;a<numRCs;a++){//Read RCs
                        st = new StringTokenizer(br.readLine()," ");

                        String AAName = st.nextToken();
                        RCInfo[a][0] = curStrandRCs.rl.getAARotamerIndex(AAName);
                        RCInfo[a][1] = Integer.valueOf(st.nextToken());
                        RCInfo[a][2] = Integer.valueOf(st.nextToken());
                        resRCCount[RCInfo[a][0]]++;
                        numRCsRead++;
                    }

                    if(numRCsRead != numRCs)
                        throw new Exception("Wrong number of RCs provided for residue "+inputNumber);


                    //Separate RCs according to their AA type
                    for(int curAA=0; curAA<numAATypes; curAA++){
                        RCRots[strandNumber][strResNum][curAA] = new int[resRCCount[curAA]];
                        RCPertStates[strandNumber][strResNum][curAA] = new int[resRCCount[curAA]];
                    }

                    int rcind[] = new int[numAATypes];
                    for(int a=0;a<numRCs;a++){
                        RCRots[strandNumber][strResNum][RCInfo[a][0]][rcind[RCInfo[a][0]]] = RCInfo[a][1];
                        RCPertStates[strandNumber][strResNum][RCInfo[a][0]][rcind[RCInfo[a][0]]] = RCInfo[a][2];
                        rcind[RCInfo[a][0]]++;
                    }

                    //numResRead++;
                }

            }

            //if(numResToRead != numResRead)
            //    throw new Exception("Information given for "+numResRead+" residues; should be "+numResToRead);

            for(int str=0; str<numStrands; str++){
                ((StrandRCs)strandRot[str]).RCRots = RCRots[str];
                ((StrandRCs)strandRot[str]).RCPertStates = RCPertStates[str];
            }

            //Find the successors and, from those, the predecessors to each perturbation
            m.findPerturbationSuccessors();

            for(Residue res : m.residue )//Assign each residue's affected perturbations
                m.assignAffectedPerts(res);
            
            br.close();
        }

        catch(Exception e){
            System.out.println("Error reading perturbation file:");
            System.out.println(e.getMessage());
            e.printStackTrace();
            System.exit(1);//Problems setting up the perturbations will mess up the whole run so quit
        }
    }


    
    public static Perturbation[] readPerts(BufferedReader br, Molecule m) throws Exception {

            int numPerts = Integer.valueOf(br.readLine().trim());
            Perturbation[] perts = new Perturbation[numPerts];
            //int resNumPerts[] = new int[m.residue.length];//Number of perturbations for each residue

            StringTokenizer st;

            for(int a=0;a<numPerts;a++){//Read perturbations
                String pertType = br.readLine();
                st = new StringTokenizer(br.readLine()," ");
                int numAffectedRes = st.countTokens();
                int affectedRes[] = new int[numAffectedRes];

                for(int b=0;b<numAffectedRes;b++){
                    String inputNumber = st.nextToken();
                    affectedRes[b] = m.mapPDBresNumToMolResNum(inputNumber);
                    //resNumPerts[affectedRes[b]]++;
                }

                perts[a] = Perturbation.generate(pertType, m, affectedRes, br);
                st = new StringTokenizer(br.readLine()," ");
                int numStates = Integer.valueOf(st.nextToken());
                perts[a].maxParams = new float[numStates];
                perts[a].minParams = new float[numStates];

                for(int b=0;b<numStates;b++){
                    st = new StringTokenizer(br.readLine()," ");
                    if(st.countTokens() != 2)
                        throw new java.lang.Exception("Bad formatting of perturbation "+a);
                    perts[a].minParams[b] = Float.valueOf(st.nextToken());
                    perts[a].maxParams[b] = Float.valueOf(st.nextToken());
                }
            }

            return perts;
    }

    

    public static void writePertFile(String fileName, Molecule m, PrunedRotamers<Boolean> eliminatedRCAtPos,
            StrandRotamers strandRot[], int strandMut[][], boolean screenOutput ){
        //Generates a perturbation file based on the perturbations, states, and RCs in m, sysLR, and ligROT
        //If screenOutput is true, then writes a "P" after every pruned
        //perturbation, perturbation state, residue perturbation state, or RC

        //strandMut is only needed for screenOutput (can be null otherwise)

        //Arrays indicating whether each of these things are pruned
        boolean prunedPerts[] = null;
        boolean prunedPertStates[][] = null;
        boolean prunedResPertStates[][] = null;
        boolean prunedRCs[][][] = null;

        StrandRCs curStrandRCs;

        //int ligNum = -1;
        //if( ligROT != null)
        //    ligNum = m.strand[ligROT.strandNumber].residue[0].moleculeResidueNumber;
        
        
        if(screenOutput){//Check what's pruned

            prunedPerts = new boolean[m.perts.length];
            prunedPertStates = new boolean[m.perts.length][];
            prunedResPertStates = new boolean[m.residue.length][];
            prunedRCs = new boolean[m.residue.length][][];


            for(int pertNum=0;pertNum<m.perts.length;pertNum++){//Initialize the non-RC ones to true
                prunedPerts[pertNum] = true;
                prunedPertStates[pertNum] = new boolean[m.perts[pertNum].maxParams.length];
                for(int state=0;state<prunedPertStates[pertNum].length;state++)
                    prunedPertStates[pertNum][state] = true;
            }
            for(int pos=0; pos<m.residue.length; pos++){
                Residue res = m.residue[pos];

                if(res.perts != null){
                    curStrandRCs = (StrandRCs)strandRot[res.strandNumber];
                    prunedResPertStates[pos] = new boolean[res.pertStates.length];
                    prunedRCs[pos] = new boolean[curStrandRCs.getNumAATypes()][];
                    for(int state=0;state<prunedResPertStates[pos].length;state++)
                        prunedResPertStates[pos][state] = true;
                }
            }

            int flexPos = 0;//position among flexible residues
            int flexPosInStrand = 0;//Position among flexible residues in the strand

            for(int pos=0; pos<m.residue.length; pos++){//pos is a molecule residue number
                Residue res = m.residue[pos];
                int posInStrand = res.strandResidueNumber;

                if(posInStrand==0)
                    flexPosInStrand=0;

                if(res.perts != null){
                    curStrandRCs = (StrandRCs)strandRot[res.strandNumber];

                    for(int AA=0;AA<curStrandRCs.getNumAllowable(posInStrand);AA++){
                        int curAA = curStrandRCs.getIndexOfNthAllowable(posInStrand,AA);

                        int numRCs = curStrandRCs.getNumRCs(posInStrand, curAA);
                        prunedRCs[pos][curAA] = new boolean[numRCs];

                        for(int RC=0; RC<numRCs; RC++){
                            boolean isPruned = eliminatedRCAtPos.get(flexPos, curAA, RC);
                            prunedRCs[pos][curAA][RC] = isPruned;
                            if(!isPruned){
                                int resPertState = curStrandRCs.RCPertStates[posInStrand][curAA][RC];
                                prunedResPertStates[pos][resPertState] = false;
                                for(int a=0;a<res.pertStates[resPertState].length;a++){//Loop over perturbations of this residue
                                    int pertNum = res.perts[a];
                                    int pertState = res.pertStates[resPertState][a];
                                    prunedPertStates[pertNum][pertState] = false;
                                    if(pertState > 0)//If only unperturbed states of a perturbation are left unpruned, we can call that perturbation pruned
                                        prunedPerts[pertNum] = false;
                                }
                            }
                        }
                    }
                }

                if( strandMut[res.strandNumber][flexPosInStrand] == posInStrand){
                    flexPos++;
                    flexPosInStrand++;
                }

            }
        }
        


        try{
            BufferedWriter bw=new BufferedWriter(new FileWriter(fileName));
            bw.append("PERTURBATIONS");
            bw.newLine();
            bw.append(String.valueOf(m.perts.length));
            bw.newLine();

            for(int pertNum=0;pertNum<m.perts.length;pertNum++){//Perturbation info
                Perturbation pert = m.perts[pertNum];
                bw.append(pert.type);

                if(screenOutput){
                    if(prunedPerts[pertNum])
                        bw.append(" P");
                }
                bw.newLine();

                for(int rc=0;rc<pert.resAffected.length;rc++){
                    int molResNum=pert.resAffected[rc];
                    bw.append(m.residue[molResNum].getResNumber()+" ");
                }

                bw.newLine();

                pert.writeExtraInfo(bw);

                bw.append(pert.maxParams.length+" states");
                bw.newLine();
                for(int state=0;state<pert.maxParams.length;state++){
                    bw.append(pert.minParams[state]+" "+pert.maxParams[state]);
                    if(screenOutput){
                        if(prunedPertStates[pertNum][state])
                            bw.append(" P");
                    }
                    bw.newLine();
                }
            }

            for(int pos=0;pos<m.residue.length;pos++){//Residue perturbation state, RC info
                Residue res=m.residue[pos];
                int posInStrand = res.strandResidueNumber;

                if( res.perts.length > 0 ){
                    curStrandRCs = (StrandRCs)strandRot[res.strandNumber];

                    bw.append("RES");
                    bw.newLine();
                    bw.append(String.valueOf(res.getResNumber()));
                    bw.newLine();

                    if(screenOutput)
                        bw.append(res.pertStates.length + " states " +
                            curStrandRCs.getTotNumRCs(posInStrand) + " RCs");
                    else//We will not print unperturbed RCs in this case
                        bw.append(res.pertStates.length + " states " +
                            curStrandRCs.getTotNumPerturbedRCs(posInStrand) + " RCs");
                    bw.newLine();

                    bw.append("PERTURBATIONS ");
                    for(int a=0;a<res.perts.length;a++)
                        bw.append(res.perts[a] + " ");
                    bw.newLine();

                    for(int state=0;state<res.pertStates.length;state++){
                        for(int pertInd=0;pertInd<res.pertStates[state].length;pertInd++)
                            bw.append(res.pertStates[state][pertInd] + " ");

                        if(screenOutput){
                            if(prunedResPertStates[pos][state])
                                bw.append("P");
                        }
                        bw.newLine();
                    }

                    bw.append("RCs");
                    bw.newLine();

                    for(int AA=0;AA<curStrandRCs.getNumAllowable(posInStrand);AA++){
                        int curAA = curStrandRCs.getIndexOfNthAllowable(posInStrand,AA);
                        String AAName = curStrandRCs.rl.getAAName(curAA);

                        for(int curRC=0; curRC<curStrandRCs.RCRots[posInStrand][curAA].length; curRC++){

                            if( ( curStrandRCs.RCPertStates[posInStrand][curAA][curRC] != 0 ) && ( ! screenOutput ) ){
                                //If the perturbation file is to be loaded later don't include unperturbed RCs...StrandRCs.addUnperturbedRCs will add those
                                bw.append( AAName + " " + curStrandRCs.RCRots[posInStrand][curAA][curRC] + " "
                                    + curStrandRCs.RCPertStates[posInStrand][curAA][curRC] );
                                bw.newLine();
                            }
                            
                            else if(screenOutput)
                            {//Print the information for all RCs if this is screening output
                                bw.append( AAName + " " + curStrandRCs.RCRots[posInStrand][curAA][curRC] + " "
                                    + curStrandRCs.RCPertStates[posInStrand][curAA][curRC] );
                                
                                if(prunedRCs[pos][curAA][curRC])
                                    bw.append(" P");

                                bw.newLine();
                            }
                        }
                    }
                }
            }

            bw.close();
        }

        catch(IOException e){
            System.out.println(e.getMessage());
            e.printStackTrace();
        }

    }

}
