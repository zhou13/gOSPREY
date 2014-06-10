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
//	Perturbation.java
//
//	Version:           2.1 beta
//
//
//	  authors:
// 	  initials    name                 organization                email
//	 ---------   -----------------    ------------------------    ----------------------------
//	  MAH           Mark A. Hallen	  Duke University               mah43@duke.edu
///////////////////////////////////////////////////////////////////////////////////////////////
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.util.*;


public abstract class Perturbation {
//Extended by classes like Shear, Backrub
    //Each perturbation has a single parameter of type float
    //However, discrete perturbations may admit as few as two discrete values of this parameter
    //Perturbations can move any number of residues, but must treat the sidechain as a rigid body
    //and should not mess with chi1 (the idea of perturbations is to be orthogonal to the SC-only adjustments)
    //They also should not interfere with the ligand rotations and translations, if any
    //Parameter value 0 and curState=0 both denote the unperturbed state

    String type;//e.g. "Shear"
    Molecule m;//The affected molecule
    int resAffected[];//Molecule-based numbers of affected residues

    float minParams[];//The minimum parameter value for each state of this perturbation
    //(Or the only one if this is a discrete perturbation or a discrete state: then minParams[number] and maxaParams[number] will be the same)
    float maxParams[];//The maximum parameter value for each state of this perturbation
    //The first perturbation state is always the unperturbed state: generally this means minParams[0] == maxParams[0] == 0
    //But it can also mean minParams=-a, maxParams=a for some small a (StrandRCs.applyPertState applies 0 for this state  when trying to save time for invalid perturbation states)
    //These states will be incoporated into RCs

    int curState = 0;//The current state of the perturbation (index in minParams, maxParams)
    //Starts out unperturbed.  This state will be changed by StrandRCs.applyRC and implemented in the actualCoordinates and Atom.coord array

    float curParam = 0;//The current parameter value of the perturbation.  This is changed by applyPerturbation and may be only implemented in the actualCoordinates

    //The current PMinimizer assumes that minimized perturbations commute with sidechain dihedral changes
    //(as all backbone perturbations that treat the sidechain as a rigid body do)


    int[] successors;//Perturbations that need to be applied after this one because they are affected by this one's state
    //These are indices in m.perts

    int[] predecessors;//Similarly, these indices in m.perts are perturbations that need to be applied before this one because they affect its state

    ArrayList<HashMap<String, float[]>> oldCoords;//Old coordinates to revert to upon undoing

    static boolean idealizeSC = true;//Should sidechains be idealized?



//******Applying perturbations********

//These functions apply perturbations to all affected residues, in the m.actualCoordinates

    abstract public boolean doPerturbationMotion(float param);
    //This function should apply the actual perturbation motion to the actualCoordinates
    //True indicates success; false indicates failure (the perturbation cannot be and thus is not applied)

//Unlike dihedrals, some perturbations may be complicated to adjust by small amounts
//(Easier to just apply with known parameter)


    public boolean applyPerturbation(float param) {
        //Apply the perturbation with a given value of the perturbation parameter
        //This operation can be reversed by undo(). It should be undoable even if dihedrals have been changed.
        //Because we only change the actualCoordinates here we do not update the curState in this function, just the curParam


        if( param == 0 ) { //Unperturbed state
            curParam = 0;
            return true;
        } else {
            //First store information needed to undo
            oldCoords = new ArrayList<HashMap<String, float[]>>(resAffected.length);
            double genChi1[] = getGenChi1();//Note the generalized chi1 for each residue so we can restore it after the perturbation

            for(int a=0; a<resAffected.length; a++)
                oldCoords.add(a, storeResBB(m.residue[resAffected[a]]) );

            boolean outcome = doPerturbationMotion(param);

            if(outcome) {
                if(idealizeSC)
                    outcome = idealizeSidechains();

                curParam = param;
            } else //Failed to apply perturbation
                curParam=0;

            setGenChi1(genChi1);

            return outcome;
        }

    }



    //This function undoes this perturbation's successor perturbations and then this one itself,
    //then applies the new value and reapplies the successors
    public void changePerturbationParameter(float param) {
        for(int a=successors.length-1; a>=0; a--)
            m.perts[successors[a]].undo();
        undo();
        applyPerturbation(param);
        for(int a=0; a<successors.length; a++)
            m.perts[successors[a]].applyPerturbation(m.perts[successors[a]].curParam);
    }


    public void undo() { //Undo the perturbation as it was most recently applied
        //We do not update the curParam now because we will typically reapply the perturbation shortly, often with the same parameter value

        if(curParam==0)//Nothing to do
            return;

        double[] genChi1 = getGenChi1();

        for(int a=0; a<resAffected.length; a++)
            restoreResBB( m.residue[resAffected[a]], oldCoords.get(a) );

        setGenChi1(genChi1);

    }



    //Return a step size for the minimizer to use when optimizing this perturbation
    //This should be overridden for any continuous perturbations...discrete perturbation cannot be minimized
    //so they don't need to override this (the returned value will not be used)
    public float getStepSizeForMinimizer() {
        return Float.NaN;
    }

    abstract public void setDefaultParams();//Set up the default perturbation states


    //Perturbations with extra info in the perturbation file
    //should override this function to output it
    public void writeExtraInfo(BufferedWriter bw) {
        return;
    }

    //Creating perturbations
    public static Perturbation generate(String type, Molecule m,int resList[],BufferedReader br) { //Generates a perturbation of the indicated type, affecting units in unitList
        //(Basically a constructor wrapper)
        Arrays.sort(resList);//We need the list of affected residues to be in ascending order

        if(type.equalsIgnoreCase("SHEAR"))
            return new Shear(m,resList);

        else if(type.equalsIgnoreCase("BACKRUB"))
            return new Backrub(m,resList);

        else if(type.equalsIgnoreCase("LOOP CLOSURE ADJUSTMENT"))
            return new LoopClosureAdjustment("LOOP CLOSURE ADJUSTMENT",m,resList);

        else if(type.equalsIgnoreCase("SSNE"))
            return new LoopClosureAdjustment("SSNE",m,resList);

        else if(type.equalsIgnoreCase("SSCE"))
            return new LoopClosureAdjustment("SSCE",m,resList);

        else if(type.equalsIgnoreCase("PARTIAL STRUCTURE SWITCH"))
            return new PartialStructureSwitch(m, resList, br);

        else if(type.equalsIgnoreCase("FULL STRUCTURE SWITCH"))
            return new FullStructureSwitch(m, resList, br);

        else if(type.equalsIgnoreCase("PROLINE FLIP"))
            return new ProlineFlip(m, resList);

        else { //not clear
            System.err.println("Perturbation type not identified: " + type);
            System.exit(1);
            return null;
        }
    }

    //Generate all the feasible perturbations of the given type in the molecule m
    //flexMolResNum = molecule residue numbers of flexible residues, in order
    //and resCount consecutive residues for each perturbation
    //Does not work for full structure switches, other perturbations requiring additional information
    //(i.e. the BufferedReader constructor argument)
    //Used by PerturbationSelector
    public static Perturbation[] generateAll(String type, Molecule m, int flexMolResNum[], int resCount) {

        ArrayList<Integer> startRes = new ArrayList<Integer>();//Molecule residue numbers to start perturbations at
        int numPerts = 0;

        //Look for any sets of resCount residues in a row in flexMolResNum
        for( int a=0; a<flexMolResNum.length - resCount + 1; a++ ) {
            int molResNum = flexMolResNum[a];

            boolean found = true;
            for(int b=1; b<resCount; b++) {
                if(flexMolResNum[a+b] != molResNum + b)
                    found = false;
            }

            if(found) {
                if( goodStart(type, molResNum, m) ) { //If molResNum is a good place to start a perturbation of the given type
                    startRes.add(molResNum);
                    numPerts++;
                }
            }
        }

        Perturbation[] ans = new Perturbation[numPerts];
        int place=0;

        for(int shift=0; shift<resCount; shift++) {
            //We add perturbations in order of their starting molecule residue number modulo resCount,
            //then in order of the actual molecule residue number, to minimize extra RCs from overlaps
            for(int pertNum=0; pertNum<numPerts; pertNum++) {
                if( startRes.get(pertNum) % resCount == shift ) {

                    int resList[] = new int[resCount];
                    for(int b=0; b<resCount; b++)
                        resList[b] = startRes.get(pertNum) + b;

                    Perturbation p = generate(type, m, resList, null);
                    p.setDefaultParams();
                    ans[place] = p;
                    place++;
                }

            }

        }



        return ans;
    }



    public static boolean goodStart(String type, int resNum, Molecule m) {
        //Is molecule residue number resNum of molecule m
        //a good place to start a perturbation of the given type?

        boolean good = true;

        if(type.equalsIgnoreCase("SHEAR")) { //We need at least three of this and the next 3 residues to be in a helix
            int count = 0;
            for(int shift=0; shift<4; shift++) {
                if( m.residue[resNum+shift].secondaryStruct == Residue.HELIX )
                    count++;
            }
            if(count < 3)
                good = false;
        }

        else if(type.equalsIgnoreCase("BACKRUB")) { //We need this and the next 2 residues to be in a sheet or loop
            for(int shift=0; shift<3; shift++) {
                if( m.residue[resNum+shift].secondaryStruct == Residue.HELIX )
                    good = false;
            }
        }

        else if(type.equalsIgnoreCase("LOOP CLOSURE ADJUSTMENT")) { //We need at least two of this and the next 2 residues to be in a loop
            int count = 0;
            for(int shift=0; shift<3; shift++) {
                if( m.residue[resNum+shift].secondaryStruct == Residue.LOOP )
                    count++;
            }
            if(count < 2)
                good = false;
        }

        else if(type.equalsIgnoreCase("SSCE")) { //Secondary structure C-terminal expansion
            //We need this and the next 2 residues to be in a loop, with a secondary structure element before the first one
            if( ! m.checkNBonded(resNum) )
                good = false;
            else if( m.residue[resNum-1].secondaryStruct == Residue.LOOP )
                good = false;

            for(int shift=0; shift<3; shift++) {
                if( m.residue[resNum+shift].secondaryStruct != Residue.LOOP)
                    good = false;
            }
        }

        else if(type.equalsIgnoreCase("SSNE")) { //Secondary structure N-terminal expansion
            //We need the next three residues to be in a loop, with a secondary structure element after the last one
            if( ! m.checkCBonded( resNum+3 ) )
                good = false;
            else if( m.residue[resNum+4].secondaryStruct == Residue.LOOP )
                good = false;

            for(int shift=1; shift<4; shift++) {
                if( m.residue[resNum+shift].secondaryStruct != Residue.LOOP)
                    good = false;
            }
        }

        else {
            System.err.println("Perturbation type " + type + " not recognized by goodStart: assuming it is infeasible");
            good = false;
        }


        return good;
    }





    //This function, to be used in applyPerturbation, idealizes all sidechains affected directly by the perturbation
    public boolean idealizeSidechains() {
        boolean outcome = true;

        //This modifies the actualCoordinates (like applyPerturbation)
        for(int a=0; a<resAffected.length; a++) {
            Residue res = m.residue[resAffected[a]];
            if( m.strand[res.strandNumber].isProtein )
                outcome = outcome && m.idealizeResSidechain(res);
        }

        return outcome;
    }


    public double[] getGenChi1() { //Get the generalized chi1 for each residue to make sure the perturbation does not change it
        //(This is either the real chi1 or another N-CA-CB-something dihedral to keep the SC oriented right)

        double ans[] = new double[resAffected.length];
        for(int a=0; a<resAffected.length; a++)
            ans[a] = getGenChi1(m, resAffected[a]);

        return ans;
    }


    public static double getGenChi1(Molecule m, int molResNum) { //Get the generalized chi1 given a molecule residue number
        //It's static because StrandRCs also uses it to place wild-type rotamers
        Residue res = m.residue[molResNum];

        if( res.name.equalsIgnoreCase("gly") || ( ! m.strand[res.strandNumber].isProtein ) )
            return 0;
        else {
            String lastAtom = getGenChi1LastAtom(res.name);
            int lastAtomNum = m.getMoleculeAtomNumber(lastAtom, molResNum);
            if(lastAtomNum == -1)//This might happen if there's 1HB instead of HB1 in alanine
                //Regardless, if present, 1HB is a good fourth atom for this dihedral in a strange residue/naming scheme
                lastAtomNum = m.getMoleculeAtomNumber("1HB",molResNum);
            if(lastAtomNum == -1) {
                System.err.println("Error: atom names not understood for " + res.fullName);
                System.exit(1);
            }


            return m.getTorsion( m.getMoleculeAtomNumber("N", molResNum),
                                 m.getMoleculeAtomNumber("CA", molResNum), m.getMoleculeAtomNumber("CB", molResNum),
                                 lastAtomNum );
        }
    }


    public void setGenChi1(double ang[]) { //Set each residue's generalized chi1 to the specified values

        for(int a=0; a<resAffected.length; a++)
            setGenChi1(m, resAffected[a], ang[a]);
    }


    public static void setGenChi1(Molecule m, int molResNum, double ang) {
        Residue res = m.residue[molResNum];

        if( m.strand[res.strandNumber].isProtein && ( ! res.name.equalsIgnoreCase("gly") ) && ( ! res.name.equalsIgnoreCase("pro") ) ) {
            //Glycine doesn't have a generalized chi1, and we cannot freely change proline's (sidechain idealization sets
            //a chi1 for proline that should remain in place)

            String lastAtom = getGenChi1LastAtom(res.name);
            int lastAtomNum = m.getMoleculeAtomNumber(lastAtom, molResNum);
            if(lastAtomNum == -1) {
                lastAtom="1HB";
                lastAtomNum = m.getMoleculeAtomNumber("1HB",molResNum);
            }
            if(lastAtomNum == -1) {
                System.err.println("Error: atom names not understood for " + res.fullName);
                System.exit(1);
            }

            //get a list of atoms to rotate...based roughly on Molecule.getAtomList
            int listSize = res.atom.length-8;
            if(res.getAtomByName("H3") != null)//Checking if this is an N-terminal residue
                listSize -= 2;
            else if(res.getAtomByName("OXT") != null)//Checking if C-terminal
                listSize--;


            int atomList[] = new int[listSize];
            int c=0;

            for(int b=0; b<res.atom.length; b++) {

                Atom at = res.atom[b];

                if( ! ( at.isBBatom || at.name.contains("HA") ||
                        at.name.equalsIgnoreCase("CB") || at.name.equalsIgnoreCase(lastAtom) ) ) {
                    atomList[c] = res.atom[b].moleculeAtomNumber;
                    c++;
                }
            }


            m.setTorsion( m.getMoleculeAtomNumber("N", molResNum),
                          m.getMoleculeAtomNumber("CA", molResNum), m.getMoleculeAtomNumber("CB", molResNum),
                          lastAtomNum, ang, atomList, listSize, false );
        }
    }


    public static String getGenChi1LastAtom(String resName) { //Get atom name X where generalized chi1 is N-CA-CB-X

        if( resName.equalsIgnoreCase("val") || resName.equalsIgnoreCase("ile") )
            return "CG1";
        else if( resName.equalsIgnoreCase("ser") )
            return "OG";
        else if( resName.equalsIgnoreCase("thr") )
            return "OG1";
        else if( resName.equalsIgnoreCase("cys") || resName.equalsIgnoreCase("cyx") )
            return "SG";
        else if( resName.equalsIgnoreCase("ala") )
            return "HB1";
        else
            return "CG";
    }

    public HashMap<String, float[]> storeResBB( Residue res ) {
        return storeResBB( res, m );
    }

    public HashMap<String, float[]> storeResBB(Residue res, Molecule molec) {
        //Store the residue's backbone coordinates and sidechain orientation in a hashmap
        //We allow it to store from a different molecule because structure-switching perturbation need to do that

        HashMap<String, float[]> resCoord = new HashMap<String, float[]>();


        Atom CA = res.getAtomByName("CA");

        for(int b=0; b<res.atom.length; b++) {
            Atom at = res.atom[b];
            if( at.isBBatom || at.name.contains("HA") ||        //backbone atom or alpha hydrogen
                    ( res.name.equalsIgnoreCase("PRO") && at.name.equalsIgnoreCase("CD") ) ) { //or proline CD (treated like part of the backbone, as an amide H would be)
                //Also we store coordinates for all atoms in proline because it doesn't have adjustable dihedrals
                float coords[] = new float[3];
                System.arraycopy( molec.actualCoordinates, 3*at.moleculeAtomNumber, coords, 0, 3);
                resCoord.put(at.name, coords);
            } else if( at.name.equalsIgnoreCase("CB") ) {
                float coords[] = new float[3];
                for(int c=0; c<3; c++)
                    coords[c] = molec.actualCoordinates[ 3*at.moleculeAtomNumber + c ]
                                - molec.actualCoordinates[ 3*CA.moleculeAtomNumber + c ] ;
                resCoord.put("CACB", coords);//Vector from CA to CB
            }
        }

        return resCoord;
    }




    public void restoreResBB( Residue res, HashMap<String, float[]> resCoord ) {
        //Restore the residue's backbone coordinates and sidechain orientation as indicated by the hashmap

        RotMatrix r = new RotMatrix();
        float SC_mtx[][] = new float[3][3], origCA[] = new float[3], curCA[] = new float[3];

        if( res.name.equalsIgnoreCase("gly") && resCoord.containsKey("CACB") ) { //There's been a mutation to glycine since the coordinates were stored
            //This might happen in a structure switch for example.
            float CAHA[] = r.scale( resCoord.get("CACB"), 0.716f);//This is 1.1/1.536: based on Molecule.idealizeResSidechain
            float HA3[] = r.add(CAHA, resCoord.get("CA") );
            resCoord.put("HA3", HA3 );
            resCoord.put( "HA2", resCoord.get("HA") );
            resCoord.remove("CACB");//Remove the CACB vector so the sidechain restoration routine is not called
        }

        else if ( ( ! res.name.equalsIgnoreCase("gly") ) && ( ! resCoord.containsKey("CACB") ) ) { //There's been a mutation away from glycine
            float HA[] = resCoord.get("HA3");
            float CACB[] = r.scale( r.subtract( HA, resCoord.get("CA") ) , 1.396f );//1.536/1.1
            resCoord.put("CACB", CACB);
            resCoord.put( "HA", resCoord.get("HA2") );
        }


        if( res.name.equalsIgnoreCase("PRO") && ( ! resCoord.containsKey("CD") ) ) { //There's been a mutation to proline since the coordinates were stored
            //This might happen in a structure switch for example.
            float NCD[] = r.subtract( resCoord.get("H"), resCoord.get("N") );
            NCD = r.scale( NCD, getIdealNHCD("PRO") / r.norm(NCD) );
            float CD[] = r.add( NCD, resCoord.get("N") );
            resCoord.put("CD", CD );
        }

        else if ( ( ! res.name.equalsIgnoreCase("PRO") ) && ( ! resCoord.containsKey("H") ) ) { //There's been a mutation from proline
            float NH[] = r.subtract( resCoord.get("CD"), resCoord.get("N") );
            NH = r.scale( NH, getIdealNHCD(res.name) / r.norm(NH) );
            float H[] = r.add( NH, resCoord.get("N") );
            resCoord.put("H", H);
        }


        if( resCoord.containsKey("CACB") ) { //We have a sidechain to restore
            origCA = resCoord.get("CA");//original (and thus post-undoing) CA coordinates
            curCA = m.getActualCoord( m.getMoleculeAtomNumber("CA", res.moleculeResidueNumber ) );//Current (pre-undoing) CA coordinates
            float origCACB[] = resCoord.get("CACB");
            float curCACB[] = r.subtract( m.getActualCoord( m.getMoleculeAtomNumber("CB", res.moleculeResidueNumber ) ) , curCA );

            //Create a matrix for rotating curCACB to origCACB (to superimpose the CBs)
            //Axis for this will be curCACB X origCACB (sign for this will give an angle < 180 degrees)
            //This is based on RotMatrix.getSuperposingRotMatrix
            float axis1[] = r.cross(curCACB, origCACB);
            if( r.norm(axis1) == 0 ) { //uold and unew are collinear
                if( r.dot(curCACB, origCACB) > 0 )//uold and unew are already superimposed
                    SC_mtx = r.identity();
                else { //Need a 180-degree rotation.
                    float normal[] = r.getPerpendicular(curCACB);
                    r.getRotMatrix( normal[0], normal[1], normal[2], 180, SC_mtx );
                }
            } else {
                float th = r.getAngle(curCACB, origCACB);//angle of this first rotation
                r.getRotMatrixRad(axis1[0], axis1[1], axis1[2], th, SC_mtx);//Rotation matrix for the sidechain
            }

        }

        for(int b=0; b<res.atom.length; b++) { //Now put the coordinates back
            Atom at = res.atom[b];
            if( resCoord.containsKey(at.name) ) { //Restore the old coordinates exactly
                float coords[] = resCoord.get(at.name);
                System.arraycopy( coords, 0, m.actualCoordinates, 3*at.moleculeAtomNumber, 3);
            } else { //Use SC_mtx to get it back in place
                float curCoords[] = m.getActualCoord( at.moleculeAtomNumber );
                float undoneCoords[] = r.add ( origCA, r.applyRotMatrix(SC_mtx, r.subtract(curCoords, curCA) ) );
                System.arraycopy( undoneCoords, 0, m.actualCoordinates, 3*at.moleculeAtomNumber, 3);
            }
        }

        if( res.name.equalsIgnoreCase("PRO") )
            m.idealizeProRing(res);//This will fix the proline geometry
    }



    /*
        //*******This function from Residue is probably better performed here (maybe in a different organization)
        public void makeAllPertStates(Molecule m,RotamerLibrary rl,int resIndex){//Given that perts is filled in, make pertStates and defaultParamAppl to allow all perturbation states
                //Uses a handle to the molecule (m), rotamer library (rl), and strandRCs object (sysLR)
                //and the residue number of this residue in the strand

                if(perts==null)//Don't need to do anything
                    return;

                int defnum[]=new int[perts.length];//Number of default states for individual perturbations
                int statenum=1;//Number of overall perturbation states
                for(int a=0;a<perts.length;a++){
                    defnum[a]=m.perts[perts[a]].storedParams.length;//Number of default values for the a'th perturbation affecting this unit
                    statenum*=defnum[a];
                }

                pertStates=new int[statenum][perts.length];
                int stateVec[]=new int[perts.length];//States for each of the perturbations
                for(int a=0;a<statenum;a++){//Each possible combination of default parameter values for the perturbations in perts gets its own overall perturbation state
    //Increment statevec
                    int toChange=perts.length-1;
                    while(stateVec[toChange]==defnum[toChange])
                        toChange--;
                    stateVec[toChange]++;
                    toChange++;
                    while(toChange<perts.length){
                        stateVec[toChange++]=0;
                    }

                    pertStates[a]=stateVec;
                }

                int numAAsAllowed = sRC.getNumAllowable(resIndex);
                defaultParamAppl=new boolean[statenum][numAAsAllowed][];
                int numRot;
                for(int a=0;a<statenum;a++){
                    for(int b=0;b<numAAsAllowed;b++){
                        if(sRC.getIndexOfNthAllowable(resIndex, b) == -1)//The unit is not a sidechain in the rotamer library
                            numRot=1;
                        else
                            numRot=rl.getNumRotForAAtype(sRC.getIndexOfNthAllowable(resIndex, b));
                        defaultParamAppl[a][b]=new boolean[numRot];
                        for(int c=0;c<numRot;c++)
                            defaultParamAppl[a][b][c]=true;
                    }
                }

            }
     *
     */


    private float getIdealNHCD(String resType) {
        //Returns the ideal N-H bond length for the given residue type (or N-CD for Pro)
        //Used in coordinate restoration with changes to and from proline
        Residue ideal = new Amber96PolyPeptideResidue().getResidue(resType);
        float N[] = ideal.getAtomByName("N").coord;
        float HCD[];
        if(resType.equalsIgnoreCase("PRO"))
            HCD = ideal.getAtomByName("CD").coord;
        else
            HCD = ideal.getAtomByName("H").coord;

        RotMatrix rm = new RotMatrix();
        return rm.norm( rm.subtract(N, HCD) );
    }



}