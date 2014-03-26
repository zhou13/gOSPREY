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
//	PMinimizer.java
//
//	Version:           2.1 beta
//
//
//	  authors:
// 	  initials    name                 organization                email
//	 ---------   -----------------    ------------------------    ----------------------------
//	  MAH           Mark A. Hallen	  Duke University               mah43@duke.edu
///////////////////////////////////////////////////////////////////////////////////////////////

import java.util.HashMap;
import java.util.HashSet;

/**
 * This class implements a simple energy minimization routine for the side-chains and perturbations. Handles the computation
 * of the Amber dihedral energy penalties (for minimizing away from the initial rotamer dihedrals).
 * Additionally there is a special residue, the ligand, that can translate and globally rotate.
 * The class extends SimpleMinimizer, relying on some of that class' functions for dihedral minimization
 */
public class PMinimizer extends SimpleMinimizer {

	// If the debug flag is set to true then additional debug statements are
	//  printed to standard out.

        boolean minimizePerturbations;//If true, minimize perturbations in addition to the usual sidechain dihedrals
        //The energies of the residues whose energy calculation flags will be minimized
        //with respect to the dihedrals of the flexible residues, plus the perturbations affecting them


	//sysRH and ligRH need to be StrandRCs objects


        int numPerturbations = 0;//number of perturbations; also the size of the next 5 arrays

        int pertList[];//List of perturbations (indices in m.perts)
        float pertParamValues[];//Current parameter values for the perturbations

        float pertStepSize[];//Step sizes for each perturbation
        float pertParamMax[];//Maxima for each perturbation's parameter
        float pertParamMin[];//Minima for each perturbation's parameter


        int flexResMap[];//Map from molecule residue number to index among flexible 
        //residues (and thus in the partial arrays for energy calculations)
        //(needed because perturbations' resAffected arrays contain molecule residue numbers)


        boolean checkMonotonic = false;//This activates an alternate stopping condition that detects when the minimization stops improving
        double energyTol = 0.01;//This is about the roundoff energy in many cases...when we are getting energy changes less than this
        //and the energy has stopped monotonically decreasing, we are as close to the minimum as we are likely to get, so we can stop

	PMinimizer(boolean minPerts){
		minimizePerturbations = minPerts;
	}
	


    @Override
	public void initialize(Molecule theM, int numStrands, Amber96ext theA96ff,
		StrandRotamers strandRotamers[], int curAANum[], boolean doDihedral){

                super.initialize(theM, numStrands, theA96ff, strandRotamers, curAANum, doDihedral);//Call the SimpleMinimizer no-ligand initialize() function

                if(minimizePerturbations){
                    flexResMap = new int[m.residue.length];
                    setupPerturbations(false);
                }
	}	




        public void setupPerturbations(boolean hasLig){//Set up perturbation information

                HashSet<Integer> pertSet=new HashSet<Integer>();

                int flexResCount = 0;
                numPerturbations = 0;


                //Find the perturbations affecting flexible residues
                for(int i=0;i<numberOfStrands;i++){
			for(int j=0;j<m.strand[i].numberOfResidues;j++){

                                Residue localResidue = m.strand[i].residue[j];

                                if(localResidue.flexible){
                                    for(int a=0;a<localResidue.perts.length;a++){//Add this residue's perturbations to the perturbation set, if they're not already in it

                                        if( !pertSet.contains(localResidue.perts[a]) ){
                                            pertSet.add(localResidue.perts[a]);
                                            numPerturbations++;
                                        }
                                    }

                                    flexResMap[localResidue.moleculeResidueNumber] = flexResCount++;
                                }
			}
		}


                //Finish preparing perturbation info
                pertList=new int[numPerturbations];
                pertParamValues=new float[numPerturbations];
                pertStepSize=new float[numPerturbations];
                pertParamMax=new float[numPerturbations];
                pertParamMin=new float[numPerturbations];

                int a=0;
                for(Integer pertNum : pertSet){
                    pertList[a++]=pertNum.intValue();
                }

                for(int b=0;b<numPerturbations;b++)
                    pertStepSize[b]=m.perts[pertList[b]].getStepSizeForMinimizer();
        }


        public void setupPertStates(){//Once the molecule is in the right RC,
            //this function will set up information for the perturbation states
            //within which the PMinimizer will minimize

            for(int b=0;b<numPerturbations;b++){
                    Perturbation pert=m.perts[pertList[b]];
                    int paramValueNum=pert.curState;//Which stored parameter state we are starting at
                    pertParamMax[b]=pert.maxParams[paramValueNum];
                    pertParamMin[b]=pert.minParams[paramValueNum];
                    pertParamValues[b]=pert.curParam;//This should be (pertParamMax[b] + pertParamMin[b])/2f if we're in the RC specified by curState
            }

        }
	
	
////////////////////////////////////////////////////////////////////////////////
//  Minimization Section
////////////////////////////////////////////////////////////////////////////////



        private void doPerturbationStep (int pertNum, float stepSize){//Calculates and applies the next perturbation parameter value in the steepest-descent minimization
            
                double initialEnergy, secondEnergy, thirdEnergy;
                Perturbation pert = m.perts[pertList[pertNum]];

                initialEnergy = getMultiResidueEnergy(pert.resAffected);//Just need energy for the units affected by the perturbation with
                //energy calculation flags set to true

		// Apply the perturbation with positive step, then negative step, and compute energies
                pert.changePerturbationParameter(pertParamValues[pertNum]+stepSize);
		secondEnergy = getMultiResidueEnergy(pert.resAffected);
                pert.changePerturbationParameter(pertParamValues[pertNum]-stepSize);
		thirdEnergy = getMultiResidueEnergy(pert.resAffected);

                //Now compute the change in the perturbation parameter, change pertParamValues appropriately
		if ((initialEnergy > secondEnergy)&&(initialEnergy > thirdEnergy)){
			if (secondEnergy < thirdEnergy)
				pertParamValues[pertNum]+=stepSize;
			else
				pertParamValues[pertNum]-=stepSize;
		}
		else if(initialEnergy > secondEnergy)
			pertParamValues[pertNum]+=stepSize;//Perhaps the difference in energies could be used to compute & return a more accurate estimate
		else if(initialEnergy > thirdEnergy)
			pertParamValues[pertNum]-=stepSize;
                //else leave pertParamValues[pertNum] unchanged: energy is at a local minimum with respect to this parameter


               checkCumulativePerturbation(pertNum);//Makes sure pertParamValues[pertNum] is within the allowed bounds, fix it if it isn't

               //Now re-apply the perturbation with the new parameter value
               pert.changePerturbationParameter(pertParamValues[pertNum]);

        }


        private void checkCumulativePerturbation(int pertNum){//,double maxPertParam,double minPertParam
            //Keep the parameter for the indicated perturbation within the allowed bounds

            if(pertParamValues[pertNum] > pertParamMax[pertNum])
                pertParamValues[pertNum] = pertParamMax[pertNum];
            else if(pertParamValues[pertNum] < pertParamMin[pertNum])
                pertParamValues[pertNum] =  pertParamMin [pertNum];

        }
        


	// Performs a simple steepest descent minimization
	// Assumptions:
	//  -a96ff is current
	//     all atoms of appropriate units have been assigned
	//     initializeEVCalculation has been called
	//  -m.actualCoords contains current atomic coordinates
	// The user should have set numMinimizationSteps and
	//  initialAngleStepSize (else defaults will be used)
	// Uses specific precomputed nonbonded arrays for each unit (this makes things
	//  run faster)
	// If ligandOnly is true, then only the ligand is allowed to minimize, while
	//		the system units are fixed
    @Override
        public void minimize(int numSteps){

                setupPertStates();//Make sure the correct perturbation states are assigned to each residue

                float tempPertStep[] = new float[numPerturbations];//Temporary step sizes, which will start at the pertStepSize values
                //and then will be reduced over the course of minimization
                float deltaPertStepSize[] = new float[numPerturbations];
                if(minimizePerturbations){
                    System.arraycopy(pertStepSize, 0, tempPertStep, 0, numPerturbations);
                    for(int a=0;a<numPerturbations;a++)
                       deltaPertStepSize[a] = tempPertStep[a] / numSteps;
                }

                
		float step = initialAngleStepSize;
		double lmaxMovement = maxMovement;
			// maximum degrees by which a torsion can
			//  cumulatively change
		float strRotStep = RotStep;
			// step size for the rigid rotation
			//  of the ligand
		float strTransStep = TransStep;
			// step size in � for the rigid ligand
			//  translation
		double strMaxTrans = MaxTrans;
			// the maximum ligand translation allowed

		int i=0;
		boolean done = false;
		double strTorque[] = new double[3];
		double strTrans[] = new double[3];

		strDihedDiff = new double[numberOfStrands][];
		strCumulativeDihedStep = new double[numberOfStrands][];
		for(int str=0; str<numberOfStrands; str++){
			strDihedDiff[str] = new double[numStrDihedrals[str]];
			strCumulativeDihedStep[str] = new double[numStrDihedrals[str]];
		}
		//sysDihedDiff = new double[numSysDihedrals];
		//ligDihedDiff = new double[numLigDihedrals];
		//sysCumulativeDihedStep = new double[numSysDihedrals];
		//ligCumulativeDihedStep = new double[numLigDihedrals];

		strStartCOM = new double[numberOfStrands][3];
		strCurCOM = new double[numberOfStrands][3];
		for(int str=0;str<numberOfStrands;str++){
			strStartCOM[str] = m.getStrandCOM(str);
			strCurCOM[str][0] = strStartCOM[str][0];
			strCurCOM[str][1] = strStartCOM[str][1];
			strCurCOM[str][2] = strStartCOM[str][2];
		}
		/*if(ligStrNum != -1){
			// get the staring COM
			ligStartCOM = m.getStrandCOM(ligStrNum);
			ligCurCOM[0] = ligStartCOM[0];
			ligCurCOM[1] = ligStartCOM[1];
			ligCurCOM[2] = ligStartCOM[2];
		}*/

		// Initialize the dihedral movement arrays
		for(int str=0; str<numberOfStrands;str++){
			for(int j=0; j<numStrDihedrals[str];j++){
				strDihedDiff[str][j] = 0.0;
				strCumulativeDihedStep[str][j] = 0.0;
			}
		}

		/*for(int j=0;j<numSysDihedrals;j++){
			sysDihedDiff[j] = 0.0;
			sysCumulativeDihedStep[j] = 0.0;
		}
		for(int j=0;j<numLigDihedrals;j++){
			ligDihedDiff[j] = 0.0;
			ligCumulativeDihedStep[j] = 0.0;
		}*/

		// If computing dihedral energies initialize them
		if (doDihedEnergy){
			if (!setupDihedralTerms()) //could not initialize dihed energies
				System.exit(1);
		}

		float deltaStep = step / numSteps;
		float deltaRotStep = strRotStep / numSteps;
		float deltaTransStep = strTransStep / numSteps;

			// numFlexRes, flexResAtomList, and flexResListSize include the ligand if one exists
			/*if(ligStrNum != -1)
				a96ff.setupPartialArrays(numFlexRes+2,MAX_NUM_ATOMS_DISTAL,flexResAtomList,
					flexResListSize);
			else*/
		a96ff.setupPartialArrays(totalFlexRes+totalTransRotStrands,MAX_NUM_ATOMS_DISTAL,flexResAtomList,
					flexResListSize);

               //Prepare for checkMonotonic stopping condition
               double bestEnergy = a96ff.calculateTotalEnergy(m.actualCoordinates,-1)[0];
               if(doDihedEnergy)
                   bestEnergy += computeDihedEnergy();
               boolean monotonic = true;//Indicates that energy is still descending monotonically; near convergence this might
               //stop happening (e.g. because we try to take a step to a lower energy but are blocked by the bounds on a parameter, ending up in a slightly higher-energy state)

               while(!done){

                    //If(minimizePerturbations), at the beginning of the loop, we are at the unperturbed starting conformation
                   //Else we're in the RC (perturbations are applied and rotamers are set)

			for(int str=0; str<numberOfStrands;str++){
				for(int j=0;j<numStrDihedrals[str];j++) {
					strDihedDiff[str][j] = computeDihedDiff(strDihedralAtNums[str][j],strDihedralDistal[str][j],
						strNumAtomsDistal[str][j],strDihedToResNum[str][j], step, str, j);
					updateCumulative(strCumulativeDihedStep[str],strDihedDiff[str],j,lmaxMovement);
					applyDihedStep(strDihedralAtNums[str][j],strDihedralDistal[str][j],strNumAtomsDistal[str][j],strDihedDiff[str][j]);
				}
			}

			/*for(int j=0;j<numLigDihedrals;j++) {
				ligDihedDiff[j] = computeDihedDiff(ligDihedralAtNums[j],ligDihedralDistal[j],
					ligNumAtomsDistal[j],ligResNumber,step,true,j);
				updateCumulative(ligCumulativeDihedStep,ligDihedDiff,j,lmaxMovement);
				applyDihedStep(ligDihedralAtNums[j],ligDihedralDistal[j],ligNumAtomsDistal[j],ligDihedDiff[j]);
			}*/

			//Translate and rotate the ligand
			/*if(ligStrNum != -1)
				doLigTransRot(ligTorque, ligTrans, lligRotStep, lligTransStep, lligMaxTrans);*/

			for(int str=0;str<numberOfStrands;str++){
				if(m.strand[str].rotTrans){
					doStrTransRot(str, strTorque, strTrans, strRotStep, strTransStep, strMaxTrans);
				}
			}


                        if(minimizePerturbations){

                            for(int j=0;j<numPerturbations;j++){
                                if(pertParamMax[j] > pertParamMin[j])//There is an interval of parameter values for this perturbation that we can minimize within
                                    doPerturbationStep(j, tempPertStep[j]);//Computes, checks, and applies step
                            }

                            for(int a=0;a<numPerturbations;a++)
                                tempPertStep[a] -= deltaPertStepSize[a];
                        }


			/*if(debug){
				String filename = String.format("badVdw/run%1$d/structure%2$02d.pdb",GLOBALNUM,i);
				//Write out pdbs of minimization movement
				m.saveMolecule(filename, 0.0f);
			}*/

			i++;
			if(i>=numSteps){
				done = true;
			}

			step -= deltaStep;
			strRotStep -= deltaRotStep;
			strTransStep -= deltaTransStep;


                        if(checkMonotonic){
                            double newEnergy = a96ff.calculateTotalEnergy(m.actualCoordinates,-1)[0];
                            if(doDihedEnergy)
                                newEnergy += computeDihedEnergy();

                            if( (!monotonic) && ( Math.abs(newEnergy - bestEnergy) < energyTol ) )//In this situation we're as close to convergence as we'll get
                                done = true;

                            if(newEnergy > bestEnergy)
                                monotonic = false;
                            else
                                bestEnergy = newEnergy;
                        }
		}

		clearMolGradient(); //after minimization is done, clear the molecule gradient

		if(debug){
			// Display movement
			for(int str=0; str<numberOfStrands;str++){
				for(int j=0;j<numStrDihedrals[str];j++) {
					System.out.print(strCumulativeDihedStep[str][j] + " ");
				}
			}

			GLOBALNUM++;

			System.out.println();
		}
	}

	
////////////////////////////////////////////////////////////////////////////////
//	 End of Minmization Section
////////////////////////////////////////////////////////////////////////////////

        private double getMultiResidueEnergy(int molResNum[]){
            //Given a set of residues with specified molecule residue numbers,
            //calculates the energy of those with energy calculation flags set to true,
            //using partial arrays if this subset is just one residue
            //(THIS MAY BE TOO INEFFICIENT--IF SO THE PARTIAL ARRAYS NEED TO BE REORGANIZED)
            double totEnergy[];

            int calcRes = -1;
            for(int a=0; a<molResNum.length; a++){

                if( ! m.residue[molResNum[a]].validConf )//If the perturbation change whose energy we're getting caused an invalid conformation, we apply an infinite energy penalty
                        return Float.POSITIVE_INFINITY;

                boolean needEnergy = m.residue[molResNum[a]].getEnergyEvalBB() || m.residue[molResNum[a]].getEnergyEvalSC();//Do we need the energy for this residue?
                if(needEnergy){
                    if(calcRes == -1)
                        calcRes = flexResMap[molResNum[a]];
                    else{//Need energy for more than one residue
                        totEnergy = a96ff.calculateTotalEnergy(m.actualCoordinates,-1);
                        return totEnergy[0];
                    }
                }
            }

            if(calcRes == -1)//Don't need the energy for any of the residues so return zero
                return 0;
            else{
                totEnergy = a96ff.calculateTotalEnergy(m.actualCoordinates, calcRes);
                return totEnergy[0];
            }
        }


}
