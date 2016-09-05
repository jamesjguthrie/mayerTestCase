/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    edcPisoFoam (2015aSRx) for OF 2.3.1

Description
    edcPisoFoam is a specialized application for turbulent combustion modeling.
    A fully compressible solver incorporates the conventional URANS formulation 
    as well as Large-Eddy simulation. The Eddy Dissipation Concept with a detailed chemistry approach
    is used for the turbulence-chemistry interaction. A robust implicit Runge-Kutta method (RADAU5) 
    for integrating stiff ordinary differential equations allows to calculate chemical systems of any type of complicity.

    support: dmitry.lysenko@decgroup.org 
				

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "volFields.H"
#include "bound.H"
#include "rhoCombustionModel.H"  // instead of "rhoChemistryCombustionModel.H" //


#include "turbulenceModel.H"
#include "multivariateScheme.H"
#include "radiationModel.H" 
#include "ignition.H" 

#include "fvIOoptionList.H" //added

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

    Info<< "\n edcPisoFoam build in QT IDE" << endl;


    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "readGravitationalAcceleration.H"
    #include "createFields.H"

    #include "createFvOptions.H" //added

    #include "createRadiationModel.H"     

    #include "initContinuityErrs.H"
    #include "readTimeControls.H"
    #include "compressibleCourantNo.H"
    #include "setInitialDeltaT.H"



    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    // mixture properties 
    Info<< "\nMixture properties\n" << endl;
    // Laminar diffusion coefficient for species transport
    // DY DY [0 2 -1 0 0 0 0] 2.88e-5;
    //scalar DY = 2.88e-5;

    dimensionedScalar DY
    (
	"DY",
	dimensionSet(0,2,-1,0,0,0,0),
	scalar(2.88e-5)
     ); 

    Info<< "\nLaminar diffusion coefficient\n" << DY << endl;
    // turbulent Schmidt number
    // v = 2.083 x 10^-5 (2.08e-5)  m2/s - kinematic viscosity
    // Sct = v/DY
    // Sct Sct [0 0 0 0 0 0 0] 0.723;
    scalar Sct = 0.72; 
    Info<< "\nTurbulent Schmidt number, Sc = " << Sct << endl;
     



    // check mass fraction conservation on principle intet and outlet boundaries 
    #include "checkYi.H"


	
    Info<< "\nStarting time loop\n" << endl;
    bool firstTimeStep = true;	


    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "readPISOControls.H"

        const bool edcLES = pisoDict.lookupOrDefault<bool>("edcLES", 0);
        dimensionedScalar rhoMax(pisoDict.lookup("rhoMax"));
        dimensionedScalar rhoMin(pisoDict.lookup("rhoMin"));


        /* we do not need these staff for OF 2.3.1*/
        /*
        const bool solverIterOutput  = pisoDict.lookupOrDefault<bool>("solverOutput", 1);
    	if (!solverIterOutput)
            lduMatrix::debug = 0;
        */

        /*
         * setting solverPerformance == 0 & lduMatrix = 0 in OpenFOAM-2.3.1/etc/controlDic
        */



        #include "compressibleCourantNo.H"
        #include "setDeltaT.H"

 

        runTime++;
        Info<< "Time = " << runTime.timeName() << nl << endl;

      
	

        /////////////////////////////////////////////////////////////////////////
        //
        //   we are not solving chemistry for the first time step!!!!!!!!!!!!!!
        //
        /////////////////////////////////////////////////////////////////////////
	
        if (firstTimeStep )
        {
            firstTimeStep = false;
        }
        else
        {
            #include "gammaStarEqn.H"
            combustion->correct();
        }
	
	
        rho = thermo.rho();

        #include "rhoEqn.H"


        for (label ocorr=1; ocorr <= nOuterCorr; ocorr++)
        {
            #include "UEqn.H"
            #include "multivariateConvection.H" //has to be called before YEqn.H and hsEqn.H 

	    
 
            #include "YEqn.H"
            #include "heEqn.H"

            // --- PISO loop
            for (int corr=1; corr<=nCorr; corr++)
            {
                #include "pEqn.H"
            }
        } // end for


        turbulence->correct();
        #include "ignite.H"
	
        rho = thermo.rho();
        rho = max(rho, rhoMin);
        rho = min(rho, rhoMax);
        thermo.rho() = rho;
        //rho.relax();

        Info<< "rho max/min : " << max(rho).value() << " " << min(rho).value() << endl;

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"   << nl << endl;


    } // end while

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
