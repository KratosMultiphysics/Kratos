/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    pimpleDyMFoam.C

Description
    Transient solver for incompressible, turbulent flow of Newtonian fluids
    on a moving mesh.

    Turbulence modelling is generic, i.e. laminar, RAS or LES may be selected.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "CorrectPhi.H"
#include "fvOptions.H"
#include "fixedValuePointPatchField.H"
#include "IOstream.H"
#include "Ostream.H"
#include "forces.H"

#include <stdlib.h>
#include <assert.h>
#include <set>
#include <string>
#include <map>
#include <sstream>
#include <unistd.h>
#include <fstream>
#include <iostream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

	#include "postProcess.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createControls.H"
    #include "createFields.H"
    #include "createUf.H"
    #include "createFvOptions.H"
	#include "CourantNo.H"
	#include "setInitialDeltaT.H"
	
    
	turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	
	
    while (runTime.run()) // or True?
    {
        usleep(1000); // Expressed in microseconds 
        
    	if (exists("next.coco"))
		{
    		remove("next.coco");
    		
			#include "readControls.H"
			#include "CourantNo.H"
			#include "setDeltaT.H"
    		
    		runTime++;
    		OFstream outfile ("next_ready.coco");
    		outfile << "Joris says: good job on next.coco" << endl;
			Info << "Time = " << runTime.timeName() << nl << endl; // Might be deleted when linked to CoCoNuT (which already outputs current time step)
		}
    	
    	if (exists("continue.coco"))
		{
    		remove("continue.coco");
    		// Define movement of the coupling interface
    		
    		/*label patchID = mesh.boundaryMesh().findPatchID();
    		const polyPatch& cPatch = mesh.boundaryMesh()[patchID];
    		forAll (cPatch,faceI)
    		{
    		pointDisplacement.boundaryField()[patchID][faceI] = ;
    		}*/
    		// See EMPIRE coupling code
    		
    		
    		
    		
    		// Calculate the mesh motion and update the mesh
            mesh.update();

            // Calculate absolute flux from the mapped surface velocity
            phi = mesh.Sf() & Uf;

            if (mesh.changing() && correctPhi)
            {
                #include "correctPhi.H"
            }

            // Make the flux relative to the mesh motion
            fvc::makeRelative(phi, U);

            if (mesh.changing() && checkMeshCourantNo)
            {
                #include "meshCourantNo.H"
            }

            // --- Pressure-velocity PIMPLE corrector loop
            while (pimple.loop())
            {
                #include "UEqn.H"

                // --- Pressure corrector loop
                while (pimple.correct())
                {
                    #include "pEqn.H"
                }
                
                if (pimple.turbCorr())
                {
                    laminarTransport.correct();
                    turbulence->correct();
                }
            }
            
            // Return the coupling interface output

            Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
                << "  ClockTime = " << runTime.elapsedClockTime() << " s"
                << nl << endl;
            
            IOobject controlDict_IO = IOobject("controlDict", runTime.system(),mesh,IOobject::MUST_READ,IOobject::AUTO_WRITE);
            IOdictionary controlDict(controlDict_IO);
            controlDict.Foam::regIOobject::write();
            
            OFstream outfile ("continue_ready.coco");
        	outfile << "Joris says good job on continue.coco" << endl;
			            
    		
		}
		
    	if (exists("save.coco"))
		{
    		remove("save.coco");
    		runTime.write(); // OF-command: loops over all objects and requests writing - writing is done based on the specific settings of each variable (AUTO_WRITE, NO_WRITE)
    		OFstream outfile ("save_ready.coco");
			outfile << "Joris says: good job on save.coco" << endl;
		}
    	
    	if (exists("stop.coco"))
		{
    		remove("stop.coco");
    		OFstream outfile ("stop_ready.coco"); 
    		outfile << "Joris says: good job on stop.coco" << endl;
    		break;
		}  
    }
    return 0;
}

// ************************************************************************* //
