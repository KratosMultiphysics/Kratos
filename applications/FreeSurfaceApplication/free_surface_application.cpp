//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Antonia Larese
//


// System includes


// External includes


// Project includes
#include "free_surface_application.h"
#include "free_surface_application_variables.h"


namespace Kratos {

KratosFreeSurfaceApplication::KratosFreeSurfaceApplication():
    KratosApplication("FreeSurfaceApplication")
    {}

void KratosFreeSurfaceApplication::Register() {
 	// calling base class register to register Kratos components
 	KratosApplication::Register();
 	std::cout << "Initializing KratosFreeSurfaceApplication... " << std::endl;



}
}  // namespace Kratos.
