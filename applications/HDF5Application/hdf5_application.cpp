//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//


// System includes


// External includes


// Project includes
#include "hdf5_application.h"
#include "hdf5_application_variables.h"


namespace Kratos {

KratosHDF5Application::KratosHDF5Application() : KratosApplication("HDF5Application") {}

void KratosHDF5Application::Register() {
 	// calling base class register to register Kratos components
 	KratosApplication::Register();
 	std::cout << "Initializing KratosHDF5Application... " << std::endl;


}
}  // namespace Kratos.
