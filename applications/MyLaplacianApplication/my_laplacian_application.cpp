//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand@{KRATOS_APP_AUTHOR}
//


// System includes


// External includes


// Project includes
#include "my_laplacian_application.h"
#include "my_laplacian_application_variables.h"


namespace Kratos {

KratosMyLaplacianApplication::KratosMyLaplacianApplication() {}

void KratosMyLaplacianApplication::Register() {
 	// calling base class register to register Kratos components
 	KratosApplication::Register();
 	std::cout << "Initializing KratosMyLaplacianApplication... " << std::endl;

  KRATOS_REGISTER_VARIABLE( MY_SCALAR )
  KRATOS_REGISTER_VARIABLE( MY_BOOL )
  KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( MY_VECTOR )

}
}  // namespace Kratos.
