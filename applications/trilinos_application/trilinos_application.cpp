//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:             BSD License 
//                                       Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//        

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "trilinos_application.h"
#include "Epetra_ConfigDefs.h"
#include "Epetra_DataAccess.h"

// #include "Epetra_MpiComm.h"

//#include "trilinos_solver.h"

namespace Kratos
{
// 	KRATOS_CREATE_VARIABLE( bool, IS_INACTIVE )

void KratosTrilinosApplication::Register()
{
    // calling base class register to register Kratos components
    KratosApplication::Register();
//     std::cout << "Initializing KratosTrilinosApplication... " << std::endl;
    std::cout << "     KRATOS   _____     _ _ _                 " << std::endl;
    std::cout << "             |_   _| __(_) (_)_ __   ___  ___ " << std::endl;
    std::cout << "               | || '__| | | | '_ \\ / _ \\/ __|" << std::endl;
    std::cout << "               | || |  | | | | | | | (_) \\__ \\" << std::endl;
    std::cout << "               |_||_|  |_|_|_|_| |_|\\___/|___/ APPLICATION     " << std::endl;
    // 		KRATOS_REGISTER_VARIABLE( IS_INACTIVE )
}


}  // namespace Kratos.


