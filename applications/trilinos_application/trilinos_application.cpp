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
#include "mpi.h"
// #include "Epetra_MpiComm.h"

namespace Kratos
{
// 	KRATOS_CREATE_VARIABLE( bool, IS_INACTIVE )

void KratosTrilinosApplication::Register()
{
    // calling base class register to register Kratos components
    KratosApplication::Register();

    std::stringstream banner;
    banner << "     KRATOS   _____     _ _ _                 " << std::endl;
    banner << "             |_   _| __(_) (_)_ __   ___  ___ " << std::endl;
    banner << "               | || '__| | | | '_ \\ / _ \\/ __|" << std::endl;
    banner << "               | || |  | | | | | | | (_) \\__ \\" << std::endl;
    banner << "               |_||_|  |_|_|_|_| |_|\\___/|___/ APPLICATION     " << std::endl;

    int mpi_is_initialized = 0;
    int rank = -1;
    MPI_Initialized(&mpi_is_initialized);

    if (mpi_is_initialized)
    {
        MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    }

    if (mpi_is_initialized)
    {
        if (rank == 0) std::cout << banner.str();
    }
    else
    {
        std::cout << banner.str();
    }
}


}  // namespace Kratos.


