//  KRATOS  _____     _ _ _
//         |_   _| __(_) (_)_ __   ___  ___
//           | || '__| | | | '_ \ / _ \/ __|
//           | || |  | | | | | | | (_) \__
//           |_||_|  |_|_|_|_| |_|\___/|___/ APPLICATION
//
//  License:             BSD License
//                                       Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

// System includes

// External includes
#include "mpi.h"

// Project includes
#include "trilinos_application.h"
#include "custom_factories/trilinos_linear_solver_factory.h"

namespace Kratos
{
void KratosTrilinosApplication::Register()
{
    // calling base class register to register Kratos components
    KratosApplication::Register();

    std::stringstream banner;
    banner << "    KRATOS  _____     _ _ _\n"
           << "           |_   _| __(_) (_)_ __   ___  ___\n"
           << "             | || '__| | | | '_ \\ / _ \\/ __|\n"
           << "             | || |  | | | | | | | (_) \\__ \\\n"
           << "             |_||_|  |_|_|_|_| |_|\\___/|___/\n"
           << "Initializing KratosTrilinosApplication..." << std::endl;

    int mpi_is_initialized = 0;
    int rank = -1;
    MPI_Initialized(&mpi_is_initialized);

    if (mpi_is_initialized)
    {
        MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    }

    if (mpi_is_initialized)
    {
        if (rank == 0) KRATOS_INFO("") << banner.str();
    }
    else
    {
        KRATOS_INFO("") << banner.str();
    }

    RegisterTrilinosLinearSolvers();
}


}  // namespace Kratos.


