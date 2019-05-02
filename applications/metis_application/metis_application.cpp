//
//   Project Name:        Kratos
//   Last Modified by:    $Author: rrossi $
//   Date:                $Date: 2008-07-09 15:09:17 $
//   Revision:            $Revision: 1.2 $
//
//

// System includes

// External includes
#include "mpi.h"

// Project includes
#include "metis_application.h"

namespace Kratos {

KratosMetisApplication::KratosMetisApplication()
    : KratosApplication("MetisApplication") {}

void KratosMetisApplication::Register()
{
    // calling base class register to register Kratos components
    KratosApplication::Register();

    std::stringstream banner; // TODO: use Logger once mpi-logger is implemented
    banner << "    KRATOS  __  __      _   _\n"
           << "           |  \\/  | ___| |_(_)___\n"
           << "           | |\\/| |/ _ \\ __| / __|\n"
           << "           | |  | |  __/ |_| \\__ \\\n"
           << "           |_|  |_|\\___|\\__|_|___/\n"
           << "Initializing KratosMetisApplication..." << std::endl;

    int mpi_is_initialized = 0;
    int rank = -1;
    MPI_Initialized(&mpi_is_initialized);

    if (mpi_is_initialized){
        MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    }

    if (mpi_is_initialized) {
        if (rank == 0) KRATOS_INFO("") << banner.str();
    }
    else {
        KRATOS_INFO("") << banner.str();
    }
}

}  // namespace Kratos.
