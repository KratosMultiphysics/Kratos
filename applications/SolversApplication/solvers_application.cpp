//
//   Project Name:        KratosSolversApplication    $
//   Created by:          $Author:        JMCarbonell $
//   Last modified by:    $Co-Author:                 $
//   Date:                $Date:         January 2019 $
//   Revision:            $Revision:              0.0 $
//
//

// System includes

// External includes

// Project includes
#include "solvers_application.h"
#include "solvers_application_variables.h"


namespace Kratos {

KratosSolversApplication::KratosSolversApplication():
    KratosApplication("SolversApplication")
{}

void KratosSolversApplication::Register()
{
  // calling base class register to register Kratos components
  KratosApplication::Register();

  std::stringstream banner;
  banner << "Initializing KratosSolversApplication... " << std::endl;

  int mpi_is_initialized = 0;
  int rank = -1;

#ifdef KRATOS_MPI

  MPI_Initialized(&mpi_is_initialized);

  if (mpi_is_initialized)
  {
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  }

#endif

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
