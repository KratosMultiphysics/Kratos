//
//   Project Name:        KratosSolversApplication $
//   Created by:          $Author:     JMCarbonell $
//   Last modified by:    $Co-Author:              $
//   Date:                $Date:      January 2019 $
//   Revision:            $Revision:           0.0 $
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

  banner << "            ___      _                            \n"
         << "    KRATOS / __| ___| |_ _ ___ _ _ ___            \n"
         << "           \\__ \\/ _ \\ \\ V / -_) '_|_-<            \n"
         << "           |___/\\___/_|\\_/\\___|_| /__/ APPLICATION\n"
         << "Initialize KratosSolversApplication..." << std::endl;

  // mpi initialization
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
    if (rank == 0) KRATOS_INFO("") << banner.str();
  }
  else
  {
    KRATOS_INFO("") << banner.str();
  }

  // Register Variables (variables created in solvers_application_variables.cpp)

  // time settings
  KRATOS_REGISTER_VARIABLE(MESHING_STEP_TIME)
  KRATOS_REGISTER_VARIABLE(CONTACT_STEP_TIME)
  KRATOS_REGISTER_VARIABLE(RESTART_STEP_TIME)

  // time integration methods
  KRATOS_REGISTER_VARIABLE(VECTOR_TIME_INTEGRATION_METHODS)
  KRATOS_REGISTER_VARIABLE(SCALAR_TIME_INTEGRATION_METHODS)
  KRATOS_REGISTER_VARIABLE(COMPONENT_TIME_INTEGRATION_METHODS)

  // implicit solution
  KRATOS_REGISTER_VARIABLE(CONVERGENCE_ACHIEVED)
  KRATOS_REGISTER_VARIABLE(COMPUTE_CONSISTENT_MASS_MATRIX)

  KRATOS_REGISTER_VARIABLE(SEGREGATED_STEP)
  KRATOS_REGISTER_VARIABLE(TIME_INTEGRATION_ORDER)

  KRATOS_REGISTER_VARIABLE(RAYLEIGH_ALPHA)
  KRATOS_REGISTER_VARIABLE(RAYLEIGH_BETA)

  // explicit solution
  KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(MIDDLE_VELOCITY)

  KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(EXTERNAL_MOMENT)
  KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(POSITION_MOMENTUM)
  KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(ROTATION_MOMENTUM)
  KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(RESIDUAL_LYAPUNOV)

  KRATOS_REGISTER_VARIABLE(INERTIA_DYADIC)
  KRATOS_REGISTER_VARIABLE(TANGENT_MATRIX)
  KRATOS_REGISTER_VARIABLE(TANGENT_LYAPUNOV)

  KRATOS_REGISTER_VARIABLE(ALPHA_TRAPEZOIDAL_RULE)
  KRATOS_REGISTER_VARIABLE(POSITION_UPDATE_LABEL)
  KRATOS_REGISTER_VARIABLE(ROTATION_UPDATE_LABEL)
  KRATOS_REGISTER_VARIABLE(MOMENTUM_UPDATE_LABEL)

  // eigenvalue solution
  KRATOS_REGISTER_VARIABLE(BUILD_LEVEL)
  KRATOS_REGISTER_VARIABLE(EIGENVALUE_VECTOR)
  KRATOS_REGISTER_VARIABLE(EIGENVECTOR_MATRIX)

  //Register Solver Factories

#ifdef INCLUDE_SUPERLU_MT
  KRATOS_REGISTER_LINEAR_SOLVER("superlu_direct", mSuperLUmtDirectSolverFactory);
#else
  KRATOS_REGISTER_LINEAR_SOLVER("superlu_direct", mSuperLUDirectSolverFactory);
  //KRATOS_REGISTER_LINEAR_SOLVER("super_lu_iterative", mSuperLUIterativeSolverFactory);
#endif

#ifdef INCLUDE_FEAST
  KRATOS_REGISTER_LINEAR_SOLVER("feast_eigen", mFEASTEigenValueSolverFactory);
#endif

}
}  // namespace Kratos.
