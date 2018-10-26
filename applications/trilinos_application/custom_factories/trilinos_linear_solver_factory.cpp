//  KRATOS  _____     _ _ _
//         |_   _| __(_) (_)_ __   ___  ___
//           | || '__| | | | '_ \ / _ \/ __|
//           | || |  | | | | | | | (_) \__
//           |_||_|  |_|_|_|_| |_|\___/|___/ APPLICATION
//
//  License:             BSD License
//                                       Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes
#include "mpi.h"

// External includes
/* Trilinos includes */
#include "Epetra_MpiComm.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_FECrsGraph.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_SerialDenseMatrix.h"

// Teuchos parameter list
#include "Teuchos_ParameterList.hpp"

#include "external_includes/epetra_default_utility.h"
#include "external_includes/aztec_solver.h"
#include "external_includes/amesos_solver.h"
#include "external_includes/ml_solver.h"

#include "external_includes/amgcl_mpi_solver.h"
//#include "external_includes/amgcl_deflation_solver.h"
#include "external_includes/amgcl_mpi_schur_complement_solver.h"

// Project includes
#include "trilinos_application.h"
#include "spaces/ublas_space.h"
#include "includes/model_part.h"
#include "custom_factories/trilinos_linear_solver_factory.h"

// Linear solvers
#include "linear_solvers/linear_solver.h"

namespace Kratos
{
    void RegisterTrilinosLinearSolvers()
    {
        typedef TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector> TrilinosSparseSpaceType;
        typedef UblasSpace<double, Matrix, Vector> TrilinosLocalSpaceType;
        typedef LinearSolver<TrilinosSparseSpaceType, TrilinosLocalSpaceType > TrilinosLinearSolverType;

        typedef AztecSolver<TrilinosSparseSpaceType, TrilinosLocalSpaceType > AztecSolverType;
        typedef AmesosSolver<TrilinosSparseSpaceType, TrilinosLocalSpaceType > AmesosSolverType;
        typedef MultiLevelSolver<TrilinosSparseSpaceType, TrilinosLocalSpaceType > MLSolverType;
        typedef AmgclMPISolver<TrilinosSparseSpaceType, TrilinosLocalSpaceType > AmgclMPISolverType;
        typedef AmgclMPISchurComplementSolver<TrilinosSparseSpaceType, TrilinosLocalSpaceType > AmgclMPISchurComplementSolverType;

        //NOTE: here we must create persisting objects for the linear solvers
        static auto AztecSolverFactory = TrilinosLinearSolverFactory<TrilinosSparseSpaceType,TrilinosLocalSpaceType,AztecSolverType>();
        static auto AmesosSolverFactory = TrilinosLinearSolverFactory<TrilinosSparseSpaceType,TrilinosLocalSpaceType,AmesosSolverType>();
        static auto MultiLevelSolverFactory = TrilinosLinearSolverFactory<TrilinosSparseSpaceType,TrilinosLocalSpaceType,MLSolverType>();
        static auto AmgclMPISolverFactory = TrilinosLinearSolverFactory<TrilinosSparseSpaceType,TrilinosLocalSpaceType,AmgclMPISolverType>();
        static auto AmgclMPISchurComplementSolverFactory = TrilinosLinearSolverFactory<TrilinosSparseSpaceType,TrilinosLocalSpaceType,AmgclMPISchurComplementSolverType>();

        // Registration of linear solvers
        KRATOS_REGISTER_TRILINOS_LINEAR_SOLVER("AztecSolver", AztecSolverFactory);
        KRATOS_REGISTER_TRILINOS_LINEAR_SOLVER("AmesosSolver", AmesosSolverFactory);
        KRATOS_REGISTER_TRILINOS_LINEAR_SOLVER("MultiLevelSolver", MultiLevelSolverFactory);
        KRATOS_REGISTER_TRILINOS_LINEAR_SOLVER("AmgclMPISolver", AmgclMPISolverFactory);
        KRATOS_REGISTER_TRILINOS_LINEAR_SOLVER("AmgclMPISchurComplementSolver", AmgclMPISchurComplementSolverFactory);

    };
} // Namespace Kratos

