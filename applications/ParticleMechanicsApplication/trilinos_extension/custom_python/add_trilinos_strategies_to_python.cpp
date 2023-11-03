//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Manuel Messmer
//

// Project includes
#include "custom_python/add_trilinos_strategies_to_python.h"

// Trilinos includes
#include "Epetra_FEVector.h"
#include "Epetra_FECrsMatrix.h"

// TrilinosApplication dependencies
#include "trilinos_space.h"

// KratosCore dependencies
#include "includes/model_part.h"
#include "linear_solvers/linear_solver.h"
#include "spaces/ublas_space.h"

// ParticleMechanicsApplication
#include "../custom_strategies/strategies/mpm_residual_based_newton_raphson_strategy.hpp"
#include "../custom_strategies/schemes/mpm_residual_based_bossak_scheme.hpp"
#include "../custom_builder_and_solvers/trilinos_mpm_block_builder_and_solver.h"

namespace Kratos::Python {

void AddTrilinosStrategiesToPython(pybind11::module& m){

    namespace py = pybind11;

    typedef TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector> TrilinosSparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> TrilinosLocalSpaceType;
    typedef LinearSolver<TrilinosSparseSpaceType, TrilinosLocalSpaceType > TrilinosLinearSolverType;
    typedef ConvergenceCriteria< TrilinosSparseSpaceType, TrilinosLocalSpaceType > TrilinosConvergenceCriteria;
    typedef SolvingStrategy< TrilinosSparseSpaceType, TrilinosLocalSpaceType > TrilinosBaseSolvingStrategyType;
    typedef Scheme< TrilinosSparseSpaceType, TrilinosLocalSpaceType > TrilinosBaseSchemeType;
    typedef BuilderAndSolver< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType > TrilinosBuilderAndSolverType;

    // Trilinos MPM Residual Based Strategy
    typedef MPMResidualBasedNewtonRaphsonStrategy< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType> TrilinosMPMResidualBasedNewtonRaphsonStrategyType;
    py::class_< TrilinosMPMResidualBasedNewtonRaphsonStrategyType,typename TrilinosMPMResidualBasedNewtonRaphsonStrategyType::Pointer, TrilinosBaseSolvingStrategyType >
        (m,"TrilinosMPMResidualBasedNewtonRaphsonStrategy")
        .def(py::init< ModelPart&, TrilinosBaseSchemeType::Pointer, TrilinosConvergenceCriteria::Pointer, TrilinosBuilderAndSolverType::Pointer, int, bool, bool, bool >() )
        ;

    // Trilinos MPM Residual Based Bossak Scheme Type
    typedef MPMResidualBasedBossakScheme<TrilinosSparseSpaceType, TrilinosLocalSpaceType> TrilinosMPMResidualBasedBossakSchemeType;
    py::class_< TrilinosMPMResidualBasedBossakSchemeType,typename TrilinosMPMResidualBasedBossakSchemeType::Pointer, TrilinosBaseSchemeType >
        (m,"TrilinosMPMResidualBasedBossakScheme")
        .def(py::init < ModelPart&, unsigned int, unsigned int, double, double, bool>())
        .def("Initialize", &TrilinosMPMResidualBasedBossakSchemeType::Initialize)
        ;

    // Trilinos MPM Block Builder and Solver
    typedef TrilinosBlockBuilderAndSolver< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType > TrilinosBlockBuilderAndSolverType;
    typedef TrilinosMPMBlockBuilderAndSolver< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType > TrilinosMPMBlockBuilderAndSolverType;
    py::class_< TrilinosMPMBlockBuilderAndSolverType, typename TrilinosMPMBlockBuilderAndSolverType::Pointer, TrilinosBlockBuilderAndSolverType  >
        (m, "TrilinosMPMBlockBuilderAndSolver")
        .def(py::init<Epetra_MpiComm&, int, TrilinosLinearSolverType::Pointer > () )
        ;

}

} // namespace Kratos::Python
