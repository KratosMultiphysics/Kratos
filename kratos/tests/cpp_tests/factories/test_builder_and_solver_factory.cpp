//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes
#include <limits>

// External includes

// Project includes
#include "testing/testing.h"
#include "containers/model.h"

// Utility includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "factories/factory.h"
#include "spaces/ublas_space.h"

namespace Kratos
{
    namespace Testing
    {
        /// Tests

        // Spaces
        using SparseSpaceType = UblasSpace<double, CompressedMatrix, Vector>;
        using LocalSpaceType = UblasSpace<double, Matrix, Vector>;
        using LinearSolverType = LinearSolver<SparseSpaceType, LocalSpaceType>;

        /// The definition of the p_builder_and_solver
        using BuilderAndSolverType = BuilderAndSolver<SparseSpaceType,LocalSpaceType, LinearSolverType>;

        /// The definition of the factory
        using BuilderAndSolverFactoryType = Factory<BuilderAndSolverType, LinearSolverType>;

        /**
         * Checks if the ResidualBasedEliminationBuilderAndSolver performs correctly the Factory
         */
        KRATOS_TEST_CASE_IN_SUITE(ResidualBasedEliminationBuilderAndSolverFactory, KratosCoreFastSuite)
        {
            LinearSolverType::Pointer p_solver = nullptr;
            Parameters this_parameters = Parameters(R"({"name" : "elimination_builder_and_solver"})");
            BuilderAndSolverType::Pointer p_builder_and_solver = BuilderAndSolverFactoryType().Create(p_solver, this_parameters);
            KRATOS_CHECK_STRING_EQUAL(p_builder_and_solver->Info(), "ResidualBasedEliminationBuilderAndSolver");
        }

        /**
         * Checks if the ResidualBasedEliminationBuilderAndSolverWithConstraints performs correctly the Factory
         */
        KRATOS_TEST_CASE_IN_SUITE(ResidualBasedEliminationBuilderAndSolverWithConstraintsFactory, KratosCoreFastSuite)
        {
            LinearSolverType::Pointer p_solver = nullptr;
            Parameters this_parameters = Parameters(R"({"name" : "elimination_builder_and_solver_with_constraints"})");
            BuilderAndSolverType::Pointer p_builder_and_solver = BuilderAndSolverFactoryType().Create(p_solver, this_parameters);
            KRATOS_CHECK_STRING_EQUAL(p_builder_and_solver->Info(), "ResidualBasedEliminationBuilderAndSolverWithConstraints");
        }

        /**
         * Checks if the ResidualBasedBlockBuilderAndSolver performs correctly the Factory
         */
        KRATOS_TEST_CASE_IN_SUITE(ResidualBasedBlockBuilderAndSolverFactory, KratosCoreFastSuite)
        {
            LinearSolverType::Pointer p_solver = nullptr;
            Parameters this_parameters = Parameters(R"({"name" : "block_builder_and_solver"})");
            BuilderAndSolverType::Pointer p_builder_and_solver = BuilderAndSolverFactoryType().Create(p_solver, this_parameters);
            KRATOS_CHECK_STRING_EQUAL(p_builder_and_solver->Info(), "ResidualBasedBlockBuilderAndSolver");
        }

    } // namespace Testing
}  // namespace Kratos.

