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
#include "factories/builder_and_solver_factory.h"
#include "spaces/ublas_space.h"

namespace Kratos 
{
    namespace Testing 
    {
        /// Tests
       
        typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
        typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
        typedef LinearSolver<SparseSpaceType, LocalSpaceType> LinearSolverType;
        
        /// The definition of the p_builder_ans_solver
        typedef BuilderAndSolver<SparseSpaceType,LocalSpaceType, LinearSolverType> BuilderAndSolverType;

        /// The definition of the factory
        typedef BuilderAndSolverFactory<SparseSpaceType,LocalSpaceType, LinearSolverType> FactoryType;
     
        /**
         * Checks if the ResidualBasedEliminationBuilderAndSolver performs correctly the Factory
         */
        KRATOS_TEST_CASE_IN_SUITE(ResidualBasedEliminationBuilderAndSolverFactory, KratosCoreFastSuite)
        {
            LinearSolverType::Pointer p_solver = nullptr;
            Parameters this_parameters = Parameters(R"({"builder_and_solver_type" : "ResidualBasedEliminationBuilderAndSolver"})");
            BuilderAndSolverType::Pointer p_builder_ans_solver = FactoryType().Create(p_solver, this_parameters);
            KRATOS_CHECK_C_STRING_EQUAL((p_builder_ans_solver->Info()).c_str(), "ResidualBasedEliminationBuilderAndSolver");
        }

        /**
         * Checks if the ResidualBasedBlockBuilderAndSolver performs correctly the Factory
         */
        KRATOS_TEST_CASE_IN_SUITE(ResidualBasedBlockBuilderAndSolverFactory, KratosCoreFastSuite)
        {
            LinearSolverType::Pointer p_solver = nullptr;
            Parameters this_parameters = Parameters(R"({"builder_and_solver_type" : "ResidualBasedBlockBuilderAndSolver"})");
            BuilderAndSolverType::Pointer p_builder_ans_solver = FactoryType().Create(p_solver, this_parameters);
            KRATOS_CHECK_C_STRING_EQUAL((p_builder_ans_solver->Info()).c_str(), "ResidualBasedBlockBuilderAndSolver");
        }

        /**
         * Checks if the ResidualBasedBlockBuilderAndSolverWithConstraints performs correctly the Factory
         */
        KRATOS_TEST_CASE_IN_SUITE(ResidualBasedBlockBuilderAndSolverWithConstraintsFactory, KratosCoreFastSuite)
        {
            LinearSolverType::Pointer p_solver = nullptr;
            Parameters this_parameters = Parameters(R"({"builder_and_solver_type" : "ResidualBasedBlockBuilderAndSolverWithConstraints"})");
            BuilderAndSolverType::Pointer p_builder_ans_solver = FactoryType().Create(p_solver, this_parameters);
            KRATOS_CHECK_C_STRING_EQUAL((p_builder_ans_solver->Info()).c_str(), "ResidualBasedBlockBuilderAndSolverWithConstraints");
        }
        
    } // namespace Testing
}  // namespace Kratos.

