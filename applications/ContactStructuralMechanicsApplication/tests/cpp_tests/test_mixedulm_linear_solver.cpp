// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:             BSD License
//                                       license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes
#include <limits>

/* External includes */

/* Project includes */
#include "testing/testing.h"
#include "containers/model.h"

/* Utility includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "spaces/ublas_space.h"
#include "includes/matrix_market_interface.h"
#include "includes/kratos_filesystem.h"

// Linear solvers
#include "linear_solvers/reorderer.h"
#include "linear_solvers/preconditioner.h"
#include "linear_solvers/direct_solver.h"
#include "linear_solvers/linear_solver.h"
#include "linear_solvers/skyline_lu_factorization_solver.h"
#include "linear_solvers/amgcl_solver.h"
#include "custom_linear_solvers/mixedulm_linear_solver.h"

namespace Kratos
{
    namespace Testing
    {
        /// Tests
        typedef Node<3> NodeType;
        typedef Geometry<NodeType> GeometryType;
        typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
        typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;

        // The direct solver
        typedef Reorderer<SparseSpaceType,  LocalSpaceType > ReordererType;
        typedef DirectSolver<SparseSpaceType,  LocalSpaceType, ReordererType > DirectSolverType;
        typedef LinearSolver<SparseSpaceType,LocalSpaceType> LinearSolverType;
        typedef AMGCLSolver<SparseSpaceType,  LocalSpaceType, ReordererType > AMGCLSolverType;
        typedef SkylineLUFactorizationSolver<SparseSpaceType,  LocalSpaceType, ReordererType > SkylineLUFactorizationSolverType;
        typedef Preconditioner<SparseSpaceType, LocalSpaceType> PreconditionerType;
        typedef MixedULMLinearSolver<SparseSpaceType,  LocalSpaceType, PreconditionerType, ReordererType> MixedULMLinearSolverType;

        // Dof arrays
        typedef ModelPart::DofsArrayType DofsArrayType;

        /**
         * @brief This creates the auxiliar files
         */
        void CreateAuxiliarFiles()
        {
            const std::string lhs = "%%MatrixMarket matrix coordinate real general\n16 16 184\n1 1 4.00969e+10\n1 2 0\n1 3 0\n1 4 0\n1 7 0\n1 8 0\n2 1 0\n2 2 1.35566e+11\n2 3 0\n2 4 0\n2 7 0\n2 8 0\n3 1 0\n3 2 0\n3 3 1.4559e+11\n3 4 0\n3 5 0\n3 6 0\n3 7 0\n3 8 0\n3 9 0\n3 10 0\n3 11 0\n3 12 0\n3 13 0\n3 14 0\n4 1 0\n4 2 0\n4 3 0\n4 4 2.22633e+10\n4 5 0\n4 6 -1.0345e+11\n4 7 0\n4 8 -6.20547e+09\n4 9 0\n4 10 -1.14853e-05\n4 11 0\n4 12 0.000525451\n4 13 0\n4 14 5.1725e+10\n5 3 0\n5 4 0\n5 5 2.069e+11\n5 6 0\n5 7 0\n5 8 0\n5 9 0\n5 10 0\n5 11 0\n5 12 0\n5 13 0\n5 14 0\n6 3 0\n6 4 -1.0345e+11\n6 5 0\n6 6 0\n6 7 0\n6 8 0\n6 9 0\n6 10 0\n6 11 0\n6 12 0.00103942\n6 13 0\n6 14 1.0345e+11\n7 1 0\n7 2 0\n7 3 0\n7 4 0\n7 5 0\n7 6 0\n7 7 1.4559e+11\n7 8 0\n7 9 0\n7 10 0\n7 11 0\n7 12 0\n7 13 0\n7 14 0\n8 1 0\n8 2 0\n8 3 0\n8 4 -6.20547e+09\n8 5 0\n8 6 0\n8 7 0\n8 8 2.22633e+10\n8 9 0\n8 10 -1.0345e+11\n8 11 0\n8 12 5.1725e+10\n8 13 0\n8 14 5.74263e-06\n9 3 0\n9 4 0\n9 5 0\n9 6 0\n9 7 0\n9 8 0\n9 9 2.069e+11\n9 10 0\n9 11 0\n9 12 0\n9 13 0\n9 14 0\n10 3 0\n10 4 -1.14853e-05\n10 5 0\n10 6 0\n10 7 0\n10 8 -1.0345e+11\n10 9 0\n10 10 0\n10 11 0\n10 12 1.0345e+11\n10 13 0\n10 14 1.14853e-05\n11 3 0\n11 4 0\n11 5 0\n11 6 0\n11 7 0\n11 8 0\n11 9 0\n11 10 0\n11 11 1.45464e+11\n11 12 0\n11 13 0\n11 14 0\n11 15 0\n11 16 0\n12 3 0\n12 4 0.000525451\n12 5 0\n12 6 0.00103942\n12 7 0\n12 8 5.1725e+10\n12 9 0\n12 10 1.0345e+11\n12 11 0\n12 12 2.22572e+10\n12 13 0\n12 14 -6.13145e+09\n12 15 0\n12 16 0\n13 3 0\n13 4 0\n13 5 0\n13 6 0\n13 7 0\n13 8 0\n13 9 0\n13 10 0\n13 11 0\n13 12 0\n13 13 1.45464e+11\n13 14 0\n13 15 0\n13 16 0\n14 3 0\n14 4 5.1725e+10\n14 5 0\n14 6 1.0345e+11\n14 7 0\n14 8 5.74263e-06\n14 9 0\n14 10 1.14853e-05\n14 11 0\n14 12 -6.13145e+09\n14 13 0\n14 14 2.22572e+10\n14 15 0\n14 16 0\n15 11 0\n15 12 0\n15 13 0\n15 14 0\n15 15 4.0137e+10\n15 16 0\n16 11 0\n16 12 0\n16 13 0\n16 14 0\n16 15 0\n16 16 1.35701e+11\n";

            std::filebuf buffer_lhs;
            buffer_lhs.open(FilesystemExtensions::JoinPaths({FilesystemExtensions::CurrentWorkingDirectory(), "A_testing_condensation.mm"}),std::ios::out);
            std::ostream os_lhs(&buffer_lhs);
            os_lhs << lhs;
            buffer_lhs.close();

            const std::string rhs = "%%MatrixMarket matrix array real general\n16 1\n0\n0\n0\n-5.1725e+07\n0\n-1.0345e+08\n0\n-5.1725e+07\n0\n-1.0345e+08\n0\n-6.23386e+08\n0\n-6.23386e+08\n0\n0\n";

            std::filebuf buffer_rhs;
            buffer_rhs.open(FilesystemExtensions::JoinPaths({FilesystemExtensions::CurrentWorkingDirectory(), "b_testing_condensation.rhs"}),std::ios::out);
            std::ostream os_rhs(&buffer_rhs);
            os_rhs << rhs;
            buffer_rhs.close();
        }

        /**
         * Checks if the MixedULMLinear solver performs correctly the resolution of the system
         */
        KRATOS_TEST_CASE_IN_SUITE(MixedULMLinearSolverSimplestSystem, KratosContactStructuralMechanicsFastSuite)
        {
            constexpr double tolerance = 1e-6;

            Model this_model;
            ModelPart& r_model_part = this_model.CreateModelPart("Main", 3);

            LinearSolverType::Pointer psolver = LinearSolverType::Pointer( new SkylineLUFactorizationSolverType() );
//             Parameters empty_parameters =  Parameters(R"({})");
//             LinearSolverType::Pointer psolver = LinearSolverType::Pointer( new AMGCLSolverType(empty_parameters) );
            LinearSolverType::Pointer pmixed_solver = LinearSolverType::Pointer( new MixedULMLinearSolverType(psolver) );

            r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
            r_model_part.AddNodalSolutionStepVariable(VECTOR_LAGRANGE_MULTIPLIER);

            NodeType::Pointer pnode1 = r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
            NodeType::Pointer pnode2 = r_model_part.CreateNewNode(2, 0.0, 0.0, 0.0);
            pnode2->Set(INTERFACE, true);
            pnode2->Set(MASTER, true);
            pnode2->Set(SLAVE, false);
            NodeType::Pointer pnode3 = r_model_part.CreateNewNode(3, 0.0, 0.0, 0.0);
            pnode3->Set(INTERFACE, true);
            pnode3->Set(ACTIVE, true);
            pnode3->Set(MASTER, false);
            pnode3->Set(SLAVE, true);

            pnode1->AddDof(DISPLACEMENT_X);
            pnode2->AddDof(DISPLACEMENT_X);
            pnode3->AddDof(DISPLACEMENT_X);
            pnode3->AddDof(VECTOR_LAGRANGE_MULTIPLIER_X);

            std::vector< Dof<double>::Pointer > DoF;
            DoF.reserve(4);
            DoF.push_back(pnode1->pGetDof(DISPLACEMENT_X));
            DoF.push_back(pnode2->pGetDof(DISPLACEMENT_X));
            DoF.push_back(pnode3->pGetDof(DISPLACEMENT_X));
            DoF.push_back(pnode3->pGetDof(VECTOR_LAGRANGE_MULTIPLIER_X));

            // Set initial solution
            (pnode1->FastGetSolutionStepValue(DISPLACEMENT)).clear();
            (pnode2->FastGetSolutionStepValue(DISPLACEMENT)).clear();
            (pnode3->FastGetSolutionStepValue(DISPLACEMENT)).clear();
            (pnode3->FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER)).clear();

            DofsArrayType Doftemp;
            Doftemp.reserve(DoF.size());
            for (auto it= DoF.begin(); it!= DoF.end(); it++)
                Doftemp.push_back( *it );

            const std::size_t system_size = 4;
            CompressedMatrix A(system_size, system_size);
            Vector ref_Dx = ZeroVector(system_size);
            Vector Dx = ZeroVector(system_size);
            Vector b = ZeroVector(system_size);
            double count = 0.0;
            for (std::size_t i = 0; i < system_size; ++i) {
                for (std::size_t j = 0; j < system_size; ++j) {
                    if (((i == 0 && j == system_size - 1) || (j == 0 && i == system_size - 1)) == false) {
                        count += 1.0;
                        A.push_back(i, j, std::sqrt(count));
                    }
                }
            }
            count = 0.0;
            for (std::size_t i = 0; i < system_size; ++i) {
                count += 1.0;
                b[i] = count;
            }

            // We solve the reference system
            psolver->Solve(A, ref_Dx, b);

            // We solve the block system
            pmixed_solver->ProvideAdditionalData(A, Dx, b, Doftemp, r_model_part);
            pmixed_solver->Solve(A, Dx, b);

            for (std::size_t i = 0; i < system_size; ++i) {
                KRATOS_CHECK_NEAR(ref_Dx[i], Dx[i], tolerance);
            }
        }

        /**
         * Checks if the MixedULMLinear solver performs correctly the resolution of the system with inactive nodes
         */
        KRATOS_TEST_CASE_IN_SUITE(MixedULMLinearSolverSimplestWithInactiveSystem, KratosContactStructuralMechanicsFastSuite)
        {
            constexpr double tolerance = 1e-6;

            Model this_model;
            ModelPart& r_model_part = this_model.CreateModelPart("Main", 3);

            LinearSolverType::Pointer psolver = LinearSolverType::Pointer( new SkylineLUFactorizationSolverType() );
//             Parameters empty_parameters =  Parameters(R"({})");
//             LinearSolverType::Pointer psolver = LinearSolverType::Pointer( new AMGCLSolverType(empty_parameters) );
            LinearSolverType::Pointer pmixed_solver = LinearSolverType::Pointer( new MixedULMLinearSolverType(psolver) );

            r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
            r_model_part.AddNodalSolutionStepVariable(VECTOR_LAGRANGE_MULTIPLIER);

            NodeType::Pointer pnode1 = r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
            NodeType::Pointer pnode2 = r_model_part.CreateNewNode(2, 0.0, 0.0, 0.0);

            pnode2->Set(INTERFACE, true);
            pnode2->Set(MASTER, true);
            pnode2->Set(SLAVE, false);
            NodeType::Pointer pnode3 = r_model_part.CreateNewNode(3, 0.0, 0.0, 0.0);
            pnode3->Set(INTERFACE, true);
            pnode3->Set(ACTIVE, false);
            pnode3->Set(MASTER, false);
            pnode3->Set(SLAVE, true);
            NodeType::Pointer pnode4 = r_model_part.CreateNewNode(4, 0.0, 0.0, 0.0);
            pnode4->Set(INTERFACE, true);
            pnode4->Set(ACTIVE, true);
            pnode4->Set(MASTER, false);
            pnode4->Set(SLAVE, true);

            pnode1->AddDof(DISPLACEMENT_X);
            pnode2->AddDof(DISPLACEMENT_X);
            pnode3->AddDof(DISPLACEMENT_X);
            pnode4->AddDof(DISPLACEMENT_X);
            pnode4->AddDof(VECTOR_LAGRANGE_MULTIPLIER_X);

            std::vector< Dof<double>::Pointer > DoF;
            DoF.reserve(5);
            DoF.push_back(pnode1->pGetDof(DISPLACEMENT_X));
            DoF.push_back(pnode2->pGetDof(DISPLACEMENT_X));
            DoF.push_back(pnode3->pGetDof(DISPLACEMENT_X));
            DoF.push_back(pnode4->pGetDof(DISPLACEMENT_X));
            DoF.push_back(pnode4->pGetDof(VECTOR_LAGRANGE_MULTIPLIER_X));

            // Set initial solution
            (pnode1->FastGetSolutionStepValue(DISPLACEMENT)).clear();
            (pnode2->FastGetSolutionStepValue(DISPLACEMENT)).clear();
            (pnode4->FastGetSolutionStepValue(DISPLACEMENT)).clear();
            (pnode4->FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER)).clear();

            DofsArrayType Doftemp;
            Doftemp.reserve(DoF.size());
            for (auto it= DoF.begin(); it!= DoF.end(); it++)
                Doftemp.push_back( *it );

            const std::size_t system_size = 5;
            CompressedMatrix A(system_size, system_size);
            Vector ref_Dx = ZeroVector(system_size);
            Vector Dx = ZeroVector(system_size);
            Vector b = ZeroVector(system_size);
            double count = 0.0;
            for (std::size_t i = 0; i < system_size; ++i) {
                for (std::size_t j = 0; j < system_size; ++j) {
                    if (((i ==2 && j == system_size - 1) || (i == 0 && j == system_size - 1) || (j == 0 && i == system_size - 1)) == false) {
                        count += 1.0;
                        A.push_back(i, j, std::sqrt(count));
                    }
                }
            }
            count = 0.0;
            for (std::size_t i = 0; i < system_size; ++i) {
                count += 1.0;
                b[i] = count;
            }

            // We solve the reference system
            psolver->Solve(A, ref_Dx, b);

            // We solve the block system
            pmixed_solver->ProvideAdditionalData(A, Dx, b, Doftemp, r_model_part);
            pmixed_solver->Solve(A, Dx, b);

            for (std::size_t i = 0; i < system_size; ++i) {
                KRATOS_CHECK_NEAR(ref_Dx[i], Dx[i], tolerance);
            }
        }

        /**
         * Checks if the MixedULMLinear solver performs correctly the resolution of the system. Unordered case
         */
        KRATOS_TEST_CASE_IN_SUITE(MixedULMLinearSolverSimplestUnorderedSystem, KratosContactStructuralMechanicsFastSuite)
        {
            constexpr double tolerance = 1e-6;

            Model this_model;
            ModelPart& r_model_part = this_model.CreateModelPart("Main", 3);

            LinearSolverType::Pointer psolver = LinearSolverType::Pointer( new SkylineLUFactorizationSolverType() );
//             Parameters empty_parameters =  Parameters(R"({})");
//             LinearSolverType::Pointer psolver = LinearSolverType::Pointer( new AMGCLSolverType(empty_parameters) );
            LinearSolverType::Pointer pmixed_solver = LinearSolverType::Pointer( new MixedULMLinearSolverType(psolver) );

            r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
            r_model_part.AddNodalSolutionStepVariable(VECTOR_LAGRANGE_MULTIPLIER);

            NodeType::Pointer pnode1 = r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
            NodeType::Pointer pnode3 = r_model_part.CreateNewNode(3, 0.0, 0.0, 0.0);
            pnode3->Set(INTERFACE, true);
            pnode3->Set(ACTIVE, true);
            pnode3->Set(MASTER, false);
            pnode3->Set(SLAVE, true);
            NodeType::Pointer pnode2 = r_model_part.CreateNewNode(2, 0.0, 0.0, 0.0);
            pnode2->Set(INTERFACE, true);
            pnode2->Set(MASTER, true);
            pnode2->Set(SLAVE, false);

            pnode1->AddDof(DISPLACEMENT_X);
            pnode3->AddDof(DISPLACEMENT_X);
            pnode3->AddDof(VECTOR_LAGRANGE_MULTIPLIER_X);
            pnode2->AddDof(DISPLACEMENT_X);

            std::vector< Dof<double>::Pointer > DoF;
            DoF.reserve(4);
            DoF.push_back(pnode1->pGetDof(DISPLACEMENT_X));
            DoF.push_back(pnode3->pGetDof(DISPLACEMENT_X));
            DoF.push_back(pnode3->pGetDof(VECTOR_LAGRANGE_MULTIPLIER_X));
            DoF.push_back(pnode2->pGetDof(DISPLACEMENT_X));

            // Set initial solution
            (pnode1->FastGetSolutionStepValue(DISPLACEMENT)).clear();
            (pnode3->FastGetSolutionStepValue(DISPLACEMENT)).clear();
            (pnode3->FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER)).clear();
            (pnode2->FastGetSolutionStepValue(DISPLACEMENT)).clear();

            DofsArrayType Doftemp;
            Doftemp.reserve(DoF.size());
            for (auto it= DoF.begin(); it!= DoF.end(); it++)
                Doftemp.push_back( *it );

            const std::size_t system_size = 4;
            CompressedMatrix A(system_size, system_size);
            CompressedMatrix Aaux(system_size, system_size);
            Vector ref_Dx = ZeroVector(system_size);
            Vector Dx = ZeroVector(system_size);
            Vector b = ZeroVector(system_size);
            Vector baux = ZeroVector(system_size);
            double count = 0.0;
            for (std::size_t i = 0; i < system_size; ++i) {
                for (std::size_t j = 0; j < system_size; ++j) {
                    if (((i == 0 && j == system_size - 1) || (j == 0 && i == system_size - 1)) == false) {
                        count += 1.0;
                        Aaux.push_back(i, j, std::sqrt(count));
                    }
                }
            }

            for (std::size_t i = 0; i < system_size; ++i) {
                std::size_t iaux = i;
                if (i == 1)
                    iaux = 2;
                else if (i == 2)
                    iaux = 3;
                else if (i == 3)
                    iaux = 1;
                for (std::size_t j = 0; j < system_size; ++j) {
                    std::size_t jaux = j;
                    if (j == 1)
                        jaux = 2;
                    else if (j == 2)
                        jaux = 3;
                    else if (j == 3)
                        jaux = 1;
                    A.push_back(i, j, Aaux(iaux, jaux));
                }
            }

            count = 0.0;
            for (std::size_t i = 0; i < system_size; ++i) {
                count += 1.0;
                baux[i] = count;
            }

            b[0] = baux[0];
            b[1] = baux[2];
            b[2] = baux[3];
            b[3] = baux[1];

            // We solve the reference system
            psolver->Solve(A, ref_Dx, b);

            // We solve the block system
            pmixed_solver->ProvideAdditionalData(A, Dx, b, Doftemp, r_model_part);
            pmixed_solver->Solve(A, Dx, b);

            for (std::size_t i = 0; i < system_size; ++i) {
                KRATOS_CHECK_NEAR(ref_Dx[i], Dx[i], tolerance);
            }
        }

        /**
         * Checks if the MixedULMLinear solver performs correctly the resolution of the system. Multiple dofs
         */

        KRATOS_TEST_CASE_IN_SUITE(MixedULMLinearSolverTwoDoFSystem, KratosContactStructuralMechanicsFastSuite)
        {
            constexpr double tolerance = 1e-5;

            Model this_model;
            ModelPart& r_model_part = this_model.CreateModelPart("Main", 3);

            LinearSolverType::Pointer psolver = LinearSolverType::Pointer( new SkylineLUFactorizationSolverType() );
//             Parameters empty_parameters =  Parameters(R"({})");
//             LinearSolverType::Pointer psolver = LinearSolverType::Pointer( new AMGCLSolverType(empty_parameters) );
            LinearSolverType::Pointer pmixed_solver = LinearSolverType::Pointer( new MixedULMLinearSolverType(psolver) );

            r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
            r_model_part.AddNodalSolutionStepVariable(VECTOR_LAGRANGE_MULTIPLIER);

            NodeType::Pointer pnode1 = r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
            NodeType::Pointer pnode2 = r_model_part.CreateNewNode(2, 0.0, 0.0, 0.0);
            pnode2->Set(INTERFACE, true);
            pnode2->Set(MASTER, true);
            pnode2->Set(SLAVE, false);
            NodeType::Pointer pnode3 = r_model_part.CreateNewNode(3, 0.0, 0.0, 0.0);
            pnode3->Set(INTERFACE, true);
            pnode3->Set(ACTIVE, true);
            pnode3->Set(MASTER, false);
            pnode3->Set(SLAVE, true);

            pnode1->AddDof(DISPLACEMENT_X);
            pnode1->AddDof(DISPLACEMENT_Y);
            pnode2->AddDof(DISPLACEMENT_X);
            pnode2->AddDof(DISPLACEMENT_Y);
            pnode3->AddDof(DISPLACEMENT_X);
            pnode3->AddDof(DISPLACEMENT_Y);
            pnode3->AddDof(VECTOR_LAGRANGE_MULTIPLIER_X);
            pnode3->AddDof(VECTOR_LAGRANGE_MULTIPLIER_Y);

            std::vector< Dof<double>::Pointer > DoF;
            DoF.reserve(4);
            DoF.push_back(pnode1->pGetDof(DISPLACEMENT_X));
            DoF.push_back(pnode1->pGetDof(DISPLACEMENT_Y));
            DoF.push_back(pnode2->pGetDof(DISPLACEMENT_X));
            DoF.push_back(pnode2->pGetDof(DISPLACEMENT_Y));
            DoF.push_back(pnode3->pGetDof(DISPLACEMENT_X));
            DoF.push_back(pnode3->pGetDof(DISPLACEMENT_Y));
            DoF.push_back(pnode3->pGetDof(VECTOR_LAGRANGE_MULTIPLIER_X));
            DoF.push_back(pnode3->pGetDof(VECTOR_LAGRANGE_MULTIPLIER_Y));

            // Set initial solution
            (pnode1->FastGetSolutionStepValue(DISPLACEMENT)).clear();
            (pnode2->FastGetSolutionStepValue(DISPLACEMENT)).clear();
            (pnode3->FastGetSolutionStepValue(DISPLACEMENT)).clear();
            (pnode3->FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER)).clear();

            DofsArrayType Doftemp;
            Doftemp.reserve(DoF.size());
            for (auto it= DoF.begin(); it!= DoF.end(); it++)
                Doftemp.push_back( *it );

            const std::size_t system_size = 8;
            CompressedMatrix A(system_size, system_size);
            Vector ref_Dx = ZeroVector(system_size);
            Vector Dx = ZeroVector(system_size);
            Vector b = ZeroVector(system_size);
            double count = 0.0;
            for (std::size_t i = 0; i < system_size; ++i) {
                for (std::size_t j = 0; j < system_size; ++j) {
                    if ((((i == 0 || i == 1)&& ((j == system_size - 1) || (j == system_size - 2))) || ((j == 0 || j==1) && ((i == system_size - 1) || (i == system_size - 2))) || (i == 4 && j == 7) || (i == 5 && j == 6)) == false) {
                        count += 1.0;
                        A.push_back(i, j, std::sqrt(count));
                    }
                }
            }

            count = 0.0;
            for (std::size_t i = 0; i < system_size; ++i) {
                count += 1.0;
                b[i] = count;
            }

            // We solve the reference system
            psolver->Solve(A, ref_Dx, b);

            // We solve the block system
            pmixed_solver->ProvideAdditionalData(A, Dx, b, Doftemp, r_model_part);
            pmixed_solver->Solve(A, Dx, b);

            for (std::size_t i = 0; i < system_size; ++i) {
                KRATOS_CHECK_NEAR(ref_Dx[i], Dx[i], tolerance);
            }
        }

        /**
         * Checks if the MixedULMLinear solver performs correctly the resolution of the system. Multiple dof Unordered case
         */
        KRATOS_TEST_CASE_IN_SUITE(MixedULMLinearSolverTwoDoFUnorderedSystem, KratosContactStructuralMechanicsFastSuite)
        {
            constexpr double tolerance = 1e-5;

            Model this_model;
            ModelPart& r_model_part = this_model.CreateModelPart("Main", 3);

            LinearSolverType::Pointer psolver = LinearSolverType::Pointer( new SkylineLUFactorizationSolverType() );
//             Parameters empty_parameters =  Parameters(R"({})");
//             LinearSolverType::Pointer psolver = LinearSolverType::Pointer( new AMGCLSolverType(empty_parameters) );
            LinearSolverType::Pointer pmixed_solver = LinearSolverType::Pointer( new MixedULMLinearSolverType(psolver) );

            r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
            r_model_part.AddNodalSolutionStepVariable(VECTOR_LAGRANGE_MULTIPLIER);

            NodeType::Pointer pnode1 = r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
            NodeType::Pointer pnode3 = r_model_part.CreateNewNode(3, 0.0, 0.0, 0.0);
            pnode3->Set(INTERFACE, true);
            pnode3->Set(ACTIVE, true);
            pnode3->Set(MASTER, false);
            pnode3->Set(SLAVE, true);
            NodeType::Pointer pnode2 = r_model_part.CreateNewNode(2, 0.0, 0.0, 0.0);
            pnode2->Set(INTERFACE, true);
            pnode2->Set(MASTER, true);
            pnode2->Set(SLAVE, false);

            pnode1->AddDof(DISPLACEMENT_X);
            pnode1->AddDof(DISPLACEMENT_Y);
            pnode3->AddDof(DISPLACEMENT_X);
            pnode3->AddDof(DISPLACEMENT_Y);
            pnode3->AddDof(VECTOR_LAGRANGE_MULTIPLIER_X);
            pnode3->AddDof(VECTOR_LAGRANGE_MULTIPLIER_Y);
            pnode2->AddDof(DISPLACEMENT_X);
            pnode2->AddDof(DISPLACEMENT_Y);

            std::vector< Dof<double>::Pointer > DoF;
            DoF.reserve(8);
            DoF.push_back(pnode1->pGetDof(DISPLACEMENT_X));
            DoF.push_back(pnode1->pGetDof(DISPLACEMENT_Y));
            DoF.push_back(pnode3->pGetDof(DISPLACEMENT_X));
            DoF.push_back(pnode3->pGetDof(DISPLACEMENT_Y));
            DoF.push_back(pnode3->pGetDof(VECTOR_LAGRANGE_MULTIPLIER_X));
            DoF.push_back(pnode3->pGetDof(VECTOR_LAGRANGE_MULTIPLIER_Y));
            DoF.push_back(pnode2->pGetDof(DISPLACEMENT_X));
            DoF.push_back(pnode2->pGetDof(DISPLACEMENT_Y));

            // Set initial solution
            (pnode1->FastGetSolutionStepValue(DISPLACEMENT)).clear();
            (pnode3->FastGetSolutionStepValue(DISPLACEMENT)).clear();
            (pnode3->FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER)).clear();
            (pnode2->FastGetSolutionStepValue(DISPLACEMENT)).clear();

            DofsArrayType Doftemp;
            Doftemp.reserve(DoF.size());
            for (auto it= DoF.begin(); it!= DoF.end(); it++)
                Doftemp.push_back( *it );

            const std::size_t system_size = 8;
            CompressedMatrix A(system_size, system_size);
            CompressedMatrix Aaux(system_size, system_size);
            Vector ref_Dx = ZeroVector(system_size);
            Vector Dx = ZeroVector(system_size);
            Vector b = ZeroVector(system_size);
            Vector baux = ZeroVector(system_size);
            double count = 0.0;
            for (std::size_t i = 0; i < system_size; ++i) {
                for (std::size_t j = 0; j < system_size; ++j) {
                    if ((((i == 0 || i == 1)&& ((j == system_size - 1) || (j == system_size - 2))) || ((j == 0 || j==1) && ((i == system_size - 1) || (i == system_size - 2))) || (i == 4 && j == 7) || (i == 5 && j == 6)) == false) {
                        count += 1.0;
                        Aaux.push_back(i, j, std::sqrt(count));
                    }
                }
            }

            for (std::size_t i = 0; i < system_size; ++i) {
                std::size_t iaux = i;
                if (i == 2)
                    iaux = 4;
                else if (i == 3)
                    iaux = 5;
                else if (i == 4)
                    iaux = 6;
                else if (i == 5)
                    iaux = 7;
                else if (i == 6)
                    iaux = 2;
                else if (i == 7)
                    iaux = 3;
                for (std::size_t j = 0; j < system_size; ++j) {
                    std::size_t jaux = j;
                    if (j == 2)
                        jaux = 4;
                    else if (j == 3)
                        jaux = 5;
                    else if (j == 4)
                        jaux = 6;
                    else if (j == 5)
                        jaux = 7;
                    else if (j == 6)
                        jaux = 2;
                    else if (j == 7)
                        jaux = 3;
                    A.push_back(i, j, Aaux(iaux, jaux));
                }
            }

            count = 0.0;
            for (std::size_t i = 0; i < system_size; ++i) {
                count += 1.0;
                baux[i] = count;
            }

            b[0] = baux[0];
            b[1] = baux[1];
            b[2] = baux[4];
            b[3] = baux[5];
            b[4] = baux[6];
            b[5] = baux[7];
            b[6] = baux[2];
            b[7] = baux[3];

            // We solve the reference system
            psolver->Solve(A, ref_Dx, b);

            // We solve the block system
            pmixed_solver->ProvideAdditionalData(A, Dx, b, Doftemp, r_model_part);
            pmixed_solver->Solve(A, Dx, b);

            for (std::size_t i = 0; i < system_size; ++i) {
                KRATOS_CHECK_NEAR(ref_Dx[i], Dx[i], tolerance);
            }
        }

        /**
         * Checks if the MixedULMLinear solver performs correctly the resolution of the system. Multiple dofs (ii)
         */

        KRATOS_TEST_CASE_IN_SUITE(MixedULMLinearSolverThreeDoFSystem, KratosContactStructuralMechanicsFastSuite)
        {
            constexpr double tolerance = 5e-2;

            Model this_model;
            ModelPart& r_model_part = this_model.CreateModelPart("Main", 3);

            LinearSolverType::Pointer psolver = LinearSolverType::Pointer( new SkylineLUFactorizationSolverType() );
//             Parameters empty_parameters =  Parameters(R"({})");
//             LinearSolverType::Pointer psolver = LinearSolverType::Pointer( new AMGCLSolverType(empty_parameters) );
            LinearSolverType::Pointer pmixed_solver = LinearSolverType::Pointer( new MixedULMLinearSolverType(psolver) );

            r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
            r_model_part.AddNodalSolutionStepVariable(VECTOR_LAGRANGE_MULTIPLIER);

            NodeType::Pointer pnode1 = r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
            NodeType::Pointer pnode2 = r_model_part.CreateNewNode(2, 0.0, 0.0, 0.0);

            pnode2->Set(INTERFACE, true);
            pnode2->Set(MASTER, true);
            pnode2->Set(SLAVE, false);
            NodeType::Pointer pnode3 = r_model_part.CreateNewNode(3, 0.0, 0.0, 0.0);
            pnode3->Set(INTERFACE, true);
            pnode3->Set(ACTIVE, true);
            pnode3->Set(MASTER, false);
            pnode3->Set(SLAVE, true);

            pnode1->AddDof(DISPLACEMENT_X);
            pnode1->AddDof(DISPLACEMENT_Y);
            pnode1->AddDof(DISPLACEMENT_Z);
            pnode2->AddDof(DISPLACEMENT_X);
            pnode2->AddDof(DISPLACEMENT_Y);
            pnode2->AddDof(DISPLACEMENT_Z);
            pnode3->AddDof(DISPLACEMENT_X);
            pnode3->AddDof(DISPLACEMENT_Y);
            pnode3->AddDof(DISPLACEMENT_Z);
            pnode3->AddDof(VECTOR_LAGRANGE_MULTIPLIER_X);
            pnode3->AddDof(VECTOR_LAGRANGE_MULTIPLIER_Y);
            pnode3->AddDof(VECTOR_LAGRANGE_MULTIPLIER_Z);

            std::vector< Dof<double>::Pointer > DoF;
            DoF.reserve(12);
            DoF.push_back(pnode1->pGetDof(DISPLACEMENT_X));
            DoF.push_back(pnode1->pGetDof(DISPLACEMENT_Y));
            DoF.push_back(pnode1->pGetDof(DISPLACEMENT_Z));
            DoF.push_back(pnode2->pGetDof(DISPLACEMENT_X));
            DoF.push_back(pnode2->pGetDof(DISPLACEMENT_Y));
            DoF.push_back(pnode2->pGetDof(DISPLACEMENT_Z));
            DoF.push_back(pnode3->pGetDof(DISPLACEMENT_X));
            DoF.push_back(pnode3->pGetDof(DISPLACEMENT_Y));
            DoF.push_back(pnode3->pGetDof(DISPLACEMENT_Z));
            DoF.push_back(pnode3->pGetDof(VECTOR_LAGRANGE_MULTIPLIER_X));
            DoF.push_back(pnode3->pGetDof(VECTOR_LAGRANGE_MULTIPLIER_Y));
            DoF.push_back(pnode3->pGetDof(VECTOR_LAGRANGE_MULTIPLIER_Z));

            // Set initial solution
            (pnode1->FastGetSolutionStepValue(DISPLACEMENT)).clear();
            (pnode2->FastGetSolutionStepValue(DISPLACEMENT)).clear();
            (pnode3->FastGetSolutionStepValue(DISPLACEMENT)).clear();
            (pnode3->FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER)).clear();

            DofsArrayType Doftemp;
            Doftemp.reserve(DoF.size());
            for (auto it= DoF.begin(); it!= DoF.end(); it++)
                Doftemp.push_back( *it );

            const std::size_t system_size = 12;
            CompressedMatrix A(system_size, system_size);
            Vector ref_Dx = ZeroVector(system_size);
            Vector Dx = ZeroVector(system_size);
            Vector b = ZeroVector(system_size);
            double count = 0.0;
            for (std::size_t i = 0; i < system_size; ++i) {
                for (std::size_t j = 0; j < system_size; ++j) {
                    if ((((i == 0 || i == 1 || i == 2) && ((j == system_size - 1) || (j == system_size - 2) || (j == system_size - 3))) || ((j == 0 || j==1 || j==2) && ((i == system_size - 1) || (i == system_size - 2) || (i == system_size - 3))) || (i == 6 && (j == 10 || j == 11)) || (i == 7 && (j == 9 || j == 11)) || (i == 8 && (j == 9 || j == 10))) == false) {
                        count += 1.0;
                        A.push_back(i, j, 5.0 + std::sqrt(count));
                    }
                }
            }

            count = 0.0;
            for (std::size_t i = 0; i < system_size; ++i) {
                count += 1.0;
                b[i] = count;
            }

            // We solve the reference system
            psolver->Solve(A, ref_Dx, b);

            // We solve the block system
            pmixed_solver->ProvideAdditionalData(A, Dx, b, Doftemp, r_model_part);
            pmixed_solver->Solve(A, Dx, b);

            for (std::size_t i = 0; i < system_size; ++i) {
                KRATOS_CHECK_NEAR(std::abs(ref_Dx[i] - Dx[i])/std::abs(ref_Dx[i]), 0.0, tolerance);
            }
        }

        /**
         * Checks if the MixedULMLinear solver performs correctly the resolution of the system. Multiple dof Unordered case (II)
         */
        KRATOS_TEST_CASE_IN_SUITE(MixedULMLinearSolverThreeDoFUnorderedSystem, KratosContactStructuralMechanicsFastSuite)
        {
            constexpr double tolerance = 1e-3;

            Model this_model;
            ModelPart& r_model_part = this_model.CreateModelPart("Main", 3);

            LinearSolverType::Pointer psolver = LinearSolverType::Pointer( new SkylineLUFactorizationSolverType() );
//             Parameters empty_parameters =  Parameters(R"({})");
//             LinearSolverType::Pointer psolver = LinearSolverType::Pointer( new AMGCLSolverType(empty_parameters) );
            LinearSolverType::Pointer pmixed_solver = LinearSolverType::Pointer( new MixedULMLinearSolverType(psolver) );

            r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
            r_model_part.AddNodalSolutionStepVariable(VECTOR_LAGRANGE_MULTIPLIER);

            NodeType::Pointer pnode1 = r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
            NodeType::Pointer pnode3 = r_model_part.CreateNewNode(3, 0.0, 0.0, 0.0);
            pnode3->Set(INTERFACE, true);
            pnode3->Set(ACTIVE, true);
            pnode3->Set(MASTER, false);
            pnode3->Set(SLAVE, true);
            NodeType::Pointer pnode2 = r_model_part.CreateNewNode(2, 0.0, 0.0, 0.0);
            pnode2->Set(INTERFACE, true);
            pnode2->Set(MASTER, true);
            pnode2->Set(SLAVE, false);

            pnode1->AddDof(DISPLACEMENT_X);
            pnode1->AddDof(DISPLACEMENT_Y);
            pnode1->AddDof(DISPLACEMENT_Z);
            pnode3->AddDof(DISPLACEMENT_X);
            pnode3->AddDof(DISPLACEMENT_Y);
            pnode3->AddDof(DISPLACEMENT_Z);
            pnode3->AddDof(VECTOR_LAGRANGE_MULTIPLIER_X);
            pnode3->AddDof(VECTOR_LAGRANGE_MULTIPLIER_Y);
            pnode3->AddDof(VECTOR_LAGRANGE_MULTIPLIER_Z);
            pnode2->AddDof(DISPLACEMENT_X);
            pnode2->AddDof(DISPLACEMENT_Y);
            pnode2->AddDof(DISPLACEMENT_Z);

            std::vector< Dof<double>::Pointer > DoF;
            DoF.reserve(12);
            DoF.push_back(pnode1->pGetDof(DISPLACEMENT_X));
            DoF.push_back(pnode1->pGetDof(DISPLACEMENT_Y));
            DoF.push_back(pnode1->pGetDof(DISPLACEMENT_Z));
            DoF.push_back(pnode3->pGetDof(DISPLACEMENT_X));
            DoF.push_back(pnode3->pGetDof(DISPLACEMENT_Y));
            DoF.push_back(pnode3->pGetDof(DISPLACEMENT_Z));
            DoF.push_back(pnode3->pGetDof(VECTOR_LAGRANGE_MULTIPLIER_X));
            DoF.push_back(pnode3->pGetDof(VECTOR_LAGRANGE_MULTIPLIER_Y));
            DoF.push_back(pnode3->pGetDof(VECTOR_LAGRANGE_MULTIPLIER_Z));
            DoF.push_back(pnode2->pGetDof(DISPLACEMENT_X));
            DoF.push_back(pnode2->pGetDof(DISPLACEMENT_Y));
            DoF.push_back(pnode2->pGetDof(DISPLACEMENT_Z));

            // Set initial solution
            (pnode1->FastGetSolutionStepValue(DISPLACEMENT)).clear();
            (pnode3->FastGetSolutionStepValue(DISPLACEMENT)).clear();
            (pnode3->FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER)).clear();
            (pnode2->FastGetSolutionStepValue(DISPLACEMENT)).clear();

            DofsArrayType Doftemp;
            Doftemp.reserve(DoF.size());
            for (auto it= DoF.begin(); it!= DoF.end(); it++)
                Doftemp.push_back( *it );

            const std::size_t system_size = 12;
            CompressedMatrix A(system_size, system_size);
            CompressedMatrix Aaux(system_size, system_size);
            Vector ref_Dx = ZeroVector(system_size);
            Vector Dx = ZeroVector(system_size);
            Vector b = ZeroVector(system_size);
            Vector baux = ZeroVector(system_size);
            double count = 0.0;
            for (std::size_t i = 0; i < system_size; ++i) {
                for (std::size_t j = 0; j < system_size; ++j) {
                    if ((((i == 0 || i == 1 || i == 2) && ((j == system_size - 1) || (j == system_size - 2) || (j == system_size - 3))) || ((j == 0 || j==1 || j==2) && ((i == system_size - 1) || (i == system_size - 2) || (i == system_size - 3))) || (i == 6 && (j == 10 || j == 11)) || (i == 7 && (j == 9 || j == 11)) || (i == 8 && (j == 9 || j == 10))) == false) {
                        count += 1.0;
                        Aaux.push_back(i, j, std::sqrt(count));
                    }
                }
            }

            for (std::size_t i = 0; i < system_size; ++i) {
                std::size_t iaux = i;
                if (i == 3)
                    iaux = 6;
                else if (i == 4)
                    iaux = 7;
                else if (i == 5)
                    iaux = 8;
                else if (i == 6)
                    iaux = 9;
                else if (i == 7)
                    iaux = 10;
                else if (i == 8)
                    iaux = 11;
                else if (i == 9)
                    iaux = 3;
                else if (i == 10)
                    iaux = 4;
                else if (i == 11)
                    iaux = 5;
                for (std::size_t j = 0; j < system_size; ++j) {
                    std::size_t jaux = j;
                    if (j == 3)
                        jaux = 6;
                    else if (j == 4)
                        jaux = 7;
                    else if (j == 5)
                        jaux = 8;
                    else if (j == 6)
                        jaux = 9;
                    else if (j == 7)
                        jaux = 10;
                    else if (j == 8)
                        jaux = 11;
                    else if (j == 9)
                        jaux = 3;
                    else if (j == 10)
                        jaux = 4;
                    else if (j == 11)
                        jaux = 5;
                    A.push_back(i, j, Aaux(iaux, jaux));
                }
            }

            count = 0.0;
            for (std::size_t i = 0; i < system_size; ++i) {
                count += 1.0;
                baux[i] = count;
            }

            b[0] = baux[0];
            b[1] = baux[1];
            b[2] = baux[2];
            b[3] = baux[6];
            b[4] = baux[7];
            b[5] = baux[8];
            b[6] = baux[9];
            b[7] = baux[10];
            b[8] = baux[11];
            b[9] = baux[3];
            b[10] = baux[4];
            b[11] = baux[5];

            // We solve the reference system
            psolver->Solve(A, ref_Dx, b);

            // We solve the block system
            pmixed_solver->ProvideAdditionalData(A, Dx, b, Doftemp, r_model_part);
            pmixed_solver->Solve(A, Dx, b);

            for (std::size_t i = 0; i < system_size; ++i) {
                KRATOS_CHECK_NEAR(std::abs(ref_Dx[i] - Dx[i])/std::abs(ref_Dx[i]), 0.0, tolerance);
            }
        }

        /**
         * Checks if the MixedULMLinear solver performs correctly the resolution of the system. Multiple dof Unordered case (II)
         */
        KRATOS_TEST_CASE_IN_SUITE(MixedULMLinearSolverRealSystem, KratosContactStructuralMechanicsFastSuite)
        {
            constexpr double tolerance = 1e-3;

            Model this_model;
            ModelPart& r_model_part = this_model.CreateModelPart("Main", 3);

            LinearSolverType::Pointer psolver = LinearSolverType::Pointer( new SkylineLUFactorizationSolverType() );
//             Parameters empty_parameters =  Parameters(R"({})");
//             LinearSolverType::Pointer psolver = LinearSolverType::Pointer( new AMGCLSolverType(empty_parameters) );
            LinearSolverType::Pointer pmixed_solver = LinearSolverType::Pointer( new MixedULMLinearSolverType(psolver) );

            r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
            r_model_part.AddNodalSolutionStepVariable(VECTOR_LAGRANGE_MULTIPLIER);

            NodeType::Pointer pnode1 = r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
            NodeType::Pointer pnode2 = r_model_part.CreateNewNode(2, 0.0, 0.0, 0.0);
            pnode2->Set(INTERFACE, true);
            pnode2->Set(ACTIVE, true);
            pnode2->Set(MASTER, false);
            pnode2->Set(SLAVE, true);
            NodeType::Pointer pnode3 = r_model_part.CreateNewNode(3, 0.0, 0.0, 0.0);
            pnode3->Set(INTERFACE, true);
            pnode3->Set(ACTIVE, true);
            pnode3->Set(MASTER, false);
            pnode3->Set(SLAVE, true);
            NodeType::Pointer pnode4 = r_model_part.CreateNewNode(4, 0.0, 0.0, 0.0);
            pnode4->Set(INTERFACE, true);
            pnode4->Set(MASTER, true);
            pnode4->Set(SLAVE, false);
            NodeType::Pointer pnode5 = r_model_part.CreateNewNode(5, 0.0, 0.0, 0.0);
            pnode5->Set(INTERFACE, true);
            pnode5->Set(MASTER, true);
            pnode5->Set(SLAVE, false);
            NodeType::Pointer pnode6 = r_model_part.CreateNewNode(6, 0.0, 0.0, 0.0);

            pnode1->AddDof(DISPLACEMENT_X);
            pnode1->AddDof(DISPLACEMENT_Y);
            pnode2->AddDof(DISPLACEMENT_X);
            pnode2->AddDof(DISPLACEMENT_Y);
            pnode2->AddDof(VECTOR_LAGRANGE_MULTIPLIER_X);
            pnode2->AddDof(VECTOR_LAGRANGE_MULTIPLIER_Y);
            pnode3->AddDof(DISPLACEMENT_X);
            pnode3->AddDof(DISPLACEMENT_Y);
            pnode3->AddDof(VECTOR_LAGRANGE_MULTIPLIER_X);
            pnode3->AddDof(VECTOR_LAGRANGE_MULTIPLIER_Y);
            pnode4->AddDof(DISPLACEMENT_X);
            pnode4->AddDof(DISPLACEMENT_Y);
            pnode5->AddDof(DISPLACEMENT_X);
            pnode5->AddDof(DISPLACEMENT_Y);
            pnode6->AddDof(DISPLACEMENT_X);
            pnode6->AddDof(DISPLACEMENT_Y);

            std::vector< Dof<double>::Pointer > DoF;
            DoF.reserve(16);
            DoF.push_back(pnode1->pGetDof(DISPLACEMENT_X));
            DoF.push_back(pnode1->pGetDof(DISPLACEMENT_Y));
            DoF.push_back(pnode2->pGetDof(DISPLACEMENT_X));
            DoF.push_back(pnode2->pGetDof(DISPLACEMENT_Y));
            DoF.push_back(pnode2->pGetDof(VECTOR_LAGRANGE_MULTIPLIER_X));
            DoF.push_back(pnode2->pGetDof(VECTOR_LAGRANGE_MULTIPLIER_Y));
            DoF.push_back(pnode3->pGetDof(DISPLACEMENT_X));
            DoF.push_back(pnode3->pGetDof(DISPLACEMENT_Y));
            DoF.push_back(pnode3->pGetDof(VECTOR_LAGRANGE_MULTIPLIER_X));
            DoF.push_back(pnode3->pGetDof(VECTOR_LAGRANGE_MULTIPLIER_Y));
            DoF.push_back(pnode4->pGetDof(DISPLACEMENT_X));
            DoF.push_back(pnode4->pGetDof(DISPLACEMENT_Y));
            DoF.push_back(pnode5->pGetDof(DISPLACEMENT_X));
            DoF.push_back(pnode5->pGetDof(DISPLACEMENT_Y));
            DoF.push_back(pnode6->pGetDof(DISPLACEMENT_X));
            DoF.push_back(pnode6->pGetDof(DISPLACEMENT_Y));

            // Set initial solution
            (pnode1->FastGetSolutionStepValue(DISPLACEMENT)).clear();
            (pnode2->FastGetSolutionStepValue(DISPLACEMENT)).clear();
            (pnode2->FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER)).clear();
            (pnode3->FastGetSolutionStepValue(DISPLACEMENT)).clear();
            (pnode3->FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER)).clear();
            (pnode4->FastGetSolutionStepValue(DISPLACEMENT)).clear();
            (pnode5->FastGetSolutionStepValue(DISPLACEMENT)).clear();
            (pnode6->FastGetSolutionStepValue(DISPLACEMENT)).clear();

            DofsArrayType Doftemp;
            Doftemp.reserve(DoF.size());
            for (auto it= DoF.begin(); it!= DoF.end(); it++)
                Doftemp.push_back( *it );

            const std::size_t system_size = 16;
            CompressedMatrix A(system_size, system_size);
            Vector ref_Dx = ZeroVector(system_size);
            Vector Dx = ZeroVector(system_size);
            Vector b = ZeroVector(system_size);

            // CHANGE THIS TO ADAPT TO YOUR PROBLEM
            CreateAuxiliarFiles();
            const bool read_a = Kratos::ReadMatrixMarketMatrix(FilesystemExtensions::JoinPaths({FilesystemExtensions::CurrentWorkingDirectory(), "A_testing_condensation.mm"}).c_str(), A);
            const bool read_b = Kratos::ReadMatrixMarketVector(FilesystemExtensions::JoinPaths({FilesystemExtensions::CurrentWorkingDirectory(), "b_testing_condensation.rhs"}).c_str(), b);

            // Removing files
            Kratos::filesystem::remove(FilesystemExtensions::JoinPaths({FilesystemExtensions::CurrentWorkingDirectory(), "A_testing_condensation.mm"}));
            Kratos::filesystem::remove(FilesystemExtensions::JoinPaths({FilesystemExtensions::CurrentWorkingDirectory(), "b_testing_condensation.rhs"}));

            if (read_a && read_b) {
                // We solve the reference system
                psolver->Solve(A, ref_Dx, b);

                // We solve the block system
                pmixed_solver->ProvideAdditionalData(A, Dx, b, Doftemp, r_model_part);
                pmixed_solver->Solve(A, Dx, b);

                for (std::size_t i = 0; i < system_size; ++i) {
                    KRATOS_CHECK_NEAR(std::abs(ref_Dx[i] - Dx[i])/std::abs(ref_Dx[i]), 0.0, tolerance);
                }
            }
        }
    } // namespace Testing
}  // namespace Kratos.
