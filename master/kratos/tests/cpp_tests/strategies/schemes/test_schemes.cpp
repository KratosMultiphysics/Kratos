//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
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
#include "spaces/ublas_space.h"
#include "solving_strategies/schemes/residual_based_newmark_displacement_scheme.hpp"
#include "solving_strategies/schemes/residual_based_bossak_displacement_scheme.hpp"
#include "solving_strategies/schemes/residual_based_bdf_displacement_scheme.h"
#include "solving_strategies/schemes/residual_based_bdf_custom_scheme.h"

namespace Kratos
{
    namespace Testing
    {
        /// Tests

        typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
        typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;

        typedef Scheme< SparseSpaceType, LocalSpaceType >  SchemeType;
        typedef ModelPart::DofsArrayType DofsArrayType;
        /**
         * @brief Common method to set the displacement scheme
         */
        static DofsArrayType BasicTestSchemeDisplacement(
            ModelPart& rModelPart,
            SchemeType::Pointer pScheme,
            std::vector< Dof<double>::Pointer >& rDoF,
            const double DeltaTime = 1.0e-4
            )
        {
            const std::size_t buffer_size = 3;

            rModelPart.SetBufferSize(buffer_size);
            rModelPart.GetProcessInfo().SetValue(DOMAIN_SIZE, 3);

            rModelPart.AddNodalSolutionStepVariable(DISPLACEMENT);
            rModelPart.AddNodalSolutionStepVariable(VELOCITY);
            rModelPart.AddNodalSolutionStepVariable(ACCELERATION);

            auto pnode = rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);

            pnode->AddDof(DISPLACEMENT_X);
            pnode->AddDof(DISPLACEMENT_Y);
            pnode->AddDof(DISPLACEMENT_Z);
            pnode->AddDof(VELOCITY_X);
            pnode->AddDof(VELOCITY_Y);
            pnode->AddDof(VELOCITY_Z);
            pnode->AddDof(ACCELERATION_X);
            pnode->AddDof(ACCELERATION_Y);
            pnode->AddDof(ACCELERATION_Z);

            rModelPart.GetProcessInfo()[DELTA_TIME] = DeltaTime;

            rDoF.reserve(9);
            rDoF.push_back(pnode->pGetDof(DISPLACEMENT_X));
            rDoF.push_back(pnode->pGetDof(DISPLACEMENT_Y));
            rDoF.push_back(pnode->pGetDof(DISPLACEMENT_Z));
            rDoF.push_back(pnode->pGetDof(VELOCITY_X));
            rDoF.push_back(pnode->pGetDof(VELOCITY_Y));
            rDoF.push_back(pnode->pGetDof(VELOCITY_Z));
            rDoF.push_back(pnode->pGetDof(ACCELERATION_X));
            rDoF.push_back(pnode->pGetDof(ACCELERATION_Y));
            rDoF.push_back(pnode->pGetDof(ACCELERATION_Z));

            // Set initial solution
            const array_1d<double, 3> zero_vector = ZeroVector(3);
            for (std::size_t i = 0; i < buffer_size; ++i) {
                pnode->FastGetSolutionStepValue(DISPLACEMENT, i) = zero_vector;
                pnode->FastGetSolutionStepValue(VELOCITY, i) = zero_vector;
                pnode->FastGetSolutionStepValue(ACCELERATION, i) = zero_vector;
            }

            DofsArrayType Doftemp;
            Doftemp.reserve(rDoF.size());
            for (auto it= rDoF.begin(); it!= rDoF.end(); ++it) {
                Doftemp.push_back( *it );
            }
            Doftemp.Sort();

            CompressedMatrix A(ZeroMatrix(3, 3));
            Vector Dx = ZeroVector(3);
            Vector b = ZeroVector(3);

            pScheme->Initialize(rModelPart);

            return Doftemp;
        }

        /**
         * @brief Common method to test the schemes
         */
        static void TestScheme(
            SchemeType::Pointer pScheme,
            const double DeltaTime = 1.0e-4,
            const std::string TestType = "DISPLACEMENT",
            const bool TestPredict = false,
            const bool TestUpdatePredict = false
            )
        {
            Model current_model;

            constexpr double tolerance = 1e-6;

            ModelPart& r_model_part = current_model.CreateModelPart("Main");

            std::vector< Dof<double>::Pointer > DoF;
            DofsArrayType Doftemp;
            if (TestType == "DISPLACEMENT")
                Doftemp = BasicTestSchemeDisplacement(r_model_part, pScheme, DoF, DeltaTime);
            else
                KRATOS_ERROR << "Case not implemented" << std::endl;

            CompressedMatrix A = ZeroMatrix(3, 3);
            Vector Dx = ZeroVector(3);
            Vector b = ZeroVector(3);

            // Check InitializeSolutionStep and Update
            double time = 0;

            auto pnode = r_model_part.pGetNode(1);

            pnode->FastGetSolutionStepValue(DISPLACEMENT_X) = std::cos(time);
            pnode->FastGetSolutionStepValue(DISPLACEMENT_X, 1) = std::cos(time - DeltaTime);
            pnode->FastGetSolutionStepValue(DISPLACEMENT_X, 2) = std::cos(time - 2.0 * DeltaTime);
            pnode->FastGetSolutionStepValue(VELOCITY_X) = -std::sin(time);
            pnode->FastGetSolutionStepValue(VELOCITY_X, 1) = -std::sin(time - DeltaTime);
            pnode->FastGetSolutionStepValue(VELOCITY_X, 2) = -std::sin(time - 2.0 * DeltaTime);
            pnode->FastGetSolutionStepValue(ACCELERATION_X) = -std::cos(time);
            pnode->FastGetSolutionStepValue(ACCELERATION_X, 1) = -std::cos(time - DeltaTime);
            pnode->FastGetSolutionStepValue(ACCELERATION_X, 2) = -std::cos(time - 2.0 * DeltaTime);

            pScheme->Initialize(r_model_part);

            const std::size_t number_iterations = 10;
            for (std::size_t iter = 0; iter < number_iterations; ++iter) {
                time += DeltaTime;

                r_model_part.CloneTimeStep(time);

                Dx[0] = std::cos(time) - std::cos(time - DeltaTime);

                pScheme->InitializeSolutionStep(r_model_part, A, Dx, b);
                pScheme->Update(r_model_part, Doftemp, A, Dx, b);

                const double x = pnode->FastGetSolutionStepValue(DISPLACEMENT_X);
                const double v = pnode->FastGetSolutionStepValue(VELOCITY_X);
                const double a = pnode->FastGetSolutionStepValue(ACCELERATION_X);

//                 // Debug
//                 std::cout << time << "\t" << x << "\t" << v << "\t" << a << std::endl;

                KRATOS_EXPECT_LE(std::abs(x - std::cos(time)), tolerance);
                KRATOS_EXPECT_LE(std::abs(v + std::sin(time)), tolerance);
                KRATOS_EXPECT_LE(std::abs(a + std::cos(time)), tolerance);
            }

            // Test updates with fixed velocities and accelerations
            if (TestUpdatePredict) {
                Dx[0] = 0.0;

                // Check Update (fix velocity)
                pScheme->Clear();
                pnode->pGetDof(DISPLACEMENT_X)->FreeDof();
                pnode->pGetDof(VELOCITY_X)->FixDof();
                time = -DeltaTime;
                r_model_part.CloneTimeStep(time);
                time = 0.0;
                r_model_part.CloneTimeStep(time);

                pnode->FastGetSolutionStepValue(DISPLACEMENT_X) = std::cos(time);
                pnode->FastGetSolutionStepValue(DISPLACEMENT_X, 1) = std::cos(time - DeltaTime);
                pnode->FastGetSolutionStepValue(DISPLACEMENT_X, 2) = std::cos(time - 2.0 * DeltaTime);
                pnode->FastGetSolutionStepValue(VELOCITY_X) = -std::sin(time);
                pnode->FastGetSolutionStepValue(VELOCITY_X, 1) = -std::sin(time - DeltaTime);
                pnode->FastGetSolutionStepValue(VELOCITY_X, 2) = -std::sin(time - 2.0 * DeltaTime);
                pnode->FastGetSolutionStepValue(ACCELERATION_X) = -std::cos(time);
                pnode->FastGetSolutionStepValue(ACCELERATION_X, 1) = -std::cos(time - DeltaTime);
                pnode->FastGetSolutionStepValue(ACCELERATION_X, 2) = -std::cos(time - 2.0 * DeltaTime);

                for (std::size_t iter = 0; iter < number_iterations; ++iter) {
                    time += DeltaTime;

                    r_model_part.CloneTimeStep(time);

                    pnode->FastGetSolutionStepValue(VELOCITY_X) = - std::sin(time);

                    pScheme->InitializeSolutionStep(r_model_part, A, Dx, b);
                    pScheme->Predict(r_model_part, Doftemp, A, Dx, b);
                    pScheme->Update(r_model_part, Doftemp, A, Dx, b);

                    const double x = pnode->FastGetSolutionStepValue(DISPLACEMENT_X);
                    const double v = pnode->FastGetSolutionStepValue(VELOCITY_X);
                    const double a = pnode->FastGetSolutionStepValue(ACCELERATION_X);

//                     // Debug
//                     std::cout << time << "\t" << x << "\t" << v << "\t" << a << std::endl;

                    KRATOS_EXPECT_LE(std::abs(x - std::cos(time)), tolerance);
                    KRATOS_EXPECT_LE(std::abs(v + std::sin(time)), tolerance);
                    KRATOS_EXPECT_LE(std::abs(a + std::cos(time)), tolerance);
                }

                // Check Update (fix acceleration)
                pScheme->Clear();
                pnode->pGetDof(DISPLACEMENT_X)->FreeDof();
                pnode->pGetDof(VELOCITY_X)->FreeDof();
                pnode->pGetDof(ACCELERATION_X)->FixDof();
                time = -DeltaTime;
                r_model_part.CloneTimeStep(time);
                time = 0.0;
                r_model_part.CloneTimeStep(time);

                pnode->FastGetSolutionStepValue(DISPLACEMENT_X) = std::cos(time);
                pnode->FastGetSolutionStepValue(DISPLACEMENT_X, 1) = std::cos(time - DeltaTime);
                pnode->FastGetSolutionStepValue(DISPLACEMENT_X, 2) = std::cos(time - 2.0 * DeltaTime);
                pnode->FastGetSolutionStepValue(VELOCITY_X) = -std::sin(time);
                pnode->FastGetSolutionStepValue(VELOCITY_X, 1) = -std::sin(time - DeltaTime);
                pnode->FastGetSolutionStepValue(VELOCITY_X, 2) = -std::sin(time - 2.0 * DeltaTime);
                pnode->FastGetSolutionStepValue(ACCELERATION_X) = -std::cos(time);
                pnode->FastGetSolutionStepValue(ACCELERATION_X, 1) = -std::cos(time - DeltaTime);
                pnode->FastGetSolutionStepValue(ACCELERATION_X, 2) = -std::cos(time - 2.0 * DeltaTime);

                for (std::size_t iter = 0; iter < number_iterations; ++iter) {
                    time += DeltaTime;

                    r_model_part.CloneTimeStep(time);

                    pnode->FastGetSolutionStepValue(ACCELERATION_X) = - std::cos(time);

                    pScheme->InitializeSolutionStep(r_model_part, A, Dx, b);
                    pScheme->Predict(r_model_part, Doftemp, A, Dx, b);
                    pScheme->Update(r_model_part, Doftemp, A, Dx, b);

                    const double x = pnode->FastGetSolutionStepValue(DISPLACEMENT_X);
                    const double v = pnode->FastGetSolutionStepValue(VELOCITY_X);
                    const double a = pnode->FastGetSolutionStepValue(ACCELERATION_X);

//                     // Debug
//                     std::cout << time << "\t" << x << "\t" << v << "\t" << a << std::endl;

                    KRATOS_EXPECT_LE(std::abs(x - std::cos(time)), tolerance);
                    KRATOS_EXPECT_LE(std::abs(v + std::sin(time)), tolerance);
                    KRATOS_EXPECT_LE(std::abs(a + std::cos(time)), tolerance);
                }
            }

            // Check Predict (displacement)
            if (TestPredict) {
                pnode->pGetDof(DISPLACEMENT_X)->FixDof();
                time = 0;

                pnode->FastGetSolutionStepValue(DISPLACEMENT_X) = std::cos(time);
                pnode->FastGetSolutionStepValue(DISPLACEMENT_X, 1) = std::cos(time - DeltaTime);
                pnode->FastGetSolutionStepValue(DISPLACEMENT_X, 2) = std::cos(time - 2.0 * DeltaTime);
                pnode->FastGetSolutionStepValue(VELOCITY_X) = -std::sin(time);
                pnode->FastGetSolutionStepValue(VELOCITY_X, 1) = -std::sin(time - DeltaTime);
                pnode->FastGetSolutionStepValue(VELOCITY_X, 2) = -std::sin(time - 2.0 * DeltaTime);
                pnode->FastGetSolutionStepValue(ACCELERATION_X) = -std::cos(time);
                pnode->FastGetSolutionStepValue(ACCELERATION_X, 1) = -std::cos(time - DeltaTime);
                pnode->FastGetSolutionStepValue(ACCELERATION_X, 2) = -std::cos(time - 2.0 * DeltaTime);

                for (std::size_t iter = 0; iter < number_iterations; ++iter) {
                    time += DeltaTime;

                    r_model_part.CloneTimeStep(time);

                    pnode->FastGetSolutionStepValue(DISPLACEMENT_X) = std::cos(time);

                    pScheme->Predict(r_model_part, Doftemp, A, Dx, b);

                    const double x = pnode->FastGetSolutionStepValue(DISPLACEMENT_X);
                    const double v = pnode->FastGetSolutionStepValue(VELOCITY_X);
                    const double a = pnode->FastGetSolutionStepValue(ACCELERATION_X);

//                     // Debug
//                     std::cout << time << "\t" << x << "\t" << v << "\t" << a << std::endl;

                    KRATOS_EXPECT_LE(std::abs(x - std::cos(time)), tolerance);
                    KRATOS_EXPECT_LE(std::abs(v + std::sin(time)), tolerance);
                    KRATOS_EXPECT_LE(std::abs(a + std::cos(time)), tolerance);
                }

                // Check Predict (velocity)
                pnode->pGetDof(DISPLACEMENT_X)->FreeDof();
                pnode->pGetDof(VELOCITY_X)->FixDof();
                time = 0.0;

                pnode->FastGetSolutionStepValue(DISPLACEMENT_X) = std::cos(time);
                pnode->FastGetSolutionStepValue(DISPLACEMENT_X, 1) = std::cos(time - DeltaTime);
                pnode->FastGetSolutionStepValue(DISPLACEMENT_X, 2) = std::cos(time - 2.0 * DeltaTime);
                pnode->FastGetSolutionStepValue(VELOCITY_X) = -std::sin(time);
                pnode->FastGetSolutionStepValue(VELOCITY_X, 1) = -std::sin(time - DeltaTime);
                pnode->FastGetSolutionStepValue(VELOCITY_X, 2) = -std::sin(time - 2.0 * DeltaTime);
                pnode->FastGetSolutionStepValue(ACCELERATION_X) = -std::cos(time);
                pnode->FastGetSolutionStepValue(ACCELERATION_X, 1) = -std::cos(time - DeltaTime);
                pnode->FastGetSolutionStepValue(ACCELERATION_X, 2) = -std::cos(time - 2.0 * DeltaTime);

                for (std::size_t iter = 0; iter < number_iterations; ++iter) {
                    time += DeltaTime;

                    r_model_part.CloneTimeStep(time);

                    pnode->FastGetSolutionStepValue(VELOCITY_X) = - std::sin(time);

                    pScheme->Predict(r_model_part, Doftemp, A, Dx, b);

                    const double x = pnode->FastGetSolutionStepValue(DISPLACEMENT_X);
                    const double v = pnode->FastGetSolutionStepValue(VELOCITY_X);
                    const double a = pnode->FastGetSolutionStepValue(ACCELERATION_X);

//                     // Debug
//                     std::cout << time << "\t" << x << "\t" << v << "\t" << a << std::endl;

                    KRATOS_EXPECT_LE(std::abs(x - std::cos(time)), tolerance);
                    KRATOS_EXPECT_LE(std::abs(v + std::sin(time)), tolerance);
                    KRATOS_EXPECT_LE(std::abs(a + std::cos(time)), tolerance);
                }

                // Check Predict (acceleration)
                pnode->pGetDof(VELOCITY_X)->FreeDof();
                pnode->pGetDof(ACCELERATION_X)->FixDof();
                time = 0;

                pnode->FastGetSolutionStepValue(DISPLACEMENT_X) = std::cos(time);
                pnode->FastGetSolutionStepValue(DISPLACEMENT_X, 1) = std::cos(time - DeltaTime);
                pnode->FastGetSolutionStepValue(DISPLACEMENT_X, 2) = std::cos(time - 2.0 * DeltaTime);
                pnode->FastGetSolutionStepValue(VELOCITY_X) = -std::sin(time);
                pnode->FastGetSolutionStepValue(VELOCITY_X, 1) = -std::sin(time - DeltaTime);
                pnode->FastGetSolutionStepValue(VELOCITY_X, 2) = -std::sin(time - 2.0 * DeltaTime);
                pnode->FastGetSolutionStepValue(ACCELERATION_X) = -std::cos(time);
                pnode->FastGetSolutionStepValue(ACCELERATION_X, 1) = -std::cos(time - DeltaTime);
                pnode->FastGetSolutionStepValue(ACCELERATION_X, 2) = -std::cos(time - 2.0 * DeltaTime);

                for (std::size_t iter = 0; iter < number_iterations; ++iter) {
                    time += DeltaTime;

                    r_model_part.CloneTimeStep(time);

                    pnode->FastGetSolutionStepValue(ACCELERATION_X) = - std::cos(time);

                    pScheme->Predict(r_model_part, Doftemp, A, Dx, b);

                    const double x = pnode->FastGetSolutionStepValue(DISPLACEMENT_X);
                    const double v = pnode->FastGetSolutionStepValue(VELOCITY_X);
                    const double a = pnode->FastGetSolutionStepValue(ACCELERATION_X);

//                     // Debug
//                     std::cout << time << "\t" << x << "\t" << v << "\t" << a << std::endl;

                    KRATOS_EXPECT_LE(std::abs(x - std::cos(time)), tolerance);
                    KRATOS_EXPECT_LE(std::abs(v + std::sin(time)), tolerance);
                    KRATOS_EXPECT_LE(std::abs(a + std::cos(time)), tolerance);
                }
            }
        }

        /**
         * Checks if the Newmark scheme performs correctly the integration
         */
        KRATOS_TEST_CASE_IN_SUITE(DisplacementNewmarkSchemeTest, KratosCoreFastSuite)
        {
            typedef ResidualBasedNewmarkDisplacementScheme< SparseSpaceType, LocalSpaceType >  ResidualBasedNewmarkDisplacementSchemeType;
            SchemeType::Pointer p_scheme = SchemeType::Pointer( new ResidualBasedNewmarkDisplacementSchemeType() );

            const double delta_time = 1.0e-4;

            TestScheme(p_scheme, delta_time, "DISPLACEMENT", true, true);
        }

        /**
         * Checks if the Bossak scheme performs correctly the integration
         */
        KRATOS_TEST_CASE_IN_SUITE(DisplacementBossakSchemeTest, KratosCoreFastSuite)
        {
            typedef ResidualBasedBossakDisplacementScheme< SparseSpaceType, LocalSpaceType >  ResidualBasedBossakDisplacementSchemeType;
            SchemeType::Pointer p_scheme = SchemeType::Pointer( new ResidualBasedBossakDisplacementSchemeType() );

            const double delta_time = 1.0e-4;

            TestScheme(p_scheme, delta_time, "DISPLACEMENT", true, true);
        }

        /**
         * Checks if the BDF2 scheme performs correctly the integration
         */
        KRATOS_TEST_CASE_IN_SUITE(DisplacementBDF2SchemeTest, KratosCoreFastSuite)
        {
            typedef ResidualBasedBDFDisplacementScheme< SparseSpaceType, LocalSpaceType >  ResidualBasedBDFDisplacementSchemeType;
            SchemeType::Pointer p_scheme = SchemeType::Pointer( new ResidualBasedBDFDisplacementSchemeType(2) );

            const double delta_time = 1.0e-4;

            TestScheme(p_scheme, delta_time, "DISPLACEMENT", true, true);
        }

        /**
         * Checks if the custom BDF2 scheme performs correctly the integration
         */
        KRATOS_TEST_CASE_IN_SUITE(CustomBDF2SchemeTest, KratosCoreFastSuite)
        {
            typedef ResidualBasedBDFCustomScheme< SparseSpaceType, LocalSpaceType >  ResidualBasedBDFCustomSchemeType;
            SchemeType::Pointer p_scheme = SchemeType::Pointer( new ResidualBasedBDFCustomSchemeType(2) );

            const double delta_time = 1.0e-4;

            TestScheme(p_scheme, delta_time, "DISPLACEMENT", true, true);
        }

    } // namespace Testing
}  // namespace Kratos.

