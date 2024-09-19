// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Aron Noordam
//

#include <string>

// Project includes
#include "containers/model.h"
#include "custom_conditions/Pw_normal_flux_condition.hpp"
#include "custom_constitutive/linear_elastic_2D_interface_law.h"
#include "geo_mechanics_fast_suite.h"

namespace Kratos::Testing
{
        KRATOS_TEST_CASE_IN_SUITE(CalculateHorizontalNormalFlux, KratosGeoMechanicsFastSuiteWithoutKernel)
        {

            // initialize modelpart
            Model current_model;
            auto &r_model_part = current_model.CreateModelPart("ModelPart", 1);
            r_model_part.GetProcessInfo().SetValue(DOMAIN_SIZE, 2);
            r_model_part.AddNodalSolutionStepVariable(WATER_PRESSURE);
            r_model_part.AddNodalSolutionStepVariable(NORMAL_FLUID_FLUX);

            // Set the element properties
            auto cond_prop = r_model_part.CreateNewProperties(0);

            // set constitutive law
            cond_prop->SetValue(CONSTITUTIVE_LAW, LinearElastic2DInterfaceLaw().Clone());

            // Create the test piping element nodes
            auto node_1 = r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
            auto node_2 = r_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);

            // set normal fluid flux
            node_1->SetLock();
            node_2->SetLock();

            node_1->FastGetSolutionStepValue(NORMAL_FLUID_FLUX) = 100;
            node_2->FastGetSolutionStepValue(NORMAL_FLUID_FLUX) = 100;

            node_1->UnSetLock();
            node_2->UnSetLock();

            unsigned int conditionSize = 2;

            // Create the test piping element
            std::vector<ModelPart::IndexType> cond_nodes{1, 2};
            Condition::Pointer                cond = Kratos::make_intrusive<PwNormalFluxCondition<2, 2>>(
                1, make_shared<Line2D2<Node>>(node_1, node_2), cond_prop);
            r_model_part.AddCondition(cond);

            // Initialize the condition to initialize the constitutive law
            const auto &r_process_info = r_model_part.GetProcessInfo();
            cond->Initialize(r_process_info);

            // get element geometry
            auto Geom = cond->GetGeometry();
            Matrix rLeftHandSideMatrix = ZeroMatrix(conditionSize, conditionSize);
            Vector rRightHandSideVector = ZeroVector(conditionSize);

            noalias(rLeftHandSideMatrix) = ZeroMatrix(conditionSize, conditionSize);

            cond->CalculateLocalSystem(rLeftHandSideMatrix, rRightHandSideVector, r_process_info);

            // set expected_results
            Vector expected_vector = ZeroVector(conditionSize);
            expected_vector[0] = -50;
            expected_vector[1] = -50;

            for (unsigned int i = 0; i < conditionSize; ++i)
            {
                KRATOS_EXPECT_NEAR(
                    rRightHandSideVector[i],
                    expected_vector[i],
                    1.0e-6);
            }
        }

        KRATOS_TEST_CASE_IN_SUITE(CalculateInclinedNormalFlux, KratosGeoMechanicsFastSuiteWithoutKernel)
        {

            // initialize modelpart
            Model current_model;
            auto &r_model_part = current_model.CreateModelPart("ModelPart", 1);
            r_model_part.GetProcessInfo().SetValue(DOMAIN_SIZE, 2);
            r_model_part.AddNodalSolutionStepVariable(WATER_PRESSURE);
            r_model_part.AddNodalSolutionStepVariable(NORMAL_FLUID_FLUX);

            // Set the element properties
            auto cond_prop = r_model_part.CreateNewProperties(0);

            // set constitutive law
            cond_prop->SetValue(CONSTITUTIVE_LAW, LinearElastic2DInterfaceLaw().Clone());

            // Create the test piping element nodes
            auto node_1 = r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
            auto node_2 = r_model_part.CreateNewNode(2, 3.0, 4.0, 0.0);

            // set normal fluid flux
            node_1->SetLock();
            node_2->SetLock();

            node_1->FastGetSolutionStepValue(NORMAL_FLUID_FLUX) = 100;
            node_2->FastGetSolutionStepValue(NORMAL_FLUID_FLUX) = 100;

            node_1->UnSetLock();
            node_2->UnSetLock();

            unsigned int conditionSize = 2;

            // Create the test piping element
            std::vector<ModelPart::IndexType> cond_nodes{1, 2};
            Condition::Pointer                cond = Kratos::make_intrusive<PwNormalFluxCondition<2, 2>>(
                1, make_shared<Line2D2<Node>>(node_1, node_2), cond_prop);
            r_model_part.AddCondition(cond);            // Initialize the element to initialize the constitutive law
            const auto &r_process_info = r_model_part.GetProcessInfo();
            cond->Initialize(r_process_info);

            // get element geometry
            auto Geom = cond->GetGeometry();
            Matrix rLeftHandSideMatrix = ZeroMatrix(conditionSize, conditionSize);
            Vector rRightHandSideVector = ZeroVector(conditionSize);

            noalias(rLeftHandSideMatrix) = ZeroMatrix(conditionSize, conditionSize);

            cond->CalculateLocalSystem(rLeftHandSideMatrix, rRightHandSideVector, r_process_info);

            // set expected_results
            Vector expected_vector = ZeroVector(conditionSize);
            expected_vector[0] = -250;
            expected_vector[1] = -250;

            for (unsigned int i = 0; i < conditionSize; ++i)
            {
                KRATOS_EXPECT_NEAR(
                    rRightHandSideVector[i],
                    expected_vector[i],
                    1.0e-6);
            }
        }
}
