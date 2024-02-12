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
#include "testing/testing.h"
#include "custom_elements/steady_state_Pw_piping_element.hpp"

namespace Kratos::Testing
{
        KRATOS_TEST_CASE_IN_SUITE(CalculateEquilibriumPipeHeight, KratosGeoMechanicsFastSuite)
        {

            // initialize modelpart
            Model current_model;
            auto &r_model_part = current_model.CreateModelPart("ModelPart", 1);
            r_model_part.GetProcessInfo().SetValue(DOMAIN_SIZE, 2);
            r_model_part.AddNodalSolutionStepVariable(WATER_PRESSURE);
            r_model_part.AddNodalSolutionStepVariable(VOLUME_ACCELERATION);

            // Set the element properties
            auto p_elem_prop = r_model_part.CreateNewProperties(0);
            p_elem_prop->SetValue(PIPE_D_70, 2e-4);
            p_elem_prop->SetValue(PIPE_ETA, 0.25);
            p_elem_prop->SetValue(PIPE_THETA, 30);
            p_elem_prop->SetValue(DENSITY_SOLID, 2650);
            p_elem_prop->SetValue(DENSITY_WATER, 1000);

            p_elem_prop->SetValue(PIPE_ELEMENT_LENGTH, 0.5);
            p_elem_prop->SetValue(PIPE_MODIFIED_D, false);
            p_elem_prop->SetValue(PIPE_MODEL_FACTOR, 1);

            // set constitutive law
            const auto &r_clone_cl = KratosComponents<ConstitutiveLaw>::Get("LinearElastic2DInterfaceLaw");
            p_elem_prop->SetValue(CONSTITUTIVE_LAW, r_clone_cl.Clone());

            // Create the test piping element nodes
            auto p_node_1 = r_model_part.CreateNewNode(1, 0.0, -0.1, 0.0);
            auto p_node_2 = r_model_part.CreateNewNode(2, 1.0, -0.1, 0.0);
            auto p_node_3 = r_model_part.CreateNewNode(3, 1.0, 0.1, 0.0);
            auto p_node_4 = r_model_part.CreateNewNode(4, 0.0, 0.1, 0.0);

            array_1d<double, 3> gravity_array;
            gravity_array[0] = 0;
            gravity_array[1] = 9.81;
            gravity_array[2] = 0;

            // set water pressure values to nodes
            p_node_1->SetLock();
            p_node_2->SetLock();
            p_node_3->SetLock();
            p_node_4->SetLock();

            p_node_1->FastGetSolutionStepValue(WATER_PRESSURE) = 0;
            p_node_2->FastGetSolutionStepValue(WATER_PRESSURE) = 500;
            p_node_3->FastGetSolutionStepValue(WATER_PRESSURE) = 500;
            p_node_4->FastGetSolutionStepValue(WATER_PRESSURE) = 0;

            p_node_1->FastGetSolutionStepValue(VOLUME_ACCELERATION) = gravity_array;
            p_node_2->FastGetSolutionStepValue(VOLUME_ACCELERATION) = gravity_array;
            p_node_3->FastGetSolutionStepValue(VOLUME_ACCELERATION) = gravity_array;
            p_node_4->FastGetSolutionStepValue(VOLUME_ACCELERATION) = gravity_array;

            p_node_1->UnSetLock();
            p_node_2->UnSetLock();
            p_node_3->UnSetLock();
            p_node_4->UnSetLock();

            // Create the test piping element
            std::vector<ModelPart::IndexType> element_nodes{1, 2, 3, 4};
            auto p_element = r_model_part.CreateNewElement("SteadyStatePwPipingElement2D4N", 1, element_nodes, p_elem_prop);

            // Initialize the element to initialize the constitutive law
            const auto &r_process_info = r_model_part.GetProcessInfo();
            p_element->Initialize(r_process_info);

            // get element geometry
            auto Geom = p_element->GetGeometry();

            // cast to piping element
            auto PipeEl = dynamic_cast<SteadyStatePwPipingElement<2, 4> *>(p_element.get());

            // calculate equilibrium height
            double expected_eq_height = PipeEl->CalculateEquilibriumPipeHeight(*p_elem_prop, Geom, p_elem_prop->GetValue(PIPE_ELEMENT_LENGTH));

            KRATOS_EXPECT_NEAR(
                expected_eq_height,
                0.000489,
                1.0e-6);
        }

        KRATOS_TEST_CASE_IN_SUITE(CalculateHeadGradient, KratosGeoMechanicsFastSuite)
        {
            // initialize modelpart
            Model current_model;
            auto &r_model_part = current_model.CreateModelPart("ModelPart", 1);
            r_model_part.GetProcessInfo().SetValue(DOMAIN_SIZE, 2);
            r_model_part.AddNodalSolutionStepVariable(WATER_PRESSURE);
            r_model_part.AddNodalSolutionStepVariable(VOLUME_ACCELERATION);

            // Set the element properties
            auto p_elem_prop = r_model_part.CreateNewProperties(0);
            p_elem_prop->SetValue(PIPE_ELEMENT_LENGTH, 1);
            p_elem_prop->SetValue(DENSITY_WATER, 1000);

            // set constitutive law
            const auto &r_clone_cl = KratosComponents<ConstitutiveLaw>::Get("LinearElastic2DInterfaceLaw");
            p_elem_prop->SetValue(CONSTITUTIVE_LAW, r_clone_cl.Clone());

            // Create the test piping element nodes
            auto p_node_1 = r_model_part.CreateNewNode(1, 0.0, -0.1, 0.0);
            auto p_node_2 = r_model_part.CreateNewNode(2, 1.0, -0.1, 0.0);
            auto p_node_3 = r_model_part.CreateNewNode(3, 1.0, 0.1, 0.0);
            auto p_node_4 = r_model_part.CreateNewNode(4, 0.0, 0.1, 0.0);

            array_1d<double, 3> gravity_array;
            gravity_array[0] = 0;
            gravity_array[1] = 10;
            gravity_array[2] = 0;

            // set water pressure values to nodes
            p_node_1->SetLock();
            p_node_2->SetLock();
            p_node_3->SetLock();
            p_node_4->SetLock();

            p_node_1->FastGetSolutionStepValue(WATER_PRESSURE) = 0;
            p_node_2->FastGetSolutionStepValue(WATER_PRESSURE) = 20000;
            p_node_3->FastGetSolutionStepValue(WATER_PRESSURE) = 20000;
            p_node_4->FastGetSolutionStepValue(WATER_PRESSURE) = 0;

            p_node_1->FastGetSolutionStepValue(VOLUME_ACCELERATION) = gravity_array;
            p_node_2->FastGetSolutionStepValue(VOLUME_ACCELERATION) = gravity_array;
            p_node_3->FastGetSolutionStepValue(VOLUME_ACCELERATION) = gravity_array;
            p_node_4->FastGetSolutionStepValue(VOLUME_ACCELERATION) = gravity_array;

            p_node_1->UnSetLock();
            p_node_2->UnSetLock();
            p_node_3->UnSetLock();
            p_node_4->UnSetLock();

            // Create the test piping element
            std::vector<ModelPart::IndexType> element_nodes{1, 2, 3, 4};
            auto p_element = r_model_part.CreateNewElement("SteadyStatePwPipingElement2D4N", 1, element_nodes, p_elem_prop);

            // Initialize the element to initialize the constitutive law
            const auto &r_process_info = r_model_part.GetProcessInfo();
            p_element->Initialize(r_process_info);

            // get element geometry
            auto Geom = p_element->GetGeometry();

            // cast to piping element
            auto PipeEl = dynamic_cast<SteadyStatePwPipingElement<2, 4> *>(p_element.get());

            // calculate water pressure gradient
            double expected_gradient = PipeEl->CalculateHeadGradient(*p_elem_prop, Geom, p_elem_prop->GetValue(PIPE_ELEMENT_LENGTH));

            // assert gradient
            // expected gradient should be 2. Test is failing on purpose to check CI
            KRATOS_EXPECT_NEAR(
                expected_gradient,
                2,
                1.0e-10);
        }
}