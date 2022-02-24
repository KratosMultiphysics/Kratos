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
//#include "stdafx.h"
#include <string>
#include <iostream>
// Project includes
#include "containers/model.h"
#include "testing/testing.h"
//#include "structural_mechanics_application_variables.h"
//#include "custom_elements/total_lagrangian.h"
//#include "custom_processes/set_cartesian_local_axes_process.h"
//#include "custom_processes/set_cylindrical_local_axes_process.h"
//#include "custom_processes/set_spherical_local_axes_process.h"
#include "custom_elements/steady_state_Pw_piping_element.hpp"
//#include "custom_elements/steady_state_Pw_piping_element.cpp"




using namespace std;

namespace Kratos
{
    namespace Testing
    {

        void AddWaterPressureDofs(ModelPart& rModelPart) {
            for (auto& r_node : rModelPart.Nodes()) {
                r_node.AddDof(WATER_PRESSURE);
            }
        }


        KRATOS_TEST_CASE_IN_SUITE(CalculateEquilibriumPipeHeight, KratosGeoMechanicsFastSuite)
        {
            // todo make actual test

            Model current_model;
            auto& r_model_part = current_model.CreateModelPart("ModelPart", 1);
            r_model_part.GetProcessInfo().SetValue(DOMAIN_SIZE, 2);

            r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
            r_model_part.AddNodalSolutionStepVariable(VOLUME_ACCELERATION);

            // Set the element properties
            auto p_elem_prop = r_model_part.CreateNewProperties(0);
            p_elem_prop->SetValue(PIPE_D_70, 2e-4);
            p_elem_prop->SetValue(PIPE_ETA, 0.25);
            p_elem_prop->SetValue(PIPE_THETA, 30);
            p_elem_prop->SetValue(DENSITY_SOLID, 2650);
            p_elem_prop->SetValue(DENSITY_WATER, 1000);


            p_elem_prop->SetValue(THICKNESS, 0.01);

            
            /*KRATOS_CHECK_NEAR(
                1,
                5,
                1.0e-10);*/

        }

        KRATOS_TEST_CASE_IN_SUITE(CalculateWaterPressureGradient, KratosGeoMechanicsFastSuite)
        {
            // initialize modelpart
            Model current_model;
            auto& r_model_part = current_model.CreateModelPart("ModelPart", 1);
            r_model_part.GetProcessInfo().SetValue(DOMAIN_SIZE, 2);
            r_model_part.AddNodalSolutionStepVariable(WATER_PRESSURE);


            // Set the element properties
            auto p_elem_prop = r_model_part.CreateNewProperties(0);
            p_elem_prop->SetValue(PIPE_ELEMENT_LENGTH, 1);

            // set constitutive law
            const auto& r_clone_cl = KratosComponents<ConstitutiveLaw>::Get("LinearElastic2DInterfaceLaw");
            p_elem_prop->SetValue(CONSTITUTIVE_LAW, r_clone_cl.Clone());

            // Create the test piping element nodes
            auto p_node_1 = r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
            auto p_node_2 = r_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
            auto p_node_3 = r_model_part.CreateNewNode(3, 1.0, 0.1, 0.0);
            auto p_node_4 = r_model_part.CreateNewNode(4, 0.0, 0.1, 0.0);

            // set water pressure values to node 1 and 2
            p_node_1->SetLock();
            p_node_2->SetLock();

            p_node_1->FastGetSolutionStepValue(WATER_PRESSURE) = 0;
            p_node_2->FastGetSolutionStepValue(WATER_PRESSURE) = 2;

            p_node_1->UnSetLock();
            p_node_2->UnSetLock();

            // Create the test piping element
            std::vector<ModelPart::IndexType> element_nodes{ 1,2,3,4};
            auto p_element = r_model_part.CreateNewElement("SteadyStatePwPipingElement2D4N", 1, element_nodes, p_elem_prop);

            // Initialize the element to initialize the constitutive law
            const auto& r_process_info = r_model_part.GetProcessInfo();
            p_element->Initialize(r_process_info); 

            // get element geometry
            auto Geom  = p_element->GetGeometry();

            // cast to piping element
            typedef typename SteadyStatePwPipingElement<2, 4>                    SteadyStatePwPipingElement;
            SteadyStatePwPipingElement PipeEl = *p_element;

            // calculate water pressure gradient
            double expected_gradient = PipeEl.CalculateWaterPressureGradient(*p_elem_prop, Geom);

            // assert gradient
            // expected gradient should be 2. Test is failing on purpose to check CI
            KRATOS_CHECK_NEAR(
                expected_gradient,
                5,
                1.0e-10);

        }
    }
}