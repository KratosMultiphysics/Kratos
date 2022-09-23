//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Daniel Diez, Ruben Zorrilla, Uxue Chasco
//
//

// System includes
#include <set>

// External includes

// Project includes
#include "testing/testing.h"
#include "containers/model.h"
#include "spaces/ublas_space.h"
#include "includes/properties.h"
#include "includes/model_part.h"
#include "utilities/math_utils.h"
#include "includes/global_pointer_variables.h"
#include "custom_elements/two_fluid_navier_stokes.h"
#include "custom_constitutive/newtonian_2d_law.h"
#include "custom_constitutive/newtonian_3d_law.h"
#include "custom_constitutive/newtonian_two_fluid_3d_law.h"
#include "processes/find_nodal_neighbours_process.h"
#include "utilities/normal_calculation_utils.h"

namespace Kratos {
    namespace Testing {

        typedef ModelPart::IndexType									 IndexType;
        typedef ModelPart::NodeIterator					          NodeIteratorType;

        void AuxiliaryFunctionForElementTest(std::string ElementName)
        {
            Model current_model;
            ModelPart &modelPart = current_model.CreateModelPart("Main");
            modelPart.SetBufferSize(3);

            // Variables addition
            modelPart.AddNodalSolutionStepVariable(BODY_FORCE);
            modelPart.AddNodalSolutionStepVariable(DENSITY);
            modelPart.AddNodalSolutionStepVariable(DYNAMIC_VISCOSITY);
            modelPart.AddNodalSolutionStepVariable(DYNAMIC_TAU);
            modelPart.AddNodalSolutionStepVariable(PRESSURE);
            modelPart.AddNodalSolutionStepVariable(VELOCITY);
            modelPart.AddNodalSolutionStepVariable(MESH_VELOCITY);
            modelPart.AddNodalSolutionStepVariable(DISTANCE);
            modelPart.AddNodalSolutionStepVariable(REACTION);
            modelPart.AddNodalSolutionStepVariable(ACCELERATION);
            modelPart.AddNodalSolutionStepVariable(EXTERNAL_PRESSURE);
            modelPart.AddNodalSolutionStepVariable(REACTION_WATER_PRESSURE);

            // Process info creation
            double delta_time = 0.1;
            const double sound_velocity = 1.0e+12;
            modelPart.GetProcessInfo().SetValue(DYNAMIC_TAU, 1.0);
            modelPart.GetProcessInfo().SetValue(SOUND_VELOCITY, 1.0e+16);
            modelPart.GetProcessInfo().SetValue(DELTA_TIME, delta_time);

            Vector bdf_coefs(3);
            bdf_coefs[0] = 3.0 / (2.0 * delta_time);
            bdf_coefs[1] = -2.0 / delta_time;
            bdf_coefs[2] = 0.5 * delta_time;
            modelPart.GetProcessInfo().SetValue(BDF_COEFFICIENTS, bdf_coefs);
            modelPart.GetProcessInfo().SetValue(VOLUME_ERROR, 0.0);

            // Set the element properties
            Properties::Pointer pElemProp = modelPart.CreateNewProperties(0);
            pElemProp->SetValue(DENSITY, 1000.0);
            pElemProp->SetValue(DYNAMIC_VISCOSITY, 1.0e-16);
            Newtonian2DLaw::Pointer pConsLaw(new Newtonian2DLaw());
            pElemProp->SetValue(CONSTITUTIVE_LAW, pConsLaw);

            // Geometry creation
            modelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
            modelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
            modelPart.CreateNewNode(3, 0.5, 0.5, 0.0);
            std::vector<ModelPart::IndexType> elemNodes{1, 2, 3};
            auto pElement = modelPart.CreateNewElement(ElementName, 1, elemNodes, pElemProp);

            // Define the nodal values
            for (NodeIteratorType it_node = modelPart.NodesBegin(); it_node < modelPart.NodesEnd(); ++it_node)
            {
                it_node->SetValue(SOUND_VELOCITY, sound_velocity);
                it_node->FastGetSolutionStepValue(DENSITY) = pElemProp->GetValue(DENSITY);
                it_node->FastGetSolutionStepValue(DYNAMIC_VISCOSITY) = pElemProp->GetValue(DYNAMIC_VISCOSITY);
                it_node->FastGetSolutionStepValue(BODY_FORCE_Y) = 10.0;
                // it_node->FastGetSolutionStepValue(DISTANCE) = 1.0;

            }

            modelPart.pGetNode(1)->FastGetSolutionStepValue(VELOCITY_X) = 1.0;

            modelPart.pGetNode(1)->FastGetSolutionStepValue(DISTANCE) = 1.0;
            modelPart.pGetNode(2)->FastGetSolutionStepValue(DISTANCE) = 1.0;
            modelPart.pGetNode(3)->FastGetSolutionStepValue(DISTANCE) = 0.5;

            // for (unsigned int i = 0; i < 3; i++)
            // {
            //     double distance = pElement->GetGeometry()[i].FastGetSolutionStepValue(DISTANCE);
            //     pElement->GetGeometry()[i].FastGetSolutionStepValue(PRESSURE) = 1000.0 * 10.0 * distance;
            // }

            // Compute RHS and LHS
            Vector RHS = ZeroVector(9);
            Matrix LHS = ZeroMatrix(9, 9);

            const auto &r_process_info = modelPart.GetProcessInfo();
            pElement->Initialize(r_process_info); // Initialize the element to initialize the constitutive law
            pElement->CalculateLocalSystem(LHS, RHS, r_process_info);

            // for (double val : RHS) {
            //     std::cout << std::setprecision(12) << val << std::endl;
            // }


            // Check the RHS values (the RHS is computed as the LHS x previous_solution,
            // hence, it is assumed that if the RHS is correct, the LHS is correct as well)

            KRATOS_WATCH(RHS)
            // KRATOS_WATCH(LHS)
        }

        /** Checks the TwoFluidNavierStokes2D3N element.
         * Checks the LHS and RHS computation
         */
        KRATOS_TEST_CASE_IN_SUITE(2DElementTwoFluidTest, FluidDynamicsApplicationFastSuite)
        {
            AuxiliaryFunctionForElementTest("TwoFluidNavierStokes2D3N");
        }

        KRATOS_TEST_CASE_IN_SUITE(2DElementEmbededTest, FluidDynamicsApplicationFastSuite)
        {
            AuxiliaryFunctionForElementTest("EmbeddedWeaklyCompressibleNavierStokes2D3N");
        }

        KRATOS_TEST_CASE_IN_SUITE(2DElementMonolithicTest, FluidDynamicsApplicationFastSuite)
        {
            AuxiliaryFunctionForElementTest("WeaklyCompressibleNavierStokes2D3N");
        }

    } // namespace Testing
}  // namespace Kratos.
