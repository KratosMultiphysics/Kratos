//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//
//

// System includes
#include <set>

// External includes

// Project includes
#include "testing/testing.h"
#include "includes/model_part.h"
#include "includes/cfd_variables.h"

// Application includes
#include "custom_elements/navier_stokes.h"
#include "custom_utilities/drag_utilities.h"
#include "custom_constitutive/newtonian_2d_law.h"
#include "fluid_dynamics_application_variables.h"

namespace Kratos {
	namespace Testing {

        /**
	     * Auxiliar function to generate a triangular element to be tested.
	     */
        void GenerateTestModelPart(
            ModelPart& rModelPart,
            bool is_embedded = false) {

            // Set buffer size
            rModelPart.SetBufferSize(3);

            // Variables addition
            rModelPart.AddNodalSolutionStepVariable(DENSITY);
            rModelPart.AddNodalSolutionStepVariable(DISTANCE);
            rModelPart.AddNodalSolutionStepVariable(REACTION);
            rModelPart.AddNodalSolutionStepVariable(PRESSURE);
            rModelPart.AddNodalSolutionStepVariable(VELOCITY);
            rModelPart.AddNodalSolutionStepVariable(BODY_FORCE);
            rModelPart.AddNodalSolutionStepVariable(DYNAMIC_TAU);
            rModelPart.AddNodalSolutionStepVariable(MESH_VELOCITY);
            rModelPart.AddNodalSolutionStepVariable(SOUND_VELOCITY);
            rModelPart.AddNodalSolutionStepVariable(DYNAMIC_VISCOSITY);

            // Process info creation
            double delta_time = 0.1;
            rModelPart.GetProcessInfo().SetValue(DOMAIN_SIZE, 2);
            rModelPart.GetProcessInfo().SetValue(DYNAMIC_TAU, 1.0);
            rModelPart.GetProcessInfo().SetValue(SOUND_VELOCITY, 1.0e+12);
            rModelPart.GetProcessInfo().SetValue(DELTA_TIME, delta_time);
            Vector bdf_coefs(3);
            bdf_coefs[0] = 3.0 / (2.0 * delta_time);
            bdf_coefs[1] = -2.0 / delta_time;
            bdf_coefs[2] = 0.5 * delta_time;
            rModelPart.GetProcessInfo().SetValue(BDF_COEFFICIENTS, bdf_coefs);

            // Set the element properties
            Newtonian2DLaw::Pointer p_cons_law(new Newtonian2DLaw());
            Properties::Pointer p_elem_prop = rModelPart.pGetProperties(0);
            p_elem_prop->SetValue(DENSITY, 1.0e+00);
            p_elem_prop->SetValue(DYNAMIC_VISCOSITY, 3.0e-02);
            p_elem_prop->SetValue(CONSTITUTIVE_LAW, p_cons_law);

            // Element creation
            rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
            rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
            rModelPart.CreateNewNode(3, 0.0, 1.0, 0.0);
            std::vector<ModelPart::IndexType> cond_nodes {1, 2};
            std::vector<ModelPart::IndexType> elem_nodes {1, 2, 3};
            if (is_embedded){
                rModelPart.CreateNewElement("EmbeddedNavierStokes2D3N", 1, elem_nodes, p_elem_prop);
            } else {
                rModelPart.CreateNewElement("NavierStokes2D3N", 1, elem_nodes, p_elem_prop);
            }
            rModelPart.CreateNewCondition("NavierStokesWallCondition2D2N", 1, cond_nodes, p_elem_prop);

            // Set the drag computation submodelpart
            ModelPart* p_sub_model_part = &rModelPart.CreateSubModelPart("DragModelPart");
            std::vector<ModelPart::IndexType> sub_model_part_nodes = {1, 2};
            std::vector<ModelPart::IndexType> sub_model_part_conds = {1};
            p_sub_model_part->AddNodes(sub_model_part_nodes);
            p_sub_model_part->AddConditions(sub_model_part_conds);

            // Add DOFs
            auto nodes_begin = rModelPart.NodesBegin();
            for (unsigned int i_node = 0; i_node < rModelPart.NumberOfNodes(); ++i_node){
                auto it_node = nodes_begin + i_node;
                it_node->AddDof(VELOCITY_X,REACTION_X);
                it_node->AddDof(VELOCITY_Y,REACTION_Y);
                it_node->AddDof(VELOCITY_Z,REACTION_Z);
                it_node->AddDof(PRESSURE,REACTION_WATER_PRESSURE);
            }

            // Set the VELOCITY and PRESSURE nodal values
            const double p_1 = 1.5;
            const double p_2 = 1.0;
            const double p_3 = 0.5;
            array_1d<double, 3> v_1 = ZeroVector(3);
            array_1d<double, 3> v_2 = ZeroVector(3);
            array_1d<double, 3> v_3 = ZeroVector(3);
            v_1[0] = 1.0;
            v_2[0] = 2.0;
            v_3[0] = 3.0;
            v_3[1] = 0.5;
            (rModelPart.GetNode(1)).GetSolutionStepValue(VELOCITY) = v_1;
            (rModelPart.GetNode(1)).GetSolutionStepValue(PRESSURE) = p_1;
            (rModelPart.GetNode(2)).GetSolutionStepValue(VELOCITY) = v_2;
            (rModelPart.GetNode(2)).GetSolutionStepValue(PRESSURE) = p_2;
            (rModelPart.GetNode(3)).GetSolutionStepValue(VELOCITY) = v_3;
            (rModelPart.GetNode(3)).GetSolutionStepValue(PRESSURE) = p_3;

            // Set the DENSITY and DYNAMIC_VISCOSITY nodal values
            for (ModelPart::NodeIterator it_node = rModelPart.NodesBegin(); it_node < rModelPart.NodesEnd(); ++it_node){
                it_node->FastGetSolutionStepValue(DENSITY) = p_elem_prop->GetValue(DENSITY);
                it_node->FastGetSolutionStepValue(DYNAMIC_VISCOSITY) = p_elem_prop->GetValue(DYNAMIC_VISCOSITY);
            }

            // If proceeds, set the DISTANCE function
            if (is_embedded){
                (rModelPart.GetNode(1)).GetSolutionStepValue(DISTANCE) = 1.0;
                (rModelPart.GetNode(2)).GetSolutionStepValue(DISTANCE) = 1.0;
                (rModelPart.GetNode(3)).GetSolutionStepValue(DISTANCE) = -1.0;
            }

        }

	    /**
	     * Checks the body fitted drag computation utility.
	     */
	    KRATOS_TEST_CASE_IN_SUITE(ComputeBodyFittedDrag, FluidDynamicsApplicationFastSuite)
		{
            // Create a test element inside a modelpart
            Model model;
            ModelPart& model_part = model.CreateModelPart("Main", 3);
            GenerateTestModelPart(model_part);
            Element::Pointer p_element = model_part.pGetElement(1);

            // Initialize the fluid element
            p_element->Initialize();

            // Set the reaction values manually. Note that the body fitted drag utilities assume
            // that the REACTION has been already computed. Since this is assumed to be done by
            // the builder and solver, which is out of the scope of this test, we do it manually.
            model_part.GetNode(1).GetDof(VELOCITY_X).GetSolutionStepReactionValue() = 5.0;
            model_part.GetNode(1).GetDof(VELOCITY_Y).GetSolutionStepReactionValue() = 10.0;
            model_part.GetNode(2).GetDof(VELOCITY_X).GetSolutionStepReactionValue() = -20.0;
            model_part.GetNode(2).GetDof(VELOCITY_Y).GetSolutionStepReactionValue() = -40.0;

            // Call the body fitted drag utility
            DragUtilities drag_utilities;
            array_1d<double, 3> drag_force = drag_utilities.CalculateBodyFittedDrag(model_part.GetSubModelPart("DragModelPart"));

            // Check computed values
            KRATOS_CHECK_NEAR(drag_force[0], 15.0, 1e-6);
            KRATOS_CHECK_NEAR(drag_force[1], 30.0, 1e-6);
            KRATOS_CHECK_NEAR(drag_force[2], 0.0, 1e-6);
	    }

	    /**
	     * Checks the embedded drag computation utility.
	     */
	    KRATOS_TEST_CASE_IN_SUITE(ComputeEmbeddedDrag, FluidDynamicsApplicationFastSuite)
		{
            bool is_embedded = true;

            // Create a test element inside a modelpart
            Model model;
            ModelPart& model_part = model.CreateModelPart("Main", 3);
            GenerateTestModelPart(model_part, is_embedded);
            Element::Pointer p_element = model_part.pGetElement(1);

            // Initialize the fluid element
            p_element->Initialize();

            // Call the embedded drag utility
            DragUtilities drag_utilities;
            array_1d<double, 3> drag_force = drag_utilities.CalculateEmbeddedDrag(model_part);

            // Check computed values
            KRATOS_CHECK_NEAR(drag_force[0], -0.03, 1e-2);
            KRATOS_CHECK_NEAR(drag_force[1], 0.4375, 1e-4);
            KRATOS_CHECK_NEAR(drag_force[2], 0.0, 1e-6);
	    }
    } // namespace Testing
}  // namespace Kratos.
