// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Riccardo Rossi
//

// Project includes
#include "containers/model.h"
#include "testing/testing.h"
#include "structural_mechanics_application_variables.h"
#include "custom_elements/total_lagrangian.h"
#include "custom_processes/set_cartesian_local_axes_process.h"
#include "custom_processes/set_cylindrical_local_axes_process.h"
#include "custom_processes/set_spherical_local_axes_process.h"

namespace Kratos
{
namespace Testing
{

    KRATOS_TEST_CASE_IN_SUITE(RotatedElementCartesian2D3N, KratosStructuralMechanicsFastSuite)
    {
        Model current_model;
        auto &r_model_part = current_model.CreateModelPart("ModelPart",1);
        r_model_part.GetProcessInfo().SetValue(DOMAIN_SIZE, 2);

        r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
        r_model_part.AddNodalSolutionStepVariable(VOLUME_ACCELERATION);

        // Set the element properties
        auto p_elem_prop = r_model_part.CreateNewProperties(0);
        p_elem_prop->SetValue(YOUNG_MODULUS, 2.0e+06);
        p_elem_prop->SetValue(POISSON_RATIO, 0.3);
        p_elem_prop->SetValue(THICKNESS, 0.01);
        const auto &r_clone_cl = KratosComponents<ConstitutiveLaw>::Get("LinearElasticPlaneStress2DLaw");
        p_elem_prop->SetValue(CONSTITUTIVE_LAW, r_clone_cl.Clone());

        // Constants for the computation of the stress
        const double E = p_elem_prop->GetValue(YOUNG_MODULUS);
        const double NU = p_elem_prop->GetValue(POISSON_RATIO);
        const double c1 = E / (1.00 - NU * NU);
        const double c2 = c1 * NU;
        const double c3 = 0.5* E / (1 + NU);


        // Create the test element
        auto p_node_1 = r_model_part.CreateNewNode(1, 0.0 , 0.0 , 0.0);
        auto p_node_2 = r_model_part.CreateNewNode(2, 1.0 , 0.0 , 0.0);
        auto p_node_3 = r_model_part.CreateNewNode(3, 0.0 , 1.0 , 0.0);

        for (auto& r_node : r_model_part.Nodes()){
            r_node.AddDof(DISPLACEMENT_X);
            r_node.AddDof(DISPLACEMENT_Y);
            r_node.AddDof(DISPLACEMENT_Z);
        }

        std::vector<ModelPart::IndexType> element_nodes {1,2,3};
        auto p_element = r_model_part.CreateNewElement("SmallDisplacementElement2D3N", 1, element_nodes, p_elem_prop);

        // Now we create the rotation process
        Parameters parameters = Parameters(R"(
        {
            "cartesian_local_axis"          : [[0.0,1.0,0.0],[1.0,0.0,0.0]]
        })");
        auto &r_cartesian_orientation_process = SetCartesianLocalAxesProcess(r_model_part, parameters);
        r_cartesian_orientation_process.ExecuteInitialize();

        array_1d<double, 3> local_axis_1 = ZeroVector(3);
        local_axis_1[1] = 1.0;
        const array_1d<double, 3> computed_local_axis_1 = p_element->GetValue(LOCAL_AXIS_1);
        KRATOS_CHECK_VECTOR_EQUAL(computed_local_axis_1, local_axis_1);
    }




}
}
