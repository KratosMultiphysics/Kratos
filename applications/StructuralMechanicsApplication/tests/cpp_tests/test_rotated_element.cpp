// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Alejandro Cornejo
//

// Project includes
#include "containers/model.h"
#include "structural_mechanics_fast_suite.h"
#include "structural_mechanics_application_variables.h"
#include "custom_elements/solid_elements/total_lagrangian.h"
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

        // Create the test element
        auto p_node_1 = r_model_part.CreateNewNode(1, 0.0 , 0.0 , 0.0);
        auto p_node_2 = r_model_part.CreateNewNode(2, 1.0 , 0.0 , 0.0);
        auto p_node_3 = r_model_part.CreateNewNode(3, 0.0 , 1.0 , 0.0);

        std::vector<ModelPart::IndexType> element_nodes {1,2,3};
        auto p_element = r_model_part.CreateNewElement("SmallDisplacementElement2D3N", 1, element_nodes, p_elem_prop);

        // Now we create the rotation process
        Parameters parameters = Parameters(R"(
        {
            "cartesian_local_axis"          : [0.0,1.0,0.0]
        })");
        SetCartesianLocalAxesProcess(r_model_part, parameters).ExecuteInitialize();

        array_1d<double, 3> local_axis_1 = ZeroVector(3);
        local_axis_1[1] = 1.0;
        const array_1d<double, 3>& r_computed_local_axis_1 = p_element->GetValue(LOCAL_AXIS_1);
        KRATOS_EXPECT_VECTOR_EQ(r_computed_local_axis_1, local_axis_1);
    }

    KRATOS_TEST_CASE_IN_SUITE(RotatedElementCartesian3D4N, KratosStructuralMechanicsFastSuite)
    {
        Model current_model;
        auto &r_model_part = current_model.CreateModelPart("ModelPart",1);
        r_model_part.GetProcessInfo().SetValue(DOMAIN_SIZE, 3);

        r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
        r_model_part.AddNodalSolutionStepVariable(VOLUME_ACCELERATION);

        // Set the element properties
        auto p_elem_prop = r_model_part.CreateNewProperties(0);
        p_elem_prop->SetValue(YOUNG_MODULUS, 2.0e+06);
        p_elem_prop->SetValue(POISSON_RATIO, 0.3);
        const auto &r_clone_cl = KratosComponents<ConstitutiveLaw>::Get("LinearElastic3DLaw");
        p_elem_prop->SetValue(CONSTITUTIVE_LAW, r_clone_cl.Clone());

        // Create the test element
        auto p_node_1 = r_model_part.CreateNewNode(1, 0.0 , 0.0 , 0.0);
        auto p_node_2 = r_model_part.CreateNewNode(2, 1.0 , 0.0 , 0.0);
        auto p_node_3 = r_model_part.CreateNewNode(3, 0.0 , 1.0 , 0.0);
        auto p_node_4 = r_model_part.CreateNewNode(4, 0.0 , 0.0 , 1.0);

        std::vector<ModelPart::IndexType> element_nodes {1,2,3,4};
        auto p_element = r_model_part.CreateNewElement("SmallDisplacementElement3D4N", 1, element_nodes, p_elem_prop);

        // Now we create the rotation process
        Parameters parameters = Parameters(R"(
        {
            "cartesian_local_axis" : [[0.0,1.0,0.0],[0.0,0.0,1.0]]
        })");
        SetCartesianLocalAxesProcess(r_model_part, parameters).ExecuteInitialize();

        array_1d<double, 3> local_axis_1 = ZeroVector(3);
        array_1d<double, 3> local_axis_2 = ZeroVector(3);
        local_axis_1[1] = 1.0;
        local_axis_2[2] = 1.0;
        const array_1d<double, 3>& r_computed_local_axis_1 = p_element->GetValue(LOCAL_AXIS_1);
        const array_1d<double, 3>& r_computed_local_axis_2 = p_element->GetValue(LOCAL_AXIS_2);
        KRATOS_EXPECT_VECTOR_EQ(r_computed_local_axis_1, local_axis_1);
        KRATOS_EXPECT_VECTOR_EQ(r_computed_local_axis_2, local_axis_2);
    }

    KRATOS_TEST_CASE_IN_SUITE(RotatedElementCylindrical2D3N, KratosStructuralMechanicsFastSuite)
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

        // Create the test element
        auto p_node_1 = r_model_part.CreateNewNode(1, 0.0 , 0.0 , 0.0);
        auto p_node_2 = r_model_part.CreateNewNode(2, 1.0 , 0.0 , 0.0);
        auto p_node_3 = r_model_part.CreateNewNode(3, 0.0 , 1.0 , 0.0);

        std::vector<ModelPart::IndexType> element_nodes {1,2,3};
        auto p_element = r_model_part.CreateNewElement("SmallDisplacementElement2D3N", 1, element_nodes, p_elem_prop);

        // Now we create the rotation process
        Parameters parameters = Parameters(R"(
        {
            "cylindrical_generatrix_axis"   : [0.0,0.0,1.0],
            "cylindrical_generatrix_point"  : [0.0,0.0,0.0]
        })");
        SetCylindricalLocalAxesProcess(r_model_part, parameters).ExecuteInitialize();

        array_1d<double, 3> local_axis_1 = ZeroVector(3);
        const double number = 0.5 * sqrt(2.0);
        local_axis_1[0] = number;
        local_axis_1[1] = number;
        const array_1d<double, 3>& r_computed_local_axis_1 = p_element->GetValue(LOCAL_AXIS_1);
        KRATOS_EXPECT_VECTOR_EQ(r_computed_local_axis_1, local_axis_1);
    }

    KRATOS_TEST_CASE_IN_SUITE(RotatedElementCylindrical3D4N, KratosStructuralMechanicsFastSuite)
    {
        Model current_model;
        auto &r_model_part = current_model.CreateModelPart("ModelPart",1);
        r_model_part.GetProcessInfo().SetValue(DOMAIN_SIZE, 3);

        r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
        r_model_part.AddNodalSolutionStepVariable(VOLUME_ACCELERATION);

        // Set the element properties
        auto p_elem_prop = r_model_part.CreateNewProperties(0);
        p_elem_prop->SetValue(YOUNG_MODULUS, 2.0e+06);
        p_elem_prop->SetValue(POISSON_RATIO, 0.3);
        const auto &r_clone_cl = KratosComponents<ConstitutiveLaw>::Get("LinearElastic3DLaw");
        p_elem_prop->SetValue(CONSTITUTIVE_LAW, r_clone_cl.Clone());

        // Create the test element
        auto p_node_1 = r_model_part.CreateNewNode(1, 0.0 , 0.0 , 0.0);
        auto p_node_2 = r_model_part.CreateNewNode(2, 1.0 , 0.0 , 0.0);
        auto p_node_3 = r_model_part.CreateNewNode(3, 0.0 , 1.0 , 0.0);
        auto p_node_4 = r_model_part.CreateNewNode(4, 0.0 , 0.0 , 1.0);

        std::vector<ModelPart::IndexType> element_nodes {1,2,3,4};
        auto p_element = r_model_part.CreateNewElement("SmallDisplacementElement3D4N", 1, element_nodes, p_elem_prop);

        // Now we create the rotation process
        Parameters parameters = Parameters(R"(
        {
            "cylindrical_generatrix_axis"   : [0.0,0.0,1.0],
            "cylindrical_generatrix_point"  : [0.0,0.0,0.0]
        })");
        SetCylindricalLocalAxesProcess(r_model_part, parameters).ExecuteInitialize();

        array_1d<double, 3> local_axis_1 = ZeroVector(3);
        array_1d<double, 3> local_axis_2 = ZeroVector(3);
        const double number = 0.5 * sqrt(2.0);
        local_axis_1[0] = number;
        local_axis_1[1] = number;
        local_axis_2[2] = 1.0;
        const array_1d<double, 3>& r_computed_local_axis_1 = p_element->GetValue(LOCAL_AXIS_1);
        const array_1d<double, 3>& r_computed_local_axis_2 = p_element->GetValue(LOCAL_AXIS_2);
        KRATOS_EXPECT_VECTOR_EQ(r_computed_local_axis_1, local_axis_1);
        KRATOS_EXPECT_VECTOR_EQ(r_computed_local_axis_2, local_axis_2);
    }

    KRATOS_TEST_CASE_IN_SUITE(RotatedElementSpherical2D3N, KratosStructuralMechanicsFastSuite)
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

        // Create the test element
        auto p_node_1 = r_model_part.CreateNewNode(1, 0.0 , 0.0 , 0.0);
        auto p_node_2 = r_model_part.CreateNewNode(2, 1.0 , 0.0 , 0.0);
        auto p_node_3 = r_model_part.CreateNewNode(3, 0.0 , 1.0 , 0.0);

        std::vector<ModelPart::IndexType> element_nodes {1,2,3};
        auto p_element = r_model_part.CreateNewElement("SmallDisplacementElement2D3N", 1, element_nodes, p_elem_prop);

        // Now we create the rotation process
        Parameters parameters = Parameters(R"(
        {
            "spherical_reference_axis"   : [0.0,0.0,1.0],
            "spherical_central_point"    : [0.0,0.0,0.0]
        })");
        SetSphericalLocalAxesProcess(r_model_part, parameters).ExecuteInitialize();

        array_1d<double, 3> local_axis_1 = ZeroVector(3);
        const double number = 0.5 * sqrt(2.0);
        local_axis_1[0] = number;
        local_axis_1[1] = number;
        const array_1d<double, 3>& r_computed_local_axis_1 = p_element->GetValue(LOCAL_AXIS_1);
        KRATOS_EXPECT_VECTOR_EQ(r_computed_local_axis_1, local_axis_1);
    }

    KRATOS_TEST_CASE_IN_SUITE(RotatedElementSpherical3D4N, KratosStructuralMechanicsFastSuite)
    {
        Model current_model;
        auto &r_model_part = current_model.CreateModelPart("ModelPart",1);
        r_model_part.GetProcessInfo().SetValue(DOMAIN_SIZE, 3);

        r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
        r_model_part.AddNodalSolutionStepVariable(VOLUME_ACCELERATION);

        // Set the element properties
        auto p_elem_prop = r_model_part.CreateNewProperties(0);
        p_elem_prop->SetValue(YOUNG_MODULUS, 2.0e+06);
        p_elem_prop->SetValue(POISSON_RATIO, 0.3);
        const auto &r_clone_cl = KratosComponents<ConstitutiveLaw>::Get("LinearElastic3DLaw");
        p_elem_prop->SetValue(CONSTITUTIVE_LAW, r_clone_cl.Clone());

        // Create the test element
        auto p_node_1 = r_model_part.CreateNewNode(1, 0.0 , 0.0 , 0.0);
        auto p_node_2 = r_model_part.CreateNewNode(2, 1.0 , 0.0 , 0.0);
        auto p_node_3 = r_model_part.CreateNewNode(3, 0.0 , 1.0 , 0.0);
        auto p_node_4 = r_model_part.CreateNewNode(4, 0.0 , 0.0 , 1.0);

        std::vector<ModelPart::IndexType> element_nodes {1,2,3,4};
        auto p_element = r_model_part.CreateNewElement("SmallDisplacementElement3D4N", 1, element_nodes, p_elem_prop);

        // Now we create the rotation process
        Parameters parameters = Parameters(R"(
        {
            "spherical_reference_axis"   : [0.0,0.0,1.0],
            "spherical_central_point"    : [0.0,0.0,0.0]
        })");
        SetSphericalLocalAxesProcess(r_model_part, parameters).ExecuteInitialize();

        array_1d<double, 3> local_axis_1 = ZeroVector(3);
        array_1d<double, 3> local_axis_2 = ZeroVector(3);
        const double number = 1.0 / sqrt(3.0);
        const double number2 = 0.5 * sqrt(2.0);
        local_axis_1[0] = number;
        local_axis_1[1] = number;
        local_axis_1[2] = number;
        local_axis_2[0] = number2;
        local_axis_2[1] = -number2;
        const array_1d<double, 3>& r_computed_local_axis_1 = p_element->GetValue(LOCAL_AXIS_1);
        const array_1d<double, 3>& r_computed_local_axis_2 = p_element->GetValue(LOCAL_AXIS_2);
        KRATOS_EXPECT_VECTOR_EQ(r_computed_local_axis_1, local_axis_1);
        KRATOS_EXPECT_VECTOR_EQ(r_computed_local_axis_2, local_axis_2);
    }
}
}
