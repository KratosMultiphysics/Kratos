// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Gennady Markelov
//

// Project includes
#include "containers/model.h"
#include "structural_mechanics_fast_suite.h"
#include "structural_mechanics_application_variables.h"

#include <utility>
#include <boost/numeric/ublas/assignment.hpp>

namespace Kratos::Testing
{
template<SizeType TNNodes>
void FillModelPartWithVariablesNodesAndDoF(ModelPart& rModelPart, double XEnd, double YEnd)
{
    rModelPart.GetProcessInfo().SetValue(DOMAIN_SIZE, 3);
    rModelPart.AddNodalSolutionStepVariable(DISPLACEMENT);
    rModelPart.AddNodalSolutionStepVariable(ROTATION_Z);

    // Create the test geometry
    rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    rModelPart.CreateNewNode(2, XEnd, YEnd, 0.0);
    if constexpr (TNNodes==3) {
        rModelPart.CreateNewNode(3, XEnd/2, YEnd/2, 0.0);
    }

    for (auto& r_node : rModelPart.Nodes()){
        r_node.AddDof(DISPLACEMENT_X);
        r_node.AddDof(DISPLACEMENT_Y);
        r_node.AddDof(ROTATION_Z);
    }
}

std::vector<ModelPart::IndexType> GetElementNodesFromModelPart(ModelPart& rModelPart)
{
    std::vector<ModelPart::IndexType> element_node_ids;
    const auto& r_nodes = rModelPart.Nodes();
    element_node_ids.reserve(r_nodes.size());
    std::transform(r_nodes.ptr_begin(), r_nodes.ptr_end(), std::back_inserter(element_node_ids),
                   [](const auto& rNodePtr){ return rNodePtr->Id(); });
    return element_node_ids;
}

template<SizeType TNNodes>
void Create2DBeamModel_and_CheckPK2Stress(const std::string & TimoshenkoBeamElementName, const std::vector<double>& rExpectedShearStress, const std::vector<double>& rExpectedBendingMoment)
{
    Model current_model;
    auto &r_model_part = current_model.CreateModelPart("ModelPart",1);
    constexpr double directional_length = 2.0;
    FillModelPartWithVariablesNodesAndDoF<TNNodes>(r_model_part, directional_length, directional_length);

    // Set the element properties
    auto p_elem_prop = r_model_part.CreateNewProperties(0);
    constexpr auto youngs_modulus = 2.0e+06;
    p_elem_prop->SetValue(YOUNG_MODULUS, youngs_modulus);
    p_elem_prop->SetValue(CROSS_AREA, 1.0);
    p_elem_prop->SetValue(I33, 1.0);
    p_elem_prop->SetValue(AREA_EFFECTIVE_Y, 5.0/6.0);

    const auto &r_clone_cl = KratosComponents<ConstitutiveLaw>::Get("TimoshenkoBeamElasticConstitutiveLaw");
    p_elem_prop->SetValue(CONSTITUTIVE_LAW, r_clone_cl.Clone());

    auto element_node_ids = GetElementNodesFromModelPart(r_model_part);
    auto p_element = r_model_part.CreateNewElement(TimoshenkoBeamElementName, 1, element_node_ids, p_elem_prop);
    const auto& r_process_info = r_model_part.GetProcessInfo();
    p_element->Initialize(r_process_info); // Initialize the element to initialize the constitutive law

    //constexpr auto induced_strain = 0.1;
    constexpr auto induced_strain = 0.1;
    p_element->GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT) += ScalarVector(3, induced_strain * directional_length);
    p_element->GetGeometry()[1].FastGetSolutionStepValue(ROTATION_Z) += 0.1;
    if constexpr(TNNodes==3) {
        p_element->GetGeometry()[2].FastGetSolutionStepValue(DISPLACEMENT) += ScalarVector(3, induced_strain * directional_length/2.0);
        p_element->GetGeometry()[2].FastGetSolutionStepValue(ROTATION_Z) += 0.1/2.0;
    }

    std::vector<Vector> stress_vectors;
    p_element->CalculateOnIntegrationPoints(PK2_STRESS_VECTOR, stress_vectors, r_process_info);

    constexpr auto expected_stress = induced_strain * youngs_modulus;
    constexpr auto tolerance = 1.0e-5;
    Vector expected_stress_vector(3);
    expected_stress_vector <<= expected_stress, rExpectedBendingMoment[0], rExpectedShearStress[0];
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(expected_stress_vector, stress_vectors[0], tolerance);
    expected_stress_vector <<= expected_stress, rExpectedBendingMoment[1], rExpectedShearStress[1];
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(expected_stress_vector, stress_vectors[1], tolerance);
    if (stress_vectors.size()>2) {
        expected_stress_vector <<= expected_stress, rExpectedBendingMoment[2], rExpectedShearStress[2];
        KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(expected_stress_vector, stress_vectors[2], tolerance);
    }

    expected_stress_vector <<= expected_stress, rExpectedBendingMoment[0], rExpectedShearStress[0];
    Vector pre_stress(3);
    pre_stress <<= 1.0e5, 1.0e4, 1.0e3;
    p_element->GetProperties().SetValue(BEAM_PRESTRESS_PK2, pre_stress);
    p_element->CalculateOnIntegrationPoints(PK2_STRESS_VECTOR, stress_vectors, r_process_info);
    expected_stress_vector += pre_stress;
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(expected_stress_vector, stress_vectors[0], tolerance);
}

template<SizeType TNNodes>
void Create2DPlaneStrainBeamModel_and_CheckPK2Stress(const std::string& TimoshenkoBeamElementName)
{
    Model current_model;
    auto &r_model_part = current_model.CreateModelPart("ModelPart",1);
    constexpr double directional_length = 2.0;
    FillModelPartWithVariablesNodesAndDoF<TNNodes>(r_model_part, directional_length, 0.0);

    // Set the element properties
    auto p_elem_prop = r_model_part.CreateNewProperties(0);
    constexpr auto youngs_modulus = 2.0e+06;
    constexpr auto poissons_ratio = 0.2;
    constexpr auto thickness = 1.0;
    p_elem_prop->SetValue(YOUNG_MODULUS, youngs_modulus);
    p_elem_prop->SetValue(POISSON_RATIO, poissons_ratio);
    p_elem_prop->SetValue(THICKNESS, thickness);
    constexpr auto effective_shear_thickness = 5.0 / 6.0;
    p_elem_prop->SetValue(THICKNESS_EFFECTIVE_Y, effective_shear_thickness);

    const auto &r_clone_cl = KratosComponents<ConstitutiveLaw>::Get("TimoshenkoBeamPlaneStrainElasticConstitutiveLaw");
    p_elem_prop->SetValue(CONSTITUTIVE_LAW, r_clone_cl.Clone());

    // Create the test element
    auto element_node_ids = GetElementNodesFromModelPart(r_model_part);
    auto p_element = r_model_part.CreateNewElement(TimoshenkoBeamElementName, 1, element_node_ids, p_elem_prop);
    const auto& r_process_info = r_model_part.GetProcessInfo();
    p_element->Initialize(r_process_info); // Initialize the element to initialize the constitutive law
    constexpr auto induced_strain = 0.2;
    constexpr auto induced_rotation = 0.1;
    constexpr auto transverse_bending_displacement = induced_rotation * directional_length / 2.0;
    constexpr auto shear_displacement = 0.03;
    Vector end_displacement(3);
    end_displacement <<= induced_strain * directional_length, transverse_bending_displacement + shear_displacement, 0.0;
    p_element->GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT) += end_displacement;
    p_element->GetGeometry()[1].FastGetSolutionStepValue(ROTATION_Z) += induced_rotation;
    constexpr auto mid_induced_rotation = induced_rotation / 2.0;
    constexpr auto mid_transverse_bending_displacement = mid_induced_rotation * directional_length / 4.0;
    Vector mid_displacement(3);
    mid_displacement <<= induced_strain * directional_length / 2.0, mid_transverse_bending_displacement + shear_displacement / 2.0, 0.0;
    p_element->GetGeometry()[2].FastGetSolutionStepValue(DISPLACEMENT) += mid_displacement;
    p_element->GetGeometry()[2].FastGetSolutionStepValue(ROTATION_Z) += mid_induced_rotation;

    std::vector<Vector> stress_vectors;
    p_element->CalculateOnIntegrationPoints(PK2_STRESS_VECTOR, stress_vectors, r_process_info);
    KRATOS_INFO("Hallo 3") << std::endl;

    constexpr auto expected_distributed_normal_force = induced_strain * youngs_modulus * thickness / (1.0 - poissons_ratio *poissons_ratio);
    constexpr auto expected_distributed_moment =  induced_rotation * youngs_modulus * (thickness * thickness * thickness / 12.0) / (directional_length * ( 1.0 - poissons_ratio *poissons_ratio));
    constexpr auto expected_distributed_shear_force = -(youngs_modulus * effective_shear_thickness / (2.0 * (1.0 + poissons_ratio) )) * (shear_displacement / directional_length);
    constexpr auto tolerance = 1.0e-5;
    Vector expected_stress_vector(5);
    expected_stress_vector <<= expected_distributed_normal_force, expected_distributed_moment, expected_distributed_shear_force, poissons_ratio*expected_distributed_normal_force, poissons_ratio*expected_distributed_moment;
    KRATOS_EXPECT_VECTOR_NEAR(expected_stress_vector, stress_vectors[0], tolerance);
    KRATOS_EXPECT_VECTOR_NEAR(expected_stress_vector, stress_vectors[1], tolerance);

    Vector pre_stress(5);
    pre_stress <<= 1.0e5, 1.0e4, 1.0e3, 2.e4, 2.e3;
    p_element->GetProperties().SetValue(BEAM_PRESTRESS_PK2, pre_stress);
    p_element->CalculateOnIntegrationPoints(PK2_STRESS_VECTOR, stress_vectors, r_process_info);
    expected_stress_vector += pre_stress;
    KRATOS_EXPECT_VECTOR_NEAR(expected_stress_vector, stress_vectors[0], tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(LinearTimoshenkoBeam2D2N_CalculatesPK2Stress, KratosStructuralMechanicsFastSuite)
{
    const std::vector expected_shear_stress{-32608.7, -32608.7, -32608.7};
    const std::vector expected_bending_moment{34989.6, 70710.7, 106432.0};

    Create2DBeamModel_and_CheckPK2Stress<2>("LinearTimoshenkoBeamElement2D2N", expected_shear_stress, expected_bending_moment);
}

KRATOS_TEST_CASE_IN_SUITE(LinearTimoshenkodBeam2D3N_CalculatesPK2Stress, KratosStructuralMechanicsFastSuite)
{
    const std::vector expected_shear_stress{2604.74, -21673.6, -48084.9};
    const std::vector expected_bending_moment{86720.4, 68516.6, 61527.5};
    Create2DBeamModel_and_CheckPK2Stress<3>("LinearTimoshenkoBeamElement2D3N", expected_shear_stress, expected_bending_moment);
}

KRATOS_TEST_CASE_IN_SUITE(LinearTimoshenkodCurvedBeam2D3N_CalculatesPK2Stress, KratosStructuralMechanicsFastSuite)
{
    const std::vector expected_shear_stress{17610.4, 65722.9};
    const std::vector expected_bending_moment{70710.7, 70710.7};
    Create2DBeamModel_and_CheckPK2Stress<3>("LinearTimoshenkoCurvedBeamElement2D3N", expected_shear_stress, expected_bending_moment);
}

KRATOS_TEST_CASE_IN_SUITE(LinearTimoshenkodCurvedBeam2D3N_CalculatesPK2StressPlaneStrain, KratosStructuralMechanicsFastSuite)
{
    Create2DPlaneStrainBeamModel_and_CheckPK2Stress<3>("LinearTimoshenkoCurvedBeamElement2D3N");
}

}
