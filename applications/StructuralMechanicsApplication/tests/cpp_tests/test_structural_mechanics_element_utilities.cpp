// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Philipp Bucher
//

// System includes

// External includes

// Project includes
#include "structural_mechanics_fast_suite.h"
#include "custom_utilities/structural_mechanics_element_utilities.h"
#include "custom_utilities/constitutive_law_utilities.h"
#include "structural_mechanics_application_variables.h"
#include "containers/model.h"

namespace Kratos {
namespace Testing {

// testing the selection for the mass-matrix (has to be consistent if nothing else
// is specified through Properties or ProcessInfo)
KRATOS_TEST_CASE_IN_SUITE(MassMatrixSelection, KratosStructuralMechanicsFastSuite)
{
    Properties aux_props;
    ProcessInfo aux_process_info;

    KRATOS_EXPECT_FALSE(StructuralMechanicsElementUtilities::ComputeLumpedMassMatrix(aux_props, aux_process_info));

    aux_props[COMPUTE_LUMPED_MASS_MATRIX] = true;
    KRATOS_EXPECT_TRUE(StructuralMechanicsElementUtilities::ComputeLumpedMassMatrix(aux_props, aux_process_info));

    aux_props[COMPUTE_LUMPED_MASS_MATRIX] = true;
    aux_process_info[COMPUTE_LUMPED_MASS_MATRIX] = true;
    KRATOS_EXPECT_TRUE(StructuralMechanicsElementUtilities::ComputeLumpedMassMatrix(aux_props, aux_process_info));

    // setting provided through ProcessInfo has priority!
    aux_props[COMPUTE_LUMPED_MASS_MATRIX] = false;
    aux_process_info[COMPUTE_LUMPED_MASS_MATRIX] = true;
    KRATOS_EXPECT_TRUE(StructuralMechanicsElementUtilities::ComputeLumpedMassMatrix(aux_props, aux_process_info));

    // setting provided through ProcessInfo has priority!
    aux_props[COMPUTE_LUMPED_MASS_MATRIX] = true;
    aux_process_info[COMPUTE_LUMPED_MASS_MATRIX] = false;
    KRATOS_EXPECT_FALSE(StructuralMechanicsElementUtilities::ComputeLumpedMassMatrix(aux_props, aux_process_info));
}

// testing the selection for rayleigh-damping (has to be consistent if nothing else
// is specified through Properties or ProcessInfo)
KRATOS_TEST_CASE_IN_SUITE(RayleighDampingSelection, KratosStructuralMechanicsFastSuite)
{
    Properties aux_props_1;
    ProcessInfo aux_process_info_1;

    KRATOS_EXPECT_FALSE(StructuralMechanicsElementUtilities::HasRayleighDamping(aux_props_1, aux_process_info_1));

    const double alpha_0 = StructuralMechanicsElementUtilities::GetRayleighAlpha(aux_props_1, aux_process_info_1);
    KRATOS_EXPECT_DOUBLE_EQ(0.00, alpha_0);

    const double beta_0 = StructuralMechanicsElementUtilities::GetRayleighBeta(aux_props_1, aux_process_info_1);
    KRATOS_EXPECT_DOUBLE_EQ(0.0, beta_0);

    const double val_alpha = 0.01;
    aux_props_1[RAYLEIGH_ALPHA] = val_alpha;
    KRATOS_EXPECT_TRUE(StructuralMechanicsElementUtilities::HasRayleighDamping(aux_props_1, aux_process_info_1));
    aux_process_info_1[RAYLEIGH_ALPHA] = 0.05;
    KRATOS_EXPECT_TRUE(StructuralMechanicsElementUtilities::HasRayleighDamping(aux_props_1, aux_process_info_1));

    Properties aux_props_2;
    ProcessInfo aux_process_info_2;

    const double val_beta = 0.025;
    aux_props_2[RAYLEIGH_BETA] = val_beta;
    KRATOS_EXPECT_TRUE(StructuralMechanicsElementUtilities::HasRayleighDamping(aux_props_2, aux_process_info_2));
    aux_process_info_2[RAYLEIGH_BETA] = 0.06;
    KRATOS_EXPECT_TRUE(StructuralMechanicsElementUtilities::HasRayleighDamping(aux_props_2, aux_process_info_2));

    Properties aux_props_3;
    ProcessInfo aux_process_info_3;
    aux_process_info_3[RAYLEIGH_BETA] = 0.07;
    KRATOS_EXPECT_TRUE(StructuralMechanicsElementUtilities::HasRayleighDamping(aux_props_3, aux_process_info_3));

    // testing if the value defined in the Properties has priority over the value defined in the ProcessInfo
    const double alpha_1 = StructuralMechanicsElementUtilities::GetRayleighAlpha(aux_props_1, aux_process_info_1);
    KRATOS_EXPECT_DOUBLE_EQ(val_alpha, alpha_1);

    const double beta_1 = StructuralMechanicsElementUtilities::GetRayleighBeta(aux_props_2, aux_process_info_2);
    KRATOS_EXPECT_DOUBLE_EQ(val_beta, beta_1);
}

KRATOS_TEST_CASE_IN_SUITE(CalculateElementLength, KratosStructuralMechanicsFastSuite)
{
    // create model part
    Model current_model;

    // Generate a model part with the previous
    ModelPart& model_part = current_model.CreateModelPart("test_model_part");
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT);

    // Fill the model part geometry data
    model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    auto p_node_2 = model_part.CreateNewNode(2, 1.0, 1.0, 1.0);
    model_part.CreateNewNode(3, 0.0, 0.0, 0.0);

    // Set the updated node-position (due to MOVE-MESH)
    p_node_2->X() = 0.5;
    p_node_2->Y() = 0.5;
    p_node_2->Z() = 0.5;

    // Set the DISPLACEMENT
    p_node_2->FastGetSolutionStepValue(DISPLACEMENT_X) = 1.0;
    p_node_2->FastGetSolutionStepValue(DISPLACEMENT_Y) = 1.0;
    p_node_2->FastGetSolutionStepValue(DISPLACEMENT_Z) = 1.0;

    Properties::Pointer p_properties(new Properties(0));
    const std::vector<ModelPart::IndexType> nodes1{1,2};
    const std::vector<ModelPart::IndexType> nodes2{1,3};
    auto p_element_1 = model_part.CreateNewElement("Element3D2N", 1, nodes1, p_properties);
    auto p_element_2 = model_part.CreateNewElement("Element3D2N", 2, nodes2, p_properties);
    auto p_element_3 = model_part.CreateNewElement("Element2D2N", 3, nodes1, p_properties);
    auto p_element_4 = model_part.CreateNewElement("Element2D2N", 4, nodes2, p_properties);

    // length 3D
    KRATOS_EXPECT_DOUBLE_EQ(StructuralMechanicsElementUtilities::CalculateReferenceLength3D2N(*p_element_1), std::sqrt(3.0));
    KRATOS_EXPECT_DOUBLE_EQ(StructuralMechanicsElementUtilities::CalculateCurrentLength3D2N(*p_element_1), std::sqrt(12.0));

    // 0 length throws ERROR  3D
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        StructuralMechanicsElementUtilities::CalculateCurrentLength3D2N(*p_element_2),
        "Error: Element #2 has a current length of zero!");

    // length 2D
    KRATOS_EXPECT_DOUBLE_EQ(StructuralMechanicsElementUtilities::CalculateReferenceLength2D2N(*p_element_3), std::sqrt(2.0));
    KRATOS_EXPECT_DOUBLE_EQ(StructuralMechanicsElementUtilities::CalculateCurrentLength2D2N(*p_element_3), std::sqrt(8.0));

    // 0 length throws ERROR 2D
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        StructuralMechanicsElementUtilities::CalculateCurrentLength2D2N(*p_element_4),
        "Error: Element #4 has a current length of zero!");
}

KRATOS_TEST_CASE_IN_SUITE(CalculateShearPhivalue, KratosStructuralMechanicsFastSuite)
{
    auto props = Properties(1);
    const double L = 0.25;
    props.SetValue(YOUNG_MODULUS, 2.0E9);
    props.SetValue(I33, 1.0e-4);
    props.SetValue(AREA_EFFECTIVE_Y, 0.1);
    const double G = ConstitutiveLawUtilities<3>::CalculateShearModulus(props);
    double phi = StructuralMechanicsElementUtilities::CalculatePhi(props, L);
    KRATOS_EXPECT_DOUBLE_EQ(phi, 12.0 * 2.0e9 * 1.0e-4 / (G * 0.1 * std::pow(L, 2)));

    props.SetValue(AREA_EFFECTIVE_Y, 0.0);
    phi = StructuralMechanicsElementUtilities::CalculatePhi(props, L);
    KRATOS_EXPECT_DOUBLE_EQ(phi, 0.0);
}

} // namespace Testing
} // namespace Kratos
