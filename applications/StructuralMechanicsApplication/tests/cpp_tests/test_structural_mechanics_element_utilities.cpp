// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Philipp Bucher
//

// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "custom_utilities/structural_mechanics_element_utilities.h"
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

    KRATOS_CHECK_IS_FALSE(StructuralMechanicsElementUtilities::ComputeLumpedMassMatrix(aux_props, aux_process_info));

    aux_props[COMPUTE_LUMPED_MASS_MATRIX] = true;
    KRATOS_CHECK(StructuralMechanicsElementUtilities::ComputeLumpedMassMatrix(aux_props, aux_process_info));

    aux_props[COMPUTE_LUMPED_MASS_MATRIX] = true;
    aux_process_info[COMPUTE_LUMPED_MASS_MATRIX] = true;
    KRATOS_CHECK(StructuralMechanicsElementUtilities::ComputeLumpedMassMatrix(aux_props, aux_process_info));

    // setting provided through ProcessInfo has priority!
    aux_props[COMPUTE_LUMPED_MASS_MATRIX] = false;
    aux_process_info[COMPUTE_LUMPED_MASS_MATRIX] = true;
    KRATOS_CHECK(StructuralMechanicsElementUtilities::ComputeLumpedMassMatrix(aux_props, aux_process_info));

    // setting provided through ProcessInfo has priority!
    aux_props[COMPUTE_LUMPED_MASS_MATRIX] = true;
    aux_process_info[COMPUTE_LUMPED_MASS_MATRIX] = false;
    KRATOS_CHECK_IS_FALSE(StructuralMechanicsElementUtilities::ComputeLumpedMassMatrix(aux_props, aux_process_info));
}

// testing the selection for rayleigh-damping (has to be consistent if nothing else
// is specified through Properties or ProcessInfo)
KRATOS_TEST_CASE_IN_SUITE(RayleighDampingSelection, KratosStructuralMechanicsFastSuite)
{
    Properties aux_props_1;
    ProcessInfo aux_process_info_1;

    KRATOS_CHECK_IS_FALSE(StructuralMechanicsElementUtilities::HasRayleighDamping(aux_props_1, aux_process_info_1));

    const double alpha_0 = StructuralMechanicsElementUtilities::GetRayleighAlpha(aux_props_1, aux_process_info_1);
    KRATOS_CHECK_DOUBLE_EQUAL(0.00, alpha_0);

    const double beta_0 = StructuralMechanicsElementUtilities::GetRayleighBeta(aux_props_1, aux_process_info_1);
    KRATOS_CHECK_DOUBLE_EQUAL(0.0, beta_0);

    const double val_alpha = 0.01;
    aux_props_1[RAYLEIGH_ALPHA] = val_alpha;
    KRATOS_CHECK(StructuralMechanicsElementUtilities::HasRayleighDamping(aux_props_1, aux_process_info_1));
    aux_process_info_1[RAYLEIGH_ALPHA] = 0.05;
    KRATOS_CHECK(StructuralMechanicsElementUtilities::HasRayleighDamping(aux_props_1, aux_process_info_1));

    Properties aux_props_2;
    ProcessInfo aux_process_info_2;

    const double val_beta = 0.025;
    aux_props_2[RAYLEIGH_BETA] = val_beta;
    KRATOS_CHECK(StructuralMechanicsElementUtilities::HasRayleighDamping(aux_props_2, aux_process_info_2));
    aux_process_info_2[RAYLEIGH_BETA] = 0.06;
    KRATOS_CHECK(StructuralMechanicsElementUtilities::HasRayleighDamping(aux_props_2, aux_process_info_2));

    Properties aux_props_3;
    ProcessInfo aux_process_info_3;
    aux_process_info_3[RAYLEIGH_BETA] = 0.07;
    KRATOS_CHECK(StructuralMechanicsElementUtilities::HasRayleighDamping(aux_props_3, aux_process_info_3));

    // testing if the value defined in the Properties has priority over the value defined in the ProcessInfo
    const double alpha_1 = StructuralMechanicsElementUtilities::GetRayleighAlpha(aux_props_1, aux_process_info_1);
    KRATOS_CHECK_DOUBLE_EQUAL(val_alpha, alpha_1);

    const double beta_1 = StructuralMechanicsElementUtilities::GetRayleighBeta(aux_props_2, aux_process_info_2);
    KRATOS_CHECK_DOUBLE_EQUAL(val_beta, beta_1);
}

KRATOS_TEST_CASE_IN_SUITE(CalculateElementLength, KratosStructuralMechanicsFastSuite)
{
    // create model part
    Model current_model;

    // Generate a model part with the previous
    ModelPart& model_part = current_model.CreateModelPart("Tetrahedra");
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT);

    // Fill the model part geometry data
    model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    model_part.CreateNewNode(2, 1.0, 1.0, 1.0);
    model_part.CreateNewNode(3, 0.0, 0.0, 0.0);


    // Set the DISTANCE field
    model_part.Nodes()[2].FastGetSolutionStepValue(DISPLACEMENT_X) = 1.0;
    model_part.Nodes()[2].FastGetSolutionStepValue(DISPLACEMENT_Y) = 1.0;
    model_part.Nodes()[2].FastGetSolutionStepValue(DISPLACEMENT_Z) = 1.0;

    Properties::Pointer p_properties(new Properties(0));
    std::vector<ModelPart::IndexType> nodes1{1,2};
    std::vector<ModelPart::IndexType> nodes2{1,3};
    model_part.CreateNewElement("Element3D2N", 1, nodes1, p_properties);
    model_part.CreateNewElement("Element3D2N", 2, nodes2, p_properties);
    model_part.CreateNewElement("Element2D2N", 3, nodes1, p_properties);
    model_part.CreateNewElement("Element2D2N", 4, nodes2, p_properties);


    // length 3D
    Element& r_element_1 = model_part.GetElement(1);
    KRATOS_CHECK_DOUBLE_EQUAL(StructuralMechanicsElementUtilities::CalculateReferenceLength3D2N(r_element_1), sqrt(3.0));
    KRATOS_CHECK_DOUBLE_EQUAL(StructuralMechanicsElementUtilities::CalculateCurrentLength3D2N(r_element_1), sqrt(12.0));

    // 0 length throws ERROR  3D
    Element& r_element_2 = model_part.GetElement(2);
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        StructuralMechanicsElementUtilities::CalculateReferenceLength3D2N(r_element_2),
        "Error: Reference length of element 2 ~0");
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        StructuralMechanicsElementUtilities::CalculateCurrentLength3D2N(r_element_2),
        "Error: Current length of element 2 ~0");

    // length 2D
    Element& r_element_3 = model_part.GetElement(3);
    KRATOS_CHECK_DOUBLE_EQUAL(StructuralMechanicsElementUtilities::CalculateReferenceLength2D2N(r_element_3), sqrt(2.0));
    KRATOS_CHECK_DOUBLE_EQUAL(StructuralMechanicsElementUtilities::CalculateCurrentLength2D2N(r_element_3), sqrt(8.0));

    // 0 length throws ERROR 2D
    Element& r_element_4 = model_part.GetElement(4);
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        StructuralMechanicsElementUtilities::CalculateReferenceLength2D2N(r_element_4),
        "Error: Reference length of element 4 ~0");
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        StructuralMechanicsElementUtilities::CalculateCurrentLength2D2N(r_element_4),
        "Error: Current length of element 4 ~0");
}

} // namespace Testing
} // namespace Kratos
