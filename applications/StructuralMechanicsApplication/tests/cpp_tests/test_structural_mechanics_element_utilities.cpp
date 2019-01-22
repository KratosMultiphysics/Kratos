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

} // namespace Testing
} // namespace Kratos
