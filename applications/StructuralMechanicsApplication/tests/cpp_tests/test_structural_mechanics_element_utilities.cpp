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

} // namespace Testing
} // namespace Kratos
