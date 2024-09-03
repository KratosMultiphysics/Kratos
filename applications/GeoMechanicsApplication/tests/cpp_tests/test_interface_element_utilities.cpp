// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Richard Faasse
//

#include "geo_mechanics_fast_suite.h"

#include <custom_utilities/interface_element_utilities.hpp>

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(InterfaceElementUtilities_RotationMatrixDoesNotChangeLength, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    Vector vector = ScalarVector{2, 1.0};

    // Act
    Matrix rotation_matrix = InterfaceElementUtilities::Calculate2DRotationMatrix();
    Vector rotated_vector = prod(rotation_matrix, vector);

    // Assert
    KRATOS_CHECK_NEAR(norm_2(rotated_vector), norm_2(vector), 1e-6);
}


}