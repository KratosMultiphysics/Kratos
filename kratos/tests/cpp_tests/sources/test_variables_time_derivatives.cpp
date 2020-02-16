//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//

// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "includes/variables.h"
#include "includes/variables_time_derivatives.h"


namespace Kratos {
namespace Testing {

KRATOS_TEST_CASE_IN_SUITE(VariableDerivativeHasMethods, KratosCoreFastSuite)
{
    const bool checktrue = VariablesTimeDerivatives<Variable<array_1d<double, 3>>>::Has(DISPLACEMENT);
    KRATOS_CHECK(checktrue);
    const bool checkfalse = VariablesTimeDerivatives<Variable<array_1d<double, 3>>>::Has(VECTOR_LAGRANGE_MULTIPLIER);
    KRATOS_CHECK_IS_FALSE(checkfalse);
}

KRATOS_TEST_CASE_IN_SUITE(VariableDerivativeGetMethods, KratosCoreFastSuite)
{
    const auto& r_velocity = VariablesTimeDerivatives<Variable<array_1d<double, 3>>>::GetFirstDerivative(DISPLACEMENT);
    KRATOS_CHECK(r_velocity.Name() == "VELOCITY");
    const auto& r_acceleration = VariablesTimeDerivatives<Variable<array_1d<double, 3>>>::GetSecondDerivative(DISPLACEMENT);
    KRATOS_CHECK(r_acceleration.Name() == "ACCELERATION");
}

}
}  // namespace Kratos.
