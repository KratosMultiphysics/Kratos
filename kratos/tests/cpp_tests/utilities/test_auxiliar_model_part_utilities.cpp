//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:   BSD License
//      Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "testing/testing.h"
#include "includes/model_part.h"
#include "utilities/auxiliar_model_part_utilities.h"

// Utilities
#include "utilities/cpp_tests_utilities.h"

namespace Kratos {
namespace Testing {

KRATOS_TEST_CASE_IN_SUITE(AuxiliarModelPartUtilities_GetData_Node_historical, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& this_model_part = current_model.CreateModelPart("Main");
    this_model_part.SetBufferSize(2);

    CppTestsUtilities::Create2DGeometry(this_model_part, "Element2D3N");

    // ...
    // check other tests and also kratos/includes/checks.h

    std::cout << "Welcome to testing!" << std::endl;

    auto data = AuxiliarModelPartUtilities(this_model_part).GetVariableData(DISPLACEMENT, DataLocation::NodeHistorical);
}

} // namespace Testing
} // namespace Kratos.
