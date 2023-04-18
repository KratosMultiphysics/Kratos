//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// System includes
#include <functional>
#include <sstream>

// External includes

// Project includes
#include "containers/model.h"
#include "containers/global_pointers_vector.h"
#include "includes/checks.h"
#include "includes/model_part.h"
#include "includes/ublas_interface.h"
#include "includes/cfd_variables.h"

// Application includes
#include "custom_utilities/fluid_test_utilities.h"
#include "custom_utilities/rans_calculation_utilities.h"
#include "custom_utilities/rans_variable_utilities.h"

// Include base h
#include "test_utilities.h"

namespace Kratos
{
namespace RansApplicationTestUtilities
{

ModelPart& CreateScalarVariableTestModelPart(
    Model& rModel,
    const std::string& rElementName,
    const std::string& rConditionName,
    const std::function<void(Properties&)>& rSetElementProperties,
    const std::function<void(Properties&)>& rSetConditionProperties,
    const std::function<void(ModelPart& rModelPart)>& rAddNodalSolutionStepVariablesFuncion,
    const Variable<double>& rDofVariable,
    const int BufferSize,
    const bool DoInitializeElements,
    const bool DoInitializeConditions)
{
    auto& r_model_part = FluidTestUtilities::CreateTestModelPart(
        rModel, "test", rElementName, rConditionName, rSetElementProperties,
        rSetConditionProperties, rAddNodalSolutionStepVariablesFuncion,
        [&rDofVariable](ModelPart::NodeType& rNode) {
            rNode.AddDof(rDofVariable).SetEquationId(rNode.Id());
        },
        BufferSize);

    if (DoInitializeElements) {
        r_model_part.Elements().front().Initialize(r_model_part.GetProcessInfo());
    }

    if (DoInitializeConditions) {
        r_model_part.Conditions().front().Initialize(r_model_part.GetProcessInfo());
    }

    return r_model_part;
}

void CheckElementsAndConditions(
    const ModelPart& rModelPart)
{
    const auto& r_process_info = rModelPart.GetProcessInfo();
    for (const auto& r_element : rModelPart.Elements()) {
        r_element.Check(r_process_info);
    }

    for (const auto& r_condition : rModelPart.Conditions()) {
        r_condition.Check(r_process_info);
    }
}

} // namespace RansApplicationTestUtilities
} // namespace Kratos