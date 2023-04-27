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
#include <numeric>

// External includes
#include <pybind11/stl.h>
#include <pybind11/operators.h>
#include <pybind11/numpy.h>

// Project includes
#include "containers/container_expression/container_data_io.h"
#include "containers/container_expression/container_expression.h"
#include "containers/container_expression/specialized_container_expression.h"
#include "add_container_expression_to_python_utils.h"

// Include base h
#include "add_container_expression_to_python.h"

namespace Kratos::Python
{

void  AddContainerExpressionToPython(pybind11::module& m)
{
    auto sub_module = m.def_submodule("ContainerExpression");

    AddContainerExpressionToPython<ModelPart::NodesContainerType>(sub_module, "NodalExpression");
    AddContainerExpressionToPython<ModelPart::ConditionsContainerType>(sub_module, "ConditionExpression");
    AddContainerExpressionToPython<ModelPart::ElementsContainerType>(sub_module, "ElementExpression");

    AddSpecializedContainerExpressionToPython<ModelPart::NodesContainerType, ContainerDataIOTags::Historical>(sub_module, "HistoricalExpression");
    AddSpecializedContainerExpressionToPython<ModelPart::NodesContainerType, ContainerDataIOTags::NonHistorical>(sub_module, "NodalNonHistoricalExpression");
    AddSpecializedContainerExpressionToPython<ModelPart::ConditionsContainerType, ContainerDataIOTags::NonHistorical>(sub_module, "ConditionNonHistoricalExpression");
    AddSpecializedContainerExpressionToPython<ModelPart::ElementsContainerType, ContainerDataIOTags::NonHistorical>(sub_module, "ElementNonHistoricalExpression");
}

} // namespace Kratos::Python
