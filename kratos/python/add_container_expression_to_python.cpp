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
#include "includes/define_python.h"
#include "expression/container_data_io.h"
#include "expression/container_expression.h"
#include "expression/specialized_container_expression.h"
#include "add_container_expression_to_python_utils.h"
#include "add_expression_io_to_python.h"

// Include base h
#include "add_container_expression_to_python.h"

namespace Kratos::Python
{


void  AddContainerExpressionToPython(pybind11::module& m)
{
    auto container_exp_sub_module = m.def_submodule("Expression");

    AddContainerExpressionToPython<ModelPart::NodesContainerType>(container_exp_sub_module, "NodalExpression");
    AddContainerExpressionToPython<ModelPart::ConditionsContainerType>(container_exp_sub_module, "ConditionExpression");
    AddContainerExpressionToPython<ModelPart::ElementsContainerType>(container_exp_sub_module, "ElementExpression");

    AddSpecializedContainerExpressionToPython<ModelPart::NodesContainerType, ContainerDataIOTags::Historical>(container_exp_sub_module, "HistoricalExpression");
    AddSpecializedContainerExpressionToPython<ModelPart::NodesContainerType, ContainerDataIOTags::NonHistorical>(container_exp_sub_module, "NodalNonHistoricalExpression");
    AddSpecializedContainerExpressionToPython<ModelPart::ConditionsContainerType, ContainerDataIOTags::NonHistorical>(container_exp_sub_module, "ConditionNonHistoricalExpression");
    AddSpecializedContainerExpressionToPython<ModelPart::ElementsContainerType, ContainerDataIOTags::NonHistorical>(container_exp_sub_module, "ElementNonHistoricalExpression");

    AddExpressionIOToPython(container_exp_sub_module);
}

} // namespace Kratos::Python
