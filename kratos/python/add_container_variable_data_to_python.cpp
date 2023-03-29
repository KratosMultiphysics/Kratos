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
#include "containers/container_variable_data/container_data_io.h"
#include "containers/container_variable_data/container_variable_data.h"
#include "containers/container_variable_data/specialized_container_variable_data.h"
#include "add_container_variable_data_to_python_utils.h"

// Include base h
#include "add_container_variable_data_to_python.h"

namespace Kratos::Python
{

void  AddContainerVariableDataToPython(pybind11::module& m)
{
    auto sub_module = m.def_submodule("ContainerVariableData");

    AddContainerVariableDataToPython<ModelPart::NodesContainerType>(sub_module, "NodalVariableData");
    AddContainerVariableDataToPython<ModelPart::ConditionsContainerType>(sub_module, "ConditionVariableData");
    AddContainerVariableDataToPython<ModelPart::ElementsContainerType>(sub_module, "ElementVariableData");

    AddSpecializedContainerVariableDataToPython<ModelPart::NodesContainerType, ContainerDataIOTags::Historical>(sub_module, "HistoricalVariableData");
    AddSpecializedContainerVariableDataToPython<ModelPart::NodesContainerType, ContainerDataIOTags::NonHistorical>(sub_module, "NodalNonHistoricalVariableData");
    AddSpecializedContainerVariableDataToPython<ModelPart::ConditionsContainerType, ContainerDataIOTags::NonHistorical>(sub_module, "ConditionNonHistoricalVariableData");
    AddSpecializedContainerVariableDataToPython<ModelPart::ElementsContainerType, ContainerDataIOTags::NonHistorical>(sub_module, "ElementNonHistoricalVariableData");
}

} // namespace Kratos::Python
