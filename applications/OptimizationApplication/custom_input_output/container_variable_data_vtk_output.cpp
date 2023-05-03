//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main author:     Fabian Meister
//

#include "input_output/vtk_output.h"

// Include base h
#include "container_variable_data_vtk_output.h"

namespace Kratos
{

        explicit ContainerVariableDataVtkOutput(ModelPart &rModelPart, Parameters Parameters = Parameters(R"({})")) : VtkOutput(rModelPart, Parameters)
        {
        }

        void ContainerVariableDataVtkOutput::WriteContainerDataToFile(const ContainerVariableData<TContainerType> &rContainerVariableData, const std::string &rDataName, std::ofstream &rFileStream) const
        {
                rFileStream << "Test";
        }

        template void ContainerVariableDataVtkOutput::WriteContainerDataToFile(const ModelPart::NodesContainerType &rContainerVariableData, const std::string &rDataName, std::ofstream &rFileStream);
        template void ContainerVariableDataVtkOutput::WriteContainerDataToFile(const ModelPart::ConditionsContainerType &rContainerVariableData, const std::string &rDataName, std::ofstream &rFileStream);
        template void ContainerVariableDataVtkOutput::WriteContainerDataToFile(const ModelPart::ElementsContainerType &rContainerVariableData, const std::string &rDataName, std::ofstream &rFileStream);
}