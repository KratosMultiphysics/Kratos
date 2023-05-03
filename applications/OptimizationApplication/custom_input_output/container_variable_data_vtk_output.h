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

// System includes

// External includes

// Project includes
#include "input_output/vtk_output.h"

#include "custom_utilities/container_variable_data/container_variable_data.h"

namespace Kratos
{
        /** \brief VtkOutput
         * A simple class that has functionality to write vtk outputs for container variable data - extention of the kratos core vtk output
         * @see : https://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf
         */

        /**
         * @brief Print the given rContainerVariableData as a VTK file
         * @param rContainerVariableData ContainerVariableData which is outputted
         * @param rDataName name of the data in the vtk file
         */
        class ContainerVariableDataVtkOutput : public VtkOutput
        {
        public:
                explicit ContainerVariableDataVtkOutput(ModelPart &rModelPart, Parameters Parameters = Parameters(R"({})"));

                template <class TContainerType>
                void WriteContainerDataToFile(const ContainerVariableData<TContainerType> &rContainerVariableData, const std::string &rDataName, std::ofstream &rFileStream);
        };
}