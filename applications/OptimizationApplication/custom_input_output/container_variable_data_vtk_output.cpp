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

        ContainerVariableDataVtkOutput::ContainerVariableDataVtkOutput(ModelPart &rModelPart, Parameters parameters = Parameters(R"({})")) : VtkOutput(rModelPart, parameters){

                                                                                                                                             };

        void ContainerVariableDataVtkOutput::TestFunction()
        {
                std::fstream my_file;
                my_file.open("test.txt");
                my_file << "Test";
                my_file.close();
        };

        template <class TContainerType>
        void ContainerVariableDataVtkOutput::WriteContainerDataToFile(const ContainerVariableData<TContainerType> &rContainerVariableData, const std::string &rDataName){

        };

        template <class number>
        void ContainerVariableDataVtkOutput::printNumberType(const number num){
                std::cout<<"some number: "<<num<<std::endl;
        };

        template <>
        void ContainerVariableDataVtkOutput::printNumberType(const int num){
                std::cout<<"an int: "<<num<<std::endl;
        };

        template void ContainerVariableDataVtkOutput::printNumberType(double);
}