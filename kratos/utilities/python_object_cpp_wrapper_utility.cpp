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

// System includes

// External includes
#include <pybind11/operators.h>

// Project includes
#include "utilities/python_object_cpp_wrapper_utility.h"
#include "containers/model.h"

namespace Kratos
{
PythonObjectCppWrapperUtility::PythonObjectCppWrapperUtility(const std::string& rNameFile)
{
    const char * name = rNameFile.c_str();
    pybind11::object object = pybind11::module::import(name);
    mListPythonObjects.resize(1);
    mListPythonObjects[0] = object;
}

/***********************************************************************************/
/***********************************************************************************/

PythonObjectCppWrapperUtility::PythonObjectCppWrapperUtility(ListType& rObjectList)
{
    const std::size_t size_objects = len(rObjectList);
    mListPythonObjects.resize(size_objects);
    for (std::size_t i_object = 0; i_object < size_objects; ++i_object)
        mListPythonObjects[i_object] = pybind11::cast<ObjectType>(rObjectList[i_object]);
}

/***********************************************************************************/
/***********************************************************************************/

PythonObjectCppWrapperUtility::PythonObjectCppWrapperUtility(ObjectType& rObject)
{
    mListPythonObjects.resize(1);
    mListPythonObjects[0] = rObject;
}

/***********************************************************************************/
/***********************************************************************************/

void PythonObjectCppWrapperUtility::AddObject(ObjectType& rObject)
{
    mListPythonObjects.push_back(rObject);
}

/***********************************************************************************/
/***********************************************************************************/

void PythonObjectCppWrapperUtility::AddObjects(ListType& rObjectList)
{
    const std::size_t size_objects = len(rObjectList);
    for (std::size_t i_object = 0; i_object < size_objects; ++i_object)
        mListPythonObjects.push_back(pybind11::cast<ObjectType>(rObjectList[i_object]));
}

/***********************************************************************************/
/***********************************************************************************/

void PythonObjectCppWrapperUtility::Execute(const std::string& rNameMethod)
{
    for (auto& process : mListPythonObjects)
        process.attr(rNameMethod.c_str())();
}

/***********************************************************************************/
/***********************************************************************************/

void PythonObjectCppWrapperUtility::RunStructuralAnalysisStage(const std::string& rProjectParametersFile)
{
    Parameters parameters = ReadParameters(rProjectParametersFile);

    Model model = Model();

    pybind11::object kratos = pybind11::module::import("KratosMultiphysics");
    pybind11::object structural = pybind11::module::import("KratosMultiphysics.StructuralMechanicsApplication");
    pybind11::object structura_analysis = pybind11::module::import("structural_mechanics_analysis").attr("StructuralMechanicsAnalysis");
//     pybind11::object init = structura_analysis.attr("__init__");
//     init(structura_analysis, model, parameters);
//     structura_analysis.attr("Run");
}

/***********************************************************************************/
/***********************************************************************************/

Parameters PythonObjectCppWrapperUtility::ReadParameters(const std::string& rProjectParametersFile)
{
    std::ifstream infile(rProjectParametersFile);
    KRATOS_ERROR_IF_NOT(infile.good()) << "Materials file: " << rProjectParametersFile << " cannot be found" << std::endl;
    std::stringstream buffer;
    buffer << infile.rdbuf();
    Parameters parameters(buffer.str());

    return parameters;
}
}  // namespace Kratos.
