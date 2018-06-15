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

void PythonObjectCppWrapperUtility::RunAnalysisStage(
        const std::string& rProjectParametersFile, 
        Parameters StageParameters
        )
{
    Parameters parameters = ReadParameters(rProjectParametersFile);

    Parameters default_stage_parameters = Parameters(R"(
    {
        "list_applications_to_import"          : ["StructuralMechanicsApplication"],
        "analysis_stage_name"                  : "structural_mechanics_analysis.StructuralMechanicsAnalysis"
    })" );
    
    StageParameters.ValidateAndAssignDefaults(default_stage_parameters);

    Model model = Model();

    pybind11::object kratos = pybind11::module::import("KratosMultiphysics");
    if (StageParameters["list_applications_to_import"].IsArray()) {
        auto list_applications_to_import = StageParameters["list_applications_to_import"];
        for (IndexType i_app = 0; i_app < list_applications_to_import.size(); ++i_app) {
            const std::string& app_name = list_applications_to_import[i_app].GetString();
            const char * char_name = ("KratosMultiphysics."+app_name).c_str();
            pybind11::object application = pybind11::module::import(char_name);
        }
    }
    std::string first, second;
    TrimComponentName(StageParameters["analysis_stage_name"].GetString(), first, second);
    const char * char_first = first.c_str();
    const char * char_second = second.c_str();
    pybind11::object analysis = pybind11::module::import(char_first).attr(char_second);
    pybind11::object my_analysis = analysis(model, parameters);
    auto run = my_analysis.attr("Run");
    run();
}

/***********************************************************************************/
/***********************************************************************************/

void PythonObjectCppWrapperUtility::RunStructuralAnalysisStage(const std::string& rProjectParametersFile)
{
    Parameters parameters = ReadParameters(rProjectParametersFile);

    Model model = Model();

    pybind11::object kratos = pybind11::module::import("KratosMultiphysics");
    pybind11::object structural = pybind11::module::import("KratosMultiphysics.StructuralMechanicsApplication");
    pybind11::object structural_analysis = pybind11::module::import("structural_mechanics_analysis").attr("StructuralMechanicsAnalysis");
    pybind11::object my_structural_analysis = structural_analysis(model, parameters);
    auto run = my_structural_analysis.attr("Run");
    run();
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

/***********************************************************************************/
/***********************************************************************************/

void PythonObjectCppWrapperUtility::TrimComponentName(
    const std::string& rCompleteName,
    std::string& rFirstPart,
    std::string& rSecondPart
    )
{
    std::string::size_type found = rCompleteName.find(".");
    rFirstPart = rCompleteName.substr(0, found - 1);
    rSecondPart = rCompleteName.substr(found + 1);
}
}  // namespace Kratos.
