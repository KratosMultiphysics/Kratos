// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Ignasi de Pouplana,
//                   Vahid Galavi
//

// External includes

// Project includes
#include "custom_python/add_custom_utilities_to_python.h"

#include "custom_utilities/node_utilities.h"
#include "custom_utilities/process_utilities.h"
#include "custom_workflows/custom_workflow_factory.h"
#include "custom_workflows/dgeoflow.h"
#include "custom_workflows/dgeosettlement.hpp"

namespace Kratos::Python
{

void AddCustomUtilitiesToPython(const pybind11::module& rModule)
{
    pybind11::class_<NodeUtilities>(rModule, "NodeUtilities")
        .def("AssignUpdatedVectorVariableToNodes", &NodeUtilities::AssignUpdatedVectorVariableToNodes);

    pybind11::class_<CustomWorkflowFactory>(rModule, "CustomWorkflowFactory")
        .def_static("CreateKratosGeoSettlement", &CustomWorkflowFactory::CreateKratosGeoSettlement,
                    pybind11::return_value_policy::take_ownership)
        .def_static("CreateKratosGeoFlow", &CustomWorkflowFactory::CreateKratosGeoFlow,
                    pybind11::return_value_policy::take_ownership);

    pybind11::class_<KratosGeoSettlement>(rModule, "KratosGeoSettlement").def("RunStage", &KratosGeoSettlement::RunStage);

    pybind11::class_<KratosExecute>(rModule, "KratosExecute").def("ExecuteFlowAnalysis", &KratosExecute::ExecuteFlowAnalysis);
    pybind11::class_<KratosExecute::CriticalHeadInfo>(rModule, "KratosExecuteCriticalHeadInfo")
        .def(pybind11::init<double, double, double>());
    pybind11::class_<KratosExecute::CallBackFunctions>(rModule, "KratosExecuteCallBackFunctions")
        .def(pybind11::init<std::function<void(const char*)>, std::function<void(double)>,
                            std::function<void(const char*)>, std::function<bool()>>());

    pybind11::class_<ProcessUtilities>(rModule, "ProcessUtilities")
        .def_static("AddProcessesSubModelPartListToSolverSettings",
                    &ProcessUtilities::AddProcessesSubModelPartListToSolverSettings,
                    pybind11::arg("project_parameters"), pybind11::arg("solver_settings"));
}

} // Namespace Kratos::Python.
