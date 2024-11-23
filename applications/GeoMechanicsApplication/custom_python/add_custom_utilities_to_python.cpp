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
#include "custom_workflows/custom_workflow_factory.h"
#include "custom_workflows/dgeosettlement.h"

namespace Kratos::Python
{

void AddCustomUtilitiesToPython(const pybind11::module& rModule)
{
    pybind11::class_<NodeUtilities>(rModule, "NodeUtilities")
        .def("AssignUpdatedVectorVariableToNonFixedComponentsOfNodes",
             &NodeUtilities::AssignUpdatedVectorVariableToNonFixedComponentsOfNodes);

    pybind11::class_<CustomWorkflowFactory>(rModule, "CustomWorkflowFactory")
        .def("CreateKratosGeoSettlement", &CustomWorkflowFactory::CreateKratosGeoSettlement,
             pybind11::return_value_policy::take_ownership);

    pybind11::class_<KratosGeoSettlement>(rModule, "KratosGeoSettlement").def("RunStage", &KratosGeoSettlement::RunStage);
}

} // Namespace Kratos::Python.
