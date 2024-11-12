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

namespace Kratos::Python
{

void AddCustomUtilitiesToPython(pybind11::module& m)
{
    pybind11::class_<NodeUtilities>(m, "NodeUtilities")
        .def("AssignUpdatedVectorVariableToNonFixedComponentsOfNodes",
             &NodeUtilities::AssignUpdatedVectorVariableToNonFixedComponentsOfNodes);
}

} // Namespace Kratos::Python.
