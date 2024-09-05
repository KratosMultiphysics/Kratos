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
#include "includes/kratos_parameters.h"

#include "custom_utilities/condition_utilities.hpp"
#include "custom_utilities/element_utilities.hpp"
#include "custom_utilities/interface_element_utilities.h"

namespace Kratos
{
namespace Python
{

void AddCustomUtilitiesToPython(pybind11::module& m) { namespace py = pybind11; }

} // Namespace Python.
} // Namespace Kratos
