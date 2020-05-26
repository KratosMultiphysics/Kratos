// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//

// System includes

// External includes

// Project includes
#include "custom_python/add_custom_utilities_to_python.h"

//Utilities
#include "custom_utilities/intersection_utilities.h"

namespace Kratos {
namespace Python {

void  AddCustomUtilitiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    // IntersectionUtilities
    m.def("FindIntersection1DGeometries2D", &IntersectionUtilities::FindIntersection1DGeometries2D);
    m.def("CreateQuadraturePointsCoupling1DGeometries2D", &IntersectionUtilities::CreateQuadraturePointsCoupling1DGeometries2D);

}

}  // namespace Python.
} // Namespace Kratos

