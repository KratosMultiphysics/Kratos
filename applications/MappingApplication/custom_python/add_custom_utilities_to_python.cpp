//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Peter Wilson
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

