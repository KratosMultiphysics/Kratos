
//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher, Jordi Cotela
//
// See Master-Thesis P.Bucher
// "Development and Implementation of a Parallel
//  Framework for Non-Matching Grid Mapping"

// System includes

// External includes

// Project includes
#include "includes/define_python.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_utilities/mapper_utilities.h"


namespace Kratos {
namespace Python {

void  AddCustomUtilitiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    auto m_mapper_utilities = m.def_submodule("MapperUtilities");

    m_mapper_utilities.def("SaveCurrentConfiguration", &MapperUtilities::SaveCurrentConfiguration);
    m_mapper_utilities.def("RestoreCurrentConfiguration", &MapperUtilities::RestoreCurrentConfiguration);
}

}  // namespace Python.
} // Namespace Kratos
