//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Riccardo Rossi
//



// System includes

// External includes

// Project includes
#include "includes/define_python.h"
#include "python/add_search_strategies_to_python.h"
#include "spatial_containers/spatial_search.h"

namespace Kratos
{

namespace Python
{

void  AddSearchStrategiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    py::class_<SpatialSearch, SpatialSearch::Pointer>(m, "SpatialSearch")
        .def(py::init< >())
        ;
}

}  // namespace Python.

} // Namespace Kratos

