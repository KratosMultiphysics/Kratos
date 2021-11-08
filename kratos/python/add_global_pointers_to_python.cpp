//
//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//


// System includes

// External includes

// Project includes
#include "includes/define_python.h"
#include "processes/process.h"
#include "python/add_global_pointers_to_python.h"
#include "containers/global_pointers_vector.h"
#include "utilities/pointer_communicator.h"
#include "processes/process.h"
#include "processes/find_global_nodal_neighbours_process.h"

namespace Kratos {
namespace Python {


void AddGlobalPointersToPython(pybind11::module& m)
{
    namespace py = pybind11;

    // NOTE: this function is special in that it accepts a "pyObject" - this is the reason for which it is defined in this same file
    py::class_< GlobalPointer<Node<3>>  >(m,"GlobalNodePointer");
    py::class_< GlobalPointer<Element>  >(m,"GlobalElementPointer");
    py::class_< GlobalPointer<Condition>  >(m,"GlobalConditionPointer");

    py::class_< GlobalPointersVector<Node<3>> >(m,"GlobalNodePointersVector");

}

} // namespace Python.
} // Namespace Kratos
