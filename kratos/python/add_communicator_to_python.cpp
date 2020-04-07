//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main author:     Philipp Bucher
//

// System includes

// External includes

// Project includes
#include "includes/define_python.h"
#include "add_communicator_to_python.h"
#include "includes/communicator.h"

namespace Kratos {
namespace Python {

void AddCommunicatorToPython(pybind11::module &m)
{
    namespace py = pybind11;
}

} // namespace Python.
} // Namespace Kratos
