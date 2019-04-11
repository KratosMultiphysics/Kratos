//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Guillermo Casas (gcasas@cimne.upc.edu)
//
//


// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "custom_python/add_custom_processes_to_python.h"
#include "includes/define_python.h"
#include "includes/model_part.h"
#include "processes/process.h"

namespace Kratos
{

namespace Python
{

void AddCustomProcessesToPython(pybind11::module& m)
{
    namespace py = pybind11;

}

} // namespace Python.

} // Namespace Kratos
