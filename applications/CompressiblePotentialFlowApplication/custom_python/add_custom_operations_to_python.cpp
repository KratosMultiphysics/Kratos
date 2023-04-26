//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Marco Antonio Zuñiga Perez
//

// System includes

// External includes

// Project includes
#include "includes/define_python.h"
#include "operations/operation.h"
#include "custom_python/add_custom_operations_to_python.h"
#include "custom_operations/potential_to_compressible_navier_stokes_operation.h"

namespace Kratos::Python
{

void AddCustomOperationsToPython(pybind11::module& m)
{
    namespace py = pybind11;

    py::class_<PotentialToCompressibleNavierStokesOperation, PotentialToCompressibleNavierStokesOperation::Pointer, Operation > (m,"PotentialToCompressibleNavierStokesOperation")
    .def(py::init<Model&, Parameters>())
    ;
}

} // Namespace Kratos::Python
