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
#include "custom_operations/primal_to_adjoint_operation.h"
#include "custom_operations/define_2d_wake_operation.h"
#include "custom_operations/define_3d_wake_operation.h"

namespace Kratos::Python
{

void AddCustomOperationsToPython(pybind11::module& m)
{
    namespace py = pybind11;

    py::class_<PotentialToCompressibleNavierStokesOperation, PotentialToCompressibleNavierStokesOperation::Pointer, Operation > (m,"PotentialToCompressibleNavierStokesOperation")
    .def(py::init<Model&, Parameters>())
    ;

    py::class_<PrimalToAdjointOperation, PrimalToAdjointOperation::Pointer, Operation > (m,"PrimalToAdjointOperation")
    .def(py::init<Model&, Parameters>())
    ;

    py::class_<Define2DWakeOperation, Define2DWakeOperation::Pointer, Operation > (m,"Define2DWakeOperation")
    .def(py::init< Model&, Parameters >())
    ;

    py::class_<Define3DWakeOperation, Define3DWakeOperation::Pointer, Operation > (m,"Define3DWakeOperation")
    .def(py::init< Model&, Parameters >())
    ;
}

} // Namespace Kratos::Python
