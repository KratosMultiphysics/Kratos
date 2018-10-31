//  KRATOS  _____     _ _ _
//         |_   _| __(_) (_)_ __   ___  ___
//           | || '__| | | | '_ \ / _ \/ __|
//           | || |  | | | | | | | (_) \__
//           |_||_|  |_|_|_|_| |_|\___/|___/ APPLICATION
//
//  License:             BSD License
//                                       Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

// System includes

#if defined(KRATOS_PYTHON)
// External includes
#include <pybind11/pybind11.h>

// Project includes
#include "includes/define_python.h"

//Trilinos includes
#include "mpi.h"

// Project includes
#include "custom_utilities/zoltan_partition_utility.h"
#include "custom_python/add_zoltan_processes_to_python.h"

namespace Kratos
{
namespace Python
{
namespace py = pybind11;

void AddZoltanProcessesToPython(pybind11::module& m)
{
    py::class_<ZoltanPartitionUtility >(m,"ZoltanPartitionUtility")
    .def(py::init< >() )
    .def("CalculatePartition", &ZoltanPartitionUtility::CalculatePartition )
    ;
}

} // namespace Python.

} // namespace Kratos.

#endif // KRATOS_PYTHON defined
