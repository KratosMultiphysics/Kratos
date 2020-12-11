//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Manuel Messmer
//

// System includes
#include <pybind11/stl.h>

// External includes

// Project includes
#include "custom_utilities/mpi/mpi_utilities.h"

namespace Kratos{
namespace Python{

    void  AddCustomMPIUtilitiesToPython(pybind11::module& m)
    {
        namespace py = pybind11;

        py::class_<MPM_MPI_Utilities, MPM_MPI_Utilities::Pointer>(m, "MPM_MPI_Utilities")
            .def_static("TransferElements", &MPM_MPI_Utilities::TransferElements)
            .def_static("TransferConditions", &MPM_MPI_Utilities::TransferConditions);
    }

}  // namespace Python.
} // Namespace Kratos