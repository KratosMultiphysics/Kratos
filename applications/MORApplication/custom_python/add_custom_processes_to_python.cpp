//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt

// System includes
#include "custom_python/add_custom_processes_to_python.h"

// External includes

// Project includes
// #include "includes/define.h"
// #include "includes/define_python.h"
#include "spaces/ublas_space.h"

#include "processes/process.h"
#include "custom_processes/monolithic_mapping_process.hpp"

namespace Kratos
{

namespace Python
{

void  AddCustomProcessesToPython(pybind11::module& m)
{

    namespace py = pybind11;
    typedef TUblasSparseSpace<double> SparseSpaceType;
    typedef SparseSpaceType::MatrixType SparseMatrixType;

    py::class_<MonolithicMappingProcess<SparseMatrixType>, typename MonolithicMappingProcess<SparseMatrixType>::Pointer, Process>(m, "MonolithicMappingProcess")
        .def(py::init<ModelPart&, ModelPart&, SparseMatrixType&>())
        ;

}

}  // namespace Python.

} // Namespace Kratos
