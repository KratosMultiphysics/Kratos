//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main author:     Riccardo Rossi
//

#ifndef KRATOS_ADD_DISTRIBUTED_SPARSE_MATRICES_TO_PYTHON_H_INCLUDED
#define KRATOS_ADD_DISTRIBUTED_SPARSE_MATRICES_TO_PYTHON_H_INCLUDED

// System includes

// External includes
#include "includes/define_python.h"

// Project includes

namespace Kratos {
namespace Python {

void AddDistributedSparseMatricesToPython(pybind11::module& m);

} // namespace Python
} // namespace Kratos

#endif // KRATOS_ADD_DISTRIBUTED_SPARSE_MATRICES_TO_PYTHON_H_INCLUDED

