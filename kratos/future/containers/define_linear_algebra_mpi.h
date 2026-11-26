//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

#pragma once

// System includes

// External includes

// Project includes
#include "containers/distributed_csr_matrix.h"
#include "containers/distributed_sparse_graph.h"
#include "containers/distributed_system_vector.h"

namespace Kratos::Future
{

///@addtogroup KratosCore
///@{

///@name Kratos Classes
///@{

struct DistributedLinearAlgebraTraits
{
    using DataType = double;

    using IndexType = std::size_t;

    using MatrixType = DistributedCsrMatrix<DataType, IndexType>;

    using VectorType = DistributedSystemVector<DataType, IndexType>;

    using SparseGraphType = DistributedSparseGraph<IndexType>;
};

///@}
///@} addtogroup block

} // namespace Kratos::Future
