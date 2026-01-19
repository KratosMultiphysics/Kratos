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
#include "containers/csr_matrix.h"
#include "containers/sparse_contiguous_row_graph.h"
#include "containers/system_vector.h"

namespace Kratos::Future
{

///@addtogroup KratosCore
///@{

///@name Kratos Classes
///@{

struct SerialLinearAlgebraTraits
{
    using DataType = double;

    using IndexType = std::size_t;

    using MatrixType = CsrMatrix<DataType, IndexType>;

    using VectorType = SystemVector<DataType, IndexType>;

    using DenseMatrixType = DenseMatrix<DataType>; //TODO: think about this one (could it be std::vector<VectorType> ?)

    using SparseGraphType = SparseContiguousRowGraph<IndexType>;
};

///@}
///@} addtogroup block

} // namespace Kratos::Future
