//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#pragma once

// System includes
#include <cstddef>
#include <memory>
#include <tuple>
#include <type_traits>

// External includes
#include <amgcl/adapter/crs_tuple.hpp>
#include <amgcl/adapter/zero_copy.hpp>
#include <amgcl/backend/builtin.hpp>
#include <span/span.hpp>

namespace Kratos {

/**
 * @brief amgcl CSR adapter tolerant of the index width of the sparse backend.
 * @details amgcl::adapter::zero_copy reinterprets the index arrays as
 * std::ptrdiff_t and therefore statically rejects narrower index types — such
 * as the (default) 32-bit indices of the Eigen sparse backend. For
 * ptrdiff_t-sized indices this is the plain zero-copy view; narrower indices
 * are widened while copying through amgcl's generic range-tuple interface,
 * a copy amgcl's hierarchy/factorization construction performs anyway.
 * @tparam TPtr Type of the row pointer entries.
 * @tparam TCol Type of the column index entries.
 * @tparam TVal Type of the value entries.
 * @param NRows Number of rows of the matrix.
 * @param pRowExtents Row pointer array of the CSR matrix, with NRows + 1 entries.
 * @param pColumnIndices Column index array of the nonzero entries.
 * @param pValues Value array of the nonzero entries.
 * @return Shared pointer to the amgcl CSR matrix: a zero-copy view of the arrays when the index type has the width
 * of std::ptrdiff_t, otherwise a copy with widened indices.
 */
template<typename TPtr, typename TCol, typename TVal>
std::shared_ptr<amgcl::backend::crs<TVal>> MakeAmgclZeroCopyAdapter(
    const std::size_t NRows,
    const TPtr* pRowExtents,
    const TCol* pColumnIndices,
    const TVal* pValues
    )
{
    if constexpr (sizeof(TPtr) == sizeof(std::ptrdiff_t) && sizeof(TCol) == sizeof(std::ptrdiff_t)) {
        return amgcl::adapter::zero_copy(NRows, pRowExtents, pColumnIndices, pValues);
    } else {
        const std::size_t nnz = NRows ? static_cast<std::size_t>(pRowExtents[NRows]) : 0;
        return std::make_shared<amgcl::backend::crs<TVal>>(std::make_tuple(
            NRows,
            Kratos::span<const TPtr>(pRowExtents, NRows + 1),
            Kratos::span<const TCol>(pColumnIndices, nnz),
            Kratos::span<const TVal>(pValues, nnz)));
    }
}

} // namespace Kratos
