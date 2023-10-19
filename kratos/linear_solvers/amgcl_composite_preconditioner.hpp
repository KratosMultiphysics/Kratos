//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Máté Kelemen
//

#pragma once

// External includes
#include "amgcl/adapter/crs_tuple.hpp"
#include "amgcl/adapter/ublas.hpp"
#include "amgcl/adapter/zero_copy.hpp"
#include "amgcl/backend/builtin.hpp"
#include "amgcl/value_type/static_matrix.hpp"
#include "amgcl/make_solver.hpp"
#include "amgcl/make_block_solver.hpp"
#include "amgcl/solver/runtime.hpp"

// Project includes
#include "utilities/parallel_utilities.h"

// STL includes
#include <array> // array
#include <vector> // vector
#include <algorithm> // count, count_if


namespace Kratos {


namespace detail {


/// Composite preconditioner traits
struct CPTraits
{
    using Mask = std::vector<char>;

    struct AMGCL
    {
        using ScalarBackend = amgcl::backend::builtin<double>;
    }; // struct amgcl
}; // struct CPTraits


/** @brief Generate four submatrices from a square matrix.
 *
 *  @details Let @f(A@f) be the square input matrix, and
 *           @f(A_{ll}@f), @f(A_{lq}@f), @f(A_{ql}@f), @f(A_{qq}@f) the
 *           output matrices in this exact order. Then this function resizes
 *           and populates the output matrices such that
 *           @f[
 *              A = \begin{bmatrix}
 *                  A_{ll} & A_{lq} \\
 *                  A_{ql} & A_{qq}
 *              \end{bmatrix}
 *           @f]
 *
 *          @f(A_{ll}@f) is defined by a boolean mask, represented by a @c char
 *          array that evaluate to true at row indices where @f(A_{ll}@f) has
 *          components.
 */
inline void MakeSubblocks(const CPTraits::AMGCL::ScalarBackend::matrix& rRootMatrix,
                          const CPTraits::Mask& rMask,
                          const std::array<CPTraits::AMGCL::ScalarBackend::matrix*,4>& rOutput)
{
    KRATOS_ERROR_IF_NOT(rMask.size() == rRootMatrix.nrows)
        << "Mask size mismatch: mask of size " << rMask.size()
        << " provided for matrix of size " << rRootMatrix.nrows << "x" << rRootMatrix.ncols;

    KRATOS_ERROR_IF_NOT(rRootMatrix.nrows == rRootMatrix.ncols)
        << "Expecting a square matrix, but got " << rRootMatrix.nrows << "x" << rRootMatrix.ncols;

    KRATOS_ERROR_IF_NOT(std::count(rOutput.begin(), rOutput.end(), nullptr) == 0)
        << "Output pointers must point to existing matrices";

    const std::size_t system_size = rMask.size();

    // Compute the index of each component local to the submatrix it belongs to.
    std::vector<std::size_t> block_local_indices(system_size);
    std::size_t masked_count = 0;
    std::size_t unmasked_count = 0;
    for (std::size_t i_mask=0; i_mask<system_size; ++i_mask) {
        const auto mask_value = rMask[i_mask];
        if (mask_value == 1) {
            block_local_indices[i_mask] = masked_count++;
        } else if (mask_value == 0) {
            block_local_indices[i_mask] = unmasked_count++;
        } else {
            KRATOS_ERROR << "Invalid mask value '" << static_cast<int>(mask_value)
                         << "' at index " << i_mask << ". The mask must only consist of "
                         << " 0s or 1s.";
        }
    }

    // Resize submatrices.
    KRATOS_TRY
    rOutput[0]->set_size(masked_count, masked_count, true);
    rOutput[1]->set_size(masked_count, unmasked_count, true);
    rOutput[2]->set_size(unmasked_count, masked_count, true);
    rOutput[3]->set_size(unmasked_count, unmasked_count, true);
    KRATOS_CATCH("")

    // Compute sparsity patterns.
    KRATOS_TRY
    IndexPartition<std::size_t>(system_size).for_each([&rOutput, &rMask, &rRootMatrix, &block_local_indices] (std::size_t i_root) {
        const std::size_t i_block_local = block_local_indices[i_root];
        const char row_mask = rMask[i_root];

        for (auto it_column=rRootMatrix.row_begin(i_root); it_column; ++it_column) {
            const char column_mask = rMask[it_column.col()];
            const char i_block = (~(column_mask | (row_mask << 1))) & 0b11; // <== flattened index of the block
            ++rOutput[i_block]->ptr[i_block_local + 1];
        }
    });

    for (CPTraits::AMGCL::ScalarBackend::matrix* p_output : rOutput) {
        p_output->set_nonzeros(p_output->scan_row_sizes());
    }
    KRATOS_CATCH("")

    // Copy data to submatrices
    KRATOS_TRY
    IndexPartition<std::size_t>(system_size).for_each([&rRootMatrix, &rMask, &rOutput, &block_local_indices](std::size_t i_root) -> void {
        std::array<std::size_t,4> head_indices {0ul, 0ul, 0ul, 0ul};
        const std::size_t i_block_row = block_local_indices[i_root];
        const char row_mask = rMask[i_root];

        const std::size_t i_head = row_mask << 1;
        head_indices[i_head] = rOutput[i_head]->ptr[i_block_row];
        head_indices[i_head + 1] = rOutput[i_head + 1]->ptr[i_block_row];

        for (auto it_column=rRootMatrix.row_begin(i_root); it_column; ++it_column) {
            const std::size_t i_column = static_cast<std::size_t>(it_column.col());
            const std::size_t i_block_column = block_local_indices[i_column];
            const double value = it_column.value();
            const char column_mask = rMask[i_column];

            const char i_block = (~(column_mask | (row_mask << 1))) & 0b11; // <== flattened index of the block
            rOutput[i_block]->col[head_indices[i_block]] = i_block_column;
            rOutput[i_block]->val[head_indices[i_block]] = value;
            ++head_indices[i_block];
        }
    });
    KRATOS_CATCH("")
}


} // namespace detail


} // namespace Kratos
