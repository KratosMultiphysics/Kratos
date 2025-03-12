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

// Project includes
#include "solving_strategies/builder_and_solvers/p_multigrid/diagonal_scaling.hpp" // DiagonalScaling, ParseDiagonalScaling, GetDiagonalScaleFactor
#include "spaces/ublas_space.h" // TUblasSparseSpace
#include "utilities/reduction_utilities.h" // MaxReduction, AbsMaxReduction
#include "utilities/profiler.h" // KRATOS_PROFILE_SCOPE

// STL includes
#include <optional> // std::optional


namespace Kratos {


DiagonalScaling ParseDiagonalScaling(Parameters Settings)
{
    KRATOS_TRY
    const auto diagonal_scaling_strategy = Settings["diagonal_scaling"].Get<std::string>();
    if (diagonal_scaling_strategy == "none") {
        return DiagonalScaling::None;
    } else if (diagonal_scaling_strategy == "abs_max") {
        return DiagonalScaling::AbsMax;
    } else if (diagonal_scaling_strategy == "norm") {
        return DiagonalScaling::Norm;
    } else {
        KRATOS_ERROR << "unsupported setting for \"diagonal_scaling\": "
                     << diagonal_scaling_strategy << ". Options are:\n"
                     << "- \"none\"\n"
                     << "- \"abs_max\"\n"
                     << "- \"norm\"\n";
    }
    KRATOS_CATCH("")
}


namespace detail {


template <class TSparse>
std::optional<typename TSparse::DataType> FindDiagonal(const typename TSparse::MatrixType& rMatrix,
                                                       const std::size_t iRow) noexcept
{
    const auto i_entry_begin = rMatrix.index1_data()[iRow];
    const auto i_entry_end   = rMatrix.index1_data()[iRow + 1];

    const auto it_column_begin = rMatrix.index2_data().begin() + i_entry_begin;
    const auto it_column_end = rMatrix.index2_data().begin() + i_entry_end;

    // Look for the diagonal entry in the current row.
    const auto it_column = std::lower_bound(it_column_begin, it_column_end, iRow);

    return (it_column == it_column_end or *it_column != iRow)
         ? std::optional<typename TSparse::DataType>()
         : rMatrix.value_data()[rMatrix.index1_data()[iRow] + std::distance(it_column_begin, it_column)];
}


} // namespace detail


template <class TSparse>
typename TSparse::DataType
GetDiagonalScaleFactor(const typename TSparse::MatrixType& rMatrix,
                       const DiagonalScaling ScalingStrategy)
{
    KRATOS_PROFILE_SCOPE(KRATOS_CODE_LOCATION);
    KRATOS_TRY
    switch (ScalingStrategy) {
        case DiagonalScaling::None:
            return 1;

        case DiagonalScaling::AbsMax:
            return IndexPartition(rMatrix.size1()).template for_each<AbsMaxReduction<typename TSparse::DataType>>(
                [&rMatrix](std::size_t i_row) -> typename TSparse::DataType {
                    const auto maybe_diagonal_entry = detail::FindDiagonal<TSparse>(rMatrix, i_row);
                    KRATOS_ERROR_IF_NOT(maybe_diagonal_entry.has_value()) << "row " << i_row << " has no diagonal entry";
                    return maybe_diagonal_entry.value();
            });

        case DiagonalScaling::Norm: {
            typename TSparse::DataType output = IndexPartition(rMatrix.size1()).template for_each<SumReduction<typename TSparse::DataType>>(
                [&rMatrix](std::size_t i_row) -> typename TSparse::DataType {
                    const auto maybe_diagonal_entry = detail::FindDiagonal<TSparse>(rMatrix, i_row);
                    KRATOS_ERROR_IF_NOT(maybe_diagonal_entry.has_value()) << "row " << i_row << " has no diagonal entry";
                    const auto diagonal_entry = maybe_diagonal_entry.value();
                    return diagonal_entry * diagonal_entry;
            });
            return std::sqrt(output) / rMatrix.size1();
        }

        default: {
            KRATOS_ERROR << "unsupported diagonal scaling (" << (int)ScalingStrategy << ')';
        }
    } // switch ScalingStrategy
    KRATOS_CATCH("")
}


template <class TSparse>
void NormalizeRows(typename TSparse::MatrixType& rLhs,
                   typename TSparse::VectorType& rRhs)
{
    KRATOS_TRY
    IndexPartition<std::size_t>(rLhs.size1()).for_each([&rLhs, &rRhs](const std::size_t i_row){
        const auto maybe_diagonal_entry = detail::FindDiagonal<TSparse>(rLhs, i_row);
        KRATOS_ERROR_IF_NOT(maybe_diagonal_entry.has_value()) << "row " << i_row << " has no diagonal entry";
        const auto diagonal_entry = maybe_diagonal_entry.value();

        if (diagonal_entry) {
            const auto scale = static_cast<typename TSparse::DataType>(1) / diagonal_entry;

            const auto i_entry_begin = rLhs.index1_data()[i_row];
            const auto i_entry_end = rLhs.index1_data()[i_row + 1];
            for (auto i_entry=i_entry_begin; i_entry<i_entry_end; ++i_entry) {
                rLhs.value_data()[i_entry] *= scale;
            } // for i_entry in range(i_entry_begin, i_entry_end)

            rRhs[i_row] *= scale;
        } /*if diagonal_entry*/ else {
            const auto i_entry_begin = rLhs.index1_data()[i_row];
            const auto i_entry_end = rLhs.index1_data()[i_row + 1];
            for (auto i_entry=i_entry_begin; i_entry<i_entry_end; ++i_entry) {
                const auto i_column = rLhs.index2_data()[i_entry];
                const auto value = rLhs.value_data()[i_entry];
                KRATOS_ERROR_IF(i_column != i_row and value) << "the diagonal of row " << i_row << " vanishes, but has off-diagonal components";
                //if (i_column == i_row) rLhs.value_data()[i_entry] = static_cast<typename TSparse::DataType>(1);
            } // for i_entry in range(i_entry_begin, i_entry_end)
            //KRATOS_ERROR_IF(rRhs[i_row]) << "row " << i_row << " is empty, but the corresponding RHS component is " << rRhs[i_row];
        } // if not diagonal_entry
    });
    KRATOS_CATCH("")
}


template <class TSparse>
void NormalizeSystem(typename TSparse::MatrixType& rLhs,
                     typename TSparse::VectorType& rRhs,
                     typename TSparse::DataType Coefficient)
{
    block_for_each(&*rLhs.value_data().begin(),
                   (&*rLhs.value_data().begin()) + rLhs.value_data().size(),
                  [Coefficient](typename TSparse::DataType& r_entry){r_entry *= Coefficient;});
    block_for_each(&*rRhs.begin(),
                   (&*rRhs.begin()) + rRhs.size(),
                  [Coefficient](typename TSparse::DataType& r_entry){r_entry *= Coefficient;});
}


// Explicit instantiations for UBLAS.
#define KRATOS_DEFINE_DIAGONAL_SCALE_FUNCTIONS(TSparseSpace)                                                        \
    template TSparseSpace::DataType GetDiagonalScaleFactor<TSparseSpace>(const typename TSparseSpace::MatrixType&,  \
                                                                         const DiagonalScaling);                    \
    template void NormalizeRows<TSparseSpace>(typename TSparseSpace::MatrixType&,                                   \
                                              typename TSparseSpace::VectorType&);                                  \
    template void NormalizeSystem<TSparseSpace>(typename TSparseSpace::MatrixType&,                                 \
                                                typename TSparseSpace::VectorType&,                                 \
                                                typename TSparseSpace::DataType)


KRATOS_DEFINE_DIAGONAL_SCALE_FUNCTIONS(TUblasSparseSpace<double>);

KRATOS_DEFINE_DIAGONAL_SCALE_FUNCTIONS(TUblasSparseSpace<float>);

#undef KRATOS_DEFINE_DIAGONAL_SCALE_FUNCTIONS


} // namespace Kratos
