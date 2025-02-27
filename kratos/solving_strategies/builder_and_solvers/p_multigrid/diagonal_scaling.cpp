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
#include "includes/variables.h" // BUILD_SCALE_FACTOR
#include "spaces/ublas_space.h" // TUblasSparseSpace
#include "utilities/reduction_utilities.h" // MaxReduction, AbsMaxReduction
#include "utilities/profiler.h" // KRATOS_PROFILE_SCOPE


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
                    const auto i_entry_begin = rMatrix.index1_data()[i_row];
                    const auto i_entry_end   = rMatrix.index1_data()[i_row + 1];

                    const auto it_column_begin = rMatrix.index2_data().begin() + i_entry_begin;
                    const auto it_column_end = rMatrix.index2_data().begin() + i_entry_end;

                    // Look for the diagonal entry in the current row.
                    const auto it_column = std::lower_bound(it_column_begin, it_column_end, i_row);
                    KRATOS_ERROR_IF(it_column == it_column_end or *it_column != i_row)
                        << "row " << i_row << " has no diagonal entry";

                    const auto diagonal_entry = rMatrix.value_data()[rMatrix.index1_data()[i_row] + std::distance(it_column_begin, it_column)];
                    return diagonal_entry;
            });

        case DiagonalScaling::Norm: {
            typename TSparse::DataType output = IndexPartition(rMatrix.size1()).template for_each<SumReduction<typename TSparse::DataType>>(
                [&rMatrix](std::size_t i_row) -> typename TSparse::DataType {
                    const auto i_entry_begin = rMatrix.index1_data()[i_row];
                    const auto i_entry_end   = rMatrix.index1_data()[i_row + 1];

                    const auto it_column_begin = rMatrix.index2_data().begin() + i_entry_begin;
                    const auto it_column_end = rMatrix.index2_data().begin() + i_entry_end;

                    // Look for the diagonal entry in the current row.
                    const auto it_column = std::lower_bound(it_column_begin, it_column_end, i_row);
                    KRATOS_ERROR_IF(it_column == it_column_end or *it_column != i_row)
                        << "row " << i_row << " has no diagonal entry";

                    const auto diagonal_entry = rMatrix.value_data()[rMatrix.index1_data()[i_row] + std::distance(it_column_begin, it_column)];
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


// Explicit instantiations for UBLAS.
#define KRATOS_DEFINE_DIAGONAL_SCALE_FACTOR_GETTER(TSparseSpace)                                                    \
    template TSparseSpace::DataType GetDiagonalScaleFactor<TSparseSpace>(const typename TSparseSpace::MatrixType&,  \
                                                                         const DiagonalScaling)

KRATOS_DEFINE_DIAGONAL_SCALE_FACTOR_GETTER(TUblasSparseSpace<double>);

KRATOS_DEFINE_DIAGONAL_SCALE_FACTOR_GETTER(TUblasSparseSpace<float>);

#undef KRATOS_DEFINE_DIAGONAL_SCALE_FACTOR_GETTER


} // namespace Kratos
