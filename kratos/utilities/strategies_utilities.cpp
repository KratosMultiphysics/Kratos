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

// System includes

// External includes

// Project includes
#include "spaces/ublas_space.h"
#include "utilities/strategies_utilities.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"

namespace Kratos
{

/// The Ublas sparse space type
typedef UblasSpace<double, CompressedMatrix, Vector> UblasSparseSpaceType;

/***********************************************************************************/
/***********************************************************************************/

template<>
double StrategiesUtilities::CheckAndCorrectZeroDiagonalValues<UblasSparseSpaceType>(
    ModelPart& rModelPart,
    UblasSparseSpaceType::MatrixType& rA,
    UblasSparseSpaceType::VectorType& rb,
    const SCALING_DIAGONAL ScalingDiagonal
    )
{
    const std::size_t system_size = rA.size1();

    double* Avalues = rA.value_data().begin();
    std::size_t* Arow_indices = rA.index1_data().begin();

    // The diagonal considered
    const double scale_factor = GetScaleNorm<UblasSparseSpaceType>(rModelPart, rA, ScalingDiagonal);

    // Detect if there is a line of all zeros and set the diagonal to a 1 if this happens
    IndexPartition<std::size_t>(system_size).for_each([&](std::size_t Index){
        bool empty = true;

        const std::size_t col_begin = Arow_indices[Index];
        const std::size_t col_end = Arow_indices[Index + 1];

        for (std::size_t j = col_begin; j < col_end; ++j) {
            if(Avalues[j] != 0.0) {
                empty = false;
                break;
            }
        }

        if(empty) {
            rA(Index, Index) = scale_factor;
            rb[Index] = 0.0;
        }
    });

    return scale_factor;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
double StrategiesUtilities::GetDiagonalNorm<UblasSparseSpaceType>(UblasSparseSpaceType::MatrixType& rA)
{
    double* Avalues = rA.value_data().begin();
    std::size_t* Arow_indices = rA.index1_data().begin();
    std::size_t* Acol_indices = rA.index2_data().begin();

    const double diagonal_norm = IndexPartition<std::size_t>(UblasSparseSpaceType::Size1(rA)).for_each<SumReduction<double>>([&](std::size_t Index){
        const std::size_t col_begin = Arow_indices[Index];
        const std::size_t col_end = Arow_indices[Index+1];
        for (std::size_t j = col_begin; j < col_end; ++j) {
            if (Acol_indices[j] == Index ) {
                return std::pow(Avalues[j], 2);
            }
        }
        return 0.0;
    });

    return std::sqrt(diagonal_norm);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
double StrategiesUtilities::GetMaxDiagonal<UblasSparseSpaceType>(UblasSparseSpaceType::MatrixType& rA)
{
    double* Avalues = rA.value_data().begin();
    std::size_t* Arow_indices = rA.index1_data().begin();
    std::size_t* Acol_indices = rA.index2_data().begin();

    return IndexPartition<std::size_t>(UblasSparseSpaceType::Size1(rA)).for_each<MaxReduction<double>>([&](std::size_t Index){
        const std::size_t col_begin = Arow_indices[Index];
        const std::size_t col_end = Arow_indices[Index+1];
        for (std::size_t j = col_begin; j < col_end; ++j) {
            if (Acol_indices[j] == Index ) {
                return std::abs(Avalues[j]);
            }
        }
        return 0.0;
    });
}

/***********************************************************************************/
/***********************************************************************************/

template<>
double StrategiesUtilities::GetMinDiagonal<UblasSparseSpaceType>(UblasSparseSpaceType::MatrixType& rA)
{
    double* Avalues = rA.value_data().begin();
    std::size_t* Arow_indices = rA.index1_data().begin();
    std::size_t* Acol_indices = rA.index2_data().begin();

    return IndexPartition<std::size_t>(UblasSparseSpaceType::Size1(rA)).for_each<MinReduction<double>>([&](std::size_t Index){
        const std::size_t col_begin = Arow_indices[Index];
        const std::size_t col_end = Arow_indices[Index+1];
        for (std::size_t j = col_begin; j < col_end; ++j) {
            if (Acol_indices[j] == Index ) {
                return std::abs(Avalues[j]);
            }
        }
        return std::numeric_limits<double>::max();
    });
}

}  // namespace Kratos
