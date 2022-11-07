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

// External includes

// Project includes
#include "includes/model_part.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

// Scaling enum
enum class SCALING_DIAGONAL {NO_SCALING = 0, CONSIDER_NORM_DIAGONAL = 1, CONSIDER_MAX_DIAGONAL = 2, CONSIDER_PRESCRIBED_DIAGONAL = 3};

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{
    
/**
 * @class StrategiesUtilities
 * @ingroup KratosCore
 * @brief This namespace includes several utilities necessaries for the strategies
 * @author Vicente Mataix Ferrandiz
 */
class StrategiesUtilities
{
public:
    ///@name Operations
    ///@{

    /**
     * @brief This method checks and corrects the zero diagonal values
     * @details This method returns the scale norm considering for scaling the diagonal
     * @param rModelPart The problem model part
     * @param rA The LHS matrix
     * @param rb The RHS vector
     * @param ScalingDiagonal The type of caling diagonal considered
     * @return The scale norm
     */
    template<class TSparseSpace>
    static double CheckAndCorrectZeroDiagonalValues(
        ModelPart& rModelPart,
        typename TSparseSpace::MatrixType& rA,
        typename TSparseSpace::VectorType& rb,
        const SCALING_DIAGONAL ScalingDiagonal = SCALING_DIAGONAL::NO_SCALING
        )
    {
        const std::size_t system_size = rA.size1();

        double* Avalues = rA.value_data().begin();
        std::size_t* Arow_indices = rA.index1_data().begin();

        // Define  zero value tolerance
        const double zero_tolerance = std::numeric_limits<double>::epsilon();

        // The diagonal considered
        const double scale_factor = GetScaleNorm<TSparseSpace>(rModelPart.GetProcessInfo(), rA, ScalingDiagonal);

        // Detect if there is a line of all zeros and set the diagonal to a 1 if this happens
        IndexPartition<std::size_t>(system_size).for_each([&](std::size_t Index){
            bool empty = true;

            const std::size_t col_begin = Arow_indices[Index];
            const std::size_t col_end = Arow_indices[Index + 1];

            for (std::size_t j = col_begin; j < col_end; ++j) {
                if(std::abs(Avalues[j]) < zero_tolerance) {
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

    /**
     * @brief This method returns the scale norm considering for scaling the diagonal
     * @param rProcessInfo The problem process info
     * @param rA The LHS matrix
     * @param ScalingDiagonal The type of caling diagonal considered
     * @return The scale norm
     */
    template<class TSparseSpace>
    static double GetScaleNorm(
        const ProcessInfo& rProcessInfo,
        typename TSparseSpace::MatrixType& rA,
        const SCALING_DIAGONAL ScalingDiagonal = SCALING_DIAGONAL::NO_SCALING
        )
    {
        switch (ScalingDiagonal) {
            case SCALING_DIAGONAL::NO_SCALING:
                return 1.0;
            case SCALING_DIAGONAL::CONSIDER_PRESCRIBED_DIAGONAL: {
                KRATOS_ERROR_IF_NOT(rProcessInfo.Has(BUILD_SCALE_FACTOR)) << "Scale factor not defined at process info" << std::endl;
                return rProcessInfo.GetValue(BUILD_SCALE_FACTOR);
            }
            case SCALING_DIAGONAL::CONSIDER_NORM_DIAGONAL:
                return GetDiagonalNorm<TSparseSpace>(rA)/static_cast<double>(rA.size1());
            case SCALING_DIAGONAL::CONSIDER_MAX_DIAGONAL:
                return GetMaxDiagonal<TSparseSpace>(rA);
//                 return TSparseSpace::TwoNorm(rA)/static_cast<double>(rA.size1());
            default:
                return GetMaxDiagonal<TSparseSpace>(rA);
        }
    }

    /**
     * @brief This method returns the diagonal norm considering for scaling the diagonal
     * @param rA The LHS matrix
     * @return The diagonal norm
     */
    template<class TSparseSpace>
    static double GetDiagonalNorm(typename TSparseSpace::MatrixType& rA)
    {
        const double* Avalues = rA.value_data().begin();
        const std::size_t* Arow_indices = rA.index1_data().begin();
        const std::size_t* Acol_indices = rA.index2_data().begin();

        const double diagonal_norm = IndexPartition<std::size_t>(TSparseSpace::Size1(rA)).for_each<SumReduction<double>>([&](std::size_t Index){
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

    /**
     * @brief This method returns the diagonal max value
     * @param rA The LHS matrix
     * @return The diagonal  max value
     */
    template<class TSparseSpace>
    static double GetAveragevalueDiagonal(typename TSparseSpace::MatrixType& rA)
    {
        return 0.5 * (GetMaxDiagonal<TSparseSpace>(rA) + GetMinDiagonal<TSparseSpace>(rA));
    }

    /**
     * @brief This method returns the diagonal max value
     * @param rA The LHS matrix
     * @return The diagonal  max value
     */
    template<class TSparseSpace>
    static double GetMaxDiagonal(typename TSparseSpace::MatrixType& rA)
    {
        double* Avalues = rA.value_data().begin();
        std::size_t* Arow_indices = rA.index1_data().begin();
        std::size_t* Acol_indices = rA.index2_data().begin();

        return IndexPartition<std::size_t>(TSparseSpace::Size1(rA)).for_each<MaxReduction<double>>([&](std::size_t Index){
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

    /**
     * @brief This method returns the diagonal min value
     * @param rA The LHS matrix
     * @return The diagonal min value
     */
    template<class TSparseSpace>
    static double GetMinDiagonal(typename TSparseSpace::MatrixType& rA)
    {
        double* Avalues = rA.value_data().begin();
        std::size_t* Arow_indices = rA.index1_data().begin();
        std::size_t* Acol_indices = rA.index2_data().begin();

        return IndexPartition<std::size_t>(TSparseSpace::Size1(rA)).for_each<MinReduction<double>>([&](std::size_t Index){
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

    ///@}

}; // namespace StrategiesUtilities
}  // namespace Kratos
