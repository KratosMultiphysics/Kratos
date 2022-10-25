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
        );

    /**
     * @brief This method returns the scale norm considering for scaling the diagonal
     * @param rModelPart The problem model part
     * @param rA The LHS matrix
     * @param ScalingDiagonal The type of caling diagonal considered
     * @return The scale norm
     */
    template<class TSparseSpace>
    static double GetScaleNorm(
        ModelPart& rModelPart,
        typename TSparseSpace::MatrixType& rA,
        const SCALING_DIAGONAL ScalingDiagonal = SCALING_DIAGONAL::NO_SCALING
        )
    {
        switch (ScalingDiagonal) {
            case SCALING_DIAGONAL::NO_SCALING:
                return 1.0;
            case SCALING_DIAGONAL::CONSIDER_PRESCRIBED_DIAGONAL: {
                const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();
                KRATOS_ERROR_IF_NOT(r_current_process_info.Has(BUILD_SCALE_FACTOR)) << "Scale factor not defined at process info" << std::endl;
                return r_current_process_info.GetValue(BUILD_SCALE_FACTOR);
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
    static double GetDiagonalNorm(typename TSparseSpace::MatrixType& rA);

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
    static double GetMaxDiagonal(typename TSparseSpace::MatrixType& rA);

    /**
     * @brief This method returns the diagonal min value
     * @param rA The LHS matrix
     * @return The diagonal min value
     */
    template<class TSparseSpace>
    static double GetMinDiagonal(typename TSparseSpace::MatrixType& rA);

    ///@}

}; // Class StrategiesUtilities
}  // namespace Kratos
