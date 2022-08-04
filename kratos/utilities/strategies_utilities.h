//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_STRATEGIES_UTILITIES)
#define KRATOS_STRATEGIES_UTILITIES

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
     * @brief This method returns the scale norm considering for scaling the diagonal
     * @param rModelPart The problem model part
     * @param rA The LHS matrix
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
    static double GetDiagonalNorm(typename TSparseSpace::MatrixType& rA)
    {
        double diagonal_norm = 0.0;
        diagonal_norm = IndexPartition<std::size_t>(TSparseSpace::Size1(rA)).for_each<SumReduction<double>>([&](std::size_t Index){
            return std::pow(rA(Index,Index), 2);
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
        const std::size_t size = TSparseSpace::Size1(rA);
        std::vector<double> abs_diagonal_values(size);
        IndexPartition<std::size_t>(size).for_each([&](std::size_t i) {
            abs_diagonal_values[i] = std::abs(rA(i,i));
        });
        return block_for_each<MaxReduction<double>>(abs_diagonal_values, [&](double& rValue) { return rValue; });
    }

    /**
     * @brief This method returns the diagonal min value
     * @param rA The LHS matrix
     * @return The diagonal min value
     */
    template<class TSparseSpace>
    static double GetMinDiagonal(typename TSparseSpace::MatrixType& rA)
    {
        const std::size_t size = TSparseSpace::Size1(rA);
        std::vector<double> abs_diagonal_values(size);
        IndexPartition<std::size_t>(size).for_each([&](std::size_t i) {
            abs_diagonal_values[i] = std::abs(rA(i,i));
        });
        return block_for_each<MinReduction<double>>(abs_diagonal_values, [&](double& rValue) { return rValue; });
    }

    ///@}

}; // namespace StrategiesUtilities
}  // namespace Kratos
#endif /* KRATOS_STRATEGIES_UTILITIES defined */
