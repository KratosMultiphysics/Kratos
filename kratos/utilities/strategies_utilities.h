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
                return GetDiagonalNorm(rA)/static_cast<double>(rA.size1());
            case SCALING_DIAGONAL::CONSIDER_MAX_DIAGONAL:
                return GetMaxDiagonal(rA);
//                 return TSparseSpace::TwoNorm(rA)/static_cast<double>(rA.size1());
            default:
                return GetMaxDiagonal(rA);
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
        return 0.5 * (GetMaxDiagonal(rA) + GetMinDiagonal(rA));
    }

    /**
     * @brief This method returns the diagonal max value
     * @param rA The LHS matrix
     * @return The diagonal  max value
     */
    template<class TSparseSpace>
    static double GetMaxDiagonal(typename TSparseSpace::MatrixType& rA)
    {
//         // NOTE: Reduction failing in MSVC
//         double max_diag = 0.0;
//         #pragma omp parallel for reduction(max:max_diag)
//         for(int i = 0; i < static_cast<int>(TSparseSpace::Size1(rA)); ++i) {
//             max_diag = std::max(max_diag, std::abs(rA(i,i)));
//         }
//         return max_diag;

        // Creating a buffer for parallel vector fill
        const int num_threads = ParallelUtilities::GetNumThreads();
        Vector max_vector(num_threads, 0.0);

        #pragma omp parallel for
        for(int i = 0; i < static_cast<int>(TSparseSpace::Size1(rA)); ++i) {
            const int id = OpenMPUtils::ThisThread();
            const double abs_value_ii = std::abs(rA(i,i));
            if (abs_value_ii > max_vector[id])
                max_vector[id] = abs_value_ii;
        }

        double max_diag = 0.0;
        for(int i = 0; i < num_threads; ++i) {
            max_diag = std::max(max_diag, max_vector[i]);
        }
        return max_diag;
    }

    /**
     * @brief This method returns the diagonal min value
     * @param rA The LHS matrix
     * @return The diagonal min value
     */
    template<class TSparseSpace>
    static double GetMinDiagonal(typename TSparseSpace::MatrixType& rA)
    {
//         // NOTE: Reduction failing in MSVC
//         double min_diag = std::numeric_limits<double>::max();
//         #pragma omp parallel for reduction(min:min_diag)
//         for(int i = 0; i < static_cast<int>(TSparseSpace::Size1(rA)); ++i) {
//             min_diag = std::min(min_diag, std::abs(rA(i,i)));
//         }
//         return min_diag;

        // Creating a buffer for parallel vector fill
        const int num_threads = ParallelUtilities::GetNumThreads();
        Vector min_vector(num_threads, std::numeric_limits<double>::max());

        #pragma omp parallel for
        for(int i = 0; i < static_cast<int>(TSparseSpace::Size1(rA)); ++i) {
            const int id = OpenMPUtils::ThisThread();
            const double abs_value_ii = std::abs(rA(i,i));
            if (abs_value_ii < min_vector[id])
                min_vector[id] = abs_value_ii;
        }

        double min_diag = std::numeric_limits<double>::max();
        for(int i = 0; i < num_threads; ++i) {
            min_diag = std::min(min_diag, min_vector[i]);
        }
        return min_diag;
    }

    ///@}

}; // namespace StrategiesUtilities
}  // namespace Kratos
#endif /* KRATOS_STRATEGIES_UTILITIES defined */
