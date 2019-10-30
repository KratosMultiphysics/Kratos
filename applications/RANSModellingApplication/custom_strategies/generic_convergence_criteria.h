//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya (https://github.com/sunethwarna)
//

#ifndef KRATOS_GENERIC_CONVERGENCE_CRITERIA_H
#define KRATOS_GENERIC_CONVERGENCE_CRITERIA_H

// System includes
#include <string>

/* Project includes */
#include "includes/model_part.h"
#include "solving_strategies/convergencecriterias/convergence_criteria.h"

namespace Kratos
{
///@addtogroup IncompressibleFluidApplication
///@{

///@name Kratos Classes
///@{

/// Convergence criteria for fluid problems.
/**
 This class implements a convergence control based on nodal velocity and
 pressure values. The error is evaluated separately for each of them, and
 relative and absolute tolerances for both must be specified.
 */
template <class TSparseSpace, class TDenseSpace>
class GenericConvergenceCriteria : public ConvergenceCriteria<TSparseSpace, TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(GenericConvergenceCriteria);

    using BaseType = ConvergenceCriteria<TSparseSpace, TDenseSpace>;

    using SparseSpaceType = TSparseSpace;

    using TDataType = typename BaseType::TDataType;

    using DofsArrayType = typename BaseType::DofsArrayType;

    using TSystemMatrixType = typename BaseType::TSystemMatrixType;

    using TSystemVectorType = typename BaseType::TSystemVectorType;

    using PartitionVector = OpenMPUtils::PartitionVector;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    /**
     * @param VelRatioTolerance Relative tolerance for velocity error
     * @param VelAbsTolerance Absolute tolerance for velocity error
     * @param PrsRatioTolerance Relative tolerance for presssure error
     * @param PrsAbsTolerance Absolute tolerance for presssure error
     */
    GenericConvergenceCriteria(double rRatioTolerance, double rAbsTolerance);

    /// Destructor.
    ~GenericConvergenceCriteria() override;

    ///@}
    ///@name Operators
    ///@{

    /// Compute relative and absoute error.
    /**
     * @param rModelPart Reference to the ModelPart containing the fluid problem.
     * @param rDofSet Reference to the container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param A System matrix (unused)
     * @param Dx Vector of results (variations on nodal variables)
     * @param b RHS vector (residual)
     * @return true if convergence is achieved, false otherwise
     */
    bool PostCriteria(ModelPart& rModelPart,
                      DofsArrayType& rDofSet,
                      const TSystemMatrixType& A,
                      const TSystemVectorType& Dx,
                      const TSystemVectorType& b) override;

    /// Initialize this class before using it
    /**
     * @param rModelPart Reference to the ModelPart containing the fluid problem. (unused)
     */
    void Initialize(ModelPart& rModelPart) override;

    void InitializeSolutionStep(ModelPart& rModelPart,
                                DofsArrayType& rDofSet,
                                const TSystemMatrixType& A,
                                const TSystemVectorType& Dx,
                                const TSystemVectorType& b) override;

    /// Turn back information as a string.
    std::string Info() const override;

    ///@} // Operations

private:
    double mRatioTolerance;
    double mAbsTolerance;

    std::string mSolutionVariables;

    void CalculateConvergenceCheckNorms(double& rSolutionNorm,
                                        double& rIncreaseNorm,
                                        double& rDofSize,
                                        ModelPart& rModelPart,
                                        DofsArrayType& rDofSet,
                                        const TSystemMatrixType& A,
                                        const TSystemVectorType& Dx,
                                        const TSystemVectorType& b);

}; // namespace Kratos

///@} // Kratos classes

///@} // Application group
} // namespace Kratos

#endif /* _VEL_PR_CRITERIA_H */
