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
#include <unordered_map>
#include <cmath>

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
class KRATOS_API(RANS_APPLICATION) GenericConvergenceCriteria : public ConvergenceCriteria<TSparseSpace, TDenseSpace>
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
    GenericConvergenceCriteria(double rRatioTolerance, double rAbsTolerance)
        : BaseType(), mRatioTolerance(rRatioTolerance), mAbsTolerance(rAbsTolerance)
    {
        mSolutionVariables = "";
    }

    /// Destructor.
    ~GenericConvergenceCriteria() override = default;

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
                      const TSystemVectorType& b) override
    {
        if (SparseSpaceType::Size(Dx) != 0) // if we are solving for something
        {
            double solution_norm, increase_norm, dof_size;

            this->CalculateConvergenceCheckNorms(solution_norm, increase_norm, dof_size,
                                                 rModelPart, rDofSet, A, Dx, b);

            if (solution_norm == 0.0)
                solution_norm = 1.0;

            const double ratio = increase_norm / solution_norm;
            const double ratio_abs = increase_norm / dof_size;

            const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();
            const unsigned int iteration = r_current_process_info[NL_ITERATION_NUMBER];

            if (this->GetEchoLevel() > 0)
            {
                std::stringstream msg;
                msg << "[" << iteration << "] CONVERGENCE CHECK: ";
                msg << mSolutionVariables;
                msg << ": ratio = " << std::scientific << ratio
                    << "; exp.ratio = " << std::scientific << mRatioTolerance;
                msg << ": abs = " << std::scientific << ratio_abs
                    << "; exp.abs = " << std::scientific << mAbsTolerance << std::endl;
                KRATOS_INFO(this->Info()) << msg.str();
            }

            if ((std::abs(ratio) > mRatioTolerance) && (std::abs(ratio_abs) > mAbsTolerance))
                return false;

            KRATOS_INFO_IF(this->Info(), this->GetEchoLevel() > 0)
                << "CONVERGENCE CHECK: " << mSolutionVariables
                << ": *** CONVERGENCE IS ACHIEVED ***\n";

            return true;
        }
        else // in this case all the displacements are imposed!
        {
            return true;
        }
    }

    /// Initialize this class before using it
    /**
     * @param rModelPart Reference to the ModelPart containing the fluid problem. (unused)
     */
    void Initialize(ModelPart& rModelPart) override
    {
        BaseType::mConvergenceCriteriaIsInitialized = true;
    }

    void InitializeSolutionStep(ModelPart& rModelPart,
                                DofsArrayType& rDofSet,
                                const TSystemMatrixType& A,
                                const TSystemVectorType& Dx,
                                const TSystemVectorType& b) override
    {
        if (mSolutionVariables == "")
        {
            const int number_of_dofs = rDofSet.size();

            // create the variable list
            std::unordered_map<std::string, int> variables_map;
            for (int i_dof = 0; i_dof < number_of_dofs; ++i_dof)
            {
                variables_map.insert(std::make_pair(
                    (rDofSet.begin() + i_dof)->GetVariable().Name(), 1));
            }

            for (const auto& pair : variables_map)
            {
                mSolutionVariables += pair.first + ", ";
            }

            // remove last two characters
            mSolutionVariables.pop_back();
            mSolutionVariables.pop_back();
        }
    }

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "GenericConvergenceCriteria";
    }

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
