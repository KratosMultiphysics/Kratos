// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/model_part.h"
#include "includes/define.h"
#include "solving_strategies/convergencecriterias/convergence_criteria.h"

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

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
 * @class ResidualDisplacementAndOtherDoFCriteria
 * @ingroup StructuralMechanicsApplication
 * @brief This is a convergence criteria that employes the residual as criteria, it divides the problem in two dofs, displacement and another one
 * @details The reactions from the RHS are not computed in the residual. Can be used with for example rotations or pressure
 * @author Vicente Mataix Ferrandiz
*/
template<class TSparseSpace,
         class TDenseSpace
         >
class ResidualDisplacementAndOtherDoFCriteria
    : public ConvergenceCriteria< TSparseSpace, TDenseSpace >
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION( ResidualDisplacementAndOtherDoFCriteria );

    typedef ConvergenceCriteria< TSparseSpace, TDenseSpace > BaseType;

    typedef TSparseSpace SparseSpaceType;

    typedef typename BaseType::TDataType TDataType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    typedef std::size_t IndexType;

    typedef std::size_t SizeType;

    ///@}
    ///@name Life Cycle

    ///@{

    /** Constructor.
     * @param RatioTolerance: Relative tolerance for error
     * @param AbsoluteTolerance: Absolute tolerance for error
     * @param OtherDoFName: The name of the other DoF
     */
    ResidualDisplacementAndOtherDoFCriteria(
        TDataType RatioTolerance,
        TDataType AbsoluteTolerance,
        const std::string& OtherDoFName = "ROTATION"
        )
        : ConvergenceCriteria< TSparseSpace, TDenseSpace >(),
          mOtherDoFName(OtherDoFName),
          mRatioTolerance(RatioTolerance),
          mAbsoluteTolerance(AbsoluteTolerance)
    {
        this->mActualizeRHSIsNeeded = true;
    }

    //* Copy constructor.

    ResidualDisplacementAndOtherDoFCriteria( ResidualDisplacementAndOtherDoFCriteria const& rOther )
      :BaseType(rOther)
      ,mOtherDoFName(rOther.mOtherDoFName)
      ,mRatioTolerance(rOther.mRatioTolerance)
      ,mAbsoluteTolerance(rOther.mAbsoluteTolerance)
      ,mInitialResidualDispNorm(rOther.mInitialResidualDispNorm)
      ,mCurrentResidualDispNorm(rOther.mCurrentResidualDispNorm)
      ,mInitialResidualOtherDoFNorm(rOther.mInitialResidualOtherDoFNorm)
      ,mCurrentResidualOtherDoFNorm(rOther.mCurrentResidualOtherDoFNorm)
    {
        this->mActualizeRHSIsNeeded = true;
    }

    //* Destructor.

    ~ResidualDisplacementAndOtherDoFCriteria() override {}


    ///@}
    ///@name Operators

    ///@{

    /**
     * Compute relative and absolute error.
     * @param rModelPart Reference to the ModelPart containing the contact problem.
     * @param rDofSet Reference to the container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param rA System matrix (unused)
     * @param rDx Vector of results (variations on nodal variables)
     * @param rb RHS vector (residual)
     * @return true if convergence is achieved, false otherwise
     */

    bool PostCriteria(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& rA,
        const TSystemVectorType& rDx,
        const TSystemVectorType& rb
    ) override
    {
        if (TSparseSpace::Size(rb) != 0) { //if we are solving for something
            TDataType ratio_displacement = 0.0;
            TDataType ratio_other_dof     = 0.0;

            SizeType disp_size;
            CalculateResidualNorm(rModelPart, mCurrentResidualDispNorm, mCurrentResidualOtherDoFNorm, disp_size, rDofSet, rb);

            if (mInitialResidualDispNorm == 0.0) {
                ratio_displacement = 0.0;
            } else {
                ratio_displacement = mCurrentResidualDispNorm/mInitialResidualDispNorm;
            }

            if (mInitialResidualOtherDoFNorm == 0.0) {
                ratio_other_dof = 0.0;
            } else {
                ratio_other_dof = mCurrentResidualOtherDoFNorm/mInitialResidualOtherDoFNorm;
            }

            const std::size_t system_size = TSparseSpace::Size(rb);
            const TDataType absolute_norm_disp      = (mCurrentResidualDispNorm/static_cast<TDataType>(disp_size));
            const TDataType absolute_norm_other_dof = (mCurrentResidualOtherDoFNorm/static_cast<TDataType>(system_size - disp_size));

            KRATOS_INFO_IF("ResidualDisplacementAndOtherDoFCriteria", this->GetEchoLevel() > 0) << "RESIDUAL DISPLACEMENT CRITERION :: Ratio = "<< ratio_displacement  << ";  Norm = " << absolute_norm_disp << std::endl;
            KRATOS_INFO_IF("ResidualDisplacementAndOtherDoFCriteria", this->GetEchoLevel() > 0) << "RESIDUAL " << mOtherDoFName << " CRITERION :: Ratio = "<< ratio_other_dof  << ";  Norm = " << absolute_norm_other_dof << std::endl;

            rModelPart.GetProcessInfo()[CONVERGENCE_RATIO] = ratio_displacement;
            rModelPart.GetProcessInfo()[RESIDUAL_NORM] = absolute_norm_disp;

            if ((ratio_displacement <= mRatioTolerance || absolute_norm_disp < mAbsoluteTolerance) && (ratio_other_dof <= mRatioTolerance || absolute_norm_other_dof < mAbsoluteTolerance)) {
                KRATOS_INFO_IF("ResidualDisplacementAndOtherDoFCriteria", this->GetEchoLevel() > 0) << "Convergence is achieved" << std::endl;
                return true;
            } else {
                return false;
            }
        } else {
            return true;
        }
    }

    /**
     * This function initialize the convergence criteria
     * @param rModelPart Reference to the ModelPart containing the contact problem. (unused)
     */

    void Initialize(
        ModelPart& rModelPart
    ) override
    {
        BaseType::mConvergenceCriteriaIsInitialized = true;
    }

    /**
     * This function initializes the solution step
     * @param rModelPart Reference to the ModelPart containing the contact problem.
     * @param rDofSet Reference to the container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param A System matrix (unused)
     * @param Dx Vector of results (variations on nodal variables)
     * @param b RHS vector (residual)
     */

    void InitializeSolutionStep(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& rA,
        const TSystemVectorType& rDx,
        const TSystemVectorType& rb
    ) override
    {
        BaseType::InitializeSolutionStep(rModelPart, rDofSet, rA, rDx, rb);

        // Filling mActiveDofs when MPC exist
        if (rModelPart.NumberOfMasterSlaveConstraints() > 0) {
            mActiveDofs.resize(rDofSet.size());

            #pragma omp parallel for
            for(int i=0; i<static_cast<int>(mActiveDofs.size()); ++i) {
                mActiveDofs[i] = true;
            }

            #pragma omp parallel for
            for (int i=0; i<static_cast<int>(rDofSet.size()); ++i) {
                const auto it_dof = rDofSet.begin() + i;
                if (it_dof->IsFixed()) {
                    mActiveDofs[it_dof->EquationId()] = false;
                }
            }

            for (const auto& r_mpc : rModelPart.MasterSlaveConstraints()) {
                for (const auto& r_dof : r_mpc.GetMasterDofsVector()) {
                    mActiveDofs[r_dof->EquationId()] = false;
                }
                for (const auto& r_dof : r_mpc.GetSlaveDofsVector()) {
                    mActiveDofs[r_dof->EquationId()] = false;
                }
            }
        }

        SizeType size_residual;
        CalculateResidualNorm(rModelPart, mInitialResidualDispNorm, mInitialResidualOtherDoFNorm, size_residual, rDofSet, rb);
    }

    ///@}
    ///@name Operations
    ///@{

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    std::string mOtherDoFName;                // The name of the other DoF

    TDataType mInitialResidualDispNorm;       // The initial residual norm for displacements
    TDataType mCurrentResidualDispNorm;       // The current residual norm for displacements
    TDataType mInitialResidualOtherDoFNorm;   // The initial residual norm for displacements
    TDataType mCurrentResidualOtherDoFNorm;   // The current residual norm for displacements

    TDataType mRatioTolerance;                // The tolerance admited in the ratio
    TDataType mAbsoluteTolerance;             // The tolerance admited in the absolutte value

    std::vector<bool> mActiveDofs;  /// This vector contains the dofs that are active

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief This method computes the norm of the residual
     * @details It checks if the dof is fixed
     * @param rModelPart Reference to the ModelPart containing the problem.
     * @param rResidualSolutionNorm The norm of the residual
     * @param rDofNum The number of DoFs
     * @param rDofSet Reference to the container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param rb RHS vector (residual + reactions)
     */
    virtual void CalculateResidualNorm(
        ModelPart& rModelPart,
        TDataType& rResidualSolutionNormDisp,
        TDataType& rResidualSolutionNormOtherDof,
        SizeType& rDofNumDisp,
        DofsArrayType& rDofSet,
        const TSystemVectorType& rb
        )
    {
        // Initialize
        TDataType residual_solution_norm_disp = TDataType();
        TDataType residual_solution_norm_other_dof = TDataType();
        SizeType disp_dof_num = 0;

        // Auxiliary values
        TDataType residual_dof_value = 0.0;
        const auto it_dof_begin = rDofSet.begin();
        const int number_of_dof = static_cast<int>(rDofSet.size());

        // Loop over Dofs
        if (rModelPart.NumberOfMasterSlaveConstraints() > 0) {
            #pragma omp parallel for firstprivate(residual_dof_value) reduction(+:residual_solution_norm_disp, residual_solution_norm_other_dof, disp_dof_num)
            for (int i = 0; i < number_of_dof; i++) {
                auto it_dof = it_dof_begin + i;

                const IndexType dof_id = it_dof->EquationId();
                residual_dof_value = TSparseSpace::GetValue(rb,dof_id);

                if (mActiveDofs[dof_id]) {
                    if (it_dof->GetVariable() == DISPLACEMENT_X || it_dof->GetVariable() == DISPLACEMENT_Y || it_dof->GetVariable() == DISPLACEMENT_Z) {
                        residual_solution_norm_disp += std::pow(residual_dof_value, 2);
                        disp_dof_num++;
                    } else {
                        residual_solution_norm_other_dof += std::pow(residual_dof_value, 2);
                    }
                }
            }
        } else {
            #pragma omp parallel for firstprivate(residual_dof_value) reduction(+:residual_solution_norm_disp, residual_solution_norm_other_dof, disp_dof_num)
            for (int i = 0; i < number_of_dof; i++) {
                auto it_dof = it_dof_begin + i;

                if (!it_dof->IsFixed()) {
                    const IndexType dof_id = it_dof->EquationId();
                    residual_dof_value = TSparseSpace::GetValue(rb,dof_id);

                    if (it_dof->GetVariable() == DISPLACEMENT_X || it_dof->GetVariable() == DISPLACEMENT_Y || it_dof->GetVariable() == DISPLACEMENT_Z) {
                        residual_solution_norm_disp += std::pow(residual_dof_value, 2);
                        disp_dof_num++;
                    } else {
                        residual_solution_norm_other_dof += std::pow(residual_dof_value, 2);
                    }
                }
            }
        }

        rDofNumDisp = disp_dof_num;
        rResidualSolutionNormDisp = std::sqrt(residual_solution_norm_disp);
        rResidualSolutionNormOtherDof = std::sqrt(residual_solution_norm_other_dof);
    }

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

}; // Class ClassName

///@}

///@name Type Definitions
///@{


///@}

}  // namespace Kratos.

