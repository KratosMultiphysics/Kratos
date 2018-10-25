//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//
//

#if !defined(KRATOS_DISPLACEMENT_CRITERIA )
#define  KRATOS_DISPLACEMENT_CRITERIA

/* System includes */

/* External includes */

/* Project includes */
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
 * @class DisplacementCriteria
 * @ingroup KratosCore
 * @brief This is a convergence criteria that employes the increment on the solution as criteria
 * @details The reactions from the RHS are not computed in the solution
 * @author Riccardo Rossi
*/
template<class TSparseSpace,
         class TDenseSpace
         >
class DisplacementCriteria
    : public ConvergenceCriteria< TSparseSpace, TDenseSpace >
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION( DisplacementCriteria );

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
    */
    explicit DisplacementCriteria(Kratos::Parameters Settings)
        : BaseType()
    {
        if (Settings.Has("displacement_absolute_tolerance")) {
            mAlwaysConvergedNorm = Settings["displacement_absolute_tolerance"].GetDouble();
        } else if (Settings.Has("absolute_tolerance")) {
            mAlwaysConvergedNorm = Settings["absolute_tolerance"].GetDouble();
        } else {
            KRATOS_WARNING("DisplacementCriteria") << "displacement_absolute_tolerance or absolute_tolerance nor defined on settings. Using default 1.0e-9" << std::endl;
            mAlwaysConvergedNorm = 1.0e-9;
        }
        if (Settings.Has("displacement_relative_tolerance")) {
            mRatioTolerance = Settings["displacement_relative_tolerance"].GetDouble();
        } else if (Settings.Has("relative_tolerance")) {
            mRatioTolerance = Settings["relative_tolerance"].GetDouble();
        } else {
            KRATOS_WARNING("DisplacementCriteria") << "displacement_relative_tolerance or relative_tolerance nor defined on settings. Using default 1.0e-4" << std::endl;
            mRatioTolerance = 1.0e-4;
        }
    }

    /** Constructor.
    */
    explicit DisplacementCriteria(
        TDataType NewRatioTolerance,
        TDataType AlwaysConvergedNorm)
        : BaseType(),
          mAlwaysConvergedNorm(AlwaysConvergedNorm),
          mRatioTolerance(NewRatioTolerance)
    {
    }

    /** Copy constructor.
    */
    explicit DisplacementCriteria( DisplacementCriteria const& rOther )
      :BaseType(rOther)
      ,mRatioTolerance(rOther.mRatioTolerance)
      ,mAlwaysConvergedNorm(rOther.mAlwaysConvergedNorm)
      ,mReferenceDispNorm(rOther.mReferenceDispNorm)
    {
    }

    /** Destructor.
    */
    ~DisplacementCriteria() override {}


    ///@}
    ///@name Operators
    ///@{

    /**
     * Compute relative and absolute error.
     * @param rModelPart Reference to the ModelPart containing the problem.
     * @param rDofSet Reference to the container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param A System matrix (unused)
     * @param Dx Vector of results (variations on nodal variables)
     * @param b RHS vector (residual + reactions)
     * @return true if convergence is achieved, false otherwise
     */
    bool PostCriteria(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& A,
        const TSystemVectorType& Dx,
        const TSystemVectorType& b
        ) override
    {
        const TDataType approx_zero_tolerance = std::numeric_limits<TDataType>::epsilon();
        const SizeType size_Dx = Dx.size();
        if (size_Dx != 0) { //if we are solving for something
            SizeType size_solution;
            TDataType final_correction_norm = CalculateFinalCorrectionNorm(size_solution, rDofSet, Dx);

            TDataType ratio = 0.0;

            CalculateReferenceNorm(rDofSet);
            if (mReferenceDispNorm < approx_zero_tolerance) {
                KRATOS_WARNING("DisplacementCriteria") << "NaN norm is detected. Setting reference to convergence criteria" << std::endl;
                mReferenceDispNorm = final_correction_norm;
            }

            if(final_correction_norm < approx_zero_tolerance) {
                ratio = 0.0;
            } else {
                ratio = final_correction_norm/mReferenceDispNorm;
            }

            const TDataType float_size_solution = static_cast<TDataType>(size_solution);

            const TDataType absolute_norm = (final_correction_norm/std::sqrt(float_size_solution));

            KRATOS_INFO_IF("DISPLACEMENT CRITERION", this->GetEchoLevel() > 0 && rModelPart.GetCommunicator().MyPID() == 0) << " :: [ Obtained ratio = " << ratio << "; Expected ratio = " << mRatioTolerance << "; Absolute norm = " << absolute_norm << "; Expected norm =  " << mAlwaysConvergedNorm << "]" << std::endl;

            rModelPart.GetProcessInfo()[CONVERGENCE_RATIO] = ratio;
            rModelPart.GetProcessInfo()[RESIDUAL_NORM] = absolute_norm;

            if ( ratio <= mRatioTolerance  ||  absolute_norm<mAlwaysConvergedNorm )  { //  || (final_correction_norm/x.size())<=1e-7)
                KRATOS_INFO_IF("DISPLACEMENT CRITERION", this->GetEchoLevel() > 0 && rModelPart.GetCommunicator().MyPID() == 0) << "Convergence is achieved" << std::endl;
                return true;
            } else {
                return false;
            }
        } else { //in this case all the displacements are imposed!
            return true;
        }
    }

    /**
     * This function initialize the convergence criteria
     * @param rModelPart Reference to the ModelPart containing the problem. (unused)
     */
    void Initialize(
        ModelPart& rModelPart
        ) override
    {
        BaseType::mConvergenceCriteriaIsInitialized = true;
    }

    /**
     * This function initializes the solution step
     * @param rModelPart Reference to the ModelPart containing the problem.
     * @param rDofSet Reference to the container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param A System matrix (unused)
     * @param Dx Vector of results (variations on nodal variables)
     * @param b RHS vector (residual + reactions)
     */
    void InitializeSolutionStep(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& A,
        const TSystemVectorType& Dx,
        const TSystemVectorType& b
        ) override
    {
        BaseType::InitializeSolutionStep(rModelPart, rDofSet, A, Dx, b);
    }

    /**
     * This function finalizes the solution step
     * @param rModelPart Reference to the ModelPart containing the problem.
     * @param rDofSet Reference to the container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param A System matrix (unused)
     * @param Dx Vector of results (variations on nodal variables)
     * @param b RHS vector (residual + reactions)
     */
    void FinalizeSolutionStep(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& A,
        const TSystemVectorType& Dx,
        const TSystemVectorType& b
        ) override
    {
        BaseType::FinalizeSolutionStep(rModelPart, rDofSet, A, Dx, b);
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

    TDataType mRatioTolerance;      /// The ratio threshold for the norm of the residual

    TDataType mAlwaysConvergedNorm; /// The absolute value threshold for the norm of the residual

    TDataType mReferenceDispNorm;   /// The norm at the beginning of the iterations


    ///@}
    ///@name Private Operators
    ///@{

    /**
     * @brief This method computes the reference norm
     * @details It checks if the dof is fixed
     * @param rDofSet Reference to the container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     */
    void CalculateReferenceNorm(DofsArrayType& rDofSet)
    {
        TDataType reference_disp_norm = TDataType();
        TDataType dof_value;

        #pragma omp parallel for reduction(+:reference_disp_norm)
        for (int i = 0; i < static_cast<int>(rDofSet.size()); i++) {
            auto it_dof = rDofSet.begin() + i;

            if(it_dof->IsFree()) {
                dof_value = it_dof->GetSolutionStepValue();
                reference_disp_norm += dof_value * dof_value;
            }
        }
        mReferenceDispNorm = std::sqrt(reference_disp_norm);
    }

    /**
     * @brief This method computes the final norm
     * @details It checks if the dof is fixed
     * @param rDofNum The number of DoFs
     * @param rDofSet Reference to the container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param Dx Vector of results (variations on nodal variables)
     */
    TDataType CalculateFinalCorrectionNorm(
        SizeType& rDofNum,
        DofsArrayType& rDofSet,
        const TSystemVectorType& Dx
        )
    {
        // Initialize
        TDataType final_correction_norm = TDataType();
        SizeType dof_num = 0;

        // Loop over Dofs
        #pragma omp parallel for reduction(+:final_correction_norm,dof_num)
        for (int i = 0; i < static_cast<int>(rDofSet.size()); i++) {
            auto it_dof = rDofSet.begin() + i;

            IndexType dof_id;
            TDataType variation_dof_value;

            if (it_dof->IsFree()) {
                dof_id = it_dof->EquationId();
                variation_dof_value = Dx[dof_id];
                final_correction_norm += variation_dof_value * variation_dof_value;
                dof_num++;
            }
        }

        rDofNum = dof_num;
        return final_correction_norm;
    }

    ///@}
    ///@name Private Operations
    ///@{


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

}; /* Class DisplacementCriteria */

///@}

///@name Type Definitions
///@{


///@}

}  /* namespace Kratos.*/

#endif /* KRATOS_DISPLACEMENT_CRITERIA  defined */

