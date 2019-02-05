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

#if !defined(KRATOS_RESIDUAL_DISPLACEMENT_AND_OTHER_DOF_CRITERIA )
#define  KRATOS_RESIDUAL_DISPLACEMENT_AND_OTHER_DOF_CRITERIA

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

/** Short class definition.
Detail class definition.

\URL[Example of use html]{ extended_documentation/no_ex_of_use.html}

\URL[Example of use pdf]{ extended_documentation/no_ex_of_use.pdf}

\URL[Example of use doc]{ extended_documentation/no_ex_of_use.doc}

\URL[Example of use ps]{ extended_documentation/no_ex_of_use.ps}

\URL[Extended documentation html]{ extended_documentation/no_ext_doc.html}

\URL[Extended documentation pdf]{ extended_documentation/no_ext_doc.pdf}

\URL[Extended documentation doc]{ extended_documentation/no_ext_doc.doc}

\URL[Extended documentation ps]{ extended_documentation/no_ext_doc.ps}
*/

template<class TSparseSpace,
         class TDenseSpace
         >
class ResidualDisplacementAndOtherDoFCriteria : public ConvergenceCriteria< TSparseSpace, TDenseSpace >
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
        std::string OtherDoFName = "ROTATION"
        )
        : ConvergenceCriteria< TSparseSpace, TDenseSpace >(),
          mOtherDoFName(OtherDoFName),
          mRatioTolerance(RatioTolerance),
          mAbsoluteTolerance(AbsoluteTolerance)
    {
        mInitialResidualIsSet = false;
    }

    //* Copy constructor.

    ResidualDisplacementAndOtherDoFCriteria( ResidualDisplacementAndOtherDoFCriteria const& rOther )
      :BaseType(rOther)
      ,mOtherDoFName(rOther.mOtherDoFName)
      ,mInitialResidualIsSet(rOther.mInitialResidualIsSet)
      ,mRatioTolerance(rOther.mRatioTolerance)
      ,mAbsoluteTolerance(rOther.mAbsoluteTolerance)
      ,mInitialResidualDispNorm(rOther.mInitialResidualDispNorm)
      ,mCurrentResidualDispNorm(rOther.mCurrentResidualDispNorm)
      ,mInitialResidualOtherDoFNorm(rOther.mInitialResidualOtherDoFNorm)
      ,mCurrentResidualOtherDoFNorm(rOther.mCurrentResidualOtherDoFNorm)
      ,mReferenceResidualDispNorm(rOther.mReferenceResidualDispNorm)
      ,mReferenceResidualOtherDoFNorm(rOther.mReferenceResidualOtherDoFNorm)
    {
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
     * @param A System matrix (unused)
     * @param Dx Vector of results (variations on nodal variables)
     * @param b RHS vector (residual)
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
        if (TSparseSpace::Size(b) != 0) //if we are solving for something
        {

            if (mInitialResidualIsSet == false)
            {
                CalculateResidualNorm(rDofSet, b, mInitialResidualDispNorm, mInitialResidualOtherDoFNorm);
                mInitialResidualIsSet = true;
            }

            TDataType RatioDisplacement = 0.0;
            TDataType RatioOtherDoF     = 0.0;

            const std::size_t DispSize = CalculateResidualNorm(rDofSet, b, mCurrentResidualDispNorm, mCurrentResidualOtherDoFNorm);

            if(mInitialResidualDispNorm == 0.0)
            {
                RatioDisplacement = 0.0;
            }
            else
            {
                RatioDisplacement = mCurrentResidualDispNorm/mInitialResidualDispNorm;
            }

            if(mInitialResidualOtherDoFNorm == 0.0)
            {
                RatioOtherDoF = 0.0;
            }
            else
            {
                RatioOtherDoF = mCurrentResidualOtherDoFNorm/mInitialResidualOtherDoFNorm;
            }

            const std::size_t SystemSize = TSparseSpace::Size(b);
            const TDataType AbsoluteNormDisp     = (mCurrentResidualDispNorm/static_cast<TDataType>(DispSize));
            const TDataType AbsoluteNormOtherDoF = (mCurrentResidualOtherDoFNorm/static_cast<TDataType>(SystemSize - DispSize));

            if (rModelPart.GetCommunicator().MyPID() == 0 && this->GetEchoLevel() > 0)
            {
                std::cout << "RESIDUAL DISPLACEMENT CRITERION :: Ratio = "<< RatioDisplacement  << ";  Norm = " << AbsoluteNormDisp << std::endl;
                std::cout << "RESIDUAL " << mOtherDoFName << " CRITERION :: Ratio = "<< RatioOtherDoF  << ";  Norm = " << AbsoluteNormOtherDoF << std::endl;
            }

            rModelPart.GetProcessInfo()[CONVERGENCE_RATIO] = RatioDisplacement;
            rModelPart.GetProcessInfo()[RESIDUAL_NORM] = AbsoluteNormDisp;

            if ((RatioDisplacement <= mRatioTolerance || AbsoluteNormDisp < mAbsoluteTolerance) && (RatioOtherDoF <= mRatioTolerance || AbsoluteNormOtherDoF < mAbsoluteTolerance))
            {
                if (rModelPart.GetCommunicator().MyPID() == 0 && this->GetEchoLevel() > 0)
                {
                    std::cout << "Convergence is achieved" << std::endl;
                }
                return true;
            }
            else
            {
                return false;
            }
        }
        else
        {
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
        const TSystemMatrixType& A,
        const TSystemVectorType& Dx,
        const TSystemVectorType& b
    ) override
    {
        mInitialResidualIsSet = false;
    }

    /**
     * This function finalizes the solution step
     * @param rModelPart Reference to the ModelPart containing the contact problem.
     * @param rDofSet Reference to the container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param A System matrix (unused)
     * @param Dx Vector of results (variations on nodal variables)
     * @param b RHS vector (residual)
     */

    void FinalizeSolutionStep(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& A,
        const TSystemVectorType& Dx,
        const TSystemVectorType& b
    ) override
    {

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

    bool mInitialResidualIsSet;               // The flag to check if the residual has been initialized

    TDataType mInitialResidualDispNorm;       // The initial residual norm for displacements
    TDataType mCurrentResidualDispNorm;       // The current residual norm for displacements
    TDataType mInitialResidualOtherDoFNorm;   // The initial residual norm for displacements
    TDataType mCurrentResidualOtherDoFNorm;   // The current residual norm for displacements

    TDataType mRatioTolerance;                // The tolerance admited in the ratio
    TDataType mAbsoluteTolerance;             // The tolerance admited in the absolutte value

    TDataType mReferenceResidualDispNorm;     // The norm of reference for the displacement
    TDataType mReferenceResidualOtherDoFNorm; // The norm of reference for the other DoF

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * This function
     * @param rDofSet Reference to the container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param b RHS vector (residual)
     * @param DispNorm: The norm of the concerning part to residual of displacement
     * @param OtherDoFNorm: The norm of the concerning part to the residual of other Dof
     * @return SizeDeltaDisp: The number of components concerning to the displacement
     */

    std::size_t CalculateResidualNorm(
        DofsArrayType& rDofSet,
        const TSystemVectorType& b,
        TDataType& DispNorm,
        TDataType& OtherDoFNorm
        )
    {
        // Initialize variables
        std::size_t SizeDisp = 0;

        DispNorm     = 0.0;
        OtherDoFNorm = 0.0;

        TDataType AuxValue;

        for(typename DofsArrayType::iterator itDof = rDofSet.begin() ; itDof != rDofSet.end() ; ++itDof)
        {
            std::size_t DofId;

            if(itDof->IsFree())
            {
                DofId = itDof->EquationId();
                AuxValue = b[DofId];
                if (itDof->GetVariable() == DISPLACEMENT_X || itDof->GetVariable() == DISPLACEMENT_Y || itDof->GetVariable() == DISPLACEMENT_Z)
                {
                    DispNorm += AuxValue*AuxValue;
                    SizeDisp += 1;
                }
                else
                {
                    OtherDoFNorm += AuxValue*AuxValue;
                }
            }
        }

        DispNorm     = std::sqrt(DispNorm);
        OtherDoFNorm = std::sqrt(OtherDoFNorm);

        return SizeDisp;
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

#endif // KRATOS_NEW_DISPLACEMENT_CRITERIA  defined

