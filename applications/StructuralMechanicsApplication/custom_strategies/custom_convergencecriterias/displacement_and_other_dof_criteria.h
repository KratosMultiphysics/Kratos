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
//

#if !defined(KRATOS_DISPLACEMENT_AND_OTHER_DOF_CRITERIA )
#define  KRATOS_DISPLACEMENT_AND_OTHER_DOF_CRITERIA

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
class DisplacementAndOtherDoFCriteria : public ConvergenceCriteria< TSparseSpace, TDenseSpace >
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION( DisplacementAndOtherDoFCriteria );

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
    DisplacementAndOtherDoFCriteria(
        TDataType RatioTolerance,
        TDataType AbsoluteTolerance,
        std::string OtherDoFName = "ROTATION"
        )
        : ConvergenceCriteria< TSparseSpace, TDenseSpace >(),
          mOtherDoFName(OtherDoFName),
          mRatioTolerance(RatioTolerance),
          mAbsoluteTolerance(AbsoluteTolerance)
    {
    }

    /** Copy constructor.
    */
    DisplacementAndOtherDoFCriteria( DisplacementAndOtherDoFCriteria const& rOther )
      :BaseType(rOther)
      ,mRatioTolerance(rOther.mRatioTolerance)
      ,mAbsoluteTolerance(rOther.mAbsoluteTolerance)
      ,mReferenceDispNorm(rOther.mReferenceDispNorm)
      ,mReferenceOtherDoFNorm(rOther.mReferenceOtherDoFNorm)
      ,mOtherDoFName(rOther.mOtherDoFName)
    {
    }

    /** Destructor.
    */
    ~DisplacementAndOtherDoFCriteria() override {}


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
        if (SparseSpaceType::Size(Dx) != 0) //if we are solving for something
        {
            TDataType DeltaDisplacenentNorm = 0.0;
            TDataType DeltaOtherDoFNorm     = 0.0;

            const std::size_t SizeDx = SparseSpaceType::Size(Dx);
            const std::size_t SizeDisp = CalculateDeltaxNorm(rDofSet, Dx, DeltaDisplacenentNorm, DeltaOtherDoFNorm);

            TDataType RatioDisplacement = 0.0;
            TDataType RatioOtherDoF     = 0.0;

            CalculateReferenceNorm(rDofSet);

            if(DeltaDisplacenentNorm == 0)
            {
                RatioDisplacement = 0.0;
            }
            else
            {
                if(mReferenceDispNorm == 0.0)
                {
                    KRATOS_ERROR << "NaN norm is detected" << std::endl;
                }

                RatioDisplacement = DeltaDisplacenentNorm/mReferenceDispNorm;
            }

            if(DeltaOtherDoFNorm == 0)
            {
                RatioOtherDoF = 0.0;
            }
            else
            {
                if(mReferenceOtherDoFNorm == 0.0)
                {
                    KRATOS_ERROR << "NaN norm is detected" << std::endl;
                }

                RatioOtherDoF = DeltaOtherDoFNorm/mReferenceOtherDoFNorm;
            }

            const TDataType AbsoluteDisplacementNorm = (DeltaDisplacenentNorm/std::sqrt(static_cast<TDataType>(SizeDisp)));
            const TDataType AbsoluteOtherDoFNorm     = (DeltaOtherDoFNorm/std::sqrt(static_cast<TDataType>(SizeDx - SizeDisp)));

            if (rModelPart.GetCommunicator().MyPID() == 0 && this->GetEchoLevel() > 0)
            {
                std::cout << "DISPLACEMENT CRITERION :: [ Obtained ratio = " << RatioDisplacement << "; Expected ratio = " << mRatioTolerance << "; Absolute norm = " << AbsoluteDisplacementNorm << "; ]" << std::endl;
                std::cout << "OTHER "<< mOtherDoFName << " CRITERION :: [ Obtained ratio = " << RatioOtherDoF << "; Expected ratio = " << mRatioTolerance << "; Absolute norm = " << AbsoluteOtherDoFNorm << "; ]" << std::endl;
            }

            rModelPart.GetProcessInfo()[CONVERGENCE_RATIO] = RatioDisplacement;
            rModelPart.GetProcessInfo()[RESIDUAL_NORM] = AbsoluteDisplacementNorm;

            if ( (RatioDisplacement <= mRatioTolerance  ||  AbsoluteDisplacementNorm < mAbsoluteTolerance) || RatioOtherDoF <= mRatioTolerance  || AbsoluteOtherDoFNorm < mAbsoluteTolerance  )
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
        else //in this case all the displacements are imposed!
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

    std::string mOtherDoFName;        // The name of the other DoF

    TDataType mRatioTolerance;        // The tolerance admited in the ratio
    TDataType mAbsoluteTolerance;     // The tolerance admited in the absolutte value

    TDataType mReferenceDispNorm;     // The norm of reference for the displacement
    TDataType mReferenceOtherDoFNorm; // The norm of reference for the other DoF

    ///@}
    ///@name Private Operators
    ///@{

    /**
     * This function calculates the reference norm
     * @param rDofSet Reference to the container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     */

    void CalculateReferenceNorm(DofsArrayType& rDofSet)
    {
        mReferenceDispNorm     = TDataType();
        mReferenceOtherDoFNorm = TDataType();
        TDataType AuxValue;

        for(typename DofsArrayType::iterator itDof = rDofSet.begin() ; itDof != rDofSet.end() ; ++itDof)
        {
            if(itDof->IsFree())
            {
                AuxValue = itDof->GetSolutionStepValue();
                if (itDof->GetVariable() == DISPLACEMENT_X || itDof->GetVariable() == DISPLACEMENT_Y || itDof->GetVariable() == DISPLACEMENT_Z)
                {
                    mReferenceDispNorm += AuxValue*AuxValue;
                }
                else
                {
                    mReferenceOtherDoFNorm += AuxValue*AuxValue;
                }
            }
        }

        mReferenceDispNorm     = std::sqrt(mReferenceDispNorm);
        mReferenceOtherDoFNorm = std::sqrt(mReferenceOtherDoFNorm);
    }

    /**
     * This function
     * @param rDofSet Reference to the container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param Dx Vector of results (variations on nodal variables)
     * @param DeltaDisplacenentNorm: The norm of the concerning part to increment of displacement
     * @param DeltaOtherDoFNorm: The norm of the concerning part to the increment of other Dof
     * @return SizeDeltaDisp: The number of components concerning to the displacement
     */

    std::size_t CalculateDeltaxNorm(
        DofsArrayType& rDofSet,
        const TSystemVectorType& Dx,
        TDataType& DeltaDisplacenentNorm,
        TDataType& DeltaOtherDoFNorm
        )
    {
        std::size_t SizeDeltaDisp = 0;

        TDataType AuxValue;

        for(typename DofsArrayType::iterator itDof = rDofSet.begin() ; itDof != rDofSet.end() ; ++itDof)
        {
            std::size_t DofId;

            if(itDof->IsFree())
            {
                DofId = itDof->EquationId();
                AuxValue = Dx[DofId];
                if (itDof->GetVariable() == DISPLACEMENT_X || itDof->GetVariable() == DISPLACEMENT_Y || itDof->GetVariable() == DISPLACEMENT_Z)
                {
                    DeltaDisplacenentNorm += AuxValue*AuxValue;
                    SizeDeltaDisp += 1;
                }
                else
                {
                    mReferenceOtherDoFNorm += AuxValue*AuxValue;
                }
            }
        }

        DeltaDisplacenentNorm = std::sqrt(DeltaDisplacenentNorm);
        DeltaOtherDoFNorm     = std::sqrt(DeltaOtherDoFNorm);

        return SizeDeltaDisp;
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

}; /* Class ClassName */

///@}

///@name Type Definitions
///@{


///@}

}  /* namespace Kratos.*/

#endif /* KRATOS_DISPLACEMENT_AND_OTHER_DOF_CRITERIA  defined */

