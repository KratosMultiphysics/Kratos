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

#if !defined(KRATOS_NEW_CONVERGENCE_CRITERIA )
#define  KRATOS_NEW_CONVERGENCE_CRITERIA


/* System includes */


/* External includes */


/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"

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
         class TDenseSpace //= DenseSpace<double>
         >
class ConvergenceCriteria
{
public:
    ///@name Type Definitions
    ///@{

    typedef typename TSparseSpace::DataType TDataType;
    typedef typename TSparseSpace::MatrixType TSystemMatrixType;
    typedef typename TSparseSpace::VectorType TSystemVectorType;

    typedef typename TDenseSpace::MatrixType LocalSystemMatrixType;
    typedef typename TDenseSpace::VectorType LocalSystemVectorType;

    //typedef Dof<double> TDofType;
    typedef ModelPart::DofsArrayType DofsArrayType;

    /** Counted pointer of ConvergenceCriteria */
    KRATOS_CLASS_POINTER_DEFINITION(ConvergenceCriteria);
    ///@}
    ///@name Life Cycle
    ///@{

    /** Constructor.
     */
    ConvergenceCriteria()
    {
        mActualizeRHSIsNeeded = false;
        mConvergenceCriteriaIsInitialized = false;
        SetEchoLevel(1);
    }

    /** Copy constructor.
     */
    ConvergenceCriteria( ConvergenceCriteria const& rOther)
      :mActualizeRHSIsNeeded(rOther.mActualizeRHSIsNeeded)
      ,mConvergenceCriteriaIsInitialized(rOther.mConvergenceCriteriaIsInitialized)
      ,mEchoLevel(rOther.mEchoLevel)
    {
    }

    /** Destructor.
     */
    virtual ~ConvergenceCriteria()
    {
    }

    ///@}
    ///@name Member Variables
    ///@{



    ///@}
    ///@name Operators
    ///@{

    /**
     * Get component wise element components
     */
    virtual std::vector<TSystemVectorType>&  GetRHS_Element_Components()
    {
      KRATOS_ERROR <<"Asking for Global Components to the CONVERGENCE CRITERION base class which is not component wise and not contains this member variable" << std::endl;
    }

    /**
     * Get component wise element variables
     */
    virtual std::vector< Variable< LocalSystemVectorType > >&  GetRHS_Element_Variables()
    {
      KRATOS_ERROR <<"Asking for Global Components to the CONVERGENCE CRITERION base class which is not component wise and not contains this member variable" << std::endl;
    }

    /**
     * Get component wise condition components
     */
    virtual std::vector<TSystemVectorType>&  GetRHS_Condition_Components()
    {
      KRATOS_ERROR <<"Asking for Global Components to the CONVERGENCE CRITERION base class which is not component wise and not contains this member variable" << std::endl;
    }

    /**
     * Get component wise condition variables
     */
    virtual std::vector< Variable< LocalSystemVectorType > >&  GetRHS_Condition_Variables()
    {
      KRATOS_ERROR <<"Asking for Global Components to the CONVERGENCE CRITERION base class which is not component wise and not contains this member variable" << std::endl;
    }

    //*********************************************************************************

    /**level of echo for the convergence criterion
    0 -> mute... no echo at all
    1 -> print basic informations
    2 -> print extra informations
     */
    virtual void SetEchoLevel(int Level)
    {
        mEchoLevel = Level;
    }

    int GetEchoLevel()
    {
        return mEchoLevel;
    }


    void SetActualizeRHSFlag(bool flag)
    {
        mActualizeRHSIsNeeded = flag;
    }

    bool GetActualizeRHSflag()
    {
        return mActualizeRHSIsNeeded;
    }

    /*Criterias that need to be called before getting the solution */
    virtual bool PreCriteria(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& A,
        const TSystemVectorType& Dx,
        const TSystemVectorType& b
    )
    {
        return true;
    }

    /*Criterias that need to be called after getting the solution */
    virtual bool PostCriteria(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& A,
        const TSystemVectorType& Dx,
        const TSystemVectorType& b
    )
    {
        return true;
    }

    virtual void Initialize(
        ModelPart& rModelPart
        )
    {
        mConvergenceCriteriaIsInitialized = true;
    }

    virtual bool IsInitialized()
    {return mConvergenceCriteriaIsInitialized;}



    virtual void InitializeSolutionStep(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& A,
        const TSystemVectorType& Dx,
        const TSystemVectorType& b
    )
    {
    }

    virtual void InitializeNonLinearIteration(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& A,
        const TSystemVectorType& Dx,
        const TSystemVectorType& b
    )
    {
    }

    virtual void FinalizeSolutionStep(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& A,
        const TSystemVectorType& Dx,
        const TSystemVectorType& b
    )
    {
    }

    virtual void FinalizeNonLinearIteration(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& A,
        const TSystemVectorType& Dx,
        const TSystemVectorType& b
    )
    {
    }

    /**
     * This function is designed to be called once to perform all the checks needed
     * on the input provided. Checks can be "expensive" as the function is designed
     * to catch user's errors.
     * @param rModelPart
     * @return 0 all ok
     */
    virtual int Check(ModelPart& rModelPart)
    {
        KRATOS_TRY

        return 0;
        KRATOS_CATCH("");
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
    bool mActualizeRHSIsNeeded = false;
    bool mConvergenceCriteriaIsInitialized = false  ;
    int  mEchoLevel;

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

    ///@}
    ///@name Private Operators
    ///@{

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

}; /* Class ConvergenceCriteria */

///@}

///@name Type Definitions */
///@{

///@}

} /* namespace Kratos.*/

#endif /* KRATOS_NEW_CONVERGENCE_CRITERIA  defined */

