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

#if !defined(KRATOS_NEW_RESIDUAL_CRITERIA )
#define  KRATOS_NEW_RESIDUAL_CRITERIA


/* System includes */


/* External includes */


/* Project includes */
#include "includes/model_part.h"
#include "includes/define.h"
#include "solving_strategies/convergencecriterias/convergence_criteria.h"

namespace Kratos
{

/**@name Kratos Globals */
/*@{ */


/*@} */
/**@name Type Definitions */
/*@{ */

/*@} */


/**@name  Enum's */
/*@{ */


/*@} */
/**@name  Functions */
/*@{ */



/*@} */
/**@name Kratos Classes */
/*@{ */

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
class ResidualCriteria : public virtual  ConvergenceCriteria< TSparseSpace, TDenseSpace >
{
public:
    /**@name Type Definitions */
    /*@{ */

    KRATOS_CLASS_POINTER_DEFINITION( ResidualCriteria );

    typedef ConvergenceCriteria< TSparseSpace, TDenseSpace > BaseType;

    typedef TSparseSpace SparseSpaceType;

    typedef typename BaseType::TDataType TDataType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    /*@} */
    /**@name Life Cycle
    */
    /*@{ */

    /** Constructor.
    */
    ResidualCriteria(
        TDataType NewRatioTolerance,
        TDataType AlwaysConvergedNorm)
        : ConvergenceCriteria< TSparseSpace, TDenseSpace >()
    {
        mRatioTolerance       = NewRatioTolerance;
        mAlwaysConvergedNorm  = AlwaysConvergedNorm;
        mInitialResidualIsSet = false;
    }

    /** Copy constructor.
    */
    ResidualCriteria( ResidualCriteria const& rOther )
      :BaseType(rOther) 
      ,mInitialResidualIsSet(rOther.mInitialResidualIsSet)
      ,mRatioTolerance(rOther.mRatioTolerance)
      ,mInitialResidualNorm(rOther.mInitialResidualNorm)
      ,mCurrentResidualNorm(rOther.mCurrentResidualNorm)
      ,mAlwaysConvergedNorm(rOther.mAlwaysConvergedNorm)
      ,mReferenceDispNorm(rOther.mReferenceDispNorm)
    {
    }

    /** Destructor.
    */
    virtual ~ResidualCriteria() {}


    /*@} */
    /**@name Operators
    */
    /*@{ */

    /*Criterias that need to be called after getting the solution */
    bool PostCriteria(
        ModelPart& r_model_part,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& A,
        const TSystemVectorType& Dx,
        const TSystemVectorType& b
    )
    {
        if (TSparseSpace::Size(b) != 0) //if we are solving for something
        {

            if (mInitialResidualIsSet == false)
            {
                mInitialResidualNorm = TSparseSpace::TwoNorm(b);
                mInitialResidualIsSet = true;
            }

            TDataType ratio;
            mCurrentResidualNorm = TSparseSpace::TwoNorm(b);

            double b_size = TSparseSpace::Size(b);


            if(mInitialResidualNorm == 0.00) ratio = 0.00;

            else ratio = mCurrentResidualNorm/mInitialResidualNorm;

            if (r_model_part.GetCommunicator().MyPID() == 0)
	      if (this->GetEchoLevel() == 1)
                std::cout << "RESIDUAL CRITERION :: Ratio = " << ratio  << ";  Norm   = " << mCurrentResidualNorm/b_size << std::endl;

	    r_model_part.GetProcessInfo()[CONVERGENCE_RATIO] = ratio;
	    r_model_part.GetProcessInfo()[RESIDUAL_NORM] = mCurrentResidualNorm/b_size;

            if (
                ratio <= mRatioTolerance
                ||
                (mCurrentResidualNorm/b_size) <mAlwaysConvergedNorm
            )
            {
                if (r_model_part.GetCommunicator().MyPID() == 0)
		  if (this->GetEchoLevel() == 1)
                    std::cout << "Convergence is achieved" << std::endl;
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




    void Initialize(
        ModelPart& r_model_part
    )
    {
        BaseType::mConvergenceCriteriaIsInitialized = true;
    }

    void InitializeSolutionStep(
        ModelPart& r_model_part,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& A,
        const TSystemVectorType& Dx,
        const TSystemVectorType& b
    )
    {
        mInitialResidualIsSet = false;
    }

    void FinalizeSolutionStep(
        ModelPart& r_model_part,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& A,
        const TSystemVectorType& Dx,
        const TSystemVectorType& b
    ) {}



    /*@} */
    /**@name Operations */
    /*@{ */


    /*@} */
    /**@name Access */
    /*@{ */


    /*@} */
    /**@name Inquiry */
    /*@{ */


    /*@} */
    /**@name Friends */
    /*@{ */


    /*@} */

protected:
    /**@name Protected static Member Variables */
    /*@{ */


    /*@} */
    /**@name Protected member Variables */
    /*@{ */


    /*@} */
    /**@name Protected Operators*/
    /*@{ */


    /*@} */
    /**@name Protected Operations*/
    /*@{ */


    /*@} */
    /**@name Protected  Access */
    /*@{ */


    /*@} */
    /**@name Protected Inquiry */
    /*@{ */


    /*@} */
    /**@name Protected LifeCycle */
    /*@{ */



    /*@} */

private:
    /**@name Static Member Variables */
    /*@{ */


    /*@} */
    /**@name Member Variables */
    /*@{ */


    bool mInitialResidualIsSet;

    TDataType mRatioTolerance;

    TDataType mInitialResidualNorm;

    TDataType mCurrentResidualNorm;

    TDataType mAlwaysConvergedNorm;

    TDataType mReferenceDispNorm;


    /*@} */
    /**@name Private Operators*/
    /*@{ */

    /*@} */
    /**@name Private Operations*/
    /*@{ */


    /*@} */
    /**@name Private  Access */
    /*@{ */


    /*@} */
    /**@name Private Inquiry */
    /*@{ */


    /*@} */
    /**@name Un accessible methods */
    /*@{ */


    /*@} */

}; /* Class ClassName */

/*@} */

/**@name Type Definitions */
/*@{ */


/*@} */

}  /* namespace Kratos.*/

#endif /* KRATOS_NEW_DISPLACEMENT_CRITERIA  defined */

