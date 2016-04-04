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

#if !defined(KRATOS_NEW_DISPLACEMENT_CRITERIA )
#define  KRATOS_NEW_DISPLACEMENT_CRITERIA


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
class DisplacementCriteria : virtual public ConvergenceCriteria< TSparseSpace, TDenseSpace >
{
public:
    /**@name Type Definitions */
    /*@{ */

    KRATOS_CLASS_POINTER_DEFINITION( DisplacementCriteria );

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
    DisplacementCriteria(
        TDataType NewRatioTolerance,
        TDataType AlwaysConvergedNorm)
        : ConvergenceCriteria< TSparseSpace, TDenseSpace >()
    {
        mRatioTolerance = NewRatioTolerance;
        mAlwaysConvergedNorm = AlwaysConvergedNorm;
    }

    /** Copy constructor.
    */
    DisplacementCriteria( DisplacementCriteria const& rOther )
      :BaseType(rOther)
      ,mRatioTolerance(rOther.mRatioTolerance)
      ,mAlwaysConvergedNorm(rOther.mAlwaysConvergedNorm)
      ,mReferenceDispNorm(rOther.mReferenceDispNorm)
    {
    }

    /** Destructor.
    */
    virtual ~DisplacementCriteria() {}


    /*@} */
    /**@name Operators
    */
    /*@{ */

    /*Criterias that need to be called before getting the solution */
    virtual bool PreCriteria(
        ModelPart& r_model_part,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& A,
        const TSystemVectorType& Dx,
        const TSystemVectorType& b
    )
    {
        return false;
    }

    /*Criterias that need to be called after getting the solution */
    bool PostCriteria(
        ModelPart& r_model_part,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& A,
        const TSystemVectorType& Dx,
        const TSystemVectorType& b
    )
    {
        if (SparseSpaceType::Size(Dx) != 0) //if we are solving for something
        {
            //TDataType mFinalCorrectionNorm = sqrt(std::inner_product(Dx.begin(),Dx.end(),Dx.begin(),TDataType()));
            //TDataType mFinalCorrectionNorm = sqrt(Dot(Dx,Dx));
            TDataType mFinalCorrectionNorm = TSparseSpace::TwoNorm(Dx);

            TDataType ratio = 0.00;

            CalculateReferenceNorm(rDofSet);

            if(mFinalCorrectionNorm == 0)
                ratio = 0.0;
            else
            {
                if(mReferenceDispNorm == 0)
                    KRATOS_THROW_ERROR(std::logic_error, "NaN norm is detected", "")
                ratio = mFinalCorrectionNorm/mReferenceDispNorm;
            }
KRATOS_WATCH(mFinalCorrectionNorm)
KRATOS_WATCH(mReferenceDispNorm)
            double aaa = SparseSpaceType::Size(Dx);

            double AbsoluteNorm = (mFinalCorrectionNorm/sqrt(aaa));
//KRATOS_WATCH(AbsoluteNorm)
//KRATOS_WATCH(mAlwaysConvergedNorm)
//KRATOS_WATCH(mRatioTolerance)
	    if (this->GetEchoLevel() == 1)
	      std::cout << "DISPLACEMENT CRITERION :: [ Obtained tol = " << ratio << "; Expected ratio = " << mRatioTolerance << "; Absolute tol = " << AbsoluteNorm << "; ]" << std::endl;

	    r_model_part.GetProcessInfo()[CONVERGENCE_RATIO] = ratio;
	    r_model_part.GetProcessInfo()[RESIDUAL_NORM] = AbsoluteNorm;

            if ( ratio <= mRatioTolerance  ||  AbsoluteNorm<mAlwaysConvergedNorm )  //  || (mFinalCorrectionNorm/x.size())<=1e-7)
            {
	      if (this->GetEchoLevel() == 1)
		std::cout << "Convergence is achieved" << std::endl;
	      
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

    TDataType mRatioTolerance;

    TDataType mAlwaysConvergedNorm;

    TDataType mReferenceDispNorm;

    /*@} */
    /**@name Private Operators*/
    /*@{ */

    void CalculateReferenceNorm(DofsArrayType& rDofSet)
    {
        mReferenceDispNorm = TDataType();
        TDataType temp;

        //// x is set to its value at the beginning of the time step
        //typename DofsArrayType::iterator it2;
        //for (it2=rDofSet.begin();it2 != rDofSet.end(); ++it2)
        //{
        //	if ( !(*it2)->IsFixed()  )
        //	{
        //		temp = mpModel->Value(*it2);
        //		mReferenceDispNorm += temp*temp;
        //	}
        //}
        for(typename DofsArrayType::iterator i_dof = rDofSet.begin() ; i_dof != rDofSet.end() ; ++i_dof)
        {
            if(i_dof->IsFree())
            {
                temp = i_dof->GetSolutionStepValue();
                mReferenceDispNorm += temp*temp;
            }
        }
        mReferenceDispNorm = sqrt(mReferenceDispNorm);
    }

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

