//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_DISPLACEMENT_CONVERGENCE_CRITERION )
#define  KRATOS_DISPLACEMENT_CONVERGENCE_CRITERION


/* System includes */


/* External includes */


/* Project includes */
#include "includes/model_part.h"
#include "includes/define.h"

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


template<class TSparseSpace,
         class TDenseSpace
         >
class DisplacementConvergenceCriterion : public ConvergenceCriteria< TSparseSpace, TDenseSpace >
{
public:
    /**@name Type Definitions */
    /*@{ */

    KRATOS_CLASS_POINTER_DEFINITION( DisplacementConvergenceCriterion );

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
    DisplacementConvergenceCriterion(TDataType NewRatioTolerance,
                                    TDataType AlwaysConvergedNorm)
        : ConvergenceCriteria< TSparseSpace, TDenseSpace >()
    {
        mRatioTolerance = NewRatioTolerance;
        mAlwaysConvergedNorm = AlwaysConvergedNorm;
    }

    /** Copy constructor.	
    */
    DisplacementConvergenceCriterion( DisplacementConvergenceCriterion const& rOther )
      :BaseType(rOther) 
      ,mRatioTolerance(rOther.mRatioTolerance)
      ,mAlwaysConvergedNorm(rOther.mAlwaysConvergedNorm)
    {
    }

    /** Destructor.
    */
    virtual ~DisplacementConvergenceCriterion() {}


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
        if (SparseSpaceType::Size(Dx) != 0) //if we are solving for something
        {
            TDataType disp_norm       = 0;
            TDataType delta_disp_norm = 0;
            TDataType ratio           = 0;

            this->CalculateDispNorms(rDofSet,Dx,disp_norm,delta_disp_norm);


            if( disp_norm!=0 )
                ratio = delta_disp_norm/disp_norm;

	    if( ratio == 0 && int(delta_disp_norm-int(disp_norm))==0 )
	      ratio = 1;

            //std::cout << "delta_disp_norm = " << delta_disp_norm << ";  disp_norm = " << disp_norm << std::endl;

	    TDataType Dx_size = SparseSpaceType::Size(Dx);
            TDataType absolute_norm = (delta_disp_norm/sqrt(Dx_size));

	    if( disp_norm < mAlwaysConvergedNorm && delta_disp_norm <= disp_norm )
	      ratio = absolute_norm;
	    
	    if (r_model_part.GetCommunicator().MyPID() == 0)
            {
                if (this->GetEchoLevel() >= 1)
                {
                    std::cout << "DISPLACEMENT CRITERION :: Ratio = "<< ratio  << ";  Norm = " << absolute_norm << std::endl;
                }
            }
	    
	    r_model_part.GetProcessInfo()[CONVERGENCE_RATIO] = ratio;
	    r_model_part.GetProcessInfo()[RESIDUAL_NORM] = absolute_norm;

            if (ratio <= mRatioTolerance || absolute_norm < mAlwaysConvergedNorm)
            {
	      if (r_model_part.GetCommunicator().MyPID() == 0)
	      {
                    if (this->GetEchoLevel() >= 1)
                    {
                        std::cout << "Convergence is achieved" << std::endl;
                    }
              }
	      return true;
            }
            else
            {
                if( int(delta_disp_norm-int(disp_norm))==0 && int(ratio)==1 && disp_norm <= mRatioTolerance * 1e-2)
		{
		    if (r_model_part.GetCommunicator().MyPID() == 0)
		      {
			if (this->GetEchoLevel() >= 1)
			  {
			    std::cout << "convergence is achieved : - no movement - " << std::endl;
			  }
		      }		  
		    return true;
		}
                else
                {
                  return false;
                }
            }
        }
        else  //nothing is moving
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

    /*@} */
    /**@name Private Operators*/
    /*@{ */


    void CalculateDispNorms(DofsArrayType& rDofSet,const TSystemVectorType& Dx,TDataType& disp_norm,TDataType& delta_disp_norm)
    {
        for(typename DofsArrayType::iterator i_dof = rDofSet.begin() ; i_dof != rDofSet.end() ; ++i_dof)
        {
            if(i_dof->IsFree())
            {
                if (i_dof->GetVariable() == DISPLACEMENT_X || i_dof->GetVariable() == DISPLACEMENT_Y || i_dof->GetVariable() == DISPLACEMENT_Z)
                {
                    //here we do: d_n+1^it-d_n
                    delta_disp_norm += Dx[i_dof->EquationId()] * Dx[i_dof->EquationId()];
                    disp_norm       += pow( i_dof->GetSolutionStepValue() - i_dof->GetSolutionStepValue(1) , 2);
                }
            }
        }

        delta_disp_norm = sqrt(delta_disp_norm);
        disp_norm = sqrt(disp_norm+1e-20);   //to avoid 0
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

#endif /* KRATOS_DISPLACEMENT_CONVERGENCE_CRITERION  defined */

