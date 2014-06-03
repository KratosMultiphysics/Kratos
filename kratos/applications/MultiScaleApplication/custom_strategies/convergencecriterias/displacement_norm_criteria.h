/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/

/* *********************************************************
*
*   Last Modified by:    $Author: Massimo Petracca $
*   Date:                $Date: 2013-10-03 19:00:00 $
*   Revision:            $Revision: 1.0 $
*
* ***********************************************************/


#if !defined(KRATOS_DISPLACEMENT_NORM_CRITERIA )
#define  KRATOS_DISPLACEMENT_NORM_CRITERIA


/* System includes */
#include <iomanip>

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

/** Residual Norm Criteria.
*/
template<class TSparseSpace,class TDenseSpace>
class DisplacementNormCriteria : public virtual  ConvergenceCriteria< TSparseSpace, TDenseSpace >
{
public:
    /**@name Type Definitions */
    /*@{ */

    KRATOS_CLASS_POINTER_DEFINITION( DisplacementNormCriteria );

    typedef ConvergenceCriteria< TSparseSpace, TDenseSpace > BaseType;

    typedef TSparseSpace SparseSpaceType;

	typedef typename SparseSpaceType::SizeType SizeType;

    typedef typename BaseType::TDataType TDataType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    /*@} */
	
    /**@name Life Cycle
    */
    /*@{ */

    DisplacementNormCriteria(TDataType NewRatioTolerance, TDataType AlwaysConvergedNorm)
        : ConvergenceCriteria< TSparseSpace, TDenseSpace >()
		, mRatioTolerance(NewRatioTolerance)
		, mAlwaysConvergedNorm(AlwaysConvergedNorm)
		, mVerbose(true)
    {
    }

	DisplacementNormCriteria(TDataType NewRatioTolerance, TDataType AlwaysConvergedNorm, bool Verbose)
        : ConvergenceCriteria< TSparseSpace, TDenseSpace >()
		, mRatioTolerance(NewRatioTolerance)
		, mAlwaysConvergedNorm(AlwaysConvergedNorm)
		, mVerbose(Verbose)
    {
    }

    virtual ~DisplacementNormCriteria() 
	{
	}

    /*@} */
	
    /**@name Operators
    */
    /*@{ */

    bool PostCriteria(ModelPart& r_model_part,
                      DofsArrayType& rDofSet,
                      const TSystemMatrixType& A,
                      const TSystemVectorType& Dx,
                      const TSystemVectorType& b)
    {
		SizeType vsize = TSparseSpace::Size(Dx);
		if(vsize == 0) return true;

		TDataType disp_norm       = 0;
        TDataType delta_disp_norm = 0;
        TDataType ratio           = 0;

        CalculateDispNorms(rDofSet, Dx, disp_norm, delta_disp_norm);

		if(disp_norm > 0.0)
			ratio = delta_disp_norm / disp_norm;
		else
			ratio = 1.0;

		TDataType absoluteNorm = delta_disp_norm / (TDataType)vsize;

		if (r_model_part.GetCommunicator().MyPID() == 0 && mVerbose)
		{
			std::stringstream ss;
			ss.precision(3);
			ss << "   DISPLACEMENT NORM | Ratio : " 
			   << std::scientific << ratio << 
			   " | Norm : " 
			   << std::scientific << absoluteNorm 
			   << std::endl;
			std::cout << ss.str();
		}

		return (ratio <= mRatioTolerance || absoluteNorm < mAlwaysConvergedNorm);
    }

    void Initialize(ModelPart& r_model_part)
    {
        BaseType::mConvergenceCriteriaIsInitialized = true;
    }

    void InitializeSolutionStep(
        ModelPart& r_model_part,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& A,
        const TSystemVectorType& Dx,
        const TSystemVectorType& b)
    {
    }

    void FinalizeSolutionStep(
        ModelPart& r_model_part,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& A,
        const TSystemVectorType& Dx,
        const TSystemVectorType& b) 
	{
	}



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
	bool mVerbose;

    /*@} */
    /**@name Private Operators*/
    /*@{ */

    /*@} */
    /**@name Private Operations*/
    /*@{ */

	void CalculateDispNorms(DofsArrayType& rDofSet,const TSystemVectorType& Dx,TDataType& disp_norm,TDataType& delta_disp_norm)
    {
        for(typename DofsArrayType::iterator i_dof = rDofSet.begin() ; i_dof != rDofSet.end() ; ++i_dof)
        {
            if(i_dof->IsFree())
            {
				//here we do: d_n+1^it-d_n
				TDataType dxi = Dx[i_dof->EquationId()];
                delta_disp_norm += dxi * dxi;
				TDataType xi = i_dof->GetSolutionStepValue() - i_dof->GetSolutionStepValue(1);
                disp_norm       += xi * xi;
            }
        }

        delta_disp_norm = sqrt(delta_disp_norm);
        disp_norm = sqrt(disp_norm);
    }

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

