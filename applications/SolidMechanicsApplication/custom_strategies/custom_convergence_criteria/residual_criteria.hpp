//
//   Project Name:        KratosSolidMechanicsApplication $
//   Last modified by:    $Author:            JMCarbonell $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_STEP_RESIDUAL_CRITERIA )
#define  KRATOS_STEP_RESIDUAL_CRITERIA


/* System includes */


/* External includes */


/* Project includes */
#include "includes/model_part.h"
#include "includes/define.h"
#include "solid_mechanics_application.h"

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
class ResidualConvergenceCriteria : public virtual  ConvergenceCriteria< TSparseSpace, TDenseSpace >
{
public:
    /**@name Type Definitions */
    /*@{ */

    //typedef boost::shared_ptr< DisplacementCriteria< TSparseSpace, TDenseSpace > > Pointer;
    KRATOS_CLASS_POINTER_DEFINITION( ResidualConvergenceCriteria );

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


protected:
    /**@name Type Definitions */
    /*@{ */
    struct ForceModulus
    {
        TDataType  internal;
        TDataType  external;
        TDataType  reaction;
        TDataType  dynamic;
        TDataType  contact;
        TDataType  residual;
    };

public:

    /** Constructor.
    */
    ResidualConvergenceCriteria(
        TDataType NewRatioTolerance,
        TDataType AlwaysConvergedNorm)
        : ConvergenceCriteria< TSparseSpace, TDenseSpace >()
    {
        mRatioTolerance       = NewRatioTolerance;
        mAlwaysConvergedNorm  = AlwaysConvergedNorm;
    }

    /** Destructor.
    */
    virtual ~ResidualConvergenceCriteria() {}


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

            TDataType residual_norm  = 0;
            residual_norm = TSparseSpace::TwoNorm(b);

            TDataType delta_disp_norm  = 0;
            delta_disp_norm = TSparseSpace::TwoNorm(Dx);

            ForceModulus force_norm;
            CalculateForceNorms(r_model_part,force_norm,residual_norm,delta_disp_norm);

            TDataType ratio_norm=residual_norm;
            CalculateMaximum(force_norm,ratio_norm);

            TDataType ratio=residual_norm;
            if(residual_norm/ratio_norm !=1 && ratio_norm!=0)
                ratio/=ratio_norm;

            std::cout<<" SIZES (A:"<<A.size1()<<"x"<<A.size2()<<", b:"<<b.size()<<")"<<std::endl;

            std::cout <<" residual_norm = " << residual_norm << ";  ratio_norm = " << ratio_norm  << std::endl;

            std::cout <<" RESIDUAL FORCES :: ratio = --"<< ratio <<"-- ;  Expected ratio = " << mRatioTolerance <<" Absolute tol reached = " << residual_norm << std::endl;

            if (ratio <= mRatioTolerance  ||  residual_norm<= mAlwaysConvergedNorm )
            {
                std::cout <<" Convergence is achieved." << std::endl;
                return true;
            }
            else
            {
                if(residual_norm-int(ratio_norm)==0 && int(residual_norm)==1)
                {
                    std::cout << "Convergence is achieved." << std::endl;
                    return true;
                }
                else
                {
                    return false;
                }
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


    void CalculateForceNorms(ModelPart& r_model_part,ForceModulus& rforce_norm,TDataType & rresidual_norm,TDataType & rdelta_disp_norm)
    {
        //initialize
        rforce_norm.internal=0;
        rforce_norm.external=0;
        rforce_norm.dynamic =0;
        rforce_norm.reaction=0;

        rforce_norm.residual=0;

        //set norms from nodal values
        ModelPart::NodesContainerType& rNodes = r_model_part.Nodes();
        for(ModelPart::NodesContainerType::iterator in = rNodes.begin(); in!=rNodes.end(); in++)
        {
	    in->SetLock();

            array_1d<double, 3 > & InternalForce  = in->FastGetSolutionStepValue(FORCE_INTERNAL);
            array_1d<double, 3 > & ExternalForce  = in->FastGetSolutionStepValue(FORCE_EXTERNAL);
            array_1d<double, 3 > & DynamicForce   = in->FastGetSolutionStepValue(FORCE_DYNAMIC);


            array_1d<double, 3 >   Residual       = ExternalForce+InternalForce-DynamicForce;

            for(typename DofsArrayType::iterator i_dof = in->GetDofs().begin(); i_dof!=in->GetDofs().end(); i_dof++)
            {
                if(i_dof->IsFree())
                {
                    if (i_dof->GetVariable() == DISPLACEMENT_X )
                    {
                        rforce_norm.internal += InternalForce[0]*InternalForce[0];
                        rforce_norm.external += ExternalForce[0]*ExternalForce[0];
                        rforce_norm.dynamic  += DynamicForce[0]*DynamicForce[0];

                        rforce_norm.residual += Residual[0]*Residual[0];


                    }
                    if (i_dof->GetVariable() == DISPLACEMENT_Y )
                    {
                        rforce_norm.internal += InternalForce[1]*InternalForce[1];
                        rforce_norm.external += ExternalForce[1]*ExternalForce[1];
                        rforce_norm.dynamic  += DynamicForce[1]*DynamicForce[1];

                        rforce_norm.residual += Residual[1]*Residual[1];

                    }
                    if (i_dof->GetVariable() == DISPLACEMENT_Z )
                    {
                        rforce_norm.internal += InternalForce[2]*InternalForce[2];
                        rforce_norm.external += ExternalForce[2]*ExternalForce[2];
                        rforce_norm.dynamic  += DynamicForce[2]*DynamicForce[2];

                        rforce_norm.residual += Residual[2]*Residual[2];
                    }

                }
                else
                {

                    {
                        if (i_dof->GetVariable() == DISPLACEMENT_X )
                        {
                            rforce_norm.reaction += InternalForce[0]*InternalForce[0];
                        }
                        if (i_dof->GetVariable() == DISPLACEMENT_Y )
                        {
                            rforce_norm.reaction += InternalForce[1]*InternalForce[1];
                        }
                        if (i_dof->GetVariable() == DISPLACEMENT_Z )
                        {
                            rforce_norm.reaction += InternalForce[2]*InternalForce[2];
                        }

                    }
                }
            }

	    in->UnSetLock();
        }


        rforce_norm.internal = sqrt(rforce_norm.internal);
        rforce_norm.external = sqrt(rforce_norm.external);
        rforce_norm.dynamic  = sqrt(rforce_norm.dynamic);
        rforce_norm.reaction = sqrt(rforce_norm.reaction);
        rforce_norm.residual = sqrt(rforce_norm.residual);

        std::cout<<" NORMS (internal:"<<rforce_norm.internal<<", external:"<<rforce_norm.external<<", dynamic:"<<rforce_norm.dynamic<<", reaction:"<<rforce_norm.reaction<<", residual:"<<rforce_norm.residual<<", total_residual:"<<rresidual_norm<<", delta_disp_residual:"<<rdelta_disp_norm<<" ) "<<std::endl;


        if(rdelta_disp_norm==rresidual_norm) //is not the norm of the residual,is the norm of the Dx (SuperLU)
            rresidual_norm = rforce_norm.residual;


        // if(rresidual_norm==0 && rresidual_norm < rforce_norm.residual)
        //   rresidual_norm = rforce_norm.residual;

    }


    void CalculateMaximum(ForceModulus& rforce_norm, TDataType& ratio_norm)
    {
        //ratio_norm is residual norm initially
        TDataType residual_norm = ratio_norm;

        if(ratio_norm<rforce_norm.internal)
            ratio_norm=rforce_norm.internal;

        if(ratio_norm<rforce_norm.external)
            ratio_norm=rforce_norm.external;

        if(ratio_norm<rforce_norm.dynamic)
            ratio_norm=rforce_norm.dynamic;


        bool use_of_reactions=false;
        if(rforce_norm.external == 0 && rforce_norm.dynamic == 0)
            use_of_reactions = true;

        if(ratio_norm == residual_norm )
            use_of_reactions = true;

        //only in this case
        if(use_of_reactions && rforce_norm.external==0)
            ratio_norm=rforce_norm.reaction;

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

#endif /* KRATOS_STEP_RESIDUAL_CRITERIA  defined */

