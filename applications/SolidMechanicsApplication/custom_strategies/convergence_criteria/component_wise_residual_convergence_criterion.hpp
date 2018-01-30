//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_COMPONENT_WISE_RESIDUAL_CONVERGENCE_CRITERION)
#define  KRATOS_COMPONENT_WISE_RESIDUAL_CONVERGENCE_CRITERION


/* System includes */


/* External includes */


/* Project includes */
#include "includes/model_part.h"
#include "includes/define.h"
#include "solid_mechanics_application_variables.h"

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
class ComponentWiseResidualConvergenceCriterion : public ConvergenceCriteria< TSparseSpace, TDenseSpace >
{
public:
    /**@name Type Definitions */
    /*@{ */

    KRATOS_CLASS_POINTER_DEFINITION( ComponentWiseResidualConvergenceCriterion );

    typedef ConvergenceCriteria< TSparseSpace, TDenseSpace > BaseType;

    typedef TSparseSpace SparseSpaceType;

    typedef typename BaseType::TDataType TDataType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

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
    ComponentWiseResidualConvergenceCriterion(
        TDataType NewRatioTolerance,
        TDataType AlwaysConvergedNorm)
        : ConvergenceCriteria< TSparseSpace, TDenseSpace >()
    {
        mRatioTolerance       = NewRatioTolerance;
        mAlwaysConvergedNorm  = AlwaysConvergedNorm;

	//components asked to the elements
	mRHS_Element_Variables.push_back(EXTERNAL_FORCES_VECTOR);
	mRHS_Element_Variables.push_back(INTERNAL_FORCES_VECTOR);
	
	//components asked to the conditions
	mRHS_Condition_Variables.push_back(EXTERNAL_FORCES_VECTOR);
	mRHS_Condition_Variables.push_back(CONTACT_FORCES_VECTOR);
    }

    /** Destructor.
    */
    virtual ~ComponentWiseResidualConvergenceCriterion() {}


    /*@} */
    /**@name Operators
    */
    /*@{ */

    /*Criterion that needs to be called after getting the solution */
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

            //std::cout<<" SIZES (A:"<<A.size1()<<"x"<<A.size2()<<", b:"<<b.size()<<")"<<std::endl;

            //std::cout <<" residual_norm = " << residual_norm << ";  ratio_norm = " << ratio_norm  << std::endl;

	    if (this->GetEchoLevel() == 1)
	      std::cout <<" RESIDUAL FORCES :: ratio = --"<< ratio <<"-- ;  (Expected ratio = " << mRatioTolerance <<", Absolute tol reached = " << residual_norm <<")"<< std::endl;

	    r_model_part.GetProcessInfo()[CONVERGENCE_RATIO] = ratio;
	    r_model_part.GetProcessInfo()[RESIDUAL_NORM] = residual_norm;


            if (ratio <= mRatioTolerance  ||  residual_norm<= mAlwaysConvergedNorm )
            {
	      if (this->GetEchoLevel() == 1)
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
    ) {}

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

    /**
     * Get component wise element components
     */
    std::vector<TSystemVectorType>&  GetRHS_Element_Components()
    { 
      return mRHS_Element_Components;
    } 

    /**
     * Get component wise element variables
     */
    std::vector< Variable< LocalSystemVectorType > >&  GetRHS_Element_Variables()
    { 
      return mRHS_Element_Variables;
    } 

    /**
     * Get component wise condition components
     */
    std::vector<TSystemVectorType>&  GetRHS_Condition_Components()
    { 
      return mRHS_Condition_Components;
    } 

    /**
     * Get component wise condition variables
     */
    std::vector< Variable< LocalSystemVectorType > >&  GetRHS_Condition_Variables()
    { 
      return mRHS_Condition_Variables;
    } 

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

    std::vector<TSystemVectorType> mRHS_Element_Components;

    std::vector< Variable< LocalSystemVectorType > > mRHS_Element_Variables;

    std::vector<TSystemVectorType> mRHS_Condition_Components;

    std::vector< Variable< LocalSystemVectorType > > mRHS_Condition_Variables;
      
    /*@} */
    /**@name Private Operators*/
    /*@{ */

    //Calculate Force Norms component wise

    void CalculateForceNorms(ModelPart& r_model_part,ForceModulus& rforce_norm,TDataType & rresidual_norm,TDataType & rdelta_disp_norm)
    {
        //initialize
        rforce_norm.internal=0;
        rforce_norm.external=0;
        rforce_norm.dynamic =0;
        rforce_norm.reaction=0;

        rforce_norm.residual=0;

        //set norms from component vectors
	rforce_norm.dynamic  = TSparseSpace::TwoNorm(mRHS_Element_Components[0]);   // external forces from elements weight
	//std::cout<<" element external "<<mRHS_Element_Components[0]<<std::endl;

	rforce_norm.internal = TSparseSpace::TwoNorm(mRHS_Element_Components[1]);   // interal forces from elements

	//std::cout<<" element internal "<<mRHS_Element_Components[1]<<std::endl;

	rforce_norm.external = TSparseSpace::TwoNorm(mRHS_Condition_Components[0]); // external forces from conditions
	//std::cout<<" condition external "<<mRHS_Condition_Components[0]<<std::endl;

	rforce_norm.reaction = TSparseSpace::TwoNorm(mRHS_Condition_Components[1]); // contact forces from conditions

	//std::cout<<" condition contact "<<mRHS_Condition_Components[1]<<std::endl;

	rforce_norm.residual = rresidual_norm;

        //std::cout<<" NORMS (internal:"<<rforce_norm.internal<<", external:"<<rforce_norm.external<<", dynamic:"<<rforce_norm.dynamic<<", reaction:"<<rforce_norm.reaction<<", residual:"<<rforce_norm.residual<<", total_residual:"<<rresidual_norm<<", delta_disp_residual:"<<rdelta_disp_norm<<" ) "<<std::endl;


        if(rdelta_disp_norm == rresidual_norm) //is not the norm of the residual,is the norm of the Dx (SuperLU)
	  std::cout<<" something is wrong ResidualNorm and DeltaDisplacementNorm are the same:: check linear solver return "<<std::endl;

        // if(rresidual_norm==0 && rresidual_norm < rforce_norm.residual)
        //   rresidual_norm = rforce_norm.residual;

    }

    //Calculate Force Norms nodally

    void CalculateNodalForceNorms(ModelPart& r_model_part,ForceModulus& rforce_norm,TDataType & rresidual_norm,TDataType & rdelta_disp_norm)
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

            array_1d<double, 3 > & InternalForce  = in->FastGetSolutionStepValue(INTERNAL_FORCE);
            array_1d<double, 3 > & ExternalForce  = in->FastGetSolutionStepValue(EXTERNAL_FORCE);
            array_1d<double, 3 > & DynamicForce   = in->FastGetSolutionStepValue(CONTACT_FORCE);


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

#endif /* KRATOS_COMPONENT_WISE_RESIDUAL_CONVERGENCE_CRITERION  defined */

