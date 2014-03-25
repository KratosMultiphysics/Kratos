/* *********************************************************
*
*   Last Modified by:    $Author: rrossi $
*   Date:                $Date: 2007-03-06 10:30:32 $
*   Revision:            $Revision: 1.2 $
*
* ***********************************************************/


#if !defined(KRATOS_RESIDUALBASED_CONVECTION_DIFFUSION_STRATEGY )
#define  KRATOS_RESIDUALBASED_CONVECTION_DIFFUSION_STRATEGY


/* System includes */


/* External includes */
#include "boost/smart_ptr.hpp"


/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "kElectrostatic.h"
#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver.h"



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

/// Short class definition.
/**   Detail class definition.

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
         class TDenseSpace,
         class TLinearSolver
         >
class ResidualBasedConvectionDiffusionStrategy
    : public SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver>
{
public:
    /**@name Type Definitions */
    /*@{ */

    /** Counted pointer of ClassName */
    KRATOS_CLASS_POINTER_DEFINITION( ResidualBasedConvectionDiffusionStrategy );

    typedef SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver> BaseType;

    typedef typename BaseType::TDataType TDataType;

    //typedef typename BaseType::DofSetType DofSetType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;



    /*@} */
    /**@name Life Cycle
    */
    /*@{ */

    /** Constructor.
    */


    ResidualBasedConvectionDiffusionStrategy(
        ModelPart& model_part,
        typename TLinearSolver::Pointer pNewLinearSolver,
        bool ReformDofAtEachIteration = true,
        int time_order = 2,
        int prediction_order = 2
    )
        : SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver>(model_part,false)
    {
        KRATOS_TRY

        mtime_order = time_order;
        mOldDt = 0.00;
        mprediction_order = prediction_order;


        //initializing fractional velocity solution step
        typedef Scheme< TSparseSpace,  TDenseSpace > SchemeType;
        typename SchemeType::Pointer pscheme = typename SchemeType::Pointer
                                               ( new ResidualBasedIncrementalUpdateStaticScheme< TSparseSpace,  TDenseSpace >() );

        bool CalculateReactions = false;
        bool CalculateNormDxFlag = true;

        //choosing the solving strategy
        mstep1 = typename BaseType::Pointer(
                     new ResidualBasedLinearStrategy<TSparseSpace,  TDenseSpace, TLinearSolver >
                     (model_part,pscheme,pNewLinearSolver,CalculateReactions,ReformDofAtEachIteration,CalculateNormDxFlag)  );
        mstep1->SetEchoLevel(2);


        typedef typename BuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver>::Pointer BuilderSolverTypePointer;

        //BuilderSolverTypePointer build = BuilderSolverTypePointer(new	ResidualBasedEliminationBuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver > (pNewLinearSolver) );
        //mstep1 = typename BaseType::Pointer( new ResidualBasedLinearStrategy<TSparseSpace,  TDenseSpace, TLinearSolver > 				(model_part,pscheme,pNewLinearSolver,build,CalculateReactions,ReformDofAtEachIteration,CalculateNormDxFlag)  );
        //mstep1->SetEchoLevel(2);

        KRATOS_CATCH("")
    }



    /** Destructor.
    */
    virtual ~ResidualBasedConvectionDiffusionStrategy() {}

    /** Destructor.
    */

    //*********************************************************************************
    //**********************************************************************
    double Solve()
    {
        KRATOS_TRY

        //calculate the BDF coefficients
        ProcessInfo& rCurrentProcessInfo = BaseType::GetModelPart().GetProcessInfo();
        double Dt = rCurrentProcessInfo[DELTA_TIME];

        if(mOldDt == 0.00) //needed for the first step
            mOldDt = Dt;
        if(mtime_order == 2)
        {
            if(BaseType::GetModelPart().GetBufferSize() < 3)
                KRATOS_ERROR(std::logic_error,"insufficient buffer size for BDF2","")

                rCurrentProcessInfo[BDF_COEFFICIENTS].resize(3);
            Vector& BDFcoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];
            BDFcoeffs[0] =	1.5 / Dt;	//coefficient for step n+1
            BDFcoeffs[1] =	-2.0 / Dt;//coefficient for step n
            BDFcoeffs[2] =	0.5 / Dt;//coefficient for step n-1
        }
        else
        {
            rCurrentProcessInfo[BDF_COEFFICIENTS].resize(2);
            Vector& BDFcoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];
            BDFcoeffs[0] =	1.0 / Dt;	//coefficient for step n+1
            BDFcoeffs[1] =	-1.0 / Dt;//coefficient for step n
        }

        //second order prediction for the velocity
        if(mprediction_order == 2)
        {
            for(ModelPart::NodeIterator i = BaseType::GetModelPart().NodesBegin() ;
                    i != BaseType::GetModelPart().NodesEnd() ; ++i)
            {
                i->FastGetSolutionStepValue(ELECTRIC_POTENTIAL) = 2.00*i->FastGetSolutionStepValue(ELECTRIC_POTENTIAL,1) - i->FastGetSolutionStepValue(ELECTRIC_POTENTIAL,2);
            }
            CalculateProjection();
        }

        //SOLVING THE PROBLEM
        rCurrentProcessInfo[FRACTIONAL_STEP] = 1;

        double Dp_norm = mstep1->Solve();
        CalculateProjection();


//Dp_norm = mstep1->Solve();
//CalculateProjection();
//Dp_norm = mstep1->Solve();
//CalculateProjection();

        return Dp_norm;
        KRATOS_CATCH("")
    }





    //******************************************************************************************************
    //******************************************************************************************************
    //calculation of projection
    void CalculateProjection()
    {
        KRATOS_TRY;

        ProcessInfo& rCurrentProcessInfo = BaseType::GetModelPart().GetProcessInfo();

        //first of all set to zero the nodal variables to be updated nodally
        for(ModelPart::NodeIterator i = BaseType::GetModelPart().NodesBegin() ;
                i != BaseType::GetModelPart().NodesEnd() ; ++i)
        {
            (i)->GetSolutionStepValue(TEMP_CONV_PROJ) = 0.00;
            (i)->GetSolutionStepValue(NODAL_AREA) = 0.00;
        }

        //add the elemental contributions for the calculation of the velocity
        //and the determination of the nodal area
        rCurrentProcessInfo[FRACTIONAL_STEP] = 2;
        for(ModelPart::ElementIterator i = BaseType::GetModelPart().ElementsBegin() ;
                i != BaseType::GetModelPart().ElementsEnd() ; ++i)
        {
            (i)->InitializeSolutionStep(rCurrentProcessInfo);
        }

        //solve nodally for the velocity
        for(ModelPart::NodeIterator i = BaseType::GetModelPart().NodesBegin() ;
                i != BaseType::GetModelPart().NodesEnd() ; ++i)
        {
            double& conv_proj = (i)->GetSolutionStepValue(TEMP_CONV_PROJ);
            double temp = 1.00 / (i)->GetSolutionStepValue(NODAL_AREA);
            conv_proj *= temp;
        }

        KRATOS_CATCH("")
    }

    virtual void SetEchoLevel(int Level)
    {
        mstep1->SetEchoLevel(Level);
    }

    virtual void Clear()
    {
        mstep1->Clear();
    }

    /*@} */
    /**@name Operators
    */
    /*@{ */

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
    typename BaseType::Pointer mstep1;
    double mOldDt;
    int mtime_order;
    int mprediction_order;





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

    /** Copy constructor.
    */
    ResidualBasedConvectionDiffusionStrategy(const ResidualBasedConvectionDiffusionStrategy& Other);


    /*@} */

}; /* Class ResidualBasedConvectionDiffusionStrategy */

/*@} */

/**@name Type Definitions */
/*@{ */


/*@} */

}  /* namespace Kratos.*/

#endif /* KRATOS_RESIDUALBASED_CONVECTION_DIFFUSION_STRATEGY  defined */

