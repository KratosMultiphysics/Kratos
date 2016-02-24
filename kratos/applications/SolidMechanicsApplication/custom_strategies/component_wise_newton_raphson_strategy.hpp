//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_COMPONENT_WISE_NEWTON_RAPHSON_STRATEGY )
#define  KRATOS_COMPONENT_WISE_NEWTON_RAPHSON_STRATEGY

/* System includes */


/* External includes */
#include "boost/smart_ptr.hpp"


/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/convergencecriterias/convergence_criteria.h"

//default builder and solver
#include "custom_strategies/custom_builders_and_solvers/component_wise_builder_and_solver.hpp"

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
         class TDenseSpace, // = DenseSpace<double>,
         class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
         >
class ComponentWiseNewtonRaphsonStrategy
    : public ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
public:
    /**@name Type Definitions */
    /*@{ */
    typedef ConvergenceCriteria<TSparseSpace, TDenseSpace> TConvergenceCriteriaType;

    /** Counted pointer of ClassName */
    KRATOS_CLASS_POINTER_DEFINITION( ComponentWiseNewtonRaphsonStrategy );

    typedef SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

    typedef typename BaseType::TBuilderAndSolverType TBuilderAndSolverType;

    typedef typename BaseType::TDataType TDataType;

    typedef TSparseSpace SparseSpaceType;

    typedef typename BaseType::TSchemeType TSchemeType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

    typedef typename BaseType::TSystemMatrixPointerType TSystemMatrixPointerType;
 
    typedef typename BaseType::TSystemVectorPointerType TSystemVectorPointerType;

    typedef typename TBuilderAndSolverType::GlobalSystemComponents GlobalSystemComponentsType;

    /*@} */
    /**@name Life Cycle
     */
    /*@{ */

    /** Constructor.
     */
    ComponentWiseNewtonRaphsonStrategy(
        ModelPart& model_part,
        typename TSchemeType::Pointer pScheme,
        typename TLinearSolver::Pointer pNewLinearSolver,
        typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
        int MaxIterations = 30,
        bool CalculateReactions = false,
        bool ReformDofSetAtEachStep = false,
        bool MoveMeshFlag = false
    )
      : ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(model_part, MoveMeshFlag)
    {
        KRATOS_TRY

	this->mKeepSystemConstantDuringIterations = false;

        //set flags to default values
        this->SetMaxIterationNumber(MaxIterations);
        this->mCalculateReactionsFlag = CalculateReactions;


        this->mReformDofSetAtEachStep = ReformDofSetAtEachStep;

        //saving the convergence criteria to be used
        this->mpConvergenceCriteria = pNewConvergenceCriteria;

        //saving the scheme
        this->mpScheme = pScheme;

        //saving the linear solver
        this->mpLinearSolver = pNewLinearSolver;

        //setting up the default builder and solver
        this->mpBuilderAndSolver = typename TBuilderAndSolverType::Pointer
	  (
	   new ComponentWiseBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver > (this->mpLinearSolver)
	   );

        //set flags to start correcty the calculations
        this->mSolutionStepIsInitialized = false;

        this->mInitializeWasPerformed = false;

	this->mFinalizeSolutionStep = true;

        //tells to the builder and solver if the reactions have to be Calculated or not
        this->GetBuilderAndSolver()->SetCalculateReactionsFlag(this->mCalculateReactionsFlag);

        //tells to the Builder And Solver if the system matrix and vectors need to
        //be reshaped at each step or not
        this->GetBuilderAndSolver()->SetReshapeMatrixFlag(this->mReformDofSetAtEachStep);

        //set EchoLevel to the default value (only time is displayed)
        this->SetEchoLevel(1);

        //by default the matrices are rebuilt at each iteration
        this->SetRebuildLevel(2);


	//component-wise options
	GlobalSystemComponentsType& rGlobalSystem = this->mpBuilderAndSolver->GetGlobalSystemComponents();
	
	rGlobalSystem.SetRHS_Element_Components( this->mpConvergenceCriteria->GetRHS_Element_Components() );
	rGlobalSystem.SetRHS_Element_Variables( this->mpConvergenceCriteria->GetRHS_Element_Variables() );

	rGlobalSystem.SetRHS_Condition_Components( this->mpConvergenceCriteria->GetRHS_Condition_Components() );
	rGlobalSystem.SetRHS_Condition_Variables( this->mpConvergenceCriteria->GetRHS_Condition_Variables() );
	//component-wise options


        KRATOS_CATCH( "" )
    }

    ComponentWiseNewtonRaphsonStrategy(
        ModelPart& model_part,
        typename TSchemeType::Pointer pScheme,
        typename TLinearSolver::Pointer pNewLinearSolver,
        typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
        typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver,
        int MaxIterations = 30,
        bool CalculateReactions = false,
        bool ReformDofSetAtEachStep = false,
        bool MoveMeshFlag = false
    )
      : ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(model_part, pScheme, pNewLinearSolver, pNewConvergenceCriteria, pNewBuilderAndSolver, MaxIterations, CalculateReactions, ReformDofSetAtEachStep, MoveMeshFlag)
    {
        KRATOS_TRY
	  
	// std::cout<<" STRATEGY: ComponentWiseNewtonRaphsonStrategy "<<std::endl;
       
	//component-wise options
	GlobalSystemComponentsType& rGlobalSystem = this->mpBuilderAndSolver->GetGlobalSystemComponents();
	
	rGlobalSystem.SetRHS_Element_Components( this->mpConvergenceCriteria->GetRHS_Element_Components() );
	rGlobalSystem.SetRHS_Element_Variables( this->mpConvergenceCriteria->GetRHS_Element_Variables() );

	rGlobalSystem.SetRHS_Condition_Components( this->mpConvergenceCriteria->GetRHS_Condition_Components() );
	rGlobalSystem.SetRHS_Condition_Variables( this->mpConvergenceCriteria->GetRHS_Condition_Variables() );
	//component-wise options

        KRATOS_CATCH( "" )
    }

    /** Destructor.
     */
    virtual ~ComponentWiseNewtonRaphsonStrategy()
    {
    }

    /*@} */
    /**@name Operators */
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

private:
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
protected:
    /**@name Static Member Variables */
    /*@{ */
    /*@} */
    /**@name Member Variables */
    /*@{ */
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
    ComponentWiseNewtonRaphsonStrategy(const ComponentWiseNewtonRaphsonStrategy& Other)
    {
    };


    /*@} */

}; /* Class ComponentWiseNewtonRaphsonStrategy */

/*@} */

/**@name Type Definitions */
/*@{ */


/*@} */

} /* namespace Kratos.*/

#endif /* KRATOS_COMPONENT_WISE_NEWTON_RAPHSON_STRATEGY  defined */

