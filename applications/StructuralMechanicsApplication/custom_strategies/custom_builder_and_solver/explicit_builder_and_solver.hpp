//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Klaus B Sautter (based on the work of JMCarbonell)
//					 

#if !defined(KRATOS_EXPLICIT_BUILDER_AND_SOLVER )
#define  KRATOS_EXPLICIT_BUILDER_AND_SOLVER


/* System includes */


#ifdef _OPENMP
#include <omp.h>
#endif

/* External includes */


/* Project includes */
#include "includes/define.h"
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"
#include "includes/model_part.h"

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
         class TDenseSpace, //= DenseSpace<double>,
         class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
         >
class ExplicitBuilderAndSolver
    : public BuilderAndSolver< TSparseSpace, TDenseSpace, TLinearSolver >
{
public:
    /**@name Type Definitions */
    /*@{ */
    KRATOS_CLASS_POINTER_DEFINITION( ExplicitBuilderAndSolver );

    typedef BuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

    typedef typename BaseType::TSchemeType TSchemeType;


    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;


    typedef typename BaseType::NodesArrayType NodesArrayType;

    typedef typename BaseType::ElementsArrayType ElementsArrayType;

    typedef typename BaseType::ConditionsArrayType ConditionsArrayType;


    /*@} */
    /**@name Life Cycle
     */
    /*@{ */

    /** Constructor.
     */
    ExplicitBuilderAndSolver(
        typename TLinearSolver::Pointer pNewLinearSystemSolver)
        : BuilderAndSolver< TSparseSpace, TDenseSpace, TLinearSolver >(pNewLinearSystemSolver)
    {
    }

    /** Destructor.
     */
    virtual ~ExplicitBuilderAndSolver()
    {
    }


    /*@} */
    /**@name Operators
     */
    /*@{ */


    /**
     * This does not build a matrix but is used to create the LHS
     * in explicit time integration each node is updated and thus
     * does not need to assemble any matrix
     */

    void BuildLHS(
            typename TSchemeType::Pointer pScheme,
            ModelPart& rModelPart,
            TSystemMatrixType& A)
    {
        KRATOS_TRY

        //Set Nodal Mass to zero
        NodesArrayType& r_nodes             = rModelPart.Nodes();
        ElementsArrayType& r_elements       = rModelPart.Elements();
        ProcessInfo& r_current_process_info   = rModelPart.GetProcessInfo();


        auto it_node = rModelPart.NodesBegin();
        #pragma omp parallel for 
        for(int i=0;i<static_cast<int>(r_nodes.size());++i)
        {
            (it_node+i)->SetValue(NODAL_MASS,0.00);
        }

        if (rModelPart.NodesBegin()->HasDofFor(ROTATION_Z))
        {
            #pragma omp parallel for 
            for(int i=0;i<static_cast<int>(r_nodes.size());++i)
                {
                    array_1d<double,3>& r_nodal_inertia = (it_node+i)->GetValue(NODAL_INERTIA);
                    noalias(r_nodal_inertia) = ZeroVector(3);
                }
        }

        auto it_elem = rModelPart.ElementsBegin();
        #pragma omp parallel for 
        for (int i=0;i<static_cast<int>(r_elements.size());++i)  
        {
            //Getting nodal mass and inertia from element
            Vector dummy_vector;
            // this function needs to be implemented in the respective 
            // element to provide inertias and nodal masses 
            (it_elem+i)->AddExplicitContribution(dummy_vector,RESIDUAL_VECTOR,NODAL_INERTIA,r_current_process_info);
        }
        
        
        KRATOS_CATCH( "" )
    }

    //**************************************************************************
    //**************************************************************************

    void BuildRHS(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemVectorType& b)
    {
        KRATOS_TRY
        
        // Compute condition contributions to RHS.
        CalculateAndAddConditionsRHS(pScheme, rModelPart);

        // Compute element contributions to RHS.
        CalculateAndAddElementsRHS(pScheme, rModelPart);

        
        KRATOS_CATCH( "" )

    }


    
    
    //***************************************************************************
    //***************************************************************************

   
    void CalculateAndAddConditionsRHS(typename TSchemeType::Pointer pScheme, ModelPart& rModelPart )
    {
        KRATOS_TRY

        ProcessInfo& rCurrentProcessInfo      = rModelPart.GetProcessInfo();
        ConditionsArrayType& pConditions      = rModelPart.Conditions();

        typename ConditionsArrayType::ptr_iterator it_begin=pConditions.ptr_begin();    
        #pragma omp parallel for 
        for (int i =0;i<static_cast<int>(pConditions.size());++i)
        {
            LocalSystemVectorType RHS_Condition_Contribution = LocalSystemVectorType(0);
            Element::EquationIdVectorType equation_id_vector_dummy; //Dummy

            pScheme->Condition_Calculate_RHS_Contribution(*(it_begin+i), RHS_Condition_Contribution,
            equation_id_vector_dummy, rCurrentProcessInfo);
        }







        KRATOS_CATCH("")
    }

    
    //***************************************************************************
    //***************************************************************************


    void CalculateAndAddElementsRHS(typename TSchemeType::Pointer pScheme, ModelPart& rModelPart )
    {

        KRATOS_TRY
        
        ProcessInfo& rCurrentProcessInfo    = rModelPart.GetProcessInfo();
        ElementsArrayType& pElements        = rModelPart.Elements();
        
        typename ElementsArrayType::ptr_iterator it_begin=pElements.ptr_begin();  
        #pragma omp parallel for 
        for (int i =0;i<static_cast<int>(pElements.size());++i)
        {
            LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);
            Element::EquationIdVectorType equation_id_vector_dummy; //Dummy

            pScheme->Calculate_RHS_Contribution(*(it_begin+i), RHS_Contribution,
                equation_id_vector_dummy, rCurrentProcessInfo);
        }


        KRATOS_CATCH("")
    }
    
    
    //**************************************************************************
    //**************************************************************************

    
    void InitializeSolutionStep(
        ModelPart& rModelPart,
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb)
    {
        KRATOS_TRY
	KRATOS_CATCH( "" )
    }

    //**************************************************************************
    //**************************************************************************

    void FinalizeSolutionStep(
        ModelPart& rModelPart,
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb)
    {
        KRATOS_TRY
        KRATOS_CATCH( "" )
    }

    //**************************************************************************
    //**************************************************************************

    void ApplyDirichletConditions(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb)
    {
    }

    //**************************************************************************
    //**************************************************************************

    void ApplyPointLoads(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemVectorType& rb)
    {
    }

    /**
    this function is intended to be called at the end of the solution step to clean up memory
    storage not needed
     */
    void Clear()
    {
        this->mDofSet.clear(); // = DofsArrayType();

        if (this->mpReactionsVector != nullptr)
            TSparseSpace::Clear((this->mpReactionsVector));
        //      this->mReactionsVector = TSystemVectorType();

        this->mpLinearSystemSolver->Clear();

        if (this->GetEchoLevel() > 1)
        {
            std::cout << "ExplicitBuilderAndSolver Clear Function called" << std::endl;
        }
    }

    /**
     * This function is designed to be called once to perform all the checks needed
     * on the input provided. Checks can be "expensive" as the function is designed
     * to catch user's errors.
     * @param r_model_part
     * @return 0 all ok
     */
    virtual int Check(ModelPart& rModelPart)
    {
        KRATOS_TRY

        return 0;
        KRATOS_CATCH( "" )
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

    //******************************************************************************************
    //******************************************************************************************


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

}; /* Class ExplicitBuilderAndSolver */

/*@} */

/**@name Type Definitions */
/*@{ */
/*@} */

} /* namespace Kratos.*/

#endif /* KRATOS_EXPLICIT_BUILDER_AND_SOLVER  defined */

