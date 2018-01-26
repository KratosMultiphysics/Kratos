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

        #ifdef _OPENMP
            int number_of_threads = omp_get_max_threads();
        #else
            int number_of_threads = 1;
        #endif

        vector<unsigned int> node_partition;
        OpenMPUtils::CreatePartition(number_of_threads, r_nodes.size(), node_partition);

        vector<unsigned int> element_partition;
        OpenMPUtils::CreatePartition(number_of_threads, r_elements.size(), element_partition);


        #pragma omp parallel
        {

            #pragma omp for

                for(int k=0; k<number_of_threads; k++)
                {
                    typename NodesArrayType::iterator i_begin=r_nodes.ptr_begin()+node_partition[k];
                    typename NodesArrayType::iterator i_end=r_nodes.ptr_begin()+node_partition[k+1];

                    for(ModelPart::NodeIterator i=i_begin; i!= i_end; ++i)
                    {
                        double& nodal_mass    =  i->FastGetSolutionStepValue(NODAL_MASS);
                        nodal_mass = 0.0;

                        if (i->HasDofFor(ROTATION_X))
                        {
                            array_1d<double,3>& nodal_inertia = i->FastGetSolutionStepValue(NODAL_INERTIA);
                            noalias(nodal_inertia) = ZeroVector(3);
                        }

                    }
                }
        }
    
        #pragma omp parallel
        {
            int k = OpenMPUtils::ThisThread();
            typename ElementsArrayType::iterator ElemBegin = r_elements.begin() + element_partition[k];
            typename ElementsArrayType::iterator ElemEnd = r_elements.begin() + element_partition[k + 1];

            for (typename ElementsArrayType::iterator itElem = ElemBegin; itElem != ElemEnd; itElem++)  
            {
                //Getting nodal mass and inertia from element
                Vector Testtemp;

                // this function needs to be implemented in the respective 
                // element to provide inertias and nodal masses 
                itElem->AddExplicitContribution(Testtemp,RESIDUAL_VECTOR,NODAL_INERTIA,r_current_process_info);
            }
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

        #ifdef _OPENMP
            int number_of_threads = omp_get_max_threads();
        #else
            int number_of_threads = 1;
        #endif

        vector<unsigned int> condition_partition;
        OpenMPUtils::CreatePartition(number_of_threads, pConditions.size(), condition_partition);


        #pragma omp parallel for 
        for(int k=0; k<number_of_threads; k++)
        {
        typename ConditionsArrayType::ptr_iterator it_begin=pConditions.ptr_begin()+condition_partition[k];
        typename ConditionsArrayType::ptr_iterator it_end=pConditions.ptr_begin()+condition_partition[k+1];

        for (typename ConditionsArrayType::ptr_iterator it= it_begin; it!=it_end; ++it)
        {

            LocalSystemVectorType RHS_Condition_Contribution = LocalSystemVectorType(0);

            Element::EquationIdVectorType equation_id_vector_dummy; //Dummy

            pScheme->Condition_Calculate_RHS_Contribution(*it, RHS_Condition_Contribution, equation_id_vector_dummy, rCurrentProcessInfo);
        }
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

        #ifdef _OPENMP
            int number_of_threads = omp_get_max_threads();
        #else
            int number_of_threads = 1;
        #endif

        vector<unsigned int> element_partition;
        OpenMPUtils::CreatePartition(number_of_threads, pElements.size(), element_partition);

        #pragma omp parallel for 
        for(int k=0; k<number_of_threads; k++)
        {
            typename ElementsArrayType::ptr_iterator it_begin=pElements.ptr_begin()+element_partition[k];
            typename ElementsArrayType::ptr_iterator it_end=pElements.ptr_begin()+element_partition[k+1];
            for (typename ElementsArrayType::ptr_iterator it= it_begin; it!=it_end; ++it)
            {

                LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);
                Element::EquationIdVectorType equation_id_vector_dummy; //Dummy

                pScheme->Calculate_RHS_Contribution(*it, RHS_Contribution, equation_id_vector_dummy, rCurrentProcessInfo);

            }
        }

        KRATOS_CATCH("")
    }
    
    
    //**************************************************************************
    //**************************************************************************

    
    void InitializeSolutionStep(
        ModelPart& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b)
    {
        KRATOS_TRY
	KRATOS_CATCH( "" )
    }

    //**************************************************************************
    //**************************************************************************

    void FinalizeSolutionStep(
        ModelPart& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b)
    {
        KRATOS_TRY
        KRATOS_CATCH( "" )
    }

    //**************************************************************************
    //**************************************************************************

    void ApplyDirichletConditions(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b)
    {
    }

    //**************************************************************************
    //**************************************************************************

    void ApplyPointLoads(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemVectorType& b)
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

