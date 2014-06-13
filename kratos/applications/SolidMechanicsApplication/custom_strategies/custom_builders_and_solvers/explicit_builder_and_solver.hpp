//
//   Project Name:        KratosSolidMechanicsApplication $
//   Last modified by:    $Author:            JMCarbonell $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_EXPLICIT_BUILDER_AND_SOLVER )
#define  KRATOS_EXPLICIT_BUILDER_AND_SOLVER


/* System includes */
#include <set>

#ifdef _OPENMP
#include <omp.h>
#endif

/* External includes */
#include "boost/smart_ptr.hpp"
#include "utilities/timer.h"

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
    //typedef boost::shared_ptr< ExplicitBuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver> > Pointer;
    KRATOS_CLASS_POINTER_DEFINITION( ExplicitBuilderAndSolver );

    typedef BuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

    typedef typename BaseType::TSchemeType TSchemeType;

    typedef typename BaseType::TDataType TDataType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

    typedef typename BaseType::TSystemMatrixPointerType TSystemMatrixPointerType;

    typedef typename BaseType::TSystemVectorPointerType TSystemVectorPointerType;

    typedef typename BaseType::NodesArrayType NodesArrayType;

    typedef typename BaseType::ElementsArrayType ElementsArrayType;

    typedef typename BaseType::ConditionsArrayType ConditionsArrayType;

    typedef typename BaseType::ElementsContainerType ElementsContainerType;

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


    //**************************************************************************
    //**************************************************************************

    void BuildLHS(
        typename TSchemeType::Pointer pScheme,
        ModelPart& r_model_part,
        TSystemMatrixType& A)
    {
        KRATOS_TRY

         //Set Nodal Mass to zero
         NodesArrayType& pNodes             = r_model_part.Nodes();
         ElementsArrayType& pElements       = r_model_part.Elements();
         ProcessInfo& rCurrentProcessInfo   = r_model_part.GetProcessInfo();

         #ifdef _OPENMP
                 int number_of_threads = omp_get_max_threads();
         #else
                 int number_of_threads = 1;
         #endif

         vector<unsigned int> node_partition;
         OpenMPUtils::CreatePartition(number_of_threads, pNodes.size(), node_partition);

         vector<unsigned int> element_partition;
         OpenMPUtils::CreatePartition(number_of_threads, pElements.size(), element_partition);


         #pragma omp parallel
         {

            #pragma omp for

             for(int k=0; k<number_of_threads; k++)
             {
                 typename NodesArrayType::iterator i_begin=pNodes.ptr_begin()+node_partition[k];
                 typename NodesArrayType::iterator i_end=pNodes.ptr_begin()+node_partition[k+1];

                 for(ModelPart::NodeIterator i=i_begin; i!= i_end; ++i)
                 {
                     double& nodal_mass    =  i->FastGetSolutionStepValue(NODAL_MASS);
                     nodal_mass = 0.0;
                 }
             }

         }

         //Calculate and assemble Mass Matrix on nodes

         unsigned int index = 0;

         #pragma omp parallel
         {
             int k = OpenMPUtils::ThisThread();
             typename ElementsArrayType::iterator ElemBegin = pElements.begin() + element_partition[k];
             typename ElementsArrayType::iterator ElemEnd = pElements.begin() + element_partition[k + 1];

             for (typename ElementsArrayType::iterator itElem = ElemBegin; itElem != ElemEnd; itElem++)  //MSI: LOOP SOBRE ENTER i ABANS...not sure this is well parallelized. pragma omp for doesnt compile.
             {
                 Matrix MassMatrix;

                 Element::GeometryType& geom = itElem->GetGeometry(); // Nodos del elemento

                 (itElem)->CalculateMassMatrix(MassMatrix, rCurrentProcessInfo);

                 const unsigned int dim   = geom.WorkingSpaceDimension();

                 index = 0;
                 for (unsigned int i = 0; i <geom.size(); i++)
                 {
                     double& mass = geom(i)->FastGetSolutionStepValue(NODAL_MASS);

                     geom(i)->SetLock();
                     index = i*dim;

                     mass += MassMatrix(index,index);
                     geom(i)->UnSetLock();
                 }
             }
         }

        KRATOS_CATCH( "" )

    }

    //**************************************************************************
    //**************************************************************************

    void BuildRHS(
        typename TSchemeType::Pointer pScheme,
        ModelPart& r_model_part,
        TSystemVectorType& b)
    {
        KRATOS_TRY
        
          /// Set to zero de RHS
        SetToZeroRHS(pScheme, r_model_part);

        /// Compute the global external nodal force. //Damping included
        CalculateAndAddConditionsRHS(pScheme, r_model_part);

        /// Compute the stress and body force of the element. ( No lineal analysis) //Damping included
        CalculateAndAddElementsRHS(pScheme, r_model_part);

        
        KRATOS_CATCH( "" )

    }


    //**************************************************************************
    //**************************************************************************
    
    void SetToZeroRHS(typename TSchemeType::Pointer pScheme, ModelPart& r_model_part )

    {
        KRATOS_TRY

        NodesArrayType& pNodes   = r_model_part.Nodes();

#ifdef _OPENMP
        int number_of_threads = omp_get_max_threads();
#else
        int number_of_threads = 1;
#endif

        vector<unsigned int> node_partition;
        OpenMPUtils::CreatePartition(number_of_threads, pNodes.size(), node_partition);

        #pragma omp parallel for
        for(int k=0; k<number_of_threads; k++)
        {
            typename NodesArrayType::iterator i_begin=pNodes.ptr_begin()+node_partition[k];
            typename NodesArrayType::iterator i_end=pNodes.ptr_begin()+node_partition[k+1];

            for(ModelPart::NodeIterator i=i_begin; i!= i_end; ++i)
            {
                array_1d<double,3>& node_rhs  = (i->FastGetSolutionStepValue(RHS));
                noalias(node_rhs)             = ZeroVector(3);

            }
        }

        KRATOS_CATCH("")
    }

    
    
//***************************************************************************
//***************************************************************************

   
    void CalculateAndAddConditionsRHS(typename TSchemeType::Pointer pScheme, ModelPart& r_model_part )
    {

    KRATOS_TRY

    ProcessInfo& rCurrentProcessInfo      = r_model_part.GetProcessInfo();
    ConditionsArrayType& pConditions      = r_model_part.Conditions();

#ifdef _OPENMP
    int number_of_threads = omp_get_max_threads();
#else
    int number_of_threads = 1;
#endif

    vector<unsigned int> condition_partition;
    OpenMPUtils::CreatePartition(number_of_threads, pConditions.size(), condition_partition);

    unsigned int index;
    #pragma omp parallel for private (index)
    for(int k=0; k<number_of_threads; k++)
    {
       typename ConditionsArrayType::ptr_iterator it_begin=pConditions.ptr_begin()+condition_partition[k];
       typename ConditionsArrayType::ptr_iterator it_end=pConditions.ptr_begin()+condition_partition[k+1];

       for (typename ConditionsArrayType::ptr_iterator it= it_begin; it!=it_end; ++it)
       {

           Condition::GeometryType& geom = (*it)->GetGeometry();

           const unsigned int& dim = (*it)->GetGeometry().WorkingSpaceDimension();

           LocalSystemVectorType RHS_Condition_Contribution = LocalSystemVectorType(0);
           Element::EquationIdVectorType EquationId; //Dummy

           pScheme->Condition_Calculate_RHS_Contribution(*it, RHS_Condition_Contribution, EquationId, rCurrentProcessInfo);

           for (unsigned int i = 0; i <geom.size(); i++)
           {
               index = i*dim;
               array_1d<double,3>& node_rhs = geom(i)->FastGetSolutionStepValue(RHS);
               for(unsigned int kk=0; kk<dim; kk++)
               {
                   geom(i)->SetLock();
                   node_rhs[kk] += RHS_Condition_Contribution[index+kk];
                   geom(i)->UnSetLock();
               }

           }
       }
    }

    KRATOS_CATCH("")
    }

    
//***************************************************************************
//***************************************************************************


    void CalculateAndAddElementsRHS(typename TSchemeType::Pointer pScheme, ModelPart& r_model_part )
    {

        KRATOS_TRY
        
        ProcessInfo& rCurrentProcessInfo    = r_model_part.GetProcessInfo();
        ElementsArrayType& pElements        = r_model_part.Elements();

#ifdef _OPENMP
        int number_of_threads = omp_get_max_threads();
#else
        int number_of_threads = 1;
#endif

        vector<unsigned int> element_partition;
        OpenMPUtils::CreatePartition(number_of_threads, pElements.size(), element_partition);

        unsigned int index;
        #pragma omp parallel for private (index)
        for(int k=0; k<number_of_threads; k++)
        {
            typename ElementsArrayType::ptr_iterator it_begin=pElements.ptr_begin()+element_partition[k];
            typename ElementsArrayType::ptr_iterator it_end=pElements.ptr_begin()+element_partition[k+1];
            for (typename ElementsArrayType::ptr_iterator it= it_begin; it!=it_end; ++it)
            {

                Element::GeometryType& geom = (*it)->GetGeometry();
                const unsigned int& dim = (*it)->GetGeometry().WorkingSpaceDimension();

                LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);
                Element::EquationIdVectorType EquationId; //Dummy

                pScheme->Calculate_RHS_Contribution(*it, RHS_Contribution, EquationId, rCurrentProcessInfo);

                for (unsigned int i = 0; i <geom.size(); i++)
                {
                    index = i*dim;
                    array_1d<double,3>& node_rhs = geom(i)->FastGetSolutionStepValue(RHS);
                    for(unsigned int kk=0; kk<dim; kk++)
                    {
                        geom(i)->SetLock();


                        node_rhs[kk] += RHS_Contribution[index+kk];


                        geom(i)->UnSetLock();
                    }

                }
            }
        }

        KRATOS_CATCH("")
    }
    
    
    //**************************************************************************
    //**************************************************************************

    
    void InitializeSolutionStep(
        ModelPart& r_model_part,
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
        ModelPart& r_model_part,
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
        ModelPart& r_model_part,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b)
    {
    }

    //**************************************************************************
    //**************************************************************************

    void ApplyPointLoads(
        typename TSchemeType::Pointer pScheme,
        ModelPart& r_model_part,
        TSystemVectorType& b)
    {
    }

    /**
    this function is intended to be called at the end of the solution step to clean up memory
    storage not needed
     */
    void Clear()
    {
        this->mDofSet = DofsArrayType();

        if (this->mpReactionsVector != NULL)
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
    virtual int Check(ModelPart& r_model_part)
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

