//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
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
#include "includes/element.h"

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

                    if (i->HasDofFor(ROTATION_X))
                    {
                        array_1d<double,3>& nodal_inertia = i->FastGetSolutionStepValue(NODAL_INERTIA);
                        noalias(nodal_inertia) = ZeroVector(3);
                    }

                 }
             }

         }

         //Calculate and assemble Mass Matrix on nodes

         unsigned int index = 0;

	 bool CalculateLumpedMassMatrix = false;
	 if( rCurrentProcessInfo.Has(COMPUTE_LUMPED_MASS_MATRIX) ){
	   CalculateLumpedMassMatrix = rCurrentProcessInfo[COMPUTE_LUMPED_MASS_MATRIX];	   
	 }
     
         #pragma omp parallel
         {
             int k = OpenMPUtils::ThisThread();
             typename ElementsArrayType::iterator ElemBegin = pElements.begin() + element_partition[k];
             typename ElementsArrayType::iterator ElemEnd = pElements.begin() + element_partition[k + 1];

             for (typename ElementsArrayType::iterator itElem = ElemBegin; itElem != ElemEnd; itElem++)  //MSI: To be parallelized
             {
                 Matrix MassMatrix;

                 Element::GeometryType& geometry = itElem->GetGeometry();
		 
                 (itElem)->CalculateMassMatrix(MassMatrix, rCurrentProcessInfo); 

                 const bool check_rot = itElem->GetGeometry()[0].HasDofFor(ROTATION_X);
                 const unsigned int dimension   = geometry.WorkingSpaceDimension();
                 const int MatSize = MassMatrix.size1();

                 Element::EquationIdVectorType testVector;
                 itElem->EquationIdVector(testVector,rCurrentProcessInfo);
                 KRATOS_WATCH(testVector[0]);

                 index = 0;
                 for (unsigned int i = 0; i <geometry.size(); i++)
                 {
                     index = i*dimension;
                     if (MassMatrix.size1() == (dimension*geometry.size()*2)) index *= 2;

                     double& mass = geometry(i)->FastGetSolutionStepValue(NODAL_MASS);
                     array_1d<double,3>& nodal_inertia = geometry(i)->FastGetSolutionStepValue(NODAL_INERTIA);                     

                    KRATOS_WATCH(geometry(i)->GetDofPosition(DISPLACEMENT_X));
                    KRATOS_WATCH(geometry(i)->GetDofPosition(DISPLACEMENT_Y));
                    KRATOS_WATCH(geometry(i)->GetDofPosition(DISPLACEMENT_Z));
                    KRATOS_WATCH(geometry(i)->GetDofPosition(ROTATION_X));
                    KRATOS_WATCH(geometry(i)->GetDofPosition(ROTATION_Y));
                    KRATOS_WATCH(geometry(i)->GetDofPosition(ROTATION_Z));


                     geometry(i)->SetLock();

                     for (int j = 0;j<MatSize;++j)
                     {
                        // claculate consistent mass matrix
                        // and condense it down to one diag entry
                        mass += MassMatrix(index,j);  
                        if(check_rot & ((index+dimension+2)<= MatSize))
                        {
                        nodal_inertia[0] += MassMatrix(index+dimension,j);   
                        nodal_inertia[1] += MassMatrix(index+dimension+1,j);   
                        nodal_inertia[2] += MassMatrix(index+dimension+2,j); 
                        }
                     }
                     if (check_rot) for (int j =0;j<3;++j) nodal_inertia[j] = std::abs(nodal_inertia[j]);

                     KRATOS_WATCH(mass);
                     KRATOS_WATCH(nodal_inertia);
                     KRATOS_WATCH(MassMatrix);

                     geometry(i)->UnSetLock();
                 }
             }
         }

	    rCurrentProcessInfo[COMPUTE_LUMPED_MASS_MATRIX] = CalculateLumpedMassMatrix;
        
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
        
        // Compute condition contributions to RHS.
        CalculateAndAddConditionsRHS(pScheme, r_model_part);

        // Compute element contributions to RHS.
        CalculateAndAddElementsRHS(pScheme, r_model_part);

        
        KRATOS_CATCH( "" )

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


    #pragma omp parallel for 
    for(int k=0; k<number_of_threads; k++)
    {
       typename ConditionsArrayType::ptr_iterator it_begin=pConditions.ptr_begin()+condition_partition[k];
       typename ConditionsArrayType::ptr_iterator it_end=pConditions.ptr_begin()+condition_partition[k+1];

       for (typename ConditionsArrayType::ptr_iterator it= it_begin; it!=it_end; ++it)
       {

           LocalSystemVectorType RHS_Condition_Contribution = LocalSystemVectorType(0);

           Element::EquationIdVectorType EquationId; //Dummy

           pScheme->Condition_Calculate_RHS_Contribution(*it, RHS_Condition_Contribution, EquationId, rCurrentProcessInfo);


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

        #pragma omp parallel for 
        for(int k=0; k<number_of_threads; k++)
        {
            typename ElementsArrayType::ptr_iterator it_begin=pElements.ptr_begin()+element_partition[k];
            typename ElementsArrayType::ptr_iterator it_end=pElements.ptr_begin()+element_partition[k+1];
            for (typename ElementsArrayType::ptr_iterator it= it_begin; it!=it_end; ++it)
            {

                LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);
                Element::EquationIdVectorType EquationId; //Dummy

                pScheme->Calculate_RHS_Contribution(*it, RHS_Contribution, EquationId, rCurrentProcessInfo);

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
        this->mDofSet.clear(); // = DofsArrayType();

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

