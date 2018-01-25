//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_COMPONENT_WISE_BUILDER_AND_SOLVER )
#define  KRATOS_COMPONENT_WISE_BUILDER_AND_SOLVER


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
class ComponentWiseBuilderAndSolver
    : public BuilderAndSolver< TSparseSpace, TDenseSpace, TLinearSolver >
{
public:

    /**@name Type Definitions */
    /*@{ */

    KRATOS_CLASS_POINTER_DEFINITION( ComponentWiseBuilderAndSolver );

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

    typedef typename BaseType::GlobalSystemComponents GlobalSystemComponentsType;

    typedef typename TSchemeType::LocalSystemComponents LocalSystemComponentsType;

    /*@} */
    /**@name Life Cycle
     */
    /*@{ */

    /** Constructor.
     */
    ComponentWiseBuilderAndSolver(
        typename TLinearSolver::Pointer pNewLinearSystemSolver)
        : BuilderAndSolver< TSparseSpace, TDenseSpace, TLinearSolver >(pNewLinearSystemSolver)
    {
      mGlobalSystem.Initialize();
      std::cout<<" component wise builder and solver selected "<<std::endl;
    }

    /** Destructor.
     */
    virtual ~ComponentWiseBuilderAndSolver()
    {
    }


    /*@} */
    /**@name Operators
     */
    /*@{ */

    GlobalSystemComponentsType& GetGlobalSystemComponents()
    {
      return mGlobalSystem;
    }

    //**************************************************************************
    //**************************************************************************

    void Build(
        typename TSchemeType::Pointer pScheme,
        ModelPart& r_model_part,
        TSystemMatrixType& A,
        TSystemVectorType& b)
    {
        KRATOS_TRY

        if (!pScheme)
            KRATOS_THROW_ERROR( std::runtime_error, "No scheme provided!", "" )

        //getting the elements from the model
        ElementsArrayType& pElements = r_model_part.Elements();

        //std::cout<<" Elements "<<r_model_part.NumberOfElements()<<std::endl;

        //getting the array of the conditions
        ConditionsArrayType& ConditionsArray = r_model_part.Conditions();

        //std::cout<<" Conditions "<<r_model_part.NumberOfConditions()<<std::endl;

        //resetting to zero the vector of reactions
        TSparseSpace::SetToZero(*(BaseType::mpReactionsVector));

        //resetting to false the reactions flag
        bool CalculateReactionsFlag = BaseType::mCalculateReactionsFlag;
	BaseType::mCalculateReactionsFlag = false;

        //double StartTime = GetTickCount();

 	//resize system components and set local variables
	this->ResizeAndInitializeComponents(r_model_part, pScheme);

        // assemble all elements

#ifndef _OPENMP


	//contributions to the element local system
	LocalSystemComponentsType& ElementLocalSystem = pScheme->GetLocalSystemComponents();
	
	LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0, 0);
        LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

       	std::vector<LocalSystemMatrixType> rLHS_LocalElementComponents; 
	ElementLocalSystem.SetLHS_Element_Components(rLHS_LocalElementComponents);

	std::vector<LocalSystemVectorType> rRHS_LocalElementComponents; 
	ElementLocalSystem.SetRHS_Element_Components(rRHS_LocalElementComponents);

        //contributions to the global system
	bool LHS_Element_Components_Set = mGlobalSystem.Are_LHS_Element_Components_Set();
	std::vector<TSystemMatrixType>& rLHS_GlobalElementComponents = mGlobalSystem.GetLHS_Element_Components();

	if( LHS_Element_Components_Set ){
	  
	  if( ElementLocalSystem.GetLHS_Element_Variables().size() != rLHS_LocalElementComponents.size() )
	    rLHS_LocalElementComponents.resize( ElementLocalSystem.GetLHS_Element_Variables().size() );

	  for( unsigned int i=0; i<rLHS_LocalElementComponents.size(); i++ )
	    {
	      rLHS_LocalElementComponents[i] = LocalSystemMatrixType(0, 0);
	    }
	}

	bool RHS_Element_Components_Set = mGlobalSystem.Are_RHS_Element_Components_Set();
	std::vector<TSystemVectorType>& rRHS_GlobalElementComponents = mGlobalSystem.GetRHS_Element_Components();

	if( RHS_Element_Components_Set ){

	  if( ElementLocalSystem.GetRHS_Element_Variables().size() != rRHS_LocalElementComponents.size() )
	    rRHS_LocalElementComponents.resize( ElementLocalSystem.GetRHS_Element_Variables().size() );

	  for( unsigned int i=0; i<rRHS_LocalElementComponents.size(); i++ )
	    {
	      rRHS_LocalElementComponents[i] = LocalSystemVectorType(0);
	    }
	}

        //vector containing the localization in the system of the different terms
	Element::EquationIdVectorType EquationId;

        ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

        for (typename ElementsArrayType::ptr_iterator it = pElements.ptr_begin(); it != pElements.ptr_end(); ++it)
        {
            //calculate elemental contribution
            pScheme->CalculateSystemContributions(*it, LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);

            //assemble the elemental contribution
            this->AssembleLHS(A, LHS_Contribution, EquationId);
            this->AssembleRHS(b, RHS_Contribution, EquationId);

	    if( LHS_Element_Components_Set ){
		  
	      for( unsigned int i=0; i<rLHS_GlobalElementComponents.size(); i++ )
		{	    
		  this->AssembleLHS(rLHS_GlobalElementComponents[i], rLHS_LocalElementComponents[i], EquationId);
		}

	    }

	    if( RHS_Element_Components_Set ){

	      for( unsigned int i=0; i<rRHS_GlobalElementComponents.size(); i++ )
		{
		  this->AssembleRHS(rRHS_GlobalElementComponents[i], rRHS_LocalElementComponents[i], EquationId);
		}

	    }

            // clean local elemental memory
            pScheme->CleanMemory(*it);
        }
        //double EndTime = GetTickCount();

        //std::cout << "total time " << EndTime - StartTime << std::endl;
        //std::cout << "writing in the system matrix " << ccc << std::endl;
        //std::cout << "calculating the elemental contrib " << ddd << std::endl;


        //contributions to the condition local system
	LocalSystemComponentsType& ConditionLocalSystem = pScheme->GetLocalSystemComponents();

        LHS_Contribution.resize(0, 0, false);
        RHS_Contribution.resize(0, false);

	std::vector<LocalSystemMatrixType> rLHS_LocalConditionComponents; 
	ConditionLocalSystem.SetLHS_Condition_Components(rLHS_LocalConditionComponents);

	std::vector<LocalSystemVectorType> rRHS_LocalConditionComponents; 
	ConditionLocalSystem.SetRHS_Condition_Components(rRHS_LocalConditionComponents);

        //contributions to the global system
	bool LHS_Condition_Components_Set = mGlobalSystem.Are_LHS_Condition_Components_Set();
	std::vector<TSystemMatrixType>& rLHS_GlobalConditionComponents = mGlobalSystem.GetLHS_Condition_Components();

	if( LHS_Condition_Components_Set ){

	  if( ConditionLocalSystem.GetLHS_Element_Variables().size() != rLHS_LocalConditionComponents.size() )
	    rLHS_LocalConditionComponents.resize( ConditionLocalSystem.GetLHS_Condition_Variables().size() );

	  for( unsigned int i=0; i<rLHS_LocalConditionComponents.size(); i++ )
	    {
	      rLHS_LocalConditionComponents[i] = LocalSystemMatrixType(0, 0);
	    }

	}

	bool RHS_Condition_Components_Set = mGlobalSystem.Are_RHS_Condition_Components_Set();
	std::vector<TSystemVectorType>& rRHS_GlobalConditionComponents = mGlobalSystem.GetRHS_Condition_Components();

	if( RHS_Condition_Components_Set ){

	  if( ConditionLocalSystem.GetRHS_Element_Variables().size() != rRHS_LocalConditionComponents.size() )
	    rRHS_LocalConditionComponents.resize( ConditionLocalSystem.GetRHS_Condition_Variables().size() );

	  for( unsigned int i=0; i<rRHS_LocalConditionComponents.size(); i++ )
	    {
	      rRHS_LocalConditionComponents[i] = LocalSystemVectorType(0);
	    }

	}
        // assemble all conditions
        for (typename ConditionsArrayType::ptr_iterator it = ConditionsArray.ptr_begin(); it != ConditionsArray.ptr_end(); ++it)
        {
            //calculate condition contribution
  	    pScheme->Condition_CalculateSystemContributions(*it, LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);

            //assemble the condition contribution
            this->AssembleLHS(A, LHS_Contribution, EquationId);
            this->AssembleRHS(b, RHS_Contribution, EquationId);

	    if( LHS_Condition_Components_Set ){

	      for( unsigned int i=0; i<rLHS_GlobalConditionComponents.size(); i++ )
		{	    
		  this->AssembleLHS(rLHS_GlobalConditionComponents[i], rLHS_LocalConditionComponents[i], EquationId);
		}

	    }

	    if( RHS_Condition_Components_Set ){

	      for( unsigned int i=0; i<rRHS_GlobalConditionComponents.size(); i++ )
		{
		  this->AssembleRHS(rRHS_GlobalConditionComponents[i], rRHS_LocalConditionComponents[i], EquationId);
		}

	    }
        }

#else
        //creating an array of lock variables of the size of the system matrix
	std::vector< omp_lock_t > lock_array(b.size());

        int b_size = b.size();
        for (int i = 0; i < b_size; i++)
	  omp_init_lock(&lock_array[i]);

        //create a partition of the element array
        int number_of_threads = omp_get_max_threads();

        vector<unsigned int> element_partition;
        this->CreatePartition(number_of_threads, pElements.size(), element_partition);
        if( this->GetEchoLevel() > 2 && r_model_part.GetCommunicator().MyPID() == 0)
        {
            KRATOS_WATCH( number_of_threads )
            KRATOS_WATCH( element_partition )
        }

	
        double start_prod = omp_get_wtime();

        #pragma omp parallel for
        for (int k = 0; k < number_of_threads; k++)
        {

	    //Parallelism problems, accesing the same member variables of the scheme
	    typename TSchemeType::Pointer pLocalScheme = pScheme->Clone();

	    //contributions to the element local system
	    LocalSystemComponentsType& ElementLocalSystem = pLocalScheme->GetLocalSystemComponents();
	  
            LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0, 0);
            LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);


	    std::vector<LocalSystemMatrixType> rLHS_LocalElementComponents; 
	    ElementLocalSystem.SetLHS_Element_Components(rLHS_LocalElementComponents);
	    
	    std::vector<LocalSystemVectorType> rRHS_LocalElementComponents; 
	    ElementLocalSystem.SetRHS_Element_Components(rRHS_LocalElementComponents);
	    
    
	    //contributions to the global system
	    bool LHS_Element_Components_Set = mGlobalSystem.Are_LHS_Element_Components_Set();
	    std::vector<TSystemMatrixType>& rLHS_GlobalElementComponents = mGlobalSystem.GetLHS_Element_Components();

	    if( LHS_Element_Components_Set ){
	    	    
	      if( ElementLocalSystem.GetLHS_Element_Variables().size() != rLHS_LocalElementComponents.size() )
		rLHS_LocalElementComponents.resize( ElementLocalSystem.GetLHS_Element_Variables().size() );

	      for( unsigned int i=0; i<rLHS_LocalElementComponents.size(); i++ )
		{
		  rLHS_LocalElementComponents[i] = LocalSystemMatrixType(0, 0);
		}

	    }

	    bool RHS_Element_Components_Set = mGlobalSystem.Are_RHS_Element_Components_Set();	
	    std::vector<TSystemVectorType>& rRHS_GlobalElementComponents = mGlobalSystem.GetRHS_Element_Components();

	    if( RHS_Element_Components_Set ){	    
	    
	      if( ElementLocalSystem.GetRHS_Element_Variables().size() != rRHS_LocalElementComponents.size() )
		rRHS_LocalElementComponents.resize( ElementLocalSystem.GetRHS_Element_Variables().size() );

	      for( unsigned int i=0; i<rRHS_LocalElementComponents.size(); i++ )
		{
		  rRHS_LocalElementComponents[i] = LocalSystemVectorType(0);
		}
	    
	    }

	    //vector containing the localization in the system of the different terms
	    Element::EquationIdVectorType EquationId;

            ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

            typename ElementsArrayType::ptr_iterator it_begin = pElements.ptr_begin() + element_partition[k];
            typename ElementsArrayType::ptr_iterator it_end = pElements.ptr_begin() + element_partition[k + 1];

            // assemble all elements
            for (typename ElementsArrayType::ptr_iterator it = it_begin; it != it_end; ++it)
            {

                //calculate elemental contribution
	        pLocalScheme->CalculateSystemContributions(*it, LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);

                //assemble the elemental contribution
                this->Assemble(A, b, LHS_Contribution, RHS_Contribution, EquationId, lock_array);

		//assemble the elemental contribution
		//AssembleLHS(mGlobalSystem.LeftHandSideContribution, ElementLocalSystem.LeftHandSideContribution, EquationId);
		//AssembleRHS(mGlobalSystem.RightHandSideContribution, ElementLocalSystem.RightHandSideContribution, EquationId);
		if( LHS_Element_Components_Set ){

		  for( unsigned int i=0; i<rLHS_GlobalElementComponents.size(); i++ )
		    {	    
		      this->AssembleLHS(rLHS_GlobalElementComponents[i], rLHS_LocalElementComponents[i], EquationId, lock_array);
		    }
		}

		if( RHS_Element_Components_Set ){

		  for( unsigned int i=0; i<rRHS_GlobalElementComponents.size(); i++ )
		    {
		      this->AssembleRHS(rRHS_GlobalElementComponents[i], rRHS_LocalElementComponents[i], EquationId, lock_array);
		    }
		}

                // clean local elemental memory
                pLocalScheme->CleanMemory(*it);
		
            }
	    
        }


        vector<unsigned int> condition_partition;
        this->CreatePartition(number_of_threads, ConditionsArray.size(), condition_partition);

        #pragma omp parallel for
        for (int k = 0; k < number_of_threads; k++)
        {
	   	  
	    //Parallelism problems, accesing the same member variables of the scheme
  	    typename TSchemeType::Pointer pLocalScheme = pScheme->Clone();
	  
	    //contributions to the condition local system
	    LocalSystemComponentsType& ConditionLocalSystem = pLocalScheme->GetLocalSystemComponents();
	    
	    LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0, 0);
            LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

	    std::vector<LocalSystemMatrixType> rLHS_LocalConditionComponents; 
	    ConditionLocalSystem.SetLHS_Condition_Components(rLHS_LocalConditionComponents);
	    
	    std::vector<LocalSystemVectorType> rRHS_LocalConditionComponents; 
	    ConditionLocalSystem.SetRHS_Condition_Components(rRHS_LocalConditionComponents);
	    
 	    //contributions to the global system
	    bool LHS_Condition_Components_Set = mGlobalSystem.Are_LHS_Condition_Components_Set();
	    std::vector<TSystemMatrixType>& rLHS_GlobalConditionComponents = mGlobalSystem.GetLHS_Condition_Components();

	    if( LHS_Condition_Components_Set ){
	      
	      if( ConditionLocalSystem.GetLHS_Condition_Variables().size() != rLHS_LocalConditionComponents.size() )
		rLHS_LocalConditionComponents.resize( ConditionLocalSystem.GetLHS_Condition_Variables().size() );

	      for( unsigned int i=0; i<rLHS_LocalConditionComponents.size(); i++ )
		{
		  rLHS_LocalConditionComponents[i] = LocalSystemMatrixType(0, 0);
		}
	    
	    }

	    
	    bool RHS_Condition_Components_Set = mGlobalSystem.Are_RHS_Condition_Components_Set();
	    std::vector<TSystemVectorType>& rRHS_GlobalConditionComponents = mGlobalSystem.GetRHS_Condition_Components();

	    if( RHS_Condition_Components_Set ){
    
	      if( ConditionLocalSystem.GetRHS_Condition_Variables().size() != rRHS_LocalConditionComponents.size() )
		rRHS_LocalConditionComponents.resize( ConditionLocalSystem.GetRHS_Condition_Variables().size() );

	      for( unsigned int i=0; i<rRHS_LocalConditionComponents.size(); i++ )
		{
		  rRHS_LocalConditionComponents[i] = LocalSystemVectorType(0);
		}

	    }
	    //vector containing the localization in the system of the different terms
	    Element::EquationIdVectorType EquationId;

            ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

            typename ConditionsArrayType::ptr_iterator it_begin = ConditionsArray.ptr_begin() + condition_partition[k];
            typename ConditionsArrayType::ptr_iterator it_end = ConditionsArray.ptr_begin() + condition_partition[k + 1];

            // assemble all elements
            for (typename ConditionsArrayType::ptr_iterator it = it_begin; it != it_end; ++it)
            {
                //calculate condition contribution
	        pLocalScheme->Condition_CalculateSystemContributions(*it, LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);

                //assemble the condition contribution
                this->Assemble(A, b, LHS_Contribution, RHS_Contribution, EquationId, lock_array);

 		//assemble the elemental contribution
		//AssembleLHS(mGlobalSystem.LeftHandSideContribution, ConditionLocalSystem.LeftHandSideContribution, EquationId);
		//AssembleRHS(mGlobalSystem.RightHandSideContribution, ConditionLocalSystem.RightHandSideContribution, EquationId);
		
		if( LHS_Condition_Components_Set ){

		  for( unsigned int i=0; i<rLHS_GlobalConditionComponents.size(); i++ )
		    {	    
		      this->AssembleLHS(rLHS_GlobalConditionComponents[i], rLHS_LocalConditionComponents[i], EquationId, lock_array);
		    }
		}
		
		if( RHS_Condition_Components_Set ){
	      
		  for( unsigned int i=0; i<rRHS_GlobalConditionComponents.size(); i++ )
		    {
		      this->AssembleRHS(rRHS_GlobalConditionComponents[i], rRHS_LocalConditionComponents[i], EquationId, lock_array);
		    }
		}
            }

        }

	double stop_prod = omp_get_wtime();
        if (this->GetEchoLevel() > 2 && r_model_part.GetCommunicator().MyPID() == 0)
            std::cout << "time: " << stop_prod - start_prod << std::endl;

        for (int i = 0; i < b_size; i++)
            omp_destroy_lock(&lock_array[i]);
        if( this->GetEchoLevel() > 2 && r_model_part.GetCommunicator().MyPID() == 0)
        {
            KRATOS_WATCH( "finished parallel building" )
        }

        //to ensure that all the threads are syncronized here
        // #pragma omp barrier
#endif

	//recovering the reactions flag
	BaseType::mCalculateReactionsFlag = CalculateReactionsFlag;


        KRATOS_CATCH( "" )

    }

    //**************************************************************************
    //**************************************************************************

    void BuildLHS(
        typename TSchemeType::Pointer pScheme,
        ModelPart& r_model_part,
	TSystemMatrixType& A)
    {
        KRATOS_TRY

        //getting the elements from the model
        ElementsArrayType& pElements = r_model_part.Elements();

        //getting the array of the conditions
        ConditionsArrayType& ConditionsArray = r_model_part.Conditions();

        //resetting to zero the vector of reactions
        TSparseSpace::SetToZero(*(BaseType::mpReactionsVector));

 	//resize system components and set local variables
	this->ResizeAndInitializeComponents(r_model_part, pScheme);

        //contributions to the element local system
	LocalSystemComponentsType& ElementLocalLHS = pScheme->GetLocalSystemComponents();
	
        LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0, 0);

	std::vector<LocalSystemMatrixType> rLHS_LocalElementComponents; 
	ElementLocalLHS.SetLHS_Element_Components(rLHS_LocalElementComponents);

        //contributions to the global system
	bool LHS_Element_Components_Set = mGlobalSystem.Are_LHS_Element_Components_Set();
	std::vector<TSystemMatrixType>& rLHS_GlobalElementComponents = mGlobalSystem.GetLHS_Element_Components();
	  
	if( LHS_Element_Components_Set ){
	      
	  if( mGlobalSystem.GetLHS_Element_Variables().size() != rLHS_LocalElementComponents.size() )
	    rLHS_LocalElementComponents.resize( mGlobalSystem.GetLHS_Element_Variables().size() );

	  for( unsigned int i=0; i<rLHS_LocalElementComponents.size(); i++ )
	    {
	      rLHS_LocalElementComponents[i] = LocalSystemMatrixType(0, 0);
	    }

	}
        //vector containing the localization in the system of the different
        //terms
        Element::EquationIdVectorType EquationId;

        ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

        // assemble all elements
        for (typename ElementsArrayType::ptr_iterator it = pElements.ptr_begin(); it != pElements.ptr_end(); ++it)
        {
            //calculate elemental contribution
            pScheme->Calculate_LHS_Contribution(*it, LHS_Contribution, EquationId, CurrentProcessInfo);

            //assemble the elemental contribution
            this->AssembleLHS(A, LHS_Contribution, EquationId);

	    if( LHS_Element_Components_Set ){

	      for( unsigned int i=0; i<rLHS_GlobalElementComponents.size(); i++ )
		{	    
		  this->AssembleLHS(rLHS_GlobalElementComponents[i], rLHS_LocalElementComponents[i], EquationId);
		}

	    }

            // clean local elemental memory
            pScheme->CleanMemory(*it);
        }


        //contributions to the condition local system
	LocalSystemComponentsType& ConditionLocalLHS = pScheme->GetLocalSystemComponents();

        LHS_Contribution.resize(0, 0, false);

	std::vector<LocalSystemMatrixType> rLHS_LocalConditionComponents; 
	ConditionLocalLHS.SetLHS_Condition_Components(rLHS_LocalConditionComponents);

        //contributions to the global system
	bool LHS_Condition_Components_Set = mGlobalSystem.Are_LHS_Condition_Components_Set();
	std::vector<TSystemMatrixType>& rLHS_GlobalConditionComponents = mGlobalSystem.GetLHS_Condition_Components();

	if( LHS_Condition_Components_Set ){

	  if( mGlobalSystem.GetLHS_Condition_Variables().size() != rLHS_LocalConditionComponents.size() )
	    rLHS_LocalConditionComponents.resize( mGlobalSystem.GetLHS_Condition_Variables().size() );

	  for( unsigned int i=0; i<rLHS_LocalConditionComponents.size(); i++ )
	    {
	      rLHS_LocalConditionComponents[i] = LocalSystemMatrixType(0, 0);
	    }

	}
        // assemble all conditions
        for (typename ConditionsArrayType::ptr_iterator it = ConditionsArray.ptr_begin(); it != ConditionsArray.ptr_end(); ++it)
        {
            //calculate elemental contribution
            pScheme->Condition_Calculate_LHS_Contribution(*it, LHS_Contribution, EquationId, CurrentProcessInfo);

            //assemble the elemental contribution
            this->AssembleLHS(A, LHS_Contribution, EquationId);

	    if( LHS_Condition_Components_Set ){

	      for( unsigned int i=0; i<rLHS_GlobalConditionComponents.size(); i++ )
		{	    
		  this->AssembleLHS(rLHS_GlobalConditionComponents[i], rLHS_LocalConditionComponents[i], EquationId);
		}

	    }

        }

        KRATOS_CATCH( "" )

    }

    //**************************************************************************
    //**************************************************************************

    void BuildLHS_CompleteOnFreeRows(
        typename TSchemeType::Pointer pScheme,
        ModelPart& r_model_part,
        TSystemMatrixType& A)
    {
        KRATOS_TRY

        //getting the elements from the model
        ElementsArrayType& pElements = r_model_part.Elements();

        //getting the array of the conditions
        ConditionsArrayType& ConditionsArray = r_model_part.Conditions();

        ProcessInfo& rCurrentProcessInfo = r_model_part.GetProcessInfo();

        //resetting to zero the vector of reactions
        TSparseSpace::SetToZero(*(BaseType::mpReactionsVector));

        //contributions to the system
        LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0, 0);

        //vector containing the localization in the system of the different
        //terms
        Element::EquationIdVectorType EquationId;

        // assemble all elements
        for (typename ElementsArrayType::ptr_iterator it = pElements.ptr_begin(); it != pElements.ptr_end(); ++it)
        {
            //calculate elemental contribution
            pScheme->Calculate_LHS_Contribution(*it, LHS_Contribution, EquationId, rCurrentProcessInfo);

            //assemble the elemental contribution
            AssembleLHS_CompleteOnFreeRows(A, LHS_Contribution, EquationId);

            // clean local elemental memory
            pScheme->CleanMemory(*it);
        }

        LHS_Contribution.resize(0, 0, false);
        // assemble all conditions
        for (typename ConditionsArrayType::ptr_iterator it = ConditionsArray.ptr_begin(); it != ConditionsArray.ptr_end(); ++it)
        {
            //calculate elemental contribution
            pScheme->Condition_Calculate_LHS_Contribution(*it, LHS_Contribution, EquationId, rCurrentProcessInfo);

            //assemble the elemental contribution
            AssembleLHS_CompleteOnFreeRows(A, LHS_Contribution, EquationId);
        }


        KRATOS_CATCH( "" )

    }

    //**************************************************************************
    //**************************************************************************

    void SystemSolve(
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b
    )
    {
        KRATOS_TRY

        double norm_b;
        if (TSparseSpace::Size(b) != 0)
            norm_b = TSparseSpace::TwoNorm(b);
        else
            norm_b = 0.00;

        if (norm_b != 0.00)
        {
            //do solve
            BaseType::mpLinearSystemSolver->Solve(A, Dx, b);
        }
        else
            TSparseSpace::SetToZero(Dx);

        //prints informations about the current time
        if (this->GetEchoLevel() > 1)
        {
            std::cout << *(BaseType::mpLinearSystemSolver) << std::endl;
        }

        KRATOS_CATCH( "" )

    }

    void SystemSolveWithPhysics(
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b,
        ModelPart& r_model_part
    )
    {
        KRATOS_TRY

        double norm_b;
        if (TSparseSpace::Size(b) != 0)
            norm_b = TSparseSpace::TwoNorm(b);
        else
            norm_b = 0.00;

        if (norm_b != 0.00)
        {
            //provide physical data as needed
            if(BaseType::mpLinearSystemSolver->AdditionalPhysicalDataIsNeeded() )
                BaseType::mpLinearSystemSolver->ProvideAdditionalData(A, Dx, b, BaseType::mDofSet, r_model_part);

            //do solve
            BaseType::mpLinearSystemSolver->Solve(A, Dx, b);
        }
        else
        {
            TSparseSpace::SetToZero(Dx);
            std::cout << "  -[WARNING: system RHS is zero]-" << std::endl;
        }

        //prints informations about the current time
        if (this->GetEchoLevel() > 1)
        {
            std::cout << *(BaseType::mpLinearSystemSolver) << std::endl;
        }

        KRATOS_CATCH( "" )

    }

    //**************************************************************************
    //**************************************************************************

    void BuildAndSolve(
        typename TSchemeType::Pointer pScheme,
        ModelPart& r_model_part,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b)
    {
        KRATOS_TRY

        Timer::Start("Build");

        //resetting to false the reactions flag
        bool CalculateReactionsFlag = BaseType::mCalculateReactionsFlag;
	BaseType::mCalculateReactionsFlag = false;


        this->Build(pScheme, r_model_part, A, b);

	//recovering the reactions flag
	BaseType::mCalculateReactionsFlag = CalculateReactionsFlag;


        Timer::Stop("Build");


        //does nothing...dirichlet conditions are naturally dealt with in defining the residual
        ApplyDirichletConditions(pScheme, r_model_part, A, Dx, b);

        if (this->GetEchoLevel() >= 2)
        {
            std::cout << "Before solving the system :" << std::endl;
	    if (this->GetEchoLevel() == 3){
	      std::cout << "System Matrix   = " << A << std::endl;
	      std::cout << "Unknowns vector = " << Dx << std::endl;
	    }
            std::cout << "RHS vector      = " << b << std::endl;
        }

        Timer::Start("Solve");

        SystemSolveWithPhysics(A, Dx, b, r_model_part);

        Timer::Stop("Solve");


        if (this->GetEchoLevel() >= 2)
        {
            std::cout << "After solving the system:" << std::endl;
	    if (this->GetEchoLevel() == 3){
	      std::cout << "System Matrix   = " << A << std::endl;
	    }
	    std::cout << "Results vector  = " << Dx << std::endl;
            std::cout << "RHS vector      = " << b << std::endl;
        }

        KRATOS_CATCH( "" )
    }


    //**************************************************************************
    //**************************************************************************

    void BuildRHSAndSolve(
        typename TSchemeType::Pointer pScheme,
        ModelPart& r_model_part,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b)
    {
        KRATOS_TRY
	  
	//resetting to false the reactions flag
        bool CalculateReactionsFlag = BaseType::mCalculateReactionsFlag;
	BaseType::mCalculateReactionsFlag = false;

	this->BuildRHS(pScheme, r_model_part, b);

	//recovering the reactions flag
	BaseType::mCalculateReactionsFlag = CalculateReactionsFlag;
		

        SystemSolve(A, Dx, b);

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

        //Getting the Elements
        ElementsArrayType& pElements = r_model_part.Elements();

        //getting the array of the conditions
        ConditionsArrayType& ConditionsArray = r_model_part.Conditions();

        ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

        //resetting to zero the vector of reactions
        TSparseSpace::SetToZero(*(BaseType::mpReactionsVector));

 	//resize system components and set local variables
	this->ResizeAndInitializeComponents(r_model_part, pScheme);

        //contributions to the system
	LocalSystemComponentsType& ElementLocalRHS = pScheme->GetLocalSystemComponents();
	
        LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);
	
	std::vector<LocalSystemVectorType> rRHS_LocalElementComponents; 
	ElementLocalRHS.SetRHS_Element_Components(rRHS_LocalElementComponents);
	
	//contributions to the global system
	bool RHS_Element_Components_Set = mGlobalSystem.Are_RHS_Element_Components_Set();
	std::vector<TSystemVectorType>& rRHS_GlobalElementComponents = mGlobalSystem.GetRHS_Element_Components();
	  
	if( RHS_Element_Components_Set ){

	  if( mGlobalSystem.GetRHS_Element_Variables().size() != rRHS_LocalElementComponents.size() )
	    rRHS_LocalElementComponents.resize( mGlobalSystem.GetRHS_Element_Variables().size() );

	  for( unsigned int i=0; i<rRHS_LocalElementComponents.size(); i++ )
	    {
	      rRHS_LocalElementComponents[i] = LocalSystemVectorType(0);
	    }
	    
	}
        //vector containing the localization in the system of the different
        //terms
        Element::EquationIdVectorType EquationId;

        // assemble all elements
        for (typename ElementsArrayType::ptr_iterator it = pElements.ptr_begin(); it != pElements.ptr_end(); ++it)
        {
            //calculate elemental Right Hand Side Contribution
            pScheme->Calculate_RHS_Contribution(*it, RHS_Contribution, EquationId, CurrentProcessInfo);

            //assemble the elemental contribution
            this->AssembleRHS(b, RHS_Contribution, EquationId);

	    if( RHS_Element_Components_Set ){
	  
	      for( unsigned int i=0; i<rRHS_GlobalElementComponents.size(); i++ )
		{
		  this->AssembleRHS(rRHS_GlobalElementComponents[i], rRHS_LocalElementComponents[i], EquationId);
		}

	    }
        }

        //contributions to the condition local system
	LocalSystemComponentsType& ConditionLocalRHS = pScheme->GetLocalSystemComponents();

        RHS_Contribution.resize(0, false);

	std::vector<LocalSystemVectorType> rRHS_LocalConditionComponents; 
	ConditionLocalRHS.SetRHS_Condition_Components(rRHS_LocalConditionComponents);

        //contributions to the global system
	bool RHS_Condition_Components_Set = mGlobalSystem.Are_RHS_Condition_Components_Set();
	std::vector<TSystemVectorType>& rRHS_GlobalConditionComponents = mGlobalSystem.GetRHS_Condition_Components();
	  
	if( RHS_Condition_Components_Set ){

	  if( mGlobalSystem.GetRHS_Condition_Variables().size() != rRHS_LocalConditionComponents.size() )
	    rRHS_LocalConditionComponents.resize( mGlobalSystem.GetRHS_Condition_Variables().size() );

	  for( unsigned int i=0; i<rRHS_LocalConditionComponents.size(); i++ )
	    {
	      rRHS_LocalConditionComponents[i] = LocalSystemVectorType(0);
	    }

	}

        // assemble all conditions
        for (typename ConditionsArrayType::ptr_iterator it = ConditionsArray.ptr_begin(); it != ConditionsArray.ptr_end(); ++it)
        {
            //calculate elemental contribution
            pScheme->Condition_Calculate_RHS_Contribution(*it, RHS_Contribution, EquationId, CurrentProcessInfo);

            //assemble the elemental contribution
            this->AssembleRHS(b, RHS_Contribution, EquationId);

	    if( RHS_Condition_Components_Set ){

	      for( unsigned int i=0; i<rRHS_GlobalConditionComponents.size(); i++ )
		{
		  this->AssembleRHS(rRHS_GlobalConditionComponents[i], rRHS_LocalConditionComponents[i], EquationId);
		}

	    }

        }

        KRATOS_CATCH( "" )

    }

    //**************************************************************************
    //**************************************************************************

    void SetUpDofSet(
        typename TSchemeType::Pointer pScheme,
        ModelPart& r_model_part
    )
    {
        KRATOS_TRY;

        if( this->GetEchoLevel() > 1 && r_model_part.GetCommunicator().MyPID() == 0)
        {
            std::cout << " Setting up the dofs " << std::endl;
        }

        //Gets the array of elements from the modeler
        ElementsArrayType& pElements = r_model_part.Elements();

        Element::DofsVectorType ElementalDofList;

        ProcessInfo& rCurrentProcessInfo = r_model_part.GetProcessInfo();

        DofsArrayType Doftemp;
        BaseType::mDofSet.clear(); // = DofsArrayType();
        //mDofSet.clear();

        //double StartTime = GetTickCount();
        for (typename ElementsArrayType::ptr_iterator it = pElements.ptr_begin(); it != pElements.ptr_end(); ++it)
        {
            // gets list of Dof involved on every element
            //aaa = GetTickCount();
            pScheme->GetElementalDofList(*it, ElementalDofList, rCurrentProcessInfo);
            //bbb += GetTickCount() - aaa;
            // KRATOS_WATCH( (*it)->Id() )
            // std::cout << "node ids" << std::endl;
            // for(unsigned int i=0; i<((*it)->GetGeometry()).size(); i++)
            //         std::cout << ((*it)->GetGeometry())[i].Id() << " ";
            // std::cout << std::endl;
            // for(unsigned int i=0; i<ElementalDofList.size(); i++)
            //         std::cout << (ElementalDofList[i]->Id()) << " ";
            // std::cout << std::endl;

            //ccc = GetTickCount();
            for (typename Element::DofsVectorType::iterator i = ElementalDofList.begin(); i != ElementalDofList.end(); ++i)
            {
                Doftemp.push_back(i->get());
                //mDofSet.push_back(*i);
            }
            //ddd += GetTickCount() - ccc;
        }

        //std::cout << "searching " << bbb << std::endl;
        //std::cout << "inserting " << ddd << std::endl;

        //taking in account conditions
        ConditionsArrayType& pConditions = r_model_part.Conditions();
        for (typename ConditionsArrayType::ptr_iterator it = pConditions.ptr_begin(); it != pConditions.ptr_end(); ++it)
        {
            // gets list of Dof involved on every element
            pScheme->GetConditionDofList(*it, ElementalDofList, rCurrentProcessInfo);

            //ccc = GetTickCount();
            for (typename Element::DofsVectorType::iterator i = ElementalDofList.begin(); i != ElementalDofList.end(); ++i)
            {
                //mDofSet.push_back(*i);
                Doftemp.push_back(i->get());
            }
            //ddd += GetTickCount() - ccc;
        }
        //std::cout << "searching " << bbb << std::endl;
        //std::cout << "inserting " << ddd << std::endl;
        /*for (typename DofsArrayType::iterator dof_iterator = Doftemp.begin(); dof_iterator != Doftemp.end(); ++dof_iterator)
        {
                KRATOS_WATCH( *dof_iterator )
        }
        std::cout << "DofTemp before Unique" << Doftemp.size() << std::endl;
         */
        //ccc = GetTickCount();
        Doftemp.Unique();
        //std::cout << "DofTemp after Unique" << Doftemp.size() << std::endl;
        BaseType::mDofSet = Doftemp;

        //ddd = GetTickCount() - ccc;
        //std::cout << "Unique " << ddd << std::endl;

        //throws an execption if there are no Degrees of freedom involved in the analysis
        if (BaseType::mDofSet.size() == 0)
            KRATOS_THROW_ERROR( std::logic_error, "No degrees of freedom!", "" )

        BaseType::mDofSetIsInitialized = true;
        if( this->GetEchoLevel() > 2 && r_model_part.GetCommunicator().MyPID() == 0)
        {
            std::cout << "finished setting up the dofs" << std::endl;
        }

        KRATOS_CATCH( "" )
    }

    //**************************************************************************
    //**************************************************************************

    void SetUpSystem(
        ModelPart& r_model_part
    )
    {
        // Set equation id for degrees of freedom
        // the free degrees of freedom are positioned at the beginning of the system,
        // while the fixed one are at the end (in opposite order).
        //
        // that means that if the EquationId is greater than "mEquationSystemSize"
        // the pointed degree of freedom is restrained
        //
        int free_id = 0;
        int fix_id = BaseType::mDofSet.size();

        for (typename DofsArrayType::iterator dof_iterator = BaseType::mDofSet.begin(); dof_iterator != BaseType::mDofSet.end(); ++dof_iterator)
            if (dof_iterator->IsFixed())
                dof_iterator->SetEquationId(--fix_id);
            else
                dof_iterator->SetEquationId(free_id++);

        BaseType::mEquationSystemSize = fix_id;

    }

    //**************************************************************************
    //**************************************************************************

    void ResizeAndInitializeVectors( typename TSchemeType::Pointer pScheme,
        TSystemMatrixPointerType& pA,
        TSystemVectorPointerType& pDx,
        TSystemVectorPointerType& pb,
        ElementsArrayType& rElements,
        ConditionsArrayType& rConditions,
        ProcessInfo& rCurrentProcessInfo
    )
    {
        KRATOS_TRY
        if (pA == NULL) //if the pointer is not initialized initialize it to an empty matrix
        {
            TSystemMatrixPointerType pNewA = TSystemMatrixPointerType(new TSystemMatrixType(0, 0));
            pA.swap(pNewA);
        }
        if (pDx == NULL) //if the pointer is not initialized initialize it to an empty matrix
        {
            TSystemVectorPointerType pNewDx = TSystemVectorPointerType(new TSystemVectorType(0));
            pDx.swap(pNewDx);
        }
        if (pb == NULL) //if the pointer is not initialized initialize it to an empty matrix
        {
            TSystemVectorPointerType pNewb = TSystemVectorPointerType(new TSystemVectorType(0));
            pb.swap(pNewb);
        }
        if (BaseType::mpReactionsVector == NULL) //if the pointer is not initialized initialize it to an empty matrix
        {
            TSystemVectorPointerType pNewReactionsVector = TSystemVectorPointerType(new TSystemVectorType(0));
            BaseType::mpReactionsVector.swap(pNewReactionsVector);
        }

        TSystemMatrixType& A = *pA;
        TSystemVectorType& Dx = *pDx;
        TSystemVectorType& b = *pb;

        //resizing the system vectors and matrix
        if (A.size1() == 0 || BaseType::GetReshapeMatrixFlag() == true) //if the matrix is not initialized
        {
            A.resize(BaseType::mEquationSystemSize, BaseType::mEquationSystemSize, false);
            ConstructMatrixStructure(pScheme, A, rElements, rConditions, rCurrentProcessInfo);
        }
        else
        {
            if (A.size1() != BaseType::mEquationSystemSize || A.size2() != BaseType::mEquationSystemSize)
            {
                KRATOS_WATCH( "it should not come here!!!!!!!! ... this is SLOW" )
                A.resize(BaseType::mEquationSystemSize, BaseType::mEquationSystemSize, true);
                ConstructMatrixStructure(pScheme, A, rElements, rConditions, rCurrentProcessInfo);
            }
        }
        if (Dx.size() != BaseType::mEquationSystemSize)
            Dx.resize(BaseType::mEquationSystemSize, false);
        if (b.size() != BaseType::mEquationSystemSize)
            b.resize(BaseType::mEquationSystemSize, false);

        //


        //if needed resize the vector for the calculation of reactions
        if (BaseType::mCalculateReactionsFlag == true)
        {
            unsigned int ReactionsVectorSize = BaseType::mDofSet.size() - BaseType::mEquationSystemSize;
            if (BaseType::mpReactionsVector->size() != ReactionsVectorSize)
                BaseType::mpReactionsVector->resize(ReactionsVectorSize, false);
        }

        KRATOS_CATCH( "" )

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

    void CalculateReactions(
        typename TSchemeType::Pointer pScheme,
        ModelPart& r_model_part,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b)
    {

        //refresh RHS to have the correct reactions
        this->BuildRHS(pScheme, r_model_part, b);

        int i;
        int systemsize = BaseType::mDofSet.size() - TSparseSpace::Size(*BaseType::mpReactionsVector);

        typename DofsArrayType::ptr_iterator it2;
        // KRATOS_WATCH( *BaseType::mpReactionsVector )
        //updating variables
        TSystemVectorType& ReactionsVector = *BaseType::mpReactionsVector;
        int num=1;
        for (it2 = BaseType::mDofSet.ptr_begin(); it2 != BaseType::mDofSet.ptr_end(); ++it2)
        {
            if ((*it2)->IsFixed())
            {
                i = (*it2)->EquationId();
                i -= systemsize;
                // KRATOS_WATCH( (*it2)->GetSolutionStepReactionValue() )
                // KRATOS_WATCH( ReactionsVector[i] )
                (*it2)->GetSolutionStepReactionValue() = ReactionsVector[i];
            }
	    else
	    {
	        (*it2)->GetSolutionStepReactionValue() = 0;
	    }

            num++;
        }

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
        // 			this->mReactionsVector = TSystemVectorType();

        this->mpLinearSystemSolver->Clear();

        if (this->GetEchoLevel() > 1)
        {
            std::cout << "ResidualBasedBuilderAndSolver Clear Function called" << std::endl;
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
    /*@} */
    /**@name Protected Operations*/
    /*@{ */


    //******************************************************************************************
    //******************************************************************************************

    virtual void ConstructMatrixStructure( typename TSchemeType::Pointer pScheme,
        TSystemMatrixType& A,
        ElementsContainerType& rElements,
        ConditionsArrayType& rConditions,
        ProcessInfo& rCurrentProcessInfo)
    {

        std::size_t equation_size = A.size1();
        std::vector<std::vector<std::size_t> > indices(equation_size);
        //std::vector<std::vector<std::size_t> > dirichlet_indices(TSystemSpaceType::Size1(mDirichletMatrix));

        Element::EquationIdVectorType ids(3, 0);
        for (typename ElementsContainerType::iterator i_element = rElements.begin(); i_element != rElements.end(); i_element++)
        {
             pScheme->EquationId( *(i_element.base()) , ids, rCurrentProcessInfo);

            for (std::size_t i = 0; i < ids.size(); i++)
                if (ids[i] < equation_size)
                {
                    std::vector<std::size_t>& row_indices = indices[ids[i]];
                    for (std::size_t j = 0; j < ids.size(); j++)
                        if (ids[j] < equation_size)
                        {
                            AddUnique(row_indices, ids[j]);
                            //indices[ids[i]].push_back(ids[j]);
                        }
                }

        }

        for (typename ConditionsArrayType::iterator i_condition = rConditions.begin(); i_condition != rConditions.end(); i_condition++)
        {
            pScheme->Condition_EquationId( *(i_condition.base()), ids, rCurrentProcessInfo);
            for (std::size_t i = 0; i < ids.size(); i++)
                if (ids[i] < equation_size)
                {
                    std::vector<std::size_t>& row_indices = indices[ids[i]];
                    for (std::size_t j = 0; j < ids.size(); j++)
                        if (ids[j] < equation_size)
                        {
                            AddUnique(row_indices, ids[j]);
                            //	indices[ids[i]].push_back(ids[j]);
                        }
                }
        }

        //allocating the memory needed
        int data_size = 0;
        for (std::size_t i = 0; i < indices.size(); i++)
        {
            data_size += indices[i].size();
        }
        A.reserve(data_size, false);

        //filling with zero the matrix (creating the structure)
        Timer::Start("MatrixStructure");
#ifndef _OPENMP
        for (std::size_t i = 0; i < indices.size(); i++)
        {
            std::vector<std::size_t>& row_indices = indices[i];
            std::sort(row_indices.begin(), row_indices.end());

            for (std::vector<std::size_t>::iterator it = row_indices.begin(); it != row_indices.end(); it++)
            {
                A.push_back(i, *it, 0.00);
            }
            row_indices.clear();
        }
#else
        int number_of_threads = omp_get_max_threads();
        vector<unsigned int> matrix_partition;
        CreatePartition(number_of_threads, indices.size(), matrix_partition);
        if (this->GetEchoLevel() > 2)
        {
            KRATOS_WATCH( matrix_partition )
        }
        for (int k = 0; k < number_of_threads; k++)
        {
            #pragma omp parallel
            if (omp_get_thread_num() == k)
            {
                for (std::size_t i = matrix_partition[k]; i < matrix_partition[k + 1]; i++)
                {
                    std::vector<std::size_t>& row_indices = indices[i];
                    std::sort(row_indices.begin(), row_indices.end());

                    for (std::vector<std::size_t>::iterator it = row_indices.begin(); it != row_indices.end(); it++)
                    {
                        A.push_back(i, *it, 0.00);
                    }
                    row_indices.clear();
                }
            }
        }
#endif
        Timer::Stop("MatrixStructure");
    }
    
    
    
    //******************************************************************************************
    //******************************************************************************************

    virtual void ConstructMatrixStructure(
        TSystemMatrixType& A,
        ElementsContainerType& rElements,
        ConditionsArrayType& rConditions,
        ProcessInfo& rCurrentProcessInfo)
    {

        std::size_t equation_size = A.size1();
        std::vector<std::vector<std::size_t> > indices(equation_size);
        //std::vector<std::vector<std::size_t> > dirichlet_indices(TSystemSpaceType::Size1(mDirichletMatrix));

        Element::EquationIdVectorType ids(3, 0);
        for (typename ElementsContainerType::iterator i_element = rElements.begin(); i_element != rElements.end(); i_element++)
        {
            (i_element)->EquationIdVector(ids, rCurrentProcessInfo);

            for (std::size_t i = 0; i < ids.size(); i++)
                if (ids[i] < equation_size)
                {
                    std::vector<std::size_t>& row_indices = indices[ids[i]];
                    for (std::size_t j = 0; j < ids.size(); j++)
                        if (ids[j] < equation_size)
                        {
                            AddUnique(row_indices, ids[j]);
                            //indices[ids[i]].push_back(ids[j]);
                        }
                }

        }

        for (typename ConditionsArrayType::iterator i_condition = rConditions.begin(); i_condition != rConditions.end(); i_condition++)
        {
            (i_condition)->EquationIdVector(ids, rCurrentProcessInfo);
            for (std::size_t i = 0; i < ids.size(); i++)
                if (ids[i] < equation_size)
                {
                    std::vector<std::size_t>& row_indices = indices[ids[i]];
                    for (std::size_t j = 0; j < ids.size(); j++)
                        if (ids[j] < equation_size)
                        {
                            AddUnique(row_indices, ids[j]);
                            //  indices[ids[i]].push_back(ids[j]);
                        }
                }
        }

        //allocating the memory needed
        int data_size = 0;
        for (std::size_t i = 0; i < indices.size(); i++)
        {
            data_size += indices[i].size();
        }
        A.reserve(data_size, false);

        //filling with zero the matrix (creating the structure)
        Timer::Start("MatrixStructure");
#ifndef _OPENMP
        for (std::size_t i = 0; i < indices.size(); i++)
        {
            std::vector<std::size_t>& row_indices = indices[i];
            std::sort(row_indices.begin(), row_indices.end());

            for (std::vector<std::size_t>::iterator it = row_indices.begin(); it != row_indices.end(); it++)
            {
                A.push_back(i, *it, 0.00);
            }
            row_indices.clear();
        }
#else
        int number_of_threads = omp_get_max_threads();
        vector<unsigned int> matrix_partition;
        CreatePartition(number_of_threads, indices.size(), matrix_partition);
        if (this->GetEchoLevel() > 2)
        {
            KRATOS_WATCH( matrix_partition )
        }
        for (int k = 0; k < number_of_threads; k++)
        {
            #pragma omp parallel
            if (omp_get_thread_num() == k)
            {
                for (std::size_t i = matrix_partition[k]; i < matrix_partition[k + 1]; i++)
                {
                    std::vector<std::size_t>& row_indices = indices[i];
                    std::sort(row_indices.begin(), row_indices.end());

                    for (std::vector<std::size_t>::iterator it = row_indices.begin(); it != row_indices.end(); it++)
                    {
                        A.push_back(i, *it, 0.00);
                    }
                    row_indices.clear();
                }
            }
        }
#endif
        Timer::Stop("MatrixStructure");
    }    
    
    

    //******************************************************************************************
    //******************************************************************************************

    void AssembleLHS(
        TSystemMatrixType& A,
        LocalSystemMatrixType& LHS_Contribution,
        Element::EquationIdVectorType& EquationId
    )
    {
        unsigned int local_size = LHS_Contribution.size1();

        for (unsigned int i_local = 0; i_local < local_size; i_local++)
        {
            unsigned int i_global = EquationId[i_local];
            if (i_global < BaseType::mEquationSystemSize)
            {
                for (unsigned int j_local = 0; j_local < local_size; j_local++)
                {
                    unsigned int j_global = EquationId[j_local];
                    if (j_global < BaseType::mEquationSystemSize)
                        A(i_global, j_global) += LHS_Contribution(i_local, j_local);
                }
            }
        }
    }


    //******************************************************************************************
    //******************************************************************************************

    void AssembleRHS(
        TSystemVectorType& b,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId
    )
    {
        unsigned int local_size = RHS_Contribution.size();

        if (BaseType::mCalculateReactionsFlag == false) //if we don't need to calculate reactions
        {
            for (unsigned int i_local = 0; i_local < local_size; i_local++)
            {
                unsigned int i_global = EquationId[i_local];
                if (i_global < BaseType::mEquationSystemSize) //on "free" DOFs
                {
                    // ASSEMBLING THE SYSTEM VECTOR
                    b[i_global] += RHS_Contribution[i_local];
                }
            }
        }
        else   //when the calculation of reactions is needed
        {
            TSystemVectorType& ReactionsVector = *BaseType::mpReactionsVector;
            for (unsigned int i_local = 0; i_local < local_size; i_local++)
            {
                unsigned int i_global = EquationId[i_local];
                if (i_global < BaseType::mEquationSystemSize) //on "free" DOFs
                {
                    // ASSEMBLING THE SYSTEM VECTOR
                    b[i_global] += RHS_Contribution[i_local];
                }
                else   //on "fixed" DOFs
                {
                    // Assembling the Vector of REACTIONS
                    ReactionsVector[i_global - BaseType::mEquationSystemSize] -= RHS_Contribution[i_local];
                }
            }
        }
    }


    //******************************************************************************************
    //******************************************************************************************

    void AssembleLHS_CompleteOnFreeRows(
        TSystemMatrixType& A,
        LocalSystemMatrixType& LHS_Contribution,
        Element::EquationIdVectorType& EquationId
    )
    {
        unsigned int local_size = LHS_Contribution.size1();
        for (unsigned int i_local = 0; i_local < local_size; i_local++)
        {
            unsigned int i_global = EquationId[i_local];
            if (i_global < BaseType::mEquationSystemSize)
            {
                for (unsigned int j_local = 0; j_local < local_size; j_local++)
                {
                    int j_global = EquationId[j_local];

                    A(i_global, j_global) += LHS_Contribution(i_local, j_local);
                }
            }
        }
    }


    //******************************************************************************************
    //******************************************************************************************

    inline void AddUnique(std::vector<std::size_t>& v, const std::size_t& candidate)
    {
        std::vector<std::size_t>::iterator i = v.begin();
        std::vector<std::size_t>::iterator endit = v.end();
        while (i != endit && (*i) != candidate)
        {
            i++;
        }
        if (i == endit)
        {
            v.push_back(candidate);
        }

    }

    //******************************************************************************************
    //******************************************************************************************

    inline void CreatePartition(unsigned int number_of_threads, const int number_of_rows, vector<unsigned int>& partitions)
    {
        partitions.resize(number_of_threads + 1);
        int partition_size = number_of_rows / number_of_threads;
        partitions[0] = 0;
        partitions[number_of_threads] = number_of_rows;
        for (unsigned int i = 1; i < number_of_threads; i++)
            partitions[i] = partitions[i - 1] + partition_size;
    }


    //******************************************************************************************
    //******************************************************************************************


#ifdef _OPENMP

    void Assemble(
        TSystemMatrixType& A,
        TSystemVectorType& b,
        const LocalSystemMatrixType& LHS_Contribution,
        const LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        std::vector< omp_lock_t >& lock_array
    )
    {
        unsigned int local_size = RHS_Contribution.size();

        for (unsigned int i_local = 0; i_local < local_size; i_local++)
        {
            unsigned int i_global = EquationId[i_local];

            if (i_global < BaseType::mEquationSystemSize)
            {
                omp_set_lock(&lock_array[i_global]);

                b[i_global] += RHS_Contribution(i_local);

                for (unsigned int j_local = 0; j_local < local_size; j_local++)
                {
                    unsigned int j_global = EquationId[j_local];
                    if (j_global < BaseType::mEquationSystemSize)
                    {
                        A(i_global, j_global) += LHS_Contribution(i_local, j_local);
                    }
                }

                omp_unset_lock(&lock_array[i_global]);


            }
            //note that computation of reactions is not performed here!
        }
    }

    //******************************************************************************************
    //******************************************************************************************

    void AssembleLHS(
        TSystemMatrixType& A,
        const LocalSystemMatrixType& LHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        std::vector< omp_lock_t >& lock_array
    )
    {
        unsigned int local_size = LHS_Contribution.size1();

        for (unsigned int i_local = 0; i_local < local_size; i_local++)
        {
            unsigned int i_global = EquationId[i_local];

            if (i_global < BaseType::mEquationSystemSize)
            {
                omp_set_lock(&lock_array[i_global]);

                for (unsigned int j_local = 0; j_local < local_size; j_local++)
                {
                    unsigned int j_global = EquationId[j_local];
                    if (j_global < BaseType::mEquationSystemSize)
                    {
                        A(i_global, j_global) += LHS_Contribution(i_local, j_local);
                    }
                }

                omp_unset_lock(&lock_array[i_global]);


            }
            //note that computation of reactions is not performed here!
        }
    }

    //******************************************************************************************
    //******************************************************************************************
    void AssembleRHS(
        TSystemVectorType& b,
        const LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        std::vector< omp_lock_t >& lock_array
    )
    {

        unsigned int local_size = RHS_Contribution.size();

        if (BaseType::mCalculateReactionsFlag == false) //if we don't need to calculate reactions
	    {
			for (unsigned int i_local = 0; i_local < local_size; i_local++)
			{
				unsigned int i_global = EquationId[i_local];

				if (i_global < BaseType::mEquationSystemSize)
				{
					omp_set_lock(&lock_array[i_global]);

		            b[i_global] += RHS_Contribution(i_local);
			        
		            omp_unset_lock(&lock_array[i_global]);
				}
		        //note that computation of reactions is not performed here!
            }

	    }
        else   //when the calculation of reactions is needed
		{
			TSystemVectorType& ReactionsVector = *BaseType::mpReactionsVector;

            for (unsigned int i_local = 0; i_local < local_size; i_local++)
			{
				unsigned int i_global = EquationId[i_local];

                if (i_global < BaseType::mEquationSystemSize) //on "free" DOFs
				{
					omp_set_lock(&lock_array[i_global]);

                    b[i_global] += RHS_Contribution[i_local];

                    omp_unset_lock(&lock_array[i_global]);
				}
                else   //on "fixed" DOFs
				{
					/*std::cout << " - LK.size() = " << lock_array.size();
					std::cout << " - REAC.size() = " << ReactionsVector.size();
					std::cout << " - TARGET ID = " << i_global - BaseType::mEquationSystemSize;
					std::cout << " -system = " << BaseType::mEquationSystemSize << std::endl;*/

					// bug fixed: lock_array now has dimension(=num_free + num_fixed) if 
					// calculate reaction is needed
					omp_set_lock(&lock_array[i_global/* - BaseType::mEquationSystemSize*/]);

                    ReactionsVector[i_global - BaseType::mEquationSystemSize] -= RHS_Contribution[i_local];
					
					omp_unset_lock(&lock_array[i_global/* - BaseType::mEquationSystemSize*/]);
				}
			}
		}
	}


#endif

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
    GlobalSystemComponentsType mGlobalSystem;
    /*@} */
    /**@name Private Operators*/
    /*@{ */
    /*@} */
    /**@name Private Operations*/
    /*@{ */

    //**************************************************************************
    //**************************************************************************
    void ResizeAndInitializeComponents(ModelPart& r_model_part, typename TSchemeType::Pointer pScheme)
    {

        // Element components resize:
      
        //contributions to the element local system
	LocalSystemComponentsType& ElementLocalSystem = pScheme->GetLocalSystemComponents();

	bool LHS_Element_Components_Set = mGlobalSystem.Are_LHS_Element_Components_Set();
	if( LHS_Element_Components_Set ){

	  std::vector<TSystemMatrixType>& rLHS_GlobalElementComponents = mGlobalSystem.GetLHS_Element_Components();
	
	  ElementLocalSystem.SetLHS_Element_Variables( mGlobalSystem.GetLHS_Element_Variables() );

	  if( mGlobalSystem.GetLHS_Element_Variables().size() != rLHS_GlobalElementComponents.size() )
	    rLHS_GlobalElementComponents.resize(mGlobalSystem.GetLHS_Element_Variables().size());
	  
	  this->ResizeAndInitializeMatrices(r_model_part, rLHS_GlobalElementComponents);
	}


	bool RHS_Element_Components_Set = mGlobalSystem.Are_RHS_Element_Components_Set();
	if( RHS_Element_Components_Set ){

	  std::vector<TSystemVectorType>& rRHS_GlobalElementComponents = mGlobalSystem.GetRHS_Element_Components();

	  ElementLocalSystem.SetRHS_Element_Variables( mGlobalSystem.GetRHS_Element_Variables() );

	  if( mGlobalSystem.GetRHS_Element_Variables().size() != rRHS_GlobalElementComponents.size() )
	    rRHS_GlobalElementComponents.resize(mGlobalSystem.GetRHS_Element_Variables().size());

	  this->ResizeAndInitializeVectorsByComponents(r_model_part, rRHS_GlobalElementComponents);

	}

	// Condition components resize:
	
	//contributions to the condition local system
	LocalSystemComponentsType& ConditionLocalSystem = pScheme->GetLocalSystemComponents();

	bool LHS_Condition_Components_Set = mGlobalSystem.Are_LHS_Condition_Components_Set();
	if( LHS_Condition_Components_Set ){
	  
	  std::vector<TSystemMatrixType>& rLHS_GlobalConditionComponents = mGlobalSystem.GetLHS_Condition_Components();

	  ConditionLocalSystem.SetLHS_Condition_Variables( mGlobalSystem.GetLHS_Condition_Variables() );
	  
	  if( mGlobalSystem.GetLHS_Condition_Variables().size() != rLHS_GlobalConditionComponents.size() )
	    rLHS_GlobalConditionComponents.resize(mGlobalSystem.GetLHS_Condition_Variables().size());
	  
	  this->ResizeAndInitializeMatrices(r_model_part, rLHS_GlobalConditionComponents);
	  
	}
	
	
	bool RHS_Condition_Components_Set = mGlobalSystem.Are_RHS_Condition_Components_Set();
	if( RHS_Condition_Components_Set ){
	  
	  std::vector<TSystemVectorType>& rRHS_GlobalConditionComponents = mGlobalSystem.GetRHS_Condition_Components();
	  ConditionLocalSystem.SetRHS_Condition_Variables( mGlobalSystem.GetRHS_Condition_Variables() );
	  
	  if( mGlobalSystem.GetRHS_Condition_Variables().size() != rRHS_GlobalConditionComponents.size() )
	    rRHS_GlobalConditionComponents.resize(mGlobalSystem.GetRHS_Condition_Variables().size());
	  
	  this->ResizeAndInitializeVectorsByComponents(r_model_part, rRHS_GlobalConditionComponents);
	}

    }

    //**************************************************************************
    //**************************************************************************

    void ResizeAndInitializeMatrices(ModelPart& r_model_part, std::vector<TSystemMatrixType>& rGlobalSystemMatrices)
    {

      ElementsArrayType& rElements     = r_model_part.Elements(); 
      ConditionsArrayType& rConditions = r_model_part.Conditions(); 
      ProcessInfo& rCurrentProcessInfo = r_model_part.GetProcessInfo();
      
      for( unsigned int i=0; i< rGlobalSystemMatrices.size(); i++ )
	{
	  //resizing the system vectors and matrix
	  if( rGlobalSystemMatrices[i].size1() == 0 || BaseType::GetReshapeMatrixFlag() == true ) //if the matrix is not initialized
	    {
	      rGlobalSystemMatrices[i].resize( BaseType::mEquationSystemSize, BaseType::mEquationSystemSize, false );
	      this->ConstructMatrixStructure(rGlobalSystemMatrices[i], rElements, rConditions, rCurrentProcessInfo);
	    }
	  else
	    {
	      if( rGlobalSystemMatrices[i].size1() != BaseType::mEquationSystemSize || rGlobalSystemMatrices[i].size2() != BaseType::mEquationSystemSize)
		{
		  KRATOS_WATCH( "it should not come here!!!!!!!! ... this is SLOW" )
		  rGlobalSystemMatrices[i].resize(BaseType::mEquationSystemSize, BaseType::mEquationSystemSize, true);
		  this->ConstructMatrixStructure(rGlobalSystemMatrices[i], rElements, rConditions, rCurrentProcessInfo);
		}
	    }
	}
      
    }


    //**************************************************************************
    //**************************************************************************

    void ResizeAndInitializeVectorsByComponents(ModelPart& r_model_part, std::vector<TSystemVectorType>& rGlobalSystemVectors)
    {
      
      for( unsigned int i=0; i< rGlobalSystemVectors.size(); i++ )
	{
	  if( rGlobalSystemVectors[i].size() != BaseType::mEquationSystemSize )
	    rGlobalSystemVectors[i].resize(BaseType::mEquationSystemSize, false);

	  rGlobalSystemVectors[i].clear(); // set components to zero
	}
      
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

}; /* Class ComponentWiseBuilderAndSolver */

/*@} */

/**@name Type Definitions */
/*@{ */
/*@} */

} /* namespace Kratos.*/

#endif /* KRATOS_COMPONENT_WISE_BUILDER_AND_SOLVER  defined */

