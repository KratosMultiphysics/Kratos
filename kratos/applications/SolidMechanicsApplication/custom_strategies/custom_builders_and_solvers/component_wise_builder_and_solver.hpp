//
//   Project Name:        KratosSolidMechanicsApplication $
//   Last modified by:    $Author:            JMCarbonell $
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
#include "custom_strategies/custom_builders_and_solvers/residual_based_builder_and_solver.hpp"

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
    : public ResidualBasedBuilderAndSolver< TSparseSpace, TDenseSpace, TLinearSolver >
{
public:

    /**@name Type Definitions */
    /*@{ */
    //typedef boost::shared_ptr< ComponentWiseBuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver> > Pointer;
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
        : ResidualBasedBuilderAndSolver< TSparseSpace, TDenseSpace, TLinearSolver >(pNewLinearSystemSolver)
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

