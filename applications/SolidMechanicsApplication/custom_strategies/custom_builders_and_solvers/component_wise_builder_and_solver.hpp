//
//   Project Name:        KratosSolidMechanicsApplication $
//   Last modified by:    $Author:            JMCarbonell $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_COMPONENT_BASED_BUILDER_AND_SOLVER )
#define  KRATOS_COMPONENT_BASED_BUILDER_AND_SOLVER


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
class ComponentBasedBuilderAndSolver
    : public ResidualBasedBuilderAndSolver< TSparseSpace, TDenseSpace, TLinearSolver >
{
public:

    struct GlobalSystemContributions
    {
    private:

      TSystemMatrixType *mpLeftHandSideContribution;
      TSystemVectorType *mpRightHandSideContribution;

      //elements
      std::vector<TSystemMatrixType> *mpLHS_Element_Components;
      const std::vector< Variable< TSystemMatrixType > > *mpLHS_Element_Variables;

      std::vector<TSystemVectorType> *mpRHS_Element_Components;
      const std::vector< Variable< TSystemVectorType > > *mpRHS_Element_Variables;
      
      //conditions
      std::vector<TSystemMatrixType> *mpLHS_Condition_Components;
      const std::vector< Variable< TSystemMatrixType > > *mpLHS_Condition_Variables;

      std::vector<TSystemVectorType> *mpRHS_Condition_Components;
      const std::vector< Variable< TSystemVectorType > > *mpRHS_Condition_Variables;
      
    public:
      
      //setting pointer variables
      void SetLeftHandSideContribution  ( TSystemMatrixType& rLeftHandSideContribution )  { mpLeftHandSideContribution = &rLeftHandSideContribution; };
      void SetRightHandSideContribution ( TSystemVectorType& rRightHandSideContribution ) { mpRightHandSideContribution = &rRightHandSideContribution; };

      //elements
      void SetLHS_Element_Components ( std::vector<TSystemMatrixType>& rLHS_Element_Components ) { mpLHS_Element_Components = &rLHS_Element_Components; };
      void SetLHS_Element_Variables     ( std::vector< Variable< LocalSystemMatrixType > >& rLHS_Element_Variables ) { mpLHS_Element_Variables = &rLHS_Element_Variables; };
      void SetRHS_Element_Components ( std::vector<TSystemVectorType>& rRHS_Element_Components ) { mpRHS_Element_Components = &rRHS_Element_Components; };
      void SetRHS_Element_Variables     ( std::vector< Variable< LocalSystemVectorType > >& rRHS_Element_Variables ) { mpRHS_Element_Variables = &rRHS_Element_Variables; };

      //conditions
      void SetLHS_Condition_Components ( std::vector<TSystemMatrixType>& rLHS_Condition_Components ) { mpLHS_Condition_Components = &rLHS_Condition_Components; };
      void SetLHS_Condition_Variables     ( std::vector< Variable< LocalSystemMatrixType > >& rLHS_Condition_Variables ) { mpLHS_Condition_Variables = &rLHS_Condition_Variables; };
      void SetRHS_Condition_Components ( std::vector<TSystemVectorType>& rRHS_Condition_Components ) { mpRHS_Condition_Components = &rRHS_Condition_Components; };
      void SetRHS_Condition_Variables     ( std::vector< Variable< LocalSystemVectorType > >& rRHS_Condition_Variables ) { mpRHS_Condition_Variables = &rRHS_Condition_Variables; };

      //getting pointer variables
      TSystemMatrixType& GetLeftHandSideContribution  () { return *mpLeftHandSideContribution; };
      TSystemVectorType& GetRightHandSideContribution () { return *mpRightHandSideContribution; };

      //elements
      std::vector<TSystemMatrixType>& GetLHS_Element_Components() { return *mpLHS_Element_Components; };
      std::vector< Variable< LocalSystemMatrixType > >& GetLHS_Element_Variables() { return *mpLHS_Element_Variables; };
      std::vector<TSystemVectorType>& GetRHS_Element_Components() { return *mpRHS_Element_Components; };
      std::vector< Variable< LocalSystemVectorType > >& GetRHS_Element_Variables() { return *mpRHS_Element_Variables; };

      //conditions
      std::vector<TSystemMatrixType>& GetLHS_Condition_Components() { return *mpLHS_Condition_Components; };
      std::vector< Variable< LocalSystemMatrixType > >& GetLHS_Condition_Variables() { return *mpLHS_Condition_Variables; };
      std::vector<TSystemVectorType>& GetRHS_Condition_Components() { return *mpRHS_Condition_Components; };
      std::vector< Variable< LocalSystemVectorType > >& GetRHS_Condition_Variables() { return *mpRHS_Condition_Variables; };

    }

    /**@name Type Definitions */
    /*@{ */
    //typedef boost::shared_ptr< ComponentBasedBuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver> > Pointer;
    KRATOS_CLASS_POINTER_DEFINITION(ComponentBasedBuilderAndSolver);

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
    ComponentBasedBuilderAndSolver(
        typename TLinearSolver::Pointer pNewLinearSystemSolver)
        : ResidualBasedBuilderAndSolver< TSparseSpace, TDenseSpace, TLinearSolver >(pNewLinearSystemSolver)
    {
    }

    /** Destructor.
     */
    virtual ~ComponentBasedBuilderAndSolver()
    {
    }


    /*@} */
    /**@name Operators
     */
    /*@{ */

    //**************************************************************************
    //**************************************************************************

    void Build(
        typename TSchemeType::Pointer pScheme,
        ModelPart& r_model_part,
	GlobalSystemContributions& rGlobalSystem)
    {
        KRATOS_TRY

        if (!pScheme)
            KRATOS_ERROR(std::runtime_error, "No scheme provided!", "");

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

        //contributions to the element local system
	ComponentBasedBossakScheme::SystemContributions ElementLocalSystem;
	
        LocalSystemMatrixType LeftHandSideElementContribution  = LocalSystemMatrixType(0, 0);
        LocalSystemVectorType RightHandSideElementContribution = LocalSystemVectorType(0);

	ElementLocalSystem.SetLeftHandSideContribution( LeftHandSideElementContribution );
	ElementLocalSystem.SetRightHandSideContribution( RightHandSideElementContribution );

	std::vector<LocalSystemMatrixType> rLHS_LocalElementComponents; 
	ElementLocalSystem.SetLHS_Components(rLHS_LocalElementComponents);

	std::vector<LocalSystemVectorType> rRHS_LocalElementComponents; 
	ElementLocalSystem.SetRHS_Components(rRHS_LocalElementComponents);

        //contributions to the global system
	std::vector<TSystemMatrixType>& rLHS_GlobalElementComponents = rGlobalSystem.GetLHS_Element_Components();
	ElementLocalSystem.SetLHS_Element_Variables( rGlobalSystem.GetLHS_Element_Variables() );

	if( rGlobalSystem.GetLHS_Element_Variables().size() != rLHS_GlobalElementComponents.size() )
	  rLHS_GlobalElementComponents.resize(rGlobalSystem.GetLHS_Element_Variables().size());

	for( unsigned int i=0; i<rLHS_GlobalElementComponents.size(); i++ )
	  {
	    rLHS_GlobalElementComponents[i] = LocalSystemMatrixType(0, 0);
	  }

	std::vector<TSystemVectorType>& rRHS_GlobalElementComponents = rGlobalSystem.GetRHS_Element_Components();
	ElementLocalSystem.SetRHS_Element_Variables( rGlobalSystem.GetRHS_Element_Variables() );

	if( rGlobalSystem.GetRHS_Element_Variables().size() != rRHS_GlobalElementComponents.size() )
	  rRHS_GlobalElementComponents.resize(rGlobalSystem.GetRHS_Element_Variables().size());

	for( unsigned int i=0; i<rRHS_GlobalElementComponents.size(); i++ )
	  {
	    rRHS_GlobalElementComponents[i] = LocalSystemVectorType(0);
	  }


	TSystemMatrixType& rGlobalLeftHandSideContribution  = rGlobalSystem.GetLeftHandSideContribution();
	TSystemVectorType& rGlobalRightHandSideContribution = rGlobalSystem.GetRightHandSideContribution();

        //vector containing the localization in the system of the different terms
	Element::EquationIdVectorType EquationId;

        //double StartTime = GetTickCount();

        // assemble all elements
#ifndef _OPENMP
        ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

        for (typename ElementsArrayType::ptr_iterator it = pElements.ptr_begin(); it != pElements.ptr_end(); ++it)
        {
            //calculate elemental contribution
            pScheme->CalculateSystemContributions(*it, ElementLocalSystem, EquationId, CurrentProcessInfo);

            //assemble the elemental contribution
            AssembleLHS(rGlobalLeftHandSideContribution, LeftHandSideElementContribution, EquationId);
            AssembleRHS(rGlobalRightHandSideContribution, RightHandSideElementContribution, EquationId);

	    for( unsigned int i=0; i<rLHS_GlobalComponents.size(); i++ )
	      {	    
		AssembleLHS(rLHS_GlobalElementElementComponents[i], rLHS_LocalElementComponents[i], EquationId);
	      }

	    for( unsigned int i=0; i<rRHS_GlobalElementComponents.size(); i++ )
	      {
		AssembleRHS(rRHS_GlobalElementComponents[i], rRHS_LocalElementComponents[i], EquationId);
	      }

            // clean local elemental memory
            pScheme->CleanMemory(*it);
        }
        //double EndTime = GetTickCount();

        //std::cout << "total time " << EndTime - StartTime << std::endl;
        //std::cout << "writing in the system matrix " << ccc << std::endl;
        //std::cout << "calculating the elemental contrib " << ddd << std::endl;


        //contributions to the condition local system
	ComponentBasedBossakScheme::SystemContributions ConditionLocalSystem;

        LocalSystemMatrixType LeftHandSideConditionContribution  = LocalSystemMatrixType(0, 0);
        LocalSystemVectorType RightHandSideConditionContribution = LocalSystemVectorType(0);

	ConditionLocalSystem.SetLeftHandSideContribution( LeftHandSideConditionContribution );
	ConditionLocalSystem.SetRightHandSideContribution( RightHandSideConditionContribution );

	std::vector<LocalSystemMatrixType> rLHS_LocalConditionComponents; 
	ConditionLocalSystem.SetLHS_Components(rLHS_LocalConditionComponents);

	std::vector<LocalSystemVectorType> rRHS_LocalConditionComponents; 
	ConditionLocalSystem.SetRHS_Components(rRHS_LocalConditionComponents);

        //contributions to the global system
	std::vector<TSystemMatrixType>& rLHS_GlobalConditionComponents = rGlobalSystem.GetLHS_Condition_Components();
	ConditionLocalSystem.SetLHS_Condition_Variables( rGlobalSystem.GetLHS_Condition_Variables() );

	if( rGlobalSystem.GetLHS_Condition_Variables().size() != rLHS_GlobalConditionComponents.size() )
	  rLHS_GlobalConditionComponents.resize(rGlobalSystem.GetLHS_Condition_Variables().size());

	for( unsigned int i=0; i<rLHS_GlobalConditionComponents.size(); i++ )
	  {
	    rLHS_GlobalConditionComponents[i] = LocalSystemMatrixType(0, 0);
	  }

	std::vector<TSystemVectorType>& rRHS_GlobalConditionComponents = rGlobalSystem.GetRHS_Condition_Components();
	ConditionLocalSystem.SetRHS_Condition_Variables( rGlobalSystem.GetRHS_Condition_Variables() );

	if( rGlobalSystem.GetRHS_Condition_Variables().size() != rRHS_GlobalConditionComponents.size() )
	  rRHS_GlobalConditionComponents.resize(rGlobalSystem.GetRHS_Condition_Variables().size());

	for( unsigned int i=0; i<rRHS_GlobalConditionComponents.size(); i++ )
	  {
	    rRHS_GlobalConditionComponents[i] = LocalSystemVectorType(0);
	  }


        // assemble all conditions
        for (typename ConditionsArrayType::ptr_iterator it = ConditionsArray.ptr_begin(); it != ConditionsArray.ptr_end(); ++it)
        {
            //calculate condition contribution
  	    pScheme->Condition_CalculateSystemContributions(*it, ConditionLocalSytem, EquationId, CurrentProcessInfo);

            //assemble the condition contribution
            AssembleLHS(rGlobalLeftHandSideContribution, LeftHandSideConditionContribution, EquationId);
            AssembleRHS(rGlobalRightHandSideContribution, RightHandSideConditionContribution, EquationId);

	    for( unsigned int i=0; i<rLHS_GlobalComponents.size(); i++ )
	      {	    
		AssembleLHS(rLHS_GlobalConditionComponents[i], rLHS_LocalConditionComponents[i], EquationId);
	      }

	    for( unsigned int i=0; i<rRHS_GlobalComponents.size(); i++ )
	      {
		AssembleRHS(rRHS_GlobalConditionComponents[i], rRHS_LocalConditionComponents[i], EquationId);
	      }
        }

#else
        //creating an array of lock variables of the size of the system matrix
        std::vector< omp_lock_t > lock_array(rGlobalLeftHandSideContribution.size1());

        int A_size = rGlobalLeftHandSideContribution.size1();
        for (int i = 0; i < A_size; i++)
            omp_init_lock(&lock_array[i]);

        //create a partition of the element array
        int number_of_threads = omp_get_max_threads();

        vector<unsigned int> element_partition;
        CreatePartition(number_of_threads, pElements.size(), element_partition);
        if( this->GetEchoLevel() > 2 && r_model_part.GetCommunicator().MyPID() == 0)
        {
            KRATOS_WATCH(number_of_threads);
            KRATOS_WATCH(element_partition);
        }


        double start_prod = omp_get_wtime();

        #pragma omp parallel for
        for (int k = 0; k < number_of_threads; k++)
        {

	    //contributions to the element local system
	    ComponentBasedBossakScheme::SystemContributions ElementLocalSystem;
	  
	    LocalSystemMatrixType LeftHandSideElementContribution  = LocalSystemMatrixType(0, 0);
	    LocalSystemVectorType RightHandSideElementContribution = LocalSystemVectorType(0);
	    
	    ElementLocalSystem.SetLeftHandSideContribution( LeftHandSideElementContribution );
	    ElementLocalSystem.SetRightHandSideContribution( RightHandSideElementContribution );

	    std::vector<LocalSystemMatrixType> rLHS_LocalElementComponents; 
	    ElementLocalSystem.SetLHS_Components(rLHS_LocalElementComponents);
	    
	    std::vector<LocalSystemVectorType> rRHS_LocalElementComponents; 
	    ElementLocalSystem.SetRHS_Components(rRHS_LocalElementComponents);
	    
    
	    //contributions to the global system
	    std::vector<TSystemMatrixType>& rLHS_GlobalElementComponents = rGlobalSystem.GetLHS_Element_Components();
	    ElementLocalSystem.SetLHS_Element_Variables( rGlobalSystem.GetLHS_Element_Variables() );
	    
	    if( rGlobalSystem.GetLHS_Element_Variables().size() != rLHS_GlobalElementComponents.size() )
	      rLHS_GlobalElementComponents.resize(rGlobalSystem.GetLHS_Element_Variables().size());
	    
	    for( unsigned int i=0; i<rLHS_GlobalElementComponents.size(); i++ )
	      {
		rLHS_GlobalElementComponents[i] = LocalSystemMatrixType(0, 0);
	      }
	    
	    std::vector<TSystemVectorType>& rRHS_GlobalElementComponents = rGlobalSystem.GetRHS_Element_Components();
	    ElementLocalSystem.SetRHS_Element_Variables( rGlobalSystem.GetRHS_Element_Variables() );
	    
	    if( rGlobalSystem.GetRHS_Element_Variables().size() != rRHS_GlobalElementComponents.size() )
	      rRHS_GlobalElementComponents.resize(rGlobalSystem.GetRHS_Element_Variables().size());
	    
	    for( unsigned int i=0; i<rRHS_GlobalElementComponents.size(); i++ )
	      {
		rRHS_GlobalElementComponents[i] = LocalSystemVectorType(0);
	      }
	    

            ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();
            typename ElementsArrayType::ptr_iterator it_begin = pElements.ptr_begin() + element_partition[k];
            typename ElementsArrayType::ptr_iterator it_end = pElements.ptr_begin() + element_partition[k + 1];

            // assemble all elements
            for (typename ElementsArrayType::ptr_iterator it = it_begin; it != it_end; ++it)
            {

                //calculate elemental contribution
	        pScheme->CalculateSystemContributions(*it, ElementLocalSystem, EquationId, CurrentProcessInfo);

                //assemble the elemental contribution
                Assemble(rGlobalLeftHandSideContribution, rGlobalRightHandSideContribution, 
			 LeftHandSideElementContribution, RightHandSideElementContribution, EquationId, lock_array);

		//assemble the elemental contribution
		//AssembleLHS(rGlobalSystem.LeftHandSideContribution, ElementLocalSystem.LeftHandSideContribution, EquationId);
		//AssembleRHS(rGlobalSystem.RightHandSideContribution, ElementLocalSystem.RightHandSideContribution, EquationId);
		
		for( unsigned int i=0; i<rLHS_GlobalElementComponents.size(); i++ )
		  {	    
		    AssembleLHS(rLHS_GlobalElementComponents[i], rLHS_LocalElementComponents[i], EquationId, lock_array);
		  }
		
		for( unsigned int i=0; i<rRHS_GlobalElementComponents.size(); i++ )
		  {
		    AssembleRHS(rRHS_GlobalElementComponents[i], rRHS_LocalElementComponents[i], EquationId, lock_array);
		  }


                // clean local elemental memory
                pScheme->CleanMemory(*it);


            }
        }

        vector<unsigned int> condition_partition;
        CreatePartition(number_of_threads, ConditionsArray.size(), condition_partition);

        #pragma omp parallel for
        for (int k = 0; k < number_of_threads; k++)
        {

	    //contributions to the condition local system
	    ComponentBasedBossakScheme::SystemContributions ConditionLocalSystem;

	    LocalSystemMatrixType LeftHandSideConditionContribution  = LocalSystemMatrixType(0, 0);
	    LocalSystemVectorType RightHandSideConditionContribution = LocalSystemVectorType(0);
	    
	    ConditionLocalSystem.SetLeftHandSideContribution( LeftHandSideConditionContribution );
	    ConditionLocalSystem.SetRightHandSideContribution( RightHandSideConditionContribution );

	    std::vector<LocalSystemMatrixType> rLHS_LocalConditionComponents; 
	    ConditionLocalSystem.SetLHS_Components(rLHS_LocalConditionComponents);
	    
	    std::vector<LocalSystemVectorType> rRHS_LocalConditionComponents; 
	    ConditionLocalSystem.SetRHS_Components(rRHS_LocalConditionComponents);
	    
	    ConditionLocalSystem.LeftHandSideContribution  = LocalSystemMatrixType(0, 0);
	    ConditionLocalSystem.RightHandSideContribution = LocalSystemVectorType(0); 
	    
	    //contributions to the global system
	    std::vector<TSystemMatrixType>& rLHS_GlobalConditionComponents = rGlobalSystem.GetLHS_Condition_Components();
	    ConditionLocalSystem.SetLHS_Condition_Variables( rGlobalSystem.GetLHS_Condition_Variables() );
	    
	    if( rGlobalSystem.GetLHS_Condition_Variables().size() != rLHS_GlobalConditionComponents.size() )
	      rLHS_GlobalConditionComponents.resize(rGlobalSystem.GetLHS_Condition_Variables().size());
	    
	    for( unsigned int i=0; i<rLHS_GlobalConditionComponents.size(); i++ )
	      {
		rLHS_GlobalConditionComponents[i] = LocalSystemMatrixType(0, 0);
	      }
	    
	    std::vector<TSystemVectorType>& rRHS_GlobalConditionComponents = rGlobalSystem.GetRHS_Condition_Components();
	    ConditionLocalSystem.SetRHS_Condition_Variables( rGlobalSystem.GetRHS_Condition_Variables() );
	    
	    if( rGlobalSystem.GetRHS_Condition_Variables().size() != rRHS_GlobalConditionComponents.size() )
	      rRHS_GlobalConditionComponents.resize(rGlobalSystem.GetRHS_Condition_Variables().size());
	    
	    for( unsigned int i=0; i<rRHS_GlobalConditionComponents.size(); i++ )
	      {
		rRHS_GlobalConditionComponents[i] = LocalSystemVectorType(0);
	      }

            ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();
            typename ConditionsArrayType::ptr_iterator it_begin = ConditionsArray.ptr_begin() + condition_partition[k];
            typename ConditionsArrayType::ptr_iterator it_end = ConditionsArray.ptr_begin() + condition_partition[k + 1];

            // assemble all elements
            for (typename ConditionsArrayType::ptr_iterator it = it_begin; it != it_end; ++it)
            {
                //calculate condition contribution
                pScheme->Condition_CalculateSystemContributions(*it, ConditionLocalSytem, EquationId, CurrentProcessInfo);

                //assemble the condition contribution
                Assemble(rGlobalLeftHandSideContribution, rGlobalRightHandSideContribution, 
			 LeftHandSideConditionContribution, RightHandSideConditionContribution, EquationId, lock_array);

 		//assemble the elemental contribution
		//AssembleLHS(rGlobalSystem.LeftHandSideContribution, ConditionLocalSystem.LeftHandSideContribution, EquationId);
		//AssembleRHS(rGlobalSystem.RightHandSideContribution, ConditionLocalSystem.RightHandSideContribution, EquationId);
		
		for( unsigned int i=0; i<rLHS_GlobalConditionComponents.size(); i++ )
		  {	    
		    AssembleLHS(rLHS_GlobalConditionComponents[i], rLHS_LocalConditionComponents[i], EquationId, lock_array);
		  }
		
		for( unsigned int i=0; i<rRHS_GlobalConditionComponents.size(); i++ )
		  {
		    AssembleRHS(rRHS_GlobalConditionComponents[i], rRHS_LocalConditionComponents[i], EquationId, lock_array);
		  }

            }
        }

	double stop_prod = omp_get_wtime();
        if (this->GetEchoLevel() > 2 && r_model_part.GetCommunicator().MyPID() == 0)
            std::cout << "time: " << stop_prod - start_prod << std::endl;

        for (int i = 0; i < A_size; i++)
            omp_destroy_lock(&lock_array[i]);
        if( this->GetEchoLevel() > 2 && r_model_part.GetCommunicator().MyPID() == 0)
        {
            KRATOS_WATCH("finished parallel building");
        }

        //to ensure that all the threads are syncronized here
        // #pragma omp barrier
#endif

	//recovering the reactions flag
	BaseType::mCalculateReactionsFlag = CalculateReactionsFlag;

        KRATOS_CATCH("")

    }

    //**************************************************************************
    //**************************************************************************

    void BuildLHS(
        typename TSchemeType::Pointer pScheme,
        ModelPart& r_model_part,
	GlobalSystemContributions& rGlobalSystem)
    {
        KRATOS_TRY

        //getting the elements from the model
        ElementsArrayType& pElements = r_model_part.Elements();

        //getting the array of the conditions
        ConditionsArrayType& ConditionsArray = r_model_part.Conditions();

        //resetting to zero the vector of reactions
        TSparseSpace::SetToZero(*(BaseType::mpReactionsVector));

        //contributions to the element local system
	ComponentBasedBossakScheme::LHS_SystemContributions ElementLocalLHS;
	
        LocalSystemMatrixType LeftHandSideElementContribution  = LocalSystemMatrixType(0, 0);
	ElementLocalLHS.SetLeftHandSideContribution( LeftHandSideElementContribution );

	std::vector<LocalSystemMatrixType> rLHS_LocalElementComponents; 
	ElementLocalLHS.SetLHS_Components(rLHS_LocalElementComponents);

        //contributions to the global system
	std::vector<TSystemMatrixType>& rLHS_GlobalElementComponents = rGlobalSystem.GetLHS_Element_Components();
	ElementLocalLHS.SetLHS_Element_Variables( rGlobalSystem.GetLHS_Element_Variables() );

	if( rGlobalSystem.GetLHS_Element_Variables().size() != rLHS_GlobalElementComponents.size() )
	  rLHS_GlobalElementComponents.resize(rGlobalSystem.GetLHS_Element_Variables().size());

	for( unsigned int i=0; i<rLHS_GlobalElementComponents.size(); i++ )
	  {
	    rLHS_GlobalElementComponents[i] = LocalSystemMatrixType(0, 0);
	  }


	TSystemMatrixType& rGlobalLeftHandSideContribution  = rGlobalSystem.GetLeftHandSideContribution();

        //vector containing the localization in the system of the different
        //terms
        Element::EquationIdVectorType EquationId;

        ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

        // assemble all elements
        for (typename ElementsArrayType::ptr_iterator it = pElements.ptr_begin(); it != pElements.ptr_end(); ++it)
        {
            //calculate elemental contribution
            pScheme->Calculate_LHS_Contribution(*it, ElementLocalLHS, EquationId, CurrentProcessInfo);

            //assemble the elemental contribution
            AssembleLHS(rGlobalLeftHandSideContribution, LeftHandSideElementContribution, EquationId);

	    for( unsigned int i=0; i<rLHS_GlobalElementComponents.size(); i++ )
	      {	    
		AssembleLHS(rLHS_GlobalElementComponents[i], rLHS_LocalElementComponents[i], EquationId);
	      }

            // clean local elemental memory
            pScheme->CleanMemory(*it);
        }


        //contributions to the condition local system
	ComponentBasedBossakScheme::LHS_SystemContributions ConditionLocalLHS;

        LocalSystemMatrixType LeftHandSideConditionContribution  = LocalSystemMatrixType(0, 0);
	ConditionLocalLHS.SetLeftHandSideContribution( LeftHandSideConditionContribution );

	std::vector<LocalSystemMatrixType> rLHS_LocalConditionComponents; 
	ConditionLocalRHS.SetLHS_Components(rLHS_LocalConditionComponents);

        //contributions to the global system
	std::vector<TSystemMatrixType>& rLHS_GlobalConditionComponents = rGlobalSystem.GetLHS_Condition_Components();
	ConditionLocalLHS.SetLHS_Condition_Variables( rGlobalSystem.GetLHS_Condition_Variables() );

	if( rGlobalSystem.GetLHS_Condition_Variables().size() != rLHS_GlobalConditionComponents.size() )
	  rLHS_GlobalConditionComponents.resize(rGlobalSystem.GetLHS_Condition_Variables().size());

	for( unsigned int i=0; i<rLHS_GlobalConditionComponents.size(); i++ )
	  {
	    rLHS_GlobalConditionComponents[i] = LocalSystemMatrixType(0, 0);
	  }

        // assemble all conditions
        for (typename ConditionsArrayType::ptr_iterator it = ConditionsArray.ptr_begin(); it != ConditionsArray.ptr_end(); ++it)
        {
            //calculate elemental contribution
            pScheme->Condition_Calculate_LHS_Contribution(*it, ConditionalLocalLHS, EquationId, CurrentProcessInfo);

            //assemble the elemental contribution
            AssembleLHS(rGlobalLeftHandSideContribution, LeftHandSideConditionContribution, EquationId);

	    for( unsigned int i=0; i<rLHS_GlobalConditionComponents.size(); i++ )
	      {	    
		AssembleLHS(rLHS_GlobalConditionComponents[i], rLHS_LocalConditionComponents[i], EquationId);
	      }

        }

        KRATOS_CATCH("")

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

        ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

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
            pScheme->Calculate_LHS_Contribution(*it, LHS_Contribution, EquationId, CurrentProcessInfo);

            //assemble the elemental contribution
            AssembleLHS_CompleteOnFreeRows(A, LHS_Contribution, EquationId);

            //clean local elemental memory
            pScheme->CleanMemory(*it);
        }

        LHS_Contribution.resize(0, 0, false);
        // assemble all conditions
        for (typename ConditionsArrayType::ptr_iterator it = ConditionsArray.ptr_begin(); it != ConditionsArray.ptr_end(); ++it)
        {
            //calculate elemental contribution
            pScheme->Condition_Calculate_LHS_Contribution(*it, LHS_Contribution, EquationId, CurrentProcessInfo);

            //assemble the elemental contribution
            AssembleLHS_CompleteOnFreeRows(A, LHS_Contribution, EquationId);
        }


        KRATOS_CATCH("")

    }


    //**************************************************************************
    //**************************************************************************

    void BuildAndSolve(
        typename TSchemeType::Pointer pScheme,
        ModelPart& r_model_part,
	GlobalSystemContributions& rGlobalSystem,
	TSystemVectorType& rDx)
    {
        KRATOS_TRY

        Timer::Start("Build");

        Build(pScheme, r_model_part, rGlobalSystem);

        Timer::Stop("Build");

	TSystemMatrixType& rA = rGlobalSystem.LeftHandSideSystemContribution;
	TSystemVectorType& rb = rGlobalSystem.RightHandSideSystemContribution;


        //does nothing...dirichlet conditions are naturally dealt with in defining the residual
        ApplyDirichletConditions(pScheme, r_model_part, rA, rDx, rb);

        if (this->GetEchoLevel() == 3)
        {
            std::cout << "Before solving the system :" << std::endl;
            std::cout << "System Matrix   = " << rA << std::endl;
            std::cout << "Unknowns vector = " << rDx << std::endl;
            std::cout << "RHS vector      = " << rb << std::endl;
        }

        Timer::Start("Solve");

        SystemSolveWithPhysics(rA, rDx, rb, r_model_part);

        Timer::Stop("Solve");

        if (this->GetEchoLevel() == 3)
        {
            std::cout << "After solving the system:" << std::endl;
            std::cout << "System Matrix   = " << rA << std::endl;
            std::cout << "Unknowns vector = " << rDx << std::endl;
            std::cout << "RHS vector      = " << rb << std::endl;
        }

        KRATOS_CATCH("")
    }

    //**************************************************************************
    //**************************************************************************

    void BuildRHSAndSolve(
        typename TSchemeType::Pointer pScheme,
        ModelPart& r_model_part,
	GlobalSystemContributions& rGlobalSystem,
        TSystemVectorType& rDx)
    {
        KRATOS_TRY

	TSystemMatrixType& rA = rGlobalSystem.LeftHandSideSystemContribution;
	TSystemVectorType& rb = rGlobalSystem.RightHandSideSystemContribution;

        bool CalculateReactionsFlag = BaseType::mCalculateReactionsFlag;
	
	BaseType::mCalculateReactionsFlag = false;
        BuildRHS(pScheme, r_model_part, rGlobalSystem);
	BaseType::mCalculateReactionsFlag = CalculateReactionsFlag;

        SystemSolve(rA, rDx, rb);

        KRATOS_CATCH("")
    }

    //**************************************************************************
    //**************************************************************************

    void BuildRHS(
        typename TSchemeType::Pointer pScheme,
        ModelPart& r_model_part,
	GlobalSystemContributions& rGlobalSystem)
    {
        KRATOS_TRY

        //Getting the Elements
        ElementsArrayType& pElements = r_model_part.Elements();

        //getting the array of the conditions
        ConditionsArrayType& ConditionsArray = r_model_part.Conditions();

        ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

        //resetting to zero the vector of reactions
        TSparseSpace::SetToZero(*(BaseType::mpReactionsVector));

        //contributions to the system
	ComponentBasedBossakScheme::RHS_SystemContributions ElementLocalRHS;
	
	LocalSystemVectorType RightHandSideElementContribution = LocalSystemVectorType(0);
	ElementLocalRHS.SetRightHandSideContribution( RightHandSideElementContribution );
	
	std::vector<LocalSystemVectorType> rRHS_LocalElementComponents; 
	ElementLocalRHS.SetRHS_Components(rRHS_LocalElementComponents);
	
	    
	//contributions to the global system   
	std::vector<TSystemVectorType>& rRHS_GlobalElementComponents = rGlobalSystem.GetRHS_Element_Components();
	ElementLocalRHS.SetRHS_Element_Variables( rGlobalSystem.GetRHS_Element_Variables() );
	
	if( rGlobalSystem.GetRHS_Element_Variables().size() != rRHS_GlobalElementComponents.size() )
	  rRHS_GlobalElementComponents.resize(rGlobalSystem.GetRHS_Element_Variables().size());
	
	for( unsigned int i=0; i<rRHS_GlobalElementComponents.size(); i++ )
	  {
	    rRHS_GlobalElementComponents[i] = LocalSystemVectorType(0);
	  }
	    

	TSystemVectorType& rGlobalRightHandSideContribution = rGlobalSystem.GetRightHandSideContribution();

        LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0, 0);
        LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

        //vector containing the localization in the system of the different
        //terms
        Element::EquationIdVectorType EquationId;

        // assemble all elements
        for (typename ElementsArrayType::ptr_iterator it = pElements.ptr_begin(); it != pElements.ptr_end(); ++it)
        {
            //calculate elemental Right Hand Side Contribution
            pScheme->Calculate_RHS_Contribution(*it, ElementLocalRHS, EquationId, CurrentProcessInfo);

            //assemble the elemental contribution
            AssembleRHS(rGlobalRightHandSideContribution, RightHandSideElementContribution, EquationId);

	    for( unsigned int i=0; i<rRHS_GlobalElementComponents.size(); i++ )
	      {
		AssembleRHS(rRHS_GlobalElementComponents[i], rRHS_LocalElementComponents[i], EquationId);
	      }
        }

        //contributions to the condition local system
	ComponentBasedBossakScheme::RHS_SystemContributions ConditionLocalRHS;

        LocalSystemVectorType RightHandSideConditionContribution = LocalSystemVectorType(0);
	ConditionLocalRHS.SetRightHandSideContribution( RightHandSideConditionContribution );

	std::vector<LocalSystemVectorType> rRHS_LocalConditionComponents; 
	ConditionLocaRHS.SetRHS_Components(rRHS_LocalConditionComponents);

        //contributions to the global system
	std::vector<TSystemVectorType>& rRHS_GlobalConditionComponents = rGlobalSystem.GetRHS_Condition_Components();
	ConditionLocalRHS.SetRHS_Condition_Variables( rGlobalSystem.GetRHS_Condition_Variables() );

	if( rGlobalSystem.GetRHS_Condition_Variables().size() != rRHS_GlobalConditionComponents.size() )
	  rRHS_GlobalConditionComponents.resize(rGlobalSystem.GetRHS_Condition_Variables().size());

	for( unsigned int i=0; i<rRHS_GlobalConditionComponents.size(); i++ )
	  {
	    rRHS_GlobalConditionComponents[i] = LocalSystemVectorType(0);
	  }


        // assemble all conditions
        for (typename ConditionsArrayType::ptr_iterator it = ConditionsArray.ptr_begin(); it != ConditionsArray.ptr_end(); ++it)
        {
            //calculate elemental contribution
            pScheme->Condition_Calculate_RHS_Contribution(*it, ConditionLocalRHS, EquationId, CurrentProcessInfo);

            //assemble the elemental contribution
            AssembleRHS(rGlobalRightHandSideContribution, RightHandSideConditionContribution, EquationId);

	    for( unsigned int i=0; i<rRHS_GlobalConditionComponents.size(); i++ )
	      {
		AssembleRHS(rRHS_GlobalConditionComponents[i], rRHS_LocalConditionComponents[i], EquationId);
	      }
        }

        KRATOS_CATCH("")

    }


    //**************************************************************************
    //**************************************************************************

    void ResizeAndInitializeVectors(
        TSystemMatrixPointerType& pA,
        TSystemVectorPointerType& pDx,
        TSystemVectorPointerType& pb,
        ElementsArrayType& rElements,
        ConditionsArrayType& rConditions,
        ProcessInfo& CurrentProcessInfo
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
            ConstructMatrixStructure(A, rElements, rConditions, CurrentProcessInfo);
        }
        else
        {
            if (A.size1() != BaseType::mEquationSystemSize || A.size2() != BaseType::mEquationSystemSize)
            {
                KRATOS_WATCH("it should not come here!!!!!!!! ... this is SLOW");
                A.resize(BaseType::mEquationSystemSize, BaseType::mEquationSystemSize, true);
                ConstructMatrixStructure(A, rElements, rConditions, CurrentProcessInfo);
            }
        }
        if (Dx.size() != BaseType::mEquationSystemSize)
            Dx.resize(BaseType::mEquationSystemSize, false);
        if (b.size() != BaseType::mEquationSystemSize)
            b.resize(BaseType::mEquationSystemSize, false);

 
        //if needed resize the vector for the calculation of reactions
        if (BaseType::mCalculateReactionsFlag == true)
        {
            unsigned int ReactionsVectorSize = BaseType::mDofSet.size() - BaseType::mEquationSystemSize;
            if (BaseType::mpReactionsVector->size() != ReactionsVectorSize)
                BaseType::mpReactionsVector->resize(ReactionsVectorSize, false);
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
        KRATOS_CATCH("")
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
        KRATOS_CATCH("")
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
        KRATOS_CATCH("");
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

}; /* Class ComponentBasedBuilderAndSolver */

/*@} */

/**@name Type Definitions */
/*@{ */
/*@} */

} /* namespace Kratos.*/

#endif /* KRATOS_COMPONENT_BASED_BUILDER_AND_SOLVER  defined */

