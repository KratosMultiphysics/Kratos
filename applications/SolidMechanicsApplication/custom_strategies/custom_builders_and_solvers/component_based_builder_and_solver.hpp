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
      std::vector<TSystemMatrixType> *mpLHS_Element_Contributions;
      const std::vector< Variable< TSystemMatrixType > > *mpLHS_Element_Variables;

      std::vector<TSystemVectorType> *mpRHS_Element_Contributions;
      const std::vector< Variable< TSystemVectorType > > *mpRHS_Element_Variables;
      
      //conditions
      std::vector<TSystemMatrixType> *mpLHS_Condition_Contributions;
      const std::vector< Variable< TSystemMatrixType > > *mpLHS_Condition_Variables;

      std::vector<TSystemVectorType> *mpRHS_Condition_Contributions;
      const std::vector< Variable< TSystemVectorType > > *mpRHS_Condition_Variables;
      
    public:
      
      //setting pointer variables
      void SetLeftHandSideContribution  ( TSystemMatrixType& rLeftHandSideContribution )  { mpLeftHandSideContribution = &rLeftHandSideContribution; };
      void SetRightHandSideContribution ( TSystemVectorType& rRightHandSideContribution ) { mpRightHandSideContribution = &rRightHandSideContribution; };

      //elements
      void SetLHS_Element_Contributions ( std::vector<TSystemMatrixType>& rLHS_Element_Contributions ) { mpLHS_Element_Contributions = &rLHS_Element_Contributions; };
      void SetLHS_Element_Variables     ( std::vector< Variable< LocalSystemMatrixType > >& rLHS_Element_Variables ) { mpLHS_Element_Variables = &rLHS_Element_Variables; };
      void SetRHS_Element_Contributions ( std::vector<TSystemVectorType>& rRHS_Element_Contributions ) { mpRHS_Element_Contributions = &rRHS_Element_Contributions; };
      void SetRHS_Element_Variables     ( std::vector< Variable< LocalSystemVectorType > >& rRHS_Element_Variables ) { mpRHS_Element_Variables = &rRHS_Element_Variables; };

      //conditions
      void SetLHS_Condition_Contributions ( std::vector<TSystemMatrixType>& rLHS_Condition_Contributions ) { mpLHS_Condition_Contributions = &rLHS_Condition_Contributions; };
      void SetLHS_Condition_Variables     ( std::vector< Variable< LocalSystemMatrixType > >& rLHS_Condition_Variables ) { mpLHS_Condition_Variables = &rLHS_Condition_Variables; };
      void SetRHS_Condition_Contributions ( std::vector<TSystemVectorType>& rRHS_Condition_Contributions ) { mpRHS_Condition_Contributions = &rRHS_Condition_Contributions; };
      void SetRHS_Condition_Variables     ( std::vector< Variable< LocalSystemVectorType > >& rRHS_Condition_Variables ) { mpRHS_Condition_Variables = &rRHS_Condition_Variables; };

      //getting pointer variables
      TSystemMatrixType& GetLeftHandSideContribution  () { return *mpLeftHandSideContribution; };
      TSystemVectorType& GetRightHandSideContribution () { return *mpRightHandSideContribution; };

      //elements
      std::vector<TSystemMatrixType>& GetLHS_Element_Contributions() { return *mpLHS_Element_Contributions; };
      std::vector< Variable< LocalSystemMatrixType > >& GetLHS_Element_Variables() { return *mpLHS_Element_Variables; };
      std::vector<TSystemVectorType>& GetRHS_Element_Contributions() { return *mpRHS_Element_Contributions; };
      std::vector< Variable< LocalSystemVectorType > >& GetRHS_Element_Variables() { return *mpRHS_Element_Variables; };

      //conditions
      std::vector<TSystemMatrixType>& GetLHS_Condition_Contributions() { return *mpLHS_Condition_Contributions; };
      std::vector< Variable< LocalSystemMatrixType > >& GetLHS_Condition_Variables() { return *mpLHS_Condition_Variables; };
      std::vector<TSystemVectorType>& GetRHS_Condition_Contributions() { return *mpRHS_Condition_Contributions; };
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
	GlobalSystemMatrixType& rGlobalSystem)
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

        //contributions to the element local system
	ComponentBasedBossakScheme::SystemContributions ElementLocalSystem;
	
	std::vector<LocalSystemMatrixType> rLHS_LocalElementContributions; 
	ElementLocalSystem.SetLHS_Contributions(rLHS_LocalElementContributions);

	std::vector<LocalSystemVectorType> rRHS_LocalElementContributions; 
	ElementLocalSystem.SetRHS_Contributions(rRHS_LocalElementContributions);

        ElementLocalSystem.LeftHandSideContribution  = LocalSystemMatrixType(0, 0);
        ElementLocalSystem.RightHandSideContribution = LocalSystemVectorType(0);

        //contributions to the global system
	std::vector<TSystemMatrixType>& rLHS_GlobalElementContributions = rGlobalSystem.GetLHS_Element_Contributions();
	ElementLocalSystem.SetLHS_Element_Variables( rGlobalSystem.GetLHS_Element_Variables() );

	if( rGlobalSystem.GetLHS_Element_Variables().size() != rLHS_GlobalElementContributions.size() )
	  rLHS_GlobalElementContributions.resize(rGlobalSystem.GetLHS_Element_Variables().size());

	for( unsigned int i=0; i<rLHS_GlobalElementContributions.size(); i++ )
	  {
	    rLHS_GlobalElementContributions[i] = LocalSystemMatrixType(0, 0);
	  }

	std::vector<TSystemVectorType>& rRHS_GlobalElementContributions = rGlobalSystem.GetRHS_Element_Contributions();
	ElementLocalSystem.SetRHS_Element_Variables( rGlobalSystem.GetRHS_Element_Variables() );

	if( rGlobalSystem.GetRHS_Element_Variables().size() != rRHS_GlobalElementContributions.size() )
	  rRHS_GlobalElementContributions.resize(rGlobalSystem.GetRHS_Element_Variables().size());

	for( unsigned int i=0; i<rRHS_GlobalElementContributions.size(); i++ )
	  {
	    rRHS_GlobalElementContributions[i] = LocalSystemVectorType(0);
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
            AssembleLHS(rGlobalLeftHandSideContribution, ElementLocalSystem.LeftHandSideContribution, EquationId);
            AssembleRHS(rGlobalRightHandSideContribution, ElementLocalSystem.RightHandSideContribution, EquationId);

	    for( unsigned int i=0; i<rLHS_GlobalContributions.size(); i++ )
	      {	    
		AssembleLHS(rLHS_GlobalElementElementContributions[i], rLHS_LocalElementContributions[i], EquationId);
	      }

	    for( unsigned int i=0; i<rRHS_GlobalElementContributions.size(); i++ )
	      {
		AssembleRHS(rRHS_GlobalElementContributions[i], rRHS_LocalElementContributions[i], EquationId);
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

	std::vector<LocalSystemMatrixType> rLHS_LocalConditionContributions; 
	ConditionLocalSystem.SetLHS_Contributions(rLHS_LocalConditionContributions);

	std::vector<LocalSystemVectorType> rRHS_LocalConditionContributions; 
	ConditionLocalSystem.SetRHS_Contributions(rRHS_LocalConditionContributions);

	ConditionLocalSystem.LeftHandSideContribution  = LocalSystemMatrixType(0, 0);
        ConditionLocalSystem.RightHandSideContribution = LocalSystemVectorType(0); 

        //contributions to the global system
	std::vector<TSystemMatrixType>& rLHS_GlobalConditionContributions = rGlobalSystem.GetLHS_Condition_Contributions();
	ConditionLocalSystem.SetLHS_Condition_Variables( rGlobalSystem.GetLHS_Condition_Variables() );

	if( rGlobalSystem.GetLHS_Condition_Variables().size() != rLHS_GlobalConditionContributions.size() )
	  rLHS_GlobalConditionContributions.resize(rGlobalSystem.GetLHS_Condition_Variables().size());

	for( unsigned int i=0; i<rLHS_GlobalConditionContributions.size(); i++ )
	  {
	    rLHS_GlobalConditionContributions[i] = LocalSystemMatrixType(0, 0);
	  }

	std::vector<TSystemVectorType>& rRHS_GlobalConditionContributions = rGlobalSystem.GetRHS_Condition_Contributions();
	ConditionLocalSystem.SetRHS_Condition_Variables( rGlobalSystem.GetRHS_Condition_Variables() );

	if( rGlobalSystem.GetRHS_Condition_Variables().size() != rRHS_GlobalConditionContributions.size() )
	  rRHS_GlobalConditionContributions.resize(rGlobalSystem.GetRHS_Condition_Variables().size());

	for( unsigned int i=0; i<rRHS_GlobalConditionContributions.size(); i++ )
	  {
	    rRHS_GlobalConditionContributions[i] = LocalSystemVectorType(0);
	  }


        // assemble all conditions
        for (typename ConditionsArrayType::ptr_iterator it = ConditionsArray.ptr_begin(); it != ConditionsArray.ptr_end(); ++it)
        {
            //calculate condition contribution
  	    pScheme->Condition_CalculateSystemContributions(*it, ConditionLocalSytem, EquationId, CurrentProcessInfo);

            //assemble the condition contribution
            AssembleLHS(rGlobalLeftHandSideContribution, ConditionLocalSystem.LeftHandSideContribution, EquationId);
            AssembleRHS(rGlobalRightHandSideContribution, ConditionLocalSystem.RightHandSideContribution, EquationId);

	    for( unsigned int i=0; i<rLHS_GlobalContributions.size(); i++ )
	      {	    
		AssembleLHS(rLHS_GlobalConditionContributions[i], rLHS_LocalConditionContributions[i], EquationId);
	      }

	    for( unsigned int i=0; i<rRHS_GlobalContributions.size(); i++ )
	      {
		AssembleRHS(rRHS_GlobalConditionContributions[i], rRHS_LocalConditionContributions[i], EquationId);
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
	  
	    std::vector<LocalSystemMatrixType> rLHS_LocalElementContributions; 
	    ElementLocalSystem.SetLHS_Contributions(rLHS_LocalElementContributions);
	    
	    std::vector<LocalSystemVectorType> rRHS_LocalElementContributions; 
	    ElementLocalSystem.SetRHS_Contributions(rRHS_LocalElementContributions);
	    
	    ElementLocalSystem.LeftHandSideContribution  = LocalSystemMatrixType(0, 0);
	    ElementLocalSystem.RightHandSideContribution = LocalSystemVectorType(0);
	    
	    //contributions to the global system
	    std::vector<TSystemMatrixType>& rLHS_GlobalElementContributions = rGlobalSystem.GetLHS_Element_Contributions();
	    ElementLocalSystem.SetLHS_Element_Variables( rGlobalSystem.GetLHS_Element_Variables() );
	    
	    if( rGlobalSystem.GetLHS_Element_Variables().size() != rLHS_GlobalElementContributions.size() )
	      rLHS_GlobalElementContributions.resize(rGlobalSystem.GetLHS_Element_Variables().size());
	    
	    for( unsigned int i=0; i<rLHS_GlobalElementContributions.size(); i++ )
	      {
		rLHS_GlobalElementContributions[i] = LocalSystemMatrixType(0, 0);
	      }
	    
	    std::vector<TSystemVectorType>& rRHS_GlobalElementContributions = rGlobalSystem.GetRHS_Element_Contributions();
	    ElementLocalSystem.SetRHS_Element_Variables( rGlobalSystem.GetRHS_Element_Variables() );
	    
	    if( rGlobalSystem.GetRHS_Element_Variables().size() != rRHS_GlobalElementContributions.size() )
	      rRHS_GlobalElementContributions.resize(rGlobalSystem.GetRHS_Element_Variables().size());
	    
	    for( unsigned int i=0; i<rRHS_GlobalElementContributions.size(); i++ )
	      {
		rRHS_GlobalElementContributions[i] = LocalSystemVectorType(0);
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
                Assemble(rGlobalSystem.LeftHandSideContribution, rGlobalSystem.RightHandSideContribution, 
			 ElementLocalSystem.LeftHandSideContribution, ElementLocalSystem.RightHandSideContribution, EquationId, lock_array);

		//assemble the elemental contribution
		//AssembleLHS(rGlobalSystem.LeftHandSideContribution, ElementLocalSystem.LeftHandSideContribution, EquationId);
		//AssembleRHS(rGlobalSystem.RightHandSideContribution, ElementLocalSystem.RightHandSideContribution, EquationId);
		
		for( unsigned int i=0; i<rLHS_GlobalContributions.size(); i++ )
		  {	    
		    AssembleLHS(rLHS_GlobalElementElementContributions[i], rLHS_LocalElementContributions[i], EquationId, lock_array);
		  }
		
		for( unsigned int i=0; i<rRHS_GlobalElementContributions.size(); i++ )
		  {
		    AssembleRHS(rRHS_GlobalElementContributions[i], rRHS_LocalElementContributions[i], EquationId, lock_array);
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

	    std::vector<LocalSystemMatrixType> rLHS_LocalConditionContributions; 
	    ConditionLocalSystem.SetLHS_Contributions(rLHS_LocalConditionContributions);
	    
	    std::vector<LocalSystemVectorType> rRHS_LocalConditionContributions; 
	    ConditionLocalSystem.SetRHS_Contributions(rRHS_LocalConditionContributions);
	    
	    ConditionLocalSystem.LeftHandSideContribution  = LocalSystemMatrixType(0, 0);
	    ConditionLocalSystem.RightHandSideContribution = LocalSystemVectorType(0); 
	    
	    //contributions to the global system
	    std::vector<TSystemMatrixType>& rLHS_GlobalConditionContributions = rGlobalSystem.GetLHS_Condition_Contributions();
	    ConditionLocalSystem.SetLHS_Condition_Variables( rGlobalSystem.GetLHS_Condition_Variables() );
	    
	    if( rGlobalSystem.GetLHS_Condition_Variables().size() != rLHS_GlobalConditionContributions.size() )
	      rLHS_GlobalConditionContributions.resize(rGlobalSystem.GetLHS_Condition_Variables().size());
	    
	    for( unsigned int i=0; i<rLHS_GlobalConditionContributions.size(); i++ )
	      {
		rLHS_GlobalConditionContributions[i] = LocalSystemMatrixType(0, 0);
	      }
	    
	    std::vector<TSystemVectorType>& rRHS_GlobalConditionContributions = rGlobalSystem.GetRHS_Condition_Contributions();
	    ConditionLocalSystem.SetRHS_Condition_Variables( rGlobalSystem.GetRHS_Condition_Variables() );
	    
	    if( rGlobalSystem.GetRHS_Condition_Variables().size() != rRHS_GlobalConditionContributions.size() )
	      rRHS_GlobalConditionContributions.resize(rGlobalSystem.GetRHS_Condition_Variables().size());
	    
	    for( unsigned int i=0; i<rRHS_GlobalConditionContributions.size(); i++ )
	      {
		rRHS_GlobalConditionContributions[i] = LocalSystemVectorType(0);
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
                Assemble(rGlobalSystem.LeftHandSideContribution, rGlobalSystem.RightHandSideContribution, 
			 ConditionLocalSystem.LeftHandSideContribution, ConditionLocalSystem.RightHandSideContribution, EquationId, lock_array);

 		//assemble the elemental contribution
		//AssembleLHS(rGlobalSystem.LeftHandSideContribution, ConditionLocalSystem.LeftHandSideContribution, EquationId);
		//AssembleRHS(rGlobalSystem.RightHandSideContribution, ConditionLocalSystem.RightHandSideContribution, EquationId);
		
		for( unsigned int i=0; i<rLHS_GlobalContributions.size(); i++ )
		  {	    
		    AssembleLHS(rLHS_GlobalConditionContributions[i], rLHS_LocalConditionContributions[i], EquationId, lock_array);
		  }
		
		for( unsigned int i=0; i<rRHS_GlobalContributions.size(); i++ )
		  {
		    AssembleRHS(rRHS_GlobalConditionContributions[i], rRHS_LocalConditionContributions[i], EquationId, lock_array);
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


        KRATOS_CATCH("")

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

        //contributions to the system
        LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0, 0);

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
            AssembleLHS(A, LHS_Contribution, EquationId);

            // clean local elemental memory
            pScheme->CleanMemory(*it);
        }

        LHS_Contribution.resize(0, 0, false);

        // assemble all conditions
        for (typename ConditionsArrayType::ptr_iterator it = ConditionsArray.ptr_begin(); it != ConditionsArray.ptr_end(); ++it)
        {
            //calculate elemental contribution
            pScheme->Condition_Calculate_LHS_Contribution(*it, LHS_Contribution, EquationId, CurrentProcessInfo);

            //assemble the elemental contribution
            AssembleLHS(A, LHS_Contribution, EquationId);
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

            // clean local elemental memory
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
	GlobalSystemMatrixType& rGlobalSystem,
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
	GlobalSystemMatrixType& rGlobalSystem,
        TSystemVectorType& rDx)
    {
        KRATOS_TRY

	TSystemMatrixType& rA = rGlobalSystem.LeftHandSideSystemContribution;
	TSystemVectorType& rb = rGlobalSystem.RightHandSideSystemContribution;

        BuildRHS(pScheme, r_model_part, rb);
        SystemSolve(rA, rDx, rb);

        KRATOS_CATCH("")
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

        //contributions to the system
        LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0, 0);
        LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

        //vector containing the localization in the system of the different
        //terms
        Element::EquationIdVectorType EquationId;

        // assemble all elements
        for (typename ElementsArrayType::ptr_iterator it = pElements.ptr_begin(); it != pElements.ptr_end(); ++it)
        {
            //calculate elemental Right Hand Side Contribution
            pScheme->Calculate_RHS_Contribution(*it, RHS_Contribution, EquationId, CurrentProcessInfo);

            //assemble the elemental contribution
            AssembleRHS(b, RHS_Contribution, EquationId);
        }

        LHS_Contribution.resize(0, 0, false);
        RHS_Contribution.resize(0, false);

        // assemble all conditions
        for (typename ConditionsArrayType::ptr_iterator it = ConditionsArray.ptr_begin(); it != ConditionsArray.ptr_end(); ++it)
        {
            //calculate elemental contribution
            pScheme->Condition_Calculate_RHS_Contribution(*it, RHS_Contribution, EquationId, CurrentProcessInfo);

            //assemble the elemental contribution
            AssembleRHS(b, RHS_Contribution, EquationId);
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
        BuildRHS(pScheme, r_model_part, b);

        int i;
        int systemsize = BaseType::mDofSet.size() - TSparseSpace::Size(*BaseType::mpReactionsVector);

        typename DofsArrayType::ptr_iterator it2;
        // KRATOS_WATCH(*BaseType::mpReactionsVector);
        //updating variables
        TSystemVectorType& ReactionsVector = *BaseType::mpReactionsVector;
        int num=1;
        for (it2 = BaseType::mDofSet.ptr_begin(); it2 != BaseType::mDofSet.ptr_end(); ++it2)
        {
            if ((*it2)->IsFixed())
            {
                i = (*it2)->EquationId();
                i -= systemsize;
                /*KRATOS_WATCH((*it2)->GetSolutionStepReactionValue());
                KRATOS_WATCH(ReactionsVector[i]);*/
                (*it2)->GetSolutionStepReactionValue() = ReactionsVector[i];
            }
            num++;
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

