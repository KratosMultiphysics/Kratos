//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//                    
//

#if !defined(KRATOS_RESIDUAL_BASED_ELIMINATION_BUILDER_AND_SOLVER )
#define  KRATOS_RESIDUAL_BASED_ELIMINATION_BUILDER_AND_SOLVER


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

// #define USE_GOOGLE_HASH
#ifdef USE_GOOGLE_HASH
#include "sparsehash/dense_hash_set" //included in external libraries
#else
#include <unordered_set>
#endif


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

/** Short class definition.

Detail class definition.

Current class provides an implementation for standard builder and solving operations.

the RHS is constituted by the unbalanced loads (residual)

Degrees of freedom are reordered putting the restrained degrees of freedom at
the end of the system ordered in reverse order with respect to the DofSet.

Imposition of the dirichlet conditions is naturally dealt with as the residual already contains
this information.

Calculation of the reactions involves a cost very similiar to the calculation of the total residual

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
         class TDenseSpace, //= DenseSpace<double>,
         class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
         >
class ResidualBasedEliminationBuilderAndSolver
    : public BuilderAndSolver< TSparseSpace, TDenseSpace, TLinearSolver >
{
public:
    /**@name Type Definitions */
    /*@{ */
    KRATOS_CLASS_POINTER_DEFINITION(ResidualBasedEliminationBuilderAndSolver);

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

	struct dof_iterator_hash
	{
		size_t operator()(const Node<3>::DofType::Pointer& it) const
		{
			std::size_t seed = 0;
			boost::hash_combine(seed, it->Id());
			boost::hash_combine(seed, (it->GetVariable()).Key());
			return seed;
		}
	};

	struct dof_iterator_equal
	{
		size_t operator()(const Node<3>::DofType::Pointer& it1, const Node<3>::DofType::Pointer& it2) const
		{
			return (((it1->Id() == it2->Id() && (it1->GetVariable()).Key()) == (it2->GetVariable()).Key()));
		}
	};


    /*@} */
    /**@name Life Cycle
     */
    /*@{ */

    /** Constructor.
     */
    ResidualBasedEliminationBuilderAndSolver(
        typename TLinearSolver::Pointer pNewLinearSystemSolver)
        : BuilderAndSolver< TSparseSpace, TDenseSpace, TLinearSolver >(pNewLinearSystemSolver)
    {

        /* 			std::cout << "using the standard builder and solver " << std::endl; */

    }

    /** Destructor.
     */
    virtual ~ResidualBasedEliminationBuilderAndSolver()
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
        TSystemMatrixType& A,
        TSystemVectorType& b)
    {
        KRATOS_TRY
			if (!pScheme)
				KRATOS_THROW_ERROR(std::runtime_error, "No scheme provided!", "");

		//getting the elements from the model
		const int nelements = static_cast<int>(r_model_part.Elements().size());

		//getting the array of the conditions
		const int nconditions = static_cast<int>(r_model_part.Conditions().size());

		ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();
		ModelPart::ElementsContainerType::iterator el_begin = r_model_part.ElementsBegin();
		ModelPart::ConditionsContainerType::iterator cond_begin = r_model_part.ConditionsBegin();

		//contributions to the system
		LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0, 0);
		LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

		//vector containing the localization in the system of the different
		//terms
		Element::EquationIdVectorType EquationId;

		// assemble all elements
		double start_build = OpenMPUtils::GetCurrentTime();
		
                #pragma omp parallel firstprivate(nelements, nconditions,  LHS_Contribution, RHS_Contribution, EquationId )
                {
                    #pragma omp  for schedule(guided, 512) nowait
                    for (int k = 0; k < nelements; k++)
                    {
                            ModelPart::ElementsContainerType::iterator it = el_begin + k;

                            //detect if the element is active or not. If the user did not make any choice the element
                            //is active by default
                            bool element_is_active = true;
                            if ((it)->IsDefined(ACTIVE))
                                    element_is_active = (it)->Is(ACTIVE);

                            if (element_is_active)
                            {
                                    //calculate elemental contribution
                                    pScheme->CalculateSystemContributions(*(it.base()), LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);

                                    //assemble the elemental contribution
    #ifdef _OPENMP
                                    Assemble(A, b, LHS_Contribution, RHS_Contribution, EquationId, mlock_array);
    #else
                                    Assemble(A, b, LHS_Contribution, RHS_Contribution, EquationId);
    #endif
                                    // clean local elemental memory
                                    pScheme->CleanMemory(*(it.base()));
                                    
                            }

                    }

                    #pragma omp  for schedule(guided, 512)
                    for (int k = 0; k < nconditions; k++)
                    {
                            ModelPart::ConditionsContainerType::iterator it = cond_begin + k;

                            //detect if the element is active or not. If the user did not make any choice the element
                            //is active by default
                            bool condition_is_active = true;
                            if ((it)->IsDefined(ACTIVE))
                                    condition_is_active = (it)->Is(ACTIVE);

                            if (condition_is_active)
                            {
                                    //calculate elemental contribution
                                    pScheme->Condition_CalculateSystemContributions(*(it.base()), LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);

    #ifdef _OPENMP
                                    Assemble(A, b, LHS_Contribution, RHS_Contribution, EquationId, mlock_array);
    #else
                                    Assemble(A, b, LHS_Contribution, RHS_Contribution, EquationId);
    #endif

                                    // clean local elemental memory
                                    pScheme->CleanMemory(*(it.base()));
                            }
                    }
                }

		double stop_build = OpenMPUtils::GetCurrentTime();
		if (this->GetEchoLevel() > 1 && r_model_part.GetCommunicator().MyPID() == 0)
			std::cout << "build time: " << stop_build - start_build << std::endl;

		if (this->GetEchoLevel() > 2 && r_model_part.GetCommunicator().MyPID() == 0)
		{
			KRATOS_WATCH("finished building");
		}


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
        ElementsArrayType& rElements = r_model_part.Elements();

        //getting the array of the conditions
        ConditionsArrayType& rConditions = r_model_part.Conditions();

        //resetting to zero the vector of reactions
        TSparseSpace::SetToZero(*(BaseType::mpReactionsVector));

        //contributions to the system
        LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0, 0);

        //vector containing the localization in the system of the different
        //terms
        Element::EquationIdVectorType EquationId;

        ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

        // assemble all elements
        for (typename ElementsArrayType::ptr_iterator it = rElements.ptr_begin(); it != rElements.ptr_end(); ++it)
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
        for (typename ConditionsArrayType::ptr_iterator it = rConditions.ptr_begin(); it != rConditions.ptr_end(); ++it)
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
        ElementsArrayType& rElements = r_model_part.Elements();

        //getting the array of the conditions
        ConditionsArrayType& rConditions = r_model_part.Conditions();

        ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

        //resetting to zero the vector of reactions
        TSparseSpace::SetToZero(*(BaseType::mpReactionsVector));

        //contributions to the system
        LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0, 0);

        //vector containing the localization in the system of the different
        //terms
        Element::EquationIdVectorType EquationId;

        // assemble all elements
        for (typename ElementsArrayType::ptr_iterator it = rElements.ptr_begin(); it != rElements.ptr_end(); ++it)
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
        for (typename ConditionsArrayType::ptr_iterator it = rConditions.ptr_begin(); it != rConditions.ptr_end(); ++it)
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

        KRATOS_CATCH("")

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
            std::cout << "ATTENTION! setting the RHS to zero!" << std::endl;
        }

        //prints informations about the current time
        if (this->GetEchoLevel() > 1)
        {
            std::cout << *(BaseType::mpLinearSystemSolver) << std::endl;
        }

        KRATOS_CATCH("")

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

        // 			boost::timer building_time;

        Timer::Start("Build");

        Build(pScheme, r_model_part, A, b);

        Timer::Stop("Build");


        // 			if(this->GetEchoLevel()>0)
        // 			{
        // 				std::cout << "Building Time : " << building_time.elapsed() << std::endl;
        // 			}

        //			ApplyPointLoads(pScheme,r_model_part,b);

        //does nothing...dirichlet conditions are naturally dealt with in defining the residual
        ApplyDirichletConditions(pScheme, r_model_part, A, Dx, b);

        if (this->GetEchoLevel() == 3)
        {
            std::cout << "before the solution of the system" << std::endl;
            std::cout << "System Matrix = " << A << std::endl;
            std::cout << "unknowns vector = " << Dx << std::endl;
            std::cout << "RHS vector = " << b << std::endl;
        }

        // 			boost::timer solve_time;
        Timer::Start("Solve");

        SystemSolveWithPhysics(A, Dx, b, r_model_part);

        Timer::Stop("Solve");

        // 			if(this->GetEchoLevel()>0)
        // 			{
        // 				std::cout << "System Solve Time : " << solve_time.elapsed() << std::endl;
        // 			}
        if (this->GetEchoLevel() == 3)
        {
            std::cout << "after the solution of the system" << std::endl;
            std::cout << "System Matrix = " << A << std::endl;
            std::cout << "unknowns vector = " << Dx << std::endl;
            std::cout << "RHS vector = " << b << std::endl;
        }

        KRATOS_CATCH("")
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

        BuildRHS(pScheme, r_model_part, b);
        SystemSolve(A, Dx, b);

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
                    
        //resetting to zero the vector of reactions
        
        if(BaseType::mCalculateReactionsFlag)
        {
            TSparseSpace::SetToZero(*(BaseType::mpReactionsVector));
        }
        

		//Getting the Elements
		ElementsArrayType& pElements = r_model_part.Elements();

		//getting the array of the conditions
		ConditionsArrayType& pConditions = r_model_part.Conditions();

		ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

		//contributions to the system
		LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

		//vector containing the localization in the system of the different terms
		Element::EquationIdVectorType EquationId;

		// assemble all elements
		
                #pragma omp parallel firstprivate( RHS_Contribution, EquationId)
                {
                    const int nelements = static_cast<int>(pElements.size());
                    #pragma omp for schedule(guided, 512) nowait
                    for (int i = 0; i<nelements; i++)
                    {
                            typename ElementsArrayType::iterator it = pElements.begin() + i;
                            //detect if the element is active or not. If the user did not make any choice the element
                            //is active by default
                            bool element_is_active = true;
                            if ((it)->IsDefined(ACTIVE))
                                    element_is_active = (it)->Is(ACTIVE);

                            if (element_is_active)
                            {
                                    //calculate elemental Right Hand Side Contribution
                                    pScheme->Calculate_RHS_Contribution(*(it.base()), RHS_Contribution, EquationId, CurrentProcessInfo);

                                    //assemble the elemental contribution
                                    AssembleRHS(b, RHS_Contribution, EquationId);
                            }
                    }

                    // assemble all conditions
                    const int nconditions = static_cast<int>(pConditions.size());
                    #pragma omp  for schedule(guided, 512)
                    for (int i = 0; i<nconditions; i++)
                    {
                            auto it = pConditions.begin() + i;
                            //detect if the element is active or not. If the user did not make any choice the element
                            //is active by default
                            bool condition_is_active = true;
                            if ((it)->IsDefined(ACTIVE))
                                    condition_is_active = (it)->Is(ACTIVE);

                            if (condition_is_active)
                            {

                                    //calculate elemental contribution
                                    pScheme->Condition_Calculate_RHS_Contribution(*(it.base()), RHS_Contribution, EquationId, CurrentProcessInfo);

                                    //assemble the elemental contribution
                                    AssembleRHS(b, RHS_Contribution, EquationId);
                            }
                    }
                }


        KRATOS_CATCH("")
    }
    //**************************************************************************
    //**************************************************************************

    void SetUpDofSet(
        typename TSchemeType::Pointer pScheme,
        ModelPart& r_model_part
    )
    {
        KRATOS_TRY;

		if (this->GetEchoLevel() > 1 && r_model_part.GetCommunicator().MyPID() == 0)
		{
			std::cout << "Setting up the dofs" << std::endl;
		}

		//Gets the array of elements from the modeler
		ElementsArrayType& pElements = r_model_part.Elements();
		const int nelements = static_cast<int>(pElements.size());

		Element::DofsVectorType ElementalDofList;

		ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();



		unsigned int nthreads = OpenMPUtils::GetNumThreads();


		//        typedef boost::fast_pool_allocator< Node<3>::DofType::Pointer > allocator_type;
		// 		typedef std::unordered_set < Node<3>::DofType::Pointer,
		// 			dof_iterator_hash,
		// 			dof_iterator_equal,
		// 			allocator_type 	>  set_type;

#ifdef USE_GOOGLE_HASH
		typedef google::dense_hash_set < Node<3>::DofType::Pointer, dof_iterator_hash>  set_type;
#else
		typedef std::unordered_set < Node<3>::DofType::Pointer, dof_iterator_hash>  set_type;
#endif
		//         


		std::vector<set_type> dofs_aux_list(nthreads);
		// 		std::vector<allocator_type> allocators(nthreads);

			for (int i = 0; i < static_cast<int>(nthreads); i++)
			{
#ifdef USE_GOOGLE_HASH
				dofs_aux_list[i].set_empty_key(Node<3>::DofType::Pointer());
#else
				// 			dofs_aux_list[i] = set_type( allocators[i]);
				dofs_aux_list[i].reserve(nelements);
#endif
			}
		
#pragma omp parallel for firstprivate(nelements, ElementalDofList)
			for (int i = 0; i < static_cast<int>(nelements); i++)
			{
				typename ElementsArrayType::iterator it = pElements.begin() + i;
				const unsigned int this_thread_id = OpenMPUtils::ThisThread();

				// gets list of Dof involved on every element
				pScheme->GetElementalDofList(*(it.base()), ElementalDofList, CurrentProcessInfo);

				dofs_aux_list[this_thread_id].insert(ElementalDofList.begin(), ElementalDofList.end());
			}

		ConditionsArrayType& pConditions = r_model_part.Conditions();
		const int nconditions = static_cast<int>(pConditions.size());
#pragma omp parallel for firstprivate(nconditions, ElementalDofList)
		for (int i = 0; i < nconditions; i++)
		{
			typename ConditionsArrayType::iterator it = pConditions.begin() + i;
			const unsigned int this_thread_id = OpenMPUtils::ThisThread();

			// gets list of Dof involved on every element
			pScheme->GetConditionDofList(*(it.base()), ElementalDofList, CurrentProcessInfo);
			dofs_aux_list[this_thread_id].insert(ElementalDofList.begin(), ElementalDofList.end());

		}

			//here we do a reduction in a tree so to have everything on thread 0
			unsigned int old_max = nthreads;
		unsigned int new_max = ceil(0.5*static_cast<double>(old_max));
		while (new_max >= 1 && new_max != old_max)
		{
			//just for debugging
		/*	std::cout << "old_max" << old_max << " new_max:" << new_max << std::endl;
			for (int i = 0; i < new_max; i++)
			{
				if (i + new_max < old_max)
				{
					std::cout << i << " - " << i + new_max << std::endl;
				}
			}
			std::cout << "********************" << std::endl;
*/
#pragma omp parallel for
			for (int i = 0; i < static_cast<int>(new_max); i++)
			{
				if (i + new_max < old_max)
				{
					dofs_aux_list[i].insert(dofs_aux_list[i + new_max].begin(), dofs_aux_list[i + new_max].end());
					dofs_aux_list[i + new_max].clear();
				}
			}

			old_max = new_max;
			new_max = ceil(0.5*static_cast<double>(old_max));

		}

			DofsArrayType Doftemp;
		BaseType::mDofSet = DofsArrayType();

		Doftemp.reserve(dofs_aux_list[0].size());
		for (auto it = dofs_aux_list[0].begin(); it != dofs_aux_list[0].end(); it++)
		{
			Doftemp.push_back(it->get());
		}
		Doftemp.Sort();

		BaseType::mDofSet = Doftemp;

		//throws an execption if there are no Degrees of freedom involved in the analysis
		if (BaseType::mDofSet.size() == 0)
			KRATOS_THROW_ERROR(std::logic_error, "No degrees of freedom!", "");

			BaseType::mDofSetIsInitialized = true;
		if (this->GetEchoLevel() > 2 && r_model_part.GetCommunicator().MyPID() == 0)
		{
			std::cout << "finished setting up the dofs" << std::endl;
		}

#ifdef _OPENMP
			if (mlock_array.size() != 0)
			{
				for (int i = 0; i < static_cast<int>(mlock_array.size()); i++)
					omp_destroy_lock(&mlock_array[i]);
			}

		mlock_array.resize(BaseType::mDofSet.size());

		for (int i = 0; i < static_cast<int>(mlock_array.size()); i++)
			omp_init_lock(&mlock_array[i]);
#endif

        KRATOS_CATCH("");
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

        //


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
        for (it2 = BaseType::mDofSet.ptr_begin(); it2 != BaseType::mDofSet.ptr_end(); ++it2)
        {
            if ((*it2)->IsFixed())
            {
                i = (*it2)->EquationId();
                i -= systemsize;
                /*KRATOS_WATCH((*it2)->GetSolutionStepReactionValue());
                KRATOS_WATCH(ReactionsVector[i]);*/
                (*it2)->GetSolutionStepReactionValue() = -ReactionsVector[i];
            }
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

    /**
    this function is intended to be called at the end of the solution step to clean up memory
    storage not needed
     */
    void Clear()
    {
        this->mDofSet = DofsArrayType();

        if (this->mpReactionsVector != NULL)
            TSparseSpace::Clear((this->mpReactionsVector));
        // 			this->mReactionsVector = TSystemVectorType();

        this->mpLinearSystemSolver->Clear();

        if (this->GetEchoLevel() > 1)
        {
            std::cout << "ResidualBasedEliminationBuilderAndSolver Clear Function called" << std::endl;
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
    
    void Assemble(
        TSystemMatrixType& A,
        TSystemVectorType& b,
        const LocalSystemMatrixType& LHS_Contribution,
        const LocalSystemVectorType& RHS_Contribution,
        const Element::EquationIdVectorType& EquationId
#ifdef _OPENMP
        ,std::vector< omp_lock_t >& lock_array
#endif
    )
    {
        unsigned int local_size = LHS_Contribution.size1();

        for (unsigned int i_local = 0; i_local < local_size; i_local++)
        {
            unsigned int i_global = EquationId[i_local];

            if (i_global < BaseType::mEquationSystemSize)
            {
#ifdef _OPENMP
                omp_set_lock(&lock_array[i_global]);
#endif
                b[i_global] += RHS_Contribution(i_local);
                for (unsigned int j_local = 0; j_local < local_size; j_local++)
                {
                    unsigned int j_global = EquationId[j_local];
                    if (j_global < BaseType::mEquationSystemSize)
                    {
                        A(i_global, j_global) += LHS_Contribution(i_local, j_local);
                    }
                }
#ifdef _OPENMP
                omp_unset_lock(&lock_array[i_global]);
#endif

            }
            //note that assembly on fixed rows is not performed here
        }
    }

    
    //**************************************************************************
	virtual void ConstructMatrixStructure(
		TSystemMatrixType& A,
		ElementsContainerType& rElements,
		ConditionsArrayType& rConditions,
		ProcessInfo& CurrentProcessInfo)
	{
		//filling with zero the matrix (creating the structure)
		Timer::Start("MatrixStructure");

		const std::size_t equation_size = BaseType::mEquationSystemSize;

#ifdef USE_GOOGLE_HASH
		std::vector<google::dense_hash_set<std::size_t> > indices(equation_size);
		const std::size_t empty_key = 2 * equation_size + 10;
#else
		std::vector<std::unordered_set<std::size_t> > indices(equation_size);
#endif

#pragma omp parallel for firstprivate(equation_size)
		for (int iii = 0; iii < static_cast<int>(equation_size); iii++)
		{
#ifdef USE_GOOGLE_HASH
			indices[iii].set_empty_key(empty_key);
#else
			indices[iii].reserve(40);
#endif            
		}

		Element::EquationIdVectorType ids(3, 0);

		const int nelements = static_cast<int>(rElements.size());
#pragma omp parallel for firstprivate(nelements, ids)
		for (int iii = 0; iii<nelements; iii++)
		{
			typename ElementsContainerType::iterator i_element = rElements.begin() + iii;
			(i_element)->EquationIdVector(ids, CurrentProcessInfo);

			for (std::size_t i = 0; i < ids.size(); i++)
			{
				if (ids[i] < BaseType::mEquationSystemSize)
				{
#ifdef _OPENMP
                                    omp_set_lock(&mlock_array[ids[i]]);
#endif
					auto& row_indices = indices[ids[i]];
					for (auto it = ids.begin(); it != ids.end(); it++)
					{
						if (*it < BaseType::mEquationSystemSize)
							row_indices.insert(*it);
					}
#ifdef _OPENMP
					omp_unset_lock(&mlock_array[ids[i]]);
#endif
				}
			}

		}

		const int nconditions = static_cast<int>(rConditions.size());
#pragma omp parallel for firstprivate(nconditions, ids)
		for (int iii = 0; iii<nconditions; iii++)
		{
			typename ConditionsArrayType::iterator i_condition = rConditions.begin() + iii;
			(i_condition)->EquationIdVector(ids, CurrentProcessInfo);
			for (std::size_t i = 0; i < ids.size(); i++)
			{
				if (ids[i] < BaseType::mEquationSystemSize)
				{
#ifdef _OPENMP
					omp_set_lock(&mlock_array[ids[i]]);
#endif
					auto& row_indices = indices[ids[i]];
					for (auto it = ids.begin(); it != ids.end(); it++)
					{
						if (*it < BaseType::mEquationSystemSize)
							row_indices.insert(*it);
					}
#ifdef _OPENMP
					omp_unset_lock(&mlock_array[ids[i]]);
#endif
				}
			}
		}

		//count the row sizes
		unsigned int nnz = 0;
		for (unsigned int i = 0; i < indices.size(); i++)
			nnz += indices[i].size();

		A = boost::numeric::ublas::compressed_matrix<double>(indices.size(), indices.size(), nnz);

		double* Avalues = A.value_data().begin();
		std::size_t* Arow_indices = A.index1_data().begin();
		std::size_t* Acol_indices = A.index2_data().begin();

		//filling the index1 vector - DO NOT MAKE PARALLEL THE FOLLOWING LOOP!
		Arow_indices[0] = 0;
		for (int i = 0; i < static_cast<int>(A.size1()); i++)
			Arow_indices[i + 1] = Arow_indices[i] + indices[i].size();



#pragma omp parallel for
		for (int i = 0; i < static_cast<int>(A.size1()); i++)
		{
			const unsigned int row_begin = Arow_indices[i];
			const unsigned int row_end = Arow_indices[i + 1];
			unsigned int k = row_begin;
			for (auto it = indices[i].begin(); it != indices[i].end(); it++)
			{
				Acol_indices[k] = *it;
				Avalues[k] = 0.0;
				k++;
			}

			std::sort(&Acol_indices[row_begin], &Acol_indices[row_end]);

		}

		A.set_filled(indices.size() + 1, nnz);

		Timer::Stop("MatrixStructure");
	}

//    virtual void ConstructMatrixStructure(
//        TSystemMatrixType& A,
//        ElementsContainerType& rElements,
//        ConditionsArrayType& rConditions,
//        ProcessInfo& CurrentProcessInfo)
//    {
//
//        std::size_t equation_size = A.size1();
//        std::vector<std::vector<std::size_t> > indices(equation_size);
//        //				std::vector<std::vector<std::size_t> > dirichlet_indices(TSystemSpaceType::Size1(mDirichletMatrix));
//
//        Element::EquationIdVectorType ids(3, 0);
//        for (typename ElementsContainerType::iterator i_element = rElements.begin(); i_element != rElements.end(); i_element++)
//        {
//            (i_element)->EquationIdVector(ids, CurrentProcessInfo);
//
//            for (std::size_t i = 0; i < ids.size(); i++)
//                if (ids[i] < equation_size)
//                {
//                    std::vector<std::size_t>& row_indices = indices[ids[i]];
//                    for (std::size_t j = 0; j < ids.size(); j++)
//                        if (ids[j] < equation_size)
//                        {
//                            AddUnique(row_indices, ids[j]);
//                            //indices[ids[i]].push_back(ids[j]);
//                        }
//                }
//
//        }
//
//        for (typename ConditionsArrayType::iterator i_condition = rConditions.begin(); i_condition != rConditions.end(); i_condition++)
//        {
//            (i_condition)->EquationIdVector(ids, CurrentProcessInfo);
//            for (std::size_t i = 0; i < ids.size(); i++)
//                if (ids[i] < equation_size)
//                {
//                    std::vector<std::size_t>& row_indices = indices[ids[i]];
//                    for (std::size_t j = 0; j < ids.size(); j++)
//                        if (ids[j] < equation_size)
//                        {
//                            AddUnique(row_indices, ids[j]);
//                            //	indices[ids[i]].push_back(ids[j]);
//                        }
//                }
//        }
//
//        //allocating the memory needed
//        int data_size = 0;
//        for (std::size_t i = 0; i < indices.size(); i++)
//        {
//            data_size += indices[i].size();
//        }
//        A.reserve(data_size, false);
//
//        //filling with zero the matrix (creating the structure)
//        Timer::Start("MatrixStructure");
//#ifndef _OPENMP
//        for (std::size_t i = 0; i < indices.size(); i++)
//        {
//            std::vector<std::size_t>& row_indices = indices[i];
//            std::sort(row_indices.begin(), row_indices.end());
//
//            for (std::vector<std::size_t>::iterator it = row_indices.begin(); it != row_indices.end(); it++)
//            {
//                A.push_back(i, *it, 0.00);
//            }
//            row_indices.clear();
//        }
//#else
//        int number_of_threads = omp_get_max_threads();
//        vector<unsigned int> matrix_partition;
//        CreatePartition(number_of_threads, indices.size(), matrix_partition);
//        if (this->GetEchoLevel() > 2)
//        {
//            KRATOS_WATCH(matrix_partition);
//        }
//        for (int k = 0; k < number_of_threads; k++)
//        {
//            #pragma omp parallel
//            if (omp_get_thread_num() == k)
//            {
//                for (std::size_t i = matrix_partition[k]; i < matrix_partition[k + 1]; i++)
//                {
//                    std::vector<std::size_t>& row_indices = indices[i];
//                    std::sort(row_indices.begin(), row_indices.end());
//
//                    for (std::vector<std::size_t>::iterator it = row_indices.begin(); it != row_indices.end(); it++)
//                    {
//                        A.push_back(i, *it, 0.00);
//                    }
//                    row_indices.clear();
//                }
//            }
//        }
//#endif
//        Timer::Stop("MatrixStructure");
//    }

    //**************************************************************************

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
#ifdef _OPENMP
	std::vector< omp_lock_t > mlock_array;
#endif 
    /*@} */
    /**@name Private Operators*/
    /*@{ */


    /*@} */
    /**@name Private Operations*/
    /*@{ */



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


	void AssembleRHS(
		TSystemVectorType& b,
		const LocalSystemVectorType& RHS_Contribution,
		const Element::EquationIdVectorType& EquationId
	)
	{
		unsigned int local_size = RHS_Contribution.size();

                if (BaseType::mCalculateReactionsFlag == false)
                {
                    for (unsigned int i_local = 0; i_local < local_size; i_local++)
                    {
                            const unsigned int i_global = EquationId[i_local];

                            if (i_global < BaseType::mEquationSystemSize) //free dof
                            {
                                    // ASSEMBLING THE SYSTEM VECTOR
                                    double& b_value = b[i_global];
                                    const double& rhs_value = RHS_Contribution[i_local];

                                    #pragma omp atomic
                                    b_value += rhs_value;
                            }
                    }
                }
                else
                {
                    TSystemVectorType& ReactionsVector = *BaseType::mpReactionsVector;
                    for (unsigned int i_local = 0; i_local < local_size; i_local++)
                    {
                            const unsigned int i_global = EquationId[i_local];

                            if (i_global < BaseType::mEquationSystemSize) //free dof
                            {
                                    // ASSEMBLING THE SYSTEM VECTOR
                                    double& b_value = b[i_global];
                                    const double& rhs_value = RHS_Contribution[i_local];

                                    #pragma omp atomic
                                    b_value += rhs_value;
                            }
                            else //fixed dof
                            {
                                    double& b_value = ReactionsVector[i_global - BaseType::mEquationSystemSize];
                                    const double& rhs_value = RHS_Contribution[i_local];

                                    #pragma omp atomic
                                    b_value += rhs_value;
                            }
                    }
                }
	}
	
    //**************************************************************************

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

}; /* Class ResidualBasedEliminationBuilderAndSolver */

/*@} */

/**@name Type Definitions */
/*@{ */


/*@} */

} /* namespace Kratos.*/

#endif /* KRATOS_RESIDUAL_BASED_ELIMINATION_BUILDER_AND_SOLVER  defined */

