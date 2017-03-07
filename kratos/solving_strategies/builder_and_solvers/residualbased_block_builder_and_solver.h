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
#if !defined(KRATOS_RESIDUAL_BASED_BLOCK_BUILDER_AND_SOLVER )
#define  KRATOS_RESIDUAL_BASED_BLOCK_BUILDER_AND_SOLVER


/* System includes */

#include "utilities/openmp_utils.h"


#include <unordered_set>
/* External includes */
#include "boost/smart_ptr.hpp"

#include "utilities/timer.h"

/* Project includes */
#include "includes/define.h"
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"
#include "includes/model_part.h"
#include "includes/kratos_flags.h"

// #define USE_GOOGLE_HASH
#ifdef USE_GOOGLE_HASH
    #include "sparsehash/dense_hash_set" //included in external libraries
#endif

// #include <iostream>
// #include <fstream>

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
class ResidualBasedBlockBuilderAndSolver
    : public BuilderAndSolver< TSparseSpace, TDenseSpace, TLinearSolver >
{
public:
    /**@name Type Definitions */
    /*@{ */
    KRATOS_CLASS_POINTER_DEFINITION(ResidualBasedBlockBuilderAndSolver);


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
            return (((it1->Id() == it2->Id() && (it1->GetVariable()).Key())==(it2->GetVariable()).Key()));
        }
    };
    /*@} */
    /**@name Life Cycle
     */
    /*@{ */

    /** Constructor.
     */
    ResidualBasedBlockBuilderAndSolver(
        typename TLinearSolver::Pointer pNewLinearSystemSolver)
        : BuilderAndSolver< TSparseSpace, TDenseSpace, TLinearSolver >(pNewLinearSystemSolver)
    {

        /* 			std::cout << "using the standard builder and solver " << std::endl; */

    }

    /** Destructor.
     */
    virtual ~ResidualBasedBlockBuilderAndSolver()
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
        TSystemVectorType& b) override
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

        #pragma omp parallel for firstprivate(nelements, LHS_Contribution, RHS_Contribution, EquationId )
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


        #pragma omp parallel for firstprivate(nconditions, LHS_Contribution, RHS_Contribution, EquationId )
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


        // equation_ids.close();
        // KRATOS_THROW_ERROR(std::logic_error,"i want to stop here :-D","")

        double stop_build = OpenMPUtils::GetCurrentTime();
        if (this->GetEchoLevel() >= 1 && r_model_part.GetCommunicator().MyPID() == 0)
            std::cout << "build time: " << stop_build - start_build << std::endl;

        //for (int i = 0; i < A_size; i++)
        //    omp_destroy_lock(&lock_array[i]);
        if (this->GetEchoLevel() > 2 && r_model_part.GetCommunicator().MyPID() == 0)
        {
            KRATOS_WATCH("finished parallel building");
        }




        KRATOS_CATCH("")

    }

    //**************************************************************************
    //**************************************************************************

    void BuildLHS(
        typename TSchemeType::Pointer pScheme,
        ModelPart& r_model_part,
        TSystemMatrixType& A) override
    {
		KRATOS_TRY

			TSystemVectorType tmp(A.size1(), 0.0);
		this->Build(pScheme, r_model_part, A, tmp);

        KRATOS_CATCH("")

    }

    //**************************************************************************
    //**************************************************************************

    void BuildLHS_CompleteOnFreeRows(
        typename TSchemeType::Pointer pScheme,
        ModelPart& r_model_part,
        TSystemMatrixType& A) override
    {
        KRATOS_TRY

			TSystemVectorType tmp(A.size1(), 0.0);
		this->Build(pScheme, r_model_part, A, tmp);

        KRATOS_CATCH("")

    }

    //**************************************************************************
    //**************************************************************************

    void SystemSolve(
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b
    ) override
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
        TSystemVectorType& b) override
    {
        KRATOS_TRY

        Timer::Start("Build");

        Build(pScheme, r_model_part, A, b);

        Timer::Stop("Build");


        ApplyDirichletConditions(pScheme, r_model_part, A, Dx, b);

        if (this->GetEchoLevel() == 3)
        {
            std::cout << "before the solution of the system" << std::endl;
            std::cout << "System Matrix = " << A << std::endl;
            std::cout << "unknowns vector = " << Dx << std::endl;
            std::cout << "RHS vector = " << b << std::endl;
        }

        double start_solve = OpenMPUtils::GetCurrentTime();
        Timer::Start("Solve");

        SystemSolveWithPhysics(A, Dx, b, r_model_part);

        Timer::Stop("Solve");
        double stop_solve = OpenMPUtils::GetCurrentTime();
        if (this->GetEchoLevel() >=1 && r_model_part.GetCommunicator().MyPID() == 0)
            std::cout << "system solve time: " << stop_solve - start_solve << std::endl;

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
        TSystemVectorType& b) override
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
        TSystemVectorType& b) override
    {
        KRATOS_TRY

        BuildRHSNoDirichlet(pScheme,r_model_part,b);

		const int ndofs = static_cast<int>(BaseType::mDofSet.size());

		//NOTE: dofs are assumed to be numbered consecutively in the BlockBuilderAndSolver
		#pragma omp parallel for firstprivate(ndofs)
		for (int k = 0; k<ndofs; k++)
		{
			typename DofsArrayType::iterator dof_iterator = BaseType::mDofSet.begin() + k;
			const std::size_t i = dof_iterator->EquationId();

			if (dof_iterator->IsFixed())
				b[i] = 0.0f;
		}

        KRATOS_CATCH("")

    }
    //**************************************************************************
    //**************************************************************************

    void SetUpDofSet(
        typename TSchemeType::Pointer pScheme,
        ModelPart& r_model_part
    ) override
    {
        KRATOS_TRY;

        if( this->GetEchoLevel() > 1 && r_model_part.GetCommunicator().MyPID() == 0)
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

        if( this->GetEchoLevel() > 2)
        {
            std::cout << "Number of threads" << nthreads << "\n" << std::endl;
        }
        for (int i = 0; i < static_cast<int>(nthreads); i++)
        {
#ifdef USE_GOOGLE_HASH
            dofs_aux_list[i].set_empty_key(Node<3>::DofType::Pointer());
#else
// 			dofs_aux_list[i] = set_type( allocators[i]);
            dofs_aux_list[i].reserve(nelements);
#endif
        }
        if( this->GetEchoLevel() > 2)
        {
            KRATOS_WATCH("Initializing element loop")
        }
        #pragma omp parallel for firstprivate(nelements, ElementalDofList)
        for (int i = 0; i < nelements; i++)
        {
            typename ElementsArrayType::iterator it = pElements.begin() + i;
            const unsigned int this_thread_id = OpenMPUtils::ThisThread();

            // gets list of Dof involved on every element
            pScheme->GetElementalDofList(*(it.base()), ElementalDofList, CurrentProcessInfo);

            dofs_aux_list[this_thread_id].insert(ElementalDofList.begin(), ElementalDofList.end());
        }
        if( this->GetEchoLevel() > 2)
        {
            std::cout << "Initializing condition loop\n" << std::endl;
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

        if( this->GetEchoLevel() > 2)
        {
            std::cout << "Initializing tree reduction\n" << std::endl;
        }
        // Here we do a reduction in a tree so to have everything on thread 0
        unsigned int old_max = nthreads;
        unsigned int new_max = ceil(0.5*static_cast<double>(old_max));
        while (new_max>=1 && new_max != old_max)
        {
            if( this->GetEchoLevel() > 2)
            {
                //just for debugging
                std::cout << "old_max" << old_max << " new_max:" << new_max << std::endl;
                for (int i = 0; i < static_cast<int>(new_max); i++)
                {
                    if (i + new_max < old_max)
                    {
                        std::cout << i << " - " << i+new_max << std::endl;
                    }
                }
                std::cout << "********************" << std::endl;
            }
            
            #pragma omp parallel for
            for (int i = 0; i < static_cast<int>(new_max); i++)
            {
                if (i + new_max < old_max)
                {
                    dofs_aux_list[i].insert(dofs_aux_list[i+new_max].begin(), dofs_aux_list[i+new_max].end());
                    dofs_aux_list[i+new_max].clear();
                }
            }
            
            old_max = new_max;
            new_max = ceil(0.5*static_cast<double>(old_max));
            
        }

        if( this->GetEchoLevel() > 2)
        {
            std::cout << "Initializing ordered array filling\n" << std::endl;
        }

        DofsArrayType Doftemp;
        BaseType::mDofSet = DofsArrayType();

        Doftemp.reserve(dofs_aux_list[0].size());
        for (auto it= dofs_aux_list[0].begin(); it!= dofs_aux_list[0].end(); it++)
        {
            Doftemp.push_back( it->get() );
        }
        Doftemp.Sort();

        BaseType::mDofSet = Doftemp;

        //Throws an exception if there are no Degrees Of Freedom involved in the analysis
        if (BaseType::mDofSet.size() == 0)
        {
            KRATOS_THROW_ERROR(std::logic_error, "No degrees of freedom!", "");
        }
        if( this->GetEchoLevel() > 2)
        {
            KRATOS_WATCH(BaseType::mDofSet.size())
        }
        BaseType::mDofSetIsInitialized = true;
        if( this->GetEchoLevel() > 2 && r_model_part.GetCommunicator().MyPID() == 0)
        {
            std::cout << "Finished setting up the dofs" << std::endl;
        }

        if( this->GetEchoLevel() > 2)
        {
            KRATOS_WATCH("Initializing lock array")
        }
        
#ifdef _OPENMP
        if (mlock_array.size() != 0)
        {
            for (int i = 0; i < static_cast<int>(mlock_array.size()); i++)
            {
                omp_destroy_lock(&mlock_array[i]);
            }
        }
        mlock_array.resize(BaseType::mDofSet.size());
        
        for (int i = 0; i < static_cast<int>(mlock_array.size()); i++)
        {
            omp_init_lock(&mlock_array[i]);
        }
#endif
        if( this->GetEchoLevel() > 2)
        {
            std::cout << "End of setupdofset\n" << std::endl;
        }
        
        KRATOS_CATCH("");
    }

    //**************************************************************************
    //**************************************************************************

    void SetUpSystem(
        ModelPart& r_model_part
    ) override
    {

        //int free_id = 0;
        BaseType::mEquationSystemSize = BaseType::mDofSet.size();
        int ndofs = static_cast<int>(BaseType::mDofSet.size());

        #pragma omp parallel for firstprivate(ndofs)
        for (int i = 0; i < static_cast<int>(ndofs); i++)
        {
            typename DofsArrayType::iterator dof_iterator = BaseType::mDofSet.begin() + i;
            dof_iterator->SetEquationId(i);

        }

        //for (typename DofsArrayType::iterator dof_iterator = BaseType::mDofSet.begin(); dof_iterator != BaseType::mDofSet.end(); ++dof_iterator)
        //    dof_iterator->SetEquationId(free_id++);



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
    ) override
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




        KRATOS_CATCH("")

    }



    //**************************************************************************
    //**************************************************************************

    void InitializeSolutionStep(
        ModelPart& r_model_part,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
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
        TSystemVectorType& b) override
    {
    }


    //**************************************************************************
    //**************************************************************************

    void CalculateReactions(
        typename TSchemeType::Pointer pScheme,
        ModelPart& r_model_part,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
        TSparseSpace::SetToZero(b);

        //refresh RHS to have the correct reactions
        BuildRHSNoDirichlet(pScheme, r_model_part, b);

		const int ndofs = static_cast<int>(BaseType::mDofSet.size());

		//NOTE: dofs are assumed to be numbered consecutively in the BlockBuilderAndSolver
		#pragma omp parallel for firstprivate(ndofs)
		for (int k = 0; k<ndofs; k++)
		{
			typename DofsArrayType::iterator dof_iterator = BaseType::mDofSet.begin() + k;

			if (dof_iterator->IsFixed())
			{
				const int i = (dof_iterator)->EquationId();
				(dof_iterator)->GetSolutionStepReactionValue() = -b[i];
			}
		}

        //KRATOS_WATCH(__LINE__)
    }

    //**************************************************************************
    //**************************************************************************

    void ApplyDirichletConditions(
        typename TSchemeType::Pointer pScheme,
        ModelPart& r_model_part,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
        std::size_t system_size = A.size1();
        std::vector<double> scaling_factors (system_size, 0.0f);

        const int ndofs = static_cast<int>(BaseType::mDofSet.size());

        //NOTE: dofs are assumed to be numbered consecutively in the BlockBuilderAndSolver
        #pragma omp parallel for firstprivate(ndofs)
        for(int k = 0; k<ndofs; k++)
        {
            typename DofsArrayType::iterator dof_iterator = BaseType::mDofSet.begin() + k;
            if(dof_iterator->IsFixed())
                scaling_factors[k] = 0.0f;
            else
                scaling_factors[k] = 1.0f;

        }

        double* Avalues = A.value_data().begin();
        std::size_t* Arow_indices = A.index1_data().begin();
        std::size_t* Acol_indices = A.index2_data().begin();

        //detect if there is a line of all zeros and set the diagonal to a 1 if this happens
        #pragma omp parallel for firstprivate(system_size)
        for(int k = 0; k < static_cast<int>(system_size); ++k)
        {
            std::size_t col_begin = Arow_indices[k];
            std::size_t col_end = Arow_indices[k+1];
            bool empty = true;
            for (std::size_t j = col_begin; j < col_end; ++j)
            {
                if(Avalues[j] != 0.0)
                {
                    empty = false;
                    break;
                }
            }

            if(empty == true)
            {
                A(k,k) = 1.0;
                b[k] = 0.0;
            }
        }

        #pragma omp parallel for
        for (int k = 0; k < static_cast<int>(system_size); ++k)
        {
            std::size_t col_begin = Arow_indices[k];
            std::size_t col_end = Arow_indices[k+1];
            double k_factor = scaling_factors[k];
            if (k_factor == 0)
            {
                // zero out the whole row, except the diagonal
                for (std::size_t j = col_begin; j < col_end; ++j)
                    if (static_cast<int>(Acol_indices[j]) != k )
                        Avalues[j] = 0.0;

                // zero out the RHS
                b[k] = 0.0;
            }
            else
            {
                // zero out the column which is associated with the zero'ed row
                for (std::size_t j = col_begin; j < col_end; ++j)
                    if(scaling_factors[ Acol_indices[j] ] == 0 )
                        Avalues[j] = 0.0;
            }
        }
    }



    /**
    this function is intended to be called at the end of the solution step to clean up memory
    storage not needed
     */
    void Clear() override
    {
        this->mDofSet = DofsArrayType();

        this->mpLinearSystemSolver->Clear();

        if (this->GetEchoLevel() > 0)
        {
            std::cout << "ResidualBasedBlockBuilderAndSolver Clear Function called" << std::endl;
        }

#ifdef _OPENMP
        for (int i = 0; i < static_cast<int>(mlock_array.size()); i++)
            omp_destroy_lock(&mlock_array[i]);
        mlock_array.resize(0);
#endif 
        
    }

    /**
     * This function is designed to be called once to perform all the checks needed
     * on the input provided. Checks can be "expensive" as the function is designed
     * to catch user's errors.
     * @param r_model_part
     * @return 0 all ok
     */
    virtual int Check(ModelPart& r_mUSE_GOOGLE_HASHodel_part) override
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
#ifdef _OPENMP
    std::vector< omp_lock_t > mlock_array;
#endif

    /*@} */
    /**@name Protected Operators*/
    /*@{ */
    //**************************************************************************

    virtual void ConstructMatrixStructure(
        TSystemMatrixType& A,
        ElementsContainerType& rElements,
        ConditionsArrayType& rConditions,
        ProcessInfo& CurrentProcessInfo)
    {
        //filling with zero the matrix (creating the structure)
        Timer::Start("MatrixStructure");

        const std::size_t equation_size = BaseType::mDofSet.size();
        
#ifdef USE_GOOGLE_HASH
        std::vector<google::dense_hash_set<std::size_t> > indices(equation_size);
        const std::size_t empty_key = 2*equation_size + 10;
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
        for(int iii=0; iii<nelements; iii++)
        {
            typename ElementsContainerType::iterator i_element = rElements.begin() + iii;
            (i_element)->EquationIdVector(ids, CurrentProcessInfo);

            for (std::size_t i = 0; i < ids.size(); i++)
            {
#ifdef _OPENMP
                omp_set_lock(&mlock_array[ids[i]]);
#endif 
                auto& row_indices = indices[ids[i]];
                row_indices.insert(ids.begin(), ids.end());

#ifdef _OPENMP
                omp_unset_lock(&mlock_array[ids[i]]);
#endif 
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
#ifdef _OPENMP
                omp_set_lock(&mlock_array[ids[i]]);
#endif
                auto& row_indices = indices[ids[i]];
                row_indices.insert(ids.begin(), ids.end());
#ifdef _OPENMP
                omp_unset_lock(&mlock_array[ids[i]]);
#endif
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
            Arow_indices[i+1] = Arow_indices[i] + indices[i].size();



        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(A.size1()); i++)
        {
            const unsigned int row_begin = Arow_indices[i];
            const unsigned int row_end = Arow_indices[i+1];
            unsigned int k = row_begin;
            for (auto it = indices[i].begin(); it != indices[i].end(); it++)
            {
                Acol_indices[k] = *it;
                Avalues[k] = 0.0;
                k++;
            }

			indices[i].clear(); //deallocating the memory 

            std::sort(&Acol_indices[row_begin], &Acol_indices[row_end]);

        }

        A.set_filled(indices.size()+1, nnz);

        Timer::Stop("MatrixStructure");
    }

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

            for (unsigned int j_local = 0; j_local < local_size; j_local++)
            {
                unsigned int j_global = EquationId[j_local];

                A(i_global, j_global) += LHS_Contribution(i_local, j_local);
            }
        }

    }


    void Assemble(
        TSystemMatrixType& A,
        TSystemVectorType& b,
        const LocalSystemMatrixType& LHS_Contribution,
        const LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId
#ifdef _OPENMP
        ,std::vector< omp_lock_t >& lock_array
#endif
    )
    {
        unsigned int local_size = LHS_Contribution.size1();

        for (unsigned int i_local = 0; i_local < local_size; i_local++)
        {
            unsigned int i_global = EquationId[i_local];

#ifdef _OPENMP
            omp_set_lock(&lock_array[i_global]);
#endif
            b[i_global] += RHS_Contribution(i_local);

            AssembleRowContribution(A, LHS_Contribution, i_global, i_local, EquationId);

#ifdef _OPENMP
            omp_unset_lock(&lock_array[i_global]);
#endif


            //note that computation of reactions is not performed here!
        }
    }

    
    //**************************************************************************

    void AssembleRHS(
        TSystemVectorType& b,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId
    )
    {
        unsigned int local_size = RHS_Contribution.size();


        for (unsigned int i_local = 0; i_local < local_size; i_local++)
        {
            unsigned int i_global = EquationId[i_local];

            // ASSEMBLING THE SYSTEM VECTOR
            double& b_value = b[i_global];
            const double& rhs_value = RHS_Contribution[i_local];

            #pragma omp atomic
            b_value += rhs_value;


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

    void BuildRHSNoDirichlet(
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

        //contributions to the system
        LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0, 0);
        LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

        //vector containing the localization in the system of the different
        //terms
        Element::EquationIdVectorType EquationId;

        // assemble all elements
        //for (typename ElementsArrayType::ptr_iterator it = pElements.ptr_begin(); it != pElements.ptr_end(); ++it)

        int nelements = static_cast<int>(pElements.size());
        #pragma omp parallel for firstprivate(nelements, RHS_Contribution, EquationId)
        for(int i=0; i<nelements; i++)
        {
            typename ElementsArrayType::iterator it = pElements.begin() + i;
            //detect if the element is active or not. If the user did not make any choice the element
            //is active by default
            bool element_is_active = true;
            if( (it)->IsDefined(ACTIVE) )
                element_is_active = (it)->Is(ACTIVE);

            if(element_is_active)
            {
                //calculate elemental Right Hand Side Contribution
                pScheme->Calculate_RHS_Contribution(*(it.base()), RHS_Contribution, EquationId, CurrentProcessInfo);

                //assemble the elemental contribution
                AssembleRHS(b, RHS_Contribution, EquationId);
            }
        }

        LHS_Contribution.resize(0, 0, false);
        RHS_Contribution.resize(0, false);

        // assemble all conditions
        //for (typename ConditionsArrayType::ptr_iterator it = ConditionsArray.ptr_begin(); it != ConditionsArray.ptr_end(); ++it)
        int nconditions = static_cast<int>(ConditionsArray.size());
        #pragma omp parallel for firstprivate(nconditions, RHS_Contribution, EquationId)
        for (int i = 0; i<nconditions; i++)
        {
            auto it = ConditionsArray.begin() + i;
            //detect if the element is active or not. If the user did not make any choice the element
            //is active by default
            bool condition_is_active = true;
            if( (it)->IsDefined(ACTIVE) )
                condition_is_active = (it)->Is(ACTIVE);

            if(condition_is_active)
            {

                //calculate elemental contribution
                pScheme->Condition_Calculate_RHS_Contribution(*(it.base()), RHS_Contribution, EquationId, CurrentProcessInfo);

                //assemble the elemental contribution
                AssembleRHS(b, RHS_Contribution, EquationId);
            }
        }

        KRATOS_CATCH("")

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

    inline void AssembleRowContribution(TSystemMatrixType& A, const Matrix& Alocal, const unsigned int i, const unsigned int i_local, Element::EquationIdVectorType& EquationId)
    {
        double* values_vector = A.value_data().begin();
        std::size_t* index1_vector = A.index1_data().begin();
        std::size_t* index2_vector = A.index2_data().begin();

        size_t left_limit = index1_vector[i];
//	size_t right_limit = index1_vector[i+1];

        //find the first entry
        size_t last_pos = ForwardFind(EquationId[0],left_limit,index2_vector);
        size_t last_found = EquationId[0];
        values_vector[last_pos] += Alocal(i_local,0);

        //now find all of the other entries
        size_t pos = 0;
        for(unsigned int j=1; j<EquationId.size(); j++)
        {
            unsigned int id_to_find = EquationId[j];
            if(id_to_find > last_found)
                pos = ForwardFind(id_to_find,last_pos+1,index2_vector);
            else
                pos = BackwardFind(id_to_find,last_pos-1,index2_vector);

            values_vector[pos] += Alocal(i_local,j);
            last_found = id_to_find;
            last_pos = pos;
        }
    }




    inline unsigned int ForwardFind(const unsigned int id_to_find,
                                    const unsigned int start,
                                    const size_t* index_vector)
    {
        unsigned int pos = start;
        while(id_to_find != index_vector[pos]) pos++;
        return pos;
    }

    inline unsigned int BackwardFind(const unsigned int id_to_find,
                                     const unsigned int start,
                                     const size_t* index_vector)
    {
        unsigned int pos = start;
        while(id_to_find != index_vector[pos]) pos--;
        return pos;
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

}; /* Class ResidualBasedBlockBuilderAndSolver */

/*@} */

/**@name Type Definitions */
/*@{ */


/*@} */

} /* namespace Kratos.*/

#endif /* KRATOS_RESIDUAL_BASED_BLOCK_BUILDER_AND_SOLVER  defined */
