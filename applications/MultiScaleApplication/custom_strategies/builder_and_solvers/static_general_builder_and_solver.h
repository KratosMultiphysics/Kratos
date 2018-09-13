//
//   Project Name:        Kratos
//   Last Modified by:    $Author: mpetracca $
//   Date:                $Date: 2013-06-06 10:37:00 $
//   Revision:            $Revision: 1.00 $
//
//

#if !defined(KRATOS_STATIC_GENERAL_BUILDER_AND_SOLVER )
#define  KRATOS_STATIC_GENERAL_BUILDER_AND_SOLVER


/* System includes */
#include <set>
#include <iomanip>

#ifdef _OPENMP
#define STATIC_GENERAL_BUILDER_AND_SOLVER_USES_OPENMP
#endif

#ifdef STATIC_GENERAL_BUILDER_AND_SOLVER_USES_OPENMP
#include <omp.h>
#endif

#define STATIC_GENERAL_BUILDER_AND_SOLVER_OPTIMIZE_SOLUTION

/* External includes */
#include "boost/smart_ptr.hpp"
#include "utilities/timer.h"

/* Project includes */
#include "multiscale_application_variables.h"
#include "includes/define.h"
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"
#include "includes/model_part.h"

namespace Kratos
{

/** Short class definition.

Detail class definition.

Current class provides an implementation for standard builder and solving operations.

the RHS is constituted by the unbalanced loads (residual)

Degrees of freedom are reordered putting the restrained degrees of freedom at
the end of the system ordered in reverse order with respect to the DofSet.

Imposition of the dirichlet conditions is naturally dealt with as the residual already contains
this information.

Calculation of the reactions involves a cost very similiar to the calculation of the total residual

 */
template<class TSparseSpace,class TDenseSpace,class TLinearSolver>
class StaticGeneralBuilderAndSolver
    : public BuilderAndSolver< TSparseSpace, TDenseSpace, TLinearSolver >
{
public:

    KRATOS_CLASS_POINTER_DEFINITION(StaticGeneralBuilderAndSolver);

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

public:

    StaticGeneralBuilderAndSolver(typename TLinearSolver::Pointer pNewLinearSystemSolver)
        : BuilderAndSolver< TSparseSpace, TDenseSpace, TLinearSolver >(pNewLinearSystemSolver)
    {
		mTotalBuildTime = 0.0;
        mTotalSolutionTime = 0.0;
		mTotalBuildTime_accum = 0.0;
        mTotalSolutionTime_accum = 0.0;
#ifdef STATIC_GENERAL_BUILDER_AND_SOLVER_OPTIMIZE_SOLUTION
		mDoFactorization = false;
#endif // STATIC_GENERAL_BUILDER_AND_SOLVER_OPTIMIZE_SOLUTION
    }

    virtual ~StaticGeneralBuilderAndSolver()
    {
    }

public:

    void Build(typename TSchemeType::Pointer pScheme, ModelPart& r_model_part, TSystemMatrixType& A, TSystemVectorType& b)
    {
		double time_begin = Timer::GetTime();

        if (!pScheme)
            KRATOS_THROW_ERROR( std::runtime_error, "No scheme provided!", "" )

        //getting the elements from the model
        ElementsArrayType& pElements = r_model_part.Elements();

        //getting the array of the conditions
        ConditionsArrayType& ConditionsArray = r_model_part.Conditions();

        //resetting to zero the vector of reactions
        TSparseSpace::SetToZero(*(BaseType::mpReactionsVector));

        // assemble all elements
#ifndef STATIC_GENERAL_BUILDER_AND_SOLVER_USES_OPENMP

		//vector containing the localization in the system of the different terms
        Element::EquationIdVectorType EquationId;

        //contributions to the system
        LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0, 0);
        LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

        ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

        for (typename ElementsArrayType::ptr_iterator it = pElements.ptr_begin(); it != pElements.ptr_end(); ++it)
        {
            //calculate elemental contribution
            pScheme->CalculateSystemContributions(*it, LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);

            //assemble the elemental contribution
            AssembleLHS(A, LHS_Contribution, EquationId);
            AssembleRHS(b, RHS_Contribution, EquationId);

            // clean local elemental memory
            pScheme->CleanMemory(*it);
        }

        LHS_Contribution.resize(0, 0, false);
        RHS_Contribution.resize(0, false);

        // assemble all conditions
        for (typename ConditionsArrayType::ptr_iterator it = ConditionsArray.ptr_begin(); it != ConditionsArray.ptr_end(); ++it)
        {
            //calculate elemental contribution
            pScheme->Condition_CalculateSystemContributions(*it, LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);

            //assemble the elemental contribution
            AssembleLHS(A, LHS_Contribution, EquationId);
            AssembleRHS(b, RHS_Contribution, EquationId);
        }

#else
        //creating an array of lock variables of the size of the system matrix
        std::vector< omp_lock_t > lock_array(A.size1());

        int A_size = A.size1();
        for (int i = 0; i < A_size; i++)
            omp_init_lock(&lock_array[i]);

        //create a partition of the element array
        int number_of_threads = omp_get_max_threads();

        vector<unsigned int> element_partition;
        CreatePartition(number_of_threads, pElements.size(), element_partition);

        if( this->GetEchoLevel() > 2 && r_model_part.GetCommunicator().MyPID() == 0)
        {
            KRATOS_WATCH( number_of_threads )
            KRATOS_WATCH( element_partition )
        }

        double start_prod = omp_get_wtime();

        #pragma omp parallel for
        for (int k = 0; k < number_of_threads; k++)
        {
            //contributions to the system
            LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0, 0);
            LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

            //vector containing the localization in the system of the different
            //terms
            Element::EquationIdVectorType EquationId;

            ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

            typename ElementsArrayType::ptr_iterator it_begin = pElements.ptr_begin() + element_partition[k];
            typename ElementsArrayType::ptr_iterator it_end = pElements.ptr_begin() + element_partition[k + 1];

            // assemble all elements
            for (typename ElementsArrayType::ptr_iterator it = it_begin; it != it_end; ++it)
            {
                //calculate elemental contribution
				#pragma omp critical
				{
					pScheme->CalculateSystemContributions(*it, LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);
				}
                //assemble the elemental contribution
                Assemble(A, b, LHS_Contribution, RHS_Contribution, EquationId, lock_array);

                // clean local elemental memory
                pScheme->CleanMemory(*it);
            }
        }

        vector<unsigned int> condition_partition;
        CreatePartition(number_of_threads, ConditionsArray.size(), condition_partition);

        #pragma omp parallel for
        for (int k = 0; k < number_of_threads; k++)
        {
            //contributions to the system
            LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0, 0);
            LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

            Condition::EquationIdVectorType EquationId;

            ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

            typename ConditionsArrayType::ptr_iterator it_begin = ConditionsArray.ptr_begin() + condition_partition[k];
            typename ConditionsArrayType::ptr_iterator it_end = ConditionsArray.ptr_begin() + condition_partition[k + 1];

            // assemble all elements
            for (typename ConditionsArrayType::ptr_iterator it = it_begin; it != it_end; ++it)
            {
                //calculate elemental contribution
                pScheme->Condition_CalculateSystemContributions(*it, LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);

                //assemble the elemental contribution
                Assemble(A, b, LHS_Contribution, RHS_Contribution, EquationId, lock_array);
            }
        }

        double stop_prod = omp_get_wtime();
        if (this->GetEchoLevel() > 2 && r_model_part.GetCommunicator().MyPID() == 0)
            std::cout << "time: " << stop_prod - start_prod << std::endl;

        for (int i = 0; i < A_size; i++)
            omp_destroy_lock(&lock_array[i]);

        if( this->GetEchoLevel() > 2 && r_model_part.GetCommunicator().MyPID() == 0)
            KRATOS_WATCH( "finished parallel building" )

#endif
		mTotalBuildTime += Timer::GetTime() - time_begin;

#ifdef STATIC_GENERAL_BUILDER_AND_SOLVER_OPTIMIZE_SOLUTION
		mDoFactorization = true;
#endif // STATIC_GENERAL_BUILDER_AND_SOLVER_OPTIMIZE_SOLUTION
    }

    void BuildLHS(typename TSchemeType::Pointer pScheme, ModelPart& r_model_part, TSystemMatrixType& A)
    {
		double time_begin = Timer::GetTime();

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

		mTotalBuildTime += Timer::GetTime() - time_begin;

#ifdef STATIC_GENERAL_BUILDER_AND_SOLVER_OPTIMIZE_SOLUTION
		mDoFactorization = true;
#endif // STATIC_GENERAL_BUILDER_AND_SOLVER_OPTIMIZE_SOLUTION
    }

    void BuildLHS_CompleteOnFreeRows(typename TSchemeType::Pointer pScheme, ModelPart& r_model_part, TSystemMatrixType& A)
    {
		double time_begin = Timer::GetTime();

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

		mTotalBuildTime += Timer::GetTime() - time_begin;

#ifdef STATIC_GENERAL_BUILDER_AND_SOLVER_OPTIMIZE_SOLUTION
		mDoFactorization = true;
#endif // STATIC_GENERAL_BUILDER_AND_SOLVER_OPTIMIZE_SOLUTION
    }

    void SystemSolve(TSystemMatrixType& A, TSystemVectorType& Dx, TSystemVectorType& b)
    {
		double time_begin = Timer::GetTime();

        double norm_b;
        if (TSparseSpace::Size(b) != 0)
            norm_b = TSparseSpace::TwoNorm(b);
        else
            norm_b = 0.00;

        if (norm_b != 0.00)
        {
#ifdef STATIC_GENERAL_BUILDER_AND_SOLVER_OPTIMIZE_SOLUTION
			if(mDoFactorization) {
				BaseType::mpLinearSystemSolver->InitializeSolutionStep(A, Dx, b);
				BaseType::mpLinearSystemSolver->PerformSolutionStep(A, Dx, b);
				mDoFactorization = false;
			}
			else {
				BaseType::mpLinearSystemSolver->PerformSolutionStep(A, Dx, b);
			}
#else
			BaseType::mpLinearSystemSolver->Solve(A, Dx, b);
#endif // STATIC_GENERAL_BUILDER_AND_SOLVER_OPTIMIZE_SOLUTION
        }
        else
            TSparseSpace::SetToZero(Dx);

		mTotalSolutionTime += Timer::GetTime() - time_begin;
    }

    void SystemSolveWithPhysics(TSystemMatrixType& A, TSystemVectorType& Dx, TSystemVectorType& b, ModelPart& r_model_part)
    {
		double time_begin = Timer::GetTime();

        double norm_b;
        if (TSparseSpace::Size(b) != 0)
            norm_b = TSparseSpace::TwoNorm(b);
        else
            norm_b = 0.00;

        if (norm_b != 0.00)
        {
            if(BaseType::mpLinearSystemSolver->AdditionalPhysicalDataIsNeeded() )
                BaseType::mpLinearSystemSolver->ProvideAdditionalData(A, Dx, b, BaseType::mDofSet, r_model_part);

#ifdef STATIC_GENERAL_BUILDER_AND_SOLVER_OPTIMIZE_SOLUTION
			if(mDoFactorization) {
				BaseType::mpLinearSystemSolver->InitializeSolutionStep(A, Dx, b);
				BaseType::mpLinearSystemSolver->PerformSolutionStep(A, Dx, b);
				mDoFactorization = false;
			}
			else {
				BaseType::mpLinearSystemSolver->PerformSolutionStep(A, Dx, b);
			}
#else
			BaseType::mpLinearSystemSolver->Solve(A, Dx, b);
#endif // STATIC_GENERAL_BUILDER_AND_SOLVER_OPTIMIZE_SOLUTION

        }
        else
        {
            TSparseSpace::SetToZero(Dx);
			if( this->GetEchoLevel() > 0 && r_model_part.GetCommunicator().MyPID() == 0)
				std::cout << "ATTENTION! setting the RHS to zero!" << std::endl;
        }

		mTotalSolutionTime += Timer::GetTime() - time_begin;
    }

    void BuildAndSolve(typename TSchemeType::Pointer pScheme, ModelPart& r_model_part, TSystemMatrixType& A, TSystemVectorType& Dx, TSystemVectorType& b)
    {
        bool CalculateReactionsFlag = BaseType::mCalculateReactionsFlag;
		BaseType::mCalculateReactionsFlag = false;
        this->Build(pScheme, r_model_part, A, b);
		BaseType::mCalculateReactionsFlag = CalculateReactionsFlag;
        SystemSolveWithPhysics(A, Dx, b, r_model_part);
    }

    void BuildRHSAndSolve(typename TSchemeType::Pointer pScheme, ModelPart& r_model_part, TSystemMatrixType& A, TSystemVectorType& Dx, TSystemVectorType& b)
    {
        bool CalculateReactionsFlag = BaseType::mCalculateReactionsFlag;
		BaseType::mCalculateReactionsFlag = false;
		this->BuildRHS(pScheme, r_model_part, b);
		BaseType::mCalculateReactionsFlag = CalculateReactionsFlag;
        SystemSolve(A, Dx, b);
    }

    void BuildRHS(typename TSchemeType::Pointer pScheme, ModelPart& r_model_part, TSystemVectorType& b)
    {
		double time_begin = Timer::GetTime();

        if (!pScheme)
            KRATOS_THROW_ERROR( std::runtime_error, "No scheme provided!", "" )

        //getting the elements from the model
        ElementsArrayType& pElements = r_model_part.Elements();

        //getting the array of the conditions
        ConditionsArrayType& ConditionsArray = r_model_part.Conditions();

         //resetting to zero the vector of reactions
        TSparseSpace::SetToZero(*(BaseType::mpReactionsVector));


        // assemble all elements
#ifndef STATIC_GENERAL_BUILDER_AND_SOLVER_USES_OPENMP

        //contributions to the system
        LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0, 0);
        LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

		ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

        //vector containing the localization in the system of the different terms
        Element::EquationIdVectorType EquationId;

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

#else
        //creating an array of lock variables of the size of the system vector
        /*std::vector< omp_lock_t > lock_array(b.size());*/
	    int total_size = b.size();
	    if (BaseType::mCalculateReactionsFlag)
			total_size += (*BaseType::mpReactionsVector).size();
	    std::vector< omp_lock_t > lock_array(total_size);

        //int b_size = b.size();
        for (int i = 0; i < total_size; i++)
            omp_init_lock(&lock_array[i]);

        //create a partition of the element array
        int number_of_threads = omp_get_max_threads();

        vector<unsigned int> element_partition;
        CreatePartition(number_of_threads, pElements.size(), element_partition);


        if( this->GetEchoLevel() > 2 && r_model_part.GetCommunicator().MyPID() == 0)
        {
            KRATOS_WATCH( number_of_threads )
            KRATOS_WATCH( element_partition )
        }

        double start_prod = omp_get_wtime();

        #pragma omp parallel for
        for (int k = 0; k < number_of_threads; k++)
        {
            //contributions to the system
            LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

            //vector containing the localization in the system of the different
            //terms
            Element::EquationIdVectorType EquationId;

            ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

            typename ElementsArrayType::ptr_iterator it_begin = pElements.ptr_begin() + element_partition[k];
            typename ElementsArrayType::ptr_iterator it_end = pElements.ptr_begin() + element_partition[k + 1];

            // assemble all elements
            for (typename ElementsArrayType::ptr_iterator it = it_begin; it != it_end; ++it)
            {
                //calculate elemental contribution
				#pragma omp critical
				{
					pScheme->Calculate_RHS_Contribution(*it, RHS_Contribution, EquationId, CurrentProcessInfo);
				}

                //assemble the elemental contribution
                AssembleRHS(b, RHS_Contribution, EquationId, lock_array);

                // clean local elemental memory
                pScheme->CleanMemory(*it);
            }
        }

        vector<unsigned int> condition_partition;
        CreatePartition(number_of_threads, ConditionsArray.size(), condition_partition);

        #pragma omp parallel for
        for (int k = 0; k < number_of_threads; k++)
        {
            //contributions to the system
            LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

            Condition::EquationIdVectorType EquationId;

            ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

            typename ConditionsArrayType::ptr_iterator it_begin = ConditionsArray.ptr_begin() + condition_partition[k];
            typename ConditionsArrayType::ptr_iterator it_end = ConditionsArray.ptr_begin() + condition_partition[k + 1];

            // assemble all elements
            for (typename ConditionsArrayType::ptr_iterator it = it_begin; it != it_end; ++it)
            {
                //calculate elemental contribution
                pScheme->Condition_Calculate_RHS_Contribution(*it, RHS_Contribution, EquationId, CurrentProcessInfo);

                //assemble the elemental contribution
                AssembleRHS(b, RHS_Contribution, EquationId, lock_array);

            }
        }


        double stop_prod = omp_get_wtime();
        if (this->GetEchoLevel() > 2 && r_model_part.GetCommunicator().MyPID() == 0)
            std::cout << "time: " << stop_prod - start_prod << std::endl;

        for (int i = 0; i < total_size; i++)
            omp_destroy_lock(&lock_array[i]);

        if( this->GetEchoLevel() > 2 && r_model_part.GetCommunicator().MyPID() == 0)
        {
            KRATOS_WATCH( "finished parallel building" )
        }

#endif

		mTotalBuildTime += Timer::GetTime() - time_begin;
    }

    void SetUpDofSet(typename TSchemeType::Pointer pScheme, ModelPart& r_model_part)
    {
        KRATOS_TRY;

        if( this->GetEchoLevel() > 2 && r_model_part.GetCommunicator().MyPID() == 0)
            std::cout << " Setting up the dofs " << std::endl;

        ElementsArrayType& pElements = r_model_part.Elements();

        Element::DofsVectorType ElementalDofList;

        ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

        DofsArrayType Doftemp;
        BaseType::mDofSet = DofsArrayType();

        for (typename ElementsArrayType::ptr_iterator it = pElements.ptr_begin(); it != pElements.ptr_end(); ++it)
        {
            pScheme->GetElementalDofList(*it, ElementalDofList, CurrentProcessInfo);
            for (typename Element::DofsVectorType::iterator i = ElementalDofList.begin(); i != ElementalDofList.end(); ++i)
                Doftemp.push_back(i->get());
        }

        ConditionsArrayType& pConditions = r_model_part.Conditions();
        for (typename ConditionsArrayType::ptr_iterator it = pConditions.ptr_begin(); it != pConditions.ptr_end(); ++it)
        {
            pScheme->GetConditionDofList(*it, ElementalDofList, CurrentProcessInfo);
            for (typename Element::DofsVectorType::iterator i = ElementalDofList.begin(); i != ElementalDofList.end(); ++i)
                Doftemp.push_back(i->get());
        }

        Doftemp.Unique();
        BaseType::mDofSet = Doftemp;

        if (BaseType::mDofSet.size() == 0)
            KRATOS_THROW_ERROR( std::logic_error, "No degrees of freedom!", "" )

        BaseType::mDofSetIsInitialized = true;
        if( this->GetEchoLevel() > 2 && r_model_part.GetCommunicator().MyPID() == 0)
            std::cout << "finished setting up the dofs" << std::endl;

        KRATOS_CATCH( "" )
    }

    void SetUpSystem(ModelPart& r_model_part)
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

    void ResizeAndInitializeVectors( typename TSchemeType::Pointer pScheme,TSystemMatrixPointerType& pA,
                                    TSystemVectorPointerType& pDx,
                                    TSystemVectorPointerType& pb,
                                    ModelPart& rModelPart)
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
            ConstructMatrixStructure(pScheme, A, rModelPart.Elements(), rModelPart.Conditions(), rModelPart.GetProcessInfo());
        }
        else
        {
            if (A.size1() != BaseType::mEquationSystemSize || A.size2() != BaseType::mEquationSystemSize)
            {
                KRATOS_WATCH( "it should not come here!!!!!!!! ... this is SLOW" )
                A.resize(BaseType::mEquationSystemSize, BaseType::mEquationSystemSize, true);
                ConstructMatrixStructure(pScheme, A, rModelPart.Elements(), rModelPart.Conditions(), rModelPart.GetProcessInfo());
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

        KRATOS_CATCH( "" )
    }

    void InitializeSolutionStep(ModelPart& r_model_part, TSystemMatrixType& A, TSystemVectorType& Dx, TSystemVectorType& b)
    {
        mTotalBuildTime = 0.0;
        mTotalSolutionTime = 0.0;
    }

    void FinalizeSolutionStep(ModelPart& r_model_part, TSystemMatrixType& A, TSystemVectorType& Dx, TSystemVectorType& b)
    {
		mTotalBuildTime_accum += mTotalBuildTime;
		mTotalSolutionTime_accum += mTotalSolutionTime;
        if( this->GetEchoLevel() > 0 && r_model_part.GetCommunicator().MyPID() == 0)
            PromptStatistics(r_model_part.GetProcessInfo());
    }

    void CalculateReactions(typename TSchemeType::Pointer pScheme, ModelPart& r_model_part, TSystemMatrixType& A, TSystemVectorType& Dx, TSystemVectorType& b)
    {
        ////refresh RHS to have the correct reactions
        //this->BuildRHS(pScheme, r_model_part, b);

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
                (*it2)->GetSolutionStepReactionValue() = ReactionsVector[i];
            }
            num++;
        }
    }

    void Clear()
    {
        this->mDofSet = DofsArrayType();

        if (this->mpReactionsVector != NULL)
            TSparseSpace::Clear((this->mpReactionsVector));

        this->mpLinearSystemSolver->Clear();
    }

    virtual int Check(ModelPart& r_model_part)
    {
        return 0;
    }

protected:

    virtual void ConstructMatrixStructure( typename TSchemeType::Pointer pScheme,TSystemMatrixType& A,
                                          ModelPart& rModelPart)
    {
        std::size_t equation_size = A.size1();
        std::vector<std::vector<std::size_t> > indices(equation_size);

        Element::EquationIdVectorType ids(3, 0);
        for (typename ElementsContainerType::iterator i_element = rModelPart.ElementsBegin(); i_element != rModelPart.ElementsEnd(); i_element++)
        {

            pScheme->EquationId( *(i_element.base()) , ids, CurrentProcessInfo);
            for (std::size_t i = 0; i < ids.size(); i++)
                if (ids[i] < equation_size)
                {
                    std::vector<std::size_t>& row_indices = indices[ids[i]];
                    for (std::size_t j = 0; j < ids.size(); j++)
                        if (ids[j] < equation_size)
                        {
                            AddUnique(row_indices, ids[j]);
                        }
                }

        }

        for (typename ConditionsArrayType::iterator i_condition = rModelPart.ConditionsBegin(); i_condition != rModelPart.ConditionsEnd(); i_condition++)
        {

            pScheme->Condition_EquationId( *(i_condition.base()), ids, CurrentProcessInfo);
            for (std::size_t i = 0; i < ids.size(); i++)
                if (ids[i] < equation_size)
                {
                    std::vector<std::size_t>& row_indices = indices[ids[i]];
                    for (std::size_t j = 0; j < ids.size(); j++)
                        if (ids[j] < equation_size)
                        {
                            AddUnique(row_indices, ids[j]);
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
#ifndef STATIC_GENERAL_BUILDER_AND_SOLVER_USES_OPENMP
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

    void AssembleLHS(TSystemMatrixType& A,
                     LocalSystemMatrixType& LHS_Contribution,
                     Element::EquationIdVectorType& EquationId)
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

    void AssembleRHS(TSystemVectorType& b,
                     LocalSystemVectorType& RHS_Contribution,
                     Element::EquationIdVectorType& EquationId)
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

	void AssembleLHS_CompleteOnFreeRows(TSystemMatrixType& A,
                                        LocalSystemMatrixType& LHS_Contribution,
                                        Element::EquationIdVectorType& EquationId)
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

    inline void CreatePartition(unsigned int number_of_threads, const int number_of_rows, vector<unsigned int>& partitions)
    {
        partitions.resize(number_of_threads + 1);
        int partition_size = number_of_rows / number_of_threads;
        partitions[0] = 0;
        partitions[number_of_threads] = number_of_rows;
        for (unsigned int i = 1; i < number_of_threads; i++)
            partitions[i] = partitions[i - 1] + partition_size;
    }

private:

    double mTotalBuildTime;
    double mTotalSolutionTime;
	double mTotalBuildTime_accum;
    double mTotalSolutionTime_accum;
    unsigned int mNumberOfFreeDofs;

#ifdef STATIC_GENERAL_BUILDER_AND_SOLVER_OPTIMIZE_SOLUTION
	bool mDoFactorization;
#endif // STATIC_GENERAL_BUILDER_AND_SOLVER_OPTIMIZE_SOLUTION

private:

    inline void PromptStatistics(const ProcessInfo & rCurrentProcessInfo)
    {
        std::stringstream ss;
        ss << " Total Build Time:    \t" << mTotalBuildTime << " seconds" << std::endl;
        ss << " Total Solution Time: \t" << mTotalSolutionTime << " seconds" << std::endl;
        int nit = rCurrentProcessInfo[NL_ITERATION_NUMBER];
		double dnit = nit == 0.0 ? 1.0 : (double)nit;
        ss << " Number of Iterations:\t" << nit << std::endl;
        ss << " Mean Build Time:    \t" << mTotalBuildTime / dnit << " seconds" << std::endl;
        ss << " Mean Solution Time: \t" << mTotalSolutionTime / dnit << " seconds" << std::endl;
		ss << " Total Build Time from start: " << mTotalBuildTime_accum << " seconds" << std::endl;
		ss << " Total Solution Time from start: " << mTotalSolutionTime_accum << " seconds" << std::endl;
		ss << " Total Time from start: " << mTotalSolutionTime_accum + mTotalBuildTime_accum << " seconds" << std::endl;
        std::cout << ss.str();
    }

#ifdef STATIC_GENERAL_BUILDER_AND_SOLVER_USES_OPENMP

    void Assemble(TSystemMatrixType& A,
                  TSystemVectorType& b,
                  const LocalSystemMatrixType& LHS_Contribution,
                  const LocalSystemVectorType& RHS_Contribution,
                  Element::EquationIdVectorType& EquationId,
                  std::vector< omp_lock_t >& lock_array)
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
        }
    }

	void AssembleLHS(TSystemMatrixType& A,
                     const LocalSystemMatrixType& LHS_Contribution,
                     Element::EquationIdVectorType& EquationId,
                     std::vector< omp_lock_t >& lock_array)
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
        }
    }

	void AssembleRHS(TSystemVectorType& b,
                     const LocalSystemVectorType& RHS_Contribution,
                     Element::EquationIdVectorType& EquationId,
                     std::vector< omp_lock_t >& lock_array)
    {
        unsigned int local_size = RHS_Contribution.size();

        if (BaseType::mCalculateReactionsFlag == false)
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
            }
	    }
        else
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

};

}

#endif /* KRATOS_STATIC_GENERAL_BUILDER_AND_SOLVER  defined */

