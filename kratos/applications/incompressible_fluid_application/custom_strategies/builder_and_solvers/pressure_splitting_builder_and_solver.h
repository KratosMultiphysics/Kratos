/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/

/*
 * File:   pressure_splitting_builder_and_solver.h
 * Author: jcotela
 *
 * Created on 18 February 2010, 15:07
 */

#if !defined(KRATOS_PRESSURE_SPLITTING_BUILDER_AND_SOLVER_H)
#define	KRATOS_PRESSURE_SPLITTING_BUILDER_AND_SOLVER_H

/* System includes */
#include <set>
#ifdef _OPENMP
#include <fstream>
#include <omp.h>
#endif

/* External includes */
#include <boost/smart_ptr.hpp>
#include "utilities/timer.h"
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

/* Project includes */
#include "includes/define.h"
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"

namespace Kratos
{

    // PressureSplittingBuilderAndSolver class
    template< class TSparseSpace,
	class TDenseSpace , //= DenseSpace<double>,
	class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
	>
    class PressureSplittingBuilderAndSolver:
        public BuilderAndSolver< TSparseSpace,TDenseSpace,TLinearSolver >
    {
        /* Please note that this builder and solver splits the system matrix
         * in 4 blocks and stores them separately. The system matrix A used as
         * input (to preserve the same function signstures as other builder and
         * solvers) is used here to store the matrix of the pressure equation.
         */
    public:

        KRATOS_CLASS_POINTER_DEFINITION( PressureSplittingBuilderAndSolver );

        // Type Definitions
        typedef BuilderAndSolver<TSparseSpace,TDenseSpace, TLinearSolver> BaseType;

        typedef typename BaseType::TSchemeType TSchemeType;

        typedef typename BaseType::TDataType TDataType;

        typedef Dof<TDataType> TDofType;
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

        typedef std::size_t KeyType; // For Dof->GetVariable().Key()

        typedef boost::numeric::ublas::vector<int> IndexVector;
        typedef std::vector<unsigned int> PartitionVector;

        // UBLAS matrix access types
        typedef typename TSystemMatrixType::iterator1 OuterIt;
        typedef typename TSystemMatrixType::iterator2 InnerIt;
        typedef typename boost::numeric::ublas::matrix_row< TSystemMatrixType > RowType;

        // Life Cycle

        /* Constructor */
        PressureSplittingBuilderAndSolver(
                typename TLinearSolver::Pointer pNewLinearSystemSolver,
                unsigned int RebuildLevel, // If > 0 Do not re-check the shape of all matrices each step (see definition in private atribute variables)
                unsigned int VelocityCorrection, // If > 0, explicitly solve the velocity to be divergence-free at each step
                bool UseInexactNewton, // If true, dynamically set the linear iterative solver tolerance for the pressure system
                double NonLinearTol = 1e-3, // Only used if InexactNewton == true, otherwise the solver will use it's own tolerance
                double MaxTolFactor = 0.1,
                double Gamma = 0.9):
            BuilderAndSolver< TSparseSpace,TDenseSpace,TLinearSolver >(pNewLinearSystemSolver),
            mRebuildLevel(RebuildLevel),
            mInitializedMatrices(false),
            mVelocityCorrection(VelocityCorrection),
            mInexactNewton(UseInexactNewton),
            mMaxTolFactor(MaxTolFactor),
            mGamma(Gamma),
            mFirstIteration(true),
//            mVelTolFactor(0),
            mPressTolFactor(0),
//            mLastVelRHSNorm(0),
            mLastPressRHSNorm(0)
        {
            mSmallTol = 0.5*NonLinearTol;
        }

        /* Destructor */
        virtual ~PressureSplittingBuilderAndSolver() {}

        // Operations

        /* Build System */
        void Build(
                typename TSchemeType::Pointer pScheme,
                ModelPart& rModelPart,
                TSystemMatrixType& A,
                TSystemVectorType& b)
        {
            KRATOS_TRY
            if (!pScheme)
                KRATOS_ERROR(std::runtime_error, "No scheme provided!", "");

            // Get elements and conditions
            ElementsArrayType& rElements = rModelPart.Elements();
            ConditionsArrayType& rConditions = rModelPart.Conditions();

            // resetting to zero the vector of reactions
            TSparseSpace::SetToZero( *(BaseType::mpReactionsVector) );

            // Reset internally stored matrices
            TSparseSpace::SetToZero( *mpS );
            TSparseSpace::SetToZero( *mpG );
            TSparseSpace::SetToZero( *mpD );
            TSparseSpace::SetToZero( *mpL );

            // Define contributions to the system
            LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0,0);
            LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

            // Will store the position of each Dof in the system
            Element::EquationIdVectorType EquationIds;

            ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();

#ifndef _OPENMP
            // Assemble contributions from elements
            for( typename ElementsArrayType::ptr_iterator pElem = rElements.ptr_begin();
                    pElem != rElements.ptr_end(); pElem++)
            {
                // Get Elemental Contributions
                pScheme->CalculateSystemContributions(*pElem,LHS_Contribution,
                        RHS_Contribution,EquationIds,CurrentProcessInfo);

                //assemble the elemental contribution
                Assemble(A,b,LHS_Contribution,RHS_Contribution,EquationIds);

                // clean local elemental memory
                pScheme->CleanMemory(*pElem);
            }

            LHS_Contribution.resize(0,0,false);
            RHS_Contribution.resize(0,false);

            // assemble all conditions
            for ( typename ConditionsArrayType::ptr_iterator pCond = rConditions.ptr_begin();
                    pCond != rConditions.ptr_end(); pCond++)
            {
                // Get condition Contributions
                pScheme->Condition_CalculateSystemContributions(*pCond,LHS_Contribution,
                        RHS_Contribution,EquationIds,CurrentProcessInfo);

                // Assemble condition contribution
                Assemble(A,b,LHS_Contribution,RHS_Contribution,EquationIds);
            }

#else
            //creating an array of lock variables of the size of the system matrix
            int TotalSize = BaseType::mEquationSystemSize;
            std::vector< omp_lock_t > lock_array(TotalSize);

            for(unsigned int i = 0; i<TotalSize; i++)
                omp_init_lock( &lock_array[i] );

            //create a partition of the element array
            unsigned int NumThreads = GetNumThreads();
            PartitionVector ElementPartition;

            DivideInPartitions(rElements.size(),ElementPartition,NumThreads);

            double start_prod = omp_get_wtime();

            #pragma omp parallel for
            for( int k=0; k < NumThreads; k++)
            {
                //contributions to the system
                LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0,0);
                LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

                //vector containing the localization in the system of the different
                //terms
                Element::EquationIdVectorType EquationId;
                ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();
                typename ElementsArrayType::ptr_iterator ElemBegin = rElements.ptr_begin()+ElementPartition[k];
                typename ElementsArrayType::ptr_iterator ElemEnd = rElements.ptr_begin()+ElementPartition[k+1];

                // assemble all elements
                for (typename ElementsArrayType::ptr_iterator pElem = ElemBegin; pElem != ElemEnd; pElem++)
                {
                    //calculate elemental contribution
                    pScheme->CalculateSystemContributions(*pElem,LHS_Contribution,RHS_Contribution,EquationId,CurrentProcessInfo);

                    //assemble the elemental contribution
                    Parallel_Assemble(A,b,LHS_Contribution,RHS_Contribution,EquationId,lock_array);

                    // clean local elemental memory
                    pScheme->CleanMemory(*pElem);
                }
            }

            PartitionVector ConditionPartition;

            DivideInPartitions(rConditions.size(),ConditionPartition,NumThreads);

            #pragma omp parallel
            {
                unsigned int k = ThisThread();

                //contributions to the system
                LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0,0);
                LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

                Condition::EquationIdVectorType EquationId;

                ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();

                typename ConditionsArrayType::ptr_iterator CondBegin = rConditions.ptr_begin()+ConditionPartition[k];
                typename ConditionsArrayType::ptr_iterator CondEnd = rConditions.ptr_begin()+ConditionPartition[k+1];

                // assemble all conditions
                for (typename ConditionsArrayType::ptr_iterator pCond = CondBegin; pCond != CondEnd; pCond++)
                {
                    //calculate condition contribution
                    pScheme->Condition_CalculateSystemContributions(*pCond,LHS_Contribution,RHS_Contribution,EquationId,CurrentProcessInfo);

                    //assemble condition contribution
                    Parallel_Assemble(A,b,LHS_Contribution,RHS_Contribution,EquationId,lock_array);
                }
            }


            double stop_prod = omp_get_wtime();
            std::cout << "time: " << stop_prod - start_prod << std::endl;

            for(int i = 0; i< TotalSize; i++)
                omp_destroy_lock(&lock_array[i]);

            KRATOS_WATCH("finished parallel building");
#endif

            /* Build the pressure system matrix */
            if (mRebuildLevel == 0 || mInitializedMatrices == false)
            {
                CalculateSystemMatrix(A);
                mInitializedMatrices = true;
            }

            KRATOS_CATCH("");
        }

        /* Build LHS only */
        void BuildLHS(
                typename TSchemeType::Pointer pScheme,
                ModelPart& rModelPart,
                TSystemMatrixType& A)
        {
            KRATOS_TRY;
            if (!pScheme)
                KRATOS_ERROR(std::runtime_error, "No scheme provided!", "");

            // Get elements and conditions
            ElementsArrayType& rElements = rModelPart.Elements();
            ConditionsArrayType& rConditions = rModelPart.Conditions();

            // resetting to zero the vector of reactions
            TSparseSpace::SetToZero(*(BaseType::mpReactionsVector));

            // Reset internally stored matrices
            TSparseSpace::SetToZero( *mpS );
            TSparseSpace::SetToZero( *mpG );
            TSparseSpace::SetToZero( *mpD );
            TSparseSpace::SetToZero( *mpL );

            // Define contributions to the system
            LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0, 0);

            // Will store the position of each Dof in the system
            Element::EquationIdVectorType EquationIds;

            ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();

//#ifndef _OPENMP
            // Assemble contributions from elements
            for (typename ElementsArrayType::ptr_iterator pElem = rElements.ptr_begin();
                    pElem != rElements.ptr_end(); pElem++)
            {
                // Get Elemental Contributions
                pScheme->Calculate_LHS_Contribution(*pElem, LHS_Contribution,
                        EquationIds, CurrentProcessInfo);

                //assemble the elemental contribution
                AssembleLHS(A, LHS_Contribution, EquationIds);

                // clean local elemental memory
                pScheme->CleanMemory(*pElem);
            }

            LHS_Contribution.resize(0, 0, false);

            // assemble all conditions
            for ( typename ConditionsArrayType::ptr_iterator pCond = rConditions.ptr_begin();
                    pCond != rConditions.ptr_end(); pCond++)
            {
                // Get Elemental Contributions
                pScheme->Condition_Calculate_LHS_Contribution(*pCond, LHS_Contribution,
                        EquationIds, CurrentProcessInfo);

                //assemble the elemental contribution
                AssembleLHS(A, LHS_Contribution, EquationIds);
            }

            /* Build the pressure system matrix */
            if (mRebuildLevel == 0 || mInitializedMatrices == false)
            {
                CalculateSystemMatrix(A);
                mInitializedMatrices = true;
            }
//#else
//            // DO MP VERSION HERE
//#endif

            KRATOS_CATCH("");
        }

        /* Build LHS with extra columns, including fixed Dofs */
        void BuildLHS_CompleteOnFreeRows(
                typename TSchemeType::Pointer pScheme,
                ModelPart& rModelPart,
                TSystemMatrixType& A)
        {
            KRATOS_TRY;
            if (!pScheme)
                KRATOS_ERROR(std::runtime_error, "No scheme provided!", "");

            // Get elements and conditions
            ElementsArrayType& rElements = rModelPart.Elements();
            ConditionsArrayType& rConditions = rModelPart.Conditions();

            // resetting to zero the vector of reactions
            TSparseSpace::SetToZero(*(BaseType::mpReactionsVector));

            // Reset internally stored matrices
            TSparseSpace::SetToZero( *mpS );
            TSparseSpace::SetToZero( *mpG );
            TSparseSpace::SetToZero( *mpD );
            TSparseSpace::SetToZero( *mpL );

            // Define contributions to the system
            LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0, 0);

            // Will store the position of each Dof in the system
            Element::EquationIdVectorType EquationIds;

            ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();

//#ifndef _OPENMP
            // Assemble contributions from elements
            for (typename ElementsArrayType::ptr_iterator pElem = rElements.ptr_begin();
                    pElem != rElements.ptr_end(); pElem++)
            {
                // Get Elemental Contributions
                pScheme->Calculate_LHS_Contribution(*pElem, LHS_Contribution,
                        EquationIds, CurrentProcessInfo);

                //assemble the elemental contribution
                AssembleLHS_CompleteOnFreeRows(A, LHS_Contribution, EquationIds);

                // clean local elemental memory
                pScheme->CleanMemory(*pElem);
            }

            LHS_Contribution.resize(0, 0, false);

            // assemble all conditions
            for (typename ConditionsArrayType::ptr_iterator pCond = rConditions.ptr_begin();
                    pCond != rConditions.ptr_end(); pCond++)
            {
                // Get Elemental Contributions
                pScheme->Condition_Calculate_LHS_Contribution(*pCond, LHS_Contribution,
                        EquationIds, CurrentProcessInfo);

                //assemble the elemental contribution
                AssembleLHS_CompleteOnFreeRows(A, LHS_Contribution, EquationIds);
            }

            /* Build the pressure system matrix */
            if (mRebuildLevel == 0 || mInitializedMatrices == false)
            {
                CalculateSystemMatrix(A);
                mInitializedMatrices = true;
            }
//#else
//            // DO MP VERSION HERE
//#endif

            KRATOS_CATCH("");
        }

        /* Call Solver */
        void SystemSolve(
                TSystemMatrixType& A,
                TSystemVectorType& Dx,
                TSystemVectorType& b)
        {
            KRATOS_TRY;

            double NormB;

            if (TSparseSpace::Size(b) != 0)
                NormB = TSparseSpace::TwoNorm(b);
            else
                NormB = 0;

            if (NormB != 0.0)
            {
                /* Solve current iteration in several steps */

                // Initialize required variables and pointers
                TSystemVectorPointerType pDVel(new TSystemVectorType(mVelFreeDofs));
                TSystemVectorType& rDVel = *pDVel;
                TSystemVectorPointerType pVelRHS(new TSystemVectorType(mVelFreeDofs));
                TSystemVectorType& rVelRHS = *pVelRHS;
                for (unsigned int i = 0; i < mVelFreeDofs; i++)
                {
                    rDVel[i] = 0.0;
                    rVelRHS[i] = b[i];
                }

                TSystemVectorPointerType pDPress(new TSystemVectorType(mPressFreeDofs));
                TSystemVectorType& rDPress = *pDPress;
                TSystemVectorPointerType pPressRHS(new TSystemVectorType(mPressFreeDofs));
                TSystemVectorType& rPressRHS = *pPressRHS;
                for (unsigned int i = 0, j=mVelFreeDofs; i < mPressFreeDofs; i++,j++)
                {
                    rDPress[i] = 0.0;
                    rPressRHS[i] = b[j];
                }

                TSystemMatrixType& rS = *mpS; // Create a reference to the velocity system matrix

                // 1. Compute intermediate velocity
//                axpy_prod(*mpG, -rDPress, rVelRHS, false);

//                if (mInexactNewton == true)
//                {
//                    double VelRHSNorm = TSparseSpace::TwoNorm(rVelRHS);
//                    if(mFirstIteration == true) {
//                        SetInitialTolerance(VelRHSNorm,mVelTolFactor);
//                    } else {
//                        UpdateTolerance(mLastVelRHSNorm,VelRHSNorm,mVelTolFactor);
//                    }
//                    mLastVelRHSNorm = VelRHSNorm;
//                }

                BaseType::mpLinearSystemSolver->Solve(rS, rDVel, rVelRHS);

                // 2. Compute Pressure Variation
#ifndef _OPENMP
                axpy_prod(*mpD, -rDVel, rPressRHS, false);
#else
                Parallel_ProductAdd(*mpD, -rDVel, rPressRHS);
#endif

                if (mInexactNewton == true)
                {
                    double PressRHSNorm = TSparseSpace::TwoNorm(rPressRHS);
                    if(mFirstIteration == true) {
                        SetInitialTolerance(PressRHSNorm,mPressTolFactor);
                    } else {
                        UpdateTolerance(mLastPressRHSNorm,PressRHSNorm,mPressTolFactor);
                    }
                    mLastPressRHSNorm = PressRHSNorm;
                }

                BaseType::mpLinearSystemSolver->Solve(A, rDPress, rPressRHS);

                // 3. Determine End of Step velocity
                if ( mVelocityCorrection > 0)
                {
                    TSparseSpace::Mult(*mpG, rDPress,rVelRHS);

                    if ( mVelocityCorrection == 1 ) {
                        for (unsigned int m = 0; m < mVelFreeDofs; m++)
                            rDVel[m] -= (*mpIDiagS)[m]*rVelRHS[m];
                    } else if ( mVelocityCorrection == 2 ) {
                        TSystemVectorPointerType pVelUpdate(new TSystemVectorType(mVelFreeDofs));
                        TSystemVectorType& rVelUpdate = *pVelUpdate;
                        for (unsigned int i = 0; i < mVelFreeDofs; i++) rVelUpdate[i] = 0.0;

                        BaseType::mpLinearSystemSolver->Solve(rS, rVelUpdate, rVelRHS);
                        noalias(rDVel) -= rVelUpdate;
                    }
                }

                // Preconditioner
//                A = *mpL - A;
//                noalias(rPressRHS) = prod(A,rPressRHS);

                // Copy the solution to output variable
                for (unsigned int i = 0; i < mVelFreeDofs; i++) Dx[i] = rDVel[i];
                for (unsigned int i = 0; i < mPressFreeDofs; i++) Dx[mVelFreeDofs + i] = rDPress[i];

                if (mFirstIteration == true) mFirstIteration = false;
            }
            else
                TSparseSpace::SetToZero(Dx);

            //prints informations about the current time
            if (this->GetEchoLevel() > 1)
            {
                std::cout << *(BaseType::mpLinearSystemSolver) << std::endl;
            }

            KRATOS_CATCH("");
        }

        /* Build and Solve sytem */
        void BuildAndSolve(
                typename TSchemeType::Pointer pScheme,
                ModelPart& rModelPart,
                TSystemMatrixType& A,
                TSystemVectorType& Dx,
                TSystemVectorType& b)
        {
            KRATOS_TRY

#ifdef _OPENMP
            std::ofstream results;
            results.open("results.csv", std::ios::out | std::ios::app );

            double t0 = omp_get_wtime();
#endif

            Timer::Start("Build");

            Build(pScheme, rModelPart, A, b);

            Timer::Stop("Build");
#ifdef _OPENMP
            double t1 = omp_get_wtime();
#endif

//        //does nothing...dirichlet conditions are naturally dealt with in defining the residual
//        ApplyDirichletConditions(pScheme,rModelPart,A,Dx,b);

            if (this->GetEchoLevel() == 3) {
                std::cout << "before the solution of the system" << std::endl;
                std::cout << "System Matrix = " << A << std::endl;
                std::cout << "unknowns vector = " << Dx << std::endl;
                std::cout << "RHS vector = " << b << std::endl;
            }

            Timer::Start("Solve");

            SystemSolve(A, Dx, b);

            Timer::Stop("Solve");
#ifdef _OPENMP
            double t2 = omp_get_wtime();

            if (this->GetEchoLevel() == 3) {
                std::cout << "after the solution of the system" << std::endl;
                std::cout << "System Matrix = " << A << std::endl;
                std::cout << "unknowns vector = " << Dx << std::endl;
                std::cout << "RHS vector = " << b << std::endl;
            }

            results << (t1-t0) << " " << (t2-t1) << "\n";
            results.close();
#endif
            KRATOS_CATCH("");
        }

        /* Solve System for updated RHS */
        void BuildRHSAndSolve(
                typename TSchemeType::Pointer pScheme,
                ModelPart& rModelPart,
                TSystemMatrixType& A,
                TSystemVectorType& Dx,
                TSystemVectorType& b)
        {
            KRATOS_TRY

            BuildRHS(pScheme,rModelPart,b);
            SystemSolve(A,Dx,b);

            KRATOS_CATCH("");
        }

        /* Build RHS only */
        void BuildRHS(
                typename TSchemeType::Pointer pScheme,
                ModelPart& rModelPart,
                TSystemVectorType& b)
        {
            KRATOS_TRY;
            if (!pScheme)
                KRATOS_ERROR(std::runtime_error, "No scheme provided!", "");

            // Get elements and conditions
            ElementsArrayType& rElements = rModelPart.Elements();
            ConditionsArrayType& rConditions = rModelPart.Conditions();

            // resetting to zero the vector of reactions
            TSparseSpace::SetToZero(*(BaseType::mpReactionsVector));

            // Define contributions to the system
            LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

            // Will store the position of each Dof in the system
            Element::EquationIdVectorType EquationIds;

            ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();

//#ifndef _OPENMP

            // Assemble contributions from elements
            for (typename ElementsArrayType::ptr_iterator pElem = rElements.ptr_begin();
                    pElem != rElements.ptr_end(); pElem++)
            {
                // Get Elemental Contributions
                pScheme->Calculate_RHS_Contribution(*pElem, RHS_Contribution,
                        EquationIds, CurrentProcessInfo);

                //assemble the elemental contribution
                AssembleRHS(b, RHS_Contribution, EquationIds);

                // clean local elemental memory
                pScheme->CleanMemory(*pElem);
            }

            RHS_Contribution.resize(0, false);

            // assemble all conditions
            for (typename ConditionsArrayType::ptr_iterator pCond = rConditions.ptr_begin();
                    pCond != rConditions.ptr_end(); pCond++)
            {
                // Get Elemental Contributions
                pScheme->Condition_Calculate_RHS_Contribution(*pCond,
                        RHS_Contribution, EquationIds, CurrentProcessInfo);

                //assemble the elemental contribution
                AssembleRHS(b, RHS_Contribution, EquationIds);
            }

//#else
//            // DO MP VERSION HERE
//#endif

            KRATOS_CATCH("");
        }

        /* Identify Dofs and store pointers to them */
        /* Single list version */
        void SetUpDofSet(
                typename TSchemeType::Pointer pScheme,
                ModelPart& rModelPart)
        {
            KRATOS_TRY;

            KRATOS_WATCH("setting up the dofs");
            //Gets the array of elements from the modeler
            ElementsArrayType& pElements = rModelPart.Elements();
            ConditionsArrayType& pConditions = rModelPart.Conditions();

            BaseType::mDofSet = DofsArrayType();

            ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();

            unsigned int NumThreads = GetNumThreads();
            PartitionVector ElemPartition;
            PartitionVector CondPartition;

            DivideInPartitions(pElements.size(),ElemPartition,NumThreads);
            DivideInPartitions(pConditions.size(),CondPartition,NumThreads);

            std::vector< DofsArrayType > DofContainer(NumThreads);

            #pragma omp parallel
            {
                unsigned int k = ThisThread();
                Element::DofsVectorType ElementalDofList;

                for (typename ElementsArrayType::ptr_iterator it = pElements.ptr_begin() + ElemPartition[k];
                        it != pElements.ptr_begin() + ElemPartition[k+1]; ++it) {
                    // gets list of Dof involved on every element
                    pScheme->GetElementalDofList(*it, ElementalDofList, CurrentProcessInfo);

                    for (typename Element::DofsVectorType::iterator i = ElementalDofList.begin();
                            i != ElementalDofList.end(); ++i) {
                        DofContainer[k].push_back(*i);
                    }
                }

                //taking into account conditions
                for (typename ConditionsArrayType::ptr_iterator it = pConditions.ptr_begin() + CondPartition[k];
                        it != pConditions.ptr_begin() + CondPartition[k+1]; ++it) {
                    // gets list of Dof involved on every element
                    pScheme->GetConditionDofList(*it, ElementalDofList, CurrentProcessInfo);

                    for (typename Element::DofsVectorType::iterator i = ElementalDofList.begin();
                            i != ElementalDofList.end(); ++i) {
                        DofContainer[k].push_back(*i);
                    }
                }

                // Remove duplicates in the partial list
                // (try to do as much work as possible in the parallel region)
                DofContainer[k].Unique();
            }

            // Generate a single list
            for (unsigned int k = 0; k < NumThreads ; k++)
                for( typename DofsArrayType::ptr_iterator itDof = DofContainer[k].ptr_begin();
                        itDof != DofContainer[k].ptr_end(); itDof++) {
                    BaseType::mDofSet.push_back(*itDof);
                }

            BaseType::mDofSet.Unique();

            //throws an execption if there are no Degrees of freedom involved in the analysis
            if (BaseType::mDofSet.size() == 0)
                KRATOS_ERROR(std::logic_error, "No degrees of freedom!", "");

            BaseType::mDofSetIsInitialized = true;

            KRATOS_CATCH("");
        }

        /* Organise Dofs, separating fixed and free nodes */
        void SetUpSystem(ModelPart& rModelPart)
        {
            KRATOS_TRY;

            unsigned int FreeId = 0;
            unsigned int FixedId = BaseType::mDofSet.size();

            for (typename DofsArrayType::iterator itDof = BaseType::mDofSet.begin();
                    itDof != BaseType::mDofSet.end(); itDof++)
            {
                KeyType CurrVar = itDof->GetVariable().Key();
                if ((CurrVar == VELOCITY_X) || (CurrVar == VELOCITY_Y)
                        || (CurrVar == VELOCITY_Z))
                {
                    if (itDof->IsFree())
                        itDof->SetEquationId(FreeId++);
                }
                else
                {
                    if (itDof->IsFixed())
                        itDof->SetEquationId(--FixedId);
                }
            }

            mVelFreeDofs = FreeId;
            mVelFixedDofsEnd = FixedId;

            for (typename DofsArrayType::iterator itDof = BaseType::mDofSet.begin();
                    itDof != BaseType::mDofSet.end(); itDof++)
            {
                KeyType CurrVar = itDof->GetVariable().Key();
                if ((CurrVar == VELOCITY_X) || (CurrVar == VELOCITY_Y)
                        || (CurrVar == VELOCITY_Z))
                {
                    if (itDof->IsFixed())
                        itDof->SetEquationId(--FixedId);
                }
                else
                {
                    if (itDof->IsFree())
                        itDof->SetEquationId(FreeId++);
                }
            }

//            mVelDofsNum += (mVelFixedDofsEnd - FreeId);
            mPressFreeDofs = FreeId - mVelFreeDofs;
            BaseType::mEquationSystemSize = mVelFreeDofs + mPressFreeDofs;

            KRATOS_CATCH("");
        }

        /* Initialize pointers to the different vectors involved, check that they
         * have the correct size
         */
        void ResizeAndInitializeVectors(
                TSystemMatrixPointerType& pA,
                TSystemVectorPointerType& pDx,
                TSystemVectorPointerType& pb,
                ElementsArrayType& rElements,
                ConditionsArrayType& rConditions,
                ProcessInfo& CurrentProcessInfo)
        {
            KRATOS_TRY;

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

            // Member Matrices
            if (mpS == NULL)
            {
                TSystemMatrixPointerType pNewS(new TSystemMatrixType(0, 0));
                mpS.swap(pNewS);
            }
            if (mpD == NULL)
            {
                TSystemMatrixPointerType pNewD(new TSystemMatrixType(0, 0));
                mpD.swap(pNewD);
            }
            if (mpG == NULL)
            {
                TSystemMatrixPointerType pNewG(new TSystemMatrixType(0, 0));
                mpG.swap(pNewG);
            }
            if (mpL == NULL)
            {
                TSystemMatrixPointerType pNewL(new TSystemMatrixType(0, 0));
                mpL.swap(pNewL);
            }

            if (mpIDiagS == NULL) //if the pointer is not initialized initialize it to an empty matrix
            {
                TSystemVectorPointerType pNewIDiagS = TSystemVectorPointerType(new TSystemVectorType(0));
                mpIDiagS.swap(pNewIDiagS);
            }

            TSystemMatrixType& A  = *pA;
            TSystemVectorType& Dx = *pDx;
            TSystemVectorType& b = *pb;

            TSystemMatrixType& S = *mpS;
            TSystemMatrixType& G = *mpG;
            TSystemMatrixType& D = *mpD;
            TSystemMatrixType& L = *mpL;

            TSystemVectorType& IDiagS = *mpIDiagS;

            //resizing the system vectors and matrix
            if (BaseType::GetReshapeMatrixFlag() == true || // if we must remesh
                    S.size1() == 0 || D.size1() == 0 || G.size1() == 0 ||
                     L.size1() == 0 || A.size1() == 0 ) //if the matrices are not initialized
            {
                S.resize(mVelFreeDofs, mVelFreeDofs, false);
                L.resize(mPressFreeDofs, mPressFreeDofs, false);

                if( mRebuildLevel < 2 || mInitializedMatrices == false)
                {
                    G.resize(mVelFreeDofs, mPressFreeDofs, false);
                    D.resize(mPressFreeDofs, mVelFreeDofs, false);
                }

                ConstructMatrixStructure(S, D, G, L, rElements, rConditions, CurrentProcessInfo);

                if( mRebuildLevel == 0 || mInitializedMatrices == false)
                {
                    A.resize(mPressFreeDofs, mPressFreeDofs, false);
                    IDiagS.resize(mVelFreeDofs);

                    AllocateSystemMatrix(A);
                    ConstructSystemMatrix(A);
                }

                //mInitializedMatrices = true; Set flag once the matrix is BUILT (not just resized)

            } else {
                if (S.size1() != mVelFreeDofs || S.size2() != mVelFreeDofs ||
                        D.size1() != mVelFreeDofs || D.size2() != mPressFreeDofs ||
                        G.size1() != mPressFreeDofs || G.size2() != mVelFreeDofs ||
                        L.size1() != mPressFreeDofs || L.size2() != mPressFreeDofs ||
                        A.size1() != mPressFreeDofs || A.size2() != mPressFreeDofs)
                {
                    KRATOS_WATCH("it should not come here!!!!!!!! ... this is SLOW");

                    S.resize(mVelFreeDofs, mVelFreeDofs, false);
                    L.resize(mPressFreeDofs, mPressFreeDofs, false);

                    if( mRebuildLevel < 2 || mInitializedMatrices == false)
                    {
                        G.resize(mVelFreeDofs, mPressFreeDofs, false);
                        D.resize(mPressFreeDofs, mVelFreeDofs, false);
                    }

                    ConstructMatrixStructure(S, D, G, L, rElements, rConditions, CurrentProcessInfo);

                    if( mRebuildLevel == 0 || mInitializedMatrices == false)
                    {
                        A.resize(mPressFreeDofs, mPressFreeDofs, false);
                        IDiagS.resize(mVelFreeDofs);

                        AllocateSystemMatrix(A);
                        ConstructSystemMatrix(A);
                    }

                    //mInitializedMatrices = true;
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

            KRATOS_CATCH("");
        }

        void InitializeSolutionStep(
                ModelPart& rModelPart,
                TSystemMatrixType& A,
                TSystemVectorType& Dx,
                TSystemVectorType& b)
        {
            KRATOS_TRY;
            mFirstIteration = true;
            mInitializedMatrices = false;
            KRATOS_CATCH("");
        }

        void FinalizeSolutionStep(
                ModelPart& rModelPart,
                TSystemMatrixType& A,
                TSystemVectorType& Dx,
                TSystemVectorType& b)
        {}

        /* Calculate Reactions */
        void CalculateReactions(
                typename TSchemeType::Pointer pScheme,
                ModelPart& rModelPart,
                TSystemMatrixType& A,
                TSystemVectorType& Dx,
                TSystemVectorType& b)
        {
            //refresh RHS to have the correct reactions
            BuildRHS(pScheme, rModelPart, b);

            int i;
            int systemsize = BaseType::mDofSet.size() - TSparseSpace::Size(*BaseType::mpReactionsVector);

            typename DofsArrayType::ptr_iterator it2;

            //updating variables
            TSystemVectorType& ReactionsVector = *BaseType::mpReactionsVector;
            for (it2 = BaseType::mDofSet.ptr_begin(); it2 != BaseType::mDofSet.ptr_end(); ++it2)
            {
                if ((*it2)->IsFixed())
                {
                    i = (*it2)->EquationId();
                    i -= systemsize;

                    (*it2)->GetSolutionStepReactionValue() = ReactionsVector[i];
                }
            }
        }

        void ApplyDirichletConditions(
                typename TSchemeType::Pointer pScheme,
                ModelPart& rModelPart,
                TSystemMatrixType& A,
                TSystemVectorType& Dx,
                TSystemVectorType& b)
        {}



        void ApplyPointLoads(
                typename TSchemeType::Pointer pScheme,
                ModelPart& rModelPart,
                TSystemVectorType& b)
        {}

        /* this function is intended to be called at the end of the solution
         * step to clean up memory storage not needed
         */
        void Clear()
        {
            this->mDofSet = DofsArrayType();

            if (this->mpReactionsVector != NULL)
                TSparseSpace::Clear((this->mpReactionsVector));

            if (this->mpS != NULL)
                TSparseSpace::Clear((this->mpS));

            if (this->mpL != NULL)
                TSparseSpace::Clear((this->mpL));

            if (mRebuildLevel < 2)
            {
                if (this->mpG != NULL)
                    TSparseSpace::Clear((this->mpG));

                if (this->mpD != NULL)
                    TSparseSpace::Clear((this->mpD));

                if (mRebuildLevel == 0)
                {
                    if (this->mpIDiagS != NULL)
                        TSparseSpace::Clear((this->mpIDiagS));
                }
            }

            if (this->GetEchoLevel() > 0) {
                KRATOS_WATCH("ResidualBasedEliminationBuilderAndSolver Clear Function called");
            }
        }

        void SetRebuildSystemMatricesFlag(bool MustRebuild = true)
        {
            mInitializedMatrices = MustRebuild;
        }

    protected:

        virtual void ConstructMatrixStructure(
                TSystemMatrixType& S,
                TSystemMatrixType& D,
                TSystemMatrixType& G,
                TSystemMatrixType& L,
                const ElementsContainerType& rElements,
                const ConditionsArrayType& rConditions,
                ProcessInfo& CurrentProcessInfo)
        {
            Timer::Start("ConstructMatrices");
            std::vector< std::vector<std::size_t> > indicesS(mVelFreeDofs);
            std::vector< std::vector<std::size_t> > indicesG(mVelFreeDofs);
            std::vector< std::vector<std::size_t> > indicesD(mPressFreeDofs);
            std::vector< std::vector<std::size_t> > indicesL(mPressFreeDofs);

            Element::EquationIdVectorType ids;
            ids.reserve(16); // 16 as initial capacity: 4 Dofs per node assumed

            // Identify and collect the indices of non-zero terms in each matrix
            for (typename ElementsContainerType::const_iterator itElem = rElements.begin();
                    itElem != rElements.end(); itElem++)
            {
                itElem->EquationIdVector(ids, CurrentProcessInfo);

                for (std::size_t i = 0; i < ids.size(); i++)
                    if (ids[i] < mVelFreeDofs)
                    {
                        std::vector<std::size_t>& RowS = indicesS[ids[i]];
                        std::vector<std::size_t>& RowG = indicesG[ids[i]];

                        for (std::size_t j = 0; j < ids.size(); j++)
                        {
                            if (ids[j] < mVelFreeDofs)
                                AddUnique(RowS, ids[j]);
                            else if (ids[j] < BaseType::mEquationSystemSize)
                                AddUnique(RowG, ids[j] - mVelFreeDofs);
                        }
                    }
                    else if (ids[i] < BaseType::mEquationSystemSize)
                    {
                        std::vector<std::size_t>& RowD = indicesD[ids[i] - mVelFreeDofs];
                        std::vector<std::size_t>& RowL = indicesL[ids[i] - mVelFreeDofs];

                        for (std::size_t j = 0; j < ids.size(); j++)
                        {
                            if (ids[j] < mVelFreeDofs)
                                AddUnique(RowD, ids[j]);
                            else if (ids[j] < BaseType::mEquationSystemSize)
                                AddUnique(RowL, ids[j] - mVelFreeDofs);
                        }
                    }
            }

            // Do the same for conditions
            for (typename ConditionsArrayType::const_iterator itCond = rConditions.begin();
                    itCond != rConditions.end(); itCond++)
            {
                itCond->EquationIdVector(ids, CurrentProcessInfo);

                for (std::size_t i = 0; i < ids.size(); i++)
                    if (ids[i] < mVelFreeDofs)
                    {
                        std::vector<std::size_t>& RowS = indicesS[ids[i]];
                        std::vector<std::size_t>& RowG = indicesG[ids[i]];

                        for (std::size_t j = 0; j < ids.size(); j++)
                        {
                            if (ids[j] < mVelFreeDofs)
                                AddUnique(RowS, ids[j]);
                            else if (ids[j] < BaseType::mEquationSystemSize)
                                AddUnique(RowG, ids[j] - mVelFreeDofs);
                        }
                    }
                    else if (ids[i] < BaseType::mEquationSystemSize)
                    {
                        std::vector<std::size_t>& RowD = indicesD[ids[i] - mVelFreeDofs];
                        std::vector<std::size_t>& RowL = indicesL[ids[i] - mVelFreeDofs];

                        for (std::size_t j = 0; j < ids.size(); j++)
                        {
                            if (ids[j] < mVelFreeDofs)
                                AddUnique(RowD, ids[j]);
                            else if (ids[j] < BaseType::mEquationSystemSize)
                                AddUnique(RowL, ids[j] - mVelFreeDofs);
                        }
                    }
            }

            // Allocate memory and initialize matrices with zeros
            int NumTermsS = 0; // Counters for non-zero terms
            int NumTermsL = 0;

            for (std::size_t i = 0; i < indicesS.size(); i++)
                NumTermsS += indicesS[i].size();

            S.reserve(NumTermsS, false);

            for (std::size_t i = 0; i < indicesL.size(); i++)
                NumTermsL += indicesL[i].size();

            L.reserve(NumTermsL, false);

            if (mRebuildLevel < 2 || mInitializedMatrices == false)
            {
                int NumTermsG = 0;
                int NumTermsD = 0;

                for (std::size_t i = 0; i < indicesG.size(); i++)
                    NumTermsG += indicesG[i].size();

                G.reserve(NumTermsG, false);

                for (std::size_t i = 0; i < indicesD.size(); i++)
                    NumTermsD += indicesD[i].size();

                D.reserve(NumTermsD, false);
            }
            Timer::Stop("ConstructMatrices");

            // Fill with zeros the matrix (create the structure)
            Timer::Start("MatrixStructure");

            // Note that AllocateSpace has different definitions depending on _OPENMP
            AllocateSpace(S, indicesS);
            AllocateSpace(L, indicesL);

            if (mRebuildLevel < 2 || mInitializedMatrices == false)
            {
                AllocateSpace(G, indicesG);
                AllocateSpace(D, indicesD);
            }

            Timer::Stop("MatrixStructure");
        }

    private:

        /* Add a matrix position to the list (check that it hasn't been used before) */
        void AddUnique(
                std::vector<std::size_t>& Row,
                const std::size_t& Candidate)
        {
            std::vector<std::size_t>::iterator i = Row.begin();
            std::vector<std::size_t>::iterator endit = Row.end();

            while ( i != endit && (*i) != Candidate)
                i++;
            if( i == endit )
                Row.push_back(Candidate);
        }



        /*
         * Assembly functions
         */

        // ALL ASSEMBLY FUNCTIONS NEED OPENMP. CONSIDER ADDING Assemble_OnFreeRows
        void Assemble(
                TSystemMatrixType& A,
                TSystemVectorType& b,
                LocalSystemMatrixType& LHS_Contribution,
                LocalSystemVectorType& RHS_Contribution,
                Element::EquationIdVectorType& EquationId)
        {
            unsigned int ContributionSize = EquationId.size();
            for (unsigned int i = 0; i < ContributionSize; i++)
            {
                unsigned int Global_i = EquationId[i];

                if (Global_i < mVelFreeDofs)
                {
                    for (unsigned int j = 0; j < ContributionSize; j++)
                    {
                        unsigned int Global_j = EquationId[j];
                        if (Global_j < mVelFreeDofs)
                        {
                            mpS->operator()(Global_i, Global_j) += LHS_Contribution(i, j);
                        }
                        else if (Global_j < BaseType::mEquationSystemSize)
                        {
                            mpG->operator()(Global_i, Global_j - mVelFreeDofs) += LHS_Contribution(i, j);
                        }
                    }
                    b[Global_i] += RHS_Contribution[i];
                }
                else if (Global_i < BaseType::mEquationSystemSize)
                {
                    for (unsigned int j = 0; j < ContributionSize; j++)
                    {
                        unsigned int Global_j = EquationId[j];
                        if (Global_j < mVelFreeDofs)
                        {
                            mpD->operator()(Global_i - mVelFreeDofs, Global_j) += LHS_Contribution(i, j);
                        }
                        else if (Global_j < BaseType::mEquationSystemSize)
                        {
                            mpL->operator()(Global_i - mVelFreeDofs, Global_j - mVelFreeDofs) += LHS_Contribution(i, j);
                        }
                    }
                    b[Global_i] += RHS_Contribution[i];
                }
            }
        }

#ifdef _OPENMP

        void Parallel_Assemble(
                TSystemMatrixType& A,
                TSystemVectorType& b,
                const LocalSystemMatrixType& LHS_Contribution,
                const LocalSystemVectorType& RHS_Contribution,
                Element::EquationIdVectorType& EquationId,
                std::vector< omp_lock_t >& lock_array)
        {

            unsigned int ContributionSize = EquationId.size();

            for (unsigned int i = 0; i < ContributionSize; i++)
            {
                unsigned int Global_i = EquationId[i];

                if ( Global_i < mVelFreeDofs )
                {
                    omp_set_lock(&lock_array[Global_i]);

                    for (unsigned int j = 0; j < ContributionSize; j++)
                    {
                        unsigned int Global_j = EquationId[j];
                        if (Global_j < mVelFreeDofs)
                        {
                            mpS->operator()(Global_i, Global_j) += LHS_Contribution(i, j);
                        }
                        else if (Global_j < BaseType::mEquationSystemSize)
                        {
                            mpG->operator()(Global_i, Global_j - mVelFreeDofs) += LHS_Contribution(i, j);
                        }
                    }

                    b[Global_i] += RHS_Contribution[i];

                    omp_unset_lock(&lock_array[Global_i]);

                } else if (Global_i < BaseType::mEquationSystemSize) {

                    omp_set_lock(&lock_array[Global_i]);

                    for (unsigned int j = 0; j < ContributionSize; j++)
                    {
                        unsigned int Global_j = EquationId[j];
                        if (Global_j < mVelFreeDofs)
                        {
                            mpD->operator()(Global_i - mVelFreeDofs, Global_j) += LHS_Contribution(i, j);
                        }
                        else if (Global_j < BaseType::mEquationSystemSize)
                        {
                            mpL->operator()(Global_i - mVelFreeDofs, Global_j - mVelFreeDofs) += LHS_Contribution(i, j);
                        }
                    }

                    b[Global_i] += RHS_Contribution[i];

                    omp_unset_lock(&lock_array[Global_i]);
                }
            }
        }

#endif

        void AssembleLHS(
                TSystemMatrixType& A,
                LocalSystemMatrixType& LHS_Contribution,
                Element::EquationIdVectorType& EquationId)
        {
            TSystemMatrixType& rS = *mpS;
            TSystemMatrixType& rD = *mpD;
            TSystemMatrixType& rG = *mpG;
            TSystemMatrixType& rL = *mpL;

            unsigned int ContributionSize =  EquationId.size();

            for( unsigned int i = 0; i < ContributionSize; i++)
            {
                unsigned int Global_i = EquationId[i];
                if ( Global_i < mVelFreeDofs )
                {
                    for( unsigned int j = 0; j < ContributionSize; j++)
                    {
                        unsigned int Global_j = EquationId[j];
                        if( Global_j < mVelFreeDofs)
                        {
                            rS(Global_i,Global_j) += LHS_Contribution(i,j);
                        }
                        else if( Global_j < BaseType::mEquationSystemSize)
                        {
                            rG(Global_i,Global_j - mVelFreeDofs) += LHS_Contribution(i,j);
                        }
                    }
                }
                else if( Global_i < BaseType::mEquationSystemSize )
                {
                    for( unsigned int j = 0; j < ContributionSize; j++)
                    {
                        unsigned int Global_j = EquationId[j];
                        if( Global_j < mVelFreeDofs)
                        {
                            rD(Global_i - mVelFreeDofs, Global_j) += LHS_Contribution(i,j);
                        }
                        else if( Global_j < BaseType::mEquationSystemSize)
                        {
                            rL(Global_i - mVelFreeDofs, Global_j - mVelFreeDofs) += LHS_Contribution(i,j);
                        }
                    }
                }
            }
        }

        void AssembleLHS_CompleteOnFreeRows(
                TSystemMatrixType& A,
                LocalSystemMatrixType& LHS_Contribution,
                Element::EquationIdVectorType& EquationId)
        {
            unsigned int ContributionSize =  EquationId.size();

            for( unsigned int i = 0; i < ContributionSize; i++)
            {
                unsigned int Global_i = EquationId[i];
                if ( Global_i < mVelFreeDofs )
                {
                    for( unsigned int j = 0; j < ContributionSize; j++)
                    {
                        unsigned int Global_j = EquationId[j];
                        if( Global_j < mVelFreeDofs)
                        {
                            mpS->operator()(Global_i,Global_j) += LHS_Contribution(i,j);
                        }
                        else
                        {
                            mpG->operator()(Global_i,Global_j - mVelFreeDofs) += LHS_Contribution(i,j);
                        }
                    }
                }
                else if( Global_i < BaseType::mEquationSystemSize )
                {
                    for( unsigned int j = 0; j < ContributionSize; j++)
                    {
                        unsigned int Global_j = EquationId[j];
                        if( Global_j < mVelFreeDofs)
                        {
                            mpD->operator()(Global_i - mVelFreeDofs, Global_j) += LHS_Contribution(i,j);
                        }
                        else
                        {
                            mpL->operator()(Global_i - mVelFreeDofs, Global_j - mVelFreeDofs) += LHS_Contribution(i,j);
                        }
                    }
                }
            }
        }


        void AssembleRHS(
                TSystemVectorType& b,
                LocalSystemVectorType& RHS_Contribution,
                Element::EquationIdVectorType& EquationId)
        {
            unsigned int ContributionSize =  EquationId.size();

            for( unsigned int i = 0; i < ContributionSize; i++)
            {
                unsigned int Global_i = EquationId[i];
                if ( Global_i < BaseType::mEquationSystemSize )
                {
                    b[Global_i] += RHS_Contribution[i];
                }
            }
        }


        /*
         * Pressure System matrix functions
         * Compute L - D*Inv(Diag(S))*G
         */
#ifndef _OPENMP

        void AllocateSystemMatrix(TSystemMatrixType& A)
        {
            /* All non-zero positions in A = L - D*Inv(Diag(S))*G have to be stored.
             * This method allocates the required memory based on the shapes of
             * member matrices mpD (Divergence operator), mpG (Gradient Operator)
             * and mpL (stabilization term)
             */

            TSystemMatrixType& rG = *mpG;
            TSystemMatrixType& rD = *mpD;
            TSystemMatrixType& rL = *mpL;

            std::size_t NumTerms = 0;
            std::vector<int> mask(rL.size2(),-1);
            // Keeps track of used cols in a given row.
            // When a col is used, mask[col] is filled with row num.

            for( OuterIt RowD = rD.begin1(), RowL = rL.begin1() ;
                    RowD != rD.end1();
                    RowD++, RowL++)
            {
                // Find terms filled by the matrix product
                for( InnerIt ItD =  RowD.begin(); ItD != RowD.end() ; ItD++ )
                {
                    RowType RowG(rG,ItD.index2());
                    for( typename RowType::iterator ItG = RowG.begin(); ItG != RowG.end(); ItG++ )
                    {
                        if( mask[ItG.index()] != int (ItD.index1()) )
                            // Cast to int to avoid a compilation warning, as index1() returns an unsigned int
                        {
                            mask[ItG.index()] = ItD.index1();
                            NumTerms++;
                        }
                    }
                }

                // Find extra terms introduced by matrix difference
                for( InnerIt ItL = RowL.begin(); ItL != RowL.end(); ItL++ )
                {
                    if( mask[ItL.index2()] != int (ItL.index1()) )
                        // Cast to int to avoid a compilation warning, as index1() returns an unsigned int
                    {
                        mask[ItL.index2()] = ItL.index1();
                        NumTerms++;
                    }
                }
            }

            A.reserve(NumTerms);
        }

#else
        // we can't allocate in parallel!!
        void AllocateSystemMatrix(TSystemMatrixType& A)
        {}

#endif


        /* Identify non-zero tems in the system matrix */
        void ConstructSystemMatrix(
                TSystemMatrixType& A)
        {
            // Retrieve matrices
            TSystemMatrixType& rG = *mpG;
            TSystemMatrixType& rD = *mpD;
            TSystemMatrixType& rL = *mpL;

            PartitionVector Partition;
            unsigned int NumThreads = GetNumThreads();

            DivideInPartitions(A.size1(),Partition,NumThreads);

            for( unsigned int k = 0 ; k < NumThreads ; k++)
            {
                #pragma omp parallel
                if( ThisThread() == k) {
                    boost::shared_ptr< IndexVector > pNext( new IndexVector(mPressFreeDofs) );
                    IndexVector& Next = *pNext; // Keeps track of which columns were filled
                    for(unsigned int m = 0; m < mPressFreeDofs; m++) Next[m] = -1;

                    std::size_t NumTerms = 0; // Full positions in a row
                    boost::shared_ptr< std::vector<unsigned int> > pUsedCols( new std::vector<unsigned int>);
                    std::vector<unsigned int>& UsedCols = *pUsedCols;
                    UsedCols.reserve(mPressFreeDofs);

                    for( std::size_t RowIndex = Partition[k] ;
                            RowIndex != Partition[k+1] ; RowIndex++ ) {
                        RowType RowD(rD,RowIndex);
                        RowType RowL(rL,RowIndex);

                        int head = -2;
                        std::size_t Length = 0;

                        // Terms filled by L
                        for( typename RowType::iterator ItL = RowL.begin();
                                ItL != RowL.end(); ItL++ ) {
                            if( Next[ItL.index()] == -1) {
                                Next[ItL.index()] = head;
                                head = ItL.index();
                                Length++;
                            }
                        }

                        // Additional terms due to D*Inv(Diag(S))*G
                        for( typename RowType::iterator ItD = RowD.begin();
                                ItD != RowD.end(); ItD++ ) {
                            RowType RowG(rG,ItD.index());

                            for( typename RowType::iterator ItG = RowG.begin();
                                    ItG != RowG.end(); ItG++ ) {
                                if( Next[ItG.index()] == -1) {
                                    Next[ItG.index()] = head;
                                    head = ItG.index();
                                    Length++;
                                }
                            }
                        }

                        // Identify full terms for ordering
                        for( std::size_t i = 0; i < Length; i++) {
                            if( Next[head] != -1 ) {
                                UsedCols.push_back(head);
                                NumTerms++;
                            }

                            int temp = head;
                            head = Next[head];

                            // Clear 'Next' for next iteration
                            Next[temp] = -1;
                        }

                        // Sort Column indices
                        SortCols(UsedCols,NumTerms);

                        // Store row in matrix, clean temporary variables
                        for( unsigned int i = 0; i < NumTerms; i++) {
                            A.push_back(RowIndex,UsedCols[i],0);
                        }
                        NumTerms = 0;
                        UsedCols.resize(0,false);
                    }
                }
            }
        }

        /* Compute the System Matrix A = L - D*Inv(Diag(S))*G. The multiplication
         * is performed in random order, so each row will be stored in a temporary
         * variable, ordered and copied in input matrix A.
         */
        void CalculateSystemMatrix(
                TSystemMatrixType& A)
        {
            Timer::Start("CalculateA");
            // Retrieve matrices
            TSystemMatrixType& rS = *mpS;
            TSystemMatrixType& rG = *mpG;
            TSystemMatrixType& rD = *mpD;
            TSystemMatrixType& rL = *mpL;

            // Compute Inv(Diag(S))
            TSystemVectorType& rIDiagS = *mpIDiagS;
            int DiagSize = int(mVelFreeDofs); // to avoid comparison between int & unsigned int
            #pragma omp parallel for
            for( int i = 0; i < DiagSize; i++)
                rIDiagS[i] = 1/rS(i,i);

            PartitionVector Partition;
            unsigned int NumThreads = GetNumThreads();

            DivideInPartitions(A.size1(),Partition,NumThreads);

            #pragma omp parallel
            {
                unsigned int k = ThisThread();
                TSystemVectorPointerType pCurrentRow( new TSystemVectorType(mPressFreeDofs) );
                TSystemVectorType& CurrentRow = *pCurrentRow; // Values for current row
                for (unsigned int i = 0; i < mPressFreeDofs; i++) CurrentRow[i] = 0.0;

                boost::shared_ptr< IndexVector > pNext( new IndexVector(mPressFreeDofs) );
                IndexVector& Next = *pNext; // Keeps track of which columns were filled
                for(unsigned int m=0; m < mPressFreeDofs; m++) Next[m] = -1;

                for( std::size_t RowIndex = Partition[k] ;
                        RowIndex != Partition[k+1] ; RowIndex++ ) {
                    RowType RowD(rD,RowIndex);
                    RowType RowL(rL,RowIndex);

                    int head = -2;
                    std::size_t Length = 0;

                    // Write L in A
                    for( typename RowType::iterator ItL = RowL.begin();
                            ItL != RowL.end(); ItL++ ) {
                        CurrentRow(ItL.index()) = *ItL;

                        if( Next[ItL.index()] == -1)
                        {
                            Next[ItL.index()] = head;
                            head = ItL.index();
                            Length++;
                        }
                    }

                    // Substract D*Inv(Diag(S))*G
                    for( typename RowType::iterator ItD = RowD.begin();
                            ItD != RowD.end(); ItD++ ) {
                        RowType RowG(rG,ItD.index());

                        for( typename RowType::iterator ItG = RowG.begin();
                                ItG != RowG.end(); ItG++ ) {
                            CurrentRow[ItG.index()] -= (*ItD) * rIDiagS[ItD.index()] * (*ItG);

                            if( Next[ItG.index()] == -1)
                            {
                                Next[ItG.index()] = head;
                                head = ItG.index();
                                Length++;
                            }
                        }
                    }

                    // Fill matrix row (Random access, as ConstructSystemMatrix
                    // guarantees that required matrix terms already exist),
                    // then clean temporary variables.
                    for( std::size_t i = 0; i < Length; i++) {
                        if( CurrentRow[head] != 0 ) {
                            A(RowIndex,head) = CurrentRow[head];
                            CurrentRow[head] = 0;
                        }

                        // Clear variables for next iteration
                        int temp = head;
                        head = Next[head];
                        Next[temp] = -1;
                    }
                }
            }
            Timer::Stop("CalculateA");
        }


        /* Helper function for Sytem matrix functions */
        void SortCols(
                std::vector<unsigned int>& ColList,
                std::size_t& NumCols)
        {
            bool swap = true;
            unsigned int d = NumCols;
            int temp;

            while( swap || d > 1 )
            {
                swap = false;
                d = (d+1)/2;
                for( unsigned int i=0; i<(NumCols - d); i++)
                    if( ColList[i+d] < ColList[i] )
                    {
                        temp = ColList[i+d];
                        ColList[i+d] = ColList[i];
                        ColList[i] = temp;
                        swap = true;
                    }
            }
        }


        /* allocate space for member matrices */

        void AllocateSpace(
                TSystemMatrixType& A,
                std::vector< std::vector<std::size_t> >& indices)
        {
            Timer::Start("AllocateOther");
            unsigned int NumThreads = GetNumThreads();
            PartitionVector MatrixPartition;

            DivideInPartitions(indices.size(),MatrixPartition,NumThreads);

            for( unsigned int k=0; k < NumThreads; k++ )
            {
                #pragma omp parallel
                if( ThisThread() == k )
                {
                    for( std::size_t i = MatrixPartition[k]; i < MatrixPartition[k+1]; i++ )
                    {
                        std::vector<std::size_t>& Row = indices[i];
                        std::sort(Row.begin(), Row.end());

                        for(std::vector<std::size_t>::iterator it= Row.begin(); it != Row.end() ; it++)
                        {
                            A.push_back(i,*it,0.00);
                        }
                        Row.clear();
                    }
                }
            }
            Timer::Start("AllocateOther");
        }

#ifdef _OPENMP

        /* y += A*x in parallel */
        void Parallel_ProductAdd(
                const TSystemMatrixType& A,
                const TSystemVectorType& in,
                TSystemVectorType& out)
        {
            Timer::Start("CustomProduct");
            typedef  unsigned int size_type;
            typedef  double value_type;

            //create partition
            PartitionVector partition;
            unsigned int number_of_threads = GetNumThreads();
            DivideInPartitions(A.size1(),partition,number_of_threads);

            #pragma omp parallel
            {
                int thread_id = omp_get_thread_num();
                int number_of_rows = partition[thread_id+1] - partition[thread_id];
                typename compressed_matrix<TDataType>::index_array_type::const_iterator row_iter_begin = A.index1_data().begin()+partition[thread_id];
                 typename compressed_matrix<TDataType>::index_array_type::const_iterator index_2_begin = A.index2_data().begin()+*row_iter_begin;
                 typename compressed_matrix<TDataType>::value_array_type::const_iterator value_begin = A.value_data().begin()+*row_iter_begin;


                partial_product_add(number_of_rows,
                                    row_iter_begin,
                                    index_2_begin,
                                    value_begin,
                                    in,
                                    partition[thread_id],
				    out);
            }
            Timer::Stop("CustomProduct");
        }

	/**
	 * calculates partial product
	 */
        void partial_product_add(
                int number_of_rows,
                typename compressed_matrix<TDataType>::index_array_type::const_iterator row_begin,
                typename compressed_matrix<TDataType>::index_array_type::const_iterator index2_begin,
                typename compressed_matrix<TDataType>::value_array_type::const_iterator value_begin,
                const TSystemVectorType& input_vec,
		unsigned int output_begin_index,
		TSystemVectorType& output_vec)
        {
            int row_size;
	    int kkk = output_begin_index;
            typename TSystemMatrixType::index_array_type::const_iterator row_it = row_begin;
            for(int k = 0; k < number_of_rows; k++)
            {
                row_size= *(row_it+1)-*row_it;
                row_it++;
                TDataType t = TDataType();

                for(int i = 0; i<row_size; i++)
                    t += *value_begin++ *   ( input_vec[*index2_begin++]);

		output_vec[kkk++] += t;

            }
        }

#endif

        /* Set iterative solver tolerance using inexact Newton criteria */
        void SetInitialTolerance(
                double RHSNorm,
                double& TolFactor)
        {
            TolFactor = mMaxTolFactor;
//            double Tolerance = TolFactor*RHSNorm;
//            (BaseType::mpLinearSystemSolver)->SetTolerance(Tolerance);
            (BaseType::mpLinearSystemSolver)->SetTolerance(mSmallTol);
//            (BaseType::mpLinearSystemSolver)->SetTolerance(TolFactor);
            std::cout << "Set iterative solver tolerance to " << TolFactor << std::endl;
        }

        void UpdateTolerance(
                double OldRHSNorm,
                double NewRHSNorm,
                double& TolFactor)
        {
            const double MaxDecreaseFactor = 0.1;

            double CandidateFactor = mGamma*(NewRHSNorm*NewRHSNorm)/(OldRHSNorm*OldRHSNorm);
            std::cout << "Norm Ratio: " << NewRHSNorm/OldRHSNorm << std::endl;
            double CandidateFactor_LimitedDecrease = mGamma*TolFactor*TolFactor;

            if (CandidateFactor_LimitedDecrease < MaxDecreaseFactor) {
                TolFactor = (CandidateFactor < mMaxTolFactor) ? CandidateFactor : mMaxTolFactor;
            } else {
                double Temp = (CandidateFactor > CandidateFactor_LimitedDecrease) ? CandidateFactor : CandidateFactor_LimitedDecrease;
                TolFactor = (Temp < mMaxTolFactor) ? Temp : mMaxTolFactor;
            }

//            double Tolerance = TolFactor * NewRHSNorm;
//            if (Tolerance < mSmallTol)
//            {
//                Tolerance = mSmallTol;
//                TolFactor = mSmallTol / NewRHSNorm;
//            }
//            (BaseType::mpLinearSystemSolver)->SetTolerance(Tolerance);
            if (TolFactor < mSmallTol) TolFactor = mSmallTol;
            (BaseType::mpLinearSystemSolver)->SetTolerance(TolFactor);
            std::cout << "Corrected iterative solver tolerance to " << TolFactor << std::endl;
        }

        /* The following functions retrieve basic OpenMP information or a
         * suitable alternative when they are compilied without OpenMP
         */

        inline unsigned int GetNumThreads()
        {
            #ifdef _OPENMP
            return omp_get_max_threads();
            #else
            return 1;
            #endif
        }

        inline unsigned int ThisThread()
        {
            #ifdef _OPENMP
            return omp_get_thread_num();
            #else
            return 0;
            #endif
        }

        /* Divide the matrix in groups of contiguous rows.
         * Each group will be assigned to a different OMP thread.
         */

        inline void DivideInPartitions(
                const int NumTerms,
                PartitionVector& Partitions,
                unsigned int& NumThreads)
        {
            #ifdef _OPENMP
            Partitions.resize(NumThreads + 1);
            int PartitionSize = NumTerms / NumThreads;
            Partitions[0] = 0;
            Partitions[NumThreads] = NumTerms;
            for(int i = 1; i < NumThreads; i++)
                Partitions[i] = Partitions[i-1] + PartitionSize ;
            #else
            Partitions.resize(2);
            Partitions[0] = 0;
            Partitions[1] = NumTerms;
            #endif
        }


        // Total number of Dofs: BaseType::mEquationSystemSize;
//        unsigned int mVelDofsNum; // Number of Velocity Dofs
        unsigned int mVelFreeDofs; // Position of 1st Free Pressure Dof
        unsigned int mVelFixedDofsEnd; // Position of 1st Fixed Pressure Dof

        unsigned int mPressFreeDofs;

        TSystemMatrixPointerType mpS; // System Matrix for the momentum equation
        TSystemMatrixPointerType mpD; // Discrete Divergence operator
        TSystemMatrixPointerType mpG; // Discrete Gradient operator
        TSystemMatrixPointerType mpL; // Stabilization term

        TSystemVectorPointerType mpIDiagS; // Inv(Diag(S)), stored as a vector

        // Flags for matrix reconstruction
        unsigned int mRebuildLevel;
        /* Decide what to rebuil each time new matrices are written
         * 0: Rebuild all matrices
         * 1: Keep A (pressure system matrix) constant
         * 2: Keep D, G and A constant
         */
        bool mInitializedMatrices;

        unsigned int mVelocityCorrection;
        /* Ensure that velocity is divergence free after each iteration
         * 0: Do not
         * 1: Solve the system Diag(S)*dv = -G*dp
         * 2: Solve the full system S*dv = -G*dp
         */

        // Variables for inexact Newton tolerance control
        bool mInexactNewton; // Use Inexact Newton method
        double mMaxTolFactor; // used to compute maximum solution tolerance
        double mSmallTol; // Minimum solution tolerance: 0.5*NonLinearTol
        double mGamma;

        bool mFirstIteration; // Keeps track of the first iteration in each step

//        double mVelTolFactor;
        double mPressTolFactor;

//        double mLastVelRHSNorm;
        double mLastPressRHSNorm;
    };
}

#endif	/* KRATOS_PRESSURE_SPLITTING_BUILDER_AND_SOLVER_H */
