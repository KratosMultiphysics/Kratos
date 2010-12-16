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
//#include <set>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "utilities/openmp_utils.h"

/* External includes */
#include <boost/smart_ptr.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

/* Project includes */
#include "includes/define.h"
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"

namespace Kratos
{

    /// A builder and solver scheme based on separating pressure and velocity Dofs
    /**
     The PressureSplittingBuilderAndSolver class is intended to solve the
     velocity and pressure degrees of freedom separately, thanks to a Pressure
     Schur complement type decomposition. Given a system equation of the type
     \f[
      \left[ \begin{array}{cc}
       S & G \\
       D & L
      \end{array} \right]
      \left[ \begin{array}{c}
       \delta u \\
       \delta p
      \end{array} \right]
      =
      \left[ \begin{array}{c}
       r_u \\
       r_p
      \end{array} \right]
     \f]
     This class divides the system and solves it as
     \f[
      S \delta u = r_u - G \, \delta p
     \f]
     \f[
      \left( L - D S^{-1} G \right) \, \delta p = r_p - D {S}^{-1} \, r_u
     \f]
     where \f$ S^{-1} \f$ is approximated by \f$ \left( Diag \left( S \right) \right)^{-1} \f$
     */
    template< class TSparseSpace,
	class TDenseSpace , //= DenseSpace<double>,
	class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
	>
    class PressureSplittingBuilderAndSolver:
        public BuilderAndSolver< TSparseSpace,TDenseSpace,TLinearSolver >
    {
        /* Please note that this builder and solver splits the system matrix
         * in 4 blocks and stores them separately. The system matrix A given as
         * input (to preserve the same function signatures as other builder and
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

        typedef OpenMPUtils::PartitionVector PartitionVector;

        // UBLAS matrix access types
        typedef typename TSystemMatrixType::iterator1 OuterIt;
        typedef typename TSystemMatrixType::iterator2 InnerIt;
        typedef typename boost::numeric::ublas::matrix_row< TSystemMatrixType > RowType;

        // Life Cycle

        /// Constructor
        PressureSplittingBuilderAndSolver(
                typename TLinearSolver::Pointer pNewVelLinearSystemSolver, // Velocity solver, stored internally
                typename TLinearSolver::Pointer pNewPressLinearSystemSolver, // Pressure solver, stored by the base class
                unsigned int VelocityCorrection, // If > 0, explicitly solve the velocity to be divergence-free at each step
                bool UseInexactNewton, // If true, dynamically set the linear iterative solver tolerance for the pressure system
                double NonLinearTol = 1e-3, // Only used if InexactNewton == true, otherwise the solver will use it's own tolerance
                double MaxTolFactor = 0.1,
                double Gamma = 0.9):
            BuilderAndSolver< TSparseSpace,TDenseSpace,TLinearSolver >(pNewPressLinearSystemSolver),
            mpVelLinearSystemSolver(pNewVelLinearSystemSolver),
            mDofSetChanged(true),
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

        /// Destructor
        virtual ~PressureSplittingBuilderAndSolver() {}

        // Operations

        /// Build System
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

#ifndef _OPENMP
            // Define contributions to the system
            LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0,0);
            LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

            // Will store the position of each Dof in the system
            Element::EquationIdVectorType EquationIds;

            ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();

            // Assemble contributions from elements
            for( typename ElementsArrayType::ptr_iterator pElem = rElements.ptr_begin();
                    pElem != rElements.ptr_end(); pElem++)
            {
                // Get Elemental Contributions
                pScheme->CalculateSystemContributions(*pElem,LHS_Contribution,
                        RHS_Contribution,EquationIds,CurrentProcessInfo);

                //assemble the elemental contribution
                Assemble(b,LHS_Contribution,RHS_Contribution,EquationIds);

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
                Assemble(b,LHS_Contribution,RHS_Contribution,EquationIds);
            }

#else
            //creating an array of lock variables of the size of the system matrix
            int TotalSize = BaseType::mEquationSystemSize;
            std::vector< omp_lock_t > lock_array(TotalSize);

            for(int i = 0; i<TotalSize; i++)
                omp_init_lock( &lock_array[i] );

            //create a partition of the element array
            int NumThreads = OpenMPUtils::GetNumThreads();
            PartitionVector ElementPartition;

            OpenMPUtils::DivideInPartitions(rElements.size(),NumThreads,ElementPartition);

            double start_prod = omp_get_wtime();

            #pragma omp parallel
            {
                int k = OpenMPUtils::ThisThread();

                //contributions to the system
                LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0,0);
                LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

                //vector containing the localization in the system of the different terms
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
                    Assemble(b,LHS_Contribution,RHS_Contribution,EquationId,lock_array);

                    // clean local elemental memory
                    pScheme->CleanMemory(*pElem);
                }
            }

            PartitionVector ConditionPartition;

            OpenMPUtils::DivideInPartitions(rConditions.size(),NumThreads,ConditionPartition);

            #pragma omp parallel
            {
                int k = OpenMPUtils::ThisThread();

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
                    Assemble(b,LHS_Contribution,RHS_Contribution,EquationId,lock_array);
                }
            }


            double stop_prod = omp_get_wtime();
            std::cout << "time: " << stop_prod - start_prod << std::endl;

            for(int i = 0; i< TotalSize; i++)
                omp_destroy_lock(&lock_array[i]);

#endif

            /* Build the pressure system matrix */
            CalculateSystemMatrix(A);

            KRATOS_CATCH("");
        }

        /// Build Left Hand Side only
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
                AssembleLHS(LHS_Contribution, EquationIds);

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
                AssembleLHS(LHS_Contribution, EquationIds);
            }
//#else
//            // DO MP VERSION HERE
//#endif

            /* Build the pressure system matrix */
            CalculateSystemMatrix(A);

            KRATOS_CATCH("");
        }

        /// Build Left Hand Side with extra columns, including fixed Dofs
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
                AssembleLHS_CompleteOnFreeRows(LHS_Contribution, EquationIds);

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
                AssembleLHS_CompleteOnFreeRows(LHS_Contribution, EquationIds);
            }
//#else
//            // DO MP VERSION HERE
//#endif

            /* Build the pressure system matrix */
            CalculateSystemMatrix(A);

            KRATOS_CATCH("");
        }

        /// Solve one iteration
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

                mpVelLinearSystemSolver->Solve(rS, rDVel, rVelRHS);
// std::cout << "velocity solution" << *(BaseType::mpLinearSystemSolver) << std::endl;

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
/*for(unsigned int i=0; i<A.size1(); i++)
	std::cout << i << " " << A(i,i) << std::endl;*/
                BaseType::mpLinearSystemSolver->Solve(A, rDPress, rPressRHS);
// std::cout << "pressure solution" << *(BaseType::mpLinearSystemSolver) << std::endl;

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

                        mpVelLinearSystemSolver->Solve(rS, rVelUpdate, rVelRHS);
// std::cout << "second velocity solution" << *(BaseType::mpLinearSystemSolver) << std::endl;
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

        /// Build and Solve system
        void BuildAndSolve(
                typename TSchemeType::Pointer pScheme,
                ModelPart& rModelPart,
                TSystemMatrixType& A,
                TSystemVectorType& Dx,
                TSystemVectorType& b)
        {
            KRATOS_TRY

            std::cout << "Building system" << std::endl;

            Build(pScheme, rModelPart, A, b);

//        //does nothing...dirichlet conditions are naturally dealt with in defining the residual
//        ApplyDirichletConditions(pScheme,rModelPart,A,Dx,b);

            if (this->GetEchoLevel() == 3) {
                std::cout << "before the solution of the system" << std::endl;
                std::cout << "System Matrix = " << A << std::endl;
                std::cout << "unknowns vector = " << Dx << std::endl;
                std::cout << "RHS vector = " << b << std::endl;
            }

            std::cout << "Solving System" << std::endl;

            SystemSolve(A, Dx, b);


            if (this->GetEchoLevel() == 3) {
                std::cout << "after the solution of the system" << std::endl;
                std::cout << "System Matrix = " << A << std::endl;
                std::cout << "unknowns vector = " << Dx << std::endl;
                std::cout << "RHS vector = " << b << std::endl;
            }

            KRATOS_CATCH("");
        }

        /// Solve System for updated Right Hand Side
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

        /// Build Right Hand Side only
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

            //create a partition of the element array
            int NumThreads = OpenMPUtils::GetNumThreads();
            PartitionVector ElementPartition;

            OpenMPUtils::DivideInPartitions(rElements.size(),NumThreads,ElementPartition);

            // Note that to assemble the LHS we use locks to avoid data races in OpenMP
            // Here we use #pragma omp atomic, as we need to write a double as oposed to
            // a row in a matrix
            #pragma omp parallel
            {
                int k = OpenMPUtils::ThisThread();

                //contributions to the system
                LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

                //vector containing the localization in the system of the different terms
                Element::EquationIdVectorType EquationId;
                ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();
                typename ElementsArrayType::ptr_iterator ElemBegin = rElements.ptr_begin()+ElementPartition[k];
                typename ElementsArrayType::ptr_iterator ElemEnd = rElements.ptr_begin()+ElementPartition[k+1];

                // assemble all elements
                for (typename ElementsArrayType::ptr_iterator pElem = ElemBegin; pElem != ElemEnd; pElem++)
                {
                    //calculate elemental contribution
                    pScheme->Calculate_RHS_Contribution(*pElem,RHS_Contribution,EquationId,CurrentProcessInfo);

                    //assemble the elemental contribution
                    AssembleRHS(b,RHS_Contribution,EquationId/*,lock_array*/);

                    // clean local elemental memory
                    pScheme->CleanMemory(*pElem);
                }
            }

            PartitionVector ConditionPartition;

            OpenMPUtils::DivideInPartitions(rConditions.size(),NumThreads,ConditionPartition);

            #pragma omp parallel
            {
                int k = OpenMPUtils::ThisThread();

                //contributions to the system
                LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

                Condition::EquationIdVectorType EquationId;

                ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();

                typename ConditionsArrayType::ptr_iterator CondBegin = rConditions.ptr_begin()+ConditionPartition[k];
                typename ConditionsArrayType::ptr_iterator CondEnd = rConditions.ptr_begin()+ConditionPartition[k+1];

                // assemble all conditions
                for (typename ConditionsArrayType::ptr_iterator pCond = CondBegin; pCond != CondEnd; pCond++)
                {
                    //calculate condition contribution
                    pScheme->Condition_Calculate_RHS_Contribution(*pCond,RHS_Contribution,EquationId,CurrentProcessInfo);

                    //assemble condition contribution
                    AssembleRHS(b,RHS_Contribution,EquationId/*,lock_array*/);
                }
            }

            KRATOS_CATCH("");
        }

        /// Identify Dofs and store pointers to them
        /**
         * Generates a list of the degrees of freedom in the model.
         * @param pScheme A pointer to the solution scheme
         * @param rModelPart A reference to the model part that contains the dofs
         */
        void SetUpDofSet(
                typename TSchemeType::Pointer pScheme,
                ModelPart& rModelPart)
        {
            KRATOS_TRY;

            std::cout << "Setting up degrees of freedom" << std::endl;
            //Gets the array of elements from the modeler
            ElementsArrayType& pElements = rModelPart.Elements();
            ConditionsArrayType& pConditions = rModelPart.Conditions();

            BaseType::mDofSet = DofsArrayType();

            ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();

            int NumThreads = OpenMPUtils::GetNumThreads();
            PartitionVector ElemPartition;
            PartitionVector CondPartition;

            OpenMPUtils::DivideInPartitions(pElements.size(),NumThreads,ElemPartition);
            OpenMPUtils::DivideInPartitions(pConditions.size(),NumThreads,CondPartition);

            std::vector< DofsArrayType > DofContainer(NumThreads);

            #pragma omp parallel
            {
                int k = OpenMPUtils::ThisThread();
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
            for (int k = 0; k < NumThreads ; k++)
                for( typename DofsArrayType::ptr_iterator itDof = DofContainer[k].ptr_begin();
                        itDof != DofContainer[k].ptr_end(); itDof++) {
                    BaseType::mDofSet.push_back(*itDof);
                }

            BaseType::mDofSet.Unique();

            //throws an execption if there are no Degrees of freedom involved in the analysis
            if (BaseType::mDofSet.size() == 0)
                KRATOS_ERROR(std::logic_error, "No degrees of freedom!", "");

            BaseType::mDofSetIsInitialized = true;
            mDofSetChanged = true;

            KRATOS_CATCH("");
        }

        /// Organise Dofs, separating fixed and free nodes
        /**
         * This function orders the degrees of freedom of the system.
         * To prepare for uncoupled solution, the numeration is assigned as follows:
         * free velocity dofs, free pressure dofs, fixed velocity dofs and fixed pressure dofs
         * This ordering allows us to set up separated dof counters for each dof type.
         * They will be used to keep track of the type (variable and fixity) of each dof
         *  during dof loops.
         * @param rModelPart A reference to the model part that contains the dofs
         */
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

        /// Initialize pointers to system vectors and matrices, check sizes
        /**
         *  Initialize pointers to the different vectors involved, check that they
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
                    mDofSetChanged == true || // if the dof set has changed
                    S.size1() == 0 || D.size1() == 0 || G.size1() == 0 ||
                     L.size1() == 0 || A.size1() == 0 ) //if the matrices are not initialized
            {
                S.resize(mVelFreeDofs, mVelFreeDofs, false);
                G.resize(mVelFreeDofs, mPressFreeDofs, false);
                D.resize(mPressFreeDofs, mVelFreeDofs, false);
                L.resize(mPressFreeDofs, mPressFreeDofs, false);

                ConstructMatrixStructure(S, D, G, L, rElements, rConditions, CurrentProcessInfo);

                A.resize(mPressFreeDofs, mPressFreeDofs, false);
                IDiagS.resize(mVelFreeDofs);

                AllocateSystemMatrix(A);
                ConstructSystemMatrix(A);
                mDofSetChanged = false;
            }
            else
            {
                // I do the check only for A, as the remaining matrices are private.
                // They are managed by this class, so they shouldn't change size spontaneously
                if (A.size1() != mPressFreeDofs || A.size2() != mPressFreeDofs )
                {
                    KRATOS_WATCH("it should not come here!!!!!!!! ... this is SLOW");

                    A.resize(mPressFreeDofs, mPressFreeDofs, false);
                    IDiagS.resize(mVelFreeDofs);

                    AllocateSystemMatrix(A);
                    ConstructSystemMatrix(A);
                    mDofSetChanged = false;
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
            KRATOS_CATCH("");
        }

        void FinalizeSolutionStep(
                ModelPart& rModelPart,
                TSystemMatrixType& A,
                TSystemVectorType& Dx,
                TSystemVectorType& b)
        {}

        /// Calculate Reactions
        void CalculateReactions(
                typename TSchemeType::Pointer pScheme,
                ModelPart& rModelPart,
                TSystemMatrixType& A,
                TSystemVectorType& Dx,
                TSystemVectorType& b)
        {
            //refresh RHS to have the correct reactions
            BuildRHS(pScheme, rModelPart, b);

            int systemsize = BaseType::mDofSet.size() - TSparseSpace::Size(*BaseType::mpReactionsVector);
            TSystemVectorType& ReactionsVector = *BaseType::mpReactionsVector;

            int NumThreads = OpenMPUtils::GetNumThreads();
            OpenMPUtils::PartitionVector DofPartition;
            OpenMPUtils::DivideInPartitions(BaseType::mDofSet.size(),NumThreads,DofPartition);

            //updating variables
            #pragma omp parallel
            {
                int i; // Position of the Dof in the Reaction vector
                int k = OpenMPUtils::ThisThread();
                typename DofsArrayType::ptr_iterator DofsBegin = BaseType::mDofSet.ptr_begin() + DofPartition[k];
                typename DofsArrayType::ptr_iterator DofsEnd = BaseType::mDofSet.ptr_begin() + DofPartition[k+1];

                for (typename DofsArrayType::ptr_iterator itDof = DofsBegin; itDof != DofsEnd; itDof++)
                {
                    if ( (*itDof)->IsFixed() )
                    {
                        i = (*itDof)->EquationId() - systemsize;
                        (*itDof)->GetSolutionStepReactionValue() = ReactionsVector[i];
                    }
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

        /// Free memory used by class members once no longer needed
        /**
         * this function is intended to be called at the end of the solution
         * step to clean up memory storage not needed
         */
        void Clear()
        {
            this->mDofSet = DofsArrayType();

            if (this->mpReactionsVector != NULL)
                TSparseSpace::Clear((this->mpReactionsVector));

            if (this->mpS != NULL)
                TSparseSpace::Clear((this->mpS));

            if (this->mpG != NULL)
                    TSparseSpace::Clear((this->mpG));

            if (this->mpD != NULL)
                    TSparseSpace::Clear((this->mpD));

            if (this->mpL != NULL)
                TSparseSpace::Clear((this->mpL));

            if (this->mpIDiagS != NULL)
                        TSparseSpace::Clear((this->mpIDiagS));

            if (this->GetEchoLevel() > 0) {
                std::cout << "ResidualBasedEliminationBuilderAndSolver Clear Function called" << std::endl;
            }
        }

        /// Force a recalculation of the graph for the pressure system matrix
        /**
         * Provides a way to modify the system matrix from outside the buider and solver.
         * Use carefully. Remember that this should be called before Build(), where
         * the matrix will be filled. If the DofSet changes between time steps, the
         * system matrix will be reformed internally, so this function should only be used
         * if the system matrix structure changes inside the time step.
         */
        inline void ReshapeSystemMatrix(TSystemMatrixType& A)
        {
            AllocateSystemMatrix(A);
            ConstructSystemMatrix(A);
        }

    protected:

        /// Compute graphs for the different matrices involved in the problem
        virtual void ConstructMatrixStructure(
                TSystemMatrixType& S,
                TSystemMatrixType& D,
                TSystemMatrixType& G,
                TSystemMatrixType& L,
                const ElementsContainerType& rElements,
                const ConditionsArrayType& rConditions,
                ProcessInfo& CurrentProcessInfo)
        {
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
            int NumTermsG = 0;
            int NumTermsD = 0;
            int NumTermsL = 0;

            for (std::size_t i = 0; i < indicesS.size(); i++)
                NumTermsS += indicesS[i].size();

            S.reserve(NumTermsS, false);

            for (std::size_t i = 0; i < indicesG.size(); i++)
                NumTermsG += indicesG[i].size();

            G.reserve(NumTermsG, false);

            for (std::size_t i = 0; i < indicesD.size(); i++)
                NumTermsD += indicesD[i].size();

            D.reserve(NumTermsD, false);

            for (std::size_t i = 0; i < indicesL.size(); i++)
                NumTermsL += indicesL[i].size();

            L.reserve(NumTermsL, false);

            // Create the matrix structure, filling it with zeros
            AllocateSpace(S, indicesS);
            AllocateSpace(G, indicesG);
            AllocateSpace(D, indicesD);
            AllocateSpace(L, indicesL);
        }

    private:

        /// Add a matrix position to the list (check that it hasn't been used before)
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

        /// Add terms from an elemental matrix to system matrices
#ifndef _OPENMP
        void Assemble(
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

#else
        void Assemble(
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

        /// Add terms from an elemental matrix to system matrices (LHS only)
        void AssembleLHS(
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

        /// Asseble local contribution to the right hand side vector
        /**
         * Asseble local contribution to the right hand side vector
         *
         * @param b Reference to RHS vector
         * @param RHS_Contribution local RHS contibution
         * @param EquationId Global ids of the local contribution
         * @param lock_array (only for OpenMP compilations) Vector containing a lock for each row of b
         */
        void AssembleRHS(
                TSystemVectorType& b,
                LocalSystemVectorType& RHS_Contribution,
                Element::EquationIdVectorType& EquationId)
        {
            unsigned int ContributionSize =  EquationId.size();

            for (unsigned int i = 0; i < ContributionSize; i++)
            {
                unsigned int Global_i = EquationId[i];
                if ( Global_i < BaseType::mEquationSystemSize )
                {
                    #pragma omp atomic
                    b[Global_i] += RHS_Contribution[i];
                }
                else if (BaseType::mCalculateReactionsFlag) // == true
                {
                    #pragma omp atomic
                    BaseType::mpReactionsVector->operator[](Global_i - BaseType::mEquationSystemSize) -= RHS_Contribution[i];
                }
            }
        }

        /**
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


        /// Identify non-zero tems in the system matrix
        void ConstructSystemMatrix(
                TSystemMatrixType& A)
        {
            // Retrieve matrices
            TSystemMatrixType& rG = *mpG;
            TSystemMatrixType& rD = *mpD;
            TSystemMatrixType& rL = *mpL;

            PartitionVector Partition;
            int NumThreads = OpenMPUtils::GetNumThreads();

            OpenMPUtils::DivideInPartitions(A.size1(),NumThreads,Partition);

            for( int k = 0 ; k < NumThreads ; k++)
            {
                // This code is serial, the pragma is here to ensure that each
                // row block is assigned to the processor that will fill it
                #pragma omp parallel
                if( OpenMPUtils::ThisThread() == k) {
                    boost::shared_ptr< IndexVector > pNext( new IndexVector(mPressFreeDofs) );
                    IndexVector& Next = *pNext; // Keeps track of which columns were filled
                    for(unsigned int m = 0; m < mPressFreeDofs; m++) Next[m] = -1;

                    std::size_t NumTerms = 0; // Full positions in a row
                    boost::shared_ptr< std::vector<unsigned int> > pUsedCols( new std::vector<unsigned int>);
                    std::vector<unsigned int>& UsedCols = *pUsedCols;
                    UsedCols.reserve(mPressFreeDofs);

                    for( int RowIndex = Partition[k] ;
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

        /// Compute the Pressure System Matrix
        /**
         *  Compute the System Matrix A = L - D*Inv(Diag(S))*G. The multiplication
         * is performed in random order, so each row will be stored in a temporary
         * variable, ordered and copied in input matrix A.
         */
        void CalculateSystemMatrix(
                TSystemMatrixType& A)
        {
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
                rIDiagS[i] = 1.0/rS(i,i);

            PartitionVector Partition;
            int NumThreads = OpenMPUtils::GetNumThreads();

            OpenMPUtils::DivideInPartitions(A.size1(),NumThreads,Partition);

            #pragma omp parallel
            {
                int k = OpenMPUtils::ThisThread();
                TSystemVectorPointerType pCurrentRow( new TSystemVectorType(mPressFreeDofs) );
                TSystemVectorType& CurrentRow = *pCurrentRow; // Values for current row
                for (unsigned int i = 0; i < mPressFreeDofs; i++) CurrentRow[i] = 0.0;

                boost::shared_ptr< IndexVector > pNext( new IndexVector(mPressFreeDofs) );
                IndexVector& Next = *pNext; // Keeps track of which columns were filled
                for(unsigned int m=0; m < mPressFreeDofs; m++) Next[m] = -1;

                std::size_t NumTerms = 0; // Full positions in a row
                boost::shared_ptr< std::vector<unsigned int> > pUsedCols( new std::vector<unsigned int>);
                std::vector<unsigned int>& UsedCols = *pUsedCols;
                UsedCols.reserve(mPressFreeDofs);

                for( int RowIndex = Partition[k] ;
                        RowIndex != Partition[k+1] ; RowIndex++ )
                {
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

                    // Identify full terms for ordering
                    for( std::size_t i = 0; i < Length; i++)
                    {
                        if( Next[head] != -1 )
                        {
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

                    // Fill matrix row, then clean temporary variables.
                    RowType RowA(A,RowIndex);
                    std::size_t n = 0;
                    unsigned int Col;

                    for( typename RowType::iterator ItA = RowA.begin(); ItA != RowA.end(); ItA++)
                    {
                        Col = UsedCols[n++];
                        *ItA = CurrentRow[Col];
                        CurrentRow[Col] = 0;
                    }
                    NumTerms = 0;
                    UsedCols.resize(0,false);
                }
            }

        }


        /// Helper function for Sytem matrix functions
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


        /// allocate space for member matrices
        void AllocateSpace(
//        unsigned int AllocateSpace( // Commented version checks for empty rows in the matrix
                TSystemMatrixType& A,
                std::vector< std::vector<std::size_t> >& indices)
        {
//            unsigned int NumEmptyRows = 0; // Keeps track of empty rows for error checking (they can be a symptom of insuficient boundary conditions)
            int NumThreads = OpenMPUtils::GetNumThreads();
            PartitionVector MatrixPartition;

            OpenMPUtils::DivideInPartitions(indices.size(),NumThreads,MatrixPartition);

            for( int k=0; k < NumThreads; k++ )
            {
                // First touch: Make the thread that will manipulate each partition
                // be the one that initializes it, so the relevant variables will
                // belong to it.
                #pragma omp parallel
                if( OpenMPUtils::ThisThread() == k )
                {
                    for( int i = MatrixPartition[k]; i < MatrixPartition[k+1]; i++ )
                    {
                        std::vector<std::size_t>& Row = indices[i];
                        std::sort(Row.begin(), Row.end());

                        for(std::vector<std::size_t>::iterator it= Row.begin(); it != Row.end() ; it++)
                        {
                            A.push_back(i,*it,0.00);
                        }
//                        if(Row.begin() == Row.end()) NumEmptyRows++;
                        Row.clear();
                    }
                }
            }
//            return NumEmptyRows;
        }

#ifdef _OPENMP

        /// y += A*x in parallel
        void Parallel_ProductAdd(
                const TSystemMatrixType& A,
                const TSystemVectorType& in,
                TSystemVectorType& out)
        {
            //create partition
            PartitionVector partition;
            int number_of_threads = OpenMPUtils::GetNumThreads();
            OpenMPUtils::DivideInPartitions(A.filled1()-1,number_of_threads,partition);

            #pragma omp parallel
            {
                int thread_id = omp_get_thread_num();
                unsigned int number_of_rows = partition[thread_id+1] - partition[thread_id];
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
        }

	/// calculates partial product for Parallel_ProductAdd
        void partial_product_add(
                unsigned int number_of_rows,
                typename compressed_matrix<TDataType>::index_array_type::const_iterator row_begin,
                typename compressed_matrix<TDataType>::index_array_type::const_iterator index2_begin,
                typename compressed_matrix<TDataType>::value_array_type::const_iterator value_begin,
                const TSystemVectorType& input_vec,
		unsigned int output_begin_index,
		TSystemVectorType& output_vec)
        {
	    unsigned int kkk = output_begin_index;
            typename TSystemMatrixType::index_array_type::const_iterator row_it = row_begin;
            for(unsigned int k = 0; k < number_of_rows; k++)
            {
                unsigned int RowBegin = *row_it;
                unsigned int RowEnd = *(++row_it);
                TDataType t = TDataType();

                for(unsigned int i = RowBegin; i < RowEnd; i++)
                    t += *value_begin++ *   ( input_vec[*index2_begin++]);

		output_vec[kkk++] += t;

            }
        }

#endif

        /// Set initial iterative solver tolerance using inexact Newton criteria
        void SetInitialTolerance(
                double RHSNorm,
                double& TolFactor)
        {
            TolFactor = mMaxTolFactor;
            (BaseType::mpLinearSystemSolver)->SetTolerance(mMaxTolFactor);
            std::cout << "Set iterative solver tolerance to " << TolFactor << std::endl;
        }

        /// Correct iterative solver tolerance using inexact Newton criteria
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

            if (TolFactor < mSmallTol) TolFactor = mSmallTol;
            (BaseType::mpLinearSystemSolver)->SetTolerance(TolFactor);
            std::cout << "Corrected iterative solver tolerance to " << TolFactor << std::endl;
        }

        // Total number of Dofs: BaseType::mEquationSystemSize;
        /// Position of 1st Free Pressure Dof
        unsigned int mVelFreeDofs;
        /// Position of 1st Fixed Pressure Dof
        unsigned int mVelFixedDofsEnd;

        unsigned int mPressFreeDofs;

        /// System Matrix for the momentum equation
        TSystemMatrixPointerType mpS;
        /// Discrete Divergence operator
        TSystemMatrixPointerType mpD;
        /// Discrete Gradient operator
        TSystemMatrixPointerType mpG;
        /// Stabilization term
        TSystemMatrixPointerType mpL;

        /// Inv(Diag(S)), stored as a vector
        TSystemVectorPointerType mpIDiagS;

        /// Linear solver for velocity system (pressure solver stored in base class)
        typename TLinearSolver::Pointer mpVelLinearSystemSolver;

        /// Flag for matrix reconstruction
        bool mDofSetChanged;

        /// Choice of method to ensure that velocity solution is divergence free
        /** Ensure that velocity is divergence free after each iteration
         * 0: Do not force (iteration process will eventually reach correct solution)
         * 1: By solving the system Diag(S)*dv = -G*dp
         * 2: By solving the full system S*dv = -G*dp
         */
        unsigned int mVelocityCorrection;

        // Variables for inexact Newton tolerance control
        /// Use Inexact Newton method
        bool mInexactNewton;
        /// used to compute maximum solution tolerance in Inexact Newton
        double mMaxTolFactor;
        /// Minimum solution tolerance for inexact Newton: 0.5*NonLinearTol
        double mSmallTol;
        /// Inexact Newton parameter
        double mGamma;

        /// Keeps track of the first iteration in each step
        bool mFirstIteration;

//        double mVelTolFactor;
        /// Inexact Newton pressure tolerance factor
        double mPressTolFactor;

//        double mLastVelRHSNorm;
        /// Inexact Newton pressure tolerance factor
        double mLastPressRHSNorm;
    };
}

#endif	/* KRATOS_PRESSURE_SPLITTING_BUILDER_AND_SOLVER_H */
