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
 * File:   trilinos_pressure_splitting_builder_and_solver.h
 * Author: jcotela
 *
 * Created on 19 April 2010, 16:02
 */

#ifndef KRATOS_TRILINOS_PRESSURE_SPLITTING_BUILDER_AND_SOLVER_H
#define	KRATOS_TRILINOS_PRESSURE_SPLITTING_BUILDER_AND_SOLVER_H

/* System includes */
#include <set>
#include <unistd.h>
#include <sys/types.h>

/* External includes */
#include "boost/smart_ptr.hpp"
#include <mpi.h>
// Trilinos Epetra
#include "Epetra_MpiComm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_FECrsGraph.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_LinearProblem.h"
// Trilinos EpetraExt
#include "EpetraExt_MatrixMatrix.h"
// Trilinos AztecOO
#include "AztecOO.h"
// Trilinos Amesos
#include "Amesos.h"

/* Project includes */
#include "includes/define.h"
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"

namespace Kratos
{
    /// An MPI builder and solver for the navier-stokes equations

    /**
     * A reimplementation of PressureSplittingBuilderAndSolver using Trilinos.
     * This builder and solver takes the monolithic linear system resulting from
     * the discretization and linearization of the incompressible Navier-Stokes
     * equations and divides it between the velocity and pressure degrees of
     * freedom. In the end, the problem is reduced to three smaller linear systems,
     * one of which can be explicitly inverted if a lumped mass matrix is used.
     * @see PressureSplittingBuilderAndSolver
     */
    template< class TSparseSpace,
    class TDenseSpace, //= DenseSpace<double>,
    class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
    >
    class TrilinosPressureSplittingBuilderAndSolver:
        public BuilderAndSolver< TSparseSpace,TDenseSpace,TLinearSolver >
    {
        /* Please note that this builder and solver splits the system matrix
         * in 4 blocks and stores them separately. The system matrix A used as
         * input (to preserve the same function signatures as other builder and
         * solvers) is used here to store the matrix of the pressure system.
         */
    public:

        KRATOS_CLASS_POINTER_DEFINITION(TrilinosPressureSplittingBuilderAndSolver);

        // Type Definitions
        typedef BuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

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

        // Life Cycle

        /// Constructor
        /**
         * Creates a TrilinosPressureSplittingBuilderAndSolver instance
         * @param Comm An MPI communicator for the model
         * @param guess_row_size Estimate of the number of non-zero terms in each matrix row
         * @param pVelocityLinearSystemSolver Pointer to the linear solver for the velocity system
         * @param pPressureLinearSystemSolver Pointer to the linear solver for the pressure system
         * @param VelocityCorrection If > 0, explicitly solve the velocity to be divergence-free at each step (other
         * @param UseInexactNewton If true, dynamically set the linear iterative solver tolerance for the pressure syste
         * @param NonLinearTol Only used if InexactNewton == true, otherwise the solver will use it's own tolerance
         * @param MaxTolFactor Inexact Newton parameter (used to set an upper bound for tolerance)
         * @param Gamma Inexact Newton parameter
         */
        TrilinosPressureSplittingBuilderAndSolver(Epetra_MpiComm& Comm,
                                                  int guess_row_size,
                                                  typename TLinearSolver::Pointer pVelocityLinearSystemSolver,
                                                  typename TLinearSolver::Pointer pPressureLinearSystemSolver,
                                                  unsigned int VelocityCorrection, // If > 0, explicitly solve the velocity to be divergence-free at each step
                                                  bool UseInexactNewton, // If true, dynamically set the linear iterative solver tolerance for the pressure system
                                                  double NonLinearTol = 1e-3, // Only used if InexactNewton == true, otherwise the solver will use it's own tolerance
                                                  double MaxTolFactor = 0.1,
                                                  double Gamma = 0.9) :
        BuilderAndSolver< TSparseSpace, TDenseSpace, TLinearSolver >(pPressureLinearSystemSolver),
        mpVelocityLinearSystemSolver(pVelocityLinearSystemSolver),
        mrComm(Comm),
        mRowSizeGuess(guess_row_size),
        mDofSetChanged(true),
        mVelocityCorrection(VelocityCorrection),
        mInexactNewton(UseInexactNewton),
        mMaxTolFactor(MaxTolFactor),
        mGamma(Gamma),
        mFirstIteration(true),
        mPressTolFactor(0),
        mLastPressRHSNorm(0)
        {
            mSmallTol = 0.5 * NonLinearTol;
        }

        /// Destructor
        virtual ~TrilinosPressureSplittingBuilderAndSolver() { }

        /// Build System
        void Build(typename TSchemeType::Pointer pScheme,
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
            TSparseSpace::SetToZero(*(BaseType::mpReactionsVector));

            // Reset internally stored matrices
            TSparseSpace::SetToZero(*mpS);
            TSparseSpace::SetToZero(*mpG);
            TSparseSpace::SetToZero(*mpD);
            TSparseSpace::SetToZero(*mpL);

            // Define contributions to the system
            LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0, 0);
            LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

            // Will store the position of each Dof in the system
            Element::EquationIdVectorType EquationIds;

            ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();

            // Assemble contributions from elements
            for (typename ElementsArrayType::ptr_iterator pElem = rElements.ptr_begin();
                 pElem != rElements.ptr_end(); pElem++)
            {
                // Get Elemental Contributions
                pScheme->CalculateSystemContributions(*pElem, LHS_Contribution,
                                                      RHS_Contribution, EquationIds, CurrentProcessInfo);

                //assemble the elemental contribution
                Assemble(b, LHS_Contribution, RHS_Contribution, EquationIds);

                // clean local elemental memory
                pScheme->CleanMemory(*pElem);
            }

            LHS_Contribution.resize(0, 0, false);
            RHS_Contribution.resize(0, false);

            // assemble contributions from conditions
            for (typename ConditionsArrayType::ptr_iterator pCond = rConditions.ptr_begin();
                 pCond != rConditions.ptr_end(); pCond++)
            {
                // Get condition Contributions
                pScheme->Condition_CalculateSystemContributions(*pCond, LHS_Contribution,
                                                                RHS_Contribution, EquationIds, CurrentProcessInfo);

                // Assemble condition contribution
                Assemble(b, LHS_Contribution, RHS_Contribution, EquationIds);
            }


            // Finalise Assembly
            b.GlobalAssemble();
            mpS->GlobalAssemble();
            mpL->GlobalAssemble();
            // Non-square matrices have to be given their DomainMap and RangeMap explicitly
            mpG->GlobalAssemble(mpL->RowMap(), mpS->RowMap());
            mpD->GlobalAssemble(mpS->RowMap(), mpL->RowMap());

            // Assemble the pressure system matrix. Assumptions:
            // 1- A was created with the proper graph and has Filled() == true (true if created using ResizeAndInitializeVectors())
            // 2- A is filled with zeros (we rely on the strategy to do that)
            Epetra_CrsMatrix* pScaledG = new Epetra_CrsMatrix(*mpG);
            Epetra_Vector* pDiagS = new Epetra_Vector(mpS->RowMap(), false);
            mpS->ExtractDiagonalCopy(*pDiagS);

            // A = L - D*(1/Diag(S))*G
            // Start by computing ScaledG = (1/Diag(S))*G
            pDiagS->Reciprocal(*pDiagS); // Inv{Diag(S)}[i] <- 1/Diag(S)[i]

            // Store Inv{Diag(S)} * G
            pScaledG->LeftScale(*pDiagS);

            // A <- D*ScaledG
            EpetraExt::MatrixMatrix::Multiply(*mpD, false, *pScaledG, false, A, false);

            // A <- L - A
            EpetraExt::MatrixMatrix::Add(*mpL, false, 1.0, A, -1.0); // bool is for transpose of L

            // If we intend to use Inv{Diag(S)} again in the solve phase, store it
            if (mVelocityCorrection == 1)
            {
                boost::shared_ptr<Epetra_Vector> pNewDiagS(new Epetra_Vector(*pDiagS));
                mpIDiagS.swap(pNewDiagS);
            }

            delete pDiagS;
            delete pScaledG;
            pDiagS = 0;
            pScaledG = 0;

            KRATOS_CATCH("")
        }

        void BuildLHS(typename TSchemeType::Pointer pScheme,
                      ModelPart& rModelPart,
                      TSystemMatrixType& A)
        {
            KRATOS_TRY
            if (!pScheme)
                KRATOS_ERROR(std::runtime_error, "No scheme provided!", "");

            // Get elements and conditions
            ElementsArrayType& rElements = rModelPart.Elements();
            ConditionsArrayType& rConditions = rModelPart.Conditions();

            //            // resetting to zero the vector of reactions
            //            TSparseSpace::SetToZero( *(BaseType::mpReactionsVector) );

            // Reset internally stored matrices
            TSparseSpace::SetToZero(*mpS);
            TSparseSpace::SetToZero(*mpG);
            TSparseSpace::SetToZero(*mpD);
            TSparseSpace::SetToZero(*mpL);

            // Define contributions to the system
            LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0, 0);

            // Will store the position of each Dof in the system
            Element::EquationIdVectorType EquationIds;

            ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();

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

            // assemble contributions from conditions
            for (typename ConditionsArrayType::ptr_iterator pCond = rConditions.ptr_begin();
                 pCond != rConditions.ptr_end(); pCond++)
            {
                // Get condition Contributions
                pScheme->Condition_Calculate_LHS_Contribution(*pCond, LHS_Contribution,
                                                              EquationIds, CurrentProcessInfo);

                // Assemble condition contribution
                AssembleLHS(LHS_Contribution, EquationIds);
            }


            // Finalize Assembly
            mpS->GlobalAssemble();
            mpL->GlobalAssemble();
            // Non-square matrices have to be given their DomainMap and RangeMap explicitly
            mpG->GlobalAssemble(mpL->RowMap(), mpS->RowMap());
            mpD->GlobalAssemble(mpS->RowMap(), mpL->RowMap());

            // Assemble the pressure system matrix. Assumptions:
            // 1- A was created with the proper graph and has Filled() == true (true if created using ResizeAndInitializeVectors())
            // 2- A is filled with zeros (we rely on the strategy to do that)
            Epetra_CrsMatrix* pScaledG = new Epetra_CrsMatrix(*mpG);
            Epetra_Vector* pDiagS = new Epetra_Vector(mpS->RowMap(), false);
            mpS->ExtractDiagonalCopy(*pDiagS);

            // A = L - D*(1/Diag(S))*G
            // Start by computing ScaledG = (1/Diag(S))*G
            pDiagS->Reciprocal(*pDiagS); // Inv{Diag(S)}[i] <- 1/Diag(S)[i]

            // Store Inv{Diag(S)} * G
            pScaledG->LeftScale(*pDiagS);

            // A <- D*ScaledG
            EpetraExt::MatrixMatrix::Multiply(*mpD, false, *pScaledG, false, A, false);

            // A <- L - A
            EpetraExt::MatrixMatrix::Add(*mpL, false, 1.0, A, -1.0); // bool is for transpose of L

            // If we intend to use Inv{Diag(S)} again in the solve phase, store it
            if (mVelocityCorrection == 1)
            {
                boost::shared_ptr<Epetra_Vector> pNewDiagS(new Epetra_Vector(*pDiagS));
                mpIDiagS.swap(pNewDiagS);
            }

            delete pDiagS;
            delete pScaledG;
            pDiagS = 0;
            pScaledG = 0;

            KRATOS_CATCH("")
        }

        /// Call Solver
        void SystemSolve(TSystemMatrixType& A,
                         TSystemVectorType& Dx,
                         TSystemVectorType& b)
        {
            KRATOS_TRY;

            double NormB;

            if (TSparseSpace::Size(b) != 0)
            {
                NormB = TSparseSpace::TwoNorm(b);
            }
            else
                NormB = 0;

            if (NormB != 0.0)
            {

                /* Solve current iteration in several steps */

                // Maps for temporary distributed vectors
                const Epetra_Map& UMap = mpS->RowMap();
                const Epetra_Map& PMap = mpL->RowMap();

                // Local variable counts and indices
                int LocalDofNum = mLocalVelDofNum + mLocalPressDofNum;

                int* VelIndices = new int[mLocalVelDofNum];
                for (int i = 0; i < mLocalVelDofNum; i++)
                    VelIndices[i] = mLocalVelBegin + i;

                int* PressIndices = new int[mLocalPressDofNum];
                for (int i = 0; i < mLocalPressDofNum; i++)
                    PressIndices[i] = mLocalPressBegin + i;

                // Extract the coupled RHS to build its velocity and pressure parts
                double* RHSValues = new double[LocalDofNum];
                b.ExtractCopy(RHSValues, LocalDofNum);

                // Build decoupled velocity and pressure vectors
                TSystemVectorPointerType pDVel(new TSystemVectorType(UMap));
                TSystemVectorPointerType pVelRHS(new TSystemVectorType(UMap));
                pVelRHS->ReplaceGlobalValues(mLocalVelDofNum, VelIndices, RHSValues, 0);

                TSystemVectorPointerType pDPress(new TSystemVectorType(PMap));
                TSystemVectorPointerType pPressRHS(new TSystemVectorType(PMap));

                TSystemMatrixType& rS = *mpS; // Create a reference to the velocity system matrix

                // Solution Process

                // 1. Compute intermediate velocity
                mpVelocityLinearSystemSolver->Solve(rS, *pDVel, *pVelRHS);

                // 2. Compute Pressure Variation
                TSystemVectorPointerType pTemp(new TSystemVectorType(PMap));
                TSparseSpace::Mult(*mpD, *pDVel, *pTemp);
                double* pPressBuffer = new double[mLocalPressDofNum];
                pTemp->ExtractCopy(pPressBuffer, mLocalPressDofNum);
                for (int i = 0; i < mLocalPressDofNum; i++)
                    pPressBuffer[i] = RHSValues[mLocalVelDofNum + i] - pPressBuffer[i];
                pPressRHS->ReplaceGlobalValues(mLocalPressDofNum, PressIndices, pPressBuffer, 0);

                // Update linear tolerance (for Inexact Newton-Raphson)
                if (mInexactNewton == true)
                {
                    double PressRHSNorm = TSparseSpace::TwoNorm(*pPressRHS);
                    if (mFirstIteration == true)
                    {
                        SetInitialTolerance(PressRHSNorm, mPressTolFactor);
                    }
                    else
                    {
                        UpdateTolerance(mLastPressRHSNorm, PressRHSNorm, mPressTolFactor);
                    }
                    mLastPressRHSNorm = PressRHSNorm;
                }

                // Solve the system
                BaseType::mpLinearSystemSolver->Solve(A, *pDPress, *pPressRHS);

                // 3. Determine End of Step velocity
                double* pVelBuffer = new double[mLocalVelDofNum];
                if (mVelocityCorrection == 0)
                {
                    pDVel->ExtractCopy(pVelBuffer, mLocalVelDofNum);
                }
                else
                {
                    TSparseSpace::Mult(*mpG, *pDPress, *pVelRHS);

                    if (mVelocityCorrection == 1)
                    {
                        // DVel = DVel - Inv{Diag(S)}*VelRHS
                        pDVel->Multiply(-1.0, *mpIDiagS, *pVelRHS, 1.0);

                        pDVel->ExtractCopy(pVelBuffer, mLocalVelDofNum);
                    }
                    else if (mVelocityCorrection == 2)
                    {
                        TSystemVectorPointerType pVelUpdate(new TSystemVectorType(UMap));

                        mpVelocityLinearSystemSolver->Solve(rS, *pVelUpdate, *pVelRHS);
                        pDVel->Update(-1.0, *pVelUpdate, 1.0); // DVel = 1.0*DVel -1.0*pVelUpdate

                        pDVel->ExtractCopy(pVelBuffer, mLocalVelDofNum);
                    }
                }

                // Preconditioner
                //                A = *mpL - A;
                //                noalias(rPressRHS) = prod(A,rPressRHS);


                int* GlobalIndices = new int[LocalDofNum];
                for (int i = 0; i < mLocalVelDofNum; i++)
                {
                    RHSValues[i] = pVelBuffer[i];
                    GlobalIndices[i] = VelIndices[i];
                }
                pDPress->ExtractCopy(pPressBuffer, mLocalPressDofNum);
                for (int i = 0, j = mLocalVelDofNum; i < mLocalPressDofNum; i++, j++)
                {
                    RHSValues[j] = pPressBuffer[i];
                    GlobalIndices[j] = mVelFreeDofs + PressIndices[i];
                }

                // Copy the solution to output variable
                Dx.ReplaceGlobalValues(LocalDofNum, GlobalIndices, RHSValues, 0);

                if (mFirstIteration == true) mFirstIteration = false;

                // Clean allocated space
                delete [] VelIndices;
                delete [] PressIndices;
                delete [] RHSValues;
                delete [] pPressBuffer;
                delete [] pVelBuffer;
                delete [] GlobalIndices;
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

        void BuildAndSolve(typename TSchemeType::Pointer pScheme,
                           ModelPart& rModelPart,
                           TSystemMatrixType& A,
                           TSystemVectorType& Dx,
                           TSystemVectorType& b)
        {
            KRATOS_TRY

            boost::timer building_time;

            int rank = rModelPart.GetCommunicator().MyPID();

            Build(pScheme, rModelPart, A, b);

            if (BaseType::GetEchoLevel() > 0)
            {
                if (rank == 0) std::cout << "Building Time : " << building_time.elapsed() << std::endl;
            }

            //does nothing...dirichlet conditions are naturally dealt with in defining the residual
            ApplyDirichletConditions(pScheme, rModelPart, A, Dx, b);

            if (BaseType::GetEchoLevel() == 3)
            {
                if (rank == 0)
                {
                    std::cout << "before the solution of the system" << std::endl;
                    std::cout << "System Matrix = " << A << std::endl;
                    std::cout << "unknowns vector = " << Dx << std::endl;
                    std::cout << "RHS vector = " << b << std::endl;
                }
            }

            boost::timer solve_time;

            SystemSolve(A, Dx, b);

            if (BaseType::GetEchoLevel() > 0)
            {
                if (rank == 0) std::cout << "System Solve Time : " << solve_time.elapsed() << std::endl;
            }
            if (BaseType::GetEchoLevel() == 3)
            {
                if (rank == 0)
                {
                    std::cout << "after the solution of the system" << std::endl;
                    std::cout << "System Matrix = " << A << std::endl;
                    std::cout << "unknowns vector = " << Dx << std::endl;
                    std::cout << "RHS vector = " << b << std::endl;
                }
            }

            KRATOS_CATCH("")
        }

        /// Solve System for updated RHS
        void BuildRHSAndSolve(typename TSchemeType::Pointer pScheme,
                              ModelPart& rModelPart,
                              TSystemMatrixType& A,
                              TSystemVectorType& Dx,
                              TSystemVectorType& b)
        {
            KRATOS_TRY

            BuildRHS(pScheme, rModelPart, b);
            SystemSolve(A, Dx, b);

            KRATOS_CATCH("");
        }

        /// Build RHS only
        void BuildRHS(typename TSchemeType::Pointer pScheme,
                      ModelPart& rModelPart,
                      TSystemVectorType& b)
        {
            KRATOS_TRY
            if (!pScheme)
                KRATOS_ERROR(std::runtime_error, "No scheme provided!", "");

            // Get elements and conditions
            ElementsArrayType& rElements = rModelPart.Elements();
            ConditionsArrayType& rConditions = rModelPart.Conditions();

            // resetting to zero the vector of reactions
            TSparseSpace::SetToZero(*(BaseType::mpReactionsVector));

            LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

            // Will store the position of each Dof in the system
            Element::EquationIdVectorType EquationIds;

            ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();

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

            // assemble contributions from conditions
            for (typename ConditionsArrayType::ptr_iterator pCond = rConditions.ptr_begin();
                 pCond != rConditions.ptr_end(); pCond++)
            {
                // Get condition Contributions
                pScheme->Condition_Calculate_RHS_Contribution(*pCond, RHS_Contribution,
                                                              EquationIds, CurrentProcessInfo);

                // Assemble condition contribution
                AssembleRHS(b, RHS_Contribution, EquationIds);
            }


            // Finalize Assembly
            b.GlobalAssemble();

            KRATOS_CATCH("")
        }

        /// Identify Dofs and store pointers to them
        void SetUpDofSet(typename TSchemeType::Pointer pScheme,
                         ModelPart& rModelPart)
        {
            KRATOS_TRY;

            //Gets the array of elements from the modeler
            ElementsArrayType& pElements = rModelPart.GetCommunicator().LocalMesh().Elements();
            ConditionsArrayType& pConditions = rModelPart.GetCommunicator().LocalMesh().Conditions();

            DofsArrayType Doftemp;
            BaseType::mDofSet = DofsArrayType();

            Element::DofsVectorType ElementalDofList;

            ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();

            for (typename ElementsArrayType::ptr_iterator it = pElements.ptr_begin(); it != pElements.ptr_end(); ++it)
            {
                // gets list of Dof involved on every element
                pScheme->GetElementalDofList(*it, ElementalDofList, CurrentProcessInfo);

                for (typename Element::DofsVectorType::iterator i = ElementalDofList.begin(); i != ElementalDofList.end(); ++i)
                {
                    Doftemp.push_back(*i);
                }
            }

            //taking in account conditions
            for (typename ConditionsArrayType::ptr_iterator it = pConditions.ptr_begin(); it != pConditions.ptr_end(); ++it)
            {
                // gets list of Dof involved on every element
                pScheme->GetConditionDofList(*it, ElementalDofList, CurrentProcessInfo);

                for (typename Element::DofsVectorType::iterator i = ElementalDofList.begin(); i != ElementalDofList.end(); ++i)
                {
                    //mDofSet.push_back(*i);
                    Doftemp.push_back(*i);
                }
            }

            Doftemp.Unique();


            BaseType::mDofSet = Doftemp;

            //throws an execption if there are no Degrees of freedom involved in the analysis
            if (BaseType::mDofSet.size() == 0)
                KRATOS_ERROR(std::logic_error, "No degrees of freedom!", "");

            BaseType::mDofSetIsInitialized = true;
            mDofSetChanged = true;

            KRATOS_CATCH("")
        }

        /// Organise Dofs, separating fixed and free nodes
        void SetUpSystem(ModelPart& rModelPart)
        {
            KRATOS_TRY;

            int Rank = mrComm.MyPID();

            // Count local Local Dofs
            mLocalVelDofNum = 0;
            mLocalPressDofNum = 0;
            int VelFixedCount = 0;
            int PressFixedCount = 0;

            for (typename DofsArrayType::iterator itDof = BaseType::mDofSet.begin();
                 itDof != BaseType::mDofSet.end(); itDof++)
            {
                // I kwow all Dofs in my partition plus some dofs in the shared boundary.
                // Those nodes belong to other partitions, so I won't count them
                if (itDof->GetSolutionStepValue(PARTITION_INDEX) == Rank)
                {
                    KeyType CurrVar = itDof->GetVariable().Key(); // Get the Dof's variable
                    if ((CurrVar == VELOCITY_X) || (CurrVar == VELOCITY_Y)
                        || (CurrVar == VELOCITY_Z))
                    {
                        if (itDof->IsFree())
                            mLocalVelDofNum++;
                        else
                            VelFixedCount++;
                    }
                    else // "Not-a-velocity" Dof (it is assumed to represent pressure)
                    {
                        if (itDof->IsFree())
                            mLocalPressDofNum++;
                        else
                            PressFixedCount++;
                    }
                }
            }

            // Compute ofsets for Dof indexes
            int LocalCounts[4] = {mLocalVelDofNum, VelFixedCount, mLocalPressDofNum, PressFixedCount};
            int GlobalCounts[4]; // Total count of each Dof type
            int GlobalOffsets[4]; // First position of each Dof type for this process

            MPI_Allreduce(&LocalCounts, &GlobalCounts, 4, MPI_INT, MPI_SUM, mrComm.Comm());
            MPI_Scan(&LocalCounts, &GlobalOffsets, 4, MPI_INT, MPI_SUM, mrComm.Comm());
            // NOTE: The output of Exscan is undefined for process 0.
            // For us, its proper output would be an array of zeros

            // Store some useful quantities
            mVelFreeDofs = GlobalCounts[0];
            mPressFreeDofs = GlobalCounts[2];
            BaseType::mEquationSystemSize = mVelFreeDofs + mPressFreeDofs;

            unsigned int VelFreeIndex, PressFreeIndex, VelFixedIndex, PressFixedIndex;

            // Each process will number its Dofs starting from:
            VelFreeIndex = GlobalOffsets[0]-LocalCounts[0];
            PressFreeIndex = GlobalOffsets[2]-LocalCounts[2] + mVelFreeDofs;
            VelFixedIndex = GlobalOffsets[1]-LocalCounts[1] + BaseType::mEquationSystemSize;
            PressFixedIndex = GlobalOffsets[3]-LocalCounts[3] + BaseType::mEquationSystemSize + GlobalCounts[2];

            mLocalVelBegin = VelFreeIndex;
            mLocalPressBegin = PressFreeIndex - mVelFreeDofs;

            // Assign Ids
            for (typename DofsArrayType::iterator itDof = BaseType::mDofSet.begin();
                 itDof != BaseType::mDofSet.end(); itDof++)
                if (itDof->GetSolutionStepValue(PARTITION_INDEX) == Rank)
                {
                    KeyType CurrVar = itDof->GetVariable().Key(); // Get the Dof's variable
                    if ((CurrVar == VELOCITY_X) || (CurrVar == VELOCITY_Y)
                        || (CurrVar == VELOCITY_Z))
                    {
                        if (itDof->IsFree())
                            itDof->SetEquationId(VelFreeIndex++);
                        else
                            itDof->SetEquationId(VelFixedIndex++);
                    }
                    else // "Not-a-velocity" Dof (it is assumed to represent pressure)
                    {
                        if (itDof->IsFree())
                            itDof->SetEquationId(PressFreeIndex++);
                        else
                            itDof->SetEquationId(PressFixedIndex++);
                    }
                }

            rModelPart.GetCommunicator().SynchronizeDofs();

            KRATOS_CATCH("");
        }

        void ResizeAndInitializeVectors(TSystemMatrixPointerType& pA,
                                        TSystemVectorPointerType& pDx,
                                        TSystemVectorPointerType& pb,
                                        ElementsArrayType& rElements,
                                        ConditionsArrayType& rConditions,
                                        ProcessInfo& CurrentProcessInfo)
        {
            KRATOS_TRY
            if (this->GetEchoLevel() > 1)
                std::cout << "entering ResizeAndInitializeVectors" << std::endl;

            // Resizing system matrices
            if (BaseType::GetReshapeMatrixFlag() == true ||
                mDofSetChanged == true ||
                pA == NULL || TSparseSpace::Size1(*pA) == 0 ||
                mpS == NULL || TSparseSpace::Size1(*mpS) == 0 ||
                mpG == NULL || TSparseSpace::Size1(*mpG) == 0 ||
                mpD == NULL || TSparseSpace::Size1(*mpD) == 0 ||
                mpL == NULL || TSparseSpace::Size1(*mpL) == 0)
            {
                int TempVelSize = (mLocalVelDofNum > 1000) ? mLocalVelDofNum : 1000;
                int TempPressSize = (mLocalPressDofNum > 500) ? mLocalPressDofNum : 500;

                int* LocalVelIndices = new int[TempVelSize];
                int* LocalPressIndices = new int[TempPressSize];

                // Genrate a map for velocity Dofs and another one for pressure
                for (int i = 0, j = mLocalVelBegin; i < mLocalVelDofNum; i++, j++)
                    LocalVelIndices[i] = j;
                Epetra_Map UMap(mVelFreeDofs, mLocalVelDofNum, LocalVelIndices, 0, mrComm);

                for (int i = 0, j = mLocalPressBegin; i < mLocalPressDofNum; i++, j++)
                    LocalPressIndices[i] = j;
                Epetra_Map PMap(mPressFreeDofs, mLocalPressDofNum, LocalPressIndices, 0, mrComm);

                // Generate two extra maps containing all Dofs of the same type in a single
                // processor. They will be used as ColMaps for rectangular matrices
                int* ColVelIndices = new int[mVelFreeDofs];
                int* ColPressIndices = new int[mPressFreeDofs];

                for (unsigned int i = 0; i < mVelFreeDofs; i++)
                    ColVelIndices[i] = i;

                for (unsigned int i = 0; i < mPressFreeDofs; i++)
                    ColPressIndices[i] = i;

                Epetra_Map ColUMap(mVelFreeDofs, mVelFreeDofs, ColVelIndices, 0, mrComm);
                Epetra_Map ColPMap(mPressFreeDofs, mPressFreeDofs, ColPressIndices, 0, mrComm);

                // Create and fill the graph for the matrices
                Epetra_FECrsGraph SGraph(Copy, UMap, mRowSizeGuess);
                Epetra_FECrsGraph GGraph(Copy, UMap, ColPMap, mRowSizeGuess);
                Epetra_FECrsGraph DGraph(Copy, PMap, ColUMap, mRowSizeGuess);
                Epetra_FECrsGraph LGraph(Copy, PMap, mRowSizeGuess);

                delete [] LocalVelIndices; LocalVelIndices = 0;
                delete [] LocalPressIndices; LocalPressIndices = 0;
                delete [] ColVelIndices; ColVelIndices = 0;
                delete [] ColPressIndices; ColPressIndices = 0;

                // Generate Epetra_Graphs containing the structure of each matrix
                int* VelIds = new int[100];
                int* PressIds = new int[100];

                Element::EquationIdVectorType EquationId;

                // Fill graphs with element contributions
                for (typename ElementsArrayType::ptr_iterator pElem = rElements.ptr_begin();
                     pElem != rElements.ptr_end(); pElem++)
                {
                    // This should go through the Scheme !
                    (*pElem)->EquationIdVector(EquationId, CurrentProcessInfo);

                    unsigned int NumVelTerms(0), NumPressTerms(0);

                    for (unsigned int k = 0; k < EquationId.size(); k++)
                    {
                        if (EquationId[k] < mVelFreeDofs)
                        {
                            VelIds[NumVelTerms] = EquationId[k];
                            NumVelTerms++;
                        }
                        else if (EquationId[k] < BaseType::mEquationSystemSize)
                        {
                            PressIds[NumPressTerms] = EquationId[k] - mVelFreeDofs;
                            NumPressTerms++;
                        }
                    }

                    // Write terms
                    if (NumVelTerms != 0)
                    {
                        SGraph.InsertGlobalIndices(NumVelTerms, VelIds, NumVelTerms, VelIds);
                    }
                    if (NumVelTerms != 0 && NumPressTerms != 0)
                    {
                        GGraph.InsertGlobalIndices(NumVelTerms, VelIds, NumPressTerms, PressIds);
                        DGraph.InsertGlobalIndices(NumPressTerms, PressIds, NumVelTerms, VelIds);
                    }
                    if (NumPressTerms != 0)
                    {
                        LGraph.InsertGlobalIndices(NumPressTerms, PressIds, NumPressTerms, PressIds);
                    }
                }

                // Fill graphs with condition contributions
                for (typename ConditionsArrayType::ptr_iterator pCond = rConditions.ptr_begin();
                     pCond != rConditions.ptr_end(); pCond++)
                {
                    // Get condition Contributions
                    (*pCond)->EquationIdVector(EquationId, CurrentProcessInfo);

                    unsigned int NumVelTerms(0), NumPressTerms(0);

                    for (unsigned int k = 0; k < EquationId.size(); k++)
                    {
                        if (EquationId[k] < mVelFreeDofs)
                        {
                            VelIds[NumVelTerms] = EquationId[k];
                            NumVelTerms++;
                        }
                        else if (EquationId[k] < BaseType::mEquationSystemSize)
                        {
                            PressIds[NumPressTerms] = EquationId[k] - mVelFreeDofs;
                            NumPressTerms++;
                        }
                    }

                    // Write terms
                    if (NumVelTerms != 0)
                    {
                        SGraph.InsertGlobalIndices(NumVelTerms, VelIds, NumVelTerms, VelIds);
                    }
                    if (NumVelTerms != 0 && NumPressTerms != 0)
                    {
                        GGraph.InsertGlobalIndices(NumVelTerms, VelIds, NumPressTerms, PressIds);
                        DGraph.InsertGlobalIndices(NumPressTerms, PressIds, NumVelTerms, VelIds);
                    }
                    if (NumPressTerms != 0)
                    {
                        LGraph.InsertGlobalIndices(NumPressTerms, PressIds, NumPressTerms, PressIds);
                    }
                }

                // Finish completed graphs
                // Note: Non-square matrices need their DomainMap and RangeMap as input
                int GraphError = 0;
                GraphError += SGraph.GlobalAssemble(true);
                GraphError += GGraph.GlobalAssemble(PMap, UMap, true);
                GraphError += DGraph.GlobalAssemble(UMap, PMap, true);
                GraphError += LGraph.GlobalAssemble(true);

                if (GraphError != 0) KRATOS_ERROR(std::logic_error, "Epetra failure during matrix inicialization", "");

                // Create & store matrices
                TSystemMatrixPointerType pNewS = TSystemMatrixPointerType(new TSystemMatrixType(Copy, SGraph));
                TSystemMatrixPointerType pNewG = TSystemMatrixPointerType(new TSystemMatrixType(Copy, GGraph));
                TSystemMatrixPointerType pNewD = TSystemMatrixPointerType(new TSystemMatrixType(Copy, DGraph));
                TSystemMatrixPointerType pNewL = TSystemMatrixPointerType(new TSystemMatrixType(Copy, LGraph));

                mpS.swap(pNewS);
                mpG.swap(pNewG);
                mpD.swap(pNewD);
                mpL.swap(pNewL);

                // Finally, define the pressure system matrix's graph from the shapes of the different matrices
                boost::shared_ptr<Epetra_CrsGraph> pGraphA;

                ConstructSystemMatrixGraph(pGraphA);
                mDofSetChanged = false;

                // Create an empty pressure system matrix at the location pointed by pA
                TSystemMatrixPointerType pNewA = TSystemMatrixPointerType(new TSystemMatrixType(Copy, *pGraphA));
                pA.swap(pNewA);

                delete [] VelIds;
                VelIds = 0;
                delete [] PressIds;
                PressIds = 0;

            }
            else if (TSparseSpace::Size1(*pA) == 0 || TSparseSpace::Size1(*pA) != BaseType::mEquationSystemSize ||
                     TSparseSpace::Size2(*pA) != BaseType::mEquationSystemSize)
                KRATOS_ERROR(std::logic_error, "ResizeAndInitialize Error: Unexpected change of matrix dimensions!", "");


            // System Vectors
            if (pb == NULL || TSparseSpace::Size(*pb) != BaseType::mEquationSystemSize ||
                pDx == NULL || TSparseSpace::Size(*pDx) != BaseType::mEquationSystemSize ||
                BaseType::mpReactionsVector == NULL)
            {
                KRATOS_WATCH("System Vectors")
                int LocalDofNum = mLocalVelDofNum + mLocalPressDofNum;

                int* LocalIndices = new int[LocalDofNum];

                for (int i = 0, j = mLocalVelBegin; i < mLocalVelDofNum; i++, j++)
                    LocalIndices[i] = j;
                for (int i = mLocalVelDofNum, j = mVelFreeDofs + mLocalPressBegin; i < LocalDofNum; i++, j++)
                    LocalIndices[i] = j;

                Epetra_Map VectorMap(BaseType::mEquationSystemSize, LocalDofNum, LocalIndices, 0, mrComm);

                TSystemVectorPointerType pNewb = TSystemVectorPointerType(new TSystemVectorType(VectorMap));
                pb.swap(pNewb);
                KRATOS_WATCH(pb->MyLength())

                TSystemVectorPointerType pNewDx = TSystemVectorPointerType(new TSystemVectorType(VectorMap));
                pDx.swap(pNewDx);

                TSystemVectorPointerType pNewReactionsVector = TSystemVectorPointerType(new TSystemVectorType(VectorMap));
                BaseType::mpReactionsVector.swap(pNewReactionsVector);

                delete [] LocalIndices;
                LocalIndices = 0;
            }

            //if needed resize the vector for the calculation of reactions
            if (BaseType::mCalculateReactionsFlag == true)
            {
                KRATOS_ERROR(std::logic_error, "calculation of reactions not yet implemented with Trilinos", "");
            }

            KRATOS_CATCH("")
        }

        void InitializeSolutionStep(ModelPart& r_model_part,
                                    TSystemMatrixType& A,
                                    TSystemVectorType& Dx,
                                    TSystemVectorType& b)
        {
            KRATOS_TRY
            mFirstIteration = true;
            KRATOS_CATCH("")
        }

        void FinalizeSolutionStep(ModelPart& r_model_part,
                                  TSystemMatrixType& A,
                                  TSystemVectorType& Dx,
                                  TSystemVectorType& b)
        {}

        void CalculateReactions(typename TSchemeType::Pointer pScheme,
                                ModelPart& r_model_part,
                                TSystemMatrixType& A,
                                TSystemVectorType& Dx,
                                TSystemVectorType& b)
        {
            KRATOS_ERROR(std::logic_error, "method CalculateReactions not implemented in Trilinos Builder And Solver ", "")
        }

        void BuildLHS_CompleteOnFreeRows(typename TSchemeType::Pointer pScheme,
                                         ModelPart& r_model_part,
                                         TSystemMatrixType& A)
        {
            KRATOS_ERROR(std::logic_error, "method BuildLHS_CompleteOnFreeRows not implemented in Trilinos Builder And Solver ", "");
        }

        void ApplyDirichletConditions(typename TSchemeType::Pointer pScheme,
                                      ModelPart& r_model_part,
                                      TSystemMatrixType& A,
                                      TSystemVectorType& Dx,
                                      TSystemVectorType& b)
        {}

        void ApplyPointLoads(typename TSchemeType::Pointer pScheme,
                             ModelPart& r_model_part,
                             TSystemVectorType& b)
        {}

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

            if (mVelocityCorrection == 1 && this->mpIDiagS != NULL)
            {
                // Note that this vector is not cleaned via TSparseSpace because
                // TSystemVectorType can be an Epetra_FEVector, but this has to
                // be a regular Epetra_Vector
                int global_elems = 0;
                Epetra_Map Map(global_elems, 0, mrComm);
                boost::shared_ptr<Epetra_Vector> pNewEmptyX(new Epetra_Vector(Map));
                mpIDiagS.swap(pNewEmptyX);
            }

            if (this->GetEchoLevel() > 0)
            {
                KRATOS_WATCH("ResidualBasedEliminationBuilderAndSolver Clear Function called");
            }
        }

    protected:

    private:

        inline void Assemble(TSystemVectorType& b,
                             LocalSystemMatrixType& LHS_Contribution,
                             LocalSystemVectorType& RHS_Contribution,
                             std::vector<std::size_t>& EquationId)
        {
            unsigned int Size = EquationId.size();

            int * pVGlobalIds(new int[Size]), * pPGlobalIds(new int[Size]);
            int * pVLocalIds(new int[Size]), * pPLocalIds(new int[Size]);
            int VTermNum(0), PTermNum(0);

            int RHSTermNum(0);
            int* pRHSIndices(new int[Size]);
            double* pRHSValues(new double[Size]);

            for (unsigned int i = 0; i < Size; i++)
            {
                if (EquationId[i] < mVelFreeDofs)
                {
                    pVLocalIds[VTermNum] = i;
                    pVGlobalIds[VTermNum++] = EquationId[i];
                    pRHSValues[RHSTermNum] = RHS_Contribution[i];
                    pRHSIndices[RHSTermNum++] = EquationId[i];
                }
                else if (EquationId[i] < BaseType::mEquationSystemSize)
                {
                    pPLocalIds[PTermNum] = i;
                    pPGlobalIds[PTermNum++] = EquationId[i] - mVelFreeDofs;
                    pRHSValues[RHSTermNum] = RHS_Contribution[i];
                    pRHSIndices[RHSTermNum++] = EquationId[i];
                }
            }

            SumMatrixContribution(mpS, VTermNum, pVGlobalIds, pVLocalIds, VTermNum, pVGlobalIds, pVLocalIds, LHS_Contribution);
            SumMatrixContribution(mpG, VTermNum, pVGlobalIds, pVLocalIds, PTermNum, pPGlobalIds, pPLocalIds, LHS_Contribution);
            SumMatrixContribution(mpD, PTermNum, pPGlobalIds, pPLocalIds, VTermNum, pVGlobalIds, pVLocalIds, LHS_Contribution);
            SumMatrixContribution(mpL, PTermNum, pPGlobalIds, pPLocalIds, PTermNum, pPGlobalIds, pPLocalIds, LHS_Contribution);

            b.SumIntoGlobalValues(RHSTermNum, pRHSIndices, pRHSValues);

            delete [] pVGlobalIds;
            delete [] pPGlobalIds;
            delete [] pVLocalIds;
            delete [] pPLocalIds;
            delete [] pRHSIndices;
            delete [] pRHSValues;
        }

        inline void SumMatrixContribution(TSystemMatrixPointerType pB, // pointer to one of the 4 matrix blocks
                                          const int RowNum, // Number of Rows
                                          const int* pGlobalRowIndices, // Global row indices
                                          const int* pElementalRowIndices, // Row indices in this element's contribution
                                          const int ColNum, // Number of columns
                                          const int* pGlobalColIndices, // Global column indices
                                          const int* pElementalColIndices, // Column indices in this element's contribution
                                          const LocalSystemMatrixType& Contribution) // 'Monolithic' elemental matrix
        {
            double pValues[RowNum * ColNum];

            for (int i = 0; i < RowNum; i++)
                for (int j = 0; j < ColNum; j++)
                    pValues[ColNum * i + j] = Contribution(pElementalRowIndices[i], pElementalColIndices[j]);

            pB->SumIntoGlobalValues(RowNum, pGlobalRowIndices, ColNum, pGlobalColIndices, &pValues[0], Epetra_FECrsMatrix::ROW_MAJOR);

            //delete [] pValues;
        }

        inline void AssembleLHS(LocalSystemMatrixType& LHS_Contribution,
                                std::vector<std::size_t>& EquationId)
        {
            unsigned int Size = EquationId.size();

            int * pVGlobalIds(new int[Size]), * pPGlobalIds(new int[Size]);
            int * pVLocalIds(new int[Size]), * pPLocalIds(new int[Size]);
            int VTermNum(0), PTermNum(0);

            for (unsigned int i = 0; i < Size; i++)
            {
                if (EquationId[i] < mVelFreeDofs)
                {
                    pVLocalIds[VTermNum] = i;
                    pVGlobalIds[VTermNum++] = EquationId[i];
                }
                else if (EquationId[i] < BaseType::mEquationSystemSize)
                {
                    pPLocalIds[PTermNum] = i;
                    pPGlobalIds[PTermNum++] = EquationId[i] - mVelFreeDofs;
                }
            }

            SumMatrixContribution(mpS, VTermNum, pVGlobalIds, pVLocalIds, VTermNum, pVGlobalIds, pVLocalIds, LHS_Contribution);
            SumMatrixContribution(mpG, VTermNum, pVGlobalIds, pVLocalIds, PTermNum, pPGlobalIds, pPLocalIds, LHS_Contribution);
            SumMatrixContribution(mpD, PTermNum, pPGlobalIds, pPLocalIds, VTermNum, pVGlobalIds, pVLocalIds, LHS_Contribution);
            SumMatrixContribution(mpL, PTermNum, pPGlobalIds, pPLocalIds, PTermNum, pPGlobalIds, pPLocalIds, LHS_Contribution);

            delete [] pVGlobalIds;
            delete [] pPGlobalIds;
            delete [] pVLocalIds;
            delete [] pPLocalIds;
        }

        inline void AssembleRHS(TSystemVectorType& b,
                                LocalSystemVectorType& RHS_Contribution,
                                std::vector<std::size_t>& EquationId)
        {
            unsigned int Size = EquationId.size();

            int RHSTermNum(0);
            int* pRHSIndices(new int[Size]);
            double* pRHSValues(new double[Size]);

            for (unsigned int i = 0; i < Size; i++)
            {
                if (EquationId[i] < BaseType::mEquationSystemSize)
                {
                    pRHSValues[RHSTermNum] = RHS_Contribution[i];
                    pRHSIndices[RHSTermNum++] = EquationId[i];
                }
            }

            b.SumIntoGlobalValues(RHSTermNum, pRHSIndices, pRHSValues);

            delete [] pRHSIndices;
            delete [] pRHSValues;
        }

        void AssembleLHS_CompleteOnFreeRows(TSystemMatrixType& A,
                                            LocalSystemMatrixType& LHS_Contribution,
                                            Element::EquationIdVectorType& EquationId)
        {
            KRATOS_ERROR(std::logic_error, "AssembleLHS_CompleteOnFreeRows() method is not implemented for Trilinos", "");
        }

        /// Compute the Graph of A from the shapes of D, G and L
        void ConstructSystemMatrixGraph(boost::shared_ptr<Epetra_CrsGraph>& pGraphA)
        {
            KRATOS_TRY// Fill the matrices with ones to compute the result. I'm counting on the next build() to clean them
            mpD->PutScalar(1.0);
            mpG->PutScalar(1.0);
            mpL->PutScalar(1.0);

            TSystemMatrixPointerType pA(new TSystemMatrixType(Copy, mpL->RowMap(), 0));

            // Get the shape of D*Inv(Diag(S))*G
            EpetraExt::MatrixMatrix::Multiply(*mpD, false, *mpG, false, *pA, false);
            // Get the shape of L - D*Inv(Diag(S))*G
            EpetraExt::MatrixMatrix::Add(*mpL, false, 1.0, *pA, 1.0);
            // Finalize the result, providing a DomainMap and a RangeMap. We provide the map because
            // the matrix was obtained from a product of rectangular matices, and it is not clear that
            // it has the proper ones. Fortunately, we know which maps we want for the result.
            pA->GlobalAssemble(mpL->DomainMap(), mpL->RangeMap(), true);

            // Store the graph using given variable
            boost::shared_ptr<Epetra_CrsGraph> pNewGraphA(new Epetra_CrsGraph(pA->Graph()));
            pGraphA.swap(pNewGraphA);
            KRATOS_CATCH("")
        }

        /// Set iterative solver tolerance using inexact Newton criteria
        void SetInitialTolerance(double RHSNorm,
                                 double& TolFactor)
        {
            TolFactor = mMaxTolFactor;
            (BaseType::mpLinearSystemSolver)->SetTolerance(mSmallTol);
            std::cout << "Set iterative solver tolerance to " << TolFactor << std::endl;
        }

        void UpdateTolerance(double OldRHSNorm,
                             double NewRHSNorm,
                             double& TolFactor)
        {
            const double MaxDecreaseFactor = 0.1;

            double CandidateFactor = mGamma * (NewRHSNorm * NewRHSNorm) / (OldRHSNorm * OldRHSNorm);
            std::cout << "Norm Ratio: " << NewRHSNorm / OldRHSNorm << std::endl;
            double CandidateFactor_LimitedDecrease = mGamma * TolFactor*TolFactor;

            if (CandidateFactor_LimitedDecrease < MaxDecreaseFactor)
            {
                TolFactor = (CandidateFactor < mMaxTolFactor) ? CandidateFactor : mMaxTolFactor;
            }
            else
            {
                double Temp = (CandidateFactor > CandidateFactor_LimitedDecrease) ? CandidateFactor : CandidateFactor_LimitedDecrease;
                TolFactor = (Temp < mMaxTolFactor) ? Temp : mMaxTolFactor;
            }

            if (TolFactor < mSmallTol) TolFactor = mSmallTol;
            (BaseType::mpLinearSystemSolver)->SetTolerance(TolFactor);
            std::cout << "Corrected iterative solver tolerance to " << TolFactor << std::endl;
        }

        /// Position of 1st Free Pressure Dof
        unsigned int mVelFreeDofs;
        // Position of 1st Fixed Pressure Dof
        unsigned int mVelFixedDofsEnd;

        /// Number of Pressure Dofs
        unsigned int mPressFreeDofs;

        /// System Matrix for the momentum equation (assumning a given pressure)
        TSystemMatrixPointerType mpS;
        /// Discrete Divergence operator
        TSystemMatrixPointerType mpD;
        /// Discrete Gradient operator
        TSystemMatrixPointerType mpG;
        /// Stabilization term
        TSystemMatrixPointerType mpL;

        /// Inv(Diag(S)), stored as a vector
        boost::shared_ptr<Epetra_Vector> mpIDiagS;

        /// Pointer to the velocity system solver (The pressure one is stored in BaseType::mpLinearSystemSolver)
        typename TLinearSolver::Pointer mpVelocityLinearSystemSolver;

        /// MPI Communicator
        Epetra_MpiComm& mrComm;
        /// Approximate row size
        int mRowSizeGuess;

        int mLocalVelBegin;
        int mLocalPressBegin;
        int mLocalVelDofNum;
        int mLocalPressDofNum;

        /// Flags for matrix reconstruction
        bool mDofSetChanged;

        /** Ensure that velocity is divergence free after each iteration
         * 0: Do not
         * 1: Solve the system Diag(S)*dv = -G*dp
         * 2: Solve the full system S*dv = -G*dp
         */
        unsigned int mVelocityCorrection;

        // Variables for inexact Newton tolerance control
        /// Use Inexact Newton method flag
        bool mInexactNewton;
        /// used to compute maximum solution tolerance
        double mMaxTolFactor;
        /// Minimum solution tolerance: 0.5*NonLinearTol
        double mSmallTol;
        double mGamma;

        /// Keeps track of the first iteration in each step
        bool mFirstIteration;
        double mPressTolFactor;
        double mLastPressRHSNorm;
    };
}

#endif	/* KRATOS_TRILINOS_PRESSURE_SPLITTING_BUILDER_AND_SOLVER_H */
