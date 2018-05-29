//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//   License:        BSD License
//   Kratos default license: kratos/license.txt
//
//   Project Name:        $StructuralMechanicsApplication $
//   Last modified by:    $Author: michael.andre@tum.de   $
//   Date:                $Date:         September 2016   $
//   Revision:            $Revision:                0.0   $

#if !defined(KRATOS_EIGENSOLVER_STRATEGY )
#define  KRATOS_EIGENSOLVER_STRATEGY

// System includes

// External includes

// Project includes
#include "solving_strategies/strategies/solving_strategy.h"
#include "utilities/builtin_timer.h"

// Application includes
#include "structural_mechanics_application_variables.h"

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Strategy for solving generalized eigenvalue problems.
template<class TSparseSpace,
         class TDenseSpace,
         class TLinearSolver
         >
class EigensolverStrategy
    : public SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(EigensolverStrategy);

    typedef SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

    typedef typename BaseType::TSchemeType::Pointer SchemePointerType;

    typedef typename BaseType::TBuilderAndSolverType::Pointer BuilderAndSolverPointerType;

    typedef typename TDenseSpace::VectorType DenseVectorType;

    typedef typename TDenseSpace::MatrixType DenseMatrixType;

    typedef TSparseSpace SparseSpaceType;

    typedef typename TSparseSpace::VectorPointerType SparseVectorPointerType;

    typedef typename TSparseSpace::MatrixPointerType SparseMatrixPointerType;

    typedef typename TSparseSpace::MatrixType SparseMatrixType;

    typedef typename TSparseSpace::VectorType SparseVectorType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    EigensolverStrategy(
        ModelPart& rModelPart,
        SchemePointerType pScheme,
        BuilderAndSolverPointerType pBuilderAndSolver
        )
        : SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(rModelPart)
    {
        KRATOS_TRY

        mpScheme = pScheme;

        mpBuilderAndSolver = pBuilderAndSolver;

        // ensure initialization of system matrices in InitializeSolutionStep()
        mpBuilderAndSolver->SetDofSetIsInitializedFlag(false);

        // default echo level (mute)
        this->SetEchoLevel(0);

        // default rebuild level (build at each solution step)
        this->SetRebuildLevel(1);

        SparseMatrixType* AuxMassMatrix = new SparseMatrixType;
        mpMassMatrix = Kratos::shared_ptr<SparseMatrixType>(AuxMassMatrix);
        SparseMatrixType* AuxStiffnessMatrix = new SparseMatrixType;
        mpStiffnessMatrix = Kratos::shared_ptr<SparseMatrixType>(AuxStiffnessMatrix);

        KRATOS_CATCH("")
    }

    /// Deleted copy constructor.
    EigensolverStrategy(const EigensolverStrategy& Other) = delete;

    /// Destructor.
    ~EigensolverStrategy() override
    {
        // Clear() controls order of deallocation to avoid invalid memory access
        // in some special cases.
        // warning: BaseType::GetModelPart() may be invalid here.
        this->Clear();
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    void SetIsInitialized(bool val)
    {
        mInitializeWasPerformed = val;
    }

    bool GetIsInitialized() const
    {
        return mInitializeWasPerformed;
    }

    void SetScheme(SchemePointerType pScheme)
    {
        mpScheme = pScheme;
    };

    SchemePointerType& pGetScheme()
    {
        return mpScheme;
    };

    void SetBuilderAndSolver(BuilderAndSolverPointerType pNewBuilderAndSolver)
    {
        mpBuilderAndSolver = pNewBuilderAndSolver;
    };

    BuilderAndSolverPointerType& pGetBuilderAndSolver()
    {
        return mpBuilderAndSolver;
    };

    SparseMatrixType& GetMassMatrix()
    {
//         if (mpMassMatrix == nullptr)
//         {
//             KRATOS_ERROR << "CHEK YOUR MASS MATRIX PLEASE IS INITIALIZED" << std::endl;
//         }
        return *mpMassMatrix;
    }

    SparseMatrixType& GetStiffnessMatrix()
    {
//         if (mpStiffnessMatrix == nullptr)
//         {
//             KRATOS_ERROR << "CHEK YOUR STIFFNESS MATRIX PLEASE IS INITIALIZED" << std::endl;
//         }
        return *mpStiffnessMatrix;
    }

    SparseMatrixPointerType& pGetMassMatrix()
    {
//         if (mpMassMatrix == nullptr)
//         {
//             KRATOS_ERROR << "CHEK YOUR MASS MATRIX PLEASE IS INITIALIZED" << std::endl;
//         }
        return mpMassMatrix;
    }

    SparseMatrixPointerType& pGetStiffnessMatrix()
    {
//         if (mpStiffnessMatrix == nullptr)
//         {
//             KRATOS_ERROR << "CHEK YOUR STIFFNESS MATRIX PLEASE IS INITIALIZED" << std::endl;
//         }
        return mpStiffnessMatrix;
    }

    void SetReformDofSetAtEachStepFlag(bool flag)
    {
        this->pGetBuilderAndSolver()->SetReshapeMatrixFlag(flag);
    }

    bool GetReformDofSetAtEachStepFlag() const
    {
        return this->pGetBuilderAndSolver()->GetReshapeMatrixFlag();
    }

    /// Set verbosity level of the solving strategy.
    /**
     * - 0 -> mute... no echo at all
     * - 1 -> print time and basic information
     * - 2 -> print linear solver data
     * - 3 -> print debug information
     */
    void SetEchoLevel(int Level) override
    {
        BaseType::SetEchoLevel(Level);
        this->pGetBuilderAndSolver()->SetEchoLevel(Level);
    }

    /**
     * Initialization to be performed once before using the strategy.
     */
    void Initialize() override
    {
        KRATOS_TRY

        ModelPart& rModelPart = BaseType::GetModelPart();
        const int rank = rModelPart.GetCommunicator().MyPID();

        KRATOS_INFO_IF("EigensolverStrategy", BaseType::GetEchoLevel() > 2 && rank == 0)
            <<  "Entering Initialize" << std::endl;

        if (mInitializeWasPerformed == false)
        {
            SchemePointerType& pScheme = this->pGetScheme();

            if (pScheme->SchemeIsInitialized() == false)
                pScheme->Initialize(rModelPart);

            if (pScheme->ElementsAreInitialized() == false)
                pScheme->InitializeElements(rModelPart);

            if (pScheme->ConditionsAreInitialized() == false)
                pScheme->InitializeConditions(rModelPart);
        }

        KRATOS_INFO_IF("EigensolverStrategy", BaseType::GetEchoLevel() > 2 && rank == 0)
            <<  "Exiting Initialize" << std::endl;

        KRATOS_CATCH("")
    }

    /**
     * Clears the internal storage
     */
    void Clear() override
    {
        KRATOS_TRY

        // if the preconditioner is saved between solves, it should be cleared here
        BuilderAndSolverPointerType& pBuilderAndSolver = this->pGetBuilderAndSolver();
        pBuilderAndSolver->GetLinearSystemSolver()->Clear();

        if (this->pGetMassMatrix() != nullptr)
            this->pGetMassMatrix() = nullptr;

        if (this->pGetStiffnessMatrix() != nullptr)
            this->pGetStiffnessMatrix() = nullptr;

        // Re-setting internal flag to ensure that the dof sets are recalculated
        pBuilderAndSolver->SetDofSetIsInitializedFlag(false);

        pBuilderAndSolver->Clear();

        this->pGetScheme()->Clear();

        mInitializeWasPerformed = false;

        KRATOS_CATCH("")
    }

    /**
     * Performs all the required operations that should be done (for each step)
     * before solving the solution step.
     * A member variable should be used as a flag to make sure this function is called only once per step.
     */
    void InitializeSolutionStep() override
    {
        KRATOS_TRY

        ModelPart& rModelPart = BaseType::GetModelPart();
        const int rank = rModelPart.GetCommunicator().MyPID();

        KRATOS_INFO_IF("EigensolverStrategy", BaseType::GetEchoLevel() > 2 && rank == 0)
            <<  "Entering InitializeSolutionStep" << std::endl;

        BuilderAndSolverPointerType& pBuilderAndSolver = this->pGetBuilderAndSolver();
        SchemePointerType& pScheme = this->pGetScheme();
        SparseMatrixPointerType& pStiffnessMatrix = this->pGetStiffnessMatrix();
        SparseMatrixType& rStiffnessMatrix = this->GetStiffnessMatrix();

        // Initialize dummy vectors
        SparseVectorPointerType pDx = SparseSpaceType::CreateEmptyVectorPointer();
        SparseVectorPointerType pb = SparseSpaceType::CreateEmptyVectorPointer();
        auto& rDx = *pDx;
        auto& rb = *pb;

        // Reset solution dofs
        BuiltinTimer system_construction_time;
        if (pBuilderAndSolver->GetDofSetIsInitializedFlag() == false ||
            pBuilderAndSolver->GetReshapeMatrixFlag() == true)
        {
            // Set up list of dofs
            BuiltinTimer setup_dofs_time;
            pBuilderAndSolver->SetUpDofSet(pScheme, rModelPart);

            KRATOS_INFO_IF("Setup Dofs Time", BaseType::GetEchoLevel() > 0 && rank == 0)
                << setup_dofs_time.ElapsedSeconds() << std::endl;

            // Set global equation ids
            BuiltinTimer setup_system_time;
            pBuilderAndSolver->SetUpSystem(rModelPart);

            KRATOS_INFO_IF("Setup System Time", BaseType::GetEchoLevel() > 0 && rank == 0)
                << setup_system_time.ElapsedSeconds() << std::endl;

            // Resize and initialize system matrices
            BuiltinTimer system_matrix_resize_time;
            SparseMatrixPointerType& pMassMatrix = this->pGetMassMatrix();

            // Mass matrix
            pBuilderAndSolver->ResizeAndInitializeVectors(
                pScheme, pMassMatrix, pDx, pb, rModelPart);

            // Stiffness matrix
            pBuilderAndSolver->ResizeAndInitializeVectors(
                pScheme, pStiffnessMatrix, pDx, pb, rModelPart);

            KRATOS_INFO_IF("System Matrix Resize Time", BaseType::GetEchoLevel() > 0 && rank == 0)
                << system_matrix_resize_time.ElapsedSeconds() << std::endl;
        }
        else
        {
            SparseSpaceType::Resize(rb, SparseSpaceType::Size1(rStiffnessMatrix));
            SparseSpaceType::Set(rb, 0.0);
            SparseSpaceType::Resize(rDx, SparseSpaceType::Size1(rStiffnessMatrix));
            SparseSpaceType::Set(rDx, 0.0);
        }

        KRATOS_INFO_IF("System Construction Time", BaseType::GetEchoLevel() > 0 && rank == 0)
            << system_construction_time.ElapsedSeconds() << std::endl;

        // Initial operations ... things that are constant over the solution
        // step
        pBuilderAndSolver->InitializeSolutionStep(BaseType::GetModelPart(),
                                                  rStiffnessMatrix, rDx, rb);

        // Initial operations ... things that are constant over the solution
        // step
        pScheme->InitializeSolutionStep(BaseType::GetModelPart(), rStiffnessMatrix, rDx, rb);

        KRATOS_INFO_IF("EigensolverStrategy", BaseType::GetEchoLevel() > 2 && rank == 0)
            <<  "Exiting InitializeSolutionStep" << std::endl;

        KRATOS_CATCH("")
    }

    bool SolveSolutionStep() override
    {
        KRATOS_TRY;

        ModelPart& rModelPart = BaseType::GetModelPart();
        const int rank = rModelPart.GetCommunicator().MyPID();

        SchemePointerType& pScheme = this->pGetScheme();
        SparseMatrixType& rMassMatrix = this->GetMassMatrix();
        SparseMatrixType& rStiffnessMatrix = this->GetStiffnessMatrix();

        // Initialize dummy rhs vector
        SparseVectorType b;
        SparseSpaceType::Resize(b,SparseSpaceType::Size1(rMassMatrix));
        SparseSpaceType::Set(b,0.0);

        // Generate lhs matrix. the factor 1 is chosen to preserve
        // SPD property
        rModelPart.GetProcessInfo()[BUILD_LEVEL] = 1;
        TSparseSpace::SetToZero(rMassMatrix);
        this->pGetBuilderAndSolver()->Build(pScheme,rModelPart,rMassMatrix,b);
        this->ApplyDirichletConditions(rMassMatrix, 1.0);

        // Generate rhs matrix. the factor -1 is chosen to make
        // Eigenvalues corresponding to fixed dofs negative
        rModelPart.GetProcessInfo()[BUILD_LEVEL] = 2;
        TSparseSpace::SetToZero(rStiffnessMatrix);
        this->pGetBuilderAndSolver()->Build(pScheme,rModelPart,rStiffnessMatrix,b);
        ApplyDirichletConditions(rStiffnessMatrix,-1.0);

        // Eigenvector matrix and eigenvalue vector are initialized by the solver
        DenseVectorType Eigenvalues;
        DenseMatrixType Eigenvectors;

        // Solve for eigenvalues and eigenvectors
        BuiltinTimer system_solve_time;
        this->pGetBuilderAndSolver()->GetLinearSystemSolver()->Solve(
                rStiffnessMatrix,
                rMassMatrix,
                Eigenvalues,
                Eigenvectors);

        KRATOS_INFO_IF("System Solve Time", BaseType::GetEchoLevel() > 0 && rank == 0)
                << system_solve_time.ElapsedSeconds() << std::endl;

        this->AssignVariables(Eigenvalues,Eigenvectors);

        return true;
        KRATOS_CATCH("")
    }

    void FinalizeSolutionStep() override
    {
        KRATOS_TRY;

        const int rank = BaseType::GetModelPart().GetCommunicator().MyPID();
        KRATOS_INFO_IF("EigensolverStrategy", BaseType::GetEchoLevel() > 2 && rank == 0)
            <<  "Entering FinalizeSolutionStep" << std::endl;

        SparseMatrixType& rStiffnessMatrix = this->GetStiffnessMatrix();
        SparseVectorPointerType pDx = SparseSpaceType::CreateEmptyVectorPointer();
        SparseVectorPointerType pb = SparseSpaceType::CreateEmptyVectorPointer();
        pGetBuilderAndSolver()->FinalizeSolutionStep(
            BaseType::GetModelPart(), rStiffnessMatrix, *pDx, *pb);
        pGetScheme()->FinalizeSolutionStep(BaseType::GetModelPart(),
                                           rStiffnessMatrix, *pDx, *pb);
        KRATOS_INFO_IF("EigensolverStrategy", BaseType::GetEchoLevel() > 2 && rank == 0)
            <<  "Exiting FinalizeSolutionStep" << std::endl;

        KRATOS_CATCH("");
    }

    /**
     * Function to perform expensive checks.
     * It is designed to be called ONCE to verify that the input is correct.
     */
    int Check() override
    {
        KRATOS_TRY

        ModelPart& rModelPart = BaseType::GetModelPart();
        const int rank = rModelPart.GetCommunicator().MyPID();

        KRATOS_INFO_IF("EigensolverStrategy", BaseType::GetEchoLevel() > 2 && rank == 0)
            <<  "Entering Check" << std::endl;

        // check the model part
        BaseType::Check();

        // check the scheme
        this->pGetScheme()->Check(rModelPart);

        // check the builder and solver
        this->pGetBuilderAndSolver()->Check(rModelPart);

        KRATOS_INFO_IF("EigensolverStrategy", BaseType::GetEchoLevel() > 2 && rank == 0)
            <<  "Exiting Check" << std::endl;

        return 0;

        KRATOS_CATCH("")
    }

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    SchemePointerType mpScheme;

    BuilderAndSolverPointerType mpBuilderAndSolver;

    SparseMatrixPointerType mpMassMatrix;

    SparseMatrixPointerType mpStiffnessMatrix;

    bool mInitializeWasPerformed = false;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /// Apply Dirichlet boundary conditions without modifying dof pattern.
    /**
     *  The dof pattern is preserved to support algebraic multigrid solvers with
     *  component-wise aggregation. Rows and columns of the fixed dofs are replaced
     *  with zeros on the off-diagonal and the diagonal is scaled by factor.
     */
    void ApplyDirichletConditions(
        SparseMatrixType& rA,
        double Factor)
    {
        KRATOS_TRY

        const int rank = BaseType::GetModelPart().GetCommunicator().MyPID();

        KRATOS_INFO_IF("EigensolverStrategy", BaseType::GetEchoLevel() > 2 && rank == 0)
            <<  "Entering ApplyDirichletConditions" << std::endl;

        const std::size_t SystemSize = rA.size1();
        std::vector<double> ScalingFactors(SystemSize);
        auto& rDofSet = this->pGetBuilderAndSolver()->GetDofSet();
        const int NumDofs = static_cast<int>(rDofSet.size());

        // NOTE: dofs are assumed to be numbered consecutively
        #pragma omp parallel for firstprivate(NumDofs)
        for(int k = 0; k<NumDofs; k++)
        {
            auto dof_iterator = std::begin(rDofSet) + k;
            ScalingFactors[k] = (dof_iterator->IsFixed()) ? 0.0 : 1.0;
        }

        double* AValues = std::begin(rA.value_data());
        std::size_t* ARowIndices = std::begin(rA.index1_data());
        std::size_t* AColIndices = std::begin(rA.index2_data());

        // if there is a line of all zeros, put one on the diagonal
        // #pragma omp parallel for firstprivate(SystemSize)
        // for(int k = 0; k < static_cast<int>(SystemSize); ++k)
        // {
        //     std::size_t ColBegin = ARowIndices[k];
        //     std::size_t ColEnd = ARowIndices[k+1];
        //     bool empty = true;
        //     for (auto j = ColBegin; j < ColEnd; ++j)
        //         if(AValues[j] != 0.0)
        //         {
        //             empty = false;
        //             break;
        //         }
        //     if(empty == true)
        //         rA(k,k) = 1.0;
        // }

        #pragma omp parallel for
        for (int k = 0; k < static_cast<int>(SystemSize); ++k)
        {
            std::size_t ColBegin = ARowIndices[k];
            std::size_t ColEnd = ARowIndices[k+1];
            if (ScalingFactors[k] == 0.0)
            {
                // row dof is fixed. zero off-diagonal columns and factor diagonal
                for (std::size_t j = ColBegin; j < ColEnd; ++j)
                {
                    if (static_cast<int>(AColIndices[j]) != k)
                    {
                        AValues[j] = 0.0;
                    }
                    else
                    {
                        AValues[j] *= Factor;
                    }
                }
            }
            else
            {
                // row dof is not fixed. zero columns associated with fixed dofs
                for (std::size_t j = ColBegin; j < ColEnd; ++j)
                {
                    AValues[j] *= ScalingFactors[AColIndices[j]];
                }
            }
        }

        KRATOS_INFO_IF("EigensolverStrategy", BaseType::GetEchoLevel() > 2 && rank == 0)
            <<  "Exiting ApplyDirichletConditions" << std::endl;

        KRATOS_CATCH("")
    }

    /// Assign eigenvalues and eigenvectors to kratos variables.
    void AssignVariables(DenseVectorType& rEigenvalues, DenseMatrixType& rEigenvectors)
    {
        ModelPart& rModelPart = BaseType::GetModelPart();
        const std::size_t NumEigenvalues = rEigenvalues.size();

        // store eigenvalues in process info
        rModelPart.GetProcessInfo()[EIGENVALUE_VECTOR] = rEigenvalues;

        for (ModelPart::NodeIterator itNode = rModelPart.NodesBegin(); itNode!= rModelPart.NodesEnd(); itNode++)
        {
            ModelPart::NodeType::DofsContainerType& NodeDofs = itNode->GetDofs();
            const std::size_t NumNodeDofs = NodeDofs.size();
            Matrix& rNodeEigenvectors = itNode->GetValue(EIGENVECTOR_MATRIX);
            if (rNodeEigenvectors.size1() != NumEigenvalues || rNodeEigenvectors.size2() != NumNodeDofs)
            {
                rNodeEigenvectors.resize(NumEigenvalues,NumNodeDofs,false);
            }

            // the jth column index of EIGENVECTOR_MATRIX corresponds to the jth nodal dof. therefore,
            // the dof ordering must not change.
            if (NodeDofs.IsSorted() == false)
            {
                NodeDofs.Sort();
            }

            // fill the EIGENVECTOR_MATRIX
            for (std::size_t i = 0; i < NumEigenvalues; i++)
                for (std::size_t j = 0; j < NumNodeDofs; j++)
                {
                    auto itDof = std::begin(NodeDofs) + j;
                    rNodeEigenvectors(i,j) = rEigenvectors(i,itDof->EquationId());
                }
        }
    }

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}

}; /* Class EigensolverStrategy */

///@}

///@name Type Definitions
///@{


///@}

} /* namespace Kratos */

#endif /* KRATOS_EIGENSOLVER_STRATEGY  defined */

