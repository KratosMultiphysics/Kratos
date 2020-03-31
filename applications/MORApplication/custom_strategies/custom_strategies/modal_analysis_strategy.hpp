//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Author1 Fullname
//                   Author2 Fullname
//

#if !defined(KRATOS_MODAL_ANALYSIS_STRATEGY )
#define  KRATOS_MODAL_ANALYSIS_STRATEGY

// System includes
#include "includes/define.h"

// External includes
#include "includes/model_part.h"

// Project includes
#include "solving_strategies/strategies/solving_strategy.h"
#include "utilities/builtin_timer.h"
#include "spaces/ublas_space.h"

// Application includes
#include "mor_application_variables.h"
#include "custom_utilities/dirichlet_utility.hpp"

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
template<
    class TSparseSpace,
    class TDenseSpace,
    class TLinearSolver, //= LinearSolver<TSparseSpace, TDenseSpace>,
    class TSparseSpaceOut,
    class TDenseSpaceOut
    >
class ModalAnalysisStrategy
    : public SolvingStrategy<TSparseSpace, TDenseSpace, LinearSolver<TSparseSpace, TDenseSpace>>
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(ModalAnalysisStrategy);

    typedef SolvingStrategy<TSparseSpace, TDenseSpace, LinearSolver<TSparseSpace, TDenseSpace>> BaseType;

    typedef typename BaseType::TSchemeType::Pointer SchemePointerType;

    typedef typename BaseType::TBuilderAndSolverType::Pointer BuilderAndSolverPointerType;

    typedef typename TDenseSpaceOut::VectorType DenseVectorType;

    typedef typename TDenseSpaceOut::MatrixType DenseMatrixType;

    typedef typename TDenseSpaceOut::VectorPointerType DenseVectorPointerType;

    typedef typename TDenseSpaceOut::MatrixPointerType DenseMatrixPointerType;

    typedef TSparseSpace SparseSpaceType;

    typedef typename TSparseSpace::VectorPointerType SparseVectorPointerType;

    typedef typename TSparseSpace::MatrixPointerType SparseMatrixPointerType;

    typedef typename TSparseSpace::MatrixType SparseMatrixType;

    typedef typename TSparseSpace::VectorType SparseVectorType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    ModalAnalysisStrategy(
        ModelPart& rModelPart,
        SchemePointerType pScheme,
        BuilderAndSolverPointerType pBuilderAndSolver,
        typename TLinearSolver::Pointer pLinearSolver
        // bool ComputeModalDecomposition = false
        )
        : SolvingStrategy<TSparseSpace, TDenseSpace, LinearSolver<TSparseSpace, TDenseSpace>>(rModelPart),
            mpScheme(pScheme),
            mpBuilderAndSolver(pBuilderAndSolver),
            mpLinearSolver(pLinearSolver)
    {
        KRATOS_TRY

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

        mpEigenvalues = TDenseSpaceOut::CreateEmptyVectorPointer();
        mpEigenvectors = TDenseSpaceOut::CreateEmptyMatrixPointer();

        KRATOS_CATCH("")
    }

    /// Deleted copy constructor.
    ModalAnalysisStrategy(const ModalAnalysisStrategy& Other) = delete;

    /// Destructor.
    ~ModalAnalysisStrategy() override
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
        return *mpMassMatrix;
    }

    SparseMatrixType& GetStiffnessMatrix()
    {
        return *mpStiffnessMatrix;
    }

    SparseMatrixPointerType& pGetMassMatrix()
    {
        return mpMassMatrix;
    }

    SparseMatrixPointerType& pGetStiffnessMatrix()
    {
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

        KRATOS_INFO_IF("ModalAnalysisStrategy", BaseType::GetEchoLevel() > 2 && rank == 0)
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

        KRATOS_INFO_IF("ModalAnalysisStrategy", BaseType::GetEchoLevel() > 2 && rank == 0)
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

        KRATOS_INFO_IF("ModalAnalysisStrategy", BaseType::GetEchoLevel() > 2 && rank == 0)
            <<  "Entering InitializeSolutionStep" << std::endl;

        BuilderAndSolverPointerType& pBuilderAndSolver = this->pGetBuilderAndSolver();
        SchemePointerType& pScheme = this->pGetScheme();
        SparseMatrixPointerType& pStiffnessMatrix = this->pGetStiffnessMatrix();
        SparseMatrixType& rStiffnessMatrix = this->GetStiffnessMatrix();
        SparseMatrixType& rMassMatrix = this->GetMassMatrix();

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

            //set up system matrices
            std::vector<unsigned int> fixed_dofs;
            DirichletUtility::GetDirichletConstraints<typename BaseType::TBuilderAndSolverType::DofsArrayType>(pBuilderAndSolver->GetDofSet(), fixed_dofs);

            // set up mass matrix
            rModelPart.GetProcessInfo()[BUILD_LEVEL] = 201;
            // TSparseSpace::SetToZero(rMassMatrix);
            this->pGetBuilderAndSolver()->Build(pScheme,rModelPart,rMassMatrix,rb);
            DirichletUtility::ApplyDirichletConditions<TSparseSpace>(rMassMatrix, rb, fixed_dofs, 1.0);

            if (BaseType::GetEchoLevel() == 4) {
                TSparseSpace::WriteMatrixMarketMatrix("MassMatrix.mm", rMassMatrix, false);
            }

            // Generate stiffness matrix. the factor -1 is chosen to make
            // Eigenvalues corresponding to fixed dofs negative
            rModelPart.GetProcessInfo()[BUILD_LEVEL] = 1;
            // TSparseSpace::SetToZero(rStiffnessMatrix);
            this->pGetBuilderAndSolver()->Build(pScheme,rModelPart,rStiffnessMatrix,rb);
            DirichletUtility::ApplyDirichletConditions<TSparseSpace>(rStiffnessMatrix, rb, fixed_dofs, -1.0);
            // ApplyDirichletConditions(rStiffnessMatrix,-1.0);

            if (BaseType::GetEchoLevel() == 4) {
                TSparseSpace::WriteMatrixMarketMatrix("StiffnessMatrix.mm", rStiffnessMatrix, false);
            }
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

        KRATOS_INFO_IF("ModalAnalysisStrategy", BaseType::GetEchoLevel() > 2 && rank == 0)
            <<  "Exiting InitializeSolutionStep" << std::endl;

        KRATOS_CATCH("")
    }

    bool SolveSolutionStep() override
    {
        KRATOS_TRY;

        ModelPart& rModelPart = BaseType::GetModelPart();
        const int rank = rModelPart.GetCommunicator().MyPID();

        // SchemePointerType& pScheme = this->pGetScheme();
        SparseMatrixType& rMassMatrix = this->GetMassMatrix();
        SparseMatrixType& rStiffnessMatrix = this->GetStiffnessMatrix();

        // Eigenvector matrix and eigenvalue vector are initialized by the solver
        DenseVectorType& Eigenvalues = this->GetEigenvalueVector();
        DenseMatrixType& Eigenvectors = this->GetEigenvectorMatrix();

        // Solve for eigenvalues and eigenvectors
        BuiltinTimer system_solve_time;
        mpLinearSolver->Solve(
                rStiffnessMatrix,
                rMassMatrix,
                Eigenvalues,
                Eigenvectors);

        KRATOS_INFO_IF("System Solve Time", BaseType::GetEchoLevel() > 0 && rank == 0)
                << system_solve_time.ElapsedSeconds() << std::endl;


        this->AssignVariables<typename DenseMatrixType::value_type>(Eigenvalues,Eigenvectors);

        return true;
        KRATOS_CATCH("")
    }

    void FinalizeSolutionStep() override
    {
        KRATOS_TRY;

        const int rank = BaseType::GetModelPart().GetCommunicator().MyPID();
        KRATOS_INFO_IF("ModalAnalysisStrategy", BaseType::GetEchoLevel() > 2 && rank == 0)
            <<  "Entering FinalizeSolutionStep" << std::endl;

        SparseMatrixType& rStiffnessMatrix = this->GetStiffnessMatrix();
        SparseVectorPointerType pDx = SparseSpaceType::CreateEmptyVectorPointer();
        SparseVectorPointerType pb = SparseSpaceType::CreateEmptyVectorPointer();
        pGetBuilderAndSolver()->FinalizeSolutionStep(
            BaseType::GetModelPart(), rStiffnessMatrix, *pDx, *pb);
        pGetScheme()->FinalizeSolutionStep(BaseType::GetModelPart(),
                                           rStiffnessMatrix, *pDx, *pb);
        KRATOS_INFO_IF("ModalAnalysisStrategy", BaseType::GetEchoLevel() > 2 && rank == 0)
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

        KRATOS_INFO_IF("ModalAnalysisStrategy", BaseType::GetEchoLevel() > 2 && rank == 0)
            <<  "Entering Check" << std::endl;

        // check the model part
        BaseType::Check();

        // check the scheme
        this->pGetScheme()->Check(rModelPart);

        // check the builder and solver
        this->pGetBuilderAndSolver()->Check(rModelPart);

        KRATOS_INFO_IF("ModalAnalysisStrategy", BaseType::GetEchoLevel() > 2 && rank == 0)
            <<  "Exiting Check" << std::endl;

        return 0;

        KRATOS_CATCH("")
    }

    ///@}
    ///@name Access
    ///@{

    virtual DenseMatrixType& GetEigenvectorMatrix()
    {
        DenseMatrixType& eigenvectors = *mpEigenvectors;

        return eigenvectors;
    }

    virtual DenseVectorType& GetEigenvalueVector()
    {
        DenseVectorType& eigenvalues = *mpEigenvalues;

        return eigenvalues;
    }

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

    typename TLinearSolver::Pointer mpLinearSolver;

    SparseMatrixPointerType mpMassMatrix;

    SparseMatrixPointerType mpStiffnessMatrix;

    DenseVectorPointerType mpEigenvalues;

    DenseMatrixPointerType mpEigenvectors;

    bool mInitializeWasPerformed = false;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    // Assign eigenvalues and eigenvectors to kratos variables.
    template<typename TScalar, typename std::enable_if<std::is_same<double, TScalar>::value, int>::type = 0>
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

            // TO BE VERIFIED!! In the current implmentation of Dofs there are nor reordered and only pushec back.
            // // the jth column index of EIGENVECTOR_MATRIX corresponds to the jth nodal dof. therefore,
            // // the dof ordering must not change.
            // if (NodeDofs.IsSorted() == false)
            // {
            //     NodeDofs.Sort();
            // }

            // fill the EIGENVECTOR_MATRIX
            for (std::size_t i = 0; i < NumEigenvalues; i++) {
                for (std::size_t j = 0; j < NumNodeDofs; j++) {
                    auto itDof = std::begin(NodeDofs) + j;
                    rNodeEigenvectors(i,j) = rEigenvectors(i,(*itDof)->EquationId());
                }
            }
        }
    }

    // Assign eigenvalues and eigenvectors to kratos variables.
    template<typename TScalar, typename std::enable_if<std::is_same<std::complex<double>, TScalar>::value, int>::type = 0>
    void AssignVariables(DenseVectorType& rEigenvalues, DenseMatrixType& rEigenvectors)
    {
        ModelPart& rModelPart = BaseType::GetModelPart();
        const std::size_t NumEigenvalues = rEigenvalues.size();

        Vector ev_r(NumEigenvalues);
        Vector ev_i(NumEigenvalues);
        for( size_t i=0; i<NumEigenvalues; ++i ) {
            ev_r[i] = std::real(rEigenvalues[i]);
            ev_i[i] = std::imag(rEigenvalues[i]);
        }

        // store eigenvalues in process info
        rModelPart.GetProcessInfo()[EIGENVALUE_VECTOR] = ev_r;
        rModelPart.GetProcessInfo()[IMAG_EIGENVALUE_VECTOR] = ev_i;

        for (ModelPart::NodeIterator itNode = rModelPart.NodesBegin(); itNode!= rModelPart.NodesEnd(); itNode++)
        {
            ModelPart::NodeType::DofsContainerType& NodeDofs = itNode->GetDofs();
            const std::size_t NumNodeDofs = NodeDofs.size();
            Matrix& rNodeEigenvectorsReal = itNode->GetValue(EIGENVECTOR_MATRIX);
            Matrix& rNodeEigenvectorsImag = itNode->GetValue(IMAG_EIGENVECTOR_MATRIX);
            if (rNodeEigenvectorsReal.size1() != NumEigenvalues || rNodeEigenvectorsReal.size2() != NumNodeDofs)
            {
                rNodeEigenvectorsReal.resize(NumEigenvalues,NumNodeDofs,false);
                rNodeEigenvectorsImag.resize(NumEigenvalues,NumNodeDofs,false);
            }

            // TO BE VERIFIED!! In the current implmentation of Dofs there are nor reordered and only pushec back.
            // // the jth column index of EIGENVECTOR_MATRIX corresponds to the jth nodal dof. therefore,
            // // the dof ordering must not change.
            // if (NodeDofs.IsSorted() == false)
            // {
            //     NodeDofs.Sort();
            // }

            // fill the EIGENVECTOR_MATRIX
            for (std::size_t i = 0; i < NumEigenvalues; i++) {
                for (std::size_t j = 0; j < NumNodeDofs; j++) {
                    auto itDof = std::begin(NodeDofs) + j;
                    rNodeEigenvectorsReal(i,j) = std::real(rEigenvectors(i,(*itDof)->EquationId()));
                    rNodeEigenvectorsImag(i,j) = std::imag(rEigenvectors(i,(*itDof)->EquationId()));
                }
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

}; /* Class ModalAnalysisStrategy */

///@}

///@name Type Definitions
///@{


///@}

} /* namespace Kratos */

#endif /* KRATOS_MODAL_ANALYSIS_STRATEGY  defined */

