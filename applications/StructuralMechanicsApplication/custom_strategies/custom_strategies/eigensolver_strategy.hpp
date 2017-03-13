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
#include<iostream>
#include<vector>
#include<iterator>

// External includes
#include<boost/timer.hpp>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/ublas_interface.h"
#include "solving_strategies/strategies/solving_strategy.h"

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

    typedef typename TDenseSpace::VectorPointerType DenseVectorPointerType;

    typedef typename TDenseSpace::MatrixPointerType DenseMatrixPointerType;

    typedef typename TDenseSpace::VectorType DenseVectorType;

    typedef typename TDenseSpace::MatrixType DenseMatrixType;

    typedef TSparseSpace SparseSpaceType;

    typedef typename TSparseSpace::MatrixPointerType SparseMatrixPointerType;

    typedef typename TSparseSpace::MatrixType SparseMatrixType;

    typedef typename TSparseSpace::VectorType SparseVectorType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    EigensolverStrategy(
        ModelPart& model_part,
        SchemePointerType pScheme,
        BuilderAndSolverPointerType pBuilderAndSolver
        )
        : SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(model_part)
    {
        KRATOS_TRY

        mpScheme = pScheme;

        mpBuilderAndSolver = pBuilderAndSolver;

        // ensure initialization of system matrices in InitializeSolutionStep()
        mpBuilderAndSolver->SetDofSetIsInitializedFlag(false);

        mInitializeWasPerformed = false;

        // default echo level (mute)
        this->SetEchoLevel(0);

        // default rebuild level (build at each solution step)
        this->SetRebuildLevel(1);

        KRATOS_CATCH("")
    }

    /// Deleted copy constructor.
    EigensolverStrategy(const EigensolverStrategy& Other) = delete;

    /// Destructor.
    virtual ~EigensolverStrategy()
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
    void SetEchoLevel(int Level)
    {
        BaseType::SetEchoLevel(Level);
        this->pGetBuilderAndSolver()->SetEchoLevel(Level);
    }

    /// Initialization to be performed once before using the strategy.
    virtual void Initialize()
    {
        KRATOS_TRY

        auto& rModelPart = BaseType::GetModelPart();
        const auto rank = rModelPart.GetCommunicator().MyPID();

        if (BaseType::GetEchoLevel() > 2 && rank == 0)
            std::cout << "Entering Initialize() of EigensolverStrategy." << std::endl;

        this->Check();

        auto& pScheme = this->pGetScheme();

        if (pScheme->SchemeIsInitialized() == false)
            pScheme->Initialize(rModelPart);

        if (pScheme->ElementsAreInitialized() == false)
            pScheme->InitializeElements(rModelPart);

        if (pScheme->ConditionsAreInitialized() == false)
            pScheme->InitializeConditions(rModelPart);

        if (BaseType::GetEchoLevel() > 2 && rank == 0)
            std::cout << "Exiting Initialize() of EigensolverStrategy." << std::endl;

        KRATOS_CATCH("")
    }

    double Solve()
    {
        KRATOS_TRY

        auto& rModelPart = BaseType::GetModelPart();
        const auto rank = rModelPart.GetCommunicator().MyPID();

        // operations to be done once
        if (this->GetIsInitialized() == false)
        {
            Initialize();
            this->SetIsInitialized(true);
        }

        this->InitializeSolutionStep();

        auto& pScheme = this->pGetScheme();
        auto& rMassMatrix = this->GetMassMatrix();
        auto& rStiffnessMatrix = this->GetStiffnessMatrix();

        // initialize dummy rhs vector
        SparseVectorType b;
        SparseSpaceType::Resize(b,SparseSpaceType::Size1(rMassMatrix));
        SparseSpaceType::Set(b,0.0);

        // generate lhs matrix. the factor 1 is chosen to preserve
        // spd property
        rModelPart.GetProcessInfo()[BUILD_LEVEL] = 1;
        this->pGetBuilderAndSolver()->Build(pScheme,rModelPart,rMassMatrix,b);
        this->ApplyDirichletConditions(rMassMatrix, 1.0);

        // generate rhs matrix. the factor -1 is chosen to make
        // eigenvalues corresponding to fixed dofs negative
        rModelPart.GetProcessInfo()[BUILD_LEVEL] = 2;
        this->pGetBuilderAndSolver()->Build(pScheme,rModelPart,rStiffnessMatrix,b);
        ApplyDirichletConditions(rStiffnessMatrix,-1.0);

        // eigenvector matrix and eigenvalue vector are initialized by the solver
        DenseVectorType Eigenvalues;
        DenseMatrixType Eigenvectors;

        // solve for eigenvalues and eigenvectors
        boost::timer system_solve_time;
        this->pGetBuilderAndSolver()->GetLinearSystemSolver()->Solve(
                rStiffnessMatrix,
                rMassMatrix,
                Eigenvalues,
                Eigenvectors);
        if (BaseType::GetEchoLevel() > 0 && rank == 0)
            std::cout << "system_solve_time : " << system_solve_time.elapsed() << std::endl;

        this->AssignVariables(Eigenvalues,Eigenvectors);

        this->FinalizeSolutionStep();

        return 0.0;

        KRATOS_CATCH("")
    }

    /// Clear the strategy.
    virtual void Clear()
    {
        KRATOS_TRY

        // if the preconditioner is saved between solves, it should be cleared here
        auto& pBuilderAndSolver = this->pGetBuilderAndSolver();
        pBuilderAndSolver->GetLinearSystemSolver()->Clear();

        if (this->pGetMassMatrix() != nullptr)
            this->pGetMassMatrix() = nullptr;

        if (this->pGetStiffnessMatrix() != nullptr)
            this->pGetStiffnessMatrix() = nullptr;


        // re-setting internal flag to ensure that the dof sets are recalculated
        pBuilderAndSolver->SetDofSetIsInitializedFlag(false);

        pBuilderAndSolver->Clear();

        this->pGetScheme()->Clear();

        mInitializeWasPerformed = false;

        KRATOS_CATCH("")
    }

    /// Initialization to be performed before every solve.
    virtual void InitializeSolutionStep()
    {
        KRATOS_TRY

        auto& rModelPart = BaseType::GetModelPart();
        const auto rank = rModelPart.GetCommunicator().MyPID();

        if (BaseType::GetEchoLevel() > 2 && rank == 0)
            std::cout << "Entering InitializeSolutionStep() of EigensolverStrategy" << std::endl;

        auto& pBuilderAndSolver = this->pGetBuilderAndSolver();
        auto& pScheme = this->pGetScheme();
        auto& pStiffnessMatrix = this->pGetStiffnessMatrix();
        auto& rStiffnessMatrix = *pStiffnessMatrix;

        // initialize dummy vectors
        auto pDx = SparseSpaceType::CreateEmptyVectorPointer();
        auto pb = SparseSpaceType::CreateEmptyVectorPointer();
        auto& rDx = *pDx;
        auto& rb = *pb;

        // reset solution dofs
        boost::timer system_construction_time;
        if (pBuilderAndSolver->GetDofSetIsInitializedFlag() == false ||
                pBuilderAndSolver->GetReshapeMatrixFlag() == true)
        {
            // set up list of dofs
            boost::timer setup_dofs_time;
            pBuilderAndSolver->SetUpDofSet(pScheme, rModelPart);
            if (BaseType::GetEchoLevel() > 0 && rank == 0)
                std::cout << "setup_dofs_time : " << setup_dofs_time.elapsed() << std::endl;

            // set global equation ids
            boost::timer setup_system_time;
            pBuilderAndSolver->SetUpSystem(rModelPart);
            if (BaseType::GetEchoLevel() > 0 && rank == 0)
                std::cout << "setup_system_time : " << setup_system_time.elapsed() << std::endl;

            // resize and initialize system matrices
            boost::timer system_matrix_resize_time;
            auto& pMassMatrix = this->pGetMassMatrix();

            // mass matrix
            pBuilderAndSolver->ResizeAndInitializeVectors(
                    pMassMatrix,
                    pDx,
                    pb,
                    rModelPart.Elements(),
                    rModelPart.Conditions(),
                    rModelPart.GetProcessInfo());

            // stiffness matrix
            pBuilderAndSolver->ResizeAndInitializeVectors(
                    pStiffnessMatrix,
                    pDx,
                    pb,
                    rModelPart.Elements(),
                    rModelPart.Conditions(),
                    rModelPart.GetProcessInfo());

            if (BaseType::GetEchoLevel() > 0 && rank == 0)
                std::cout << "system_matrix_resize_time : " << system_matrix_resize_time.elapsed() << std::endl;
        }
        else
        {
            SparseSpaceType::Resize(rb,SparseSpaceType::Size1(rStiffnessMatrix));
            SparseSpaceType::Set(rb,0.0);
            SparseSpaceType::Resize(rDx,SparseSpaceType::Size1(rStiffnessMatrix));
            SparseSpaceType::Set(rDx,0.0);
        }
        if (BaseType::GetEchoLevel() > 0 && rank == 0)
            std::cout << "system_construction_time : " << system_construction_time.elapsed() << std::endl;

        // initial operations ... things that are constant over the solution step
        pBuilderAndSolver->InitializeSolutionStep(BaseType::GetModelPart(),rStiffnessMatrix,rDx,rb);

        // initial operations ... things that are constant over the solution step
        pScheme->InitializeSolutionStep(BaseType::GetModelPart(),rStiffnessMatrix,rDx,rb);

        if (BaseType::GetEchoLevel() > 2 && rank == 0)
            std::cout << "Exiting InitializeSolutionStep() of EigensolverStrategy" << std::endl;

        KRATOS_CATCH("")
    }

    /// Check whether initial input is valid.
    virtual int Check()
    {
        KRATOS_TRY

        auto& rModelPart = BaseType::GetModelPart();
        const auto rank = rModelPart.GetCommunicator().MyPID();

        if (BaseType::GetEchoLevel() > 2 && rank == 0)
            std::cout << "Entering Check() of EigensolverStrategy" << std::endl;

        // check the model part
        BaseType::Check();

        // check the scheme
        this->pGetScheme()->Check(rModelPart);

        // check the builder and solver
        this->pGetBuilderAndSolver()->Check(rModelPart);

        if (BaseType::GetEchoLevel() > 2 && rank == 0)
            std::cout << "Exiting Check() of EigensolverStrategy" << std::endl;

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

    bool mInitializeWasPerformed;

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
    void ApplyDirichletConditions(SparseMatrixType& rA, double Factor)
    {
        KRATOS_TRY

        const auto rank = BaseType::GetModelPart().GetCommunicator().MyPID();

        if (BaseType::GetEchoLevel() > 2 && rank == 0)
            std::cout << "Entering ApplyDirichletConditions() of EigensolverStrategy" << std::endl;

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
        #pragma omp parallel for firstprivate(SystemSize)
        for(int k = 0; k < static_cast<int>(SystemSize); ++k)
        {
            std::size_t ColBegin = ARowIndices[k];
            std::size_t ColEnd = ARowIndices[k+1];
            bool empty = true;
            for (auto j = ColBegin; j < ColEnd; ++j)
                if(AValues[j] != 0.0)
                {
                    empty = false;
                    break;
                }
            if(empty == true)
                rA(k,k) = 1.0;
        }

        #pragma omp parallel for
        for (int k = 0; k < static_cast<int>(SystemSize); ++k)
        {
            std::size_t ColBegin = ARowIndices[k];
            std::size_t ColEnd = ARowIndices[k+1];
            if (ScalingFactors[k] == 0.0)
            {
                // row dof is fixed. zero off-diagonal columns and factor diagonal
                for (std::size_t j = ColBegin; j < ColEnd; ++j)
                    if (static_cast<int>(AColIndices[j]) != k)
                        AValues[j] = 0.0;
                    else
                        AValues[j] *= Factor;
            }
            else
            {
                // row dof is not fixed. zero columns associated with fixed dofs
                for (std::size_t j = ColBegin; j < ColEnd; ++j)
                    AValues[j] *= ScalingFactors[AColIndices[j]];
            }
        }

        if (BaseType::GetEchoLevel() > 2 && rank == 0)
            std::cout << "Exiting ApplyDirichletConditions() of EigensolverStrategy" << std::endl;

        KRATOS_CATCH("")
    }

    /// Assign eigenvalues and eigenvectors to kratos variables.
    void AssignVariables(DenseVectorType& rEigenvalues, DenseMatrixType& rEigenvectors)
    {
        auto& rModelPart = BaseType::GetModelPart();
        const auto NumEigenvalues = rEigenvalues.size();

        // store eigenvalues in process info
        rModelPart.GetProcessInfo()[EIGENVALUE_VECTOR] = rEigenvalues;

        for (ModelPart::NodeIterator itNode = rModelPart.NodesBegin(); itNode!= rModelPart.NodesEnd(); itNode++)
        {
            ModelPart::NodeType::DofsContainerType& NodeDofs = itNode->GetDofs();
            const auto NumNodeDofs = NodeDofs.size();
            Matrix& rNodeEigenvectors = itNode->GetValue(EIGENVECTOR_MATRIX);
            if (rNodeEigenvectors.size1() != NumEigenvalues || rNodeEigenvectors.size2() != NumNodeDofs)
                rNodeEigenvectors.resize(NumEigenvalues,NumNodeDofs,false);

            // the jth column index of EIGENVECTOR_MATRIX corresponds to the jth nodal dof. therefore,
            // the dof ordering must not change.
            if (NodeDofs.IsSorted() == false)
                NodeDofs.Sort();

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

