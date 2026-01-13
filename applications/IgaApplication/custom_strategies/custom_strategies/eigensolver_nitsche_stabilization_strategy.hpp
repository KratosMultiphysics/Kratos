//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//   License:        BSD License
//   Kratos default license: kratos/license.txt

#if !defined(KRATOS_EIGENSOLVER_NITSCHE_STABILIZATION_STRATEGY )
#define  KRATOS_EIGENSOLVER_NITSCHE_STABILIZATION_STRATEGY

// System includes

// External includes

// Project includes
#include "solving_strategies/strategies/implicit_solving_strategy.h"
#include "utilities/builtin_timer.h"
#include "spaces/ublas_space.h"

#include "containers/model.h"

// Application includes
#include "custom_utilities/iga_flags.h"
#include "iga_application_variables.h"

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

/// Strategy for solving generalized eigenvalue problems to obtain Nitsche stabilization factor.
template<class TSparseSpace,
         class TDenseSpace,
         class TLinearSolver
         >
class EigensolverNitscheStabilizationStrategy
    : public ImplicitSolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(EigensolverNitscheStabilizationStrategy);

    typedef ImplicitSolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

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
    EigensolverNitscheStabilizationStrategy(
        ModelPart& rModelPart,
        SchemePointerType pScheme,
        BuilderAndSolverPointerType pBuilderAndSolver
        )
        : BaseType(rModelPart),
            mpScheme(pScheme), mpBuilderAndSolver(pBuilderAndSolver)
    {
        KRATOS_TRY

        // ensure initialization of system matrices in InitializeSolutionStep()
        mpBuilderAndSolver->SetDofSetIsInitializedFlag(false);

        // default echo level (mute)
        this->SetEchoLevel(0);

        // default rebuild level (build at each solution step)
        this->SetRebuildLevel(1);

        SparseMatrixType* AuxStiffnessMatrix = new SparseMatrixType;
        mpStiffnessMatrix = Kratos::shared_ptr<SparseMatrixType>(AuxStiffnessMatrix);
        SparseMatrixType* AuxStabilizationMatrix = new SparseMatrixType;
        mpStabilizationMatrix = Kratos::shared_ptr<SparseMatrixType>(AuxStabilizationMatrix);

        KRATOS_CATCH("")
    }

    /// Deleted copy constructor.
    EigensolverNitscheStabilizationStrategy(const EigensolverNitscheStabilizationStrategy& Other) = delete;

    /// Destructor.
    ~EigensolverNitscheStabilizationStrategy() override
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

    SparseMatrixType& GetStiffnessMatrix()
    {
        return *mpStiffnessMatrix;
    }

    SparseMatrixType& GetStabilizationMatrix()
    {
        return *mpStabilizationMatrix;
    }

    SparseMatrixPointerType& pGetStiffnessMatrix()
    {
        return mpStiffnessMatrix;
    }

    SparseMatrixPointerType& pGetStabilizationMatrix()
    {
        return mpStabilizationMatrix;
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

        KRATOS_INFO_IF("EigensolverNitscheStabilizationStrategy", BaseType::GetEchoLevel() > 2 && rank == 0)
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

        KRATOS_INFO_IF("EigensolverNitscheStabilizationStrategy", BaseType::GetEchoLevel() > 2 && rank == 0)
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

        if (this->pGetStiffnessMatrix() != nullptr)
            this->pGetStiffnessMatrix() = nullptr;

        if (this->pGetStabilizationMatrix() != nullptr)
            this->pGetStabilizationMatrix() = nullptr;

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

        KRATOS_INFO_IF("EigensolverNitscheStabilizationStrategy", BaseType::GetEchoLevel() > 2 && rank == 0)
            <<  "Entering InitializeSolutionStep" << std::endl;

        BuilderAndSolverPointerType& pBuilderAndSolver = this->pGetBuilderAndSolver();
        SchemePointerType& pScheme = this->pGetScheme();
        SparseMatrixPointerType& pStabilizationMatrix = this->pGetStabilizationMatrix();
        SparseMatrixType& rStabilizationMatrix = this->GetStabilizationMatrix();

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
            SparseMatrixPointerType& pStiffnessMatrix = this->pGetStiffnessMatrix();

            // Stiffness matrix
            pBuilderAndSolver->ResizeAndInitializeVectors(
                pScheme, pStiffnessMatrix, pDx, pb, rModelPart);

            // Stabilization matrix
            pBuilderAndSolver->ResizeAndInitializeVectors(
                pScheme, pStabilizationMatrix, pDx, pb, rModelPart);

            KRATOS_INFO_IF("System Matrix Resize Time", BaseType::GetEchoLevel() > 0 && rank == 0)
                << system_matrix_resize_time.ElapsedSeconds() << std::endl;
        }
        else
        {
            SparseSpaceType::Resize(rb, SparseSpaceType::Size1(rStabilizationMatrix));
            SparseSpaceType::Set(rb, 0.0);
            SparseSpaceType::Resize(rDx, SparseSpaceType::Size1(rStabilizationMatrix));
            SparseSpaceType::Set(rDx, 0.0);
        }

        KRATOS_INFO_IF("System Construction Time", BaseType::GetEchoLevel() > 0 && rank == 0)
            << system_construction_time.ElapsedSeconds() << std::endl;

        // Initial operations ... things that are constant over the solution
        // step
        pBuilderAndSolver->InitializeSolutionStep(BaseType::GetModelPart(),
                                                  rStabilizationMatrix, rDx, rb);

        // Initial operations ... things that are constant over the solution
        // step
        pScheme->InitializeSolutionStep(BaseType::GetModelPart(), rStabilizationMatrix, rDx, rb);

        KRATOS_INFO_IF("EigensolverNitscheStabilizationStrategy", BaseType::GetEchoLevel() > 2 && rank == 0)
            <<  "Exiting InitializeSolutionStep" << std::endl;

        KRATOS_CATCH("")
    }

    /**
     * Solve all the required operations that should be done (for each step).
     * Solve for eigenvalues for the Nitsche stabilization factor.
     */
    bool SolveSolutionStep() override
    {
        ModelPart& rModelPart = BaseType::GetModelPart();

        SchemePointerType& pScheme = this->pGetScheme();
        SparseMatrixType& rStiffnessMatrix = this->GetStiffnessMatrix();
        SparseMatrixType& rStabilizationMatrix = this->GetStabilizationMatrix();

        // Initialize dummy rhs vector
        SparseVectorType b;
        SparseSpaceType::Resize(b,SparseSpaceType::Size1(rStiffnessMatrix));
        SparseSpaceType::Set(b,0.0);

        // Get the global stiffness matrix of the respective model part, linked with the Nitsche stabilization scheme.
        rModelPart.GetProcessInfo()[BUILD_LEVEL] = 1;
        TSparseSpace::SetToZero(rStiffnessMatrix);
        this->pGetBuilderAndSolver()->Build(pScheme,rModelPart,rStiffnessMatrix,b);

        if (BaseType::GetEchoLevel() == 4) {
            TSparseSpace::WriteMatrixMarketMatrix("StiffnessMatrix.mm", rStiffnessMatrix, false);
        }

        // Get the global stabilization matrix of the respective model part, linked with the Nitsche stabilization scheme.
        rModelPart.GetProcessInfo()[BUILD_LEVEL] = 2;
        TSparseSpace::SetToZero(rStabilizationMatrix);
        this->pGetBuilderAndSolver()->Build(pScheme,rModelPart,rStabilizationMatrix,b);

        if (BaseType::GetEchoLevel() == 4) {
            TSparseSpace::WriteMatrixMarketMatrix("StabilizationMatrix.mm", rStabilizationMatrix, false);
        }

        // Obtain the stifness and stabilization matrices on the current interface boundary:
        SparseMatrixType reduced_stabilization_matrix;
        SparseMatrixType reduced_stiffness_matrix;

        Vector rResult;
        std::size_t number_of_nodes_master, number_of_nodes_slave, number_of_nodes;
        if(rModelPart.ConditionsBegin()->GetGeometry().NumberOfGeometryParts() != 0) // Coupling Nitsche condition
        {
            // 1. find the DOFs on the current interface boundary
            Model reduced_model_master;
            ModelPart& reduced_model_part_master = reduced_model_master.CreateModelPart("new_model");

            for (auto& r_cond : rModelPart.Conditions()) {
                auto& r_geom_master = r_cond.GetGeometry().GetGeometryPart(0);
                auto& r_N_master = r_geom_master.ShapeFunctionsValues();

                for (IndexType i = 0; i<r_N_master.size2();++i)
                {
                    if(r_N_master(0,i) > 1e-6)
                    {
                        reduced_model_part_master.AddNode(r_geom_master.pGetPoint(i));
                    }
                }
            }
            number_of_nodes_master = reduced_model_part_master.NumberOfNodes();

            Model reduced_model_slave;
            ModelPart& reduced_model_part_slave = reduced_model_slave.CreateModelPart("new_model");

            for (auto& r_cond : rModelPart.Conditions()) {
                auto& r_geom_slave = r_cond.GetGeometry().GetGeometryPart(1);
                auto& r_N_slave = r_geom_slave.ShapeFunctionsValues();

                for (IndexType i = 0; i<r_N_slave.size2();++i)
                {
                    if(r_N_slave(0,i) > 1e-6)
                    {
                        reduced_model_part_slave.AddNode(r_geom_slave.pGetPoint(i));
                    }
                }
            }
            number_of_nodes_slave = reduced_model_part_slave.NumberOfNodes();

            // 2. create the result vector
            rResult.resize((number_of_nodes_master+number_of_nodes_slave)*3);

            IndexType i_master = 0;
            for (auto& r_node : reduced_model_part_master.Nodes()) {
                const IndexType index = i_master * 3;

                rResult[index]     = r_node.GetDof(DISPLACEMENT_X).EquationId();
                rResult[index + 1] = r_node.GetDof(DISPLACEMENT_Y).EquationId();
                rResult[index + 2] = r_node.GetDof(DISPLACEMENT_Z).EquationId();

                i_master++;
            }

            IndexType i_slave = 0;
            for (auto& r_node : reduced_model_part_slave.Nodes()) {
                const IndexType index = i_slave * 3 + 3 * number_of_nodes_master;

                rResult[index]     = r_node.GetDof(DISPLACEMENT_X).EquationId();
                rResult[index + 1] = r_node.GetDof(DISPLACEMENT_Y).EquationId();
                rResult[index + 2] = r_node.GetDof(DISPLACEMENT_Z).EquationId();

                i_slave++;
            }

            // 3. assigned the value of the reduced stifness and stabilization matrices based on result vector
            reduced_stabilization_matrix = ZeroMatrix((number_of_nodes_master+number_of_nodes_slave)*3,(number_of_nodes_master+number_of_nodes_slave)*3);

            for (IndexType i = 0; i < (number_of_nodes_master+number_of_nodes_slave)*3; i++)
            {
                for (IndexType j = 0; j <= i; j++)
                {
                    auto value = rStabilizationMatrix(rResult(i), rResult(j));
                    reduced_stabilization_matrix(i,j) = value;
                    if (i != j)
                    {
                        reduced_stabilization_matrix(j,i) = value;
                    }
                }
            }

            reduced_stiffness_matrix = ZeroMatrix((number_of_nodes_master+number_of_nodes_slave)*3,(number_of_nodes_master+number_of_nodes_slave)*3);

            for (IndexType i = 0; i < (number_of_nodes_master+number_of_nodes_slave)*3; i++)
            {
                for (IndexType j = 0; j <= i; j++)
                {
                    reduced_stiffness_matrix(i,j) = rStiffnessMatrix(rResult(i),rResult(j));
                    if (i != j)
                    {
                        reduced_stiffness_matrix(j,i) = rStiffnessMatrix(rResult(i),rResult(j));
                    }
                }
            }
        }
        else // Support Nitsche condition
        {
            // 1. find the DOFs on the current interface boundary
            Model reduced_model;
            ModelPart& reduced_model_part = reduced_model.CreateModelPart("new_model");

            for (auto& r_cond : rModelPart.Conditions()) {
                auto& r_geom = r_cond.GetGeometry();
                auto& r_N = r_geom.ShapeFunctionsValues();

                for (IndexType i = 0; i<r_N.size2();++i)
                {
                    if(r_N(0,i) > 1e-6)
                    {
                        reduced_model_part.AddNode(r_geom.pGetPoint(i));
                    }
                }
            }
            number_of_nodes = reduced_model_part.NumberOfNodes();

            // 2. create the result vector
            rResult.resize((number_of_nodes)*3);

            IndexType i = 0;
            for (auto& r_node : reduced_model_part.Nodes()) {
                const IndexType index = i * 3;

                rResult[index]     = r_node.GetDof(DISPLACEMENT_X).EquationId();
                rResult[index + 1] = r_node.GetDof(DISPLACEMENT_Y).EquationId();
                rResult[index + 2] = r_node.GetDof(DISPLACEMENT_Z).EquationId();

                i++;
            }

            // 3. assigned the value of the reduced stifness and stabilization matrices based on result vector
            reduced_stabilization_matrix = ZeroMatrix((number_of_nodes)*3,(number_of_nodes)*3);

            for (IndexType i = 0; i < (number_of_nodes)*3; i++)
            {
                for (IndexType j = 0; j <= i; j++)
                {
                    reduced_stabilization_matrix(i,j) = rStabilizationMatrix(rResult(i),rResult(j));
                    if (i != j)
                        reduced_stabilization_matrix(j,i) = rStabilizationMatrix(rResult(i),rResult(j));
                }
            }

            reduced_stiffness_matrix = ZeroMatrix((number_of_nodes)*3,(number_of_nodes)*3);

            for (IndexType i = 0; i < (number_of_nodes)*3; i++)
            {
                for (IndexType j = 0; j <= i; j++)
                {
                    reduced_stiffness_matrix(i,j) = rStiffnessMatrix(rResult(i),rResult(j));
                    if (i != j)
                        reduced_stiffness_matrix(j,i) = rStiffnessMatrix(rResult(i),rResult(j));
                }
            }
        }

        // Eigenvector matrix and eigenvalue vector are initialized by the solver
        DenseVectorType Eigenvalues;
        DenseMatrixType Eigenvectors;

        // Solve for eigenvalues and eigenvectors
        BuiltinTimer system_solve_time;
        this->pGetBuilderAndSolver()->GetLinearSystemSolver()->Solve(
                reduced_stabilization_matrix,
                reduced_stiffness_matrix,
                Eigenvalues,
                Eigenvectors);

        KRATOS_INFO_IF("System Solve Time", BaseType::GetEchoLevel() > 0)
                << system_solve_time.ElapsedSeconds() << std::endl;

        rModelPart.GetProcessInfo()[EIGENVALUE_NITSCHE_STABILIZATION_VECTOR] = Eigenvalues;

        if (rModelPart.Conditions().begin()->Is(IgaFlags::FIX_ROTATION_X))
        {
            SparseMatrixType& rStabilizationRotationMatrix = this->GetStabilizationMatrix();

            // Get the global stabilization matrix of the respective model part, linked with the Nitsche stabilization scheme.
            rModelPart.GetProcessInfo()[BUILD_LEVEL] = 3;
            TSparseSpace::SetToZero(rStabilizationRotationMatrix);

            this->pGetBuilderAndSolver()->Build(pScheme,rModelPart,rStabilizationRotationMatrix,b);

            if (BaseType::GetEchoLevel() == 4) {
                TSparseSpace::WriteMatrixMarketMatrix("StabilizationRotationMatrix.mm", rStabilizationRotationMatrix, false);
            }

            SparseMatrixType reduced_stabilization_rotation_matrix = ZeroMatrix((number_of_nodes_master+number_of_nodes_slave)*3,(number_of_nodes_master+number_of_nodes_slave)*3);
            for (IndexType i = 0; i < (number_of_nodes_master+number_of_nodes_slave)*3; i++)
            {
                for (IndexType j = 0; j <= i; j++)
                {
                    reduced_stabilization_rotation_matrix(i,j) = rStabilizationRotationMatrix(rResult(i),rResult(j));
                    if (i != j)
                    {
                        reduced_stabilization_rotation_matrix(j,i) = rStabilizationRotationMatrix(rResult(i),rResult(j));
                    }
                }
            }
       
            // Eigenvector matrix and eigenvalue vector are initialized by the solver
            DenseVectorType Eigenvaluesrotation;
            DenseMatrixType Eigenvectorsrotation;

            // // Solve for eigenvalues and eigenvectors
            BuiltinTimer system_solve_time;
            this->pGetBuilderAndSolver()->GetLinearSystemSolver()->Solve(
                reduced_stabilization_rotation_matrix,
                reduced_stiffness_matrix,
                Eigenvaluesrotation,
                Eigenvectorsrotation);

            KRATOS_INFO_IF("System Solve Time", BaseType::GetEchoLevel() > 0)
                << system_solve_time.ElapsedSeconds() << std::endl;

            rModelPart.GetProcessInfo()[EIGENVALUE_NITSCHE_STABILIZATION_ROTATION_VECTOR] = Eigenvaluesrotation;
        }
        else
        {
            DenseVectorType Eigenvaluesrotation(1, 0.0);
            rModelPart.GetProcessInfo()[EIGENVALUE_NITSCHE_STABILIZATION_ROTATION_VECTOR] = Eigenvaluesrotation;
        }

        return true;
    }

    void FinalizeSolutionStep() override
    {
        KRATOS_TRY;

        const int rank = BaseType::GetModelPart().GetCommunicator().MyPID();
        KRATOS_INFO_IF("EigensolverNitscheStabilizationStrategy", BaseType::GetEchoLevel() > 2 && rank == 0)
            <<  "Entering FinalizeSolutionStep" << std::endl;

        SparseMatrixType& rStabilizationMatrix = this->GetStabilizationMatrix();
        SparseVectorPointerType pDx = SparseSpaceType::CreateEmptyVectorPointer();
        SparseVectorPointerType pb = SparseSpaceType::CreateEmptyVectorPointer();
        pGetBuilderAndSolver()->FinalizeSolutionStep(
            BaseType::GetModelPart(), rStabilizationMatrix, *pDx, *pb);
        pGetScheme()->FinalizeSolutionStep(BaseType::GetModelPart(),
                                           rStabilizationMatrix, *pDx, *pb);
        KRATOS_INFO_IF("EigensolverNitscheStabilizationStrategy", BaseType::GetEchoLevel() > 2 && rank == 0)
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

        KRATOS_INFO_IF("EigensolverNitscheStabilizationStrategy", BaseType::GetEchoLevel() > 2 && rank == 0)
            <<  "Entering Check" << std::endl;

        // check the model part
        BaseType::Check();

        // check the scheme
        this->pGetScheme()->Check(rModelPart);

        // check the builder and solver
        this->pGetBuilderAndSolver()->Check(rModelPart);

        KRATOS_INFO_IF("EigensolverNitscheStabilizationStrategy", BaseType::GetEchoLevel() > 2 && rank == 0)
            <<  "Exiting Check" << std::endl;

        return 0;

        KRATOS_CATCH("")
    }

    ///@}

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    SchemePointerType mpScheme;

    BuilderAndSolverPointerType mpBuilderAndSolver;

    SparseMatrixPointerType mpStiffnessMatrix;

    SparseMatrixPointerType mpStabilizationMatrix;

    bool mInitializeWasPerformed = false;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}

}; /* Class EigensolverNitscheStabilizationStrategy */

///@}

///@name Type Definitions
///@{


///@}

} /* namespace Kratos */

#endif /* KRATOS_EIGENSOLVER_NITSCHE_STABILIZATION_STRATEGY  defined */

