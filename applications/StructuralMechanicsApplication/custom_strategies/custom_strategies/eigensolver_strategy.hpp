// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main author:    Michael Andre, https://github.com/msandre
//

#pragma once

// System includes

// External includes

// Project includes
#include "solving_strategies/strategies/implicit_solving_strategy.h"
#include "utilities/builtin_timer.h"
#include "utilities/atomic_utilities.h"
#include "utilities/entities_utilities.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"

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
    : public ImplicitSolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(EigensolverStrategy);

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
    EigensolverStrategy(
        ModelPart& rModelPart,
        SchemePointerType pScheme,
        BuilderAndSolverPointerType pBuilderAndSolver,
        double MassMatrixDiagonalValue,
        double StiffnessMatrixDiagonalValue,
        bool ComputeModalDecomposition = false,
        bool NormalizeEigenvectorsWithMassMatrix = false
        )
        : ImplicitSolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(rModelPart),
            mMassMatrixDiagonalValue(MassMatrixDiagonalValue),
            mStiffnessMatrixDiagonalValue(StiffnessMatrixDiagonalValue),
            mComputeModalDecompostion(ComputeModalDecomposition),
            mNormalizeEigenvectorsWithMassMatrix(NormalizeEigenvectorsWithMassMatrix)
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

        KRATOS_INFO_IF("EigensolverStrategy", BaseType::GetEchoLevel() > 2)
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

        KRATOS_INFO_IF("EigensolverStrategy", BaseType::GetEchoLevel() > 2)
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

        KRATOS_INFO_IF("EigensolverStrategy", BaseType::GetEchoLevel() > 2)
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

            KRATOS_INFO_IF("Setup Dofs Time", BaseType::GetEchoLevel() > 0)
                << setup_dofs_time << std::endl;

            // Set global equation ids
            BuiltinTimer setup_system_time;
            pBuilderAndSolver->SetUpSystem(rModelPart);

            KRATOS_INFO_IF("Setup System Time", BaseType::GetEchoLevel() > 0)
                << setup_system_time << std::endl;

            // Resize and initialize system matrices
            BuiltinTimer system_matrix_resize_time;
            SparseMatrixPointerType& pMassMatrix = this->pGetMassMatrix();

            // Mass matrix
            pBuilderAndSolver->ResizeAndInitializeVectors(
                pScheme, pMassMatrix, pDx, pb, rModelPart);

            // Stiffness matrix
            pBuilderAndSolver->ResizeAndInitializeVectors(
                pScheme, pStiffnessMatrix, pDx, pb, rModelPart);

            KRATOS_INFO_IF("System Matrix Resize Time", BaseType::GetEchoLevel() > 0)
                << system_matrix_resize_time << std::endl;
        }
        else
        {
            SparseSpaceType::Resize(rb, SparseSpaceType::Size1(rStiffnessMatrix));
            SparseSpaceType::Set(rb, 0.0);
            SparseSpaceType::Resize(rDx, SparseSpaceType::Size1(rStiffnessMatrix));
            SparseSpaceType::Set(rDx, 0.0);
        }

        KRATOS_INFO_IF("System Construction Time", BaseType::GetEchoLevel() > 0)
            << system_construction_time << std::endl;

        // Initial operations ... things that are constant over the solution
        // step
        pBuilderAndSolver->InitializeSolutionStep(BaseType::GetModelPart(),
                                                  rStiffnessMatrix, rDx, rb);

        // Initial operations ... things that are constant over the solution
        // step
        pScheme->InitializeSolutionStep(BaseType::GetModelPart(), rStiffnessMatrix, rDx, rb);

        KRATOS_INFO_IF("EigensolverStrategy", BaseType::GetEchoLevel() > 2)
            <<  "Exiting InitializeSolutionStep" << std::endl;

        KRATOS_CATCH("")
    }

    bool SolveSolutionStep() override
    {
        KRATOS_TRY;

        ModelPart& rModelPart = BaseType::GetModelPart();

        SchemePointerType& pScheme = this->pGetScheme();
        SparseMatrixType& rMassMatrix = this->GetMassMatrix();
        SparseMatrixType& rStiffnessMatrix = this->GetStiffnessMatrix();

        // Initialize dummy rhs vector
        SparseVectorType b;
        SparseSpaceType::Resize(b,SparseSpaceType::Size1(rMassMatrix));
        SparseSpaceType::Set(b,0.0);

        // if the size of the Matrix is the same as the size of the dofset, then also the dirichlet-dofs are in the matrix
        // i.e. a BlockBuilder is used
        const bool matrix_contains_dirichlet_dofs = SparseSpaceType::Size1(rMassMatrix) == this->pGetBuilderAndSolver()->GetDofSet().size();

        const bool master_slave_constraints_defined = rModelPart.NumberOfMasterSlaveConstraints() != 0;

        rModelPart.GetProcessInfo()[BUILD_LEVEL] = 1;
        TSparseSpace::SetToZero(rMassMatrix);

        EntitiesUtilities::InitializeNonLinearIterationAllEntities(rModelPart);

        this->pGetBuilderAndSolver()->Build(pScheme,rModelPart,rMassMatrix,b);
        if (master_slave_constraints_defined) {
            this->pGetBuilderAndSolver()->ApplyConstraints(pScheme, rModelPart, rMassMatrix, b);
        }
        if (matrix_contains_dirichlet_dofs) {
            this->ApplyDirichletConditions(rMassMatrix, mMassMatrixDiagonalValue);
        }

        if (BaseType::GetEchoLevel() == 4) {
            TSparseSpace::WriteMatrixMarketMatrix("MassMatrix.mm", rMassMatrix, false);
        }

        rModelPart.GetProcessInfo()[BUILD_LEVEL] = 2;
        TSparseSpace::SetToZero(rStiffnessMatrix);
        this->pGetBuilderAndSolver()->Build(pScheme,rModelPart,rStiffnessMatrix,b);
        if (master_slave_constraints_defined) {
            this->pGetBuilderAndSolver()->ApplyConstraints(pScheme, rModelPart, rStiffnessMatrix, b);
        }

        if (matrix_contains_dirichlet_dofs) {
            this->ApplyDirichletConditions(rStiffnessMatrix, mStiffnessMatrixDiagonalValue);
        }

        EntitiesUtilities::FinalizeNonLinearIterationAllEntities(rModelPart);

        if (BaseType::GetEchoLevel() == 4) {
            TSparseSpace::WriteMatrixMarketMatrix("StiffnessMatrix.mm", rStiffnessMatrix, false);
        }

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

        KRATOS_INFO_IF("System Solve Time", BaseType::GetEchoLevel() > 0)
                << system_solve_time << std::endl;

        if (master_slave_constraints_defined){
            this->ReconstructSlaveSolution(Eigenvectors);
        }

        if (mNormalizeEigenvectorsWithMassMatrix){
            MassNormalizeEigenvectors(Eigenvectors);
        }

        this->AssignVariables(Eigenvalues,Eigenvectors);


        if (mComputeModalDecompostion) {
            ComputeModalDecomposition(Eigenvectors);
        }

        return true;
        KRATOS_CATCH("")
    }

    void FinalizeSolutionStep() override
    {
        KRATOS_TRY;

        KRATOS_INFO_IF("EigensolverStrategy", BaseType::GetEchoLevel() > 2)
            <<  "Entering FinalizeSolutionStep" << std::endl;

        SparseMatrixType& rStiffnessMatrix = this->GetStiffnessMatrix();
        SparseVectorPointerType pDx = SparseSpaceType::CreateEmptyVectorPointer();
        SparseVectorPointerType pb = SparseSpaceType::CreateEmptyVectorPointer();
        pGetBuilderAndSolver()->FinalizeSolutionStep(
            BaseType::GetModelPart(), rStiffnessMatrix, *pDx, *pb);
        pGetScheme()->FinalizeSolutionStep(BaseType::GetModelPart(),
                                           rStiffnessMatrix, *pDx, *pb);
        KRATOS_INFO_IF("EigensolverStrategy", BaseType::GetEchoLevel() > 2)
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

        KRATOS_INFO_IF("EigensolverStrategy", BaseType::GetEchoLevel() > 2)
            <<  "Entering Check" << std::endl;

        // check the model part
        BaseType::Check();

        // check the scheme
        this->pGetScheme()->Check(rModelPart);

        // check the builder and solver
        this->pGetBuilderAndSolver()->Check(rModelPart);

        KRATOS_INFO_IF("EigensolverStrategy", BaseType::GetEchoLevel() > 2)
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

    double mMassMatrixDiagonalValue = 0.0;
    double mStiffnessMatrixDiagonalValue = 1.0;

    bool mComputeModalDecompostion = false;

    bool mNormalizeEigenvectorsWithMassMatrix = false; 

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief Recover the solution related to the slave dofs in case of master-slave constraints.
     * @details This function is an adaptation of the implementation in ConstraintUtilities,
     *  since there the variables are assumed to be stored in SolutionStepValue.
     *  Beware that this implementation is only valid for Block B&S, since the master-slave constraints
     *  don't work with Elimination B&S yet.
     */
    void ReconstructSlaveSolution(
        DenseMatrixType& rEigenvectors
    )
    {
        KRATOS_TRY

        auto& r_model_part = BaseType::GetModelPart();
        const std::size_t number_of_eigenvalues = rEigenvectors.size1();

        struct TLS{
            Matrix relation_matrix;
            Vector constant_vector;
            Vector master_dofs_values;
        };

        for (std::size_t i_eigenvalue = 0; i_eigenvalue < number_of_eigenvalues; ++i_eigenvalue){

            // Reset slave dofs
            block_for_each(r_model_part.MasterSlaveConstraints(), [&i_eigenvalue, &rEigenvectors](const MasterSlaveConstraint& r_master_slave_constraint){
                const auto& r_slave_dofs_vector = r_master_slave_constraint.GetSlaveDofsVector();
                for (const auto& r_slave_dof: r_slave_dofs_vector){
                    AtomicMult(rEigenvectors(i_eigenvalue, r_slave_dof->EquationId()), 0.0);
                }
            });

            // Apply constraints
            block_for_each(r_model_part.MasterSlaveConstraints(), TLS(), [&i_eigenvalue, &rEigenvectors, &r_model_part](const MasterSlaveConstraint& r_master_slave_constraint, TLS& rTLS){
                // Detect if the constraint is active or not. If the user did not make any choice the constraint
                // It is active by default
                bool constraint_is_active = true;
                if (r_master_slave_constraint.IsDefined(ACTIVE))
                    constraint_is_active = r_master_slave_constraint.Is(ACTIVE);
                if (constraint_is_active) {
                    // Saving the master dofs values
                    const auto& r_master_dofs_vector = r_master_slave_constraint.GetMasterDofsVector();
                    const auto& r_slave_dofs_vector = r_master_slave_constraint.GetSlaveDofsVector();
                    rTLS.master_dofs_values.resize(r_master_dofs_vector.size());
                    for (IndexType i = 0; i < r_master_dofs_vector.size(); ++i) {
                        rTLS.master_dofs_values[i] = rEigenvectors(i_eigenvalue, r_master_dofs_vector[i]->EquationId());
                    }
                    // Apply the constraint to the slave dofs
                    r_master_slave_constraint.GetLocalSystem(rTLS.relation_matrix, rTLS.constant_vector, r_model_part.GetProcessInfo());
                    double aux;
                    for (IndexType i = 0; i < rTLS.relation_matrix.size1(); ++i) {
                        aux = rTLS.constant_vector[i];
                        for(IndexType j = 0; j < rTLS.relation_matrix.size2(); ++j) {
                            aux += rTLS.relation_matrix(i,j) * rTLS.master_dofs_values[j];
                        }
                        AtomicAdd(rEigenvectors(i_eigenvalue, r_slave_dofs_vector[i]->EquationId()), aux);
                    }
                }
            });
        }

        KRATOS_CATCH("")
    }

    /// Normalize the eigenmodes with respect to the mass matrix
    void MassNormalizeEigenvectors(DenseMatrixType& rEigenvectors)
    {
        const SparseMatrixType& rMassMatrix = this->GetMassMatrix();
        const std::size_t num_modes = rEigenvectors.size1();
        const std::size_t num_dofs  = rEigenvectors.size2();

        DenseVectorType phi(num_dofs);
        DenseVectorType tmp(num_dofs);

        for (std::size_t j = 0; j < num_modes; ++j) {
             // copy row j -> phi
            IndexPartition<std::size_t>(num_dofs).for_each([&phi, &rEigenvectors, j](std::size_t i){
                phi[i] = rEigenvectors(j, i);
            });

            // tmp = M * phi_j
            TSparseSpace::Mult(rMassMatrix, phi, tmp);

            // Compute modal_mass = phi^T * tmp
            const double modal_mass = IndexPartition<std::size_t>(num_dofs).for_each<SumReduction<double>>([&phi, &tmp](std::size_t i)->double {
                return phi[i] * tmp[i];
            });
            const auto global_modal_mass = BaseType::GetModelPart().GetCommunicator().GetDataCommunicator().SumAll(modal_mass);

            if (!(global_modal_mass > 0.0)) {
                KRATOS_WARNING("MassNormalizeEigenvectors")
                    << "Modal mass for mode " << j << " is non-positive (" << global_modal_mass << "). Skipping normalization.\n";
                continue;
            }

            const double scale = 1.0 / std::sqrt(global_modal_mass);

           // scale and write back into rEigenvectors 
            IndexPartition<std::size_t>(num_dofs).for_each([&rEigenvectors, &phi, scale, j](std::size_t i){
                rEigenvectors(j, i) = phi[i] * scale;
            });
        }
    }

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

        KRATOS_INFO_IF("EigensolverStrategy", BaseType::GetEchoLevel() > 2)
            <<  "Entering ApplyDirichletConditions" << std::endl;

        const std::size_t SystemSize = rA.size1();
        std::vector<double> ScalingFactors(SystemSize);
        auto& rDofSet = this->pGetBuilderAndSolver()->GetDofSet();
        const int NumDofs = static_cast<int>(rDofSet.size());

        // NOTE: dofs are assumed to be numbered consecutively
        IndexPartition(NumDofs).for_each([&rDofSet, &ScalingFactors](std::size_t k){
            auto dof_iterator = std::begin(rDofSet) + k;
            ScalingFactors[k] = (dof_iterator->IsFixed()) ? 0.0 : 1.0;
        });

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

        IndexPartition(SystemSize).for_each([&](std::size_t k){
            std::size_t ColBegin = ARowIndices[k];
            std::size_t ColEnd = ARowIndices[k+1];
            if (ScalingFactors[k] == 0.0)
            {
                // row dof is fixed. zero off-diagonal columns and factor diagonal
                for (std::size_t j = ColBegin; j < ColEnd; ++j)
                {
                    if (AColIndices[j] != k)
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
        });

        KRATOS_INFO_IF("EigensolverStrategy", BaseType::GetEchoLevel() > 2)
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

        const auto& r_dof_set = this->pGetBuilderAndSolver()->GetDofSet();

        for (ModelPart::NodeIterator itNode = rModelPart.NodesBegin(); itNode!= rModelPart.NodesEnd(); itNode++) {
            ModelPart::NodeType::DofsContainerType& NodeDofs = itNode->GetDofs();
            const std::size_t NumNodeDofs = NodeDofs.size();
            Matrix& rNodeEigenvectors = itNode->GetValue(EIGENVECTOR_MATRIX);
            if (rNodeEigenvectors.size1() != NumEigenvalues || rNodeEigenvectors.size2() != NumNodeDofs) {
                rNodeEigenvectors.resize(NumEigenvalues,NumNodeDofs,false);
            }

            // fill the EIGENVECTOR_MATRIX
            for (std::size_t i = 0; i < NumEigenvalues; i++) {
                for (std::size_t j = 0; j < NumNodeDofs; j++)
                {
                    const auto itDof = std::begin(NodeDofs) + j;
                    bool is_active = !(r_dof_set.find(**itDof) == r_dof_set.end());
                    if ((*itDof)->IsFree() && is_active) {
                       rNodeEigenvectors(i,j) = rEigenvectors(i,(*itDof)->EquationId());
                    }
                    else {
                       rNodeEigenvectors(i,j) = 0.0;
                    }
                }
            }
        }
    }

    ///
     /**
     * Computes the modal decomposition depending on the number of eigenvalues
     * chosen and stores them in the corresponding variables. Can be activated by setting
     * bool variable exposed to the python interface.
     */
    void ComputeModalDecomposition(const DenseMatrixType& rEigenvectors)
    {
        const SparseMatrixType& rMassMatrix = this->GetMassMatrix();
        SparseMatrixType m_temp = ZeroMatrix(rEigenvectors.size1(),rEigenvectors.size2());
        boost::numeric::ublas::axpy_prod(rEigenvectors,rMassMatrix,m_temp,true);
        Matrix modal_mass_matrix = ZeroMatrix(m_temp.size1(),m_temp.size1());
        boost::numeric::ublas::axpy_prod(m_temp,trans(rEigenvectors),modal_mass_matrix);

        const SparseMatrixType& rStiffnessMatrix = this->GetStiffnessMatrix();
        SparseMatrixType k_temp = ZeroMatrix(rEigenvectors.size1(),rEigenvectors.size2());
        boost::numeric::ublas::axpy_prod(rEigenvectors,rStiffnessMatrix,k_temp,true);
        Matrix modal_stiffness_matrix = ZeroMatrix(k_temp.size1(),k_temp.size1());
        boost::numeric::ublas::axpy_prod(k_temp,trans(rEigenvectors),modal_stiffness_matrix);

        ModelPart& rModelPart = BaseType::GetModelPart();
        rModelPart.GetProcessInfo()[MODAL_MASS_MATRIX] = modal_mass_matrix;
        rModelPart.GetProcessInfo()[MODAL_STIFFNESS_MATRIX] = modal_stiffness_matrix;

        KRATOS_INFO("ModalMassMatrix")      << modal_mass_matrix << std::endl;
        KRATOS_INFO("ModalStiffnessMatrix") << modal_stiffness_matrix << std::endl;
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
