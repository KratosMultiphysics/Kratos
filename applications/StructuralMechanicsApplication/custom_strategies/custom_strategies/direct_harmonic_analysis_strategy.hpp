#pragma once

// System includes
#include <cmath>
#include <complex>
#include <vector>

// External includes
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

// Project includes
#include "solving_strategies/strategies/implicit_solving_strategy.h"
#include "linear_solvers/linear_solver.h"
#include "spaces/ublas_space.h"
#include "utilities/builtin_timer.h"
#include "utilities/atomic_utilities.h"
#include "utilities/entities_utilities.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"

// Application includes
#include "structural_mechanics_application_variables.h"

namespace Kratos
{

template<class TSparseSpace, class TDenseSpace, class TLinearSolver>
class DirectHarmonicAnalysisStrategy
    : public ImplicitSolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(DirectHarmonicAnalysisStrategy);

    using BaseType = ImplicitSolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>;
    using SchemePointerType = typename BaseType::TSchemeType::Pointer;
    using BuilderAndSolverPointerType = typename BaseType::TBuilderAndSolverType::Pointer;

    using SparseSpaceType = TSparseSpace;
    using DenseSpaceType  = TDenseSpace;

    using SparseMatrixType        = typename TSparseSpace::MatrixType;
    using SparseMatrixPointerType = typename TSparseSpace::MatrixPointerType;
    using SparseVectorType        = typename TSparseSpace::VectorType;
    using SparseVectorPointerType = typename TSparseSpace::VectorPointerType;

    using DenseVectorType = typename TDenseSpace::VectorType;
    using DenseMatrixType = typename TDenseSpace::MatrixType;

    using RealType = double;
    using ComplexType = std::complex<double>;

    using ComplexSparseMatrixType =
        boost::numeric::ublas::compressed_matrix<ComplexType>;

    using ComplexVectorType =
        boost::numeric::ublas::vector<ComplexType>;

    using ComplexDenseMatrixType =
        boost::numeric::ublas::matrix<ComplexType>;

    using ComplexSparseSpaceType =
        UblasSpace<ComplexType, ComplexSparseMatrixType, ComplexVectorType>;

    using ComplexDenseSpaceType =
        UblasSpace<ComplexType, ComplexDenseMatrixType, ComplexVectorType>;

    using ComplexLinearSolverType =
        LinearSolver<ComplexSparseSpaceType, ComplexDenseSpaceType>;

    using ComplexLinearSolverPointerType =
        typename ComplexLinearSolverType::Pointer;

    using ComplexSparseMatrixPointerType =
        Kratos::shared_ptr<ComplexSparseMatrixType>;

    using ComplexVectorPointerType =
        Kratos::shared_ptr<ComplexVectorType>;

    using LocalSystemMatrixType = typename BaseType::LocalSystemMatrixType;
    using LocalSystemVectorType = typename BaseType::LocalSystemVectorType;

    /// Constructor
    DirectHarmonicAnalysisStrategy(
        ModelPart& rModelPart,
        SchemePointerType pScheme,
        BuilderAndSolverPointerType pBuilderAndSolver,
        ComplexLinearSolverPointerType pComplexLinearSolver,
        const double MassMatrixDiagonalValue = 1.0,
        const double StiffnessMatrixDiagonalValue = 1.0,
        const double DampingMatrixDiagonalValue = 1.0,
        const bool ReformDofSetAtEachStep = false,
        const bool AssembleDampingMatrix = false,
        const std::string RealLoadSubModelPart = "",
        const std::string ImaginaryLoadSubModelPart = ""
    )
        : BaseType(rModelPart),
          mMassMatrixDiagonalValue(MassMatrixDiagonalValue),
          mStiffnessMatrixDiagonalValue(StiffnessMatrixDiagonalValue),
          mDampingMatrixDiagonalValue(DampingMatrixDiagonalValue),
          mAssembleDampingMatrix(AssembleDampingMatrix)
    {
        KRATOS_TRY

        mpScheme = pScheme;
        mpBuilderAndSolver = pBuilderAndSolver;
        mpComplexLinearSolver = pComplexLinearSolver;

        KRATOS_ERROR_IF(mpComplexLinearSolver == nullptr)
            << "Complex linear solver is null." << std::endl;

        mpBuilderAndSolver->SetDofSetIsInitializedFlag(false);
        mpBuilderAndSolver->SetReshapeMatrixFlag(ReformDofSetAtEachStep);

        mpMassMatrix = Kratos::make_shared<SparseMatrixType>();
        mpStiffnessMatrix = Kratos::make_shared<SparseMatrixType>();
        mpDampingMatrix = Kratos::make_shared<SparseMatrixType>();

        mpRealLoadVector = SparseSpaceType::CreateEmptyVectorPointer();
        mpImaginaryLoadVector = SparseSpaceType::CreateEmptyVectorPointer();

        mpDynamicMatrix = Kratos::make_shared<ComplexSparseMatrixType>();
        mpComplexRHS = Kratos::make_shared<ComplexVectorType>();
        mpComplexSolution = Kratos::make_shared<ComplexVectorType>();

        mRealLoadSubModelPart = RealLoadSubModelPart;
        mImaginaryLoadSubModelPart = ImaginaryLoadSubModelPart;

        this->SetEchoLevel(0);
        this->SetRebuildLevel(0);

        KRATOS_CATCH("")
    }

    /// Deleted copy constructor
    DirectHarmonicAnalysisStrategy(const DirectHarmonicAnalysisStrategy& Other) = delete;

    ~DirectHarmonicAnalysisStrategy() override
    {
        this->Clear();
    }

    void SetEchoLevel(int Level) override
    {
        BaseType::SetEchoLevel(Level);
        mpBuilderAndSolver->SetEchoLevel(Level);
    }

    SparseMatrixPointerType& pGetMassMatrix()
    {
        return mpMassMatrix;
    }

    SparseMatrixPointerType& pGetStiffnessMatrix()
    {
        return mpStiffnessMatrix;
    }

    SparseMatrixPointerType& pGetDampingMatrix()
    {
        return mpDampingMatrix;
    }

    void SetAssembleDampingMatrixFlag(const bool AssembleDampingMatrix)
    {
        mAssembleDampingMatrix = AssembleDampingMatrix;
    }

    bool GetAssembleDampingMatrixFlag() const
    {
        return mAssembleDampingMatrix;
    }

    void Initialize() override
    {
        KRATOS_TRY

        ModelPart& r_model_part = BaseType::GetModelPart();

        KRATOS_INFO_IF("DirectHarmonicAnalysisStrategy", BaseType::GetEchoLevel() > 2)
            << "Entering Initialize" << std::endl;

        if (!mInitializeWasPerformed) {
            auto& p_scheme = mpScheme;

            if (!p_scheme->SchemeIsInitialized()) {
                p_scheme->Initialize(r_model_part);
            }

            if (!p_scheme->ElementsAreInitialized()) {
                p_scheme->InitializeElements(r_model_part);
            }

            if (!p_scheme->ConditionsAreInitialized()) {
                p_scheme->InitializeConditions(r_model_part);
            }

            // Read Rayleigh coefficients from properties, if available
            for (auto& r_prop : r_model_part.PropertiesArray()) {
                if (r_prop->Has(RAYLEIGH_ALPHA)) {
                    mRayleighAlpha = r_prop->GetValue(RAYLEIGH_ALPHA);
                }
                if (r_prop->Has(RAYLEIGH_BETA)) {
                    mRayleighBeta = r_prop->GetValue(RAYLEIGH_BETA);
                }
            }

            mInitializeWasPerformed = true;
        }

        KRATOS_INFO_IF("DirectHarmonicAnalysisStrategy", BaseType::GetEchoLevel() > 2)
            << "Exiting Initialize" << std::endl;

        KRATOS_CATCH("")
    }

    void InitializeSolutionStep() override
    {
        KRATOS_TRY

        ModelPart& r_model_part = BaseType::GetModelPart();

        KRATOS_INFO_IF("DirectHarmonicAnalysisStrategy", BaseType::GetEchoLevel() > 2)
            << "Entering InitializeSolutionStep" << std::endl;

        auto& p_builder_and_solver = mpBuilderAndSolver;
        auto& p_scheme = mpScheme;

        auto& p_mass_matrix = this->pGetMassMatrix();
        auto& p_stiffness_matrix = this->pGetStiffnessMatrix();
        auto& p_damping_matrix = this->pGetDampingMatrix();

        SparseVectorPointerType pDx = SparseSpaceType::CreateEmptyVectorPointer();
        SparseVectorPointerType pb  = SparseSpaceType::CreateEmptyVectorPointer();
        auto& rDx = *pDx;
        auto& rb = *pb;

        if (!p_builder_and_solver->GetDofSetIsInitializedFlag() ||
            p_builder_and_solver->GetReshapeMatrixFlag())
        {
            BuiltinTimer setup_dofs_time;
            p_builder_and_solver->SetUpDofSet(p_scheme, r_model_part);

            KRATOS_INFO_IF("Setup Dofs Time", BaseType::GetEchoLevel() > 0)
                << setup_dofs_time << std::endl;

            BuiltinTimer setup_system_time;
            p_builder_and_solver->SetUpSystem(r_model_part);

            KRATOS_INFO_IF("Setup System Time", BaseType::GetEchoLevel() > 0)
                << setup_system_time << std::endl;

            BuiltinTimer resize_time;
            p_builder_and_solver->ResizeAndInitializeVectors(
                p_scheme, p_mass_matrix, pDx, pb, r_model_part);

            p_builder_and_solver->ResizeAndInitializeVectors(
                p_scheme, p_stiffness_matrix, pDx, pb, r_model_part);

            p_builder_and_solver->ResizeAndInitializeVectors(
                p_scheme, p_damping_matrix, pDx, pb, r_model_part);

            const std::size_t system_size = SparseSpaceType::Size1(*p_stiffness_matrix);

            if (mpRealLoadVector->size() != system_size) {
                SparseSpaceType::Resize(*mpRealLoadVector, system_size);
            }
            SparseSpaceType::Set(*mpRealLoadVector, 0.0);
            if (mpImaginaryLoadVector->size() != system_size) {
                SparseSpaceType::Resize(*mpImaginaryLoadVector, system_size);
            }
            SparseSpaceType::Set(*mpImaginaryLoadVector, 0.0);

            if (mpComplexRHS->size() != system_size) {
                mpComplexRHS->resize(system_size, false);
            }
            if (mpComplexSolution->size() != system_size) {
                mpComplexSolution->resize(system_size, false);
            }

            for (std::size_t i = 0; i < system_size; ++i) {
                (*mpComplexRHS)[i] = ComplexType(0.0, 0.0);
                (*mpComplexSolution)[i] = ComplexType(0.0, 0.0);
            }

            if (mpDynamicMatrix->size1() != system_size || mpDynamicMatrix->size2() != system_size) {
                mpDynamicMatrix->resize(system_size, system_size, false);
            }

            KRATOS_INFO_IF("System Resize Time", BaseType::GetEchoLevel() > 0)
                << resize_time << std::endl;
        }
        else
        {
            SparseSpaceType::Resize(rb, SparseSpaceType::Size1(*p_stiffness_matrix));
            SparseSpaceType::Set(rb, 0.0);
            SparseSpaceType::Resize(rDx, SparseSpaceType::Size1(*p_stiffness_matrix));
            SparseSpaceType::Set(rDx, 0.0);
        }

        p_builder_and_solver->InitializeSolutionStep(
            r_model_part, *p_stiffness_matrix, rDx, rb);

        p_scheme->InitializeSolutionStep(
            r_model_part, *p_stiffness_matrix, rDx, rb);

        KRATOS_INFO_IF("DirectHarmonicAnalysisStrategy", BaseType::GetEchoLevel() > 2)
            << "Exiting InitializeSolutionStep" << std::endl;

        KRATOS_CATCH("")
    }

    bool SolveSolutionStep() override
    {
        KRATOS_TRY

        ModelPart& r_model_part = BaseType::GetModelPart();
        auto& r_process_info = r_model_part.GetProcessInfo();

        auto& r_mass_matrix = *(this->pGetMassMatrix());
        auto& r_stiffness_matrix = *(this->pGetStiffnessMatrix());
        auto& r_damping_matrix = *(this->pGetDampingMatrix());
        auto& r_real_load_vector = *mpRealLoadVector;
        auto& r_imaginary_load_vector = *mpImaginaryLoadVector;

        // For now the angular frequency is stored inside TIME
        const double omega = r_process_info[TIME];

        SparseVectorType b;
        SparseSpaceType::Resize(b, SparseSpaceType::Size1(r_mass_matrix));
        SparseSpaceType::Set(b, 0.0);

        const bool master_slave_constraints_defined =
            r_model_part.NumberOfMasterSlaveConstraints() != 0;

        // Assemble M
        r_process_info[BUILD_LEVEL] = 1;
        TSparseSpace::SetToZero(r_mass_matrix);

        EntitiesUtilities::InitializeNonLinearIterationAllEntities(r_model_part);
        mpBuilderAndSolver->Build(mpScheme, r_model_part, r_mass_matrix, b);

        if (master_slave_constraints_defined) {
            mpBuilderAndSolver->ApplyConstraints(mpScheme, r_model_part, r_mass_matrix, b);
        }

        // Assemble K
        r_process_info[BUILD_LEVEL] = 2;
        TSparseSpace::SetToZero(r_stiffness_matrix);
        TSparseSpace::SetToZero(b);

        mpBuilderAndSolver->Build(mpScheme, r_model_part, r_stiffness_matrix, b);

        if (master_slave_constraints_defined) {
            mpBuilderAndSolver->ApplyConstraints(mpScheme, r_model_part, r_stiffness_matrix, b);
        }

        // Assemble C
        TSparseSpace::SetToZero(r_damping_matrix);

        if (mAssembleDampingMatrix) {
            r_process_info[BUILD_LEVEL] = 3;
            TSparseSpace::SetToZero(b);

            mpBuilderAndSolver->Build(mpScheme, r_model_part, r_damping_matrix, b);

            if (master_slave_constraints_defined) {
                mpBuilderAndSolver->ApplyConstraints(mpScheme, r_model_part, r_damping_matrix, b);
            }
        }

        // Assemble RHS  manually to avoid the residual formulation
        TSparseSpace::SetToZero(r_real_load_vector);
        TSparseSpace::SetToZero(r_imaginary_load_vector);

        // Build real part of the RHS (we do it here to avoid the residual formulation)
        if (mRealLoadSubModelPart != ""){
            BuildRHS(r_model_part.GetSubModelPart(mRealLoadSubModelPart), r_real_load_vector);
        }

        // Build imaginary part of the RHS (we do it here to avoid the residual formulation)
        if (mImaginaryLoadSubModelPart != ""){
            BuildRHS(r_model_part.GetSubModelPart(mImaginaryLoadSubModelPart), r_imaginary_load_vector);
        }
        
        // Keep this consistent with your current infrastructure.
        if (master_slave_constraints_defined) {
            mpBuilderAndSolver->ApplyConstraints(mpScheme, r_model_part, r_stiffness_matrix, r_real_load_vector);
            mpBuilderAndSolver->ApplyConstraints(mpScheme, r_model_part, r_stiffness_matrix, r_imaginary_load_vector);
        }

        EntitiesUtilities::FinalizeNonLinearIterationAllEntities(r_model_part);

        // Build complex dynamic system D(om) = -om^2 * M + j * om * C + K
        BuildComplexDynamicSystem(omega);
        ApplyDirichletConditionsComplex(*mpDynamicMatrix, *mpComplexRHS, ComplexType(1.0, 0.0));

        // Solve complex system
        BuiltinTimer solve_time;
        mpComplexLinearSolver->Solve(
            *mpDynamicMatrix,
            *mpComplexSolution,
            *mpComplexRHS);

        KRATOS_INFO_IF("DirectHarmonicAnalysisStrategy", BaseType::GetEchoLevel() > 0)
            << "Direct harmonic solve time: " << solve_time << std::endl;

        // Assign results to DISPLACEMENT and DISPLACEMENT_IMAGINARY
        AssignVariables(*mpComplexSolution);

        return true;

        KRATOS_CATCH("")
    }

    void FinalizeSolutionStep() override
    {
        KRATOS_TRY

        auto& r_stiffness_matrix = *mpStiffnessMatrix;
        SparseVectorPointerType pDx = SparseSpaceType::CreateEmptyVectorPointer();
        SparseVectorPointerType pb  = SparseSpaceType::CreateEmptyVectorPointer();

        mpBuilderAndSolver->FinalizeSolutionStep(
            BaseType::GetModelPart(), r_stiffness_matrix, *pDx, *pb);

        mpScheme->FinalizeSolutionStep(
            BaseType::GetModelPart(), r_stiffness_matrix, *pDx, *pb);

        KRATOS_CATCH("")
    }

    void Clear() override
    {
        KRATOS_TRY

        auto& p_builder_and_solver = mpBuilderAndSolver;
        if (p_builder_and_solver != nullptr) {
            p_builder_and_solver->GetLinearSystemSolver()->Clear();
        }

        if (mpComplexLinearSolver != nullptr) {
            mpComplexLinearSolver->Clear();
        }

        if (mpMassMatrix != nullptr) {
            mpMassMatrix = nullptr;
        }
        if (mpStiffnessMatrix != nullptr) {
            mpStiffnessMatrix = nullptr;
        }
        if (mpDampingMatrix != nullptr) {
            mpDampingMatrix = nullptr;
        }
        if (mpDynamicMatrix != nullptr) {
            mpDynamicMatrix = nullptr;
        }

        SparseSpaceType::Clear(mpRealLoadVector);
        SparseSpaceType::Clear(mpImaginaryLoadVector);

        if (mpComplexRHS != nullptr) {
            mpComplexRHS = nullptr;
        }
        if (mpComplexSolution != nullptr) {
            mpComplexSolution = nullptr;
        }

        if (p_builder_and_solver != nullptr) {
            p_builder_and_solver->SetDofSetIsInitializedFlag(false);
            p_builder_and_solver->Clear();
        }

        if (mpScheme != nullptr) {
            mpScheme->Clear();
        }

        mInitializeWasPerformed = false;
        mRayleighAlpha = 0.0;
        mRayleighBeta = 0.0;

        KRATOS_CATCH("")
    }

    int Check() override
    {
        KRATOS_TRY

        ModelPart& r_model_part = BaseType::GetModelPart();

        BaseType::Check();
        mpScheme->Check(r_model_part);
        mpBuilderAndSolver->Check(r_model_part);

        KRATOS_ERROR_IF(mpComplexLinearSolver == nullptr)
            << "DirectHarmonicAnalysisStrategy: complex linear solver is null in Check()." << std::endl;

        KRATOS_ERROR_IF(mRealLoadSubModelPart.empty())
            << "DirectHarmonicAnalysisStrategy: \"real_load_sub_model_part\" is not defined." << std::endl;

        KRATOS_ERROR_IF_NOT(r_model_part.HasSubModelPart(mRealLoadSubModelPart))
            << "DirectHarmonicAnalysisStrategy: real load sub model part \""
            << mRealLoadSubModelPart << "\" was not found in model part \""
            << r_model_part.Name() << "\"." << std::endl;

        if (!mImaginaryLoadSubModelPart.empty()) {
            KRATOS_ERROR_IF_NOT(r_model_part.HasSubModelPart(mImaginaryLoadSubModelPart))
                << "DirectHarmonicAnalysisStrategy: imaginary load sub model part \""
                << mImaginaryLoadSubModelPart << "\" was not found in model part \""
                << r_model_part.Name() << "\"." << std::endl;
        }

        return 0;

        KRATOS_CATCH("")
    }

private:
    SchemePointerType mpScheme;
    BuilderAndSolverPointerType mpBuilderAndSolver;
    ComplexLinearSolverPointerType mpComplexLinearSolver;

    SparseMatrixPointerType mpMassMatrix;
    SparseMatrixPointerType mpStiffnessMatrix;
    SparseMatrixPointerType mpDampingMatrix;

    SparseVectorPointerType mpRealLoadVector;
    SparseVectorPointerType mpImaginaryLoadVector;

    ComplexSparseMatrixPointerType mpDynamicMatrix;
    ComplexVectorPointerType mpComplexRHS;
    ComplexVectorPointerType mpComplexSolution;

    std::string mRealLoadSubModelPart;
    std::string mImaginaryLoadSubModelPart;

    bool mInitializeWasPerformed = false;

    double mMassMatrixDiagonalValue = 1.0;
    double mStiffnessMatrixDiagonalValue = 1.0;
    double mDampingMatrixDiagonalValue = 1.0;

    bool mAssembleDampingMatrix = false;

    double mRayleighAlpha = 0.0;
    double mRayleighBeta = 0.0;

private:
    void BuildComplexDynamicSystem(const double omega)
    {
        auto& rK = *mpStiffnessMatrix;
        auto& rM = *mpMassMatrix;
        auto& rC = *mpDampingMatrix;
        auto& rF_real = *mpRealLoadVector;
        auto& rF_imaginary = *mpImaginaryLoadVector;

        auto& rA = *mpDynamicMatrix;
        auto& rb = *mpComplexRHS;
        auto& rx = *mpComplexSolution;

        const std::size_t n = SparseSpaceType::Size1(rK);

        // Reset
        rA.clear();
        rA.resize(n, n, false);

        for (std::size_t i = 0; i < n; ++i) {
            rb[i] = ComplexType(rF_real[i], rF_imaginary[i]);
            rx[i] = ComplexType(0.0, 0.0);
        }

        // Real part: K - omega^2 M
        AddRealMatrixToComplex(rA, rK, ComplexType(1.0, 0.0));
        AddRealMatrixToComplex(rA, rM, ComplexType(-omega * omega, 0.0));

        // Imaginary part: i*omega*C
        if (mAssembleDampingMatrix) {
            AddRealMatrixToComplex(rA, rC, ComplexType(0.0, omega));
        }

        // Rayleigh damping contribution: i*omega*(alpha*M + beta*K)
        if (std::abs(mRayleighAlpha) > 0.0) {
            AddRealMatrixToComplex(rA, rM, ComplexType(0.0, omega * mRayleighAlpha));
        }

        if (std::abs(mRayleighBeta) > 0.0) {
            AddRealMatrixToComplex(rA, rK, ComplexType(0.0, omega * mRayleighBeta));
        }
    }

    void AddRealMatrixToComplex(
        ComplexSparseMatrixType& rA,
        const SparseMatrixType& rB,
        const ComplexType Coeff)
    {
        for (typename SparseMatrixType::const_iterator1 it1 = rB.begin1(); it1 != rB.end1(); ++it1) {
            for (typename SparseMatrixType::const_iterator2 it2 = it1.begin(); it2 != it1.end(); ++it2) {
                const std::size_t I = it2.index1();
                const std::size_t J = it2.index2();
                rA(I, J) += Coeff * (*it2);
            }
        }
    }

    void ApplyDirichletConditionsComplex(
        ComplexSparseMatrixType& rA,
        ComplexVectorType& rb,
        const ComplexType DiagonalValue)
    {
        KRATOS_TRY

        const std::size_t system_size = rA.size1();
        std::vector<double> scaling_factors(system_size, 1.0);

        auto& r_dof_set = mpBuilderAndSolver->GetDofSet();
        const int num_dofs = static_cast<int>(r_dof_set.size());

        IndexPartition(num_dofs).for_each([&r_dof_set, &scaling_factors](std::size_t k){
            auto dof_iterator = std::begin(r_dof_set) + k;
            scaling_factors[k] = dof_iterator->IsFixed() ? 0.0 : 1.0;
        });

        // Ensure diagonal exists on fixed rows
        for (std::size_t k = 0; k < system_size; ++k) {
            if (scaling_factors[k] == 0.0) {
                bool has_diagonal = false;

                const std::size_t col_begin = rA.index1_data()[k];
                const std::size_t col_end   = rA.index1_data()[k + 1];

                for (std::size_t j = col_begin; j < col_end; ++j) {
                    if (rA.index2_data()[j] == k) {
                        has_diagonal = true;
                        break;
                    }
                }

                if (!has_diagonal) {
                    rA(k, k) = DiagonalValue;
                }
            }
        }

        auto* AValues = std::begin(rA.value_data());
        auto* ARowIndices = std::begin(rA.index1_data());
        auto* AColIndices = std::begin(rA.index2_data());

        IndexPartition(system_size).for_each([&](std::size_t k){
            const std::size_t col_begin = ARowIndices[k];
            const std::size_t col_end   = ARowIndices[k + 1];

            if (scaling_factors[k] == 0.0) {
                for (std::size_t j = col_begin; j < col_end; ++j) {
                    if (AColIndices[j] != k) {
                        AValues[j] = ComplexType(0.0, 0.0);
                    } else {
                        AValues[j] = DiagonalValue;
                    }
                }
                rb[k] = ComplexType(0.0, 0.0);
            } else {
                for (std::size_t j = col_begin; j < col_end; ++j) {
                    if (scaling_factors[AColIndices[j]] == 0.0) {
                        AValues[j] = ComplexType(0.0, 0.0);
                    }
                }
            }
        });

        KRATOS_CATCH("")
    }
    
    void AssignVariables(const ComplexVectorType& rSolution)
    {
        ModelPart& r_model_part = BaseType::GetModelPart();

        for (auto& r_node : r_model_part.Nodes()) {

            // initialize the full vector to zero first
            auto& r_displacement_imaginary = r_node.FastGetSolutionStepValue(DISPLACEMENT_IMAGINARY);
            r_displacement_imaginary[0] = 0.0;
            r_displacement_imaginary[1] = 0.0;
            r_displacement_imaginary[2] = 0.0;

            auto& r_node_dofs = r_node.GetDofs();

            for (auto it_dof = r_node_dofs.begin(); it_dof != r_node_dofs.end(); ++it_dof) {
                auto& p_dof = *it_dof;
                const auto& r_var = p_dof->GetVariable();

                if (!p_dof->IsFixed()) {
                    const std::size_t eq_id = p_dof->EquationId();
                    const ComplexType u = rSolution[eq_id];

                    const double u_real = std::real(u);
                    const double u_imag = std::imag(u);

                    // real part -> actual displacement DOF value
                    p_dof->GetSolutionStepValue() = u_real;

                    // imaginary part -> MESH_DISPLACEMENT vector
                    if (r_var == DISPLACEMENT_X) {
                        r_displacement_imaginary[0] = u_imag;
                    } else if (r_var == DISPLACEMENT_Y) {
                        r_displacement_imaginary[1] = u_imag;
                    } else if (r_var == DISPLACEMENT_Z) {
                        r_displacement_imaginary[2] = u_imag;
                    }
                } else {
                    p_dof->GetSolutionStepValue() = 0.0;

                    if (r_var == DISPLACEMENT_X) {
                        r_displacement_imaginary[0] = 0.0;
                    } else if (r_var == DISPLACEMENT_Y) {
                        r_displacement_imaginary[1] = 0.0;
                    } else if (r_var == DISPLACEMENT_Z) {
                        r_displacement_imaginary[2] = 0.0;
                    }
                }
            }
        }
    }

    void BuildRHS(ModelPart& rModelPart, Vector& rRHS)
    {
        KRATOS_TRY

        SparseSpaceType::SetToZero(rRHS);

        auto& r_process_info = rModelPart.GetProcessInfo();
        const std::size_t number_of_conditions = rModelPart.NumberOfConditions();

        IndexPartition<std::size_t>(number_of_conditions).for_each(
            std::tuple<Vector, Element::EquationIdVectorType>{},
            [&rModelPart, &r_process_info, &rRHS](const std::size_t Index, auto& rTLS) {
                auto& r_condition = *(rModelPart.ConditionsBegin() + Index);

                if (!r_condition.IsActive()) {
                    return;
                }

                auto& r_local_rhs = std::get<0>(rTLS);
                auto& r_equation_id = std::get<1>(rTLS);

                r_condition.CalculateRightHandSide(r_local_rhs, r_process_info);
                r_condition.EquationIdVector(r_equation_id, r_process_info);

                for (IndexType i = 0; i < r_equation_id.size(); ++i) {
                    const auto eq_id = r_equation_id[i];
                    if (eq_id < rRHS.size()) {
                        AtomicAdd(rRHS[eq_id], r_local_rhs[i]);
                    }
                }
            });

        KRATOS_CATCH("")
    }
};

} // namespace Kratos