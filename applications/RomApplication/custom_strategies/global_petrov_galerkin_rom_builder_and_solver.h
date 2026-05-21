//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:        BSD License
//                  Kratos default license: kratos/license.txt
//
//  Main authors:   Sebastian Ares de Parga Regalado
//

#pragma once

/* System includes */

/* External includes */
#include <Eigen/Core>
#include <Eigen/Dense>

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/schemes/scheme.h"
#include "custom_strategies/global_rom_builder_and_solver.h"
#include "utilities/builtin_timer.h"
#include "utilities/dense_householder_qr_decomposition.h"

/* Application includes */
#include "rom_application_variables.h"
#include "custom_utilities/rom_auxiliary_utilities.h"
#include "custom_utilities/rom_residuals_utility.h"

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

/**
 * @class GlobalPetrovGalerkinROMBuilderAndSolver
 * @ingroup RomApplication
 * @brief This class provides an implementation for the GlobalPetrovGalerkinROM builder and solver operations.
 * This Builder and Solver (B&S) now inherits from the GlobalROMBuilderAndSolver, which in turn inherits
 * from the ResidualBasedBlockBuilderAndSolver. The Right-Hand Side (RHS) is composed of unbalanced loads
 * (residual) and is constructed using the ResidualBasedBlockBuilderAndSolver. Similarly, the
 * Left-Hand Side (LHS) is constructed using the ResidualBasedBlockBuilderAndSolver and is then multiplied
 * by the ROM RIGHT BASIS. We then project it onto the ROM LEFT BASIS, yielding a rectangular system (ROM size) that is then
 * solved using the QR decomposition. The degrees of freedom are rearranged so that the restrained ones are
 * placed at the end of the system, ordered inversely to the DofSet, mirroring the arrangement
 * in the Full Order Model (FOM).
 * @tparam TSparseSpace Defines the sparse system under consideration
 * @tparam TDenseSpace Defines the dense system under consideration
 * @tparam TLinearSolver Specifies the linear solver being utilized
 * @author Sebastian Ares de Parga Regalado
 */


template <class TSparseSpace, class TDenseSpace, class TLinearSolver>
class GlobalPetrovGalerkinROMBuilderAndSolver : public GlobalROMBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>
{
public:

    ///@name Type Definitions
    ///@{

    // Class pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(GlobalPetrovGalerkinROMBuilderAndSolver);

    // The size_t types
    using SizeType = std::size_t;
    using IndexType = std::size_t;

    /// The definition of the current class
    using ClassType = GlobalPetrovGalerkinROMBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>;

    /// Definition of the classes from the base class
    using BaseBuilderAndSolverType = BuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>;
    using BaseType = GlobalROMBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>;
    using TSchemeType = typename BaseType::TSchemeType;
    using TSystemMatrixType = typename BaseType::TSystemMatrixType;
    using TSystemVectorType = typename BaseType::TSystemVectorType;
    using ElementsArrayType = typename BaseType::ElementsArrayType;
    using ConditionsArrayType = typename BaseType::ConditionsArrayType;
    using LocalSystemVectorType = typename BaseBuilderAndSolverType::LocalSystemVectorType;
    using LocalSystemMatrixType = typename BaseBuilderAndSolverType::LocalSystemMatrixType;

    // Eigen definitions
    using EigenDynamicMatrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    using EigenDynamicVector = Eigen::Matrix<double, Eigen::Dynamic, 1>;
    using EigenSparseMatrix = Eigen::SparseMatrix<double, Eigen::RowMajor, int>;

    /// DoF types definition
    using DofType = typename Node::DofType;
    using DofPointerType = typename DofType::Pointer;
    using DofQueue = moodycamel::ConcurrentQueue<DofType::Pointer>;

    // Bring the base class function into scope
    using BaseType::ProjectROM;

    ///@}
    ///@name Life cycle
    ///@{

    explicit GlobalPetrovGalerkinROMBuilderAndSolver(
        typename TLinearSolver::Pointer pNewLinearSystemSolver,
        Parameters ThisParameters) : BaseType(pNewLinearSystemSolver)
    {
        // Validate and assign defaults
        Parameters this_parameters_copy = ThisParameters.Clone();
        this_parameters_copy = BaseType::ValidateAndAssignParameters(this_parameters_copy, GetDefaultParameters());
        AssignSettings(this_parameters_copy);
    }

    ~GlobalPetrovGalerkinROMBuilderAndSolver() = default;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /**
     * Builds and projects the reduced system of equations
     */
    virtual void BuildAndProjectROM(
        typename TSchemeType::Pointer pScheme,
        ModelPart &rModelPart,
        TSystemMatrixType &rA,
        TSystemVectorType &rb,
        TSystemVectorType &rDx) override
    {
        KRATOS_TRY

        KRATOS_ERROR_IF(!pScheme) << "No scheme provided!" << std::endl;

        const auto assembling_timer = BuiltinTimer();

        if (rA.size1() != BaseBuilderAndSolverType::mEquationSystemSize || rA.size2() != BaseBuilderAndSolverType::mEquationSystemSize) {
            rA.resize(BaseBuilderAndSolverType::mEquationSystemSize, BaseBuilderAndSolverType::mEquationSystemSize, false);
            BaseType::ConstructMatrixStructure(pScheme, rA, rModelPart);
        }

        BaseType::Build(pScheme, rModelPart, rA, rb);

        if (BaseType::GetMonotonicityPreservingFlag()) {
            BaseType::MonotonicityPreserving(rA, rb);
        }

        ProjectROM(rModelPart, rA, rb);

        double time = assembling_timer.ElapsedSeconds();
        KRATOS_INFO_IF("GlobalPetrovGalerkinROMBuilderAndSolver", (BaseType::GetEchoLevel() > 0)) << "Build and project time: " << time << std::endl;

        KRATOS_CATCH("")
    }

    void GetRightROMBasis(
        const ModelPart& rModelPart,
        Matrix& rPhiGlobal
    )
    {
        BaseType::BuildRightROMBasis(rModelPart, rPhiGlobal);
    }

    /**
     * Projects the reduced system of equations
     */
    virtual void ProjectROM(
        ModelPart &rModelPart,
        TSystemMatrixType &rA,
        TSystemVectorType &rb
        ) override
    {
        KRATOS_TRY

        if (mRightRomBasisInitialized==false){
            mPhiGlobal = ZeroMatrix(BaseBuilderAndSolverType::mEquationSystemSize, BaseType::GetNumberOfROMModes());
            mRightRomBasisInitialized = true;
        }

        if (mLeftRomBasisInitialized==false){
            mPsiGlobal = ZeroMatrix(BaseBuilderAndSolverType::mEquationSystemSize, mNumberOfPetrovGalerkinRomModes);
            mLeftRomBasisInitialized = true;
        }

        BaseType::BuildRightROMBasis(rModelPart, mPhiGlobal);
        BuildLeftROMBasis(rModelPart, mPsiGlobal);
        auto a_wrapper = UblasWrapper<double>(rA);
        const auto& eigen_rA = a_wrapper.matrix();
        Eigen::Map<EigenDynamicVector> eigen_rb(rb.data().begin(), rb.size());
        Eigen::Map<EigenDynamicMatrix> eigen_mPhiGlobal(mPhiGlobal.data().begin(), mPhiGlobal.size1(), mPhiGlobal.size2());
        Eigen::Map<EigenDynamicMatrix> eigen_mPsiGlobal(mPsiGlobal.data().begin(), mPsiGlobal.size1(), mPsiGlobal.size2());

        EigenDynamicMatrix eigen_rA_times_mPhiGlobal = eigen_rA * eigen_mPhiGlobal; //TODO: Make it in parallel.

        // Compute the matrix multiplication
        mEigenRomA = eigen_mPsiGlobal.transpose() * eigen_rA_times_mPhiGlobal; //TODO: Make it in parallel.
        mEigenRomB = eigen_mPsiGlobal.transpose() * eigen_rb; //TODO: Make it in parallel.

        KRATOS_CATCH("")
    }

    void BuildLeftROMBasis(
        const ModelPart& rModelPart,
        Matrix& rPsiGlobal)
    {
        const auto& r_dof_set = BaseBuilderAndSolverType::GetDofSet();
        block_for_each(r_dof_set, [&](const DofType& r_dof)
        {
            const auto& r_node = rModelPart.GetNode(r_dof.Id());
            const Matrix& r_rom_nodal_basis = r_node.GetValue(ROM_LEFT_BASIS);
            const Matrix::size_type row_id = BaseType::mMapPhi.at(r_dof.GetVariable().Key());
            if (r_dof.IsFixed())
            {
                noalias(row(rPsiGlobal, r_dof.EquationId())) = ZeroVector(r_rom_nodal_basis.size2());
            }
            else{
                noalias(row(rPsiGlobal, r_dof.EquationId())) = row(r_rom_nodal_basis, row_id);
            }
        });
    }

    void BuildAndSolve(
        typename TSchemeType::Pointer pScheme,
        ModelPart &rModelPart,
        TSystemMatrixType &A,
        TSystemVectorType &Dx,
        TSystemVectorType &b) override
    {
        KRATOS_TRY

        BuildAndProjectROM(pScheme, rModelPart, A, b, Dx);

        BaseType::SolveROM(rModelPart, mEigenRomA, mEigenRomB, Dx);

        KRATOS_CATCH("")
    }

    Parameters GetDefaultParameters() const override
    {
        Parameters default_parameters = Parameters(R"(
        {
            "name" : "global_petrov_galerkin_rom_builder_and_solver",
            "nodal_unknowns" : [],
            "number_of_rom_dofs" : 10,
            "petrov_galerkin_number_of_rom_dofs" : 10,
            "rom_bns_settings": {
                "monotonicity_preserving" : false
            }
        })");
        default_parameters.AddMissingParameters(BaseType::GetDefaultParameters());

        return default_parameters;
    }

    static std::string Name()
    {
        return "global_petrov_galerkin_rom_builder_and_solver";
    }

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const override
    {
        return "GlobalPetrovGalerkinROMBuilderAndSolver";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream &rOStream) const override
    {
        rOStream << Info();
    }

    /// Print object's
    virtual void PrintData(std::ostream &rOStream) const override
    {
        rOStream << Info();
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}
protected:

    ///@}
    ///@name Protected static member variables
    ///@{


    ///@}
    ///@name Protected member variables
    ///@{

    ///@}
    ///@name Protected operators
    ///@{


    ///@}
    ///@name Protected operations
    ///@{

    void AssignSettings(const Parameters ThisParameters) override
    {
        BaseType::AssignSettings(ThisParameters);

        // // Set member variables
        mNumberOfPetrovGalerkinRomModes = ThisParameters["petrov_galerkin_number_of_rom_dofs"].GetInt();
    }

    ///@}
    ///@name Protected access
    ///@{


    ///@}
    ///@name Protected inquiry
    ///@{


    ///@}
    ///@name Protected life cycle
    ///@{

private:
    SizeType mNumberOfPetrovGalerkinRomModes;
    EigenDynamicMatrix mEigenRomA;
    EigenDynamicVector mEigenRomB;
    Matrix mPhiGlobal;
    Matrix mPsiGlobal;
    bool mRightRomBasisInitialized = false;
    bool mLeftRomBasisInitialized = false;
    std::unordered_set<std::size_t> mSelectedDofs;
    bool mIsSelectedDofsInitialized = false;

    ///@}
    ///@name Private operations
    ///@{

    ///@}
}; /* Class GlobalPetrovGalerkinROMBuilderAndSolver */

///@}
///@name Type Definitions
///@{


///@}

} /* namespace Kratos.*/