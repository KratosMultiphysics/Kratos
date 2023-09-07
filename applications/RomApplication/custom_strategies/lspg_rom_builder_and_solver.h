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
//                  Eduard Gomez Escandell
//

#pragma once

/* System includes */

/* External includes */
#include "concurrentqueue/concurrentqueue.h"
#include <Eigen/Core>
#include <Eigen/Dense>

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/schemes/scheme.h"
#include "custom_strategies/global_rom_builder_and_solver.h"
#include "utilities/builtin_timer.h"
#include "utilities/reduction_utilities.h"
#include "utilities/dense_householder_qr_decomposition.h"
#include "processes/find_nodal_neighbours_process.h"

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
 * @class LeastSquaresPetrovGalerkinROMBuilderAndSolver
 * @ingroup RomApplication
 * @brief This class provides an implementation for the LeastSquaresPetrovGalerkinROM builder and solver operations.
 * This B&S now inherits from the GlobalROMBuilderAndSolver, which in turn inherits 
 * from the ResidualBasedBlockBuilderAndSolver. The RHS is composed of unbalanced loads 
 * (residual) and is constructed using the ResidualBasedBlockBuilderAndSolver. Similarly, the 
 * LHS is constructed using the ResidualBasedBlockBuilderAndSolver and is then multiplied 
 * by the ROM RIGHT BASIS. This results in a rectangular system with dimensions of FOM size 
 * by ROM size. This system can be solved using either the normal equations or the QR 
 * decomposition. The degrees of freedom are rearranged so that the restrained ones are 
 * placed at the end of the system, ordered inversely to the DofSet, mirroring the arrangement 
 * in the FOM.
 * @tparam TSparseSpace Defines the sparse system under consideration
 * @tparam TDenseSpace Defines the dense system under consideration
 * @tparam TLinearSolver Specifies the linear solver being utilized
 * @author Sebastian Ares de Parga Regalado
 */

template <class TSparseSpace, class TDenseSpace, class TLinearSolver>
class LeastSquaresPetrovGalerkinROMBuilderAndSolver : public GlobalROMBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>
{
public:

    ///@name Type Definitions
    ///@{

    // Class pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(LeastSquaresPetrovGalerkinROMBuilderAndSolver);

    // The size_t types
    using SizeType = std::size_t;
    using IndexType = std::size_t;

    /// The definition of the current class
    using ClassType = LeastSquaresPetrovGalerkinROMBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>;

    /// Definition of the classes from the base class
    using BaseBuilderAndSolverType = BuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>;
    using BaseType = GlobalROMBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>;
    using TSchemeType = typename BaseType::TSchemeType;
    using TSystemMatrixType = typename BaseType::TSystemMatrixType;
    using TSystemVectorType = typename BaseType::TSystemVectorType;
    using LocalSystemVectorType = typename BaseType::LocalSystemVectorType;
    using LocalSystemMatrixType = typename BaseType::LocalSystemMatrixType;
    using TSystemMatrixPointerType = typename BaseType::TSystemMatrixPointerType;
    using TSystemVectorPointerType = typename BaseType::TSystemVectorPointerType;
    using ElementsArrayType = typename BaseType::ElementsArrayType;
    using ConditionsArrayType = typename BaseType::ConditionsArrayType;

    /// Additional definitions
    using MasterSlaveConstraintContainerType = typename ModelPart::MasterSlaveConstraintContainerType;
    using EquationIdVectorType = Element::EquationIdVectorType;

    // Eigen definitions
    using EigenDynamicMatrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    using EigenDynamicVector = Eigen::Matrix<double, Eigen::Dynamic, 1>;
    using EigenSparseMatrix = Eigen::SparseMatrix<double, Eigen::RowMajor, int>;

    /// DoF types definition
    using NodeType = Node;
    using DofType = typename NodeType::DofType;
    using DofPointerType = typename DofType::Pointer;
    using DofQueue = moodycamel::ConcurrentQueue<DofType::Pointer>;
    using DofsVectorType = Element::DofsVectorType;


    ///@}
    ///@name Life cycle
    ///@{

    explicit LeastSquaresPetrovGalerkinROMBuilderAndSolver(
        typename TLinearSolver::Pointer pNewLinearSystemSolver,
        Parameters ThisParameters) : BaseType(pNewLinearSystemSolver) 
    {
        // Validate and assign defaults
        Parameters this_parameters_copy = ThisParameters.Clone();
        this_parameters_copy = this->ValidateAndAssignParameters(this_parameters_copy, this->GetDefaultParameters());
        this->AssignSettings(this_parameters_copy);
    } 

    ~LeastSquaresPetrovGalerkinROMBuilderAndSolver() = default;

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

        BuildAndApplyDirichletConditions(pScheme, rModelPart, rA, rb, rDx);

        // if (rA.size1() != BaseType::mEquationSystemSize || rA.size2() != BaseType::mEquationSystemSize) {
        //     rA.resize(BaseType::mEquationSystemSize, BaseType::mEquationSystemSize, false);
        //     BaseType::ConstructMatrixStructure(pScheme, rA, rModelPart);
        // }

        // this->Build(pScheme, rModelPart, rA, rb);

        // BaseType::ApplyDirichletConditions(pScheme, rModelPart, rA, rDx, rb);

        // Initialize the mask vector with zeros
        Vector hrom_dof_mask_vector = ZeroVector(BaseType::GetEquationSystemSize());

        ProjectROM(rModelPart, rA, rb);

        double time = assembling_timer.ElapsedSeconds();
        KRATOS_INFO_IF("LeastSquaresPetrovGalerkinROMBuilderAndSolver", (this->GetEchoLevel() > 0)) << "Build and project time: " << time << std::endl;

        KRATOS_CATCH("")
    }

    void BuildAndApplyDirichletConditions(
        typename TSchemeType::Pointer pScheme,
        ModelPart &rModelPart,
        TSystemMatrixType &rA,
        TSystemVectorType &rb,
        TSystemVectorType &rDx
    ){
        if (rA.size1() != BaseType::mEquationSystemSize || rA.size2() != BaseType::mEquationSystemSize) {
            rA.resize(BaseType::mEquationSystemSize, BaseType::mEquationSystemSize, false);
            BaseType::ConstructMatrixStructure(pScheme, rA, rModelPart);
        }

        this->Build(pScheme, rModelPart, rA, rb);

        BaseType::ApplyDirichletConditions(pScheme, rModelPart, rA, rDx, rb);
    }

    void BuildRightRomBasis(
        const ModelPart& rModelPart,
        Matrix& rPhiGlobal
    )
    {
        this->BuildRightROMBasis(rModelPart, rPhiGlobal);
    }

    /**
     * Projects the reduced system of equations 
     */
    virtual void ProjectROM(
        ModelPart &rModelPart,
        TSystemMatrixType &rA,
        TSystemVectorType &rb) override
    {
        KRATOS_TRY

        if (mRightRomBasisInitialized==false){
            mPhiGlobal = ZeroMatrix(BaseBuilderAndSolverType::GetEquationSystemSize(), this->GetNumberOfROMModes());
            mRightRomBasisInitialized = true;
        }

        this->BuildRightROMBasis(rModelPart, mPhiGlobal);

        auto a_wrapper = UblasWrapper<double>(rA);
        const auto& eigen_rA = a_wrapper.matrix();
        Eigen::Map<EigenDynamicVector> eigen_rb(rb.data().begin(), rb.size());
        Eigen::Map<EigenDynamicMatrix> eigen_mPhiGlobal(mPhiGlobal.data().begin(), mPhiGlobal.size1(), mPhiGlobal.size2());
        
        EigenDynamicMatrix eigen_rA_times_mPhiGlobal = eigen_rA * eigen_mPhiGlobal; //TODO: Make it in parallel.

        if (mSolvingTechnique == "normal_equations"){
            // Compute the matrix multiplication
            this->mEigenRomA = eigen_rA_times_mPhiGlobal.transpose() * eigen_rA_times_mPhiGlobal; //TODO: Make it in parallel.
            this->mEigenRomB = eigen_rA_times_mPhiGlobal.transpose() * eigen_rb; //TODO: Make it in parallel.
        }
        else if (mSolvingTechnique == "qr_decomposition") {
            this->mEigenRomA = eigen_rA_times_mPhiGlobal; 
            this->mEigenRomB = eigen_rb;
        }

        KRATOS_CATCH("")
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

        // Obtain the assembled residuals vector (To build a basis for Petrov-Galerkin)
        if (mTrainPetrovGalerkinFlag & mBasisStrategy=="residuals"){
            std::stringstream matrix_market_vector_name;
            matrix_market_vector_name << "R_" << rModelPart.GetProcessInfo()[TIME] << "_" << rModelPart.GetProcessInfo()[NL_ITERATION_NUMBER] <<  ".res.mm";
            SparseSpaceType::WriteMatrixMarketVector((matrix_market_vector_name.str()).c_str(), b);
        }
        
        if (mSolvingTechnique == "normal_equations")
            this->SolveROM(rModelPart, this->mEigenRomA, this->mEigenRomB, Dx);
        else if (mSolvingTechnique == "qr_decomposition")
            SolveROM(rModelPart, this->mEigenRomA, this->mEigenRomB, Dx);

        KRATOS_CATCH("")
    }

    Parameters GetDefaultParameters() const override
    {
        Parameters default_parameters = Parameters(R"(
        {
            "name" : "lspg_rom_builder_and_solver",
            "nodal_unknowns" : [],
            "number_of_rom_dofs" : 10,
            "inner_rom_settings": {
                "train_petrov_galerkin" : false,
                "solving_technique" : "normal_equations",
                "basis_strategy" : "residuals",
                "monotonicity_preserving" : false
            }
        })");
        default_parameters.AddMissingParameters(BaseType::GetDefaultParameters());

        return default_parameters;
    }

    static std::string Name() 
    {
        return "lspg_rom_builder_and_solver";
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
        return "LeastSquaresPetrovGalerkinROMBuilderAndSolver";
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
    
    SizeType mNodalDofs;
    ElementsArrayType mComplementaryElements;
    ConditionsArrayType mComplementaryConditions;
    ElementsArrayType mNeighbouringAndSelectedElements;
    ConditionsArrayType mNeighbouringAndSelectedConditions;
    EigenDynamicMatrix mEigenRomA;
    EigenDynamicVector mEigenRomB;
    Matrix mPhiGlobal;
    bool mRightRomBasisInitialized = false;

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
        mTrainPetrovGalerkinFlag = ThisParameters["inner_rom_settings"]["train_petrov_galerkin"].GetBool();
        mBasisStrategy = ThisParameters["inner_rom_settings"]["basis_strategy"].GetString();
        mSolvingTechnique = ThisParameters["inner_rom_settings"]["solving_technique"].GetString();
    }

    /**
     * Resizes a Matrix if it's not the right size
     */
    template<typename TMatrix>
    static void ResizeIfNeeded(TMatrix& rMat, const SizeType Rows, const SizeType Cols)

    {
        if(rMat.size1() != Rows || rMat.size2() != Cols) {
            rMat.resize(Rows, Cols, false);
        }
    };

    /**
     * Solves reduced system of equations and broadcasts it
     */
    void SolveROM(
        ModelPart &rModelPart,
        EigenDynamicMatrix &rEigenRomA,
        EigenDynamicVector &rEigenRomB,
        TSystemVectorType &rDx) override
    {
        KRATOS_TRY
        
        const auto solving_timer = BuiltinTimer();

        LocalSystemVectorType dxrom(this->GetNumberOfROMModes());

        using EigenDynamicVector = Eigen::Matrix<double, Eigen::Dynamic, 1>;
        Eigen::Map<EigenDynamicVector> dxrom_eigen(dxrom.data().begin(), dxrom.size());
        dxrom_eigen = rEigenRomA.colPivHouseholderQr().solve(rEigenRomB);
        KRATOS_WATCH(rEigenRomA)
        KRATOS_WATCH(rEigenRomB)
        KRATOS_WATCH(dxrom_eigen)
        KRATOS_INFO_IF("LeastSquaresPetrovGalerkinROMBuilderAndSolver", (this->GetEchoLevel() > 0)) << "Solve reduced system time: " << solving_timer.ElapsedSeconds() << std::endl;

        // Save the ROM solution increment in the root modelpart database
        auto& r_root_mp = rModelPart.GetRootModelPart();
        noalias(r_root_mp.GetValue(ROM_SOLUTION_INCREMENT)) += dxrom;

        // project reduced solution back to full order model
        const auto backward_projection_timer = BuiltinTimer();
        this->ProjectToFineBasis(dxrom, rModelPart, rDx);
        KRATOS_WATCH(rDx)
        KRATOS_INFO_IF("LeastSquaresPetrovGalerkinROMBuilderAndSolver", (this->GetEchoLevel() > 0)) << "Project to fine basis time: " << backward_projection_timer.ElapsedSeconds() << std::endl;

        KRATOS_CATCH("")
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
    bool mTrainPetrovGalerkinFlag = false;
    std::string mBasisStrategy;
    std::string mSolvingTechnique;

    ///@}
    ///@name Private operations 
    ///@{

    ///@}
}; /* Class LeastSquaresPetrovGalerkinROMBuilderAndSolver */

///@}
///@name Type Definitions
///@{


///@}

} /* namespace Kratos.*/

