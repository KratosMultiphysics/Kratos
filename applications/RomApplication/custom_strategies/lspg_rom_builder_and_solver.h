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

    // Eigen definitions
    using EigenDynamicMatrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    using EigenDynamicVector = Eigen::Matrix<double, Eigen::Dynamic, 1>;
    using EigenSparseMatrix = Eigen::SparseMatrix<double, Eigen::RowMajor, int>;

    /// DoF types definition
    using NodeType = Node;


    ///@}
    ///@name Life cycle
    ///@{

    explicit LeastSquaresPetrovGalerkinROMBuilderAndSolver(
        typename TLinearSolver::Pointer pNewLinearSystemSolver,
        Parameters ThisParameters) : BaseType(pNewLinearSystemSolver) 
    {
        // Validate and assign defaults
        Parameters this_parameters_copy = ThisParameters.Clone();
        this_parameters_copy = BaseType::ValidateAndAssignParameters(this_parameters_copy, GetDefaultParameters());
        AssignSettings(this_parameters_copy);
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

        if (BaseType::GetMonotonicityPreservingFlag()) {
            BaseType::MonotonicityPreserving(rA, rb);
        }

        ProjectROM(rModelPart, rA, rb);

        double time = assembling_timer.ElapsedSeconds();
        KRATOS_INFO_IF("LeastSquaresPetrovGalerkinROMBuilderAndSolver", (BaseType::GetEchoLevel() > 0)) << "Build and project time: " << time << std::endl;

        KRATOS_CATCH("")
    }

    void BuildAndApplyDirichletConditions(
        typename TSchemeType::Pointer pScheme,
        ModelPart &rModelPart,
        TSystemMatrixType &rA,
        TSystemVectorType &rb,
        TSystemVectorType &rDx
    ){
        if (rA.size1() != BaseBuilderAndSolverType::mEquationSystemSize || rA.size2() != BaseBuilderAndSolverType::mEquationSystemSize) {
            rA.resize(BaseBuilderAndSolverType::mEquationSystemSize, BaseBuilderAndSolverType::mEquationSystemSize, false);
            BaseType::ConstructMatrixStructure(pScheme, rA, rModelPart);
        }

        BaseType::Build(pScheme, rModelPart, rA, rb);

        BaseType::ApplyDirichletConditions(pScheme, rModelPart, rA, rDx, rb);
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
        TSystemVectorType &rb) override
    {
        KRATOS_TRY

        if (mRightRomBasisInitialized==false){
            mPhiGlobal = ZeroMatrix(BaseBuilderAndSolverType::mEquationSystemSize, BaseType::GetNumberOfROMModes());
            mRightRomBasisInitialized = true;
        }

        BaseType::BuildRightROMBasis(rModelPart, mPhiGlobal);
        auto a_wrapper = UblasWrapper<double>(rA);
        const auto& eigen_rA = a_wrapper.matrix();
        Eigen::Map<EigenDynamicVector> eigen_rb(rb.data().begin(), rb.size());
        Eigen::Map<EigenDynamicMatrix> eigen_mPhiGlobal(mPhiGlobal.data().begin(), mPhiGlobal.size1(), mPhiGlobal.size2());

        EigenDynamicMatrix eigen_rA_times_mPhiGlobal = eigen_rA * eigen_mPhiGlobal; //TODO: Make it in parallel.

        if (mSolvingTechnique == "normal_equations"){
            // Compute the matrix multiplication
            mEigenRomA = eigen_rA_times_mPhiGlobal.transpose() * eigen_rA_times_mPhiGlobal; //TODO: Make it in parallel.
            mEigenRomB = eigen_rA_times_mPhiGlobal.transpose() * eigen_rb; //TODO: Make it in parallel.
        }
        else if (mSolvingTechnique == "qr_decomposition") {
            mEigenRomA = eigen_rA_times_mPhiGlobal; 
            mEigenRomB = eigen_rb;
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
        
        if (mSolvingTechnique == "normal_equations"){
            BaseType::SolveROM(rModelPart, mEigenRomA, mEigenRomB, Dx);
        }
        else if (mSolvingTechnique == "qr_decomposition"){
            BaseType::SolveROM(rModelPart, mEigenRomA, mEigenRomB, Dx);
        }

        KRATOS_CATCH("")
    }

    Parameters GetDefaultParameters() const override
    {
        Parameters default_parameters = Parameters(R"(
        {
            "name" : "lspg_rom_builder_and_solver",
            "nodal_unknowns" : [],
            "number_of_rom_dofs" : 10,
            "rom_bns_settings": {
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
        mTrainPetrovGalerkinFlag = ThisParameters["rom_bns_settings"]["train_petrov_galerkin"].GetBool();
        mBasisStrategy = ThisParameters["rom_bns_settings"]["basis_strategy"].GetString();
        mSolvingTechnique = ThisParameters["rom_bns_settings"]["solving_technique"].GetString();
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
    EigenDynamicMatrix mEigenRomA;
    EigenDynamicVector mEigenRomB;
    Matrix mPhiGlobal;
    bool mRightRomBasisInitialized = false;

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

