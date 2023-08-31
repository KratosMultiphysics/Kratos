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

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/schemes/scheme.h"
#include "custom_strategies/rom_builder_and_solver.h"
#include "utilities/builtin_timer.h"
#include "utilities/reduction_utilities.h"
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
 * @brief Current class provides an implementation for LeastSquaresPetrovGalerkinROM builder and solving operations.
 * @details The RHS is constituted by the unbalanced loads (residual).
 * The LHS is constituted by multiplying the Jacobian or its approximation with the ROM RIGHT BASIS, 
 * yielding a rectangular system (FOM size) that is then solved using the QR decomposition.
 * Degrees of freedom are reordered putting the restrained degrees of freedom at
 * the end of the system ordered in reverse order with respect to the DofSet (as for the FOM).
 * Imposition of the dirichlet conditions is naturally dealt with as the residual already contains
 * this information (as for the FOM).
 * @tparam TSparseSpace The sparse system considered
 * @tparam TDenseSpace The dense system considered
 * @tparam TLinearSolver The linear solver considered
 * @author SebastianAres
 */

template <class TSparseSpace, class TDenseSpace, class TLinearSolver>
class LeastSquaresPetrovGalerkinROMBuilderAndSolver : public ROMBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>
{
public:

    ///@name Type Definitions
    ///@{

    // Class pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(LeastSquaresPetrovGalerkinROMBuilderAndSolver);

    // The size_t types
    typedef std::size_t SizeType;
    typedef std::size_t IndexType;

    /// The definition of the current class
    typedef LeastSquaresPetrovGalerkinROMBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver> ClassType;

    /// Definition of the classes from the base class
    typedef ROMBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;
    typedef typename BaseType::TSchemeType TSchemeType;
    typedef typename BaseType::DofsArrayType DofsArrayType;
    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;
    typedef typename BaseType::TSystemVectorType TSystemVectorType;
    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;
    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;
    typedef typename BaseType::TSystemMatrixPointerType TSystemMatrixPointerType;
    typedef typename BaseType::TSystemVectorPointerType TSystemVectorPointerType;
    typedef typename BaseType::ElementsArrayType ElementsArrayType;
    typedef typename BaseType::ConditionsArrayType ConditionsArrayType;

    /// Definition of the classes from the ROM base class
    // typedef ROMBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver> ROMBuilderAndSolverType;

    /// Additional definitions
    typedef typename ModelPart::MasterSlaveConstraintContainerType MasterSlaveConstraintContainerType;
    typedef Element::EquationIdVectorType EquationIdVectorType;
    typedef Element::DofsVectorType DofsVectorType;

    // Non-distributed, dense:
    typedef LocalSystemMatrixType RomSystemMatrixType;
    typedef LocalSystemVectorType RomSystemVectorType;

    //Distributed, dense
    typedef RomSystemMatrixType LSPGSystemMatrixType; 
    typedef RomSystemVectorType LSPGSystemVectorType;
    //      ^ Change this to a distributed dense type

    /// DoF types definition
    typedef Node NodeType;
    typedef typename NodeType::DofType DofType;
    typedef typename DofType::Pointer DofPointerType;
    typedef moodycamel::ConcurrentQueue<DofType::Pointer> DofQueue;

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

    void SetUpDofSet(
        typename TSchemeType::Pointer pScheme,
        ModelPart &rModelPart) override
    {
        KRATOS_TRY;

        KRATOS_INFO_IF("LeastSquaresPetrovGalerkinROMBuilderAndSolver", (this->GetEchoLevel() > 1)) << "Setting up the dofs" << std::endl;
        KRATOS_INFO_IF("LeastSquaresPetrovGalerkinROMBuilderAndSolver", (this->GetEchoLevel() > 2)) << "Number of threads" << ParallelUtilities::GetNumThreads() << "\n" << std::endl;
        KRATOS_INFO_IF("LeastSquaresPetrovGalerkinROMBuilderAndSolver", (this->GetEchoLevel() > 2)) << "Initializing element loop" << std::endl;

        // Get model part data
        if (this->mHromWeightsInitialized == false) {
            this->InitializeHROMWeights(rModelPart);
        }

        auto dof_queue = this->ExtractDofSet(pScheme, rModelPart);

        // Fill a sorted auxiliary array of with the DOFs set
        KRATOS_INFO_IF("LeastSquaresPetrovGalerkinROMBuilderAndSolver", (this->GetEchoLevel() > 2)) << "Initializing ordered array filling\n" << std::endl;
        auto dof_array = this->SortAndRemoveDuplicateDofs(dof_queue);

        // Update base builder and solver DOFs array and set corresponding flag
        BaseType::GetDofSet().swap(dof_array);
        BaseType::SetDofSetIsInitializedFlag(true);

        // Throw an exception if there are no DOFs involved in the analysis
        KRATOS_ERROR_IF(BaseType::GetDofSet().size() == 0) << "No degrees of freedom!" << std::endl;
        KRATOS_INFO_IF("LeastSquaresPetrovGalerkinROMBuilderAndSolver", (this->GetEchoLevel() > 2)) << "Number of degrees of freedom:" << BaseType::GetDofSet().size() << std::endl;
        KRATOS_INFO_IF("LeastSquaresPetrovGalerkinROMBuilderAndSolver", (this->GetEchoLevel() > 2)) << "Finished setting up the dofs" << std::endl;

#ifdef KRATOS_DEBUG
        // If reactions are to be calculated, we check if all the dofs have reactions defined
        if (BaseType::GetCalculateReactionsFlag())
        {
            for (const auto& r_dof: BaseType::GetDofSet())
            {
                KRATOS_ERROR_IF_NOT(r_dof.HasReaction())
                    << "Reaction variable not set for the following :\n"
                    << "Node : " << r_dof.Id() << '\n'
                    << "Dof  : " << r_dof      << '\n'
                    << "Not possible to calculate reactions." << std::endl;
            }
        }
#endif
        KRATOS_CATCH("");
    } 
    
    void BuildAndSolve(
        typename TSchemeType::Pointer pScheme,
        ModelPart &rModelPart,
        TSystemMatrixType &A,
        TSystemVectorType &Dx,
        TSystemVectorType &b) override
    {
        KRATOS_TRY
        LSPGSystemMatrixType Arom = ZeroMatrix(BaseType::GetEquationSystemSize(), this->GetNumberOfROMModes());
        LSPGSystemVectorType brom = ZeroVector(BaseType::GetEquationSystemSize());
        BuildROM(pScheme, rModelPart, Arom, brom);

        // Obtain the assembled residuals vector (To build a basis for Petrov-Galerkin)
        if (mTrainPetrovGalerkinFlag){
            TSystemVectorType r_residual;
            GetAssembledResiduals(pScheme, rModelPart, r_residual);
        }

        SolveROM(rModelPart, Arom, brom, Dx);


        KRATOS_CATCH("")
    }

    Parameters GetDefaultParameters() const override
    {
        Parameters default_parameters = Parameters(R"(
        {
            "name" : "lspg_rom_builder_and_solver",
            "nodal_unknowns" : [],
            "number_of_rom_dofs" : 10,
            "train_petrov_galerkin" : false
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

    SizeType mNodalDofs;
    bool mTrainPetrovGalerkinFlag = false;

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
        mTrainPetrovGalerkinFlag = ThisParameters["train_petrov_galerkin"].GetBool();
    }

    /**
    * Thread Local Storage containing dynamically allocated structures to avoid reallocating each iteration.
    */
    struct AssemblyTLS
    {
        Matrix phiE = {};                // Elemental Phi
        LocalSystemMatrixType lhs = {};  // Elemental LHS
        EquationIdVectorType eq_id = {}; // Elemental equation ID vector
        DofsVectorType dofs = {};        // Elemental dof vector
        RomSystemMatrixType romA;        // reduced LHS
        RomSystemVectorType romB;        // reduced RHS
    };

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
     * Builds the reduced system of equations on rank 0 
     */
    void BuildROM(
        typename TSchemeType::Pointer pScheme,
        ModelPart &rModelPart,
        LSPGSystemMatrixType &rA,
        LSPGSystemVectorType &rb) override
    {
        KRATOS_TRY
        // Define a dense matrix to hold the reduced problem
        rA = ZeroMatrix(BaseType::GetEquationSystemSize(), this->GetNumberOfROMModes());
        rb = ZeroVector(BaseType::GetEquationSystemSize());

        // Build the system matrix by looping over elements and conditions and assembling to A
        KRATOS_ERROR_IF(!pScheme) << "No scheme provided!" << std::endl;

        // Get ProcessInfo from main model part
        const auto& r_current_process_info = rModelPart.GetProcessInfo();


        // Assemble all entities
        const auto assembling_timer = BuiltinTimer();

        AssemblyTLS assembly_tls_container;

        const auto& r_elements = rModelPart.Elements();

        if(!r_elements.empty())
        {
            block_for_each(r_elements, assembly_tls_container, 
                [&](Element& r_element, AssemblyTLS& r_thread_prealloc)
            {
                CalculateLocalContributionLSPG(r_element, rA, rb, r_thread_prealloc, *pScheme, r_current_process_info);
            });
        }


        const auto& r_conditions = rModelPart.Conditions();

        if(!r_conditions.empty())
        {
            block_for_each(r_conditions, assembly_tls_container, 
                [&](Condition& r_condition, AssemblyTLS& r_thread_prealloc)
            {
                CalculateLocalContributionLSPG(r_condition, rA, rb, r_thread_prealloc, *pScheme, r_current_process_info);
            });
        }

        KRATOS_INFO_IF("LeastSquaresPetrovGalerkinROMBuilderAndSolver", (this->GetEchoLevel() > 0)) << "Build time: " << assembling_timer.ElapsedSeconds() << std::endl;
        KRATOS_INFO_IF("LeastSquaresPetrovGalerkinROMBuilderAndSolver", (this->GetEchoLevel() > 2)) << "Finished parallel building" << std::endl;

        KRATOS_CATCH("")
    }

    /**
     * Solves reduced system of equations and broadcasts it
     */
    void SolveROM(
        ModelPart &rModelPart,
        LSPGSystemMatrixType &rA,
        LSPGSystemVectorType &rb,
        TSystemVectorType &rDx) override
    {
        KRATOS_TRY

        LSPGSystemVectorType dxrom(this->GetNumberOfROMModes());
        
        const auto solving_timer = BuiltinTimer();
        // Calculate the QR decomposition
        DenseHouseholderQRDecomposition<TDenseSpace> qr_decomposition;
        //                              ^Correct after properly defining LSPGSystemMatrixType
        qr_decomposition.Compute(rA);
        qr_decomposition.Solve(rb, dxrom);
        KRATOS_INFO_IF("LeastSquaresPetrovGalerkinROMBuilderAndSolver", (this->GetEchoLevel() > 0)) << "Solve reduced system time: " << solving_timer.ElapsedSeconds() << std::endl;

        // Save the ROM solution increment in the root modelpart database
        auto& r_root_mp = rModelPart.GetRootModelPart();
        noalias(r_root_mp.GetValue(ROM_SOLUTION_INCREMENT)) += dxrom;

        // project reduced solution back to full order model
        const auto backward_projection_timer = BuiltinTimer();
        this->ProjectToFineBasis(dxrom, rModelPart, rDx);
        KRATOS_INFO_IF("LeastSquaresPetrovGalerkinROMBuilderAndSolver", (this->GetEchoLevel() > 0)) << "Project to fine basis time: " << backward_projection_timer.ElapsedSeconds() << std::endl;

        KRATOS_CATCH("")
    }

    /**
     * Builds the reduced system of equations on rank 0 
     */
    void GetAssembledResiduals(
        typename TSchemeType::Pointer pScheme,
        ModelPart &rModelPart,
        LSPGSystemVectorType &rb) 
    {
        KRATOS_TRY
        const auto residual_writing_timer = BuiltinTimer();
        // Define a dense matrix to hold the reduced problem
        rb = ZeroVector(BaseType::GetEquationSystemSize());

        // Build the system matrix by looping over elements and conditions and assembling to A
        KRATOS_ERROR_IF(!pScheme) << "No scheme provided!" << std::endl;

        // Get ProcessInfo from main model part
        const auto& r_current_process_info = rModelPart.GetProcessInfo();

        AssemblyTLS assembly_tls_container;

        const auto& r_elements = rModelPart.Elements();

        if(!r_elements.empty())
        {
            block_for_each(r_elements, assembly_tls_container, 
                [&](Element& r_element, AssemblyTLS& r_thread_prealloc)
            {
                CalculateAssembledResiduals(r_element, rb, r_thread_prealloc, *pScheme, r_current_process_info);
            });
        }


        const auto& r_conditions = rModelPart.Conditions();

        if(!r_conditions.empty())
        {
            block_for_each(r_conditions, assembly_tls_container, 
                [&](Condition& r_condition, AssemblyTLS& r_thread_prealloc)
            {
                CalculateAssembledResiduals(r_condition, rb, r_thread_prealloc, *pScheme, r_current_process_info);
            });
        }

        std::stringstream matrix_market_vector_name;
        matrix_market_vector_name << "R_" << rModelPart.GetProcessInfo()[TIME] << "_" << rModelPart.GetProcessInfo()[NL_ITERATION_NUMBER] <<  ".res.mm";
        SparseSpaceType::WriteMatrixMarketVector((matrix_market_vector_name.str()).c_str(), rb);

        KRATOS_INFO_IF("LeastSquaresPetrovGalerkinROMBuilderAndSolver", (this->GetEchoLevel() > 0)) << "Write residuals to train Petrov Galerkin time: " << residual_writing_timer.ElapsedSeconds() << std::endl;
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

    SizeType mNumberOfRomModes;

    ///@}
    ///@name Private operations 
    ///@{

    /**
     * Computes the local contribution of an element or condition for LSPG
     */
    template<typename TEntity>
    void CalculateLocalContributionLSPG(
        TEntity& rEntity,
        LSPGSystemMatrixType& rAglobal,
        LSPGSystemVectorType& rBglobal,
        AssemblyTLS& rPreAlloc,
        TSchemeType& rScheme,
        const ProcessInfo& rCurrentProcessInfo)
    {
        if (rEntity.IsDefined(ACTIVE) && rEntity.IsNot(ACTIVE)) return;

        // Calculate elemental contribution
        rScheme.CalculateSystemContributions(rEntity, rPreAlloc.lhs, rPreAlloc.romB, rPreAlloc.eq_id, rCurrentProcessInfo);
        rEntity.GetDofList(rPreAlloc.dofs, rCurrentProcessInfo);

        const SizeType ndofs = rPreAlloc.dofs.size();
        ResizeIfNeeded(rPreAlloc.phiE, ndofs, this->GetNumberOfROMModes());
        ResizeIfNeeded(rPreAlloc.romA, ndofs, this->GetNumberOfROMModes());

        const auto &r_geom = rEntity.GetGeometry();
        RomAuxiliaryUtilities::GetPhiElemental(rPreAlloc.phiE, rPreAlloc.dofs, r_geom, this->mMapPhi);

        noalias(rPreAlloc.romA) = prod(rPreAlloc.lhs, rPreAlloc.phiE);


        // Assembly
        for(SizeType row=0; row < ndofs; ++row)
        {
            const SizeType global_row = rPreAlloc.eq_id[row];

            double& r_bi = rBglobal(global_row);
            AtomicAdd(r_bi, rPreAlloc.romB[row]);

            if(rPreAlloc.dofs[row]->IsFree())
            {
                for(SizeType col=0; col < this->GetNumberOfROMModes(); ++col)
                {
                    // const SizeType global_col = rPreAlloc.eq_id[col];
                    const SizeType global_col = col;
                    double& r_aij = rAglobal(global_row, global_col);
                    AtomicAdd(r_aij, rPreAlloc.romA(row, col));
                }
            }
        }
    }

    /**
     * Computes the global assembled residual 
     */
    template<typename TEntity>
    void CalculateAssembledResiduals(
        TEntity& rEntity,
        LSPGSystemVectorType& rBglobal,
        AssemblyTLS& rPreAlloc,
        TSchemeType& rScheme,
        const ProcessInfo& rCurrentProcessInfo)
    {
        if (rEntity.IsDefined(ACTIVE) && rEntity.IsNot(ACTIVE)) return;

        // Calculate elemental contribution
        rScheme.CalculateRHSContribution(rEntity, rPreAlloc.romB, rPreAlloc.eq_id, rCurrentProcessInfo);
        rEntity.GetDofList(rPreAlloc.dofs, rCurrentProcessInfo);

        const SizeType ndofs = rPreAlloc.dofs.size();

        // Assembly
        for(SizeType row=0; row < ndofs; ++row)
        {
            const SizeType global_row = rPreAlloc.eq_id[row];

            if(rPreAlloc.dofs[row]->IsFree())
            {
                double& r_bi = rBglobal(global_row);
                AtomicAdd(r_bi, rPreAlloc.romB[row]);
            }
        }
    }




    ///@}
}; /* Class LeastSquaresPetrovGalerkinROMBuilderAndSolver */

///@}
///@name Type Definitions
///@{


///@}

} /* namespace Kratos.*/

