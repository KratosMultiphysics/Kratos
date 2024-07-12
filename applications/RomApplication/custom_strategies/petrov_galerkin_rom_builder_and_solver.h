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
 * @class PetrovGalerkinROMBuilderAndSolver
 * @ingroup RomApplication
 * @brief Current class provides an implementation for PetrovGalerkinROM builder and solving operations.
 * @details The RHS is constituted by the unbalanced loads (residual) and projected onto the ROM LEFT BASIS.
 * The LHS is constituted by first multiplying the Jacobian or its approximation with the ROM RIGHT BASIS
 * and then projecting it onto the ROM LEFT BASIS, yielding a rectangular system (ROM size) that is then 
 * solved using the QR decomposition.
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
class PetrovGalerkinROMBuilderAndSolver : public ROMBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>
{
public:


    ///@name Type Definitions
    ///@{

    // Class pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(PetrovGalerkinROMBuilderAndSolver);

    // The size_t types
    typedef std::size_t SizeType;
    typedef std::size_t IndexType;

    /// The definition of the current class
    typedef PetrovGalerkinROMBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver> ClassType;

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

    /// Additional definitions
    typedef typename ModelPart::MasterSlaveConstraintContainerType MasterSlaveConstraintContainerType;
    typedef Element::EquationIdVectorType EquationIdVectorType;
    typedef Element::DofsVectorType DofsVectorType;

    // Non-distributed, dense:
    typedef LocalSystemMatrixType RomSystemMatrixType;
    typedef LocalSystemVectorType RomSystemVectorType;

    //Distributed, dense
    typedef RomSystemMatrixType PetrovGalerkinSystemMatrixType; 
    typedef RomSystemVectorType PetrovGalerkinSystemVectorType;
    //      ^ Change this to a distributed dense type

    /// DoF types definition
    typedef Node NodeType;
    typedef typename NodeType::DofType DofType;
    typedef typename DofType::Pointer DofPointerType;
    typedef moodycamel::ConcurrentQueue<DofType::Pointer> DofQueue;

    ///@}
    ///@name Life cycle
    ///@{

    explicit PetrovGalerkinROMBuilderAndSolver(
        typename TLinearSolver::Pointer pNewLinearSystemSolver,
        Parameters ThisParameters)
        : ROMBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>(pNewLinearSystemSolver) 
    {
        // Validate and assign defaults
        Parameters this_parameters_copy = ThisParameters.Clone();
        this_parameters_copy = this->ValidateAndAssignParameters(this_parameters_copy, this->GetDefaultParameters());
        this->AssignSettings(this_parameters_copy);
    }

    ~PetrovGalerkinROMBuilderAndSolver() = default;

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

        KRATOS_INFO_IF("PetrovGalerkinROMBuilderAndSolver", (this->GetEchoLevel() > 1)) << "Setting up the dofs" << std::endl;
        KRATOS_INFO_IF("PetrovGalerkinROMBuilderAndSolver", (this->GetEchoLevel() > 2)) << "Number of threads" << ParallelUtilities::GetNumThreads() << "\n" << std::endl;
        KRATOS_INFO_IF("PetrovGalerkinROMBuilderAndSolver", (this->GetEchoLevel() > 2)) << "Initializing element loop" << std::endl;

        // Get model part data
        if (this->mHromWeightsInitialized == false) {
            this->InitializeHROMWeights(rModelPart);
        }

        auto dof_queue = this->ExtractDofSet(pScheme, rModelPart);

        // Fill a sorted auxiliary array of with the DOFs set
        KRATOS_INFO_IF("PetrovGalerkinROMBuilderAndSolver", (this->GetEchoLevel() > 2)) << "Initializing ordered array filling\n" << std::endl;
        auto dof_array = this->SortAndRemoveDuplicateDofs(dof_queue);

        // Update base builder and solver DOFs array and set corresponding flag
        BaseType::GetDofSet().swap(dof_array);
        BaseType::SetDofSetIsInitializedFlag(true);

        // Throw an exception if there are no DOFs involved in the analysis
        KRATOS_ERROR_IF(BaseType::GetDofSet().size() == 0) << "No degrees of freedom!" << std::endl;
        KRATOS_INFO_IF("PetrovGalerkinROMBuilderAndSolver", (this->GetEchoLevel() > 2)) << "Number of degrees of freedom:" << BaseType::GetDofSet().size() << std::endl;
        KRATOS_INFO_IF("PetrovGalerkinROMBuilderAndSolver", (this->GetEchoLevel() > 2)) << "Finished setting up the dofs" << std::endl;

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
        // KRATOS_TRY
        PetrovGalerkinSystemMatrixType Arom = ZeroMatrix(mNumberOfPetrovGalerkinRomModes, this->GetNumberOfROMModes());
        PetrovGalerkinSystemVectorType brom = ZeroVector(mNumberOfPetrovGalerkinRomModes);
        BuildROM(pScheme, rModelPart, Arom, brom);
        SolveROM(rModelPart, Arom, brom, Dx);


        // KRATOS_CATCH("")
    }

    Parameters GetDefaultParameters() const override
    {
        Parameters default_parameters = Parameters(R"(
        {
            "name" : "petrov_galerkin_rom_builder_and_solver",
            "nodal_unknowns" : [],
            "number_of_rom_dofs" : 10,
            "petrov_galerkin_number_of_rom_dofs" : 10
        })");
        default_parameters.AddMissingParameters(BaseType::GetDefaultParameters());

        return default_parameters;
    }

    static std::string Name() 
    {
        return "petrov_galerkin_rom_builder_and_solver";
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
        return "PetrovGalerkinROMBuilderAndSolver";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream &rOStream) const override
    {
        rOStream << Info();
    }

    /// Print object's data.GetNumberOfROMModes()
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

    SizeType mNumberOfPetrovGalerkinRomModes;

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

    /**
    * Thread Local Storage containing dynamically allocated structures to avoid reallocating each iteration.
    */
    struct AssemblyTLS
    {
        AssemblyTLS(SizeType NRomModes, SizeType NPetrovGalerkinRomModes)
            : romA(ZeroMatrix(NPetrovGalerkinRomModes, NRomModes)),
              romB(ZeroVector(NPetrovGalerkinRomModes))
        { }
        AssemblyTLS() = delete;

        Matrix phiE = {};                // Elemental Phi
        Matrix psiE = {};                // Elemental Psi
        LocalSystemMatrixType lhs = {};  // Elemental LHS
        LocalSystemVectorType rhs = {};  // Elemental RHS
        EquationIdVectorType eq_id = {}; // Elemental equation ID vector
        DofsVectorType dofs = {};        // Elemental dof vector
        PetrovGalerkinSystemMatrixType romA;        // reduced LHS
        PetrovGalerkinSystemVectorType romB;        // reduced RHS
        RomSystemMatrixType aux = {};    // Auxiliary: romA = psi.t * (LHS * phi) := psi.t * aux
    };

    /**
     * Class to sum-reduce matrices and vectors.
     */
    template<typename T>
    struct NonTrivialSumReduction
    {
        typedef T value_type;
        typedef T return_type;

        T mValue;
        bool mInitialized = false;

        void Init(const value_type& first_value)
        {
            mValue = first_value;
            mInitialized = true;
        }

        /// access to reduced value
        return_type GetValue() const
        {
            return mValue;
        }

        void LocalReduce(const value_type& value)
        {
            if(!mInitialized) {
                Init(value);
            } else {
                noalias(mValue) += value;
            }
        }

        void ThreadSafeReduce(const NonTrivialSumReduction& rOther)
        {
            if(!rOther.mInitialized) return;

            const std::lock_guard<LockObject> scope_lock(ParallelUtilities::GetGlobalLock());
            LocalReduce(rOther.mValue);
        }
    };

    /**
     * Resizes a Matrix if it's not the right size
     */
    template<typename TMatrix>
    static void ResizeIfNeeded(TMatrix& mat, const SizeType rows, const SizeType cols)
    {
        if(mat.size1() != rows || mat.size2() != cols) {
            mat.resize(rows, cols, false);
        }
    };

    /**
     * Builds the reduced system of equations on rank 0 
     */
    void BuildROM(
        typename TSchemeType::Pointer pScheme,
        ModelPart &rModelPart,
        PetrovGalerkinSystemMatrixType &rA,
        PetrovGalerkinSystemVectorType &rb) override
    {
        // KRATOS_TRY
        // Define a dense matrix to hold the reduced problem
        rA = ZeroMatrix(mNumberOfPetrovGalerkinRomModes, this->GetNumberOfROMModes());
        rb = ZeroVector(mNumberOfPetrovGalerkinRomModes);

        // Build the system matrix by looping over elements and conditions and assembling to A
        KRATOS_ERROR_IF(!pScheme) << "No scheme provided!" << std::endl;

        // Get ProcessInfo from main model part
        const auto& r_current_process_info = rModelPart.GetProcessInfo();


        // Assemble all entities
        const auto assembling_timer = BuiltinTimer();

        using SystemSumReducer = CombinedReduction<NonTrivialSumReduction<PetrovGalerkinSystemMatrixType>, NonTrivialSumReduction<PetrovGalerkinSystemVectorType>>;
        AssemblyTLS assembly_tls_container(this->GetNumberOfROMModes(), mNumberOfPetrovGalerkinRomModes);

        const auto& r_elements = this->mHromSimulation ? this->mSelectedElements : rModelPart.Elements();

        if(!r_elements.empty())
        {
            std::tie(rA, rb) =
            block_for_each<SystemSumReducer>(r_elements, assembly_tls_container, 
                [&](Element& r_element, AssemblyTLS& r_thread_prealloc)
            {
                return CalculateLocalContributionPetrovGalerkin(r_element, rA, rb, r_thread_prealloc, *pScheme, r_current_process_info);
            });
        }


        const auto& r_conditions = this->mHromSimulation ? this->mSelectedConditions : rModelPart.Conditions();

        if(!r_conditions.empty())
        {
            PetrovGalerkinSystemMatrixType aconditions;
            PetrovGalerkinSystemVectorType bconditions;

            std::tie(aconditions, bconditions) =
            block_for_each<SystemSumReducer>(r_conditions, assembly_tls_container, 
                [&](Condition& r_condition, AssemblyTLS& r_thread_prealloc)
            {
                return CalculateLocalContributionPetrovGalerkin(r_condition, rA, rb, r_thread_prealloc, *pScheme, r_current_process_info);
            });

            rA += aconditions;
            rb += bconditions;
        }

        KRATOS_INFO_IF("PetrovGalerkinROMBuilderAndSolver", (this->GetEchoLevel() > 0)) << "Build time: " << assembling_timer.ElapsedSeconds() << std::endl;
        KRATOS_INFO_IF("PetrovGalerkinROMBuilderAndSolver", (this->GetEchoLevel() > 2)) << "Finished parallel building" << std::endl;

        // KRATOS_CATCH("")
    }

    /**
     * Solves reduced system of equations and broadcasts it
     */
    void SolveROM(
        ModelPart &rModelPart,
        PetrovGalerkinSystemMatrixType &rA,
        PetrovGalerkinSystemVectorType &rb,
        TSystemVectorType &rDx) override
    {
        KRATOS_TRY

        RomSystemVectorType dxrom(this->GetNumberOfROMModes());
        
        const auto solving_timer = BuiltinTimer();
        // Calculate the QR decomposition
        DenseHouseholderQRDecomposition<TDenseSpace> qr_decomposition;
        //                              ^Correct after properly defining PetrovGalerkinSystemMatrixType
        qr_decomposition.Compute(rA);
        qr_decomposition.Solve(rb, dxrom);
        KRATOS_INFO_IF("PetrovGalerkinROMBuilderAndSolver", (this->GetEchoLevel() > 0)) << "Solve reduced system time: " << solving_timer.ElapsedSeconds() << std::endl;

        // Save the ROM solution increment in the root modelpart database
        auto& r_root_mp = rModelPart.GetRootModelPart();
        noalias(r_root_mp.GetValue(ROM_CURRENT_SOLUTION_TOTAL      )) += dxrom;
        noalias(r_root_mp.GetValue(ROM_CURRENT_NON_LINEAR_INCREMENT))  = dxrom;

        // project reduced solution back to full order model
        const auto backward_projection_timer = BuiltinTimer();
        this->ProjectToFineBasis(dxrom, rModelPart, rDx);
        KRATOS_INFO_IF("PetrovGalerkinROMBuilderAndSolver", (this->GetEchoLevel() > 0)) << "Project to fine basis time: " << backward_projection_timer.ElapsedSeconds() << std::endl;

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
    ///@}
    ///@name Private operations 
    ///@{

    /**
     * Computes the local contribution of an element or condition for PetrovGalerkin
     */
    template<typename TEntity>
    std::tuple<LocalSystemMatrixType, LocalSystemVectorType> CalculateLocalContributionPetrovGalerkin(
        TEntity& rEntity,
        PetrovGalerkinSystemMatrixType& rAglobal,
        PetrovGalerkinSystemVectorType& rBglobal,
        AssemblyTLS& rPreAlloc,
        TSchemeType& rScheme,
        const ProcessInfo& rCurrentProcessInfo)
    {
        if (rEntity.IsDefined(ACTIVE) && rEntity.IsNot(ACTIVE))
        {
            rPreAlloc.romA = ZeroMatrix(mNumberOfPetrovGalerkinRomModes, this->GetNumberOfROMModes());
            rPreAlloc.romB = ZeroVector(mNumberOfPetrovGalerkinRomModes);
            return std::tie(rPreAlloc.romA, rPreAlloc.romB);
        }

        // Calculate elemental contribution
        rScheme.CalculateSystemContributions(rEntity, rPreAlloc.lhs, rPreAlloc.rhs, rPreAlloc.eq_id, rCurrentProcessInfo);
        rEntity.GetDofList(rPreAlloc.dofs, rCurrentProcessInfo);

        const SizeType ndofs = rPreAlloc.dofs.size();
        ResizeIfNeeded(rPreAlloc.phiE, ndofs, this->GetNumberOfROMModes());
        ResizeIfNeeded(rPreAlloc.psiE, ndofs, mNumberOfPetrovGalerkinRomModes);
        ResizeIfNeeded(rPreAlloc.aux, ndofs, this->GetNumberOfROMModes());

        const auto &r_geom = rEntity.GetGeometry();
        RomAuxiliaryUtilities::GetPhiElemental(rPreAlloc.phiE, rPreAlloc.dofs, r_geom, this->mMapPhi);
        RomAuxiliaryUtilities::GetPsiElemental(rPreAlloc.psiE, rPreAlloc.dofs, r_geom, this->mMapPhi);

        const double h_rom_weight = this->mHromSimulation ? rEntity.GetValue(HROM_WEIGHT) : 1.0;

        noalias(rPreAlloc.aux) = prod(rPreAlloc.lhs, rPreAlloc.phiE);
        noalias(rPreAlloc.romA) = prod(trans(rPreAlloc.psiE), rPreAlloc.aux) * h_rom_weight;
        noalias(rPreAlloc.romB) = prod(trans(rPreAlloc.psiE), rPreAlloc.rhs) * h_rom_weight;

        return std::tie(rPreAlloc.romA, rPreAlloc.romB);
    }


    ///@}
}; /* Class PetrovGalerkinROMBuilderAndSolver */

///@}
///@name Type Definitions
///@{


///@}

} /* namespace Kratos.*/

