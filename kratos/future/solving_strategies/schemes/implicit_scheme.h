//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//                   Riccardo Rossi
//

#pragma once

// System includes

// External includes

// Project includes
#include "containers/csr_matrix.h"
#include "containers/system_vector.h"
#include "includes/kratos_parameters.h"
#include "includes/master_slave_constraint.h"
#include "includes/model_part.h"
#include "spaces/kratos_space.h"
#include "utilities/builtin_timer.h"
#include "utilities/dof_utilities/dof_array_utilities.h"
#include "utilities/entities_utilities.h"
#include "utilities/openmp_utils.h" //TODO: SOME FILES INCLUDING scheme.h RELY ON THIS. LEAVING AS FUTURE TODO.
#include "utilities/parallel_utilities.h"
#include "utilities/timer.h"

#ifdef KRATOS_USE_FUTURE
#include "future/linear_solvers/amgcl_solver.h"
#include "future/solving_strategies/builders/builder.h"
#include "future/solving_strategies/builders/block_builder.h"
#endif

namespace Kratos::Future
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

/// TLS type definition

/**
 * @brief Implicit scheme TLS type definition
 * Thread Local Storage container to be used in the parallel assembly of implicit problems
 * @tparam DataType data type of the problem to be solved
 */
template<class TDataType = double >
struct ImplicitThreadLocalStorage
{
    // Local LHS contribution
    DenseMatrix<TDataType> LocalMatrix;

    // Local RHS constribution
    DenseVector<TDataType> LocalVector;

    // Vector containing the localization in the system of the different terms
    Element::EquationIdVectorType LocalEqIds;

    // Vector containing the slave equation ids
    MasterSlaveConstraint::EquationIdVectorType SlaveEqIds;

    // Vector containing the master equation ids
    MasterSlaveConstraint::EquationIdVectorType MasterEqIds;
};

/**
 * @class ImplicitScheme
 * @ingroup KratosCore
 * @brief This class provides the implementation of the basic tasks that are needed by the solution strategy.
 * @details It is intended to be the place for tailoring the solution strategies to problem specific tasks.
 * @author Ruben Zorrilla
 */
//TODO: Think about the template parameters

template<class TSparseMatrixType, class TSparseVectorType, class TSparseGraphType>
class ImplicitScheme
{
public:

    // FIXME: Does not work... ask @Charlie
    // /// Add scheme to Kratos registry
    // KRATOS_REGISTRY_ADD_TEMPLATE_PROTOTYPE("Schemes.KratosMultiphysics", ImplicitScheme, ImplicitScheme, TSparseMatrixType, TSparseVectorType)
    // KRATOS_REGISTRY_ADD_TEMPLATE_PROTOTYPE("Schemes.All", ImplicitScheme, ImplicitScheme, TSparseMatrixType, TSparseVectorType)

    ///@name Type Definitions
    ///@{

    /// Pointer definition of ImplicitScheme
    KRATOS_CLASS_POINTER_DEFINITION(ImplicitScheme);

    /// Size type definition
    using SizeType = std::size_t;

    /// Index type definition
    using IndexType = typename TSparseMatrixType::IndexType;

    /// Data type definition
    using DataType = typename TSparseMatrixType::DataType;

    /// TLS type
    using TLSType = ImplicitThreadLocalStorage<DataType>;

    /// Assembly helper type
    //FIXME: This should be set in the constructor
    using BuilderType = Future::BlockBuilder<TLSType, TSparseMatrixType, TSparseVectorType, TSparseGraphType>;

    /// DoF type definition
    using DofType = Dof<DataType>;

    /// DoF array type definition
    using DofsArrayType = ModelPart::DofsArrayType;

    /// Effective DOFs map type definition
    using EffectiveDofsMapType = typename BuilderType::EffectiveDofsMapType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default Constructor
     * @details Initializes the flags
     */
    explicit ImplicitScheme() = default;

    /**
     * @brief Constructor with Parameters
     */
    explicit ImplicitScheme(
        ModelPart& rModelPart,
        Parameters ThisParameters)
        : mpModelPart(&rModelPart)
    {
        // Validate default parameters
        ThisParameters = this->ValidateAndAssignParameters(ThisParameters, this->GetDefaultParameters());
        this->AssignSettings(ThisParameters);

        //TODO: User-definable reshaping stuff

        // Set up the assembly helper
        Parameters build_settings = ThisParameters["build_settings"];
        build_settings.AddInt("echo_level", ThisParameters["echo_level"].GetInt());
        mpBuilder = Kratos::make_unique<BuilderType>(rModelPart, build_settings);
    }

    /** Copy Constructor.
     */
    explicit ImplicitScheme(ImplicitScheme& rOther)
      : mSchemeIsInitialized(rOther.mSchemeIsInitialized)
      , mSchemeSolutionStepIsInitialized(rOther.mSchemeSolutionStepIsInitialized)
    {
        //TODO: Check this... particularly the mpBuilder pointer
    }

    /** Destructor.
     */
    virtual ~ImplicitScheme() = default;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Create method
     * @param ThisParameters The configuration parameters
     */
    virtual typename ImplicitScheme<TSparseMatrixType, TSparseVectorType, TSparseGraphType>::Pointer Create(
        ModelPart& rModelPart,
        Parameters ThisParameters) const
    {
        return Kratos::make_shared<ImplicitScheme<TSparseMatrixType, TSparseVectorType, TSparseGraphType>>(rModelPart, ThisParameters);
    }

    /**
     * @brief Clone method
     * @return The pointer of the cloned ImplicitScheme
     */
    virtual typename ImplicitScheme<TSparseMatrixType, TSparseVectorType, TSparseGraphType>::Pointer Clone()
    {
        return Kratos::make_shared<ImplicitScheme<TSparseMatrixType, TSparseVectorType, TSparseGraphType>>(*this) ;
    }

    /**
     * @brief This is the place to initialize the ImplicitScheme.
     * @details This is intended to be called just once when the strategy is initialized
     */
    virtual void Initialize()
    {
        KRATOS_TRY

        // Check if the Initialize has been already performed
        if (!mSchemeIsInitialized) {
            // Initialize elements, conditions and constraints
            EntitiesUtilities::InitializeAllEntities(*mpModelPart);

            // Set the flag to avoid calling this twice
            mSchemeIsInitialized = true; //TODO: Discuss with the KTC if these should remain or not
        }

        KRATOS_CATCH("")
    }

    /**
     * @brief Function called once at the beginning of each solution step
     * The basic operations to be carried out in here are the following:
     * 1) Set up the DOF array from the element and conditions DOFs (SetUpDofArray)
     * 2) Set up the system ids (i.e., the DOFs equation id) for the element and condition DOF array
     * 3) Construct the master slave constraints structure (this includes the effective DOF array, the effective DOF id map and the relation and constant arrays)
     * 4) Allocate the memory for the system arrays (note that this implies building the sparse matrix graph)
     * 5) Call the InitializeSolutionStep of all entities
     * Further operations might be required depending on the time integration scheme
     * Note that steps from 1 to 4 can be done once if the DOF set does not change (i.e., the mesh and the constraints active/inactive status do not change in time)
     * @param rDofSet The array of DOFs from elements and conditions
     * @param rEffectiveDofSet The array of DOFs to be solved after the application of constraints
     * @param rEffectiveDofIdMap A map relating each effective DOF to its effective id
     * @param rpA The system left hand side matrix
     * @param rpEffectiveLhs The effective left hand side matrix
     * @param rpB The system right hand side vector
     * @param rpEffectiveRhs The effective right hand side vector
     * @param rpDx The solution update vector
     * @param rpEffectiveDx The effective solution update vector
     * @param rConstraintsRelationMatrix The assembled constraints relation matrix (i.e. T)
     * @param rConstraintsConstantVector The assembled constraints constant vector
     * @param ReformDofSet Flag to indicate if the DOFs have changed and need to be updated
     */
    virtual void InitializeSolutionStep(
        DofsArrayType& rDofSet,
        DofsArrayType& rEffectiveDofSet,
        EffectiveDofsMapType& rEffectiveDofIdMap,
        typename TSparseMatrixType::Pointer& rpA,
        typename TSparseMatrixType::Pointer& rpEffectiveLhs,
        typename TSparseVectorType::Pointer& rpB,
        typename TSparseVectorType::Pointer& rpEffectiveRhs,
        typename TSparseVectorType::Pointer& rpDx,
        typename TSparseVectorType::Pointer& rpEffectiveDx,
        TSparseMatrixType& rConstraintsRelationMatrix,
        TSparseVectorType& rConstraintsConstantVector,
        const bool ReformDofSets = true)
    {
        KRATOS_ERROR << "\'ImplicitScheme\' does not implement \'InitializeSolutionStep\' method. Call derived class one." << std::endl;
    }

    /**
     * @brief Function called once at the end of a solution step, after convergence is reached if an iterative process is needed
     * @param rModelPart The model part of the problem to solve
     * @param A LHS matrix
     * @param Dx Incremental update of primary variables
     * @param b RHS Vector
     */
    //TODO: Think on the arguments of this one (I'd pass all in order to provide maximum flexibility in derived classes)
    virtual void FinalizeSolutionStep(
        TSparseMatrixType& A,
        TSparseVectorType& Dx,
        TSparseVectorType& b)
    {
        KRATOS_TRY

        // Finalizes solution step for all of the elements, conditions and constraints
        EntitiesUtilities::FinalizeSolutionStepAllEntities(*mpModelPart);

        // Reset flags for next step
        mSchemeSolutionStepIsInitialized = false;

        KRATOS_CATCH("")
    }

    /**
     * @brief unction to be called when it is needed to initialize an iteration. It is designed to be called at the beginning of each non linear iteration
     * @note Take care: the elemental function with the same name is NOT called here.
     * @warning Must be defined in derived classes
     * @details The function is called in the builder for memory efficiency
     * @param rModelPart The model part of the problem to solve
     * @param A LHS matrix
     * @param Dx Incremental update of primary variables
     * @param b RHS Vector
     */
    //TODO: Think on the arguments of this one (I'd pass all in order to provide maximum flexibility in derived classes)
    virtual void InitializeNonLinIteration(
        TSparseMatrixType& A,
        TSparseVectorType& Dx,
        TSparseVectorType& b)
    {
        KRATOS_TRY

        // Finalizes non-linear iteration for all of the elements, conditions and constraints
        EntitiesUtilities::InitializeNonLinearIterationAllEntities(*mpModelPart);

        KRATOS_CATCH("")
    }

    /**
     * @brief Function to be called when it is needed to finalize an iteration. It is designed to be called at the end of each non linear iteration
     * @param rModelPart The model part of the problem to solve
     * @param A LHS matrix
     * @param Dx Incremental update of primary variables
     * @param b RHS Vector
     */
    //TODO: Think on the arguments of this one (I'd pass all in order to provide maximum flexibility in derived classes)
    virtual void FinalizeNonLinIteration(
        TSparseMatrixType& A,
        TSparseVectorType& Dx,
        TSparseVectorType& b)
    {
        KRATOS_TRY

        // Finalizes non-linear iteration for all of the elements, conditions and constraints
        EntitiesUtilities::FinalizeNonLinearIterationAllEntities(*mpModelPart);

        KRATOS_CATCH("")
    }

    virtual void SetUpDofArray(DofsArrayType& rDofSet)
    {
        // Call the external utility to set up the DOFs array
        DofArrayUtilities::SetUpDofArray(*mpModelPart, rDofSet, mEchoLevel);

        // Set the corresponding flag
        mDofSetIsInitialized = true;

        KRATOS_INFO_IF("ImplicitScheme", mEchoLevel >= 2) << "Finished DOFs array set up." << std::endl;
    }

    SizeType SetUpSystemIds(DofsArrayType& rDofSet)
    {
        KRATOS_TRY

        KRATOS_ERROR_IF(rDofSet.empty()) << "DOFs set is empty. Call the 'SetUpDofArray' first." << std::endl;

        const SizeType equation_system_size = mpBuilder->SetUpSystemIds(rDofSet);

        KRATOS_INFO_IF("ImplicitScheme", mEchoLevel >= 2) << "Finished system set up." << std::endl;

        return equation_system_size;

        KRATOS_CATCH("")
    }

    //TODO: Think on the overloads for the mass and damping matrices
    virtual void ResizeAndInitializeVectors(
        const DofsArrayType& rDofSet,
        const DofsArrayType& rEffectiveDofSet,
        typename TSparseMatrixType::Pointer& rpLhs,
        typename TSparseMatrixType::Pointer& rpEffectiveLhs,
        typename TSparseVectorType::Pointer& rpRhs,
        typename TSparseVectorType::Pointer& rpEffectiveRhs,
        typename TSparseVectorType::Pointer& rpDx,
        typename TSparseVectorType::Pointer& rpEffectiveDx,
        const bool CalculateReactions = false)
    {
        KRATOS_TRY

        // Call the assembly helper to allocate and initialize the required vectors
        // Note that this also allocates the required reaction vectors (e.g., elimination build)
        (this->GetBuilder()).ResizeAndInitializeVectors(rDofSet, rEffectiveDofSet, rpLhs, rpEffectiveLhs, rpRhs, rpEffectiveRhs, rpDx, rpEffectiveDx, CalculateReactions);

        KRATOS_INFO_IF("StaticScheme", this->GetEchoLevel() >= 2) << "Finished system initialization." << std::endl;

        KRATOS_CATCH("")
    }

    virtual void ConstructMasterSlaveConstraintsStructure(
        const DofsArrayType& rDofSet,
        DofsArrayType& rEffectiveDofSet,
        EffectiveDofsMapType& rEffectiveDofIdMap,
        TSparseMatrixType& rConstraintsRelationMatrix,
        TSparseVectorType& rConstraintsConstantVector)
    {
        KRATOS_TRY

        // Call the assembly helper to set the master-slave constraints
        (this->GetBuilder()).ConstructMasterSlaveConstraintsStructure(*mpModelPart, rDofSet, rEffectiveDofSet, rEffectiveDofIdMap, rConstraintsRelationMatrix, rConstraintsConstantVector);

        KRATOS_INFO_IF("StaticScheme", this->GetEchoLevel() >= 2) << "Finished constraints initialization." << std::endl;

        KRATOS_CATCH("")
    }

    virtual void Build(
        TSparseMatrixType& rLHS,
        TSparseVectorType& rRHS)
    {
        Timer::Start("Build");

        const auto timer = BuiltinTimer();

        // const auto elem_func = [](ModelPart::ElementConstantIterator ItElem, const ProcessInfo& rProcessInfo, TLSType& rTLS){
        //     if (ItElem->Is(ACTIVE)) {
        //         // Calculate local LHS and RHS contributions
        //         ItElem->CalculateLocalSystem(rTLS.LocalMatrix, rTLS.LocalVector, rProcessInfo);

        //         // Get the positions in the global system
        //         ItElem->EquationIdVector(rTLS.LocalEqIds, rProcessInfo);

        //         // The element is active and is to be assembled
        //         return true;
        //     } else {
        //         // Clear current TLS values
        //         rTLS.LocalEqIds.clear();
        //         rTLS.LocalVector.clear();
        //         rTLS.LocalMatrix.clear();

        //         // The element is inactive and is not to be assembled
        //         return false;
        //     }
        // };

        // const auto cond_func = [](ModelPart::ConditionConstantIterator ItCond, const ProcessInfo& rProcessInfo, TLSType& rTLS){
        //     if (ItCond->Is(ACTIVE)) {
        //         // Calculate local LHS and RHS contributions
        //         ItCond->CalculateLocalSystem(rTLS.LocalMatrix, rTLS.LocalVector, rProcessInfo);

        //         // Get the positions in the global system
        //         ItCond->EquationIdVector(rTLS.LocalEqIds, rProcessInfo);

        //         // The condition is active and is to be assembled
        //         return true;
        //     } else {
        //         // Clear current TLS values
        //         rTLS.LocalEqIds.clear();
        //         rTLS.LocalVector.clear();
        //         rTLS.LocalMatrix.clear();

        //         // The condition is inactive and is not to be assembled
        //         return false;
        //     }
        // };

        // TLSType aux_tls;
        // auto& r_assembly_helper = GetBuilder();
        // r_assembly_helper.SetElementAssemblyFunction(elem_func);
        // r_assembly_helper.SetConditionAssemblyFunction(cond_func);
        // r_assembly_helper.Assemble(rLHS, rRHS, aux_tls);

        // Getting conditions and elements to be assembled
        const auto& r_elems = mpModelPart->Elements();
        const auto& r_conds = mpModelPart->Conditions();
        const auto& r_process_info = mpModelPart->GetProcessInfo();

        // Getting entities container data
        auto elems_begin = r_elems.begin();
        auto conds_begin = r_conds.begin();
        const std::size_t n_elems = r_elems.size();
        const std::size_t n_conds = r_conds.size();

        // Initialize RHS and LHS assembly
        rRHS.BeginAssemble();
        rLHS.BeginAssemble();

        // Assemble entities
        TLSType aux_tls;
        #pragma omp parallel firstprivate(n_elems, n_conds, elems_begin, conds_begin, r_process_info, aux_tls)
        {
            // Assemble elements
            # pragma omp for schedule(guided, 512) nowait
            for (int k = 0; k < n_elems; ++k) {
                // Calculate local LHS and RHS contributions
                auto it_elem = elems_begin + k;
                const bool assemble = CalculateLocalSystemContribution(*it_elem, aux_tls, r_process_info);

                // Assemble the local contributions to the global system
                if (assemble) {
                    rRHS.Assemble(aux_tls.LocalVector, aux_tls.LocalEqIds); // RHS contributions assembly FIXME: Do the ij-test
                    rLHS.Assemble(aux_tls.LocalMatrix, aux_tls.LocalEqIds); // LHS contributions assembly FIXME: Do the ij-test
                }
            }

            // Assemble conditions
            # pragma omp for schedule(guided, 512)
            for (int k = 0; k < n_conds; ++k) {
                // Calculate local LHS and RHS contributions
                auto it_cond = conds_begin + k;
                const bool assemble = CalculateLocalSystemContribution(*it_cond, aux_tls, r_process_info);

                // Assemble the local contributions to the global system
                if (assemble) {
                    rRHS.Assemble(aux_tls.LocalVector, aux_tls.LocalEqIds); // RHS contributions assembly FIXME: Do the ij-test
                    rLHS.Assemble(aux_tls.LocalMatrix, aux_tls.LocalEqIds); // LHS contributions assembly FIXME: Do the ij-test
                }
            }
        }

        // Finalize RHS and LHS assembly
        rRHS.FinalizeAssemble();
        rLHS.FinalizeAssemble();

        KRATOS_INFO_IF("ImplicitScheme", mEchoLevel >= 1) << "Build time: " << timer << std::endl;
        KRATOS_INFO_IF("ImplicitScheme", mEchoLevel >= 2) << "Finished parallel building" << std::endl;

        Timer::Stop("Build");
    }

    virtual void Build(TSparseVectorType& rRHS)
    {
        Timer::Start("BuildRightHandSide");

        const auto timer = BuiltinTimer();

        // const auto elem_func = [](ModelPart::ElementConstantIterator ItElem, const ProcessInfo& rProcessInfo, TLSType& rTLS){
        //     if (ItElem->Is(ACTIVE)) {
        //         // Calculate the RHS contribution
        //         ItElem->CalculateRightHandSide(rTLS.LocalVector, rProcessInfo);

        //         // Get the positions in the global system
        //         ItElem->EquationIdVector(rTLS.LocalEqIds, rProcessInfo);

        //         // The element is active and is to be assembled
        //         return true;
        //     } else {
        //         // Clear current TLS values
        //         rTLS.LocalEqIds.clear();
        //         rTLS.LocalVector.clear();

        //         // The element is inactive and is not to be assembled
        //         return false;
        //     }
        // };

        // const auto cond_func = [](ModelPart::ConditionConstantIterator ItCond, const ProcessInfo& rProcessInfo, TLSType& rTLS){
        //     if (ItCond->Is(ACTIVE)) {
        //         // Calculate the RHS contribution
        //         ItCond->CalculateRightHandSide(rTLS.LocalVector, rProcessInfo);

        //         // Get the positions in the global system
        //         ItCond->EquationIdVector(rTLS.LocalEqIds, rProcessInfo);

        //         // The condition is active and is to be assembled
        //         return true;
        //     } else {
        //         // Clear current TLS values
        //         rTLS.LocalEqIds.clear();
        //         rTLS.LocalVector.clear();

        //         // The condition is inactive and is not to be assembled
        //         return false;
        //     }
        // };

        // TLSType aux_tls;
        // auto& r_assembly_helper = GetBuilder();
        // r_assembly_helper.SetElementAssemblyFunction(elem_func);
        // r_assembly_helper.SetConditionAssemblyFunction(cond_func);
        // r_assembly_helper.Assemble(rRHS, aux_tls);

        // Getting conditions and elements to be assembled
        const auto& r_elems = mpModelPart->Elements();
        const auto& r_conds = mpModelPart->Conditions();
        const auto& r_process_info = mpModelPart->GetProcessInfo();

        // Getting entities container data
        auto elems_begin = r_elems.begin();
        auto conds_begin = r_conds.begin();
        const std::size_t n_elems = r_elems.size();
        const std::size_t n_conds = r_conds.size();

        // Initialize RHS assembly
        rRHS.BeginAssemble();

        // Assemble entities
        TLSType aux_tls;
        #pragma omp parallel firstprivate(n_elems, n_conds, elems_begin, conds_begin, r_process_info, aux_tls)
        {
            // Assemble elements
            # pragma omp for schedule(guided, 512) nowait
            for (int k = 0; k < n_elems; ++k) {
                // Calculate local RHS contribution
                auto it_elem = elems_begin + k;
                const bool assemble = CalculateRightHandSideContribution(*it_elem, aux_tls, r_process_info);

                // Assemble the local contributions to the global system
                if (assemble) {
                    rRHS.Assemble(aux_tls.LocalVector, aux_tls.LocalEqIds); // RHS contributions assembly FIXME: Do the ij-test
                }
            }

            // Assemble conditions
            # pragma omp for schedule(guided, 512)
            for (int k = 0; k < n_conds; ++k) {
                // Calculate local RHS contribution
                auto it_cond = conds_begin + k;
                const bool assemble = CalculateRightHandSideContribution(*it_cond, aux_tls, r_process_info);

                // Assemble the local contributions to the global system
                if (assemble) {
                    rRHS.Assemble(aux_tls.LocalVector, aux_tls.LocalEqIds); // RHS contributions assembly FIXME: Do the ij-test
                }
            }
        }

        // Finalize RHS assembly
        rRHS.FinalizeAssemble();

        KRATOS_INFO_IF("ImplicitScheme", mEchoLevel >= 1) << "Build time (RightHandSide only): " << timer << std::endl;
        KRATOS_INFO_IF("ImplicitScheme", mEchoLevel >= 2) << "Finished parallel building (RightHandSide only)" << std::endl;

        Timer::Stop("BuildRightHandSide");
    }

    virtual void Build(TSparseMatrixType& rLHS)
    {
        Timer::Start("BuildLeftHandSide");

        const auto timer = BuiltinTimer();

        // const auto elem_func = [](ModelPart::ElementConstantIterator ItElem, const ProcessInfo& rProcessInfo, TLSType& rTLS){
        //     if (ItElem->Is(ACTIVE)) {
        //         // Calculate local LHS contribution
        //         ItElem->CalculateLeftHandSide(rTLS.LocalMatrix, rProcessInfo);

        //         // Get the positions in the global system
        //         ItElem->EquationIdVector(rTLS.LocalEqIds, rProcessInfo);

        //         // The element is active and is to be assembled
        //         return true;
        //     } else {
        //         // Clear current TLS values
        //         rTLS.LocalEqIds.clear();
        //         rTLS.LocalMatrix.clear();

        //         // The element is inactive and is not to be assembled
        //         return false;
        //     }
        // };

        // const auto cond_func = [](ModelPart::ConditionConstantIterator ItCond, const ProcessInfo& rProcessInfo, TLSType& rTLS){
        //     if (ItCond->Is(ACTIVE)) {
        //         // Calculate local LHS contribution
        //         ItCond->CalculateLeftHandSide(rTLS.LocalMatrix, rProcessInfo);

        //         // Get the positions in the global system
        //         ItCond->EquationIdVector(rTLS.LocalEqIds, rProcessInfo);

        //         // The condition is active and is to be assembled
        //         return true;
        //     } else {
        //         // Clear current TLS values
        //         rTLS.LocalEqIds.clear();
        //         rTLS.LocalMatrix.clear();

        //         // The condition is inactive and is not to be assembled
        //         return false;
        //     }
        // };

        // TLSType aux_tls;
        // auto& r_assembly_helper = GetBuilder();
        // r_assembly_helper.SetElementAssemblyFunction(elem_func);
        // r_assembly_helper.SetConditionAssemblyFunction(cond_func);
        // r_assembly_helper.Assemble(rLHS, aux_tls);

        // Getting conditions and elements to be assembled
        const auto& r_elems = mpModelPart->Elements();
        const auto& r_conds = mpModelPart->Conditions();
        const auto& r_process_info = mpModelPart->GetProcessInfo();

        // Getting entities container data
        auto elems_begin = r_elems.begin();
        auto conds_begin = r_conds.begin();
        const std::size_t n_elems = r_elems.size();
        const std::size_t n_conds = r_conds.size();

        // Initialize LHS assembly
        rLHS.BeginAssemble();

        // Assemble entities
        TLSType aux_tls;
        #pragma omp parallel firstprivate(n_elems, n_conds, elems_begin, conds_begin, r_process_info, aux_tls)
        {
            // Assemble elements
            # pragma omp for schedule(guided, 512) nowait
            for (int k = 0; k < n_elems; ++k) {
                // Calculate local LHS contribution
                auto it_elem = elems_begin + k;
                const bool assemble = CalculateLeftHandSideContribution(*it_elem, aux_tls, r_process_info);

                // Assemble the local contributions to the global system
                if (assemble) {
                    rLHS.Assemble(aux_tls.LocalMatrix, aux_tls.LocalEqIds); // LHS contributions assembly FIXME: Do the ij-test
                }
            }

            // Assemble conditions
            # pragma omp for schedule(guided, 512)
            for (int k = 0; k < n_conds; ++k) {
                // Calculate local LHS contribution
                auto it_cond = conds_begin + k;
                const bool assemble = CalculateLeftHandSideContribution(*it_cond, aux_tls, r_process_info);

                // Assemble the local contributions to the global system
                if (assemble) {
                    rLHS.Assemble(aux_tls.LocalMatrix, aux_tls.LocalEqIds); // LHS contributions assembly FIXME: Do the ij-test
                }
            }
        }

        // Finalize LHS assembly
        rLHS.FinalizeAssemble();

        KRATOS_INFO_IF("ImplicitScheme", mEchoLevel >= 1) << "Build time (LeftHandSide only): " << timer << std::endl;
        KRATOS_INFO_IF("ImplicitScheme", mEchoLevel >= 2) << "Finished parallel building (LeftHandSide only)" << std::endl;

        Timer::Stop("BuildLeftHandSide");
    }

    virtual void BuildMassMatrix(TSparseMatrixType& rMassMatrix)
    {
        Timer::Start("BuildMassMatrix");

        const auto timer = BuiltinTimer();

        // const auto elem_func = [](ModelPart::ElementConstantIterator ItElem, const ProcessInfo& rProcessInfo, TLSType& rTLS){
        //     if (ItElem->Is(ACTIVE)) {
        //         // Calculate local mass matrix contribution
        //         ItElem->CalculateMassMatrix(rTLS.LocalMatrix, rProcessInfo);

        //         // Get the positions in the global system
        //         ItElem->EquationIdVector(rTLS.LocalEqIds, rProcessInfo);

        //         // The element is active and is to be assembled
        //         return true;
        //     } else {
        //         // Clear current TLS values
        //         rTLS.LocalEqIds.clear();
        //         rTLS.LocalMatrix.clear();

        //         // The element is inactive and is not to be assembled
        //         return false;
        //     }
        // };

        // const auto cond_func = [](ModelPart::ConditionConstantIterator ItCond, const ProcessInfo& rProcessInfo, TLSType& rTLS){
        //     if (ItCond->Is(ACTIVE)) {
        //         // Calculate local mass matrix contribution
        //         ItCond->CalculateMassMatrix(rTLS.LocalMatrix, rProcessInfo);

        //         // Get the positions in the global system
        //         ItCond->EquationIdVector(rTLS.LocalEqIds, rProcessInfo);

        //         // The condition is active and is to be assembled
        //         return true;
        //     } else {
        //         // Clear current TLS values
        //         rTLS.LocalEqIds.clear();
        //         rTLS.LocalMatrix.clear();

        //         // The condition is inactive and is not to be assembled
        //         return false;
        //     }
        // };

        // TLSType aux_tls;
        // auto& r_assembly_helper = GetBuilder();
        // r_assembly_helper.SetElementAssemblyFunction(elem_func);
        // r_assembly_helper.SetConditionAssemblyFunction(cond_func);
        // r_assembly_helper.Assemble(rMassMatrix, aux_tls);

        // Getting conditions and elements to be assembled
        const auto& r_elems = mpModelPart->Elements();
        const auto& r_conds = mpModelPart->Conditions();
        const auto& r_process_info = mpModelPart->GetProcessInfo();

        // Getting entities container data
        auto elems_begin = r_elems.begin();
        auto conds_begin = r_conds.begin();
        const std::size_t n_elems = r_elems.size();
        const std::size_t n_conds = r_conds.size();

        // Initialize LHS assembly
        rMassMatrix.BeginAssemble();

        // Assemble entities
        TLSType aux_tls;
        #pragma omp parallel firstprivate(n_elems, n_conds, elems_begin, conds_begin, r_process_info, aux_tls)
        {
            // Assemble elements
            # pragma omp for schedule(guided, 512) nowait
            for (int k = 0; k < n_elems; ++k) {
                // Calculate local LHS contribution
                auto it_elem = elems_begin + k;
                const bool assemble = CalculateMassMatrixContribution(*it_elem, aux_tls, r_process_info);

                // Assemble the local contributions to the global system
                if (assemble) {
                    rMassMatrix.Assemble(aux_tls.LocalMatrix, aux_tls.LocalEqIds); // LHS contributions assembly FIXME: Do the ij-test
                }
            }

            // Assemble conditions
            # pragma omp for schedule(guided, 512)
            for (int k = 0; k < n_conds; ++k) {
                // Calculate local LHS contribution
                auto it_cond = conds_begin + k;
                const bool assemble = CalculateMassMatrixContribution(*it_cond, aux_tls, r_process_info);

                // Assemble the local contributions to the global system
                if (assemble) {
                    rMassMatrix.Assemble(aux_tls.LocalMatrix, aux_tls.LocalEqIds); // LHS contributions assembly FIXME: Do the ij-test
                }
            }
        }

        // Finalize LHS assembly
        rMassMatrix.FinalizeAssemble();

        KRATOS_INFO_IF("ImplicitScheme", mEchoLevel >= 1) << "Build mass matrix time: " << timer << std::endl;
        KRATOS_INFO_IF("ImplicitScheme", mEchoLevel >= 2) << "Finished mass matrix parallel building" << std::endl;

        Timer::Stop("BuildMassMatrix");
    }

    virtual void BuildDampingMatrix(TSparseMatrixType& rDampingMatrix)
    {
        Timer::Start("BuildDampingMatrix");

        const auto timer = BuiltinTimer();

        // const auto elem_func = [](ModelPart::ElementConstantIterator ItElem, const ProcessInfo& rProcessInfo, TLSType& rTLS){
        //     if (ItElem->Is(ACTIVE)) {
        //         // Calculate local damping matrix contribution
        //         ItElem->CalculateDampingMatrix(rTLS.LocalMatrix, rProcessInfo);

        //         // Get the positions in the global system
        //         ItElem->EquationIdVector(rTLS.LocalEqIds, rProcessInfo);

        //         // The element is active and is to be assembled
        //         return true;
        //     } else {
        //         // Clear current TLS values
        //         rTLS.LocalEqIds.clear();
        //         rTLS.LocalMatrix.clear();

        //         // The element is inactive and is not to be assembled
        //         return false;
        //     }
        // };

        // const auto cond_func = [](ModelPart::ConditionConstantIterator ItCond, const ProcessInfo& rProcessInfo, TLSType& rTLS){
        //     if (ItCond->Is(ACTIVE)) {
        //         // Calculate local damping matrix contribution
        //         ItCond->CalculateDampingMatrix(rTLS.LocalMatrix, rProcessInfo);

        //         // Get the positions in the global system
        //         ItCond->EquationIdVector(rTLS.LocalEqIds, rProcessInfo);

        //         // The condition is active and is to be assembled
        //         return true;
        //     } else {
        //         // Clear current TLS values
        //         rTLS.LocalEqIds.clear();
        //         rTLS.LocalMatrix.clear();

        //         // The condition is inactive and is not to be assembled
        //         return false;
        //     }
        // };

        // TLSType aux_tls;
        // auto& r_assembly_helper = GetBuilder();
        // r_assembly_helper.SetElementAssemblyFunction(elem_func);
        // r_assembly_helper.SetConditionAssemblyFunction(cond_func);
        // r_assembly_helper.Assemble(rDampingMatrix, aux_tls);

        // Getting conditions and elements to be assembled
        const auto& r_elems = mpModelPart->Elements();
        const auto& r_conds = mpModelPart->Conditions();
        const auto& r_process_info = mpModelPart->GetProcessInfo();

        // Getting entities container data
        auto elems_begin = r_elems.begin();
        auto conds_begin = r_conds.begin();
        const std::size_t n_elems = r_elems.size();
        const std::size_t n_conds = r_conds.size();

        // Initialize LHS assembly
        rDampingMatrix.BeginAssemble();

        // Assemble entities
        TLSType aux_tls;
        #pragma omp parallel firstprivate(n_elems, n_conds, elems_begin, conds_begin, r_process_info, aux_tls)
        {
            // Assemble elements
            # pragma omp for schedule(guided, 512) nowait
            for (int k = 0; k < n_elems; ++k) {
                // Calculate local LHS contribution
                auto it_elem = elems_begin + k;
                const bool assemble = CalculateDampingMatrixContribution(*it_elem, aux_tls, r_process_info);

                // Assemble the local contributions to the global system
                if (assemble) {
                    rDampingMatrix.Assemble(aux_tls.LocalMatrix, aux_tls.LocalEqIds); // LHS contributions assembly FIXME: Do the ij-test
                }
            }

            // Assemble conditions
            # pragma omp for schedule(guided, 512)
            for (int k = 0; k < n_conds; ++k) {
                // Calculate local LHS contribution
                auto it_cond = conds_begin + k;
                const bool assemble = CalculateDampingMatrixContribution(*it_cond, aux_tls, r_process_info);

                // Assemble the local contributions to the global system
                if (assemble) {
                    rDampingMatrix.Assemble(aux_tls.LocalMatrix, aux_tls.LocalEqIds); // LHS contributions assembly FIXME: Do the ij-test
                }
            }
        }

        // Finalize LHS assembly
        rDampingMatrix.FinalizeAssemble();

        KRATOS_INFO_IF("ImplicitScheme", mEchoLevel >= 1) << "Build damping matrix time: " << timer << std::endl;
        KRATOS_INFO_IF("ImplicitScheme", mEchoLevel >= 2) << "Finished damping matrix parallel building" << std::endl;

        Timer::Stop("BuildDampingMatrix");
    }

    virtual void BuildMasterSlaveConstraints(
        const DofsArrayType& rDofSet,
        const EffectiveDofsMapType& rDofIdMap,
        TSparseMatrixType& rConstraintsRelationMatrix,
        TSparseVectorType& rConstraintsConstantVector)
    {
        if (mpModelPart->NumberOfMasterSlaveConstraints() != 0) {
            Timer::Start("BuildConstraints");

            const auto timer_constraints = BuiltinTimer();

            // const auto const_func = [](ModelPart::MasterSlaveConstraintConstantIteratorType ItConst, const ProcessInfo& rProcessInfo, TLSType& rTLS){
            //     if (ItConst->Is(ACTIVE)) {
            //         // Calculate local relation matrix and constant vector contributions
            //         ItConst->CalculateLocalSystem(rTLS.LocalMatrix, rTLS.LocalVector, rProcessInfo);

            //         // The constraint is active and is to be assembled
            //         return true;
            //     } else {
            //         // Reset the local relation matrix and constant vector contributions
            //         rTLS.LocalVector.resize(0, false);
            //         rTLS.LocalMatrix.resize(0, 0, false);

            //         // The constraint is inactive and is not to be assembled
            //         return false;
            //     }
            // };

            // TLSType aux_tls;
            // auto& r_assembly_helper = GetBuilder();
            // r_assembly_helper.SetConstraintAssemblyFunction(const_func);
            // r_assembly_helper.AssembleMasterSlaveConstraints(rDofSet, rDofIdMap, rConstraintsRelationMatrix, rConstraintsConstantVector, aux_tls);

            // Getting constraints to be assembled
            const auto& r_consts = mpModelPart->MasterSlaveConstraints();
            const auto& r_process_info = mpModelPart->GetProcessInfo();

            // Getting constraints container data
            auto consts_begin = r_consts.begin();
            const std::size_t n_consts = r_consts.size();

            // Initialize constraints arrays
            rConstraintsRelationMatrix.SetValue(0.0);
            rConstraintsConstantVector.SetValue(0.0);

            // We clear the inactive DOFs set
            std::unordered_set<IndexType> inactive_slave_dofs;

            rConstraintsRelationMatrix.BeginAssemble();
            rConstraintsConstantVector.BeginAssemble();

            TLSType aux_tls;
            #pragma omp parallel firstprivate(rDofIdMap, consts_begin, r_process_info)
            {
                // Auxiliary set to store the inactive constraints slave DOFs (required by the block build)
                std::unordered_set<IndexType> auxiliar_inactive_slave_dofs;

                // Assemble constraints
                # pragma omp for schedule(guided, 512) nowait
                for (int k = 0; k < n_consts; ++k) {
                    // Calculate local contributions
                    auto it_const = consts_begin + k;
                    const bool assemble_const = CalculateConstraintContribution(*it_const, aux_tls, r_process_info);

                    // Set the master and slave equation ids
                    // Note that the slaves follow the system equation ids while the masters use the effective map ones
                    const auto& r_slave_dofs = it_const->GetSlaveDofsVector();
                    auto& r_slave_eq_ids = aux_tls.SlaveEqIds;
                    const std::size_t n_slaves = r_slave_dofs.size();
                    if (r_slave_eq_ids.size() != n_slaves) {
                        r_slave_eq_ids.resize(n_slaves);
                    }
                    for (IndexType i_slave = 0; i_slave < n_slaves; ++i_slave) {
                        r_slave_eq_ids[i_slave] = (*(r_slave_dofs.begin() + i_slave))->EquationId();
                    }

                    const auto& r_master_dofs = it_const->GetMasterDofsVector();
                    auto& r_master_eq_ids = aux_tls.MasterEqIds;
                    const std::size_t n_masters = r_master_dofs.size();
                    if (r_master_eq_ids.size() != n_masters) {
                        r_master_eq_ids.resize(n_masters);
                    }
                    for (IndexType i_master = 0; i_master < n_masters; ++i_master) {
                        auto p_master = *(r_master_dofs.begin() + i_master);
                        auto p_master_find = rDofIdMap.find(p_master);
                        KRATOS_ERROR_IF(p_master_find == rDofIdMap.end()) << "Master DOF cannot be found in DOF ids map." << std::endl;
                        r_master_eq_ids[i_master] = p_master_find->second;
                    }

                    // Assemble the constraints local contributions to the global system
                    if (assemble_const) {
                        // Assemble relation matrix contribution
                        rConstraintsRelationMatrix.Assemble(aux_tls.LocalMatrix, r_slave_eq_ids, r_master_eq_ids);

                        // Assemble constant vector contribution
                        rConstraintsConstantVector.Assemble(aux_tls.LocalVector, r_slave_eq_ids);
                    } else {
                        auxiliar_inactive_slave_dofs.insert(r_slave_eq_ids.begin(), r_slave_eq_ids.end());
                    }
                }

                // We merge all the sets in one thread
                #pragma omp critical
                {
                    inactive_slave_dofs.insert(auxiliar_inactive_slave_dofs.begin(), auxiliar_inactive_slave_dofs.end());
                }
            }

            rConstraintsRelationMatrix.FinalizeAssemble();
            rConstraintsConstantVector.FinalizeAssemble();

            // Setting the missing effective but not constrain-related DOFs into the T and C system
            // For doing so we loop the standard DOF array (the one from elements and conditions)
            // We search for each DOF in the effective DOF ids map, if present it means its effective
            IndexPartition<IndexType>(rDofSet.size()).for_each([&](IndexType Index){
                const auto p_dof = *(rDofSet.ptr_begin() + Index);
                const auto p_dof_find = rDofIdMap.find(p_dof);
                if (p_dof_find != rDofIdMap.end()) {
                    rConstraintsConstantVector[p_dof->EquationId()] = 0.0;
                    rConstraintsRelationMatrix(p_dof->EquationId(), p_dof_find->second) = 1.0;
                }
            });

            // Setting inactive slave dofs in the T and C system
            //TODO: Can't this be parallel?
            for (auto eq_id : inactive_slave_dofs) {
                rConstraintsConstantVector[eq_id] = 0.0;
                rConstraintsRelationMatrix(eq_id, eq_id) = 1.0;
            }

            KRATOS_INFO_IF("ImplicitScheme", mEchoLevel >= 1) << "Constraints build time: " << timer_constraints << std::endl;

            Timer::Stop("BuildConstraints");
        } else {
            KRATOS_INFO_IF("ImplicitScheme", mEchoLevel >= 1) << "There are no constraints to build." << std::endl;
        }
    }

    virtual void ApplyMasterSlaveConstraints(
        typename TSparseMatrixType::Pointer& rpLhs,
        typename TSparseMatrixType::Pointer& rpEffectiveLhs,
        typename TSparseVectorType::Pointer& rpRhs,
        typename TSparseVectorType::Pointer& rpEffectiveRhs,
        typename TSparseVectorType::Pointer& rpDx,
        typename TSparseVectorType::Pointer& rpEffectiveDx,
        const TSparseMatrixType& rConstraintsRelationMatrix,
        const TSparseVectorType& rConstraintsConstantVector)
    {
        if (mpModelPart->NumberOfMasterSlaveConstraints() != 0) {
            Timer::Start("ApplyConstraints");

            const auto timer_constraints = BuiltinTimer();

            auto& r_assembly_helper = GetBuilder();
            r_assembly_helper.ApplyMasterSlaveConstraints(rpLhs, rpEffectiveLhs, rpRhs, *rpEffectiveRhs, *rpDx, *rpEffectiveDx, rConstraintsRelationMatrix, rConstraintsConstantVector);

            KRATOS_INFO_IF("ImplicitScheme", mEchoLevel >= 1) << "Constraints apply time: " << timer_constraints << std::endl;

            Timer::Stop("ApplyConstraints");
        } else {
            // If there are no constraints the effective arrays are the same as the input ones
            // Note that we avoid duplicating the memory by making the effective pointers to point to the same object
            rpEffectiveLhs = rpLhs;
            rpEffectiveRhs = rpRhs;
            rpEffectiveDx = rpDx;
        }
    }

    virtual void ApplyMasterSlaveConstraints(
        typename TSparseVectorType::Pointer& rpRhs,
        typename TSparseVectorType::Pointer& rpEffectiveRhs,
        typename TSparseVectorType::Pointer& rpDx,
        typename TSparseVectorType::Pointer& rpEffectiveDx,
        const TSparseMatrixType& rConstraintsRelationMatrix,
        const TSparseVectorType& rConstraintsConstantVector)
    {
        if (mpModelPart->NumberOfMasterSlaveConstraints() != 0) {
            Timer::Start("ApplyConstraintsRightHandSide");

            const auto timer_constraints = BuiltinTimer();

            auto& r_assembly_helper = GetBuilder();
            r_assembly_helper.ApplyMasterSlaveConstraints(rpRhs, *rpEffectiveRhs, *rpDx, *rpEffectiveDx, rConstraintsRelationMatrix, rConstraintsConstantVector);

            KRATOS_INFO_IF("ImplicitScheme", mEchoLevel >= 1) << "Constraints apply time (RightHandSide only): " << timer_constraints << std::endl;

            Timer::Stop("ApplyConstraintsRightHandSide");
        } else {
            // If there are no constraints the effective arrays are the same as the input ones
            // Note that we avoid duplicating the memory by making the effective pointers to point to the same object
            rpEffectiveRhs = rpRhs;
            rpEffectiveDx = rpDx;
        }
    }

    //TODO: Think about the dynamic case and the mass and damping matrices!!
    virtual void ApplyDirichletConditions(
        const DofsArrayType& rDofArray,
        const EffectiveDofsMapType& rDofIdMap,
        TSparseMatrixType& rLHS,
        TSparseVectorType& rRHS)
    {
        GetBuilder().ApplyDirichletConditions(rDofArray, rDofIdMap, rLHS, rRHS);
    }

    virtual void ApplyDirichletConditions(
        const DofsArrayType& rDofArray,
        const EffectiveDofsMapType& rDofIdMap,
        TSparseVectorType& rRHS)
    {
        GetBuilder().ApplyDirichletConditions(rDofArray, rDofIdMap, rRHS);
    }

    //TODO: Think about the dynamic case and the mass and damping matrices!!
    virtual void CalculateReactions(
        const DofsArrayType& rDofSet,
        TSparseVectorType& rRHS)
    {
        //TODO: To be implemented

        // const auto elem_func = [](ModelPart::ElementConstantIterator ItElem, const ProcessInfo& rProcessInfo, TLSType& rTLS){
        //     if (ItElem->Is(ACTIVE)) {
        //         // Calculate the RHS contributions
        //         ItElem->CalculateRightHandSide(rTLS.LocalVector, rProcessInfo);

        //         // Get the positions in the global system
        //         ItElem->EquationIdVector(rTLS.LocalEqIds, rProcessInfo);

        //         // The element is active and is to be assembled
        //         return true;
        //     } else {
        //         // Clear current TLS values
        //         rTLS.LocalEqIds.clear();
        //         rTLS.LocalVector.clear();

        //         // The element is inactive and is not to be assembled
        //         return false;
        //     }
        // };

        // const auto cond_func = [](ModelPart::ConditionConstantIterator ItCond, const ProcessInfo& rProcessInfo, TLSType& rTLS){
        //     if (ItCond->Is(ACTIVE)) {
        //         // Calculate the RHS contributions
        //         ItCond->CalculateRightHandSide(rTLS.LocalVector, rProcessInfo);

        //         // Get the positions in the global system
        //         ItCond->EquationIdVector(rTLS.LocalEqIds, rProcessInfo);

        //         // The condition is active and is to be assembled
        //         return true;
        //     } else {
        //         // Clear current TLS values
        //         rTLS.LocalEqIds.clear();
        //         rTLS.LocalVector.clear();

        //         // The condition is inactive and is not to be assembled
        //         return false;
        //     }
        // };

        // TLSType aux_tls;
        // auto& r_assembly_helper = GetBuilder();
        // r_assembly_helper.SetElementAssemblyFunction(elem_func);
        // r_assembly_helper.SetConditionAssemblyFunction(cond_func);
        // r_assembly_helper.CalculateReactionsRightHandSide(rDofSet, rRHS, aux_tls);
    }

    /**
     * @brief Performing the prediction of the solution.
     * @warning Must be defined in derived classes
     * @param A LHS matrix
     * @param Dx Incremental update of primary variables
     * @param b RHS Vector
     */
    virtual void Predict(
        DofsArrayType& rDofSet,
        DofsArrayType& rEffectiveDofSet,
        EffectiveDofsMapType& rEffectiveDofIdMap,
        TSparseMatrixType& rA,
        TSparseVectorType& rb,
        TSparseVectorType& rDx,
        TSparseVectorType& rEffectiveDx,
        TSparseMatrixType& rConstraintsRelationMatrix,
        TSparseVectorType& rConstraintsConstantVector)
    {
        KRATOS_ERROR << "\'ImplicitScheme\' does not implement \'Predict\' method. Call derived class one." << std::endl;
    }

    /**
     * @brief Performing the update of the solution.
     * @warning Must be defined in derived classes
     * @param rDofSet Set of all primary variables
     * @param A LHS matrix
     * @param Dx Incremental update of primary variables
     * @param b RHS Vector
     */
    virtual void Update(
        DofsArrayType& rDofSet,
        EffectiveDofsMapType& rEffectiveDofIdMap,
        TSparseMatrixType& rA,
        TSparseVectorType& rb,
        TSparseVectorType& rDx,
        const TSparseVectorType& rEffectiveDx,
        const TSparseMatrixType& rConstraintsRelationMatrix)
    {
        KRATOS_ERROR << "\'ImplicitScheme\' does not implement \'Update\' method. Call derived class one." << std::endl;
    }

    void UpdateConstraintsLooseDofs(
        const TSparseVectorType& rEffectiveDx,
        DofsArrayType& rDofSet,
        EffectiveDofsMapType& rEffectiveDofIdMap)
    {
        if (mpModelPart->NumberOfMasterSlaveConstraints() != 0) {
            // Create a temporary std::unordered_set to make the search efficient
            std::unordered_set<typename DofType::Pointer> aux_dof_set(rDofSet.GetContainer().begin(), rDofSet.GetContainer().end());

            // Find the loose constraint DOFs by comparing the standard DOF set to the effective one
            // Those nodes appearing in the effective DOFs map but not in the standard DOF set are loose DOFs
            std::vector<std::pair<typename DofType::Pointer, IndexType>> loose_dofs;
            loose_dofs.reserve(rEffectiveDofIdMap.size());
            for (auto& rEffectiveDofPair : rEffectiveDofIdMap) {
                auto p_dof = rEffectiveDofPair.first;
                if (aux_dof_set.find(p_dof) == aux_dof_set.end()) {
                    loose_dofs.push_back(rEffectiveDofPair);
                }
            }

            // Update the constraint loose DOFs with the effective solution increment values
            IndexPartition<IndexType>(loose_dofs.size()).for_each([&](IndexType Index){
                auto loose_dof_info = loose_dofs[Index];
                auto p_loose_dof = loose_dof_info.first;
                if (p_loose_dof->IsFree()) {
                    p_loose_dof->GetSolutionStepValue() += rEffectiveDx[loose_dof_info.second];
                }
            });
        }
    }

    /**
     * @brief Calculates the update vector
     * This method computes the solution update vector from the effective one
     * @param rConstraintsRelationMatrix The constraints relation matrix (i.e., T)
     * @param rEffectiveDx The effective solution update vector
     * @param rDx The solution update vector
     */
    void CalculateUpdateVector(
        const TSparseMatrixType& rConstraintsRelationMatrix,
        const TSparseVectorType& rEffectiveDx,
        TSparseVectorType& rDx)
    {
        if (mpModelPart->NumberOfMasterSlaveConstraints() != 0) {
            rDx.SetValue(0.0);
            rConstraintsRelationMatrix.SpMV(rEffectiveDx, rDx);
        } else {
            rDx = rEffectiveDx;
        }
    }

    /**
     * @brief Functions to be called to prepare the data needed for the output of results.
     * @warning Must be defined in derived classes
     * @param A LHS matrix
     * @param Dx Incremental update of primary variables
     * @param b RHS Vector
     */
    virtual void CalculateOutputData(
        TSparseMatrixType& A,
        TSparseVectorType& Dx,
        TSparseVectorType& b)
    {
        KRATOS_TRY

        //TODO: Think about creating a datavalue container in here that we can access from outside.

        KRATOS_CATCH("")
    }

    /**
     * @brief Liberate internal storage.
     * @warning Must be implemented in the derived classes
     */
    virtual void Clear()
    {
        KRATOS_TRY

        // Reset initialization flags
        mSchemeIsInitialized = false;

        mSchemeSolutionStepIsInitialized = false;

        // Clear the assembly helper
        GetBuilder().Clear();

        KRATOS_CATCH("")
    }

    /**
     * @brief This function is designed to be called once to perform all the checks needed
     * on the input provided. Checks can be "expensive" as the function is designed
     * to catch user's errors.
     * @details Checks can be "expensive" as the function is designed
     * @return 0 all OK, 1 otherwise
     */
    virtual int Check() const
    {
        KRATOS_TRY

        return 0;

        KRATOS_CATCH("");
    }

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     */
    virtual Parameters GetDefaultParameters() const
    {
        const Parameters default_parameters = Parameters(R"({
            "name" : "implicit_scheme",
            "build_settings" : {
                "build_type" : "block",
                "scaling_type" : "max_diagonal"
            },
            "echo_level" : 0,
            "move_mesh" : false
        })");

        return default_parameters;
    }

    /**
     * @brief Returns the name of the class as used in the settings (snake_case format)
     * @return The name of the class
     */
    static std::string Name()
    {
        return "implicit_scheme";
    }

    ///@}
    ///@name Access
    ///@{

    /**
     * @brief This method sets the value of mMoveMesh
     * @param MoveMesh If the flag must be set to true or false
     */
    void SetMoveMesh(const bool MoveMesh)
    {
        mMoveMesh = MoveMesh;
    }

    /**
     * @brief This method sets the value of mEchoLevel
     * @param EchoLevel The value to set
     */
    void SetEchoLevel(const int EchoLevel)
    {
        mEchoLevel = EchoLevel;
    }

    /**
     * @brief This method sets the value of mDofSetIsInitialized
     * @param DofSetIsInitialized The value to set
     */
    void SetDofSetIsInitialized(const bool DofSetIsInitialized)
    {
        mDofSetIsInitialized = DofSetIsInitialized;
    }

    /**
     * @brief This method sets the value of mSchemeIsInitialized
     * @param SchemeIsInitialized The value to set
     */
    void SetSchemeIsInitialized(const bool SchemeIsInitialized)
    {
        mSchemeIsInitialized = SchemeIsInitialized;
    }

    /**
     * @brief This method sets the value of mSchemeIsInitialized
     * @param SchemeIsInitialized The value to set
     */
    void SetSchemeSolutionStepIsInitialized(const bool SchemeSolutionStepIsInitialized)
    {
        mSchemeSolutionStepIsInitialized = SchemeSolutionStepIsInitialized;
    }


    ModelPart& GetModelPart()
    {
        return *mpModelPart;
    }

    ModelPart& GetModelPart() const
    {
        return *mpModelPart;
    }

    ///@}
    ///@name Inquiry
    ///@{

    /**
     * @brief This method returns if the mesh has to be updated
     * @return bool True if to be moved, false otherwise
     */
    int GetMoveMesh() const
    {
        return mMoveMesh;
    }

    /**
     * @brief This method returns the echo level value (verbosity level)
     * @return int Echo level value
     */
    int GetEchoLevel() const
    {
        return mEchoLevel;
    }

    /**
     * @brief This method returns if the DOF set is initialized
     * @return bool True if initialized, false otherwise
     */
    bool GetDofSetIsInitialized() const
    {
        return mDofSetIsInitialized;
    }

    /**
     * @brief This method returns if the scheme is initialized
     * @return bool True if initialized, false otherwise
     */
    bool GetSchemeIsInitialized() const
    {
        return mSchemeIsInitialized;
    }

    /**
     * @brief This method returns if the scheme is initialized
     * @return bool True if initialized, false otherwise
     */
    bool GetSchemeSolutionStepIsInitialized() const
    {
        return mSchemeSolutionStepIsInitialized;
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "ImplicitScheme";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

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

    template<class TEntityType>
    bool CalculateLocalSystemContribution(
        TEntityType& rEntity,
        TLSType& rTLS,
        const ProcessInfo &rProcessInfo)
    {
        if (rEntity.Is(ACTIVE)) {
            // Calculate local RHS contribution
            rEntity.CalculateLocalSystem(rTLS.LocalMatrix, rTLS.LocalVector, rProcessInfo);

            // Get the positions in the global system
            rEntity.EquationIdVector(rTLS.LocalEqIds, rProcessInfo);

            // The element is active and is to be assembled
            return true;
        } else {
            // Clear current TLS values
            rTLS.LocalEqIds.clear();
            rTLS.LocalVector.clear();
            rTLS.LocalMatrix.clear();

            // The element is inactive and is not to be assembled
            return false;
        }
    }

    template<class TEntityType>
    bool CalculateRightHandSideContribution(
        TEntityType& rEntity,
        TLSType& rTLS,
        const ProcessInfo &rProcessInfo)
    {
        if (rEntity.Is(ACTIVE)) {
            // Calculate local RHS contribution
            rEntity.CalculateRightHandSide(rTLS.LocalVector, rProcessInfo);

            // Get the positions in the global system
            rEntity.EquationIdVector(rTLS.LocalEqIds, rProcessInfo);

            // The element is active and is to be assembled
            return true;
        } else {
            // Clear current TLS values
            rTLS.LocalEqIds.clear();
            rTLS.LocalVector.clear();

            // The element is inactive and is not to be assembled
            return false;
        }
    }

    template<class TEntityType>
    bool CalculateLeftHandSideContribution(
        TEntityType& rEntity,
        TLSType& rTLS,
        const ProcessInfo &rProcessInfo)
    {
        if (rEntity.Is(ACTIVE)) {
            // Calculate local RHS contribution
            rEntity.CalculateLeftHandSide(rTLS.LocalMatrix, rProcessInfo);

            // Get the positions in the global system
            rEntity.EquationIdVector(rTLS.LocalEqIds, rProcessInfo);

            // The element is active and is to be assembled
            return true;
        } else {
            // Clear current TLS values
            rTLS.LocalEqIds.clear();
            rTLS.LocalMatrix.clear();

            // The element is inactive and is not to be assembled
            return false;
        }
    }

    template<class TEntityType>
    bool CalculateMassMatrixContribution(
        TEntityType& rEntity,
        TLSType& rTLS,
        const ProcessInfo &rProcessInfo)
    {
        if (rEntity.Is(ACTIVE)) {
            // Calculate local RHS contribution
            rEntity.CalculateMassMatrix(rTLS.LocalMatrix, rProcessInfo);

            // Get the positions in the global system
            rEntity.EquationIdVector(rTLS.LocalEqIds, rProcessInfo);

            // The element is active and is to be assembled
            return true;
        } else {
            // Clear current TLS values
            rTLS.LocalEqIds.clear();
            rTLS.LocalMatrix.clear();

            // The element is inactive and is not to be assembled
            return false;
        }
    }

    template<class TEntityType>
    bool CalculateDampingMatrixContribution(
        TEntityType& rEntity,
        TLSType& rTLS,
        const ProcessInfo &rProcessInfo)
    {
        if (rEntity.Is(ACTIVE)) {
            // Calculate local RHS contribution
            rEntity.CalculateDampingMatrix(rTLS.LocalMatrix, rProcessInfo);

            // Get the positions in the global system
            rEntity.EquationIdVector(rTLS.LocalEqIds, rProcessInfo);

            // The element is active and is to be assembled
            return true;
        } else {
            // Clear current TLS values
            rTLS.LocalEqIds.clear();
            rTLS.LocalMatrix.clear();

            // The element is inactive and is not to be assembled
            return false;
        }
    }

    bool CalculateConstraintContribution(
        ModelPart::MasterSlaveConstraintType& rConstraint,
        TLSType& rTLS,
        const ProcessInfo &rProcessInfo)
    {
        if (rConstraint.Is(ACTIVE)) {
            // Calculate local RHS contribution
            rConstraint.CalculateLocalSystem(rTLS.LocalMatrix, rTLS.LocalVector, rProcessInfo);

            // The constraint is active and is to be assembled
            return true;
        } else {
            // Clear current TLS values
            rTLS.LocalVector.clear();
            rTLS.LocalMatrix.clear();

            // The constraint is inactive and is not to be assembled
            return false;
        }
    }

    /**
     * @brief This function is designed to move the mesh
     * @note It considers DISPLACEMENT as the variable storing the motion. Derive it to adapt to your own strategies.
     */
    virtual void MoveMesh()
    {
        KRATOS_TRY

        KRATOS_ERROR_IF_NOT(mpModelPart->HasNodalSolutionStepVariable(DISPLACEMENT_X))
            << "It is impossible to move the mesh since the DISPLACEMENT variable is not in the ModelPart. Either use SetMoveMeshFlag(False) or add DISPLACEMENT to the list of variables" << std::endl;

        block_for_each(mpModelPart->Nodes(), [](Node& rNode){
            noalias(rNode.Coordinates()) = rNode.GetInitialPosition().Coordinates();
            noalias(rNode.Coordinates()) += rNode.FastGetSolutionStepValue(DISPLACEMENT);
        });

        KRATOS_INFO_IF("SolvingStrategy", mEchoLevel != 0) << "Mesh moved." << std::endl;

        KRATOS_CATCH("")
    }

    /**
     * @brief This method validate and assign default parameters
     * @param rParameters Parameters to be validated
     * @param DefaultParameters The default parameters
     * @return Returns validated Parameters
     */
    Parameters ValidateAndAssignParameters(
        Parameters ThisParameters,
        const Parameters DefaultParameters) const
    {
        ThisParameters.ValidateAndAssignDefaults(DefaultParameters);
        return ThisParameters;
    }

    /**
     * @brief This method assigns settings to member variables
     * @param ThisParameters Parameters that are assigned to the member variables
     */
    virtual void AssignSettings(const Parameters ThisParameters)
    {
        mMoveMesh = ThisParameters["move_mesh"].GetBool();
        mEchoLevel = ThisParameters["echo_level"].GetInt();
    }

    ///@}
    ///@name Protected  Access
    ///@{

    BuilderType& GetBuilder()
    {
        return *mpBuilder;
    }

    BuilderType& GetBuilder() const
    {
        return *mpBuilder;
    }

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

    int mEchoLevel = 0;

    bool mMoveMesh = false; /// Flag to activate the mesh motion from the DISPLACEMENT variable

    bool mDofSetIsInitialized = false; /// Flag to be used in controlling if the DOF set has been already set

    bool mSchemeIsInitialized = false; /// Flag to be used in controlling if the Scheme has been initialized or not

    bool mSchemeSolutionStepIsInitialized = false; /// Flag to be used in controlling if the Scheme solution step has been initialized or not

    ModelPart* mpModelPart = nullptr;

    typename BuilderType::UniquePointer mpBuilder = nullptr;

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
    ///@name Un accessible methods
    ///@{

    ///@}
}; // Class Scheme

} // namespace Kratos::Future.
