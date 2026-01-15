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
#include <algorithm>
#include <execution>

// Project includes
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"
#include "utilities/builtin_timer.h"
#include "utilities/dof_utilities/dof_array_utilities.h"
#include "utilities/entities_utilities.h"
#include "utilities/openmp_utils.h" //TODO: SOME FILES INCLUDING scheme.h RELY ON THIS. LEAVING AS FUTURE TODO.
#include "utilities/parallel_utilities.h"
#include "utilities/timer.h"

#ifdef KRATOS_USE_FUTURE
#include "future/solving_strategies/builders/builder.h"
#include "future/solving_strategies/builders/block_builder.h"
#include "future/solving_strategies/builders/elimination_builder.h"
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

//FIXME: We need to find a way to redefine the TLS in derived classes

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
 * @brief This class provides the implementation of the basic tasks that are needed by all implicit schemes
 * @author Ruben Zorrilla
 */
//TODO: Think about the template parameters
template<class TLinearAlgebra>
class ImplicitScheme
{
public:

    // FIXME: Does not work... ask @Charlie
    // /// Add scheme to Kratos registry
    // KRATOS_REGISTRY_ADD_TEMPLATE_PROTOTYPE("Schemes.KratosMultiphysics", ImplicitScheme, ImplicitScheme, MatrixType, VectorType, TSparseGraphType)
    // KRATOS_REGISTRY_ADD_TEMPLATE_PROTOTYPE("Schemes.All", ImplicitScheme, ImplicitScheme, MatrixType, VectorType, TSparseGraphType)

    ///@name Type Definitions
    ///@{

    /// Pointer definition of ImplicitScheme
    KRATOS_CLASS_POINTER_DEFINITION(ImplicitScheme);

    /// Index type definition
    using IndexType = typename TLinearAlgebra::IndexType;

    /// Data type definition
    using DataType = typename TLinearAlgebra::DataType;

    /// Matrix type definition
    using MatrixType = typename TLinearAlgebra::MatrixType;

    /// Vector type definition
    using VectorType = typename TLinearAlgebra::VectorType;

    /// TLS type
    using TLSType = ImplicitThreadLocalStorage<DataType>;

    /// Block builder type
    using BuilderType = Future::Builder<TLinearAlgebra>;

    /// Block builder type
    using BlockBuilderType = Future::BlockBuilder<TLinearAlgebra>;

    /// Elimination builder type
    using EliminationBuilderType = Future::EliminationBuilder<TLinearAlgebra>;

    /// DoF type definition
    using DofType = Dof<DataType>;

    /// DoF array type definition
    using DofsArrayType = ModelPart::DofsArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// @brief Default constructor
    explicit ImplicitScheme() = default;

    /// @brief Constructor with parameters
    /// @param rModelPart Reference to the model part
    /// @param ThisParameters Parameters object encapsulating the settings
    explicit ImplicitScheme(
        ModelPart& rModelPart,
        Parameters ThisParameters)
        : mpModelPart(&rModelPart)
    {
        // Validate default parameters
        ThisParameters = this->ValidateAndAssignParameters(ThisParameters, this->GetDefaultParameters());
        this->AssignSettings(ThisParameters);

        // Set up the assembly helper
        Parameters build_settings = ThisParameters["build_settings"];
        build_settings.AddInt("echo_level", ThisParameters["echo_level"].GetInt());
        if (build_settings.Has("name")) {
            const std::string builder_type = build_settings["name"].GetString();
            if (builder_type == "block_builder") {
                mpBuilder = Kratos::make_unique<BlockBuilderType>(rModelPart, build_settings); // TODO: Use the registry in here
            } else if (builder_type == "elimination_builder") {
                mpBuilder = Kratos::make_unique<EliminationBuilderType>(rModelPart, build_settings); // TODO: Use the registry in here
            } else {
                KRATOS_ERROR << "Wrong builder type \'" << builder_type << "\'. Available options are \'block_builder\' and \'elimination_builder\'." << std::endl;
            }
        } else {
            KRATOS_WARNING("ImplicitScheme") << "Builder type not provided. Defaulting to \'block_builder\'." << std::endl;
            mpBuilder = Kratos::make_unique<BlockBuilderType>(rModelPart, build_settings); // TODO: Use the registry in here
        }
    }

    /// @brief Copy constructor
    /// @param rOther Other ImplicitScheme
    explicit ImplicitScheme(ImplicitScheme& rOther)
    {
        //TODO: Check this... particularly the mpBuilder pointer
    }

    /// @brief Destructor
    virtual ~ImplicitScheme() = default;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Create method
     * @param rModelPart Reference to the model part
     * @param ThisParameters The configuration parameters
     */
    virtual typename ImplicitScheme<TLinearAlgebra>::Pointer Create(
        ModelPart& rModelPart,
        Parameters ThisParameters) const
    {
        return Kratos::make_shared<ImplicitScheme<TLinearAlgebra>>(rModelPart, ThisParameters);
    }

    /**
     * @brief Clone method
     * @return The pointer of the cloned ImplicitScheme
     */
    virtual typename ImplicitScheme<TLinearAlgebra>::Pointer Clone()
    {
        return Kratos::make_shared<ImplicitScheme<TLinearAlgebra>>(*this) ;
    }

    /**
     * @brief This is the place to initialize the ImplicitScheme.
     * @details This is intended to be called just once when the strategy is initialized
     * This method sets up the linear system of equations (DOF sets and allocation) and calls the Initialize of all entities
     * Further operations might be required depending on the time integration scheme
     * Note that steps from 1 to 4 can be done once if the DOF set does not change (i.e., the mesh and the constraints active/inactive status do not change in time)
     * @param rImplicitStrategyDataContainer Auxiliary container with the linear system arrays
     */
    virtual void Initialize(ImplicitStrategyDataContainer<TLinearAlgebra> &rImplicitStrategyDataContainer)
    {
        KRATOS_TRY

        // Set up the system
        InitializeLinearSystem(rImplicitStrategyDataContainer);

        // Initialize elements, conditions and constraints
        EntitiesUtilities::InitializeAllEntities(*mpModelPart);

        KRATOS_CATCH("")
    }

    /**
     * @brief Function called once at the beginning of each solution step
     * @param rImplicitStrategyDataContainer Auxiliary container with the linear system arrays
     */
    virtual void InitializeSolutionStep(ImplicitStrategyDataContainer<TLinearAlgebra>& rImplicitStrategyDataContainer)
    {
        // Initializes solution step for all of the elements, conditions and constraints
        EntitiesUtilities::InitializeSolutionStepAllEntities(*mpModelPart);
    }

    /**
     * @brief Performing the prediction of the solution.
     * @warning Must be defined in derived classes
     * @param rImplicitStrategyDataContainer Auxiliary container with the linear system arrays
     */
    virtual void Predict(ImplicitStrategyDataContainer<TLinearAlgebra> &rImplicitStrategyDataContainer)
    {
        KRATOS_ERROR << "\'ImplicitScheme\' does not implement \'Predict\' method. Call derived class one." << std::endl;
    }

    /**
     * @brief Function called once at the end of a solution step, after convergence is reached if an iterative process is needed
     * @param rImplicitStrategyDataContainer Auxiliary container with the linear system arrays
     */
    //TODO: Think on the arguments of this one (I'd pass all in order to provide maximum flexibility in derived classes)
    virtual void FinalizeSolutionStep(ImplicitStrategyDataContainer<TLinearAlgebra> &rImplicitStrategyDataContainer)
    {
        KRATOS_TRY

        // Finalizes solution step for all of the elements, conditions and constraints
        EntitiesUtilities::FinalizeSolutionStepAllEntities(*mpModelPart);

        KRATOS_CATCH("")
    }

    /**
     * @brief Function to be called when it is needed to initialize an iteration. It is designed to be called at the beginning of each non linear iteration
     * @param rImplicitStrategyDataContainer Auxiliary container with the linear system arrays
     */
    virtual void InitializeNonLinIteration(ImplicitStrategyDataContainer<TLinearAlgebra> &rImplicitStrategyDataContainer)
    {
        KRATOS_TRY

        // Finalizes non-linear iteration for all of the elements, conditions and constraints
        EntitiesUtilities::InitializeNonLinearIterationAllEntities(*mpModelPart);

        KRATOS_CATCH("")
    }

    /**
     * @brief Function to be called when it is needed to finalize an iteration. It is designed to be called at the end of each non linear iteration
     * @param rImplicitStrategyDataContainer Auxiliary container with the linear system arrays
     */
    virtual void FinalizeNonLinIteration(ImplicitStrategyDataContainer<TLinearAlgebra> &rImplicitStrategyDataContainer)
    {
        KRATOS_TRY

        // Finalizes non-linear iteration for all of the elements, conditions and constraints
        EntitiesUtilities::FinalizeNonLinearIterationAllEntities(*mpModelPart);

        KRATOS_CATCH("")
    }

    /**
     * @brief Set the Up Dof Arrays
     * This method sets the standard and effective DOF sets
     * @param pDofSet Pointer to the standard DOF set
     * @param pEffectiveDofSet Pointer to the effective DOF set
     * @param rSlaveToMasterDofsMap The map containing the corresponding master(s) for each slave DOF
     * @return std::pair<std::size_t, std::size_t> Sizes of the standard and effective DOF sets
     */
    virtual std::pair<std::size_t, std::size_t> SetUpDofArrays(
        typename DofsArrayType::Pointer pDofSet,
        typename DofsArrayType::Pointer pEffectiveDofSet)
    {
        // Call the external utility to set up the DOFs array
        DofArrayUtilities::SetUpDofArray(*mpModelPart, *pDofSet, mEchoLevel);
        KRATOS_INFO_IF("ImplicitScheme", mEchoLevel >= 2) << "Finished DOFs array set up." << std::endl;

        // Call the external utility to set up the DOFs array
        DofArrayUtilities::SetUpEffectiveDofArray(*mpModelPart, *pDofSet, *pEffectiveDofSet, mEchoLevel);
        KRATOS_INFO_IF("ImplicitScheme", mEchoLevel >= 2) << "Finished effective DOFs array set up." << std::endl;

        // Return the sizes of the two DOF sets
        return std::make_pair(pDofSet->size(), pEffectiveDofSet->size());
    }

    /**
     * @brief Set the Up System Ids
     * This method sets the standard and effective DOF ids
     * @param pDofSet Pointer to the standard DOF set
     * @param pEffectiveDofSet Pointer to the effective DOF set
     */
    virtual void SetUpSystemIds(
        typename DofsArrayType::Pointer pDofSet,
        typename DofsArrayType::Pointer pEffectiveDofSet)
    {
        KRATOS_TRY

        // Check if the provided DOF arrays have been already set
        KRATOS_ERROR_IF(pDofSet->empty()) << "DOFs set is empty. Call the 'SetUpDofArray' first." << std::endl;
        KRATOS_ERROR_IF(pEffectiveDofSet->empty()) << "Effective DOFs set is empty. Call the 'SetUpDofArray' first." << std::endl;

        // Call the external utility to set up the DOF ids
        DofArrayUtilities::SetDofEquationIds(*pDofSet);
        KRATOS_INFO_IF("ImplicitScheme", mEchoLevel >= 2) << "Finished setting the DOF ids." << std::endl;

        // Call the external utility to set up the effective DOF ids
        DofArrayUtilities::SetEffectiveDofEquationIds(*pDofSet, *pEffectiveDofSet);
        KRATOS_INFO_IF("ImplicitScheme", mEchoLevel >= 2) << "Finished setting the effective DOF ids." << std::endl;

        KRATOS_CATCH("")
    }

    virtual void Build(
        MatrixType& rLHS,
        VectorType& rRHS)
    {
        Timer::Start("Build");

        const auto timer = BuiltinTimer();

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
        #pragma omp parallel
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

    virtual void BuildWithSafeAssemble(
        MatrixType& rLHS,
        VectorType& rRHS)
    {
        Timer::Start("BuildWithSafeAssemble");

        const auto timer = BuiltinTimer();

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
        #pragma omp parallel
        {
            // Assemble elements
            # pragma omp for schedule(guided, 512) nowait
            for (int k = 0; k < n_elems; ++k) {
                // Calculate local LHS and RHS contributions
                auto it_elem = elems_begin + k;
                const bool assemble = CalculateLocalSystemContribution(*it_elem, aux_tls, r_process_info);

                // Assemble the local contributions to the global system
                if (assemble) {
                    rRHS.SafeAssemble(aux_tls.LocalVector, aux_tls.LocalEqIds); // RHS contributions assembly FIXME: Do the ij-test
                    rLHS.SafeAssemble(aux_tls.LocalMatrix, aux_tls.LocalEqIds); // LHS contributions assembly FIXME: Do the ij-test
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
                    rRHS.SafeAssemble(aux_tls.LocalVector, aux_tls.LocalEqIds); // RHS contributions assembly FIXME: Do the ij-test
                    rLHS.SafeAssemble(aux_tls.LocalMatrix, aux_tls.LocalEqIds); // LHS contributions assembly FIXME: Do the ij-test
                }
            }
        }

        // Finalize RHS and LHS assembly
        rRHS.FinalizeAssemble();
        rLHS.FinalizeAssemble();

        KRATOS_INFO_IF("ImplicitScheme", mEchoLevel >= 1) << "Build w/ safe assemble time: " << timer << std::endl;
        KRATOS_INFO_IF("ImplicitScheme", mEchoLevel >= 2) << "Finished parallel building" << std::endl;

        Timer::Stop("BuildWithSafeAssemble");
    }

#ifdef KRATOS_USE_TBB
    virtual void BuildWithThreadLocal(
        MatrixType& rLHS,
        VectorType& rRHS)
    {
        Timer::Start("BuildWithThreadLocal");

        const auto timer = BuiltinTimer();

        // Getting conditions and elements to be assembled
        const auto& r_process_info = mpModelPart->GetProcessInfo();

        // Initialize RHS and LHS assembly
        rRHS.BeginAssemble();
        rLHS.BeginAssemble();

        thread_local Matrix LocalMatrix;
        thread_local Vector LocalVector;
        thread_local Element::EquationIdVectorType LocalEqIds;

        const auto& r_elems = mpModelPart->Elements();
        auto elems_begin = std::begin(r_elems.GetContainer());
        auto elems_end = std::end(r_elems.GetContainer());
        // auto elems_begin = r_elems.begin();
        // auto elems_end = r_elems.end();

        std::for_each(std::execution::par_unseq, elems_begin, elems_end,
            [&](auto& it_elem) {
                // const bool assemble = CalculateLocalSystemContribution(*it_elem, aux_tls, r_process_info);

                if (it_elem->IsActive()) {
                    // Calculate local RHS contribution
                    it_elem->CalculateLocalSystem(LocalMatrix, LocalVector, r_process_info);

                    // Get the positions in the global system
                    it_elem->EquationIdVector(LocalEqIds, r_process_info);

                    rRHS.Assemble(LocalVector, LocalEqIds);
                    // rRHS.SafeAssemble(LocalVector, LocalEqIds);
                    rLHS.Assemble(LocalMatrix, LocalEqIds);
                    // rLHS.SafeAssemble(LocalMatrix, LocalEqIds);
                } else {
                    LocalEqIds.clear();
                    LocalVector.clear();
                    LocalMatrix.clear();
                }
            }
        );

        const auto& r_conds = mpModelPart->Conditions();
        auto conds_begin = std::begin(r_conds.GetContainer());
        auto conds_end = std::end(r_conds.GetContainer());
        std::for_each(std::execution::par_unseq, conds_begin, conds_end,
            [&](auto& it_cond) {
                // const bool assemble = CalculateLocalSystemContribution(*it_cond, aux_tls, r_process_info);

                if (it_cond->IsActive()) {
                    // Calculate local RHS contribution
                    it_cond->CalculateLocalSystem(LocalMatrix, LocalVector, r_process_info);

                    // Get the positions in the global system
                    it_cond->EquationIdVector(LocalEqIds, r_process_info);

                    // rRHS.SafeAssemble(LocalVector, LocalEqIds);
                    rRHS.Assemble(LocalVector, LocalEqIds);
                    // rLHS.SafeAssemble(LocalMatrix, LocalEqIds);
                    rLHS.Assemble(LocalMatrix, LocalEqIds);
                } else {
                    LocalEqIds.clear();
                    LocalVector.clear();
                    LocalMatrix.clear();
                }
            }
        );

        // Finalize RHS and LHS assembly
        rRHS.FinalizeAssemble();
        rLHS.FinalizeAssemble();

        KRATOS_INFO_IF("ImplicitScheme", mEchoLevel >= 1) << "Build w/ thread local time: " << timer << std::endl;
        KRATOS_INFO_IF("ImplicitScheme", mEchoLevel >= 2) << "Finished parallel building" << std::endl;

        Timer::Stop("BuildWithThreadLocal");
    }

    virtual void BuildWithLocalAllocation(
        MatrixType& rLHS,
        VectorType& rRHS)
    {
        Timer::Start("BuildWithLocalAllocation");

        const auto timer = BuiltinTimer();

        // Getting conditions and elements to be assembled
        const auto& r_process_info = mpModelPart->GetProcessInfo();

        // Initialize RHS and LHS assembly
        rRHS.BeginAssemble();
        rLHS.BeginAssemble();

        const auto& r_elems = mpModelPart->Elements();
        auto elems_begin = std::begin(r_elems.GetContainer());
        auto elems_end = std::end(r_elems.GetContainer());
        // auto elems_begin = r_elems.begin();
        // auto elems_end = r_elems.end();

        std::for_each(std::execution::par_unseq, elems_begin, elems_end,
            [&](auto& it_elem) {
                // const bool assemble = CalculateLocalSystemContribution(*it_elem, aux_tls, r_process_info);

                if (it_elem->IsActive()) {

                    Matrix LocalMatrix;
                    Vector LocalVector;
                    Element::EquationIdVectorType LocalEqIds;

                    // Calculate local RHS contribution
                    it_elem->CalculateLocalSystem(LocalMatrix, LocalVector, r_process_info);

                    // Get the positions in the global system
                    it_elem->EquationIdVector(LocalEqIds, r_process_info);

                    rRHS.Assemble(LocalVector, LocalEqIds);
                    // rRHS.SafeAssemble(LocalVector, LocalEqIds);
                    rLHS.Assemble(LocalMatrix, LocalEqIds);
                    // rLHS.SafeAssemble(LocalMatrix, LocalEqIds);
                }
            }
        );

        const auto& r_conds = mpModelPart->Conditions();
        auto conds_begin = std::begin(r_conds.GetContainer());
        auto conds_end = std::end(r_conds.GetContainer());
        std::for_each(std::execution::par_unseq, conds_begin, conds_end,
            [&](auto& it_cond) {
                // const bool assemble = CalculateLocalSystemContribution(*it_cond, aux_tls, r_process_info);

                if (it_cond->IsActive()) {
                    Matrix LocalMatrix;
                    Vector LocalVector;
                    Element::EquationIdVectorType LocalEqIds;

                    // Calculate local RHS contribution
                    it_cond->CalculateLocalSystem(LocalMatrix, LocalVector, r_process_info);

                    // Get the positions in the global system
                    it_cond->EquationIdVector(LocalEqIds, r_process_info);

                    // rRHS.SafeAssemble(LocalVector, LocalEqIds);
                    rRHS.Assemble(LocalVector, LocalEqIds);
                    // rLHS.SafeAssemble(LocalMatrix, LocalEqIds);
                    rLHS.Assemble(LocalMatrix, LocalEqIds);
                }
            }
        );

        // Finalize RHS and LHS assembly
        rRHS.FinalizeAssemble();
        rLHS.FinalizeAssemble();

        KRATOS_INFO_IF("ImplicitScheme", mEchoLevel >= 1) << "Build w/ local allocation time: " << timer << std::endl;
        KRATOS_INFO_IF("ImplicitScheme", mEchoLevel >= 2) << "Finished parallel building" << std::endl;

        Timer::Stop("BuildWithLocalAllocation");
    }
#endif

    virtual void Build(VectorType& rRHS)
    {
        Timer::Start("BuildRightHandSide");

        const auto timer = BuiltinTimer();

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
        #pragma omp parallel
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

    virtual void Build(MatrixType& rLHS)
    {
        Timer::Start("BuildLeftHandSide");

        const auto timer = BuiltinTimer();

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
        #pragma omp parallel
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

    virtual void BuildMassMatrix(MatrixType& rMassMatrix)
    {
        Timer::Start("BuildMassMatrix");

        const auto timer = BuiltinTimer();

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
        #pragma omp parallel
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

    virtual void BuildDampingMatrix(MatrixType& rDampingMatrix)
    {
        Timer::Start("BuildDampingMatrix");

        const auto timer = BuiltinTimer();

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
        #pragma omp parallel
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

    virtual void BuildMasterSlaveConstraints(ImplicitStrategyDataContainer<TLinearAlgebra>& rImplicitStrategyDataContainer)
    {
        if (mpModelPart->NumberOfMasterSlaveConstraints() != 0) {
            Timer::Start("BuildConstraints");

            const auto timer_constraints = BuiltinTimer();

            // Getting constraints to be assembled
            const auto& r_consts = mpModelPart->MasterSlaveConstraints();
            const auto& r_process_info = mpModelPart->GetProcessInfo();

            // Getting constraints container data
            auto consts_begin = r_consts.begin();
            const std::size_t n_consts = r_consts.size();

            // Get constraints arrays from the linear system container
            auto& r_constraints_T = *(rImplicitStrategyDataContainer.pConstraintsT);
            auto& r_constraints_q = *(rImplicitStrategyDataContainer.pConstraintsQ);

            // Initialize constraints arrays
            r_constraints_T.SetValue(0.0);
            r_constraints_q.SetValue(0.0);

            // We clear the inactive DOFs set
            std::unordered_set<IndexType> inactive_slave_dofs;

            r_constraints_T.BeginAssemble();
            r_constraints_q.BeginAssemble();

            TLSType aux_tls;
            auto& r_eff_dof_set = *(rImplicitStrategyDataContainer.pEffectiveDofSet);
            #pragma omp parallel
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
                        auto p_master_find = r_eff_dof_set.find(*p_master);
                        KRATOS_ERROR_IF(p_master_find == r_eff_dof_set.end()) << "Master DOF cannot be found in effective DOF set." << std::endl;
                        r_master_eq_ids[i_master] = p_master_find->EffectiveEquationId();
                    }

                    // Assemble the constraints local contributions to the global system
                    if (assemble_const) {
                        // Assemble relation matrix contribution
                        r_constraints_T.Assemble(aux_tls.LocalMatrix, r_slave_eq_ids, r_master_eq_ids);

                        // Assemble constant vector contribution
                        r_constraints_q.Assemble(aux_tls.LocalVector, r_slave_eq_ids);
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

            r_constraints_T.FinalizeAssemble();
            r_constraints_q.FinalizeAssemble();

            // Setting the missing effective but not constrain-related DOFs into the T and C system
            // For doing so we loop the standard DOF array (the one from elements and conditions)
            // We search for each DOF in the effective DOF ids map, if present it means its effective
            auto& r_dof_set = *(rImplicitStrategyDataContainer.pDofSet);
            IndexPartition<IndexType>(r_dof_set.size()).for_each([&](IndexType Index){
                const auto p_dof = *(r_dof_set.ptr_begin() + Index);
                const auto p_dof_find = r_eff_dof_set.find(*p_dof);
                if (p_dof_find != r_eff_dof_set.end()) {
                    r_constraints_q[p_dof->EquationId()] = 0.0;
                    r_constraints_T(p_dof->EquationId(), p_dof_find->EffectiveEquationId()) = 1.0;
                }
            });

            // Setting inactive slave dofs in the T and C system
            //TODO: Can't this be parallel?
            for (auto eq_id : inactive_slave_dofs) {
                r_constraints_q[eq_id] = 0.0;
                r_constraints_T(eq_id, eq_id) = 1.0;
            }

            KRATOS_INFO_IF("ImplicitScheme", mEchoLevel >= 1) << "Constraints build time: " << timer_constraints << std::endl;

            Timer::Stop("BuildConstraints");
        } else {
            KRATOS_INFO_IF("ImplicitScheme", mEchoLevel >= 1) << "There are no constraints to build." << std::endl;
        }
    }

    /**
     * @brief Builds the linear system constraints
     * This method builds the linear system constraints, that is, the master-slave and eventual Dirichlet constraints
     * The master-slave constraints are build according to the scheme implementation while the Dirichlet ones depend on the build type
     * @param rImplicitStrategyDataContainer Auxiliary container with the linear system arrays
     */
    virtual void BuildLinearSystemConstraints(ImplicitStrategyDataContainer<TLinearAlgebra> &rImplicitStrategyDataContainer)
    {
        BuiltinTimer build_linear_system_constraints_time;

        // Build the master-slave constraints relation matrix and constant vector
        BuildMasterSlaveConstraints(rImplicitStrategyDataContainer);

        KRATOS_INFO_IF("ImplicitScheme", this->GetEchoLevel() > 0) << "Build linear system constraints time: " << build_linear_system_constraints_time << std::endl;
    }

    /**
     * @brief Applies the linear system constraints
     * This method applies the linear system constraints, that is the master-slave and the Dirichlet constraints
     * @param rEffectiveDofSet The effective DOFs array (i.e., those that are not slaves)
     * @param rImplicitStrategyDataContainer Auxiliary container with the linear system arrays
     */
    virtual void ApplyLinearSystemConstraints(ImplicitStrategyDataContainer<TLinearAlgebra> &rImplicitStrategyDataContainer)
    {
        BuiltinTimer apply_linear_system_constraints_time;

        GetBuilder().ApplyLinearSystemConstraints(rImplicitStrategyDataContainer);

        KRATOS_INFO_IF("ImplicitScheme", this->GetEchoLevel() > 0) << "Apply linear system constraints time: " << apply_linear_system_constraints_time << std::endl;
    }

    //TODO: Think about the dynamic case and the mass and damping matrices!!
    virtual void CalculateReactions(
        const DofsArrayType& rDofSet,
        VectorType& rRHS)
    {
        //TODO: To be implemented
        KRATOS_ERROR << "Not implemented yet." << std::endl;

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
     * @brief Performs the update of the solution
     * @warning Must be defined in derived classes
     * @param rDofSet Set of all primary variables
     * @param A LHS matrix
     * @param Dx Incremental update of primary variables
     * @param b RHS Vector
     */
    virtual void Update(ImplicitStrategyDataContainer<TLinearAlgebra>& rImplicitStrategyDataContainer)
    {
        KRATOS_ERROR << "\'ImplicitScheme\' does not implement \'Update\' method. Call derived class one." << std::endl;
    }

    /**
     * @brief Performs the update of the constraints only DOFs
     * This function performs the update of the constraints only DOFs (i.e., those that are not associated to elements/conditions)
     * @param rEffectiveDx Effective solution update vector
     * @param rDofSet The array of DOFs from elements and conditions
     * @param rEffectiveDofSet The effective DOFs array (i.e., those that are not slaves)
     */
    void UpdateConstraintsOnlyDofs(
        const VectorType& rEffectiveDx,
        DofsArrayType& rDofSet,
        DofsArrayType& rEffectiveDofSet)
    {
        if (mpModelPart->NumberOfMasterSlaveConstraints() != 0) {
            // Create a temporary std::unordered_set to make the search efficient
            std::unordered_set<typename DofType::Pointer> aux_dof_set(rDofSet.GetContainer().begin(), rDofSet.GetContainer().end());

            // Find the constraints only DOFs by comparing the standard DOF set to the effective one
            // Those nodes appearing in the effective DOFs set but not in the standard one are constraints only DOFs
            std::vector<typename DofType::Pointer> constr_only_dofs;
            constr_only_dofs.reserve(rEffectiveDofSet.size());
            for (unsigned int i = 0; i < rEffectiveDofSet.size(); ++i) {
                auto p_eff_dof = rEffectiveDofSet.ptr_begin() + i;
                if (aux_dof_set.find(*p_eff_dof) == aux_dof_set.end()) {
                    constr_only_dofs.push_back(*p_eff_dof);
                }
            }

            // Update the constraints only DOFs with the effective solution increment values
            IndexPartition<IndexType>(constr_only_dofs.size()).for_each([&](IndexType Index){
                auto p_constr_only_dof = constr_only_dofs[Index];
                if (p_constr_only_dof->IsFree()) {
                    p_constr_only_dof->GetSolutionStepValue() += rEffectiveDx[p_constr_only_dof->EffectiveEquationId()];
                }
            });
        }
    }

    /**
     * @brief Calculates the update vector
     * This method computes the solution update vector from the effective one
     * @param rImplicitStrategyDataContainer Auxiliary container with the linear system arrays
     */
    void CalculateUpdateVector(ImplicitStrategyDataContainer<TLinearAlgebra>& rImplicitStrategyDataContainer)
    {
        // Check if the effective relation matrix is set
        auto p_eff_T = rImplicitStrategyDataContainer.pEffectiveT;
        KRATOS_ERROR_IF(mpModelPart->NumberOfMasterSlaveConstraints() != 0 && p_eff_T == nullptr) <<
            "There are constraints but effective relation matrix is not set. Solution update vector cannot be computed." << std::endl;

        // Compute the solution vector from the effective one
        auto& r_dx = *rImplicitStrategyDataContainer.pDx;
        auto& r_eff_dx = *rImplicitStrategyDataContainer.pEffectiveDx;
        if (p_eff_T != nullptr) {
            r_dx.SetValue(0.0);
            p_eff_T->SpMV(r_eff_dx, r_dx);
        } else {
            r_dx = r_eff_dx;
        }
    }

    /**
     * @brief Liberate internal storage
     */
    virtual void Clear()
    {
        KRATOS_TRY

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
                "name" : "block_builder"
            },
            "echo_level" : 0,
            "move_mesh" : false,
            "reform_dofs_at_each_step" : false
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
     * @brief This method sets the value of mReformDofsAtEachStep
     * @param ReformDofsAtEachStep If the flag must be set to true or false
     */
    void SetReformDofsAtEachStep(const bool ReformDofsAtEachStep)
    {
        mReformDofsAtEachStep = ReformDofsAtEachStep;
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
     * @brief Get the Model Part object
     * Returns a reference to the model part the scheme is referring to
     * @return ModelPart& Reference to the scheme model part
     */
    ModelPart& GetModelPart()
    {
        return *mpModelPart;
    }

    /**
     * @brief Get the Model Part object
     * Returns a reference to the model part the scheme is referring to
     * @return const ModelPart& Reference to the scheme model part
     */
    const ModelPart& GetModelPart() const
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
    bool GetMoveMesh() const
    {
        return mMoveMesh;
    }

    /**
     * @brief This method returns if DOF sets have to be updated at each time step
     * @return bool True if to be updated, false otherwise
     */
    bool GetReformDofsAtEachStep() const
    {
        return mReformDofsAtEachStep;
    }

    /**
     * @brief This method returns the echo level value (verbosity level)
     * @return int Echo level value
     */
    int GetEchoLevel() const
    {
        return mEchoLevel;
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

    /**
     * @brief Auxiliary function to set up the implicit linear system of equations
     * The basic operations to be carried out in here are the following:
     * 1) Set up the DOF arrays from the element and conditions DOFs and the corresponding effective ones accounting for the constraints
     * 2) Set up the system ids (i.e., the DOFs equation ids), including the effective DOF ids, which may not match the "standard" ones
     * 3) Allocate the memory for the linear system constraints arrays (note that the operations done in here may depend on the build type)
     * 4) Allocate the memory for the system arrays (note that this implies building the sparse matrix graph)
     * @param rImplicitStrategyDataContainer Auxiliary container with the linear system arrays
     */
    void InitializeLinearSystem(ImplicitStrategyDataContainer<TLinearAlgebra>& rImplicitStrategyDataContainer)
    {
        KRATOS_TRY

        // Setting up the DOFs list
        BuiltinTimer setup_dofs_time;
        auto p_dof_set = rImplicitStrategyDataContainer.pDofSet;
        auto p_eff_dof_set = rImplicitStrategyDataContainer.pEffectiveDofSet;
        auto [eq_system_size, eff_eq_system_size] = this->SetUpDofArrays(p_dof_set, p_eff_dof_set);
        KRATOS_INFO_IF("ImplicitScheme", this->GetEchoLevel() > 0) << "Setup DOFs Time: " << setup_dofs_time << std::endl;

        // Set up the equation ids
        BuiltinTimer setup_system_ids_time;
        this->SetUpSystemIds(p_dof_set, p_eff_dof_set);
        KRATOS_INFO_IF("ImplicitScheme", this->GetEchoLevel() > 0) << "Set up system time: " << setup_system_ids_time << std::endl;
        KRATOS_INFO_IF("ImplicitScheme", this->GetEchoLevel() > 0) << "Equation system size: " << eq_system_size << std::endl;
        KRATOS_INFO_IF("ImplicitScheme", this->GetEchoLevel() > 0) << "Effective equation system size: " << eff_eq_system_size << std::endl;

        // Allocating the system constraints arrays
        BuiltinTimer constraints_allocation_time;
        (this->GetBuilder()).AllocateLinearSystemConstraints(rImplicitStrategyDataContainer);
        KRATOS_INFO_IF("ImplicitScheme", this->GetEchoLevel() > 0) << "Linear system constraints allocation time: " << constraints_allocation_time << std::endl;

        // Call the builder to allocate and initialize the system vectors
        BuiltinTimer linear_system_allocation_time;
        (this->GetBuilder()).AllocateLinearSystem(rImplicitStrategyDataContainer);
        KRATOS_INFO_IF("ImplicitScheme", this->GetEchoLevel() > 0) << "Linear system allocation time: " << linear_system_allocation_time << std::endl;

        KRATOS_CATCH("")
    }

    template<class TEntityType>
    bool CalculateLocalSystemContribution(
        TEntityType& rEntity,
        TLSType& rTLS,
        const ProcessInfo &rProcessInfo)
    {
        if (rEntity.IsActive()) {
            // Calculate local RHS contribution
            rEntity.CalculateLocalSystem(rTLS.LocalMatrix, rTLS.LocalVector, rProcessInfo);

            // Get the positions in the global system
            rEntity.EquationIdVector(rTLS.LocalEqIds, rProcessInfo);

            // The element is active and is to be assembled
            return true;
        } else {
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
        if (rEntity.IsActive()) {
            // Calculate local RHS contribution
            rEntity.CalculateRightHandSide(rTLS.LocalVector, rProcessInfo);

            // Get the positions in the global system
            rEntity.EquationIdVector(rTLS.LocalEqIds, rProcessInfo);

            // The element is active and is to be assembled
            return true;
        } else {
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
        if (rEntity.IsActive()) {
            // Calculate local RHS contribution
            rEntity.CalculateLeftHandSide(rTLS.LocalMatrix, rProcessInfo);

            // Get the positions in the global system
            rEntity.EquationIdVector(rTLS.LocalEqIds, rProcessInfo);

            // The element is active and is to be assembled
            return true;
        } else {
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
        if (rEntity.IsActive()) {
            // Calculate local RHS contribution
            rEntity.CalculateMassMatrix(rTLS.LocalMatrix, rProcessInfo);

            // Get the positions in the global system
            rEntity.EquationIdVector(rTLS.LocalEqIds, rProcessInfo);

            // The element is active and is to be assembled
            return true;
        } else {
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
        if (rEntity.IsActive()) {
            // Calculate local RHS contribution
            rEntity.CalculateDampingMatrix(rTLS.LocalMatrix, rProcessInfo);

            // Get the positions in the global system
            rEntity.EquationIdVector(rTLS.LocalEqIds, rProcessInfo);

            // The element is active and is to be assembled
            return true;
        } else {
            // The element is inactive and is not to be assembled
            return false;
        }
    }

    bool CalculateConstraintContribution(
        ModelPart::MasterSlaveConstraintType& rConstraint,
        TLSType& rTLS,
        const ProcessInfo &rProcessInfo)
    {
        if (rConstraint.IsActive()) {
            // Calculate local RHS contribution
            rConstraint.CalculateLocalSystem(rTLS.LocalMatrix, rTLS.LocalVector, rProcessInfo);

            // The constraint is active and is to be assembled
            return true;
        } else {
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
        mReformDofsAtEachStep = ThisParameters["reform_dofs_at_each_step"].GetBool();
    }

    ///@}
    ///@name Protected  Access
    ///@{

    /**
     * @brief Get the Builder object
     * @return BuilderType& Reference to the builder class
     */
    BuilderType& GetBuilder()
    {
        return *mpBuilder;
    }

    /**
     * @brief Get the Builder object
     * @return const BuilderType& Reference to the builder class
     */
    const BuilderType& GetBuilder() const
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

    int mEchoLevel = 0; /// The level of verbosity

    bool mMoveMesh = false; /// Flag to activate the mesh motion from the DISPLACEMENT variable

    bool mReformDofsAtEachStep = false; /// Flag to indicate if the DOF sets are required to be computed at each time step

    ModelPart* mpModelPart = nullptr; /// Pointer to the ModelPart the scheme refers to

    typename BuilderType::UniquePointer mpBuilder = nullptr; /// Pointer to the corresponding builder

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

} // namespace Kratos::Future

