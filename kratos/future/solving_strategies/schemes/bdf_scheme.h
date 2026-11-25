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
//

#pragma once

// System includes
#include <array>
#include <vector>

// External includes

// Project includes
#include "includes/checks.h"
#include "includes/process_info.h"
#include "includes/variables.h"
#include "utilities/time_discretization.h"

#ifdef KRATOS_USE_FUTURE
#include "future/solving_strategies/schemes/implicit_scheme.h"
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

/**
 * @brief BDF scheme TLS type definition
 * Thread Local Storage container used in the parallel assembly of the dynamic (BDF) implicit problems
 * @tparam DataType data type of the problem to be solved
 */
template<class TDataType = double >
struct BdfThreadLocalStorage
{
    // Local LHS contribution
    DenseMatrix<TDataType> LocalMatrix;

    // Local RHS contribution
    DenseVector<TDataType> LocalVector;

    // Scratch for the local mass matrix
    DenseMatrix<TDataType> MassMatrix;

    // Scratch for the local damping matrix
    DenseMatrix<TDataType> DampingMatrix;

    // Scratch for the first derivatives (velocity) at the reference step
    DenseVector<TDataType> FirstDerivatives;

    // Scratch for the second derivatives (acceleration) at the reference step
    DenseVector<TDataType> SecondDerivatives;

    // Vector containing the localization in the system of the different terms
    Element::EquationIdVectorType LocalEqIds;

    // Vector containing the slave equation ids
    MasterSlaveConstraint::EquationIdVectorType SlaveEqIds;

    // Vector containing the master equation ids
    MasterSlaveConstraint::EquationIdVectorType MasterEqIds;
};

/**
 * @class BDFScheme
 * @ingroup KratosCore
 * @brief Backward Differentiation Formula (BDF) time integration scheme for dynamic implicit problems
 * @tparam TLinearAlgebra The linear algebra type
 * @author Ruben Zorrilla
 */
template<class TLinearAlgebra>
class BDFScheme : public ImplicitScheme<TLinearAlgebra>
{
public:

    // FIXME: Does not work... ask @Charlie
    // /// Add scheme to Kratos registry
    // KRATOS_REGISTRY_ADD_TEMPLATE_PROTOTYPE("Schemes.KratosMultiphysics", BDFScheme, BDFScheme, TLinearAlgebra)
    // KRATOS_REGISTRY_ADD_TEMPLATE_PROTOTYPE("Schemes.All", BDFScheme, BDFScheme, TLinearAlgebra)

    ///@name Type Definitions
    ///@{

    /// Pointer definition of BDFScheme
    KRATOS_CLASS_POINTER_DEFINITION(BDFScheme);

    /// The definition of the current class
    using BaseType = ImplicitScheme<TLinearAlgebra>;

    /// Index type definition
    using IndexType = typename TLinearAlgebra::IndexType;

    /// Data type definition
    using DataType = typename TLinearAlgebra::DataType;

    /// Matrix type definition
    using MatrixType = typename TLinearAlgebra::MatrixType;

    /// Vector type definition
    using VectorType = typename TLinearAlgebra::VectorType;

    /// TLS type
    using TLSType = BdfThreadLocalStorage<DataType>;

    /// DoF type definition
    using DofType = Dof<DataType>;

    /// DoF array type definition
    using DofsArrayType = ModelPart::DofsArrayType;

    /// BDF coefficients type
    using BdfCoefficientsType = Vector;

    ///@}
    ///@name Life Cycle
    ///@{

    /// @brief Default constructor
    explicit BDFScheme() = default;

    /// @brief Constructor with parameters
    /// @param rModelPart Reference to the model part
    /// @param ThisParameters Parameters object encapsulating the settings
    explicit BDFScheme(
        ModelPart& rModelPart,
        Parameters ThisParameters)
        : BaseType(rModelPart)
    {
        // Validate the parameters against the default ones and assign the settings
        this->ValidateAndAssignSettings(ThisParameters);
    }

    /// @brief Copy constructor
    /// @param rOther Other BDFScheme
    explicit BDFScheme(BDFScheme& rOther)
        : BaseType(rOther)
        , mOrder(rOther.mOrder)
        , mBDFCoefficients(rOther.mBDFCoefficients)
        , mpBDF(Kratos::make_unique<TimeDiscretization::BDF>(rOther.mOrder))
        , mSolutionVariables(rOther.mSolutionVariables)
        , mFirstDerivativeVariables(rOther.mFirstDerivativeVariables)
        , mSecondDerivativeVariables(rOther.mSecondDerivativeVariables)
    {
    }

    /// @brief Destructor
    virtual ~BDFScheme() = default;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    typename BaseType::Pointer Create(
        ModelPart& rModelPart,
        Parameters ThisParameters) const override
    {
        return Kratos::make_shared<BDFScheme<TLinearAlgebra>>(rModelPart, ThisParameters);
    }

    typename BaseType::Pointer Clone() override
    {
        return Kratos::make_shared<BDFScheme<TLinearAlgebra>>(*this);
    }

    /**
     * @brief Function called once at the beginning of each solution step
     * @details In addition to the base class operations (initializing the entities), it computes and stores the
     * BDF coefficients for the current time step. It also fixes the solution DOFs whose first or second time
     * derivative DOF is fixed.
     * @param rImplicitStrategyData Auxiliary container with the linear system arrays
     */
    void InitializeSolutionStep(ImplicitStrategyData<TLinearAlgebra>& rImplicitStrategyData) override
    {
        KRATOS_TRY

        // Initialize the entities (base class behavior)
        BaseType::InitializeSolutionStep(rImplicitStrategyData);

        // Get the model part process info
        auto& r_model_part = this->GetModelPart();
        auto& r_process_info = r_model_part.GetProcessInfo();

        // Check that the time step is available
        KRATOS_ERROR_IF_NOT(r_process_info.Has(DELTA_TIME))
            << "DELTA_TIME is not set in the ProcessInfo. The BDF coefficients cannot be computed." << std::endl;

        const double dt_0 = r_process_info[DELTA_TIME];

        // For orders higher than 2 a constant time step is required
        if (mOrder > 2) {
            const double dt_1 = r_process_info.GetPreviousTimeStepInfo(1)[DELTA_TIME];
            KRATOS_ERROR_IF(std::abs(dt_0 - dt_1) > 1.0e-10 * (dt_0 + dt_1))
                << "For BDF orders higher than 2 a constant time step is required. "
                << "Previous time step: " << dt_1 << ", current time step: " << dt_0 << std::endl;
        }

        // Compute and save the BDF coefficients in the model part ProcessInfo
        mpBDF->ComputeAndSaveBDFCoefficients(r_process_info);

        // Update the internal copy of the BDF coefficients
        mBDFCoefficients = r_process_info[BDF_COEFFICIENTS];

        // If a first or second time derivative DOF is fixed, the corresponding solution DOF is fixed too
        block_for_each(r_model_part.Nodes(), [this](Node& rNode){
            for (std::size_t i_var = 0; i_var < mSolutionVariables.size(); ++i_var) {
                const auto& r_solution_variable = *mSolutionVariables[i_var];
                const auto& r_first_derivative_variable = *mFirstDerivativeVariables[i_var];
                const auto& r_second_derivative_variable = *mSecondDerivativeVariables[i_var];

                if (rNode.HasDofFor(r_second_derivative_variable) && rNode.IsFixed(r_second_derivative_variable)) {
                    rNode.Fix(r_solution_variable);
                } else if (rNode.HasDofFor(r_first_derivative_variable) && rNode.IsFixed(r_first_derivative_variable)) {
                    rNode.Fix(r_solution_variable);
                }
            }
        });

        KRATOS_CATCH("")
    }

    /**
     * @brief Builds the full linear system (LHS and RHS) including the dynamic (BDF) contributions
     * @details The LHS reads K + c0^2*M + c0*D and the RHS reads (f - K*u) - M*a - D*v, where all the quantities
     * are evaluated at the predicted state of the current step.
     * @param rLHS Left-hand side matrix
     * @param rRHS Right-hand side vector
     */
    void Build(
        MatrixType& rLHS,
        VectorType& rRHS) override
    {
        Timer::Start("Build");

        const auto timer = BuiltinTimer();

        // Getting conditions and elements to be assembled
        const auto& r_model_part = this->GetModelPart();
        const auto& r_elems = r_model_part.Elements();
        const auto& r_conds = r_model_part.Conditions();
        const auto& r_process_info = r_model_part.GetProcessInfo();

        // Getting entities container data
        auto elems_begin = r_elems.begin();
        auto conds_begin = r_conds.begin();
        const std::size_t n_elems = r_elems.size();
        const std::size_t n_conds = r_conds.size();

        // Initialize RHS and LHS assembly
        rRHS.BeginAssemble();
        rLHS.BeginAssemble();

        // Assemble entities
        #pragma omp parallel
        {
            TLSType aux_tls;

            // Assemble elements
            # pragma omp for schedule(guided, 512) nowait
            for (int k = 0; k < n_elems; ++k) {
                // Calculate the local LHS and RHS contributions including the dynamics
                auto it_elem = elems_begin + k;
                const bool assemble = this->CalculateLocalSystemContributionWithDynamics(*it_elem, aux_tls, r_process_info);

                // Assemble the local contributions to the global system
                if (assemble) {
                    rRHS.Assemble(aux_tls.LocalVector, aux_tls.LocalEqIds);
                    rLHS.Assemble(aux_tls.LocalMatrix, aux_tls.LocalEqIds);
                }
            }

            // Assemble conditions
            # pragma omp for schedule(guided, 512)
            for (int k = 0; k < n_conds; ++k) {
                // Calculate the local LHS and RHS contributions including the dynamics
                auto it_cond = conds_begin + k;
                const bool assemble = this->CalculateLocalSystemContributionWithDynamics(*it_cond, aux_tls, r_process_info);

                // Assemble the local contributions to the global system
                if (assemble) {
                    rRHS.Assemble(aux_tls.LocalVector, aux_tls.LocalEqIds);
                    rLHS.Assemble(aux_tls.LocalMatrix, aux_tls.LocalEqIds);
                }
            }
        }

        // Finalize RHS and LHS assembly
        rRHS.FinalizeAssemble();
        rLHS.FinalizeAssemble();

        KRATOS_INFO_IF("BDFScheme", this->GetEchoLevel() >= 1) << "Build time: " << timer << std::endl;
        KRATOS_INFO_IF("BDFScheme", this->GetEchoLevel() >= 2) << "Finished parallel building" << std::endl;

        Timer::Stop("Build");
    }

    /**
     * @brief Builds only the right-hand side (RHS) including the dynamic (BDF) contributions
     * @details The RHS reads (f - K*u) - M*a - D*v, where all the quantities are evaluated at the predicted state
     * of the current step.
     * @param rRHS Right-hand side vector
     */
    void Build(VectorType& rRHS) override
    {
        Timer::Start("BuildRightHandSide");

        const auto timer = BuiltinTimer();

        // Getting conditions and elements to be assembled
        const auto& r_model_part = this->GetModelPart();
        const auto& r_elems = r_model_part.Elements();
        const auto& r_conds = r_model_part.Conditions();
        const auto& r_process_info = r_model_part.GetProcessInfo();

        // Getting entities container data
        auto elems_begin = r_elems.begin();
        auto conds_begin = r_conds.begin();
        const std::size_t n_elems = r_elems.size();
        const std::size_t n_conds = r_conds.size();

        // Initialize RHS assembly
        rRHS.BeginAssemble();

        // Assemble entities
        #pragma omp parallel
        {
            TLSType aux_tls;

            // Assemble elements
            # pragma omp for schedule(guided, 512) nowait
            for (int k = 0; k < n_elems; ++k) {
                // Calculate the local RHS contribution including the dynamics
                auto it_elem = elems_begin + k;
                const bool assemble = this->CalculateRightHandSideContributionWithDynamics(*it_elem, aux_tls, r_process_info);

                // Assemble the local contributions to the global system
                if (assemble) {
                    rRHS.Assemble(aux_tls.LocalVector, aux_tls.LocalEqIds);
                }
            }

            // Assemble conditions
            # pragma omp for schedule(guided, 512)
            for (int k = 0; k < n_conds; ++k) {
                // Calculate the local RHS contribution including the dynamics
                auto it_cond = conds_begin + k;
                const bool assemble = this->CalculateRightHandSideContributionWithDynamics(*it_cond, aux_tls, r_process_info);

                // Assemble the local contributions to the global system
                if (assemble) {
                    rRHS.Assemble(aux_tls.LocalVector, aux_tls.LocalEqIds);
                }
            }
        }

        // Finalize RHS assembly
        rRHS.FinalizeAssemble();

        KRATOS_INFO_IF("BDFScheme", this->GetEchoLevel() >= 1) << "Build time (RightHandSide only): " << timer << std::endl;
        KRATOS_INFO_IF("BDFScheme", this->GetEchoLevel() >= 2) << "Finished parallel building (RightHandSide only)" << std::endl;

        Timer::Stop("BuildRightHandSide");
    }

    /**
     * @brief Performs the update of the solution
     * @details Updates the displacement (primary variable) DOFs with the solution increment and then recomputes the
     * time derivatives (velocity and acceleration) from the updated displacement history with the BDF formula.
     * @param rImplicitStrategyData Auxiliary container with the linear system arrays
     */
    void Update(ImplicitStrategyData<TLinearAlgebra>& rImplicitStrategyData) override
    {
        KRATOS_TRY

        // Get linear system arrays
        auto& r_dx = *(rImplicitStrategyData.pGetLinearSystem()->pGetVector(LinearSystemTags::DenseVectorTag::Dx));
        auto& r_eff_dx = *(rImplicitStrategyData.pGetEffectiveLinearSystem()->pGetVector(LinearSystemTags::DenseVectorTag::Dx));
        auto& r_dof_set = *(rImplicitStrategyData.pGetDofSet());
        auto& r_eff_dof_set = *(rImplicitStrategyData.pGetEffectiveDofSet());

        // First update the constraints only DOFs with the effective solution vector
        this->UpdateConstraintsOnlyDofs(r_eff_dx, r_dof_set, r_eff_dof_set);

        // Get the solution update vector from the effective one
        this->CalculateUpdateVector(rImplicitStrategyData);

        // Update DOFs with solution values (note that we solve for the increments)
        block_for_each(r_dof_set, [&r_dx](DofType& rDof){
            if (rDof.IsFree()) {
                rDof.GetSolutionStepValue() += r_dx[rDof.EquationId()];
            }
        });

        // Update the time derivatives (velocity and acceleration) from the updated displacement history
        this->UpdateDerivatives();

        // If the mesh is to be updated, call the MoveMesh() method
        if (this->GetMoveMesh()) {
            this->MoveMesh();
        }

        KRATOS_CATCH("")
    }

    /**
     * @brief This function is designed to be called once to perform all the checks needed on the input provided
     * @details In addition to the base class checks, it verifies that the model part buffer size is large enough
     * to hold the history required by the selected BDF order.
     * @return 0 all OK, 1 otherwise
     */
    int Check() const override
    {
        KRATOS_TRY

        int check = BaseType::Check();

        // Check that the buffer size is large enough for the selected BDF order
        KRATOS_ERROR_IF(this->GetModelPart().GetBufferSize() < mOrder + 1)
            << "Insufficient buffer size for a BDF scheme of order " << mOrder << ". "
            << "The buffer size should be greater than " << (mOrder + 1) << ". Current size is "
            << this->GetModelPart().GetBufferSize() << std::endl;

        // Check that the configured variable lists are consistent
        KRATOS_ERROR_IF(mSolutionVariables.empty())
            << "The BDF scheme requires at least one solution variable." << std::endl;
        KRATOS_ERROR_IF(mFirstDerivativeVariables.empty())
            << "The BDF scheme requires at least one first derivative variable." << std::endl;
        KRATOS_ERROR_IF(mSecondDerivativeVariables.empty())
            << "The BDF scheme requires at least one second derivative variable." << std::endl;
        KRATOS_ERROR_IF(mSolutionVariables.size() != mFirstDerivativeVariables.size())
            << "The number of solution variables and first derivative variables must be equal." << std::endl;
        KRATOS_ERROR_IF(mSolutionVariables.size() != mSecondDerivativeVariables.size())
            << "The number of solution variables and second derivative variables must be equal." << std::endl;

        // Check that the configured variables are allocated in the nodal solution step data
        // and that the solution variables have a DOF in every node
        for (const auto& r_node : this->GetModelPart().Nodes()) {
            for (const auto p_variable : mSolutionVariables) {
                KRATOS_CHECK_VARIABLE_IN_NODAL_DATA((*p_variable), r_node);
                KRATOS_CHECK_DOF_IN_NODE((*p_variable), r_node);
            }
            for (const auto p_variable : mFirstDerivativeVariables) {
                KRATOS_CHECK_VARIABLE_IN_NODAL_DATA((*p_variable), r_node);
            }
            for (const auto p_variable : mSecondDerivativeVariables) {
                KRATOS_CHECK_VARIABLE_IN_NODAL_DATA((*p_variable), r_node);
            }
        }

        return check;

        KRATOS_CATCH("");
    }

    /**
     * @brief This method provides the default parameters to avoid conflicts between the different constructors
     * @return The default parameters
     */
    Parameters GetDefaultParameters() const override
    {
        // Current class default parameters
        Parameters default_parameters = Parameters(R"({
            "name" : "bdf_scheme",
            "integration_order" : 2,
            "solution_variables" : ["DISPLACEMENT"]
        })");

        // Add base class default parameters
        default_parameters.RecursivelyAddMissingParameters(BaseType::GetDefaultParameters());

        return default_parameters;
    }

    /**
     * @brief Returns the name of the class as used in the settings (snake_case format)
     * @return The name of the class
     */
    static std::string Name()
    {
        return "bdf_scheme";
    }

    ///@}
    ///@name Access
    ///@{

    /**
     * @brief Get the BDF integration order
     * @return std::size_t The BDF integration order
     */
    std::size_t GetOrder() const
    {
        return mOrder;
    }

    /**
     * @brief Get the BDF coefficients of the current time step
     * @return const BdfCoefficientsType& The BDF coefficients
     */
    const BdfCoefficientsType& GetBDFCoefficients() const
    {
        return mBDFCoefficients;
    }

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "BDFScheme (order " + std::to_string(mOrder) + ")";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
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
     * @brief This method assigns settings to member variables
     * @param ThisParameters Parameters that are assigned to the member variables
     */
    void AssignSettings(const Parameters ThisParameters) override
    {
        // Assign base scheme settings
        BaseType::AssignSettings(ThisParameters);

        // Assign the BDF settings
        mOrder = static_cast<std::size_t>(ThisParameters["integration_order"].GetInt());
        KRATOS_ERROR_IF(mOrder < 1) << "The BDF integration order must be at least 1. Current value: " << mOrder << std::endl;
        KRATOS_ERROR_IF(mOrder > 6) << "The BDF integration order must be at most 6. Current value: " << mOrder << std::endl;

        // Set up the BDF coefficients utility and reserve the coefficients
        mpBDF = Kratos::make_unique<TimeDiscretization::BDF>(mOrder);
        mBDFCoefficients = Vector(mOrder + 1, 0.0);

        // Create the solution and derivative variable lists from the settings
        this->CreateVariablesList(ThisParameters);
    }

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

    std::size_t mOrder = 2; /// The BDF integration order

    BdfCoefficientsType mBDFCoefficients; /// The BDF coefficients of the current time step

    TimeDiscretization::BDF::UniquePointer mpBDF = nullptr; /// Utility to compute the BDF coefficients

    std::vector<const Variable<double>*> mSolutionVariables; /// The solution variables

    std::vector<const Variable<double>*> mFirstDerivativeVariables; /// The first derivative variables

    std::vector<const Variable<double>*> mSecondDerivativeVariables; /// The second derivative variables

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    /**
     * @brief Creates the lists of solution and derivative variables from the scheme settings
     * @param ThisParameters Parameters containing the variable names
     */
    void CreateVariablesList(const Parameters ThisParameters)
    {
        KRATOS_TRY

        mSolutionVariables.clear();
        mFirstDerivativeVariables.clear();
        mSecondDerivativeVariables.clear();

        // Get the domain size from the model part process info
        const auto& r_process_info = this->GetModelPart().GetProcessInfo();
        KRATOS_ERROR_IF_NOT(r_process_info.Has(DOMAIN_SIZE))
            << "DOMAIN_SIZE is not set in the ProcessInfo. The BDF scheme cannot expand the vector variable components." << std::endl;
        const std::size_t domain_size = static_cast<std::size_t>(r_process_info.GetValue(DOMAIN_SIZE));
        KRATOS_ERROR_IF_NOT(domain_size == 2 || domain_size == 3)
            << "DOMAIN_SIZE must be 2 or 3. Current value: " << domain_size << std::endl;

        // Read the solution variable names from the settings
        const auto solution_variable_names = ThisParameters["solution_variables"].GetStringArray();

        // Expand each variable name into its components and add the corresponding time derivatives
        for (const auto& r_var_name : solution_variable_names) {
            if (!r_var_name.empty()) {
                AddVariableWithDerivatives(
                    r_var_name,
                    domain_size);
            }
        }

        // Validate the variable lists
        KRATOS_ERROR_IF(mSolutionVariables.empty())
            << "The BDF scheme requires at least one solution variable." << std::endl;
        KRATOS_ERROR_IF(mFirstDerivativeVariables.empty())
            << "The BDF scheme requires at least one first derivative variable." << std::endl;
        KRATOS_ERROR_IF(mSecondDerivativeVariables.empty())
            << "The BDF scheme requires at least one second derivative variable." << std::endl;
        KRATOS_ERROR_IF(mSolutionVariables.size() != mFirstDerivativeVariables.size())
            << "The number of solution variables and first derivative variables must be equal. "
            << "Solution variables: " << mSolutionVariables.size()
            << ", first derivative variables: " << mFirstDerivativeVariables.size() << std::endl;
        KRATOS_ERROR_IF(mSolutionVariables.size() != mSecondDerivativeVariables.size())
            << "The number of solution variables and second derivative variables must be equal. "
            << "Solution variables: " << mSolutionVariables.size()
            << ", second derivative variables: " << mSecondDerivativeVariables.size() << std::endl;

        KRATOS_CATCH("")
    }

    /**
     * @brief Adds the component variables of a scalar or vector variable to the lists
     * @details For each component the first and second time derivatives are also added to
     * the corresponding lists. Note that the solution variables are assumed to be second
     * order, that is, both time derivatives must be defined.
     * @param rName Name of the variable to expand
     * @param DomainSize Number of components to use
     */
    void AddVariableWithDerivatives(
        const std::string& rName,
        const std::size_t DomainSize)
    {
        // A scalar variable (including a component such as DISPLACEMENT_X)
        if (KratosComponents<Variable<double>>::Has(rName)) {
            AddScalarVariableWithDerivatives(KratosComponents<Variable<double>>::Get(rName));
            return;
        }

        // A vector variable whose components are expanded according to the domain size
        if (KratosComponents<Variable<array_1d<double, 3>>>::Has(rName)) {
            const std::array<std::string, 3> component_suffixes = {"_X", "_Y", "_Z"};
            for (std::size_t i_component = 0; i_component < DomainSize; ++i_component) {
                const std::string component_name = rName + component_suffixes[i_component];
                AddScalarVariableWithDerivatives(KratosComponents<Variable<double>>::Get(component_name));
            }
            return;
        }

        KRATOS_ERROR << "The variable " << rName << " is not registered." << std::endl;
    }

    /**
     * @brief Adds a variable and its first and second time derivatives to the lists
     * @param rVariable The solution variable
     */
    void AddScalarVariableWithDerivatives(const Variable<double>& rVariable)
    {
        const auto& r_first_derivative_var = rVariable.GetTimeDerivative();
        const auto& r_second_derivative_var = r_first_derivative_var.GetTimeDerivative();
        mSolutionVariables.push_back(&rVariable);
        mFirstDerivativeVariables.push_back(&r_first_derivative_var);
        mSecondDerivativeVariables.push_back(&r_second_derivative_var);
    }

    /**
     * @brief Calculates the local LHS and RHS contributions of an entity including the dynamic (BDF) terms
     * @details The LHS is K + c0^2*M + c0*D and the RHS is (f - K*u) - M*a - D*v, evaluated at the predicted state.
     * @tparam TEntityType The entity type (Element or Condition)
     * @param rEntity The entity to compute
     * @param rTLS The thread local storage holding the scratch arrays
     * @param rProcessInfo The current process info
     * @return bool True if the entity is active and its contributions are to be assembled
     */
    template<class TEntityType>
    bool CalculateLocalSystemContributionWithDynamics(
        TEntityType& rEntity,
        TLSType& rTLS,
        const ProcessInfo& rProcessInfo)
    {
        if (rEntity.IsActive()) {
            // Calculate the static (mechanical) local system: LHS = K, RHS = f - K*u
            rEntity.CalculateLocalSystem(rTLS.LocalMatrix, rTLS.LocalVector, rProcessInfo);

            // Get the positions in the global system
            rEntity.EquationIdVector(rTLS.LocalEqIds, rProcessInfo);

            // Calculate the local mass and damping matrices
            rEntity.CalculateMassMatrix(rTLS.MassMatrix, rProcessInfo);
            rEntity.CalculateDampingMatrix(rTLS.DampingMatrix, rProcessInfo);

            // Add the dynamic LHS contributions: LHS += c0^2 * M + c0 * D
            const double c0 = mBDFCoefficients[0];
            if (rTLS.MassMatrix.size1() != 0) {
                noalias(rTLS.LocalMatrix) += rTLS.MassMatrix * (c0 * c0);
            }
            if (rTLS.DampingMatrix.size1() != 0) {
                noalias(rTLS.LocalMatrix) += rTLS.DampingMatrix * c0;
            }

            // Add the dynamic RHS contributions: RHS -= M * a + D * v (evaluated at the predicted state)
            if (rTLS.MassMatrix.size1() != 0) {
                rEntity.GetSecondDerivativesVector(rTLS.SecondDerivatives, 0);
                noalias(rTLS.LocalVector) -= prod(rTLS.MassMatrix, rTLS.SecondDerivatives);
            }
            if (rTLS.DampingMatrix.size1() != 0) {
                rEntity.GetFirstDerivativesVector(rTLS.FirstDerivatives, 0);
                noalias(rTLS.LocalVector) -= prod(rTLS.DampingMatrix, rTLS.FirstDerivatives);
            }

            // The entity is active and is to be assembled
            return true;
        } else {
            // Clear current TLS values
            rTLS.LocalEqIds.clear();
            rTLS.LocalVector.clear();
            rTLS.LocalMatrix.clear();

            // The entity is inactive and is not to be assembled
            return false;
        }
    }

    /**
     * @brief Calculates the local RHS contribution of an entity including the dynamic (BDF) terms
     * @details The RHS is (f - K*u) - M*a - D*v, evaluated at the predicted state.
     * @tparam TEntityType The entity type (Element or Condition)
     * @param rEntity The entity to compute
     * @param rTLS The thread local storage holding the scratch arrays
     * @param rProcessInfo The current process info
     * @return bool True if the entity is active and its contributions are to be assembled
     */
    template<class TEntityType>
    bool CalculateRightHandSideContributionWithDynamics(
        TEntityType& rEntity,
        TLSType& rTLS,
        const ProcessInfo& rProcessInfo)
    {
        if (rEntity.IsActive()) {
            // Calculate the static (mechanical) RHS: RHS = f - K*u
            rEntity.CalculateRightHandSide(rTLS.LocalVector, rProcessInfo);

            // Get the positions in the global system
            rEntity.EquationIdVector(rTLS.LocalEqIds, rProcessInfo);

            // Calculate the local mass and damping matrices
            rEntity.CalculateMassMatrix(rTLS.MassMatrix, rProcessInfo);
            rEntity.CalculateDampingMatrix(rTLS.DampingMatrix, rProcessInfo);

            // Add the dynamic RHS contributions: RHS -= M * a + D * v (evaluated at the predicted state)
            if (rTLS.MassMatrix.size1() != 0) {
                rEntity.GetSecondDerivativesVector(rTLS.SecondDerivatives, 0);
                noalias(rTLS.LocalVector) -= prod(rTLS.MassMatrix, rTLS.SecondDerivatives);
            }
            if (rTLS.DampingMatrix.size1() != 0) {
                rEntity.GetFirstDerivativesVector(rTLS.FirstDerivatives, 0);
                noalias(rTLS.LocalVector) -= prod(rTLS.DampingMatrix, rTLS.FirstDerivatives);
            }

            // The entity is active and is to be assembled
            return true;
        } else {
            // Clear current TLS values
            rTLS.LocalEqIds.clear();
            rTLS.LocalVector.clear();

            // The entity is inactive and is not to be assembled
            return false;
        }
    }

    /**
     * @brief Performing the prediction of the solution
     * @details The generic implicit-scheme setup (initializing the linear systems and building the
     * master-slave constraints) is carried out by the base Predict through
     * PrepareLinearSystemAndConstraints. This method performs the scheme-specific BDF time
     * extrapolation:
     *  - The solution variables are extrapolated to the current step with a quadratic formula.
     *  - The master-slave constraints are then applied to the extrapolated solution.
     *  - Finally, the velocity and acceleration are recomputed from the constrained solution history.
     * @param rImplicitStrategyData Auxiliary container with the linear system arrays
     */
    void PredictImpl(ImplicitStrategyData<TLinearAlgebra>& rImplicitStrategyData) override
    {
        KRATOS_TRY

        // Get the model part and the current time step size
        const auto& r_model_part = this->GetModelPart();
        const auto& r_process_info = r_model_part.GetProcessInfo();
        const double delta_time = r_process_info[DELTA_TIME];

        // Extrapolate the solution variables to the current step (nodally for efficiency)
        // Note that the solution variables are assumed to be second order, so both the first and
        // second derivative variables are always available for each solution variable
        block_for_each(r_model_part.Nodes(), [this, delta_time](Node& rNode){
            for (std::size_t i_var = 0; i_var < mSolutionVariables.size(); ++i_var) {
                const auto& r_solution_variable = *mSolutionVariables[i_var];
                const auto& r_first_derivative_variable = *mFirstDerivativeVariables[i_var];
                const auto& r_second_derivative_variable = *mSecondDerivativeVariables[i_var];

                double& u0 = rNode.FastGetSolutionStepValue(r_solution_variable);
                double& v0 = rNode.FastGetSolutionStepValue(r_first_derivative_variable);

                if (rNode.HasDofFor(r_second_derivative_variable) &&
                    rNode.IsFixed(r_second_derivative_variable)) {
                    // The acceleration is fixed: compute the velocity and the solution from the
                    // fixed acceleration value using the BDF formula
                    const double a0 = rNode.FastGetSolutionStepValue(r_second_derivative_variable);
                    v0 = a0;
                    for (std::size_t i_order = 1; i_order < mOrder + 1; ++i_order) {
                        v0 -= mBDFCoefficients[i_order] *
                              rNode.FastGetSolutionStepValue(r_first_derivative_variable, i_order);
                    }
                    v0 /= mBDFCoefficients[0];

                    u0 = v0;
                    for (std::size_t i_order = 1; i_order < mOrder + 1; ++i_order) {
                        u0 -= mBDFCoefficients[i_order] *
                              rNode.FastGetSolutionStepValue(r_solution_variable, i_order);
                    }
                    u0 /= mBDFCoefficients[0];
                } else if (rNode.HasDofFor(r_first_derivative_variable) &&
                           rNode.IsFixed(r_first_derivative_variable)) {
                    // The velocity is fixed: compute the solution from the fixed velocity
                    // value using the BDF formula
                    u0 = v0;
                    for (std::size_t i_order = 1; i_order < mOrder + 1; ++i_order) {
                        u0 -= mBDFCoefficients[i_order] *
                              rNode.FastGetSolutionStepValue(r_solution_variable, i_order);
                    }
                    u0 /= mBDFCoefficients[0];
                } else if (!rNode.IsFixed(r_solution_variable)) {
                    // Free DOF: quadratic extrapolation from the previous step values
                    const double u1 = rNode.FastGetSolutionStepValue(r_solution_variable, 1);
                    const double v1 = rNode.FastGetSolutionStepValue(r_first_derivative_variable, 1);
                    const double a1 = rNode.FastGetSolutionStepValue(r_second_derivative_variable, 1);
                    u0 = u1 + delta_time * v1 + 0.5 * delta_time * delta_time * a1;
                }
            }
        });

        // Apply the master-slave constraints to the prediction
        const auto& r_comm = r_model_part.GetCommunicator().GetDataCommunicator();
        const unsigned int n_constraints_loc = r_model_part.MasterSlaveConstraints().size();
        const unsigned int n_constraints_glob = r_comm.SumAll(n_constraints_loc);

        auto& r_dof_set = *(rImplicitStrategyData.pGetDofSet());
        auto& r_eff_dof_set = *(rImplicitStrategyData.pGetEffectiveDofSet());

        if (n_constraints_glob != 0) {
            // Check that constraints have been built
            auto p_constraints_T = rImplicitStrategyData.pGetConstraintsT();
            auto p_constraints_Q = rImplicitStrategyData.pGetConstraintsQ();
            KRATOS_ERROR_IF(p_constraints_T == nullptr || p_constraints_Q == nullptr) << "The master-slave constraints have not been built yet." << std::endl;

            // Set the prediction effective vector from the effective DOFs values
            VectorType y_p(r_eff_dof_set.size());
            block_for_each(r_eff_dof_set, [&y_p](DofType& rEffDof){
                y_p[rEffDof.EffectiveEquationId()] = rEffDof.GetSolutionStepValue();
            });

            // Apply the prediction to the constraints only DOFs
            IndexPartition<IndexType>(this->mConstraintsOnlyDofs.size()).for_each([&](IndexType Index){
                auto p_constr_only_dof = this->mConstraintsOnlyDofs[Index];
                if (p_constr_only_dof->IsFree()) {
                    p_constr_only_dof->GetSolutionStepValue() = y_p[p_constr_only_dof->EffectiveEquationId()];
                }
            });

            // Apply the master-slave constraints to the prediction effective vector
            VectorType x_p(r_dof_set.size());
            (this->GetBuilder()).CalculateSolutionVector(y_p, *p_constraints_T, *p_constraints_Q, x_p); // x_p = T*y_p + q

            // Update the nodal DOFs with the prediction
            block_for_each(r_dof_set, [&x_p](DofType& rDof){
                if (rDof.IsFree()) {
                    rDof.GetSolutionStepValue() = x_p[rDof.EquationId()];
                }
            });
        }

        // Update the time derivatives (velocity and acceleration) from the predicted solution
        this->UpdateDerivatives();

        KRATOS_CATCH("")
    }

    /**
     * @brief Updates the time derivatives (velocity and acceleration) of all the nodes
     * @details It recomputes the velocity and the acceleration from the (updated) displacement history with the BDF formula.
     */
    void UpdateDerivatives()
    {
        const auto& r_model_part = this->GetModelPart();
        block_for_each(r_model_part.Nodes(), [this](Node& rNode){

            for (std::size_t i_var = 0; i_var < mSolutionVariables.size(); ++i_var) {
                const auto& r_solution_variable = *mSolutionVariables[i_var];
                const auto& r_first_derivative_variable = *mFirstDerivativeVariables[i_var];

                // Update the first derivative (velocity)
                double first_derivative = mBDFCoefficients[0] * rNode.FastGetSolutionStepValue(r_solution_variable);
                for (std::size_t i_order = 1; i_order < mOrder + 1; ++i_order) {
                    first_derivative += mBDFCoefficients[i_order] * rNode.FastGetSolutionStepValue(r_solution_variable, i_order);
                }
                rNode.FastGetSolutionStepValue(r_first_derivative_variable) = first_derivative;

                // Update the second derivative (acceleration)
                const auto& r_second_derivative_variable = *mSecondDerivativeVariables[i_var];

                double second_derivative = mBDFCoefficients[0] * rNode.FastGetSolutionStepValue(r_first_derivative_variable);
                for (std::size_t i_order = 1; i_order < mOrder + 1; ++i_order) {
                    second_derivative += mBDFCoefficients[i_order] * rNode.FastGetSolutionStepValue(r_first_derivative_variable, i_order);
                }
                rNode.FastGetSolutionStepValue(r_second_derivative_variable) = second_derivative;
            }
        });
    }

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
}; // Class BDFScheme

} // namespace Kratos::Future.
