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
//                   Eduard Gomez
//
//

#pragma once

// System includes
#include <numeric>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "factories/factory.h"
#include "solving_strategies/strategies/explicit_solving_strategy.h"
#include "butcher_tableau.h"

namespace Kratos {

///@name Kratos Classes
///@{

/** @brief Family of explicit Runge-Kutta schemes
 *
 * @details:
 * Formulation:
 *
 * - The differential equation is u_t = f(t, u)
 *
 * The Runge-Kutta method is, for i = 0...N substeps:
 * - k^(i) = f(t + c_i*dt, u^(i-1))              -> Where u^(0) := u^n
 * - u^(i) = u^n + \sum_{j=1}^{i-1} A_{ij} k^(j) -> Intermediate steps. u^(N) is not needed, therefore neither is A[N, :]
 * - u^{n+1} = u^n + \sum_{i} b_{i} k(u^(i))     -> Solution
 *
 * @tparam TSparseSpace
 * @tparam TDenseSpace
 * @tparam TButcherTableau specifies
 *  - The sets of coefficients A, b and c
 *  - The number of Runge-Kutta substeps
 *
 * @see: ButcherTableau
 */
template <class TSparseSpace, class TDenseSpace, class TButcherTableau>
class ExplicitSolvingStrategyRungeKutta : public ExplicitSolvingStrategy<TSparseSpace, TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{

    /// The base solving strategy class definition
    using SolvingStrategyType = SolvingStrategy<TSparseSpace, TDenseSpace>;

    /// The base class definition
    using BaseType = ExplicitSolvingStrategy<TSparseSpace, TDenseSpace>;

    /// The definition of the current class
    using ClassType = ExplicitSolvingStrategyRungeKutta<TSparseSpace, TDenseSpace, TButcherTableau>;

    /// The explicit builder and solver definition
    using ExplicitBuilderType = typename BaseType::ExplicitBuilderType;

    /// The DOF type
    using DofType = typename BaseType::DofType;

    /// The local vector definition
    using LocalSystemVectorType = typename TDenseSpace::VectorType;
    using LocalSystemMatrixType = typename TDenseSpace::MatrixType;

    /** Counted pointer of ClassName */
    KRATOS_CLASS_POINTER_DEFINITION(ExplicitSolvingStrategyRungeKutta);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor. (empty)
     */
    explicit ExplicitSolvingStrategyRungeKutta()
        : BaseType()
    {
    }

    /**
     * @brief Default constructor. (with parameters)
     * @param rModelPart The model part of the problem
     * @param ThisParameters The configuration parameters
     */
    explicit ExplicitSolvingStrategyRungeKutta(
        ModelPart &rModelPart,
        Parameters ThisParameters)
        : BaseType(rModelPart)
    {
        // Validate and assign defaults
        ThisParameters = this->ValidateAndAssignParameters(ThisParameters, this->GetDefaultParameters());
        this->AssignSettings(ThisParameters);
    }

    /**
     * @brief Default constructor.
     * @param rModelPart The model part to be computed
     * @param pExplicitBuilder The pointer to the explicit builder and solver
     * @param MoveMeshFlag The flag to set if the mesh is moved or not
     */
    explicit ExplicitSolvingStrategyRungeKutta(
        ModelPart &rModelPart,
        typename ExplicitBuilderType::Pointer pExplicitBuilder,
        bool MoveMeshFlag = false,
        int RebuildLevel = 0)
        : BaseType(rModelPart, pExplicitBuilder, MoveMeshFlag, RebuildLevel)
    {
    }

    /**
     * @brief Default constructor.
     * @param rModelPart The model part to be computed
     * @param MoveMeshFlag The flag to set if the mesh is moved or not
     */
    explicit ExplicitSolvingStrategyRungeKutta(
        ModelPart &rModelPart,
        bool MoveMeshFlag = false,
        int RebuildLevel = 0)
        : BaseType(rModelPart, MoveMeshFlag, RebuildLevel)
    {
    }

    /**
     * @brief Create method
     * @param rModelPart The model part to be computed
     * @param ThisParameters The configuration parameters
     */
    typename SolvingStrategyType::Pointer Create(
        ModelPart& rModelPart,
        Parameters ThisParameters
        ) const override
    {
        return Kratos::make_shared<ClassType>(rModelPart, ThisParameters);
    }

    /** Copy constructor.
     */
    ExplicitSolvingStrategyRungeKutta(const ExplicitSolvingStrategyRungeKutta &Other) = delete;

    /** Destructor.
     */
    ~ExplicitSolvingStrategyRungeKutta() override = default;

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     * @return The default parameters
     */
    Parameters GetDefaultParameters() const override
    {
        Parameters default_parameters = Parameters(R"(
        {
            "name" : "explicit_solving_strategy_runge_kutta"
        })");

        // Getting base class default parameters
        const Parameters base_default_parameters = BaseType::GetDefaultParameters();
        default_parameters.RecursivelyAddMissingParameters(base_default_parameters);
        return default_parameters;
    }

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    void Initialize() override
    {
        BaseType::Initialize();
        BaseType::GetModelPart().GetProcessInfo().SetValue(TIME_INTEGRATION_THETA, 0.0);
    }

    ///@}
    ///@name Input and output
    ///@{

    /**
     * @brief Returns the name of the class as used in the settings (snake_case format)
     * @return The name of the class
     */
    static std::string Name()
    {
        std::stringstream s;
        s << "explicit_solving_strategy_runge_kutta_" << TButcherTableau::Name();
        return s.str();
    }

    /// Return information as a string.
    std::string Info() const override
    {
        std::stringstream ss;
        ss << "ExplicitSolvingStrategyRungeKutta<" << mButcherTableau.Info() << ">";
        return ss.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream &rOStream) const override
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream &rOStream) const override
    {
        rOStream << Info();
    }

    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{
    const TButcherTableau mButcherTableau;


    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{

    void SolveWithLumpedMassMatrix() override
    {
        KRATOS_TRY

        // Get the required data from the explicit builder and solver
        auto& r_explicit_bs = BaseType::GetExplicitBuilder();
        auto& r_dof_set = r_explicit_bs.GetDofSet();
        const unsigned int dof_size = r_explicit_bs.GetEquationSystemSize();
        const auto& r_lumped_mass_vector = r_explicit_bs.GetLumpedMassMatrixVector();

        // Set the auxiliary RK vectors
        LocalSystemVectorType u_n(dof_size); // TODO: THIS IS INEFICCIENT. CREATE A UNORDERED_SET WITH THE IDOF AND VALUE AS ENTRIES. THIS HAS TO BE OPTIONAL
        LocalSystemMatrixType rk_K(dof_size, TButcherTableau::SubstepCount());

        // Perform the RK update
        const double dt = BaseType::GetDeltaTime();
        KRATOS_ERROR_IF(dt < 1.0e-12) << "ProcessInfo DELTA_TIME is close to zero." << std::endl;

        // Set the previous step solution in the current buffer position
        // Note that we set the 0 position of the buffer to be equal to the values in step n (not n+1)
        // Additionally, we save in an auxiliary vector the value of the fixed DOFs, which is also taken from the previous time step
        IndexPartition<int>(r_dof_set.size()).for_each(
            [&](int i_dof){
                auto it_dof = r_dof_set.begin() + i_dof;
                double& r_u_0 = it_dof->GetSolutionStepValue(0);
                const double& r_u_1 = it_dof->GetSolutionStepValue(1);
                if (it_dof->IsFixed()) {
                    u_n(i_dof) = r_u_0;
                }
                r_u_0 = r_u_1;
            }
        );

        // Calculate the RK intermediate sub steps
        for(unsigned int i=1; i<TButcherTableau::SubstepCount(); ++i) {
            PerformRungeKuttaIntermediateSubStep(i, u_n, rk_K);
        }
        PerformRungeKuttaLastSubStep(rk_K);

        // Do the final solution update
        const auto& weights = mButcherTableau.GetWeights();
        IndexPartition<int>(r_dof_set.size()).for_each(
            [&](int i_dof){
                auto it_dof = r_dof_set.begin() + i_dof;
                // Do the DOF update
                double& r_u = it_dof->GetSolutionStepValue(0);
                const double& r_u_old = it_dof->GetSolutionStepValue(1);
                if (!it_dof->IsFixed()) {
                    const double mass = r_lumped_mass_vector(i_dof);
                    const MatrixRow<LocalSystemMatrixType> substeps_k = row(rk_K, i_dof);
                    r_u = r_u_old + (dt / mass) * std::inner_product(weights.begin(), weights.end(), substeps_k.begin(), 0.0);
                } else {
                    r_u = u_n(i_dof);
                }
            }
        );

        KRATOS_CATCH("");
    }

    /**
     * @brief Initialize the Runge-Kutta intermediate substep
     * This method is intended to implement all the operations required before each Runge-Kutta intermediate substep
     */
    virtual void InitializeRungeKuttaIntermediateSubStep() {};

    /**
     * @brief Finalize the Runge-Kutta intermediate substep
     * This method is intended to implement all the operations required after each Runge-Kutta intermediate substep
     */
    virtual void FinalizeRungeKuttaIntermediateSubStep() {};

    /**
     * @brief Initialize the Runge-Kutta last substep
     * This method is intended to implement all the operations required before each Runge-Kutta last substep
     */
    virtual void InitializeRungeKuttaLastSubStep() {};

    /**
     * @brief Finalize the Runge-Kutta last substep
     * This method is intended to implement all the operations required after each Runge-Kutta last substep
     */
    virtual void FinalizeRungeKuttaLastSubStep() {};

    /**
     * @brief Performs an intermediate RK4 step
     * This functions performs all the operations required in an intermediate RK4 sub step
     * @param SubStepIndex The sub step index
     * @param SubStepCoefficients The sub step coefficients (these are saved as member variables)
     * @param rFixedDofsValues The vector containing the step n+1 values of the fixed DOFs
     * @param rIntermediateStepResidualVector The vector to store the intermediate sub step residual
     */
    virtual void PerformRungeKuttaIntermediateSubStep(
        const IndexType SubStepIndex,
        const LocalSystemVectorType& rFixedDofsValues,
        LocalSystemMatrixType& rIntermediateStepResidualVectors)
    {
        KRATOS_TRY

        // Get the required data from the explicit builder and solver
        auto& r_explicit_bs = BaseType::GetExplicitBuilder();
        auto& r_dof_set = r_explicit_bs.GetDofSet();
        const auto& r_lumped_mass_vector = r_explicit_bs.GetLumpedMassMatrixVector();

        // Get model part and information
        const double dt = BaseType::GetDeltaTime();
        KRATOS_ERROR_IF(dt < 1.0e-12) << "ProcessInfo DELTA_TIME is close to zero." << std::endl;
        auto& r_model_part = BaseType::GetModelPart();
        auto& r_process_info = r_model_part.GetProcessInfo();

        // Fetch this substeps's values from tableau
        const double integration_theta = mButcherTableau.GetIntegrationTheta(SubStepIndex);
        const auto alphas_begin = mButcherTableau.GetMatrixRowBegin(SubStepIndex); // Runge kutta matrix row begin
        const auto alphas_end = mButcherTableau.GetMatrixRowEnd(SubStepIndex); // Runge kutta matrix row end

        // Set the RUNGE_KUTTA_STEP value. This has to be done prior to the InitializeRungeKuttaStep()
        r_process_info.GetValue(RUNGE_KUTTA_STEP) = SubStepIndex;
        r_process_info.GetValue(TIME_INTEGRATION_THETA) = integration_theta;

        // Perform the intermidate sub step update
        InitializeRungeKuttaIntermediateSubStep();
        r_explicit_bs.BuildRHS(r_model_part);

        IndexPartition<int>(r_dof_set.size()).for_each(
            [&](int i_dof){
                auto it_dof = r_dof_set.begin() + i_dof;
                // Save current value in the corresponding vector
                const double& r_res = it_dof->GetSolutionStepReactionValue();
                rIntermediateStepResidualVectors(i_dof, SubStepIndex-1) = r_res;
                // Do the DOF update
                double& r_u = it_dof->GetSolutionStepValue(0);
                const double& r_u_old = it_dof->GetSolutionStepValue(1);
                if (!it_dof->IsFixed()) {
                    const double mass = r_lumped_mass_vector(i_dof);
                    const auto k = row(rIntermediateStepResidualVectors, i_dof);
                    r_u = r_u_old + (dt / mass) * std::inner_product(alphas_begin, alphas_end, k.begin(), 0.0);
                    /*                            ^~~~~~~~~~~~~~~~~~
                     * Using std::inner_product instead of boost's inner_prod because it allows us to
                     * chose a begin and, more importantly, an end.
                     *
                     * This is useful because alpha_ij = 0 for j >= i, hence the tail end of this scalar product
                     * can be ignored by setting the past-the-end iterator at j=i
                     */
                } else {
                    const double delta_u = rFixedDofsValues(i_dof) - r_u_old;
                    r_u = r_u_old + integration_theta * delta_u;
                }
            }
        );

        FinalizeRungeKuttaIntermediateSubStep();

        KRATOS_CATCH("SubstepIndex = " + std::to_string(SubStepIndex));
    }

    /**
     * @brief Performs the last RK4 step
     * This functions performs all the operations required in the last RK4 sub step
     * @param rLastStepResidualVector The vector to store the last sub step residual
     */
    virtual void PerformRungeKuttaLastSubStep(LocalSystemMatrixType& rLastStepResidualVector)
    {
        KRATOS_TRY

        // Get the required data from the explicit builder and solver
        auto& r_explicit_bs = BaseType::GetExplicitBuilder();
        auto& r_dof_set = r_explicit_bs.GetDofSet();

        // Get model part
        auto& r_model_part = BaseType::GetModelPart();
        auto& r_process_info = r_model_part.GetProcessInfo();
        constexpr unsigned int substep_index = TButcherTableau::SubstepCount();

        // Set the RUNGE_KUTTA_STEP value. This has to be done prior to the InitializeRungeKuttaStep()
        r_process_info.GetValue(RUNGE_KUTTA_STEP) = TButcherTableau::SubstepCount();
        r_process_info.GetValue(TIME_INTEGRATION_THETA) = mButcherTableau.GetIntegrationTheta(substep_index);

        // Perform the last sub step residual calculation
        InitializeRungeKuttaLastSubStep();
        r_explicit_bs.BuildRHS(r_model_part);

        IndexPartition<int>(r_dof_set.size()).for_each(
            [&](int i_dof){
                const auto it_dof = r_dof_set.begin() + i_dof;
                // Save current value in the corresponding vector
                const double& r_res = it_dof->GetSolutionStepReactionValue();
                rLastStepResidualVector(i_dof, substep_index - 1) = r_res;
            }
        );

        FinalizeRungeKuttaLastSubStep();

        KRATOS_CATCH("");
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
};

///@}

///@name Type Definitions
///@{

template <class TSparseSpace, class TDenseSpace>
using ExplicitSolvingStrategyRungeKutta4 = ExplicitSolvingStrategyRungeKutta<TSparseSpace, TDenseSpace, ButcherTableauRK4>;

template <class TSparseSpace, class TDenseSpace>
using ExplicitSolvingStrategyRungeKutta3TVD = ExplicitSolvingStrategyRungeKutta<TSparseSpace, TDenseSpace, ButcherTableauRK3TVD>;

template <class TSparseSpace, class TDenseSpace>
using ExplicitSolvingStrategyRungeKutta2 = ExplicitSolvingStrategyRungeKutta<TSparseSpace, TDenseSpace, ButcherTableauMidPointMethod>;

template <class TSparseSpace, class TDenseSpace>
using ExplicitSolvingStrategyRungeKutta1 = ExplicitSolvingStrategyRungeKutta<TSparseSpace, TDenseSpace, ButcherTableauForwardEuler>;

///@}

} /* namespace Kratos.*/
