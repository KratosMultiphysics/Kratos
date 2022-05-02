//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla, Edurd Gómez
//
//

#if !defined(KRATOS_EXPLICIT_SOLVING_STRATEGY_BFECC)
#define KRATOS_EXPLICIT_SOLVING_STRATEGY_BFECC

/* System includes */
#include <numeric>

/* External includes */

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "factories/factory.h"
#include "solving_strategies/strategies/explicit_solving_strategy.h"

namespace Kratos
{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/** @brief Family of explicit Runge-Kutta schemes
 *
 * @details:
 * Formulation:
 *
 * WIP
 *
 * @tparam TSparseSpace
 * @tparam TDenseSpace
 */
template <class TSparseSpace, class TDenseSpace>
class ExplicitSolvingStrategyBFECC : public ExplicitSolvingStrategy<TSparseSpace, TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{

    enum Substep { FORWARD, BACKWARD, FINAL };

    typedef ModelPart::SizeType SizeType;

    // The base solving strategy class definition
    typedef SolvingStrategy<TSparseSpace, TDenseSpace> SolvingStrategyType;

    // The base class definition
    typedef ExplicitSolvingStrategy<TSparseSpace, TDenseSpace> BaseType;

    /// The definition of the current class
    typedef ExplicitSolvingStrategyBFECC<TSparseSpace, TDenseSpace> ClassType;

    // The explicit builder and solver definition
    typedef typename BaseType::ExplicitBuilderType ExplicitBuilderType;

    /// The DOF type
    typedef typename BaseType::DofType DofType;

    /// The local vector definition
    typedef typename TDenseSpace::VectorType LocalSystemVectorType;
    typedef typename TDenseSpace::MatrixType LocalSystemMatrixType;

    /** Counted pointer of ClassName */
    KRATOS_CLASS_POINTER_DEFINITION(ExplicitSolvingStrategyBFECC);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor. (empty)
     */
    explicit ExplicitSolvingStrategyBFECC()
        : BaseType()
    {
    }

    /**
     * @brief Default constructor. (with parameters)
     * @param rModelPart The model part of the problem
     * @param ThisParameters The configuration parameters
     */
    explicit ExplicitSolvingStrategyBFECC(
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
    explicit ExplicitSolvingStrategyBFECC(
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
    explicit ExplicitSolvingStrategyBFECC(
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
    ExplicitSolvingStrategyBFECC(const ExplicitSolvingStrategyBFECC &Other) = delete;

    /** Destructor.
     */
    ~ExplicitSolvingStrategyBFECC() override = default;

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     * @return The default parameters
     */
    virtual Parameters GetDefaultParameters() const override
    {
        Parameters default_parameters = Parameters(R"(
        {
            "name" : "explicit_solving_strategy_bfecc"
        })");

        // Getting base class default parameters
        const Parameters base_default_parameters = BaseType::GetDefaultParameters();
        default_parameters.RecursivelyAddMissingParameters(base_default_parameters);
        return default_parameters;
    }

    /**
     * @brief This method assigns settings to member variables
     * @param ThisParameters Parameters that are assigned to the member variables
     */
    virtual void AssignSettings(const Parameters ThisParameters) override
    {
        BaseType::AssignSettings(ThisParameters);
    }

    /**
     * @brief Returns the name of the class as used in the settings (snake_case format)
     * @return The name of the class
     */
    static std::string Name()
    {
        std::stringstream s;
        s << "explicit_solving_strategy_bfecc";
        return s.str();
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

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream ss;
        ss << "ExplicitSolvingStrategyBFECC";
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

        // Validate time step size
        const double dt = BaseType::GetDeltaTime();
        KRATOS_ERROR_IF(dt < 1.0e-12) << "ProcessInfo DELTA_TIME is close to zero." << std::endl;

        // Stash the prev step solution
        const auto original_starting_values = ExtractSolutionStepData(1);

        // Take step forward
        {
            const auto u_fixed = CopySolutionStepData(1, 0);
            PerformSubstep(u_fixed, FORWARD);
        }

        // Overwrite previous step solution with backwards step
        {
            const auto u_fixed = CopySolutionStepData(0, 1);
            PerformSubstep(u_fixed, BACKWARD);
        }

        // Correct error with the previous step solution
        CorrectErrorAfterForwardsAndBackwards(original_starting_values);

        // Final update
        {
            const auto u_fixed = CopySolutionStepData(1, 0);
            PerformSubstep(u_fixed, FINAL);
        }

        KRATOS_CATCH("");
    }

    /**
     * Ovrwrite the destination buffer position with data from the source buffer position.
     * Additionally, we save in an auxiliary vector the value of the fixed DOFs in the destination buffer position.
     */
    LocalSystemVectorType CopySolutionStepData(const SizeType source, const SizeType destination)
    {
        // Get the required data from the explicit builder and solver
        auto& explicit_bs = BaseType::GetExplicitBuilder();
        auto& r_dof_set = explicit_bs.GetDofSet();
        const SizeType dof_size = explicit_bs.GetEquationSystemSize();

        LocalSystemVectorType u_fixed(dof_size);
        IndexPartition<SizeType>(r_dof_set.size()).for_each(
            [&](const SizeType i_dof){
                auto it_dof = r_dof_set.begin() + i_dof;
                double& r_u_0 = it_dof->GetSolutionStepValue(destination);
                if (it_dof->IsFixed()) {
                    u_fixed(i_dof) = r_u_0;
                }
                r_u_0 = it_dof->GetSolutionStepValue(source);
            }
        );

        return u_fixed;
    }

    LocalSystemVectorType ExtractSolutionStepData(const SizeType BufferPosition) const
    {
        // Get the required data from the explicit builder and solver
        auto& explicit_bs = BaseType::GetExplicitBuilder();
        auto& r_dof_set = explicit_bs.GetDofSet();
        const SizeType dof_size = explicit_bs.GetEquationSystemSize();

        LocalSystemVectorType u_copy(dof_size);
        IndexPartition<SizeType>(r_dof_set.size()).for_each(
            [&](const SizeType i_dof){
                const auto it_dof = r_dof_set.cbegin() + i_dof;
                u_copy[i_dof] = it_dof->GetSolutionStepValue(BufferPosition);
            }
        );

        return u_copy;
    }

    /**
     * @brief Initialize the BFECC initial forward substep
     * This method is intended to implement all the operations required before each BFECC initial forward substep
     */
    virtual void InitializeBFECCForwardSubstep() {};

    /**
     * @brief Finalize the BFECC initial forward substep
     * This method is intended to implement all the operations required after each BFECC initial forward substep
     */
    virtual void FinalizeBFECCForwardSubstep() {};

    /**
     * @brief Initialize the BFECC backward substep
     * This method is intended to implement all the operations required before each BFECC backward substep
     */
    virtual void InitializeBFECCBackwardSubstep() {};

    /**
     * @brief Finalize the BFECC backward substep
     * This method is intended to implement all the operations required after each BFECC backward substep
     */
    virtual void FinalizeBFECCBackwardSubstep() {};

    /**
     * @brief Initialize the BFECC final substep
     * This method is intended to implement all the operations required before each BFECC final substep
     */
    virtual void InitializeBFECCFinalSubstep() {};

    /**
     * @brief Finalize the BFECC final substep
     * This method is intended to implement all the operations required after each BFECC final substep
     */
    virtual void FinalizeBFECCFinalSubstep() {};

    /**
     * @brief Performs an intermediate RK4 step
     * This functions performs all the operations required in an intermediate RK4 sub step
     * @param rFixedDofsValues The vector containing the step n+1 values of the fixed DOFs
     * @param Substep The type of substep it is
     */
    virtual void PerformSubstep(const LocalSystemVectorType& rFixedDofsValues, const Substep SubstepType)
    {
        KRATOS_TRY

        // Get the required data from the explicit builder and solver
        auto& explicit_bs = BaseType::GetExplicitBuilder();
        auto& r_dof_set = explicit_bs.GetDofSet();
        const auto& r_lumped_mass_vector = explicit_bs.GetLumpedMassMatrixVector();

        // Get model part and information
        const double dt = BaseType::GetDeltaTime();
        KRATOS_ERROR_IF(dt < 1.0e-12) << "ProcessInfo DELTA_TIME is close to zero." << std::endl;
        auto& r_model_part = BaseType::GetModelPart();
        auto& r_process_info = r_model_part.GetProcessInfo();

        // Perform the intermidate sub step update

        double time_integration_theta = -1.0;
        double time_direction = 0.0;
        switch(SubstepType)
        {
            case FORWARD:
                InitializeBFECCForwardSubstep();
                time_integration_theta = 0.0;
                time_direction = 1.0;
                break;
            case BACKWARD:
                InitializeBFECCBackwardSubstep();
                time_integration_theta = 1.0;
                time_direction =-1.0;
                break;
            case FINAL:
                InitializeBFECCFinalSubstep();
                time_integration_theta = 0.0;
                time_direction = 1.0;
                break;
            default:
                KRATOS_ERROR << "Invalid value for Substep" << std::endl;
        }

        r_process_info.GetValue(TIME_INTEGRATION_THETA) = time_integration_theta;
        explicit_bs.BuildRHS(r_model_part);

        IndexPartition<SizeType>(r_dof_set.size()).for_each(
            [&](SizeType i_dof){
                auto it_dof = r_dof_set.begin() + i_dof;

                // Save current value in the corresponding vector
                const double residual = it_dof->GetSolutionStepReactionValue();

                // Do the DOF update
                double& r_u = it_dof->GetSolutionStepValue(0);
                double& r_u_prev_step = it_dof->GetSolutionStepValue(1);

                if (!it_dof->IsFixed()) {
                    const double mass = r_lumped_mass_vector[i_dof];
                    r_u += time_direction * (dt / mass) * residual;
                } else {
                    r_u = time_integration_theta*rFixedDofsValues[i_dof] + (1 - time_integration_theta)*r_u_prev_step;
                }
            }
        );

        switch(SubstepType) {
            case FORWARD:  FinalizeBFECCForwardSubstep();   break;
            case BACKWARD: FinalizeBFECCBackwardSubstep();  break;
            case FINAL:    FinalizeBFECCFinalSubstep();     break;
            default: KRATOS_ERROR << "Invalid value for Substep" << std::endl;
        }

        KRATOS_CATCH("Substep type: " + [](Substep substep) -> std::string {
            switch(substep) {
            case FORWARD:  return "FORWARD";
            case BACKWARD: return "BACKWARD";
            case FINAL:    return "FINAL";
            default: return std::to_string(static_cast<int>(substep));
            }
        }(SubstepType));
    }

    /**
     * Computes the error by comparing the buffer position 0 and the previous step
     * solution. Then stores the corrected starting point in the buffer position 1.
     * @param rPrevStepSolution
     */
    void CorrectErrorAfterForwardsAndBackwards(const LocalSystemVectorType& rPrevStepSolution)
    {
        auto& explicit_bs = BaseType::GetExplicitBuilder();
        auto& r_dof_set = explicit_bs.GetDofSet();

        IndexPartition<SizeType>(r_dof_set.size()).for_each(
            [&](const SizeType dof_index)
            {
                const double prev_step_solution = rPrevStepSolution[dof_index];
                auto& r_dof = *(r_dof_set.begin() + dof_index);

                const double error = (r_dof.GetSolutionStepValue(0) - prev_step_solution)/2.0;

                r_dof.GetSolutionStepValue(1) = prev_step_solution - error;

            }
        );
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

///@}

} /* namespace Kratos.*/

#endif /* KRATOS_EXPLICIT_SOLVING_STRATEGY_BFECC  defined */
