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
    Parameters GetDefaultParameters() const override
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
     * @brief Returns the name of the class as used in the settings (snake_case format)
     * @return The name of the class
     */
    static std::string Name()
    {
        std::stringstream s;
        s << "explicit_solving_strategy_bfecc"
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

        const auto u_fixed = CloneSolutionStepData();

        // Calculate the intermediate sub steps
        SubstepForward(u_fixed);
        SubstepBackward(u_fixed);

        // Final update
        CorrectErrorAfterForwardsAndBackwards();
        FinalSubstep(u_fixed);

        KRATOS_CATCH("");
    }

    /**
     * Set the previous step solution in the current buffer position
     * Note that we set the 0 position of the buffer to be equal to the values in step n (not n+1)
     * Additionally, we save in an auxiliary vector the value of the fixed DOFs, which is also taken from the previous time step
     */
    LocalSystemVectorType CloneSolutionStepData()
    {
        // Get the required data from the explicit builder and solver
        auto& explicit_bs = BaseType::pGetExplicitBuilder();
        auto& r_dof_set = explicit_bs.GetDofSet();
        const SizeType dof_size = explicit_bs.GetEquationSystemSize();

        LocalSystemVectorType u_fixed(dof_size);
        IndexPartition<SizeType>(r_dof_set.size()).for_each(
            [&](const SizeType i_dof){
                auto it_dof = r_dof_set.begin() + i_dof;
                double& r_u_0 = it_dof->GetSolutionStepValue(0);
                if (it_dof->IsFixed()) {
                    u_fixed(i_dof) = r_u_0;
                }
                r_u_0 = it_dof->GetSolutionStepValue(1);
            }
        );

        return u_fixed;
    }

    /**
     * @brief Initialize the BFECC intermediate substep
     * This method is intended to implement all the operations required before each BFECC intermediate substep
     */
    virtual void InitializeBFECCIntermediateSubStep(SizeType SubstepIndex) {};

    /**
     * @brief Finalize the BFECC intermediate substep
     * This method is intended to implement all the operations required after each BFECC intermediate substep
     */
    virtual void FinalizeBFECCIntermediateSubStep(SizeType SubstepIndex) {};

    /**
     * @brief Initialize the BFECC last substep
     * This method is intended to implement all the operations required before each BFECC last substep
     */
    virtual void InitializeBFECCLastSubStep(SizeType SubstepIndex) {};

    /**
     * @brief Finalize the BFECC last substep
     * This method is intended to implement all the operations required after each BFECC last substep
     */
    virtual void FinalizeBFECCLastSubStep(SizeType SubstepIndex) {};

    void SubstepForward(const LocalSystemVectorType& rFixedDofsValues)
    {
        return PerformSubStep(1, 0.0, rFixedDofsValues, 1.0, false);
    }

    void SubstepBackward(const LocalSystemVectorType& rFixedDofsValues)
    {
        return PerformSubStep(2, 1.0, rFixedDofsValues,-1.0, false);
    }

    void FinalSubstep(const LocalSystemVectorType& rFixedDofsValues)
    {
        return PerformSubStep(3, 0.0, rFixedDofsValues, 1.0, true);
    }

    /**
     * @brief Performs an intermediate RK4 step
     * This functions performs all the operations required in an intermediate RK4 sub step
     * @param SubStepIndex The sub step index
     * @param TimeIntegrationTheta The point in time-step to evaulate in
     * @param rFixedDofsValues The vector containing the step n+1 values of the fixed DOFs
     * @param TimeDirection Whether going forawrds or backwards
     * @param LastSubstep Whether this is the last substep or not
     */
    virtual void PerformSubStep(
        const IndexType SubstepIndex,
        const double TimeIntegrationTheta,
        const LocalSystemVectorType& rFixedDofsValues,
        const double TimeDirection,
        bool LastSubstep)
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
        r_process_info.GetValue(TIME_INTEGRATION_THETA) = TimeIntegrationTheta;

        if(LastSubstep) {
            InitializeBFECCLastSubStep(SubstepIndex);
        } else {
            InitializeBFECCIntermediateSubStep(SubstepIndex);
        }

        explicit_bs.BuildRHS(r_model_part);

        IndexPartition<int>(r_dof_set.size()).for_each(
            [&](int i_dof){
                auto it_dof = r_dof_set.begin() + i_dof;

                // Save current value in the corresponding vector
                const double residual = it_dof->GetSolutionStepReactionValue();

                // Do the DOF update
                double& r_u = it_dof->GetSolutionStepValue(0);
                double& r_u_prev_step = it_dof->GetSolutionStepValue(0);

                if (!it_dof->IsFixed()) {
                    const double mass = r_lumped_mass_vector[i_dof];
                    r_u += TimeDirection * (dt / mass) * residual;
                } else {
                    r_u = integration_theta*rFixedDofsValues[i_dof] +(1-integration_theta)*r_u_prev_step;
                }
            }
        );

        if(LastSubstep) {
            FinalizeBFECCLastSubStep(SubstepIndex);
        } else {
            FinalizeBFECCIntermediateSubStep(SubstepIndex);
        }

        KRATOS_CATCH("SubstepIndex = " + std::to_string(SubStepIndex));
    }

    void CorrectErrorAfterForwardsAndBackwards()
    {
        auto& explicit_bs = BaseType::GetExplicitBuilder();
        auto& r_dof_set = explicit_bs.GetDofSet();

        block_for_each(r_dof_set, [&](Dof<double>& r_dof)
        {
            const double old_value = r_dof.GetSolutionStepValue(1);
            double& solution_step_value = r_dof.GetSolutionStepValue(0);
            const double error = (solution_step_value - old_value)/2.0;
            solution_step_value = old_value - error;
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
