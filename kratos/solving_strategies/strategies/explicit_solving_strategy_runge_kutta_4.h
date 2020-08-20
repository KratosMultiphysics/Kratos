//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//
//

#if !defined(KRATOS_EXPLICIT_SOLVING_STRATEGY_RUNGE_KUTTA_4)
#define KRATOS_EXPLICIT_SOLVING_STRATEGY_RUNGE_KUTTA_4

/* System includes */

/* External includes */

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/strategies/explicit_solving_strategy.h"

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

/** @brief Explicit solving strategy base class
 * @details This is the base class from which we will derive all the explicit strategies (FE, RK4, ...)
 */
template <class TSparseSpace, class TDenseSpace>
class ExplicitSolvingStrategyRungeKutta4 : public ExplicitSolvingStrategy<TSparseSpace, TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{

    // The base class definition
    typedef ExplicitSolvingStrategy<TSparseSpace, TDenseSpace> BaseType;

    /// The definition of the current class
    typedef ExplicitSolvingStrategyRungeKutta4<TSparseSpace, TDenseSpace> ClassType;

    // The explicit builder and solver definition
    typedef typename BaseType::ExplicitBuilderType ExplicitBuilderType;

    /// The DOF type
    typedef typename BaseType::DofType DofType;

    /// The local vector definition
    typedef typename TDenseSpace::VectorType LocalSystemVectorType;

    /** Counted pointer of ClassName */
    KRATOS_CLASS_POINTER_DEFINITION(ExplicitSolvingStrategyRungeKutta4);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor. (empty)
     */
    explicit ExplicitSolvingStrategyRungeKutta4()
        : BaseType()
    {
    }

    /**
     * @brief Default constructor. (with parameters)
     * @param rModelPart The model part of the problem
     * @param ThisParameters The configuration parameters
     */
    explicit ExplicitSolvingStrategyRungeKutta4(
        ModelPart &rModelPart,
        Parameters ThisParameters)
        : BaseType(rModelPart, ThisParameters)
    {
    }

    /**
     * @brief Default constructor.
     * @param rModelPart The model part to be computed
     * @param pExplicitBuilder The pointer to the explicit builder and solver
     * @param MoveMeshFlag The flag to set if the mesh is moved or not
     */
    explicit ExplicitSolvingStrategyRungeKutta4(
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
    explicit ExplicitSolvingStrategyRungeKutta4(
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
    typename BaseType::Pointer Create(
        ModelPart& rModelPart,
        Parameters ThisParameters
        ) const override
    {
        return Kratos::make_shared<ClassType>(rModelPart, ThisParameters);
    }

    /** Copy constructor.
     */
    ExplicitSolvingStrategyRungeKutta4(const ExplicitSolvingStrategyRungeKutta4 &Other) = delete;

    /** Destructor.
     */
    ~ExplicitSolvingStrategyRungeKutta4() override = default;

    /**
     * @brief Returns the name of the class as used in the settings (snake_case format)
     * @return The name of the class
     */
    static std::string Name()
    {
        return "explicit_solving_strategy_runge_kutta_4";
    }

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "ExplicitSolvingStrategyRungeKutta4";
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
        // Get the required data from the explicit builder and solver
        const auto p_explicit_bs = BaseType::pGetExplicitBuilder();
        auto& r_dof_set = p_explicit_bs->GetDofSet();
        const unsigned int dof_size = p_explicit_bs->GetEquationSystemSize();
        const auto& r_lumped_mass_vector = p_explicit_bs->GetLumpedMassMatrixVector();

        // Set the auxiliary RK vectors
        LocalSystemVectorType u_n(dof_size); // TODO: THIS IS INEFICCIENT. CREATE A UNORDERED_SET WITH THE IDOF AND VALUE AS ENTRIES. THIS HAS TO BE OPTIONAL
        LocalSystemVectorType rk_k1(dof_size);
        LocalSystemVectorType rk_k2(dof_size);
        LocalSystemVectorType rk_k3(dof_size);
        LocalSystemVectorType rk_k4(dof_size);

        // Perform the RK 4 update
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
                    u_n(i_dof) = r_u_1;
                }
                r_u_0 = r_u_1;
            }
        );

        // Calculate the RK4 intermediate sub steps
        PerformRungeKuttaIntermediateSubStep(1, mA21, u_n, rk_k1);
        PerformRungeKuttaIntermediateSubStep(2, mA32, u_n, rk_k2);
        PerformRungeKuttaIntermediateSubStep(3, mA43, u_n, rk_k3);
        PerformRungeKuttaLastSubStep(rk_k4);

        // Do the final solution update
        IndexPartition<int>(r_dof_set.size()).for_each(
            [&](int i_dof){
                auto it_dof = r_dof_set.begin() + i_dof;
                // Do the DOF update
                double& r_u = it_dof->GetSolutionStepValue(0);
                const double& r_u_old = it_dof->GetSolutionStepValue(1);
                if (!it_dof->IsFixed()) {
                    const double mass = r_lumped_mass_vector(i_dof);
                    r_u = r_u_old + (dt / mass) * (mB1 * rk_k1(i_dof) + mB2 * rk_k2(i_dof) + mB3 * rk_k3(i_dof) + mB4 * rk_k4(i_dof));
                } else {
                    r_u = u_n(i_dof);
                }
            }
        );
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
     * @param SubStepCoefficient The sub step coefficient (these are saved as member variables)
     * @param rFixedDofsValues The vector containing the step n+1 values of the fixed DOFs
     * @param rIntermediateStepResidualVector The vector to store the intermediate sub step residual
     */
    virtual void PerformRungeKuttaIntermediateSubStep(
        const IndexType SubStepIndex,
        const double SubStepCoefficient,
        const LocalSystemVectorType& rFixedDofsValues,
        LocalSystemVectorType& rIntermediateStepResidualVector)
    {
        // Get the required data from the explicit builder and solver
        const auto p_explicit_bs = BaseType::pGetExplicitBuilder();
        auto& r_dof_set = p_explicit_bs->GetDofSet();
        const auto& r_lumped_mass_vector = p_explicit_bs->GetLumpedMassMatrixVector();

        // Get model part and information
        const double dt = BaseType::GetDeltaTime();
        KRATOS_ERROR_IF(dt < 1.0e-12) << "ProcessInfo DELTA_TIME is close to zero." << std::endl;
        auto& r_model_part = BaseType::GetModelPart();
        auto& r_process_info = r_model_part.GetProcessInfo();

        // Set the RUNGE_KUTTA_STEP value. This has to be done prior to the InitializeRungeKuttaStep()
        r_process_info.GetValue(RUNGE_KUTTA_STEP) = SubStepIndex;

        // Perform the intermidate sub step update
        InitializeRungeKuttaIntermediateSubStep();
        p_explicit_bs->BuildRHS(r_model_part);

        IndexPartition<int>(r_dof_set.size()).for_each(
            [&](int i_dof){
                auto it_dof = r_dof_set.begin() + i_dof;
                // Save current value in the corresponding vector
                const double& r_res = it_dof->GetSolutionStepReactionValue();
                rIntermediateStepResidualVector(i_dof) = r_res;
                // Do the DOF update
                double& r_u = it_dof->GetSolutionStepValue(0);
                const double& r_u_old = it_dof->GetSolutionStepValue(1);
                if (!it_dof->IsFixed()) {
                    const double mass = r_lumped_mass_vector(i_dof);
                    r_u = r_u_old + SubStepCoefficient * (dt / mass) * r_res;
                } else {
                    const double delta_u = rFixedDofsValues(i_dof) - r_u_old;
                    r_u = r_u_old + SubStepCoefficient * delta_u;
                }
            }
        );

        FinalizeRungeKuttaIntermediateSubStep();
    }

    /**
     * @brief Performs the last RK4 step
     * This functions performs all the operations required in the last RK4 sub step
     * @param rLastStepResidualVector The vector to store the last sub step residual
     */
    virtual void PerformRungeKuttaLastSubStep(LocalSystemVectorType& rLastStepResidualVector)
    {
        // Get the required data from the explicit builder and solver
        const auto p_explicit_bs = BaseType::pGetExplicitBuilder();
        auto& r_dof_set = p_explicit_bs->GetDofSet();

        // Get model part
        auto& r_model_part = BaseType::GetModelPart();
        auto& r_process_info = r_model_part.GetProcessInfo();

        // Set the RUNGE_KUTTA_STEP value. This has to be done prior to the InitializeRungeKuttaStep()
        r_process_info.GetValue(RUNGE_KUTTA_STEP) = 4;

        // Perform the last sub step residual calculation
        InitializeRungeKuttaLastSubStep();
        p_explicit_bs->BuildRHS(r_model_part);

        IndexPartition<int>(r_dof_set.size()).for_each(
            [&](int i_dof){
                const auto it_dof = r_dof_set.begin() + i_dof;
                // Save current value in the corresponding vector
                const double& r_res = it_dof->GetSolutionStepReactionValue();
                rLastStepResidualVector(i_dof) = r_res;
            }
        );

        FinalizeRungeKuttaLastSubStep();
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

    const double mA21 = 0.5; // RK4 a_21 coefficient
    const double mA32 = 0.5; // RK4 a_32 coefficient
    const double mA43 = 1.0; // RK4 a_43 coefficient
    const double mB1 = 1.0/6.0; // RK4 b_1 coefficient
    const double mB2 = 1.0/3.0; // RK4 b_2 coefficient
    const double mB3 = 1.0/3.0; // RK4 b_3 coefficient
    const double mB4 = 1.0/6.0; // RK4 b_4 coefficient

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
}; /* Class NewExplicitSolvingStrategyRungeKutta4 */

///@}

///@name Type Definitions
///@{

///@}

} /* namespace Kratos.*/

#endif /* KRATOS_EXPLICIT_SOLVING_STRATEGY_RUNGE_KUTTA_4  defined */
