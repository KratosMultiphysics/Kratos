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

    // The explicit builder and solver definition
    typedef typename BaseType::ExplicitBuilderAndSolverType ExplicitBuilderAndSolverType;

    /// The local vector definition
    typedef typename TDenseSpace::VectorType LocalSystemVectorType;

    // typedef typename DofsArrayType::iterator DofIteratorType;

    // typedef typename DofsArrayType::const_iterator DofConstantIteratorType;

    // typedef ModelPart::NodesContainerType NodesArrayType;

    // typedef ModelPart::ElementsContainerType ElementsArrayType;

    // typedef ModelPart::ConditionsContainerType ConditionsArrayType;

    /** Counted pointer of ClassName */
    KRATOS_CLASS_POINTER_DEFINITION(ExplicitSolvingStrategyRungeKutta4);

    ///@}
    ///@name Life Cycle
    ///@{

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
     * @param pExplicitBuilderAndSolver The pointer to the explicit builder and solver
     * @param MoveMeshFlag The flag to set if the mesh is moved or not
     */
    explicit ExplicitSolvingStrategyRungeKutta4(
        ModelPart &rModelPart,
        typename ExplicitBuilderAndSolverType::Pointer pExplicitBuilderAndSolver,
        bool MoveMeshFlag = false,
        int RebuildLevel = 0)
        : BaseType(rModelPart, pExplicitBuilderAndSolver, MoveMeshFlag, RebuildLevel)
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

    /** Copy constructor.
     */
    ExplicitSolvingStrategyRungeKutta4(const ExplicitSolvingStrategyRungeKutta4 &Other) = delete;

    /** Destructor.
     */
    virtual ~ExplicitSolvingStrategyRungeKutta4() = default;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    // /**
    //  * @brief Initialization of member variables and prior operations
    //  */
    // virtual void Initialize()
    // {
    //     // Base class initialize
    //     BaseType::Initialize();

    //     // Initialize the RK4 auxiliary vectors
    //     InitializeRungeKuttaAuxiliaryVectors();
    // }

    // /**
    //  * @brief Performs all the required operations that should be done (for each step) before solving the solution step.
    //  * @details A member variable should be used as a flag to make sure this function is called only once per step.
    //  */
    // virtual void InitializeSolutionStep()
    // {
    //     // Base class initialize solution step
    //     BaseType::Initialize();

    //     // If the DOFs have been rebuilt, reinitialize the RK auxiliary vectors
    //     if (BaseType::GetRebuildLevel != 0) {
    //         mpExplicitBuilderAndSolver->SetResetDofSetFlag(true);
    //     }
    // }


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
        const auto p_explicit_bs = BaseType::pGetExplicitBuilderAndSolver();
        auto& r_dof_set = p_explicit_bs->GetDofSet();
        const unsigned int dof_size = p_explicit_bs->GetEquationSystemSize();
        const auto& r_lumped_mass_vector = p_explicit_bs->GetLumpedMassMatrixVector();

        // Set the auxiliary RK vectors
        LocalSystemVectorType rk_k1(dof_size);
        LocalSystemVectorType rk_k2(dof_size);
        LocalSystemVectorType rk_k3(dof_size);
        LocalSystemVectorType rk_k4(dof_size);

        // Perform the RK 4 update
        const double dt = BaseType::GetDeltaTime();
        KRATOS_ERROR_IF(dt < 1.0e-12) << "ProcessInfo DELTA_TIME is close to zero." << std::endl;
        auto& r_model_part = BaseType::GetModelPart();

        // Set the previous step solution in the current buffer position
#pragma omp parallel for firstprivate(dof_size)
        for (int i_dof = 0; i_dof < dof_size; ++i_dof) {
            auto it_dof = r_dof_set.begin() + i_dof;
            if (!it_dof->IsFixed()) {
                it_dof->GetSolutionStepValue(0) = it_dof->GetSolutionStepValue(1);
            }
        }

        // 1st RK step
        p_explicit_bs->BuildRHS(r_model_part);
#pragma omp parallel for firstprivate(dof_size)
        for (int i_dof = 0; i_dof < dof_size; ++i_dof) {
            auto it_dof = r_dof_set.begin() + i_dof;
            // Save current value in the corresponding vector
            const double& r_res = it_dof->GetSolutionStepReactionValue();
            rk_k1(i_dof) = r_res;
            // Update the current DOF values
            if (!it_dof->IsFixed()) {
                const double mass = r_lumped_mass_vector(i_dof);
                it_dof->GetSolutionStepValue(0) = it_dof->GetSolutionStepValue(1) + mA21 * (dt / mass) * r_res;
            }
        }

        // 2nd RK STEP
        p_explicit_bs->BuildRHS(r_model_part);
#pragma omp parallel for firstprivate(dof_size)
        for (int i_dof = 0; i_dof < dof_size; ++i_dof) {
            auto it_dof = r_dof_set.begin() + i_dof;
            // Save current value in the corresponding vector
            const double& r_res = it_dof->GetSolutionStepReactionValue();
            rk_k2(i_dof) = r_res;
            // Do the DOF update
            if (!it_dof->IsFixed()) {
                const double mass = r_lumped_mass_vector(i_dof);
                it_dof->GetSolutionStepValue(0) = it_dof->GetSolutionStepValue(1) + mA32 * (dt / mass) * r_res;
            }
        }

        // 3rd RK STEP
        p_explicit_bs->BuildRHS(r_model_part);
#pragma omp parallel for firstprivate(dof_size)
        for (int i_dof = 0; i_dof < dof_size; ++i_dof) {
            auto it_dof = r_dof_set.begin() + i_dof;
            // Save current value in the corresponding vector
            const double& r_res = it_dof->GetSolutionStepReactionValue();
            rk_k3(i_dof) = r_res;
            // Do the DOF update
            if (!it_dof->IsFixed()) {
                const double mass = r_lumped_mass_vector(i_dof);
                it_dof->GetSolutionStepValue(0) = it_dof->GetSolutionStepValue(1) + mA43 * (dt / mass) * r_res;
            }
        }

        // 4rd RK STEP
        p_explicit_bs->BuildRHS(r_model_part);
#pragma omp parallel for firstprivate(dof_size)
        for (int i_dof = 0; i_dof < dof_size; ++i_dof) {
            auto it_dof = r_dof_set.begin() + i_dof;
            // Save current value in the corresponding vector
            const double& r_res = it_dof->GetSolutionStepReactionValue();
            rk_k4(i_dof) = r_res;
        }

        // Do the final solution update
#pragma omp parallel for firstprivate(dof_size)
        for (int i_dof = 0; i_dof < dof_size; ++i_dof) {
            auto it_dof = r_dof_set.begin() + i_dof;
            // Do the DOF update
            if (!it_dof->IsFixed()) {
                const double mass = r_lumped_mass_vector(i_dof);
                it_dof->GetSolutionStepValue(0) = it_dof->GetSolutionStepValue(1) + (dt / mass) * (mB1 * rk_k1(i_dof) + mB2 * rk_k2(i_dof) + mB3 * rk_k3(i_dof) + mB4 * rk_k4(i_dof));
            }
        }
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
