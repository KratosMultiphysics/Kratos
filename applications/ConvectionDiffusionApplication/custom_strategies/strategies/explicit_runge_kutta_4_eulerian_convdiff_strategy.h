//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Tosi
//
//

#if !defined(KRATOS_EXPLICIT_RUNGE_KUTTA_4_EULERIAN_CONVDIFF_STRATEGY)
#define KRATOS_EXPLICIT_RUNGE_KUTTA_4_EULERIAN_CONVDIFF_STRATEGY

/* System includes */

/* External includes */

/* Project includes */
#include "solving_strategies/strategies/explicit_solving_strategy_runge_kutta_4.h"
#include "includes/convection_diffusion_settings.h"

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

/** @brief Explicit solving strategy for convection diffusion explicit eulerian element
 *  @details This strategy adds the OSS step
 */
template <class TSparseSpace, class TDenseSpace>
class ExplicitSolvingStrategyRungeKutta4ConvectionDiffusion : public ExplicitSolvingStrategyRungeKutta4<TSparseSpace, TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{

    // The base class definition
    typedef ExplicitSolvingStrategyRungeKutta4<TSparseSpace, TDenseSpace> BaseType;

    // The explicit builder and solver definition
    typedef typename BaseType::ExplicitBuilderAndSolverType ExplicitBuilderAndSolverType;

    /// The local vector definition
    typedef typename TDenseSpace::VectorType LocalSystemVectorType;

    /** Counted pointer of ClassName */
    KRATOS_CLASS_POINTER_DEFINITION(ExplicitSolvingStrategyRungeKutta4ConvectionDiffusion);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor. (with parameters)
     * @param rModelPart The model part of the problem
     * @param ThisParameters The configuration parameters
     */
    explicit ExplicitSolvingStrategyRungeKutta4ConvectionDiffusion(
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
    explicit ExplicitSolvingStrategyRungeKutta4ConvectionDiffusion(
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
    explicit ExplicitSolvingStrategyRungeKutta4ConvectionDiffusion(
        ModelPart &rModelPart,
        bool MoveMeshFlag = false,
        int RebuildLevel = 0)
        : BaseType(rModelPart, MoveMeshFlag, RebuildLevel)
    {
    }

    /** Copy constructor.
     */
    ExplicitSolvingStrategyRungeKutta4ConvectionDiffusion(const ExplicitSolvingStrategyRungeKutta4ConvectionDiffusion &Other) = delete;

    /** Destructor.
     */
    virtual ~ExplicitSolvingStrategyRungeKutta4ConvectionDiffusion() = default;

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
        return "ExplicitSolvingStrategyRungeKutta4ConvectionDiffusion";
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

    /**
     * @brief Initialize the Runge-Kutta substep
     * The Orthogonal Subgrid Scale step is run
     */
    virtual void InitializeRungeKuttaSubStep() override
    {
        // Get the required data from the explicit builder and solver
        const auto p_explicit_bs = BaseType::pGetExplicitBuilderAndSolver();
        auto& r_dof_set = p_explicit_bs->GetDofSet();
        const auto& r_lumped_mass_vector = p_explicit_bs->GetLumpedMassMatrixVector();

        // Perform Orthogonal Subgrid Scale step if USE_OSS is active
        auto& r_model_part = BaseType::GetModelPart();
        auto& r_process_info = r_model_part.GetProcessInfo();
        if (r_process_info.GetValue(USE_OSS) == 1)
        {
            // Activate OSS flag used inside the element
            r_process_info.GetValue(OSS_SWITCH) = 1;
            p_explicit_bs->BuildRHS(r_model_part);

            ConvectionDiffusionSettings::Pointer p_settings = r_process_info[CONVECTION_DIFFUSION_SETTINGS];
            auto& r_settings = *p_settings;
            const auto number_of_nodes = r_model_part.NumberOfNodes();
#pragma omp parallel for firstprivate(number_of_nodes)
            for (unsigned int i_node = 0; i_node < number_of_nodes; i_node++)
            {
                auto& current_node = r_model_part.GetNode(i_node+1);
                const double mass = r_lumped_mass_vector(i_node);
                current_node.FastGetSolutionStepValue(r_settings.GetProjectionVariable()) = current_node.FastGetSolutionStepValue(r_settings.GetReactionVariable()) / mass;
            }
            // End OSS step
        }
        // Deactivate OSS flag used inside the element
        r_process_info.GetValue(OSS_SWITCH) = 0;
    };

    /**
     * @brief Initialize the Runge-Kutta substep
     * The Orthogonal Subgrid Scale step is run
     */
    virtual void FinalizeRungeKuttaStep() override
    {
        // Get the required data from the explicit builder and solver
        const auto p_explicit_bs = BaseType::pGetExplicitBuilderAndSolver();
        auto& r_dof_set = p_explicit_bs->GetDofSet();
        const auto& r_lumped_mass_vector = p_explicit_bs->GetLumpedMassMatrixVector();

        // Perform Orthogonal Subgrid Scale step if USE_OSS is active
        auto& r_model_part = BaseType::GetModelPart();
        auto& r_process_info = r_model_part.GetProcessInfo();
        if (r_process_info.GetValue(USE_OSS) == 1)
        {
            if (r_process_info.GetValue(RUNGE_KUTTA_STEP) == 4)
            {
                // Activate OSS flag used inside the element
                r_process_info.GetValue(OSS_SWITCH) = 1;
                p_explicit_bs->BuildRHS(r_model_part);

                ConvectionDiffusionSettings::Pointer p_settings = r_process_info[CONVECTION_DIFFUSION_SETTINGS];
                auto& r_settings = *p_settings;
                const auto number_of_nodes = r_model_part.NumberOfNodes();
#pragma omp parallel for firstprivate(number_of_nodes)
                for (unsigned int i_node = 0; i_node < number_of_nodes; i_node++)
                {
                    auto& current_node = r_model_part.GetNode(i_node+1);
                    const double mass = r_lumped_mass_vector(i_node);
                    current_node.FastGetSolutionStepValue(r_settings.GetProjectionVariable()) = current_node.FastGetSolutionStepValue(r_settings.GetReactionVariable()) / mass;
                }
                // End OSS step
            }
        }
        // Deactivate OSS flag used inside the element
        r_process_info.GetValue(OSS_SWITCH) = 0;
    };

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
}; /* Class ExplicitSolvingStrategyRungeKutta4ConvectionDiffusion */

///@}

///@name Type Definitions
///@{

///@}

} /* namespace Kratos.*/

#endif /* KRATOS_EXPLICIT_RUNGE_KUTTA_4_EULERIAN_CONVDIFF_STRATEGY  defined */
