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

/**
 * @class ExplicitSolvingStrategyRungeKutta4ConvectionDiffusion
 * @ingroup ConvectionDiffusionApplication
 * @brief This strategy adds the orthogonal subgrid projections computation to the base explicit runge kutta 4 integration method.
 * @details The orthogonal subgrid scale projections are computed before each Runge-Kutta 4 step,
 * and after the final update, before updating the dynamic subscales.
 * @author Riccardo Tosi
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
    typedef typename BaseType::ExplicitBuilderType ExplicitBuilderType;

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
     * @param pExplicitBuilder The pointer to the explicit builder and solver
     * @param MoveMeshFlag The flag to set if the mesh is moved or not
     */
    explicit ExplicitSolvingStrategyRungeKutta4ConvectionDiffusion(
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

    /**
     * @brief Initialization of variables.
     * @details In this method, we call the base strategy initialize and initialize the projection variable.
     * This is required to prevent OpenMP errors as the projection variable is stored in the non-historical database.
     */
    void Initialize() override
    {
        auto& r_model_part = BaseType::GetModelPart();
        const auto& r_process_info = r_model_part.GetProcessInfo();
        // Call the base method
        BaseType::Initialize();
        // If required, initialize the OSS projection variables
        if (r_process_info[OSS_SWITCH]) {
            for (auto& r_node : r_model_part.GetCommunicator().LocalMesh().Nodes()) {
                r_node.SetValue(UNKNOWN_PROJECTION, 0.0);
            }
        }
    }

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
     * @brief Initialize the Runge-Kutta substep.
     * @details Calculate the orthogonal subscale projections if required.
     */
    virtual void InitializeRungeKuttaIntermediateSubStep() override
    {
        BaseType::InitializeRungeKuttaIntermediateSubStep();
        auto& r_model_part = BaseType::GetModelPart();
        const auto& r_process_info = r_model_part.GetProcessInfo();
        if (r_process_info[OSS_SWITCH] == 1) {
            ExecuteOSSStep();
        }
    };

    /**
     * @brief Initialize the Runge-Kutta substep.
     * @details Calculate the orthogonal subscale projections if required.
     */
    virtual void InitializeRungeKuttaLastSubStep() override
    {
        BaseType::InitializeRungeKuttaLastSubStep();
        auto& r_model_part = BaseType::GetModelPart();
        const auto& r_process_info = r_model_part.GetProcessInfo();
        if (r_process_info[OSS_SWITCH] == 1) {
            ExecuteOSSStep();
        }
    };

    /**
     * @brief Finalize the Runge Kutta explicit solver step.
     * @details Calculate the orthogonal subscale projections if required.
     */
    virtual void FinalizeSolutionStep() override
    {
        auto& r_model_part = BaseType::GetModelPart();
        const auto& r_process_info = r_model_part.GetProcessInfo();
        if (r_process_info[OSS_SWITCH] == 1) {
            ExecuteOSSStep();
        }
        BaseType::FinalizeSolutionStep();
    };

    /**
     * @brief Execute the OSS step.
     */
    virtual void ExecuteOSSStep()
    {
        // Get the required data from the explicit builder and solver
        const auto p_explicit_bs = BaseType::pGetExplicitBuilder();
        const auto& r_lumped_mass_vector = p_explicit_bs->GetLumpedMassMatrixVector();

        // Get model part data
        auto& r_model_part = BaseType::GetModelPart();
        auto& r_process_info = r_model_part.GetProcessInfo();
        const auto n_nodes = r_model_part.NumberOfNodes();
        const int n_elem = r_model_part.NumberOfElements();

        ConvectionDiffusionSettings::Pointer p_settings = r_process_info[CONVECTION_DIFFUSION_SETTINGS];
        auto& r_settings = *p_settings;

        // Initialize the projection value
#pragma omp parallel for
        for (std::size_t i_node = 0; i_node < n_nodes; ++i_node) {
            auto it_node = r_model_part.NodesBegin() + i_node;
            it_node->GetValue(UNKNOWN_PROJECTION) = 0.0;
        }

        // Calculate the unknown projection
        double unknown_proj;
#pragma omp parallel for
        for (int i_elem = 0; i_elem < n_elem; ++i_elem) {
            auto it_elem = r_model_part.ElementsBegin() + i_elem;
            it_elem->Calculate(UNKNOWN_PROJECTION, unknown_proj, r_process_info);
        }
#pragma omp parallel for
        for (std::size_t i_node = 0; i_node < n_nodes; ++i_node) {
            auto it_node = r_model_part.NodesBegin() + i_node;
            const double mass = r_lumped_mass_vector(i_node);
            it_node->FastGetSolutionStepValue(r_settings.GetProjectionVariable()) = it_node->GetValue(UNKNOWN_PROJECTION) / mass;
        }
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
