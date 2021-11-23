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

#if !defined(KRATOS_COMPRESSIBLE_NAVIER_STOKES_EXPLICIT_SOLVING_STRATEGY_RUNGE_KUTTA)
#define KRATOS_COMPRESSIBLE_NAVIER_STOKES_EXPLICIT_SOLVING_STRATEGY_RUNGE_KUTTA

// System includes
#include <functional>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/global_variables.h"
#include "includes/model_part.h"
#include "solving_strategies/strategies/explicit_solving_strategy_runge_kutta.h"
#include "utilities/atomic_utilities.h"
#include "utilities/math_utils.h"
#include "utilities/parallel_utilities.h"

// Application includes
#include "fluid_dynamics_application_variables.h"
#include "custom_processes/shock_capturing_physics_based_process.h"
#include "custom_processes/shock_capturing_entropy_viscosity_process.h"

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
template <class TSparseSpace, class TDenseSpace, class TButcherTableau>
class CompressibleNavierStokesExplicitSolvingStrategyRungeKutta
: public ExplicitSolvingStrategyRungeKutta<TSparseSpace, TDenseSpace, TButcherTableau>
{
public:
    ///@name Type Definitions
    ///@{

    /// The base class definition
    typedef ExplicitSolvingStrategyRungeKutta<TSparseSpace, TDenseSpace, TButcherTableau> BaseType;

    /// The explicit builder and solver definition
    typedef typename BaseType::ExplicitBuilderType ExplicitBuilderType;

    /// The local vector definition
    typedef typename TDenseSpace::VectorType LocalSystemVectorType;

    /// Pointer definition of CompressibleNavierStokesExplicitSolvingStrategyRungeKutta
    KRATOS_CLASS_POINTER_DEFINITION(CompressibleNavierStokesExplicitSolvingStrategyRungeKutta);

    /// Local Flags
    KRATOS_DEFINE_LOCAL_FLAG(SHOCK_CAPTURING);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor. (with parameters)
     * @param rModelPart The model part of the problem
     * @param ThisParameters The configuration parameters
     */
    explicit CompressibleNavierStokesExplicitSolvingStrategyRungeKutta(
        ModelPart &rModelPart,
        Parameters ThisParameters)
        : BaseType(rModelPart)
    {
        KRATOS_TRY

        // Validate and assign defaults
        ThisParameters = this->ValidateAndAssignParameters(ThisParameters, this->GetDefaultParameters());
        this->AssignSettings(ThisParameters);

        KRATOS_CATCH("")
    }

    /**
     * @brief Default constructor.
     * @param rModelPart The model part to be computed
     * @param pExplicitBuilder The pointer to the explicit builder and solver
     * @param MoveMeshFlag The flag to set if the mesh is moved or not
     */
    explicit CompressibleNavierStokesExplicitSolvingStrategyRungeKutta(
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
    explicit CompressibleNavierStokesExplicitSolvingStrategyRungeKutta(
        ModelPart &rModelPart,
        bool MoveMeshFlag = false,
        int RebuildLevel = 0)
        : BaseType(rModelPart, MoveMeshFlag, RebuildLevel)
    {
    }

    /** Copy constructor.
     */
    CompressibleNavierStokesExplicitSolvingStrategyRungeKutta(const CompressibleNavierStokesExplicitSolvingStrategyRungeKutta &Other) = delete;

    /** Destructor.
     */
    virtual ~CompressibleNavierStokesExplicitSolvingStrategyRungeKutta() = default;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    int Check() override
    {
        int err_code = BaseType::Check();

        // Check that the process info already contains the OSS activation variable
        const auto& r_process_info = BaseType::GetModelPart().GetProcessInfo();
        KRATOS_ERROR_IF_NOT(r_process_info.Has(OSS_SWITCH)) << "OSS_SWITCH variable has not been set in the ProcessInfo container. Please set it in the strategy \'Initialize\'." << std::endl;

        // Shock capturing process check
        if (mpShockCapturingProcess) {
            mpShockCapturingProcess->Check();
        }

        return err_code;
    }

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     * @return The default parameters
     */
    Parameters GetDefaultParameters() const override
    {
        KRATOS_TRY

        Parameters default_parameters = Parameters(R"(
        {
            "name" : "compressible_navier_stokes_explicit_explicit_solving_strategy_runge_kutta_4",
            "rebuild_level" : 0,
            "move_mesh_flag": false,
            "calculate_non_conservative_magnitudes" : true,
            "shock_capturing_settings" : { }
        })");

        // Getting base class default parameters
        const Parameters base_default_parameters = BaseType::GetDefaultParameters();
        default_parameters.RecursivelyAddMissingParameters(base_default_parameters);
        return default_parameters;

        KRATOS_CATCH("")
    }

    /**
     * @brief Returns the name of the class as used in the settings (snake_case format)
     * @return The name of the class
     */
    static std::string Name()
    {
        return "compressible_navier_stokes_explicit_solving_strategy_runge_kutta_4";
    }

    /**
     * @brief This method assigns settings to member variables
     * @param ThisParameters Parameters that are assigned to the member variables
     */
    void AssignSettings(const Parameters ThisParameters) override
    {
        // Base class assign settings call
        BaseType::AssignSettings(ThisParameters);
        mCalculateNonConservativeMagnitudes = ThisParameters["calculate_non_conservative_magnitudes"].GetBool();

        SetUpShockCapturing(ThisParameters["shock_capturing_settings"]);

        if (mpShockCapturingProcess && !mCalculateNonConservativeMagnitudes) {
            KRATOS_WARNING("CompressibleNavierStokesExplicitSolvingStrategyRungeKutta4")
                << "\'shock_capturing\' requires \'calculate_non_conservative_magnitudes\' to be active. Activating it." << std::endl;
            mCalculateNonConservativeMagnitudes = true;
        }
    }

    void SetUpShockCapturing(Parameters ShockCapturingParameters)
    {
        KRATOS_TRY

        auto defaults = Parameters(R"(
            {
                "type" : "physics_based",
                "Parameters" : {
                    "model_part_name" : ""
                }
            }
        )");
        defaults["Parameters"]["model_part_name"].SetString(BaseType::GetModelPart().Name());

        ShockCapturingParameters.ValidateAndAssignDefaults(defaults);
        ShockCapturingParameters.RecursivelyAddMissingParameters(defaults);

        using ShockCapturingFactoryType = Process::UniquePointer (*)(ModelPart&, Parameters);

        const std::map<const std::string, ShockCapturingFactoryType> shock_capturing_factory_map
        {
            {"none"         , [](ModelPart& m, Parameters p) -> Process::UniquePointer {return nullptr;}},
            {"physics_based", [](ModelPart& m, Parameters p) -> Process::UniquePointer {return Kratos::make_unique<ShockCapturingPhysicsBasedProcess>(m, p);}},
            {"entropy_based", [](ModelPart& m, Parameters p) -> Process::UniquePointer {return Kratos::make_unique<ShockCapturingEntropyViscosityProcess>(m, p);}}
        };

        const auto sc_type = ShockCapturingParameters["type"].GetString();

        const auto it = shock_capturing_factory_map.find(sc_type);
        if(it != shock_capturing_factory_map.end())
        {
            mpShockCapturingProcess = it->second(BaseType::GetModelPart(), ShockCapturingParameters["Parameters"]);
        }
        else
        {
            std::stringstream msg;
            msg << "Provided shock capturing type \""<< sc_type <<"\" does not match any of the available ones.\n";
            msg << "Please choose one from the following list:\n";
            for(const auto& keyvaluepair: shock_capturing_factory_map)
            {
                msg <<" - " << keyvaluepair.first << "\n";
            }
            msg << std::endl;

            KRATOS_ERROR << msg.str();
        }

        KRATOS_CATCH("")
    }

    /**
     * @brief Initialization of member variables and prior operations
     * In this method we call the base strategy initialize and initialize the time derivatives
     * This is required to prevent OpenMP errors as the time derivatives are stored in the non-historical database
     */
    void Initialize() override
    {
        auto& r_model_part = BaseType::GetModelPart();
        const auto& r_process_info = r_model_part.GetProcessInfo();

        // Call the base RK4 finalize substep method
        BaseType::Initialize();

        // Initialize the postprocess non-historical variables
        block_for_each(r_model_part.Nodes(), [](Node<3>& rNode){
            rNode.SetValue(MACH, 0.0);
            rNode.SetValue(SOUND_VELOCITY, 0.0);
        });

        // If required, initialize the OSS projection variables
        if (r_process_info[OSS_SWITCH]) {
            block_for_each(r_model_part.Nodes(), [](Node<3>& rNode){
                rNode.SetValue(DENSITY_PROJECTION, 0.0);
                rNode.SetValue(TOTAL_ENERGY_PROJECTION, 0.0);
                rNode.SetValue(MOMENTUM_PROJECTION, ZeroVector(3));
            });
        }

        // If required, initialize the physics-based shock capturing variables
        if (mpShockCapturingProcess) {
            mpShockCapturingProcess->ExecuteInitialize();
        }
    }

    void InitializeSolutionStep() override
    {
        BaseType::InitializeSolutionStep();

        if (mCalculateNonConservativeMagnitudes) {
            CalculateNonConservativeMagnitudes();
        }

        if (mpShockCapturingProcess) {
            mpShockCapturingProcess->ExecuteInitializeSolutionStep();
        }
    }

    /**
     * @brief Finalize the Runge-Kutta step
     * In this method we calculate the final linearised time derivatives after the final update
     * These will be the time derivatives employed in the first RK4 sub step of the next time step
     */
    void FinalizeSolutionStep() override
    {
        // Call the base RK4 finalize substep method
        BaseType::FinalizeSolutionStep();

        // Apply the momentum slip condition
        if (mApplySlipCondition) {
            ApplySlipCondition();
        }

        // Postprocess the non-conservative magnitudes
        // This needs to be done before the shock capturing as it is based on these
        if (mCalculateNonConservativeMagnitudes) {
            CalculateNonConservativeMagnitudes();
        }

        // Perform the shock capturing detection and artificial values calculation
        // This needs to be done at the end of the step in order to include the future shock
        // capturing magnitudes in the next automatic dt calculation
        if (mpShockCapturingProcess) {
            mpShockCapturingProcess->ExecuteFinalizeSolutionStep();
        }
    }

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "CompressibleNavierStokesExplicitSolvingStrategyRungeKutta";
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

    void InitializeRungeKuttaIntermediateSubStep() override
    {
        // Call the base RK4 to perform the initialize intermediate RK sub step
        BaseType::InitializeRungeKuttaIntermediateSubStep();

        // Approximate the unknowns time derivatives with a FE scheme
        // These will be used in the next RK substep residual calculation to compute the subscales
        auto& r_model_part = BaseType::GetModelPart();
        auto& r_process_info = r_model_part.GetProcessInfo();

        // Calculate the Orthogonal SubsScales projections
        if (r_process_info[OSS_SWITCH]) {
            CalculateOrthogonalSubScalesProjection();
        }
    }

    void InitializeRungeKuttaLastSubStep() override
    {
        // Call the base RK4 to perform the initialize intermediate RK sub step
        BaseType::InitializeRungeKuttaLastSubStep();

        // Calculate the Orthogonal SubsScales projections
        auto& r_model_part = BaseType::GetModelPart();
        auto& r_process_info = r_model_part.GetProcessInfo();

        // Calculate the Orthogonal SubsScales projections
        if (r_process_info[OSS_SWITCH]) {
            CalculateOrthogonalSubScalesProjection();
        }
    }

    /**
     * @brief Finalize the Runge-Kutta intermediate substep
     * In this method we calculate the linearised time derivatives after the intemediate substep
     */
    void FinalizeRungeKuttaIntermediateSubStep() override
    {
        // Call the base RK4 finalize substep method
        BaseType::FinalizeRungeKuttaIntermediateSubStep();

        // Apply the momentum slip condition
        if (mApplySlipCondition) {
            ApplySlipCondition();
        }
    }

    void FinalizeRungeKuttaLastSubStep() override
    {
        // Call the base RK4 finalize substep method
        BaseType::FinalizeRungeKuttaLastSubStep();

        // Apply the momentum slip condition
        //TODO: THIS SHOULDN'T BE REQUIRED --> DOING IT AFTER THE FINAL UPDATE MUST BE ENOUGH
        if (mApplySlipCondition) {
            ApplySlipCondition();
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

    bool mApplySlipCondition = true;
    bool mCalculateNonConservativeMagnitudes = true;

    Process::UniquePointer mpShockCapturingProcess = nullptr;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    void CalculateOrthogonalSubScalesProjection()
    {
        // Get model part data
        auto& r_model_part = BaseType::GetModelPart();
        const int n_nodes = r_model_part.NumberOfNodes();
        const auto& r_process_info = r_model_part.GetProcessInfo();
        const unsigned int block_size = r_process_info[DOMAIN_SIZE] + 2;

        // Get the required data from the explicit builder and solver
        // The lumped mass vector will be used to get the NODAL_AREA for the residuals projection
        const auto p_explicit_bs = BaseType::pGetExplicitBuilder();
        const auto& r_lumped_mass_vector = p_explicit_bs->GetLumpedMassMatrixVector();

        // Initialize the projection values
        block_for_each(r_model_part.Nodes(), [](Node<3>& rNode){
            rNode.GetValue(DENSITY_PROJECTION) = 0.0;
            rNode.GetValue(MOMENTUM_PROJECTION) = ZeroVector(3);
            rNode.GetValue(TOTAL_ENERGY_PROJECTION) = 0.0;
        });

        // Calculate the residuals projection
        std::tuple<double, double, array_1d<double,3>> oss_proj_tls;
        block_for_each(r_model_part.Elements(), oss_proj_tls, [&](Element& rElement, std::tuple<double, double, array_1d<double,3>>& rOssProjTLS){
            double& rho_proj = std::get<0>(rOssProjTLS);
            double& tot_ener_proj = std::get<1>(rOssProjTLS);
            array_1d<double,3>& mom_proj = std::get<2>(rOssProjTLS);
            rElement.Calculate(DENSITY_PROJECTION, rho_proj, r_process_info);
            rElement.Calculate(MOMENTUM_PROJECTION, mom_proj, r_process_info);
            rElement.Calculate(TOTAL_ENERGY_PROJECTION, tot_ener_proj, r_process_info);
        });

        // Do thhe nodal weighting
        // Note that to avoid calculating the NODAL_AREA we took from the first DOF of the lumped mass vector
        IndexPartition<>(n_nodes).for_each([&](const ModelPart::SizeType iNode){
            auto it_node = r_model_part.NodesBegin() + iNode;
            const double nodal_area = r_lumped_mass_vector[iNode * block_size];
            it_node->GetValue(DENSITY_PROJECTION) /= nodal_area;
            it_node->GetValue(MOMENTUM_PROJECTION) /= nodal_area;
            it_node->GetValue(TOTAL_ENERGY_PROJECTION) /= nodal_area;
        });

    }

    /**
     * @brief Calculates the non-conservative magnitudes
     * This function calculates the non-conservative magnitudes from the obtained conservative ones
     * The set of non-conservative magnitudes (VELOCITY, PRESSURE and TEMPERATURE) are saved in the
     * historical database. The speed of sound (SOUND_VELOCITY) and Mach number (MACH) are also
     * computed and stored in the non-historical database.
     */
    void CalculateNonConservativeMagnitudes()
    {
        // Get fluid properties from first element in mesh
        // Note that in here we assume that the fluid is single-phase
        auto& r_model_part = BaseType::GetModelPart();
        const auto& r_prop = r_model_part.ElementsBegin()->GetProperties();
        const double c_v = r_prop.GetValue(SPECIFIC_HEAT); // Specific heat at constant volume
        const double gamma = r_prop.GetValue(HEAT_CAPACITY_RATIO); // Heat capacity ratio
        const double R = (gamma - 1.0) * c_v; // Ideal gas constant

        // Loop the nodes to calculate the non-conservative magnitudes
        array_1d<double,3> aux_vel;
        block_for_each(r_model_part.Nodes(), aux_vel,[&] (Node<3> &rNode, array_1d<double,3>&rVelocity) {
            const auto& r_mom = rNode.FastGetSolutionStepValue(MOMENTUM);
            const double& r_rho = rNode.FastGetSolutionStepValue(DENSITY);
            const double& r_tot_ener = rNode.FastGetSolutionStepValue(TOTAL_ENERGY);
            rVelocity = r_mom / r_rho;
            const double temp = (r_tot_ener / r_rho - 0.5 * inner_prod(rVelocity, rVelocity)) / c_v;
            const double sound_velocity = std::sqrt(gamma * R * temp);
            rNode.FastGetSolutionStepValue(VELOCITY) = rVelocity;
            rNode.FastGetSolutionStepValue(TEMPERATURE) = temp;
            rNode.FastGetSolutionStepValue(PRESSURE) = r_rho * R * temp;
            rNode.GetValue(SOUND_VELOCITY) = sound_velocity;
            rNode.GetValue(MACH) = norm_2(rVelocity) / sound_velocity;
        });
    }

    /**
     * @brief Appy the slip condition
     * This method substracts the normal projection of the momentum for all the nodes flagged as SLIP
     * The correction is computed as m_slip = m - (m·n) x n. It is intended to be called after each RK substep
     */
    void ApplySlipCondition()
    {
        auto &r_model_part = BaseType::GetModelPart();

        block_for_each(r_model_part.Nodes(), [](Node<3>& rNode){
            if (rNode.Is(SLIP)) {
                auto unit_normal = rNode.FastGetSolutionStepValue(NORMAL);
                unit_normal /= norm_2(unit_normal);
                auto& r_mom = rNode.FastGetSolutionStepValue(MOMENTUM);
                const double r_mom_n = inner_prod(r_mom, unit_normal);
                noalias(r_mom) -= r_mom_n * unit_normal;
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
}; /* Class CompressibleNavierStokesExplicitSolvingStrategyRungeKutta */

///@}

///@name Type Definitions
///@{

///@}

template<class TSparseSpace, class TDenseSpace>
using CompressibleNavierStokesExplicitSolvingStrategyRungeKutta4 = CompressibleNavierStokesExplicitSolvingStrategyRungeKutta<TSparseSpace, TDenseSpace, ButcherTableauRK4>;

template<class TSparseSpace, class TDenseSpace>
using CompressibleNavierStokesExplicitSolvingStrategyRungeKutta3TVD = CompressibleNavierStokesExplicitSolvingStrategyRungeKutta<TSparseSpace, TDenseSpace, ButcherTableauRK3TVD>;

} /* namespace Kratos.*/

#endif /* KRATOS_COMPRESSIBLE_NAVIER_STOKES_EXPLICIT_SOLVING_STRATEGY_RUNGE_KUTTA  defined */
