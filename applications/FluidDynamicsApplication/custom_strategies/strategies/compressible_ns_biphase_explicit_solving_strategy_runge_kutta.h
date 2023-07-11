//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla (adapted by Andrea Montanino)
//
//

#if !defined(KRATOS_COMPRESSIBLE_NS_BIPHASE_EXPLICIT_SOLVING_STRATEGY_RUNGE_KUTTA)
#define KRATOS_COMPRESSIBLE_NS_BIPHASE_EXPLICIT_SOLVING_STRATEGY_RUNGE_KUTTA

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
#include "processes/find_nodal_neighbours_process.h" 
#include "processes/find_global_nodal_neighbours_process.h"

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
class CompressibleNSBiphaseExplicitSolvingStrategyRungeKutta
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

    /// Pointer definition of CompressibleNSBiphaseExplicitSolvingStrategyRungeKutta
    KRATOS_CLASS_POINTER_DEFINITION(CompressibleNSBiphaseExplicitSolvingStrategyRungeKutta);

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
    explicit CompressibleNSBiphaseExplicitSolvingStrategyRungeKutta(
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
    explicit CompressibleNSBiphaseExplicitSolvingStrategyRungeKutta(
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
    explicit CompressibleNSBiphaseExplicitSolvingStrategyRungeKutta(
        ModelPart &rModelPart,
        bool MoveMeshFlag = false,
        int RebuildLevel = 0)
        : BaseType(rModelPart, MoveMeshFlag, RebuildLevel)
    {
    }

    /** Copy constructor.
     */
    CompressibleNSBiphaseExplicitSolvingStrategyRungeKutta(const CompressibleNSBiphaseExplicitSolvingStrategyRungeKutta &Other) = delete;

    /** Destructor.
     */
    virtual ~CompressibleNSBiphaseExplicitSolvingStrategyRungeKutta() = default;

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
            "name" : "compressible_ns_biphase_explicit_solving_strategy_runge_kutta_4",
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
        return "compressible_navier_stokes_explicit_solving_strategy_runge_kutta_4";  // CheckAndreaMontanino
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
            KRATOS_WARNING("CompressibleNSBiphaseExplicitSolvingStrategyRungeKutta4")
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

        FindNodalNeighboursProcess process(r_model_part);
            process.Execute();

        // Call the base RK4 finalize substep method
        BaseType::Initialize();

        // Initialize the postprocess non-historical variables
        block_for_each(r_model_part.Nodes(), [](Node<3>& rNode){
            rNode.SetValue(MACH, 0.0);
            rNode.SetValue(SOUND_VELOCITY, 0.0);            
            rNode.SetValue(DYNAMIC_PRESSURE, 0.0);
            rNode.SetValue(SOLID_CONCENTRATION, 0.0);
            rNode.SetValue(GAS_PRESSURE, 0.0);
        });
        // If required, initialize the OSS projection variables
        if (r_process_info[OSS_SWITCH]) {
            block_for_each(r_model_part.Nodes(), [](Node<3>& rNode){
                rNode.SetValue(DENSITY_PROJECTION, 0.0);
                rNode.SetValue(DENSITY_SOLID_PROJECTION, 0.0);   // Added by Andrea
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
        return "CompressibleNSBiphaseExplicitSolvingStrategyRungeKutta";
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
        /*
        if (r_process_info[OSS_SWITCH]) {
            CalculateOrthogonalSubScalesProjection();
        }
        */
    }

    void InitializeRungeKuttaLastSubStep() override
    {
        // Call the base RK4 to perform the initialize intermediate RK sub step
        BaseType::InitializeRungeKuttaLastSubStep();

        // Calculate the Orthogonal SubsScales projections
        auto& r_model_part = BaseType::GetModelPart();
        auto& r_process_info = r_model_part.GetProcessInfo();

        // Calculate the Orthogonal SubsScales projections
        /*
        if (r_process_info[OSS_SWITCH]) {
            CalculateOrthogonalSubScalesProjection();
        }
        */
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
/*
    void CalculateOrthogonalSubScalesProjection()
    {
        // Get model part data
        auto& r_model_part = BaseType::GetModelPart();
        const int n_nodes = r_model_part.NumberOfNodes();
        const auto& r_process_info = r_model_part.GetProcessInfo();
        const unsigned int block_size = r_process_info[DOMAIN_SIZE] + 3;     // From +2 to +3 (because it is biphase)

        // Get the required data from the explicit builder and solver
        // The lumped mass vector will be used to get the NODAL_AREA for the residuals projection
        const auto p_explicit_bs = BaseType::pGetExplicitBuilder();
        const auto& r_lumped_mass_vector = p_explicit_bs->GetLumpedMassMatrixVector();

        // Initialize the projection values
        block_for_each(r_model_part.Nodes(), [](Node<3>& rNode){
            rNode.GetValue(DENSITY_PROJECTION) = 0.0;
            //rNode.GetValue(DENSITY_SOLID_PROJECTION) = 0.0;
            rNode.GetValue(MOMENTUM_PROJECTION) = ZeroVector(3);
            rNode.GetValue(TOTAL_ENERGY_PROJECTION) = 0.0;
        });

        // Calculate the residuals projection
        std::tuple<double, double, array_1d<double,3>> oss_proj_tls;
        block_for_each(r_model_part.Elements(), oss_proj_tls, [&](Element& rElement, std::tuple<double, double, array_1d<double,3>>& rOssProjTLS){
            double& rho_proj = std::get<0>(rOssProjTLS);
            double& rho_solid_proj = std::get<1>(rOssProjTLS);
            double& tot_ener_proj = std::get<2>(rOssProjTLS);
            array_1d<double,3>& mom_proj = std::get<3>(rOssProjTLS);   //Check this order
            rElement.Calculate(DENSITY_PROJECTION, rho_proj, r_process_info);
            //rElement.Calculate(DENSITY_SOLID_PROJECTION, rho_solid_proj, r_process_info);
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
*/
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
        double tol = 1e-16;
        
        auto& r_model_part = BaseType::GetModelPart();
        const auto& r_prop = r_model_part.ElementsBegin()->GetProperties();
        const double c_v = r_prop.GetValue(SPECIFIC_HEAT); // Specific heat at constant volume
        const double gamma = r_prop.GetValue(HEAT_CAPACITY_RATIO); // Heat capacity ratio
        const double R = (gamma - 1.0) * c_v; // Ideal gas constant
        const double rho_s = r_prop.GetValue(SOLID_MATERIAL_DENSITY); // Density of the solid part
        const double c_s = r_prop.GetValue(SOLID_MATERIAL_SPECIFIC_HEAT); // Specific heat of the solid part

        // Cleaning DS

        double N1 = 0, N3 = 0;

        double M1 = 0, M3 = 1;

        while (M3 > tol){

            M1 = 0;
            M3 = 0;
            N1 = 0;
            N3 = 0;

            #pragma omp parallel for
            for (int i = 0; i < static_cast<int>(r_model_part.NumberOfNodes()); ++i)
            {   
                auto it_node    = r_model_part.NodesBegin() + i;

                double denS     = it_node->FastGetSolutionStepValue(DENSITY_SOLID);
                
                double A = it_node->FastGetSolutionStepValue(NODAL_AREA);
                
                if (denS > 0)    {
                    M1 += denS*A;
                    N1 += A;
                }
                if (denS < 0)    {
                    M3 += fabs(denS)*A;
                    N3 += A;
                }
            }

            #pragma omp parallel for
            for (int i = 0; i < static_cast<int>(r_model_part.NumberOfNodes()); ++i)
            {   
                auto it_node    = r_model_part.NodesBegin() + i;

                double denS = it_node->FastGetSolutionStepValue(DENSITY_SOLID); 
                
                if (denS > 0)   denS -= M3/N1;
                if (denS < 0)   denS  = 0.0;

                it_node->FastGetSolutionStepValue(DENSITY_SOLID) = denS;

            }

            printf("M3 = %.3e \n", M3);

        }

// Cleaning ETOT

        N1 = 0;
        N3 = 0;

        M1 = 0;
        M3 = 0;
        tol = 1;

        while (M3 > tol){

            M1 = 0;
            M3 = 0;
            N1 = 0;
            N3 = 0;

            #pragma omp parallel for
            for (int i = 0; i < static_cast<int>(r_model_part.NumberOfNodes()); ++i)
            {   
                auto it_node    = r_model_part.NodesBegin() + i;

                double denE     = it_node->FastGetSolutionStepValue(TOTAL_ENERGY);
                
                double A = it_node->FastGetSolutionStepValue(NODAL_AREA);
                
                if (denE > 0)    {
                    M1 += denE*A;
                    N1 += A;
                }
                if (denE < 0)    {
                    M3 += fabs(denE)*A;
                    N3 += A;
                }
            }

            #pragma omp parallel for
            for (int i = 0; i < static_cast<int>(r_model_part.NumberOfNodes()); ++i)
            {   
                auto it_node    = r_model_part.NodesBegin() + i;

                double denE = it_node->FastGetSolutionStepValue(TOTAL_ENERGY); 
                
                if (denE > 0)   denE -= M3/N1;
                if (denE < 0)   denE  = -0.1*tol;

                it_node->FastGetSolutionStepValue(TOTAL_ENERGY) = denE;

            }

            printf("M3 = %.3e \n", M3);
            M3 = 0;

        }        

// Smoothing phase         
        
        
        double  alpha_dt = 1*4e-3;
        double  alpha_ds = 1*4e-3;
        double  alpha_dm = 1*8e-4;
        double  alpha_de = 1*8e-3;


        int check = 1;
        int count = 0;

        while (check == 1 && count < 1)
        {
            
            check = 0;

            #pragma omp parallel for
            for (int i = 0; i < static_cast<int>(r_model_part.NumberOfNodes()); ++i)
            {
                auto it_node    = r_model_part.NodesBegin() + i;

                double dt_i = it_node->FastGetSolutionStepValue(DENSITY);
                double ds_i = it_node->FastGetSolutionStepValue(DENSITY_SOLID);
                array_1d<double,3> mom_i = it_node->FastGetSolutionStepValue(MOMENTUM);
                double ene_i = it_node->FastGetSolutionStepValue(TOTAL_ENERGY);

                double iArea     = it_node->FastGetSolutionStepValue(NODAL_AREA);

                double iAreaTOT  = iArea;
                
                double dt_corr   = 0.0;
                double ds_corr   = 0.0;
                array_1d<double,3> dm_corr;
                dm_corr[0] = 0.0;
                dm_corr[1] = 0.0;
                dm_corr[2] = 0.0;
                double ene_corr  = 0.0;
                
                GlobalPointersVector< Node<3> >& rneigh = it_node->GetValue(NEIGHBOUR_NODES);
                for( GlobalPointersVector<Node<3> >::iterator jt_node = rneigh.begin(); jt_node!=rneigh.end(); jt_node++){
                    
                    double jArea    = jt_node->FastGetSolutionStepValue(NODAL_AREA);
                    
                    double dt_j     = jt_node->FastGetSolutionStepValue(DENSITY);
                    double ds_j     = jt_node->FastGetSolutionStepValue(DENSITY_SOLID);
                    array_1d<double,3> mom_j = jt_node->FastGetSolutionStepValue(MOMENTUM);
                    double ene_j    = jt_node->FastGetSolutionStepValue(TOTAL_ENERGY);
                    

                    double Areaij   = 0.5*(iArea + jArea);

                    iAreaTOT += Areaij;

                    dt_corr  += (dt_i - dt_j)*Areaij;
                    ds_corr  += (ds_i - ds_j)*Areaij;
                    dm_corr  += (mom_i - mom_j)*Areaij; 
                    ene_corr += (ene_i - ene_j)*Areaij;
                    
                }

                it_node->FastGetSolutionStepValue(DENSITY_RHS)       =  dt_i - alpha_dt*dt_corr/iAreaTOT;
                it_node->FastGetSolutionStepValue(DENSITY_SOLID_RHS) =  ds_i - alpha_ds*ds_corr/iAreaTOT;
                it_node->FastGetSolutionStepValue(MOMENTUM_RHS)      =  mom_i - alpha_dm*dm_corr/iAreaTOT;
                it_node->FastGetSolutionStepValue(TOTAL_ENERGY_RHS)  =  ene_i - alpha_de*ene_corr/iAreaTOT;
                
 //               if(fabs( - alpha_de*ene_corr/iAreaTOT) > 1e-3)
 //                   printf("corr_ene = %.3e\n",  - alpha_de*ene_corr/iAreaTOT);


            }

            if (check == 1) {
                count++;
                printf("count = %d\n", count);
            }

            printf("alpha = %.3e \n", alpha_de);
            alpha_de *= 2;
            alpha_ds *= 2;

        }
        
        // Loop the nodes to calculate the non-conservative magnitudes
        array_1d<double,3> aux_vel;
        block_for_each(r_model_part.Nodes(), aux_vel,[&] (Node<3> &rNode, array_1d<double,3>&rVelocity) {
            
            rNode.FastGetSolutionStepValue(MOMENTUM)        = rNode.FastGetSolutionStepValue(MOMENTUM_RHS);
            rNode.FastGetSolutionStepValue(DENSITY)         = rNode.FastGetSolutionStepValue(DENSITY_RHS);
            rNode.FastGetSolutionStepValue(DENSITY_SOLID)   = rNode.FastGetSolutionStepValue(DENSITY_SOLID_RHS);
            rNode.FastGetSolutionStepValue(TOTAL_ENERGY)    = rNode.FastGetSolutionStepValue(TOTAL_ENERGY_RHS);
            

            const auto& r_mom = rNode.FastGetSolutionStepValue(MOMENTUM);
            const double& r_rho = rNode.FastGetSolutionStepValue(DENSITY);
            const double& r_rho_solid = rNode.FastGetSolutionStepValue(DENSITY_SOLID);
            const double& r_tot_ener = rNode.FastGetSolutionStepValue(TOTAL_ENERGY);
        
            rVelocity = r_mom / r_rho;
            double c_mixed = c_v*(r_rho - r_rho_solid) + c_s*r_rho_solid;
            double Cp_mixed = (c_v + R)*(r_rho - r_rho_solid) + c_s*r_rho_solid;
            double sol_conc = r_rho_solid/rho_s;
            
            const double temp = (r_tot_ener - 0.5 * r_rho * inner_prod(rVelocity, rVelocity)) / c_mixed;   // Modify for biphase flow
            const double pressure = (r_rho - r_rho_solid) * R * temp;
            const double sound_velocity = std::sqrt(r_rho - r_rho_solid)*R*Cp_mixed*(2*r_tot_ener - r_rho*inner_prod(rVelocity, rVelocity))/(2*r_rho*c_mixed*c_mixed); 
            rNode.FastGetSolutionStepValue(VELOCITY) = rVelocity;
            rNode.FastGetSolutionStepValue(SOLID_CONCENTRATION) = sol_conc; 
            rNode.FastGetSolutionStepValue(TEMPERATURE) = temp;
            rNode.FastGetSolutionStepValue(PRESSURE) = pressure;
            rNode.FastGetSolutionStepValue(GAS_PRESSURE) = pressure/(1 - sol_conc); 
            rNode.FastGetSolutionStepValue(SOUND_VELOCITY) = sound_velocity;
            rNode.FastGetSolutionStepValue(MACH) = norm_2(rVelocity) / sound_velocity;
            rNode.FastGetSolutionStepValue(DYNAMIC_PRESSURE) = 0.5 * r_rho * inner_prod(rVelocity, rVelocity);
        
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
}; /* Class CompressibleNSBiphaseExplicitSolvingStrategyRungeKutta */

///@}

///@name Type Definitions
///@{

///@}

template<class TSparseSpace, class TDenseSpace>
using CompressibleNSBiphaseExplicitSolvingStrategyRungeKutta4 = CompressibleNSBiphaseExplicitSolvingStrategyRungeKutta<TSparseSpace, TDenseSpace, ButcherTableauRK4>;

template<class TSparseSpace, class TDenseSpace>
using CompressibleNSBiphaseExplicitSolvingStrategyRungeKutta3TVD = CompressibleNSBiphaseExplicitSolvingStrategyRungeKutta<TSparseSpace, TDenseSpace, ButcherTableauRK3TVD>;

template<class TSparseSpace, class TDenseSpace>
using CompressibleNSBiphaseExplicitSolvingStrategyRungeKuttaForwardEuler = CompressibleNSBiphaseExplicitSolvingStrategyRungeKutta<TSparseSpace, TDenseSpace, ButcherTableauForwardEuler>;

} /* namespace Kratos.*/

#endif /* KRATOS_COMPRESSIBLE_NS_BIPHASE_EXPLICIT_SOLVING_STRATEGY_RUNGE_KUTTA  defined */

    