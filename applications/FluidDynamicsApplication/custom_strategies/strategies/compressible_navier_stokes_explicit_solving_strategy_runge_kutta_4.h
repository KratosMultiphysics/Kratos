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

#if !defined(KRATOS_COMPRESSIBLE_NAVIER_STOKES_EXPLICIT_SOLVING_STRATEGY_RUNGE_KUTTA_4)
#define KRATOS_COMPRESSIBLE_NAVIER_STOKES_EXPLICIT_SOLVING_STRATEGY_RUNGE_KUTTA_4

// System includes
#include <functional>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/global_variables.h"
#include "includes/model_part.h"
#include "solving_strategies/strategies/explicit_solving_strategy_runge_kutta_4.h"
#include "utilities/math_utils.h"

// Application includes
#include "fluid_dynamics_application_variables.h"

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
class CompressibleNavierStokesExplicitSolvingStrategyRungeKutta4 : public ExplicitSolvingStrategyRungeKutta4<TSparseSpace, TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{

    /// The base class definition
    typedef ExplicitSolvingStrategyRungeKutta4<TSparseSpace, TDenseSpace> BaseType;

    /// The explicit builder and solver definition
    typedef typename BaseType::ExplicitBuilderType ExplicitBuilderType;

    /// The local vector definition
    typedef typename TDenseSpace::VectorType LocalSystemVectorType;

    /// Pointer definition of CompressibleNavierStokesExplicitSolvingStrategyRungeKutta4
    KRATOS_CLASS_POINTER_DEFINITION(CompressibleNavierStokesExplicitSolvingStrategyRungeKutta4);

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
    explicit CompressibleNavierStokesExplicitSolvingStrategyRungeKutta4(
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
    explicit CompressibleNavierStokesExplicitSolvingStrategyRungeKutta4(
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
    explicit CompressibleNavierStokesExplicitSolvingStrategyRungeKutta4(
        ModelPart &rModelPart,
        bool MoveMeshFlag = false,
        int RebuildLevel = 0)
        : BaseType(rModelPart, MoveMeshFlag, RebuildLevel)
    {
    }

    /** Copy constructor.
     */
    CompressibleNavierStokesExplicitSolvingStrategyRungeKutta4(const CompressibleNavierStokesExplicitSolvingStrategyRungeKutta4 &Other) = delete;

    /** Destructor.
     */
    virtual ~CompressibleNavierStokesExplicitSolvingStrategyRungeKutta4() = default;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     * @return The default parameters
     */
    Parameters GetDefaultParameters() const override
    {
        Parameters default_parameters = Parameters(R"(
        {
            "name" : "compressible_navier_stokes_explicit_explicit_solving_strategy_runge_kutta_4",
            "rebuild_level" : 0,
            "move_mesh_flag": false,
            "shock_capturing" : true
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

        // Set the specific compressible NS settings
        mShockCapturing = ThisParameters["shock_capturing"].GetBool();
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

        // Initialize the non-historical database values
        for (auto& r_node : r_model_part.GetCommunicator().LocalMesh().Nodes()) {
            // Initialize the unknowns time derivatives to zero
            r_node.SetValue(DENSITY_TIME_DERIVATIVE, 0.0);
            r_node.SetValue(MOMENTUM_TIME_DERIVATIVE, MOMENTUM_TIME_DERIVATIVE.Zero());
            r_node.SetValue(TOTAL_ENERGY_TIME_DERIVATIVE, 0.0);
        }

        // If required, initialize the OSS projection variables
        if (r_process_info[OSS_SWITCH]) {
            for (auto& r_node : r_model_part.GetCommunicator().LocalMesh().Nodes()) {
                r_node.SetValue(NODAL_AREA, 0.0);
                r_node.SetValue(DENSITY_PROJECTION, 0.0);
                r_node.SetValue(TOTAL_ENERGY_PROJECTION, 0.0);
                r_node.SetValue(MOMENTUM_PROJECTION, ZeroVector(3));
            }
        }

        // If required, initialize the orthogonal projection shock capturing variables
        if (mShockCapturing) {
            // Initialize nodal values
            for (auto& r_node : r_model_part.GetCommunicator().LocalMesh().Nodes()) {
                r_node.SetValue(NODAL_AREA, 0.0);
                r_node.SetValue(ARTIFICIAL_CONDUCTIVITY, 0.0);
                r_node.SetValue(ARTIFICIAL_BULK_VISCOSITY, 0.0);
                r_node.SetValue(ARTIFICIAL_DYNAMIC_VISCOSITY, 0.0);
                r_node.SetValue(DENSITY_GRADIENT, ZeroVector(3));
            }

            // Initialize elemental values
            for (auto& r_elem : r_model_part.GetCommunicator().LocalMesh().Elements()) {
                r_elem.SetValue(DENSITY_GRADIENT, ZeroVector(3));
            }
        }
    }

    /**
     * @brief Initialize the Runge-Kutta step
     * In case the mesh has been updated in the previous step we need to reinitialize the shock capturing
     * This includes the calculation of the nodal element size, nodal area and nodal neighbours
     */
    void InitializeSolutionStep() override
    {
        // Call the base RK4 initialize substep method
        BaseType::InitializeSolutionStep();

        // Calculate the magnitudes time derivatives
        UpdateUnknownsTimeDerivatives(1.0);
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

        // Perform the shock capturing detection and artificial values calculation
        // This needs to be done at the end of the step in order to include the future shock 
        // capturing magnitudes in the next automatic dt calculation
        if (mShockCapturing) {
            CalculatePhysicsBasedShockCapturing();
        }
    }

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "CompressibleNavierStokesExplicitSolvingStrategyRungeKutta4";
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
        const auto& r_process_info = r_model_part.GetProcessInfo();

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

    bool mShockCapturing = true;
    bool mApplySlipCondition = true;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief Update the compressible Navier-Stokes unknowns time derivatives
     * This method approximates the compressible Navier-Stokes unknowns time derivatives
     * These are required to calculate the inertial stabilization terms in the compressible NS element
     * To that purpose a linear Forward-Euler interpolation is used
     */
    void UpdateUnknownsTimeDerivatives(const double SubStepAccCoefficient)
    {
        const double dt = BaseType::GetDeltaTime();
        KRATOS_ERROR_IF(dt < 1.0e-12) << "ProcessInfo DELTA_TIME is close to zero." << std::endl;
        auto &r_model_part = BaseType::GetModelPart();

#pragma omp parallel for
        for (int i_node = 0; i_node < static_cast<int>(r_model_part.NumberOfNodes()); ++i_node) {
            auto it_node = r_model_part.NodesBegin() + i_node;

            // Density DOF time derivative
            const double& r_rho = it_node->FastGetSolutionStepValue(DENSITY);
            const double& r_rho_old = it_node->FastGetSolutionStepValue(DENSITY, 1);
            it_node->GetValue(DENSITY_TIME_DERIVATIVE) = SubStepAccCoefficient * (r_rho - r_rho_old) / dt;

            // Momentum DOF time derivative
            const auto& r_mom = it_node->FastGetSolutionStepValue(MOMENTUM);
            const auto& r_mom_old = it_node->FastGetSolutionStepValue(MOMENTUM, 1);
            it_node->GetValue(MOMENTUM_TIME_DERIVATIVE) = SubStepAccCoefficient * (r_mom - r_mom_old) / dt;

            // Total energy DOF time derivative
            const double& r_tot_enr = it_node->FastGetSolutionStepValue(TOTAL_ENERGY);
            const double& r_tot_enr_old = it_node->FastGetSolutionStepValue(TOTAL_ENERGY, 1);
            it_node->GetValue(TOTAL_ENERGY_TIME_DERIVATIVE) = SubStepAccCoefficient * (r_tot_enr - r_tot_enr_old) / dt;
        }
    }

    void CalculateOrthogonalSubScalesProjection()
    {
        // Get model part data
        auto& r_model_part = BaseType::GetModelPart();
        const int n_nodes = r_model_part.NumberOfNodes();
        const int n_elem = r_model_part.NumberOfElements();
        const auto& r_process_info = r_model_part.GetProcessInfo();

        // Get the required data from the explicit builder and solver
        const auto p_explicit_bs = BaseType::pGetExplicitBuilder();
        // TODO_ USE THE LUMPED MASS VECTOR AS NODAL_AREA
        // TODO_ USE THE LUMPED MASS VECTOR AS NODAL_AREA
        // TODO_ USE THE LUMPED MASS VECTOR AS NODAL_AREA
        // TODO_ USE THE LUMPED MASS VECTOR AS NODAL_AREA
        const auto& r_lumped_mass_vector = p_explicit_bs->GetLumpedMassMatrixVector();

        // Initialize the projection values
#pragma omp parallel for
        for (int i_node = 0; i_node < n_nodes; ++i_node) {
            auto it_node = r_model_part.NodesBegin() + i_node;
            it_node->GetValue(NODAL_AREA) = 0.0;
            it_node->GetValue(DENSITY_PROJECTION) = 0.0;
            it_node->GetValue(MOMENTUM_PROJECTION) = ZeroVector(3);
            it_node->GetValue(TOTAL_ENERGY_PROJECTION) = 0.0;
        }

        // Calculate the residuals projection
        double dens_proj;
        double tot_ener_proj;
        array_1d<double,3> mom_proj;
#pragma omp parallel for
        for (int i_elem = 0; i_elem < n_elem; ++i_elem) {
            auto it_elem = r_model_part.ElementsBegin() + i_elem;
            // Calculate the projection values
            it_elem->Calculate(DENSITY_PROJECTION, dens_proj, r_process_info);
            it_elem->Calculate(MOMENTUM_PROJECTION, mom_proj, r_process_info);
            it_elem->Calculate(TOTAL_ENERGY_PROJECTION, tot_ener_proj, r_process_info);
            // Calculate the NODAL_AREA
            // TODO: This is not probably required each time
            auto& r_geom = it_elem->GetGeometry();
            const unsigned int n_nodes = r_geom.PointsNumber();
            const double geom_domain_size = r_geom.DomainSize();
            const double aux_weight = geom_domain_size / static_cast<double>(n_nodes);
            for (auto& r_node : r_geom) {
#pragma omp atomic
                r_node.GetValue(NODAL_AREA) += aux_weight;
            }
        }

#pragma omp parallel for
        for (int i_node = 0; i_node < n_nodes; ++i_node) {
            auto it_node = r_model_part.NodesBegin() + i_node;
            const double nodal_area = it_node->GetValue(NODAL_AREA);
            it_node->GetValue(DENSITY_PROJECTION) /= nodal_area;
            it_node->GetValue(MOMENTUM_PROJECTION) /= nodal_area;
            it_node->GetValue(TOTAL_ENERGY_PROJECTION) /= nodal_area;
        }
    }

    void CalculatePhysicsBasedShockCapturing()
    {
        // Calculate the model part data
        auto& r_model_part = BaseType::GetModelPart();
        const auto r_process_info = r_model_part.GetProcessInfo();
        const int n_nodes = r_model_part.NumberOfNodes();
        const int n_elems = r_model_part.NumberOfElements();

        // Initialize the values to zero
#pragma omp parallel for
        for (int i_node = 0; i_node < n_nodes; ++i_node) {
            auto it_node = r_model_part.NodesBegin() + i_node;
            it_node->GetValue(NODAL_AREA) = 0.0;
            it_node->GetValue(DENSITY_GRADIENT) = ZeroVector(3);
            it_node->GetValue(ARTIFICIAL_CONDUCTIVITY) = 0.0;
            it_node->GetValue(ARTIFICIAL_BULK_VISCOSITY) = 0.0;
            it_node->GetValue(ARTIFICIAL_DYNAMIC_VISCOSITY) = 0.0;
        }

        // Set the functor to calculate the element size
        // Note that this assumes a unique geometry in the computational mesh
        std::function<std::tuple<double, double, Matrix>(Geometry<Node<3>>& rGeometry)> elem_metric_function;
        const GeometryData::KratosGeometryType geometry_type = (r_model_part.ElementsBegin()->GetGeometry()).GetGeometryType();
        switch (geometry_type) {
            case GeometryData::KratosGeometryType::Kratos_Triangle2D3:
                elem_metric_function = [&](const Geometry<Node<3>>& rGeometry){return CalculateTriangleMetricTensor(rGeometry);};
                break;
            case GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4:
                elem_metric_function = [&](const Geometry<Node<3>>& rGeometry){return CalculateTetrahedraMetricTensor(rGeometry);};
                break;
            default:
                KRATOS_ERROR << "Asking for a non-implemented geometry.";
        }

        // Loop the elements to project the gradients
        // Note that it is assumed that the gradient is constant within the element
        // Hence, only one Gauss point is used
        const double eps = 1.0e-7;

        double div_v;
        Matrix grad_vel;
        array_1d<double,3> rot_v;
        array_1d<double,3> grad_rho;
        array_1d<double,3> grad_temp;
        Geometry<Node<3>>::ShapeFunctionsGradientsType dNdX_container;
#pragma omp parallel for private(dNdX_container, div_v, rot_v, grad_rho, grad_temp, grad_vel)
        for (int i_elem = 0; i_elem < n_elems; ++i_elem) {
            auto it_elem = r_model_part.ElementsBegin() + i_elem;
            auto& r_geom = it_elem->GetGeometry();
            const unsigned int n_nodes = r_geom.PointsNumber();

            // Get fluid properties
            const auto p_prop = it_elem->pGetProperties();
            const double c_v = p_prop->GetValue(SPECIFIC_HEAT);
            const double gamma = p_prop->GetValue(HEAT_CAPACITY_RATIO);

            // Calculate elemental magnitudes
            const double k = 1.0; // Polynomial order of the numerical simulation
            double c_ref; // Reference speed of sound (Ma = 1.0)
            // TODO: CALLING THE CALCULATES IS NOT THE MOST EFFICIENT WAY... THINK ABOUT THIS...
            it_elem->Calculate(SOUND_VELOCITY, c_ref, r_process_info);             
            it_elem->Calculate(DENSITY_GRADIENT, grad_rho, r_process_info);
            it_elem->Calculate(VELOCITY_DIVERGENCE, div_v, r_process_info);
            it_elem->Calculate(VELOCITY_ROTATIONAL, rot_v, r_process_info);
            it_elem->Calculate(VELOCITY_GRADIENT, grad_vel, r_process_info);
            it_elem->Calculate(TEMPERATURE_GRADIENT, grad_temp, r_process_info);

            // Calculate midpoint values
            double midpoint_rho = 0.0;
            double midpoint_tot_ener = 0.0;
            array_1d<double, 3> midpoint_v = ZeroVector(3);
            for (unsigned int i_node = 0; i_node < n_nodes; ++i_node) {
                const double& r_rho = r_geom[i_node].FastGetSolutionStepValue(DENSITY);
                const double& r_tot_ener = r_geom[i_node].FastGetSolutionStepValue(TOTAL_ENERGY);
                midpoint_rho += r_rho;
                midpoint_tot_ener += r_tot_ener;
                midpoint_v += r_geom[i_node].FastGetSolutionStepValue(MOMENTUM) / r_rho;
            }
            midpoint_rho /= static_cast<double>(n_nodes);
            midpoint_v /= static_cast<double>(n_nodes);
            const double v_norm_pow = midpoint_v[0]*midpoint_v[0] + midpoint_v[1]*midpoint_v[1] + midpoint_v[2]*midpoint_v[2];

            // Inverse metric tensor calculation
            // const double h_ref = avg_h_function(r_geom); // Reference element size used in the metric tensor
            // const auto metric_tensor = elem_metric_function(r_geom, h_ref); // Metric tensor relative to the reference element size
            const auto metric_data = elem_metric_function(r_geom);
            const double h_ref = std::get<0>(metric_data); // Reference element size used in the metric tensor
            const double metric_tensor_inf = std::get<1>(metric_data); // Metric tensor infimum norm (smallest eigenvalue)
            const auto metric_tensor = std::get<2>(metric_data); // Metric tensor relative to the reference element size
            double aux_det;
            Matrix inv_metric_tensor;
            MathUtils<double>::InvertMatrix(metric_tensor, inv_metric_tensor, aux_det);

            // Characteristic element sizes
            array_1d<double,3> inv_metric_grad_rho = ZeroVector(3);
            array_1d<double,3> inv_metric_grad_temp = ZeroVector(3);
            for (unsigned int i = 0; i < inv_metric_tensor.size1(); ++i) {
                for (unsigned int j = 0; j < inv_metric_tensor.size2(); ++j) {
                    inv_metric_grad_rho(i) += inv_metric_tensor(i,j) * grad_rho(j);
                    inv_metric_grad_temp(i) += inv_metric_tensor(i,j) * grad_temp(j);
                }
            }
            const double h_beta = h_ref * norm_2(grad_rho) / std::sqrt(inner_prod(grad_rho, inv_metric_grad_rho) + eps); // Characteristic element size along the direction of the density gradient
            const double h_kappa = h_ref * norm_2(grad_temp) / std::sqrt(inner_prod(grad_temp, inv_metric_grad_temp) + eps); // Characteristic element size along the direction of the temperature gradient
            const double h_mu = h_ref * metric_tensor_inf;

            // Dilatation sensor (activates in shock waves)
            const double s_omega = - h_beta * div_v / k / c_ref;

            // Vorticity sensor (vanishes in vorticity dominated regions)
            const double div_v_pow = std::pow(div_v, 2);
            const double rot_v_norm_pow = rot_v[0] * rot_v[0] + rot_v[1] * rot_v[1] + rot_v[2] * rot_v[2];
            const double s_w = div_v_pow / (div_v_pow + rot_v_norm_pow + eps);

            // Calculate limited shock sensor
            const double s_beta_0 = 0.01;
            const double s_beta_max = 2.0 / std::sqrt(std::pow(gamma, 2) - 1.0);
            const double s_beta = s_omega * s_w;
            // const double s_beta_hat = LimitingFunction(s_beta, s_beta_0, s_beta_max);
            const double s_beta_hat = SmoothedLimitingFunction(s_beta, s_beta_0, s_beta_max);
            it_elem->GetValue(SHOCK_SENSOR) = s_beta_hat;

            // Thermal sensor (detect thermal gradients that are larger than possible with the grid resolution)
            Matrix mid_pt_jacobian;
            r_geom.Jacobian(mid_pt_jacobian, 0, GeometryData::GI_GAUSS_1);
            array_1d<double,3> local_grad_temp = ZeroVector(3);
            for (unsigned int i = 0; i < mid_pt_jacobian.size1(); ++i) {
                for (unsigned int j = 0; j < mid_pt_jacobian.size2(); ++j) {
                    local_grad_temp(i) += mid_pt_jacobian(j,i) * grad_temp(j);
                }
            }
            const double stagnation_temp = midpoint_tot_ener / midpoint_rho / c_v;

            const double s_kappa_0 = 1.0;
            const double s_kappa_max = 2.0;
            const double s_kappa = h_ref * norm_2(local_grad_temp) / k / stagnation_temp;
            const double s_kappa_hat = SmoothedLimitingFunction(s_kappa, s_kappa_0, s_kappa_max);
            it_elem->GetValue(THERMAL_SENSOR) = s_kappa_hat;

            // Shear sensor (detect velocity gradients that are larger than possible with the grid resolution)
            const unsigned int dim = r_geom.WorkingSpaceDimension();
            Matrix shear_grad_vel(dim, dim);
            for (unsigned int d1 = 0; d1 < dim; ++d1) {
                for (unsigned int d2 = 0; d2 < dim; ++d2) {
                    shear_grad_vel(d1, d2) = d1 == d2 ? 0.0 : grad_vel(d1, d2);
                }
            }
            const Matrix local_shear_grad_vel = prod(shear_grad_vel, Matrix(trans(mid_pt_jacobian)));
            Matrix eigen_vect_mat, eigen_val_mat;
            MathUtils<double>::GaussSeidelEigenSystem(local_shear_grad_vel, eigen_vect_mat, eigen_val_mat);
            double shear_spect_norm = 0.0; 
            for (unsigned int d = 0; d < eigen_val_mat.size1(); ++d){
                if (eigen_val_mat(d,d) > shear_spect_norm) {
                    shear_spect_norm = eigen_val_mat(d,d);
                }
            }
            const double isentropic_max_vel = std::sqrt(v_norm_pow + (2.0 / (gamma - 1.0)) *  std::pow(c_ref, 2)); // TODO: THIS ISN'T c_ref ACCORDING TO PERAIRE 

            const double s_mu_0 = 1.0;
            const double s_mu_max = 2.0;
            const double s_mu = h_ref * shear_spect_norm / isentropic_max_vel / k;
            // const double s_mu_hat = LimitingFunction(s_mu, s_mu_0, s_mu_max);
            const double s_mu_hat = SmoothedLimitingFunction(s_mu, s_mu_0, s_mu_max);
            it_elem->GetValue(SHEAR_SENSOR) = s_mu_hat;

            // Calculate artificial magnitudes
            const double ref_mom_norm = midpoint_rho * std::sqrt(v_norm_pow + std::pow(c_ref,2));

            // Calculate elemental artificial bulk viscosity
            const double k_beta = 1.5;
            const double elem_b_star =  (k_beta * h_beta / k) * ref_mom_norm * s_beta_hat; 
            it_elem->GetValue(ARTIFICIAL_BULK_VISCOSITY) = elem_b_star;

            // Calculate elemental artificial conductivity (dilatancy)
            const double Pr_beta_min = 0.9;
            const double alpha_pr_beta = 2.0;
            const double Mach_threshold = 3.0;
            const double Mach = norm_2(midpoint_v) / c_ref;
            const double Pr_beta = Pr_beta_min * (1.0 + std::exp(-2.0 * alpha_pr_beta * (Mach - Mach_threshold)));
            const double elem_k1_star = (gamma * c_v / Pr_beta) * elem_b_star;

            // Calculate elemental artificial conductivity (thermal sensor)
            const double k_kappa = 1.0;
            const double elem_k2_star = (gamma * c_v) * (k_kappa * h_kappa / k) * ref_mom_norm * s_kappa_hat;
            it_elem->GetValue(ARTIFICIAL_CONDUCTIVITY) = elem_k1_star + elem_k2_star;

            // Calculate elemental artificial dynamic viscosity
            const double k_mu = 1.0;
            const double elem_mu_star = (k_mu * h_mu / k) * ref_mom_norm * s_mu_hat;
            it_elem->GetValue(ARTIFICIAL_DYNAMIC_VISCOSITY) = elem_mu_star;

            // Project the shock capturing magnitudes to the nodes
            const double geom_domain_size = r_geom.DomainSize();
            const double aux_weight = geom_domain_size / static_cast<double>(n_nodes);
            for (unsigned int i_node = 0; i_node < n_nodes; ++i_node) {
                auto& r_node = r_geom[i_node];
#pragma omp atomic
                r_node.GetValue(ARTIFICIAL_BULK_VISCOSITY) += aux_weight * elem_b_star;
#pragma omp atomic
                r_node.GetValue(ARTIFICIAL_CONDUCTIVITY) += aux_weight * (elem_k1_star + elem_k2_star);
#pragma omp atomic
                r_node.GetValue(ARTIFICIAL_DYNAMIC_VISCOSITY) += aux_weight * elem_mu_star;
#pragma omp atomic
                r_node.GetValue(NODAL_AREA) += aux_weight;
            }
        }

        // Nodal smoothing of the shock capturing magnitudes
#pragma omp parallel for
        for (int i_node = 0; i_node < n_nodes; ++i_node) {
            auto it_node = r_model_part.NodesBegin() + i_node;
            const double weight = it_node->GetValue(NODAL_AREA);
            it_node->GetValue(ARTIFICIAL_CONDUCTIVITY) /= weight;
            it_node->GetValue(ARTIFICIAL_BULK_VISCOSITY) /= weight;
            it_node->GetValue(ARTIFICIAL_DYNAMIC_VISCOSITY) /= weight;
        }
    }

    double LimitingFunction(
        const double s,
        const double s_0,
        const double s_max,
        const double s_min = 0.0)
    {
        const double aux_1 = std::max(s - s_0, s_min);
        const double aux_2 = std::min(aux_1 - s_max, s_min);
        return aux_2 + s_max;
    }

    double SmoothedLimitingFunction(
        const double s,
        const double s_0,
        const double s_max)
    {
        const double aux_1 = SmoothedMaxFunction(s - s_0);
        const double aux_2 = SmoothedMinFunction(aux_1 - s_max);
        return aux_2 + s_max;
    }

    // Smooth approximation of the max(s,0) function
    double SmoothedMaxFunction(const double s)
    {
        const double b = 100;
        const double l_max = (s / Globals::Pi) * std::atan(b * s) + 0.5 * s - (1.0 / Globals::Pi) * std::atan(b) + 0.5;
        return l_max;
    }

    // Smooth approximation of the min(s,0) function
    double SmoothedMinFunction(const double s)
    {
        return s - SmoothedMaxFunction(s);
    }

    // https://es.wikipedia.org/wiki/Circunelipse_de_Steiner
    std::tuple<double, double, Matrix> CalculateTriangleMetricTensor(const Geometry<Node<3>>& rGeometry)
    {
        const array_1d<double, 3> p_1 = rGeometry[0].Coordinates();
        const array_1d<double, 3> p_2 = rGeometry[1].Coordinates();
        const array_1d<double, 3> p_3 = rGeometry[2].Coordinates();

        // Solve the metric problem trans(e)*M*e = 1
        // This means, find the coefficients of the matrix M such that all the edges have unit length
        Vector sol;
        array_1d<double,3> aux_vect;
        BoundedMatrix<double,3,3> aux_mat;
        aux_mat(0,0) = std::pow(p_1[0]-p_2[0], 2); aux_mat(0,1) = 2.0*(p_1[0]-p_2[0])*(p_1[1]-p_2[1]); aux_mat(0,2) = std::pow(p_1[1]-p_2[1], 2);
        aux_mat(1,0) = std::pow(p_1[0]-p_3[0], 2); aux_mat(1,1) = 2.0*(p_1[0]-p_3[0])*(p_1[1]-p_3[1]); aux_mat(1,2) = std::pow(p_1[1]-p_3[1], 2);
        aux_mat(2,0) = std::pow(p_2[0]-p_3[0], 2); aux_mat(2,1) = 2.0*(p_2[0]-p_3[0])*(p_2[1]-p_3[1]); aux_mat(2,2) = std::pow(p_2[1]-p_3[1], 2);
        aux_vect[0] = 1.0;
        aux_vect[1] = 1.0;
        aux_vect[2] = 1.0;        
        MathUtils<double>::Solve(aux_mat, sol, aux_vect);

        // Set the metric tensor
        Matrix metric(2,2);
        metric(0,0) = sol[0]; metric(0,1) = sol[1];
        metric(1,0) = sol[1]; metric(1,1) = sol[2];

        // Calculate the eigenvalues of the metric tensor to obtain the ellipsis of inertia axes lengths
        BoundedMatrix<double,2,2> eigenvects, eigenvals;
        MathUtils<double>::GaussSeidelEigenSystem(metric, eigenvects, eigenvals);
        const double h_1 = std::sqrt(1.0 / eigenvals(0,0));
        const double h_2 = std::sqrt(1.0 / eigenvals(1,1));

        // Calculate the reference element size as the average of the ellipsis of intertia axes lengths
        const double h_ref = 0.5 * (h_1 + h_2);
        
        // Make the metric dimensionless
        metric *= std::pow(h_ref,2);

        // Calculate metric infimum norm
        const double metric_inf = std::min(eigenvals(0,0), eigenvals(1,1));

        return std::make_tuple(h_ref, metric_inf, metric);
    }

    // https://es.wikipedia.org/wiki/Circunelipse_de_Steiner --> 3D extension
    std::tuple<double, double, Matrix> CalculateTetrahedraMetricTensor(const Geometry<Node<3>>& rGeometry)
    {
        // Solve the metric problem trans(e)*M*e = 1
        // This means, find the coefficients of the matrix M such that all the edges have unit length
        Vector sol;
        array_1d<double, 6> aux_vect;
        BoundedMatrix<double, 6, 6> aux_mat;
        unsigned int row = 0;
        for (unsigned int i = 0; i < 3; ++i) {
            const auto& i_coord = rGeometry[i].Coordinates();
            for (unsigned int j = i + 1; j < 4; ++j) {
                const auto& j_coord = rGeometry[j].Coordinates();
                aux_mat(row, 0) = std::pow(i_coord[0]-j_coord[0], 2);
                aux_mat(row, 1) = 2.0*(i_coord[0]-j_coord[0])*(i_coord[1]-j_coord[1]);
                aux_mat(row, 2) = 2.0*(i_coord[0]-j_coord[0])*(i_coord[2]-j_coord[2]);
                aux_mat(row, 3) = std::pow(i_coord[1]-j_coord[1], 2);
                aux_mat(row, 4) = 2.0*(i_coord[1]-j_coord[1])*(i_coord[2]-j_coord[2]);
                aux_mat(row, 5) = std::pow(i_coord[2]-j_coord[2], 2);
                aux_vect(row) = 1.0;
                row++;
            }
        }
        MathUtils<double>::Solve(aux_mat, sol, aux_vect);

        // Set the metric tensor
        Matrix metric(3,3);
        metric(0,0) = sol[0]; metric(0,1) = sol[1]; metric(0,2) = sol[2];
        metric(1,0) = sol[1]; metric(1,1) = sol[3]; metric(1,2) = sol[4];
        metric(2,0) = sol[2]; metric(2,1) = sol[4]; metric(2,2) = sol[5];

        // Calculate the eigenvalues of the metric tensor to obtain the ellipsis of inertia axes lengths
        BoundedMatrix<double,3,3> eigenvects, eigenvals;
        MathUtils<double>::GaussSeidelEigenSystem(metric, eigenvects, eigenvals);
        const double h_1 = std::sqrt(1.0 / eigenvals(0,0));
        const double h_2 = std::sqrt(1.0 / eigenvals(1,1));
        const double h_3 = std::sqrt(1.0 / eigenvals(2,2));

        // Calculate the reference element size as the average of the ellipsis of intertia axes lengths
        const double h_ref = (h_1 + h_2 + h_3) / 3.0;
        
        // Make the metric dimensionless
        metric *= std::pow(h_ref,2);

        // Calculate metric infimum norm
        const double metric_inf = std::min(eigenvals(0,0), std::min(eigenvals(1,1), eigenvals(2,2)));

        return std::make_tuple(h_ref, metric_inf, metric);
    }

    void ApplySlipCondition()
    {
        // Calculate the model part data
        auto &r_model_part = BaseType::GetModelPart();
        const int n_nodes = r_model_part.NumberOfNodes();

        // Calculate and substract the normal contribution
#pragma omp parallel for
        for (unsigned int i_node = 0; i_node < n_nodes; ++i_node) {
            auto it_node = r_model_part.NodesBegin() + i_node;
            if (it_node->Is(SLIP)) {
                auto unit_normal = it_node->FastGetSolutionStepValue(NORMAL);
                unit_normal /= norm_2(unit_normal);
                auto& r_mom = it_node->FastGetSolutionStepValue(MOMENTUM);
                const double r_mom_n = inner_prod(r_mom, unit_normal);
                r_mom -= r_mom_n * unit_normal;
            }
        }
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
}; /* Class CompressibleNavierStokesExplicitSolvingStrategyRungeKutta4 */

///@}

///@name Type Definitions
///@{

///@}

} /* namespace Kratos.*/

#endif /* KRATOS_COMPRESSIBLE_NAVIER_STOKES_EXPLICIT_SOLVING_STRATEGY_RUNGE_KUTTA_4  defined */
