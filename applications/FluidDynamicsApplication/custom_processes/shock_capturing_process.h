//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main author:     Ruben Zorrilla
//

#if !defined(KRATOS_SHOCK_CAPTURING_UTILITIES_H_INCLUDED)
#define  KRATOS_SHOCK_CAPTURING_UTILITIES_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "containers/model.h"
#include "includes/define.h"
#include "includes/model_part.h"
#include "processes/process.h"
#include "utilities/atomic_utilities.h"

// Application includes
#include "fluid_dynamics_application_variables.h"

namespace Kratos
{
///@addtogroup FluidDynamicsApplication
///@{

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

class KRATOS_API(FLUID_DYNAMICS_APPLICATION) ShockCapturingProcess : public Process
{
public:

    ///@name Type Definitions
    ///@{

    /// The base process type
    typedef Process BaseType;

    /// Type for the metric calculation function
    typedef std::function<std::tuple<double, double, Matrix>(Geometry<Node<3>> &rGeometry)> ElementMetricFunctionType;

    /// Type for the 2D (linear triangle) TLS geometry data
    typedef std::tuple<double, BoundedMatrix<double, 3, 2>, array_1d<double, 3>> ShockCapturingTLSType2D3N;

    /// Type for the 3D (linear tetrahedra) TLS geometry data
    typedef std::tuple<double, BoundedMatrix<double, 4, 3>, array_1d<double, 4>> ShockCapturingTLSType3D4N;

    /// Pointer definition of ShockCapturingProcess
    KRATOS_CLASS_POINTER_DEFINITION(ShockCapturingProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor with model
    ShockCapturingProcess(
        Model& rModel,
        Parameters& rParameters)
        : mrModelPart(rModel.GetModelPart(rParameters["model_part_name"].GetString()))
    {
        ValidateAndAssignParameters(rParameters);
    };

    /// Constructor with model part
    ShockCapturingProcess(
        ModelPart& rModelPart,
        Parameters& rParameters)
        : mrModelPart(rModelPart)
    {
        ValidateAndAssignParameters(rParameters);
    };

    /// Destructor.
    ~ShockCapturingProcess() {};

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    void Execute() override;

    void ExecuteInitialize() override;

    void ExecuteFinalizeSolutionStep() override;

    int Check() override;

    const Parameters GetDefaultParameters() const override;

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "ShockCapturingProcess";
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream &rOStream) const override { rOStream << "ShockCapturingProcess"; }

    /// Print object's data.
    void PrintData(std::ostream &rOStream) const override {}

    ///@}
    ///@name Friends
    ///@{


    ///@}
private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    ModelPart& mrModelPart;

    bool mUpdateNodalArea;
    bool mShockSensor;
    bool mShearSensor;
    bool mThermalSensor;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    void ValidateAndAssignParameters(Parameters &rParameters);

    /**
     * @brief Set the element metric function object
     * This method checks the first element type and sets the function to calculate
     * the metric tensor accordingly. Note that, for the sake of performance, it is
     * assumed that the element geometry type in the mesh is unique.
     * Such metric tensor calculation function returns a tuple that contains
     * 0 : The reference element size used to build the metric tensor
     * 1 : The infimum norm of the metric tensor
     * 2 : The metric tensor relative to the reference element size stored in 0
     * @return ElementMetricFunctionType Function to calculate the metric tensor
     */
    ElementMetricFunctionType SetElementMetricFunction();

    ElementMetricFunctionType SetShockCapturingTLSContainer();

    /**
     * @brief Physics-based shock capturing
     * This function calculates the artificial magnitudes using a physics-based shock capturing method.
     * References https://arc.aiaa.org/doi/abs/10.2514/6.2018-0062
     */
    void CalculatePhysicsBasedShockCapturing();

    template<class TTLSContainerType>
    void CalculatePhysicsBasedShockCapturingElementContribution(
        Element &rElement,
        TTLSContainerType &rShockCapturingTLS)
    {
        auto &r_geom = rElement.GetGeometry();
        const unsigned int n_nodes = r_geom.PointsNumber();

        // Get TLS values and calculate geometry data
        double& r_vol = std::get<0>(rShockCapturingTLS);
        auto& r_DN_DX = std::get<1>(rShockCapturingTLS);
        auto& r_N = std::get<2>(rShockCapturingTLS);
        GeometryUtils::CalculateGeometryData(r_geom, r_DN_DX, r_N, r_vol);

        // Get fluid physical properties
        const auto& r_prop = rElement.GetProperties();
        const double c_v = r_prop.GetValue(SPECIFIC_HEAT);
        const double gamma = r_prop.GetValue(HEAT_CAPACITY_RATIO);

        // Calculate elemental magnitudes
        const double k = 1.0; // Polynomial order of the numerical simulation
        double c_ref;         // Elemental speed of sound
        // TODO: CALLING THE CALCULATES IS NOT THE MOST EFFICIENT WAY... THINK ABOUT THIS...


        // rElement.Calculate(DENSITY_GRADIENT, grad_rho, r_process_info); // Midpoint density gradient --> Shock sensor (h_beta)
        // rElement.Calculate(VELOCITY_DIVERGENCE, div_v, r_process_info); // Migpoint velocity divergence --> Shock sensor
        // rElement.Calculate(VELOCITY_ROTATIONAL, rot_v, r_process_info); // Midpoint velocity rotational --> Shock sensor


        // rElement.Calculate(SOUND_VELOCITY, c_ref, r_process_info);       // Midpoint sound velocity --> Shear sensor
        // rElement.Calculate(VELOCITY_GRADIENT, grad_vel, r_process_info); // Midpoint velocity gradient --> Shear sensor


        // rElement.Calculate(TEMPERATURE_GRADIENT, grad_temp, r_process_info); // Temperature gradient --> Thermal sensor



        // Calculate midpoint values
        double midpoint_rho = 0.0;
        double midpoint_tot_ener = 0.0;
        array_1d<double, 3> midpoint_v = ZeroVector(3);
        for (unsigned int i_node = 0; i_node < n_nodes; ++i_node)
        {
            const double &r_rho = r_geom[i_node].FastGetSolutionStepValue(DENSITY);
            const double &r_tot_ener = r_geom[i_node].FastGetSolutionStepValue(TOTAL_ENERGY);
            midpoint_rho += r_rho;
            midpoint_tot_ener += r_tot_ener;
            midpoint_v += r_geom[i_node].FastGetSolutionStepValue(MOMENTUM) / r_rho;
        }
        midpoint_rho /= static_cast<double>(n_nodes);
        midpoint_v /= static_cast<double>(n_nodes);
        const double v_norm_pow = midpoint_v[0] * midpoint_v[0] + midpoint_v[1] * midpoint_v[1] + midpoint_v[2] * midpoint_v[2];
        const double stagnation_temp = midpoint_tot_ener / midpoint_rho / c_v;
        const double c_star = std::sqrt(gamma * (gamma - 1.0) * c_v * stagnation_temp * (2.0 / (gamma + 1.0))); // Critical speed of sound

        // Inverse metric tensor calculation
        // const double h_ref = avg_h_function(r_geom); // Reference element size used in the metric tensor
        // const auto metric_tensor = elem_metric_function(r_geom, h_ref); // Metric tensor relative to the reference element size
        const auto metric_data = elem_metric_function(r_geom);
        const double h_ref = std::get<0>(metric_data);             // Reference element size used in the metric tensor
        const double metric_tensor_inf = std::get<1>(metric_data); // Metric tensor infimum norm (smallest eigenvalue)
        const auto metric_tensor = std::get<2>(metric_data);       // Metric tensor relative to the reference element size
        double aux_det;
        Matrix inv_metric_tensor;
        MathUtils<double>::InvertMatrix(metric_tensor, inv_metric_tensor, aux_det);

        // Characteristic element sizes
        array_1d<double, 3> inv_metric_grad_rho = ZeroVector(3);
        array_1d<double, 3> inv_metric_grad_temp = ZeroVector(3);
        for (unsigned int i = 0; i < inv_metric_tensor.size1(); ++i)
        {
            for (unsigned int j = 0; j < inv_metric_tensor.size2(); ++j)
            {
                inv_metric_grad_rho(i) += inv_metric_tensor(i, j) * grad_rho(j);
                inv_metric_grad_temp(i) += inv_metric_tensor(i, j) * grad_temp(j);
            }
        }
        const double h_beta = h_ref * norm_2(grad_rho) / std::sqrt(inner_prod(grad_rho, inv_metric_grad_rho) + eps);     // Characteristic element size along the direction of the density gradient
        const double h_kappa = h_ref * norm_2(grad_temp) / std::sqrt(inner_prod(grad_temp, inv_metric_grad_temp) + eps); // Characteristic element size along the direction of the temperature gradient
        const double h_mu = h_ref * metric_tensor_inf;

        // Dilatation sensor (activates in shock waves)
        const double s_omega = -h_beta * div_v / k / c_star;

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
        rElement.GetValue(SHOCK_SENSOR) = s_beta_hat;

        // Thermal sensor (detect thermal gradients that are larger than possible with the grid resolution)
        Matrix mid_pt_jacobian;
        r_geom.Jacobian(mid_pt_jacobian, 0, GeometryData::GI_GAUSS_1);
        array_1d<double, 3> local_grad_temp = ZeroVector(3);
        for (unsigned int i = 0; i < mid_pt_jacobian.size1(); ++i)
        {
            for (unsigned int j = 0; j < mid_pt_jacobian.size2(); ++j)
            {
                local_grad_temp(i) += mid_pt_jacobian(j, i) * grad_temp(j);
            }
        }

        const double s_kappa_0 = 1.0;
        const double s_kappa_max = 2.0;
        const double s_kappa = h_ref * norm_2(local_grad_temp) / k / stagnation_temp;
        const double s_kappa_hat = SmoothedLimitingFunction(s_kappa, s_kappa_0, s_kappa_max);
        rElement.GetValue(THERMAL_SENSOR) = s_kappa_hat;

        // Shear sensor (detect velocity gradients that are larger than possible with the grid resolution)
        const unsigned int dim = r_geom.WorkingSpaceDimension();
        Matrix shear_grad_vel(dim, dim);
        for (unsigned int d1 = 0; d1 < dim; ++d1)
        {
            for (unsigned int d2 = 0; d2 < dim; ++d2)
            {
                shear_grad_vel(d1, d2) = d1 == d2 ? 0.0 : grad_vel(d1, d2);
            }
        }
        const Matrix local_shear_grad_vel = prod(shear_grad_vel, trans(mid_pt_jacobian));
        Matrix eigen_vect_mat, eigen_val_mat;
        MathUtils<double>::GaussSeidelEigenSystem(local_shear_grad_vel, eigen_vect_mat, eigen_val_mat);
        double shear_spect_norm = 0.0;
        for (unsigned int d = 0; d < eigen_val_mat.size1(); ++d)
        {
            if (eigen_val_mat(d, d) > shear_spect_norm)
            {
                shear_spect_norm = eigen_val_mat(d, d);
            }
        }
        const double isentropic_max_vel = std::sqrt(v_norm_pow + (2.0 / (gamma - 1.0)) * std::pow(c_ref, 2));

        const double s_mu_0 = 1.0;
        const double s_mu_max = 2.0;
        const double s_mu = h_ref * shear_spect_norm / isentropic_max_vel / k;
        // const double s_mu_hat = LimitingFunction(s_mu, s_mu_0, s_mu_max);
        const double s_mu_hat = SmoothedLimitingFunction(s_mu, s_mu_0, s_mu_max);
        rElement.GetValue(SHEAR_SENSOR) = s_mu_hat;

        // Calculate artificial magnitudes
        const double ref_mom_norm = midpoint_rho * std::sqrt(v_norm_pow + std::pow(c_star, 2));

        // Calculate elemental artificial bulk viscosity
        const double k_beta = 1.5;
        const double elem_b_star = (k_beta * h_beta / k) * ref_mom_norm * s_beta_hat;
        rElement.GetValue(ARTIFICIAL_BULK_VISCOSITY) = elem_b_star;

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
        rElement.GetValue(ARTIFICIAL_CONDUCTIVITY) = elem_k1_star + elem_k2_star;

        // Calculate elemental artificial dynamic viscosity
        const double k_mu = 1.0;
        const double elem_mu_star = (k_mu * h_mu / k) * ref_mom_norm * s_mu_hat;
        rElement.GetValue(ARTIFICIAL_DYNAMIC_VISCOSITY) = elem_mu_star;

        // Project the shock capturing magnitudes to the nodes
        const double geom_domain_size = r_geom.DomainSize();
        const double aux_weight = geom_domain_size / static_cast<double>(n_nodes);
        for (unsigned int i_node = 0; i_node < n_nodes; ++i_node)
        {
            auto &r_node = r_geom[i_node];
            AtomicAdd(r_node.GetValue(ARTIFICIAL_BULK_VISCOSITY), aux_weight * elem_b_star);
            AtomicAdd(r_node.GetValue(ARTIFICIAL_CONDUCTIVITY), aux_weight * (elem_k1_star + elem_k2_star));
            AtomicAdd(r_node.GetValue(ARTIFICIAL_DYNAMIC_VISCOSITY), aux_weight * elem_mu_star);
        }
    }

    double LimitingFunction(
        const double s,
        const double s_0,
        const double s_max,
        const double s_min = 0.0);

    double SmoothedLimitingFunction(
        const double s,
        const double s_0,
        const double s_max);

    // Smooth approximation of the max(s,0) function
    double SmoothedMaxFunction(const double s);

    // Smooth approximation of the min(s,0) function
    double SmoothedMinFunction(const double s);

    // TODO: SHALL WE MOVE THIS TO THE FLUID_ELEMENT_UTILITIES?
    // https://es.wikipedia.org/wiki/Circunelipse_de_Steiner
    std::tuple<double, double, Matrix> CalculateTriangleMetricTensor(const Geometry<Node<3>> &rGeometry);

    // TODO: SHALL WE MOVE THIS TO THE FLUID_ELEMENT_UTILITIES?
    // https://es.wikipedia.org/wiki/Circunelipse_de_Steiner --> 3D extension
    std::tuple<double, double, Matrix> CalculateTetrahedraMetricTensor(const Geometry<Node<3>> &rGeometry);

    // TODO: REMOVE AFTER SUNETH'S PR
    template <class TDataType>
    void UpdateValue(TDataType &rOutput, const TDataType &rInput);

    // TODO: REMOVE AFTER SUNETH'S PR
    template <class... TRefVariableValuePairArgs>
    void EvaluateInPoint(
        const Geometry<Node<3>> &rGeometry,
        const Vector &rShapeFunction,
        const int Step,
        const TRefVariableValuePairArgs &... rValueVariablePairs)
    {
        KRATOS_TRY

        const int number_of_nodes = rGeometry.PointsNumber();

        const auto &r_node = rGeometry[0];
        const double shape_function_value = rShapeFunction[0];

        int dummy[sizeof...(TRefVariableValuePairArgs)] = {(
            std::get<0>(rValueVariablePairs) =
                r_node.FastGetSolutionStepValue(std::get<1>(rValueVariablePairs), Step) * shape_function_value,
            0)...};

        // this can be removed with fold expressions in c++17
        *dummy = 0;

        for (int c = 1; c < number_of_nodes; ++c)
        {
            const auto &r_node = rGeometry[c];
            const double shape_function_value = rShapeFunction[c];

            int dummy[sizeof...(TRefVariableValuePairArgs)] = {(
                UpdateValue<typename std::remove_reference<typename std::tuple_element<0, TRefVariableValuePairArgs>::type>::type>(
                    std::get<0>(rValueVariablePairs),
                    r_node.FastGetSolutionStepValue(std::get<1>(rValueVariablePairs), Step) * shape_function_value),
                0)...};

            // this can be removed with fold expressions in c++17
            *dummy = 0;
        }

        KRATOS_CATCH("");
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

    /// Assignment operator.
    ShockCapturingProcess& operator=(ShockCapturingProcess const& rOther);

    /// Copy constructor.
    ShockCapturingProcess(ShockCapturingProcess const& rOther);

    ///@}

}; // Class ShockCapturingProcess

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator << (
    std::ostream& rOStream,
    const ShockCapturingProcess& rThis);

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_SHOCK_CAPTURING_UTILITIES_H_INCLUDED  defined
