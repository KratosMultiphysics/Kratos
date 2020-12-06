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
#include "utilities/geometry_metric_calculator.h"

// Application includes
#include "custom_utilities/fluid_calculation_utilities.h"
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
    typedef std::tuple<
        double,
        BoundedMatrix<double,3,2>,
        array_1d<double,3>,
        BoundedMatrix<double,2,2>,
        array_1d<double,3>,
        array_1d<double,3>,
        array_1d<double,3>,
        array_1d<double,3>,
        array_1d<double,3>,
        BoundedMatrix<double,2,2>> ShockCapturingTLSType2D3N;

    /// Type for the 3D (linear tetrahedra) TLS geometry data
    typedef std::tuple<
        double,
        BoundedMatrix<double,4,3>,
        array_1d<double,4>,
        BoundedMatrix<double,3,3>,
        array_1d<double,3>,
        array_1d<double,3>,
        array_1d<double,3>,
        array_1d<double,3>,
        array_1d<double,3>,
        BoundedMatrix<double,3,3>> ShockCapturingTLSType3D4N;

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
    bool mThermallyCoupledFormulation;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    void ValidateAndAssignParameters(Parameters &rParameters);

    /**
     * @brief Physics-based shock capturing
     * This function calculates the artificial magnitudes using a physics-based shock capturing method.
     * References https://arc.aiaa.org/doi/abs/10.2514/6.2018-0062
     */
    void CalculatePhysicsBasedShockCapturing();

    /**
     * @brief Set the tetrahedra TLS container
     * This method sets the TLS container for a 2D3N triangle geometry
     * @return ShockCapturingTLSType2D3N The TLS container
     */
    ShockCapturingTLSType2D3N SetTLSContainer2D3N();

    /**
     * @brief Set the tetrahedra TLS container
     * This method sets the TLS container for a 3D4N tetrahedra geometry
     * @return ShockCapturingTLSType3D4N The TLS container
     */
    ShockCapturingTLSType3D4N SetTLSContainer3D4N();

    /**
     * @brief Calculate elemental shock capturing contribution
     * This method calculates the elemental physics-based shock capturing contribution
     * It is intended to be called from the lambda function used in the parallel region
     * @tparam TDim Problem dimension
     * @tparam TNumNodes Geometry number of nodes
     * @tparam TTLSContainerType The TLS container type
     * @param rElement Reference to the element of interest
     * @param rShockCapturingTLS Reference to the TLS container
     */
    template<std::size_t TDim, std::size_t TNumNodes, class TTLSContainerType>
    void KRATOS_API(FLUID_DYNAMICS_APPLICATION) CalculatePhysicsBasedShockCapturingElementContribution(
        Element &rElement,
        TTLSContainerType &rShockCapturingTLS)
    {
        auto &r_geom = rElement.GetGeometry();
        const unsigned int n_nodes = r_geom.PointsNumber();

        // Initialize artificial magnitudes values
        if (mShockSensor) {rElement.GetValue(ARTIFICIAL_BULK_VISCOSITY) = 0.0;}
        if (mShearSensor) {rElement.GetValue(ARTIFICIAL_DYNAMIC_VISCOSITY) = 0.0;}
        if (mThermalSensor || mShockSensor) {rElement.GetValue(ARTIFICIAL_CONDUCTIVITY) = 0.0;} // Note that conductivity is modified in both sensors

        // Get TLS values geometry values
        double& r_vol = std::get<0>(rShockCapturingTLS);
        auto& r_DN_DX = std::get<1>(rShockCapturingTLS);
        auto& r_N = std::get<2>(rShockCapturingTLS);
        auto& r_metric_tensor = std::get<3>(rShockCapturingTLS);

        // Calculate geometry data
        GeometryUtils::CalculateGeometryData(r_geom, r_DN_DX, r_N, r_vol);

        // Calculate the geometry metric
        double h_ref, metric_tensor_inf, metric_tensor_sup;
        GeometryMetricCalculator<TDim,TNumNodes>::CalculateMetricTensorDimensionless(
            r_geom,
            r_metric_tensor,
            h_ref,
            metric_tensor_inf,
            metric_tensor_sup);

        // Inverse metric tensor calculation
        double aux_det;
        Matrix inv_metric_tensor; //FIXME: We should use a bounded matrix in here
        MathUtils<double>::InvertMatrix(r_metric_tensor, inv_metric_tensor, aux_det);

        // Get fluid physical properties
        const auto& r_prop = rElement.GetProperties();
        const double c_v = r_prop.GetValue(SPECIFIC_HEAT);
        const double gamma = r_prop.GetValue(HEAT_CAPACITY_RATIO);

        // Set auxiliary constants
        const double k = 1.0; // Polynomial order of the numerical simulation
        const double eps = 1.0e-7; // Small constant to avoid division by 0

        // Calculate the midpoint values
        double midpoint_rho, midpoint_tot_ener;
        auto& r_midpoint_v = std::get<4>(rShockCapturingTLS);
        if (mThermallyCoupledFormulation) {
            // Get required midpoint values
            double midpoint_temp;
            FluidCalculationUtilities::EvaluateInPoint(r_geom, r_N, std::tie(r_midpoint_v, VELOCITY), std::tie(midpoint_rho, DENSITY), std::tie(midpoint_temp, TEMPERATURE));
            // If the formulation is thermally coupled, the total energy is the summation of the thermal and kinetic ones
            midpoint_tot_ener = midpoint_rho * (c_v * midpoint_temp + 0.5 * midpoint_rho * inner_prod(r_midpoint_v, r_midpoint_v));
        } else {
            // Get required midpoint values
            FluidCalculationUtilities::EvaluateInPoint(r_geom, r_N, std::tie(r_midpoint_v, VELOCITY), std::tie(midpoint_rho, DENSITY));
            // If the formulation is not energy coupled, the total energy equals the kinetic one
            midpoint_tot_ener = 0.5 * midpoint_rho * inner_prod(r_midpoint_v, r_midpoint_v);
        }

        // Calculate common values
        const double v_norm_pow = SquaredArrayNorm(r_midpoint_v);
        const double stagnation_temp = midpoint_tot_ener / midpoint_rho / c_v;
        const double c_star = std::sqrt(gamma * (gamma - 1.0) * c_v * stagnation_temp * (2.0 / (gamma + 1.0))); // Critical speed of sound
        const double ref_mom_norm = midpoint_rho * std::sqrt(v_norm_pow + std::pow(c_star, 2));

        // Shock sensor values
        if (mShockSensor) {
            // Calculate the required differential operators
            double div_v, mach;
            auto& r_grad_rho = std::get<5>(rShockCapturingTLS);
            auto& r_rot_v = std::get<6>(rShockCapturingTLS);
            CalculateShockSensorValues<TDim,TNumNodes>(r_geom, r_N, r_DN_DX, mach, div_v, r_grad_rho, r_rot_v);

            // Characteristic element size along the direction of the density gradient
            const double h_beta = h_ref * norm_2(r_grad_rho) / std::sqrt(CalculateProjectedInverseMetricElementSize(inv_metric_tensor, r_grad_rho) + eps);

            // Dilatation sensor (activates in shock waves)
            const double s_omega = -h_beta * div_v / k / c_star;

            // Vorticity sensor (vanishes in vorticity dominated regions)
            const double div_v_pow = std::pow(div_v, 2);
            const double rot_v_norm_pow = SquaredArrayNorm(r_rot_v);
            const double s_w = div_v_pow / (div_v_pow + rot_v_norm_pow + eps);

            // Calculate limited shock sensor
            const double s_beta_0 = 0.01;
            const double s_beta_max = 2.0 / std::sqrt(std::pow(gamma, 2) - 1.0);
            const double s_beta = s_omega * s_w;
            // const double s_beta_hat = LimitingFunction(s_beta, s_beta_0, s_beta_max);
            const double s_beta_hat = SmoothedLimitingFunction(s_beta, s_beta_0, s_beta_max);
            rElement.GetValue(SHOCK_SENSOR) = s_beta_hat;

            // Calculate elemental artificial bulk viscosity
            const double k_beta = 1.5;
            const double elem_b_star = (k_beta * h_beta / k) * ref_mom_norm * s_beta_hat;
            rElement.GetValue(ARTIFICIAL_BULK_VISCOSITY) = elem_b_star;

            // Calculate elemental artificial conductivity (dilatancy)
            const double Pr_beta_min = 0.9;
            const double alpha_pr_beta = 2.0;
            const double mach_threshold = 3.0;
            const double Pr_beta = Pr_beta_min * (1.0 + std::exp(-2.0 * alpha_pr_beta * (mach - mach_threshold)));
            const double elem_k1_star = (gamma * c_v / Pr_beta) * elem_b_star;
            rElement.GetValue(ARTIFICIAL_CONDUCTIVITY) += elem_k1_star;
        }

        if (mThermalSensor || mShearSensor) {
            // Calculate Jacobian matrix (non-required for the shock sensor)
            Matrix mid_pt_jacobian; //FIXME: We should use a bounded matrix in here
            r_geom.Jacobian(mid_pt_jacobian, 0, GeometryData::GI_GAUSS_1);

            // Thermal sensor values
            if (mThermalSensor) {
                // Calculate temperature gradients
                auto& r_grad_temp = std::get<7>(rShockCapturingTLS);
                auto& r_grad_temp_local = std::get<8>(rShockCapturingTLS);
                CalculateTemperatureGradients(r_geom, r_DN_DX, mid_pt_jacobian, r_grad_temp, r_grad_temp_local);

                // Characteristic element size along the direction of the temperature gradient
                const double h_kappa = h_ref * norm_2(r_grad_temp) / std::sqrt(CalculateProjectedInverseMetricElementSize(inv_metric_tensor, r_grad_temp) + eps);

                // Thermal sensor (detect thermal gradients that are larger than possible with the grid resolution)
                const double s_kappa_0 = 1.0;
                const double s_kappa_max = 2.0;
                const double s_kappa = h_ref * norm_2(r_grad_temp_local) / k / stagnation_temp;
                const double s_kappa_hat = SmoothedLimitingFunction(s_kappa, s_kappa_0, s_kappa_max);
                rElement.GetValue(THERMAL_SENSOR) = s_kappa_hat;

                // Calculate elemental artificial conductivity (thermal sensor)
                const double k_kappa = 1.0;
                const double elem_k2_star = (gamma * c_v) * (k_kappa * h_kappa / k) * ref_mom_norm * s_kappa_hat;
                rElement.GetValue(ARTIFICIAL_CONDUCTIVITY) += elem_k2_star;
            }

            // Shear sensor values
            if (mShearSensor) {
                // Calculate shear sensor values
                double r_c;
                auto& r_local_shear_grad_v = std::get<9>(rShockCapturingTLS);
                CalculateShearSensorValues(r_geom, r_N, r_DN_DX, mid_pt_jacobian, r_local_shear_grad_v, r_c);
                BoundedMatrix<double,TDim,TDim> eigen_vect_mat, eigen_val_mat;
                MathUtils<double>::GaussSeidelEigenSystem(r_local_shear_grad_v, eigen_vect_mat, eigen_val_mat);
                double shear_spect_norm = 0.0;
                for (unsigned int d = 0; d < eigen_val_mat.size1(); ++d) {
                    if (eigen_val_mat(d, d) > shear_spect_norm) {
                        shear_spect_norm = eigen_val_mat(d, d);
                    }
                }
                shear_spect_norm = std::sqrt(shear_spect_norm);

                // Characteristic element size for the shear sensor
                const double h_mu = h_ref * metric_tensor_inf;

                // Shear sensor (detect velocity gradients that are larger than possible with the grid resolution)
                const double isentropic_max_vel = std::sqrt(v_norm_pow + (2.0 / (gamma - 1.0)) * std::pow(r_c, 2));
                const double s_mu_0 = 1.0;
                const double s_mu_max = 2.0;
                const double s_mu = h_ref * shear_spect_norm / isentropic_max_vel / k;
                // const double s_mu_hat = LimitingFunction(s_mu, s_mu_0, s_mu_max);
                const double s_mu_hat = SmoothedLimitingFunction(s_mu, s_mu_0, s_mu_max);
                rElement.GetValue(SHEAR_SENSOR) = s_mu_hat;

                // Calculate elemental artificial dynamic viscosity
                const double k_mu = 1.0;
                const double elem_mu_star = (k_mu * h_mu / k) * ref_mom_norm * s_mu_hat;
                rElement.GetValue(ARTIFICIAL_DYNAMIC_VISCOSITY) = elem_mu_star;
            }
        }

        // Project the shock capturing magnitudes to the nodes
        const double aux_weight = r_vol / static_cast<double>(n_nodes);
        for (unsigned int i_node = 0; i_node < n_nodes; ++i_node) {
            auto &r_node = r_geom[i_node];
            if (mShockSensor) {AtomicAdd(r_node.GetValue(ARTIFICIAL_BULK_VISCOSITY), aux_weight * rElement.GetValue(ARTIFICIAL_BULK_VISCOSITY));}
            if (mShearSensor) {AtomicAdd(r_node.GetValue(ARTIFICIAL_DYNAMIC_VISCOSITY), aux_weight * rElement.GetValue(ARTIFICIAL_DYNAMIC_VISCOSITY));}
            if (mShockSensor || mThermalSensor) {AtomicAdd(r_node.GetValue(ARTIFICIAL_CONDUCTIVITY), aux_weight * rElement.GetValue(ARTIFICIAL_CONDUCTIVITY));}
        }
    }

    /**
     * @brief Calculate the squared norm of an array_1d<double,3>
     * Auxiliary inline function to calculate the squared norm of an array_1d<double,3>
     * @param rArray The array to calculate the squared norm
     * @return double The squared norm of the given array
     */
    inline double SquaredArrayNorm(const array_1d<double,3>& rArray)
    {
        return rArray[0]*rArray[0] + rArray[1]*rArray[1] + rArray[2]*rArray[2];
    }

    /**
     * @brief Limiting function
     * Function to limit a given value between two bounds
     * @param s Value to limit
     * @param s_0 Initial value
     * @param s_max Maximum value
     * @param s_min Minimum values
     * @return double Provided value or the corresponding minimum or maximum one
     */
    double LimitingFunction(
        const double s,
        const double s_0,
        const double s_max,
        const double s_min = 0.0);

    /**
     * @brief Smoothed limiting function
     * Function to smoothly limit a given value between two bounds
     * @param s Value to limit
     * @param s_0 Initial value
     * @param s_max Maximum value
     * @return double Provided value or the corresponding minimum or maximum one
     */
    double SmoothedLimitingFunction(
        const double s,
        const double s_0,
        const double s_max);

    /**
     * @brief Smoothed max function
     * Smooth approximation of the max function
     * @param s Value to check
     * @return double Provided value or the maximum one
     */
    double SmoothedMaxFunction(const double s);

    /**
     * @brief Smoothed min function
     * Smooth approximation of the min(s,0) function
     * @param s Value to check
     * @return double Provided value or the minimum one
     */
    double SmoothedMinFunction(const double s);

    /**
     * @brief Calculate projected metric-based element size
     * This method calculates the element size as trans(grad(a))*inv(M)*grad(a)
     * that is the element size along the projection of the provided gradient
     * @param rInverseMetricTensor Inverse of the geometry metric tensor
     * @param rScalarGradient Scalar magnitude gradient to be projected
     * @return double Projected element size
     */
    double CalculateProjectedInverseMetricElementSize(
        const Matrix& rInverseMetricTensor,
        const array_1d<double,3>& rScalarGradient);

    /**
     * @brief Calculate the shock sensor values
     * This method calculates the magnitudes required for the shock sensor
     * @tparam TDim Problem dimension
     * @tparam TNumNodes Geometry number of nodes
     * @param rGeometry Element geometry
     * @param rN Shape function values
     * @param rDN_DX Shape function gradients
     * @param rMachNumber Computed Mach number
     * @param rVelocityDivergence Computed velocity divergence
     * @param rDensityGradient Computed density gradient
     * @param rVelocityRotational Computed velocity rotational
     */
    template<std::size_t TDim, std::size_t TNumNodes>
    void KRATOS_API(FLUID_DYNAMICS_APPLICATION) CalculateShockSensorValues(
        const Geometry<Node<3>>& rGeometry,
        const array_1d<double,TNumNodes>& rN,
        const BoundedMatrix<double,TNumNodes,TDim>& rDN_DX,
        double& rMachNumber,
        double& rVelocityDivergence,
        array_1d<double,3>& rDensityGradient,
        array_1d<double,3>& rVelocityRotational);

    /**
     * @brief Calculate the temperature gradients
     * This method calculates the temperature gradients required by the thermal sensor
     * @tparam TDim Problem dimension
     * @tparam TNumNodes Geometry number of nodes
     * @param rGeometry Element geometry
     * @param rDN_DX Shape function gradients
     * @param rJacobianMatrix Geometry Jacobian matrix
     * @param rTemperatureGradient Computed temperature gradient
     * @param rTemperatureLocalGradient Computed temperature local gradient
     */
    template<std::size_t TDim, std::size_t TNumNodes>
    void KRATOS_API(FLUID_DYNAMICS_APPLICATION) CalculateTemperatureGradients(
        const Geometry<Node<3>>& rGeometry,
        const BoundedMatrix<double,TNumNodes,TDim>& rDN_DX,
        const Matrix& rJacobianMatrix,
        array_1d<double,3>& rTemperatureGradient,
        array_1d<double,3>& rTemperatureLocalGradient);

    /**
     * @brief Calculate the shear sensor values
     * This method calculates the values required by the shear sensor
     * @tparam TDim Problem dimension
     * @tparam TNumNodes Geometry number of nodes
     * @param rGeometry Element geometry
     * @param rN Shape function values
     * @param rDN_DX Shape function gradients
     * @param rJacobianMatrix Geometry Jacobian matrix
     * @param rLocalVelocityShearGradient Computed velocity shear local gradient
     * @param rSoundVelocity Computed sound velocity
     */
    template<std::size_t TDim, std::size_t TNumNodes>
    void KRATOS_API(FLUID_DYNAMICS_APPLICATION) CalculateShearSensorValues(
        const Geometry<Node<3>>& rGeometry,
        const array_1d<double,TNumNodes>& rN,
        const BoundedMatrix<double,TNumNodes,TDim>& rDN_DX,
        const Matrix& rJacobianMatrix,
        BoundedMatrix<double,TDim,TDim>& rLocalVelocityShearGradient,
        double& rSoundVelocity);

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
