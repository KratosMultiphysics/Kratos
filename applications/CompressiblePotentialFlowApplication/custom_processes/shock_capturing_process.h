//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main author:     
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
///@addtogroup CompressiblePotentialFlowApplication
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

class KRATOS_API(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION) ShockCapturingProcess : public Process
{
public:

    ///@name Type Definitions
    ///@{

    /// The base process type
    typedef Process BaseType;

    /// Type for the metric calculation function
    typedef std::function<std::tuple<double, double, Matrix>(Geometry<Node> &rGeometry)> ElementMetricFunctionType;

    /// Type for the 2D (linear triangle) TLS geometry data
    struct ShockCapturingTLSContainer2D3N
    {
        double Vol;
        array_1d<double,3> N;
        BoundedMatrix<double,3,2> DN_DX;
        BoundedMatrix<double,2,2> MetricTensor;
        BoundedMatrix<double,2,2> InverseMetricTensor;
        array_1d<double,3> MidpointVelocity;
        array_1d<double,3> DensityGradient;
        array_1d<double,3> VelocityRotational;
    };

    /// Type for the 3D (linear tetrahedra) TLS geometry data
    struct ShockCapturingTLSContainer3D4N
    {
        double Vol;
        array_1d<double,4> N;
        BoundedMatrix<double,4,3> DN_DX;
        BoundedMatrix<double,3,3> MetricTensor;
        BoundedMatrix<double,3,3> InverseMetricTensor;
        array_1d<double,3> MidpointVelocity;
        array_1d<double,3> DensityGradient;
        array_1d<double,3> VelocityRotational;
    };

    /// Pointer definition of ShockCapturingProcess
    KRATOS_CLASS_POINTER_DEFINITION(ShockCapturingProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor with model
    ShockCapturingProcess(
        Model& rModel,
        Parameters rParameters)
        : ShockCapturingProcess(rModel.GetModelPart(rParameters["model_part_name"].GetString()), rParameters) {};

    /// Constructor with model part
    ShockCapturingProcess(
        ModelPart& rModelPart,
        Parameters rParameters)
        : Process()
        , mrModelPart(rModelPart)
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

    void ExecuteInitializeSolutionStep() override;

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
    void CalculatePhysicsBasedShockCapturingElementContribution(
        Element &rElement,
        TTLSContainerType &rShockCapturingTLS)
    {
        auto &r_geom = rElement.GetGeometry();

        // Get TLS values geometry values
        double& r_vol = rShockCapturingTLS.Vol;
        auto& r_DN_DX = rShockCapturingTLS.DN_DX;
        auto& r_N = rShockCapturingTLS.N;
        auto& r_metric_tensor = rShockCapturingTLS.MetricTensor;
        auto& r_inv_metric_tensor = rShockCapturingTLS.InverseMetricTensor;

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
        MathUtils<double>::InvertMatrix(r_metric_tensor, r_inv_metric_tensor, aux_det);

        // Get fluid physical properties
        const double c_v = 722.14;
        const double gamma = 1.4;

        // Set auxiliary constants
        const double k = 1.0; // Polynomial order of the numerical simulation
        const double eps = std::numeric_limits<double>::epsilon(); // Small constant to avoid division by 0

        // Calculate the midpoint values
        double midpoint_cp, midpoint_rho, formulation_temp;
        auto& r_midpoint_v = rShockCapturingTLS.MidpointVelocity;
        auto& r_grad_rho = rShockCapturingTLS.DensityGradient;
        auto& r_rot_v = rShockCapturingTLS.VelocityRotational;

        // Calculate the required differential operators
        double div_v;

        CalculateShockSensorValues<TDim,TNumNodes>(r_geom, r_N, r_DN_DX, r_midpoint_v, midpoint_cp, midpoint_rho, div_v, r_grad_rho, r_rot_v);

        formulation_temp = (0.5 * midpoint_rho * inner_prod(r_midpoint_v,r_midpoint_v) * midpoint_cp) / (midpoint_rho * (gamma - 1.0) * c_v);

        // Calculate common values
        const double stagnation_temp = formulation_temp + 0.5 * inner_prod(r_midpoint_v,r_midpoint_v)/(gamma * c_v);
        const double critical_temp = stagnation_temp * (2.0 / (gamma + 1.0));
        const double c_star = std::sqrt(gamma * (gamma - 1.0) * c_v * critical_temp); // Critical speed of sound

        // Characteristic element size along the direction of the density gradient
        const double h_beta = h_ref * norm_2(r_grad_rho) / std::sqrt(CalculateProjectedInverseMetricElementSize(r_inv_metric_tensor, r_grad_rho) + eps);

        // Dilatation sensor (activates in shock waves)
        // Critical speed of sound might be zero if zero initial condition is prescribed (i.e. thermally uncoupled framework)
        const double s_omega = c_star > 0.0 ? -h_beta * div_v / (k * c_star) : 0.0;

        // Vorticity sensor (vanishes in vorticity dominated regions)
        const double div_v_pow = std::pow(div_v, 2);
        const double rot_v_norm_pow = SquaredArrayNorm(r_rot_v);
        const double s_w = div_v_pow / (div_v_pow + rot_v_norm_pow + eps);

        // Calculate limited shock sensor
        const double s_beta_0 = 0.01;
        const double s_beta_max = 2.0 / std::sqrt(std::pow(gamma, 2) - 1.0);
        const double s_beta = s_omega * s_w;

        // KRATOS_INFO("ShockCapturing")  << "::::::::::::::::::::::::::::::::::::::::::::::::" << std::endl;  
        // KRATOS_INFO("ShockCapturing")  << "s_omega = " << s_omega << std::endl; 
        // KRATOS_INFO("ShockCapturing")  << "::::::::::::::::::::::::::::::::::::::::::::::::" << std::endl;       
        
        const double s_beta_hat = SmoothedLimitingFunction(s_beta, s_beta_0, s_beta_max);
        rElement.GetValue(SHOCK_SENSOR) = s_beta_hat;
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
     * @param rVelocity Computed velocity
     * @param rPressureCoefficient Computed pressure coefficient
     * @param rDensity Computed density
     * @param rVelocityDivergence Computed velocity divergence
     * @param rDensityGradient Computed density gradient
     * @param rVelocityRotational Computed velocity rotational
     */
    template<std::size_t TDim, std::size_t TNumNodes>
    void KRATOS_API(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION) CalculateShockSensorValues(
        const Geometry<Node>& rGeometry,
        const array_1d<double,TNumNodes>& rN,
        const BoundedMatrix<double,TNumNodes,TDim>& rDN_DX,
        array_1d<double,3>& rVelocity,
        double& rPressureCoefficient,
        double& rDensity,
        double& rVelocityDivergence,
        array_1d<double,3>& rDensityGradient,
        array_1d<double,3>& rVelocityRotational);

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
