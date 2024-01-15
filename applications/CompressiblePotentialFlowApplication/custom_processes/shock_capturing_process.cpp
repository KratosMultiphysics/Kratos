//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    
//
//

// System includes

// External includes

// Project includes
#include "geometries/geometry.h"
#include "geometries/geometry_data.h"
#include "utilities/geometry_utilities.h"
#include "utilities/parallel_utilities.h"

// Application includes
#include "shock_capturing_process.h"
#include "fluid_dynamics_application_variables.h"

namespace Kratos
{
    /* Public functions *******************************************************/

    void ShockCapturingProcess::Execute()
    {
        ExecuteInitialize();
        ExecuteInitializeSolutionStep();
    }

    void ShockCapturingProcess::ExecuteInitialize()
    {
        // Initialize elemental values
        block_for_each(mrModelPart.Elements(), [&](Element& rElement){
            rElement.SetValue(SHOCK_SENSOR, 0.0);
        });
    }

    void ShockCapturingProcess::ExecuteInitializeSolutionStep()
    {
        CalculatePhysicsBasedShockCapturing();
    }

    const Parameters ShockCapturingProcess::GetDefaultParameters() const
    {
        const Parameters default_parameters = Parameters(R"({
            "model_part_name" : ""
        })");
        return default_parameters;
    }

    /* Private functions *****************************************************/

    void ShockCapturingProcess::ValidateAndAssignParameters(Parameters& rParameters)
    {
        // Validate and assign defaults
        rParameters.ValidateAndAssignDefaults(GetDefaultParameters());
    }

    /**
     * @brief Physics-based shock capturing
     * This function calculates the artificial magnitudes using a physics-based shock capturing method.
     * References https://arc.aiaa.org/doi/abs/10.2514/6.2018-0062
     */
    void ShockCapturingProcess::CalculatePhysicsBasedShockCapturing()
    {
        // Calculate the elemental contributions of the shock capturing
        const auto geometry_type = (mrModelPart.ElementsBegin()->GetGeometry()).GetGeometryType();
        if (geometry_type == GeometryData::KratosGeometryType::Kratos_Triangle2D3) {
            // Set auxiliary TLS container and elemental function
            ShockCapturingTLSContainer2D3N tls_container_2D3N;
            auto aux_function_2D3N = [&, this] (Element &rElement, ShockCapturingTLSContainer2D3N &rShockCapturingTLS) {this->CalculatePhysicsBasedShockCapturingElementContribution<2,3>(rElement, rShockCapturingTLS);};
            // Perform the elemental loop
            block_for_each(mrModelPart.Elements(), tls_container_2D3N, aux_function_2D3N);
        } else if (geometry_type == GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4) {
            // Set auxiliary TLS container and elemental function
            ShockCapturingTLSContainer3D4N tls_container_3D4N;
            auto aux_function_3D4N = [&, this] (Element &rElement, ShockCapturingTLSContainer3D4N &rShockCapturingTLS) {this->CalculatePhysicsBasedShockCapturingElementContribution<3,4>(rElement, rShockCapturingTLS);};
            // Perform the elemental loop
            block_for_each(mrModelPart.Elements(), tls_container_3D4N, aux_function_3D4N);
        } else {
            KRATOS_ERROR << "Asking for a non-supported geometry. Shock capturing only supports \'Triangle2D3\' and \'Tetrahedra3D4\' geometries.";
        }
    }

    double ShockCapturingProcess::SmoothedLimitingFunction(
        const double s,
        const double s_0,
        const double s_max)
    {
        const double aux_1 = SmoothedMaxFunction(s - s_0);
        const double aux_2 = SmoothedMinFunction(aux_1 - s_max);
        return aux_2 + s_max;
    }

    // Smooth approximation of the max(s,0) function
    double ShockCapturingProcess::SmoothedMaxFunction(const double s)
    {
        const double b = 100;
        const double l_max = (s / Globals::Pi) * std::atan(b * s) + 0.5 * s - (1.0 / Globals::Pi) * std::atan(b) + 0.5;
        return l_max;
    }

    // Smooth approximation of the min(s,0) function
    double ShockCapturingProcess::SmoothedMinFunction(const double s)
    {
        return s - SmoothedMaxFunction(s);
    }


    double ShockCapturingProcess::CalculateProjectedInverseMetricElementSize(
        const Matrix& rInverseMetricTensor,
        const array_1d<double,3>& rScalarGradient)
    {
        double h_proj = 0.0;
        for (std::size_t i = 0; i < rInverseMetricTensor.size1(); ++i) {
            for (std::size_t j = 0; j < rInverseMetricTensor.size2(); ++j) {
                h_proj += rScalarGradient(i) * rInverseMetricTensor(i,j) * rScalarGradient(j);
            }
        }
        return h_proj;
    }

    template<>
    void ShockCapturingProcess::CalculateShockSensorValues<2,3>(
        const Geometry<Node>& rGeometry,
        const array_1d<double,3>& rN,
        const BoundedMatrix<double,3,2>& rDN_DX,
        array_1d<double,3>& rVelocity,
        double& rPressureCoefficient,
        double& rDensity,
        double& rVelocityDivergence,
        array_1d<double,3>& rDensityGradient,
        array_1d<double,3>& rVelocityRotational)
    {
        rVelocityDivergence = 0.0;
        rDensityGradient = ZeroVector(3);
        rVelocityRotational = ZeroVector(3);
        double dvy_dx = 0.0;
        double dvx_dy = 0.0;
        for (std::size_t j = 0; j < 3; ++j) {
            rVelocity = rGeometry[j].GetValue(VELOCITY);
            rPressureCoefficient = rGeometry[j].GetValue(PRESSURE_COEFFICIENT);
            rDensity = rGeometry[j].GetValue(DENSITY);
            dvy_dx += rVelocity(1) * rDN_DX(j,0);
            dvx_dy += rVelocity(0) * rDN_DX(j,1);
            for (std::size_t i = 0; i < 2; ++i) {
                rDensityGradient(i) += rDN_DX(j,i) * rDensity;
                rVelocityDivergence += rDN_DX(j,i) * rVelocity(i);
            }
        }
        rVelocityRotational(2) = dvy_dx - dvx_dy;
    }

    template<>
    void ShockCapturingProcess::CalculateShockSensorValues<3,4>(
        const Geometry<Node>& rGeometry,
        const array_1d<double,4>& rN,
        const BoundedMatrix<double,4,3>& rDN_DX,
        array_1d<double,3>& rVelocity,
        double& rPressureCoefficient,
        double& rDensity,
        double& rVelocityDivergence,
        array_1d<double,3>& rDensityGradient,
        array_1d<double,3>& rVelocityRotational)
    {
        rVelocityDivergence = 0.0;
        rDensityGradient = ZeroVector(3);
        rVelocityRotational = ZeroVector(3);
        double dvx_dy = 0.0;
        double dvx_dz = 0.0;
        double dvy_dx = 0.0;
        double dvy_dz = 0.0;
        double dvz_dx = 0.0;
        double dvz_dy = 0.0;
        for (std::size_t j = 0; j < 4; ++j) {
            rVelocity = rGeometry[j].GetValue(VELOCITY);
            rPressureCoefficient = rGeometry[j].GetValue(PRESSURE_COEFFICIENT);
            rDensity = rGeometry[j].GetValue(DENSITY);
            dvx_dy += rVelocity(0) * rDN_DX(j,1);
            dvx_dz += rVelocity(0) * rDN_DX(j,2);
            dvy_dx += rVelocity(1) * rDN_DX(j,0);
            dvy_dz += rVelocity(1) * rDN_DX(j,2);
            dvz_dx += rVelocity(2) * rDN_DX(j,0);
            dvz_dy += rVelocity(2) * rDN_DX(j,1);
            for (std::size_t i = 0; i < 3; ++i) {
                rDensityGradient(i) += rDN_DX(j,i) * rDensity;
                rVelocityDivergence += rDN_DX(j,i) * rVelocity(i);
            }
        }
        rVelocityRotational(0) = dvz_dy - dvy_dz;
        rVelocityRotational(1) = dvx_dz - dvz_dx;
        rVelocityRotational(2) = dvy_dx - dvx_dy;
    }

    /* External functions *****************************************************/
    /// output stream function
    inline std::ostream& operator << (
        std::ostream& rOStream,
        const ShockCapturingProcess& rThis) {

        rThis.PrintData(rOStream);
        return rOStream;
    }
}
