//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//
//

// System includes

// External includes

// Project includes
#include "geometries/geometry.h"
#include "geometries/geometry_data.h"
#include "processes/calculate_nodal_area_process.h"
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
        Check();
        ExecuteInitialize();
        ExecuteFinalizeSolutionStep();
    }

    void ShockCapturingProcess::ExecuteInitialize()
    {
        // Initialize nodal values
        block_for_each(mrModelPart.Nodes(), [&](Node<3>& rNode) {
            if (mShockSensor) {rNode.SetValue(ARTIFICIAL_BULK_VISCOSITY, 0.0);}
            if (mShearSensor) {rNode.SetValue(ARTIFICIAL_DYNAMIC_VISCOSITY, 0.0);}
            if (mShockSensor || mThermalSensor) {rNode.SetValue(ARTIFICIAL_CONDUCTIVITY, 0.0);}
        });

        // Initialize elemental values
        block_for_each(mrModelPart.Elements(), [&](Element& rElement){
            if (mShockSensor || mThermalSensor) {
                // Artificial conductivity is used in both sensors
                rElement.SetValue(ARTIFICIAL_CONDUCTIVITY, 0.0);
                // Shock sensor specific values
                if (mShockSensor) {
                    rElement.SetValue(SHOCK_SENSOR, 0.0);
                    rElement.SetValue(ARTIFICIAL_BULK_VISCOSITY, 0.0);
                }
                // Thermal sensor specific values
                if (mThermalSensor) {
                    rElement.SetValue(THERMAL_SENSOR, 0.0);
                }
            }
            // Shear sensor specific values
            if (mShearSensor) {
                rElement.SetValue(SHEAR_SENSOR, 0.0);
                rElement.SetValue(ARTIFICIAL_DYNAMIC_VISCOSITY, 0.0);
            }
        });

        // Calculate the NODAL_AREA
        CalculateNodalAreaProcess<false> nodal_area_process(mrModelPart);
        nodal_area_process.Execute();
    }

    void ShockCapturingProcess::ExecuteFinalizeSolutionStep()
    {
        CalculatePhysicsBasedShockCapturing();
    }

    const Parameters ShockCapturingProcess::GetDefaultParameters() const
    {
        const Parameters default_parameters = Parameters(R"({
            "model_part_name" : "",
            "calculate_nodal_area_at_each_step" : false,
            "shock_sensor" : true,
            "shear_sensor" : true,
            "thermal_sensor" : true,
            "thermally_coupled_formulation" : true
        })");

        return default_parameters;
    }

    int ShockCapturingProcess::Check()
    {
        // Base process check
        int err_code = BaseType::Check();

        // Check that the required variables are in the nodal data
        for (const auto& rNode : mrModelPart.Nodes()) {
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DENSITY, rNode)
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY, rNode)
            if (mThermalSensor || mThermallyCoupledFormulation) {
                KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TEMPERATURE, rNode)
            }
        }

        // Check that the required material properties are in the elemental properties
        for (const auto& rElement : mrModelPart.Elements()) {
            const auto& r_prop = rElement.GetProperties();
            KRATOS_ERROR_IF_NOT(r_prop.Has(SPECIFIC_HEAT)) << "Element " << rElement.Id() << " properties " << r_prop.Id() <<  " has no SPECIFIC_HEAT." << std::endl;
            KRATOS_ERROR_IF_NOT(r_prop.Has(HEAT_CAPACITY_RATIO)) << "Element " << rElement.Id() << " properties " << r_prop.Id() << " has no HEAT_CAPACITY_RATIO." << std::endl;
        }

        return err_code;
    }

    /* Private functions *****************************************************/

    void ShockCapturingProcess::ValidateAndAssignParameters(Parameters& rParameters)
    {
        // Validate and assign defaults
        rParameters.ValidateAndAssignDefaults(GetDefaultParameters());

        // Assign settings
        mUpdateNodalArea = rParameters["calculate_nodal_area_at_each_step"].GetBool();
        mShockSensor = rParameters["shock_sensor"].GetBool();
        mShearSensor = rParameters["shear_sensor"].GetBool();
        mThermalSensor = rParameters["thermal_sensor"].GetBool();
        mThermallyCoupledFormulation = rParameters["thermally_coupled_formulation"].GetBool();

        // Check user-provided assigned settings
        KRATOS_ERROR_IF(mThermalSensor && !mThermallyCoupledFormulation)
            << "Thermal sensor cannot be computed using a thermally non-coupled formulation. Check \'thermal_sensor\' and \'thermally_coupled_formulation\' provided settings." << std::endl;
    }

    /**
     * @brief Physics-based shock capturing
     * This function calculates the artificial magnitudes using a physics-based shock capturing method.
     * References https://arc.aiaa.org/doi/abs/10.2514/6.2018-0062
     */
    void ShockCapturingProcess::CalculatePhysicsBasedShockCapturing()
    {
        // Initialize the values to zero
        block_for_each(mrModelPart.Nodes(), [](Node<3> &rNode) {
            rNode.GetValue(ARTIFICIAL_CONDUCTIVITY) = 0.0;
            rNode.GetValue(ARTIFICIAL_BULK_VISCOSITY) = 0.0;
            rNode.GetValue(ARTIFICIAL_DYNAMIC_VISCOSITY) = 0.0;
        });

        // If required, update the NODAL_AREA
        if (mUpdateNodalArea) {
            CalculateNodalAreaProcess<false> nodal_area_process(mrModelPart);
            nodal_area_process.Execute();
        }

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
            KRATOS_ERROR << "Asking for a non-supported geometry. Physics-based shock capturing only supports \'Triangle2D3\' and \'Tetrahedra3D4\' geometries.";
        }

        // Nodal smoothing of the shock capturing magnitudes
        block_for_each(mrModelPart.Nodes(), [](Node<3>& rNode) {
            const double nodal_area = rNode.GetValue(NODAL_AREA);
            rNode.GetValue(ARTIFICIAL_CONDUCTIVITY) /= nodal_area;
            rNode.GetValue(ARTIFICIAL_BULK_VISCOSITY) /= nodal_area;
            rNode.GetValue(ARTIFICIAL_DYNAMIC_VISCOSITY) /= nodal_area;
        });
    }

    double ShockCapturingProcess::LimitingFunction(
        const double s,
        const double s_0,
        const double s_max,
        const double s_min)
    {
        const double aux_1 = std::max(s - s_0, s_min);
        const double aux_2 = std::min(aux_1 - s_max, s_min);
        return aux_2 + s_max;
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
        const Geometry<Node<3>>& rGeometry,
        const array_1d<double,3>& rN,
        const BoundedMatrix<double,3,2>& rDN_DX,
        double& rMachNumber,
        double& rVelocityDivergence,
        array_1d<double,3>& rDensityGradient,
        array_1d<double,3>& rVelocityRotational)
    {
        rMachNumber = 0.0;
        rVelocityDivergence = 0.0;
        rDensityGradient = ZeroVector(3);
        rVelocityRotational = ZeroVector(3);
        double dvy_dx = 0.0;
        double dvx_dy = 0.0;
        for (std::size_t j = 0; j < 3; ++j) {
            rMachNumber = rN(j) * rGeometry[j].GetValue(MACH);
            const auto& r_v_j = rGeometry[j].FastGetSolutionStepValue(VELOCITY);
            const double& r_rho_j = rGeometry[j].FastGetSolutionStepValue(DENSITY);
            dvy_dx += r_v_j(1) * rDN_DX(j,0);
            dvx_dy += r_v_j(0) * rDN_DX(j,1);
            for (std::size_t i = 0; i < 2; ++i) {
                rDensityGradient(i) += rDN_DX(j,i) * r_rho_j;
                rVelocityDivergence += rDN_DX(j,i) * r_v_j(i);
            }
        }
        rVelocityRotational(2) = dvy_dx - dvx_dy;
    }

    template<>
    void ShockCapturingProcess::CalculateShockSensorValues<3,4>(
        const Geometry<Node<3>>& rGeometry,
        const array_1d<double,4>& rN,
        const BoundedMatrix<double,4,3>& rDN_DX,
        double& rMachNumber,
        double& rVelocityDivergence,
        array_1d<double,3>& rDensityGradient,
        array_1d<double,3>& rVelocityRotational)
    {
        rMachNumber = 0.0;
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
            rMachNumber = rN(j) * rGeometry[j].GetValue(MACH);
            const auto& r_v_j = rGeometry[j].FastGetSolutionStepValue(VELOCITY);
            const double& r_rho_j = rGeometry[j].FastGetSolutionStepValue(DENSITY);
            dvx_dy += r_v_j(0) * rDN_DX(j,1);
            dvx_dz += r_v_j(0) * rDN_DX(j,2);
            dvy_dx += r_v_j(1) * rDN_DX(j,0);
            dvy_dz += r_v_j(1) * rDN_DX(j,2);
            dvz_dx += r_v_j(2) * rDN_DX(j,0);
            dvz_dy += r_v_j(2) * rDN_DX(j,1);
            for (std::size_t i = 0; i < 3; ++i) {
                rDensityGradient(i) += rDN_DX(j,i) * r_rho_j;
                rVelocityDivergence += rDN_DX(j,i) * r_v_j(i);
            }
        }
        rVelocityRotational(0) = dvz_dy - dvy_dz;
        rVelocityRotational(1) = dvx_dz - dvz_dx;
        rVelocityRotational(2) = dvy_dx - dvx_dy;
    }

    template<std::size_t TDim, std::size_t TNumNodes>
    void ShockCapturingProcess::CalculateTemperatureGradients(
        const Geometry<Node<3>>& rGeometry,
        const BoundedMatrix<double,TNumNodes,TDim>& rDN_DX,
        const Matrix& rJacobianMatrix,
        array_1d<double,3>& rTemperatureGradient,
        array_1d<double,3>& rTemperatureLocalGradient)
    {
        // Calculate temperature gradient
        array_1d<double,3> grad_temp = ZeroVector(3);
        for (std::size_t j = 0; j < TNumNodes; ++j) {
            const double& r_temp_j = rGeometry[j].FastGetSolutionStepValue(TEMPERATURE);
            for (std::size_t i = 0; i < TDim; ++i) {
                grad_temp(i) += rDN_DX(j,i) * r_temp_j;
            }
        }

        // Calculate temperature local gradient
        rTemperatureLocalGradient = ZeroVector(3);
        for (unsigned int i = 0; i < TDim; ++i) {
            for (unsigned int j = 0; j < TDim; ++j) {
                rTemperatureLocalGradient(i) += rJacobianMatrix(j, i) * grad_temp(j);
            }
        }
    }

    template<std::size_t TDim, std::size_t TNumNodes>
    void ShockCapturingProcess::CalculateShearSensorValues(
        const Geometry<Node<3>>& rGeometry,
        const array_1d<double,TNumNodes>& rN,
        const BoundedMatrix<double,TNumNodes,TDim>& rDN_DX,
        const Matrix& rJacobianMatrix,
        BoundedMatrix<double,TDim,TDim>& rLocalVelocityShearGradient,
        double& rSoundVelocity)
    {
        rSoundVelocity = 0.0;
        BoundedMatrix<double,TDim,TDim> shear_grad = ZeroMatrix(TDim,TDim);
        for (std::size_t k = 0; k < TNumNodes; ++k) {
            const auto& r_vel = rGeometry[k].FastGetSolutionStepValue(VELOCITY);
            rSoundVelocity += rN(k) * rGeometry[k].GetValue(SOUND_VELOCITY);
            for (std::size_t i = 0; i < TDim; ++i) {
                for (std::size_t j = 0; j < TDim; ++j) {
                    shear_grad(i,j) += rDN_DX(k,j) * r_vel(i);
                }
            }
        }

        // Only keep the shear component of the velocity gradient
        for (unsigned int d = 0; d < TDim; ++d) {
            shear_grad(d, d) = 0.0;
        }

        // Multiply by the Jacobian matrix to obtain the local gradient
        rLocalVelocityShearGradient = prod(shear_grad, trans(rJacobianMatrix));
    }

    /* External functions *****************************************************/

    /// output stream function
    inline std::ostream& operator << (
        std::ostream& rOStream,
        const ShockCapturingProcess& rThis) {

        rThis.PrintData(rOStream);
        return rOStream;
    }

    /* Explicit template instantiation ****************************************/
    template void KRATOS_API(FLUID_DYNAMICS_APPLICATION) ShockCapturingProcess::CalculateTemperatureGradients<2,3>(
        const Geometry<Node<3>>& rGeometry,
        const BoundedMatrix<double,3,2>& rDN_DX,
        const Matrix& rJacobianMatrix,
        array_1d<double,3>& rTemperatureGradient,
        array_1d<double,3>& rTemperatureLocalGradient);

    template void KRATOS_API(FLUID_DYNAMICS_APPLICATION) ShockCapturingProcess::CalculateTemperatureGradients<3,4>(
        const Geometry<Node<3>>& rGeometry,
        const BoundedMatrix<double,4,3>& rDN_DX,
        const Matrix& rJacobianMatrix,
        array_1d<double,3>& rTemperatureGradient,
        array_1d<double,3>& rTemperatureLocalGradient);

    template void KRATOS_API(FLUID_DYNAMICS_APPLICATION) ShockCapturingProcess::CalculateShearSensorValues<2,3>(
        const Geometry<Node<3>>& rGeometry,
        const array_1d<double,3>& rN,
        const BoundedMatrix<double,3,2>& rDN_DX,
        const Matrix& rJacobianMatrix,
        BoundedMatrix<double,2,2>& rLocalVelocityShearGradient,
        double& rSoundVelocity);

    template void KRATOS_API(FLUID_DYNAMICS_APPLICATION) ShockCapturingProcess::CalculateShearSensorValues<3,4>(
        const Geometry<Node<3>>& rGeometry,
        const array_1d<double,4>& rN,
        const BoundedMatrix<double,4,3>& rDN_DX,
        const Matrix& rJacobianMatrix,
        BoundedMatrix<double,3,3>& rLocalVelocityShearGradient,
        double& rSoundVelocity);

}
