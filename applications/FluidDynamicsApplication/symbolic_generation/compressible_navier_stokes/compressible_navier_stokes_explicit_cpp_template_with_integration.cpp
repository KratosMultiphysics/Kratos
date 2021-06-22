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

// System includes


// External includes


// Project includes
#include "includes/checks.h"
#include "utilities/element_size_calculator.h"

// Application includes
#include "custom_elements/compressible_navier_stokes_explicit.h"


namespace Kratos {

template <>
void CompressibleNavierStokesExplicit<2,3>::EquationIdVector(
    EquationIdVectorType &rResult,
    const ProcessInfo &rCurrentProcessInfo) const
{
    KRATOS_TRY

    constexpr unsigned int n_nodes = 3;
    constexpr unsigned int block_size = 4;
    constexpr unsigned int dof_size = n_nodes * block_size;

    if (rResult.size() != dof_size) {
        rResult.resize(dof_size);
    }

    unsigned int local_index = 0;
    const auto& r_geometry = GetGeometry();
    const unsigned int den_pos = r_geometry[0].GetDofPosition(DENSITY);
    const unsigned int mom_pos = r_geometry[0].GetDofPosition(MOMENTUM);
    const unsigned int enr_pos = r_geometry[0].GetDofPosition(TOTAL_ENERGY);
    for (unsigned int i_node = 0; i_node < n_nodes; ++i_node) {
        rResult[local_index++] = r_geometry[i_node].GetDof(DENSITY, den_pos).EquationId();
        rResult[local_index++] = r_geometry[i_node].GetDof(MOMENTUM_X, mom_pos).EquationId();
        rResult[local_index++] = r_geometry[i_node].GetDof(MOMENTUM_Y, mom_pos + 1).EquationId();
        rResult[local_index++] = r_geometry[i_node].GetDof(TOTAL_ENERGY, enr_pos).EquationId();
    }

    KRATOS_CATCH("");
}

template <>
void CompressibleNavierStokesExplicit<3,4>::EquationIdVector(
    EquationIdVectorType &rResult,
    const ProcessInfo &rCurrentProcessInfo) const
{
    KRATOS_TRY

    constexpr unsigned int n_nodes = 4;
    constexpr unsigned int block_size = 5;
    constexpr unsigned int dof_size = n_nodes * block_size;

    if (rResult.size() != dof_size) {
        rResult.resize(dof_size);
    }

    unsigned int local_index = 0;
    const auto& r_geometry = GetGeometry();
    const unsigned int den_pos = r_geometry[0].GetDofPosition(DENSITY);
    const unsigned int mom_pos = r_geometry[0].GetDofPosition(MOMENTUM);
    const unsigned int enr_pos = r_geometry[0].GetDofPosition(TOTAL_ENERGY);
    for (unsigned int i_node = 0; i_node < n_nodes; ++i_node) {
        rResult[local_index++] = r_geometry[i_node].GetDof(DENSITY, den_pos).EquationId();
        rResult[local_index++] = r_geometry[i_node].GetDof(MOMENTUM_X, mom_pos).EquationId();
        rResult[local_index++] = r_geometry[i_node].GetDof(MOMENTUM_Y, mom_pos + 1).EquationId();
        rResult[local_index++] = r_geometry[i_node].GetDof(MOMENTUM_Z, mom_pos + 2).EquationId();
        rResult[local_index++] = r_geometry[i_node].GetDof(TOTAL_ENERGY, enr_pos).EquationId();
    }

    KRATOS_CATCH("");
}

template <>
void CompressibleNavierStokesExplicit<2,3>::GetDofList(
    DofsVectorType &ElementalDofList,
    const ProcessInfo &rCurrentProcessInfo) const
{
    KRATOS_TRY

    constexpr unsigned int n_nodes = 3;
    constexpr unsigned int block_size = 4;
    constexpr unsigned int dof_size = n_nodes * block_size;

    if (ElementalDofList.size() != dof_size) {
        ElementalDofList.resize(dof_size);
    }

    unsigned int local_index = 0;
    const auto& r_geometry = GetGeometry();
    const unsigned int den_pos = r_geometry[0].GetDofPosition(DENSITY);
    const unsigned int mom_pos = r_geometry[0].GetDofPosition(MOMENTUM);
    const unsigned int enr_pos = r_geometry[0].GetDofPosition(TOTAL_ENERGY);
    for (unsigned int i_node = 0; i_node < n_nodes; ++i_node) {
        ElementalDofList[local_index++] = r_geometry[i_node].pGetDof(DENSITY, den_pos);
        ElementalDofList[local_index++] = r_geometry[i_node].pGetDof(MOMENTUM_X, mom_pos);
        ElementalDofList[local_index++] = r_geometry[i_node].pGetDof(MOMENTUM_Y, mom_pos + 1);
        ElementalDofList[local_index++] = r_geometry[i_node].pGetDof(TOTAL_ENERGY, enr_pos);
    }

    KRATOS_CATCH("");
}

template <>
void CompressibleNavierStokesExplicit<3,4>::GetDofList(
    DofsVectorType &ElementalDofList,
    const ProcessInfo &rCurrentProcessInfo) const
{
    KRATOS_TRY

    constexpr unsigned int n_nodes = 4;
    constexpr unsigned int block_size = 5;
    unsigned int dof_size = n_nodes * block_size;

    if (ElementalDofList.size() != dof_size) {
        ElementalDofList.resize(dof_size);
    }

    unsigned int local_index = 0;
    const auto &r_geometry = GetGeometry();
    const unsigned int den_pos = r_geometry[0].GetDofPosition(DENSITY);
    const unsigned int mom_pos = r_geometry[0].GetDofPosition(MOMENTUM);
    const unsigned int enr_pos = r_geometry[0].GetDofPosition(TOTAL_ENERGY);
    for (unsigned int i_node = 0; i_node < n_nodes; ++i_node) {
        ElementalDofList[local_index++] = this->GetGeometry()[i_node].pGetDof(DENSITY, den_pos);
        ElementalDofList[local_index++] = this->GetGeometry()[i_node].pGetDof(MOMENTUM_X, mom_pos);
        ElementalDofList[local_index++] = this->GetGeometry()[i_node].pGetDof(MOMENTUM_Y, mom_pos + 1);
        ElementalDofList[local_index++] = this->GetGeometry()[i_node].pGetDof(MOMENTUM_Z, mom_pos + 2);
        ElementalDofList[local_index++] = this->GetGeometry()[i_node].pGetDof(TOTAL_ENERGY, enr_pos);
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
int CompressibleNavierStokesExplicit<TDim, TNumNodes>::Check(const ProcessInfo &rCurrentProcessInfo) const
{
    KRATOS_TRY

    // Perform basic element checks
    int ErrorCode = Kratos::Element::Check(rCurrentProcessInfo);
    if (ErrorCode != 0) {
        return ErrorCode;
    }

    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    for (unsigned int i = 0; i < TNumNodes; ++i) {
        KRATOS_ERROR_IF_NOT(this->GetGeometry()[i].SolutionStepsDataHas(DENSITY)) << "Missing DENSITY variable on solution step data for node " << this->GetGeometry()[i].Id();
        KRATOS_ERROR_IF_NOT(this->GetGeometry()[i].SolutionStepsDataHas(MOMENTUM)) << "Missing MOMENTUM variable on solution step data for node " << this->GetGeometry()[i].Id();
        KRATOS_ERROR_IF_NOT(this->GetGeometry()[i].SolutionStepsDataHas(TOTAL_ENERGY)) << "Missing TOTAL_ENERGY variable on solution step data for node " << this->GetGeometry()[i].Id();
        KRATOS_ERROR_IF_NOT(this->GetGeometry()[i].SolutionStepsDataHas(BODY_FORCE)) << "Missing BODY_FORCE variable on solution step data for node " << this->GetGeometry()[i].Id();
        KRATOS_ERROR_IF_NOT(this->GetGeometry()[i].SolutionStepsDataHas(HEAT_SOURCE)) << "Missing HEAT_SOURCE variable on solution step data for node " << this->GetGeometry()[i].Id();

        // Activate as soon as we start using the explicit DOF based strategy
        KRATOS_ERROR_IF_NOT(this->GetGeometry()[i].HasDofFor(DENSITY)) << "Missing DENSITY DOF in node ", this->GetGeometry()[i].Id();
        KRATOS_ERROR_IF_NOT(this->GetGeometry()[i].HasDofFor(MOMENTUM_X) || this->GetGeometry()[i].HasDofFor(MOMENTUM_Y)) << "Missing MOMENTUM component DOF in node ", this->GetGeometry()[i].Id();
        if (TDim == 3) {
            KRATOS_ERROR_IF_NOT(this->GetGeometry()[i].HasDofFor(MOMENTUM_Z)) << "Missing MOMENTUM component DOF in node ", this->GetGeometry()[i].Id();
        }
        KRATOS_ERROR_IF_NOT(this->GetGeometry()[i].HasDofFor(TOTAL_ENERGY)) << "Missing TOTAL_ENERGY DOF in node ", this->GetGeometry()[i].Id();
    }

    return 0;

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void CompressibleNavierStokesExplicit<TDim, TNumNodes>::Calculate(
    const Variable<double>& rVariable,
    double& Output,
    const ProcessInfo& rCurrentProcessInfo)
{
    // Lumped projection terms
    if (rVariable == DENSITY_PROJECTION) {
        CalculateDensityProjection(rCurrentProcessInfo);
    } else if (rVariable == TOTAL_ENERGY_PROJECTION) {
        CalculateTotalEnergyProjection(rCurrentProcessInfo);
    } else if (rVariable == VELOCITY_DIVERGENCE) {
        Output = CalculateMidPointVelocityDivergence();
    } else if (rVariable == SOUND_VELOCITY) {
        Output = CalculateMidPointSoundVelocity();
    } else {
        KRATOS_ERROR << "Variable not implemented." << std::endl;
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void CompressibleNavierStokesExplicit<TDim, TNumNodes>::Calculate(
    const Variable<array_1d<double, 3 > >& rVariable,
    array_1d<double, 3 > & Output,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rVariable == DENSITY_GRADIENT) {
        Output = CalculateMidPointDensityGradient();
    } else if (rVariable == TEMPERATURE_GRADIENT) {
        Output = CalculateMidPointTemperatureGradient();
    } else if (rVariable == VELOCITY_ROTATIONAL) {
        Output = CalculateMidPointVelocityRotational();
    } else if (rVariable == MOMENTUM_PROJECTION) {
        CalculateMomentumProjection(rCurrentProcessInfo);
    } else {
        KRATOS_ERROR << "Variable not implemented." << std::endl;
    }
}


template <unsigned int TDim, unsigned int TNumNodes>
void CompressibleNavierStokesExplicit<TDim, TNumNodes>::Calculate(
    const Variable<Matrix>& rVariable,
    Matrix & Output,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rVariable == VELOCITY_GRADIENT) {
        Output = CalculateMidPointVelocityGradient();
    } else {
        KRATOS_ERROR << "Variable not implemented." << std::endl;
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void CompressibleNavierStokesExplicit<TDim, TNumNodes>::CalculateOnIntegrationPoints(
    const Variable<double>& rVariable,
    std::vector<double>& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    const auto& r_geometry = GetGeometry();
    const auto& r_integration_points = r_geometry.IntegrationPoints();
    if (rOutput.size() != r_integration_points.size()) {
        rOutput.resize( r_integration_points.size() );
    }

    if (rVariable == SHOCK_SENSOR) {
        const double sc = this->GetValue(SHOCK_SENSOR);
        for (unsigned int i_gauss = 0; i_gauss < r_integration_points.size(); ++i_gauss) {
            rOutput[i_gauss] = sc;
        }
    } else if (rVariable == SHEAR_SENSOR) {
        const double sc = this->GetValue(SHEAR_SENSOR);
        for (unsigned int i_gauss = 0; i_gauss < r_integration_points.size(); ++i_gauss) {
            rOutput[i_gauss] = sc;
        }
    } else if (rVariable == THERMAL_SENSOR) {
        const double sc = this->GetValue(THERMAL_SENSOR);
        for (unsigned int i_gauss = 0; i_gauss < r_integration_points.size(); ++i_gauss) {
            rOutput[i_gauss] = sc;
        }
    } else if (rVariable == ARTIFICIAL_CONDUCTIVITY) {
        const double k_star = this->GetValue(ARTIFICIAL_CONDUCTIVITY);
        for (unsigned int i_gauss = 0; i_gauss < r_integration_points.size(); ++i_gauss) {
            rOutput[i_gauss] = k_star;
        }
    } else if (rVariable == ARTIFICIAL_BULK_VISCOSITY) {
        const double beta_star = this->GetValue(ARTIFICIAL_BULK_VISCOSITY);
        for (unsigned int i_gauss = 0; i_gauss < r_integration_points.size(); ++i_gauss) {
            rOutput[i_gauss] = beta_star;
        }
    } else if (rVariable == VELOCITY_DIVERGENCE) {
        const double div_v = CalculateMidPointVelocityDivergence();
        for (unsigned int i_gauss = 0; i_gauss < r_integration_points.size(); ++i_gauss) {
            rOutput[i_gauss] = div_v;
        }
    } else {
        KRATOS_ERROR << "Variable not implemented." << std::endl;
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void CompressibleNavierStokesExplicit<TDim, TNumNodes>::CalculateOnIntegrationPoints(
    const Variable<array_1d<double,3>>& rVariable,
    std::vector<array_1d<double,3>>& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    const auto& r_geometry = GetGeometry();
    const auto& r_integration_points = r_geometry.IntegrationPoints();
    if (rOutput.size() != r_integration_points.size()) {
        rOutput.resize( r_integration_points.size() );
    }

    if (rVariable == DENSITY_GRADIENT) {
        const array_1d<double,3> rho_grad = CalculateMidPointDensityGradient();
        for (unsigned int i_gauss = 0; i_gauss < r_integration_points.size(); ++i_gauss) {
            rOutput[i_gauss] = rho_grad;
        }
    } else if (rVariable == TEMPERATURE_GRADIENT) {
        const array_1d<double,3> temp_grad = CalculateMidPointTemperatureGradient();

        for (unsigned int i_gauss = 0; i_gauss < r_integration_points.size(); ++i_gauss) {
            rOutput[i_gauss] = temp_grad;
        }
    } else if (rVariable == VELOCITY_ROTATIONAL) {
        const array_1d<double,3> rot_v = CalculateMidPointVelocityRotational();
        for (unsigned int i_gauss = 0; i_gauss < r_integration_points.size(); ++i_gauss) {
            rOutput[i_gauss] = rot_v;
        }
    } else {
        KRATOS_ERROR << "Variable not implemented." << std::endl;
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void CompressibleNavierStokesExplicit<TDim, TNumNodes>::FillElementData(
    ElementDataStruct &rData,
    const ProcessInfo &rCurrentProcessInfo)
{
    // Getting data for the given geometry
    const auto& r_geometry = GetGeometry();
    GeometryUtils::CalculateGeometryData(r_geometry, rData.DN_DX, rData.N, rData.volume);

    // Compute element size
    rData.h = ElementSizeCalculator<TDim, TNumNodes>::GradientsElementSize(rData.DN_DX);

    // Database access to all of the variables needed
    Properties &r_properties = this->GetProperties();
    rData.mu = r_properties.GetValue(DYNAMIC_VISCOSITY);
    rData.lambda = r_properties.GetValue(CONDUCTIVITY);
    rData.c_v = r_properties.GetValue(SPECIFIC_HEAT); // TODO: WE SHOULD SPECIFY WHICH ONE --> CREATE SPECIFIC_HEAT_CONSTANT_VOLUME
    rData.gamma = r_properties.GetValue(HEAT_CAPACITY_RATIO);

    rData.UseOSS = rCurrentProcessInfo[OSS_SWITCH];
    rData.ShockCapturing = rCurrentProcessInfo[SHOCK_CAPTURING_SWITCH];

    // Magnitudes to calculate the time derivatives
    const double time_step = rCurrentProcessInfo[DELTA_TIME];
    const double theta = rCurrentProcessInfo[TIME_INTEGRATION_THETA];
    const double aux_theta = theta > 0 ? 1.0 / (theta * time_step) : 0.0;

    // Get nodal values
    if (rData.UseOSS) {
        for (unsigned int i = 0; i < TNumNodes; ++i) {
            const auto& r_node = r_geometry[i];
            // Vector data
            const array_1d<double,3>& r_momentum = r_node.FastGetSolutionStepValue(MOMENTUM);
            const array_1d<double,3>& r_momentum_old = r_node.FastGetSolutionStepValue(MOMENTUM, 1);
            const array_1d<double,3>& r_momentum_projection = r_node.GetValue(MOMENTUM_PROJECTION);
            const array_1d<double,3> mom_inc = r_momentum - r_momentum_old;
            const auto& r_body_force = r_node.FastGetSolutionStepValue(BODY_FORCE);
            for (unsigned int k = 0; k < TDim; ++k) {
                rData.U(i, k + 1) = r_momentum[k];
                rData.dUdt(i, k + 1) = aux_theta * mom_inc[k];
                rData.ResProj(i, k + 1) = r_momentum_projection[k];
                rData.f_ext(i, k) = r_body_force[k];
            }
            // Density data
            const double& r_rho = r_node.FastGetSolutionStepValue(DENSITY);
            const double& r_rho_old = r_node.FastGetSolutionStepValue(DENSITY, 1);
            const double rho_inc = r_rho - r_rho_old;
            rData.U(i, 0) = r_rho;
            rData.dUdt(i, 0) = aux_theta * rho_inc;
            rData.ResProj(i, 0) = r_node.GetValue(DENSITY_PROJECTION);
            // Total energy data
            const double& r_tot_ener = r_node.FastGetSolutionStepValue(TOTAL_ENERGY);
            const double& r_tot_ener_old = r_node.FastGetSolutionStepValue(TOTAL_ENERGY, 1);
            const double tot_ener_inc = r_tot_ener - r_tot_ener_old;
            rData.U(i, TDim + 1) = r_tot_ener;
            rData.dUdt(i, TDim + 1) = aux_theta * tot_ener_inc;
            rData.ResProj(i, TDim + 1) = r_node.GetValue(TOTAL_ENERGY_PROJECTION);
            // Source data
            rData.r_ext(i) = r_node.FastGetSolutionStepValue(HEAT_SOURCE);
            rData.m_ext(i) = r_node.FastGetSolutionStepValue(MASS_SOURCE);
            // Shock capturing data
            rData.mu_sc_nodes(i) = r_node.GetValue(ARTIFICIAL_DYNAMIC_VISCOSITY);
            rData.beta_sc_nodes(i) = r_node.GetValue(ARTIFICIAL_BULK_VISCOSITY);
            rData.lamb_sc_nodes(i) = r_node.GetValue(ARTIFICIAL_CONDUCTIVITY);
        }
    } else {
        for (unsigned int i = 0; i < TNumNodes; ++i) {
            const auto& r_node = r_geometry[i];
            // Vector data
            const array_1d<double,3>& r_momentum = r_node.FastGetSolutionStepValue(MOMENTUM);
            const array_1d<double,3>& r_momentum_old = r_node.FastGetSolutionStepValue(MOMENTUM, 1);
            const array_1d<double,3> mom_inc = r_momentum - r_momentum_old;
            const auto& r_body_force = r_node.FastGetSolutionStepValue(BODY_FORCE);
            for (unsigned int k = 0; k < TDim; ++k) {
                rData.U(i, k + 1) = r_momentum[k];
                rData.dUdt(i, k + 1) = aux_theta * mom_inc[k];
                rData.f_ext(i, k) = r_body_force[k];
            }
            // Density data
            const double& r_rho = r_node.FastGetSolutionStepValue(DENSITY);
            const double& r_rho_old = r_node.FastGetSolutionStepValue(DENSITY, 1);
            rData.U(i, 0) = r_rho;
            rData.dUdt(i, 0) = aux_theta * (r_rho - r_rho_old);
            // Total energy data
            const double& r_tot_ener = r_node.FastGetSolutionStepValue(TOTAL_ENERGY);
            const double& r_tot_ener_old = r_node.FastGetSolutionStepValue(TOTAL_ENERGY, 1);
            rData.U(i, TDim + 1) = r_tot_ener;
            rData.dUdt(i, TDim + 1) = aux_theta * (r_tot_ener - r_tot_ener_old);
            // Source data
            rData.r_ext(i) = r_node.FastGetSolutionStepValue(HEAT_SOURCE);
            rData.m_ext(i) = r_node.FastGetSolutionStepValue(MASS_SOURCE);
            // Shock capturing data
            rData.mu_sc_nodes(i) = r_node.GetValue(ARTIFICIAL_DYNAMIC_VISCOSITY);
            rData.beta_sc_nodes(i) = r_node.GetValue(ARTIFICIAL_BULK_VISCOSITY);
            rData.lamb_sc_nodes(i) = r_node.GetValue(ARTIFICIAL_CONDUCTIVITY);
        }
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
array_1d<double,3> CompressibleNavierStokesExplicit<TDim, TNumNodes>::CalculateMidPointDensityGradient() const
{
    // Get geometry data
    const auto& r_geom = GetGeometry();
    const unsigned int n_nodes = r_geom.PointsNumber();
    Geometry<Node<3>>::ShapeFunctionsGradientsType dNdX_container;
    r_geom.ShapeFunctionsIntegrationPointsGradients(dNdX_container, GeometryData::GI_GAUSS_1);
    const auto& r_dNdX = dNdX_container[0];

    // Calculate midpoint magnitudes
    array_1d<double,3> midpoint_grad_rho = ZeroVector(3);
    for (unsigned int i_node = 0; i_node < n_nodes; ++i_node) {
        auto& r_node = r_geom[i_node];
        const auto node_dNdX = row(r_dNdX, i_node);
        const double& r_rho = r_node.FastGetSolutionStepValue(DENSITY);
        for (unsigned int d1 = 0; d1 < TDim; ++d1) {
            midpoint_grad_rho[d1] += node_dNdX(d1) * r_rho;
        }
    }

    return midpoint_grad_rho;
}

// TODO: IN HERE WE HAVE LINEARIZED THE DERIVATIVE... CHECK IF WE SHOULD PROPERLY COMPUTE IT
// TODO: IN HERE WE HAVE LINEARIZED THE DERIVATIVE... CHECK IF WE SHOULD PROPERLY COMPUTE IT
// TODO: IN HERE WE HAVE LINEARIZED THE DERIVATIVE... CHECK IF WE SHOULD PROPERLY COMPUTE IT
// TODO: IN HERE WE HAVE LINEARIZED THE DERIVATIVE... CHECK IF WE SHOULD PROPERLY COMPUTE IT
// TODO: IN HERE WE HAVE LINEARIZED THE DERIVATIVE... CHECK IF WE SHOULD PROPERLY COMPUTE IT
// TODO: IN HERE WE HAVE LINEARIZED THE DERIVATIVE... CHECK IF WE SHOULD PROPERLY COMPUTE IT
// TODO: IN HERE WE HAVE LINEARIZED THE DERIVATIVE... CHECK IF WE SHOULD PROPERLY COMPUTE IT
template <unsigned int TDim, unsigned int TNumNodes>
array_1d<double,3> CompressibleNavierStokesExplicit<TDim, TNumNodes>::CalculateMidPointTemperatureGradient() const
{
    // Get geometry data
    const auto& r_geom = GetGeometry();
    const unsigned int n_nodes = r_geom.PointsNumber();
    Geometry<Node<3>>::ShapeFunctionsGradientsType dNdX_container;
    r_geom.ShapeFunctionsIntegrationPointsGradients(dNdX_container, GeometryData::GI_GAUSS_1);
    const auto& r_dNdX = dNdX_container[0];

    // Calculate midpoint magnitudes
    const double c_v = GetProperties()[SPECIFIC_HEAT];
    array_1d<double,3> midpoint_grad_temp = ZeroVector(3);
    for (unsigned int i_node = 0; i_node < n_nodes; ++i_node) {
        auto& r_node = r_geom[i_node];
        const auto node_dNdX = row(r_dNdX, i_node);
        const auto& r_mom = r_node.FastGetSolutionStepValue(MOMENTUM);
        const double& r_rho = r_node.FastGetSolutionStepValue(DENSITY);
        const double& r_tot_ener = r_node.FastGetSolutionStepValue(TOTAL_ENERGY);
        const array_1d<double, 3> vel = r_mom / r_rho;
        const double temp = (r_tot_ener / r_rho - 0.5 * inner_prod(vel, vel)) / c_v;
        for (unsigned int d1 = 0; d1 < TDim; ++d1) {
            midpoint_grad_temp[d1] += node_dNdX(d1) * temp;
        }
    }

    return midpoint_grad_temp;
}

template <unsigned int TDim, unsigned int TNumNodes>
double CompressibleNavierStokesExplicit<TDim, TNumNodes>::CalculateMidPointSoundVelocity() const
{
    // Get geometry data
    const auto& r_geom = GetGeometry();
    const unsigned int n_nodes = r_geom.PointsNumber();

    // Calculate midpoint magnitudes
    double midpoint_rho = 0.0;
    double midpoint_tot_ener = 0.0;
    array_1d<double,TDim> midpoint_mom = ZeroVector(TDim);
    for (unsigned int i_node = 0; i_node < n_nodes; ++i_node) {
        auto& r_node = r_geom[i_node];
        const auto& r_mom = r_node.FastGetSolutionStepValue(MOMENTUM);
        const double& r_rho = r_node.FastGetSolutionStepValue(DENSITY);
        const double& r_tot_ener = r_node.FastGetSolutionStepValue(TOTAL_ENERGY);
        midpoint_rho += r_rho;
        midpoint_tot_ener += r_tot_ener;
        for (unsigned int d1 = 0; d1 < TDim; ++d1) {
            midpoint_mom[d1] += r_mom(d1);
        }
    }
    midpoint_rho /= n_nodes;
    midpoint_mom /= n_nodes;
    midpoint_tot_ener /= n_nodes;

    // Calculate midpoint speed of sound
    const auto& r_prop = GetProperties();
    const double c_v = r_prop.GetValue(SPECIFIC_HEAT);
    const double gamma = r_prop.GetValue(HEAT_CAPACITY_RATIO);
    const double temp = (midpoint_tot_ener / midpoint_rho - inner_prod(midpoint_mom, midpoint_mom) / (2 * std::pow(midpoint_rho, 2))) / c_v;
    double midpoint_c = std::sqrt(gamma * (gamma - 1.0) * c_v * temp);
    return midpoint_c;
}

template <unsigned int TDim, unsigned int TNumNodes>
double CompressibleNavierStokesExplicit<TDim, TNumNodes>::CalculateMidPointVelocityDivergence() const
{
    // Get geometry data
    const auto& r_geom = GetGeometry();
    const unsigned int n_nodes = r_geom.PointsNumber();
    Geometry<Node<3>>::ShapeFunctionsGradientsType dNdX_container;
    r_geom.ShapeFunctionsIntegrationPointsGradients(dNdX_container, GeometryData::GI_GAUSS_1);
    const auto& r_dNdX = dNdX_container[0];

    // Calculate midpoint magnitudes
    double midpoint_rho = 0.0;
    double midpoint_div_mom = 0.0;
    array_1d<double,TDim> midpoint_mom = ZeroVector(TDim);
    array_1d<double,TDim> midpoint_grad_rho = ZeroVector(TDim);
    for (unsigned int i_node = 0; i_node < n_nodes; ++i_node) {
        auto& r_node = r_geom[i_node];
        const auto node_dNdX = row(r_dNdX, i_node);
        const auto& r_mom = r_node.FastGetSolutionStepValue(MOMENTUM);
        const double& r_rho = r_node.FastGetSolutionStepValue(DENSITY);
        midpoint_rho += r_rho;
        for (unsigned int d1 = 0; d1 < TDim; ++d1) {
            midpoint_mom[d1] += r_mom(d1);
            midpoint_div_mom += node_dNdX(d1) * r_mom(d1);
            midpoint_grad_rho[d1] += node_dNdX(d1) * r_rho;
        }
    }
    midpoint_rho /= n_nodes;
    midpoint_mom /= n_nodes;

    // Calculate velocity divergence
    // Note that the formulation is written in conservative variables. Hence we do div(mom/rho).
    double midpoint_div_v = (midpoint_rho * midpoint_div_mom - inner_prod(midpoint_mom, midpoint_grad_rho)) / std::pow(midpoint_rho, 2);
    return midpoint_div_v;
}

template <>
array_1d<double,3> CompressibleNavierStokesExplicit<2,3>::CalculateMidPointVelocityRotational() const
{
    // Get geometry data
    const auto& r_geom = GetGeometry();
    const unsigned int n_nodes = r_geom.PointsNumber();
    Geometry<Node<3>>::ShapeFunctionsGradientsType dNdX_container;
    r_geom.ShapeFunctionsIntegrationPointsGradients(dNdX_container, GeometryData::GI_GAUSS_1);
    const auto& r_dNdX = dNdX_container[0];

    // Calculate midpoint magnitudes
    double midpoint_rho = 0.0;
    double midpoint_dmy_dx = 0.0;
    double midpoint_dmx_dy = 0.0;
    double midpoint_rho_dx = 0.0;
    double midpoint_rho_dy = 0.0;
    array_1d<double,3> midpoint_mom = ZeroVector(3);
    for (unsigned int i_node = 0; i_node < n_nodes; ++i_node) {
        auto& r_node = r_geom[i_node];
        const auto node_dNdX = row(r_dNdX, i_node);
        const auto& r_mom = r_node.FastGetSolutionStepValue(MOMENTUM);
        const double& r_rho = r_node.FastGetSolutionStepValue(DENSITY);
        midpoint_rho += r_rho;
        midpoint_mom += r_mom;
        midpoint_dmy_dx += r_mom[1] * node_dNdX[0];
        midpoint_dmx_dy += r_mom[0] * node_dNdX[1];
        midpoint_rho_dx += r_rho * node_dNdX[0];
        midpoint_rho_dy += r_rho * node_dNdX[1];
    }
    midpoint_rho /= n_nodes;
    midpoint_mom /= n_nodes;

    // Calculate velocity rotational
    // Note that the formulation is written in conservative variables. Hence we do rot(mom/rho).
    const double dvy_dx = (midpoint_dmy_dx * midpoint_rho - midpoint_mom[1] * midpoint_rho_dx) / std::pow(midpoint_rho, 2);
    const double dvx_dy = (midpoint_dmx_dy * midpoint_rho - midpoint_mom[0] * midpoint_rho_dy) / std::pow(midpoint_rho, 2);
    array_1d<double,3> midpoint_rot_v;
    midpoint_rot_v[0] = 0.0;
    midpoint_rot_v[1] = 0.0;
    midpoint_rot_v[2] = dvy_dx - dvx_dy;
    return midpoint_rot_v;
}

template <>
array_1d<double,3> CompressibleNavierStokesExplicit<3,4>::CalculateMidPointVelocityRotational() const
{
    // Get geometry data
    const auto& r_geom = GetGeometry();
    const unsigned int n_nodes = r_geom.PointsNumber();
    Geometry<Node<3>>::ShapeFunctionsGradientsType dNdX_container;
    r_geom.ShapeFunctionsIntegrationPointsGradients(dNdX_container, GeometryData::GI_GAUSS_1);
    const auto& r_dNdX = dNdX_container[0];

    // Calculate midpoint magnitudes
    double midpoint_rho = 0.0;
    double midpoint_dmx_dy = 0.0;
    double midpoint_dmx_dz = 0.0;
    double midpoint_dmy_dx = 0.0;
    double midpoint_dmy_dz = 0.0;
    double midpoint_dmz_dx = 0.0;
    double midpoint_dmz_dy = 0.0;
    double midpoint_rho_dx = 0.0;
    double midpoint_rho_dy = 0.0;
    double midpoint_rho_dz = 0.0;
    array_1d<double,3> midpoint_mom = ZeroVector(3);
    for (unsigned int i_node = 0; i_node < n_nodes; ++i_node) {
        auto& r_node = r_geom[i_node];
        const auto node_dNdX = row(r_dNdX, i_node);
        const auto& r_mom = r_node.FastGetSolutionStepValue(MOMENTUM);
        const double& r_rho = r_node.FastGetSolutionStepValue(DENSITY);
        midpoint_rho += r_rho;
        midpoint_mom += r_mom;
        midpoint_dmx_dy += r_mom[0] * node_dNdX[1];
        midpoint_dmx_dz += r_mom[0] * node_dNdX[2];
        midpoint_dmy_dx += r_mom[1] * node_dNdX[0];
        midpoint_dmy_dz += r_mom[1] * node_dNdX[2];
        midpoint_dmz_dx += r_mom[2] * node_dNdX[0];
        midpoint_dmz_dy += r_mom[2] * node_dNdX[1];
        midpoint_rho_dx += r_rho * node_dNdX[0];
        midpoint_rho_dy += r_rho * node_dNdX[1];
        midpoint_rho_dz += r_rho * node_dNdX[2];
    }
    midpoint_rho /= n_nodes;
    midpoint_mom /= n_nodes;

    // Calculate velocity rotational
    // Note that the formulation is written in conservative variables. Hence we do rot(mom/rho).
    const double rho_pow = std::pow(midpoint_rho, 2);
    const double dvz_dy = (midpoint_dmz_dy * midpoint_rho - midpoint_mom[2] * midpoint_rho_dy) / rho_pow;
    const double dvy_dz = (midpoint_dmy_dz * midpoint_rho - midpoint_mom[1] * midpoint_rho_dz) / rho_pow;
    const double dvx_dz = (midpoint_dmx_dz * midpoint_rho - midpoint_mom[0] * midpoint_rho_dz) / rho_pow;
    const double dvz_dx = (midpoint_dmz_dx * midpoint_rho - midpoint_mom[2] * midpoint_rho_dx) / rho_pow;
    const double dvy_dx = (midpoint_dmy_dx * midpoint_rho - midpoint_mom[1] * midpoint_rho_dx) / rho_pow;
    const double dvx_dy = (midpoint_dmx_dy * midpoint_rho - midpoint_mom[0] * midpoint_rho_dy) / rho_pow;
    array_1d<double,3> midpoint_rot_v;
    midpoint_rot_v[0] = dvz_dy - dvy_dz;
    midpoint_rot_v[1] = dvx_dz - dvz_dx;
    midpoint_rot_v[2] = dvy_dx - dvx_dy;
    return midpoint_rot_v;
}

template <>
BoundedMatrix<double, 3, 3> CompressibleNavierStokesExplicit<2, 3>::CalculateMidPointVelocityGradient() const
{
    KRATOS_TRY

    // Get geometry data
    const auto& r_geom = GetGeometry();
    const unsigned int n_nodes = r_geom.PointsNumber();
    Geometry<Node<3>>::ShapeFunctionsGradientsType dNdX_container;
    r_geom.ShapeFunctionsIntegrationPointsGradients(dNdX_container, GeometryData::GI_GAUSS_1);
    const auto& r_dNdX = dNdX_container[0];

    // Calculate midpoint magnitudes
    double midpoint_rho = 0.0;
    double midpoint_dmx_dx = 0.0;
    double midpoint_dmx_dy = 0.0;
    double midpoint_dmy_dx = 0.0;
    double midpoint_dmy_dy = 0.0;
    double midpoint_rho_dx = 0.0;
    double midpoint_rho_dy = 0.0;
    array_1d<double,3> midpoint_mom = ZeroVector(3);
    for (unsigned int i_node = 0; i_node < n_nodes; ++i_node) {
        auto& r_node = r_geom[i_node];
        const auto node_dNdX = row(r_dNdX, i_node);
        const auto& r_mom = r_node.FastGetSolutionStepValue(MOMENTUM);
        const double& r_rho = r_node.FastGetSolutionStepValue(DENSITY);
        midpoint_rho += r_rho;
        midpoint_mom += r_mom;
        midpoint_dmx_dx += r_mom[0] * node_dNdX[0];
        midpoint_dmx_dy += r_mom[0] * node_dNdX[1];
        midpoint_dmy_dx += r_mom[1] * node_dNdX[0];
        midpoint_dmy_dy += r_mom[1] * node_dNdX[1];
        midpoint_rho_dx += r_rho * node_dNdX[0];
        midpoint_rho_dy += r_rho * node_dNdX[1];
    }
    midpoint_rho /= n_nodes;
    midpoint_mom /= n_nodes;

    // Calculate velocity gradient
    // Note that the formulation is written in conservative variables. Hence we do grad(mom/rho).
    BoundedMatrix<double, 3, 3> midpoint_grad_v = ZeroMatrix(3, 3);
    midpoint_grad_v(0,0) = (midpoint_dmx_dx * midpoint_rho - midpoint_mom[0] * midpoint_rho_dx);
    midpoint_grad_v(0,1) = (midpoint_dmx_dy * midpoint_rho - midpoint_mom[0] * midpoint_rho_dy);
    midpoint_grad_v(1,0) = (midpoint_dmy_dx * midpoint_rho - midpoint_mom[1] * midpoint_rho_dx);
    midpoint_grad_v(1,1) = (midpoint_dmy_dy * midpoint_rho - midpoint_mom[1] * midpoint_rho_dy);
    midpoint_grad_v /= std::pow(midpoint_rho, 2);

    return midpoint_grad_v;

    KRATOS_CATCH("")
}

template <>
BoundedMatrix<double, 3, 3> CompressibleNavierStokesExplicit<3, 4>::CalculateMidPointVelocityGradient() const
{
    KRATOS_TRY

    // Get geometry data
    const auto& r_geom = GetGeometry();
    const unsigned int n_nodes = r_geom.PointsNumber();
    Geometry<Node<3>>::ShapeFunctionsGradientsType dNdX_container;
    r_geom.ShapeFunctionsIntegrationPointsGradients(dNdX_container, GeometryData::GI_GAUSS_1);
    const auto& r_dNdX = dNdX_container[0];

    // Calculate midpoint magnitudes
    double midpoint_rho = 0.0;
    double midpoint_dmx_dx = 0.0;
    double midpoint_dmx_dy = 0.0;
    double midpoint_dmx_dz = 0.0;
    double midpoint_dmy_dx = 0.0;
    double midpoint_dmy_dy = 0.0;
    double midpoint_dmy_dz = 0.0;
    double midpoint_dmz_dx = 0.0;
    double midpoint_dmz_dy = 0.0;
    double midpoint_dmz_dz = 0.0;
    double midpoint_rho_dx = 0.0;
    double midpoint_rho_dy = 0.0;
    double midpoint_rho_dz = 0.0;
    array_1d<double,3> midpoint_mom = ZeroVector(3);
    for (unsigned int i_node = 0; i_node < n_nodes; ++i_node) {
        auto& r_node = r_geom[i_node];
        const auto node_dNdX = row(r_dNdX, i_node);
        const auto& r_mom = r_node.FastGetSolutionStepValue(MOMENTUM);
        const double& r_rho = r_node.FastGetSolutionStepValue(DENSITY);
        midpoint_rho += r_rho;
        midpoint_mom += r_mom;
        midpoint_dmx_dx += r_mom[0] * node_dNdX[0];
        midpoint_dmx_dy += r_mom[0] * node_dNdX[1];
        midpoint_dmx_dz += r_mom[0] * node_dNdX[2];
        midpoint_dmy_dx += r_mom[1] * node_dNdX[0];
        midpoint_dmy_dy += r_mom[1] * node_dNdX[1];
        midpoint_dmy_dz += r_mom[1] * node_dNdX[2];
        midpoint_dmz_dx += r_mom[2] * node_dNdX[0];
        midpoint_dmz_dy += r_mom[2] * node_dNdX[1];
        midpoint_dmz_dz += r_mom[2] * node_dNdX[2];
        midpoint_rho_dx += r_rho * node_dNdX[0];
        midpoint_rho_dy += r_rho * node_dNdX[1];
        midpoint_rho_dz += r_rho * node_dNdX[2];
    }
    midpoint_rho /= n_nodes;
    midpoint_mom /= n_nodes;

    // Calculate velocity gradient
    // Note that the formulation is written in conservative variables. Hence we do grad(mom/rho).
    BoundedMatrix<double, 3, 3> midpoint_grad_v;
    midpoint_grad_v(0,0) = (midpoint_dmx_dx * midpoint_rho - midpoint_mom[0] * midpoint_rho_dx);
    midpoint_grad_v(0,1) = (midpoint_dmx_dy * midpoint_rho - midpoint_mom[0] * midpoint_rho_dy);
    midpoint_grad_v(0,2) = (midpoint_dmx_dz * midpoint_rho - midpoint_mom[0] * midpoint_rho_dz);
    midpoint_grad_v(1,0) = (midpoint_dmy_dx * midpoint_rho - midpoint_mom[1] * midpoint_rho_dx);
    midpoint_grad_v(1,1) = (midpoint_dmy_dy * midpoint_rho - midpoint_mom[1] * midpoint_rho_dy);
    midpoint_grad_v(1,2) = (midpoint_dmy_dz * midpoint_rho - midpoint_mom[1] * midpoint_rho_dz);
    midpoint_grad_v(2,0) = (midpoint_dmz_dx * midpoint_rho - midpoint_mom[2] * midpoint_rho_dx);
    midpoint_grad_v(2,1) = (midpoint_dmz_dy * midpoint_rho - midpoint_mom[2] * midpoint_rho_dy);
    midpoint_grad_v(2,2) = (midpoint_dmz_dz * midpoint_rho - midpoint_mom[2] * midpoint_rho_dz);
    midpoint_grad_v /= std::pow(midpoint_rho, 2);

    return midpoint_grad_v;

    KRATOS_CATCH("")
}

template <>
void CompressibleNavierStokesExplicit<2,3>::CalculateMomentumProjection(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    constexpr unsigned int dim = 2;
    constexpr unsigned int n_nodes = 3;

    // Struct to pass around the data
    ElementDataStruct data;
    this->FillElementData(data, rCurrentProcessInfo);

    // Substitute the formulation symbols by the data structure values
    const BoundedMatrix<double, n_nodes, 2> &f_ext = data.f_ext;
    const double gamma = data.gamma;

    // Solution vector values and time derivatives from nodal data
    // This is intentionally done in this way to limit the matrix acceses
    // The notation U_i_j DOF j value in node i
    const double &U_0_0 = data.U(0, 0);
    const double &U_0_1 = data.U(0, 1);
    const double &U_0_2 = data.U(0, 2);
    const double &U_0_3 = data.U(0, 3);
    const double &U_1_0 = data.U(1, 0);
    const double &U_1_1 = data.U(1, 1);
    const double &U_1_2 = data.U(1, 2);
    const double &U_1_3 = data.U(1, 3);
    const double &U_2_0 = data.U(2, 0);
    const double &U_2_1 = data.U(2, 1);
    const double &U_2_2 = data.U(2, 2);
    const double &U_2_3 = data.U(2, 3);

    const double &dUdt_0_1 = data.dUdt(0, 1);
    const double &dUdt_0_2 = data.dUdt(0, 2);
    const double &dUdt_1_1 = data.dUdt(1, 1);
    const double &dUdt_1_2 = data.dUdt(1, 2);
    const double &dUdt_2_1 = data.dUdt(2, 1);
    const double &dUdt_2_2 = data.dUdt(2, 2);

    // Hardcoded shape functions gradients for linear triangular element
    // This is explicitly done to minimize the matrix acceses
    // The notation DN_i_j means shape function for node i in dimension j
    const double &DN_DX_0_0 = data.DN_DX(0, 0);
    const double &DN_DX_0_1 = data.DN_DX(0, 1);
    const double &DN_DX_1_0 = data.DN_DX(1, 0);
    const double &DN_DX_1_1 = data.DN_DX(1, 1);
    const double &DN_DX_2_0 = data.DN_DX(2, 0);
    const double &DN_DX_2_1 = data.DN_DX(2, 1);

    // Calculate shock capturing values
    BoundedVector<double, 6> mom_proj;

    //substitute_mom_proj_2D
    // Here we assume that all the weights of the gauss points are the same so we multiply at the end by Volume/n_nodes
    mom_proj *= data.volume / static_cast<double>(n_nodes);

    // Assembly the projection contributions
    auto& r_geometry = GetGeometry();
    for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
        const IndexType aux = i_node * dim;
        auto& r_mom_proj = r_geometry[i_node].GetValue(MOMENTUM_PROJECTION);
        for (IndexType d = 0; d < dim; ++d) {
#pragma omp atomic
            r_mom_proj[d] += mom_proj[aux + d];
        }
    }

    KRATOS_CATCH("")
}

template <>
void CompressibleNavierStokesExplicit<3,4>::CalculateMomentumProjection(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    constexpr unsigned int dim = 3;
    constexpr unsigned int n_nodes = 4;

    // Struct to pass around the data
    ElementDataStruct data;
    this->FillElementData(data, rCurrentProcessInfo);

    // Substitute the formulation symbols by the data structure values
    const BoundedMatrix<double, n_nodes, 3> &f_ext = data.f_ext;
    const double gamma = data.gamma;

    // Solution vector values and time derivatives from nodal data
    // This is intentionally done in this way to limit the matrix acceses
    // The notation U_i_j DOF j value in node i
    const double &U_0_0 = data.U(0, 0);
    const double &U_0_1 = data.U(0, 1);
    const double &U_0_2 = data.U(0, 2);
    const double &U_0_3 = data.U(0, 3);
    const double &U_0_4 = data.U(0, 4);
    const double &U_1_0 = data.U(1, 0);
    const double &U_1_1 = data.U(1, 1);
    const double &U_1_2 = data.U(1, 2);
    const double &U_1_3 = data.U(1, 3);
    const double &U_1_4 = data.U(1, 4);
    const double &U_2_0 = data.U(2, 0);
    const double &U_2_1 = data.U(2, 1);
    const double &U_2_2 = data.U(2, 2);
    const double &U_2_3 = data.U(2, 3);
    const double &U_2_4 = data.U(2, 4);
    const double &U_3_0 = data.U(3, 0);
    const double &U_3_1 = data.U(3, 1);
    const double &U_3_2 = data.U(3, 2);
    const double &U_3_3 = data.U(3, 3);
    const double &U_3_4 = data.U(3, 4);

    const double &dUdt_0_1 = data.dUdt(0, 1);
    const double &dUdt_0_2 = data.dUdt(0, 2);
    const double &dUdt_0_3 = data.dUdt(0, 3);
    const double &dUdt_1_1 = data.dUdt(1, 1);
    const double &dUdt_1_2 = data.dUdt(1, 2);
    const double &dUdt_1_3 = data.dUdt(1, 3);
    const double &dUdt_2_1 = data.dUdt(2, 1);
    const double &dUdt_2_2 = data.dUdt(2, 2);
    const double &dUdt_2_3 = data.dUdt(2, 3);
    const double &dUdt_3_1 = data.dUdt(3, 1);
    const double &dUdt_3_2 = data.dUdt(3, 2);
    const double &dUdt_3_3 = data.dUdt(3, 3);

    // Hardcoded shape functions gradients for linear tetrahedra element
    // This is explicitly done to minimize the matrix acceses
    // The notation DN_i_j means shape function for node i in dimension j
    const double &DN_DX_0_0 = data.DN_DX(0, 0);
    const double &DN_DX_0_1 = data.DN_DX(0, 1);
    const double &DN_DX_0_2 = data.DN_DX(0, 2);
    const double &DN_DX_1_0 = data.DN_DX(1, 0);
    const double &DN_DX_1_1 = data.DN_DX(1, 1);
    const double &DN_DX_1_2 = data.DN_DX(1, 2);
    const double &DN_DX_2_0 = data.DN_DX(2, 0);
    const double &DN_DX_2_1 = data.DN_DX(2, 1);
    const double &DN_DX_2_2 = data.DN_DX(2, 2);
    const double &DN_DX_3_0 = data.DN_DX(3, 0);
    const double &DN_DX_3_1 = data.DN_DX(3, 1);
    const double &DN_DX_3_2 = data.DN_DX(3, 2);

    // Calculate shock capturing values
    BoundedVector<double, 12> mom_proj;

    //substitute_mom_proj_3D
    // Here we assume that all the weights of the gauss points are the same so we multiply at the end by Volume/n_nodes
    mom_proj *= data.volume / static_cast<double>(n_nodes);

    // Assembly the projection contributions
    auto& r_geometry = GetGeometry();
    for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
        const IndexType aux = i_node * dim;
        auto& r_mom_proj = r_geometry[i_node].GetValue(MOMENTUM_PROJECTION);
        for (IndexType d = 0; d < dim; ++d) {
#pragma omp atomic
            r_mom_proj[d] += mom_proj[aux + d];
        }
    }

    KRATOS_CATCH("")
}

template <>
void CompressibleNavierStokesExplicit<2,3>::CalculateDensityProjection(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    constexpr unsigned int n_nodes = 3;

    // Struct to pass around the data
    ElementDataStruct data;
    this->FillElementData(data, rCurrentProcessInfo);

    // Substitute the formulation symbols by the data structure values
    const array_1d<double, n_nodes> &m_ext = data.m_ext;

    // Solution vector values and time derivatives from nodal data
    // This is intentionally done in this way to limit the matrix acceses
    // The notation U_i_j DOF j value in node i
    const double &U_0_1 = data.U(0, 1);
    const double &U_0_2 = data.U(0, 2);
    const double &U_1_1 = data.U(1, 1);
    const double &U_1_2 = data.U(1, 2);
    const double &U_2_1 = data.U(2, 1);
    const double &U_2_2 = data.U(2, 2);

    const double &dUdt_0_0 = data.dUdt(0, 0);
    const double &dUdt_1_0 = data.dUdt(1, 0);
    const double &dUdt_2_0 = data.dUdt(2, 0);

    // Hardcoded shape functions gradients for linear triangular element
    // This is explicitly done to minimize the matrix acceses
    // The notation DN_i_j means shape function for node i in dimension j
    const double &DN_DX_0_0 = data.DN_DX(0, 0);
    const double &DN_DX_0_1 = data.DN_DX(0, 1);
    const double &DN_DX_1_0 = data.DN_DX(1, 0);
    const double &DN_DX_1_1 = data.DN_DX(1, 1);
    const double &DN_DX_2_0 = data.DN_DX(2, 0);
    const double &DN_DX_2_1 = data.DN_DX(2, 1);

    // Calculate shock capturing values
    BoundedVector<double, 3> rho_proj;

    //substitute_rho_proj_2D
    // Here we assume that all the weights of the gauss points are the same so we multiply at the end by Volume/n_nodes
    rho_proj *= data.volume / static_cast<double>(n_nodes);

    // Assembly the projection contributions
    auto& r_geometry = GetGeometry();
    for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
#pragma omp atomic
        r_geometry[i_node].GetValue(DENSITY_PROJECTION) += rho_proj[i_node];
    }

    KRATOS_CATCH("")
}

template <>
void CompressibleNavierStokesExplicit<3,4>::CalculateDensityProjection(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    constexpr unsigned int n_nodes = 4;

    // Struct to pass around the data
    ElementDataStruct data;
    this->FillElementData(data, rCurrentProcessInfo);

    // Substitute the formulation symbols by the data structure values
    const array_1d<double, n_nodes> &m_ext = data.m_ext;

    // Solution vector values and time derivatives from nodal data
    // This is intentionally done in this way to limit the matrix acceses
    // The notation U_i_j DOF j value in node i
    const double &U_0_1 = data.U(0, 1);
    const double &U_0_2 = data.U(0, 2);
    const double &U_0_3 = data.U(0, 3);
    const double &U_1_1 = data.U(1, 1);
    const double &U_1_2 = data.U(1, 2);
    const double &U_1_3 = data.U(1, 3);
    const double &U_2_1 = data.U(2, 1);
    const double &U_2_2 = data.U(2, 2);
    const double &U_2_3 = data.U(2, 3);
    const double &U_3_1 = data.U(3, 1);
    const double &U_3_2 = data.U(3, 2);
    const double &U_3_3 = data.U(3, 3);

    const double &dUdt_0_0 = data.dUdt(0, 0);
    const double &dUdt_1_0 = data.dUdt(1, 0);
    const double &dUdt_2_0 = data.dUdt(2, 0);
    const double &dUdt_3_0 = data.dUdt(3, 0);

    // Hardcoded shape functions gradients for linear tetrahedra element
    // This is explicitly done to minimize the matrix acceses
    // The notation DN_i_j means shape function for node i in dimension j
    const double &DN_DX_0_0 = data.DN_DX(0, 0);
    const double &DN_DX_0_1 = data.DN_DX(0, 1);
    const double &DN_DX_0_2 = data.DN_DX(0, 2);
    const double &DN_DX_1_0 = data.DN_DX(1, 0);
    const double &DN_DX_1_1 = data.DN_DX(1, 1);
    const double &DN_DX_1_2 = data.DN_DX(1, 2);
    const double &DN_DX_2_0 = data.DN_DX(2, 0);
    const double &DN_DX_2_1 = data.DN_DX(2, 1);
    const double &DN_DX_2_2 = data.DN_DX(2, 2);
    const double &DN_DX_3_0 = data.DN_DX(3, 0);
    const double &DN_DX_3_1 = data.DN_DX(3, 1);
    const double &DN_DX_3_2 = data.DN_DX(3, 2);

    // Calculate shock capturing values
    BoundedVector<double, 4> rho_proj;

    //substitute_rho_proj_3D
    // Here we assume that all the weights of the gauss points are the same so we multiply at the end by Volume/n_nodes
    rho_proj *= data.volume / static_cast<double>(n_nodes);

    // Assembly the projection contributions
    auto& r_geometry = GetGeometry();
    for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
#pragma omp atomic
        r_geometry[i_node].GetValue(DENSITY_PROJECTION) += rho_proj[i_node];
    }

    KRATOS_CATCH("")
}

template <>
void CompressibleNavierStokesExplicit<2,3>::CalculateTotalEnergyProjection(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    constexpr unsigned int n_nodes = 3;

    // Struct to pass around the data
    ElementDataStruct data;
    this->FillElementData(data, rCurrentProcessInfo);

    // Substitute the formulation symbols by the data structure values
    const array_1d<double, n_nodes> &r_ext = data.r_ext;
    const BoundedMatrix<double, n_nodes, 2> &f_ext = data.f_ext;
    const double gamma = data.gamma;

    // Solution vector values and time derivatives from nodal data
    // This is intentionally done in this way to limit the matrix acceses
    // The notation U_i_j DOF j value in node i
    const double &U_0_0 = data.U(0, 0);
    const double &U_0_1 = data.U(0, 1);
    const double &U_0_2 = data.U(0, 2);
    const double &U_0_3 = data.U(0, 3);
    const double &U_1_0 = data.U(1, 0);
    const double &U_1_1 = data.U(1, 1);
    const double &U_1_2 = data.U(1, 2);
    const double &U_1_3 = data.U(1, 3);
    const double &U_2_0 = data.U(2, 0);
    const double &U_2_1 = data.U(2, 1);
    const double &U_2_2 = data.U(2, 2);
    const double &U_2_3 = data.U(2, 3);

    const double &dUdt_0_3 = data.dUdt(0, 3);
    const double &dUdt_1_3 = data.dUdt(1, 3);
    const double &dUdt_2_3 = data.dUdt(2, 3);

    // Hardcoded shape functions gradients for linear triangular element
    // This is explicitly done to minimize the matrix acceses
    // The notation DN_i_j means shape function for node i in dimension j
    const double &DN_DX_0_0 = data.DN_DX(0, 0);
    const double &DN_DX_0_1 = data.DN_DX(0, 1);
    const double &DN_DX_1_0 = data.DN_DX(1, 0);
    const double &DN_DX_1_1 = data.DN_DX(1, 1);
    const double &DN_DX_2_0 = data.DN_DX(2, 0);
    const double &DN_DX_2_1 = data.DN_DX(2, 1);

    // Calculate shock capturing values
    BoundedVector<double, 3> tot_ener_proj;

    //substitute_tot_ener_proj_2D
    // Here we assume that all the weights of the gauss points are the same so we multiply at the end by Volume/n_nodes
    tot_ener_proj *= data.volume / static_cast<double>(n_nodes);

    // Assembly the projection contributions
    auto& r_geometry = GetGeometry();
    for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
#pragma omp atomic
        r_geometry[i_node].GetValue(TOTAL_ENERGY_PROJECTION) += tot_ener_proj[i_node];
    }

    KRATOS_CATCH("")
}

template <>
void CompressibleNavierStokesExplicit<3,4>::CalculateTotalEnergyProjection(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    constexpr unsigned int n_nodes = 4;

    // Struct to pass around the data
    ElementDataStruct data;
    this->FillElementData(data, rCurrentProcessInfo);

    // Substitute the formulation symbols by the data structure values
    const array_1d<double, n_nodes> &r_ext = data.r_ext;
    const BoundedMatrix<double, n_nodes, 3> &f_ext = data.f_ext;
    const double gamma = data.gamma;

    // Solution vector values and time derivatives from nodal data
    // This is intentionally done in this way to limit the matrix acceses
    // The notation U_i_j DOF j value in node i
    const double &U_0_0 = data.U(0, 0);
    const double &U_0_1 = data.U(0, 1);
    const double &U_0_2 = data.U(0, 2);
    const double &U_0_3 = data.U(0, 3);
    const double &U_0_4 = data.U(0, 4);
    const double &U_1_0 = data.U(1, 0);
    const double &U_1_1 = data.U(1, 1);
    const double &U_1_2 = data.U(1, 2);
    const double &U_1_3 = data.U(1, 3);
    const double &U_1_4 = data.U(1, 4);
    const double &U_2_0 = data.U(2, 0);
    const double &U_2_1 = data.U(2, 1);
    const double &U_2_2 = data.U(2, 2);
    const double &U_2_3 = data.U(2, 3);
    const double &U_2_4 = data.U(2, 4);
    const double &U_3_0 = data.U(3, 0);
    const double &U_3_1 = data.U(3, 1);
    const double &U_3_2 = data.U(3, 2);
    const double &U_3_3 = data.U(3, 3);
    const double &U_3_4 = data.U(3, 4);

    const double &dUdt_0_4 = data.dUdt(0, 4);
    const double &dUdt_1_4 = data.dUdt(1, 4);
    const double &dUdt_2_4 = data.dUdt(2, 4);
    const double &dUdt_3_4 = data.dUdt(3, 4);

    // Hardcoded shape functions gradients for linear tetrahedra element
    // This is explicitly done to minimize the matrix acceses
    // The notation DN_i_j means shape function for node i in dimension j
    const double &DN_DX_0_0 = data.DN_DX(0, 0);
    const double &DN_DX_0_1 = data.DN_DX(0, 1);
    const double &DN_DX_0_2 = data.DN_DX(0, 2);
    const double &DN_DX_1_0 = data.DN_DX(1, 0);
    const double &DN_DX_1_1 = data.DN_DX(1, 1);
    const double &DN_DX_1_2 = data.DN_DX(1, 2);
    const double &DN_DX_2_0 = data.DN_DX(2, 0);
    const double &DN_DX_2_1 = data.DN_DX(2, 1);
    const double &DN_DX_2_2 = data.DN_DX(2, 2);
    const double &DN_DX_3_0 = data.DN_DX(3, 0);
    const double &DN_DX_3_1 = data.DN_DX(3, 1);
    const double &DN_DX_3_2 = data.DN_DX(3, 2);

    // Calculate shock capturing values
    BoundedVector<double, 4> tot_ener_proj;

    //substitute_tot_ener_proj_3D
    // Here we assume that all the weights of the gauss points are the same so we multiply at the end by Volume/n_nodes
    tot_ener_proj *= data.volume / static_cast<double>(n_nodes);

    // Assembly the projection contributions
    auto& r_geometry = GetGeometry();
    for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
#pragma omp atomic
        r_geometry[i_node].GetValue(TOTAL_ENERGY_PROJECTION) += tot_ener_proj[i_node];
    }

    KRATOS_CATCH("")
}

template <>
void CompressibleNavierStokesExplicit<2,3>::CalculateRightHandSideInternal(
    BoundedVector<double, 12> &rRightHandSideBoundedVector,
    const ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY

    constexpr unsigned int n_nodes = 3;

    // Struct to pass around the data
    ElementDataStruct data;
    this->FillElementData(data, rCurrentProcessInfo);

    // Substitute the formulation symbols by the data structure values
    const double h = data.h;
    const array_1d<double, n_nodes>& r_ext = data.r_ext;
    const array_1d<double, n_nodes>& m_ext = data.m_ext;
    const array_1d<double, n_nodes>& mu_sc_nodes = data.mu_sc_nodes;
    const array_1d<double, n_nodes>& beta_sc_nodes = data.beta_sc_nodes;
    const array_1d<double, n_nodes>& lamb_sc_nodes = data.lamb_sc_nodes;
    const BoundedMatrix<double, n_nodes, 2>& f_ext = data.f_ext;
    const double mu = data.mu;
    const double c_v = data.c_v;
    const double gamma = data.gamma;
    const double lambda = data.lambda;

    // Stabilization parameters
    const double stab_c1 = 12.0;
    const double stab_c2 = 2.0;
    const double stab_c3 = 1.0;

    // Solution vector values and time derivatives from nodal data
    // This is intentionally done in this way to limit the matrix acceses
    // The notation U_i_j DOF j value in node i
    const double& U_0_0 = data.U(0, 0);
    const double& U_0_1 = data.U(0, 1);
    const double& U_0_2 = data.U(0, 2);
    const double& U_0_3 = data.U(0, 3);
    const double& U_1_0 = data.U(1, 0);
    const double& U_1_1 = data.U(1, 1);
    const double& U_1_2 = data.U(1, 2);
    const double& U_1_3 = data.U(1, 3);
    const double& U_2_0 = data.U(2, 0);
    const double& U_2_1 = data.U(2, 1);
    const double& U_2_2 = data.U(2, 2);
    const double& U_2_3 = data.U(2, 3);

    const double& dUdt_0_0 = data.dUdt(0, 0);
    const double& dUdt_0_1 = data.dUdt(0, 1);
    const double& dUdt_0_2 = data.dUdt(0, 2);
    const double& dUdt_0_3 = data.dUdt(0, 3);
    const double& dUdt_1_0 = data.dUdt(1, 0);
    const double& dUdt_1_1 = data.dUdt(1, 1);
    const double& dUdt_1_2 = data.dUdt(1, 2);
    const double& dUdt_1_3 = data.dUdt(1, 3);
    const double& dUdt_2_0 = data.dUdt(2, 0);
    const double& dUdt_2_1 = data.dUdt(2, 1);
    const double& dUdt_2_2 = data.dUdt(2, 2);
    const double& dUdt_2_3 = data.dUdt(2, 3);

    // Hardcoded shape functions gradients for linear triangular element
    // This is explicitly done to minimize the matrix acceses
    // The notation DN_i_j means shape function for node i in dimension j
    const double& DN_DX_0_0 = data.DN_DX(0, 0);
    const double& DN_DX_0_1 = data.DN_DX(0, 1);
    const double& DN_DX_1_0 = data.DN_DX(1, 0);
    const double& DN_DX_1_1 = data.DN_DX(1, 1);
    const double& DN_DX_2_0 = data.DN_DX(2, 0);
    const double& DN_DX_2_1 = data.DN_DX(2, 1);

    if (data.UseOSS) {
        // Projections container accesses
        const double& ResProj_0_0 = data.ResProj(0, 0);
        const double& ResProj_0_1 = data.ResProj(0, 1);
        const double& ResProj_0_2 = data.ResProj(0, 2);
        const double& ResProj_0_3 = data.ResProj(0, 3);
        const double& ResProj_1_0 = data.ResProj(1, 0);
        const double& ResProj_1_1 = data.ResProj(1, 1);
        const double& ResProj_1_2 = data.ResProj(1, 2);
        const double& ResProj_1_3 = data.ResProj(1, 3);
        const double& ResProj_2_0 = data.ResProj(2, 0);
        const double& ResProj_2_1 = data.ResProj(2, 1);
        const double& ResProj_2_2 = data.ResProj(2, 2);
        const double& ResProj_2_3 = data.ResProj(2, 3);

        //substitute_rhs_2D_OSS
    } else {
        //substitute_rhs_2D_ASGS
    }

    // Here we assume that all the weights of the gauss points are the same so we multiply at the end by Volume/n_nodes
    rRightHandSideBoundedVector *= data.volume / static_cast<double>(n_nodes);

    KRATOS_CATCH("")
}

template<>
void CompressibleNavierStokesExplicit<3,4>::CalculateRightHandSideInternal(
    BoundedVector<double, 20> &rRightHandSideBoundedVector,
    const ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY

    constexpr unsigned int n_nodes = 4;

    // Struct to pass around the data
    ElementDataStruct data;
    this->FillElementData(data, rCurrentProcessInfo);

    // Substitute the formulation symbols by the data structure values
    const double h = data.h;
    const array_1d<double, n_nodes>& r_ext = data.r_ext;
    const array_1d<double, n_nodes>& m_ext = data.m_ext;
    const array_1d<double, n_nodes>& mu_sc_nodes = data.mu_sc_nodes;
    const array_1d<double, n_nodes>& beta_sc_nodes = data.beta_sc_nodes;
    const array_1d<double, n_nodes>& lamb_sc_nodes = data.lamb_sc_nodes;
    const BoundedMatrix<double, n_nodes, 3>& f_ext = data.f_ext;
    const double mu = data.mu;
    const double c_v = data.c_v;
    const double gamma = data.gamma;
    const double lambda = data.lambda;

    // Stabilization parameters
    const double stab_c1 = 12.0;
    const double stab_c2 = 2.0;
    const double stab_c3 = 1.0;

    // Solution vector values and time derivatives from nodal data
    // This is intentionally done in this way to limit the matrix acceses
    // The notation U_i_j DOF j value in node i
    const double& U_0_0 = data.U(0, 0);
    const double& U_0_1 = data.U(0, 1);
    const double& U_0_2 = data.U(0, 2);
    const double& U_0_3 = data.U(0, 3);
    const double& U_0_4 = data.U(0, 4);
    const double& U_1_0 = data.U(1, 0);
    const double& U_1_1 = data.U(1, 1);
    const double& U_1_2 = data.U(1, 2);
    const double& U_1_3 = data.U(1, 3);
    const double& U_1_4 = data.U(1, 4);
    const double& U_2_0 = data.U(2, 0);
    const double& U_2_1 = data.U(2, 1);
    const double& U_2_2 = data.U(2, 2);
    const double& U_2_3 = data.U(2, 3);
    const double& U_2_4 = data.U(2, 4);
    const double& U_3_0 = data.U(3, 0);
    const double& U_3_1 = data.U(3, 1);
    const double& U_3_2 = data.U(3, 2);
    const double& U_3_3 = data.U(3, 3);
    const double& U_3_4 = data.U(3, 4);

    const double& dUdt_0_0 = data.dUdt(0, 0);
    const double& dUdt_0_1 = data.dUdt(0, 1);
    const double& dUdt_0_2 = data.dUdt(0, 2);
    const double& dUdt_0_3 = data.dUdt(0, 3);
    const double& dUdt_0_4 = data.dUdt(0, 4);
    const double& dUdt_1_0 = data.dUdt(1, 0);
    const double& dUdt_1_1 = data.dUdt(1, 1);
    const double& dUdt_1_2 = data.dUdt(1, 2);
    const double& dUdt_1_3 = data.dUdt(1, 3);
    const double& dUdt_1_4 = data.dUdt(1, 4);
    const double& dUdt_2_0 = data.dUdt(2, 0);
    const double& dUdt_2_1 = data.dUdt(2, 1);
    const double& dUdt_2_2 = data.dUdt(2, 2);
    const double& dUdt_2_3 = data.dUdt(2, 3);
    const double& dUdt_2_4 = data.dUdt(2, 4);
    const double& dUdt_3_0 = data.dUdt(3, 0);
    const double& dUdt_3_1 = data.dUdt(3, 1);
    const double& dUdt_3_2 = data.dUdt(3, 2);
    const double& dUdt_3_3 = data.dUdt(3, 3);
    const double& dUdt_3_4 = data.dUdt(3, 4);

    // Hardcoded shape functions gradients for linear tetrahedra element
    // This is explicitly done to minimize the matrix acceses
    // The notation DN_i_j means shape function for node i in dimension j
    const double& DN_DX_0_0 = data.DN_DX(0, 0);
    const double& DN_DX_0_1 = data.DN_DX(0, 1);
    const double& DN_DX_0_2 = data.DN_DX(0, 2);
    const double& DN_DX_1_0 = data.DN_DX(1, 0);
    const double& DN_DX_1_1 = data.DN_DX(1, 1);
    const double& DN_DX_1_2 = data.DN_DX(1, 2);
    const double& DN_DX_2_0 = data.DN_DX(2, 0);
    const double& DN_DX_2_1 = data.DN_DX(2, 1);
    const double& DN_DX_2_2 = data.DN_DX(2, 2);
    const double& DN_DX_3_0 = data.DN_DX(3, 0);
    const double& DN_DX_3_1 = data.DN_DX(3, 1);
    const double& DN_DX_3_2 = data.DN_DX(3, 2);

    if (data.UseOSS) {
        // Projections container accesses
        const double& ResProj_0_0 = data.ResProj(0, 0);
        const double& ResProj_0_1 = data.ResProj(0, 1);
        const double& ResProj_0_2 = data.ResProj(0, 2);
        const double& ResProj_0_3 = data.ResProj(0, 3);
        const double& ResProj_0_4 = data.ResProj(0, 4);
        const double& ResProj_1_0 = data.ResProj(1, 0);
        const double& ResProj_1_1 = data.ResProj(1, 1);
        const double& ResProj_1_2 = data.ResProj(1, 2);
        const double& ResProj_1_3 = data.ResProj(1, 3);
        const double& ResProj_1_4 = data.ResProj(1, 4);
        const double& ResProj_2_0 = data.ResProj(2, 0);
        const double& ResProj_2_1 = data.ResProj(2, 1);
        const double& ResProj_2_2 = data.ResProj(2, 2);
        const double& ResProj_2_3 = data.ResProj(2, 3);
        const double& ResProj_2_4 = data.ResProj(2, 4);
        const double& ResProj_3_0 = data.ResProj(3, 0);
        const double& ResProj_3_1 = data.ResProj(3, 1);
        const double& ResProj_3_2 = data.ResProj(3, 2);
        const double& ResProj_3_3 = data.ResProj(3, 3);
        const double& ResProj_3_4 = data.ResProj(3, 4);

        //substitute_rhs_3D_OSS
    } else {
        //substitute_rhs_3D_ASGS
    }

    // Here we assume that all the weights of the gauss points are the same so we multiply at the end by Volume/n_nodes
    rRightHandSideBoundedVector *= data.volume / static_cast<double>(n_nodes);

    KRATOS_CATCH("")
}

template <>
void CompressibleNavierStokesExplicit<2,3>::AddExplicitContribution(const ProcessInfo &rCurrentProcessInfo)
{
    constexpr IndexType dim = 2;
    constexpr IndexType n_nodes = 3;
    constexpr IndexType block_size = 4;

    // Calculate the explicit residual vector
    BoundedVector<double, 12> rhs;
    CalculateRightHandSideInternal(rhs, rCurrentProcessInfo);

    // Add the residual contribution
    // Note that the reaction is indeed the formulation residual
    auto& r_geometry = GetGeometry();
    for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
        const IndexType aux = i_node * block_size;
#pragma omp atomic
        r_geometry[i_node].FastGetSolutionStepValue(REACTION_DENSITY) += rhs[aux];
        auto& r_mom = r_geometry[i_node].FastGetSolutionStepValue(REACTION);
        for (IndexType d = 0; d < dim; ++d) {
#pragma omp atomic
            r_mom[d] += rhs[aux + (d + 1)];
        }
#pragma omp atomic
        r_geometry[i_node].FastGetSolutionStepValue(REACTION_ENERGY) += rhs[aux + 3];
    }
}

template <>
void CompressibleNavierStokesExplicit<3,4>::AddExplicitContribution(const ProcessInfo &rCurrentProcessInfo)
{
    constexpr IndexType dim = 3;
    constexpr IndexType n_nodes = 4;
    constexpr IndexType block_size = 5;

    // Calculate the explicit residual vector
    BoundedVector<double, 20> rhs;
    CalculateRightHandSideInternal(rhs, rCurrentProcessInfo);

    // Add the residual contribution
    // Note that the reaction is indeed the formulation residual
    auto& r_geometry = GetGeometry();
    for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
        const IndexType aux = i_node * block_size;
#pragma omp atomic
        r_geometry[i_node].FastGetSolutionStepValue(REACTION_DENSITY) += rhs[aux];
        auto& r_mom = r_geometry[i_node].FastGetSolutionStepValue(REACTION);
        for (IndexType d = 0; d < dim; ++d) {
#pragma omp atomic
            r_mom[d] += rhs[aux + (d + 1)];
        }
#pragma omp atomic
        r_geometry[i_node].FastGetSolutionStepValue(REACTION_ENERGY) += rhs[aux + 4];
    }
}

template <>
void CompressibleNavierStokesExplicit<2,3>::CalculateMassMatrix(
    MatrixType &rMassMatrix,
    const ProcessInfo &rCurrentProcessInfo)
{
    constexpr IndexType n_nodes = 3;
    constexpr IndexType block_size = 4;

    // Initialize and fill the mass matrix values
    const double one_six = 1.0 / 6.0;
    const double one_twelve = 1.0 / 12.0;
    const unsigned int size = n_nodes * block_size;
    rMassMatrix = ZeroMatrix(size, size);
    rMassMatrix(0, 0) = one_six; rMassMatrix(0, 4) = one_twelve; rMassMatrix(0, 8) = one_twelve;
    rMassMatrix(1, 1) = one_six; rMassMatrix(1, 5) = one_twelve; rMassMatrix(1, 9) = one_twelve;
    rMassMatrix(2, 2) = one_six; rMassMatrix(2, 6) = one_twelve; rMassMatrix(2, 10) = one_twelve;
    rMassMatrix(3, 3) = one_six; rMassMatrix(3, 7) = one_twelve; rMassMatrix(3, 11) = one_twelve;
    rMassMatrix(4, 0) = one_twelve; rMassMatrix(4, 4) = one_six; rMassMatrix(4, 8) = one_twelve;
    rMassMatrix(5, 1) = one_twelve; rMassMatrix(5, 5) = one_six; rMassMatrix(5, 9) = one_twelve;
    rMassMatrix(6, 2) = one_twelve; rMassMatrix(6, 6) = one_six; rMassMatrix(6, 10) = one_twelve;
    rMassMatrix(7, 3) = one_twelve; rMassMatrix(7, 7) = one_six; rMassMatrix(7, 11) = one_twelve;
    rMassMatrix(8, 0) = one_twelve; rMassMatrix(8, 4) = one_twelve; rMassMatrix(8, 8) = one_six;
    rMassMatrix(9, 1) = one_twelve; rMassMatrix(9, 5) = one_twelve; rMassMatrix(9, 9) = one_six;
    rMassMatrix(10, 2) = one_twelve; rMassMatrix(10, 6) = one_twelve; rMassMatrix(10, 10) = one_six;
    rMassMatrix(11, 3) = one_twelve; rMassMatrix(11, 7) = one_twelve; rMassMatrix(11, 11) = one_six;

    // Here we assume that all the Gauss pt. have the same weight so we multiply by the volume
    rMassMatrix *= GetGeometry().Area();
}

template <>
void CompressibleNavierStokesExplicit<3,4>::CalculateMassMatrix(
    MatrixType &rMassMatrix,
    const ProcessInfo &rCurrentProcessInfo)
{
    constexpr IndexType n_nodes = 4;
    constexpr IndexType block_size = 5;

    // Initialize and fill the mass matrix values
    const double one_ten = 0.1;
    const double one_twenty = 0.05;
    const unsigned int size = n_nodes * block_size;
    rMassMatrix = ZeroMatrix(size, size);
    rMassMatrix(0, 0) = one_ten; rMassMatrix(0, 5) = one_twenty; rMassMatrix(0, 10) = one_twenty; rMassMatrix(0,15) = one_twenty;
    rMassMatrix(1, 1) = one_ten; rMassMatrix(1, 6) = one_twenty; rMassMatrix(1, 11) = one_twenty; rMassMatrix(1,16) = one_twenty;
    rMassMatrix(2, 2) = one_ten; rMassMatrix(2, 7) = one_twenty; rMassMatrix(2, 12) = one_twenty; rMassMatrix(2,17) = one_twenty;
    rMassMatrix(3, 3) = one_ten; rMassMatrix(3, 8) = one_twenty; rMassMatrix(3, 13) = one_twenty; rMassMatrix(3,18) = one_twenty;
    rMassMatrix(4, 4) = one_ten; rMassMatrix(4, 9) = one_twenty; rMassMatrix(4, 14) = one_twenty; rMassMatrix(4,19) = one_twenty;
    rMassMatrix(5, 0) = one_twenty; rMassMatrix(5, 5) = one_ten; rMassMatrix(5, 10) = one_twenty; rMassMatrix(5,15) = one_twenty;
    rMassMatrix(6, 1) = one_twenty; rMassMatrix(6, 6) = one_ten; rMassMatrix(6, 11) = one_twenty; rMassMatrix(6,16) = one_twenty;
    rMassMatrix(7, 2) = one_twenty; rMassMatrix(7, 7) = one_ten; rMassMatrix(7, 12) = one_twenty; rMassMatrix(7,17) = one_twenty;
    rMassMatrix(8, 3) = one_twenty; rMassMatrix(8, 8) = one_ten; rMassMatrix(8, 13) = one_twenty; rMassMatrix(8,18) = one_twenty;
    rMassMatrix(9, 4) = one_twenty; rMassMatrix(9, 9) = one_ten; rMassMatrix(9, 14) = one_twenty; rMassMatrix(9,19) = one_twenty;
    rMassMatrix(10, 0) = one_twenty; rMassMatrix(10, 5) = one_twenty; rMassMatrix(10, 10) = one_ten; rMassMatrix(10,15) = one_twenty;
    rMassMatrix(11, 1) = one_twenty; rMassMatrix(11, 6) = one_twenty; rMassMatrix(11, 11) = one_ten; rMassMatrix(11,16) = one_twenty;
    rMassMatrix(12, 2) = one_twenty; rMassMatrix(12, 7) = one_twenty; rMassMatrix(12, 12) = one_ten; rMassMatrix(12,17) = one_twenty;
    rMassMatrix(13, 3) = one_twenty; rMassMatrix(13, 8) = one_twenty; rMassMatrix(13, 13) = one_ten; rMassMatrix(13,18) = one_twenty;
    rMassMatrix(14, 4) = one_twenty; rMassMatrix(14, 9) = one_twenty; rMassMatrix(14, 14) = one_ten; rMassMatrix(14,19) = one_twenty;
    rMassMatrix(15, 0) = one_twenty; rMassMatrix(15, 5) = one_twenty; rMassMatrix(15, 10) = one_twenty; rMassMatrix(15,15) = one_ten;
    rMassMatrix(16, 1) = one_twenty; rMassMatrix(16, 6) = one_twenty; rMassMatrix(16, 11) = one_twenty; rMassMatrix(16,16) = one_ten;
    rMassMatrix(17, 2) = one_twenty; rMassMatrix(17, 7) = one_twenty; rMassMatrix(17, 12) = one_twenty; rMassMatrix(17,17) = one_ten;
    rMassMatrix(18, 3) = one_twenty; rMassMatrix(18, 8) = one_twenty; rMassMatrix(18, 13) = one_twenty; rMassMatrix(18,18) = one_ten;
    rMassMatrix(19, 4) = one_twenty; rMassMatrix(19, 9) = one_twenty; rMassMatrix(19, 14) = one_twenty; rMassMatrix(19,19) = one_ten;

    // Here we assume that all the Gauss pt. have the same weight so we multiply by the volume
    rMassMatrix *= GetGeometry().Volume();
}

template <unsigned int TDim, unsigned int TNumNodes>
void CompressibleNavierStokesExplicit<TDim, TNumNodes>::CalculateLumpedMassVector(
    VectorType& rLumpedMassVector,
    const ProcessInfo& rCurrentProcessInfo) const
{
    // Initialize the lumped mass vector
    constexpr IndexType size = TNumNodes * BlockSize;
    if (rLumpedMassVector.size() != BlockSize) {
        rLumpedMassVector.resize(size, false);
    }

    // Fill the lumped mass vector
    const double nodal_mass = GetGeometry().DomainSize() / TNumNodes;
    std::fill(rLumpedMassVector.begin(),rLumpedMassVector.end(),nodal_mass);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation

template class CompressibleNavierStokesExplicit<2,3>;
template class CompressibleNavierStokesExplicit<3,4>;

}
