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

// Application includes
#include "custom_elements/compressible_navier_stokes_explicit.h"

namespace Kratos {

template <>
void CompressibleNavierStokesExplicit<2>::EquationIdVector(
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
void CompressibleNavierStokesExplicit<3>::EquationIdVector(
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
void CompressibleNavierStokesExplicit<2>::GetDofList(
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
void CompressibleNavierStokesExplicit<3>::GetDofList(
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

template <unsigned int TDim, unsigned int TNumNodes, unsigned int TBlockSize>
int CompressibleNavierStokesExplicit<TDim, TNumNodes, TBlockSize>::Check(const ProcessInfo &rCurrentProcessInfo)
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

template <unsigned int TDim, unsigned int TNumNodes, unsigned int TBlockSize>
void CompressibleNavierStokesExplicit<TDim, TNumNodes, TBlockSize>::Calculate(
    const Variable<double>& rVariable,
    double& Output,
    const ProcessInfo& rCurrentProcessInfo)
{
    // Lumped projection terms
    if (rVariable == DENSITY_PROJECTION) {
        CalculateDensityProjection(rCurrentProcessInfo);
    } else if (rVariable == TOTAL_ENERGY_PROJECTION) {
        CalculateTotalEnergyProjection(rCurrentProcessInfo);
    } else {
        KRATOS_ERROR << "Variable not implemented." << std::endl;
    }
}

template <unsigned int TDim, unsigned int TNumNodes, unsigned int TBlockSize>
void CompressibleNavierStokesExplicit<TDim, TNumNodes, TBlockSize>::Calculate(
    const Variable<array_1d<double, 3 > >& rVariable,
    array_1d<double, 3 > & Output,
    const ProcessInfo& rCurrentProcessInfo)
{
    // Lumped projection terms
    if (rVariable == MOMENTUM_PROJECTION) {
        CalculateMomentumProjection(rCurrentProcessInfo);
    } else {
        KRATOS_ERROR << "Variable not implemented." << std::endl;
    }
}

template <unsigned int TDim, unsigned int TNumNodes, unsigned int TBlockSize>
void CompressibleNavierStokesExplicit<TDim, TNumNodes, TBlockSize>::CalculateOnIntegrationPoints(
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
    } else if (rVariable == DENSITY_SHOCK_SENSOR) {
        const double sc = this->GetValue(DENSITY_SHOCK_SENSOR);
        for (unsigned int i_gauss = 0; i_gauss < r_integration_points.size(); ++i_gauss) {
            rOutput[i_gauss] = sc;
        }
    } else if (rVariable == MOMENTUM_SHOCK_SENSOR) {
        const double sc = this->GetValue(MOMENTUM_SHOCK_SENSOR);
        for (unsigned int i_gauss = 0; i_gauss < r_integration_points.size(); ++i_gauss) {
            rOutput[i_gauss] = sc;
        }
    } else if (rVariable == TOTAL_ENERGY_SHOCK_SENSOR) {
        const double sc = this->GetValue(TOTAL_ENERGY_SHOCK_SENSOR);
        for (unsigned int i_gauss = 0; i_gauss < r_integration_points.size(); ++i_gauss) {
            rOutput[i_gauss] = sc;
        }
    } else if (rVariable == SHOCK_CAPTURING_VISCOSITY) {
        const double nu_sc = this->GetValue(SHOCK_CAPTURING_VISCOSITY);
        for (unsigned int i_gauss = 0; i_gauss < r_integration_points.size(); ++i_gauss) {
            rOutput[i_gauss] = nu_sc;
        }
    } else if (rVariable == SHOCK_CAPTURING_CONDUCTIVITY) {
        const double lambda_sc = this->GetValue(SHOCK_CAPTURING_CONDUCTIVITY);
        for (unsigned int i_gauss = 0; i_gauss < r_integration_points.size(); ++i_gauss) {
            rOutput[i_gauss] = lambda_sc;
        }
    } else {
        KRATOS_ERROR << "Variable not implemented." << std::endl;
    }
}

template <unsigned int TDim, unsigned int TNumNodes, unsigned int TBlockSize>
void CompressibleNavierStokesExplicit<TDim, TNumNodes, TBlockSize>::CalculateOnIntegrationPoints(
    const Variable<array_1d<double,3>>& rVariable,
    std::vector<array_1d<double,3>>& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    const auto& r_geometry = GetGeometry();
    const auto& r_integration_points = r_geometry.IntegrationPoints();
    if (rOutput.size() != r_integration_points.size()) {
        rOutput.resize( r_integration_points.size() );
    }

    if (rVariable == TOTAL_ENERGY_GRADIENT) {
        const auto& tot_ener_grad = this->GetValue(TOTAL_ENERGY_GRADIENT);
        for (unsigned int i_gauss = 0; i_gauss < r_integration_points.size(); ++i_gauss) {
            rOutput[i_gauss] = tot_ener_grad;
        }
    } else if (rVariable == DENSITY_GRADIENT) {
        const auto& rho_grad = this->GetValue(DENSITY_GRADIENT);
        for (unsigned int i_gauss = 0; i_gauss < r_integration_points.size(); ++i_gauss) {
            rOutput[i_gauss] = rho_grad;
        }
    } else if (rVariable == PRESSURE_GRADIENT) {
        const auto& pres_grad = this->GetValue(PRESSURE_GRADIENT);
        for (unsigned int i_gauss = 0; i_gauss < r_integration_points.size(); ++i_gauss) {
            rOutput[i_gauss] = pres_grad;
        }
    } else {
        KRATOS_ERROR << "Variable not implemented." << std::endl;
    }
}

template <unsigned int TDim, unsigned int TNumNodes, unsigned int TBlockSize>
void CompressibleNavierStokesExplicit<TDim, TNumNodes, TBlockSize>::CalculateOnIntegrationPoints(
    const Variable<Matrix>& rVariable,
    std::vector<Matrix>& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    const auto& r_geometry = GetGeometry();
    const auto& r_integration_points = r_geometry.IntegrationPoints();
    if (rOutput.size() != r_integration_points.size()) {
        rOutput.resize( r_integration_points.size() );
    }

    if (rVariable == MOMENTUM_GRADIENT) {
        const auto& mom_grad = this->GetValue(MOMENTUM_GRADIENT);
        for (unsigned int i_gauss = 0; i_gauss < r_integration_points.size(); ++i_gauss) {
            rOutput[i_gauss] = mom_grad;
        }
    } else {
        KRATOS_ERROR << "Variable not implemented." << std::endl;
    }
}

template <unsigned int TDim, unsigned int TNumNodes, unsigned int TBlockSize>
void CompressibleNavierStokesExplicit<TDim, TNumNodes, TBlockSize>::FillElementData(
    ElementDataStruct &rData,
    const ProcessInfo &rCurrentProcessInfo)
{
    // Getting data for the given geometry
    const auto& r_geometry = GetGeometry();
    GeometryUtils::CalculateGeometryData(r_geometry, rData.DN_DX, rData.N, rData.volume);

    // Compute element size
    rData.h = CalculateElementSize(rData.DN_DX);

    // Database access to all of the variables needed
    Properties &r_properties = this->GetProperties();
    rData.mu = r_properties.GetValue(DYNAMIC_VISCOSITY);
    rData.lambda = r_properties.GetValue(CONDUCTIVITY);
    rData.c_v = r_properties.GetValue(SPECIFIC_HEAT); // TODO: WE SHOULD SPECIFY WHICH ONE --> CREATE SPECIFIC_HEAT_CONSTANT_VOLUME
    rData.gamma = r_properties.GetValue(HEAT_CAPACITY_RATIO);

    rData.UseOSS = rCurrentProcessInfo[OSS_SWITCH];
    rData.ShockCapturing = rCurrentProcessInfo[SHOCK_CAPTURING_SWITCH];

    // Get nodal values
    if (rData.UseOSS) {
        for (unsigned int i = 0; i < TNumNodes; ++i) {
            const auto& r_node = r_geometry[i];
            const auto& r_momentum = r_node.FastGetSolutionStepValue(MOMENTUM);
            const auto& r_momentum_projection = r_node.GetValue(MOMENTUM_PROJECTION);
            const auto& r_momentum_time_derivative = r_node.GetValue(MOMENTUM_TIME_DERIVATIVE);
            const auto& r_body_force = r_node.FastGetSolutionStepValue(BODY_FORCE);

            for (unsigned int k = 0; k < TDim; ++k) {
                rData.U(i, k + 1) = r_momentum[k];
                rData.ResProj(i, k + 1) = r_momentum_projection[k];
                rData.dUdt(i, k + 1) = r_momentum_time_derivative[k];
                rData.f_ext(i, k) = r_body_force[k];
            }
            rData.U(i, 0) = r_node.FastGetSolutionStepValue(DENSITY);
            rData.ResProj(i, 0) = r_node.GetValue(DENSITY_PROJECTION);
            rData.dUdt(i, 0) = r_node.GetValue(DENSITY_TIME_DERIVATIVE);
            rData.U(i, TDim + 1) = r_node.FastGetSolutionStepValue(TOTAL_ENERGY);
            rData.dUdt(i, TDim + 1) = r_node.GetValue(TOTAL_ENERGY_TIME_DERIVATIVE);
            rData.r_ext(i) = r_node.FastGetSolutionStepValue(HEAT_SOURCE);
            rData.ResProj(i, TDim + 1) = r_node.GetValue(TOTAL_ENERGY_PROJECTION);
        }
    } else {
        for (unsigned int i = 0; i < TNumNodes; ++i) {
            const auto& r_node = r_geometry[i];
            const auto& r_momentum = r_node.FastGetSolutionStepValue(MOMENTUM);
            const auto& r_momentum_time_derivative = r_node.GetValue(MOMENTUM_TIME_DERIVATIVE);
            const auto& r_body_force = r_node.FastGetSolutionStepValue(BODY_FORCE);

            for (unsigned int k = 0; k < TDim; ++k) {
                rData.U(i, k + 1) = r_momentum[k];
                rData.dUdt(i, k + 1) = r_momentum_time_derivative[k];
                rData.f_ext(i, k) = r_body_force[k];
            }
            rData.U(i, 0) = r_node.FastGetSolutionStepValue(DENSITY);
            rData.dUdt(i, 0) = r_node.GetValue(DENSITY_TIME_DERIVATIVE);
            rData.U(i, TDim + 1) = r_node.FastGetSolutionStepValue(TOTAL_ENERGY);
            rData.dUdt(i, TDim + 1) = r_node.GetValue(TOTAL_ENERGY_TIME_DERIVATIVE);
            rData.r_ext(i) = r_node.FastGetSolutionStepValue(HEAT_SOURCE);
        }

    }

    // Shock capturing values
    rData.nu_sc = this->GetValue(SHOCK_CAPTURING_VISCOSITY);
    rData.lambda_sc = this->GetValue(SHOCK_CAPTURING_CONDUCTIVITY);
}

template <unsigned int TDim, unsigned int TNumNodes, unsigned int TBlockSize>
double CompressibleNavierStokesExplicit<TDim, TNumNodes, TBlockSize>::CalculateElementSize(const BoundedMatrix<double,TNumNodes, TDim>& rDN_DX)
{
    double h = 0.0;
    for (unsigned int i = 0; i < TNumNodes; ++i) {
        double h_inv = 0.0;
        for (unsigned int k = 0; k < TDim; ++k) {
            h_inv += rDN_DX(i,k) * rDN_DX(i,k);
        }
        h += 1.0/h_inv;
    }
    h = sqrt(h) / static_cast<double>(TNumNodes);
    return h;
}

template <>
void CompressibleNavierStokesExplicit<2>::CalculateMomentumProjection(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    constexpr unsigned int dim = 2;
    constexpr unsigned int n_nodes = 3;

    // Struct to pass around the data
    ElementDataStruct data;
    this->FillElementData(data, rCurrentProcessInfo);

    // Substitute the formulation symbols by the data structure values
    const BoundedMatrix<double, n_nodes, 2> &f_ext = data.f_ext;
    const double c_v = data.c_v;
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

    const double &dUdt_0_0 = data.dUdt(0, 0);
    const double &dUdt_0_1 = data.dUdt(0, 1);
    const double &dUdt_0_2 = data.dUdt(0, 2);
    const double &dUdt_0_3 = data.dUdt(0, 3);
    const double &dUdt_1_0 = data.dUdt(1, 0);
    const double &dUdt_1_1 = data.dUdt(1, 1);
    const double &dUdt_1_2 = data.dUdt(1, 2);
    const double &dUdt_1_3 = data.dUdt(1, 3);
    const double &dUdt_2_0 = data.dUdt(2, 0);
    const double &dUdt_2_1 = data.dUdt(2, 1);
    const double &dUdt_2_2 = data.dUdt(2, 2);
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
    BoundedVector<double, 6> mom_proj;

    const double cmom_proj0 =             0.16666666666666666*f_ext(1,0);
const double cmom_proj1 =             0.16666666666666666*f_ext(2,0);
const double cmom_proj2 =             cmom_proj0 + cmom_proj1 + 0.66666666666666663*f_ext(0,0);
const double cmom_proj3 =             0.16666666666666666*U_1_0;
const double cmom_proj4 =             0.16666666666666666*U_2_0;
const double cmom_proj5 =             0.66666666666666663*U_0_0 + cmom_proj3 + cmom_proj4;
const double cmom_proj6 =             0.66666666666666663*cmom_proj5;
const double cmom_proj7 =             DN_DX_0_0*U_0_0 + DN_DX_1_0*U_1_0 + DN_DX_2_0*U_2_0;
const double cmom_proj8 =             0.25*U_1_1;
const double cmom_proj9 =             0.25*U_2_1;
const double cmom_proj10 =             pow(U_0_1 + cmom_proj8 + cmom_proj9, 2);
const double cmom_proj11 =             0.25*U_1_2;
const double cmom_proj12 =             0.25*U_2_2;
const double cmom_proj13 =             pow(U_0_2 + cmom_proj11 + cmom_proj12, 2);
const double cmom_proj14 =             gamma - 1;
const double cmom_proj15 =             0.50000000000000011*cmom_proj14;
const double cmom_proj16 =             cmom_proj15*(cmom_proj10 + cmom_proj13);
const double cmom_proj17 =             -1.0000000000000002*cmom_proj10 + cmom_proj16;
const double cmom_proj18 =             0.25*U_1_0;
const double cmom_proj19 =             0.25*U_2_0;
const double cmom_proj20 =             pow(U_0_0 + cmom_proj18 + cmom_proj19, -2);
const double cmom_proj21 =             0.66666666666666663*cmom_proj20;
const double cmom_proj22 =             DN_DX_0_1*U_0_2 + DN_DX_1_1*U_1_2 + DN_DX_2_1*U_2_2;
const double cmom_proj23 =             0.16666666666666666*U_1_1;
const double cmom_proj24 =             0.16666666666666666*U_2_1;
const double cmom_proj25 =             0.66666666666666663*U_0_1 + cmom_proj23 + cmom_proj24;
const double cmom_proj26 =             1.0/cmom_proj5;
const double cmom_proj27 =             0.66666666666666663*cmom_proj26;
const double cmom_proj28 =             cmom_proj25*cmom_proj27;
const double cmom_proj29 =             DN_DX_0_1*U_0_1 + DN_DX_1_1*U_1_1 + DN_DX_2_1*U_2_1;
const double cmom_proj30 =             0.16666666666666666*U_1_2;
const double cmom_proj31 =             0.16666666666666666*U_2_2;
const double cmom_proj32 =             0.66666666666666663*U_0_2 + cmom_proj30 + cmom_proj31;
const double cmom_proj33 =             cmom_proj27*cmom_proj32;
const double cmom_proj34 =             DN_DX_0_0*U_0_2 + DN_DX_1_0*U_1_2 + DN_DX_2_0*U_2_2;
const double cmom_proj35 =             cmom_proj14*cmom_proj34;
const double cmom_proj36 =             DN_DX_0_0*U_0_1 + DN_DX_1_0*U_1_1 + DN_DX_2_0*U_2_1;
const double cmom_proj37 =             1.0*gamma - 3.0;
const double cmom_proj38 =             cmom_proj36*cmom_proj37;
const double cmom_proj39 =             DN_DX_0_1*U_0_0 + DN_DX_1_1*U_1_0 + DN_DX_2_1*U_2_0;
const double cmom_proj40 =             1.5000000000000002*cmom_proj39;
const double cmom_proj41 =             cmom_proj20*cmom_proj25*cmom_proj32;
const double cmom_proj42 =             0.16666666666666666*f_ext(0,0);
const double cmom_proj43 =             cmom_proj1 + cmom_proj42 + 0.66666666666666663*f_ext(1,0);
const double cmom_proj44 =             0.16666666666666666*U_0_0;
const double cmom_proj45 =             0.66666666666666663*U_1_0 + cmom_proj4 + cmom_proj44;
const double cmom_proj46 =             0.16666666666666666*cmom_proj45;
const double cmom_proj47 =             0.16666666666666666*cmom_proj7;
const double cmom_proj48 =             0.25*U_0_0;
const double cmom_proj49 =             pow(U_1_0 + cmom_proj19 + cmom_proj48, -2);
const double cmom_proj50 =             0.25*U_0_1;
const double cmom_proj51 =             pow(U_1_1 + cmom_proj50 + cmom_proj9, 2);
const double cmom_proj52 =             0.25*U_0_2;
const double cmom_proj53 =             pow(U_1_2 + cmom_proj12 + cmom_proj52, 2);
const double cmom_proj54 =             cmom_proj15*(cmom_proj51 + cmom_proj53);
const double cmom_proj55 =             cmom_proj49*(-1.0000000000000002*cmom_proj51 + cmom_proj54);
const double cmom_proj56 =             0.16666666666666666*U_0_1;
const double cmom_proj57 =             0.66666666666666663*U_1_1 + cmom_proj24 + cmom_proj56;
const double cmom_proj58 =             1.0/cmom_proj45;
const double cmom_proj59 =             0.16666666666666666*cmom_proj58;
const double cmom_proj60 =             cmom_proj57*cmom_proj59;
const double cmom_proj61 =             0.16666666666666666*U_0_2;
const double cmom_proj62 =             0.66666666666666663*U_1_2 + cmom_proj31 + cmom_proj61;
const double cmom_proj63 =             cmom_proj59*cmom_proj62;
const double cmom_proj64 =             0.37500000000000006*cmom_proj39;
const double cmom_proj65 =             cmom_proj49*cmom_proj57*cmom_proj62;
const double cmom_proj66 =             -cmom_proj22*cmom_proj60 - cmom_proj29*cmom_proj63 + cmom_proj35*cmom_proj63 + cmom_proj38*cmom_proj60 + cmom_proj43*cmom_proj46 - cmom_proj47*cmom_proj55 + cmom_proj64*cmom_proj65 - 0.25*dUdt_1_1;
const double cmom_proj67 =             cmom_proj0 + cmom_proj42 + 0.66666666666666663*f_ext(2,0);
const double cmom_proj68 =             0.66666666666666663*U_2_0 + cmom_proj3 + cmom_proj44;
const double cmom_proj69 =             0.16666666666666666*cmom_proj68;
const double cmom_proj70 =             DN_DX_0_0*U_0_3 + DN_DX_1_0*U_1_3 + DN_DX_2_0*U_2_3;
const double cmom_proj71 =             0.99999999999999989*cmom_proj14;
const double cmom_proj72 =             pow(U_2_0 + cmom_proj18 + cmom_proj48, -2);
const double cmom_proj73 =             pow(U_2_1 + cmom_proj50 + cmom_proj8, 2);
const double cmom_proj74 =             pow(U_2_2 + cmom_proj11 + cmom_proj52, 2);
const double cmom_proj75 =             cmom_proj15*(cmom_proj73 + cmom_proj74);
const double cmom_proj76 =             cmom_proj72*(-1.0000000000000002*cmom_proj73 + cmom_proj75);
const double cmom_proj77 =             0.66666666666666663*U_2_1 + cmom_proj23 + cmom_proj56;
const double cmom_proj78 =             1.0/cmom_proj68;
const double cmom_proj79 =             0.16666666666666666*cmom_proj78;
const double cmom_proj80 =             cmom_proj77*cmom_proj79;
const double cmom_proj81 =             0.66666666666666663*U_2_2 + cmom_proj30 + cmom_proj61;
const double cmom_proj82 =             cmom_proj79*cmom_proj81;
const double cmom_proj83 =             cmom_proj72*cmom_proj77*cmom_proj81;
const double cmom_proj84 =             -cmom_proj22*cmom_proj80 - cmom_proj29*cmom_proj82 + cmom_proj35*cmom_proj82 + cmom_proj38*cmom_proj80 - cmom_proj47*cmom_proj76 + cmom_proj64*cmom_proj83 + cmom_proj67*cmom_proj69 - cmom_proj70*cmom_proj71 - 0.25*dUdt_2_1;
const double cmom_proj85 =             0.16666666666666666*f_ext(1,1);
const double cmom_proj86 =             0.16666666666666666*f_ext(2,1);
const double cmom_proj87 =             cmom_proj85 + cmom_proj86 + 0.66666666666666663*f_ext(0,1);
const double cmom_proj88 =             -1.0000000000000002*cmom_proj13 + cmom_proj16;
const double cmom_proj89 =             cmom_proj14*cmom_proj29;
const double cmom_proj90 =             cmom_proj22*cmom_proj37;
const double cmom_proj91 =             1.5000000000000002*cmom_proj7;
const double cmom_proj92 =             0.16666666666666666*f_ext(0,1);
const double cmom_proj93 =             cmom_proj86 + cmom_proj92 + 0.66666666666666663*f_ext(1,1);
const double cmom_proj94 =             0.16666666666666666*cmom_proj39;
const double cmom_proj95 =             cmom_proj49*(-1.0000000000000002*cmom_proj53 + cmom_proj54);
const double cmom_proj96 =             0.37500000000000006*cmom_proj7;
const double cmom_proj97 =             -cmom_proj34*cmom_proj60 - cmom_proj36*cmom_proj63 + cmom_proj46*cmom_proj93 + cmom_proj60*cmom_proj89 + cmom_proj63*cmom_proj90 + cmom_proj65*cmom_proj96 - cmom_proj94*cmom_proj95 - 0.25*dUdt_1_2;
const double cmom_proj98 =             cmom_proj85 + cmom_proj92 + 0.66666666666666663*f_ext(2,1);
const double cmom_proj99 =             DN_DX_0_1*U_0_3 + DN_DX_1_1*U_1_3 + DN_DX_2_1*U_2_3;
const double cmom_proj100 =             cmom_proj72*(-1.0000000000000002*cmom_proj74 + cmom_proj75);
const double cmom_proj101 =             -cmom_proj100*cmom_proj94 - cmom_proj34*cmom_proj80 - cmom_proj36*cmom_proj82 + cmom_proj69*cmom_proj98 - cmom_proj71*cmom_proj99 + cmom_proj80*cmom_proj89 + cmom_proj82*cmom_proj90 + cmom_proj83*cmom_proj96 - 0.25*dUdt_2_2;
const double cmom_proj102 =             0.66666666666666663*cmom_proj45;
const double cmom_proj103 =             0.66666666666666663*cmom_proj7;
const double cmom_proj104 =             0.66666666666666663*cmom_proj58;
const double cmom_proj105 =             cmom_proj104*cmom_proj57;
const double cmom_proj106 =             cmom_proj104*cmom_proj62;
const double cmom_proj107 =             0.16666666666666666*cmom_proj5;
const double cmom_proj108 =             0.16666666666666666*cmom_proj26;
const double cmom_proj109 =             cmom_proj108*cmom_proj25;
const double cmom_proj110 =             cmom_proj108*cmom_proj32;
const double cmom_proj111 =             cmom_proj107*cmom_proj2 - cmom_proj109*cmom_proj22 + cmom_proj109*cmom_proj38 - cmom_proj110*cmom_proj29 + cmom_proj110*cmom_proj35 - cmom_proj17*cmom_proj20*cmom_proj47 + cmom_proj41*cmom_proj64 - 0.25*dUdt_0_1;
const double cmom_proj112 =             0.66666666666666663*cmom_proj39;
const double cmom_proj113 =             cmom_proj107*cmom_proj87 - cmom_proj109*cmom_proj34 + cmom_proj109*cmom_proj89 - cmom_proj110*cmom_proj36 + cmom_proj110*cmom_proj90 - cmom_proj20*cmom_proj88*cmom_proj94 + cmom_proj41*cmom_proj96 - 0.25*dUdt_0_2;
const double cmom_proj114 =             0.66666666666666663*cmom_proj68;
const double cmom_proj115 =             1.0*cmom_proj14;
const double cmom_proj116 =             0.66666666666666663*cmom_proj78;
const double cmom_proj117 =             cmom_proj116*cmom_proj77;
const double cmom_proj118 =             cmom_proj116*cmom_proj81;
            mom_proj[0]=-cmom_proj17*cmom_proj21*cmom_proj7 + cmom_proj2*cmom_proj6 - cmom_proj22*cmom_proj28 + cmom_proj28*cmom_proj38 - cmom_proj29*cmom_proj33 + cmom_proj33*cmom_proj35 + cmom_proj40*cmom_proj41 + cmom_proj66 + cmom_proj84 - 0.5*dUdt_0_1;
            mom_proj[1]=cmom_proj101 - cmom_proj21*cmom_proj39*cmom_proj88 - cmom_proj28*cmom_proj34 + cmom_proj28*cmom_proj89 - cmom_proj33*cmom_proj36 + cmom_proj33*cmom_proj90 + cmom_proj41*cmom_proj91 + cmom_proj6*cmom_proj87 + cmom_proj97 - 0.5*dUdt_0_2;
            mom_proj[2]=cmom_proj102*cmom_proj43 - cmom_proj103*cmom_proj55 - cmom_proj105*cmom_proj22 + cmom_proj105*cmom_proj38 - cmom_proj106*cmom_proj29 + cmom_proj106*cmom_proj35 + cmom_proj111 + cmom_proj40*cmom_proj65 + cmom_proj84 - 0.5*dUdt_1_1;
            mom_proj[3]=cmom_proj101 + cmom_proj102*cmom_proj93 - cmom_proj105*cmom_proj34 + cmom_proj105*cmom_proj89 - cmom_proj106*cmom_proj36 + cmom_proj106*cmom_proj90 - cmom_proj112*cmom_proj95 + cmom_proj113 + cmom_proj65*cmom_proj91 - 0.5*dUdt_1_2;
            mom_proj[4]=-cmom_proj103*cmom_proj76 + cmom_proj111 + cmom_proj114*cmom_proj67 - cmom_proj115*cmom_proj70 - cmom_proj117*cmom_proj22 + cmom_proj117*cmom_proj38 - cmom_proj118*cmom_proj29 + cmom_proj118*cmom_proj35 + cmom_proj40*cmom_proj83 + cmom_proj66 - 0.5*dUdt_2_1;
            mom_proj[5]=-cmom_proj100*cmom_proj112 + cmom_proj113 + cmom_proj114*cmom_proj98 - cmom_proj115*cmom_proj99 - cmom_proj117*cmom_proj34 + cmom_proj117*cmom_proj89 - cmom_proj118*cmom_proj36 + cmom_proj118*cmom_proj90 + cmom_proj83*cmom_proj91 + cmom_proj97 - 0.5*dUdt_2_2;

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
void CompressibleNavierStokesExplicit<3>::CalculateMomentumProjection(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    constexpr unsigned int dim = 3;
    constexpr unsigned int n_nodes = 4;

    // Struct to pass around the data
    ElementDataStruct data;
    this->FillElementData(data, rCurrentProcessInfo);

    // Substitute the formulation symbols by the data structure values
    const BoundedMatrix<double, n_nodes, 3> &f_ext = data.f_ext;
    const double c_v = data.c_v;
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

    const double &dUdt_0_0 = data.dUdt(0, 0);
    const double &dUdt_0_1 = data.dUdt(0, 1);
    const double &dUdt_0_2 = data.dUdt(0, 2);
    const double &dUdt_0_3 = data.dUdt(0, 3);
    const double &dUdt_0_4 = data.dUdt(0, 4);
    const double &dUdt_1_0 = data.dUdt(1, 0);
    const double &dUdt_1_1 = data.dUdt(1, 1);
    const double &dUdt_1_2 = data.dUdt(1, 2);
    const double &dUdt_1_3 = data.dUdt(1, 3);
    const double &dUdt_1_4 = data.dUdt(1, 4);
    const double &dUdt_2_0 = data.dUdt(2, 0);
    const double &dUdt_2_1 = data.dUdt(2, 1);
    const double &dUdt_2_2 = data.dUdt(2, 2);
    const double &dUdt_2_3 = data.dUdt(2, 3);
    const double &dUdt_2_4 = data.dUdt(2, 4);
    const double &dUdt_3_0 = data.dUdt(3, 0);
    const double &dUdt_3_1 = data.dUdt(3, 1);
    const double &dUdt_3_2 = data.dUdt(3, 2);
    const double &dUdt_3_3 = data.dUdt(3, 3);
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
void CompressibleNavierStokesExplicit<2>::CalculateDensityProjection(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    constexpr unsigned int n_nodes = 3;

    // Struct to pass around the data
    ElementDataStruct data;
    this->FillElementData(data, rCurrentProcessInfo);

    // Substitute the formulation symbols by the data structure values
    const double h = data.h;
    const array_1d<double, n_nodes> &r_ext = data.r_ext;
    const BoundedMatrix<double, n_nodes, 2> &f_ext = data.f_ext;
    const double mu = data.mu;
    const double lambda = data.lambda;
    const double c_v = data.c_v;
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

    const double &dUdt_0_0 = data.dUdt(0, 0);
    const double &dUdt_0_1 = data.dUdt(0, 1);
    const double &dUdt_0_2 = data.dUdt(0, 2);
    const double &dUdt_0_3 = data.dUdt(0, 3);
    const double &dUdt_1_0 = data.dUdt(1, 0);
    const double &dUdt_1_1 = data.dUdt(1, 1);
    const double &dUdt_1_2 = data.dUdt(1, 2);
    const double &dUdt_1_3 = data.dUdt(1, 3);
    const double &dUdt_2_0 = data.dUdt(2, 0);
    const double &dUdt_2_1 = data.dUdt(2, 1);
    const double &dUdt_2_2 = data.dUdt(2, 2);
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
    BoundedVector<double, 3> rho_proj;

    const double crho_proj0 =             0.25*dUdt_1_0;
const double crho_proj1 =             DN_DX_0_0*U_0_1;
const double crho_proj2 =             DN_DX_0_1*U_0_2;
const double crho_proj3 =             DN_DX_1_0*U_1_1;
const double crho_proj4 =             DN_DX_1_1*U_1_2;
const double crho_proj5 =             DN_DX_2_0*U_2_1;
const double crho_proj6 =             DN_DX_2_1*U_2_2;
const double crho_proj7 =             0.99999999999999989*crho_proj1 + 0.99999999999999989*crho_proj2 + 0.99999999999999989*crho_proj3 + 0.99999999999999989*crho_proj4 + 0.99999999999999989*crho_proj5 + 0.99999999999999989*crho_proj6 + 0.25*dUdt_2_0;
const double crho_proj8 =             0.25*dUdt_0_0;
            rho_proj[0]=-crho_proj0 - crho_proj7 - 0.5*dUdt_0_0;
            rho_proj[1]=-crho_proj7 - crho_proj8 - 0.5*dUdt_1_0;
            rho_proj[2]=-crho_proj0 - 1.0*crho_proj1 - 1.0*crho_proj2 - 1.0*crho_proj3 - 1.0*crho_proj4 - 1.0*crho_proj5 - 1.0*crho_proj6 - crho_proj8 - 0.5*dUdt_2_0;

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
void CompressibleNavierStokesExplicit<3>::CalculateDensityProjection(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    constexpr unsigned int dim = 3;
    constexpr unsigned int n_nodes = 4;

    // Struct to pass around the data
    ElementDataStruct data;
    this->FillElementData(data, rCurrentProcessInfo);

    // Substitute the formulation symbols by the data structure values
    const double h = data.h;
    const array_1d<double, n_nodes> &r_ext = data.r_ext;
    const BoundedMatrix<double, n_nodes, 3> &f_ext = data.f_ext;
    const double mu = data.mu;
    const double lambda = data.lambda;
    const double c_v = data.c_v;
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

    const double &dUdt_0_0 = data.dUdt(0, 0);
    const double &dUdt_0_1 = data.dUdt(0, 1);
    const double &dUdt_0_2 = data.dUdt(0, 2);
    const double &dUdt_0_3 = data.dUdt(0, 3);
    const double &dUdt_0_4 = data.dUdt(0, 4);
    const double &dUdt_1_0 = data.dUdt(1, 0);
    const double &dUdt_1_1 = data.dUdt(1, 1);
    const double &dUdt_1_2 = data.dUdt(1, 2);
    const double &dUdt_1_3 = data.dUdt(1, 3);
    const double &dUdt_1_4 = data.dUdt(1, 4);
    const double &dUdt_2_0 = data.dUdt(2, 0);
    const double &dUdt_2_1 = data.dUdt(2, 1);
    const double &dUdt_2_2 = data.dUdt(2, 2);
    const double &dUdt_2_3 = data.dUdt(2, 3);
    const double &dUdt_2_4 = data.dUdt(2, 4);
    const double &dUdt_3_0 = data.dUdt(3, 0);
    const double &dUdt_3_1 = data.dUdt(3, 1);
    const double &dUdt_3_2 = data.dUdt(3, 2);
    const double &dUdt_3_3 = data.dUdt(3, 3);
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
void CompressibleNavierStokesExplicit<2>::CalculateTotalEnergyProjection(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    constexpr unsigned int n_nodes = 3;

    // Struct to pass around the data
    ElementDataStruct data;
    this->FillElementData(data, rCurrentProcessInfo);

    // Substitute the formulation symbols by the data structure values
    const double h = data.h;
    const array_1d<double, n_nodes> &r_ext = data.r_ext;
    const BoundedMatrix<double, n_nodes, 2> &f_ext = data.f_ext;
    const double mu = data.mu;
    const double lambda = data.lambda;
    const double c_v = data.c_v;
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

    const double &dUdt_0_0 = data.dUdt(0, 0);
    const double &dUdt_0_1 = data.dUdt(0, 1);
    const double &dUdt_0_2 = data.dUdt(0, 2);
    const double &dUdt_0_3 = data.dUdt(0, 3);
    const double &dUdt_1_0 = data.dUdt(1, 0);
    const double &dUdt_1_1 = data.dUdt(1, 1);
    const double &dUdt_1_2 = data.dUdt(1, 2);
    const double &dUdt_1_3 = data.dUdt(1, 3);
    const double &dUdt_2_0 = data.dUdt(2, 0);
    const double &dUdt_2_1 = data.dUdt(2, 1);
    const double &dUdt_2_2 = data.dUdt(2, 2);
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

    const double ctot_ener_proj0 =             0.25*U_1_1;
const double ctot_ener_proj1 =             0.25*U_2_1;
const double ctot_ener_proj2 =             pow(U_0_1 + ctot_ener_proj0 + ctot_ener_proj1, 2);
const double ctot_ener_proj3 =             0.25*U_1_0;
const double ctot_ener_proj4 =             0.25*U_2_0;
const double ctot_ener_proj5 =             U_0_0 + ctot_ener_proj3 + ctot_ener_proj4;
const double ctot_ener_proj6 =             pow(ctot_ener_proj5, -2);
const double ctot_ener_proj7 =             gamma - 1;
const double ctot_ener_proj8 =             1.0000000000000002*ctot_ener_proj7;
const double ctot_ener_proj9 =             ctot_ener_proj6*ctot_ener_proj8;
const double ctot_ener_proj10 =             0.16666666666666666*U_1_0;
const double ctot_ener_proj11 =             0.16666666666666666*U_2_0;
const double ctot_ener_proj12 =             0.66666666666666663*U_0_0 + ctot_ener_proj10 + ctot_ener_proj11;
const double ctot_ener_proj13 =             1.0/ctot_ener_proj12;
const double ctot_ener_proj14 =             0.25*U_1_2;
const double ctot_ener_proj15 =             0.25*U_2_2;
const double ctot_ener_proj16 =             pow(U_0_2 + ctot_ener_proj14 + ctot_ener_proj15, 2);
const double ctot_ener_proj17 =             0.16666666666666666*U_1_3;
const double ctot_ener_proj18 =             0.16666666666666666*U_2_3;
const double ctot_ener_proj19 =             0.66666666666666663*U_0_3 + ctot_ener_proj17 + ctot_ener_proj18;
const double ctot_ener_proj20 =             ctot_ener_proj7*(-ctot_ener_proj13*(0.22222222222222221*ctot_ener_proj16 + 0.22222222222222221*ctot_ener_proj2) + ctot_ener_proj19);
const double ctot_ener_proj21 =             -ctot_ener_proj13*(ctot_ener_proj19 + ctot_ener_proj20);
const double ctot_ener_proj22 =             ctot_ener_proj2*ctot_ener_proj9 + ctot_ener_proj21;
const double ctot_ener_proj23 =             DN_DX_0_0*U_0_1 + DN_DX_1_0*U_1_1 + DN_DX_2_0*U_2_1;
const double ctot_ener_proj24 =             0.66666666666666663*ctot_ener_proj23;
const double ctot_ener_proj25 =             ctot_ener_proj16*ctot_ener_proj9 + ctot_ener_proj21;
const double ctot_ener_proj26 =             DN_DX_0_1*U_0_2 + DN_DX_1_1*U_1_2 + DN_DX_2_1*U_2_2;
const double ctot_ener_proj27 =             0.66666666666666663*ctot_ener_proj26;
const double ctot_ener_proj28 =             0.16666666666666666*r_ext[1];
const double ctot_ener_proj29 =             0.16666666666666666*r_ext[2];
const double ctot_ener_proj30 =             ctot_ener_proj12*(ctot_ener_proj28 + ctot_ener_proj29 + 0.66666666666666663*r_ext[0]);
const double ctot_ener_proj31 =             0.16666666666666666*U_1_1;
const double ctot_ener_proj32 =             0.16666666666666666*U_2_1;
const double ctot_ener_proj33 =             0.66666666666666663*U_0_1 + ctot_ener_proj31 + ctot_ener_proj32;
const double ctot_ener_proj34 =             0.16666666666666666*f_ext(1,0);
const double ctot_ener_proj35 =             0.16666666666666666*f_ext(2,0);
const double ctot_ener_proj36 =             ctot_ener_proj33*(ctot_ener_proj34 + ctot_ener_proj35 + 0.66666666666666663*f_ext(0,0));
const double ctot_ener_proj37 =             0.16666666666666666*U_1_2;
const double ctot_ener_proj38 =             0.16666666666666666*U_2_2;
const double ctot_ener_proj39 =             0.66666666666666663*U_0_2 + ctot_ener_proj37 + ctot_ener_proj38;
const double ctot_ener_proj40 =             0.16666666666666666*f_ext(1,1);
const double ctot_ener_proj41 =             0.16666666666666666*f_ext(2,1);
const double ctot_ener_proj42 =             ctot_ener_proj39*(ctot_ener_proj40 + ctot_ener_proj41 + 0.66666666666666663*f_ext(0,1));
const double ctot_ener_proj43 =             DN_DX_0_0*U_0_3 + DN_DX_1_0*U_1_3 + DN_DX_2_0*U_2_3;
const double ctot_ener_proj44 =             0.66666666666666663*ctot_ener_proj33;
const double ctot_ener_proj45 =             ctot_ener_proj13*gamma;
const double ctot_ener_proj46 =             DN_DX_0_1*U_0_3 + DN_DX_1_1*U_1_3 + DN_DX_2_1*U_2_3;
const double ctot_ener_proj47 =             0.66666666666666663*ctot_ener_proj39;
const double ctot_ener_proj48 =             DN_DX_0_0*U_0_0 + DN_DX_1_0*U_1_0 + DN_DX_2_0*U_2_0;
const double ctot_ener_proj49 =             -0.37500000000000006*U_1_3;
const double ctot_ener_proj50 =             -0.37500000000000006*U_2_3;
const double ctot_ener_proj51 =             0.75000000000000011*ctot_ener_proj7;
const double ctot_ener_proj52 =             ctot_ener_proj6*(-1.5000000000000002*U_0_3 - 2.2500000000000004*ctot_ener_proj20 + ctot_ener_proj49 + ctot_ener_proj50 + ctot_ener_proj51*(ctot_ener_proj16 + ctot_ener_proj2)/ctot_ener_proj5);
const double ctot_ener_proj53 =             DN_DX_0_1*U_0_0 + DN_DX_1_1*U_1_0 + DN_DX_2_1*U_2_0;
const double ctot_ener_proj54 =             ctot_ener_proj33*ctot_ener_proj39*ctot_ener_proj6;
const double ctot_ener_proj55 =             DN_DX_0_0*U_0_2 + DN_DX_1_0*U_1_2 + DN_DX_2_0*U_2_2;
const double ctot_ener_proj56 =             1.5000000000000002*ctot_ener_proj7;
const double ctot_ener_proj57 =             ctot_ener_proj55*ctot_ener_proj56;
const double ctot_ener_proj58 =             DN_DX_0_1*U_0_1 + DN_DX_1_1*U_1_1 + DN_DX_2_1*U_2_1;
const double ctot_ener_proj59 =             ctot_ener_proj56*ctot_ener_proj58;
const double ctot_ener_proj60 =             0.25*U_0_1;
const double ctot_ener_proj61 =             pow(U_2_1 + ctot_ener_proj0 + ctot_ener_proj60, 2);
const double ctot_ener_proj62 =             0.25*U_0_0;
const double ctot_ener_proj63 =             U_2_0 + ctot_ener_proj3 + ctot_ener_proj62;
const double ctot_ener_proj64 =             pow(ctot_ener_proj63, -2);
const double ctot_ener_proj65 =             ctot_ener_proj64*ctot_ener_proj8;
const double ctot_ener_proj66 =             0.16666666666666666*U_0_0;
const double ctot_ener_proj67 =             0.66666666666666663*U_2_0 + ctot_ener_proj10 + ctot_ener_proj66;
const double ctot_ener_proj68 =             1.0/ctot_ener_proj67;
const double ctot_ener_proj69 =             0.25*U_0_2;
const double ctot_ener_proj70 =             pow(U_2_2 + ctot_ener_proj14 + ctot_ener_proj69, 2);
const double ctot_ener_proj71 =             0.16666666666666666*U_0_3;
const double ctot_ener_proj72 =             0.66666666666666663*U_2_3 + ctot_ener_proj17 + ctot_ener_proj71;
const double ctot_ener_proj73 =             ctot_ener_proj7*(-ctot_ener_proj68*(0.22222222222222221*ctot_ener_proj61 + 0.22222222222222221*ctot_ener_proj70) + ctot_ener_proj72);
const double ctot_ener_proj74 =             -ctot_ener_proj68*(ctot_ener_proj72 + ctot_ener_proj73);
const double ctot_ener_proj75 =             ctot_ener_proj61*ctot_ener_proj65 + ctot_ener_proj74;
const double ctot_ener_proj76 =             0.16666666666666666*ctot_ener_proj23;
const double ctot_ener_proj77 =             ctot_ener_proj65*ctot_ener_proj70 + ctot_ener_proj74;
const double ctot_ener_proj78 =             0.16666666666666666*ctot_ener_proj26;
const double ctot_ener_proj79 =             0.16666666666666666*r_ext[0];
const double ctot_ener_proj80 =             ctot_ener_proj67*(ctot_ener_proj28 + ctot_ener_proj79 + 0.66666666666666663*r_ext[2]);
const double ctot_ener_proj81 =             0.16666666666666666*U_0_1;
const double ctot_ener_proj82 =             0.66666666666666663*U_2_1 + ctot_ener_proj31 + ctot_ener_proj81;
const double ctot_ener_proj83 =             0.16666666666666666*f_ext(0,0);
const double ctot_ener_proj84 =             ctot_ener_proj82*(ctot_ener_proj34 + ctot_ener_proj83 + 0.66666666666666663*f_ext(2,0));
const double ctot_ener_proj85 =             0.16666666666666666*U_0_2;
const double ctot_ener_proj86 =             0.66666666666666663*U_2_2 + ctot_ener_proj37 + ctot_ener_proj85;
const double ctot_ener_proj87 =             0.16666666666666666*f_ext(0,1);
const double ctot_ener_proj88 =             ctot_ener_proj86*(ctot_ener_proj40 + ctot_ener_proj87 + 0.66666666666666663*f_ext(2,1));
const double ctot_ener_proj89 =             0.16666666666666666*gamma;
const double ctot_ener_proj90 =             ctot_ener_proj68*ctot_ener_proj89;
const double ctot_ener_proj91 =             ctot_ener_proj43*ctot_ener_proj82;
const double ctot_ener_proj92 =             ctot_ener_proj46*ctot_ener_proj86;
const double ctot_ener_proj93 =             0.16666666666666666*ctot_ener_proj48;
const double ctot_ener_proj94 =             -0.37500000000000006*U_0_3;
const double ctot_ener_proj95 =             ctot_ener_proj64*(-1.5000000000000002*U_2_3 + ctot_ener_proj49 + ctot_ener_proj51*(ctot_ener_proj61 + ctot_ener_proj70)/ctot_ener_proj63 - 2.2500000000000004*ctot_ener_proj73 + ctot_ener_proj94);
const double ctot_ener_proj96 =             ctot_ener_proj82*ctot_ener_proj95;
const double ctot_ener_proj97 =             0.16666666666666666*ctot_ener_proj53;
const double ctot_ener_proj98 =             ctot_ener_proj86*ctot_ener_proj95;
const double ctot_ener_proj99 =             0.37500000000000006*ctot_ener_proj7;
const double ctot_ener_proj100 =             ctot_ener_proj55*ctot_ener_proj99;
const double ctot_ener_proj101 =             ctot_ener_proj64*ctot_ener_proj82*ctot_ener_proj86;
const double ctot_ener_proj102 =             ctot_ener_proj58*ctot_ener_proj99;
const double ctot_ener_proj103 =             ctot_ener_proj100*ctot_ener_proj101 + ctot_ener_proj101*ctot_ener_proj102 + ctot_ener_proj75*ctot_ener_proj76 + ctot_ener_proj77*ctot_ener_proj78 + 0.16666666666666666*ctot_ener_proj80 + 0.16666666666666666*ctot_ener_proj84 + 0.16666666666666666*ctot_ener_proj88 - ctot_ener_proj90*ctot_ener_proj91 - ctot_ener_proj90*ctot_ener_proj92 - ctot_ener_proj93*ctot_ener_proj96 - ctot_ener_proj97*ctot_ener_proj98 - 0.25*dUdt_2_3;
const double ctot_ener_proj104 =             pow(U_1_1 + ctot_ener_proj1 + ctot_ener_proj60, 2);
const double ctot_ener_proj105 =             U_1_0 + ctot_ener_proj4 + ctot_ener_proj62;
const double ctot_ener_proj106 =             pow(ctot_ener_proj105, -2);
const double ctot_ener_proj107 =             ctot_ener_proj106*ctot_ener_proj8;
const double ctot_ener_proj108 =             0.66666666666666663*U_1_0 + ctot_ener_proj11 + ctot_ener_proj66;
const double ctot_ener_proj109 =             1.0/ctot_ener_proj108;
const double ctot_ener_proj110 =             pow(U_1_2 + ctot_ener_proj15 + ctot_ener_proj69, 2);
const double ctot_ener_proj111 =             0.66666666666666663*U_1_3 + ctot_ener_proj18 + ctot_ener_proj71;
const double ctot_ener_proj112 =             ctot_ener_proj7*(-ctot_ener_proj109*(0.22222222222222221*ctot_ener_proj104 + 0.22222222222222221*ctot_ener_proj110) + ctot_ener_proj111);
const double ctot_ener_proj113 =             -ctot_ener_proj109*(ctot_ener_proj111 + ctot_ener_proj112);
const double ctot_ener_proj114 =             ctot_ener_proj104*ctot_ener_proj107 + ctot_ener_proj113;
const double ctot_ener_proj115 =             ctot_ener_proj107*ctot_ener_proj110 + ctot_ener_proj113;
const double ctot_ener_proj116 =             ctot_ener_proj108*(ctot_ener_proj29 + ctot_ener_proj79 + 0.66666666666666663*r_ext[1]);
const double ctot_ener_proj117 =             0.66666666666666663*U_1_1 + ctot_ener_proj32 + ctot_ener_proj81;
const double ctot_ener_proj118 =             ctot_ener_proj117*(ctot_ener_proj35 + ctot_ener_proj83 + 0.66666666666666663*f_ext(1,0));
const double ctot_ener_proj119 =             0.66666666666666663*U_1_2 + ctot_ener_proj38 + ctot_ener_proj85;
const double ctot_ener_proj120 =             ctot_ener_proj119*(ctot_ener_proj41 + ctot_ener_proj87 + 0.66666666666666663*f_ext(1,1));
const double ctot_ener_proj121 =             ctot_ener_proj109*ctot_ener_proj89;
const double ctot_ener_proj122 =             ctot_ener_proj117*ctot_ener_proj43;
const double ctot_ener_proj123 =             ctot_ener_proj119*ctot_ener_proj46;
const double ctot_ener_proj124 =             ctot_ener_proj106*(-1.5000000000000002*U_1_3 - 2.2500000000000004*ctot_ener_proj112 + ctot_ener_proj50 + ctot_ener_proj94 + ctot_ener_proj51*(ctot_ener_proj104 + ctot_ener_proj110)/ctot_ener_proj105);
const double ctot_ener_proj125 =             ctot_ener_proj117*ctot_ener_proj124;
const double ctot_ener_proj126 =             ctot_ener_proj119*ctot_ener_proj124;
const double ctot_ener_proj127 =             ctot_ener_proj106*ctot_ener_proj117*ctot_ener_proj119;
const double ctot_ener_proj128 =             ctot_ener_proj100*ctot_ener_proj127 + ctot_ener_proj102*ctot_ener_proj127 + ctot_ener_proj114*ctot_ener_proj76 + ctot_ener_proj115*ctot_ener_proj78 + 0.16666666666666666*ctot_ener_proj116 + 0.16666666666666666*ctot_ener_proj118 + 0.16666666666666666*ctot_ener_proj120 - ctot_ener_proj121*ctot_ener_proj122 - ctot_ener_proj121*ctot_ener_proj123 - ctot_ener_proj125*ctot_ener_proj93 - ctot_ener_proj126*ctot_ener_proj97 - 0.25*dUdt_1_3;
const double ctot_ener_proj129 =             0.66666666666666663*gamma;
const double ctot_ener_proj130 =             ctot_ener_proj109*ctot_ener_proj129;
const double ctot_ener_proj131 =             0.66666666666666663*ctot_ener_proj48;
const double ctot_ener_proj132 =             0.66666666666666663*ctot_ener_proj53;
const double ctot_ener_proj133 =             ctot_ener_proj13*ctot_ener_proj89;
const double ctot_ener_proj134 =             ctot_ener_proj100*ctot_ener_proj54 + ctot_ener_proj102*ctot_ener_proj54 - ctot_ener_proj133*ctot_ener_proj33*ctot_ener_proj43 - ctot_ener_proj133*ctot_ener_proj39*ctot_ener_proj46 + ctot_ener_proj22*ctot_ener_proj76 + ctot_ener_proj25*ctot_ener_proj78 + 0.16666666666666666*ctot_ener_proj30 - ctot_ener_proj33*ctot_ener_proj52*ctot_ener_proj93 + 0.16666666666666666*ctot_ener_proj36 - ctot_ener_proj39*ctot_ener_proj52*ctot_ener_proj97 + 0.16666666666666666*ctot_ener_proj42 - 0.25*dUdt_0_3;
const double ctot_ener_proj135 =             ctot_ener_proj129*ctot_ener_proj68;
            tot_ener_proj[0]=ctot_ener_proj103 + ctot_ener_proj128 + ctot_ener_proj22*ctot_ener_proj24 + ctot_ener_proj25*ctot_ener_proj27 + 0.66666666666666663*ctot_ener_proj30 + 0.66666666666666663*ctot_ener_proj36 + 0.66666666666666663*ctot_ener_proj42 - ctot_ener_proj43*ctot_ener_proj44*ctot_ener_proj45 - ctot_ener_proj44*ctot_ener_proj48*ctot_ener_proj52 - ctot_ener_proj45*ctot_ener_proj46*ctot_ener_proj47 - ctot_ener_proj47*ctot_ener_proj52*ctot_ener_proj53 + ctot_ener_proj54*ctot_ener_proj57 + ctot_ener_proj54*ctot_ener_proj59 - 0.5*dUdt_0_3;
            tot_ener_proj[1]=ctot_ener_proj103 + ctot_ener_proj114*ctot_ener_proj24 + ctot_ener_proj115*ctot_ener_proj27 + 0.66666666666666663*ctot_ener_proj116 + 0.66666666666666663*ctot_ener_proj118 + 0.66666666666666663*ctot_ener_proj120 - ctot_ener_proj122*ctot_ener_proj130 - ctot_ener_proj123*ctot_ener_proj130 - ctot_ener_proj125*ctot_ener_proj131 - ctot_ener_proj126*ctot_ener_proj132 + ctot_ener_proj127*ctot_ener_proj57 + ctot_ener_proj127*ctot_ener_proj59 + ctot_ener_proj134 - 0.5*dUdt_1_3;
            tot_ener_proj[2]=ctot_ener_proj101*ctot_ener_proj57 + ctot_ener_proj101*ctot_ener_proj59 + ctot_ener_proj128 - ctot_ener_proj131*ctot_ener_proj96 - ctot_ener_proj132*ctot_ener_proj98 + ctot_ener_proj134 - ctot_ener_proj135*ctot_ener_proj91 - ctot_ener_proj135*ctot_ener_proj92 + ctot_ener_proj24*ctot_ener_proj75 + ctot_ener_proj27*ctot_ener_proj77 + 0.66666666666666663*ctot_ener_proj80 + 0.66666666666666663*ctot_ener_proj84 + 0.66666666666666663*ctot_ener_proj88 - 0.5*dUdt_2_3;

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
void CompressibleNavierStokesExplicit<3>::CalculateTotalEnergyProjection(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    constexpr unsigned int n_nodes = 4;

    // Struct to pass around the data
    ElementDataStruct data;
    this->FillElementData(data, rCurrentProcessInfo);

    // Substitute the formulation symbols by the data structure values
    const double h = data.h;
    const array_1d<double, n_nodes> &r_ext = data.r_ext;
    const BoundedMatrix<double, n_nodes, 3> &f_ext = data.f_ext;
    const double mu = data.mu;
    const double lambda = data.lambda;
    const double c_v = data.c_v;
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

    const double &dUdt_0_0 = data.dUdt(0, 0);
    const double &dUdt_0_1 = data.dUdt(0, 1);
    const double &dUdt_0_2 = data.dUdt(0, 2);
    const double &dUdt_0_3 = data.dUdt(0, 3);
    const double &dUdt_0_4 = data.dUdt(0, 4);
    const double &dUdt_1_0 = data.dUdt(1, 0);
    const double &dUdt_1_1 = data.dUdt(1, 1);
    const double &dUdt_1_2 = data.dUdt(1, 2);
    const double &dUdt_1_3 = data.dUdt(1, 3);
    const double &dUdt_1_4 = data.dUdt(1, 4);
    const double &dUdt_2_0 = data.dUdt(2, 0);
    const double &dUdt_2_1 = data.dUdt(2, 1);
    const double &dUdt_2_2 = data.dUdt(2, 2);
    const double &dUdt_2_3 = data.dUdt(2, 3);
    const double &dUdt_2_4 = data.dUdt(2, 4);
    const double &dUdt_3_0 = data.dUdt(3, 0);
    const double &dUdt_3_1 = data.dUdt(3, 1);
    const double &dUdt_3_2 = data.dUdt(3, 2);
    const double &dUdt_3_3 = data.dUdt(3, 3);
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
void CompressibleNavierStokesExplicit<2>::CalculateRightHandSideInternal(
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
    const array_1d<double, n_nodes> &r_ext = data.r_ext;
    const BoundedMatrix<double, n_nodes, 2> &f_ext = data.f_ext;
    const double mu = data.mu;
    const double lambda = data.lambda;
    const double c_v = data.c_v;
    const double gamma = data.gamma;

    // Stabilization parameters
    const double stab_c1 = 12.0;
    const double stab_c2 = 2.0;
    const double stab_c3 = 1.0;

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

    const double &dUdt_0_0 = data.dUdt(0, 0);
    const double &dUdt_0_1 = data.dUdt(0, 1);
    const double &dUdt_0_2 = data.dUdt(0, 2);
    const double &dUdt_0_3 = data.dUdt(0, 3);
    const double &dUdt_1_0 = data.dUdt(1, 0);
    const double &dUdt_1_1 = data.dUdt(1, 1);
    const double &dUdt_1_2 = data.dUdt(1, 2);
    const double &dUdt_1_3 = data.dUdt(1, 3);
    const double &dUdt_2_0 = data.dUdt(2, 0);
    const double &dUdt_2_1 = data.dUdt(2, 1);
    const double &dUdt_2_2 = data.dUdt(2, 2);
    const double &dUdt_2_3 = data.dUdt(2, 3);

    // Shock capturing parameters
    double nu_sc = 0.0;
    double nu_st = 0.0;
    double k_sc = 0.0;
    double k_st = 0.0;
    double lin_m_norm = 1.0; // This is intentionally set to a non-zero number to avoid dividing by zero
    array_1d<double, 2> lin_m = ZeroVector(2);
    if (data.ShockCapturing) {
        nu_sc = data.nu_sc;
        k_sc = data.lambda_sc;

        const double rho_avg = (U_0_0 + U_1_0 + U_2_0) / 3.0;
        array_1d<double, 2> v_avg;
        v_avg[0] = (U_0_1 / U_0_0 + U_1_1 / U_1_0 + U_2_1 / U_2_0) / 3.0;
        v_avg[1] = (U_0_2 / U_0_0 + U_1_2 / U_1_0 + U_2_2 / U_2_0) / 3.0;
        const double v_avg_norm = norm_2(v_avg);
        const double tot_ener_avg = (U_0_3 + U_1_3 + U_2_3) / 3.0;
        const double c_avg = gamma * (gamma - 1.0) * ((tot_ener_avg / rho_avg) - 0.5 * std::pow(v_avg_norm, 2));

        // Source terms midpoint values
        double r_ext_avg = 0.0;
        array_1d<double, 2> f_ext_avg = ZeroVector(2);
        for (unsigned int i_node = 0; i_node < n_nodes; ++i_node) {
            f_ext_avg[0] += f_ext(i_node, 0);
            f_ext_avg[1] += f_ext(i_node, 1);
            r_ext_avg += r_ext[i_node];
        }
        f_ext_avg /= n_nodes;
        r_ext_avg /= n_nodes;
        const double f_ext_avg_norm = norm_2(f_ext_avg);

        const double alpha = lambda / (rho_avg * gamma * c_v);
        const double aux_1 = std::pow(r_ext_avg, 2) + 2.0 * std::pow(c_avg, 2) * std::pow(f_ext_avg_norm, 2) + std::sqrt(std::pow(r_ext_avg, 4) + 4.0 * std::pow(c_avg, 2) * std::pow(f_ext_avg_norm, 2) * std::pow(r_ext_avg, 2));
        const double aux_2 = 2.0 * std::pow(c_avg, 4);
        const double tau_rho = (stab_c2 * (v_avg_norm + c_avg)/ h) + (stab_c3 * std::sqrt(aux_1 / aux_2));
        const double tau_m_avg = 1.0 / ((4.0 * stab_c1 * mu / 3.0 / rho_avg / std::pow(h, 2)) + tau_rho);
        const double tau_et_avg = 1.0 / ((stab_c1 * alpha / std::pow(h, 2)) + tau_rho);

        nu_st = std::max(0.0, nu_sc - tau_m_avg * std::pow(v_avg_norm, 2));
        k_st = std::max(0.0, k_sc - tau_et_avg * std::pow(v_avg_norm, 2));

        lin_m[0] = (U_0_1 + U_1_1 + U_2_1) / 3.0;
        lin_m[1] = (U_0_2 + U_1_2 + U_2_2) / 3.0;
        const double zero_tol = 1.0e-12;
        const double aux = norm_2(lin_m);
        lin_m_norm = aux > zero_tol ? aux : 1.0;
    }

    // Hardcoded shape functions gradients for linear triangular element
    // This is explicitly done to minimize the matrix acceses
    // The notation DN_i_j means shape function for node i in dimension j
    const double &DN_DX_0_0 = data.DN_DX(0, 0);
    const double &DN_DX_0_1 = data.DN_DX(0, 1);
    const double &DN_DX_1_0 = data.DN_DX(1, 0);
    const double &DN_DX_1_1 = data.DN_DX(1, 1);
    const double &DN_DX_2_0 = data.DN_DX(2, 0);
    const double &DN_DX_2_1 = data.DN_DX(2, 1);

    if (data.UseOSS) {
        // Projections container accesses
        const double &ResProj_0_0 = data.ResProj(0, 0);
        const double &ResProj_0_1 = data.ResProj(0, 1);
        const double &ResProj_0_2 = data.ResProj(0, 2);
        const double &ResProj_0_3 = data.ResProj(0, 3);
        const double &ResProj_1_0 = data.ResProj(1, 0);
        const double &ResProj_1_1 = data.ResProj(1, 1);
        const double &ResProj_1_2 = data.ResProj(1, 2);
        const double &ResProj_1_3 = data.ResProj(1, 3);
        const double &ResProj_2_0 = data.ResProj(2, 0);
        const double &ResProj_2_1 = data.ResProj(2, 1);
        const double &ResProj_2_2 = data.ResProj(2, 2);
        const double &ResProj_2_3 = data.ResProj(2, 3);

        const double crRightHandSideBoundedVector0 =             0.16666666666666666*U_1_0;
const double crRightHandSideBoundedVector1 =             0.16666666666666666*U_2_0;
const double crRightHandSideBoundedVector2 =             0.66666666666666663*U_0_0 + crRightHandSideBoundedVector0 + crRightHandSideBoundedVector1;
const double crRightHandSideBoundedVector3 =             1.0/crRightHandSideBoundedVector2;
const double crRightHandSideBoundedVector4 =             stab_c1/pow(h, 2);
const double crRightHandSideBoundedVector5 =             crRightHandSideBoundedVector3*crRightHandSideBoundedVector4;
const double crRightHandSideBoundedVector6 =             1.3333333333333333*mu;
const double crRightHandSideBoundedVector7 =             0.25*U_1_0;
const double crRightHandSideBoundedVector8 =             0.25*U_2_0;
const double crRightHandSideBoundedVector9 =             U_0_0 + crRightHandSideBoundedVector7 + crRightHandSideBoundedVector8;
const double crRightHandSideBoundedVector10 =             pow(crRightHandSideBoundedVector9, -2);
const double crRightHandSideBoundedVector11 =             0.25*U_1_1;
const double crRightHandSideBoundedVector12 =             0.25*U_2_1;
const double crRightHandSideBoundedVector13 =             pow(U_0_1 + crRightHandSideBoundedVector11 + crRightHandSideBoundedVector12, 2);
const double crRightHandSideBoundedVector14 =             0.25*U_1_2;
const double crRightHandSideBoundedVector15 =             0.25*U_2_2;
const double crRightHandSideBoundedVector16 =             pow(U_0_2 + crRightHandSideBoundedVector14 + crRightHandSideBoundedVector15, 2);
const double crRightHandSideBoundedVector17 =             crRightHandSideBoundedVector13 + crRightHandSideBoundedVector16;
const double crRightHandSideBoundedVector18 =             sqrt(gamma);
const double crRightHandSideBoundedVector19 =             gamma - 1;
const double crRightHandSideBoundedVector20 =             0.5*crRightHandSideBoundedVector10;
const double crRightHandSideBoundedVector21 =             0.16666666666666666*U_1_3;
const double crRightHandSideBoundedVector22 =             0.16666666666666666*U_2_3;
const double crRightHandSideBoundedVector23 =             0.66666666666666663*U_0_3 + crRightHandSideBoundedVector21 + crRightHandSideBoundedVector22;
const double crRightHandSideBoundedVector24 =             crRightHandSideBoundedVector13*crRightHandSideBoundedVector20 + crRightHandSideBoundedVector16*crRightHandSideBoundedVector20 - crRightHandSideBoundedVector23*crRightHandSideBoundedVector3;
const double crRightHandSideBoundedVector25 =             crRightHandSideBoundedVector19*crRightHandSideBoundedVector24;
const double crRightHandSideBoundedVector26 =             stab_c2/h;
const double crRightHandSideBoundedVector27 =             pow(crRightHandSideBoundedVector19, -2);
const double crRightHandSideBoundedVector28 =             0.25*r_ext[1];
const double crRightHandSideBoundedVector29 =             0.25*r_ext[2];
const double crRightHandSideBoundedVector30 =             pow(crRightHandSideBoundedVector28 + crRightHandSideBoundedVector29 + r_ext[0], 2);
const double crRightHandSideBoundedVector31 =             0.25*f_ext(1,0);
const double crRightHandSideBoundedVector32 =             0.25*f_ext(2,0);
const double crRightHandSideBoundedVector33 =             0.25*f_ext(1,1);
const double crRightHandSideBoundedVector34 =             0.25*f_ext(2,1);
const double crRightHandSideBoundedVector35 =             crRightHandSideBoundedVector25*(pow(crRightHandSideBoundedVector31 + crRightHandSideBoundedVector32 + f_ext(0,0), 2) + pow(crRightHandSideBoundedVector33 + crRightHandSideBoundedVector34 + f_ext(0,1), 2));
const double crRightHandSideBoundedVector36 =             0.88888888888888884*gamma;
const double crRightHandSideBoundedVector37 =             0.44444444444444442*gamma;
const double crRightHandSideBoundedVector38 =             1.0/gamma;
const double crRightHandSideBoundedVector39 =             0.70710678118654757*crRightHandSideBoundedVector38*stab_c3;
const double crRightHandSideBoundedVector40 =             crRightHandSideBoundedVector26*(crRightHandSideBoundedVector18*sqrt(-crRightHandSideBoundedVector25) + 1.0*sqrt(crRightHandSideBoundedVector10*crRightHandSideBoundedVector17)) + crRightHandSideBoundedVector39*sqrt(crRightHandSideBoundedVector27*(0.44444444444444442*crRightHandSideBoundedVector30 - crRightHandSideBoundedVector35*crRightHandSideBoundedVector36 + 1.3333333333333333*sqrt(crRightHandSideBoundedVector30*(0.11111111111111109*crRightHandSideBoundedVector30 - crRightHandSideBoundedVector35*crRightHandSideBoundedVector37)))/pow(crRightHandSideBoundedVector24, 2));
const double crRightHandSideBoundedVector41 =             1.0/(crRightHandSideBoundedVector40 + crRightHandSideBoundedVector5*crRightHandSideBoundedVector6);
const double crRightHandSideBoundedVector42 =             0.16666666666666666*f_ext(1,0);
const double crRightHandSideBoundedVector43 =             0.16666666666666666*f_ext(2,0);
const double crRightHandSideBoundedVector44 =             crRightHandSideBoundedVector42 + crRightHandSideBoundedVector43 + 0.66666666666666663*f_ext(0,0);
const double crRightHandSideBoundedVector45 =             crRightHandSideBoundedVector2*crRightHandSideBoundedVector44;
const double crRightHandSideBoundedVector46 =             DN_DX_0_0*U_0_0 + DN_DX_1_0*U_1_0 + DN_DX_2_0*U_2_0;
const double crRightHandSideBoundedVector47 =             1.0000000000000002*crRightHandSideBoundedVector13;
const double crRightHandSideBoundedVector48 =             0.50000000000000011*crRightHandSideBoundedVector19;
const double crRightHandSideBoundedVector49 =             crRightHandSideBoundedVector17*crRightHandSideBoundedVector48;
const double crRightHandSideBoundedVector50 =             crRightHandSideBoundedVector10*(-crRightHandSideBoundedVector47 + crRightHandSideBoundedVector49);
const double crRightHandSideBoundedVector51 =             crRightHandSideBoundedVector46*crRightHandSideBoundedVector50;
const double crRightHandSideBoundedVector52 =             DN_DX_0_1*U_0_2;
const double crRightHandSideBoundedVector53 =             DN_DX_1_1*U_1_2;
const double crRightHandSideBoundedVector54 =             DN_DX_2_1*U_2_2;
const double crRightHandSideBoundedVector55 =             crRightHandSideBoundedVector52 + crRightHandSideBoundedVector53 + crRightHandSideBoundedVector54;
const double crRightHandSideBoundedVector56 =             0.66666666666666663*U_0_1;
const double crRightHandSideBoundedVector57 =             0.16666666666666666*U_1_1;
const double crRightHandSideBoundedVector58 =             0.16666666666666666*U_2_1;
const double crRightHandSideBoundedVector59 =             crRightHandSideBoundedVector56 + crRightHandSideBoundedVector57 + crRightHandSideBoundedVector58;
const double crRightHandSideBoundedVector60 =             crRightHandSideBoundedVector3*crRightHandSideBoundedVector59;
const double crRightHandSideBoundedVector61 =             crRightHandSideBoundedVector55*crRightHandSideBoundedVector60;
const double crRightHandSideBoundedVector62 =             DN_DX_0_1*U_0_1;
const double crRightHandSideBoundedVector63 =             DN_DX_1_1*U_1_1;
const double crRightHandSideBoundedVector64 =             DN_DX_2_1*U_2_1;
const double crRightHandSideBoundedVector65 =             crRightHandSideBoundedVector62 + crRightHandSideBoundedVector63 + crRightHandSideBoundedVector64;
const double crRightHandSideBoundedVector66 =             0.66666666666666663*U_0_2;
const double crRightHandSideBoundedVector67 =             0.16666666666666666*U_1_2;
const double crRightHandSideBoundedVector68 =             0.16666666666666666*U_2_2;
const double crRightHandSideBoundedVector69 =             crRightHandSideBoundedVector66 + crRightHandSideBoundedVector67 + crRightHandSideBoundedVector68;
const double crRightHandSideBoundedVector70 =             crRightHandSideBoundedVector3*crRightHandSideBoundedVector69;
const double crRightHandSideBoundedVector71 =             crRightHandSideBoundedVector65*crRightHandSideBoundedVector70;
const double crRightHandSideBoundedVector72 =             DN_DX_0_0*U_0_1;
const double crRightHandSideBoundedVector73 =             DN_DX_1_0*U_1_1;
const double crRightHandSideBoundedVector74 =             DN_DX_2_0*U_2_1;
const double crRightHandSideBoundedVector75 =             crRightHandSideBoundedVector72 + crRightHandSideBoundedVector73 + crRightHandSideBoundedVector74;
const double crRightHandSideBoundedVector76 =             1.0*gamma;
const double crRightHandSideBoundedVector77 =             3.0 - crRightHandSideBoundedVector76;
const double crRightHandSideBoundedVector78 =             crRightHandSideBoundedVector75*crRightHandSideBoundedVector77;
const double crRightHandSideBoundedVector79 =             DN_DX_0_0*U_0_2;
const double crRightHandSideBoundedVector80 =             DN_DX_1_0*U_1_2;
const double crRightHandSideBoundedVector81 =             DN_DX_2_0*U_2_2;
const double crRightHandSideBoundedVector82 =             crRightHandSideBoundedVector79 + crRightHandSideBoundedVector80 + crRightHandSideBoundedVector81;
const double crRightHandSideBoundedVector83 =             crRightHandSideBoundedVector19*crRightHandSideBoundedVector82;
const double crRightHandSideBoundedVector84 =             1.0*crRightHandSideBoundedVector83;
const double crRightHandSideBoundedVector85 =             DN_DX_0_1*U_0_0 + DN_DX_1_1*U_1_0 + DN_DX_2_1*U_2_0;
const double crRightHandSideBoundedVector86 =             crRightHandSideBoundedVector69*crRightHandSideBoundedVector85;
const double crRightHandSideBoundedVector87 =             2.2500000000000004*crRightHandSideBoundedVector59;
const double crRightHandSideBoundedVector88 =             crRightHandSideBoundedVector10*crRightHandSideBoundedVector87;
const double crRightHandSideBoundedVector89 =             DN_DX_0_0*U_0_3 + DN_DX_1_0*U_1_3 + DN_DX_2_0*U_2_3;
const double crRightHandSideBoundedVector90 =             crRightHandSideBoundedVector19*crRightHandSideBoundedVector89;
const double crRightHandSideBoundedVector91 =             0.16666666666666666*ResProj_2_1 + crRightHandSideBoundedVector90 + 0.16666666666666666*dUdt_2_1;
const double crRightHandSideBoundedVector92 =             0.16666666666666666*ResProj_1_1 + 0.16666666666666666*dUdt_1_1;
const double crRightHandSideBoundedVector93 =             crRightHandSideBoundedVector41*(0.66666666666666663*ResProj_0_1 - crRightHandSideBoundedVector45 + crRightHandSideBoundedVector51 + crRightHandSideBoundedVector60*crRightHandSideBoundedVector78 + crRightHandSideBoundedVector61 - crRightHandSideBoundedVector70*crRightHandSideBoundedVector84 + crRightHandSideBoundedVector71 - crRightHandSideBoundedVector86*crRightHandSideBoundedVector88 + crRightHandSideBoundedVector91 + crRightHandSideBoundedVector92 + 0.66666666666666663*dUdt_0_1);
const double crRightHandSideBoundedVector94 =             0.16666666666666666*U_0_0;
const double crRightHandSideBoundedVector95 =             0.66666666666666663*U_1_0 + crRightHandSideBoundedVector1 + crRightHandSideBoundedVector94;
const double crRightHandSideBoundedVector96 =             1.0/crRightHandSideBoundedVector95;
const double crRightHandSideBoundedVector97 =             crRightHandSideBoundedVector4*crRightHandSideBoundedVector6;
const double crRightHandSideBoundedVector98 =             0.25*U_0_0;
const double crRightHandSideBoundedVector99 =             U_1_0 + crRightHandSideBoundedVector8 + crRightHandSideBoundedVector98;
const double crRightHandSideBoundedVector100 =             pow(crRightHandSideBoundedVector99, -2);
const double crRightHandSideBoundedVector101 =             0.25*U_0_1;
const double crRightHandSideBoundedVector102 =             pow(U_1_1 + crRightHandSideBoundedVector101 + crRightHandSideBoundedVector12, 2);
const double crRightHandSideBoundedVector103 =             0.25*U_0_2;
const double crRightHandSideBoundedVector104 =             pow(U_1_2 + crRightHandSideBoundedVector103 + crRightHandSideBoundedVector15, 2);
const double crRightHandSideBoundedVector105 =             crRightHandSideBoundedVector102 + crRightHandSideBoundedVector104;
const double crRightHandSideBoundedVector106 =             0.5*crRightHandSideBoundedVector100;
const double crRightHandSideBoundedVector107 =             0.16666666666666666*U_0_3;
const double crRightHandSideBoundedVector108 =             0.66666666666666663*U_1_3 + crRightHandSideBoundedVector107 + crRightHandSideBoundedVector22;
const double crRightHandSideBoundedVector109 =             crRightHandSideBoundedVector102*crRightHandSideBoundedVector106 + crRightHandSideBoundedVector104*crRightHandSideBoundedVector106 - crRightHandSideBoundedVector108*crRightHandSideBoundedVector96;
const double crRightHandSideBoundedVector110 =             crRightHandSideBoundedVector109*crRightHandSideBoundedVector19;
const double crRightHandSideBoundedVector111 =             0.25*r_ext[0];
const double crRightHandSideBoundedVector112 =             pow(crRightHandSideBoundedVector111 + crRightHandSideBoundedVector29 + r_ext[1], 2);
const double crRightHandSideBoundedVector113 =             0.25*f_ext(0,0);
const double crRightHandSideBoundedVector114 =             0.25*f_ext(0,1);
const double crRightHandSideBoundedVector115 =             crRightHandSideBoundedVector110*(pow(crRightHandSideBoundedVector113 + crRightHandSideBoundedVector32 + f_ext(1,0), 2) + pow(crRightHandSideBoundedVector114 + crRightHandSideBoundedVector34 + f_ext(1,1), 2));
const double crRightHandSideBoundedVector116 =             crRightHandSideBoundedVector26*(crRightHandSideBoundedVector18*sqrt(-crRightHandSideBoundedVector110) + 1.0*sqrt(crRightHandSideBoundedVector100*crRightHandSideBoundedVector105)) + crRightHandSideBoundedVector39*sqrt(crRightHandSideBoundedVector27*(0.44444444444444442*crRightHandSideBoundedVector112 - crRightHandSideBoundedVector115*crRightHandSideBoundedVector36 + 1.3333333333333333*sqrt(crRightHandSideBoundedVector112*(0.11111111111111109*crRightHandSideBoundedVector112 - crRightHandSideBoundedVector115*crRightHandSideBoundedVector37)))/pow(crRightHandSideBoundedVector109, 2));
const double crRightHandSideBoundedVector117 =             1.0/(crRightHandSideBoundedVector116 + crRightHandSideBoundedVector96*crRightHandSideBoundedVector97);
const double crRightHandSideBoundedVector118 =             0.16666666666666666*f_ext(0,0);
const double crRightHandSideBoundedVector119 =             crRightHandSideBoundedVector118 + crRightHandSideBoundedVector43 + 0.66666666666666663*f_ext(1,0);
const double crRightHandSideBoundedVector120 =             crRightHandSideBoundedVector119*crRightHandSideBoundedVector95;
const double crRightHandSideBoundedVector121 =             1.0000000000000002*crRightHandSideBoundedVector102;
const double crRightHandSideBoundedVector122 =             crRightHandSideBoundedVector105*crRightHandSideBoundedVector48;
const double crRightHandSideBoundedVector123 =             crRightHandSideBoundedVector100*(-crRightHandSideBoundedVector121 + crRightHandSideBoundedVector122);
const double crRightHandSideBoundedVector124 =             crRightHandSideBoundedVector123*crRightHandSideBoundedVector46;
const double crRightHandSideBoundedVector125 =             0.16666666666666666*U_0_1;
const double crRightHandSideBoundedVector126 =             0.66666666666666663*U_1_1 + crRightHandSideBoundedVector125 + crRightHandSideBoundedVector58;
const double crRightHandSideBoundedVector127 =             crRightHandSideBoundedVector126*crRightHandSideBoundedVector96;
const double crRightHandSideBoundedVector128 =             crRightHandSideBoundedVector127*crRightHandSideBoundedVector55;
const double crRightHandSideBoundedVector129 =             0.16666666666666666*U_0_2;
const double crRightHandSideBoundedVector130 =             0.66666666666666663*U_1_2 + crRightHandSideBoundedVector129 + crRightHandSideBoundedVector68;
const double crRightHandSideBoundedVector131 =             crRightHandSideBoundedVector130*crRightHandSideBoundedVector96;
const double crRightHandSideBoundedVector132 =             crRightHandSideBoundedVector131*crRightHandSideBoundedVector65;
const double crRightHandSideBoundedVector133 =             crRightHandSideBoundedVector130*crRightHandSideBoundedVector85;
const double crRightHandSideBoundedVector134 =             2.2500000000000004*crRightHandSideBoundedVector126;
const double crRightHandSideBoundedVector135 =             crRightHandSideBoundedVector100*crRightHandSideBoundedVector134;
const double crRightHandSideBoundedVector136 =             0.16666666666666666*ResProj_0_1 + 0.16666666666666666*dUdt_0_1;
const double crRightHandSideBoundedVector137 =             crRightHandSideBoundedVector117*(0.66666666666666663*ResProj_1_1 - crRightHandSideBoundedVector120 + crRightHandSideBoundedVector124 + crRightHandSideBoundedVector127*crRightHandSideBoundedVector78 + crRightHandSideBoundedVector128 - crRightHandSideBoundedVector131*crRightHandSideBoundedVector84 + crRightHandSideBoundedVector132 - crRightHandSideBoundedVector133*crRightHandSideBoundedVector135 + crRightHandSideBoundedVector136 + crRightHandSideBoundedVector91 + 0.66666666666666663*dUdt_1_1);
const double crRightHandSideBoundedVector138 =             0.66666666666666663*U_2_0 + crRightHandSideBoundedVector0 + crRightHandSideBoundedVector94;
const double crRightHandSideBoundedVector139 =             1.0/crRightHandSideBoundedVector138;
const double crRightHandSideBoundedVector140 =             U_2_0 + crRightHandSideBoundedVector7 + crRightHandSideBoundedVector98;
const double crRightHandSideBoundedVector141 =             pow(crRightHandSideBoundedVector140, -2);
const double crRightHandSideBoundedVector142 =             pow(U_2_1 + crRightHandSideBoundedVector101 + crRightHandSideBoundedVector11, 2);
const double crRightHandSideBoundedVector143 =             pow(U_2_2 + crRightHandSideBoundedVector103 + crRightHandSideBoundedVector14, 2);
const double crRightHandSideBoundedVector144 =             crRightHandSideBoundedVector142 + crRightHandSideBoundedVector143;
const double crRightHandSideBoundedVector145 =             0.5*crRightHandSideBoundedVector141;
const double crRightHandSideBoundedVector146 =             0.66666666666666663*U_2_3 + crRightHandSideBoundedVector107 + crRightHandSideBoundedVector21;
const double crRightHandSideBoundedVector147 =             -crRightHandSideBoundedVector139*crRightHandSideBoundedVector146 + crRightHandSideBoundedVector142*crRightHandSideBoundedVector145 + crRightHandSideBoundedVector143*crRightHandSideBoundedVector145;
const double crRightHandSideBoundedVector148 =             crRightHandSideBoundedVector147*crRightHandSideBoundedVector19;
const double crRightHandSideBoundedVector149 =             pow(crRightHandSideBoundedVector111 + crRightHandSideBoundedVector28 + r_ext[2], 2);
const double crRightHandSideBoundedVector150 =             crRightHandSideBoundedVector148*(pow(crRightHandSideBoundedVector113 + crRightHandSideBoundedVector31 + f_ext(2,0), 2) + pow(crRightHandSideBoundedVector114 + crRightHandSideBoundedVector33 + f_ext(2,1), 2));
const double crRightHandSideBoundedVector151 =             crRightHandSideBoundedVector26*(crRightHandSideBoundedVector18*sqrt(-crRightHandSideBoundedVector148) + 1.0*sqrt(crRightHandSideBoundedVector141*crRightHandSideBoundedVector144)) + crRightHandSideBoundedVector39*sqrt(crRightHandSideBoundedVector27*(0.44444444444444442*crRightHandSideBoundedVector149 - crRightHandSideBoundedVector150*crRightHandSideBoundedVector36 + 1.3333333333333333*sqrt(crRightHandSideBoundedVector149*(0.11111111111111109*crRightHandSideBoundedVector149 - crRightHandSideBoundedVector150*crRightHandSideBoundedVector37)))/pow(crRightHandSideBoundedVector147, 2));
const double crRightHandSideBoundedVector152 =             1.0/(crRightHandSideBoundedVector139*crRightHandSideBoundedVector97 + crRightHandSideBoundedVector151);
const double crRightHandSideBoundedVector153 =             crRightHandSideBoundedVector118 + crRightHandSideBoundedVector42 + 0.66666666666666663*f_ext(2,0);
const double crRightHandSideBoundedVector154 =             crRightHandSideBoundedVector138*crRightHandSideBoundedVector153;
const double crRightHandSideBoundedVector155 =             1.0000000000000002*crRightHandSideBoundedVector142;
const double crRightHandSideBoundedVector156 =             crRightHandSideBoundedVector144*crRightHandSideBoundedVector48;
const double crRightHandSideBoundedVector157 =             crRightHandSideBoundedVector141*(-crRightHandSideBoundedVector155 + crRightHandSideBoundedVector156);
const double crRightHandSideBoundedVector158 =             crRightHandSideBoundedVector157*crRightHandSideBoundedVector46;
const double crRightHandSideBoundedVector159 =             0.66666666666666663*U_2_1 + crRightHandSideBoundedVector125 + crRightHandSideBoundedVector57;
const double crRightHandSideBoundedVector160 =             crRightHandSideBoundedVector139*crRightHandSideBoundedVector159;
const double crRightHandSideBoundedVector161 =             crRightHandSideBoundedVector160*crRightHandSideBoundedVector55;
const double crRightHandSideBoundedVector162 =             0.66666666666666663*U_2_2 + crRightHandSideBoundedVector129 + crRightHandSideBoundedVector67;
const double crRightHandSideBoundedVector163 =             crRightHandSideBoundedVector139*crRightHandSideBoundedVector162;
const double crRightHandSideBoundedVector164 =             crRightHandSideBoundedVector163*crRightHandSideBoundedVector65;
const double crRightHandSideBoundedVector165 =             crRightHandSideBoundedVector162*crRightHandSideBoundedVector85;
const double crRightHandSideBoundedVector166 =             2.2500000000000004*crRightHandSideBoundedVector159;
const double crRightHandSideBoundedVector167 =             crRightHandSideBoundedVector141*crRightHandSideBoundedVector166;
const double crRightHandSideBoundedVector168 =             crRightHandSideBoundedVector152*(0.66666666666666663*ResProj_2_1 + crRightHandSideBoundedVector136 - crRightHandSideBoundedVector154 + crRightHandSideBoundedVector158 + crRightHandSideBoundedVector160*crRightHandSideBoundedVector78 + crRightHandSideBoundedVector161 - crRightHandSideBoundedVector163*crRightHandSideBoundedVector84 + crRightHandSideBoundedVector164 - crRightHandSideBoundedVector165*crRightHandSideBoundedVector167 + crRightHandSideBoundedVector90 + crRightHandSideBoundedVector92 + 0.66666666666666663*dUdt_2_1);
const double crRightHandSideBoundedVector169 =             0.16666666666666666*f_ext(1,1);
const double crRightHandSideBoundedVector170 =             0.16666666666666666*f_ext(2,1);
const double crRightHandSideBoundedVector171 =             crRightHandSideBoundedVector169 + crRightHandSideBoundedVector170 + 0.66666666666666663*f_ext(0,1);
const double crRightHandSideBoundedVector172 =             crRightHandSideBoundedVector171*crRightHandSideBoundedVector2;
const double crRightHandSideBoundedVector173 =             1.0000000000000002*crRightHandSideBoundedVector16;
const double crRightHandSideBoundedVector174 =             crRightHandSideBoundedVector10*(-crRightHandSideBoundedVector173 + crRightHandSideBoundedVector49);
const double crRightHandSideBoundedVector175 =             crRightHandSideBoundedVector174*crRightHandSideBoundedVector85;
const double crRightHandSideBoundedVector176 =             crRightHandSideBoundedVector60*crRightHandSideBoundedVector82;
const double crRightHandSideBoundedVector177 =             crRightHandSideBoundedVector70*crRightHandSideBoundedVector75;
const double crRightHandSideBoundedVector178 =             crRightHandSideBoundedVector55*crRightHandSideBoundedVector77;
const double crRightHandSideBoundedVector179 =             crRightHandSideBoundedVector19*crRightHandSideBoundedVector65;
const double crRightHandSideBoundedVector180 =             1.0*crRightHandSideBoundedVector179;
const double crRightHandSideBoundedVector181 =             2.2500000000000004*crRightHandSideBoundedVector69;
const double crRightHandSideBoundedVector182 =             crRightHandSideBoundedVector46*crRightHandSideBoundedVector59;
const double crRightHandSideBoundedVector183 =             crRightHandSideBoundedVector10*crRightHandSideBoundedVector182;
const double crRightHandSideBoundedVector184 =             DN_DX_0_1*U_0_3 + DN_DX_1_1*U_1_3 + DN_DX_2_1*U_2_3;
const double crRightHandSideBoundedVector185 =             crRightHandSideBoundedVector184*crRightHandSideBoundedVector19;
const double crRightHandSideBoundedVector186 =             0.16666666666666666*ResProj_2_2 + crRightHandSideBoundedVector185 + 0.16666666666666666*dUdt_2_2;
const double crRightHandSideBoundedVector187 =             0.16666666666666666*ResProj_1_2 + 0.16666666666666666*dUdt_1_2;
const double crRightHandSideBoundedVector188 =             crRightHandSideBoundedVector41*(0.66666666666666663*ResProj_0_2 - crRightHandSideBoundedVector172 + crRightHandSideBoundedVector175 + crRightHandSideBoundedVector176 + crRightHandSideBoundedVector177 + crRightHandSideBoundedVector178*crRightHandSideBoundedVector70 - crRightHandSideBoundedVector180*crRightHandSideBoundedVector60 - crRightHandSideBoundedVector181*crRightHandSideBoundedVector183 + crRightHandSideBoundedVector186 + crRightHandSideBoundedVector187 + 0.66666666666666663*dUdt_0_2);
const double crRightHandSideBoundedVector189 =             0.16666666666666666*f_ext(0,1);
const double crRightHandSideBoundedVector190 =             crRightHandSideBoundedVector170 + crRightHandSideBoundedVector189 + 0.66666666666666663*f_ext(1,1);
const double crRightHandSideBoundedVector191 =             crRightHandSideBoundedVector190*crRightHandSideBoundedVector95;
const double crRightHandSideBoundedVector192 =             1.0000000000000002*crRightHandSideBoundedVector104;
const double crRightHandSideBoundedVector193 =             crRightHandSideBoundedVector100*(crRightHandSideBoundedVector122 - crRightHandSideBoundedVector192);
const double crRightHandSideBoundedVector194 =             crRightHandSideBoundedVector193*crRightHandSideBoundedVector85;
const double crRightHandSideBoundedVector195 =             crRightHandSideBoundedVector127*crRightHandSideBoundedVector82;
const double crRightHandSideBoundedVector196 =             crRightHandSideBoundedVector131*crRightHandSideBoundedVector75;
const double crRightHandSideBoundedVector197 =             2.2500000000000004*crRightHandSideBoundedVector130;
const double crRightHandSideBoundedVector198 =             crRightHandSideBoundedVector126*crRightHandSideBoundedVector46;
const double crRightHandSideBoundedVector199 =             crRightHandSideBoundedVector100*crRightHandSideBoundedVector198;
const double crRightHandSideBoundedVector200 =             0.16666666666666666*ResProj_0_2 + 0.16666666666666666*dUdt_0_2;
const double crRightHandSideBoundedVector201 =             crRightHandSideBoundedVector117*(0.66666666666666663*ResProj_1_2 - crRightHandSideBoundedVector127*crRightHandSideBoundedVector180 + crRightHandSideBoundedVector131*crRightHandSideBoundedVector178 + crRightHandSideBoundedVector186 - crRightHandSideBoundedVector191 + crRightHandSideBoundedVector194 + crRightHandSideBoundedVector195 + crRightHandSideBoundedVector196 - crRightHandSideBoundedVector197*crRightHandSideBoundedVector199 + crRightHandSideBoundedVector200 + 0.66666666666666663*dUdt_1_2);
const double crRightHandSideBoundedVector202 =             crRightHandSideBoundedVector169 + crRightHandSideBoundedVector189 + 0.66666666666666663*f_ext(2,1);
const double crRightHandSideBoundedVector203 =             crRightHandSideBoundedVector138*crRightHandSideBoundedVector202;
const double crRightHandSideBoundedVector204 =             1.0000000000000002*crRightHandSideBoundedVector143;
const double crRightHandSideBoundedVector205 =             crRightHandSideBoundedVector141*(crRightHandSideBoundedVector156 - crRightHandSideBoundedVector204);
const double crRightHandSideBoundedVector206 =             crRightHandSideBoundedVector205*crRightHandSideBoundedVector85;
const double crRightHandSideBoundedVector207 =             crRightHandSideBoundedVector160*crRightHandSideBoundedVector82;
const double crRightHandSideBoundedVector208 =             crRightHandSideBoundedVector163*crRightHandSideBoundedVector75;
const double crRightHandSideBoundedVector209 =             2.2500000000000004*crRightHandSideBoundedVector162;
const double crRightHandSideBoundedVector210 =             crRightHandSideBoundedVector159*crRightHandSideBoundedVector46;
const double crRightHandSideBoundedVector211 =             crRightHandSideBoundedVector141*crRightHandSideBoundedVector210;
const double crRightHandSideBoundedVector212 =             crRightHandSideBoundedVector152*(0.66666666666666663*ResProj_2_2 - crRightHandSideBoundedVector160*crRightHandSideBoundedVector180 + crRightHandSideBoundedVector163*crRightHandSideBoundedVector178 + crRightHandSideBoundedVector185 + crRightHandSideBoundedVector187 + crRightHandSideBoundedVector200 - crRightHandSideBoundedVector203 + crRightHandSideBoundedVector206 + crRightHandSideBoundedVector207 + crRightHandSideBoundedVector208 - crRightHandSideBoundedVector209*crRightHandSideBoundedVector211 + 0.66666666666666663*dUdt_2_2);
const double crRightHandSideBoundedVector213 =             crRightHandSideBoundedVector55 + crRightHandSideBoundedVector75;
const double crRightHandSideBoundedVector214 =             crRightHandSideBoundedVector2*crRightHandSideBoundedVector75;
const double crRightHandSideBoundedVector215 =             -crRightHandSideBoundedVector182 + crRightHandSideBoundedVector214;
const double crRightHandSideBoundedVector216 =             crRightHandSideBoundedVector10*mu;
const double crRightHandSideBoundedVector217 =             4.5000000000000009*crRightHandSideBoundedVector216;
const double crRightHandSideBoundedVector218 =             crRightHandSideBoundedVector2*crRightHandSideBoundedVector55;
const double crRightHandSideBoundedVector219 =             crRightHandSideBoundedVector218 - crRightHandSideBoundedVector86;
const double crRightHandSideBoundedVector220 =             -1.5000000000000002*crRightHandSideBoundedVector216*(crRightHandSideBoundedVector215 + crRightHandSideBoundedVector219);
const double crRightHandSideBoundedVector221 =             0.66666666666666663*crRightHandSideBoundedVector182;
const double crRightHandSideBoundedVector222 =             crRightHandSideBoundedVector3*(-0.66666666666666663*crRightHandSideBoundedVector214 + 1.3333333333333335*crRightHandSideBoundedVector218 + crRightHandSideBoundedVector221 - 1.3333333333333335*crRightHandSideBoundedVector86);
const double crRightHandSideBoundedVector223 =             1.0/lin_m_norm;
const double crRightHandSideBoundedVector224 =             crRightHandSideBoundedVector223*lin_m[0]*lin_m[1];
const double crRightHandSideBoundedVector225 =             crRightHandSideBoundedVector224*nu_st;
const double crRightHandSideBoundedVector226 =             0.66666666666666663*crRightHandSideBoundedVector86;
const double crRightHandSideBoundedVector227 =             crRightHandSideBoundedVector3*(-1.3333333333333335*crRightHandSideBoundedVector182 + 1.3333333333333335*crRightHandSideBoundedVector214 - 0.66666666666666663*crRightHandSideBoundedVector218 + crRightHandSideBoundedVector226);
const double crRightHandSideBoundedVector228 =             crRightHandSideBoundedVector223*pow(lin_m[0], 2);
const double crRightHandSideBoundedVector229 =             crRightHandSideBoundedVector228*nu_st;
const double crRightHandSideBoundedVector230 =             crRightHandSideBoundedVector224*nu_sc;
const double crRightHandSideBoundedVector231 =             1 - crRightHandSideBoundedVector228;
const double crRightHandSideBoundedVector232 =             crRightHandSideBoundedVector231*nu_sc;
const double crRightHandSideBoundedVector233 =             crRightHandSideBoundedVector215*crRightHandSideBoundedVector217 + crRightHandSideBoundedVector220 + crRightHandSideBoundedVector222*crRightHandSideBoundedVector225 - crRightHandSideBoundedVector222*crRightHandSideBoundedVector230 + crRightHandSideBoundedVector227*crRightHandSideBoundedVector229 + crRightHandSideBoundedVector227*crRightHandSideBoundedVector232;
const double crRightHandSideBoundedVector234 =             crRightHandSideBoundedVector75*crRightHandSideBoundedVector95;
const double crRightHandSideBoundedVector235 =             -crRightHandSideBoundedVector198 + crRightHandSideBoundedVector234;
const double crRightHandSideBoundedVector236 =             crRightHandSideBoundedVector100*mu;
const double crRightHandSideBoundedVector237 =             4.5000000000000009*crRightHandSideBoundedVector236;
const double crRightHandSideBoundedVector238 =             crRightHandSideBoundedVector55*crRightHandSideBoundedVector95;
const double crRightHandSideBoundedVector239 =             -crRightHandSideBoundedVector133 + crRightHandSideBoundedVector238;
const double crRightHandSideBoundedVector240 =             -1.5000000000000002*crRightHandSideBoundedVector236*(crRightHandSideBoundedVector235 + crRightHandSideBoundedVector239);
const double crRightHandSideBoundedVector241 =             0.66666666666666663*crRightHandSideBoundedVector198;
const double crRightHandSideBoundedVector242 =             crRightHandSideBoundedVector96*(-1.3333333333333335*crRightHandSideBoundedVector133 - 0.66666666666666663*crRightHandSideBoundedVector234 + 1.3333333333333335*crRightHandSideBoundedVector238 + crRightHandSideBoundedVector241);
const double crRightHandSideBoundedVector243 =             0.66666666666666663*crRightHandSideBoundedVector133;
const double crRightHandSideBoundedVector244 =             crRightHandSideBoundedVector96*(-1.3333333333333335*crRightHandSideBoundedVector198 + 1.3333333333333335*crRightHandSideBoundedVector234 - 0.66666666666666663*crRightHandSideBoundedVector238 + crRightHandSideBoundedVector243);
const double crRightHandSideBoundedVector245 =             crRightHandSideBoundedVector225*crRightHandSideBoundedVector242 + crRightHandSideBoundedVector229*crRightHandSideBoundedVector244 - crRightHandSideBoundedVector230*crRightHandSideBoundedVector242 + crRightHandSideBoundedVector232*crRightHandSideBoundedVector244 + crRightHandSideBoundedVector235*crRightHandSideBoundedVector237 + crRightHandSideBoundedVector240;
const double crRightHandSideBoundedVector246 =             crRightHandSideBoundedVector138*crRightHandSideBoundedVector75;
const double crRightHandSideBoundedVector247 =             -crRightHandSideBoundedVector210 + crRightHandSideBoundedVector246;
const double crRightHandSideBoundedVector248 =             crRightHandSideBoundedVector141*mu;
const double crRightHandSideBoundedVector249 =             4.5000000000000009*crRightHandSideBoundedVector248;
const double crRightHandSideBoundedVector250 =             crRightHandSideBoundedVector138*crRightHandSideBoundedVector55;
const double crRightHandSideBoundedVector251 =             -crRightHandSideBoundedVector165 + crRightHandSideBoundedVector250;
const double crRightHandSideBoundedVector252 =             -1.5000000000000002*crRightHandSideBoundedVector248*(crRightHandSideBoundedVector247 + crRightHandSideBoundedVector251);
const double crRightHandSideBoundedVector253 =             0.66666666666666663*crRightHandSideBoundedVector210;
const double crRightHandSideBoundedVector254 =             crRightHandSideBoundedVector139*(-1.3333333333333335*crRightHandSideBoundedVector165 - 0.66666666666666663*crRightHandSideBoundedVector246 + 1.3333333333333335*crRightHandSideBoundedVector250 + crRightHandSideBoundedVector253);
const double crRightHandSideBoundedVector255 =             0.66666666666666663*crRightHandSideBoundedVector165;
const double crRightHandSideBoundedVector256 =             crRightHandSideBoundedVector139*(-1.3333333333333335*crRightHandSideBoundedVector210 + 1.3333333333333335*crRightHandSideBoundedVector246 - 0.66666666666666663*crRightHandSideBoundedVector250 + crRightHandSideBoundedVector255);
const double crRightHandSideBoundedVector257 =             crRightHandSideBoundedVector225*crRightHandSideBoundedVector254 + crRightHandSideBoundedVector229*crRightHandSideBoundedVector256 - crRightHandSideBoundedVector230*crRightHandSideBoundedVector254 + crRightHandSideBoundedVector232*crRightHandSideBoundedVector256 + crRightHandSideBoundedVector247*crRightHandSideBoundedVector249 + crRightHandSideBoundedVector252;
const double crRightHandSideBoundedVector258 =             DN_DX_0_1*crRightHandSideBoundedVector60;
const double crRightHandSideBoundedVector259 =             DN_DX_0_1*crRightHandSideBoundedVector56 + 0.66666666666666663*crRightHandSideBoundedVector63 + 0.66666666666666663*crRightHandSideBoundedVector64;
const double crRightHandSideBoundedVector260 =             0.66666666666666663*crRightHandSideBoundedVector3;
const double crRightHandSideBoundedVector261 =             DN_DX_0_0*crRightHandSideBoundedVector70;
const double crRightHandSideBoundedVector262 =             1.0*crRightHandSideBoundedVector19;
const double crRightHandSideBoundedVector263 =             crRightHandSideBoundedVector59*crRightHandSideBoundedVector85;
const double crRightHandSideBoundedVector264 =             1.5000000000000002*crRightHandSideBoundedVector10;
const double crRightHandSideBoundedVector265 =             crRightHandSideBoundedVector46*crRightHandSideBoundedVector69;
const double crRightHandSideBoundedVector266 =             1.5000000000000002*crRightHandSideBoundedVector19;
const double crRightHandSideBoundedVector267 =             crRightHandSideBoundedVector10*crRightHandSideBoundedVector266;
const double crRightHandSideBoundedVector268 =             1.0*crRightHandSideBoundedVector188;
const double crRightHandSideBoundedVector269 =             DN_DX_0_1*crRightHandSideBoundedVector127;
const double crRightHandSideBoundedVector270 =             DN_DX_0_0*crRightHandSideBoundedVector131;
const double crRightHandSideBoundedVector271 =             DN_DX_1_1*crRightHandSideBoundedVector57 + DN_DX_2_1*crRightHandSideBoundedVector58 + 0.16666666666666666*crRightHandSideBoundedVector62;
const double crRightHandSideBoundedVector272 =             0.16666666666666666*crRightHandSideBoundedVector83;
const double crRightHandSideBoundedVector273 =             crRightHandSideBoundedVector126*crRightHandSideBoundedVector85;
const double crRightHandSideBoundedVector274 =             0.37500000000000006*crRightHandSideBoundedVector100;
const double crRightHandSideBoundedVector275 =             crRightHandSideBoundedVector273*crRightHandSideBoundedVector274;
const double crRightHandSideBoundedVector276 =             crRightHandSideBoundedVector130*crRightHandSideBoundedVector46;
const double crRightHandSideBoundedVector277 =             crRightHandSideBoundedVector274*crRightHandSideBoundedVector276;
const double crRightHandSideBoundedVector278 =             -crRightHandSideBoundedVector19*crRightHandSideBoundedVector277 - crRightHandSideBoundedVector271*crRightHandSideBoundedVector96 + crRightHandSideBoundedVector272*crRightHandSideBoundedVector96 + crRightHandSideBoundedVector275;
const double crRightHandSideBoundedVector279 =             1.0*crRightHandSideBoundedVector201;
const double crRightHandSideBoundedVector280 =             DN_DX_0_1*crRightHandSideBoundedVector160;
const double crRightHandSideBoundedVector281 =             DN_DX_0_0*crRightHandSideBoundedVector163;
const double crRightHandSideBoundedVector282 =             crRightHandSideBoundedVector159*crRightHandSideBoundedVector85;
const double crRightHandSideBoundedVector283 =             0.37500000000000006*crRightHandSideBoundedVector141;
const double crRightHandSideBoundedVector284 =             crRightHandSideBoundedVector282*crRightHandSideBoundedVector283;
const double crRightHandSideBoundedVector285 =             crRightHandSideBoundedVector162*crRightHandSideBoundedVector46;
const double crRightHandSideBoundedVector286 =             crRightHandSideBoundedVector283*crRightHandSideBoundedVector285;
const double crRightHandSideBoundedVector287 =             -crRightHandSideBoundedVector139*crRightHandSideBoundedVector271 + crRightHandSideBoundedVector139*crRightHandSideBoundedVector272 - crRightHandSideBoundedVector19*crRightHandSideBoundedVector286 + crRightHandSideBoundedVector284;
const double crRightHandSideBoundedVector288 =             1.0*crRightHandSideBoundedVector212;
const double crRightHandSideBoundedVector289 =             pow(crRightHandSideBoundedVector9, -3);
const double crRightHandSideBoundedVector290 =             0.66666666666666663*crRightHandSideBoundedVector289;
const double crRightHandSideBoundedVector291 =             3.0000000000000004*crRightHandSideBoundedVector13;
const double crRightHandSideBoundedVector292 =             crRightHandSideBoundedVector17*crRightHandSideBoundedVector266;
const double crRightHandSideBoundedVector293 =             crRightHandSideBoundedVector46*(-crRightHandSideBoundedVector291 + crRightHandSideBoundedVector292);
const double crRightHandSideBoundedVector294 =             crRightHandSideBoundedVector264*crRightHandSideBoundedVector59;
const double crRightHandSideBoundedVector295 =             crRightHandSideBoundedVector264*crRightHandSideBoundedVector69;
const double crRightHandSideBoundedVector296 =             4.5*crRightHandSideBoundedVector289;
const double crRightHandSideBoundedVector297 =             crRightHandSideBoundedVector59*crRightHandSideBoundedVector86;
const double crRightHandSideBoundedVector298 =             crRightHandSideBoundedVector296*crRightHandSideBoundedVector297;
const double crRightHandSideBoundedVector299 =             0.66666666666666663*crRightHandSideBoundedVector75;
const double crRightHandSideBoundedVector300 =             2.2500000000000004*gamma - 6.7500000000000018;
const double crRightHandSideBoundedVector301 =             crRightHandSideBoundedVector10*crRightHandSideBoundedVector300;
const double crRightHandSideBoundedVector302 =             crRightHandSideBoundedVector301*crRightHandSideBoundedVector59;
const double crRightHandSideBoundedVector303 =             crRightHandSideBoundedVector267*crRightHandSideBoundedVector69;
const double crRightHandSideBoundedVector304 =             crRightHandSideBoundedVector303*crRightHandSideBoundedVector82;
const double crRightHandSideBoundedVector305 =             DN_DX_0_1*crRightHandSideBoundedVector10;
const double crRightHandSideBoundedVector306 =             crRightHandSideBoundedVector305*crRightHandSideBoundedVector69;
const double crRightHandSideBoundedVector307 =             crRightHandSideBoundedVector306*crRightHandSideBoundedVector87;
const double crRightHandSideBoundedVector308 =             0.1111111111111111*f_ext(1,0);
const double crRightHandSideBoundedVector309 =             0.1111111111111111*f_ext(2,0);
const double crRightHandSideBoundedVector310 =             crRightHandSideBoundedVector308 + crRightHandSideBoundedVector309 + 0.44444444444444442*f_ext(0,0);
const double crRightHandSideBoundedVector311 =             0.16666666666666666*ResProj_2_0 + crRightHandSideBoundedVector213 + 0.16666666666666666*dUdt_2_0;
const double crRightHandSideBoundedVector312 =             0.16666666666666666*ResProj_1_0 + 0.16666666666666666*dUdt_1_0;
const double crRightHandSideBoundedVector313 =             1.0*(0.66666666666666663*ResProj_0_0 + crRightHandSideBoundedVector311 + crRightHandSideBoundedVector312 + 0.66666666666666663*dUdt_0_0)/crRightHandSideBoundedVector40;
const double crRightHandSideBoundedVector314 =             DN_DX_0_1*crRightHandSideBoundedVector100;
const double crRightHandSideBoundedVector315 =             crRightHandSideBoundedVector130*crRightHandSideBoundedVector314;
const double crRightHandSideBoundedVector316 =             crRightHandSideBoundedVector134*crRightHandSideBoundedVector315;
const double crRightHandSideBoundedVector317 =             0.16666666666666666*crRightHandSideBoundedVector46;
const double crRightHandSideBoundedVector318 =             pow(crRightHandSideBoundedVector99, -3);
const double crRightHandSideBoundedVector319 =             3.0000000000000004*crRightHandSideBoundedVector102;
const double crRightHandSideBoundedVector320 =             crRightHandSideBoundedVector105*crRightHandSideBoundedVector266;
const double crRightHandSideBoundedVector321 =             crRightHandSideBoundedVector318*(-crRightHandSideBoundedVector319 + crRightHandSideBoundedVector320);
const double crRightHandSideBoundedVector322 =             crRightHandSideBoundedVector126*crRightHandSideBoundedVector274;
const double crRightHandSideBoundedVector323 =             crRightHandSideBoundedVector322*crRightHandSideBoundedVector55;
const double crRightHandSideBoundedVector324 =             crRightHandSideBoundedVector130*crRightHandSideBoundedVector274;
const double crRightHandSideBoundedVector325 =             crRightHandSideBoundedVector324*crRightHandSideBoundedVector65;
const double crRightHandSideBoundedVector326 =             1.125*crRightHandSideBoundedVector318;
const double crRightHandSideBoundedVector327 =             crRightHandSideBoundedVector126*crRightHandSideBoundedVector133;
const double crRightHandSideBoundedVector328 =             crRightHandSideBoundedVector326*crRightHandSideBoundedVector327;
const double crRightHandSideBoundedVector329 =             0.16666666666666666*crRightHandSideBoundedVector75;
const double crRightHandSideBoundedVector330 =             crRightHandSideBoundedVector300*crRightHandSideBoundedVector329;
const double crRightHandSideBoundedVector331 =             crRightHandSideBoundedVector100*crRightHandSideBoundedVector126;
const double crRightHandSideBoundedVector332 =             crRightHandSideBoundedVector324*crRightHandSideBoundedVector83;
const double crRightHandSideBoundedVector333 =             0.027777777777777776*f_ext(0,0);
const double crRightHandSideBoundedVector334 =             0.027777777777777776*f_ext(2,0);
const double crRightHandSideBoundedVector335 =             crRightHandSideBoundedVector308 + crRightHandSideBoundedVector333 + crRightHandSideBoundedVector334;
const double crRightHandSideBoundedVector336 =             -crRightHandSideBoundedVector317*crRightHandSideBoundedVector321 - crRightHandSideBoundedVector323 - crRightHandSideBoundedVector325 + crRightHandSideBoundedVector328 + crRightHandSideBoundedVector330*crRightHandSideBoundedVector331 + crRightHandSideBoundedVector332 + crRightHandSideBoundedVector335;
const double crRightHandSideBoundedVector337 =             0.16666666666666666*ResProj_0_0 + 0.16666666666666666*dUdt_0_0;
const double crRightHandSideBoundedVector338 =             1.0*(0.66666666666666663*ResProj_1_0 + crRightHandSideBoundedVector311 + crRightHandSideBoundedVector337 + 0.66666666666666663*dUdt_1_0)/crRightHandSideBoundedVector116;
const double crRightHandSideBoundedVector339 =             DN_DX_0_1*crRightHandSideBoundedVector141;
const double crRightHandSideBoundedVector340 =             crRightHandSideBoundedVector162*crRightHandSideBoundedVector339;
const double crRightHandSideBoundedVector341 =             crRightHandSideBoundedVector166*crRightHandSideBoundedVector340;
const double crRightHandSideBoundedVector342 =             pow(crRightHandSideBoundedVector140, -3);
const double crRightHandSideBoundedVector343 =             3.0000000000000004*crRightHandSideBoundedVector142;
const double crRightHandSideBoundedVector344 =             crRightHandSideBoundedVector144*crRightHandSideBoundedVector266;
const double crRightHandSideBoundedVector345 =             crRightHandSideBoundedVector342*(-crRightHandSideBoundedVector343 + crRightHandSideBoundedVector344);
const double crRightHandSideBoundedVector346 =             crRightHandSideBoundedVector159*crRightHandSideBoundedVector283;
const double crRightHandSideBoundedVector347 =             crRightHandSideBoundedVector346*crRightHandSideBoundedVector55;
const double crRightHandSideBoundedVector348 =             crRightHandSideBoundedVector162*crRightHandSideBoundedVector283;
const double crRightHandSideBoundedVector349 =             crRightHandSideBoundedVector348*crRightHandSideBoundedVector65;
const double crRightHandSideBoundedVector350 =             1.125*crRightHandSideBoundedVector342;
const double crRightHandSideBoundedVector351 =             crRightHandSideBoundedVector159*crRightHandSideBoundedVector165;
const double crRightHandSideBoundedVector352 =             crRightHandSideBoundedVector350*crRightHandSideBoundedVector351;
const double crRightHandSideBoundedVector353 =             crRightHandSideBoundedVector141*crRightHandSideBoundedVector159;
const double crRightHandSideBoundedVector354 =             crRightHandSideBoundedVector348*crRightHandSideBoundedVector83;
const double crRightHandSideBoundedVector355 =             0.027777777777777776*f_ext(1,0);
const double crRightHandSideBoundedVector356 =             crRightHandSideBoundedVector309 + crRightHandSideBoundedVector333 + crRightHandSideBoundedVector355;
const double crRightHandSideBoundedVector357 =             -crRightHandSideBoundedVector317*crRightHandSideBoundedVector345 + crRightHandSideBoundedVector330*crRightHandSideBoundedVector353 - crRightHandSideBoundedVector347 - crRightHandSideBoundedVector349 + crRightHandSideBoundedVector352 + crRightHandSideBoundedVector354 + crRightHandSideBoundedVector356;
const double crRightHandSideBoundedVector358 =             1.0*(0.66666666666666663*ResProj_2_0 + crRightHandSideBoundedVector213 + crRightHandSideBoundedVector312 + crRightHandSideBoundedVector337 + 0.66666666666666663*dUdt_2_0)/crRightHandSideBoundedVector151;
const double crRightHandSideBoundedVector359 =             crRightHandSideBoundedVector76 - 3.0;
const double crRightHandSideBoundedVector360 =             crRightHandSideBoundedVector3*crRightHandSideBoundedVector359;
const double crRightHandSideBoundedVector361 =             DN_DX_0_0*crRightHandSideBoundedVector60;
const double crRightHandSideBoundedVector362 =             0.66666666666666663*crRightHandSideBoundedVector52 + 0.66666666666666663*crRightHandSideBoundedVector53 + 0.66666666666666663*crRightHandSideBoundedVector54;
const double crRightHandSideBoundedVector363 =             DN_DX_0_1*crRightHandSideBoundedVector70;
const double crRightHandSideBoundedVector364 =             1.5000000000000002*crRightHandSideBoundedVector86;
const double crRightHandSideBoundedVector365 =             crRightHandSideBoundedVector10*crRightHandSideBoundedVector364;
const double crRightHandSideBoundedVector366 =             crRightHandSideBoundedVector3*crRightHandSideBoundedVector362 + crRightHandSideBoundedVector363 - crRightHandSideBoundedVector365;
const double crRightHandSideBoundedVector367 =             1.0*crRightHandSideBoundedVector93;
const double crRightHandSideBoundedVector368 =             DN_DX_0_0*crRightHandSideBoundedVector127;
const double crRightHandSideBoundedVector369 =             DN_DX_0_1*crRightHandSideBoundedVector131;
const double crRightHandSideBoundedVector370 =             0.16666666666666666*crRightHandSideBoundedVector52 + 0.16666666666666666*crRightHandSideBoundedVector53 + 0.16666666666666666*crRightHandSideBoundedVector54;
const double crRightHandSideBoundedVector371 =             -crRightHandSideBoundedVector133*crRightHandSideBoundedVector274 + crRightHandSideBoundedVector370*crRightHandSideBoundedVector96;
const double crRightHandSideBoundedVector372 =             crRightHandSideBoundedVector369 + crRightHandSideBoundedVector371;
const double crRightHandSideBoundedVector373 =             crRightHandSideBoundedVector329*crRightHandSideBoundedVector359;
const double crRightHandSideBoundedVector374 =             0.16666666666666666*crRightHandSideBoundedVector300;
const double crRightHandSideBoundedVector375 =             crRightHandSideBoundedVector199*crRightHandSideBoundedVector374 - crRightHandSideBoundedVector373*crRightHandSideBoundedVector96;
const double crRightHandSideBoundedVector376 =             1.0*crRightHandSideBoundedVector137;
const double crRightHandSideBoundedVector377 =             DN_DX_0_0*crRightHandSideBoundedVector160;
const double crRightHandSideBoundedVector378 =             DN_DX_0_1*crRightHandSideBoundedVector163;
const double crRightHandSideBoundedVector379 =             crRightHandSideBoundedVector139*crRightHandSideBoundedVector370 - crRightHandSideBoundedVector165*crRightHandSideBoundedVector283;
const double crRightHandSideBoundedVector380 =             crRightHandSideBoundedVector378 + crRightHandSideBoundedVector379;
const double crRightHandSideBoundedVector381 =             -crRightHandSideBoundedVector139*crRightHandSideBoundedVector373 + crRightHandSideBoundedVector211*crRightHandSideBoundedVector374;
const double crRightHandSideBoundedVector382 =             1.0*crRightHandSideBoundedVector168;
const double crRightHandSideBoundedVector383 =             nu_sc*(1 - crRightHandSideBoundedVector224);
const double crRightHandSideBoundedVector384 =             crRightHandSideBoundedVector2*crRightHandSideBoundedVector82;
const double crRightHandSideBoundedVector385 =             crRightHandSideBoundedVector2*crRightHandSideBoundedVector65;
const double crRightHandSideBoundedVector386 =             (crRightHandSideBoundedVector2*crRightHandSideBoundedVector225 + crRightHandSideBoundedVector2*crRightHandSideBoundedVector383 + mu)*(-2.2500000000000004*crRightHandSideBoundedVector263 - 2.2500000000000004*crRightHandSideBoundedVector265 + 2.2500000000000004*crRightHandSideBoundedVector384 + 2.2500000000000004*crRightHandSideBoundedVector385);
const double crRightHandSideBoundedVector387 =             crRightHandSideBoundedVector82*crRightHandSideBoundedVector95;
const double crRightHandSideBoundedVector388 =             crRightHandSideBoundedVector65*crRightHandSideBoundedVector95;
const double crRightHandSideBoundedVector389 =             (crRightHandSideBoundedVector225*crRightHandSideBoundedVector95 + crRightHandSideBoundedVector383*crRightHandSideBoundedVector95 + mu)*(-2.2500000000000004*crRightHandSideBoundedVector273 - 2.2500000000000004*crRightHandSideBoundedVector276 + 2.2500000000000004*crRightHandSideBoundedVector387 + 2.2500000000000004*crRightHandSideBoundedVector388);
const double crRightHandSideBoundedVector390 =             crRightHandSideBoundedVector138*crRightHandSideBoundedVector82;
const double crRightHandSideBoundedVector391 =             crRightHandSideBoundedVector138*crRightHandSideBoundedVector65;
const double crRightHandSideBoundedVector392 =             (crRightHandSideBoundedVector138*crRightHandSideBoundedVector225 + crRightHandSideBoundedVector138*crRightHandSideBoundedVector383 + mu)*(-2.2500000000000004*crRightHandSideBoundedVector282 - 2.2500000000000004*crRightHandSideBoundedVector285 + 2.2500000000000004*crRightHandSideBoundedVector390 + 2.2500000000000004*crRightHandSideBoundedVector391);
const double crRightHandSideBoundedVector393 =             0.66666666666666663*crRightHandSideBoundedVector83;
const double crRightHandSideBoundedVector394 =             crRightHandSideBoundedVector299*crRightHandSideBoundedVector359;
const double crRightHandSideBoundedVector395 =             DN_DX_0_0*crRightHandSideBoundedVector262;
const double crRightHandSideBoundedVector396 =             lambda/c_v;
const double crRightHandSideBoundedVector397 =             crRightHandSideBoundedVector38*crRightHandSideBoundedVector396;
const double crRightHandSideBoundedVector398 =             crRightHandSideBoundedVector10*crRightHandSideBoundedVector19;
const double crRightHandSideBoundedVector399 =             crRightHandSideBoundedVector398*crRightHandSideBoundedVector47;
const double crRightHandSideBoundedVector400 =             crRightHandSideBoundedVector19*(crRightHandSideBoundedVector23 - crRightHandSideBoundedVector3*(0.22222222222222221*crRightHandSideBoundedVector13 + 0.22222222222222221*crRightHandSideBoundedVector16));
const double crRightHandSideBoundedVector401 =             crRightHandSideBoundedVector3*(crRightHandSideBoundedVector23 + crRightHandSideBoundedVector400);
const double crRightHandSideBoundedVector402 =             crRightHandSideBoundedVector173*crRightHandSideBoundedVector398;
const double crRightHandSideBoundedVector403 =             0.16666666666666666*r_ext[1];
const double crRightHandSideBoundedVector404 =             0.16666666666666666*r_ext[2];
const double crRightHandSideBoundedVector405 =             crRightHandSideBoundedVector2*(crRightHandSideBoundedVector403 + crRightHandSideBoundedVector404 + 0.66666666666666663*r_ext[0]);
const double crRightHandSideBoundedVector406 =             crRightHandSideBoundedVector44*crRightHandSideBoundedVector59;
const double crRightHandSideBoundedVector407 =             crRightHandSideBoundedVector171*crRightHandSideBoundedVector69;
const double crRightHandSideBoundedVector408 =             crRightHandSideBoundedVector89*gamma;
const double crRightHandSideBoundedVector409 =             crRightHandSideBoundedVector408*crRightHandSideBoundedVector60;
const double crRightHandSideBoundedVector410 =             crRightHandSideBoundedVector184*gamma;
const double crRightHandSideBoundedVector411 =             crRightHandSideBoundedVector410*crRightHandSideBoundedVector70;
const double crRightHandSideBoundedVector412 =             -0.37500000000000006*U_1_3;
const double crRightHandSideBoundedVector413 =             -0.37500000000000006*U_2_3;
const double crRightHandSideBoundedVector414 =             0.75000000000000011*crRightHandSideBoundedVector19;
const double crRightHandSideBoundedVector415 =             1.0/crRightHandSideBoundedVector9;
const double crRightHandSideBoundedVector416 =             crRightHandSideBoundedVector17*crRightHandSideBoundedVector415;
const double crRightHandSideBoundedVector417 =             -1.5000000000000002*U_0_3 - 2.2500000000000004*crRightHandSideBoundedVector400 + crRightHandSideBoundedVector412 + crRightHandSideBoundedVector413 + crRightHandSideBoundedVector414*crRightHandSideBoundedVector416;
const double crRightHandSideBoundedVector418 =             crRightHandSideBoundedVector10*crRightHandSideBoundedVector417;
const double crRightHandSideBoundedVector419 =             crRightHandSideBoundedVector182*crRightHandSideBoundedVector418;
const double crRightHandSideBoundedVector420 =             crRightHandSideBoundedVector418*crRightHandSideBoundedVector86;
const double crRightHandSideBoundedVector421 =             crRightHandSideBoundedVector69*crRightHandSideBoundedVector88;
const double crRightHandSideBoundedVector422 =             crRightHandSideBoundedVector65*crRightHandSideBoundedVector87;
const double crRightHandSideBoundedVector423 =             crRightHandSideBoundedVector398*crRightHandSideBoundedVector69;
const double crRightHandSideBoundedVector424 =             0.16666666666666666*ResProj_2_3 + 0.16666666666666666*dUdt_2_3;
const double crRightHandSideBoundedVector425 =             0.16666666666666666*ResProj_1_3 + 0.16666666666666666*dUdt_1_3;
const double crRightHandSideBoundedVector426 =             (0.66666666666666663*ResProj_0_3 - crRightHandSideBoundedVector405 - crRightHandSideBoundedVector406 - crRightHandSideBoundedVector407 + crRightHandSideBoundedVector409 + crRightHandSideBoundedVector411 + crRightHandSideBoundedVector419 + crRightHandSideBoundedVector420 - crRightHandSideBoundedVector421*crRightHandSideBoundedVector83 - crRightHandSideBoundedVector422*crRightHandSideBoundedVector423 + crRightHandSideBoundedVector424 + crRightHandSideBoundedVector425 + crRightHandSideBoundedVector55*(crRightHandSideBoundedVector401 - crRightHandSideBoundedVector402) + crRightHandSideBoundedVector75*(-crRightHandSideBoundedVector399 + crRightHandSideBoundedVector401) + 0.66666666666666663*dUdt_0_3)/(crRightHandSideBoundedVector397*crRightHandSideBoundedVector5 + crRightHandSideBoundedVector40);
const double crRightHandSideBoundedVector427 =             crRightHandSideBoundedVector397*crRightHandSideBoundedVector4;
const double crRightHandSideBoundedVector428 =             crRightHandSideBoundedVector100*crRightHandSideBoundedVector19;
const double crRightHandSideBoundedVector429 =             crRightHandSideBoundedVector121*crRightHandSideBoundedVector428;
const double crRightHandSideBoundedVector430 =             crRightHandSideBoundedVector19*(crRightHandSideBoundedVector108 - crRightHandSideBoundedVector96*(0.22222222222222221*crRightHandSideBoundedVector102 + 0.22222222222222221*crRightHandSideBoundedVector104));
const double crRightHandSideBoundedVector431 =             crRightHandSideBoundedVector96*(crRightHandSideBoundedVector108 + crRightHandSideBoundedVector430);
const double crRightHandSideBoundedVector432 =             crRightHandSideBoundedVector192*crRightHandSideBoundedVector428;
const double crRightHandSideBoundedVector433 =             0.16666666666666666*r_ext[0];
const double crRightHandSideBoundedVector434 =             crRightHandSideBoundedVector95*(crRightHandSideBoundedVector404 + crRightHandSideBoundedVector433 + 0.66666666666666663*r_ext[1]);
const double crRightHandSideBoundedVector435 =             crRightHandSideBoundedVector119*crRightHandSideBoundedVector126;
const double crRightHandSideBoundedVector436 =             crRightHandSideBoundedVector130*crRightHandSideBoundedVector190;
const double crRightHandSideBoundedVector437 =             crRightHandSideBoundedVector127*crRightHandSideBoundedVector408;
const double crRightHandSideBoundedVector438 =             crRightHandSideBoundedVector131*crRightHandSideBoundedVector410;
const double crRightHandSideBoundedVector439 =             -0.37500000000000006*U_0_3;
const double crRightHandSideBoundedVector440 =             1.0/crRightHandSideBoundedVector99;
const double crRightHandSideBoundedVector441 =             crRightHandSideBoundedVector105*crRightHandSideBoundedVector440;
const double crRightHandSideBoundedVector442 =             -1.5000000000000002*U_1_3 + crRightHandSideBoundedVector413 + crRightHandSideBoundedVector414*crRightHandSideBoundedVector441 - 2.2500000000000004*crRightHandSideBoundedVector430 + crRightHandSideBoundedVector439;
const double crRightHandSideBoundedVector443 =             crRightHandSideBoundedVector100*crRightHandSideBoundedVector442;
const double crRightHandSideBoundedVector444 =             crRightHandSideBoundedVector198*crRightHandSideBoundedVector443;
const double crRightHandSideBoundedVector445 =             crRightHandSideBoundedVector133*crRightHandSideBoundedVector443;
const double crRightHandSideBoundedVector446 =             crRightHandSideBoundedVector130*crRightHandSideBoundedVector135;
const double crRightHandSideBoundedVector447 =             crRightHandSideBoundedVector134*crRightHandSideBoundedVector65;
const double crRightHandSideBoundedVector448 =             crRightHandSideBoundedVector130*crRightHandSideBoundedVector428;
const double crRightHandSideBoundedVector449 =             0.16666666666666666*ResProj_0_3 + 0.16666666666666666*dUdt_0_3;
const double crRightHandSideBoundedVector450 =             (0.66666666666666663*ResProj_1_3 + crRightHandSideBoundedVector424 - crRightHandSideBoundedVector434 - crRightHandSideBoundedVector435 - crRightHandSideBoundedVector436 + crRightHandSideBoundedVector437 + crRightHandSideBoundedVector438 + crRightHandSideBoundedVector444 + crRightHandSideBoundedVector445 - crRightHandSideBoundedVector446*crRightHandSideBoundedVector83 - crRightHandSideBoundedVector447*crRightHandSideBoundedVector448 + crRightHandSideBoundedVector449 + crRightHandSideBoundedVector55*(crRightHandSideBoundedVector431 - crRightHandSideBoundedVector432) + crRightHandSideBoundedVector75*(-crRightHandSideBoundedVector429 + crRightHandSideBoundedVector431) + 0.66666666666666663*dUdt_1_3)/(crRightHandSideBoundedVector116 + crRightHandSideBoundedVector427*crRightHandSideBoundedVector96);
const double crRightHandSideBoundedVector451 =             crRightHandSideBoundedVector141*crRightHandSideBoundedVector19;
const double crRightHandSideBoundedVector452 =             crRightHandSideBoundedVector155*crRightHandSideBoundedVector451;
const double crRightHandSideBoundedVector453 =             crRightHandSideBoundedVector19*(-crRightHandSideBoundedVector139*(0.22222222222222221*crRightHandSideBoundedVector142 + 0.22222222222222221*crRightHandSideBoundedVector143) + crRightHandSideBoundedVector146);
const double crRightHandSideBoundedVector454 =             crRightHandSideBoundedVector139*(crRightHandSideBoundedVector146 + crRightHandSideBoundedVector453);
const double crRightHandSideBoundedVector455 =             crRightHandSideBoundedVector204*crRightHandSideBoundedVector451;
const double crRightHandSideBoundedVector456 =             crRightHandSideBoundedVector138*(crRightHandSideBoundedVector403 + crRightHandSideBoundedVector433 + 0.66666666666666663*r_ext[2]);
const double crRightHandSideBoundedVector457 =             crRightHandSideBoundedVector153*crRightHandSideBoundedVector159;
const double crRightHandSideBoundedVector458 =             crRightHandSideBoundedVector162*crRightHandSideBoundedVector202;
const double crRightHandSideBoundedVector459 =             crRightHandSideBoundedVector160*crRightHandSideBoundedVector408;
const double crRightHandSideBoundedVector460 =             crRightHandSideBoundedVector163*crRightHandSideBoundedVector410;
const double crRightHandSideBoundedVector461 =             1.0/crRightHandSideBoundedVector140;
const double crRightHandSideBoundedVector462 =             crRightHandSideBoundedVector144*crRightHandSideBoundedVector461;
const double crRightHandSideBoundedVector463 =             -1.5000000000000002*U_2_3 + crRightHandSideBoundedVector412 + crRightHandSideBoundedVector414*crRightHandSideBoundedVector462 + crRightHandSideBoundedVector439 - 2.2500000000000004*crRightHandSideBoundedVector453;
const double crRightHandSideBoundedVector464 =             crRightHandSideBoundedVector141*crRightHandSideBoundedVector463;
const double crRightHandSideBoundedVector465 =             crRightHandSideBoundedVector210*crRightHandSideBoundedVector464;
const double crRightHandSideBoundedVector466 =             crRightHandSideBoundedVector165*crRightHandSideBoundedVector464;
const double crRightHandSideBoundedVector467 =             crRightHandSideBoundedVector162*crRightHandSideBoundedVector167;
const double crRightHandSideBoundedVector468 =             crRightHandSideBoundedVector166*crRightHandSideBoundedVector65;
const double crRightHandSideBoundedVector469 =             crRightHandSideBoundedVector162*crRightHandSideBoundedVector451;
const double crRightHandSideBoundedVector470 =             (0.66666666666666663*ResProj_2_3 + crRightHandSideBoundedVector425 + crRightHandSideBoundedVector449 - crRightHandSideBoundedVector456 - crRightHandSideBoundedVector457 - crRightHandSideBoundedVector458 + crRightHandSideBoundedVector459 + crRightHandSideBoundedVector460 + crRightHandSideBoundedVector465 + crRightHandSideBoundedVector466 - crRightHandSideBoundedVector467*crRightHandSideBoundedVector83 - crRightHandSideBoundedVector468*crRightHandSideBoundedVector469 + crRightHandSideBoundedVector55*(crRightHandSideBoundedVector454 - crRightHandSideBoundedVector455) + crRightHandSideBoundedVector75*(-crRightHandSideBoundedVector452 + crRightHandSideBoundedVector454) + 0.66666666666666663*dUdt_2_3)/(crRightHandSideBoundedVector139*crRightHandSideBoundedVector427 + crRightHandSideBoundedVector151);
const double crRightHandSideBoundedVector471 =             0.16666666666666666*crRightHandSideBoundedVector154 - 0.16666666666666666*crRightHandSideBoundedVector158 + crRightHandSideBoundedVector160*crRightHandSideBoundedVector373 - 0.16666666666666666*crRightHandSideBoundedVector161 + crRightHandSideBoundedVector163*crRightHandSideBoundedVector272 - 0.16666666666666666*crRightHandSideBoundedVector164 + crRightHandSideBoundedVector165*crRightHandSideBoundedVector346;
const double crRightHandSideBoundedVector472 =             0.16666666666666666*crRightHandSideBoundedVector120 - 0.16666666666666666*crRightHandSideBoundedVector124 + crRightHandSideBoundedVector127*crRightHandSideBoundedVector373 - 0.16666666666666666*crRightHandSideBoundedVector128 + crRightHandSideBoundedVector131*crRightHandSideBoundedVector272 - 0.16666666666666666*crRightHandSideBoundedVector132 + crRightHandSideBoundedVector133*crRightHandSideBoundedVector322;
const double crRightHandSideBoundedVector473 =             crRightHandSideBoundedVector223*pow(lin_m[1], 2);
const double crRightHandSideBoundedVector474 =             crRightHandSideBoundedVector473*nu_st;
const double crRightHandSideBoundedVector475 =             1 - crRightHandSideBoundedVector473;
const double crRightHandSideBoundedVector476 =             crRightHandSideBoundedVector475*nu_sc;
const double crRightHandSideBoundedVector477 =             crRightHandSideBoundedVector217*crRightHandSideBoundedVector219 + crRightHandSideBoundedVector220 + crRightHandSideBoundedVector222*crRightHandSideBoundedVector474 + crRightHandSideBoundedVector222*crRightHandSideBoundedVector476 + crRightHandSideBoundedVector225*crRightHandSideBoundedVector227 - crRightHandSideBoundedVector227*crRightHandSideBoundedVector230;
const double crRightHandSideBoundedVector478 =             crRightHandSideBoundedVector225*crRightHandSideBoundedVector244 - crRightHandSideBoundedVector230*crRightHandSideBoundedVector244 + crRightHandSideBoundedVector237*crRightHandSideBoundedVector239 + crRightHandSideBoundedVector240 + crRightHandSideBoundedVector242*crRightHandSideBoundedVector474 + crRightHandSideBoundedVector242*crRightHandSideBoundedVector476;
const double crRightHandSideBoundedVector479 =             crRightHandSideBoundedVector225*crRightHandSideBoundedVector256 - crRightHandSideBoundedVector230*crRightHandSideBoundedVector256 + crRightHandSideBoundedVector249*crRightHandSideBoundedVector251 + crRightHandSideBoundedVector252 + crRightHandSideBoundedVector254*crRightHandSideBoundedVector474 + crRightHandSideBoundedVector254*crRightHandSideBoundedVector476;
const double crRightHandSideBoundedVector480 =             3.0000000000000004*crRightHandSideBoundedVector16;
const double crRightHandSideBoundedVector481 =             crRightHandSideBoundedVector85*(crRightHandSideBoundedVector292 - crRightHandSideBoundedVector480);
const double crRightHandSideBoundedVector482 =             crRightHandSideBoundedVector182*crRightHandSideBoundedVector69;
const double crRightHandSideBoundedVector483 =             crRightHandSideBoundedVector296*crRightHandSideBoundedVector482;
const double crRightHandSideBoundedVector484 =             0.66666666666666663*crRightHandSideBoundedVector55;
const double crRightHandSideBoundedVector485 =             crRightHandSideBoundedVector301*crRightHandSideBoundedVector69;
const double crRightHandSideBoundedVector486 =             crRightHandSideBoundedVector267*crRightHandSideBoundedVector59;
const double crRightHandSideBoundedVector487 =             crRightHandSideBoundedVector486*crRightHandSideBoundedVector65;
const double crRightHandSideBoundedVector488 =             0.1111111111111111*f_ext(1,1);
const double crRightHandSideBoundedVector489 =             0.1111111111111111*f_ext(2,1);
const double crRightHandSideBoundedVector490 =             crRightHandSideBoundedVector488 + crRightHandSideBoundedVector489 + 0.44444444444444442*f_ext(0,1);
const double crRightHandSideBoundedVector491 =             0.16666666666666666*crRightHandSideBoundedVector85;
const double crRightHandSideBoundedVector492 =             3.0000000000000004*crRightHandSideBoundedVector104;
const double crRightHandSideBoundedVector493 =             crRightHandSideBoundedVector318*(crRightHandSideBoundedVector320 - crRightHandSideBoundedVector492);
const double crRightHandSideBoundedVector494 =             crRightHandSideBoundedVector324*crRightHandSideBoundedVector75;
const double crRightHandSideBoundedVector495 =             crRightHandSideBoundedVector130*crRightHandSideBoundedVector198;
const double crRightHandSideBoundedVector496 =             crRightHandSideBoundedVector326*crRightHandSideBoundedVector495;
const double crRightHandSideBoundedVector497 =             0.16666666666666666*crRightHandSideBoundedVector55;
const double crRightHandSideBoundedVector498 =             crRightHandSideBoundedVector300*crRightHandSideBoundedVector497;
const double crRightHandSideBoundedVector499 =             crRightHandSideBoundedVector100*crRightHandSideBoundedVector130;
const double crRightHandSideBoundedVector500 =             crRightHandSideBoundedVector179*crRightHandSideBoundedVector322;
const double crRightHandSideBoundedVector501 =             0.027777777777777776*f_ext(0,1);
const double crRightHandSideBoundedVector502 =             0.027777777777777776*f_ext(2,1);
const double crRightHandSideBoundedVector503 =             crRightHandSideBoundedVector488 + crRightHandSideBoundedVector501 + crRightHandSideBoundedVector502;
const double crRightHandSideBoundedVector504 =             -crRightHandSideBoundedVector322*crRightHandSideBoundedVector82 - crRightHandSideBoundedVector491*crRightHandSideBoundedVector493 - crRightHandSideBoundedVector494 + crRightHandSideBoundedVector496 + crRightHandSideBoundedVector498*crRightHandSideBoundedVector499 + crRightHandSideBoundedVector500 + crRightHandSideBoundedVector503;
const double crRightHandSideBoundedVector505 =             3.0000000000000004*crRightHandSideBoundedVector143;
const double crRightHandSideBoundedVector506 =             crRightHandSideBoundedVector342*(crRightHandSideBoundedVector344 - crRightHandSideBoundedVector505);
const double crRightHandSideBoundedVector507 =             crRightHandSideBoundedVector348*crRightHandSideBoundedVector75;
const double crRightHandSideBoundedVector508 =             crRightHandSideBoundedVector162*crRightHandSideBoundedVector210;
const double crRightHandSideBoundedVector509 =             crRightHandSideBoundedVector350*crRightHandSideBoundedVector508;
const double crRightHandSideBoundedVector510 =             crRightHandSideBoundedVector141*crRightHandSideBoundedVector162;
const double crRightHandSideBoundedVector511 =             crRightHandSideBoundedVector179*crRightHandSideBoundedVector346;
const double crRightHandSideBoundedVector512 =             0.027777777777777776*f_ext(1,1);
const double crRightHandSideBoundedVector513 =             crRightHandSideBoundedVector489 + crRightHandSideBoundedVector501 + crRightHandSideBoundedVector512;
const double crRightHandSideBoundedVector514 =             -crRightHandSideBoundedVector346*crRightHandSideBoundedVector82 - crRightHandSideBoundedVector491*crRightHandSideBoundedVector506 + crRightHandSideBoundedVector498*crRightHandSideBoundedVector510 - crRightHandSideBoundedVector507 + crRightHandSideBoundedVector509 + crRightHandSideBoundedVector511 + crRightHandSideBoundedVector513;
const double crRightHandSideBoundedVector515 =             0.66666666666666663*crRightHandSideBoundedVector72 + 0.66666666666666663*crRightHandSideBoundedVector73 + 0.66666666666666663*crRightHandSideBoundedVector74;
const double crRightHandSideBoundedVector516 =             1.5000000000000002*crRightHandSideBoundedVector182;
const double crRightHandSideBoundedVector517 =             crRightHandSideBoundedVector10*crRightHandSideBoundedVector516;
const double crRightHandSideBoundedVector518 =             crRightHandSideBoundedVector3*crRightHandSideBoundedVector515 + crRightHandSideBoundedVector361 - crRightHandSideBoundedVector517;
const double crRightHandSideBoundedVector519 =             DN_DX_0_0*crRightHandSideBoundedVector66 + 0.66666666666666663*crRightHandSideBoundedVector80 + 0.66666666666666663*crRightHandSideBoundedVector81;
const double crRightHandSideBoundedVector520 =             0.16666666666666666*crRightHandSideBoundedVector72 + 0.16666666666666666*crRightHandSideBoundedVector73 + 0.16666666666666666*crRightHandSideBoundedVector74;
const double crRightHandSideBoundedVector521 =             -crRightHandSideBoundedVector198*crRightHandSideBoundedVector274 + crRightHandSideBoundedVector520*crRightHandSideBoundedVector96;
const double crRightHandSideBoundedVector522 =             crRightHandSideBoundedVector368 + crRightHandSideBoundedVector521;
const double crRightHandSideBoundedVector523 =             crRightHandSideBoundedVector359*crRightHandSideBoundedVector497;
const double crRightHandSideBoundedVector524 =             crRightHandSideBoundedVector100*crRightHandSideBoundedVector133*crRightHandSideBoundedVector374 - crRightHandSideBoundedVector523*crRightHandSideBoundedVector96;
const double crRightHandSideBoundedVector525 =             DN_DX_1_0*crRightHandSideBoundedVector67 + DN_DX_2_0*crRightHandSideBoundedVector68 + 0.16666666666666666*crRightHandSideBoundedVector79;
const double crRightHandSideBoundedVector526 =             0.16666666666666666*crRightHandSideBoundedVector179;
const double crRightHandSideBoundedVector527 =             crRightHandSideBoundedVector19*crRightHandSideBoundedVector275 - crRightHandSideBoundedVector277 + crRightHandSideBoundedVector525*crRightHandSideBoundedVector96 - crRightHandSideBoundedVector526*crRightHandSideBoundedVector96;
const double crRightHandSideBoundedVector528 =             crRightHandSideBoundedVector139*crRightHandSideBoundedVector520 - crRightHandSideBoundedVector210*crRightHandSideBoundedVector283;
const double crRightHandSideBoundedVector529 =             crRightHandSideBoundedVector377 + crRightHandSideBoundedVector528;
const double crRightHandSideBoundedVector530 =             -crRightHandSideBoundedVector139*crRightHandSideBoundedVector523 + crRightHandSideBoundedVector141*crRightHandSideBoundedVector165*crRightHandSideBoundedVector374;
const double crRightHandSideBoundedVector531 =             crRightHandSideBoundedVector139*crRightHandSideBoundedVector525 - crRightHandSideBoundedVector139*crRightHandSideBoundedVector526 + crRightHandSideBoundedVector19*crRightHandSideBoundedVector284 - crRightHandSideBoundedVector286;
const double crRightHandSideBoundedVector532 =             crRightHandSideBoundedVector10*crRightHandSideBoundedVector386;
const double crRightHandSideBoundedVector533 =             crRightHandSideBoundedVector100*crRightHandSideBoundedVector389;
const double crRightHandSideBoundedVector534 =             crRightHandSideBoundedVector141*crRightHandSideBoundedVector392;
const double crRightHandSideBoundedVector535 =             0.66666666666666663*crRightHandSideBoundedVector179;
const double crRightHandSideBoundedVector536 =             crRightHandSideBoundedVector359*crRightHandSideBoundedVector484;
const double crRightHandSideBoundedVector537 =             DN_DX_0_1*crRightHandSideBoundedVector262;
const double crRightHandSideBoundedVector538 =             crRightHandSideBoundedVector160*crRightHandSideBoundedVector526 + crRightHandSideBoundedVector163*crRightHandSideBoundedVector523 + 0.16666666666666666*crRightHandSideBoundedVector203 - 0.16666666666666666*crRightHandSideBoundedVector206 - 0.16666666666666666*crRightHandSideBoundedVector207 - 0.16666666666666666*crRightHandSideBoundedVector208 + crRightHandSideBoundedVector210*crRightHandSideBoundedVector348;
const double crRightHandSideBoundedVector539 =             crRightHandSideBoundedVector127*crRightHandSideBoundedVector526 + crRightHandSideBoundedVector131*crRightHandSideBoundedVector523 + 0.16666666666666666*crRightHandSideBoundedVector191 - 0.16666666666666666*crRightHandSideBoundedVector194 - 0.16666666666666666*crRightHandSideBoundedVector195 - 0.16666666666666666*crRightHandSideBoundedVector196 + crRightHandSideBoundedVector198*crRightHandSideBoundedVector324;
const double crRightHandSideBoundedVector540 =             -crRightHandSideBoundedVector401;
const double crRightHandSideBoundedVector541 =             crRightHandSideBoundedVector399 + crRightHandSideBoundedVector540;
const double crRightHandSideBoundedVector542 =             crRightHandSideBoundedVector402 + crRightHandSideBoundedVector540;
const double crRightHandSideBoundedVector543 =             crRightHandSideBoundedVector70*mu;
const double crRightHandSideBoundedVector544 =             -2.2500000000000004*crRightHandSideBoundedVector263 - 2.2500000000000004*crRightHandSideBoundedVector265 + 2.2500000000000004*crRightHandSideBoundedVector384 + 2.2500000000000004*crRightHandSideBoundedVector385;
const double crRightHandSideBoundedVector545 =             crRightHandSideBoundedVector60*mu;
const double crRightHandSideBoundedVector546 =             2.2500000000000004*crRightHandSideBoundedVector2;
const double crRightHandSideBoundedVector547 =             2.2500000000000004*crRightHandSideBoundedVector23;
const double crRightHandSideBoundedVector548 =             1.5000000000000002*crRightHandSideBoundedVector415;
const double crRightHandSideBoundedVector549 =             crRightHandSideBoundedVector46*crRightHandSideBoundedVector548;
const double crRightHandSideBoundedVector550 =             crRightHandSideBoundedVector13*crRightHandSideBoundedVector549 + crRightHandSideBoundedVector16*crRightHandSideBoundedVector549 - crRightHandSideBoundedVector181*crRightHandSideBoundedVector82 - crRightHandSideBoundedVector46*crRightHandSideBoundedVector547 + crRightHandSideBoundedVector546*crRightHandSideBoundedVector89 - crRightHandSideBoundedVector75*crRightHandSideBoundedVector87;
const double crRightHandSideBoundedVector551 =             gamma*k_st;
const double crRightHandSideBoundedVector552 =             crRightHandSideBoundedVector548*crRightHandSideBoundedVector85;
const double crRightHandSideBoundedVector553 =             crRightHandSideBoundedVector13*crRightHandSideBoundedVector552 + crRightHandSideBoundedVector16*crRightHandSideBoundedVector552 - crRightHandSideBoundedVector181*crRightHandSideBoundedVector55 + crRightHandSideBoundedVector184*crRightHandSideBoundedVector546 - crRightHandSideBoundedVector422 - crRightHandSideBoundedVector547*crRightHandSideBoundedVector85;
const double crRightHandSideBoundedVector554 =             crRightHandSideBoundedVector2*crRightHandSideBoundedVector553;
const double crRightHandSideBoundedVector555 =             crRightHandSideBoundedVector224*crRightHandSideBoundedVector554;
const double crRightHandSideBoundedVector556 =             crRightHandSideBoundedVector2*crRightHandSideBoundedVector550;
const double crRightHandSideBoundedVector557 =             crRightHandSideBoundedVector228*crRightHandSideBoundedVector551;
const double crRightHandSideBoundedVector558 =             gamma*k_sc;
const double crRightHandSideBoundedVector559 =             crRightHandSideBoundedVector231*crRightHandSideBoundedVector558;
const double crRightHandSideBoundedVector560 =             crRightHandSideBoundedVector10*(crRightHandSideBoundedVector396*crRightHandSideBoundedVector550 + crRightHandSideBoundedVector543*crRightHandSideBoundedVector544 + crRightHandSideBoundedVector545*(-3.0000000000000009*crRightHandSideBoundedVector182 + 3.0000000000000009*crRightHandSideBoundedVector214 - 1.5000000000000002*crRightHandSideBoundedVector218 + crRightHandSideBoundedVector364) + crRightHandSideBoundedVector551*crRightHandSideBoundedVector555 - crRightHandSideBoundedVector555*crRightHandSideBoundedVector558 + crRightHandSideBoundedVector556*crRightHandSideBoundedVector557 + crRightHandSideBoundedVector556*crRightHandSideBoundedVector559);
const double crRightHandSideBoundedVector561 =             crRightHandSideBoundedVector131*mu;
const double crRightHandSideBoundedVector562 =             -2.2500000000000004*crRightHandSideBoundedVector273 - 2.2500000000000004*crRightHandSideBoundedVector276 + 2.2500000000000004*crRightHandSideBoundedVector387 + 2.2500000000000004*crRightHandSideBoundedVector388;
const double crRightHandSideBoundedVector563 =             1.5000000000000002*crRightHandSideBoundedVector133;
const double crRightHandSideBoundedVector564 =             crRightHandSideBoundedVector127*mu;
const double crRightHandSideBoundedVector565 =             2.2500000000000004*crRightHandSideBoundedVector95;
const double crRightHandSideBoundedVector566 =             2.2500000000000004*crRightHandSideBoundedVector108;
const double crRightHandSideBoundedVector567 =             1.5000000000000002*crRightHandSideBoundedVector440;
const double crRightHandSideBoundedVector568 =             crRightHandSideBoundedVector46*crRightHandSideBoundedVector567;
const double crRightHandSideBoundedVector569 =             crRightHandSideBoundedVector102*crRightHandSideBoundedVector568 + crRightHandSideBoundedVector104*crRightHandSideBoundedVector568 - crRightHandSideBoundedVector134*crRightHandSideBoundedVector75 - crRightHandSideBoundedVector197*crRightHandSideBoundedVector82 - crRightHandSideBoundedVector46*crRightHandSideBoundedVector566 + crRightHandSideBoundedVector565*crRightHandSideBoundedVector89;
const double crRightHandSideBoundedVector570 =             crRightHandSideBoundedVector567*crRightHandSideBoundedVector85;
const double crRightHandSideBoundedVector571 =             crRightHandSideBoundedVector102*crRightHandSideBoundedVector570 + crRightHandSideBoundedVector104*crRightHandSideBoundedVector570 + crRightHandSideBoundedVector184*crRightHandSideBoundedVector565 - crRightHandSideBoundedVector197*crRightHandSideBoundedVector55 - crRightHandSideBoundedVector447 - crRightHandSideBoundedVector566*crRightHandSideBoundedVector85;
const double crRightHandSideBoundedVector572 =             crRightHandSideBoundedVector571*crRightHandSideBoundedVector95;
const double crRightHandSideBoundedVector573 =             crRightHandSideBoundedVector224*crRightHandSideBoundedVector572;
const double crRightHandSideBoundedVector574 =             crRightHandSideBoundedVector569*crRightHandSideBoundedVector95;
const double crRightHandSideBoundedVector575 =             crRightHandSideBoundedVector100*(crRightHandSideBoundedVector396*crRightHandSideBoundedVector569 + crRightHandSideBoundedVector551*crRightHandSideBoundedVector573 + crRightHandSideBoundedVector557*crRightHandSideBoundedVector574 - crRightHandSideBoundedVector558*crRightHandSideBoundedVector573 + crRightHandSideBoundedVector559*crRightHandSideBoundedVector574 + crRightHandSideBoundedVector561*crRightHandSideBoundedVector562 + crRightHandSideBoundedVector564*(-3.0000000000000009*crRightHandSideBoundedVector198 + 3.0000000000000009*crRightHandSideBoundedVector234 - 1.5000000000000002*crRightHandSideBoundedVector238 + crRightHandSideBoundedVector563));
const double crRightHandSideBoundedVector576 =             crRightHandSideBoundedVector163*mu;
const double crRightHandSideBoundedVector577 =             -2.2500000000000004*crRightHandSideBoundedVector282 - 2.2500000000000004*crRightHandSideBoundedVector285 + 2.2500000000000004*crRightHandSideBoundedVector390 + 2.2500000000000004*crRightHandSideBoundedVector391;
const double crRightHandSideBoundedVector578 =             1.5000000000000002*crRightHandSideBoundedVector165;
const double crRightHandSideBoundedVector579 =             crRightHandSideBoundedVector160*mu;
const double crRightHandSideBoundedVector580 =             2.2500000000000004*crRightHandSideBoundedVector138;
const double crRightHandSideBoundedVector581 =             2.2500000000000004*crRightHandSideBoundedVector146;
const double crRightHandSideBoundedVector582 =             1.5000000000000002*crRightHandSideBoundedVector461;
const double crRightHandSideBoundedVector583 =             crRightHandSideBoundedVector46*crRightHandSideBoundedVector582;
const double crRightHandSideBoundedVector584 =             crRightHandSideBoundedVector142*crRightHandSideBoundedVector583 + crRightHandSideBoundedVector143*crRightHandSideBoundedVector583 - crRightHandSideBoundedVector166*crRightHandSideBoundedVector75 - crRightHandSideBoundedVector209*crRightHandSideBoundedVector82 - crRightHandSideBoundedVector46*crRightHandSideBoundedVector581 + crRightHandSideBoundedVector580*crRightHandSideBoundedVector89;
const double crRightHandSideBoundedVector585 =             crRightHandSideBoundedVector582*crRightHandSideBoundedVector85;
const double crRightHandSideBoundedVector586 =             crRightHandSideBoundedVector142*crRightHandSideBoundedVector585 + crRightHandSideBoundedVector143*crRightHandSideBoundedVector585 + crRightHandSideBoundedVector184*crRightHandSideBoundedVector580 - crRightHandSideBoundedVector209*crRightHandSideBoundedVector55 - crRightHandSideBoundedVector468 - crRightHandSideBoundedVector581*crRightHandSideBoundedVector85;
const double crRightHandSideBoundedVector587 =             crRightHandSideBoundedVector138*crRightHandSideBoundedVector586;
const double crRightHandSideBoundedVector588 =             crRightHandSideBoundedVector224*crRightHandSideBoundedVector587;
const double crRightHandSideBoundedVector589 =             crRightHandSideBoundedVector138*crRightHandSideBoundedVector584;
const double crRightHandSideBoundedVector590 =             crRightHandSideBoundedVector141*(crRightHandSideBoundedVector396*crRightHandSideBoundedVector584 + crRightHandSideBoundedVector551*crRightHandSideBoundedVector588 + crRightHandSideBoundedVector557*crRightHandSideBoundedVector589 - crRightHandSideBoundedVector558*crRightHandSideBoundedVector588 + crRightHandSideBoundedVector559*crRightHandSideBoundedVector589 + crRightHandSideBoundedVector576*crRightHandSideBoundedVector577 + crRightHandSideBoundedVector579*(-3.0000000000000009*crRightHandSideBoundedVector210 + 3.0000000000000009*crRightHandSideBoundedVector246 - 1.5000000000000002*crRightHandSideBoundedVector250 + crRightHandSideBoundedVector578));
const double crRightHandSideBoundedVector591 =             crRightHandSideBoundedVector224*crRightHandSideBoundedVector556;
const double crRightHandSideBoundedVector592 =             crRightHandSideBoundedVector473*crRightHandSideBoundedVector551;
const double crRightHandSideBoundedVector593 =             crRightHandSideBoundedVector475*crRightHandSideBoundedVector558;
const double crRightHandSideBoundedVector594 =             crRightHandSideBoundedVector10*(crRightHandSideBoundedVector396*crRightHandSideBoundedVector553 + crRightHandSideBoundedVector543*(-1.5000000000000002*crRightHandSideBoundedVector214 + 3.0000000000000009*crRightHandSideBoundedVector218 + crRightHandSideBoundedVector516 - 3.0000000000000009*crRightHandSideBoundedVector86) + crRightHandSideBoundedVector544*crRightHandSideBoundedVector545 + crRightHandSideBoundedVector551*crRightHandSideBoundedVector591 + crRightHandSideBoundedVector554*crRightHandSideBoundedVector592 + crRightHandSideBoundedVector554*crRightHandSideBoundedVector593 - crRightHandSideBoundedVector558*crRightHandSideBoundedVector591);
const double crRightHandSideBoundedVector595 =             1.5000000000000002*crRightHandSideBoundedVector198;
const double crRightHandSideBoundedVector596 =             crRightHandSideBoundedVector224*crRightHandSideBoundedVector574;
const double crRightHandSideBoundedVector597 =             crRightHandSideBoundedVector100*(crRightHandSideBoundedVector396*crRightHandSideBoundedVector571 + crRightHandSideBoundedVector551*crRightHandSideBoundedVector596 - crRightHandSideBoundedVector558*crRightHandSideBoundedVector596 + crRightHandSideBoundedVector561*(-3.0000000000000009*crRightHandSideBoundedVector133 - 1.5000000000000002*crRightHandSideBoundedVector234 + 3.0000000000000009*crRightHandSideBoundedVector238 + crRightHandSideBoundedVector595) + crRightHandSideBoundedVector562*crRightHandSideBoundedVector564 + crRightHandSideBoundedVector572*crRightHandSideBoundedVector592 + crRightHandSideBoundedVector572*crRightHandSideBoundedVector593);
const double crRightHandSideBoundedVector598 =             1.5000000000000002*crRightHandSideBoundedVector210;
const double crRightHandSideBoundedVector599 =             crRightHandSideBoundedVector224*crRightHandSideBoundedVector589;
const double crRightHandSideBoundedVector600 =             crRightHandSideBoundedVector141*(crRightHandSideBoundedVector396*crRightHandSideBoundedVector586 + crRightHandSideBoundedVector551*crRightHandSideBoundedVector599 - crRightHandSideBoundedVector558*crRightHandSideBoundedVector599 + crRightHandSideBoundedVector576*(-3.0000000000000009*crRightHandSideBoundedVector165 - 1.5000000000000002*crRightHandSideBoundedVector246 + 3.0000000000000009*crRightHandSideBoundedVector250 + crRightHandSideBoundedVector598) + crRightHandSideBoundedVector577*crRightHandSideBoundedVector579 + crRightHandSideBoundedVector587*crRightHandSideBoundedVector592 + crRightHandSideBoundedVector587*crRightHandSideBoundedVector593);
const double crRightHandSideBoundedVector601 =             0.1111111111111111*r_ext[1];
const double crRightHandSideBoundedVector602 =             0.1111111111111111*r_ext[2];
const double crRightHandSideBoundedVector603 =             crRightHandSideBoundedVector59*crRightHandSideBoundedVector69;
const double crRightHandSideBoundedVector604 =             crRightHandSideBoundedVector296*crRightHandSideBoundedVector603;
const double crRightHandSideBoundedVector605 =             crRightHandSideBoundedVector418*crRightHandSideBoundedVector59;
const double crRightHandSideBoundedVector606 =             -1.125*U_1_3;
const double crRightHandSideBoundedVector607 =             -1.125*U_2_3;
const double crRightHandSideBoundedVector608 =             4.5000000000000009*crRightHandSideBoundedVector19;
const double crRightHandSideBoundedVector609 =             -4.5*U_0_3 - 6.7500000000000009*crRightHandSideBoundedVector400 + crRightHandSideBoundedVector416*crRightHandSideBoundedVector608 + crRightHandSideBoundedVector606 + crRightHandSideBoundedVector607;
const double crRightHandSideBoundedVector610 =             crRightHandSideBoundedVector289*crRightHandSideBoundedVector609;
const double crRightHandSideBoundedVector611 =             crRightHandSideBoundedVector19*crRightHandSideBoundedVector415;
const double crRightHandSideBoundedVector612 =             crRightHandSideBoundedVector10*(crRightHandSideBoundedVector291*crRightHandSideBoundedVector611 + crRightHandSideBoundedVector417);
const double crRightHandSideBoundedVector613 =             crRightHandSideBoundedVector10*(crRightHandSideBoundedVector417 + crRightHandSideBoundedVector480*crRightHandSideBoundedVector611);
const double crRightHandSideBoundedVector614 =             crRightHandSideBoundedVector126*crRightHandSideBoundedVector443;
const double crRightHandSideBoundedVector615 =             0.027777777777777776*r_ext[0];
const double crRightHandSideBoundedVector616 =             0.027777777777777776*r_ext[2];
const double crRightHandSideBoundedVector617 =             crRightHandSideBoundedVector19*crRightHandSideBoundedVector440;
const double crRightHandSideBoundedVector618 =             crRightHandSideBoundedVector100*(crRightHandSideBoundedVector319*crRightHandSideBoundedVector617 + crRightHandSideBoundedVector442);
const double crRightHandSideBoundedVector619 =             crRightHandSideBoundedVector100*(crRightHandSideBoundedVector442 + crRightHandSideBoundedVector492*crRightHandSideBoundedVector617);
const double crRightHandSideBoundedVector620 =             -1.125*U_0_3;
const double crRightHandSideBoundedVector621 =             crRightHandSideBoundedVector318*(-4.5*U_1_3 - 6.7500000000000009*crRightHandSideBoundedVector430 + crRightHandSideBoundedVector441*crRightHandSideBoundedVector608 + crRightHandSideBoundedVector607 + crRightHandSideBoundedVector620);
const double crRightHandSideBoundedVector622 =             0.16666666666666666*crRightHandSideBoundedVector621;
const double crRightHandSideBoundedVector623 =             crRightHandSideBoundedVector126*crRightHandSideBoundedVector130;
const double crRightHandSideBoundedVector624 =             crRightHandSideBoundedVector326*crRightHandSideBoundedVector623;
const double crRightHandSideBoundedVector625 =             -crRightHandSideBoundedVector133*crRightHandSideBoundedVector622 + crRightHandSideBoundedVector179*crRightHandSideBoundedVector624 - crRightHandSideBoundedVector198*crRightHandSideBoundedVector622 - crRightHandSideBoundedVector322*crRightHandSideBoundedVector408 - crRightHandSideBoundedVector324*crRightHandSideBoundedVector410 + crRightHandSideBoundedVector329*crRightHandSideBoundedVector618 + crRightHandSideBoundedVector497*crRightHandSideBoundedVector619 + crRightHandSideBoundedVector601 + crRightHandSideBoundedVector615 + crRightHandSideBoundedVector616 + crRightHandSideBoundedVector624*crRightHandSideBoundedVector83;
const double crRightHandSideBoundedVector626 =             crRightHandSideBoundedVector159*crRightHandSideBoundedVector464;
const double crRightHandSideBoundedVector627 =             0.027777777777777776*r_ext[1];
const double crRightHandSideBoundedVector628 =             crRightHandSideBoundedVector19*crRightHandSideBoundedVector461;
const double crRightHandSideBoundedVector629 =             crRightHandSideBoundedVector141*(crRightHandSideBoundedVector343*crRightHandSideBoundedVector628 + crRightHandSideBoundedVector463);
const double crRightHandSideBoundedVector630 =             crRightHandSideBoundedVector141*(crRightHandSideBoundedVector463 + crRightHandSideBoundedVector505*crRightHandSideBoundedVector628);
const double crRightHandSideBoundedVector631 =             crRightHandSideBoundedVector342*(-4.5*U_2_3 - 6.7500000000000009*crRightHandSideBoundedVector453 + crRightHandSideBoundedVector462*crRightHandSideBoundedVector608 + crRightHandSideBoundedVector606 + crRightHandSideBoundedVector620);
const double crRightHandSideBoundedVector632 =             0.16666666666666666*crRightHandSideBoundedVector631;
const double crRightHandSideBoundedVector633 =             crRightHandSideBoundedVector159*crRightHandSideBoundedVector162;
const double crRightHandSideBoundedVector634 =             crRightHandSideBoundedVector350*crRightHandSideBoundedVector633;
const double crRightHandSideBoundedVector635 =             -crRightHandSideBoundedVector165*crRightHandSideBoundedVector632 + crRightHandSideBoundedVector179*crRightHandSideBoundedVector634 - crRightHandSideBoundedVector210*crRightHandSideBoundedVector632 + crRightHandSideBoundedVector329*crRightHandSideBoundedVector629 - crRightHandSideBoundedVector346*crRightHandSideBoundedVector408 - crRightHandSideBoundedVector348*crRightHandSideBoundedVector410 + crRightHandSideBoundedVector497*crRightHandSideBoundedVector630 + crRightHandSideBoundedVector602 + crRightHandSideBoundedVector615 + crRightHandSideBoundedVector627 + crRightHandSideBoundedVector634*crRightHandSideBoundedVector83;
const double crRightHandSideBoundedVector636 =             0.66666666666666663*crRightHandSideBoundedVector46;
const double crRightHandSideBoundedVector637 =             4.5000000000000009*crRightHandSideBoundedVector398;
const double crRightHandSideBoundedVector638 =             crRightHandSideBoundedVector59*crRightHandSideBoundedVector75;
const double crRightHandSideBoundedVector639 =             crRightHandSideBoundedVector303*crRightHandSideBoundedVector65;
const double crRightHandSideBoundedVector640 =             0.66666666666666663*crRightHandSideBoundedVector85;
const double crRightHandSideBoundedVector641 =             crRightHandSideBoundedVector55*crRightHandSideBoundedVector69;
const double crRightHandSideBoundedVector642 =             crRightHandSideBoundedVector423*crRightHandSideBoundedVector87;
const double crRightHandSideBoundedVector643 =             -crRightHandSideBoundedVector431;
const double crRightHandSideBoundedVector644 =             crRightHandSideBoundedVector429 + crRightHandSideBoundedVector643;
const double crRightHandSideBoundedVector645 =             0.16666666666666666*crRightHandSideBoundedVector96;
const double crRightHandSideBoundedVector646 =             1.1250000000000002*crRightHandSideBoundedVector428;
const double crRightHandSideBoundedVector647 =             crRightHandSideBoundedVector126*crRightHandSideBoundedVector75;
const double crRightHandSideBoundedVector648 =             -crRightHandSideBoundedVector19*crRightHandSideBoundedVector323 - crRightHandSideBoundedVector19*crRightHandSideBoundedVector325 + crRightHandSideBoundedVector19*crRightHandSideBoundedVector328 + crRightHandSideBoundedVector317*crRightHandSideBoundedVector618 - crRightHandSideBoundedVector332 + crRightHandSideBoundedVector335 + crRightHandSideBoundedVector408*crRightHandSideBoundedVector645 - crRightHandSideBoundedVector646*crRightHandSideBoundedVector647;
const double crRightHandSideBoundedVector649 =             crRightHandSideBoundedVector432 + crRightHandSideBoundedVector643;
const double crRightHandSideBoundedVector650 =             crRightHandSideBoundedVector134*crRightHandSideBoundedVector448;
const double crRightHandSideBoundedVector651 =             crRightHandSideBoundedVector130*crRightHandSideBoundedVector55;
const double crRightHandSideBoundedVector652 =             crRightHandSideBoundedVector322*crRightHandSideBoundedVector83;
const double crRightHandSideBoundedVector653 =             -crRightHandSideBoundedVector19*crRightHandSideBoundedVector494 + crRightHandSideBoundedVector19*crRightHandSideBoundedVector496 + crRightHandSideBoundedVector410*crRightHandSideBoundedVector645 + crRightHandSideBoundedVector491*crRightHandSideBoundedVector619 - crRightHandSideBoundedVector500 + crRightHandSideBoundedVector503 - crRightHandSideBoundedVector646*crRightHandSideBoundedVector651 - crRightHandSideBoundedVector652;
const double crRightHandSideBoundedVector654 =             -crRightHandSideBoundedVector454;
const double crRightHandSideBoundedVector655 =             crRightHandSideBoundedVector452 + crRightHandSideBoundedVector654;
const double crRightHandSideBoundedVector656 =             0.16666666666666666*crRightHandSideBoundedVector139;
const double crRightHandSideBoundedVector657 =             1.1250000000000002*crRightHandSideBoundedVector451;
const double crRightHandSideBoundedVector658 =             crRightHandSideBoundedVector159*crRightHandSideBoundedVector75;
const double crRightHandSideBoundedVector659 =             -crRightHandSideBoundedVector19*crRightHandSideBoundedVector347 - crRightHandSideBoundedVector19*crRightHandSideBoundedVector349 + crRightHandSideBoundedVector19*crRightHandSideBoundedVector352 + crRightHandSideBoundedVector317*crRightHandSideBoundedVector629 - crRightHandSideBoundedVector354 + crRightHandSideBoundedVector356 + crRightHandSideBoundedVector408*crRightHandSideBoundedVector656 - crRightHandSideBoundedVector657*crRightHandSideBoundedVector658;
const double crRightHandSideBoundedVector660 =             crRightHandSideBoundedVector455 + crRightHandSideBoundedVector654;
const double crRightHandSideBoundedVector661 =             crRightHandSideBoundedVector166*crRightHandSideBoundedVector469;
const double crRightHandSideBoundedVector662 =             crRightHandSideBoundedVector162*crRightHandSideBoundedVector55;
const double crRightHandSideBoundedVector663 =             crRightHandSideBoundedVector346*crRightHandSideBoundedVector83;
const double crRightHandSideBoundedVector664 =             -crRightHandSideBoundedVector19*crRightHandSideBoundedVector507 + crRightHandSideBoundedVector19*crRightHandSideBoundedVector509 + crRightHandSideBoundedVector410*crRightHandSideBoundedVector656 + crRightHandSideBoundedVector491*crRightHandSideBoundedVector630 - crRightHandSideBoundedVector511 + crRightHandSideBoundedVector513 - crRightHandSideBoundedVector657*crRightHandSideBoundedVector662 - crRightHandSideBoundedVector663;
const double crRightHandSideBoundedVector665 =             crRightHandSideBoundedVector426*crRightHandSideBoundedVector76;
const double crRightHandSideBoundedVector666 =             crRightHandSideBoundedVector450*crRightHandSideBoundedVector76;
const double crRightHandSideBoundedVector667 =             crRightHandSideBoundedVector470*crRightHandSideBoundedVector76;
const double crRightHandSideBoundedVector668 =             crRightHandSideBoundedVector162*crRightHandSideBoundedVector511 + crRightHandSideBoundedVector162*crRightHandSideBoundedVector663 + crRightHandSideBoundedVector329*crRightHandSideBoundedVector655 + 0.16666666666666666*crRightHandSideBoundedVector456 + 0.16666666666666666*crRightHandSideBoundedVector457 + 0.16666666666666666*crRightHandSideBoundedVector458 - 0.16666666666666666*crRightHandSideBoundedVector459 - 0.16666666666666666*crRightHandSideBoundedVector460 - 0.16666666666666666*crRightHandSideBoundedVector465 - 0.16666666666666666*crRightHandSideBoundedVector466 + crRightHandSideBoundedVector497*crRightHandSideBoundedVector660;
const double crRightHandSideBoundedVector669 =             crRightHandSideBoundedVector130*crRightHandSideBoundedVector500 + crRightHandSideBoundedVector130*crRightHandSideBoundedVector652 + crRightHandSideBoundedVector329*crRightHandSideBoundedVector644 + 0.16666666666666666*crRightHandSideBoundedVector434 + 0.16666666666666666*crRightHandSideBoundedVector435 + 0.16666666666666666*crRightHandSideBoundedVector436 - 0.16666666666666666*crRightHandSideBoundedVector437 - 0.16666666666666666*crRightHandSideBoundedVector438 - 0.16666666666666666*crRightHandSideBoundedVector444 - 0.16666666666666666*crRightHandSideBoundedVector445 + crRightHandSideBoundedVector497*crRightHandSideBoundedVector649;
const double crRightHandSideBoundedVector670 =             0.99999999999999989*crRightHandSideBoundedVector52 + 0.99999999999999989*crRightHandSideBoundedVector53 + 0.99999999999999989*crRightHandSideBoundedVector54 + 0.99999999999999989*crRightHandSideBoundedVector72 + 0.99999999999999989*crRightHandSideBoundedVector73 + 0.99999999999999989*crRightHandSideBoundedVector74;
const double crRightHandSideBoundedVector671 =             DN_DX_1_1*crRightHandSideBoundedVector60;
const double crRightHandSideBoundedVector672 =             DN_DX_1_0*crRightHandSideBoundedVector70;
const double crRightHandSideBoundedVector673 =             0.37500000000000006*crRightHandSideBoundedVector10;
const double crRightHandSideBoundedVector674 =             0.37500000000000006*crRightHandSideBoundedVector398;
const double crRightHandSideBoundedVector675 =             crRightHandSideBoundedVector263*crRightHandSideBoundedVector673 - crRightHandSideBoundedVector265*crRightHandSideBoundedVector674 - crRightHandSideBoundedVector271*crRightHandSideBoundedVector3 + crRightHandSideBoundedVector272*crRightHandSideBoundedVector3;
const double crRightHandSideBoundedVector676 =             DN_DX_1_1*crRightHandSideBoundedVector127;
const double crRightHandSideBoundedVector677 =             0.66666666666666663*crRightHandSideBoundedVector96;
const double crRightHandSideBoundedVector678 =             DN_DX_1_0*crRightHandSideBoundedVector131;
const double crRightHandSideBoundedVector679 =             1.5000000000000002*crRightHandSideBoundedVector100;
const double crRightHandSideBoundedVector680 =             crRightHandSideBoundedVector100*crRightHandSideBoundedVector266;
const double crRightHandSideBoundedVector681 =             DN_DX_1_1*crRightHandSideBoundedVector160;
const double crRightHandSideBoundedVector682 =             DN_DX_1_0*crRightHandSideBoundedVector163;
const double crRightHandSideBoundedVector683 =             0.16666666666666666*crRightHandSideBoundedVector289;
const double crRightHandSideBoundedVector684 =             crRightHandSideBoundedVector59*crRightHandSideBoundedVector673;
const double crRightHandSideBoundedVector685 =             crRightHandSideBoundedVector673*crRightHandSideBoundedVector69;
const double crRightHandSideBoundedVector686 =             1.125*crRightHandSideBoundedVector289;
const double crRightHandSideBoundedVector687 =             crRightHandSideBoundedVector297*crRightHandSideBoundedVector686;
const double crRightHandSideBoundedVector688 =             crRightHandSideBoundedVector685*crRightHandSideBoundedVector83;
const double crRightHandSideBoundedVector689 =             0.1111111111111111*f_ext(0,0);
const double crRightHandSideBoundedVector690 =             crRightHandSideBoundedVector334 + crRightHandSideBoundedVector355 + crRightHandSideBoundedVector689;
const double crRightHandSideBoundedVector691 =             -crRightHandSideBoundedVector293*crRightHandSideBoundedVector683 + crRightHandSideBoundedVector302*crRightHandSideBoundedVector329 - crRightHandSideBoundedVector55*crRightHandSideBoundedVector684 - crRightHandSideBoundedVector65*crRightHandSideBoundedVector685 + crRightHandSideBoundedVector687 + crRightHandSideBoundedVector688 + crRightHandSideBoundedVector690;
const double crRightHandSideBoundedVector692 =             crRightHandSideBoundedVector126*crRightHandSideBoundedVector679;
const double crRightHandSideBoundedVector693 =             crRightHandSideBoundedVector130*crRightHandSideBoundedVector679;
const double crRightHandSideBoundedVector694 =             4.5*crRightHandSideBoundedVector318;
const double crRightHandSideBoundedVector695 =             crRightHandSideBoundedVector327*crRightHandSideBoundedVector694;
const double crRightHandSideBoundedVector696 =             crRightHandSideBoundedVector299*crRightHandSideBoundedVector300;
const double crRightHandSideBoundedVector697 =             crRightHandSideBoundedVector266*crRightHandSideBoundedVector499;
const double crRightHandSideBoundedVector698 =             crRightHandSideBoundedVector697*crRightHandSideBoundedVector82;
const double crRightHandSideBoundedVector699 =             crRightHandSideBoundedVector309 + crRightHandSideBoundedVector689 + 0.44444444444444442*f_ext(1,0);
const double crRightHandSideBoundedVector700 =             DN_DX_1_0*crRightHandSideBoundedVector60;
const double crRightHandSideBoundedVector701 =             DN_DX_1_1*crRightHandSideBoundedVector70;
const double crRightHandSideBoundedVector702 =             crRightHandSideBoundedVector3*crRightHandSideBoundedVector370 - crRightHandSideBoundedVector673*crRightHandSideBoundedVector86;
const double crRightHandSideBoundedVector703 =             crRightHandSideBoundedVector701 + crRightHandSideBoundedVector702;
const double crRightHandSideBoundedVector704 =             crRightHandSideBoundedVector183*crRightHandSideBoundedVector374 - crRightHandSideBoundedVector329*crRightHandSideBoundedVector360;
const double crRightHandSideBoundedVector705 =             crRightHandSideBoundedVector359*crRightHandSideBoundedVector96;
const double crRightHandSideBoundedVector706 =             DN_DX_1_0*crRightHandSideBoundedVector127;
const double crRightHandSideBoundedVector707 =             crRightHandSideBoundedVector100*crRightHandSideBoundedVector300;
const double crRightHandSideBoundedVector708 =             DN_DX_1_1*crRightHandSideBoundedVector131;
const double crRightHandSideBoundedVector709 =             crRightHandSideBoundedVector100*crRightHandSideBoundedVector563;
const double crRightHandSideBoundedVector710 =             crRightHandSideBoundedVector362*crRightHandSideBoundedVector96 + crRightHandSideBoundedVector708 - crRightHandSideBoundedVector709;
const double crRightHandSideBoundedVector711 =             DN_DX_1_0*crRightHandSideBoundedVector160;
const double crRightHandSideBoundedVector712 =             DN_DX_1_1*crRightHandSideBoundedVector163;
const double crRightHandSideBoundedVector713 =             crRightHandSideBoundedVector379 + crRightHandSideBoundedVector712;
const double crRightHandSideBoundedVector714 =             DN_DX_1_0*crRightHandSideBoundedVector262;
const double crRightHandSideBoundedVector715 =             crRightHandSideBoundedVector272*crRightHandSideBoundedVector70 + crRightHandSideBoundedVector373*crRightHandSideBoundedVector60 + 0.16666666666666666*crRightHandSideBoundedVector45 - 0.16666666666666666*crRightHandSideBoundedVector51 - 0.16666666666666666*crRightHandSideBoundedVector61 + crRightHandSideBoundedVector684*crRightHandSideBoundedVector86 - 0.16666666666666666*crRightHandSideBoundedVector71 - 0.99999999999999989*crRightHandSideBoundedVector90;
const double crRightHandSideBoundedVector716 =             crRightHandSideBoundedVector482*crRightHandSideBoundedVector686;
const double crRightHandSideBoundedVector717 =             crRightHandSideBoundedVector59*crRightHandSideBoundedVector674;
const double crRightHandSideBoundedVector718 =             crRightHandSideBoundedVector65*crRightHandSideBoundedVector717;
const double crRightHandSideBoundedVector719 =             0.1111111111111111*f_ext(0,1);
const double crRightHandSideBoundedVector720 =             crRightHandSideBoundedVector502 + crRightHandSideBoundedVector512 + crRightHandSideBoundedVector719;
const double crRightHandSideBoundedVector721 =             -crRightHandSideBoundedVector481*crRightHandSideBoundedVector683 + crRightHandSideBoundedVector485*crRightHandSideBoundedVector497 - crRightHandSideBoundedVector684*crRightHandSideBoundedVector82 - crRightHandSideBoundedVector685*crRightHandSideBoundedVector75 + crRightHandSideBoundedVector716 + crRightHandSideBoundedVector718 + crRightHandSideBoundedVector720;
const double crRightHandSideBoundedVector722 =             crRightHandSideBoundedVector495*crRightHandSideBoundedVector694;
const double crRightHandSideBoundedVector723 =             crRightHandSideBoundedVector300*crRightHandSideBoundedVector484;
const double crRightHandSideBoundedVector724 =             crRightHandSideBoundedVector266*crRightHandSideBoundedVector331;
const double crRightHandSideBoundedVector725 =             crRightHandSideBoundedVector65*crRightHandSideBoundedVector724;
const double crRightHandSideBoundedVector726 =             crRightHandSideBoundedVector489 + crRightHandSideBoundedVector719 + 0.44444444444444442*f_ext(1,1);
const double crRightHandSideBoundedVector727 =             -crRightHandSideBoundedVector182*crRightHandSideBoundedVector673 + crRightHandSideBoundedVector3*crRightHandSideBoundedVector520;
const double crRightHandSideBoundedVector728 =             crRightHandSideBoundedVector700 + crRightHandSideBoundedVector727;
const double crRightHandSideBoundedVector729 =             0.16666666666666666*crRightHandSideBoundedVector301*crRightHandSideBoundedVector86 - crRightHandSideBoundedVector360*crRightHandSideBoundedVector497;
const double crRightHandSideBoundedVector730 =             crRightHandSideBoundedVector263*crRightHandSideBoundedVector674 - crRightHandSideBoundedVector265*crRightHandSideBoundedVector673 + crRightHandSideBoundedVector3*crRightHandSideBoundedVector525 - crRightHandSideBoundedVector3*crRightHandSideBoundedVector526;
const double crRightHandSideBoundedVector731 =             crRightHandSideBoundedVector100*crRightHandSideBoundedVector595;
const double crRightHandSideBoundedVector732 =             crRightHandSideBoundedVector515*crRightHandSideBoundedVector96 + crRightHandSideBoundedVector706 - crRightHandSideBoundedVector731;
const double crRightHandSideBoundedVector733 =             crRightHandSideBoundedVector528 + crRightHandSideBoundedVector711;
const double crRightHandSideBoundedVector734 =             DN_DX_1_1*crRightHandSideBoundedVector262;
const double crRightHandSideBoundedVector735 =             0.16666666666666666*crRightHandSideBoundedVector172 - 0.16666666666666666*crRightHandSideBoundedVector175 - 0.16666666666666666*crRightHandSideBoundedVector176 - 0.16666666666666666*crRightHandSideBoundedVector177 + crRightHandSideBoundedVector182*crRightHandSideBoundedVector685 - 0.99999999999999989*crRightHandSideBoundedVector185 + crRightHandSideBoundedVector523*crRightHandSideBoundedVector70 + crRightHandSideBoundedVector526*crRightHandSideBoundedVector60;
const double crRightHandSideBoundedVector736 =             crRightHandSideBoundedVector418*crRightHandSideBoundedVector69;
const double crRightHandSideBoundedVector737 =             0.1111111111111111*r_ext[0];
const double crRightHandSideBoundedVector738 =             crRightHandSideBoundedVector609*crRightHandSideBoundedVector683;
const double crRightHandSideBoundedVector739 =             crRightHandSideBoundedVector603*crRightHandSideBoundedVector686;
const double crRightHandSideBoundedVector740 =             crRightHandSideBoundedVector179*crRightHandSideBoundedVector739 - crRightHandSideBoundedVector182*crRightHandSideBoundedVector738 + crRightHandSideBoundedVector329*crRightHandSideBoundedVector612 - crRightHandSideBoundedVector408*crRightHandSideBoundedVector684 - crRightHandSideBoundedVector410*crRightHandSideBoundedVector685 + crRightHandSideBoundedVector497*crRightHandSideBoundedVector613 + crRightHandSideBoundedVector616 + crRightHandSideBoundedVector627 + crRightHandSideBoundedVector737 - crRightHandSideBoundedVector738*crRightHandSideBoundedVector86 + crRightHandSideBoundedVector739*crRightHandSideBoundedVector83;
const double crRightHandSideBoundedVector741 =             crRightHandSideBoundedVector623*crRightHandSideBoundedVector694;
const double crRightHandSideBoundedVector742 =             crRightHandSideBoundedVector130*crRightHandSideBoundedVector443;
const double crRightHandSideBoundedVector743 =             crRightHandSideBoundedVector162*crRightHandSideBoundedVector464;
const double crRightHandSideBoundedVector744 =             0.16666666666666666*crRightHandSideBoundedVector3;
const double crRightHandSideBoundedVector745 =             1.1250000000000002*crRightHandSideBoundedVector398;
const double crRightHandSideBoundedVector746 =             crRightHandSideBoundedVector674*crRightHandSideBoundedVector69;
const double crRightHandSideBoundedVector747 =             crRightHandSideBoundedVector19*crRightHandSideBoundedVector687 + crRightHandSideBoundedVector317*crRightHandSideBoundedVector612 + crRightHandSideBoundedVector408*crRightHandSideBoundedVector744 - crRightHandSideBoundedVector55*crRightHandSideBoundedVector717 - crRightHandSideBoundedVector638*crRightHandSideBoundedVector745 - crRightHandSideBoundedVector65*crRightHandSideBoundedVector746 - crRightHandSideBoundedVector688 + crRightHandSideBoundedVector690;
const double crRightHandSideBoundedVector748 =             crRightHandSideBoundedVector684*crRightHandSideBoundedVector83;
const double crRightHandSideBoundedVector749 =             crRightHandSideBoundedVector19*crRightHandSideBoundedVector716 + crRightHandSideBoundedVector410*crRightHandSideBoundedVector744 + crRightHandSideBoundedVector491*crRightHandSideBoundedVector613 - crRightHandSideBoundedVector641*crRightHandSideBoundedVector745 - crRightHandSideBoundedVector718 + crRightHandSideBoundedVector720 - crRightHandSideBoundedVector746*crRightHandSideBoundedVector75 - crRightHandSideBoundedVector748;
const double crRightHandSideBoundedVector750 =             4.5000000000000009*crRightHandSideBoundedVector428;
const double crRightHandSideBoundedVector751 =             crRightHandSideBoundedVector724*crRightHandSideBoundedVector82;
const double crRightHandSideBoundedVector752 =             crRightHandSideBoundedVector329*crRightHandSideBoundedVector541 + 0.16666666666666666*crRightHandSideBoundedVector405 + 0.16666666666666666*crRightHandSideBoundedVector406 + 0.16666666666666666*crRightHandSideBoundedVector407 - 0.16666666666666666*crRightHandSideBoundedVector409 - 0.16666666666666666*crRightHandSideBoundedVector411 - 0.16666666666666666*crRightHandSideBoundedVector419 - 0.16666666666666666*crRightHandSideBoundedVector420 + crRightHandSideBoundedVector497*crRightHandSideBoundedVector542 + crRightHandSideBoundedVector69*crRightHandSideBoundedVector718 + crRightHandSideBoundedVector69*crRightHandSideBoundedVector748;
const double crRightHandSideBoundedVector753 =             DN_DX_2_1*crRightHandSideBoundedVector60;
const double crRightHandSideBoundedVector754 =             DN_DX_2_0*crRightHandSideBoundedVector70;
const double crRightHandSideBoundedVector755 =             DN_DX_2_1*crRightHandSideBoundedVector127;
const double crRightHandSideBoundedVector756 =             DN_DX_2_0*crRightHandSideBoundedVector131;
const double crRightHandSideBoundedVector757 =             DN_DX_2_1*crRightHandSideBoundedVector160;
const double crRightHandSideBoundedVector758 =             0.66666666666666663*crRightHandSideBoundedVector139;
const double crRightHandSideBoundedVector759 =             DN_DX_2_0*crRightHandSideBoundedVector163;
const double crRightHandSideBoundedVector760 =             1.5000000000000002*crRightHandSideBoundedVector141;
const double crRightHandSideBoundedVector761 =             crRightHandSideBoundedVector141*crRightHandSideBoundedVector266;
const double crRightHandSideBoundedVector762 =             crRightHandSideBoundedVector159*crRightHandSideBoundedVector760;
const double crRightHandSideBoundedVector763 =             crRightHandSideBoundedVector162*crRightHandSideBoundedVector760;
const double crRightHandSideBoundedVector764 =             4.5*crRightHandSideBoundedVector342;
const double crRightHandSideBoundedVector765 =             crRightHandSideBoundedVector351*crRightHandSideBoundedVector764;
const double crRightHandSideBoundedVector766 =             crRightHandSideBoundedVector266*crRightHandSideBoundedVector510;
const double crRightHandSideBoundedVector767 =             crRightHandSideBoundedVector766*crRightHandSideBoundedVector82;
const double crRightHandSideBoundedVector768 =             crRightHandSideBoundedVector308 + crRightHandSideBoundedVector689 + 0.44444444444444442*f_ext(2,0);
const double crRightHandSideBoundedVector769 =             DN_DX_2_0*crRightHandSideBoundedVector60;
const double crRightHandSideBoundedVector770 =             DN_DX_2_1*crRightHandSideBoundedVector70;
const double crRightHandSideBoundedVector771 =             crRightHandSideBoundedVector702 + crRightHandSideBoundedVector770;
const double crRightHandSideBoundedVector772 =             DN_DX_2_0*crRightHandSideBoundedVector127;
const double crRightHandSideBoundedVector773 =             DN_DX_2_1*crRightHandSideBoundedVector131;
const double crRightHandSideBoundedVector774 =             crRightHandSideBoundedVector371 + crRightHandSideBoundedVector773;
const double crRightHandSideBoundedVector775 =             crRightHandSideBoundedVector139*crRightHandSideBoundedVector359;
const double crRightHandSideBoundedVector776 =             DN_DX_2_0*crRightHandSideBoundedVector160;
const double crRightHandSideBoundedVector777 =             crRightHandSideBoundedVector141*crRightHandSideBoundedVector300;
const double crRightHandSideBoundedVector778 =             DN_DX_2_1*crRightHandSideBoundedVector163;
const double crRightHandSideBoundedVector779 =             crRightHandSideBoundedVector141*crRightHandSideBoundedVector578;
const double crRightHandSideBoundedVector780 =             crRightHandSideBoundedVector139*crRightHandSideBoundedVector362 + crRightHandSideBoundedVector778 - crRightHandSideBoundedVector779;
const double crRightHandSideBoundedVector781 =             DN_DX_2_0*crRightHandSideBoundedVector262;
const double crRightHandSideBoundedVector782 =             crRightHandSideBoundedVector508*crRightHandSideBoundedVector764;
const double crRightHandSideBoundedVector783 =             crRightHandSideBoundedVector266*crRightHandSideBoundedVector353;
const double crRightHandSideBoundedVector784 =             crRightHandSideBoundedVector65*crRightHandSideBoundedVector783;
const double crRightHandSideBoundedVector785 =             crRightHandSideBoundedVector488 + crRightHandSideBoundedVector719 + 0.44444444444444442*f_ext(2,1);
const double crRightHandSideBoundedVector786 =             crRightHandSideBoundedVector727 + crRightHandSideBoundedVector769;
const double crRightHandSideBoundedVector787 =             crRightHandSideBoundedVector521 + crRightHandSideBoundedVector772;
const double crRightHandSideBoundedVector788 =             crRightHandSideBoundedVector141*crRightHandSideBoundedVector598;
const double crRightHandSideBoundedVector789 =             crRightHandSideBoundedVector139*crRightHandSideBoundedVector515 + crRightHandSideBoundedVector776 - crRightHandSideBoundedVector788;
const double crRightHandSideBoundedVector790 =             DN_DX_2_1*crRightHandSideBoundedVector262;
const double crRightHandSideBoundedVector791 =             crRightHandSideBoundedVector633*crRightHandSideBoundedVector764;
const double crRightHandSideBoundedVector792 =             4.5000000000000009*crRightHandSideBoundedVector451;
const double crRightHandSideBoundedVector793 =             crRightHandSideBoundedVector783*crRightHandSideBoundedVector82;
            rRightHandSideBoundedVector[0]=-1.0*DN_DX_0_0*crRightHandSideBoundedVector137 - 1.0*DN_DX_0_0*crRightHandSideBoundedVector168 - 1.0*DN_DX_0_0*crRightHandSideBoundedVector93 - 1.0*DN_DX_0_1*crRightHandSideBoundedVector188 - 1.0*DN_DX_0_1*crRightHandSideBoundedVector201 - 1.0*DN_DX_0_1*crRightHandSideBoundedVector212 - 1.0*crRightHandSideBoundedVector213;
            rRightHandSideBoundedVector[1]=-DN_DX_0_0*crRightHandSideBoundedVector233 - DN_DX_0_0*crRightHandSideBoundedVector245 - DN_DX_0_0*crRightHandSideBoundedVector257 + crRightHandSideBoundedVector268*(-crRightHandSideBoundedVector258 - crRightHandSideBoundedVector259*crRightHandSideBoundedVector3 + crRightHandSideBoundedVector260*crRightHandSideBoundedVector83 + crRightHandSideBoundedVector261*crRightHandSideBoundedVector262 + crRightHandSideBoundedVector263*crRightHandSideBoundedVector264 - crRightHandSideBoundedVector265*crRightHandSideBoundedVector267) + crRightHandSideBoundedVector279*(crRightHandSideBoundedVector262*crRightHandSideBoundedVector270 - crRightHandSideBoundedVector269 + crRightHandSideBoundedVector278) + crRightHandSideBoundedVector288*(crRightHandSideBoundedVector262*crRightHandSideBoundedVector281 - crRightHandSideBoundedVector280 + crRightHandSideBoundedVector287) - crRightHandSideBoundedVector305*crRightHandSideBoundedVector386 - crRightHandSideBoundedVector313*(DN_DX_0_0*crRightHandSideBoundedVector50 - crRightHandSideBoundedVector290*crRightHandSideBoundedVector293 - crRightHandSideBoundedVector294*crRightHandSideBoundedVector55 - crRightHandSideBoundedVector295*crRightHandSideBoundedVector65 + crRightHandSideBoundedVector298 + crRightHandSideBoundedVector299*crRightHandSideBoundedVector302 + crRightHandSideBoundedVector304 - crRightHandSideBoundedVector307 + crRightHandSideBoundedVector310) - crRightHandSideBoundedVector314*crRightHandSideBoundedVector389 - crRightHandSideBoundedVector338*(DN_DX_0_0*crRightHandSideBoundedVector123 - crRightHandSideBoundedVector316 + crRightHandSideBoundedVector336) - crRightHandSideBoundedVector339*crRightHandSideBoundedVector392 - crRightHandSideBoundedVector358*(DN_DX_0_0*crRightHandSideBoundedVector157 - crRightHandSideBoundedVector341 + crRightHandSideBoundedVector357) + crRightHandSideBoundedVector365*crRightHandSideBoundedVector59 - crRightHandSideBoundedVector367*(crRightHandSideBoundedVector221*crRightHandSideBoundedVector301 - crRightHandSideBoundedVector299*crRightHandSideBoundedVector360 - crRightHandSideBoundedVector359*crRightHandSideBoundedVector361 + crRightHandSideBoundedVector366) - crRightHandSideBoundedVector376*(-crRightHandSideBoundedVector359*crRightHandSideBoundedVector368 + crRightHandSideBoundedVector372 + crRightHandSideBoundedVector375) - crRightHandSideBoundedVector382*(-crRightHandSideBoundedVector359*crRightHandSideBoundedVector377 + crRightHandSideBoundedVector380 + crRightHandSideBoundedVector381) + crRightHandSideBoundedVector393*crRightHandSideBoundedVector70 + crRightHandSideBoundedVector394*crRightHandSideBoundedVector60 - crRightHandSideBoundedVector395*crRightHandSideBoundedVector426 - crRightHandSideBoundedVector395*crRightHandSideBoundedVector450 - crRightHandSideBoundedVector395*crRightHandSideBoundedVector470 + 0.66666666666666663*crRightHandSideBoundedVector45 + crRightHandSideBoundedVector471 + crRightHandSideBoundedVector472 - 0.66666666666666663*crRightHandSideBoundedVector51 - 0.66666666666666663*crRightHandSideBoundedVector61 - 0.66666666666666663*crRightHandSideBoundedVector71 - 1.0*crRightHandSideBoundedVector90;
            rRightHandSideBoundedVector[2]=-DN_DX_0_0*crRightHandSideBoundedVector532 - DN_DX_0_0*crRightHandSideBoundedVector533 - DN_DX_0_0*crRightHandSideBoundedVector534 - DN_DX_0_1*crRightHandSideBoundedVector477 - DN_DX_0_1*crRightHandSideBoundedVector478 - DN_DX_0_1*crRightHandSideBoundedVector479 + 0.66666666666666663*crRightHandSideBoundedVector172 - 0.66666666666666663*crRightHandSideBoundedVector175 - 0.66666666666666663*crRightHandSideBoundedVector176 - 0.66666666666666663*crRightHandSideBoundedVector177 - 1.0*crRightHandSideBoundedVector185 - crRightHandSideBoundedVector268*(crRightHandSideBoundedVector226*crRightHandSideBoundedVector301 - crRightHandSideBoundedVector359*crRightHandSideBoundedVector363 - crRightHandSideBoundedVector360*crRightHandSideBoundedVector484 + crRightHandSideBoundedVector518) - crRightHandSideBoundedVector279*(-crRightHandSideBoundedVector359*crRightHandSideBoundedVector369 + crRightHandSideBoundedVector522 + crRightHandSideBoundedVector524) - crRightHandSideBoundedVector288*(-crRightHandSideBoundedVector359*crRightHandSideBoundedVector378 + crRightHandSideBoundedVector529 + crRightHandSideBoundedVector530) - crRightHandSideBoundedVector313*(-DN_DX_0_0*crRightHandSideBoundedVector421 + DN_DX_0_1*crRightHandSideBoundedVector174 - crRightHandSideBoundedVector290*crRightHandSideBoundedVector481 - crRightHandSideBoundedVector294*crRightHandSideBoundedVector82 - crRightHandSideBoundedVector295*crRightHandSideBoundedVector75 + crRightHandSideBoundedVector483 + crRightHandSideBoundedVector484*crRightHandSideBoundedVector485 + crRightHandSideBoundedVector487 + crRightHandSideBoundedVector490) - crRightHandSideBoundedVector338*(-DN_DX_0_0*crRightHandSideBoundedVector446 + DN_DX_0_1*crRightHandSideBoundedVector193 + crRightHandSideBoundedVector504) - crRightHandSideBoundedVector358*(-DN_DX_0_0*crRightHandSideBoundedVector467 + DN_DX_0_1*crRightHandSideBoundedVector205 + crRightHandSideBoundedVector514) - crRightHandSideBoundedVector367*(-crRightHandSideBoundedVector179*crRightHandSideBoundedVector260 - crRightHandSideBoundedVector258*crRightHandSideBoundedVector262 + crRightHandSideBoundedVector261 + crRightHandSideBoundedVector263*crRightHandSideBoundedVector267 - crRightHandSideBoundedVector264*crRightHandSideBoundedVector265 + crRightHandSideBoundedVector3*crRightHandSideBoundedVector519) - crRightHandSideBoundedVector376*(-crRightHandSideBoundedVector262*crRightHandSideBoundedVector269 + crRightHandSideBoundedVector270 + crRightHandSideBoundedVector527) - crRightHandSideBoundedVector382*(-crRightHandSideBoundedVector262*crRightHandSideBoundedVector280 + crRightHandSideBoundedVector281 + crRightHandSideBoundedVector531) - crRightHandSideBoundedVector426*crRightHandSideBoundedVector537 - crRightHandSideBoundedVector450*crRightHandSideBoundedVector537 - crRightHandSideBoundedVector470*crRightHandSideBoundedVector537 + crRightHandSideBoundedVector517*crRightHandSideBoundedVector69 + crRightHandSideBoundedVector535*crRightHandSideBoundedVector60 + crRightHandSideBoundedVector536*crRightHandSideBoundedVector70 + crRightHandSideBoundedVector538 + crRightHandSideBoundedVector539;
            rRightHandSideBoundedVector[3]=-DN_DX_0_0*crRightHandSideBoundedVector560 - DN_DX_0_0*crRightHandSideBoundedVector575 - DN_DX_0_0*crRightHandSideBoundedVector590 - DN_DX_0_1*crRightHandSideBoundedVector594 - DN_DX_0_1*crRightHandSideBoundedVector597 - DN_DX_0_1*crRightHandSideBoundedVector600 - crRightHandSideBoundedVector221*crRightHandSideBoundedVector418 - crRightHandSideBoundedVector226*crRightHandSideBoundedVector418 - crRightHandSideBoundedVector268*(-DN_DX_0_0*crRightHandSideBoundedVector642 - DN_DX_0_1*crRightHandSideBoundedVector542 + crRightHandSideBoundedVector19*crRightHandSideBoundedVector483 + crRightHandSideBoundedVector260*crRightHandSideBoundedVector410 - crRightHandSideBoundedVector303*crRightHandSideBoundedVector75 - crRightHandSideBoundedVector486*crRightHandSideBoundedVector82 - crRightHandSideBoundedVector487 + crRightHandSideBoundedVector490 + crRightHandSideBoundedVector613*crRightHandSideBoundedVector640 - crRightHandSideBoundedVector637*crRightHandSideBoundedVector641) - crRightHandSideBoundedVector279*(-DN_DX_0_0*crRightHandSideBoundedVector650 - DN_DX_0_1*crRightHandSideBoundedVector649 + crRightHandSideBoundedVector653) - crRightHandSideBoundedVector288*(-DN_DX_0_0*crRightHandSideBoundedVector661 - DN_DX_0_1*crRightHandSideBoundedVector660 + crRightHandSideBoundedVector664) + crRightHandSideBoundedVector299*crRightHandSideBoundedVector541 + crRightHandSideBoundedVector304*crRightHandSideBoundedVector59 - crRightHandSideBoundedVector313*(DN_DX_0_0*crRightHandSideBoundedVector605 + crRightHandSideBoundedVector179*crRightHandSideBoundedVector604 - crRightHandSideBoundedVector221*crRightHandSideBoundedVector610 - crRightHandSideBoundedVector226*crRightHandSideBoundedVector610 - crRightHandSideBoundedVector294*crRightHandSideBoundedVector408 - crRightHandSideBoundedVector295*crRightHandSideBoundedVector410 + crRightHandSideBoundedVector299*crRightHandSideBoundedVector612 + crRightHandSideBoundedVector306*crRightHandSideBoundedVector417 + crRightHandSideBoundedVector484*crRightHandSideBoundedVector613 + crRightHandSideBoundedVector601 + crRightHandSideBoundedVector602 + crRightHandSideBoundedVector604*crRightHandSideBoundedVector83 + 0.44444444444444442*r_ext[0]) - crRightHandSideBoundedVector338*(DN_DX_0_0*crRightHandSideBoundedVector614 + crRightHandSideBoundedVector315*crRightHandSideBoundedVector442 + crRightHandSideBoundedVector625) - crRightHandSideBoundedVector358*(DN_DX_0_0*crRightHandSideBoundedVector626 + crRightHandSideBoundedVector340*crRightHandSideBoundedVector463 + crRightHandSideBoundedVector635) - crRightHandSideBoundedVector367*(-DN_DX_0_0*crRightHandSideBoundedVector541 + crRightHandSideBoundedVector19*crRightHandSideBoundedVector298 - crRightHandSideBoundedVector19*crRightHandSideBoundedVector307 + crRightHandSideBoundedVector260*crRightHandSideBoundedVector408 - crRightHandSideBoundedVector304 + crRightHandSideBoundedVector310 - crRightHandSideBoundedVector486*crRightHandSideBoundedVector55 + crRightHandSideBoundedVector612*crRightHandSideBoundedVector636 - crRightHandSideBoundedVector637*crRightHandSideBoundedVector638 - crRightHandSideBoundedVector639) - crRightHandSideBoundedVector376*(-DN_DX_0_0*crRightHandSideBoundedVector644 - crRightHandSideBoundedVector19*crRightHandSideBoundedVector316 + crRightHandSideBoundedVector648) - crRightHandSideBoundedVector382*(-DN_DX_0_0*crRightHandSideBoundedVector655 - crRightHandSideBoundedVector19*crRightHandSideBoundedVector341 + crRightHandSideBoundedVector659) + 0.66666666666666663*crRightHandSideBoundedVector405 + 0.66666666666666663*crRightHandSideBoundedVector406 + 0.66666666666666663*crRightHandSideBoundedVector407 - 0.66666666666666663*crRightHandSideBoundedVector409 - 0.66666666666666663*crRightHandSideBoundedVector411 + crRightHandSideBoundedVector484*crRightHandSideBoundedVector542 + crRightHandSideBoundedVector59*crRightHandSideBoundedVector639 - crRightHandSideBoundedVector665*(crRightHandSideBoundedVector366 + crRightHandSideBoundedVector518) - crRightHandSideBoundedVector666*(crRightHandSideBoundedVector372 + crRightHandSideBoundedVector522) - crRightHandSideBoundedVector667*(crRightHandSideBoundedVector380 + crRightHandSideBoundedVector529) + crRightHandSideBoundedVector668 + crRightHandSideBoundedVector669;
            rRightHandSideBoundedVector[4]=-DN_DX_1_0*crRightHandSideBoundedVector367 - DN_DX_1_0*crRightHandSideBoundedVector376 - DN_DX_1_0*crRightHandSideBoundedVector382 - DN_DX_1_1*crRightHandSideBoundedVector268 - DN_DX_1_1*crRightHandSideBoundedVector279 - DN_DX_1_1*crRightHandSideBoundedVector288 - crRightHandSideBoundedVector670;
            rRightHandSideBoundedVector[5]=-DN_DX_1_0*crRightHandSideBoundedVector233 - DN_DX_1_0*crRightHandSideBoundedVector245 - DN_DX_1_0*crRightHandSideBoundedVector257 - DN_DX_1_1*crRightHandSideBoundedVector532 - DN_DX_1_1*crRightHandSideBoundedVector533 - DN_DX_1_1*crRightHandSideBoundedVector534 + 0.66666666666666663*crRightHandSideBoundedVector120 - 0.66666666666666663*crRightHandSideBoundedVector124 + crRightHandSideBoundedVector126*crRightHandSideBoundedVector709 + crRightHandSideBoundedVector127*crRightHandSideBoundedVector394 - 0.66666666666666663*crRightHandSideBoundedVector128 + crRightHandSideBoundedVector131*crRightHandSideBoundedVector393 - 0.66666666666666663*crRightHandSideBoundedVector132 + crRightHandSideBoundedVector268*(crRightHandSideBoundedVector262*crRightHandSideBoundedVector672 - crRightHandSideBoundedVector671 + crRightHandSideBoundedVector675) + crRightHandSideBoundedVector279*(-crRightHandSideBoundedVector259*crRightHandSideBoundedVector96 + crRightHandSideBoundedVector262*crRightHandSideBoundedVector678 + crRightHandSideBoundedVector273*crRightHandSideBoundedVector679 - crRightHandSideBoundedVector276*crRightHandSideBoundedVector680 - crRightHandSideBoundedVector676 + crRightHandSideBoundedVector677*crRightHandSideBoundedVector83) + crRightHandSideBoundedVector288*(crRightHandSideBoundedVector262*crRightHandSideBoundedVector682 + crRightHandSideBoundedVector287 - crRightHandSideBoundedVector681) - crRightHandSideBoundedVector313*(DN_DX_1_0*crRightHandSideBoundedVector50 - DN_DX_1_1*crRightHandSideBoundedVector421 + crRightHandSideBoundedVector691) - crRightHandSideBoundedVector338*(DN_DX_1_0*crRightHandSideBoundedVector123 - DN_DX_1_1*crRightHandSideBoundedVector446 - crRightHandSideBoundedVector321*crRightHandSideBoundedVector636 + crRightHandSideBoundedVector331*crRightHandSideBoundedVector696 - crRightHandSideBoundedVector55*crRightHandSideBoundedVector692 - crRightHandSideBoundedVector65*crRightHandSideBoundedVector693 + crRightHandSideBoundedVector695 + crRightHandSideBoundedVector698 + crRightHandSideBoundedVector699) - crRightHandSideBoundedVector358*(DN_DX_1_0*crRightHandSideBoundedVector157 - DN_DX_1_1*crRightHandSideBoundedVector467 + crRightHandSideBoundedVector357) - crRightHandSideBoundedVector367*(-crRightHandSideBoundedVector359*crRightHandSideBoundedVector700 + crRightHandSideBoundedVector703 + crRightHandSideBoundedVector704) - crRightHandSideBoundedVector376*(crRightHandSideBoundedVector241*crRightHandSideBoundedVector707 - crRightHandSideBoundedVector299*crRightHandSideBoundedVector705 - crRightHandSideBoundedVector359*crRightHandSideBoundedVector706 + crRightHandSideBoundedVector710) - crRightHandSideBoundedVector382*(-crRightHandSideBoundedVector359*crRightHandSideBoundedVector711 + crRightHandSideBoundedVector381 + crRightHandSideBoundedVector713) - crRightHandSideBoundedVector426*crRightHandSideBoundedVector714 - crRightHandSideBoundedVector450*crRightHandSideBoundedVector714 - crRightHandSideBoundedVector470*crRightHandSideBoundedVector714 + crRightHandSideBoundedVector471 + crRightHandSideBoundedVector715;
            rRightHandSideBoundedVector[6]=-DN_DX_1_0*crRightHandSideBoundedVector532 - DN_DX_1_0*crRightHandSideBoundedVector533 - DN_DX_1_0*crRightHandSideBoundedVector534 - DN_DX_1_1*crRightHandSideBoundedVector477 - DN_DX_1_1*crRightHandSideBoundedVector478 - DN_DX_1_1*crRightHandSideBoundedVector479 + crRightHandSideBoundedVector127*crRightHandSideBoundedVector535 + crRightHandSideBoundedVector130*crRightHandSideBoundedVector731 + crRightHandSideBoundedVector131*crRightHandSideBoundedVector536 + 0.66666666666666663*crRightHandSideBoundedVector191 - 0.66666666666666663*crRightHandSideBoundedVector194 - 0.66666666666666663*crRightHandSideBoundedVector195 - 0.66666666666666663*crRightHandSideBoundedVector196 - crRightHandSideBoundedVector268*(-crRightHandSideBoundedVector359*crRightHandSideBoundedVector701 + crRightHandSideBoundedVector728 + crRightHandSideBoundedVector729) - crRightHandSideBoundedVector279*(crRightHandSideBoundedVector243*crRightHandSideBoundedVector707 - crRightHandSideBoundedVector359*crRightHandSideBoundedVector708 - crRightHandSideBoundedVector484*crRightHandSideBoundedVector705 + crRightHandSideBoundedVector732) - crRightHandSideBoundedVector288*(-crRightHandSideBoundedVector359*crRightHandSideBoundedVector712 + crRightHandSideBoundedVector530 + crRightHandSideBoundedVector733) - crRightHandSideBoundedVector313*(-DN_DX_1_0*crRightHandSideBoundedVector421 + DN_DX_1_1*crRightHandSideBoundedVector174 + crRightHandSideBoundedVector721) - crRightHandSideBoundedVector338*(-DN_DX_1_0*crRightHandSideBoundedVector446 + DN_DX_1_1*crRightHandSideBoundedVector193 - crRightHandSideBoundedVector493*crRightHandSideBoundedVector640 + crRightHandSideBoundedVector499*crRightHandSideBoundedVector723 - crRightHandSideBoundedVector692*crRightHandSideBoundedVector82 - crRightHandSideBoundedVector693*crRightHandSideBoundedVector75 + crRightHandSideBoundedVector722 + crRightHandSideBoundedVector725 + crRightHandSideBoundedVector726) - crRightHandSideBoundedVector358*(-DN_DX_1_0*crRightHandSideBoundedVector467 + DN_DX_1_1*crRightHandSideBoundedVector205 + crRightHandSideBoundedVector514) - crRightHandSideBoundedVector367*(-crRightHandSideBoundedVector262*crRightHandSideBoundedVector671 + crRightHandSideBoundedVector672 + crRightHandSideBoundedVector730) - crRightHandSideBoundedVector376*(-crRightHandSideBoundedVector179*crRightHandSideBoundedVector677 - crRightHandSideBoundedVector262*crRightHandSideBoundedVector676 + crRightHandSideBoundedVector273*crRightHandSideBoundedVector680 - crRightHandSideBoundedVector276*crRightHandSideBoundedVector679 + crRightHandSideBoundedVector519*crRightHandSideBoundedVector96 + crRightHandSideBoundedVector678) - crRightHandSideBoundedVector382*(-crRightHandSideBoundedVector262*crRightHandSideBoundedVector681 + crRightHandSideBoundedVector531 + crRightHandSideBoundedVector682) - crRightHandSideBoundedVector426*crRightHandSideBoundedVector734 - crRightHandSideBoundedVector450*crRightHandSideBoundedVector734 - crRightHandSideBoundedVector470*crRightHandSideBoundedVector734 + crRightHandSideBoundedVector538 + crRightHandSideBoundedVector735;
            rRightHandSideBoundedVector[7]=-DN_DX_1_0*crRightHandSideBoundedVector560 - DN_DX_1_0*crRightHandSideBoundedVector575 - DN_DX_1_0*crRightHandSideBoundedVector590 - DN_DX_1_1*crRightHandSideBoundedVector594 - DN_DX_1_1*crRightHandSideBoundedVector597 - DN_DX_1_1*crRightHandSideBoundedVector600 + crRightHandSideBoundedVector130*crRightHandSideBoundedVector725 + crRightHandSideBoundedVector130*crRightHandSideBoundedVector751 - crRightHandSideBoundedVector241*crRightHandSideBoundedVector443 - crRightHandSideBoundedVector243*crRightHandSideBoundedVector443 - crRightHandSideBoundedVector268*(-DN_DX_1_0*crRightHandSideBoundedVector642 - DN_DX_1_1*crRightHandSideBoundedVector542 + crRightHandSideBoundedVector749) - crRightHandSideBoundedVector279*(-DN_DX_1_0*crRightHandSideBoundedVector650 - DN_DX_1_1*crRightHandSideBoundedVector649 + crRightHandSideBoundedVector19*crRightHandSideBoundedVector722 + crRightHandSideBoundedVector410*crRightHandSideBoundedVector677 + crRightHandSideBoundedVector619*crRightHandSideBoundedVector640 - crRightHandSideBoundedVector651*crRightHandSideBoundedVector750 - crRightHandSideBoundedVector697*crRightHandSideBoundedVector75 - crRightHandSideBoundedVector725 + crRightHandSideBoundedVector726 - crRightHandSideBoundedVector751) - crRightHandSideBoundedVector288*(-DN_DX_1_0*crRightHandSideBoundedVector661 - DN_DX_1_1*crRightHandSideBoundedVector660 + crRightHandSideBoundedVector664) + crRightHandSideBoundedVector299*crRightHandSideBoundedVector644 - crRightHandSideBoundedVector313*(DN_DX_1_0*crRightHandSideBoundedVector605 + DN_DX_1_1*crRightHandSideBoundedVector736 + crRightHandSideBoundedVector740) - crRightHandSideBoundedVector338*(DN_DX_1_0*crRightHandSideBoundedVector614 + DN_DX_1_1*crRightHandSideBoundedVector742 + crRightHandSideBoundedVector179*crRightHandSideBoundedVector741 - crRightHandSideBoundedVector241*crRightHandSideBoundedVector621 - crRightHandSideBoundedVector243*crRightHandSideBoundedVector621 + crRightHandSideBoundedVector299*crRightHandSideBoundedVector618 - crRightHandSideBoundedVector408*crRightHandSideBoundedVector692 - crRightHandSideBoundedVector410*crRightHandSideBoundedVector693 + crRightHandSideBoundedVector484*crRightHandSideBoundedVector619 + crRightHandSideBoundedVector602 + crRightHandSideBoundedVector737 + crRightHandSideBoundedVector741*crRightHandSideBoundedVector83 + 0.44444444444444442*r_ext[1]) - crRightHandSideBoundedVector358*(DN_DX_1_0*crRightHandSideBoundedVector626 + DN_DX_1_1*crRightHandSideBoundedVector743 + crRightHandSideBoundedVector635) - crRightHandSideBoundedVector367*(-DN_DX_1_0*crRightHandSideBoundedVector541 - DN_DX_1_1*crRightHandSideBoundedVector642 + crRightHandSideBoundedVector747) - crRightHandSideBoundedVector376*(-DN_DX_1_0*crRightHandSideBoundedVector644 - DN_DX_1_1*crRightHandSideBoundedVector650 + crRightHandSideBoundedVector19*crRightHandSideBoundedVector695 + crRightHandSideBoundedVector408*crRightHandSideBoundedVector677 - crRightHandSideBoundedVector55*crRightHandSideBoundedVector724 + crRightHandSideBoundedVector618*crRightHandSideBoundedVector636 - crRightHandSideBoundedVector647*crRightHandSideBoundedVector750 - crRightHandSideBoundedVector65*crRightHandSideBoundedVector697 - crRightHandSideBoundedVector698 + crRightHandSideBoundedVector699) - crRightHandSideBoundedVector382*(-DN_DX_1_0*crRightHandSideBoundedVector655 - DN_DX_1_1*crRightHandSideBoundedVector661 + crRightHandSideBoundedVector659) + 0.66666666666666663*crRightHandSideBoundedVector434 + 0.66666666666666663*crRightHandSideBoundedVector435 + 0.66666666666666663*crRightHandSideBoundedVector436 - 0.66666666666666663*crRightHandSideBoundedVector437 - 0.66666666666666663*crRightHandSideBoundedVector438 + crRightHandSideBoundedVector484*crRightHandSideBoundedVector649 - crRightHandSideBoundedVector665*(crRightHandSideBoundedVector703 + crRightHandSideBoundedVector728) - crRightHandSideBoundedVector666*(crRightHandSideBoundedVector710 + crRightHandSideBoundedVector732) - crRightHandSideBoundedVector667*(crRightHandSideBoundedVector713 + crRightHandSideBoundedVector733) + crRightHandSideBoundedVector668 + crRightHandSideBoundedVector752;
            rRightHandSideBoundedVector[8]=-DN_DX_2_0*crRightHandSideBoundedVector367 - DN_DX_2_0*crRightHandSideBoundedVector376 - DN_DX_2_0*crRightHandSideBoundedVector382 - DN_DX_2_1*crRightHandSideBoundedVector268 - DN_DX_2_1*crRightHandSideBoundedVector279 - DN_DX_2_1*crRightHandSideBoundedVector288 - crRightHandSideBoundedVector670;
            rRightHandSideBoundedVector[9]=-DN_DX_2_0*crRightHandSideBoundedVector233 - DN_DX_2_0*crRightHandSideBoundedVector245 - DN_DX_2_0*crRightHandSideBoundedVector257 - DN_DX_2_1*crRightHandSideBoundedVector532 - DN_DX_2_1*crRightHandSideBoundedVector533 - DN_DX_2_1*crRightHandSideBoundedVector534 + 0.66666666666666663*crRightHandSideBoundedVector154 - 0.66666666666666663*crRightHandSideBoundedVector158 + crRightHandSideBoundedVector159*crRightHandSideBoundedVector779 + crRightHandSideBoundedVector160*crRightHandSideBoundedVector394 - 0.66666666666666663*crRightHandSideBoundedVector161 + crRightHandSideBoundedVector163*crRightHandSideBoundedVector393 - 0.66666666666666663*crRightHandSideBoundedVector164 + crRightHandSideBoundedVector268*(crRightHandSideBoundedVector262*crRightHandSideBoundedVector754 + crRightHandSideBoundedVector675 - crRightHandSideBoundedVector753) + crRightHandSideBoundedVector279*(crRightHandSideBoundedVector262*crRightHandSideBoundedVector756 + crRightHandSideBoundedVector278 - crRightHandSideBoundedVector755) + crRightHandSideBoundedVector288*(-crRightHandSideBoundedVector139*crRightHandSideBoundedVector259 + crRightHandSideBoundedVector262*crRightHandSideBoundedVector759 + crRightHandSideBoundedVector282*crRightHandSideBoundedVector760 - crRightHandSideBoundedVector285*crRightHandSideBoundedVector761 - crRightHandSideBoundedVector757 + crRightHandSideBoundedVector758*crRightHandSideBoundedVector83) - crRightHandSideBoundedVector313*(DN_DX_2_0*crRightHandSideBoundedVector50 - DN_DX_2_1*crRightHandSideBoundedVector421 + crRightHandSideBoundedVector691) - crRightHandSideBoundedVector338*(DN_DX_2_0*crRightHandSideBoundedVector123 - DN_DX_2_1*crRightHandSideBoundedVector446 + crRightHandSideBoundedVector336) - crRightHandSideBoundedVector358*(DN_DX_2_0*crRightHandSideBoundedVector157 - DN_DX_2_1*crRightHandSideBoundedVector467 - crRightHandSideBoundedVector345*crRightHandSideBoundedVector636 + crRightHandSideBoundedVector353*crRightHandSideBoundedVector696 - crRightHandSideBoundedVector55*crRightHandSideBoundedVector762 - crRightHandSideBoundedVector65*crRightHandSideBoundedVector763 + crRightHandSideBoundedVector765 + crRightHandSideBoundedVector767 + crRightHandSideBoundedVector768) - crRightHandSideBoundedVector367*(-crRightHandSideBoundedVector359*crRightHandSideBoundedVector769 + crRightHandSideBoundedVector704 + crRightHandSideBoundedVector771) - crRightHandSideBoundedVector376*(-crRightHandSideBoundedVector359*crRightHandSideBoundedVector772 + crRightHandSideBoundedVector375 + crRightHandSideBoundedVector774) - crRightHandSideBoundedVector382*(crRightHandSideBoundedVector253*crRightHandSideBoundedVector777 - crRightHandSideBoundedVector299*crRightHandSideBoundedVector775 - crRightHandSideBoundedVector359*crRightHandSideBoundedVector776 + crRightHandSideBoundedVector780) - crRightHandSideBoundedVector426*crRightHandSideBoundedVector781 - crRightHandSideBoundedVector450*crRightHandSideBoundedVector781 - crRightHandSideBoundedVector470*crRightHandSideBoundedVector781 + crRightHandSideBoundedVector472 + crRightHandSideBoundedVector715;
            rRightHandSideBoundedVector[10]=-DN_DX_2_0*crRightHandSideBoundedVector532 - DN_DX_2_0*crRightHandSideBoundedVector533 - DN_DX_2_0*crRightHandSideBoundedVector534 - DN_DX_2_1*crRightHandSideBoundedVector477 - DN_DX_2_1*crRightHandSideBoundedVector478 - DN_DX_2_1*crRightHandSideBoundedVector479 + crRightHandSideBoundedVector160*crRightHandSideBoundedVector535 + crRightHandSideBoundedVector162*crRightHandSideBoundedVector788 + crRightHandSideBoundedVector163*crRightHandSideBoundedVector536 + 0.66666666666666663*crRightHandSideBoundedVector203 - 0.66666666666666663*crRightHandSideBoundedVector206 - 0.66666666666666663*crRightHandSideBoundedVector207 - 0.66666666666666663*crRightHandSideBoundedVector208 - crRightHandSideBoundedVector268*(-crRightHandSideBoundedVector359*crRightHandSideBoundedVector770 + crRightHandSideBoundedVector729 + crRightHandSideBoundedVector786) - crRightHandSideBoundedVector279*(-crRightHandSideBoundedVector359*crRightHandSideBoundedVector773 + crRightHandSideBoundedVector524 + crRightHandSideBoundedVector787) - crRightHandSideBoundedVector288*(crRightHandSideBoundedVector255*crRightHandSideBoundedVector777 - crRightHandSideBoundedVector359*crRightHandSideBoundedVector778 - crRightHandSideBoundedVector484*crRightHandSideBoundedVector775 + crRightHandSideBoundedVector789) - crRightHandSideBoundedVector313*(-DN_DX_2_0*crRightHandSideBoundedVector421 + DN_DX_2_1*crRightHandSideBoundedVector174 + crRightHandSideBoundedVector721) - crRightHandSideBoundedVector338*(-DN_DX_2_0*crRightHandSideBoundedVector446 + DN_DX_2_1*crRightHandSideBoundedVector193 + crRightHandSideBoundedVector504) - crRightHandSideBoundedVector358*(-DN_DX_2_0*crRightHandSideBoundedVector467 + DN_DX_2_1*crRightHandSideBoundedVector205 - crRightHandSideBoundedVector506*crRightHandSideBoundedVector640 + crRightHandSideBoundedVector510*crRightHandSideBoundedVector723 - crRightHandSideBoundedVector75*crRightHandSideBoundedVector763 - crRightHandSideBoundedVector762*crRightHandSideBoundedVector82 + crRightHandSideBoundedVector782 + crRightHandSideBoundedVector784 + crRightHandSideBoundedVector785) - crRightHandSideBoundedVector367*(-crRightHandSideBoundedVector262*crRightHandSideBoundedVector753 + crRightHandSideBoundedVector730 + crRightHandSideBoundedVector754) - crRightHandSideBoundedVector376*(-crRightHandSideBoundedVector262*crRightHandSideBoundedVector755 + crRightHandSideBoundedVector527 + crRightHandSideBoundedVector756) - crRightHandSideBoundedVector382*(crRightHandSideBoundedVector139*crRightHandSideBoundedVector519 - crRightHandSideBoundedVector179*crRightHandSideBoundedVector758 - crRightHandSideBoundedVector262*crRightHandSideBoundedVector757 + crRightHandSideBoundedVector282*crRightHandSideBoundedVector761 - crRightHandSideBoundedVector285*crRightHandSideBoundedVector760 + crRightHandSideBoundedVector759) - crRightHandSideBoundedVector426*crRightHandSideBoundedVector790 - crRightHandSideBoundedVector450*crRightHandSideBoundedVector790 - crRightHandSideBoundedVector470*crRightHandSideBoundedVector790 + crRightHandSideBoundedVector539 + crRightHandSideBoundedVector735;
            rRightHandSideBoundedVector[11]=-DN_DX_2_0*crRightHandSideBoundedVector560 - DN_DX_2_0*crRightHandSideBoundedVector575 - DN_DX_2_0*crRightHandSideBoundedVector590 - DN_DX_2_1*crRightHandSideBoundedVector594 - DN_DX_2_1*crRightHandSideBoundedVector597 - DN_DX_2_1*crRightHandSideBoundedVector600 + crRightHandSideBoundedVector162*crRightHandSideBoundedVector784 + crRightHandSideBoundedVector162*crRightHandSideBoundedVector793 - crRightHandSideBoundedVector253*crRightHandSideBoundedVector464 - crRightHandSideBoundedVector255*crRightHandSideBoundedVector464 - crRightHandSideBoundedVector268*(-DN_DX_2_0*crRightHandSideBoundedVector642 - DN_DX_2_1*crRightHandSideBoundedVector542 + crRightHandSideBoundedVector749) - crRightHandSideBoundedVector279*(-DN_DX_2_0*crRightHandSideBoundedVector650 - DN_DX_2_1*crRightHandSideBoundedVector649 + crRightHandSideBoundedVector653) - crRightHandSideBoundedVector288*(-DN_DX_2_0*crRightHandSideBoundedVector661 - DN_DX_2_1*crRightHandSideBoundedVector660 + crRightHandSideBoundedVector19*crRightHandSideBoundedVector782 + crRightHandSideBoundedVector410*crRightHandSideBoundedVector758 + crRightHandSideBoundedVector630*crRightHandSideBoundedVector640 - crRightHandSideBoundedVector662*crRightHandSideBoundedVector792 - crRightHandSideBoundedVector75*crRightHandSideBoundedVector766 - crRightHandSideBoundedVector784 + crRightHandSideBoundedVector785 - crRightHandSideBoundedVector793) + crRightHandSideBoundedVector299*crRightHandSideBoundedVector655 - crRightHandSideBoundedVector313*(DN_DX_2_0*crRightHandSideBoundedVector605 + DN_DX_2_1*crRightHandSideBoundedVector736 + crRightHandSideBoundedVector740) - crRightHandSideBoundedVector338*(DN_DX_2_0*crRightHandSideBoundedVector614 + DN_DX_2_1*crRightHandSideBoundedVector742 + crRightHandSideBoundedVector625) - crRightHandSideBoundedVector358*(DN_DX_2_0*crRightHandSideBoundedVector626 + DN_DX_2_1*crRightHandSideBoundedVector743 + crRightHandSideBoundedVector179*crRightHandSideBoundedVector791 - crRightHandSideBoundedVector253*crRightHandSideBoundedVector631 - crRightHandSideBoundedVector255*crRightHandSideBoundedVector631 + crRightHandSideBoundedVector299*crRightHandSideBoundedVector629 - crRightHandSideBoundedVector408*crRightHandSideBoundedVector762 - crRightHandSideBoundedVector410*crRightHandSideBoundedVector763 + crRightHandSideBoundedVector484*crRightHandSideBoundedVector630 + crRightHandSideBoundedVector601 + crRightHandSideBoundedVector737 + crRightHandSideBoundedVector791*crRightHandSideBoundedVector83 + 0.44444444444444442*r_ext[2]) - crRightHandSideBoundedVector367*(-DN_DX_2_0*crRightHandSideBoundedVector541 - DN_DX_2_1*crRightHandSideBoundedVector642 + crRightHandSideBoundedVector747) - crRightHandSideBoundedVector376*(-DN_DX_2_0*crRightHandSideBoundedVector644 - DN_DX_2_1*crRightHandSideBoundedVector650 + crRightHandSideBoundedVector648) - crRightHandSideBoundedVector382*(-DN_DX_2_0*crRightHandSideBoundedVector655 - DN_DX_2_1*crRightHandSideBoundedVector661 + crRightHandSideBoundedVector19*crRightHandSideBoundedVector765 + crRightHandSideBoundedVector408*crRightHandSideBoundedVector758 - crRightHandSideBoundedVector55*crRightHandSideBoundedVector783 + crRightHandSideBoundedVector629*crRightHandSideBoundedVector636 - crRightHandSideBoundedVector65*crRightHandSideBoundedVector766 - crRightHandSideBoundedVector658*crRightHandSideBoundedVector792 - crRightHandSideBoundedVector767 + crRightHandSideBoundedVector768) + 0.66666666666666663*crRightHandSideBoundedVector456 + 0.66666666666666663*crRightHandSideBoundedVector457 + 0.66666666666666663*crRightHandSideBoundedVector458 - 0.66666666666666663*crRightHandSideBoundedVector459 - 0.66666666666666663*crRightHandSideBoundedVector460 + crRightHandSideBoundedVector484*crRightHandSideBoundedVector660 - crRightHandSideBoundedVector665*(crRightHandSideBoundedVector771 + crRightHandSideBoundedVector786) - crRightHandSideBoundedVector666*(crRightHandSideBoundedVector774 + crRightHandSideBoundedVector787) - crRightHandSideBoundedVector667*(crRightHandSideBoundedVector780 + crRightHandSideBoundedVector789) + crRightHandSideBoundedVector669 + crRightHandSideBoundedVector752;

    } else {
        const double crRightHandSideBoundedVector0 =             0.16666666666666666*U_1_0;
const double crRightHandSideBoundedVector1 =             0.16666666666666666*U_2_0;
const double crRightHandSideBoundedVector2 =             0.66666666666666663*U_0_0 + crRightHandSideBoundedVector0 + crRightHandSideBoundedVector1;
const double crRightHandSideBoundedVector3 =             1.0/crRightHandSideBoundedVector2;
const double crRightHandSideBoundedVector4 =             stab_c1/pow(h, 2);
const double crRightHandSideBoundedVector5 =             crRightHandSideBoundedVector3*crRightHandSideBoundedVector4;
const double crRightHandSideBoundedVector6 =             1.3333333333333333*mu;
const double crRightHandSideBoundedVector7 =             0.25*U_1_0;
const double crRightHandSideBoundedVector8 =             0.25*U_2_0;
const double crRightHandSideBoundedVector9 =             U_0_0 + crRightHandSideBoundedVector7 + crRightHandSideBoundedVector8;
const double crRightHandSideBoundedVector10 =             pow(crRightHandSideBoundedVector9, -2);
const double crRightHandSideBoundedVector11 =             0.25*U_1_1;
const double crRightHandSideBoundedVector12 =             0.25*U_2_1;
const double crRightHandSideBoundedVector13 =             pow(U_0_1 + crRightHandSideBoundedVector11 + crRightHandSideBoundedVector12, 2);
const double crRightHandSideBoundedVector14 =             0.25*U_1_2;
const double crRightHandSideBoundedVector15 =             0.25*U_2_2;
const double crRightHandSideBoundedVector16 =             pow(U_0_2 + crRightHandSideBoundedVector14 + crRightHandSideBoundedVector15, 2);
const double crRightHandSideBoundedVector17 =             crRightHandSideBoundedVector13 + crRightHandSideBoundedVector16;
const double crRightHandSideBoundedVector18 =             sqrt(gamma);
const double crRightHandSideBoundedVector19 =             gamma - 1;
const double crRightHandSideBoundedVector20 =             0.5*crRightHandSideBoundedVector10;
const double crRightHandSideBoundedVector21 =             0.16666666666666666*U_1_3;
const double crRightHandSideBoundedVector22 =             0.16666666666666666*U_2_3;
const double crRightHandSideBoundedVector23 =             0.66666666666666663*U_0_3 + crRightHandSideBoundedVector21 + crRightHandSideBoundedVector22;
const double crRightHandSideBoundedVector24 =             crRightHandSideBoundedVector13*crRightHandSideBoundedVector20 + crRightHandSideBoundedVector16*crRightHandSideBoundedVector20 - crRightHandSideBoundedVector23*crRightHandSideBoundedVector3;
const double crRightHandSideBoundedVector25 =             crRightHandSideBoundedVector19*crRightHandSideBoundedVector24;
const double crRightHandSideBoundedVector26 =             stab_c2/h;
const double crRightHandSideBoundedVector27 =             pow(crRightHandSideBoundedVector19, -2);
const double crRightHandSideBoundedVector28 =             0.25*r_ext[1];
const double crRightHandSideBoundedVector29 =             0.25*r_ext[2];
const double crRightHandSideBoundedVector30 =             pow(crRightHandSideBoundedVector28 + crRightHandSideBoundedVector29 + r_ext[0], 2);
const double crRightHandSideBoundedVector31 =             0.25*f_ext(1,0);
const double crRightHandSideBoundedVector32 =             0.25*f_ext(2,0);
const double crRightHandSideBoundedVector33 =             0.25*f_ext(1,1);
const double crRightHandSideBoundedVector34 =             0.25*f_ext(2,1);
const double crRightHandSideBoundedVector35 =             crRightHandSideBoundedVector25*(pow(crRightHandSideBoundedVector31 + crRightHandSideBoundedVector32 + f_ext(0,0), 2) + pow(crRightHandSideBoundedVector33 + crRightHandSideBoundedVector34 + f_ext(0,1), 2));
const double crRightHandSideBoundedVector36 =             0.88888888888888884*gamma;
const double crRightHandSideBoundedVector37 =             0.44444444444444442*gamma;
const double crRightHandSideBoundedVector38 =             1.0/gamma;
const double crRightHandSideBoundedVector39 =             0.70710678118654757*crRightHandSideBoundedVector38*stab_c3;
const double crRightHandSideBoundedVector40 =             crRightHandSideBoundedVector26*(crRightHandSideBoundedVector18*sqrt(-crRightHandSideBoundedVector25) + 1.0*sqrt(crRightHandSideBoundedVector10*crRightHandSideBoundedVector17)) + crRightHandSideBoundedVector39*sqrt(crRightHandSideBoundedVector27*(0.44444444444444442*crRightHandSideBoundedVector30 - crRightHandSideBoundedVector35*crRightHandSideBoundedVector36 + 1.3333333333333333*sqrt(crRightHandSideBoundedVector30*(0.11111111111111109*crRightHandSideBoundedVector30 - crRightHandSideBoundedVector35*crRightHandSideBoundedVector37)))/pow(crRightHandSideBoundedVector24, 2));
const double crRightHandSideBoundedVector41 =             1.0/(crRightHandSideBoundedVector40 + crRightHandSideBoundedVector5*crRightHandSideBoundedVector6);
const double crRightHandSideBoundedVector42 =             0.16666666666666666*dUdt_1_1;
const double crRightHandSideBoundedVector43 =             0.16666666666666666*f_ext(1,0);
const double crRightHandSideBoundedVector44 =             0.16666666666666666*f_ext(2,0);
const double crRightHandSideBoundedVector45 =             crRightHandSideBoundedVector43 + crRightHandSideBoundedVector44 + 0.66666666666666663*f_ext(0,0);
const double crRightHandSideBoundedVector46 =             crRightHandSideBoundedVector2*crRightHandSideBoundedVector45;
const double crRightHandSideBoundedVector47 =             DN_DX_0_0*U_0_0 + DN_DX_1_0*U_1_0 + DN_DX_2_0*U_2_0;
const double crRightHandSideBoundedVector48 =             1.0000000000000002*crRightHandSideBoundedVector13;
const double crRightHandSideBoundedVector49 =             0.50000000000000011*crRightHandSideBoundedVector19;
const double crRightHandSideBoundedVector50 =             crRightHandSideBoundedVector17*crRightHandSideBoundedVector49;
const double crRightHandSideBoundedVector51 =             crRightHandSideBoundedVector10*(-crRightHandSideBoundedVector48 + crRightHandSideBoundedVector50);
const double crRightHandSideBoundedVector52 =             crRightHandSideBoundedVector47*crRightHandSideBoundedVector51;
const double crRightHandSideBoundedVector53 =             DN_DX_0_1*U_0_2;
const double crRightHandSideBoundedVector54 =             DN_DX_1_1*U_1_2;
const double crRightHandSideBoundedVector55 =             DN_DX_2_1*U_2_2;
const double crRightHandSideBoundedVector56 =             crRightHandSideBoundedVector53 + crRightHandSideBoundedVector54 + crRightHandSideBoundedVector55;
const double crRightHandSideBoundedVector57 =             0.66666666666666663*U_0_1;
const double crRightHandSideBoundedVector58 =             0.16666666666666666*U_1_1;
const double crRightHandSideBoundedVector59 =             0.16666666666666666*U_2_1;
const double crRightHandSideBoundedVector60 =             crRightHandSideBoundedVector57 + crRightHandSideBoundedVector58 + crRightHandSideBoundedVector59;
const double crRightHandSideBoundedVector61 =             crRightHandSideBoundedVector3*crRightHandSideBoundedVector60;
const double crRightHandSideBoundedVector62 =             crRightHandSideBoundedVector56*crRightHandSideBoundedVector61;
const double crRightHandSideBoundedVector63 =             DN_DX_0_1*U_0_1;
const double crRightHandSideBoundedVector64 =             DN_DX_1_1*U_1_1;
const double crRightHandSideBoundedVector65 =             DN_DX_2_1*U_2_1;
const double crRightHandSideBoundedVector66 =             crRightHandSideBoundedVector63 + crRightHandSideBoundedVector64 + crRightHandSideBoundedVector65;
const double crRightHandSideBoundedVector67 =             0.66666666666666663*U_0_2;
const double crRightHandSideBoundedVector68 =             0.16666666666666666*U_1_2;
const double crRightHandSideBoundedVector69 =             0.16666666666666666*U_2_2;
const double crRightHandSideBoundedVector70 =             crRightHandSideBoundedVector67 + crRightHandSideBoundedVector68 + crRightHandSideBoundedVector69;
const double crRightHandSideBoundedVector71 =             crRightHandSideBoundedVector3*crRightHandSideBoundedVector70;
const double crRightHandSideBoundedVector72 =             crRightHandSideBoundedVector66*crRightHandSideBoundedVector71;
const double crRightHandSideBoundedVector73 =             DN_DX_0_0*U_0_1;
const double crRightHandSideBoundedVector74 =             DN_DX_1_0*U_1_1;
const double crRightHandSideBoundedVector75 =             DN_DX_2_0*U_2_1;
const double crRightHandSideBoundedVector76 =             crRightHandSideBoundedVector73 + crRightHandSideBoundedVector74 + crRightHandSideBoundedVector75;
const double crRightHandSideBoundedVector77 =             1.0*gamma;
const double crRightHandSideBoundedVector78 =             3.0 - crRightHandSideBoundedVector77;
const double crRightHandSideBoundedVector79 =             crRightHandSideBoundedVector76*crRightHandSideBoundedVector78;
const double crRightHandSideBoundedVector80 =             DN_DX_0_0*U_0_2;
const double crRightHandSideBoundedVector81 =             DN_DX_1_0*U_1_2;
const double crRightHandSideBoundedVector82 =             DN_DX_2_0*U_2_2;
const double crRightHandSideBoundedVector83 =             crRightHandSideBoundedVector80 + crRightHandSideBoundedVector81 + crRightHandSideBoundedVector82;
const double crRightHandSideBoundedVector84 =             crRightHandSideBoundedVector19*crRightHandSideBoundedVector83;
const double crRightHandSideBoundedVector85 =             1.0*crRightHandSideBoundedVector84;
const double crRightHandSideBoundedVector86 =             DN_DX_0_1*U_0_0 + DN_DX_1_1*U_1_0 + DN_DX_2_1*U_2_0;
const double crRightHandSideBoundedVector87 =             crRightHandSideBoundedVector70*crRightHandSideBoundedVector86;
const double crRightHandSideBoundedVector88 =             2.2500000000000004*crRightHandSideBoundedVector60;
const double crRightHandSideBoundedVector89 =             crRightHandSideBoundedVector10*crRightHandSideBoundedVector88;
const double crRightHandSideBoundedVector90 =             DN_DX_0_0*U_0_3 + DN_DX_1_0*U_1_3 + DN_DX_2_0*U_2_3;
const double crRightHandSideBoundedVector91 =             crRightHandSideBoundedVector19*crRightHandSideBoundedVector90;
const double crRightHandSideBoundedVector92 =             crRightHandSideBoundedVector91 + 0.16666666666666666*dUdt_2_1;
const double crRightHandSideBoundedVector93 =             crRightHandSideBoundedVector41*(crRightHandSideBoundedVector42 - crRightHandSideBoundedVector46 + crRightHandSideBoundedVector52 + crRightHandSideBoundedVector61*crRightHandSideBoundedVector79 + crRightHandSideBoundedVector62 - crRightHandSideBoundedVector71*crRightHandSideBoundedVector85 + crRightHandSideBoundedVector72 - crRightHandSideBoundedVector87*crRightHandSideBoundedVector89 + crRightHandSideBoundedVector92 + 0.66666666666666663*dUdt_0_1);
const double crRightHandSideBoundedVector94 =             0.16666666666666666*U_0_0;
const double crRightHandSideBoundedVector95 =             0.66666666666666663*U_1_0 + crRightHandSideBoundedVector1 + crRightHandSideBoundedVector94;
const double crRightHandSideBoundedVector96 =             1.0/crRightHandSideBoundedVector95;
const double crRightHandSideBoundedVector97 =             crRightHandSideBoundedVector4*crRightHandSideBoundedVector6;
const double crRightHandSideBoundedVector98 =             0.25*U_0_0;
const double crRightHandSideBoundedVector99 =             U_1_0 + crRightHandSideBoundedVector8 + crRightHandSideBoundedVector98;
const double crRightHandSideBoundedVector100 =             pow(crRightHandSideBoundedVector99, -2);
const double crRightHandSideBoundedVector101 =             0.25*U_0_1;
const double crRightHandSideBoundedVector102 =             pow(U_1_1 + crRightHandSideBoundedVector101 + crRightHandSideBoundedVector12, 2);
const double crRightHandSideBoundedVector103 =             0.25*U_0_2;
const double crRightHandSideBoundedVector104 =             pow(U_1_2 + crRightHandSideBoundedVector103 + crRightHandSideBoundedVector15, 2);
const double crRightHandSideBoundedVector105 =             crRightHandSideBoundedVector102 + crRightHandSideBoundedVector104;
const double crRightHandSideBoundedVector106 =             0.5*crRightHandSideBoundedVector100;
const double crRightHandSideBoundedVector107 =             0.16666666666666666*U_0_3;
const double crRightHandSideBoundedVector108 =             0.66666666666666663*U_1_3 + crRightHandSideBoundedVector107 + crRightHandSideBoundedVector22;
const double crRightHandSideBoundedVector109 =             crRightHandSideBoundedVector102*crRightHandSideBoundedVector106 + crRightHandSideBoundedVector104*crRightHandSideBoundedVector106 - crRightHandSideBoundedVector108*crRightHandSideBoundedVector96;
const double crRightHandSideBoundedVector110 =             crRightHandSideBoundedVector109*crRightHandSideBoundedVector19;
const double crRightHandSideBoundedVector111 =             0.25*r_ext[0];
const double crRightHandSideBoundedVector112 =             pow(crRightHandSideBoundedVector111 + crRightHandSideBoundedVector29 + r_ext[1], 2);
const double crRightHandSideBoundedVector113 =             0.25*f_ext(0,0);
const double crRightHandSideBoundedVector114 =             0.25*f_ext(0,1);
const double crRightHandSideBoundedVector115 =             crRightHandSideBoundedVector110*(pow(crRightHandSideBoundedVector113 + crRightHandSideBoundedVector32 + f_ext(1,0), 2) + pow(crRightHandSideBoundedVector114 + crRightHandSideBoundedVector34 + f_ext(1,1), 2));
const double crRightHandSideBoundedVector116 =             crRightHandSideBoundedVector26*(crRightHandSideBoundedVector18*sqrt(-crRightHandSideBoundedVector110) + 1.0*sqrt(crRightHandSideBoundedVector100*crRightHandSideBoundedVector105)) + crRightHandSideBoundedVector39*sqrt(crRightHandSideBoundedVector27*(0.44444444444444442*crRightHandSideBoundedVector112 - crRightHandSideBoundedVector115*crRightHandSideBoundedVector36 + 1.3333333333333333*sqrt(crRightHandSideBoundedVector112*(0.11111111111111109*crRightHandSideBoundedVector112 - crRightHandSideBoundedVector115*crRightHandSideBoundedVector37)))/pow(crRightHandSideBoundedVector109, 2));
const double crRightHandSideBoundedVector117 =             1.0/(crRightHandSideBoundedVector116 + crRightHandSideBoundedVector96*crRightHandSideBoundedVector97);
const double crRightHandSideBoundedVector118 =             0.16666666666666666*dUdt_0_1;
const double crRightHandSideBoundedVector119 =             0.16666666666666666*f_ext(0,0);
const double crRightHandSideBoundedVector120 =             crRightHandSideBoundedVector119 + crRightHandSideBoundedVector44 + 0.66666666666666663*f_ext(1,0);
const double crRightHandSideBoundedVector121 =             crRightHandSideBoundedVector120*crRightHandSideBoundedVector95;
const double crRightHandSideBoundedVector122 =             1.0000000000000002*crRightHandSideBoundedVector102;
const double crRightHandSideBoundedVector123 =             crRightHandSideBoundedVector105*crRightHandSideBoundedVector49;
const double crRightHandSideBoundedVector124 =             crRightHandSideBoundedVector100*(-crRightHandSideBoundedVector122 + crRightHandSideBoundedVector123);
const double crRightHandSideBoundedVector125 =             crRightHandSideBoundedVector124*crRightHandSideBoundedVector47;
const double crRightHandSideBoundedVector126 =             0.16666666666666666*U_0_1;
const double crRightHandSideBoundedVector127 =             0.66666666666666663*U_1_1 + crRightHandSideBoundedVector126 + crRightHandSideBoundedVector59;
const double crRightHandSideBoundedVector128 =             crRightHandSideBoundedVector127*crRightHandSideBoundedVector96;
const double crRightHandSideBoundedVector129 =             crRightHandSideBoundedVector128*crRightHandSideBoundedVector56;
const double crRightHandSideBoundedVector130 =             0.16666666666666666*U_0_2;
const double crRightHandSideBoundedVector131 =             0.66666666666666663*U_1_2 + crRightHandSideBoundedVector130 + crRightHandSideBoundedVector69;
const double crRightHandSideBoundedVector132 =             crRightHandSideBoundedVector131*crRightHandSideBoundedVector96;
const double crRightHandSideBoundedVector133 =             crRightHandSideBoundedVector132*crRightHandSideBoundedVector66;
const double crRightHandSideBoundedVector134 =             crRightHandSideBoundedVector131*crRightHandSideBoundedVector86;
const double crRightHandSideBoundedVector135 =             2.2500000000000004*crRightHandSideBoundedVector127;
const double crRightHandSideBoundedVector136 =             crRightHandSideBoundedVector100*crRightHandSideBoundedVector135;
const double crRightHandSideBoundedVector137 =             crRightHandSideBoundedVector117*(crRightHandSideBoundedVector118 - crRightHandSideBoundedVector121 + crRightHandSideBoundedVector125 + crRightHandSideBoundedVector128*crRightHandSideBoundedVector79 + crRightHandSideBoundedVector129 - crRightHandSideBoundedVector132*crRightHandSideBoundedVector85 + crRightHandSideBoundedVector133 - crRightHandSideBoundedVector134*crRightHandSideBoundedVector136 + crRightHandSideBoundedVector92 + 0.66666666666666663*dUdt_1_1);
const double crRightHandSideBoundedVector138 =             0.66666666666666663*U_2_0 + crRightHandSideBoundedVector0 + crRightHandSideBoundedVector94;
const double crRightHandSideBoundedVector139 =             1.0/crRightHandSideBoundedVector138;
const double crRightHandSideBoundedVector140 =             U_2_0 + crRightHandSideBoundedVector7 + crRightHandSideBoundedVector98;
const double crRightHandSideBoundedVector141 =             pow(crRightHandSideBoundedVector140, -2);
const double crRightHandSideBoundedVector142 =             pow(U_2_1 + crRightHandSideBoundedVector101 + crRightHandSideBoundedVector11, 2);
const double crRightHandSideBoundedVector143 =             pow(U_2_2 + crRightHandSideBoundedVector103 + crRightHandSideBoundedVector14, 2);
const double crRightHandSideBoundedVector144 =             crRightHandSideBoundedVector142 + crRightHandSideBoundedVector143;
const double crRightHandSideBoundedVector145 =             0.5*crRightHandSideBoundedVector141;
const double crRightHandSideBoundedVector146 =             0.66666666666666663*U_2_3 + crRightHandSideBoundedVector107 + crRightHandSideBoundedVector21;
const double crRightHandSideBoundedVector147 =             -crRightHandSideBoundedVector139*crRightHandSideBoundedVector146 + crRightHandSideBoundedVector142*crRightHandSideBoundedVector145 + crRightHandSideBoundedVector143*crRightHandSideBoundedVector145;
const double crRightHandSideBoundedVector148 =             crRightHandSideBoundedVector147*crRightHandSideBoundedVector19;
const double crRightHandSideBoundedVector149 =             pow(crRightHandSideBoundedVector111 + crRightHandSideBoundedVector28 + r_ext[2], 2);
const double crRightHandSideBoundedVector150 =             crRightHandSideBoundedVector148*(pow(crRightHandSideBoundedVector113 + crRightHandSideBoundedVector31 + f_ext(2,0), 2) + pow(crRightHandSideBoundedVector114 + crRightHandSideBoundedVector33 + f_ext(2,1), 2));
const double crRightHandSideBoundedVector151 =             crRightHandSideBoundedVector26*(crRightHandSideBoundedVector18*sqrt(-crRightHandSideBoundedVector148) + 1.0*sqrt(crRightHandSideBoundedVector141*crRightHandSideBoundedVector144)) + crRightHandSideBoundedVector39*sqrt(crRightHandSideBoundedVector27*(0.44444444444444442*crRightHandSideBoundedVector149 - crRightHandSideBoundedVector150*crRightHandSideBoundedVector36 + 1.3333333333333333*sqrt(crRightHandSideBoundedVector149*(0.11111111111111109*crRightHandSideBoundedVector149 - crRightHandSideBoundedVector150*crRightHandSideBoundedVector37)))/pow(crRightHandSideBoundedVector147, 2));
const double crRightHandSideBoundedVector152 =             1.0/(crRightHandSideBoundedVector139*crRightHandSideBoundedVector97 + crRightHandSideBoundedVector151);
const double crRightHandSideBoundedVector153 =             crRightHandSideBoundedVector119 + crRightHandSideBoundedVector43 + 0.66666666666666663*f_ext(2,0);
const double crRightHandSideBoundedVector154 =             crRightHandSideBoundedVector138*crRightHandSideBoundedVector153;
const double crRightHandSideBoundedVector155 =             0.66666666666666663*U_2_1 + crRightHandSideBoundedVector126 + crRightHandSideBoundedVector58;
const double crRightHandSideBoundedVector156 =             crRightHandSideBoundedVector139*crRightHandSideBoundedVector155;
const double crRightHandSideBoundedVector157 =             crRightHandSideBoundedVector156*crRightHandSideBoundedVector56;
const double crRightHandSideBoundedVector158 =             0.66666666666666663*U_2_2 + crRightHandSideBoundedVector130 + crRightHandSideBoundedVector68;
const double crRightHandSideBoundedVector159 =             crRightHandSideBoundedVector139*crRightHandSideBoundedVector158;
const double crRightHandSideBoundedVector160 =             crRightHandSideBoundedVector159*crRightHandSideBoundedVector66;
const double crRightHandSideBoundedVector161 =             crRightHandSideBoundedVector158*crRightHandSideBoundedVector86;
const double crRightHandSideBoundedVector162 =             2.2500000000000004*crRightHandSideBoundedVector155;
const double crRightHandSideBoundedVector163 =             crRightHandSideBoundedVector141*crRightHandSideBoundedVector162;
const double crRightHandSideBoundedVector164 =             1.0000000000000002*crRightHandSideBoundedVector142;
const double crRightHandSideBoundedVector165 =             crRightHandSideBoundedVector144*crRightHandSideBoundedVector49;
const double crRightHandSideBoundedVector166 =             crRightHandSideBoundedVector141*(-crRightHandSideBoundedVector164 + crRightHandSideBoundedVector165);
const double crRightHandSideBoundedVector167 =             crRightHandSideBoundedVector166*crRightHandSideBoundedVector47;
const double crRightHandSideBoundedVector168 =             crRightHandSideBoundedVector152*(crRightHandSideBoundedVector118 - crRightHandSideBoundedVector154 + crRightHandSideBoundedVector156*crRightHandSideBoundedVector79 + crRightHandSideBoundedVector157 - crRightHandSideBoundedVector159*crRightHandSideBoundedVector85 + crRightHandSideBoundedVector160 - crRightHandSideBoundedVector161*crRightHandSideBoundedVector163 + crRightHandSideBoundedVector167 + crRightHandSideBoundedVector42 + crRightHandSideBoundedVector91 + 0.66666666666666663*dUdt_2_1);
const double crRightHandSideBoundedVector169 =             0.16666666666666666*dUdt_1_2;
const double crRightHandSideBoundedVector170 =             0.16666666666666666*f_ext(1,1);
const double crRightHandSideBoundedVector171 =             0.16666666666666666*f_ext(2,1);
const double crRightHandSideBoundedVector172 =             crRightHandSideBoundedVector170 + crRightHandSideBoundedVector171 + 0.66666666666666663*f_ext(0,1);
const double crRightHandSideBoundedVector173 =             crRightHandSideBoundedVector172*crRightHandSideBoundedVector2;
const double crRightHandSideBoundedVector174 =             1.0000000000000002*crRightHandSideBoundedVector16;
const double crRightHandSideBoundedVector175 =             crRightHandSideBoundedVector10*(-crRightHandSideBoundedVector174 + crRightHandSideBoundedVector50);
const double crRightHandSideBoundedVector176 =             crRightHandSideBoundedVector175*crRightHandSideBoundedVector86;
const double crRightHandSideBoundedVector177 =             crRightHandSideBoundedVector61*crRightHandSideBoundedVector83;
const double crRightHandSideBoundedVector178 =             crRightHandSideBoundedVector71*crRightHandSideBoundedVector76;
const double crRightHandSideBoundedVector179 =             crRightHandSideBoundedVector56*crRightHandSideBoundedVector78;
const double crRightHandSideBoundedVector180 =             crRightHandSideBoundedVector19*crRightHandSideBoundedVector66;
const double crRightHandSideBoundedVector181 =             1.0*crRightHandSideBoundedVector180;
const double crRightHandSideBoundedVector182 =             2.2500000000000004*crRightHandSideBoundedVector70;
const double crRightHandSideBoundedVector183 =             crRightHandSideBoundedVector47*crRightHandSideBoundedVector60;
const double crRightHandSideBoundedVector184 =             crRightHandSideBoundedVector10*crRightHandSideBoundedVector183;
const double crRightHandSideBoundedVector185 =             DN_DX_0_1*U_0_3 + DN_DX_1_1*U_1_3 + DN_DX_2_1*U_2_3;
const double crRightHandSideBoundedVector186 =             crRightHandSideBoundedVector185*crRightHandSideBoundedVector19;
const double crRightHandSideBoundedVector187 =             crRightHandSideBoundedVector186 + 0.16666666666666666*dUdt_2_2;
const double crRightHandSideBoundedVector188 =             crRightHandSideBoundedVector41*(crRightHandSideBoundedVector169 - crRightHandSideBoundedVector173 + crRightHandSideBoundedVector176 + crRightHandSideBoundedVector177 + crRightHandSideBoundedVector178 + crRightHandSideBoundedVector179*crRightHandSideBoundedVector71 - crRightHandSideBoundedVector181*crRightHandSideBoundedVector61 - crRightHandSideBoundedVector182*crRightHandSideBoundedVector184 + crRightHandSideBoundedVector187 + 0.66666666666666663*dUdt_0_2);
const double crRightHandSideBoundedVector189 =             0.16666666666666666*dUdt_0_2;
const double crRightHandSideBoundedVector190 =             0.16666666666666666*f_ext(0,1);
const double crRightHandSideBoundedVector191 =             crRightHandSideBoundedVector171 + crRightHandSideBoundedVector190 + 0.66666666666666663*f_ext(1,1);
const double crRightHandSideBoundedVector192 =             crRightHandSideBoundedVector191*crRightHandSideBoundedVector95;
const double crRightHandSideBoundedVector193 =             1.0000000000000002*crRightHandSideBoundedVector104;
const double crRightHandSideBoundedVector194 =             crRightHandSideBoundedVector100*(crRightHandSideBoundedVector123 - crRightHandSideBoundedVector193);
const double crRightHandSideBoundedVector195 =             crRightHandSideBoundedVector194*crRightHandSideBoundedVector86;
const double crRightHandSideBoundedVector196 =             crRightHandSideBoundedVector128*crRightHandSideBoundedVector83;
const double crRightHandSideBoundedVector197 =             crRightHandSideBoundedVector132*crRightHandSideBoundedVector76;
const double crRightHandSideBoundedVector198 =             2.2500000000000004*crRightHandSideBoundedVector131;
const double crRightHandSideBoundedVector199 =             crRightHandSideBoundedVector127*crRightHandSideBoundedVector47;
const double crRightHandSideBoundedVector200 =             crRightHandSideBoundedVector100*crRightHandSideBoundedVector199;
const double crRightHandSideBoundedVector201 =             crRightHandSideBoundedVector117*(-crRightHandSideBoundedVector128*crRightHandSideBoundedVector181 + crRightHandSideBoundedVector132*crRightHandSideBoundedVector179 + crRightHandSideBoundedVector187 + crRightHandSideBoundedVector189 - crRightHandSideBoundedVector192 + crRightHandSideBoundedVector195 + crRightHandSideBoundedVector196 + crRightHandSideBoundedVector197 - crRightHandSideBoundedVector198*crRightHandSideBoundedVector200 + 0.66666666666666663*dUdt_1_2);
const double crRightHandSideBoundedVector202 =             crRightHandSideBoundedVector170 + crRightHandSideBoundedVector190 + 0.66666666666666663*f_ext(2,1);
const double crRightHandSideBoundedVector203 =             crRightHandSideBoundedVector138*crRightHandSideBoundedVector202;
const double crRightHandSideBoundedVector204 =             crRightHandSideBoundedVector156*crRightHandSideBoundedVector83;
const double crRightHandSideBoundedVector205 =             crRightHandSideBoundedVector159*crRightHandSideBoundedVector76;
const double crRightHandSideBoundedVector206 =             2.2500000000000004*crRightHandSideBoundedVector158;
const double crRightHandSideBoundedVector207 =             crRightHandSideBoundedVector155*crRightHandSideBoundedVector47;
const double crRightHandSideBoundedVector208 =             crRightHandSideBoundedVector141*crRightHandSideBoundedVector207;
const double crRightHandSideBoundedVector209 =             1.0000000000000002*crRightHandSideBoundedVector143;
const double crRightHandSideBoundedVector210 =             crRightHandSideBoundedVector141*(crRightHandSideBoundedVector165 - crRightHandSideBoundedVector209);
const double crRightHandSideBoundedVector211 =             crRightHandSideBoundedVector210*crRightHandSideBoundedVector86;
const double crRightHandSideBoundedVector212 =             crRightHandSideBoundedVector152*(-crRightHandSideBoundedVector156*crRightHandSideBoundedVector181 + crRightHandSideBoundedVector159*crRightHandSideBoundedVector179 + crRightHandSideBoundedVector169 + crRightHandSideBoundedVector186 + crRightHandSideBoundedVector189 - crRightHandSideBoundedVector203 + crRightHandSideBoundedVector204 + crRightHandSideBoundedVector205 - crRightHandSideBoundedVector206*crRightHandSideBoundedVector208 + crRightHandSideBoundedVector211 + 0.66666666666666663*dUdt_2_2);
const double crRightHandSideBoundedVector213 =             crRightHandSideBoundedVector56 + crRightHandSideBoundedVector76;
const double crRightHandSideBoundedVector214 =             crRightHandSideBoundedVector2*crRightHandSideBoundedVector76;
const double crRightHandSideBoundedVector215 =             -crRightHandSideBoundedVector183 + crRightHandSideBoundedVector214;
const double crRightHandSideBoundedVector216 =             crRightHandSideBoundedVector10*mu;
const double crRightHandSideBoundedVector217 =             4.5000000000000009*crRightHandSideBoundedVector216;
const double crRightHandSideBoundedVector218 =             crRightHandSideBoundedVector2*crRightHandSideBoundedVector56;
const double crRightHandSideBoundedVector219 =             crRightHandSideBoundedVector218 - crRightHandSideBoundedVector87;
const double crRightHandSideBoundedVector220 =             -1.5000000000000002*crRightHandSideBoundedVector216*(crRightHandSideBoundedVector215 + crRightHandSideBoundedVector219);
const double crRightHandSideBoundedVector221 =             0.66666666666666663*crRightHandSideBoundedVector183;
const double crRightHandSideBoundedVector222 =             crRightHandSideBoundedVector3*(-0.66666666666666663*crRightHandSideBoundedVector214 + 1.3333333333333335*crRightHandSideBoundedVector218 + crRightHandSideBoundedVector221 - 1.3333333333333335*crRightHandSideBoundedVector87);
const double crRightHandSideBoundedVector223 =             1.0/lin_m_norm;
const double crRightHandSideBoundedVector224 =             crRightHandSideBoundedVector223*lin_m[0]*lin_m[1];
const double crRightHandSideBoundedVector225 =             crRightHandSideBoundedVector224*nu_st;
const double crRightHandSideBoundedVector226 =             0.66666666666666663*crRightHandSideBoundedVector87;
const double crRightHandSideBoundedVector227 =             crRightHandSideBoundedVector3*(-1.3333333333333335*crRightHandSideBoundedVector183 + 1.3333333333333335*crRightHandSideBoundedVector214 - 0.66666666666666663*crRightHandSideBoundedVector218 + crRightHandSideBoundedVector226);
const double crRightHandSideBoundedVector228 =             crRightHandSideBoundedVector223*pow(lin_m[0], 2);
const double crRightHandSideBoundedVector229 =             crRightHandSideBoundedVector228*nu_st;
const double crRightHandSideBoundedVector230 =             crRightHandSideBoundedVector224*nu_sc;
const double crRightHandSideBoundedVector231 =             1 - crRightHandSideBoundedVector228;
const double crRightHandSideBoundedVector232 =             crRightHandSideBoundedVector231*nu_sc;
const double crRightHandSideBoundedVector233 =             crRightHandSideBoundedVector215*crRightHandSideBoundedVector217 + crRightHandSideBoundedVector220 + crRightHandSideBoundedVector222*crRightHandSideBoundedVector225 - crRightHandSideBoundedVector222*crRightHandSideBoundedVector230 + crRightHandSideBoundedVector227*crRightHandSideBoundedVector229 + crRightHandSideBoundedVector227*crRightHandSideBoundedVector232;
const double crRightHandSideBoundedVector234 =             crRightHandSideBoundedVector76*crRightHandSideBoundedVector95;
const double crRightHandSideBoundedVector235 =             -crRightHandSideBoundedVector199 + crRightHandSideBoundedVector234;
const double crRightHandSideBoundedVector236 =             crRightHandSideBoundedVector100*mu;
const double crRightHandSideBoundedVector237 =             4.5000000000000009*crRightHandSideBoundedVector236;
const double crRightHandSideBoundedVector238 =             crRightHandSideBoundedVector56*crRightHandSideBoundedVector95;
const double crRightHandSideBoundedVector239 =             -crRightHandSideBoundedVector134 + crRightHandSideBoundedVector238;
const double crRightHandSideBoundedVector240 =             -1.5000000000000002*crRightHandSideBoundedVector236*(crRightHandSideBoundedVector235 + crRightHandSideBoundedVector239);
const double crRightHandSideBoundedVector241 =             0.66666666666666663*crRightHandSideBoundedVector199;
const double crRightHandSideBoundedVector242 =             crRightHandSideBoundedVector96*(-1.3333333333333335*crRightHandSideBoundedVector134 - 0.66666666666666663*crRightHandSideBoundedVector234 + 1.3333333333333335*crRightHandSideBoundedVector238 + crRightHandSideBoundedVector241);
const double crRightHandSideBoundedVector243 =             0.66666666666666663*crRightHandSideBoundedVector134;
const double crRightHandSideBoundedVector244 =             crRightHandSideBoundedVector96*(-1.3333333333333335*crRightHandSideBoundedVector199 + 1.3333333333333335*crRightHandSideBoundedVector234 - 0.66666666666666663*crRightHandSideBoundedVector238 + crRightHandSideBoundedVector243);
const double crRightHandSideBoundedVector245 =             crRightHandSideBoundedVector225*crRightHandSideBoundedVector242 + crRightHandSideBoundedVector229*crRightHandSideBoundedVector244 - crRightHandSideBoundedVector230*crRightHandSideBoundedVector242 + crRightHandSideBoundedVector232*crRightHandSideBoundedVector244 + crRightHandSideBoundedVector235*crRightHandSideBoundedVector237 + crRightHandSideBoundedVector240;
const double crRightHandSideBoundedVector246 =             crRightHandSideBoundedVector138*crRightHandSideBoundedVector76;
const double crRightHandSideBoundedVector247 =             -crRightHandSideBoundedVector207 + crRightHandSideBoundedVector246;
const double crRightHandSideBoundedVector248 =             crRightHandSideBoundedVector141*mu;
const double crRightHandSideBoundedVector249 =             4.5000000000000009*crRightHandSideBoundedVector248;
const double crRightHandSideBoundedVector250 =             crRightHandSideBoundedVector138*crRightHandSideBoundedVector56;
const double crRightHandSideBoundedVector251 =             -crRightHandSideBoundedVector161 + crRightHandSideBoundedVector250;
const double crRightHandSideBoundedVector252 =             -1.5000000000000002*crRightHandSideBoundedVector248*(crRightHandSideBoundedVector247 + crRightHandSideBoundedVector251);
const double crRightHandSideBoundedVector253 =             0.66666666666666663*crRightHandSideBoundedVector207;
const double crRightHandSideBoundedVector254 =             crRightHandSideBoundedVector139*(-1.3333333333333335*crRightHandSideBoundedVector161 - 0.66666666666666663*crRightHandSideBoundedVector246 + 1.3333333333333335*crRightHandSideBoundedVector250 + crRightHandSideBoundedVector253);
const double crRightHandSideBoundedVector255 =             0.66666666666666663*crRightHandSideBoundedVector161;
const double crRightHandSideBoundedVector256 =             crRightHandSideBoundedVector139*(-1.3333333333333335*crRightHandSideBoundedVector207 + 1.3333333333333335*crRightHandSideBoundedVector246 - 0.66666666666666663*crRightHandSideBoundedVector250 + crRightHandSideBoundedVector255);
const double crRightHandSideBoundedVector257 =             crRightHandSideBoundedVector225*crRightHandSideBoundedVector254 + crRightHandSideBoundedVector229*crRightHandSideBoundedVector256 - crRightHandSideBoundedVector230*crRightHandSideBoundedVector254 + crRightHandSideBoundedVector232*crRightHandSideBoundedVector256 + crRightHandSideBoundedVector247*crRightHandSideBoundedVector249 + crRightHandSideBoundedVector252;
const double crRightHandSideBoundedVector258 =             DN_DX_0_1*crRightHandSideBoundedVector61;
const double crRightHandSideBoundedVector259 =             DN_DX_0_1*crRightHandSideBoundedVector57 + 0.66666666666666663*crRightHandSideBoundedVector64 + 0.66666666666666663*crRightHandSideBoundedVector65;
const double crRightHandSideBoundedVector260 =             0.66666666666666663*crRightHandSideBoundedVector3;
const double crRightHandSideBoundedVector261 =             DN_DX_0_0*crRightHandSideBoundedVector71;
const double crRightHandSideBoundedVector262 =             1.0*crRightHandSideBoundedVector19;
const double crRightHandSideBoundedVector263 =             crRightHandSideBoundedVector60*crRightHandSideBoundedVector86;
const double crRightHandSideBoundedVector264 =             1.5000000000000002*crRightHandSideBoundedVector10;
const double crRightHandSideBoundedVector265 =             crRightHandSideBoundedVector47*crRightHandSideBoundedVector70;
const double crRightHandSideBoundedVector266 =             1.5000000000000002*crRightHandSideBoundedVector19;
const double crRightHandSideBoundedVector267 =             crRightHandSideBoundedVector10*crRightHandSideBoundedVector266;
const double crRightHandSideBoundedVector268 =             1.0*crRightHandSideBoundedVector188;
const double crRightHandSideBoundedVector269 =             DN_DX_0_1*crRightHandSideBoundedVector128;
const double crRightHandSideBoundedVector270 =             DN_DX_0_0*crRightHandSideBoundedVector132;
const double crRightHandSideBoundedVector271 =             DN_DX_1_1*crRightHandSideBoundedVector58 + DN_DX_2_1*crRightHandSideBoundedVector59 + 0.16666666666666666*crRightHandSideBoundedVector63;
const double crRightHandSideBoundedVector272 =             0.16666666666666666*crRightHandSideBoundedVector84;
const double crRightHandSideBoundedVector273 =             crRightHandSideBoundedVector127*crRightHandSideBoundedVector86;
const double crRightHandSideBoundedVector274 =             0.37500000000000006*crRightHandSideBoundedVector100;
const double crRightHandSideBoundedVector275 =             crRightHandSideBoundedVector273*crRightHandSideBoundedVector274;
const double crRightHandSideBoundedVector276 =             crRightHandSideBoundedVector131*crRightHandSideBoundedVector47;
const double crRightHandSideBoundedVector277 =             crRightHandSideBoundedVector274*crRightHandSideBoundedVector276;
const double crRightHandSideBoundedVector278 =             -crRightHandSideBoundedVector19*crRightHandSideBoundedVector277 - crRightHandSideBoundedVector271*crRightHandSideBoundedVector96 + crRightHandSideBoundedVector272*crRightHandSideBoundedVector96 + crRightHandSideBoundedVector275;
const double crRightHandSideBoundedVector279 =             1.0*crRightHandSideBoundedVector201;
const double crRightHandSideBoundedVector280 =             DN_DX_0_1*crRightHandSideBoundedVector156;
const double crRightHandSideBoundedVector281 =             DN_DX_0_0*crRightHandSideBoundedVector159;
const double crRightHandSideBoundedVector282 =             crRightHandSideBoundedVector155*crRightHandSideBoundedVector86;
const double crRightHandSideBoundedVector283 =             0.37500000000000006*crRightHandSideBoundedVector141;
const double crRightHandSideBoundedVector284 =             crRightHandSideBoundedVector282*crRightHandSideBoundedVector283;
const double crRightHandSideBoundedVector285 =             crRightHandSideBoundedVector158*crRightHandSideBoundedVector47;
const double crRightHandSideBoundedVector286 =             crRightHandSideBoundedVector283*crRightHandSideBoundedVector285;
const double crRightHandSideBoundedVector287 =             -crRightHandSideBoundedVector139*crRightHandSideBoundedVector271 + crRightHandSideBoundedVector139*crRightHandSideBoundedVector272 - crRightHandSideBoundedVector19*crRightHandSideBoundedVector286 + crRightHandSideBoundedVector284;
const double crRightHandSideBoundedVector288 =             1.0*crRightHandSideBoundedVector212;
const double crRightHandSideBoundedVector289 =             pow(crRightHandSideBoundedVector9, -3);
const double crRightHandSideBoundedVector290 =             0.66666666666666663*crRightHandSideBoundedVector289;
const double crRightHandSideBoundedVector291 =             3.0000000000000004*crRightHandSideBoundedVector13;
const double crRightHandSideBoundedVector292 =             crRightHandSideBoundedVector17*crRightHandSideBoundedVector266;
const double crRightHandSideBoundedVector293 =             crRightHandSideBoundedVector47*(-crRightHandSideBoundedVector291 + crRightHandSideBoundedVector292);
const double crRightHandSideBoundedVector294 =             crRightHandSideBoundedVector264*crRightHandSideBoundedVector60;
const double crRightHandSideBoundedVector295 =             crRightHandSideBoundedVector264*crRightHandSideBoundedVector70;
const double crRightHandSideBoundedVector296 =             4.5*crRightHandSideBoundedVector289;
const double crRightHandSideBoundedVector297 =             crRightHandSideBoundedVector60*crRightHandSideBoundedVector87;
const double crRightHandSideBoundedVector298 =             crRightHandSideBoundedVector296*crRightHandSideBoundedVector297;
const double crRightHandSideBoundedVector299 =             0.66666666666666663*crRightHandSideBoundedVector76;
const double crRightHandSideBoundedVector300 =             2.2500000000000004*gamma - 6.7500000000000018;
const double crRightHandSideBoundedVector301 =             crRightHandSideBoundedVector10*crRightHandSideBoundedVector300;
const double crRightHandSideBoundedVector302 =             crRightHandSideBoundedVector301*crRightHandSideBoundedVector60;
const double crRightHandSideBoundedVector303 =             crRightHandSideBoundedVector267*crRightHandSideBoundedVector70;
const double crRightHandSideBoundedVector304 =             crRightHandSideBoundedVector303*crRightHandSideBoundedVector83;
const double crRightHandSideBoundedVector305 =             DN_DX_0_1*crRightHandSideBoundedVector10;
const double crRightHandSideBoundedVector306 =             crRightHandSideBoundedVector305*crRightHandSideBoundedVector70;
const double crRightHandSideBoundedVector307 =             crRightHandSideBoundedVector306*crRightHandSideBoundedVector88;
const double crRightHandSideBoundedVector308 =             0.1111111111111111*f_ext(1,0);
const double crRightHandSideBoundedVector309 =             0.1111111111111111*f_ext(2,0);
const double crRightHandSideBoundedVector310 =             crRightHandSideBoundedVector308 + crRightHandSideBoundedVector309 + 0.44444444444444442*f_ext(0,0);
const double crRightHandSideBoundedVector311 =             0.16666666666666666*dUdt_1_0;
const double crRightHandSideBoundedVector312 =             crRightHandSideBoundedVector213 + 0.16666666666666666*dUdt_2_0;
const double crRightHandSideBoundedVector313 =             1.0*(crRightHandSideBoundedVector311 + crRightHandSideBoundedVector312 + 0.66666666666666663*dUdt_0_0)/crRightHandSideBoundedVector40;
const double crRightHandSideBoundedVector314 =             DN_DX_0_1*crRightHandSideBoundedVector100;
const double crRightHandSideBoundedVector315 =             crRightHandSideBoundedVector131*crRightHandSideBoundedVector314;
const double crRightHandSideBoundedVector316 =             crRightHandSideBoundedVector135*crRightHandSideBoundedVector315;
const double crRightHandSideBoundedVector317 =             0.16666666666666666*crRightHandSideBoundedVector47;
const double crRightHandSideBoundedVector318 =             pow(crRightHandSideBoundedVector99, -3);
const double crRightHandSideBoundedVector319 =             3.0000000000000004*crRightHandSideBoundedVector102;
const double crRightHandSideBoundedVector320 =             crRightHandSideBoundedVector105*crRightHandSideBoundedVector266;
const double crRightHandSideBoundedVector321 =             crRightHandSideBoundedVector318*(-crRightHandSideBoundedVector319 + crRightHandSideBoundedVector320);
const double crRightHandSideBoundedVector322 =             crRightHandSideBoundedVector127*crRightHandSideBoundedVector274;
const double crRightHandSideBoundedVector323 =             crRightHandSideBoundedVector322*crRightHandSideBoundedVector56;
const double crRightHandSideBoundedVector324 =             crRightHandSideBoundedVector131*crRightHandSideBoundedVector274;
const double crRightHandSideBoundedVector325 =             crRightHandSideBoundedVector324*crRightHandSideBoundedVector66;
const double crRightHandSideBoundedVector326 =             1.125*crRightHandSideBoundedVector318;
const double crRightHandSideBoundedVector327 =             crRightHandSideBoundedVector127*crRightHandSideBoundedVector134;
const double crRightHandSideBoundedVector328 =             crRightHandSideBoundedVector326*crRightHandSideBoundedVector327;
const double crRightHandSideBoundedVector329 =             0.16666666666666666*crRightHandSideBoundedVector76;
const double crRightHandSideBoundedVector330 =             crRightHandSideBoundedVector300*crRightHandSideBoundedVector329;
const double crRightHandSideBoundedVector331 =             crRightHandSideBoundedVector100*crRightHandSideBoundedVector127;
const double crRightHandSideBoundedVector332 =             crRightHandSideBoundedVector324*crRightHandSideBoundedVector84;
const double crRightHandSideBoundedVector333 =             0.027777777777777776*f_ext(0,0);
const double crRightHandSideBoundedVector334 =             0.027777777777777776*f_ext(2,0);
const double crRightHandSideBoundedVector335 =             crRightHandSideBoundedVector308 + crRightHandSideBoundedVector333 + crRightHandSideBoundedVector334;
const double crRightHandSideBoundedVector336 =             -crRightHandSideBoundedVector317*crRightHandSideBoundedVector321 - crRightHandSideBoundedVector323 - crRightHandSideBoundedVector325 + crRightHandSideBoundedVector328 + crRightHandSideBoundedVector330*crRightHandSideBoundedVector331 + crRightHandSideBoundedVector332 + crRightHandSideBoundedVector335;
const double crRightHandSideBoundedVector337 =             0.16666666666666666*dUdt_0_0;
const double crRightHandSideBoundedVector338 =             1.0*(crRightHandSideBoundedVector312 + crRightHandSideBoundedVector337 + 0.66666666666666663*dUdt_1_0)/crRightHandSideBoundedVector116;
const double crRightHandSideBoundedVector339 =             DN_DX_0_1*crRightHandSideBoundedVector141;
const double crRightHandSideBoundedVector340 =             crRightHandSideBoundedVector158*crRightHandSideBoundedVector339;
const double crRightHandSideBoundedVector341 =             crRightHandSideBoundedVector162*crRightHandSideBoundedVector340;
const double crRightHandSideBoundedVector342 =             pow(crRightHandSideBoundedVector140, -3);
const double crRightHandSideBoundedVector343 =             3.0000000000000004*crRightHandSideBoundedVector142;
const double crRightHandSideBoundedVector344 =             crRightHandSideBoundedVector144*crRightHandSideBoundedVector266;
const double crRightHandSideBoundedVector345 =             crRightHandSideBoundedVector342*(-crRightHandSideBoundedVector343 + crRightHandSideBoundedVector344);
const double crRightHandSideBoundedVector346 =             crRightHandSideBoundedVector155*crRightHandSideBoundedVector283;
const double crRightHandSideBoundedVector347 =             crRightHandSideBoundedVector346*crRightHandSideBoundedVector56;
const double crRightHandSideBoundedVector348 =             crRightHandSideBoundedVector158*crRightHandSideBoundedVector283;
const double crRightHandSideBoundedVector349 =             crRightHandSideBoundedVector348*crRightHandSideBoundedVector66;
const double crRightHandSideBoundedVector350 =             1.125*crRightHandSideBoundedVector342;
const double crRightHandSideBoundedVector351 =             crRightHandSideBoundedVector155*crRightHandSideBoundedVector161;
const double crRightHandSideBoundedVector352 =             crRightHandSideBoundedVector350*crRightHandSideBoundedVector351;
const double crRightHandSideBoundedVector353 =             crRightHandSideBoundedVector141*crRightHandSideBoundedVector155;
const double crRightHandSideBoundedVector354 =             crRightHandSideBoundedVector348*crRightHandSideBoundedVector84;
const double crRightHandSideBoundedVector355 =             0.027777777777777776*f_ext(1,0);
const double crRightHandSideBoundedVector356 =             crRightHandSideBoundedVector309 + crRightHandSideBoundedVector333 + crRightHandSideBoundedVector355;
const double crRightHandSideBoundedVector357 =             -crRightHandSideBoundedVector317*crRightHandSideBoundedVector345 + crRightHandSideBoundedVector330*crRightHandSideBoundedVector353 - crRightHandSideBoundedVector347 - crRightHandSideBoundedVector349 + crRightHandSideBoundedVector352 + crRightHandSideBoundedVector354 + crRightHandSideBoundedVector356;
const double crRightHandSideBoundedVector358 =             1.0*(crRightHandSideBoundedVector213 + crRightHandSideBoundedVector311 + crRightHandSideBoundedVector337 + 0.66666666666666663*dUdt_2_0)/crRightHandSideBoundedVector151;
const double crRightHandSideBoundedVector359 =             crRightHandSideBoundedVector77 - 3.0;
const double crRightHandSideBoundedVector360 =             crRightHandSideBoundedVector3*crRightHandSideBoundedVector359;
const double crRightHandSideBoundedVector361 =             DN_DX_0_0*crRightHandSideBoundedVector61;
const double crRightHandSideBoundedVector362 =             0.66666666666666663*crRightHandSideBoundedVector53 + 0.66666666666666663*crRightHandSideBoundedVector54 + 0.66666666666666663*crRightHandSideBoundedVector55;
const double crRightHandSideBoundedVector363 =             DN_DX_0_1*crRightHandSideBoundedVector71;
const double crRightHandSideBoundedVector364 =             1.5000000000000002*crRightHandSideBoundedVector87;
const double crRightHandSideBoundedVector365 =             crRightHandSideBoundedVector10*crRightHandSideBoundedVector364;
const double crRightHandSideBoundedVector366 =             crRightHandSideBoundedVector3*crRightHandSideBoundedVector362 + crRightHandSideBoundedVector363 - crRightHandSideBoundedVector365;
const double crRightHandSideBoundedVector367 =             1.0*crRightHandSideBoundedVector93;
const double crRightHandSideBoundedVector368 =             DN_DX_0_0*crRightHandSideBoundedVector128;
const double crRightHandSideBoundedVector369 =             DN_DX_0_1*crRightHandSideBoundedVector132;
const double crRightHandSideBoundedVector370 =             0.16666666666666666*crRightHandSideBoundedVector53 + 0.16666666666666666*crRightHandSideBoundedVector54 + 0.16666666666666666*crRightHandSideBoundedVector55;
const double crRightHandSideBoundedVector371 =             -crRightHandSideBoundedVector134*crRightHandSideBoundedVector274 + crRightHandSideBoundedVector370*crRightHandSideBoundedVector96;
const double crRightHandSideBoundedVector372 =             crRightHandSideBoundedVector369 + crRightHandSideBoundedVector371;
const double crRightHandSideBoundedVector373 =             crRightHandSideBoundedVector329*crRightHandSideBoundedVector359;
const double crRightHandSideBoundedVector374 =             0.16666666666666666*crRightHandSideBoundedVector300;
const double crRightHandSideBoundedVector375 =             crRightHandSideBoundedVector200*crRightHandSideBoundedVector374 - crRightHandSideBoundedVector373*crRightHandSideBoundedVector96;
const double crRightHandSideBoundedVector376 =             1.0*crRightHandSideBoundedVector137;
const double crRightHandSideBoundedVector377 =             DN_DX_0_0*crRightHandSideBoundedVector156;
const double crRightHandSideBoundedVector378 =             DN_DX_0_1*crRightHandSideBoundedVector159;
const double crRightHandSideBoundedVector379 =             crRightHandSideBoundedVector139*crRightHandSideBoundedVector370 - crRightHandSideBoundedVector161*crRightHandSideBoundedVector283;
const double crRightHandSideBoundedVector380 =             crRightHandSideBoundedVector378 + crRightHandSideBoundedVector379;
const double crRightHandSideBoundedVector381 =             -crRightHandSideBoundedVector139*crRightHandSideBoundedVector373 + crRightHandSideBoundedVector208*crRightHandSideBoundedVector374;
const double crRightHandSideBoundedVector382 =             1.0*crRightHandSideBoundedVector168;
const double crRightHandSideBoundedVector383 =             nu_sc*(1 - crRightHandSideBoundedVector224);
const double crRightHandSideBoundedVector384 =             crRightHandSideBoundedVector2*crRightHandSideBoundedVector83;
const double crRightHandSideBoundedVector385 =             crRightHandSideBoundedVector2*crRightHandSideBoundedVector66;
const double crRightHandSideBoundedVector386 =             (crRightHandSideBoundedVector2*crRightHandSideBoundedVector225 + crRightHandSideBoundedVector2*crRightHandSideBoundedVector383 + mu)*(-2.2500000000000004*crRightHandSideBoundedVector263 - 2.2500000000000004*crRightHandSideBoundedVector265 + 2.2500000000000004*crRightHandSideBoundedVector384 + 2.2500000000000004*crRightHandSideBoundedVector385);
const double crRightHandSideBoundedVector387 =             crRightHandSideBoundedVector83*crRightHandSideBoundedVector95;
const double crRightHandSideBoundedVector388 =             crRightHandSideBoundedVector66*crRightHandSideBoundedVector95;
const double crRightHandSideBoundedVector389 =             (crRightHandSideBoundedVector225*crRightHandSideBoundedVector95 + crRightHandSideBoundedVector383*crRightHandSideBoundedVector95 + mu)*(-2.2500000000000004*crRightHandSideBoundedVector273 - 2.2500000000000004*crRightHandSideBoundedVector276 + 2.2500000000000004*crRightHandSideBoundedVector387 + 2.2500000000000004*crRightHandSideBoundedVector388);
const double crRightHandSideBoundedVector390 =             crRightHandSideBoundedVector138*crRightHandSideBoundedVector83;
const double crRightHandSideBoundedVector391 =             crRightHandSideBoundedVector138*crRightHandSideBoundedVector66;
const double crRightHandSideBoundedVector392 =             (crRightHandSideBoundedVector138*crRightHandSideBoundedVector225 + crRightHandSideBoundedVector138*crRightHandSideBoundedVector383 + mu)*(-2.2500000000000004*crRightHandSideBoundedVector282 - 2.2500000000000004*crRightHandSideBoundedVector285 + 2.2500000000000004*crRightHandSideBoundedVector390 + 2.2500000000000004*crRightHandSideBoundedVector391);
const double crRightHandSideBoundedVector393 =             0.66666666666666663*crRightHandSideBoundedVector84;
const double crRightHandSideBoundedVector394 =             crRightHandSideBoundedVector299*crRightHandSideBoundedVector359;
const double crRightHandSideBoundedVector395 =             DN_DX_0_0*crRightHandSideBoundedVector262;
const double crRightHandSideBoundedVector396 =             lambda/c_v;
const double crRightHandSideBoundedVector397 =             crRightHandSideBoundedVector38*crRightHandSideBoundedVector396;
const double crRightHandSideBoundedVector398 =             0.16666666666666666*dUdt_1_3;
const double crRightHandSideBoundedVector399 =             0.16666666666666666*dUdt_2_3;
const double crRightHandSideBoundedVector400 =             0.16666666666666666*r_ext[1];
const double crRightHandSideBoundedVector401 =             0.16666666666666666*r_ext[2];
const double crRightHandSideBoundedVector402 =             crRightHandSideBoundedVector2*(crRightHandSideBoundedVector400 + crRightHandSideBoundedVector401 + 0.66666666666666663*r_ext[0]);
const double crRightHandSideBoundedVector403 =             crRightHandSideBoundedVector45*crRightHandSideBoundedVector60;
const double crRightHandSideBoundedVector404 =             crRightHandSideBoundedVector172*crRightHandSideBoundedVector70;
const double crRightHandSideBoundedVector405 =             crRightHandSideBoundedVector90*gamma;
const double crRightHandSideBoundedVector406 =             crRightHandSideBoundedVector405*crRightHandSideBoundedVector61;
const double crRightHandSideBoundedVector407 =             crRightHandSideBoundedVector185*gamma;
const double crRightHandSideBoundedVector408 =             crRightHandSideBoundedVector407*crRightHandSideBoundedVector71;
const double crRightHandSideBoundedVector409 =             crRightHandSideBoundedVector70*crRightHandSideBoundedVector89;
const double crRightHandSideBoundedVector410 =             crRightHandSideBoundedVector66*crRightHandSideBoundedVector88;
const double crRightHandSideBoundedVector411 =             crRightHandSideBoundedVector10*crRightHandSideBoundedVector19;
const double crRightHandSideBoundedVector412 =             crRightHandSideBoundedVector411*crRightHandSideBoundedVector70;
const double crRightHandSideBoundedVector413 =             crRightHandSideBoundedVector411*crRightHandSideBoundedVector48;
const double crRightHandSideBoundedVector414 =             crRightHandSideBoundedVector19*(crRightHandSideBoundedVector23 - crRightHandSideBoundedVector3*(0.22222222222222221*crRightHandSideBoundedVector13 + 0.22222222222222221*crRightHandSideBoundedVector16));
const double crRightHandSideBoundedVector415 =             crRightHandSideBoundedVector3*(crRightHandSideBoundedVector23 + crRightHandSideBoundedVector414);
const double crRightHandSideBoundedVector416 =             crRightHandSideBoundedVector174*crRightHandSideBoundedVector411;
const double crRightHandSideBoundedVector417 =             -0.37500000000000006*U_1_3;
const double crRightHandSideBoundedVector418 =             -0.37500000000000006*U_2_3;
const double crRightHandSideBoundedVector419 =             0.75000000000000011*crRightHandSideBoundedVector19;
const double crRightHandSideBoundedVector420 =             1.0/crRightHandSideBoundedVector9;
const double crRightHandSideBoundedVector421 =             crRightHandSideBoundedVector17*crRightHandSideBoundedVector420;
const double crRightHandSideBoundedVector422 =             -1.5000000000000002*U_0_3 - 2.2500000000000004*crRightHandSideBoundedVector414 + crRightHandSideBoundedVector417 + crRightHandSideBoundedVector418 + crRightHandSideBoundedVector419*crRightHandSideBoundedVector421;
const double crRightHandSideBoundedVector423 =             crRightHandSideBoundedVector10*crRightHandSideBoundedVector422;
const double crRightHandSideBoundedVector424 =             crRightHandSideBoundedVector183*crRightHandSideBoundedVector423;
const double crRightHandSideBoundedVector425 =             crRightHandSideBoundedVector423*crRightHandSideBoundedVector87;
const double crRightHandSideBoundedVector426 =             (crRightHandSideBoundedVector398 + crRightHandSideBoundedVector399 - crRightHandSideBoundedVector402 - crRightHandSideBoundedVector403 - crRightHandSideBoundedVector404 + crRightHandSideBoundedVector406 + crRightHandSideBoundedVector408 - crRightHandSideBoundedVector409*crRightHandSideBoundedVector84 - crRightHandSideBoundedVector410*crRightHandSideBoundedVector412 + crRightHandSideBoundedVector424 + crRightHandSideBoundedVector425 + crRightHandSideBoundedVector56*(crRightHandSideBoundedVector415 - crRightHandSideBoundedVector416) + crRightHandSideBoundedVector76*(-crRightHandSideBoundedVector413 + crRightHandSideBoundedVector415) + 0.66666666666666663*dUdt_0_3)/(crRightHandSideBoundedVector397*crRightHandSideBoundedVector5 + crRightHandSideBoundedVector40);
const double crRightHandSideBoundedVector427 =             crRightHandSideBoundedVector397*crRightHandSideBoundedVector4;
const double crRightHandSideBoundedVector428 =             0.16666666666666666*dUdt_0_3;
const double crRightHandSideBoundedVector429 =             0.16666666666666666*r_ext[0];
const double crRightHandSideBoundedVector430 =             crRightHandSideBoundedVector95*(crRightHandSideBoundedVector401 + crRightHandSideBoundedVector429 + 0.66666666666666663*r_ext[1]);
const double crRightHandSideBoundedVector431 =             crRightHandSideBoundedVector120*crRightHandSideBoundedVector127;
const double crRightHandSideBoundedVector432 =             crRightHandSideBoundedVector131*crRightHandSideBoundedVector191;
const double crRightHandSideBoundedVector433 =             crRightHandSideBoundedVector128*crRightHandSideBoundedVector405;
const double crRightHandSideBoundedVector434 =             crRightHandSideBoundedVector132*crRightHandSideBoundedVector407;
const double crRightHandSideBoundedVector435 =             crRightHandSideBoundedVector131*crRightHandSideBoundedVector136;
const double crRightHandSideBoundedVector436 =             crRightHandSideBoundedVector135*crRightHandSideBoundedVector66;
const double crRightHandSideBoundedVector437 =             crRightHandSideBoundedVector100*crRightHandSideBoundedVector19;
const double crRightHandSideBoundedVector438 =             crRightHandSideBoundedVector131*crRightHandSideBoundedVector437;
const double crRightHandSideBoundedVector439 =             crRightHandSideBoundedVector122*crRightHandSideBoundedVector437;
const double crRightHandSideBoundedVector440 =             crRightHandSideBoundedVector19*(crRightHandSideBoundedVector108 - crRightHandSideBoundedVector96*(0.22222222222222221*crRightHandSideBoundedVector102 + 0.22222222222222221*crRightHandSideBoundedVector104));
const double crRightHandSideBoundedVector441 =             crRightHandSideBoundedVector96*(crRightHandSideBoundedVector108 + crRightHandSideBoundedVector440);
const double crRightHandSideBoundedVector442 =             crRightHandSideBoundedVector193*crRightHandSideBoundedVector437;
const double crRightHandSideBoundedVector443 =             -0.37500000000000006*U_0_3;
const double crRightHandSideBoundedVector444 =             1.0/crRightHandSideBoundedVector99;
const double crRightHandSideBoundedVector445 =             crRightHandSideBoundedVector105*crRightHandSideBoundedVector444;
const double crRightHandSideBoundedVector446 =             -1.5000000000000002*U_1_3 + crRightHandSideBoundedVector418 + crRightHandSideBoundedVector419*crRightHandSideBoundedVector445 - 2.2500000000000004*crRightHandSideBoundedVector440 + crRightHandSideBoundedVector443;
const double crRightHandSideBoundedVector447 =             crRightHandSideBoundedVector100*crRightHandSideBoundedVector446;
const double crRightHandSideBoundedVector448 =             crRightHandSideBoundedVector199*crRightHandSideBoundedVector447;
const double crRightHandSideBoundedVector449 =             crRightHandSideBoundedVector134*crRightHandSideBoundedVector447;
const double crRightHandSideBoundedVector450 =             (crRightHandSideBoundedVector399 + crRightHandSideBoundedVector428 - crRightHandSideBoundedVector430 - crRightHandSideBoundedVector431 - crRightHandSideBoundedVector432 + crRightHandSideBoundedVector433 + crRightHandSideBoundedVector434 - crRightHandSideBoundedVector435*crRightHandSideBoundedVector84 - crRightHandSideBoundedVector436*crRightHandSideBoundedVector438 + crRightHandSideBoundedVector448 + crRightHandSideBoundedVector449 + crRightHandSideBoundedVector56*(crRightHandSideBoundedVector441 - crRightHandSideBoundedVector442) + crRightHandSideBoundedVector76*(-crRightHandSideBoundedVector439 + crRightHandSideBoundedVector441) + 0.66666666666666663*dUdt_1_3)/(crRightHandSideBoundedVector116 + crRightHandSideBoundedVector427*crRightHandSideBoundedVector96);
const double crRightHandSideBoundedVector451 =             crRightHandSideBoundedVector138*(crRightHandSideBoundedVector400 + crRightHandSideBoundedVector429 + 0.66666666666666663*r_ext[2]);
const double crRightHandSideBoundedVector452 =             crRightHandSideBoundedVector153*crRightHandSideBoundedVector155;
const double crRightHandSideBoundedVector453 =             crRightHandSideBoundedVector158*crRightHandSideBoundedVector202;
const double crRightHandSideBoundedVector454 =             crRightHandSideBoundedVector156*crRightHandSideBoundedVector405;
const double crRightHandSideBoundedVector455 =             crRightHandSideBoundedVector159*crRightHandSideBoundedVector407;
const double crRightHandSideBoundedVector456 =             crRightHandSideBoundedVector158*crRightHandSideBoundedVector163;
const double crRightHandSideBoundedVector457 =             crRightHandSideBoundedVector162*crRightHandSideBoundedVector66;
const double crRightHandSideBoundedVector458 =             crRightHandSideBoundedVector141*crRightHandSideBoundedVector19;
const double crRightHandSideBoundedVector459 =             crRightHandSideBoundedVector158*crRightHandSideBoundedVector458;
const double crRightHandSideBoundedVector460 =             crRightHandSideBoundedVector164*crRightHandSideBoundedVector458;
const double crRightHandSideBoundedVector461 =             crRightHandSideBoundedVector19*(-crRightHandSideBoundedVector139*(0.22222222222222221*crRightHandSideBoundedVector142 + 0.22222222222222221*crRightHandSideBoundedVector143) + crRightHandSideBoundedVector146);
const double crRightHandSideBoundedVector462 =             crRightHandSideBoundedVector139*(crRightHandSideBoundedVector146 + crRightHandSideBoundedVector461);
const double crRightHandSideBoundedVector463 =             crRightHandSideBoundedVector209*crRightHandSideBoundedVector458;
const double crRightHandSideBoundedVector464 =             1.0/crRightHandSideBoundedVector140;
const double crRightHandSideBoundedVector465 =             crRightHandSideBoundedVector144*crRightHandSideBoundedVector464;
const double crRightHandSideBoundedVector466 =             -1.5000000000000002*U_2_3 + crRightHandSideBoundedVector417 + crRightHandSideBoundedVector419*crRightHandSideBoundedVector465 + crRightHandSideBoundedVector443 - 2.2500000000000004*crRightHandSideBoundedVector461;
const double crRightHandSideBoundedVector467 =             crRightHandSideBoundedVector141*crRightHandSideBoundedVector466;
const double crRightHandSideBoundedVector468 =             crRightHandSideBoundedVector207*crRightHandSideBoundedVector467;
const double crRightHandSideBoundedVector469 =             crRightHandSideBoundedVector161*crRightHandSideBoundedVector467;
const double crRightHandSideBoundedVector470 =             (crRightHandSideBoundedVector398 + crRightHandSideBoundedVector428 - crRightHandSideBoundedVector451 - crRightHandSideBoundedVector452 - crRightHandSideBoundedVector453 + crRightHandSideBoundedVector454 + crRightHandSideBoundedVector455 - crRightHandSideBoundedVector456*crRightHandSideBoundedVector84 - crRightHandSideBoundedVector457*crRightHandSideBoundedVector459 + crRightHandSideBoundedVector468 + crRightHandSideBoundedVector469 + crRightHandSideBoundedVector56*(crRightHandSideBoundedVector462 - crRightHandSideBoundedVector463) + crRightHandSideBoundedVector76*(-crRightHandSideBoundedVector460 + crRightHandSideBoundedVector462) + 0.66666666666666663*dUdt_2_3)/(crRightHandSideBoundedVector139*crRightHandSideBoundedVector427 + crRightHandSideBoundedVector151);
const double crRightHandSideBoundedVector471 =             0.16666666666666666*crRightHandSideBoundedVector154 + crRightHandSideBoundedVector156*crRightHandSideBoundedVector373 - 0.16666666666666666*crRightHandSideBoundedVector157 + crRightHandSideBoundedVector159*crRightHandSideBoundedVector272 - 0.16666666666666666*crRightHandSideBoundedVector160 + crRightHandSideBoundedVector161*crRightHandSideBoundedVector346 - 0.16666666666666666*crRightHandSideBoundedVector167;
const double crRightHandSideBoundedVector472 =             0.16666666666666666*crRightHandSideBoundedVector121 - 0.16666666666666666*crRightHandSideBoundedVector125 + crRightHandSideBoundedVector128*crRightHandSideBoundedVector373 - 0.16666666666666666*crRightHandSideBoundedVector129 + crRightHandSideBoundedVector132*crRightHandSideBoundedVector272 - 0.16666666666666666*crRightHandSideBoundedVector133 + crRightHandSideBoundedVector134*crRightHandSideBoundedVector322;
const double crRightHandSideBoundedVector473 =             crRightHandSideBoundedVector223*pow(lin_m[1], 2);
const double crRightHandSideBoundedVector474 =             crRightHandSideBoundedVector473*nu_st;
const double crRightHandSideBoundedVector475 =             1 - crRightHandSideBoundedVector473;
const double crRightHandSideBoundedVector476 =             crRightHandSideBoundedVector475*nu_sc;
const double crRightHandSideBoundedVector477 =             crRightHandSideBoundedVector217*crRightHandSideBoundedVector219 + crRightHandSideBoundedVector220 + crRightHandSideBoundedVector222*crRightHandSideBoundedVector474 + crRightHandSideBoundedVector222*crRightHandSideBoundedVector476 + crRightHandSideBoundedVector225*crRightHandSideBoundedVector227 - crRightHandSideBoundedVector227*crRightHandSideBoundedVector230;
const double crRightHandSideBoundedVector478 =             crRightHandSideBoundedVector225*crRightHandSideBoundedVector244 - crRightHandSideBoundedVector230*crRightHandSideBoundedVector244 + crRightHandSideBoundedVector237*crRightHandSideBoundedVector239 + crRightHandSideBoundedVector240 + crRightHandSideBoundedVector242*crRightHandSideBoundedVector474 + crRightHandSideBoundedVector242*crRightHandSideBoundedVector476;
const double crRightHandSideBoundedVector479 =             crRightHandSideBoundedVector225*crRightHandSideBoundedVector256 - crRightHandSideBoundedVector230*crRightHandSideBoundedVector256 + crRightHandSideBoundedVector249*crRightHandSideBoundedVector251 + crRightHandSideBoundedVector252 + crRightHandSideBoundedVector254*crRightHandSideBoundedVector474 + crRightHandSideBoundedVector254*crRightHandSideBoundedVector476;
const double crRightHandSideBoundedVector480 =             3.0000000000000004*crRightHandSideBoundedVector16;
const double crRightHandSideBoundedVector481 =             crRightHandSideBoundedVector86*(crRightHandSideBoundedVector292 - crRightHandSideBoundedVector480);
const double crRightHandSideBoundedVector482 =             crRightHandSideBoundedVector183*crRightHandSideBoundedVector70;
const double crRightHandSideBoundedVector483 =             crRightHandSideBoundedVector296*crRightHandSideBoundedVector482;
const double crRightHandSideBoundedVector484 =             0.66666666666666663*crRightHandSideBoundedVector56;
const double crRightHandSideBoundedVector485 =             crRightHandSideBoundedVector301*crRightHandSideBoundedVector70;
const double crRightHandSideBoundedVector486 =             crRightHandSideBoundedVector267*crRightHandSideBoundedVector60;
const double crRightHandSideBoundedVector487 =             crRightHandSideBoundedVector486*crRightHandSideBoundedVector66;
const double crRightHandSideBoundedVector488 =             0.1111111111111111*f_ext(1,1);
const double crRightHandSideBoundedVector489 =             0.1111111111111111*f_ext(2,1);
const double crRightHandSideBoundedVector490 =             crRightHandSideBoundedVector488 + crRightHandSideBoundedVector489 + 0.44444444444444442*f_ext(0,1);
const double crRightHandSideBoundedVector491 =             0.16666666666666666*crRightHandSideBoundedVector86;
const double crRightHandSideBoundedVector492 =             3.0000000000000004*crRightHandSideBoundedVector104;
const double crRightHandSideBoundedVector493 =             crRightHandSideBoundedVector318*(crRightHandSideBoundedVector320 - crRightHandSideBoundedVector492);
const double crRightHandSideBoundedVector494 =             crRightHandSideBoundedVector324*crRightHandSideBoundedVector76;
const double crRightHandSideBoundedVector495 =             crRightHandSideBoundedVector131*crRightHandSideBoundedVector199;
const double crRightHandSideBoundedVector496 =             crRightHandSideBoundedVector326*crRightHandSideBoundedVector495;
const double crRightHandSideBoundedVector497 =             0.16666666666666666*crRightHandSideBoundedVector56;
const double crRightHandSideBoundedVector498 =             crRightHandSideBoundedVector300*crRightHandSideBoundedVector497;
const double crRightHandSideBoundedVector499 =             crRightHandSideBoundedVector100*crRightHandSideBoundedVector131;
const double crRightHandSideBoundedVector500 =             crRightHandSideBoundedVector180*crRightHandSideBoundedVector322;
const double crRightHandSideBoundedVector501 =             0.027777777777777776*f_ext(0,1);
const double crRightHandSideBoundedVector502 =             0.027777777777777776*f_ext(2,1);
const double crRightHandSideBoundedVector503 =             crRightHandSideBoundedVector488 + crRightHandSideBoundedVector501 + crRightHandSideBoundedVector502;
const double crRightHandSideBoundedVector504 =             -crRightHandSideBoundedVector322*crRightHandSideBoundedVector83 - crRightHandSideBoundedVector491*crRightHandSideBoundedVector493 - crRightHandSideBoundedVector494 + crRightHandSideBoundedVector496 + crRightHandSideBoundedVector498*crRightHandSideBoundedVector499 + crRightHandSideBoundedVector500 + crRightHandSideBoundedVector503;
const double crRightHandSideBoundedVector505 =             3.0000000000000004*crRightHandSideBoundedVector143;
const double crRightHandSideBoundedVector506 =             crRightHandSideBoundedVector342*(crRightHandSideBoundedVector344 - crRightHandSideBoundedVector505);
const double crRightHandSideBoundedVector507 =             crRightHandSideBoundedVector348*crRightHandSideBoundedVector76;
const double crRightHandSideBoundedVector508 =             crRightHandSideBoundedVector158*crRightHandSideBoundedVector207;
const double crRightHandSideBoundedVector509 =             crRightHandSideBoundedVector350*crRightHandSideBoundedVector508;
const double crRightHandSideBoundedVector510 =             crRightHandSideBoundedVector141*crRightHandSideBoundedVector158;
const double crRightHandSideBoundedVector511 =             crRightHandSideBoundedVector180*crRightHandSideBoundedVector346;
const double crRightHandSideBoundedVector512 =             0.027777777777777776*f_ext(1,1);
const double crRightHandSideBoundedVector513 =             crRightHandSideBoundedVector489 + crRightHandSideBoundedVector501 + crRightHandSideBoundedVector512;
const double crRightHandSideBoundedVector514 =             -crRightHandSideBoundedVector346*crRightHandSideBoundedVector83 - crRightHandSideBoundedVector491*crRightHandSideBoundedVector506 + crRightHandSideBoundedVector498*crRightHandSideBoundedVector510 - crRightHandSideBoundedVector507 + crRightHandSideBoundedVector509 + crRightHandSideBoundedVector511 + crRightHandSideBoundedVector513;
const double crRightHandSideBoundedVector515 =             0.66666666666666663*crRightHandSideBoundedVector73 + 0.66666666666666663*crRightHandSideBoundedVector74 + 0.66666666666666663*crRightHandSideBoundedVector75;
const double crRightHandSideBoundedVector516 =             1.5000000000000002*crRightHandSideBoundedVector183;
const double crRightHandSideBoundedVector517 =             crRightHandSideBoundedVector10*crRightHandSideBoundedVector516;
const double crRightHandSideBoundedVector518 =             crRightHandSideBoundedVector3*crRightHandSideBoundedVector515 + crRightHandSideBoundedVector361 - crRightHandSideBoundedVector517;
const double crRightHandSideBoundedVector519 =             DN_DX_0_0*crRightHandSideBoundedVector67 + 0.66666666666666663*crRightHandSideBoundedVector81 + 0.66666666666666663*crRightHandSideBoundedVector82;
const double crRightHandSideBoundedVector520 =             0.16666666666666666*crRightHandSideBoundedVector73 + 0.16666666666666666*crRightHandSideBoundedVector74 + 0.16666666666666666*crRightHandSideBoundedVector75;
const double crRightHandSideBoundedVector521 =             -crRightHandSideBoundedVector199*crRightHandSideBoundedVector274 + crRightHandSideBoundedVector520*crRightHandSideBoundedVector96;
const double crRightHandSideBoundedVector522 =             crRightHandSideBoundedVector368 + crRightHandSideBoundedVector521;
const double crRightHandSideBoundedVector523 =             crRightHandSideBoundedVector359*crRightHandSideBoundedVector497;
const double crRightHandSideBoundedVector524 =             crRightHandSideBoundedVector100*crRightHandSideBoundedVector134*crRightHandSideBoundedVector374 - crRightHandSideBoundedVector523*crRightHandSideBoundedVector96;
const double crRightHandSideBoundedVector525 =             DN_DX_1_0*crRightHandSideBoundedVector68 + DN_DX_2_0*crRightHandSideBoundedVector69 + 0.16666666666666666*crRightHandSideBoundedVector80;
const double crRightHandSideBoundedVector526 =             0.16666666666666666*crRightHandSideBoundedVector180;
const double crRightHandSideBoundedVector527 =             crRightHandSideBoundedVector19*crRightHandSideBoundedVector275 - crRightHandSideBoundedVector277 + crRightHandSideBoundedVector525*crRightHandSideBoundedVector96 - crRightHandSideBoundedVector526*crRightHandSideBoundedVector96;
const double crRightHandSideBoundedVector528 =             crRightHandSideBoundedVector139*crRightHandSideBoundedVector520 - crRightHandSideBoundedVector207*crRightHandSideBoundedVector283;
const double crRightHandSideBoundedVector529 =             crRightHandSideBoundedVector377 + crRightHandSideBoundedVector528;
const double crRightHandSideBoundedVector530 =             -crRightHandSideBoundedVector139*crRightHandSideBoundedVector523 + crRightHandSideBoundedVector141*crRightHandSideBoundedVector161*crRightHandSideBoundedVector374;
const double crRightHandSideBoundedVector531 =             crRightHandSideBoundedVector139*crRightHandSideBoundedVector525 - crRightHandSideBoundedVector139*crRightHandSideBoundedVector526 + crRightHandSideBoundedVector19*crRightHandSideBoundedVector284 - crRightHandSideBoundedVector286;
const double crRightHandSideBoundedVector532 =             crRightHandSideBoundedVector10*crRightHandSideBoundedVector386;
const double crRightHandSideBoundedVector533 =             crRightHandSideBoundedVector100*crRightHandSideBoundedVector389;
const double crRightHandSideBoundedVector534 =             crRightHandSideBoundedVector141*crRightHandSideBoundedVector392;
const double crRightHandSideBoundedVector535 =             0.66666666666666663*crRightHandSideBoundedVector180;
const double crRightHandSideBoundedVector536 =             crRightHandSideBoundedVector359*crRightHandSideBoundedVector484;
const double crRightHandSideBoundedVector537 =             DN_DX_0_1*crRightHandSideBoundedVector262;
const double crRightHandSideBoundedVector538 =             crRightHandSideBoundedVector156*crRightHandSideBoundedVector526 + crRightHandSideBoundedVector159*crRightHandSideBoundedVector523 + 0.16666666666666666*crRightHandSideBoundedVector203 - 0.16666666666666666*crRightHandSideBoundedVector204 - 0.16666666666666666*crRightHandSideBoundedVector205 + crRightHandSideBoundedVector207*crRightHandSideBoundedVector348 - 0.16666666666666666*crRightHandSideBoundedVector211;
const double crRightHandSideBoundedVector539 =             crRightHandSideBoundedVector128*crRightHandSideBoundedVector526 + crRightHandSideBoundedVector132*crRightHandSideBoundedVector523 + 0.16666666666666666*crRightHandSideBoundedVector192 - 0.16666666666666666*crRightHandSideBoundedVector195 - 0.16666666666666666*crRightHandSideBoundedVector196 - 0.16666666666666666*crRightHandSideBoundedVector197 + crRightHandSideBoundedVector199*crRightHandSideBoundedVector324;
const double crRightHandSideBoundedVector540 =             -crRightHandSideBoundedVector415;
const double crRightHandSideBoundedVector541 =             crRightHandSideBoundedVector413 + crRightHandSideBoundedVector540;
const double crRightHandSideBoundedVector542 =             crRightHandSideBoundedVector416 + crRightHandSideBoundedVector540;
const double crRightHandSideBoundedVector543 =             crRightHandSideBoundedVector71*mu;
const double crRightHandSideBoundedVector544 =             -2.2500000000000004*crRightHandSideBoundedVector263 - 2.2500000000000004*crRightHandSideBoundedVector265 + 2.2500000000000004*crRightHandSideBoundedVector384 + 2.2500000000000004*crRightHandSideBoundedVector385;
const double crRightHandSideBoundedVector545 =             crRightHandSideBoundedVector61*mu;
const double crRightHandSideBoundedVector546 =             2.2500000000000004*crRightHandSideBoundedVector2;
const double crRightHandSideBoundedVector547 =             2.2500000000000004*crRightHandSideBoundedVector23;
const double crRightHandSideBoundedVector548 =             1.5000000000000002*crRightHandSideBoundedVector420;
const double crRightHandSideBoundedVector549 =             crRightHandSideBoundedVector47*crRightHandSideBoundedVector548;
const double crRightHandSideBoundedVector550 =             crRightHandSideBoundedVector13*crRightHandSideBoundedVector549 + crRightHandSideBoundedVector16*crRightHandSideBoundedVector549 - crRightHandSideBoundedVector182*crRightHandSideBoundedVector83 - crRightHandSideBoundedVector47*crRightHandSideBoundedVector547 + crRightHandSideBoundedVector546*crRightHandSideBoundedVector90 - crRightHandSideBoundedVector76*crRightHandSideBoundedVector88;
const double crRightHandSideBoundedVector551 =             gamma*k_st;
const double crRightHandSideBoundedVector552 =             crRightHandSideBoundedVector548*crRightHandSideBoundedVector86;
const double crRightHandSideBoundedVector553 =             crRightHandSideBoundedVector13*crRightHandSideBoundedVector552 + crRightHandSideBoundedVector16*crRightHandSideBoundedVector552 - crRightHandSideBoundedVector182*crRightHandSideBoundedVector56 + crRightHandSideBoundedVector185*crRightHandSideBoundedVector546 - crRightHandSideBoundedVector410 - crRightHandSideBoundedVector547*crRightHandSideBoundedVector86;
const double crRightHandSideBoundedVector554 =             crRightHandSideBoundedVector2*crRightHandSideBoundedVector553;
const double crRightHandSideBoundedVector555 =             crRightHandSideBoundedVector224*crRightHandSideBoundedVector554;
const double crRightHandSideBoundedVector556 =             crRightHandSideBoundedVector2*crRightHandSideBoundedVector550;
const double crRightHandSideBoundedVector557 =             crRightHandSideBoundedVector228*crRightHandSideBoundedVector551;
const double crRightHandSideBoundedVector558 =             gamma*k_sc;
const double crRightHandSideBoundedVector559 =             crRightHandSideBoundedVector231*crRightHandSideBoundedVector558;
const double crRightHandSideBoundedVector560 =             crRightHandSideBoundedVector10*(crRightHandSideBoundedVector396*crRightHandSideBoundedVector550 + crRightHandSideBoundedVector543*crRightHandSideBoundedVector544 + crRightHandSideBoundedVector545*(-3.0000000000000009*crRightHandSideBoundedVector183 + 3.0000000000000009*crRightHandSideBoundedVector214 - 1.5000000000000002*crRightHandSideBoundedVector218 + crRightHandSideBoundedVector364) + crRightHandSideBoundedVector551*crRightHandSideBoundedVector555 - crRightHandSideBoundedVector555*crRightHandSideBoundedVector558 + crRightHandSideBoundedVector556*crRightHandSideBoundedVector557 + crRightHandSideBoundedVector556*crRightHandSideBoundedVector559);
const double crRightHandSideBoundedVector561 =             crRightHandSideBoundedVector132*mu;
const double crRightHandSideBoundedVector562 =             -2.2500000000000004*crRightHandSideBoundedVector273 - 2.2500000000000004*crRightHandSideBoundedVector276 + 2.2500000000000004*crRightHandSideBoundedVector387 + 2.2500000000000004*crRightHandSideBoundedVector388;
const double crRightHandSideBoundedVector563 =             1.5000000000000002*crRightHandSideBoundedVector134;
const double crRightHandSideBoundedVector564 =             crRightHandSideBoundedVector128*mu;
const double crRightHandSideBoundedVector565 =             2.2500000000000004*crRightHandSideBoundedVector95;
const double crRightHandSideBoundedVector566 =             2.2500000000000004*crRightHandSideBoundedVector108;
const double crRightHandSideBoundedVector567 =             1.5000000000000002*crRightHandSideBoundedVector444;
const double crRightHandSideBoundedVector568 =             crRightHandSideBoundedVector47*crRightHandSideBoundedVector567;
const double crRightHandSideBoundedVector569 =             crRightHandSideBoundedVector102*crRightHandSideBoundedVector568 + crRightHandSideBoundedVector104*crRightHandSideBoundedVector568 - crRightHandSideBoundedVector135*crRightHandSideBoundedVector76 - crRightHandSideBoundedVector198*crRightHandSideBoundedVector83 - crRightHandSideBoundedVector47*crRightHandSideBoundedVector566 + crRightHandSideBoundedVector565*crRightHandSideBoundedVector90;
const double crRightHandSideBoundedVector570 =             crRightHandSideBoundedVector567*crRightHandSideBoundedVector86;
const double crRightHandSideBoundedVector571 =             crRightHandSideBoundedVector102*crRightHandSideBoundedVector570 + crRightHandSideBoundedVector104*crRightHandSideBoundedVector570 + crRightHandSideBoundedVector185*crRightHandSideBoundedVector565 - crRightHandSideBoundedVector198*crRightHandSideBoundedVector56 - crRightHandSideBoundedVector436 - crRightHandSideBoundedVector566*crRightHandSideBoundedVector86;
const double crRightHandSideBoundedVector572 =             crRightHandSideBoundedVector571*crRightHandSideBoundedVector95;
const double crRightHandSideBoundedVector573 =             crRightHandSideBoundedVector224*crRightHandSideBoundedVector572;
const double crRightHandSideBoundedVector574 =             crRightHandSideBoundedVector569*crRightHandSideBoundedVector95;
const double crRightHandSideBoundedVector575 =             crRightHandSideBoundedVector100*(crRightHandSideBoundedVector396*crRightHandSideBoundedVector569 + crRightHandSideBoundedVector551*crRightHandSideBoundedVector573 + crRightHandSideBoundedVector557*crRightHandSideBoundedVector574 - crRightHandSideBoundedVector558*crRightHandSideBoundedVector573 + crRightHandSideBoundedVector559*crRightHandSideBoundedVector574 + crRightHandSideBoundedVector561*crRightHandSideBoundedVector562 + crRightHandSideBoundedVector564*(-3.0000000000000009*crRightHandSideBoundedVector199 + 3.0000000000000009*crRightHandSideBoundedVector234 - 1.5000000000000002*crRightHandSideBoundedVector238 + crRightHandSideBoundedVector563));
const double crRightHandSideBoundedVector576 =             crRightHandSideBoundedVector159*mu;
const double crRightHandSideBoundedVector577 =             -2.2500000000000004*crRightHandSideBoundedVector282 - 2.2500000000000004*crRightHandSideBoundedVector285 + 2.2500000000000004*crRightHandSideBoundedVector390 + 2.2500000000000004*crRightHandSideBoundedVector391;
const double crRightHandSideBoundedVector578 =             1.5000000000000002*crRightHandSideBoundedVector161;
const double crRightHandSideBoundedVector579 =             crRightHandSideBoundedVector156*mu;
const double crRightHandSideBoundedVector580 =             2.2500000000000004*crRightHandSideBoundedVector138;
const double crRightHandSideBoundedVector581 =             2.2500000000000004*crRightHandSideBoundedVector146;
const double crRightHandSideBoundedVector582 =             1.5000000000000002*crRightHandSideBoundedVector464;
const double crRightHandSideBoundedVector583 =             crRightHandSideBoundedVector47*crRightHandSideBoundedVector582;
const double crRightHandSideBoundedVector584 =             crRightHandSideBoundedVector142*crRightHandSideBoundedVector583 + crRightHandSideBoundedVector143*crRightHandSideBoundedVector583 - crRightHandSideBoundedVector162*crRightHandSideBoundedVector76 - crRightHandSideBoundedVector206*crRightHandSideBoundedVector83 - crRightHandSideBoundedVector47*crRightHandSideBoundedVector581 + crRightHandSideBoundedVector580*crRightHandSideBoundedVector90;
const double crRightHandSideBoundedVector585 =             crRightHandSideBoundedVector582*crRightHandSideBoundedVector86;
const double crRightHandSideBoundedVector586 =             crRightHandSideBoundedVector142*crRightHandSideBoundedVector585 + crRightHandSideBoundedVector143*crRightHandSideBoundedVector585 + crRightHandSideBoundedVector185*crRightHandSideBoundedVector580 - crRightHandSideBoundedVector206*crRightHandSideBoundedVector56 - crRightHandSideBoundedVector457 - crRightHandSideBoundedVector581*crRightHandSideBoundedVector86;
const double crRightHandSideBoundedVector587 =             crRightHandSideBoundedVector138*crRightHandSideBoundedVector586;
const double crRightHandSideBoundedVector588 =             crRightHandSideBoundedVector224*crRightHandSideBoundedVector587;
const double crRightHandSideBoundedVector589 =             crRightHandSideBoundedVector138*crRightHandSideBoundedVector584;
const double crRightHandSideBoundedVector590 =             crRightHandSideBoundedVector141*(crRightHandSideBoundedVector396*crRightHandSideBoundedVector584 + crRightHandSideBoundedVector551*crRightHandSideBoundedVector588 + crRightHandSideBoundedVector557*crRightHandSideBoundedVector589 - crRightHandSideBoundedVector558*crRightHandSideBoundedVector588 + crRightHandSideBoundedVector559*crRightHandSideBoundedVector589 + crRightHandSideBoundedVector576*crRightHandSideBoundedVector577 + crRightHandSideBoundedVector579*(-3.0000000000000009*crRightHandSideBoundedVector207 + 3.0000000000000009*crRightHandSideBoundedVector246 - 1.5000000000000002*crRightHandSideBoundedVector250 + crRightHandSideBoundedVector578));
const double crRightHandSideBoundedVector591 =             crRightHandSideBoundedVector224*crRightHandSideBoundedVector556;
const double crRightHandSideBoundedVector592 =             crRightHandSideBoundedVector473*crRightHandSideBoundedVector551;
const double crRightHandSideBoundedVector593 =             crRightHandSideBoundedVector475*crRightHandSideBoundedVector558;
const double crRightHandSideBoundedVector594 =             crRightHandSideBoundedVector10*(crRightHandSideBoundedVector396*crRightHandSideBoundedVector553 + crRightHandSideBoundedVector543*(-1.5000000000000002*crRightHandSideBoundedVector214 + 3.0000000000000009*crRightHandSideBoundedVector218 + crRightHandSideBoundedVector516 - 3.0000000000000009*crRightHandSideBoundedVector87) + crRightHandSideBoundedVector544*crRightHandSideBoundedVector545 + crRightHandSideBoundedVector551*crRightHandSideBoundedVector591 + crRightHandSideBoundedVector554*crRightHandSideBoundedVector592 + crRightHandSideBoundedVector554*crRightHandSideBoundedVector593 - crRightHandSideBoundedVector558*crRightHandSideBoundedVector591);
const double crRightHandSideBoundedVector595 =             1.5000000000000002*crRightHandSideBoundedVector199;
const double crRightHandSideBoundedVector596 =             crRightHandSideBoundedVector224*crRightHandSideBoundedVector574;
const double crRightHandSideBoundedVector597 =             crRightHandSideBoundedVector100*(crRightHandSideBoundedVector396*crRightHandSideBoundedVector571 + crRightHandSideBoundedVector551*crRightHandSideBoundedVector596 - crRightHandSideBoundedVector558*crRightHandSideBoundedVector596 + crRightHandSideBoundedVector561*(-3.0000000000000009*crRightHandSideBoundedVector134 - 1.5000000000000002*crRightHandSideBoundedVector234 + 3.0000000000000009*crRightHandSideBoundedVector238 + crRightHandSideBoundedVector595) + crRightHandSideBoundedVector562*crRightHandSideBoundedVector564 + crRightHandSideBoundedVector572*crRightHandSideBoundedVector592 + crRightHandSideBoundedVector572*crRightHandSideBoundedVector593);
const double crRightHandSideBoundedVector598 =             1.5000000000000002*crRightHandSideBoundedVector207;
const double crRightHandSideBoundedVector599 =             crRightHandSideBoundedVector224*crRightHandSideBoundedVector589;
const double crRightHandSideBoundedVector600 =             crRightHandSideBoundedVector141*(crRightHandSideBoundedVector396*crRightHandSideBoundedVector586 + crRightHandSideBoundedVector551*crRightHandSideBoundedVector599 - crRightHandSideBoundedVector558*crRightHandSideBoundedVector599 + crRightHandSideBoundedVector576*(-3.0000000000000009*crRightHandSideBoundedVector161 - 1.5000000000000002*crRightHandSideBoundedVector246 + 3.0000000000000009*crRightHandSideBoundedVector250 + crRightHandSideBoundedVector598) + crRightHandSideBoundedVector577*crRightHandSideBoundedVector579 + crRightHandSideBoundedVector587*crRightHandSideBoundedVector592 + crRightHandSideBoundedVector587*crRightHandSideBoundedVector593);
const double crRightHandSideBoundedVector601 =             0.1111111111111111*r_ext[1];
const double crRightHandSideBoundedVector602 =             0.1111111111111111*r_ext[2];
const double crRightHandSideBoundedVector603 =             crRightHandSideBoundedVector60*crRightHandSideBoundedVector70;
const double crRightHandSideBoundedVector604 =             crRightHandSideBoundedVector296*crRightHandSideBoundedVector603;
const double crRightHandSideBoundedVector605 =             crRightHandSideBoundedVector423*crRightHandSideBoundedVector60;
const double crRightHandSideBoundedVector606 =             -1.125*U_1_3;
const double crRightHandSideBoundedVector607 =             -1.125*U_2_3;
const double crRightHandSideBoundedVector608 =             4.5000000000000009*crRightHandSideBoundedVector19;
const double crRightHandSideBoundedVector609 =             -4.5*U_0_3 - 6.7500000000000009*crRightHandSideBoundedVector414 + crRightHandSideBoundedVector421*crRightHandSideBoundedVector608 + crRightHandSideBoundedVector606 + crRightHandSideBoundedVector607;
const double crRightHandSideBoundedVector610 =             crRightHandSideBoundedVector289*crRightHandSideBoundedVector609;
const double crRightHandSideBoundedVector611 =             crRightHandSideBoundedVector19*crRightHandSideBoundedVector420;
const double crRightHandSideBoundedVector612 =             crRightHandSideBoundedVector10*(crRightHandSideBoundedVector291*crRightHandSideBoundedVector611 + crRightHandSideBoundedVector422);
const double crRightHandSideBoundedVector613 =             crRightHandSideBoundedVector10*(crRightHandSideBoundedVector422 + crRightHandSideBoundedVector480*crRightHandSideBoundedVector611);
const double crRightHandSideBoundedVector614 =             crRightHandSideBoundedVector127*crRightHandSideBoundedVector447;
const double crRightHandSideBoundedVector615 =             0.027777777777777776*r_ext[0];
const double crRightHandSideBoundedVector616 =             0.027777777777777776*r_ext[2];
const double crRightHandSideBoundedVector617 =             crRightHandSideBoundedVector19*crRightHandSideBoundedVector444;
const double crRightHandSideBoundedVector618 =             crRightHandSideBoundedVector100*(crRightHandSideBoundedVector319*crRightHandSideBoundedVector617 + crRightHandSideBoundedVector446);
const double crRightHandSideBoundedVector619 =             crRightHandSideBoundedVector100*(crRightHandSideBoundedVector446 + crRightHandSideBoundedVector492*crRightHandSideBoundedVector617);
const double crRightHandSideBoundedVector620 =             -1.125*U_0_3;
const double crRightHandSideBoundedVector621 =             crRightHandSideBoundedVector318*(-4.5*U_1_3 - 6.7500000000000009*crRightHandSideBoundedVector440 + crRightHandSideBoundedVector445*crRightHandSideBoundedVector608 + crRightHandSideBoundedVector607 + crRightHandSideBoundedVector620);
const double crRightHandSideBoundedVector622 =             0.16666666666666666*crRightHandSideBoundedVector621;
const double crRightHandSideBoundedVector623 =             crRightHandSideBoundedVector127*crRightHandSideBoundedVector131;
const double crRightHandSideBoundedVector624 =             crRightHandSideBoundedVector326*crRightHandSideBoundedVector623;
const double crRightHandSideBoundedVector625 =             -crRightHandSideBoundedVector134*crRightHandSideBoundedVector622 + crRightHandSideBoundedVector180*crRightHandSideBoundedVector624 - crRightHandSideBoundedVector199*crRightHandSideBoundedVector622 - crRightHandSideBoundedVector322*crRightHandSideBoundedVector405 - crRightHandSideBoundedVector324*crRightHandSideBoundedVector407 + crRightHandSideBoundedVector329*crRightHandSideBoundedVector618 + crRightHandSideBoundedVector497*crRightHandSideBoundedVector619 + crRightHandSideBoundedVector601 + crRightHandSideBoundedVector615 + crRightHandSideBoundedVector616 + crRightHandSideBoundedVector624*crRightHandSideBoundedVector84;
const double crRightHandSideBoundedVector626 =             crRightHandSideBoundedVector155*crRightHandSideBoundedVector467;
const double crRightHandSideBoundedVector627 =             0.027777777777777776*r_ext[1];
const double crRightHandSideBoundedVector628 =             crRightHandSideBoundedVector19*crRightHandSideBoundedVector464;
const double crRightHandSideBoundedVector629 =             crRightHandSideBoundedVector141*(crRightHandSideBoundedVector343*crRightHandSideBoundedVector628 + crRightHandSideBoundedVector466);
const double crRightHandSideBoundedVector630 =             crRightHandSideBoundedVector141*(crRightHandSideBoundedVector466 + crRightHandSideBoundedVector505*crRightHandSideBoundedVector628);
const double crRightHandSideBoundedVector631 =             crRightHandSideBoundedVector342*(-4.5*U_2_3 - 6.7500000000000009*crRightHandSideBoundedVector461 + crRightHandSideBoundedVector465*crRightHandSideBoundedVector608 + crRightHandSideBoundedVector606 + crRightHandSideBoundedVector620);
const double crRightHandSideBoundedVector632 =             0.16666666666666666*crRightHandSideBoundedVector631;
const double crRightHandSideBoundedVector633 =             crRightHandSideBoundedVector155*crRightHandSideBoundedVector158;
const double crRightHandSideBoundedVector634 =             crRightHandSideBoundedVector350*crRightHandSideBoundedVector633;
const double crRightHandSideBoundedVector635 =             -crRightHandSideBoundedVector161*crRightHandSideBoundedVector632 + crRightHandSideBoundedVector180*crRightHandSideBoundedVector634 - crRightHandSideBoundedVector207*crRightHandSideBoundedVector632 + crRightHandSideBoundedVector329*crRightHandSideBoundedVector629 - crRightHandSideBoundedVector346*crRightHandSideBoundedVector405 - crRightHandSideBoundedVector348*crRightHandSideBoundedVector407 + crRightHandSideBoundedVector497*crRightHandSideBoundedVector630 + crRightHandSideBoundedVector602 + crRightHandSideBoundedVector615 + crRightHandSideBoundedVector627 + crRightHandSideBoundedVector634*crRightHandSideBoundedVector84;
const double crRightHandSideBoundedVector636 =             0.66666666666666663*crRightHandSideBoundedVector47;
const double crRightHandSideBoundedVector637 =             4.5000000000000009*crRightHandSideBoundedVector411;
const double crRightHandSideBoundedVector638 =             crRightHandSideBoundedVector60*crRightHandSideBoundedVector76;
const double crRightHandSideBoundedVector639 =             crRightHandSideBoundedVector303*crRightHandSideBoundedVector66;
const double crRightHandSideBoundedVector640 =             0.66666666666666663*crRightHandSideBoundedVector86;
const double crRightHandSideBoundedVector641 =             crRightHandSideBoundedVector56*crRightHandSideBoundedVector70;
const double crRightHandSideBoundedVector642 =             crRightHandSideBoundedVector412*crRightHandSideBoundedVector88;
const double crRightHandSideBoundedVector643 =             -crRightHandSideBoundedVector441;
const double crRightHandSideBoundedVector644 =             crRightHandSideBoundedVector439 + crRightHandSideBoundedVector643;
const double crRightHandSideBoundedVector645 =             0.16666666666666666*crRightHandSideBoundedVector96;
const double crRightHandSideBoundedVector646 =             1.1250000000000002*crRightHandSideBoundedVector437;
const double crRightHandSideBoundedVector647 =             crRightHandSideBoundedVector127*crRightHandSideBoundedVector76;
const double crRightHandSideBoundedVector648 =             -crRightHandSideBoundedVector19*crRightHandSideBoundedVector323 - crRightHandSideBoundedVector19*crRightHandSideBoundedVector325 + crRightHandSideBoundedVector19*crRightHandSideBoundedVector328 + crRightHandSideBoundedVector317*crRightHandSideBoundedVector618 - crRightHandSideBoundedVector332 + crRightHandSideBoundedVector335 + crRightHandSideBoundedVector405*crRightHandSideBoundedVector645 - crRightHandSideBoundedVector646*crRightHandSideBoundedVector647;
const double crRightHandSideBoundedVector649 =             crRightHandSideBoundedVector442 + crRightHandSideBoundedVector643;
const double crRightHandSideBoundedVector650 =             crRightHandSideBoundedVector135*crRightHandSideBoundedVector438;
const double crRightHandSideBoundedVector651 =             crRightHandSideBoundedVector131*crRightHandSideBoundedVector56;
const double crRightHandSideBoundedVector652 =             crRightHandSideBoundedVector322*crRightHandSideBoundedVector84;
const double crRightHandSideBoundedVector653 =             -crRightHandSideBoundedVector19*crRightHandSideBoundedVector494 + crRightHandSideBoundedVector19*crRightHandSideBoundedVector496 + crRightHandSideBoundedVector407*crRightHandSideBoundedVector645 + crRightHandSideBoundedVector491*crRightHandSideBoundedVector619 - crRightHandSideBoundedVector500 + crRightHandSideBoundedVector503 - crRightHandSideBoundedVector646*crRightHandSideBoundedVector651 - crRightHandSideBoundedVector652;
const double crRightHandSideBoundedVector654 =             -crRightHandSideBoundedVector462;
const double crRightHandSideBoundedVector655 =             crRightHandSideBoundedVector460 + crRightHandSideBoundedVector654;
const double crRightHandSideBoundedVector656 =             0.16666666666666666*crRightHandSideBoundedVector139;
const double crRightHandSideBoundedVector657 =             1.1250000000000002*crRightHandSideBoundedVector458;
const double crRightHandSideBoundedVector658 =             crRightHandSideBoundedVector155*crRightHandSideBoundedVector76;
const double crRightHandSideBoundedVector659 =             -crRightHandSideBoundedVector19*crRightHandSideBoundedVector347 - crRightHandSideBoundedVector19*crRightHandSideBoundedVector349 + crRightHandSideBoundedVector19*crRightHandSideBoundedVector352 + crRightHandSideBoundedVector317*crRightHandSideBoundedVector629 - crRightHandSideBoundedVector354 + crRightHandSideBoundedVector356 + crRightHandSideBoundedVector405*crRightHandSideBoundedVector656 - crRightHandSideBoundedVector657*crRightHandSideBoundedVector658;
const double crRightHandSideBoundedVector660 =             crRightHandSideBoundedVector463 + crRightHandSideBoundedVector654;
const double crRightHandSideBoundedVector661 =             crRightHandSideBoundedVector162*crRightHandSideBoundedVector459;
const double crRightHandSideBoundedVector662 =             crRightHandSideBoundedVector158*crRightHandSideBoundedVector56;
const double crRightHandSideBoundedVector663 =             crRightHandSideBoundedVector346*crRightHandSideBoundedVector84;
const double crRightHandSideBoundedVector664 =             -crRightHandSideBoundedVector19*crRightHandSideBoundedVector507 + crRightHandSideBoundedVector19*crRightHandSideBoundedVector509 + crRightHandSideBoundedVector407*crRightHandSideBoundedVector656 + crRightHandSideBoundedVector491*crRightHandSideBoundedVector630 - crRightHandSideBoundedVector511 + crRightHandSideBoundedVector513 - crRightHandSideBoundedVector657*crRightHandSideBoundedVector662 - crRightHandSideBoundedVector663;
const double crRightHandSideBoundedVector665 =             crRightHandSideBoundedVector426*crRightHandSideBoundedVector77;
const double crRightHandSideBoundedVector666 =             crRightHandSideBoundedVector450*crRightHandSideBoundedVector77;
const double crRightHandSideBoundedVector667 =             crRightHandSideBoundedVector470*crRightHandSideBoundedVector77;
const double crRightHandSideBoundedVector668 =             crRightHandSideBoundedVector158*crRightHandSideBoundedVector511 + crRightHandSideBoundedVector158*crRightHandSideBoundedVector663 + crRightHandSideBoundedVector329*crRightHandSideBoundedVector655 + 0.16666666666666666*crRightHandSideBoundedVector451 + 0.16666666666666666*crRightHandSideBoundedVector452 + 0.16666666666666666*crRightHandSideBoundedVector453 - 0.16666666666666666*crRightHandSideBoundedVector454 - 0.16666666666666666*crRightHandSideBoundedVector455 - 0.16666666666666666*crRightHandSideBoundedVector468 - 0.16666666666666666*crRightHandSideBoundedVector469 + crRightHandSideBoundedVector497*crRightHandSideBoundedVector660;
const double crRightHandSideBoundedVector669 =             crRightHandSideBoundedVector131*crRightHandSideBoundedVector500 + crRightHandSideBoundedVector131*crRightHandSideBoundedVector652 + crRightHandSideBoundedVector329*crRightHandSideBoundedVector644 + 0.16666666666666666*crRightHandSideBoundedVector430 + 0.16666666666666666*crRightHandSideBoundedVector431 + 0.16666666666666666*crRightHandSideBoundedVector432 - 0.16666666666666666*crRightHandSideBoundedVector433 - 0.16666666666666666*crRightHandSideBoundedVector434 - 0.16666666666666666*crRightHandSideBoundedVector448 - 0.16666666666666666*crRightHandSideBoundedVector449 + crRightHandSideBoundedVector497*crRightHandSideBoundedVector649;
const double crRightHandSideBoundedVector670 =             0.99999999999999989*crRightHandSideBoundedVector53 + 0.99999999999999989*crRightHandSideBoundedVector54 + 0.99999999999999989*crRightHandSideBoundedVector55 + 0.99999999999999989*crRightHandSideBoundedVector73 + 0.99999999999999989*crRightHandSideBoundedVector74 + 0.99999999999999989*crRightHandSideBoundedVector75;
const double crRightHandSideBoundedVector671 =             DN_DX_1_1*crRightHandSideBoundedVector61;
const double crRightHandSideBoundedVector672 =             DN_DX_1_0*crRightHandSideBoundedVector71;
const double crRightHandSideBoundedVector673 =             0.37500000000000006*crRightHandSideBoundedVector10;
const double crRightHandSideBoundedVector674 =             0.37500000000000006*crRightHandSideBoundedVector411;
const double crRightHandSideBoundedVector675 =             crRightHandSideBoundedVector263*crRightHandSideBoundedVector673 - crRightHandSideBoundedVector265*crRightHandSideBoundedVector674 - crRightHandSideBoundedVector271*crRightHandSideBoundedVector3 + crRightHandSideBoundedVector272*crRightHandSideBoundedVector3;
const double crRightHandSideBoundedVector676 =             DN_DX_1_1*crRightHandSideBoundedVector128;
const double crRightHandSideBoundedVector677 =             0.66666666666666663*crRightHandSideBoundedVector96;
const double crRightHandSideBoundedVector678 =             DN_DX_1_0*crRightHandSideBoundedVector132;
const double crRightHandSideBoundedVector679 =             1.5000000000000002*crRightHandSideBoundedVector100;
const double crRightHandSideBoundedVector680 =             crRightHandSideBoundedVector100*crRightHandSideBoundedVector266;
const double crRightHandSideBoundedVector681 =             DN_DX_1_1*crRightHandSideBoundedVector156;
const double crRightHandSideBoundedVector682 =             DN_DX_1_0*crRightHandSideBoundedVector159;
const double crRightHandSideBoundedVector683 =             0.16666666666666666*crRightHandSideBoundedVector289;
const double crRightHandSideBoundedVector684 =             crRightHandSideBoundedVector60*crRightHandSideBoundedVector673;
const double crRightHandSideBoundedVector685 =             crRightHandSideBoundedVector673*crRightHandSideBoundedVector70;
const double crRightHandSideBoundedVector686 =             1.125*crRightHandSideBoundedVector289;
const double crRightHandSideBoundedVector687 =             crRightHandSideBoundedVector297*crRightHandSideBoundedVector686;
const double crRightHandSideBoundedVector688 =             crRightHandSideBoundedVector685*crRightHandSideBoundedVector84;
const double crRightHandSideBoundedVector689 =             0.1111111111111111*f_ext(0,0);
const double crRightHandSideBoundedVector690 =             crRightHandSideBoundedVector334 + crRightHandSideBoundedVector355 + crRightHandSideBoundedVector689;
const double crRightHandSideBoundedVector691 =             -crRightHandSideBoundedVector293*crRightHandSideBoundedVector683 + crRightHandSideBoundedVector302*crRightHandSideBoundedVector329 - crRightHandSideBoundedVector56*crRightHandSideBoundedVector684 - crRightHandSideBoundedVector66*crRightHandSideBoundedVector685 + crRightHandSideBoundedVector687 + crRightHandSideBoundedVector688 + crRightHandSideBoundedVector690;
const double crRightHandSideBoundedVector692 =             crRightHandSideBoundedVector127*crRightHandSideBoundedVector679;
const double crRightHandSideBoundedVector693 =             crRightHandSideBoundedVector131*crRightHandSideBoundedVector679;
const double crRightHandSideBoundedVector694 =             4.5*crRightHandSideBoundedVector318;
const double crRightHandSideBoundedVector695 =             crRightHandSideBoundedVector327*crRightHandSideBoundedVector694;
const double crRightHandSideBoundedVector696 =             crRightHandSideBoundedVector299*crRightHandSideBoundedVector300;
const double crRightHandSideBoundedVector697 =             crRightHandSideBoundedVector266*crRightHandSideBoundedVector499;
const double crRightHandSideBoundedVector698 =             crRightHandSideBoundedVector697*crRightHandSideBoundedVector83;
const double crRightHandSideBoundedVector699 =             crRightHandSideBoundedVector309 + crRightHandSideBoundedVector689 + 0.44444444444444442*f_ext(1,0);
const double crRightHandSideBoundedVector700 =             DN_DX_1_0*crRightHandSideBoundedVector61;
const double crRightHandSideBoundedVector701 =             DN_DX_1_1*crRightHandSideBoundedVector71;
const double crRightHandSideBoundedVector702 =             crRightHandSideBoundedVector3*crRightHandSideBoundedVector370 - crRightHandSideBoundedVector673*crRightHandSideBoundedVector87;
const double crRightHandSideBoundedVector703 =             crRightHandSideBoundedVector701 + crRightHandSideBoundedVector702;
const double crRightHandSideBoundedVector704 =             crRightHandSideBoundedVector184*crRightHandSideBoundedVector374 - crRightHandSideBoundedVector329*crRightHandSideBoundedVector360;
const double crRightHandSideBoundedVector705 =             crRightHandSideBoundedVector359*crRightHandSideBoundedVector96;
const double crRightHandSideBoundedVector706 =             DN_DX_1_0*crRightHandSideBoundedVector128;
const double crRightHandSideBoundedVector707 =             crRightHandSideBoundedVector100*crRightHandSideBoundedVector300;
const double crRightHandSideBoundedVector708 =             DN_DX_1_1*crRightHandSideBoundedVector132;
const double crRightHandSideBoundedVector709 =             crRightHandSideBoundedVector100*crRightHandSideBoundedVector563;
const double crRightHandSideBoundedVector710 =             crRightHandSideBoundedVector362*crRightHandSideBoundedVector96 + crRightHandSideBoundedVector708 - crRightHandSideBoundedVector709;
const double crRightHandSideBoundedVector711 =             DN_DX_1_0*crRightHandSideBoundedVector156;
const double crRightHandSideBoundedVector712 =             DN_DX_1_1*crRightHandSideBoundedVector159;
const double crRightHandSideBoundedVector713 =             crRightHandSideBoundedVector379 + crRightHandSideBoundedVector712;
const double crRightHandSideBoundedVector714 =             DN_DX_1_0*crRightHandSideBoundedVector262;
const double crRightHandSideBoundedVector715 =             crRightHandSideBoundedVector272*crRightHandSideBoundedVector71 + crRightHandSideBoundedVector373*crRightHandSideBoundedVector61 + 0.16666666666666666*crRightHandSideBoundedVector46 - 0.16666666666666666*crRightHandSideBoundedVector52 - 0.16666666666666666*crRightHandSideBoundedVector62 + crRightHandSideBoundedVector684*crRightHandSideBoundedVector87 - 0.16666666666666666*crRightHandSideBoundedVector72 - 0.99999999999999989*crRightHandSideBoundedVector91;
const double crRightHandSideBoundedVector716 =             crRightHandSideBoundedVector482*crRightHandSideBoundedVector686;
const double crRightHandSideBoundedVector717 =             crRightHandSideBoundedVector60*crRightHandSideBoundedVector674;
const double crRightHandSideBoundedVector718 =             crRightHandSideBoundedVector66*crRightHandSideBoundedVector717;
const double crRightHandSideBoundedVector719 =             0.1111111111111111*f_ext(0,1);
const double crRightHandSideBoundedVector720 =             crRightHandSideBoundedVector502 + crRightHandSideBoundedVector512 + crRightHandSideBoundedVector719;
const double crRightHandSideBoundedVector721 =             -crRightHandSideBoundedVector481*crRightHandSideBoundedVector683 + crRightHandSideBoundedVector485*crRightHandSideBoundedVector497 - crRightHandSideBoundedVector684*crRightHandSideBoundedVector83 - crRightHandSideBoundedVector685*crRightHandSideBoundedVector76 + crRightHandSideBoundedVector716 + crRightHandSideBoundedVector718 + crRightHandSideBoundedVector720;
const double crRightHandSideBoundedVector722 =             crRightHandSideBoundedVector495*crRightHandSideBoundedVector694;
const double crRightHandSideBoundedVector723 =             crRightHandSideBoundedVector300*crRightHandSideBoundedVector484;
const double crRightHandSideBoundedVector724 =             crRightHandSideBoundedVector266*crRightHandSideBoundedVector331;
const double crRightHandSideBoundedVector725 =             crRightHandSideBoundedVector66*crRightHandSideBoundedVector724;
const double crRightHandSideBoundedVector726 =             crRightHandSideBoundedVector489 + crRightHandSideBoundedVector719 + 0.44444444444444442*f_ext(1,1);
const double crRightHandSideBoundedVector727 =             -crRightHandSideBoundedVector183*crRightHandSideBoundedVector673 + crRightHandSideBoundedVector3*crRightHandSideBoundedVector520;
const double crRightHandSideBoundedVector728 =             crRightHandSideBoundedVector700 + crRightHandSideBoundedVector727;
const double crRightHandSideBoundedVector729 =             0.16666666666666666*crRightHandSideBoundedVector301*crRightHandSideBoundedVector87 - crRightHandSideBoundedVector360*crRightHandSideBoundedVector497;
const double crRightHandSideBoundedVector730 =             crRightHandSideBoundedVector263*crRightHandSideBoundedVector674 - crRightHandSideBoundedVector265*crRightHandSideBoundedVector673 + crRightHandSideBoundedVector3*crRightHandSideBoundedVector525 - crRightHandSideBoundedVector3*crRightHandSideBoundedVector526;
const double crRightHandSideBoundedVector731 =             crRightHandSideBoundedVector100*crRightHandSideBoundedVector595;
const double crRightHandSideBoundedVector732 =             crRightHandSideBoundedVector515*crRightHandSideBoundedVector96 + crRightHandSideBoundedVector706 - crRightHandSideBoundedVector731;
const double crRightHandSideBoundedVector733 =             crRightHandSideBoundedVector528 + crRightHandSideBoundedVector711;
const double crRightHandSideBoundedVector734 =             DN_DX_1_1*crRightHandSideBoundedVector262;
const double crRightHandSideBoundedVector735 =             0.16666666666666666*crRightHandSideBoundedVector173 - 0.16666666666666666*crRightHandSideBoundedVector176 - 0.16666666666666666*crRightHandSideBoundedVector177 - 0.16666666666666666*crRightHandSideBoundedVector178 + crRightHandSideBoundedVector183*crRightHandSideBoundedVector685 - 0.99999999999999989*crRightHandSideBoundedVector186 + crRightHandSideBoundedVector523*crRightHandSideBoundedVector71 + crRightHandSideBoundedVector526*crRightHandSideBoundedVector61;
const double crRightHandSideBoundedVector736 =             crRightHandSideBoundedVector423*crRightHandSideBoundedVector70;
const double crRightHandSideBoundedVector737 =             0.1111111111111111*r_ext[0];
const double crRightHandSideBoundedVector738 =             crRightHandSideBoundedVector609*crRightHandSideBoundedVector683;
const double crRightHandSideBoundedVector739 =             crRightHandSideBoundedVector603*crRightHandSideBoundedVector686;
const double crRightHandSideBoundedVector740 =             crRightHandSideBoundedVector180*crRightHandSideBoundedVector739 - crRightHandSideBoundedVector183*crRightHandSideBoundedVector738 + crRightHandSideBoundedVector329*crRightHandSideBoundedVector612 - crRightHandSideBoundedVector405*crRightHandSideBoundedVector684 - crRightHandSideBoundedVector407*crRightHandSideBoundedVector685 + crRightHandSideBoundedVector497*crRightHandSideBoundedVector613 + crRightHandSideBoundedVector616 + crRightHandSideBoundedVector627 + crRightHandSideBoundedVector737 - crRightHandSideBoundedVector738*crRightHandSideBoundedVector87 + crRightHandSideBoundedVector739*crRightHandSideBoundedVector84;
const double crRightHandSideBoundedVector741 =             crRightHandSideBoundedVector623*crRightHandSideBoundedVector694;
const double crRightHandSideBoundedVector742 =             crRightHandSideBoundedVector131*crRightHandSideBoundedVector447;
const double crRightHandSideBoundedVector743 =             crRightHandSideBoundedVector158*crRightHandSideBoundedVector467;
const double crRightHandSideBoundedVector744 =             0.16666666666666666*crRightHandSideBoundedVector3;
const double crRightHandSideBoundedVector745 =             1.1250000000000002*crRightHandSideBoundedVector411;
const double crRightHandSideBoundedVector746 =             crRightHandSideBoundedVector674*crRightHandSideBoundedVector70;
const double crRightHandSideBoundedVector747 =             crRightHandSideBoundedVector19*crRightHandSideBoundedVector687 + crRightHandSideBoundedVector317*crRightHandSideBoundedVector612 + crRightHandSideBoundedVector405*crRightHandSideBoundedVector744 - crRightHandSideBoundedVector56*crRightHandSideBoundedVector717 - crRightHandSideBoundedVector638*crRightHandSideBoundedVector745 - crRightHandSideBoundedVector66*crRightHandSideBoundedVector746 - crRightHandSideBoundedVector688 + crRightHandSideBoundedVector690;
const double crRightHandSideBoundedVector748 =             crRightHandSideBoundedVector684*crRightHandSideBoundedVector84;
const double crRightHandSideBoundedVector749 =             crRightHandSideBoundedVector19*crRightHandSideBoundedVector716 + crRightHandSideBoundedVector407*crRightHandSideBoundedVector744 + crRightHandSideBoundedVector491*crRightHandSideBoundedVector613 - crRightHandSideBoundedVector641*crRightHandSideBoundedVector745 - crRightHandSideBoundedVector718 + crRightHandSideBoundedVector720 - crRightHandSideBoundedVector746*crRightHandSideBoundedVector76 - crRightHandSideBoundedVector748;
const double crRightHandSideBoundedVector750 =             4.5000000000000009*crRightHandSideBoundedVector437;
const double crRightHandSideBoundedVector751 =             crRightHandSideBoundedVector724*crRightHandSideBoundedVector83;
const double crRightHandSideBoundedVector752 =             crRightHandSideBoundedVector329*crRightHandSideBoundedVector541 + 0.16666666666666666*crRightHandSideBoundedVector402 + 0.16666666666666666*crRightHandSideBoundedVector403 + 0.16666666666666666*crRightHandSideBoundedVector404 - 0.16666666666666666*crRightHandSideBoundedVector406 - 0.16666666666666666*crRightHandSideBoundedVector408 - 0.16666666666666666*crRightHandSideBoundedVector424 - 0.16666666666666666*crRightHandSideBoundedVector425 + crRightHandSideBoundedVector497*crRightHandSideBoundedVector542 + crRightHandSideBoundedVector70*crRightHandSideBoundedVector718 + crRightHandSideBoundedVector70*crRightHandSideBoundedVector748;
const double crRightHandSideBoundedVector753 =             DN_DX_2_1*crRightHandSideBoundedVector61;
const double crRightHandSideBoundedVector754 =             DN_DX_2_0*crRightHandSideBoundedVector71;
const double crRightHandSideBoundedVector755 =             DN_DX_2_1*crRightHandSideBoundedVector128;
const double crRightHandSideBoundedVector756 =             DN_DX_2_0*crRightHandSideBoundedVector132;
const double crRightHandSideBoundedVector757 =             DN_DX_2_1*crRightHandSideBoundedVector156;
const double crRightHandSideBoundedVector758 =             0.66666666666666663*crRightHandSideBoundedVector139;
const double crRightHandSideBoundedVector759 =             DN_DX_2_0*crRightHandSideBoundedVector159;
const double crRightHandSideBoundedVector760 =             1.5000000000000002*crRightHandSideBoundedVector141;
const double crRightHandSideBoundedVector761 =             crRightHandSideBoundedVector141*crRightHandSideBoundedVector266;
const double crRightHandSideBoundedVector762 =             crRightHandSideBoundedVector155*crRightHandSideBoundedVector760;
const double crRightHandSideBoundedVector763 =             crRightHandSideBoundedVector158*crRightHandSideBoundedVector760;
const double crRightHandSideBoundedVector764 =             4.5*crRightHandSideBoundedVector342;
const double crRightHandSideBoundedVector765 =             crRightHandSideBoundedVector351*crRightHandSideBoundedVector764;
const double crRightHandSideBoundedVector766 =             crRightHandSideBoundedVector266*crRightHandSideBoundedVector510;
const double crRightHandSideBoundedVector767 =             crRightHandSideBoundedVector766*crRightHandSideBoundedVector83;
const double crRightHandSideBoundedVector768 =             crRightHandSideBoundedVector308 + crRightHandSideBoundedVector689 + 0.44444444444444442*f_ext(2,0);
const double crRightHandSideBoundedVector769 =             DN_DX_2_0*crRightHandSideBoundedVector61;
const double crRightHandSideBoundedVector770 =             DN_DX_2_1*crRightHandSideBoundedVector71;
const double crRightHandSideBoundedVector771 =             crRightHandSideBoundedVector702 + crRightHandSideBoundedVector770;
const double crRightHandSideBoundedVector772 =             DN_DX_2_0*crRightHandSideBoundedVector128;
const double crRightHandSideBoundedVector773 =             DN_DX_2_1*crRightHandSideBoundedVector132;
const double crRightHandSideBoundedVector774 =             crRightHandSideBoundedVector371 + crRightHandSideBoundedVector773;
const double crRightHandSideBoundedVector775 =             crRightHandSideBoundedVector139*crRightHandSideBoundedVector359;
const double crRightHandSideBoundedVector776 =             DN_DX_2_0*crRightHandSideBoundedVector156;
const double crRightHandSideBoundedVector777 =             crRightHandSideBoundedVector141*crRightHandSideBoundedVector300;
const double crRightHandSideBoundedVector778 =             DN_DX_2_1*crRightHandSideBoundedVector159;
const double crRightHandSideBoundedVector779 =             crRightHandSideBoundedVector141*crRightHandSideBoundedVector578;
const double crRightHandSideBoundedVector780 =             crRightHandSideBoundedVector139*crRightHandSideBoundedVector362 + crRightHandSideBoundedVector778 - crRightHandSideBoundedVector779;
const double crRightHandSideBoundedVector781 =             DN_DX_2_0*crRightHandSideBoundedVector262;
const double crRightHandSideBoundedVector782 =             crRightHandSideBoundedVector508*crRightHandSideBoundedVector764;
const double crRightHandSideBoundedVector783 =             crRightHandSideBoundedVector266*crRightHandSideBoundedVector353;
const double crRightHandSideBoundedVector784 =             crRightHandSideBoundedVector66*crRightHandSideBoundedVector783;
const double crRightHandSideBoundedVector785 =             crRightHandSideBoundedVector488 + crRightHandSideBoundedVector719 + 0.44444444444444442*f_ext(2,1);
const double crRightHandSideBoundedVector786 =             crRightHandSideBoundedVector727 + crRightHandSideBoundedVector769;
const double crRightHandSideBoundedVector787 =             crRightHandSideBoundedVector521 + crRightHandSideBoundedVector772;
const double crRightHandSideBoundedVector788 =             crRightHandSideBoundedVector141*crRightHandSideBoundedVector598;
const double crRightHandSideBoundedVector789 =             crRightHandSideBoundedVector139*crRightHandSideBoundedVector515 + crRightHandSideBoundedVector776 - crRightHandSideBoundedVector788;
const double crRightHandSideBoundedVector790 =             DN_DX_2_1*crRightHandSideBoundedVector262;
const double crRightHandSideBoundedVector791 =             crRightHandSideBoundedVector633*crRightHandSideBoundedVector764;
const double crRightHandSideBoundedVector792 =             4.5000000000000009*crRightHandSideBoundedVector458;
const double crRightHandSideBoundedVector793 =             crRightHandSideBoundedVector783*crRightHandSideBoundedVector83;
            rRightHandSideBoundedVector[0]=-1.0*DN_DX_0_0*crRightHandSideBoundedVector137 - 1.0*DN_DX_0_0*crRightHandSideBoundedVector168 - 1.0*DN_DX_0_0*crRightHandSideBoundedVector93 - 1.0*DN_DX_0_1*crRightHandSideBoundedVector188 - 1.0*DN_DX_0_1*crRightHandSideBoundedVector201 - 1.0*DN_DX_0_1*crRightHandSideBoundedVector212 - 1.0*crRightHandSideBoundedVector213;
            rRightHandSideBoundedVector[1]=-DN_DX_0_0*crRightHandSideBoundedVector233 - DN_DX_0_0*crRightHandSideBoundedVector245 - DN_DX_0_0*crRightHandSideBoundedVector257 + crRightHandSideBoundedVector268*(-crRightHandSideBoundedVector258 - crRightHandSideBoundedVector259*crRightHandSideBoundedVector3 + crRightHandSideBoundedVector260*crRightHandSideBoundedVector84 + crRightHandSideBoundedVector261*crRightHandSideBoundedVector262 + crRightHandSideBoundedVector263*crRightHandSideBoundedVector264 - crRightHandSideBoundedVector265*crRightHandSideBoundedVector267) + crRightHandSideBoundedVector279*(crRightHandSideBoundedVector262*crRightHandSideBoundedVector270 - crRightHandSideBoundedVector269 + crRightHandSideBoundedVector278) + crRightHandSideBoundedVector288*(crRightHandSideBoundedVector262*crRightHandSideBoundedVector281 - crRightHandSideBoundedVector280 + crRightHandSideBoundedVector287) - crRightHandSideBoundedVector305*crRightHandSideBoundedVector386 - crRightHandSideBoundedVector313*(DN_DX_0_0*crRightHandSideBoundedVector51 - crRightHandSideBoundedVector290*crRightHandSideBoundedVector293 - crRightHandSideBoundedVector294*crRightHandSideBoundedVector56 - crRightHandSideBoundedVector295*crRightHandSideBoundedVector66 + crRightHandSideBoundedVector298 + crRightHandSideBoundedVector299*crRightHandSideBoundedVector302 + crRightHandSideBoundedVector304 - crRightHandSideBoundedVector307 + crRightHandSideBoundedVector310) - crRightHandSideBoundedVector314*crRightHandSideBoundedVector389 - crRightHandSideBoundedVector338*(DN_DX_0_0*crRightHandSideBoundedVector124 - crRightHandSideBoundedVector316 + crRightHandSideBoundedVector336) - crRightHandSideBoundedVector339*crRightHandSideBoundedVector392 - crRightHandSideBoundedVector358*(DN_DX_0_0*crRightHandSideBoundedVector166 - crRightHandSideBoundedVector341 + crRightHandSideBoundedVector357) + crRightHandSideBoundedVector365*crRightHandSideBoundedVector60 - crRightHandSideBoundedVector367*(crRightHandSideBoundedVector221*crRightHandSideBoundedVector301 - crRightHandSideBoundedVector299*crRightHandSideBoundedVector360 - crRightHandSideBoundedVector359*crRightHandSideBoundedVector361 + crRightHandSideBoundedVector366) - crRightHandSideBoundedVector376*(-crRightHandSideBoundedVector359*crRightHandSideBoundedVector368 + crRightHandSideBoundedVector372 + crRightHandSideBoundedVector375) - crRightHandSideBoundedVector382*(-crRightHandSideBoundedVector359*crRightHandSideBoundedVector377 + crRightHandSideBoundedVector380 + crRightHandSideBoundedVector381) + crRightHandSideBoundedVector393*crRightHandSideBoundedVector71 + crRightHandSideBoundedVector394*crRightHandSideBoundedVector61 - crRightHandSideBoundedVector395*crRightHandSideBoundedVector426 - crRightHandSideBoundedVector395*crRightHandSideBoundedVector450 - crRightHandSideBoundedVector395*crRightHandSideBoundedVector470 + 0.66666666666666663*crRightHandSideBoundedVector46 + crRightHandSideBoundedVector471 + crRightHandSideBoundedVector472 - 0.66666666666666663*crRightHandSideBoundedVector52 - 0.66666666666666663*crRightHandSideBoundedVector62 - 0.66666666666666663*crRightHandSideBoundedVector72 - 1.0*crRightHandSideBoundedVector91;
            rRightHandSideBoundedVector[2]=-DN_DX_0_0*crRightHandSideBoundedVector532 - DN_DX_0_0*crRightHandSideBoundedVector533 - DN_DX_0_0*crRightHandSideBoundedVector534 - DN_DX_0_1*crRightHandSideBoundedVector477 - DN_DX_0_1*crRightHandSideBoundedVector478 - DN_DX_0_1*crRightHandSideBoundedVector479 + 0.66666666666666663*crRightHandSideBoundedVector173 - 0.66666666666666663*crRightHandSideBoundedVector176 - 0.66666666666666663*crRightHandSideBoundedVector177 - 0.66666666666666663*crRightHandSideBoundedVector178 - 1.0*crRightHandSideBoundedVector186 - crRightHandSideBoundedVector268*(crRightHandSideBoundedVector226*crRightHandSideBoundedVector301 - crRightHandSideBoundedVector359*crRightHandSideBoundedVector363 - crRightHandSideBoundedVector360*crRightHandSideBoundedVector484 + crRightHandSideBoundedVector518) - crRightHandSideBoundedVector279*(-crRightHandSideBoundedVector359*crRightHandSideBoundedVector369 + crRightHandSideBoundedVector522 + crRightHandSideBoundedVector524) - crRightHandSideBoundedVector288*(-crRightHandSideBoundedVector359*crRightHandSideBoundedVector378 + crRightHandSideBoundedVector529 + crRightHandSideBoundedVector530) - crRightHandSideBoundedVector313*(-DN_DX_0_0*crRightHandSideBoundedVector409 + DN_DX_0_1*crRightHandSideBoundedVector175 - crRightHandSideBoundedVector290*crRightHandSideBoundedVector481 - crRightHandSideBoundedVector294*crRightHandSideBoundedVector83 - crRightHandSideBoundedVector295*crRightHandSideBoundedVector76 + crRightHandSideBoundedVector483 + crRightHandSideBoundedVector484*crRightHandSideBoundedVector485 + crRightHandSideBoundedVector487 + crRightHandSideBoundedVector490) - crRightHandSideBoundedVector338*(-DN_DX_0_0*crRightHandSideBoundedVector435 + DN_DX_0_1*crRightHandSideBoundedVector194 + crRightHandSideBoundedVector504) - crRightHandSideBoundedVector358*(-DN_DX_0_0*crRightHandSideBoundedVector456 + DN_DX_0_1*crRightHandSideBoundedVector210 + crRightHandSideBoundedVector514) - crRightHandSideBoundedVector367*(-crRightHandSideBoundedVector180*crRightHandSideBoundedVector260 - crRightHandSideBoundedVector258*crRightHandSideBoundedVector262 + crRightHandSideBoundedVector261 + crRightHandSideBoundedVector263*crRightHandSideBoundedVector267 - crRightHandSideBoundedVector264*crRightHandSideBoundedVector265 + crRightHandSideBoundedVector3*crRightHandSideBoundedVector519) - crRightHandSideBoundedVector376*(-crRightHandSideBoundedVector262*crRightHandSideBoundedVector269 + crRightHandSideBoundedVector270 + crRightHandSideBoundedVector527) - crRightHandSideBoundedVector382*(-crRightHandSideBoundedVector262*crRightHandSideBoundedVector280 + crRightHandSideBoundedVector281 + crRightHandSideBoundedVector531) - crRightHandSideBoundedVector426*crRightHandSideBoundedVector537 - crRightHandSideBoundedVector450*crRightHandSideBoundedVector537 - crRightHandSideBoundedVector470*crRightHandSideBoundedVector537 + crRightHandSideBoundedVector517*crRightHandSideBoundedVector70 + crRightHandSideBoundedVector535*crRightHandSideBoundedVector61 + crRightHandSideBoundedVector536*crRightHandSideBoundedVector71 + crRightHandSideBoundedVector538 + crRightHandSideBoundedVector539;
            rRightHandSideBoundedVector[3]=-DN_DX_0_0*crRightHandSideBoundedVector560 - DN_DX_0_0*crRightHandSideBoundedVector575 - DN_DX_0_0*crRightHandSideBoundedVector590 - DN_DX_0_1*crRightHandSideBoundedVector594 - DN_DX_0_1*crRightHandSideBoundedVector597 - DN_DX_0_1*crRightHandSideBoundedVector600 - crRightHandSideBoundedVector221*crRightHandSideBoundedVector423 - crRightHandSideBoundedVector226*crRightHandSideBoundedVector423 - crRightHandSideBoundedVector268*(-DN_DX_0_0*crRightHandSideBoundedVector642 - DN_DX_0_1*crRightHandSideBoundedVector542 + crRightHandSideBoundedVector19*crRightHandSideBoundedVector483 + crRightHandSideBoundedVector260*crRightHandSideBoundedVector407 - crRightHandSideBoundedVector303*crRightHandSideBoundedVector76 - crRightHandSideBoundedVector486*crRightHandSideBoundedVector83 - crRightHandSideBoundedVector487 + crRightHandSideBoundedVector490 + crRightHandSideBoundedVector613*crRightHandSideBoundedVector640 - crRightHandSideBoundedVector637*crRightHandSideBoundedVector641) - crRightHandSideBoundedVector279*(-DN_DX_0_0*crRightHandSideBoundedVector650 - DN_DX_0_1*crRightHandSideBoundedVector649 + crRightHandSideBoundedVector653) - crRightHandSideBoundedVector288*(-DN_DX_0_0*crRightHandSideBoundedVector661 - DN_DX_0_1*crRightHandSideBoundedVector660 + crRightHandSideBoundedVector664) + crRightHandSideBoundedVector299*crRightHandSideBoundedVector541 + crRightHandSideBoundedVector304*crRightHandSideBoundedVector60 - crRightHandSideBoundedVector313*(DN_DX_0_0*crRightHandSideBoundedVector605 + crRightHandSideBoundedVector180*crRightHandSideBoundedVector604 - crRightHandSideBoundedVector221*crRightHandSideBoundedVector610 - crRightHandSideBoundedVector226*crRightHandSideBoundedVector610 - crRightHandSideBoundedVector294*crRightHandSideBoundedVector405 - crRightHandSideBoundedVector295*crRightHandSideBoundedVector407 + crRightHandSideBoundedVector299*crRightHandSideBoundedVector612 + crRightHandSideBoundedVector306*crRightHandSideBoundedVector422 + crRightHandSideBoundedVector484*crRightHandSideBoundedVector613 + crRightHandSideBoundedVector601 + crRightHandSideBoundedVector602 + crRightHandSideBoundedVector604*crRightHandSideBoundedVector84 + 0.44444444444444442*r_ext[0]) - crRightHandSideBoundedVector338*(DN_DX_0_0*crRightHandSideBoundedVector614 + crRightHandSideBoundedVector315*crRightHandSideBoundedVector446 + crRightHandSideBoundedVector625) - crRightHandSideBoundedVector358*(DN_DX_0_0*crRightHandSideBoundedVector626 + crRightHandSideBoundedVector340*crRightHandSideBoundedVector466 + crRightHandSideBoundedVector635) - crRightHandSideBoundedVector367*(-DN_DX_0_0*crRightHandSideBoundedVector541 + crRightHandSideBoundedVector19*crRightHandSideBoundedVector298 - crRightHandSideBoundedVector19*crRightHandSideBoundedVector307 + crRightHandSideBoundedVector260*crRightHandSideBoundedVector405 - crRightHandSideBoundedVector304 + crRightHandSideBoundedVector310 - crRightHandSideBoundedVector486*crRightHandSideBoundedVector56 + crRightHandSideBoundedVector612*crRightHandSideBoundedVector636 - crRightHandSideBoundedVector637*crRightHandSideBoundedVector638 - crRightHandSideBoundedVector639) - crRightHandSideBoundedVector376*(-DN_DX_0_0*crRightHandSideBoundedVector644 - crRightHandSideBoundedVector19*crRightHandSideBoundedVector316 + crRightHandSideBoundedVector648) - crRightHandSideBoundedVector382*(-DN_DX_0_0*crRightHandSideBoundedVector655 - crRightHandSideBoundedVector19*crRightHandSideBoundedVector341 + crRightHandSideBoundedVector659) + 0.66666666666666663*crRightHandSideBoundedVector402 + 0.66666666666666663*crRightHandSideBoundedVector403 + 0.66666666666666663*crRightHandSideBoundedVector404 - 0.66666666666666663*crRightHandSideBoundedVector406 - 0.66666666666666663*crRightHandSideBoundedVector408 + crRightHandSideBoundedVector484*crRightHandSideBoundedVector542 + crRightHandSideBoundedVector60*crRightHandSideBoundedVector639 - crRightHandSideBoundedVector665*(crRightHandSideBoundedVector366 + crRightHandSideBoundedVector518) - crRightHandSideBoundedVector666*(crRightHandSideBoundedVector372 + crRightHandSideBoundedVector522) - crRightHandSideBoundedVector667*(crRightHandSideBoundedVector380 + crRightHandSideBoundedVector529) + crRightHandSideBoundedVector668 + crRightHandSideBoundedVector669;
            rRightHandSideBoundedVector[4]=-DN_DX_1_0*crRightHandSideBoundedVector367 - DN_DX_1_0*crRightHandSideBoundedVector376 - DN_DX_1_0*crRightHandSideBoundedVector382 - DN_DX_1_1*crRightHandSideBoundedVector268 - DN_DX_1_1*crRightHandSideBoundedVector279 - DN_DX_1_1*crRightHandSideBoundedVector288 - crRightHandSideBoundedVector670;
            rRightHandSideBoundedVector[5]=-DN_DX_1_0*crRightHandSideBoundedVector233 - DN_DX_1_0*crRightHandSideBoundedVector245 - DN_DX_1_0*crRightHandSideBoundedVector257 - DN_DX_1_1*crRightHandSideBoundedVector532 - DN_DX_1_1*crRightHandSideBoundedVector533 - DN_DX_1_1*crRightHandSideBoundedVector534 + 0.66666666666666663*crRightHandSideBoundedVector121 - 0.66666666666666663*crRightHandSideBoundedVector125 + crRightHandSideBoundedVector127*crRightHandSideBoundedVector709 + crRightHandSideBoundedVector128*crRightHandSideBoundedVector394 - 0.66666666666666663*crRightHandSideBoundedVector129 + crRightHandSideBoundedVector132*crRightHandSideBoundedVector393 - 0.66666666666666663*crRightHandSideBoundedVector133 + crRightHandSideBoundedVector268*(crRightHandSideBoundedVector262*crRightHandSideBoundedVector672 - crRightHandSideBoundedVector671 + crRightHandSideBoundedVector675) + crRightHandSideBoundedVector279*(-crRightHandSideBoundedVector259*crRightHandSideBoundedVector96 + crRightHandSideBoundedVector262*crRightHandSideBoundedVector678 + crRightHandSideBoundedVector273*crRightHandSideBoundedVector679 - crRightHandSideBoundedVector276*crRightHandSideBoundedVector680 - crRightHandSideBoundedVector676 + crRightHandSideBoundedVector677*crRightHandSideBoundedVector84) + crRightHandSideBoundedVector288*(crRightHandSideBoundedVector262*crRightHandSideBoundedVector682 + crRightHandSideBoundedVector287 - crRightHandSideBoundedVector681) - crRightHandSideBoundedVector313*(DN_DX_1_0*crRightHandSideBoundedVector51 - DN_DX_1_1*crRightHandSideBoundedVector409 + crRightHandSideBoundedVector691) - crRightHandSideBoundedVector338*(DN_DX_1_0*crRightHandSideBoundedVector124 - DN_DX_1_1*crRightHandSideBoundedVector435 - crRightHandSideBoundedVector321*crRightHandSideBoundedVector636 + crRightHandSideBoundedVector331*crRightHandSideBoundedVector696 - crRightHandSideBoundedVector56*crRightHandSideBoundedVector692 - crRightHandSideBoundedVector66*crRightHandSideBoundedVector693 + crRightHandSideBoundedVector695 + crRightHandSideBoundedVector698 + crRightHandSideBoundedVector699) - crRightHandSideBoundedVector358*(DN_DX_1_0*crRightHandSideBoundedVector166 - DN_DX_1_1*crRightHandSideBoundedVector456 + crRightHandSideBoundedVector357) - crRightHandSideBoundedVector367*(-crRightHandSideBoundedVector359*crRightHandSideBoundedVector700 + crRightHandSideBoundedVector703 + crRightHandSideBoundedVector704) - crRightHandSideBoundedVector376*(crRightHandSideBoundedVector241*crRightHandSideBoundedVector707 - crRightHandSideBoundedVector299*crRightHandSideBoundedVector705 - crRightHandSideBoundedVector359*crRightHandSideBoundedVector706 + crRightHandSideBoundedVector710) - crRightHandSideBoundedVector382*(-crRightHandSideBoundedVector359*crRightHandSideBoundedVector711 + crRightHandSideBoundedVector381 + crRightHandSideBoundedVector713) - crRightHandSideBoundedVector426*crRightHandSideBoundedVector714 - crRightHandSideBoundedVector450*crRightHandSideBoundedVector714 - crRightHandSideBoundedVector470*crRightHandSideBoundedVector714 + crRightHandSideBoundedVector471 + crRightHandSideBoundedVector715;
            rRightHandSideBoundedVector[6]=-DN_DX_1_0*crRightHandSideBoundedVector532 - DN_DX_1_0*crRightHandSideBoundedVector533 - DN_DX_1_0*crRightHandSideBoundedVector534 - DN_DX_1_1*crRightHandSideBoundedVector477 - DN_DX_1_1*crRightHandSideBoundedVector478 - DN_DX_1_1*crRightHandSideBoundedVector479 + crRightHandSideBoundedVector128*crRightHandSideBoundedVector535 + crRightHandSideBoundedVector131*crRightHandSideBoundedVector731 + crRightHandSideBoundedVector132*crRightHandSideBoundedVector536 + 0.66666666666666663*crRightHandSideBoundedVector192 - 0.66666666666666663*crRightHandSideBoundedVector195 - 0.66666666666666663*crRightHandSideBoundedVector196 - 0.66666666666666663*crRightHandSideBoundedVector197 - crRightHandSideBoundedVector268*(-crRightHandSideBoundedVector359*crRightHandSideBoundedVector701 + crRightHandSideBoundedVector728 + crRightHandSideBoundedVector729) - crRightHandSideBoundedVector279*(crRightHandSideBoundedVector243*crRightHandSideBoundedVector707 - crRightHandSideBoundedVector359*crRightHandSideBoundedVector708 - crRightHandSideBoundedVector484*crRightHandSideBoundedVector705 + crRightHandSideBoundedVector732) - crRightHandSideBoundedVector288*(-crRightHandSideBoundedVector359*crRightHandSideBoundedVector712 + crRightHandSideBoundedVector530 + crRightHandSideBoundedVector733) - crRightHandSideBoundedVector313*(-DN_DX_1_0*crRightHandSideBoundedVector409 + DN_DX_1_1*crRightHandSideBoundedVector175 + crRightHandSideBoundedVector721) - crRightHandSideBoundedVector338*(-DN_DX_1_0*crRightHandSideBoundedVector435 + DN_DX_1_1*crRightHandSideBoundedVector194 - crRightHandSideBoundedVector493*crRightHandSideBoundedVector640 + crRightHandSideBoundedVector499*crRightHandSideBoundedVector723 - crRightHandSideBoundedVector692*crRightHandSideBoundedVector83 - crRightHandSideBoundedVector693*crRightHandSideBoundedVector76 + crRightHandSideBoundedVector722 + crRightHandSideBoundedVector725 + crRightHandSideBoundedVector726) - crRightHandSideBoundedVector358*(-DN_DX_1_0*crRightHandSideBoundedVector456 + DN_DX_1_1*crRightHandSideBoundedVector210 + crRightHandSideBoundedVector514) - crRightHandSideBoundedVector367*(-crRightHandSideBoundedVector262*crRightHandSideBoundedVector671 + crRightHandSideBoundedVector672 + crRightHandSideBoundedVector730) - crRightHandSideBoundedVector376*(-crRightHandSideBoundedVector180*crRightHandSideBoundedVector677 - crRightHandSideBoundedVector262*crRightHandSideBoundedVector676 + crRightHandSideBoundedVector273*crRightHandSideBoundedVector680 - crRightHandSideBoundedVector276*crRightHandSideBoundedVector679 + crRightHandSideBoundedVector519*crRightHandSideBoundedVector96 + crRightHandSideBoundedVector678) - crRightHandSideBoundedVector382*(-crRightHandSideBoundedVector262*crRightHandSideBoundedVector681 + crRightHandSideBoundedVector531 + crRightHandSideBoundedVector682) - crRightHandSideBoundedVector426*crRightHandSideBoundedVector734 - crRightHandSideBoundedVector450*crRightHandSideBoundedVector734 - crRightHandSideBoundedVector470*crRightHandSideBoundedVector734 + crRightHandSideBoundedVector538 + crRightHandSideBoundedVector735;
            rRightHandSideBoundedVector[7]=-DN_DX_1_0*crRightHandSideBoundedVector560 - DN_DX_1_0*crRightHandSideBoundedVector575 - DN_DX_1_0*crRightHandSideBoundedVector590 - DN_DX_1_1*crRightHandSideBoundedVector594 - DN_DX_1_1*crRightHandSideBoundedVector597 - DN_DX_1_1*crRightHandSideBoundedVector600 + crRightHandSideBoundedVector131*crRightHandSideBoundedVector725 + crRightHandSideBoundedVector131*crRightHandSideBoundedVector751 - crRightHandSideBoundedVector241*crRightHandSideBoundedVector447 - crRightHandSideBoundedVector243*crRightHandSideBoundedVector447 - crRightHandSideBoundedVector268*(-DN_DX_1_0*crRightHandSideBoundedVector642 - DN_DX_1_1*crRightHandSideBoundedVector542 + crRightHandSideBoundedVector749) - crRightHandSideBoundedVector279*(-DN_DX_1_0*crRightHandSideBoundedVector650 - DN_DX_1_1*crRightHandSideBoundedVector649 + crRightHandSideBoundedVector19*crRightHandSideBoundedVector722 + crRightHandSideBoundedVector407*crRightHandSideBoundedVector677 + crRightHandSideBoundedVector619*crRightHandSideBoundedVector640 - crRightHandSideBoundedVector651*crRightHandSideBoundedVector750 - crRightHandSideBoundedVector697*crRightHandSideBoundedVector76 - crRightHandSideBoundedVector725 + crRightHandSideBoundedVector726 - crRightHandSideBoundedVector751) - crRightHandSideBoundedVector288*(-DN_DX_1_0*crRightHandSideBoundedVector661 - DN_DX_1_1*crRightHandSideBoundedVector660 + crRightHandSideBoundedVector664) + crRightHandSideBoundedVector299*crRightHandSideBoundedVector644 - crRightHandSideBoundedVector313*(DN_DX_1_0*crRightHandSideBoundedVector605 + DN_DX_1_1*crRightHandSideBoundedVector736 + crRightHandSideBoundedVector740) - crRightHandSideBoundedVector338*(DN_DX_1_0*crRightHandSideBoundedVector614 + DN_DX_1_1*crRightHandSideBoundedVector742 + crRightHandSideBoundedVector180*crRightHandSideBoundedVector741 - crRightHandSideBoundedVector241*crRightHandSideBoundedVector621 - crRightHandSideBoundedVector243*crRightHandSideBoundedVector621 + crRightHandSideBoundedVector299*crRightHandSideBoundedVector618 - crRightHandSideBoundedVector405*crRightHandSideBoundedVector692 - crRightHandSideBoundedVector407*crRightHandSideBoundedVector693 + crRightHandSideBoundedVector484*crRightHandSideBoundedVector619 + crRightHandSideBoundedVector602 + crRightHandSideBoundedVector737 + crRightHandSideBoundedVector741*crRightHandSideBoundedVector84 + 0.44444444444444442*r_ext[1]) - crRightHandSideBoundedVector358*(DN_DX_1_0*crRightHandSideBoundedVector626 + DN_DX_1_1*crRightHandSideBoundedVector743 + crRightHandSideBoundedVector635) - crRightHandSideBoundedVector367*(-DN_DX_1_0*crRightHandSideBoundedVector541 - DN_DX_1_1*crRightHandSideBoundedVector642 + crRightHandSideBoundedVector747) - crRightHandSideBoundedVector376*(-DN_DX_1_0*crRightHandSideBoundedVector644 - DN_DX_1_1*crRightHandSideBoundedVector650 + crRightHandSideBoundedVector19*crRightHandSideBoundedVector695 + crRightHandSideBoundedVector405*crRightHandSideBoundedVector677 - crRightHandSideBoundedVector56*crRightHandSideBoundedVector724 + crRightHandSideBoundedVector618*crRightHandSideBoundedVector636 - crRightHandSideBoundedVector647*crRightHandSideBoundedVector750 - crRightHandSideBoundedVector66*crRightHandSideBoundedVector697 - crRightHandSideBoundedVector698 + crRightHandSideBoundedVector699) - crRightHandSideBoundedVector382*(-DN_DX_1_0*crRightHandSideBoundedVector655 - DN_DX_1_1*crRightHandSideBoundedVector661 + crRightHandSideBoundedVector659) + 0.66666666666666663*crRightHandSideBoundedVector430 + 0.66666666666666663*crRightHandSideBoundedVector431 + 0.66666666666666663*crRightHandSideBoundedVector432 - 0.66666666666666663*crRightHandSideBoundedVector433 - 0.66666666666666663*crRightHandSideBoundedVector434 + crRightHandSideBoundedVector484*crRightHandSideBoundedVector649 - crRightHandSideBoundedVector665*(crRightHandSideBoundedVector703 + crRightHandSideBoundedVector728) - crRightHandSideBoundedVector666*(crRightHandSideBoundedVector710 + crRightHandSideBoundedVector732) - crRightHandSideBoundedVector667*(crRightHandSideBoundedVector713 + crRightHandSideBoundedVector733) + crRightHandSideBoundedVector668 + crRightHandSideBoundedVector752;
            rRightHandSideBoundedVector[8]=-DN_DX_2_0*crRightHandSideBoundedVector367 - DN_DX_2_0*crRightHandSideBoundedVector376 - DN_DX_2_0*crRightHandSideBoundedVector382 - DN_DX_2_1*crRightHandSideBoundedVector268 - DN_DX_2_1*crRightHandSideBoundedVector279 - DN_DX_2_1*crRightHandSideBoundedVector288 - crRightHandSideBoundedVector670;
            rRightHandSideBoundedVector[9]=-DN_DX_2_0*crRightHandSideBoundedVector233 - DN_DX_2_0*crRightHandSideBoundedVector245 - DN_DX_2_0*crRightHandSideBoundedVector257 - DN_DX_2_1*crRightHandSideBoundedVector532 - DN_DX_2_1*crRightHandSideBoundedVector533 - DN_DX_2_1*crRightHandSideBoundedVector534 + 0.66666666666666663*crRightHandSideBoundedVector154 + crRightHandSideBoundedVector155*crRightHandSideBoundedVector779 + crRightHandSideBoundedVector156*crRightHandSideBoundedVector394 - 0.66666666666666663*crRightHandSideBoundedVector157 + crRightHandSideBoundedVector159*crRightHandSideBoundedVector393 - 0.66666666666666663*crRightHandSideBoundedVector160 - 0.66666666666666663*crRightHandSideBoundedVector167 + crRightHandSideBoundedVector268*(crRightHandSideBoundedVector262*crRightHandSideBoundedVector754 + crRightHandSideBoundedVector675 - crRightHandSideBoundedVector753) + crRightHandSideBoundedVector279*(crRightHandSideBoundedVector262*crRightHandSideBoundedVector756 + crRightHandSideBoundedVector278 - crRightHandSideBoundedVector755) + crRightHandSideBoundedVector288*(-crRightHandSideBoundedVector139*crRightHandSideBoundedVector259 + crRightHandSideBoundedVector262*crRightHandSideBoundedVector759 + crRightHandSideBoundedVector282*crRightHandSideBoundedVector760 - crRightHandSideBoundedVector285*crRightHandSideBoundedVector761 - crRightHandSideBoundedVector757 + crRightHandSideBoundedVector758*crRightHandSideBoundedVector84) - crRightHandSideBoundedVector313*(DN_DX_2_0*crRightHandSideBoundedVector51 - DN_DX_2_1*crRightHandSideBoundedVector409 + crRightHandSideBoundedVector691) - crRightHandSideBoundedVector338*(DN_DX_2_0*crRightHandSideBoundedVector124 - DN_DX_2_1*crRightHandSideBoundedVector435 + crRightHandSideBoundedVector336) - crRightHandSideBoundedVector358*(DN_DX_2_0*crRightHandSideBoundedVector166 - DN_DX_2_1*crRightHandSideBoundedVector456 - crRightHandSideBoundedVector345*crRightHandSideBoundedVector636 + crRightHandSideBoundedVector353*crRightHandSideBoundedVector696 - crRightHandSideBoundedVector56*crRightHandSideBoundedVector762 - crRightHandSideBoundedVector66*crRightHandSideBoundedVector763 + crRightHandSideBoundedVector765 + crRightHandSideBoundedVector767 + crRightHandSideBoundedVector768) - crRightHandSideBoundedVector367*(-crRightHandSideBoundedVector359*crRightHandSideBoundedVector769 + crRightHandSideBoundedVector704 + crRightHandSideBoundedVector771) - crRightHandSideBoundedVector376*(-crRightHandSideBoundedVector359*crRightHandSideBoundedVector772 + crRightHandSideBoundedVector375 + crRightHandSideBoundedVector774) - crRightHandSideBoundedVector382*(crRightHandSideBoundedVector253*crRightHandSideBoundedVector777 - crRightHandSideBoundedVector299*crRightHandSideBoundedVector775 - crRightHandSideBoundedVector359*crRightHandSideBoundedVector776 + crRightHandSideBoundedVector780) - crRightHandSideBoundedVector426*crRightHandSideBoundedVector781 - crRightHandSideBoundedVector450*crRightHandSideBoundedVector781 - crRightHandSideBoundedVector470*crRightHandSideBoundedVector781 + crRightHandSideBoundedVector472 + crRightHandSideBoundedVector715;
            rRightHandSideBoundedVector[10]=-DN_DX_2_0*crRightHandSideBoundedVector532 - DN_DX_2_0*crRightHandSideBoundedVector533 - DN_DX_2_0*crRightHandSideBoundedVector534 - DN_DX_2_1*crRightHandSideBoundedVector477 - DN_DX_2_1*crRightHandSideBoundedVector478 - DN_DX_2_1*crRightHandSideBoundedVector479 + crRightHandSideBoundedVector156*crRightHandSideBoundedVector535 + crRightHandSideBoundedVector158*crRightHandSideBoundedVector788 + crRightHandSideBoundedVector159*crRightHandSideBoundedVector536 + 0.66666666666666663*crRightHandSideBoundedVector203 - 0.66666666666666663*crRightHandSideBoundedVector204 - 0.66666666666666663*crRightHandSideBoundedVector205 - 0.66666666666666663*crRightHandSideBoundedVector211 - crRightHandSideBoundedVector268*(-crRightHandSideBoundedVector359*crRightHandSideBoundedVector770 + crRightHandSideBoundedVector729 + crRightHandSideBoundedVector786) - crRightHandSideBoundedVector279*(-crRightHandSideBoundedVector359*crRightHandSideBoundedVector773 + crRightHandSideBoundedVector524 + crRightHandSideBoundedVector787) - crRightHandSideBoundedVector288*(crRightHandSideBoundedVector255*crRightHandSideBoundedVector777 - crRightHandSideBoundedVector359*crRightHandSideBoundedVector778 - crRightHandSideBoundedVector484*crRightHandSideBoundedVector775 + crRightHandSideBoundedVector789) - crRightHandSideBoundedVector313*(-DN_DX_2_0*crRightHandSideBoundedVector409 + DN_DX_2_1*crRightHandSideBoundedVector175 + crRightHandSideBoundedVector721) - crRightHandSideBoundedVector338*(-DN_DX_2_0*crRightHandSideBoundedVector435 + DN_DX_2_1*crRightHandSideBoundedVector194 + crRightHandSideBoundedVector504) - crRightHandSideBoundedVector358*(-DN_DX_2_0*crRightHandSideBoundedVector456 + DN_DX_2_1*crRightHandSideBoundedVector210 - crRightHandSideBoundedVector506*crRightHandSideBoundedVector640 + crRightHandSideBoundedVector510*crRightHandSideBoundedVector723 - crRightHandSideBoundedVector76*crRightHandSideBoundedVector763 - crRightHandSideBoundedVector762*crRightHandSideBoundedVector83 + crRightHandSideBoundedVector782 + crRightHandSideBoundedVector784 + crRightHandSideBoundedVector785) - crRightHandSideBoundedVector367*(-crRightHandSideBoundedVector262*crRightHandSideBoundedVector753 + crRightHandSideBoundedVector730 + crRightHandSideBoundedVector754) - crRightHandSideBoundedVector376*(-crRightHandSideBoundedVector262*crRightHandSideBoundedVector755 + crRightHandSideBoundedVector527 + crRightHandSideBoundedVector756) - crRightHandSideBoundedVector382*(crRightHandSideBoundedVector139*crRightHandSideBoundedVector519 - crRightHandSideBoundedVector180*crRightHandSideBoundedVector758 - crRightHandSideBoundedVector262*crRightHandSideBoundedVector757 + crRightHandSideBoundedVector282*crRightHandSideBoundedVector761 - crRightHandSideBoundedVector285*crRightHandSideBoundedVector760 + crRightHandSideBoundedVector759) - crRightHandSideBoundedVector426*crRightHandSideBoundedVector790 - crRightHandSideBoundedVector450*crRightHandSideBoundedVector790 - crRightHandSideBoundedVector470*crRightHandSideBoundedVector790 + crRightHandSideBoundedVector539 + crRightHandSideBoundedVector735;
            rRightHandSideBoundedVector[11]=-DN_DX_2_0*crRightHandSideBoundedVector560 - DN_DX_2_0*crRightHandSideBoundedVector575 - DN_DX_2_0*crRightHandSideBoundedVector590 - DN_DX_2_1*crRightHandSideBoundedVector594 - DN_DX_2_1*crRightHandSideBoundedVector597 - DN_DX_2_1*crRightHandSideBoundedVector600 + crRightHandSideBoundedVector158*crRightHandSideBoundedVector784 + crRightHandSideBoundedVector158*crRightHandSideBoundedVector793 - crRightHandSideBoundedVector253*crRightHandSideBoundedVector467 - crRightHandSideBoundedVector255*crRightHandSideBoundedVector467 - crRightHandSideBoundedVector268*(-DN_DX_2_0*crRightHandSideBoundedVector642 - DN_DX_2_1*crRightHandSideBoundedVector542 + crRightHandSideBoundedVector749) - crRightHandSideBoundedVector279*(-DN_DX_2_0*crRightHandSideBoundedVector650 - DN_DX_2_1*crRightHandSideBoundedVector649 + crRightHandSideBoundedVector653) - crRightHandSideBoundedVector288*(-DN_DX_2_0*crRightHandSideBoundedVector661 - DN_DX_2_1*crRightHandSideBoundedVector660 + crRightHandSideBoundedVector19*crRightHandSideBoundedVector782 + crRightHandSideBoundedVector407*crRightHandSideBoundedVector758 + crRightHandSideBoundedVector630*crRightHandSideBoundedVector640 - crRightHandSideBoundedVector662*crRightHandSideBoundedVector792 - crRightHandSideBoundedVector76*crRightHandSideBoundedVector766 - crRightHandSideBoundedVector784 + crRightHandSideBoundedVector785 - crRightHandSideBoundedVector793) + crRightHandSideBoundedVector299*crRightHandSideBoundedVector655 - crRightHandSideBoundedVector313*(DN_DX_2_0*crRightHandSideBoundedVector605 + DN_DX_2_1*crRightHandSideBoundedVector736 + crRightHandSideBoundedVector740) - crRightHandSideBoundedVector338*(DN_DX_2_0*crRightHandSideBoundedVector614 + DN_DX_2_1*crRightHandSideBoundedVector742 + crRightHandSideBoundedVector625) - crRightHandSideBoundedVector358*(DN_DX_2_0*crRightHandSideBoundedVector626 + DN_DX_2_1*crRightHandSideBoundedVector743 + crRightHandSideBoundedVector180*crRightHandSideBoundedVector791 - crRightHandSideBoundedVector253*crRightHandSideBoundedVector631 - crRightHandSideBoundedVector255*crRightHandSideBoundedVector631 + crRightHandSideBoundedVector299*crRightHandSideBoundedVector629 - crRightHandSideBoundedVector405*crRightHandSideBoundedVector762 - crRightHandSideBoundedVector407*crRightHandSideBoundedVector763 + crRightHandSideBoundedVector484*crRightHandSideBoundedVector630 + crRightHandSideBoundedVector601 + crRightHandSideBoundedVector737 + crRightHandSideBoundedVector791*crRightHandSideBoundedVector84 + 0.44444444444444442*r_ext[2]) - crRightHandSideBoundedVector367*(-DN_DX_2_0*crRightHandSideBoundedVector541 - DN_DX_2_1*crRightHandSideBoundedVector642 + crRightHandSideBoundedVector747) - crRightHandSideBoundedVector376*(-DN_DX_2_0*crRightHandSideBoundedVector644 - DN_DX_2_1*crRightHandSideBoundedVector650 + crRightHandSideBoundedVector648) - crRightHandSideBoundedVector382*(-DN_DX_2_0*crRightHandSideBoundedVector655 - DN_DX_2_1*crRightHandSideBoundedVector661 + crRightHandSideBoundedVector19*crRightHandSideBoundedVector765 + crRightHandSideBoundedVector405*crRightHandSideBoundedVector758 - crRightHandSideBoundedVector56*crRightHandSideBoundedVector783 + crRightHandSideBoundedVector629*crRightHandSideBoundedVector636 - crRightHandSideBoundedVector658*crRightHandSideBoundedVector792 - crRightHandSideBoundedVector66*crRightHandSideBoundedVector766 - crRightHandSideBoundedVector767 + crRightHandSideBoundedVector768) + 0.66666666666666663*crRightHandSideBoundedVector451 + 0.66666666666666663*crRightHandSideBoundedVector452 + 0.66666666666666663*crRightHandSideBoundedVector453 - 0.66666666666666663*crRightHandSideBoundedVector454 - 0.66666666666666663*crRightHandSideBoundedVector455 + crRightHandSideBoundedVector484*crRightHandSideBoundedVector660 - crRightHandSideBoundedVector665*(crRightHandSideBoundedVector771 + crRightHandSideBoundedVector786) - crRightHandSideBoundedVector666*(crRightHandSideBoundedVector774 + crRightHandSideBoundedVector787) - crRightHandSideBoundedVector667*(crRightHandSideBoundedVector780 + crRightHandSideBoundedVector789) + crRightHandSideBoundedVector669 + crRightHandSideBoundedVector752;

    }

    // Here we assume that all the weights of the gauss points are the same so we multiply at the end by Volume/n_nodes
    rRightHandSideBoundedVector *= data.volume / static_cast<double>(n_nodes);

    KRATOS_CATCH("")
}

template<>
void CompressibleNavierStokesExplicit<3>::CalculateRightHandSideInternal(
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
    const array_1d<double, n_nodes> &r_ext = data.r_ext;
    const BoundedMatrix<double, n_nodes, 3> &f_ext = data.f_ext;
    const double mu = data.mu;
    const double lambda = data.lambda;
    const double c_v = data.c_v;
    const double gamma = data.gamma;

    // Stabilization parameters
    const double stab_c1 = 12.0;
    const double stab_c2 = 2.0;
    const double stab_c3 = 1.0;

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

    const double &dUdt_0_0 = data.dUdt(0, 0);
    const double &dUdt_0_1 = data.dUdt(0, 1);
    const double &dUdt_0_2 = data.dUdt(0, 2);
    const double &dUdt_0_3 = data.dUdt(0, 3);
    const double &dUdt_0_4 = data.dUdt(0, 4);
    const double &dUdt_1_0 = data.dUdt(1, 0);
    const double &dUdt_1_1 = data.dUdt(1, 1);
    const double &dUdt_1_2 = data.dUdt(1, 2);
    const double &dUdt_1_3 = data.dUdt(1, 3);
    const double &dUdt_1_4 = data.dUdt(1, 4);
    const double &dUdt_2_0 = data.dUdt(2, 0);
    const double &dUdt_2_1 = data.dUdt(2, 1);
    const double &dUdt_2_2 = data.dUdt(2, 2);
    const double &dUdt_2_3 = data.dUdt(2, 3);
    const double &dUdt_2_4 = data.dUdt(2, 4);
    const double &dUdt_3_0 = data.dUdt(3, 0);
    const double &dUdt_3_1 = data.dUdt(3, 1);
    const double &dUdt_3_2 = data.dUdt(3, 2);
    const double &dUdt_3_3 = data.dUdt(3, 3);
    const double &dUdt_3_4 = data.dUdt(3, 4);

    // Shock capturing parameters
    double nu_sc = 0.0;
    double nu_st = 0.0;
    double k_sc = 0.0;
    double k_st = 0.0;
    double lin_m_norm = 1.0; // This is intentionally set to a non-zero number to avoid dividing by zero
    array_1d<double, 3> lin_m = ZeroVector(3);
    if (data.ShockCapturing) {
        nu_sc = data.nu_sc;
        k_sc = data.lambda_sc;

        const double rho_avg = (U_0_0 + U_1_0 + U_2_0 + U_3_0) / 4.0;
        array_1d<double, 3> v_avg;
        v_avg[0] = (U_0_1 / U_0_0 + U_1_1 / U_1_0 + U_2_1 / U_2_0 + U_3_1 / U_3_0 ) / 4.0;
        v_avg[1] = (U_0_2 / U_0_0 + U_1_2 / U_1_0 + U_2_2 / U_2_0 + U_3_2 / U_3_0 ) / 4.0;
        v_avg[2] = (U_0_3 / U_0_0 + U_1_3 / U_1_0 + U_2_3 / U_2_0 + U_3_3 / U_3_0 ) / 4.0;
        const double v_avg_norm = norm_2(v_avg);
        const double tot_ener_avg = (U_0_4 + U_1_4 + U_2_4 + U_3_4) / 4.0;
        const double c_avg = gamma * (gamma - 1.0) * ((tot_ener_avg / rho_avg) - 0.5 * std::pow(v_avg_norm, 2));

        // Source terms midpoint values
        double r_ext_avg = 0.0;
        array_1d<double, 3> f_ext_avg = ZeroVector(3);
        for (unsigned int i_node = 0; i_node < n_nodes; ++i_node) {
            f_ext_avg[0] += f_ext(i_node, 0);
            f_ext_avg[1] += f_ext(i_node, 1);
            f_ext_avg[2] += f_ext(i_node, 2);
            r_ext_avg += r_ext[i_node];
        }
        f_ext_avg /= n_nodes;
        r_ext_avg /= n_nodes;
        const double f_ext_avg_norm = norm_2(f_ext_avg);

        const double alpha = lambda / (rho_avg * gamma * c_v);
        const double aux_1 = std::pow(r_ext_avg, 2) + 2.0 * std::pow(c_avg, 2) * std::pow(f_ext_avg_norm, 2) + std::sqrt(std::pow(r_ext_avg, 4) + 4.0 * std::pow(c_avg, 2) * std::pow(f_ext_avg_norm, 2) * std::pow(r_ext_avg, 2));
        const double aux_2 = 2.0 * std::pow(c_avg, 4);
        const double tau_rho = (stab_c2 * (v_avg_norm + c_avg)/ h) + (stab_c3 * std::sqrt(aux_1 / aux_2));
        const double tau_m_avg = 1.0 / ((4.0 * stab_c1 * mu / 3.0 / rho_avg / std::pow(h, 2)) + tau_rho);
        const double tau_et_avg = 1.0 / ((stab_c1 * alpha / std::pow(h, 2)) + tau_rho);

        nu_st = std::max(0.0, nu_sc - tau_m_avg * std::pow(v_avg_norm, 2));
        k_st = std::max(0.0, k_sc - tau_et_avg * std::pow(v_avg_norm, 2));

        lin_m[0] = (U_0_1 + U_1_1 + U_2_1 + U_3_1) / 4.0;
        lin_m[1] = (U_0_2 + U_1_2 + U_2_2 + U_3_2) / 4.0;
        lin_m[2] = (U_0_3 + U_1_3 + U_2_3 + U_3_3) / 4.0;
        const double zero_tol = 1.0e-12;
        const double aux = norm_2(lin_m);
        lin_m_norm = aux > zero_tol ? aux : 1.0;
    }

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

    if (data.UseOSS) {
        // Projections container accesses
        const double &ResProj_0_0 = data.ResProj(0, 0);
        const double &ResProj_0_1 = data.ResProj(0, 1);
        const double &ResProj_0_2 = data.ResProj(0, 2);
        const double &ResProj_0_3 = data.ResProj(0, 3);
        const double &ResProj_0_4 = data.ResProj(0, 4);
        const double &ResProj_1_0 = data.ResProj(1, 0);
        const double &ResProj_1_1 = data.ResProj(1, 1);
        const double &ResProj_1_2 = data.ResProj(1, 2);
        const double &ResProj_1_3 = data.ResProj(1, 3);
        const double &ResProj_1_4 = data.ResProj(1, 4);
        const double &ResProj_2_0 = data.ResProj(2, 0);
        const double &ResProj_2_1 = data.ResProj(2, 1);
        const double &ResProj_2_2 = data.ResProj(2, 2);
        const double &ResProj_2_3 = data.ResProj(2, 3);
        const double &ResProj_2_4 = data.ResProj(2, 4);
        const double &ResProj_3_0 = data.ResProj(3, 0);
        const double &ResProj_3_1 = data.ResProj(3, 1);
        const double &ResProj_3_2 = data.ResProj(3, 2);
        const double &ResProj_3_3 = data.ResProj(3, 3);
        const double &ResProj_3_4 = data.ResProj(3, 4);

        //substitute_rhs_3D_OSS
    } else {
        //substitute_rhs_3D_ASGS
    }

    // Here we assume that all the weights of the gauss points are the same so we multiply at the end by Volume/n_nodes
    rRightHandSideBoundedVector *= data.volume / static_cast<double>(n_nodes);

    KRATOS_CATCH("")
}

template <>
void CompressibleNavierStokesExplicit<2>::AddExplicitContribution(const ProcessInfo &rCurrentProcessInfo)
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
void CompressibleNavierStokesExplicit<3>::AddExplicitContribution(const ProcessInfo &rCurrentProcessInfo)
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
void CompressibleNavierStokesExplicit<2>::CalculateMassMatrix(
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
void CompressibleNavierStokesExplicit<3>::CalculateMassMatrix(
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

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation

template class CompressibleNavierStokesExplicit<2>;
template class CompressibleNavierStokesExplicit<3>;

}
