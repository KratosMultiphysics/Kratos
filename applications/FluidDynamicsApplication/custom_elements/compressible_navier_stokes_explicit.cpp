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
        KRATOS_ERROR_IF_NOT(this->GetGeometry()[i].SolutionStepsDataHas(EXTERNAL_PRESSURE)) << "Missing EXTERNAL_PRESSURE variable on solution step data for node " << this->GetGeometry()[i].Id();

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
            rData.r(i) = r_node.FastGetSolutionStepValue(EXTERNAL_PRESSURE);
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
            rData.r(i) = r_node.FastGetSolutionStepValue(EXTERNAL_PRESSURE);
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
    const array_1d<double, n_nodes> &r = data.r;
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
    const array_1d<double, n_nodes> &r = data.r;
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
    const array_1d<double, n_nodes> &r = data.r;
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
const double ctot_ener_proj28 =             0.16666666666666666*r[1];
const double ctot_ener_proj29 =             0.16666666666666666*r[2];
const double ctot_ener_proj30 =             ctot_ener_proj12*(ctot_ener_proj28 + ctot_ener_proj29 + 0.66666666666666663*r[0]);
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
const double ctot_ener_proj79 =             0.16666666666666666*r[0];
const double ctot_ener_proj80 =             ctot_ener_proj67*(ctot_ener_proj28 + ctot_ener_proj79 + 0.66666666666666663*r[2]);
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
const double ctot_ener_proj116 =             ctot_ener_proj108*(ctot_ener_proj29 + ctot_ener_proj79 + 0.66666666666666663*r[1]);
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
    const array_1d<double, n_nodes> &r = data.r;
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
    const array_1d<double, n_nodes> &r = data.r;
    const BoundedMatrix<double, n_nodes, 2> &f_ext = data.f_ext;
    const double mu = data.mu;
    const double lambda = data.lambda;
    const double c_v = data.c_v;
    const double gamma = data.gamma;

    // Stabilization parameters
    const double stab_c1 = 12.0;
    const double stab_c2 = 2.0;

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

        const double alpha = lambda / (rho_avg * gamma * c_v);
        const double tau_m_avg = 1.0 / ((4.0 * stab_c1 * mu / 3.0 / rho_avg / std::pow(h, 2)) + (stab_c2 * (v_avg_norm + c_avg)/ h));
        const double tau_et_avg = 1.0 / ((stab_c1 * alpha / std::pow(h, 2)) + (stab_c2 * (v_avg_norm + c_avg)/ h));

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

        const double crRightHandSideBoundedVector0 =             DN_DX_0_0*h;
const double crRightHandSideBoundedVector1 =             0.16666666666666666*U_1_0;
const double crRightHandSideBoundedVector2 =             0.16666666666666666*U_2_0;
const double crRightHandSideBoundedVector3 =             0.66666666666666663*U_0_0 + crRightHandSideBoundedVector1 + crRightHandSideBoundedVector2;
const double crRightHandSideBoundedVector4 =             1.0/crRightHandSideBoundedVector3;
const double crRightHandSideBoundedVector5 =             stab_c1/h;
const double crRightHandSideBoundedVector6 =             crRightHandSideBoundedVector4*crRightHandSideBoundedVector5;
const double crRightHandSideBoundedVector7 =             1.3333333333333333*mu;
const double crRightHandSideBoundedVector8 =             0.25*U_1_0;
const double crRightHandSideBoundedVector9 =             0.25*U_2_0;
const double crRightHandSideBoundedVector10 =             U_0_0 + crRightHandSideBoundedVector8 + crRightHandSideBoundedVector9;
const double crRightHandSideBoundedVector11 =             pow(crRightHandSideBoundedVector10, -2);
const double crRightHandSideBoundedVector12 =             0.25*U_1_1;
const double crRightHandSideBoundedVector13 =             0.25*U_2_1;
const double crRightHandSideBoundedVector14 =             pow(U_0_1 + crRightHandSideBoundedVector12 + crRightHandSideBoundedVector13, 2);
const double crRightHandSideBoundedVector15 =             0.25*U_1_2;
const double crRightHandSideBoundedVector16 =             0.25*U_2_2;
const double crRightHandSideBoundedVector17 =             pow(U_0_2 + crRightHandSideBoundedVector15 + crRightHandSideBoundedVector16, 2);
const double crRightHandSideBoundedVector18 =             crRightHandSideBoundedVector14 + crRightHandSideBoundedVector17;
const double crRightHandSideBoundedVector19 =             sqrt(gamma);
const double crRightHandSideBoundedVector20 =             gamma - 1;
const double crRightHandSideBoundedVector21 =             0.5*crRightHandSideBoundedVector11;
const double crRightHandSideBoundedVector22 =             0.16666666666666666*U_1_3;
const double crRightHandSideBoundedVector23 =             0.16666666666666666*U_2_3;
const double crRightHandSideBoundedVector24 =             0.66666666666666663*U_0_3 + crRightHandSideBoundedVector22 + crRightHandSideBoundedVector23;
const double crRightHandSideBoundedVector25 =             crRightHandSideBoundedVector19*sqrt(-crRightHandSideBoundedVector20*(crRightHandSideBoundedVector14*crRightHandSideBoundedVector21 + crRightHandSideBoundedVector17*crRightHandSideBoundedVector21 - crRightHandSideBoundedVector24*crRightHandSideBoundedVector4)) + 1.0*sqrt(crRightHandSideBoundedVector11*crRightHandSideBoundedVector18);
const double crRightHandSideBoundedVector26 =             crRightHandSideBoundedVector25*stab_c2;
const double crRightHandSideBoundedVector27 =             1.0/(crRightHandSideBoundedVector26 + crRightHandSideBoundedVector6*crRightHandSideBoundedVector7);
const double crRightHandSideBoundedVector28 =             0.16666666666666666*f_ext(1,0);
const double crRightHandSideBoundedVector29 =             0.16666666666666666*f_ext(2,0);
const double crRightHandSideBoundedVector30 =             crRightHandSideBoundedVector28 + crRightHandSideBoundedVector29 + 0.66666666666666663*f_ext(0,0);
const double crRightHandSideBoundedVector31 =             crRightHandSideBoundedVector3*crRightHandSideBoundedVector30;
const double crRightHandSideBoundedVector32 =             DN_DX_0_0*U_0_0 + DN_DX_1_0*U_1_0 + DN_DX_2_0*U_2_0;
const double crRightHandSideBoundedVector33 =             1.0000000000000002*crRightHandSideBoundedVector14;
const double crRightHandSideBoundedVector34 =             0.50000000000000011*crRightHandSideBoundedVector20;
const double crRightHandSideBoundedVector35 =             crRightHandSideBoundedVector18*crRightHandSideBoundedVector34;
const double crRightHandSideBoundedVector36 =             crRightHandSideBoundedVector11*(-crRightHandSideBoundedVector33 + crRightHandSideBoundedVector35);
const double crRightHandSideBoundedVector37 =             crRightHandSideBoundedVector32*crRightHandSideBoundedVector36;
const double crRightHandSideBoundedVector38 =             DN_DX_0_1*U_0_2;
const double crRightHandSideBoundedVector39 =             DN_DX_1_1*U_1_2;
const double crRightHandSideBoundedVector40 =             DN_DX_2_1*U_2_2;
const double crRightHandSideBoundedVector41 =             crRightHandSideBoundedVector38 + crRightHandSideBoundedVector39 + crRightHandSideBoundedVector40;
const double crRightHandSideBoundedVector42 =             0.66666666666666663*U_0_1;
const double crRightHandSideBoundedVector43 =             0.16666666666666666*U_1_1;
const double crRightHandSideBoundedVector44 =             0.16666666666666666*U_2_1;
const double crRightHandSideBoundedVector45 =             crRightHandSideBoundedVector42 + crRightHandSideBoundedVector43 + crRightHandSideBoundedVector44;
const double crRightHandSideBoundedVector46 =             crRightHandSideBoundedVector4*crRightHandSideBoundedVector45;
const double crRightHandSideBoundedVector47 =             crRightHandSideBoundedVector41*crRightHandSideBoundedVector46;
const double crRightHandSideBoundedVector48 =             DN_DX_0_1*U_0_1;
const double crRightHandSideBoundedVector49 =             DN_DX_1_1*U_1_1;
const double crRightHandSideBoundedVector50 =             DN_DX_2_1*U_2_1;
const double crRightHandSideBoundedVector51 =             crRightHandSideBoundedVector48 + crRightHandSideBoundedVector49 + crRightHandSideBoundedVector50;
const double crRightHandSideBoundedVector52 =             0.66666666666666663*U_0_2;
const double crRightHandSideBoundedVector53 =             0.16666666666666666*U_1_2;
const double crRightHandSideBoundedVector54 =             0.16666666666666666*U_2_2;
const double crRightHandSideBoundedVector55 =             crRightHandSideBoundedVector52 + crRightHandSideBoundedVector53 + crRightHandSideBoundedVector54;
const double crRightHandSideBoundedVector56 =             crRightHandSideBoundedVector4*crRightHandSideBoundedVector55;
const double crRightHandSideBoundedVector57 =             crRightHandSideBoundedVector51*crRightHandSideBoundedVector56;
const double crRightHandSideBoundedVector58 =             DN_DX_0_0*U_0_1;
const double crRightHandSideBoundedVector59 =             DN_DX_1_0*U_1_1;
const double crRightHandSideBoundedVector60 =             DN_DX_2_0*U_2_1;
const double crRightHandSideBoundedVector61 =             crRightHandSideBoundedVector58 + crRightHandSideBoundedVector59 + crRightHandSideBoundedVector60;
const double crRightHandSideBoundedVector62 =             1.0*gamma;
const double crRightHandSideBoundedVector63 =             3.0 - crRightHandSideBoundedVector62;
const double crRightHandSideBoundedVector64 =             crRightHandSideBoundedVector61*crRightHandSideBoundedVector63;
const double crRightHandSideBoundedVector65 =             DN_DX_0_0*U_0_2;
const double crRightHandSideBoundedVector66 =             DN_DX_1_0*U_1_2;
const double crRightHandSideBoundedVector67 =             DN_DX_2_0*U_2_2;
const double crRightHandSideBoundedVector68 =             crRightHandSideBoundedVector65 + crRightHandSideBoundedVector66 + crRightHandSideBoundedVector67;
const double crRightHandSideBoundedVector69 =             crRightHandSideBoundedVector20*crRightHandSideBoundedVector68;
const double crRightHandSideBoundedVector70 =             1.0*crRightHandSideBoundedVector69;
const double crRightHandSideBoundedVector71 =             DN_DX_0_1*U_0_0 + DN_DX_1_1*U_1_0 + DN_DX_2_1*U_2_0;
const double crRightHandSideBoundedVector72 =             crRightHandSideBoundedVector55*crRightHandSideBoundedVector71;
const double crRightHandSideBoundedVector73 =             2.2500000000000004*crRightHandSideBoundedVector45;
const double crRightHandSideBoundedVector74 =             crRightHandSideBoundedVector11*crRightHandSideBoundedVector73;
const double crRightHandSideBoundedVector75 =             DN_DX_0_0*U_0_3 + DN_DX_1_0*U_1_3 + DN_DX_2_0*U_2_3;
const double crRightHandSideBoundedVector76 =             crRightHandSideBoundedVector20*crRightHandSideBoundedVector75;
const double crRightHandSideBoundedVector77 =             0.16666666666666666*ResProj_2_1 + crRightHandSideBoundedVector76 + 0.16666666666666666*dUdt_2_1;
const double crRightHandSideBoundedVector78 =             0.16666666666666666*ResProj_1_1 + 0.16666666666666666*dUdt_1_1;
const double crRightHandSideBoundedVector79 =             crRightHandSideBoundedVector27*(0.66666666666666663*ResProj_0_1 - crRightHandSideBoundedVector31 + crRightHandSideBoundedVector37 + crRightHandSideBoundedVector46*crRightHandSideBoundedVector64 + crRightHandSideBoundedVector47 - crRightHandSideBoundedVector56*crRightHandSideBoundedVector70 + crRightHandSideBoundedVector57 - crRightHandSideBoundedVector72*crRightHandSideBoundedVector74 + crRightHandSideBoundedVector77 + crRightHandSideBoundedVector78 + 0.66666666666666663*dUdt_0_1);
const double crRightHandSideBoundedVector80 =             0.16666666666666666*U_0_0;
const double crRightHandSideBoundedVector81 =             0.66666666666666663*U_1_0 + crRightHandSideBoundedVector2 + crRightHandSideBoundedVector80;
const double crRightHandSideBoundedVector82 =             1.0/crRightHandSideBoundedVector81;
const double crRightHandSideBoundedVector83 =             crRightHandSideBoundedVector5*crRightHandSideBoundedVector7;
const double crRightHandSideBoundedVector84 =             0.25*U_0_0;
const double crRightHandSideBoundedVector85 =             U_1_0 + crRightHandSideBoundedVector84 + crRightHandSideBoundedVector9;
const double crRightHandSideBoundedVector86 =             pow(crRightHandSideBoundedVector85, -2);
const double crRightHandSideBoundedVector87 =             0.25*U_0_1;
const double crRightHandSideBoundedVector88 =             pow(U_1_1 + crRightHandSideBoundedVector13 + crRightHandSideBoundedVector87, 2);
const double crRightHandSideBoundedVector89 =             0.25*U_0_2;
const double crRightHandSideBoundedVector90 =             pow(U_1_2 + crRightHandSideBoundedVector16 + crRightHandSideBoundedVector89, 2);
const double crRightHandSideBoundedVector91 =             crRightHandSideBoundedVector88 + crRightHandSideBoundedVector90;
const double crRightHandSideBoundedVector92 =             0.5*crRightHandSideBoundedVector86;
const double crRightHandSideBoundedVector93 =             0.16666666666666666*U_0_3;
const double crRightHandSideBoundedVector94 =             0.66666666666666663*U_1_3 + crRightHandSideBoundedVector23 + crRightHandSideBoundedVector93;
const double crRightHandSideBoundedVector95 =             crRightHandSideBoundedVector19*sqrt(-crRightHandSideBoundedVector20*(-crRightHandSideBoundedVector82*crRightHandSideBoundedVector94 + crRightHandSideBoundedVector88*crRightHandSideBoundedVector92 + crRightHandSideBoundedVector90*crRightHandSideBoundedVector92)) + 1.0*sqrt(crRightHandSideBoundedVector86*crRightHandSideBoundedVector91);
const double crRightHandSideBoundedVector96 =             crRightHandSideBoundedVector95*stab_c2;
const double crRightHandSideBoundedVector97 =             1.0/(crRightHandSideBoundedVector82*crRightHandSideBoundedVector83 + crRightHandSideBoundedVector96);
const double crRightHandSideBoundedVector98 =             0.16666666666666666*f_ext(0,0);
const double crRightHandSideBoundedVector99 =             crRightHandSideBoundedVector29 + crRightHandSideBoundedVector98 + 0.66666666666666663*f_ext(1,0);
const double crRightHandSideBoundedVector100 =             crRightHandSideBoundedVector81*crRightHandSideBoundedVector99;
const double crRightHandSideBoundedVector101 =             1.0000000000000002*crRightHandSideBoundedVector88;
const double crRightHandSideBoundedVector102 =             crRightHandSideBoundedVector34*crRightHandSideBoundedVector91;
const double crRightHandSideBoundedVector103 =             crRightHandSideBoundedVector86*(-crRightHandSideBoundedVector101 + crRightHandSideBoundedVector102);
const double crRightHandSideBoundedVector104 =             crRightHandSideBoundedVector103*crRightHandSideBoundedVector32;
const double crRightHandSideBoundedVector105 =             0.16666666666666666*U_0_1;
const double crRightHandSideBoundedVector106 =             0.66666666666666663*U_1_1 + crRightHandSideBoundedVector105 + crRightHandSideBoundedVector44;
const double crRightHandSideBoundedVector107 =             crRightHandSideBoundedVector106*crRightHandSideBoundedVector82;
const double crRightHandSideBoundedVector108 =             crRightHandSideBoundedVector107*crRightHandSideBoundedVector41;
const double crRightHandSideBoundedVector109 =             0.16666666666666666*U_0_2;
const double crRightHandSideBoundedVector110 =             0.66666666666666663*U_1_2 + crRightHandSideBoundedVector109 + crRightHandSideBoundedVector54;
const double crRightHandSideBoundedVector111 =             crRightHandSideBoundedVector110*crRightHandSideBoundedVector82;
const double crRightHandSideBoundedVector112 =             crRightHandSideBoundedVector111*crRightHandSideBoundedVector51;
const double crRightHandSideBoundedVector113 =             2.2500000000000004*crRightHandSideBoundedVector106;
const double crRightHandSideBoundedVector114 =             crRightHandSideBoundedVector110*crRightHandSideBoundedVector71;
const double crRightHandSideBoundedVector115 =             crRightHandSideBoundedVector114*crRightHandSideBoundedVector86;
const double crRightHandSideBoundedVector116 =             0.16666666666666666*ResProj_0_1 + 0.16666666666666666*dUdt_0_1;
const double crRightHandSideBoundedVector117 =             crRightHandSideBoundedVector97*(0.66666666666666663*ResProj_1_1 - crRightHandSideBoundedVector100 + crRightHandSideBoundedVector104 + crRightHandSideBoundedVector107*crRightHandSideBoundedVector64 + crRightHandSideBoundedVector108 - crRightHandSideBoundedVector111*crRightHandSideBoundedVector70 + crRightHandSideBoundedVector112 - crRightHandSideBoundedVector113*crRightHandSideBoundedVector115 + crRightHandSideBoundedVector116 + crRightHandSideBoundedVector77 + 0.66666666666666663*dUdt_1_1);
const double crRightHandSideBoundedVector118 =             0.66666666666666663*U_2_0 + crRightHandSideBoundedVector1 + crRightHandSideBoundedVector80;
const double crRightHandSideBoundedVector119 =             1.0/crRightHandSideBoundedVector118;
const double crRightHandSideBoundedVector120 =             U_2_0 + crRightHandSideBoundedVector8 + crRightHandSideBoundedVector84;
const double crRightHandSideBoundedVector121 =             pow(crRightHandSideBoundedVector120, -2);
const double crRightHandSideBoundedVector122 =             pow(U_2_1 + crRightHandSideBoundedVector12 + crRightHandSideBoundedVector87, 2);
const double crRightHandSideBoundedVector123 =             pow(U_2_2 + crRightHandSideBoundedVector15 + crRightHandSideBoundedVector89, 2);
const double crRightHandSideBoundedVector124 =             crRightHandSideBoundedVector122 + crRightHandSideBoundedVector123;
const double crRightHandSideBoundedVector125 =             0.5*crRightHandSideBoundedVector121;
const double crRightHandSideBoundedVector126 =             0.66666666666666663*U_2_3 + crRightHandSideBoundedVector22 + crRightHandSideBoundedVector93;
const double crRightHandSideBoundedVector127 =             crRightHandSideBoundedVector19*sqrt(-crRightHandSideBoundedVector20*(-crRightHandSideBoundedVector119*crRightHandSideBoundedVector126 + crRightHandSideBoundedVector122*crRightHandSideBoundedVector125 + crRightHandSideBoundedVector123*crRightHandSideBoundedVector125)) + 1.0*sqrt(crRightHandSideBoundedVector121*crRightHandSideBoundedVector124);
const double crRightHandSideBoundedVector128 =             crRightHandSideBoundedVector127*stab_c2;
const double crRightHandSideBoundedVector129 =             1.0/(crRightHandSideBoundedVector119*crRightHandSideBoundedVector83 + crRightHandSideBoundedVector128);
const double crRightHandSideBoundedVector130 =             crRightHandSideBoundedVector28 + crRightHandSideBoundedVector98 + 0.66666666666666663*f_ext(2,0);
const double crRightHandSideBoundedVector131 =             crRightHandSideBoundedVector118*crRightHandSideBoundedVector130;
const double crRightHandSideBoundedVector132 =             1.0000000000000002*crRightHandSideBoundedVector122;
const double crRightHandSideBoundedVector133 =             crRightHandSideBoundedVector124*crRightHandSideBoundedVector34;
const double crRightHandSideBoundedVector134 =             crRightHandSideBoundedVector121*(-crRightHandSideBoundedVector132 + crRightHandSideBoundedVector133);
const double crRightHandSideBoundedVector135 =             crRightHandSideBoundedVector134*crRightHandSideBoundedVector32;
const double crRightHandSideBoundedVector136 =             0.66666666666666663*U_2_1 + crRightHandSideBoundedVector105 + crRightHandSideBoundedVector43;
const double crRightHandSideBoundedVector137 =             crRightHandSideBoundedVector119*crRightHandSideBoundedVector136;
const double crRightHandSideBoundedVector138 =             crRightHandSideBoundedVector137*crRightHandSideBoundedVector41;
const double crRightHandSideBoundedVector139 =             0.66666666666666663*U_2_2 + crRightHandSideBoundedVector109 + crRightHandSideBoundedVector53;
const double crRightHandSideBoundedVector140 =             crRightHandSideBoundedVector119*crRightHandSideBoundedVector139;
const double crRightHandSideBoundedVector141 =             crRightHandSideBoundedVector140*crRightHandSideBoundedVector51;
const double crRightHandSideBoundedVector142 =             2.2500000000000004*crRightHandSideBoundedVector136;
const double crRightHandSideBoundedVector143 =             crRightHandSideBoundedVector139*crRightHandSideBoundedVector71;
const double crRightHandSideBoundedVector144 =             crRightHandSideBoundedVector121*crRightHandSideBoundedVector143;
const double crRightHandSideBoundedVector145 =             crRightHandSideBoundedVector129*(0.66666666666666663*ResProj_2_1 + crRightHandSideBoundedVector116 - crRightHandSideBoundedVector131 + crRightHandSideBoundedVector135 + crRightHandSideBoundedVector137*crRightHandSideBoundedVector64 + crRightHandSideBoundedVector138 - crRightHandSideBoundedVector140*crRightHandSideBoundedVector70 + crRightHandSideBoundedVector141 - crRightHandSideBoundedVector142*crRightHandSideBoundedVector144 + crRightHandSideBoundedVector76 + crRightHandSideBoundedVector78 + 0.66666666666666663*dUdt_2_1);
const double crRightHandSideBoundedVector146 =             DN_DX_0_1*h;
const double crRightHandSideBoundedVector147 =             0.16666666666666666*f_ext(1,1);
const double crRightHandSideBoundedVector148 =             0.16666666666666666*f_ext(2,1);
const double crRightHandSideBoundedVector149 =             crRightHandSideBoundedVector147 + crRightHandSideBoundedVector148 + 0.66666666666666663*f_ext(0,1);
const double crRightHandSideBoundedVector150 =             crRightHandSideBoundedVector149*crRightHandSideBoundedVector3;
const double crRightHandSideBoundedVector151 =             1.0000000000000002*crRightHandSideBoundedVector17;
const double crRightHandSideBoundedVector152 =             crRightHandSideBoundedVector11*(-crRightHandSideBoundedVector151 + crRightHandSideBoundedVector35);
const double crRightHandSideBoundedVector153 =             crRightHandSideBoundedVector152*crRightHandSideBoundedVector71;
const double crRightHandSideBoundedVector154 =             crRightHandSideBoundedVector46*crRightHandSideBoundedVector68;
const double crRightHandSideBoundedVector155 =             crRightHandSideBoundedVector56*crRightHandSideBoundedVector61;
const double crRightHandSideBoundedVector156 =             crRightHandSideBoundedVector41*crRightHandSideBoundedVector63;
const double crRightHandSideBoundedVector157 =             crRightHandSideBoundedVector20*crRightHandSideBoundedVector51;
const double crRightHandSideBoundedVector158 =             1.0*crRightHandSideBoundedVector157;
const double crRightHandSideBoundedVector159 =             2.2500000000000004*crRightHandSideBoundedVector55;
const double crRightHandSideBoundedVector160 =             crRightHandSideBoundedVector32*crRightHandSideBoundedVector45;
const double crRightHandSideBoundedVector161 =             crRightHandSideBoundedVector11*crRightHandSideBoundedVector160;
const double crRightHandSideBoundedVector162 =             DN_DX_0_1*U_0_3 + DN_DX_1_1*U_1_3 + DN_DX_2_1*U_2_3;
const double crRightHandSideBoundedVector163 =             crRightHandSideBoundedVector162*crRightHandSideBoundedVector20;
const double crRightHandSideBoundedVector164 =             0.16666666666666666*ResProj_2_2 + crRightHandSideBoundedVector163 + 0.16666666666666666*dUdt_2_2;
const double crRightHandSideBoundedVector165 =             0.16666666666666666*ResProj_1_2 + 0.16666666666666666*dUdt_1_2;
const double crRightHandSideBoundedVector166 =             crRightHandSideBoundedVector27*(0.66666666666666663*ResProj_0_2 - crRightHandSideBoundedVector150 + crRightHandSideBoundedVector153 + crRightHandSideBoundedVector154 + crRightHandSideBoundedVector155 + crRightHandSideBoundedVector156*crRightHandSideBoundedVector56 - crRightHandSideBoundedVector158*crRightHandSideBoundedVector46 - crRightHandSideBoundedVector159*crRightHandSideBoundedVector161 + crRightHandSideBoundedVector164 + crRightHandSideBoundedVector165 + 0.66666666666666663*dUdt_0_2);
const double crRightHandSideBoundedVector167 =             0.16666666666666666*f_ext(0,1);
const double crRightHandSideBoundedVector168 =             crRightHandSideBoundedVector148 + crRightHandSideBoundedVector167 + 0.66666666666666663*f_ext(1,1);
const double crRightHandSideBoundedVector169 =             crRightHandSideBoundedVector168*crRightHandSideBoundedVector81;
const double crRightHandSideBoundedVector170 =             1.0000000000000002*crRightHandSideBoundedVector90;
const double crRightHandSideBoundedVector171 =             crRightHandSideBoundedVector86*(crRightHandSideBoundedVector102 - crRightHandSideBoundedVector170);
const double crRightHandSideBoundedVector172 =             crRightHandSideBoundedVector171*crRightHandSideBoundedVector71;
const double crRightHandSideBoundedVector173 =             crRightHandSideBoundedVector107*crRightHandSideBoundedVector68;
const double crRightHandSideBoundedVector174 =             crRightHandSideBoundedVector111*crRightHandSideBoundedVector61;
const double crRightHandSideBoundedVector175 =             2.2500000000000004*crRightHandSideBoundedVector110;
const double crRightHandSideBoundedVector176 =             crRightHandSideBoundedVector106*crRightHandSideBoundedVector32;
const double crRightHandSideBoundedVector177 =             crRightHandSideBoundedVector176*crRightHandSideBoundedVector86;
const double crRightHandSideBoundedVector178 =             0.16666666666666666*ResProj_0_2 + 0.16666666666666666*dUdt_0_2;
const double crRightHandSideBoundedVector179 =             crRightHandSideBoundedVector97*(0.66666666666666663*ResProj_1_2 - crRightHandSideBoundedVector107*crRightHandSideBoundedVector158 + crRightHandSideBoundedVector111*crRightHandSideBoundedVector156 + crRightHandSideBoundedVector164 - crRightHandSideBoundedVector169 + crRightHandSideBoundedVector172 + crRightHandSideBoundedVector173 + crRightHandSideBoundedVector174 - crRightHandSideBoundedVector175*crRightHandSideBoundedVector177 + crRightHandSideBoundedVector178 + 0.66666666666666663*dUdt_1_2);
const double crRightHandSideBoundedVector180 =             crRightHandSideBoundedVector147 + crRightHandSideBoundedVector167 + 0.66666666666666663*f_ext(2,1);
const double crRightHandSideBoundedVector181 =             crRightHandSideBoundedVector118*crRightHandSideBoundedVector180;
const double crRightHandSideBoundedVector182 =             1.0000000000000002*crRightHandSideBoundedVector123;
const double crRightHandSideBoundedVector183 =             crRightHandSideBoundedVector121*(crRightHandSideBoundedVector133 - crRightHandSideBoundedVector182);
const double crRightHandSideBoundedVector184 =             crRightHandSideBoundedVector183*crRightHandSideBoundedVector71;
const double crRightHandSideBoundedVector185 =             crRightHandSideBoundedVector137*crRightHandSideBoundedVector68;
const double crRightHandSideBoundedVector186 =             crRightHandSideBoundedVector140*crRightHandSideBoundedVector61;
const double crRightHandSideBoundedVector187 =             2.2500000000000004*crRightHandSideBoundedVector139;
const double crRightHandSideBoundedVector188 =             crRightHandSideBoundedVector136*crRightHandSideBoundedVector32;
const double crRightHandSideBoundedVector189 =             crRightHandSideBoundedVector121*crRightHandSideBoundedVector188;
const double crRightHandSideBoundedVector190 =             crRightHandSideBoundedVector129*(0.66666666666666663*ResProj_2_2 - crRightHandSideBoundedVector137*crRightHandSideBoundedVector158 + crRightHandSideBoundedVector140*crRightHandSideBoundedVector156 + crRightHandSideBoundedVector163 + crRightHandSideBoundedVector165 + crRightHandSideBoundedVector178 - crRightHandSideBoundedVector181 + crRightHandSideBoundedVector184 + crRightHandSideBoundedVector185 + crRightHandSideBoundedVector186 - crRightHandSideBoundedVector187*crRightHandSideBoundedVector189 + 0.66666666666666663*dUdt_2_2);
const double crRightHandSideBoundedVector191 =             crRightHandSideBoundedVector41 + crRightHandSideBoundedVector61;
const double crRightHandSideBoundedVector192 =             crRightHandSideBoundedVector3*crRightHandSideBoundedVector61;
const double crRightHandSideBoundedVector193 =             -crRightHandSideBoundedVector160 + crRightHandSideBoundedVector192;
const double crRightHandSideBoundedVector194 =             crRightHandSideBoundedVector11*mu;
const double crRightHandSideBoundedVector195 =             4.5000000000000009*crRightHandSideBoundedVector194;
const double crRightHandSideBoundedVector196 =             crRightHandSideBoundedVector3*crRightHandSideBoundedVector41;
const double crRightHandSideBoundedVector197 =             crRightHandSideBoundedVector196 - crRightHandSideBoundedVector72;
const double crRightHandSideBoundedVector198 =             -1.5000000000000002*crRightHandSideBoundedVector194*(crRightHandSideBoundedVector193 + crRightHandSideBoundedVector197);
const double crRightHandSideBoundedVector199 =             0.66666666666666663*crRightHandSideBoundedVector160;
const double crRightHandSideBoundedVector200 =             crRightHandSideBoundedVector4*(-0.66666666666666663*crRightHandSideBoundedVector192 + 1.3333333333333335*crRightHandSideBoundedVector196 + crRightHandSideBoundedVector199 - 1.3333333333333335*crRightHandSideBoundedVector72);
const double crRightHandSideBoundedVector201 =             1.0/lin_m_norm;
const double crRightHandSideBoundedVector202 =             crRightHandSideBoundedVector201*lin_m[0]*lin_m[1];
const double crRightHandSideBoundedVector203 =             crRightHandSideBoundedVector202*nu_st;
const double crRightHandSideBoundedVector204 =             0.66666666666666663*crRightHandSideBoundedVector72;
const double crRightHandSideBoundedVector205 =             crRightHandSideBoundedVector4*(-1.3333333333333335*crRightHandSideBoundedVector160 + 1.3333333333333335*crRightHandSideBoundedVector192 - 0.66666666666666663*crRightHandSideBoundedVector196 + crRightHandSideBoundedVector204);
const double crRightHandSideBoundedVector206 =             crRightHandSideBoundedVector201*pow(lin_m[0], 2);
const double crRightHandSideBoundedVector207 =             crRightHandSideBoundedVector206*nu_st;
const double crRightHandSideBoundedVector208 =             crRightHandSideBoundedVector202*nu_sc;
const double crRightHandSideBoundedVector209 =             1 - crRightHandSideBoundedVector206;
const double crRightHandSideBoundedVector210 =             crRightHandSideBoundedVector209*nu_sc;
const double crRightHandSideBoundedVector211 =             crRightHandSideBoundedVector193*crRightHandSideBoundedVector195 + crRightHandSideBoundedVector198 + crRightHandSideBoundedVector200*crRightHandSideBoundedVector203 - crRightHandSideBoundedVector200*crRightHandSideBoundedVector208 + crRightHandSideBoundedVector205*crRightHandSideBoundedVector207 + crRightHandSideBoundedVector205*crRightHandSideBoundedVector210;
const double crRightHandSideBoundedVector212 =             crRightHandSideBoundedVector61*crRightHandSideBoundedVector81;
const double crRightHandSideBoundedVector213 =             -crRightHandSideBoundedVector176 + crRightHandSideBoundedVector212;
const double crRightHandSideBoundedVector214 =             crRightHandSideBoundedVector86*mu;
const double crRightHandSideBoundedVector215 =             4.5000000000000009*crRightHandSideBoundedVector214;
const double crRightHandSideBoundedVector216 =             crRightHandSideBoundedVector41*crRightHandSideBoundedVector81;
const double crRightHandSideBoundedVector217 =             -crRightHandSideBoundedVector114 + crRightHandSideBoundedVector216;
const double crRightHandSideBoundedVector218 =             -1.5000000000000002*crRightHandSideBoundedVector214*(crRightHandSideBoundedVector213 + crRightHandSideBoundedVector217);
const double crRightHandSideBoundedVector219 =             0.66666666666666663*crRightHandSideBoundedVector176;
const double crRightHandSideBoundedVector220 =             crRightHandSideBoundedVector82*(-1.3333333333333335*crRightHandSideBoundedVector114 - 0.66666666666666663*crRightHandSideBoundedVector212 + 1.3333333333333335*crRightHandSideBoundedVector216 + crRightHandSideBoundedVector219);
const double crRightHandSideBoundedVector221 =             0.66666666666666663*crRightHandSideBoundedVector114;
const double crRightHandSideBoundedVector222 =             crRightHandSideBoundedVector82*(-1.3333333333333335*crRightHandSideBoundedVector176 + 1.3333333333333335*crRightHandSideBoundedVector212 - 0.66666666666666663*crRightHandSideBoundedVector216 + crRightHandSideBoundedVector221);
const double crRightHandSideBoundedVector223 =             crRightHandSideBoundedVector203*crRightHandSideBoundedVector220 + crRightHandSideBoundedVector207*crRightHandSideBoundedVector222 - crRightHandSideBoundedVector208*crRightHandSideBoundedVector220 + crRightHandSideBoundedVector210*crRightHandSideBoundedVector222 + crRightHandSideBoundedVector213*crRightHandSideBoundedVector215 + crRightHandSideBoundedVector218;
const double crRightHandSideBoundedVector224 =             crRightHandSideBoundedVector118*crRightHandSideBoundedVector61;
const double crRightHandSideBoundedVector225 =             -crRightHandSideBoundedVector188 + crRightHandSideBoundedVector224;
const double crRightHandSideBoundedVector226 =             crRightHandSideBoundedVector121*mu;
const double crRightHandSideBoundedVector227 =             4.5000000000000009*crRightHandSideBoundedVector226;
const double crRightHandSideBoundedVector228 =             crRightHandSideBoundedVector118*crRightHandSideBoundedVector41;
const double crRightHandSideBoundedVector229 =             -crRightHandSideBoundedVector143 + crRightHandSideBoundedVector228;
const double crRightHandSideBoundedVector230 =             -1.5000000000000002*crRightHandSideBoundedVector226*(crRightHandSideBoundedVector225 + crRightHandSideBoundedVector229);
const double crRightHandSideBoundedVector231 =             0.66666666666666663*crRightHandSideBoundedVector188;
const double crRightHandSideBoundedVector232 =             crRightHandSideBoundedVector119*(-1.3333333333333335*crRightHandSideBoundedVector143 - 0.66666666666666663*crRightHandSideBoundedVector224 + 1.3333333333333335*crRightHandSideBoundedVector228 + crRightHandSideBoundedVector231);
const double crRightHandSideBoundedVector233 =             0.66666666666666663*crRightHandSideBoundedVector143;
const double crRightHandSideBoundedVector234 =             crRightHandSideBoundedVector119*(-1.3333333333333335*crRightHandSideBoundedVector188 + 1.3333333333333335*crRightHandSideBoundedVector224 - 0.66666666666666663*crRightHandSideBoundedVector228 + crRightHandSideBoundedVector233);
const double crRightHandSideBoundedVector235 =             crRightHandSideBoundedVector203*crRightHandSideBoundedVector232 + crRightHandSideBoundedVector207*crRightHandSideBoundedVector234 - crRightHandSideBoundedVector208*crRightHandSideBoundedVector232 + crRightHandSideBoundedVector210*crRightHandSideBoundedVector234 + crRightHandSideBoundedVector225*crRightHandSideBoundedVector227 + crRightHandSideBoundedVector230;
const double crRightHandSideBoundedVector236 =             DN_DX_0_1*crRightHandSideBoundedVector11;
const double crRightHandSideBoundedVector237 =             nu_sc*(1 - crRightHandSideBoundedVector202);
const double crRightHandSideBoundedVector238 =             crRightHandSideBoundedVector3*crRightHandSideBoundedVector68;
const double crRightHandSideBoundedVector239 =             crRightHandSideBoundedVector3*crRightHandSideBoundedVector51;
const double crRightHandSideBoundedVector240 =             crRightHandSideBoundedVector45*crRightHandSideBoundedVector71;
const double crRightHandSideBoundedVector241 =             crRightHandSideBoundedVector32*crRightHandSideBoundedVector55;
const double crRightHandSideBoundedVector242 =             (crRightHandSideBoundedVector203*crRightHandSideBoundedVector3 + crRightHandSideBoundedVector237*crRightHandSideBoundedVector3 + mu)*(2.2500000000000004*crRightHandSideBoundedVector238 + 2.2500000000000004*crRightHandSideBoundedVector239 - 2.2500000000000004*crRightHandSideBoundedVector240 - 2.2500000000000004*crRightHandSideBoundedVector241);
const double crRightHandSideBoundedVector243 =             DN_DX_0_1*crRightHandSideBoundedVector86;
const double crRightHandSideBoundedVector244 =             crRightHandSideBoundedVector68*crRightHandSideBoundedVector81;
const double crRightHandSideBoundedVector245 =             crRightHandSideBoundedVector51*crRightHandSideBoundedVector81;
const double crRightHandSideBoundedVector246 =             crRightHandSideBoundedVector106*crRightHandSideBoundedVector71;
const double crRightHandSideBoundedVector247 =             crRightHandSideBoundedVector110*crRightHandSideBoundedVector32;
const double crRightHandSideBoundedVector248 =             (crRightHandSideBoundedVector203*crRightHandSideBoundedVector81 + crRightHandSideBoundedVector237*crRightHandSideBoundedVector81 + mu)*(2.2500000000000004*crRightHandSideBoundedVector244 + 2.2500000000000004*crRightHandSideBoundedVector245 - 2.2500000000000004*crRightHandSideBoundedVector246 - 2.2500000000000004*crRightHandSideBoundedVector247);
const double crRightHandSideBoundedVector249 =             DN_DX_0_1*crRightHandSideBoundedVector121;
const double crRightHandSideBoundedVector250 =             crRightHandSideBoundedVector118*crRightHandSideBoundedVector68;
const double crRightHandSideBoundedVector251 =             crRightHandSideBoundedVector118*crRightHandSideBoundedVector51;
const double crRightHandSideBoundedVector252 =             crRightHandSideBoundedVector136*crRightHandSideBoundedVector71;
const double crRightHandSideBoundedVector253 =             crRightHandSideBoundedVector139*crRightHandSideBoundedVector32;
const double crRightHandSideBoundedVector254 =             (crRightHandSideBoundedVector118*crRightHandSideBoundedVector203 + crRightHandSideBoundedVector118*crRightHandSideBoundedVector237 + mu)*(2.2500000000000004*crRightHandSideBoundedVector250 + 2.2500000000000004*crRightHandSideBoundedVector251 - 2.2500000000000004*crRightHandSideBoundedVector252 - 2.2500000000000004*crRightHandSideBoundedVector253);
const double crRightHandSideBoundedVector255 =             DN_DX_0_1*crRightHandSideBoundedVector46;
const double crRightHandSideBoundedVector256 =             DN_DX_0_1*crRightHandSideBoundedVector42 + 0.66666666666666663*crRightHandSideBoundedVector49 + 0.66666666666666663*crRightHandSideBoundedVector50;
const double crRightHandSideBoundedVector257 =             0.66666666666666663*crRightHandSideBoundedVector4;
const double crRightHandSideBoundedVector258 =             DN_DX_0_0*crRightHandSideBoundedVector56;
const double crRightHandSideBoundedVector259 =             1.0*crRightHandSideBoundedVector20;
const double crRightHandSideBoundedVector260 =             1.5000000000000002*crRightHandSideBoundedVector11;
const double crRightHandSideBoundedVector261 =             1.5000000000000002*crRightHandSideBoundedVector20;
const double crRightHandSideBoundedVector262 =             crRightHandSideBoundedVector11*crRightHandSideBoundedVector261;
const double crRightHandSideBoundedVector263 =             1.0*h;
const double crRightHandSideBoundedVector264 =             crRightHandSideBoundedVector166*crRightHandSideBoundedVector263;
const double crRightHandSideBoundedVector265 =             DN_DX_0_1*crRightHandSideBoundedVector107;
const double crRightHandSideBoundedVector266 =             DN_DX_0_0*crRightHandSideBoundedVector111;
const double crRightHandSideBoundedVector267 =             DN_DX_1_1*crRightHandSideBoundedVector43 + DN_DX_2_1*crRightHandSideBoundedVector44 + 0.16666666666666666*crRightHandSideBoundedVector48;
const double crRightHandSideBoundedVector268 =             0.16666666666666666*crRightHandSideBoundedVector69;
const double crRightHandSideBoundedVector269 =             0.37500000000000006*crRightHandSideBoundedVector86;
const double crRightHandSideBoundedVector270 =             crRightHandSideBoundedVector246*crRightHandSideBoundedVector269;
const double crRightHandSideBoundedVector271 =             crRightHandSideBoundedVector247*crRightHandSideBoundedVector269;
const double crRightHandSideBoundedVector272 =             -crRightHandSideBoundedVector20*crRightHandSideBoundedVector271 - crRightHandSideBoundedVector267*crRightHandSideBoundedVector82 + crRightHandSideBoundedVector268*crRightHandSideBoundedVector82 + crRightHandSideBoundedVector270;
const double crRightHandSideBoundedVector273 =             crRightHandSideBoundedVector179*crRightHandSideBoundedVector263;
const double crRightHandSideBoundedVector274 =             DN_DX_0_1*crRightHandSideBoundedVector137;
const double crRightHandSideBoundedVector275 =             DN_DX_0_0*crRightHandSideBoundedVector140;
const double crRightHandSideBoundedVector276 =             0.37500000000000006*crRightHandSideBoundedVector121;
const double crRightHandSideBoundedVector277 =             crRightHandSideBoundedVector252*crRightHandSideBoundedVector276;
const double crRightHandSideBoundedVector278 =             crRightHandSideBoundedVector253*crRightHandSideBoundedVector276;
const double crRightHandSideBoundedVector279 =             -crRightHandSideBoundedVector119*crRightHandSideBoundedVector267 + crRightHandSideBoundedVector119*crRightHandSideBoundedVector268 - crRightHandSideBoundedVector20*crRightHandSideBoundedVector278 + crRightHandSideBoundedVector277;
const double crRightHandSideBoundedVector280 =             crRightHandSideBoundedVector190*crRightHandSideBoundedVector263;
const double crRightHandSideBoundedVector281 =             0.66666666666666663*crRightHandSideBoundedVector69;
const double crRightHandSideBoundedVector282 =             crRightHandSideBoundedVector62 - 3.0;
const double crRightHandSideBoundedVector283 =             0.66666666666666663*crRightHandSideBoundedVector61;
const double crRightHandSideBoundedVector284 =             crRightHandSideBoundedVector282*crRightHandSideBoundedVector283;
const double crRightHandSideBoundedVector285 =             1.5000000000000002*crRightHandSideBoundedVector72;
const double crRightHandSideBoundedVector286 =             crRightHandSideBoundedVector11*crRightHandSideBoundedVector285;
const double crRightHandSideBoundedVector287 =             crRightHandSideBoundedVector282*crRightHandSideBoundedVector4;
const double crRightHandSideBoundedVector288 =             DN_DX_0_0*crRightHandSideBoundedVector46;
const double crRightHandSideBoundedVector289 =             2.2500000000000004*gamma - 6.7500000000000018;
const double crRightHandSideBoundedVector290 =             crRightHandSideBoundedVector11*crRightHandSideBoundedVector289;
const double crRightHandSideBoundedVector291 =             0.66666666666666663*crRightHandSideBoundedVector38 + 0.66666666666666663*crRightHandSideBoundedVector39 + 0.66666666666666663*crRightHandSideBoundedVector40;
const double crRightHandSideBoundedVector292 =             DN_DX_0_1*crRightHandSideBoundedVector56;
const double crRightHandSideBoundedVector293 =             -crRightHandSideBoundedVector286 + crRightHandSideBoundedVector291*crRightHandSideBoundedVector4 + crRightHandSideBoundedVector292;
const double crRightHandSideBoundedVector294 =             crRightHandSideBoundedVector263*crRightHandSideBoundedVector79;
const double crRightHandSideBoundedVector295 =             DN_DX_0_0*crRightHandSideBoundedVector107;
const double crRightHandSideBoundedVector296 =             DN_DX_0_1*crRightHandSideBoundedVector111;
const double crRightHandSideBoundedVector297 =             0.16666666666666666*crRightHandSideBoundedVector38 + 0.16666666666666666*crRightHandSideBoundedVector39 + 0.16666666666666666*crRightHandSideBoundedVector40;
const double crRightHandSideBoundedVector298 =             crRightHandSideBoundedVector114*crRightHandSideBoundedVector269;
const double crRightHandSideBoundedVector299 =             crRightHandSideBoundedVector297*crRightHandSideBoundedVector82 - crRightHandSideBoundedVector298;
const double crRightHandSideBoundedVector300 =             crRightHandSideBoundedVector296 + crRightHandSideBoundedVector299;
const double crRightHandSideBoundedVector301 =             0.16666666666666666*crRightHandSideBoundedVector61;
const double crRightHandSideBoundedVector302 =             crRightHandSideBoundedVector282*crRightHandSideBoundedVector301;
const double crRightHandSideBoundedVector303 =             0.16666666666666666*crRightHandSideBoundedVector289;
const double crRightHandSideBoundedVector304 =             crRightHandSideBoundedVector177*crRightHandSideBoundedVector303 - crRightHandSideBoundedVector302*crRightHandSideBoundedVector82;
const double crRightHandSideBoundedVector305 =             crRightHandSideBoundedVector117*crRightHandSideBoundedVector263;
const double crRightHandSideBoundedVector306 =             DN_DX_0_0*crRightHandSideBoundedVector137;
const double crRightHandSideBoundedVector307 =             DN_DX_0_1*crRightHandSideBoundedVector140;
const double crRightHandSideBoundedVector308 =             crRightHandSideBoundedVector143*crRightHandSideBoundedVector276;
const double crRightHandSideBoundedVector309 =             crRightHandSideBoundedVector119*crRightHandSideBoundedVector297 - crRightHandSideBoundedVector308;
const double crRightHandSideBoundedVector310 =             crRightHandSideBoundedVector307 + crRightHandSideBoundedVector309;
const double crRightHandSideBoundedVector311 =             -crRightHandSideBoundedVector119*crRightHandSideBoundedVector302 + crRightHandSideBoundedVector189*crRightHandSideBoundedVector303;
const double crRightHandSideBoundedVector312 =             crRightHandSideBoundedVector145*crRightHandSideBoundedVector263;
const double crRightHandSideBoundedVector313 =             crRightHandSideBoundedVector0*crRightHandSideBoundedVector259;
const double crRightHandSideBoundedVector314 =             lambda/c_v;
const double crRightHandSideBoundedVector315 =             crRightHandSideBoundedVector314/gamma;
const double crRightHandSideBoundedVector316 =             crRightHandSideBoundedVector11*crRightHandSideBoundedVector20;
const double crRightHandSideBoundedVector317 =             crRightHandSideBoundedVector316*crRightHandSideBoundedVector33;
const double crRightHandSideBoundedVector318 =             crRightHandSideBoundedVector20*(crRightHandSideBoundedVector24 - crRightHandSideBoundedVector4*(0.22222222222222221*crRightHandSideBoundedVector14 + 0.22222222222222221*crRightHandSideBoundedVector17));
const double crRightHandSideBoundedVector319 =             crRightHandSideBoundedVector4*(crRightHandSideBoundedVector24 + crRightHandSideBoundedVector318);
const double crRightHandSideBoundedVector320 =             crRightHandSideBoundedVector151*crRightHandSideBoundedVector316;
const double crRightHandSideBoundedVector321 =             0.16666666666666666*r[1];
const double crRightHandSideBoundedVector322 =             0.16666666666666666*r[2];
const double crRightHandSideBoundedVector323 =             crRightHandSideBoundedVector3*(crRightHandSideBoundedVector321 + crRightHandSideBoundedVector322 + 0.66666666666666663*r[0]);
const double crRightHandSideBoundedVector324 =             crRightHandSideBoundedVector30*crRightHandSideBoundedVector45;
const double crRightHandSideBoundedVector325 =             crRightHandSideBoundedVector149*crRightHandSideBoundedVector55;
const double crRightHandSideBoundedVector326 =             crRightHandSideBoundedVector75*gamma;
const double crRightHandSideBoundedVector327 =             crRightHandSideBoundedVector326*crRightHandSideBoundedVector46;
const double crRightHandSideBoundedVector328 =             crRightHandSideBoundedVector162*gamma;
const double crRightHandSideBoundedVector329 =             crRightHandSideBoundedVector328*crRightHandSideBoundedVector56;
const double crRightHandSideBoundedVector330 =             -0.37500000000000006*U_1_3;
const double crRightHandSideBoundedVector331 =             -0.37500000000000006*U_2_3;
const double crRightHandSideBoundedVector332 =             0.75000000000000011*crRightHandSideBoundedVector20;
const double crRightHandSideBoundedVector333 =             1.0/crRightHandSideBoundedVector10;
const double crRightHandSideBoundedVector334 =             crRightHandSideBoundedVector18*crRightHandSideBoundedVector333;
const double crRightHandSideBoundedVector335 =             -1.5000000000000002*U_0_3 - 2.2500000000000004*crRightHandSideBoundedVector318 + crRightHandSideBoundedVector330 + crRightHandSideBoundedVector331 + crRightHandSideBoundedVector332*crRightHandSideBoundedVector334;
const double crRightHandSideBoundedVector336 =             crRightHandSideBoundedVector11*crRightHandSideBoundedVector335;
const double crRightHandSideBoundedVector337 =             crRightHandSideBoundedVector160*crRightHandSideBoundedVector336;
const double crRightHandSideBoundedVector338 =             crRightHandSideBoundedVector336*crRightHandSideBoundedVector72;
const double crRightHandSideBoundedVector339 =             crRightHandSideBoundedVector55*crRightHandSideBoundedVector74;
const double crRightHandSideBoundedVector340 =             crRightHandSideBoundedVector51*crRightHandSideBoundedVector73;
const double crRightHandSideBoundedVector341 =             crRightHandSideBoundedVector316*crRightHandSideBoundedVector55;
const double crRightHandSideBoundedVector342 =             0.16666666666666666*ResProj_2_3 + 0.16666666666666666*dUdt_2_3;
const double crRightHandSideBoundedVector343 =             0.16666666666666666*ResProj_1_3 + 0.16666666666666666*dUdt_1_3;
const double crRightHandSideBoundedVector344 =             (0.66666666666666663*ResProj_0_3 - crRightHandSideBoundedVector323 - crRightHandSideBoundedVector324 - crRightHandSideBoundedVector325 + crRightHandSideBoundedVector327 + crRightHandSideBoundedVector329 + crRightHandSideBoundedVector337 + crRightHandSideBoundedVector338 - crRightHandSideBoundedVector339*crRightHandSideBoundedVector69 - crRightHandSideBoundedVector340*crRightHandSideBoundedVector341 + crRightHandSideBoundedVector342 + crRightHandSideBoundedVector343 + crRightHandSideBoundedVector41*(crRightHandSideBoundedVector319 - crRightHandSideBoundedVector320) + crRightHandSideBoundedVector61*(-crRightHandSideBoundedVector317 + crRightHandSideBoundedVector319) + 0.66666666666666663*dUdt_0_3)/(crRightHandSideBoundedVector26 + crRightHandSideBoundedVector315*crRightHandSideBoundedVector6);
const double crRightHandSideBoundedVector345 =             crRightHandSideBoundedVector315*crRightHandSideBoundedVector5;
const double crRightHandSideBoundedVector346 =             crRightHandSideBoundedVector20*crRightHandSideBoundedVector86;
const double crRightHandSideBoundedVector347 =             crRightHandSideBoundedVector101*crRightHandSideBoundedVector346;
const double crRightHandSideBoundedVector348 =             crRightHandSideBoundedVector20*(-crRightHandSideBoundedVector82*(0.22222222222222221*crRightHandSideBoundedVector88 + 0.22222222222222221*crRightHandSideBoundedVector90) + crRightHandSideBoundedVector94);
const double crRightHandSideBoundedVector349 =             crRightHandSideBoundedVector82*(crRightHandSideBoundedVector348 + crRightHandSideBoundedVector94);
const double crRightHandSideBoundedVector350 =             crRightHandSideBoundedVector170*crRightHandSideBoundedVector346;
const double crRightHandSideBoundedVector351 =             0.16666666666666666*r[0];
const double crRightHandSideBoundedVector352 =             crRightHandSideBoundedVector81*(crRightHandSideBoundedVector322 + crRightHandSideBoundedVector351 + 0.66666666666666663*r[1]);
const double crRightHandSideBoundedVector353 =             crRightHandSideBoundedVector106*crRightHandSideBoundedVector99;
const double crRightHandSideBoundedVector354 =             crRightHandSideBoundedVector110*crRightHandSideBoundedVector168;
const double crRightHandSideBoundedVector355 =             crRightHandSideBoundedVector107*crRightHandSideBoundedVector326;
const double crRightHandSideBoundedVector356 =             crRightHandSideBoundedVector111*crRightHandSideBoundedVector328;
const double crRightHandSideBoundedVector357 =             -0.37500000000000006*U_0_3;
const double crRightHandSideBoundedVector358 =             1.0/crRightHandSideBoundedVector85;
const double crRightHandSideBoundedVector359 =             crRightHandSideBoundedVector358*crRightHandSideBoundedVector91;
const double crRightHandSideBoundedVector360 =             -1.5000000000000002*U_1_3 + crRightHandSideBoundedVector331 + crRightHandSideBoundedVector332*crRightHandSideBoundedVector359 - 2.2500000000000004*crRightHandSideBoundedVector348 + crRightHandSideBoundedVector357;
const double crRightHandSideBoundedVector361 =             crRightHandSideBoundedVector360*crRightHandSideBoundedVector86;
const double crRightHandSideBoundedVector362 =             crRightHandSideBoundedVector176*crRightHandSideBoundedVector361;
const double crRightHandSideBoundedVector363 =             crRightHandSideBoundedVector114*crRightHandSideBoundedVector361;
const double crRightHandSideBoundedVector364 =             crRightHandSideBoundedVector110*crRightHandSideBoundedVector86;
const double crRightHandSideBoundedVector365 =             crRightHandSideBoundedVector113*crRightHandSideBoundedVector364;
const double crRightHandSideBoundedVector366 =             crRightHandSideBoundedVector113*crRightHandSideBoundedVector51;
const double crRightHandSideBoundedVector367 =             crRightHandSideBoundedVector110*crRightHandSideBoundedVector346;
const double crRightHandSideBoundedVector368 =             0.16666666666666666*ResProj_0_3 + 0.16666666666666666*dUdt_0_3;
const double crRightHandSideBoundedVector369 =             (0.66666666666666663*ResProj_1_3 + crRightHandSideBoundedVector342 - crRightHandSideBoundedVector352 - crRightHandSideBoundedVector353 - crRightHandSideBoundedVector354 + crRightHandSideBoundedVector355 + crRightHandSideBoundedVector356 + crRightHandSideBoundedVector362 + crRightHandSideBoundedVector363 - crRightHandSideBoundedVector365*crRightHandSideBoundedVector69 - crRightHandSideBoundedVector366*crRightHandSideBoundedVector367 + crRightHandSideBoundedVector368 + crRightHandSideBoundedVector41*(crRightHandSideBoundedVector349 - crRightHandSideBoundedVector350) + crRightHandSideBoundedVector61*(-crRightHandSideBoundedVector347 + crRightHandSideBoundedVector349) + 0.66666666666666663*dUdt_1_3)/(crRightHandSideBoundedVector345*crRightHandSideBoundedVector82 + crRightHandSideBoundedVector96);
const double crRightHandSideBoundedVector370 =             crRightHandSideBoundedVector121*crRightHandSideBoundedVector20;
const double crRightHandSideBoundedVector371 =             crRightHandSideBoundedVector132*crRightHandSideBoundedVector370;
const double crRightHandSideBoundedVector372 =             crRightHandSideBoundedVector20*(-crRightHandSideBoundedVector119*(0.22222222222222221*crRightHandSideBoundedVector122 + 0.22222222222222221*crRightHandSideBoundedVector123) + crRightHandSideBoundedVector126);
const double crRightHandSideBoundedVector373 =             crRightHandSideBoundedVector119*(crRightHandSideBoundedVector126 + crRightHandSideBoundedVector372);
const double crRightHandSideBoundedVector374 =             crRightHandSideBoundedVector182*crRightHandSideBoundedVector370;
const double crRightHandSideBoundedVector375 =             crRightHandSideBoundedVector118*(crRightHandSideBoundedVector321 + crRightHandSideBoundedVector351 + 0.66666666666666663*r[2]);
const double crRightHandSideBoundedVector376 =             crRightHandSideBoundedVector130*crRightHandSideBoundedVector136;
const double crRightHandSideBoundedVector377 =             crRightHandSideBoundedVector139*crRightHandSideBoundedVector180;
const double crRightHandSideBoundedVector378 =             crRightHandSideBoundedVector137*crRightHandSideBoundedVector326;
const double crRightHandSideBoundedVector379 =             crRightHandSideBoundedVector140*crRightHandSideBoundedVector328;
const double crRightHandSideBoundedVector380 =             1.0/crRightHandSideBoundedVector120;
const double crRightHandSideBoundedVector381 =             crRightHandSideBoundedVector124*crRightHandSideBoundedVector380;
const double crRightHandSideBoundedVector382 =             -1.5000000000000002*U_2_3 + crRightHandSideBoundedVector330 + crRightHandSideBoundedVector332*crRightHandSideBoundedVector381 + crRightHandSideBoundedVector357 - 2.2500000000000004*crRightHandSideBoundedVector372;
const double crRightHandSideBoundedVector383 =             crRightHandSideBoundedVector121*crRightHandSideBoundedVector382;
const double crRightHandSideBoundedVector384 =             crRightHandSideBoundedVector188*crRightHandSideBoundedVector383;
const double crRightHandSideBoundedVector385 =             crRightHandSideBoundedVector143*crRightHandSideBoundedVector383;
const double crRightHandSideBoundedVector386 =             crRightHandSideBoundedVector121*crRightHandSideBoundedVector139;
const double crRightHandSideBoundedVector387 =             crRightHandSideBoundedVector142*crRightHandSideBoundedVector386;
const double crRightHandSideBoundedVector388 =             crRightHandSideBoundedVector142*crRightHandSideBoundedVector51;
const double crRightHandSideBoundedVector389 =             crRightHandSideBoundedVector139*crRightHandSideBoundedVector370;
const double crRightHandSideBoundedVector390 =             (0.66666666666666663*ResProj_2_3 + crRightHandSideBoundedVector343 + crRightHandSideBoundedVector368 - crRightHandSideBoundedVector375 - crRightHandSideBoundedVector376 - crRightHandSideBoundedVector377 + crRightHandSideBoundedVector378 + crRightHandSideBoundedVector379 + crRightHandSideBoundedVector384 + crRightHandSideBoundedVector385 - crRightHandSideBoundedVector387*crRightHandSideBoundedVector69 - crRightHandSideBoundedVector388*crRightHandSideBoundedVector389 + crRightHandSideBoundedVector41*(crRightHandSideBoundedVector373 - crRightHandSideBoundedVector374) + crRightHandSideBoundedVector61*(-crRightHandSideBoundedVector371 + crRightHandSideBoundedVector373) + 0.66666666666666663*dUdt_2_3)/(crRightHandSideBoundedVector119*crRightHandSideBoundedVector345 + crRightHandSideBoundedVector128);
const double crRightHandSideBoundedVector391 =             pow(crRightHandSideBoundedVector10, -3);
const double crRightHandSideBoundedVector392 =             0.66666666666666663*crRightHandSideBoundedVector391;
const double crRightHandSideBoundedVector393 =             3.0000000000000004*crRightHandSideBoundedVector14;
const double crRightHandSideBoundedVector394 =             crRightHandSideBoundedVector18*crRightHandSideBoundedVector261;
const double crRightHandSideBoundedVector395 =             crRightHandSideBoundedVector32*(-crRightHandSideBoundedVector393 + crRightHandSideBoundedVector394);
const double crRightHandSideBoundedVector396 =             crRightHandSideBoundedVector260*crRightHandSideBoundedVector45;
const double crRightHandSideBoundedVector397 =             crRightHandSideBoundedVector260*crRightHandSideBoundedVector55;
const double crRightHandSideBoundedVector398 =             4.5*crRightHandSideBoundedVector391;
const double crRightHandSideBoundedVector399 =             crRightHandSideBoundedVector45*crRightHandSideBoundedVector72;
const double crRightHandSideBoundedVector400 =             crRightHandSideBoundedVector398*crRightHandSideBoundedVector399;
const double crRightHandSideBoundedVector401 =             crRightHandSideBoundedVector290*crRightHandSideBoundedVector45;
const double crRightHandSideBoundedVector402 =             crRightHandSideBoundedVector262*crRightHandSideBoundedVector55;
const double crRightHandSideBoundedVector403 =             crRightHandSideBoundedVector402*crRightHandSideBoundedVector68;
const double crRightHandSideBoundedVector404 =             crRightHandSideBoundedVector236*crRightHandSideBoundedVector55;
const double crRightHandSideBoundedVector405 =             crRightHandSideBoundedVector404*crRightHandSideBoundedVector73;
const double crRightHandSideBoundedVector406 =             0.1111111111111111*f_ext(1,0);
const double crRightHandSideBoundedVector407 =             0.1111111111111111*f_ext(2,0);
const double crRightHandSideBoundedVector408 =             crRightHandSideBoundedVector406 + crRightHandSideBoundedVector407 + 0.44444444444444442*f_ext(0,0);
const double crRightHandSideBoundedVector409 =             0.16666666666666666*ResProj_2_0 + crRightHandSideBoundedVector191 + 0.16666666666666666*dUdt_2_0;
const double crRightHandSideBoundedVector410 =             0.16666666666666666*ResProj_1_0 + 0.16666666666666666*dUdt_1_0;
const double crRightHandSideBoundedVector411 =             crRightHandSideBoundedVector263/stab_c2;
const double crRightHandSideBoundedVector412 =             crRightHandSideBoundedVector411*(0.66666666666666663*ResProj_0_0 + crRightHandSideBoundedVector409 + crRightHandSideBoundedVector410 + 0.66666666666666663*dUdt_0_0)/crRightHandSideBoundedVector25;
const double crRightHandSideBoundedVector413 =             crRightHandSideBoundedVector110*crRightHandSideBoundedVector243;
const double crRightHandSideBoundedVector414 =             crRightHandSideBoundedVector113*crRightHandSideBoundedVector413;
const double crRightHandSideBoundedVector415 =             0.16666666666666666*crRightHandSideBoundedVector32;
const double crRightHandSideBoundedVector416 =             pow(crRightHandSideBoundedVector85, -3);
const double crRightHandSideBoundedVector417 =             3.0000000000000004*crRightHandSideBoundedVector88;
const double crRightHandSideBoundedVector418 =             crRightHandSideBoundedVector261*crRightHandSideBoundedVector91;
const double crRightHandSideBoundedVector419 =             crRightHandSideBoundedVector416*(-crRightHandSideBoundedVector417 + crRightHandSideBoundedVector418);
const double crRightHandSideBoundedVector420 =             crRightHandSideBoundedVector106*crRightHandSideBoundedVector269;
const double crRightHandSideBoundedVector421 =             crRightHandSideBoundedVector41*crRightHandSideBoundedVector420;
const double crRightHandSideBoundedVector422 =             crRightHandSideBoundedVector110*crRightHandSideBoundedVector269;
const double crRightHandSideBoundedVector423 =             crRightHandSideBoundedVector422*crRightHandSideBoundedVector51;
const double crRightHandSideBoundedVector424 =             1.125*crRightHandSideBoundedVector416;
const double crRightHandSideBoundedVector425 =             crRightHandSideBoundedVector106*crRightHandSideBoundedVector114;
const double crRightHandSideBoundedVector426 =             crRightHandSideBoundedVector424*crRightHandSideBoundedVector425;
const double crRightHandSideBoundedVector427 =             crRightHandSideBoundedVector289*crRightHandSideBoundedVector301;
const double crRightHandSideBoundedVector428 =             crRightHandSideBoundedVector106*crRightHandSideBoundedVector86;
const double crRightHandSideBoundedVector429 =             crRightHandSideBoundedVector422*crRightHandSideBoundedVector69;
const double crRightHandSideBoundedVector430 =             0.027777777777777776*f_ext(0,0);
const double crRightHandSideBoundedVector431 =             0.027777777777777776*f_ext(2,0);
const double crRightHandSideBoundedVector432 =             crRightHandSideBoundedVector406 + crRightHandSideBoundedVector430 + crRightHandSideBoundedVector431;
const double crRightHandSideBoundedVector433 =             -crRightHandSideBoundedVector415*crRightHandSideBoundedVector419 - crRightHandSideBoundedVector421 - crRightHandSideBoundedVector423 + crRightHandSideBoundedVector426 + crRightHandSideBoundedVector427*crRightHandSideBoundedVector428 + crRightHandSideBoundedVector429 + crRightHandSideBoundedVector432;
const double crRightHandSideBoundedVector434 =             0.16666666666666666*ResProj_0_0 + 0.16666666666666666*dUdt_0_0;
const double crRightHandSideBoundedVector435 =             crRightHandSideBoundedVector411*(0.66666666666666663*ResProj_1_0 + crRightHandSideBoundedVector409 + crRightHandSideBoundedVector434 + 0.66666666666666663*dUdt_1_0)/crRightHandSideBoundedVector95;
const double crRightHandSideBoundedVector436 =             crRightHandSideBoundedVector139*crRightHandSideBoundedVector249;
const double crRightHandSideBoundedVector437 =             crRightHandSideBoundedVector142*crRightHandSideBoundedVector436;
const double crRightHandSideBoundedVector438 =             pow(crRightHandSideBoundedVector120, -3);
const double crRightHandSideBoundedVector439 =             3.0000000000000004*crRightHandSideBoundedVector122;
const double crRightHandSideBoundedVector440 =             crRightHandSideBoundedVector124*crRightHandSideBoundedVector261;
const double crRightHandSideBoundedVector441 =             crRightHandSideBoundedVector438*(-crRightHandSideBoundedVector439 + crRightHandSideBoundedVector440);
const double crRightHandSideBoundedVector442 =             crRightHandSideBoundedVector136*crRightHandSideBoundedVector276;
const double crRightHandSideBoundedVector443 =             crRightHandSideBoundedVector41*crRightHandSideBoundedVector442;
const double crRightHandSideBoundedVector444 =             crRightHandSideBoundedVector139*crRightHandSideBoundedVector276;
const double crRightHandSideBoundedVector445 =             crRightHandSideBoundedVector444*crRightHandSideBoundedVector51;
const double crRightHandSideBoundedVector446 =             1.125*crRightHandSideBoundedVector438;
const double crRightHandSideBoundedVector447 =             crRightHandSideBoundedVector136*crRightHandSideBoundedVector143;
const double crRightHandSideBoundedVector448 =             crRightHandSideBoundedVector446*crRightHandSideBoundedVector447;
const double crRightHandSideBoundedVector449 =             crRightHandSideBoundedVector121*crRightHandSideBoundedVector136;
const double crRightHandSideBoundedVector450 =             crRightHandSideBoundedVector444*crRightHandSideBoundedVector69;
const double crRightHandSideBoundedVector451 =             0.027777777777777776*f_ext(1,0);
const double crRightHandSideBoundedVector452 =             crRightHandSideBoundedVector407 + crRightHandSideBoundedVector430 + crRightHandSideBoundedVector451;
const double crRightHandSideBoundedVector453 =             -crRightHandSideBoundedVector415*crRightHandSideBoundedVector441 + crRightHandSideBoundedVector427*crRightHandSideBoundedVector449 - crRightHandSideBoundedVector443 - crRightHandSideBoundedVector445 + crRightHandSideBoundedVector448 + crRightHandSideBoundedVector450 + crRightHandSideBoundedVector452;
const double crRightHandSideBoundedVector454 =             crRightHandSideBoundedVector411*(0.66666666666666663*ResProj_2_0 + crRightHandSideBoundedVector191 + crRightHandSideBoundedVector410 + crRightHandSideBoundedVector434 + 0.66666666666666663*dUdt_2_0)/crRightHandSideBoundedVector127;
const double crRightHandSideBoundedVector455 =             0.16666666666666666*crRightHandSideBoundedVector131 - 0.16666666666666666*crRightHandSideBoundedVector135 + crRightHandSideBoundedVector136*crRightHandSideBoundedVector308 + crRightHandSideBoundedVector137*crRightHandSideBoundedVector302 - 0.16666666666666666*crRightHandSideBoundedVector138 + crRightHandSideBoundedVector140*crRightHandSideBoundedVector268 - 0.16666666666666666*crRightHandSideBoundedVector141;
const double crRightHandSideBoundedVector456 =             0.16666666666666666*crRightHandSideBoundedVector100 - 0.16666666666666666*crRightHandSideBoundedVector104 + crRightHandSideBoundedVector106*crRightHandSideBoundedVector298 + crRightHandSideBoundedVector107*crRightHandSideBoundedVector302 - 0.16666666666666666*crRightHandSideBoundedVector108 + crRightHandSideBoundedVector111*crRightHandSideBoundedVector268 - 0.16666666666666666*crRightHandSideBoundedVector112;
const double crRightHandSideBoundedVector457 =             crRightHandSideBoundedVector201*pow(lin_m[1], 2);
const double crRightHandSideBoundedVector458 =             crRightHandSideBoundedVector457*nu_st;
const double crRightHandSideBoundedVector459 =             1 - crRightHandSideBoundedVector457;
const double crRightHandSideBoundedVector460 =             crRightHandSideBoundedVector459*nu_sc;
const double crRightHandSideBoundedVector461 =             crRightHandSideBoundedVector195*crRightHandSideBoundedVector197 + crRightHandSideBoundedVector198 + crRightHandSideBoundedVector200*crRightHandSideBoundedVector458 + crRightHandSideBoundedVector200*crRightHandSideBoundedVector460 + crRightHandSideBoundedVector203*crRightHandSideBoundedVector205 - crRightHandSideBoundedVector205*crRightHandSideBoundedVector208;
const double crRightHandSideBoundedVector462 =             crRightHandSideBoundedVector203*crRightHandSideBoundedVector222 - crRightHandSideBoundedVector208*crRightHandSideBoundedVector222 + crRightHandSideBoundedVector215*crRightHandSideBoundedVector217 + crRightHandSideBoundedVector218 + crRightHandSideBoundedVector220*crRightHandSideBoundedVector458 + crRightHandSideBoundedVector220*crRightHandSideBoundedVector460;
const double crRightHandSideBoundedVector463 =             crRightHandSideBoundedVector203*crRightHandSideBoundedVector234 - crRightHandSideBoundedVector208*crRightHandSideBoundedVector234 + crRightHandSideBoundedVector227*crRightHandSideBoundedVector229 + crRightHandSideBoundedVector230 + crRightHandSideBoundedVector232*crRightHandSideBoundedVector458 + crRightHandSideBoundedVector232*crRightHandSideBoundedVector460;
const double crRightHandSideBoundedVector464 =             crRightHandSideBoundedVector11*crRightHandSideBoundedVector242;
const double crRightHandSideBoundedVector465 =             crRightHandSideBoundedVector248*crRightHandSideBoundedVector86;
const double crRightHandSideBoundedVector466 =             crRightHandSideBoundedVector121*crRightHandSideBoundedVector254;
const double crRightHandSideBoundedVector467 =             0.66666666666666663*crRightHandSideBoundedVector157;
const double crRightHandSideBoundedVector468 =             0.66666666666666663*crRightHandSideBoundedVector41;
const double crRightHandSideBoundedVector469 =             crRightHandSideBoundedVector282*crRightHandSideBoundedVector468;
const double crRightHandSideBoundedVector470 =             1.5000000000000002*crRightHandSideBoundedVector160;
const double crRightHandSideBoundedVector471 =             crRightHandSideBoundedVector11*crRightHandSideBoundedVector470;
const double crRightHandSideBoundedVector472 =             0.66666666666666663*crRightHandSideBoundedVector58 + 0.66666666666666663*crRightHandSideBoundedVector59 + 0.66666666666666663*crRightHandSideBoundedVector60;
const double crRightHandSideBoundedVector473 =             crRightHandSideBoundedVector288 + crRightHandSideBoundedVector4*crRightHandSideBoundedVector472 - crRightHandSideBoundedVector471;
const double crRightHandSideBoundedVector474 =             DN_DX_0_0*crRightHandSideBoundedVector52 + 0.66666666666666663*crRightHandSideBoundedVector66 + 0.66666666666666663*crRightHandSideBoundedVector67;
const double crRightHandSideBoundedVector475 =             0.16666666666666666*crRightHandSideBoundedVector58 + 0.16666666666666666*crRightHandSideBoundedVector59 + 0.16666666666666666*crRightHandSideBoundedVector60;
const double crRightHandSideBoundedVector476 =             -crRightHandSideBoundedVector176*crRightHandSideBoundedVector269 + crRightHandSideBoundedVector475*crRightHandSideBoundedVector82;
const double crRightHandSideBoundedVector477 =             crRightHandSideBoundedVector295 + crRightHandSideBoundedVector476;
const double crRightHandSideBoundedVector478 =             0.16666666666666666*crRightHandSideBoundedVector41;
const double crRightHandSideBoundedVector479 =             crRightHandSideBoundedVector282*crRightHandSideBoundedVector478;
const double crRightHandSideBoundedVector480 =             crRightHandSideBoundedVector115*crRightHandSideBoundedVector303 - crRightHandSideBoundedVector479*crRightHandSideBoundedVector82;
const double crRightHandSideBoundedVector481 =             DN_DX_1_0*crRightHandSideBoundedVector53 + DN_DX_2_0*crRightHandSideBoundedVector54 + 0.16666666666666666*crRightHandSideBoundedVector65;
const double crRightHandSideBoundedVector482 =             0.16666666666666666*crRightHandSideBoundedVector157;
const double crRightHandSideBoundedVector483 =             crRightHandSideBoundedVector20*crRightHandSideBoundedVector270 - crRightHandSideBoundedVector271 + crRightHandSideBoundedVector481*crRightHandSideBoundedVector82 - crRightHandSideBoundedVector482*crRightHandSideBoundedVector82;
const double crRightHandSideBoundedVector484 =             crRightHandSideBoundedVector119*crRightHandSideBoundedVector475 - crRightHandSideBoundedVector188*crRightHandSideBoundedVector276;
const double crRightHandSideBoundedVector485 =             crRightHandSideBoundedVector306 + crRightHandSideBoundedVector484;
const double crRightHandSideBoundedVector486 =             -crRightHandSideBoundedVector119*crRightHandSideBoundedVector479 + crRightHandSideBoundedVector144*crRightHandSideBoundedVector303;
const double crRightHandSideBoundedVector487 =             crRightHandSideBoundedVector119*crRightHandSideBoundedVector481 - crRightHandSideBoundedVector119*crRightHandSideBoundedVector482 + crRightHandSideBoundedVector20*crRightHandSideBoundedVector277 - crRightHandSideBoundedVector278;
const double crRightHandSideBoundedVector488 =             crRightHandSideBoundedVector146*crRightHandSideBoundedVector259;
const double crRightHandSideBoundedVector489 =             3.0000000000000004*crRightHandSideBoundedVector17;
const double crRightHandSideBoundedVector490 =             crRightHandSideBoundedVector71*(crRightHandSideBoundedVector394 - crRightHandSideBoundedVector489);
const double crRightHandSideBoundedVector491 =             crRightHandSideBoundedVector160*crRightHandSideBoundedVector55;
const double crRightHandSideBoundedVector492 =             crRightHandSideBoundedVector398*crRightHandSideBoundedVector491;
const double crRightHandSideBoundedVector493 =             crRightHandSideBoundedVector290*crRightHandSideBoundedVector55;
const double crRightHandSideBoundedVector494 =             crRightHandSideBoundedVector262*crRightHandSideBoundedVector45;
const double crRightHandSideBoundedVector495 =             crRightHandSideBoundedVector494*crRightHandSideBoundedVector51;
const double crRightHandSideBoundedVector496 =             0.1111111111111111*f_ext(1,1);
const double crRightHandSideBoundedVector497 =             0.1111111111111111*f_ext(2,1);
const double crRightHandSideBoundedVector498 =             crRightHandSideBoundedVector496 + crRightHandSideBoundedVector497 + 0.44444444444444442*f_ext(0,1);
const double crRightHandSideBoundedVector499 =             0.16666666666666666*crRightHandSideBoundedVector71;
const double crRightHandSideBoundedVector500 =             3.0000000000000004*crRightHandSideBoundedVector90;
const double crRightHandSideBoundedVector501 =             crRightHandSideBoundedVector416*(crRightHandSideBoundedVector418 - crRightHandSideBoundedVector500);
const double crRightHandSideBoundedVector502 =             crRightHandSideBoundedVector422*crRightHandSideBoundedVector61;
const double crRightHandSideBoundedVector503 =             crRightHandSideBoundedVector110*crRightHandSideBoundedVector176;
const double crRightHandSideBoundedVector504 =             crRightHandSideBoundedVector424*crRightHandSideBoundedVector503;
const double crRightHandSideBoundedVector505 =             crRightHandSideBoundedVector289*crRightHandSideBoundedVector478;
const double crRightHandSideBoundedVector506 =             crRightHandSideBoundedVector157*crRightHandSideBoundedVector420;
const double crRightHandSideBoundedVector507 =             0.027777777777777776*f_ext(0,1);
const double crRightHandSideBoundedVector508 =             0.027777777777777776*f_ext(2,1);
const double crRightHandSideBoundedVector509 =             crRightHandSideBoundedVector496 + crRightHandSideBoundedVector507 + crRightHandSideBoundedVector508;
const double crRightHandSideBoundedVector510 =             crRightHandSideBoundedVector364*crRightHandSideBoundedVector505 - crRightHandSideBoundedVector420*crRightHandSideBoundedVector68 - crRightHandSideBoundedVector499*crRightHandSideBoundedVector501 - crRightHandSideBoundedVector502 + crRightHandSideBoundedVector504 + crRightHandSideBoundedVector506 + crRightHandSideBoundedVector509;
const double crRightHandSideBoundedVector511 =             3.0000000000000004*crRightHandSideBoundedVector123;
const double crRightHandSideBoundedVector512 =             crRightHandSideBoundedVector438*(crRightHandSideBoundedVector440 - crRightHandSideBoundedVector511);
const double crRightHandSideBoundedVector513 =             crRightHandSideBoundedVector444*crRightHandSideBoundedVector61;
const double crRightHandSideBoundedVector514 =             crRightHandSideBoundedVector139*crRightHandSideBoundedVector188;
const double crRightHandSideBoundedVector515 =             crRightHandSideBoundedVector446*crRightHandSideBoundedVector514;
const double crRightHandSideBoundedVector516 =             crRightHandSideBoundedVector157*crRightHandSideBoundedVector442;
const double crRightHandSideBoundedVector517 =             0.027777777777777776*f_ext(1,1);
const double crRightHandSideBoundedVector518 =             crRightHandSideBoundedVector497 + crRightHandSideBoundedVector507 + crRightHandSideBoundedVector517;
const double crRightHandSideBoundedVector519 =             crRightHandSideBoundedVector386*crRightHandSideBoundedVector505 - crRightHandSideBoundedVector442*crRightHandSideBoundedVector68 - crRightHandSideBoundedVector499*crRightHandSideBoundedVector512 - crRightHandSideBoundedVector513 + crRightHandSideBoundedVector515 + crRightHandSideBoundedVector516 + crRightHandSideBoundedVector518;
const double crRightHandSideBoundedVector520 =             crRightHandSideBoundedVector137*crRightHandSideBoundedVector482 + crRightHandSideBoundedVector140*crRightHandSideBoundedVector479 + 0.16666666666666666*crRightHandSideBoundedVector181 - 0.16666666666666666*crRightHandSideBoundedVector184 - 0.16666666666666666*crRightHandSideBoundedVector185 - 0.16666666666666666*crRightHandSideBoundedVector186 + crRightHandSideBoundedVector188*crRightHandSideBoundedVector444;
const double crRightHandSideBoundedVector521 =             crRightHandSideBoundedVector107*crRightHandSideBoundedVector482 + crRightHandSideBoundedVector111*crRightHandSideBoundedVector479 + 0.16666666666666666*crRightHandSideBoundedVector169 - 0.16666666666666666*crRightHandSideBoundedVector172 - 0.16666666666666666*crRightHandSideBoundedVector173 - 0.16666666666666666*crRightHandSideBoundedVector174 + crRightHandSideBoundedVector176*crRightHandSideBoundedVector422;
const double crRightHandSideBoundedVector522 =             -crRightHandSideBoundedVector319;
const double crRightHandSideBoundedVector523 =             crRightHandSideBoundedVector317 + crRightHandSideBoundedVector522;
const double crRightHandSideBoundedVector524 =             crRightHandSideBoundedVector320 + crRightHandSideBoundedVector522;
const double crRightHandSideBoundedVector525 =             crRightHandSideBoundedVector56*mu;
const double crRightHandSideBoundedVector526 =             2.2500000000000004*crRightHandSideBoundedVector238 + 2.2500000000000004*crRightHandSideBoundedVector239 - 2.2500000000000004*crRightHandSideBoundedVector240 - 2.2500000000000004*crRightHandSideBoundedVector241;
const double crRightHandSideBoundedVector527 =             crRightHandSideBoundedVector46*mu;
const double crRightHandSideBoundedVector528 =             2.2500000000000004*crRightHandSideBoundedVector3;
const double crRightHandSideBoundedVector529 =             2.2500000000000004*crRightHandSideBoundedVector24;
const double crRightHandSideBoundedVector530 =             1.5000000000000002*crRightHandSideBoundedVector333;
const double crRightHandSideBoundedVector531 =             crRightHandSideBoundedVector32*crRightHandSideBoundedVector530;
const double crRightHandSideBoundedVector532 =             crRightHandSideBoundedVector14*crRightHandSideBoundedVector531 - crRightHandSideBoundedVector159*crRightHandSideBoundedVector68 + crRightHandSideBoundedVector17*crRightHandSideBoundedVector531 - crRightHandSideBoundedVector32*crRightHandSideBoundedVector529 + crRightHandSideBoundedVector528*crRightHandSideBoundedVector75 - crRightHandSideBoundedVector61*crRightHandSideBoundedVector73;
const double crRightHandSideBoundedVector533 =             gamma*k_st;
const double crRightHandSideBoundedVector534 =             crRightHandSideBoundedVector530*crRightHandSideBoundedVector71;
const double crRightHandSideBoundedVector535 =             crRightHandSideBoundedVector14*crRightHandSideBoundedVector534 - crRightHandSideBoundedVector159*crRightHandSideBoundedVector41 + crRightHandSideBoundedVector162*crRightHandSideBoundedVector528 + crRightHandSideBoundedVector17*crRightHandSideBoundedVector534 - crRightHandSideBoundedVector340 - crRightHandSideBoundedVector529*crRightHandSideBoundedVector71;
const double crRightHandSideBoundedVector536 =             crRightHandSideBoundedVector3*crRightHandSideBoundedVector535;
const double crRightHandSideBoundedVector537 =             crRightHandSideBoundedVector202*crRightHandSideBoundedVector536;
const double crRightHandSideBoundedVector538 =             crRightHandSideBoundedVector3*crRightHandSideBoundedVector532;
const double crRightHandSideBoundedVector539 =             crRightHandSideBoundedVector206*crRightHandSideBoundedVector533;
const double crRightHandSideBoundedVector540 =             gamma*k_sc;
const double crRightHandSideBoundedVector541 =             crRightHandSideBoundedVector209*crRightHandSideBoundedVector540;
const double crRightHandSideBoundedVector542 =             crRightHandSideBoundedVector11*(crRightHandSideBoundedVector314*crRightHandSideBoundedVector532 + crRightHandSideBoundedVector525*crRightHandSideBoundedVector526 + crRightHandSideBoundedVector527*(-3.0000000000000009*crRightHandSideBoundedVector160 + 3.0000000000000009*crRightHandSideBoundedVector192 - 1.5000000000000002*crRightHandSideBoundedVector196 + crRightHandSideBoundedVector285) + crRightHandSideBoundedVector533*crRightHandSideBoundedVector537 - crRightHandSideBoundedVector537*crRightHandSideBoundedVector540 + crRightHandSideBoundedVector538*crRightHandSideBoundedVector539 + crRightHandSideBoundedVector538*crRightHandSideBoundedVector541);
const double crRightHandSideBoundedVector543 =             crRightHandSideBoundedVector111*mu;
const double crRightHandSideBoundedVector544 =             2.2500000000000004*crRightHandSideBoundedVector244 + 2.2500000000000004*crRightHandSideBoundedVector245 - 2.2500000000000004*crRightHandSideBoundedVector246 - 2.2500000000000004*crRightHandSideBoundedVector247;
const double crRightHandSideBoundedVector545 =             1.5000000000000002*crRightHandSideBoundedVector114;
const double crRightHandSideBoundedVector546 =             crRightHandSideBoundedVector107*mu;
const double crRightHandSideBoundedVector547 =             2.2500000000000004*crRightHandSideBoundedVector81;
const double crRightHandSideBoundedVector548 =             2.2500000000000004*crRightHandSideBoundedVector94;
const double crRightHandSideBoundedVector549 =             1.5000000000000002*crRightHandSideBoundedVector358;
const double crRightHandSideBoundedVector550 =             crRightHandSideBoundedVector32*crRightHandSideBoundedVector549;
const double crRightHandSideBoundedVector551 =             -crRightHandSideBoundedVector113*crRightHandSideBoundedVector61 - crRightHandSideBoundedVector175*crRightHandSideBoundedVector68 - crRightHandSideBoundedVector32*crRightHandSideBoundedVector548 + crRightHandSideBoundedVector547*crRightHandSideBoundedVector75 + crRightHandSideBoundedVector550*crRightHandSideBoundedVector88 + crRightHandSideBoundedVector550*crRightHandSideBoundedVector90;
const double crRightHandSideBoundedVector552 =             crRightHandSideBoundedVector549*crRightHandSideBoundedVector71;
const double crRightHandSideBoundedVector553 =             crRightHandSideBoundedVector162*crRightHandSideBoundedVector547 - crRightHandSideBoundedVector175*crRightHandSideBoundedVector41 - crRightHandSideBoundedVector366 - crRightHandSideBoundedVector548*crRightHandSideBoundedVector71 + crRightHandSideBoundedVector552*crRightHandSideBoundedVector88 + crRightHandSideBoundedVector552*crRightHandSideBoundedVector90;
const double crRightHandSideBoundedVector554 =             crRightHandSideBoundedVector553*crRightHandSideBoundedVector81;
const double crRightHandSideBoundedVector555 =             crRightHandSideBoundedVector202*crRightHandSideBoundedVector554;
const double crRightHandSideBoundedVector556 =             crRightHandSideBoundedVector551*crRightHandSideBoundedVector81;
const double crRightHandSideBoundedVector557 =             crRightHandSideBoundedVector86*(crRightHandSideBoundedVector314*crRightHandSideBoundedVector551 + crRightHandSideBoundedVector533*crRightHandSideBoundedVector555 + crRightHandSideBoundedVector539*crRightHandSideBoundedVector556 - crRightHandSideBoundedVector540*crRightHandSideBoundedVector555 + crRightHandSideBoundedVector541*crRightHandSideBoundedVector556 + crRightHandSideBoundedVector543*crRightHandSideBoundedVector544 + crRightHandSideBoundedVector546*(-3.0000000000000009*crRightHandSideBoundedVector176 + 3.0000000000000009*crRightHandSideBoundedVector212 - 1.5000000000000002*crRightHandSideBoundedVector216 + crRightHandSideBoundedVector545));
const double crRightHandSideBoundedVector558 =             crRightHandSideBoundedVector140*mu;
const double crRightHandSideBoundedVector559 =             2.2500000000000004*crRightHandSideBoundedVector250 + 2.2500000000000004*crRightHandSideBoundedVector251 - 2.2500000000000004*crRightHandSideBoundedVector252 - 2.2500000000000004*crRightHandSideBoundedVector253;
const double crRightHandSideBoundedVector560 =             1.5000000000000002*crRightHandSideBoundedVector143;
const double crRightHandSideBoundedVector561 =             crRightHandSideBoundedVector137*mu;
const double crRightHandSideBoundedVector562 =             2.2500000000000004*crRightHandSideBoundedVector118;
const double crRightHandSideBoundedVector563 =             2.2500000000000004*crRightHandSideBoundedVector126;
const double crRightHandSideBoundedVector564 =             1.5000000000000002*crRightHandSideBoundedVector380;
const double crRightHandSideBoundedVector565 =             crRightHandSideBoundedVector32*crRightHandSideBoundedVector564;
const double crRightHandSideBoundedVector566 =             crRightHandSideBoundedVector122*crRightHandSideBoundedVector565 + crRightHandSideBoundedVector123*crRightHandSideBoundedVector565 - crRightHandSideBoundedVector142*crRightHandSideBoundedVector61 - crRightHandSideBoundedVector187*crRightHandSideBoundedVector68 - crRightHandSideBoundedVector32*crRightHandSideBoundedVector563 + crRightHandSideBoundedVector562*crRightHandSideBoundedVector75;
const double crRightHandSideBoundedVector567 =             crRightHandSideBoundedVector564*crRightHandSideBoundedVector71;
const double crRightHandSideBoundedVector568 =             crRightHandSideBoundedVector122*crRightHandSideBoundedVector567 + crRightHandSideBoundedVector123*crRightHandSideBoundedVector567 + crRightHandSideBoundedVector162*crRightHandSideBoundedVector562 - crRightHandSideBoundedVector187*crRightHandSideBoundedVector41 - crRightHandSideBoundedVector388 - crRightHandSideBoundedVector563*crRightHandSideBoundedVector71;
const double crRightHandSideBoundedVector569 =             crRightHandSideBoundedVector118*crRightHandSideBoundedVector568;
const double crRightHandSideBoundedVector570 =             crRightHandSideBoundedVector202*crRightHandSideBoundedVector569;
const double crRightHandSideBoundedVector571 =             crRightHandSideBoundedVector118*crRightHandSideBoundedVector566;
const double crRightHandSideBoundedVector572 =             crRightHandSideBoundedVector121*(crRightHandSideBoundedVector314*crRightHandSideBoundedVector566 + crRightHandSideBoundedVector533*crRightHandSideBoundedVector570 + crRightHandSideBoundedVector539*crRightHandSideBoundedVector571 - crRightHandSideBoundedVector540*crRightHandSideBoundedVector570 + crRightHandSideBoundedVector541*crRightHandSideBoundedVector571 + crRightHandSideBoundedVector558*crRightHandSideBoundedVector559 + crRightHandSideBoundedVector561*(-3.0000000000000009*crRightHandSideBoundedVector188 + 3.0000000000000009*crRightHandSideBoundedVector224 - 1.5000000000000002*crRightHandSideBoundedVector228 + crRightHandSideBoundedVector560));
const double crRightHandSideBoundedVector573 =             crRightHandSideBoundedVector202*crRightHandSideBoundedVector538;
const double crRightHandSideBoundedVector574 =             crRightHandSideBoundedVector457*crRightHandSideBoundedVector533;
const double crRightHandSideBoundedVector575 =             crRightHandSideBoundedVector459*crRightHandSideBoundedVector540;
const double crRightHandSideBoundedVector576 =             crRightHandSideBoundedVector11*(crRightHandSideBoundedVector314*crRightHandSideBoundedVector535 + crRightHandSideBoundedVector525*(-1.5000000000000002*crRightHandSideBoundedVector192 + 3.0000000000000009*crRightHandSideBoundedVector196 + crRightHandSideBoundedVector470 - 3.0000000000000009*crRightHandSideBoundedVector72) + crRightHandSideBoundedVector526*crRightHandSideBoundedVector527 + crRightHandSideBoundedVector533*crRightHandSideBoundedVector573 + crRightHandSideBoundedVector536*crRightHandSideBoundedVector574 + crRightHandSideBoundedVector536*crRightHandSideBoundedVector575 - crRightHandSideBoundedVector540*crRightHandSideBoundedVector573);
const double crRightHandSideBoundedVector577 =             1.5000000000000002*crRightHandSideBoundedVector176;
const double crRightHandSideBoundedVector578 =             crRightHandSideBoundedVector202*crRightHandSideBoundedVector556;
const double crRightHandSideBoundedVector579 =             crRightHandSideBoundedVector86*(crRightHandSideBoundedVector314*crRightHandSideBoundedVector553 + crRightHandSideBoundedVector533*crRightHandSideBoundedVector578 - crRightHandSideBoundedVector540*crRightHandSideBoundedVector578 + crRightHandSideBoundedVector543*(-3.0000000000000009*crRightHandSideBoundedVector114 - 1.5000000000000002*crRightHandSideBoundedVector212 + 3.0000000000000009*crRightHandSideBoundedVector216 + crRightHandSideBoundedVector577) + crRightHandSideBoundedVector544*crRightHandSideBoundedVector546 + crRightHandSideBoundedVector554*crRightHandSideBoundedVector574 + crRightHandSideBoundedVector554*crRightHandSideBoundedVector575);
const double crRightHandSideBoundedVector580 =             1.5000000000000002*crRightHandSideBoundedVector188;
const double crRightHandSideBoundedVector581 =             crRightHandSideBoundedVector202*crRightHandSideBoundedVector571;
const double crRightHandSideBoundedVector582 =             crRightHandSideBoundedVector121*(crRightHandSideBoundedVector314*crRightHandSideBoundedVector568 + crRightHandSideBoundedVector533*crRightHandSideBoundedVector581 - crRightHandSideBoundedVector540*crRightHandSideBoundedVector581 + crRightHandSideBoundedVector558*(-3.0000000000000009*crRightHandSideBoundedVector143 - 1.5000000000000002*crRightHandSideBoundedVector224 + 3.0000000000000009*crRightHandSideBoundedVector228 + crRightHandSideBoundedVector580) + crRightHandSideBoundedVector559*crRightHandSideBoundedVector561 + crRightHandSideBoundedVector569*crRightHandSideBoundedVector574 + crRightHandSideBoundedVector569*crRightHandSideBoundedVector575);
const double crRightHandSideBoundedVector583 =             crRightHandSideBoundedVector20*crRightHandSideBoundedVector333;
const double crRightHandSideBoundedVector584 =             crRightHandSideBoundedVector335 + crRightHandSideBoundedVector393*crRightHandSideBoundedVector583;
const double crRightHandSideBoundedVector585 =             0.66666666666666663*crRightHandSideBoundedVector11;
const double crRightHandSideBoundedVector586 =             4.5000000000000009*crRightHandSideBoundedVector316;
const double crRightHandSideBoundedVector587 =             crRightHandSideBoundedVector45*crRightHandSideBoundedVector61;
const double crRightHandSideBoundedVector588 =             crRightHandSideBoundedVector402*crRightHandSideBoundedVector51;
const double crRightHandSideBoundedVector589 =             crRightHandSideBoundedVector335 + crRightHandSideBoundedVector489*crRightHandSideBoundedVector583;
const double crRightHandSideBoundedVector590 =             crRightHandSideBoundedVector41*crRightHandSideBoundedVector55;
const double crRightHandSideBoundedVector591 =             crRightHandSideBoundedVector341*crRightHandSideBoundedVector73;
const double crRightHandSideBoundedVector592 =             -crRightHandSideBoundedVector349;
const double crRightHandSideBoundedVector593 =             crRightHandSideBoundedVector347 + crRightHandSideBoundedVector592;
const double crRightHandSideBoundedVector594 =             0.16666666666666666*crRightHandSideBoundedVector82;
const double crRightHandSideBoundedVector595 =             crRightHandSideBoundedVector20*crRightHandSideBoundedVector358;
const double crRightHandSideBoundedVector596 =             crRightHandSideBoundedVector86*(crRightHandSideBoundedVector360 + crRightHandSideBoundedVector417*crRightHandSideBoundedVector595);
const double crRightHandSideBoundedVector597 =             1.1250000000000002*crRightHandSideBoundedVector346;
const double crRightHandSideBoundedVector598 =             crRightHandSideBoundedVector106*crRightHandSideBoundedVector61;
const double crRightHandSideBoundedVector599 =             -crRightHandSideBoundedVector20*crRightHandSideBoundedVector421 - crRightHandSideBoundedVector20*crRightHandSideBoundedVector423 + crRightHandSideBoundedVector20*crRightHandSideBoundedVector426 + crRightHandSideBoundedVector326*crRightHandSideBoundedVector594 + crRightHandSideBoundedVector415*crRightHandSideBoundedVector596 - crRightHandSideBoundedVector429 + crRightHandSideBoundedVector432 - crRightHandSideBoundedVector597*crRightHandSideBoundedVector598;
const double crRightHandSideBoundedVector600 =             crRightHandSideBoundedVector350 + crRightHandSideBoundedVector592;
const double crRightHandSideBoundedVector601 =             crRightHandSideBoundedVector113*crRightHandSideBoundedVector367;
const double crRightHandSideBoundedVector602 =             crRightHandSideBoundedVector86*(crRightHandSideBoundedVector360 + crRightHandSideBoundedVector500*crRightHandSideBoundedVector595);
const double crRightHandSideBoundedVector603 =             crRightHandSideBoundedVector110*crRightHandSideBoundedVector41;
const double crRightHandSideBoundedVector604 =             crRightHandSideBoundedVector420*crRightHandSideBoundedVector69;
const double crRightHandSideBoundedVector605 =             -crRightHandSideBoundedVector20*crRightHandSideBoundedVector502 + crRightHandSideBoundedVector20*crRightHandSideBoundedVector504 + crRightHandSideBoundedVector328*crRightHandSideBoundedVector594 + crRightHandSideBoundedVector499*crRightHandSideBoundedVector602 - crRightHandSideBoundedVector506 + crRightHandSideBoundedVector509 - crRightHandSideBoundedVector597*crRightHandSideBoundedVector603 - crRightHandSideBoundedVector604;
const double crRightHandSideBoundedVector606 =             -crRightHandSideBoundedVector373;
const double crRightHandSideBoundedVector607 =             crRightHandSideBoundedVector371 + crRightHandSideBoundedVector606;
const double crRightHandSideBoundedVector608 =             0.16666666666666666*crRightHandSideBoundedVector119;
const double crRightHandSideBoundedVector609 =             crRightHandSideBoundedVector20*crRightHandSideBoundedVector380;
const double crRightHandSideBoundedVector610 =             crRightHandSideBoundedVector121*(crRightHandSideBoundedVector382 + crRightHandSideBoundedVector439*crRightHandSideBoundedVector609);
const double crRightHandSideBoundedVector611 =             1.1250000000000002*crRightHandSideBoundedVector370;
const double crRightHandSideBoundedVector612 =             crRightHandSideBoundedVector136*crRightHandSideBoundedVector61;
const double crRightHandSideBoundedVector613 =             -crRightHandSideBoundedVector20*crRightHandSideBoundedVector443 - crRightHandSideBoundedVector20*crRightHandSideBoundedVector445 + crRightHandSideBoundedVector20*crRightHandSideBoundedVector448 + crRightHandSideBoundedVector326*crRightHandSideBoundedVector608 + crRightHandSideBoundedVector415*crRightHandSideBoundedVector610 - crRightHandSideBoundedVector450 + crRightHandSideBoundedVector452 - crRightHandSideBoundedVector611*crRightHandSideBoundedVector612;
const double crRightHandSideBoundedVector614 =             crRightHandSideBoundedVector374 + crRightHandSideBoundedVector606;
const double crRightHandSideBoundedVector615 =             crRightHandSideBoundedVector142*crRightHandSideBoundedVector389;
const double crRightHandSideBoundedVector616 =             crRightHandSideBoundedVector121*(crRightHandSideBoundedVector382 + crRightHandSideBoundedVector511*crRightHandSideBoundedVector609);
const double crRightHandSideBoundedVector617 =             crRightHandSideBoundedVector139*crRightHandSideBoundedVector41;
const double crRightHandSideBoundedVector618 =             crRightHandSideBoundedVector442*crRightHandSideBoundedVector69;
const double crRightHandSideBoundedVector619 =             -crRightHandSideBoundedVector20*crRightHandSideBoundedVector513 + crRightHandSideBoundedVector20*crRightHandSideBoundedVector515 + crRightHandSideBoundedVector328*crRightHandSideBoundedVector608 + crRightHandSideBoundedVector499*crRightHandSideBoundedVector616 - crRightHandSideBoundedVector516 + crRightHandSideBoundedVector518 - crRightHandSideBoundedVector611*crRightHandSideBoundedVector617 - crRightHandSideBoundedVector618;
const double crRightHandSideBoundedVector620 =             crRightHandSideBoundedVector62*h;
const double crRightHandSideBoundedVector621 =             crRightHandSideBoundedVector344*crRightHandSideBoundedVector620;
const double crRightHandSideBoundedVector622 =             crRightHandSideBoundedVector369*crRightHandSideBoundedVector620;
const double crRightHandSideBoundedVector623 =             crRightHandSideBoundedVector390*crRightHandSideBoundedVector620;
const double crRightHandSideBoundedVector624 =             0.1111111111111111*r[1];
const double crRightHandSideBoundedVector625 =             0.1111111111111111*r[2];
const double crRightHandSideBoundedVector626 =             crRightHandSideBoundedVector45*crRightHandSideBoundedVector55;
const double crRightHandSideBoundedVector627 =             crRightHandSideBoundedVector398*crRightHandSideBoundedVector626;
const double crRightHandSideBoundedVector628 =             crRightHandSideBoundedVector336*crRightHandSideBoundedVector45;
const double crRightHandSideBoundedVector629 =             -1.125*U_1_3;
const double crRightHandSideBoundedVector630 =             -1.125*U_2_3;
const double crRightHandSideBoundedVector631 =             4.5000000000000009*crRightHandSideBoundedVector20;
const double crRightHandSideBoundedVector632 =             -4.5*U_0_3 - 6.7500000000000009*crRightHandSideBoundedVector318 + crRightHandSideBoundedVector334*crRightHandSideBoundedVector631 + crRightHandSideBoundedVector629 + crRightHandSideBoundedVector630;
const double crRightHandSideBoundedVector633 =             crRightHandSideBoundedVector391*crRightHandSideBoundedVector632;
const double crRightHandSideBoundedVector634 =             crRightHandSideBoundedVector11*crRightHandSideBoundedVector584;
const double crRightHandSideBoundedVector635 =             crRightHandSideBoundedVector11*crRightHandSideBoundedVector589;
const double crRightHandSideBoundedVector636 =             crRightHandSideBoundedVector106*crRightHandSideBoundedVector361;
const double crRightHandSideBoundedVector637 =             0.027777777777777776*r[0];
const double crRightHandSideBoundedVector638 =             0.027777777777777776*r[2];
const double crRightHandSideBoundedVector639 =             -1.125*U_0_3;
const double crRightHandSideBoundedVector640 =             crRightHandSideBoundedVector416*(-4.5*U_1_3 - 6.7500000000000009*crRightHandSideBoundedVector348 + crRightHandSideBoundedVector359*crRightHandSideBoundedVector631 + crRightHandSideBoundedVector630 + crRightHandSideBoundedVector639);
const double crRightHandSideBoundedVector641 =             0.16666666666666666*crRightHandSideBoundedVector640;
const double crRightHandSideBoundedVector642 =             crRightHandSideBoundedVector106*crRightHandSideBoundedVector110;
const double crRightHandSideBoundedVector643 =             crRightHandSideBoundedVector424*crRightHandSideBoundedVector642;
const double crRightHandSideBoundedVector644 =             -crRightHandSideBoundedVector114*crRightHandSideBoundedVector641 + crRightHandSideBoundedVector157*crRightHandSideBoundedVector643 - crRightHandSideBoundedVector176*crRightHandSideBoundedVector641 + crRightHandSideBoundedVector301*crRightHandSideBoundedVector596 - crRightHandSideBoundedVector326*crRightHandSideBoundedVector420 - crRightHandSideBoundedVector328*crRightHandSideBoundedVector422 + crRightHandSideBoundedVector478*crRightHandSideBoundedVector602 + crRightHandSideBoundedVector624 + crRightHandSideBoundedVector637 + crRightHandSideBoundedVector638 + crRightHandSideBoundedVector643*crRightHandSideBoundedVector69;
const double crRightHandSideBoundedVector645 =             crRightHandSideBoundedVector136*crRightHandSideBoundedVector383;
const double crRightHandSideBoundedVector646 =             0.027777777777777776*r[1];
const double crRightHandSideBoundedVector647 =             crRightHandSideBoundedVector438*(-4.5*U_2_3 - 6.7500000000000009*crRightHandSideBoundedVector372 + crRightHandSideBoundedVector381*crRightHandSideBoundedVector631 + crRightHandSideBoundedVector629 + crRightHandSideBoundedVector639);
const double crRightHandSideBoundedVector648 =             0.16666666666666666*crRightHandSideBoundedVector647;
const double crRightHandSideBoundedVector649 =             crRightHandSideBoundedVector136*crRightHandSideBoundedVector139;
const double crRightHandSideBoundedVector650 =             crRightHandSideBoundedVector446*crRightHandSideBoundedVector649;
const double crRightHandSideBoundedVector651 =             -crRightHandSideBoundedVector143*crRightHandSideBoundedVector648 + crRightHandSideBoundedVector157*crRightHandSideBoundedVector650 - crRightHandSideBoundedVector188*crRightHandSideBoundedVector648 + crRightHandSideBoundedVector301*crRightHandSideBoundedVector610 - crRightHandSideBoundedVector326*crRightHandSideBoundedVector442 - crRightHandSideBoundedVector328*crRightHandSideBoundedVector444 + crRightHandSideBoundedVector478*crRightHandSideBoundedVector616 + crRightHandSideBoundedVector625 + crRightHandSideBoundedVector637 + crRightHandSideBoundedVector646 + crRightHandSideBoundedVector650*crRightHandSideBoundedVector69;
const double crRightHandSideBoundedVector652 =             crRightHandSideBoundedVector139*crRightHandSideBoundedVector516 + crRightHandSideBoundedVector139*crRightHandSideBoundedVector618 + crRightHandSideBoundedVector301*crRightHandSideBoundedVector607 + 0.16666666666666666*crRightHandSideBoundedVector375 + 0.16666666666666666*crRightHandSideBoundedVector376 + 0.16666666666666666*crRightHandSideBoundedVector377 - 0.16666666666666666*crRightHandSideBoundedVector378 - 0.16666666666666666*crRightHandSideBoundedVector379 - 0.16666666666666666*crRightHandSideBoundedVector384 - 0.16666666666666666*crRightHandSideBoundedVector385 + crRightHandSideBoundedVector478*crRightHandSideBoundedVector614;
const double crRightHandSideBoundedVector653 =             crRightHandSideBoundedVector110*crRightHandSideBoundedVector506 + crRightHandSideBoundedVector110*crRightHandSideBoundedVector604 + crRightHandSideBoundedVector301*crRightHandSideBoundedVector593 + 0.16666666666666666*crRightHandSideBoundedVector352 + 0.16666666666666666*crRightHandSideBoundedVector353 + 0.16666666666666666*crRightHandSideBoundedVector354 - 0.16666666666666666*crRightHandSideBoundedVector355 - 0.16666666666666666*crRightHandSideBoundedVector356 - 0.16666666666666666*crRightHandSideBoundedVector362 - 0.16666666666666666*crRightHandSideBoundedVector363 + crRightHandSideBoundedVector478*crRightHandSideBoundedVector600;
const double crRightHandSideBoundedVector654 =             0.99999999999999989*crRightHandSideBoundedVector38 + 0.99999999999999989*crRightHandSideBoundedVector39 + 0.99999999999999989*crRightHandSideBoundedVector40 + 0.99999999999999989*crRightHandSideBoundedVector58 + 0.99999999999999989*crRightHandSideBoundedVector59 + 0.99999999999999989*crRightHandSideBoundedVector60;
const double crRightHandSideBoundedVector655 =             DN_DX_1_1*crRightHandSideBoundedVector46;
const double crRightHandSideBoundedVector656 =             DN_DX_1_0*crRightHandSideBoundedVector56;
const double crRightHandSideBoundedVector657 =             0.37500000000000006*crRightHandSideBoundedVector11;
const double crRightHandSideBoundedVector658 =             0.37500000000000006*crRightHandSideBoundedVector316;
const double crRightHandSideBoundedVector659 =             crRightHandSideBoundedVector240*crRightHandSideBoundedVector657 - crRightHandSideBoundedVector241*crRightHandSideBoundedVector658 - crRightHandSideBoundedVector267*crRightHandSideBoundedVector4 + crRightHandSideBoundedVector268*crRightHandSideBoundedVector4;
const double crRightHandSideBoundedVector660 =             DN_DX_1_1*crRightHandSideBoundedVector107;
const double crRightHandSideBoundedVector661 =             0.66666666666666663*crRightHandSideBoundedVector82;
const double crRightHandSideBoundedVector662 =             DN_DX_1_0*crRightHandSideBoundedVector111;
const double crRightHandSideBoundedVector663 =             1.5000000000000002*crRightHandSideBoundedVector86;
const double crRightHandSideBoundedVector664 =             crRightHandSideBoundedVector261*crRightHandSideBoundedVector86;
const double crRightHandSideBoundedVector665 =             DN_DX_1_1*crRightHandSideBoundedVector137;
const double crRightHandSideBoundedVector666 =             DN_DX_1_0*crRightHandSideBoundedVector140;
const double crRightHandSideBoundedVector667 =             crRightHandSideBoundedVector545*crRightHandSideBoundedVector86;
const double crRightHandSideBoundedVector668 =             DN_DX_1_0*crRightHandSideBoundedVector46;
const double crRightHandSideBoundedVector669 =             DN_DX_1_1*crRightHandSideBoundedVector56;
const double crRightHandSideBoundedVector670 =             crRightHandSideBoundedVector657*crRightHandSideBoundedVector72;
const double crRightHandSideBoundedVector671 =             crRightHandSideBoundedVector297*crRightHandSideBoundedVector4 - crRightHandSideBoundedVector670;
const double crRightHandSideBoundedVector672 =             crRightHandSideBoundedVector669 + crRightHandSideBoundedVector671;
const double crRightHandSideBoundedVector673 =             crRightHandSideBoundedVector161*crRightHandSideBoundedVector303 - crRightHandSideBoundedVector287*crRightHandSideBoundedVector301;
const double crRightHandSideBoundedVector674 =             crRightHandSideBoundedVector282*crRightHandSideBoundedVector82;
const double crRightHandSideBoundedVector675 =             DN_DX_1_0*crRightHandSideBoundedVector107;
const double crRightHandSideBoundedVector676 =             crRightHandSideBoundedVector289*crRightHandSideBoundedVector86;
const double crRightHandSideBoundedVector677 =             DN_DX_1_1*crRightHandSideBoundedVector111;
const double crRightHandSideBoundedVector678 =             crRightHandSideBoundedVector291*crRightHandSideBoundedVector82 - crRightHandSideBoundedVector667 + crRightHandSideBoundedVector677;
const double crRightHandSideBoundedVector679 =             DN_DX_1_0*crRightHandSideBoundedVector137;
const double crRightHandSideBoundedVector680 =             DN_DX_1_1*crRightHandSideBoundedVector140;
const double crRightHandSideBoundedVector681 =             crRightHandSideBoundedVector309 + crRightHandSideBoundedVector680;
const double crRightHandSideBoundedVector682 =             crRightHandSideBoundedVector259*h;
const double crRightHandSideBoundedVector683 =             DN_DX_1_0*crRightHandSideBoundedVector682;
const double crRightHandSideBoundedVector684 =             0.16666666666666666*crRightHandSideBoundedVector391;
const double crRightHandSideBoundedVector685 =             crRightHandSideBoundedVector45*crRightHandSideBoundedVector657;
const double crRightHandSideBoundedVector686 =             crRightHandSideBoundedVector55*crRightHandSideBoundedVector657;
const double crRightHandSideBoundedVector687 =             1.125*crRightHandSideBoundedVector391;
const double crRightHandSideBoundedVector688 =             crRightHandSideBoundedVector399*crRightHandSideBoundedVector687;
const double crRightHandSideBoundedVector689 =             crRightHandSideBoundedVector686*crRightHandSideBoundedVector69;
const double crRightHandSideBoundedVector690 =             0.1111111111111111*f_ext(0,0);
const double crRightHandSideBoundedVector691 =             crRightHandSideBoundedVector431 + crRightHandSideBoundedVector451 + crRightHandSideBoundedVector690;
const double crRightHandSideBoundedVector692 =             crRightHandSideBoundedVector301*crRightHandSideBoundedVector401 - crRightHandSideBoundedVector395*crRightHandSideBoundedVector684 - crRightHandSideBoundedVector41*crRightHandSideBoundedVector685 - crRightHandSideBoundedVector51*crRightHandSideBoundedVector686 + crRightHandSideBoundedVector688 + crRightHandSideBoundedVector689 + crRightHandSideBoundedVector691;
const double crRightHandSideBoundedVector693 =             0.66666666666666663*crRightHandSideBoundedVector32;
const double crRightHandSideBoundedVector694 =             crRightHandSideBoundedVector106*crRightHandSideBoundedVector663;
const double crRightHandSideBoundedVector695 =             crRightHandSideBoundedVector110*crRightHandSideBoundedVector663;
const double crRightHandSideBoundedVector696 =             4.5*crRightHandSideBoundedVector416;
const double crRightHandSideBoundedVector697 =             crRightHandSideBoundedVector425*crRightHandSideBoundedVector696;
const double crRightHandSideBoundedVector698 =             crRightHandSideBoundedVector283*crRightHandSideBoundedVector289;
const double crRightHandSideBoundedVector699 =             crRightHandSideBoundedVector261*crRightHandSideBoundedVector364;
const double crRightHandSideBoundedVector700 =             crRightHandSideBoundedVector68*crRightHandSideBoundedVector699;
const double crRightHandSideBoundedVector701 =             crRightHandSideBoundedVector407 + crRightHandSideBoundedVector690 + 0.44444444444444442*f_ext(1,0);
const double crRightHandSideBoundedVector702 =             crRightHandSideBoundedVector268*crRightHandSideBoundedVector56 + crRightHandSideBoundedVector302*crRightHandSideBoundedVector46 + 0.16666666666666666*crRightHandSideBoundedVector31 - 0.16666666666666666*crRightHandSideBoundedVector37 + crRightHandSideBoundedVector45*crRightHandSideBoundedVector670 - 0.16666666666666666*crRightHandSideBoundedVector47 - 0.16666666666666666*crRightHandSideBoundedVector57 - 0.99999999999999989*crRightHandSideBoundedVector76;
const double crRightHandSideBoundedVector703 =             crRightHandSideBoundedVector577*crRightHandSideBoundedVector86;
const double crRightHandSideBoundedVector704 =             -crRightHandSideBoundedVector160*crRightHandSideBoundedVector657 + crRightHandSideBoundedVector4*crRightHandSideBoundedVector475;
const double crRightHandSideBoundedVector705 =             crRightHandSideBoundedVector668 + crRightHandSideBoundedVector704;
const double crRightHandSideBoundedVector706 =             -crRightHandSideBoundedVector287*crRightHandSideBoundedVector478 + 0.16666666666666666*crRightHandSideBoundedVector290*crRightHandSideBoundedVector72;
const double crRightHandSideBoundedVector707 =             crRightHandSideBoundedVector240*crRightHandSideBoundedVector658 - crRightHandSideBoundedVector241*crRightHandSideBoundedVector657 + crRightHandSideBoundedVector4*crRightHandSideBoundedVector481 - crRightHandSideBoundedVector4*crRightHandSideBoundedVector482;
const double crRightHandSideBoundedVector708 =             crRightHandSideBoundedVector472*crRightHandSideBoundedVector82 + crRightHandSideBoundedVector675 - crRightHandSideBoundedVector703;
const double crRightHandSideBoundedVector709 =             crRightHandSideBoundedVector484 + crRightHandSideBoundedVector679;
const double crRightHandSideBoundedVector710 =             DN_DX_1_1*crRightHandSideBoundedVector682;
const double crRightHandSideBoundedVector711 =             crRightHandSideBoundedVector491*crRightHandSideBoundedVector687;
const double crRightHandSideBoundedVector712 =             crRightHandSideBoundedVector45*crRightHandSideBoundedVector658;
const double crRightHandSideBoundedVector713 =             crRightHandSideBoundedVector51*crRightHandSideBoundedVector712;
const double crRightHandSideBoundedVector714 =             0.1111111111111111*f_ext(0,1);
const double crRightHandSideBoundedVector715 =             crRightHandSideBoundedVector508 + crRightHandSideBoundedVector517 + crRightHandSideBoundedVector714;
const double crRightHandSideBoundedVector716 =             crRightHandSideBoundedVector478*crRightHandSideBoundedVector493 - crRightHandSideBoundedVector490*crRightHandSideBoundedVector684 - crRightHandSideBoundedVector61*crRightHandSideBoundedVector686 - crRightHandSideBoundedVector68*crRightHandSideBoundedVector685 + crRightHandSideBoundedVector711 + crRightHandSideBoundedVector713 + crRightHandSideBoundedVector715;
const double crRightHandSideBoundedVector717 =             0.66666666666666663*crRightHandSideBoundedVector71;
const double crRightHandSideBoundedVector718 =             crRightHandSideBoundedVector503*crRightHandSideBoundedVector696;
const double crRightHandSideBoundedVector719 =             crRightHandSideBoundedVector289*crRightHandSideBoundedVector468;
const double crRightHandSideBoundedVector720 =             crRightHandSideBoundedVector261*crRightHandSideBoundedVector428;
const double crRightHandSideBoundedVector721 =             crRightHandSideBoundedVector51*crRightHandSideBoundedVector720;
const double crRightHandSideBoundedVector722 =             crRightHandSideBoundedVector497 + crRightHandSideBoundedVector714 + 0.44444444444444442*f_ext(1,1);
const double crRightHandSideBoundedVector723 =             0.16666666666666666*crRightHandSideBoundedVector150 - 0.16666666666666666*crRightHandSideBoundedVector153 - 0.16666666666666666*crRightHandSideBoundedVector154 - 0.16666666666666666*crRightHandSideBoundedVector155 + crRightHandSideBoundedVector160*crRightHandSideBoundedVector686 - 0.99999999999999989*crRightHandSideBoundedVector163 + crRightHandSideBoundedVector46*crRightHandSideBoundedVector482 + crRightHandSideBoundedVector479*crRightHandSideBoundedVector56;
const double crRightHandSideBoundedVector724 =             0.16666666666666666*crRightHandSideBoundedVector4;
const double crRightHandSideBoundedVector725 =             1.1250000000000002*crRightHandSideBoundedVector316;
const double crRightHandSideBoundedVector726 =             crRightHandSideBoundedVector55*crRightHandSideBoundedVector658;
const double crRightHandSideBoundedVector727 =             crRightHandSideBoundedVector20*crRightHandSideBoundedVector688 + crRightHandSideBoundedVector326*crRightHandSideBoundedVector724 - crRightHandSideBoundedVector41*crRightHandSideBoundedVector712 + crRightHandSideBoundedVector415*crRightHandSideBoundedVector634 - crRightHandSideBoundedVector51*crRightHandSideBoundedVector726 - crRightHandSideBoundedVector587*crRightHandSideBoundedVector725 - crRightHandSideBoundedVector689 + crRightHandSideBoundedVector691;
const double crRightHandSideBoundedVector728 =             crRightHandSideBoundedVector685*crRightHandSideBoundedVector69;
const double crRightHandSideBoundedVector729 =             crRightHandSideBoundedVector20*crRightHandSideBoundedVector711 + crRightHandSideBoundedVector328*crRightHandSideBoundedVector724 + crRightHandSideBoundedVector499*crRightHandSideBoundedVector635 - crRightHandSideBoundedVector590*crRightHandSideBoundedVector725 - crRightHandSideBoundedVector61*crRightHandSideBoundedVector726 - crRightHandSideBoundedVector713 + crRightHandSideBoundedVector715 - crRightHandSideBoundedVector728;
const double crRightHandSideBoundedVector730 =             4.5000000000000009*crRightHandSideBoundedVector346;
const double crRightHandSideBoundedVector731 =             crRightHandSideBoundedVector68*crRightHandSideBoundedVector720;
const double crRightHandSideBoundedVector732 =             crRightHandSideBoundedVector336*crRightHandSideBoundedVector55;
const double crRightHandSideBoundedVector733 =             0.1111111111111111*r[0];
const double crRightHandSideBoundedVector734 =             crRightHandSideBoundedVector632*crRightHandSideBoundedVector684;
const double crRightHandSideBoundedVector735 =             crRightHandSideBoundedVector626*crRightHandSideBoundedVector687;
const double crRightHandSideBoundedVector736 =             crRightHandSideBoundedVector157*crRightHandSideBoundedVector735 - crRightHandSideBoundedVector160*crRightHandSideBoundedVector734 + crRightHandSideBoundedVector301*crRightHandSideBoundedVector634 - crRightHandSideBoundedVector326*crRightHandSideBoundedVector685 - crRightHandSideBoundedVector328*crRightHandSideBoundedVector686 + crRightHandSideBoundedVector478*crRightHandSideBoundedVector635 + crRightHandSideBoundedVector638 + crRightHandSideBoundedVector646 + crRightHandSideBoundedVector69*crRightHandSideBoundedVector735 - crRightHandSideBoundedVector72*crRightHandSideBoundedVector734 + crRightHandSideBoundedVector733;
const double crRightHandSideBoundedVector737 =             crRightHandSideBoundedVector642*crRightHandSideBoundedVector696;
const double crRightHandSideBoundedVector738 =             crRightHandSideBoundedVector110*crRightHandSideBoundedVector361;
const double crRightHandSideBoundedVector739 =             crRightHandSideBoundedVector139*crRightHandSideBoundedVector383;
const double crRightHandSideBoundedVector740 =             crRightHandSideBoundedVector301*crRightHandSideBoundedVector523 + 0.16666666666666666*crRightHandSideBoundedVector323 + 0.16666666666666666*crRightHandSideBoundedVector324 + 0.16666666666666666*crRightHandSideBoundedVector325 - 0.16666666666666666*crRightHandSideBoundedVector327 - 0.16666666666666666*crRightHandSideBoundedVector329 - 0.16666666666666666*crRightHandSideBoundedVector337 - 0.16666666666666666*crRightHandSideBoundedVector338 + crRightHandSideBoundedVector478*crRightHandSideBoundedVector524 + crRightHandSideBoundedVector55*crRightHandSideBoundedVector713 + crRightHandSideBoundedVector55*crRightHandSideBoundedVector728;
const double crRightHandSideBoundedVector741 =             DN_DX_2_1*crRightHandSideBoundedVector46;
const double crRightHandSideBoundedVector742 =             DN_DX_2_0*crRightHandSideBoundedVector56;
const double crRightHandSideBoundedVector743 =             DN_DX_2_1*crRightHandSideBoundedVector107;
const double crRightHandSideBoundedVector744 =             DN_DX_2_0*crRightHandSideBoundedVector111;
const double crRightHandSideBoundedVector745 =             DN_DX_2_1*crRightHandSideBoundedVector137;
const double crRightHandSideBoundedVector746 =             0.66666666666666663*crRightHandSideBoundedVector119;
const double crRightHandSideBoundedVector747 =             DN_DX_2_0*crRightHandSideBoundedVector140;
const double crRightHandSideBoundedVector748 =             1.5000000000000002*crRightHandSideBoundedVector121;
const double crRightHandSideBoundedVector749 =             crRightHandSideBoundedVector121*crRightHandSideBoundedVector261;
const double crRightHandSideBoundedVector750 =             crRightHandSideBoundedVector121*crRightHandSideBoundedVector560;
const double crRightHandSideBoundedVector751 =             DN_DX_2_0*crRightHandSideBoundedVector46;
const double crRightHandSideBoundedVector752 =             DN_DX_2_1*crRightHandSideBoundedVector56;
const double crRightHandSideBoundedVector753 =             crRightHandSideBoundedVector671 + crRightHandSideBoundedVector752;
const double crRightHandSideBoundedVector754 =             DN_DX_2_0*crRightHandSideBoundedVector107;
const double crRightHandSideBoundedVector755 =             DN_DX_2_1*crRightHandSideBoundedVector111;
const double crRightHandSideBoundedVector756 =             crRightHandSideBoundedVector299 + crRightHandSideBoundedVector755;
const double crRightHandSideBoundedVector757 =             crRightHandSideBoundedVector119*crRightHandSideBoundedVector282;
const double crRightHandSideBoundedVector758 =             DN_DX_2_0*crRightHandSideBoundedVector137;
const double crRightHandSideBoundedVector759 =             crRightHandSideBoundedVector121*crRightHandSideBoundedVector289;
const double crRightHandSideBoundedVector760 =             DN_DX_2_1*crRightHandSideBoundedVector140;
const double crRightHandSideBoundedVector761 =             crRightHandSideBoundedVector119*crRightHandSideBoundedVector291 - crRightHandSideBoundedVector750 + crRightHandSideBoundedVector760;
const double crRightHandSideBoundedVector762 =             DN_DX_2_0*crRightHandSideBoundedVector682;
const double crRightHandSideBoundedVector763 =             crRightHandSideBoundedVector136*crRightHandSideBoundedVector748;
const double crRightHandSideBoundedVector764 =             crRightHandSideBoundedVector139*crRightHandSideBoundedVector748;
const double crRightHandSideBoundedVector765 =             4.5*crRightHandSideBoundedVector438;
const double crRightHandSideBoundedVector766 =             crRightHandSideBoundedVector447*crRightHandSideBoundedVector765;
const double crRightHandSideBoundedVector767 =             crRightHandSideBoundedVector261*crRightHandSideBoundedVector386;
const double crRightHandSideBoundedVector768 =             crRightHandSideBoundedVector68*crRightHandSideBoundedVector767;
const double crRightHandSideBoundedVector769 =             crRightHandSideBoundedVector406 + crRightHandSideBoundedVector690 + 0.44444444444444442*f_ext(2,0);
const double crRightHandSideBoundedVector770 =             crRightHandSideBoundedVector121*crRightHandSideBoundedVector580;
const double crRightHandSideBoundedVector771 =             crRightHandSideBoundedVector704 + crRightHandSideBoundedVector751;
const double crRightHandSideBoundedVector772 =             crRightHandSideBoundedVector476 + crRightHandSideBoundedVector754;
const double crRightHandSideBoundedVector773 =             crRightHandSideBoundedVector119*crRightHandSideBoundedVector472 + crRightHandSideBoundedVector758 - crRightHandSideBoundedVector770;
const double crRightHandSideBoundedVector774 =             DN_DX_2_1*crRightHandSideBoundedVector682;
const double crRightHandSideBoundedVector775 =             crRightHandSideBoundedVector514*crRightHandSideBoundedVector765;
const double crRightHandSideBoundedVector776 =             crRightHandSideBoundedVector261*crRightHandSideBoundedVector449;
const double crRightHandSideBoundedVector777 =             crRightHandSideBoundedVector51*crRightHandSideBoundedVector776;
const double crRightHandSideBoundedVector778 =             crRightHandSideBoundedVector496 + crRightHandSideBoundedVector714 + 0.44444444444444442*f_ext(2,1);
const double crRightHandSideBoundedVector779 =             4.5000000000000009*crRightHandSideBoundedVector370;
const double crRightHandSideBoundedVector780 =             crRightHandSideBoundedVector68*crRightHandSideBoundedVector776;
const double crRightHandSideBoundedVector781 =             crRightHandSideBoundedVector649*crRightHandSideBoundedVector765;
            rRightHandSideBoundedVector[0]=-1.0*crRightHandSideBoundedVector0*crRightHandSideBoundedVector117 - 1.0*crRightHandSideBoundedVector0*crRightHandSideBoundedVector145 - 1.0*crRightHandSideBoundedVector0*crRightHandSideBoundedVector79 - 1.0*crRightHandSideBoundedVector146*crRightHandSideBoundedVector166 - 1.0*crRightHandSideBoundedVector146*crRightHandSideBoundedVector179 - 1.0*crRightHandSideBoundedVector146*crRightHandSideBoundedVector190 - 1.0*crRightHandSideBoundedVector191;
            rRightHandSideBoundedVector[1]=-DN_DX_0_0*crRightHandSideBoundedVector211 - DN_DX_0_0*crRightHandSideBoundedVector223 - DN_DX_0_0*crRightHandSideBoundedVector235 - crRightHandSideBoundedVector236*crRightHandSideBoundedVector242 - crRightHandSideBoundedVector243*crRightHandSideBoundedVector248 - crRightHandSideBoundedVector249*crRightHandSideBoundedVector254 + crRightHandSideBoundedVector264*(crRightHandSideBoundedVector240*crRightHandSideBoundedVector260 - crRightHandSideBoundedVector241*crRightHandSideBoundedVector262 - crRightHandSideBoundedVector255 - crRightHandSideBoundedVector256*crRightHandSideBoundedVector4 + crRightHandSideBoundedVector257*crRightHandSideBoundedVector69 + crRightHandSideBoundedVector258*crRightHandSideBoundedVector259) + crRightHandSideBoundedVector273*(crRightHandSideBoundedVector259*crRightHandSideBoundedVector266 - crRightHandSideBoundedVector265 + crRightHandSideBoundedVector272) + crRightHandSideBoundedVector280*(crRightHandSideBoundedVector259*crRightHandSideBoundedVector275 - crRightHandSideBoundedVector274 + crRightHandSideBoundedVector279) + crRightHandSideBoundedVector281*crRightHandSideBoundedVector56 + crRightHandSideBoundedVector284*crRightHandSideBoundedVector46 + crRightHandSideBoundedVector286*crRightHandSideBoundedVector45 - crRightHandSideBoundedVector294*(crRightHandSideBoundedVector199*crRightHandSideBoundedVector290 - crRightHandSideBoundedVector282*crRightHandSideBoundedVector288 - crRightHandSideBoundedVector283*crRightHandSideBoundedVector287 + crRightHandSideBoundedVector293) - crRightHandSideBoundedVector305*(-crRightHandSideBoundedVector282*crRightHandSideBoundedVector295 + crRightHandSideBoundedVector300 + crRightHandSideBoundedVector304) + 0.66666666666666663*crRightHandSideBoundedVector31 - crRightHandSideBoundedVector312*(-crRightHandSideBoundedVector282*crRightHandSideBoundedVector306 + crRightHandSideBoundedVector310 + crRightHandSideBoundedVector311) - crRightHandSideBoundedVector313*crRightHandSideBoundedVector344 - crRightHandSideBoundedVector313*crRightHandSideBoundedVector369 - crRightHandSideBoundedVector313*crRightHandSideBoundedVector390 - 0.66666666666666663*crRightHandSideBoundedVector37 - crRightHandSideBoundedVector412*(DN_DX_0_0*crRightHandSideBoundedVector36 + crRightHandSideBoundedVector283*crRightHandSideBoundedVector401 - crRightHandSideBoundedVector392*crRightHandSideBoundedVector395 - crRightHandSideBoundedVector396*crRightHandSideBoundedVector41 - crRightHandSideBoundedVector397*crRightHandSideBoundedVector51 + crRightHandSideBoundedVector400 + crRightHandSideBoundedVector403 - crRightHandSideBoundedVector405 + crRightHandSideBoundedVector408) - crRightHandSideBoundedVector435*(DN_DX_0_0*crRightHandSideBoundedVector103 - crRightHandSideBoundedVector414 + crRightHandSideBoundedVector433) - crRightHandSideBoundedVector454*(DN_DX_0_0*crRightHandSideBoundedVector134 - crRightHandSideBoundedVector437 + crRightHandSideBoundedVector453) + crRightHandSideBoundedVector455 + crRightHandSideBoundedVector456 - 0.66666666666666663*crRightHandSideBoundedVector47 - 0.66666666666666663*crRightHandSideBoundedVector57 - 1.0*crRightHandSideBoundedVector76;
            rRightHandSideBoundedVector[2]=-DN_DX_0_0*crRightHandSideBoundedVector464 - DN_DX_0_0*crRightHandSideBoundedVector465 - DN_DX_0_0*crRightHandSideBoundedVector466 - DN_DX_0_1*crRightHandSideBoundedVector461 - DN_DX_0_1*crRightHandSideBoundedVector462 - DN_DX_0_1*crRightHandSideBoundedVector463 + 0.66666666666666663*crRightHandSideBoundedVector150 - 0.66666666666666663*crRightHandSideBoundedVector153 - 0.66666666666666663*crRightHandSideBoundedVector154 - 0.66666666666666663*crRightHandSideBoundedVector155 - 1.0*crRightHandSideBoundedVector163 - crRightHandSideBoundedVector264*(crRightHandSideBoundedVector204*crRightHandSideBoundedVector290 - crRightHandSideBoundedVector282*crRightHandSideBoundedVector292 - crRightHandSideBoundedVector287*crRightHandSideBoundedVector468 + crRightHandSideBoundedVector473) - crRightHandSideBoundedVector273*(-crRightHandSideBoundedVector282*crRightHandSideBoundedVector296 + crRightHandSideBoundedVector477 + crRightHandSideBoundedVector480) - crRightHandSideBoundedVector280*(-crRightHandSideBoundedVector282*crRightHandSideBoundedVector307 + crRightHandSideBoundedVector485 + crRightHandSideBoundedVector486) - crRightHandSideBoundedVector294*(-crRightHandSideBoundedVector157*crRightHandSideBoundedVector257 + crRightHandSideBoundedVector240*crRightHandSideBoundedVector262 - crRightHandSideBoundedVector241*crRightHandSideBoundedVector260 - crRightHandSideBoundedVector255*crRightHandSideBoundedVector259 + crRightHandSideBoundedVector258 + crRightHandSideBoundedVector4*crRightHandSideBoundedVector474) - crRightHandSideBoundedVector305*(-crRightHandSideBoundedVector259*crRightHandSideBoundedVector265 + crRightHandSideBoundedVector266 + crRightHandSideBoundedVector483) - crRightHandSideBoundedVector312*(-crRightHandSideBoundedVector259*crRightHandSideBoundedVector274 + crRightHandSideBoundedVector275 + crRightHandSideBoundedVector487) - crRightHandSideBoundedVector344*crRightHandSideBoundedVector488 - crRightHandSideBoundedVector369*crRightHandSideBoundedVector488 - crRightHandSideBoundedVector390*crRightHandSideBoundedVector488 - crRightHandSideBoundedVector412*(-DN_DX_0_0*crRightHandSideBoundedVector339 + DN_DX_0_1*crRightHandSideBoundedVector152 - crRightHandSideBoundedVector392*crRightHandSideBoundedVector490 - crRightHandSideBoundedVector396*crRightHandSideBoundedVector68 - crRightHandSideBoundedVector397*crRightHandSideBoundedVector61 + crRightHandSideBoundedVector468*crRightHandSideBoundedVector493 + crRightHandSideBoundedVector492 + crRightHandSideBoundedVector495 + crRightHandSideBoundedVector498) - crRightHandSideBoundedVector435*(-DN_DX_0_0*crRightHandSideBoundedVector365 + DN_DX_0_1*crRightHandSideBoundedVector171 + crRightHandSideBoundedVector510) - crRightHandSideBoundedVector454*(-DN_DX_0_0*crRightHandSideBoundedVector387 + DN_DX_0_1*crRightHandSideBoundedVector183 + crRightHandSideBoundedVector519) + crRightHandSideBoundedVector46*crRightHandSideBoundedVector467 + crRightHandSideBoundedVector469*crRightHandSideBoundedVector56 + crRightHandSideBoundedVector471*crRightHandSideBoundedVector55 + crRightHandSideBoundedVector520 + crRightHandSideBoundedVector521;
            rRightHandSideBoundedVector[3]=-DN_DX_0_0*crRightHandSideBoundedVector542 - DN_DX_0_0*crRightHandSideBoundedVector557 - DN_DX_0_0*crRightHandSideBoundedVector572 - DN_DX_0_1*crRightHandSideBoundedVector576 - DN_DX_0_1*crRightHandSideBoundedVector579 - DN_DX_0_1*crRightHandSideBoundedVector582 - crRightHandSideBoundedVector199*crRightHandSideBoundedVector336 - crRightHandSideBoundedVector204*crRightHandSideBoundedVector336 - crRightHandSideBoundedVector264*(-DN_DX_0_0*crRightHandSideBoundedVector591 - DN_DX_0_1*crRightHandSideBoundedVector524 + crRightHandSideBoundedVector20*crRightHandSideBoundedVector492 + crRightHandSideBoundedVector257*crRightHandSideBoundedVector328 - crRightHandSideBoundedVector402*crRightHandSideBoundedVector61 - crRightHandSideBoundedVector494*crRightHandSideBoundedVector68 - crRightHandSideBoundedVector495 + crRightHandSideBoundedVector498 + crRightHandSideBoundedVector585*crRightHandSideBoundedVector589*crRightHandSideBoundedVector71 - crRightHandSideBoundedVector586*crRightHandSideBoundedVector590) - crRightHandSideBoundedVector273*(-DN_DX_0_0*crRightHandSideBoundedVector601 - DN_DX_0_1*crRightHandSideBoundedVector600 + crRightHandSideBoundedVector605) - crRightHandSideBoundedVector280*(-DN_DX_0_0*crRightHandSideBoundedVector615 - DN_DX_0_1*crRightHandSideBoundedVector614 + crRightHandSideBoundedVector619) + crRightHandSideBoundedVector283*crRightHandSideBoundedVector523 - crRightHandSideBoundedVector294*(-DN_DX_0_0*crRightHandSideBoundedVector523 + crRightHandSideBoundedVector20*crRightHandSideBoundedVector400 - crRightHandSideBoundedVector20*crRightHandSideBoundedVector405 + crRightHandSideBoundedVector257*crRightHandSideBoundedVector326 + crRightHandSideBoundedVector32*crRightHandSideBoundedVector584*crRightHandSideBoundedVector585 - crRightHandSideBoundedVector403 + crRightHandSideBoundedVector408 - crRightHandSideBoundedVector41*crRightHandSideBoundedVector494 - crRightHandSideBoundedVector586*crRightHandSideBoundedVector587 - crRightHandSideBoundedVector588) - crRightHandSideBoundedVector305*(-DN_DX_0_0*crRightHandSideBoundedVector593 - crRightHandSideBoundedVector20*crRightHandSideBoundedVector414 + crRightHandSideBoundedVector599) - crRightHandSideBoundedVector312*(-DN_DX_0_0*crRightHandSideBoundedVector607 - crRightHandSideBoundedVector20*crRightHandSideBoundedVector437 + crRightHandSideBoundedVector613) + 0.66666666666666663*crRightHandSideBoundedVector323 + 0.66666666666666663*crRightHandSideBoundedVector324 + 0.66666666666666663*crRightHandSideBoundedVector325 - 0.66666666666666663*crRightHandSideBoundedVector327 - 0.66666666666666663*crRightHandSideBoundedVector329 + crRightHandSideBoundedVector403*crRightHandSideBoundedVector45 - crRightHandSideBoundedVector412*(DN_DX_0_0*crRightHandSideBoundedVector628 + crRightHandSideBoundedVector157*crRightHandSideBoundedVector627 - crRightHandSideBoundedVector199*crRightHandSideBoundedVector633 - crRightHandSideBoundedVector204*crRightHandSideBoundedVector633 + crRightHandSideBoundedVector283*crRightHandSideBoundedVector634 - crRightHandSideBoundedVector326*crRightHandSideBoundedVector396 - crRightHandSideBoundedVector328*crRightHandSideBoundedVector397 + crRightHandSideBoundedVector335*crRightHandSideBoundedVector404 + crRightHandSideBoundedVector468*crRightHandSideBoundedVector635 + crRightHandSideBoundedVector624 + crRightHandSideBoundedVector625 + crRightHandSideBoundedVector627*crRightHandSideBoundedVector69 + 0.44444444444444442*r[0]) - crRightHandSideBoundedVector435*(DN_DX_0_0*crRightHandSideBoundedVector636 + crRightHandSideBoundedVector360*crRightHandSideBoundedVector413 + crRightHandSideBoundedVector644) + crRightHandSideBoundedVector45*crRightHandSideBoundedVector588 - crRightHandSideBoundedVector454*(DN_DX_0_0*crRightHandSideBoundedVector645 + crRightHandSideBoundedVector382*crRightHandSideBoundedVector436 + crRightHandSideBoundedVector651) + crRightHandSideBoundedVector468*crRightHandSideBoundedVector524 - crRightHandSideBoundedVector621*(crRightHandSideBoundedVector293 + crRightHandSideBoundedVector473) - crRightHandSideBoundedVector622*(crRightHandSideBoundedVector300 + crRightHandSideBoundedVector477) - crRightHandSideBoundedVector623*(crRightHandSideBoundedVector310 + crRightHandSideBoundedVector485) + crRightHandSideBoundedVector652 + crRightHandSideBoundedVector653;
            rRightHandSideBoundedVector[4]=-DN_DX_1_0*crRightHandSideBoundedVector294 - DN_DX_1_0*crRightHandSideBoundedVector305 - DN_DX_1_0*crRightHandSideBoundedVector312 - DN_DX_1_1*crRightHandSideBoundedVector264 - DN_DX_1_1*crRightHandSideBoundedVector273 - DN_DX_1_1*crRightHandSideBoundedVector280 - crRightHandSideBoundedVector654;
            rRightHandSideBoundedVector[5]=-DN_DX_1_0*crRightHandSideBoundedVector211 - DN_DX_1_0*crRightHandSideBoundedVector223 - DN_DX_1_0*crRightHandSideBoundedVector235 - DN_DX_1_1*crRightHandSideBoundedVector464 - DN_DX_1_1*crRightHandSideBoundedVector465 - DN_DX_1_1*crRightHandSideBoundedVector466 + 0.66666666666666663*crRightHandSideBoundedVector100 - 0.66666666666666663*crRightHandSideBoundedVector104 + crRightHandSideBoundedVector106*crRightHandSideBoundedVector667 + crRightHandSideBoundedVector107*crRightHandSideBoundedVector284 - 0.66666666666666663*crRightHandSideBoundedVector108 + crRightHandSideBoundedVector111*crRightHandSideBoundedVector281 - 0.66666666666666663*crRightHandSideBoundedVector112 + crRightHandSideBoundedVector264*(crRightHandSideBoundedVector259*crRightHandSideBoundedVector656 - crRightHandSideBoundedVector655 + crRightHandSideBoundedVector659) + crRightHandSideBoundedVector273*(crRightHandSideBoundedVector246*crRightHandSideBoundedVector663 - crRightHandSideBoundedVector247*crRightHandSideBoundedVector664 - crRightHandSideBoundedVector256*crRightHandSideBoundedVector82 + crRightHandSideBoundedVector259*crRightHandSideBoundedVector662 - crRightHandSideBoundedVector660 + crRightHandSideBoundedVector661*crRightHandSideBoundedVector69) + crRightHandSideBoundedVector280*(crRightHandSideBoundedVector259*crRightHandSideBoundedVector666 + crRightHandSideBoundedVector279 - crRightHandSideBoundedVector665) - crRightHandSideBoundedVector294*(-crRightHandSideBoundedVector282*crRightHandSideBoundedVector668 + crRightHandSideBoundedVector672 + crRightHandSideBoundedVector673) - crRightHandSideBoundedVector305*(crRightHandSideBoundedVector219*crRightHandSideBoundedVector676 - crRightHandSideBoundedVector282*crRightHandSideBoundedVector675 - crRightHandSideBoundedVector283*crRightHandSideBoundedVector674 + crRightHandSideBoundedVector678) - crRightHandSideBoundedVector312*(-crRightHandSideBoundedVector282*crRightHandSideBoundedVector679 + crRightHandSideBoundedVector311 + crRightHandSideBoundedVector681) - crRightHandSideBoundedVector344*crRightHandSideBoundedVector683 - crRightHandSideBoundedVector369*crRightHandSideBoundedVector683 - crRightHandSideBoundedVector390*crRightHandSideBoundedVector683 - crRightHandSideBoundedVector412*(DN_DX_1_0*crRightHandSideBoundedVector36 - DN_DX_1_1*crRightHandSideBoundedVector339 + crRightHandSideBoundedVector692) - crRightHandSideBoundedVector435*(DN_DX_1_0*crRightHandSideBoundedVector103 - DN_DX_1_1*crRightHandSideBoundedVector365 - crRightHandSideBoundedVector41*crRightHandSideBoundedVector694 - crRightHandSideBoundedVector419*crRightHandSideBoundedVector693 + crRightHandSideBoundedVector428*crRightHandSideBoundedVector698 - crRightHandSideBoundedVector51*crRightHandSideBoundedVector695 + crRightHandSideBoundedVector697 + crRightHandSideBoundedVector700 + crRightHandSideBoundedVector701) - crRightHandSideBoundedVector454*(DN_DX_1_0*crRightHandSideBoundedVector134 - DN_DX_1_1*crRightHandSideBoundedVector387 + crRightHandSideBoundedVector453) + crRightHandSideBoundedVector455 + crRightHandSideBoundedVector702;
            rRightHandSideBoundedVector[6]=-DN_DX_1_0*crRightHandSideBoundedVector464 - DN_DX_1_0*crRightHandSideBoundedVector465 - DN_DX_1_0*crRightHandSideBoundedVector466 - DN_DX_1_1*crRightHandSideBoundedVector461 - DN_DX_1_1*crRightHandSideBoundedVector462 - DN_DX_1_1*crRightHandSideBoundedVector463 + crRightHandSideBoundedVector107*crRightHandSideBoundedVector467 + crRightHandSideBoundedVector110*crRightHandSideBoundedVector703 + crRightHandSideBoundedVector111*crRightHandSideBoundedVector469 + 0.66666666666666663*crRightHandSideBoundedVector169 - 0.66666666666666663*crRightHandSideBoundedVector172 - 0.66666666666666663*crRightHandSideBoundedVector173 - 0.66666666666666663*crRightHandSideBoundedVector174 - crRightHandSideBoundedVector264*(-crRightHandSideBoundedVector282*crRightHandSideBoundedVector669 + crRightHandSideBoundedVector705 + crRightHandSideBoundedVector706) - crRightHandSideBoundedVector273*(crRightHandSideBoundedVector221*crRightHandSideBoundedVector676 - crRightHandSideBoundedVector282*crRightHandSideBoundedVector677 - crRightHandSideBoundedVector468*crRightHandSideBoundedVector674 + crRightHandSideBoundedVector708) - crRightHandSideBoundedVector280*(-crRightHandSideBoundedVector282*crRightHandSideBoundedVector680 + crRightHandSideBoundedVector486 + crRightHandSideBoundedVector709) - crRightHandSideBoundedVector294*(-crRightHandSideBoundedVector259*crRightHandSideBoundedVector655 + crRightHandSideBoundedVector656 + crRightHandSideBoundedVector707) - crRightHandSideBoundedVector305*(-crRightHandSideBoundedVector157*crRightHandSideBoundedVector661 + crRightHandSideBoundedVector246*crRightHandSideBoundedVector664 - crRightHandSideBoundedVector247*crRightHandSideBoundedVector663 - crRightHandSideBoundedVector259*crRightHandSideBoundedVector660 + crRightHandSideBoundedVector474*crRightHandSideBoundedVector82 + crRightHandSideBoundedVector662) - crRightHandSideBoundedVector312*(-crRightHandSideBoundedVector259*crRightHandSideBoundedVector665 + crRightHandSideBoundedVector487 + crRightHandSideBoundedVector666) - crRightHandSideBoundedVector344*crRightHandSideBoundedVector710 - crRightHandSideBoundedVector369*crRightHandSideBoundedVector710 - crRightHandSideBoundedVector390*crRightHandSideBoundedVector710 - crRightHandSideBoundedVector412*(-DN_DX_1_0*crRightHandSideBoundedVector339 + DN_DX_1_1*crRightHandSideBoundedVector152 + crRightHandSideBoundedVector716) - crRightHandSideBoundedVector435*(-DN_DX_1_0*crRightHandSideBoundedVector365 + DN_DX_1_1*crRightHandSideBoundedVector171 + crRightHandSideBoundedVector364*crRightHandSideBoundedVector719 - crRightHandSideBoundedVector501*crRightHandSideBoundedVector717 - crRightHandSideBoundedVector61*crRightHandSideBoundedVector695 - crRightHandSideBoundedVector68*crRightHandSideBoundedVector694 + crRightHandSideBoundedVector718 + crRightHandSideBoundedVector721 + crRightHandSideBoundedVector722) - crRightHandSideBoundedVector454*(-DN_DX_1_0*crRightHandSideBoundedVector387 + DN_DX_1_1*crRightHandSideBoundedVector183 + crRightHandSideBoundedVector519) + crRightHandSideBoundedVector520 + crRightHandSideBoundedVector723;
            rRightHandSideBoundedVector[7]=-DN_DX_1_0*crRightHandSideBoundedVector542 - DN_DX_1_0*crRightHandSideBoundedVector557 - DN_DX_1_0*crRightHandSideBoundedVector572 - DN_DX_1_1*crRightHandSideBoundedVector576 - DN_DX_1_1*crRightHandSideBoundedVector579 - DN_DX_1_1*crRightHandSideBoundedVector582 + crRightHandSideBoundedVector110*crRightHandSideBoundedVector721 + crRightHandSideBoundedVector110*crRightHandSideBoundedVector731 - crRightHandSideBoundedVector219*crRightHandSideBoundedVector361 - crRightHandSideBoundedVector221*crRightHandSideBoundedVector361 - crRightHandSideBoundedVector264*(-DN_DX_1_0*crRightHandSideBoundedVector591 - DN_DX_1_1*crRightHandSideBoundedVector524 + crRightHandSideBoundedVector729) - crRightHandSideBoundedVector273*(-DN_DX_1_0*crRightHandSideBoundedVector601 - DN_DX_1_1*crRightHandSideBoundedVector600 + crRightHandSideBoundedVector20*crRightHandSideBoundedVector718 + crRightHandSideBoundedVector328*crRightHandSideBoundedVector661 + crRightHandSideBoundedVector602*crRightHandSideBoundedVector717 - crRightHandSideBoundedVector603*crRightHandSideBoundedVector730 - crRightHandSideBoundedVector61*crRightHandSideBoundedVector699 - crRightHandSideBoundedVector721 + crRightHandSideBoundedVector722 - crRightHandSideBoundedVector731) - crRightHandSideBoundedVector280*(-DN_DX_1_0*crRightHandSideBoundedVector615 - DN_DX_1_1*crRightHandSideBoundedVector614 + crRightHandSideBoundedVector619) + crRightHandSideBoundedVector283*crRightHandSideBoundedVector593 - crRightHandSideBoundedVector294*(-DN_DX_1_0*crRightHandSideBoundedVector523 - DN_DX_1_1*crRightHandSideBoundedVector591 + crRightHandSideBoundedVector727) - crRightHandSideBoundedVector305*(-DN_DX_1_0*crRightHandSideBoundedVector593 - DN_DX_1_1*crRightHandSideBoundedVector601 + crRightHandSideBoundedVector20*crRightHandSideBoundedVector697 + crRightHandSideBoundedVector326*crRightHandSideBoundedVector661 - crRightHandSideBoundedVector41*crRightHandSideBoundedVector720 - crRightHandSideBoundedVector51*crRightHandSideBoundedVector699 + crRightHandSideBoundedVector596*crRightHandSideBoundedVector693 - crRightHandSideBoundedVector598*crRightHandSideBoundedVector730 - crRightHandSideBoundedVector700 + crRightHandSideBoundedVector701) - crRightHandSideBoundedVector312*(-DN_DX_1_0*crRightHandSideBoundedVector607 - DN_DX_1_1*crRightHandSideBoundedVector615 + crRightHandSideBoundedVector613) + 0.66666666666666663*crRightHandSideBoundedVector352 + 0.66666666666666663*crRightHandSideBoundedVector353 + 0.66666666666666663*crRightHandSideBoundedVector354 - 0.66666666666666663*crRightHandSideBoundedVector355 - 0.66666666666666663*crRightHandSideBoundedVector356 - crRightHandSideBoundedVector412*(DN_DX_1_0*crRightHandSideBoundedVector628 + DN_DX_1_1*crRightHandSideBoundedVector732 + crRightHandSideBoundedVector736) - crRightHandSideBoundedVector435*(DN_DX_1_0*crRightHandSideBoundedVector636 + DN_DX_1_1*crRightHandSideBoundedVector738 + crRightHandSideBoundedVector157*crRightHandSideBoundedVector737 - crRightHandSideBoundedVector219*crRightHandSideBoundedVector640 - crRightHandSideBoundedVector221*crRightHandSideBoundedVector640 + crRightHandSideBoundedVector283*crRightHandSideBoundedVector596 - crRightHandSideBoundedVector326*crRightHandSideBoundedVector694 - crRightHandSideBoundedVector328*crRightHandSideBoundedVector695 + crRightHandSideBoundedVector468*crRightHandSideBoundedVector602 + crRightHandSideBoundedVector625 + crRightHandSideBoundedVector69*crRightHandSideBoundedVector737 + crRightHandSideBoundedVector733 + 0.44444444444444442*r[1]) - crRightHandSideBoundedVector454*(DN_DX_1_0*crRightHandSideBoundedVector645 + DN_DX_1_1*crRightHandSideBoundedVector739 + crRightHandSideBoundedVector651) + crRightHandSideBoundedVector468*crRightHandSideBoundedVector600 - crRightHandSideBoundedVector621*(crRightHandSideBoundedVector672 + crRightHandSideBoundedVector705) - crRightHandSideBoundedVector622*(crRightHandSideBoundedVector678 + crRightHandSideBoundedVector708) - crRightHandSideBoundedVector623*(crRightHandSideBoundedVector681 + crRightHandSideBoundedVector709) + crRightHandSideBoundedVector652 + crRightHandSideBoundedVector740;
            rRightHandSideBoundedVector[8]=-DN_DX_2_0*crRightHandSideBoundedVector294 - DN_DX_2_0*crRightHandSideBoundedVector305 - DN_DX_2_0*crRightHandSideBoundedVector312 - DN_DX_2_1*crRightHandSideBoundedVector264 - DN_DX_2_1*crRightHandSideBoundedVector273 - DN_DX_2_1*crRightHandSideBoundedVector280 - crRightHandSideBoundedVector654;
            rRightHandSideBoundedVector[9]=-DN_DX_2_0*crRightHandSideBoundedVector211 - DN_DX_2_0*crRightHandSideBoundedVector223 - DN_DX_2_0*crRightHandSideBoundedVector235 - DN_DX_2_1*crRightHandSideBoundedVector464 - DN_DX_2_1*crRightHandSideBoundedVector465 - DN_DX_2_1*crRightHandSideBoundedVector466 + 0.66666666666666663*crRightHandSideBoundedVector131 - 0.66666666666666663*crRightHandSideBoundedVector135 + crRightHandSideBoundedVector136*crRightHandSideBoundedVector750 + crRightHandSideBoundedVector137*crRightHandSideBoundedVector284 - 0.66666666666666663*crRightHandSideBoundedVector138 + crRightHandSideBoundedVector140*crRightHandSideBoundedVector281 - 0.66666666666666663*crRightHandSideBoundedVector141 + crRightHandSideBoundedVector264*(crRightHandSideBoundedVector259*crRightHandSideBoundedVector742 + crRightHandSideBoundedVector659 - crRightHandSideBoundedVector741) + crRightHandSideBoundedVector273*(crRightHandSideBoundedVector259*crRightHandSideBoundedVector744 + crRightHandSideBoundedVector272 - crRightHandSideBoundedVector743) + crRightHandSideBoundedVector280*(-crRightHandSideBoundedVector119*crRightHandSideBoundedVector256 + crRightHandSideBoundedVector252*crRightHandSideBoundedVector748 - crRightHandSideBoundedVector253*crRightHandSideBoundedVector749 + crRightHandSideBoundedVector259*crRightHandSideBoundedVector747 + crRightHandSideBoundedVector69*crRightHandSideBoundedVector746 - crRightHandSideBoundedVector745) - crRightHandSideBoundedVector294*(-crRightHandSideBoundedVector282*crRightHandSideBoundedVector751 + crRightHandSideBoundedVector673 + crRightHandSideBoundedVector753) - crRightHandSideBoundedVector305*(-crRightHandSideBoundedVector282*crRightHandSideBoundedVector754 + crRightHandSideBoundedVector304 + crRightHandSideBoundedVector756) - crRightHandSideBoundedVector312*(crRightHandSideBoundedVector231*crRightHandSideBoundedVector759 - crRightHandSideBoundedVector282*crRightHandSideBoundedVector758 - crRightHandSideBoundedVector283*crRightHandSideBoundedVector757 + crRightHandSideBoundedVector761) - crRightHandSideBoundedVector344*crRightHandSideBoundedVector762 - crRightHandSideBoundedVector369*crRightHandSideBoundedVector762 - crRightHandSideBoundedVector390*crRightHandSideBoundedVector762 - crRightHandSideBoundedVector412*(DN_DX_2_0*crRightHandSideBoundedVector36 - DN_DX_2_1*crRightHandSideBoundedVector339 + crRightHandSideBoundedVector692) - crRightHandSideBoundedVector435*(DN_DX_2_0*crRightHandSideBoundedVector103 - DN_DX_2_1*crRightHandSideBoundedVector365 + crRightHandSideBoundedVector433) - crRightHandSideBoundedVector454*(DN_DX_2_0*crRightHandSideBoundedVector134 - DN_DX_2_1*crRightHandSideBoundedVector387 - crRightHandSideBoundedVector41*crRightHandSideBoundedVector763 - crRightHandSideBoundedVector441*crRightHandSideBoundedVector693 + crRightHandSideBoundedVector449*crRightHandSideBoundedVector698 - crRightHandSideBoundedVector51*crRightHandSideBoundedVector764 + crRightHandSideBoundedVector766 + crRightHandSideBoundedVector768 + crRightHandSideBoundedVector769) + crRightHandSideBoundedVector456 + crRightHandSideBoundedVector702;
            rRightHandSideBoundedVector[10]=-DN_DX_2_0*crRightHandSideBoundedVector464 - DN_DX_2_0*crRightHandSideBoundedVector465 - DN_DX_2_0*crRightHandSideBoundedVector466 - DN_DX_2_1*crRightHandSideBoundedVector461 - DN_DX_2_1*crRightHandSideBoundedVector462 - DN_DX_2_1*crRightHandSideBoundedVector463 + crRightHandSideBoundedVector137*crRightHandSideBoundedVector467 + crRightHandSideBoundedVector139*crRightHandSideBoundedVector770 + crRightHandSideBoundedVector140*crRightHandSideBoundedVector469 + 0.66666666666666663*crRightHandSideBoundedVector181 - 0.66666666666666663*crRightHandSideBoundedVector184 - 0.66666666666666663*crRightHandSideBoundedVector185 - 0.66666666666666663*crRightHandSideBoundedVector186 - crRightHandSideBoundedVector264*(-crRightHandSideBoundedVector282*crRightHandSideBoundedVector752 + crRightHandSideBoundedVector706 + crRightHandSideBoundedVector771) - crRightHandSideBoundedVector273*(-crRightHandSideBoundedVector282*crRightHandSideBoundedVector755 + crRightHandSideBoundedVector480 + crRightHandSideBoundedVector772) - crRightHandSideBoundedVector280*(crRightHandSideBoundedVector233*crRightHandSideBoundedVector759 - crRightHandSideBoundedVector282*crRightHandSideBoundedVector760 - crRightHandSideBoundedVector468*crRightHandSideBoundedVector757 + crRightHandSideBoundedVector773) - crRightHandSideBoundedVector294*(-crRightHandSideBoundedVector259*crRightHandSideBoundedVector741 + crRightHandSideBoundedVector707 + crRightHandSideBoundedVector742) - crRightHandSideBoundedVector305*(-crRightHandSideBoundedVector259*crRightHandSideBoundedVector743 + crRightHandSideBoundedVector483 + crRightHandSideBoundedVector744) - crRightHandSideBoundedVector312*(crRightHandSideBoundedVector119*crRightHandSideBoundedVector474 - crRightHandSideBoundedVector157*crRightHandSideBoundedVector746 + crRightHandSideBoundedVector252*crRightHandSideBoundedVector749 - crRightHandSideBoundedVector253*crRightHandSideBoundedVector748 - crRightHandSideBoundedVector259*crRightHandSideBoundedVector745 + crRightHandSideBoundedVector747) - crRightHandSideBoundedVector344*crRightHandSideBoundedVector774 - crRightHandSideBoundedVector369*crRightHandSideBoundedVector774 - crRightHandSideBoundedVector390*crRightHandSideBoundedVector774 - crRightHandSideBoundedVector412*(-DN_DX_2_0*crRightHandSideBoundedVector339 + DN_DX_2_1*crRightHandSideBoundedVector152 + crRightHandSideBoundedVector716) - crRightHandSideBoundedVector435*(-DN_DX_2_0*crRightHandSideBoundedVector365 + DN_DX_2_1*crRightHandSideBoundedVector171 + crRightHandSideBoundedVector510) - crRightHandSideBoundedVector454*(-DN_DX_2_0*crRightHandSideBoundedVector387 + DN_DX_2_1*crRightHandSideBoundedVector183 + crRightHandSideBoundedVector386*crRightHandSideBoundedVector719 - crRightHandSideBoundedVector512*crRightHandSideBoundedVector717 - crRightHandSideBoundedVector61*crRightHandSideBoundedVector764 - crRightHandSideBoundedVector68*crRightHandSideBoundedVector763 + crRightHandSideBoundedVector775 + crRightHandSideBoundedVector777 + crRightHandSideBoundedVector778) + crRightHandSideBoundedVector521 + crRightHandSideBoundedVector723;
            rRightHandSideBoundedVector[11]=-DN_DX_2_0*crRightHandSideBoundedVector542 - DN_DX_2_0*crRightHandSideBoundedVector557 - DN_DX_2_0*crRightHandSideBoundedVector572 - DN_DX_2_1*crRightHandSideBoundedVector576 - DN_DX_2_1*crRightHandSideBoundedVector579 - DN_DX_2_1*crRightHandSideBoundedVector582 + crRightHandSideBoundedVector139*crRightHandSideBoundedVector777 + crRightHandSideBoundedVector139*crRightHandSideBoundedVector780 - crRightHandSideBoundedVector231*crRightHandSideBoundedVector383 - crRightHandSideBoundedVector233*crRightHandSideBoundedVector383 - crRightHandSideBoundedVector264*(-DN_DX_2_0*crRightHandSideBoundedVector591 - DN_DX_2_1*crRightHandSideBoundedVector524 + crRightHandSideBoundedVector729) - crRightHandSideBoundedVector273*(-DN_DX_2_0*crRightHandSideBoundedVector601 - DN_DX_2_1*crRightHandSideBoundedVector600 + crRightHandSideBoundedVector605) - crRightHandSideBoundedVector280*(-DN_DX_2_0*crRightHandSideBoundedVector615 - DN_DX_2_1*crRightHandSideBoundedVector614 + crRightHandSideBoundedVector20*crRightHandSideBoundedVector775 + crRightHandSideBoundedVector328*crRightHandSideBoundedVector746 - crRightHandSideBoundedVector61*crRightHandSideBoundedVector767 + crRightHandSideBoundedVector616*crRightHandSideBoundedVector717 - crRightHandSideBoundedVector617*crRightHandSideBoundedVector779 - crRightHandSideBoundedVector777 + crRightHandSideBoundedVector778 - crRightHandSideBoundedVector780) + crRightHandSideBoundedVector283*crRightHandSideBoundedVector607 - crRightHandSideBoundedVector294*(-DN_DX_2_0*crRightHandSideBoundedVector523 - DN_DX_2_1*crRightHandSideBoundedVector591 + crRightHandSideBoundedVector727) - crRightHandSideBoundedVector305*(-DN_DX_2_0*crRightHandSideBoundedVector593 - DN_DX_2_1*crRightHandSideBoundedVector601 + crRightHandSideBoundedVector599) - crRightHandSideBoundedVector312*(-DN_DX_2_0*crRightHandSideBoundedVector607 - DN_DX_2_1*crRightHandSideBoundedVector615 + crRightHandSideBoundedVector20*crRightHandSideBoundedVector766 + crRightHandSideBoundedVector326*crRightHandSideBoundedVector746 - crRightHandSideBoundedVector41*crRightHandSideBoundedVector776 - crRightHandSideBoundedVector51*crRightHandSideBoundedVector767 + crRightHandSideBoundedVector610*crRightHandSideBoundedVector693 - crRightHandSideBoundedVector612*crRightHandSideBoundedVector779 - crRightHandSideBoundedVector768 + crRightHandSideBoundedVector769) + 0.66666666666666663*crRightHandSideBoundedVector375 + 0.66666666666666663*crRightHandSideBoundedVector376 + 0.66666666666666663*crRightHandSideBoundedVector377 - 0.66666666666666663*crRightHandSideBoundedVector378 - 0.66666666666666663*crRightHandSideBoundedVector379 - crRightHandSideBoundedVector412*(DN_DX_2_0*crRightHandSideBoundedVector628 + DN_DX_2_1*crRightHandSideBoundedVector732 + crRightHandSideBoundedVector736) - crRightHandSideBoundedVector435*(DN_DX_2_0*crRightHandSideBoundedVector636 + DN_DX_2_1*crRightHandSideBoundedVector738 + crRightHandSideBoundedVector644) - crRightHandSideBoundedVector454*(DN_DX_2_0*crRightHandSideBoundedVector645 + DN_DX_2_1*crRightHandSideBoundedVector739 + crRightHandSideBoundedVector157*crRightHandSideBoundedVector781 - crRightHandSideBoundedVector231*crRightHandSideBoundedVector647 - crRightHandSideBoundedVector233*crRightHandSideBoundedVector647 + crRightHandSideBoundedVector283*crRightHandSideBoundedVector610 - crRightHandSideBoundedVector326*crRightHandSideBoundedVector763 - crRightHandSideBoundedVector328*crRightHandSideBoundedVector764 + crRightHandSideBoundedVector468*crRightHandSideBoundedVector616 + crRightHandSideBoundedVector624 + crRightHandSideBoundedVector69*crRightHandSideBoundedVector781 + crRightHandSideBoundedVector733 + 0.44444444444444442*r[2]) + crRightHandSideBoundedVector468*crRightHandSideBoundedVector614 - crRightHandSideBoundedVector621*(crRightHandSideBoundedVector753 + crRightHandSideBoundedVector771) - crRightHandSideBoundedVector622*(crRightHandSideBoundedVector756 + crRightHandSideBoundedVector772) - crRightHandSideBoundedVector623*(crRightHandSideBoundedVector761 + crRightHandSideBoundedVector773) + crRightHandSideBoundedVector653 + crRightHandSideBoundedVector740;

    } else {
        const double crRightHandSideBoundedVector0 =             DN_DX_0_0*h;
const double crRightHandSideBoundedVector1 =             0.16666666666666666*U_1_0;
const double crRightHandSideBoundedVector2 =             0.16666666666666666*U_2_0;
const double crRightHandSideBoundedVector3 =             0.66666666666666663*U_0_0 + crRightHandSideBoundedVector1 + crRightHandSideBoundedVector2;
const double crRightHandSideBoundedVector4 =             1.0/crRightHandSideBoundedVector3;
const double crRightHandSideBoundedVector5 =             stab_c1/h;
const double crRightHandSideBoundedVector6 =             crRightHandSideBoundedVector4*crRightHandSideBoundedVector5;
const double crRightHandSideBoundedVector7 =             1.3333333333333333*mu;
const double crRightHandSideBoundedVector8 =             0.25*U_1_0;
const double crRightHandSideBoundedVector9 =             0.25*U_2_0;
const double crRightHandSideBoundedVector10 =             U_0_0 + crRightHandSideBoundedVector8 + crRightHandSideBoundedVector9;
const double crRightHandSideBoundedVector11 =             pow(crRightHandSideBoundedVector10, -2);
const double crRightHandSideBoundedVector12 =             0.25*U_1_1;
const double crRightHandSideBoundedVector13 =             0.25*U_2_1;
const double crRightHandSideBoundedVector14 =             pow(U_0_1 + crRightHandSideBoundedVector12 + crRightHandSideBoundedVector13, 2);
const double crRightHandSideBoundedVector15 =             0.25*U_1_2;
const double crRightHandSideBoundedVector16 =             0.25*U_2_2;
const double crRightHandSideBoundedVector17 =             pow(U_0_2 + crRightHandSideBoundedVector15 + crRightHandSideBoundedVector16, 2);
const double crRightHandSideBoundedVector18 =             crRightHandSideBoundedVector14 + crRightHandSideBoundedVector17;
const double crRightHandSideBoundedVector19 =             sqrt(gamma);
const double crRightHandSideBoundedVector20 =             gamma - 1;
const double crRightHandSideBoundedVector21 =             0.5*crRightHandSideBoundedVector11;
const double crRightHandSideBoundedVector22 =             0.16666666666666666*U_1_3;
const double crRightHandSideBoundedVector23 =             0.16666666666666666*U_2_3;
const double crRightHandSideBoundedVector24 =             0.66666666666666663*U_0_3 + crRightHandSideBoundedVector22 + crRightHandSideBoundedVector23;
const double crRightHandSideBoundedVector25 =             crRightHandSideBoundedVector19*sqrt(-crRightHandSideBoundedVector20*(crRightHandSideBoundedVector14*crRightHandSideBoundedVector21 + crRightHandSideBoundedVector17*crRightHandSideBoundedVector21 - crRightHandSideBoundedVector24*crRightHandSideBoundedVector4)) + 1.0*sqrt(crRightHandSideBoundedVector11*crRightHandSideBoundedVector18);
const double crRightHandSideBoundedVector26 =             crRightHandSideBoundedVector25*stab_c2;
const double crRightHandSideBoundedVector27 =             1.0/(crRightHandSideBoundedVector26 + crRightHandSideBoundedVector6*crRightHandSideBoundedVector7);
const double crRightHandSideBoundedVector28 =             0.16666666666666666*dUdt_1_1;
const double crRightHandSideBoundedVector29 =             0.16666666666666666*f_ext(1,0);
const double crRightHandSideBoundedVector30 =             0.16666666666666666*f_ext(2,0);
const double crRightHandSideBoundedVector31 =             crRightHandSideBoundedVector29 + crRightHandSideBoundedVector30 + 0.66666666666666663*f_ext(0,0);
const double crRightHandSideBoundedVector32 =             crRightHandSideBoundedVector3*crRightHandSideBoundedVector31;
const double crRightHandSideBoundedVector33 =             DN_DX_0_0*U_0_0 + DN_DX_1_0*U_1_0 + DN_DX_2_0*U_2_0;
const double crRightHandSideBoundedVector34 =             1.0000000000000002*crRightHandSideBoundedVector14;
const double crRightHandSideBoundedVector35 =             0.50000000000000011*crRightHandSideBoundedVector20;
const double crRightHandSideBoundedVector36 =             crRightHandSideBoundedVector18*crRightHandSideBoundedVector35;
const double crRightHandSideBoundedVector37 =             crRightHandSideBoundedVector11*(-crRightHandSideBoundedVector34 + crRightHandSideBoundedVector36);
const double crRightHandSideBoundedVector38 =             crRightHandSideBoundedVector33*crRightHandSideBoundedVector37;
const double crRightHandSideBoundedVector39 =             DN_DX_0_1*U_0_2;
const double crRightHandSideBoundedVector40 =             DN_DX_1_1*U_1_2;
const double crRightHandSideBoundedVector41 =             DN_DX_2_1*U_2_2;
const double crRightHandSideBoundedVector42 =             crRightHandSideBoundedVector39 + crRightHandSideBoundedVector40 + crRightHandSideBoundedVector41;
const double crRightHandSideBoundedVector43 =             0.66666666666666663*U_0_1;
const double crRightHandSideBoundedVector44 =             0.16666666666666666*U_1_1;
const double crRightHandSideBoundedVector45 =             0.16666666666666666*U_2_1;
const double crRightHandSideBoundedVector46 =             crRightHandSideBoundedVector43 + crRightHandSideBoundedVector44 + crRightHandSideBoundedVector45;
const double crRightHandSideBoundedVector47 =             crRightHandSideBoundedVector4*crRightHandSideBoundedVector46;
const double crRightHandSideBoundedVector48 =             crRightHandSideBoundedVector42*crRightHandSideBoundedVector47;
const double crRightHandSideBoundedVector49 =             DN_DX_0_1*U_0_1;
const double crRightHandSideBoundedVector50 =             DN_DX_1_1*U_1_1;
const double crRightHandSideBoundedVector51 =             DN_DX_2_1*U_2_1;
const double crRightHandSideBoundedVector52 =             crRightHandSideBoundedVector49 + crRightHandSideBoundedVector50 + crRightHandSideBoundedVector51;
const double crRightHandSideBoundedVector53 =             0.66666666666666663*U_0_2;
const double crRightHandSideBoundedVector54 =             0.16666666666666666*U_1_2;
const double crRightHandSideBoundedVector55 =             0.16666666666666666*U_2_2;
const double crRightHandSideBoundedVector56 =             crRightHandSideBoundedVector53 + crRightHandSideBoundedVector54 + crRightHandSideBoundedVector55;
const double crRightHandSideBoundedVector57 =             crRightHandSideBoundedVector4*crRightHandSideBoundedVector56;
const double crRightHandSideBoundedVector58 =             crRightHandSideBoundedVector52*crRightHandSideBoundedVector57;
const double crRightHandSideBoundedVector59 =             DN_DX_0_0*U_0_1;
const double crRightHandSideBoundedVector60 =             DN_DX_1_0*U_1_1;
const double crRightHandSideBoundedVector61 =             DN_DX_2_0*U_2_1;
const double crRightHandSideBoundedVector62 =             crRightHandSideBoundedVector59 + crRightHandSideBoundedVector60 + crRightHandSideBoundedVector61;
const double crRightHandSideBoundedVector63 =             1.0*gamma;
const double crRightHandSideBoundedVector64 =             3.0 - crRightHandSideBoundedVector63;
const double crRightHandSideBoundedVector65 =             crRightHandSideBoundedVector62*crRightHandSideBoundedVector64;
const double crRightHandSideBoundedVector66 =             DN_DX_0_0*U_0_2;
const double crRightHandSideBoundedVector67 =             DN_DX_1_0*U_1_2;
const double crRightHandSideBoundedVector68 =             DN_DX_2_0*U_2_2;
const double crRightHandSideBoundedVector69 =             crRightHandSideBoundedVector66 + crRightHandSideBoundedVector67 + crRightHandSideBoundedVector68;
const double crRightHandSideBoundedVector70 =             crRightHandSideBoundedVector20*crRightHandSideBoundedVector69;
const double crRightHandSideBoundedVector71 =             1.0*crRightHandSideBoundedVector70;
const double crRightHandSideBoundedVector72 =             DN_DX_0_1*U_0_0 + DN_DX_1_1*U_1_0 + DN_DX_2_1*U_2_0;
const double crRightHandSideBoundedVector73 =             crRightHandSideBoundedVector56*crRightHandSideBoundedVector72;
const double crRightHandSideBoundedVector74 =             2.2500000000000004*crRightHandSideBoundedVector46;
const double crRightHandSideBoundedVector75 =             crRightHandSideBoundedVector11*crRightHandSideBoundedVector74;
const double crRightHandSideBoundedVector76 =             DN_DX_0_0*U_0_3 + DN_DX_1_0*U_1_3 + DN_DX_2_0*U_2_3;
const double crRightHandSideBoundedVector77 =             crRightHandSideBoundedVector20*crRightHandSideBoundedVector76;
const double crRightHandSideBoundedVector78 =             crRightHandSideBoundedVector77 + 0.16666666666666666*dUdt_2_1;
const double crRightHandSideBoundedVector79 =             crRightHandSideBoundedVector27*(crRightHandSideBoundedVector28 - crRightHandSideBoundedVector32 + crRightHandSideBoundedVector38 + crRightHandSideBoundedVector47*crRightHandSideBoundedVector65 + crRightHandSideBoundedVector48 - crRightHandSideBoundedVector57*crRightHandSideBoundedVector71 + crRightHandSideBoundedVector58 - crRightHandSideBoundedVector73*crRightHandSideBoundedVector75 + crRightHandSideBoundedVector78 + 0.66666666666666663*dUdt_0_1);
const double crRightHandSideBoundedVector80 =             0.16666666666666666*U_0_0;
const double crRightHandSideBoundedVector81 =             0.66666666666666663*U_1_0 + crRightHandSideBoundedVector2 + crRightHandSideBoundedVector80;
const double crRightHandSideBoundedVector82 =             1.0/crRightHandSideBoundedVector81;
const double crRightHandSideBoundedVector83 =             crRightHandSideBoundedVector5*crRightHandSideBoundedVector7;
const double crRightHandSideBoundedVector84 =             0.25*U_0_0;
const double crRightHandSideBoundedVector85 =             U_1_0 + crRightHandSideBoundedVector84 + crRightHandSideBoundedVector9;
const double crRightHandSideBoundedVector86 =             pow(crRightHandSideBoundedVector85, -2);
const double crRightHandSideBoundedVector87 =             0.25*U_0_1;
const double crRightHandSideBoundedVector88 =             pow(U_1_1 + crRightHandSideBoundedVector13 + crRightHandSideBoundedVector87, 2);
const double crRightHandSideBoundedVector89 =             0.25*U_0_2;
const double crRightHandSideBoundedVector90 =             pow(U_1_2 + crRightHandSideBoundedVector16 + crRightHandSideBoundedVector89, 2);
const double crRightHandSideBoundedVector91 =             crRightHandSideBoundedVector88 + crRightHandSideBoundedVector90;
const double crRightHandSideBoundedVector92 =             0.5*crRightHandSideBoundedVector86;
const double crRightHandSideBoundedVector93 =             0.16666666666666666*U_0_3;
const double crRightHandSideBoundedVector94 =             0.66666666666666663*U_1_3 + crRightHandSideBoundedVector23 + crRightHandSideBoundedVector93;
const double crRightHandSideBoundedVector95 =             crRightHandSideBoundedVector19*sqrt(-crRightHandSideBoundedVector20*(-crRightHandSideBoundedVector82*crRightHandSideBoundedVector94 + crRightHandSideBoundedVector88*crRightHandSideBoundedVector92 + crRightHandSideBoundedVector90*crRightHandSideBoundedVector92)) + 1.0*sqrt(crRightHandSideBoundedVector86*crRightHandSideBoundedVector91);
const double crRightHandSideBoundedVector96 =             crRightHandSideBoundedVector95*stab_c2;
const double crRightHandSideBoundedVector97 =             1.0/(crRightHandSideBoundedVector82*crRightHandSideBoundedVector83 + crRightHandSideBoundedVector96);
const double crRightHandSideBoundedVector98 =             0.16666666666666666*dUdt_0_1;
const double crRightHandSideBoundedVector99 =             0.16666666666666666*f_ext(0,0);
const double crRightHandSideBoundedVector100 =             crRightHandSideBoundedVector30 + crRightHandSideBoundedVector99 + 0.66666666666666663*f_ext(1,0);
const double crRightHandSideBoundedVector101 =             crRightHandSideBoundedVector100*crRightHandSideBoundedVector81;
const double crRightHandSideBoundedVector102 =             1.0000000000000002*crRightHandSideBoundedVector88;
const double crRightHandSideBoundedVector103 =             crRightHandSideBoundedVector35*crRightHandSideBoundedVector91;
const double crRightHandSideBoundedVector104 =             crRightHandSideBoundedVector86*(-crRightHandSideBoundedVector102 + crRightHandSideBoundedVector103);
const double crRightHandSideBoundedVector105 =             crRightHandSideBoundedVector104*crRightHandSideBoundedVector33;
const double crRightHandSideBoundedVector106 =             0.16666666666666666*U_0_1;
const double crRightHandSideBoundedVector107 =             0.66666666666666663*U_1_1 + crRightHandSideBoundedVector106 + crRightHandSideBoundedVector45;
const double crRightHandSideBoundedVector108 =             crRightHandSideBoundedVector107*crRightHandSideBoundedVector82;
const double crRightHandSideBoundedVector109 =             crRightHandSideBoundedVector108*crRightHandSideBoundedVector42;
const double crRightHandSideBoundedVector110 =             0.16666666666666666*U_0_2;
const double crRightHandSideBoundedVector111 =             0.66666666666666663*U_1_2 + crRightHandSideBoundedVector110 + crRightHandSideBoundedVector55;
const double crRightHandSideBoundedVector112 =             crRightHandSideBoundedVector111*crRightHandSideBoundedVector82;
const double crRightHandSideBoundedVector113 =             crRightHandSideBoundedVector112*crRightHandSideBoundedVector52;
const double crRightHandSideBoundedVector114 =             2.2500000000000004*crRightHandSideBoundedVector107;
const double crRightHandSideBoundedVector115 =             crRightHandSideBoundedVector111*crRightHandSideBoundedVector72;
const double crRightHandSideBoundedVector116 =             crRightHandSideBoundedVector115*crRightHandSideBoundedVector86;
const double crRightHandSideBoundedVector117 =             crRightHandSideBoundedVector97*(-crRightHandSideBoundedVector101 + crRightHandSideBoundedVector105 + crRightHandSideBoundedVector108*crRightHandSideBoundedVector65 + crRightHandSideBoundedVector109 - crRightHandSideBoundedVector112*crRightHandSideBoundedVector71 + crRightHandSideBoundedVector113 - crRightHandSideBoundedVector114*crRightHandSideBoundedVector116 + crRightHandSideBoundedVector78 + crRightHandSideBoundedVector98 + 0.66666666666666663*dUdt_1_1);
const double crRightHandSideBoundedVector118 =             0.66666666666666663*U_2_0 + crRightHandSideBoundedVector1 + crRightHandSideBoundedVector80;
const double crRightHandSideBoundedVector119 =             1.0/crRightHandSideBoundedVector118;
const double crRightHandSideBoundedVector120 =             U_2_0 + crRightHandSideBoundedVector8 + crRightHandSideBoundedVector84;
const double crRightHandSideBoundedVector121 =             pow(crRightHandSideBoundedVector120, -2);
const double crRightHandSideBoundedVector122 =             pow(U_2_1 + crRightHandSideBoundedVector12 + crRightHandSideBoundedVector87, 2);
const double crRightHandSideBoundedVector123 =             pow(U_2_2 + crRightHandSideBoundedVector15 + crRightHandSideBoundedVector89, 2);
const double crRightHandSideBoundedVector124 =             crRightHandSideBoundedVector122 + crRightHandSideBoundedVector123;
const double crRightHandSideBoundedVector125 =             0.5*crRightHandSideBoundedVector121;
const double crRightHandSideBoundedVector126 =             0.66666666666666663*U_2_3 + crRightHandSideBoundedVector22 + crRightHandSideBoundedVector93;
const double crRightHandSideBoundedVector127 =             crRightHandSideBoundedVector19*sqrt(-crRightHandSideBoundedVector20*(-crRightHandSideBoundedVector119*crRightHandSideBoundedVector126 + crRightHandSideBoundedVector122*crRightHandSideBoundedVector125 + crRightHandSideBoundedVector123*crRightHandSideBoundedVector125)) + 1.0*sqrt(crRightHandSideBoundedVector121*crRightHandSideBoundedVector124);
const double crRightHandSideBoundedVector128 =             crRightHandSideBoundedVector127*stab_c2;
const double crRightHandSideBoundedVector129 =             1.0/(crRightHandSideBoundedVector119*crRightHandSideBoundedVector83 + crRightHandSideBoundedVector128);
const double crRightHandSideBoundedVector130 =             crRightHandSideBoundedVector29 + crRightHandSideBoundedVector99 + 0.66666666666666663*f_ext(2,0);
const double crRightHandSideBoundedVector131 =             crRightHandSideBoundedVector118*crRightHandSideBoundedVector130;
const double crRightHandSideBoundedVector132 =             0.66666666666666663*U_2_1 + crRightHandSideBoundedVector106 + crRightHandSideBoundedVector44;
const double crRightHandSideBoundedVector133 =             crRightHandSideBoundedVector119*crRightHandSideBoundedVector132;
const double crRightHandSideBoundedVector134 =             crRightHandSideBoundedVector133*crRightHandSideBoundedVector42;
const double crRightHandSideBoundedVector135 =             0.66666666666666663*U_2_2 + crRightHandSideBoundedVector110 + crRightHandSideBoundedVector54;
const double crRightHandSideBoundedVector136 =             crRightHandSideBoundedVector119*crRightHandSideBoundedVector135;
const double crRightHandSideBoundedVector137 =             crRightHandSideBoundedVector136*crRightHandSideBoundedVector52;
const double crRightHandSideBoundedVector138 =             2.2500000000000004*crRightHandSideBoundedVector132;
const double crRightHandSideBoundedVector139 =             crRightHandSideBoundedVector135*crRightHandSideBoundedVector72;
const double crRightHandSideBoundedVector140 =             crRightHandSideBoundedVector121*crRightHandSideBoundedVector139;
const double crRightHandSideBoundedVector141 =             1.0000000000000002*crRightHandSideBoundedVector122;
const double crRightHandSideBoundedVector142 =             crRightHandSideBoundedVector124*crRightHandSideBoundedVector35;
const double crRightHandSideBoundedVector143 =             crRightHandSideBoundedVector121*(-crRightHandSideBoundedVector141 + crRightHandSideBoundedVector142);
const double crRightHandSideBoundedVector144 =             crRightHandSideBoundedVector143*crRightHandSideBoundedVector33;
const double crRightHandSideBoundedVector145 =             crRightHandSideBoundedVector129*(-crRightHandSideBoundedVector131 + crRightHandSideBoundedVector133*crRightHandSideBoundedVector65 + crRightHandSideBoundedVector134 - crRightHandSideBoundedVector136*crRightHandSideBoundedVector71 + crRightHandSideBoundedVector137 - crRightHandSideBoundedVector138*crRightHandSideBoundedVector140 + crRightHandSideBoundedVector144 + crRightHandSideBoundedVector28 + crRightHandSideBoundedVector77 + crRightHandSideBoundedVector98 + 0.66666666666666663*dUdt_2_1);
const double crRightHandSideBoundedVector146 =             DN_DX_0_1*h;
const double crRightHandSideBoundedVector147 =             0.16666666666666666*dUdt_1_2;
const double crRightHandSideBoundedVector148 =             0.16666666666666666*f_ext(1,1);
const double crRightHandSideBoundedVector149 =             0.16666666666666666*f_ext(2,1);
const double crRightHandSideBoundedVector150 =             crRightHandSideBoundedVector148 + crRightHandSideBoundedVector149 + 0.66666666666666663*f_ext(0,1);
const double crRightHandSideBoundedVector151 =             crRightHandSideBoundedVector150*crRightHandSideBoundedVector3;
const double crRightHandSideBoundedVector152 =             1.0000000000000002*crRightHandSideBoundedVector17;
const double crRightHandSideBoundedVector153 =             crRightHandSideBoundedVector11*(-crRightHandSideBoundedVector152 + crRightHandSideBoundedVector36);
const double crRightHandSideBoundedVector154 =             crRightHandSideBoundedVector153*crRightHandSideBoundedVector72;
const double crRightHandSideBoundedVector155 =             crRightHandSideBoundedVector47*crRightHandSideBoundedVector69;
const double crRightHandSideBoundedVector156 =             crRightHandSideBoundedVector57*crRightHandSideBoundedVector62;
const double crRightHandSideBoundedVector157 =             crRightHandSideBoundedVector42*crRightHandSideBoundedVector64;
const double crRightHandSideBoundedVector158 =             crRightHandSideBoundedVector20*crRightHandSideBoundedVector52;
const double crRightHandSideBoundedVector159 =             1.0*crRightHandSideBoundedVector158;
const double crRightHandSideBoundedVector160 =             2.2500000000000004*crRightHandSideBoundedVector56;
const double crRightHandSideBoundedVector161 =             crRightHandSideBoundedVector33*crRightHandSideBoundedVector46;
const double crRightHandSideBoundedVector162 =             crRightHandSideBoundedVector11*crRightHandSideBoundedVector161;
const double crRightHandSideBoundedVector163 =             DN_DX_0_1*U_0_3 + DN_DX_1_1*U_1_3 + DN_DX_2_1*U_2_3;
const double crRightHandSideBoundedVector164 =             crRightHandSideBoundedVector163*crRightHandSideBoundedVector20;
const double crRightHandSideBoundedVector165 =             crRightHandSideBoundedVector164 + 0.16666666666666666*dUdt_2_2;
const double crRightHandSideBoundedVector166 =             crRightHandSideBoundedVector27*(crRightHandSideBoundedVector147 - crRightHandSideBoundedVector151 + crRightHandSideBoundedVector154 + crRightHandSideBoundedVector155 + crRightHandSideBoundedVector156 + crRightHandSideBoundedVector157*crRightHandSideBoundedVector57 - crRightHandSideBoundedVector159*crRightHandSideBoundedVector47 - crRightHandSideBoundedVector160*crRightHandSideBoundedVector162 + crRightHandSideBoundedVector165 + 0.66666666666666663*dUdt_0_2);
const double crRightHandSideBoundedVector167 =             0.16666666666666666*dUdt_0_2;
const double crRightHandSideBoundedVector168 =             0.16666666666666666*f_ext(0,1);
const double crRightHandSideBoundedVector169 =             crRightHandSideBoundedVector149 + crRightHandSideBoundedVector168 + 0.66666666666666663*f_ext(1,1);
const double crRightHandSideBoundedVector170 =             crRightHandSideBoundedVector169*crRightHandSideBoundedVector81;
const double crRightHandSideBoundedVector171 =             1.0000000000000002*crRightHandSideBoundedVector90;
const double crRightHandSideBoundedVector172 =             crRightHandSideBoundedVector86*(crRightHandSideBoundedVector103 - crRightHandSideBoundedVector171);
const double crRightHandSideBoundedVector173 =             crRightHandSideBoundedVector172*crRightHandSideBoundedVector72;
const double crRightHandSideBoundedVector174 =             crRightHandSideBoundedVector108*crRightHandSideBoundedVector69;
const double crRightHandSideBoundedVector175 =             crRightHandSideBoundedVector112*crRightHandSideBoundedVector62;
const double crRightHandSideBoundedVector176 =             2.2500000000000004*crRightHandSideBoundedVector111;
const double crRightHandSideBoundedVector177 =             crRightHandSideBoundedVector107*crRightHandSideBoundedVector33;
const double crRightHandSideBoundedVector178 =             crRightHandSideBoundedVector177*crRightHandSideBoundedVector86;
const double crRightHandSideBoundedVector179 =             crRightHandSideBoundedVector97*(-crRightHandSideBoundedVector108*crRightHandSideBoundedVector159 + crRightHandSideBoundedVector112*crRightHandSideBoundedVector157 + crRightHandSideBoundedVector165 + crRightHandSideBoundedVector167 - crRightHandSideBoundedVector170 + crRightHandSideBoundedVector173 + crRightHandSideBoundedVector174 + crRightHandSideBoundedVector175 - crRightHandSideBoundedVector176*crRightHandSideBoundedVector178 + 0.66666666666666663*dUdt_1_2);
const double crRightHandSideBoundedVector180 =             crRightHandSideBoundedVector148 + crRightHandSideBoundedVector168 + 0.66666666666666663*f_ext(2,1);
const double crRightHandSideBoundedVector181 =             crRightHandSideBoundedVector118*crRightHandSideBoundedVector180;
const double crRightHandSideBoundedVector182 =             crRightHandSideBoundedVector133*crRightHandSideBoundedVector69;
const double crRightHandSideBoundedVector183 =             crRightHandSideBoundedVector136*crRightHandSideBoundedVector62;
const double crRightHandSideBoundedVector184 =             2.2500000000000004*crRightHandSideBoundedVector135;
const double crRightHandSideBoundedVector185 =             crRightHandSideBoundedVector132*crRightHandSideBoundedVector33;
const double crRightHandSideBoundedVector186 =             crRightHandSideBoundedVector121*crRightHandSideBoundedVector185;
const double crRightHandSideBoundedVector187 =             1.0000000000000002*crRightHandSideBoundedVector123;
const double crRightHandSideBoundedVector188 =             crRightHandSideBoundedVector121*(crRightHandSideBoundedVector142 - crRightHandSideBoundedVector187);
const double crRightHandSideBoundedVector189 =             crRightHandSideBoundedVector188*crRightHandSideBoundedVector72;
const double crRightHandSideBoundedVector190 =             crRightHandSideBoundedVector129*(-crRightHandSideBoundedVector133*crRightHandSideBoundedVector159 + crRightHandSideBoundedVector136*crRightHandSideBoundedVector157 + crRightHandSideBoundedVector147 + crRightHandSideBoundedVector164 + crRightHandSideBoundedVector167 - crRightHandSideBoundedVector181 + crRightHandSideBoundedVector182 + crRightHandSideBoundedVector183 - crRightHandSideBoundedVector184*crRightHandSideBoundedVector186 + crRightHandSideBoundedVector189 + 0.66666666666666663*dUdt_2_2);
const double crRightHandSideBoundedVector191 =             crRightHandSideBoundedVector42 + crRightHandSideBoundedVector62;
const double crRightHandSideBoundedVector192 =             crRightHandSideBoundedVector3*crRightHandSideBoundedVector62;
const double crRightHandSideBoundedVector193 =             -crRightHandSideBoundedVector161 + crRightHandSideBoundedVector192;
const double crRightHandSideBoundedVector194 =             crRightHandSideBoundedVector11*mu;
const double crRightHandSideBoundedVector195 =             4.5000000000000009*crRightHandSideBoundedVector194;
const double crRightHandSideBoundedVector196 =             crRightHandSideBoundedVector3*crRightHandSideBoundedVector42;
const double crRightHandSideBoundedVector197 =             crRightHandSideBoundedVector196 - crRightHandSideBoundedVector73;
const double crRightHandSideBoundedVector198 =             -1.5000000000000002*crRightHandSideBoundedVector194*(crRightHandSideBoundedVector193 + crRightHandSideBoundedVector197);
const double crRightHandSideBoundedVector199 =             0.66666666666666663*crRightHandSideBoundedVector161;
const double crRightHandSideBoundedVector200 =             crRightHandSideBoundedVector4*(-0.66666666666666663*crRightHandSideBoundedVector192 + 1.3333333333333335*crRightHandSideBoundedVector196 + crRightHandSideBoundedVector199 - 1.3333333333333335*crRightHandSideBoundedVector73);
const double crRightHandSideBoundedVector201 =             1.0/lin_m_norm;
const double crRightHandSideBoundedVector202 =             crRightHandSideBoundedVector201*lin_m[0]*lin_m[1];
const double crRightHandSideBoundedVector203 =             crRightHandSideBoundedVector202*nu_st;
const double crRightHandSideBoundedVector204 =             0.66666666666666663*crRightHandSideBoundedVector73;
const double crRightHandSideBoundedVector205 =             crRightHandSideBoundedVector4*(-1.3333333333333335*crRightHandSideBoundedVector161 + 1.3333333333333335*crRightHandSideBoundedVector192 - 0.66666666666666663*crRightHandSideBoundedVector196 + crRightHandSideBoundedVector204);
const double crRightHandSideBoundedVector206 =             crRightHandSideBoundedVector201*pow(lin_m[0], 2);
const double crRightHandSideBoundedVector207 =             crRightHandSideBoundedVector206*nu_st;
const double crRightHandSideBoundedVector208 =             crRightHandSideBoundedVector202*nu_sc;
const double crRightHandSideBoundedVector209 =             1 - crRightHandSideBoundedVector206;
const double crRightHandSideBoundedVector210 =             crRightHandSideBoundedVector209*nu_sc;
const double crRightHandSideBoundedVector211 =             crRightHandSideBoundedVector193*crRightHandSideBoundedVector195 + crRightHandSideBoundedVector198 + crRightHandSideBoundedVector200*crRightHandSideBoundedVector203 - crRightHandSideBoundedVector200*crRightHandSideBoundedVector208 + crRightHandSideBoundedVector205*crRightHandSideBoundedVector207 + crRightHandSideBoundedVector205*crRightHandSideBoundedVector210;
const double crRightHandSideBoundedVector212 =             crRightHandSideBoundedVector62*crRightHandSideBoundedVector81;
const double crRightHandSideBoundedVector213 =             -crRightHandSideBoundedVector177 + crRightHandSideBoundedVector212;
const double crRightHandSideBoundedVector214 =             crRightHandSideBoundedVector86*mu;
const double crRightHandSideBoundedVector215 =             4.5000000000000009*crRightHandSideBoundedVector214;
const double crRightHandSideBoundedVector216 =             crRightHandSideBoundedVector42*crRightHandSideBoundedVector81;
const double crRightHandSideBoundedVector217 =             -crRightHandSideBoundedVector115 + crRightHandSideBoundedVector216;
const double crRightHandSideBoundedVector218 =             -1.5000000000000002*crRightHandSideBoundedVector214*(crRightHandSideBoundedVector213 + crRightHandSideBoundedVector217);
const double crRightHandSideBoundedVector219 =             0.66666666666666663*crRightHandSideBoundedVector177;
const double crRightHandSideBoundedVector220 =             crRightHandSideBoundedVector82*(-1.3333333333333335*crRightHandSideBoundedVector115 - 0.66666666666666663*crRightHandSideBoundedVector212 + 1.3333333333333335*crRightHandSideBoundedVector216 + crRightHandSideBoundedVector219);
const double crRightHandSideBoundedVector221 =             0.66666666666666663*crRightHandSideBoundedVector115;
const double crRightHandSideBoundedVector222 =             crRightHandSideBoundedVector82*(-1.3333333333333335*crRightHandSideBoundedVector177 + 1.3333333333333335*crRightHandSideBoundedVector212 - 0.66666666666666663*crRightHandSideBoundedVector216 + crRightHandSideBoundedVector221);
const double crRightHandSideBoundedVector223 =             crRightHandSideBoundedVector203*crRightHandSideBoundedVector220 + crRightHandSideBoundedVector207*crRightHandSideBoundedVector222 - crRightHandSideBoundedVector208*crRightHandSideBoundedVector220 + crRightHandSideBoundedVector210*crRightHandSideBoundedVector222 + crRightHandSideBoundedVector213*crRightHandSideBoundedVector215 + crRightHandSideBoundedVector218;
const double crRightHandSideBoundedVector224 =             crRightHandSideBoundedVector118*crRightHandSideBoundedVector62;
const double crRightHandSideBoundedVector225 =             -crRightHandSideBoundedVector185 + crRightHandSideBoundedVector224;
const double crRightHandSideBoundedVector226 =             crRightHandSideBoundedVector121*mu;
const double crRightHandSideBoundedVector227 =             4.5000000000000009*crRightHandSideBoundedVector226;
const double crRightHandSideBoundedVector228 =             crRightHandSideBoundedVector118*crRightHandSideBoundedVector42;
const double crRightHandSideBoundedVector229 =             -crRightHandSideBoundedVector139 + crRightHandSideBoundedVector228;
const double crRightHandSideBoundedVector230 =             -1.5000000000000002*crRightHandSideBoundedVector226*(crRightHandSideBoundedVector225 + crRightHandSideBoundedVector229);
const double crRightHandSideBoundedVector231 =             0.66666666666666663*crRightHandSideBoundedVector185;
const double crRightHandSideBoundedVector232 =             crRightHandSideBoundedVector119*(-1.3333333333333335*crRightHandSideBoundedVector139 - 0.66666666666666663*crRightHandSideBoundedVector224 + 1.3333333333333335*crRightHandSideBoundedVector228 + crRightHandSideBoundedVector231);
const double crRightHandSideBoundedVector233 =             0.66666666666666663*crRightHandSideBoundedVector139;
const double crRightHandSideBoundedVector234 =             crRightHandSideBoundedVector119*(-1.3333333333333335*crRightHandSideBoundedVector185 + 1.3333333333333335*crRightHandSideBoundedVector224 - 0.66666666666666663*crRightHandSideBoundedVector228 + crRightHandSideBoundedVector233);
const double crRightHandSideBoundedVector235 =             crRightHandSideBoundedVector203*crRightHandSideBoundedVector232 + crRightHandSideBoundedVector207*crRightHandSideBoundedVector234 - crRightHandSideBoundedVector208*crRightHandSideBoundedVector232 + crRightHandSideBoundedVector210*crRightHandSideBoundedVector234 + crRightHandSideBoundedVector225*crRightHandSideBoundedVector227 + crRightHandSideBoundedVector230;
const double crRightHandSideBoundedVector236 =             DN_DX_0_1*crRightHandSideBoundedVector11;
const double crRightHandSideBoundedVector237 =             nu_sc*(1 - crRightHandSideBoundedVector202);
const double crRightHandSideBoundedVector238 =             crRightHandSideBoundedVector3*crRightHandSideBoundedVector69;
const double crRightHandSideBoundedVector239 =             crRightHandSideBoundedVector3*crRightHandSideBoundedVector52;
const double crRightHandSideBoundedVector240 =             crRightHandSideBoundedVector46*crRightHandSideBoundedVector72;
const double crRightHandSideBoundedVector241 =             crRightHandSideBoundedVector33*crRightHandSideBoundedVector56;
const double crRightHandSideBoundedVector242 =             (crRightHandSideBoundedVector203*crRightHandSideBoundedVector3 + crRightHandSideBoundedVector237*crRightHandSideBoundedVector3 + mu)*(2.2500000000000004*crRightHandSideBoundedVector238 + 2.2500000000000004*crRightHandSideBoundedVector239 - 2.2500000000000004*crRightHandSideBoundedVector240 - 2.2500000000000004*crRightHandSideBoundedVector241);
const double crRightHandSideBoundedVector243 =             DN_DX_0_1*crRightHandSideBoundedVector86;
const double crRightHandSideBoundedVector244 =             crRightHandSideBoundedVector69*crRightHandSideBoundedVector81;
const double crRightHandSideBoundedVector245 =             crRightHandSideBoundedVector52*crRightHandSideBoundedVector81;
const double crRightHandSideBoundedVector246 =             crRightHandSideBoundedVector107*crRightHandSideBoundedVector72;
const double crRightHandSideBoundedVector247 =             crRightHandSideBoundedVector111*crRightHandSideBoundedVector33;
const double crRightHandSideBoundedVector248 =             (crRightHandSideBoundedVector203*crRightHandSideBoundedVector81 + crRightHandSideBoundedVector237*crRightHandSideBoundedVector81 + mu)*(2.2500000000000004*crRightHandSideBoundedVector244 + 2.2500000000000004*crRightHandSideBoundedVector245 - 2.2500000000000004*crRightHandSideBoundedVector246 - 2.2500000000000004*crRightHandSideBoundedVector247);
const double crRightHandSideBoundedVector249 =             DN_DX_0_1*crRightHandSideBoundedVector121;
const double crRightHandSideBoundedVector250 =             crRightHandSideBoundedVector118*crRightHandSideBoundedVector69;
const double crRightHandSideBoundedVector251 =             crRightHandSideBoundedVector118*crRightHandSideBoundedVector52;
const double crRightHandSideBoundedVector252 =             crRightHandSideBoundedVector132*crRightHandSideBoundedVector72;
const double crRightHandSideBoundedVector253 =             crRightHandSideBoundedVector135*crRightHandSideBoundedVector33;
const double crRightHandSideBoundedVector254 =             (crRightHandSideBoundedVector118*crRightHandSideBoundedVector203 + crRightHandSideBoundedVector118*crRightHandSideBoundedVector237 + mu)*(2.2500000000000004*crRightHandSideBoundedVector250 + 2.2500000000000004*crRightHandSideBoundedVector251 - 2.2500000000000004*crRightHandSideBoundedVector252 - 2.2500000000000004*crRightHandSideBoundedVector253);
const double crRightHandSideBoundedVector255 =             DN_DX_0_1*crRightHandSideBoundedVector47;
const double crRightHandSideBoundedVector256 =             DN_DX_0_1*crRightHandSideBoundedVector43 + 0.66666666666666663*crRightHandSideBoundedVector50 + 0.66666666666666663*crRightHandSideBoundedVector51;
const double crRightHandSideBoundedVector257 =             0.66666666666666663*crRightHandSideBoundedVector4;
const double crRightHandSideBoundedVector258 =             DN_DX_0_0*crRightHandSideBoundedVector57;
const double crRightHandSideBoundedVector259 =             1.0*crRightHandSideBoundedVector20;
const double crRightHandSideBoundedVector260 =             1.5000000000000002*crRightHandSideBoundedVector11;
const double crRightHandSideBoundedVector261 =             1.5000000000000002*crRightHandSideBoundedVector20;
const double crRightHandSideBoundedVector262 =             crRightHandSideBoundedVector11*crRightHandSideBoundedVector261;
const double crRightHandSideBoundedVector263 =             1.0*h;
const double crRightHandSideBoundedVector264 =             crRightHandSideBoundedVector166*crRightHandSideBoundedVector263;
const double crRightHandSideBoundedVector265 =             DN_DX_0_1*crRightHandSideBoundedVector108;
const double crRightHandSideBoundedVector266 =             DN_DX_0_0*crRightHandSideBoundedVector112;
const double crRightHandSideBoundedVector267 =             DN_DX_1_1*crRightHandSideBoundedVector44 + DN_DX_2_1*crRightHandSideBoundedVector45 + 0.16666666666666666*crRightHandSideBoundedVector49;
const double crRightHandSideBoundedVector268 =             0.16666666666666666*crRightHandSideBoundedVector70;
const double crRightHandSideBoundedVector269 =             0.37500000000000006*crRightHandSideBoundedVector86;
const double crRightHandSideBoundedVector270 =             crRightHandSideBoundedVector246*crRightHandSideBoundedVector269;
const double crRightHandSideBoundedVector271 =             crRightHandSideBoundedVector247*crRightHandSideBoundedVector269;
const double crRightHandSideBoundedVector272 =             -crRightHandSideBoundedVector20*crRightHandSideBoundedVector271 - crRightHandSideBoundedVector267*crRightHandSideBoundedVector82 + crRightHandSideBoundedVector268*crRightHandSideBoundedVector82 + crRightHandSideBoundedVector270;
const double crRightHandSideBoundedVector273 =             crRightHandSideBoundedVector179*crRightHandSideBoundedVector263;
const double crRightHandSideBoundedVector274 =             DN_DX_0_1*crRightHandSideBoundedVector133;
const double crRightHandSideBoundedVector275 =             DN_DX_0_0*crRightHandSideBoundedVector136;
const double crRightHandSideBoundedVector276 =             0.37500000000000006*crRightHandSideBoundedVector121;
const double crRightHandSideBoundedVector277 =             crRightHandSideBoundedVector252*crRightHandSideBoundedVector276;
const double crRightHandSideBoundedVector278 =             crRightHandSideBoundedVector253*crRightHandSideBoundedVector276;
const double crRightHandSideBoundedVector279 =             -crRightHandSideBoundedVector119*crRightHandSideBoundedVector267 + crRightHandSideBoundedVector119*crRightHandSideBoundedVector268 - crRightHandSideBoundedVector20*crRightHandSideBoundedVector278 + crRightHandSideBoundedVector277;
const double crRightHandSideBoundedVector280 =             crRightHandSideBoundedVector190*crRightHandSideBoundedVector263;
const double crRightHandSideBoundedVector281 =             0.66666666666666663*crRightHandSideBoundedVector70;
const double crRightHandSideBoundedVector282 =             crRightHandSideBoundedVector63 - 3.0;
const double crRightHandSideBoundedVector283 =             0.66666666666666663*crRightHandSideBoundedVector62;
const double crRightHandSideBoundedVector284 =             crRightHandSideBoundedVector282*crRightHandSideBoundedVector283;
const double crRightHandSideBoundedVector285 =             1.5000000000000002*crRightHandSideBoundedVector73;
const double crRightHandSideBoundedVector286 =             crRightHandSideBoundedVector11*crRightHandSideBoundedVector285;
const double crRightHandSideBoundedVector287 =             crRightHandSideBoundedVector282*crRightHandSideBoundedVector4;
const double crRightHandSideBoundedVector288 =             DN_DX_0_0*crRightHandSideBoundedVector47;
const double crRightHandSideBoundedVector289 =             2.2500000000000004*gamma - 6.7500000000000018;
const double crRightHandSideBoundedVector290 =             crRightHandSideBoundedVector11*crRightHandSideBoundedVector289;
const double crRightHandSideBoundedVector291 =             0.66666666666666663*crRightHandSideBoundedVector39 + 0.66666666666666663*crRightHandSideBoundedVector40 + 0.66666666666666663*crRightHandSideBoundedVector41;
const double crRightHandSideBoundedVector292 =             DN_DX_0_1*crRightHandSideBoundedVector57;
const double crRightHandSideBoundedVector293 =             -crRightHandSideBoundedVector286 + crRightHandSideBoundedVector291*crRightHandSideBoundedVector4 + crRightHandSideBoundedVector292;
const double crRightHandSideBoundedVector294 =             crRightHandSideBoundedVector263*crRightHandSideBoundedVector79;
const double crRightHandSideBoundedVector295 =             DN_DX_0_0*crRightHandSideBoundedVector108;
const double crRightHandSideBoundedVector296 =             DN_DX_0_1*crRightHandSideBoundedVector112;
const double crRightHandSideBoundedVector297 =             0.16666666666666666*crRightHandSideBoundedVector39 + 0.16666666666666666*crRightHandSideBoundedVector40 + 0.16666666666666666*crRightHandSideBoundedVector41;
const double crRightHandSideBoundedVector298 =             crRightHandSideBoundedVector115*crRightHandSideBoundedVector269;
const double crRightHandSideBoundedVector299 =             crRightHandSideBoundedVector297*crRightHandSideBoundedVector82 - crRightHandSideBoundedVector298;
const double crRightHandSideBoundedVector300 =             crRightHandSideBoundedVector296 + crRightHandSideBoundedVector299;
const double crRightHandSideBoundedVector301 =             0.16666666666666666*crRightHandSideBoundedVector62;
const double crRightHandSideBoundedVector302 =             crRightHandSideBoundedVector282*crRightHandSideBoundedVector301;
const double crRightHandSideBoundedVector303 =             0.16666666666666666*crRightHandSideBoundedVector289;
const double crRightHandSideBoundedVector304 =             crRightHandSideBoundedVector178*crRightHandSideBoundedVector303 - crRightHandSideBoundedVector302*crRightHandSideBoundedVector82;
const double crRightHandSideBoundedVector305 =             crRightHandSideBoundedVector117*crRightHandSideBoundedVector263;
const double crRightHandSideBoundedVector306 =             DN_DX_0_0*crRightHandSideBoundedVector133;
const double crRightHandSideBoundedVector307 =             DN_DX_0_1*crRightHandSideBoundedVector136;
const double crRightHandSideBoundedVector308 =             crRightHandSideBoundedVector139*crRightHandSideBoundedVector276;
const double crRightHandSideBoundedVector309 =             crRightHandSideBoundedVector119*crRightHandSideBoundedVector297 - crRightHandSideBoundedVector308;
const double crRightHandSideBoundedVector310 =             crRightHandSideBoundedVector307 + crRightHandSideBoundedVector309;
const double crRightHandSideBoundedVector311 =             -crRightHandSideBoundedVector119*crRightHandSideBoundedVector302 + crRightHandSideBoundedVector186*crRightHandSideBoundedVector303;
const double crRightHandSideBoundedVector312 =             crRightHandSideBoundedVector145*crRightHandSideBoundedVector263;
const double crRightHandSideBoundedVector313 =             crRightHandSideBoundedVector0*crRightHandSideBoundedVector259;
const double crRightHandSideBoundedVector314 =             lambda/c_v;
const double crRightHandSideBoundedVector315 =             crRightHandSideBoundedVector314/gamma;
const double crRightHandSideBoundedVector316 =             0.16666666666666666*dUdt_1_3;
const double crRightHandSideBoundedVector317 =             0.16666666666666666*dUdt_2_3;
const double crRightHandSideBoundedVector318 =             0.16666666666666666*r[1];
const double crRightHandSideBoundedVector319 =             0.16666666666666666*r[2];
const double crRightHandSideBoundedVector320 =             crRightHandSideBoundedVector3*(crRightHandSideBoundedVector318 + crRightHandSideBoundedVector319 + 0.66666666666666663*r[0]);
const double crRightHandSideBoundedVector321 =             crRightHandSideBoundedVector31*crRightHandSideBoundedVector46;
const double crRightHandSideBoundedVector322 =             crRightHandSideBoundedVector150*crRightHandSideBoundedVector56;
const double crRightHandSideBoundedVector323 =             crRightHandSideBoundedVector76*gamma;
const double crRightHandSideBoundedVector324 =             crRightHandSideBoundedVector323*crRightHandSideBoundedVector47;
const double crRightHandSideBoundedVector325 =             crRightHandSideBoundedVector163*gamma;
const double crRightHandSideBoundedVector326 =             crRightHandSideBoundedVector325*crRightHandSideBoundedVector57;
const double crRightHandSideBoundedVector327 =             crRightHandSideBoundedVector56*crRightHandSideBoundedVector75;
const double crRightHandSideBoundedVector328 =             crRightHandSideBoundedVector52*crRightHandSideBoundedVector74;
const double crRightHandSideBoundedVector329 =             crRightHandSideBoundedVector11*crRightHandSideBoundedVector20;
const double crRightHandSideBoundedVector330 =             crRightHandSideBoundedVector329*crRightHandSideBoundedVector56;
const double crRightHandSideBoundedVector331 =             crRightHandSideBoundedVector329*crRightHandSideBoundedVector34;
const double crRightHandSideBoundedVector332 =             crRightHandSideBoundedVector20*(crRightHandSideBoundedVector24 - crRightHandSideBoundedVector4*(0.22222222222222221*crRightHandSideBoundedVector14 + 0.22222222222222221*crRightHandSideBoundedVector17));
const double crRightHandSideBoundedVector333 =             crRightHandSideBoundedVector4*(crRightHandSideBoundedVector24 + crRightHandSideBoundedVector332);
const double crRightHandSideBoundedVector334 =             crRightHandSideBoundedVector152*crRightHandSideBoundedVector329;
const double crRightHandSideBoundedVector335 =             -0.37500000000000006*U_1_3;
const double crRightHandSideBoundedVector336 =             -0.37500000000000006*U_2_3;
const double crRightHandSideBoundedVector337 =             0.75000000000000011*crRightHandSideBoundedVector20;
const double crRightHandSideBoundedVector338 =             1.0/crRightHandSideBoundedVector10;
const double crRightHandSideBoundedVector339 =             crRightHandSideBoundedVector18*crRightHandSideBoundedVector338;
const double crRightHandSideBoundedVector340 =             -1.5000000000000002*U_0_3 - 2.2500000000000004*crRightHandSideBoundedVector332 + crRightHandSideBoundedVector335 + crRightHandSideBoundedVector336 + crRightHandSideBoundedVector337*crRightHandSideBoundedVector339;
const double crRightHandSideBoundedVector341 =             crRightHandSideBoundedVector11*crRightHandSideBoundedVector340;
const double crRightHandSideBoundedVector342 =             crRightHandSideBoundedVector161*crRightHandSideBoundedVector341;
const double crRightHandSideBoundedVector343 =             crRightHandSideBoundedVector341*crRightHandSideBoundedVector73;
const double crRightHandSideBoundedVector344 =             (crRightHandSideBoundedVector316 + crRightHandSideBoundedVector317 - crRightHandSideBoundedVector320 - crRightHandSideBoundedVector321 - crRightHandSideBoundedVector322 + crRightHandSideBoundedVector324 + crRightHandSideBoundedVector326 - crRightHandSideBoundedVector327*crRightHandSideBoundedVector70 - crRightHandSideBoundedVector328*crRightHandSideBoundedVector330 + crRightHandSideBoundedVector342 + crRightHandSideBoundedVector343 + crRightHandSideBoundedVector42*(crRightHandSideBoundedVector333 - crRightHandSideBoundedVector334) + crRightHandSideBoundedVector62*(-crRightHandSideBoundedVector331 + crRightHandSideBoundedVector333) + 0.66666666666666663*dUdt_0_3)/(crRightHandSideBoundedVector26 + crRightHandSideBoundedVector315*crRightHandSideBoundedVector6);
const double crRightHandSideBoundedVector345 =             crRightHandSideBoundedVector315*crRightHandSideBoundedVector5;
const double crRightHandSideBoundedVector346 =             0.16666666666666666*dUdt_0_3;
const double crRightHandSideBoundedVector347 =             0.16666666666666666*r[0];
const double crRightHandSideBoundedVector348 =             crRightHandSideBoundedVector81*(crRightHandSideBoundedVector319 + crRightHandSideBoundedVector347 + 0.66666666666666663*r[1]);
const double crRightHandSideBoundedVector349 =             crRightHandSideBoundedVector100*crRightHandSideBoundedVector107;
const double crRightHandSideBoundedVector350 =             crRightHandSideBoundedVector111*crRightHandSideBoundedVector169;
const double crRightHandSideBoundedVector351 =             crRightHandSideBoundedVector108*crRightHandSideBoundedVector323;
const double crRightHandSideBoundedVector352 =             crRightHandSideBoundedVector112*crRightHandSideBoundedVector325;
const double crRightHandSideBoundedVector353 =             crRightHandSideBoundedVector111*crRightHandSideBoundedVector86;
const double crRightHandSideBoundedVector354 =             crRightHandSideBoundedVector114*crRightHandSideBoundedVector353;
const double crRightHandSideBoundedVector355 =             crRightHandSideBoundedVector114*crRightHandSideBoundedVector52;
const double crRightHandSideBoundedVector356 =             crRightHandSideBoundedVector20*crRightHandSideBoundedVector86;
const double crRightHandSideBoundedVector357 =             crRightHandSideBoundedVector111*crRightHandSideBoundedVector356;
const double crRightHandSideBoundedVector358 =             crRightHandSideBoundedVector102*crRightHandSideBoundedVector356;
const double crRightHandSideBoundedVector359 =             crRightHandSideBoundedVector20*(-crRightHandSideBoundedVector82*(0.22222222222222221*crRightHandSideBoundedVector88 + 0.22222222222222221*crRightHandSideBoundedVector90) + crRightHandSideBoundedVector94);
const double crRightHandSideBoundedVector360 =             crRightHandSideBoundedVector82*(crRightHandSideBoundedVector359 + crRightHandSideBoundedVector94);
const double crRightHandSideBoundedVector361 =             crRightHandSideBoundedVector171*crRightHandSideBoundedVector356;
const double crRightHandSideBoundedVector362 =             -0.37500000000000006*U_0_3;
const double crRightHandSideBoundedVector363 =             1.0/crRightHandSideBoundedVector85;
const double crRightHandSideBoundedVector364 =             crRightHandSideBoundedVector363*crRightHandSideBoundedVector91;
const double crRightHandSideBoundedVector365 =             -1.5000000000000002*U_1_3 + crRightHandSideBoundedVector336 + crRightHandSideBoundedVector337*crRightHandSideBoundedVector364 - 2.2500000000000004*crRightHandSideBoundedVector359 + crRightHandSideBoundedVector362;
const double crRightHandSideBoundedVector366 =             crRightHandSideBoundedVector365*crRightHandSideBoundedVector86;
const double crRightHandSideBoundedVector367 =             crRightHandSideBoundedVector177*crRightHandSideBoundedVector366;
const double crRightHandSideBoundedVector368 =             crRightHandSideBoundedVector115*crRightHandSideBoundedVector366;
const double crRightHandSideBoundedVector369 =             (crRightHandSideBoundedVector317 + crRightHandSideBoundedVector346 - crRightHandSideBoundedVector348 - crRightHandSideBoundedVector349 - crRightHandSideBoundedVector350 + crRightHandSideBoundedVector351 + crRightHandSideBoundedVector352 - crRightHandSideBoundedVector354*crRightHandSideBoundedVector70 - crRightHandSideBoundedVector355*crRightHandSideBoundedVector357 + crRightHandSideBoundedVector367 + crRightHandSideBoundedVector368 + crRightHandSideBoundedVector42*(crRightHandSideBoundedVector360 - crRightHandSideBoundedVector361) + crRightHandSideBoundedVector62*(-crRightHandSideBoundedVector358 + crRightHandSideBoundedVector360) + 0.66666666666666663*dUdt_1_3)/(crRightHandSideBoundedVector345*crRightHandSideBoundedVector82 + crRightHandSideBoundedVector96);
const double crRightHandSideBoundedVector370 =             crRightHandSideBoundedVector118*(crRightHandSideBoundedVector318 + crRightHandSideBoundedVector347 + 0.66666666666666663*r[2]);
const double crRightHandSideBoundedVector371 =             crRightHandSideBoundedVector130*crRightHandSideBoundedVector132;
const double crRightHandSideBoundedVector372 =             crRightHandSideBoundedVector135*crRightHandSideBoundedVector180;
const double crRightHandSideBoundedVector373 =             crRightHandSideBoundedVector133*crRightHandSideBoundedVector323;
const double crRightHandSideBoundedVector374 =             crRightHandSideBoundedVector136*crRightHandSideBoundedVector325;
const double crRightHandSideBoundedVector375 =             crRightHandSideBoundedVector121*crRightHandSideBoundedVector135;
const double crRightHandSideBoundedVector376 =             crRightHandSideBoundedVector138*crRightHandSideBoundedVector375;
const double crRightHandSideBoundedVector377 =             crRightHandSideBoundedVector138*crRightHandSideBoundedVector52;
const double crRightHandSideBoundedVector378 =             crRightHandSideBoundedVector121*crRightHandSideBoundedVector20;
const double crRightHandSideBoundedVector379 =             crRightHandSideBoundedVector135*crRightHandSideBoundedVector378;
const double crRightHandSideBoundedVector380 =             crRightHandSideBoundedVector141*crRightHandSideBoundedVector378;
const double crRightHandSideBoundedVector381 =             crRightHandSideBoundedVector20*(-crRightHandSideBoundedVector119*(0.22222222222222221*crRightHandSideBoundedVector122 + 0.22222222222222221*crRightHandSideBoundedVector123) + crRightHandSideBoundedVector126);
const double crRightHandSideBoundedVector382 =             crRightHandSideBoundedVector119*(crRightHandSideBoundedVector126 + crRightHandSideBoundedVector381);
const double crRightHandSideBoundedVector383 =             crRightHandSideBoundedVector187*crRightHandSideBoundedVector378;
const double crRightHandSideBoundedVector384 =             1.0/crRightHandSideBoundedVector120;
const double crRightHandSideBoundedVector385 =             crRightHandSideBoundedVector124*crRightHandSideBoundedVector384;
const double crRightHandSideBoundedVector386 =             -1.5000000000000002*U_2_3 + crRightHandSideBoundedVector335 + crRightHandSideBoundedVector337*crRightHandSideBoundedVector385 + crRightHandSideBoundedVector362 - 2.2500000000000004*crRightHandSideBoundedVector381;
const double crRightHandSideBoundedVector387 =             crRightHandSideBoundedVector121*crRightHandSideBoundedVector386;
const double crRightHandSideBoundedVector388 =             crRightHandSideBoundedVector185*crRightHandSideBoundedVector387;
const double crRightHandSideBoundedVector389 =             crRightHandSideBoundedVector139*crRightHandSideBoundedVector387;
const double crRightHandSideBoundedVector390 =             (crRightHandSideBoundedVector316 + crRightHandSideBoundedVector346 - crRightHandSideBoundedVector370 - crRightHandSideBoundedVector371 - crRightHandSideBoundedVector372 + crRightHandSideBoundedVector373 + crRightHandSideBoundedVector374 - crRightHandSideBoundedVector376*crRightHandSideBoundedVector70 - crRightHandSideBoundedVector377*crRightHandSideBoundedVector379 + crRightHandSideBoundedVector388 + crRightHandSideBoundedVector389 + crRightHandSideBoundedVector42*(crRightHandSideBoundedVector382 - crRightHandSideBoundedVector383) + crRightHandSideBoundedVector62*(-crRightHandSideBoundedVector380 + crRightHandSideBoundedVector382) + 0.66666666666666663*dUdt_2_3)/(crRightHandSideBoundedVector119*crRightHandSideBoundedVector345 + crRightHandSideBoundedVector128);
const double crRightHandSideBoundedVector391 =             pow(crRightHandSideBoundedVector10, -3);
const double crRightHandSideBoundedVector392 =             0.66666666666666663*crRightHandSideBoundedVector391;
const double crRightHandSideBoundedVector393 =             3.0000000000000004*crRightHandSideBoundedVector14;
const double crRightHandSideBoundedVector394 =             crRightHandSideBoundedVector18*crRightHandSideBoundedVector261;
const double crRightHandSideBoundedVector395 =             crRightHandSideBoundedVector33*(-crRightHandSideBoundedVector393 + crRightHandSideBoundedVector394);
const double crRightHandSideBoundedVector396 =             crRightHandSideBoundedVector260*crRightHandSideBoundedVector46;
const double crRightHandSideBoundedVector397 =             crRightHandSideBoundedVector260*crRightHandSideBoundedVector56;
const double crRightHandSideBoundedVector398 =             4.5*crRightHandSideBoundedVector391;
const double crRightHandSideBoundedVector399 =             crRightHandSideBoundedVector46*crRightHandSideBoundedVector73;
const double crRightHandSideBoundedVector400 =             crRightHandSideBoundedVector398*crRightHandSideBoundedVector399;
const double crRightHandSideBoundedVector401 =             crRightHandSideBoundedVector290*crRightHandSideBoundedVector46;
const double crRightHandSideBoundedVector402 =             crRightHandSideBoundedVector262*crRightHandSideBoundedVector56;
const double crRightHandSideBoundedVector403 =             crRightHandSideBoundedVector402*crRightHandSideBoundedVector69;
const double crRightHandSideBoundedVector404 =             crRightHandSideBoundedVector236*crRightHandSideBoundedVector56;
const double crRightHandSideBoundedVector405 =             crRightHandSideBoundedVector404*crRightHandSideBoundedVector74;
const double crRightHandSideBoundedVector406 =             0.1111111111111111*f_ext(1,0);
const double crRightHandSideBoundedVector407 =             0.1111111111111111*f_ext(2,0);
const double crRightHandSideBoundedVector408 =             crRightHandSideBoundedVector406 + crRightHandSideBoundedVector407 + 0.44444444444444442*f_ext(0,0);
const double crRightHandSideBoundedVector409 =             0.16666666666666666*dUdt_1_0;
const double crRightHandSideBoundedVector410 =             crRightHandSideBoundedVector191 + 0.16666666666666666*dUdt_2_0;
const double crRightHandSideBoundedVector411 =             crRightHandSideBoundedVector263/stab_c2;
const double crRightHandSideBoundedVector412 =             crRightHandSideBoundedVector411*(crRightHandSideBoundedVector409 + crRightHandSideBoundedVector410 + 0.66666666666666663*dUdt_0_0)/crRightHandSideBoundedVector25;
const double crRightHandSideBoundedVector413 =             crRightHandSideBoundedVector111*crRightHandSideBoundedVector243;
const double crRightHandSideBoundedVector414 =             crRightHandSideBoundedVector114*crRightHandSideBoundedVector413;
const double crRightHandSideBoundedVector415 =             0.16666666666666666*crRightHandSideBoundedVector33;
const double crRightHandSideBoundedVector416 =             pow(crRightHandSideBoundedVector85, -3);
const double crRightHandSideBoundedVector417 =             3.0000000000000004*crRightHandSideBoundedVector88;
const double crRightHandSideBoundedVector418 =             crRightHandSideBoundedVector261*crRightHandSideBoundedVector91;
const double crRightHandSideBoundedVector419 =             crRightHandSideBoundedVector416*(-crRightHandSideBoundedVector417 + crRightHandSideBoundedVector418);
const double crRightHandSideBoundedVector420 =             crRightHandSideBoundedVector107*crRightHandSideBoundedVector269;
const double crRightHandSideBoundedVector421 =             crRightHandSideBoundedVector42*crRightHandSideBoundedVector420;
const double crRightHandSideBoundedVector422 =             crRightHandSideBoundedVector111*crRightHandSideBoundedVector269;
const double crRightHandSideBoundedVector423 =             crRightHandSideBoundedVector422*crRightHandSideBoundedVector52;
const double crRightHandSideBoundedVector424 =             1.125*crRightHandSideBoundedVector416;
const double crRightHandSideBoundedVector425 =             crRightHandSideBoundedVector107*crRightHandSideBoundedVector115;
const double crRightHandSideBoundedVector426 =             crRightHandSideBoundedVector424*crRightHandSideBoundedVector425;
const double crRightHandSideBoundedVector427 =             crRightHandSideBoundedVector289*crRightHandSideBoundedVector301;
const double crRightHandSideBoundedVector428 =             crRightHandSideBoundedVector107*crRightHandSideBoundedVector86;
const double crRightHandSideBoundedVector429 =             crRightHandSideBoundedVector422*crRightHandSideBoundedVector70;
const double crRightHandSideBoundedVector430 =             0.027777777777777776*f_ext(0,0);
const double crRightHandSideBoundedVector431 =             0.027777777777777776*f_ext(2,0);
const double crRightHandSideBoundedVector432 =             crRightHandSideBoundedVector406 + crRightHandSideBoundedVector430 + crRightHandSideBoundedVector431;
const double crRightHandSideBoundedVector433 =             -crRightHandSideBoundedVector415*crRightHandSideBoundedVector419 - crRightHandSideBoundedVector421 - crRightHandSideBoundedVector423 + crRightHandSideBoundedVector426 + crRightHandSideBoundedVector427*crRightHandSideBoundedVector428 + crRightHandSideBoundedVector429 + crRightHandSideBoundedVector432;
const double crRightHandSideBoundedVector434 =             0.16666666666666666*dUdt_0_0;
const double crRightHandSideBoundedVector435 =             crRightHandSideBoundedVector411*(crRightHandSideBoundedVector410 + crRightHandSideBoundedVector434 + 0.66666666666666663*dUdt_1_0)/crRightHandSideBoundedVector95;
const double crRightHandSideBoundedVector436 =             crRightHandSideBoundedVector135*crRightHandSideBoundedVector249;
const double crRightHandSideBoundedVector437 =             crRightHandSideBoundedVector138*crRightHandSideBoundedVector436;
const double crRightHandSideBoundedVector438 =             pow(crRightHandSideBoundedVector120, -3);
const double crRightHandSideBoundedVector439 =             3.0000000000000004*crRightHandSideBoundedVector122;
const double crRightHandSideBoundedVector440 =             crRightHandSideBoundedVector124*crRightHandSideBoundedVector261;
const double crRightHandSideBoundedVector441 =             crRightHandSideBoundedVector438*(-crRightHandSideBoundedVector439 + crRightHandSideBoundedVector440);
const double crRightHandSideBoundedVector442 =             crRightHandSideBoundedVector132*crRightHandSideBoundedVector276;
const double crRightHandSideBoundedVector443 =             crRightHandSideBoundedVector42*crRightHandSideBoundedVector442;
const double crRightHandSideBoundedVector444 =             crRightHandSideBoundedVector135*crRightHandSideBoundedVector276;
const double crRightHandSideBoundedVector445 =             crRightHandSideBoundedVector444*crRightHandSideBoundedVector52;
const double crRightHandSideBoundedVector446 =             1.125*crRightHandSideBoundedVector438;
const double crRightHandSideBoundedVector447 =             crRightHandSideBoundedVector132*crRightHandSideBoundedVector139;
const double crRightHandSideBoundedVector448 =             crRightHandSideBoundedVector446*crRightHandSideBoundedVector447;
const double crRightHandSideBoundedVector449 =             crRightHandSideBoundedVector121*crRightHandSideBoundedVector132;
const double crRightHandSideBoundedVector450 =             crRightHandSideBoundedVector444*crRightHandSideBoundedVector70;
const double crRightHandSideBoundedVector451 =             0.027777777777777776*f_ext(1,0);
const double crRightHandSideBoundedVector452 =             crRightHandSideBoundedVector407 + crRightHandSideBoundedVector430 + crRightHandSideBoundedVector451;
const double crRightHandSideBoundedVector453 =             -crRightHandSideBoundedVector415*crRightHandSideBoundedVector441 + crRightHandSideBoundedVector427*crRightHandSideBoundedVector449 - crRightHandSideBoundedVector443 - crRightHandSideBoundedVector445 + crRightHandSideBoundedVector448 + crRightHandSideBoundedVector450 + crRightHandSideBoundedVector452;
const double crRightHandSideBoundedVector454 =             crRightHandSideBoundedVector411*(crRightHandSideBoundedVector191 + crRightHandSideBoundedVector409 + crRightHandSideBoundedVector434 + 0.66666666666666663*dUdt_2_0)/crRightHandSideBoundedVector127;
const double crRightHandSideBoundedVector455 =             0.16666666666666666*crRightHandSideBoundedVector131 + crRightHandSideBoundedVector132*crRightHandSideBoundedVector308 + crRightHandSideBoundedVector133*crRightHandSideBoundedVector302 - 0.16666666666666666*crRightHandSideBoundedVector134 + crRightHandSideBoundedVector136*crRightHandSideBoundedVector268 - 0.16666666666666666*crRightHandSideBoundedVector137 - 0.16666666666666666*crRightHandSideBoundedVector144;
const double crRightHandSideBoundedVector456 =             0.16666666666666666*crRightHandSideBoundedVector101 - 0.16666666666666666*crRightHandSideBoundedVector105 + crRightHandSideBoundedVector107*crRightHandSideBoundedVector298 + crRightHandSideBoundedVector108*crRightHandSideBoundedVector302 - 0.16666666666666666*crRightHandSideBoundedVector109 + crRightHandSideBoundedVector112*crRightHandSideBoundedVector268 - 0.16666666666666666*crRightHandSideBoundedVector113;
const double crRightHandSideBoundedVector457 =             crRightHandSideBoundedVector201*pow(lin_m[1], 2);
const double crRightHandSideBoundedVector458 =             crRightHandSideBoundedVector457*nu_st;
const double crRightHandSideBoundedVector459 =             1 - crRightHandSideBoundedVector457;
const double crRightHandSideBoundedVector460 =             crRightHandSideBoundedVector459*nu_sc;
const double crRightHandSideBoundedVector461 =             crRightHandSideBoundedVector195*crRightHandSideBoundedVector197 + crRightHandSideBoundedVector198 + crRightHandSideBoundedVector200*crRightHandSideBoundedVector458 + crRightHandSideBoundedVector200*crRightHandSideBoundedVector460 + crRightHandSideBoundedVector203*crRightHandSideBoundedVector205 - crRightHandSideBoundedVector205*crRightHandSideBoundedVector208;
const double crRightHandSideBoundedVector462 =             crRightHandSideBoundedVector203*crRightHandSideBoundedVector222 - crRightHandSideBoundedVector208*crRightHandSideBoundedVector222 + crRightHandSideBoundedVector215*crRightHandSideBoundedVector217 + crRightHandSideBoundedVector218 + crRightHandSideBoundedVector220*crRightHandSideBoundedVector458 + crRightHandSideBoundedVector220*crRightHandSideBoundedVector460;
const double crRightHandSideBoundedVector463 =             crRightHandSideBoundedVector203*crRightHandSideBoundedVector234 - crRightHandSideBoundedVector208*crRightHandSideBoundedVector234 + crRightHandSideBoundedVector227*crRightHandSideBoundedVector229 + crRightHandSideBoundedVector230 + crRightHandSideBoundedVector232*crRightHandSideBoundedVector458 + crRightHandSideBoundedVector232*crRightHandSideBoundedVector460;
const double crRightHandSideBoundedVector464 =             crRightHandSideBoundedVector11*crRightHandSideBoundedVector242;
const double crRightHandSideBoundedVector465 =             crRightHandSideBoundedVector248*crRightHandSideBoundedVector86;
const double crRightHandSideBoundedVector466 =             crRightHandSideBoundedVector121*crRightHandSideBoundedVector254;
const double crRightHandSideBoundedVector467 =             0.66666666666666663*crRightHandSideBoundedVector158;
const double crRightHandSideBoundedVector468 =             0.66666666666666663*crRightHandSideBoundedVector42;
const double crRightHandSideBoundedVector469 =             crRightHandSideBoundedVector282*crRightHandSideBoundedVector468;
const double crRightHandSideBoundedVector470 =             1.5000000000000002*crRightHandSideBoundedVector161;
const double crRightHandSideBoundedVector471 =             crRightHandSideBoundedVector11*crRightHandSideBoundedVector470;
const double crRightHandSideBoundedVector472 =             0.66666666666666663*crRightHandSideBoundedVector59 + 0.66666666666666663*crRightHandSideBoundedVector60 + 0.66666666666666663*crRightHandSideBoundedVector61;
const double crRightHandSideBoundedVector473 =             crRightHandSideBoundedVector288 + crRightHandSideBoundedVector4*crRightHandSideBoundedVector472 - crRightHandSideBoundedVector471;
const double crRightHandSideBoundedVector474 =             DN_DX_0_0*crRightHandSideBoundedVector53 + 0.66666666666666663*crRightHandSideBoundedVector67 + 0.66666666666666663*crRightHandSideBoundedVector68;
const double crRightHandSideBoundedVector475 =             0.16666666666666666*crRightHandSideBoundedVector59 + 0.16666666666666666*crRightHandSideBoundedVector60 + 0.16666666666666666*crRightHandSideBoundedVector61;
const double crRightHandSideBoundedVector476 =             -crRightHandSideBoundedVector177*crRightHandSideBoundedVector269 + crRightHandSideBoundedVector475*crRightHandSideBoundedVector82;
const double crRightHandSideBoundedVector477 =             crRightHandSideBoundedVector295 + crRightHandSideBoundedVector476;
const double crRightHandSideBoundedVector478 =             0.16666666666666666*crRightHandSideBoundedVector42;
const double crRightHandSideBoundedVector479 =             crRightHandSideBoundedVector282*crRightHandSideBoundedVector478;
const double crRightHandSideBoundedVector480 =             crRightHandSideBoundedVector116*crRightHandSideBoundedVector303 - crRightHandSideBoundedVector479*crRightHandSideBoundedVector82;
const double crRightHandSideBoundedVector481 =             DN_DX_1_0*crRightHandSideBoundedVector54 + DN_DX_2_0*crRightHandSideBoundedVector55 + 0.16666666666666666*crRightHandSideBoundedVector66;
const double crRightHandSideBoundedVector482 =             0.16666666666666666*crRightHandSideBoundedVector158;
const double crRightHandSideBoundedVector483 =             crRightHandSideBoundedVector20*crRightHandSideBoundedVector270 - crRightHandSideBoundedVector271 + crRightHandSideBoundedVector481*crRightHandSideBoundedVector82 - crRightHandSideBoundedVector482*crRightHandSideBoundedVector82;
const double crRightHandSideBoundedVector484 =             crRightHandSideBoundedVector119*crRightHandSideBoundedVector475 - crRightHandSideBoundedVector185*crRightHandSideBoundedVector276;
const double crRightHandSideBoundedVector485 =             crRightHandSideBoundedVector306 + crRightHandSideBoundedVector484;
const double crRightHandSideBoundedVector486 =             -crRightHandSideBoundedVector119*crRightHandSideBoundedVector479 + crRightHandSideBoundedVector140*crRightHandSideBoundedVector303;
const double crRightHandSideBoundedVector487 =             crRightHandSideBoundedVector119*crRightHandSideBoundedVector481 - crRightHandSideBoundedVector119*crRightHandSideBoundedVector482 + crRightHandSideBoundedVector20*crRightHandSideBoundedVector277 - crRightHandSideBoundedVector278;
const double crRightHandSideBoundedVector488 =             crRightHandSideBoundedVector146*crRightHandSideBoundedVector259;
const double crRightHandSideBoundedVector489 =             3.0000000000000004*crRightHandSideBoundedVector17;
const double crRightHandSideBoundedVector490 =             crRightHandSideBoundedVector72*(crRightHandSideBoundedVector394 - crRightHandSideBoundedVector489);
const double crRightHandSideBoundedVector491 =             crRightHandSideBoundedVector161*crRightHandSideBoundedVector56;
const double crRightHandSideBoundedVector492 =             crRightHandSideBoundedVector398*crRightHandSideBoundedVector491;
const double crRightHandSideBoundedVector493 =             crRightHandSideBoundedVector290*crRightHandSideBoundedVector56;
const double crRightHandSideBoundedVector494 =             crRightHandSideBoundedVector262*crRightHandSideBoundedVector46;
const double crRightHandSideBoundedVector495 =             crRightHandSideBoundedVector494*crRightHandSideBoundedVector52;
const double crRightHandSideBoundedVector496 =             0.1111111111111111*f_ext(1,1);
const double crRightHandSideBoundedVector497 =             0.1111111111111111*f_ext(2,1);
const double crRightHandSideBoundedVector498 =             crRightHandSideBoundedVector496 + crRightHandSideBoundedVector497 + 0.44444444444444442*f_ext(0,1);
const double crRightHandSideBoundedVector499 =             0.16666666666666666*crRightHandSideBoundedVector72;
const double crRightHandSideBoundedVector500 =             3.0000000000000004*crRightHandSideBoundedVector90;
const double crRightHandSideBoundedVector501 =             crRightHandSideBoundedVector416*(crRightHandSideBoundedVector418 - crRightHandSideBoundedVector500);
const double crRightHandSideBoundedVector502 =             crRightHandSideBoundedVector422*crRightHandSideBoundedVector62;
const double crRightHandSideBoundedVector503 =             crRightHandSideBoundedVector111*crRightHandSideBoundedVector177;
const double crRightHandSideBoundedVector504 =             crRightHandSideBoundedVector424*crRightHandSideBoundedVector503;
const double crRightHandSideBoundedVector505 =             crRightHandSideBoundedVector289*crRightHandSideBoundedVector478;
const double crRightHandSideBoundedVector506 =             crRightHandSideBoundedVector158*crRightHandSideBoundedVector420;
const double crRightHandSideBoundedVector507 =             0.027777777777777776*f_ext(0,1);
const double crRightHandSideBoundedVector508 =             0.027777777777777776*f_ext(2,1);
const double crRightHandSideBoundedVector509 =             crRightHandSideBoundedVector496 + crRightHandSideBoundedVector507 + crRightHandSideBoundedVector508;
const double crRightHandSideBoundedVector510 =             crRightHandSideBoundedVector353*crRightHandSideBoundedVector505 - crRightHandSideBoundedVector420*crRightHandSideBoundedVector69 - crRightHandSideBoundedVector499*crRightHandSideBoundedVector501 - crRightHandSideBoundedVector502 + crRightHandSideBoundedVector504 + crRightHandSideBoundedVector506 + crRightHandSideBoundedVector509;
const double crRightHandSideBoundedVector511 =             3.0000000000000004*crRightHandSideBoundedVector123;
const double crRightHandSideBoundedVector512 =             crRightHandSideBoundedVector438*(crRightHandSideBoundedVector440 - crRightHandSideBoundedVector511);
const double crRightHandSideBoundedVector513 =             crRightHandSideBoundedVector444*crRightHandSideBoundedVector62;
const double crRightHandSideBoundedVector514 =             crRightHandSideBoundedVector135*crRightHandSideBoundedVector185;
const double crRightHandSideBoundedVector515 =             crRightHandSideBoundedVector446*crRightHandSideBoundedVector514;
const double crRightHandSideBoundedVector516 =             crRightHandSideBoundedVector158*crRightHandSideBoundedVector442;
const double crRightHandSideBoundedVector517 =             0.027777777777777776*f_ext(1,1);
const double crRightHandSideBoundedVector518 =             crRightHandSideBoundedVector497 + crRightHandSideBoundedVector507 + crRightHandSideBoundedVector517;
const double crRightHandSideBoundedVector519 =             crRightHandSideBoundedVector375*crRightHandSideBoundedVector505 - crRightHandSideBoundedVector442*crRightHandSideBoundedVector69 - crRightHandSideBoundedVector499*crRightHandSideBoundedVector512 - crRightHandSideBoundedVector513 + crRightHandSideBoundedVector515 + crRightHandSideBoundedVector516 + crRightHandSideBoundedVector518;
const double crRightHandSideBoundedVector520 =             crRightHandSideBoundedVector133*crRightHandSideBoundedVector482 + crRightHandSideBoundedVector136*crRightHandSideBoundedVector479 + 0.16666666666666666*crRightHandSideBoundedVector181 - 0.16666666666666666*crRightHandSideBoundedVector182 - 0.16666666666666666*crRightHandSideBoundedVector183 + crRightHandSideBoundedVector185*crRightHandSideBoundedVector444 - 0.16666666666666666*crRightHandSideBoundedVector189;
const double crRightHandSideBoundedVector521 =             crRightHandSideBoundedVector108*crRightHandSideBoundedVector482 + crRightHandSideBoundedVector112*crRightHandSideBoundedVector479 + 0.16666666666666666*crRightHandSideBoundedVector170 - 0.16666666666666666*crRightHandSideBoundedVector173 - 0.16666666666666666*crRightHandSideBoundedVector174 - 0.16666666666666666*crRightHandSideBoundedVector175 + crRightHandSideBoundedVector177*crRightHandSideBoundedVector422;
const double crRightHandSideBoundedVector522 =             -crRightHandSideBoundedVector333;
const double crRightHandSideBoundedVector523 =             crRightHandSideBoundedVector331 + crRightHandSideBoundedVector522;
const double crRightHandSideBoundedVector524 =             crRightHandSideBoundedVector334 + crRightHandSideBoundedVector522;
const double crRightHandSideBoundedVector525 =             crRightHandSideBoundedVector57*mu;
const double crRightHandSideBoundedVector526 =             2.2500000000000004*crRightHandSideBoundedVector238 + 2.2500000000000004*crRightHandSideBoundedVector239 - 2.2500000000000004*crRightHandSideBoundedVector240 - 2.2500000000000004*crRightHandSideBoundedVector241;
const double crRightHandSideBoundedVector527 =             crRightHandSideBoundedVector47*mu;
const double crRightHandSideBoundedVector528 =             2.2500000000000004*crRightHandSideBoundedVector3;
const double crRightHandSideBoundedVector529 =             2.2500000000000004*crRightHandSideBoundedVector24;
const double crRightHandSideBoundedVector530 =             1.5000000000000002*crRightHandSideBoundedVector338;
const double crRightHandSideBoundedVector531 =             crRightHandSideBoundedVector33*crRightHandSideBoundedVector530;
const double crRightHandSideBoundedVector532 =             crRightHandSideBoundedVector14*crRightHandSideBoundedVector531 - crRightHandSideBoundedVector160*crRightHandSideBoundedVector69 + crRightHandSideBoundedVector17*crRightHandSideBoundedVector531 - crRightHandSideBoundedVector33*crRightHandSideBoundedVector529 + crRightHandSideBoundedVector528*crRightHandSideBoundedVector76 - crRightHandSideBoundedVector62*crRightHandSideBoundedVector74;
const double crRightHandSideBoundedVector533 =             gamma*k_st;
const double crRightHandSideBoundedVector534 =             crRightHandSideBoundedVector530*crRightHandSideBoundedVector72;
const double crRightHandSideBoundedVector535 =             crRightHandSideBoundedVector14*crRightHandSideBoundedVector534 - crRightHandSideBoundedVector160*crRightHandSideBoundedVector42 + crRightHandSideBoundedVector163*crRightHandSideBoundedVector528 + crRightHandSideBoundedVector17*crRightHandSideBoundedVector534 - crRightHandSideBoundedVector328 - crRightHandSideBoundedVector529*crRightHandSideBoundedVector72;
const double crRightHandSideBoundedVector536 =             crRightHandSideBoundedVector3*crRightHandSideBoundedVector535;
const double crRightHandSideBoundedVector537 =             crRightHandSideBoundedVector202*crRightHandSideBoundedVector536;
const double crRightHandSideBoundedVector538 =             crRightHandSideBoundedVector3*crRightHandSideBoundedVector532;
const double crRightHandSideBoundedVector539 =             crRightHandSideBoundedVector206*crRightHandSideBoundedVector533;
const double crRightHandSideBoundedVector540 =             gamma*k_sc;
const double crRightHandSideBoundedVector541 =             crRightHandSideBoundedVector209*crRightHandSideBoundedVector540;
const double crRightHandSideBoundedVector542 =             crRightHandSideBoundedVector11*(crRightHandSideBoundedVector314*crRightHandSideBoundedVector532 + crRightHandSideBoundedVector525*crRightHandSideBoundedVector526 + crRightHandSideBoundedVector527*(-3.0000000000000009*crRightHandSideBoundedVector161 + 3.0000000000000009*crRightHandSideBoundedVector192 - 1.5000000000000002*crRightHandSideBoundedVector196 + crRightHandSideBoundedVector285) + crRightHandSideBoundedVector533*crRightHandSideBoundedVector537 - crRightHandSideBoundedVector537*crRightHandSideBoundedVector540 + crRightHandSideBoundedVector538*crRightHandSideBoundedVector539 + crRightHandSideBoundedVector538*crRightHandSideBoundedVector541);
const double crRightHandSideBoundedVector543 =             crRightHandSideBoundedVector112*mu;
const double crRightHandSideBoundedVector544 =             2.2500000000000004*crRightHandSideBoundedVector244 + 2.2500000000000004*crRightHandSideBoundedVector245 - 2.2500000000000004*crRightHandSideBoundedVector246 - 2.2500000000000004*crRightHandSideBoundedVector247;
const double crRightHandSideBoundedVector545 =             1.5000000000000002*crRightHandSideBoundedVector115;
const double crRightHandSideBoundedVector546 =             crRightHandSideBoundedVector108*mu;
const double crRightHandSideBoundedVector547 =             2.2500000000000004*crRightHandSideBoundedVector81;
const double crRightHandSideBoundedVector548 =             2.2500000000000004*crRightHandSideBoundedVector94;
const double crRightHandSideBoundedVector549 =             1.5000000000000002*crRightHandSideBoundedVector363;
const double crRightHandSideBoundedVector550 =             crRightHandSideBoundedVector33*crRightHandSideBoundedVector549;
const double crRightHandSideBoundedVector551 =             -crRightHandSideBoundedVector114*crRightHandSideBoundedVector62 - crRightHandSideBoundedVector176*crRightHandSideBoundedVector69 - crRightHandSideBoundedVector33*crRightHandSideBoundedVector548 + crRightHandSideBoundedVector547*crRightHandSideBoundedVector76 + crRightHandSideBoundedVector550*crRightHandSideBoundedVector88 + crRightHandSideBoundedVector550*crRightHandSideBoundedVector90;
const double crRightHandSideBoundedVector552 =             crRightHandSideBoundedVector549*crRightHandSideBoundedVector72;
const double crRightHandSideBoundedVector553 =             crRightHandSideBoundedVector163*crRightHandSideBoundedVector547 - crRightHandSideBoundedVector176*crRightHandSideBoundedVector42 - crRightHandSideBoundedVector355 - crRightHandSideBoundedVector548*crRightHandSideBoundedVector72 + crRightHandSideBoundedVector552*crRightHandSideBoundedVector88 + crRightHandSideBoundedVector552*crRightHandSideBoundedVector90;
const double crRightHandSideBoundedVector554 =             crRightHandSideBoundedVector553*crRightHandSideBoundedVector81;
const double crRightHandSideBoundedVector555 =             crRightHandSideBoundedVector202*crRightHandSideBoundedVector554;
const double crRightHandSideBoundedVector556 =             crRightHandSideBoundedVector551*crRightHandSideBoundedVector81;
const double crRightHandSideBoundedVector557 =             crRightHandSideBoundedVector86*(crRightHandSideBoundedVector314*crRightHandSideBoundedVector551 + crRightHandSideBoundedVector533*crRightHandSideBoundedVector555 + crRightHandSideBoundedVector539*crRightHandSideBoundedVector556 - crRightHandSideBoundedVector540*crRightHandSideBoundedVector555 + crRightHandSideBoundedVector541*crRightHandSideBoundedVector556 + crRightHandSideBoundedVector543*crRightHandSideBoundedVector544 + crRightHandSideBoundedVector546*(-3.0000000000000009*crRightHandSideBoundedVector177 + 3.0000000000000009*crRightHandSideBoundedVector212 - 1.5000000000000002*crRightHandSideBoundedVector216 + crRightHandSideBoundedVector545));
const double crRightHandSideBoundedVector558 =             crRightHandSideBoundedVector136*mu;
const double crRightHandSideBoundedVector559 =             2.2500000000000004*crRightHandSideBoundedVector250 + 2.2500000000000004*crRightHandSideBoundedVector251 - 2.2500000000000004*crRightHandSideBoundedVector252 - 2.2500000000000004*crRightHandSideBoundedVector253;
const double crRightHandSideBoundedVector560 =             1.5000000000000002*crRightHandSideBoundedVector139;
const double crRightHandSideBoundedVector561 =             crRightHandSideBoundedVector133*mu;
const double crRightHandSideBoundedVector562 =             2.2500000000000004*crRightHandSideBoundedVector118;
const double crRightHandSideBoundedVector563 =             2.2500000000000004*crRightHandSideBoundedVector126;
const double crRightHandSideBoundedVector564 =             1.5000000000000002*crRightHandSideBoundedVector384;
const double crRightHandSideBoundedVector565 =             crRightHandSideBoundedVector33*crRightHandSideBoundedVector564;
const double crRightHandSideBoundedVector566 =             crRightHandSideBoundedVector122*crRightHandSideBoundedVector565 + crRightHandSideBoundedVector123*crRightHandSideBoundedVector565 - crRightHandSideBoundedVector138*crRightHandSideBoundedVector62 - crRightHandSideBoundedVector184*crRightHandSideBoundedVector69 - crRightHandSideBoundedVector33*crRightHandSideBoundedVector563 + crRightHandSideBoundedVector562*crRightHandSideBoundedVector76;
const double crRightHandSideBoundedVector567 =             crRightHandSideBoundedVector564*crRightHandSideBoundedVector72;
const double crRightHandSideBoundedVector568 =             crRightHandSideBoundedVector122*crRightHandSideBoundedVector567 + crRightHandSideBoundedVector123*crRightHandSideBoundedVector567 + crRightHandSideBoundedVector163*crRightHandSideBoundedVector562 - crRightHandSideBoundedVector184*crRightHandSideBoundedVector42 - crRightHandSideBoundedVector377 - crRightHandSideBoundedVector563*crRightHandSideBoundedVector72;
const double crRightHandSideBoundedVector569 =             crRightHandSideBoundedVector118*crRightHandSideBoundedVector568;
const double crRightHandSideBoundedVector570 =             crRightHandSideBoundedVector202*crRightHandSideBoundedVector569;
const double crRightHandSideBoundedVector571 =             crRightHandSideBoundedVector118*crRightHandSideBoundedVector566;
const double crRightHandSideBoundedVector572 =             crRightHandSideBoundedVector121*(crRightHandSideBoundedVector314*crRightHandSideBoundedVector566 + crRightHandSideBoundedVector533*crRightHandSideBoundedVector570 + crRightHandSideBoundedVector539*crRightHandSideBoundedVector571 - crRightHandSideBoundedVector540*crRightHandSideBoundedVector570 + crRightHandSideBoundedVector541*crRightHandSideBoundedVector571 + crRightHandSideBoundedVector558*crRightHandSideBoundedVector559 + crRightHandSideBoundedVector561*(-3.0000000000000009*crRightHandSideBoundedVector185 + 3.0000000000000009*crRightHandSideBoundedVector224 - 1.5000000000000002*crRightHandSideBoundedVector228 + crRightHandSideBoundedVector560));
const double crRightHandSideBoundedVector573 =             crRightHandSideBoundedVector202*crRightHandSideBoundedVector538;
const double crRightHandSideBoundedVector574 =             crRightHandSideBoundedVector457*crRightHandSideBoundedVector533;
const double crRightHandSideBoundedVector575 =             crRightHandSideBoundedVector459*crRightHandSideBoundedVector540;
const double crRightHandSideBoundedVector576 =             crRightHandSideBoundedVector11*(crRightHandSideBoundedVector314*crRightHandSideBoundedVector535 + crRightHandSideBoundedVector525*(-1.5000000000000002*crRightHandSideBoundedVector192 + 3.0000000000000009*crRightHandSideBoundedVector196 + crRightHandSideBoundedVector470 - 3.0000000000000009*crRightHandSideBoundedVector73) + crRightHandSideBoundedVector526*crRightHandSideBoundedVector527 + crRightHandSideBoundedVector533*crRightHandSideBoundedVector573 + crRightHandSideBoundedVector536*crRightHandSideBoundedVector574 + crRightHandSideBoundedVector536*crRightHandSideBoundedVector575 - crRightHandSideBoundedVector540*crRightHandSideBoundedVector573);
const double crRightHandSideBoundedVector577 =             1.5000000000000002*crRightHandSideBoundedVector177;
const double crRightHandSideBoundedVector578 =             crRightHandSideBoundedVector202*crRightHandSideBoundedVector556;
const double crRightHandSideBoundedVector579 =             crRightHandSideBoundedVector86*(crRightHandSideBoundedVector314*crRightHandSideBoundedVector553 + crRightHandSideBoundedVector533*crRightHandSideBoundedVector578 - crRightHandSideBoundedVector540*crRightHandSideBoundedVector578 + crRightHandSideBoundedVector543*(-3.0000000000000009*crRightHandSideBoundedVector115 - 1.5000000000000002*crRightHandSideBoundedVector212 + 3.0000000000000009*crRightHandSideBoundedVector216 + crRightHandSideBoundedVector577) + crRightHandSideBoundedVector544*crRightHandSideBoundedVector546 + crRightHandSideBoundedVector554*crRightHandSideBoundedVector574 + crRightHandSideBoundedVector554*crRightHandSideBoundedVector575);
const double crRightHandSideBoundedVector580 =             1.5000000000000002*crRightHandSideBoundedVector185;
const double crRightHandSideBoundedVector581 =             crRightHandSideBoundedVector202*crRightHandSideBoundedVector571;
const double crRightHandSideBoundedVector582 =             crRightHandSideBoundedVector121*(crRightHandSideBoundedVector314*crRightHandSideBoundedVector568 + crRightHandSideBoundedVector533*crRightHandSideBoundedVector581 - crRightHandSideBoundedVector540*crRightHandSideBoundedVector581 + crRightHandSideBoundedVector558*(-3.0000000000000009*crRightHandSideBoundedVector139 - 1.5000000000000002*crRightHandSideBoundedVector224 + 3.0000000000000009*crRightHandSideBoundedVector228 + crRightHandSideBoundedVector580) + crRightHandSideBoundedVector559*crRightHandSideBoundedVector561 + crRightHandSideBoundedVector569*crRightHandSideBoundedVector574 + crRightHandSideBoundedVector569*crRightHandSideBoundedVector575);
const double crRightHandSideBoundedVector583 =             crRightHandSideBoundedVector20*crRightHandSideBoundedVector338;
const double crRightHandSideBoundedVector584 =             crRightHandSideBoundedVector340 + crRightHandSideBoundedVector393*crRightHandSideBoundedVector583;
const double crRightHandSideBoundedVector585 =             0.66666666666666663*crRightHandSideBoundedVector11;
const double crRightHandSideBoundedVector586 =             4.5000000000000009*crRightHandSideBoundedVector329;
const double crRightHandSideBoundedVector587 =             crRightHandSideBoundedVector46*crRightHandSideBoundedVector62;
const double crRightHandSideBoundedVector588 =             crRightHandSideBoundedVector402*crRightHandSideBoundedVector52;
const double crRightHandSideBoundedVector589 =             crRightHandSideBoundedVector340 + crRightHandSideBoundedVector489*crRightHandSideBoundedVector583;
const double crRightHandSideBoundedVector590 =             crRightHandSideBoundedVector42*crRightHandSideBoundedVector56;
const double crRightHandSideBoundedVector591 =             crRightHandSideBoundedVector330*crRightHandSideBoundedVector74;
const double crRightHandSideBoundedVector592 =             -crRightHandSideBoundedVector360;
const double crRightHandSideBoundedVector593 =             crRightHandSideBoundedVector358 + crRightHandSideBoundedVector592;
const double crRightHandSideBoundedVector594 =             0.16666666666666666*crRightHandSideBoundedVector82;
const double crRightHandSideBoundedVector595 =             crRightHandSideBoundedVector20*crRightHandSideBoundedVector363;
const double crRightHandSideBoundedVector596 =             crRightHandSideBoundedVector86*(crRightHandSideBoundedVector365 + crRightHandSideBoundedVector417*crRightHandSideBoundedVector595);
const double crRightHandSideBoundedVector597 =             1.1250000000000002*crRightHandSideBoundedVector356;
const double crRightHandSideBoundedVector598 =             crRightHandSideBoundedVector107*crRightHandSideBoundedVector62;
const double crRightHandSideBoundedVector599 =             -crRightHandSideBoundedVector20*crRightHandSideBoundedVector421 - crRightHandSideBoundedVector20*crRightHandSideBoundedVector423 + crRightHandSideBoundedVector20*crRightHandSideBoundedVector426 + crRightHandSideBoundedVector323*crRightHandSideBoundedVector594 + crRightHandSideBoundedVector415*crRightHandSideBoundedVector596 - crRightHandSideBoundedVector429 + crRightHandSideBoundedVector432 - crRightHandSideBoundedVector597*crRightHandSideBoundedVector598;
const double crRightHandSideBoundedVector600 =             crRightHandSideBoundedVector361 + crRightHandSideBoundedVector592;
const double crRightHandSideBoundedVector601 =             crRightHandSideBoundedVector114*crRightHandSideBoundedVector357;
const double crRightHandSideBoundedVector602 =             crRightHandSideBoundedVector86*(crRightHandSideBoundedVector365 + crRightHandSideBoundedVector500*crRightHandSideBoundedVector595);
const double crRightHandSideBoundedVector603 =             crRightHandSideBoundedVector111*crRightHandSideBoundedVector42;
const double crRightHandSideBoundedVector604 =             crRightHandSideBoundedVector420*crRightHandSideBoundedVector70;
const double crRightHandSideBoundedVector605 =             -crRightHandSideBoundedVector20*crRightHandSideBoundedVector502 + crRightHandSideBoundedVector20*crRightHandSideBoundedVector504 + crRightHandSideBoundedVector325*crRightHandSideBoundedVector594 + crRightHandSideBoundedVector499*crRightHandSideBoundedVector602 - crRightHandSideBoundedVector506 + crRightHandSideBoundedVector509 - crRightHandSideBoundedVector597*crRightHandSideBoundedVector603 - crRightHandSideBoundedVector604;
const double crRightHandSideBoundedVector606 =             -crRightHandSideBoundedVector382;
const double crRightHandSideBoundedVector607 =             crRightHandSideBoundedVector380 + crRightHandSideBoundedVector606;
const double crRightHandSideBoundedVector608 =             0.16666666666666666*crRightHandSideBoundedVector119;
const double crRightHandSideBoundedVector609 =             crRightHandSideBoundedVector20*crRightHandSideBoundedVector384;
const double crRightHandSideBoundedVector610 =             crRightHandSideBoundedVector121*(crRightHandSideBoundedVector386 + crRightHandSideBoundedVector439*crRightHandSideBoundedVector609);
const double crRightHandSideBoundedVector611 =             1.1250000000000002*crRightHandSideBoundedVector378;
const double crRightHandSideBoundedVector612 =             crRightHandSideBoundedVector132*crRightHandSideBoundedVector62;
const double crRightHandSideBoundedVector613 =             -crRightHandSideBoundedVector20*crRightHandSideBoundedVector443 - crRightHandSideBoundedVector20*crRightHandSideBoundedVector445 + crRightHandSideBoundedVector20*crRightHandSideBoundedVector448 + crRightHandSideBoundedVector323*crRightHandSideBoundedVector608 + crRightHandSideBoundedVector415*crRightHandSideBoundedVector610 - crRightHandSideBoundedVector450 + crRightHandSideBoundedVector452 - crRightHandSideBoundedVector611*crRightHandSideBoundedVector612;
const double crRightHandSideBoundedVector614 =             crRightHandSideBoundedVector383 + crRightHandSideBoundedVector606;
const double crRightHandSideBoundedVector615 =             crRightHandSideBoundedVector138*crRightHandSideBoundedVector379;
const double crRightHandSideBoundedVector616 =             crRightHandSideBoundedVector121*(crRightHandSideBoundedVector386 + crRightHandSideBoundedVector511*crRightHandSideBoundedVector609);
const double crRightHandSideBoundedVector617 =             crRightHandSideBoundedVector135*crRightHandSideBoundedVector42;
const double crRightHandSideBoundedVector618 =             crRightHandSideBoundedVector442*crRightHandSideBoundedVector70;
const double crRightHandSideBoundedVector619 =             -crRightHandSideBoundedVector20*crRightHandSideBoundedVector513 + crRightHandSideBoundedVector20*crRightHandSideBoundedVector515 + crRightHandSideBoundedVector325*crRightHandSideBoundedVector608 + crRightHandSideBoundedVector499*crRightHandSideBoundedVector616 - crRightHandSideBoundedVector516 + crRightHandSideBoundedVector518 - crRightHandSideBoundedVector611*crRightHandSideBoundedVector617 - crRightHandSideBoundedVector618;
const double crRightHandSideBoundedVector620 =             crRightHandSideBoundedVector63*h;
const double crRightHandSideBoundedVector621 =             crRightHandSideBoundedVector344*crRightHandSideBoundedVector620;
const double crRightHandSideBoundedVector622 =             crRightHandSideBoundedVector369*crRightHandSideBoundedVector620;
const double crRightHandSideBoundedVector623 =             crRightHandSideBoundedVector390*crRightHandSideBoundedVector620;
const double crRightHandSideBoundedVector624 =             0.1111111111111111*r[1];
const double crRightHandSideBoundedVector625 =             0.1111111111111111*r[2];
const double crRightHandSideBoundedVector626 =             crRightHandSideBoundedVector46*crRightHandSideBoundedVector56;
const double crRightHandSideBoundedVector627 =             crRightHandSideBoundedVector398*crRightHandSideBoundedVector626;
const double crRightHandSideBoundedVector628 =             crRightHandSideBoundedVector341*crRightHandSideBoundedVector46;
const double crRightHandSideBoundedVector629 =             -1.125*U_1_3;
const double crRightHandSideBoundedVector630 =             -1.125*U_2_3;
const double crRightHandSideBoundedVector631 =             4.5000000000000009*crRightHandSideBoundedVector20;
const double crRightHandSideBoundedVector632 =             -4.5*U_0_3 - 6.7500000000000009*crRightHandSideBoundedVector332 + crRightHandSideBoundedVector339*crRightHandSideBoundedVector631 + crRightHandSideBoundedVector629 + crRightHandSideBoundedVector630;
const double crRightHandSideBoundedVector633 =             crRightHandSideBoundedVector391*crRightHandSideBoundedVector632;
const double crRightHandSideBoundedVector634 =             crRightHandSideBoundedVector11*crRightHandSideBoundedVector584;
const double crRightHandSideBoundedVector635 =             crRightHandSideBoundedVector11*crRightHandSideBoundedVector589;
const double crRightHandSideBoundedVector636 =             crRightHandSideBoundedVector107*crRightHandSideBoundedVector366;
const double crRightHandSideBoundedVector637 =             0.027777777777777776*r[0];
const double crRightHandSideBoundedVector638 =             0.027777777777777776*r[2];
const double crRightHandSideBoundedVector639 =             -1.125*U_0_3;
const double crRightHandSideBoundedVector640 =             crRightHandSideBoundedVector416*(-4.5*U_1_3 - 6.7500000000000009*crRightHandSideBoundedVector359 + crRightHandSideBoundedVector364*crRightHandSideBoundedVector631 + crRightHandSideBoundedVector630 + crRightHandSideBoundedVector639);
const double crRightHandSideBoundedVector641 =             0.16666666666666666*crRightHandSideBoundedVector640;
const double crRightHandSideBoundedVector642 =             crRightHandSideBoundedVector107*crRightHandSideBoundedVector111;
const double crRightHandSideBoundedVector643 =             crRightHandSideBoundedVector424*crRightHandSideBoundedVector642;
const double crRightHandSideBoundedVector644 =             -crRightHandSideBoundedVector115*crRightHandSideBoundedVector641 + crRightHandSideBoundedVector158*crRightHandSideBoundedVector643 - crRightHandSideBoundedVector177*crRightHandSideBoundedVector641 + crRightHandSideBoundedVector301*crRightHandSideBoundedVector596 - crRightHandSideBoundedVector323*crRightHandSideBoundedVector420 - crRightHandSideBoundedVector325*crRightHandSideBoundedVector422 + crRightHandSideBoundedVector478*crRightHandSideBoundedVector602 + crRightHandSideBoundedVector624 + crRightHandSideBoundedVector637 + crRightHandSideBoundedVector638 + crRightHandSideBoundedVector643*crRightHandSideBoundedVector70;
const double crRightHandSideBoundedVector645 =             crRightHandSideBoundedVector132*crRightHandSideBoundedVector387;
const double crRightHandSideBoundedVector646 =             0.027777777777777776*r[1];
const double crRightHandSideBoundedVector647 =             crRightHandSideBoundedVector438*(-4.5*U_2_3 - 6.7500000000000009*crRightHandSideBoundedVector381 + crRightHandSideBoundedVector385*crRightHandSideBoundedVector631 + crRightHandSideBoundedVector629 + crRightHandSideBoundedVector639);
const double crRightHandSideBoundedVector648 =             0.16666666666666666*crRightHandSideBoundedVector647;
const double crRightHandSideBoundedVector649 =             crRightHandSideBoundedVector132*crRightHandSideBoundedVector135;
const double crRightHandSideBoundedVector650 =             crRightHandSideBoundedVector446*crRightHandSideBoundedVector649;
const double crRightHandSideBoundedVector651 =             -crRightHandSideBoundedVector139*crRightHandSideBoundedVector648 + crRightHandSideBoundedVector158*crRightHandSideBoundedVector650 - crRightHandSideBoundedVector185*crRightHandSideBoundedVector648 + crRightHandSideBoundedVector301*crRightHandSideBoundedVector610 - crRightHandSideBoundedVector323*crRightHandSideBoundedVector442 - crRightHandSideBoundedVector325*crRightHandSideBoundedVector444 + crRightHandSideBoundedVector478*crRightHandSideBoundedVector616 + crRightHandSideBoundedVector625 + crRightHandSideBoundedVector637 + crRightHandSideBoundedVector646 + crRightHandSideBoundedVector650*crRightHandSideBoundedVector70;
const double crRightHandSideBoundedVector652 =             crRightHandSideBoundedVector135*crRightHandSideBoundedVector516 + crRightHandSideBoundedVector135*crRightHandSideBoundedVector618 + crRightHandSideBoundedVector301*crRightHandSideBoundedVector607 + 0.16666666666666666*crRightHandSideBoundedVector370 + 0.16666666666666666*crRightHandSideBoundedVector371 + 0.16666666666666666*crRightHandSideBoundedVector372 - 0.16666666666666666*crRightHandSideBoundedVector373 - 0.16666666666666666*crRightHandSideBoundedVector374 - 0.16666666666666666*crRightHandSideBoundedVector388 - 0.16666666666666666*crRightHandSideBoundedVector389 + crRightHandSideBoundedVector478*crRightHandSideBoundedVector614;
const double crRightHandSideBoundedVector653 =             crRightHandSideBoundedVector111*crRightHandSideBoundedVector506 + crRightHandSideBoundedVector111*crRightHandSideBoundedVector604 + crRightHandSideBoundedVector301*crRightHandSideBoundedVector593 + 0.16666666666666666*crRightHandSideBoundedVector348 + 0.16666666666666666*crRightHandSideBoundedVector349 + 0.16666666666666666*crRightHandSideBoundedVector350 - 0.16666666666666666*crRightHandSideBoundedVector351 - 0.16666666666666666*crRightHandSideBoundedVector352 - 0.16666666666666666*crRightHandSideBoundedVector367 - 0.16666666666666666*crRightHandSideBoundedVector368 + crRightHandSideBoundedVector478*crRightHandSideBoundedVector600;
const double crRightHandSideBoundedVector654 =             0.99999999999999989*crRightHandSideBoundedVector39 + 0.99999999999999989*crRightHandSideBoundedVector40 + 0.99999999999999989*crRightHandSideBoundedVector41 + 0.99999999999999989*crRightHandSideBoundedVector59 + 0.99999999999999989*crRightHandSideBoundedVector60 + 0.99999999999999989*crRightHandSideBoundedVector61;
const double crRightHandSideBoundedVector655 =             DN_DX_1_1*crRightHandSideBoundedVector47;
const double crRightHandSideBoundedVector656 =             DN_DX_1_0*crRightHandSideBoundedVector57;
const double crRightHandSideBoundedVector657 =             0.37500000000000006*crRightHandSideBoundedVector11;
const double crRightHandSideBoundedVector658 =             0.37500000000000006*crRightHandSideBoundedVector329;
const double crRightHandSideBoundedVector659 =             crRightHandSideBoundedVector240*crRightHandSideBoundedVector657 - crRightHandSideBoundedVector241*crRightHandSideBoundedVector658 - crRightHandSideBoundedVector267*crRightHandSideBoundedVector4 + crRightHandSideBoundedVector268*crRightHandSideBoundedVector4;
const double crRightHandSideBoundedVector660 =             DN_DX_1_1*crRightHandSideBoundedVector108;
const double crRightHandSideBoundedVector661 =             0.66666666666666663*crRightHandSideBoundedVector82;
const double crRightHandSideBoundedVector662 =             DN_DX_1_0*crRightHandSideBoundedVector112;
const double crRightHandSideBoundedVector663 =             1.5000000000000002*crRightHandSideBoundedVector86;
const double crRightHandSideBoundedVector664 =             crRightHandSideBoundedVector261*crRightHandSideBoundedVector86;
const double crRightHandSideBoundedVector665 =             DN_DX_1_1*crRightHandSideBoundedVector133;
const double crRightHandSideBoundedVector666 =             DN_DX_1_0*crRightHandSideBoundedVector136;
const double crRightHandSideBoundedVector667 =             crRightHandSideBoundedVector545*crRightHandSideBoundedVector86;
const double crRightHandSideBoundedVector668 =             DN_DX_1_0*crRightHandSideBoundedVector47;
const double crRightHandSideBoundedVector669 =             DN_DX_1_1*crRightHandSideBoundedVector57;
const double crRightHandSideBoundedVector670 =             crRightHandSideBoundedVector657*crRightHandSideBoundedVector73;
const double crRightHandSideBoundedVector671 =             crRightHandSideBoundedVector297*crRightHandSideBoundedVector4 - crRightHandSideBoundedVector670;
const double crRightHandSideBoundedVector672 =             crRightHandSideBoundedVector669 + crRightHandSideBoundedVector671;
const double crRightHandSideBoundedVector673 =             crRightHandSideBoundedVector162*crRightHandSideBoundedVector303 - crRightHandSideBoundedVector287*crRightHandSideBoundedVector301;
const double crRightHandSideBoundedVector674 =             crRightHandSideBoundedVector282*crRightHandSideBoundedVector82;
const double crRightHandSideBoundedVector675 =             DN_DX_1_0*crRightHandSideBoundedVector108;
const double crRightHandSideBoundedVector676 =             crRightHandSideBoundedVector289*crRightHandSideBoundedVector86;
const double crRightHandSideBoundedVector677 =             DN_DX_1_1*crRightHandSideBoundedVector112;
const double crRightHandSideBoundedVector678 =             crRightHandSideBoundedVector291*crRightHandSideBoundedVector82 - crRightHandSideBoundedVector667 + crRightHandSideBoundedVector677;
const double crRightHandSideBoundedVector679 =             DN_DX_1_0*crRightHandSideBoundedVector133;
const double crRightHandSideBoundedVector680 =             DN_DX_1_1*crRightHandSideBoundedVector136;
const double crRightHandSideBoundedVector681 =             crRightHandSideBoundedVector309 + crRightHandSideBoundedVector680;
const double crRightHandSideBoundedVector682 =             crRightHandSideBoundedVector259*h;
const double crRightHandSideBoundedVector683 =             DN_DX_1_0*crRightHandSideBoundedVector682;
const double crRightHandSideBoundedVector684 =             0.16666666666666666*crRightHandSideBoundedVector391;
const double crRightHandSideBoundedVector685 =             crRightHandSideBoundedVector46*crRightHandSideBoundedVector657;
const double crRightHandSideBoundedVector686 =             crRightHandSideBoundedVector56*crRightHandSideBoundedVector657;
const double crRightHandSideBoundedVector687 =             1.125*crRightHandSideBoundedVector391;
const double crRightHandSideBoundedVector688 =             crRightHandSideBoundedVector399*crRightHandSideBoundedVector687;
const double crRightHandSideBoundedVector689 =             crRightHandSideBoundedVector686*crRightHandSideBoundedVector70;
const double crRightHandSideBoundedVector690 =             0.1111111111111111*f_ext(0,0);
const double crRightHandSideBoundedVector691 =             crRightHandSideBoundedVector431 + crRightHandSideBoundedVector451 + crRightHandSideBoundedVector690;
const double crRightHandSideBoundedVector692 =             crRightHandSideBoundedVector301*crRightHandSideBoundedVector401 - crRightHandSideBoundedVector395*crRightHandSideBoundedVector684 - crRightHandSideBoundedVector42*crRightHandSideBoundedVector685 - crRightHandSideBoundedVector52*crRightHandSideBoundedVector686 + crRightHandSideBoundedVector688 + crRightHandSideBoundedVector689 + crRightHandSideBoundedVector691;
const double crRightHandSideBoundedVector693 =             0.66666666666666663*crRightHandSideBoundedVector33;
const double crRightHandSideBoundedVector694 =             crRightHandSideBoundedVector107*crRightHandSideBoundedVector663;
const double crRightHandSideBoundedVector695 =             crRightHandSideBoundedVector111*crRightHandSideBoundedVector663;
const double crRightHandSideBoundedVector696 =             4.5*crRightHandSideBoundedVector416;
const double crRightHandSideBoundedVector697 =             crRightHandSideBoundedVector425*crRightHandSideBoundedVector696;
const double crRightHandSideBoundedVector698 =             crRightHandSideBoundedVector283*crRightHandSideBoundedVector289;
const double crRightHandSideBoundedVector699 =             crRightHandSideBoundedVector261*crRightHandSideBoundedVector353;
const double crRightHandSideBoundedVector700 =             crRightHandSideBoundedVector69*crRightHandSideBoundedVector699;
const double crRightHandSideBoundedVector701 =             crRightHandSideBoundedVector407 + crRightHandSideBoundedVector690 + 0.44444444444444442*f_ext(1,0);
const double crRightHandSideBoundedVector702 =             crRightHandSideBoundedVector268*crRightHandSideBoundedVector57 + crRightHandSideBoundedVector302*crRightHandSideBoundedVector47 + 0.16666666666666666*crRightHandSideBoundedVector32 - 0.16666666666666666*crRightHandSideBoundedVector38 + crRightHandSideBoundedVector46*crRightHandSideBoundedVector670 - 0.16666666666666666*crRightHandSideBoundedVector48 - 0.16666666666666666*crRightHandSideBoundedVector58 - 0.99999999999999989*crRightHandSideBoundedVector77;
const double crRightHandSideBoundedVector703 =             crRightHandSideBoundedVector577*crRightHandSideBoundedVector86;
const double crRightHandSideBoundedVector704 =             -crRightHandSideBoundedVector161*crRightHandSideBoundedVector657 + crRightHandSideBoundedVector4*crRightHandSideBoundedVector475;
const double crRightHandSideBoundedVector705 =             crRightHandSideBoundedVector668 + crRightHandSideBoundedVector704;
const double crRightHandSideBoundedVector706 =             -crRightHandSideBoundedVector287*crRightHandSideBoundedVector478 + 0.16666666666666666*crRightHandSideBoundedVector290*crRightHandSideBoundedVector73;
const double crRightHandSideBoundedVector707 =             crRightHandSideBoundedVector240*crRightHandSideBoundedVector658 - crRightHandSideBoundedVector241*crRightHandSideBoundedVector657 + crRightHandSideBoundedVector4*crRightHandSideBoundedVector481 - crRightHandSideBoundedVector4*crRightHandSideBoundedVector482;
const double crRightHandSideBoundedVector708 =             crRightHandSideBoundedVector472*crRightHandSideBoundedVector82 + crRightHandSideBoundedVector675 - crRightHandSideBoundedVector703;
const double crRightHandSideBoundedVector709 =             crRightHandSideBoundedVector484 + crRightHandSideBoundedVector679;
const double crRightHandSideBoundedVector710 =             DN_DX_1_1*crRightHandSideBoundedVector682;
const double crRightHandSideBoundedVector711 =             crRightHandSideBoundedVector491*crRightHandSideBoundedVector687;
const double crRightHandSideBoundedVector712 =             crRightHandSideBoundedVector46*crRightHandSideBoundedVector658;
const double crRightHandSideBoundedVector713 =             crRightHandSideBoundedVector52*crRightHandSideBoundedVector712;
const double crRightHandSideBoundedVector714 =             0.1111111111111111*f_ext(0,1);
const double crRightHandSideBoundedVector715 =             crRightHandSideBoundedVector508 + crRightHandSideBoundedVector517 + crRightHandSideBoundedVector714;
const double crRightHandSideBoundedVector716 =             crRightHandSideBoundedVector478*crRightHandSideBoundedVector493 - crRightHandSideBoundedVector490*crRightHandSideBoundedVector684 - crRightHandSideBoundedVector62*crRightHandSideBoundedVector686 - crRightHandSideBoundedVector685*crRightHandSideBoundedVector69 + crRightHandSideBoundedVector711 + crRightHandSideBoundedVector713 + crRightHandSideBoundedVector715;
const double crRightHandSideBoundedVector717 =             0.66666666666666663*crRightHandSideBoundedVector72;
const double crRightHandSideBoundedVector718 =             crRightHandSideBoundedVector503*crRightHandSideBoundedVector696;
const double crRightHandSideBoundedVector719 =             crRightHandSideBoundedVector289*crRightHandSideBoundedVector468;
const double crRightHandSideBoundedVector720 =             crRightHandSideBoundedVector261*crRightHandSideBoundedVector428;
const double crRightHandSideBoundedVector721 =             crRightHandSideBoundedVector52*crRightHandSideBoundedVector720;
const double crRightHandSideBoundedVector722 =             crRightHandSideBoundedVector497 + crRightHandSideBoundedVector714 + 0.44444444444444442*f_ext(1,1);
const double crRightHandSideBoundedVector723 =             0.16666666666666666*crRightHandSideBoundedVector151 - 0.16666666666666666*crRightHandSideBoundedVector154 - 0.16666666666666666*crRightHandSideBoundedVector155 - 0.16666666666666666*crRightHandSideBoundedVector156 + crRightHandSideBoundedVector161*crRightHandSideBoundedVector686 - 0.99999999999999989*crRightHandSideBoundedVector164 + crRightHandSideBoundedVector47*crRightHandSideBoundedVector482 + crRightHandSideBoundedVector479*crRightHandSideBoundedVector57;
const double crRightHandSideBoundedVector724 =             0.16666666666666666*crRightHandSideBoundedVector4;
const double crRightHandSideBoundedVector725 =             1.1250000000000002*crRightHandSideBoundedVector329;
const double crRightHandSideBoundedVector726 =             crRightHandSideBoundedVector56*crRightHandSideBoundedVector658;
const double crRightHandSideBoundedVector727 =             crRightHandSideBoundedVector20*crRightHandSideBoundedVector688 + crRightHandSideBoundedVector323*crRightHandSideBoundedVector724 + crRightHandSideBoundedVector415*crRightHandSideBoundedVector634 - crRightHandSideBoundedVector42*crRightHandSideBoundedVector712 - crRightHandSideBoundedVector52*crRightHandSideBoundedVector726 - crRightHandSideBoundedVector587*crRightHandSideBoundedVector725 - crRightHandSideBoundedVector689 + crRightHandSideBoundedVector691;
const double crRightHandSideBoundedVector728 =             crRightHandSideBoundedVector685*crRightHandSideBoundedVector70;
const double crRightHandSideBoundedVector729 =             crRightHandSideBoundedVector20*crRightHandSideBoundedVector711 + crRightHandSideBoundedVector325*crRightHandSideBoundedVector724 + crRightHandSideBoundedVector499*crRightHandSideBoundedVector635 - crRightHandSideBoundedVector590*crRightHandSideBoundedVector725 - crRightHandSideBoundedVector62*crRightHandSideBoundedVector726 - crRightHandSideBoundedVector713 + crRightHandSideBoundedVector715 - crRightHandSideBoundedVector728;
const double crRightHandSideBoundedVector730 =             4.5000000000000009*crRightHandSideBoundedVector356;
const double crRightHandSideBoundedVector731 =             crRightHandSideBoundedVector69*crRightHandSideBoundedVector720;
const double crRightHandSideBoundedVector732 =             crRightHandSideBoundedVector341*crRightHandSideBoundedVector56;
const double crRightHandSideBoundedVector733 =             0.1111111111111111*r[0];
const double crRightHandSideBoundedVector734 =             crRightHandSideBoundedVector632*crRightHandSideBoundedVector684;
const double crRightHandSideBoundedVector735 =             crRightHandSideBoundedVector626*crRightHandSideBoundedVector687;
const double crRightHandSideBoundedVector736 =             crRightHandSideBoundedVector158*crRightHandSideBoundedVector735 - crRightHandSideBoundedVector161*crRightHandSideBoundedVector734 + crRightHandSideBoundedVector301*crRightHandSideBoundedVector634 - crRightHandSideBoundedVector323*crRightHandSideBoundedVector685 - crRightHandSideBoundedVector325*crRightHandSideBoundedVector686 + crRightHandSideBoundedVector478*crRightHandSideBoundedVector635 + crRightHandSideBoundedVector638 + crRightHandSideBoundedVector646 + crRightHandSideBoundedVector70*crRightHandSideBoundedVector735 - crRightHandSideBoundedVector73*crRightHandSideBoundedVector734 + crRightHandSideBoundedVector733;
const double crRightHandSideBoundedVector737 =             crRightHandSideBoundedVector642*crRightHandSideBoundedVector696;
const double crRightHandSideBoundedVector738 =             crRightHandSideBoundedVector111*crRightHandSideBoundedVector366;
const double crRightHandSideBoundedVector739 =             crRightHandSideBoundedVector135*crRightHandSideBoundedVector387;
const double crRightHandSideBoundedVector740 =             crRightHandSideBoundedVector301*crRightHandSideBoundedVector523 + 0.16666666666666666*crRightHandSideBoundedVector320 + 0.16666666666666666*crRightHandSideBoundedVector321 + 0.16666666666666666*crRightHandSideBoundedVector322 - 0.16666666666666666*crRightHandSideBoundedVector324 - 0.16666666666666666*crRightHandSideBoundedVector326 - 0.16666666666666666*crRightHandSideBoundedVector342 - 0.16666666666666666*crRightHandSideBoundedVector343 + crRightHandSideBoundedVector478*crRightHandSideBoundedVector524 + crRightHandSideBoundedVector56*crRightHandSideBoundedVector713 + crRightHandSideBoundedVector56*crRightHandSideBoundedVector728;
const double crRightHandSideBoundedVector741 =             DN_DX_2_1*crRightHandSideBoundedVector47;
const double crRightHandSideBoundedVector742 =             DN_DX_2_0*crRightHandSideBoundedVector57;
const double crRightHandSideBoundedVector743 =             DN_DX_2_1*crRightHandSideBoundedVector108;
const double crRightHandSideBoundedVector744 =             DN_DX_2_0*crRightHandSideBoundedVector112;
const double crRightHandSideBoundedVector745 =             DN_DX_2_1*crRightHandSideBoundedVector133;
const double crRightHandSideBoundedVector746 =             0.66666666666666663*crRightHandSideBoundedVector119;
const double crRightHandSideBoundedVector747 =             DN_DX_2_0*crRightHandSideBoundedVector136;
const double crRightHandSideBoundedVector748 =             1.5000000000000002*crRightHandSideBoundedVector121;
const double crRightHandSideBoundedVector749 =             crRightHandSideBoundedVector121*crRightHandSideBoundedVector261;
const double crRightHandSideBoundedVector750 =             crRightHandSideBoundedVector121*crRightHandSideBoundedVector560;
const double crRightHandSideBoundedVector751 =             DN_DX_2_0*crRightHandSideBoundedVector47;
const double crRightHandSideBoundedVector752 =             DN_DX_2_1*crRightHandSideBoundedVector57;
const double crRightHandSideBoundedVector753 =             crRightHandSideBoundedVector671 + crRightHandSideBoundedVector752;
const double crRightHandSideBoundedVector754 =             DN_DX_2_0*crRightHandSideBoundedVector108;
const double crRightHandSideBoundedVector755 =             DN_DX_2_1*crRightHandSideBoundedVector112;
const double crRightHandSideBoundedVector756 =             crRightHandSideBoundedVector299 + crRightHandSideBoundedVector755;
const double crRightHandSideBoundedVector757 =             crRightHandSideBoundedVector119*crRightHandSideBoundedVector282;
const double crRightHandSideBoundedVector758 =             DN_DX_2_0*crRightHandSideBoundedVector133;
const double crRightHandSideBoundedVector759 =             crRightHandSideBoundedVector121*crRightHandSideBoundedVector289;
const double crRightHandSideBoundedVector760 =             DN_DX_2_1*crRightHandSideBoundedVector136;
const double crRightHandSideBoundedVector761 =             crRightHandSideBoundedVector119*crRightHandSideBoundedVector291 - crRightHandSideBoundedVector750 + crRightHandSideBoundedVector760;
const double crRightHandSideBoundedVector762 =             DN_DX_2_0*crRightHandSideBoundedVector682;
const double crRightHandSideBoundedVector763 =             crRightHandSideBoundedVector132*crRightHandSideBoundedVector748;
const double crRightHandSideBoundedVector764 =             crRightHandSideBoundedVector135*crRightHandSideBoundedVector748;
const double crRightHandSideBoundedVector765 =             4.5*crRightHandSideBoundedVector438;
const double crRightHandSideBoundedVector766 =             crRightHandSideBoundedVector447*crRightHandSideBoundedVector765;
const double crRightHandSideBoundedVector767 =             crRightHandSideBoundedVector261*crRightHandSideBoundedVector375;
const double crRightHandSideBoundedVector768 =             crRightHandSideBoundedVector69*crRightHandSideBoundedVector767;
const double crRightHandSideBoundedVector769 =             crRightHandSideBoundedVector406 + crRightHandSideBoundedVector690 + 0.44444444444444442*f_ext(2,0);
const double crRightHandSideBoundedVector770 =             crRightHandSideBoundedVector121*crRightHandSideBoundedVector580;
const double crRightHandSideBoundedVector771 =             crRightHandSideBoundedVector704 + crRightHandSideBoundedVector751;
const double crRightHandSideBoundedVector772 =             crRightHandSideBoundedVector476 + crRightHandSideBoundedVector754;
const double crRightHandSideBoundedVector773 =             crRightHandSideBoundedVector119*crRightHandSideBoundedVector472 + crRightHandSideBoundedVector758 - crRightHandSideBoundedVector770;
const double crRightHandSideBoundedVector774 =             DN_DX_2_1*crRightHandSideBoundedVector682;
const double crRightHandSideBoundedVector775 =             crRightHandSideBoundedVector514*crRightHandSideBoundedVector765;
const double crRightHandSideBoundedVector776 =             crRightHandSideBoundedVector261*crRightHandSideBoundedVector449;
const double crRightHandSideBoundedVector777 =             crRightHandSideBoundedVector52*crRightHandSideBoundedVector776;
const double crRightHandSideBoundedVector778 =             crRightHandSideBoundedVector496 + crRightHandSideBoundedVector714 + 0.44444444444444442*f_ext(2,1);
const double crRightHandSideBoundedVector779 =             4.5000000000000009*crRightHandSideBoundedVector378;
const double crRightHandSideBoundedVector780 =             crRightHandSideBoundedVector69*crRightHandSideBoundedVector776;
const double crRightHandSideBoundedVector781 =             crRightHandSideBoundedVector649*crRightHandSideBoundedVector765;
            rRightHandSideBoundedVector[0]=-1.0*crRightHandSideBoundedVector0*crRightHandSideBoundedVector117 - 1.0*crRightHandSideBoundedVector0*crRightHandSideBoundedVector145 - 1.0*crRightHandSideBoundedVector0*crRightHandSideBoundedVector79 - 1.0*crRightHandSideBoundedVector146*crRightHandSideBoundedVector166 - 1.0*crRightHandSideBoundedVector146*crRightHandSideBoundedVector179 - 1.0*crRightHandSideBoundedVector146*crRightHandSideBoundedVector190 - 1.0*crRightHandSideBoundedVector191;
            rRightHandSideBoundedVector[1]=-DN_DX_0_0*crRightHandSideBoundedVector211 - DN_DX_0_0*crRightHandSideBoundedVector223 - DN_DX_0_0*crRightHandSideBoundedVector235 - crRightHandSideBoundedVector236*crRightHandSideBoundedVector242 - crRightHandSideBoundedVector243*crRightHandSideBoundedVector248 - crRightHandSideBoundedVector249*crRightHandSideBoundedVector254 + crRightHandSideBoundedVector264*(crRightHandSideBoundedVector240*crRightHandSideBoundedVector260 - crRightHandSideBoundedVector241*crRightHandSideBoundedVector262 - crRightHandSideBoundedVector255 - crRightHandSideBoundedVector256*crRightHandSideBoundedVector4 + crRightHandSideBoundedVector257*crRightHandSideBoundedVector70 + crRightHandSideBoundedVector258*crRightHandSideBoundedVector259) + crRightHandSideBoundedVector273*(crRightHandSideBoundedVector259*crRightHandSideBoundedVector266 - crRightHandSideBoundedVector265 + crRightHandSideBoundedVector272) + crRightHandSideBoundedVector280*(crRightHandSideBoundedVector259*crRightHandSideBoundedVector275 - crRightHandSideBoundedVector274 + crRightHandSideBoundedVector279) + crRightHandSideBoundedVector281*crRightHandSideBoundedVector57 + crRightHandSideBoundedVector284*crRightHandSideBoundedVector47 + crRightHandSideBoundedVector286*crRightHandSideBoundedVector46 - crRightHandSideBoundedVector294*(crRightHandSideBoundedVector199*crRightHandSideBoundedVector290 - crRightHandSideBoundedVector282*crRightHandSideBoundedVector288 - crRightHandSideBoundedVector283*crRightHandSideBoundedVector287 + crRightHandSideBoundedVector293) - crRightHandSideBoundedVector305*(-crRightHandSideBoundedVector282*crRightHandSideBoundedVector295 + crRightHandSideBoundedVector300 + crRightHandSideBoundedVector304) - crRightHandSideBoundedVector312*(-crRightHandSideBoundedVector282*crRightHandSideBoundedVector306 + crRightHandSideBoundedVector310 + crRightHandSideBoundedVector311) - crRightHandSideBoundedVector313*crRightHandSideBoundedVector344 - crRightHandSideBoundedVector313*crRightHandSideBoundedVector369 - crRightHandSideBoundedVector313*crRightHandSideBoundedVector390 + 0.66666666666666663*crRightHandSideBoundedVector32 - 0.66666666666666663*crRightHandSideBoundedVector38 - crRightHandSideBoundedVector412*(DN_DX_0_0*crRightHandSideBoundedVector37 + crRightHandSideBoundedVector283*crRightHandSideBoundedVector401 - crRightHandSideBoundedVector392*crRightHandSideBoundedVector395 - crRightHandSideBoundedVector396*crRightHandSideBoundedVector42 - crRightHandSideBoundedVector397*crRightHandSideBoundedVector52 + crRightHandSideBoundedVector400 + crRightHandSideBoundedVector403 - crRightHandSideBoundedVector405 + crRightHandSideBoundedVector408) - crRightHandSideBoundedVector435*(DN_DX_0_0*crRightHandSideBoundedVector104 - crRightHandSideBoundedVector414 + crRightHandSideBoundedVector433) - crRightHandSideBoundedVector454*(DN_DX_0_0*crRightHandSideBoundedVector143 - crRightHandSideBoundedVector437 + crRightHandSideBoundedVector453) + crRightHandSideBoundedVector455 + crRightHandSideBoundedVector456 - 0.66666666666666663*crRightHandSideBoundedVector48 - 0.66666666666666663*crRightHandSideBoundedVector58 - 1.0*crRightHandSideBoundedVector77;
            rRightHandSideBoundedVector[2]=-DN_DX_0_0*crRightHandSideBoundedVector464 - DN_DX_0_0*crRightHandSideBoundedVector465 - DN_DX_0_0*crRightHandSideBoundedVector466 - DN_DX_0_1*crRightHandSideBoundedVector461 - DN_DX_0_1*crRightHandSideBoundedVector462 - DN_DX_0_1*crRightHandSideBoundedVector463 + 0.66666666666666663*crRightHandSideBoundedVector151 - 0.66666666666666663*crRightHandSideBoundedVector154 - 0.66666666666666663*crRightHandSideBoundedVector155 - 0.66666666666666663*crRightHandSideBoundedVector156 - 1.0*crRightHandSideBoundedVector164 - crRightHandSideBoundedVector264*(crRightHandSideBoundedVector204*crRightHandSideBoundedVector290 - crRightHandSideBoundedVector282*crRightHandSideBoundedVector292 - crRightHandSideBoundedVector287*crRightHandSideBoundedVector468 + crRightHandSideBoundedVector473) - crRightHandSideBoundedVector273*(-crRightHandSideBoundedVector282*crRightHandSideBoundedVector296 + crRightHandSideBoundedVector477 + crRightHandSideBoundedVector480) - crRightHandSideBoundedVector280*(-crRightHandSideBoundedVector282*crRightHandSideBoundedVector307 + crRightHandSideBoundedVector485 + crRightHandSideBoundedVector486) - crRightHandSideBoundedVector294*(-crRightHandSideBoundedVector158*crRightHandSideBoundedVector257 + crRightHandSideBoundedVector240*crRightHandSideBoundedVector262 - crRightHandSideBoundedVector241*crRightHandSideBoundedVector260 - crRightHandSideBoundedVector255*crRightHandSideBoundedVector259 + crRightHandSideBoundedVector258 + crRightHandSideBoundedVector4*crRightHandSideBoundedVector474) - crRightHandSideBoundedVector305*(-crRightHandSideBoundedVector259*crRightHandSideBoundedVector265 + crRightHandSideBoundedVector266 + crRightHandSideBoundedVector483) - crRightHandSideBoundedVector312*(-crRightHandSideBoundedVector259*crRightHandSideBoundedVector274 + crRightHandSideBoundedVector275 + crRightHandSideBoundedVector487) - crRightHandSideBoundedVector344*crRightHandSideBoundedVector488 - crRightHandSideBoundedVector369*crRightHandSideBoundedVector488 - crRightHandSideBoundedVector390*crRightHandSideBoundedVector488 - crRightHandSideBoundedVector412*(-DN_DX_0_0*crRightHandSideBoundedVector327 + DN_DX_0_1*crRightHandSideBoundedVector153 - crRightHandSideBoundedVector392*crRightHandSideBoundedVector490 - crRightHandSideBoundedVector396*crRightHandSideBoundedVector69 - crRightHandSideBoundedVector397*crRightHandSideBoundedVector62 + crRightHandSideBoundedVector468*crRightHandSideBoundedVector493 + crRightHandSideBoundedVector492 + crRightHandSideBoundedVector495 + crRightHandSideBoundedVector498) - crRightHandSideBoundedVector435*(-DN_DX_0_0*crRightHandSideBoundedVector354 + DN_DX_0_1*crRightHandSideBoundedVector172 + crRightHandSideBoundedVector510) - crRightHandSideBoundedVector454*(-DN_DX_0_0*crRightHandSideBoundedVector376 + DN_DX_0_1*crRightHandSideBoundedVector188 + crRightHandSideBoundedVector519) + crRightHandSideBoundedVector467*crRightHandSideBoundedVector47 + crRightHandSideBoundedVector469*crRightHandSideBoundedVector57 + crRightHandSideBoundedVector471*crRightHandSideBoundedVector56 + crRightHandSideBoundedVector520 + crRightHandSideBoundedVector521;
            rRightHandSideBoundedVector[3]=-DN_DX_0_0*crRightHandSideBoundedVector542 - DN_DX_0_0*crRightHandSideBoundedVector557 - DN_DX_0_0*crRightHandSideBoundedVector572 - DN_DX_0_1*crRightHandSideBoundedVector576 - DN_DX_0_1*crRightHandSideBoundedVector579 - DN_DX_0_1*crRightHandSideBoundedVector582 - crRightHandSideBoundedVector199*crRightHandSideBoundedVector341 - crRightHandSideBoundedVector204*crRightHandSideBoundedVector341 - crRightHandSideBoundedVector264*(-DN_DX_0_0*crRightHandSideBoundedVector591 - DN_DX_0_1*crRightHandSideBoundedVector524 + crRightHandSideBoundedVector20*crRightHandSideBoundedVector492 + crRightHandSideBoundedVector257*crRightHandSideBoundedVector325 - crRightHandSideBoundedVector402*crRightHandSideBoundedVector62 - crRightHandSideBoundedVector494*crRightHandSideBoundedVector69 - crRightHandSideBoundedVector495 + crRightHandSideBoundedVector498 + crRightHandSideBoundedVector585*crRightHandSideBoundedVector589*crRightHandSideBoundedVector72 - crRightHandSideBoundedVector586*crRightHandSideBoundedVector590) - crRightHandSideBoundedVector273*(-DN_DX_0_0*crRightHandSideBoundedVector601 - DN_DX_0_1*crRightHandSideBoundedVector600 + crRightHandSideBoundedVector605) - crRightHandSideBoundedVector280*(-DN_DX_0_0*crRightHandSideBoundedVector615 - DN_DX_0_1*crRightHandSideBoundedVector614 + crRightHandSideBoundedVector619) + crRightHandSideBoundedVector283*crRightHandSideBoundedVector523 - crRightHandSideBoundedVector294*(-DN_DX_0_0*crRightHandSideBoundedVector523 + crRightHandSideBoundedVector20*crRightHandSideBoundedVector400 - crRightHandSideBoundedVector20*crRightHandSideBoundedVector405 + crRightHandSideBoundedVector257*crRightHandSideBoundedVector323 + crRightHandSideBoundedVector33*crRightHandSideBoundedVector584*crRightHandSideBoundedVector585 - crRightHandSideBoundedVector403 + crRightHandSideBoundedVector408 - crRightHandSideBoundedVector42*crRightHandSideBoundedVector494 - crRightHandSideBoundedVector586*crRightHandSideBoundedVector587 - crRightHandSideBoundedVector588) - crRightHandSideBoundedVector305*(-DN_DX_0_0*crRightHandSideBoundedVector593 - crRightHandSideBoundedVector20*crRightHandSideBoundedVector414 + crRightHandSideBoundedVector599) - crRightHandSideBoundedVector312*(-DN_DX_0_0*crRightHandSideBoundedVector607 - crRightHandSideBoundedVector20*crRightHandSideBoundedVector437 + crRightHandSideBoundedVector613) + 0.66666666666666663*crRightHandSideBoundedVector320 + 0.66666666666666663*crRightHandSideBoundedVector321 + 0.66666666666666663*crRightHandSideBoundedVector322 - 0.66666666666666663*crRightHandSideBoundedVector324 - 0.66666666666666663*crRightHandSideBoundedVector326 + crRightHandSideBoundedVector403*crRightHandSideBoundedVector46 - crRightHandSideBoundedVector412*(DN_DX_0_0*crRightHandSideBoundedVector628 + crRightHandSideBoundedVector158*crRightHandSideBoundedVector627 - crRightHandSideBoundedVector199*crRightHandSideBoundedVector633 - crRightHandSideBoundedVector204*crRightHandSideBoundedVector633 + crRightHandSideBoundedVector283*crRightHandSideBoundedVector634 - crRightHandSideBoundedVector323*crRightHandSideBoundedVector396 - crRightHandSideBoundedVector325*crRightHandSideBoundedVector397 + crRightHandSideBoundedVector340*crRightHandSideBoundedVector404 + crRightHandSideBoundedVector468*crRightHandSideBoundedVector635 + crRightHandSideBoundedVector624 + crRightHandSideBoundedVector625 + crRightHandSideBoundedVector627*crRightHandSideBoundedVector70 + 0.44444444444444442*r[0]) - crRightHandSideBoundedVector435*(DN_DX_0_0*crRightHandSideBoundedVector636 + crRightHandSideBoundedVector365*crRightHandSideBoundedVector413 + crRightHandSideBoundedVector644) - crRightHandSideBoundedVector454*(DN_DX_0_0*crRightHandSideBoundedVector645 + crRightHandSideBoundedVector386*crRightHandSideBoundedVector436 + crRightHandSideBoundedVector651) + crRightHandSideBoundedVector46*crRightHandSideBoundedVector588 + crRightHandSideBoundedVector468*crRightHandSideBoundedVector524 - crRightHandSideBoundedVector621*(crRightHandSideBoundedVector293 + crRightHandSideBoundedVector473) - crRightHandSideBoundedVector622*(crRightHandSideBoundedVector300 + crRightHandSideBoundedVector477) - crRightHandSideBoundedVector623*(crRightHandSideBoundedVector310 + crRightHandSideBoundedVector485) + crRightHandSideBoundedVector652 + crRightHandSideBoundedVector653;
            rRightHandSideBoundedVector[4]=-DN_DX_1_0*crRightHandSideBoundedVector294 - DN_DX_1_0*crRightHandSideBoundedVector305 - DN_DX_1_0*crRightHandSideBoundedVector312 - DN_DX_1_1*crRightHandSideBoundedVector264 - DN_DX_1_1*crRightHandSideBoundedVector273 - DN_DX_1_1*crRightHandSideBoundedVector280 - crRightHandSideBoundedVector654;
            rRightHandSideBoundedVector[5]=-DN_DX_1_0*crRightHandSideBoundedVector211 - DN_DX_1_0*crRightHandSideBoundedVector223 - DN_DX_1_0*crRightHandSideBoundedVector235 - DN_DX_1_1*crRightHandSideBoundedVector464 - DN_DX_1_1*crRightHandSideBoundedVector465 - DN_DX_1_1*crRightHandSideBoundedVector466 + 0.66666666666666663*crRightHandSideBoundedVector101 - 0.66666666666666663*crRightHandSideBoundedVector105 + crRightHandSideBoundedVector107*crRightHandSideBoundedVector667 + crRightHandSideBoundedVector108*crRightHandSideBoundedVector284 - 0.66666666666666663*crRightHandSideBoundedVector109 + crRightHandSideBoundedVector112*crRightHandSideBoundedVector281 - 0.66666666666666663*crRightHandSideBoundedVector113 + crRightHandSideBoundedVector264*(crRightHandSideBoundedVector259*crRightHandSideBoundedVector656 - crRightHandSideBoundedVector655 + crRightHandSideBoundedVector659) + crRightHandSideBoundedVector273*(crRightHandSideBoundedVector246*crRightHandSideBoundedVector663 - crRightHandSideBoundedVector247*crRightHandSideBoundedVector664 - crRightHandSideBoundedVector256*crRightHandSideBoundedVector82 + crRightHandSideBoundedVector259*crRightHandSideBoundedVector662 - crRightHandSideBoundedVector660 + crRightHandSideBoundedVector661*crRightHandSideBoundedVector70) + crRightHandSideBoundedVector280*(crRightHandSideBoundedVector259*crRightHandSideBoundedVector666 + crRightHandSideBoundedVector279 - crRightHandSideBoundedVector665) - crRightHandSideBoundedVector294*(-crRightHandSideBoundedVector282*crRightHandSideBoundedVector668 + crRightHandSideBoundedVector672 + crRightHandSideBoundedVector673) - crRightHandSideBoundedVector305*(crRightHandSideBoundedVector219*crRightHandSideBoundedVector676 - crRightHandSideBoundedVector282*crRightHandSideBoundedVector675 - crRightHandSideBoundedVector283*crRightHandSideBoundedVector674 + crRightHandSideBoundedVector678) - crRightHandSideBoundedVector312*(-crRightHandSideBoundedVector282*crRightHandSideBoundedVector679 + crRightHandSideBoundedVector311 + crRightHandSideBoundedVector681) - crRightHandSideBoundedVector344*crRightHandSideBoundedVector683 - crRightHandSideBoundedVector369*crRightHandSideBoundedVector683 - crRightHandSideBoundedVector390*crRightHandSideBoundedVector683 - crRightHandSideBoundedVector412*(DN_DX_1_0*crRightHandSideBoundedVector37 - DN_DX_1_1*crRightHandSideBoundedVector327 + crRightHandSideBoundedVector692) - crRightHandSideBoundedVector435*(DN_DX_1_0*crRightHandSideBoundedVector104 - DN_DX_1_1*crRightHandSideBoundedVector354 - crRightHandSideBoundedVector419*crRightHandSideBoundedVector693 - crRightHandSideBoundedVector42*crRightHandSideBoundedVector694 + crRightHandSideBoundedVector428*crRightHandSideBoundedVector698 - crRightHandSideBoundedVector52*crRightHandSideBoundedVector695 + crRightHandSideBoundedVector697 + crRightHandSideBoundedVector700 + crRightHandSideBoundedVector701) - crRightHandSideBoundedVector454*(DN_DX_1_0*crRightHandSideBoundedVector143 - DN_DX_1_1*crRightHandSideBoundedVector376 + crRightHandSideBoundedVector453) + crRightHandSideBoundedVector455 + crRightHandSideBoundedVector702;
            rRightHandSideBoundedVector[6]=-DN_DX_1_0*crRightHandSideBoundedVector464 - DN_DX_1_0*crRightHandSideBoundedVector465 - DN_DX_1_0*crRightHandSideBoundedVector466 - DN_DX_1_1*crRightHandSideBoundedVector461 - DN_DX_1_1*crRightHandSideBoundedVector462 - DN_DX_1_1*crRightHandSideBoundedVector463 + crRightHandSideBoundedVector108*crRightHandSideBoundedVector467 + crRightHandSideBoundedVector111*crRightHandSideBoundedVector703 + crRightHandSideBoundedVector112*crRightHandSideBoundedVector469 + 0.66666666666666663*crRightHandSideBoundedVector170 - 0.66666666666666663*crRightHandSideBoundedVector173 - 0.66666666666666663*crRightHandSideBoundedVector174 - 0.66666666666666663*crRightHandSideBoundedVector175 - crRightHandSideBoundedVector264*(-crRightHandSideBoundedVector282*crRightHandSideBoundedVector669 + crRightHandSideBoundedVector705 + crRightHandSideBoundedVector706) - crRightHandSideBoundedVector273*(crRightHandSideBoundedVector221*crRightHandSideBoundedVector676 - crRightHandSideBoundedVector282*crRightHandSideBoundedVector677 - crRightHandSideBoundedVector468*crRightHandSideBoundedVector674 + crRightHandSideBoundedVector708) - crRightHandSideBoundedVector280*(-crRightHandSideBoundedVector282*crRightHandSideBoundedVector680 + crRightHandSideBoundedVector486 + crRightHandSideBoundedVector709) - crRightHandSideBoundedVector294*(-crRightHandSideBoundedVector259*crRightHandSideBoundedVector655 + crRightHandSideBoundedVector656 + crRightHandSideBoundedVector707) - crRightHandSideBoundedVector305*(-crRightHandSideBoundedVector158*crRightHandSideBoundedVector661 + crRightHandSideBoundedVector246*crRightHandSideBoundedVector664 - crRightHandSideBoundedVector247*crRightHandSideBoundedVector663 - crRightHandSideBoundedVector259*crRightHandSideBoundedVector660 + crRightHandSideBoundedVector474*crRightHandSideBoundedVector82 + crRightHandSideBoundedVector662) - crRightHandSideBoundedVector312*(-crRightHandSideBoundedVector259*crRightHandSideBoundedVector665 + crRightHandSideBoundedVector487 + crRightHandSideBoundedVector666) - crRightHandSideBoundedVector344*crRightHandSideBoundedVector710 - crRightHandSideBoundedVector369*crRightHandSideBoundedVector710 - crRightHandSideBoundedVector390*crRightHandSideBoundedVector710 - crRightHandSideBoundedVector412*(-DN_DX_1_0*crRightHandSideBoundedVector327 + DN_DX_1_1*crRightHandSideBoundedVector153 + crRightHandSideBoundedVector716) - crRightHandSideBoundedVector435*(-DN_DX_1_0*crRightHandSideBoundedVector354 + DN_DX_1_1*crRightHandSideBoundedVector172 + crRightHandSideBoundedVector353*crRightHandSideBoundedVector719 - crRightHandSideBoundedVector501*crRightHandSideBoundedVector717 - crRightHandSideBoundedVector62*crRightHandSideBoundedVector695 - crRightHandSideBoundedVector69*crRightHandSideBoundedVector694 + crRightHandSideBoundedVector718 + crRightHandSideBoundedVector721 + crRightHandSideBoundedVector722) - crRightHandSideBoundedVector454*(-DN_DX_1_0*crRightHandSideBoundedVector376 + DN_DX_1_1*crRightHandSideBoundedVector188 + crRightHandSideBoundedVector519) + crRightHandSideBoundedVector520 + crRightHandSideBoundedVector723;
            rRightHandSideBoundedVector[7]=-DN_DX_1_0*crRightHandSideBoundedVector542 - DN_DX_1_0*crRightHandSideBoundedVector557 - DN_DX_1_0*crRightHandSideBoundedVector572 - DN_DX_1_1*crRightHandSideBoundedVector576 - DN_DX_1_1*crRightHandSideBoundedVector579 - DN_DX_1_1*crRightHandSideBoundedVector582 + crRightHandSideBoundedVector111*crRightHandSideBoundedVector721 + crRightHandSideBoundedVector111*crRightHandSideBoundedVector731 - crRightHandSideBoundedVector219*crRightHandSideBoundedVector366 - crRightHandSideBoundedVector221*crRightHandSideBoundedVector366 - crRightHandSideBoundedVector264*(-DN_DX_1_0*crRightHandSideBoundedVector591 - DN_DX_1_1*crRightHandSideBoundedVector524 + crRightHandSideBoundedVector729) - crRightHandSideBoundedVector273*(-DN_DX_1_0*crRightHandSideBoundedVector601 - DN_DX_1_1*crRightHandSideBoundedVector600 + crRightHandSideBoundedVector20*crRightHandSideBoundedVector718 + crRightHandSideBoundedVector325*crRightHandSideBoundedVector661 + crRightHandSideBoundedVector602*crRightHandSideBoundedVector717 - crRightHandSideBoundedVector603*crRightHandSideBoundedVector730 - crRightHandSideBoundedVector62*crRightHandSideBoundedVector699 - crRightHandSideBoundedVector721 + crRightHandSideBoundedVector722 - crRightHandSideBoundedVector731) - crRightHandSideBoundedVector280*(-DN_DX_1_0*crRightHandSideBoundedVector615 - DN_DX_1_1*crRightHandSideBoundedVector614 + crRightHandSideBoundedVector619) + crRightHandSideBoundedVector283*crRightHandSideBoundedVector593 - crRightHandSideBoundedVector294*(-DN_DX_1_0*crRightHandSideBoundedVector523 - DN_DX_1_1*crRightHandSideBoundedVector591 + crRightHandSideBoundedVector727) - crRightHandSideBoundedVector305*(-DN_DX_1_0*crRightHandSideBoundedVector593 - DN_DX_1_1*crRightHandSideBoundedVector601 + crRightHandSideBoundedVector20*crRightHandSideBoundedVector697 + crRightHandSideBoundedVector323*crRightHandSideBoundedVector661 - crRightHandSideBoundedVector42*crRightHandSideBoundedVector720 - crRightHandSideBoundedVector52*crRightHandSideBoundedVector699 + crRightHandSideBoundedVector596*crRightHandSideBoundedVector693 - crRightHandSideBoundedVector598*crRightHandSideBoundedVector730 - crRightHandSideBoundedVector700 + crRightHandSideBoundedVector701) - crRightHandSideBoundedVector312*(-DN_DX_1_0*crRightHandSideBoundedVector607 - DN_DX_1_1*crRightHandSideBoundedVector615 + crRightHandSideBoundedVector613) + 0.66666666666666663*crRightHandSideBoundedVector348 + 0.66666666666666663*crRightHandSideBoundedVector349 + 0.66666666666666663*crRightHandSideBoundedVector350 - 0.66666666666666663*crRightHandSideBoundedVector351 - 0.66666666666666663*crRightHandSideBoundedVector352 - crRightHandSideBoundedVector412*(DN_DX_1_0*crRightHandSideBoundedVector628 + DN_DX_1_1*crRightHandSideBoundedVector732 + crRightHandSideBoundedVector736) - crRightHandSideBoundedVector435*(DN_DX_1_0*crRightHandSideBoundedVector636 + DN_DX_1_1*crRightHandSideBoundedVector738 + crRightHandSideBoundedVector158*crRightHandSideBoundedVector737 - crRightHandSideBoundedVector219*crRightHandSideBoundedVector640 - crRightHandSideBoundedVector221*crRightHandSideBoundedVector640 + crRightHandSideBoundedVector283*crRightHandSideBoundedVector596 - crRightHandSideBoundedVector323*crRightHandSideBoundedVector694 - crRightHandSideBoundedVector325*crRightHandSideBoundedVector695 + crRightHandSideBoundedVector468*crRightHandSideBoundedVector602 + crRightHandSideBoundedVector625 + crRightHandSideBoundedVector70*crRightHandSideBoundedVector737 + crRightHandSideBoundedVector733 + 0.44444444444444442*r[1]) - crRightHandSideBoundedVector454*(DN_DX_1_0*crRightHandSideBoundedVector645 + DN_DX_1_1*crRightHandSideBoundedVector739 + crRightHandSideBoundedVector651) + crRightHandSideBoundedVector468*crRightHandSideBoundedVector600 - crRightHandSideBoundedVector621*(crRightHandSideBoundedVector672 + crRightHandSideBoundedVector705) - crRightHandSideBoundedVector622*(crRightHandSideBoundedVector678 + crRightHandSideBoundedVector708) - crRightHandSideBoundedVector623*(crRightHandSideBoundedVector681 + crRightHandSideBoundedVector709) + crRightHandSideBoundedVector652 + crRightHandSideBoundedVector740;
            rRightHandSideBoundedVector[8]=-DN_DX_2_0*crRightHandSideBoundedVector294 - DN_DX_2_0*crRightHandSideBoundedVector305 - DN_DX_2_0*crRightHandSideBoundedVector312 - DN_DX_2_1*crRightHandSideBoundedVector264 - DN_DX_2_1*crRightHandSideBoundedVector273 - DN_DX_2_1*crRightHandSideBoundedVector280 - crRightHandSideBoundedVector654;
            rRightHandSideBoundedVector[9]=-DN_DX_2_0*crRightHandSideBoundedVector211 - DN_DX_2_0*crRightHandSideBoundedVector223 - DN_DX_2_0*crRightHandSideBoundedVector235 - DN_DX_2_1*crRightHandSideBoundedVector464 - DN_DX_2_1*crRightHandSideBoundedVector465 - DN_DX_2_1*crRightHandSideBoundedVector466 + 0.66666666666666663*crRightHandSideBoundedVector131 + crRightHandSideBoundedVector132*crRightHandSideBoundedVector750 + crRightHandSideBoundedVector133*crRightHandSideBoundedVector284 - 0.66666666666666663*crRightHandSideBoundedVector134 + crRightHandSideBoundedVector136*crRightHandSideBoundedVector281 - 0.66666666666666663*crRightHandSideBoundedVector137 - 0.66666666666666663*crRightHandSideBoundedVector144 + crRightHandSideBoundedVector264*(crRightHandSideBoundedVector259*crRightHandSideBoundedVector742 + crRightHandSideBoundedVector659 - crRightHandSideBoundedVector741) + crRightHandSideBoundedVector273*(crRightHandSideBoundedVector259*crRightHandSideBoundedVector744 + crRightHandSideBoundedVector272 - crRightHandSideBoundedVector743) + crRightHandSideBoundedVector280*(-crRightHandSideBoundedVector119*crRightHandSideBoundedVector256 + crRightHandSideBoundedVector252*crRightHandSideBoundedVector748 - crRightHandSideBoundedVector253*crRightHandSideBoundedVector749 + crRightHandSideBoundedVector259*crRightHandSideBoundedVector747 + crRightHandSideBoundedVector70*crRightHandSideBoundedVector746 - crRightHandSideBoundedVector745) - crRightHandSideBoundedVector294*(-crRightHandSideBoundedVector282*crRightHandSideBoundedVector751 + crRightHandSideBoundedVector673 + crRightHandSideBoundedVector753) - crRightHandSideBoundedVector305*(-crRightHandSideBoundedVector282*crRightHandSideBoundedVector754 + crRightHandSideBoundedVector304 + crRightHandSideBoundedVector756) - crRightHandSideBoundedVector312*(crRightHandSideBoundedVector231*crRightHandSideBoundedVector759 - crRightHandSideBoundedVector282*crRightHandSideBoundedVector758 - crRightHandSideBoundedVector283*crRightHandSideBoundedVector757 + crRightHandSideBoundedVector761) - crRightHandSideBoundedVector344*crRightHandSideBoundedVector762 - crRightHandSideBoundedVector369*crRightHandSideBoundedVector762 - crRightHandSideBoundedVector390*crRightHandSideBoundedVector762 - crRightHandSideBoundedVector412*(DN_DX_2_0*crRightHandSideBoundedVector37 - DN_DX_2_1*crRightHandSideBoundedVector327 + crRightHandSideBoundedVector692) - crRightHandSideBoundedVector435*(DN_DX_2_0*crRightHandSideBoundedVector104 - DN_DX_2_1*crRightHandSideBoundedVector354 + crRightHandSideBoundedVector433) - crRightHandSideBoundedVector454*(DN_DX_2_0*crRightHandSideBoundedVector143 - DN_DX_2_1*crRightHandSideBoundedVector376 - crRightHandSideBoundedVector42*crRightHandSideBoundedVector763 - crRightHandSideBoundedVector441*crRightHandSideBoundedVector693 + crRightHandSideBoundedVector449*crRightHandSideBoundedVector698 - crRightHandSideBoundedVector52*crRightHandSideBoundedVector764 + crRightHandSideBoundedVector766 + crRightHandSideBoundedVector768 + crRightHandSideBoundedVector769) + crRightHandSideBoundedVector456 + crRightHandSideBoundedVector702;
            rRightHandSideBoundedVector[10]=-DN_DX_2_0*crRightHandSideBoundedVector464 - DN_DX_2_0*crRightHandSideBoundedVector465 - DN_DX_2_0*crRightHandSideBoundedVector466 - DN_DX_2_1*crRightHandSideBoundedVector461 - DN_DX_2_1*crRightHandSideBoundedVector462 - DN_DX_2_1*crRightHandSideBoundedVector463 + crRightHandSideBoundedVector133*crRightHandSideBoundedVector467 + crRightHandSideBoundedVector135*crRightHandSideBoundedVector770 + crRightHandSideBoundedVector136*crRightHandSideBoundedVector469 + 0.66666666666666663*crRightHandSideBoundedVector181 - 0.66666666666666663*crRightHandSideBoundedVector182 - 0.66666666666666663*crRightHandSideBoundedVector183 - 0.66666666666666663*crRightHandSideBoundedVector189 - crRightHandSideBoundedVector264*(-crRightHandSideBoundedVector282*crRightHandSideBoundedVector752 + crRightHandSideBoundedVector706 + crRightHandSideBoundedVector771) - crRightHandSideBoundedVector273*(-crRightHandSideBoundedVector282*crRightHandSideBoundedVector755 + crRightHandSideBoundedVector480 + crRightHandSideBoundedVector772) - crRightHandSideBoundedVector280*(crRightHandSideBoundedVector233*crRightHandSideBoundedVector759 - crRightHandSideBoundedVector282*crRightHandSideBoundedVector760 - crRightHandSideBoundedVector468*crRightHandSideBoundedVector757 + crRightHandSideBoundedVector773) - crRightHandSideBoundedVector294*(-crRightHandSideBoundedVector259*crRightHandSideBoundedVector741 + crRightHandSideBoundedVector707 + crRightHandSideBoundedVector742) - crRightHandSideBoundedVector305*(-crRightHandSideBoundedVector259*crRightHandSideBoundedVector743 + crRightHandSideBoundedVector483 + crRightHandSideBoundedVector744) - crRightHandSideBoundedVector312*(crRightHandSideBoundedVector119*crRightHandSideBoundedVector474 - crRightHandSideBoundedVector158*crRightHandSideBoundedVector746 + crRightHandSideBoundedVector252*crRightHandSideBoundedVector749 - crRightHandSideBoundedVector253*crRightHandSideBoundedVector748 - crRightHandSideBoundedVector259*crRightHandSideBoundedVector745 + crRightHandSideBoundedVector747) - crRightHandSideBoundedVector344*crRightHandSideBoundedVector774 - crRightHandSideBoundedVector369*crRightHandSideBoundedVector774 - crRightHandSideBoundedVector390*crRightHandSideBoundedVector774 - crRightHandSideBoundedVector412*(-DN_DX_2_0*crRightHandSideBoundedVector327 + DN_DX_2_1*crRightHandSideBoundedVector153 + crRightHandSideBoundedVector716) - crRightHandSideBoundedVector435*(-DN_DX_2_0*crRightHandSideBoundedVector354 + DN_DX_2_1*crRightHandSideBoundedVector172 + crRightHandSideBoundedVector510) - crRightHandSideBoundedVector454*(-DN_DX_2_0*crRightHandSideBoundedVector376 + DN_DX_2_1*crRightHandSideBoundedVector188 + crRightHandSideBoundedVector375*crRightHandSideBoundedVector719 - crRightHandSideBoundedVector512*crRightHandSideBoundedVector717 - crRightHandSideBoundedVector62*crRightHandSideBoundedVector764 - crRightHandSideBoundedVector69*crRightHandSideBoundedVector763 + crRightHandSideBoundedVector775 + crRightHandSideBoundedVector777 + crRightHandSideBoundedVector778) + crRightHandSideBoundedVector521 + crRightHandSideBoundedVector723;
            rRightHandSideBoundedVector[11]=-DN_DX_2_0*crRightHandSideBoundedVector542 - DN_DX_2_0*crRightHandSideBoundedVector557 - DN_DX_2_0*crRightHandSideBoundedVector572 - DN_DX_2_1*crRightHandSideBoundedVector576 - DN_DX_2_1*crRightHandSideBoundedVector579 - DN_DX_2_1*crRightHandSideBoundedVector582 + crRightHandSideBoundedVector135*crRightHandSideBoundedVector777 + crRightHandSideBoundedVector135*crRightHandSideBoundedVector780 - crRightHandSideBoundedVector231*crRightHandSideBoundedVector387 - crRightHandSideBoundedVector233*crRightHandSideBoundedVector387 - crRightHandSideBoundedVector264*(-DN_DX_2_0*crRightHandSideBoundedVector591 - DN_DX_2_1*crRightHandSideBoundedVector524 + crRightHandSideBoundedVector729) - crRightHandSideBoundedVector273*(-DN_DX_2_0*crRightHandSideBoundedVector601 - DN_DX_2_1*crRightHandSideBoundedVector600 + crRightHandSideBoundedVector605) - crRightHandSideBoundedVector280*(-DN_DX_2_0*crRightHandSideBoundedVector615 - DN_DX_2_1*crRightHandSideBoundedVector614 + crRightHandSideBoundedVector20*crRightHandSideBoundedVector775 + crRightHandSideBoundedVector325*crRightHandSideBoundedVector746 + crRightHandSideBoundedVector616*crRightHandSideBoundedVector717 - crRightHandSideBoundedVector617*crRightHandSideBoundedVector779 - crRightHandSideBoundedVector62*crRightHandSideBoundedVector767 - crRightHandSideBoundedVector777 + crRightHandSideBoundedVector778 - crRightHandSideBoundedVector780) + crRightHandSideBoundedVector283*crRightHandSideBoundedVector607 - crRightHandSideBoundedVector294*(-DN_DX_2_0*crRightHandSideBoundedVector523 - DN_DX_2_1*crRightHandSideBoundedVector591 + crRightHandSideBoundedVector727) - crRightHandSideBoundedVector305*(-DN_DX_2_0*crRightHandSideBoundedVector593 - DN_DX_2_1*crRightHandSideBoundedVector601 + crRightHandSideBoundedVector599) - crRightHandSideBoundedVector312*(-DN_DX_2_0*crRightHandSideBoundedVector607 - DN_DX_2_1*crRightHandSideBoundedVector615 + crRightHandSideBoundedVector20*crRightHandSideBoundedVector766 + crRightHandSideBoundedVector323*crRightHandSideBoundedVector746 - crRightHandSideBoundedVector42*crRightHandSideBoundedVector776 - crRightHandSideBoundedVector52*crRightHandSideBoundedVector767 + crRightHandSideBoundedVector610*crRightHandSideBoundedVector693 - crRightHandSideBoundedVector612*crRightHandSideBoundedVector779 - crRightHandSideBoundedVector768 + crRightHandSideBoundedVector769) + 0.66666666666666663*crRightHandSideBoundedVector370 + 0.66666666666666663*crRightHandSideBoundedVector371 + 0.66666666666666663*crRightHandSideBoundedVector372 - 0.66666666666666663*crRightHandSideBoundedVector373 - 0.66666666666666663*crRightHandSideBoundedVector374 - crRightHandSideBoundedVector412*(DN_DX_2_0*crRightHandSideBoundedVector628 + DN_DX_2_1*crRightHandSideBoundedVector732 + crRightHandSideBoundedVector736) - crRightHandSideBoundedVector435*(DN_DX_2_0*crRightHandSideBoundedVector636 + DN_DX_2_1*crRightHandSideBoundedVector738 + crRightHandSideBoundedVector644) - crRightHandSideBoundedVector454*(DN_DX_2_0*crRightHandSideBoundedVector645 + DN_DX_2_1*crRightHandSideBoundedVector739 + crRightHandSideBoundedVector158*crRightHandSideBoundedVector781 - crRightHandSideBoundedVector231*crRightHandSideBoundedVector647 - crRightHandSideBoundedVector233*crRightHandSideBoundedVector647 + crRightHandSideBoundedVector283*crRightHandSideBoundedVector610 - crRightHandSideBoundedVector323*crRightHandSideBoundedVector763 - crRightHandSideBoundedVector325*crRightHandSideBoundedVector764 + crRightHandSideBoundedVector468*crRightHandSideBoundedVector616 + crRightHandSideBoundedVector624 + crRightHandSideBoundedVector70*crRightHandSideBoundedVector781 + crRightHandSideBoundedVector733 + 0.44444444444444442*r[2]) + crRightHandSideBoundedVector468*crRightHandSideBoundedVector614 - crRightHandSideBoundedVector621*(crRightHandSideBoundedVector753 + crRightHandSideBoundedVector771) - crRightHandSideBoundedVector622*(crRightHandSideBoundedVector756 + crRightHandSideBoundedVector772) - crRightHandSideBoundedVector623*(crRightHandSideBoundedVector761 + crRightHandSideBoundedVector773) + crRightHandSideBoundedVector653 + crRightHandSideBoundedVector740;

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
    const array_1d<double, n_nodes> &r = data.r;
    const BoundedMatrix<double, n_nodes, 3> &f_ext = data.f_ext;
    const double mu = data.mu;
    const double lambda = data.lambda;
    const double c_v = data.c_v;
    const double gamma = data.gamma;

    // Stabilization parameters
    const double stab_c1 = 12.0;
    const double stab_c2 = 2.0;

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

        const double alpha = lambda / (rho_avg * gamma * c_v);
        const double tau_m_avg = 1.0 / ((4.0 * stab_c1 * mu / 3.0 / rho_avg / std::pow(h, 2)) + (stab_c2 * (v_avg_norm + c_avg)/ h));
        const double tau_et_avg = 1.0 / ((stab_c1 * alpha / std::pow(h, 2)) + (stab_c2 * (v_avg_norm + c_avg)/ h));

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
