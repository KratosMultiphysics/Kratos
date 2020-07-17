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

    if (rVariable == SHOCK_CAPTURING_VISCOSITY) {
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
    rData.c_v = r_properties.GetValue(SPECIFIC_HEAT);
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

    const double cmom_proj0 =             -0.25*dUdt_1_1;
const double cmom_proj1 =             0.166666666666667*U_0_0;
const double cmom_proj2 =             0.166666666666667*U_2_0;
const double cmom_proj3 =             0.666666666666667*U_1_0 + cmom_proj1 + cmom_proj2;
const double cmom_proj4 =             0.166666666666667*f_ext(0,0);
const double cmom_proj5 =             0.166666666666667*f_ext(2,0);
const double cmom_proj6 =             cmom_proj3*(cmom_proj4 + cmom_proj5 + 0.666666666666667*f_ext(1,0));
const double cmom_proj7 =             0.166666666666667*cmom_proj6;
const double cmom_proj8 =             0.166666666666667*U_0_1;
const double cmom_proj9 =             0.166666666666667*U_2_1;
const double cmom_proj10 =             0.666666666666667*U_1_1 + cmom_proj8 + cmom_proj9;
const double cmom_proj11 =             DN_DX_0_1*U_0_2 + DN_DX_1_1*U_1_2 + DN_DX_2_1*U_2_2;
const double cmom_proj12 =             1.0/cmom_proj3;
const double cmom_proj13 =             cmom_proj10*cmom_proj11*cmom_proj12;
const double cmom_proj14 =             -0.166666666666667*cmom_proj13;
const double cmom_proj15 =             0.166666666666667*U_0_2;
const double cmom_proj16 =             0.166666666666667*U_2_2;
const double cmom_proj17 =             0.666666666666667*U_1_2 + cmom_proj15 + cmom_proj16;
const double cmom_proj18 =             DN_DX_0_1*U_0_1 + DN_DX_1_1*U_1_1 + DN_DX_2_1*U_2_1;
const double cmom_proj19 =             cmom_proj12*cmom_proj17*cmom_proj18;
const double cmom_proj20 =             -0.166666666666667*cmom_proj19;
const double cmom_proj21 =             gamma - 1;
const double cmom_proj22 =             DN_DX_0_0*U_0_2 + DN_DX_1_0*U_1_2 + DN_DX_2_0*U_2_2;
const double cmom_proj23 =             cmom_proj12*cmom_proj17*cmom_proj21*cmom_proj22;
const double cmom_proj24 =             0.166666666666667*cmom_proj23;
const double cmom_proj25 =             1.0*gamma - 3.0;
const double cmom_proj26 =             DN_DX_0_0*U_0_1 + DN_DX_1_0*U_1_1 + DN_DX_2_0*U_2_1;
const double cmom_proj27 =             cmom_proj10*cmom_proj12*cmom_proj25*cmom_proj26;
const double cmom_proj28 =             0.166666666666667*cmom_proj27;
const double cmom_proj29 =             DN_DX_0_1*U_0_0 + DN_DX_1_1*U_1_0 + DN_DX_2_1*U_2_0;
const double cmom_proj30 =             pow(cmom_proj3, -2);
const double cmom_proj31 =             cmom_proj10*cmom_proj17*cmom_proj29*cmom_proj30;
const double cmom_proj32 =             0.166666666666667*cmom_proj31;
const double cmom_proj33 =             DN_DX_0_0*U_0_0 + DN_DX_1_0*U_1_0 + DN_DX_2_0*U_2_0;
const double cmom_proj34 =             pow(cmom_proj10, 2);
const double cmom_proj35 =             0.5*gamma - 0.5;
const double cmom_proj36 =             pow(cmom_proj17, 2);
const double cmom_proj37 =             cmom_proj35*(cmom_proj34 + cmom_proj36);
const double cmom_proj38 =             cmom_proj30*cmom_proj33*(-cmom_proj34 + cmom_proj37);
const double cmom_proj39 =             -0.166666666666667*cmom_proj38;
const double cmom_proj40 =             -0.25*dUdt_2_1;
const double cmom_proj41 =             cmom_proj21*(DN_DX_0_0*U_0_3 + DN_DX_1_0*U_1_3 + DN_DX_2_0*U_2_3);
const double cmom_proj42 =             -1.0*cmom_proj41;
const double cmom_proj43 =             0.166666666666667*U_1_0;
const double cmom_proj44 =             0.666666666666667*U_2_0 + cmom_proj1 + cmom_proj43;
const double cmom_proj45 =             0.166666666666667*f_ext(1,0);
const double cmom_proj46 =             cmom_proj44*(cmom_proj4 + cmom_proj45 + 0.666666666666667*f_ext(2,0));
const double cmom_proj47 =             0.166666666666667*cmom_proj46;
const double cmom_proj48 =             0.666666666666667*U_0_0 + cmom_proj2 + cmom_proj43;
const double cmom_proj49 =             cmom_proj48*(cmom_proj45 + cmom_proj5 + 0.666666666666667*f_ext(0,0));
const double cmom_proj50 =             0.166666666666667*U_1_1;
const double cmom_proj51 =             0.666666666666667*U_2_1 + cmom_proj50 + cmom_proj8;
const double cmom_proj52 =             1.0/cmom_proj44;
const double cmom_proj53 =             cmom_proj11*cmom_proj51*cmom_proj52;
const double cmom_proj54 =             -0.166666666666667*cmom_proj53;
const double cmom_proj55 =             0.166666666666667*U_1_2;
const double cmom_proj56 =             0.666666666666667*U_2_2 + cmom_proj15 + cmom_proj55;
const double cmom_proj57 =             cmom_proj18*cmom_proj52*cmom_proj56;
const double cmom_proj58 =             -0.166666666666667*cmom_proj57;
const double cmom_proj59 =             0.666666666666667*U_0_1 + cmom_proj50 + cmom_proj9;
const double cmom_proj60 =             1.0/cmom_proj48;
const double cmom_proj61 =             cmom_proj11*cmom_proj59*cmom_proj60;
const double cmom_proj62 =             0.666666666666667*U_0_2 + cmom_proj16 + cmom_proj55;
const double cmom_proj63 =             cmom_proj18*cmom_proj60*cmom_proj62;
const double cmom_proj64 =             cmom_proj21*cmom_proj22*cmom_proj52*cmom_proj56;
const double cmom_proj65 =             0.166666666666667*cmom_proj64;
const double cmom_proj66 =             cmom_proj21*cmom_proj22*cmom_proj60*cmom_proj62;
const double cmom_proj67 =             cmom_proj25*cmom_proj26*cmom_proj51*cmom_proj52;
const double cmom_proj68 =             0.166666666666667*cmom_proj67;
const double cmom_proj69 =             cmom_proj25*cmom_proj26*cmom_proj59*cmom_proj60;
const double cmom_proj70 =             pow(cmom_proj44, -2);
const double cmom_proj71 =             cmom_proj29*cmom_proj51*cmom_proj56*cmom_proj70;
const double cmom_proj72 =             0.166666666666667*cmom_proj71;
const double cmom_proj73 =             pow(cmom_proj48, -2);
const double cmom_proj74 =             cmom_proj29*cmom_proj59*cmom_proj62*cmom_proj73;
const double cmom_proj75 =             pow(cmom_proj51, 2);
const double cmom_proj76 =             pow(cmom_proj56, 2);
const double cmom_proj77 =             cmom_proj35*(cmom_proj75 + cmom_proj76);
const double cmom_proj78 =             cmom_proj33*cmom_proj70*(-cmom_proj75 + cmom_proj77);
const double cmom_proj79 =             -0.166666666666667*cmom_proj78;
const double cmom_proj80 =             pow(cmom_proj59, 2);
const double cmom_proj81 =             pow(cmom_proj62, 2);
const double cmom_proj82 =             cmom_proj35*(cmom_proj80 + cmom_proj81);
const double cmom_proj83 =             cmom_proj33*cmom_proj73*(-cmom_proj80 + cmom_proj82);
const double cmom_proj84 =             -0.25*dUdt_1_2;
const double cmom_proj85 =             0.166666666666667*f_ext(0,1);
const double cmom_proj86 =             0.166666666666667*f_ext(2,1);
const double cmom_proj87 =             cmom_proj3*(cmom_proj85 + cmom_proj86 + 0.666666666666667*f_ext(1,1));
const double cmom_proj88 =             0.166666666666667*cmom_proj87;
const double cmom_proj89 =             cmom_proj10*cmom_proj12*cmom_proj22;
const double cmom_proj90 =             -0.166666666666667*cmom_proj89;
const double cmom_proj91 =             cmom_proj12*cmom_proj17*cmom_proj26;
const double cmom_proj92 =             -0.166666666666667*cmom_proj91;
const double cmom_proj93 =             cmom_proj10*cmom_proj12*cmom_proj18*cmom_proj21;
const double cmom_proj94 =             0.166666666666667*cmom_proj93;
const double cmom_proj95 =             cmom_proj11*cmom_proj12*cmom_proj17*cmom_proj25;
const double cmom_proj96 =             0.166666666666667*cmom_proj95;
const double cmom_proj97 =             cmom_proj10*cmom_proj17*cmom_proj30*cmom_proj33;
const double cmom_proj98 =             0.166666666666667*cmom_proj97;
const double cmom_proj99 =             cmom_proj29*cmom_proj30*(-cmom_proj36 + cmom_proj37);
const double cmom_proj100 =             -0.166666666666667*cmom_proj99;
const double cmom_proj101 =             -0.25*dUdt_2_2;
const double cmom_proj102 =             cmom_proj21*(DN_DX_0_1*U_0_3 + DN_DX_1_1*U_1_3 + DN_DX_2_1*U_2_3);
const double cmom_proj103 =             -1.0*cmom_proj102;
const double cmom_proj104 =             0.166666666666667*f_ext(1,1);
const double cmom_proj105 =             cmom_proj44*(cmom_proj104 + cmom_proj85 + 0.666666666666667*f_ext(2,1));
const double cmom_proj106 =             0.166666666666667*cmom_proj105;
const double cmom_proj107 =             cmom_proj48*(cmom_proj104 + cmom_proj86 + 0.666666666666667*f_ext(0,1));
const double cmom_proj108 =             cmom_proj22*cmom_proj51*cmom_proj52;
const double cmom_proj109 =             -0.166666666666667*cmom_proj108;
const double cmom_proj110 =             cmom_proj26*cmom_proj52*cmom_proj56;
const double cmom_proj111 =             -0.166666666666667*cmom_proj110;
const double cmom_proj112 =             cmom_proj22*cmom_proj59*cmom_proj60;
const double cmom_proj113 =             cmom_proj26*cmom_proj60*cmom_proj62;
const double cmom_proj114 =             cmom_proj18*cmom_proj21*cmom_proj51*cmom_proj52;
const double cmom_proj115 =             0.166666666666667*cmom_proj114;
const double cmom_proj116 =             cmom_proj18*cmom_proj21*cmom_proj59*cmom_proj60;
const double cmom_proj117 =             cmom_proj11*cmom_proj25*cmom_proj52*cmom_proj56;
const double cmom_proj118 =             0.166666666666667*cmom_proj117;
const double cmom_proj119 =             cmom_proj11*cmom_proj25*cmom_proj60*cmom_proj62;
const double cmom_proj120 =             cmom_proj33*cmom_proj51*cmom_proj56*cmom_proj70;
const double cmom_proj121 =             0.166666666666667*cmom_proj120;
const double cmom_proj122 =             cmom_proj33*cmom_proj59*cmom_proj62*cmom_proj73;
const double cmom_proj123 =             cmom_proj29*cmom_proj70*(-cmom_proj76 + cmom_proj77);
const double cmom_proj124 =             -0.166666666666667*cmom_proj123;
const double cmom_proj125 =             cmom_proj29*cmom_proj73*(-cmom_proj81 + cmom_proj82);
const double cmom_proj126 =             0.166666666666667*cmom_proj49 - 0.166666666666667*cmom_proj61 - 0.166666666666667*cmom_proj63 + 0.166666666666667*cmom_proj66 + 0.166666666666667*cmom_proj69 + 0.166666666666667*cmom_proj74 - 0.166666666666667*cmom_proj83 - 0.25*dUdt_0_1;
const double cmom_proj127 =             0.166666666666667*cmom_proj107 - 0.166666666666667*cmom_proj112 - 0.166666666666667*cmom_proj113 + 0.166666666666667*cmom_proj116 + 0.166666666666667*cmom_proj119 + 0.166666666666667*cmom_proj122 - 0.166666666666667*cmom_proj125 - 0.25*dUdt_0_2;
            mom_proj[0]=cmom_proj0 + cmom_proj14 + cmom_proj20 + cmom_proj24 + cmom_proj28 + cmom_proj32 + cmom_proj39 + cmom_proj40 + cmom_proj42 + cmom_proj47 + 0.666666666666667*cmom_proj49 + cmom_proj54 + cmom_proj58 - 0.666666666666667*cmom_proj61 - 0.666666666666667*cmom_proj63 + cmom_proj65 + 0.666666666666667*cmom_proj66 + cmom_proj68 + 0.666666666666667*cmom_proj69 + cmom_proj7 + cmom_proj72 + 0.666666666666667*cmom_proj74 + cmom_proj79 - 0.666666666666667*cmom_proj83 - 0.5*dUdt_0_1;
            mom_proj[1]=cmom_proj100 + cmom_proj101 + cmom_proj103 + cmom_proj106 + 0.666666666666667*cmom_proj107 + cmom_proj109 + cmom_proj111 - 0.666666666666667*cmom_proj112 - 0.666666666666667*cmom_proj113 + cmom_proj115 + 0.666666666666667*cmom_proj116 + cmom_proj118 + 0.666666666666667*cmom_proj119 + cmom_proj121 + 0.666666666666667*cmom_proj122 + cmom_proj124 - 0.666666666666667*cmom_proj125 + cmom_proj84 + cmom_proj88 + cmom_proj90 + cmom_proj92 + cmom_proj94 + cmom_proj96 + cmom_proj98 - 0.5*dUdt_0_2;
            mom_proj[2]=cmom_proj126 - 0.666666666666667*cmom_proj13 - 0.666666666666667*cmom_proj19 + 0.666666666666667*cmom_proj23 + 0.666666666666667*cmom_proj27 + 0.666666666666667*cmom_proj31 - 0.666666666666667*cmom_proj38 + cmom_proj40 + cmom_proj42 + cmom_proj47 + cmom_proj54 + cmom_proj58 + 0.666666666666667*cmom_proj6 + cmom_proj65 + cmom_proj68 + cmom_proj72 + cmom_proj79 - 0.5*dUdt_1_1;
            mom_proj[3]=cmom_proj101 + cmom_proj103 + cmom_proj106 + cmom_proj109 + cmom_proj111 + cmom_proj115 + cmom_proj118 + cmom_proj121 + cmom_proj124 + cmom_proj127 + 0.666666666666667*cmom_proj87 - 0.666666666666667*cmom_proj89 - 0.666666666666667*cmom_proj91 + 0.666666666666667*cmom_proj93 + 0.666666666666667*cmom_proj95 + 0.666666666666667*cmom_proj97 - 0.666666666666667*cmom_proj99 - 0.5*dUdt_1_2;
            mom_proj[4]=cmom_proj0 + cmom_proj126 + cmom_proj14 + cmom_proj20 + cmom_proj24 + cmom_proj28 + cmom_proj32 + cmom_proj39 - 1.0*cmom_proj41 + 0.666666666666667*cmom_proj46 - 0.666666666666667*cmom_proj53 - 0.666666666666667*cmom_proj57 + 0.666666666666667*cmom_proj64 + 0.666666666666667*cmom_proj67 + cmom_proj7 + 0.666666666666667*cmom_proj71 - 0.666666666666667*cmom_proj78 - 0.5*dUdt_2_1;
            mom_proj[5]=cmom_proj100 - 1.0*cmom_proj102 + 0.666666666666667*cmom_proj105 - 0.666666666666667*cmom_proj108 - 0.666666666666667*cmom_proj110 + 0.666666666666667*cmom_proj114 + 0.666666666666667*cmom_proj117 + 0.666666666666667*cmom_proj120 - 0.666666666666667*cmom_proj123 + cmom_proj127 + cmom_proj84 + cmom_proj88 + cmom_proj90 + cmom_proj92 + cmom_proj94 + cmom_proj96 + cmom_proj98 - 0.5*dUdt_2_2;

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
const double crho_proj7 =             1.0*crho_proj1 + 1.0*crho_proj2 + 1.0*crho_proj3 + 1.0*crho_proj4 + 1.0*crho_proj5 + 1.0*crho_proj6 + 0.25*dUdt_2_0;
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

    const double ctot_ener_proj0 =             -0.25*dUdt_1_3;
const double ctot_ener_proj1 =             0.166666666666667*U_0_0;
const double ctot_ener_proj2 =             0.166666666666667*U_2_0;
const double ctot_ener_proj3 =             0.666666666666667*U_1_0 + ctot_ener_proj1 + ctot_ener_proj2;
const double ctot_ener_proj4 =             0.166666666666667*r[0];
const double ctot_ener_proj5 =             0.166666666666667*r[2];
const double ctot_ener_proj6 =             ctot_ener_proj3*(ctot_ener_proj4 + ctot_ener_proj5 + 0.666666666666667*r[1]);
const double ctot_ener_proj7 =             0.166666666666667*ctot_ener_proj6;
const double ctot_ener_proj8 =             0.166666666666667*U_0_1;
const double ctot_ener_proj9 =             0.166666666666667*U_2_1;
const double ctot_ener_proj10 =             0.666666666666667*U_1_1 + ctot_ener_proj8 + ctot_ener_proj9;
const double ctot_ener_proj11 =             0.166666666666667*f_ext(0,0);
const double ctot_ener_proj12 =             0.166666666666667*f_ext(2,0);
const double ctot_ener_proj13 =             ctot_ener_proj10*(ctot_ener_proj11 + ctot_ener_proj12 + 0.666666666666667*f_ext(1,0));
const double ctot_ener_proj14 =             0.166666666666667*ctot_ener_proj13;
const double ctot_ener_proj15 =             0.166666666666667*U_0_2;
const double ctot_ener_proj16 =             0.166666666666667*U_2_2;
const double ctot_ener_proj17 =             0.666666666666667*U_1_2 + ctot_ener_proj15 + ctot_ener_proj16;
const double ctot_ener_proj18 =             0.166666666666667*f_ext(0,1);
const double ctot_ener_proj19 =             0.166666666666667*f_ext(2,1);
const double ctot_ener_proj20 =             ctot_ener_proj17*(ctot_ener_proj18 + ctot_ener_proj19 + 0.666666666666667*f_ext(1,1));
const double ctot_ener_proj21 =             0.166666666666667*ctot_ener_proj20;
const double ctot_ener_proj22 =             DN_DX_0_0*U_0_3 + DN_DX_1_0*U_1_3 + DN_DX_2_0*U_2_3;
const double ctot_ener_proj23 =             1.0/ctot_ener_proj3;
const double ctot_ener_proj24 =             ctot_ener_proj10*ctot_ener_proj22*ctot_ener_proj23*gamma;
const double ctot_ener_proj25 =             -0.166666666666667*ctot_ener_proj24;
const double ctot_ener_proj26 =             DN_DX_0_1*U_0_3 + DN_DX_1_1*U_1_3 + DN_DX_2_1*U_2_3;
const double ctot_ener_proj27 =             ctot_ener_proj17*ctot_ener_proj23*ctot_ener_proj26*gamma;
const double ctot_ener_proj28 =             -0.166666666666667*ctot_ener_proj27;
const double ctot_ener_proj29 =             gamma - 1;
const double ctot_ener_proj30 =             DN_DX_0_0*U_0_2 + DN_DX_1_0*U_1_2 + DN_DX_2_0*U_2_2;
const double ctot_ener_proj31 =             pow(ctot_ener_proj3, -2);
const double ctot_ener_proj32 =             ctot_ener_proj10*ctot_ener_proj17*ctot_ener_proj29*ctot_ener_proj30*ctot_ener_proj31;
const double ctot_ener_proj33 =             0.166666666666667*ctot_ener_proj32;
const double ctot_ener_proj34 =             DN_DX_0_1*U_0_1 + DN_DX_1_1*U_1_1 + DN_DX_2_1*U_2_1;
const double ctot_ener_proj35 =             ctot_ener_proj10*ctot_ener_proj17*ctot_ener_proj29*ctot_ener_proj31*ctot_ener_proj34;
const double ctot_ener_proj36 =             0.166666666666667*ctot_ener_proj35;
const double ctot_ener_proj37 =             DN_DX_0_0*U_0_1 + DN_DX_1_0*U_1_1 + DN_DX_2_0*U_2_1;
const double ctot_ener_proj38 =             pow(ctot_ener_proj10, 2);
const double ctot_ener_proj39 =             1.0*ctot_ener_proj23*ctot_ener_proj29;
const double ctot_ener_proj40 =             0.166666666666667*U_0_3;
const double ctot_ener_proj41 =             -ctot_ener_proj40;
const double ctot_ener_proj42 =             0.666666666666667*U_1_3;
const double ctot_ener_proj43 =             0.166666666666667*U_2_3;
const double ctot_ener_proj44 =             -ctot_ener_proj43;
const double ctot_ener_proj45 =             pow(ctot_ener_proj17, 2);
const double ctot_ener_proj46 =             -ctot_ener_proj29*(-ctot_ener_proj23*(0.5*ctot_ener_proj38 + 0.5*ctot_ener_proj45) + ctot_ener_proj40 + ctot_ener_proj42 + ctot_ener_proj43) + ctot_ener_proj41 - ctot_ener_proj42 + ctot_ener_proj44;
const double ctot_ener_proj47 =             ctot_ener_proj23*ctot_ener_proj37*(ctot_ener_proj38*ctot_ener_proj39 + ctot_ener_proj46);
const double ctot_ener_proj48 =             0.166666666666667*ctot_ener_proj47;
const double ctot_ener_proj49 =             DN_DX_0_1*U_0_2 + DN_DX_1_1*U_1_2 + DN_DX_2_1*U_2_2;
const double ctot_ener_proj50 =             ctot_ener_proj23*ctot_ener_proj49*(ctot_ener_proj39*ctot_ener_proj45 + ctot_ener_proj46);
const double ctot_ener_proj51 =             0.166666666666667*ctot_ener_proj50;
const double ctot_ener_proj52 =             DN_DX_0_0*U_0_0 + DN_DX_1_0*U_1_0 + DN_DX_2_0*U_2_0;
const double ctot_ener_proj53 =             0.5*gamma - 0.5;
const double ctot_ener_proj54 =             ctot_ener_proj23*ctot_ener_proj53*(ctot_ener_proj38 + ctot_ener_proj45) + ctot_ener_proj46;
const double ctot_ener_proj55 =             ctot_ener_proj10*ctot_ener_proj31*ctot_ener_proj52*ctot_ener_proj54;
const double ctot_ener_proj56 =             -0.166666666666667*ctot_ener_proj55;
const double ctot_ener_proj57 =             DN_DX_0_1*U_0_0 + DN_DX_1_1*U_1_0 + DN_DX_2_1*U_2_0;
const double ctot_ener_proj58 =             ctot_ener_proj17*ctot_ener_proj31*ctot_ener_proj54*ctot_ener_proj57;
const double ctot_ener_proj59 =             -0.166666666666667*ctot_ener_proj58;
const double ctot_ener_proj60 =             -0.25*dUdt_2_3;
const double ctot_ener_proj61 =             0.166666666666667*U_1_0;
const double ctot_ener_proj62 =             0.666666666666667*U_2_0 + ctot_ener_proj1 + ctot_ener_proj61;
const double ctot_ener_proj63 =             0.166666666666667*r[1];
const double ctot_ener_proj64 =             ctot_ener_proj62*(ctot_ener_proj4 + ctot_ener_proj63 + 0.666666666666667*r[2]);
const double ctot_ener_proj65 =             0.166666666666667*ctot_ener_proj64;
const double ctot_ener_proj66 =             0.666666666666667*U_0_0 + ctot_ener_proj2 + ctot_ener_proj61;
const double ctot_ener_proj67 =             ctot_ener_proj66*(ctot_ener_proj5 + ctot_ener_proj63 + 0.666666666666667*r[0]);
const double ctot_ener_proj68 =             0.166666666666667*U_1_1;
const double ctot_ener_proj69 =             0.666666666666667*U_2_1 + ctot_ener_proj68 + ctot_ener_proj8;
const double ctot_ener_proj70 =             0.166666666666667*f_ext(1,0);
const double ctot_ener_proj71 =             ctot_ener_proj69*(ctot_ener_proj11 + ctot_ener_proj70 + 0.666666666666667*f_ext(2,0));
const double ctot_ener_proj72 =             0.166666666666667*ctot_ener_proj71;
const double ctot_ener_proj73 =             0.666666666666667*U_0_1 + ctot_ener_proj68 + ctot_ener_proj9;
const double ctot_ener_proj74 =             ctot_ener_proj73*(ctot_ener_proj12 + ctot_ener_proj70 + 0.666666666666667*f_ext(0,0));
const double ctot_ener_proj75 =             0.166666666666667*U_1_2;
const double ctot_ener_proj76 =             0.666666666666667*U_2_2 + ctot_ener_proj15 + ctot_ener_proj75;
const double ctot_ener_proj77 =             0.166666666666667*f_ext(1,1);
const double ctot_ener_proj78 =             ctot_ener_proj76*(ctot_ener_proj18 + ctot_ener_proj77 + 0.666666666666667*f_ext(2,1));
const double ctot_ener_proj79 =             0.166666666666667*ctot_ener_proj78;
const double ctot_ener_proj80 =             0.666666666666667*U_0_2 + ctot_ener_proj16 + ctot_ener_proj75;
const double ctot_ener_proj81 =             ctot_ener_proj80*(ctot_ener_proj19 + ctot_ener_proj77 + 0.666666666666667*f_ext(0,1));
const double ctot_ener_proj82 =             1.0/ctot_ener_proj62;
const double ctot_ener_proj83 =             ctot_ener_proj22*ctot_ener_proj69*ctot_ener_proj82*gamma;
const double ctot_ener_proj84 =             -0.166666666666667*ctot_ener_proj83;
const double ctot_ener_proj85 =             ctot_ener_proj26*ctot_ener_proj76*ctot_ener_proj82*gamma;
const double ctot_ener_proj86 =             -0.166666666666667*ctot_ener_proj85;
const double ctot_ener_proj87 =             1.0/ctot_ener_proj66;
const double ctot_ener_proj88 =             ctot_ener_proj22*ctot_ener_proj73*ctot_ener_proj87*gamma;
const double ctot_ener_proj89 =             ctot_ener_proj26*ctot_ener_proj80*ctot_ener_proj87*gamma;
const double ctot_ener_proj90 =             pow(ctot_ener_proj62, -2);
const double ctot_ener_proj91 =             ctot_ener_proj29*ctot_ener_proj30*ctot_ener_proj69*ctot_ener_proj76*ctot_ener_proj90;
const double ctot_ener_proj92 =             0.166666666666667*ctot_ener_proj91;
const double ctot_ener_proj93 =             ctot_ener_proj29*ctot_ener_proj34*ctot_ener_proj69*ctot_ener_proj76*ctot_ener_proj90;
const double ctot_ener_proj94 =             0.166666666666667*ctot_ener_proj93;
const double ctot_ener_proj95 =             pow(ctot_ener_proj66, -2);
const double ctot_ener_proj96 =             ctot_ener_proj29*ctot_ener_proj30*ctot_ener_proj73*ctot_ener_proj80*ctot_ener_proj95;
const double ctot_ener_proj97 =             ctot_ener_proj29*ctot_ener_proj34*ctot_ener_proj73*ctot_ener_proj80*ctot_ener_proj95;
const double ctot_ener_proj98 =             pow(ctot_ener_proj69, 2);
const double ctot_ener_proj99 =             1.0*ctot_ener_proj29*ctot_ener_proj82;
const double ctot_ener_proj100 =             0.166666666666667*U_1_3;
const double ctot_ener_proj101 =             -ctot_ener_proj100;
const double ctot_ener_proj102 =             0.666666666666667*U_2_3;
const double ctot_ener_proj103 =             pow(ctot_ener_proj76, 2);
const double ctot_ener_proj104 =             ctot_ener_proj101 - ctot_ener_proj102 - ctot_ener_proj29*(ctot_ener_proj100 + ctot_ener_proj102 + ctot_ener_proj40 - ctot_ener_proj82*(0.5*ctot_ener_proj103 + 0.5*ctot_ener_proj98)) + ctot_ener_proj41;
const double ctot_ener_proj105 =             ctot_ener_proj37*ctot_ener_proj82*(ctot_ener_proj104 + ctot_ener_proj98*ctot_ener_proj99);
const double ctot_ener_proj106 =             0.166666666666667*ctot_ener_proj105;
const double ctot_ener_proj107 =             ctot_ener_proj49*ctot_ener_proj82*(ctot_ener_proj103*ctot_ener_proj99 + ctot_ener_proj104);
const double ctot_ener_proj108 =             0.166666666666667*ctot_ener_proj107;
const double ctot_ener_proj109 =             pow(ctot_ener_proj73, 2);
const double ctot_ener_proj110 =             1.0*ctot_ener_proj29*ctot_ener_proj87;
const double ctot_ener_proj111 =             0.666666666666667*U_0_3;
const double ctot_ener_proj112 =             pow(ctot_ener_proj80, 2);
const double ctot_ener_proj113 =             ctot_ener_proj101 - ctot_ener_proj111 - ctot_ener_proj29*(ctot_ener_proj100 + ctot_ener_proj111 + ctot_ener_proj43 - ctot_ener_proj87*(0.5*ctot_ener_proj109 + 0.5*ctot_ener_proj112)) + ctot_ener_proj44;
const double ctot_ener_proj114 =             ctot_ener_proj37*ctot_ener_proj87*(ctot_ener_proj109*ctot_ener_proj110 + ctot_ener_proj113);
const double ctot_ener_proj115 =             ctot_ener_proj49*ctot_ener_proj87*(ctot_ener_proj110*ctot_ener_proj112 + ctot_ener_proj113);
const double ctot_ener_proj116 =             ctot_ener_proj104 + ctot_ener_proj53*ctot_ener_proj82*(ctot_ener_proj103 + ctot_ener_proj98);
const double ctot_ener_proj117 =             ctot_ener_proj116*ctot_ener_proj52*ctot_ener_proj69*ctot_ener_proj90;
const double ctot_ener_proj118 =             -0.166666666666667*ctot_ener_proj117;
const double ctot_ener_proj119 =             ctot_ener_proj116*ctot_ener_proj57*ctot_ener_proj76*ctot_ener_proj90;
const double ctot_ener_proj120 =             -0.166666666666667*ctot_ener_proj119;
const double ctot_ener_proj121 =             ctot_ener_proj113 + ctot_ener_proj53*ctot_ener_proj87*(ctot_ener_proj109 + ctot_ener_proj112);
const double ctot_ener_proj122 =             ctot_ener_proj121*ctot_ener_proj52*ctot_ener_proj73*ctot_ener_proj95;
const double ctot_ener_proj123 =             ctot_ener_proj121*ctot_ener_proj57*ctot_ener_proj80*ctot_ener_proj95;
const double ctot_ener_proj124 =             0.166666666666667*ctot_ener_proj114 + 0.166666666666667*ctot_ener_proj115 - 0.166666666666667*ctot_ener_proj122 - 0.166666666666667*ctot_ener_proj123 + 0.166666666666667*ctot_ener_proj67 + 0.166666666666667*ctot_ener_proj74 + 0.166666666666667*ctot_ener_proj81 - 0.166666666666667*ctot_ener_proj88 - 0.166666666666667*ctot_ener_proj89 + 0.166666666666667*ctot_ener_proj96 + 0.166666666666667*ctot_ener_proj97 - 0.25*dUdt_0_3;
            tot_ener_proj[0]=ctot_ener_proj0 + ctot_ener_proj106 + ctot_ener_proj108 + 0.666666666666667*ctot_ener_proj114 + 0.666666666666667*ctot_ener_proj115 + ctot_ener_proj118 + ctot_ener_proj120 - 0.666666666666667*ctot_ener_proj122 - 0.666666666666667*ctot_ener_proj123 + ctot_ener_proj14 + ctot_ener_proj21 + ctot_ener_proj25 + ctot_ener_proj28 + ctot_ener_proj33 + ctot_ener_proj36 + ctot_ener_proj48 + ctot_ener_proj51 + ctot_ener_proj56 + ctot_ener_proj59 + ctot_ener_proj60 + ctot_ener_proj65 + 0.666666666666667*ctot_ener_proj67 + ctot_ener_proj7 + ctot_ener_proj72 + 0.666666666666667*ctot_ener_proj74 + ctot_ener_proj79 + 0.666666666666667*ctot_ener_proj81 + ctot_ener_proj84 + ctot_ener_proj86 - 0.666666666666667*ctot_ener_proj88 - 0.666666666666667*ctot_ener_proj89 + ctot_ener_proj92 + ctot_ener_proj94 + 0.666666666666667*ctot_ener_proj96 + 0.666666666666667*ctot_ener_proj97 - 0.5*dUdt_0_3;
            tot_ener_proj[1]=ctot_ener_proj106 + ctot_ener_proj108 + ctot_ener_proj118 + ctot_ener_proj120 + ctot_ener_proj124 + 0.666666666666667*ctot_ener_proj13 + 0.666666666666667*ctot_ener_proj20 - 0.666666666666667*ctot_ener_proj24 - 0.666666666666667*ctot_ener_proj27 + 0.666666666666667*ctot_ener_proj32 + 0.666666666666667*ctot_ener_proj35 + 0.666666666666667*ctot_ener_proj47 + 0.666666666666667*ctot_ener_proj50 - 0.666666666666667*ctot_ener_proj55 - 0.666666666666667*ctot_ener_proj58 + 0.666666666666667*ctot_ener_proj6 + ctot_ener_proj60 + ctot_ener_proj65 + ctot_ener_proj72 + ctot_ener_proj79 + ctot_ener_proj84 + ctot_ener_proj86 + ctot_ener_proj92 + ctot_ener_proj94 - 0.5*dUdt_1_3;
            tot_ener_proj[2]=ctot_ener_proj0 + 0.666666666666667*ctot_ener_proj105 + 0.666666666666667*ctot_ener_proj107 - 0.666666666666667*ctot_ener_proj117 - 0.666666666666667*ctot_ener_proj119 + ctot_ener_proj124 + ctot_ener_proj14 + ctot_ener_proj21 + ctot_ener_proj25 + ctot_ener_proj28 + ctot_ener_proj33 + ctot_ener_proj36 + ctot_ener_proj48 + ctot_ener_proj51 + ctot_ener_proj56 + ctot_ener_proj59 + 0.666666666666667*ctot_ener_proj64 + ctot_ener_proj7 + 0.666666666666667*ctot_ener_proj71 + 0.666666666666667*ctot_ener_proj78 - 0.666666666666667*ctot_ener_proj83 - 0.666666666666667*ctot_ener_proj85 + 0.666666666666667*ctot_ener_proj91 + 0.666666666666667*ctot_ener_proj93 - 0.5*dUdt_2_3;

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

        const double tau_m_avg = 1.0 / ((4.0 * stab_c1 * mu / 3.0 / rho_avg / std::pow(h, 2)) + (stab_c2 * (v_avg_norm + c_avg)/ h));
        const double tau_et_avg = 1.0 / ((stab_c1 * lambda / std::pow(h, 2)) + (stab_c2 * (v_avg_norm + c_avg)/ h));

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

        const double crRightHandSideBoundedVector0 =             DN_DX_0_0*U_0_1;
const double crRightHandSideBoundedVector1 =             DN_DX_0_1*U_0_2;
const double crRightHandSideBoundedVector2 =             DN_DX_1_0*U_1_1;
const double crRightHandSideBoundedVector3 =             DN_DX_1_1*U_1_2;
const double crRightHandSideBoundedVector4 =             DN_DX_2_0*U_2_1;
const double crRightHandSideBoundedVector5 =             DN_DX_2_1*U_2_2;
const double crRightHandSideBoundedVector6 =             crRightHandSideBoundedVector0 + crRightHandSideBoundedVector1 + crRightHandSideBoundedVector2 + crRightHandSideBoundedVector3 + crRightHandSideBoundedVector4 + crRightHandSideBoundedVector5;
const double crRightHandSideBoundedVector7 =             DN_DX_0_0*h;
const double crRightHandSideBoundedVector8 =             1.0/h;
const double crRightHandSideBoundedVector9 =             1.33333333333333*crRightHandSideBoundedVector8*mu*stab_c1;
const double crRightHandSideBoundedVector10 =             0.166666666666667*U_0_0;
const double crRightHandSideBoundedVector11 =             0.166666666666667*U_1_0;
const double crRightHandSideBoundedVector12 =             0.666666666666667*U_2_0 + crRightHandSideBoundedVector10 + crRightHandSideBoundedVector11;
const double crRightHandSideBoundedVector13 =             1.0/crRightHandSideBoundedVector12;
const double crRightHandSideBoundedVector14 =             pow(crRightHandSideBoundedVector12, -2);
const double crRightHandSideBoundedVector15 =             0.166666666666667*U_0_1;
const double crRightHandSideBoundedVector16 =             0.166666666666667*U_1_1;
const double crRightHandSideBoundedVector17 =             0.666666666666667*U_2_1;
const double crRightHandSideBoundedVector18 =             crRightHandSideBoundedVector15 + crRightHandSideBoundedVector16 + crRightHandSideBoundedVector17;
const double crRightHandSideBoundedVector19 =             pow(crRightHandSideBoundedVector18, 2);
const double crRightHandSideBoundedVector20 =             0.166666666666667*U_0_2;
const double crRightHandSideBoundedVector21 =             0.166666666666667*U_1_2;
const double crRightHandSideBoundedVector22 =             0.666666666666667*U_2_2;
const double crRightHandSideBoundedVector23 =             crRightHandSideBoundedVector20 + crRightHandSideBoundedVector21 + crRightHandSideBoundedVector22;
const double crRightHandSideBoundedVector24 =             pow(crRightHandSideBoundedVector23, 2);
const double crRightHandSideBoundedVector25 =             crRightHandSideBoundedVector19 + crRightHandSideBoundedVector24;
const double crRightHandSideBoundedVector26 =             sqrt(gamma);
const double crRightHandSideBoundedVector27 =             gamma - 1;
const double crRightHandSideBoundedVector28 =             crRightHandSideBoundedVector13*crRightHandSideBoundedVector27;
const double crRightHandSideBoundedVector29 =             0.166666666666667*U_0_3;
const double crRightHandSideBoundedVector30 =             -crRightHandSideBoundedVector29;
const double crRightHandSideBoundedVector31 =             0.166666666666667*U_1_3;
const double crRightHandSideBoundedVector32 =             -crRightHandSideBoundedVector31;
const double crRightHandSideBoundedVector33 =             0.666666666666667*U_2_3;
const double crRightHandSideBoundedVector34 =             -crRightHandSideBoundedVector33;
const double crRightHandSideBoundedVector35 =             0.5*crRightHandSideBoundedVector19;
const double crRightHandSideBoundedVector36 =             0.5*crRightHandSideBoundedVector24;
const double crRightHandSideBoundedVector37 =             crRightHandSideBoundedVector26*sqrt(-crRightHandSideBoundedVector28*(crRightHandSideBoundedVector13*crRightHandSideBoundedVector35 + crRightHandSideBoundedVector13*crRightHandSideBoundedVector36 + crRightHandSideBoundedVector30 + crRightHandSideBoundedVector32 + crRightHandSideBoundedVector34)) + sqrt(crRightHandSideBoundedVector14*crRightHandSideBoundedVector25);
const double crRightHandSideBoundedVector38 =             crRightHandSideBoundedVector37*stab_c2;
const double crRightHandSideBoundedVector39 =             1.0/(crRightHandSideBoundedVector13*crRightHandSideBoundedVector9 + crRightHandSideBoundedVector38);
const double crRightHandSideBoundedVector40 =             0.166666666666667*ResProj_1_1;
const double crRightHandSideBoundedVector41 =             0.166666666666667*dUdt_1_1;
const double crRightHandSideBoundedVector42 =             0.166666666666667*ResProj_0_1;
const double crRightHandSideBoundedVector43 =             0.166666666666667*dUdt_0_1;
const double crRightHandSideBoundedVector44 =             DN_DX_0_0*U_0_3 + DN_DX_1_0*U_1_3 + DN_DX_2_0*U_2_3;
const double crRightHandSideBoundedVector45 =             crRightHandSideBoundedVector27*crRightHandSideBoundedVector44;
const double crRightHandSideBoundedVector46 =             0.166666666666667*f_ext(0,0);
const double crRightHandSideBoundedVector47 =             0.166666666666667*f_ext(1,0);
const double crRightHandSideBoundedVector48 =             crRightHandSideBoundedVector46 + crRightHandSideBoundedVector47 + 0.666666666666667*f_ext(2,0);
const double crRightHandSideBoundedVector49 =             crRightHandSideBoundedVector12*crRightHandSideBoundedVector48;
const double crRightHandSideBoundedVector50 =             crRightHandSideBoundedVector1 + crRightHandSideBoundedVector3 + crRightHandSideBoundedVector5;
const double crRightHandSideBoundedVector51 =             crRightHandSideBoundedVector13*crRightHandSideBoundedVector18;
const double crRightHandSideBoundedVector52 =             crRightHandSideBoundedVector50*crRightHandSideBoundedVector51;
const double crRightHandSideBoundedVector53 =             DN_DX_0_1*U_0_1 + DN_DX_1_1*U_1_1 + DN_DX_2_1*U_2_1;
const double crRightHandSideBoundedVector54 =             crRightHandSideBoundedVector13*crRightHandSideBoundedVector23;
const double crRightHandSideBoundedVector55 =             crRightHandSideBoundedVector53*crRightHandSideBoundedVector54;
const double crRightHandSideBoundedVector56 =             DN_DX_0_0*U_0_2 + DN_DX_1_0*U_1_2 + DN_DX_2_0*U_2_2;
const double crRightHandSideBoundedVector57 =             1.0*crRightHandSideBoundedVector27*crRightHandSideBoundedVector56;
const double crRightHandSideBoundedVector58 =             1.0*gamma;
const double crRightHandSideBoundedVector59 =             -crRightHandSideBoundedVector58 + 3.0;
const double crRightHandSideBoundedVector60 =             crRightHandSideBoundedVector13*crRightHandSideBoundedVector59;
const double crRightHandSideBoundedVector61 =             crRightHandSideBoundedVector0 + crRightHandSideBoundedVector2 + crRightHandSideBoundedVector4;
const double crRightHandSideBoundedVector62 =             crRightHandSideBoundedVector18*crRightHandSideBoundedVector61;
const double crRightHandSideBoundedVector63 =             DN_DX_0_1*U_0_0 + DN_DX_1_1*U_1_0 + DN_DX_2_1*U_2_0;
const double crRightHandSideBoundedVector64 =             crRightHandSideBoundedVector14*crRightHandSideBoundedVector18*crRightHandSideBoundedVector23;
const double crRightHandSideBoundedVector65 =             crRightHandSideBoundedVector63*crRightHandSideBoundedVector64;
const double crRightHandSideBoundedVector66 =             DN_DX_0_0*U_0_0 + DN_DX_1_0*U_1_0 + DN_DX_2_0*U_2_0;
const double crRightHandSideBoundedVector67 =             0.5*gamma - 0.5;
const double crRightHandSideBoundedVector68 =             crRightHandSideBoundedVector25*crRightHandSideBoundedVector67;
const double crRightHandSideBoundedVector69 =             -crRightHandSideBoundedVector19 + crRightHandSideBoundedVector68;
const double crRightHandSideBoundedVector70 =             crRightHandSideBoundedVector14*crRightHandSideBoundedVector69;
const double crRightHandSideBoundedVector71 =             crRightHandSideBoundedVector66*crRightHandSideBoundedVector70;
const double crRightHandSideBoundedVector72 =             0.666666666666667*ResProj_2_1 + crRightHandSideBoundedVector40 + crRightHandSideBoundedVector41 + crRightHandSideBoundedVector42 + crRightHandSideBoundedVector43 + crRightHandSideBoundedVector45 - crRightHandSideBoundedVector49 + crRightHandSideBoundedVector52 - crRightHandSideBoundedVector54*crRightHandSideBoundedVector57 + crRightHandSideBoundedVector55 + crRightHandSideBoundedVector60*crRightHandSideBoundedVector62 - crRightHandSideBoundedVector65 + crRightHandSideBoundedVector71 + 0.666666666666667*dUdt_2_1;
const double crRightHandSideBoundedVector73 =             crRightHandSideBoundedVector39*crRightHandSideBoundedVector72;
const double crRightHandSideBoundedVector74 =             0.166666666666667*U_2_0;
const double crRightHandSideBoundedVector75 =             0.666666666666667*U_1_0 + crRightHandSideBoundedVector10 + crRightHandSideBoundedVector74;
const double crRightHandSideBoundedVector76 =             1.0/crRightHandSideBoundedVector75;
const double crRightHandSideBoundedVector77 =             pow(crRightHandSideBoundedVector75, -2);
const double crRightHandSideBoundedVector78 =             0.666666666666667*U_1_1;
const double crRightHandSideBoundedVector79 =             0.166666666666667*U_2_1;
const double crRightHandSideBoundedVector80 =             crRightHandSideBoundedVector15 + crRightHandSideBoundedVector78 + crRightHandSideBoundedVector79;
const double crRightHandSideBoundedVector81 =             pow(crRightHandSideBoundedVector80, 2);
const double crRightHandSideBoundedVector82 =             0.666666666666667*U_1_2;
const double crRightHandSideBoundedVector83 =             0.166666666666667*U_2_2;
const double crRightHandSideBoundedVector84 =             crRightHandSideBoundedVector20 + crRightHandSideBoundedVector82 + crRightHandSideBoundedVector83;
const double crRightHandSideBoundedVector85 =             pow(crRightHandSideBoundedVector84, 2);
const double crRightHandSideBoundedVector86 =             crRightHandSideBoundedVector81 + crRightHandSideBoundedVector85;
const double crRightHandSideBoundedVector87 =             crRightHandSideBoundedVector27*crRightHandSideBoundedVector76;
const double crRightHandSideBoundedVector88 =             0.666666666666667*U_1_3;
const double crRightHandSideBoundedVector89 =             -crRightHandSideBoundedVector88;
const double crRightHandSideBoundedVector90 =             0.166666666666667*U_2_3;
const double crRightHandSideBoundedVector91 =             -crRightHandSideBoundedVector90;
const double crRightHandSideBoundedVector92 =             0.5*crRightHandSideBoundedVector81;
const double crRightHandSideBoundedVector93 =             0.5*crRightHandSideBoundedVector85;
const double crRightHandSideBoundedVector94 =             crRightHandSideBoundedVector26*sqrt(-crRightHandSideBoundedVector87*(crRightHandSideBoundedVector30 + crRightHandSideBoundedVector76*crRightHandSideBoundedVector92 + crRightHandSideBoundedVector76*crRightHandSideBoundedVector93 + crRightHandSideBoundedVector89 + crRightHandSideBoundedVector91)) + sqrt(crRightHandSideBoundedVector77*crRightHandSideBoundedVector86);
const double crRightHandSideBoundedVector95 =             crRightHandSideBoundedVector94*stab_c2;
const double crRightHandSideBoundedVector96 =             1.0/(crRightHandSideBoundedVector76*crRightHandSideBoundedVector9 + crRightHandSideBoundedVector95);
const double crRightHandSideBoundedVector97 =             0.166666666666667*ResProj_2_1 + 0.166666666666667*dUdt_2_1;
const double crRightHandSideBoundedVector98 =             0.166666666666667*f_ext(2,0);
const double crRightHandSideBoundedVector99 =             crRightHandSideBoundedVector46 + crRightHandSideBoundedVector98 + 0.666666666666667*f_ext(1,0);
const double crRightHandSideBoundedVector100 =             crRightHandSideBoundedVector75*crRightHandSideBoundedVector99;
const double crRightHandSideBoundedVector101 =             crRightHandSideBoundedVector76*crRightHandSideBoundedVector80;
const double crRightHandSideBoundedVector102 =             crRightHandSideBoundedVector101*crRightHandSideBoundedVector50;
const double crRightHandSideBoundedVector103 =             crRightHandSideBoundedVector76*crRightHandSideBoundedVector84;
const double crRightHandSideBoundedVector104 =             crRightHandSideBoundedVector103*crRightHandSideBoundedVector53;
const double crRightHandSideBoundedVector105 =             crRightHandSideBoundedVector59*crRightHandSideBoundedVector76;
const double crRightHandSideBoundedVector106 =             crRightHandSideBoundedVector61*crRightHandSideBoundedVector80;
const double crRightHandSideBoundedVector107 =             crRightHandSideBoundedVector77*crRightHandSideBoundedVector80*crRightHandSideBoundedVector84;
const double crRightHandSideBoundedVector108 =             crRightHandSideBoundedVector107*crRightHandSideBoundedVector63;
const double crRightHandSideBoundedVector109 =             crRightHandSideBoundedVector67*crRightHandSideBoundedVector86;
const double crRightHandSideBoundedVector110 =             crRightHandSideBoundedVector109 - crRightHandSideBoundedVector81;
const double crRightHandSideBoundedVector111 =             crRightHandSideBoundedVector110*crRightHandSideBoundedVector77;
const double crRightHandSideBoundedVector112 =             crRightHandSideBoundedVector111*crRightHandSideBoundedVector66;
const double crRightHandSideBoundedVector113 =             0.666666666666667*ResProj_1_1 - crRightHandSideBoundedVector100 + crRightHandSideBoundedVector102 - crRightHandSideBoundedVector103*crRightHandSideBoundedVector57 + crRightHandSideBoundedVector104 + crRightHandSideBoundedVector105*crRightHandSideBoundedVector106 - crRightHandSideBoundedVector108 + crRightHandSideBoundedVector112 + crRightHandSideBoundedVector42 + crRightHandSideBoundedVector43 + crRightHandSideBoundedVector45 + crRightHandSideBoundedVector97 + 0.666666666666667*dUdt_1_1;
const double crRightHandSideBoundedVector114 =             crRightHandSideBoundedVector113*crRightHandSideBoundedVector96;
const double crRightHandSideBoundedVector115 =             0.666666666666667*U_0_0 + crRightHandSideBoundedVector11 + crRightHandSideBoundedVector74;
const double crRightHandSideBoundedVector116 =             1.0/crRightHandSideBoundedVector115;
const double crRightHandSideBoundedVector117 =             pow(crRightHandSideBoundedVector115, -2);
const double crRightHandSideBoundedVector118 =             0.666666666666667*U_0_1;
const double crRightHandSideBoundedVector119 =             crRightHandSideBoundedVector118 + crRightHandSideBoundedVector16 + crRightHandSideBoundedVector79;
const double crRightHandSideBoundedVector120 =             pow(crRightHandSideBoundedVector119, 2);
const double crRightHandSideBoundedVector121 =             0.666666666666667*U_0_2;
const double crRightHandSideBoundedVector122 =             crRightHandSideBoundedVector121 + crRightHandSideBoundedVector21 + crRightHandSideBoundedVector83;
const double crRightHandSideBoundedVector123 =             pow(crRightHandSideBoundedVector122, 2);
const double crRightHandSideBoundedVector124 =             crRightHandSideBoundedVector120 + crRightHandSideBoundedVector123;
const double crRightHandSideBoundedVector125 =             crRightHandSideBoundedVector116*crRightHandSideBoundedVector27;
const double crRightHandSideBoundedVector126 =             0.666666666666667*U_0_3;
const double crRightHandSideBoundedVector127 =             -crRightHandSideBoundedVector126;
const double crRightHandSideBoundedVector128 =             0.5*crRightHandSideBoundedVector120;
const double crRightHandSideBoundedVector129 =             0.5*crRightHandSideBoundedVector123;
const double crRightHandSideBoundedVector130 =             crRightHandSideBoundedVector26*sqrt(-crRightHandSideBoundedVector125*(crRightHandSideBoundedVector116*crRightHandSideBoundedVector128 + crRightHandSideBoundedVector116*crRightHandSideBoundedVector129 + crRightHandSideBoundedVector127 + crRightHandSideBoundedVector32 + crRightHandSideBoundedVector91)) + sqrt(crRightHandSideBoundedVector117*crRightHandSideBoundedVector124);
const double crRightHandSideBoundedVector131 =             crRightHandSideBoundedVector130*stab_c2;
const double crRightHandSideBoundedVector132 =             1.0/(crRightHandSideBoundedVector116*crRightHandSideBoundedVector9 + crRightHandSideBoundedVector131);
const double crRightHandSideBoundedVector133 =             crRightHandSideBoundedVector47 + crRightHandSideBoundedVector98 + 0.666666666666667*f_ext(0,0);
const double crRightHandSideBoundedVector134 =             crRightHandSideBoundedVector115*crRightHandSideBoundedVector133;
const double crRightHandSideBoundedVector135 =             crRightHandSideBoundedVector116*crRightHandSideBoundedVector119;
const double crRightHandSideBoundedVector136 =             crRightHandSideBoundedVector135*crRightHandSideBoundedVector50;
const double crRightHandSideBoundedVector137 =             crRightHandSideBoundedVector116*crRightHandSideBoundedVector122;
const double crRightHandSideBoundedVector138 =             crRightHandSideBoundedVector137*crRightHandSideBoundedVector53;
const double crRightHandSideBoundedVector139 =             crRightHandSideBoundedVector116*crRightHandSideBoundedVector59;
const double crRightHandSideBoundedVector140 =             crRightHandSideBoundedVector119*crRightHandSideBoundedVector61;
const double crRightHandSideBoundedVector141 =             crRightHandSideBoundedVector117*crRightHandSideBoundedVector119*crRightHandSideBoundedVector122;
const double crRightHandSideBoundedVector142 =             crRightHandSideBoundedVector141*crRightHandSideBoundedVector63;
const double crRightHandSideBoundedVector143 =             crRightHandSideBoundedVector124*crRightHandSideBoundedVector67;
const double crRightHandSideBoundedVector144 =             -crRightHandSideBoundedVector120 + crRightHandSideBoundedVector143;
const double crRightHandSideBoundedVector145 =             crRightHandSideBoundedVector117*crRightHandSideBoundedVector144;
const double crRightHandSideBoundedVector146 =             crRightHandSideBoundedVector145*crRightHandSideBoundedVector66;
const double crRightHandSideBoundedVector147 =             0.666666666666667*ResProj_0_1 - crRightHandSideBoundedVector134 + crRightHandSideBoundedVector136 - crRightHandSideBoundedVector137*crRightHandSideBoundedVector57 + crRightHandSideBoundedVector138 + crRightHandSideBoundedVector139*crRightHandSideBoundedVector140 - crRightHandSideBoundedVector142 + crRightHandSideBoundedVector146 + crRightHandSideBoundedVector40 + crRightHandSideBoundedVector41 + crRightHandSideBoundedVector45 + crRightHandSideBoundedVector97 + 0.666666666666667*dUdt_0_1;
const double crRightHandSideBoundedVector148 =             crRightHandSideBoundedVector132*crRightHandSideBoundedVector147;
const double crRightHandSideBoundedVector149 =             DN_DX_0_1*h;
const double crRightHandSideBoundedVector150 =             0.166666666666667*ResProj_1_2;
const double crRightHandSideBoundedVector151 =             0.166666666666667*dUdt_1_2;
const double crRightHandSideBoundedVector152 =             0.166666666666667*ResProj_0_2;
const double crRightHandSideBoundedVector153 =             0.166666666666667*dUdt_0_2;
const double crRightHandSideBoundedVector154 =             DN_DX_0_1*U_0_3 + DN_DX_1_1*U_1_3 + DN_DX_2_1*U_2_3;
const double crRightHandSideBoundedVector155 =             crRightHandSideBoundedVector154*crRightHandSideBoundedVector27;
const double crRightHandSideBoundedVector156 =             0.166666666666667*f_ext(0,1);
const double crRightHandSideBoundedVector157 =             0.166666666666667*f_ext(1,1);
const double crRightHandSideBoundedVector158 =             crRightHandSideBoundedVector156 + crRightHandSideBoundedVector157 + 0.666666666666667*f_ext(2,1);
const double crRightHandSideBoundedVector159 =             crRightHandSideBoundedVector12*crRightHandSideBoundedVector158;
const double crRightHandSideBoundedVector160 =             crRightHandSideBoundedVector51*crRightHandSideBoundedVector56;
const double crRightHandSideBoundedVector161 =             crRightHandSideBoundedVector54*crRightHandSideBoundedVector61;
const double crRightHandSideBoundedVector162 =             1.0*crRightHandSideBoundedVector27*crRightHandSideBoundedVector53;
const double crRightHandSideBoundedVector163 =             crRightHandSideBoundedVector23*crRightHandSideBoundedVector50;
const double crRightHandSideBoundedVector164 =             crRightHandSideBoundedVector64*crRightHandSideBoundedVector66;
const double crRightHandSideBoundedVector165 =             -crRightHandSideBoundedVector24 + crRightHandSideBoundedVector68;
const double crRightHandSideBoundedVector166 =             crRightHandSideBoundedVector14*crRightHandSideBoundedVector165;
const double crRightHandSideBoundedVector167 =             crRightHandSideBoundedVector166*crRightHandSideBoundedVector63;
const double crRightHandSideBoundedVector168 =             0.666666666666667*ResProj_2_2 + crRightHandSideBoundedVector150 + crRightHandSideBoundedVector151 + crRightHandSideBoundedVector152 + crRightHandSideBoundedVector153 + crRightHandSideBoundedVector155 - crRightHandSideBoundedVector159 + crRightHandSideBoundedVector160 + crRightHandSideBoundedVector161 - crRightHandSideBoundedVector162*crRightHandSideBoundedVector51 + crRightHandSideBoundedVector163*crRightHandSideBoundedVector60 - crRightHandSideBoundedVector164 + crRightHandSideBoundedVector167 + 0.666666666666667*dUdt_2_2;
const double crRightHandSideBoundedVector169 =             crRightHandSideBoundedVector168*crRightHandSideBoundedVector39;
const double crRightHandSideBoundedVector170 =             0.166666666666667*ResProj_2_2 + 0.166666666666667*dUdt_2_2;
const double crRightHandSideBoundedVector171 =             0.166666666666667*f_ext(2,1);
const double crRightHandSideBoundedVector172 =             crRightHandSideBoundedVector156 + crRightHandSideBoundedVector171 + 0.666666666666667*f_ext(1,1);
const double crRightHandSideBoundedVector173 =             crRightHandSideBoundedVector172*crRightHandSideBoundedVector75;
const double crRightHandSideBoundedVector174 =             crRightHandSideBoundedVector101*crRightHandSideBoundedVector56;
const double crRightHandSideBoundedVector175 =             crRightHandSideBoundedVector103*crRightHandSideBoundedVector61;
const double crRightHandSideBoundedVector176 =             crRightHandSideBoundedVector50*crRightHandSideBoundedVector84;
const double crRightHandSideBoundedVector177 =             crRightHandSideBoundedVector107*crRightHandSideBoundedVector66;
const double crRightHandSideBoundedVector178 =             crRightHandSideBoundedVector109 - crRightHandSideBoundedVector85;
const double crRightHandSideBoundedVector179 =             crRightHandSideBoundedVector178*crRightHandSideBoundedVector77;
const double crRightHandSideBoundedVector180 =             crRightHandSideBoundedVector179*crRightHandSideBoundedVector63;
const double crRightHandSideBoundedVector181 =             0.666666666666667*ResProj_1_2 - crRightHandSideBoundedVector101*crRightHandSideBoundedVector162 + crRightHandSideBoundedVector105*crRightHandSideBoundedVector176 + crRightHandSideBoundedVector152 + crRightHandSideBoundedVector153 + crRightHandSideBoundedVector155 + crRightHandSideBoundedVector170 - crRightHandSideBoundedVector173 + crRightHandSideBoundedVector174 + crRightHandSideBoundedVector175 - crRightHandSideBoundedVector177 + crRightHandSideBoundedVector180 + 0.666666666666667*dUdt_1_2;
const double crRightHandSideBoundedVector182 =             crRightHandSideBoundedVector181*crRightHandSideBoundedVector96;
const double crRightHandSideBoundedVector183 =             crRightHandSideBoundedVector157 + crRightHandSideBoundedVector171 + 0.666666666666667*f_ext(0,1);
const double crRightHandSideBoundedVector184 =             crRightHandSideBoundedVector115*crRightHandSideBoundedVector183;
const double crRightHandSideBoundedVector185 =             crRightHandSideBoundedVector135*crRightHandSideBoundedVector56;
const double crRightHandSideBoundedVector186 =             crRightHandSideBoundedVector137*crRightHandSideBoundedVector61;
const double crRightHandSideBoundedVector187 =             crRightHandSideBoundedVector122*crRightHandSideBoundedVector50;
const double crRightHandSideBoundedVector188 =             crRightHandSideBoundedVector141*crRightHandSideBoundedVector66;
const double crRightHandSideBoundedVector189 =             -crRightHandSideBoundedVector123 + crRightHandSideBoundedVector143;
const double crRightHandSideBoundedVector190 =             crRightHandSideBoundedVector117*crRightHandSideBoundedVector189;
const double crRightHandSideBoundedVector191 =             crRightHandSideBoundedVector190*crRightHandSideBoundedVector63;
const double crRightHandSideBoundedVector192 =             0.666666666666667*ResProj_0_2 - crRightHandSideBoundedVector135*crRightHandSideBoundedVector162 + crRightHandSideBoundedVector139*crRightHandSideBoundedVector187 + crRightHandSideBoundedVector150 + crRightHandSideBoundedVector151 + crRightHandSideBoundedVector155 + crRightHandSideBoundedVector170 - crRightHandSideBoundedVector184 + crRightHandSideBoundedVector185 + crRightHandSideBoundedVector186 - crRightHandSideBoundedVector188 + crRightHandSideBoundedVector191 + 0.666666666666667*dUdt_0_2;
const double crRightHandSideBoundedVector193 =             crRightHandSideBoundedVector132*crRightHandSideBoundedVector192;
const double crRightHandSideBoundedVector194 =             0.166666666666667*crRightHandSideBoundedVector100;
const double crRightHandSideBoundedVector195 =             -0.166666666666667*crRightHandSideBoundedVector102;
const double crRightHandSideBoundedVector196 =             -0.166666666666667*crRightHandSideBoundedVector104;
const double crRightHandSideBoundedVector197 =             crRightHandSideBoundedVector27*crRightHandSideBoundedVector56;
const double crRightHandSideBoundedVector198 =             0.166666666666667*crRightHandSideBoundedVector197;
const double crRightHandSideBoundedVector199 =             crRightHandSideBoundedVector103*crRightHandSideBoundedVector198;
const double crRightHandSideBoundedVector200 =             crRightHandSideBoundedVector58 - 3.0;
const double crRightHandSideBoundedVector201 =             crRightHandSideBoundedVector200*crRightHandSideBoundedVector61;
const double crRightHandSideBoundedVector202 =             0.166666666666667*crRightHandSideBoundedVector201;
const double crRightHandSideBoundedVector203 =             crRightHandSideBoundedVector101*crRightHandSideBoundedVector202;
const double crRightHandSideBoundedVector204 =             0.166666666666667*crRightHandSideBoundedVector108;
const double crRightHandSideBoundedVector205 =             -0.166666666666667*crRightHandSideBoundedVector112;
const double crRightHandSideBoundedVector206 =             0.166666666666667*crRightHandSideBoundedVector49;
const double crRightHandSideBoundedVector207 =             -0.166666666666667*crRightHandSideBoundedVector52;
const double crRightHandSideBoundedVector208 =             -0.166666666666667*crRightHandSideBoundedVector55;
const double crRightHandSideBoundedVector209 =             crRightHandSideBoundedVector198*crRightHandSideBoundedVector54;
const double crRightHandSideBoundedVector210 =             0.666666666666667*crRightHandSideBoundedVector197;
const double crRightHandSideBoundedVector211 =             crRightHandSideBoundedVector202*crRightHandSideBoundedVector51;
const double crRightHandSideBoundedVector212 =             0.666666666666667*crRightHandSideBoundedVector201;
const double crRightHandSideBoundedVector213 =             0.166666666666667*crRightHandSideBoundedVector65;
const double crRightHandSideBoundedVector214 =             crRightHandSideBoundedVector117*crRightHandSideBoundedVector119;
const double crRightHandSideBoundedVector215 =             0.666666666666667*DN_DX_0_1*U_0_0 + 0.666666666666667*DN_DX_1_1*U_1_0 + 0.666666666666667*DN_DX_2_1*U_2_0;
const double crRightHandSideBoundedVector216 =             crRightHandSideBoundedVector122*crRightHandSideBoundedVector215;
const double crRightHandSideBoundedVector217 =             -0.166666666666667*crRightHandSideBoundedVector71;
const double crRightHandSideBoundedVector218 =             0.666666666666667*DN_DX_0_0*U_0_0 + 0.666666666666667*DN_DX_1_0*U_1_0 + 0.666666666666667*DN_DX_2_0*U_2_0;
const double crRightHandSideBoundedVector219 =             crRightHandSideBoundedVector53 + crRightHandSideBoundedVector56;
const double crRightHandSideBoundedVector220 =             DN_DX_0_1*crRightHandSideBoundedVector219*mu;
const double crRightHandSideBoundedVector221 =             1.0/lin_m_norm;
const double crRightHandSideBoundedVector222 =             1.0/mu;
const double crRightHandSideBoundedVector223 =             crRightHandSideBoundedVector221*crRightHandSideBoundedVector222*lin_m[0]*lin_m[1]*nu_st;
const double crRightHandSideBoundedVector224 =             crRightHandSideBoundedVector12*crRightHandSideBoundedVector223 + 1;
const double crRightHandSideBoundedVector225 =             crRightHandSideBoundedVector221*lin_m[0]*lin_m[1];
const double crRightHandSideBoundedVector226 =             crRightHandSideBoundedVector222*nu_sc*(crRightHandSideBoundedVector225 - 1);
const double crRightHandSideBoundedVector227 =             -crRightHandSideBoundedVector12*crRightHandSideBoundedVector226 + crRightHandSideBoundedVector224;
const double crRightHandSideBoundedVector228 =             crRightHandSideBoundedVector223*crRightHandSideBoundedVector75 + 1;
const double crRightHandSideBoundedVector229 =             -crRightHandSideBoundedVector226*crRightHandSideBoundedVector75 + crRightHandSideBoundedVector228;
const double crRightHandSideBoundedVector230 =             crRightHandSideBoundedVector115*crRightHandSideBoundedVector223 + 1;
const double crRightHandSideBoundedVector231 =             -crRightHandSideBoundedVector115*crRightHandSideBoundedVector226 + crRightHandSideBoundedVector230;
const double crRightHandSideBoundedVector232 =             -nu_sc + nu_st;
const double crRightHandSideBoundedVector233 =             crRightHandSideBoundedVector12*crRightHandSideBoundedVector221*crRightHandSideBoundedVector232*lin_m[0]*lin_m[1];
const double crRightHandSideBoundedVector234 =             2*crRightHandSideBoundedVector1 + 2*crRightHandSideBoundedVector3 + 2*crRightHandSideBoundedVector5;
const double crRightHandSideBoundedVector235 =             0.111111111111111*U_0_0;
const double crRightHandSideBoundedVector236 =             0.111111111111111*U_1_0;
const double crRightHandSideBoundedVector237 =             0.444444444444444*U_2_0 + crRightHandSideBoundedVector235 + crRightHandSideBoundedVector236;
const double crRightHandSideBoundedVector238 =             -crRightHandSideBoundedVector14*(-crRightHandSideBoundedVector18*crRightHandSideBoundedVector218 - crRightHandSideBoundedVector215*crRightHandSideBoundedVector23 + crRightHandSideBoundedVector237*crRightHandSideBoundedVector50 + crRightHandSideBoundedVector237*crRightHandSideBoundedVector61);
const double crRightHandSideBoundedVector239 =             crRightHandSideBoundedVector234 + crRightHandSideBoundedVector238;
const double crRightHandSideBoundedVector240 =             pow(lin_m[0], 2);
const double crRightHandSideBoundedVector241 =             crRightHandSideBoundedVector221*crRightHandSideBoundedVector222*crRightHandSideBoundedVector240*nu_st;
const double crRightHandSideBoundedVector242 =             -crRightHandSideBoundedVector221*crRightHandSideBoundedVector240 + 1;
const double crRightHandSideBoundedVector243 =             crRightHandSideBoundedVector222*crRightHandSideBoundedVector242*nu_sc;
const double crRightHandSideBoundedVector244 =             2*crRightHandSideBoundedVector0 + 2*crRightHandSideBoundedVector2 + 2*crRightHandSideBoundedVector4;
const double crRightHandSideBoundedVector245 =             crRightHandSideBoundedVector238 + crRightHandSideBoundedVector244;
const double crRightHandSideBoundedVector246 =             crRightHandSideBoundedVector233*crRightHandSideBoundedVector239 + crRightHandSideBoundedVector245*mu*(crRightHandSideBoundedVector12*crRightHandSideBoundedVector241 + crRightHandSideBoundedVector12*crRightHandSideBoundedVector243 + 1);
const double crRightHandSideBoundedVector247 =             crRightHandSideBoundedVector221*crRightHandSideBoundedVector232*crRightHandSideBoundedVector75*lin_m[0]*lin_m[1];
const double crRightHandSideBoundedVector248 =             0.111111111111111*U_2_0;
const double crRightHandSideBoundedVector249 =             0.444444444444444*U_1_0 + crRightHandSideBoundedVector235 + crRightHandSideBoundedVector248;
const double crRightHandSideBoundedVector250 =             -crRightHandSideBoundedVector77*(-crRightHandSideBoundedVector215*crRightHandSideBoundedVector84 - crRightHandSideBoundedVector218*crRightHandSideBoundedVector80 + crRightHandSideBoundedVector249*crRightHandSideBoundedVector50 + crRightHandSideBoundedVector249*crRightHandSideBoundedVector61);
const double crRightHandSideBoundedVector251 =             crRightHandSideBoundedVector234 + crRightHandSideBoundedVector250;
const double crRightHandSideBoundedVector252 =             crRightHandSideBoundedVector244 + crRightHandSideBoundedVector250;
const double crRightHandSideBoundedVector253 =             crRightHandSideBoundedVector247*crRightHandSideBoundedVector251 + crRightHandSideBoundedVector252*mu*(crRightHandSideBoundedVector241*crRightHandSideBoundedVector75 + crRightHandSideBoundedVector243*crRightHandSideBoundedVector75 + 1);
const double crRightHandSideBoundedVector254 =             crRightHandSideBoundedVector115*crRightHandSideBoundedVector221*crRightHandSideBoundedVector232*lin_m[0]*lin_m[1];
const double crRightHandSideBoundedVector255 =             0.444444444444444*U_0_0 + crRightHandSideBoundedVector236 + crRightHandSideBoundedVector248;
const double crRightHandSideBoundedVector256 =             -crRightHandSideBoundedVector117*(-crRightHandSideBoundedVector119*crRightHandSideBoundedVector218 - crRightHandSideBoundedVector122*crRightHandSideBoundedVector215 + crRightHandSideBoundedVector255*crRightHandSideBoundedVector50 + crRightHandSideBoundedVector255*crRightHandSideBoundedVector61);
const double crRightHandSideBoundedVector257 =             crRightHandSideBoundedVector234 + crRightHandSideBoundedVector256;
const double crRightHandSideBoundedVector258 =             crRightHandSideBoundedVector244 + crRightHandSideBoundedVector256;
const double crRightHandSideBoundedVector259 =             crRightHandSideBoundedVector254*crRightHandSideBoundedVector257 + crRightHandSideBoundedVector258*mu*(crRightHandSideBoundedVector115*crRightHandSideBoundedVector241 + crRightHandSideBoundedVector115*crRightHandSideBoundedVector243 + 1);
const double crRightHandSideBoundedVector260 =             1.0/stab_c2;
const double crRightHandSideBoundedVector261 =             0.166666666666667*ResProj_1_0;
const double crRightHandSideBoundedVector262 =             0.166666666666667*dUdt_1_0;
const double crRightHandSideBoundedVector263 =             0.166666666666667*ResProj_0_0;
const double crRightHandSideBoundedVector264 =             0.166666666666667*dUdt_0_0;
const double crRightHandSideBoundedVector265 =             1.0*crRightHandSideBoundedVector260*h*(0.666666666666667*ResProj_2_0 + crRightHandSideBoundedVector0 + crRightHandSideBoundedVector1 + crRightHandSideBoundedVector2 + crRightHandSideBoundedVector261 + crRightHandSideBoundedVector262 + crRightHandSideBoundedVector263 + crRightHandSideBoundedVector264 + crRightHandSideBoundedVector3 + crRightHandSideBoundedVector4 + crRightHandSideBoundedVector5 + 0.666666666666667*dUdt_2_0)/crRightHandSideBoundedVector37;
const double crRightHandSideBoundedVector266 =             DN_DX_0_1*crRightHandSideBoundedVector18;
const double crRightHandSideBoundedVector267 =             crRightHandSideBoundedVector14*crRightHandSideBoundedVector23;
const double crRightHandSideBoundedVector268 =             0.0277777777777778*f_ext(0,0);
const double crRightHandSideBoundedVector269 =             0.0277777777777778*f_ext(1,0);
const double crRightHandSideBoundedVector270 =             0.111111111111111*f_ext(2,0);
const double crRightHandSideBoundedVector271 =             crRightHandSideBoundedVector14*crRightHandSideBoundedVector18*crRightHandSideBoundedVector50;
const double crRightHandSideBoundedVector272 =             0.166666666666667*crRightHandSideBoundedVector271;
const double crRightHandSideBoundedVector273 =             crRightHandSideBoundedVector14*crRightHandSideBoundedVector23*crRightHandSideBoundedVector53;
const double crRightHandSideBoundedVector274 =             crRightHandSideBoundedVector198*crRightHandSideBoundedVector267;
const double crRightHandSideBoundedVector275 =             crRightHandSideBoundedVector14*crRightHandSideBoundedVector18;
const double crRightHandSideBoundedVector276 =             pow(crRightHandSideBoundedVector12, -3);
const double crRightHandSideBoundedVector277 =             crRightHandSideBoundedVector18*crRightHandSideBoundedVector23*crRightHandSideBoundedVector276*crRightHandSideBoundedVector63;
const double crRightHandSideBoundedVector278 =             0.333333333333333*crRightHandSideBoundedVector277;
const double crRightHandSideBoundedVector279 =             crRightHandSideBoundedVector276*crRightHandSideBoundedVector66*crRightHandSideBoundedVector69;
const double crRightHandSideBoundedVector280 =             crRightHandSideBoundedVector202*crRightHandSideBoundedVector275 + crRightHandSideBoundedVector268 + crRightHandSideBoundedVector269 + crRightHandSideBoundedVector270 - crRightHandSideBoundedVector272 - 0.166666666666667*crRightHandSideBoundedVector273 + crRightHandSideBoundedVector274 + crRightHandSideBoundedVector278 - 0.333333333333333*crRightHandSideBoundedVector279;
const double crRightHandSideBoundedVector281 =             0.166666666666667*ResProj_2_0;
const double crRightHandSideBoundedVector282 =             0.166666666666667*dUdt_2_0;
const double crRightHandSideBoundedVector283 =             1.0*crRightHandSideBoundedVector260*h*(0.666666666666667*ResProj_1_0 + crRightHandSideBoundedVector0 + crRightHandSideBoundedVector1 + crRightHandSideBoundedVector2 + crRightHandSideBoundedVector263 + crRightHandSideBoundedVector264 + crRightHandSideBoundedVector281 + crRightHandSideBoundedVector282 + crRightHandSideBoundedVector3 + crRightHandSideBoundedVector4 + crRightHandSideBoundedVector5 + 0.666666666666667*dUdt_1_0)/crRightHandSideBoundedVector94;
const double crRightHandSideBoundedVector284 =             DN_DX_0_1*crRightHandSideBoundedVector80;
const double crRightHandSideBoundedVector285 =             crRightHandSideBoundedVector77*crRightHandSideBoundedVector84;
const double crRightHandSideBoundedVector286 =             0.111111111111111*f_ext(1,0);
const double crRightHandSideBoundedVector287 =             0.0277777777777778*f_ext(2,0);
const double crRightHandSideBoundedVector288 =             crRightHandSideBoundedVector50*crRightHandSideBoundedVector77*crRightHandSideBoundedVector80;
const double crRightHandSideBoundedVector289 =             0.166666666666667*crRightHandSideBoundedVector288;
const double crRightHandSideBoundedVector290 =             crRightHandSideBoundedVector53*crRightHandSideBoundedVector77*crRightHandSideBoundedVector84;
const double crRightHandSideBoundedVector291 =             crRightHandSideBoundedVector198*crRightHandSideBoundedVector285;
const double crRightHandSideBoundedVector292 =             crRightHandSideBoundedVector77*crRightHandSideBoundedVector80;
const double crRightHandSideBoundedVector293 =             pow(crRightHandSideBoundedVector75, -3);
const double crRightHandSideBoundedVector294 =             crRightHandSideBoundedVector293*crRightHandSideBoundedVector63*crRightHandSideBoundedVector80*crRightHandSideBoundedVector84;
const double crRightHandSideBoundedVector295 =             0.333333333333333*crRightHandSideBoundedVector294;
const double crRightHandSideBoundedVector296 =             crRightHandSideBoundedVector110*crRightHandSideBoundedVector293*crRightHandSideBoundedVector66;
const double crRightHandSideBoundedVector297 =             crRightHandSideBoundedVector202*crRightHandSideBoundedVector292 + crRightHandSideBoundedVector268 + crRightHandSideBoundedVector286 + crRightHandSideBoundedVector287 - crRightHandSideBoundedVector289 - 0.166666666666667*crRightHandSideBoundedVector290 + crRightHandSideBoundedVector291 + crRightHandSideBoundedVector295 - 0.333333333333333*crRightHandSideBoundedVector296;
const double crRightHandSideBoundedVector298 =             1.0*crRightHandSideBoundedVector260*h*(0.666666666666667*ResProj_0_0 + crRightHandSideBoundedVector261 + crRightHandSideBoundedVector262 + crRightHandSideBoundedVector281 + crRightHandSideBoundedVector282 + crRightHandSideBoundedVector6 + 0.666666666666667*dUdt_0_0)/crRightHandSideBoundedVector130;
const double crRightHandSideBoundedVector299 =             crRightHandSideBoundedVector270 + crRightHandSideBoundedVector286 + 0.444444444444444*f_ext(0,0);
const double crRightHandSideBoundedVector300 =             crRightHandSideBoundedVector117*crRightHandSideBoundedVector119*crRightHandSideBoundedVector50;
const double crRightHandSideBoundedVector301 =             0.666666666666667*crRightHandSideBoundedVector300;
const double crRightHandSideBoundedVector302 =             crRightHandSideBoundedVector117*crRightHandSideBoundedVector122*crRightHandSideBoundedVector53;
const double crRightHandSideBoundedVector303 =             DN_DX_0_1*crRightHandSideBoundedVector119;
const double crRightHandSideBoundedVector304 =             crRightHandSideBoundedVector117*crRightHandSideBoundedVector122;
const double crRightHandSideBoundedVector305 =             crRightHandSideBoundedVector210*crRightHandSideBoundedVector304;
const double crRightHandSideBoundedVector306 =             pow(crRightHandSideBoundedVector115, -3);
const double crRightHandSideBoundedVector307 =             crRightHandSideBoundedVector119*crRightHandSideBoundedVector122*crRightHandSideBoundedVector306*crRightHandSideBoundedVector63;
const double crRightHandSideBoundedVector308 =             1.33333333333333*crRightHandSideBoundedVector307;
const double crRightHandSideBoundedVector309 =             crRightHandSideBoundedVector144*crRightHandSideBoundedVector306*crRightHandSideBoundedVector66;
const double crRightHandSideBoundedVector310 =             crRightHandSideBoundedVector58 - 1.0;
const double crRightHandSideBoundedVector311 =             DN_DX_0_0*crRightHandSideBoundedVector23;
const double crRightHandSideBoundedVector312 =             -DN_DX_0_1*crRightHandSideBoundedVector15;
const double crRightHandSideBoundedVector313 =             -DN_DX_1_1*crRightHandSideBoundedVector16;
const double crRightHandSideBoundedVector314 =             -DN_DX_2_1*crRightHandSideBoundedVector79;
const double crRightHandSideBoundedVector315 =             0.166666666666667*crRightHandSideBoundedVector13*crRightHandSideBoundedVector18;
const double crRightHandSideBoundedVector316 =             crRightHandSideBoundedVector315*crRightHandSideBoundedVector63;
const double crRightHandSideBoundedVector317 =             0.166666666666667*crRightHandSideBoundedVector13*crRightHandSideBoundedVector23;
const double crRightHandSideBoundedVector318 =             crRightHandSideBoundedVector317*crRightHandSideBoundedVector66;
const double crRightHandSideBoundedVector319 =             crRightHandSideBoundedVector198 - crRightHandSideBoundedVector27*crRightHandSideBoundedVector318 + crRightHandSideBoundedVector312 + crRightHandSideBoundedVector313 + crRightHandSideBoundedVector314 + crRightHandSideBoundedVector316;
const double crRightHandSideBoundedVector320 =             1.0*crRightHandSideBoundedVector13*crRightHandSideBoundedVector168*crRightHandSideBoundedVector39*h;
const double crRightHandSideBoundedVector321 =             DN_DX_0_0*crRightHandSideBoundedVector84;
const double crRightHandSideBoundedVector322 =             0.166666666666667*crRightHandSideBoundedVector76*crRightHandSideBoundedVector80;
const double crRightHandSideBoundedVector323 =             crRightHandSideBoundedVector322*crRightHandSideBoundedVector63;
const double crRightHandSideBoundedVector324 =             0.166666666666667*crRightHandSideBoundedVector76*crRightHandSideBoundedVector84;
const double crRightHandSideBoundedVector325 =             crRightHandSideBoundedVector324*crRightHandSideBoundedVector66;
const double crRightHandSideBoundedVector326 =             crRightHandSideBoundedVector198 - crRightHandSideBoundedVector27*crRightHandSideBoundedVector325 + crRightHandSideBoundedVector312 + crRightHandSideBoundedVector313 + crRightHandSideBoundedVector314 + crRightHandSideBoundedVector323;
const double crRightHandSideBoundedVector327 =             1.0*crRightHandSideBoundedVector181*crRightHandSideBoundedVector76*crRightHandSideBoundedVector96*h;
const double crRightHandSideBoundedVector328 =             -DN_DX_0_1*crRightHandSideBoundedVector118 - DN_DX_1_1*crRightHandSideBoundedVector78 - DN_DX_2_1*crRightHandSideBoundedVector17 + crRightHandSideBoundedVector210;
const double crRightHandSideBoundedVector329 =             DN_DX_0_0*crRightHandSideBoundedVector122;
const double crRightHandSideBoundedVector330 =             1.0*crRightHandSideBoundedVector116*crRightHandSideBoundedVector132*crRightHandSideBoundedVector192*h;
const double crRightHandSideBoundedVector331 =             DN_DX_0_1*crRightHandSideBoundedVector23;
const double crRightHandSideBoundedVector332 =             DN_DX_0_0*crRightHandSideBoundedVector18;
const double crRightHandSideBoundedVector333 =             DN_DX_0_1*crRightHandSideBoundedVector20;
const double crRightHandSideBoundedVector334 =             DN_DX_1_1*crRightHandSideBoundedVector21;
const double crRightHandSideBoundedVector335 =             DN_DX_2_1*crRightHandSideBoundedVector83;
const double crRightHandSideBoundedVector336 =             -crRightHandSideBoundedVector202;
const double crRightHandSideBoundedVector337 =             crRightHandSideBoundedVector317*crRightHandSideBoundedVector63;
const double crRightHandSideBoundedVector338 =             -crRightHandSideBoundedVector337;
const double crRightHandSideBoundedVector339 =             crRightHandSideBoundedVector315*crRightHandSideBoundedVector66;
const double crRightHandSideBoundedVector340 =             crRightHandSideBoundedVector200*crRightHandSideBoundedVector339 + crRightHandSideBoundedVector333 + crRightHandSideBoundedVector334 + crRightHandSideBoundedVector335 + crRightHandSideBoundedVector336 + crRightHandSideBoundedVector338;
const double crRightHandSideBoundedVector341 =             1.0*crRightHandSideBoundedVector13*crRightHandSideBoundedVector39*crRightHandSideBoundedVector72*h;
const double crRightHandSideBoundedVector342 =             DN_DX_0_1*crRightHandSideBoundedVector84;
const double crRightHandSideBoundedVector343 =             DN_DX_0_0*crRightHandSideBoundedVector80;
const double crRightHandSideBoundedVector344 =             crRightHandSideBoundedVector324*crRightHandSideBoundedVector63;
const double crRightHandSideBoundedVector345 =             -crRightHandSideBoundedVector344;
const double crRightHandSideBoundedVector346 =             crRightHandSideBoundedVector322*crRightHandSideBoundedVector66;
const double crRightHandSideBoundedVector347 =             crRightHandSideBoundedVector200*crRightHandSideBoundedVector346 + crRightHandSideBoundedVector333 + crRightHandSideBoundedVector334 + crRightHandSideBoundedVector335 + crRightHandSideBoundedVector336 + crRightHandSideBoundedVector345;
const double crRightHandSideBoundedVector348 =             1.0*crRightHandSideBoundedVector113*crRightHandSideBoundedVector76*crRightHandSideBoundedVector96*h;
const double crRightHandSideBoundedVector349 =             DN_DX_0_1*crRightHandSideBoundedVector121;
const double crRightHandSideBoundedVector350 =             DN_DX_1_1*crRightHandSideBoundedVector82;
const double crRightHandSideBoundedVector351 =             DN_DX_2_1*crRightHandSideBoundedVector22;
const double crRightHandSideBoundedVector352 =             -crRightHandSideBoundedVector212 + crRightHandSideBoundedVector349 + crRightHandSideBoundedVector350 + crRightHandSideBoundedVector351;
const double crRightHandSideBoundedVector353 =             DN_DX_0_1*crRightHandSideBoundedVector122;
const double crRightHandSideBoundedVector354 =             DN_DX_0_0*crRightHandSideBoundedVector119;
const double crRightHandSideBoundedVector355 =             -crRightHandSideBoundedVector116*crRightHandSideBoundedVector216;
const double crRightHandSideBoundedVector356 =             1.0*crRightHandSideBoundedVector116*crRightHandSideBoundedVector132*crRightHandSideBoundedVector147*h;
const double crRightHandSideBoundedVector357 =             DN_DX_0_0*crRightHandSideBoundedVector310*h;
const double crRightHandSideBoundedVector358 =             1.0/c_v;
const double crRightHandSideBoundedVector359 =             crRightHandSideBoundedVector358*crRightHandSideBoundedVector8*lambda*stab_c1/gamma;
const double crRightHandSideBoundedVector360 =             1.0/(crRightHandSideBoundedVector13*crRightHandSideBoundedVector359 + crRightHandSideBoundedVector38);
const double crRightHandSideBoundedVector361 =             0.166666666666667*ResProj_1_3;
const double crRightHandSideBoundedVector362 =             0.166666666666667*dUdt_1_3;
const double crRightHandSideBoundedVector363 =             0.166666666666667*ResProj_0_3;
const double crRightHandSideBoundedVector364 =             0.166666666666667*dUdt_0_3;
const double crRightHandSideBoundedVector365 =             0.166666666666667*r[0];
const double crRightHandSideBoundedVector366 =             0.166666666666667*r[1];
const double crRightHandSideBoundedVector367 =             crRightHandSideBoundedVector12*(crRightHandSideBoundedVector365 + crRightHandSideBoundedVector366 + 0.666666666666667*r[2]);
const double crRightHandSideBoundedVector368 =             crRightHandSideBoundedVector18*crRightHandSideBoundedVector48;
const double crRightHandSideBoundedVector369 =             crRightHandSideBoundedVector158*crRightHandSideBoundedVector23;
const double crRightHandSideBoundedVector370 =             crRightHandSideBoundedVector13*crRightHandSideBoundedVector44*gamma;
const double crRightHandSideBoundedVector371 =             crRightHandSideBoundedVector13*crRightHandSideBoundedVector154*gamma;
const double crRightHandSideBoundedVector372 =             crRightHandSideBoundedVector13*crRightHandSideBoundedVector19;
const double crRightHandSideBoundedVector373 =             crRightHandSideBoundedVector310*crRightHandSideBoundedVector372;
const double crRightHandSideBoundedVector374 =             crRightHandSideBoundedVector29 + crRightHandSideBoundedVector31 + crRightHandSideBoundedVector33;
const double crRightHandSideBoundedVector375 =             crRightHandSideBoundedVector27*(-crRightHandSideBoundedVector13*(crRightHandSideBoundedVector35 + crRightHandSideBoundedVector36) + crRightHandSideBoundedVector374);
const double crRightHandSideBoundedVector376 =             crRightHandSideBoundedVector29 + crRightHandSideBoundedVector31 + crRightHandSideBoundedVector33 + crRightHandSideBoundedVector375;
const double crRightHandSideBoundedVector377 =             crRightHandSideBoundedVector13*crRightHandSideBoundedVector24;
const double crRightHandSideBoundedVector378 =             crRightHandSideBoundedVector310*crRightHandSideBoundedVector377;
const double crRightHandSideBoundedVector379 =             crRightHandSideBoundedVector30 + crRightHandSideBoundedVector32 + crRightHandSideBoundedVector34 - crRightHandSideBoundedVector375;
const double crRightHandSideBoundedVector380 =             crRightHandSideBoundedVector13*crRightHandSideBoundedVector25*crRightHandSideBoundedVector67 + crRightHandSideBoundedVector379;
const double crRightHandSideBoundedVector381 =             crRightHandSideBoundedVector275*crRightHandSideBoundedVector380*crRightHandSideBoundedVector66;
const double crRightHandSideBoundedVector382 =             crRightHandSideBoundedVector267*crRightHandSideBoundedVector380*crRightHandSideBoundedVector63;
const double crRightHandSideBoundedVector383 =             0.666666666666667*ResProj_2_3 + crRightHandSideBoundedVector13*crRightHandSideBoundedVector50*(crRightHandSideBoundedVector376 - crRightHandSideBoundedVector378) + crRightHandSideBoundedVector13*crRightHandSideBoundedVector61*(-crRightHandSideBoundedVector373 + crRightHandSideBoundedVector376) - crRightHandSideBoundedVector162*crRightHandSideBoundedVector64 + crRightHandSideBoundedVector18*crRightHandSideBoundedVector370 + crRightHandSideBoundedVector23*crRightHandSideBoundedVector371 + crRightHandSideBoundedVector361 + crRightHandSideBoundedVector362 + crRightHandSideBoundedVector363 + crRightHandSideBoundedVector364 - crRightHandSideBoundedVector367 - crRightHandSideBoundedVector368 - crRightHandSideBoundedVector369 + crRightHandSideBoundedVector381 + crRightHandSideBoundedVector382 - crRightHandSideBoundedVector57*crRightHandSideBoundedVector64 + 0.666666666666667*dUdt_2_3;
const double crRightHandSideBoundedVector384 =             crRightHandSideBoundedVector360*crRightHandSideBoundedVector383;
const double crRightHandSideBoundedVector385 =             1.0/(crRightHandSideBoundedVector359*crRightHandSideBoundedVector76 + crRightHandSideBoundedVector95);
const double crRightHandSideBoundedVector386 =             0.166666666666667*ResProj_2_3 + 0.166666666666667*dUdt_2_3;
const double crRightHandSideBoundedVector387 =             0.166666666666667*r[2];
const double crRightHandSideBoundedVector388 =             crRightHandSideBoundedVector75*(crRightHandSideBoundedVector365 + crRightHandSideBoundedVector387 + 0.666666666666667*r[1]);
const double crRightHandSideBoundedVector389 =             crRightHandSideBoundedVector80*crRightHandSideBoundedVector99;
const double crRightHandSideBoundedVector390 =             crRightHandSideBoundedVector172*crRightHandSideBoundedVector84;
const double crRightHandSideBoundedVector391 =             crRightHandSideBoundedVector44*crRightHandSideBoundedVector76*gamma;
const double crRightHandSideBoundedVector392 =             crRightHandSideBoundedVector154*crRightHandSideBoundedVector76*gamma;
const double crRightHandSideBoundedVector393 =             crRightHandSideBoundedVector76*crRightHandSideBoundedVector81;
const double crRightHandSideBoundedVector394 =             crRightHandSideBoundedVector310*crRightHandSideBoundedVector393;
const double crRightHandSideBoundedVector395 =             crRightHandSideBoundedVector29 + crRightHandSideBoundedVector88 + crRightHandSideBoundedVector90;
const double crRightHandSideBoundedVector396 =             crRightHandSideBoundedVector27*(crRightHandSideBoundedVector395 - crRightHandSideBoundedVector76*(crRightHandSideBoundedVector92 + crRightHandSideBoundedVector93));
const double crRightHandSideBoundedVector397 =             crRightHandSideBoundedVector29 + crRightHandSideBoundedVector396 + crRightHandSideBoundedVector88 + crRightHandSideBoundedVector90;
const double crRightHandSideBoundedVector398 =             crRightHandSideBoundedVector76*crRightHandSideBoundedVector85;
const double crRightHandSideBoundedVector399 =             crRightHandSideBoundedVector310*crRightHandSideBoundedVector398;
const double crRightHandSideBoundedVector400 =             crRightHandSideBoundedVector30 - crRightHandSideBoundedVector396 + crRightHandSideBoundedVector89 + crRightHandSideBoundedVector91;
const double crRightHandSideBoundedVector401 =             crRightHandSideBoundedVector400 + crRightHandSideBoundedVector67*crRightHandSideBoundedVector76*crRightHandSideBoundedVector86;
const double crRightHandSideBoundedVector402 =             crRightHandSideBoundedVector292*crRightHandSideBoundedVector401*crRightHandSideBoundedVector66;
const double crRightHandSideBoundedVector403 =             crRightHandSideBoundedVector285*crRightHandSideBoundedVector401*crRightHandSideBoundedVector63;
const double crRightHandSideBoundedVector404 =             0.666666666666667*ResProj_1_3 - crRightHandSideBoundedVector107*crRightHandSideBoundedVector162 - crRightHandSideBoundedVector107*crRightHandSideBoundedVector57 + crRightHandSideBoundedVector363 + crRightHandSideBoundedVector364 + crRightHandSideBoundedVector386 - crRightHandSideBoundedVector388 - crRightHandSideBoundedVector389 - crRightHandSideBoundedVector390 + crRightHandSideBoundedVector391*crRightHandSideBoundedVector80 + crRightHandSideBoundedVector392*crRightHandSideBoundedVector84 + crRightHandSideBoundedVector402 + crRightHandSideBoundedVector403 + crRightHandSideBoundedVector50*crRightHandSideBoundedVector76*(crRightHandSideBoundedVector397 - crRightHandSideBoundedVector399) + crRightHandSideBoundedVector61*crRightHandSideBoundedVector76*(-crRightHandSideBoundedVector394 + crRightHandSideBoundedVector397) + 0.666666666666667*dUdt_1_3;
const double crRightHandSideBoundedVector405 =             crRightHandSideBoundedVector385*crRightHandSideBoundedVector404;
const double crRightHandSideBoundedVector406 =             1.0/(crRightHandSideBoundedVector116*crRightHandSideBoundedVector359 + crRightHandSideBoundedVector131);
const double crRightHandSideBoundedVector407 =             crRightHandSideBoundedVector115*(crRightHandSideBoundedVector366 + crRightHandSideBoundedVector387 + 0.666666666666667*r[0]);
const double crRightHandSideBoundedVector408 =             crRightHandSideBoundedVector119*crRightHandSideBoundedVector133;
const double crRightHandSideBoundedVector409 =             crRightHandSideBoundedVector122*crRightHandSideBoundedVector183;
const double crRightHandSideBoundedVector410 =             crRightHandSideBoundedVector116*crRightHandSideBoundedVector44*gamma;
const double crRightHandSideBoundedVector411 =             crRightHandSideBoundedVector116*crRightHandSideBoundedVector154*gamma;
const double crRightHandSideBoundedVector412 =             crRightHandSideBoundedVector116*crRightHandSideBoundedVector120;
const double crRightHandSideBoundedVector413 =             crRightHandSideBoundedVector310*crRightHandSideBoundedVector412;
const double crRightHandSideBoundedVector414 =             crRightHandSideBoundedVector126 + crRightHandSideBoundedVector31 + crRightHandSideBoundedVector90;
const double crRightHandSideBoundedVector415 =             crRightHandSideBoundedVector27*(-crRightHandSideBoundedVector116*(crRightHandSideBoundedVector128 + crRightHandSideBoundedVector129) + crRightHandSideBoundedVector414);
const double crRightHandSideBoundedVector416 =             crRightHandSideBoundedVector126 + crRightHandSideBoundedVector31 + crRightHandSideBoundedVector415 + crRightHandSideBoundedVector90;
const double crRightHandSideBoundedVector417 =             crRightHandSideBoundedVector116*crRightHandSideBoundedVector123;
const double crRightHandSideBoundedVector418 =             crRightHandSideBoundedVector310*crRightHandSideBoundedVector417;
const double crRightHandSideBoundedVector419 =             crRightHandSideBoundedVector127 + crRightHandSideBoundedVector32 - crRightHandSideBoundedVector415 + crRightHandSideBoundedVector91;
const double crRightHandSideBoundedVector420 =             crRightHandSideBoundedVector116*crRightHandSideBoundedVector124*crRightHandSideBoundedVector67 + crRightHandSideBoundedVector419;
const double crRightHandSideBoundedVector421 =             crRightHandSideBoundedVector214*crRightHandSideBoundedVector420*crRightHandSideBoundedVector66;
const double crRightHandSideBoundedVector422 =             crRightHandSideBoundedVector304*crRightHandSideBoundedVector420*crRightHandSideBoundedVector63;
const double crRightHandSideBoundedVector423 =             0.666666666666667*ResProj_0_3 + crRightHandSideBoundedVector116*crRightHandSideBoundedVector50*(crRightHandSideBoundedVector416 - crRightHandSideBoundedVector418) + crRightHandSideBoundedVector116*crRightHandSideBoundedVector61*(-crRightHandSideBoundedVector413 + crRightHandSideBoundedVector416) + crRightHandSideBoundedVector119*crRightHandSideBoundedVector410 + crRightHandSideBoundedVector122*crRightHandSideBoundedVector411 - crRightHandSideBoundedVector141*crRightHandSideBoundedVector162 - crRightHandSideBoundedVector141*crRightHandSideBoundedVector57 + crRightHandSideBoundedVector361 + crRightHandSideBoundedVector362 + crRightHandSideBoundedVector386 - crRightHandSideBoundedVector407 - crRightHandSideBoundedVector408 - crRightHandSideBoundedVector409 + crRightHandSideBoundedVector421 + crRightHandSideBoundedVector422 + 0.666666666666667*dUdt_0_3;
const double crRightHandSideBoundedVector424 =             crRightHandSideBoundedVector406*crRightHandSideBoundedVector423;
const double crRightHandSideBoundedVector425 =             0.166666666666667*crRightHandSideBoundedVector173;
const double crRightHandSideBoundedVector426 =             -0.166666666666667*crRightHandSideBoundedVector174;
const double crRightHandSideBoundedVector427 =             -0.166666666666667*crRightHandSideBoundedVector175;
const double crRightHandSideBoundedVector428 =             crRightHandSideBoundedVector27*crRightHandSideBoundedVector53;
const double crRightHandSideBoundedVector429 =             0.166666666666667*crRightHandSideBoundedVector428;
const double crRightHandSideBoundedVector430 =             crRightHandSideBoundedVector101*crRightHandSideBoundedVector429;
const double crRightHandSideBoundedVector431 =             crRightHandSideBoundedVector200*crRightHandSideBoundedVector50;
const double crRightHandSideBoundedVector432 =             0.166666666666667*crRightHandSideBoundedVector431;
const double crRightHandSideBoundedVector433 =             crRightHandSideBoundedVector103*crRightHandSideBoundedVector432;
const double crRightHandSideBoundedVector434 =             0.166666666666667*crRightHandSideBoundedVector177;
const double crRightHandSideBoundedVector435 =             -0.166666666666667*crRightHandSideBoundedVector180;
const double crRightHandSideBoundedVector436 =             0.166666666666667*crRightHandSideBoundedVector159;
const double crRightHandSideBoundedVector437 =             -0.166666666666667*crRightHandSideBoundedVector160;
const double crRightHandSideBoundedVector438 =             -0.166666666666667*crRightHandSideBoundedVector161;
const double crRightHandSideBoundedVector439 =             crRightHandSideBoundedVector429*crRightHandSideBoundedVector51;
const double crRightHandSideBoundedVector440 =             0.666666666666667*crRightHandSideBoundedVector428;
const double crRightHandSideBoundedVector441 =             crRightHandSideBoundedVector432*crRightHandSideBoundedVector54;
const double crRightHandSideBoundedVector442 =             0.666666666666667*crRightHandSideBoundedVector431;
const double crRightHandSideBoundedVector443 =             0.166666666666667*crRightHandSideBoundedVector164;
const double crRightHandSideBoundedVector444 =             crRightHandSideBoundedVector119*crRightHandSideBoundedVector218;
const double crRightHandSideBoundedVector445 =             -0.166666666666667*crRightHandSideBoundedVector167;
const double crRightHandSideBoundedVector446 =             DN_DX_0_0*crRightHandSideBoundedVector219*mu;
const double crRightHandSideBoundedVector447 =             pow(lin_m[1], 2);
const double crRightHandSideBoundedVector448 =             crRightHandSideBoundedVector221*crRightHandSideBoundedVector222*crRightHandSideBoundedVector447*nu_st;
const double crRightHandSideBoundedVector449 =             -crRightHandSideBoundedVector221*crRightHandSideBoundedVector447 + 1;
const double crRightHandSideBoundedVector450 =             crRightHandSideBoundedVector222*crRightHandSideBoundedVector449*nu_sc;
const double crRightHandSideBoundedVector451 =             crRightHandSideBoundedVector233*crRightHandSideBoundedVector245 + crRightHandSideBoundedVector239*mu*(crRightHandSideBoundedVector12*crRightHandSideBoundedVector448 + crRightHandSideBoundedVector12*crRightHandSideBoundedVector450 + 1);
const double crRightHandSideBoundedVector452 =             crRightHandSideBoundedVector247*crRightHandSideBoundedVector252 + crRightHandSideBoundedVector251*mu*(crRightHandSideBoundedVector448*crRightHandSideBoundedVector75 + crRightHandSideBoundedVector450*crRightHandSideBoundedVector75 + 1);
const double crRightHandSideBoundedVector453 =             crRightHandSideBoundedVector254*crRightHandSideBoundedVector258 + crRightHandSideBoundedVector257*mu*(crRightHandSideBoundedVector115*crRightHandSideBoundedVector448 + crRightHandSideBoundedVector115*crRightHandSideBoundedVector450 + 1);
const double crRightHandSideBoundedVector454 =             crRightHandSideBoundedVector267*crRightHandSideBoundedVector332;
const double crRightHandSideBoundedVector455 =             0.0277777777777778*f_ext(0,1);
const double crRightHandSideBoundedVector456 =             0.0277777777777778*f_ext(1,1);
const double crRightHandSideBoundedVector457 =             0.111111111111111*f_ext(2,1);
const double crRightHandSideBoundedVector458 =             crRightHandSideBoundedVector14*crRightHandSideBoundedVector18*crRightHandSideBoundedVector56;
const double crRightHandSideBoundedVector459 =             crRightHandSideBoundedVector14*crRightHandSideBoundedVector23*crRightHandSideBoundedVector61;
const double crRightHandSideBoundedVector460 =             0.166666666666667*crRightHandSideBoundedVector459;
const double crRightHandSideBoundedVector461 =             crRightHandSideBoundedVector275*crRightHandSideBoundedVector429;
const double crRightHandSideBoundedVector462 =             crRightHandSideBoundedVector18*crRightHandSideBoundedVector23*crRightHandSideBoundedVector276*crRightHandSideBoundedVector66;
const double crRightHandSideBoundedVector463 =             0.333333333333333*crRightHandSideBoundedVector462;
const double crRightHandSideBoundedVector464 =             crRightHandSideBoundedVector165*crRightHandSideBoundedVector276*crRightHandSideBoundedVector63;
const double crRightHandSideBoundedVector465 =             crRightHandSideBoundedVector267*crRightHandSideBoundedVector432 + crRightHandSideBoundedVector455 + crRightHandSideBoundedVector456 + crRightHandSideBoundedVector457 - 0.166666666666667*crRightHandSideBoundedVector458 - crRightHandSideBoundedVector460 + crRightHandSideBoundedVector461 + crRightHandSideBoundedVector463 - 0.333333333333333*crRightHandSideBoundedVector464;
const double crRightHandSideBoundedVector466 =             crRightHandSideBoundedVector285*crRightHandSideBoundedVector343;
const double crRightHandSideBoundedVector467 =             0.111111111111111*f_ext(1,1);
const double crRightHandSideBoundedVector468 =             0.0277777777777778*f_ext(2,1);
const double crRightHandSideBoundedVector469 =             crRightHandSideBoundedVector56*crRightHandSideBoundedVector77*crRightHandSideBoundedVector80;
const double crRightHandSideBoundedVector470 =             crRightHandSideBoundedVector61*crRightHandSideBoundedVector77*crRightHandSideBoundedVector84;
const double crRightHandSideBoundedVector471 =             0.166666666666667*crRightHandSideBoundedVector470;
const double crRightHandSideBoundedVector472 =             crRightHandSideBoundedVector292*crRightHandSideBoundedVector429;
const double crRightHandSideBoundedVector473 =             crRightHandSideBoundedVector293*crRightHandSideBoundedVector66*crRightHandSideBoundedVector80*crRightHandSideBoundedVector84;
const double crRightHandSideBoundedVector474 =             0.333333333333333*crRightHandSideBoundedVector473;
const double crRightHandSideBoundedVector475 =             crRightHandSideBoundedVector178*crRightHandSideBoundedVector293*crRightHandSideBoundedVector63;
const double crRightHandSideBoundedVector476 =             crRightHandSideBoundedVector285*crRightHandSideBoundedVector432 + crRightHandSideBoundedVector455 + crRightHandSideBoundedVector467 + crRightHandSideBoundedVector468 - 0.166666666666667*crRightHandSideBoundedVector469 - crRightHandSideBoundedVector471 + crRightHandSideBoundedVector472 + crRightHandSideBoundedVector474 - 0.333333333333333*crRightHandSideBoundedVector475;
const double crRightHandSideBoundedVector477 =             crRightHandSideBoundedVector457 + crRightHandSideBoundedVector467 + 0.444444444444444*f_ext(0,1);
const double crRightHandSideBoundedVector478 =             crRightHandSideBoundedVector117*crRightHandSideBoundedVector119*crRightHandSideBoundedVector56;
const double crRightHandSideBoundedVector479 =             crRightHandSideBoundedVector117*crRightHandSideBoundedVector122*crRightHandSideBoundedVector61;
const double crRightHandSideBoundedVector480 =             0.666666666666667*crRightHandSideBoundedVector479;
const double crRightHandSideBoundedVector481 =             crRightHandSideBoundedVector304*crRightHandSideBoundedVector354;
const double crRightHandSideBoundedVector482 =             crRightHandSideBoundedVector214*crRightHandSideBoundedVector440;
const double crRightHandSideBoundedVector483 =             crRightHandSideBoundedVector119*crRightHandSideBoundedVector122*crRightHandSideBoundedVector306*crRightHandSideBoundedVector66;
const double crRightHandSideBoundedVector484 =             1.33333333333333*crRightHandSideBoundedVector483;
const double crRightHandSideBoundedVector485 =             crRightHandSideBoundedVector189*crRightHandSideBoundedVector306*crRightHandSideBoundedVector63;
const double crRightHandSideBoundedVector486 =             DN_DX_0_0*crRightHandSideBoundedVector20;
const double crRightHandSideBoundedVector487 =             DN_DX_1_0*crRightHandSideBoundedVector21;
const double crRightHandSideBoundedVector488 =             DN_DX_2_0*crRightHandSideBoundedVector83;
const double crRightHandSideBoundedVector489 =             -crRightHandSideBoundedVector429;
const double crRightHandSideBoundedVector490 =             crRightHandSideBoundedVector27*crRightHandSideBoundedVector316 - crRightHandSideBoundedVector318 + crRightHandSideBoundedVector486 + crRightHandSideBoundedVector487 + crRightHandSideBoundedVector488 + crRightHandSideBoundedVector489;
const double crRightHandSideBoundedVector491 =             crRightHandSideBoundedVector27*crRightHandSideBoundedVector323 - crRightHandSideBoundedVector325 + crRightHandSideBoundedVector486 + crRightHandSideBoundedVector487 + crRightHandSideBoundedVector488 + crRightHandSideBoundedVector489;
const double crRightHandSideBoundedVector492 =             DN_DX_0_0*crRightHandSideBoundedVector121 + DN_DX_1_0*crRightHandSideBoundedVector82 + DN_DX_2_0*crRightHandSideBoundedVector22 - crRightHandSideBoundedVector440;
const double crRightHandSideBoundedVector493 =             DN_DX_0_0*crRightHandSideBoundedVector15;
const double crRightHandSideBoundedVector494 =             DN_DX_1_0*crRightHandSideBoundedVector16;
const double crRightHandSideBoundedVector495 =             DN_DX_2_0*crRightHandSideBoundedVector79;
const double crRightHandSideBoundedVector496 =             -crRightHandSideBoundedVector432;
const double crRightHandSideBoundedVector497 =             -crRightHandSideBoundedVector339;
const double crRightHandSideBoundedVector498 =             crRightHandSideBoundedVector200*crRightHandSideBoundedVector337 + crRightHandSideBoundedVector493 + crRightHandSideBoundedVector494 + crRightHandSideBoundedVector495 + crRightHandSideBoundedVector496 + crRightHandSideBoundedVector497;
const double crRightHandSideBoundedVector499 =             -crRightHandSideBoundedVector346;
const double crRightHandSideBoundedVector500 =             crRightHandSideBoundedVector200*crRightHandSideBoundedVector344 + crRightHandSideBoundedVector493 + crRightHandSideBoundedVector494 + crRightHandSideBoundedVector495 + crRightHandSideBoundedVector496 + crRightHandSideBoundedVector499;
const double crRightHandSideBoundedVector501 =             DN_DX_0_0*crRightHandSideBoundedVector118;
const double crRightHandSideBoundedVector502 =             DN_DX_1_0*crRightHandSideBoundedVector78;
const double crRightHandSideBoundedVector503 =             DN_DX_2_0*crRightHandSideBoundedVector17;
const double crRightHandSideBoundedVector504 =             -crRightHandSideBoundedVector442 + crRightHandSideBoundedVector501 + crRightHandSideBoundedVector502 + crRightHandSideBoundedVector503;
const double crRightHandSideBoundedVector505 =             -crRightHandSideBoundedVector116*crRightHandSideBoundedVector444;
const double crRightHandSideBoundedVector506 =             DN_DX_0_1*crRightHandSideBoundedVector310*h;
const double crRightHandSideBoundedVector507 =             0.166666666666667*crRightHandSideBoundedVector388;
const double crRightHandSideBoundedVector508 =             0.166666666666667*crRightHandSideBoundedVector389;
const double crRightHandSideBoundedVector509 =             0.166666666666667*crRightHandSideBoundedVector390;
const double crRightHandSideBoundedVector510 =             0.166666666666667*crRightHandSideBoundedVector391;
const double crRightHandSideBoundedVector511 =             -crRightHandSideBoundedVector510*crRightHandSideBoundedVector80;
const double crRightHandSideBoundedVector512 =             0.166666666666667*crRightHandSideBoundedVector392;
const double crRightHandSideBoundedVector513 =             -crRightHandSideBoundedVector512*crRightHandSideBoundedVector84;
const double crRightHandSideBoundedVector514 =             crRightHandSideBoundedVector291*crRightHandSideBoundedVector80;
const double crRightHandSideBoundedVector515 =             crRightHandSideBoundedVector285*crRightHandSideBoundedVector429;
const double crRightHandSideBoundedVector516 =             crRightHandSideBoundedVector515*crRightHandSideBoundedVector80;
const double crRightHandSideBoundedVector517 =             crRightHandSideBoundedVector394 + crRightHandSideBoundedVector400;
const double crRightHandSideBoundedVector518 =             crRightHandSideBoundedVector517*crRightHandSideBoundedVector61*crRightHandSideBoundedVector76;
const double crRightHandSideBoundedVector519 =             0.166666666666667*crRightHandSideBoundedVector518;
const double crRightHandSideBoundedVector520 =             crRightHandSideBoundedVector399 + crRightHandSideBoundedVector400;
const double crRightHandSideBoundedVector521 =             crRightHandSideBoundedVector50*crRightHandSideBoundedVector520*crRightHandSideBoundedVector76;
const double crRightHandSideBoundedVector522 =             0.166666666666667*crRightHandSideBoundedVector521;
const double crRightHandSideBoundedVector523 =             -0.166666666666667*crRightHandSideBoundedVector402;
const double crRightHandSideBoundedVector524 =             -0.166666666666667*crRightHandSideBoundedVector403;
const double crRightHandSideBoundedVector525 =             0.166666666666667*crRightHandSideBoundedVector367;
const double crRightHandSideBoundedVector526 =             0.166666666666667*crRightHandSideBoundedVector368;
const double crRightHandSideBoundedVector527 =             0.166666666666667*crRightHandSideBoundedVector369;
const double crRightHandSideBoundedVector528 =             0.166666666666667*crRightHandSideBoundedVector370;
const double crRightHandSideBoundedVector529 =             -crRightHandSideBoundedVector18*crRightHandSideBoundedVector528;
const double crRightHandSideBoundedVector530 =             0.166666666666667*crRightHandSideBoundedVector371;
const double crRightHandSideBoundedVector531 =             -crRightHandSideBoundedVector23*crRightHandSideBoundedVector530;
const double crRightHandSideBoundedVector532 =             0.666666666666667*crRightHandSideBoundedVector410;
const double crRightHandSideBoundedVector533 =             0.666666666666667*crRightHandSideBoundedVector411;
const double crRightHandSideBoundedVector534 =             crRightHandSideBoundedVector18*crRightHandSideBoundedVector274;
const double crRightHandSideBoundedVector535 =             crRightHandSideBoundedVector267*crRightHandSideBoundedVector429;
const double crRightHandSideBoundedVector536 =             crRightHandSideBoundedVector18*crRightHandSideBoundedVector535;
const double crRightHandSideBoundedVector537 =             crRightHandSideBoundedVector304*crRightHandSideBoundedVector440;
const double crRightHandSideBoundedVector538 =             crRightHandSideBoundedVector373 + crRightHandSideBoundedVector379;
const double crRightHandSideBoundedVector539 =             crRightHandSideBoundedVector13*crRightHandSideBoundedVector538*crRightHandSideBoundedVector61;
const double crRightHandSideBoundedVector540 =             0.166666666666667*crRightHandSideBoundedVector539;
const double crRightHandSideBoundedVector541 =             crRightHandSideBoundedVector378 + crRightHandSideBoundedVector379;
const double crRightHandSideBoundedVector542 =             crRightHandSideBoundedVector13*crRightHandSideBoundedVector50*crRightHandSideBoundedVector541;
const double crRightHandSideBoundedVector543 =             0.166666666666667*crRightHandSideBoundedVector542;
const double crRightHandSideBoundedVector544 =             crRightHandSideBoundedVector413 + crRightHandSideBoundedVector419;
const double crRightHandSideBoundedVector545 =             crRightHandSideBoundedVector116*crRightHandSideBoundedVector544*crRightHandSideBoundedVector61;
const double crRightHandSideBoundedVector546 =             crRightHandSideBoundedVector418 + crRightHandSideBoundedVector419;
const double crRightHandSideBoundedVector547 =             crRightHandSideBoundedVector116*crRightHandSideBoundedVector50*crRightHandSideBoundedVector546;
const double crRightHandSideBoundedVector548 =             -0.166666666666667*crRightHandSideBoundedVector381;
const double crRightHandSideBoundedVector549 =             -0.166666666666667*crRightHandSideBoundedVector382;
const double crRightHandSideBoundedVector550 =             crRightHandSideBoundedVector117*crRightHandSideBoundedVector420;
const double crRightHandSideBoundedVector551 =             crRightHandSideBoundedVector222*nu_sc*(-crRightHandSideBoundedVector225 + 1);
const double crRightHandSideBoundedVector552 =             crRightHandSideBoundedVector219*mu*(crRightHandSideBoundedVector12*crRightHandSideBoundedVector551 + crRightHandSideBoundedVector224);
const double crRightHandSideBoundedVector553 =             crRightHandSideBoundedVector221*lin_m[0]*lin_m[1]*(-k_sc + k_st);
const double crRightHandSideBoundedVector554 =             crRightHandSideBoundedVector12*crRightHandSideBoundedVector154 - crRightHandSideBoundedVector13*(crRightHandSideBoundedVector163 + crRightHandSideBoundedVector18*crRightHandSideBoundedVector53 - crRightHandSideBoundedVector25*crRightHandSideBoundedVector63) - crRightHandSideBoundedVector374*crRightHandSideBoundedVector63;
const double crRightHandSideBoundedVector555 =             crRightHandSideBoundedVector13*crRightHandSideBoundedVector358*lambda;
const double crRightHandSideBoundedVector556 =             1.0/lambda;
const double crRightHandSideBoundedVector557 =             c_v*crRightHandSideBoundedVector221*crRightHandSideBoundedVector240*crRightHandSideBoundedVector556*k_st;
const double crRightHandSideBoundedVector558 =             c_v*crRightHandSideBoundedVector242*crRightHandSideBoundedVector556*k_sc;
const double crRightHandSideBoundedVector559 =             crRightHandSideBoundedVector12*crRightHandSideBoundedVector44 - crRightHandSideBoundedVector13*(crRightHandSideBoundedVector23*crRightHandSideBoundedVector56 - crRightHandSideBoundedVector25*crRightHandSideBoundedVector66 + crRightHandSideBoundedVector62) - crRightHandSideBoundedVector374*crRightHandSideBoundedVector66;
const double crRightHandSideBoundedVector560 =             crRightHandSideBoundedVector13*(crRightHandSideBoundedVector18*crRightHandSideBoundedVector246 + crRightHandSideBoundedVector23*crRightHandSideBoundedVector552 - crRightHandSideBoundedVector553*crRightHandSideBoundedVector554 - crRightHandSideBoundedVector555*crRightHandSideBoundedVector559*(crRightHandSideBoundedVector12*crRightHandSideBoundedVector557 + crRightHandSideBoundedVector12*crRightHandSideBoundedVector558 + 1));
const double crRightHandSideBoundedVector561 =             crRightHandSideBoundedVector219*mu*(crRightHandSideBoundedVector228 + crRightHandSideBoundedVector551*crRightHandSideBoundedVector75);
const double crRightHandSideBoundedVector562 =             crRightHandSideBoundedVector154*crRightHandSideBoundedVector75 - crRightHandSideBoundedVector395*crRightHandSideBoundedVector63 - crRightHandSideBoundedVector76*(crRightHandSideBoundedVector176 + crRightHandSideBoundedVector53*crRightHandSideBoundedVector80 - crRightHandSideBoundedVector63*crRightHandSideBoundedVector86);
const double crRightHandSideBoundedVector563 =             crRightHandSideBoundedVector358*crRightHandSideBoundedVector76*lambda;
const double crRightHandSideBoundedVector564 =             -crRightHandSideBoundedVector395*crRightHandSideBoundedVector66 + crRightHandSideBoundedVector44*crRightHandSideBoundedVector75 - crRightHandSideBoundedVector76*(crRightHandSideBoundedVector106 + crRightHandSideBoundedVector56*crRightHandSideBoundedVector84 - crRightHandSideBoundedVector66*crRightHandSideBoundedVector86);
const double crRightHandSideBoundedVector565 =             crRightHandSideBoundedVector76*(crRightHandSideBoundedVector253*crRightHandSideBoundedVector80 - crRightHandSideBoundedVector553*crRightHandSideBoundedVector562 + crRightHandSideBoundedVector561*crRightHandSideBoundedVector84 - crRightHandSideBoundedVector563*crRightHandSideBoundedVector564*(crRightHandSideBoundedVector557*crRightHandSideBoundedVector75 + crRightHandSideBoundedVector558*crRightHandSideBoundedVector75 + 1));
const double crRightHandSideBoundedVector566 =             crRightHandSideBoundedVector219*mu*(crRightHandSideBoundedVector115*crRightHandSideBoundedVector551 + crRightHandSideBoundedVector230);
const double crRightHandSideBoundedVector567 =             crRightHandSideBoundedVector115*crRightHandSideBoundedVector154 - crRightHandSideBoundedVector116*(crRightHandSideBoundedVector119*crRightHandSideBoundedVector53 - crRightHandSideBoundedVector124*crRightHandSideBoundedVector63 + crRightHandSideBoundedVector187) - crRightHandSideBoundedVector414*crRightHandSideBoundedVector63;
const double crRightHandSideBoundedVector568 =             crRightHandSideBoundedVector116*crRightHandSideBoundedVector358*lambda;
const double crRightHandSideBoundedVector569 =             crRightHandSideBoundedVector115*crRightHandSideBoundedVector44 - crRightHandSideBoundedVector116*(crRightHandSideBoundedVector122*crRightHandSideBoundedVector56 - crRightHandSideBoundedVector124*crRightHandSideBoundedVector66 + crRightHandSideBoundedVector140) - crRightHandSideBoundedVector414*crRightHandSideBoundedVector66;
const double crRightHandSideBoundedVector570 =             crRightHandSideBoundedVector116*(crRightHandSideBoundedVector119*crRightHandSideBoundedVector259 + crRightHandSideBoundedVector122*crRightHandSideBoundedVector566 - crRightHandSideBoundedVector553*crRightHandSideBoundedVector567 - crRightHandSideBoundedVector568*crRightHandSideBoundedVector569*(crRightHandSideBoundedVector115*crRightHandSideBoundedVector557 + crRightHandSideBoundedVector115*crRightHandSideBoundedVector558 + 1));
const double crRightHandSideBoundedVector571 =             c_v*crRightHandSideBoundedVector221*crRightHandSideBoundedVector447*crRightHandSideBoundedVector556*k_st;
const double crRightHandSideBoundedVector572 =             c_v*crRightHandSideBoundedVector449*crRightHandSideBoundedVector556*k_sc;
const double crRightHandSideBoundedVector573 =             crRightHandSideBoundedVector13*(crRightHandSideBoundedVector18*crRightHandSideBoundedVector552 + crRightHandSideBoundedVector23*crRightHandSideBoundedVector451 - crRightHandSideBoundedVector553*crRightHandSideBoundedVector559 - crRightHandSideBoundedVector554*crRightHandSideBoundedVector555*(crRightHandSideBoundedVector12*crRightHandSideBoundedVector571 + crRightHandSideBoundedVector12*crRightHandSideBoundedVector572 + 1));
const double crRightHandSideBoundedVector574 =             crRightHandSideBoundedVector76*(crRightHandSideBoundedVector452*crRightHandSideBoundedVector84 - crRightHandSideBoundedVector553*crRightHandSideBoundedVector564 + crRightHandSideBoundedVector561*crRightHandSideBoundedVector80 - crRightHandSideBoundedVector562*crRightHandSideBoundedVector563*(crRightHandSideBoundedVector571*crRightHandSideBoundedVector75 + crRightHandSideBoundedVector572*crRightHandSideBoundedVector75 + 1));
const double crRightHandSideBoundedVector575 =             crRightHandSideBoundedVector116*(crRightHandSideBoundedVector119*crRightHandSideBoundedVector566 + crRightHandSideBoundedVector122*crRightHandSideBoundedVector453 - crRightHandSideBoundedVector553*crRightHandSideBoundedVector569 - crRightHandSideBoundedVector567*crRightHandSideBoundedVector568*(crRightHandSideBoundedVector115*crRightHandSideBoundedVector571 + crRightHandSideBoundedVector115*crRightHandSideBoundedVector572 + 1));
const double crRightHandSideBoundedVector576 =             1.0*crRightHandSideBoundedVector39*crRightHandSideBoundedVector72*h;
const double crRightHandSideBoundedVector577 =             crRightHandSideBoundedVector13*crRightHandSideBoundedVector538;
const double crRightHandSideBoundedVector578 =             crRightHandSideBoundedVector14*crRightHandSideBoundedVector18*crRightHandSideBoundedVector61;
const double crRightHandSideBoundedVector579 =             2.0*gamma - 2.0;
const double crRightHandSideBoundedVector580 =             crRightHandSideBoundedVector372*crRightHandSideBoundedVector579 + crRightHandSideBoundedVector380;
const double crRightHandSideBoundedVector581 =             0.166666666666667*crRightHandSideBoundedVector14*crRightHandSideBoundedVector580;
const double crRightHandSideBoundedVector582 =             crRightHandSideBoundedVector268 + crRightHandSideBoundedVector269 - crRightHandSideBoundedVector27*crRightHandSideBoundedVector272 + crRightHandSideBoundedVector27*crRightHandSideBoundedVector278 + crRightHandSideBoundedVector270 - crRightHandSideBoundedVector274 + crRightHandSideBoundedVector528 - crRightHandSideBoundedVector535 - crRightHandSideBoundedVector578*crRightHandSideBoundedVector67 + crRightHandSideBoundedVector581*crRightHandSideBoundedVector66;
const double crRightHandSideBoundedVector583 =             1.0*crRightHandSideBoundedVector168*crRightHandSideBoundedVector39*h;
const double crRightHandSideBoundedVector584 =             crRightHandSideBoundedVector13*crRightHandSideBoundedVector541;
const double crRightHandSideBoundedVector585 =             crRightHandSideBoundedVector14*crRightHandSideBoundedVector23*crRightHandSideBoundedVector50;
const double crRightHandSideBoundedVector586 =             crRightHandSideBoundedVector377*crRightHandSideBoundedVector579 + crRightHandSideBoundedVector380;
const double crRightHandSideBoundedVector587 =             0.166666666666667*crRightHandSideBoundedVector14*crRightHandSideBoundedVector586;
const double crRightHandSideBoundedVector588 =             -crRightHandSideBoundedVector198*crRightHandSideBoundedVector275 - crRightHandSideBoundedVector27*crRightHandSideBoundedVector460 + crRightHandSideBoundedVector27*crRightHandSideBoundedVector463 + crRightHandSideBoundedVector455 + crRightHandSideBoundedVector456 + crRightHandSideBoundedVector457 - crRightHandSideBoundedVector461 + crRightHandSideBoundedVector530 - crRightHandSideBoundedVector585*crRightHandSideBoundedVector67 + crRightHandSideBoundedVector587*crRightHandSideBoundedVector63;
const double crRightHandSideBoundedVector589 =             1.0*crRightHandSideBoundedVector113*crRightHandSideBoundedVector96*h;
const double crRightHandSideBoundedVector590 =             crRightHandSideBoundedVector517*crRightHandSideBoundedVector76;
const double crRightHandSideBoundedVector591 =             crRightHandSideBoundedVector61*crRightHandSideBoundedVector77*crRightHandSideBoundedVector80;
const double crRightHandSideBoundedVector592 =             crRightHandSideBoundedVector393*crRightHandSideBoundedVector579 + crRightHandSideBoundedVector401;
const double crRightHandSideBoundedVector593 =             0.166666666666667*crRightHandSideBoundedVector592*crRightHandSideBoundedVector77;
const double crRightHandSideBoundedVector594 =             crRightHandSideBoundedVector268 - crRightHandSideBoundedVector27*crRightHandSideBoundedVector289 + crRightHandSideBoundedVector27*crRightHandSideBoundedVector295 + crRightHandSideBoundedVector286 + crRightHandSideBoundedVector287 - crRightHandSideBoundedVector291 + crRightHandSideBoundedVector510 - crRightHandSideBoundedVector515 - crRightHandSideBoundedVector591*crRightHandSideBoundedVector67 + crRightHandSideBoundedVector593*crRightHandSideBoundedVector66;
const double crRightHandSideBoundedVector595 =             1.0*crRightHandSideBoundedVector181*crRightHandSideBoundedVector96*h;
const double crRightHandSideBoundedVector596 =             crRightHandSideBoundedVector520*crRightHandSideBoundedVector76;
const double crRightHandSideBoundedVector597 =             crRightHandSideBoundedVector50*crRightHandSideBoundedVector77*crRightHandSideBoundedVector84;
const double crRightHandSideBoundedVector598 =             crRightHandSideBoundedVector398*crRightHandSideBoundedVector579 + crRightHandSideBoundedVector401;
const double crRightHandSideBoundedVector599 =             0.166666666666667*crRightHandSideBoundedVector598*crRightHandSideBoundedVector77;
const double crRightHandSideBoundedVector600 =             -crRightHandSideBoundedVector198*crRightHandSideBoundedVector292 - crRightHandSideBoundedVector27*crRightHandSideBoundedVector471 + crRightHandSideBoundedVector27*crRightHandSideBoundedVector474 + crRightHandSideBoundedVector455 + crRightHandSideBoundedVector467 + crRightHandSideBoundedVector468 - crRightHandSideBoundedVector472 + crRightHandSideBoundedVector512 - crRightHandSideBoundedVector597*crRightHandSideBoundedVector67 + crRightHandSideBoundedVector599*crRightHandSideBoundedVector63;
const double crRightHandSideBoundedVector601 =             1.0*crRightHandSideBoundedVector132*crRightHandSideBoundedVector147*h;
const double crRightHandSideBoundedVector602 =             crRightHandSideBoundedVector117*crRightHandSideBoundedVector119*crRightHandSideBoundedVector61;
const double crRightHandSideBoundedVector603 =             crRightHandSideBoundedVector116*crRightHandSideBoundedVector544;
const double crRightHandSideBoundedVector604 =             crRightHandSideBoundedVector412*crRightHandSideBoundedVector579 + crRightHandSideBoundedVector420;
const double crRightHandSideBoundedVector605 =             crRightHandSideBoundedVector117*crRightHandSideBoundedVector604;
const double crRightHandSideBoundedVector606 =             1.0*crRightHandSideBoundedVector132*crRightHandSideBoundedVector192*h;
const double crRightHandSideBoundedVector607 =             crRightHandSideBoundedVector117*crRightHandSideBoundedVector122*crRightHandSideBoundedVector50;
const double crRightHandSideBoundedVector608 =             crRightHandSideBoundedVector116*crRightHandSideBoundedVector546;
const double crRightHandSideBoundedVector609 =             crRightHandSideBoundedVector417*crRightHandSideBoundedVector579 + crRightHandSideBoundedVector420;
const double crRightHandSideBoundedVector610 =             crRightHandSideBoundedVector117*crRightHandSideBoundedVector609;
const double crRightHandSideBoundedVector611 =             crRightHandSideBoundedVector333 + crRightHandSideBoundedVector334 + crRightHandSideBoundedVector335 + crRightHandSideBoundedVector338 + crRightHandSideBoundedVector493 + crRightHandSideBoundedVector494 + crRightHandSideBoundedVector495 + crRightHandSideBoundedVector497;
const double crRightHandSideBoundedVector612 =             1.0*crRightHandSideBoundedVector13*crRightHandSideBoundedVector360*crRightHandSideBoundedVector383*gamma*h;
const double crRightHandSideBoundedVector613 =             crRightHandSideBoundedVector333 + crRightHandSideBoundedVector334 + crRightHandSideBoundedVector335 + crRightHandSideBoundedVector345 + crRightHandSideBoundedVector493 + crRightHandSideBoundedVector494 + crRightHandSideBoundedVector495 + crRightHandSideBoundedVector499;
const double crRightHandSideBoundedVector614 =             1.0*crRightHandSideBoundedVector385*crRightHandSideBoundedVector404*crRightHandSideBoundedVector76*gamma*h;
const double crRightHandSideBoundedVector615 =             1.0*crRightHandSideBoundedVector116*crRightHandSideBoundedVector406*crRightHandSideBoundedVector423*gamma*h;
const double crRightHandSideBoundedVector616 =             crRightHandSideBoundedVector14*crRightHandSideBoundedVector380;
const double crRightHandSideBoundedVector617 =             0.0277777777777778*r[0];
const double crRightHandSideBoundedVector618 =             0.0277777777777778*r[1];
const double crRightHandSideBoundedVector619 =             0.111111111111111*r[2];
const double crRightHandSideBoundedVector620 =             crRightHandSideBoundedVector14*crRightHandSideBoundedVector18*crRightHandSideBoundedVector44*gamma;
const double crRightHandSideBoundedVector621 =             crRightHandSideBoundedVector14*crRightHandSideBoundedVector154*crRightHandSideBoundedVector23*gamma;
const double crRightHandSideBoundedVector622 =             crRightHandSideBoundedVector18*crRightHandSideBoundedVector23*crRightHandSideBoundedVector27*crRightHandSideBoundedVector276*crRightHandSideBoundedVector56;
const double crRightHandSideBoundedVector623 =             crRightHandSideBoundedVector18*crRightHandSideBoundedVector23*crRightHandSideBoundedVector27*crRightHandSideBoundedVector276*crRightHandSideBoundedVector53;
const double crRightHandSideBoundedVector624 =             crRightHandSideBoundedVector25*crRightHandSideBoundedVector28 + crRightHandSideBoundedVector379;
const double crRightHandSideBoundedVector625 =             crRightHandSideBoundedVector18*crRightHandSideBoundedVector276*crRightHandSideBoundedVector624*crRightHandSideBoundedVector66;
const double crRightHandSideBoundedVector626 =             crRightHandSideBoundedVector23*crRightHandSideBoundedVector276*crRightHandSideBoundedVector624*crRightHandSideBoundedVector63;
const double crRightHandSideBoundedVector627 =             crRightHandSideBoundedVector50*crRightHandSideBoundedVector587 + crRightHandSideBoundedVector581*crRightHandSideBoundedVector61 + crRightHandSideBoundedVector617 + crRightHandSideBoundedVector618 + crRightHandSideBoundedVector619 - 0.166666666666667*crRightHandSideBoundedVector620 - 0.166666666666667*crRightHandSideBoundedVector621 + 0.333333333333333*crRightHandSideBoundedVector622 + 0.333333333333333*crRightHandSideBoundedVector623 - 0.333333333333333*crRightHandSideBoundedVector625 - 0.333333333333333*crRightHandSideBoundedVector626;
const double crRightHandSideBoundedVector628 =             crRightHandSideBoundedVector401*crRightHandSideBoundedVector77;
const double crRightHandSideBoundedVector629 =             0.111111111111111*r[1];
const double crRightHandSideBoundedVector630 =             0.0277777777777778*r[2];
const double crRightHandSideBoundedVector631 =             crRightHandSideBoundedVector44*crRightHandSideBoundedVector77*crRightHandSideBoundedVector80*gamma;
const double crRightHandSideBoundedVector632 =             crRightHandSideBoundedVector154*crRightHandSideBoundedVector77*crRightHandSideBoundedVector84*gamma;
const double crRightHandSideBoundedVector633 =             crRightHandSideBoundedVector27*crRightHandSideBoundedVector293*crRightHandSideBoundedVector56*crRightHandSideBoundedVector80*crRightHandSideBoundedVector84;
const double crRightHandSideBoundedVector634 =             crRightHandSideBoundedVector27*crRightHandSideBoundedVector293*crRightHandSideBoundedVector53*crRightHandSideBoundedVector80*crRightHandSideBoundedVector84;
const double crRightHandSideBoundedVector635 =             crRightHandSideBoundedVector400 + crRightHandSideBoundedVector86*crRightHandSideBoundedVector87;
const double crRightHandSideBoundedVector636 =             crRightHandSideBoundedVector293*crRightHandSideBoundedVector635*crRightHandSideBoundedVector66*crRightHandSideBoundedVector80;
const double crRightHandSideBoundedVector637 =             crRightHandSideBoundedVector293*crRightHandSideBoundedVector63*crRightHandSideBoundedVector635*crRightHandSideBoundedVector84;
const double crRightHandSideBoundedVector638 =             crRightHandSideBoundedVector50*crRightHandSideBoundedVector599 + crRightHandSideBoundedVector593*crRightHandSideBoundedVector61 + crRightHandSideBoundedVector617 + crRightHandSideBoundedVector629 + crRightHandSideBoundedVector630 - 0.166666666666667*crRightHandSideBoundedVector631 - 0.166666666666667*crRightHandSideBoundedVector632 + 0.333333333333333*crRightHandSideBoundedVector633 + 0.333333333333333*crRightHandSideBoundedVector634 - 0.333333333333333*crRightHandSideBoundedVector636 - 0.333333333333333*crRightHandSideBoundedVector637;
const double crRightHandSideBoundedVector639 =             crRightHandSideBoundedVector117*crRightHandSideBoundedVector119*crRightHandSideBoundedVector44*gamma;
const double crRightHandSideBoundedVector640 =             crRightHandSideBoundedVector117*crRightHandSideBoundedVector122*crRightHandSideBoundedVector154*gamma;
const double crRightHandSideBoundedVector641 =             crRightHandSideBoundedVector119*crRightHandSideBoundedVector122*crRightHandSideBoundedVector27*crRightHandSideBoundedVector306*crRightHandSideBoundedVector56;
const double crRightHandSideBoundedVector642 =             crRightHandSideBoundedVector119*crRightHandSideBoundedVector122*crRightHandSideBoundedVector27*crRightHandSideBoundedVector306*crRightHandSideBoundedVector53;
const double crRightHandSideBoundedVector643 =             crRightHandSideBoundedVector124*crRightHandSideBoundedVector125 + crRightHandSideBoundedVector419;
const double crRightHandSideBoundedVector644 =             crRightHandSideBoundedVector119*crRightHandSideBoundedVector306*crRightHandSideBoundedVector643*crRightHandSideBoundedVector66;
const double crRightHandSideBoundedVector645 =             crRightHandSideBoundedVector122*crRightHandSideBoundedVector306*crRightHandSideBoundedVector63*crRightHandSideBoundedVector643;
const double crRightHandSideBoundedVector646 =             crRightHandSideBoundedVector501 + crRightHandSideBoundedVector502 + crRightHandSideBoundedVector503;
const double crRightHandSideBoundedVector647 =             crRightHandSideBoundedVector349 + crRightHandSideBoundedVector350 + crRightHandSideBoundedVector351;
const double crRightHandSideBoundedVector648 =             1.0*crRightHandSideBoundedVector0 + 1.0*crRightHandSideBoundedVector1 + 1.0*crRightHandSideBoundedVector2 + 1.0*crRightHandSideBoundedVector3 + 1.0*crRightHandSideBoundedVector4 + 1.0*crRightHandSideBoundedVector5;
const double crRightHandSideBoundedVector649 =             1.0*DN_DX_1_0*h;
const double crRightHandSideBoundedVector650 =             1.0*DN_DX_1_1*h;
const double crRightHandSideBoundedVector651 =             0.166666666666667*crRightHandSideBoundedVector134 + crRightHandSideBoundedVector135*crRightHandSideBoundedVector202 - 0.166666666666667*crRightHandSideBoundedVector136 + crRightHandSideBoundedVector137*crRightHandSideBoundedVector198 - 0.166666666666667*crRightHandSideBoundedVector138 + 0.166666666666667*crRightHandSideBoundedVector142 - 0.166666666666667*crRightHandSideBoundedVector146 - 1.0*crRightHandSideBoundedVector45;
const double crRightHandSideBoundedVector652 =             crRightHandSideBoundedVector215*crRightHandSideBoundedVector84;
const double crRightHandSideBoundedVector653 =             DN_DX_1_1*crRightHandSideBoundedVector219*mu;
const double crRightHandSideBoundedVector654 =             DN_DX_1_1*crRightHandSideBoundedVector18;
const double crRightHandSideBoundedVector655 =             0.111111111111111*f_ext(0,0);
const double crRightHandSideBoundedVector656 =             crRightHandSideBoundedVector270 + crRightHandSideBoundedVector655 + 0.444444444444444*f_ext(1,0);
const double crRightHandSideBoundedVector657 =             0.666666666666667*crRightHandSideBoundedVector288;
const double crRightHandSideBoundedVector658 =             DN_DX_1_1*crRightHandSideBoundedVector80;
const double crRightHandSideBoundedVector659 =             crRightHandSideBoundedVector210*crRightHandSideBoundedVector285;
const double crRightHandSideBoundedVector660 =             1.33333333333333*crRightHandSideBoundedVector294;
const double crRightHandSideBoundedVector661 =             DN_DX_1_1*crRightHandSideBoundedVector119;
const double crRightHandSideBoundedVector662 =             0.166666666666667*crRightHandSideBoundedVector300;
const double crRightHandSideBoundedVector663 =             crRightHandSideBoundedVector198*crRightHandSideBoundedVector304;
const double crRightHandSideBoundedVector664 =             0.333333333333333*crRightHandSideBoundedVector307;
const double crRightHandSideBoundedVector665 =             crRightHandSideBoundedVector202*crRightHandSideBoundedVector214 + crRightHandSideBoundedVector269 + crRightHandSideBoundedVector287 - 0.166666666666667*crRightHandSideBoundedVector302 - 0.333333333333333*crRightHandSideBoundedVector309 + crRightHandSideBoundedVector655 - crRightHandSideBoundedVector662 + crRightHandSideBoundedVector663 + crRightHandSideBoundedVector664;
const double crRightHandSideBoundedVector666 =             DN_DX_1_0*crRightHandSideBoundedVector23;
const double crRightHandSideBoundedVector667 =             DN_DX_1_0*crRightHandSideBoundedVector84;
const double crRightHandSideBoundedVector668 =             DN_DX_1_0*crRightHandSideBoundedVector122;
const double crRightHandSideBoundedVector669 =             0.166666666666667*crRightHandSideBoundedVector116*crRightHandSideBoundedVector119;
const double crRightHandSideBoundedVector670 =             crRightHandSideBoundedVector63*crRightHandSideBoundedVector669;
const double crRightHandSideBoundedVector671 =             0.166666666666667*crRightHandSideBoundedVector116*crRightHandSideBoundedVector122;
const double crRightHandSideBoundedVector672 =             crRightHandSideBoundedVector66*crRightHandSideBoundedVector671;
const double crRightHandSideBoundedVector673 =             crRightHandSideBoundedVector198 - crRightHandSideBoundedVector27*crRightHandSideBoundedVector672 + crRightHandSideBoundedVector312 + crRightHandSideBoundedVector313 + crRightHandSideBoundedVector314 + crRightHandSideBoundedVector670;
const double crRightHandSideBoundedVector674 =             DN_DX_1_1*crRightHandSideBoundedVector23;
const double crRightHandSideBoundedVector675 =             DN_DX_1_0*crRightHandSideBoundedVector18;
const double crRightHandSideBoundedVector676 =             DN_DX_1_1*crRightHandSideBoundedVector84;
const double crRightHandSideBoundedVector677 =             DN_DX_1_0*crRightHandSideBoundedVector80;
const double crRightHandSideBoundedVector678 =             -crRightHandSideBoundedVector652*crRightHandSideBoundedVector76;
const double crRightHandSideBoundedVector679 =             DN_DX_1_1*crRightHandSideBoundedVector122;
const double crRightHandSideBoundedVector680 =             DN_DX_1_0*crRightHandSideBoundedVector119;
const double crRightHandSideBoundedVector681 =             crRightHandSideBoundedVector63*crRightHandSideBoundedVector671;
const double crRightHandSideBoundedVector682 =             -crRightHandSideBoundedVector681;
const double crRightHandSideBoundedVector683 =             crRightHandSideBoundedVector66*crRightHandSideBoundedVector669;
const double crRightHandSideBoundedVector684 =             crRightHandSideBoundedVector200*crRightHandSideBoundedVector683 + crRightHandSideBoundedVector333 + crRightHandSideBoundedVector334 + crRightHandSideBoundedVector335 + crRightHandSideBoundedVector336 + crRightHandSideBoundedVector682;
const double crRightHandSideBoundedVector685 =             DN_DX_1_0*crRightHandSideBoundedVector310*h;
const double crRightHandSideBoundedVector686 =             crRightHandSideBoundedVector135*crRightHandSideBoundedVector429 + crRightHandSideBoundedVector137*crRightHandSideBoundedVector432 - 1.0*crRightHandSideBoundedVector155 + 0.166666666666667*crRightHandSideBoundedVector184 - 0.166666666666667*crRightHandSideBoundedVector185 - 0.166666666666667*crRightHandSideBoundedVector186 + 0.166666666666667*crRightHandSideBoundedVector188 - 0.166666666666667*crRightHandSideBoundedVector191;
const double crRightHandSideBoundedVector687 =             crRightHandSideBoundedVector218*crRightHandSideBoundedVector80;
const double crRightHandSideBoundedVector688 =             DN_DX_1_0*crRightHandSideBoundedVector219*mu;
const double crRightHandSideBoundedVector689 =             crRightHandSideBoundedVector267*crRightHandSideBoundedVector675;
const double crRightHandSideBoundedVector690 =             0.111111111111111*f_ext(0,1);
const double crRightHandSideBoundedVector691 =             crRightHandSideBoundedVector457 + crRightHandSideBoundedVector690 + 0.444444444444444*f_ext(1,1);
const double crRightHandSideBoundedVector692 =             0.666666666666667*crRightHandSideBoundedVector470;
const double crRightHandSideBoundedVector693 =             crRightHandSideBoundedVector285*crRightHandSideBoundedVector677;
const double crRightHandSideBoundedVector694 =             crRightHandSideBoundedVector292*crRightHandSideBoundedVector440;
const double crRightHandSideBoundedVector695 =             1.33333333333333*crRightHandSideBoundedVector473;
const double crRightHandSideBoundedVector696 =             crRightHandSideBoundedVector304*crRightHandSideBoundedVector680;
const double crRightHandSideBoundedVector697 =             0.166666666666667*crRightHandSideBoundedVector479;
const double crRightHandSideBoundedVector698 =             crRightHandSideBoundedVector214*crRightHandSideBoundedVector429;
const double crRightHandSideBoundedVector699 =             0.333333333333333*crRightHandSideBoundedVector483;
const double crRightHandSideBoundedVector700 =             crRightHandSideBoundedVector304*crRightHandSideBoundedVector432 + crRightHandSideBoundedVector456 + crRightHandSideBoundedVector468 - 0.166666666666667*crRightHandSideBoundedVector478 - 0.333333333333333*crRightHandSideBoundedVector485 + crRightHandSideBoundedVector690 - crRightHandSideBoundedVector697 + crRightHandSideBoundedVector698 + crRightHandSideBoundedVector699;
const double crRightHandSideBoundedVector701 =             crRightHandSideBoundedVector27*crRightHandSideBoundedVector670 + crRightHandSideBoundedVector486 + crRightHandSideBoundedVector487 + crRightHandSideBoundedVector488 + crRightHandSideBoundedVector489 - crRightHandSideBoundedVector672;
const double crRightHandSideBoundedVector702 =             -crRightHandSideBoundedVector687*crRightHandSideBoundedVector76;
const double crRightHandSideBoundedVector703 =             -crRightHandSideBoundedVector683;
const double crRightHandSideBoundedVector704 =             crRightHandSideBoundedVector200*crRightHandSideBoundedVector681 + crRightHandSideBoundedVector493 + crRightHandSideBoundedVector494 + crRightHandSideBoundedVector495 + crRightHandSideBoundedVector496 + crRightHandSideBoundedVector703;
const double crRightHandSideBoundedVector705 =             DN_DX_1_1*crRightHandSideBoundedVector310*h;
const double crRightHandSideBoundedVector706 =             0.166666666666667*crRightHandSideBoundedVector410;
const double crRightHandSideBoundedVector707 =             0.166666666666667*crRightHandSideBoundedVector411;
const double crRightHandSideBoundedVector708 =             crRightHandSideBoundedVector304*crRightHandSideBoundedVector429;
const double crRightHandSideBoundedVector709 =             crRightHandSideBoundedVector119*crRightHandSideBoundedVector663 - crRightHandSideBoundedVector119*crRightHandSideBoundedVector706 + crRightHandSideBoundedVector119*crRightHandSideBoundedVector708 - crRightHandSideBoundedVector122*crRightHandSideBoundedVector707 + 0.166666666666667*crRightHandSideBoundedVector407 + 0.166666666666667*crRightHandSideBoundedVector408 + 0.166666666666667*crRightHandSideBoundedVector409 - 0.166666666666667*crRightHandSideBoundedVector421 - 0.166666666666667*crRightHandSideBoundedVector422 + 0.166666666666667*crRightHandSideBoundedVector545 + 0.166666666666667*crRightHandSideBoundedVector547;
const double crRightHandSideBoundedVector710 =             0.666666666666667*crRightHandSideBoundedVector391;
const double crRightHandSideBoundedVector711 =             0.666666666666667*crRightHandSideBoundedVector392;
const double crRightHandSideBoundedVector712 =             crRightHandSideBoundedVector285*crRightHandSideBoundedVector440;
const double crRightHandSideBoundedVector713 =             crRightHandSideBoundedVector592*crRightHandSideBoundedVector77;
const double crRightHandSideBoundedVector714 =             crRightHandSideBoundedVector598*crRightHandSideBoundedVector77;
const double crRightHandSideBoundedVector715 =             0.166666666666667*crRightHandSideBoundedVector117*crRightHandSideBoundedVector604;
const double crRightHandSideBoundedVector716 =             crRightHandSideBoundedVector269 - crRightHandSideBoundedVector27*crRightHandSideBoundedVector662 + crRightHandSideBoundedVector27*crRightHandSideBoundedVector664 + crRightHandSideBoundedVector287 - crRightHandSideBoundedVector602*crRightHandSideBoundedVector67 + crRightHandSideBoundedVector655 + crRightHandSideBoundedVector66*crRightHandSideBoundedVector715 - crRightHandSideBoundedVector663 + crRightHandSideBoundedVector706 - crRightHandSideBoundedVector708;
const double crRightHandSideBoundedVector717 =             0.166666666666667*crRightHandSideBoundedVector117*crRightHandSideBoundedVector609;
const double crRightHandSideBoundedVector718 =             -crRightHandSideBoundedVector198*crRightHandSideBoundedVector214 - crRightHandSideBoundedVector27*crRightHandSideBoundedVector697 + crRightHandSideBoundedVector27*crRightHandSideBoundedVector699 + crRightHandSideBoundedVector456 + crRightHandSideBoundedVector468 - crRightHandSideBoundedVector607*crRightHandSideBoundedVector67 + crRightHandSideBoundedVector63*crRightHandSideBoundedVector717 + crRightHandSideBoundedVector690 - crRightHandSideBoundedVector698 + crRightHandSideBoundedVector707;
const double crRightHandSideBoundedVector719 =             crRightHandSideBoundedVector333 + crRightHandSideBoundedVector334 + crRightHandSideBoundedVector335 + crRightHandSideBoundedVector493 + crRightHandSideBoundedVector494 + crRightHandSideBoundedVector495 + crRightHandSideBoundedVector682 + crRightHandSideBoundedVector703;
const double crRightHandSideBoundedVector720 =             0.111111111111111*r[0];
const double crRightHandSideBoundedVector721 =             crRightHandSideBoundedVector50*crRightHandSideBoundedVector717 + crRightHandSideBoundedVector61*crRightHandSideBoundedVector715 + crRightHandSideBoundedVector618 + crRightHandSideBoundedVector630 - 0.166666666666667*crRightHandSideBoundedVector639 - 0.166666666666667*crRightHandSideBoundedVector640 + 0.333333333333333*crRightHandSideBoundedVector641 + 0.333333333333333*crRightHandSideBoundedVector642 - 0.333333333333333*crRightHandSideBoundedVector644 - 0.333333333333333*crRightHandSideBoundedVector645 + crRightHandSideBoundedVector720;
const double crRightHandSideBoundedVector722 =             1.0*DN_DX_2_0*h;
const double crRightHandSideBoundedVector723 =             1.0*DN_DX_2_1*h;
const double crRightHandSideBoundedVector724 =             crRightHandSideBoundedVector215*crRightHandSideBoundedVector23;
const double crRightHandSideBoundedVector725 =             DN_DX_2_1*crRightHandSideBoundedVector219*mu;
const double crRightHandSideBoundedVector726 =             crRightHandSideBoundedVector286 + crRightHandSideBoundedVector655 + 0.444444444444444*f_ext(2,0);
const double crRightHandSideBoundedVector727 =             0.666666666666667*crRightHandSideBoundedVector271;
const double crRightHandSideBoundedVector728 =             DN_DX_2_1*crRightHandSideBoundedVector18;
const double crRightHandSideBoundedVector729 =             crRightHandSideBoundedVector210*crRightHandSideBoundedVector267;
const double crRightHandSideBoundedVector730 =             1.33333333333333*crRightHandSideBoundedVector277;
const double crRightHandSideBoundedVector731 =             DN_DX_2_1*crRightHandSideBoundedVector80;
const double crRightHandSideBoundedVector732 =             DN_DX_2_1*crRightHandSideBoundedVector119;
const double crRightHandSideBoundedVector733 =             DN_DX_2_0*crRightHandSideBoundedVector23;
const double crRightHandSideBoundedVector734 =             DN_DX_2_0*crRightHandSideBoundedVector84;
const double crRightHandSideBoundedVector735 =             DN_DX_2_0*crRightHandSideBoundedVector122;
const double crRightHandSideBoundedVector736 =             DN_DX_2_1*crRightHandSideBoundedVector23;
const double crRightHandSideBoundedVector737 =             DN_DX_2_0*crRightHandSideBoundedVector18;
const double crRightHandSideBoundedVector738 =             -crRightHandSideBoundedVector13*crRightHandSideBoundedVector724;
const double crRightHandSideBoundedVector739 =             DN_DX_2_1*crRightHandSideBoundedVector84;
const double crRightHandSideBoundedVector740 =             DN_DX_2_0*crRightHandSideBoundedVector80;
const double crRightHandSideBoundedVector741 =             DN_DX_2_1*crRightHandSideBoundedVector122;
const double crRightHandSideBoundedVector742 =             DN_DX_2_0*crRightHandSideBoundedVector119;
const double crRightHandSideBoundedVector743 =             DN_DX_2_0*crRightHandSideBoundedVector310*h;
const double crRightHandSideBoundedVector744 =             crRightHandSideBoundedVector18*crRightHandSideBoundedVector218;
const double crRightHandSideBoundedVector745 =             DN_DX_2_0*crRightHandSideBoundedVector219*mu;
const double crRightHandSideBoundedVector746 =             crRightHandSideBoundedVector467 + crRightHandSideBoundedVector690 + 0.444444444444444*f_ext(2,1);
const double crRightHandSideBoundedVector747 =             0.666666666666667*crRightHandSideBoundedVector459;
const double crRightHandSideBoundedVector748 =             crRightHandSideBoundedVector267*crRightHandSideBoundedVector737;
const double crRightHandSideBoundedVector749 =             crRightHandSideBoundedVector275*crRightHandSideBoundedVector440;
const double crRightHandSideBoundedVector750 =             1.33333333333333*crRightHandSideBoundedVector462;
const double crRightHandSideBoundedVector751 =             crRightHandSideBoundedVector285*crRightHandSideBoundedVector740;
const double crRightHandSideBoundedVector752 =             crRightHandSideBoundedVector304*crRightHandSideBoundedVector742;
const double crRightHandSideBoundedVector753 =             -crRightHandSideBoundedVector13*crRightHandSideBoundedVector744;
const double crRightHandSideBoundedVector754 =             DN_DX_2_1*crRightHandSideBoundedVector310*h;
const double crRightHandSideBoundedVector755 =             0.666666666666667*crRightHandSideBoundedVector370;
const double crRightHandSideBoundedVector756 =             0.666666666666667*crRightHandSideBoundedVector371;
const double crRightHandSideBoundedVector757 =             crRightHandSideBoundedVector267*crRightHandSideBoundedVector440;
const double crRightHandSideBoundedVector758 =             crRightHandSideBoundedVector14*crRightHandSideBoundedVector580;
const double crRightHandSideBoundedVector759 =             crRightHandSideBoundedVector14*crRightHandSideBoundedVector586;
            rRightHandSideBoundedVector[0]=-1.0*crRightHandSideBoundedVector114*crRightHandSideBoundedVector7 - 1.0*crRightHandSideBoundedVector148*crRightHandSideBoundedVector7 - 1.0*crRightHandSideBoundedVector149*crRightHandSideBoundedVector169 - 1.0*crRightHandSideBoundedVector149*crRightHandSideBoundedVector182 - 1.0*crRightHandSideBoundedVector149*crRightHandSideBoundedVector193 - 1.0*crRightHandSideBoundedVector6 - 1.0*crRightHandSideBoundedVector7*crRightHandSideBoundedVector73;
            rRightHandSideBoundedVector[1]=-DN_DX_0_0*crRightHandSideBoundedVector246 - DN_DX_0_0*crRightHandSideBoundedVector253 - DN_DX_0_0*crRightHandSideBoundedVector259 + 0.666666666666667*crRightHandSideBoundedVector134 + crRightHandSideBoundedVector135*crRightHandSideBoundedVector212 - 0.666666666666667*crRightHandSideBoundedVector136 + crRightHandSideBoundedVector137*crRightHandSideBoundedVector210 - 0.666666666666667*crRightHandSideBoundedVector138 - crRightHandSideBoundedVector145*crRightHandSideBoundedVector218 + crRightHandSideBoundedVector194 + crRightHandSideBoundedVector195 + crRightHandSideBoundedVector196 + crRightHandSideBoundedVector199 + crRightHandSideBoundedVector203 + crRightHandSideBoundedVector204 + crRightHandSideBoundedVector205 + crRightHandSideBoundedVector206 + crRightHandSideBoundedVector207 + crRightHandSideBoundedVector208 + crRightHandSideBoundedVector209 + crRightHandSideBoundedVector211 + crRightHandSideBoundedVector213 + crRightHandSideBoundedVector214*crRightHandSideBoundedVector216 + crRightHandSideBoundedVector217 - crRightHandSideBoundedVector220*crRightHandSideBoundedVector227 - crRightHandSideBoundedVector220*crRightHandSideBoundedVector229 - crRightHandSideBoundedVector220*crRightHandSideBoundedVector231 - crRightHandSideBoundedVector265*(DN_DX_0_0*crRightHandSideBoundedVector70 - crRightHandSideBoundedVector266*crRightHandSideBoundedVector267 + crRightHandSideBoundedVector280) - crRightHandSideBoundedVector283*(DN_DX_0_0*crRightHandSideBoundedVector111 - crRightHandSideBoundedVector284*crRightHandSideBoundedVector285 + crRightHandSideBoundedVector297) - crRightHandSideBoundedVector298*(DN_DX_0_0*crRightHandSideBoundedVector145 + crRightHandSideBoundedVector212*crRightHandSideBoundedVector214 + crRightHandSideBoundedVector299 - crRightHandSideBoundedVector301 - 0.666666666666667*crRightHandSideBoundedVector302 - crRightHandSideBoundedVector303*crRightHandSideBoundedVector304 + crRightHandSideBoundedVector305 + crRightHandSideBoundedVector308 - 1.33333333333333*crRightHandSideBoundedVector309) + crRightHandSideBoundedVector320*(-crRightHandSideBoundedVector266 + crRightHandSideBoundedVector310*crRightHandSideBoundedVector311 + crRightHandSideBoundedVector319) + crRightHandSideBoundedVector327*(-crRightHandSideBoundedVector284 + crRightHandSideBoundedVector310*crRightHandSideBoundedVector321 + crRightHandSideBoundedVector326) + crRightHandSideBoundedVector330*(-crRightHandSideBoundedVector116*crRightHandSideBoundedVector122*crRightHandSideBoundedVector218*crRightHandSideBoundedVector27 + crRightHandSideBoundedVector135*crRightHandSideBoundedVector215 - crRightHandSideBoundedVector303 + crRightHandSideBoundedVector310*crRightHandSideBoundedVector329 + crRightHandSideBoundedVector328) - crRightHandSideBoundedVector341*(-crRightHandSideBoundedVector200*crRightHandSideBoundedVector332 + crRightHandSideBoundedVector331 + crRightHandSideBoundedVector340) - crRightHandSideBoundedVector348*(-crRightHandSideBoundedVector200*crRightHandSideBoundedVector343 + crRightHandSideBoundedVector342 + crRightHandSideBoundedVector347) - crRightHandSideBoundedVector356*(crRightHandSideBoundedVector116*crRightHandSideBoundedVector119*crRightHandSideBoundedVector200*crRightHandSideBoundedVector218 - crRightHandSideBoundedVector200*crRightHandSideBoundedVector354 + crRightHandSideBoundedVector352 + crRightHandSideBoundedVector353 + crRightHandSideBoundedVector355) - crRightHandSideBoundedVector357*crRightHandSideBoundedVector384 - crRightHandSideBoundedVector357*crRightHandSideBoundedVector405 - crRightHandSideBoundedVector357*crRightHandSideBoundedVector424 - 1.0*crRightHandSideBoundedVector45;
            rRightHandSideBoundedVector[2]=-DN_DX_0_1*crRightHandSideBoundedVector451 - DN_DX_0_1*crRightHandSideBoundedVector452 - DN_DX_0_1*crRightHandSideBoundedVector453 + crRightHandSideBoundedVector135*crRightHandSideBoundedVector440 + crRightHandSideBoundedVector137*crRightHandSideBoundedVector442 - 1.0*crRightHandSideBoundedVector155 + 0.666666666666667*crRightHandSideBoundedVector184 - 0.666666666666667*crRightHandSideBoundedVector185 - 0.666666666666667*crRightHandSideBoundedVector186 - crRightHandSideBoundedVector190*crRightHandSideBoundedVector215 - crRightHandSideBoundedVector227*crRightHandSideBoundedVector446 - crRightHandSideBoundedVector229*crRightHandSideBoundedVector446 - crRightHandSideBoundedVector231*crRightHandSideBoundedVector446 - crRightHandSideBoundedVector265*(DN_DX_0_1*crRightHandSideBoundedVector166 - crRightHandSideBoundedVector454 + crRightHandSideBoundedVector465) - crRightHandSideBoundedVector283*(DN_DX_0_1*crRightHandSideBoundedVector179 - crRightHandSideBoundedVector466 + crRightHandSideBoundedVector476) - crRightHandSideBoundedVector298*(DN_DX_0_1*crRightHandSideBoundedVector190 + crRightHandSideBoundedVector304*crRightHandSideBoundedVector442 + crRightHandSideBoundedVector477 - 0.666666666666667*crRightHandSideBoundedVector478 - crRightHandSideBoundedVector480 - crRightHandSideBoundedVector481 + crRightHandSideBoundedVector482 + crRightHandSideBoundedVector484 - 1.33333333333333*crRightHandSideBoundedVector485) + crRightHandSideBoundedVector304*crRightHandSideBoundedVector444 - crRightHandSideBoundedVector320*(-crRightHandSideBoundedVector200*crRightHandSideBoundedVector331 + crRightHandSideBoundedVector332 + crRightHandSideBoundedVector498) - crRightHandSideBoundedVector327*(-crRightHandSideBoundedVector200*crRightHandSideBoundedVector342 + crRightHandSideBoundedVector343 + crRightHandSideBoundedVector500) - crRightHandSideBoundedVector330*(crRightHandSideBoundedVector116*crRightHandSideBoundedVector122*crRightHandSideBoundedVector200*crRightHandSideBoundedVector215 - crRightHandSideBoundedVector200*crRightHandSideBoundedVector353 + crRightHandSideBoundedVector354 + crRightHandSideBoundedVector504 + crRightHandSideBoundedVector505) - crRightHandSideBoundedVector341*(-crRightHandSideBoundedVector266*crRightHandSideBoundedVector310 + crRightHandSideBoundedVector311 + crRightHandSideBoundedVector490) - crRightHandSideBoundedVector348*(-crRightHandSideBoundedVector284*crRightHandSideBoundedVector310 + crRightHandSideBoundedVector321 + crRightHandSideBoundedVector491) - crRightHandSideBoundedVector356*(crRightHandSideBoundedVector116*crRightHandSideBoundedVector119*crRightHandSideBoundedVector215*crRightHandSideBoundedVector27 - crRightHandSideBoundedVector137*crRightHandSideBoundedVector218 - crRightHandSideBoundedVector303*crRightHandSideBoundedVector310 + crRightHandSideBoundedVector329 + crRightHandSideBoundedVector492) - crRightHandSideBoundedVector384*crRightHandSideBoundedVector506 - crRightHandSideBoundedVector405*crRightHandSideBoundedVector506 - crRightHandSideBoundedVector424*crRightHandSideBoundedVector506 + crRightHandSideBoundedVector425 + crRightHandSideBoundedVector426 + crRightHandSideBoundedVector427 + crRightHandSideBoundedVector430 + crRightHandSideBoundedVector433 + crRightHandSideBoundedVector434 + crRightHandSideBoundedVector435 + crRightHandSideBoundedVector436 + crRightHandSideBoundedVector437 + crRightHandSideBoundedVector438 + crRightHandSideBoundedVector439 + crRightHandSideBoundedVector441 + crRightHandSideBoundedVector443 + crRightHandSideBoundedVector445;
            rRightHandSideBoundedVector[3]=-DN_DX_0_0*crRightHandSideBoundedVector560 - DN_DX_0_0*crRightHandSideBoundedVector565 - DN_DX_0_0*crRightHandSideBoundedVector570 - DN_DX_0_1*crRightHandSideBoundedVector573 - DN_DX_0_1*crRightHandSideBoundedVector574 - DN_DX_0_1*crRightHandSideBoundedVector575 + crRightHandSideBoundedVector119*crRightHandSideBoundedVector305 - crRightHandSideBoundedVector119*crRightHandSideBoundedVector532 + crRightHandSideBoundedVector119*crRightHandSideBoundedVector537 - crRightHandSideBoundedVector122*crRightHandSideBoundedVector533 - crRightHandSideBoundedVector216*crRightHandSideBoundedVector550 - crRightHandSideBoundedVector265*(crRightHandSideBoundedVector331*crRightHandSideBoundedVector616 + crRightHandSideBoundedVector332*crRightHandSideBoundedVector616 + crRightHandSideBoundedVector627) - crRightHandSideBoundedVector283*(crRightHandSideBoundedVector342*crRightHandSideBoundedVector628 + crRightHandSideBoundedVector343*crRightHandSideBoundedVector628 + crRightHandSideBoundedVector638) - crRightHandSideBoundedVector298*(crRightHandSideBoundedVector353*crRightHandSideBoundedVector550 + crRightHandSideBoundedVector354*crRightHandSideBoundedVector550 + crRightHandSideBoundedVector605*crRightHandSideBoundedVector646 + crRightHandSideBoundedVector610*crRightHandSideBoundedVector647 + crRightHandSideBoundedVector619 + crRightHandSideBoundedVector629 - 0.666666666666667*crRightHandSideBoundedVector639 - 0.666666666666667*crRightHandSideBoundedVector640 + 1.33333333333333*crRightHandSideBoundedVector641 + 1.33333333333333*crRightHandSideBoundedVector642 - 1.33333333333333*crRightHandSideBoundedVector644 - 1.33333333333333*crRightHandSideBoundedVector645 + 0.444444444444444*r[0]) + 0.666666666666667*crRightHandSideBoundedVector407 + 0.666666666666667*crRightHandSideBoundedVector408 + 0.666666666666667*crRightHandSideBoundedVector409 - crRightHandSideBoundedVector444*crRightHandSideBoundedVector550 + crRightHandSideBoundedVector507 + crRightHandSideBoundedVector508 + crRightHandSideBoundedVector509 + crRightHandSideBoundedVector511 + crRightHandSideBoundedVector513 + crRightHandSideBoundedVector514 + crRightHandSideBoundedVector516 + crRightHandSideBoundedVector519 + crRightHandSideBoundedVector522 + crRightHandSideBoundedVector523 + crRightHandSideBoundedVector524 + crRightHandSideBoundedVector525 + crRightHandSideBoundedVector526 + crRightHandSideBoundedVector527 + crRightHandSideBoundedVector529 + crRightHandSideBoundedVector531 + crRightHandSideBoundedVector534 + crRightHandSideBoundedVector536 + crRightHandSideBoundedVector540 + crRightHandSideBoundedVector543 + 0.666666666666667*crRightHandSideBoundedVector545 + 0.666666666666667*crRightHandSideBoundedVector547 + crRightHandSideBoundedVector548 + crRightHandSideBoundedVector549 - crRightHandSideBoundedVector576*(-DN_DX_0_0*crRightHandSideBoundedVector577 - DN_DX_0_1*crRightHandSideBoundedVector18*crRightHandSideBoundedVector267*crRightHandSideBoundedVector310 + crRightHandSideBoundedVector582) - crRightHandSideBoundedVector583*(-DN_DX_0_1*crRightHandSideBoundedVector584 - crRightHandSideBoundedVector310*crRightHandSideBoundedVector454 + crRightHandSideBoundedVector588) - crRightHandSideBoundedVector589*(-DN_DX_0_0*crRightHandSideBoundedVector590 - DN_DX_0_1*crRightHandSideBoundedVector285*crRightHandSideBoundedVector310*crRightHandSideBoundedVector80 + crRightHandSideBoundedVector594) - crRightHandSideBoundedVector595*(-DN_DX_0_1*crRightHandSideBoundedVector596 - crRightHandSideBoundedVector310*crRightHandSideBoundedVector466 + crRightHandSideBoundedVector600) - crRightHandSideBoundedVector601*(-DN_DX_0_0*crRightHandSideBoundedVector603 - DN_DX_0_1*crRightHandSideBoundedVector119*crRightHandSideBoundedVector304*crRightHandSideBoundedVector310 + crRightHandSideBoundedVector218*crRightHandSideBoundedVector605 - crRightHandSideBoundedVector27*crRightHandSideBoundedVector301 + crRightHandSideBoundedVector27*crRightHandSideBoundedVector308 + crRightHandSideBoundedVector299 - crRightHandSideBoundedVector305 + crRightHandSideBoundedVector532 - crRightHandSideBoundedVector537 - crRightHandSideBoundedVector579*crRightHandSideBoundedVector602) - crRightHandSideBoundedVector606*(-DN_DX_0_1*crRightHandSideBoundedVector608 - crRightHandSideBoundedVector210*crRightHandSideBoundedVector214 + crRightHandSideBoundedVector215*crRightHandSideBoundedVector610 - crRightHandSideBoundedVector27*crRightHandSideBoundedVector480 + crRightHandSideBoundedVector27*crRightHandSideBoundedVector484 - crRightHandSideBoundedVector310*crRightHandSideBoundedVector481 + crRightHandSideBoundedVector477 - crRightHandSideBoundedVector482 + crRightHandSideBoundedVector533 - crRightHandSideBoundedVector579*crRightHandSideBoundedVector607) - crRightHandSideBoundedVector612*(crRightHandSideBoundedVector331 + crRightHandSideBoundedVector332 + crRightHandSideBoundedVector611) - crRightHandSideBoundedVector614*(crRightHandSideBoundedVector342 + crRightHandSideBoundedVector343 + crRightHandSideBoundedVector613) - crRightHandSideBoundedVector615*(crRightHandSideBoundedVector349 + crRightHandSideBoundedVector350 + crRightHandSideBoundedVector351 + crRightHandSideBoundedVector353 + crRightHandSideBoundedVector354 + crRightHandSideBoundedVector355 + crRightHandSideBoundedVector501 + crRightHandSideBoundedVector502 + crRightHandSideBoundedVector503 + crRightHandSideBoundedVector505);
            rRightHandSideBoundedVector[4]=-crRightHandSideBoundedVector114*crRightHandSideBoundedVector649 - crRightHandSideBoundedVector148*crRightHandSideBoundedVector649 - crRightHandSideBoundedVector169*crRightHandSideBoundedVector650 - crRightHandSideBoundedVector182*crRightHandSideBoundedVector650 - crRightHandSideBoundedVector193*crRightHandSideBoundedVector650 - crRightHandSideBoundedVector648 - crRightHandSideBoundedVector649*crRightHandSideBoundedVector73;
            rRightHandSideBoundedVector[5]=-DN_DX_1_0*crRightHandSideBoundedVector246 - DN_DX_1_0*crRightHandSideBoundedVector253 - DN_DX_1_0*crRightHandSideBoundedVector259 + 0.666666666666667*crRightHandSideBoundedVector100 + crRightHandSideBoundedVector101*crRightHandSideBoundedVector212 - 0.666666666666667*crRightHandSideBoundedVector102 + crRightHandSideBoundedVector103*crRightHandSideBoundedVector210 - 0.666666666666667*crRightHandSideBoundedVector104 - crRightHandSideBoundedVector111*crRightHandSideBoundedVector218 + crRightHandSideBoundedVector206 + crRightHandSideBoundedVector207 + crRightHandSideBoundedVector208 + crRightHandSideBoundedVector209 + crRightHandSideBoundedVector211 + crRightHandSideBoundedVector213 + crRightHandSideBoundedVector217 - crRightHandSideBoundedVector227*crRightHandSideBoundedVector653 - crRightHandSideBoundedVector229*crRightHandSideBoundedVector653 - crRightHandSideBoundedVector231*crRightHandSideBoundedVector653 - crRightHandSideBoundedVector265*(DN_DX_1_0*crRightHandSideBoundedVector70 - crRightHandSideBoundedVector267*crRightHandSideBoundedVector654 + crRightHandSideBoundedVector280) - crRightHandSideBoundedVector283*(DN_DX_1_0*crRightHandSideBoundedVector111 + crRightHandSideBoundedVector212*crRightHandSideBoundedVector292 - crRightHandSideBoundedVector285*crRightHandSideBoundedVector658 - 0.666666666666667*crRightHandSideBoundedVector290 - 1.33333333333333*crRightHandSideBoundedVector296 + crRightHandSideBoundedVector656 - crRightHandSideBoundedVector657 + crRightHandSideBoundedVector659 + crRightHandSideBoundedVector660) + crRightHandSideBoundedVector292*crRightHandSideBoundedVector652 - crRightHandSideBoundedVector298*(DN_DX_1_0*crRightHandSideBoundedVector145 - crRightHandSideBoundedVector304*crRightHandSideBoundedVector661 + crRightHandSideBoundedVector665) + crRightHandSideBoundedVector320*(crRightHandSideBoundedVector310*crRightHandSideBoundedVector666 + crRightHandSideBoundedVector319 - crRightHandSideBoundedVector654) + crRightHandSideBoundedVector327*(crRightHandSideBoundedVector101*crRightHandSideBoundedVector215 - crRightHandSideBoundedVector218*crRightHandSideBoundedVector27*crRightHandSideBoundedVector76*crRightHandSideBoundedVector84 + crRightHandSideBoundedVector310*crRightHandSideBoundedVector667 + crRightHandSideBoundedVector328 - crRightHandSideBoundedVector658) + crRightHandSideBoundedVector330*(crRightHandSideBoundedVector310*crRightHandSideBoundedVector668 - crRightHandSideBoundedVector661 + crRightHandSideBoundedVector673) - crRightHandSideBoundedVector341*(-crRightHandSideBoundedVector200*crRightHandSideBoundedVector675 + crRightHandSideBoundedVector340 + crRightHandSideBoundedVector674) - crRightHandSideBoundedVector348*(crRightHandSideBoundedVector200*crRightHandSideBoundedVector218*crRightHandSideBoundedVector76*crRightHandSideBoundedVector80 - crRightHandSideBoundedVector200*crRightHandSideBoundedVector677 + crRightHandSideBoundedVector352 + crRightHandSideBoundedVector676 + crRightHandSideBoundedVector678) - crRightHandSideBoundedVector356*(-crRightHandSideBoundedVector200*crRightHandSideBoundedVector680 + crRightHandSideBoundedVector679 + crRightHandSideBoundedVector684) - crRightHandSideBoundedVector384*crRightHandSideBoundedVector685 - crRightHandSideBoundedVector405*crRightHandSideBoundedVector685 - crRightHandSideBoundedVector424*crRightHandSideBoundedVector685 + crRightHandSideBoundedVector651;
            rRightHandSideBoundedVector[6]=-DN_DX_1_1*crRightHandSideBoundedVector451 - DN_DX_1_1*crRightHandSideBoundedVector452 - DN_DX_1_1*crRightHandSideBoundedVector453 + crRightHandSideBoundedVector101*crRightHandSideBoundedVector440 + crRightHandSideBoundedVector103*crRightHandSideBoundedVector442 + 0.666666666666667*crRightHandSideBoundedVector173 - 0.666666666666667*crRightHandSideBoundedVector174 - 0.666666666666667*crRightHandSideBoundedVector175 - crRightHandSideBoundedVector179*crRightHandSideBoundedVector215 - crRightHandSideBoundedVector227*crRightHandSideBoundedVector688 - crRightHandSideBoundedVector229*crRightHandSideBoundedVector688 - crRightHandSideBoundedVector231*crRightHandSideBoundedVector688 - crRightHandSideBoundedVector265*(DN_DX_1_1*crRightHandSideBoundedVector166 + crRightHandSideBoundedVector465 - crRightHandSideBoundedVector689) - crRightHandSideBoundedVector283*(DN_DX_1_1*crRightHandSideBoundedVector179 + crRightHandSideBoundedVector285*crRightHandSideBoundedVector442 - 0.666666666666667*crRightHandSideBoundedVector469 - 1.33333333333333*crRightHandSideBoundedVector475 + crRightHandSideBoundedVector691 - crRightHandSideBoundedVector692 - crRightHandSideBoundedVector693 + crRightHandSideBoundedVector694 + crRightHandSideBoundedVector695) + crRightHandSideBoundedVector285*crRightHandSideBoundedVector687 - crRightHandSideBoundedVector298*(DN_DX_1_1*crRightHandSideBoundedVector190 - crRightHandSideBoundedVector696 + crRightHandSideBoundedVector700) - crRightHandSideBoundedVector320*(-crRightHandSideBoundedVector200*crRightHandSideBoundedVector674 + crRightHandSideBoundedVector498 + crRightHandSideBoundedVector675) - crRightHandSideBoundedVector327*(crRightHandSideBoundedVector200*crRightHandSideBoundedVector215*crRightHandSideBoundedVector76*crRightHandSideBoundedVector84 - crRightHandSideBoundedVector200*crRightHandSideBoundedVector676 + crRightHandSideBoundedVector504 + crRightHandSideBoundedVector677 + crRightHandSideBoundedVector702) - crRightHandSideBoundedVector330*(-crRightHandSideBoundedVector200*crRightHandSideBoundedVector679 + crRightHandSideBoundedVector680 + crRightHandSideBoundedVector704) - crRightHandSideBoundedVector341*(-crRightHandSideBoundedVector310*crRightHandSideBoundedVector654 + crRightHandSideBoundedVector490 + crRightHandSideBoundedVector666) - crRightHandSideBoundedVector348*(-crRightHandSideBoundedVector103*crRightHandSideBoundedVector218 + crRightHandSideBoundedVector215*crRightHandSideBoundedVector27*crRightHandSideBoundedVector76*crRightHandSideBoundedVector80 - crRightHandSideBoundedVector310*crRightHandSideBoundedVector658 + crRightHandSideBoundedVector492 + crRightHandSideBoundedVector667) - crRightHandSideBoundedVector356*(-crRightHandSideBoundedVector310*crRightHandSideBoundedVector661 + crRightHandSideBoundedVector668 + crRightHandSideBoundedVector701) - crRightHandSideBoundedVector384*crRightHandSideBoundedVector705 - crRightHandSideBoundedVector405*crRightHandSideBoundedVector705 - crRightHandSideBoundedVector424*crRightHandSideBoundedVector705 + crRightHandSideBoundedVector436 + crRightHandSideBoundedVector437 + crRightHandSideBoundedVector438 + crRightHandSideBoundedVector439 + crRightHandSideBoundedVector441 + crRightHandSideBoundedVector443 + crRightHandSideBoundedVector445 + crRightHandSideBoundedVector686;
            rRightHandSideBoundedVector[7]=-DN_DX_1_0*crRightHandSideBoundedVector560 - DN_DX_1_0*crRightHandSideBoundedVector565 - DN_DX_1_0*crRightHandSideBoundedVector570 - DN_DX_1_1*crRightHandSideBoundedVector573 - DN_DX_1_1*crRightHandSideBoundedVector574 - DN_DX_1_1*crRightHandSideBoundedVector575 - crRightHandSideBoundedVector265*(crRightHandSideBoundedVector616*crRightHandSideBoundedVector674 + crRightHandSideBoundedVector616*crRightHandSideBoundedVector675 + crRightHandSideBoundedVector627) - crRightHandSideBoundedVector283*(crRightHandSideBoundedVector619 + crRightHandSideBoundedVector628*crRightHandSideBoundedVector676 + crRightHandSideBoundedVector628*crRightHandSideBoundedVector677 - 0.666666666666667*crRightHandSideBoundedVector631 - 0.666666666666667*crRightHandSideBoundedVector632 + 1.33333333333333*crRightHandSideBoundedVector633 + 1.33333333333333*crRightHandSideBoundedVector634 - 1.33333333333333*crRightHandSideBoundedVector636 - 1.33333333333333*crRightHandSideBoundedVector637 + crRightHandSideBoundedVector646*crRightHandSideBoundedVector713 + crRightHandSideBoundedVector647*crRightHandSideBoundedVector714 + crRightHandSideBoundedVector720 + 0.444444444444444*r[1]) - crRightHandSideBoundedVector298*(crRightHandSideBoundedVector550*crRightHandSideBoundedVector679 + crRightHandSideBoundedVector550*crRightHandSideBoundedVector680 + crRightHandSideBoundedVector721) + 0.666666666666667*crRightHandSideBoundedVector388 + 0.666666666666667*crRightHandSideBoundedVector389 + 0.666666666666667*crRightHandSideBoundedVector390 + 0.666666666666667*crRightHandSideBoundedVector518 + 0.666666666666667*crRightHandSideBoundedVector521 + crRightHandSideBoundedVector525 + crRightHandSideBoundedVector526 + crRightHandSideBoundedVector527 + crRightHandSideBoundedVector529 + crRightHandSideBoundedVector531 + crRightHandSideBoundedVector534 + crRightHandSideBoundedVector536 + crRightHandSideBoundedVector540 + crRightHandSideBoundedVector543 + crRightHandSideBoundedVector548 + crRightHandSideBoundedVector549 - crRightHandSideBoundedVector576*(-DN_DX_1_0*crRightHandSideBoundedVector577 - DN_DX_1_1*crRightHandSideBoundedVector18*crRightHandSideBoundedVector267*crRightHandSideBoundedVector310 + crRightHandSideBoundedVector582) - crRightHandSideBoundedVector583*(-DN_DX_1_1*crRightHandSideBoundedVector584 - crRightHandSideBoundedVector310*crRightHandSideBoundedVector689 + crRightHandSideBoundedVector588) - crRightHandSideBoundedVector589*(-DN_DX_1_0*crRightHandSideBoundedVector590 - DN_DX_1_1*crRightHandSideBoundedVector285*crRightHandSideBoundedVector310*crRightHandSideBoundedVector80 + crRightHandSideBoundedVector218*crRightHandSideBoundedVector713 - crRightHandSideBoundedVector27*crRightHandSideBoundedVector657 + crRightHandSideBoundedVector27*crRightHandSideBoundedVector660 - crRightHandSideBoundedVector579*crRightHandSideBoundedVector591 + crRightHandSideBoundedVector656 - crRightHandSideBoundedVector659 + crRightHandSideBoundedVector710 - crRightHandSideBoundedVector712) - crRightHandSideBoundedVector595*(-DN_DX_1_1*crRightHandSideBoundedVector596 - crRightHandSideBoundedVector210*crRightHandSideBoundedVector292 + crRightHandSideBoundedVector215*crRightHandSideBoundedVector714 - crRightHandSideBoundedVector27*crRightHandSideBoundedVector692 + crRightHandSideBoundedVector27*crRightHandSideBoundedVector695 - crRightHandSideBoundedVector310*crRightHandSideBoundedVector693 - crRightHandSideBoundedVector579*crRightHandSideBoundedVector597 + crRightHandSideBoundedVector691 - crRightHandSideBoundedVector694 + crRightHandSideBoundedVector711) - crRightHandSideBoundedVector601*(-DN_DX_1_0*crRightHandSideBoundedVector603 - DN_DX_1_1*crRightHandSideBoundedVector119*crRightHandSideBoundedVector304*crRightHandSideBoundedVector310 + crRightHandSideBoundedVector716) - crRightHandSideBoundedVector606*(-DN_DX_1_1*crRightHandSideBoundedVector608 - crRightHandSideBoundedVector310*crRightHandSideBoundedVector696 + crRightHandSideBoundedVector718) - crRightHandSideBoundedVector612*(crRightHandSideBoundedVector611 + crRightHandSideBoundedVector674 + crRightHandSideBoundedVector675) - crRightHandSideBoundedVector614*(crRightHandSideBoundedVector349 + crRightHandSideBoundedVector350 + crRightHandSideBoundedVector351 + crRightHandSideBoundedVector501 + crRightHandSideBoundedVector502 + crRightHandSideBoundedVector503 + crRightHandSideBoundedVector676 + crRightHandSideBoundedVector677 + crRightHandSideBoundedVector678 + crRightHandSideBoundedVector702) - crRightHandSideBoundedVector615*(crRightHandSideBoundedVector679 + crRightHandSideBoundedVector680 + crRightHandSideBoundedVector719) - crRightHandSideBoundedVector628*crRightHandSideBoundedVector652 - crRightHandSideBoundedVector628*crRightHandSideBoundedVector687 + crRightHandSideBoundedVector659*crRightHandSideBoundedVector80 + crRightHandSideBoundedVector709 - crRightHandSideBoundedVector710*crRightHandSideBoundedVector80 - crRightHandSideBoundedVector711*crRightHandSideBoundedVector84 + crRightHandSideBoundedVector712*crRightHandSideBoundedVector80;
            rRightHandSideBoundedVector[8]=-crRightHandSideBoundedVector114*crRightHandSideBoundedVector722 - crRightHandSideBoundedVector148*crRightHandSideBoundedVector722 - crRightHandSideBoundedVector169*crRightHandSideBoundedVector723 - crRightHandSideBoundedVector182*crRightHandSideBoundedVector723 - crRightHandSideBoundedVector193*crRightHandSideBoundedVector723 - crRightHandSideBoundedVector648 - crRightHandSideBoundedVector722*crRightHandSideBoundedVector73;
            rRightHandSideBoundedVector[9]=-DN_DX_2_0*crRightHandSideBoundedVector246 - DN_DX_2_0*crRightHandSideBoundedVector253 - DN_DX_2_0*crRightHandSideBoundedVector259 + crRightHandSideBoundedVector194 + crRightHandSideBoundedVector195 + crRightHandSideBoundedVector196 + crRightHandSideBoundedVector199 + crRightHandSideBoundedVector203 + crRightHandSideBoundedVector204 + crRightHandSideBoundedVector205 + crRightHandSideBoundedVector210*crRightHandSideBoundedVector54 + crRightHandSideBoundedVector212*crRightHandSideBoundedVector51 - crRightHandSideBoundedVector218*crRightHandSideBoundedVector70 - crRightHandSideBoundedVector227*crRightHandSideBoundedVector725 - crRightHandSideBoundedVector229*crRightHandSideBoundedVector725 - crRightHandSideBoundedVector231*crRightHandSideBoundedVector725 - crRightHandSideBoundedVector265*(DN_DX_2_0*crRightHandSideBoundedVector70 + crRightHandSideBoundedVector212*crRightHandSideBoundedVector275 - crRightHandSideBoundedVector267*crRightHandSideBoundedVector728 - 0.666666666666667*crRightHandSideBoundedVector273 - 1.33333333333333*crRightHandSideBoundedVector279 + crRightHandSideBoundedVector726 - crRightHandSideBoundedVector727 + crRightHandSideBoundedVector729 + crRightHandSideBoundedVector730) + crRightHandSideBoundedVector275*crRightHandSideBoundedVector724 - crRightHandSideBoundedVector283*(DN_DX_2_0*crRightHandSideBoundedVector111 - crRightHandSideBoundedVector285*crRightHandSideBoundedVector731 + crRightHandSideBoundedVector297) - crRightHandSideBoundedVector298*(DN_DX_2_0*crRightHandSideBoundedVector145 - crRightHandSideBoundedVector304*crRightHandSideBoundedVector732 + crRightHandSideBoundedVector665) + crRightHandSideBoundedVector320*(-crRightHandSideBoundedVector13*crRightHandSideBoundedVector218*crRightHandSideBoundedVector23*crRightHandSideBoundedVector27 + crRightHandSideBoundedVector215*crRightHandSideBoundedVector51 + crRightHandSideBoundedVector310*crRightHandSideBoundedVector733 + crRightHandSideBoundedVector328 - crRightHandSideBoundedVector728) + crRightHandSideBoundedVector327*(crRightHandSideBoundedVector310*crRightHandSideBoundedVector734 + crRightHandSideBoundedVector326 - crRightHandSideBoundedVector731) + crRightHandSideBoundedVector330*(crRightHandSideBoundedVector310*crRightHandSideBoundedVector735 + crRightHandSideBoundedVector673 - crRightHandSideBoundedVector732) - crRightHandSideBoundedVector341*(crRightHandSideBoundedVector13*crRightHandSideBoundedVector18*crRightHandSideBoundedVector200*crRightHandSideBoundedVector218 - crRightHandSideBoundedVector200*crRightHandSideBoundedVector737 + crRightHandSideBoundedVector352 + crRightHandSideBoundedVector736 + crRightHandSideBoundedVector738) - crRightHandSideBoundedVector348*(-crRightHandSideBoundedVector200*crRightHandSideBoundedVector740 + crRightHandSideBoundedVector347 + crRightHandSideBoundedVector739) - crRightHandSideBoundedVector356*(-crRightHandSideBoundedVector200*crRightHandSideBoundedVector742 + crRightHandSideBoundedVector684 + crRightHandSideBoundedVector741) - crRightHandSideBoundedVector384*crRightHandSideBoundedVector743 - crRightHandSideBoundedVector405*crRightHandSideBoundedVector743 - crRightHandSideBoundedVector424*crRightHandSideBoundedVector743 + 0.666666666666667*crRightHandSideBoundedVector49 - 0.666666666666667*crRightHandSideBoundedVector52 - 0.666666666666667*crRightHandSideBoundedVector55 + crRightHandSideBoundedVector651;
            rRightHandSideBoundedVector[10]=-DN_DX_2_1*crRightHandSideBoundedVector451 - DN_DX_2_1*crRightHandSideBoundedVector452 - DN_DX_2_1*crRightHandSideBoundedVector453 + 0.666666666666667*crRightHandSideBoundedVector159 - 0.666666666666667*crRightHandSideBoundedVector160 - 0.666666666666667*crRightHandSideBoundedVector161 - crRightHandSideBoundedVector166*crRightHandSideBoundedVector215 - crRightHandSideBoundedVector227*crRightHandSideBoundedVector745 - crRightHandSideBoundedVector229*crRightHandSideBoundedVector745 - crRightHandSideBoundedVector231*crRightHandSideBoundedVector745 - crRightHandSideBoundedVector265*(DN_DX_2_1*crRightHandSideBoundedVector166 + crRightHandSideBoundedVector267*crRightHandSideBoundedVector442 - 0.666666666666667*crRightHandSideBoundedVector458 - 1.33333333333333*crRightHandSideBoundedVector464 + crRightHandSideBoundedVector746 - crRightHandSideBoundedVector747 - crRightHandSideBoundedVector748 + crRightHandSideBoundedVector749 + crRightHandSideBoundedVector750) + crRightHandSideBoundedVector267*crRightHandSideBoundedVector744 - crRightHandSideBoundedVector283*(DN_DX_2_1*crRightHandSideBoundedVector179 + crRightHandSideBoundedVector476 - crRightHandSideBoundedVector751) - crRightHandSideBoundedVector298*(DN_DX_2_1*crRightHandSideBoundedVector190 + crRightHandSideBoundedVector700 - crRightHandSideBoundedVector752) - crRightHandSideBoundedVector320*(crRightHandSideBoundedVector13*crRightHandSideBoundedVector200*crRightHandSideBoundedVector215*crRightHandSideBoundedVector23 - crRightHandSideBoundedVector200*crRightHandSideBoundedVector736 + crRightHandSideBoundedVector504 + crRightHandSideBoundedVector737 + crRightHandSideBoundedVector753) - crRightHandSideBoundedVector327*(-crRightHandSideBoundedVector200*crRightHandSideBoundedVector739 + crRightHandSideBoundedVector500 + crRightHandSideBoundedVector740) - crRightHandSideBoundedVector330*(-crRightHandSideBoundedVector200*crRightHandSideBoundedVector741 + crRightHandSideBoundedVector704 + crRightHandSideBoundedVector742) - crRightHandSideBoundedVector341*(crRightHandSideBoundedVector13*crRightHandSideBoundedVector18*crRightHandSideBoundedVector215*crRightHandSideBoundedVector27 - crRightHandSideBoundedVector218*crRightHandSideBoundedVector54 - crRightHandSideBoundedVector310*crRightHandSideBoundedVector728 + crRightHandSideBoundedVector492 + crRightHandSideBoundedVector733) - crRightHandSideBoundedVector348*(-crRightHandSideBoundedVector310*crRightHandSideBoundedVector731 + crRightHandSideBoundedVector491 + crRightHandSideBoundedVector734) - crRightHandSideBoundedVector356*(-crRightHandSideBoundedVector310*crRightHandSideBoundedVector732 + crRightHandSideBoundedVector701 + crRightHandSideBoundedVector735) - crRightHandSideBoundedVector384*crRightHandSideBoundedVector754 - crRightHandSideBoundedVector405*crRightHandSideBoundedVector754 - crRightHandSideBoundedVector424*crRightHandSideBoundedVector754 + crRightHandSideBoundedVector425 + crRightHandSideBoundedVector426 + crRightHandSideBoundedVector427 + crRightHandSideBoundedVector430 + crRightHandSideBoundedVector433 + crRightHandSideBoundedVector434 + crRightHandSideBoundedVector435 + crRightHandSideBoundedVector440*crRightHandSideBoundedVector51 + crRightHandSideBoundedVector442*crRightHandSideBoundedVector54 + crRightHandSideBoundedVector686;
            rRightHandSideBoundedVector[11]=-DN_DX_2_0*crRightHandSideBoundedVector560 - DN_DX_2_0*crRightHandSideBoundedVector565 - DN_DX_2_0*crRightHandSideBoundedVector570 - DN_DX_2_1*crRightHandSideBoundedVector573 - DN_DX_2_1*crRightHandSideBoundedVector574 - DN_DX_2_1*crRightHandSideBoundedVector575 + crRightHandSideBoundedVector18*crRightHandSideBoundedVector729 - crRightHandSideBoundedVector18*crRightHandSideBoundedVector755 + crRightHandSideBoundedVector18*crRightHandSideBoundedVector757 - crRightHandSideBoundedVector23*crRightHandSideBoundedVector756 - crRightHandSideBoundedVector265*(crRightHandSideBoundedVector616*crRightHandSideBoundedVector736 + crRightHandSideBoundedVector616*crRightHandSideBoundedVector737 - 0.666666666666667*crRightHandSideBoundedVector620 - 0.666666666666667*crRightHandSideBoundedVector621 + 1.33333333333333*crRightHandSideBoundedVector622 + 1.33333333333333*crRightHandSideBoundedVector623 - 1.33333333333333*crRightHandSideBoundedVector625 - 1.33333333333333*crRightHandSideBoundedVector626 + crRightHandSideBoundedVector629 + crRightHandSideBoundedVector646*crRightHandSideBoundedVector758 + crRightHandSideBoundedVector647*crRightHandSideBoundedVector759 + crRightHandSideBoundedVector720 + 0.444444444444444*r[2]) - crRightHandSideBoundedVector283*(crRightHandSideBoundedVector628*crRightHandSideBoundedVector739 + crRightHandSideBoundedVector628*crRightHandSideBoundedVector740 + crRightHandSideBoundedVector638) - crRightHandSideBoundedVector298*(crRightHandSideBoundedVector550*crRightHandSideBoundedVector741 + crRightHandSideBoundedVector550*crRightHandSideBoundedVector742 + crRightHandSideBoundedVector721) + 0.666666666666667*crRightHandSideBoundedVector367 + 0.666666666666667*crRightHandSideBoundedVector368 + 0.666666666666667*crRightHandSideBoundedVector369 + crRightHandSideBoundedVector507 + crRightHandSideBoundedVector508 + crRightHandSideBoundedVector509 + crRightHandSideBoundedVector511 + crRightHandSideBoundedVector513 + crRightHandSideBoundedVector514 + crRightHandSideBoundedVector516 + crRightHandSideBoundedVector519 + crRightHandSideBoundedVector522 + crRightHandSideBoundedVector523 + crRightHandSideBoundedVector524 + 0.666666666666667*crRightHandSideBoundedVector539 + 0.666666666666667*crRightHandSideBoundedVector542 - crRightHandSideBoundedVector576*(-DN_DX_2_0*crRightHandSideBoundedVector577 - DN_DX_2_1*crRightHandSideBoundedVector18*crRightHandSideBoundedVector267*crRightHandSideBoundedVector310 + crRightHandSideBoundedVector218*crRightHandSideBoundedVector758 - crRightHandSideBoundedVector27*crRightHandSideBoundedVector727 + crRightHandSideBoundedVector27*crRightHandSideBoundedVector730 - crRightHandSideBoundedVector578*crRightHandSideBoundedVector579 + crRightHandSideBoundedVector726 - crRightHandSideBoundedVector729 + crRightHandSideBoundedVector755 - crRightHandSideBoundedVector757) - crRightHandSideBoundedVector583*(-DN_DX_2_1*crRightHandSideBoundedVector584 - crRightHandSideBoundedVector210*crRightHandSideBoundedVector275 + crRightHandSideBoundedVector215*crRightHandSideBoundedVector759 - crRightHandSideBoundedVector27*crRightHandSideBoundedVector747 + crRightHandSideBoundedVector27*crRightHandSideBoundedVector750 - crRightHandSideBoundedVector310*crRightHandSideBoundedVector748 - crRightHandSideBoundedVector579*crRightHandSideBoundedVector585 + crRightHandSideBoundedVector746 - crRightHandSideBoundedVector749 + crRightHandSideBoundedVector756) - crRightHandSideBoundedVector589*(-DN_DX_2_0*crRightHandSideBoundedVector590 - DN_DX_2_1*crRightHandSideBoundedVector285*crRightHandSideBoundedVector310*crRightHandSideBoundedVector80 + crRightHandSideBoundedVector594) - crRightHandSideBoundedVector595*(-DN_DX_2_1*crRightHandSideBoundedVector596 - crRightHandSideBoundedVector310*crRightHandSideBoundedVector751 + crRightHandSideBoundedVector600) - crRightHandSideBoundedVector601*(-DN_DX_2_0*crRightHandSideBoundedVector603 - DN_DX_2_1*crRightHandSideBoundedVector119*crRightHandSideBoundedVector304*crRightHandSideBoundedVector310 + crRightHandSideBoundedVector716) - crRightHandSideBoundedVector606*(-DN_DX_2_1*crRightHandSideBoundedVector608 - crRightHandSideBoundedVector310*crRightHandSideBoundedVector752 + crRightHandSideBoundedVector718) - crRightHandSideBoundedVector612*(crRightHandSideBoundedVector349 + crRightHandSideBoundedVector350 + crRightHandSideBoundedVector351 + crRightHandSideBoundedVector501 + crRightHandSideBoundedVector502 + crRightHandSideBoundedVector503 + crRightHandSideBoundedVector736 + crRightHandSideBoundedVector737 + crRightHandSideBoundedVector738 + crRightHandSideBoundedVector753) - crRightHandSideBoundedVector614*(crRightHandSideBoundedVector613 + crRightHandSideBoundedVector739 + crRightHandSideBoundedVector740) - crRightHandSideBoundedVector615*(crRightHandSideBoundedVector719 + crRightHandSideBoundedVector741 + crRightHandSideBoundedVector742) - crRightHandSideBoundedVector616*crRightHandSideBoundedVector724 - crRightHandSideBoundedVector616*crRightHandSideBoundedVector744 + crRightHandSideBoundedVector709;

    } else {
        const double crRightHandSideBoundedVector0 =             DN_DX_0_0*U_0_1;
const double crRightHandSideBoundedVector1 =             DN_DX_0_1*U_0_2;
const double crRightHandSideBoundedVector2 =             DN_DX_1_0*U_1_1;
const double crRightHandSideBoundedVector3 =             DN_DX_1_1*U_1_2;
const double crRightHandSideBoundedVector4 =             DN_DX_2_0*U_2_1;
const double crRightHandSideBoundedVector5 =             DN_DX_2_1*U_2_2;
const double crRightHandSideBoundedVector6 =             crRightHandSideBoundedVector0 + crRightHandSideBoundedVector1 + crRightHandSideBoundedVector2 + crRightHandSideBoundedVector3 + crRightHandSideBoundedVector4 + crRightHandSideBoundedVector5;
const double crRightHandSideBoundedVector7 =             DN_DX_0_0*h;
const double crRightHandSideBoundedVector8 =             1.0/h;
const double crRightHandSideBoundedVector9 =             1.33333333333333*crRightHandSideBoundedVector8*mu*stab_c1;
const double crRightHandSideBoundedVector10 =             0.166666666666667*U_0_0;
const double crRightHandSideBoundedVector11 =             0.166666666666667*U_1_0;
const double crRightHandSideBoundedVector12 =             0.666666666666667*U_2_0 + crRightHandSideBoundedVector10 + crRightHandSideBoundedVector11;
const double crRightHandSideBoundedVector13 =             1.0/crRightHandSideBoundedVector12;
const double crRightHandSideBoundedVector14 =             pow(crRightHandSideBoundedVector12, -2);
const double crRightHandSideBoundedVector15 =             0.166666666666667*U_0_1;
const double crRightHandSideBoundedVector16 =             0.166666666666667*U_1_1;
const double crRightHandSideBoundedVector17 =             0.666666666666667*U_2_1;
const double crRightHandSideBoundedVector18 =             crRightHandSideBoundedVector15 + crRightHandSideBoundedVector16 + crRightHandSideBoundedVector17;
const double crRightHandSideBoundedVector19 =             pow(crRightHandSideBoundedVector18, 2);
const double crRightHandSideBoundedVector20 =             0.166666666666667*U_0_2;
const double crRightHandSideBoundedVector21 =             0.166666666666667*U_1_2;
const double crRightHandSideBoundedVector22 =             0.666666666666667*U_2_2;
const double crRightHandSideBoundedVector23 =             crRightHandSideBoundedVector20 + crRightHandSideBoundedVector21 + crRightHandSideBoundedVector22;
const double crRightHandSideBoundedVector24 =             pow(crRightHandSideBoundedVector23, 2);
const double crRightHandSideBoundedVector25 =             crRightHandSideBoundedVector19 + crRightHandSideBoundedVector24;
const double crRightHandSideBoundedVector26 =             sqrt(gamma);
const double crRightHandSideBoundedVector27 =             gamma - 1;
const double crRightHandSideBoundedVector28 =             crRightHandSideBoundedVector13*crRightHandSideBoundedVector27;
const double crRightHandSideBoundedVector29 =             0.166666666666667*U_0_3;
const double crRightHandSideBoundedVector30 =             -crRightHandSideBoundedVector29;
const double crRightHandSideBoundedVector31 =             0.166666666666667*U_1_3;
const double crRightHandSideBoundedVector32 =             -crRightHandSideBoundedVector31;
const double crRightHandSideBoundedVector33 =             0.666666666666667*U_2_3;
const double crRightHandSideBoundedVector34 =             -crRightHandSideBoundedVector33;
const double crRightHandSideBoundedVector35 =             0.5*crRightHandSideBoundedVector19;
const double crRightHandSideBoundedVector36 =             0.5*crRightHandSideBoundedVector24;
const double crRightHandSideBoundedVector37 =             crRightHandSideBoundedVector26*sqrt(-crRightHandSideBoundedVector28*(crRightHandSideBoundedVector13*crRightHandSideBoundedVector35 + crRightHandSideBoundedVector13*crRightHandSideBoundedVector36 + crRightHandSideBoundedVector30 + crRightHandSideBoundedVector32 + crRightHandSideBoundedVector34)) + sqrt(crRightHandSideBoundedVector14*crRightHandSideBoundedVector25);
const double crRightHandSideBoundedVector38 =             crRightHandSideBoundedVector37*stab_c2;
const double crRightHandSideBoundedVector39 =             1.0/(crRightHandSideBoundedVector13*crRightHandSideBoundedVector9 + crRightHandSideBoundedVector38);
const double crRightHandSideBoundedVector40 =             DN_DX_0_0*U_0_3 + DN_DX_1_0*U_1_3 + DN_DX_2_0*U_2_3;
const double crRightHandSideBoundedVector41 =             crRightHandSideBoundedVector27*crRightHandSideBoundedVector40;
const double crRightHandSideBoundedVector42 =             crRightHandSideBoundedVector41 + 0.166666666666667*dUdt_0_1;
const double crRightHandSideBoundedVector43 =             0.166666666666667*dUdt_1_1;
const double crRightHandSideBoundedVector44 =             0.166666666666667*f_ext(0,0);
const double crRightHandSideBoundedVector45 =             0.166666666666667*f_ext(1,0);
const double crRightHandSideBoundedVector46 =             crRightHandSideBoundedVector44 + crRightHandSideBoundedVector45 + 0.666666666666667*f_ext(2,0);
const double crRightHandSideBoundedVector47 =             crRightHandSideBoundedVector12*crRightHandSideBoundedVector46;
const double crRightHandSideBoundedVector48 =             crRightHandSideBoundedVector1 + crRightHandSideBoundedVector3 + crRightHandSideBoundedVector5;
const double crRightHandSideBoundedVector49 =             crRightHandSideBoundedVector13*crRightHandSideBoundedVector18;
const double crRightHandSideBoundedVector50 =             crRightHandSideBoundedVector48*crRightHandSideBoundedVector49;
const double crRightHandSideBoundedVector51 =             DN_DX_0_1*U_0_1 + DN_DX_1_1*U_1_1 + DN_DX_2_1*U_2_1;
const double crRightHandSideBoundedVector52 =             crRightHandSideBoundedVector13*crRightHandSideBoundedVector23;
const double crRightHandSideBoundedVector53 =             crRightHandSideBoundedVector51*crRightHandSideBoundedVector52;
const double crRightHandSideBoundedVector54 =             DN_DX_0_0*U_0_2 + DN_DX_1_0*U_1_2 + DN_DX_2_0*U_2_2;
const double crRightHandSideBoundedVector55 =             1.0*crRightHandSideBoundedVector27*crRightHandSideBoundedVector54;
const double crRightHandSideBoundedVector56 =             1.0*gamma;
const double crRightHandSideBoundedVector57 =             -crRightHandSideBoundedVector56 + 3.0;
const double crRightHandSideBoundedVector58 =             crRightHandSideBoundedVector13*crRightHandSideBoundedVector57;
const double crRightHandSideBoundedVector59 =             crRightHandSideBoundedVector0 + crRightHandSideBoundedVector2 + crRightHandSideBoundedVector4;
const double crRightHandSideBoundedVector60 =             crRightHandSideBoundedVector18*crRightHandSideBoundedVector59;
const double crRightHandSideBoundedVector61 =             DN_DX_0_1*U_0_0 + DN_DX_1_1*U_1_0 + DN_DX_2_1*U_2_0;
const double crRightHandSideBoundedVector62 =             crRightHandSideBoundedVector14*crRightHandSideBoundedVector18*crRightHandSideBoundedVector23;
const double crRightHandSideBoundedVector63 =             crRightHandSideBoundedVector61*crRightHandSideBoundedVector62;
const double crRightHandSideBoundedVector64 =             DN_DX_0_0*U_0_0 + DN_DX_1_0*U_1_0 + DN_DX_2_0*U_2_0;
const double crRightHandSideBoundedVector65 =             0.5*gamma - 0.5;
const double crRightHandSideBoundedVector66 =             crRightHandSideBoundedVector25*crRightHandSideBoundedVector65;
const double crRightHandSideBoundedVector67 =             -crRightHandSideBoundedVector19 + crRightHandSideBoundedVector66;
const double crRightHandSideBoundedVector68 =             crRightHandSideBoundedVector14*crRightHandSideBoundedVector67;
const double crRightHandSideBoundedVector69 =             crRightHandSideBoundedVector64*crRightHandSideBoundedVector68;
const double crRightHandSideBoundedVector70 =             crRightHandSideBoundedVector42 + crRightHandSideBoundedVector43 - crRightHandSideBoundedVector47 + crRightHandSideBoundedVector50 - crRightHandSideBoundedVector52*crRightHandSideBoundedVector55 + crRightHandSideBoundedVector53 + crRightHandSideBoundedVector58*crRightHandSideBoundedVector60 - crRightHandSideBoundedVector63 + crRightHandSideBoundedVector69 + 0.666666666666667*dUdt_2_1;
const double crRightHandSideBoundedVector71 =             crRightHandSideBoundedVector39*crRightHandSideBoundedVector70;
const double crRightHandSideBoundedVector72 =             0.166666666666667*U_2_0;
const double crRightHandSideBoundedVector73 =             0.666666666666667*U_1_0 + crRightHandSideBoundedVector10 + crRightHandSideBoundedVector72;
const double crRightHandSideBoundedVector74 =             1.0/crRightHandSideBoundedVector73;
const double crRightHandSideBoundedVector75 =             pow(crRightHandSideBoundedVector73, -2);
const double crRightHandSideBoundedVector76 =             0.666666666666667*U_1_1;
const double crRightHandSideBoundedVector77 =             0.166666666666667*U_2_1;
const double crRightHandSideBoundedVector78 =             crRightHandSideBoundedVector15 + crRightHandSideBoundedVector76 + crRightHandSideBoundedVector77;
const double crRightHandSideBoundedVector79 =             pow(crRightHandSideBoundedVector78, 2);
const double crRightHandSideBoundedVector80 =             0.666666666666667*U_1_2;
const double crRightHandSideBoundedVector81 =             0.166666666666667*U_2_2;
const double crRightHandSideBoundedVector82 =             crRightHandSideBoundedVector20 + crRightHandSideBoundedVector80 + crRightHandSideBoundedVector81;
const double crRightHandSideBoundedVector83 =             pow(crRightHandSideBoundedVector82, 2);
const double crRightHandSideBoundedVector84 =             crRightHandSideBoundedVector79 + crRightHandSideBoundedVector83;
const double crRightHandSideBoundedVector85 =             crRightHandSideBoundedVector27*crRightHandSideBoundedVector74;
const double crRightHandSideBoundedVector86 =             0.666666666666667*U_1_3;
const double crRightHandSideBoundedVector87 =             -crRightHandSideBoundedVector86;
const double crRightHandSideBoundedVector88 =             0.166666666666667*U_2_3;
const double crRightHandSideBoundedVector89 =             -crRightHandSideBoundedVector88;
const double crRightHandSideBoundedVector90 =             0.5*crRightHandSideBoundedVector79;
const double crRightHandSideBoundedVector91 =             0.5*crRightHandSideBoundedVector83;
const double crRightHandSideBoundedVector92 =             crRightHandSideBoundedVector26*sqrt(-crRightHandSideBoundedVector85*(crRightHandSideBoundedVector30 + crRightHandSideBoundedVector74*crRightHandSideBoundedVector90 + crRightHandSideBoundedVector74*crRightHandSideBoundedVector91 + crRightHandSideBoundedVector87 + crRightHandSideBoundedVector89)) + sqrt(crRightHandSideBoundedVector75*crRightHandSideBoundedVector84);
const double crRightHandSideBoundedVector93 =             crRightHandSideBoundedVector92*stab_c2;
const double crRightHandSideBoundedVector94 =             1.0/(crRightHandSideBoundedVector74*crRightHandSideBoundedVector9 + crRightHandSideBoundedVector93);
const double crRightHandSideBoundedVector95 =             0.166666666666667*dUdt_2_1;
const double crRightHandSideBoundedVector96 =             0.166666666666667*f_ext(2,0);
const double crRightHandSideBoundedVector97 =             crRightHandSideBoundedVector44 + crRightHandSideBoundedVector96 + 0.666666666666667*f_ext(1,0);
const double crRightHandSideBoundedVector98 =             crRightHandSideBoundedVector73*crRightHandSideBoundedVector97;
const double crRightHandSideBoundedVector99 =             crRightHandSideBoundedVector74*crRightHandSideBoundedVector78;
const double crRightHandSideBoundedVector100 =             crRightHandSideBoundedVector48*crRightHandSideBoundedVector99;
const double crRightHandSideBoundedVector101 =             crRightHandSideBoundedVector74*crRightHandSideBoundedVector82;
const double crRightHandSideBoundedVector102 =             crRightHandSideBoundedVector101*crRightHandSideBoundedVector51;
const double crRightHandSideBoundedVector103 =             crRightHandSideBoundedVector57*crRightHandSideBoundedVector74;
const double crRightHandSideBoundedVector104 =             crRightHandSideBoundedVector59*crRightHandSideBoundedVector78;
const double crRightHandSideBoundedVector105 =             crRightHandSideBoundedVector75*crRightHandSideBoundedVector78*crRightHandSideBoundedVector82;
const double crRightHandSideBoundedVector106 =             crRightHandSideBoundedVector105*crRightHandSideBoundedVector61;
const double crRightHandSideBoundedVector107 =             crRightHandSideBoundedVector65*crRightHandSideBoundedVector84;
const double crRightHandSideBoundedVector108 =             crRightHandSideBoundedVector107 - crRightHandSideBoundedVector79;
const double crRightHandSideBoundedVector109 =             crRightHandSideBoundedVector108*crRightHandSideBoundedVector75;
const double crRightHandSideBoundedVector110 =             crRightHandSideBoundedVector109*crRightHandSideBoundedVector64;
const double crRightHandSideBoundedVector111 =             crRightHandSideBoundedVector100 - crRightHandSideBoundedVector101*crRightHandSideBoundedVector55 + crRightHandSideBoundedVector102 + crRightHandSideBoundedVector103*crRightHandSideBoundedVector104 - crRightHandSideBoundedVector106 + crRightHandSideBoundedVector110 + crRightHandSideBoundedVector42 + crRightHandSideBoundedVector95 - crRightHandSideBoundedVector98 + 0.666666666666667*dUdt_1_1;
const double crRightHandSideBoundedVector112 =             crRightHandSideBoundedVector111*crRightHandSideBoundedVector94;
const double crRightHandSideBoundedVector113 =             0.666666666666667*U_0_0 + crRightHandSideBoundedVector11 + crRightHandSideBoundedVector72;
const double crRightHandSideBoundedVector114 =             1.0/crRightHandSideBoundedVector113;
const double crRightHandSideBoundedVector115 =             pow(crRightHandSideBoundedVector113, -2);
const double crRightHandSideBoundedVector116 =             0.666666666666667*U_0_1;
const double crRightHandSideBoundedVector117 =             crRightHandSideBoundedVector116 + crRightHandSideBoundedVector16 + crRightHandSideBoundedVector77;
const double crRightHandSideBoundedVector118 =             pow(crRightHandSideBoundedVector117, 2);
const double crRightHandSideBoundedVector119 =             0.666666666666667*U_0_2;
const double crRightHandSideBoundedVector120 =             crRightHandSideBoundedVector119 + crRightHandSideBoundedVector21 + crRightHandSideBoundedVector81;
const double crRightHandSideBoundedVector121 =             pow(crRightHandSideBoundedVector120, 2);
const double crRightHandSideBoundedVector122 =             crRightHandSideBoundedVector118 + crRightHandSideBoundedVector121;
const double crRightHandSideBoundedVector123 =             crRightHandSideBoundedVector114*crRightHandSideBoundedVector27;
const double crRightHandSideBoundedVector124 =             0.666666666666667*U_0_3;
const double crRightHandSideBoundedVector125 =             -crRightHandSideBoundedVector124;
const double crRightHandSideBoundedVector126 =             0.5*crRightHandSideBoundedVector118;
const double crRightHandSideBoundedVector127 =             0.5*crRightHandSideBoundedVector121;
const double crRightHandSideBoundedVector128 =             crRightHandSideBoundedVector26*sqrt(-crRightHandSideBoundedVector123*(crRightHandSideBoundedVector114*crRightHandSideBoundedVector126 + crRightHandSideBoundedVector114*crRightHandSideBoundedVector127 + crRightHandSideBoundedVector125 + crRightHandSideBoundedVector32 + crRightHandSideBoundedVector89)) + sqrt(crRightHandSideBoundedVector115*crRightHandSideBoundedVector122);
const double crRightHandSideBoundedVector129 =             crRightHandSideBoundedVector128*stab_c2;
const double crRightHandSideBoundedVector130 =             1.0/(crRightHandSideBoundedVector114*crRightHandSideBoundedVector9 + crRightHandSideBoundedVector129);
const double crRightHandSideBoundedVector131 =             crRightHandSideBoundedVector45 + crRightHandSideBoundedVector96 + 0.666666666666667*f_ext(0,0);
const double crRightHandSideBoundedVector132 =             crRightHandSideBoundedVector113*crRightHandSideBoundedVector131;
const double crRightHandSideBoundedVector133 =             crRightHandSideBoundedVector114*crRightHandSideBoundedVector117;
const double crRightHandSideBoundedVector134 =             crRightHandSideBoundedVector133*crRightHandSideBoundedVector48;
const double crRightHandSideBoundedVector135 =             crRightHandSideBoundedVector114*crRightHandSideBoundedVector120;
const double crRightHandSideBoundedVector136 =             crRightHandSideBoundedVector135*crRightHandSideBoundedVector51;
const double crRightHandSideBoundedVector137 =             crRightHandSideBoundedVector114*crRightHandSideBoundedVector57;
const double crRightHandSideBoundedVector138 =             crRightHandSideBoundedVector117*crRightHandSideBoundedVector59;
const double crRightHandSideBoundedVector139 =             crRightHandSideBoundedVector115*crRightHandSideBoundedVector117*crRightHandSideBoundedVector120;
const double crRightHandSideBoundedVector140 =             crRightHandSideBoundedVector139*crRightHandSideBoundedVector61;
const double crRightHandSideBoundedVector141 =             crRightHandSideBoundedVector122*crRightHandSideBoundedVector65;
const double crRightHandSideBoundedVector142 =             -crRightHandSideBoundedVector118 + crRightHandSideBoundedVector141;
const double crRightHandSideBoundedVector143 =             crRightHandSideBoundedVector115*crRightHandSideBoundedVector142;
const double crRightHandSideBoundedVector144 =             crRightHandSideBoundedVector143*crRightHandSideBoundedVector64;
const double crRightHandSideBoundedVector145 =             -crRightHandSideBoundedVector132 + crRightHandSideBoundedVector134 - crRightHandSideBoundedVector135*crRightHandSideBoundedVector55 + crRightHandSideBoundedVector136 + crRightHandSideBoundedVector137*crRightHandSideBoundedVector138 - crRightHandSideBoundedVector140 + crRightHandSideBoundedVector144 + crRightHandSideBoundedVector41 + crRightHandSideBoundedVector43 + crRightHandSideBoundedVector95 + 0.666666666666667*dUdt_0_1;
const double crRightHandSideBoundedVector146 =             crRightHandSideBoundedVector130*crRightHandSideBoundedVector145;
const double crRightHandSideBoundedVector147 =             DN_DX_0_1*h;
const double crRightHandSideBoundedVector148 =             DN_DX_0_1*U_0_3 + DN_DX_1_1*U_1_3 + DN_DX_2_1*U_2_3;
const double crRightHandSideBoundedVector149 =             crRightHandSideBoundedVector148*crRightHandSideBoundedVector27;
const double crRightHandSideBoundedVector150 =             crRightHandSideBoundedVector149 + 0.166666666666667*dUdt_0_2;
const double crRightHandSideBoundedVector151 =             0.166666666666667*dUdt_1_2;
const double crRightHandSideBoundedVector152 =             0.166666666666667*f_ext(0,1);
const double crRightHandSideBoundedVector153 =             0.166666666666667*f_ext(1,1);
const double crRightHandSideBoundedVector154 =             crRightHandSideBoundedVector152 + crRightHandSideBoundedVector153 + 0.666666666666667*f_ext(2,1);
const double crRightHandSideBoundedVector155 =             crRightHandSideBoundedVector12*crRightHandSideBoundedVector154;
const double crRightHandSideBoundedVector156 =             crRightHandSideBoundedVector49*crRightHandSideBoundedVector54;
const double crRightHandSideBoundedVector157 =             crRightHandSideBoundedVector52*crRightHandSideBoundedVector59;
const double crRightHandSideBoundedVector158 =             1.0*crRightHandSideBoundedVector27*crRightHandSideBoundedVector51;
const double crRightHandSideBoundedVector159 =             crRightHandSideBoundedVector23*crRightHandSideBoundedVector48;
const double crRightHandSideBoundedVector160 =             crRightHandSideBoundedVector62*crRightHandSideBoundedVector64;
const double crRightHandSideBoundedVector161 =             -crRightHandSideBoundedVector24 + crRightHandSideBoundedVector66;
const double crRightHandSideBoundedVector162 =             crRightHandSideBoundedVector14*crRightHandSideBoundedVector161;
const double crRightHandSideBoundedVector163 =             crRightHandSideBoundedVector162*crRightHandSideBoundedVector61;
const double crRightHandSideBoundedVector164 =             crRightHandSideBoundedVector150 + crRightHandSideBoundedVector151 - crRightHandSideBoundedVector155 + crRightHandSideBoundedVector156 + crRightHandSideBoundedVector157 - crRightHandSideBoundedVector158*crRightHandSideBoundedVector49 + crRightHandSideBoundedVector159*crRightHandSideBoundedVector58 - crRightHandSideBoundedVector160 + crRightHandSideBoundedVector163 + 0.666666666666667*dUdt_2_2;
const double crRightHandSideBoundedVector165 =             crRightHandSideBoundedVector164*crRightHandSideBoundedVector39;
const double crRightHandSideBoundedVector166 =             0.166666666666667*dUdt_2_2;
const double crRightHandSideBoundedVector167 =             0.166666666666667*f_ext(2,1);
const double crRightHandSideBoundedVector168 =             crRightHandSideBoundedVector152 + crRightHandSideBoundedVector167 + 0.666666666666667*f_ext(1,1);
const double crRightHandSideBoundedVector169 =             crRightHandSideBoundedVector168*crRightHandSideBoundedVector73;
const double crRightHandSideBoundedVector170 =             crRightHandSideBoundedVector54*crRightHandSideBoundedVector99;
const double crRightHandSideBoundedVector171 =             crRightHandSideBoundedVector101*crRightHandSideBoundedVector59;
const double crRightHandSideBoundedVector172 =             crRightHandSideBoundedVector48*crRightHandSideBoundedVector82;
const double crRightHandSideBoundedVector173 =             crRightHandSideBoundedVector105*crRightHandSideBoundedVector64;
const double crRightHandSideBoundedVector174 =             crRightHandSideBoundedVector107 - crRightHandSideBoundedVector83;
const double crRightHandSideBoundedVector175 =             crRightHandSideBoundedVector174*crRightHandSideBoundedVector75;
const double crRightHandSideBoundedVector176 =             crRightHandSideBoundedVector175*crRightHandSideBoundedVector61;
const double crRightHandSideBoundedVector177 =             crRightHandSideBoundedVector103*crRightHandSideBoundedVector172 + crRightHandSideBoundedVector150 - crRightHandSideBoundedVector158*crRightHandSideBoundedVector99 + crRightHandSideBoundedVector166 - crRightHandSideBoundedVector169 + crRightHandSideBoundedVector170 + crRightHandSideBoundedVector171 - crRightHandSideBoundedVector173 + crRightHandSideBoundedVector176 + 0.666666666666667*dUdt_1_2;
const double crRightHandSideBoundedVector178 =             crRightHandSideBoundedVector177*crRightHandSideBoundedVector94;
const double crRightHandSideBoundedVector179 =             crRightHandSideBoundedVector153 + crRightHandSideBoundedVector167 + 0.666666666666667*f_ext(0,1);
const double crRightHandSideBoundedVector180 =             crRightHandSideBoundedVector113*crRightHandSideBoundedVector179;
const double crRightHandSideBoundedVector181 =             crRightHandSideBoundedVector133*crRightHandSideBoundedVector54;
const double crRightHandSideBoundedVector182 =             crRightHandSideBoundedVector135*crRightHandSideBoundedVector59;
const double crRightHandSideBoundedVector183 =             crRightHandSideBoundedVector120*crRightHandSideBoundedVector48;
const double crRightHandSideBoundedVector184 =             crRightHandSideBoundedVector139*crRightHandSideBoundedVector64;
const double crRightHandSideBoundedVector185 =             -crRightHandSideBoundedVector121 + crRightHandSideBoundedVector141;
const double crRightHandSideBoundedVector186 =             crRightHandSideBoundedVector115*crRightHandSideBoundedVector185;
const double crRightHandSideBoundedVector187 =             crRightHandSideBoundedVector186*crRightHandSideBoundedVector61;
const double crRightHandSideBoundedVector188 =             -crRightHandSideBoundedVector133*crRightHandSideBoundedVector158 + crRightHandSideBoundedVector137*crRightHandSideBoundedVector183 + crRightHandSideBoundedVector149 + crRightHandSideBoundedVector151 + crRightHandSideBoundedVector166 - crRightHandSideBoundedVector180 + crRightHandSideBoundedVector181 + crRightHandSideBoundedVector182 - crRightHandSideBoundedVector184 + crRightHandSideBoundedVector187 + 0.666666666666667*dUdt_0_2;
const double crRightHandSideBoundedVector189 =             crRightHandSideBoundedVector130*crRightHandSideBoundedVector188;
const double crRightHandSideBoundedVector190 =             0.166666666666667*crRightHandSideBoundedVector98;
const double crRightHandSideBoundedVector191 =             -0.166666666666667*crRightHandSideBoundedVector100;
const double crRightHandSideBoundedVector192 =             -0.166666666666667*crRightHandSideBoundedVector102;
const double crRightHandSideBoundedVector193 =             crRightHandSideBoundedVector27*crRightHandSideBoundedVector54;
const double crRightHandSideBoundedVector194 =             0.166666666666667*crRightHandSideBoundedVector193;
const double crRightHandSideBoundedVector195 =             crRightHandSideBoundedVector101*crRightHandSideBoundedVector194;
const double crRightHandSideBoundedVector196 =             crRightHandSideBoundedVector56 - 3.0;
const double crRightHandSideBoundedVector197 =             crRightHandSideBoundedVector196*crRightHandSideBoundedVector59;
const double crRightHandSideBoundedVector198 =             0.166666666666667*crRightHandSideBoundedVector197;
const double crRightHandSideBoundedVector199 =             crRightHandSideBoundedVector198*crRightHandSideBoundedVector99;
const double crRightHandSideBoundedVector200 =             0.166666666666667*crRightHandSideBoundedVector106;
const double crRightHandSideBoundedVector201 =             -0.166666666666667*crRightHandSideBoundedVector110;
const double crRightHandSideBoundedVector202 =             0.166666666666667*crRightHandSideBoundedVector47;
const double crRightHandSideBoundedVector203 =             -0.166666666666667*crRightHandSideBoundedVector50;
const double crRightHandSideBoundedVector204 =             -0.166666666666667*crRightHandSideBoundedVector53;
const double crRightHandSideBoundedVector205 =             crRightHandSideBoundedVector194*crRightHandSideBoundedVector52;
const double crRightHandSideBoundedVector206 =             0.666666666666667*crRightHandSideBoundedVector193;
const double crRightHandSideBoundedVector207 =             crRightHandSideBoundedVector198*crRightHandSideBoundedVector49;
const double crRightHandSideBoundedVector208 =             0.666666666666667*crRightHandSideBoundedVector197;
const double crRightHandSideBoundedVector209 =             0.166666666666667*crRightHandSideBoundedVector63;
const double crRightHandSideBoundedVector210 =             crRightHandSideBoundedVector115*crRightHandSideBoundedVector117;
const double crRightHandSideBoundedVector211 =             0.666666666666667*DN_DX_0_1*U_0_0 + 0.666666666666667*DN_DX_1_1*U_1_0 + 0.666666666666667*DN_DX_2_1*U_2_0;
const double crRightHandSideBoundedVector212 =             crRightHandSideBoundedVector120*crRightHandSideBoundedVector211;
const double crRightHandSideBoundedVector213 =             -0.166666666666667*crRightHandSideBoundedVector69;
const double crRightHandSideBoundedVector214 =             0.666666666666667*DN_DX_0_0*U_0_0 + 0.666666666666667*DN_DX_1_0*U_1_0 + 0.666666666666667*DN_DX_2_0*U_2_0;
const double crRightHandSideBoundedVector215 =             crRightHandSideBoundedVector51 + crRightHandSideBoundedVector54;
const double crRightHandSideBoundedVector216 =             DN_DX_0_1*crRightHandSideBoundedVector215*mu;
const double crRightHandSideBoundedVector217 =             1.0/lin_m_norm;
const double crRightHandSideBoundedVector218 =             1.0/mu;
const double crRightHandSideBoundedVector219 =             crRightHandSideBoundedVector217*crRightHandSideBoundedVector218*lin_m[0]*lin_m[1]*nu_st;
const double crRightHandSideBoundedVector220 =             crRightHandSideBoundedVector12*crRightHandSideBoundedVector219 + 1;
const double crRightHandSideBoundedVector221 =             crRightHandSideBoundedVector217*lin_m[0]*lin_m[1];
const double crRightHandSideBoundedVector222 =             crRightHandSideBoundedVector218*nu_sc*(crRightHandSideBoundedVector221 - 1);
const double crRightHandSideBoundedVector223 =             -crRightHandSideBoundedVector12*crRightHandSideBoundedVector222 + crRightHandSideBoundedVector220;
const double crRightHandSideBoundedVector224 =             crRightHandSideBoundedVector219*crRightHandSideBoundedVector73 + 1;
const double crRightHandSideBoundedVector225 =             -crRightHandSideBoundedVector222*crRightHandSideBoundedVector73 + crRightHandSideBoundedVector224;
const double crRightHandSideBoundedVector226 =             crRightHandSideBoundedVector113*crRightHandSideBoundedVector219 + 1;
const double crRightHandSideBoundedVector227 =             -crRightHandSideBoundedVector113*crRightHandSideBoundedVector222 + crRightHandSideBoundedVector226;
const double crRightHandSideBoundedVector228 =             -nu_sc + nu_st;
const double crRightHandSideBoundedVector229 =             crRightHandSideBoundedVector12*crRightHandSideBoundedVector217*crRightHandSideBoundedVector228*lin_m[0]*lin_m[1];
const double crRightHandSideBoundedVector230 =             2*crRightHandSideBoundedVector1 + 2*crRightHandSideBoundedVector3 + 2*crRightHandSideBoundedVector5;
const double crRightHandSideBoundedVector231 =             0.111111111111111*U_0_0;
const double crRightHandSideBoundedVector232 =             0.111111111111111*U_1_0;
const double crRightHandSideBoundedVector233 =             0.444444444444444*U_2_0 + crRightHandSideBoundedVector231 + crRightHandSideBoundedVector232;
const double crRightHandSideBoundedVector234 =             -crRightHandSideBoundedVector14*(-crRightHandSideBoundedVector18*crRightHandSideBoundedVector214 - crRightHandSideBoundedVector211*crRightHandSideBoundedVector23 + crRightHandSideBoundedVector233*crRightHandSideBoundedVector48 + crRightHandSideBoundedVector233*crRightHandSideBoundedVector59);
const double crRightHandSideBoundedVector235 =             crRightHandSideBoundedVector230 + crRightHandSideBoundedVector234;
const double crRightHandSideBoundedVector236 =             pow(lin_m[0], 2);
const double crRightHandSideBoundedVector237 =             crRightHandSideBoundedVector217*crRightHandSideBoundedVector218*crRightHandSideBoundedVector236*nu_st;
const double crRightHandSideBoundedVector238 =             -crRightHandSideBoundedVector217*crRightHandSideBoundedVector236 + 1;
const double crRightHandSideBoundedVector239 =             crRightHandSideBoundedVector218*crRightHandSideBoundedVector238*nu_sc;
const double crRightHandSideBoundedVector240 =             2*crRightHandSideBoundedVector0 + 2*crRightHandSideBoundedVector2 + 2*crRightHandSideBoundedVector4;
const double crRightHandSideBoundedVector241 =             crRightHandSideBoundedVector234 + crRightHandSideBoundedVector240;
const double crRightHandSideBoundedVector242 =             crRightHandSideBoundedVector229*crRightHandSideBoundedVector235 + crRightHandSideBoundedVector241*mu*(crRightHandSideBoundedVector12*crRightHandSideBoundedVector237 + crRightHandSideBoundedVector12*crRightHandSideBoundedVector239 + 1);
const double crRightHandSideBoundedVector243 =             crRightHandSideBoundedVector217*crRightHandSideBoundedVector228*crRightHandSideBoundedVector73*lin_m[0]*lin_m[1];
const double crRightHandSideBoundedVector244 =             0.111111111111111*U_2_0;
const double crRightHandSideBoundedVector245 =             0.444444444444444*U_1_0 + crRightHandSideBoundedVector231 + crRightHandSideBoundedVector244;
const double crRightHandSideBoundedVector246 =             -crRightHandSideBoundedVector75*(-crRightHandSideBoundedVector211*crRightHandSideBoundedVector82 - crRightHandSideBoundedVector214*crRightHandSideBoundedVector78 + crRightHandSideBoundedVector245*crRightHandSideBoundedVector48 + crRightHandSideBoundedVector245*crRightHandSideBoundedVector59);
const double crRightHandSideBoundedVector247 =             crRightHandSideBoundedVector230 + crRightHandSideBoundedVector246;
const double crRightHandSideBoundedVector248 =             crRightHandSideBoundedVector240 + crRightHandSideBoundedVector246;
const double crRightHandSideBoundedVector249 =             crRightHandSideBoundedVector243*crRightHandSideBoundedVector247 + crRightHandSideBoundedVector248*mu*(crRightHandSideBoundedVector237*crRightHandSideBoundedVector73 + crRightHandSideBoundedVector239*crRightHandSideBoundedVector73 + 1);
const double crRightHandSideBoundedVector250 =             crRightHandSideBoundedVector113*crRightHandSideBoundedVector217*crRightHandSideBoundedVector228*lin_m[0]*lin_m[1];
const double crRightHandSideBoundedVector251 =             0.444444444444444*U_0_0 + crRightHandSideBoundedVector232 + crRightHandSideBoundedVector244;
const double crRightHandSideBoundedVector252 =             -crRightHandSideBoundedVector115*(-crRightHandSideBoundedVector117*crRightHandSideBoundedVector214 - crRightHandSideBoundedVector120*crRightHandSideBoundedVector211 + crRightHandSideBoundedVector251*crRightHandSideBoundedVector48 + crRightHandSideBoundedVector251*crRightHandSideBoundedVector59);
const double crRightHandSideBoundedVector253 =             crRightHandSideBoundedVector230 + crRightHandSideBoundedVector252;
const double crRightHandSideBoundedVector254 =             crRightHandSideBoundedVector240 + crRightHandSideBoundedVector252;
const double crRightHandSideBoundedVector255 =             crRightHandSideBoundedVector250*crRightHandSideBoundedVector253 + crRightHandSideBoundedVector254*mu*(crRightHandSideBoundedVector113*crRightHandSideBoundedVector237 + crRightHandSideBoundedVector113*crRightHandSideBoundedVector239 + 1);
const double crRightHandSideBoundedVector256 =             1.0/stab_c2;
const double crRightHandSideBoundedVector257 =             0.166666666666667*dUdt_1_0;
const double crRightHandSideBoundedVector258 =             crRightHandSideBoundedVector0 + crRightHandSideBoundedVector1 + crRightHandSideBoundedVector2 + crRightHandSideBoundedVector3 + crRightHandSideBoundedVector4 + crRightHandSideBoundedVector5 + 0.166666666666667*dUdt_0_0;
const double crRightHandSideBoundedVector259 =             1.0*crRightHandSideBoundedVector256*h*(crRightHandSideBoundedVector257 + crRightHandSideBoundedVector258 + 0.666666666666667*dUdt_2_0)/crRightHandSideBoundedVector37;
const double crRightHandSideBoundedVector260 =             DN_DX_0_1*crRightHandSideBoundedVector18;
const double crRightHandSideBoundedVector261 =             crRightHandSideBoundedVector14*crRightHandSideBoundedVector23;
const double crRightHandSideBoundedVector262 =             0.0277777777777778*f_ext(0,0);
const double crRightHandSideBoundedVector263 =             0.0277777777777778*f_ext(1,0);
const double crRightHandSideBoundedVector264 =             0.111111111111111*f_ext(2,0);
const double crRightHandSideBoundedVector265 =             crRightHandSideBoundedVector14*crRightHandSideBoundedVector18*crRightHandSideBoundedVector48;
const double crRightHandSideBoundedVector266 =             0.166666666666667*crRightHandSideBoundedVector265;
const double crRightHandSideBoundedVector267 =             crRightHandSideBoundedVector14*crRightHandSideBoundedVector23*crRightHandSideBoundedVector51;
const double crRightHandSideBoundedVector268 =             crRightHandSideBoundedVector194*crRightHandSideBoundedVector261;
const double crRightHandSideBoundedVector269 =             crRightHandSideBoundedVector14*crRightHandSideBoundedVector18;
const double crRightHandSideBoundedVector270 =             pow(crRightHandSideBoundedVector12, -3);
const double crRightHandSideBoundedVector271 =             crRightHandSideBoundedVector18*crRightHandSideBoundedVector23*crRightHandSideBoundedVector270*crRightHandSideBoundedVector61;
const double crRightHandSideBoundedVector272 =             0.333333333333333*crRightHandSideBoundedVector271;
const double crRightHandSideBoundedVector273 =             crRightHandSideBoundedVector270*crRightHandSideBoundedVector64*crRightHandSideBoundedVector67;
const double crRightHandSideBoundedVector274 =             crRightHandSideBoundedVector198*crRightHandSideBoundedVector269 + crRightHandSideBoundedVector262 + crRightHandSideBoundedVector263 + crRightHandSideBoundedVector264 - crRightHandSideBoundedVector266 - 0.166666666666667*crRightHandSideBoundedVector267 + crRightHandSideBoundedVector268 + crRightHandSideBoundedVector272 - 0.333333333333333*crRightHandSideBoundedVector273;
const double crRightHandSideBoundedVector275 =             0.166666666666667*dUdt_2_0;
const double crRightHandSideBoundedVector276 =             1.0*crRightHandSideBoundedVector256*h*(crRightHandSideBoundedVector258 + crRightHandSideBoundedVector275 + 0.666666666666667*dUdt_1_0)/crRightHandSideBoundedVector92;
const double crRightHandSideBoundedVector277 =             DN_DX_0_1*crRightHandSideBoundedVector78;
const double crRightHandSideBoundedVector278 =             crRightHandSideBoundedVector75*crRightHandSideBoundedVector82;
const double crRightHandSideBoundedVector279 =             0.111111111111111*f_ext(1,0);
const double crRightHandSideBoundedVector280 =             0.0277777777777778*f_ext(2,0);
const double crRightHandSideBoundedVector281 =             crRightHandSideBoundedVector48*crRightHandSideBoundedVector75*crRightHandSideBoundedVector78;
const double crRightHandSideBoundedVector282 =             0.166666666666667*crRightHandSideBoundedVector281;
const double crRightHandSideBoundedVector283 =             crRightHandSideBoundedVector51*crRightHandSideBoundedVector75*crRightHandSideBoundedVector82;
const double crRightHandSideBoundedVector284 =             crRightHandSideBoundedVector194*crRightHandSideBoundedVector278;
const double crRightHandSideBoundedVector285 =             crRightHandSideBoundedVector75*crRightHandSideBoundedVector78;
const double crRightHandSideBoundedVector286 =             pow(crRightHandSideBoundedVector73, -3);
const double crRightHandSideBoundedVector287 =             crRightHandSideBoundedVector286*crRightHandSideBoundedVector61*crRightHandSideBoundedVector78*crRightHandSideBoundedVector82;
const double crRightHandSideBoundedVector288 =             0.333333333333333*crRightHandSideBoundedVector287;
const double crRightHandSideBoundedVector289 =             crRightHandSideBoundedVector108*crRightHandSideBoundedVector286*crRightHandSideBoundedVector64;
const double crRightHandSideBoundedVector290 =             crRightHandSideBoundedVector198*crRightHandSideBoundedVector285 + crRightHandSideBoundedVector262 + crRightHandSideBoundedVector279 + crRightHandSideBoundedVector280 - crRightHandSideBoundedVector282 - 0.166666666666667*crRightHandSideBoundedVector283 + crRightHandSideBoundedVector284 + crRightHandSideBoundedVector288 - 0.333333333333333*crRightHandSideBoundedVector289;
const double crRightHandSideBoundedVector291 =             1.0*crRightHandSideBoundedVector256*h*(crRightHandSideBoundedVector257 + crRightHandSideBoundedVector275 + crRightHandSideBoundedVector6 + 0.666666666666667*dUdt_0_0)/crRightHandSideBoundedVector128;
const double crRightHandSideBoundedVector292 =             crRightHandSideBoundedVector264 + crRightHandSideBoundedVector279 + 0.444444444444444*f_ext(0,0);
const double crRightHandSideBoundedVector293 =             crRightHandSideBoundedVector115*crRightHandSideBoundedVector117*crRightHandSideBoundedVector48;
const double crRightHandSideBoundedVector294 =             0.666666666666667*crRightHandSideBoundedVector293;
const double crRightHandSideBoundedVector295 =             crRightHandSideBoundedVector115*crRightHandSideBoundedVector120*crRightHandSideBoundedVector51;
const double crRightHandSideBoundedVector296 =             DN_DX_0_1*crRightHandSideBoundedVector117;
const double crRightHandSideBoundedVector297 =             crRightHandSideBoundedVector115*crRightHandSideBoundedVector120;
const double crRightHandSideBoundedVector298 =             crRightHandSideBoundedVector206*crRightHandSideBoundedVector297;
const double crRightHandSideBoundedVector299 =             pow(crRightHandSideBoundedVector113, -3);
const double crRightHandSideBoundedVector300 =             crRightHandSideBoundedVector117*crRightHandSideBoundedVector120*crRightHandSideBoundedVector299*crRightHandSideBoundedVector61;
const double crRightHandSideBoundedVector301 =             1.33333333333333*crRightHandSideBoundedVector300;
const double crRightHandSideBoundedVector302 =             crRightHandSideBoundedVector142*crRightHandSideBoundedVector299*crRightHandSideBoundedVector64;
const double crRightHandSideBoundedVector303 =             crRightHandSideBoundedVector56 - 1.0;
const double crRightHandSideBoundedVector304 =             DN_DX_0_0*crRightHandSideBoundedVector23;
const double crRightHandSideBoundedVector305 =             -DN_DX_0_1*crRightHandSideBoundedVector15;
const double crRightHandSideBoundedVector306 =             -DN_DX_1_1*crRightHandSideBoundedVector16;
const double crRightHandSideBoundedVector307 =             -DN_DX_2_1*crRightHandSideBoundedVector77;
const double crRightHandSideBoundedVector308 =             0.166666666666667*crRightHandSideBoundedVector13*crRightHandSideBoundedVector18;
const double crRightHandSideBoundedVector309 =             crRightHandSideBoundedVector308*crRightHandSideBoundedVector61;
const double crRightHandSideBoundedVector310 =             0.166666666666667*crRightHandSideBoundedVector13*crRightHandSideBoundedVector23;
const double crRightHandSideBoundedVector311 =             crRightHandSideBoundedVector310*crRightHandSideBoundedVector64;
const double crRightHandSideBoundedVector312 =             crRightHandSideBoundedVector194 - crRightHandSideBoundedVector27*crRightHandSideBoundedVector311 + crRightHandSideBoundedVector305 + crRightHandSideBoundedVector306 + crRightHandSideBoundedVector307 + crRightHandSideBoundedVector309;
const double crRightHandSideBoundedVector313 =             1.0*crRightHandSideBoundedVector13*crRightHandSideBoundedVector164*crRightHandSideBoundedVector39*h;
const double crRightHandSideBoundedVector314 =             DN_DX_0_0*crRightHandSideBoundedVector82;
const double crRightHandSideBoundedVector315 =             0.166666666666667*crRightHandSideBoundedVector74*crRightHandSideBoundedVector78;
const double crRightHandSideBoundedVector316 =             crRightHandSideBoundedVector315*crRightHandSideBoundedVector61;
const double crRightHandSideBoundedVector317 =             0.166666666666667*crRightHandSideBoundedVector74*crRightHandSideBoundedVector82;
const double crRightHandSideBoundedVector318 =             crRightHandSideBoundedVector317*crRightHandSideBoundedVector64;
const double crRightHandSideBoundedVector319 =             crRightHandSideBoundedVector194 - crRightHandSideBoundedVector27*crRightHandSideBoundedVector318 + crRightHandSideBoundedVector305 + crRightHandSideBoundedVector306 + crRightHandSideBoundedVector307 + crRightHandSideBoundedVector316;
const double crRightHandSideBoundedVector320 =             1.0*crRightHandSideBoundedVector177*crRightHandSideBoundedVector74*crRightHandSideBoundedVector94*h;
const double crRightHandSideBoundedVector321 =             -DN_DX_0_1*crRightHandSideBoundedVector116 - DN_DX_1_1*crRightHandSideBoundedVector76 - DN_DX_2_1*crRightHandSideBoundedVector17 + crRightHandSideBoundedVector206;
const double crRightHandSideBoundedVector322 =             DN_DX_0_0*crRightHandSideBoundedVector120;
const double crRightHandSideBoundedVector323 =             1.0*crRightHandSideBoundedVector114*crRightHandSideBoundedVector130*crRightHandSideBoundedVector188*h;
const double crRightHandSideBoundedVector324 =             DN_DX_0_1*crRightHandSideBoundedVector23;
const double crRightHandSideBoundedVector325 =             DN_DX_0_0*crRightHandSideBoundedVector18;
const double crRightHandSideBoundedVector326 =             DN_DX_0_1*crRightHandSideBoundedVector20;
const double crRightHandSideBoundedVector327 =             DN_DX_1_1*crRightHandSideBoundedVector21;
const double crRightHandSideBoundedVector328 =             DN_DX_2_1*crRightHandSideBoundedVector81;
const double crRightHandSideBoundedVector329 =             -crRightHandSideBoundedVector198;
const double crRightHandSideBoundedVector330 =             crRightHandSideBoundedVector310*crRightHandSideBoundedVector61;
const double crRightHandSideBoundedVector331 =             -crRightHandSideBoundedVector330;
const double crRightHandSideBoundedVector332 =             crRightHandSideBoundedVector308*crRightHandSideBoundedVector64;
const double crRightHandSideBoundedVector333 =             crRightHandSideBoundedVector196*crRightHandSideBoundedVector332 + crRightHandSideBoundedVector326 + crRightHandSideBoundedVector327 + crRightHandSideBoundedVector328 + crRightHandSideBoundedVector329 + crRightHandSideBoundedVector331;
const double crRightHandSideBoundedVector334 =             1.0*crRightHandSideBoundedVector13*crRightHandSideBoundedVector39*crRightHandSideBoundedVector70*h;
const double crRightHandSideBoundedVector335 =             DN_DX_0_1*crRightHandSideBoundedVector82;
const double crRightHandSideBoundedVector336 =             DN_DX_0_0*crRightHandSideBoundedVector78;
const double crRightHandSideBoundedVector337 =             crRightHandSideBoundedVector317*crRightHandSideBoundedVector61;
const double crRightHandSideBoundedVector338 =             -crRightHandSideBoundedVector337;
const double crRightHandSideBoundedVector339 =             crRightHandSideBoundedVector315*crRightHandSideBoundedVector64;
const double crRightHandSideBoundedVector340 =             crRightHandSideBoundedVector196*crRightHandSideBoundedVector339 + crRightHandSideBoundedVector326 + crRightHandSideBoundedVector327 + crRightHandSideBoundedVector328 + crRightHandSideBoundedVector329 + crRightHandSideBoundedVector338;
const double crRightHandSideBoundedVector341 =             1.0*crRightHandSideBoundedVector111*crRightHandSideBoundedVector74*crRightHandSideBoundedVector94*h;
const double crRightHandSideBoundedVector342 =             DN_DX_0_1*crRightHandSideBoundedVector119;
const double crRightHandSideBoundedVector343 =             DN_DX_1_1*crRightHandSideBoundedVector80;
const double crRightHandSideBoundedVector344 =             DN_DX_2_1*crRightHandSideBoundedVector22;
const double crRightHandSideBoundedVector345 =             -crRightHandSideBoundedVector208 + crRightHandSideBoundedVector342 + crRightHandSideBoundedVector343 + crRightHandSideBoundedVector344;
const double crRightHandSideBoundedVector346 =             DN_DX_0_1*crRightHandSideBoundedVector120;
const double crRightHandSideBoundedVector347 =             DN_DX_0_0*crRightHandSideBoundedVector117;
const double crRightHandSideBoundedVector348 =             -crRightHandSideBoundedVector114*crRightHandSideBoundedVector212;
const double crRightHandSideBoundedVector349 =             1.0*crRightHandSideBoundedVector114*crRightHandSideBoundedVector130*crRightHandSideBoundedVector145*h;
const double crRightHandSideBoundedVector350 =             DN_DX_0_0*crRightHandSideBoundedVector303*h;
const double crRightHandSideBoundedVector351 =             1.0/c_v;
const double crRightHandSideBoundedVector352 =             crRightHandSideBoundedVector351*crRightHandSideBoundedVector8*lambda*stab_c1/gamma;
const double crRightHandSideBoundedVector353 =             1.0/(crRightHandSideBoundedVector13*crRightHandSideBoundedVector352 + crRightHandSideBoundedVector38);
const double crRightHandSideBoundedVector354 =             0.166666666666667*dUdt_0_3;
const double crRightHandSideBoundedVector355 =             0.166666666666667*dUdt_1_3;
const double crRightHandSideBoundedVector356 =             0.166666666666667*r[0];
const double crRightHandSideBoundedVector357 =             0.166666666666667*r[1];
const double crRightHandSideBoundedVector358 =             crRightHandSideBoundedVector12*(crRightHandSideBoundedVector356 + crRightHandSideBoundedVector357 + 0.666666666666667*r[2]);
const double crRightHandSideBoundedVector359 =             crRightHandSideBoundedVector18*crRightHandSideBoundedVector46;
const double crRightHandSideBoundedVector360 =             crRightHandSideBoundedVector154*crRightHandSideBoundedVector23;
const double crRightHandSideBoundedVector361 =             crRightHandSideBoundedVector13*crRightHandSideBoundedVector40*gamma;
const double crRightHandSideBoundedVector362 =             crRightHandSideBoundedVector13*crRightHandSideBoundedVector148*gamma;
const double crRightHandSideBoundedVector363 =             crRightHandSideBoundedVector13*crRightHandSideBoundedVector19;
const double crRightHandSideBoundedVector364 =             crRightHandSideBoundedVector303*crRightHandSideBoundedVector363;
const double crRightHandSideBoundedVector365 =             crRightHandSideBoundedVector29 + crRightHandSideBoundedVector31 + crRightHandSideBoundedVector33;
const double crRightHandSideBoundedVector366 =             crRightHandSideBoundedVector27*(-crRightHandSideBoundedVector13*(crRightHandSideBoundedVector35 + crRightHandSideBoundedVector36) + crRightHandSideBoundedVector365);
const double crRightHandSideBoundedVector367 =             crRightHandSideBoundedVector29 + crRightHandSideBoundedVector31 + crRightHandSideBoundedVector33 + crRightHandSideBoundedVector366;
const double crRightHandSideBoundedVector368 =             crRightHandSideBoundedVector13*crRightHandSideBoundedVector24;
const double crRightHandSideBoundedVector369 =             crRightHandSideBoundedVector303*crRightHandSideBoundedVector368;
const double crRightHandSideBoundedVector370 =             crRightHandSideBoundedVector30 + crRightHandSideBoundedVector32 + crRightHandSideBoundedVector34 - crRightHandSideBoundedVector366;
const double crRightHandSideBoundedVector371 =             crRightHandSideBoundedVector13*crRightHandSideBoundedVector25*crRightHandSideBoundedVector65 + crRightHandSideBoundedVector370;
const double crRightHandSideBoundedVector372 =             crRightHandSideBoundedVector269*crRightHandSideBoundedVector371*crRightHandSideBoundedVector64;
const double crRightHandSideBoundedVector373 =             crRightHandSideBoundedVector261*crRightHandSideBoundedVector371*crRightHandSideBoundedVector61;
const double crRightHandSideBoundedVector374 =             crRightHandSideBoundedVector13*crRightHandSideBoundedVector48*(crRightHandSideBoundedVector367 - crRightHandSideBoundedVector369) + crRightHandSideBoundedVector13*crRightHandSideBoundedVector59*(-crRightHandSideBoundedVector364 + crRightHandSideBoundedVector367) - crRightHandSideBoundedVector158*crRightHandSideBoundedVector62 + crRightHandSideBoundedVector18*crRightHandSideBoundedVector361 + crRightHandSideBoundedVector23*crRightHandSideBoundedVector362 + crRightHandSideBoundedVector354 + crRightHandSideBoundedVector355 - crRightHandSideBoundedVector358 - crRightHandSideBoundedVector359 - crRightHandSideBoundedVector360 + crRightHandSideBoundedVector372 + crRightHandSideBoundedVector373 - crRightHandSideBoundedVector55*crRightHandSideBoundedVector62 + 0.666666666666667*dUdt_2_3;
const double crRightHandSideBoundedVector375 =             crRightHandSideBoundedVector353*crRightHandSideBoundedVector374;
const double crRightHandSideBoundedVector376 =             1.0/(crRightHandSideBoundedVector352*crRightHandSideBoundedVector74 + crRightHandSideBoundedVector93);
const double crRightHandSideBoundedVector377 =             0.166666666666667*dUdt_2_3;
const double crRightHandSideBoundedVector378 =             0.166666666666667*r[2];
const double crRightHandSideBoundedVector379 =             crRightHandSideBoundedVector73*(crRightHandSideBoundedVector356 + crRightHandSideBoundedVector378 + 0.666666666666667*r[1]);
const double crRightHandSideBoundedVector380 =             crRightHandSideBoundedVector78*crRightHandSideBoundedVector97;
const double crRightHandSideBoundedVector381 =             crRightHandSideBoundedVector168*crRightHandSideBoundedVector82;
const double crRightHandSideBoundedVector382 =             crRightHandSideBoundedVector40*crRightHandSideBoundedVector74*gamma;
const double crRightHandSideBoundedVector383 =             crRightHandSideBoundedVector148*crRightHandSideBoundedVector74*gamma;
const double crRightHandSideBoundedVector384 =             crRightHandSideBoundedVector74*crRightHandSideBoundedVector79;
const double crRightHandSideBoundedVector385 =             crRightHandSideBoundedVector303*crRightHandSideBoundedVector384;
const double crRightHandSideBoundedVector386 =             crRightHandSideBoundedVector29 + crRightHandSideBoundedVector86 + crRightHandSideBoundedVector88;
const double crRightHandSideBoundedVector387 =             crRightHandSideBoundedVector27*(crRightHandSideBoundedVector386 - crRightHandSideBoundedVector74*(crRightHandSideBoundedVector90 + crRightHandSideBoundedVector91));
const double crRightHandSideBoundedVector388 =             crRightHandSideBoundedVector29 + crRightHandSideBoundedVector387 + crRightHandSideBoundedVector86 + crRightHandSideBoundedVector88;
const double crRightHandSideBoundedVector389 =             crRightHandSideBoundedVector74*crRightHandSideBoundedVector83;
const double crRightHandSideBoundedVector390 =             crRightHandSideBoundedVector303*crRightHandSideBoundedVector389;
const double crRightHandSideBoundedVector391 =             crRightHandSideBoundedVector30 - crRightHandSideBoundedVector387 + crRightHandSideBoundedVector87 + crRightHandSideBoundedVector89;
const double crRightHandSideBoundedVector392 =             crRightHandSideBoundedVector391 + crRightHandSideBoundedVector65*crRightHandSideBoundedVector74*crRightHandSideBoundedVector84;
const double crRightHandSideBoundedVector393 =             crRightHandSideBoundedVector285*crRightHandSideBoundedVector392*crRightHandSideBoundedVector64;
const double crRightHandSideBoundedVector394 =             crRightHandSideBoundedVector278*crRightHandSideBoundedVector392*crRightHandSideBoundedVector61;
const double crRightHandSideBoundedVector395 =             -crRightHandSideBoundedVector105*crRightHandSideBoundedVector158 - crRightHandSideBoundedVector105*crRightHandSideBoundedVector55 + crRightHandSideBoundedVector354 + crRightHandSideBoundedVector377 - crRightHandSideBoundedVector379 - crRightHandSideBoundedVector380 - crRightHandSideBoundedVector381 + crRightHandSideBoundedVector382*crRightHandSideBoundedVector78 + crRightHandSideBoundedVector383*crRightHandSideBoundedVector82 + crRightHandSideBoundedVector393 + crRightHandSideBoundedVector394 + crRightHandSideBoundedVector48*crRightHandSideBoundedVector74*(crRightHandSideBoundedVector388 - crRightHandSideBoundedVector390) + crRightHandSideBoundedVector59*crRightHandSideBoundedVector74*(-crRightHandSideBoundedVector385 + crRightHandSideBoundedVector388) + 0.666666666666667*dUdt_1_3;
const double crRightHandSideBoundedVector396 =             crRightHandSideBoundedVector376*crRightHandSideBoundedVector395;
const double crRightHandSideBoundedVector397 =             1.0/(crRightHandSideBoundedVector114*crRightHandSideBoundedVector352 + crRightHandSideBoundedVector129);
const double crRightHandSideBoundedVector398 =             crRightHandSideBoundedVector113*(crRightHandSideBoundedVector357 + crRightHandSideBoundedVector378 + 0.666666666666667*r[0]);
const double crRightHandSideBoundedVector399 =             crRightHandSideBoundedVector117*crRightHandSideBoundedVector131;
const double crRightHandSideBoundedVector400 =             crRightHandSideBoundedVector120*crRightHandSideBoundedVector179;
const double crRightHandSideBoundedVector401 =             crRightHandSideBoundedVector114*crRightHandSideBoundedVector40*gamma;
const double crRightHandSideBoundedVector402 =             crRightHandSideBoundedVector114*crRightHandSideBoundedVector148*gamma;
const double crRightHandSideBoundedVector403 =             crRightHandSideBoundedVector114*crRightHandSideBoundedVector118;
const double crRightHandSideBoundedVector404 =             crRightHandSideBoundedVector303*crRightHandSideBoundedVector403;
const double crRightHandSideBoundedVector405 =             crRightHandSideBoundedVector124 + crRightHandSideBoundedVector31 + crRightHandSideBoundedVector88;
const double crRightHandSideBoundedVector406 =             crRightHandSideBoundedVector27*(-crRightHandSideBoundedVector114*(crRightHandSideBoundedVector126 + crRightHandSideBoundedVector127) + crRightHandSideBoundedVector405);
const double crRightHandSideBoundedVector407 =             crRightHandSideBoundedVector124 + crRightHandSideBoundedVector31 + crRightHandSideBoundedVector406 + crRightHandSideBoundedVector88;
const double crRightHandSideBoundedVector408 =             crRightHandSideBoundedVector114*crRightHandSideBoundedVector121;
const double crRightHandSideBoundedVector409 =             crRightHandSideBoundedVector303*crRightHandSideBoundedVector408;
const double crRightHandSideBoundedVector410 =             crRightHandSideBoundedVector125 + crRightHandSideBoundedVector32 - crRightHandSideBoundedVector406 + crRightHandSideBoundedVector89;
const double crRightHandSideBoundedVector411 =             crRightHandSideBoundedVector114*crRightHandSideBoundedVector122*crRightHandSideBoundedVector65 + crRightHandSideBoundedVector410;
const double crRightHandSideBoundedVector412 =             crRightHandSideBoundedVector210*crRightHandSideBoundedVector411*crRightHandSideBoundedVector64;
const double crRightHandSideBoundedVector413 =             crRightHandSideBoundedVector297*crRightHandSideBoundedVector411*crRightHandSideBoundedVector61;
const double crRightHandSideBoundedVector414 =             crRightHandSideBoundedVector114*crRightHandSideBoundedVector48*(crRightHandSideBoundedVector407 - crRightHandSideBoundedVector409) + crRightHandSideBoundedVector114*crRightHandSideBoundedVector59*(-crRightHandSideBoundedVector404 + crRightHandSideBoundedVector407) + crRightHandSideBoundedVector117*crRightHandSideBoundedVector401 + crRightHandSideBoundedVector120*crRightHandSideBoundedVector402 - crRightHandSideBoundedVector139*crRightHandSideBoundedVector158 - crRightHandSideBoundedVector139*crRightHandSideBoundedVector55 + crRightHandSideBoundedVector355 + crRightHandSideBoundedVector377 - crRightHandSideBoundedVector398 - crRightHandSideBoundedVector399 - crRightHandSideBoundedVector400 + crRightHandSideBoundedVector412 + crRightHandSideBoundedVector413 + 0.666666666666667*dUdt_0_3;
const double crRightHandSideBoundedVector415 =             crRightHandSideBoundedVector397*crRightHandSideBoundedVector414;
const double crRightHandSideBoundedVector416 =             0.166666666666667*crRightHandSideBoundedVector169;
const double crRightHandSideBoundedVector417 =             -0.166666666666667*crRightHandSideBoundedVector170;
const double crRightHandSideBoundedVector418 =             -0.166666666666667*crRightHandSideBoundedVector171;
const double crRightHandSideBoundedVector419 =             crRightHandSideBoundedVector27*crRightHandSideBoundedVector51;
const double crRightHandSideBoundedVector420 =             0.166666666666667*crRightHandSideBoundedVector419;
const double crRightHandSideBoundedVector421 =             crRightHandSideBoundedVector420*crRightHandSideBoundedVector99;
const double crRightHandSideBoundedVector422 =             crRightHandSideBoundedVector196*crRightHandSideBoundedVector48;
const double crRightHandSideBoundedVector423 =             0.166666666666667*crRightHandSideBoundedVector422;
const double crRightHandSideBoundedVector424 =             crRightHandSideBoundedVector101*crRightHandSideBoundedVector423;
const double crRightHandSideBoundedVector425 =             0.166666666666667*crRightHandSideBoundedVector173;
const double crRightHandSideBoundedVector426 =             -0.166666666666667*crRightHandSideBoundedVector176;
const double crRightHandSideBoundedVector427 =             0.166666666666667*crRightHandSideBoundedVector155;
const double crRightHandSideBoundedVector428 =             -0.166666666666667*crRightHandSideBoundedVector156;
const double crRightHandSideBoundedVector429 =             -0.166666666666667*crRightHandSideBoundedVector157;
const double crRightHandSideBoundedVector430 =             crRightHandSideBoundedVector420*crRightHandSideBoundedVector49;
const double crRightHandSideBoundedVector431 =             0.666666666666667*crRightHandSideBoundedVector419;
const double crRightHandSideBoundedVector432 =             crRightHandSideBoundedVector423*crRightHandSideBoundedVector52;
const double crRightHandSideBoundedVector433 =             0.666666666666667*crRightHandSideBoundedVector422;
const double crRightHandSideBoundedVector434 =             0.166666666666667*crRightHandSideBoundedVector160;
const double crRightHandSideBoundedVector435 =             crRightHandSideBoundedVector117*crRightHandSideBoundedVector214;
const double crRightHandSideBoundedVector436 =             -0.166666666666667*crRightHandSideBoundedVector163;
const double crRightHandSideBoundedVector437 =             DN_DX_0_0*crRightHandSideBoundedVector215*mu;
const double crRightHandSideBoundedVector438 =             pow(lin_m[1], 2);
const double crRightHandSideBoundedVector439 =             crRightHandSideBoundedVector217*crRightHandSideBoundedVector218*crRightHandSideBoundedVector438*nu_st;
const double crRightHandSideBoundedVector440 =             -crRightHandSideBoundedVector217*crRightHandSideBoundedVector438 + 1;
const double crRightHandSideBoundedVector441 =             crRightHandSideBoundedVector218*crRightHandSideBoundedVector440*nu_sc;
const double crRightHandSideBoundedVector442 =             crRightHandSideBoundedVector229*crRightHandSideBoundedVector241 + crRightHandSideBoundedVector235*mu*(crRightHandSideBoundedVector12*crRightHandSideBoundedVector439 + crRightHandSideBoundedVector12*crRightHandSideBoundedVector441 + 1);
const double crRightHandSideBoundedVector443 =             crRightHandSideBoundedVector243*crRightHandSideBoundedVector248 + crRightHandSideBoundedVector247*mu*(crRightHandSideBoundedVector439*crRightHandSideBoundedVector73 + crRightHandSideBoundedVector441*crRightHandSideBoundedVector73 + 1);
const double crRightHandSideBoundedVector444 =             crRightHandSideBoundedVector250*crRightHandSideBoundedVector254 + crRightHandSideBoundedVector253*mu*(crRightHandSideBoundedVector113*crRightHandSideBoundedVector439 + crRightHandSideBoundedVector113*crRightHandSideBoundedVector441 + 1);
const double crRightHandSideBoundedVector445 =             crRightHandSideBoundedVector261*crRightHandSideBoundedVector325;
const double crRightHandSideBoundedVector446 =             0.0277777777777778*f_ext(0,1);
const double crRightHandSideBoundedVector447 =             0.0277777777777778*f_ext(1,1);
const double crRightHandSideBoundedVector448 =             0.111111111111111*f_ext(2,1);
const double crRightHandSideBoundedVector449 =             crRightHandSideBoundedVector14*crRightHandSideBoundedVector18*crRightHandSideBoundedVector54;
const double crRightHandSideBoundedVector450 =             crRightHandSideBoundedVector14*crRightHandSideBoundedVector23*crRightHandSideBoundedVector59;
const double crRightHandSideBoundedVector451 =             0.166666666666667*crRightHandSideBoundedVector450;
const double crRightHandSideBoundedVector452 =             crRightHandSideBoundedVector269*crRightHandSideBoundedVector420;
const double crRightHandSideBoundedVector453 =             crRightHandSideBoundedVector18*crRightHandSideBoundedVector23*crRightHandSideBoundedVector270*crRightHandSideBoundedVector64;
const double crRightHandSideBoundedVector454 =             0.333333333333333*crRightHandSideBoundedVector453;
const double crRightHandSideBoundedVector455 =             crRightHandSideBoundedVector161*crRightHandSideBoundedVector270*crRightHandSideBoundedVector61;
const double crRightHandSideBoundedVector456 =             crRightHandSideBoundedVector261*crRightHandSideBoundedVector423 + crRightHandSideBoundedVector446 + crRightHandSideBoundedVector447 + crRightHandSideBoundedVector448 - 0.166666666666667*crRightHandSideBoundedVector449 - crRightHandSideBoundedVector451 + crRightHandSideBoundedVector452 + crRightHandSideBoundedVector454 - 0.333333333333333*crRightHandSideBoundedVector455;
const double crRightHandSideBoundedVector457 =             crRightHandSideBoundedVector278*crRightHandSideBoundedVector336;
const double crRightHandSideBoundedVector458 =             0.111111111111111*f_ext(1,1);
const double crRightHandSideBoundedVector459 =             0.0277777777777778*f_ext(2,1);
const double crRightHandSideBoundedVector460 =             crRightHandSideBoundedVector54*crRightHandSideBoundedVector75*crRightHandSideBoundedVector78;
const double crRightHandSideBoundedVector461 =             crRightHandSideBoundedVector59*crRightHandSideBoundedVector75*crRightHandSideBoundedVector82;
const double crRightHandSideBoundedVector462 =             0.166666666666667*crRightHandSideBoundedVector461;
const double crRightHandSideBoundedVector463 =             crRightHandSideBoundedVector285*crRightHandSideBoundedVector420;
const double crRightHandSideBoundedVector464 =             crRightHandSideBoundedVector286*crRightHandSideBoundedVector64*crRightHandSideBoundedVector78*crRightHandSideBoundedVector82;
const double crRightHandSideBoundedVector465 =             0.333333333333333*crRightHandSideBoundedVector464;
const double crRightHandSideBoundedVector466 =             crRightHandSideBoundedVector174*crRightHandSideBoundedVector286*crRightHandSideBoundedVector61;
const double crRightHandSideBoundedVector467 =             crRightHandSideBoundedVector278*crRightHandSideBoundedVector423 + crRightHandSideBoundedVector446 + crRightHandSideBoundedVector458 + crRightHandSideBoundedVector459 - 0.166666666666667*crRightHandSideBoundedVector460 - crRightHandSideBoundedVector462 + crRightHandSideBoundedVector463 + crRightHandSideBoundedVector465 - 0.333333333333333*crRightHandSideBoundedVector466;
const double crRightHandSideBoundedVector468 =             crRightHandSideBoundedVector448 + crRightHandSideBoundedVector458 + 0.444444444444444*f_ext(0,1);
const double crRightHandSideBoundedVector469 =             crRightHandSideBoundedVector115*crRightHandSideBoundedVector117*crRightHandSideBoundedVector54;
const double crRightHandSideBoundedVector470 =             crRightHandSideBoundedVector115*crRightHandSideBoundedVector120*crRightHandSideBoundedVector59;
const double crRightHandSideBoundedVector471 =             0.666666666666667*crRightHandSideBoundedVector470;
const double crRightHandSideBoundedVector472 =             crRightHandSideBoundedVector297*crRightHandSideBoundedVector347;
const double crRightHandSideBoundedVector473 =             crRightHandSideBoundedVector210*crRightHandSideBoundedVector431;
const double crRightHandSideBoundedVector474 =             crRightHandSideBoundedVector117*crRightHandSideBoundedVector120*crRightHandSideBoundedVector299*crRightHandSideBoundedVector64;
const double crRightHandSideBoundedVector475 =             1.33333333333333*crRightHandSideBoundedVector474;
const double crRightHandSideBoundedVector476 =             crRightHandSideBoundedVector185*crRightHandSideBoundedVector299*crRightHandSideBoundedVector61;
const double crRightHandSideBoundedVector477 =             DN_DX_0_0*crRightHandSideBoundedVector20;
const double crRightHandSideBoundedVector478 =             DN_DX_1_0*crRightHandSideBoundedVector21;
const double crRightHandSideBoundedVector479 =             DN_DX_2_0*crRightHandSideBoundedVector81;
const double crRightHandSideBoundedVector480 =             -crRightHandSideBoundedVector420;
const double crRightHandSideBoundedVector481 =             crRightHandSideBoundedVector27*crRightHandSideBoundedVector309 - crRightHandSideBoundedVector311 + crRightHandSideBoundedVector477 + crRightHandSideBoundedVector478 + crRightHandSideBoundedVector479 + crRightHandSideBoundedVector480;
const double crRightHandSideBoundedVector482 =             crRightHandSideBoundedVector27*crRightHandSideBoundedVector316 - crRightHandSideBoundedVector318 + crRightHandSideBoundedVector477 + crRightHandSideBoundedVector478 + crRightHandSideBoundedVector479 + crRightHandSideBoundedVector480;
const double crRightHandSideBoundedVector483 =             DN_DX_0_0*crRightHandSideBoundedVector119 + DN_DX_1_0*crRightHandSideBoundedVector80 + DN_DX_2_0*crRightHandSideBoundedVector22 - crRightHandSideBoundedVector431;
const double crRightHandSideBoundedVector484 =             DN_DX_0_0*crRightHandSideBoundedVector15;
const double crRightHandSideBoundedVector485 =             DN_DX_1_0*crRightHandSideBoundedVector16;
const double crRightHandSideBoundedVector486 =             DN_DX_2_0*crRightHandSideBoundedVector77;
const double crRightHandSideBoundedVector487 =             -crRightHandSideBoundedVector423;
const double crRightHandSideBoundedVector488 =             -crRightHandSideBoundedVector332;
const double crRightHandSideBoundedVector489 =             crRightHandSideBoundedVector196*crRightHandSideBoundedVector330 + crRightHandSideBoundedVector484 + crRightHandSideBoundedVector485 + crRightHandSideBoundedVector486 + crRightHandSideBoundedVector487 + crRightHandSideBoundedVector488;
const double crRightHandSideBoundedVector490 =             -crRightHandSideBoundedVector339;
const double crRightHandSideBoundedVector491 =             crRightHandSideBoundedVector196*crRightHandSideBoundedVector337 + crRightHandSideBoundedVector484 + crRightHandSideBoundedVector485 + crRightHandSideBoundedVector486 + crRightHandSideBoundedVector487 + crRightHandSideBoundedVector490;
const double crRightHandSideBoundedVector492 =             DN_DX_0_0*crRightHandSideBoundedVector116;
const double crRightHandSideBoundedVector493 =             DN_DX_1_0*crRightHandSideBoundedVector76;
const double crRightHandSideBoundedVector494 =             DN_DX_2_0*crRightHandSideBoundedVector17;
const double crRightHandSideBoundedVector495 =             -crRightHandSideBoundedVector433 + crRightHandSideBoundedVector492 + crRightHandSideBoundedVector493 + crRightHandSideBoundedVector494;
const double crRightHandSideBoundedVector496 =             -crRightHandSideBoundedVector114*crRightHandSideBoundedVector435;
const double crRightHandSideBoundedVector497 =             DN_DX_0_1*crRightHandSideBoundedVector303*h;
const double crRightHandSideBoundedVector498 =             0.166666666666667*crRightHandSideBoundedVector379;
const double crRightHandSideBoundedVector499 =             0.166666666666667*crRightHandSideBoundedVector380;
const double crRightHandSideBoundedVector500 =             0.166666666666667*crRightHandSideBoundedVector381;
const double crRightHandSideBoundedVector501 =             0.166666666666667*crRightHandSideBoundedVector382;
const double crRightHandSideBoundedVector502 =             -crRightHandSideBoundedVector501*crRightHandSideBoundedVector78;
const double crRightHandSideBoundedVector503 =             0.166666666666667*crRightHandSideBoundedVector383;
const double crRightHandSideBoundedVector504 =             -crRightHandSideBoundedVector503*crRightHandSideBoundedVector82;
const double crRightHandSideBoundedVector505 =             crRightHandSideBoundedVector284*crRightHandSideBoundedVector78;
const double crRightHandSideBoundedVector506 =             crRightHandSideBoundedVector278*crRightHandSideBoundedVector420;
const double crRightHandSideBoundedVector507 =             crRightHandSideBoundedVector506*crRightHandSideBoundedVector78;
const double crRightHandSideBoundedVector508 =             crRightHandSideBoundedVector385 + crRightHandSideBoundedVector391;
const double crRightHandSideBoundedVector509 =             crRightHandSideBoundedVector508*crRightHandSideBoundedVector59*crRightHandSideBoundedVector74;
const double crRightHandSideBoundedVector510 =             0.166666666666667*crRightHandSideBoundedVector509;
const double crRightHandSideBoundedVector511 =             crRightHandSideBoundedVector390 + crRightHandSideBoundedVector391;
const double crRightHandSideBoundedVector512 =             crRightHandSideBoundedVector48*crRightHandSideBoundedVector511*crRightHandSideBoundedVector74;
const double crRightHandSideBoundedVector513 =             0.166666666666667*crRightHandSideBoundedVector512;
const double crRightHandSideBoundedVector514 =             -0.166666666666667*crRightHandSideBoundedVector393;
const double crRightHandSideBoundedVector515 =             -0.166666666666667*crRightHandSideBoundedVector394;
const double crRightHandSideBoundedVector516 =             0.166666666666667*crRightHandSideBoundedVector358;
const double crRightHandSideBoundedVector517 =             0.166666666666667*crRightHandSideBoundedVector359;
const double crRightHandSideBoundedVector518 =             0.166666666666667*crRightHandSideBoundedVector360;
const double crRightHandSideBoundedVector519 =             0.166666666666667*crRightHandSideBoundedVector361;
const double crRightHandSideBoundedVector520 =             -crRightHandSideBoundedVector18*crRightHandSideBoundedVector519;
const double crRightHandSideBoundedVector521 =             0.166666666666667*crRightHandSideBoundedVector362;
const double crRightHandSideBoundedVector522 =             -crRightHandSideBoundedVector23*crRightHandSideBoundedVector521;
const double crRightHandSideBoundedVector523 =             0.666666666666667*crRightHandSideBoundedVector401;
const double crRightHandSideBoundedVector524 =             0.666666666666667*crRightHandSideBoundedVector402;
const double crRightHandSideBoundedVector525 =             crRightHandSideBoundedVector18*crRightHandSideBoundedVector268;
const double crRightHandSideBoundedVector526 =             crRightHandSideBoundedVector261*crRightHandSideBoundedVector420;
const double crRightHandSideBoundedVector527 =             crRightHandSideBoundedVector18*crRightHandSideBoundedVector526;
const double crRightHandSideBoundedVector528 =             crRightHandSideBoundedVector297*crRightHandSideBoundedVector431;
const double crRightHandSideBoundedVector529 =             crRightHandSideBoundedVector364 + crRightHandSideBoundedVector370;
const double crRightHandSideBoundedVector530 =             crRightHandSideBoundedVector13*crRightHandSideBoundedVector529*crRightHandSideBoundedVector59;
const double crRightHandSideBoundedVector531 =             0.166666666666667*crRightHandSideBoundedVector530;
const double crRightHandSideBoundedVector532 =             crRightHandSideBoundedVector369 + crRightHandSideBoundedVector370;
const double crRightHandSideBoundedVector533 =             crRightHandSideBoundedVector13*crRightHandSideBoundedVector48*crRightHandSideBoundedVector532;
const double crRightHandSideBoundedVector534 =             0.166666666666667*crRightHandSideBoundedVector533;
const double crRightHandSideBoundedVector535 =             crRightHandSideBoundedVector404 + crRightHandSideBoundedVector410;
const double crRightHandSideBoundedVector536 =             crRightHandSideBoundedVector114*crRightHandSideBoundedVector535*crRightHandSideBoundedVector59;
const double crRightHandSideBoundedVector537 =             crRightHandSideBoundedVector409 + crRightHandSideBoundedVector410;
const double crRightHandSideBoundedVector538 =             crRightHandSideBoundedVector114*crRightHandSideBoundedVector48*crRightHandSideBoundedVector537;
const double crRightHandSideBoundedVector539 =             -0.166666666666667*crRightHandSideBoundedVector372;
const double crRightHandSideBoundedVector540 =             -0.166666666666667*crRightHandSideBoundedVector373;
const double crRightHandSideBoundedVector541 =             crRightHandSideBoundedVector115*crRightHandSideBoundedVector411;
const double crRightHandSideBoundedVector542 =             crRightHandSideBoundedVector218*nu_sc*(-crRightHandSideBoundedVector221 + 1);
const double crRightHandSideBoundedVector543 =             crRightHandSideBoundedVector215*mu*(crRightHandSideBoundedVector12*crRightHandSideBoundedVector542 + crRightHandSideBoundedVector220);
const double crRightHandSideBoundedVector544 =             crRightHandSideBoundedVector217*lin_m[0]*lin_m[1]*(-k_sc + k_st);
const double crRightHandSideBoundedVector545 =             crRightHandSideBoundedVector12*crRightHandSideBoundedVector148 - crRightHandSideBoundedVector13*(crRightHandSideBoundedVector159 + crRightHandSideBoundedVector18*crRightHandSideBoundedVector51 - crRightHandSideBoundedVector25*crRightHandSideBoundedVector61) - crRightHandSideBoundedVector365*crRightHandSideBoundedVector61;
const double crRightHandSideBoundedVector546 =             crRightHandSideBoundedVector13*crRightHandSideBoundedVector351*lambda;
const double crRightHandSideBoundedVector547 =             1.0/lambda;
const double crRightHandSideBoundedVector548 =             c_v*crRightHandSideBoundedVector217*crRightHandSideBoundedVector236*crRightHandSideBoundedVector547*k_st;
const double crRightHandSideBoundedVector549 =             c_v*crRightHandSideBoundedVector238*crRightHandSideBoundedVector547*k_sc;
const double crRightHandSideBoundedVector550 =             crRightHandSideBoundedVector12*crRightHandSideBoundedVector40 - crRightHandSideBoundedVector13*(crRightHandSideBoundedVector23*crRightHandSideBoundedVector54 - crRightHandSideBoundedVector25*crRightHandSideBoundedVector64 + crRightHandSideBoundedVector60) - crRightHandSideBoundedVector365*crRightHandSideBoundedVector64;
const double crRightHandSideBoundedVector551 =             crRightHandSideBoundedVector13*(crRightHandSideBoundedVector18*crRightHandSideBoundedVector242 + crRightHandSideBoundedVector23*crRightHandSideBoundedVector543 - crRightHandSideBoundedVector544*crRightHandSideBoundedVector545 - crRightHandSideBoundedVector546*crRightHandSideBoundedVector550*(crRightHandSideBoundedVector12*crRightHandSideBoundedVector548 + crRightHandSideBoundedVector12*crRightHandSideBoundedVector549 + 1));
const double crRightHandSideBoundedVector552 =             crRightHandSideBoundedVector215*mu*(crRightHandSideBoundedVector224 + crRightHandSideBoundedVector542*crRightHandSideBoundedVector73);
const double crRightHandSideBoundedVector553 =             crRightHandSideBoundedVector148*crRightHandSideBoundedVector73 - crRightHandSideBoundedVector386*crRightHandSideBoundedVector61 - crRightHandSideBoundedVector74*(crRightHandSideBoundedVector172 + crRightHandSideBoundedVector51*crRightHandSideBoundedVector78 - crRightHandSideBoundedVector61*crRightHandSideBoundedVector84);
const double crRightHandSideBoundedVector554 =             crRightHandSideBoundedVector351*crRightHandSideBoundedVector74*lambda;
const double crRightHandSideBoundedVector555 =             -crRightHandSideBoundedVector386*crRightHandSideBoundedVector64 + crRightHandSideBoundedVector40*crRightHandSideBoundedVector73 - crRightHandSideBoundedVector74*(crRightHandSideBoundedVector104 + crRightHandSideBoundedVector54*crRightHandSideBoundedVector82 - crRightHandSideBoundedVector64*crRightHandSideBoundedVector84);
const double crRightHandSideBoundedVector556 =             crRightHandSideBoundedVector74*(crRightHandSideBoundedVector249*crRightHandSideBoundedVector78 - crRightHandSideBoundedVector544*crRightHandSideBoundedVector553 + crRightHandSideBoundedVector552*crRightHandSideBoundedVector82 - crRightHandSideBoundedVector554*crRightHandSideBoundedVector555*(crRightHandSideBoundedVector548*crRightHandSideBoundedVector73 + crRightHandSideBoundedVector549*crRightHandSideBoundedVector73 + 1));
const double crRightHandSideBoundedVector557 =             crRightHandSideBoundedVector215*mu*(crRightHandSideBoundedVector113*crRightHandSideBoundedVector542 + crRightHandSideBoundedVector226);
const double crRightHandSideBoundedVector558 =             crRightHandSideBoundedVector113*crRightHandSideBoundedVector148 - crRightHandSideBoundedVector114*(crRightHandSideBoundedVector117*crRightHandSideBoundedVector51 - crRightHandSideBoundedVector122*crRightHandSideBoundedVector61 + crRightHandSideBoundedVector183) - crRightHandSideBoundedVector405*crRightHandSideBoundedVector61;
const double crRightHandSideBoundedVector559 =             crRightHandSideBoundedVector114*crRightHandSideBoundedVector351*lambda;
const double crRightHandSideBoundedVector560 =             crRightHandSideBoundedVector113*crRightHandSideBoundedVector40 - crRightHandSideBoundedVector114*(crRightHandSideBoundedVector120*crRightHandSideBoundedVector54 - crRightHandSideBoundedVector122*crRightHandSideBoundedVector64 + crRightHandSideBoundedVector138) - crRightHandSideBoundedVector405*crRightHandSideBoundedVector64;
const double crRightHandSideBoundedVector561 =             crRightHandSideBoundedVector114*(crRightHandSideBoundedVector117*crRightHandSideBoundedVector255 + crRightHandSideBoundedVector120*crRightHandSideBoundedVector557 - crRightHandSideBoundedVector544*crRightHandSideBoundedVector558 - crRightHandSideBoundedVector559*crRightHandSideBoundedVector560*(crRightHandSideBoundedVector113*crRightHandSideBoundedVector548 + crRightHandSideBoundedVector113*crRightHandSideBoundedVector549 + 1));
const double crRightHandSideBoundedVector562 =             c_v*crRightHandSideBoundedVector217*crRightHandSideBoundedVector438*crRightHandSideBoundedVector547*k_st;
const double crRightHandSideBoundedVector563 =             c_v*crRightHandSideBoundedVector440*crRightHandSideBoundedVector547*k_sc;
const double crRightHandSideBoundedVector564 =             crRightHandSideBoundedVector13*(crRightHandSideBoundedVector18*crRightHandSideBoundedVector543 + crRightHandSideBoundedVector23*crRightHandSideBoundedVector442 - crRightHandSideBoundedVector544*crRightHandSideBoundedVector550 - crRightHandSideBoundedVector545*crRightHandSideBoundedVector546*(crRightHandSideBoundedVector12*crRightHandSideBoundedVector562 + crRightHandSideBoundedVector12*crRightHandSideBoundedVector563 + 1));
const double crRightHandSideBoundedVector565 =             crRightHandSideBoundedVector74*(crRightHandSideBoundedVector443*crRightHandSideBoundedVector82 - crRightHandSideBoundedVector544*crRightHandSideBoundedVector555 + crRightHandSideBoundedVector552*crRightHandSideBoundedVector78 - crRightHandSideBoundedVector553*crRightHandSideBoundedVector554*(crRightHandSideBoundedVector562*crRightHandSideBoundedVector73 + crRightHandSideBoundedVector563*crRightHandSideBoundedVector73 + 1));
const double crRightHandSideBoundedVector566 =             crRightHandSideBoundedVector114*(crRightHandSideBoundedVector117*crRightHandSideBoundedVector557 + crRightHandSideBoundedVector120*crRightHandSideBoundedVector444 - crRightHandSideBoundedVector544*crRightHandSideBoundedVector560 - crRightHandSideBoundedVector558*crRightHandSideBoundedVector559*(crRightHandSideBoundedVector113*crRightHandSideBoundedVector562 + crRightHandSideBoundedVector113*crRightHandSideBoundedVector563 + 1));
const double crRightHandSideBoundedVector567 =             1.0*crRightHandSideBoundedVector39*crRightHandSideBoundedVector70*h;
const double crRightHandSideBoundedVector568 =             crRightHandSideBoundedVector13*crRightHandSideBoundedVector529;
const double crRightHandSideBoundedVector569 =             crRightHandSideBoundedVector14*crRightHandSideBoundedVector18*crRightHandSideBoundedVector59;
const double crRightHandSideBoundedVector570 =             2.0*gamma - 2.0;
const double crRightHandSideBoundedVector571 =             crRightHandSideBoundedVector363*crRightHandSideBoundedVector570 + crRightHandSideBoundedVector371;
const double crRightHandSideBoundedVector572 =             0.166666666666667*crRightHandSideBoundedVector14*crRightHandSideBoundedVector571;
const double crRightHandSideBoundedVector573 =             crRightHandSideBoundedVector262 + crRightHandSideBoundedVector263 + crRightHandSideBoundedVector264 - crRightHandSideBoundedVector266*crRightHandSideBoundedVector27 - crRightHandSideBoundedVector268 + crRightHandSideBoundedVector27*crRightHandSideBoundedVector272 + crRightHandSideBoundedVector519 - crRightHandSideBoundedVector526 - crRightHandSideBoundedVector569*crRightHandSideBoundedVector65 + crRightHandSideBoundedVector572*crRightHandSideBoundedVector64;
const double crRightHandSideBoundedVector574 =             1.0*crRightHandSideBoundedVector164*crRightHandSideBoundedVector39*h;
const double crRightHandSideBoundedVector575 =             crRightHandSideBoundedVector13*crRightHandSideBoundedVector532;
const double crRightHandSideBoundedVector576 =             crRightHandSideBoundedVector14*crRightHandSideBoundedVector23*crRightHandSideBoundedVector48;
const double crRightHandSideBoundedVector577 =             crRightHandSideBoundedVector368*crRightHandSideBoundedVector570 + crRightHandSideBoundedVector371;
const double crRightHandSideBoundedVector578 =             0.166666666666667*crRightHandSideBoundedVector14*crRightHandSideBoundedVector577;
const double crRightHandSideBoundedVector579 =             -crRightHandSideBoundedVector194*crRightHandSideBoundedVector269 - crRightHandSideBoundedVector27*crRightHandSideBoundedVector451 + crRightHandSideBoundedVector27*crRightHandSideBoundedVector454 + crRightHandSideBoundedVector446 + crRightHandSideBoundedVector447 + crRightHandSideBoundedVector448 - crRightHandSideBoundedVector452 + crRightHandSideBoundedVector521 - crRightHandSideBoundedVector576*crRightHandSideBoundedVector65 + crRightHandSideBoundedVector578*crRightHandSideBoundedVector61;
const double crRightHandSideBoundedVector580 =             1.0*crRightHandSideBoundedVector111*crRightHandSideBoundedVector94*h;
const double crRightHandSideBoundedVector581 =             crRightHandSideBoundedVector508*crRightHandSideBoundedVector74;
const double crRightHandSideBoundedVector582 =             crRightHandSideBoundedVector59*crRightHandSideBoundedVector75*crRightHandSideBoundedVector78;
const double crRightHandSideBoundedVector583 =             crRightHandSideBoundedVector384*crRightHandSideBoundedVector570 + crRightHandSideBoundedVector392;
const double crRightHandSideBoundedVector584 =             0.166666666666667*crRightHandSideBoundedVector583*crRightHandSideBoundedVector75;
const double crRightHandSideBoundedVector585 =             crRightHandSideBoundedVector262 - crRightHandSideBoundedVector27*crRightHandSideBoundedVector282 + crRightHandSideBoundedVector27*crRightHandSideBoundedVector288 + crRightHandSideBoundedVector279 + crRightHandSideBoundedVector280 - crRightHandSideBoundedVector284 + crRightHandSideBoundedVector501 - crRightHandSideBoundedVector506 - crRightHandSideBoundedVector582*crRightHandSideBoundedVector65 + crRightHandSideBoundedVector584*crRightHandSideBoundedVector64;
const double crRightHandSideBoundedVector586 =             1.0*crRightHandSideBoundedVector177*crRightHandSideBoundedVector94*h;
const double crRightHandSideBoundedVector587 =             crRightHandSideBoundedVector511*crRightHandSideBoundedVector74;
const double crRightHandSideBoundedVector588 =             crRightHandSideBoundedVector48*crRightHandSideBoundedVector75*crRightHandSideBoundedVector82;
const double crRightHandSideBoundedVector589 =             crRightHandSideBoundedVector389*crRightHandSideBoundedVector570 + crRightHandSideBoundedVector392;
const double crRightHandSideBoundedVector590 =             0.166666666666667*crRightHandSideBoundedVector589*crRightHandSideBoundedVector75;
const double crRightHandSideBoundedVector591 =             -crRightHandSideBoundedVector194*crRightHandSideBoundedVector285 - crRightHandSideBoundedVector27*crRightHandSideBoundedVector462 + crRightHandSideBoundedVector27*crRightHandSideBoundedVector465 + crRightHandSideBoundedVector446 + crRightHandSideBoundedVector458 + crRightHandSideBoundedVector459 - crRightHandSideBoundedVector463 + crRightHandSideBoundedVector503 - crRightHandSideBoundedVector588*crRightHandSideBoundedVector65 + crRightHandSideBoundedVector590*crRightHandSideBoundedVector61;
const double crRightHandSideBoundedVector592 =             1.0*crRightHandSideBoundedVector130*crRightHandSideBoundedVector145*h;
const double crRightHandSideBoundedVector593 =             crRightHandSideBoundedVector115*crRightHandSideBoundedVector117*crRightHandSideBoundedVector59;
const double crRightHandSideBoundedVector594 =             crRightHandSideBoundedVector114*crRightHandSideBoundedVector535;
const double crRightHandSideBoundedVector595 =             crRightHandSideBoundedVector403*crRightHandSideBoundedVector570 + crRightHandSideBoundedVector411;
const double crRightHandSideBoundedVector596 =             crRightHandSideBoundedVector115*crRightHandSideBoundedVector595;
const double crRightHandSideBoundedVector597 =             1.0*crRightHandSideBoundedVector130*crRightHandSideBoundedVector188*h;
const double crRightHandSideBoundedVector598 =             crRightHandSideBoundedVector115*crRightHandSideBoundedVector120*crRightHandSideBoundedVector48;
const double crRightHandSideBoundedVector599 =             crRightHandSideBoundedVector114*crRightHandSideBoundedVector537;
const double crRightHandSideBoundedVector600 =             crRightHandSideBoundedVector408*crRightHandSideBoundedVector570 + crRightHandSideBoundedVector411;
const double crRightHandSideBoundedVector601 =             crRightHandSideBoundedVector115*crRightHandSideBoundedVector600;
const double crRightHandSideBoundedVector602 =             crRightHandSideBoundedVector326 + crRightHandSideBoundedVector327 + crRightHandSideBoundedVector328 + crRightHandSideBoundedVector331 + crRightHandSideBoundedVector484 + crRightHandSideBoundedVector485 + crRightHandSideBoundedVector486 + crRightHandSideBoundedVector488;
const double crRightHandSideBoundedVector603 =             1.0*crRightHandSideBoundedVector13*crRightHandSideBoundedVector353*crRightHandSideBoundedVector374*gamma*h;
const double crRightHandSideBoundedVector604 =             crRightHandSideBoundedVector326 + crRightHandSideBoundedVector327 + crRightHandSideBoundedVector328 + crRightHandSideBoundedVector338 + crRightHandSideBoundedVector484 + crRightHandSideBoundedVector485 + crRightHandSideBoundedVector486 + crRightHandSideBoundedVector490;
const double crRightHandSideBoundedVector605 =             1.0*crRightHandSideBoundedVector376*crRightHandSideBoundedVector395*crRightHandSideBoundedVector74*gamma*h;
const double crRightHandSideBoundedVector606 =             1.0*crRightHandSideBoundedVector114*crRightHandSideBoundedVector397*crRightHandSideBoundedVector414*gamma*h;
const double crRightHandSideBoundedVector607 =             crRightHandSideBoundedVector14*crRightHandSideBoundedVector371;
const double crRightHandSideBoundedVector608 =             0.0277777777777778*r[0];
const double crRightHandSideBoundedVector609 =             0.0277777777777778*r[1];
const double crRightHandSideBoundedVector610 =             0.111111111111111*r[2];
const double crRightHandSideBoundedVector611 =             crRightHandSideBoundedVector14*crRightHandSideBoundedVector18*crRightHandSideBoundedVector40*gamma;
const double crRightHandSideBoundedVector612 =             crRightHandSideBoundedVector14*crRightHandSideBoundedVector148*crRightHandSideBoundedVector23*gamma;
const double crRightHandSideBoundedVector613 =             crRightHandSideBoundedVector18*crRightHandSideBoundedVector23*crRightHandSideBoundedVector27*crRightHandSideBoundedVector270*crRightHandSideBoundedVector54;
const double crRightHandSideBoundedVector614 =             crRightHandSideBoundedVector18*crRightHandSideBoundedVector23*crRightHandSideBoundedVector27*crRightHandSideBoundedVector270*crRightHandSideBoundedVector51;
const double crRightHandSideBoundedVector615 =             crRightHandSideBoundedVector25*crRightHandSideBoundedVector28 + crRightHandSideBoundedVector370;
const double crRightHandSideBoundedVector616 =             crRightHandSideBoundedVector18*crRightHandSideBoundedVector270*crRightHandSideBoundedVector615*crRightHandSideBoundedVector64;
const double crRightHandSideBoundedVector617 =             crRightHandSideBoundedVector23*crRightHandSideBoundedVector270*crRightHandSideBoundedVector61*crRightHandSideBoundedVector615;
const double crRightHandSideBoundedVector618 =             crRightHandSideBoundedVector48*crRightHandSideBoundedVector578 + crRightHandSideBoundedVector572*crRightHandSideBoundedVector59 + crRightHandSideBoundedVector608 + crRightHandSideBoundedVector609 + crRightHandSideBoundedVector610 - 0.166666666666667*crRightHandSideBoundedVector611 - 0.166666666666667*crRightHandSideBoundedVector612 + 0.333333333333333*crRightHandSideBoundedVector613 + 0.333333333333333*crRightHandSideBoundedVector614 - 0.333333333333333*crRightHandSideBoundedVector616 - 0.333333333333333*crRightHandSideBoundedVector617;
const double crRightHandSideBoundedVector619 =             crRightHandSideBoundedVector392*crRightHandSideBoundedVector75;
const double crRightHandSideBoundedVector620 =             0.111111111111111*r[1];
const double crRightHandSideBoundedVector621 =             0.0277777777777778*r[2];
const double crRightHandSideBoundedVector622 =             crRightHandSideBoundedVector40*crRightHandSideBoundedVector75*crRightHandSideBoundedVector78*gamma;
const double crRightHandSideBoundedVector623 =             crRightHandSideBoundedVector148*crRightHandSideBoundedVector75*crRightHandSideBoundedVector82*gamma;
const double crRightHandSideBoundedVector624 =             crRightHandSideBoundedVector27*crRightHandSideBoundedVector286*crRightHandSideBoundedVector54*crRightHandSideBoundedVector78*crRightHandSideBoundedVector82;
const double crRightHandSideBoundedVector625 =             crRightHandSideBoundedVector27*crRightHandSideBoundedVector286*crRightHandSideBoundedVector51*crRightHandSideBoundedVector78*crRightHandSideBoundedVector82;
const double crRightHandSideBoundedVector626 =             crRightHandSideBoundedVector391 + crRightHandSideBoundedVector84*crRightHandSideBoundedVector85;
const double crRightHandSideBoundedVector627 =             crRightHandSideBoundedVector286*crRightHandSideBoundedVector626*crRightHandSideBoundedVector64*crRightHandSideBoundedVector78;
const double crRightHandSideBoundedVector628 =             crRightHandSideBoundedVector286*crRightHandSideBoundedVector61*crRightHandSideBoundedVector626*crRightHandSideBoundedVector82;
const double crRightHandSideBoundedVector629 =             crRightHandSideBoundedVector48*crRightHandSideBoundedVector590 + crRightHandSideBoundedVector584*crRightHandSideBoundedVector59 + crRightHandSideBoundedVector608 + crRightHandSideBoundedVector620 + crRightHandSideBoundedVector621 - 0.166666666666667*crRightHandSideBoundedVector622 - 0.166666666666667*crRightHandSideBoundedVector623 + 0.333333333333333*crRightHandSideBoundedVector624 + 0.333333333333333*crRightHandSideBoundedVector625 - 0.333333333333333*crRightHandSideBoundedVector627 - 0.333333333333333*crRightHandSideBoundedVector628;
const double crRightHandSideBoundedVector630 =             crRightHandSideBoundedVector115*crRightHandSideBoundedVector117*crRightHandSideBoundedVector40*gamma;
const double crRightHandSideBoundedVector631 =             crRightHandSideBoundedVector115*crRightHandSideBoundedVector120*crRightHandSideBoundedVector148*gamma;
const double crRightHandSideBoundedVector632 =             crRightHandSideBoundedVector117*crRightHandSideBoundedVector120*crRightHandSideBoundedVector27*crRightHandSideBoundedVector299*crRightHandSideBoundedVector54;
const double crRightHandSideBoundedVector633 =             crRightHandSideBoundedVector117*crRightHandSideBoundedVector120*crRightHandSideBoundedVector27*crRightHandSideBoundedVector299*crRightHandSideBoundedVector51;
const double crRightHandSideBoundedVector634 =             crRightHandSideBoundedVector122*crRightHandSideBoundedVector123 + crRightHandSideBoundedVector410;
const double crRightHandSideBoundedVector635 =             crRightHandSideBoundedVector117*crRightHandSideBoundedVector299*crRightHandSideBoundedVector634*crRightHandSideBoundedVector64;
const double crRightHandSideBoundedVector636 =             crRightHandSideBoundedVector120*crRightHandSideBoundedVector299*crRightHandSideBoundedVector61*crRightHandSideBoundedVector634;
const double crRightHandSideBoundedVector637 =             crRightHandSideBoundedVector492 + crRightHandSideBoundedVector493 + crRightHandSideBoundedVector494;
const double crRightHandSideBoundedVector638 =             crRightHandSideBoundedVector342 + crRightHandSideBoundedVector343 + crRightHandSideBoundedVector344;
const double crRightHandSideBoundedVector639 =             1.0*crRightHandSideBoundedVector0 + 1.0*crRightHandSideBoundedVector1 + 1.0*crRightHandSideBoundedVector2 + 1.0*crRightHandSideBoundedVector3 + 1.0*crRightHandSideBoundedVector4 + 1.0*crRightHandSideBoundedVector5;
const double crRightHandSideBoundedVector640 =             1.0*DN_DX_1_0*h;
const double crRightHandSideBoundedVector641 =             1.0*DN_DX_1_1*h;
const double crRightHandSideBoundedVector642 =             0.166666666666667*crRightHandSideBoundedVector132 + crRightHandSideBoundedVector133*crRightHandSideBoundedVector198 - 0.166666666666667*crRightHandSideBoundedVector134 + crRightHandSideBoundedVector135*crRightHandSideBoundedVector194 - 0.166666666666667*crRightHandSideBoundedVector136 + 0.166666666666667*crRightHandSideBoundedVector140 - 0.166666666666667*crRightHandSideBoundedVector144 - 1.0*crRightHandSideBoundedVector41;
const double crRightHandSideBoundedVector643 =             crRightHandSideBoundedVector211*crRightHandSideBoundedVector82;
const double crRightHandSideBoundedVector644 =             DN_DX_1_1*crRightHandSideBoundedVector215*mu;
const double crRightHandSideBoundedVector645 =             DN_DX_1_1*crRightHandSideBoundedVector18;
const double crRightHandSideBoundedVector646 =             0.111111111111111*f_ext(0,0);
const double crRightHandSideBoundedVector647 =             crRightHandSideBoundedVector264 + crRightHandSideBoundedVector646 + 0.444444444444444*f_ext(1,0);
const double crRightHandSideBoundedVector648 =             0.666666666666667*crRightHandSideBoundedVector281;
const double crRightHandSideBoundedVector649 =             DN_DX_1_1*crRightHandSideBoundedVector78;
const double crRightHandSideBoundedVector650 =             crRightHandSideBoundedVector206*crRightHandSideBoundedVector278;
const double crRightHandSideBoundedVector651 =             1.33333333333333*crRightHandSideBoundedVector287;
const double crRightHandSideBoundedVector652 =             DN_DX_1_1*crRightHandSideBoundedVector117;
const double crRightHandSideBoundedVector653 =             0.166666666666667*crRightHandSideBoundedVector293;
const double crRightHandSideBoundedVector654 =             crRightHandSideBoundedVector194*crRightHandSideBoundedVector297;
const double crRightHandSideBoundedVector655 =             0.333333333333333*crRightHandSideBoundedVector300;
const double crRightHandSideBoundedVector656 =             crRightHandSideBoundedVector198*crRightHandSideBoundedVector210 + crRightHandSideBoundedVector263 + crRightHandSideBoundedVector280 - 0.166666666666667*crRightHandSideBoundedVector295 - 0.333333333333333*crRightHandSideBoundedVector302 + crRightHandSideBoundedVector646 - crRightHandSideBoundedVector653 + crRightHandSideBoundedVector654 + crRightHandSideBoundedVector655;
const double crRightHandSideBoundedVector657 =             DN_DX_1_0*crRightHandSideBoundedVector23;
const double crRightHandSideBoundedVector658 =             DN_DX_1_0*crRightHandSideBoundedVector82;
const double crRightHandSideBoundedVector659 =             DN_DX_1_0*crRightHandSideBoundedVector120;
const double crRightHandSideBoundedVector660 =             0.166666666666667*crRightHandSideBoundedVector114*crRightHandSideBoundedVector117;
const double crRightHandSideBoundedVector661 =             crRightHandSideBoundedVector61*crRightHandSideBoundedVector660;
const double crRightHandSideBoundedVector662 =             0.166666666666667*crRightHandSideBoundedVector114*crRightHandSideBoundedVector120;
const double crRightHandSideBoundedVector663 =             crRightHandSideBoundedVector64*crRightHandSideBoundedVector662;
const double crRightHandSideBoundedVector664 =             crRightHandSideBoundedVector194 - crRightHandSideBoundedVector27*crRightHandSideBoundedVector663 + crRightHandSideBoundedVector305 + crRightHandSideBoundedVector306 + crRightHandSideBoundedVector307 + crRightHandSideBoundedVector661;
const double crRightHandSideBoundedVector665 =             DN_DX_1_1*crRightHandSideBoundedVector23;
const double crRightHandSideBoundedVector666 =             DN_DX_1_0*crRightHandSideBoundedVector18;
const double crRightHandSideBoundedVector667 =             DN_DX_1_1*crRightHandSideBoundedVector82;
const double crRightHandSideBoundedVector668 =             DN_DX_1_0*crRightHandSideBoundedVector78;
const double crRightHandSideBoundedVector669 =             -crRightHandSideBoundedVector643*crRightHandSideBoundedVector74;
const double crRightHandSideBoundedVector670 =             DN_DX_1_1*crRightHandSideBoundedVector120;
const double crRightHandSideBoundedVector671 =             DN_DX_1_0*crRightHandSideBoundedVector117;
const double crRightHandSideBoundedVector672 =             crRightHandSideBoundedVector61*crRightHandSideBoundedVector662;
const double crRightHandSideBoundedVector673 =             -crRightHandSideBoundedVector672;
const double crRightHandSideBoundedVector674 =             crRightHandSideBoundedVector64*crRightHandSideBoundedVector660;
const double crRightHandSideBoundedVector675 =             crRightHandSideBoundedVector196*crRightHandSideBoundedVector674 + crRightHandSideBoundedVector326 + crRightHandSideBoundedVector327 + crRightHandSideBoundedVector328 + crRightHandSideBoundedVector329 + crRightHandSideBoundedVector673;
const double crRightHandSideBoundedVector676 =             DN_DX_1_0*crRightHandSideBoundedVector303*h;
const double crRightHandSideBoundedVector677 =             crRightHandSideBoundedVector133*crRightHandSideBoundedVector420 + crRightHandSideBoundedVector135*crRightHandSideBoundedVector423 - 1.0*crRightHandSideBoundedVector149 + 0.166666666666667*crRightHandSideBoundedVector180 - 0.166666666666667*crRightHandSideBoundedVector181 - 0.166666666666667*crRightHandSideBoundedVector182 + 0.166666666666667*crRightHandSideBoundedVector184 - 0.166666666666667*crRightHandSideBoundedVector187;
const double crRightHandSideBoundedVector678 =             crRightHandSideBoundedVector214*crRightHandSideBoundedVector78;
const double crRightHandSideBoundedVector679 =             DN_DX_1_0*crRightHandSideBoundedVector215*mu;
const double crRightHandSideBoundedVector680 =             crRightHandSideBoundedVector261*crRightHandSideBoundedVector666;
const double crRightHandSideBoundedVector681 =             0.111111111111111*f_ext(0,1);
const double crRightHandSideBoundedVector682 =             crRightHandSideBoundedVector448 + crRightHandSideBoundedVector681 + 0.444444444444444*f_ext(1,1);
const double crRightHandSideBoundedVector683 =             0.666666666666667*crRightHandSideBoundedVector461;
const double crRightHandSideBoundedVector684 =             crRightHandSideBoundedVector278*crRightHandSideBoundedVector668;
const double crRightHandSideBoundedVector685 =             crRightHandSideBoundedVector285*crRightHandSideBoundedVector431;
const double crRightHandSideBoundedVector686 =             1.33333333333333*crRightHandSideBoundedVector464;
const double crRightHandSideBoundedVector687 =             crRightHandSideBoundedVector297*crRightHandSideBoundedVector671;
const double crRightHandSideBoundedVector688 =             0.166666666666667*crRightHandSideBoundedVector470;
const double crRightHandSideBoundedVector689 =             crRightHandSideBoundedVector210*crRightHandSideBoundedVector420;
const double crRightHandSideBoundedVector690 =             0.333333333333333*crRightHandSideBoundedVector474;
const double crRightHandSideBoundedVector691 =             crRightHandSideBoundedVector297*crRightHandSideBoundedVector423 + crRightHandSideBoundedVector447 + crRightHandSideBoundedVector459 - 0.166666666666667*crRightHandSideBoundedVector469 - 0.333333333333333*crRightHandSideBoundedVector476 + crRightHandSideBoundedVector681 - crRightHandSideBoundedVector688 + crRightHandSideBoundedVector689 + crRightHandSideBoundedVector690;
const double crRightHandSideBoundedVector692 =             crRightHandSideBoundedVector27*crRightHandSideBoundedVector661 + crRightHandSideBoundedVector477 + crRightHandSideBoundedVector478 + crRightHandSideBoundedVector479 + crRightHandSideBoundedVector480 - crRightHandSideBoundedVector663;
const double crRightHandSideBoundedVector693 =             -crRightHandSideBoundedVector678*crRightHandSideBoundedVector74;
const double crRightHandSideBoundedVector694 =             -crRightHandSideBoundedVector674;
const double crRightHandSideBoundedVector695 =             crRightHandSideBoundedVector196*crRightHandSideBoundedVector672 + crRightHandSideBoundedVector484 + crRightHandSideBoundedVector485 + crRightHandSideBoundedVector486 + crRightHandSideBoundedVector487 + crRightHandSideBoundedVector694;
const double crRightHandSideBoundedVector696 =             DN_DX_1_1*crRightHandSideBoundedVector303*h;
const double crRightHandSideBoundedVector697 =             0.166666666666667*crRightHandSideBoundedVector401;
const double crRightHandSideBoundedVector698 =             0.166666666666667*crRightHandSideBoundedVector402;
const double crRightHandSideBoundedVector699 =             crRightHandSideBoundedVector297*crRightHandSideBoundedVector420;
const double crRightHandSideBoundedVector700 =             crRightHandSideBoundedVector117*crRightHandSideBoundedVector654 - crRightHandSideBoundedVector117*crRightHandSideBoundedVector697 + crRightHandSideBoundedVector117*crRightHandSideBoundedVector699 - crRightHandSideBoundedVector120*crRightHandSideBoundedVector698 + 0.166666666666667*crRightHandSideBoundedVector398 + 0.166666666666667*crRightHandSideBoundedVector399 + 0.166666666666667*crRightHandSideBoundedVector400 - 0.166666666666667*crRightHandSideBoundedVector412 - 0.166666666666667*crRightHandSideBoundedVector413 + 0.166666666666667*crRightHandSideBoundedVector536 + 0.166666666666667*crRightHandSideBoundedVector538;
const double crRightHandSideBoundedVector701 =             0.666666666666667*crRightHandSideBoundedVector382;
const double crRightHandSideBoundedVector702 =             0.666666666666667*crRightHandSideBoundedVector383;
const double crRightHandSideBoundedVector703 =             crRightHandSideBoundedVector278*crRightHandSideBoundedVector431;
const double crRightHandSideBoundedVector704 =             crRightHandSideBoundedVector583*crRightHandSideBoundedVector75;
const double crRightHandSideBoundedVector705 =             crRightHandSideBoundedVector589*crRightHandSideBoundedVector75;
const double crRightHandSideBoundedVector706 =             0.166666666666667*crRightHandSideBoundedVector115*crRightHandSideBoundedVector595;
const double crRightHandSideBoundedVector707 =             crRightHandSideBoundedVector263 - crRightHandSideBoundedVector27*crRightHandSideBoundedVector653 + crRightHandSideBoundedVector27*crRightHandSideBoundedVector655 + crRightHandSideBoundedVector280 - crRightHandSideBoundedVector593*crRightHandSideBoundedVector65 + crRightHandSideBoundedVector64*crRightHandSideBoundedVector706 + crRightHandSideBoundedVector646 - crRightHandSideBoundedVector654 + crRightHandSideBoundedVector697 - crRightHandSideBoundedVector699;
const double crRightHandSideBoundedVector708 =             0.166666666666667*crRightHandSideBoundedVector115*crRightHandSideBoundedVector600;
const double crRightHandSideBoundedVector709 =             -crRightHandSideBoundedVector194*crRightHandSideBoundedVector210 - crRightHandSideBoundedVector27*crRightHandSideBoundedVector688 + crRightHandSideBoundedVector27*crRightHandSideBoundedVector690 + crRightHandSideBoundedVector447 + crRightHandSideBoundedVector459 - crRightHandSideBoundedVector598*crRightHandSideBoundedVector65 + crRightHandSideBoundedVector61*crRightHandSideBoundedVector708 + crRightHandSideBoundedVector681 - crRightHandSideBoundedVector689 + crRightHandSideBoundedVector698;
const double crRightHandSideBoundedVector710 =             crRightHandSideBoundedVector326 + crRightHandSideBoundedVector327 + crRightHandSideBoundedVector328 + crRightHandSideBoundedVector484 + crRightHandSideBoundedVector485 + crRightHandSideBoundedVector486 + crRightHandSideBoundedVector673 + crRightHandSideBoundedVector694;
const double crRightHandSideBoundedVector711 =             0.111111111111111*r[0];
const double crRightHandSideBoundedVector712 =             crRightHandSideBoundedVector48*crRightHandSideBoundedVector708 + crRightHandSideBoundedVector59*crRightHandSideBoundedVector706 + crRightHandSideBoundedVector609 + crRightHandSideBoundedVector621 - 0.166666666666667*crRightHandSideBoundedVector630 - 0.166666666666667*crRightHandSideBoundedVector631 + 0.333333333333333*crRightHandSideBoundedVector632 + 0.333333333333333*crRightHandSideBoundedVector633 - 0.333333333333333*crRightHandSideBoundedVector635 - 0.333333333333333*crRightHandSideBoundedVector636 + crRightHandSideBoundedVector711;
const double crRightHandSideBoundedVector713 =             1.0*DN_DX_2_0*h;
const double crRightHandSideBoundedVector714 =             1.0*DN_DX_2_1*h;
const double crRightHandSideBoundedVector715 =             crRightHandSideBoundedVector211*crRightHandSideBoundedVector23;
const double crRightHandSideBoundedVector716 =             DN_DX_2_1*crRightHandSideBoundedVector215*mu;
const double crRightHandSideBoundedVector717 =             crRightHandSideBoundedVector279 + crRightHandSideBoundedVector646 + 0.444444444444444*f_ext(2,0);
const double crRightHandSideBoundedVector718 =             0.666666666666667*crRightHandSideBoundedVector265;
const double crRightHandSideBoundedVector719 =             DN_DX_2_1*crRightHandSideBoundedVector18;
const double crRightHandSideBoundedVector720 =             crRightHandSideBoundedVector206*crRightHandSideBoundedVector261;
const double crRightHandSideBoundedVector721 =             1.33333333333333*crRightHandSideBoundedVector271;
const double crRightHandSideBoundedVector722 =             DN_DX_2_1*crRightHandSideBoundedVector78;
const double crRightHandSideBoundedVector723 =             DN_DX_2_1*crRightHandSideBoundedVector117;
const double crRightHandSideBoundedVector724 =             DN_DX_2_0*crRightHandSideBoundedVector23;
const double crRightHandSideBoundedVector725 =             DN_DX_2_0*crRightHandSideBoundedVector82;
const double crRightHandSideBoundedVector726 =             DN_DX_2_0*crRightHandSideBoundedVector120;
const double crRightHandSideBoundedVector727 =             DN_DX_2_1*crRightHandSideBoundedVector23;
const double crRightHandSideBoundedVector728 =             DN_DX_2_0*crRightHandSideBoundedVector18;
const double crRightHandSideBoundedVector729 =             -crRightHandSideBoundedVector13*crRightHandSideBoundedVector715;
const double crRightHandSideBoundedVector730 =             DN_DX_2_1*crRightHandSideBoundedVector82;
const double crRightHandSideBoundedVector731 =             DN_DX_2_0*crRightHandSideBoundedVector78;
const double crRightHandSideBoundedVector732 =             DN_DX_2_1*crRightHandSideBoundedVector120;
const double crRightHandSideBoundedVector733 =             DN_DX_2_0*crRightHandSideBoundedVector117;
const double crRightHandSideBoundedVector734 =             DN_DX_2_0*crRightHandSideBoundedVector303*h;
const double crRightHandSideBoundedVector735 =             crRightHandSideBoundedVector18*crRightHandSideBoundedVector214;
const double crRightHandSideBoundedVector736 =             DN_DX_2_0*crRightHandSideBoundedVector215*mu;
const double crRightHandSideBoundedVector737 =             crRightHandSideBoundedVector458 + crRightHandSideBoundedVector681 + 0.444444444444444*f_ext(2,1);
const double crRightHandSideBoundedVector738 =             0.666666666666667*crRightHandSideBoundedVector450;
const double crRightHandSideBoundedVector739 =             crRightHandSideBoundedVector261*crRightHandSideBoundedVector728;
const double crRightHandSideBoundedVector740 =             crRightHandSideBoundedVector269*crRightHandSideBoundedVector431;
const double crRightHandSideBoundedVector741 =             1.33333333333333*crRightHandSideBoundedVector453;
const double crRightHandSideBoundedVector742 =             crRightHandSideBoundedVector278*crRightHandSideBoundedVector731;
const double crRightHandSideBoundedVector743 =             crRightHandSideBoundedVector297*crRightHandSideBoundedVector733;
const double crRightHandSideBoundedVector744 =             -crRightHandSideBoundedVector13*crRightHandSideBoundedVector735;
const double crRightHandSideBoundedVector745 =             DN_DX_2_1*crRightHandSideBoundedVector303*h;
const double crRightHandSideBoundedVector746 =             0.666666666666667*crRightHandSideBoundedVector361;
const double crRightHandSideBoundedVector747 =             0.666666666666667*crRightHandSideBoundedVector362;
const double crRightHandSideBoundedVector748 =             crRightHandSideBoundedVector261*crRightHandSideBoundedVector431;
const double crRightHandSideBoundedVector749 =             crRightHandSideBoundedVector14*crRightHandSideBoundedVector571;
const double crRightHandSideBoundedVector750 =             crRightHandSideBoundedVector14*crRightHandSideBoundedVector577;
            rRightHandSideBoundedVector[0]=-1.0*crRightHandSideBoundedVector112*crRightHandSideBoundedVector7 - 1.0*crRightHandSideBoundedVector146*crRightHandSideBoundedVector7 - 1.0*crRightHandSideBoundedVector147*crRightHandSideBoundedVector165 - 1.0*crRightHandSideBoundedVector147*crRightHandSideBoundedVector178 - 1.0*crRightHandSideBoundedVector147*crRightHandSideBoundedVector189 - 1.0*crRightHandSideBoundedVector6 - 1.0*crRightHandSideBoundedVector7*crRightHandSideBoundedVector71;
            rRightHandSideBoundedVector[1]=-DN_DX_0_0*crRightHandSideBoundedVector242 - DN_DX_0_0*crRightHandSideBoundedVector249 - DN_DX_0_0*crRightHandSideBoundedVector255 + 0.666666666666667*crRightHandSideBoundedVector132 + crRightHandSideBoundedVector133*crRightHandSideBoundedVector208 - 0.666666666666667*crRightHandSideBoundedVector134 + crRightHandSideBoundedVector135*crRightHandSideBoundedVector206 - 0.666666666666667*crRightHandSideBoundedVector136 - crRightHandSideBoundedVector143*crRightHandSideBoundedVector214 + crRightHandSideBoundedVector190 + crRightHandSideBoundedVector191 + crRightHandSideBoundedVector192 + crRightHandSideBoundedVector195 + crRightHandSideBoundedVector199 + crRightHandSideBoundedVector200 + crRightHandSideBoundedVector201 + crRightHandSideBoundedVector202 + crRightHandSideBoundedVector203 + crRightHandSideBoundedVector204 + crRightHandSideBoundedVector205 + crRightHandSideBoundedVector207 + crRightHandSideBoundedVector209 + crRightHandSideBoundedVector210*crRightHandSideBoundedVector212 + crRightHandSideBoundedVector213 - crRightHandSideBoundedVector216*crRightHandSideBoundedVector223 - crRightHandSideBoundedVector216*crRightHandSideBoundedVector225 - crRightHandSideBoundedVector216*crRightHandSideBoundedVector227 - crRightHandSideBoundedVector259*(DN_DX_0_0*crRightHandSideBoundedVector68 - crRightHandSideBoundedVector260*crRightHandSideBoundedVector261 + crRightHandSideBoundedVector274) - crRightHandSideBoundedVector276*(DN_DX_0_0*crRightHandSideBoundedVector109 - crRightHandSideBoundedVector277*crRightHandSideBoundedVector278 + crRightHandSideBoundedVector290) - crRightHandSideBoundedVector291*(DN_DX_0_0*crRightHandSideBoundedVector143 + crRightHandSideBoundedVector208*crRightHandSideBoundedVector210 + crRightHandSideBoundedVector292 - crRightHandSideBoundedVector294 - 0.666666666666667*crRightHandSideBoundedVector295 - crRightHandSideBoundedVector296*crRightHandSideBoundedVector297 + crRightHandSideBoundedVector298 + crRightHandSideBoundedVector301 - 1.33333333333333*crRightHandSideBoundedVector302) + crRightHandSideBoundedVector313*(-crRightHandSideBoundedVector260 + crRightHandSideBoundedVector303*crRightHandSideBoundedVector304 + crRightHandSideBoundedVector312) + crRightHandSideBoundedVector320*(-crRightHandSideBoundedVector277 + crRightHandSideBoundedVector303*crRightHandSideBoundedVector314 + crRightHandSideBoundedVector319) + crRightHandSideBoundedVector323*(-crRightHandSideBoundedVector114*crRightHandSideBoundedVector120*crRightHandSideBoundedVector214*crRightHandSideBoundedVector27 + crRightHandSideBoundedVector133*crRightHandSideBoundedVector211 - crRightHandSideBoundedVector296 + crRightHandSideBoundedVector303*crRightHandSideBoundedVector322 + crRightHandSideBoundedVector321) - crRightHandSideBoundedVector334*(-crRightHandSideBoundedVector196*crRightHandSideBoundedVector325 + crRightHandSideBoundedVector324 + crRightHandSideBoundedVector333) - crRightHandSideBoundedVector341*(-crRightHandSideBoundedVector196*crRightHandSideBoundedVector336 + crRightHandSideBoundedVector335 + crRightHandSideBoundedVector340) - crRightHandSideBoundedVector349*(crRightHandSideBoundedVector114*crRightHandSideBoundedVector117*crRightHandSideBoundedVector196*crRightHandSideBoundedVector214 - crRightHandSideBoundedVector196*crRightHandSideBoundedVector347 + crRightHandSideBoundedVector345 + crRightHandSideBoundedVector346 + crRightHandSideBoundedVector348) - crRightHandSideBoundedVector350*crRightHandSideBoundedVector375 - crRightHandSideBoundedVector350*crRightHandSideBoundedVector396 - crRightHandSideBoundedVector350*crRightHandSideBoundedVector415 - 1.0*crRightHandSideBoundedVector41;
            rRightHandSideBoundedVector[2]=-DN_DX_0_1*crRightHandSideBoundedVector442 - DN_DX_0_1*crRightHandSideBoundedVector443 - DN_DX_0_1*crRightHandSideBoundedVector444 + crRightHandSideBoundedVector133*crRightHandSideBoundedVector431 + crRightHandSideBoundedVector135*crRightHandSideBoundedVector433 - 1.0*crRightHandSideBoundedVector149 + 0.666666666666667*crRightHandSideBoundedVector180 - 0.666666666666667*crRightHandSideBoundedVector181 - 0.666666666666667*crRightHandSideBoundedVector182 - crRightHandSideBoundedVector186*crRightHandSideBoundedVector211 - crRightHandSideBoundedVector223*crRightHandSideBoundedVector437 - crRightHandSideBoundedVector225*crRightHandSideBoundedVector437 - crRightHandSideBoundedVector227*crRightHandSideBoundedVector437 - crRightHandSideBoundedVector259*(DN_DX_0_1*crRightHandSideBoundedVector162 - crRightHandSideBoundedVector445 + crRightHandSideBoundedVector456) - crRightHandSideBoundedVector276*(DN_DX_0_1*crRightHandSideBoundedVector175 - crRightHandSideBoundedVector457 + crRightHandSideBoundedVector467) - crRightHandSideBoundedVector291*(DN_DX_0_1*crRightHandSideBoundedVector186 + crRightHandSideBoundedVector297*crRightHandSideBoundedVector433 + crRightHandSideBoundedVector468 - 0.666666666666667*crRightHandSideBoundedVector469 - crRightHandSideBoundedVector471 - crRightHandSideBoundedVector472 + crRightHandSideBoundedVector473 + crRightHandSideBoundedVector475 - 1.33333333333333*crRightHandSideBoundedVector476) + crRightHandSideBoundedVector297*crRightHandSideBoundedVector435 - crRightHandSideBoundedVector313*(-crRightHandSideBoundedVector196*crRightHandSideBoundedVector324 + crRightHandSideBoundedVector325 + crRightHandSideBoundedVector489) - crRightHandSideBoundedVector320*(-crRightHandSideBoundedVector196*crRightHandSideBoundedVector335 + crRightHandSideBoundedVector336 + crRightHandSideBoundedVector491) - crRightHandSideBoundedVector323*(crRightHandSideBoundedVector114*crRightHandSideBoundedVector120*crRightHandSideBoundedVector196*crRightHandSideBoundedVector211 - crRightHandSideBoundedVector196*crRightHandSideBoundedVector346 + crRightHandSideBoundedVector347 + crRightHandSideBoundedVector495 + crRightHandSideBoundedVector496) - crRightHandSideBoundedVector334*(-crRightHandSideBoundedVector260*crRightHandSideBoundedVector303 + crRightHandSideBoundedVector304 + crRightHandSideBoundedVector481) - crRightHandSideBoundedVector341*(-crRightHandSideBoundedVector277*crRightHandSideBoundedVector303 + crRightHandSideBoundedVector314 + crRightHandSideBoundedVector482) - crRightHandSideBoundedVector349*(crRightHandSideBoundedVector114*crRightHandSideBoundedVector117*crRightHandSideBoundedVector211*crRightHandSideBoundedVector27 - crRightHandSideBoundedVector135*crRightHandSideBoundedVector214 - crRightHandSideBoundedVector296*crRightHandSideBoundedVector303 + crRightHandSideBoundedVector322 + crRightHandSideBoundedVector483) - crRightHandSideBoundedVector375*crRightHandSideBoundedVector497 - crRightHandSideBoundedVector396*crRightHandSideBoundedVector497 - crRightHandSideBoundedVector415*crRightHandSideBoundedVector497 + crRightHandSideBoundedVector416 + crRightHandSideBoundedVector417 + crRightHandSideBoundedVector418 + crRightHandSideBoundedVector421 + crRightHandSideBoundedVector424 + crRightHandSideBoundedVector425 + crRightHandSideBoundedVector426 + crRightHandSideBoundedVector427 + crRightHandSideBoundedVector428 + crRightHandSideBoundedVector429 + crRightHandSideBoundedVector430 + crRightHandSideBoundedVector432 + crRightHandSideBoundedVector434 + crRightHandSideBoundedVector436;
            rRightHandSideBoundedVector[3]=-DN_DX_0_0*crRightHandSideBoundedVector551 - DN_DX_0_0*crRightHandSideBoundedVector556 - DN_DX_0_0*crRightHandSideBoundedVector561 - DN_DX_0_1*crRightHandSideBoundedVector564 - DN_DX_0_1*crRightHandSideBoundedVector565 - DN_DX_0_1*crRightHandSideBoundedVector566 + crRightHandSideBoundedVector117*crRightHandSideBoundedVector298 - crRightHandSideBoundedVector117*crRightHandSideBoundedVector523 + crRightHandSideBoundedVector117*crRightHandSideBoundedVector528 - crRightHandSideBoundedVector120*crRightHandSideBoundedVector524 - crRightHandSideBoundedVector212*crRightHandSideBoundedVector541 - crRightHandSideBoundedVector259*(crRightHandSideBoundedVector324*crRightHandSideBoundedVector607 + crRightHandSideBoundedVector325*crRightHandSideBoundedVector607 + crRightHandSideBoundedVector618) - crRightHandSideBoundedVector276*(crRightHandSideBoundedVector335*crRightHandSideBoundedVector619 + crRightHandSideBoundedVector336*crRightHandSideBoundedVector619 + crRightHandSideBoundedVector629) - crRightHandSideBoundedVector291*(crRightHandSideBoundedVector346*crRightHandSideBoundedVector541 + crRightHandSideBoundedVector347*crRightHandSideBoundedVector541 + crRightHandSideBoundedVector596*crRightHandSideBoundedVector637 + crRightHandSideBoundedVector601*crRightHandSideBoundedVector638 + crRightHandSideBoundedVector610 + crRightHandSideBoundedVector620 - 0.666666666666667*crRightHandSideBoundedVector630 - 0.666666666666667*crRightHandSideBoundedVector631 + 1.33333333333333*crRightHandSideBoundedVector632 + 1.33333333333333*crRightHandSideBoundedVector633 - 1.33333333333333*crRightHandSideBoundedVector635 - 1.33333333333333*crRightHandSideBoundedVector636 + 0.444444444444444*r[0]) + 0.666666666666667*crRightHandSideBoundedVector398 + 0.666666666666667*crRightHandSideBoundedVector399 + 0.666666666666667*crRightHandSideBoundedVector400 - crRightHandSideBoundedVector435*crRightHandSideBoundedVector541 + crRightHandSideBoundedVector498 + crRightHandSideBoundedVector499 + crRightHandSideBoundedVector500 + crRightHandSideBoundedVector502 + crRightHandSideBoundedVector504 + crRightHandSideBoundedVector505 + crRightHandSideBoundedVector507 + crRightHandSideBoundedVector510 + crRightHandSideBoundedVector513 + crRightHandSideBoundedVector514 + crRightHandSideBoundedVector515 + crRightHandSideBoundedVector516 + crRightHandSideBoundedVector517 + crRightHandSideBoundedVector518 + crRightHandSideBoundedVector520 + crRightHandSideBoundedVector522 + crRightHandSideBoundedVector525 + crRightHandSideBoundedVector527 + crRightHandSideBoundedVector531 + crRightHandSideBoundedVector534 + 0.666666666666667*crRightHandSideBoundedVector536 + 0.666666666666667*crRightHandSideBoundedVector538 + crRightHandSideBoundedVector539 + crRightHandSideBoundedVector540 - crRightHandSideBoundedVector567*(-DN_DX_0_0*crRightHandSideBoundedVector568 - DN_DX_0_1*crRightHandSideBoundedVector18*crRightHandSideBoundedVector261*crRightHandSideBoundedVector303 + crRightHandSideBoundedVector573) - crRightHandSideBoundedVector574*(-DN_DX_0_1*crRightHandSideBoundedVector575 - crRightHandSideBoundedVector303*crRightHandSideBoundedVector445 + crRightHandSideBoundedVector579) - crRightHandSideBoundedVector580*(-DN_DX_0_0*crRightHandSideBoundedVector581 - DN_DX_0_1*crRightHandSideBoundedVector278*crRightHandSideBoundedVector303*crRightHandSideBoundedVector78 + crRightHandSideBoundedVector585) - crRightHandSideBoundedVector586*(-DN_DX_0_1*crRightHandSideBoundedVector587 - crRightHandSideBoundedVector303*crRightHandSideBoundedVector457 + crRightHandSideBoundedVector591) - crRightHandSideBoundedVector592*(-DN_DX_0_0*crRightHandSideBoundedVector594 - DN_DX_0_1*crRightHandSideBoundedVector117*crRightHandSideBoundedVector297*crRightHandSideBoundedVector303 + crRightHandSideBoundedVector214*crRightHandSideBoundedVector596 - crRightHandSideBoundedVector27*crRightHandSideBoundedVector294 + crRightHandSideBoundedVector27*crRightHandSideBoundedVector301 + crRightHandSideBoundedVector292 - crRightHandSideBoundedVector298 + crRightHandSideBoundedVector523 - crRightHandSideBoundedVector528 - crRightHandSideBoundedVector570*crRightHandSideBoundedVector593) - crRightHandSideBoundedVector597*(-DN_DX_0_1*crRightHandSideBoundedVector599 - crRightHandSideBoundedVector206*crRightHandSideBoundedVector210 + crRightHandSideBoundedVector211*crRightHandSideBoundedVector601 - crRightHandSideBoundedVector27*crRightHandSideBoundedVector471 + crRightHandSideBoundedVector27*crRightHandSideBoundedVector475 - crRightHandSideBoundedVector303*crRightHandSideBoundedVector472 + crRightHandSideBoundedVector468 - crRightHandSideBoundedVector473 + crRightHandSideBoundedVector524 - crRightHandSideBoundedVector570*crRightHandSideBoundedVector598) - crRightHandSideBoundedVector603*(crRightHandSideBoundedVector324 + crRightHandSideBoundedVector325 + crRightHandSideBoundedVector602) - crRightHandSideBoundedVector605*(crRightHandSideBoundedVector335 + crRightHandSideBoundedVector336 + crRightHandSideBoundedVector604) - crRightHandSideBoundedVector606*(crRightHandSideBoundedVector342 + crRightHandSideBoundedVector343 + crRightHandSideBoundedVector344 + crRightHandSideBoundedVector346 + crRightHandSideBoundedVector347 + crRightHandSideBoundedVector348 + crRightHandSideBoundedVector492 + crRightHandSideBoundedVector493 + crRightHandSideBoundedVector494 + crRightHandSideBoundedVector496);
            rRightHandSideBoundedVector[4]=-crRightHandSideBoundedVector112*crRightHandSideBoundedVector640 - crRightHandSideBoundedVector146*crRightHandSideBoundedVector640 - crRightHandSideBoundedVector165*crRightHandSideBoundedVector641 - crRightHandSideBoundedVector178*crRightHandSideBoundedVector641 - crRightHandSideBoundedVector189*crRightHandSideBoundedVector641 - crRightHandSideBoundedVector639 - crRightHandSideBoundedVector640*crRightHandSideBoundedVector71;
            rRightHandSideBoundedVector[5]=-DN_DX_1_0*crRightHandSideBoundedVector242 - DN_DX_1_0*crRightHandSideBoundedVector249 - DN_DX_1_0*crRightHandSideBoundedVector255 - 0.666666666666667*crRightHandSideBoundedVector100 + crRightHandSideBoundedVector101*crRightHandSideBoundedVector206 - 0.666666666666667*crRightHandSideBoundedVector102 - crRightHandSideBoundedVector109*crRightHandSideBoundedVector214 + crRightHandSideBoundedVector202 + crRightHandSideBoundedVector203 + crRightHandSideBoundedVector204 + crRightHandSideBoundedVector205 + crRightHandSideBoundedVector207 + crRightHandSideBoundedVector208*crRightHandSideBoundedVector99 + crRightHandSideBoundedVector209 + crRightHandSideBoundedVector213 - crRightHandSideBoundedVector223*crRightHandSideBoundedVector644 - crRightHandSideBoundedVector225*crRightHandSideBoundedVector644 - crRightHandSideBoundedVector227*crRightHandSideBoundedVector644 - crRightHandSideBoundedVector259*(DN_DX_1_0*crRightHandSideBoundedVector68 - crRightHandSideBoundedVector261*crRightHandSideBoundedVector645 + crRightHandSideBoundedVector274) - crRightHandSideBoundedVector276*(DN_DX_1_0*crRightHandSideBoundedVector109 + crRightHandSideBoundedVector208*crRightHandSideBoundedVector285 - crRightHandSideBoundedVector278*crRightHandSideBoundedVector649 - 0.666666666666667*crRightHandSideBoundedVector283 - 1.33333333333333*crRightHandSideBoundedVector289 + crRightHandSideBoundedVector647 - crRightHandSideBoundedVector648 + crRightHandSideBoundedVector650 + crRightHandSideBoundedVector651) + crRightHandSideBoundedVector285*crRightHandSideBoundedVector643 - crRightHandSideBoundedVector291*(DN_DX_1_0*crRightHandSideBoundedVector143 - crRightHandSideBoundedVector297*crRightHandSideBoundedVector652 + crRightHandSideBoundedVector656) + crRightHandSideBoundedVector313*(crRightHandSideBoundedVector303*crRightHandSideBoundedVector657 + crRightHandSideBoundedVector312 - crRightHandSideBoundedVector645) + crRightHandSideBoundedVector320*(crRightHandSideBoundedVector211*crRightHandSideBoundedVector99 - crRightHandSideBoundedVector214*crRightHandSideBoundedVector27*crRightHandSideBoundedVector74*crRightHandSideBoundedVector82 + crRightHandSideBoundedVector303*crRightHandSideBoundedVector658 + crRightHandSideBoundedVector321 - crRightHandSideBoundedVector649) + crRightHandSideBoundedVector323*(crRightHandSideBoundedVector303*crRightHandSideBoundedVector659 - crRightHandSideBoundedVector652 + crRightHandSideBoundedVector664) - crRightHandSideBoundedVector334*(-crRightHandSideBoundedVector196*crRightHandSideBoundedVector666 + crRightHandSideBoundedVector333 + crRightHandSideBoundedVector665) - crRightHandSideBoundedVector341*(crRightHandSideBoundedVector196*crRightHandSideBoundedVector214*crRightHandSideBoundedVector74*crRightHandSideBoundedVector78 - crRightHandSideBoundedVector196*crRightHandSideBoundedVector668 + crRightHandSideBoundedVector345 + crRightHandSideBoundedVector667 + crRightHandSideBoundedVector669) - crRightHandSideBoundedVector349*(-crRightHandSideBoundedVector196*crRightHandSideBoundedVector671 + crRightHandSideBoundedVector670 + crRightHandSideBoundedVector675) - crRightHandSideBoundedVector375*crRightHandSideBoundedVector676 - crRightHandSideBoundedVector396*crRightHandSideBoundedVector676 - crRightHandSideBoundedVector415*crRightHandSideBoundedVector676 + crRightHandSideBoundedVector642 + 0.666666666666667*crRightHandSideBoundedVector98;
            rRightHandSideBoundedVector[6]=-DN_DX_1_1*crRightHandSideBoundedVector442 - DN_DX_1_1*crRightHandSideBoundedVector443 - DN_DX_1_1*crRightHandSideBoundedVector444 + crRightHandSideBoundedVector101*crRightHandSideBoundedVector433 + 0.666666666666667*crRightHandSideBoundedVector169 - 0.666666666666667*crRightHandSideBoundedVector170 - 0.666666666666667*crRightHandSideBoundedVector171 - crRightHandSideBoundedVector175*crRightHandSideBoundedVector211 - crRightHandSideBoundedVector223*crRightHandSideBoundedVector679 - crRightHandSideBoundedVector225*crRightHandSideBoundedVector679 - crRightHandSideBoundedVector227*crRightHandSideBoundedVector679 - crRightHandSideBoundedVector259*(DN_DX_1_1*crRightHandSideBoundedVector162 + crRightHandSideBoundedVector456 - crRightHandSideBoundedVector680) - crRightHandSideBoundedVector276*(DN_DX_1_1*crRightHandSideBoundedVector175 + crRightHandSideBoundedVector278*crRightHandSideBoundedVector433 - 0.666666666666667*crRightHandSideBoundedVector460 - 1.33333333333333*crRightHandSideBoundedVector466 + crRightHandSideBoundedVector682 - crRightHandSideBoundedVector683 - crRightHandSideBoundedVector684 + crRightHandSideBoundedVector685 + crRightHandSideBoundedVector686) + crRightHandSideBoundedVector278*crRightHandSideBoundedVector678 - crRightHandSideBoundedVector291*(DN_DX_1_1*crRightHandSideBoundedVector186 - crRightHandSideBoundedVector687 + crRightHandSideBoundedVector691) - crRightHandSideBoundedVector313*(-crRightHandSideBoundedVector196*crRightHandSideBoundedVector665 + crRightHandSideBoundedVector489 + crRightHandSideBoundedVector666) - crRightHandSideBoundedVector320*(crRightHandSideBoundedVector196*crRightHandSideBoundedVector211*crRightHandSideBoundedVector74*crRightHandSideBoundedVector82 - crRightHandSideBoundedVector196*crRightHandSideBoundedVector667 + crRightHandSideBoundedVector495 + crRightHandSideBoundedVector668 + crRightHandSideBoundedVector693) - crRightHandSideBoundedVector323*(-crRightHandSideBoundedVector196*crRightHandSideBoundedVector670 + crRightHandSideBoundedVector671 + crRightHandSideBoundedVector695) - crRightHandSideBoundedVector334*(-crRightHandSideBoundedVector303*crRightHandSideBoundedVector645 + crRightHandSideBoundedVector481 + crRightHandSideBoundedVector657) - crRightHandSideBoundedVector341*(-crRightHandSideBoundedVector101*crRightHandSideBoundedVector214 + crRightHandSideBoundedVector211*crRightHandSideBoundedVector27*crRightHandSideBoundedVector74*crRightHandSideBoundedVector78 - crRightHandSideBoundedVector303*crRightHandSideBoundedVector649 + crRightHandSideBoundedVector483 + crRightHandSideBoundedVector658) - crRightHandSideBoundedVector349*(-crRightHandSideBoundedVector303*crRightHandSideBoundedVector652 + crRightHandSideBoundedVector659 + crRightHandSideBoundedVector692) - crRightHandSideBoundedVector375*crRightHandSideBoundedVector696 - crRightHandSideBoundedVector396*crRightHandSideBoundedVector696 - crRightHandSideBoundedVector415*crRightHandSideBoundedVector696 + crRightHandSideBoundedVector427 + crRightHandSideBoundedVector428 + crRightHandSideBoundedVector429 + crRightHandSideBoundedVector430 + crRightHandSideBoundedVector431*crRightHandSideBoundedVector99 + crRightHandSideBoundedVector432 + crRightHandSideBoundedVector434 + crRightHandSideBoundedVector436 + crRightHandSideBoundedVector677;
            rRightHandSideBoundedVector[7]=-DN_DX_1_0*crRightHandSideBoundedVector551 - DN_DX_1_0*crRightHandSideBoundedVector556 - DN_DX_1_0*crRightHandSideBoundedVector561 - DN_DX_1_1*crRightHandSideBoundedVector564 - DN_DX_1_1*crRightHandSideBoundedVector565 - DN_DX_1_1*crRightHandSideBoundedVector566 - crRightHandSideBoundedVector259*(crRightHandSideBoundedVector607*crRightHandSideBoundedVector665 + crRightHandSideBoundedVector607*crRightHandSideBoundedVector666 + crRightHandSideBoundedVector618) - crRightHandSideBoundedVector276*(crRightHandSideBoundedVector610 + crRightHandSideBoundedVector619*crRightHandSideBoundedVector667 + crRightHandSideBoundedVector619*crRightHandSideBoundedVector668 - 0.666666666666667*crRightHandSideBoundedVector622 - 0.666666666666667*crRightHandSideBoundedVector623 + 1.33333333333333*crRightHandSideBoundedVector624 + 1.33333333333333*crRightHandSideBoundedVector625 - 1.33333333333333*crRightHandSideBoundedVector627 - 1.33333333333333*crRightHandSideBoundedVector628 + crRightHandSideBoundedVector637*crRightHandSideBoundedVector704 + crRightHandSideBoundedVector638*crRightHandSideBoundedVector705 + crRightHandSideBoundedVector711 + 0.444444444444444*r[1]) - crRightHandSideBoundedVector291*(crRightHandSideBoundedVector541*crRightHandSideBoundedVector670 + crRightHandSideBoundedVector541*crRightHandSideBoundedVector671 + crRightHandSideBoundedVector712) + 0.666666666666667*crRightHandSideBoundedVector379 + 0.666666666666667*crRightHandSideBoundedVector380 + 0.666666666666667*crRightHandSideBoundedVector381 + 0.666666666666667*crRightHandSideBoundedVector509 + 0.666666666666667*crRightHandSideBoundedVector512 + crRightHandSideBoundedVector516 + crRightHandSideBoundedVector517 + crRightHandSideBoundedVector518 + crRightHandSideBoundedVector520 + crRightHandSideBoundedVector522 + crRightHandSideBoundedVector525 + crRightHandSideBoundedVector527 + crRightHandSideBoundedVector531 + crRightHandSideBoundedVector534 + crRightHandSideBoundedVector539 + crRightHandSideBoundedVector540 - crRightHandSideBoundedVector567*(-DN_DX_1_0*crRightHandSideBoundedVector568 - DN_DX_1_1*crRightHandSideBoundedVector18*crRightHandSideBoundedVector261*crRightHandSideBoundedVector303 + crRightHandSideBoundedVector573) - crRightHandSideBoundedVector574*(-DN_DX_1_1*crRightHandSideBoundedVector575 - crRightHandSideBoundedVector303*crRightHandSideBoundedVector680 + crRightHandSideBoundedVector579) - crRightHandSideBoundedVector580*(-DN_DX_1_0*crRightHandSideBoundedVector581 - DN_DX_1_1*crRightHandSideBoundedVector278*crRightHandSideBoundedVector303*crRightHandSideBoundedVector78 + crRightHandSideBoundedVector214*crRightHandSideBoundedVector704 - crRightHandSideBoundedVector27*crRightHandSideBoundedVector648 + crRightHandSideBoundedVector27*crRightHandSideBoundedVector651 - crRightHandSideBoundedVector570*crRightHandSideBoundedVector582 + crRightHandSideBoundedVector647 - crRightHandSideBoundedVector650 + crRightHandSideBoundedVector701 - crRightHandSideBoundedVector703) - crRightHandSideBoundedVector586*(-DN_DX_1_1*crRightHandSideBoundedVector587 - crRightHandSideBoundedVector206*crRightHandSideBoundedVector285 + crRightHandSideBoundedVector211*crRightHandSideBoundedVector705 - crRightHandSideBoundedVector27*crRightHandSideBoundedVector683 + crRightHandSideBoundedVector27*crRightHandSideBoundedVector686 - crRightHandSideBoundedVector303*crRightHandSideBoundedVector684 - crRightHandSideBoundedVector570*crRightHandSideBoundedVector588 + crRightHandSideBoundedVector682 - crRightHandSideBoundedVector685 + crRightHandSideBoundedVector702) - crRightHandSideBoundedVector592*(-DN_DX_1_0*crRightHandSideBoundedVector594 - DN_DX_1_1*crRightHandSideBoundedVector117*crRightHandSideBoundedVector297*crRightHandSideBoundedVector303 + crRightHandSideBoundedVector707) - crRightHandSideBoundedVector597*(-DN_DX_1_1*crRightHandSideBoundedVector599 - crRightHandSideBoundedVector303*crRightHandSideBoundedVector687 + crRightHandSideBoundedVector709) - crRightHandSideBoundedVector603*(crRightHandSideBoundedVector602 + crRightHandSideBoundedVector665 + crRightHandSideBoundedVector666) - crRightHandSideBoundedVector605*(crRightHandSideBoundedVector342 + crRightHandSideBoundedVector343 + crRightHandSideBoundedVector344 + crRightHandSideBoundedVector492 + crRightHandSideBoundedVector493 + crRightHandSideBoundedVector494 + crRightHandSideBoundedVector667 + crRightHandSideBoundedVector668 + crRightHandSideBoundedVector669 + crRightHandSideBoundedVector693) - crRightHandSideBoundedVector606*(crRightHandSideBoundedVector670 + crRightHandSideBoundedVector671 + crRightHandSideBoundedVector710) - crRightHandSideBoundedVector619*crRightHandSideBoundedVector643 - crRightHandSideBoundedVector619*crRightHandSideBoundedVector678 + crRightHandSideBoundedVector650*crRightHandSideBoundedVector78 + crRightHandSideBoundedVector700 - crRightHandSideBoundedVector701*crRightHandSideBoundedVector78 - crRightHandSideBoundedVector702*crRightHandSideBoundedVector82 + crRightHandSideBoundedVector703*crRightHandSideBoundedVector78;
            rRightHandSideBoundedVector[8]=-crRightHandSideBoundedVector112*crRightHandSideBoundedVector713 - crRightHandSideBoundedVector146*crRightHandSideBoundedVector713 - crRightHandSideBoundedVector165*crRightHandSideBoundedVector714 - crRightHandSideBoundedVector178*crRightHandSideBoundedVector714 - crRightHandSideBoundedVector189*crRightHandSideBoundedVector714 - crRightHandSideBoundedVector639 - crRightHandSideBoundedVector71*crRightHandSideBoundedVector713;
            rRightHandSideBoundedVector[9]=-DN_DX_2_0*crRightHandSideBoundedVector242 - DN_DX_2_0*crRightHandSideBoundedVector249 - DN_DX_2_0*crRightHandSideBoundedVector255 + crRightHandSideBoundedVector190 + crRightHandSideBoundedVector191 + crRightHandSideBoundedVector192 + crRightHandSideBoundedVector195 + crRightHandSideBoundedVector199 + crRightHandSideBoundedVector200 + crRightHandSideBoundedVector201 + crRightHandSideBoundedVector206*crRightHandSideBoundedVector52 + crRightHandSideBoundedVector208*crRightHandSideBoundedVector49 - crRightHandSideBoundedVector214*crRightHandSideBoundedVector68 - crRightHandSideBoundedVector223*crRightHandSideBoundedVector716 - crRightHandSideBoundedVector225*crRightHandSideBoundedVector716 - crRightHandSideBoundedVector227*crRightHandSideBoundedVector716 - crRightHandSideBoundedVector259*(DN_DX_2_0*crRightHandSideBoundedVector68 + crRightHandSideBoundedVector208*crRightHandSideBoundedVector269 - crRightHandSideBoundedVector261*crRightHandSideBoundedVector719 - 0.666666666666667*crRightHandSideBoundedVector267 - 1.33333333333333*crRightHandSideBoundedVector273 + crRightHandSideBoundedVector717 - crRightHandSideBoundedVector718 + crRightHandSideBoundedVector720 + crRightHandSideBoundedVector721) + crRightHandSideBoundedVector269*crRightHandSideBoundedVector715 - crRightHandSideBoundedVector276*(DN_DX_2_0*crRightHandSideBoundedVector109 - crRightHandSideBoundedVector278*crRightHandSideBoundedVector722 + crRightHandSideBoundedVector290) - crRightHandSideBoundedVector291*(DN_DX_2_0*crRightHandSideBoundedVector143 - crRightHandSideBoundedVector297*crRightHandSideBoundedVector723 + crRightHandSideBoundedVector656) + crRightHandSideBoundedVector313*(-crRightHandSideBoundedVector13*crRightHandSideBoundedVector214*crRightHandSideBoundedVector23*crRightHandSideBoundedVector27 + crRightHandSideBoundedVector211*crRightHandSideBoundedVector49 + crRightHandSideBoundedVector303*crRightHandSideBoundedVector724 + crRightHandSideBoundedVector321 - crRightHandSideBoundedVector719) + crRightHandSideBoundedVector320*(crRightHandSideBoundedVector303*crRightHandSideBoundedVector725 + crRightHandSideBoundedVector319 - crRightHandSideBoundedVector722) + crRightHandSideBoundedVector323*(crRightHandSideBoundedVector303*crRightHandSideBoundedVector726 + crRightHandSideBoundedVector664 - crRightHandSideBoundedVector723) - crRightHandSideBoundedVector334*(crRightHandSideBoundedVector13*crRightHandSideBoundedVector18*crRightHandSideBoundedVector196*crRightHandSideBoundedVector214 - crRightHandSideBoundedVector196*crRightHandSideBoundedVector728 + crRightHandSideBoundedVector345 + crRightHandSideBoundedVector727 + crRightHandSideBoundedVector729) - crRightHandSideBoundedVector341*(-crRightHandSideBoundedVector196*crRightHandSideBoundedVector731 + crRightHandSideBoundedVector340 + crRightHandSideBoundedVector730) - crRightHandSideBoundedVector349*(-crRightHandSideBoundedVector196*crRightHandSideBoundedVector733 + crRightHandSideBoundedVector675 + crRightHandSideBoundedVector732) - crRightHandSideBoundedVector375*crRightHandSideBoundedVector734 - crRightHandSideBoundedVector396*crRightHandSideBoundedVector734 - crRightHandSideBoundedVector415*crRightHandSideBoundedVector734 + 0.666666666666667*crRightHandSideBoundedVector47 - 0.666666666666667*crRightHandSideBoundedVector50 - 0.666666666666667*crRightHandSideBoundedVector53 + crRightHandSideBoundedVector642;
            rRightHandSideBoundedVector[10]=-DN_DX_2_1*crRightHandSideBoundedVector442 - DN_DX_2_1*crRightHandSideBoundedVector443 - DN_DX_2_1*crRightHandSideBoundedVector444 + 0.666666666666667*crRightHandSideBoundedVector155 - 0.666666666666667*crRightHandSideBoundedVector156 - 0.666666666666667*crRightHandSideBoundedVector157 - crRightHandSideBoundedVector162*crRightHandSideBoundedVector211 - crRightHandSideBoundedVector223*crRightHandSideBoundedVector736 - crRightHandSideBoundedVector225*crRightHandSideBoundedVector736 - crRightHandSideBoundedVector227*crRightHandSideBoundedVector736 - crRightHandSideBoundedVector259*(DN_DX_2_1*crRightHandSideBoundedVector162 + crRightHandSideBoundedVector261*crRightHandSideBoundedVector433 - 0.666666666666667*crRightHandSideBoundedVector449 - 1.33333333333333*crRightHandSideBoundedVector455 + crRightHandSideBoundedVector737 - crRightHandSideBoundedVector738 - crRightHandSideBoundedVector739 + crRightHandSideBoundedVector740 + crRightHandSideBoundedVector741) + crRightHandSideBoundedVector261*crRightHandSideBoundedVector735 - crRightHandSideBoundedVector276*(DN_DX_2_1*crRightHandSideBoundedVector175 + crRightHandSideBoundedVector467 - crRightHandSideBoundedVector742) - crRightHandSideBoundedVector291*(DN_DX_2_1*crRightHandSideBoundedVector186 + crRightHandSideBoundedVector691 - crRightHandSideBoundedVector743) - crRightHandSideBoundedVector313*(crRightHandSideBoundedVector13*crRightHandSideBoundedVector196*crRightHandSideBoundedVector211*crRightHandSideBoundedVector23 - crRightHandSideBoundedVector196*crRightHandSideBoundedVector727 + crRightHandSideBoundedVector495 + crRightHandSideBoundedVector728 + crRightHandSideBoundedVector744) - crRightHandSideBoundedVector320*(-crRightHandSideBoundedVector196*crRightHandSideBoundedVector730 + crRightHandSideBoundedVector491 + crRightHandSideBoundedVector731) - crRightHandSideBoundedVector323*(-crRightHandSideBoundedVector196*crRightHandSideBoundedVector732 + crRightHandSideBoundedVector695 + crRightHandSideBoundedVector733) - crRightHandSideBoundedVector334*(crRightHandSideBoundedVector13*crRightHandSideBoundedVector18*crRightHandSideBoundedVector211*crRightHandSideBoundedVector27 - crRightHandSideBoundedVector214*crRightHandSideBoundedVector52 - crRightHandSideBoundedVector303*crRightHandSideBoundedVector719 + crRightHandSideBoundedVector483 + crRightHandSideBoundedVector724) - crRightHandSideBoundedVector341*(-crRightHandSideBoundedVector303*crRightHandSideBoundedVector722 + crRightHandSideBoundedVector482 + crRightHandSideBoundedVector725) - crRightHandSideBoundedVector349*(-crRightHandSideBoundedVector303*crRightHandSideBoundedVector723 + crRightHandSideBoundedVector692 + crRightHandSideBoundedVector726) - crRightHandSideBoundedVector375*crRightHandSideBoundedVector745 - crRightHandSideBoundedVector396*crRightHandSideBoundedVector745 - crRightHandSideBoundedVector415*crRightHandSideBoundedVector745 + crRightHandSideBoundedVector416 + crRightHandSideBoundedVector417 + crRightHandSideBoundedVector418 + crRightHandSideBoundedVector421 + crRightHandSideBoundedVector424 + crRightHandSideBoundedVector425 + crRightHandSideBoundedVector426 + crRightHandSideBoundedVector431*crRightHandSideBoundedVector49 + crRightHandSideBoundedVector433*crRightHandSideBoundedVector52 + crRightHandSideBoundedVector677;
            rRightHandSideBoundedVector[11]=-DN_DX_2_0*crRightHandSideBoundedVector551 - DN_DX_2_0*crRightHandSideBoundedVector556 - DN_DX_2_0*crRightHandSideBoundedVector561 - DN_DX_2_1*crRightHandSideBoundedVector564 - DN_DX_2_1*crRightHandSideBoundedVector565 - DN_DX_2_1*crRightHandSideBoundedVector566 + crRightHandSideBoundedVector18*crRightHandSideBoundedVector720 - crRightHandSideBoundedVector18*crRightHandSideBoundedVector746 + crRightHandSideBoundedVector18*crRightHandSideBoundedVector748 - crRightHandSideBoundedVector23*crRightHandSideBoundedVector747 - crRightHandSideBoundedVector259*(crRightHandSideBoundedVector607*crRightHandSideBoundedVector727 + crRightHandSideBoundedVector607*crRightHandSideBoundedVector728 - 0.666666666666667*crRightHandSideBoundedVector611 - 0.666666666666667*crRightHandSideBoundedVector612 + 1.33333333333333*crRightHandSideBoundedVector613 + 1.33333333333333*crRightHandSideBoundedVector614 - 1.33333333333333*crRightHandSideBoundedVector616 - 1.33333333333333*crRightHandSideBoundedVector617 + crRightHandSideBoundedVector620 + crRightHandSideBoundedVector637*crRightHandSideBoundedVector749 + crRightHandSideBoundedVector638*crRightHandSideBoundedVector750 + crRightHandSideBoundedVector711 + 0.444444444444444*r[2]) - crRightHandSideBoundedVector276*(crRightHandSideBoundedVector619*crRightHandSideBoundedVector730 + crRightHandSideBoundedVector619*crRightHandSideBoundedVector731 + crRightHandSideBoundedVector629) - crRightHandSideBoundedVector291*(crRightHandSideBoundedVector541*crRightHandSideBoundedVector732 + crRightHandSideBoundedVector541*crRightHandSideBoundedVector733 + crRightHandSideBoundedVector712) + 0.666666666666667*crRightHandSideBoundedVector358 + 0.666666666666667*crRightHandSideBoundedVector359 + 0.666666666666667*crRightHandSideBoundedVector360 + crRightHandSideBoundedVector498 + crRightHandSideBoundedVector499 + crRightHandSideBoundedVector500 + crRightHandSideBoundedVector502 + crRightHandSideBoundedVector504 + crRightHandSideBoundedVector505 + crRightHandSideBoundedVector507 + crRightHandSideBoundedVector510 + crRightHandSideBoundedVector513 + crRightHandSideBoundedVector514 + crRightHandSideBoundedVector515 + 0.666666666666667*crRightHandSideBoundedVector530 + 0.666666666666667*crRightHandSideBoundedVector533 - crRightHandSideBoundedVector567*(-DN_DX_2_0*crRightHandSideBoundedVector568 - DN_DX_2_1*crRightHandSideBoundedVector18*crRightHandSideBoundedVector261*crRightHandSideBoundedVector303 + crRightHandSideBoundedVector214*crRightHandSideBoundedVector749 - crRightHandSideBoundedVector27*crRightHandSideBoundedVector718 + crRightHandSideBoundedVector27*crRightHandSideBoundedVector721 - crRightHandSideBoundedVector569*crRightHandSideBoundedVector570 + crRightHandSideBoundedVector717 - crRightHandSideBoundedVector720 + crRightHandSideBoundedVector746 - crRightHandSideBoundedVector748) - crRightHandSideBoundedVector574*(-DN_DX_2_1*crRightHandSideBoundedVector575 - crRightHandSideBoundedVector206*crRightHandSideBoundedVector269 + crRightHandSideBoundedVector211*crRightHandSideBoundedVector750 - crRightHandSideBoundedVector27*crRightHandSideBoundedVector738 + crRightHandSideBoundedVector27*crRightHandSideBoundedVector741 - crRightHandSideBoundedVector303*crRightHandSideBoundedVector739 - crRightHandSideBoundedVector570*crRightHandSideBoundedVector576 + crRightHandSideBoundedVector737 - crRightHandSideBoundedVector740 + crRightHandSideBoundedVector747) - crRightHandSideBoundedVector580*(-DN_DX_2_0*crRightHandSideBoundedVector581 - DN_DX_2_1*crRightHandSideBoundedVector278*crRightHandSideBoundedVector303*crRightHandSideBoundedVector78 + crRightHandSideBoundedVector585) - crRightHandSideBoundedVector586*(-DN_DX_2_1*crRightHandSideBoundedVector587 - crRightHandSideBoundedVector303*crRightHandSideBoundedVector742 + crRightHandSideBoundedVector591) - crRightHandSideBoundedVector592*(-DN_DX_2_0*crRightHandSideBoundedVector594 - DN_DX_2_1*crRightHandSideBoundedVector117*crRightHandSideBoundedVector297*crRightHandSideBoundedVector303 + crRightHandSideBoundedVector707) - crRightHandSideBoundedVector597*(-DN_DX_2_1*crRightHandSideBoundedVector599 - crRightHandSideBoundedVector303*crRightHandSideBoundedVector743 + crRightHandSideBoundedVector709) - crRightHandSideBoundedVector603*(crRightHandSideBoundedVector342 + crRightHandSideBoundedVector343 + crRightHandSideBoundedVector344 + crRightHandSideBoundedVector492 + crRightHandSideBoundedVector493 + crRightHandSideBoundedVector494 + crRightHandSideBoundedVector727 + crRightHandSideBoundedVector728 + crRightHandSideBoundedVector729 + crRightHandSideBoundedVector744) - crRightHandSideBoundedVector605*(crRightHandSideBoundedVector604 + crRightHandSideBoundedVector730 + crRightHandSideBoundedVector731) - crRightHandSideBoundedVector606*(crRightHandSideBoundedVector710 + crRightHandSideBoundedVector732 + crRightHandSideBoundedVector733) - crRightHandSideBoundedVector607*crRightHandSideBoundedVector715 - crRightHandSideBoundedVector607*crRightHandSideBoundedVector735 + crRightHandSideBoundedVector700;

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

        const double tau_m_avg = 1.0 / ((4.0 * stab_c1 * mu / 3.0 / rho_avg / std::pow(h, 2)) + (stab_c2 * (v_avg_norm + c_avg)/ h));
        const double tau_et_avg = 1.0 / ((stab_c1 * lambda / std::pow(h, 2)) + (stab_c2 * (v_avg_norm + c_avg)/ h));

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
