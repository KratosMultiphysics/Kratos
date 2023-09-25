// KRATOS ___ ___  _  ___   __   ___ ___ ___ ___
//       / __/ _ \| \| \ \ / /__|   \_ _| __| __|
//      | (_| (_) | .` |\ V /___| |) | || _|| _|
//       \___\___/|_|\_| \_/    |___/___|_| |_|  APPLICATION
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "custom_elements/eulerian_conv_diff.h"
#include "convection_diffusion_application.h"
#include "includes/convection_diffusion_settings.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"

// Application includes
#include "axisymmetric_eulerian_convection_diffusion.h"

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
void AxisymmetricEulerianConvectionDiffusionElement<TDim,TNumNodes>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Resize of LHS and RHS arrays
    if (rLeftHandSideMatrix.size1() != TNumNodes || rLeftHandSideMatrix.size2() != TNumNodes) {
        rLeftHandSideMatrix.resize(TNumNodes, TNumNodes, false);
    }

    if (rRightHandSideVector.size() != TNumNodes) {
        rRightHandSideVector.resize(TNumNodes, false);
    }

    // Initialize LHS and RHS arrays
    rLeftHandSideMatrix.clear();
    rRightHandSideVector.clear();

    // Initialize element data container
    typename BaseType::ElementVariables Variables;
    this->InitializeEulerianElement(Variables,rCurrentProcessInfo);

    // Fill element data container with nodal data
    this->GetNodalValues(Variables, rCurrentProcessInfo);

    // Calculate kinematics
    Vector det_J_vect;
    ShapeFunctionsGradientsType DN_DX;
    const auto& r_geom = this->GetGeometry();
    const auto N = r_geom.ShapeFunctionsValues(mIntegrationMethod);
    r_geom.ShapeFunctionsIntegrationPointsGradients(DN_DX, det_J_vect, mIntegrationMethod);

    // Gauss points loop
    double y_g;
    array_1d<double,TDim> v_g_th;
    array_1d<double,TNumNodes> N_g;
    array_1d<double,TNumNodes> v_g_th_dot_grad;
    BoundedMatrix<double, TDim, TDim> grad_v_g_th;
    BoundedMatrix<double, TNumNodes, TDim> DN_DX_g;
    const auto integration_points = r_geom.IntegrationPoints(mIntegrationMethod);
    const SizeType n_gauss = integration_points.size();
    for (IndexType g = 0; g < n_gauss; ++g) {
        // Get Gauss point data
        noalias(N_g) = row(N, g);
        noalias(DN_DX_g) = DN_DX[g];
        const double w_g = integration_points[g].Weight() * det_J_vect[g];

        // Calculate Gauss point values
        y_g = 0.0;
        noalias(v_g_th) = ZeroVector(TDim);
        noalias(grad_v_g_th) = ZeroMatrix(TDim, TDim);
        for (IndexType i = 0; i < TNumNodes; ++i) {
            // Gauss point radius
            y_g += N_g[i] * r_geom[i].Y();
            // Gauss point velocity and velocity gradient
            for (IndexType d1 = 0; d1 < TDim; ++d1) {
                // Gauss point velocity
                v_g_th[d1] += N_g[i] * (Variables.theta * Variables.v[i][d1] + (1.0 - Variables.theta) * Variables.vold[i][d1]);
                // Gauss point velocity gradient (stored as Dv_j/Dx_i)
                for (IndexType d2 = 0; d2 < TDim; ++d2) {
                    grad_v_g_th(d1, d2) += DN_DX_g(i, d1) * (Variables.theta * Variables.v[i][d2] + (1.0 - Variables.theta) * Variables.vold[i][d2]);
                }
            }
        }
        noalias(v_g_th_dot_grad) = prod(DN_DX_g, v_g_th);

        // Calculate stabilization constant
        const double norm_v = norm_2(v_g_th);
        const double h = this->ComputeH(DN_DX_g) / n_gauss;
        const double tau = this->CalculateTau(Variables, norm_v, h);

        // Assemble Gauss point LHS and RHS contributions
        for (IndexType i = 0; i < TNumNodes; ++i) {
            for (IndexType j = 0; j < TNumNodes; ++j) {
                // Source term
                const double aux_source = w_g * N_g[i] * N_g[j];
                rRightHandSideVector(i) += aux_source * Variables.theta * Variables.volumetric_source[j];
                rRightHandSideVector(i) += aux_source * (1.0 - Variables.theta) * Variables.volumetric_source[j]; // In case we make the body force time dependent

                // Dynamic term
                const double aux_dyn = w_g * Variables.density * Variables.specific_heat * Variables.dt_inv * N_g[i] * N_g[j];
                rLeftHandSideMatrix(i, j) += aux_dyn;
                rRightHandSideVector(i) -= aux_dyn * (Variables.phi[j] - Variables.phi_old[j]);

                // Convective terms
                const double aux_conv_1 = w_g * Variables.density * Variables.specific_heat * N_g[i] * v_g_th_dot_grad[j];
                rLeftHandSideMatrix(i,j) += aux_conv_1 * Variables.theta;
                rRightHandSideVector(i) -= aux_conv_1 * Variables.theta * Variables.phi[j];
                rRightHandSideVector(i) -= aux_conv_1 * (1.0 - Variables.theta) * Variables.phi_old[j];

                const double aux_conv_2 = w_g * Variables.density * Variables.specific_heat * Variables.beta * N_g[i] * (v_g_th[1]/y_g + grad_v_g_th(0,0) + grad_v_g_th(1,1));
                rLeftHandSideMatrix(i,j) += aux_conv_2 * Variables.theta;
                rRightHandSideVector(i) -= aux_conv_2 * Variables.theta * Variables.phi[j];
                rRightHandSideVector(i) -= aux_conv_2 * (1.0 - Variables.theta) * Variables.phi_old[j];

                // Diffusive terms
                const double aux_diff_1 = w_g * Variables.conductivity * N_g[i] * DN_DX_g(j,1) / y_g;
                rLeftHandSideMatrix(i,j) -= aux_diff_1 * Variables.theta;
                rRightHandSideVector(i) += aux_diff_1 * Variables.theta * Variables.phi[j];
                rRightHandSideVector(i) += aux_diff_1 * (1.0 - Variables.theta) * Variables.phi_old[j];

                const double aux_diff_2 = w_g * Variables.conductivity * (DN_DX_g(i, 0) * DN_DX_g(j, 0) + DN_DX_g(i, 1) * DN_DX_g(j, 1));
                rLeftHandSideMatrix(i, j) += aux_diff_2 * Variables.theta;
                rRightHandSideVector(i) -= aux_diff_2 * Variables.theta * Variables.phi[j];
                rRightHandSideVector(i) -= aux_diff_2 * (1.0 - Variables.theta) * Variables.phi_old[j];

                // Stabilization terms
                const double aux_div_stab_op = Variables.density * Variables.specific_heat * Variables.beta * N_g[i] * (v_g_th[1] / y_g + grad_v_g_th(0, 0) + grad_v_g_th(1, 1));
                const double aux_conv_stab_op = Variables.density * Variables.specific_heat * (DN_DX_g(i,0)*v_g_th[0] + N_g[i]*grad_v_g_th(0,0) + DN_DX_g(i,1)*v_g_th[1] + N_g[i]*grad_v_g_th(1,1));
                const double aux_diff_stab_op = Variables.conductivity * (DN_DX_g(i, 1) / y_g - N_g[i] / std::pow(y_g, 2));
                const double aux_stab = tau * w_g * (aux_conv_stab_op - aux_div_stab_op - aux_diff_stab_op);

                rRightHandSideVector(i) += aux_stab * N_g[j] * Variables.theta * Variables.volumetric_source[j];
                rRightHandSideVector(i) += aux_stab * N_g[j] * (1.0 - Variables.theta) * Variables.volumetric_source[j]; // In case we make the body force time dependent

                const double aux_dyn_stab = aux_stab * Variables.density * Variables.specific_heat * Variables.dt_inv * N_g[j];
                rLeftHandSideMatrix(i, j) += aux_dyn_stab;
                rRightHandSideVector(i) -= aux_dyn_stab * (Variables.phi[j] - Variables.phi_old[j]);

                const double aux_conv_stab = aux_stab * Variables.density * Variables.specific_heat * v_g_th_dot_grad[j];
                rLeftHandSideMatrix(i, j) += aux_conv_stab * Variables.theta;
                rRightHandSideVector(i) -= aux_conv_stab * Variables.theta * Variables.phi[j];
                rRightHandSideVector(i) -= aux_conv_stab * (1.0 - Variables.theta) * Variables.phi_old[j];

                const double aux_div_stab = aux_stab * Variables.density * Variables.specific_heat * Variables.beta * (v_g_th[1] / y_g + grad_v_g_th(0, 0) + grad_v_g_th(1, 1)) * N_g[j];
                rLeftHandSideMatrix(i, j) += aux_div_stab * Variables.theta;
                rRightHandSideVector(i) -= aux_div_stab * Variables.theta * Variables.phi[j];
                rRightHandSideVector(i) -= aux_div_stab * (1.0 - Variables.theta) * Variables.phi_old[j];

                const double aux_diff_stab = aux_stab * DN_DX_g(j, 1) * Variables.conductivity / y_g;
                rLeftHandSideMatrix(i, j) -= aux_diff_stab * Variables.theta;
                rRightHandSideVector(i) += aux_diff_stab * Variables.theta * Variables.phi[j];
                rRightHandSideVector(i) += aux_diff_stab * (1.0 - Variables.theta) * Variables.phi_old[j];
            }
        }
    }

    KRATOS_CATCH("")
}

template< unsigned int TDim, unsigned int TNumNodes >
void AxisymmetricEulerianConvectionDiffusionElement< TDim, TNumNodes >::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    // Matrix LeftHandSide;
    // this->CalculateLocalSystem(LeftHandSide,rRightHandSideVector,rCurrentProcessInfo);
}

template< unsigned int TDim, unsigned int TNumNodes >
int AxisymmetricEulerianConvectionDiffusionElement< TDim, TNumNodes >::Check(const ProcessInfo &rCurrentProcessInfo) const
{
    // Base element check
    int out = BaseType::Check(rCurrentProcessInfo);
    KRATOS_ERROR_IF_NOT(out == 0)
        << "Error in base class Check for Element " << this->Info() << std::endl
        << "Error code is " << out << std::endl;

    // Check that there are no negative y-coordinates (radius is always positive)
    const auto& r_geom = this->GetGeometry();
    for (const auto& r_node : r_geom) {
        KRATOS_ERROR_IF(r_node.Y() < 0.0) << "Negative y-coordinate found in node " << r_node.Id() << ". Axisymmetric radius must be positive." << std::endl;
    }

    return 0;
}

template class AxisymmetricEulerianConvectionDiffusionElement<2,3>;
template class AxisymmetricEulerianConvectionDiffusionElement<2,4>;

}
