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
    array_1d<double,TNumNodes> N_g;
    BoundedMatrix<double, TNumNodes, TDim> DN_DX_g;
    const auto integration_points = r_geom.IntegrationPoints(mIntegrationMethod);
    for (IndexType g = 0; g < integration_points.size(); ++g) {
        // Get Gauss point data
        noalias(N_g) = row(N, g);
        noalias(DN_DX_g) = DN_DX[g];
        const double w_g = integration_points[g].Weight() * det_J_vect[g];

        // Assemble Gauss point LHS and RHS contributions
        for (IndexType i = 0; i < TNumNodes; ++i) {
            for (IndexType j = 0; j < TNumNodes; ++j) {
                rRightHandSideVector(i) += w_g * Variables.theta * N_g[i] * N_g[j] * Variables.volumetric_source[j];
                rRightHandSideVector(i) -= w_g * (1.0 - Variables.theta) * N_g[i] * N_g[j] * Variables.volumetric_source[j];
            }
        }
    }


    // // Compute the geometry
    // BoundedMatrix<double,TNumNodes, TDim> DN_DX;
    // array_1d<double,TNumNodes > N;
    // double Volume;
    // this-> CalculateGeometry(DN_DX,Volume);

    // // Getting the values of shape functions on Integration Points
    // BoundedMatrix<double,TNumNodes, TNumNodes> Ncontainer;
    // const GeometryType& Geom = this->GetGeometry();
    // Ncontainer = Geom.ShapeFunctionsValues( GeometryData::IntegrationMethod::GI_GAUSS_2 );

    // // Getting the values of Current Process Info and computing the value of h
    // this-> GetNodalValues(Variables,rCurrentProcessInfo);
    // double h = this->ComputeH(DN_DX);

    // //Computing the divergence
    // for (unsigned int i = 0; i < TNumNodes; i++)
    // {
    //     for(unsigned int k=0; k<TDim; k++)
    //     {
    //         Variables.div_v += DN_DX(i,k)*(Variables.v[i][k]*Variables.theta + Variables.vold[i][k]*(1.0-Variables.theta));
    //     }
    // }

    // //Some auxilary definitions
    // BoundedMatrix<double,TNumNodes, TNumNodes> aux1 = ZeroMatrix(TNumNodes, TNumNodes); //terms multiplying dphi/dt
    // BoundedMatrix<double,TNumNodes, TNumNodes> aux2 = ZeroMatrix(TNumNodes, TNumNodes); //terms multiplying phi
    // bounded_matrix<double,TNumNodes, TDim> tmp;

    // // Gauss points and Number of nodes coincides in this case.
    // for(unsigned int igauss=0; igauss<TNumNodes; igauss++)
    // {
    //     noalias(N) = row(Ncontainer,igauss);

    //     //obtain the velocity in the middle of the tiem step
    //     array_1d<double, TDim > vel_gauss=ZeroVector(TDim);
    //     for (unsigned int i = 0; i < TNumNodes; i++)
    //     {
    //             for(unsigned int k=0; k<TDim; k++)
    //             vel_gauss[k] += N[i]*(Variables.v[i][k]*Variables.theta + Variables.vold[i][k]*(1.0-Variables.theta));
    //     }
    //     const double norm_vel = norm_2(vel_gauss);
    //     array_1d<double, TNumNodes > a_dot_grad = prod(DN_DX, vel_gauss);

    //     const double tau = this->CalculateTau(Variables,norm_vel,h);

    //     //terms multiplying dphi/dt (aux1)
    //     noalias(aux1) += (1.0+tau*Variables.beta*Variables.div_v)*outer_prod(N, N);
    //     noalias(aux1) +=  tau*outer_prod(a_dot_grad, N);

    //     //terms which multiply the gradient of phi
    //     noalias(aux2) += (1.0+tau*Variables.beta*Variables.div_v)*outer_prod(N, a_dot_grad);
    //     noalias(aux2) += tau*outer_prod(a_dot_grad, a_dot_grad);
    // }

    // //adding the second and third term in the formulation
    // noalias(rLeftHandSideMatrix)  = (Variables.dt_inv*Variables.density*Variables.specific_heat + Variables.theta*Variables.beta*Variables.div_v)*aux1;
    // noalias(rRightHandSideVector) = (Variables.dt_inv*Variables.density*Variables.specific_heat - (1.0-Variables.theta)*Variables.beta*Variables.div_v)*prod(aux1,Variables.phi_old);

    // //adding the diffusion
    // noalias(rLeftHandSideMatrix)  += (Variables.conductivity * Variables.theta * prod(DN_DX, trans(DN_DX)))*static_cast<double>(TNumNodes);
    // noalias(rRightHandSideVector) -= prod((Variables.conductivity * (1.0-Variables.theta) * prod(DN_DX, trans(DN_DX))),Variables.phi_old)*static_cast<double>(TNumNodes) ;

    // //terms in aux2
    // noalias(rLeftHandSideMatrix) += Variables.density*Variables.specific_heat*Variables.theta*aux2;
    // noalias(rRightHandSideVector) -= Variables.density*Variables.specific_heat*(1.0-Variables.theta)*prod(aux2,Variables.phi_old);

    // // volume source terms (affecting the RHS only)
    // noalias(rRightHandSideVector) += prod(aux1, Variables.volumetric_source);

    // //take out the dirichlet part to finish computing the residual
    // noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, Variables.phi);

    // rRightHandSideVector *= Volume/static_cast<double>(TNumNodes);
    // rLeftHandSideMatrix *= Volume/static_cast<double>(TNumNodes);

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
