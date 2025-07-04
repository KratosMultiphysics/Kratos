// KRATOS ___ ___  _  ___   __   ___ ___ ___ ___
//       / __/ _ \| \| \ \ / /__|   \_ _| __| __|
//      | (_| (_) | .` |\ V /___| |) | || _|| _|
//       \___\___/|_|\_| \_/    |___/___|_| |_|  APPLICATION
//
//  License: BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:  
//

// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "custom_elements/conservative_levelset_element.h"
#include "convection_diffusion_application.h"
#include "includes/convection_diffusion_settings.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"

namespace Kratos
{
    template< unsigned int TDim, unsigned int TNumNodes >
    void ConservativeLevelsetElement<TDim,TNumNodes>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        // Resize of the Left and Right Hand side
        if (rLeftHandSideMatrix.size1() != TNumNodes)
            rLeftHandSideMatrix.resize(TNumNodes, TNumNodes, false); //false says not to preserve existing storage!!

        if (rRightHandSideVector.size() != TNumNodes)
            rRightHandSideVector.resize(TNumNodes, false); //false says not to preserve existing storage!!

        //Element variables
        typename BaseType::ElementVariables Variables;
        this->InitializeEulerianElement(Variables,rCurrentProcessInfo);

        // Compute the geometry
        BoundedMatrix<double,TNumNodes, TDim> DN_DX;
        array_1d<double,TNumNodes > N;
        double Volume;
        this-> CalculateGeometry(DN_DX,Volume);

        // Getting the values of shape functions on Integration Points
        BoundedMatrix<double,TNumNodes, TNumNodes> Ncontainer;
        const GeometryType& Geom = this->GetGeometry();
        Ncontainer = Geom.ShapeFunctionsValues( GeometryData::IntegrationMethod::GI_GAUSS_2 );

        // Getting the values of Current Process Info and computing the value of h
        this-> GetNodalValues(Variables,rCurrentProcessInfo);
        double h = this->ComputeH(DN_DX);

        array_1d<double,TNumNodes> div_v_vector;
        array_1d< array_1d<double,3 >, TNumNodes> normal;

        //Computing the divergence
        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            div_v_vector[i] = 1.0;
            normal[i] = this->GetGeometry()[i].FastGetSolutionStepValue(NORMAL);
            for(unsigned int k=0; k<TDim; k++)
            {
                Variables.v[i][k] = normal[i][k];
                Variables.vold[i][k] = normal[i][k];
                Variables.div_v += DN_DX(i,k)*(Variables.v[i][k]*Variables.theta + Variables.vold[i][k]*(1.0-Variables.theta));
            }
        }
        div_v_vector *= Variables.div_v;
        
        //Some auxilary definitions
        BoundedMatrix<double,TNumNodes, TNumNodes> aux1 = ZeroMatrix(TNumNodes, TNumNodes); //terms multiplying dphi/dt
        BoundedMatrix<double,TNumNodes, TNumNodes> aux2 = ZeroMatrix(TNumNodes, TNumNodes); //terms multiplying phi
        BoundedMatrix<double,TNumNodes, TNumNodes> aux3 = ZeroMatrix(TNumNodes, TNumNodes); //terms multiplying volumetric_source
        bounded_matrix<double,TNumNodes, TDim> tmp;

        // Gauss points and Number of nodes coincides in this case.
        for(unsigned int igauss=0; igauss<TNumNodes; igauss++)
        {
            noalias(N) = row(Ncontainer,igauss);

            //obtain the velocity in the middle of the tiem step
            array_1d<double, TDim > vel_gauss=ZeroVector(TDim);

            //obtain the phi and other phi_related variables in the middle of the tiem step 
            double phi_gauss = 0.0;
            double oneminusphi_gauss = 0.0;
            //double oneminus2phi_gauss = 0.0;

            for (unsigned int i = 0; i < TNumNodes; i++)
            {
                phi_gauss += N[i]*(Variables.phi[i]*Variables.theta + Variables.phi_old[i]*(1.0-Variables.theta));
                 
                for(unsigned int k=0; k<TDim; k++)
                    vel_gauss[k] += N[i]*(Variables.v[i][k]*Variables.theta + Variables.vold[i][k]*(1.0-Variables.theta));
            }
       
            //KRATOS_INFO("Gauss point phi:") << phi_gauss << std::endl;
            oneminusphi_gauss = 1.0 - phi_gauss;
            //KRATOS_INFO("Gauss point 1-phi:") << oneminusphi_gauss << std::endl;
            //oneminus2phi_gauss = 1.0 - 2*phi_gauss;
            //KRATOS_INFO("Gauss point 1-2phi:") << oneminus2phi_gauss << std::endl;

            //const double norm_vel = norm_2(oneminus2phi_gauss*vel_gauss);
            //array_1d<double, TNumNodes > a_dot_grad = prod(DN_DX, oneminus2phi_gauss*vel_gauss);
            const double norm_vel = norm_2(oneminusphi_gauss*vel_gauss);
            array_1d<double, TNumNodes > a_dot_grad = prod(DN_DX, oneminusphi_gauss*vel_gauss);

            //const double tau = this->CalculateTau(Variables,norm_vel,h);
            const double tau = 0.0;

            //terms multiplying dphi/dt (aux1)
            noalias(aux1) += (1.0+tau*Variables.beta*Variables.div_v)*outer_prod(N, N);
            noalias(aux1) +=  tau*outer_prod(a_dot_grad, N);

            //terms multiplying volumetric_source (aux3)
            //noalias(aux3) += phi_gauss*oneminusphi_gauss*(1.0+tau*Variables.beta*Variables.div_v)*outer_prod(N, N);
            //noalias(aux3) +=  phi_gauss*oneminusphi_gauss*tau*outer_prod(a_dot_grad, N);

            //terms which multiply the gradient of phi
            //noalias(aux2) += (1.0+tau*Variables.beta*Variables.div_v)*outer_prod(N, a_dot_grad);
            noalias(aux2) -= (1.0+tau*Variables.beta*Variables.div_v)*outer_prod(a_dot_grad, N);
            noalias(aux2) += tau*outer_prod(a_dot_grad, a_dot_grad);
        }

        //adding the second and third term in the formulation
        noalias(rLeftHandSideMatrix)  = (Variables.dt_inv*Variables.density*Variables.specific_heat + Variables.theta*Variables.beta*Variables.div_v)*aux1;
        noalias(rRightHandSideVector) = (Variables.dt_inv*Variables.density*Variables.specific_heat - (1.0-Variables.theta)*Variables.beta*Variables.div_v)*prod(aux1,Variables.phi_old);

        //adding the diffusion
        noalias(rLeftHandSideMatrix)  += (Variables.conductivity * Variables.theta * prod(DN_DX, trans(DN_DX)))*static_cast<double>(TNumNodes);
        noalias(rRightHandSideVector) -= prod((Variables.conductivity * (1.0-Variables.theta) * prod(DN_DX, trans(DN_DX))),Variables.phi_old)*static_cast<double>(TNumNodes) ;

        //terms in aux2
        noalias(rLeftHandSideMatrix) += Variables.density*Variables.specific_heat*Variables.theta*aux2;
        noalias(rRightHandSideVector) -= Variables.density*Variables.specific_heat*(1.0-Variables.theta)*prod(aux2,Variables.phi_old);

        // volume source terms (affecting the RHS only)
        noalias(rRightHandSideVector) += prod(aux1, Variables.volumetric_source);
        //noalias(rRightHandSideVector) -= prod(aux3, div_v_vector);

        //take out the dirichlet part to finish computing the residual
        noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, Variables.phi);

        rRightHandSideVector *= Volume/static_cast<double>(TNumNodes);
        rLeftHandSideMatrix *= Volume/static_cast<double>(TNumNodes);

        KRATOS_CATCH("Error in Conservative Levelset Element")
    }

template class ConservativeLevelsetElement<2,3>;
template class ConservativeLevelsetElement<2,4>;
template class ConservativeLevelsetElement<3,4>;
template class ConservativeLevelsetElement<3,8>;

}