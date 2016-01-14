/*
==============================================================================
KratosConvectionDiffusionApplication 
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu 
rrossi@cimne.upc.edu
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
 */


// System includes 


// External includes 


// Project includes 
#include "includes/define.h"
#include "custom_elements/eulerian_conv_diff.h"
#include "convection_diffusion_application.h"
#include "includes/convection_diffusion_settings.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h" 

namespace Kratos {
   
    template< unsigned int TDim, unsigned int TNumNodes >
    void EulerianConvectionDiffusionElement< TDim, TNumNodes >::CalculateLocalSystem(Matrix& rLeftHandSideMatrix, Vector& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY
        
        if (rLeftHandSideMatrix.size1() != TNumNodes)
            rLeftHandSideMatrix.resize(TNumNodes, TNumNodes, false); //false says not to preserve existing storage!!

        if (rRightHandSideVector.size() != TNumNodes)
            rRightHandSideVector.resize(TNumNodes, false); //false says not to preserve existing storage!!


//         noalias(rLeftHandSideMatrix) = ZeroMatrix(TNumNodes, TNumNodes);
//         noalias(rRightHandSideVector) = ZeroVector(TNumNodes);

        const double theta = rCurrentProcessInfo[THETA]; //Variable defining the temporal scheme (0: Forward Euler, 1: Backward Euler, 0.5: Crank-Nicolson)
        
        const double dyn_st_beta = rCurrentProcessInfo[DYNAMIC_TAU];
        const double delta_t = rCurrentProcessInfo[DELTA_TIME];
        const double dt_inv = 1.0 / delta_t;
		const double lumping_factor = 1.00 / double(TNumNodes);

        ConvectionDiffusionSettings::Pointer my_settings = rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
        const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();

        //getting data for the given geometry
        boost::numeric::ublas::bounded_matrix<double, TNumNodes, TDim > DN_DX;
        array_1d<double, TNumNodes > N;
        double Volume;
        GeometryUtils::CalculateGeometryData(GetGeometry(), DN_DX, N, Volume);
        double h = ComputeH(DN_DX, Volume);
        

        //here we get all the variables we will need
        array_1d<double,TNumNodes> phi, phi_old;
        array_1d< array_1d<double,3 >, TNumNodes> v, vold;
        //bool active_convection=false;    // to kill some terms in case active_convection=false. For the moment it is inactive.
        //using only one Gauss Point for the material properties and volumetric heat flux
        double conductivity = 0.0;
        double specific_heat = 0.0;
        double density = 0.0;

		//storing locally the flags to avoid repeated check in the nodal loops
        const bool IsDefinedVelocityVariable = my_settings->IsDefinedVelocityVariable();
        const bool IsDefinedMeshVelocityVariable = my_settings->IsDefinedMeshVelocityVariable();
        const bool IsDefinedDensityVariable = my_settings->IsDefinedDensityVariable();
        const bool IsDefinedSpecificHeatVariableVariable = my_settings->IsDefinedSpecificHeatVariable();
        const bool IsDefinedDiffusionVariable = my_settings->IsDefinedDiffusionVariable();


        for (unsigned int i = 0; i < TNumNodes; i++)
        {	
            phi[i] = GetGeometry()[i].FastGetSolutionStepValue(rUnknownVar);
            phi_old[i] = GetGeometry()[i].FastGetSolutionStepValue(rUnknownVar,1);
//             dphi_dt[i] = dt_inv*(phi[i] - phi_old [i];

			v[i]=ZeroVector(3);
			vold[i]=ZeroVector(3);
            if (IsDefinedVelocityVariable)
            {
				  const Variable<array_1d<double, 3 > >& rVelocityVar = my_settings->GetVelocityVariable();
				  v[i] = GetGeometry()[i].FastGetSolutionStepValue(rVelocityVar);
				  vold[i] = GetGeometry()[i].FastGetSolutionStepValue(rVelocityVar,1);
				  //active_convection=true;
			}
			
			if (IsDefinedMeshVelocityVariable)
            {
				  const Variable<array_1d<double, 3 > >& rMeshVelocityVar = my_settings->GetMeshVelocityVariable();
				  v[i] -= GetGeometry()[i].FastGetSolutionStepValue(rMeshVelocityVar);
				  vold[i] -= GetGeometry()[i].FastGetSolutionStepValue(rMeshVelocityVar,1);
				  //active_convection=true;
			}
			
			if (IsDefinedDensityVariable)
			{
				const Variable<double>& rDensityVar = my_settings->GetDensityVariable();
				density += GetGeometry()[i].FastGetSolutionStepValue(rDensityVar);
			}
			else
				density += 1.0;
				
			if (IsDefinedSpecificHeatVariableVariable)
			{
				const Variable<double>& rSpecificHeatVar = my_settings->GetSpecificHeatVariable();
				specific_heat += GetGeometry()[i].FastGetSolutionStepValue(rSpecificHeatVar);
			}
			else
				specific_heat += 1.0;			
				
			if (IsDefinedDiffusionVariable)
			{
				const Variable<double>& rDiffusionVar = my_settings->GetDiffusionVariable();
				conductivity += GetGeometry()[i].FastGetSolutionStepValue(rDiffusionVar);
			}
			//if not, then the conductivity = 0   
        }
        array_1d<double,TDim> grad_phi_halfstep = prod(trans(DN_DX), 0.5*(phi+phi_old));
        //const double norm_grad = norm_2(grad_phi_halfstep);

        conductivity *= lumping_factor;
        density *= lumping_factor;
        specific_heat *= lumping_factor;
        //heat_flux *= lumping_factor;



        //here we use a term beta which takes into account a reaction term of the type "beta*div_v"
        //compute the divergence of v
        double div_v = 0.0;
        for (unsigned int i = 0; i < TNumNodes; i++)
            for(unsigned int k=0; k<TDim; k++)
                div_v += DN_DX(i,k)*(v[i][k]*theta + vold[i][k]*(1.0-theta));
            
        double beta = 0.0; //1.0;
        
//         unsigned int nneg=0;
//         for(unsigned int i=0; i<TNumNodes; i++) if(phi[i] < 0.0) nneg++;
//         if(nneg > 0) beta = 1.0; //beta = 0.1;

        boost::numeric::ublas::bounded_matrix<double,TNumNodes, TNumNodes> aux1 = ZeroMatrix(TNumNodes, TNumNodes); //terms multiplying dphi/dt
        boost::numeric::ublas::bounded_matrix<double,TNumNodes, TNumNodes> aux2 = ZeroMatrix(TNumNodes, TNumNodes); //terms multiplying phi
        bounded_matrix<double,TNumNodes, TDim> tmp;

            
        bounded_matrix<double,TNumNodes, TNumNodes> Ncontainer;
        GetShapeFunctionsOnGauss(Ncontainer);
        for(unsigned int igauss=0; igauss<TDim+1; igauss++)
        {
            noalias(N) = row(Ncontainer,igauss);

            //obtain the velocity in the middle of the tiem step
            array_1d<double, TDim > vel_gauss=ZeroVector(TDim);
            for (unsigned int i = 0; i < TNumNodes; i++)
            {
                 for(unsigned int k=0; k<TDim; k++)
                    vel_gauss[k] += N[i]*(v[i][k]*theta + vold[i][k]*(1.0-theta));
            }
            const double norm_vel = norm_2(vel_gauss);
            array_1d<double, TNumNodes > a_dot_grad = prod(DN_DX, vel_gauss);

            const double tau_denom = std::max(dyn_st_beta *dt_inv + 4.0 * conductivity/(h*h) + 2.0 *density*specific_heat* norm_vel / h + beta*div_v,  1e-2);
            const double tau = 1.0 / (tau_denom);

            //terms multiplying dphi/dt (aux1)
            noalias(aux1) += (1.0+tau*beta*div_v)*outer_prod(N, N);
            noalias(aux1) +=  tau*outer_prod(a_dot_grad, N);
            
            //terms which multiply the gradient of phi
            noalias(aux2) += (1.0+tau*beta*div_v)*outer_prod(N, a_dot_grad);
            noalias(aux2) += tau*outer_prod(a_dot_grad, a_dot_grad);
            
            //cross-wind term  
//             if(norm_grad > 1e-3 && norm_vel > 1e-9)
//             {
//                 const double C = 0.7;
//                 const double time_derivative = dt_inv*(inner_prod(N,phi)-inner_prod(N,phi_old));
//                 const double res = -time_derivative -inner_prod(vel_gauss, grad_phi_halfstep);
//                 
//                 const double disc_capturing_coeff = 0.5*C*h*fabs(res/norm_grad);
//                 bounded_matrix<double,TDim,TDim> D = disc_capturing_coeff*( IdentityMatrix(TDim,TDim));
//                 const double norm_vel_squared = norm_vel*norm_vel;
//                 D += (std::max( disc_capturing_coeff - tau*norm_vel_squared , 0.0) - disc_capturing_coeff)/(norm_vel_squared) * outer_prod(vel_gauss,vel_gauss);
// 
//                 noalias(tmp) = prod(DN_DX,D);
//                 noalias(aux2) += prod(tmp,trans(DN_DX));
//             }
        }
        
        //adding the second and third term in the formulation
        noalias(rLeftHandSideMatrix)  = (dt_inv*density*specific_heat + theta*beta*div_v)*aux1;
        noalias(rRightHandSideVector) = (dt_inv*density*specific_heat - (1.0-theta)*beta*div_v)*prod(aux1,phi_old);
        //adding the diffusion
        noalias(rLeftHandSideMatrix)  += (conductivity * theta * prod(DN_DX, trans(DN_DX)))*static_cast<double>(TNumNodes);
        noalias(rRightHandSideVector) -= prod((conductivity * (1.0-theta) * prod(DN_DX, trans(DN_DX))),phi_old)*static_cast<double>(TNumNodes) ;
        
        //terms in aux2
        noalias(rLeftHandSideMatrix) += theta*aux2;
        noalias(rRightHandSideVector) -= (1.0-theta)*prod(aux2,phi_old);
        
        //take out the dirichlet part to finish computing the residual
        noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, phi);

        rRightHandSideVector *= Volume/static_cast<double>(TNumNodes);
        rLeftHandSideMatrix *= Volume/static_cast<double>(TNumNodes);
        
        KRATOS_CATCH("Error in Eulerian ConvDiff Element")
    }

    template< unsigned int TDim, unsigned int TNumNodes >
    void EulerianConvectionDiffusionElement< TDim, TNumNodes >::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_THROW_ERROR(std::runtime_error, "CalculateRightHandSide not implemented","");
    }

    template< unsigned int TDim, unsigned int TNumNodes >
    void EulerianConvectionDiffusionElement< TDim, TNumNodes >::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY
        
        ConvectionDiffusionSettings::Pointer my_settings = rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
        const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();

        if (rResult.size() != TNumNodes)
            rResult.resize(TNumNodes, false);

        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            rResult[i] = GetGeometry()[i].GetDof(rUnknownVar).EquationId();
        }
        
        KRATOS_CATCH("")
    }
    
    template< unsigned int TDim, unsigned int TNumNodes >
    void EulerianConvectionDiffusionElement< TDim, TNumNodes >::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY
        
        ConvectionDiffusionSettings::Pointer my_settings = rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
        const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();

        if (ElementalDofList.size() != TNumNodes)
            ElementalDofList.resize(TNumNodes);

        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            ElementalDofList[i] = GetGeometry()[i].pGetDof(rUnknownVar);

        }
        
        KRATOS_CATCH("");
    }

    template class EulerianConvectionDiffusionElement<2,3>;
    template class EulerianConvectionDiffusionElement<3,4>;
}
