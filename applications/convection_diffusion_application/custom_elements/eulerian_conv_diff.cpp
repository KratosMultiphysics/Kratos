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

namespace Kratos 
{
     
//----------------------------------------------------------------------------------------

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

//---------------------------------------------------------------------------------------- 

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
        
        KRATOS_CATCH("")
    }
    
//----------------------------------------------------------------------------------------

    template< unsigned int TDim, unsigned int TNumNodes >
    void EulerianConvectionDiffusionElement<TDim,TNumNodes>::CalculateLocalSystem(Matrix& rLeftHandSideMatrix, Vector& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY
        
        // Resize of the Left and Right Hand side
        if (rLeftHandSideMatrix.size1() != TNumNodes)
            rLeftHandSideMatrix.resize(TNumNodes, TNumNodes, false); //false says not to preserve existing storage!!

        if (rRightHandSideVector.size() != TNumNodes)
            rRightHandSideVector.resize(TNumNodes, false); //false says not to preserve existing storage!!
                  
        //Element variables
        ElementVariables Variables;
        this->InitializeEulerianElement(Variables,rCurrentProcessInfo);
        
        // Compute the geometry
        boost::numeric::ublas::bounded_matrix<double,TNumNodes, TDim> DN_DX;
        array_1d<double,TNumNodes > N;
        double Volume;
        this-> CalculateGeometry(DN_DX,Volume);

        // Getting the values of shape functions on Integration Points
        boost::numeric::ublas::bounded_matrix<double,TNumNodes, TNumNodes> Ncontainer;
        const GeometryType& Geom = this->GetGeometry();
        Ncontainer = Geom.ShapeFunctionsValues( GeometryData::GI_GAUSS_2 );

        // Getting the values of Current Process Info and computing the value of h
        this-> GetNodalValues(Variables,rCurrentProcessInfo);
        double h = this->ComputeH(DN_DX);
        
        //Computing the divergence
        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            for(unsigned int k=0; k<TDim; k++)
            {
                Variables.div_v += DN_DX(i,k)*(Variables.v[i][k]*Variables.theta + Variables.vold[i][k]*(1.0-Variables.theta));
            }
        }
        
        //Some auxilary definitions
        boost::numeric::ublas::bounded_matrix<double,TNumNodes, TNumNodes> aux1 = ZeroMatrix(TNumNodes, TNumNodes); //terms multiplying dphi/dt
        boost::numeric::ublas::bounded_matrix<double,TNumNodes, TNumNodes> aux2 = ZeroMatrix(TNumNodes, TNumNodes); //terms multiplying phi
        bounded_matrix<double,TNumNodes, TDim> tmp;
        
        // Gauss points and Number of nodes coincides in this case.
        for(unsigned int igauss=0; igauss<TNumNodes; igauss++)
        {
            noalias(N) = row(Ncontainer,igauss);

            //obtain the velocity in the middle of the tiem step
            array_1d<double, TDim > vel_gauss=ZeroVector(TDim);
            for (unsigned int i = 0; i < TNumNodes; i++)
            {
                 for(unsigned int k=0; k<TDim; k++)
                    vel_gauss[k] += N[i]*(Variables.v[i][k]*Variables.theta + Variables.vold[i][k]*(1.0-Variables.theta));
            }
            const double norm_vel = norm_2(vel_gauss);
            array_1d<double, TNumNodes > a_dot_grad = prod(DN_DX, vel_gauss);

            const double tau_denom = std::max(Variables.dyn_st_beta *Variables.dt_inv + 4.0 * Variables.conductivity/(h*h) + 2.0 *Variables.density*Variables.specific_heat* norm_vel / h + Variables.beta*Variables.div_v,  1e-2);
            const double tau = 1.0 / (tau_denom);

            //terms multiplying dphi/dt (aux1)
            noalias(aux1) += (1.0+tau*Variables.beta*Variables.div_v)*outer_prod(N, N);
            noalias(aux1) +=  tau*outer_prod(a_dot_grad, N);
            
            //terms which multiply the gradient of phi
            noalias(aux2) += (1.0+tau*Variables.beta*Variables.div_v)*outer_prod(N, a_dot_grad);
            noalias(aux2) += tau*outer_prod(a_dot_grad, a_dot_grad);
        }
                 
        //adding the second and third term in the formulation
        noalias(rLeftHandSideMatrix)  = (Variables.dt_inv*Variables.density*Variables.specific_heat + Variables.theta*Variables.beta*Variables.div_v)*aux1;
        noalias(rRightHandSideVector) = (Variables.dt_inv*Variables.density*Variables.specific_heat - (1.0-Variables.theta)*Variables.beta*Variables.div_v)*prod(aux1,Variables.phi_old);
                
        //adding the diffusion
        noalias(rLeftHandSideMatrix)  += (Variables.conductivity * Variables.theta * prod(DN_DX, trans(DN_DX)))*static_cast<double>(TNumNodes);
        noalias(rRightHandSideVector) -= prod((Variables.conductivity * (1.0-Variables.theta) * prod(DN_DX, trans(DN_DX))),Variables.phi_old)*static_cast<double>(TNumNodes) ;
        
        //terms in aux2
        noalias(rLeftHandSideMatrix) += Variables.theta*aux2;
        noalias(rRightHandSideVector) -= (1.0-Variables.theta)*prod(aux2,Variables.phi_old);
        
        //take out the dirichlet part to finish computing the residual
        noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, Variables.phi);

        rRightHandSideVector *= Volume/static_cast<double>(TNumNodes);
        rLeftHandSideMatrix *= Volume/static_cast<double>(TNumNodes);
        
        KRATOS_CATCH("Error in Eulerian ConvDiff Element")
    }
 
//----------------------------------------------------------------------------------------

    template< unsigned int TDim, unsigned int TNumNodes >
    void EulerianConvectionDiffusionElement<TDim,TNumNodes>::InitializeEulerianElement(ElementVariables& rVariables, ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        rVariables.theta = rCurrentProcessInfo[THETA]; //Variable defining the temporal scheme (0: Forward Euler, 1: Backward Euler, 0.5: Crank-Nicolson)
        rVariables.dyn_st_beta = rCurrentProcessInfo[DYNAMIC_TAU];
        const double delta_t = rCurrentProcessInfo[DELTA_TIME];
        rVariables.dt_inv = 1.0 / delta_t;
		rVariables.lumping_factor = 1.00 / double(TNumNodes);
        
        rVariables.conductivity = 0.0;
        rVariables.specific_heat = 0.0;
        rVariables.density = 0.0;
        rVariables.beta = 0.0;
        rVariables.div_v = 0.0;

   
        KRATOS_CATCH( "" )
    }

//---------------------------------------------------------------------------------------- 

    template< unsigned int TDim, unsigned int TNumNodes >
    void EulerianConvectionDiffusionElement<TDim,TNumNodes>::CalculateGeometry(boost::numeric::ublas::bounded_matrix<double,TNumNodes,TDim>& rDN_DX, double& rVolume)
    {
                
        const GeometryType& Geom = this->GetGeometry();
        
        // We select GI_GAUSS_1 due to we are computing at the barycenter.
        const GeometryType::IntegrationPointsArrayType& integration_points = Geom.IntegrationPoints( GeometryData::GI_GAUSS_1 );
        const unsigned int NumGPoints = integration_points.size();
        rVolume = Geom.Area();
        GeometryType::ShapeFunctionsGradientsType DN_DXContainer(NumGPoints);
        Geom.ShapeFunctionsIntegrationPointsGradients(DN_DXContainer,GeometryData::GI_GAUSS_1);
        
        noalias( rDN_DX ) = DN_DXContainer[0];
        
    }

//----------------------------------------------------------------------------------------

    template< unsigned int TDim, unsigned int TNumNodes >
    double EulerianConvectionDiffusionElement< TDim, TNumNodes >::ComputeH(boost::numeric::ublas::bounded_matrix<double,TNumNodes,TDim >& DN_DX)
    {
        double h=0.0;
                
        for(unsigned int i=0; i<TNumNodes; i++)
        {
            double h_inv = 0.0;
            for(unsigned int k=0; k<TDim; k++)
            {
                h_inv += DN_DX(i,k)*DN_DX(i,k);
            }
            h += 1.0/h_inv;
        }
        h = sqrt(h)/static_cast<double>(TNumNodes);
        return h;
    }

//---------------------------------------------------------------------------------------- 
    
    template< unsigned int TDim, unsigned int TNumNodes >
    void EulerianConvectionDiffusionElement<TDim,TNumNodes>::GetNodalValues(ElementVariables& rVariables, ProcessInfo& rCurrentProcessInfo)
    {       
        ConvectionDiffusionSettings::Pointer my_settings = rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
        
        //////storing locally the flags to avoid repeated check in the nodal loops
        const bool IsDefinedVelocityVariable = my_settings->IsDefinedVelocityVariable();
        const bool IsDefinedMeshVelocityVariable = my_settings->IsDefinedMeshVelocityVariable();
        const bool IsDefinedDensityVariable = my_settings->IsDefinedDensityVariable();
        const bool IsDefinedSpecificHeatVariableVariable = my_settings->IsDefinedSpecificHeatVariable();
        const bool IsDefinedDiffusionVariable = my_settings->IsDefinedDiffusionVariable();

        const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();
        
        for (unsigned int i = 0; i < TNumNodes; i++)
        {	
            rVariables.phi[i] = GetGeometry()[i].FastGetSolutionStepValue(rUnknownVar);
            rVariables.phi_old[i] = GetGeometry()[i].FastGetSolutionStepValue(rUnknownVar,1);
            //dphi_dt[i] = dt_inv*(phi[i] - phi_old [i];

			rVariables.v[i]=ZeroVector(3);
			rVariables.vold[i]=ZeroVector(3);
            if (IsDefinedVelocityVariable)
            {
				  const Variable<array_1d<double, 3 > >& rVelocityVar = my_settings->GetVelocityVariable();
				  rVariables.v[i] = GetGeometry()[i].FastGetSolutionStepValue(rVelocityVar);
				  rVariables.vold[i] = GetGeometry()[i].FastGetSolutionStepValue(rVelocityVar,1);
				  //active_convection=true;
			}
			
			if (IsDefinedMeshVelocityVariable)
            {
				  const Variable<array_1d<double, 3 > >& rMeshVelocityVar = my_settings->GetMeshVelocityVariable();
				  rVariables.v[i] -= GetGeometry()[i].FastGetSolutionStepValue(rMeshVelocityVar);
				  rVariables.vold[i] -= GetGeometry()[i].FastGetSolutionStepValue(rMeshVelocityVar,1);
				  //active_convection=true;
			}
			
			if (IsDefinedDensityVariable)
			{
				const Variable<double>& rDensityVar = my_settings->GetDensityVariable();
				rVariables.density += GetGeometry()[i].FastGetSolutionStepValue(rDensityVar);
			}
			else
				rVariables.density += 1.0;
				
			if (IsDefinedSpecificHeatVariableVariable)
			{
				const Variable<double>& rSpecificHeatVar = my_settings->GetSpecificHeatVariable();
				rVariables.specific_heat += GetGeometry()[i].FastGetSolutionStepValue(rSpecificHeatVar);
			}
			else
				rVariables.specific_heat += 1.0;			
				
			if (IsDefinedDiffusionVariable)
			{
				const Variable<double>& rDiffusionVar = my_settings->GetDiffusionVariable();
				rVariables.conductivity += GetGeometry()[i].FastGetSolutionStepValue(rDiffusionVar);
			}
			//if not, then the conductivity = 0   
        }
        
        //array_1d<double,TDim> grad_phi_halfstep = prod(trans(DN_DX), 0.5*(phi+phi_old));
        //const double norm_grad = norm_2(grad_phi_halfstep);

        rVariables.conductivity *= rVariables.lumping_factor;
        rVariables.density *= rVariables.lumping_factor;
        rVariables.specific_heat *= rVariables.lumping_factor;

    }

//----------------------------------------------------------------------------------------

    template< unsigned int TDim, unsigned int TNumNodes >
    void EulerianConvectionDiffusionElement< TDim, TNumNodes >::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_THROW_ERROR(std::runtime_error, "CalculateRightHandSide not implemented","");
    }

//----------------------------------------------------------------------------------------


template class EulerianConvectionDiffusionElement<2,3>;
template class EulerianConvectionDiffusionElement<2,4>;
template class EulerianConvectionDiffusionElement<3,4>;
template class EulerianConvectionDiffusionElement<3,8>;

}
