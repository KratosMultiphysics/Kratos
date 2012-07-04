/*
==============================================================================
KratosIncompressibleFluidApplication 
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

//   
//   Project Name:        Kratos       
//   Last modified by:    $Author: kazem $
//   Date:                $Date: 2009-01-21 14:15:02 $
//   Revision:            $Revision: 1.6 $
//
//

//#define GRADPN_FORM
//#define STOKES

// System includes 


// External includes 


// Project includes 
#include "includes/define.h"
#include "custom_elements/SUPG_conv_diff_phase_change_2d.h"
#include "utilities/math_utils.h"
#include "thermo_mechanical_application.h"
#include "includes/convection_diffusion_settings.h"
#include "utilities/geometry_utilities.h" 

namespace Kratos {

//************************************************************************************
//************************************************************************************

SUPGConvDiffPhaseChange2D::SUPGConvDiffPhaseChange2D(IndexType NewId, GeometryType::Pointer pGeometry)
: SUPGConvDiff2D(NewId, pGeometry) {
    //DO NOT ADD DOFS HERE!!!
}

//************************************************************************************
//************************************************************************************

SUPGConvDiffPhaseChange2D::SUPGConvDiffPhaseChange2D(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
: SUPGConvDiff2D(NewId, pGeometry, pProperties) {

}

Element::Pointer SUPGConvDiffPhaseChange2D::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const {

    KRATOS_TRY
    return Element::Pointer(new SUPGConvDiffPhaseChange2D(NewId, GetGeometry().Create(ThisNodes), pProperties));
    KRATOS_CATCH("");
}

SUPGConvDiffPhaseChange2D::~SUPGConvDiffPhaseChange2D() {
}

//************************************************************************************
//************************************************************************************

    void SUPGConvDiffPhaseChange2D::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) {
        KRATOS_TRY


        unsigned int nodes_number = GetGeometry().size();
        unsigned int dim = 2;
        unsigned int matsize = nodes_number;
        const double lumping_factor = 1.00 / double(nodes_number);
	

        if (rLeftHandSideMatrix.size1() != matsize)
            rLeftHandSideMatrix.resize(matsize, matsize,false); //false says not to preserve existing storage!!

        if (rRightHandSideVector.size() != matsize)
            rRightHandSideVector.resize(matsize,false); //false says not to preserve existing storage!!


//         noalias(rLeftHandSideMatrix) = ZeroMatrix(matsize, matsize);
//         noalias(rRightHandSideVector) = ZeroVector(matsize);

        double delta_t = rCurrentProcessInfo[DELTA_TIME];

        boost::numeric::ublas::bounded_matrix<double, 3, 2 > DN_DX;
        array_1d<double, 3 > N; 

        //getting data for the given geometry
        double Area;
        GeometryUtils::CalculateGeometryData(GetGeometry(), DN_DX, N, Area);
        array_1d<double, 2 > ms_vel_gauss;
	
	
        //calculating viscosity
        ConvectionDiffusionSettings::Pointer my_settings = rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
	
        const Variable<double>& rDensityVar = my_settings->GetDensityVariable();
        const Variable<double>& rDiffusionVar = my_settings->GetDiffusionVariable();
        const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();
        const Variable<double>& rSourceVar = my_settings->GetVolumeSourceVariable();
        const Variable<array_1d<double, 3 > >& rMeshVelocityVar = my_settings->GetMeshVelocityVariable();
        const Variable<array_1d<double, 3 > >& rConvVar = my_settings->GetConvectionVariable();


        double conductivity = GetGeometry()[0].FastGetSolutionStepValue(rDiffusionVar);
        double specific_heat = GetGeometry()[0].FastGetSolutionStepValue(SPECIFIC_HEAT);
        double density = GetGeometry()[0].FastGetSolutionStepValue(rDensityVar);
        double heat_source = GetGeometry()[0].FastGetSolutionStepValue(rSourceVar);
        const array_1d<double, 3 > & v = GetGeometry()[0].FastGetSolutionStepValue(rConvVar); //VELOCITY
        const array_1d<double, 3 > & w = GetGeometry()[0].FastGetSolutionStepValue(rMeshVelocityVar); //

        double LL = rCurrentProcessInfo[LATENT_HEAT];
	double FF = GetGeometry()[0].FastGetSolutionStepValue(SOLID_FRACTION);
	double old_FF = GetGeometry()[0].FastGetSolutionStepValue(SOLID_FRACTION,1);	
        double DF_DT = GetGeometry()[0].FastGetSolutionStepValue(SOLID_FRACTION_RATE);
	
	double old_T = GetGeometry()[0].FastGetSolutionStepValue(rUnknownVar, 1);
	double T = GetGeometry()[0].FastGetSolutionStepValue(rUnknownVar);
	
        for (unsigned int j = 0; j < dim; j++)
            ms_vel_gauss[j] = v[j] - w[j];

        for (unsigned int i = 1; i < nodes_number; i++)
        {
            conductivity += GetGeometry()[i].FastGetSolutionStepValue(rDiffusionVar);
            density += GetGeometry()[i].FastGetSolutionStepValue(rDensityVar);
            specific_heat += GetGeometry()[i].FastGetSolutionStepValue(SPECIFIC_HEAT);
            heat_source += GetGeometry()[i].FastGetSolutionStepValue(rSourceVar);
	    FF += GetGeometry()[i].FastGetSolutionStepValue(SOLID_FRACTION);
	    old_FF += GetGeometry()[i].FastGetSolutionStepValue(SOLID_FRACTION,1);
	    DF_DT += GetGeometry()[i].FastGetSolutionStepValue(SOLID_FRACTION_RATE);
	    old_T += GetGeometry()[i].FastGetSolutionStepValue(rUnknownVar, 1);	    
	    T += GetGeometry()[i].FastGetSolutionStepValue(rUnknownVar);	 
	    
            const array_1d<double, 3 > & v = GetGeometry()[i].FastGetSolutionStepValue(rConvVar);
            const array_1d<double, 3 > & w = GetGeometry()[i].FastGetSolutionStepValue(rMeshVelocityVar);
            for (unsigned int j = 0; j < dim; j++)
                ms_vel_gauss[j] += v[j] - w[j];

        }
        conductivity *= lumping_factor;
        density *= lumping_factor;
        specific_heat *= lumping_factor;
        heat_source *= lumping_factor;
        ms_vel_gauss *= lumping_factor;	
	FF *= lumping_factor;
	old_FF *= lumping_factor;	
	DF_DT *= lumping_factor;
	old_T *= lumping_factor;
	T *= lumping_factor;
		
	//we divide conductivity by (ro*C) and heat_source by C
// 	conductivity /= (density*specific_heat);
// 	heat_source /= (specific_heat);	
	
        double tau;
	double conductivity_scaled = conductivity/(density*specific_heat);
        CalculateTau(ms_vel_gauss,tau,conductivity_scaled,delta_t, Area, rCurrentProcessInfo);
//         tau = 0.0;
	//Crank-Nicholson factor
	double cr_nk = 0.5;
	double dt_inv = 1.0/ delta_t;
	
        //INERTIA CONTRIBUTION
        boost::numeric::ublas::bounded_matrix<double, 3, 3 > msMassFactors = 1.0 / 3.0 * IdentityMatrix(3, 3);	
        noalias(rLeftHandSideMatrix) = dt_inv * density * specific_heat  * msMassFactors;
	
        //viscous term
	boost::numeric::ublas::bounded_matrix<double, 3, 3 > Laplacian_Matrix = prod(DN_DX , trans(DN_DX));
        noalias(rLeftHandSideMatrix) += (1.0-cr_nk) * conductivity * Laplacian_Matrix;	
		
        //Advective term
        array_1d<double, 3 > a_dot_grad;	
        noalias(a_dot_grad) = prod(DN_DX, ms_vel_gauss);
	boost::numeric::ublas::bounded_matrix<double, 3, 3 > Advective_Matrix = outer_prod(N, a_dot_grad);	
        noalias(rLeftHandSideMatrix) += (1.0-cr_nk) * density * specific_heat * Advective_Matrix;

        //stabilization terms
        array_1d<double, 3 > a_dot_grad_and_mass;
	a_dot_grad_and_mass = density *(dt_inv * N  +  (1.0-cr_nk) * a_dot_grad);
        noalias(rLeftHandSideMatrix) += tau * specific_heat * outer_prod(a_dot_grad, a_dot_grad_and_mass);	
	
	//Add heat_source
	noalias(rRightHandSideVector) = heat_source* N;
	
	//Take N_value terms
	array_1d<double, 3 > step_unknown;
        for (unsigned int iii = 0; iii < nodes_number; iii++)
            step_unknown[iii] =  GetGeometry()[iii].FastGetSolutionStepValue(rUnknownVar, 1);

        //Add N_mass terms
// 	noalias(rRightHandSideVector) += dt_inv * prod(msMassFactors, step_unknown);
// 	
// 	//Add N_advective terms	
// 	noalias(rRightHandSideVector) -= cr_nk * prod(Advective_Matrix, step_unknown);
// 	
// 	//Add N_Laplacian terms		
// 	noalias(rRightHandSideVector) -= cr_nk * conductivity * prod(Laplacian_Matrix, step_unknown);
	
	//Add all n_step terms 
	boost::numeric::ublas::bounded_matrix<double, 3, 3 > old_step_matrix = dt_inv * density * specific_heat * msMassFactors ;
	old_step_matrix -= ( cr_nk * density * specific_heat  * Advective_Matrix + cr_nk * conductivity * Laplacian_Matrix);
	noalias(rRightHandSideVector) += prod(old_step_matrix, step_unknown);
	
	//Add n_Stabilization terms (note thta stabilization for the phase change is just added for "n+1,i" (explicitly) and it is done below)		
	a_dot_grad_and_mass = density * specific_heat*(dt_inv * N  -  cr_nk * a_dot_grad);
	double old_res = inner_prod(a_dot_grad_and_mass, step_unknown);
	old_res += heat_source ;
	noalias(rRightHandSideVector) += tau * a_dot_grad * old_res;	

	//subtracting the dirichlet term
        // RHS -= LHS*temperatures
        for (unsigned int iii = 0; iii < nodes_number; iii++)
            step_unknown[iii] = GetGeometry()[iii].FastGetSolutionStepValue(rUnknownVar);
        noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, step_unknown);	

	
	//add RHS contribution of the phase change considering 3 Gauss points one at each node
	array_1d<double, 3 > phase_change_vec;
	double const_fac = 0.33333333333333333333333333*dt_inv * LL * density;
        for (unsigned int iii = 0; iii < nodes_number; iii++){
            double nd_FF =  GetGeometry()[iii].FastGetSolutionStepValue(SOLID_FRACTION);	
            double nd_FF_old =  GetGeometry()[iii].FastGetSolutionStepValue(SOLID_FRACTION,1);		    
            phase_change_vec[iii] =  (nd_FF - nd_FF_old) * const_fac;	    
	}

	noalias(rRightHandSideVector) += phase_change_vec;

	//add phase_change tangant matrix
        boost::numeric::ublas::bounded_matrix<double, 3, 3 > tan_phase_change =  IdentityMatrix(3, 3);		
	double solid_T = rCurrentProcessInfo[SOLID_TEMPERATURE];
	double fluid_T = rCurrentProcessInfo[FLUID_TEMPERATURE];	
	double mid_T = 0.5*(solid_T + fluid_T);
	double tangent_DF_DT = DF_DT;

        for (unsigned int iii = 0; iii < nodes_number; iii++){
            double nd_DF_DT =  GetGeometry()[iii].FastGetSolutionStepValue(SOLID_FRACTION_RATE);
	    double tangent_DF_DT = nd_DF_DT;
	    if(tangent_DF_DT == 0.0)
	    {
	      double nd_FF =  GetGeometry()[iii].FastGetSolutionStepValue(SOLID_FRACTION);	
	      double nd_FF_old =  GetGeometry()[iii].FastGetSolutionStepValue(SOLID_FRACTION,1);
	      double TT = GetGeometry()[iii].FastGetSolutionStepValue(rUnknownVar);

	      tangent_DF_DT = (nd_FF - nd_FF_old)/(T - mid_T);
	      
	    }
	    tan_phase_change(iii,iii) = 0.33333333333333333333333333 * tangent_DF_DT;
	}
	
// 	if( DF_DT == 0.0 )
// 	  tangent_DF_DT = (FF - old_FF)/(T - mid_T);
	
//  	tangent_DF_DT  = 0.0;
	noalias(rLeftHandSideMatrix) -= dt_inv * density  * LL  * tan_phase_change;  


	
/*KRATOS_WATCH(tan_phase_change);	
KRATOS_WATCH(phase_change_vec* Area);	
KRATOS_WATCH(FF);	
KRATOS_WATCH(old_FF);*/	
	rRightHandSideVector *= Area;	
        rLeftHandSideMatrix *= Area;
		

        KRATOS_CATCH("")
    }



} // Namespace Kratos


