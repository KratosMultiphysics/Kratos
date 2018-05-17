//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Kazem Kamran
//                   Riccardo Rossi
//

// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "custom_elements/SUPG_conv_diff_3d.h"
#include "utilities/math_utils.h"
#include "thermo_mechanical_application.h"
#include "includes/convection_diffusion_settings.h"
#include "utilities/geometry_utilities.h"
#include "includes/c2c_variables.h"
#include "includes/cfd_variables.h"

namespace Kratos
{

//************************************************************************************
//************************************************************************************

SUPGConvDiff3D::SUPGConvDiff3D(IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!
}

//************************************************************************************
//************************************************************************************

SUPGConvDiff3D::SUPGConvDiff3D(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
{

}

Element::Pointer SUPGConvDiff3D::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{

    KRATOS_TRY
    return Element::Pointer(new SUPGConvDiff3D(NewId, GetGeometry().Create(ThisNodes), pProperties));
    KRATOS_CATCH("");
}

SUPGConvDiff3D::~SUPGConvDiff3D()
{
}

//************************************************************************************
//************************************************************************************

void SUPGConvDiff3D::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    unsigned int nodes_number = GetGeometry().size();
    unsigned int dim = 3;
    unsigned int matsize = nodes_number;
    const double lumping_factor = 1.00 / double(nodes_number);

    if (rLeftHandSideMatrix.size1() != matsize)
        rLeftHandSideMatrix.resize(matsize, matsize,false); //false says not to preserve existing storage!!

    if (rRightHandSideVector.size() != matsize)
        rRightHandSideVector.resize(matsize,false); //false says not to preserve existing storage!!


    noalias(rLeftHandSideMatrix) = ZeroMatrix(matsize, matsize);
    noalias(rRightHandSideVector) = ZeroVector(matsize);

    double delta_t = rCurrentProcessInfo[DELTA_TIME];

    BoundedMatrix<double, 4, 3 > DN_DX;
    array_1d<double, 4 > N;

    //getting data for the given geometry
    double Volume;
    GeometryUtils::CalculateGeometryData(GetGeometry(), DN_DX, N, Volume);
    array_1d<double, 3 > ms_vel_gauss;

    //4 Gauss points coordinates
    BoundedMatrix<double, 4, 4 > GaussCrd = ZeroMatrix(4, 4);
    GaussCrd(0,0) = 0.58541020; GaussCrd(0,1) = 0.13819660; GaussCrd(0,2) = 0.13819660; GaussCrd(0,3) = 0.13819660;
    GaussCrd(1,0) = 0.13819660; GaussCrd(1,1) = 0.58541020; GaussCrd(1,2) = 0.13819660; GaussCrd(1,3) = 0.13819660;
    GaussCrd(2,0) = 0.13819660; GaussCrd(2,1) = 0.13819660; GaussCrd(2,2) = 0.58541020; GaussCrd(2,3) = 0.13819660;
    GaussCrd(3,0) = 0.13819660; GaussCrd(3,1) = 0.13819660; GaussCrd(3,2) = 0.13819660; GaussCrd(3,3) = 0.58541020;

    double wgauss = lumping_factor;


    //calculating viscosity
    ConvectionDiffusionSettings::Pointer my_settings = rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);

    const Variable<double>& rDensityVar = my_settings->GetDensityVariable();
    const Variable<double>& rDiffusionVar = my_settings->GetDiffusionVariable();
    const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();
    const Variable<double>& rSourceVar = my_settings->GetVolumeSourceVariable();
    const Variable<array_1d<double, 3 > >& rMeshVelocityVar = my_settings->GetMeshVelocityVariable();
    const Variable<array_1d<double, 3 > >& rConvVar = my_settings->GetConvectionVariable();

    //compute common matrix terms ( Laplacian, mass)
    BoundedMatrix<double, 4, 4 > Laplacian_Matrix = prod(DN_DX , trans(DN_DX));
    BoundedMatrix<double, 4, 4 > msMassFactors = 1.0 / 4.0 * IdentityMatrix(4, 4);
    //Take N_value terms
    array_1d<double, 4 > step_unknown;
    for (unsigned int iii = 0; iii < nodes_number; iii++)
	step_unknown[iii] =  GetGeometry()[iii].FastGetSolutionStepValue(rUnknownVar, 1);

	//Phase change parameters
    double LL = rCurrentProcessInfo[LATENT_HEAT];
	double solid_T = rCurrentProcessInfo[SOLID_TEMPERATURE];
    double fluid_T = rCurrentProcessInfo[FLUID_TEMPERATURE];
	array_1d<double, 4 > phase_change_vec = ZeroVector(4);
	BoundedMatrix<double, 4, 4 > tan_phase_change =ZeroMatrix(matsize, matsize);
	double mid_T = 0.5*(solid_T + fluid_T);

    //4GP integration rule
    for(unsigned int gp = 0; gp<4; ++gp)
    {
        N[0] = GaussCrd(gp,0); N[1] = GaussCrd(gp,1); N[2] = GaussCrd(gp,2); N[3] = GaussCrd(gp,3);

	double conductivity = N[0]*GetGeometry()[0].FastGetSolutionStepValue(rDiffusionVar);
	double specific_heat = N[0]*GetGeometry()[0].FastGetSolutionStepValue(SPECIFIC_HEAT);
	double density = N[0]*GetGeometry()[0].FastGetSolutionStepValue(rDensityVar);
	double heat_source = N[0]*GetGeometry()[0].FastGetSolutionStepValue(rSourceVar);
	const array_1d<double, 3 > & v = N[0]*GetGeometry()[0].FastGetSolutionStepValue(rConvVar); //VELOCITY
	const array_1d<double, 3 > & w = N[0]*GetGeometry()[0].FastGetSolutionStepValue(rMeshVelocityVar); //
	double gp_dist = N[0]*GetGeometry()[0].FastGetSolutionStepValue(DISTANCE);

	double FF = N[0]*GetGeometry()[0].FastGetSolutionStepValue(SOLIDFRACTION);
	double old_FF = N[0]*GetGeometry()[0].FastGetSolutionStepValue(SOLIDFRACTION,1);
    double DF_DT = N[0]*GetGeometry()[0].FastGetSolutionStepValue(SOLIDFRACTION_RATE);
	double old_T = N[0]*GetGeometry()[0].FastGetSolutionStepValue(rUnknownVar, 1);
	double T =  N[0]*GetGeometry()[0].FastGetSolutionStepValue(rUnknownVar);

	for (unsigned int j = 0; j < dim; j++)
	    ms_vel_gauss[j] = v[j] - w[j];

	for (unsigned int i = 1; i < nodes_number; i++)
	{
	    conductivity += N[i]*GetGeometry()[i].FastGetSolutionStepValue(rDiffusionVar);
	    density += N[i]*GetGeometry()[i].FastGetSolutionStepValue(rDensityVar);
	    specific_heat += N[i]*GetGeometry()[i].FastGetSolutionStepValue(SPECIFIC_HEAT);
	    heat_source += N[i]*GetGeometry()[i].FastGetSolutionStepValue(rSourceVar);

	    FF += N[i]*GetGeometry()[i].FastGetSolutionStepValue(SOLIDFRACTION);
	    old_FF += N[i]*GetGeometry()[i].FastGetSolutionStepValue(SOLIDFRACTION,1);
	    DF_DT += N[i]*GetGeometry()[i].FastGetSolutionStepValue(SOLIDFRACTION_RATE);
	    old_T += N[i]*GetGeometry()[i].FastGetSolutionStepValue(rUnknownVar, 1);
	    T += N[i]*GetGeometry()[i].FastGetSolutionStepValue(rUnknownVar);
		gp_dist += N[i]*GetGeometry()[i].FastGetSolutionStepValue(rUnknownVar);

	    const array_1d<double, 3 > & v = N[i]*GetGeometry()[i].FastGetSolutionStepValue(rConvVar);
	    const array_1d<double, 3 > & w = N[i]*GetGeometry()[i].FastGetSolutionStepValue(rMeshVelocityVar);
	    for (unsigned int j = 0; j < dim; j++)
		ms_vel_gauss[j] += v[j] - w[j];

	}

	//modify C for phase change
	double c_star = 0.0;
	if( gp_dist <= 0.0)
	 {

		//solidification terms LHS
		double tangent_DF_DT = 0.0;
		if(T <= fluid_T)
		{
			tangent_DF_DT = DF_DT;
			if( !(solid_T <= T && T <= fluid_T) ) //if not in the range of phase change
			{
				double aux_denom = T - mid_T;
				if(fabs(aux_denom) < 1e-6)
				 {
				   if(aux_denom >= 0) aux_denom = 1e-6;
				   else aux_denom = -1e-6;
				 }
			   tangent_DF_DT = (FF - old_FF)/aux_denom;
			}
		}
		c_star += tangent_DF_DT * LL;
	  }
	specific_heat += c_star;

	double tau;
	double conductivity_scaled = conductivity/(density*specific_heat);
	CalculateTau(ms_vel_gauss,tau,conductivity_scaled,delta_t, Volume, rCurrentProcessInfo);
	//tau = 0.0;
    //        tau *= density * specific_heat;

	//Crank-Nicholson factor
	double cr_nk = 0.0;
	double dt_inv = 1.0/ delta_t;

	//INERTIA CONTRIBUTION
	noalias(rLeftHandSideMatrix) += dt_inv * density * specific_heat * msMassFactors;

	//viscous term
	noalias(rLeftHandSideMatrix) += (1.0-cr_nk) * conductivity * Laplacian_Matrix;

	//Advective term
	array_1d<double, 4 > a_dot_grad;
	noalias(a_dot_grad) = prod(DN_DX, ms_vel_gauss);
	BoundedMatrix<double, 4, 4 > Advective_Matrix = outer_prod(N, a_dot_grad);
	noalias(rLeftHandSideMatrix) += (1.0-cr_nk) * density * specific_heat * Advective_Matrix;

	//stabilization terms
	array_1d<double, 4 > a_dot_grad_and_mass;
	a_dot_grad_and_mass = density * specific_heat *(dt_inv * N  +  (1.0-cr_nk) * a_dot_grad);
	noalias(rLeftHandSideMatrix) += tau * outer_prod(a_dot_grad, a_dot_grad_and_mass);

	//Add heat_source
	noalias(rRightHandSideVector) += heat_source * N;


	//Add N_mass terms
    // 	noalias(rRightHandSideVector) += dt_inv * prod(msMassFactors, step_unknown);
    //
    // 	//Add N_advective terms
    // 	noalias(rRightHandSideVector) -= cr_nk * prod(Advective_Matrix, step_unknown);
    //
    // 	//Add N_Laplacian terms
    // 	noalias(rRightHandSideVector) -= cr_nk * conductivity * prod(Laplacian_Matrix, step_unknown);

	//Add all n_step terms
	BoundedMatrix<double, 4, 4 > old_step_matrix = dt_inv * density * specific_heat * msMassFactors ;
	old_step_matrix -= ( cr_nk * density * specific_heat * Advective_Matrix + cr_nk * conductivity * Laplacian_Matrix);
	noalias(rRightHandSideVector) += prod(old_step_matrix, step_unknown);

	//Add n_Stabilization terms
	a_dot_grad_and_mass = density * specific_heat * (dt_inv * N  -  cr_nk * a_dot_grad);
	double old_res = inner_prod(a_dot_grad_and_mass, step_unknown);
	old_res += heat_source;
	noalias(rRightHandSideVector) +=  tau * a_dot_grad * old_res;

	//add shock capturing
	double art_visc = 0.0;
	CalculateArtifitialViscosity(art_visc, DN_DX, ms_vel_gauss,rUnknownVar,Volume,conductivity_scaled);
	noalias(rLeftHandSideMatrix) += art_visc * density * specific_heat * Laplacian_Matrix;

	//solidification terms RHS
	/*if( gp_dist <= 0.0)
	 {
		double const_fac = dt_inv * LL ;
		if(T <= fluid_T)
		 {
	       for (unsigned int bb = 0; bb < nodes_number; bb++)
			 phase_change_vec[bb] +=  N[bb]*(FF - old_FF) * const_fac * density;
		 }
		else
		 {
	       for (unsigned int bb = 0; bb < nodes_number; bb++)
			 phase_change_vec[bb] += 0.0;
		 }

		//solidification terms LHS
		double tangent_DF_DT = 0.0;
		if(T <= fluid_T)
		{
			tangent_DF_DT = DF_DT;
			if( !(solid_T <= T && T <= fluid_T) ) //if not in the range of phase change
			{
				double aux_denom = T - mid_T;
				if(fabs(aux_denom) < 1e-6)
				 {
				   if(aux_denom >= 0) aux_denom = 1e-6;
				   else aux_denom = -1e-6;
				 }
			   tangent_DF_DT = (FF - old_FF)/aux_denom;
			}
		}
	 	tan_phase_change += (const_fac * density * tangent_DF_DT * msMassFactors );// + density*LL*tangent_DF_DT*Advective_Matrix);

	  }*/
    }
    //subtracting the dirichlet term
    // RHS -= LHS*temperatures
    for (unsigned int iii = 0; iii < nodes_number; iii++)
	step_unknown[iii] = GetGeometry()[iii].FastGetSolutionStepValue(rUnknownVar);
    noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, step_unknown);

	//phase chaneg term
	//noalias(rRightHandSideVector) += phase_change_vec;
    //noalias(rLeftHandSideMatrix) -= tan_phase_change;

    rRightHandSideVector *= (wgauss * Volume);
    rLeftHandSideMatrix *= (wgauss * Volume);

    KRATOS_CATCH("")
}
//***********************************************************************************
//**************************************************************************************

void SUPGConvDiff3D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

//     KRATOS_THROW_ERROR(std::logic_error, "method not implemented", "");

    KRATOS_CATCH("")
}



//************************************************************************************
//************************************************************************************

void SUPGConvDiff3D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
{
    KRATOS_TRY

    unsigned int number_of_nodes = GetGeometry().PointsNumber();
    ConvectionDiffusionSettings::Pointer my_settings = CurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
    const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();


    if (rResult.size() != number_of_nodes )
        rResult.resize(number_of_nodes , false);

    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
        rResult[i] = GetGeometry()[i].GetDof(rUnknownVar).EquationId();

    }
    KRATOS_CATCH("")

}

//************************************************************************************
//************************************************************************************

void SUPGConvDiff3D::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo)
{
    KRATOS_TRY

    unsigned int number_of_nodes = GetGeometry().PointsNumber();
    ConvectionDiffusionSettings::Pointer my_settings = CurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
    const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();



    if (ElementalDofList.size() != number_of_nodes )
        ElementalDofList.resize(number_of_nodes );

    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
        ElementalDofList[i] = GetGeometry()[i].pGetDof(rUnknownVar);
    }
    KRATOS_CATCH("");
}


//*************************************************************************************
//*************************************************************************************

void SUPGConvDiff3D::CalculateTau(array_1d<double, 3 >& ms_adv_vel, double& tau,const double& K, const double time, const double volume, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY


    double advvel_norm = MathUtils<double>::Norm3(ms_adv_vel);

    double ele_length = pow(12.0*volume,0.333333333333333333333);
    ele_length = 0.666666667 * ele_length * 1.732;

    const double dyn_st_beta = rCurrentProcessInfo[DYNAMIC_TAU];
    tau = 1.0 / (dyn_st_beta / time + 4.0 * K / (ele_length * ele_length) + 2.0 * advvel_norm  / ele_length);


    KRATOS_CATCH("")


}
//*************************************************************************************
//*************************************************************************************
void SUPGConvDiff3D::CalculateArtifitialViscosity(double& art_visc,
						      BoundedMatrix<double, 4, 3 > DN_DX,
						      array_1d<double, 3 > ms_vel_gauss,
						      const Variable<double>& temperature,
						      const double volume,
						      const double scaled_K)
{
	KRATOS_TRY
    art_visc = 0.0;
	int not_cutted = 1;
	int negative = 0;
	double pos_dist_max = -1000.0;
	double neg_dist_max = +1000.0;

	double ele_length = pow(12.0*volume,0.333333333333333333333);
	ele_length = 0.666666667 * ele_length * 1.732;

	for( int ii = 0; ii<4; ++ii)
	{
	 double nd_dist = this->GetGeometry()[ii].FastGetSolutionStepValue(DISTANCE);
     if( nd_dist < 0.0)
 	   negative++;

     if( nd_dist > 0.0 && nd_dist > pos_dist_max)
		 pos_dist_max = nd_dist;

	 if( nd_dist < 0.0 && nd_dist < neg_dist_max)
		 neg_dist_max = nd_dist;
    }

 if( negative != 4 && negative != 0)
 	  not_cutted = 0;



	 array_1d<double, 3 > grad_t;
         unsigned int number_of_nodes = GetGeometry().PointsNumber();
	 double temp_cup = GetGeometry()[0].FastGetSolutionStepValue(temperature);
	 grad_t[0] = DN_DX(0,0) * temp_cup;
	 grad_t[1] = DN_DX(0,1) * temp_cup;
	 grad_t[2] = DN_DX(0,2) * temp_cup;

	for(unsigned int ii=1; ii < number_of_nodes; ++ii)
	{
	   temp_cup = GetGeometry()[ii].FastGetSolutionStepValue(temperature);
	   grad_t[0] += DN_DX(ii,0) * temp_cup;
	   grad_t[1] += DN_DX(ii,1) * temp_cup;
	   grad_t[2] += DN_DX(ii,2) * temp_cup;
	}

	double norm_grad_t = MathUtils<double>::Norm3(grad_t);
	if(norm_grad_t < 0.000001){
	  art_visc = 0.0;
	}
	else
	{
	  double a_parallel = (ms_vel_gauss[0]*grad_t[0] + ms_vel_gauss[1]*grad_t[1] + ms_vel_gauss[2]*grad_t[2]) / norm_grad_t;
	  a_parallel = fabs(a_parallel);


	  double Effective_K = scaled_K;
	  if(scaled_K < 0.00000000001)
	    Effective_K = 0.00000000001;

	  double Peclet_parallel = a_parallel*ele_length / (2.0*Effective_K);
	  double alpha;
	  if(Peclet_parallel == 0.0)
		 alpha = 0.0;
	  else
		 alpha =1.0 - 1.0/Peclet_parallel;

	  if (alpha < 0.0)
	       alpha = 0.0;

	  //art_visc = 100.0*0.5 * alpha * ele_length * a_parallel;
	  art_visc =  1.0 *  alpha * ele_length * a_parallel;

	}

   // for cut elements apply higher art_visc
	//if(not_cutted == 0)
      //  art_visc*=1.0;//25.0

	if( (pos_dist_max <= 1.0*ele_length && pos_dist_max > 0.0) || not_cutted == 0 )//(neg_dist_max < 0.0 && fabs(neg_dist_max) < 1.0*ele_length)
	{
		art_visc*=25.0;
	}

	this->GetValue(PR_ART_VISC) = art_visc;

	KRATOS_CATCH("")

}
//************************************************************************************
//************************************************************************************
//void SUPGConvDiff3D::GetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo)
//{
//
//    /*        double delta_t = rCurrentProcessInfo[DELTA_TIME];*/
//    BoundedMatrix<double, 4, 3 > DN_DX;
//    array_1d<double, 4 > N;
//    //getting data for the given geometry
//    double volume;
//    GeometryUtils::CalculateGeometryData(GetGeometry(), DN_DX, N, volume);
//
//    if (rVariable == PR_ART_VISC)
//    {
//	for (unsigned int PointNumber = 0;
//		PointNumber < 1; PointNumber++)
//	{
//	    //	KRATOS_WATCH(this->GetValue(IS_WATER));
//	    //	KRATOS_WATCH(this->Info());
//	    rValues[PointNumber] = this->GetValue(PR_ART_VISC);
//	}
//    }
//
//}

} // Namespace Kratos


