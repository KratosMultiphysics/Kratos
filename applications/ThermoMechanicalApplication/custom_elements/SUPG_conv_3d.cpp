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


//#define GRADPN_FORM
//#define STOKES

// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "custom_elements/SUPG_conv_3d.h"
#include "utilities/math_utils.h"
#include "thermo_mechanical_application.h"
#include "includes/convection_diffusion_settings.h"
#include "utilities/geometry_utilities.h"
#include "includes/c2c_variables.h"
#include "includes/cfd_variables.h"
#include "includes/deprecated_variables.h"


namespace Kratos
{

//************************************************************************************
//************************************************************************************

SUPGConv3D::SUPGConv3D(IndexType NewId, GeometryType::Pointer pGeometry)
    : SUPGConvDiff3D(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!
}

//************************************************************************************
//************************************************************************************

SUPGConv3D::SUPGConv3D(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : SUPGConvDiff3D(NewId, pGeometry, pProperties)
{

}

Element::Pointer SUPGConv3D::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{

    KRATOS_TRY
    return Element::Pointer(new SUPGConv3D(NewId, GetGeometry().Create(ThisNodes), pProperties));
    KRATOS_CATCH("");
}

SUPGConv3D::~SUPGConv3D()
{
}

//************************************************************************************
//************************************************************************************

void SUPGConv3D::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const ProcessInfo& rConstProcessInfo = rCurrentProcessInfo;
    unsigned int nodes_number = GetGeometry().size();
    unsigned int dim = 3;
    unsigned int matsize = nodes_number;
    const double lumping_factor = 1.00 / double(nodes_number);

    if (rLeftHandSideMatrix.size1() != matsize)
        rLeftHandSideMatrix.resize(matsize, matsize, false); //false says not to preserve existing storage!!

    if (rRightHandSideVector.size() != matsize)
        rRightHandSideVector.resize(matsize, false); //false says not to preserve existing storage!!


    noalias(rLeftHandSideMatrix) = ZeroMatrix(matsize, matsize);
    noalias(rRightHandSideVector) = ZeroVector(matsize);

    double delta_t = rConstProcessInfo[DELTA_TIME];

    BoundedMatrix<double, 4, 3 > DN_DX;
    array_1d<double, 4 > N;

    //getting data for the given geometry
    double Volume;
    GeometryUtils::CalculateGeometryData(GetGeometry(), DN_DX, N, Volume);
    array_1d<double, 3 > ms_vel_gauss;


    //calculating viscosity
    ConvectionDiffusionSettings::Pointer my_settings = rConstProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);

             //const Variable<double>& rDensityVar = my_settings->GetDensityVariable();
    //        const Variable<double>& rSourceVar = my_settings->GetVolumeSourceVariable();
    const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();
    const Variable<array_1d<double, 3 > >& rMeshVelocityVar = my_settings->GetMeshVelocityVariable();
    const Variable<array_1d<double, 3 > >& rConvVar = my_settings->GetConvectionVariable();

    // KRATOS_WATCH(rConvVar);
    //         double specific_heat = GetGeometry()[0].FastGetSolutionStepValue(SPECIFIC_HEAT);
    //         double density = GetGeometry()[0].FastGetSolutionStepValue(rDensityVar);
    //         double heat_source = GetGeometry()[0].FastGetSolutionStepValue(rSourceVar);
    const array_1d<double, 3 > & v = GetGeometry()[0].FastGetSolutionStepValue(rConvVar); //VELOCITY
    const array_1d<double, 3 > & w = GetGeometry()[0].FastGetSolutionStepValue(rMeshVelocityVar); //

	double density = rConstProcessInfo[DENSITY];
    //double air_density = GetGeometry()[0].FastGetSolutionStepValue(rDensityVar);
	double node_distance = GetGeometry()[0].FastGetSolutionStepValue(rUnknownVar);
	int gravity_switch = 	rConstProcessInfo[IS_GRAVITY_FILLING];
	double vel_fac = 1.0;
	if (node_distance > 0.0 && gravity_switch == 1)
		vel_fac = 1.0/density;

    for (unsigned int j = 0; j < dim; j++)
        ms_vel_gauss[j] = vel_fac * (v[j] - w[j]);

    for (unsigned int i = 1; i < nodes_number; i++)
    {
        //             density += GetGeometry()[i].FastGetSolutionStepValue(rDensityVar);
        //             specific_heat += GetGeometry()[i].FastGetSolutionStepValue(SPECIFIC_HEAT);
        //             heat_source += GetGeometry()[i].FastGetSolutionStepValue(rSourceVar);

		density =  rConstProcessInfo[DENSITY];
		node_distance = GetGeometry()[i].FastGetSolutionStepValue(rUnknownVar);
        //air_density = GetGeometry()[i].FastGetSolutionStepValue(rDensityVar);
		vel_fac = 1.0;
		if (node_distance > 0.0 && gravity_switch == 1)
			vel_fac = 1.0/density;

        const array_1d<double, 3 > & v = GetGeometry()[i].FastGetSolutionStepValue(rConvVar);
        const array_1d<double, 3 > & w = GetGeometry()[i].FastGetSolutionStepValue(rMeshVelocityVar);
        for (unsigned int j = 0; j < dim; j++)
            ms_vel_gauss[j] += vel_fac * (v[j] - w[j]);

    }
    //         density *= lumping_factor;
    //         specific_heat *= lumping_factor;
    //         heat_source *= lumping_factor;
    ms_vel_gauss *= lumping_factor;


    //we divide conductivity by (ro*C) and heat_source by C
    // 	heat_source /= (specific_heat);

    double tau;
    CalculateConvTau(ms_vel_gauss, tau, delta_t, Volume, rConstProcessInfo);

    //Crank-Nicholson factor
    double cr_nk = 0.5;
    double dt_inv = 1.0 / delta_t;

    //INERTIA CONTRIBUTION
    BoundedMatrix<double, 4, 4 > msMassFactors = 0.25* IdentityMatrix(4, 4);
    noalias(rLeftHandSideMatrix) = dt_inv * msMassFactors;


    //Advective term
    array_1d<double, 4 > a_dot_grad;
    noalias(a_dot_grad) = prod(DN_DX, ms_vel_gauss);
    BoundedMatrix<double, 4, 4 > Advective_Matrix = outer_prod(N, a_dot_grad);
    noalias(rLeftHandSideMatrix) += (1.0 - cr_nk) * Advective_Matrix;

    //stabilization terms
    array_1d<double, 4 > a_dot_grad_and_mass;
    a_dot_grad_and_mass = dt_inv * N + (1.0 - cr_nk) * a_dot_grad;
    noalias(rLeftHandSideMatrix) += tau * outer_prod(a_dot_grad, a_dot_grad_and_mass);

    //Add heat_source
    // 	noalias(rRightHandSideVector) = heat_source * N;

    //Take N_value terms
    array_1d<double, 4 > step_unknown;
    for (unsigned int iii = 0; iii < nodes_number; iii++)
        step_unknown[iii] = GetGeometry()[iii].FastGetSolutionStepValue(rUnknownVar, 1);

    //compute shock capturing term
    array_1d<double,4> aux_t;
    aux_t[0] = GetGeometry()[0].FastGetSolutionStepValue(rUnknownVar);
    aux_t[1] = GetGeometry()[1].FastGetSolutionStepValue(rUnknownVar);
    aux_t[2] = GetGeometry()[2].FastGetSolutionStepValue(rUnknownVar);
    aux_t[3] = GetGeometry()[3].FastGetSolutionStepValue(rUnknownVar);
    array_1d<double,3> grad_g = prod(trans(DN_DX),aux_t);

    double res = inner_prod(ms_vel_gauss,grad_g);

    double dphi_dt = 0.0;
    for (unsigned int i = 0; i < 4; i++)
        dphi_dt += N[i]*(aux_t[i] - GetGeometry()[i].FastGetSolutionStepValue(rUnknownVar,1));
    res += dt_inv * dphi_dt;

    double h = pow(12.0 * Volume, 0.333333333333333333333);
    h = 0.666666667 * h * 1.732;

//    double Kiso = 0.5*0.7*h*fabs(res)/(norm_2(grad_g) + 1e-12);
////        noalias(rLeftHandSideMatrix) += Kiso * prod(DN_DX,trans(DN_DX));
//
//     double kaniso = Kiso/(inner_prod(ms_vel_gauss,ms_vel_gauss)+1e-12);
//     BoundedMatrix<double, 3, 3 > aux33 = Kiso*IdentityMatrix(3, 3);
//     noalias(aux33) -= kaniso*outer_prod(ms_vel_gauss,ms_vel_gauss);
//
//     BoundedMatrix<double, 3, 4 > aux34 = prod(aux33,trans(DN_DX));
//     noalias(rLeftHandSideMatrix) += prod(DN_DX,aux34);

    //Add N_mass terms
    // 	noalias(rRightHandSideVector) += dt_inv * prod(msMassFactors, step_unknown);
    //
    // 	//Add N_advective terms
    // 	noalias(rRightHandSideVector) -= cr_nk * prod(Advective_Matrix, step_unknown);
    //
    // 	//Add N_Laplacian terms
    // 	noalias(rRightHandSideVector) -= cr_nk * conductivity * prod(Laplacian_Matrix, step_unknown);

    //Add all n_step terms
    BoundedMatrix<double, 4, 4 > old_step_matrix = dt_inv*msMassFactors;
    old_step_matrix -= (cr_nk * Advective_Matrix);
    noalias(rRightHandSideVector) = prod(old_step_matrix, step_unknown);

    //Add n_Stabilization terms
    a_dot_grad_and_mass = dt_inv * N - cr_nk * a_dot_grad;
    double old_res = inner_prod(a_dot_grad_and_mass, step_unknown);
    /*	old_res += heat_source;*/
    noalias(rRightHandSideVector) += tau * a_dot_grad * old_res;

	//add mass conservation term
	//double div_v=0.0;
	//for(unsigned int i=0; i<nodes_number; i++)
	//{
	//	const array_1d<double,3> vel = this->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
	//	for(unsigned int k=0; k<3; k++)
	//		div_v += DN_DX(i,k)* vel[k];
	//}
	//BoundedMatrix<double,4,4> Mconsistent;
	//for(unsigned int i=0; i<nodes_number; i++)
	//{
	//	Mconsistent(i,i) = 0.1;
	//	for(unsigned int j=i+1; j<nodes_number; j++)
	//	{
	//		Mconsistent(i,j) = 0.05;
	//		Mconsistent(j,i) = 0.05;
	//	}
	//}

	//noalias(rLeftHandSideMatrix) -= div_v*Mconsistent;

    //subtracting the dirichlet term
    // RHS -= LHS*temperatures
    for (unsigned int iii = 0; iii < nodes_number; iii++)
        step_unknown[iii] = GetGeometry()[iii].FastGetSolutionStepValue(rUnknownVar);
    noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, step_unknown);


    rRightHandSideVector *= Volume;
    rLeftHandSideMatrix *= Volume;

//        KRATOS_WATCH(this->Id())
//        KRATOS_WATCH(rLeftHandSideMatrix);
//        KRATOS_WATCH(rRightHandSideVector);

    KRATOS_CATCH("")
}
//***********************************************************************************
//**************************************************************************************

void SUPGConv3D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

//     KRATOS_THROW_ERROR(std::logic_error, "method not implemented", "");

    KRATOS_CATCH("")
}




//*************************************************************************************
//*************************************************************************************

void SUPGConv3D::CalculateConvTau(array_1d<double, 3 > & ms_adv_vel, double& tau, const double time, const double volume, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY


    double advvel_norm = MathUtils<double>::Norm3(ms_adv_vel);

    double ele_length = pow(12.0 * volume, 0.333333333333333333333);
    ele_length = 0.666666667 * ele_length * 1.732;

    const double dyn_st_beta = rCurrentProcessInfo[DYNAMIC_TAU];
    tau = 1.0 / (dyn_st_beta / time + 2.0 * advvel_norm / ele_length);


    KRATOS_CATCH("")


}

//************************************************************************************
//************************************************************************************

void SUPGConv3D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
{
    KRATOS_TRY
    unsigned int number_of_nodes = GetGeometry().PointsNumber();
    ConvectionDiffusionSettings::Pointer my_settings = CurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
    const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();

    if (rResult.size() != number_of_nodes)
        rResult.resize(number_of_nodes, false);

    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
        rResult[i] = GetGeometry()[i].GetDof(rUnknownVar).EquationId();
//            KRATOS_WATCH(GetGeometry()[i].GetDof(rUnknownVar).EquationId());
    }
    KRATOS_CATCH("")

}

//************************************************************************************
//************************************************************************************

void SUPGConv3D::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo)
{
    KRATOS_TRY
    unsigned int number_of_nodes = GetGeometry().PointsNumber();
    ConvectionDiffusionSettings::Pointer my_settings = CurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
    const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();

    if (ElementalDofList.size() != number_of_nodes)
        ElementalDofList.resize(number_of_nodes);

    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
        ElementalDofList[i] = GetGeometry()[i].pGetDof(rUnknownVar);

    }
    KRATOS_CATCH("");
}


} // Namespace Kratos


