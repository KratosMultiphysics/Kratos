//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Rubio
//

// System includes


// External includes


// Project includes
#include "custom_elements/SUPG_conv_3d_levelset.h"
#include "utilities/math_utils.h"
#include "thermo_mechanical_application.h"
#include "includes/convection_diffusion_settings.h"
#include "utilities/geometry_utilities.h"
#include "includes/c2c_variables.h"
#include "includes/deprecated_variables.h"


#define mu_parameter 0.0001

namespace Kratos {

//************************************************************************************
//************************************************************************************

SUPGConvLevelSet::SUPGConvLevelSet(IndexType NewId, GeometryType::Pointer pGeometry) : SUPGConv3D(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!
}

//************************************************************************************
//************************************************************************************

SUPGConvLevelSet::SUPGConvLevelSet(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) : SUPGConv3D(NewId, pGeometry, pProperties)
{

}

Element::Pointer SUPGConvLevelSet::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{

    KRATOS_TRY
    return Element::Pointer(new SUPGConvLevelSet(NewId, GetGeometry().Create(ThisNodes), pProperties));
    KRATOS_CATCH("");
}

SUPGConvLevelSet::~SUPGConvLevelSet()
{
}

//************************************************************************************
//************************************************************************************

void SUPGConvLevelSet::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
KRATOS_TRY
	//std::cout << "Inside Levelset Element" << std::endl;
	//SUPGConv3D::CalculateLocalSystem(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);
	//////////////////////////////////////////////////////////////////////////////////////////////////////// COPIED FROM SUPG_conv_3d.cpp
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

    double delta_t = rCurrentProcessInfo[DELTA_TIME];

    BoundedMatrix<double, 4, 3 > DN_DX;
    array_1d<double, 4 > N;

    //getting data for the given geometry
    double Volume;
    GeometryUtils::CalculateGeometryData(GetGeometry(), DN_DX, N, Volume);
    array_1d<double, 3 > ms_vel_gauss=ZeroVector(3);

    //calculating viscosity
    ConvectionDiffusionSettings::Pointer my_settings = rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
    const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();
    const Variable<array_1d<double, 3 > >& rMeshVelocityVar = my_settings->GetMeshVelocityVariable();
    const Variable<array_1d<double, 3 > >& rConvVar = my_settings->GetConvectionVariable();
	// Compute Distance Gradient
	//array_1d<double,3> grad_D;
	//CalculateDistanceGradient(grad_D);

    //const array_1d<double, 3 > & v = GetGeometry()[0].FastGetSolutionStepValue(rConvVar); //VELOCITY
    //const array_1d<double, 3 > & w = GetGeometry()[0].FastGetSolutionStepValue(rMeshVelocityVar); //
	double density = rCurrentProcessInfo[DENSITY];
    //double air_density = GetGeometry()[0].FastGetSolutionStepValue(rDensityVar);
	//double node_distance = GetGeometry()[0].FastGetSolutionStepValue(rUnknownVar);
	int gravity_switch = 	rCurrentProcessInfo[IS_GRAVITY_FILLING];
	double vel_fac = 1.0;
    for (unsigned int i = 0; i < nodes_number; i++)
    {
		double node_distance = GetGeometry()[i].FastGetSolutionStepValue(rUnknownVar);
		vel_fac = 1.0;
		if (node_distance > 0.0 && gravity_switch == 1)
			vel_fac = 1.0/density;

        const array_1d<double, 3 > & v = GetGeometry()[i].FastGetSolutionStepValue(rConvVar);
        const array_1d<double, 3 > & w = GetGeometry()[i].FastGetSolutionStepValue(rMeshVelocityVar);

		for (unsigned int j = 0; j < dim; j++)
		{
			double efec_vj=v[j];//-0.5*mu_parameter*grad_D[j]*0.0;
			ms_vel_gauss[j] += vel_fac * (efec_vj - w[j]);
		}

    }
    //         density *= lumping_factor;
    //         specific_heat *= lumping_factor;
    //         heat_source *= lumping_factor;
    ms_vel_gauss *= lumping_factor;


    //we divide conductivity by (ro*C) and heat_source by C
    // 	heat_source /= (specific_heat);

    double tau;
    CalculateConvTau(ms_vel_gauss, tau, delta_t, Volume, rCurrentProcessInfo);

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

    //Add all n_step terms
    BoundedMatrix<double, 4, 4 > old_step_matrix = dt_inv*msMassFactors;
    old_step_matrix -= (cr_nk * Advective_Matrix);
    noalias(rRightHandSideVector) = prod(old_step_matrix, step_unknown);

    //Add n_Stabilization terms
    a_dot_grad_and_mass = dt_inv * N - cr_nk * a_dot_grad;
    double old_res = inner_prod(a_dot_grad_and_mass, step_unknown);
    /*	old_res += heat_source;*/
    noalias(rRightHandSideVector) += tau * a_dot_grad * old_res;


    //subtracting the dirichlet term
    // RHS -= LHS*temperatures
    for (unsigned int iii = 0; iii < nodes_number; iii++)
        step_unknown[iii] = GetGeometry()[iii].FastGetSolutionStepValue(rUnknownVar);
    noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, step_unknown);


    rRightHandSideVector *= Volume;
    rLeftHandSideMatrix *= Volume;

	////////////////////////////////////////////////////////////////////////////////////////////////////////


	VectorType penalty=rRightHandSideVector;
	CalculatePenalty(penalty);
	for(unsigned int i=0; i<rRightHandSideVector.size(); i++)
    {
		rRightHandSideVector[i] -= penalty[i];
    }

KRATOS_CATCH("Error in Levelset Element")
}



//************************************************************************************
//************************************************************************************

void SUPGConvLevelSet::CalculatePenalty(VectorType& penalty)
{
KRATOS_TRY
	//std::cout << "Inside Calculating Penalty" << std::endl;
    //compute geometrical data of the element
    BoundedMatrix<double, 4, 3 > DN_DX;
    array_1d<double, 4 > Ncenter;
    double Volume;
    GeometryUtils::CalculateGeometryData(GetGeometry(), DN_DX, Ncenter, Volume);

	//gather variables on all nodes
    array_1d<double,4> distances;
    for(unsigned int i=0; i<4; i++)
    {
		distances[i] = GetGeometry()[i].FastGetSolutionStepValue(DISTANCE);
    }
    //compute distance gradient
    array_1d<double,3> gradD = prod(trans(DN_DX),distances);
    const double norm_gradD = norm_2(gradD);
	double tmp=0.50*Volume*(norm_gradD-1.0);
	//double tmp=0.50*Volume;
	for(unsigned int i=0; i<penalty.size(); i++)
    {
		penalty[i] = mu_parameter*tmp;
    }
KRATOS_CATCH("Error in Levelset Element, calculating penalty")
}

//************************************************************************************
//************************************************************************************

void SUPGConvLevelSet::CalculateDistanceGradient(array_1d<double,3>& grad_D)
{
KRATOS_TRY
	//grad_D.resize(3);
	//std::cout << "Inside Calculating Penalty" << std::endl;
    //compute geometrical data of the element
    BoundedMatrix<double, 4, 3 > DN_DX;
    array_1d<double, 4 > Ncenter;
    double Volume;
    GeometryUtils::CalculateGeometryData(GetGeometry(), DN_DX, Ncenter, Volume);

	//gather variables on all nodes
    array_1d<double,4> distances;
    for(unsigned int i=0; i<4; i++)
    {
		distances[i] = GetGeometry()[i].FastGetSolutionStepValue(DISTANCE);
    }
    //compute distance gradient
    grad_D = prod(trans(DN_DX),distances);

KRATOS_CATCH("Error in Distance Gradient Calculation")
}


}

