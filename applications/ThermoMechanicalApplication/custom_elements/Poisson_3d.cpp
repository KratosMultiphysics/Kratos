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
//


//#define GRADPN_FORM
//#define STOKES

// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "custom_elements/Poisson_3d.h"
#include "utilities/math_utils.h"
#include "thermo_mechanical_application.h"
#include "includes/convection_diffusion_settings.h"
#include "utilities/geometry_utilities.h"

#include <cmath>

namespace Kratos
{

//************************************************************************************
//************************************************************************************

Poisson3D::Poisson3D(IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!
}

//************************************************************************************
//************************************************************************************

Poisson3D::Poisson3D(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
{

}

Element::Pointer Poisson3D::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{

    KRATOS_TRY
    return Element::Pointer(new Poisson3D(NewId, GetGeometry().Create(ThisNodes), pProperties));
    KRATOS_CATCH("");
}

Poisson3D::~Poisson3D()
{
}

//************************************************************************************
//************************************************************************************

void Poisson3D::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
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


//         noalias(rLeftHandSideMatrix) = ZeroMatrix(matsize, matsize);
//         noalias(rRightHandSideVector) = ZeroVector(matsize);


    BoundedMatrix<double, 4, 3 > DN_DX;
    array_1d<double, 4 > N;

    //getting data for the given geometry
    double Volume;
    GeometryUtils::CalculateGeometryData(GetGeometry(), DN_DX, N, Volume);
    array_1d<double, 3 > ms_vel_gauss;


    //calculating viscosity
    ConvectionDiffusionSettings::Pointer my_settings = rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);

     //const Variable<double>& rDiffusionVar = my_settings->GetDiffusionVariable();
    const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();
//     const Variable<double>& rSourceVar = my_settings->GetVolumeSourceVariable();

    double mass_source = 0.0;

    double ele_length = pow(12.0*Volume,0.333333333333333333333);
    double conductivity = rCurrentProcessInfo[KAPPA];
    double epsilon = rCurrentProcessInfo[EPSILON];
//     conductivity = 10.0 * ele_length;
//     epsilon = 3.0 * ele_length;

    for (unsigned int i = 0; i < nodes_number; i++){
        mass_source += GetGeometry()[i].FastGetSolutionStepValue(VOLUME_FRACTION);

	double raw_dist =  GetGeometry()[i].FastGetSolutionStepValue(DISTANCE);
	raw_dist += GetGeometry()[i].FastGetSolutionStepValue(DISTANCE_CORRECTION);

	mass_source -= SmoothedHeavyside(epsilon,raw_dist);
	GetGeometry()[i].FastGetSolutionStepValue(AUX_INDEX) = GetGeometry()[i].FastGetSolutionStepValue(VOLUME_FRACTION) - SmoothedHeavyside(epsilon,raw_dist);
    }
    mass_source *= lumping_factor;


    //viscous term
    BoundedMatrix<double, 4, 4 > Laplacian_Matrix = prod(DN_DX , trans(DN_DX));
    noalias(rLeftHandSideMatrix) =  conductivity * Laplacian_Matrix;


    //Add heat_source
    noalias(rRightHandSideVector) = mass_source * N;




    //add shock capturing
/*    double art_visc = 0.0;
    CalculateArtifitialViscosity(art_visc, DN_DX, ms_vel_gauss,rUnknownVar,Volume,conductivity_scaled);

    noalias(rLeftHandSideMatrix) += art_visc * density * specific_heat  * Laplacian_Matrix;	*/


    //subtracting the dirichlet term
    // RHS -= LHS*temperatures
    array_1d<double, 4 > step_unknown;
    for (unsigned int iii = 0; iii < nodes_number; iii++)
        step_unknown[iii] = GetGeometry()[iii].FastGetSolutionStepValue(DISTANCE_CORRECTION);
    noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, step_unknown);


    rRightHandSideVector *= Volume;
    rLeftHandSideMatrix *= Volume;


    KRATOS_CATCH("")
}
//***********************************************************************************
//**************************************************************************************

void Poisson3D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

//     KRATOS_THROW_ERROR(std::logic_error, "method not implemented", "");

    KRATOS_CATCH("")
}



//************************************************************************************
//************************************************************************************

void Poisson3D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
{
    KRATOS_TRY

    unsigned int number_of_nodes = GetGeometry().PointsNumber();
    ConvectionDiffusionSettings::Pointer my_settings = CurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
    const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();


    if (rResult.size() != number_of_nodes )
        rResult.resize(number_of_nodes , false);

    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
        rResult[i] = GetGeometry()[i].GetDof(DISTANCE_CORRECTION).EquationId();

    }
    KRATOS_CATCH("")

}

//************************************************************************************
//************************************************************************************

void Poisson3D::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo)
{
    KRATOS_TRY

    unsigned int number_of_nodes = GetGeometry().PointsNumber();
    ConvectionDiffusionSettings::Pointer my_settings = CurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
    const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();



    if (ElementalDofList.size() != number_of_nodes )
        ElementalDofList.resize(number_of_nodes );

    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
        ElementalDofList[i] = GetGeometry()[i].pGetDof(DISTANCE_CORRECTION);
    }
    KRATOS_CATCH("");
}
//*************************************************************************************
//*************************************************************************************

double Poisson3D::SmoothedHeavyside( const double epsilon, const double distance)
{
    KRATOS_TRY



    if(distance <= -epsilon)
      return 1.0;
    else if(distance >= epsilon)
      return 0.0;
    else{
     double PI = 3.14159265;
     double h_eps = 0.5*(1.0 - distance/epsilon - sin(distance * PI/epsilon) / PI);
     return h_eps;
    }

    KRATOS_CATCH("")

}

//*************************************************************************************
//*************************************************************************************

void Poisson3D::CalculateTau(array_1d<double, 3 >& ms_adv_vel, double& tau,const double& K, const double time, const double volume, const ProcessInfo& rCurrentProcessInfo)
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
void Poisson3D::CalculateArtifitialViscosity(double& art_visc,
						      BoundedMatrix<double, 4, 3 > DN_DX,
						      array_1d<double, 3 > ms_vel_gauss,
						      const Variable<double>& temperature,
						      const double volume,
						      const double scaled_K)
{
	KRATOS_TRY

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
	if(norm_grad_t < 0.000001)
	  art_visc = 0.0;
	else
	{
	  double a_parallel = (ms_vel_gauss[0]*grad_t[0] + ms_vel_gauss[1]*grad_t[1] + ms_vel_gauss[2]*grad_t[2]) / norm_grad_t;
	  a_parallel = abs(a_parallel);
	  double ele_length = pow(12.0*volume,0.333333333333333333333);
	  ele_length = 0.666666667 * ele_length * 1.732;

	  double Effective_K = scaled_K;
	  if(scaled_K < 0.00000000001)
	    Effective_K = 0.00000000001;

	  double Peclet_parallel = a_parallel*ele_length / (2.0*Effective_K);
	  double alpha;
	  if(Peclet_parallel == 0.0)
		 alpha = 0.0;
	  else
		 alpha = 3.0 - 1.0/Peclet_parallel;

	  if (alpha < 0.0)
	       alpha = 0.0;

	  art_visc = 100.0*0.5 * alpha * ele_length * a_parallel;

	}

        this->GetValue(PR_ART_VISC) = art_visc;

	KRATOS_CATCH("")

}
//************************************************************************************
//************************************************************************************
//void Poisson3D::GetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo)
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


