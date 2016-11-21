//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/IGAStructuralMechanicsApplication/license.txt
//
//  Main authors:    Tobias Tescheamacher
//                   Riccardo Rossi
//


// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "custom_conditions/meshless_support_rotation_condition.h"
#include "utilities/math_utils.h"

#include "iga_structural_mechanics_application_variables.h"
#include "iga_structural_mechanics_application.h"

namespace Kratos
{

//************************************************************************************
//************************************************************************************
MeshlessSupportRotationCondition::MeshlessSupportRotationCondition(IndexType NewId, GeometryType::Pointer pGeometry)
    : MeshlessBaseCondition(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!
}


//************************************************************************************
//************************************************************************************
MeshlessSupportRotationCondition::MeshlessSupportRotationCondition(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
    : MeshlessBaseCondition(NewId, pGeometry, pProperties)
{
}

Condition::Pointer MeshlessSupportRotationCondition::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
{
    return MeshlessBaseCondition::Pointer(new MeshlessSupportRotationCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

// Destructor
MeshlessSupportRotationCondition::~MeshlessSupportRotationCondition()
{
}
//************************************************************************************
//************************************************************************************
void MeshlessSupportRotationCondition::CalculateRotation(const Matrix& ShapeFunctionDerivatives,
	Vector& Phi_r, Matrix& Phi_rs, array_1d<double, 2>& Phi, const array_1d<double, 2>& Tangents)

{
	KRATOS_TRY
	int number_of_points = ShapeFunctionDerivatives.size1();

	//basis vectors g10, g20 and g30 in initial undeformed configuration
	array_1d<double, 3> g10, g20, g30;
	this->GetInitialBasisVectors(ShapeFunctionDerivatives, g10, g20, g30);

	//basis vectors g1, g2 and g3 in current configuration
	array_1d<double, 3> g1, g2, g3;
	this->GetBasisVectors(ShapeFunctionDerivatives, g1, g2, g3);

	// t1 normal to trim, t2 tangential to trim
	array_1d<double, 3> t2 = Tangents(0)*g10 + Tangents(1)*g20;  
	array_1d<double, 3> t1;
	CrossProduct(t1, t2, g30);
	t2 = t2 / norm_2(t2);
	t1 = t1 / norm_2(t1);

	// computation of the a3 displacement
	array_1d<double, 3> w = g3 - g30;
	array_1d<double, 3> SinusOmegaVector;
	CrossProduct(SinusOmegaVector, g30, w);

	array_1d<double, 2> SinusOmega;
	SinusOmega(0) = inner_prod(SinusOmegaVector, t2);
	SinusOmega(1) = inner_prod(SinusOmegaVector, t1);

	array_1d<double, 3> Omega;
	if (SinusOmega(0)>1.0)
		SinusOmega(0) = 0.999999;
	if (SinusOmega(1)>1.0)
		SinusOmega(1) = 0.999999;
	Omega(0) = asin(SinusOmega(0));
	Omega(1) = asin(SinusOmega(1));

// 	array_1d<double, 2> Phi;
	Phi(0) = Omega(0);
	Phi(1) = Omega(1);

	//variation of the a3
	array_1d<double, 3> t3 = g3;
	array_1d<double, 3> tilde_t3; //g3
	CrossProduct(tilde_t3, g1, g2);
	double Length_t3 = norm_2(tilde_t3);

	for (unsigned int n = 0; n < number_of_points; n++)
	{
		for (unsigned int i = 0; i < 3; i++)
		{
			//variations of the basis vectors
			array_1d<double, 3> a1_r;
			a1_r.clear();
			array_1d<double, 3> a2_r;
			a2_r.clear();

			a1_r(i) = ShapeFunctionDerivatives(n, 0);
			a2_r(i) = ShapeFunctionDerivatives(n, 1);
			//variation of the non normalized local vector
			array_1d<double, 3> tilde_3_r = CrossProduct(a1_r, g2) + CrossProduct(g1, a2_r);
			double line_t3_r = inner_prod(t3, tilde_3_r);
			array_1d<double, 3> t3_r = tilde_3_r / Length_t3 - line_t3_r*t3 / Length_t3;
			array_1d<double, 3> SinusOmega_r;
			CrossProduct(SinusOmega_r, g30, t3_r);
			Phi_r(n * 3 + i) = 1.0 / sqrt(1.0 - pow(SinusOmega(0), 2))*inner_prod(SinusOmega_r, t2);
			// if needed at some point:
			//Phi_r_2(i * 3 + j) = 1.0 / sqrt(1.0 - pow(SinusOmega(1), 2))*inner_prod(SinusOmega_r, t1);
		}
	}
	//KRATOS_WATCH(Phi_r)
	for (unsigned int n = 0; n < number_of_points; n++)
	{
		for (unsigned int i = 0; i < 3; i++)
		{
			//variations of the basis vectors
			array_1d<double, 3> a1_r_n;
			a1_r_n.clear();
			array_1d<double, 3> a2_r_n;
			a2_r_n.clear(); /*[0] = 0.0; a2_r_n[1] = 0.0; a2_r_n[2] = 0.0;*/

			a1_r_n(i) = ShapeFunctionDerivatives(n, 0);
			a2_r_n(i) = ShapeFunctionDerivatives(n, 1);

			//variation of the non normalized local vector
			array_1d<double, 3> tilde_3_r_n = CrossProduct(a1_r_n, g2) + CrossProduct(g1, a2_r_n);
			double line_t3_r_n = inner_prod(t3, tilde_3_r_n);
			array_1d<double, 3> t3_r_n = tilde_3_r_n / Length_t3 - line_t3_r_n*t3 / Length_t3;
			array_1d<double, 3> SinusOmega_r_n;
			CrossProduct(SinusOmega_r_n, g30, t3_r_n);

			for (unsigned int m = 0; m < 3; m++)
			{
				for (unsigned int j = 0; j < 3; j++)
				{
					//variations of the basis vectors
					array_1d<double, 3> a1_r_m;
					a1_r_m.clear();
					array_1d<double, 3> a2_r_m;
					a2_r_n.clear();

					a1_r_m(j) = ShapeFunctionDerivatives(m, 0);
					a2_r_m(j) = ShapeFunctionDerivatives(m, 1);

					array_1d<double, 3> tilde_3_r_m = CrossProduct(a1_r_m, g2) + CrossProduct(g1, a2_r_m);
					double line_t3_r_m = inner_prod(t3, tilde_3_r_m);
					array_1d<double, 3> t3_r_m = tilde_3_r_m / Length_t3 - line_t3_r_m*t3 / Length_t3;
					array_1d<double, 3> SinusOmega_r_m = CrossProduct(g30, t3_r_m);


					array_1d<double, 3> tilde_t3_rs = CrossProduct(a1_r_n, a2_r_m) + CrossProduct(a2_r_n, a1_r_m);
					double line_t3_rs = inner_prod(t3_r_m, tilde_3_r_n) + inner_prod(t3, tilde_t3_rs);
					array_1d<double, 3> t3_rs = (tilde_t3_rs*Length_t3 - line_t3_r_m * tilde_3_r_n) / pow(Length_t3, 2)
						- line_t3_rs*t3 / Length_t3 - line_t3_r_n * (t3_r_m * Length_t3 - line_t3_r_m * t3) / pow(Length_t3, 2);
					array_1d<double, 3> SinusOmega_rs = CrossProduct(g30, t3_rs);

					Phi_rs(n * 3 + i, m * 3 + j) = inner_prod(SinusOmega_rs, t2) / sqrt(1.0 - pow(SinusOmega(0), 2))
						+ inner_prod(SinusOmega_r_m, t2)*inner_prod(SinusOmega_r_n, t2)*SinusOmega(0) / pow(1.0
							- pow(SinusOmega(0), 2), 1.5);
				}
			}
		}
	}
	//KRATOS_WATCH(Phi_rs)
	//KRATOS_WATCH(Phi)
	KRATOS_CATCH("")
}


//************************************************************************************
//************************************************************************************
void MeshlessSupportRotationCondition::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
	KRATOS_TRY

    const unsigned int number_of_points = GetGeometry().size();

    //resizing the system in case it does not have the right size
    if(rLeftHandSideMatrix.size1() != number_of_points*3)
        rLeftHandSideMatrix.resize(number_of_points*3,number_of_points*3,false);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(number_of_points*3,number_of_points*3); //resetting LHS

    //resizing as needed the RHS
    if(rRightHandSideVector.size() != number_of_points*3)
        rRightHandSideVector.resize(number_of_points*3,false);
    rRightHandSideVector = ZeroVector(number_of_points*3); //resetting RHS

    //Read in Data
    const Vector& ShapeFunctionsN = this->GetValue(SHAPE_FUNCTION_VALUES);
	const Matrix& DN_De = this->GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);

	const array_1d<double, 3>& support = this->GetValue(DISPLACEMENT);

	double Penalty = this->GetValue(PENALTY_FACTOR);

	const double integration_weight = this->GetValue(INTEGRATION_WEIGHT);
	Vector Tangents = this->GetValue(TANGENTS);

	const int displacement_rotation_fix = this->GetValue(DISPLACEMENT_ROTATION_FIX);

	//MAPPING Parameter: Geometric to Parameter Space
	double JGeometricToParameter;
	MappingGeometricToParameter(DN_De, Tangents, JGeometricToParameter);

        
	// Read out information of which elements are fixed
	// int cheaper to store than 4 bool
	int rot = displacement_rotation_fix / 1000;
	int dispX = (displacement_rotation_fix % 1000) / 100;
	int dispY = (displacement_rotation_fix % 100) / 10;
	int dispZ = (displacement_rotation_fix % 10) / 1;

        
	//For ROTATIONAL CLAMPING
	if (rot == 1)
	{
		Vector Phi_r = ZeroVector(number_of_points * 3);
		Matrix Phi_rs = ZeroMatrix(number_of_points * 3, number_of_points * 3);
		array_1d<double, 2> Phi;
		array_1d<double, 2> TrimTangents;
		TrimTangents[0] = Tangents[0];
		TrimTangents[1] = Tangents[1];
		CalculateRotation(DN_De, Phi_r, Phi_rs, Phi, TrimTangents);

		for (unsigned int n = 0; n < number_of_points*3; n++)
		{
			for (unsigned int m = 0; m < number_of_points*3; m++)
			{
				rLeftHandSideMatrix(n,m) += (Phi_r(n)*Phi_r(m) + Phi(0)*Phi_rs(n,m));
			}
			rRightHandSideVector(n) -= Phi(0)*Phi_r(n);
		}
	}

	
	// DISPLACEMENTS
	Matrix Stiffness =  ZeroMatrix(3, ShapeFunctionsN.size()*3);
	for (unsigned int i = 0; i < ShapeFunctionsN.size(); i++)
	{
		if (dispX == 1)
			Stiffness(0, 3 * i) = ShapeFunctionsN[i];

		if (dispY == 1)
			Stiffness(1, 3 * i + 1) = ShapeFunctionsN[i];

		if (dispZ == 1)
			Stiffness(2, 3 * i + 2) = ShapeFunctionsN[i];
	}

    Vector TDisplacements(number_of_points*3);
	for (unsigned int i = 0; i < number_of_points; i++)
	{
		const array_1d<double, 3> disp = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);
		int index = 3*i;
		TDisplacements[index]     = (disp[0] - support[0]);
		TDisplacements[index + 1] = (disp[1] - support[1]);
		TDisplacements[index + 2] = (disp[2] - support[2]);
	}

    noalias(rLeftHandSideMatrix) += prod(trans(Stiffness), Stiffness);
	noalias(rRightHandSideVector) -= prod(prod(trans(Stiffness), Stiffness), TDisplacements);

        
	//MAPPING
	rLeftHandSideMatrix *= integration_weight * JGeometricToParameter * Penalty;
	rRightHandSideVector *= integration_weight * JGeometricToParameter * Penalty;

        
    KRATOS_CATCH("")
} // MeshlessSupportRotationCondition::MeshlessSupportRotationCondition


//************************************************************************************
//************************************************************************************
void MeshlessSupportRotationCondition::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    MatrixType temp(0,0);
    CalculateLocalSystem(temp, rRightHandSideVector, rCurrentProcessInfo);
}


//************************************************************************************
//************************************************************************************
void MeshlessSupportRotationCondition::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
{
	KRATOS_TRY
		unsigned int number_of_nodes = GetGeometry().size();
	unsigned int dim = number_of_nodes * 3;

	if (rResult.size() != dim)
		rResult.resize(dim);

	for (unsigned int i = 0; i < number_of_nodes; i++)
	{
		int index = i * 3;
		rResult[index] = GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
		rResult[index + 1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();
		rResult[index + 2] = GetGeometry()[i].GetDof(DISPLACEMENT_Z).EquationId();
	}

	KRATOS_CATCH("")
}


//************************************************************************************
//************************************************************************************
void MeshlessSupportRotationCondition::GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
{

	ElementalDofList.resize(0);

	for (unsigned int i = 0; i < GetGeometry().size(); i++)
	{
		ElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
		ElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
		ElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Z));
	}
}


} // Namespace Kratos


