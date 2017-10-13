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
#include "custom_conditions/meshless_lagrange_coupling_condition.h"
#include "utilities/math_utils.h"
//
#include "iga_structural_mechanics_application.h"
#include "iga_structural_mechanics_application_variables.h"

namespace Kratos
{

//************************************************************************************
//************************************************************************************
MeshlessLagrangeCouplingCondition::MeshlessLagrangeCouplingCondition(IndexType NewId, GeometryType::Pointer pGeometry)
    : MeshlessBaseCouplingCondition(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!
}


//************************************************************************************
//************************************************************************************
MeshlessLagrangeCouplingCondition::MeshlessLagrangeCouplingCondition(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
    : MeshlessBaseCouplingCondition(NewId, pGeometry, pProperties)
{
}


Condition::Pointer MeshlessLagrangeCouplingCondition::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
{
    return MeshlessBaseCouplingCondition::Pointer(new MeshlessLagrangeCouplingCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
}


// Destructor
MeshlessLagrangeCouplingCondition::~MeshlessLagrangeCouplingCondition()
{
}

//************************************************************************************
//************************************************************************************
void MeshlessLagrangeCouplingCondition::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
	KRATOS_TRY

	const unsigned int number_of_points = GetGeometry().size();
	//system equation contains of [(3 displacments)*number_of_points + (3 Lagrange Multipliers)*number_of_points]

	//LHS
	if (rLeftHandSideMatrix.size1() != number_of_points * 6)
		rLeftHandSideMatrix.resize(number_of_points * 6, number_of_points * 6, false);
	noalias(rLeftHandSideMatrix) = ZeroMatrix(number_of_points * 6, number_of_points * 6); //resetting LHS

	//RHS
	if (rRightHandSideVector.size() != number_of_points * 6)
		rRightHandSideVector.resize(number_of_points * 6, false);
	rRightHandSideVector = ZeroVector(number_of_points * 6); //resetting RHS

	const double& Weighting = this->GetValue(INTEGRATION_WEIGHT);
	const Vector& localTrimTangents = this->GetValue(TANGENTS);

	const Vector& ShapeFunctionsN = this->GetValue(SHAPE_FUNCTION_VALUES);
	const Matrix& DN_DeMaster = this->GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES_MASTER);


	//For ROTATIONAL SUPPORT
	Vector Phi_r = ZeroVector(number_of_points * 3);
	Vector Phi_r_Lambda = ZeroVector(number_of_points * 3);
	Matrix Phi_rs = ZeroMatrix(number_of_points * 3, number_of_points * 3);
	array_1d<double, 2> Diff_Phi;
	Diff_Phi.clear();
	
	CaculateRotationalShapeFunctions(Phi_r, Phi_r_Lambda, Phi_rs, Diff_Phi);


	//SHAPE_FUNCTION_VALUES has first the shape functions of the displacements, 
	// then the shape functions of the Lagrange Multipliers
	for (unsigned int i = number_of_points; i < 2*number_of_points; i++) //loop over Lagrange Multipliers
	{
		for (unsigned int j = 0; j < number_of_points; j++) // lopp over shape functions of displacements
		{
			double NN = ShapeFunctionsN[j] * ShapeFunctionsN[i];
			//lambda in X
			rLeftHandSideMatrix(i * 3, 3 * j)         = NN + Phi_r_Lambda((i - number_of_points) * 3 )  * Phi_r(j * 3);//Phi_r_Lambda((i - number_of_points)*3)*Phi_r(j*3); ShapeFunctionsN[i] * Phi_r(j * 3);// 
			//lambda in Y
			rLeftHandSideMatrix(i * 3 + 1, 3 * j + 1) = NN + Phi_r_Lambda((i - number_of_points) * 3+1) * Phi_r(j * 3+1);
			//lambda in Z;
			rLeftHandSideMatrix(i * 3 + 2, 3 * j + 2) = NN + Phi_r_Lambda((i - number_of_points) * 3+2) * Phi_r(j * 3+2);
			//lambda in X
			rLeftHandSideMatrix(3 * j, i * 3)         = NN + Phi_r_Lambda((i - number_of_points) * 3 )  * Phi_r(j * 3);
			//lambda in Y
			rLeftHandSideMatrix(3 * j + 1, i * 3 + 1) = NN + Phi_r_Lambda((i - number_of_points) * 3+1) * Phi_r(j * 3+1);
			//lambda in Z;
			rLeftHandSideMatrix(3 * j + 2, i * 3 + 2) = NN + Phi_r_Lambda((i - number_of_points) * 3+2) * Phi_r(j * 3+2);
		}
	}

	Vector TDisplacementsLambda = ZeroVector(number_of_points * 6);

	for (unsigned int i = 0; i < number_of_points; i++)
	{
		const array_1d<double, 3> disp = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);
		int index = 3 * i;
		TDisplacementsLambda[index]		= disp[0];
		TDisplacementsLambda[index + 1] = disp[1];
		TDisplacementsLambda[index + 2] = disp[2];
	}
	for (unsigned int i = 0; i < number_of_points; i++)
	{
		const array_1d<double, 3> LagrangeMultiplier = GetGeometry()[i].FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER);
		int index = 3 * i + 3 * number_of_points;
		TDisplacementsLambda[index]		= LagrangeMultiplier[0];
		TDisplacementsLambda[index + 1] = LagrangeMultiplier[1];
		TDisplacementsLambda[index + 2] = LagrangeMultiplier[2];
	}

        array_1d<double,2> aux;
        aux[0] = localTrimTangents[0];
        aux[1] = localTrimTangents[1];

	double JGeometrictoParameter;
	MappingGeometricToParameterMasterElement(DN_DeMaster, aux, JGeometrictoParameter);

	rLeftHandSideMatrix *= (Weighting * JGeometrictoParameter);


	noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, TDisplacementsLambda);



	KRATOS_CATCH("")
} // MeshlessLagrangeCouplingCondition::MeshlessLagrangeCouplingCondition


  //************************************************************************************
  //************************************************************************************
void MeshlessLagrangeCouplingCondition::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
	MatrixType temp(0, 0);
	CalculateLocalSystem(temp, rRightHandSideVector, rCurrentProcessInfo);
}


//************************************************************************************
//************************************************************************************
void MeshlessLagrangeCouplingCondition::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
{
	KRATOS_TRY
	unsigned int number_of_nodes = GetGeometry().size();
	unsigned int number_of_dofs = number_of_nodes * 6;

	if (rResult.size() != number_of_dofs)
		rResult.resize(number_of_dofs);

	for (int i = 0; i < number_of_nodes; i++)
	{
		int index = 3 * i;
		rResult[index]	   = GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
		rResult[index + 1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();
		rResult[index + 2] = GetGeometry()[i].GetDof(DISPLACEMENT_Z).EquationId();
	}
	for (int i = 0; i<number_of_nodes; i++)
	{
		int index = 3 * i + 3 * number_of_nodes;
		rResult[index]     = GetGeometry()[i].GetDof(VECTOR_LAGRANGE_MULTIPLIER_X).EquationId();
		rResult[index + 1] = GetGeometry()[i].GetDof(VECTOR_LAGRANGE_MULTIPLIER_Y).EquationId();
		rResult[index + 2] = GetGeometry()[i].GetDof(VECTOR_LAGRANGE_MULTIPLIER_Z).EquationId();
	}

	KRATOS_CATCH("")
}


//************************************************************************************
//************************************************************************************
void MeshlessLagrangeCouplingCondition::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo)
{
	ElementalDofList.resize(0);

	for (unsigned int i = 0; i < GetGeometry().size(); i++)
	{
		ElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
		ElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
		ElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Z));
	}
	for (unsigned int i = 0; i<GetGeometry().size(); i++)
	{
		ElementalDofList.push_back(GetGeometry()[i].pGetDof(VECTOR_LAGRANGE_MULTIPLIER_X));
		ElementalDofList.push_back(GetGeometry()[i].pGetDof(VECTOR_LAGRANGE_MULTIPLIER_Y));
		ElementalDofList.push_back(GetGeometry()[i].pGetDof(VECTOR_LAGRANGE_MULTIPLIER_Z));
	}
}


} // Namespace Kratos


