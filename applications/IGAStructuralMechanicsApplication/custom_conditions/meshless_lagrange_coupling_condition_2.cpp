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
#include "custom_conditions/meshless_lagrange_coupling_condition_2.h"
#include "utilities/math_utils.h"
//
#include "iga_structural_mechanics_application.h"
#include "iga_structural_mechanics_application_variables.h"

namespace Kratos
{

//************************************************************************************
//************************************************************************************
MeshlessLagrangeCouplingCondition2::MeshlessLagrangeCouplingCondition2(IndexType NewId, GeometryType::Pointer pGeometry)
    : MeshlessBaseCouplingCondition(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!
}


//************************************************************************************
//************************************************************************************
MeshlessLagrangeCouplingCondition2::MeshlessLagrangeCouplingCondition2(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
    : MeshlessBaseCouplingCondition(NewId, pGeometry, pProperties)
{
}


Condition::Pointer MeshlessLagrangeCouplingCondition2::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
{
    return MeshlessBaseCouplingCondition::Pointer(new MeshlessLagrangeCouplingCondition2(NewId, GetGeometry().Create(ThisNodes), pProperties));
}


// Destructor
MeshlessLagrangeCouplingCondition2::~MeshlessLagrangeCouplingCondition2()
{
}

//************************************************************************************
//************************************************************************************
void MeshlessLagrangeCouplingCondition2::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
	KRATOS_TRY

	const unsigned int number_of_points = GetGeometry().size();
	//system equation contains of [(3 displacments)*number_of_points + (3 Lagrange Multipliers)*number_of_points]

	//LHS
	if (rLeftHandSideMatrix.size1() != number_of_points * 9)
		rLeftHandSideMatrix.resize(number_of_points * 9, number_of_points * 9, false);
	noalias(rLeftHandSideMatrix) = ZeroMatrix(number_of_points * 9, number_of_points * 9); //resetting LHS

	//RHS
	if (rRightHandSideVector.size() != number_of_points * 9)
		rRightHandSideVector.resize(number_of_points * 9, false);
	rRightHandSideVector = ZeroVector(number_of_points * 9); //resetting RHS

	const double& Weighting = this->GetValue(INTEGRATION_WEIGHT);
	const Vector& localTrimTangents = this->GetValue(TANGENTS);

	//Shape functions: first master, second slave, third lambda of master, fourth 0
	const Vector& NMaster = this->GetValue(SHAPE_FUNCTION_VALUES);
  const Vector& NSlave = this->GetValue(SHAPE_FUNCTION_VALUES_SLAVE);
  Vector ShapeFunctions = ZeroVector(NMaster.size() + NSlave.size());
  Vector ShapeFunctionsLambda = ZeroVector(NMaster.size() + NSlave.size());
  for (unsigned int i = 0; i < NMaster.size(); i++)
  {
    ShapeFunctions[i] = NMaster[i];
    ShapeFunctionsLambda[i] = NMaster[i];
  }
  for (unsigned int i = NMaster.size(); i < NMaster.size()+ NSlave.size(); i++)
  {
    ShapeFunctions[i] = -NSlave[i - NMaster.size()];
  }
  std::cout << "here 1" << std::endl;
	//For ROTATIONAL SUPPORT
	Vector Phi_r = ZeroVector(number_of_points * 3);
	Vector Phi_r_Lambda = ZeroVector(number_of_points * 3);
	Matrix Phi_rs = ZeroMatrix(number_of_points * 3, number_of_points * 3);
	array_1d<double, 2> Diff_Phi;
	Diff_Phi.clear();

  std::cout << "here 2" << std::endl;
	CaculateRotationalShapeFunctions(Phi_r, Phi_r_Lambda, Phi_rs, Diff_Phi);

  std::cout << "here 3" << std::endl;
	//SHAPE_FUNCTION_VALUES has first the shape functions of the displacements, 
	// then the shape functions of the Lagrange Multipliers

  std::cout << "here 4" << std::endl;
  for (unsigned int i = 0; i < number_of_points; i++) //loop over Lagrange Multipliers
	{
		for (unsigned int j = 0; j < number_of_points; j++) // lopp over shape functions of displacements
		{
			double NN = ShapeFunctions[j] * ShapeFunctionsLambda[i];
			
      const unsigned int ibase = i * 3 + 3*number_of_points;
      const unsigned int jbase = j*3;
      const unsigned int ilambda = (i) * 3;

			// Matrix in following shape:
			// |0 H^T|
			// |H 0  |
			//lambda in X
			rLeftHandSideMatrix(ibase,     jbase)     = NN;//Phi_r_Lambda((i - number_of_points)*3)*Phi_r(j*3); ShapeFunctionsN[i] * Phi_r(j * 3);// 
			//lambda in Y
			rLeftHandSideMatrix(ibase + 1, jbase + 1) = NN;
			//lambda in Z;
			rLeftHandSideMatrix(ibase + 2, jbase + 2) = NN;
			//lambda in X
			rLeftHandSideMatrix(jbase,     ibase)     = NN;
			//lambda in Y
			rLeftHandSideMatrix(jbase + 1, ibase + 1) = NN;
			//lambda in Z;
			rLeftHandSideMatrix(jbase + 2, ibase + 2) = NN;

			//Depending on the solver if needed. 
			//if (rLeftHandSideMatrix(3 * j, i * 3) <= 0.000001)
			//	rLeftHandSideMatrix(i * 3, i * 3) = 1e-14;
			//if (rLeftHandSideMatrix(3 * j + 1, i * 3 + 1) <= 0.000001)
			//	rLeftHandSideMatrix(i * 3 + 1, i * 3 + 1) = 1e-14;
			//if (rLeftHandSideMatrix(3 * j + 2, i * 3 + 2) <= 0.000001)
			//	rLeftHandSideMatrix(i * 3 +2, i * 3 +2) = 1e-14;
		}
	}

  std::cout << "here 5" << std::endl;
  for (unsigned int i = 0; i < number_of_points; i++) //loop over Lagrange Multipliers
  {
    for (unsigned int j = 0; j < number_of_points; j++) // lopp over shape functions of displacements
    {
      //double NN = ShapeFunctions[j] * ShapeFunctions[i];

      const unsigned int ibase = i * 3;
      const unsigned int jbase = j * 3;
      const unsigned int ilambda = (i ) * 3;

      const unsigned int ibase2 = i * 3 + 6* number_of_points;
      const unsigned int jbase2 = j * 3 + 6* number_of_points;
      // Matrix in following shape:
      // |0 H^T|
      // |H 0  |
      //lambda in X
      rLeftHandSideMatrix(ibase2, jbase)         = Phi_r_Lambda(ilambda)     * Phi_r(jbase);//Phi_r_Lambda((i - number_of_points)*3)*Phi_r(j*3); ShapeFunctionsN[i] * Phi_r(j * 3);// 
                                                                                        //lambda in Y
      rLeftHandSideMatrix(ibase2 + 1, jbase + 1) = Phi_r_Lambda(ilambda + 1) * Phi_r(jbase + 1);
      //lambda in Z;
      rLeftHandSideMatrix(ibase2 + 2, jbase + 2) = Phi_r_Lambda(ilambda + 2) * Phi_r(jbase + 2);
      //lambda in X
      rLeftHandSideMatrix(jbase, ibase2)         = Phi_r_Lambda(ilambda)     * Phi_r(jbase);
      //lambda in Y
      rLeftHandSideMatrix(jbase + 1, ibase2 + 1) = Phi_r_Lambda(ilambda + 1) * Phi_r(jbase + 1);
      //lambda in Z;
      rLeftHandSideMatrix(jbase + 2, ibase2 + 2) = Phi_r_Lambda(ilambda + 2) * Phi_r(jbase + 2);

      //Depending on the solver if needed. 
      //if (rLeftHandSideMatrix(3 * j, i * 3) <= 0.000001)
      //	rLeftHandSideMatrix(i * 3, i * 3) = 1e-14;
      //if (rLeftHandSideMatrix(3 * j + 1, i * 3 + 1) <= 0.000001)
      //	rLeftHandSideMatrix(i * 3 + 1, i * 3 + 1) = 1e-14;
      //if (rLeftHandSideMatrix(3 * j + 2, i * 3 + 2) <= 0.000001)
      //	rLeftHandSideMatrix(i * 3 +2, i * 3 +2) = 1e-14;
    }
  }

  std::cout << "here 6" << std::endl;

	Vector TDisplacementsLambda = ZeroVector(number_of_points * 9);

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

  std::cout << "here 7" << std::endl;
  for (unsigned int i = 0; i < number_of_points; i++)
  {
    const array_1d<double, 3> LagrangeMultiplier = GetGeometry()[i].FastGetSolutionStepValue(ROTATION);
    int index = 3 * i + 6 * number_of_points;
    TDisplacementsLambda[index] = LagrangeMultiplier[0];
    TDisplacementsLambda[index + 1] = LagrangeMultiplier[1];
    TDisplacementsLambda[index + 2] = LagrangeMultiplier[2];
  }
  std::cout << "here 8" << std::endl;
	//MAPPING Geometry Space to Parameter Space
	const Matrix& DN_DeMaster = this->GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);

	array_1d<double, 2> localTrimTangentsMaster;
	localTrimTangentsMaster[0] = localTrimTangents[0];
	localTrimTangentsMaster[1] = localTrimTangents[1];

	double JGeometrictoParameter = 1.0;
	MappingGeometricToParameterMasterElement(DN_DeMaster, localTrimTangentsMaster, JGeometrictoParameter);
	rLeftHandSideMatrix *= (Weighting * JGeometrictoParameter);
	noalias(rRightHandSideVector) += prod(rLeftHandSideMatrix, TDisplacementsLambda);
	//rLeftHandSideMatrix *= (Weighting * JGeometrictoParameter);
	KRATOS_CATCH("")
} // MeshlessLagrangeCouplingCondition2::MeshlessLagrangeCouplingCondition2


  //************************************************************************************
  //************************************************************************************
void MeshlessLagrangeCouplingCondition2::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
	MatrixType temp(0, 0);
	CalculateLocalSystem(temp, rRightHandSideVector, rCurrentProcessInfo);
}


//************************************************************************************
//************************************************************************************
void MeshlessLagrangeCouplingCondition2::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
{
	KRATOS_TRY
	unsigned int number_of_nodes = GetGeometry().size();
	unsigned int number_of_dofs = number_of_nodes * 9;

	if (rResult.size() != number_of_dofs)
		rResult.resize(number_of_dofs);

  std::cout << "here01" << std::endl;
	for (int i = 0; i < number_of_nodes; i++)
	{
		int index = 3 * i;
		rResult[index]	   = GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
		rResult[index + 1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();
		rResult[index + 2] = GetGeometry()[i].GetDof(DISPLACEMENT_Z).EquationId();
	}
  std::cout << "here02" << std::endl;
	for (int i = 0; i<number_of_nodes; i++)
	{
		int index = 3 * i + 3 * number_of_nodes;
		rResult[index]     = GetGeometry()[i].GetDof(VECTOR_LAGRANGE_MULTIPLIER_X).EquationId();
		rResult[index + 1] = GetGeometry()[i].GetDof(VECTOR_LAGRANGE_MULTIPLIER_Y).EquationId();
		rResult[index + 2] = GetGeometry()[i].GetDof(VECTOR_LAGRANGE_MULTIPLIER_Z).EquationId();
	}
  std::cout << "here03" << std::endl;
  for (int i = 0; i<number_of_nodes; i++)
  {
    int index = 3 * i + 6 * number_of_nodes;
    rResult[index] = GetGeometry()[i].GetDof(ROTATION_X).EquationId();
    rResult[index + 1] = GetGeometry()[i].GetDof(ROTATION_Y).EquationId();
    rResult[index + 2] = GetGeometry()[i].GetDof(ROTATION_Z).EquationId();
  }
	KRATOS_CATCH("")
}


//************************************************************************************
//************************************************************************************
void MeshlessLagrangeCouplingCondition2::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo)
{
	ElementalDofList.resize(0);

  std::cout << "here001" << std::endl;
	for (unsigned int i = 0; i < GetGeometry().size(); i++)
	{
		ElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
		ElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
		ElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Z));
	}
  std::cout << "here002" << std::endl;
	for (unsigned int i = 0; i<GetGeometry().size(); i++)
	{
		ElementalDofList.push_back(GetGeometry()[i].pGetDof(VECTOR_LAGRANGE_MULTIPLIER_X));
		ElementalDofList.push_back(GetGeometry()[i].pGetDof(VECTOR_LAGRANGE_MULTIPLIER_Y));
		ElementalDofList.push_back(GetGeometry()[i].pGetDof(VECTOR_LAGRANGE_MULTIPLIER_Z));
	}
  std::cout << "here003" << std::endl;
  for (unsigned int i = 0; i<GetGeometry().size(); i++)
  {
    ElementalDofList.push_back(GetGeometry()[i].pGetDof(ROTATION_X));
    ElementalDofList.push_back(GetGeometry()[i].pGetDof(ROTATION_Y));
    ElementalDofList.push_back(GetGeometry()[i].pGetDof(ROTATION_Z));
  }
}


} // Namespace Kratos


