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
//                   Michael Breitenberger
//                   Riccardo Rossi
//

// System includes

// External includes
#include "includes/define.h"
#include "custom_conditions/meshless_penalty_coupling_rotation_condition.h"
#include "utilities/math_utils.h"

// Project includes
#include "iga_structural_mechanics_application_variables.h"
#include "iga_structural_mechanics_application.h"

namespace Kratos
{

	//************************************************************************************
	//************************************************************************************
	MeshlessPenaltyCouplingRotationCondition::MeshlessPenaltyCouplingRotationCondition(IndexType NewId, GeometryType::Pointer pGeometry)
		: MeshlessBaseCouplingCondition(NewId, pGeometry)
	{
		//DO NOT ADD DOFS HERE!!!
	}


	//************************************************************************************
	//************************************************************************************
	MeshlessPenaltyCouplingRotationCondition::MeshlessPenaltyCouplingRotationCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
		: MeshlessBaseCouplingCondition(NewId, pGeometry, pProperties)
	{
	}


	Condition::Pointer MeshlessPenaltyCouplingRotationCondition::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
	{
		return MeshlessBaseCouplingCondition::Pointer(new MeshlessPenaltyCouplingRotationCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
	}


	// Destructor
	MeshlessPenaltyCouplingRotationCondition::~MeshlessPenaltyCouplingRotationCondition()
	{
	}
	//}
	//************************************************************************************
	//************************************************************************************
	void MeshlessPenaltyCouplingRotationCondition::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY
		const unsigned int number_of_points = GetGeometry().size();

		//resizing the system in case it does not have the right size
		if (rLeftHandSideMatrix.size1() != number_of_points * 3)
			rLeftHandSideMatrix.resize(number_of_points * 3, number_of_points * 3, false);

		noalias(rLeftHandSideMatrix) = ZeroMatrix(number_of_points * 3, number_of_points * 3); //resetting LHS

		//resizing as needed the RHS
		if (rRightHandSideVector.size() != number_of_points * 3)
			rRightHandSideVector.resize(number_of_points * 3, false);

		rRightHandSideVector = ZeroVector(number_of_points * 3); //resetting RHS

		const Vector& ShapeFunctionsN = this->GetValue(SHAPE_FUNCTION_VALUES);
    //KRATOS_WATCH(ShapeFunctionsN)
    const Vector& NSlave = this->GetValue(SHAPE_FUNCTION_VALUES_SLAVE);
    //KRATOS_WATCH(NSlave)
    Vector ShapeFunctions = ZeroVector(ShapeFunctionsN.size() + NSlave.size());
    for (unsigned int i = 0; i < ShapeFunctionsN.size(); i++)
    {
      ShapeFunctions[i] = ShapeFunctionsN[i];
    }
    for (unsigned int i = 0; i < NSlave.size(); i++)
    {
      ShapeFunctions[i + ShapeFunctionsN.size()] = -NSlave[i];
    }
    //KRATOS_WATCH(ShapeFunctions)
    const double Penalty = 10000000;// this->GetValue(PENALTY_FACTOR);
    //KRATOS_WATCH(Penalty)
		const double Weighting = this->GetValue(INTEGRATION_WEIGHT);
    //KRATOS_WATCH(Weighting)
		const Vector& localTrimTangents = this->GetValue(TANGENTS);
    //KRATOS_WATCH(localTrimTangents)
		const Matrix& ShapeFunctionDerivativesMaster = this->GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);
    //KRATOS_WATCH(ShapeFunctionDerivativesMaster)

		array_1d<double, 2> localTrimTangentsMaster;
		localTrimTangentsMaster[0] = localTrimTangents[0];
		localTrimTangentsMaster[1] = localTrimTangents[1];

		//const int displacement_rotation_fix = this->GetValue(DISPLACEMENT_ROTATION_FIX);
    //KRATOS_WATCH(displacement_rotation_fix)

		// Read out information of which elements are fixed
		// int cheaper to store than 4 bool
		//int rot = displacement_rotation_fix / 1000;
		//int dispX = (displacement_rotation_fix % 1000) / 100;
		//int dispY = (displacement_rotation_fix % 100) / 10;
		//int dispZ = (displacement_rotation_fix % 10) / 1;

		if (false)//(rot == 1)
		{
			Vector Phi_r = ZeroVector(number_of_points * 3);
			Vector Phi_r_Lambda = ZeroVector(number_of_points * 3);
			Matrix Phi_rs = ZeroMatrix(number_of_points * 3, number_of_points * 3);
			array_1d<double, 2> Diff_Phi;
			Diff_Phi.clear();

			CaculateRotationalShapeFunctions(Phi_r, Phi_r_Lambda, Phi_rs, Diff_Phi);
      //KRATOS_WATCH(Phi_r)
      //KRATOS_WATCH(Phi_r_Lambda)
      //KRATOS_WATCH(Phi_rs)
      //KRATOS_WATCH(Diff_Phi)

			for (unsigned int i = 0; i < number_of_points * 3; i++)
			{
				for (unsigned int j = 0; j < number_of_points * 3; j++)
				{
					rLeftHandSideMatrix(i, j) = Phi_r(i)*Phi_r(j) + Diff_Phi(0)*Phi_rs(i, j);
				}
				rRightHandSideVector[i] = Diff_Phi(0)*Phi_r(i);
			}
		}

		

		//FOR DISPLACEMENTS
		Matrix Hcomplete = ZeroMatrix(3, number_of_points * 3);
		for (unsigned int i = 0; i < number_of_points; i++)
		{
			if (true)//(dispX == 1)
				Hcomplete(0, 3 * i)     = ShapeFunctions[i];

			if (true)//(dispY == 1)
				Hcomplete(1, 3 * i + 1) = ShapeFunctions[i];

			if (true)//(dispZ == 1)
				Hcomplete(2, 3 * i + 2) = ShapeFunctions[i];
		}
		//KRATOS_WATCH(ShapeFunctions)
		Vector TDisplacements(number_of_points * 3);
		for (unsigned int i = 0; i < number_of_points; i++)
		{
			const array_1d<double, 3> disp = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);
			int index = 3 * i;
			TDisplacements[index]     = disp[0];
			TDisplacements[index + 1] = disp[1];
			TDisplacements[index + 2] = disp[2];
		}

		double JGeometrictoParameter;
		MappingGeometricToParameterMasterElement(ShapeFunctionDerivativesMaster, localTrimTangentsMaster, JGeometrictoParameter);
		
		noalias(rLeftHandSideMatrix) += prod(trans(Hcomplete), Hcomplete);
		noalias(rRightHandSideVector) += prod(prod(trans(Hcomplete), Hcomplete), TDisplacements);

		//Mapping:
		rLeftHandSideMatrix *= (Weighting * JGeometrictoParameter * Penalty);
		rRightHandSideVector *= (Weighting * JGeometrictoParameter * Penalty);

		KRATOS_CATCH("")
	} // MeshlessPenaltyCouplingRotationCondition::MeshlessPenaltyCouplingRotationCondition

	  //************************************************************************************
	  //************************************************************************************
	void MeshlessPenaltyCouplingRotationCondition::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		MatrixType temp(0, 0);
		CalculateLocalSystem(temp, rRightHandSideVector, rCurrentProcessInfo);
	}


	//************************************************************************************
	//************************************************************************************
	void MeshlessPenaltyCouplingRotationCondition::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
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
	void MeshlessPenaltyCouplingRotationCondition::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo)
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


