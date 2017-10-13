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
//


// System includes


// External includes


// Project includes
#include "custom_conditions/meshless_load_condition.h"
#include "utilities/math_utils.h"
#include "includes/define.h"

#include "iga_structural_mechanics_application_variables.h"
#include "iga_structural_mechanics_application.h"


namespace Kratos
{
	//************************************************************************************
	//************************************************************************************
	MeshlessLoadCondition::MeshlessLoadCondition(IndexType NewId, GeometryType::Pointer pGeometry)
		: MeshlessBaseCondition(NewId, pGeometry)
	{
		//DO NOT ADD DOFS HERE!!!
	}


	//************************************************************************************
	//************************************************************************************
	MeshlessLoadCondition::MeshlessLoadCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
		: MeshlessBaseCondition(NewId, pGeometry, pProperties)
	{
	}


	Condition::Pointer MeshlessLoadCondition::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
	{
		return MeshlessBaseCondition::Pointer(new MeshlessLoadCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
	}


	// Destructor
	MeshlessLoadCondition::~MeshlessLoadCondition()
	{
	}

	//************************************************************************************
	//************************************************************************************
	void MeshlessLoadCondition::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		MatrixType temp(0, 0);
		CalculateLocalSystem(temp, rRightHandSideVector, rCurrentProcessInfo);
	}

	//************************************************************************************
	//************************************************************************************
	/**
	* CalculateLocalSystem
	* Returns only rRightHandSideVector as there is no impact on the 
	* siffness due to loads
	*
	*/
	void MeshlessLoadCondition::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY

		const unsigned int number_of_points = GetGeometry().size();

		if (rLeftHandSideMatrix.size1() != number_of_points * 3)
			rLeftHandSideMatrix.resize(number_of_points * 3, number_of_points * 3, false);
		noalias(rLeftHandSideMatrix) = ZeroMatrix(number_of_points * 3, number_of_points * 3); //resetting LHS

		if (rRightHandSideVector.size() != number_of_points * 3)
			rRightHandSideVector.resize(number_of_points * 3, false);
		rRightHandSideVector = ZeroVector(number_of_points * 3); //resetting RHS

		Vector ShapeFunctionsN = this->GetValue(SHAPE_FUNCTION_VALUES);
		Matrix DN_De = this->GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);
		double integration_weight = this->GetValue(INTEGRATION_WEIGHT);

		double load = this->GetValue(DISTRIBUTED_LOAD_FACTOR);

		//TODO: find better way to pass condition type
		int conditiontype = this->GetValue(LOAD_TYPE);
		int SURFACE_PRESSURE = (conditiontype % 1000) / 100;
		int SURFACE_DEAD = (conditiontype % 100) / 10;
		int EDGE_LOAD = (conditiontype % 10) / 1;

		array_1d<double, 3> g1, g2, g3; 
		GetBasisVectors(DN_De, g1, g2, g3); //g3 normalized!

		array_1d<double, 3> direction;

		if (SURFACE_PRESSURE == 1)
		{
			CrossProduct(g3, g1, g2);

			double dArea = norm_2(g3);

			integration_weight *= dArea;

			direction = g3 / dArea;
		}
		else if (SURFACE_DEAD == 1)
		{
			direction = this->GetValue(DIRECTION);

			CrossProduct(g3, g1, g2);

			double dArea = norm_2(g3);

			integration_weight *= dArea;
		}
		else if (EDGE_LOAD == 1)
		{
			/** Edge Load:
			* simulates self weight.
			* needs direction of forces.
			* 
			*/
			direction = this->GetValue(DIRECTION);

			Vector Tangents = this->GetValue(TANGENTS);
			double JGeometricToParameter;
			MappingGeometricToParameter(DN_De, Tangents, JGeometricToParameter);

			integration_weight *= JGeometricToParameter;
		}


		Vector fLoads(number_of_points * 3);

		for (unsigned int i = 0; i < number_of_points; i++)
		{
			const array_1d<double, 3> pload = GetGeometry()[i].FastGetSolutionStepValue(POINT_LOAD);
			
			int index = 3 * i;
			fLoads[index]	  = - load * integration_weight * direction[0] * ShapeFunctionsN[i];
			fLoads[index + 1] = - load * integration_weight * direction[1] * ShapeFunctionsN[i];
			fLoads[index + 2] = - load * integration_weight * direction[2] * ShapeFunctionsN[i];
		}

		noalias(rRightHandSideVector) -= fLoads;
		KRATOS_CATCH("")
	}


	//***********************************************************************************
	//***********************************************************************************

	void MeshlessLoadCondition::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
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
	void MeshlessLoadCondition::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo)
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