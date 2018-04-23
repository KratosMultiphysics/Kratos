//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/IGAStructuralMechanicsApplication/license.txt
//
//  Main authors:    Tobias Teschemacher
//


// System includes


// External includes


// Project includes
#include "custom_conditions/meshless_point_condition.h"
#include "utilities/math_utils.h"
#include "includes/define.h"

#include "iga_structural_mechanics_application_variables.h"
#include "iga_structural_mechanics_application.h"


namespace Kratos
{
	//***********************************************************************************
	//***********************************************************************************
	void MeshlessPointCondition::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
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
	void MeshlessPointCondition::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo)
	{
		ElementalDofList.resize(0);

		for (unsigned int i = 0; i < GetGeometry().size(); i++)
		{
			ElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
			ElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
			ElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Z));
		}
	}

	//************************************************************************************
	//************************************************************************************
	void MeshlessPointCondition::GetFirstDerivativesVector(Vector& values, int Step)
	{
		if (values.size() != 3)
			values.resize(3, false);

		const int& number_of_nodes = GetGeometry().size();
		Vector N = this->GetValue(SHAPE_FUNCTION_VALUES);

		for (SizeType i = 0; i < number_of_nodes; i++)
		{
			const NodeType & iNode = GetGeometry()[i];
			const array_1d<double, 3>& vel = iNode.FastGetSolutionStepValue(VELOCITY, Step);

			values[0] += N[i] * vel[0];
			values[1] += N[i] * vel[1];
			values[2] += N[i] * vel[2];
		}
	}


	/***********************************************************************************/
	/***********************************************************************************/
	void MeshlessPointCondition::CalculateOnIntegrationPoints(
		const Variable<array_1d<double, 3>>& rVariable,
		std::vector<array_1d<double, 3>>& rOutput,
		const ProcessInfo& rCurrentProcessInfo
	)
	{
		if (rOutput.size() != 1)
			rOutput.resize(1);

		if (rVariable == DISPLACEMENT) {
			const int& number_of_nodes = GetGeometry().size();
			Vector N = this->GetValue(SHAPE_FUNCTION_VALUES);

			for (SizeType i = 0; i < number_of_nodes; i++)
			{
				const NodeType& iNode = GetGeometry()[i];
				const array_1d<double, 3>& disp = iNode.FastGetSolutionStepValue(DISPLACEMENT);

				rOutput[0][0] += N[i] * disp[0];
				rOutput[0][1] += N[i] * disp[1];
				rOutput[0][2] += N[i] * disp[2];
			}
		}
	}

	/***********************************************************************************/
	/***********************************************************************************/
	void MeshlessPointCondition::GetValueOnIntegrationPoints(
		const Variable<array_1d<double, 3>>& rVariable,
		std::vector<array_1d<double, 3>>& rOutput,
		const ProcessInfo& rCurrentProcessInfo
	)
	{
		if (rOutput.size() != 1)
			rOutput.resize(1);

		if (rVariable == DISPLACEMENT) {
			const int& number_of_nodes = GetGeometry().size();
			Vector N = this->GetValue(SHAPE_FUNCTION_VALUES);

			for (SizeType i = 0; i < number_of_nodes; i++)
			{
				const NodeType& iNode = GetGeometry()[i];
				const array_1d<double, 3>& disp = iNode.FastGetSolutionStepValue(DISPLACEMENT);

				rOutput[0][0] += N[i] * disp[0];
				rOutput[0][1] += N[i] * disp[1];
				rOutput[0][2] += N[i] * disp[2];
			}
		}
	}
} // Namespace Kratos