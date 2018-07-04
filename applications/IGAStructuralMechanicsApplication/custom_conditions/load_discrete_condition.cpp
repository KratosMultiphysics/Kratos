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
	LoadDiscreteCondition::LoadDiscreteCondition(IndexType NewId, GeometryType::Pointer pGeometry)
		: MeshlessBaseCondition(NewId, pGeometry)
	{
	}


	//************************************************************************************
	//************************************************************************************
	LoadDiscreteCondition::LoadDiscreteCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
		: MeshlessBaseCondition(NewId, pGeometry, pProperties)
	{
	}


	Condition::Pointer LoadDiscreteCondition::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
	{
		return MeshlessBaseCondition::Pointer(new LoadDiscreteCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
	}


	// Destructor
	LoadDiscreteCondition::~LoadDiscreteCondition()
	{
	}

	//************************************************************************************
	//************************************************************************************
	void LoadDiscreteCondition::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
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
	*/
	void LoadDiscreteCondition::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY

		const unsigned int number_of_control_points = GetGeometry().size();
		const unsinged int number_of_dofs = number_of_control_points * 3;

		if (rLeftHandSideMatrix.size1() != number_of_dofs)
			rLeftHandSideMatrix.resize(number_of_dofs, number_of_dofs, false);
		noalias(rLeftHandSideMatrix) = ZeroMatrix(number_of_dofs, number_of_dofs); //resetting LHS

		if (rRightHandSideVector.size() != number_of_dofs)
			rRightHandSideVector.resize(number_of_dofs, false);
		rRightHandSideVector = ZeroVector(number_of_dofs); //resetting RHS

		if (!Has(SHAPE_FUNCTION_VALUES))
			KRATOS_ERROR << "No SHAPE_FUNCTION_VALUES assigned!" << std::endl;
		Vector N = this->GetValue(SHAPE_FUNCTION_VALUES);

		// Edge and Surface loads
		if (this->Has(INTEGRATION_WEIGHT))
		{
			double integration_weight = this->GetValue(INTEGRATION_WEIGHT);

			if (!Has(SHAPE_FUNCTION_LOCAL_DERIVATIVES))
				KRATOS_ERROR << "No SHAPE_FUNCTION_LOCAL_DERIVATIVES assigned!" << std::endl;
			Matrix DN_De = this->GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);

			array_1d<double, 3> g1, g2, g3;
			GetBasisVectors(DN_De, g1, g2, g3); //g3 normalized!

			// Surface loads
			if (this->Has(SURFACE_LOAD))
			{
				const array_1d<double, 3 > surface_load = this->GetValue(SURFACE_LOAD);
				CrossProduct(g3, g1, g2);
				double dArea = norm_2(g3);
				double integration_weight_area = integration_weight * dArea;

				for (unsigned int i = 0; i < number_of_points; i++)
				{
					int index = 3 * i;
					fLoads[index]     = -surface_load[0] * integration_weight_area * N[i];
					fLoads[index + 1] = -surface_load[1] * integration_weight_area * N[i];
					fLoads[index + 2] = -surface_load[2] * integration_weight_area * N[i];
				}
			}

			// Pressure loads
			if (this->Has(PRESSURE))
			{
				double pressure = this->GetValue(PRESSURE);
				CrossProduct(g3, g1, g2);
				double dArea = norm_2(g3);

				// Check wether applied to edge or to surface.
				double integration_weight_area = 0.0;
				if (!this->Has(TANGENTS))
					integration_weight_area = integration_weight * dArea;
				else
				{
					integration_weight_area = integration_weight * norm_2(g1 * Tangents[0] + g2 * Tangents[1]);
				}

				array_1d<double, 3> direction = g3 / dArea;

				for (unsigned int i = 0; i < number_of_points; i++)
				{
					int index = 3 * i;
					fLoads[index]     = -direction[0] * pressure * integration_weight_area * N[i];
					fLoads[index + 1] = -direction[1] * pressure * integration_weight_area * N[i];
					fLoads[index + 2] = -direction[2] * pressure * integration_weight_area * N[i];
				}
			}

			// Edge loads
			if (this->Has(LINE_LOAD))
			{
				array_1d<double, 3> line_load = this->GetValue(LINE_LOAD);

				double integration_weight_area = integration_weight * norm_2(g1 * Tangents[0] + g2 * Tangents[1]);

				for (unsigned int i = 0; i < number_of_points; i++)
				{
					int index = 3 * i;
					fLoads[index]     = -line_load[0] * integration_weight_area * N[i];
					fLoads[index + 1] = -line_load[1] * integration_weight_area * N[i];
					fLoads[index + 2] = -line_load[2] * integration_weight_area * N[i];
				}
			}
		}

		// Point loads
		if (this->Has(POINT_LOAD))
		{
			const array_1d<double, 3 > PointLoad = this->GetValue(POINT_LOAD);

			for (unsigned int i = 0; i < number_of_points; i++)
			{
				int index = 3 * i;
				fLoads[index]     = -PointLoad[0] * N[i];
				fLoads[index + 1] = -PointLoad[1] * N[i];
				fLoads[index + 2] = -PointLoad[2] * N[i];
			}
		}

		noalias(rRightHandSideVector) -= fLoads;
		KRATOS_CATCH("")
	}


	//***********************************************************************************
	//***********************************************************************************

	void LoadDiscreteCondition::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
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
	void LoadDiscreteCondition::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo)
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