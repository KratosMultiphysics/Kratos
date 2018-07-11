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
#include "custom_conditions/meshless_force_interface_condition.h"
#include "utilities/math_utils.h"
#include "includes/define.h"

#include "iga_structural_mechanics_application_variables.h"
#include "iga_structural_mechanics_application.h"


namespace Kratos
{
	//************************************************************************************
	//************************************************************************************
	MeshlessForceInterfaceCondition::MeshlessForceInterfaceCondition(IndexType NewId, GeometryType::Pointer pGeometry)
		: MeshlessBaseCondition(NewId, pGeometry)
	{
		//DO NOT ADD DOFS HERE!!!
	}


	//************************************************************************************
	//************************************************************************************
	MeshlessForceInterfaceCondition::MeshlessForceInterfaceCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
		: MeshlessBaseCondition(NewId, pGeometry, pProperties)
	{
	}


	Condition::Pointer MeshlessForceInterfaceCondition::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
	{
		return MeshlessBaseCondition::Pointer(new MeshlessForceInterfaceCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
	}


	// Destructor
	MeshlessForceInterfaceCondition::~MeshlessForceInterfaceCondition()
	{
	}

	//************************************************************************************
	//************************************************************************************
	//void MeshlessForceInterfaceCondition::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	//{
	//	MatrixType temp(0, 0);
	//	CalculateLocalSystem(temp, rRightHandSideVector, rCurrentProcessInfo);
	//}

	void MeshlessForceInterfaceCondition::CalculateRightHandSide(VectorType& rRightHandSideVector,
		ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY
		MatrixType temp(0, 0);
		CalculateLocalSystem(temp, rRightHandSideVector, rCurrentProcessInfo);
		KRATOS_CATCH("")
	}

	//************************************************************************************
	//************************************************************************************
	/**
	* CalculateLocalSystem
	* Returns only rRightHandSideVector as there is no impact on the 
	* siffness due to loads
	*/
	void MeshlessForceInterfaceCondition::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY
		//std::cout << "start MeshlessForceInterfaceCondition" << std::endl;
		//Check(rCurrentProcessInfo);
		const unsigned int number_of_points = GetGeometry().size();

		if (rLeftHandSideMatrix.size1() != number_of_points * 3)
			rLeftHandSideMatrix.resize(number_of_points * 3, number_of_points * 3, false);
		noalias(rLeftHandSideMatrix) = ZeroMatrix(number_of_points * 3, number_of_points * 3); //resetting LHS

		if (rRightHandSideVector.size() != number_of_points * 3)
			rRightHandSideVector.resize(number_of_points * 3, false);
		rRightHandSideVector = ZeroVector(number_of_points * 3); //resetting RHS

		//std::cout << "here note" << std::endl;
		const Vector& N = this->GetValue(SHAPE_FUNCTION_VALUES);
		Vector& force_vector = this->GetValue(EXTERNAL_FORCES_VECTOR);
		//Vector force_vector = ZeroVector(3);
		//force_vector(2) = 1000;
		//KRATOS_WATCH(N)
		KRATOS_WATCH(force_vector)
		Vector fLoads = ZeroVector(number_of_points * 3);

		for (unsigned int i = 0; i < number_of_points; i++)
		{
			int index = 3 * i;
			fLoads[index]     -= force_vector[0] * N[i];
			fLoads[index + 1] -= force_vector[1] * N[i];
			fLoads[index + 2] -= force_vector[2] * N[i];
		}
		noalias(rRightHandSideVector) -= fLoads;
		//force_vector(0) = 0;
		//force_vector(1) = 0;
		//force_vector(2) = 0;
		//this->SetValue(EXTERNAL_FORCES_VECTOR, force_vector);
		std::vector<Vector> coordinates;
		GetValueOnIntegrationPoints(COORDINATES, coordinates, rCurrentProcessInfo);
		//KRATOS_WATCH(coordinates[0])
		//KRATOS_WATCH(rRightHandSideVector)
		//for (unsigned int i = 0; i < GetGeometry().size(); i++)
		//{
		//	std::cout << GetGeometry()[i].Id() << std::endl;
		//}

		//Vector condition_coords = ZeroVector(3);
		//for (SizeType i = 0; i < number_of_points; i++)
		//{
		//	NodeType & iNode = GetGeometry()[i];

		//	Vector part_external_force = N[i] * force_vector;
		//	Vector full = part_external_force + iNode.GetValue(EXTERNAL_FORCES_VECTOR);
		//	iNode.SetValue(EXTERNAL_FORCES_VECTOR, full);
		//}

		KRATOS_CATCH("")
	}
	/***********************************************************************************/
	/***********************************************************************************/

	int  MeshlessForceInterfaceCondition::Check(const ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY;
		if (!Has(SHAPE_FUNCTION_VALUES))
			KRATOS_ERROR << "No SHAPE_FUNCTION_VALUES assigned!" << std::endl;
		if (!this->Has(EXTERNAL_FORCES_VECTOR))
			KRATOS_ERROR << "EXTERNAL_FORCES_VECTOR not assigned!" << std::endl;
		KRATOS_CATCH("")
	}

	//***********************************************************************************
	//***********************************************************************************

	void MeshlessForceInterfaceCondition::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
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
	void MeshlessForceInterfaceCondition::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo)
	{

		ElementalDofList.resize(0);

		for (unsigned int i = 0; i < GetGeometry().size(); i++)
		{
			ElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
			ElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
			ElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Z));
		}
	}

	/***********************************************************************************/
	/***********************************************************************************/
	void MeshlessForceInterfaceCondition::GetValueOnIntegrationPoints(
		const Variable<Vector>& rVariable,
		std::vector<Vector>& rValues,
		const ProcessInfo& rCurrentProcessInfo
	)
	{
		const unsigned int size = GetGeometry().IntegrationPoints().size();

		if (rValues.size() != 1)
			rValues.resize(1);

		if (rVariable == COORDINATES) {
			const int& number_of_control_points = GetGeometry().size();
			Vector N = this->GetValue(SHAPE_FUNCTION_VALUES);
			Vector condition_coords = ZeroVector(3);
			for (SizeType i = 0; i < number_of_control_points; i++)
			{
				const NodeType & iNode = GetGeometry()[i];
				const array_1d<double, 3>& coords = iNode.Coordinates();

				condition_coords[0] += N[i] * coords[0];
				condition_coords[1] += N[i] * coords[1];
				condition_coords[2] += N[i] * coords[2];
			}
			rValues[0] = condition_coords;
		}
		else {
			CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
		}
	}
    /***********************************************************************************/
    /***********************************************************************************/
    void MeshlessForceInterfaceCondition::GetValueOnIntegrationPoints(
        const Variable<array_1d<double, 3>>& rVariable,
        std::vector<array_1d<double, 3>>& rValues,
        const ProcessInfo& rCurrentProcessInfo
    )
    {
        if (rValues.size() != 1)
            rValues.resize(1);

        if (rVariable == VELOCITY) {
            const int& number_of_control_points = GetGeometry().size();
            Vector N = this->GetValue(SHAPE_FUNCTION_VALUES);

            array_1d<double,3> velocity = ZeroVector(3);
            for (SizeType i = 0; i < number_of_control_points; i++)
            {
                const NodeType & iNode = GetGeometry()[i];
                const array_1d<double, 3>& vel = iNode.FastGetSolutionStepValue(VELOCITY, 0);

                velocity[0] += N[i] * vel[0];
                velocity[1] += N[i] * vel[1];
                velocity[2] += N[i] * vel[2];
            }
            rValues[0] = velocity;
        }
        else {
            CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
        }
    }

	////************************************************************************************
	////************************************************************************************
	//void MeshlessForceInterfaceCondition::GetValuesVector(Vector& values, int Step)
	//{
	//	if (values.size() != 3)
	//		values.resize(3, false);

	//	const int& number_of_nodes = GetGeometry().size();
	//	Vector N = this->GetValue(SHAPE_FUNCTION_VALUES);

	//	for (SizeType i = 0; i < number_of_nodes; i++)
	//	{
	//		const NodeType & iNode = GetGeometry()[i];
	//		const array_1d<double, 3>& disp = iNode.FastGetSolutionStepValue(DISPLACEMENT, Step);

	//		values[0] += N[i] * disp[0];
	//		values[1] += N[i] * disp[1];
	//		values[2] += N[i] * disp[2];
	//	}
	//}

	////************************************************************************************
	////************************************************************************************
	//void MeshlessForceInterfaceCondition::GetFirstDerivativesVector(Vector& values, int Step)
	//{
	//	if (values.size() != 3)
	//		values.resize(3, false);


	//	const int& number_of_nodes = GetGeometry().size();
	//	Vector N = this->GetValue(SHAPE_FUNCTION_VALUES);

	//	for (SizeType i = 0; i < number_of_nodes; i++)
	//	{
	//		const NodeType & iNode = GetGeometry()[i];
	//		const array_1d<double, 3>& vel = iNode.FastGetSolutionStepValue(VELOCITY, Step);

	//		values[0] += N[i] * vel[0];
	//		values[1] += N[i] * vel[1];
	//		values[2] += N[i] * vel[2];
	//	}
	//}

	

} // Namespace Kratos