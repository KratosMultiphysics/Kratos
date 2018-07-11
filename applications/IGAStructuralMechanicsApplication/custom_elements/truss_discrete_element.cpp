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
//                   Riccardo Rossi
//


// System includes
#include "includes/define.h"
#include "utilities/math_utils.h"
#include "geometries/geometry.h"

// External includes


// Project includes
#include "custom_elements/meshless_base_element.h"
#include "custom_elements/truss_discrete_element.h"

// Application includes
#include "iga_structural_mechanics_application.h"
#include "iga_structural_mechanics_application_variables.h"

namespace Kratos
{
	//************************************************************************************
	//************************************************************************************
	void TrussDiscreteElement::EquationIdVector(
		EquationIdVectorType& rResult,
		ProcessInfo& rCurrentProcessInfo
	)
	{
		KRATOS_TRY;
		unsigned int number_of_nodes = GetGeometry().size();
		unsigned int dim = number_of_nodes * 3;

		if (rResult.size() != dim)
			rResult.resize(dim);

		for (unsigned int i = 0; i < number_of_nodes; i++)
		{
			int index = i * 3;
			rResult[index]     = GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
			rResult[index + 1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();
			rResult[index + 2] = GetGeometry()[i].GetDof(DISPLACEMENT_Z).EquationId();
		}
		KRATOS_CATCH("")
	};

	//************************************************************************************
	//************************************************************************************
	void TrussDiscreteElement::GetDofList(
		DofsVectorType& rElementalDofList,
		ProcessInfo& rCurrentProcessInfo
	)
	{
		KRATOS_TRY;
		//if (rElementalDofList.size() != msElementSize) {
		//	rElementalDofList.resize(msElementSize);
		//}

		//for (int i = 0; i < msNumberOfNodes; ++i) {
		//	int index = i * msLocalSize;
		//	rElementalDofList[index] = this->GetGeometry()[i].pGetDof(DISPLACEMENT_X);
		//	rElementalDofList[index + 1] = this->GetGeometry()[i].pGetDof(DISPLACEMENT_Y);
		//	rElementalDofList[index + 2] = this->GetGeometry()[i].pGetDof(ROTATION_Z);
		//}

		rElementalDofList.resize(0);

		for (unsigned int i = 0; i < GetGeometry().size(); i++)
		{
			rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
			rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
			rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Z));
		}

		KRATOS_CATCH("")
	};
	//***********************************************************************************
	//***********************************************************************************
	void TrussDiscreteElement::FinalizeSolutionStep(
		ProcessInfo& rCurrentProcessInfo)
	{
		mConstitutiveLaw->FinalizeSolutionStep(GetProperties(),
			GetGeometry(),
			this->GetValue(SHAPE_FUNCTION_VALUES),
			rCurrentProcessInfo);
	}
	//************************************************************************************
	//************************************************************************************
	void TrussDiscreteElement::InitializeMaterial()
	{
		KRATOS_TRY

			if (GetProperties()[CONSTITUTIVE_LAW] != NULL)
			{
				//get shape functions for evaluation of stresses inside the constitutive law
				const Vector&  N = this->GetValue(SHAPE_FUNCTION_VALUES);
				const double integration_weight = this->GetValue(INTEGRATION_WEIGHT);

				mConstitutiveLaw = GetProperties()[CONSTITUTIVE_LAW]->Clone();
				ProcessInfo emptyProcessInfo = ProcessInfo();
				mConstitutiveLaw->InitializeMaterial(GetProperties(), GetGeometry(), N);
			}
			else
			{
				KRATOS_ERROR << "A constitutive law needs to be specified for the element with ID " << this->Id() << std::endl;
			}

		KRATOS_CATCH("");
	}
	//************************************************************************************
	//************************************************************************************
	void TrussDiscreteElement::CalculateAll(
		MatrixType& rLeftHandSideMatrix,
		VectorType& rRightHandSideVector,
		ProcessInfo& rCurrentProcessInfo,
		const bool CalculateStiffnessMatrixFlag,
		const bool CalculateResidualVectorFlag
	)
	{
		KRATOS_TRY
		// definition of problem size
		const unsigned int number_of_control_points = GetGeometry().size();
		unsigned int number_of_dofs = number_of_control_points * 3;

		//set up properties for Constitutive Law
		ConstitutiveLaw::Parameters Values(GetGeometry(), GetProperties(), rCurrentProcessInfo);

		Values.GetOptions().Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false);
		Values.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRESS);
		Values.GetOptions().Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

		//resizing as needed the LHS
		if (CalculateStiffnessMatrixFlag == true) //calculation of the matrix is required
		{
			if (rLeftHandSideMatrix.size1() != number_of_dofs)
				rLeftHandSideMatrix.resize(number_of_dofs, number_of_dofs);
			noalias(rLeftHandSideMatrix) = ZeroMatrix(number_of_dofs, number_of_dofs); //resetting LHS
		}
		//resizing as needed the RHS
		if (CalculateResidualVectorFlag == true) //calculation of the matrix is required
		{
			if (rRightHandSideVector.size() != number_of_dofs)
				rRightHandSideVector.resize(number_of_dofs);
			rRightHandSideVector = ZeroVector(number_of_dofs); //resetting RHS
		}

		//reading in of integration weight, shape function values and shape function derivatives
		const double& integration_weight = this->GetValue(INTEGRATION_WEIGHT);
		const Vector&   N     = this->GetValue(SHAPE_FUNCTION_VALUES);
		const Matrix&  DN_De  = this->GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);

		//Bending stabilization
		double E = GetProperties()[YOUNG_MODULUS];
		double area = GetProperties()[AREA];
		double prestress = GetProperties()[PRESTRESS_CAUCHY];

		Vector base_vector = ZeroVector(3);
		GetBaseVector(base_vector, DN_De);

		double a0 = norm_2(mBaseVector0);
		double a = norm_2(base_vector);

		//stresses
		double E11_membrane = 0.5 * (pow(a, 2) - pow(a0,2));   //Green Lagrange formulation (strain)
		//cfloat E11_b = (c_b - c_B);   //Green Lagrange formulation (strain)
		//cfloat E11_n = (c_n - c_N);   //Green Lagrange formulation (strain)

		double S11_membrane = prestress * area + E11_membrane * area*E / pow(a0,2);//normal force

		// 1st variation
		// variation of the axial strain 
		Vector epsilon_var_dof = ZeroVector(number_of_dofs);
		Get1stVariationsAxialStrain(epsilon_var_dof, base_vector, 3, DN_De);
		epsilon_var_dof = epsilon_var_dof / pow(a0, 2);

		// 2nd variation
		// variation of the axial strain 
		Matrix epsilon_var_2_dof = ZeroMatrix(number_of_dofs, number_of_dofs);
		Get2ndVariationsAxialStrain(epsilon_var_2_dof, 3, DN_De);
		epsilon_var_2_dof = epsilon_var_2_dof / pow(a0, 2);

		for (int r = 0; r<number_of_dofs; r++)
			for (int s = 0; s<number_of_dofs; s++)
			{
				rLeftHandSideMatrix(r, s) = E * area * epsilon_var_dof[r] * epsilon_var_dof[s] + S11_membrane * epsilon_var_2_dof(r, s);
			}

		rRightHandSideVector = -S11_membrane * epsilon_var_dof;

		rLeftHandSideMatrix = rLeftHandSideMatrix * integration_weight;
		rRightHandSideVector = rRightHandSideVector * integration_weight;

		KRATOS_CATCH("");
	}


	//************************************************************************************
	//************************************************************************************
	void TrussDiscreteElement::CalculateOnIntegrationPoints(
		const Variable<double>& rVariable,
		std::vector<double>& rOutput,
		const ProcessInfo& rCurrentProcessInfo
	)
	{
		if (rOutput.size() != 1)
			rOutput.resize(1);

		rOutput[0] = 0.0;// mConstitutiveLawVector[0]->GetValue(rVariable, rOutput[0]);
	}

	void TrussDiscreteElement::CalculateOnIntegrationPoints(
		const Variable<Vector>& rVariable,
		std::vector<Vector>& rValues,
		const ProcessInfo& rCurrentProcessInfo
	)
	{
		if (rValues.size() != 1)
		{
			rValues.resize(1);
		}

		rValues[0] = ZeroVector(3);
	}
} // Namespace Kratos


