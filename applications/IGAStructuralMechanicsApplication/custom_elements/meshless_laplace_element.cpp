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
#include "custom_elements/meshless_base_element.h"
#include "custom_elements/meshless_laplace_element.h"

#include "iga_structural_mechanics_application.h"
#include "iga_structural_mechanics_application_variables.h"

#include "utilities/math_utils.h"

//#include "geometries/geometry.h"
//#include "custom_geometries/meshless_geometry.h"

namespace Kratos
{

//************************************************************************************
//************************************************************************************
MeshlessLaplaceElement::MeshlessLaplaceElement(
	IndexType NewId,
	GeometryType::Pointer pGeometry)
    : MeshlessBaseElement(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!

}

//************************************************************************************
//************************************************************************************
MeshlessLaplaceElement::MeshlessLaplaceElement(
	IndexType NewId,
	GeometryType::Pointer pGeometry,  
	PropertiesType::Pointer pProperties)
    : MeshlessBaseElement(NewId, pGeometry, pProperties)
{
}

Element::Pointer MeshlessLaplaceElement::Create(
	IndexType NewId, 
	NodesArrayType const& ThisNodes,  
	PropertiesType::Pointer pProperties) const
{
    return MeshlessBaseElement::Pointer(new MeshlessLaplaceElement(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

MeshlessLaplaceElement::~MeshlessLaplaceElement()
{
}

//************************************************************************************
//************************************************************************************
void MeshlessLaplaceElement::GetDofList(
	DofsVectorType& ElementalDofList,
	ProcessInfo& rCurrentProcessInfo)

{
	unsigned int number_of_nodes = GetGeometry().PointsNumber();
	if (ElementalDofList.size() != number_of_nodes)
		ElementalDofList.resize(number_of_nodes);
	for (unsigned int i = 0; i<number_of_nodes; i++)
	{
		ElementalDofList[i] = GetGeometry()[i].pGetDof(TEMPERATURE);
	}
}
//************************************************************************************
//************************************************************************************
void MeshlessLaplaceElement::EquationIdVector(
	EquationIdVectorType& rResult,
	ProcessInfo& rCurrentProcessInfo)

{
	KRATOS_TRY
	unsigned int number_of_nodes = GetGeometry().PointsNumber();
	if (rResult.size() != number_of_nodes)
		rResult.resize(number_of_nodes);
	for (unsigned int i = 0; i<number_of_nodes; i++)
	{
		rResult[i] = GetGeometry()[i].GetDof(TEMPERATURE).EquationId();
	}

	KRATOS_CATCH("")
}
//************************************************************************************
//************************************************************************************
void MeshlessLaplaceElement::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
	KRATOS_TRY
	const unsigned int number_of_points = GetGeometry().size();
	const unsigned int dim = GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES).size2();//GetGeometry()[0].WorkingSpaceDimension();
	//resizing as needed the LHS
	if (rLeftHandSideMatrix.size1() != number_of_points)
		rLeftHandSideMatrix.resize(number_of_points, number_of_points, false);
	noalias(rLeftHandSideMatrix) = ZeroMatrix(number_of_points, number_of_points); //resetting LHS


																				   //resizing as needed the RHS
	if (rRightHandSideVector.size() != number_of_points)
		rRightHandSideVector.resize(number_of_points, false);
	rRightHandSideVector = ZeroVector(number_of_points); //resetting RHS

														 //reading integration points and local gradients
	double integration_weight = this->GetValue(INTEGRATION_WEIGHT);
	Vector N = this->GetValue(SHAPE_FUNCTION_VALUES);
	Matrix DN_De = this->GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);
	

	Matrix J0;
	Matrix DN_DX(number_of_points, dim);
	

	ComputeGlobalDerivatives(DN_De, DN_DX, J0);

	Matrix InvJ0(dim, dim);
	double DetJ0;
	Vector temp(number_of_points);


	//calculating inverse jacobian and jacobian determinant
	MathUtils<double>::InvertMatrix(J0, InvJ0, DetJ0);


	double IntToReferenceWeight = integration_weight * DetJ0;
	noalias(rLeftHandSideMatrix) += IntToReferenceWeight * prod(DN_DX, trans(DN_DX)); //
//calculating external forces
	noalias(rRightHandSideVector) = ZeroVector(number_of_points); //case of zero ext forces

																  // RHS = ExtForces - K*temp;
	for (unsigned int i = 0; i < number_of_points; i++)
	{
		temp[i] = GetGeometry()[i].GetSolutionStepValue(TEMPERATURE); //this includes the - sign
	}
																	  //axpy_prod(rLeftHandSideMatrix, temp, rRightHandSideVector, false);  //RHS -= K*temp
	noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, temp);

	KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************
void MeshlessLaplaceElement::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
	MatrixType temp(0, 0);
	CalculateLocalSystem(temp, rRightHandSideVector, rCurrentProcessInfo);
}

//************************************************************************************
//************************************************************************************
} // Namespace Kratos


