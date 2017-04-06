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

#include "iga_structural_mechanics_application.h"
#include "iga_structural_mechanics_application_variables.h"

#include "utilities/math_utils.h"

#include "geometries/geometry.h"

namespace Kratos
{
//************************************************************************************
//************************************************************************************
MeshlessBaseElement::MeshlessBaseElement(
	IndexType NewId,
	GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!
}
//************************************************************************************
//************************************************************************************
MeshlessBaseElement::MeshlessBaseElement(
	IndexType NewId,
	GeometryType::Pointer pGeometry,  
	PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
{
}
//************************************************************************************
//************************************************************************************
Element::Pointer MeshlessBaseElement::Create(
	IndexType NewId, 
	NodesArrayType const& ThisNodes,  
	PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new MeshlessBaseElement(NewId, GetGeometry().Create(ThisNodes), pProperties));
}
//************************************************************************************
//************************************************************************************
MeshlessBaseElement::~MeshlessBaseElement()
{
}
//************************************************************************************
/*
* GetGeometryData Gets data stored on the element.
*
* @param integration_weight the integration weight is the specific 
  integration weight of the Gauss Point of this element.
* @param N shape function N(xi, eta) evaluated on location of the Gauss Point.
* @param DN_De shape function derivatives evaluated on location of this Gauss Point. 
  Generally first and second column for first order derivatives. Additionally if 
  needed third and fourth for second order derivatives and fith for cross terms.
*/
//************************************************************************************
void MeshlessBaseElement::GetGeometryData(double& integration_weight,
	Vector& N,
	Matrix& DN_De
	)
{
	integration_weight = this->GetValue(INTEGRATION_WEIGHT);
	N = this->GetValue(SHAPE_FUNCTION_VALUES);
	DN_De = this->GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);
}

//************************************************************************************
/*
* ComputeGlobalDerivatives Transform from local coordinates to
  global coordinates.
*
* @param DN_De derivatives of shape functions in two directions.
  Needed parameter for the calculation.
* @param DN_DX derivatives of shape function in global coordinates.
  return value. Has to be zero in advance.
* @param Jacobian calculated Jacobian. Is always of size 3x2, as
  the function only allows mapping from from 3D in Geometry Space to
  2D in Parameter Space.
*/
//************************************************************************************
void MeshlessBaseElement::ComputeGlobalDerivatives(const Matrix& DN_De,
	Matrix& DN_DX,
	Matrix& Jacobian)
{
	KRATOS_TRY
	const unsigned int number_of_points = GetGeometry().size();
	// working_space_dimension = number of dimensions in Geometry Space:
	const unsigned int working_space_dimension = 3;
	// local_space_dimension = number of dimensions in Parameter Space
	// If DN_De is used only for the first and second derivatives this term can be used:
	// local_space_dimension = GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES).size2();// GetGeometry().WorkingSpaceDimension();
	const unsigned int local_space_dimension = 2; 

	DN_DX.resize(number_of_points,working_space_dimension);

	this->Jacobian(DN_De, Jacobian);

	noalias(DN_DX) = prod(DN_De, Jacobian);
	KRATOS_CATCH("")
}
//************************************************************************************
/**
* Jacobian gives the mapping for the given shape functions. This 
  function calculates the mapping from 3D in Geometry Space to 
  2D in Parameter Space.
*
* @param DN_De derivatives of shape functions in two directions.
*
* @param Jacobian calculated Jacobian. Is always of size 3x2, as
  the function only allows mapping from from 3D in Geometry Space to 
  2D in Parameter Space.
*
*/
//************************************************************************************
void MeshlessBaseElement::Jacobian(const Matrix& DN_De,
	Matrix& Jacobian)
{
	const unsigned int number_of_points = GetGeometry().size();
	const unsigned int working_space_dimension = 3;// GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES).size2();// GetGeometry().WorkingSpaceDimension();
	const unsigned int local_space_dimension = 2;// GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES).size2();// GetGeometry().LocalSpaceDimension();

	Jacobian.resize(working_space_dimension, local_space_dimension);

	Jacobian.clear();
	for (unsigned int i = 0; i < number_of_points; i++)
	{
    std::cout << "Coords[" << i << "]" << (GetGeometry()[i]).Coordinates() << std::endl;

		for (unsigned int k = 0; k<working_space_dimension; k++)
		{
			for (unsigned int m = 0; m<local_space_dimension; m++)
			{
				Jacobian(k, m) += (GetGeometry()[i]).Coordinates()[k] * DN_De(i, m);
			}
		}
	}
}
//************************************************************************************
/**
* Hessian calculates the Hessian for the given system with the 
  shape function derivatives DN_De.
*
* @param DN_De derivatives of shape functions.
*
* @param Hessian calculated Hessian. Is always of size 3x3.
*
*/
//************************************************************************************
void MeshlessBaseElement::Hessian(Matrix& Hessian, const Matrix& DDN_DDe)
{
	const unsigned int number_of_points = GetGeometry().size();
	const unsigned int working_space_dimension = 3;// GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES).size2();// GetGeometry().WorkingSpaceDimension();

	Hessian.resize(working_space_dimension, working_space_dimension);
	Hessian = ZeroMatrix(working_space_dimension, working_space_dimension);

	for (size_t k = 0; k<number_of_points; k++)
	{
		const array_1d<double, 3> coords = GetGeometry()[k].Coordinates();

		Hessian(0, 0) += DDN_DDe(k, 0)*coords[0];
		Hessian(0, 1) += DDN_DDe(k, 1)*coords[0];
		Hessian(0, 2) += DDN_DDe(k, 2)*coords[0];

		Hessian(1, 0) += DDN_DDe(k, 0)*coords[1];
		Hessian(1, 1) += DDN_DDe(k, 1)*coords[1];
		Hessian(1, 2) += DDN_DDe(k, 2)*coords[1];

		Hessian(2, 0) += DDN_DDe(k, 0)*coords[2];
		Hessian(2, 1) += DDN_DDe(k, 1)*coords[2];
		Hessian(2, 2) += DDN_DDe(k, 2)*coords[2];
	}
}

//***********************************************************************************
/**
* CrossProduct calculates the cross product of two 3d vectors. a x b = cross
*/
//***********************************************************************************
void MeshlessBaseElement::CrossProduct(
	array_1d<double, 3>& cross,
	const array_1d<double, 3>& a,
	const array_1d<double, 3>& b)
{
	cross[0] = a[1] * b[2] - a[2] * b[1];
	cross[1] = a[2] * b[0] - a[0] * b[2];
	cross[2] = a[0] * b[1] - a[1] * b[0];
}
//***********************************************************************************
/**
* CrossProduct calculates the cross product of two 3d vectors. a x b.
* @return cross = a x b
*/
//***********************************************************************************
array_1d<double, 3> MeshlessBaseElement::CrossProduct(
	const array_1d<double, 3>& a,
	const array_1d<double, 3>& b)
{
	array_1d<double, 3> cross;
	cross[0] = a[1] * b[2] - a[2] * b[1];
	cross[1] = a[2] * b[0] - a[0] * b[2];
	cross[2] = a[0] * b[1] - a[1] * b[0];
	return cross;
}
//***********************************************************************************
//***********************************************************************************
} // Namespace Kratos


