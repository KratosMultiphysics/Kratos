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
#include "custom_conditions/meshless_base_condition.h"
#include "utilities/math_utils.h"

#include "iga_structural_mechanics_application_variables.h"
#include "iga_structural_mechanics_application.h"

namespace Kratos
{

//************************************************************************************
//************************************************************************************
MeshlessBaseCondition::MeshlessBaseCondition(IndexType NewId, GeometryType::Pointer pGeometry)
    : Condition(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!
}
//************************************************************************************
//************************************************************************************
MeshlessBaseCondition::MeshlessBaseCondition(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
    : Condition(NewId, pGeometry, pProperties)
{
}
//************************************************************************************
//************************************************************************************
Condition::Pointer MeshlessBaseCondition::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new MeshlessBaseCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
}
//************************************************************************************
// Destructor
//************************************************************************************
MeshlessBaseCondition::~MeshlessBaseCondition()
{
}

/**
* Jacobian gives the mapping for the given shape functions. This
  function calculates the mapping from 3D in Geometry Space to
  2D in Parameter Space.
*
* @param DN_De derivatives of shape functions in two directions.
* @param Jacobian calculated Jacobian. Is always of size 3x2, as
the function only allows mapping from from 3D in Geometry Space to
2D in Parameter Space.
*
*/
void MeshlessBaseCondition::Jacobian(const Matrix& DN_De,
	Matrix& Jacobian)
{
	const unsigned int number_of_points = GetGeometry().size();
	const unsigned int working_space_dimension = 3;// GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES).size2();// GetGeometry().WorkingSpaceDimension();
	const unsigned int local_space_dimension = 2;// GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES).size2();// GetGeometry().LocalSpaceDimension();

	Jacobian.resize(working_space_dimension, local_space_dimension);

	Jacobian.clear();
	for (unsigned int i = 0; i < number_of_points; i++)
	{
		for (unsigned int k = 0; k<working_space_dimension; k++)
		{
			for (unsigned int m = 0; m<local_space_dimension; m++)
			{
				Jacobian(k, m) += (GetGeometry()[i]).Coordinates()[k] * DN_De(i, m);
			}
		}
	}
}

/**
* JacobianInitial gives the mapping for the given shape functions
  for the undeformed system. This function calculates the mapping 
  from 3D in Geometry Space to 2D in Parameter Space.
*
* @param DN_De derivatives of shape functions.
* @param Jacobian calculated Jacobian. Is always of size 3x2, as
  the function only allows mapping from from 3D in Geometry Space to
  2D in Parameter Space.
*/
void MeshlessBaseCondition::JacobianInitial(const Matrix& DN_De,
	Matrix& Jacobian)
{
	const unsigned int number_of_points = GetGeometry().size();
	const unsigned int working_space_dimension = 3;// GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES).size2();// GetGeometry().WorkingSpaceDimension();
	const unsigned int local_space_dimension = 2;// GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES).size2();// GetGeometry().LocalSpaceDimension();

	Jacobian.resize(working_space_dimension, local_space_dimension);
	Jacobian.clear();

	for (unsigned int i = 0; i < number_of_points; i++)
	{
		//for (unsigned int k = 0; k<working_space_dimension; k++)
		//{
			for (unsigned int m = 0; m<local_space_dimension; m++)
			{
				Jacobian(0, m) += (GetGeometry()[i]).X0() * DN_De(i, m);
				Jacobian(1, m) += (GetGeometry()[i]).Y0() * DN_De(i, m);
				Jacobian(2, m) += (GetGeometry()[i]).Z0() * DN_De(i, m);
			}
		//}
	}
}



/**
* MappingGeometricToParameter calculates the J-tilde for the mapping from 
  Geometric to Parameter Space. This paramater is needed for all condition 
  integrations on edges.
*
* @param DN_De derivatives of shape functions.
* @param JGeometricToParameter Mapping parameter for Geometric Space to 
  Parameter Space
* 
* @see Jacobian
*/
void MeshlessBaseCondition::MappingGeometricToParameter(const Matrix& DN_De,
	const array_1d<double, 2>& Tangents,
	double& JGeometricToParameter)
{
	Matrix J;
	Jacobian(DN_De, J);

	//basis vectors g1 and g2
	array_1d<double, 3> g1;
	array_1d<double, 3> g2;

	g1[0] = J(0, 0);
	g2[0] = J(0, 1);
	g1[1] = J(1, 0);
	g2[1] = J(1, 1);
	g1[2] = J(2, 0);
	g2[2] = J(2, 1);
	//basis vector g3
	//CrossProduct(g3, g1, g2);

	array_1d<double, 3> temp = g1*Tangents[0] + g2*Tangents[1];
	// g1*localDerivatives[0] + g2*localDerivatives[1];//(g1 / norm_2(g1))*localDerivatives[0] + (g2 / norm_2(g2))*localDerivatives[1];
	JGeometricToParameter = norm_2(temp);
}

/**
* GetBasisVectors calculates the basis vectors g1, g2 and the normalized basis g3
* @param DN_De first order derivatives of shape functions.
* @param g1, g2 basis vectors
* @param g3 NORMALIZED basis vectors
*/
void MeshlessBaseCondition::GetBasisVectors(
	const Matrix& DN_De,
	array_1d<double, 3>& g1,
	array_1d<double, 3>& g2,
	array_1d<double, 3>& g3)
{
	Matrix J;
	Jacobian(DN_De, J);

	g1[0] = J(0, 0);
	g2[0] = J(0, 1);
	g1[1] = J(1, 0);
	g2[1] = J(1, 1);
	g1[2] = J(2, 0);
	g2[2] = J(2, 1);

	//basis vector g3
	CrossProduct(g3, g1, g2);
	g3 = g3 / norm_2(g3);
}

/**
* GetInitialBasisVectors calculates the basis vectors g1, g2 and the normalized 
  basis g3 in the UNDEFORMED system.
* @param DN_De first order derivatives of shape functions.
* @param g10, g20 basis vectors
* @param g30 NORMALIZED basis vectors
*/
void MeshlessBaseCondition::GetInitialBasisVectors(
	const Matrix& DN_De,
	array_1d<double, 3>& g10,
	array_1d<double, 3>& g20,
	array_1d<double, 3>& g30)
{
	Matrix J;
	JacobianInitial(DN_De, J);

	g10[0] = J(0, 0);
	g20[0] = J(0, 1);
	g10[1] = J(1, 0);
	g20[1] = J(1, 1);
	g10[2] = J(2, 0);
	g20[2] = J(2, 1);

	//basis vector g30
	CrossProduct(g30, g10, g20);
	g30 = g30 / norm_2(g30);
}

/**
* CrossProduct calculates the cross product of two 3d vectors. a x b = cross
*/
void MeshlessBaseCondition::CrossProduct(
	array_1d<double, 3>& cross,
	const array_1d<double, 3>& a,
	const array_1d<double, 3>& b)
{
	cross[0] = a[1] * b[2] - a[2] * b[1];
	cross[1] = a[2] * b[0] - a[0] * b[2];
	cross[2] = a[0] * b[1] - a[1] * b[0];
}

/**
* CrossProduct calculates the cross product of two 3d vectors. a x b.
* @return cross = a x b
*/
array_1d<double, 3> MeshlessBaseCondition::CrossProduct(
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


