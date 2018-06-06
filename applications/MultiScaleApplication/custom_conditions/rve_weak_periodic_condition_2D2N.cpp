/*
==============================================================================
KratosMultiScaleApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
janosch.stascheit@rub.de
nagel@sd.rub.de
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
- Ruhr-University Bochum, Institute for Structural Mechanics, Germany


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/
//
//   Project Name:        Kratos
//   Last Modified by:    $Author: Massimo Petracca $
//   Date:                $Date: 2013-11-05 19:00:00 $
//   Revision:            $Revision: 1.00 $
//
//


// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "rve_weak_periodic_condition_2D2N.h"
#include "utilities/math_utils.h"
#include "geometries/point_3d.h"
#include "multiscale_application_variables.h"
#include "custom_utilities/math_helpers.h"

namespace Kratos
{

RveWeakPeriodicCondition2D2N::RveWeakPeriodicCondition2D2N(IndexType NewId, GeometryType::Pointer pGeometry)
	: MyBase(NewId, pGeometry)
{
	m_is_skew_symmetric_constraint = false;
}

RveWeakPeriodicCondition2D2N::RveWeakPeriodicCondition2D2N(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
	: MyBase(NewId, pGeometry, pProperties)
{
	m_is_skew_symmetric_constraint = false;
}

RveWeakPeriodicCondition2D2N::RveWeakPeriodicCondition2D2N(const RveWeakPeriodicCondition2D2N& rOther)
	: MyBase(rOther)
{
	m_is_skew_symmetric_constraint = rOther.IsSkewSymmetricConstraint();
}

Condition::Pointer RveWeakPeriodicCondition2D2N::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
{
	/*return Condition::Pointer(new RveWeakPeriodicCondition2D2N(NewId, GetGeometry().Create(ThisNodes), pProperties));*/
	RveWeakPeriodicCondition2D2N::Pointer copy_cnd(new RveWeakPeriodicCondition2D2N(NewId, GetGeometry().Create(ThisNodes), pProperties));
	copy_cnd->IsSkewSymmetricConstraint() = true;
	return copy_cnd;
}

RveWeakPeriodicCondition2D2N::~RveWeakPeriodicCondition2D2N()
{
}

void RveWeakPeriodicCondition2D2N::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
	//{
	//	if(rLeftHandSideMatrix.size1() != 7 || rLeftHandSideMatrix.size2() != 7) rLeftHandSideMatrix.resize(7, 7,false);
	//	noalias(rLeftHandSideMatrix) = ZeroMatrix(7, 7);

	//	if(rRightHandSideVector.size() != 7) rRightHandSideVector.resize(7, false);
	//	noalias( rRightHandSideVector ) = ZeroVector(7);

	//	GeometryType& geom = GetGeometry();

	//	// get current values

	//	Vector currentValues(7,0.0);
	//	//currentValues(0) = geom[0].FastGetSolutionStepValue(DISPLACEMENT_X); // node 1 UX
	//	//currentValues(1) = geom[0].FastGetSolutionStepValue(DISPLACEMENT_Y); // node 1 UY
	//	//currentValues(2) = geom[1].FastGetSolutionStepValue(DISPLACEMENT_X); // node 2 UX
	//	//currentValues(3) = geom[1].FastGetSolutionStepValue(DISPLACEMENT_Y); // node 2 UY
	//	currentValues(4) = geom[2].FastGetSolutionStepValue(RVE_WPC_LAGRANGIAN_DOF_X); // lagrangian node - DOF 1
	//	currentValues(5) = geom[2].FastGetSolutionStepValue(RVE_WPC_LAGRANGIAN_DOF_Y); // lagrangian node - DOF 2
	//	currentValues(6) = geom[2].FastGetSolutionStepValue(RVE_WPC_LAGRANGIAN_DOF_Z); // lagrangian node - DOF 3

	//	// compute the outward normal vector
	//	// Note: 2D on XY plane is assumed. boundary nodes in CCW order assumed

	//	double x0 = geom[0].X0();
	//	double y0 = geom[0].Y0();
	//	double x1 = geom[1].X0();
	//	double y1 = geom[1].Y0();

	//	double nxL = (y1 - y0)/2.0;
	//	double nyL = (x0 - x1)/2.0;

	//	// form the stiffness matrix

	//	Matrix& K = rLeftHandSideMatrix;

	//	K(4,0) = nxL;                 K(4,2) = nxL;
	//				   K(5,1) = nyL;                 K(5,3) = nyL;
	//	K(6,0) = nyL;  K(6,1) = nxL;  K(6,2) = nyL;  K(6,3) = nxL;

	//	K(0,4) = nxL;                 K(2,4) = nxL;
	//				   K(1,5) = nyL;                 K(3,5) = nyL;
	//	K(0,6) = nyL;  K(1,6) = nxL;  K(2,6) = nyL;  K(3,6) = nxL;

	//	// form residual

	//	noalias(rRightHandSideVector) -= prod( rLeftHandSideMatrix, currentValues );
	//	return;
	//}

	// get the number of required lagrangian dofs

	/*array_1d<unsigned int, 3> mask;
	unsigned int n_lag_dofs = CalculateLagrangianDofMask(mask);
	unsigned int n_dofs = 4 + n_lag_dofs;*/
	array_1d<unsigned int, 3> mask;
	unsigned int n_dofs;
	if(m_is_skew_symmetric_constraint) {
		n_dofs = 5;
	}
	else {
		unsigned int n_lag_dofs = CalculateLagrangianDofMask(mask);
		n_dofs = 4 + n_lag_dofs;
	}

	// resize system matrix and vector

	if(rLeftHandSideMatrix.size1() != n_dofs || rLeftHandSideMatrix.size2() != n_dofs) rLeftHandSideMatrix.resize(n_dofs, n_dofs,false);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(n_dofs, n_dofs);

    if(rRightHandSideVector.size() != n_dofs) rRightHandSideVector.resize(n_dofs, false);
	noalias( rRightHandSideVector ) = ZeroVector(n_dofs);

	GeometryType& geom = GetGeometry();

	// get current values

	Vector currentValues(n_dofs,0.0);
	//currentValues(0) = geom[0].FastGetSolutionStepValue(DISPLACEMENT_X); // node 1 UX
	//currentValues(1) = geom[0].FastGetSolutionStepValue(DISPLACEMENT_Y); // node 1 UY
	//currentValues(2) = geom[1].FastGetSolutionStepValue(DISPLACEMENT_X); // node 2 UX
	//currentValues(3) = geom[1].FastGetSolutionStepValue(DISPLACEMENT_Y); // node 2 UY
	if(m_is_skew_symmetric_constraint) {
		currentValues(4) = geom[2].FastGetSolutionStepValue(RVE_WPR_LAGRANGIAN_DOF); // lagrangian node - DOF 1
	}
	else {
		if(mask[0]>0)currentValues(mask[0]) = geom[2].FastGetSolutionStepValue(RVE_WPC_LAGRANGIAN_DOF_X); // lagrangian node - DOF 1
		if(mask[1]>0)currentValues(mask[1]) = geom[2].FastGetSolutionStepValue(RVE_WPC_LAGRANGIAN_DOF_Y); // lagrangian node - DOF 2
		if(mask[2]>0)currentValues(mask[2]) = geom[2].FastGetSolutionStepValue(RVE_WPC_LAGRANGIAN_DOF_Z); // lagrangian node - DOF 3
	}

	// compute the outward normal vector
	// Note: 2D on XY plane is assumed. boundary nodes in CCW order assumed

	double x0 = geom[0].X0();
	double y0 = geom[0].Y0();
	double x1 = geom[1].X0();
	double y1 = geom[1].Y0();

	double nxL = (y1 - y0)/2.0;
	double nyL = (x0 - x1)/2.0;

	// form the stiffness matrix

	Matrix& K = rLeftHandSideMatrix;


	if(m_is_skew_symmetric_constraint) {
		K(4,0) = -0.5*nyL;  K(4,1) = 0.5*nxL;  K(4,2) = -0.5*nyL;  K(4,3) = 0.5*nxL;
		K(0,4) = -0.5*nyL;  K(1,4) = 0.5*nxL;  K(2,4) = -0.5*nyL;  K(3,4) = 0.5*nxL;
	}
	else {
		if(mask[0]>0) {
			K(mask[0],0) = nxL;
			K(mask[0],2) = nxL;
			K(0,mask[0]) = nxL;
			K(2,mask[0]) = nxL;
		}
		if(mask[1]>0) {
			K(mask[1],1) = nyL;
			K(mask[1],3) = nyL;
			K(1,mask[1]) = nyL;
			K(3,mask[1]) = nyL;
		}
		if(mask[2]>0) {
			K(mask[2],0) = nyL;
			K(mask[2],1) = nxL;
			K(mask[2],2) = nyL;
			K(mask[2],3) = nxL;
			K(0,mask[2]) = nyL;
			K(1,mask[2]) = nxL;
			K(2,mask[2]) = nyL;
			K(3,mask[2]) = nxL;
		}
	}

	// form residual

	noalias(rRightHandSideVector) -= prod( rLeftHandSideMatrix, currentValues );
}

void RveWeakPeriodicCondition2D2N::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    MatrixType dummy;
	CalculateLocalSystem(dummy, rRightHandSideVector, rCurrentProcessInfo);
}

void RveWeakPeriodicCondition2D2N::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
{
	//if(rResult.size() != 7) rResult.resize(7);
	//GeometryType& geom = GetGeometry();
	//// first node
	//rResult[0] = geom[0].GetDof( DISPLACEMENT_X ).EquationId();
	//rResult[1] = geom[0].GetDof( DISPLACEMENT_Y ).EquationId();
	//// second node
	//rResult[2] = geom[1].GetDof( DISPLACEMENT_X ).EquationId();
	//rResult[3] = geom[1].GetDof( DISPLACEMENT_Y ).EquationId();
	//// lagrangian node
	//rResult[4] = geom[2].GetDof( RVE_WPC_LAGRANGIAN_DOF_X ).EquationId();
	//rResult[5] = geom[2].GetDof( RVE_WPC_LAGRANGIAN_DOF_Y ).EquationId();
	//rResult[6] = geom[2].GetDof( RVE_WPC_LAGRANGIAN_DOF_Z ).EquationId();

	if(m_is_skew_symmetric_constraint)
	{
		unsigned int n_dofs = 5;

		if(rResult.size() != n_dofs) rResult.resize(n_dofs);
		GeometryType& geom = GetGeometry();
		// first node
		rResult[0] = geom[0].GetDof( DISPLACEMENT_X ).EquationId();
		rResult[1] = geom[0].GetDof( DISPLACEMENT_Y ).EquationId();
		// second node
		rResult[2] = geom[1].GetDof( DISPLACEMENT_X ).EquationId();
		rResult[3] = geom[1].GetDof( DISPLACEMENT_Y ).EquationId();
		// lagrangian node
		rResult[4] = geom[2].GetDof( RVE_WPR_LAGRANGIAN_DOF ).EquationId();
	}
	else
	{
		array_1d<unsigned int, 3> mask;
		unsigned int n_lag_dofs = CalculateLagrangianDofMask(mask);
		unsigned int n_dofs = 4 + n_lag_dofs;

		if(rResult.size() != n_dofs) rResult.resize(n_dofs);
		GeometryType& geom = GetGeometry();
		// first node
		rResult[0] = geom[0].GetDof( DISPLACEMENT_X ).EquationId();
		rResult[1] = geom[0].GetDof( DISPLACEMENT_Y ).EquationId();
		// second node
		rResult[2] = geom[1].GetDof( DISPLACEMENT_X ).EquationId();
		rResult[3] = geom[1].GetDof( DISPLACEMENT_Y ).EquationId();
		// lagrangian node
		if(mask[0] > 0) rResult[mask[0]] = geom[2].GetDof( RVE_WPC_LAGRANGIAN_DOF_X ).EquationId();
		if(mask[1] > 0) rResult[mask[1]] = geom[2].GetDof( RVE_WPC_LAGRANGIAN_DOF_Y ).EquationId();
		if(mask[2] > 0) rResult[mask[2]] = geom[2].GetDof( RVE_WPC_LAGRANGIAN_DOF_Z ).EquationId();
	}
}

void RveWeakPeriodicCondition2D2N::GetDofList(DofsVectorType& ConditionalDofList,ProcessInfo& CurrentProcessInfo)
{
	//if(ConditionalDofList.size() != 7) ConditionalDofList.resize(7);
	//GeometryType& geom = GetGeometry();
	//// first node
	//ConditionalDofList[0] = geom[0].pGetDof( DISPLACEMENT_X );
	//ConditionalDofList[1] = geom[0].pGetDof( DISPLACEMENT_Y );
	//// second node
	//ConditionalDofList[2] = geom[1].pGetDof( DISPLACEMENT_X );
	//ConditionalDofList[3] = geom[1].pGetDof( DISPLACEMENT_Y );
	//// lagrangian node
	//ConditionalDofList[4] = geom[2].pGetDof( RVE_WPC_LAGRANGIAN_DOF_X );
	//ConditionalDofList[5] = geom[2].pGetDof( RVE_WPC_LAGRANGIAN_DOF_Y );
	//ConditionalDofList[6] = geom[2].pGetDof( RVE_WPC_LAGRANGIAN_DOF_Z );

	if(m_is_skew_symmetric_constraint)
	{
		unsigned int n_dofs = 5;

		if(ConditionalDofList.size() != n_dofs) ConditionalDofList.resize(n_dofs);
		GeometryType& geom = GetGeometry();
		// first node
		ConditionalDofList[0] = geom[0].pGetDof( DISPLACEMENT_X );
		ConditionalDofList[1] = geom[0].pGetDof( DISPLACEMENT_Y );
		// second node
		ConditionalDofList[2] = geom[1].pGetDof( DISPLACEMENT_X );
		ConditionalDofList[3] = geom[1].pGetDof( DISPLACEMENT_Y );
		// lagrangian node
		ConditionalDofList[4] = geom[2].pGetDof( RVE_WPR_LAGRANGIAN_DOF );
	}
	else
	{
		array_1d<unsigned int, 3> mask;
		unsigned int n_lag_dofs = CalculateLagrangianDofMask(mask);
		unsigned int n_dofs = 4 + n_lag_dofs;

		if(ConditionalDofList.size() != n_dofs) ConditionalDofList.resize(n_dofs);
		GeometryType& geom = GetGeometry();
		// first node
		ConditionalDofList[0] = geom[0].pGetDof( DISPLACEMENT_X );
		ConditionalDofList[1] = geom[0].pGetDof( DISPLACEMENT_Y );
		// second node
		ConditionalDofList[2] = geom[1].pGetDof( DISPLACEMENT_X );
		ConditionalDofList[3] = geom[1].pGetDof( DISPLACEMENT_Y );
		// lagrangian node
		if(mask[0] > 0) ConditionalDofList[mask[0]] = geom[2].pGetDof( RVE_WPC_LAGRANGIAN_DOF_X );
		if(mask[1] > 0) ConditionalDofList[mask[1]] = geom[2].pGetDof( RVE_WPC_LAGRANGIAN_DOF_Y );
		if(mask[2] > 0) ConditionalDofList[mask[2]] = geom[2].pGetDof( RVE_WPC_LAGRANGIAN_DOF_Z );
	}
}

int RveWeakPeriodicCondition2D2N::Check(const ProcessInfo& rCurrentProcessInfo)
{
	MyBase::Check(rCurrentProcessInfo);
	KRATOS_TRY

	GeometryType & geom = this->GetGeometry();

	if(geom.size() != 3) {
		std::stringstream ss;
		ss << "RveWeakPeriodicCondition2D2N - The Geometry should have 2 nodes + 1 lagrangian node. Condition with ID = " << this->GetId();
		KRATOS_THROW_ERROR(std::logic_error, ss.str(), "");
	}

	//verify that the variables are correctly initialized
	if(DISPLACEMENT.Key() == 0)
		KRATOS_THROW_ERROR(std::invalid_argument,"DISPLACEMENT has Key zero! (check if the application is correctly registered","");

	//verify that the dofs exist
	for(unsigned int i = 0; i < 2; i++)
	{
		NodeType & iNode = geom[i];

		if(iNode.SolutionStepsDataHas(DISPLACEMENT) == false)
			KRATOS_THROW_ERROR(std::invalid_argument,"missing variable DISPLACEMENT on node ", iNode.Id());

		if(iNode.HasDofFor(DISPLACEMENT_X) == false || iNode.HasDofFor(DISPLACEMENT_Y) == false)
			KRATOS_THROW_ERROR(std::invalid_argument,"missing one of the dofs for the variable DISPLACEMENT on node ", iNode.Id());
	}
	NodeType & lagNode = geom[2];
	if(lagNode.SolutionStepsDataHas(RVE_WPC_LAGRANGIAN_DOF) == false)
		KRATOS_THROW_ERROR(std::invalid_argument,"missing variable RVE_WPC_LAGRANGIAN_DOF on node ", lagNode.Id());
	if( lagNode.HasDofFor(RVE_WPC_LAGRANGIAN_DOF_X) == false ||
		lagNode.HasDofFor(RVE_WPC_LAGRANGIAN_DOF_Y) == false ||
		lagNode.HasDofFor(RVE_WPC_LAGRANGIAN_DOF_Z) == false) {
		KRATOS_THROW_ERROR(std::invalid_argument,"missing one of the dofs for the variable RVE_WPC_LAGRANGIAN_DOF on node ", lagNode.Id());
	}
	if(m_is_skew_symmetric_constraint)
	{
		if(lagNode.SolutionStepsDataHas(RVE_WPR_LAGRANGIAN_DOF) == false)
			KRATOS_THROW_ERROR(std::invalid_argument,"missing variable RVE_WPR_LAGRANGIAN_DOF on node ", lagNode.Id());
	}

	KRATOS_CATCH("")
	return 0;
}

unsigned int RveWeakPeriodicCondition2D2N::CalculateLagrangianDofMask(array_1d<unsigned int, 3>& mask)
{
	GeometryType& geom = GetGeometry();
	double x0 = geom[0].X0();
	double y0 = geom[0].Y0();
	double x1 = geom[1].X0();
	double y1 = geom[1].Y0();

	double tx = (x1-x0);
	double ty = (y1-y0);
	double L  = std::sqrt(tx*tx + ty*ty);

	double nx = (y1-y0)/L;
	double ny = (x0-x1)/L;

	unsigned int index = 3;
	double tol = 1.0e-8;
	if(std::abs(nx) > tol)
		mask[0] = ++index;
	else
		mask[0] = 0;
	if(std::abs(ny) > tol)
		mask[1] = ++index;
	else
		mask[1] = 0;
	mask[2] = ++index;
	return index-3;
	/*mask[0] = 4;
	mask[1] = 5;
	mask[2] = 6;
	return 3;*/
}

}