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
#include "rve_weak_rotation_constraint_condition_3D.h"
#include "utilities/math_utils.h"
#include "geometries/point_3d.h"
#include "multiscale_application_variables.h"
#include "custom_utilities/math_helpers.h"

namespace Kratos
{

RveWeakRotationCondition3D::RveWeakRotationCondition3D(IndexType NewId, GeometryType::Pointer pGeometry)
	: MyBase(NewId, pGeometry)
{
}

RveWeakRotationCondition3D::RveWeakRotationCondition3D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
	: MyBase(NewId, pGeometry, pProperties)
{
}

RveWeakRotationCondition3D::RveWeakRotationCondition3D(const RveWeakRotationCondition3D& rOther)
	: MyBase(rOther)
{
}

Condition::Pointer RveWeakRotationCondition3D::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
{
	return Condition::Pointer(new RveWeakRotationCondition3D(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

RveWeakRotationCondition3D::~RveWeakRotationCondition3D()
{
}

void RveWeakRotationCondition3D::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
	GeometryType& geom = GetGeometry();
	unsigned int num_nodes = geom.PointsNumber();
	unsigned int num_dofs_u = num_nodes*2; // note (we assume that the RVE is flat on the XY plane, so we don't care about the Z rotation)
	unsigned int num_dofs_lag = 2;
	unsigned int num_dofs = num_dofs_u + num_dofs_lag;

	if(rLeftHandSideMatrix.size1() != num_dofs || rLeftHandSideMatrix.size2() != num_dofs) rLeftHandSideMatrix.resize(num_dofs, num_dofs,false);
	noalias(rLeftHandSideMatrix) = ZeroMatrix(num_dofs, num_dofs);

	if(rRightHandSideVector.size() != num_dofs) rRightHandSideVector.resize(num_dofs, false);
	noalias( rRightHandSideVector ) = ZeroVector(num_dofs);

	double V = geom.DomainSize();
	double l1 = mLagrangianNode->FastGetSolutionStepValue(RVE_SHELL_WRC_LAGRANGIAN_DOF_X);
	double l2 = mLagrangianNode->FastGetSolutionStepValue(RVE_SHELL_WRC_LAGRANGIAN_DOF_Y);
	for(unsigned int inode = 0; inode < num_nodes; inode++)
	{
		unsigned int index = inode*2;
		// [C]
		rLeftHandSideMatrix(num_dofs-2, index  ) = V;
		rLeftHandSideMatrix(num_dofs-1, index+1) = V;
		// [C]^T
		rLeftHandSideMatrix(index  , num_dofs-2) = V;
		rLeftHandSideMatrix(index+1, num_dofs-1) = V;
		// R
		rRightHandSideVector(index)   = -V*l1;
		rRightHandSideVector(index+1) = -V*l2;
	}
}

void RveWeakRotationCondition3D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    MatrixType dummy;
	CalculateLocalSystem(dummy, rRightHandSideVector, rCurrentProcessInfo);
}

void RveWeakRotationCondition3D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
{
	GeometryType& geom = GetGeometry();
	unsigned int num_nodes = geom.PointsNumber();
	unsigned int num_dofs_u = num_nodes*2; // note (we assume that the RVE is flat on the XY plane, so we don't care about the Z rotation)
	unsigned int num_dofs_lag = 2;
	unsigned int num_dofs = num_dofs_u + num_dofs_lag;
	// resize
	if(rResult.size() != num_dofs) rResult.resize(num_dofs);
	for(unsigned int i = 0; i < num_nodes; i++)
	{
		unsigned int index = i*2;
		rResult[index  ] = geom[i].GetDof( ROTATION_X ).EquationId();
		rResult[index+1] = geom[i].GetDof( ROTATION_Y ).EquationId();
	}
	// lagrangian node
	rResult[num_dofs-2] = mLagrangianNode->GetDof( RVE_SHELL_WRC_LAGRANGIAN_DOF_X ).EquationId();
	rResult[num_dofs-1] = mLagrangianNode->GetDof( RVE_SHELL_WRC_LAGRANGIAN_DOF_Y ).EquationId();
}

void RveWeakRotationCondition3D::GetDofList(DofsVectorType& ConditionalDofList,ProcessInfo& CurrentProcessInfo)
{
	GeometryType& geom = GetGeometry();
	unsigned int num_nodes = geom.PointsNumber();
	unsigned int num_dofs_u = num_nodes*2; // note (we assume that the RVE is flat on the XY plane, so we don't care about the Z rotation)
	unsigned int num_dofs_lag = 2;
	unsigned int num_dofs = num_dofs_u + num_dofs_lag;
	// resize
	if(ConditionalDofList.size() != num_dofs) ConditionalDofList.resize(num_dofs);
	for(unsigned int i = 0; i < num_nodes; i++)
	{
		unsigned int index = i*2;
		ConditionalDofList[index  ] = geom[i].pGetDof( ROTATION_X );
		ConditionalDofList[index+1] = geom[i].pGetDof( ROTATION_Y );
	}
	// lagrangian node
	ConditionalDofList[num_dofs-2] = mLagrangianNode->pGetDof( RVE_SHELL_WRC_LAGRANGIAN_DOF_X );
	ConditionalDofList[num_dofs-1] = mLagrangianNode->pGetDof( RVE_SHELL_WRC_LAGRANGIAN_DOF_Y );
}

int RveWeakRotationCondition3D::Check(const ProcessInfo& rCurrentProcessInfo)
{
	MyBase::Check(rCurrentProcessInfo);
	KRATOS_TRY

	GeometryType & geom = this->GetGeometry();

	//verify that the variables are correctly initialized
	if(ROTATION.Key() == 0)
		KRATOS_THROW_ERROR(std::invalid_argument,"ROTATION has Key zero! (check if the application is correctly registered","");

	//verify that the dofs exist
	for(unsigned int i = 0; i < 2; i++)
	{
		NodeType & iNode = geom[i];

		if(iNode.SolutionStepsDataHas(ROTATION) == false)
			KRATOS_THROW_ERROR(std::invalid_argument,"missing variable ROTATION on node ", iNode.Id());

		if(iNode.HasDofFor(ROTATION_X) == false || iNode.HasDofFor(ROTATION_Y) == false)
			KRATOS_THROW_ERROR(std::invalid_argument,"missing one of the dofs for the variable ROTATION on node ", iNode.Id());
	}
	if(!mLagrangianNode)
		KRATOS_THROW_ERROR(std::invalid_argument,"missing Lagrangian node", "");
	NodeType & lagNode = *mLagrangianNode;
	if(lagNode.SolutionStepsDataHas(RVE_SHELL_WRC_LAGRANGIAN_DOF_X) == false)
		KRATOS_THROW_ERROR(std::invalid_argument,"missing variable RVE_SHELL_WRC_LAGRANGIAN_DOF_X on node ", lagNode.Id());
	if(lagNode.SolutionStepsDataHas(RVE_SHELL_WRC_LAGRANGIAN_DOF_Y) == false)
		KRATOS_THROW_ERROR(std::invalid_argument,"missing variable RVE_SHELL_WRC_LAGRANGIAN_DOF_Y on node ", lagNode.Id());

	KRATOS_CATCH("")
	return 0;
}

}