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
#include "rve_weak_tshear_constraint_condition_3D.h"
#include "utilities/math_utils.h"
#include "geometries/point_3d.h"
#include "multiscale_application_variables.h"
#include "custom_utilities/math_helpers.h"

namespace Kratos
{

RveWeakTShearCondition3D::RveWeakTShearCondition3D(IndexType NewId, GeometryType::Pointer pGeometry)
	: MyBase(NewId, pGeometry)
{
}

RveWeakTShearCondition3D::RveWeakTShearCondition3D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
	: MyBase(NewId, pGeometry, pProperties)
{
}

RveWeakTShearCondition3D::RveWeakTShearCondition3D(const RveWeakTShearCondition3D& rOther)
	: MyBase(rOther)
{
}

Condition::Pointer RveWeakTShearCondition3D::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
{
	return Condition::Pointer(new RveWeakTShearCondition3D(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

RveWeakTShearCondition3D::~RveWeakTShearCondition3D()
{
}

void RveWeakTShearCondition3D::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
	GeometryType& geom = GetGeometry();
	unsigned int num_nodes = geom.PointsNumber();
	unsigned int num_dofs_u = num_nodes*3; // note (we assume that the RVE is flat on the XY plane, so we don't care about the Z rotation)
	unsigned int num_dofs_lag = 2;
	unsigned int num_dofs = num_dofs_u + num_dofs_lag;

	if(rLeftHandSideMatrix.size1() != num_dofs || rLeftHandSideMatrix.size2() != num_dofs) rLeftHandSideMatrix.resize(num_dofs, num_dofs,false);
	noalias(rLeftHandSideMatrix) = ZeroMatrix(num_dofs, num_dofs);

	if(rRightHandSideVector.size() != num_dofs) rRightHandSideVector.resize(num_dofs, false);
	noalias( rRightHandSideVector ) = ZeroVector(num_dofs);

	double l1 = mLagrangianNode->FastGetSolutionStepValue(RVE_SHELL_WRC_LAGRANGIAN_DOF_X);
	double l2 = mLagrangianNode->FastGetSolutionStepValue(RVE_SHELL_WRC_LAGRANGIAN_DOF_Y);

	double x1 = geom[0].X0(); double y1 = geom[0].Y0();
	double x2 = geom[1].X0(); double y2 = geom[1].Y0();
	double x3 = geom[2].X0(); double y3 = geom[2].Y0();
	double x4 = geom[3].X0(); double y4 = geom[3].Y0();

	double xi = 1.0/std::sqrt(3.0);
	double eta = xi;
	double C1 = xi*xi;
	double C2 = eta*eta;
	double C3 = (C2*x3*y4)/8;
	double C4 = (C2*x4*y2)/8;
	double C5 = (C1*x4*y2)/8;
	double C6 = (C2*x1*y3)/8;
	double C7 = (C1*x2*y3)/8;

	Matrix& K = rLeftHandSideMatrix;
	K(12,0) = x4/2 - x2/2;
	K(0,12) = K(12,0);
	K(12,1) = C3 + C4 + C5 + C6 + C7 - (x1*y2)/8 + (x2*y1)/8 + (x1*y4)/8 - (x2*y3)/8 + (x3*y2)/8 - (x4*y1)/8 - (x3*y4)/8 + (x4*y3)/8 - (C1*x1*y3)/8 + (C1*x3*y1)/8 - (C2*x1*y2)/8 + (C2*x2*y1)/8 + (C1*x1*y4)/8 - (C1*x3*y2)/8 - (C1*x4*y1)/8 - (C2*x3*y1)/8 - (C1*x2*y4)/8 - (C2*x2*y4)/8 - (C2*x4*y3)/8;
	K(1,12) = K(12,1);
	K(12,3) = x1/2 - x3/2;
	K(3,12) = K(12,3);
	K(12,4) = C3 + C4 - C5 + C6 - C7 - (x1*y2)/8 + (x2*y1)/8 + (x1*y4)/8 - (x2*y3)/8 + (x3*y2)/8 - (x4*y1)/8 - (x3*y4)/8 + (x4*y3)/8 + (C1*x1*y3)/8 - (C1*x3*y1)/8 - (C2*x1*y2)/8 + (C2*x2*y1)/8 - (C1*x1*y4)/8 + (C1*x3*y2)/8 + (C1*x4*y1)/8 - (C2*x3*y1)/8 + (C1*x2*y4)/8 - (C2*x2*y4)/8 - (C2*x4*y3)/8;
	K(4,12) = K(12,4);
	K(12,6) = x2/2 - x4/2;
	K(6,12) = K(12,6);
	K(12,7) = (x2*y1)/8 - C4 - C5 - C6 - C7 - (x1*y2)/8 - C3 + (x1*y4)/8 - (x2*y3)/8 + (x3*y2)/8 - (x4*y1)/8 - (x3*y4)/8 + (x4*y3)/8 + (C1*x1*y3)/8 - (C1*x3*y1)/8 + (C2*x1*y2)/8 - (C2*x2*y1)/8 - (C1*x1*y4)/8 + (C1*x3*y2)/8 + (C1*x4*y1)/8 + (C2*x3*y1)/8 + (C1*x2*y4)/8 + (C2*x2*y4)/8 + (C2*x4*y3)/8;
	K(7,12) = K(12,7);
	K(12,9) = x3/2 - x1/2;
	K(9,12) = K(12,9);
	K(12,10) = C5 - C4 - C3 - C6 + C7 - (x1*y2)/8 + (x2*y1)/8 + (x1*y4)/8 - (x2*y3)/8 + (x3*y2)/8 - (x4*y1)/8 - (x3*y4)/8 + (x4*y3)/8 - (C1*x1*y3)/8 + (C1*x3*y1)/8 + (C2*x1*y2)/8 - (C2*x2*y1)/8 + (C1*x1*y4)/8 - (C1*x3*y2)/8 - (C1*x4*y1)/8 + (C2*x3*y1)/8 - (C1*x2*y4)/8 + (C2*x2*y4)/8 + (C2*x4*y3)/8;
	K(10,12) = K(12,10);
	K(13,0) = y2/2 - y4/2;
	K(0,13) = K(13,0);
	K(13,2) = (x1*y2)/8 - C4 - C5 - C6 - C7 - C3 - (x2*y1)/8 - (x1*y4)/8 + (x2*y3)/8 - (x3*y2)/8 + (x4*y1)/8 + (x3*y4)/8 - (x4*y3)/8 + (C1*x1*y3)/8 - (C1*x3*y1)/8 + (C2*x1*y2)/8 - (C2*x2*y1)/8 - (C1*x1*y4)/8 + (C1*x3*y2)/8 + (C1*x4*y1)/8 + (C2*x3*y1)/8 + (C1*x2*y4)/8 + (C2*x2*y4)/8 + (C2*x4*y3)/8;
	K(2,13) = K(13,2);
	K(13,3) = y3/2 - y1/2;
	K(3,13) = K(13,3);
	K(13,5) = C5 - C4 - C3 - C6 + C7 + (x1*y2)/8 - (x2*y1)/8 - (x1*y4)/8 + (x2*y3)/8 - (x3*y2)/8 + (x4*y1)/8 + (x3*y4)/8 - (x4*y3)/8 - (C1*x1*y3)/8 + (C1*x3*y1)/8 + (C2*x1*y2)/8 - (C2*x2*y1)/8 + (C1*x1*y4)/8 - (C1*x3*y2)/8 - (C1*x4*y1)/8 + (C2*x3*y1)/8 - (C1*x2*y4)/8 + (C2*x2*y4)/8 + (C2*x4*y3)/8;
	K(5,13) = K(13,5);
	K(13,6) = y4/2 - y2/2;
	K(6,13) = K(13,6);
	K(13,8) = C3 + C4 + C5 + C6 + C7 + (x1*y2)/8 - (x2*y1)/8 - (x1*y4)/8 + (x2*y3)/8 - (x3*y2)/8 + (x4*y1)/8 + (x3*y4)/8 - (x4*y3)/8 - (C1*x1*y3)/8 + (C1*x3*y1)/8 - (C2*x1*y2)/8 + (C2*x2*y1)/8 + (C1*x1*y4)/8 - (C1*x3*y2)/8 - (C1*x4*y1)/8 - (C2*x3*y1)/8 - (C1*x2*y4)/8 - (C2*x2*y4)/8 - (C2*x4*y3)/8;
	K(8,13) = K(13,8);
	K(13,9) = y1/2 - y3/2;
	K(9,13) = K(13,9);
	K(13,11) = C3 + C4 - C5 + C6 - C7 + (x1*y2)/8 - (x2*y1)/8 - (x1*y4)/8 + (x2*y3)/8 - (x3*y2)/8 + (x4*y1)/8 + (x3*y4)/8 - (x4*y3)/8 + (C1*x1*y3)/8 - (C1*x3*y1)/8 - (C2*x1*y2)/8 + (C2*x2*y1)/8 - (C1*x1*y4)/8 + (C1*x3*y2)/8 + (C1*x4*y1)/8 - (C2*x3*y1)/8 + (C1*x2*y4)/8 - (C2*x2*y4)/8 - (C2*x4*y3)/8;
	K(11,13) = K(13,11);

	Vector& R = rRightHandSideVector;
	R(0) = - l1*K(0, 12) - l2*K(0, 13);
	R(1) = -l1*K(1, 12);
	R(2) = -l2*K(2, 13);
	R(3) = - l1*K(3, 12) - l2*K(3, 13);
	R(4) = -l1*K(4, 12);
	R(5) = -l2*K(5, 13);
	R(6) = - l1*K(6, 12) - l2*K(6, 13);
	R(7) = -l1*K(7, 12);
	R(8) = -l2*K(8, 13);
	R(9) = - l1*K(9, 12) - l2*K(9, 13);
	R(10) = -l1*K(10, 12);
	R(11) = -l2*K(11, 13);
}

void RveWeakTShearCondition3D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    MatrixType dummy;
	CalculateLocalSystem(dummy, rRightHandSideVector, rCurrentProcessInfo);
}

void RveWeakTShearCondition3D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
{
	GeometryType& geom = GetGeometry();
	unsigned int num_nodes = geom.PointsNumber();
	unsigned int num_dofs_u = num_nodes*3; // note (we assume that the RVE is flat on the XY plane, so we don't care about the Z rotation)
	unsigned int num_dofs_lag = 2;
	unsigned int num_dofs = num_dofs_u + num_dofs_lag;
	// resize
	if(rResult.size() != num_dofs) rResult.resize(num_dofs);
	for(unsigned int i = 0; i < num_nodes; i++)
	{
		unsigned int index = i*3;
		rResult[index  ] = geom[i].GetDof( DISPLACEMENT_Z ).EquationId();
		rResult[index+1] = geom[i].GetDof( ROTATION_X ).EquationId();
		rResult[index+2] = geom[i].GetDof( ROTATION_Y ).EquationId();
	}
	// lagrangian node
	rResult[num_dofs-2] = mLagrangianNode->GetDof( RVE_SHELL_WRC_LAGRANGIAN_DOF_X ).EquationId();
	rResult[num_dofs-1] = mLagrangianNode->GetDof( RVE_SHELL_WRC_LAGRANGIAN_DOF_Y ).EquationId();
}

void RveWeakTShearCondition3D::GetDofList(DofsVectorType& ConditionalDofList,ProcessInfo& CurrentProcessInfo)
{
	GeometryType& geom = GetGeometry();
	unsigned int num_nodes = geom.PointsNumber();
	unsigned int num_dofs_u = num_nodes*3; // note (we assume that the RVE is flat on the XY plane, so we don't care about the Z rotation)
	unsigned int num_dofs_lag = 2;
	unsigned int num_dofs = num_dofs_u + num_dofs_lag;
	// resize
	if(ConditionalDofList.size() != num_dofs) ConditionalDofList.resize(num_dofs);
	for(unsigned int i = 0; i < num_nodes; i++)
	{
		unsigned int index = i*3;
		ConditionalDofList[index  ] = geom[i].pGetDof( DISPLACEMENT_Z );
		ConditionalDofList[index+1] = geom[i].pGetDof( ROTATION_X );
		ConditionalDofList[index+2] = geom[i].pGetDof( ROTATION_Y );
	}
	// lagrangian node
	ConditionalDofList[num_dofs-2] = mLagrangianNode->pGetDof( RVE_SHELL_WRC_LAGRANGIAN_DOF_X );
	ConditionalDofList[num_dofs-1] = mLagrangianNode->pGetDof( RVE_SHELL_WRC_LAGRANGIAN_DOF_Y );
}

int RveWeakTShearCondition3D::Check(const ProcessInfo& rCurrentProcessInfo)
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

		if(iNode.SolutionStepsDataHas(DISPLACEMENT) == false)
			KRATOS_THROW_ERROR(std::invalid_argument,"missing variable DISPLACEMENT on node ", iNode.Id());
		if(iNode.SolutionStepsDataHas(ROTATION) == false)
			KRATOS_THROW_ERROR(std::invalid_argument,"missing variable ROTATION on node ", iNode.Id());

		if(iNode.HasDofFor(DISPLACEMENT_Z) == false)
			KRATOS_THROW_ERROR(std::invalid_argument,"missing one of the dofs for the variable DISPLACEMENT on node ", iNode.Id());
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