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
#include "shell_edge_load_condition_3D2N.h"
#include "utilities/math_utils.h"
#include "geometries/point_3d.h"
#include "multiscale_application.h"
#include "custom_utilities/math_helpers.h"

namespace Kratos
{

ShellEdgeLoadCondition3D2N::ShellEdgeLoadCondition3D2N(IndexType NewId, GeometryType::Pointer pGeometry)
	: MyBase(NewId, pGeometry)
{
}

ShellEdgeLoadCondition3D2N::ShellEdgeLoadCondition3D2N(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
	: MyBase(NewId, pGeometry, pProperties)
{
}

ShellEdgeLoadCondition3D2N::ShellEdgeLoadCondition3D2N(const ShellEdgeLoadCondition3D2N& rOther)
	: MyBase(rOther)
{
}

Condition::Pointer ShellEdgeLoadCondition3D2N::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
{
	return Condition::Pointer(new ShellEdgeLoadCondition3D2N(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

ShellEdgeLoadCondition3D2N::~ShellEdgeLoadCondition3D2N()
{
}

void ShellEdgeLoadCondition3D2N::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
	GeometryType& geom = GetGeometry();
	unsigned int num_nodes = geom.PointsNumber();
	unsigned int num_dofs = 6*num_nodes;

	if(rLeftHandSideMatrix.size1() != num_dofs || rLeftHandSideMatrix.size2() != num_dofs) rLeftHandSideMatrix.resize(num_dofs, num_dofs,false);
	noalias(rLeftHandSideMatrix) = ZeroMatrix(num_dofs, num_dofs);

	if(rRightHandSideVector.size() != num_dofs) rRightHandSideVector.resize(num_dofs, false);
	noalias( rRightHandSideVector ) = ZeroVector(num_dofs);

	PropertiesType& props = GetProperties();
	double th = props[THICKNESS];

	// todo: generalize it for any geometry / load pattern
	double lx = geom[1].X0() - geom[0].X0();
	double ly = geom[1].X0() - geom[0].X0();
	double lz = geom[1].X0() - geom[0].X0();
	double L = std::sqrt(lx*lx + ly*ly + lz*lz);

	array_1d<double, 3> F_vector_over_2;
	array_1d<double, 3> M_vector_over_2;
	F_vector_over_2.clear();
	M_vector_over_2.clear();
	for(unsigned int i = 0; i < num_nodes; i++)
	{
		if(geom[i].SolutionStepsDataHas(SHELL_EDGE_FORCE))
		{
			const array_1d<double,3>& F = geom[i].FastGetSolutionStepValue(SHELL_EDGE_FORCE);

		}
		if(geom[i].SolutionStepsDataHas(SHELL_EDGE_TORQUE))
		{
			const array_1d<double,3>& M = geom[i].FastGetSolutionStepValue(SHELL_EDGE_TORQUE);
		}
	}
}

void ShellEdgeLoadCondition3D2N::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    MatrixType dummy;
	CalculateLocalSystem(dummy, rRightHandSideVector, rCurrentProcessInfo);
}

void ShellEdgeLoadCondition3D2N::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
{
	GeometryType& geom = GetGeometry();
	unsigned int num_nodes = geom.PointsNumber();
	unsigned int num_dofs = 6*num_nodes;
	// resize
	if(rResult.size() != num_dofs) rResult.resize(num_dofs);
	for(unsigned int i = 0; i < num_nodes; i++)
	{
		unsigned int index = i*6;
		rResult[index  ] = geom[i].GetDof( DISPLACEMENT_X ).EquationId();
		rResult[index+1] = geom[i].GetDof( DISPLACEMENT_Y ).EquationId();
		rResult[index+2] = geom[i].GetDof( DISPLACEMENT_Z ).EquationId();
		rResult[index+3] = geom[i].GetDof( ROTATION_X ).EquationId();
		rResult[index+4] = geom[i].GetDof( ROTATION_Y ).EquationId();
		rResult[index+5] = geom[i].GetDof( ROTATION_Z ).EquationId();
	}
}

void ShellEdgeLoadCondition3D2N::GetDofList(DofsVectorType& ConditionalDofList,ProcessInfo& CurrentProcessInfo)
{
	GeometryType& geom = GetGeometry();
	unsigned int num_nodes = geom.PointsNumber();
	unsigned int num_dofs = 6*num_nodes;
	// resize
	if(ConditionalDofList.size() != num_dofs) ConditionalDofList.resize(num_dofs);
	for(unsigned int i = 0; i < num_nodes; i++)
	{
		unsigned int index = i*6;
		ConditionalDofList[index  ] = geom[i].pGetDof( DISPLACEMENT_X );
		ConditionalDofList[index+1] = geom[i].pGetDof( DISPLACEMENT_Y );
		ConditionalDofList[index+2] = geom[i].pGetDof( DISPLACEMENT_Z );
		ConditionalDofList[index+3] = geom[i].pGetDof( ROTATION_X );
		ConditionalDofList[index+4] = geom[i].pGetDof( ROTATION_Y );
		ConditionalDofList[index+5] = geom[i].pGetDof( ROTATION_Z );
	}
}

int ShellEdgeLoadCondition3D2N::Check(const ProcessInfo& rCurrentProcessInfo)
{
	MyBase::Check(rCurrentProcessInfo);
	KRATOS_TRY

	GeometryType & geom = this->GetGeometry();
	unsigned int num_nodes = geom.PointsNumber();
	if(num_nodes != 2)
		KRATOS_THROW_ERROR(std::invalid_argument,"ShellEdgeLoadCondition works only with linear elements","");

	//verify that the variables are correctly initialized
	if(DISPLACEMENT.Key() == 0)
		KRATOS_THROW_ERROR(std::invalid_argument,"DISPLACEMENT has Key zero! (check if the application is correctly registered","");
	if(ROTATION.Key() == 0)
		KRATOS_THROW_ERROR(std::invalid_argument,"ROTATION has Key zero! (check if the application is correctly registered","");

	//verify that the dofs exist
	for(unsigned int i = 0; i < num_nodes; i++)
	{
		NodeType & iNode = geom[i];

		if(iNode.SolutionStepsDataHas(DISPLACEMENT) == false)
			KRATOS_THROW_ERROR(std::invalid_argument,"missing variable DISPLACEMENT on node ", iNode.Id());

		if(iNode.HasDofFor(DISPLACEMENT_X) == false || iNode.HasDofFor(DISPLACEMENT_Y) == false || iNode.HasDofFor(DISPLACEMENT_Z) == false)
			KRATOS_THROW_ERROR(std::invalid_argument,"missing one of the dofs for the variable DISPLACEMENT on node ", iNode.Id());

		if(iNode.SolutionStepsDataHas(ROTATION) == false)
			KRATOS_THROW_ERROR(std::invalid_argument,"missing variable ROTATION on node ", iNode.Id());

		if(iNode.HasDofFor(ROTATION_X) == false || iNode.HasDofFor(ROTATION_Y) == false || iNode.HasDofFor(ROTATION_Z) == false)
			KRATOS_THROW_ERROR(std::invalid_argument,"missing one of the dofs for the variable ROTATION on node ", iNode.Id());
	}

	PropertiesType & props = this->GetProperties();
	if(!props.Has(THICKNESS))
		KRATOS_THROW_ERROR(std::invalid_argument,"ShellEdgeLoadCondition needs THICKNESS property","");

	KRATOS_CATCH("")
	return 0;
}

}