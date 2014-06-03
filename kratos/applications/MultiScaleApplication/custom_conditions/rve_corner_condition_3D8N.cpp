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
#include "rve_corner_condition_3D8N.h"
#include "utilities/math_utils.h"
#include "geometries/point_3d.h"
#include "multiscale_application.h"

namespace Kratos
{
	
RveCornerCondition3D8N::RveCornerCondition3D8N(IndexType NewId, GeometryType::Pointer pGeometry)
	: MyBase(NewId, pGeometry)
{
}

RveCornerCondition3D8N::RveCornerCondition3D8N(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
	: MyBase(NewId, pGeometry, pProperties)
{
}
	
RveCornerCondition3D8N::RveCornerCondition3D8N(const RveCornerCondition3D8N& rOther)
	: MyBase(rOther)
{
}

Condition::Pointer RveCornerCondition3D8N::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
{
	return Condition::Pointer(new RveCornerCondition3D8N(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

RveCornerCondition3D8N::~RveCornerCondition3D8N()
{
}

void RveCornerCondition3D8N::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    if(rLeftHandSideMatrix.size1() != 48) rLeftHandSideMatrix.resize(48, 48,false);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(48, 48);

    if(rRightHandSideVector.size() != 48) rRightHandSideVector.resize(48, false);
	noalias( rRightHandSideVector ) = ZeroVector(48);

	GeometryType& geom = GetGeometry();

	// setup RHS vector with prescribed values

	Vector strainVector(6);
	if( this->GetMacroStrainVector(strainVector) )
	{
		double exx = strainVector[0];
		double eyy = strainVector[1];
		double ezz = strainVector[2];
		double exy = strainVector[3] * 0.5;
		double eyz = strainVector[4] * 0.5;
		double exz = strainVector[5] * 0.5;

		double lx = geom[1].X0() - geom[0].X0();
		double ly = geom[3].Y0() - geom[0].Y0();
		double lz = geom[4].Z0() - geom[0].Z0();

		// [0, 0, 0] : bottom-left corner -> always fixed

		// [1, 0, 0] : corner
		rRightHandSideVector( 9) = exx * lx;
		rRightHandSideVector(10) = exy * lx;
		rRightHandSideVector(11) = exz * lx;

		// [1, 1, 0] : corner
		rRightHandSideVector(15) = exx * lx + exy * ly;
		rRightHandSideVector(16) = exy * lx + eyy * ly;
		rRightHandSideVector(17) = exz * lx + eyz * ly;

		// [0, 1, 0] : corner
		rRightHandSideVector(21) = exy * ly;
		rRightHandSideVector(22) = eyy * ly;
		rRightHandSideVector(23) = eyz * ly;

		// [0, 0, 1] : corner
		rRightHandSideVector(27) = exz * lz;
		rRightHandSideVector(28) = eyz * lz;
		rRightHandSideVector(29) = ezz * lz;

		// [1, 0, 1] : corner
		rRightHandSideVector(33) = exx * lx + exz * lz;
		rRightHandSideVector(34) = exy * lx + eyz * lz;
		rRightHandSideVector(35) = exz * lx + ezz * lz;

		// [1, 1, 1] : corner
		rRightHandSideVector(39) = exx * lx + exy * ly + exz * lz;
		rRightHandSideVector(40) = exy * lx + eyy * ly + eyz * lz;
		rRightHandSideVector(41) = exz * lx + eyz * ly + ezz * lz;

		// [0, 1, 1] : corner
		rRightHandSideVector(45) = exy * ly + exz * lz;
		rRightHandSideVector(46) = eyy * ly + eyz * lz;
		rRightHandSideVector(47) = eyz * ly + ezz * lz;	
	}

	// get current values and form the system matrix

	Vector currentValues(48);
	for(size_t i = 0; i < geom.size(); i++)
	{
		size_t index = i * 6;
		NodeType& inode = geom[i];

		currentValues(index    ) = inode.FastGetSolutionStepValue(DISPLACEMENT_X);
		currentValues(index + 1) = inode.FastGetSolutionStepValue(DISPLACEMENT_Y);
		currentValues(index + 2) = inode.FastGetSolutionStepValue(DISPLACEMENT_Z);
		currentValues(index + 3) = inode.FastGetSolutionStepValue(DISPLACEMENT_LAGRANGE_X);
		currentValues(index + 4) = inode.FastGetSolutionStepValue(DISPLACEMENT_LAGRANGE_Y);
		currentValues(index + 5) = inode.FastGetSolutionStepValue(DISPLACEMENT_LAGRANGE_Z);

		rLeftHandSideMatrix( index + 3, index     ) = 1.0;
		rLeftHandSideMatrix( index + 4, index + 1 ) = 1.0;
		rLeftHandSideMatrix( index + 5, index + 2 ) = 1.0;

		rLeftHandSideMatrix( index    , index + 3 ) = 1.0;
		rLeftHandSideMatrix( index + 1, index + 4 ) = 1.0;
		rLeftHandSideMatrix( index + 2, index + 5 ) = 1.0;
	}

	// form residual

	noalias(rRightHandSideVector) -= prod( rLeftHandSideMatrix, currentValues );
}

void RveCornerCondition3D8N::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    MatrixType dummy;
	CalculateLocalSystem(dummy, rRightHandSideVector, rCurrentProcessInfo);
}

void RveCornerCondition3D8N::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
{
	if(rResult.size() != 48) rResult.resize(48);
	GeometryType& geom = GetGeometry();
	for(size_t i = 0; i < geom.size(); i++)
	{
		size_t index = i * 6;
		NodeType& node = geom[i];
		rResult[index    ] = node.GetDof( DISPLACEMENT_X ).EquationId();
		rResult[index + 1] = node.GetDof( DISPLACEMENT_Y ).EquationId();
		rResult[index + 2] = node.GetDof( DISPLACEMENT_Z ).EquationId();
		rResult[index + 3] = node.GetDof( DISPLACEMENT_LAGRANGE_X ).EquationId();
		rResult[index + 4] = node.GetDof( DISPLACEMENT_LAGRANGE_Y ).EquationId();
		rResult[index + 5] = node.GetDof( DISPLACEMENT_LAGRANGE_Z ).EquationId();
	}
}

void RveCornerCondition3D8N::GetDofList(DofsVectorType& ConditionalDofList,ProcessInfo& CurrentProcessInfo)
{
	if(ConditionalDofList.size() != 48) ConditionalDofList.resize(48);
	GeometryType& geom = GetGeometry();
	for(size_t i = 0; i < geom.size(); i++)
	{
		size_t index = i * 6;
		NodeType& node = geom[i];
		ConditionalDofList[index    ] = node.pGetDof( DISPLACEMENT_X );
		ConditionalDofList[index + 1] = node.pGetDof( DISPLACEMENT_Y );
		ConditionalDofList[index + 2] = node.pGetDof( DISPLACEMENT_Z );
		ConditionalDofList[index + 3] = node.pGetDof( DISPLACEMENT_LAGRANGE_X );
		ConditionalDofList[index + 4] = node.pGetDof( DISPLACEMENT_LAGRANGE_Y );
		ConditionalDofList[index + 5] = node.pGetDof( DISPLACEMENT_LAGRANGE_Z );
	}
}

int RveCornerCondition3D8N::Check(const ProcessInfo& rCurrentProcessInfo)
{
	MyBase::Check(rCurrentProcessInfo);
	KRATOS_TRY

	GeometryType & geom = this->GetGeometry();
	
	if(geom.size() != 8) {
		std::stringstream ss;
		ss << "RveCornerCondition3D8N - The Geometry should have 8 nodes. Condition with ID = " << this->GetId();
		KRATOS_ERROR(std::logic_error, ss.str(), "");
	}
	
	//verify that the variables are correctly initialized
	if(DISPLACEMENT.Key() == 0)
		KRATOS_ERROR(std::invalid_argument,"DISPLACEMENT has Key zero! (check if the application is correctly registered","");

	if(DISPLACEMENT_LAGRANGE.Key() == 0)
		KRATOS_ERROR(std::invalid_argument,"DISPLACEMENT_LAGRANGE has Key zero! (check if the application is correctly registered","");

	//verify that the dofs exist
	for(unsigned int i = 0; i < geom.size(); i++)
	{
		NodeType & iNode = geom[i];
		
		if(iNode.SolutionStepsDataHas(DISPLACEMENT) == false)
			KRATOS_ERROR(std::invalid_argument,"missing variable DISPLACEMENT on node ", iNode.Id());
			
		if(iNode.HasDofFor(DISPLACEMENT_X) == false || iNode.HasDofFor(DISPLACEMENT_Y) == false || iNode.HasDofFor(DISPLACEMENT_Z) == false)
			KRATOS_ERROR(std::invalid_argument,"missing one of the dofs for the variable DISPLACEMENT on node ", iNode.Id());
			
		if(iNode.SolutionStepsDataHas(DISPLACEMENT_LAGRANGE) == false)
			KRATOS_ERROR(std::invalid_argument,"missing variable DISPLACEMENT_LAGRANGE on node ", iNode.Id());
			
		if(iNode.HasDofFor(DISPLACEMENT_LAGRANGE_X) == false || iNode.HasDofFor(DISPLACEMENT_LAGRANGE_Y) == false || iNode.HasDofFor(DISPLACEMENT_LAGRANGE_Z) == false)
			KRATOS_ERROR(std::invalid_argument,"missing one of the dofs for the variable DISPLACEMENT_LAGRANGE on node ", iNode.Id());
	}

	KRATOS_CATCH("")
	return 0;
}

}


