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
#include "rve_corner_condition_2D4N.h"
#include "utilities/math_utils.h"
#include "geometries/point_3d.h"
#include "multiscale_application.h"

namespace Kratos
{
	
RveCornerCondition2D4N::RveCornerCondition2D4N(IndexType NewId, GeometryType::Pointer pGeometry)
	: MyBase(NewId, pGeometry)
{
}

RveCornerCondition2D4N::RveCornerCondition2D4N(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
	: MyBase(NewId, pGeometry, pProperties)
{
}
	
RveCornerCondition2D4N::RveCornerCondition2D4N(const RveCornerCondition2D4N& rOther)
	: MyBase(rOther)
{
}

Condition::Pointer RveCornerCondition2D4N::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
{
	return Condition::Pointer(new RveCornerCondition2D4N(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

RveCornerCondition2D4N::~RveCornerCondition2D4N()
{
}

void RveCornerCondition2D4N::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    if(rLeftHandSideMatrix.size1() != 16) rLeftHandSideMatrix.resize(16, 16,false);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(16, 16);

    if(rRightHandSideVector.size() != 16) rRightHandSideVector.resize(16, false);
	noalias( rRightHandSideVector ) = ZeroVector(16);

	GeometryType& geom = GetGeometry();

#ifdef RVE_TEST_MOD_CONDITIONS

	Vector strainVector(3);
	if(!this->GetMacroStrainVector(strainVector)) return;

	double exx = strainVector[0];
	double eyy = strainVector[1];
	double exy = strainVector[2] * 0.5;

	double x0 = geom[0].X0();
	double y0 = geom[0].Y0();

	// get current values and form the system matrix

	Vector currentValues(16);
	for(size_t i = 0; i < 4; i++)
	{
		size_t index = i * 4;
		NodeType& inode = geom[i];

		currentValues(index    ) = inode.FastGetSolutionStepValue(DISPLACEMENT_X);
		currentValues(index + 1) = inode.FastGetSolutionStepValue(DISPLACEMENT_Y);
		currentValues(index + 2) = inode.FastGetSolutionStepValue(DISPLACEMENT_LAGRANGE_X);
		currentValues(index + 3) = inode.FastGetSolutionStepValue(DISPLACEMENT_LAGRANGE_Y);

		rLeftHandSideMatrix( index + 2, index     ) = 1.0;
		rLeftHandSideMatrix( index + 3, index + 1 ) = 1.0;
		rLeftHandSideMatrix( index    , index + 2 ) = 1.0;
		rLeftHandSideMatrix( index + 1, index + 3 ) = 1.0;

		double xi = geom[i].X0() - x0;
		double yi = geom[i].Y0() - y0;

		rRightHandSideVector( index + 2 ) = exx * xi + exy * yi;
		rRightHandSideVector( index + 3 ) = exy * xi + eyy * yi;
	}

#else

	// setup RHS vector with prescribed values

	Vector strainVector(3);
	if( this->GetMacroStrainVector(strainVector) )
	{
		double exx = strainVector[0];
		double eyy = strainVector[1];
		double exy = strainVector[2] * 0.5;

		double lx = geom[1].X0() - geom[0].X0();
		double ly = geom[3].Y0() - geom[0].Y0();

		// [1, 0] : bottom-right corner
		rRightHandSideVector( 6) = exx * lx;
		rRightHandSideVector( 7) = exy * lx;

		// [1, 1] : top-right corner
		rRightHandSideVector(10) = exx * lx + exy * ly;
		rRightHandSideVector(11) = eyy * ly + exy * lx;

		// [0, 1] : top-left corner
		rRightHandSideVector(14) = exy * ly;
		rRightHandSideVector(15) = eyy * ly;
	}

	// get current values and form the system matrix

	Vector currentValues(16);
	for(size_t i = 0; i < 4; i++)
	{
		size_t index = i * 4;
		NodeType& inode = geom[i];

		currentValues(index    ) = inode.FastGetSolutionStepValue(DISPLACEMENT_X);
		currentValues(index + 1) = inode.FastGetSolutionStepValue(DISPLACEMENT_Y);
		currentValues(index + 2) = inode.FastGetSolutionStepValue(DISPLACEMENT_LAGRANGE_X);
		currentValues(index + 3) = inode.FastGetSolutionStepValue(DISPLACEMENT_LAGRANGE_Y);

		rLeftHandSideMatrix( index + 2, index     ) = 1.0;
		rLeftHandSideMatrix( index + 3, index + 1 ) = 1.0;
		rLeftHandSideMatrix( index    , index + 2 ) = 1.0;
		rLeftHandSideMatrix( index + 1, index + 3 ) = 1.0;
	}

#endif // RVE_TEST_MOD_CONDITIONS

	// form residual

	noalias(rRightHandSideVector) -= prod( rLeftHandSideMatrix, currentValues );
}

void RveCornerCondition2D4N::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    MatrixType dummy;
	CalculateLocalSystem(dummy, rRightHandSideVector, rCurrentProcessInfo);
}

void RveCornerCondition2D4N::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
{
	if(rResult.size() != 16) rResult.resize(16);
	GeometryType& geom = GetGeometry();
	for(size_t i = 0; i < 4; i++)
	{
		size_t index = i * 4;
		NodeType& node = geom[i];
		rResult[index    ] = node.GetDof( DISPLACEMENT_X ).EquationId();
		rResult[index + 1] = node.GetDof( DISPLACEMENT_Y ).EquationId();
		rResult[index + 2] = node.GetDof( DISPLACEMENT_LAGRANGE_X ).EquationId();
		rResult[index + 3] = node.GetDof( DISPLACEMENT_LAGRANGE_Y ).EquationId();
	}
}

void RveCornerCondition2D4N::GetDofList(DofsVectorType& ConditionalDofList,ProcessInfo& CurrentProcessInfo)
{
	if(ConditionalDofList.size() != 16) ConditionalDofList.resize(16);
	GeometryType& geom = GetGeometry();
	for(size_t i = 0; i < 4; i++)
	{
		size_t index = i * 4;
		NodeType& node = geom[i];
		ConditionalDofList[index    ] = node.pGetDof( DISPLACEMENT_X );
		ConditionalDofList[index + 1] = node.pGetDof( DISPLACEMENT_Y );
		ConditionalDofList[index + 2] = node.pGetDof( DISPLACEMENT_LAGRANGE_X );
		ConditionalDofList[index + 3] = node.pGetDof( DISPLACEMENT_LAGRANGE_Y );
	}
}

int RveCornerCondition2D4N::Check(const ProcessInfo& rCurrentProcessInfo)
{
	MyBase::Check(rCurrentProcessInfo);
	KRATOS_TRY

	GeometryType & geom = this->GetGeometry();
	
	if(geom.size() != 4) {
		std::stringstream ss;
		ss << "RveCornerCondition2D4N - The Geometry should have 4 nodes. Condition with ID = " << this->GetId();
		KRATOS_THROW_ERROR(std::logic_error, ss.str(), "");
	}
	
	//verify that the variables are correctly initialized
	if(DISPLACEMENT.Key() == 0)
		KRATOS_THROW_ERROR(std::invalid_argument,"DISPLACEMENT has Key zero! (check if the application is correctly registered","");

	if(DISPLACEMENT_LAGRANGE.Key() == 0)
		KRATOS_THROW_ERROR(std::invalid_argument,"DISPLACEMENT_LAGRANGE has Key zero! (check if the application is correctly registered","");

	//verify that the dofs exist
	for(unsigned int i = 0; i < 4; i++)
	{
		NodeType & iNode = geom[i];
		
		if(iNode.SolutionStepsDataHas(DISPLACEMENT) == false)
			KRATOS_THROW_ERROR(std::invalid_argument,"missing variable DISPLACEMENT on node ", iNode.Id());
			
		if(iNode.HasDofFor(DISPLACEMENT_X) == false || iNode.HasDofFor(DISPLACEMENT_Y) == false)
			KRATOS_THROW_ERROR(std::invalid_argument,"missing one of the dofs for the variable DISPLACEMENT on node ", iNode.Id());
			
		if(iNode.SolutionStepsDataHas(DISPLACEMENT_LAGRANGE) == false)
			KRATOS_THROW_ERROR(std::invalid_argument,"missing variable DISPLACEMENT_LAGRANGE on node ", iNode.Id());
			
		if(iNode.HasDofFor(DISPLACEMENT_LAGRANGE_X) == false || iNode.HasDofFor(DISPLACEMENT_LAGRANGE_Y) == false)
			KRATOS_THROW_ERROR(std::invalid_argument,"missing one of the dofs for the variable DISPLACEMENT_LAGRANGE on node ", iNode.Id());
	}

	KRATOS_CATCH("")
	return 0;
}

}