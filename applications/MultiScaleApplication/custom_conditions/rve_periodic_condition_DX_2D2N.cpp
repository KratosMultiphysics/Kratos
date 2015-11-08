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
#include "rve_periodic_condition_DX_2D2N.h"
#include "utilities/math_utils.h"
#include "geometries/point_3d.h"
#include "multiscale_application.h"

namespace Kratos
{
	
RvePeriodicConditionDX2D2N::RvePeriodicConditionDX2D2N(IndexType NewId, GeometryType::Pointer pGeometry)
	: MyBase(NewId, pGeometry)
{
}

RvePeriodicConditionDX2D2N::RvePeriodicConditionDX2D2N(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
	: MyBase(NewId, pGeometry, pProperties)
{
}
	
RvePeriodicConditionDX2D2N::RvePeriodicConditionDX2D2N(const RvePeriodicConditionDX2D2N& rOther)
	: MyBase(rOther)
{
}

Condition::Pointer RvePeriodicConditionDX2D2N::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
{
	return Condition::Pointer(new RvePeriodicConditionDX2D2N(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

RvePeriodicConditionDX2D2N::~RvePeriodicConditionDX2D2N()
{
}

void RvePeriodicConditionDX2D2N::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    if(rLeftHandSideMatrix.size1() != 6 || rLeftHandSideMatrix.size2() != 6) rLeftHandSideMatrix.resize(6, 6,false);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(6, 6);

    if(rRightHandSideVector.size() != 6) rRightHandSideVector.resize(6, false);
	noalias( rRightHandSideVector ) = ZeroVector(6);

	// some geometric info
	// Note:
	// to get the right sign for the distance do ALWAYS master.coordinate - slave.coordinate
	
	GeometryType& geom = GetGeometry();
	double lx = geom[1].X0() - geom[0].X0();

	// setup RHS vector with prescribed values

	Vector strainVector(3);
	if( this->GetMacroStrainVector(strainVector) )
	{
		double exx = strainVector(0) * lx;
		double exy = strainVector(2) * lx * 0.5;
	
		rRightHandSideVector(4) = -exx;
		rRightHandSideVector(5) = -exy;
	}

	// get current values and form the system matrix

	Vector currentValues(6);
	currentValues(0) = geom[0].FastGetSolutionStepValue(DISPLACEMENT_X);
	currentValues(1) = geom[0].FastGetSolutionStepValue(DISPLACEMENT_Y);
	currentValues(2) = geom[1].FastGetSolutionStepValue(DISPLACEMENT_X);
	currentValues(3) = geom[1].FastGetSolutionStepValue(DISPLACEMENT_Y);
	currentValues(4) = geom[0].FastGetSolutionStepValue(DISPLACEMENT_LAGRANGE_X);
	currentValues(5) = geom[0].FastGetSolutionStepValue(DISPLACEMENT_LAGRANGE_Y);
	
	rLeftHandSideMatrix(4,0) =  1.0;
	rLeftHandSideMatrix(4,2) = -1.0;
	rLeftHandSideMatrix(0,4) =  1.0;
	rLeftHandSideMatrix(2,4) = -1.0;
	rLeftHandSideMatrix(5,1) =  1.0;
	rLeftHandSideMatrix(5,3) = -1.0;
	rLeftHandSideMatrix(1,5) =  1.0;
	rLeftHandSideMatrix(3,5) = -1.0;

	// form residual

	noalias(rRightHandSideVector) -= prod( rLeftHandSideMatrix, currentValues );
}

void RvePeriodicConditionDX2D2N::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    MatrixType dummy;
	CalculateLocalSystem(dummy, rRightHandSideVector, rCurrentProcessInfo);
}

void RvePeriodicConditionDX2D2N::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
{
	if(rResult.size() != 6) rResult.resize(6);
	GeometryType& geom = GetGeometry();
	
	rResult[0] = geom[0].GetDof( DISPLACEMENT_X ).EquationId();
	rResult[1] = geom[0].GetDof( DISPLACEMENT_Y ).EquationId();
	
	rResult[2] = geom[1].GetDof( DISPLACEMENT_X ).EquationId();
	rResult[3] = geom[1].GetDof( DISPLACEMENT_Y ).EquationId();
	
	rResult[4] = geom[0].GetDof( DISPLACEMENT_LAGRANGE_X ).EquationId();
	rResult[5] = geom[0].GetDof( DISPLACEMENT_LAGRANGE_Y ).EquationId();
}

void RvePeriodicConditionDX2D2N::GetDofList(DofsVectorType& ConditionalDofList,ProcessInfo& CurrentProcessInfo)
{
	if(ConditionalDofList.size() != 6) ConditionalDofList.resize(6);
	GeometryType& geom = GetGeometry();
	
	ConditionalDofList[0] = geom[0].pGetDof( DISPLACEMENT_X );
	ConditionalDofList[1] = geom[0].pGetDof( DISPLACEMENT_Y );
	
	ConditionalDofList[2] = geom[1].pGetDof( DISPLACEMENT_X );
	ConditionalDofList[3] = geom[1].pGetDof( DISPLACEMENT_Y );
	
	ConditionalDofList[4] = geom[0].pGetDof( DISPLACEMENT_LAGRANGE_X );
	ConditionalDofList[5] = geom[0].pGetDof( DISPLACEMENT_LAGRANGE_Y );
}

int RvePeriodicConditionDX2D2N::Check(const ProcessInfo& rCurrentProcessInfo)
{
	MyBase::Check(rCurrentProcessInfo);
	KRATOS_TRY

	GeometryType & geom = this->GetGeometry();
	
	if(geom.size() != 2) {
		std::stringstream ss;
		ss << "RvePeriodicConditionDX2D2N - The Geometry should have 2 nodes. Condition with ID = " << this->GetId();
		KRATOS_THROW_ERROR(std::logic_error, ss.str(), "");
	}
	
	//verify that the variables are correctly initialized
	if(DISPLACEMENT.Key() == 0)
		KRATOS_THROW_ERROR(std::invalid_argument,"DISPLACEMENT has Key zero! (check if the application is correctly registered","");

	if(DISPLACEMENT_LAGRANGE.Key() == 0)
		KRATOS_THROW_ERROR(std::invalid_argument,"DISPLACEMENT_LAGRANGE has Key zero! (check if the application is correctly registered","");

	//verify that the dofs exist
	for(unsigned int i = 0; i < 2; i++)
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