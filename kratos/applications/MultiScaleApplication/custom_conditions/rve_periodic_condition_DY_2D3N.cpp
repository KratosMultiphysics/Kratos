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
#include "rve_periodic_condition_DY_2D3N.h"
#include "utilities/math_utils.h"
#include "geometries/point_3d.h"
#include "multiscale_application.h"

namespace Kratos
{
	
RvePeriodicConditionDY2D3N::RvePeriodicConditionDY2D3N(IndexType NewId, GeometryType::Pointer pGeometry, double C1, double C2)
	: MyBase(NewId, pGeometry)
	, mC1(C1), mC2(C2)
{
}

RvePeriodicConditionDY2D3N::RvePeriodicConditionDY2D3N(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties, double C1, double C2)
	: MyBase(NewId, pGeometry, pProperties)
	, mC1(C1), mC2(C2)
{
}
	
RvePeriodicConditionDY2D3N::RvePeriodicConditionDY2D3N(const RvePeriodicConditionDY2D3N& rOther)
	: MyBase(rOther)
	, mC1(rOther.mC1), mC2(rOther.mC2)
{
}

Condition::Pointer RvePeriodicConditionDY2D3N::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
{
	return Condition::Pointer(new RvePeriodicConditionDY2D3N(NewId, GetGeometry().Create(ThisNodes), pProperties, mC1, mC2));
}

RvePeriodicConditionDY2D3N::~RvePeriodicConditionDY2D3N()
{
}

void RvePeriodicConditionDY2D3N::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    if(rLeftHandSideMatrix.size1() != 8 || rLeftHandSideMatrix.size2() != 8) rLeftHandSideMatrix.resize(8, 8,false);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(8, 8);

    if(rRightHandSideVector.size() != 8) rRightHandSideVector.resize(8, false);
	noalias( rRightHandSideVector ) = ZeroVector(8);

	// some geometric info
	// Note:
	// to get the right sign for the distance do ALWAYS master.coordinate - slave.coordinate
	
	GeometryType& geom = GetGeometry();
	double ly = geom[1].Y0() - geom[0].Y0();

	// setup RHS vector with prescribed values

	Vector strainVector(3);
	if( this->GetMacroStrainVector(strainVector) )
	{
		double eyy = strainVector(1) * ly;
		double exy = strainVector(2) * ly * 0.5;
	
		rRightHandSideVector(6) = -exy;
		rRightHandSideVector(7) = -eyy;
	}

	// get current values and form the system matrix

	Vector currentValues(8);
	currentValues(0) = geom[0].FastGetSolutionStepValue(DISPLACEMENT_X);
	currentValues(1) = geom[0].FastGetSolutionStepValue(DISPLACEMENT_Y);
	currentValues(2) = geom[1].FastGetSolutionStepValue(DISPLACEMENT_X);
	currentValues(3) = geom[1].FastGetSolutionStepValue(DISPLACEMENT_Y);
	currentValues(4) = geom[2].FastGetSolutionStepValue(DISPLACEMENT_X);
	currentValues(5) = geom[2].FastGetSolutionStepValue(DISPLACEMENT_Y);
	currentValues(6) = geom[0].FastGetSolutionStepValue(DISPLACEMENT_LAGRANGE_X);
	currentValues(7) = geom[0].FastGetSolutionStepValue(DISPLACEMENT_LAGRANGE_Y);
	
	rLeftHandSideMatrix(6,0) =  1.0;
	rLeftHandSideMatrix(6,2) = -mC1;
	rLeftHandSideMatrix(6,4) = -mC2;
	rLeftHandSideMatrix(0,6) =  1.0;
	rLeftHandSideMatrix(2,6) = -mC1;
	rLeftHandSideMatrix(4,6) = -mC2;
	rLeftHandSideMatrix(7,1) =  1.0;
	rLeftHandSideMatrix(7,3) = -mC1;
	rLeftHandSideMatrix(7,5) = -mC2;
	rLeftHandSideMatrix(1,7) =  1.0;
	rLeftHandSideMatrix(3,7) = -mC1;
	rLeftHandSideMatrix(5,7) = -mC2;

	// form residual

	noalias(rRightHandSideVector) -= prod( rLeftHandSideMatrix, currentValues );
}

void RvePeriodicConditionDY2D3N::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    MatrixType dummy;
	CalculateLocalSystem(dummy, rRightHandSideVector, rCurrentProcessInfo);
}

void RvePeriodicConditionDY2D3N::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
{
	if(rResult.size() != 8) rResult.resize(8);
	GeometryType& geom = GetGeometry();
	
	rResult[0] = geom[0].GetDof( DISPLACEMENT_X ).EquationId();
	rResult[1] = geom[0].GetDof( DISPLACEMENT_Y ).EquationId();
	
	rResult[2] = geom[1].GetDof( DISPLACEMENT_X ).EquationId();
	rResult[3] = geom[1].GetDof( DISPLACEMENT_Y ).EquationId();
	
	rResult[4] = geom[2].GetDof( DISPLACEMENT_X ).EquationId();
	rResult[5] = geom[2].GetDof( DISPLACEMENT_Y ).EquationId();
	
	rResult[6] = geom[0].GetDof( DISPLACEMENT_LAGRANGE_X ).EquationId();
	rResult[7] = geom[0].GetDof( DISPLACEMENT_LAGRANGE_Y ).EquationId();
}

void RvePeriodicConditionDY2D3N::GetDofList(DofsVectorType& ConditionalDofList,ProcessInfo& CurrentProcessInfo)
{
	if(ConditionalDofList.size() != 8) ConditionalDofList.resize(8);
	GeometryType& geom = GetGeometry();
	
	ConditionalDofList[0] = geom[0].pGetDof( DISPLACEMENT_X );
	ConditionalDofList[1] = geom[0].pGetDof( DISPLACEMENT_Y );
	
	ConditionalDofList[2] = geom[1].pGetDof( DISPLACEMENT_X );
	ConditionalDofList[3] = geom[1].pGetDof( DISPLACEMENT_Y );
	
	ConditionalDofList[4] = geom[2].pGetDof( DISPLACEMENT_X );
	ConditionalDofList[5] = geom[2].pGetDof( DISPLACEMENT_Y );
	
	ConditionalDofList[6] = geom[0].pGetDof( DISPLACEMENT_LAGRANGE_X );
	ConditionalDofList[7] = geom[0].pGetDof( DISPLACEMENT_LAGRANGE_Y );
}

int RvePeriodicConditionDY2D3N::Check(const ProcessInfo& rCurrentProcessInfo)
{
	MyBase::Check(rCurrentProcessInfo);
	KRATOS_TRY

	GeometryType & geom = this->GetGeometry();
	
	if(geom.size() != 3) {
		std::stringstream ss;
		ss << "RvePeriodicConditionDY2D3N - The Geometry should have 3 nodes. Condition with ID = " << this->GetId();
		KRATOS_ERROR(std::logic_error, ss.str(), "");
	}
	
	//verify that the variables are correctly initialized
	if(DISPLACEMENT.Key() == 0)
		KRATOS_ERROR(std::invalid_argument,"DISPLACEMENT has Key zero! (check if the application is correctly registered","");

	if(DISPLACEMENT_LAGRANGE.Key() == 0)
		KRATOS_ERROR(std::invalid_argument,"DISPLACEMENT_LAGRANGE has Key zero! (check if the application is correctly registered","");

	//verify that the dofs exist
	for(unsigned int i = 0; i < 3; i++)
	{
		NodeType & iNode = geom[i];
		
		if(iNode.SolutionStepsDataHas(DISPLACEMENT) == false)
			KRATOS_ERROR(std::invalid_argument,"missing variable DISPLACEMENT on node ", iNode.Id());
			
		if(iNode.HasDofFor(DISPLACEMENT_X) == false || iNode.HasDofFor(DISPLACEMENT_Y) == false)
			KRATOS_ERROR(std::invalid_argument,"missing one of the dofs for the variable DISPLACEMENT on node ", iNode.Id());
			
		if(iNode.SolutionStepsDataHas(DISPLACEMENT_LAGRANGE) == false)
			KRATOS_ERROR(std::invalid_argument,"missing variable DISPLACEMENT_LAGRANGE on node ", iNode.Id());
			
		if(iNode.HasDofFor(DISPLACEMENT_LAGRANGE_X) == false || iNode.HasDofFor(DISPLACEMENT_LAGRANGE_Y) == false)
			KRATOS_ERROR(std::invalid_argument,"missing one of the dofs for the variable DISPLACEMENT_LAGRANGE on node ", iNode.Id());
	}

	KRATOS_CATCH("")
	return 0;
}

}


