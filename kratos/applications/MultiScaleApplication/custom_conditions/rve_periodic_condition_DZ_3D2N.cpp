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
#include "rve_periodic_condition_DZ_3D2N.h"
#include "utilities/math_utils.h"
#include "geometries/point_3d.h"
#include "multiscale_application.h"

namespace Kratos
{
	
RvePeriodicConditionDZ3D2N::RvePeriodicConditionDZ3D2N(IndexType NewId, GeometryType::Pointer pGeometry)
	: MyBase(NewId, pGeometry)
{
}

RvePeriodicConditionDZ3D2N::RvePeriodicConditionDZ3D2N(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
	: MyBase(NewId, pGeometry, pProperties)
{
}
	
RvePeriodicConditionDZ3D2N::RvePeriodicConditionDZ3D2N(const RvePeriodicConditionDZ3D2N& rOther)
	: MyBase(rOther)
{
}

Condition::Pointer RvePeriodicConditionDZ3D2N::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
{
	return Condition::Pointer(new RvePeriodicConditionDZ3D2N(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

RvePeriodicConditionDZ3D2N::~RvePeriodicConditionDZ3D2N()
{
}

void RvePeriodicConditionDZ3D2N::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    if(rLeftHandSideMatrix.size1() != 9 || rLeftHandSideMatrix.size2() != 9) rLeftHandSideMatrix.resize(9, 9,false);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(9, 9);

    if(rRightHandSideVector.size() != 9) rRightHandSideVector.resize(9, false);
	noalias( rRightHandSideVector ) = ZeroVector(9);

	// some geometric info
	// Note:
	// to get the right sign for the distance do ALWAYS master.coordinate - slave.coordinate
	
	GeometryType& geom = GetGeometry();
	double lz = geom[1].Z0() - geom[0].Z0();

	// setup RHS vector with prescribed values

	Vector strainVector(6);
	if( this->GetMacroStrainVector(strainVector) )
	{
		double ezz = strainVector(2) * lz;
		double exz = strainVector(5) * lz * 0.5;
		double eyz = strainVector(4) * lz * 0.5;
	
		rRightHandSideVector(6) = -exz;
		rRightHandSideVector(7) = -eyz;
		rRightHandSideVector(8) = -ezz;
	}

	// get current values and form the system matrix

	Vector currentValues(9);
	currentValues(0) = geom[0].FastGetSolutionStepValue(DISPLACEMENT_X);
	currentValues(1) = geom[0].FastGetSolutionStepValue(DISPLACEMENT_Y);
	currentValues(2) = geom[0].FastGetSolutionStepValue(DISPLACEMENT_Z);
	currentValues(3) = geom[1].FastGetSolutionStepValue(DISPLACEMENT_X);
	currentValues(4) = geom[1].FastGetSolutionStepValue(DISPLACEMENT_Y);
	currentValues(5) = geom[1].FastGetSolutionStepValue(DISPLACEMENT_Z);
	currentValues(6) = geom[0].FastGetSolutionStepValue(DISPLACEMENT_LAGRANGE_X);
	currentValues(7) = geom[0].FastGetSolutionStepValue(DISPLACEMENT_LAGRANGE_Y);
	currentValues(8) = geom[0].FastGetSolutionStepValue(DISPLACEMENT_LAGRANGE_Z);
	
	rLeftHandSideMatrix(6,0) =  1.0;
	rLeftHandSideMatrix(6,3) = -1.0;
	rLeftHandSideMatrix(0,6) =  1.0;
	rLeftHandSideMatrix(3,6) = -1.0;
	rLeftHandSideMatrix(7,1) =  1.0;
	rLeftHandSideMatrix(7,4) = -1.0;
	rLeftHandSideMatrix(1,7) =  1.0;
	rLeftHandSideMatrix(4,7) = -1.0;
	rLeftHandSideMatrix(8,2) =  1.0;
	rLeftHandSideMatrix(8,5) = -1.0;
	rLeftHandSideMatrix(2,8) =  1.0;
	rLeftHandSideMatrix(5,8) = -1.0;

	// form residual

	noalias(rRightHandSideVector) -= prod( rLeftHandSideMatrix, currentValues );
}

void RvePeriodicConditionDZ3D2N::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    MatrixType dummy;
	CalculateLocalSystem(dummy, rRightHandSideVector, rCurrentProcessInfo);
}

void RvePeriodicConditionDZ3D2N::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
{
	if(rResult.size() != 9) rResult.resize(9);
	GeometryType& geom = GetGeometry();
	
	rResult[0] = geom[0].GetDof( DISPLACEMENT_X ).EquationId();
	rResult[1] = geom[0].GetDof( DISPLACEMENT_Y ).EquationId();
	rResult[2] = geom[0].GetDof( DISPLACEMENT_Z ).EquationId();
	
	rResult[3] = geom[1].GetDof( DISPLACEMENT_X ).EquationId();
	rResult[4] = geom[1].GetDof( DISPLACEMENT_Y ).EquationId();
	rResult[5] = geom[1].GetDof( DISPLACEMENT_Z ).EquationId();
	
	rResult[6] = geom[0].GetDof( DISPLACEMENT_LAGRANGE_X ).EquationId();
	rResult[7] = geom[0].GetDof( DISPLACEMENT_LAGRANGE_Y ).EquationId();
	rResult[8] = geom[0].GetDof( DISPLACEMENT_LAGRANGE_Z ).EquationId();
}

void RvePeriodicConditionDZ3D2N::GetDofList(DofsVectorType& ConditionalDofList,ProcessInfo& CurrentProcessInfo)
{
	if(ConditionalDofList.size() != 9) ConditionalDofList.resize(9);
	GeometryType& geom = GetGeometry();
	
	ConditionalDofList[0] = geom[0].pGetDof( DISPLACEMENT_X );
	ConditionalDofList[1] = geom[0].pGetDof( DISPLACEMENT_Y );
	ConditionalDofList[2] = geom[0].pGetDof( DISPLACEMENT_Z );
	
	ConditionalDofList[3] = geom[1].pGetDof( DISPLACEMENT_X );
	ConditionalDofList[4] = geom[1].pGetDof( DISPLACEMENT_Y );
	ConditionalDofList[5] = geom[1].pGetDof( DISPLACEMENT_Z );
	
	ConditionalDofList[6] = geom[0].pGetDof( DISPLACEMENT_LAGRANGE_X );
	ConditionalDofList[7] = geom[0].pGetDof( DISPLACEMENT_LAGRANGE_Y );
	ConditionalDofList[8] = geom[0].pGetDof( DISPLACEMENT_LAGRANGE_Z );
}

int RvePeriodicConditionDZ3D2N::Check(const ProcessInfo& rCurrentProcessInfo)
{
	MyBase::Check(rCurrentProcessInfo);
	KRATOS_TRY

	GeometryType & geom = this->GetGeometry();
	
	if(geom.size() != 2) {
		std::stringstream ss;
		ss << "RvePeriodicConditionDZ3D2N - The Geometry should have 2 nodes. Condition with ID = " << this->GetId();
		KRATOS_ERROR(std::logic_error, ss.str(), "");
	}
	
	//verify that the variables are correctly initialized
	if(DISPLACEMENT.Key() == 0)
		KRATOS_ERROR(std::invalid_argument,"DISPLACEMENT has Key zero! (check if the application is correctly registered","");

	if(DISPLACEMENT_LAGRANGE.Key() == 0)
		KRATOS_ERROR(std::invalid_argument,"DISPLACEMENT_LAGRANGE has Key zero! (check if the application is correctly registered","");

	//verify that the dofs exist
	for(unsigned int i = 0; i < 2; i++)
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


