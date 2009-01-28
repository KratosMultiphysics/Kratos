//   
//   Project Name:        Kratos       
//   Last modified by:    $Author: rrossi $
//   Date:                $Date: 2007-01-23 15:38:00 $
//   Revision:            $Revision: 1.1 $
//
// 


// System includes 


// External includes 


// Project includes 
#include "includes/define.h"
#include "custom_conditions/pressure_cond.h"
#include "utilities/math_utils.h"
#include "PFEM_application.h"

namespace Kratos
{
	//************************************************************************************
	//************************************************************************************
	PressureCond::PressureCond(IndexType NewId, GeometryType::Pointer pGeometry)
		: Condition(NewId, pGeometry)
	{		
		//DO NOT ADD DOFS HERE!!!
	}

	//************************************************************************************
	//************************************************************************************
	PressureCond::PressureCond(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
		: Condition(NewId, pGeometry, pProperties)
	{
	}

	Condition::Pointer PressureCond::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		return Condition::Pointer(new PressureCond(NewId, GetGeometry().Create(ThisNodes), pProperties));
	}

	PressureCond::~PressureCond()
	{
	}


	//************************************************************************************
	//************************************************************************************
	void PressureCond::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{

		if(rRightHandSideVector.size() != 2)
            rRightHandSideVector.resize(2,false);
		//there can be a problem if I make the size of rHS vector = 2, coz the force has 3 components - comp. in Z =0
		const array_1d<double,3>& fp = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE_FORCE);

		rRightHandSideVector[0] = -fp[0];
		rRightHandSideVector[1] = -fp[1];
		
	}

	//************************************************************************************
	//************************************************************************************
	void PressureCond::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		if(rLeftHandSideMatrix.size1() != 2)
		{
			rLeftHandSideMatrix.resize(2,2,false);
		}
		noalias(rLeftHandSideMatrix) = ZeroMatrix(2,2); 

		if(rRightHandSideVector.size() != 2)
            rRightHandSideVector.resize(2,false);

		const array_1d<double,3>& fp = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE_FORCE);
		
		
		rRightHandSideVector[0] = -fp[0];
		rRightHandSideVector[1] = -fp[1];


	}



	//************************************************************************************
	//************************************************************************************
	void PressureCond::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
	{
		if(rResult.size() != 2)
			rResult.resize(2);

		rResult[0] = (GetGeometry()[0].GetDof(DISPLACEMENT_X)).EquationId();
		rResult[1] = (GetGeometry()[0].GetDof(DISPLACEMENT_Y).EquationId());

	}

	//************************************************************************************
	//************************************************************************************
	  void PressureCond::GetDofList(DofsVectorType& ConditionalDofList,ProcessInfo& CurrentProcessInfo)
	{
		if(ConditionalDofList.size() != 2)
			ConditionalDofList.resize(2);

		ConditionalDofList[0] = (GetGeometry()[0].pGetDof(DISPLACEMENT_X));
		ConditionalDofList[1] = (GetGeometry()[0].pGetDof(DISPLACEMENT_Y));
		
	}
} // Namespace Kratos


