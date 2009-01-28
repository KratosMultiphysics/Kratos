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
#include "custom_conditions/pressure_cond3d.h"
#include "utilities/math_utils.h"
#include "PFEM_application.h"

namespace Kratos
{
	//************************************************************************************
	//************************************************************************************
	PressureCond3D::PressureCond3D(IndexType NewId, GeometryType::Pointer pGeometry)
		: Condition(NewId, pGeometry)
	{		
		//DO NOT ADD DOFS HERE!!!
	}

	//************************************************************************************
	//************************************************************************************
	PressureCond3D::PressureCond3D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
		: Condition(NewId, pGeometry, pProperties)
	{
	}

	Condition::Pointer PressureCond3D::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		return Condition::Pointer(new PressureCond3D(NewId, GetGeometry().Create(ThisNodes), pProperties));
	}

	PressureCond3D::~PressureCond3D()
	{
	}


	//************************************************************************************
	//************************************************************************************
	void PressureCond3D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{

		if(rRightHandSideVector.size() != 3)
            rRightHandSideVector.resize(3,false);
		//there can be a problem if I make the size of rHS vector = 2, coz the force has 3 components - comp. in Z =0
		const array_1d<double,3>& fp = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE_FORCE);

		rRightHandSideVector[0] = -fp[0];
		rRightHandSideVector[1] = -fp[1];
		rRightHandSideVector[2] = -fp[2];
		
	}

	//************************************************************************************
	//************************************************************************************
	void PressureCond3D::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		if(rLeftHandSideMatrix.size1() != 3)
		{
			rLeftHandSideMatrix.resize(3,3,false);
		}
		noalias(rLeftHandSideMatrix) = ZeroMatrix(3,3); 

		if(rRightHandSideVector.size() != 3)
            rRightHandSideVector.resize(3,false);

		const array_1d<double,3>& fp = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE_FORCE);
		
		
		rRightHandSideVector[0] = -fp[0];
		rRightHandSideVector[1] = -fp[1];
		rRightHandSideVector[2] = -fp[2];

	}



	//************************************************************************************
	//************************************************************************************
	void PressureCond3D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
	{
		if(rResult.size() != 3)
			rResult.resize(3);

		rResult[0] = (GetGeometry()[0].GetDof(DISPLACEMENT_X)).EquationId();
		rResult[1] = (GetGeometry()[0].GetDof(DISPLACEMENT_Y).EquationId());
		rResult[2] = (GetGeometry()[0].GetDof(DISPLACEMENT_Z).EquationId());
	}

	//************************************************************************************
	//************************************************************************************
	  void PressureCond3D::GetDofList(DofsVectorType& ConditionalDofList,ProcessInfo& CurrentProcessInfo)
	{
		if(ConditionalDofList.size() != 3)
			ConditionalDofList.resize(3);

		ConditionalDofList[0] = (GetGeometry()[0].pGetDof(DISPLACEMENT_X));
		ConditionalDofList[1] = (GetGeometry()[0].pGetDof(DISPLACEMENT_Y));
		ConditionalDofList[2] = (GetGeometry()[0].pGetDof(DISPLACEMENT_Z));
		
	}
} // Namespace Kratos


