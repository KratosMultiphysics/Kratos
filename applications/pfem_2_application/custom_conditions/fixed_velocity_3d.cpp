// Project includes 
#include "includes/define.h"
#include "custom_conditions/fixed_velocity_3d.h"
#include "pfem_2_application_variables.h"
#include "utilities/math_utils.h"

namespace Kratos
{
	//************************************************************************************
	//************************************************************************************
	FixedVelocity3D::FixedVelocity3D(IndexType NewId, GeometryType::Pointer pGeometry)
		: Condition(NewId, pGeometry)
	{		
		//DO NOT ADD DOFS HERE!!!
	}

	//************************************************************************************
	//************************************************************************************
	FixedVelocity3D::FixedVelocity3D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
		: Condition(NewId, pGeometry, pProperties)
	{
	}
	Condition::Pointer FixedVelocity3D::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		return Condition::Pointer(new FixedVelocity3D(NewId, GetGeometry().Create(ThisNodes), pProperties));
	}

	FixedVelocity3D::~FixedVelocity3D()
	{
	}

	//************************************************************************************
	//************************************************************************************
	void FixedVelocity3D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY
		if(rCurrentProcessInfo[FRACTIONAL_STEP]==2)
		{
			KRATOS_THROW_ERROR(std::logic_error,  "method not implemented" , "");
			/*
			if(rRightHandSideVector.size() != 1)
				rRightHandSideVector.resize(1,false);

			double load = GetGeometry()[0].GetSolutionStepValue(POINT_HEAT_SOURCE);
			rRightHandSideVector[0] = load;
			*/
		}
		KRATOS_CATCH("")
	}

	//************************************************************************************
	//************************************************************************************
	void FixedVelocity3D::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY

		//we only have to calculate if 
		if(rCurrentProcessInfo[FRACTIONAL_STEP]==2)
		{
			if(rLeftHandSideMatrix.size1() != 1)
				rLeftHandSideMatrix.resize(1,1,false);
			noalias(rLeftHandSideMatrix) = ZeroMatrix(1,1);
			if(rRightHandSideVector.size() != 1)
				rRightHandSideVector.resize(1,false);
			
			array_1d<double,3> & velocity = GetGeometry()[0].GetSolutionStepValue(VELOCITY,1);
			array_1d<double,3> & normal = GetGeometry()[0].FastGetSolutionStepValue(NORMAL);
			
			double added_RHS_boundary_term=0.0;
			for (unsigned int i=0; i!=3 ; i++)
				 added_RHS_boundary_term+=velocity(i)*normal(i);
				 
			rRightHandSideVector[0] = added_RHS_boundary_term;
			
		}
		KRATOS_CATCH("")
	}


	//************************************************************************************
	//************************************************************************************
	void FixedVelocity3D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
	{
		if(rCurrentProcessInfo[FRACTIONAL_STEP]==2)
		{
			int number_of_nodes = GetGeometry().PointsNumber();
			unsigned int index;
			unsigned int dim = 1;
			rResult.resize(number_of_nodes*dim);
			for (int i=0;i<number_of_nodes;i++)
			{
				index = i*dim;
				rResult[index] = (GetGeometry()[i].GetDof(PRESSURE).EquationId());			
			}
		}
	}

	//************************************************************************************
	//************************************************************************************
	  void FixedVelocity3D::GetDofList(DofsVectorType& ConditionalDofList,ProcessInfo& rCurrentProcessInfo)
	{
		if(rCurrentProcessInfo[FRACTIONAL_STEP]==2)
		{
			unsigned int dim = 1;
			ConditionalDofList.resize(GetGeometry().size()*dim);
			unsigned int index;
			for (unsigned int i=0;i<GetGeometry().size();i++)
			{
				
				index = i*dim;
				ConditionalDofList[index] = (GetGeometry()[i].pGetDof(PRESSURE));
			}
		}
	}
} // Namespace Kratos
