// Project includes 
#include "includes/define.h"
#include "custom_conditions/lagrange_multiplier_mean_velocity_2d.h"
#include "fluid_rve_lagrange_multipliers_application.h"
#include "utilities/math_utils.h"

namespace Kratos
{
	//************************************************************************************
	//************************************************************************************
	MeanVelocityLagrangeMultiplierCondition2D::MeanVelocityLagrangeMultiplierCondition2D(IndexType NewId, GeometryType::Pointer pGeometry)
		: Condition(NewId, pGeometry)
	{		
		//DO NOT ADD DOFS HERE!!!
	}

	//************************************************************************************
	//************************************************************************************
	MeanVelocityLagrangeMultiplierCondition2D::MeanVelocityLagrangeMultiplierCondition2D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
		: Condition(NewId, pGeometry, pProperties)
	{
	}
	Condition::Pointer MeanVelocityLagrangeMultiplierCondition2D::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		return Condition::Pointer(new MeanVelocityLagrangeMultiplierCondition2D(NewId, GetGeometry().Create(ThisNodes), pProperties));
	}

	MeanVelocityLagrangeMultiplierCondition2D::~MeanVelocityLagrangeMultiplierCondition2D()
	{
	}

	//************************************************************************************
	//************************************************************************************
	void MeanVelocityLagrangeMultiplierCondition2D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY
		/*
		if(rCurrentProcessInfo[FRACTIONAL_STEP]==2)
		{
			KRATOS_THROW_ERROR(std::logic_error,  "method not implemented" , "");
		}
		*/
		KRATOS_THROW_ERROR(std::logic_error,  "method not implemented" , "");
		
		KRATOS_CATCH("")
	}

	//************************************************************************************
	//************************************************************************************
	void MeanVelocityLagrangeMultiplierCondition2D::CalculateLocalVelocityContribution(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)

	{
		KRATOS_TRY

		//we only have to calculate if 
		if(true) //(rCurrentProcessInfo[FRACTIONAL_STEP]==2)
		{
			const unsigned int dim = 2;
			if(rLeftHandSideMatrix.size1() != dim*2)
				rLeftHandSideMatrix.resize(dim*2,dim*2,false);
			noalias(rLeftHandSideMatrix) = ZeroMatrix(dim*2,dim*2);
			if(rRightHandSideVector.size() != dim*2)
				rRightHandSideVector.resize(dim*2,false);
			rRightHandSideVector = ZeroVector(dim*2);
		
		
			const array_1d<double,3> & velocity = GetGeometry()[0].GetSolutionStepValue(VELOCITY,0);
			const array_1d<double,3> & target_velocity = rCurrentProcessInfo[LAGRANGE_MULTIPLIER_VELOCITY];

			const double area = GetGeometry()[0].GetSolutionStepValue(NODAL_AREA);

			rLeftHandSideMatrix(0,2)= sqrt(area);
			rLeftHandSideMatrix(2,0)= sqrt(area);
			rRightHandSideVector(2)= sqrt(area)*target_velocity[0];
			rLeftHandSideMatrix(1,3)= sqrt(area);
			rLeftHandSideMatrix(3,1)= sqrt(area);
			rRightHandSideVector(3)= sqrt(area)*target_velocity[1];

			//substracting:
			array_1d<double,dim*2>  temp_vector=ZeroVector(dim*2);
			temp_vector[0]=velocity[0];
			temp_vector[1]=velocity[1];
			const array_1d<double,3> & velocity_aux_node = (this->GetValue(NEIGHBOUR_NODES)[0]).GetSolutionStepValue(LAGRANGE_MULTIPLIER_VELOCITY,0);
			temp_vector[2]=velocity_aux_node[0];
			temp_vector[3]=velocity_aux_node[1];

			//KRATOS_WATCH(rLeftHandSideMatrix);
			//KRATOS_WATCH(rRightHandSideVector);

			noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,temp_vector);
			
		}
		KRATOS_CATCH("")
	}


	void MeanVelocityLagrangeMultiplierCondition2D::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)

	{
		KRATOS_TRY

		//we only have to calculate if 
		if(true) //(rCurrentProcessInfo[FRACTIONAL_STEP]==2)
		{
			if(rLeftHandSideMatrix.size1() != 3)
				rLeftHandSideMatrix.resize(3,3,false);
			noalias(rLeftHandSideMatrix) = ZeroMatrix(3,3);
			if(rRightHandSideVector.size() != 3)
				rRightHandSideVector.resize(3,false);
			rRightHandSideVector = ZeroVector(3);
			

	
			
		}
		KRATOS_CATCH("")
	}


	//************************************************************************************
	//************************************************************************************
	void MeanVelocityLagrangeMultiplierCondition2D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
	{
		if(true)//rCurrentProcessInfo[FRACTIONAL_STEP]==2)
		{
            if (rResult.size()!=(3))
			    rResult.resize(3);
			rResult[0] = (GetGeometry()[0].GetDof(VELOCITY_X).EquationId());
			rResult[1] = (GetGeometry()[0].GetDof(VELOCITY_Y).EquationId());
			rResult[2] = (this->GetValue(NEIGHBOUR_NODES)[0]).GetDof(LAGRANGE_MULTIPLIER_VELOCITY_X).EquationId();			
			//rResult[3] = (this->GetValue(NEIGHBOUR_NODES)[0]).GetDof(LAGRANGE_MULTIPLIER_VELOCITY_Y).EquationId();			
		}
	}

	//************************************************************************************
	//************************************************************************************
	  void MeanVelocityLagrangeMultiplierCondition2D::GetDofList(DofsVectorType& ConditionalDofList,ProcessInfo& rCurrentProcessInfo)
	{
		if(true)//rCurrentProcessInfo[FRACTIONAL_STEP]==2)
		{
			unsigned int dim = 3;
            if (ConditionalDofList.size()!=3)
			     ConditionalDofList.resize(3);
			ConditionalDofList[0] = GetGeometry()[0].pGetDof(VELOCITY_X);
			ConditionalDofList[1] = GetGeometry()[0].pGetDof(VELOCITY_Y);
			ConditionalDofList[2] = (this->GetValue(NEIGHBOUR_NODES)[0]).pGetDof(LAGRANGE_MULTIPLIER_VELOCITY_X);
			//ConditionalDofList[3] = (this->GetValue(NEIGHBOUR_NODES)[0]).pGetDof(LAGRANGE_MULTIPLIER_VELOCITY_Y);
		}
	}
} // Namespace Kratos
