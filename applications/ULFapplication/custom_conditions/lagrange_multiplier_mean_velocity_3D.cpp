// Project includes 
#include "includes/define.h"
#include "custom_conditions/lagrange_multiplier_mean_velocity_3D.h"
//#include "fluid_rve_lagrange_multipliers_application.h"
#include "ULF_application.h"
#include "utilities/math_utils.h"
//#include <math.h> 

namespace Kratos
{
	//************************************************************************************
	//************************************************************************************
	MeanVelocityLagrangeMultiplierCondition3D::MeanVelocityLagrangeMultiplierCondition3D(IndexType NewId, GeometryType::Pointer pGeometry)
		: Condition(NewId, pGeometry)
	{		
		//DO NOT ADD DOFS HERE!!!
	}

	//************************************************************************************
	//************************************************************************************
	MeanVelocityLagrangeMultiplierCondition3D::MeanVelocityLagrangeMultiplierCondition3D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
		: Condition(NewId, pGeometry, pProperties)
	{
	}
	Condition::Pointer MeanVelocityLagrangeMultiplierCondition3D::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		return Condition::Pointer(new MeanVelocityLagrangeMultiplierCondition3D(NewId, GetGeometry().Create(ThisNodes), pProperties));
	}

	MeanVelocityLagrangeMultiplierCondition3D::~MeanVelocityLagrangeMultiplierCondition3D()
	{
	}

	//************************************************************************************
	//************************************************************************************
	void MeanVelocityLagrangeMultiplierCondition3D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
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
	void MeanVelocityLagrangeMultiplierCondition3D::CalculateLocalVelocityContribution(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)

	{
		KRATOS_TRY

		//we only have to calculate if 
		if(true) //(rCurrentProcessInfo[FRACTIONAL_STEP]==2)
		{
			const unsigned int dim = 3;
			if(rLeftHandSideMatrix.size1() != dim*2)
				rLeftHandSideMatrix.resize(dim*2,dim*2,false);
			noalias(rLeftHandSideMatrix) = ZeroMatrix(dim*2,dim*2);
			if(rRightHandSideVector.size() != dim*2)
				rRightHandSideVector.resize(dim*2,false);
			rRightHandSideVector = ZeroVector(dim*2);
		
		
			const array_1d<double,3> & velocity = GetGeometry()[0].GetSolutionStepValue(VELOCITY,0);
			const array_1d<double,3> & target_velocity = rCurrentProcessInfo[LAGRANGE_MULTIPLIER_VELOCITY];

			const double vol = GetGeometry()[0].GetSolutionStepValue(NODAL_AREA);

			rLeftHandSideMatrix(0,3)= cbrt(vol);
			rLeftHandSideMatrix(3,0)= cbrt(vol);
			rRightHandSideVector(3)= cbrt(vol)*target_velocity[0];
			rLeftHandSideMatrix(1,4)= cbrt(vol);
			rLeftHandSideMatrix(4,1)= cbrt(vol);
			rRightHandSideVector(4)= cbrt(vol)*target_velocity[1];
			rLeftHandSideMatrix(2,5)= cbrt(vol);
			rLeftHandSideMatrix(5,2)= cbrt(vol);
			rRightHandSideVector(5)= cbrt(vol)*target_velocity[2];

			//substracting:
			array_1d<double,dim*2>  temp_vector=ZeroVector(dim*2);
			temp_vector[0]=velocity[0];
			temp_vector[1]=velocity[1];
			temp_vector[2]=velocity[2];
			const array_1d<double,3> & velocity_aux_node = (this->GetValue(NEIGHBOUR_NODES)[0]).GetSolutionStepValue(LAGRANGE_MULTIPLIER_VELOCITY,0);
			temp_vector[3]=velocity_aux_node[0];
			temp_vector[4]=velocity_aux_node[1];
			temp_vector[5]=velocity_aux_node[2];

			//KRATOS_WATCH(rLeftHandSideMatrix);
			//KRATOS_WATCH(rRightHandSideVector);

			noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,temp_vector);
			
		}
		KRATOS_CATCH("")
	}


	void MeanVelocityLagrangeMultiplierCondition3D::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)

	{
		KRATOS_TRY

		//we only have to calculate if 
		if(true) //(rCurrentProcessInfo[FRACTIONAL_STEP]==2)
		{
			const unsigned int dim = 3;
			if(rLeftHandSideMatrix.size1() != dim*2)
				rLeftHandSideMatrix.resize(dim*2,dim*2,false);
			noalias(rLeftHandSideMatrix) = ZeroMatrix(dim*2,dim*2);
			if(rRightHandSideVector.size() != dim*2)
				rRightHandSideVector.resize(dim*2,false);
			rRightHandSideVector = ZeroVector(dim*2);
			

	
			
		}
		KRATOS_CATCH("")
	}


	//************************************************************************************
	//************************************************************************************
	void MeanVelocityLagrangeMultiplierCondition3D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
	{
		if(true)//rCurrentProcessInfo[FRACTIONAL_STEP]==2)
		{
			unsigned int dim = 3;
            if (rResult.size()!=(2*dim))
			    rResult.resize(2*dim);
			rResult[0] = (GetGeometry()[0].GetDof(VELOCITY_X).EquationId());
			rResult[1] = (GetGeometry()[0].GetDof(VELOCITY_Y).EquationId());
			rResult[2] = (GetGeometry()[0].GetDof(VELOCITY_Z).EquationId());
			rResult[3] = (this->GetValue(NEIGHBOUR_NODES)[0]).GetDof(LAGRANGE_MULTIPLIER_VELOCITY_X).EquationId();			
			rResult[4] = (this->GetValue(NEIGHBOUR_NODES)[0]).GetDof(LAGRANGE_MULTIPLIER_VELOCITY_Y).EquationId();
			rResult[5] = (this->GetValue(NEIGHBOUR_NODES)[0]).GetDof(LAGRANGE_MULTIPLIER_VELOCITY_Z).EquationId();						
		}
	}

	//************************************************************************************
	//************************************************************************************
	  void MeanVelocityLagrangeMultiplierCondition3D::GetDofList(DofsVectorType& ConditionalDofList,ProcessInfo& rCurrentProcessInfo)
	{
		if(true)//rCurrentProcessInfo[FRACTIONAL_STEP]==2)
		{
			unsigned int dim = 3;
            if (ConditionalDofList.size()!=2*dim)
			     ConditionalDofList.resize(2*dim);
			ConditionalDofList[0] = GetGeometry()[0].pGetDof(VELOCITY_X);
			ConditionalDofList[1] = GetGeometry()[0].pGetDof(VELOCITY_Y);
			ConditionalDofList[2] = GetGeometry()[0].pGetDof(VELOCITY_Z);
			ConditionalDofList[3] = (this->GetValue(NEIGHBOUR_NODES)[0]).pGetDof(LAGRANGE_MULTIPLIER_VELOCITY_X);
			ConditionalDofList[4] = (this->GetValue(NEIGHBOUR_NODES)[0]).pGetDof(LAGRANGE_MULTIPLIER_VELOCITY_Y);
			ConditionalDofList[5] = (this->GetValue(NEIGHBOUR_NODES)[0]).pGetDof(LAGRANGE_MULTIPLIER_VELOCITY_Z);
		}
	}
} // Namespace Kratos
