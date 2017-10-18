// Project includes 
#include "includes/define.h"
#include "custom_conditions/inverse_tangent_velocity_periodic_condition_2d.h"
#include "fluid_rve_lagrange_multipliers_application.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h" 


namespace Kratos
{
	//************************************************************************************
	//************************************************************************************
	InverseTangentVelocityPeriodicCondition2D2N::InverseTangentVelocityPeriodicCondition2D2N(IndexType NewId, GeometryType::Pointer pGeometry)
		: Condition(NewId, pGeometry)
	{		
		//DO NOT ADD DOFS HERE!!!
	}

	//************************************************************************************
	//************************************************************************************
	InverseTangentVelocityPeriodicCondition2D2N::InverseTangentVelocityPeriodicCondition2D2N(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
		: Condition(NewId, pGeometry, pProperties)
	{
	}
	Condition::Pointer InverseTangentVelocityPeriodicCondition2D2N::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		return Condition::Pointer(new InverseTangentVelocityPeriodicCondition2D2N(NewId, GetGeometry().Create(ThisNodes), pProperties));
	}

	InverseTangentVelocityPeriodicCondition2D2N::~InverseTangentVelocityPeriodicCondition2D2N()
	{
	}

	//************************************************************************************
	//************************************************************************************
	void InverseTangentVelocityPeriodicCondition2D2N::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
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
	/*
	void InverseTangentVelocityPeriodicCondition2D2N::CalculateLocalVelocityContribution(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)

	{
		KRATOS_TRY

		//we only have to calculate if 
		if(rLeftHandSideMatrix.size1() != 3) // two dofs per node + 3 multipliers
			rLeftHandSideMatrix.resize(3,3,false);
		noalias(rLeftHandSideMatrix) = ZeroMatrix(3,3);
		if(rRightHandSideVector.size() != 3)
			rRightHandSideVector.resize(3,false);
		rRightHandSideVector = ZeroVector(3);
		
		//double Area;
		//array_1d<double, 3 > N;
		//boost::numeric::ublas::bounded_matrix<double, 3, 2> DN_DX;
		boost::numeric::ublas::bounded_matrix<double, 1, 3 > LM_matrix = ZeroMatrix(1, 3) ; //(gradient)
		Geometry<Node<3> >& geom = this->GetGeometry();
		//GeometryUtils::CalculateGeometryData(geom, DN_DX, N, Area);
			

		LM_matrix(0, 0 ) = 1.0; // du/dy : we use DNi/DNy
		LM_matrix(0, 1 ) =  1.0; // dv/dx : we use DNi/DNx
			
			
		for (unsigned int i = 0; i < 2; i++) //velocity dofs
		{
			for (unsigned int j = 0; j < 1 ; j++) // LM	
			{
				rLeftHandSideMatrix(i,j+2) = LM_matrix(j,i);
				rLeftHandSideMatrix(j+2,i) = LM_matrix(j,i);
			}
		}

			
		//substracting:
		array_1d<double,3>  temp_vector=ZeroVector(3);
		for (unsigned int j = 0; j < 2; j++) //2 nodes
		{
			const array_1d<double,3> & velocity = GetGeometry()[j].GetSolutionStepValue(VELOCITY,0);
			if(fabs(GetGeometry()[0].X()-GetGeometry()[1].X())<0.00001)
				temp_vector[j]=velocity[0];
			else
				temp_vector[j]=velocity[1];

		}
		temp_vector[2]=GetGeometry()[0].GetSolutionStepValue(LAGRANGE_MULTIPLIER_TANGENT_VELOCITY,0);
		
 		noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,temp_vector);
			



		KRATOS_CATCH("")
	}
	*/
	
	void InverseTangentVelocityPeriodicCondition2D2N::CalculateLocalVelocityContribution(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)

	{
		KRATOS_TRY

		//we only have to calculate if 
		if(rLeftHandSideMatrix.size1() != 3) // two dofs per node + 3 multipliers
			rLeftHandSideMatrix.resize(3,3,false);
		noalias(rLeftHandSideMatrix) = ZeroMatrix(3,3);
		if(rRightHandSideVector.size() != 3)
			rRightHandSideVector.resize(3,false);
		rRightHandSideVector = ZeroVector(3);
		
		//double Area;
		//array_1d<double, 3 > N;
		//boost::numeric::ublas::bounded_matrix<double, 3, 2> DN_DX;
		boost::numeric::ublas::bounded_matrix<double, 1, 3 > LM_matrix = ZeroMatrix(1, 3) ; //(gradient)
		Geometry<Node<3> >& geom = this->GetGeometry();
		//GeometryUtils::CalculateGeometryData(geom, DN_DX, N, Area);
			

		LM_matrix(0, 0 ) = 1.0; // du/dy : we use DNi/DNy
		LM_matrix(0, 1 ) = -1.0; // dv/dx : we use DNi/DNx
			
			
		for (unsigned int i = 0; i < 2; i++) //velocity dofs
		{
			for (unsigned int j = 0; j < 1 ; j++) // LM	
			{
				rLeftHandSideMatrix(i,j+2) = LM_matrix(j,i);
				rLeftHandSideMatrix(j+2,i) = LM_matrix(j,i);
			}
		}

		//rRightHandSideVector[2] = -20.0;
		rRightHandSideVector[2] = -80.0;
			
		//substracting:
		array_1d<double,3>  temp_vector=ZeroVector(3);
		for (unsigned int j = 0; j < 2; j++) //2 nodes
		{
			const array_1d<double,3> & velocity = GetGeometry()[j].GetSolutionStepValue(VELOCITY,0);
			if(fabs(GetGeometry()[0].X()-GetGeometry()[1].X())<0.00001)
				temp_vector[j]=velocity[0];
			else
				temp_vector[j]=velocity[1];

		}
		temp_vector[2]=GetGeometry()[0].GetSolutionStepValue(LAGRANGE_MULTIPLIER_TANGENT_VELOCITY,0);
		
 		noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,temp_vector);
			



		KRATOS_CATCH("")
	}
	

	void InverseTangentVelocityPeriodicCondition2D2N::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)

	{
		KRATOS_TRY

		//we only have to calculate if 

		if(rLeftHandSideMatrix.size1() != 3)
			rLeftHandSideMatrix.resize(3,3,false);
		noalias(rLeftHandSideMatrix) = ZeroMatrix(3,3);
		if(rRightHandSideVector.size() != 3)
			rRightHandSideVector.resize(3,false);
		rRightHandSideVector = ZeroVector(3);

		KRATOS_CATCH("")
	}


	//************************************************************************************
	//************************************************************************************
	void InverseTangentVelocityPeriodicCondition2D2N::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
	{
		 if (rResult.size()!=(3))
			    rResult.resize(3);
	     unsigned int LocalIndex = 0;      
		 if(fabs(GetGeometry()[0].X()-GetGeometry()[1].X())<0.00001)// && GetGeometry()[0].GetSolutionStepValue(NODE_PAIR_X_COMPONENT)!=0 ) //then this is the condition relating top and bottom nodes that needs that is not in the corner.
		 { 
			 for (unsigned int iNode = 0; iNode < 2; ++iNode)
			 {
				rResult[LocalIndex++] = (GetGeometry()[iNode].GetDof(VELOCITY_X).EquationId());
			 }
			 rResult[LocalIndex++] = GetGeometry()[0].GetDof(LAGRANGE_MULTIPLIER_TANGENT_VELOCITY).EquationId();
		 }
		 else if(fabs(GetGeometry()[0].Y()-GetGeometry()[1].Y())<0.00001)// && GetGeometry()[0].GetSolutionStepValue(NODE_PAIR_Y_COMPONENT)!=0 ) //then this is the condition relating left and right nodes that needs that is not in the corner.
		 {
			for (unsigned int iNode = 0; iNode < 2; ++iNode)
			 {
				rResult[LocalIndex++] = (GetGeometry()[iNode].GetDof(VELOCITY_Y).EquationId());
			 }
			 rResult[LocalIndex++] = GetGeometry()[0].GetDof(LAGRANGE_MULTIPLIER_TANGENT_VELOCITY).EquationId();
		 }
	}

	//************************************************************************************
	//************************************************************************************
	  void InverseTangentVelocityPeriodicCondition2D2N::GetDofList(DofsVectorType& ConditionalDofList,ProcessInfo& rCurrentProcessInfo)
	{

         if (ConditionalDofList.size()!=3)
			  ConditionalDofList.resize(3);
	     unsigned int LocalIndex = 0;      
		 if(fabs(GetGeometry()[0].X()-GetGeometry()[1].X())<0.00001 )// && GetGeometry()[0].GetSolutionStepValue(NODE_PAIR_X_COMPONENT)!=0 ) //then this is the condition relating top and bottom nodes that needs that is not in the corner.
		 { 
			 for (unsigned int iNode = 0; iNode < 2; ++iNode)
			 {
				ConditionalDofList[LocalIndex++] = GetGeometry()[iNode].pGetDof(VELOCITY_X);
			 }
			 ConditionalDofList[LocalIndex++] = GetGeometry()[0].pGetDof(LAGRANGE_MULTIPLIER_TANGENT_VELOCITY);
		 }
		 else if(fabs(GetGeometry()[0].Y()-GetGeometry()[1].Y())<0.00001 )// && GetGeometry()[0].GetSolutionStepValue(NODE_PAIR_Y_COMPONENT)!=0 ) //then this is the condition relating left and right nodes that needs that is not in the corner.
		 {
			 for (unsigned int iNode = 0; iNode < 2; ++iNode)
			 {
				ConditionalDofList[LocalIndex++] = GetGeometry()[iNode].pGetDof(VELOCITY_Y);
			 }
			 ConditionalDofList[LocalIndex++] = GetGeometry()[0].pGetDof(LAGRANGE_MULTIPLIER_TANGENT_VELOCITY);
		 }
	}
} // Namespace Kratos
