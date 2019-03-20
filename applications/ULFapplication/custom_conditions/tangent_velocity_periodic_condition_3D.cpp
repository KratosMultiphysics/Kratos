// Project includes 
#include "includes/define.h"
#include "custom_conditions/tangent_velocity_periodic_condition_3D.h"
//#include "fluid_rve_lagrange_multipliers_application.h"
#include "ULF_application.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h" 

//author: Pavel Ryzhakov

namespace Kratos
{
	//************************************************************************************
	//************************************************************************************
	TangentVelocityPeriodicCondition3D2N::TangentVelocityPeriodicCondition3D2N(IndexType NewId, GeometryType::Pointer pGeometry)
		: Condition(NewId, pGeometry)
	{		
		//DO NOT ADD DOFS HERE!!!
	}

	//************************************************************************************
	//************************************************************************************
	TangentVelocityPeriodicCondition3D2N::TangentVelocityPeriodicCondition3D2N(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
		: Condition(NewId, pGeometry, pProperties)
	{
	}
	Condition::Pointer TangentVelocityPeriodicCondition3D2N::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		return Condition::Pointer(new TangentVelocityPeriodicCondition3D2N(NewId, GetGeometry().Create(ThisNodes), pProperties));
	}

	TangentVelocityPeriodicCondition3D2N::~TangentVelocityPeriodicCondition3D2N()
	{
	}

	//************************************************************************************
	//************************************************************************************
	void TangentVelocityPeriodicCondition3D2N::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY

		KRATOS_THROW_ERROR(std::logic_error,  "method not implemented" , "");
		
		KRATOS_CATCH("")
	}

	//************************************************************************************
	//************************************************************************************	
	void TangentVelocityPeriodicCondition3D2N::CalculateLocalVelocityContribution(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)

	{
		KRATOS_TRY
		//we only have to calculate if 
		if(rLeftHandSideMatrix.size1() != 6) // v1* v2* v1** v2** and lambda (these can be: v1x v2x, v1y v2y and lambda or v1x v2x v1z v2z and lambda etc... two tangential components per wall
			rLeftHandSideMatrix.resize(6,6,false);
		noalias(rLeftHandSideMatrix) = ZeroMatrix(6,6);
		if(rRightHandSideVector.size() != 6)
			rRightHandSideVector.resize(6,false);
		rRightHandSideVector = ZeroVector(6);
		
		boost::numeric::ublas::bounded_matrix<double, 6, 6 > LM_matrix = ZeroMatrix(6, 6) ; //(gradient)
		Geometry<Node<3> >& geom = this->GetGeometry();		

		//READING THE G JUMP VECTOR: CONTAINS Gxy, Gxz, Gyz // IT IS SAVED IN THE AUX_VECTOR VARIABLE
		array_1d<double,3>  jump_vector= rCurrentProcessInfo[AUX_VECTOR];  //ZeroVector(3);
		
		///////////////////////////////////////////////////////////

		//LM_matrix(0, 0 ) = 1.0; 
		//LM_matrix(0, 1 ) = -jump_vector[0]; 
		//LM is like this
		// 0   0   0   0   1   0
		// 0   0   0   0   0   1
		// 0   0   0   0  -1   0
		// 0   0   0   0   0  -1
		// 1   0  -1   0   0   0
		// 0   1   0  -1   0   0
		
		//to be multiplied by v1* v1** v2* v2** lambda

		LM_matrix(0,4)=1.0; LM_matrix(1,5)=1.0; LM_matrix(2,4)=-1.0; LM_matrix(3,5)=-1.0; LM_matrix(4,0)=1.0; LM_matrix(4,2)=-1.0; LM_matrix(5,1)=1.0;  LM_matrix(5,3)=-1.0;
		

		//because we multiply it afterwards by vi vj lambda, where vi and vj are the velocities of the "paired" nodes i and j
		rLeftHandSideMatrix=LM_matrix;

		//WE ARE FIRST TRYING THE PROBLEMS WHERE PERIODICITY WITH JUMPS IS PRESENT IN LEFT-RIGHT and UP-LOW wallss in X and Y components only. 
		const array_1d<double,3> & velocity0 = GetGeometry()[0].GetSolutionStepValue(VELOCITY,0);
		const array_1d<double,3> & velocity1 = GetGeometry()[1].GetSolutionStepValue(VELOCITY,0);
		//substracting:
		array_1d<double,6>  temp_vector=ZeroVector(6);
		//UPPER - LOWER PAIR
		//"paired" nodes of upper and lower walls: they have same x and same z... => tangentials are: x component vel[0] and z component vel[2]
		if(fabs(GetGeometry()[0].X()-GetGeometry()[1].X())<0.00001 && fabs(GetGeometry()[0].Z()-GetGeometry()[1].Z())<0.00001)
			{

			//x_component
			temp_vector[0]=velocity0[0];
			temp_vector[2]=velocity1[0];
			//z_component
			temp_vector[1]=velocity0[2];
			temp_vector[3]=velocity1[2];
			//primera componente tangente -> xy
			rRightHandSideVector[4] = -jump_vector[0];		//v1x*-v2x*=0
			//rRightHandSideVector[4] = -jump_vector[0]+(jump_vector[0]*diminish_jump_value);
			//segunda componente tangente -> yz			
			rRightHandSideVector[5] = -jump_vector[2];		//v1z**-v2z**=0			
			}

		//"paired" nodes of left and rights walls: they have same y and same z... => tangentials are: y component vel[1] and z component vel[2]
		if(fabs(GetGeometry()[0].Y()-GetGeometry()[1].Y())<0.00001 && fabs(GetGeometry()[0].Z()-GetGeometry()[1].Z())<0.00001)
			{
			//y_component
			temp_vector[0]=velocity0[1];
			temp_vector[2]=velocity1[1];
			//z_component
			temp_vector[1]=velocity0[2];
			temp_vector[3]=velocity1[2];
			//primera componente tangente -> xy
			rRightHandSideVector[4] = -jump_vector[0];		//v1y*-v2y*=0			
			//segunda componente tangente -> yz
			rRightHandSideVector[5] = -jump_vector[2];		//v1z**-v2z**=0			
			}


		//"paired" nodes of front and back walls: they have same x and same y... => tangentials are: x component vel[0] and y component vel[1]
		if(fabs(GetGeometry()[0].X()-GetGeometry()[1].X())<0.00001 && fabs(GetGeometry()[0].Y()-GetGeometry()[1].Y())<0.00001)
			{		
			//y_component
			temp_vector[0]=velocity0[0];
			temp_vector[2]=velocity1[0];
			//z_component
			temp_vector[1]=velocity0[1];
			temp_vector[3]=velocity1[1];
			//primera componente tangente -> xz			
			rRightHandSideVector[4] = -jump_vector[1];		//v1x*-v2x*=0
			//segunda componente tangente -> yz
			rRightHandSideVector[5] = -jump_vector[2];		//v1y**-v2y**=0
			//below is the one with a jump
			//rRightHandSideVector[5] = -jump_vector[......];
			}

		
		//so, depending on whether we are haveing nodes on horiz or vertical wall, we put the horiz or vert component of the paired nodes velocities 
		temp_vector[4]=GetGeometry()[0].GetSolutionStepValue(LAGRANGE_MULTIPLIER_TANGENT_VELOCITY,0);
		temp_vector[5]=GetGeometry()[0].GetSolutionStepValue(LAGRANGE_MULTIPLIER_SECOND_TANGENT_VELOCITY,0);
		//so temp_vector contains: [vxi vxj lambda] if these are horizontal (upper lower) walls and [vyi vyj lambda]  if thats a vertical wall
		//since this is a condition for tangential components
		
 		noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,temp_vector);
			



		KRATOS_CATCH("")
	}
	

	void TangentVelocityPeriodicCondition3D2N::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)

	{
		KRATOS_TRY

		if(rLeftHandSideMatrix.size1() != 6)
			rLeftHandSideMatrix.resize(6,6,false);
		noalias(rLeftHandSideMatrix) = ZeroMatrix(6,6);
		if(rRightHandSideVector.size() != 6)
			rRightHandSideVector.resize(6,false);
		rRightHandSideVector = ZeroVector(6);

		KRATOS_CATCH("")
	}


	//************************************************************************************
	//************************************************************************************
	void TangentVelocityPeriodicCondition3D2N::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
	{
		 if (rResult.size()!=(6))
			    rResult.resize(6);
	    
	     //UPPER - LOWER PAIR
             //"paired" nodes of upper and lower walls: they have same x and same z... => tangentials are: x component and z component 
	     if(fabs(GetGeometry()[0].X()-GetGeometry()[1].X())<0.00001 && fabs(GetGeometry()[0].Z()-GetGeometry()[1].Z())<0.00001)
		 {
		 //X
		 //first node of the pair
		 rResult[0]=(GetGeometry()[0].GetDof(VELOCITY_X).EquationId());
		 //second node of teh pair
		 rResult[2]=(GetGeometry()[1].GetDof(VELOCITY_X).EquationId());
                 //Z
		 //first node of the pair
		 rResult[1]=(GetGeometry()[0].GetDof(VELOCITY_Z).EquationId());
		 //second node of teh pair
		 rResult[3]=(GetGeometry()[1].GetDof(VELOCITY_Z).EquationId());
		 //Lag multiplier dof
		 rResult[4] = GetGeometry()[0].GetDof(LAGRANGE_MULTIPLIER_TANGENT_VELOCITY).EquationId();
		 rResult[5] = GetGeometry()[0].GetDof(LAGRANGE_MULTIPLIER_SECOND_TANGENT_VELOCITY).EquationId();
		 }

	     //LEFT - RIGHT PAIR
             //"paired" nodes of left and right walls: they have same y and same z... => tangentials are: y component  and z component 
	    if(fabs(GetGeometry()[0].Y()-GetGeometry()[1].Y())<0.00001 && fabs(GetGeometry()[0].Z()-GetGeometry()[1].Z())<0.00001)
		 {
		 //Y
		 //first node of the pair
		 rResult[0]=(GetGeometry()[0].GetDof(VELOCITY_Y).EquationId());
		 //second node of teh pair
		 rResult[2]=(GetGeometry()[1].GetDof(VELOCITY_Y).EquationId());
                 //Z
		 //first node of the pair
		 rResult[1]=(GetGeometry()[0].GetDof(VELOCITY_Z).EquationId());
		 //second node of teh pair
		 rResult[3]=(GetGeometry()[1].GetDof(VELOCITY_Z).EquationId());
		 //Lag multiplier dof
		 rResult[4] = GetGeometry()[0].GetDof(LAGRANGE_MULTIPLIER_TANGENT_VELOCITY).EquationId();
		 rResult[5] = GetGeometry()[0].GetDof(LAGRANGE_MULTIPLIER_SECOND_TANGENT_VELOCITY).EquationId();
		 }


	     //BACK - FRONT PAIR
             //"paired" nodes of back and front walls: they have same x and same y... => tangentials are: x component  and y component 
	    if(fabs(GetGeometry()[0].X()-GetGeometry()[1].X())<0.00001 && fabs(GetGeometry()[0].Y()-GetGeometry()[1].Y())<0.00001)
		 {
		 //X
		 //first node of the pair
		 rResult[0]=(GetGeometry()[0].GetDof(VELOCITY_X).EquationId());
		 //second node of teh pair
		 rResult[2]=(GetGeometry()[1].GetDof(VELOCITY_X).EquationId());
                 //Y
		 //first node of the pair
		 rResult[1]=(GetGeometry()[0].GetDof(VELOCITY_Y).EquationId());
		 //second node of teh pair
		 rResult[3]=(GetGeometry()[1].GetDof(VELOCITY_Y).EquationId());
		 //Lag multiplier dof
		 rResult[4] = GetGeometry()[0].GetDof(LAGRANGE_MULTIPLIER_TANGENT_VELOCITY).EquationId();
		 rResult[5] = GetGeometry()[0].GetDof(LAGRANGE_MULTIPLIER_SECOND_TANGENT_VELOCITY).EquationId();
		 }
	}

	//************************************************************************************
	//************************************************************************************
	  void TangentVelocityPeriodicCondition3D2N::GetDofList(DofsVectorType& ConditionalDofList,ProcessInfo& rCurrentProcessInfo)
	{

         if (ConditionalDofList.size()!=6)
			  ConditionalDofList.resize(6);
//UPPER - LOWER PAIR
             //"paired" nodes of upper and lower walls: they have same x and same z... => tangentials are: x component and z component 
	     if(fabs(GetGeometry()[0].X()-GetGeometry()[1].X())<0.00001 && fabs(GetGeometry()[0].Z()-GetGeometry()[1].Z())<0.00001)
		 {
		 //X
		 //first node of the pair
		 ConditionalDofList[0] = GetGeometry()[0].pGetDof(VELOCITY_X);
		 //second node of teh pair
  	         ConditionalDofList[2] = GetGeometry()[1].pGetDof(VELOCITY_X);
                 //Z
		 //first node of the pair
		 ConditionalDofList[1] = GetGeometry()[0].pGetDof(VELOCITY_Z);
		 //second node of teh pair
		 ConditionalDofList[3] = GetGeometry()[1].pGetDof(VELOCITY_Z);
		 //Lag multiplier dof
		 ConditionalDofList[4] = GetGeometry()[0].pGetDof(LAGRANGE_MULTIPLIER_TANGENT_VELOCITY);
		 ConditionalDofList[5] = GetGeometry()[0].pGetDof(LAGRANGE_MULTIPLIER_SECOND_TANGENT_VELOCITY);
		 }

	     //LEFT - RIGHT PAIR
             //"paired" nodes of left and right walls: they have same y and same z... => tangentials are: y component  and z component 
	    if(fabs(GetGeometry()[0].Y()-GetGeometry()[1].Y())<0.00001 && fabs(GetGeometry()[0].Z()-GetGeometry()[1].Z())<0.00001)
		 {
		 //Y
		 //first node of the pair
		 ConditionalDofList[0] = GetGeometry()[0].pGetDof(VELOCITY_Y);		
		 //second node of teh pair
		 ConditionalDofList[2] = GetGeometry()[1].pGetDof(VELOCITY_Y);
                 //Z
		 //first node of the pair
		 ConditionalDofList[1] = GetGeometry()[0].pGetDof(VELOCITY_Z);		
		 //second node of teh pair
		 ConditionalDofList[3] = GetGeometry()[1].pGetDof(VELOCITY_Z);		
		 //Lag multiplier dof
		 ConditionalDofList[4] = GetGeometry()[0].pGetDof(LAGRANGE_MULTIPLIER_TANGENT_VELOCITY);
		 ConditionalDofList[5] = GetGeometry()[0].pGetDof(LAGRANGE_MULTIPLIER_SECOND_TANGENT_VELOCITY);
		 }


	     //BACK - FRONT PAIR
             //"paired" nodes of back and front walls: they have same x and same y... => tangentials are: x component  and y component 
	    if(fabs(GetGeometry()[0].X()-GetGeometry()[1].X())<0.00001 && fabs(GetGeometry()[0].Y()-GetGeometry()[1].Y())<0.00001)
		 {
		  //X
		 //first node of the pair
		 ConditionalDofList[0] = GetGeometry()[0].pGetDof(VELOCITY_X);
		 //second node of teh pair
  	         ConditionalDofList[2] = GetGeometry()[1].pGetDof(VELOCITY_X);
		 //Y
		 //first node of the pair
		 ConditionalDofList[1] = GetGeometry()[0].pGetDof(VELOCITY_Y);
		 //second node of teh pair
		 ConditionalDofList[3] = GetGeometry()[1].pGetDof(VELOCITY_Y);
		 //Lag multiplier dof
		 ConditionalDofList[4] = GetGeometry()[0].pGetDof(LAGRANGE_MULTIPLIER_TANGENT_VELOCITY);
		 ConditionalDofList[5] = GetGeometry()[0].pGetDof(LAGRANGE_MULTIPLIER_SECOND_TANGENT_VELOCITY);
		 }
	     
	}
} // Namespace Kratos
