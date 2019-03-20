// Project includes 
#include "includes/define.h"
#include "custom_conditions/tangent_velocity_periodic_normal_to_edge_condition_3D.h"
//#include "fluid_rve_lagrange_multipliers_application.h"
#include "ULF_application.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h" 

//author: Pavel Ryzhakov

namespace Kratos
{
	//************************************************************************************
	//************************************************************************************
	TangentVelocityPeriodicNormalToEdgeCondition3D2N::TangentVelocityPeriodicNormalToEdgeCondition3D2N(IndexType NewId, GeometryType::Pointer pGeometry)
		: Condition(NewId, pGeometry)
	{		
		//DO NOT ADD DOFS HERE!!!
	}

	//************************************************************************************
	//************************************************************************************
	TangentVelocityPeriodicNormalToEdgeCondition3D2N::TangentVelocityPeriodicNormalToEdgeCondition3D2N(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
		: Condition(NewId, pGeometry, pProperties)
	{
	}
	Condition::Pointer TangentVelocityPeriodicNormalToEdgeCondition3D2N::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		return Condition::Pointer(new TangentVelocityPeriodicNormalToEdgeCondition3D2N(NewId, GetGeometry().Create(ThisNodes), pProperties));
	}

	TangentVelocityPeriodicNormalToEdgeCondition3D2N::~TangentVelocityPeriodicNormalToEdgeCondition3D2N()
	{
	}

	//************************************************************************************
	//************************************************************************************
	void TangentVelocityPeriodicNormalToEdgeCondition3D2N::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY

		KRATOS_THROW_ERROR(std::logic_error,  "method not implemented" , "");
		
		KRATOS_CATCH("")
	}

	//************************************************************************************
	//************************************************************************************	
	void TangentVelocityPeriodicNormalToEdgeCondition3D2N::CalculateLocalVelocityContribution(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)

	{
		KRATOS_TRY
		//we only have to calculate if 
		if(rLeftHandSideMatrix.size1() != 3) // v1* v2* v1** v2** and lambda (these can be: v1x v2x, v1y v2y and lambda or v1x v2x v1z v2z and lambda etc... two tangential components per wall
			rLeftHandSideMatrix.resize(3,3,false);
		noalias(rLeftHandSideMatrix) = ZeroMatrix(3,3);
		if(rRightHandSideVector.size() != 3)
			rRightHandSideVector.resize(3,false);
		rRightHandSideVector = ZeroVector(3);
		
		boost::numeric::ublas::bounded_matrix<double, 3, 3 > LM_matrix = ZeroMatrix(3, 3) ; //(gradient)
		Geometry<Node<3> >& geom = this->GetGeometry();

		//READING THE G JUMP VECTOR: CONTAINS Gxy, Gxz, Gyz // IT IS SAVED IN THE AUX_VECTOR VARIABLE
		array_1d<double,3>  jump_vector= rCurrentProcessInfo[AUX_VECTOR];  //ZeroVector(3);


		//KRATOS_WATCH("CONDITION===========================================================================================")
		//LM is like this
		// 0   0   1   
		// 0   0  -1 
		// 1  -1   0 
		
		//to be multiplied by v1* v2* lambda

		LM_matrix(0,2)=1.0; 		LM_matrix(1,2)=-1.0;  		LM_matrix(2,0)=1.0;    		LM_matrix(2,1)=-1.0;

		
		//because we multiply it afterwards by vi vj lambda, where vi and vj are the velocities of the "paired" nodes i and j
		rLeftHandSideMatrix=LM_matrix;

		//WE ARE FIRST TRYING THE PROBLEMS WHERE PERIODICITY WITH JUMPS IS PRESENT IN LEFT-RIGHT and UP-LOW wallss in X and Y components only. 
		const array_1d<double,3> & velocity0 = GetGeometry()[0].GetSolutionStepValue(VELOCITY,0);
		const array_1d<double,3> & velocity1 = GetGeometry()[1].GetSolutionStepValue(VELOCITY,0);
		//substracting:
		array_1d<double,3>  temp_vector=ZeroVector(3);

		//LEFT AND RIGHT WALLS: Z component
		if(fabs(GetGeometry()[0].Y()-GetGeometry()[1].Y())<0.00001 && fabs(GetGeometry()[0].Z()-GetGeometry()[1].Z())<0.00001)
			{			
			//y_component
			temp_vector[0]=velocity0[2];
			temp_vector[1]=velocity1[2];
			//rRightHandSideVector[2] = 0.0;		//v1y*-v2y*=0
			rRightHandSideVector[2] = -jump_vector[2];		//v1y*-v2y*=0
			 //TO FIND OUT WHICH OF THE TWO PAIRS THESE ARE: VERTICAL OR HORIZONTAL (in z)... 
			 if (GetGeometry()[0].Z()>0.000001 && GetGeometry()[0].Z()<0.9999999)				
				temp_vector[2]=GetGeometry()[1].GetSolutionStepValue(LAGRANGE_MULTIPLIER_SECOND_TANGENT_VELOCITY,0); //store this DOF on the right wall nodes
			 //note that due to periodicity of the normal comoments, the nodes on the left and right are already "joined".. Thus we do it only on one side
			 else if (GetGeometry()[0].Z()<0.000001)
				{
				temp_vector[2]=GetGeometry()[0].GetSolutionStepValue(LAGRANGE_MULTIPLIER_SECOND_TANGENT_VELOCITY,0);
				}
			}
		//UP DOWN WALLS... Z component
		if(fabs(GetGeometry()[0].X()-GetGeometry()[1].X())<0.00001 && fabs(GetGeometry()[0].Z()-GetGeometry()[1].Z())<0.00001)
			{			
			//y_component
			temp_vector[0]=velocity0[2];
			temp_vector[1]=velocity1[2];
			//rRightHandSideVector[2] = 0.0;		//v1y*-v2y*=0
			rRightHandSideVector[2] = -jump_vector[2];
			 //TO FIND OUT WHICH OF THE TWO PAIRS THESE ARE: VERTICAL OR HORIZONTAL (in z)... 
			 if (GetGeometry()[0].Z()>0.000001 && GetGeometry()[0].Z()<0.9999999 && GetGeometry()[0].X()<0.000001)				
				temp_vector[2]=GetGeometry()[1].GetSolutionStepValue(LAGRANGE_MULTIPLIER_TANGENT_VELOCITY,0); //store this DOF on the up wall nodes
			 //note that due to periodicity of the normal comoments, the nodes on the left and right are already "joined".. Thus we do it only on one side
			 else if (GetGeometry()[0].Z()>0.99999) 
				{
				temp_vector[2]=GetGeometry()[0].GetSolutionStepValue(LAGRANGE_MULTIPLIER_SECOND_TANGENT_VELOCITY,0);
				}
			}
		
		//FRONT BACK WALLS X component
		if(fabs(GetGeometry()[0].X()-GetGeometry()[1].X())<0.00001 && fabs(GetGeometry()[0].Y()-GetGeometry()[1].Y())<0.00001)
			{
			//KRATOS_WATCH("Im hereeeeeeeeeee")
				
			if (GetGeometry()[0].X()>0.000001 && GetGeometry()[0].X()<0.9999999  && GetGeometry()[0].Y()>0.99999)
				{
				//x_component
				temp_vector[0]=velocity0[0];
				temp_vector[1]=velocity1[0];			
				//primera componente tangente
				//NO JUMP IN THIS COMPONENT BETWEEN BACK AND FRONT WALLS
				rRightHandSideVector[2] = -jump_vector[2];		//v1x*-v2x*=0	
	
				temp_vector[2]=GetGeometry()[0].GetSolutionStepValue(LAGRANGE_MULTIPLIER_SECOND_TANGENT_VELOCITY,0);	
				
				}

			else if (GetGeometry()[0].Y()>0.000001 && GetGeometry()[0].Y()<0.9999999 && GetGeometry()[0].X()>0.99999)//back left vertical edge
				{
				//y_component
				temp_vector[0]=velocity0[0];
				temp_vector[1]=velocity1[0];			
				//primera componente tangente
				//NO JUMP IN THIS COMPONENT BETWEEN BACK AND FRONT WALLS
				rRightHandSideVector[2] = -jump_vector[2];		//v1x*-v2x*=0		
				temp_vector[2]=GetGeometry()[0].GetSolutionStepValue(LAGRANGE_MULTIPLIER_SECOND_TANGENT_VELOCITY,0);	
				//KRATOS_WATCH(GetGeometry()[0])
				//KRATOS_WATCH(GetGeometry()[1])
				}

			}

		
		
		
		
 		noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,temp_vector);
			



		KRATOS_CATCH("")
	}
	

	void TangentVelocityPeriodicNormalToEdgeCondition3D2N::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)

	{
		KRATOS_TRY

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
	void TangentVelocityPeriodicNormalToEdgeCondition3D2N::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
	{
		 if (rResult.size()!=(3))
			    rResult.resize(3);	
		//LEFT and RIGHT    
	     if(fabs(GetGeometry()[0].Y()-GetGeometry()[1].Y())<0.00001 && fabs(GetGeometry()[0].Z()-GetGeometry()[1].Z())<0.00001)
		 {
		 //Y
		 //first node of the pair
		 rResult[0]=(GetGeometry()[0].GetDof(VELOCITY_Z).EquationId());
		 //second node of teh pair
		 rResult[1]=(GetGeometry()[1].GetDof(VELOCITY_Z).EquationId());
		 if (GetGeometry()[0].Z()>0.000001 && GetGeometry()[0].Z()<0.9999999)		
			 rResult[2] = GetGeometry()[1].GetDof(LAGRANGE_MULTIPLIER_SECOND_TANGENT_VELOCITY).EquationId();//we put the DOF on the right wall not to coincide with the dof of vel x
		 else if (GetGeometry()[0].Z()<0.00001)
			 rResult[2] = GetGeometry()[0].GetDof(LAGRANGE_MULTIPLIER_SECOND_TANGENT_VELOCITY).EquationId();
		 }
	

		//UP DOWN WALLS... Z component
		if(fabs(GetGeometry()[0].X()-GetGeometry()[1].X())<0.00001 && fabs(GetGeometry()[0].Z()-GetGeometry()[1].Z())<0.00001)
			{			
			rResult[0]=(GetGeometry()[0].GetDof(VELOCITY_Z).EquationId());
		 	//second node of teh pair
			rResult[1]=(GetGeometry()[1].GetDof(VELOCITY_Z).EquationId());
			 //TO FIND OUT WHICH OF THE TWO PAIRS THESE ARE: VERTICAL OR HORIZONTAL (in z)... 
			 if (GetGeometry()[0].Z()>0.000001 && GetGeometry()[0].Z()<0.9999999 && GetGeometry()[0].X()<0.000001)			
				 rResult[2] = GetGeometry()[1].GetDof(LAGRANGE_MULTIPLIER_TANGENT_VELOCITY).EquationId(); //store this DOF on the up wall nodes
			 //note that due to periodicity of the normal comoments, the nodes on the left and right are already "joined".. Thus we do it only on one side
			 else if (GetGeometry()[0].Z()>0.99999) 
				{
				rResult[2] = GetGeometry()[0].GetDof(LAGRANGE_MULTIPLIER_SECOND_TANGENT_VELOCITY).EquationId();
				}
			}



	     //BACK - FRONT PAIR
             //"paired" nodes of back and front walls: they have same x and same y... => tangentials are: x component  and y component 
	       if(fabs(GetGeometry()[0].X()-GetGeometry()[1].X())<0.00001 && fabs(GetGeometry()[0].Y()-GetGeometry()[1].Y())<0.00001)
		 {
	 if (GetGeometry()[0].X()>0.000001 && GetGeometry()[0].X()<0.9999999 &&  GetGeometry()[0].Y()>0.99999)
			{
			 //X
			 //first node of the pair
			 rResult[0]=(GetGeometry()[0].GetDof(VELOCITY_X).EquationId());
			 //second node of teh pair
			 rResult[1]=(GetGeometry()[1].GetDof(VELOCITY_X).EquationId());
			 rResult[2] = GetGeometry()[0].GetDof(LAGRANGE_MULTIPLIER_SECOND_TANGENT_VELOCITY).EquationId();
			}

		else if (GetGeometry()[0].Y()>0.000001 && GetGeometry()[0].Y()<0.9999999 && GetGeometry()[0].X()>0.999999) //back left vertical edge
			{
			//X
			 //first node of the pair
			 rResult[0]=(GetGeometry()[0].GetDof(VELOCITY_X).EquationId());
			 //second node of teh pair
			 rResult[1]=(GetGeometry()[1].GetDof(VELOCITY_X).EquationId());
			 rResult[2] = GetGeometry()[0].GetDof(LAGRANGE_MULTIPLIER_SECOND_TANGENT_VELOCITY).EquationId();
			}

		 }
	}

	//************************************************************************************
	//************************************************************************************
	  void TangentVelocityPeriodicNormalToEdgeCondition3D2N::GetDofList(DofsVectorType& ConditionalDofList,ProcessInfo& rCurrentProcessInfo)
	{

         if (ConditionalDofList.size()!=3)
			  ConditionalDofList.resize(3);
 //LEFT - RIGHT PAIR
             //"paired" nodes of left and right walls: they have same y and same z... => tangentials are: y component  and z component 
	    if(fabs(GetGeometry()[0].Y()-GetGeometry()[1].Y())<0.00001 && fabs(GetGeometry()[0].Z()-GetGeometry()[1].Z())<0.00001)
		 {
		 //Y
		 //first node of the pair
		 ConditionalDofList[0] = GetGeometry()[0].pGetDof(VELOCITY_Z);
		 //second node of teh pair
  	         ConditionalDofList[1] = GetGeometry()[1].pGetDof(VELOCITY_Z);
		//we store it here on the right side, because the lower edge of the left side is already used for secon multiplier of X velocity
 		 if (GetGeometry()[0].Z()>0.000001 && GetGeometry()[0].Z()<0.9999999)		
		 	ConditionalDofList[2] = GetGeometry()[1].pGetDof(LAGRANGE_MULTIPLIER_SECOND_TANGENT_VELOCITY);//we put the DOF on the right wall not to coincide with the dof of vel x
		 else if (GetGeometry()[0].Z()<0.00001)//NOTE THAT SECOND_TANGENT variable in these nodes is already used for the previous condition! Thats why here TANGENT and not SECOND_TANGENT	
			ConditionalDofList[2] = GetGeometry()[0].pGetDof(LAGRANGE_MULTIPLIER_SECOND_TANGENT_VELOCITY); //these are the left wall edges oriented in z-direction
		 }


//UP DOWN WALLS... Z component
		if(fabs(GetGeometry()[0].X()-GetGeometry()[1].X())<0.00001 && fabs(GetGeometry()[0].Z()-GetGeometry()[1].Z())<0.00001)
			{			
			ConditionalDofList[0] = GetGeometry()[0].pGetDof(VELOCITY_Z);
		 	//second node of teh pair
			 ConditionalDofList[1] = GetGeometry()[1].pGetDof(VELOCITY_Z);
			 //TO FIND OUT WHICH OF THE TWO PAIRS THESE ARE: VERTICAL OR HORIZONTAL (in z)... 
			 if (GetGeometry()[0].Z()>0.000001 && GetGeometry()[0].Z()<0.9999999 && GetGeometry()[0].X()<0.000001)		
				 ConditionalDofList[2] = GetGeometry()[1].pGetDof(LAGRANGE_MULTIPLIER_TANGENT_VELOCITY);//store this DOF on the up wall nodes
			 //note that due to periodicity of the normal comoments, the nodes on the left and right are already "joined".. Thus we do it only on one side
			 else if (GetGeometry()[0].Z()>0.99999) 
				{
				ConditionalDofList[2] = GetGeometry()[0].pGetDof(LAGRANGE_MULTIPLIER_SECOND_TANGENT_VELOCITY); 
				}
			}


	     //BACK - FRONT PAIR
             //"paired" nodes of back and front walls: they have same x and same y... => tangentials are: x component  and y component 
	  if(fabs(GetGeometry()[0].X()-GetGeometry()[1].X())<0.00001 && fabs(GetGeometry()[0].Y()-GetGeometry()[1].Y())<0.00001)
		 {
		  //X
	 	if (GetGeometry()[0].X()>0.000001 && GetGeometry()[0].X()<0.9999999 && GetGeometry()[0].Y()>0.99999)
			{
			 //first node of the pair
			 ConditionalDofList[0] = GetGeometry()[0].pGetDof(VELOCITY_X);
			 //second node of teh pair
	  	         ConditionalDofList[1] = GetGeometry()[1].pGetDof(VELOCITY_X);
			 ConditionalDofList[2] = GetGeometry()[0].pGetDof(LAGRANGE_MULTIPLIER_SECOND_TANGENT_VELOCITY);
			}

	else if (GetGeometry()[0].Y()>0.000001 && GetGeometry()[0].Y()<0.9999999 && GetGeometry()[0].X()>0.999999) //back right vertical edge
			{
			 //first node of the pair
			 ConditionalDofList[0] = GetGeometry()[0].pGetDof(VELOCITY_X);
			 //second node of teh pair
	  	         ConditionalDofList[1] = GetGeometry()[1].pGetDof(VELOCITY_X);
			 ConditionalDofList[2] = GetGeometry()[0].pGetDof(LAGRANGE_MULTIPLIER_SECOND_TANGENT_VELOCITY);
			}


		 }
	     
	}
} // Namespace Kratos
