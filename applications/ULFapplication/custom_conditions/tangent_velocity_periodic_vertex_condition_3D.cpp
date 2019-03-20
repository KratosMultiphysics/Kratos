// Project includes 
#include "includes/define.h"
#include "custom_conditions/tangent_velocity_periodic_vertex_condition_3D.h"
//#include "fluid_rve_lagrange_multipliers_application.h"
#include "ULF_application.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h" 

//author: Pavel Ryzhakov

namespace Kratos
{
	//************************************************************************************
	//************************************************************************************
	TangentVelocityPeriodicVertexCondition3D2N::TangentVelocityPeriodicVertexCondition3D2N(IndexType NewId, GeometryType::Pointer pGeometry)
		: Condition(NewId, pGeometry)
	{		
		//DO NOT ADD DOFS HERE!!!
	}

	//************************************************************************************
	//************************************************************************************
	TangentVelocityPeriodicVertexCondition3D2N::TangentVelocityPeriodicVertexCondition3D2N(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
		: Condition(NewId, pGeometry, pProperties)
	{
	}
	Condition::Pointer TangentVelocityPeriodicVertexCondition3D2N::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		return Condition::Pointer(new TangentVelocityPeriodicVertexCondition3D2N(NewId, GetGeometry().Create(ThisNodes), pProperties));
	}

	TangentVelocityPeriodicVertexCondition3D2N::~TangentVelocityPeriodicVertexCondition3D2N()
	{
	}

	//************************************************************************************
	//************************************************************************************
	void TangentVelocityPeriodicVertexCondition3D2N::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY

		KRATOS_THROW_ERROR(std::logic_error,  "method not implemented" , "");
		
		KRATOS_CATCH("")
	}

	//************************************************************************************
	//************************************************************************************	
	void TangentVelocityPeriodicVertexCondition3D2N::CalculateLocalVelocityContribution(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)

	{
		KRATOS_TRY

		//LM is like this
		// 0   0   1   
		// 0   0  -1 
		// 1  -1   0 
		//to be multiplied by v1* v2* lambda

		//because we multiply it afterwards by vi vj lambda, where vi and vj are the velocities of the "paired" nodes i and j
		if(rLeftHandSideMatrix.size1() != 3) // v1* v2* v1** v2** and lambda (these can be: v1x v2x, v1y v2y and lambda or v1x v2x v1z v2z and lambda etc... two tang compon per wall	
			rLeftHandSideMatrix.resize(3,3,false);
		noalias(rLeftHandSideMatrix) = ZeroMatrix(3,3);
		if(rRightHandSideVector.size() != 3)
			rRightHandSideVector.resize(3,false);
		rRightHandSideVector = ZeroVector(3);
		array_1d<double,3>  temp_vector=ZeroVector(3);
			
		boost::numeric::ublas::bounded_matrix<double, 3, 3 > LM_matrix = ZeroMatrix(3, 3) ; //(gradient)
		Geometry<Node<3> >& geom = this->GetGeometry();

		//READING THE G JUMP VECTOR: CONTAINS Gxy, Gxz, Gyz // IT IS SAVED IN THE AUX_VECTOR VARIABLE
		array_1d<double,3>  jump_vector= rCurrentProcessInfo[AUX_VECTOR];  //ZeroVector(3);

	
		LM_matrix(0,2)=1.0; 		LM_matrix(1,2)=-1.0;  		LM_matrix(2,0)=1.0;    		LM_matrix(2,1)=-1.0;
			
		rLeftHandSideMatrix=LM_matrix;
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	

		//WE ARE FIRST TRYING THE PROBLEMS WHERE PERIODICITY WITH JUMPS IS PRESENT IN LEFT-RIGHT and UP-LOW wallss in X and Y components only. 
		const array_1d<double,3> & velocity0 = GetGeometry()[0].GetSolutionStepValue(VELOCITY,0);
		const array_1d<double,3> & velocity1 = GetGeometry()[1].GetSolutionStepValue(VELOCITY,0);
		//substracting:

		
		//"paired" nodes of up and down walls. We use right side of the walls.. and store teh dof in the upper right edge
		if(fabs(GetGeometry()[0].X()-GetGeometry()[1].X())<0.00001 && fabs(GetGeometry()[0].Z()-GetGeometry()[1].Z())<0.00001)
			{
			//x_component
			KRATOS_WATCH("Here X-vel")
			temp_vector[0]=velocity0[0];
			temp_vector[1]=velocity1[0];			
			rRightHandSideVector[2] = -jump_vector[0];		//v1x*-v2x*=0
			//rRightHandSideVector[2] = -jump_vector[0]+jump_vector[0]*diminish_jump_value;
			 //TO FIND OUT WHICH OF THE TWO PAIRS THESE ARE: VERTICAL OR HORIZONTAL... 
			temp_vector[2]=GetGeometry()[1].GetSolutionStepValue(LAGRANGE_MULTIPLIER_TANGENT_VELOCITY,0);	//we store it in the second node
			}



		//"paired" nodes of left and right walls. We use the left wall, down edge vertices
		if(fabs(GetGeometry()[0].Y()-GetGeometry()[1].Y())<0.00001 && fabs(GetGeometry()[0].Z()-GetGeometry()[1].Z())<0.00001)
			{			
			//y_component
			KRATOS_WATCH("Here Y-vel")
			temp_vector[0]=velocity0[1];
			temp_vector[1]=velocity1[1];
			rRightHandSideVector[2] = -jump_vector[0];		//v1y*-v2y*=0	
			//rRightHandSideVector[2] = -jump_vector[0]+jump_vector[0]*diminish_jump_value;
			temp_vector[2]=GetGeometry()[0].GetSolutionStepValue(LAGRANGE_MULTIPLIER_TANGENT_VELOCITY,0);								
			}

		//"paired" nodes of front and back walls . We store the LAG MULT dof in the second node
		if(fabs(GetGeometry()[0].X()-GetGeometry()[1].X())<0.00001 && fabs(GetGeometry()[0].Y()-GetGeometry()[1].Y())<0.00001)
			{
			KRATOS_WATCH("Here Y-vel")
			temp_vector[0]=velocity0[1];
			temp_vector[1]=velocity1[1];			
			//NO JUMP IN THIS COMPONENT BETWEEN BACK AND FRONT WALLS
			rRightHandSideVector[2] = -jump_vector[2];	//v1x*-v2x*=0		
			temp_vector[2]=GetGeometry()[1].GetSolutionStepValue(LAGRANGE_MULTIPLIER_TANGENT_VELOCITY,0);					
			}
		
 		noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,temp_vector);
			



		KRATOS_CATCH("")
	}
	

	void TangentVelocityPeriodicVertexCondition3D2N::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)

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
	void TangentVelocityPeriodicVertexCondition3D2N::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
	{
		 if (rResult.size()!=(3))
			    rResult.resize(3);
//KRATOS_WATCH("Im hereeeeeeeeeee22222222222")	    

	     //UPPER - LOWER PAIR
             //"paired" nodes of upper and lower walls: they have same x and same z... => tangentials are: x component and z component 
	     if(fabs(GetGeometry()[0].X()-GetGeometry()[1].X())<0.00001 && fabs(GetGeometry()[0].Z()-GetGeometry()[1].Z())<0.00001)
		 {
		 //X
		 //first node of the pair
		 rResult[0]=(GetGeometry()[0].GetDof(VELOCITY_X).EquationId());
		 //second node of teh pair
		 rResult[1]=(GetGeometry()[1].GetDof(VELOCITY_X).EquationId());
		 //TO FIND OUT WHICH OF THE TWO PAIRS THESE ARE: VERTICAL OR HORIZONTAL... 		 
		 rResult[2] = GetGeometry()[1].GetDof(LAGRANGE_MULTIPLIER_TANGENT_VELOCITY).EquationId(); //store LAG MULT DOF in in the second node
		 }


	     //LEFT - RIGHT PAIR
             //"paired" nodes of left and right walls: they have same y and same z... => tangentials are: y component  and z component 
	    if(fabs(GetGeometry()[0].Y()-GetGeometry()[1].Y())<0.00001 && fabs(GetGeometry()[0].Z()-GetGeometry()[1].Z())<0.00001)
		 {
		 //Y
		 //first node of the pair
		 rResult[0]=(GetGeometry()[0].GetDof(VELOCITY_Y).EquationId());
		 //second node of teh pair
		 rResult[1]=(GetGeometry()[1].GetDof(VELOCITY_Y).EquationId());
		 rResult[2] = GetGeometry()[0].GetDof(LAGRANGE_MULTIPLIER_TANGENT_VELOCITY).EquationId();		 
		 }
	

	     //BACK - FRONT PAIR
             //"paired" nodes of back and front walls: they have same x and same y... => tangentials are: x component  and y component 
	    if(fabs(GetGeometry()[0].X()-GetGeometry()[1].X())<0.00001 && fabs(GetGeometry()[0].Y()-GetGeometry()[1].Y())<0.00001)
		 {
		 //Y
		 //first node of the pair
		 rResult[0]=(GetGeometry()[0].GetDof(VELOCITY_Y).EquationId());
		 //second node of teh pair
		 rResult[1]=(GetGeometry()[1].GetDof(VELOCITY_Y).EquationId());
		 rResult[2] = GetGeometry()[1].GetDof(LAGRANGE_MULTIPLIER_TANGENT_VELOCITY).EquationId();		//we store DOF in the second node
		 }
	}

	//************************************************************************************
	//************************************************************************************
	  void TangentVelocityPeriodicVertexCondition3D2N::GetDofList(DofsVectorType& ConditionalDofList,ProcessInfo& rCurrentProcessInfo)
	{

         if (ConditionalDofList.size()!=3)
			  ConditionalDofList.resize(3);

//KRATOS_WATCH("Im hereeeeeeeeeee111111")
//UPPER - LOWER PAIR
             //"paired" nodes of upper and lower walls: they have same x and same z... => tangentials are: x component and z component 
	     if(fabs(GetGeometry()[0].X()-GetGeometry()[1].X())<0.00001 && fabs(GetGeometry()[0].Z()-GetGeometry()[1].Z())<0.00001)
		 {
		 //X
		 //first node of the pair
		 ConditionalDofList[0] = GetGeometry()[0].pGetDof(VELOCITY_X);
		 //second node of teh pair
  	         ConditionalDofList[1] = GetGeometry()[1].pGetDof(VELOCITY_X);
		 //If this is a uppper-lower horizontal edges
	 	ConditionalDofList[2] = GetGeometry()[1].pGetDof(LAGRANGE_MULTIPLIER_TANGENT_VELOCITY); //we store LAG dof in the second node
		 }

	     //LEFT - RIGHT PAIR
             //"paired" nodes of left and right walls: they have same y and same z... => tangentials are: y component  and z component 
	    if(fabs(GetGeometry()[0].Y()-GetGeometry()[1].Y())<0.00001 && fabs(GetGeometry()[0].Z()-GetGeometry()[1].Z())<0.00001)
		 {
		 //Y
		 //first node of the pair
		 ConditionalDofList[0] = GetGeometry()[0].pGetDof(VELOCITY_Y);
		 //second node of teh pair
  	         ConditionalDofList[1] = GetGeometry()[1].pGetDof(VELOCITY_Y);
 		 ConditionalDofList[2] = GetGeometry()[0].pGetDof(LAGRANGE_MULTIPLIER_TANGENT_VELOCITY);
		 }


	     //BACK - FRONT PAIR
             //"paired" nodes of back and front walls: they have same x and same y... => tangentials are: x component  and y component 
	    if(fabs(GetGeometry()[0].X()-GetGeometry()[1].X())<0.00001 && fabs(GetGeometry()[0].Y()-GetGeometry()[1].Y())<0.00001)
		 {
		  //Y
		 //first node of the pair
		 ConditionalDofList[0] = GetGeometry()[0].pGetDof(VELOCITY_Y);
		 //second node of teh pair
  	         ConditionalDofList[1] = GetGeometry()[1].pGetDof(VELOCITY_Y);
		 ConditionalDofList[2] = GetGeometry()[1].pGetDof(LAGRANGE_MULTIPLIER_TANGENT_VELOCITY);	//we store LAG dof in the second node


		 }
	     
	}
} // Namespace Kratos
