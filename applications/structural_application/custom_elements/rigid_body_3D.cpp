//   
//   Project Name:        Kratos       
//   Last modified by:    $Author: rrossi $
//   Date:                $Date: 2009-01-13 15:39:56 $
//   Revision:            $Revision: 1.8 $
//
//
 

// System includes 


// External includes 


// Project includes 
#include "includes/define.h"
#include "custom_elements/rigid_body_3D.h"
#include "structural_application.h"
#include "utilities/math_utils.h"

namespace Kratos
{

	//************************************************************************************
	//************************************************************************************
	RigidBody3D::RigidBody3D(IndexType NewId, GeometryType::Pointer pGeometry)
		: Element(NewId, pGeometry)
	{		
		//DO NOT ADD DOFS HERE!!!
	}

	//************************************************************************************
	//************************************************************************************

	RigidBody3D::RigidBody3D(IndexType NewId, 
				 GeometryType::Pointer pGeometry,  
     				 PropertiesType::Pointer pProperties, 
	  			 GeometryType::Pointer  rskin_nodes, 
       				 double& mass, 
	    			 Matrix& Inertia,
				 array_1d<double,3>& translational_stiffness,
				 array_1d<double,3>& rotational_stiffness
		)
		: Element(NewId, pGeometry, pProperties
)
	{
		mtranslational_stiffness = translational_stiffness;
		mrotational_stiffness = rotational_stiffness;

		m_skin_nodes = rskin_nodes;
		
		if( Inertia.size1() != 3)
			KRATOS_ERROR(std::logic_error,"wrong size of the INERTIA matrix .... should be 3*3","");
		mInertia.resize(3,3,false);
		noalias(mInertia) = Inertia;

		mmass = mass;
		if(mmass == 0.0)
			KRATOS_ERROR(std::logic_error,"the mass of the center node should be different from 0","");
			
		

	}

	Element::Pointer RigidBody3D::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
	  KRATOS_ERROR(std::logic_error,"this element does not allow the standard constuctor","");
	  //		return Element::Pointer(new RigidBody3D(NewId, GetGeometry().Create(ThisNodes), pProperties));
	}

	RigidBody3D::~RigidBody3D()
	{
	}

	//************************************************************************************
	//************************************************************************************
	void RigidBody3D::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY

		unsigned int number_of_points = 1;
		unsigned int ndofs = 6;

		//resizing as needed the LHS
		if(rLeftHandSideMatrix.size1() != number_of_points*ndofs)
			rLeftHandSideMatrix.resize(number_of_points*ndofs,number_of_points*ndofs,false);
		noalias(rLeftHandSideMatrix) = ZeroMatrix(number_of_points*ndofs,number_of_points*ndofs); //resetting LHS



		//resizing as needed the RHS
		if(rRightHandSideVector.size() != number_of_points*ndofs)
			rRightHandSideVector.resize(number_of_points*ndofs,false);
		rRightHandSideVector = ZeroVector(number_of_points*ndofs); //resetting RHS

		CalculateForces(rRightHandSideVector, rCurrentProcessInfo);

		//adding spring stiffnesses
		for(unsigned int i=0; i<3; i++)
		{
			rLeftHandSideMatrix(i,i) += mtranslational_stiffness[i];
			rLeftHandSideMatrix(3+i,3+i) += mrotational_stiffness[i];
		}

		//add contribution to the RHS
		Vector temp(6);
		GetValuesVector(temp,0); //gettign the current displacement
		for(unsigned int i=0; i<3; i++)
		{
			rRightHandSideVector[i] -= mtranslational_stiffness[i]*temp[i]; 
			rRightHandSideVector[3+i] -= mrotational_stiffness[i]*temp[3+i];		
			
		}
		//KRATOS_WATCH(mtranslational_stiffness);
		//KRATOS_WATCH(mrotational_stiffness);

		KRATOS_CATCH("")
	}

	//************************************************************************************
	//************************************************************************************
	void RigidBody3D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		//resizing as needed the RHS
		if(rRightHandSideVector.size() != 6)
			rRightHandSideVector.resize(6,false);
		rRightHandSideVector = ZeroVector(6); //resetting RHS

		CalculateForces(rRightHandSideVector, rCurrentProcessInfo);
	}


	//************************************************************************************
	//************************************************************************************
	//this subroutine calculates the forces on the center node.
	//the nodes on the external surface are in geometry at position 1..n
	//the center node at position 1
	void RigidBody3D::CalculateForces(VectorType& rExtForces, ProcessInfo& rCurrentProcessInfo)
	{
		unsigned int ndofs = 6;
		noalias(rExtForces) = ZeroVector(ndofs);

		array_1d<double,3> NodeForce, NodeMoment, TotForce, TotMoment, dist;

		double Xbase = GetGeometry()[0].X();
		double Ybase = GetGeometry()[0].Y();
		double Zbase = GetGeometry()[0].Z();

		noalias(TotForce) = ZeroVector(3);
		noalias(TotMoment) = ZeroVector(3);

		//loop on the nodes of the outer surface
		//note that everything is done as if it was in 3D
		//unsigned int origin = 1;
		
		for(GeometryType::iterator it = m_skin_nodes->begin();
				  it !=  m_skin_nodes->end(); it++)
		{
			noalias(NodeForce) = it->FastGetSolutionStepValue(FORCE);
		
			//calculating distance spring center - node
			dist[0] = it->X() - Xbase;
			dist[1] = it->Y() - Ybase;
			dist[2] = it->Z() - Zbase;

			//adding nodal contribution to the total force
			noalias(TotForce) += NodeForce;
		
			//adding contributions to moment
			MathUtils<double>::CrossProduct(NodeMoment,dist,NodeForce);
			noalias(TotMoment) += NodeMoment;
		
		}
		
		//adding the contribution of the weight
		noalias(TotForce) += mmass * GetGeometry()[0].FastGetSolutionStepValue(BODY_FORCE);

		//writing the contributions in the ext force
		rExtForces[0] = TotForce[0];
		rExtForces[1] = TotForce[1];
		rExtForces[2] = TotForce[2];
		rExtForces[3] = -TotMoment[0];
		rExtForces[4] = -TotMoment[1];
		rExtForces[5] = -TotMoment[2];

                KRATOS_WATCH(rExtForces)
KRATOS_WATCH(TotForce)
KRATOS_WATCH(TotMoment)
                noalias( GetGeometry()[0].GetValue(FORCE_CM) ) = TotForce;
                noalias( GetGeometry()[0].GetValue(MOMENTUM_CM) ) = -TotMoment;
	}

	//************************************************************************************
	//************************************************************************************
	//this subroutine is used to update the external shape once the central node is moved
	//it should be called just after the update phase
	void RigidBody3D::UpdateExtShape(ProcessInfo& rCurrentProcessInfo)
	{
		array_1d<double,3> NodeDisp, RotationalDisplacement, vel, dist;
		
		//reading displacement and rotation for the center node
		const array_1d<double,3>&  SpringDisplacement = GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT);
		const array_1d<double,3>  SpringRotation = GetGeometry()[0].FastGetSolutionStepValue(ROTATION) - 
							   GetGeometry()[0].GetValue(ROTATION);

		//calculating incremental rotation matrix
		double cx = cos(SpringRotation[0]);
		double sx = sin(SpringRotation[0]);
		double cy = cos(SpringRotation[1]);
		double sy = sin(SpringRotation[1]);
		double cz = cos(SpringRotation[2]);
		double sz = sin(SpringRotation[2]);
		boost::numeric::ublas::bounded_matrix<double,3,3> IncrementalRot;
		boost::numeric::ublas::bounded_matrix<double,3,3> Rot;
		
		IncrementalRot(0,0) = cz*cy+sz*sx*sy;	IncrementalRot(0,1) = sz*cx;	IncrementalRot(0,2) = -cz*sy+sz*sx*cy;
		IncrementalRot(1,0) = -sz*cy+cz*sx*sy;	IncrementalRot(1,1) = cz*cx;	IncrementalRot(1,2) = sz*sy+cz*sx*cy;
		IncrementalRot(2,0) = cx*sy;		IncrementalRot(2,1) = -sx;	IncrementalRot(2,2) = cx*cy;
		
		//calculating total rotation
		noalias(Rot) = prod(IncrementalRot,mRot);
		
/*KRATOS_WATCH(mRot);
KRATOS_WATCH(IncrementalRot);
KRATOS_WATCH(Rot);*/
		
		//saving the last rotation used in the database
		noalias(mRot) = Rot;
		noalias(GetGeometry()[0].GetValue(ROTATION)) = GetGeometry()[0].FastGetSolutionStepValue(ROTATION);
				
		double Xbase = GetGeometry()[0].X0();
		double Ybase = GetGeometry()[0].Y0();
		double Zbase = GetGeometry()[0].Z0();

		//loop on the nodes of the outer surface
		//note that everything is done as if it was in 3D
		for(GeometryType::iterator it = m_skin_nodes->begin();
				  it !=  m_skin_nodes->end(); it++)
		{
			//calculating distance spring center - node
			dist[0] = it->X0() - Xbase;
			dist[1] = it->Y0() - Ybase;
			dist[2] = it->Z0() - Zbase;

			//calculating the displacement due to the rotation
			noalias(RotationalDisplacement) = prod(Rot,dist);
			noalias(RotationalDisplacement) -= dist; //bug
// 			RotationalDisplacement[0] -= dist[0];
// 			RotationalDisplacement[1] -= dist[1];
// 			RotationalDisplacement[2] -= dist[2];

			//calculating the nodal total displacement
			noalias(NodeDisp) =	SpringDisplacement;
			noalias(NodeDisp) +=	RotationalDisplacement;

			//update position of the nodes 
			it->X() = it->X0() + NodeDisp[0];
			it->Y() = it->Y0() + NodeDisp[1];
			it->Z() = it->Z0() + NodeDisp[2];

			//update displacement value 
			it->FastGetSolutionStepValue(DISPLACEMENT) = NodeDisp;

			//calculate velocity of the nodes 
			//***********************************
			//ATTENTION!! like this it is not consistent!! put the consistent formulation!!
			//***********************************
			double DeltaTime = rCurrentProcessInfo[DELTA_TIME];
			noalias(vel) = NodeDisp;
			noalias(vel) -= it->FastGetSolutionStepValue(DISPLACEMENT,1);
			vel /= DeltaTime;
			it->FastGetSolutionStepValue(VELOCITY) = vel;

		}
		//rotating the inertia matrix
		boost::numeric::ublas::bounded_matrix<double,3,3> temp;
		noalias(temp) = prod( mInertia, trans(Rot) );
		noalias( mRotatedInertia ) = prod( Rot , temp);
	
// 		noalias(temp) = prod( mInertia, Rot );
// 		noalias( mRotatedInertia ) = prod( trans(Rot) , temp);

	}


	//************************************************************************************
	//************************************************************************************
	  void RigidBody3D::FinalizeNonLinearIteration(ProcessInfo& CurrentProcessInfo)
	  {UpdateExtShape(CurrentProcessInfo);}
	  void RigidBody3D::InitializeNonLinearIteration(ProcessInfo& CurrentProcessInfo)
	  {UpdateExtShape(CurrentProcessInfo);}
/*	  void RigidBody3D::FinalizeSolutionStep(ProcessInfo& CurrentProcessInfo)
		{UpdateExtShape(CurrentProcessInfo);}*/
	  
// 	  void RigidBody3D::InitializeNonLinearIteration(ProcessInfo& CurrentProcessInfo)
// 	  {		
// 		  Vector temp(6);
// 		  CalculateForces(temp, CurrentProcessInfo);
// 		  KRATOS_WATCH(temp);
// 	
// 	  }	  
// 	  
// 	  void RigidBody3D::InitializeSolutionStep(ProcessInfo& CurrentProcessInfo)
// 		{	
// 			
// 			Vector temp(6);
// 			CalculateForces(temp, CurrentProcessInfo);
// 			KRATOS_WATCH(temp);
// 	
// 		}

	//************************************************************************************
	//************************************************************************************
	void RigidBody3D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
	{
		unsigned int ndofs = 6;

		if(rResult.size() != ndofs)
			rResult.resize(ndofs,false);

		rResult[0] = GetGeometry()[0].GetDof(DISPLACEMENT_X).EquationId();
		rResult[1] = GetGeometry()[0].GetDof(DISPLACEMENT_Y).EquationId();
		rResult[2] = GetGeometry()[0].GetDof(DISPLACEMENT_Z).EquationId();
		rResult[3] = GetGeometry()[0].GetDof(ROTATION_X).EquationId();
		rResult[4] = GetGeometry()[0].GetDof(ROTATION_Y).EquationId();
		rResult[5] = GetGeometry()[0].GetDof(ROTATION_Z).EquationId();

	}

	//************************************************************************************
	//************************************************************************************
	  void RigidBody3D::GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
	{
		unsigned int ndofs = 6;

		if(ElementalDofList.size() != ndofs)
			ElementalDofList.resize(ndofs);

		ElementalDofList[0] = GetGeometry()[0].pGetDof(DISPLACEMENT_X);
		ElementalDofList[1] = GetGeometry()[0].pGetDof(DISPLACEMENT_Y);
		ElementalDofList[2] = GetGeometry()[0].pGetDof(DISPLACEMENT_Z);
		ElementalDofList[3] = GetGeometry()[0].pGetDof(ROTATION_X);
		ElementalDofList[4] = GetGeometry()[0].pGetDof(ROTATION_Y);
		ElementalDofList[5] = GetGeometry()[0].pGetDof(ROTATION_Z);
	}


	//************************************************************************************
	//************************************************************************************
	  void RigidBody3D::GetValuesVector(Vector& values, int Step)
	{
		if(values.size() != 6)	values.resize(6,false);
		values[0] = GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_X,Step);
		values[1] = GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_Y,Step);
		values[2] = GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_Z,Step);
		values[3] = GetGeometry()[0].FastGetSolutionStepValue(ROTATION_X,Step);
		values[4] = GetGeometry()[0].FastGetSolutionStepValue(ROTATION_Y,Step);
		values[5] = GetGeometry()[0].FastGetSolutionStepValue(ROTATION_Z,Step);
	}
	//************************************************************************************
	//************************************************************************************
	  void RigidBody3D::GetFirstDerivativesVector(Vector& values, int Step)
	{
		if(values.size() != 6)	values.resize(6,false);
		values[0] = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY_X,Step);
		values[1] = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY_Y,Step);
		values[2] = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY_Z,Step);
		values[3] = GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY_X,Step);
		values[4] = GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY_Y,Step);
		values[5] = GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY_Z,Step);
	}
	//************************************************************************************
	//************************************************************************************
	  void RigidBody3D::GetSecondDerivativesVector(Vector& values, int Step)
	{
		if(values.size() != 6)	values.resize(6,false);
		values[0] = GetGeometry()[0].FastGetSolutionStepValue(ACCELERATION_X,Step);
		values[1] = GetGeometry()[0].FastGetSolutionStepValue(ACCELERATION_Y,Step);
		values[2] = GetGeometry()[0].FastGetSolutionStepValue(ACCELERATION_Z,Step);
		values[3] = GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_ACCELERATION_X,Step);
		values[4] = GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_ACCELERATION_Y,Step);
		values[5] = GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_ACCELERATION_Z,Step);
	}

	//************************************************************************************
	//************************************************************************************
	  void RigidBody3D::MassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY
		if(rMassMatrix.size1() != 6)
			rMassMatrix.resize(6,6,false);
		
		noalias(rMassMatrix) = ZeroMatrix(6,6);
		
		rMassMatrix(0,0) = mmass;
		rMassMatrix(1,1) = mmass;
		rMassMatrix(2,2) = mmass;
		
		for(unsigned int i=0; i<3; i++)
			for(unsigned int j=0; j<3; j++)
				rMassMatrix(3+i,3+j) = mRotatedInertia(i,j);
		KRATOS_CATCH("")
	}
	
	//************************************************************************************
	//************************************************************************************
	void RigidBody3D::DampMatrix(MatrixType& rDampMatrix, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY
				
		if(rDampMatrix.size1() != 6)
				rDampMatrix.resize(6,6,false);
		
		noalias(rDampMatrix) = ZeroMatrix(6,6);
		
		const array_1d<double,3> omega = GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY);
		boost::numeric::ublas::bounded_matrix<double,3,3> cross_matrix;
		boost::numeric::ublas::bounded_matrix<double,3,3> aux;
		
		cross_matrix(0,0) =  0.00;
		cross_matrix(0,1) = -omega[2];
		cross_matrix(0,2) =  omega[1];
		cross_matrix(1,0) =  omega[2];
		cross_matrix(1,1) =  0.00;
		cross_matrix(1,2) = -omega[0];
		cross_matrix(2,0) = -omega[1];
		cross_matrix(2,1) =  omega[0];
		cross_matrix(2,2) =  0.00;
				
		noalias(aux) = prod(cross_matrix,mRotatedInertia);
		
		for(unsigned int i=0; i<3; i++)
			for(unsigned int j=0; j<3; j++)
				rDampMatrix(3+i,3+j) = aux(i,j);
		

		
		KRATOS_CATCH("")
	}	
	
	//************************************************************************************
	//************************************************************************************
	void RigidBody3D::Initialize()
	{
		//saving the last rotation used in the database
		noalias(GetGeometry()[0].GetValue(ROTATION)) = ZeroVector(3);
		noalias(mRot) = IdentityMatrix(3,3);
		
	}
	

		

} // Namespace Kratos


