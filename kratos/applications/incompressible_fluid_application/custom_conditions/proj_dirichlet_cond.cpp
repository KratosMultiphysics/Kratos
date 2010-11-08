//   
//   Project Name:        Kratos       
//   Last modified by:    $Author: paolo $
//   Date:                $Date: 2009-05-08 15:38:00 $
//   Revision:            $Revision: 1.0 $
//
// 


// System includes 


// External includes 


// Project includes 
#include "includes/define.h"
#include "custom_conditions/proj_dirichlet_cond.h"
#include "utilities/math_utils.h"
#include "incompressible_fluid_application.h"
#include "utilities/geometry_utilities.h" 
#include "utilities/math_utils.h"

namespace Kratos
{
	//************************************************************************************
	//************************************************************************************
	ProjDirichletCond::ProjDirichletCond(IndexType NewId, GeometryType::Pointer pGeometry)
		: Condition(NewId, pGeometry)
	{		
		//DO NOT ADD DOFS HERE!!!
	}

	//************************************************************************************
	//************************************************************************************
	ProjDirichletCond::ProjDirichletCond(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
		: Condition(NewId, pGeometry, pProperties)
	{
	}

	//************************************************************************************
	//************************************************************************************
	ProjDirichletCond::ProjDirichletCond(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties, array_1d<double,3> Point1, array_1d<double,3> Point2, array_1d<double,3> vel1, array_1d<double,3> vel2)
		: Condition(NewId, pGeometry, pProperties)
	{
		this->mPoint1=Point1;
		this->mPoint2=Point2;

		this->mVel1=vel1;
		this->mVel2=vel2;

		array_1d<double,3> msN;

		//and also fixing the velocity at the nodes, if the intersection is too close to them:
		unsigned int bad_vertex_index1=100;
		unsigned int bad_vertex_index2=100;
		//to check if two points arent on the same edge
		int which_edge1=100;
		int which_edge2=100;
		//first point contribution
		CalculateN_at_Point(GetGeometry(), mPoint1[0], mPoint1[1], msN);
		//KRATOS_WATCH(mPoint1)
		//KRATOS_WATCH(msN)
		//CHECK IF THE POINT IS NOT TOOO CLOSE TO THE VERTEX of destination - and identify if so, to which vertex
		double tol=0.2;
		if ( msN[0]<tol && msN[1]<tol)
			{
			bad_vertex_index1=2;	
			}
		else if ( msN[0]<tol && msN[2]<tol)
			{
			bad_vertex_index1=1;
			}
		else if ( msN[1]<tol && msN[2]<tol)
			{
			bad_vertex_index1=0;
			}
		else bad_vertex_index1=100;
		//watch out for what which:edge means
		for (int i=0;i<3;i++)
		{
		if (msN[i]<0.00000000000001)
			which_edge1=i;
		}


		CalculateN_at_Point(GetGeometry(), mPoint2[0], mPoint2[1], msN);
		//KRATOS_WATCH(mPoint2)
		//KRATOS_WATCH(msN)		
		if ( msN[0]<tol && msN[1]<tol)
			{
			bad_vertex_index2=2;			
			}
		else if ( msN[0]<tol && msN[2]<tol)
			{
			bad_vertex_index2=1;
			}
		else if ( msN[1]<tol && msN[2]<tol)
			{
			bad_vertex_index2=0;
			}
		else bad_vertex_index2=100;

		for (int i=0;i<3;i++)
		{
		if (msN[i]<0.00000000000001)
			which_edge2=i;
		}
		//check that intersections dont belong to the same point
		if (which_edge1==which_edge2)// || bad_vertex_index1<=2 || bad_vertex_index2<=2)
			{
			KRATOS_ERROR(std::logic_error,  "Two INTERSECTIONs on one edge!!!!!!!", "")// Or the intersection is close to a vertex " , "")
			}
		

		//here we fix the values of the nodes, that lie outside of origin fluid domain (i.e. IS_INTERFACE=0)
		
		for (unsigned int i=0;i<GetGeometry().size();i++)
		{
			
			
			if (GetGeometry()[i].FastGetSolutionStepValue(IS_INTERFACE)==0.0)
				{
				//if the intersection is far enough from the verreces of REAL part of fluid domain
				if (i!=bad_vertex_index1 && i!=bad_vertex_index2)	
					{
					GetGeometry()[i].FastGetSolutionStepValue(AUX_VEL_X)=GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_X);
					GetGeometry()[i].FastGetSolutionStepValue(AUX_VEL_Y)=GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_Y);
				
					GetGeometry()[i].Fix(AUX_VEL_X);
					GetGeometry()[i].Fix(AUX_VEL_Y);
					}
				else if (i==bad_vertex_index1)
					{
					KRATOS_WATCH("Case1")
					GetGeometry()[i].FastGetSolutionStepValue(IS_INTERFACE)=100.0;
					GetGeometry()[i].FastGetSolutionStepValue(AUX_VEL_X)=mVel1[0];
					GetGeometry()[i].FastGetSolutionStepValue(AUX_VEL_Y)=mVel1[1];
				
					GetGeometry()[i].Fix(AUX_VEL_X);
					GetGeometry()[i].Fix(AUX_VEL_Y);

					}
				else if (i==bad_vertex_index2)
					{
					KRATOS_WATCH("Case2")
					GetGeometry()[i].FastGetSolutionStepValue(IS_INTERFACE)=100.0;
					GetGeometry()[i].FastGetSolutionStepValue(AUX_VEL_X)=mVel2[0];
					GetGeometry()[i].FastGetSolutionStepValue(AUX_VEL_Y)=mVel2[1];
				
					GetGeometry()[i].Fix(AUX_VEL_X);
					GetGeometry()[i].Fix(AUX_VEL_Y);
					}
				}	
			//in the rest of the nodes (the ones that are lying on the fictitious part of the domain) we set aux_vel to zero
			else if (GetGeometry()[i].FastGetSolutionStepValue(IS_INTERFACE)==1.0)
				{
				GetGeometry()[i].FastGetSolutionStepValue(AUX_VEL_X)=0.0;
				GetGeometry()[i].FastGetSolutionStepValue(AUX_VEL_Y)=0.0;				
				}	
				
		}	
		
		//and apply special treatment to the cases where the intersection is too close to the node
		//if (bad_vertex_index1<=2.0 || bad_vertex_index2<=2.0)
		//	KRATOS_WATCH("You have intersections close to the nodes")		
		

		/*
		if (bad_vertex_index1<=2.0)
		{
			

			if (bad_vertex_index2<=2.0)
			{
				//check if they are close to the same vertex - if so, apply half of the sum
				if (bad_vertex_index2=bad_vertex_index1)
					{
					KRATOS_WATCH("Case1")
					GetGeometry()[bad_vertex_index1].FastGetSolutionStepValue(AUX_VEL)=(mVel1+mVel2)/2.0;
					GetGeometry()[bad_vertex_index1].Fix(AUX_VEL_X);
					GetGeometry()[bad_vertex_index1].Fix(AUX_VEL_Y);
					}
				else
					{
					KRATOS_WATCH("Case2")
					GetGeometry()[bad_vertex_index1].FastGetSolutionStepValue(AUX_VEL)=mVel1;
					GetGeometry()[bad_vertex_index1].Fix(AUX_VEL_X);
					GetGeometry()[bad_vertex_index1].Fix(AUX_VEL_Y);
					
					GetGeometry()[bad_vertex_index2].FastGetSolutionStepValue(AUX_VEL)=mVel2;
					GetGeometry()[bad_vertex_index2].Fix(AUX_VEL_X);
					GetGeometry()[bad_vertex_index2].Fix(AUX_VEL_Y);
					}

			}
			//if only one bad vertex
			else
			{
				KRATOS_WATCH("Case3")
				//directly apply the origin velocity
				GetGeometry()[bad_vertex_index1].FastGetSolutionStepValue(AUX_VEL)=mVel1;
				GetGeometry()[bad_vertex_index1].Fix(AUX_VEL_X);
				GetGeometry()[bad_vertex_index1].Fix(AUX_VEL_Y);
			}
			
		}
		else if (bad_vertex_index2<=2.0)
			{
			KRATOS_WATCH("Case4")
			GetGeometry()[bad_vertex_index2].FastGetSolutionStepValue(AUX_VEL)=mVel2;
			GetGeometry()[bad_vertex_index2].Fix(AUX_VEL_X);
			GetGeometry()[bad_vertex_index2].Fix(AUX_VEL_Y);
			}
		
		*/
	}

	Condition::Pointer ProjDirichletCond::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		return Condition::Pointer(new ProjDirichletCond(NewId, GetGeometry().Create(ThisNodes), pProperties));
	}

	ProjDirichletCond::~ProjDirichletCond()
	{
	}


	//************************************************************************************
	//************************************************************************************
	void ProjDirichletCond::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{

		/*
		if(rRightHandSideVector.size() != 6)
           		 rRightHandSideVector.resize(6,false);
		*/
		
		KRATOS_ERROR(std::logic_error,"Method not implemented!!!!","");
		
		
	}

	//************************************************************************************
	//************************************************************************************
	void ProjDirichletCond::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		array_1d<double,3> msN;
		
		if(rLeftHandSideMatrix.size1() != 6)
		{
			rLeftHandSideMatrix.resize(6,6,false);
		}
		noalias(rLeftHandSideMatrix) = ZeroMatrix(6,6); 

		if(rRightHandSideVector.size() != 6)
            rRightHandSideVector.resize(6,false);

		//save the velocity vec in order to add the rhs contribution originating from Mv
		array_1d<double,6> Vel_VEC;
		for (int i=0;i<3;i++)
		{
		Vel_VEC[2*i]=GetGeometry()[i].FastGetSolutionStepValue(AUX_VEL_X);
		Vel_VEC[2*i+1]=GetGeometry()[i].FastGetSolutionStepValue(AUX_VEL_Y);
		}
		//KRATOS_WATCH(Vel_VEC)
		/*
		const array_1d<double,3>& fp = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE_FORCE);
		
		
		rRightHandSideVector[0] = -fp[0];
		rRightHandSideVector[1] = -fp[1];
		*/
		/*
		KRATOS_WATCH(this->mPoint1)
		KRATOS_WATCH(this->mPoint2)
		KRATOS_WATCH(this->mVel1)
		KRATOS_WATCH(this->mVel2)
		*/
		//the length of the interface
		double interf_length;
		interf_length=sqrt((mPoint1[0]-mPoint2[0])*(mPoint1[0]-mPoint2[0])+(mPoint1[1]-mPoint2[1])*(mPoint1[1]-mPoint2[1]));
		//KRATOS_WATCH(interf_length)
		if (interf_length<0.0000000001)	
			KRATOS_ERROR(std::logic_error,"ZERO intersection length!!!!","");
		//double Area=GeometryUtils::CalculateVolume2D(GetGeometry());

		//we add the mass matrix to the left-hand side (we are solving: Ingeral_over_Gamma (omega*(u-u_prescribed))dGamma=0, where gamma is the interface)
		
		Matrix Mass(3,3);
		Matrix Mass1(3,3);
		//first point contribution
		CalculateN_at_Point(GetGeometry(), mPoint1[0], mPoint1[1], msN);
			
		Mass1=outer_prod(msN, trans(msN));//MathUtils<double>::TensorProduct3(Aux,Aux);

		//KRATOS_WATCH(msN)
		//KRATOS_WATCH(Mass1)
		//second point contribution
		CalculateN_at_Point(GetGeometry(), mPoint2[0], mPoint2[1], msN);
		

		Mass=Mass1+outer_prod(msN, trans(msN));//MathUtils<double>::TensorProduct3(Aux,Aux);
	
		Mass*=0.5*interf_length;
		//KRATOS_WATCH(msN)
		//KRATOS_WATCH(Mass)

		//checking with lumped mass
		
		for (int i=0;i<3;i++)
		{ 
		for (int k=0;k<3;k++)
			{
			if (i==k)
				Mass(i,k)=Mass(i,0)+Mass(i,1)+Mass(i,2);
			else 
				Mass(i,k)=0.0;
			}
		}
		

		for(unsigned int i=0;i<3;i++)
		{
			for (unsigned int k=0;k<3;k++)
			{
			rLeftHandSideMatrix(2*i,2*k)=Mass(i,k);
			rLeftHandSideMatrix(2*i+1,2*k+1)=Mass(i,k);
			}
				
		}
		//KRATOS_WATCH(rLeftHandSideMatrix)

		//and now we  assemble the RHS, which is computed using 2 Gauss points - we integrate right at the intersection points
		//so: rhs=0.5 * l *( (N1 N2 N3)_at_point1 * (v_prescribed1) + (N1 N2 N3)_at_point2 * (v_prescribed2))  , where l is the length of the interface		
		//first point contribution
		
		CalculateN_at_Point(GetGeometry(), mPoint1[0], mPoint1[1], msN);

		//first all the x_components
		rRightHandSideVector[0]=msN[0]*(mVel1[0]-(msN[0]*Vel_VEC[0]+msN[1]*Vel_VEC[2]+msN[2]*Vel_VEC[4]));
		rRightHandSideVector[2]=msN[1]*(mVel1[0]-(msN[0]*Vel_VEC[0]+msN[1]*Vel_VEC[2]+msN[2]*Vel_VEC[4]));
		rRightHandSideVector[4]=msN[2]*(mVel1[0]-(msN[0]*Vel_VEC[0]+msN[1]*Vel_VEC[2]+msN[2]*Vel_VEC[4]));
		//and now y_component
		rRightHandSideVector[1]=msN[0]*(mVel1[1]-(msN[0]*Vel_VEC[1]+msN[1]*Vel_VEC[3]+msN[2]*Vel_VEC[5]));
		rRightHandSideVector[3]=msN[1]*(mVel1[1]-(msN[0]*Vel_VEC[1]+msN[1]*Vel_VEC[3]+msN[2]*Vel_VEC[5]));
		rRightHandSideVector[5]=msN[2]*(mVel1[1]-(msN[0]*Vel_VEC[1]+msN[1]*Vel_VEC[3]+msN[2]*Vel_VEC[5]));

		//second point contribution
		CalculateN_at_Point(GetGeometry(), mPoint2[0], mPoint2[1], msN);
		//KRATOS_WATCH(msN)
		//first all the x_components
		rRightHandSideVector[0]+=msN[0]*(mVel2[0]-(msN[0]*Vel_VEC[0]+msN[1]*Vel_VEC[2]+msN[2]*Vel_VEC[4]));
		rRightHandSideVector[2]+=msN[1]*(mVel2[0]-(msN[0]*Vel_VEC[0]+msN[1]*Vel_VEC[2]+msN[2]*Vel_VEC[4]));
		rRightHandSideVector[4]+=msN[2]*(mVel2[0]-(msN[0]*Vel_VEC[0]+msN[1]*Vel_VEC[2]+msN[2]*Vel_VEC[4]));
		//and now y_component
		rRightHandSideVector[1]+=msN[0]*(mVel2[1]-(msN[0]*Vel_VEC[1]+msN[1]*Vel_VEC[3]+msN[2]*Vel_VEC[5]));
		rRightHandSideVector[3]+=msN[1]*(mVel2[1]-(msN[0]*Vel_VEC[1]+msN[1]*Vel_VEC[3]+msN[2]*Vel_VEC[5]));
		rRightHandSideVector[5]+=msN[2]*(mVel2[1]-(msN[0]*Vel_VEC[1]+msN[1]*Vel_VEC[3]+msN[2]*Vel_VEC[5]));

		//completing integration using two Gauss points
		rRightHandSideVector*=0.5*interf_length;

		
		//KRATOS_WATCH(mVel1)
		//KRATOS_WATCH(mVel2) 
		//KRATOS_WATCH(interf_length)
		//KRATOS_WATCH("RHS inside of condition")
		//KRATOS_WATCH(rRightHandSideVector)

		//and finally imposing the Dirichlet on the nodes, that lie inside the fluid (IS_INTERFACE=0)
		//and also onto the nodes that are too close to the interface (those are defined by bad_vertex_index)
		
		
		//first we apply Dirichlet to "good" conditions (no bad vertices)
		/*
		if ( bad_vertex_index1>2 && bad_vertex_index2>2 )
		{
		for (int i=0;i<GetGeometry().size();i++)
			{
			KRATOS_WATCH(GetGeometry()[i].FastGetSolutionStepValue(IS_INTERFACE))
			if (GetGeometry()[i].FastGetSolutionStepValue(IS_INTERFACE)==0.0)			
				{
				//KRATOS_WATCH("Fixing velocity")

				//CHECK that it is not yet fixed by the neighboring element
					GetGeometry()[i].FastGetSolutionStepValue(AUX_VEL_X)=GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_X);
					GetGeometry()[i].FastGetSolutionStepValue(AUX_VEL_Y)=GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_Y);

					KRATOS_WATCH(GetGeometry()[i].FastGetSolutionStepValue(VELOCITY))
					KRATOS_WATCH(GetGeometry()[i].FastGetSolutionStepValue(VELOCITY))
					
					KRATOS_WATCH(GetGeometry()[i].FastGetSolutionStepValue(AUX_VEL_X))
					KRATOS_WATCH(GetGeometry()[i].FastGetSolutionStepValue(AUX_VEL_Y))
					GetGeometry()[i].FastGetSolutionStepValue(AUX_VEL_X)=0.0;
					GetGeometry()[i].FastGetSolutionStepValue(AUX_VEL_Y)=0.0;
					KRATOS_WATCH(GetGeometry()[i].FastGetSolutionStepValue(AUX_VEL_X))
					KRATOS_WATCH(GetGeometry()[i].FastGetSolutionStepValue(AUX_VEL_Y))
					
					//KRATOS_WATCH(GetGeometry()[i].FastGetSolutionStepValue(AUX_VEL))

					GetGeometry()[i].Fix(AUX_VEL_X);
					GetGeometry()[i].Fix(AUX_VEL_Y);
				
				}
			}
		}
		KRATOS_WATCH("OK1")
		
		if (bad_vertex_index1<=2)
		{
		GetGeometry()[bad_vertex_index1].FastGetSolutionStepValue(AUX_VEL_X)=mVel1[0];
		GetGeometry()[bad_vertex_index1].FastGetSolutionStepValue(AUX_VEL_Y)=mVel1[1];
		GetGeometry()[bad_vertex_index1].Fix(AUX_VEL_X);
		GetGeometry()[bad_vertex_index1].Fix(AUX_VEL_Y);	
		KRATOS_WATCH("BAD VERTEX FOUNDDDDDDDDDDDD!!!")			
		}
		if (bad_vertex_index2<=2)
		{
		KRATOS_WATCH("BAD VERTEX FOUNDDDDDDDDDDDD!!!")			
		GetGeometry()[bad_vertex_index2].FastGetSolutionStepValue(AUX_VEL_X)=mVel2[0];
		GetGeometry()[bad_vertex_index2].FastGetSolutionStepValue(AUX_VEL_Y)=mVel2[1];
		GetGeometry()[bad_vertex_index2].Fix(AUX_VEL_X);
		GetGeometry()[bad_vertex_index2].Fix(AUX_VEL_Y);				
		}
		
		//////////////////////////////////////////////////////////////	
		/// CHECK THIS ONE AND MAYBE UNCOMMENT!!!!!!
		if (bad_vertex_index1<=2 || bad_vertex_index2<=2)
		{
			for (int i=0;i<3;i++)
			{
			if (i!=bad_vertex_index1 && i!=bad_vertex_index2)
				{
				if (GetGeometry()[i].FastGetSolutionStepValue(IS_INTERFACE)==0.0)
					{
					GetGeometry()[i].FastGetSolutionStepValue(AUX_VEL)=GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
			
					GetGeometry()[i].Fix(AUX_VEL_X);
					GetGeometry()[i].Fix(AUX_VEL_Y);
					}
				else 
					{
					GetGeometry()[i].FastGetSolutionStepValue(AUX_VEL)=0.5*(mVel1+mVel2);
			
					GetGeometry()[i].Fix(AUX_VEL_X);
					GetGeometry()[i].Fix(AUX_VEL_Y);
					}

				
				}
			}
		}
		*/
////////////////////////////////////////////

		











			/*

		//if we have both bad vertices
		if (bad_vertex_index2<=2 && bad_vertex_index1<=2)
		{
			GetGeometry()[bad_vertex_index1].FastGetSolutionStepValue(AUX_VEL)=mVel1;
			GetGeometry()[bad_vertex_index1].Fix(AUX_VEL_X);
			GetGeometry()[bad_vertex_index1].Fix(AUX_VEL_Y);
			KRATOS_WATCH("Fixing the vel at bad vertices")
	
			GetGeometry()[bad_vertex_index2].FastGetSolutionStepValue(AUX_VEL)=mVel2;
			GetGeometry()[bad_vertex_index2].Fix(AUX_VEL_X);
			GetGeometry()[bad_vertex_index2].Fix(AUX_VEL_Y);
			KRATOS_WATCH("Fixing the vel at bad vertices")
			if (bad_vertex_index1==bad_vertex_index2)
				{
				KRATOS_WATCH("Both intersections at close to the same node")
				KRATOS_WATCH(mVel1)
				KRATOS_WATCH(mVel2)
				}

			for (int i=0;i<GetGeometry().size();i++)
				{
				if (i!=bad_vertex_index1 && i!=bad_vertex_index2)
					{
					//if the intersections are inside the origin  - >fix the rest to the destination vel
					if (GetGeometry()[bad_vertex_index1].FastGetSolutionStepValue(IS_INTERFACE)==1.0 && GetGeometry()[bad_vertex_index2].FastGetSolutionStepValue(IS_INTERFACE)==1.0)
						{
						GetGeometry()[i].FastGetSolutionStepValue(AUX_VEL)=GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
						GetGeometry()[i].Fix(AUX_VEL_X);
						GetGeometry()[i].Fix(AUX_VEL_Y);
						}
					//check if you can fix to zero the inner nodes!!!!!!!
					else if (GetGeometry()[bad_vertex_index1].FastGetSolutionStepValue(IS_INTERFACE)==0.0 && GetGeometry()[bad_vertex_index2].FastGetSolutionStepValue(IS_INTERFACE)==0.0)
						{
						GetGeometry()[i].FastGetSolutionStepValue(AUX_VEL)=0.5*(mVel1+mVel2);
						GetGeometry()[i].Fix(AUX_VEL_X);
						GetGeometry()[i].Fix(AUX_VEL_Y);
						}
					//strange case if one is is_interface and another one - no.. then throw an error
					else if ( GetGeometry()[bad_vertex_index1].FastGetSolutionStepValue(IS_INTERFACE)+GetGeometry()[bad_vertex_index2].FastGetSolutionStepValue(IS_INTERFACE)==1.0 ) 
						{
						KRATOS_ERROR(std::logic_error,  "ONE BAD IN ONE BAD Out!!!!!!!1!!!! " , "")
						}
					
					}

				}

		}
		else if ((bad_vertex_index1<=2 && bad_vertex_index2>2)||(bad_vertex_index2<=2 && bad_vertex_index2>1))
			{
			int bad_vertex_index;
			array_1d<double,3> Vel_at_bad_vertex;
			array_1d<double,3> Vel_at_good_inters;
			if (bad_vertex_index1<=2)
				{
				bad_vertex_index=bad_vertex_index1;
				Vel_at_bad_vertex=mVel1;
				Vel_at_good_inters=mVel2;
				}
				
			else if (bad_vertex_index2<=2)
				{
				bad_vertex_index=bad_vertex_index2;
				Vel_at_bad_vertex=mVel2;
				Vel_at_good_inters=mVel1;
				}			
			GetGeometry()[bad_vertex_index].FastGetSolutionStepValue(AUX_VEL)=Vel_at_bad_vertex;
			GetGeometry()[bad_vertex_index].Fix(AUX_VEL_X);
			GetGeometry()[bad_vertex_index].Fix(AUX_VEL_Y);

			for (int i=0;i<GetGeometry().size();i++)
				{
				if (i!=bad_vertex_index)
					{
					if (GetGeometry()[i].FastGetSolutionStepValue(IS_INTERFACE)==0.0)
						{
						GetGeometry()[i].FastGetSolutionStepValue(AUX_VEL)=GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
						}
					else if (GetGeometry()[i].FastGetSolutionStepValue(IS_INTERFACE)==1.0)
						{
						//PUT HERE BETTER APPROXIMATION----
						//double x=GetGeometry()
						GetGeometry()[i].FastGetSolutionStepValue(AUX_VEL)=Vel_at_good_inters;
						KRATOS_WATCH(GetGeometry()[i].FastGetSolutionStepValue(AUX_VEL))
						}
					}
				}
			//KRATOS_ERROR(std::logic_error,  "CHEck, you just have one bad vertex!!! and it screws everything1!!!! " , "")
			}
		
		
		///////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////
		//											//
		//			THIS IS WORKING!!!!						//
		//	
		/////////////////////////////////////////////////////////////////////////////////////////
		//now the case that a condition contains just one bad vertex
		
		if (bad_vertex_index2<=2 || bad_vertex_index1<=2)
		{
		for (int i=0;i<GetGeometry().size();i++)
			{
			GetGeometry()[i].FastGetSolutionStepValue(AUX_VEL)=0.5*(mVel2+mVel1);
			GetGeometry()[i].Fix(AUX_VEL_X);
			GetGeometry()[i].Fix(AUX_VEL_Y);
			KRATOS_WATCH("Fixing the vel at bad vertices")
			}
		}
		

		///////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////
		
		//and now we are directly applying Dirichlet onto the nodes that are too close to the intersection (those are defined by bad_vertex_index)
		
		if (bad_vertex_index1<=2)
		{
			if (bad_vertex_index2<=2)
			{
				//check if they are close to the same vertex - if so, apply half of the sum
				if (bad_vertex_index2=bad_vertex_index1)
					{
					GetGeometry()[bad_vertex_index1].FastGetSolutionStepValue(AUX_VEL)=(mVel1+mVel2)/2.0;
					GetGeometry()[bad_vertex_index1].Fix(AUX_VEL_X);
					GetGeometry()[bad_vertex_index1].Fix(AUX_VEL_Y);
					}
				else
					{
					GetGeometry()[bad_vertex_index1].FastGetSolutionStepValue(AUX_VEL)=mVel1;
					GetGeometry()[bad_vertex_index1].Fix(AUX_VEL_X);
					GetGeometry()[bad_vertex_index1].Fix(AUX_VEL_Y);
					
					GetGeometry()[bad_vertex_index2].FastGetSolutionStepValue(AUX_VEL)=mVel2;
					GetGeometry()[bad_vertex_index2].Fix(AUX_VEL_X);
					GetGeometry()[bad_vertex_index2].Fix(AUX_VEL_Y);
					}

			}
			//if only one bad vertex
			else
			{
				//directly apply the origin velocity
				GetGeometry()[bad_vertex_index1].FastGetSolutionStepValue(AUX_VEL)=mVel1;
				GetGeometry()[bad_vertex_index1].Fix(AUX_VEL_X);
				GetGeometry()[bad_vertex_index1].Fix(AUX_VEL_Y);
//			}
//		}
//		else if (bad_vertex_index2<=2)
//		{
//		GetGeometry()[bad_vertex_index2].FastGetSolutionStepValue(AUX_VEL)=mVel2;
//		GetGeometry()[bad_vertex_index2].Fix(AUX_VEL_X);
//		GetGeometry()[bad_vertex_index2].Fix(AUX_VEL_Y);
//		}
		
		*/
		
		
	}
	//************************************************************************************
	//************************************************************************************
	void ProjDirichletCond::CalculateN_at_Point(Element::GeometryType& geom, const double xc, const double yc, array_1d<double,3>& N_at_c)
	{
			//first we calculate the area of the whole triangle			
			double x10 = geom[1].X() - geom[0].X();
			double y10 = geom[1].Y() - geom[0].Y();
			
			double x20 = geom[2].X() - geom[0].X();
			double y20 = geom[2].Y() - geom[0].Y();
			
			double detJ = x10 * y20-y10 * x20;
			double totArea=0.5*detJ;			
			//and now we calculate the areas of three respective triangle, that (xc,yc) divide the original one into
			// xc, 0, 1
			double x0c = geom[0].X() - xc ;
			double y0c = geom[0].Y() - yc ;
			
			double x1c = geom[1].X() - xc;
			double y1c = geom[1].Y() - yc;

			double x2c = geom[2].X() - xc;
			double y2c = geom[2].Y() - yc;
			//xc, 0, 1
			detJ= x0c * y1c - y0c * x1c;
			double Area2 = 0.5*detJ;
			//xc, 0, 2
			detJ= x0c * y2c - y0c * x2c;
			double Area1 = 0.5*detJ;
			//xc, 1, 2
			detJ= x1c * y2c - y1c * x2c;
			double Area0 = 0.5*detJ;

			if (totArea<0.00000000000000001)
				KRATOS_ERROR(std::logic_error,  "Your element Proj DIrichlet Cond has a zero area!!!! " , "");
			//and now we fill in the array of shape functions values:
			// 1 0 2
			N_at_c[0]=fabs(Area0/totArea);
			N_at_c[1]=fabs(Area1/totArea);
			N_at_c[2]=fabs(Area2/totArea);
			if (  (N_at_c[0]<0.05 && N_at_c[1]<0.05) || (N_at_c[0]<0.05 && N_at_c[2]<0.05) || (N_at_c[2]<0.05 && N_at_c[1]<0.05))			
			KRATOS_WATCH("Dangerous VERTICES!!!")		
			//KRATOS_ERROR(std::logic_error,  "Too close to the node is the INTERSECTION!!!! " , "")

	}


	//************************************************************************************
	//************************************************************************************
	void ProjDirichletCond::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
	{
		if(rResult.size() != 6)
			rResult.resize(6);

		rResult[0] = (GetGeometry()[0].GetDof(AUX_VEL_X)).EquationId();
		rResult[1] = (GetGeometry()[0].GetDof(AUX_VEL_Y)).EquationId();

		rResult[2] = (GetGeometry()[1].GetDof(AUX_VEL_X)).EquationId();
		rResult[3] = (GetGeometry()[1].GetDof(AUX_VEL_Y)).EquationId();

		rResult[4] = (GetGeometry()[2].GetDof(AUX_VEL_X)).EquationId();
		rResult[5] = (GetGeometry()[2].GetDof(AUX_VEL_Y)).EquationId();

	}

	//************************************************************************************
	//************************************************************************************
	void ProjDirichletCond::GetDofList(DofsVectorType& ConditionalDofList,ProcessInfo& CurrentProcessInfo)
	{
		if(ConditionalDofList.size() != 6)
			ConditionalDofList.resize(6);

		ConditionalDofList[0] = (GetGeometry()[0].pGetDof(AUX_VEL_X));
		ConditionalDofList[1] = (GetGeometry()[0].pGetDof(AUX_VEL_Y));
		
		ConditionalDofList[2] = (GetGeometry()[1].pGetDof(AUX_VEL_X));
		ConditionalDofList[3] = (GetGeometry()[1].pGetDof(AUX_VEL_Y));

		ConditionalDofList[4] = (GetGeometry()[2].pGetDof(AUX_VEL_X));
		ConditionalDofList[5] = (GetGeometry()[2].pGetDof(AUX_VEL_Y));

	}
	
} // Namespace Kratos


