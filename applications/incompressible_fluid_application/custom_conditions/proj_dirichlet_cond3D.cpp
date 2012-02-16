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
#include "custom_conditions/proj_dirichlet_cond3D.h"
#include "utilities/math_utils.h"
#include "incompressible_fluid_application.h"
#include "utilities/geometry_utilities.h" 
#include "utilities/math_utils.h"

namespace Kratos
{
	//************************************************************************************
	//************************************************************************************
	ProjDirichletCond3D::ProjDirichletCond3D(IndexType NewId, GeometryType::Pointer pGeometry)
		: Condition(NewId, pGeometry)
	{		
		//DO NOT ADD DOFS HERE!!!
	}

	//************************************************************************************
	//************************************************************************************
	ProjDirichletCond3D::ProjDirichletCond3D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
		: Condition(NewId, pGeometry, pProperties)
	{
	}

	//************************************************************************************
	//************************************************************************************
	//This is a constructor for the intersections consisting of 3 points
	ProjDirichletCond3D::ProjDirichletCond3D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties, array_1d<double,3> Point1, array_1d<double,3> Point2,array_1d<double,3> Point3, array_1d<double,3> 		vel1, array_1d<double,3> vel2, array_1d<double,3> vel3)
		: Condition(NewId, pGeometry, pProperties)
	{
		
		this->mPoint1=Point1;
		this->mPoint2=Point2;
		this->mPoint3=Point3;

		this->mPoint4=ZeroVector(3);

		this->mVel1=vel1;
		this->mVel2=vel2;
		this->mVel3=vel3;

		this->mVel4=ZeroVector(3);

		this->mNumber_of_intersections=3;

		array_1d<double,4> msN;

		//and also fixing the velocity at the nodes, if the intersection is too close to them:
		unsigned int bad_vertex_index1=100;
		unsigned int bad_vertex_index2=100;
		unsigned int bad_vertex_index3=100;
		
		//first point contribution
		CalculateN_at_Point(GetGeometry(), mPoint1[0], mPoint1[1], mPoint1[2], msN);
		//CHECK IF THE POINT IS NOT TOOO CLOSE TO THE VERTEX of destination - and identify if so, to which vertex
		double tol=0.1;
		if ( msN[0]<tol && msN[1]<tol && msN[2]<tol)			
			bad_vertex_index1=3;				
		else if ( msN[0]<tol && msN[2]<tol && msN[3]<tol)
			bad_vertex_index1=1;			
		else if ( msN[1]<tol && msN[2]<tol && msN[3]<tol)
			bad_vertex_index1=0;
		else if ( msN[0]<tol && msN[1]<tol && msN[3]<tol)
			bad_vertex_index1=2;			
		else bad_vertex_index1=100;
		
		CalculateN_at_Point(GetGeometry(), mPoint2[0], mPoint2[1], mPoint2[2], msN);
		if ( msN[0]<tol && msN[1]<tol && msN[2]<tol)
			bad_vertex_index2=3;	
		else if ( msN[0]<tol && msN[2]<tol && msN[3]<tol)
			bad_vertex_index2=1;
		else if ( msN[1]<tol && msN[2]<tol && msN[3]<tol)
			bad_vertex_index2=0;
		else if ( msN[0]<tol && msN[1]<tol && msN[3]<tol)
			bad_vertex_index2=2;
		else bad_vertex_index2=100;

		CalculateN_at_Point(GetGeometry(), mPoint3[0], mPoint3[1], mPoint3[2], msN);
		if ( msN[0]<tol && msN[1]<tol && msN[2]<tol)
			bad_vertex_index3=3;	
		else if ( msN[0]<tol && msN[2]<tol && msN[3]<tol)
			bad_vertex_index3=1;
		else if ( msN[1]<tol && msN[2]<tol && msN[3]<tol)
			bad_vertex_index3=0;
		else if ( msN[0]<tol && msN[1]<tol && msN[3]<tol)
			bad_vertex_index3=2;			
		else bad_vertex_index3=100;	
			
		//here we fix the values of the nodes, that lie outside of origin domain (i.e. IS_INTERFACE=0), i.e. IT IS THE REAL FLUID
		for (unsigned int i=0;i<GetGeometry().size();i++)
		{
						
			if (GetGeometry()[i].FastGetSolutionStepValue(IS_INTERFACE)==0.0)
				{
					
				//if the intersection is far enough from the verreces of REAL part of fluid domain
				if (i!=bad_vertex_index1 && i!=bad_vertex_index2 && i!=bad_vertex_index3)	
					{
					GetGeometry()[i].FastGetSolutionStepValue(AUX_VEL_X)=GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_X);
					GetGeometry()[i].FastGetSolutionStepValue(AUX_VEL_Y)=GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_Y);
					GetGeometry()[i].FastGetSolutionStepValue(AUX_VEL_Z)=GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_Z);					
					}
				else if (i==bad_vertex_index1)
					{
					KRATOS_WATCH("Case1")
					GetGeometry()[i].FastGetSolutionStepValue(IS_INTERFACE)=100.0;
					GetGeometry()[i].FastGetSolutionStepValue(AUX_VEL_X)=mVel1[0];
					GetGeometry()[i].FastGetSolutionStepValue(AUX_VEL_Y)=mVel1[1];
					GetGeometry()[i].FastGetSolutionStepValue(AUX_VEL_Z)=mVel1[2];				
					}
				else if (i==bad_vertex_index2)
					{
					KRATOS_WATCH("Case2")
					GetGeometry()[i].FastGetSolutionStepValue(IS_INTERFACE)=100.0;
					GetGeometry()[i].FastGetSolutionStepValue(AUX_VEL_X)=mVel2[0];
					GetGeometry()[i].FastGetSolutionStepValue(AUX_VEL_Y)=mVel2[1];
					GetGeometry()[i].FastGetSolutionStepValue(AUX_VEL_Z)=mVel2[2];
					}
				else if (i==bad_vertex_index3)
					{
					KRATOS_WATCH("Case3")
					GetGeometry()[i].FastGetSolutionStepValue(IS_INTERFACE)=100.0;
					GetGeometry()[i].FastGetSolutionStepValue(AUX_VEL_X)=mVel3[0];
					GetGeometry()[i].FastGetSolutionStepValue(AUX_VEL_Y)=mVel3[1];
					GetGeometry()[i].FastGetSolutionStepValue(AUX_VEL_Z)=mVel3[2];
					}
					
					GetGeometry()[i].Fix(AUX_VEL_X);
					GetGeometry()[i].Fix(AUX_VEL_Y);
					GetGeometry()[i].Fix(AUX_VEL_Z);
				}	
			
			//in the rest of the nodes (the ones that are lying in the fictitious part of the domain) we set aux_vel to zero
			else if (GetGeometry()[i].FastGetSolutionStepValue(IS_INTERFACE)==1.0)
				{
				GetGeometry()[i].FastGetSolutionStepValue(AUX_VEL_X)=0.0;
				GetGeometry()[i].FastGetSolutionStepValue(AUX_VEL_Y)=0.0;				
				GetGeometry()[i].FastGetSolutionStepValue(AUX_VEL_Z)=0.0;
				}	
				
		}	
		
		
	}
	//This is a constructor for the intersections consisting of 4 points
	ProjDirichletCond3D::ProjDirichletCond3D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties, array_1d<double,3> Point1, array_1d<double,3> Point2,array_1d<double,3> Point3, array_1d<double,3> Point4, array_1d<double,3> vel1, array_1d<double,3> vel2, array_1d<double,3> vel3, array_1d<double,3> vel4)
		: Condition(NewId, pGeometry, pProperties)
	{
		this->mPoint1=Point1;
		this->mPoint2=Point2;
		this->mPoint3=Point3;
		this->mPoint4=Point4;

		this->mVel1=vel1;
		this->mVel2=vel2;
		this->mVel3=vel3;
		this->mVel4=vel4;

		this->mNumber_of_intersections=4;

		array_1d<double,4> msN;

		//and also fixing the velocity at the nodes, if the intersection is too close to them:
		unsigned int bad_vertex_index1=100;
		unsigned int bad_vertex_index2=100;
		unsigned int bad_vertex_index3=100;
		unsigned int bad_vertex_index4=100;
		
		//first point contribution
		CalculateN_at_Point(GetGeometry(), mPoint1[0], mPoint1[1], mPoint1[2], msN);
		//KRATOS_WATCH(mPoint1)
		//KRATOS_WATCH(msN)
		//CHECK IF THE POINT IS NOT TOOO CLOSE TO THE VERTEX of destination - and identify if so, to which vertex
		double tol=0.1;
		if ( msN[0]<tol && msN[1]<tol && msN[2]<tol)
			bad_vertex_index1=3;	
		else if ( msN[0]<tol && msN[2]<tol && msN[3]<tol)			
			bad_vertex_index1=1;			
		else if ( msN[1]<tol && msN[2]<tol && msN[3]<tol)
			bad_vertex_index1=0;
		else if ( msN[0]<tol && msN[1]<tol && msN[3]<tol)
			bad_vertex_index1=2;			
		else bad_vertex_index1=100;
		
		CalculateN_at_Point(GetGeometry(), mPoint2[0], mPoint2[1], mPoint2[2], msN);
		if ( msN[0]<tol && msN[1]<tol && msN[2]<tol)
			bad_vertex_index2=3;				
		else if ( msN[0]<tol && msN[2]<tol && msN[3]<tol)
			bad_vertex_index2=1;
		else if ( msN[1]<tol && msN[2]<tol && msN[3]<tol)
			bad_vertex_index2=0;
		else if ( msN[0]<tol && msN[1]<tol && msN[3]<tol)
			bad_vertex_index2=2;		
		else bad_vertex_index2=100;
		
		CalculateN_at_Point(GetGeometry(), mPoint3[0], mPoint3[1], mPoint3[2], msN);
		if ( msN[0]<tol && msN[1]<tol && msN[2])
			bad_vertex_index3=3;	
		else if ( msN[0]<tol && msN[2]<tol && msN[3]<tol)
			bad_vertex_index3=1;
		else if ( msN[1]<tol && msN[2]<tol && msN[3]<tol)
			bad_vertex_index3=0;
		else if ( msN[0]<tol && msN[1]<tol && msN[3]<tol)
			bad_vertex_index3=2;
		else bad_vertex_index3=100;

		CalculateN_at_Point(GetGeometry(), mPoint4[0], mPoint4[1], mPoint4[2], msN);
		
		if ( msN[0]<tol && msN[1]<tol && msN[2]<tol)
			bad_vertex_index4=3;	
		else if ( msN[0]<tol && msN[2]<tol && msN[3]<tol)
			bad_vertex_index4=1;
		else if ( msN[1]<tol && msN[2]<tol && msN[3]<tol)
			bad_vertex_index4=0;
		else if ( msN[0]<tol && msN[1]<tol && msN[3]<tol)
			bad_vertex_index4=2;
		else bad_vertex_index4=100;
		
		
		
		//here we fix the values of the nodes, that lie outside of origin domain (i.e. IS_INTERFACE=0), i.e. IT IS THE REAL FLUID
		
		for (unsigned int i=0;i<GetGeometry().size();i++)
		{
						
			if (GetGeometry()[i].FastGetSolutionStepValue(IS_INTERFACE)==0.0)
				{
					
				//if the intersection is far enough from the verreces of REAL part of fluid domain
				if (i!=bad_vertex_index1 && i!=bad_vertex_index2 && i!=bad_vertex_index3)	
					{
					GetGeometry()[i].FastGetSolutionStepValue(AUX_VEL_X)=GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_X);
					GetGeometry()[i].FastGetSolutionStepValue(AUX_VEL_Y)=GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_Y);
					GetGeometry()[i].FastGetSolutionStepValue(AUX_VEL_Z)=GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_Z);					
					}
				else if (i==bad_vertex_index1)
					{
					GetGeometry()[i].FastGetSolutionStepValue(IS_INTERFACE)=100.0;
					GetGeometry()[i].FastGetSolutionStepValue(AUX_VEL_X)=mVel1[0];
					GetGeometry()[i].FastGetSolutionStepValue(AUX_VEL_Y)=mVel1[1];
					GetGeometry()[i].FastGetSolutionStepValue(AUX_VEL_Z)=mVel1[2];				
					}
				else if (i==bad_vertex_index2)
					{					
					GetGeometry()[i].FastGetSolutionStepValue(IS_INTERFACE)=100.0;
					GetGeometry()[i].FastGetSolutionStepValue(AUX_VEL_X)=mVel2[0];
					GetGeometry()[i].FastGetSolutionStepValue(AUX_VEL_Y)=mVel2[1];
					GetGeometry()[i].FastGetSolutionStepValue(AUX_VEL_Z)=mVel2[2];
					}
				else if (i==bad_vertex_index3)
					{					
					GetGeometry()[i].FastGetSolutionStepValue(IS_INTERFACE)=100.0;
					GetGeometry()[i].FastGetSolutionStepValue(AUX_VEL_X)=mVel3[0];
					GetGeometry()[i].FastGetSolutionStepValue(AUX_VEL_Y)=mVel3[1];
					GetGeometry()[i].FastGetSolutionStepValue(AUX_VEL_Z)=mVel3[2];
					}					
				else if (i==bad_vertex_index4)
					{					
					GetGeometry()[i].FastGetSolutionStepValue(IS_INTERFACE)=100.0;
					GetGeometry()[i].FastGetSolutionStepValue(AUX_VEL_X)=mVel4[0];
					GetGeometry()[i].FastGetSolutionStepValue(AUX_VEL_Y)=mVel4[1];
					GetGeometry()[i].FastGetSolutionStepValue(AUX_VEL_Z)=mVel4[2];
					}
					
					GetGeometry()[i].Fix(AUX_VEL_X);
					GetGeometry()[i].Fix(AUX_VEL_Y);
					GetGeometry()[i].Fix(AUX_VEL_Z);
				}	
			//in the rest of the nodes (the ones that are lying on the fictitious part of the domain) we set aux_vel to zero
			else if (GetGeometry()[i].FastGetSolutionStepValue(IS_INTERFACE)==1.0)
				{
				GetGeometry()[i].FastGetSolutionStepValue(AUX_VEL_X)=0.0;
				GetGeometry()[i].FastGetSolutionStepValue(AUX_VEL_Y)=0.0;				
				GetGeometry()[i].FastGetSolutionStepValue(AUX_VEL_Z)=0.0;
				}	
				
		}	
		
		
		
	}

	Condition::Pointer ProjDirichletCond3D::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		return Condition::Pointer(new ProjDirichletCond3D(NewId, GetGeometry().Create(ThisNodes), pProperties));
	}

	ProjDirichletCond3D::~ProjDirichletCond3D()
	{
	}


	//************************************************************************************
	//************************************************************************************
	void ProjDirichletCond3D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{

		/*
		if(rRightHandSideVector.size() != 12)
           		 rRightHandSideVector.resize(12,false);
		*/
		
		KRATOS_ERROR(std::logic_error,"Method not implemented!!!!","");
		
		
	}

	//************************************************************************************
	//************************************************************************************
	void ProjDirichletCond3D::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		array_1d<double,4> msN;
		
		if(rLeftHandSideMatrix.size1() != 12)
		{
			rLeftHandSideMatrix.resize(12,12,false);
		}
		noalias(rLeftHandSideMatrix) = ZeroMatrix(12,12); 

		if(rRightHandSideVector.size() != 12)
            		rRightHandSideVector.resize(12,false);

		//save the velocity vec in order to add the rhs contribution originating from Mv
		array_1d<double,12> Vel_VEC=ZeroVector(12);
		for (int i=0;i<4;i++)
		{
		Vel_VEC[3*i]=GetGeometry()[i].FastGetSolutionStepValue(AUX_VEL_X);
		Vel_VEC[3*i+1]=GetGeometry()[i].FastGetSolutionStepValue(AUX_VEL_Y);
		Vel_VEC[3*i+2]=GetGeometry()[i].FastGetSolutionStepValue(AUX_VEL_Y);
		}
		KRATOS_WATCH(Vel_VEC)
					
		double inters_area =  0.0;
		//if the intersection is a quadrilateral - its a sum of two triangles
		double inters_area2=0.0;
		
		inters_area=CalculateTriangleArea3D(mPoint1, mPoint2, mPoint3);

		if (inters_area<0.000000000000001)	
			KRATOS_ERROR(std::logic_error,"ZERO intersection AREA!!!!","");

		//this vector stores the "vertices" of the second triangle (in case there are 4 intersection points)
		//std::vector< array_1d<double,3> > PointsOfSecondTriangle(4, array_1d<double,3>(3,0));//Vector(3));
		
		std::vector<array_1d<double,3> > PointsOfSecondTriangle;
		std::vector<array_1d<double,3> > VelsOfSecondTriangle;
		PointsOfSecondTriangle.reserve(3);
		VelsOfSecondTriangle.reserve(3);
		if (this->mNumber_of_intersections==4)
			{
			//the fourth point always belongs to the secodn triangle.. the other two ones we need to find
			PointsOfSecondTriangle.push_back(mPoint4);
			VelsOfSecondTriangle.push_back(mVel4);
			//if there are 4 points, we need to understand which is the second triangle (connectivities)
			double a41=Length(mPoint4, mPoint1);
			double a42=Length(mPoint4, mPoint2);
			double a43=Length(mPoint4, mPoint3);			
			
			if (a41>a42 && a41>a43)
				{
				PointsOfSecondTriangle.push_back(mPoint2);
				PointsOfSecondTriangle.push_back(mPoint3);
				VelsOfSecondTriangle.push_back(mVel2);
				VelsOfSecondTriangle.push_back(mVel3);
				}
			else if (a42>a41 && a42>a43)	
				{
				PointsOfSecondTriangle.push_back(mPoint1);
				PointsOfSecondTriangle.push_back(mPoint3);
				VelsOfSecondTriangle.push_back(mVel1);
				VelsOfSecondTriangle.push_back(mVel3);				
				}
			else if (a43>a41 && a43>a42)
				{
				PointsOfSecondTriangle.push_back(mPoint1);
				PointsOfSecondTriangle.push_back(mPoint2);
				VelsOfSecondTriangle.push_back(mVel1);
				VelsOfSecondTriangle.push_back(mVel2);
				}
			else 
				KRATOS_ERROR(std::logic_error,"Intersection has 4 points, but something is wrong!!!! Maybe the distance between the points is the same","");

			
			inters_area2=CalculateTriangleArea3D(PointsOfSecondTriangle[0], PointsOfSecondTriangle[1], PointsOfSecondTriangle[2]);
			
			if (inters_area2==0.0)
				KRATOS_ERROR(std::logic_error,"Intersection area of the second triangle is ZERO", "");
				
			}		
		
		
		//we add the mass matrix to the left-hand side (we are solving: Ingeral_over_Gamma (omega*(u-u_prescribed))dGamma=0, where gamma is the interface)
		
		Matrix Mass(4,4);
		Matrix Mass1(4,4);
		Matrix Mass2(4,4);
		Matrix Mass3(4,4);
		//first point contribution
		CalculateN_at_Point(GetGeometry(), mPoint1[0], mPoint1[1], mPoint1[2], msN);			
		Mass1=outer_prod(msN, trans(msN));//MathUtils<double>::TensorProduct3(Aux,Aux);

		//second point contribution
		CalculateN_at_Point(GetGeometry(), mPoint2[0], mPoint2[1], mPoint2[2], msN);
		Mass2=outer_prod(msN, trans(msN));//MathUtils<double>::TensorProduct3(Aux,Aux);

		//third point contribution
		CalculateN_at_Point(GetGeometry(), mPoint3[0], mPoint3[1], mPoint3[2], msN);
		Mass3=outer_prod(msN, trans(msN));

		Mass=inters_area*0.3333333333333333333333333*(Mass1+Mass2+Mass3);
		KRATOS_WATCH(Mass)

		if (this->mNumber_of_intersections==4)
			{
			CalculateN_at_Point(GetGeometry(), PointsOfSecondTriangle[0][0], PointsOfSecondTriangle[0][1], PointsOfSecondTriangle[0][2], msN);
			Mass1=outer_prod(msN, trans(msN));
			CalculateN_at_Point(GetGeometry(), PointsOfSecondTriangle[1][0], PointsOfSecondTriangle[1][1], PointsOfSecondTriangle[1][2], msN);
			Mass2=outer_prod(msN, trans(msN));
			CalculateN_at_Point(GetGeometry(), PointsOfSecondTriangle[2][0], PointsOfSecondTriangle[2][1], PointsOfSecondTriangle[2][2], msN);
			Mass3=outer_prod(msN, trans(msN));

			Mass+=inters_area2*0.333333333333333333*(Mass1+Mass2+Mass3);

			}
		//checking with lumped mass
		
		for (int i=0;i<4;i++)
		{ 
		for (int k=0;k<4;k++)
			{
			if (i==k)
				Mass(i,k)=Mass(i,0)+Mass(i,1)+Mass(i,2)+Mass(i,3);
			else 
				Mass(i,k)=0.0;
			}
		}
		
		for(unsigned int i=0;i<4;i++)
		{
			for (unsigned int k=0;k<4;k++)
			{
			rLeftHandSideMatrix(3*i,3*k)=Mass(i,k);
			rLeftHandSideMatrix(3*i+1,3*k+1)=Mass(i,k);
			rLeftHandSideMatrix(3*i+2,3*k+2)=Mass(i,k);
			}
				
		}
		KRATOS_WATCH(rLeftHandSideMatrix)

		//and now we  assemble the RHS, which is computed using 2 Gauss points - we integrate right at the intersection points
		//so: rhs=0.5 * l *( (N1 N2 N3)_at_point1 * (v_prescribed1) + (N1 N2 N3)_at_point2 * (v_prescribed2))  , where l is the length of the interface		
		//first point contribution
		
		CalculateN_at_Point(GetGeometry(), mPoint1[0], mPoint1[1], mPoint1[2], msN);

		//first all the x_components
		double x_contr=mVel1[0]-(msN[0]*Vel_VEC[0]+msN[1]*Vel_VEC[3]+msN[2]*Vel_VEC[6]+msN[3]*Vel_VEC[9]);
		double y_contr=mVel1[1]-(msN[0]*Vel_VEC[1]+msN[1]*Vel_VEC[4]+msN[2]*Vel_VEC[7]+msN[3]*Vel_VEC[10]);
		double z_contr=mVel1[2]-(msN[0]*Vel_VEC[2]+msN[1]*Vel_VEC[5]+msN[2]*Vel_VEC[8]+msN[3]*Vel_VEC[11]);

		
		rRightHandSideVector[0]=msN[0]*x_contr;
		rRightHandSideVector[1]=msN[0]*y_contr;		
		rRightHandSideVector[2]=msN[0]*z_contr;

		rRightHandSideVector[3]=msN[1]*x_contr;
		rRightHandSideVector[4]=msN[1]*y_contr;
		rRightHandSideVector[5]=msN[1]*z_contr;

		rRightHandSideVector[6]=msN[2]*x_contr;
		rRightHandSideVector[7]=msN[2]*y_contr;
		rRightHandSideVector[8]=msN[2]*z_contr;

		rRightHandSideVector[9]=msN[3]*x_contr;
		rRightHandSideVector[10]=msN[3]*y_contr;
		rRightHandSideVector[11]=msN[3]*z_contr;

		KRATOS_WATCH(x_contr)
		KRATOS_WATCH(y_contr)
		KRATOS_WATCH(z_contr)

		



		////////////////////////////////////////////////////////////////////////////////////////////////////
		//second point contribution
		CalculateN_at_Point(GetGeometry(), mPoint2[0], mPoint2[1], mPoint2[2], msN);
		//KRATOS_WATCH(msN)
		//first all the x_components
		//first all the x_components
		x_contr=mVel2[0]-(msN[0]*Vel_VEC[0]+msN[1]*Vel_VEC[3]+msN[2]*Vel_VEC[6]+msN[3]*Vel_VEC[9]);
		y_contr=mVel2[1]-(msN[0]*Vel_VEC[1]+msN[1]*Vel_VEC[4]+msN[2]*Vel_VEC[7]+msN[3]*Vel_VEC[10]);
		z_contr=mVel2[2]-(msN[0]*Vel_VEC[2]+msN[1]*Vel_VEC[5]+msN[2]*Vel_VEC[8]+msN[3]*Vel_VEC[11]);
		
		rRightHandSideVector[0]+=msN[0]*x_contr;
		rRightHandSideVector[1]+=msN[0]*y_contr;		
		rRightHandSideVector[2]+=msN[0]*z_contr;

		rRightHandSideVector[3]+=msN[1]*x_contr;
		rRightHandSideVector[4]+=msN[1]*y_contr;
		rRightHandSideVector[5]+=msN[1]*z_contr;

		rRightHandSideVector[6]+=msN[2]*x_contr;
		rRightHandSideVector[7]+=msN[2]*y_contr;
		rRightHandSideVector[8]+=msN[2]*z_contr;

		rRightHandSideVector[9]+=msN[3]*x_contr;
		rRightHandSideVector[10]+=msN[3]*y_contr;
		rRightHandSideVector[11]+=msN[3]*z_contr;

		KRATOS_WATCH(x_contr)
		KRATOS_WATCH(y_contr)
		KRATOS_WATCH(z_contr)


		////////////////////////////////////////////////////////////////////////////////////////////////////////
		//third point contribution
		CalculateN_at_Point(GetGeometry(), mPoint3[0], mPoint3[1], mPoint3[2], msN);
		//KRATOS_WATCH(msN)
		//first all the x_components
		//first all the x_components
		x_contr=mVel3[0]-(msN[0]*Vel_VEC[0]+msN[1]*Vel_VEC[3]+msN[2]*Vel_VEC[6]+msN[3]*Vel_VEC[9]);
		y_contr=mVel3[1]-(msN[0]*Vel_VEC[1]+msN[1]*Vel_VEC[4]+msN[2]*Vel_VEC[7]+msN[3]*Vel_VEC[10]);
		z_contr=mVel3[2]-(msN[0]*Vel_VEC[2]+msN[1]*Vel_VEC[5]+msN[2]*Vel_VEC[8]+msN[3]*Vel_VEC[11]);	

		rRightHandSideVector[0]+=msN[0]*x_contr;
		rRightHandSideVector[1]+=msN[0]*y_contr;		
		rRightHandSideVector[2]+=msN[0]*z_contr;

		rRightHandSideVector[3]+=msN[1]*x_contr;
		rRightHandSideVector[4]+=msN[1]*y_contr;
		rRightHandSideVector[5]+=msN[1]*z_contr;

		rRightHandSideVector[6]+=msN[2]*x_contr;
		rRightHandSideVector[7]+=msN[2]*y_contr;
		rRightHandSideVector[8]+=msN[2]*z_contr;

		rRightHandSideVector[9]+=msN[3]*x_contr;
		rRightHandSideVector[10]+=msN[3]*y_contr;
		rRightHandSideVector[11]+=msN[3]*z_contr;

		KRATOS_WATCH(x_contr)
		KRATOS_WATCH(y_contr)
		KRATOS_WATCH(z_contr)
		

		//completing integration using two Gauss points
		///UTOCHNIT' NASCHET KONSTANY!!!!!!!!!!!!!
		rRightHandSideVector*=0.333333333333*inters_area;
		KRATOS_WATCH(inters_area)
		

		
		if (this->mNumber_of_intersections==4)
			{
			KRATOS_WATCH("SECond triangle")
			KRATOS_WATCH(inters_area2)
			array_1d<double,12> RHS_TEMP=ZeroVector(12);
			//first point contribution
			CalculateN_at_Point(GetGeometry(), PointsOfSecondTriangle[0][0], PointsOfSecondTriangle[0][1], PointsOfSecondTriangle[0][2], msN);

			//first all the x_components
			x_contr=VelsOfSecondTriangle[0][0]-(msN[0]*Vel_VEC[0]+msN[1]*Vel_VEC[3]+msN[2]*Vel_VEC[6]+msN[3]*Vel_VEC[9]);
			y_contr=VelsOfSecondTriangle[0][1]-(msN[0]*Vel_VEC[1]+msN[1]*Vel_VEC[4]+msN[2]*Vel_VEC[7]+msN[3]*Vel_VEC[10]);
			z_contr=VelsOfSecondTriangle[0][2]-(msN[0]*Vel_VEC[2]+msN[1]*Vel_VEC[5]+msN[2]*Vel_VEC[8]+msN[3]*Vel_VEC[11]);

			RHS_TEMP[0]=msN[0]*x_contr;
			RHS_TEMP[1]=msN[0]*y_contr;
			RHS_TEMP[2]=msN[0]*z_contr;

			RHS_TEMP[3]=msN[1]*x_contr;
			RHS_TEMP[4]=msN[1]*y_contr;
			RHS_TEMP[5]=msN[1]*z_contr;

			RHS_TEMP[6]=msN[2]*x_contr;
			RHS_TEMP[7]=msN[2]*y_contr;
			RHS_TEMP[8]=msN[2]*z_contr;

			RHS_TEMP[9]=msN[3]*x_contr;
			RHS_TEMP[10]=msN[3]*y_contr;
			RHS_TEMP[11]=msN[3]*z_contr;

			KRATOS_WATCH(x_contr)
			KRATOS_WATCH(y_contr)
			KRATOS_WATCH(z_contr)

			////////////////////////////////////////////////////////////////////////////////////////////////////
			//second point contribution
			CalculateN_at_Point(GetGeometry(), PointsOfSecondTriangle[1][0], PointsOfSecondTriangle[1][1], PointsOfSecondTriangle[1][2], msN);
			//KRATOS_WATCH(msN)
			//first all the x_components
			//first all the x_components
			x_contr=VelsOfSecondTriangle[1][0]-(msN[0]*Vel_VEC[0]+msN[1]*Vel_VEC[3]+msN[2]*Vel_VEC[6]+msN[3]*Vel_VEC[9]);
			y_contr=VelsOfSecondTriangle[1][1]-(msN[0]*Vel_VEC[1]+msN[1]*Vel_VEC[4]+msN[2]*Vel_VEC[7]+msN[3]*Vel_VEC[10]);
			z_contr=VelsOfSecondTriangle[1][2]-(msN[0]*Vel_VEC[2]+msN[1]*Vel_VEC[5]+msN[2]*Vel_VEC[8]+msN[3]*Vel_VEC[11]);


			RHS_TEMP[0]+=msN[0]*x_contr;
			RHS_TEMP[1]+=msN[0]*y_contr;
			RHS_TEMP[2]+=msN[0]*z_contr;

			RHS_TEMP[3]+=msN[1]*x_contr;
			RHS_TEMP[4]+=msN[1]*y_contr;
			RHS_TEMP[5]+=msN[1]*z_contr;

			RHS_TEMP[6]+=msN[2]*x_contr;
			RHS_TEMP[7]+=msN[2]*y_contr;
			RHS_TEMP[8]+=msN[2]*z_contr;

			RHS_TEMP[9]+=msN[3]*x_contr;
			RHS_TEMP[10]+=msN[3]*y_contr;
			RHS_TEMP[11]+=msN[3]*z_contr;

			KRATOS_WATCH(x_contr)
			KRATOS_WATCH(y_contr)
			KRATOS_WATCH(z_contr)
			////////////////////////////////////////////////////////////////////////////////////////////////////////
			//third point contribution
			CalculateN_at_Point(GetGeometry(), PointsOfSecondTriangle[2][0], PointsOfSecondTriangle[2][1], PointsOfSecondTriangle[2][2], msN);

			KRATOS_WATCH(x_contr)
			KRATOS_WATCH(y_contr)
			KRATOS_WATCH(z_contr)
			//KRATOS_WATCH(msN)
			//first all the x_components
			//first all the x_components
			x_contr=VelsOfSecondTriangle[2][0]-(msN[0]*Vel_VEC[0]+msN[1]*Vel_VEC[3]+msN[2]*Vel_VEC[6]+msN[3]*Vel_VEC[9]);
			y_contr=VelsOfSecondTriangle[2][1]-(msN[0]*Vel_VEC[1]+msN[1]*Vel_VEC[4]+msN[2]*Vel_VEC[7]+msN[3]*Vel_VEC[10]);
			z_contr=VelsOfSecondTriangle[2][2]-(msN[0]*Vel_VEC[2]+msN[1]*Vel_VEC[5]+msN[2]*Vel_VEC[8]+msN[3]*Vel_VEC[11]);
			
			RHS_TEMP[0]+=msN[0]*x_contr;
			RHS_TEMP[1]+=msN[0]*y_contr;
			RHS_TEMP[2]+=msN[0]*z_contr;

			RHS_TEMP[3]+=msN[1]*x_contr;
			RHS_TEMP[4]+=msN[1]*y_contr;
			RHS_TEMP[5]+=msN[1]*z_contr;

			RHS_TEMP[6]+=msN[2]*x_contr;
			RHS_TEMP[7]+=msN[2]*y_contr;
			RHS_TEMP[8]+=msN[2]*z_contr;

			RHS_TEMP[9]+=msN[3]*x_contr;
			RHS_TEMP[10]+=msN[3]*y_contr;
			RHS_TEMP[11]+=msN[3]*z_contr;

			KRATOS_WATCH(x_contr)
			
			RHS_TEMP*=0.333333333333*inters_area2;
			KRATOS_WATCH(rRightHandSideVector)
			KRATOS_WATCH(RHS_TEMP)

			rRightHandSideVector+=RHS_TEMP;
		
			
			}
		
		KRATOS_WATCH(rRightHandSideVector)

		
		
		
	}
	/////////////////////////////////////////////////////////////
	double ProjDirichletCond3D::Length(array_1d<double,3>& Point1, array_1d<double,3>& Point2)
	{
	KRATOS_WATCH("length calculation")
	return sqrt((Point1[0]-Point2[0])*(Point1[0]-Point2[0]) + (Point1[1]-Point2[1])*(Point1[1]-Point2[1]) +(Point1[2]-Point2[2])*(Point1[2]-Point2[2]));
	}
	////////////////////////////////////////////////////////////
	double ProjDirichletCond3D::CalculateTriangleArea3D(	array_1d<double,3>& Point1, array_1d<double,3>& Point2, array_1d<double,3>& Point3	)
		{
			//Heron's formula
			double a=Length(Point1, Point2);//sqrt((Point1[0]-Point2[0])*(Point1[0]-Point2[0]) + (Point1[1]-Point2[1])*(Point1[1]-Point2[1]) +(Point1[2]-Point2[2])*(Point1[2]-Point2[2]));
			double b=Length(Point1, Point3);//sqrt((Point3[0]-Point2[0])*(Point3[0]-Point2[0]) + (Point3[1]-Point2[1])*(Point3[1]-Point2[1]) +(Point3[2]-Point2[2])*(Point3[2]-Point2[2]));
			double c=Length(Point2, Point3);//sqrt((Point1[0]-Point3[0])*(Point1[0]-Point3[0]) + (Point1[1]-Point3[1])*(Point1[1]-Point3[1]) +(Point1[2]-Point3[2])*(Point1[2]-Point3[2]));
			double p=0.5*(a+b+c);
			return sqrt(p*(p-a)*(p-b)*(p-c));
		}

	double ProjDirichletCond3D::CalculateVol(	const double x0, const double y0, const double z0,
					const double x1, const double y1, const double z1,
    					const double x2, const double y2, const double z2,
    					const double x3, const double y3, const double z3
					  )
	{
			double x10 = x1 - x0;
			double y10 = y1 - y0;
			double z10 = z1 - z0;

			double x20 = x2 - x0;
			double y20 = y2 - y0;
			double z20 = z2 - z0;

			double x30 = x3 - x0;
			double y30 = y3 - y0;
			double z30 = z3 - z0;

			double detJ = (x10 * y20 * z30 - x10 * y30 * z20) + (y10 * z20 * x30 - y10 * x20 * z30) + (z10 * x20 * y30 - z10 * y20 * x30);
			return  detJ*0.1666666666666666666667;			
			
	}
	//************************************************************************************
	//************************************************************************************
	void ProjDirichletCond3D::CalculateN_at_Point(Element::GeometryType& geom, const double xc, const double yc, const double zc, array_1d<double,4>& N_at_c)
	{
			double x0=geom[0].X();
			double x1=geom[1].X();		
			double x2=geom[2].X();
			double x3=geom[3].X();		

			double y0=geom[0].Y();
			double y1=geom[1].Y();		
			double y2=geom[2].Y();
			double y3=geom[3].Y();		

			double z0=geom[0].Z();
			double z1=geom[1].Z();		
			double z2=geom[2].Z();
			double z3=geom[3].Z();		

			double vol = CalculateVol(x0,y0,z0,x1,y1,z1,x2,y2,z2,x3,y3,z3);

			double inv_vol = 0.0;
			if(vol < 0.000000000000001)
			  {
				KRATOS_ERROR(std::logic_error,"element with zero vol found","");
			  }
			else
			  {
				inv_vol = 1.0 / vol;
			  }

			N_at_c[0] = CalculateVol(x1,y1,z1,x3,y3,z3,x2,y2,z2,xc,yc,zc) * inv_vol;
			N_at_c[1] = CalculateVol(x3,y3,z3,x0,y0,z0,x2,y2,z2,xc,yc,zc) * inv_vol;
			N_at_c[2] = CalculateVol(x3,y3,z3,x1,y1,z1,x0,y0,z0,xc,yc,zc) * inv_vol;
	   		N_at_c[3] = CalculateVol(x0,y0,z0,x1,y1,z1,x2,y2,z2,xc,yc,zc) * inv_vol;

			//if the intersection is close one of the tetrahedra's vertices ->send this message
			if (  (N_at_c[0]<0.05 && N_at_c[1]<0.05 && N_at_c[2]<0.05) || 
			      (N_at_c[0]<0.05 && N_at_c[2]<0.05 && N_at_c[3]<0.05) || 
			      (N_at_c[2]<0.05 && N_at_c[1]<0.05 && N_at_c[3]<0.05) || 
			      (N_at_c[0]<0.05 && N_at_c[1]<0.05 && N_at_c[3]<0.05)     )			
				KRATOS_WATCH("Dangerous VERTICES!!! Intersection is very close to the node")		
			//KRATOS_ERROR(std::logic_error,  "Too close to the node is the INTERSECTION!!!! " , "")

	}


	//************************************************************************************
	//************************************************************************************
	void ProjDirichletCond3D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
	{
		if(rResult.size() != 12)
			rResult.resize(12);

		rResult[0] = (GetGeometry()[0].GetDof(AUX_VEL_X)).EquationId();
		rResult[1] = (GetGeometry()[0].GetDof(AUX_VEL_Y)).EquationId();
		rResult[2] = (GetGeometry()[0].GetDof(AUX_VEL_Z)).EquationId();

		rResult[3] = (GetGeometry()[1].GetDof(AUX_VEL_X)).EquationId();
		rResult[4] = (GetGeometry()[1].GetDof(AUX_VEL_Y)).EquationId();
		rResult[5] = (GetGeometry()[1].GetDof(AUX_VEL_Z)).EquationId();

		rResult[6] = (GetGeometry()[2].GetDof(AUX_VEL_X)).EquationId();
		rResult[7] = (GetGeometry()[2].GetDof(AUX_VEL_Y)).EquationId();
		rResult[8] = (GetGeometry()[2].GetDof(AUX_VEL_Z)).EquationId();

		rResult[9] = (GetGeometry()[3].GetDof(AUX_VEL_X)).EquationId();
		rResult[10] = (GetGeometry()[3].GetDof(AUX_VEL_Y)).EquationId();
		rResult[11] = (GetGeometry()[3].GetDof(AUX_VEL_Z)).EquationId();


	}

	//************************************************************************************
	//************************************************************************************
	void ProjDirichletCond3D::GetDofList(DofsVectorType& ConditionalDofList,ProcessInfo& CurrentProcessInfo)
	{
		if(ConditionalDofList.size() != 12)
			ConditionalDofList.resize(12);

		ConditionalDofList[0] = (GetGeometry()[0].pGetDof(AUX_VEL_X));
		ConditionalDofList[1] = (GetGeometry()[0].pGetDof(AUX_VEL_Y));
		ConditionalDofList[2] = (GetGeometry()[0].pGetDof(AUX_VEL_Z));
		
		ConditionalDofList[3] = (GetGeometry()[1].pGetDof(AUX_VEL_X));
		ConditionalDofList[4] = (GetGeometry()[1].pGetDof(AUX_VEL_Y));
		ConditionalDofList[5] = (GetGeometry()[1].pGetDof(AUX_VEL_Z));

		ConditionalDofList[6] = (GetGeometry()[2].pGetDof(AUX_VEL_X));
		ConditionalDofList[7] = (GetGeometry()[2].pGetDof(AUX_VEL_Y));
		ConditionalDofList[8] = (GetGeometry()[2].pGetDof(AUX_VEL_Z));

		ConditionalDofList[9] = (GetGeometry()[3].pGetDof(AUX_VEL_X));
		ConditionalDofList[10] = (GetGeometry()[3].pGetDof(AUX_VEL_Y));
		ConditionalDofList[11] = (GetGeometry()[3].pGetDof(AUX_VEL_Z));


	}
	
} // Namespace Kratos


