//
//   Project Name:        KratosContactMechanicsApplication $
//   Created by:          $Author:              JMCarbonell $
//   Last modified by:    $Co-Author:                       $
//   Date:                $Date:                  July 2016 $
//   Revision:            $Revision:                    0.0 $
//
//

// System includes

// External includes

// Project includes
#include "includes/kratos_flags.h"
#include "custom_utilities/contact_domain_utilities.hpp"

namespace Kratos
{

	KRATOS_CREATE_LOCAL_FLAG( ContactDomainUtilities, COMPUTE_RHS_VECTOR,                 0 );
	KRATOS_CREATE_LOCAL_FLAG( ContactDomainUtilities, COMPUTE_LHS_MATRIX,                 1 );

        KRATOS_CREATE_LOCAL_FLAG( ContactDomainUtilities, COMPUTE_RHS_VECTOR_WITH_COMPONENTS, 2 );
        KRATOS_CREATE_LOCAL_FLAG( ContactDomainUtilities, COMPUTE_LHS_MATRIX_WITH_COMPONENTS, 3 );

	KRATOS_CREATE_LOCAL_FLAG( ContactDomainUtilities, COMPUTE_FRICTION_FORCES,            4 );
	KRATOS_CREATE_LOCAL_FLAG( ContactDomainUtilities, COMPUTE_FRICTION_STIFFNESS,         5 );

	KRATOS_CREATE_LOCAL_FLAG( ContactDomainUtilities, COMPUTE_NODAL_CONTACT_FORCES,       6 );
       
	//*******************************************************************************************
	//*******************************************************************************************

        //GEOMETRICAL THEORY
        // Heron's Formula for the area of a triangle(Hero's Formula)
        // A method for calculating the area of a triangle when you know the lengths of all three sides.
        // Let a,b,c be the lengths of the sides of a triangle. The area is given by:

        // area = sqrt( p*(p−a)*(p−b)*(p−c) )   

        // where p is half the perimeter, or
  
        // p = 0.5 *(a+b+c)


        void ContactDomainUtilities::CalculateBaseArea (double& area,
							double& a,
							double& b,
							double& c)
	{

	  double p = 0.5*(a+b+c);
	  area = sqrt( p*(p-a)*(p-b)*(p-c) );
	    
	}

  	//*******************************************************************************************
	//*******************************************************************************************

        void ContactDomainUtilities::CalculateEdgeDistances (std::vector<BaseLengths>& BaseVector,
							     PointType& P1,
							     PointType& P2,
							     PointType& PS1,
							     PointType& PS2,
							     PointType& Normal)
	{

	  BaseVector[0].L=norm_2(P2-P1);
	  BaseVector[1].L=norm_2(PS2-PS1);

	  //the normal points from the side to the slave node:  
	  PointType V1; 
	  PointType V2;
	  
	  //projection of the slave on the master plane:  
	  PointType M1 = 0.5*(P2-P1);
	  PointType M2 = 0.5*(PS2-PS1);
   
	  //Comtutacion of the line bases:

	  //BaseVector[0]:
	  
	  V1 = P2-P1;
	  V2 = Normal;

	  V1 /= norm_2(V1);
	  V2 /= norm_2(V2);
	  
	  CalculateLineIntersection(BaseVector[0].B, P1, M2, V1, V2 );

	  BaseVector[0].A = BaseVector[0].L-BaseVector[0].B;

	  //BaseVector[1]:
	  
	  V1 = PS2-PS1;
	  V2 = -Normal;

	  V1 /= norm_2(V1);
	  V2 /= norm_2(V2);
	  
	  CalculateLineIntersection(BaseVector[1].B, PS1, M1 , V1, V2 );

	  BaseVector[1].A = BaseVector[1].L-BaseVector[1].B;


	  std::cout<<" BaseVector 0-> L: "<<BaseVector[0].L<<" A: "<<BaseVector[0].A<<" B: "<<BaseVector[0].B<<std::endl;
	  std::cout<<" BaseVector 1-> L: "<<BaseVector[1].L<<" A: "<<BaseVector[1].A<<" B: "<<BaseVector[1].B<<std::endl;
	

	}
  
  	//*******************************************************************************************
	//*******************************************************************************************

        void ContactDomainUtilities::CalculateBaseDistances (std::vector<BaseLengths>& BaseVector,
							     PointType& P1,
							     PointType& P2,
							     PointType& P3,
							     PointType& PS,
							     PointType& Normal)
	{

	  BaseVector[0].L=norm_2(P2-P1);
	  BaseVector[1].L=norm_2(P3-P2);
	  BaseVector[2].L=norm_2(P1-P3);
	  
	  //the normal points from the side to the slave node:  
	  PointType V1; 
	  PointType V2;
	  
	  // V1 = P2-P1;
	  // V2 = P3-P1;	  
	  // MathUtils<double>::CrossProduct(Normal,V1,V2);	  
	  // if( norm_2(Normal) != 0 )
	  //   Normal/=norm_2(Normal);	  

	  //projection of the slave on the master plane:  
	  PointType PPS = PS-P1;
	  PPS-=Normal*(inner_prod(PPS,Normal));
	  PPS+=P1;
    
	  //Comtutacion of the line bases:

	  //BaseVector[0]:
	  
	  V1 = P2-P1;
	  V2 = P3-P2;

	  V1 /= norm_2(V1);
	  V2 /= norm_2(V2);
	  
	  CalculateLineIntersection(BaseVector[0].B, P1, PPS, V1, V2 );

	  BaseVector[0].A = BaseVector[0].L-BaseVector[0].B;

	  //BaseVector[1]:
	  
	  V1 = P3-P2;
	  V2 = P1-P3;

	  V1 /= norm_2(V1);
	  V2 /= norm_2(V2);
	  
	  CalculateLineIntersection(BaseVector[1].B, P2, PPS, V1, V2 );

	  BaseVector[1].A = BaseVector[1].L-BaseVector[1].B;

	  //BaseVector[2]:
	  
	  V1 = P1-P3;
	  V2 = P2-P1;

	  V1 /= norm_2(V1);
	  V2 /= norm_2(V2);
	  
	  CalculateLineIntersection(BaseVector[1].B, P3, PPS, V1, V2 );

	  BaseVector[2].A = BaseVector[2].L-BaseVector[2].B;

	  std::cout<<" BaseVector 0-> L: "<<BaseVector[0].L<<" A: "<<BaseVector[0].A<<" B: "<<BaseVector[0].B<<std::endl;
	  std::cout<<" BaseVector 1-> L: "<<BaseVector[1].L<<" A: "<<BaseVector[1].A<<" B: "<<BaseVector[1].B<<std::endl;
	  std::cout<<" BaseVector 2-> L: "<<BaseVector[2].L<<" A: "<<BaseVector[2].A<<" B: "<<BaseVector[2].B<<std::endl;

	}

	//*******************************************************************************************
	//*******************************************************************************************
  
        //GEOMETRICAL THEORY
        // The method was to
        // use the two equations that represent the lines:

        //   L1 = P1 + aV1
        //   L2 = P2 + bV2

        // Where V1 and V2 are the director vectors
        // P1 point which belongs to L1
        // P2 point which belongs to L2
        
        // => P1 + aV1 = P2 + bV2 
        // => aV1 = (P2-P1) + bV2

        // applying the vectorial (x V2) in the two sides of the equation:

        // a(V1 x V2) = (P2-P1) x V2

        // once we have "a" we can replace it in the 1st equation to get the 
        // point of intersection

        void ContactDomainUtilities::CalculateLineIntersection (double& a,
								const PointType& P1,
								const PointType& P2,
								const PointType& V1,
								const PointType& V2)
	{
	  PointType P2_P1 = P2-P1;
	  PointType V1xV2;
	  MathUtils<double>::CrossProduct(V1xV2,V1,V2);

	  PointType P2_P1xV2;
	  MathUtils<double>::CrossProduct(P2_P1xV2,P2_P1,V2);

	  if( norm_2(V1xV2) != 0 )   
	    a = norm_2(P2_P1xV2)/norm_2(V1xV2);
	  else
	    a = 0;	  
	}

	//*******************************************************************************************
	//*******************************************************************************************

	void ContactDomainUtilities::CalculateBaseDistances (BaseLengths& Base,
							     PointType& P1,
							     PointType& P2,
							     PointType& PS,
							     PointType& Normal)
	{

		Base.L=norm_2(P2-P1);

		//the normal points from the side to the slave node:

		//projection of the slave on the master plane:  
		PointType Projection= PS-P1;
		Projection-=Normal*(inner_prod(Projection,Normal));
		Projection+=P1;

		double sign=1;
		PointType Pro1 = Projection-P1; 
		PointType Pro2 = P2-P1;

		if(double(inner_prod(Pro2,Pro1))<0)
			sign*=(-1);

		//signed distance to node 1
		Base.B= sign*norm_2(P1-Projection);

		sign=1;
		Pro1 = Projection-P2;
		Pro2 = P1-P2;
		if(inner_prod(Pro2,Pro1)<0)
			sign*=(-1);

		//signed distance to node 2
		Base.A= sign*norm_2(P2-Projection);

	}



        //************************************************************************************
        //************************************************************************************

	ContactDomainUtilities::PointType & ContactDomainUtilities::CalculateSurfaceNormal(PointType &Normal, PointType& D1, PointType& D2)
	{

		Normal.clear();

		MathUtils<double>::CrossProduct(Normal, D1, D2);
		
		if(norm_2(Normal)!=0)
			Normal/=norm_2(Normal);

		return Normal;
	}
  
        //************************************************************************************
        //************************************************************************************

	ContactDomainUtilities::PointType & ContactDomainUtilities::CalculateFaceNormal(PointType &Normal, PointType& P1, PointType &P2)
	{
   	        //contact element is ouside (sign criterion)
		Normal.clear();
		Normal[0] = - (P2[1] - P1[1]);
		Normal[1] =   (P2[0] - P1[0]);
		Normal[2] =    0.00;

		if(norm_2(Normal)!=0)
		  Normal/=norm_2(Normal);
		
		return Normal;
	}

        //************************************************************************************
        //************************************************************************************


	ContactDomainUtilities::PointType & ContactDomainUtilities::CalculateFaceTangent(PointType &Tangent ,PointType& P1, PointType &P2)
	{
   	        //contact element is ouside (sign criterion)
		Tangent.clear();
		Tangent[0] =  -(P2[0] - P1[0]);
		Tangent[1] =  -(P2[1] - P1[1]);
		Tangent[2] =    0.00;

		if(norm_2(Tangent)!=0)
			Tangent/=norm_2(Tangent);

		return Tangent;

	}



       //************************************************************************************
       //************************************************************************************


	ContactDomainUtilities::PointType & ContactDomainUtilities::CalculateFaceTangent(PointType &Tangent ,PointType& Normal)
	{
	        //counter clock-wise movement
		Tangent.clear();
		Tangent[0] =  - Normal[1];
		Tangent[1] =    Normal[0];
		Tangent[2] =    0.00;

		if(norm_2(Tangent)!=0)
			Tangent/=norm_2(Tangent);

		return Tangent;

	}



} // Namespace Kratos

