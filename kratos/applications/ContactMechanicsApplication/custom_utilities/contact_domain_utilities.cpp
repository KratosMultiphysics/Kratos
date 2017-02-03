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


	//*******************************************************************************************
	//*******************************************************************************************

        void ContactDomainUtilities::CalculateBaseDistances (std::vector<BaseLengths>& BaseVector,
							     PointType& P1,
							     PointType& P2,
							     PointType& P3,
							     PointType& PS,
							     PointType& Normal)
	{
	  // uint_mat lpofa(4,4);
	  // Getlpofa(lpofa);
  
	  // //std::cout<<" lpofa "<<lpofa(0,slave)<<" "<<lpofa(1,slave)<<" "<<lpofa(2,slave)<<" "<<slave<<std::endl;

	  // Base[0].L=VertexPosition[lpofa(0,slave)].DistancePoint(VertexPosition[lpofa(1,slave)]); 
	  // Base[1].L=VertexPosition[lpofa(1,slave)].DistancePoint(VertexPosition[lpofa(2,slave)]);  
	  // Base[2].L=VertexPosition[lpofa(2,slave)].DistancePoint(VertexPosition[lpofa(0,slave)]);  

	  // //because the normal points from the side to the slave node:
	  // Pt Projection=VertexPosition[slave].ProjectionPlane(VertexPosition[lpofa(0,slave)],VertexPosition[lpofa(1,slave)],VertexPosition[lpofa(2,slave)]);

	  // Pt BaseLine  = (VertexPosition[lpofa(1,slave)]-VertexPosition[lpofa(0,slave)]).direction(); //V1
	  // Pt PlaneLine = (VertexPosition[lpofa(1,slave)]-VertexPosition[lpofa(2,slave)]).direction(); //V2
  
	  // Pt Projected = Projection-VertexPosition[lpofa(0,slave)]; //P2-P1

	  // //intersection of lines:
	  // Pt PPlane=(Projected%PlaneLine); //(P2-P1)xV2
	  // Pt BPlane=(BaseLine%PlaneLine); //V1xV2
  
	  // double sign=PPlane.direction()*BPlane.direction();
    
  
	  // if(BPlane.modulus()>0){
	  //   Base[0].B =(PPlane.modulus()/BPlane.modulus())*sign;
	  // }
	  // else{
	  //   Base[0].B =0;
	  // }


	  // BaseLine*= (-1); //V1
	  // Projected=Projection-VertexPosition[lpofa(1,slave)]; //P2-P1


	  // PPlane=(Projected%PlaneLine);
	  // BPlane=(BaseLine%PlaneLine);
  
	  // sign=PPlane.direction()*BPlane.direction();
    
  
	  // if(BPlane.modulus()>0){
	  //   Base[0].A =(PPlane.modulus()/BPlane.modulus())*sign;
	  // }
	  // else{
	  //   Base[0].A =0;
	  // }


	  // BaseLine  = (VertexPosition[lpofa(2,slave)]-VertexPosition[lpofa(1,slave)]).direction();
	  // PlaneLine = (VertexPosition[lpofa(2,slave)]-VertexPosition[lpofa(0,slave)]).direction();
	  // Projected = Projection-VertexPosition[lpofa(1,slave)];

	  // PPlane=(Projected%PlaneLine);
	  // BPlane=(BaseLine%PlaneLine);
  
	  // sign=PPlane.direction()*BPlane.direction();
    
  
	  // if(BPlane.modulus()>0){
	  //   Base[1].B =(PPlane.modulus()/BPlane.modulus())*sign;
	  // }
	  // else{
	  //   Base[1].B =0;
	  // }


	  // BaseLine *= (-1);
	  // Projected = Projection-VertexPosition[lpofa(2,slave)];

	  // PPlane=(Projected%PlaneLine);
	  // BPlane=(BaseLine%PlaneLine);
  
	  // sign=PPlane.direction()*BPlane.direction();
    
  
	  // if(BPlane.modulus()>0){
	  //   Base[1].A =(PPlane.modulus()/BPlane.modulus())*sign;
	  // }
	  // else{
	  //   Base[1].A =0;
	  // }


	  // BaseLine  = (VertexPosition[lpofa(0,slave)]-VertexPosition[lpofa(2,slave)]).direction();
	  // PlaneLine = (VertexPosition[lpofa(0,slave)]-VertexPosition[lpofa(1,slave)]).direction();
	  // Projected = Projection-VertexPosition[lpofa(2,slave)];

	  // PPlane=(Projected%PlaneLine);
	  // BPlane=(BaseLine%PlaneLine);
  
	  // sign=PPlane.direction()*BPlane.direction();
    
  
	  // if(BPlane.modulus()>0){
	  //   Base[2].B =(PPlane.modulus()/BPlane.modulus())*sign;
	  // }
	  // else{
	  //   Base[2].B =0;
	  // }


	  // BaseLine *= (-1);
	  // Projected = Projection-VertexPosition[lpofa(0,slave)];

	  // PPlane=(Projected%PlaneLine);
	  // BPlane=(BaseLine%PlaneLine);
  
	  // sign=PPlane.direction()*BPlane.direction();
    
  
	  // if(BPlane.modulus()>0){
	  //   Base[2].A =(PPlane.modulus()/BPlane.modulus())*sign;
	  // }
	  // else{
	  //   Base[2].A =0;
	  // }

	  //   std::cout<<" Base 0-> L: "<<Base[0].L<<" A: "<<Base[0].A<<" B: "<<Base[0].B<<std::endl;
	  //   std::cout<<" Base 1-> L: "<<Base[1].L<<" A: "<<Base[1].A<<" B: "<<Base[1].B<<std::endl;
	  //   std::cout<<" Base 2-> L: "<<Base[2].L<<" A: "<<Base[2].A<<" B: "<<Base[2].B<<std::endl;

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


		//if the normal points from the side to the slave node:
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

		Normal.clear();
		Normal[0] =    P2[1] - P1[1];
		Normal[1] = - (P2[0] - P1[0]);
		Normal[2] =    0.00;

		if(norm_2(Normal)!=0)
			Normal/=norm_2(Normal);

		return Normal;
	}

        //************************************************************************************
        //************************************************************************************


	ContactDomainUtilities::PointType & ContactDomainUtilities::CalculateFaceTangent(PointType &Tangent ,PointType& P1, PointType &P2)
	{

		Tangent.clear();
		Tangent[0] =    (P2[0] - P1[0]);
		Tangent[1] =   -(P2[1] - P1[1]);
		Tangent[2] =    0.00;

		if(norm_2(Tangent)!=0)
			Tangent/=norm_2(Tangent);

		return Tangent;

	}



       //************************************************************************************
       //************************************************************************************


	ContactDomainUtilities::PointType & ContactDomainUtilities::CalculateFaceTangent(PointType &Tangent ,PointType& Normal)
	{

		Tangent.clear();
		Tangent[0] =  - Normal[1];
		Tangent[1] =    Normal[0];
		Tangent[2] =    0.00;

		if(norm_2(Tangent)!=0)
			Tangent/=norm_2(Tangent);

		return Tangent;

	}



} // Namespace Kratos

