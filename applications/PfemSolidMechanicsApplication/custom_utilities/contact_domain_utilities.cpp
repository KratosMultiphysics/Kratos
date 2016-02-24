//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                    July 2013 $
//   Revision:            $Revision:                      0.0 $
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

	void ContactDomainUtilities::CalculateBaseDistances (BaseLengths& Base,
							     LocalVectorType& P1,
							     LocalVectorType& P2,
							     LocalVectorType& PS,
							     LocalVectorType& Normal)
	{

		Base.L=norm_2(P2-P1);


		//if the normal points from the side to the slave node:
		LocalVectorType Projection= PS-P1;
		Projection-=Normal*(inner_prod(Projection,Normal));
		Projection+=P1;

		double sign=1;
		LocalVectorType Pro1 = Projection-P1; 
		LocalVectorType Pro2 = P2-P1;

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

	ContactDomainUtilities::LocalVectorType & ContactDomainUtilities::CalculateFaceNormal(LocalVectorType &Normal, LocalVectorType& P1, LocalVectorType &P2)
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


	ContactDomainUtilities::LocalVectorType & ContactDomainUtilities::CalculateFaceTangent(LocalVectorType &Tangent ,LocalVectorType& P1, LocalVectorType &P2)
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


	ContactDomainUtilities::LocalVectorType & ContactDomainUtilities::CalculateFaceTangent(LocalVectorType &Tangent ,LocalVectorType& Normal)
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

