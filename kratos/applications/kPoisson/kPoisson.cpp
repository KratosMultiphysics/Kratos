//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: it's me $
//   Date:                $Date: 2008-08-08 $
//   Revision:            $Revision: 1.0 $
//
// 

// System includes

// External includes 

// Project includes
#include "includes/define.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/line_2d.h"
#include "kPoisson.h"
#include "includes/variables.h"


namespace Kratos
{
	//for Poisson Application
	KRATOS_CREATE_VARIABLE(double, DUMMY_UNKNOWN)
	KRATOS_CREATE_VARIABLE(double, DUMMY_MATERIAL)
	KRATOS_CREATE_VARIABLE(double, DUMMY_POINT_SOURCE)


	KratosR1PoissonApplication::KratosR1PoissonApplication():
		mPoisson2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>()))))
	{}

		
	void KratosR1PoissonApplication::Register()
	{

		//for Poisson Application
		KRATOS_REGISTER_VARIABLE(DUMMY_UNKNOWN)
		KRATOS_REGISTER_VARIABLE(DUMMY_MATERIAL)
		KRATOS_REGISTER_VARIABLE(DUMMY_POINT_SOURCE)


		// calling base class register to register Kratos components
		KratosApplication::Register();
		std::cout << "Initializing KratosR1PoissonApplication... " << std::endl;

		// Registering elements and conditions here
		KRATOS_REGISTER_ELEMENT("Poisson2D", mPoisson2D);

	}
	
}  // namespace Kratos.


