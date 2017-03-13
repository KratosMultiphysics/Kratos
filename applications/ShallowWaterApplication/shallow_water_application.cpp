//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author:  $
//   Date:                $Date:  $
//   Revision:            $Revision: 1.3 $
//
// 



// System includes


// External includes 


// Project includes
#include "includes/define.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/line_2d.h"
#include "geometries/point_2d.h"
#include "shallow_water_application.h"
#include "includes/variables.h"


namespace Kratos
{
	KRATOS_CREATE_VARIABLE(double, POINT_HEAT_SOURCE)
	// CONDUCTIVITY and TEMPERATURE are already included in kernel

 	KratosShallowWaterApplication::KratosShallowWaterApplication():  //constructor
		mPoisson2D   ( 0, Element::GeometryType::Pointer( new Triangle2D3<Node<3> >( Element::GeometryType::PointsArrayType (3) ) ) ),
		mPointSource ( 0, Element::GeometryType::Pointer( new Point2D <Node<3>    >( Element::GeometryType::PointsArrayType (1) ) ) )
 	{}
 	
 	void KratosShallowWaterApplication::Register()
 	{
 		// calling base class register to register Kratos components
 		KratosApplication::Register();
 		std::cout << "Initializing KratosShallowWaterApplication... " << std::endl;
 
		// Registering variables
 		KRATOS_REGISTER_VARIABLE( POINT_HEAT_SOURCE );

		// Registering elements and conditions
		KRATOS_REGISTER_ELEMENT( "Poisson2D", mPoisson2D );
		KRATOS_REGISTER_CONDITION( "PointSource", mPointSource )
 
 	}

}  // namespace Kratos.


