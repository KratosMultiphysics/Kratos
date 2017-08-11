//   
//   Project Name:        Kratos       
//   Last Modified by:    Miguel Masó Sotomayor
//   Date:                April 26th 2017
//   Revision:            1.4
//
// 



// System includes


// External includes 


// Project includes
#include "geometries/triangle_2d_3.h"
#include "geometries/quadrilateral_2d_4.h"
#include "geometries/line_2d.h"
#include "geometries/point_2d.h"
#include "includes/define.h"
#include "includes/variables.h"
#include "shallow_water_application.h"


namespace Kratos
{

	KratosShallowWaterApplication::KratosShallowWaterApplication():
	
	mPrimitiveVarElement2D3N( 0, Element::GeometryType::Pointer( new Triangle2D3<Node<3>      >( Element::GeometryType::PointsArrayType (3) ) ) ),
	mPrimitiveVarElement2D4N( 0, Element::GeometryType::Pointer( new Quadrilateral2D4<Node<3> >( Element::GeometryType::PointsArrayType (4) ) ) ),
	
	mConservedVarElement2D3N( 0, Element::GeometryType::Pointer( new Triangle2D3<Node<3>      >( Element::GeometryType::PointsArrayType (3) ) ) ),
	mConservedVarElement2D4N( 0, Element::GeometryType::Pointer( new Quadrilateral2D4<Node<3> >( Element::GeometryType::PointsArrayType (4) ) ) ),
	
	mRainCondition2D3N( 0, Element::GeometryType::Pointer( new Triangle2D3<Node<3>      >( Element::GeometryType::PointsArrayType (3) ) ) ),
	mRainCondition2D4N( 0, Element::GeometryType::Pointer( new Quadrilateral2D4<Node<3> >( Element::GeometryType::PointsArrayType (4) ) ) )
	
	{}
	
	void KratosShallowWaterApplication::Register()
	{
		// Calling base class register to register Kratos components
		KratosApplication::Register();
		
		std::cout << " KRATOS      |          |   |                        " << std::endl;
		std::cout << "        __|   _ \\  _` | |   |    _ \\        /      " << std::endl;
		std::cout << "      \\__ `  |  | (   | |   |   (   |      /        " << std::endl;
		std::cout << "      ____/ _| _|\\__,_|\\__|\\__|\\___/  _/ _/ WATER" << std::endl;
		std::cout << "Initializing KratosShallowWaterApplication...        " << std::endl;

		KRATOS_REGISTER_VARIABLE(BATHYMETRY)                            // Geometric definition of the problem
		KRATOS_REGISTER_VARIABLE(RAIN)                                  // Source term
		
		KRATOS_REGISTER_VARIABLE(HEIGHT)                                // Main variable
		KRATOS_REGISTER_VARIABLE(PROJECTED_HEIGHT)                      // Convected variable
		KRATOS_REGISTER_VARIABLE(DELTA_HEIGHT)                          // Variable to uodate particles
		KRATOS_REGISTER_VARIABLE(PROJECTED_VELOCITY)
		KRATOS_REGISTER_VARIABLE(DELTA_VELOCITY)
		KRATOS_REGISTER_VARIABLE(PROJECTED_MOMENTUM)
		KRATOS_REGISTER_VARIABLE(DELTA_MOMENTUM)
		
		KRATOS_REGISTER_VARIABLE(MEAN_SIZE)                             // Specific variable for PFEM2
		KRATOS_REGISTER_VARIABLE(MEAN_VEL_OVER_ELEM_SIZE)               // Specific variable for PFEM2

		//~ KRATOS_REGISTER_VARIABLE(SCALAR_VELOCITY_X)
		//~ KRATOS_REGISTER_VARIABLE(SCALAR_VELOCITY_Y)
		//~ KRATOS_REGISTER_VARIABLE(SCALAR_PROJECTED_VELOCITY_X)
		//~ KRATOS_REGISTER_VARIABLE(SCALAR_PROJECTED_VELOCITY_Y)

		// Registering elements and conditions here
		KRATOS_REGISTER_ELEMENT("PrimitiveVarElement2D3N", mPrimitiveVarElement2D3N)   // mesh stage element
		KRATOS_REGISTER_ELEMENT("PrimitiveVarElement2D4N", mPrimitiveVarElement2D4N)   // mesh stage element
		
		KRATOS_REGISTER_ELEMENT("ConservedVarElement2D3N", mConservedVarElement2D3N)   // mesh stage element
		KRATOS_REGISTER_ELEMENT("ConservedVarElement2D4N", mConservedVarElement2D4N)   // mesh stage element
		
		KRATOS_REGISTER_CONDITION("RainCondition2D3N", mRainCondition2D3N)
		KRATOS_REGISTER_CONDITION("RainCondition2D4N", mRainCondition2D4N)

	}

}  // namespace Kratos.


