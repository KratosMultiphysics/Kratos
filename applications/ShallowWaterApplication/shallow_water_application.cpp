//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author:  Miguel Masó Sotomayor$
//   Date:                $Date:            april 26 2017$
//   Revision:            $Revision:                 1.3 $
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

	KratosShallowWaterApplication::KratosShallowWaterApplication():
	mProjectedSWE ( 0, Element::GeometryType::Pointer( new Triangle2D3<Node<3> >( Element::GeometryType::PointsArrayType (3) ) ) )
	{}
	
	void KratosShallowWaterApplication::Register()
	{
		// Calling base class register to register Kratos components
		KratosApplication::Register();
		std::cout << "Initializing KratosShallowWaterApplication... " << std::endl;

		KRATOS_REGISTER_VARIABLE(BATHYMETRY)
		
		KRATOS_REGISTER_VARIABLE(HEIGHT)
		KRATOS_REGISTER_VARIABLE(PROJECTED_HEIGHT)
		KRATOS_REGISTER_VARIABLE(DELTA_HEIGHT)
		KRATOS_REGISTER_VARIABLE(PROJECTED_VELOCITY)
		KRATOS_REGISTER_VARIABLE(DELTA_VELOCITY)
		
		KRATOS_REGISTER_VARIABLE(MEAN_SIZE)
		KRATOS_REGISTER_VARIABLE(MEAN_VEL_OVER_ELEM_SIZE)

		KRATOS_REGISTER_VARIABLE(SCALAR_VELOCITY_X)
		KRATOS_REGISTER_VARIABLE(SCALAR_VELOCITY_Y)
		KRATOS_REGISTER_VARIABLE(SCALAR_PROJECTED_VELOCITY_X)
		KRATOS_REGISTER_VARIABLE(SCALAR_PROJECTED_VELOCITY_Y)

		// Registering elements and conditions here
		KRATOS_REGISTER_ELEMENT("ProjectedSWE", mProjectedSWE)   // mesh stage element

 	}

}  // namespace Kratos.


