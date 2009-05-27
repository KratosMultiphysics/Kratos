//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: rrossi $
//   Date:                $Date: 2008-05-06 15:11:09 $
//   Revision:            $Revision: 1.6 $
//
// 



// System includes


// External includes 


// Project includes
#include "includes/define.h"

#include "PFEM_application.h"
#include "includes/variables.h"
#include "geometries/triangle_2d_3.h"


namespace Kratos
{

/*	KRATOS_CREATE_VARIABLE(double, NODAL_AREA)
	KRATOS_CREATE_VARIABLE(double, NODAL_H)
	KRATOS_CREATE_VARIABLE(double, IS_STRUCTURE)
	KRATOS_CREATE_VARIABLE(double, IS_FLUID)
	KRATOS_CREATE_VARIABLE(double, IS_BOUNDARY)
	KRATOS_CREATE_VARIABLE(double, IS_FREE_SURFACE)
	KRATOS_CREATE_VARIABLE(double, IS_FREE_SURFACE)
	KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(NORMAL_TO_WALL)
*/
//	KRATOS_CREATE_VARIABLE(double, PRESSURE_OLD_IT)
//	KRATOS_CREATE_VARIABLE(double, NODAL_MASS)
//	KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(FRACT_VEL)
	KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(VAUX)

	KratosPFEMApplication::KratosPFEMApplication():	
		mFreeSurfaceCond2d(0, Element::GeometryType::Pointer(new Geometry<Node<3> >(Element::GeometryType::PointsArrayType(2, Node<3>()))))
	{}


	void KratosPFEMApplication::Register()
	{
		// calling base class register to register Kratos components
		KratosApplication::Register();
		std::cout << "Initializing Kratos PFEMApplication... " << std::endl;

/*		KRATOS_REGISTER_VARIABLE( NODAL_AREA)
		KRATOS_REGISTER_VARIABLE( NODAL_H)
		KRATOS_REGISTER_VARIABLE( IS_STRUCTURE)
		KRATOS_REGISTER_VARIABLE( IS_FLUID)
		KRATOS_REGISTER_VARIABLE( IS_BOUNDARY)
		KRATOS_REGISTER_VARIABLE( IS_FREE_SURFACE)
		KRATOS_REGISTER_VARIABLE( IS_FREE_SURFACE)
		KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(NORMAL_TO_WALL)
*/
		KRATOS_REGISTER_VARIABLE( PRESSURE_OLD_IT)
		KRATOS_REGISTER_VARIABLE( NODAL_MASS)
		KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(FRACT_VEL)
		KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(VAUX)

		KRATOS_REGISTER_CONDITION("FreeSurfaceCond2d", mFreeSurfaceCond2d);
	}


	
}  // namespace Kratos.


