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
#include "geometries/triangle_3d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/line_2d.h"
#include "pfem_2_application.h"
#include "includes/variables.h"
#include "geometries/point_2d.h"   //add the point_2d (necessary for the pointsource)
#include "geometries/point_3d.h"   //add the point_2d (necessary for the pointsource)
//#include "utilities/enrich_2d_2dofs.h"

namespace Kratos
{
	//Example
 	KRATOS_CREATE_VARIABLE(double, PRESS_GRADIENT_JUMP)
	KRATOS_CREATE_VARIABLE(double, PRESS_DISCONTINUITY);
	KRATOS_CREATE_VARIABLE(double, INV_LAPLACIAN_ENRICH)
	KRATOS_CREATE_VARIABLE(double, ENRICH_RHS)
	KRATOS_CREATE_VARIABLE(double, G_VALUE)
	KRATOS_CREATE_VARIABLE(double, GRADIENT_DISCONTINUITY)
	KRATOS_CREATE_VARIABLE(double, PREVIOUS_ITERATION_PRESSURE)
	KRATOS_CREATE_VARIABLE(double, FIRST_ITERATION_PRESSURE)
	KRATOS_CREATE_VARIABLE(double, VELOCITY_OVER_ELEM_SIZE)
	KRATOS_CREATE_VARIABLE(double, MEAN_SIZE)
	KRATOS_CREATE_VARIABLE(Vector, ENRICH_LHS_ROW_3D)
	//KRATOS_CREATE_VARIABLE(double, IS_AIR)

	KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(ENRICH_LHS_ROW)
	KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(ENRICH_PRESS_PROJ_NEGATIVE)
	KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(ENRICH_PRESS_PROJ_POSITIVE)
	KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(SURFACE_NORMAL)
	KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(SURFACE_COORDINATES)
	KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(PRESS_PROJ_NO_RO)
	KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(DELTA_VELOCITY)
	
//	KRATOS_CREATE_VARIABLE(double, NODAL_AREA);
//

 	KratosPFEM2Application::KratosPFEM2Application():
 			mPFEM22D    ( 0, Element::GeometryType::Pointer( new Triangle2D3<Node<3> >(  Element::GeometryType::PointsArrayType (3, Node<3>() ) ) ) ),
 			mPFEM23D    ( 0, Element::GeometryType::Pointer( new Tetrahedra3D4<Node<3> >(  Element::GeometryType::PointsArrayType (4, Node<3>() ) ) ) ),
 			mFixedVelocity2D    ( 0, Condition::GeometryType::Pointer( new Point2D<Node<3> >(  Element::GeometryType::PointsArrayType (1, Node<3>() ) ) ) ),
 			mFixedVelocity3D    ( 0, Condition::GeometryType::Pointer( new Point3D<Node<3> >(  Element::GeometryType::PointsArrayType (1, Node<3>() ) ) ) ),
 			mPFEM2DiscontinuousVelocity2D    ( 0, Element::GeometryType::Pointer( new Triangle2D3<Node<3> >(  Element::GeometryType::PointsArrayType (3, Node<3>() ) ) ) )
 	{}
 	
 	void KratosPFEM2Application::Register()
 	{
 		// calling base class register to register Kratos components
 		KratosApplication::Register();
 		std::cout << "Initializing KratosPFEM2Application... " << std::endl;
 		KRATOS_REGISTER_ELEMENT("PFEM22D", mPFEM22D);
 		KRATOS_REGISTER_ELEMENT("PFEM23D", mPFEM23D);
 		
 		KRATOS_REGISTER_CONDITION("FixedVelocity2D", mFixedVelocity2D);
		KRATOS_REGISTER_CONDITION("FixedVelocity3D", mFixedVelocity3D);
 		
 		KRATOS_REGISTER_ELEMENT("PFEM2DiscontinuousVelocity2D",mPFEM2DiscontinuousVelocity2D);
 
 		KRATOS_REGISTER_VARIABLE( PRESS_GRADIENT_JUMP )
 		KRATOS_REGISTER_VARIABLE(PRESS_DISCONTINUITY);
 		KRATOS_REGISTER_VARIABLE(INV_LAPLACIAN_ENRICH)
 		KRATOS_REGISTER_VARIABLE(ENRICH_RHS)
 		KRATOS_REGISTER_VARIABLE(GRADIENT_DISCONTINUITY)
 		KRATOS_REGISTER_VARIABLE(PREVIOUS_ITERATION_PRESSURE)
 		KRATOS_REGISTER_VARIABLE(FIRST_ITERATION_PRESSURE)
 		KRATOS_REGISTER_VARIABLE(MEAN_SIZE)
 		KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(ENRICH_PRESS_PROJ_NEGATIVE)
 		KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(ENRICH_PRESS_PROJ_POSITIVE)
 		KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(SURFACE_NORMAL)
 		KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(SURFACE_COORDINATES)
 		KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(PRESS_PROJ_NO_RO)
 		KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(DELTA_VELOCITY)
 		KRATOS_REGISTER_VARIABLE(G_VALUE)
 		KRATOS_REGISTER_VARIABLE(VELOCITY_OVER_ELEM_SIZE)
 		KRATOS_REGISTER_VARIABLE(ENRICH_LHS_ROW_3D)
 		//KRATOS_REGISTER_VARIABLE(IS_AIR)
// 		KRATOS_REGISTER_VARIABLE(NODAL_AREA);

 
 	}

}  // namespace Kratos.


