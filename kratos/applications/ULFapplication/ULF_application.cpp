//   
//   Project Name:        Kratos       

//   Last Modified by:    $Author: virginia $
//   Date:                $Date: 2008-11-26 15:05:54 $
//   Revision:            $Revision: 1.11 $
//
// 


  
// System includes


// External includes 

 
// Project includes
#include "includes/define.h"

#include "ULF_application.h"
#include "includes/variables.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/tetrahedra_3d_4.h"


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
	//KRATOS_CREATE_VARIABLE(double, IS_LAGRANGIAN_INLET)

	
	KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(PRESSURE_FORCE)
	KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(DISP_FRAC)
	KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(VAUX)

KratosULFApplication::KratosULFApplication():
	mUpdatedLagrangianFluid2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>())))),
	mUpdatedLagrangianFluid3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4 <Node<3> >(Element::GeometryType::PointsArrayType(4, Node<3>())))),
	mUpdatedLagrangianFluid2Dinc(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>())))),
	mUpdatedLagrangianFluid3Dinc(0, Element::GeometryType::Pointer(new Tetrahedra3D4 <Node<3> >(Element::GeometryType::PointsArrayType(4, Node<3>())))),
	//new one - mis of frac step and ulf_inc	
	mUlfFrac2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>()))))
	{}
	

	void KratosULFApplication::Register()
	{
		// calling base class register to register Kratos components
		KratosApplication::Register();
		std::cout << "Initializing Kratos ULFApplication... " << std::endl;

/*		KRATOS_REGISTER_VARIABLE( NODAL_AREA)
		KRATOS_REGISTER_VARIABLE( NODAL_H)
		KRATOS_REGISTER_VARIABLE( IS_STRUCTURE)
		KRATOS_REGISTER_VARIABLE( IS_FLUID)
		KRATOS_REGISTER_VARIABLE( IS_BOUNDARY)
		KRATOS_REGISTER_VARIABLE( IS_FREE_SURFACE)
		KRATOS_REGISTER_VARIABLE( IS_FREE_SURFACE)
		KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(NORMAL_TO_WALL)
*/
	//	KRATOS_REGISTER_VARIABLE(IS_LAGRANGIAN_INLET)

		KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(PRESSURE_FORCE)
		KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(DISP_FRAC)
		KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(VAUX)
		
		KRATOS_REGISTER_ELEMENT("UpdatedLagrangianFluid2D", mUpdatedLagrangianFluid2D);
		KRATOS_REGISTER_ELEMENT("UpdatedLagrangianFluid3D", mUpdatedLagrangianFluid3D);
		KRATOS_REGISTER_ELEMENT("UpdatedLagrangianFluid2Dinc", mUpdatedLagrangianFluid2Dinc);
		KRATOS_REGISTER_ELEMENT("UpdatedLagrangianFluid3Dinc", mUpdatedLagrangianFluid3Dinc);
		//
		KRATOS_REGISTER_ELEMENT("UlfFrac2D", mUlfFrac2D);
		
	}


	
}  // namespace Kratos.


