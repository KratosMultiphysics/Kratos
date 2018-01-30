//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Author Julio Marti
//


// Project includes
#include "includes/define.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/line_2d.h"
#include "radiation_application.h"
#include "includes/variables.h"


namespace Kratos
{

	//KRATOS_CREATE_VARIABLE(double,  TEMP_CONV_PROJ)	
	//KRATOS_CREATE_VARIABLE(double,  EMISSIVITY)	
			
	//KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(CONVECTION_VELOCITY)

	KratosRadiationApplication::KratosRadiationApplication():
	mRad2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
        mRad3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
	mRadFace2D(0, Element::GeometryType::Pointer(new Geometry<Node<3> >(Element::GeometryType::PointsArrayType(2)))),
	mRadFace3D(0, Element::GeometryType::Pointer(new Triangle3D3<Node<3> >(Element::GeometryType::PointsArrayType(3))))
	{}


		
	void KratosRadiationApplication::Register()
	{
		// calling base class register to register Kratos components
		KratosApplication::Register();
		std::cout << "Initializing KratosRadiationApplication... " << std::endl;

		
		// Registering elements and conditions here
		KRATOS_REGISTER_ELEMENT("Rad2D", mRad2D);
		KRATOS_REGISTER_ELEMENT("Rad3D", mRad3D);

		KRATOS_REGISTER_CONDITION("RadFace2D", mRadFace2D);
                KRATOS_REGISTER_CONDITION("RadFace3D", mRadFace3D);


	}




}  // namespace Kratos.


