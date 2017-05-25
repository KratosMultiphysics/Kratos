//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
// 
//  Project Name:        Kratos       
//  Last Modified by:    $Author:  Miguel Masó Sotomayor
//  Date:                $Date:  april 26 2017
//  Revision:            $Revision: 1.2 $
//
//


// System includes


// External includes 


// Project includes
#include "includes/define.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/line_2d.h"
#include "geometries/point_2d.h"
#include "pure_diffusion_application.h"
#include "includes/variables.h"


namespace Kratos
{

	KratosPureDiffusionApplication::KratosPureDiffusionApplication():
		mPoisson2D    ( 0, Element::GeometryType::Pointer( new Triangle2D3<Node<3> >( Element::GeometryType::PointsArrayType (3) ) ) ),
		mProjectedSWE ( 0, Element::GeometryType::Pointer( new Triangle2D3<Node<3> >( Element::GeometryType::PointsArrayType (3) ) ) ),
		mPointSource  ( 0, Element::GeometryType::Pointer (new Point2D <Node<3>    >( Element::GeometryType::PointsArrayType (1) ) ) )
	{}
	
	void KratosPureDiffusionApplication::Register()
	{
		// calling base class register to register Kratos components
		KratosApplication::Register();
		std::cout << "Initializing KratosPureDiffusionApplication... " << std::endl;

		KRATOS_REGISTER_VARIABLE(POINT_HEAT_SOURCE)
		KRATOS_REGISTER_VARIABLE(BATHYMETRY)
		KRATOS_REGISTER_VARIABLE(HEIGHT)
		KRATOS_REGISTER_VARIABLE(PROJECTED_HEIGHT)
		KRATOS_REGISTER_VARIABLE(PROJECTED_VELOCITY)

		// Registering elements and conditions here
		KRATOS_REGISTER_ELEMENT("Poisson2D", mPoisson2D)         // here is our element
		KRATOS_REGISTER_ELEMENT("ProjectedSWE", mProjectedSWE)   // here is another element
		KRATOS_REGISTER_CONDITION( "PointSource", mPointSource ) // and our condition

	}

}  // namespace Kratos.


