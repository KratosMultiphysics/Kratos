//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: jcotela $
//   Date:                $Date: 2010-11-11 $
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
#include "fluid_dynamics_application.h"
#include "includes/variables.h"


namespace Kratos
{

 	KratosFluidDynamicsApplication::KratosFluidDynamicsApplication():
                mVMS2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>())))),
                mVMS3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4, Node<3>())))),
                mBinghamVMS2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>())))),
                mBinghamVMS3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4, Node<3>())))),
                mDynamicVMS2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>())))),
                mDynamicVMS3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4, Node<3>())))),
                mSpalartAllmaras2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>())))),
                mSpalartAllmaras3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4, Node<3>())))),
                mMonolithicWallCondition2D(0, Element::GeometryType::Pointer( new Line2D2<Node<3> >( Element::GeometryType::PointsArrayType( 2, Node<3>() ) ) ) ),
                mMonolithicWallCondition3D(0, Element::GeometryType::Pointer( new Triangle3D3<Node<3> >( Element::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) ),
                mPeriodicCondition(0, Element::GeometryType::Pointer( new Line2D2<Node<3> >( Element::GeometryType::PointsArrayType( 2, Node<3>() ) ) ) )//,
 	{}
 	
 	void KratosFluidDynamicsApplication::Register()
 	{
 		// calling base class register to register Kratos components
 		KratosApplication::Register();
 		std::cout << "Initializing KratosFluidDynamicsApplication... " << std::endl;
 
		// Register Variables (defined in fluid_dynamics_application_variables.h)
                KRATOS_REGISTER_VARIABLE(PATCH_INDEX);
                KRATOS_REGISTER_VARIABLE(TRACK_SUBSCALES);
                KRATOS_REGISTER_VARIABLE(TAUONE);
                KRATOS_REGISTER_VARIABLE(TAUTWO);
                KRATOS_REGISTER_VARIABLE(Y_WALL);
                KRATOS_REGISTER_VARIABLE(SUBSCALE_PRESSURE);
//                KRATOS_REGISTER_VARIABLE(C_SMAGORINSKY)
                KRATOS_REGISTER_VARIABLE(PERIODIC_VARIABLES);
                KRATOS_REGISTER_VARIABLE(SUBSCALE_VELOCITY);
                KRATOS_REGISTER_VARIABLE(VORTICITY);
                KRATOS_REGISTER_VARIABLE(COARSE_VELOCITY);

		// Register Elements
                KRATOS_REGISTER_ELEMENT("VMS2D",mVMS2D)
                KRATOS_REGISTER_ELEMENT("VMS3D",mVMS3D)
                KRATOS_REGISTER_ELEMENT("BinghamVMS2D",mBinghamVMS2D)
                KRATOS_REGISTER_ELEMENT("BinghamVMS3D",mBinghamVMS3D)
                KRATOS_REGISTER_ELEMENT("DynamicVMS2D",mDynamicVMS2D)
                KRATOS_REGISTER_ELEMENT("DynamicVMS3D",mDynamicVMS3D)

                KRATOS_REGISTER_ELEMENT("SpalartAllmaras2D",mSpalartAllmaras2D)
                KRATOS_REGISTER_ELEMENT("SpalartAllmaras3D",mSpalartAllmaras3D)

                // Register Conditions
                KRATOS_REGISTER_CONDITION("MonolithicWallCondition2D",mMonolithicWallCondition2D)
                KRATOS_REGISTER_CONDITION("MonolithicWallCondition3D",mMonolithicWallCondition3D)
                KRATOS_REGISTER_CONDITION("PeriodicCondition",mPeriodicCondition)
 	}


}  // namespace Kratos.


