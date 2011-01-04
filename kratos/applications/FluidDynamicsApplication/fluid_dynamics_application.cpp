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
                mVMS2DSmagorinsky(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>())))),
                mVMS3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4, Node<3>())))),
                mVMS3DSmagorinsky(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4, Node<3>()))))//,
//                mFluidPeriodicCondition2D( 0, Element::GeometryType::Pointer( new Line2D2<Node<3> >( Element::GeometryType::PointsArrayType( 2, Node<3>() ) ) ) )
 	{}
 	
 	void KratosFluidDynamicsApplication::Register()
 	{
 		// calling base class register to register Kratos components
 		KratosApplication::Register();
 		std::cout << "Initializing KratosFluidDynamicsApplication... " << std::endl;
 
		// Register Variables (defined in fluid_dynamics_application_variables.h)
 		KRATOS_REGISTER_VARIABLE(VORTICITY)

		// Register Elements
                KRATOS_REGISTER_ELEMENT("VMS2D",mVMS2D)
                KRATOS_REGISTER_ELEMENT("VMS2DSmagorinsky",mVMS2DSmagorinsky)
                KRATOS_REGISTER_ELEMENT("VMS3D",mVMS3D)
                KRATOS_REGISTER_ELEMENT("VMS3DSmagorinsky",mVMS3DSmagorinsky)

		// Register Conditions
 //               KRATOS_REGISTER_CONDITION("FluidPeriodicCondition2D",mFluidPeriodicCondition2D)
 	}


}  // namespace Kratos.


