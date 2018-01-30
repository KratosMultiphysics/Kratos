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
#include "pod_application.h"
#include "includes/variables.h"


namespace Kratos
{
//Example
KRATOS_CREATE_VARIABLE(Vector, POD_VELOCITY_X)
KRATOS_CREATE_VARIABLE(Vector, POD_VELOCITY_Y)
KRATOS_CREATE_VARIABLE(Vector, POD_VELOCITY_Z)
KRATOS_CREATE_VARIABLE(Vector, POD_PRESSURE)
// 	KRATOS_CREATE_VARIABLE(double, IS_INTERFACE);
// 	KRATOS_CREATE_VARIABLE(double, NODAL_AREA);


KratosPodApplication::KratosPodApplication()
/*		:
		mElem2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>())))),
		mElem3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4, Node<3>()))))*/
{}

void KratosPodApplication::Register()
{
    // calling base class register to register Kratos components
    KratosApplication::Register();
    std::cout << "Initializing KratosPodApplication... " << std::endl;

    KRATOS_REGISTER_VARIABLE( POD_VELOCITY_X)
    KRATOS_REGISTER_VARIABLE( POD_VELOCITY_Y)
    KRATOS_REGISTER_VARIABLE( POD_VELOCITY_Z)
    KRATOS_REGISTER_VARIABLE( POD_PRESSURE)

// 		KRATOS_REGISTER_ELEMENT("Elem2D", mElem2D);
// 		KRATOS_REGISTER_ELEMENT("Elemt3D", mElem3D);

}

}  // namespace Kratos.


