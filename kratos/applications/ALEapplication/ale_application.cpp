/* *********************************************************
*
*   Last Modified by:    $Author: AMini $
*   Date:                $Date: Mai 2015 $
*   Revision:            $Revision: 1.3 $
*
* ***********************************************************/



// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/line_2d.h"
#include "geometries/hexahedra_3d_8.h"
#include "ale_application.h"
#include "includes/variables.h"


namespace Kratos
{
//Create variables


//

KratosALEApplication::KratosALEApplication():
    mLaplacianMeshMovingElement2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
    mLaplacianMeshMovingElement3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    mStructuralMeshMovingElement2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
    mStructuralMeshMovingElement3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    mStructuralMeshMovingElement2DNonlinear(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
    mStructuralMeshMovingElement3DNonlinear(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4))))
{}


void KratosALEApplication::Register()
{
    // calling base class register to register Kratos components
    KratosApplication::Register();
    std::cout << "Initializing KratosALEApplication... " << std::endl;


    KRATOS_REGISTER_ELEMENT("LaplacianMeshMovingElement2D", mLaplacianMeshMovingElement2D);
    KRATOS_REGISTER_ELEMENT("LaplacianMeshMovingElemtent3D", mLaplacianMeshMovingElement3D);
    KRATOS_REGISTER_ELEMENT("StructuralMeshMovingElement2D3N", mStructuralMeshMovingElement2D3N);
    KRATOS_REGISTER_ELEMENT("StructuralMeshMovingElement3D4N", mStructuralMeshMovingElement3D4N);
    KRATOS_REGISTER_ELEMENT("StructuralMeshMovingElement2DNonlinear", mStructuralMeshMovingElement2DNonlinear);
    KRATOS_REGISTER_ELEMENT("StructuralMeshMovingElement3DNonlinear", mStructuralMeshMovingElement3DNonlinear);

}

}  // namespace Kratos.


