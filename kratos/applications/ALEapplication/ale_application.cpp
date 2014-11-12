//
//   Project Name:        Kratos
//   Last Modified by:    $Author: AMini $
//   Date:                $Date: Oct 2014 $
//   Revision:            $Revision: 1.4 $
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
#include "ale_application.h"
#include "includes/variables.h"


namespace Kratos
{
//Create variables


//

KratosALEApplication::KratosALEApplication():
    mLaplacianMeshMovingElem2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>())))),
    mLaplacianMeshMovingElem3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4, Node<3>())))),
    mLaplacianComponentwiseMeshMovingElem2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>())))),
    mLaplacianComponentwiseMeshMovingElem3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4, Node<3>())))),
    mStructuralMeshMovingElement2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>())))),
    mStructuralMeshMovingElement3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4, Node<3>())))),
    mLaplacianComponentwiseMeshMovingElem2DStrainbased(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>())))),
    //mLaplacianComponentwiseMeshMovingElem3DStrainbased(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4, Node<3>())))),
    mStructuralMeshMovingElem2DNonlin(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>())))),
    mStructuralMeshMovingElem3DNonlin(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4, Node<3>()))))
{}


void KratosALEApplication::Register()
{
    // calling base class register to register Kratos components
    KratosApplication::Register();
    std::cout << "Initializing KratosALEApplication... " << std::endl;


    KRATOS_REGISTER_ELEMENT("LaplacianMeshMovingElem2D", mLaplacianMeshMovingElem2D);
    KRATOS_REGISTER_ELEMENT("LaplacianMeshMovingElemt3D", mLaplacianMeshMovingElem3D);
    KRATOS_REGISTER_ELEMENT("LaplacianComponentwiseMeshMovingElemt2D", mLaplacianComponentwiseMeshMovingElem2D);
    KRATOS_REGISTER_ELEMENT("LaplacianComponentwiseMeshMovingElemt3D", mLaplacianComponentwiseMeshMovingElem3D);
    KRATOS_REGISTER_ELEMENT("StructuralMeshMovingElement2D", mStructuralMeshMovingElement2D);
    KRATOS_REGISTER_ELEMENT("StructuralMeshMovingElement3D", mStructuralMeshMovingElement3D);
    KRATOS_REGISTER_ELEMENT("LaplacianComponentwiseMeshMovingElemt2D", mLaplacianComponentwiseMeshMovingElem2DStrainbased);
    //KRATOS_REGISTER_ELEMENT("LaplacianComponentwiseMeshMovingElemt3D", mLaplacianComponentwiseMeshMovingElem3DStrainbased);
    KRATOS_REGISTER_ELEMENT("StructuralMeshMovingElemt2DNonlin", mStructuralMeshMovingElem2DNonlin);
    KRATOS_REGISTER_ELEMENT("StructuralMeshMovingElemt3DNonlin", mStructuralMeshMovingElem3DNonlin);

}

}  // namespace Kratos.


