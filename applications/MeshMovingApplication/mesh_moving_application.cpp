//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license:
// kratos/license.txt
//
//  Main authors:    Andreas Winterstein (a.winterstein@tum.de)
//

// System includes

// External includes

// Project includes
//#include "includes/define.h"
#include "mesh_moving_application.h"

#include "geometries/hexahedra_3d_8.h"
#include "geometries/prism_3d_15.h"
#include "geometries/prism_3d_6.h"
#include "geometries/quadrilateral_2d_4.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_3d_3.h"

namespace Kratos {
//

KratosMeshMovingApplication::KratosMeshMovingApplication()
    : KratosApplication("MeshMovingApplication"),
      mLaplacianMeshMovingElement2D3N(
          0, Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(
                 Element::GeometryType::PointsArrayType(3)))),
      mLaplacianMeshMovingElement2D4N(
          0, Element::GeometryType::Pointer(new Quadrilateral2D4<Node<3>>(
                 Element::GeometryType::PointsArrayType(4)))),
      mLaplacianMeshMovingElement3D4N(
          0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(
                 Element::GeometryType::PointsArrayType(4)))),
      mLaplacianMeshMovingElement3D8N(
          0, Element::GeometryType::Pointer(new Hexahedra3D8<Node<3>>(
                 Element::GeometryType::PointsArrayType(8)))),
      mStructuralMeshMovingElement2D3N(
          0, Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(
                 Element::GeometryType::PointsArrayType(3)))),
      mStructuralMeshMovingElement2D4N(
          0, Element::GeometryType::Pointer(new Quadrilateral2D4<Node<3>>(
                 Element::GeometryType::PointsArrayType(4)))),
      mStructuralMeshMovingElement3D4N(
          0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(
                 Element::GeometryType::PointsArrayType(4)))),
      mStructuralMeshMovingElement3D8N(
          0, Element::GeometryType::Pointer(new Hexahedra3D8<Node<3>>(
                 Element::GeometryType::PointsArrayType(8)))),
      mStructuralMeshMovingElement3D6N(
          0, Element::GeometryType::Pointer(new Prism3D6<Node<3>>(
                 Element::GeometryType::PointsArrayType(6)))),
      mStructuralMeshMovingElement3D15N(
          0, Element::GeometryType::Pointer(new Prism3D15<Node<3>>(
                 Element::GeometryType::PointsArrayType(15)))),
      // Elements with empty geometry for automatic element generation by
      // GenerateMeshPart function
      mLaplacianMeshMovingElement(
          0, Element::GeometryType::Pointer(new Geometry<Node<3>>())),
      mStructuralMeshMovingElement(
          0, Element::GeometryType::Pointer(new Geometry<Node<3>>())) {}

void KratosMeshMovingApplication::Register() {
  // calling base class register to register Kratos components
  KratosApplication::Register();

  std::cout << "KRATOS __  __        _    __  __         _                       " << std::endl;
  std::cout << "      |  \\/  |___ __| |_ |  \\/  |_____ _(_)_ _  __ _           " << std::endl;
  std::cout << "      | |\\/| / -_|_-< ' \\| |\\/| / _ \\ V / | ' \\/ _` |       " << std::endl;
  std::cout << "      |_|  |_\\___/__/_||_|_|  |_\\___/\\_/|_|_||_\\__, |        " << std::endl;
  std::cout << "                                               |___/  application" << std::endl;
  std::cout << "Initializing KratosMeshMovingApplication...                      " << std::endl;


  // register elements
  KRATOS_REGISTER_ELEMENT("LaplacianMeshMovingElement2D3N",
                          mLaplacianMeshMovingElement2D3N);
  KRATOS_REGISTER_ELEMENT("LaplacianMeshMovingElemtent3D4N",
                          mLaplacianMeshMovingElement3D4N);
  KRATOS_REGISTER_ELEMENT("LaplacianMeshMovingElement2D4N",
                          mLaplacianMeshMovingElement2D4N);
  KRATOS_REGISTER_ELEMENT("LaplacianMeshMovingElement3D8N",
                          mLaplacianMeshMovingElement3D8N);
  KRATOS_REGISTER_ELEMENT("LaplacianMeshMovingElement",
                          mLaplacianMeshMovingElement);
  KRATOS_REGISTER_ELEMENT("StructuralMeshMovingElement2D3N",
                          mStructuralMeshMovingElement2D3N);
  KRATOS_REGISTER_ELEMENT("StructuralMeshMovingElement2D4N",
                          mStructuralMeshMovingElement2D4N);
  KRATOS_REGISTER_ELEMENT("StructuralMeshMovingElement3D4N",
                          mStructuralMeshMovingElement3D4N);
  KRATOS_REGISTER_ELEMENT("StructuralMeshMovingElement3D8N",
                          mStructuralMeshMovingElement3D8N);
  KRATOS_REGISTER_ELEMENT("StructuralMeshMovingElement3D6N",
                          mStructuralMeshMovingElement3D6N);
  KRATOS_REGISTER_ELEMENT("StructuralMeshMovingElement3D15N",
                          mStructuralMeshMovingElement3D15N);
  KRATOS_REGISTER_ELEMENT("StructuralMeshMovingElement",
                          mStructuralMeshMovingElement);

}

} // namespace Kratos.
