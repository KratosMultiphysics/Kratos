//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    @{KRATOS_APP_AUTHOR}
//


// System includes


// External includes


// Project includes
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/quadrilateral_2d_4.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/hexahedra_3d_27.h"
#include "helmholtz_application.h"
#include "helmholtz_application_variables.h"


namespace Kratos {

KratosHelmholtzApplication::KratosHelmholtzApplication():
    KratosApplication("HelmholtzApplication"),
      mHelmholtz2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
      mHelmholtz3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
      mHelmholtz3D8N(0, Element::GeometryType::Pointer(new Hexahedra3D8<Node<3> >(Element::GeometryType::PointsArrayType(8)))),
      mHelmholtz3D27N(0, Element::GeometryType::Pointer(new Hexahedra3D27<Node<3> >(Element::GeometryType::PointsArrayType(27)))),
      mHelmholtzVec2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
      mHelmholtzVec3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
      mHelmholtzVec3D8N(0, Element::GeometryType::Pointer(new Hexahedra3D8<Node<3> >(Element::GeometryType::PointsArrayType(8)))),
      mHelmholtzVec3D27N(0, Element::GeometryType::Pointer(new Hexahedra3D27<Node<3> >(Element::GeometryType::PointsArrayType(27))))
    {}

void KratosHelmholtzApplication::Register()
{
  KRATOS_INFO("") << "Initializing KratosHelmholtzApplication..." << std::endl;

  KRATOS_REGISTER_VARIABLE(HELMHOLTZ_DIRECTION)
  KRATOS_REGISTER_VARIABLE(HELMHOLTZ_POISSON_RATIO)

  KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( HELMHOLTZ_VARS )
  KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( HELMHOLTZ_SOURCE )

  KRATOS_REGISTER_ELEMENT("HelmholtzElement2D3N", mHelmholtz2D3N);
  KRATOS_REGISTER_ELEMENT("HelmholtzElement3D4N", mHelmholtz3D4N);
  KRATOS_REGISTER_ELEMENT("HelmholtzElement3D8N", mHelmholtz3D8N);
  KRATOS_REGISTER_ELEMENT("HelmholtzElement3D27N", mHelmholtz3D27N);

  KRATOS_REGISTER_ELEMENT("HelmholtzVecElement2D3N", mHelmholtzVec2D3N);
  KRATOS_REGISTER_ELEMENT("HelmholtzVecElement3D4N", mHelmholtzVec3D4N);
  KRATOS_REGISTER_ELEMENT("HelmholtzVecElement3D8N", mHelmholtzVec3D8N);
  KRATOS_REGISTER_ELEMENT("HelmholtzVecElement3D27N",mHelmholtzVec3D27N);  


}
}  // namespace Kratos.
