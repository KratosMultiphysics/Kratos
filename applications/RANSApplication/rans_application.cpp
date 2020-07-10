//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya (https://github.com/sunethwarna)
//

// System includes

// External includes

// Project includes
#include "rans_application.h"
#include "geometries/line_2d_2.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/triangle_2d_3.h"
#include "rans_application_variables.h"

namespace Kratos
{
KratosRANSApplication::KratosRANSApplication()
    : KratosApplication("RANSApplication"),
      // incompressible potential flow elements
      mIncompressiblePotentialFlowVelocity2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mIncompressiblePotentialFlowVelocity3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(Element::GeometryType::PointsArrayType(4)))),
      mIncompressiblePotentialFlowPressure2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mIncompressiblePotentialFlowPressure3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(Element::GeometryType::PointsArrayType(4))))
{
}

void KratosRANSApplication::Register()
{
    KRATOS_INFO("") << "Initializing KratosRANSApplication..." << std::endl;

    // incompressible potential flow specific variables
    KRATOS_REGISTER_VARIABLE( VELOCITY_POTENTIAL )
    KRATOS_REGISTER_VARIABLE( PRESSURE_POTENTIAL )
    KRATOS_REGISTER_VARIABLE( RANS_IS_INLET )
    KRATOS_REGISTER_VARIABLE( RANS_IS_OUTLET )
    KRATOS_REGISTER_VARIABLE( RANS_IS_STRUCTURE )

    // registering elements
    // registering incompressible potential flow elements
    KRATOS_REGISTER_ELEMENT("RansIncompressiblePotentialFlowVelocity2D3N", mIncompressiblePotentialFlowVelocity2D);
    KRATOS_REGISTER_ELEMENT("RansIncompressiblePotentialFlowVelocity3D4N", mIncompressiblePotentialFlowVelocity3D);
    KRATOS_REGISTER_ELEMENT("RansIncompressiblePotentialFlowPressure2D3N", mIncompressiblePotentialFlowPressure2D);
    KRATOS_REGISTER_ELEMENT("RansIncompressiblePotentialFlowPressure3D4N", mIncompressiblePotentialFlowPressure3D);
}
} // namespace Kratos.
