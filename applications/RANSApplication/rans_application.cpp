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
#include "geometries/line_2d_2.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/triangle_2d_3.h"
#include "rans_application_variables.h"

// Include base h
#include "rans_application.h"

namespace Kratos
{
KratosRANSApplication::KratosRANSApplication()
    : KratosApplication("RANSApplication"),
      // k-epsilon turbulence model elements
      mRansEvmKEpsilonHighReKAFC2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mRansEvmKEpsilonHighReKAFC3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(Element::GeometryType::PointsArrayType(4)))),
      mRansEvmKEpsilonHighReEpsilonAFC2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mRansEvmKEpsilonHighReEpsilonAFC3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(Element::GeometryType::PointsArrayType(4)))),
      mRansEvmKEpsilonHighReKRFC2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mRansEvmKEpsilonHighReKRFC3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(Element::GeometryType::PointsArrayType(4)))),
      mRansEvmKEpsilonHighReEpsilonRFC2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mRansEvmKEpsilonHighReEpsilonRFC3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(Element::GeometryType::PointsArrayType(4)))),
      mRansEvmKEpsilonHighReKCWD2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mRansEvmKEpsilonHighReKCWD3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(Element::GeometryType::PointsArrayType(4)))),
      mRansEvmKEpsilonHighReEpsilonCWD2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mRansEvmKEpsilonHighReEpsilonCWD3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(Element::GeometryType::PointsArrayType(4))))
{
}

void KratosRANSApplication::Register()
{
    KRATOS_INFO("") << "Initializing KratosRANSApplication..." << std::endl;

    // residual based flux corrected stabilization variables
    KRATOS_REGISTER_VARIABLE( RANS_STABILIZATION_DISCRETE_UPWIND_OPERATOR_COEFFICIENT )
    KRATOS_REGISTER_VARIABLE( RANS_STABILIZATION_DIAGONAL_POSITIVITY_PRESERVING_COEFFICIENT )

    // k-epsilon-high-re turbulence modelling variables
    KRATOS_REGISTER_VARIABLE( TURBULENT_KINETIC_ENERGY )
    KRATOS_REGISTER_VARIABLE( TURBULENT_ENERGY_DISSIPATION_RATE )
    KRATOS_REGISTER_VARIABLE( TURBULENT_KINETIC_ENERGY_RATE )
    KRATOS_REGISTER_VARIABLE( TURBULENT_ENERGY_DISSIPATION_RATE_2 )
    KRATOS_REGISTER_VARIABLE( RANS_AUXILIARY_VARIABLE_1 )
    KRATOS_REGISTER_VARIABLE( RANS_AUXILIARY_VARIABLE_2 )
    KRATOS_REGISTER_VARIABLE( TURBULENT_KINETIC_ENERGY_SIGMA )
    KRATOS_REGISTER_VARIABLE( TURBULENT_ENERGY_DISSIPATION_RATE_SIGMA )
    KRATOS_REGISTER_VARIABLE( TURBULENCE_RANS_C_MU )
    KRATOS_REGISTER_VARIABLE( TURBULENCE_RANS_C1 )
    KRATOS_REGISTER_VARIABLE( TURBULENCE_RANS_C2 )

    // registering elements
    // registering k-epsilon algebraic flux correction based elements
    KRATOS_REGISTER_ELEMENT("RansEvmKEpsilonHighReKAFC2D3N", mRansEvmKEpsilonHighReKAFC2D);
    KRATOS_REGISTER_ELEMENT("RansEvmKEpsilonHighReKAFC3D4N", mRansEvmKEpsilonHighReKAFC3D);
    KRATOS_REGISTER_ELEMENT("RansEvmKEpsilonHighReEpsilonAFC2D3N", mRansEvmKEpsilonHighReEpsilonAFC2D);
    KRATOS_REGISTER_ELEMENT("RansEvmKEpsilonHighReEpsilonAFC3D4N", mRansEvmKEpsilonHighReEpsilonAFC3D);

    // registering k-epsilon residual fc based elements
    KRATOS_REGISTER_ELEMENT("RansEvmKEpsilonHighReKRFC2D3N", mRansEvmKEpsilonHighReKRFC2D);
    KRATOS_REGISTER_ELEMENT("RansEvmKEpsilonHighReKRFC3D4N", mRansEvmKEpsilonHighReKRFC3D);
    KRATOS_REGISTER_ELEMENT("RansEvmKEpsilonHighReEpsilonRFC2D3N", mRansEvmKEpsilonHighReEpsilonRFC2D);
    KRATOS_REGISTER_ELEMENT("RansEvmKEpsilonHighReEpsilonRFC3D4N", mRansEvmKEpsilonHighReEpsilonRFC3D);

    // registering k-epsilon cross wind stabilized elements
    KRATOS_REGISTER_ELEMENT("RansEvmKEpsilonHighReKCWD2D3N", mRansEvmKEpsilonHighReKCWD2D);
    KRATOS_REGISTER_ELEMENT("RansEvmKEpsilonHighReKCWD3D4N", mRansEvmKEpsilonHighReKCWD3D);
    KRATOS_REGISTER_ELEMENT("RansEvmKEpsilonHighReEpsilonCWD2D3N", mRansEvmKEpsilonHighReEpsilonCWD2D);
    KRATOS_REGISTER_ELEMENT("RansEvmKEpsilonHighReEpsilonCWD3D4N", mRansEvmKEpsilonHighReEpsilonCWD3D);
}
} // namespace Kratos.
