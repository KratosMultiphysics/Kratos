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
      mRansEvmKEpsilonHighReEpsilonCWD3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(Element::GeometryType::PointsArrayType(4)))),
      // k-omega turbulence model elements
      mRansEvmKOmegaKAFC2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mRansEvmKOmegaKAFC3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(Element::GeometryType::PointsArrayType(4)))),
      mRansEvmKOmegaOmegaAFC2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mRansEvmKOmegaOmegaAFC3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(Element::GeometryType::PointsArrayType(4)))),
      mRansEvmKOmegaKRFC2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mRansEvmKOmegaKRFC3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(Element::GeometryType::PointsArrayType(4)))),
      mRansEvmKOmegaOmegaRFC2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mRansEvmKOmegaOmegaRFC3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(Element::GeometryType::PointsArrayType(4)))),
      mRansEvmKOmegaKCWD2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mRansEvmKOmegaKCWD3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(Element::GeometryType::PointsArrayType(4)))),
      mRansEvmKOmegaOmegaCWD2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mRansEvmKOmegaOmegaCWD3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(Element::GeometryType::PointsArrayType(4))))
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

    // k-omega turbulence modelling specific additional variables
    KRATOS_REGISTER_VARIABLE( TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE )
    KRATOS_REGISTER_VARIABLE( TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_2 )
    KRATOS_REGISTER_VARIABLE( TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_SIGMA )
    KRATOS_REGISTER_VARIABLE( TURBULENCE_RANS_BETA )
    KRATOS_REGISTER_VARIABLE( TURBULENCE_RANS_GAMMA )

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

    // registering k-omega algebraic flux correction based elements
    KRATOS_REGISTER_ELEMENT("RansEvmKOmegaKAFC2D3N", mRansEvmKOmegaKAFC2D);
    KRATOS_REGISTER_ELEMENT("RansEvmKOmegaKAFC3D4N", mRansEvmKOmegaKAFC3D);
    KRATOS_REGISTER_ELEMENT("RansEvmKOmegaOmegaAFC2D3N", mRansEvmKOmegaOmegaAFC2D);
    KRATOS_REGISTER_ELEMENT("RansEvmKOmegaOmegaAFC3D4N", mRansEvmKOmegaOmegaAFC3D);

    // registering k-omega residual fc based elements
    KRATOS_REGISTER_ELEMENT("RansEvmKOmegaKRFC2D3N", mRansEvmKOmegaKRFC2D);
    KRATOS_REGISTER_ELEMENT("RansEvmKOmegaKRFC3D4N", mRansEvmKOmegaKRFC3D);
    KRATOS_REGISTER_ELEMENT("RansEvmKOmegaOmegaRFC2D3N", mRansEvmKOmegaOmegaRFC2D);
    KRATOS_REGISTER_ELEMENT("RansEvmKOmegaOmegaRFC3D4N", mRansEvmKOmegaOmegaRFC3D);

    // registering k-omega cross wind stabilized elements
    KRATOS_REGISTER_ELEMENT("RansEvmKOmegaKCWD2D3N", mRansEvmKOmegaKCWD2D);
    KRATOS_REGISTER_ELEMENT("RansEvmKOmegaKCWD3D4N", mRansEvmKOmegaKCWD3D);
    KRATOS_REGISTER_ELEMENT("RansEvmKOmegaOmegaCWD2D3N", mRansEvmKOmegaOmegaCWD2D);
    KRATOS_REGISTER_ELEMENT("RansEvmKOmegaOmegaCWD3D4N", mRansEvmKOmegaOmegaCWD3D);
}
} // namespace Kratos.
