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
      mIncompressiblePotentialFlowPressure3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(Element::GeometryType::PointsArrayType(4)))),
      // k-epsilon turbulence model elements
      mRansKEpsilonKAFC2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mRansKEpsilonKAFC3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(Element::GeometryType::PointsArrayType(4)))),
      mRansKEpsilonEpsilonAFC2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mRansKEpsilonEpsilonAFC3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(Element::GeometryType::PointsArrayType(4)))),
      mRansKEpsilonKRFC2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mRansKEpsilonKRFC3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(Element::GeometryType::PointsArrayType(4)))),
      mRansKEpsilonEpsilonRFC2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mRansKEpsilonEpsilonRFC3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(Element::GeometryType::PointsArrayType(4)))),
      mRansKEpsilonKCWD2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mRansKEpsilonKCWD3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(Element::GeometryType::PointsArrayType(4)))),
      mRansKEpsilonEpsilonCWD2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mRansKEpsilonEpsilonCWD3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(Element::GeometryType::PointsArrayType(4))))
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

    // residual based flux corrected stabilization variables
    KRATOS_REGISTER_VARIABLE( RANS_STABILIZATION_DISCRETE_UPWIND_OPERATOR_COEFFICIENT )
    KRATOS_REGISTER_VARIABLE( RANS_STABILIZATION_DIAGONAL_POSITIVITY_PRESERVING_COEFFICIENT )

    // algebraic flux corrected stabilization variables
    KRATOS_REGISTER_VARIABLE( AFC_POSITIVE_ANTI_DIFFUSIVE_FLUX )
    KRATOS_REGISTER_VARIABLE( AFC_NEGATIVE_ANTI_DIFFUSIVE_FLUX )
    KRATOS_REGISTER_VARIABLE( AFC_POSITIVE_ANTI_DIFFUSIVE_FLUX_LIMIT )
    KRATOS_REGISTER_VARIABLE( AFC_NEGATIVE_ANTI_DIFFUSIVE_FLUX_LIMIT )

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
    // registering incompressible potential flow elements
    KRATOS_REGISTER_ELEMENT("RansIncompressiblePotentialFlowVelocity2D3N", mIncompressiblePotentialFlowVelocity2D);
    KRATOS_REGISTER_ELEMENT("RansIncompressiblePotentialFlowVelocity3D4N", mIncompressiblePotentialFlowVelocity3D);
    KRATOS_REGISTER_ELEMENT("RansIncompressiblePotentialFlowPressure2D3N", mIncompressiblePotentialFlowPressure2D);
    KRATOS_REGISTER_ELEMENT("RansIncompressiblePotentialFlowPressure3D4N", mIncompressiblePotentialFlowPressure3D);

    // registering k-epsilon algebraic flux correction based elements
    KRATOS_REGISTER_ELEMENT("RansKEpsilonKAFC2D3N", mRansKEpsilonKAFC2D);
    KRATOS_REGISTER_ELEMENT("RansKEpsilonKAFC3D4N", mRansKEpsilonKAFC3D);
    KRATOS_REGISTER_ELEMENT("RansKEpsilonEpsilonAFC2D3N", mRansKEpsilonEpsilonAFC2D);
    KRATOS_REGISTER_ELEMENT("RansKEpsilonEpsilonAFC3D4N", mRansKEpsilonEpsilonAFC3D);

    // registering k-epsilon residual fc based elements
    KRATOS_REGISTER_ELEMENT("RansKEpsilonKRFC2D3N", mRansKEpsilonKRFC2D);
    KRATOS_REGISTER_ELEMENT("RansKEpsilonKRFC3D4N", mRansKEpsilonKRFC3D);
    KRATOS_REGISTER_ELEMENT("RansKEpsilonEpsilonRFC2D3N", mRansKEpsilonEpsilonRFC2D);
    KRATOS_REGISTER_ELEMENT("RansKEpsilonEpsilonRFC3D4N", mRansKEpsilonEpsilonRFC3D);

    // registering k-epsilon cross wind stabilized elements
    KRATOS_REGISTER_ELEMENT("RansKEpsilonKCWD2D3N", mRansKEpsilonKCWD2D);
    KRATOS_REGISTER_ELEMENT("RansKEpsilonKCWD3D4N", mRansKEpsilonKCWD3D);
    KRATOS_REGISTER_ELEMENT("RansKEpsilonEpsilonCWD2D3N", mRansKEpsilonEpsilonCWD2D);
    KRATOS_REGISTER_ELEMENT("RansKEpsilonEpsilonCWD3D4N", mRansKEpsilonEpsilonCWD3D);
}
} // namespace Kratos.
