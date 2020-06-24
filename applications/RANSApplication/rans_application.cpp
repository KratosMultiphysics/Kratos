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
      mRansEvmKOmegaOmegaCWD3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(Element::GeometryType::PointsArrayType(4)))),
      // k-omega-sst turbulence model elements
      mRansEvmKOmegaSSTKAFC2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mRansEvmKOmegaSSTKAFC3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(Element::GeometryType::PointsArrayType(4)))),
      mRansEvmKOmegaSSTOmegaAFC2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mRansEvmKOmegaSSTOmegaAFC3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(Element::GeometryType::PointsArrayType(4)))),
      mRansEvmKOmegaSSTKRFC2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mRansEvmKOmegaSSTKRFC3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(Element::GeometryType::PointsArrayType(4)))),
      mRansEvmKOmegaSSTOmegaRFC2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mRansEvmKOmegaSSTOmegaRFC3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(Element::GeometryType::PointsArrayType(4)))),
      mRansEvmKOmegaSSTKCWD2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mRansEvmKOmegaSSTKCWD3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(Element::GeometryType::PointsArrayType(4)))),
      mRansEvmKOmegaSSTOmegaCWD2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mRansEvmKOmegaSSTOmegaCWD3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(Element::GeometryType::PointsArrayType(4)))),
      // k-epsilon turbulence model conditions
      mRansEvmKEpsilonHighReEpsilonKBasedWall2D2N(0,Condition::GeometryType::Pointer(new Line2D2<Node<3>>(Condition::GeometryType::PointsArrayType(2)))),
      mRansEvmKEpsilonHighReEpsilonKBasedWall3D3N(0,Condition::GeometryType::Pointer(new Triangle3D3<Node<3>>(Condition::GeometryType::PointsArrayType(3)))),
      mRansEvmKEpsilonHighReEpsilonUBasedWall2D2N(0,Condition::GeometryType::Pointer(new Line2D2<Node<3>>(Condition::GeometryType::PointsArrayType(2)))),
      mRansEvmKEpsilonHighReEpsilonUBasedWall3D3N(0,Condition::GeometryType::Pointer(new Triangle3D3<Node<3>>(Condition::GeometryType::PointsArrayType(3)))),
      // k-omega turbulence model conditions
      mRansEvmKOmegaOmegaKBasedWall2D2N(0,Condition::GeometryType::Pointer(new Line2D2<Node<3>>(Condition::GeometryType::PointsArrayType(2)))),
      mRansEvmKOmegaOmegaKBasedWall3D3N(0,Condition::GeometryType::Pointer(new Triangle3D3<Node<3>>(Condition::GeometryType::PointsArrayType(3)))),
      mRansEvmKOmegaOmegaUBasedWall2D2N(0,Condition::GeometryType::Pointer(new Line2D2<Node<3>>(Condition::GeometryType::PointsArrayType(2)))),
      mRansEvmKOmegaOmegaUBasedWall3D3N(0,Condition::GeometryType::Pointer(new Triangle3D3<Node<3>>(Condition::GeometryType::PointsArrayType(3))))
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

    // k-omega-sst turbulence modelling specific additional variables
    KRATOS_REGISTER_VARIABLE( TURBULENT_KINETIC_ENERGY_SIGMA_1 )
    KRATOS_REGISTER_VARIABLE( TURBULENT_KINETIC_ENERGY_SIGMA_2 )
    KRATOS_REGISTER_VARIABLE( TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_SIGMA_1 )
    KRATOS_REGISTER_VARIABLE( TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_SIGMA_2 )
    KRATOS_REGISTER_VARIABLE( TURBULENCE_RANS_BETA_1 )
    KRATOS_REGISTER_VARIABLE( TURBULENCE_RANS_BETA_2 )
    KRATOS_REGISTER_VARIABLE( WALL_VON_KARMAN )

    // wall function condition specific additional variables
    KRATOS_REGISTER_VARIABLE( RANS_Y_PLUS )
    KRATOS_REGISTER_VARIABLE( RANS_LINEAR_LOG_LAW_Y_PLUS_LIMIT )
    KRATOS_REGISTER_VARIABLE( WALL_SMOOTHNESS_BETA )
    KRATOS_REGISTER_VARIABLE( RANS_IS_WALL )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( FRICTION_VELOCITY )

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

    // registering k-omega-sst algebraic flux correction based elements
    KRATOS_REGISTER_ELEMENT("RansEvmKOmegaSSTKAFC2D3N", mRansEvmKOmegaSSTKAFC2D);
    KRATOS_REGISTER_ELEMENT("RansEvmKOmegaSSTKAFC3D4N", mRansEvmKOmegaSSTKAFC3D);
    KRATOS_REGISTER_ELEMENT("RansEvmKOmegaSSTOmegaAFC2D3N", mRansEvmKOmegaSSTOmegaAFC2D);
    KRATOS_REGISTER_ELEMENT("RansEvmKOmegaSSTOmegaAFC3D4N", mRansEvmKOmegaSSTOmegaAFC3D);

    // registering k-omega-sst residual fc based elements
    KRATOS_REGISTER_ELEMENT("RansEvmKOmegaSSTKRFC2D3N", mRansEvmKOmegaSSTKRFC2D);
    KRATOS_REGISTER_ELEMENT("RansEvmKOmegaSSTKRFC3D4N", mRansEvmKOmegaSSTKRFC3D);
    KRATOS_REGISTER_ELEMENT("RansEvmKOmegaSSTOmegaRFC2D3N", mRansEvmKOmegaSSTOmegaRFC2D);
    KRATOS_REGISTER_ELEMENT("RansEvmKOmegaSSTOmegaRFC3D4N", mRansEvmKOmegaSSTOmegaRFC3D);

    // registering k-omega-sst cross wind stabilized elements
    KRATOS_REGISTER_ELEMENT("RansEvmKOmegaSSTKCWD2D3N", mRansEvmKOmegaSSTKCWD2D);
    KRATOS_REGISTER_ELEMENT("RansEvmKOmegaSSTKCWD3D4N", mRansEvmKOmegaSSTKCWD3D);
    KRATOS_REGISTER_ELEMENT("RansEvmKOmegaSSTOmegaCWD2D3N", mRansEvmKOmegaSSTOmegaCWD2D);
    KRATOS_REGISTER_ELEMENT("RansEvmKOmegaSSTOmegaCWD3D4N", mRansEvmKOmegaSSTOmegaCWD3D);

    // registering conditions
    // registering k-epsilon conditions
    KRATOS_REGISTER_CONDITION("RansEvmKEpsilonHighReEpsilonKBasedWall2D2N", mRansEvmKEpsilonHighReEpsilonKBasedWall2D2N);
    KRATOS_REGISTER_CONDITION("RansEvmKEpsilonHighReEpsilonKBasedWall3D3N", mRansEvmKEpsilonHighReEpsilonKBasedWall3D3N);

    KRATOS_REGISTER_CONDITION("RansEvmKEpsilonHighReEpsilonUBasedWall2D2N", mRansEvmKEpsilonHighReEpsilonUBasedWall2D2N);
    KRATOS_REGISTER_CONDITION("RansEvmKEpsilonHighReEpsilonUBasedWall3D3N", mRansEvmKEpsilonHighReEpsilonUBasedWall3D3N);

    // registering k-omega conditions
    KRATOS_REGISTER_CONDITION("RansEvmKOmegaOmegaKBasedWall2D2N", mRansEvmKOmegaOmegaKBasedWall2D2N);
    KRATOS_REGISTER_CONDITION("RansEvmKOmegaOmegaKBasedWall3D3N", mRansEvmKOmegaOmegaKBasedWall3D3N);

    KRATOS_REGISTER_CONDITION("RansEvmKOmegaOmegaUBasedWall2D2N", mRansEvmKOmegaOmegaUBasedWall2D2N);
    KRATOS_REGISTER_CONDITION("RansEvmKOmegaOmegaUBasedWall3D3N", mRansEvmKOmegaOmegaUBasedWall3D3N);
}
} // namespace Kratos.
