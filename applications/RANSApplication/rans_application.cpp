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
      mIncompressiblePotentialFlowVelocity2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mIncompressiblePotentialFlowVelocity3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(Element::GeometryType::PointsArrayType(4)))),
      mIncompressiblePotentialFlowPressure2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mIncompressiblePotentialFlowPressure3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(Element::GeometryType::PointsArrayType(4)))),
      mIncompressiblePotentialFlowVelocityCondition2D2N(0,Element::GeometryType::Pointer(new Line2D2<Node<3>>(Element::GeometryType::PointsArrayType(2)))),
      mIncompressiblePotentialFlowVelocityCondition3D3N(0,Element::GeometryType::Pointer(new Triangle3D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mIncompressiblePotentialFlowPressureCondition2D2N(0,Element::GeometryType::Pointer(new Line2D2<Node<3>>(Element::GeometryType::PointsArrayType(2)))),
      mIncompressiblePotentialFlowPressureCondition3D3N(0,Element::GeometryType::Pointer(new Triangle3D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mRansFractionalStep2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mRansFractionalStep3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(Element::GeometryType::PointsArrayType(4)))),
      mFSHighReKWallCondition2D2N(0,Element::GeometryType::Pointer(new Line2D2<Node<3>>(Element::GeometryType::PointsArrayType(2)))),
      mFSHighReKWallCondition3D3N(0,Element::GeometryType::Pointer(new Triangle3D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mRansEvmKEpsilonLowReK2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mRansEvmKEpsilonLowReK3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(Element::GeometryType::PointsArrayType(4)))),
      mRansEvmKEpsilonLowReEpsilon2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mRansEvmKEpsilonLowReEpsilon3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(Element::GeometryType::PointsArrayType(4)))),
      mRansEvmKEpsilonK2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3>>( Element::GeometryType::PointsArrayType(3)))),
      mRansEvmKEpsilonK3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>( Element::GeometryType::PointsArrayType(4)))),
      mRansEvmKEpsilonEpsilon2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3>>( Element::GeometryType::PointsArrayType(3)))),
      mRansEvmKEpsilonEpsilon3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>( Element::GeometryType::PointsArrayType(4)))),
      mRansEvmKEpsilonEpsilonWall2D2N(0,Element::GeometryType::Pointer(new Line2D2<Node<3>>(Element::GeometryType::PointsArrayType(2)))),
      mRansEvmKEpsilonEpsilonWall3D3N(0,Element::GeometryType::Pointer(new Triangle3D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mRansEvmKEpsilonVmsMonolithicWall2D2N(0,Element::GeometryType::Pointer(new Line2D2<Node<3>>(Element::GeometryType::PointsArrayType(2)))),
      mRansEvmKEpsilonVmsMonolithicWall3D3N(0,Element::GeometryType::Pointer(new Triangle3D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mRansEvmEpsilonAdjoint2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mRansEvmEpsilonAdjoint3D4N(0,Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(Element::GeometryType::PointsArrayType(4)))),
      mRansEvmKAdjoint2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mRansEvmKAdjoint3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(Element::GeometryType::PointsArrayType(4)))),
      mRansEvmKEpsilonVMSAdjoint2D3N(0,Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mRansEvmKEpsilonVMSAdjoint3D4N(0,Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(Element::GeometryType::PointsArrayType(4)))),
      mRansEvmMonolithicKEpsilonVMSAdjoint2D(0,Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mRansEvmMonolithicKEpsilonVMSAdjoint3D(0,Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(Element::GeometryType::PointsArrayType(4)))),
      mRansEvmEpsilonAdjointWallCondition2D2N(0,Element::GeometryType::Pointer(new Line2D2<Node<3>>(Element::GeometryType::PointsArrayType(2)))),
      mRansEvmEpsilonAdjointWallCondition3D3N(0,Element::GeometryType::Pointer(new Triangle3D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mRansEvmVmsMonolithicAdjointWallCondition2D2N(0,Element::GeometryType::Pointer(new Line2D2<Node<3>>(Element::GeometryType::PointsArrayType(2)))),
      mRansEvmVmsMonolithicAdjointWallCondition3D3N(0,Element::GeometryType::Pointer(new Triangle3D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mRansEvmMonolithicKEpsilonVMSAdjointWallCondition2D2N(0,Element::GeometryType::Pointer(new Line2D2<Node<3>>(Element::GeometryType::PointsArrayType(2)))),
      mRansEvmMonolithicKEpsilonVMSAdjointWallCondition3D3N(0,Element::GeometryType::Pointer(new Triangle3D3<Node<3>>(Element::GeometryType::PointsArrayType(3))))
{
}

void KratosRANSApplication::Register()
{
    KRATOS_INFO("") << "Initializing KratosRANSApplication..." << std::endl;

    KRATOS_REGISTER_VARIABLE(TURBULENT_KINETIC_ENERGY)
    KRATOS_REGISTER_VARIABLE(TURBULENT_ENERGY_DISSIPATION_RATE)
    KRATOS_REGISTER_VARIABLE(TURBULENT_KINETIC_ENERGY_RATE)
    KRATOS_REGISTER_VARIABLE(TURBULENT_ENERGY_DISSIPATION_RATE_2)
    KRATOS_REGISTER_VARIABLE(IS_CO_SOLVING_PROCESS_ACTIVE)
    KRATOS_REGISTER_VARIABLE(RANS_Y_PLUS)
    KRATOS_REGISTER_VARIABLE(RANS_AUXILIARY_VARIABLE_1)
    KRATOS_REGISTER_VARIABLE(RANS_AUXILIARY_VARIABLE_2)
    KRATOS_REGISTER_VARIABLE(WALL_SMOOTHNESS_BETA)
    KRATOS_REGISTER_VARIABLE(WALL_VON_KARMAN)
    KRATOS_REGISTER_VARIABLE(TURBULENCE_RANS_C_MU)
    KRATOS_REGISTER_VARIABLE(TURBULENCE_RANS_C1)
    KRATOS_REGISTER_VARIABLE(TURBULENCE_RANS_C2)
    KRATOS_REGISTER_VARIABLE(TURBULENT_KINETIC_ENERGY_SIGMA)
    KRATOS_REGISTER_VARIABLE(TURBULENT_ENERGY_DISSIPATION_RATE_SIGMA)
    KRATOS_REGISTER_VARIABLE(NUMBER_OF_NEIGHBOUR_CONDITIONS)
    KRATOS_REGISTER_VARIABLE(COUPLING_ITERATION)
    KRATOS_REGISTER_VARIABLE(ANALYSIS_STEPS)
    KRATOS_REGISTER_VARIABLE(VELOCITY_POTENTIAL)
    KRATOS_REGISTER_VARIABLE(PRESSURE_POTENTIAL)
    KRATOS_REGISTER_VARIABLE(RANS_Y_PLUS_LIMIT)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(FRICTION_VELOCITY)
    KRATOS_REGISTER_VARIABLE(WALL_MODEL_PART_NAMES)

    KRATOS_REGISTER_VARIABLE(RANS_STABILIZATION_DISCRETE_UPWIND_OPERATOR_COEFFICIENT)
    KRATOS_REGISTER_VARIABLE(RANS_STABILIZATION_DIAGONAL_POSITIVITY_PRESERVING_COEFFICIENT)

    KRATOS_REGISTER_VARIABLE(RANS_IS_INLET)
    KRATOS_REGISTER_VARIABLE(RANS_IS_OUTLET)
    KRATOS_REGISTER_VARIABLE(RANS_IS_WALL)

    // Register adjoint variables
    KRATOS_REGISTER_VARIABLE(RANS_NUT_SCALAR_PARTIAL_DERIVATIVES)
    KRATOS_REGISTER_VARIABLE(RANS_Y_PLUS_VELOCITY_DERIVATIVES)
    KRATOS_REGISTER_VARIABLE(RANS_VELOCITY_PRESSURE_PARTIAL_DERIVATIVE)
    KRATOS_REGISTER_VARIABLE(RANS_PRESSURE_PARTIAL_DERIVATIVE)
    KRATOS_REGISTER_VARIABLE(RANS_TURBULENT_KINETIC_ENERGY_PARTIAL_DERIVATIVE)
    KRATOS_REGISTER_VARIABLE(RANS_TURBULENT_ENERGY_DISSIPATION_RATE_PARTIAL_DERIVATIVE)
    KRATOS_REGISTER_VARIABLE(RANS_ACCELERATION_PARTIAL_DERIVATIVE)
    KRATOS_REGISTER_VARIABLE(RANS_TURBULENT_KINETIC_ENERGY_RATE_PARTIAL_DERIVATIVE)
    KRATOS_REGISTER_VARIABLE(RANS_TURBULENT_ENERGY_DISSIPATION_RATE_2_PARTIAL_DERIVATIVE)

    KRATOS_REGISTER_VARIABLE(RANS_SCALAR_1_ADJOINT_1)
    KRATOS_REGISTER_VARIABLE(RANS_SCALAR_1_ADJOINT_2)
    KRATOS_REGISTER_VARIABLE(RANS_SCALAR_1_ADJOINT_3)
    KRATOS_REGISTER_VARIABLE(RANS_AUX_ADJOINT_SCALAR_1)

    KRATOS_REGISTER_VARIABLE(RANS_SCALAR_2_ADJOINT_1)
    KRATOS_REGISTER_VARIABLE(RANS_SCALAR_2_ADJOINT_2)
    KRATOS_REGISTER_VARIABLE(RANS_SCALAR_2_ADJOINT_3)
    KRATOS_REGISTER_VARIABLE(RANS_AUX_ADJOINT_SCALAR_2)

    // Register Elements
    KRATOS_REGISTER_ELEMENT("RansIncompressiblePotentialFlowVelocity2D3N", mIncompressiblePotentialFlowVelocity2D);
    KRATOS_REGISTER_ELEMENT("RansIncompressiblePotentialFlowVelocity3D4N", mIncompressiblePotentialFlowVelocity3D);
    KRATOS_REGISTER_ELEMENT("RansIncompressiblePotentialFlowPressure2D3N", mIncompressiblePotentialFlowPressure2D);
    KRATOS_REGISTER_ELEMENT("RansIncompressiblePotentialFlowPressure3D4N", mIncompressiblePotentialFlowPressure3D);

    KRATOS_REGISTER_ELEMENT("RansEvmKEpsilonLowReK2D3N", mRansEvmKEpsilonLowReK2D);
    KRATOS_REGISTER_ELEMENT("RansEvmKEpsilonLowReK3D4N", mRansEvmKEpsilonLowReK3D);
    KRATOS_REGISTER_ELEMENT("RansEvmKEpsilonLowReEpsilon2D3N", mRansEvmKEpsilonLowReEpsilon2D);
    KRATOS_REGISTER_ELEMENT("RansEvmKEpsilonLowReEpsilon3D4N", mRansEvmKEpsilonLowReEpsilon3D);

    KRATOS_REGISTER_ELEMENT("RansEvmKEpsilonK2D3N", mRansEvmKEpsilonK2D);
    KRATOS_REGISTER_ELEMENT("RansEvmKEpsilonK3D4N", mRansEvmKEpsilonK3D);
    KRATOS_REGISTER_ELEMENT("RansEvmKEpsilonEpsilon2D3N", mRansEvmKEpsilonEpsilon2D);
    KRATOS_REGISTER_ELEMENT("RansEvmKEpsilonEpsilon3D4N", mRansEvmKEpsilonEpsilon3D);

    KRATOS_REGISTER_CONDITION("RansIncompressiblePotentialFlowVelocity2D2N", mIncompressiblePotentialFlowVelocityCondition2D2N);
    KRATOS_REGISTER_CONDITION("RansIncompressiblePotentialFlowVelocity3D3N", mIncompressiblePotentialFlowVelocityCondition3D3N);
    KRATOS_REGISTER_CONDITION("RansIncompressiblePotentialFlowPressure2D2N", mIncompressiblePotentialFlowPressureCondition2D2N);
    KRATOS_REGISTER_CONDITION("RansIncompressiblePotentialFlowPressure3D3N", mIncompressiblePotentialFlowPressureCondition3D3N);

    KRATOS_REGISTER_ELEMENT("RansFractionalStep2D3N", mRansFractionalStep2D);
    KRATOS_REGISTER_ELEMENT("RansFractionalStep3D4N", mRansFractionalStep3D);

    KRATOS_REGISTER_CONDITION("RansFSHighReKWall2D2N", mFSHighReKWallCondition2D2N);
    KRATOS_REGISTER_CONDITION("RansFSHighReKWall3D3N", mFSHighReKWallCondition3D3N);

    KRATOS_REGISTER_CONDITION("RansEvmKEpsilonEpsilonWall2D2N", mRansEvmKEpsilonEpsilonWall2D2N);
    KRATOS_REGISTER_CONDITION("RansEvmKEpsilonEpsilonWall3D3N", mRansEvmKEpsilonEpsilonWall3D3N);

    KRATOS_REGISTER_CONDITION("RansEvmKEpsilonVmsMonolithicWall2D2N",
                              mRansEvmKEpsilonVmsMonolithicWall2D2N);
    KRATOS_REGISTER_CONDITION("RansEvmKEpsilonVmsMonolithicWall3D3N",
                              mRansEvmKEpsilonVmsMonolithicWall3D3N);

    // Registering adjoint elements
    KRATOS_REGISTER_ELEMENT("RansEvmEpsilonAdjoint2D3N", mRansEvmEpsilonAdjoint2D3N);
    KRATOS_REGISTER_ELEMENT("RansEvmEpsilonAdjoint3D4N", mRansEvmEpsilonAdjoint3D4N);

    KRATOS_REGISTER_ELEMENT("RansEvmKAdjoint2D3N", mRansEvmKAdjoint2D3N);
    KRATOS_REGISTER_ELEMENT("RansEvmKAdjoint3D4N", mRansEvmKAdjoint3D4N);

    KRATOS_REGISTER_ELEMENT("RansEvmKEpsilonVMSAdjoint2D3N", mRansEvmKEpsilonVMSAdjoint2D3N);
    KRATOS_REGISTER_ELEMENT("RansEvmKEpsilonVMSAdjoint3D4N", mRansEvmKEpsilonVMSAdjoint3D4N);

    KRATOS_REGISTER_ELEMENT("RansEvmMonolithicKEpsilonVMSAdjoint2D",
                            mRansEvmMonolithicKEpsilonVMSAdjoint2D);
    KRATOS_REGISTER_ELEMENT("RansEvmMonolithicKEpsilonVMSAdjoint3D",
                            mRansEvmMonolithicKEpsilonVMSAdjoint3D);

    // Registering adjoint conditions
    KRATOS_REGISTER_CONDITION("RansEvmEpsilonAdjointWallCondition2D2N",
                              mRansEvmEpsilonAdjointWallCondition2D2N);
    KRATOS_REGISTER_CONDITION("RansEvmEpsilonAdjointWallCondition3D3N",
                              mRansEvmEpsilonAdjointWallCondition3D3N);

    KRATOS_REGISTER_CONDITION("RansEvmVmsMonolithicAdjointWallCondition2D2N",
                              mRansEvmVmsMonolithicAdjointWallCondition2D2N);
    KRATOS_REGISTER_CONDITION("RansEvmVmsMonolithicAdjointWallCondition3D3N",
                              mRansEvmVmsMonolithicAdjointWallCondition3D3N);

    KRATOS_REGISTER_CONDITION(
        "RansEvmMonolithicKEpsilonVMSAdjointWallCondition2D2N",
        mRansEvmMonolithicKEpsilonVMSAdjointWallCondition2D2N);
    KRATOS_REGISTER_CONDITION(
        "RansEvmMonolithicKEpsilonVMSAdjointWallCondition3D3N",
        mRansEvmMonolithicKEpsilonVMSAdjointWallCondition3D3N);
}
} // namespace Kratos.
