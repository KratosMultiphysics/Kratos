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
      // incompressible potential flow conditions
      mIncompressiblePotentialFlowVelocityCondition2D2N(0,Condition::GeometryType::Pointer(new Line2D2<Node<3>>(Condition::GeometryType::PointsArrayType(2)))),
      mIncompressiblePotentialFlowVelocityCondition3D3N(0,Condition::GeometryType::Pointer(new Triangle3D3<Node<3>>(Condition::GeometryType::PointsArrayType(3)))),
      mIncompressiblePotentialFlowPressureCondition2D2N(0,Condition::GeometryType::Pointer(new Line2D2<Node<3>>(Condition::GeometryType::PointsArrayType(2)))),
      mIncompressiblePotentialFlowPressureCondition3D3N(0,Condition::GeometryType::Pointer(new Triangle3D3<Node<3>>(Condition::GeometryType::PointsArrayType(3)))),
      // fractional step elemnets
      mRansFractionalStep2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mRansFractionalStep3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(Element::GeometryType::PointsArrayType(4)))),
      // fractional step conditions
      mFSHighReKWallCondition2D2N(0,Condition::GeometryType::Pointer(new Line2D2<Node<3>>(Condition::GeometryType::PointsArrayType(2)))),
      mFSHighReKWallCondition3D3N(0,Condition::GeometryType::Pointer(new Triangle3D3<Node<3>>(Condition::GeometryType::PointsArrayType(3)))),
      // monolithic conditions
      mRansEvmKEpsilonVmsMonolithicWall2D2N(0,Condition::GeometryType::Pointer(new Line2D2<Node<3>>(Condition::GeometryType::PointsArrayType(2)))),
      mRansEvmKEpsilonVmsMonolithicWall3D3N(0,Condition::GeometryType::Pointer(new Triangle3D3<Node<3>>(Condition::GeometryType::PointsArrayType(3)))),
      mRansVMSMonolithicKBasedWallCondition2D2N(0,Condition::GeometryType::Pointer(new Line2D2<Node<3>>(Condition::GeometryType::PointsArrayType(2)))),
      mRansVMSMonolithicKBasedWallCondition3D3N(0,Condition::GeometryType::Pointer(new Triangle3D3<Node<3>>(Condition::GeometryType::PointsArrayType(3)))),
      // k-epsilon turbulence model elements
      mRansEvmKEpsilonKAFC2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mRansEvmKEpsilonKAFC3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(Element::GeometryType::PointsArrayType(4)))),
      mRansEvmKEpsilonEpsilonAFC2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mRansEvmKEpsilonEpsilonAFC3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(Element::GeometryType::PointsArrayType(4)))),
      mRansEvmKEpsilonKResidualBasedFC2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mRansEvmKEpsilonKResidualBasedFC3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(Element::GeometryType::PointsArrayType(4)))),
      mRansEvmKEpsilonEpsilonResidualBasedFC2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mRansEvmKEpsilonEpsilonResidualBasedFC3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(Element::GeometryType::PointsArrayType(4)))),
      mRansEvmKEpsilonKCrossWindStabilized2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mRansEvmKEpsilonKCrossWindStabilized3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(Element::GeometryType::PointsArrayType(4)))),
      mRansEvmKEpsilonEpsilonCrossWindStabilized2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mRansEvmKEpsilonEpsilonCrossWindStabilized3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(Element::GeometryType::PointsArrayType(4)))),
      // k-epsilon turbulence model conditions
      mRansEvmKEpsilonEpsilonKBasedWallCondition2D2N(0,Condition::GeometryType::Pointer(new Line2D2<Node<3>>(Condition::GeometryType::PointsArrayType(2)))),
      mRansEvmKEpsilonEpsilonKBasedWallCondition3D3N(0,Condition::GeometryType::Pointer(new Triangle3D3<Node<3>>(Condition::GeometryType::PointsArrayType(3)))),
      mRansEvmKEpsilonEpsilonUBasedWallCondition2D2N(0,Condition::GeometryType::Pointer(new Line2D2<Node<3>>(Condition::GeometryType::PointsArrayType(2)))),
      mRansEvmKEpsilonEpsilonUBasedWallCondition3D3N(0,Condition::GeometryType::Pointer(new Triangle3D3<Node<3>>(Condition::GeometryType::PointsArrayType(3)))),
      // k-omega turbulence model elements
      mRansEvmKOmegaKAFC2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mRansEvmKOmegaKAFC3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(Element::GeometryType::PointsArrayType(4)))),
      mRansEvmKOmegaOmegaAFC2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mRansEvmKOmegaOmegaAFC3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(Element::GeometryType::PointsArrayType(4)))),
      mRansEvmKOmegaKResidualBasedFC2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mRansEvmKOmegaKResidualBasedFC3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(Element::GeometryType::PointsArrayType(4)))),
      mRansEvmKOmegaOmegaResidualBasedFC2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mRansEvmKOmegaOmegaResidualBasedFC3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(Element::GeometryType::PointsArrayType(4)))),
      mRansEvmKOmegaKCrossWindStabilized2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mRansEvmKOmegaKCrossWindStabilized3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(Element::GeometryType::PointsArrayType(4)))),
      mRansEvmKOmegaOmegaCrossWindStabilized2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mRansEvmKOmegaOmegaCrossWindStabilized3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(Element::GeometryType::PointsArrayType(4)))),
      // k-omega turbulence model conditions
      mRansEvmKOmegaOmegaKBasedWallCondition2D2N(0,Condition::GeometryType::Pointer(new Line2D2<Node<3>>(Condition::GeometryType::PointsArrayType(2)))),
      mRansEvmKOmegaOmegaKBasedWallCondition3D3N(0,Condition::GeometryType::Pointer(new Triangle3D3<Node<3>>(Condition::GeometryType::PointsArrayType(3)))),
      mRansEvmKOmegaOmegaUBasedWallCondition2D2N(0,Condition::GeometryType::Pointer(new Line2D2<Node<3>>(Condition::GeometryType::PointsArrayType(2)))),
      mRansEvmKOmegaOmegaUBasedWallCondition3D3N(0,Condition::GeometryType::Pointer(new Triangle3D3<Node<3>>(Condition::GeometryType::PointsArrayType(3)))),
      // k-omega-sst turbulence model elements
      mRansEvmKOmegaSSTKAFC2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mRansEvmKOmegaSSTKAFC3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(Element::GeometryType::PointsArrayType(4)))),
      mRansEvmKOmegaSSTOmegaAFC2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mRansEvmKOmegaSSTOmegaAFC3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(Element::GeometryType::PointsArrayType(4)))),
      mRansEvmKOmegaSSTKResidualBasedFC2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mRansEvmKOmegaSSTKResidualBasedFC3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(Element::GeometryType::PointsArrayType(4)))),
      mRansEvmKOmegaSSTOmegaResidualBasedFC2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mRansEvmKOmegaSSTOmegaResidualBasedFC3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(Element::GeometryType::PointsArrayType(4)))),
      mRansEvmKOmegaSSTKCrossWindStabilized2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mRansEvmKOmegaSSTKCrossWindStabilized3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(Element::GeometryType::PointsArrayType(4)))),
      mRansEvmKOmegaSSTOmegaCrossWindStabilized2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mRansEvmKOmegaSSTOmegaCrossWindStabilized3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(Element::GeometryType::PointsArrayType(4)))),
      // old elements
      mRansEvmKEpsilonLowReK2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mRansEvmKEpsilonLowReK3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(Element::GeometryType::PointsArrayType(4)))),
      mRansEvmKEpsilonLowReEpsilon2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mRansEvmKEpsilonLowReEpsilon3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(Element::GeometryType::PointsArrayType(4)))),
      mRansEvmKEpsilonK2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3>>( Element::GeometryType::PointsArrayType(3)))),
      mRansEvmKEpsilonK3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>( Element::GeometryType::PointsArrayType(4)))),
      mRansEvmKEpsilonEpsilon2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3>>( Element::GeometryType::PointsArrayType(3)))),
      mRansEvmKEpsilonEpsilon3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>( Element::GeometryType::PointsArrayType(4)))),
      mRansEvmKEpsilonEpsilonVelocityBasedWallCondition2D2N(0,Condition::GeometryType::Pointer(new Line2D2<Node<3>>(Condition::GeometryType::PointsArrayType(2)))),
      mRansEvmKEpsilonEpsilonVelocityBasedWallCondition3D3N(0,Condition::GeometryType::Pointer(new Triangle3D3<Node<3>>(Condition::GeometryType::PointsArrayType(3)))),
      mRansEvmKEpsilonEpsilonKBasedLHSWallCondition2D2N(0,Condition::GeometryType::Pointer(new Line2D2<Node<3>>(Condition::GeometryType::PointsArrayType(2)))),
      mRansEvmKEpsilonEpsilonKBasedLHSWallCondition3D3N(0,Condition::GeometryType::Pointer(new Triangle3D3<Node<3>>(Condition::GeometryType::PointsArrayType(3)))),
      mRansEvmKEpsilonEpsilonKBasedRHSWallCondition2D2N(0,Condition::GeometryType::Pointer(new Line2D2<Node<3>>(Condition::GeometryType::PointsArrayType(2)))),
      mRansEvmKEpsilonEpsilonKBasedRHSWallCondition3D3N(0,Condition::GeometryType::Pointer(new Triangle3D3<Node<3>>(Condition::GeometryType::PointsArrayType(3)))),
      mRansEvmKEpsilonEpsilonWall2D2N(0,Condition::GeometryType::Pointer(new Line2D2<Node<3>>(Condition::GeometryType::PointsArrayType(2)))),
      mRansEvmKEpsilonEpsilonWall3D3N(0,Condition::GeometryType::Pointer(new Triangle3D3<Node<3>>(Condition::GeometryType::PointsArrayType(3)))),
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
      mRansEvmVmsMonolithicAdjointWallCondition2D2N(0,Condition::GeometryType::Pointer(new Line2D2<Node<3>>(Condition::GeometryType::PointsArrayType(2)))),
      mRansEvmVmsMonolithicAdjointWallCondition3D3N(0,Condition::GeometryType::Pointer(new Triangle3D3<Node<3>>(Condition::GeometryType::PointsArrayType(3)))),
      mRansEvmMonolithicKEpsilonVMSAdjointWallCondition2D2N(0,Condition::GeometryType::Pointer(new Line2D2<Node<3>>(Condition::GeometryType::PointsArrayType(2)))),
      mRansEvmMonolithicKEpsilonVMSAdjointWallCondition3D3N(0,Condition::GeometryType::Pointer(new Triangle3D3<Node<3>>(Condition::GeometryType::PointsArrayType(3))))
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
    KRATOS_REGISTER_VARIABLE(RANS_Y_PLUS_LOWER_LIMIT)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(FRICTION_VELOCITY)
    KRATOS_REGISTER_VARIABLE(WALL_MODEL_PART_NAME)

    KRATOS_REGISTER_VARIABLE(TURBULENCE_RANS_BETA_ZERO)
    KRATOS_REGISTER_VARIABLE(TURBULENCE_RANS_BETA)
    KRATOS_REGISTER_VARIABLE(TURBULENCE_RANS_GAMMA)
    KRATOS_REGISTER_VARIABLE(TURBULENCE_RANS_K_C1)
    KRATOS_REGISTER_VARIABLE(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_SIGMA)
    KRATOS_REGISTER_VARIABLE(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE)
    KRATOS_REGISTER_VARIABLE(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_2)

    KRATOS_REGISTER_VARIABLE(TURBULENCE_RANS_A1)
    KRATOS_REGISTER_VARIABLE(TURBULENCE_RANS_BETA_1)
    KRATOS_REGISTER_VARIABLE(TURBULENCE_RANS_BETA_2)
    KRATOS_REGISTER_VARIABLE(TURBULENT_KINETIC_ENERGY_SIGMA_1)
    KRATOS_REGISTER_VARIABLE(TURBULENT_KINETIC_ENERGY_SIGMA_2)
    KRATOS_REGISTER_VARIABLE(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_SIGMA_1)
    KRATOS_REGISTER_VARIABLE(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_SIGMA_2)

    KRATOS_REGISTER_VARIABLE(RANS_STABILIZATION_DISCRETE_UPWIND_OPERATOR_COEFFICIENT)
    KRATOS_REGISTER_VARIABLE(RANS_STABILIZATION_DIAGONAL_POSITIVITY_PRESERVING_COEFFICIENT)

    KRATOS_REGISTER_VARIABLE(RANS_IS_INLET)
    KRATOS_REGISTER_VARIABLE(RANS_IS_OUTLET)
    KRATOS_REGISTER_VARIABLE(RANS_IS_WALL)

    KRATOS_REGISTER_VARIABLE(AFC_POSITIVE_ANTI_DIFFUSIVE_FLUX)
    KRATOS_REGISTER_VARIABLE(AFC_NEGATIVE_ANTI_DIFFUSIVE_FLUX)
    KRATOS_REGISTER_VARIABLE(AFC_POSITIVE_ANTI_DIFFUSIVE_FLUX_LIMIT)
    KRATOS_REGISTER_VARIABLE(AFC_NEGATIVE_ANTI_DIFFUSIVE_FLUX_LIMIT)

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

    // registering incompressible potential flow elements
    KRATOS_REGISTER_ELEMENT("RansIncompressiblePotentialFlowVelocity2D3N", mIncompressiblePotentialFlowVelocity2D);
    KRATOS_REGISTER_ELEMENT("RansIncompressiblePotentialFlowVelocity3D4N", mIncompressiblePotentialFlowVelocity3D);
    KRATOS_REGISTER_ELEMENT("RansIncompressiblePotentialFlowPressure2D3N", mIncompressiblePotentialFlowPressure2D);
    KRATOS_REGISTER_ELEMENT("RansIncompressiblePotentialFlowPressure3D4N", mIncompressiblePotentialFlowPressure3D);

    // registering incompressible potential flow conditions
    KRATOS_REGISTER_CONDITION("RansIncompressiblePotentialFlowVelocity2D2N", mIncompressiblePotentialFlowVelocityCondition2D2N);
    KRATOS_REGISTER_CONDITION("RansIncompressiblePotentialFlowVelocity3D3N", mIncompressiblePotentialFlowVelocityCondition3D3N);
    KRATOS_REGISTER_CONDITION("RansIncompressiblePotentialFlowPressure2D2N", mIncompressiblePotentialFlowPressureCondition2D2N);
    KRATOS_REGISTER_CONDITION("RansIncompressiblePotentialFlowPressure3D3N", mIncompressiblePotentialFlowPressureCondition3D3N);

    // registering fractional step elements
    KRATOS_REGISTER_ELEMENT("RansFractionalStep2D3N", mRansFractionalStep2D);
    KRATOS_REGISTER_ELEMENT("RansFractionalStep3D4N", mRansFractionalStep3D);

    // registering fractional step conditions
    KRATOS_REGISTER_CONDITION("RansFSHighReKWall2D2N", mFSHighReKWallCondition2D2N);
    KRATOS_REGISTER_CONDITION("RansFSHighReKWall3D3N", mFSHighReKWallCondition3D3N);

    // registering monolithic conditions
    KRATOS_REGISTER_CONDITION("RansEvmKEpsilonVmsMonolithicWall2D2N", mRansEvmKEpsilonVmsMonolithicWall2D2N);
    KRATOS_REGISTER_CONDITION("RansEvmKEpsilonVmsMonolithicWall3D3N", mRansEvmKEpsilonVmsMonolithicWall3D3N);

    KRATOS_REGISTER_CONDITION("RansVMSMonolithicKBasedWallCondition2D2N", mRansVMSMonolithicKBasedWallCondition2D2N);
    KRATOS_REGISTER_CONDITION("RansVMSMonolithicKBasedWallCondition3D3N", mRansVMSMonolithicKBasedWallCondition3D3N);

    // registering low re k-epsilon elements
    KRATOS_REGISTER_ELEMENT("RansEvmKEpsilonLowReK2D3N", mRansEvmKEpsilonLowReK2D);
    KRATOS_REGISTER_ELEMENT("RansEvmKEpsilonLowReK3D4N", mRansEvmKEpsilonLowReK3D);
    KRATOS_REGISTER_ELEMENT("RansEvmKEpsilonLowReEpsilon2D3N", mRansEvmKEpsilonLowReEpsilon2D);
    KRATOS_REGISTER_ELEMENT("RansEvmKEpsilonLowReEpsilon3D4N", mRansEvmKEpsilonLowReEpsilon3D);

    // registering k-epsilon elements - old way
    KRATOS_REGISTER_ELEMENT("RansEvmKEpsilonK2D3N", mRansEvmKEpsilonK2D);
    KRATOS_REGISTER_ELEMENT("RansEvmKEpsilonK3D4N", mRansEvmKEpsilonK3D);
    KRATOS_REGISTER_ELEMENT("RansEvmKEpsilonEpsilon2D3N", mRansEvmKEpsilonEpsilon2D);
    KRATOS_REGISTER_ELEMENT("RansEvmKEpsilonEpsilon3D4N", mRansEvmKEpsilonEpsilon3D);

    // registering k-epsilon algebraic flux correction based elements
    KRATOS_REGISTER_ELEMENT("RansEvmKEpsilonKAFC2D3N", mRansEvmKEpsilonKAFC2D);
    KRATOS_REGISTER_ELEMENT("RansEvmKEpsilonKAFC3D4N", mRansEvmKEpsilonKAFC3D);
    KRATOS_REGISTER_ELEMENT("RansEvmKEpsilonEpsilonAFC2D3N", mRansEvmKEpsilonEpsilonAFC2D);
    KRATOS_REGISTER_ELEMENT("RansEvmKEpsilonEpsilonAFC3D4N", mRansEvmKEpsilonEpsilonAFC3D);

    // registering k-epsilon residual fc based elements
    KRATOS_REGISTER_ELEMENT("RansEvmKEpsilonKResidualBasedFC2D3N", mRansEvmKEpsilonKResidualBasedFC2D);
    KRATOS_REGISTER_ELEMENT("RansEvmKEpsilonKResidualBasedFC3D4N", mRansEvmKEpsilonKResidualBasedFC3D);
    KRATOS_REGISTER_ELEMENT("RansEvmKEpsilonEpsilonResidualBasedFC2D3N", mRansEvmKEpsilonEpsilonResidualBasedFC2D);
    KRATOS_REGISTER_ELEMENT("RansEvmKEpsilonEpsilonResidualBasedFC3D4N", mRansEvmKEpsilonEpsilonResidualBasedFC3D);

    // registering k-epsilon cross wind stabilized elements
    KRATOS_REGISTER_ELEMENT("RansEvmKEpsilonKCrossWindStabilized2D3N", mRansEvmKEpsilonKCrossWindStabilized2D);
    KRATOS_REGISTER_ELEMENT("RansEvmKEpsilonKCrossWindStabilized3D4N", mRansEvmKEpsilonKCrossWindStabilized3D);
    KRATOS_REGISTER_ELEMENT("RansEvmKEpsilonEpsilonCrossWindStabilized2D3N", mRansEvmKEpsilonEpsilonCrossWindStabilized2D);
    KRATOS_REGISTER_ELEMENT("RansEvmKEpsilonEpsilonCrossWindStabilized3D4N", mRansEvmKEpsilonEpsilonCrossWindStabilized3D);

    // registering k-epsilon conditions
    KRATOS_REGISTER_CONDITION("RansEvmKEpsilonEpsilonKBasedWallCondition2D2N", mRansEvmKEpsilonEpsilonKBasedWallCondition2D2N);
    KRATOS_REGISTER_CONDITION("RansEvmKEpsilonEpsilonKBasedWallCondition3D3N", mRansEvmKEpsilonEpsilonKBasedWallCondition3D3N);

    KRATOS_REGISTER_CONDITION("RansEvmKEpsilonEpsilonUBasedWallCondition2D2N", mRansEvmKEpsilonEpsilonUBasedWallCondition2D2N);
    KRATOS_REGISTER_CONDITION("RansEvmKEpsilonEpsilonUBasedWallCondition3D3N", mRansEvmKEpsilonEpsilonUBasedWallCondition3D3N);

    KRATOS_REGISTER_CONDITION("RansEvmKEpsilonEpsilonVelocityBasedWallCondition2D2N", mRansEvmKEpsilonEpsilonVelocityBasedWallCondition2D2N);
    KRATOS_REGISTER_CONDITION("RansEvmKEpsilonEpsilonVelocityBasedWallCondition3D3N", mRansEvmKEpsilonEpsilonVelocityBasedWallCondition3D3N);

    KRATOS_REGISTER_CONDITION("RansEvmKEpsilonEpsilonKBasedLHSWallCondition2D2N", mRansEvmKEpsilonEpsilonKBasedLHSWallCondition2D2N);
    KRATOS_REGISTER_CONDITION("RansEvmKEpsilonEpsilonKBasedLHSWallCondition3D3N", mRansEvmKEpsilonEpsilonKBasedLHSWallCondition3D3N);

    KRATOS_REGISTER_CONDITION("RansEvmKEpsilonEpsilonKBasedRHSWallCondition2D2N", mRansEvmKEpsilonEpsilonKBasedRHSWallCondition2D2N);
    KRATOS_REGISTER_CONDITION("RansEvmKEpsilonEpsilonKBasedRHSWallCondition3D3N", mRansEvmKEpsilonEpsilonKBasedRHSWallCondition3D3N);

    KRATOS_REGISTER_CONDITION("RansEvmKEpsilonEpsilonWall2D2N", mRansEvmKEpsilonEpsilonWall2D2N);
    KRATOS_REGISTER_CONDITION("RansEvmKEpsilonEpsilonWall3D3N", mRansEvmKEpsilonEpsilonWall3D3N);

    // registering k-omega algebraic flux correction based elements
    KRATOS_REGISTER_ELEMENT("RansEvmKOmegaKAFC2D3N", mRansEvmKOmegaKAFC2D);
    KRATOS_REGISTER_ELEMENT("RansEvmKOmegaKAFC3D4N", mRansEvmKOmegaKAFC3D);
    KRATOS_REGISTER_ELEMENT("RansEvmKOmegaOmegaAFC2D3N", mRansEvmKOmegaOmegaAFC2D);
    KRATOS_REGISTER_ELEMENT("RansEvmKOmegaOmegaAFC3D4N", mRansEvmKOmegaOmegaAFC3D);

    // registering k-omega residual fc based elements
    KRATOS_REGISTER_ELEMENT("RansEvmKOmegaKResidualBasedFC2D3N", mRansEvmKOmegaKResidualBasedFC2D);
    KRATOS_REGISTER_ELEMENT("RansEvmKOmegaKResidualBasedFC3D4N", mRansEvmKOmegaKResidualBasedFC3D);
    KRATOS_REGISTER_ELEMENT("RansEvmKOmegaOmegaResidualBasedFC2D3N", mRansEvmKOmegaOmegaResidualBasedFC2D);
    KRATOS_REGISTER_ELEMENT("RansEvmKOmegaOmegaResidualBasedFC3D4N", mRansEvmKOmegaOmegaResidualBasedFC3D);

    // registering k-omega cross wind stabilized elements
    KRATOS_REGISTER_ELEMENT("RansEvmKOmegaKCrossWindStabilized2D3N", mRansEvmKOmegaKCrossWindStabilized2D);
    KRATOS_REGISTER_ELEMENT("RansEvmKOmegaKCrossWindStabilized3D4N", mRansEvmKOmegaKCrossWindStabilized3D);
    KRATOS_REGISTER_ELEMENT("RansEvmKOmegaOmegaCrossWindStabilized2D3N", mRansEvmKOmegaOmegaCrossWindStabilized2D);
    KRATOS_REGISTER_ELEMENT("RansEvmKOmegaOmegaCrossWindStabilized3D4N", mRansEvmKOmegaOmegaCrossWindStabilized3D);

    // registering k-omega conditions
    KRATOS_REGISTER_CONDITION("RansEvmKOmegaOmegaKBasedWallCondition2D2N", mRansEvmKOmegaOmegaKBasedWallCondition2D2N);
    KRATOS_REGISTER_CONDITION("RansEvmKOmegaOmegaKBasedWallCondition3D3N", mRansEvmKOmegaOmegaKBasedWallCondition3D3N);

    KRATOS_REGISTER_CONDITION("RansEvmKOmegaOmegaUBasedWallCondition2D2N", mRansEvmKOmegaOmegaUBasedWallCondition2D2N);
    KRATOS_REGISTER_CONDITION("RansEvmKOmegaOmegaUBasedWallCondition3D3N", mRansEvmKOmegaOmegaUBasedWallCondition3D3N);

    // registering k-omega-sst algebraic flux correction based elements
    KRATOS_REGISTER_ELEMENT("RansEvmKOmegaSSTKAFC2D3N", mRansEvmKOmegaSSTKAFC2D);
    KRATOS_REGISTER_ELEMENT("RansEvmKOmegaSSTKAFC3D4N", mRansEvmKOmegaSSTKAFC3D);
    KRATOS_REGISTER_ELEMENT("RansEvmKOmegaSSTOmegaAFC2D3N", mRansEvmKOmegaSSTOmegaAFC2D);
    KRATOS_REGISTER_ELEMENT("RansEvmKOmegaSSTOmegaAFC3D4N", mRansEvmKOmegaSSTOmegaAFC3D);

    // registering k-omega-sst residual fc based elements
    KRATOS_REGISTER_ELEMENT("RansEvmKOmegaSSTKResidualBasedFC2D3N", mRansEvmKOmegaSSTKResidualBasedFC2D);
    KRATOS_REGISTER_ELEMENT("RansEvmKOmegaSSTKResidualBasedFC3D4N", mRansEvmKOmegaSSTKResidualBasedFC3D);
    KRATOS_REGISTER_ELEMENT("RansEvmKOmegaSSTOmegaResidualBasedFC2D3N", mRansEvmKOmegaSSTOmegaResidualBasedFC2D);
    KRATOS_REGISTER_ELEMENT("RansEvmKOmegaSSTOmegaResidualBasedFC3D4N", mRansEvmKOmegaSSTOmegaResidualBasedFC3D);

    // registering k-omega-sst cross wind stabilized elements
    KRATOS_REGISTER_ELEMENT("RansEvmKOmegaSSTKCrossWindStabilized2D3N", mRansEvmKOmegaSSTKCrossWindStabilized2D);
    KRATOS_REGISTER_ELEMENT("RansEvmKOmegaSSTKCrossWindStabilized3D4N", mRansEvmKOmegaSSTKCrossWindStabilized3D);
    KRATOS_REGISTER_ELEMENT("RansEvmKOmegaSSTOmegaCrossWindStabilized2D3N", mRansEvmKOmegaSSTOmegaCrossWindStabilized2D);
    KRATOS_REGISTER_ELEMENT("RansEvmKOmegaSSTOmegaCrossWindStabilized3D4N", mRansEvmKOmegaSSTOmegaCrossWindStabilized3D);

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
