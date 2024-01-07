//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
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
      // stabilization validation elements
      mRansCircularConvectionAFC2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mRansCircularConvectionCWD2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mRansCircularConvectionRFC2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mRansBodyForceGovernedCDRAFC2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mRansBodyForceGovernedCDRCWD2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mRansBodyForceGovernedCDRRFC2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mRansDiffusionAFC2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mRansDiffusionCWD2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mRansDiffusionRFC2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      // incompressible potential flow elements
      mIncompressiblePotentialFlowVelocity2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mIncompressiblePotentialFlowVelocity3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(Element::GeometryType::PointsArrayType(4)))),
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
      mRansKEpsilonEpsilonCWD3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(Element::GeometryType::PointsArrayType(4)))),
      // k-omega turbulence model elements
      mRansKOmegaKAFC2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mRansKOmegaKAFC3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(Element::GeometryType::PointsArrayType(4)))),
      mRansKOmegaOmegaAFC2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mRansKOmegaOmegaAFC3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(Element::GeometryType::PointsArrayType(4)))),
      mRansKOmegaKRFC2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mRansKOmegaKRFC3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(Element::GeometryType::PointsArrayType(4)))),
      mRansKOmegaOmegaRFC2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mRansKOmegaOmegaRFC3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(Element::GeometryType::PointsArrayType(4)))),
      mRansKOmegaKCWD2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mRansKOmegaKCWD3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(Element::GeometryType::PointsArrayType(4)))),
      mRansKOmegaOmegaCWD2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mRansKOmegaOmegaCWD3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(Element::GeometryType::PointsArrayType(4)))),
      // k-omega-sst turbulence model elements
      mRansKOmegaSSTKAFC2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mRansKOmegaSSTKAFC3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(Element::GeometryType::PointsArrayType(4)))),
      mRansKOmegaSSTOmegaAFC2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mRansKOmegaSSTOmegaAFC3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(Element::GeometryType::PointsArrayType(4)))),
      mRansKOmegaSSTKRFC2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mRansKOmegaSSTKRFC3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(Element::GeometryType::PointsArrayType(4)))),
      mRansKOmegaSSTOmegaRFC2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mRansKOmegaSSTOmegaRFC3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(Element::GeometryType::PointsArrayType(4)))),
      mRansKOmegaSSTKCWD2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mRansKOmegaSSTKCWD3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(Element::GeometryType::PointsArrayType(4)))),
      mRansKOmegaSSTOmegaCWD2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mRansKOmegaSSTOmegaCWD3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(Element::GeometryType::PointsArrayType(4)))),
      // vms monolithic k based wall conditions
      mRansVMSMonolithicKBasedWall2D2N(0,Condition::GeometryType::Pointer(new Line2D2<Node<3>>(Condition::GeometryType::PointsArrayType(2)))),
      mRansVMSMonolithicKBasedWall3D3N(0,Condition::GeometryType::Pointer(new Triangle3D3<Node<3>>(Condition::GeometryType::PointsArrayType(3)))),
      mRansVMSMonolithicUBasedWall2D2N(0,Condition::GeometryType::Pointer(new Line2D2<Node<3>>(Condition::GeometryType::PointsArrayType(2)))),
      mRansVMSMonolithicUBasedWall3D3N(0,Condition::GeometryType::Pointer(new Triangle3D3<Node<3>>(Condition::GeometryType::PointsArrayType(3)))),
      // fractional step wall conditions
      mFractionalStepKBasedWall2D2N(0,Condition::GeometryType::Pointer(new Line2D2<Node<3>>(Condition::GeometryType::PointsArrayType(2)))),
      mFractionalStepKBasedWall3D3N(0,Condition::GeometryType::Pointer(new Triangle3D3<Node<3>>(Condition::GeometryType::PointsArrayType(3)))),
      // incompressible potential flow conditions
      mIncompressiblePotentialFlowVelocityInlet2D2N(0,Condition::GeometryType::Pointer(new Line2D2<Node<3>>(Condition::GeometryType::PointsArrayType(2)))),
      mIncompressiblePotentialFlowVelocityInlet3D3N(0,Condition::GeometryType::Pointer(new Triangle3D3<Node<3>>(Condition::GeometryType::PointsArrayType(3)))),
      // k-epsilon turbulence model conditions
      mRansKEpsilonKVisBasedWall2D2N(0,Condition::GeometryType::Pointer(new Line2D2<Node<3>>(Condition::GeometryType::PointsArrayType(2)))),
      mRansKEpsilonKVisBasedWall3D3N(0,Condition::GeometryType::Pointer(new Triangle3D3<Node<3>>(Condition::GeometryType::PointsArrayType(3)))),
      mRansKEpsilonEpsilonKBasedWall2D2N(0,Condition::GeometryType::Pointer(new Line2D2<Node<3>>(Condition::GeometryType::PointsArrayType(2)))),
      mRansKEpsilonEpsilonKBasedWall3D3N(0,Condition::GeometryType::Pointer(new Triangle3D3<Node<3>>(Condition::GeometryType::PointsArrayType(3)))),
      mRansKEpsilonEpsilonUBasedWall2D2N(0,Condition::GeometryType::Pointer(new Line2D2<Node<3>>(Condition::GeometryType::PointsArrayType(2)))),
      mRansKEpsilonEpsilonUBasedWall3D3N(0,Condition::GeometryType::Pointer(new Triangle3D3<Node<3>>(Condition::GeometryType::PointsArrayType(3)))),
      // k-omega turbulence model conditions
      mRansKOmegaOmegaKBasedWall2D2N(0,Condition::GeometryType::Pointer(new Line2D2<Node<3>>(Condition::GeometryType::PointsArrayType(2)))),
      mRansKOmegaOmegaKBasedWall3D3N(0,Condition::GeometryType::Pointer(new Triangle3D3<Node<3>>(Condition::GeometryType::PointsArrayType(3)))),
      mRansKOmegaOmegaUBasedWall2D2N(0,Condition::GeometryType::Pointer(new Line2D2<Node<3>>(Condition::GeometryType::PointsArrayType(2)))),
      mRansKOmegaOmegaUBasedWall3D3N(0,Condition::GeometryType::Pointer(new Triangle3D3<Node<3>>(Condition::GeometryType::PointsArrayType(3)))),
      mRansKOmegaOmegaVisLogBasedWall2D2N(0,Condition::GeometryType::Pointer(new Line2D2<Node<3>>(Condition::GeometryType::PointsArrayType(2)))),
      mRansKOmegaOmegaVisLogBasedWall3D3N(0,Condition::GeometryType::Pointer(new Triangle3D3<Node<3>>(Condition::GeometryType::PointsArrayType(3)))),
      // k-omega-sst turbulence model conditions
      mRansKOmegaSSTOmegaKBasedWall2D2N(0,Condition::GeometryType::Pointer(new Line2D2<Node<3>>(Condition::GeometryType::PointsArrayType(2)))),
      mRansKOmegaSSTOmegaKBasedWall3D3N(0,Condition::GeometryType::Pointer(new Triangle3D3<Node<3>>(Condition::GeometryType::PointsArrayType(3)))),
      mRansKOmegaSSTOmegaUBasedWall2D2N(0,Condition::GeometryType::Pointer(new Line2D2<Node<3>>(Condition::GeometryType::PointsArrayType(2)))),
      mRansKOmegaSSTOmegaUBasedWall3D3N(0,Condition::GeometryType::Pointer(new Triangle3D3<Node<3>>(Condition::GeometryType::PointsArrayType(3)))),
      // stabilization validation adjoint elements
      mRansCircularConvectionRFCAdjoint2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
      mRansDiffusionRFCAdjoint2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
      // stabilization adjoint validation conditions
      mRansScalarEquationAdjoint2D2N(0,Condition::GeometryType::Pointer(new Line2D2<Node<3>>(Condition::GeometryType::PointsArrayType(2)))),
      // k-epsilon adjoint elements
      mRansKEpsilonQSVMSRFCAdjoint2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
      mRansKEpsilonQSVMSRFCAdjoint3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
      // k-omega adjoint elements
      mRansKOmegaQSVMSRFCAdjoint2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
      mRansKOmegaQSVMSRFCAdjoint3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
      // k-omega-sst adjoint elements
      mRansKOmegaSSTQSVMSRFCAdjoint2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
      mRansKOmegaSSTQSVMSRFCAdjoint3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
      // k-epsilon adjoint conditions
      mRansKEpsilonVMSKBasedEpsilonKBasedWallAdjoint2D2N(0,Condition::GeometryType::Pointer(new Line2D2<Node<3>>(Condition::GeometryType::PointsArrayType(2)))),
      mRansKEpsilonVMSKBasedEpsilonKBasedWallAdjoint3D3N(0,Condition::GeometryType::Pointer(new Triangle3D3<Node<3>>(Condition::GeometryType::PointsArrayType(3)))),
      mRansKEpsilonVMSKBasedEpsilonUBasedWallAdjoint2D2N(0,Condition::GeometryType::Pointer(new Line2D2<Node<3>>(Condition::GeometryType::PointsArrayType(2)))),
      mRansKEpsilonVMSKBasedEpsilonUBasedWallAdjoint3D3N(0,Condition::GeometryType::Pointer(new Triangle3D3<Node<3>>(Condition::GeometryType::PointsArrayType(3)))),
      // k-omega adjoint conditions
      mRansKOmegaVMSKBasedOmegaKBasedWallAdjoint2D2N(0,Condition::GeometryType::Pointer(new Line2D2<Node<3>>(Condition::GeometryType::PointsArrayType(2)))),
      mRansKOmegaVMSKBasedOmegaKBasedWallAdjoint3D3N(0,Condition::GeometryType::Pointer(new Triangle3D3<Node<3>>(Condition::GeometryType::PointsArrayType(3)))),
      mRansKOmegaVMSKBasedOmegaUBasedWallAdjoint2D2N(0,Condition::GeometryType::Pointer(new Line2D2<Node<3>>(Condition::GeometryType::PointsArrayType(2)))),
      mRansKOmegaVMSKBasedOmegaUBasedWallAdjoint3D3N(0,Condition::GeometryType::Pointer(new Triangle3D3<Node<3>>(Condition::GeometryType::PointsArrayType(3)))),
      // k-omega-sst adjoint conditions
      mRansKOmegaSSTVMSKBasedOmegaKBasedWallAdjoint2D2N(0,Condition::GeometryType::Pointer(new Line2D2<Node<3>>(Condition::GeometryType::PointsArrayType(2)))),
      mRansKOmegaSSTVMSKBasedOmegaKBasedWallAdjoint3D3N(0,Condition::GeometryType::Pointer(new Triangle3D3<Node<3>>(Condition::GeometryType::PointsArrayType(3)))),
      mRansKOmegaSSTVMSKBasedOmegaUBasedWallAdjoint2D2N(0,Condition::GeometryType::Pointer(new Line2D2<Node<3>>(Condition::GeometryType::PointsArrayType(2)))),
      mRansKOmegaSSTVMSKBasedOmegaUBasedWallAdjoint3D3N(0,Condition::GeometryType::Pointer(new Triangle3D3<Node<3>>(Condition::GeometryType::PointsArrayType(3))))
{
}

void KratosRANSApplication::Register()
{
    KRATOS_INFO("") << "Initializing KratosRANSApplication..." << std::endl;

    // stabilization validation specific variables
    KRATOS_REGISTER_VARIABLE( CIRCULAR_CONVECTION_ROTATION_CLOCKWISE )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( CIRCULAR_CONVECTION_ROTATION_CENTER )
    KRATOS_REGISTER_VARIABLE( REACTION_COEFFICIENT )

    // incompressible potential flow specific variables
    KRATOS_REGISTER_VARIABLE( VELOCITY_POTENTIAL )
    KRATOS_REGISTER_VARIABLE( VELOCITY_POTENTIAL_RATE )
    KRATOS_REGISTER_VARIABLE( PRESSURE_POTENTIAL )
    KRATOS_REGISTER_VARIABLE( RANS_IS_INLET )
    KRATOS_REGISTER_VARIABLE( RANS_IS_OUTLET )
    KRATOS_REGISTER_VARIABLE( RANS_IS_STRUCTURE )
    KRATOS_REGISTER_VARIABLE( RANS_IS_STEADY )

    // residual based flux corrected stabilization variables
    KRATOS_REGISTER_VARIABLE( RANS_RESIDUAL_BASED_FLUX_CORRECTED_ACCELERATON_FACTOR )
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

    // wall function condition specific additional variables
    KRATOS_REGISTER_VARIABLE( RANS_Y_PLUS )
    KRATOS_REGISTER_VARIABLE( GAUSS_RANS_Y_PLUS )
    KRATOS_REGISTER_VARIABLE( RANS_LINEAR_LOG_LAW_Y_PLUS_LIMIT )
    KRATOS_REGISTER_VARIABLE( WALL_SMOOTHNESS_BETA )
    KRATOS_REGISTER_VARIABLE( WALL_ROUGHNESS_HEIGHT )
    KRATOS_REGISTER_VARIABLE( WALL_CORRECTION_FACTOR )
    KRATOS_REGISTER_VARIABLE( RANS_IS_WALL_FUNCTION_ACTIVE )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( FRICTION_VELOCITY )

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
    KRATOS_REGISTER_VARIABLE( TURBULENCE_RANS_A1 )
    KRATOS_REGISTER_VARIABLE( TURBULENCE_RANS_BETA_1 )
    KRATOS_REGISTER_VARIABLE( TURBULENCE_RANS_BETA_2 )
    KRATOS_REGISTER_VARIABLE( VON_KARMAN )

    // formulation specific variables
    KRATOS_REGISTER_VARIABLE( ANALYSIS_STEPS )
    KRATOS_REGISTER_VARIABLE( WALL_MODEL_PART_NAME )
    KRATOS_REGISTER_VARIABLE( NUMBER_OF_NEIGHBOUR_CONDITIONS )
    KRATOS_REGISTER_VARIABLE( NUMBER_OF_SMOOTHING_NEIGHBOUR_NODES )

    // adjoint variables
    KRATOS_REGISTER_VARIABLE( RANS_SCALAR_1_ADJOINT_1 )
    KRATOS_REGISTER_VARIABLE( RANS_SCALAR_1_ADJOINT_2 )
    KRATOS_REGISTER_VARIABLE( RANS_SCALAR_1_ADJOINT_3 )
    KRATOS_REGISTER_VARIABLE( RANS_AUX_ADJOINT_SCALAR_1 )

    KRATOS_REGISTER_VARIABLE( RANS_SCALAR_2_ADJOINT_1 )
    KRATOS_REGISTER_VARIABLE( RANS_SCALAR_2_ADJOINT_2 )
    KRATOS_REGISTER_VARIABLE( RANS_SCALAR_2_ADJOINT_3 )
    KRATOS_REGISTER_VARIABLE( RANS_AUX_ADJOINT_SCALAR_2 )

    // primal solution location storage variables
    KRATOS_REGISTER_VARIABLE(RANS_PRIMAL_SOLUTION_LOCATION_1)

    // Response function interpolation error variables for transient cases
    KRATOS_REGISTER_VARIABLE( RANS_RESPONSE_FUNCTION_DOFS_INTERPOLATION_ERROR )
    KRATOS_REGISTER_VARIABLE( RANS_RESPONSE_FUNCTION_DOFS_INTERPOLATION_ERROR_RATE )
    KRATOS_REGISTER_VARIABLE( RESPONSE_FUNCTION_INTERPOLATION_ERROR_AUXILIARY )

    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(TIME_AVERAGED_VELOCITY)
    KRATOS_REGISTER_VARIABLE(TIME_AVERAGED_VELOCITY_POTENTIAL)
    KRATOS_REGISTER_VARIABLE(TIME_AVERAGED_PRESSURE)
    KRATOS_REGISTER_VARIABLE(TIME_AVERAGED_TURBULENT_KINETIC_ENERGY)
    KRATOS_REGISTER_VARIABLE(TIME_AVERAGED_TURBULENT_ENERGY_DISSIPATION_RATE)
    KRATOS_REGISTER_VARIABLE(TIME_AVERAGED_TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE)

    // registering elements
    // registering stabilization validation elements
    KRATOS_REGISTER_ELEMENT("RansCircularConvectionAFC2D3N", mRansCircularConvectionAFC2D);
    KRATOS_REGISTER_ELEMENT("RansCircularConvectionCWD2D3N", mRansCircularConvectionCWD2D);
    KRATOS_REGISTER_ELEMENT("RansCircularConvectionRFC2D3N", mRansCircularConvectionRFC2D);
    KRATOS_REGISTER_ELEMENT("RansBodyForceGovernedCDRAFC2D3N", mRansBodyForceGovernedCDRAFC2D);
    KRATOS_REGISTER_ELEMENT("RansBodyForceGovernedCDRCWD2D3N", mRansBodyForceGovernedCDRCWD2D);
    KRATOS_REGISTER_ELEMENT("RansBodyForceGovernedCDRRFC2D3N", mRansBodyForceGovernedCDRRFC2D);
    KRATOS_REGISTER_ELEMENT("RansDiffusionAFC2D3N", mRansDiffusionAFC2D);
    KRATOS_REGISTER_ELEMENT("RansDiffusionCWD2D3N", mRansDiffusionCWD2D);
    KRATOS_REGISTER_ELEMENT("RansDiffusionRFC2D3N", mRansDiffusionRFC2D);

    // registering incompressible potential flow elements
    KRATOS_REGISTER_ELEMENT("RansIncompressiblePotentialFlowVelocity2D3N", mIncompressiblePotentialFlowVelocity2D);
    KRATOS_REGISTER_ELEMENT("RansIncompressiblePotentialFlowVelocity3D4N", mIncompressiblePotentialFlowVelocity3D);

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

    // registering k-omega algebraic flux correction based elements
    KRATOS_REGISTER_ELEMENT("RansKOmegaKAFC2D3N", mRansKOmegaKAFC2D);
    KRATOS_REGISTER_ELEMENT("RansKOmegaKAFC3D4N", mRansKOmegaKAFC3D);
    KRATOS_REGISTER_ELEMENT("RansKOmegaOmegaAFC2D3N", mRansKOmegaOmegaAFC2D);
    KRATOS_REGISTER_ELEMENT("RansKOmegaOmegaAFC3D4N", mRansKOmegaOmegaAFC3D);

    // registering k-omega residual fc based elements
    KRATOS_REGISTER_ELEMENT("RansKOmegaKRFC2D3N", mRansKOmegaKRFC2D);
    KRATOS_REGISTER_ELEMENT("RansKOmegaKRFC3D4N", mRansKOmegaKRFC3D);
    KRATOS_REGISTER_ELEMENT("RansKOmegaOmegaRFC2D3N", mRansKOmegaOmegaRFC2D);
    KRATOS_REGISTER_ELEMENT("RansKOmegaOmegaRFC3D4N", mRansKOmegaOmegaRFC3D);

    // registering k-omega cross wind stabilized elements
    KRATOS_REGISTER_ELEMENT("RansKOmegaKCWD2D3N", mRansKOmegaKCWD2D);
    KRATOS_REGISTER_ELEMENT("RansKOmegaKCWD3D4N", mRansKOmegaKCWD3D);
    KRATOS_REGISTER_ELEMENT("RansKOmegaOmegaCWD2D3N", mRansKOmegaOmegaCWD2D);
    KRATOS_REGISTER_ELEMENT("RansKOmegaOmegaCWD3D4N", mRansKOmegaOmegaCWD3D);

    // registering k-omega-sst algebraic flux correction based elements
    KRATOS_REGISTER_ELEMENT("RansKOmegaSSTKAFC2D3N", mRansKOmegaSSTKAFC2D);
    KRATOS_REGISTER_ELEMENT("RansKOmegaSSTKAFC3D4N", mRansKOmegaSSTKAFC3D);
    KRATOS_REGISTER_ELEMENT("RansKOmegaSSTOmegaAFC2D3N", mRansKOmegaSSTOmegaAFC2D);
    KRATOS_REGISTER_ELEMENT("RansKOmegaSSTOmegaAFC3D4N", mRansKOmegaSSTOmegaAFC3D);

    // registering k-omega-sst residual fc based elements
    KRATOS_REGISTER_ELEMENT("RansKOmegaSSTKRFC2D3N", mRansKOmegaSSTKRFC2D);
    KRATOS_REGISTER_ELEMENT("RansKOmegaSSTKRFC3D4N", mRansKOmegaSSTKRFC3D);
    KRATOS_REGISTER_ELEMENT("RansKOmegaSSTOmegaRFC2D3N", mRansKOmegaSSTOmegaRFC2D);
    KRATOS_REGISTER_ELEMENT("RansKOmegaSSTOmegaRFC3D4N", mRansKOmegaSSTOmegaRFC3D);

    // registering k-omega-sst cross wind stabilized elements
    KRATOS_REGISTER_ELEMENT("RansKOmegaSSTKCWD2D3N", mRansKOmegaSSTKCWD2D);
    KRATOS_REGISTER_ELEMENT("RansKOmegaSSTKCWD3D4N", mRansKOmegaSSTKCWD3D);
    KRATOS_REGISTER_ELEMENT("RansKOmegaSSTOmegaCWD2D3N", mRansKOmegaSSTOmegaCWD2D);
    KRATOS_REGISTER_ELEMENT("RansKOmegaSSTOmegaCWD3D4N", mRansKOmegaSSTOmegaCWD3D);

    // registering conditions
    // registering vms monolithic conditions
    KRATOS_REGISTER_CONDITION("RansVMSMonolithicKBasedWall2D2N", mRansVMSMonolithicKBasedWall2D2N);
    KRATOS_REGISTER_CONDITION("RansVMSMonolithicKBasedWall3D3N", mRansVMSMonolithicKBasedWall3D3N);

    KRATOS_REGISTER_CONDITION("RansVMSMonolithicUBasedWall2D2N", mRansVMSMonolithicUBasedWall2D2N);
    KRATOS_REGISTER_CONDITION("RansVMSMonolithicUBasedWall3D3N", mRansVMSMonolithicUBasedWall3D3N);

    // registering fractional step wall conditions
    KRATOS_REGISTER_CONDITION("RansFractionalStepKBasedWall2D2N", mFractionalStepKBasedWall2D2N);
    KRATOS_REGISTER_CONDITION("RansFractionalStepKBasedWall3D3N", mFractionalStepKBasedWall3D3N);

    // registering incompressible potential flow conditions
    KRATOS_REGISTER_CONDITION("RansIncompressiblePotentialFlowVelocityInlet2D2N", mIncompressiblePotentialFlowVelocityInlet2D2N);
    KRATOS_REGISTER_CONDITION("RansIncompressiblePotentialFlowVelocityInlet3D3N", mIncompressiblePotentialFlowVelocityInlet3D3N);

    // registering k-epsilon conditions
    KRATOS_REGISTER_CONDITION("RansKEpsilonKVisBasedWall2D2N", mRansKEpsilonKVisBasedWall2D2N);
    KRATOS_REGISTER_CONDITION("RansKEpsilonKVisBasedWall3D3N", mRansKEpsilonKVisBasedWall3D3N);

    KRATOS_REGISTER_CONDITION("RansKEpsilonEpsilonKBasedWall2D2N", mRansKEpsilonEpsilonKBasedWall2D2N);
    KRATOS_REGISTER_CONDITION("RansKEpsilonEpsilonKBasedWall3D3N", mRansKEpsilonEpsilonKBasedWall3D3N);

    KRATOS_REGISTER_CONDITION("RansKEpsilonEpsilonUBasedWall2D2N", mRansKEpsilonEpsilonUBasedWall2D2N);
    KRATOS_REGISTER_CONDITION("RansKEpsilonEpsilonUBasedWall3D3N", mRansKEpsilonEpsilonUBasedWall3D3N);

    // registering k-omega conditions
    KRATOS_REGISTER_CONDITION("RansKOmegaOmegaKBasedWall2D2N", mRansKOmegaOmegaKBasedWall2D2N);
    KRATOS_REGISTER_CONDITION("RansKOmegaOmegaKBasedWall3D3N", mRansKOmegaOmegaKBasedWall3D3N);

    KRATOS_REGISTER_CONDITION("RansKOmegaOmegaUBasedWall2D2N", mRansKOmegaOmegaUBasedWall2D2N);
    KRATOS_REGISTER_CONDITION("RansKOmegaOmegaUBasedWall3D3N", mRansKOmegaOmegaUBasedWall3D3N);

    KRATOS_REGISTER_CONDITION("RansKOmegaOmegaVisLogBasedWall2D2N", mRansKOmegaOmegaVisLogBasedWall2D2N);
    KRATOS_REGISTER_CONDITION("RansKOmegaOmegaVisLogBasedWall3D3N", mRansKOmegaOmegaVisLogBasedWall3D3N);

    // registering k-omega-sst conditions
    KRATOS_REGISTER_CONDITION("RansKOmegaSSTOmegaKBasedWall2D2N", mRansKOmegaSSTOmegaKBasedWall2D2N);
    KRATOS_REGISTER_CONDITION("RansKOmegaSSTOmegaKBasedWall3D3N", mRansKOmegaSSTOmegaKBasedWall3D3N);

    KRATOS_REGISTER_CONDITION("RansKOmegaSSTOmegaUBasedWall2D2N", mRansKOmegaSSTOmegaUBasedWall2D2N);
    KRATOS_REGISTER_CONDITION("RansKOmegaSSTOmegaUBasedWall3D3N", mRansKOmegaSSTOmegaUBasedWall3D3N);

    // registering constitutive laws
    KRATOS_REGISTER_CONSTITUTIVE_LAW("RansKEpsilonNewtonian2DLaw", mRansKEpsilonNewtonian2DLaw);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("RansKEpsilonNewtonian3DLaw", mRansKEpsilonNewtonian3DLaw);

    KRATOS_REGISTER_CONSTITUTIVE_LAW("RansKOmegaNewtonian2DLaw", mRansKOmegaNewtonian2DLaw);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("RansKOmegaNewtonian3DLaw", mRansKOmegaNewtonian3DLaw);

    KRATOS_REGISTER_CONSTITUTIVE_LAW("RansKOmegaSSTNewtonian2DLaw", mRansKOmegaSSTNewtonian2DLaw);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("RansKOmegaSSTNewtonian3DLaw", mRansKOmegaSSTNewtonian3DLaw);

    KRATOS_REGISTER_CONSTITUTIVE_LAW("RansFrozenTurbulenceKEpsilonNewtonian2DLaw", mRansFrozenTurbulenceKEpsilonNewtonian2DLaw);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("RansFrozenTurbulenceKEpsilonNewtonian3DLaw", mRansFrozenTurbulenceKEpsilonNewtonian3DLaw);

    KRATOS_REGISTER_CONSTITUTIVE_LAW("RansFrozenTurbulenceKOmegaNewtonian2DLaw", mRansFrozenTurbulenceKOmegaNewtonian2DLaw);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("RansFrozenTurbulenceKOmegaNewtonian3DLaw", mRansFrozenTurbulenceKOmegaNewtonian3DLaw);

    KRATOS_REGISTER_CONSTITUTIVE_LAW("RansFrozenTurbulenceKOmegaSSTNewtonian2DLaw", mRansFrozenTurbulenceKOmegaSSTNewtonian2DLaw);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("RansFrozenTurbulenceKOmegaSSTNewtonian3DLaw", mRansFrozenTurbulenceKOmegaSSTNewtonian3DLaw);

    // registering stabilization validation adjoint elements
    KRATOS_REGISTER_ELEMENT("RansCircularConvectionRFCAdjoint2D3N", mRansCircularConvectionRFCAdjoint2D3N);
    KRATOS_REGISTER_ELEMENT("RansDiffusionRFCAdjoint2D3N", mRansDiffusionRFCAdjoint2D3N);

    // registering stabilization validation adjoint conditions
    KRATOS_REGISTER_CONDITION("RansScalarEquationAdjoint2D2N", mRansScalarEquationAdjoint2D2N);

    // registering k-epsilon adjoint elements
    KRATOS_REGISTER_ELEMENT("RansKEpsilonQSVMSRFCAdjoint2D3N", mRansKEpsilonQSVMSRFCAdjoint2D3N);
    KRATOS_REGISTER_ELEMENT("RansKEpsilonQSVMSRFCAdjoint3D4N", mRansKEpsilonQSVMSRFCAdjoint3D4N);

    // registering k-omega adjoint elements
    KRATOS_REGISTER_ELEMENT("RansKOmegaQSVMSRFCAdjoint2D3N", mRansKOmegaQSVMSRFCAdjoint2D3N);
    KRATOS_REGISTER_ELEMENT("RansKOmegaQSVMSRFCAdjoint3D4N", mRansKOmegaQSVMSRFCAdjoint3D4N);

    // registering k-omega-sst adjoint elements
    KRATOS_REGISTER_ELEMENT("RansKOmegaSSTQSVMSRFCAdjoint2D3N", mRansKOmegaSSTQSVMSRFCAdjoint2D3N);
    KRATOS_REGISTER_ELEMENT("RansKOmegaSSTQSVMSRFCAdjoint3D4N", mRansKOmegaSSTQSVMSRFCAdjoint3D4N);

    // registering k-epsilon adjoint conditions
    KRATOS_REGISTER_CONDITION("RansKEpsilonVMSKBasedEpsilonKBasedWallAdjoint2D2N", mRansKEpsilonVMSKBasedEpsilonKBasedWallAdjoint2D2N);
    KRATOS_REGISTER_CONDITION("RansKEpsilonVMSKBasedEpsilonKBasedWallAdjoint3D3N", mRansKEpsilonVMSKBasedEpsilonKBasedWallAdjoint3D3N);

    KRATOS_REGISTER_CONDITION("RansKEpsilonVMSKBasedEpsilonUBasedWallAdjoint2D2N", mRansKEpsilonVMSKBasedEpsilonUBasedWallAdjoint2D2N);
    KRATOS_REGISTER_CONDITION("RansKEpsilonVMSKBasedEpsilonUBasedWallAdjoint3D3N", mRansKEpsilonVMSKBasedEpsilonUBasedWallAdjoint3D3N);

    // registering k-omega adjoint conditions
    KRATOS_REGISTER_CONDITION("RansKOmegaVMSKBasedOmegaKBasedWallAdjoint2D2N", mRansKOmegaVMSKBasedOmegaKBasedWallAdjoint2D2N);
    KRATOS_REGISTER_CONDITION("RansKOmegaVMSKBasedOmegaKBasedWallAdjoint3D3N", mRansKOmegaVMSKBasedOmegaKBasedWallAdjoint3D3N);

    KRATOS_REGISTER_CONDITION("RansKOmegaVMSKBasedOmegaUBasedWallAdjoint2D2N", mRansKOmegaVMSKBasedOmegaUBasedWallAdjoint2D2N);
    KRATOS_REGISTER_CONDITION("RansKOmegaVMSKBasedOmegaUBasedWallAdjoint3D3N", mRansKOmegaVMSKBasedOmegaUBasedWallAdjoint3D3N);

    // registering k-omega-sst adjoint conditions
    KRATOS_REGISTER_CONDITION("RansKOmegaSSTVMSKBasedOmegaKBasedWallAdjoint2D2N", mRansKOmegaSSTVMSKBasedOmegaKBasedWallAdjoint2D2N);
    KRATOS_REGISTER_CONDITION("RansKOmegaSSTVMSKBasedOmegaKBasedWallAdjoint3D3N", mRansKOmegaSSTVMSKBasedOmegaKBasedWallAdjoint3D3N);

    KRATOS_REGISTER_CONDITION("RansKOmegaSSTVMSKBasedOmegaUBasedWallAdjoint2D2N", mRansKOmegaSSTVMSKBasedOmegaUBasedWallAdjoint2D2N);
    KRATOS_REGISTER_CONDITION("RansKOmegaSSTVMSKBasedOmegaUBasedWallAdjoint3D3N", mRansKOmegaSSTVMSKBasedOmegaUBasedWallAdjoint3D3N);
}
} // namespace Kratos.
