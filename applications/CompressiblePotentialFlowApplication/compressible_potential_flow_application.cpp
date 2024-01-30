//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:        BSD License
//                  Kratos default license: kratos/license.txt
//
//
//  Main authors:    Riccardo Rossi, Inigo Lopez and Marc Nunez
//

// System includes

// External includes

// Project includes
#include "compressible_potential_flow_application.h"
#include "compressible_potential_flow_application_variables.h"
#include "geometries/line_2d_2.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/tetrahedra_3d_4.h"

namespace Kratos {

KratosCompressiblePotentialFlowApplication::KratosCompressiblePotentialFlowApplication():
    KratosApplication("CompressiblePotentialFlowApplication"),
    mIncompressiblePotentialFlowElement2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<Node >(Element::GeometryType::PointsArrayType(3)))),
    mIncompressiblePotentialFlowElement3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node >(Element::GeometryType::PointsArrayType(4)))),
    mCompressiblePotentialFlowElement2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<Node >(Element::GeometryType::PointsArrayType(3)))),
    mCompressiblePotentialFlowElement3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node >(Element::GeometryType::PointsArrayType(4)))),
    mIncompressiblePerturbationPotentialFlowElement2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<Node >(Element::GeometryType::PointsArrayType(3)))),
    mIncompressiblePerturbationPotentialFlowElement3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node >(Element::GeometryType::PointsArrayType(4)))),
    mCompressiblePerturbationPotentialFlowElement2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<Node >(Element::GeometryType::PointsArrayType(3)))),
    mCompressiblePerturbationPotentialFlowElement3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node >(Element::GeometryType::PointsArrayType(4)))),
    mTransonicPerturbationPotentialFlowElement2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<Node >(Element::GeometryType::PointsArrayType(3)))),
    mTransonicPerturbationPotentialFlowElement3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node >(Element::GeometryType::PointsArrayType(4)))),
    mAdjointAnalyticalIncompressiblePotentialFlowElement2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<Node >(Element::GeometryType::PointsArrayType(3)))),
    mAdjointIncompressiblePotentialFlowElement2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<Node >(Element::GeometryType::PointsArrayType(3)))),
    mAdjointIncompressiblePerturbationPotentialFlowElement2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<Node >(Element::GeometryType::PointsArrayType(3)))),
    mAdjointIncompressiblePerturbationPotentialFlowElement3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node >(Element::GeometryType::PointsArrayType(4)))),
    mAdjointCompressiblePotentialFlowElement2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<Node >(Element::GeometryType::PointsArrayType(3)))),
    mEmbeddedIncompressiblePotentialFlowElement2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<Node >(Element::GeometryType::PointsArrayType(3)))),
    mEmbeddedIncompressiblePotentialFlowElement3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node >(Element::GeometryType::PointsArrayType(4)))),
    mEmbeddedCompressiblePotentialFlowElement2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<Node >(Element::GeometryType::PointsArrayType(3)))),
    mEmbeddedCompressiblePotentialFlowElement3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node >(Element::GeometryType::PointsArrayType(4)))),
    mEmbeddedTransonicPerturbationPotentialFlowElement2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<Node >(Element::GeometryType::PointsArrayType(3)))),
    mEmbeddedTransonicPerturbationPotentialFlowElement3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node >(Element::GeometryType::PointsArrayType(4)))),
    mAdjointEmbeddedIncompressiblePotentialFlowElement2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<Node >(Element::GeometryType::PointsArrayType(3)))),
    mAdjointEmbeddedCompressiblePotentialFlowElement2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<Node >(Element::GeometryType::PointsArrayType(3)))),
    mPotentialWallCondition2D2N(0, Element::GeometryType::Pointer(new Line2D2<Node >(Element::GeometryType::PointsArrayType(2)))),
    mPotentialWallCondition3D3N(0, Element::GeometryType::Pointer(new Triangle3D3<Node >(Element::GeometryType::PointsArrayType(3)))),
    mAdjointPotentialWallCondition2D2N(0, Element::GeometryType::Pointer(new Line2D2<Node >(Element::GeometryType::PointsArrayType(2)))),
    mAdjointPotentialWallCondition3D3N(0, Element::GeometryType::Pointer(new Triangle3D3<Node >(Element::GeometryType::PointsArrayType(3))))
  {}

void KratosCompressiblePotentialFlowApplication::Register()
{
  KRATOS_INFO("") << "Initializing KratosCompressiblePotentialFlowApplication..." << std::endl;

  // Register Variables (defined in compressible_potential_flow_application_variables.h)
  // Degrees of freedom
  KRATOS_REGISTER_VARIABLE(VELOCITY_POTENTIAL);
  KRATOS_REGISTER_VARIABLE(AUXILIARY_VELOCITY_POTENTIAL);

  // Reaction variables (Degrees of freedom)
  KRATOS_REGISTER_VARIABLE(REACTION_VELOCITY_POTENTIAL);
  KRATOS_REGISTER_VARIABLE(REACTION_AUXILIARY_VELOCITY_POTENTIAL);

  //Embedded variables
  KRATOS_REGISTER_VARIABLE(GEOMETRY_DISTANCE);
  KRATOS_REGISTER_VARIABLE(ROTATION_ANGLE);

  // Wake variables
  KRATOS_REGISTER_VARIABLE(WAKE_DISTANCE);
  KRATOS_REGISTER_VARIABLE(WAKE_ELEMENTAL_DISTANCES);
  KRATOS_REGISTER_VARIABLE(WAKE_ORIGIN);

  // Adjoint variables
  KRATOS_REGISTER_VARIABLE(ADJOINT_VELOCITY_POTENTIAL);
  KRATOS_REGISTER_VARIABLE(ADJOINT_AUXILIARY_VELOCITY_POTENTIAL);

  // Flow field magnitudes
  KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(PERTURBATION_VELOCITY);
  KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(VELOCITY_LOWER);
  KRATOS_REGISTER_VARIABLE(PRESSURE_LOWER);
  KRATOS_REGISTER_VARIABLE(POTENTIAL_JUMP);
  KRATOS_REGISTER_VARIABLE(ENERGY_NORM_REFERENCE);
  KRATOS_REGISTER_VARIABLE(POTENTIAL_ENERGY_REFERENCE);

  // Free stream magnitudes
  KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(FREE_STREAM_VELOCITY);
  KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(FREE_STREAM_VELOCITY_DIRECTION);
  KRATOS_REGISTER_VARIABLE(FREE_STREAM_DENSITY);
  KRATOS_REGISTER_VARIABLE(FREE_STREAM_MACH);

  // Integral magnitudes
  KRATOS_REGISTER_VARIABLE(LIFT_COEFFICIENT);
  KRATOS_REGISTER_VARIABLE(MOMENT_COEFFICIENT);
  KRATOS_REGISTER_VARIABLE(LIFT_COEFFICIENT_JUMP);
  KRATOS_REGISTER_VARIABLE(LIFT_COEFFICIENT_FAR_FIELD);
  KRATOS_REGISTER_VARIABLE(DRAG_COEFFICIENT_FAR_FIELD);

  // Geometrical variables
  KRATOS_REGISTER_VARIABLE(REFERENCE_CHORD)
  KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(WAKE_NORMAL);
  KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(WING_SPAN_DIRECTION);
  KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(VECTOR_TO_UPWIND_ELEMENT);

  // Solver parameters
  KRATOS_REGISTER_VARIABLE(MACH_LIMIT)
  KRATOS_REGISTER_VARIABLE(CRITICAL_MACH)
  KRATOS_REGISTER_VARIABLE(UPWIND_FACTOR_CONSTANT)

  // Solver settings
  KRATOS_REGISTER_VARIABLE(ECHO_LEVEL)

  // Markers
  KRATOS_REGISTER_VARIABLE(WAKE);
  KRATOS_REGISTER_VARIABLE(KUTTA);
  KRATOS_REGISTER_VARIABLE(WING_TIP);
  KRATOS_REGISTER_VARIABLE(TRAILING_EDGE);
  KRATOS_REGISTER_VARIABLE(UPPER_SURFACE);
  KRATOS_REGISTER_VARIABLE(LOWER_SURFACE);
  KRATOS_REGISTER_VARIABLE(UPPER_WAKE);
  KRATOS_REGISTER_VARIABLE(LOWER_WAKE);
  KRATOS_REGISTER_VARIABLE(AIRFOIL);
  KRATOS_REGISTER_VARIABLE(FAR_FIELD);

  // To be removed
  KRATOS_REGISTER_VARIABLE(TRAILING_EDGE_ELEMENT);
  KRATOS_REGISTER_VARIABLE(DECOUPLED_TRAILING_EDGE_ELEMENT);
  KRATOS_REGISTER_VARIABLE(DEACTIVATED_WAKE);
  KRATOS_REGISTER_VARIABLE(ALL_TRAILING_EDGE);
  KRATOS_REGISTER_VARIABLE(ZERO_VELOCITY_CONDITION);

  //Register elements
  KRATOS_REGISTER_ELEMENT("IncompressiblePotentialFlowElement2D3N", mIncompressiblePotentialFlowElement2D3N);
  KRATOS_REGISTER_ELEMENT("IncompressiblePotentialFlowElement3D4N", mIncompressiblePotentialFlowElement3D4N);
  KRATOS_REGISTER_ELEMENT("CompressiblePotentialFlowElement2D3N", mCompressiblePotentialFlowElement2D3N);
  KRATOS_REGISTER_ELEMENT("CompressiblePotentialFlowElement3D4N", mCompressiblePotentialFlowElement3D4N);
  KRATOS_REGISTER_ELEMENT("IncompressiblePerturbationPotentialFlowElement2D3N", mIncompressiblePerturbationPotentialFlowElement2D3N);
  KRATOS_REGISTER_ELEMENT("IncompressiblePerturbationPotentialFlowElement3D4N", mIncompressiblePerturbationPotentialFlowElement3D4N);
  KRATOS_REGISTER_ELEMENT("CompressiblePerturbationPotentialFlowElement2D3N", mCompressiblePerturbationPotentialFlowElement2D3N);
  KRATOS_REGISTER_ELEMENT("CompressiblePerturbationPotentialFlowElement3D4N", mCompressiblePerturbationPotentialFlowElement3D4N);
  KRATOS_REGISTER_ELEMENT("TransonicPerturbationPotentialFlowElement2D3N", mTransonicPerturbationPotentialFlowElement2D3N);
  KRATOS_REGISTER_ELEMENT("TransonicPerturbationPotentialFlowElement3D4N", mTransonicPerturbationPotentialFlowElement3D4N);
  KRATOS_REGISTER_ELEMENT("AdjointAnalyticalIncompressiblePotentialFlowElement2D3N", mAdjointAnalyticalIncompressiblePotentialFlowElement2D3N);
  KRATOS_REGISTER_ELEMENT("AdjointIncompressiblePotentialFlowElement2D3N", mAdjointIncompressiblePotentialFlowElement2D3N);
  KRATOS_REGISTER_ELEMENT("AdjointIncompressiblePerturbationPotentialFlowElement2D3N", mAdjointIncompressiblePerturbationPotentialFlowElement2D3N);
  KRATOS_REGISTER_ELEMENT("AdjointIncompressiblePerturbationPotentialFlowElement3D4N", mAdjointIncompressiblePerturbationPotentialFlowElement3D4N);
  KRATOS_REGISTER_ELEMENT("AdjointCompressiblePotentialFlowElement2D3N", mAdjointCompressiblePotentialFlowElement2D3N);
  KRATOS_REGISTER_ELEMENT("EmbeddedIncompressiblePotentialFlowElement2D3N", mEmbeddedIncompressiblePotentialFlowElement2D3N);
  KRATOS_REGISTER_ELEMENT("EmbeddedIncompressiblePotentialFlowElement3D4N", mEmbeddedIncompressiblePotentialFlowElement3D4N);
  KRATOS_REGISTER_ELEMENT("EmbeddedCompressiblePotentialFlowElement2D3N", mEmbeddedCompressiblePotentialFlowElement2D3N);
  KRATOS_REGISTER_ELEMENT("EmbeddedCompressiblePotentialFlowElement3D4N", mEmbeddedCompressiblePotentialFlowElement3D4N);
  KRATOS_REGISTER_ELEMENT("EmbeddedTransonicPerturbationPotentialFlowElement2D3N", mEmbeddedTransonicPerturbationPotentialFlowElement2D3N);
  KRATOS_REGISTER_ELEMENT("EmbeddedTransonicPerturbationPotentialFlowElement3D4N", mEmbeddedTransonicPerturbationPotentialFlowElement3D4N);
  KRATOS_REGISTER_ELEMENT("AdjointEmbeddedIncompressiblePotentialFlowElement2D3N", mAdjointEmbeddedIncompressiblePotentialFlowElement2D3N);
  KRATOS_REGISTER_ELEMENT("AdjointEmbeddedCompressiblePotentialFlowElement2D3N", mAdjointEmbeddedCompressiblePotentialFlowElement2D3N);

  //Register conditions
  KRATOS_REGISTER_CONDITION("PotentialWallCondition2D2N", mPotentialWallCondition2D2N);
  KRATOS_REGISTER_CONDITION("PotentialWallCondition3D3N", mPotentialWallCondition3D3N);
  KRATOS_REGISTER_CONDITION("AdjointPotentialWallCondition2D2N", mAdjointPotentialWallCondition2D2N);
  KRATOS_REGISTER_CONDITION("AdjointPotentialWallCondition3D3N", mAdjointPotentialWallCondition3D3N);
}

}  // namespace Kratos.
