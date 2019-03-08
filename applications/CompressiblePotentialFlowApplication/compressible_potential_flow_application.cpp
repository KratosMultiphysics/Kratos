//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Inigo Lopez and Riccardo Rossi
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
    mCompressiblePotentialFlowElement2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
    mCompressiblePotentialFlowElement3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    mIncompressiblePotentialFlowElement2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
    mPotentialWallCondition2D2N(0, Element::GeometryType::Pointer(new Line2D2<Node<3> >(Element::GeometryType::PointsArrayType(2)))),
    mPotentialWallCondition3D3N(0, Element::GeometryType::Pointer(new Triangle3D3<Node<3> >(Element::GeometryType::PointsArrayType(3))))
  {}

void KratosCompressiblePotentialFlowApplication::Register()
{
    // calling base class register to register Kratos components
    KratosApplication::Register();
    std::cout << "Initializing KratosCompressiblePotentialFlowApplication... " << std::endl;

    // Register Variables (defined in compressible_potential_flow_application_variables.h)
    // Degrees of freedom
    KRATOS_REGISTER_VARIABLE(VELOCITY_POTENTIAL);
    KRATOS_REGISTER_VARIABLE(AUXILIARY_VELOCITY_POTENTIAL);

    // Flow field magnitudes
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(VELOCITY_INFINITY);
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(VELOCITY_LOWER);
    KRATOS_REGISTER_VARIABLE(PRESSURE_LOWER);
    KRATOS_REGISTER_VARIABLE(POTENTIAL_JUMP);
    KRATOS_REGISTER_VARIABLE(ENERGY_NORM_REFERENCE);
    KRATOS_REGISTER_VARIABLE(POTENTIAL_ENERGY_REFERENCE);

    // Markers
    KRATOS_REGISTER_VARIABLE(WAKE);
    KRATOS_REGISTER_VARIABLE(KUTTA);
    KRATOS_REGISTER_VARIABLE(TRAILING_EDGE);
    KRATOS_REGISTER_VARIABLE(UPPER_SURFACE);
    KRATOS_REGISTER_VARIABLE(LOWER_SURFACE);
    KRATOS_REGISTER_VARIABLE(UPPER_WAKE);
    KRATOS_REGISTER_VARIABLE(LOWER_WAKE);
    KRATOS_REGISTER_VARIABLE(AIRFOIL);

    // To be removed
    KRATOS_REGISTER_VARIABLE(TRAILING_EDGE_ELEMENT);
    KRATOS_REGISTER_VARIABLE(DECOUPLED_TRAILING_EDGE_ELEMENT);
    KRATOS_REGISTER_VARIABLE(DEACTIVATED_WAKE);
    KRATOS_REGISTER_VARIABLE(ALL_TRAILING_EDGE);
    KRATOS_REGISTER_VARIABLE(ZERO_VELOCITY_CONDITION);

    //Register elements
    KRATOS_REGISTER_ELEMENT("CompressiblePotentialFlowElement2D3N", mCompressiblePotentialFlowElement2D3N); //this is the name the element should have according to the naming convention
    KRATOS_REGISTER_ELEMENT("CompressiblePotentialFlowElement3D4N", mCompressiblePotentialFlowElement3D4N);
    KRATOS_REGISTER_ELEMENT("IncompressiblePotentialFlowElement2D3N", mIncompressiblePotentialFlowElement2D3N);

    //Register conditions
    KRATOS_REGISTER_CONDITION("PotentialWallCondition2D2N", mPotentialWallCondition2D2N);
    KRATOS_REGISTER_CONDITION("PotentialWallCondition3D3N", mPotentialWallCondition3D3N);
}

}  // namespace Kratos.
