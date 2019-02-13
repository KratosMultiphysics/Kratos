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
    mIncompressiblePotentialFlowElement2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
    mIncompressibleFullPotentialFlowElement2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
    mIncompressibleStressesPotentialFlowElement2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
    mIncompressibleAlphaPotentialFlowElement2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
    mIncompressibleAlphaFullPotentialFlowElement2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
    mIncompressibleStressesMixPotentialFlowElement2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
    // mIncompressiblePotentialFlowElement3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    mCompressiblePotentialFlowElement2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
    mCompressiblePotentialFlowElement3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    mCompressibleFullPotentialFlowElement2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
    mIncompressiblePotentialWallCondition2D2N(0, Element::GeometryType::Pointer(new Line2D2<Node<3> >(Element::GeometryType::PointsArrayType(2)))),
    mIncompressibleStressesPotentialWallCondition2D2N(0, Element::GeometryType::Pointer(new Line2D2<Node<3> >(Element::GeometryType::PointsArrayType(2)))),
    mIncompressiblePotentialWallCondition3D3N(0, Element::GeometryType::Pointer(new Triangle3D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
    mPotentialWallCondition2D2N(0, Element::GeometryType::Pointer(new Line2D2<Node<3> >(Element::GeometryType::PointsArrayType(2)))),
    mPotentialWallCondition3D3N(0, Element::GeometryType::Pointer(new Triangle3D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
    mCompressiblePotentialWallCondition2D2N(0, Element::GeometryType::Pointer(new Line2D2<Node<3> >(Element::GeometryType::PointsArrayType(2)))),
    mCompressiblePotentialWallCondition3D3N(0, Element::GeometryType::Pointer(new Triangle3D3<Node<3> >(Element::GeometryType::PointsArrayType(3))))
 
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

    // Embedded and adjoint variables
    KRATOS_REGISTER_VARIABLE( WAKE_DISTANCE );
    KRATOS_REGISTER_VARIABLE( LEVEL_SET_DISTANCE );
    
    KRATOS_REGISTER_VARIABLE( ADJOINT_VELOCITY_POTENTIAL );
    KRATOS_REGISTER_VARIABLE( ADJOINT_AUXILIARY_VELOCITY_POTENTIAL );
    KRATOS_REGISTER_VARIABLE( DISTANCE_SENSITIVITY );
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(COORDINATES_SENSITIVITY);
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(COORDINATES);

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
    KRATOS_REGISTER_ELEMENT("IncompressiblePotentialFlowElement2D3N",mIncompressiblePotentialFlowElement2D3N); //this is the name the element should have according to the naming convention
    KRATOS_REGISTER_ELEMENT("IncompressibleFullPotentialFlowElement2D3N",mIncompressibleFullPotentialFlowElement2D3N); //this is the name the element should have according to the naming convention
    KRATOS_REGISTER_ELEMENT("IncompressibleAdjointPotentialFlowElement2D3N",mIncompressibleAdjointPotentialFlowElement2D3N); //this is the name the element should have according to the naming convention
    KRATOS_REGISTER_ELEMENT("IncompressibleStressesPotentialFlowElement2D3N",mIncompressibleStressesPotentialFlowElement2D3N); //this is the name the element should have according to the naming convention
    KRATOS_REGISTER_ELEMENT("IncompressibleAlphaPotentialFlowElement2D3N",mIncompressibleAlphaPotentialFlowElement2D3N); //this is the name the element should have according to the naming convention
    KRATOS_REGISTER_ELEMENT("IncompressibleAlphaFullPotentialFlowElement2D3N",mIncompressibleAlphaFullPotentialFlowElement2D3N); //this is the name the element should have according to the naming convention
    KRATOS_REGISTER_ELEMENT("IncompressibleStressesMixPotentialFlowElement2D3N",mIncompressibleStressesMixPotentialFlowElement2D3N); //this is the name the element should have according to the naming convention
    // KRATOS_REGISTER_ELEMENT("IncompressiblePotentialFlowElement3D4N",mIncompressiblePotentialFlowElement3D4N); //this is the name the element should have according to the naming convention

    KRATOS_REGISTER_ELEMENT("CompressiblePotentialFlowElement2D3N",mCompressiblePotentialFlowElement2D3N); //this is the name the element should have according to the naming convention
    KRATOS_REGISTER_ELEMENT("CompressiblePotentialFlowElement3D4N",mCompressiblePotentialFlowElement3D4N); //this is the name the element should have according to the naming convention
    KRATOS_REGISTER_ELEMENT("CompressibleFullPotentialFlowElement2D3N",mCompressibleFullPotentialFlowElement2D3N); //this is the name the element should have according to the naming convention
    //Register conditions
    KRATOS_REGISTER_CONDITION("IncompressiblePotentialWallCondition2D2N",mIncompressiblePotentialWallCondition2D2N); //this is the name the element should have according to the naming convention
    KRATOS_REGISTER_CONDITION("IncompressibleAdjointPotentialWallCondition2D2N",mIncompressibleAdjointPotentialWallCondition2D2N); //this is the name the element should have according to the naming convention
    KRATOS_REGISTER_CONDITION("IncompressibleStressesPotentialWallCondition2D2N",mIncompressibleStressesPotentialWallCondition2D2N); //this is the name the element should have according to the naming convention
    KRATOS_REGISTER_CONDITION("IncompressiblePotentialWallCondition3D3N",mIncompressiblePotentialWallCondition3D3N); //this is the name the element should have according to the naming convention
    
    KRATOS_REGISTER_CONDITION("PotentialWallCondition2D2N", mPotentialWallCondition2D2N);
    KRATOS_REGISTER_CONDITION("PotentialWallCondition3D3N", mPotentialWallCondition3D3N);
    KRATOS_REGISTER_CONDITION("CompressiblePotentialWallCondition2D2N", mCompressiblePotentialWallCondition2D2N);
    KRATOS_REGISTER_CONDITION("CompressiblePotentialWallCondition3D3N", mCompressiblePotentialWallCondition3D3N);
}

}  // namespace Kratos.
