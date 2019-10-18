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
#include "rans_modelling_application.h"
#include "geometries/line_2d_2.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/triangle_2d_3.h"
#include "rans_modelling_application_variables.h"

namespace Kratos
{
KratosRANSModellingApplication::KratosRANSModellingApplication()
    : KratosApplication("RANSModellingApplication"),
      mRansEvmLowReK2D(0,
                       Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(
                           Element::GeometryType::PointsArrayType(3)))),
      mRansEvmLowReK3D(0,
                       Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(
                           Element::GeometryType::PointsArrayType(4)))),
      mRansEvmLowReEpsilon2D(0,
                             Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(
                                 Element::GeometryType::PointsArrayType(3)))),
      mRansEvmLowReEpsilon3D(0,
                             Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(
                                 Element::GeometryType::PointsArrayType(4)))),
      mRansEvmK2D(0,
                  Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(
                      Element::GeometryType::PointsArrayType(3)))),
      mRansEvmK3D(0,
                  Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(
                      Element::GeometryType::PointsArrayType(4)))),
      mRansEvmEpsilon2D(0,
                        Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(
                            Element::GeometryType::PointsArrayType(3)))),
      mRansEvmEpsilon3D(0,
                        Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(
                            Element::GeometryType::PointsArrayType(4)))),
      mRansEvmEpsilonWallCondition2D2N(
          0,
          Element::GeometryType::Pointer(
              new Line2D2<Node<3>>(Element::GeometryType::PointsArrayType(2)))),
      mRansEvmEpsilonWallCondition3D3N(
          0,
          Element::GeometryType::Pointer(
              new Triangle3D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mRansEvmVmsMonolithicWallCondition2D2N(
          0,
          Element::GeometryType::Pointer(
              new Line2D2<Node<3>>(Element::GeometryType::PointsArrayType(2)))),
      mRansEvmVmsMonolithicWallCondition3D3N(
          0,
          Element::GeometryType::Pointer(
              new Triangle3D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mRansEvmEpsilonAdjoint2D3N(0,
                                 Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(
                                     Element::GeometryType::PointsArrayType(3)))),
      mRansEvmEpsilonAdjoint3D4N(
          0,
          Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(
              Element::GeometryType::PointsArrayType(4)))),
      mRansEvmKAdjoint2D3N(0,
                           Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(
                               Element::GeometryType::PointsArrayType(3)))),
      mRansEvmKAdjoint3D4N(0,
                           Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(
                               Element::GeometryType::PointsArrayType(4)))),
      mRansEvmKEpsilonVMSAdjoint2D3N(
          0,
          Element::GeometryType::Pointer(
              new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mRansEvmKEpsilonVMSAdjoint3D4N(
          0,
          Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(
              Element::GeometryType::PointsArrayType(4)))),
      mRansEvmMonolithicKEpsilonVMSAdjoint2D(
          0,
          Element::GeometryType::Pointer(
              new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mRansEvmMonolithicKEpsilonVMSAdjoint3D(
          0,
          Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(
              Element::GeometryType::PointsArrayType(4)))),
      mRansEvmEpsilonAdjointWallCondition2D2N(
          0,
          Element::GeometryType::Pointer(
              new Line2D2<Node<3>>(Element::GeometryType::PointsArrayType(2)))),
      mRansEvmEpsilonAdjointWallCondition3D3N(
          0,
          Element::GeometryType::Pointer(
              new Triangle3D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mRansEvmVmsMonolithicAdjointWallCondition2D2N(
          0,
          Element::GeometryType::Pointer(
              new Line2D2<Node<3>>(Element::GeometryType::PointsArrayType(2)))),
      mRansEvmVmsMonolithicAdjointWallCondition3D3N(
          0,
          Element::GeometryType::Pointer(
              new Triangle3D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mRansEvmMonolithicKEpsilonVMSAdjointWallCondition2D2N(
          0,
          Element::GeometryType::Pointer(
              new Line2D2<Node<3>>(Element::GeometryType::PointsArrayType(2)))),
      mRansEvmMonolithicKEpsilonVMSAdjointWallCondition3D3N(
          0,
          Element::GeometryType::Pointer(
              new Triangle3D3<Node<3>>(Element::GeometryType::PointsArrayType(3))))
{
}

void KratosRANSModellingApplication::Register()
{
    // calling base class register to register Kratos components
    KratosApplication::Register();
    KRATOS_INFO("") << "Initializing KratosRANSModellingApplication..." << std::endl;

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
    KRATOS_REGISTER_VARIABLE(TURBULENT_VISCOSITY_MIN)
    KRATOS_REGISTER_VARIABLE(TURBULENT_VISCOSITY_MAX)
    KRATOS_REGISTER_VARIABLE(TURBULENT_KINETIC_ENERGY_SIGMA)
    KRATOS_REGISTER_VARIABLE(TURBULENT_ENERGY_DISSIPATION_RATE_SIGMA)

    KRATOS_REGISTER_VARIABLE(PARENT_ELEMENT)
    KRATOS_REGISTER_VARIABLE(NUMBER_OF_NEIGHBOUR_CONDITIONS)

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
    KRATOS_REGISTER_ELEMENT("RansEvmLowReK2D3N", mRansEvmLowReK2D);
    KRATOS_REGISTER_ELEMENT("RansEvmLowReK3D4N", mRansEvmLowReK3D);
    KRATOS_REGISTER_ELEMENT("RansEvmLowReEpsilon2D3N", mRansEvmLowReEpsilon2D);
    KRATOS_REGISTER_ELEMENT("RansEvmLowReEpsilon3D4N", mRansEvmLowReEpsilon3D);

    KRATOS_REGISTER_ELEMENT("RansEvmK2D3N", mRansEvmK2D);
    KRATOS_REGISTER_ELEMENT("RansEvmK3D4N", mRansEvmK3D);
    KRATOS_REGISTER_ELEMENT("RansEvmEpsilon2D3N", mRansEvmEpsilon2D);
    KRATOS_REGISTER_ELEMENT("RansEvmEpsilon3D4N", mRansEvmEpsilon3D);

    KRATOS_REGISTER_CONDITION("RansEvmEpsilonWallCondition2D2N", mRansEvmEpsilonWallCondition2D2N);
    KRATOS_REGISTER_CONDITION("RansEvmEpsilonWallCondition3D3N", mRansEvmEpsilonWallCondition3D3N);

    KRATOS_REGISTER_CONDITION("RansEvmVmsMonolithicWallCondition2D2N",
                              mRansEvmVmsMonolithicWallCondition2D2N);
    KRATOS_REGISTER_CONDITION("RansEvmVmsMonolithicWallCondition3D3N",
                              mRansEvmVmsMonolithicWallCondition3D3N);

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
