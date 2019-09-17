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
      mEVMLowReKElement2D3N(0,
                            Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(
                                Element::GeometryType::PointsArrayType(3)))),
      mEVMLowReKElement3D4N(0,
                            Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(
                                Element::GeometryType::PointsArrayType(4)))),
      mEVMLowReEpsilonElement2D3N(
          0,
          Element::GeometryType::Pointer(
              new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mEVMLowReEpsilonElement3D4N(
          0,
          Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(
              Element::GeometryType::PointsArrayType(4)))),
      mEVMKElement2D3N(0,
                       Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(
                           Element::GeometryType::PointsArrayType(3)))),
      mEVMKElement3D4N(0,
                       Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(
                           Element::GeometryType::PointsArrayType(4)))),
      mEVMEpsilonElement2D3N(0,
                             Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(
                                 Element::GeometryType::PointsArrayType(3)))),
      mEVMEpsilonElement3D4N(0,
                             Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(
                                 Element::GeometryType::PointsArrayType(4)))),
      mEVMEpsilonWallCondition2D2N(
          0,
          Element::GeometryType::Pointer(
              new Line2D2<Node<3>>(Element::GeometryType::PointsArrayType(2)))),
      mEVMEpsilonWallCondition3D3N(
          0,
          Element::GeometryType::Pointer(
              new Triangle3D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mEVMVMSMonolithicWallCondition2D2N(
          0,
          Element::GeometryType::Pointer(
              new Line2D2<Node<3>>(Element::GeometryType::PointsArrayType(2)))),
      mEVMVMSMonolithicWallCondition3D3N(
          0,
          Element::GeometryType::Pointer(
              new Triangle3D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mEVMEpsilonAdjointElement2D3N(
          0,
          Element::GeometryType::Pointer(
              new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mEVMEpsilonAdjointElement3D4N(
          0,
          Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(
              Element::GeometryType::PointsArrayType(4)))),
      mEVMKAdjointElement2D3N(0,
                              Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(
                                  Element::GeometryType::PointsArrayType(3)))),
      mEVMKAdjointElement3D4N(0,
                              Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(
                                  Element::GeometryType::PointsArrayType(4)))),
      mEVMKEpsilonVMSAdjointElement2D3N(
          0,
          Element::GeometryType::Pointer(
              new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mEVMKEpsilonVMSAdjointElement3D4N(
          0,
          Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(
              Element::GeometryType::PointsArrayType(4)))),
      mEVMMonolithicKEpsilonVMSAdjointElement2D3N(
          0,
          Element::GeometryType::Pointer(
              new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mEVMMonolithicKEpsilonVMSAdjointElement3D4N(
          0,
          Element::GeometryType::Pointer(
              new Tetrahedra3D4<Node<3>>(Element::GeometryType::PointsArrayType(4))))
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
    KRATOS_REGISTER_VARIABLE(RANS_WALL_Y_PLUS)
    KRATOS_REGISTER_VARIABLE(RANS_AUXILIARY_VARIABLE_1)
    KRATOS_REGISTER_VARIABLE(RANS_AUXILIARY_VARIABLE_2)
    KRATOS_REGISTER_VARIABLE(WALL_SMOOTHNESS_BETA)
    KRATOS_REGISTER_VARIABLE(WALL_VON_KARMAN)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(WALL_VELOCITY)
    KRATOS_REGISTER_VARIABLE(TURBULENCE_RANS_C_MU)
    KRATOS_REGISTER_VARIABLE(TURBULENCE_RANS_C1)
    KRATOS_REGISTER_VARIABLE(TURBULENCE_RANS_C2)
    KRATOS_REGISTER_VARIABLE(TURBULENT_VISCOSITY_MIN)
    KRATOS_REGISTER_VARIABLE(TURBULENT_VISCOSITY_MAX)
    KRATOS_REGISTER_VARIABLE(TURBULENT_KINETIC_ENERGY_SIGMA)
    KRATOS_REGISTER_VARIABLE(TURBULENT_ENERGY_DISSIPATION_RATE_SIGMA)

    KRATOS_REGISTER_VARIABLE(RANS_AUXILIARY_PROCESS_LIST)
    KRATOS_REGISTER_VARIABLE(PARENT_ELEMENT)
    KRATOS_REGISTER_VARIABLE(NUMBER_OF_NEIGHBOUR_CONDITIONS)

    // Register adjoint variables
    KRATOS_REGISTER_VARIABLE(RANS_NUT_PARTIAL_DERIVATIVES)
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

    // Register primal elements
    KRATOS_REGISTER_ELEMENT("EVMLowReKElement2D3N", mEVMLowReKElement2D3N);
    KRATOS_REGISTER_ELEMENT("EVMLowReKElement3D4N", mEVMLowReKElement3D4N);
    KRATOS_REGISTER_ELEMENT("EVMLowReEpsilonElement2D3N", mEVMLowReEpsilonElement2D3N);
    KRATOS_REGISTER_ELEMENT("EVMLowReEpsilonElement3D4N", mEVMLowReEpsilonElement3D4N);

    KRATOS_REGISTER_ELEMENT("EVMKElement2D3N", mEVMKElement2D3N);
    KRATOS_REGISTER_ELEMENT("EVMKElement3D4N", mEVMKElement3D4N);
    KRATOS_REGISTER_ELEMENT("EVMEpsilonElement2D3N", mEVMEpsilonElement2D3N);
    KRATOS_REGISTER_ELEMENT("EVMEpsilonElement3D4N", mEVMEpsilonElement3D4N);

    KRATOS_REGISTER_CONDITION("EVMEpsilonWallCondition2D2N", mEVMEpsilonWallCondition2D2N);
    KRATOS_REGISTER_CONDITION("EVMEpsilonWallCondition3D3N", mEVMEpsilonWallCondition3D3N);

    KRATOS_REGISTER_CONDITION("EVMVMSMonolithicWallCondition2D2N",
                              mEVMVMSMonolithicWallCondition2D2N);
    KRATOS_REGISTER_CONDITION("EVMVMSMonolithicWallCondition3D3N",
                              mEVMVMSMonolithicWallCondition3D3N);

    // Registering adjoint elements
    KRATOS_REGISTER_ELEMENT("EVMEpsilonAdjointElement2D3N", mEVMEpsilonAdjointElement2D3N);
    KRATOS_REGISTER_ELEMENT("EVMEpsilonAdjointElement3D4N", mEVMEpsilonAdjointElement3D4N);

    KRATOS_REGISTER_ELEMENT("EVMKAdjointElement2D3N", mEVMKAdjointElement2D3N);
    KRATOS_REGISTER_ELEMENT("EVMKAdjointElement3D4N", mEVMKAdjointElement3D4N);

    KRATOS_REGISTER_ELEMENT("EVMKEpsilonVMSAdjointElement2D3N", mEVMKEpsilonVMSAdjointElement2D3N);
    KRATOS_REGISTER_ELEMENT("EVMKEpsilonVMSAdjointElement3D4N", mEVMKEpsilonVMSAdjointElement3D4N);

    KRATOS_REGISTER_ELEMENT("EVMMonolithicKEpsilonVMSAdjointElement2D3N",
                            mEVMMonolithicKEpsilonVMSAdjointElement2D3N);
    KRATOS_REGISTER_ELEMENT("EVMMonolithicKEpsilonVMSAdjointElement3D4N",
                            mEVMMonolithicKEpsilonVMSAdjointElement3D4N);
}
} // namespace Kratos.
