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
      mRANSEVMLowReK2D(0,
                       Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(
                           Element::GeometryType::PointsArrayType(3)))),
      mRANSEVMLowReK3D(0,
                       Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(
                           Element::GeometryType::PointsArrayType(4)))),
      mRANSEVMLowReEpsilon2D(0,
                             Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(
                                 Element::GeometryType::PointsArrayType(3)))),
      mRANSEVMLowReEpsilon3D(0,
                             Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(
                                 Element::GeometryType::PointsArrayType(4)))),
      mRANSEVMK2D(0,
                  Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(
                      Element::GeometryType::PointsArrayType(3)))),
      mRANSEVMK3D(0,
                  Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(
                      Element::GeometryType::PointsArrayType(4)))),
      mRANSEVMEpsilon2D(0,
                        Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(
                            Element::GeometryType::PointsArrayType(3)))),
      mRANSEVMEpsilon3D(0,
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

    // Register Elements
    KRATOS_REGISTER_ELEMENT("RANSEVMLowReK2D3N", mRANSEVMLowReK2D);
    KRATOS_REGISTER_ELEMENT("RANSEVMLowReK3D4N", mRANSEVMLowReK3D);
    KRATOS_REGISTER_ELEMENT("RANSEVMLowReEpsilon2D3N", mRANSEVMLowReEpsilon2D);
    KRATOS_REGISTER_ELEMENT("RANSEVMLowReEpsilon3D4N", mRANSEVMLowReEpsilon3D);

    KRATOS_REGISTER_ELEMENT("RANSEVMK2D3N", mRANSEVMK2D);
    KRATOS_REGISTER_ELEMENT("RANSEVMK3D4N", mRANSEVMK3D);
    KRATOS_REGISTER_ELEMENT("RANSEVMEpsilon2D3N", mRANSEVMEpsilon2D);
    KRATOS_REGISTER_ELEMENT("RANSEVMEpsilon3D4N", mRANSEVMEpsilon3D);

    KRATOS_REGISTER_CONDITION("RansEvmEpsilonWallCondition2D2N", mRansEvmEpsilonWallCondition2D2N);
    KRATOS_REGISTER_CONDITION("RansEvmEpsilonWallCondition3D3N", mRansEvmEpsilonWallCondition3D3N);

    KRATOS_REGISTER_CONDITION("RansEvmVmsMonolithicWallCondition2D2N",
                              mRansEvmVmsMonolithicWallCondition2D2N);
    KRATOS_REGISTER_CONDITION("RansEvmVmsMonolithicWallCondition3D3N",
                              mRansEvmVmsMonolithicWallCondition3D3N);
}
} // namespace Kratos.
