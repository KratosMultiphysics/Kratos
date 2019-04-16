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
#include "rans_modelling_application_variables.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/tetrahedra_3d_4.h"

namespace Kratos {

KratosRANSModellingApplication::KratosRANSModellingApplication():
    KratosApplication("RANSModellingApplication"),
    mRANSEVMK2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
    mRANSEVMK3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    mRANSEVMEPSILON2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
    mRANSEVMEPSILON3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4))))
    {}

void KratosRANSModellingApplication::Register()
{
     // calling base class register to register Kratos components
     KratosApplication::Register();
     KRATOS_INFO("") << "Initializing KratosRANSModellingApplication..." << std::endl;

  KRATOS_REGISTER_VARIABLE( TURBULENT_KINETIC_ENERGY )
  KRATOS_REGISTER_VARIABLE( TURBULENT_ENERGY_DISSIPATION_RATE )
  KRATOS_REGISTER_VARIABLE( TURBULENT_KINETIC_ENERGY_RATE )
  KRATOS_REGISTER_VARIABLE( TURBULENT_ENERGY_DISSIPATION_RATE_2 )
  KRATOS_REGISTER_VARIABLE( IS_CO_SOLVING_PROCESS_ACTIVE )
  KRATOS_REGISTER_VARIABLE( OLD_CONVERGENCE_VARIABLE )
  KRATOS_REGISTER_VARIABLE( RANS_Y_PLUS )
  KRATOS_REGISTER_VARIABLE( RANS_AUXILIARY_VARIABLE_1 )
  KRATOS_REGISTER_VARIABLE( RANS_AUXILIARY_VARIABLE_2 )
  KRATOS_REGISTER_VARIABLE( WALL_SMOOTHNESS_BETA )
  KRATOS_REGISTER_VARIABLE( WALL_VON_KARMAN )
  KRATOS_REGISTER_VARIABLE( TURBULENCE_RANS_C_MU )
  KRATOS_REGISTER_VARIABLE( TURBULENCE_RANS_C1 )
  KRATOS_REGISTER_VARIABLE( TURBULENCE_RANS_C2 )
  KRATOS_REGISTER_VARIABLE( TURBULENT_VISCOSITY_MIN )
  KRATOS_REGISTER_VARIABLE( TURBULENT_VISCOSITY_MAX )
  KRATOS_REGISTER_VARIABLE( TURBULENT_KINETIC_ENERGY_SIGMA )
  KRATOS_REGISTER_VARIABLE( TURBULENT_ENERGY_DISSIPATION_RATE_SIGMA )

  // Register Elements
  KRATOS_REGISTER_ELEMENT("RANSEVMK2D3N",mRANSEVMK2D);
  KRATOS_REGISTER_ELEMENT("RANSEVMK3D4N",mRANSEVMK3D);
  KRATOS_REGISTER_ELEMENT("RANSEVMEPSILON2D3N",mRANSEVMEPSILON2D);
  KRATOS_REGISTER_ELEMENT("RANSEVMEPSILON3D4N",mRANSEVMEPSILON3D);

}
}  // namespace Kratos.
