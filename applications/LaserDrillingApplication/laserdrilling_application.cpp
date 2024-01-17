//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ignasi de Pouplana
//


// System includes

// External includes
//

// Project includes
#include "includes/define.h"
#include "laserdrilling_application.h"
#include "includes/variables.h"

namespace Kratos {



KratosLaserDrillingApplication::KratosLaserDrillingApplication()
    : KratosApplication("LaserDrillingApplication")
    //   mEulerianConvDiff2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
    //   mDLaserDrillingExplicit3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))) 
      {}

void KratosLaserDrillingApplication::Register() {
    KRATOS_INFO("") << "Initializing KratosLaserDrillingApplication... " << std::endl;

    // Registering variables
    KRATOS_REGISTER_VARIABLE(NO2)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(TRANSPORT_VELOCITY)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(VECTOR_VELOCITAT_VENT)
    KRATOS_REGISTER_VARIABLE(IS_POTENTIAL_FLOW_STEP)
    KRATOS_REGISTER_VARIABLE(VELOCITAT_VENT)
    KRATOS_REGISTER_VARIABLE(DIRECCIO_VENT)
    KRATOS_REGISTER_VARIABLE(METEO_DIRECCIO_VENT)
    KRATOS_REGISTER_VARIABLE(STARTING_DATE)
    KRATOS_REGISTER_VARIABLE(SIMULATION_DURATION_IN_DAYS)
    KRATOS_REGISTER_VARIABLE(WIND_AUTOMATIC_PROCESS)
    KRATOS_REGISTER_VARIABLE(POLLUTANT_AUTOMATIC_PROCESS)
    KRATOS_REGISTER_VARIABLE(CITY)
    KRATOS_REGISTER_VARIABLE(CASE_ID)
    KRATOS_REGISTER_VARIABLE(IN_PRODUCTION)

    // Registering elements and conditions here
    // KRATOS_REGISTER_ELEMENT("EulerianConvDiff2D", mEulerianConvDiff2D);
    // KRATOS_REGISTER_ELEMENT("DLaserDrillingExplicit3D4N", mDLaserDrillingExplicit3D4N);
}

}  // namespace Kratos.
