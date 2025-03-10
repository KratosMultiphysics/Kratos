//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Mohammad R. Hashemi
//


// System includes


// External includes


// Project includes
#include "droplet_dynamics_application.h"
#include "droplet_dynamics_application_variables.h"


namespace Kratos {

KratosDropletDynamicsApplication::KratosDropletDynamicsApplication():
    KratosApplication("DropletDynamicsApplication"),
    mDropletDynamics2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<Node >(Element::GeometryType::PointsArrayType(3)))),
    mDropletDynamics3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node >(Element::GeometryType::PointsArrayType(4))))
    {}

void KratosDropletDynamicsApplication::Register()
{
     KRATOS_INFO("") << "Initializing KratosDropletDynamicsApplication..." << std::endl;

    // External interfacial force, e.g. for including the electromagentic coupling
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(EXT_INT_FORCE)
    
    //Auxiliary variable to store maximum element size (h_{max})
    KRATOS_REGISTER_VARIABLE(NODAL_H_MAX)

    // Smoothing surface auxiliary distance
    KRATOS_REGISTER_VARIABLE( DISTANCE_AUX)
    KRATOS_REGISTER_VARIABLE( DISTANCE_AUX2)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( DISTANCE_GRADIENT_AUX)
    
    // Parallel levelset distance calculator needs an AREA_VARIABLE_AUX
    KRATOS_REGISTER_VARIABLE( AREA_VARIABLE_AUX)

    // A variable to check if node is on cut element (maybe in a layer farther for future!)
    //KRATOS_REGISTER_VARIABLE( IS_NEAR_CUT)

    // Contact line calculation
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(NORMAL_VECTOR)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(TANGENT_VECTOR)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(CONTACT_VECTOR)
    KRATOS_REGISTER_VARIABLE( CONTACT_ANGLE)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(CONTACT_VECTOR_MICRO)
    KRATOS_REGISTER_VARIABLE( CONTACT_ANGLE_MICRO)
    KRATOS_REGISTER_VARIABLE( CONTACT_VELOCITY)


    // Enriched pressure is an array of NumNodes components defined for elements. Access it using Element.GetValue()
    KRATOS_REGISTER_VARIABLE( ENRICHED_PRESSURE_1)
    KRATOS_REGISTER_VARIABLE( ENRICHED_PRESSURE_2)
    KRATOS_REGISTER_VARIABLE( ENRICHED_PRESSURE_3)
    KRATOS_REGISTER_VARIABLE( ENRICHED_PRESSURE_4)

    // Last known velocity and pressure to recalculate the last increment
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( VELOCITY_STAR)
    KRATOS_REGISTER_VARIABLE( PRESSURE_STAR)

    // Pressure gradient to calculate its jump over interface
    // KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( PRESSURE_GRADIENT_AUX)

    // Register Elements
    KRATOS_REGISTER_ELEMENT("DropletDynamics2D3N", mDropletDynamics2D3N);
    KRATOS_REGISTER_ELEMENT("DropletDynamics3D4N", mDropletDynamics3D4N);
}
}  // namespace Kratos.
