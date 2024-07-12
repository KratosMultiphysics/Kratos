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

    // Register Elements
    KRATOS_REGISTER_ELEMENT("DropletDynamics2D3N", mDropletDynamics2D3N);
    KRATOS_REGISTER_ELEMENT("DropletDynamics3D4N", mDropletDynamics3D4N);
}
}  // namespace Kratos.
