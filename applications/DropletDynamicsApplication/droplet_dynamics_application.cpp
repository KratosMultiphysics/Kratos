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
    mDropletDynamics2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
    mDropletDynamics3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4))))
    {}

void KratosDropletDynamicsApplication::Register()
{
     KRATOS_INFO("") << "Initializing KratosDropletDynamicsApplication..." << std::endl;

    // External interfacial force, e.g. for including the electromagentic coupling
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(EXT_INT_FORCE)

    KRATOS_REGISTER_VARIABLE( EPOTENTIAL )
    KRATOS_REGISTER_VARIABLE( PERMITTIVITYPOS )
    KRATOS_REGISTER_VARIABLE( PERMITTIVITYNEG )
    KRATOS_REGISTER_VARIABLE( CONDUCTIVITYPOS )
    KRATOS_REGISTER_VARIABLE( CONDUCTIVITYNEG )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( EFIELD )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( SCHARGE )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( EFORCE )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( EFIELDPOS )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( EFIELDNEG )


    KRATOS_REGISTER_VARIABLE( INV_K_ENRICH )
    KRATOS_REGISTER_VARIABLE( BIJ_ENRICH_ROW )
    KRATOS_REGISTER_VARIABLE( BJI_ENRICH_ROW )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( POS_GRAD_ENRICH )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( NEG_GRAD_ENRICH )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( RHS_ENRICH )
    
        // Register Elements
    KRATOS_REGISTER_ELEMENT("DropletDynamics2D3N", mDropletDynamics2D3N);
    KRATOS_REGISTER_ELEMENT("DropletDynamics3D4N", mDropletDynamics3D4N);
}
}  // namespace Kratos.


