//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//


// System includes


// External includes


// Project includes
#include "nurbs_brep_application.h"
#include "nurbs_brep_application_variables.h"


namespace Kratos {

KratosNurbsBrepApplication::KratosNurbsBrepApplication()
    : KratosApplication("NurbsBrepApplication")
{}

void KratosNurbsBrepApplication::Register() {
    // calling base class register to register Kratos components
    KratosApplication::Register();
    std::cout << "Initializing KratosNurbsBrepApplication... " << std::endl;
}

} // namespace Kratos
