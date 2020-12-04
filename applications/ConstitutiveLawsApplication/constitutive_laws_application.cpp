//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Alejandro Cornejo Velazquez
//                   Riccardo Rossi
//


// System includes


// External includes


// Project includes
#include "constitutive_laws_application.h"
#include "constitutive_laws_application_variables.h"


namespace Kratos {

KratosConstitutiveLawsApplication::KratosConstitutiveLawsApplication():
    KratosApplication("ConstitutiveLawsApplication")
    {}

void KratosConstitutiveLawsApplication::Register()
{
     KRATOS_INFO("") << "Initializing KratosConstitutiveLawsApplication..." << std::endl;

}
}  // namespace Kratos.
