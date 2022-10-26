//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//


// System includes


// External includes


// Project includes
#include "med_application.h"


namespace Kratos {

KratosMedApplication::KratosMedApplication():
    KratosApplication("MedApplication")
    {}

void KratosMedApplication::Register()
{
     KRATOS_INFO("") << "Initializing KratosMedApplication..." << std::endl;

}

}  // namespace Kratos.
