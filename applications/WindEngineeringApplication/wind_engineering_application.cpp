//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Máté Kelemen
//

// Project includes
#include "wind_engineering_application.h"


namespace Kratos
{


KratosWindEngineeringApplication::KratosWindEngineeringApplication()
    : KratosApplication("WindEngineeringApplication")
{
}


void KratosWindEngineeringApplication::Register()
{
    KRATOS_INFO("") << "Initializing KratosWindEngineeringApplication..." << std::endl;
    // Register custom variables
}


std::string KratosWindEngineeringApplication::Info() const
{
    return "KratosWindEngineeringApplication";
}


void KratosWindEngineeringApplication::PrintInfo(std::ostream& rStream) const
{
    rStream << this->Info();
    this->PrintData(rStream);
}


void KratosWindEngineeringApplication::PrintData(std::ostream& rStream) const
{
    KRATOS_WATCH("in my application");
    KRATOS_WATCH(KratosComponents<VariableData>::GetComponents().size() );

    rStream << "Variables:" << std::endl;
    KratosComponents<VariableData>().PrintData(rStream);
    rStream << std::endl;
    rStream << "Elements:" << std::endl;
    KratosComponents<Element>().PrintData(rStream);
    rStream << std::endl;
    rStream << "Conditions:" << std::endl;
    KratosComponents<Condition>().PrintData(rStream);
}


} // namespace Kratos