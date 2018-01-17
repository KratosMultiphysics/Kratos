//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//
#include <iostream>

#include "includes/kernel.h"
#include "includes/kratos_version.h"
#include "input_output/logger.h"

namespace Kratos {
Kernel::Kernel() {
    std::cout << " |  /           |             " << std::endl;
    std::cout << " ' /   __| _` | __|  _ \\   __|" << std::endl;
    std::cout << " . \\  |   (   | |   (   |\\__ \\ " << std::endl;
    std::cout << "_|\\_\\_|  \\__,_|\\__|\\___/ ____/" << std::endl;
    std::cout << "           Multi-Physics " << KRATOS_VERSION << std::endl;

    if (!IsImported("KratosMultiphysics")) {
        KratosApplication::Pointer pKratosApplication =
            Kratos::make_shared<KratosApplication>(
                std::string("KratosMultiphysics"));
        pKratosApplication->RegisterVariables();
        this->ImportApplication(pKratosApplication);
    }
}

std::unordered_map<std::string, KratosApplication::Pointer>&
Kernel::GetApplicationsList() {
    static std::unordered_map<std::string, KratosApplication::Pointer>
        application_list;
    return application_list;
}

bool Kernel::IsImported(std::string ApplicationName) const {
    if (GetApplicationsList().find(ApplicationName) !=
        GetApplicationsList().end())
        return true;
    else
        return false;
}

void Kernel::ImportApplication(KratosApplication::Pointer pNewApplication) {
    if (IsImported(pNewApplication->Name()))
        KRATOS_ERROR << "importing more than once the application : "
                     << pNewApplication->Name() << std::endl;

    pNewApplication->Register();
    Kernel::GetApplicationsList()[pNewApplication->Name()] = pNewApplication;
}

std::string Kernel::Info() const { return "kernel"; }

void Kernel::PrintInfo(std::ostream& rOStream) const { rOStream << "kernel"; }

/// Print object's data.
void Kernel::PrintData(std::ostream& rOStream) const {
    rOStream << "Variables:" << std::endl;
    KratosComponents<VariableData>().PrintData(rOStream);
    rOStream << std::endl;
    rOStream << "Elements:" << std::endl;
    KratosComponents<Element>().PrintData(rOStream);
    rOStream << std::endl;
    rOStream << "Conditions:" << std::endl;
    KratosComponents<Condition>().PrintData(rOStream);

    rOStream << "Loaded applications:" << std::endl;

    auto& application_list = Kernel::GetApplicationsList();
    rOStream << "number of loaded applications = " << application_list.size()
             << std::endl;
    for (auto it = application_list.begin(); it != application_list.end(); ++it)
        rOStream << "  " << it->first << std::endl;
}
}
