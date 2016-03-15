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











#include "includes/kernel.h"
#include "includes/kratos_version.h"


namespace Kratos
{
    Kernel::Kernel()
    {
	std::cout << " |  /           |             " << std::endl;
	std::cout << " ' /   __| _` | __|  _ \\   __|" << std::endl;
	std::cout << " . \\  |   (   | |   (   |\\__ \\ " << std::endl;
	std::cout << "_|\\_\\_|  \\__,_|\\__|\\___/ ____/" << std::endl;
    std::cout << "           Multi-Physics "<< KRATOS_VERSION << std::endl;

        mKratosApplication.RegisterVariables();
    }

    std::string Kernel::Info() const
    {
        return "kernel";
    }

    void Kernel::PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "kernel";
    }

    /// Print object's data.
    void Kernel::PrintData(std::ostream& rOStream) const
    {
        rOStream << "Variables:" << std::endl;
        KratosComponents<VariableData>().PrintData(rOStream);
        rOStream << std::endl;
        rOStream << "Elements:" << std::endl;
        KratosComponents<Element>().PrintData(rOStream);
        rOStream << std::endl;
        rOStream << "Conditions:" << std::endl;
        KratosComponents<Condition>().PrintData(rOStream);
    }
}



