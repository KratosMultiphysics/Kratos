//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:   Philipp Bucher
//
//

// System includes

// External includes

// Project includes
#include "includes/properties.h"
#include "utilities/properties_configuration.h"

namespace Kratos
{

Properties::ConfigurationType& Properties::GetAccessor()
{
    if (!mpConfiguration)
        mpConfiguration = Kratos::make_shared<PropertiesConfiguration>(this);
    return *mpConfiguration;
}

const Properties::ConfigurationType& Properties::GetAccessor() const
{
    KRATOS_DEBUG_ERROR_IF_NOT(mpConfiguration) << "PropertiesConfiguration "
        << "is not initialized!" << std::endl;
    return *mpConfiguration;
}

}  // namespace Kratos.
