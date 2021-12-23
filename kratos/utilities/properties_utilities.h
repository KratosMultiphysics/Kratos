//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_PROPERTIES_UTILITIES)
#define KRATOS_PROPERTIES_UTILITIES

// System includes

// External includes

// Project includes
#include "includes/properties.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{
    
/**
 * @namespace PropertiesUtilities
 * @ingroup KratosCore
 * @brief This namespace includes several utilities for the Properties
 * @author Vicente Mataix Ferrandiz
 */
namespace PropertiesUtilities
{
    /**
     * @brief This method copies the values from a properties to another
     * @param rOriginProperties The model of the problem to mesh
     */
    void KRATOS_API(KRATOS_CORE) CopyPropertiesValues(
        const Properties& rOriginProperties,
        Properties& rDestinatiionProperties
        );

}; // namespace PropertiesUtilities
}  // namespace Kratos
#endif /* KRATOS_PROPERTIES_UTILITIES defined */
