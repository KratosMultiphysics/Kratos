//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//


#if !defined(KRATOS_MODELER_FACTORY_H_INCLUDED )
#define  KRATOS_MODELER_FACTORY_H_INCLUDED


// System includes

// External includes

// Project includes
#include "modeler/modeler.h"

namespace Kratos
{
///@name Kratos Classes
///@{

/// Creates Modelers.
namespace ModelerFactory
{

    /// Checks if the modeler is registered
    bool KRATOS_API(KRATOS_CORE) Has(const std::string& ModelerName);

    /// Checks if the modeler is registered
    typename Modeler::Pointer KRATOS_API(KRATOS_CORE) Create(
        const std::string& ModelerName, Model& rModel, const Parameters ModelParameters);

}

}  // namespace Kratos.

#endif // KRATOS_MODELER_FACTORY_H_INCLUDED  defined

