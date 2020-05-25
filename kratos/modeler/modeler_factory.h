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
#include "includes/kratos_components.h"
#include "modeler/modeler.h"

namespace Kratos
{
///@name Kratos Classes
///@{

/// Creates Modelers.
namespace ModelerFactory
{

    /// Checks if the modeler is registered
    bool Has(const std::string& ModelerName)
    {
        return KratosComponents< Modeler >::Has(ModelerName);
    }

    /// Checks if the modeler is registered
    typename Modeler::Pointer Create(
        const std::string& ModelerName, Model& rModel, const Parameters ModelParameters)
    {
        KRATOS_ERROR_IF_NOT(Has(ModelerName))
            << "Trying to construct a modeler: "
            << ModelerName << "\" which does not exist.\n"
            << "The available options (for currently loaded applications) are:\n"
            << KratosComponents< Modeler >() << std::endl;

        Modeler const& r_clone_modeler = KratosComponents< Modeler >::Get(ModelerName);
        return r_clone_modeler.Create(rModel, ModelParameters);
    }

}

}  // namespace Kratos.

#endif // KRATOS_MODELER_FACTORY_H_INCLUDED  defined

