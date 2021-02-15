// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Vahid Galavi
//

#if !defined (KRATOS_RETENTION_LAW_FACTORY_H_INCLUDED)
#define  KRATOS_RETENTION_LAW_FACTORY_H_INCLUDED

// System includes
#include <string>
#include <iostream>
#include "includes/define.h"

// External includes

// Project includes
#include "includes/serializer.h"
#include "custom_retention/retention_law.h"
#include "custom_retention/van_genuchten_law.h"
#include "custom_retention/saturated_law.h"

// Application includes
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

/**
 * @class RetentionLawFactory
 * @ingroup GeoMechanicsApplication
 * @author Vahid Galavi
 */
class KRATOS_API(GEO_MECHANICS_APPLICATION) RetentionLawFactory
{
public:

    ///@name Type Definitions

    /// Counted pointer of RetentionLawFactory
    KRATOS_CLASS_POINTER_DEFINITION( RetentionLawFactory );

    ///@}
    ///@name Lyfe Cycle
    ///@{

    /**
     * @brief Default constructor.
     */
    RetentionLawFactory(){}

    /**
     * @brief Destructor.
     */
    ~RetentionLawFactory(){}

    static RetentionLaw::Pointer Clone(const Properties& rMaterialProperties)
     {
         if (rMaterialProperties.Has(RETENTION_LAW))
         {
            const std::string &RetentionLawName = rMaterialProperties[RETENTION_LAW];
            if (RetentionLawName == "VanGenuchtenLaw")
                return VanGenuchtenLaw::Clone();

            if (RetentionLawName == "SaturatedLaw")
                return SaturatedLaw::Clone();

            KRATOS_THROW_ERROR( std::invalid_argument, "undefined RETENTION_LAW!!", RetentionLawName )

            return nullptr;
         }

         // default is saturated law
         return SaturatedLaw::Clone();
     }

private:


}; // Class RetentionLawFactory
}  // namespace Kratos.
#endif // KRATOS_RETENTION_LAW_FACTORY_H_INCLUDED  defined
