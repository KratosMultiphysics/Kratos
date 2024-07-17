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

#pragma once

// System includes
#include <string>
#include <iostream>
#include "includes/define.h"

// External includes

// Project includes
#include "custom_retention/retention_law.h"
#include "custom_retention/van_genuchten_law.h"
#include "custom_retention/saturated_law.h"
#include "custom_retention/saturated_below_phreatic_level_law.h"

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
    /// Counted pointer of RetentionLawFactory
    KRATOS_CLASS_POINTER_DEFINITION( RetentionLawFactory );

    static unique_ptr<RetentionLaw> Clone(const Properties& rMaterialProperties)
     {
         if (rMaterialProperties.Has(RETENTION_LAW))
         {
            const std::string &RetentionLawName = rMaterialProperties[RETENTION_LAW];
            if (RetentionLawName == "VanGenuchtenLaw")
                return make_unique<VanGenuchtenLaw>();

            if (RetentionLawName == "SaturatedLaw")
                return make_unique<SaturatedLaw>();

            if (RetentionLawName == "SaturatedBelowPhreaticLevelLaw")
                return make_unique<SaturatedBelowPhreaticLevelLaw>();

            if (RetentionLawName == "PressureFilterLaw") 
                return make_unique<SaturatedLaw>();

            KRATOS_ERROR << "Undefined RETENTION_LAW! "
                         << RetentionLawName
                         << std::endl;

            return nullptr;
         }

         // default is saturated law
         return make_unique<SaturatedLaw>();

     }

}; // Class RetentionLawFactory
}  // namespace Kratos.
