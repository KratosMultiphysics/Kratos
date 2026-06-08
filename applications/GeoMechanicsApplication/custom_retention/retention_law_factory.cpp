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

// Project includes
#include "custom_retention/retention_law_factory.h"
#include "custom_retention/saturated_below_phreatic_level_law.h"
#include "custom_retention/saturated_law.h"
#include "custom_retention/van_genuchten_law.h"

// Application includes
#include "geo_mechanics_application_variables.h"

namespace Kratos
{
std::unique_ptr<RetentionLaw> RetentionLawFactory::Clone(const Properties& rMaterialProperties)
{
    if (rMaterialProperties.Has(RETENTION_LAW)) {
        const std::string& RetentionLawName = rMaterialProperties[RETENTION_LAW];
        if (RetentionLawName == "VanGenuchtenLaw") return std::make_unique<VanGenuchtenLaw>();

        if (RetentionLawName == "SaturatedLaw") return std::make_unique<SaturatedLaw>();

        if (RetentionLawName == "SaturatedBelowPhreaticLevelLaw")
            return std::make_unique<SaturatedBelowPhreaticLevelLaw>();

        if (RetentionLawName == "PressureFilterLaw") return std::make_unique<SaturatedLaw>();

        KRATOS_ERROR << "Undefined RETENTION_LAW! " << RetentionLawName << std::endl;

        return nullptr;
    }

    // default is saturated law
    return std::make_unique<SaturatedLaw>();
}

} // namespace Kratos.
