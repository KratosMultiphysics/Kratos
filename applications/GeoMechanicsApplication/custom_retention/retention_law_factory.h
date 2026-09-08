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
#include "custom_retention/retention_law.h"
#include "includes/define.h"

namespace Kratos
{
class Properties;

/**
 * @class RetentionLawFactory
 * @ingroup GeoMechanicsApplication
 * @author Vahid Galavi
 */
class KRATOS_API(GEO_MECHANICS_APPLICATION) RetentionLawFactory
{
public:
    /// Counted pointer of RetentionLawFactory
    KRATOS_CLASS_POINTER_DEFINITION(RetentionLawFactory);

    static std::unique_ptr<RetentionLaw> Clone(const Properties& rMaterialProperties);

}; // Class RetentionLawFactory
} // namespace Kratos.
