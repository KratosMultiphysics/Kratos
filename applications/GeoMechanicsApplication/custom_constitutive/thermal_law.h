// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Mohamed Nabi
//                   John van Esch
//                   Gennady Markelov
//

#pragma once

#include "geo_mechanics_application_variables.h"
#include "includes/serializer.h"

namespace Kratos
{

/**
 * @class GeoThermalDispersionLaw
 * @ingroup GeoMechanicsApplication
 * @brief This class defines the thermal dispersion for heat cases
 */
class KRATOS_API(GEO_MECHANICS_APPLICATION) GeoThermalLaw
{
public:
    /// Counted pointer of GeoThermalLaw
    KRATOS_CLASS_POINTER_DEFINITION(GeoThermalLaw);

    GeoThermalLaw() {}

    virtual ~GeoThermalLaw() = default;

    virtual Matrix CalculateThermalDispersionMatrix(const Properties&  rProp,
                                                    const ProcessInfo& rProcessInfo) const = 0;

private:
    friend class Serializer;

    virtual void save(Serializer& rSerializer) const {}

    virtual void load(Serializer& rSerializer) {}
};
} // namespace Kratos.
