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

class Properties;

class KRATOS_API(GEO_MECHANICS_APPLICATION) GeoThermalLaw
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(GeoThermalLaw);

    virtual ~GeoThermalLaw() = default;

    [[nodiscard]] virtual Matrix CalculateThermalDispersionMatrix(const Properties& rProp) const = 0;

private:
    friend class Serializer;

    virtual void save(Serializer& rSerializer) const = 0;

    virtual void load(Serializer& rSerializer) = 0;
};
} // namespace Kratos.
