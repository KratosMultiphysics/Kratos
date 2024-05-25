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
//

#pragma once

#include "custom_constitutive/thermal_law.h"

namespace Kratos
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) GeoThermalFilterLaw : public GeoThermalLaw
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(GeoThermalFilterLaw);

    [[nodiscard]] Matrix CalculateThermalDispersionMatrix(const Properties& rProp) const override;

private:
    friend class Serializer;

    void save(Serializer&) const override
    {
        // Nothing to serialize, since there are no data members, and its base class is abstract
    }

    void load(Serializer&) override
    {
        // Nothing to serialize, since there are no data members, and its base class is abstract
    }
};
} // namespace Kratos.
