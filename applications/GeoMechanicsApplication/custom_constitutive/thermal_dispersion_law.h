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

#include "custom_constitutive/thermal_law.h"

namespace Kratos
{

class RetentionLaw;

class KRATOS_API(GEO_MECHANICS_APPLICATION) GeoThermalDispersionLaw : public GeoThermalLaw
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(GeoThermalDispersionLaw);

    GeoThermalDispersionLaw() = default;

    explicit GeoThermalDispersionLaw(SizeType NumberOfDimensions);

    [[nodiscard]] Matrix CalculateThermalDispersionMatrix(const Properties& rProp) const override;

private:
    std::size_t mNumberOfDimensions = 2;

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        rSerializer.save("NumberOfDimensions", mNumberOfDimensions);
    }

    void load(Serializer& rSerializer) override
    {
        rSerializer.load("NumberOfDimensions", mNumberOfDimensions);
    }
};
} // namespace Kratos.
