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

/**
 * @class GeoThermalDispersionLaw
 * @ingroup GeoMechanicsApplication
 * @brief This class defines the thermal dispersion for heat cases
 */
class KRATOS_API(GEO_MECHANICS_APPLICATION) GeoThermalDispersionLaw : public GeoThermalLaw
{
public:
    /// Counted pointer of GeoThermalDispersionLaw
    KRATOS_CLASS_POINTER_DEFINITION(GeoThermalDispersionLaw);

    GeoThermalDispersionLaw();

    explicit GeoThermalDispersionLaw(SizeType NumberOfDimensions);

    ~GeoThermalDispersionLaw() override = default;

    Matrix CalculateThermalDispersionMatrix(const Properties& rProp, const ProcessInfo& rProcessInfo) const override;

private:
    std::size_t mNumberOfDimensions = 0;

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, GeoThermalLaw)
        rSerializer.save("NumberOfDimensions", mNumberOfDimensions);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, GeoThermalLaw)
        rSerializer.load("NumberOfDimensions", mNumberOfDimensions);
    }
};
} // namespace Kratos.
