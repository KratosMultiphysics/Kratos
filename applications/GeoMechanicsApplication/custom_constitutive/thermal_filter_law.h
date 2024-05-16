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

/**
 * @class GeoThermalFilterLaw
 * @ingroup GeoMechanicsApplication
 * @brief This class defines the thermal filter for heat cases
 */
class KRATOS_API(GEO_MECHANICS_APPLICATION) GeoThermalFilterLaw : public GeoThermalLaw
{
public:
    /// Counted pointer of GeoThermalFilterLaw
    KRATOS_CLASS_POINTER_DEFINITION(GeoThermalFilterLaw);

    GeoThermalFilterLaw();

    explicit GeoThermalFilterLaw(SizeType NumberOfDimensions);

    ~GeoThermalFilterLaw() override = default;

    GeoThermalLaw::Pointer Clone() const override;

    Matrix CalculateThermalDispersionMatrix(const Properties& rProp, const ProcessInfo&) const override;

private:
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
