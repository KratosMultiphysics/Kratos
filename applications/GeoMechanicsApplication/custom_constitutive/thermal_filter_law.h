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

#include "includes/constitutive_law.h"

namespace Kratos
{

/**
 * @class GeoThermalFilterLaw
 * @ingroup GeoMechanicsApplication
 * @brief This class defines the thermal filter for heat cases
 */
class KRATOS_API(GEO_MECHANICS_APPLICATION) GeoThermalFilterLaw : public ConstitutiveLaw
{
public:
    /// Counted pointer of GeoThermalFilterLaw
    KRATOS_CLASS_POINTER_DEFINITION(GeoThermalFilterLaw);

    GeoThermalFilterLaw();

    explicit GeoThermalFilterLaw(std::size_t NumberOfDimensions);

    ConstitutiveLaw::Pointer Clone() const override;

    SizeType WorkingSpaceDimension() override;

    Matrix CalculateThermalFilterMatrix(const Properties& rProp, const ProcessInfo& rProcessInfo) const;

private:
    SizeType mNumberOfDimensions;

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw)
        rSerializer.save("NumberOfDimensions", mNumberOfDimensions);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw)
        rSerializer.load("NumberOfDimensions", mNumberOfDimensions);
    }
};
} // namespace Kratos.
