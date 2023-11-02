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

#include "includes/constitutive_law.h"

namespace Kratos {

/**
 * @class GeoThermalDispersionLaw
 * @ingroup GeoMechanicsApplication
 * @brief This class defines the thermal dispersion for heat cases
 * @details This class derives from the linear elastic case on 3D
 */
class KRATOS_API(GEO_MECHANICS_APPLICATION) GeoThermalDispersionLaw : public ConstitutiveLaw {
public:
    /// Counted pointer of GeoThermalDispersionLaw
    KRATOS_CLASS_POINTER_DEFINITION(GeoThermalDispersionLaw);

    explicit GeoThermalDispersionLaw(SizeType NumberOfDimensions);

    ConstitutiveLaw::Pointer Clone() const override;

    SizeType WorkingSpaceDimension() override;

    Matrix CalculateThermalDispersionMatrix(const Properties& rValues) const;

private:
    SizeType mNumberOfDimensions;

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw)
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw)
    }
};
} // namespace Kratos.
