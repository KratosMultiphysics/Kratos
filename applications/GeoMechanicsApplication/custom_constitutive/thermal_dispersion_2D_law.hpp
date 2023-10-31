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

namespace Kratos {

/**
 * @class GeoThermalDispersion2DLaw
 * @ingroup GeoMechanicsApplication
 * @brief This class defines the thermal dispersion for heat cases
 * @details This class derives from the linear elastic case on 3D
 */
class KRATOS_API(GEO_MECHANICS_APPLICATION) GeoThermalDispersion2DLaw : public ConstitutiveLaw {
public:

    /// Counted pointer of LinearPlaneStrainK0Law
    KRATOS_CLASS_POINTER_DEFINITION(GeoThermalDispersion2DLaw);

    ConstitutiveLaw::Pointer Clone() const override;

    SizeType WorkingSpaceDimension() override;

    static void CalculateThermalDispersionMatrix(Matrix& rThermalDispersionMatrix,
                                                 const Properties& rValues);

private:
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
