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

#include "custom_retention/retention_law.h"
#include "includes/define.h"
#include "includes/serializer.h"

namespace Kratos
{
/**
 * @class SaturatedLaw
 * @ingroup GeoMechanicsApplication
 * @brief This class defines a saturated Soil Water Characteristic Curve (retention curve)
 * @details This class derives from the base retention law
 * @author Vahid Galavi
 */
class KRATOS_API(GEO_MECHANICS_APPLICATION) SaturatedLaw : public RetentionLaw
{
public:
    using GeometryType = Geometry<Node>;

    // Counted pointer of SaturatedLaw
    KRATOS_CLASS_POINTER_DEFINITION(SaturatedLaw);

    [[nodiscard]] RetentionLaw::Pointer Clone() const override;

    double CalculateSaturation(Parameters& rParameters) const override;

    double CalculateEffectiveSaturation(Parameters& rParameters) const override;

    double CalculateDerivativeOfSaturation(Parameters& rParameters) const override;

    double CalculateRelativePermeability(Parameters& rParameters) const override;

    double CalculateBishopCoefficient(Parameters& rParameters) const override;

    /**
     * @brief This function provides the place to perform checks on the completeness of the input.
     * @details It is designed to be called only once (or anyway, not often) typically at the beginning of the calculations, so to verify that nothing is missing from the input or that no common error is found.
     * @param rMaterialProperties The properties of the material
     * @param rCurrentProcessInfo The current process info instance
     * @return 0 if OK, 1 otherwise
     */
    int Check(const Properties& rMaterialProperties, const ProcessInfo& rCurrentProcessInfo) override;

private:
    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, RetentionLaw)
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, RetentionLaw)
    }

}; // Class SaturatedLaw
} // namespace Kratos.
