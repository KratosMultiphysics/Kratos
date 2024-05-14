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

namespace Kratos
{

class RetentionLaw;

/**
 * @class GeoPwFilterLaw
 * @ingroup GeoMechanicsApplication
 * @brief This class defines the Pw filter for water pressure cases
 */
class KRATOS_API(GEO_MECHANICS_APPLICATION) GeoPwFilterLaw : public ConstitutiveLaw
{
public:
    /// Counted pointer of GeoPwFilterLaw
    KRATOS_CLASS_POINTER_DEFINITION(GeoPwFilterLaw);

    GeoPwFilterLaw();

    explicit GeoPwFilterLaw(SizeType NumberOfDimensions);

    ConstitutiveLaw::Pointer Clone() const override;

    SizeType WorkingSpaceDimension() override;

    Matrix CalculatePwFilterMatrix(const Properties& rProp, const ProcessInfo& rProcessInfo) const;

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
