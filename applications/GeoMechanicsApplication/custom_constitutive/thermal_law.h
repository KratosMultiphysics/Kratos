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

/**
 * @class GeoThermalDispersionLaw
 * @ingroup GeoMechanicsApplication
 * @brief This class defines the thermal dispersion for heat cases
 */
class KRATOS_API(GEO_MECHANICS_APPLICATION) GeoThermalLaw
{
public:
    /// Counted pointer of GeoThermalLaw
    KRATOS_CLASS_POINTER_DEFINITION(GeoThermalLaw);

    GeoThermalLaw() {}

    virtual ~GeoThermalLaw() = default;

    virtual GeoThermalLaw::Pointer Clone() const = 0;

    std::size_t WorkingSpaceDimension() const { return mNumberOfDimensions; }

    virtual Matrix CalculateThermalDispersionMatrix(const Properties&  rProp,
                                                    const ProcessInfo& rProcessInfo) const = 0;

protected:
    std::size_t mNumberOfDimensions = 0;

private:
    friend class Serializer;

    virtual void save(Serializer& rSerializer) const
    {
        rSerializer.save("NumberOfDimensions", mNumberOfDimensions);
    }

    virtual void load(Serializer& rSerializer)
    {
        rSerializer.load("NumberOfDimensions", mNumberOfDimensions);
    }
};
} // namespace Kratos.
